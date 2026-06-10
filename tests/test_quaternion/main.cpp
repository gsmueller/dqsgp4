/// test_quaternion — Algebraic correctness for math::Quaternion<T>.
///
/// Implements design/audit/mathematical_correctness.md AUD-MC-4..14,
/// verifying REQ-DQ-4..14 from
/// design/specifications/dual_quaternion_algebra.md.
///
/// Each test draws a deterministic sample of inputs from a seeded
/// std::mt19937_64, checks the identity numerically, and asserts that
/// the actual deviation is dominated by the sum of the two sides'
/// reported `total_error()`. That double duty (identity holds AND the
/// reported bound is valid) implements both the algebraic-correctness
/// audit and the error-framework dominance check (REQ-EF-2).

#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>

using T = double;
using TV = math::TrackedValue<T>;
using math::Quaternion;
using math::Vector3;

namespace {

int tests_passed = 0;
int tests_failed = 0;

constexpr std::uint64_t kSeed = 0xC0FFEEull;
constexpr int kSamples = 8;

// --- Helpers ---

TV tv(double v) {
    T val = static_cast<T>(v);
    T prec = TV::representation_bound(val);
    return TV(val, T(0), prec, T(0));
}

Quaternion<T> random_unit_q(std::mt19937_64& rng) {
    std::normal_distribution<T> dist(T(0), T(1));
    T w = dist(rng), x = dist(rng), y = dist(rng), z = dist(rng);
    T mag = std::sqrt(w * w + x * x + y * y + z * z);
    return Quaternion<T>(tv(w / mag), tv(x / mag), tv(y / mag), tv(z / mag));
}

Vector3<T> random_v(std::mt19937_64& rng) {
    std::normal_distribution<T> dist(T(0), T(1));
    return Vector3<T>(tv(dist(rng)), tv(dist(rng)), tv(dist(rng)));
}

T q_diff(const Quaternion<T>& a, const Quaternion<T>& b) {
    T dw = std::abs(a.w.value - b.w.value);
    T dx = std::abs(a.x.value - b.x.value);
    T dy = std::abs(a.y.value - b.y.value);
    T dz = std::abs(a.z.value - b.z.value);
    return std::max({dw, dx, dy, dz});
}

T v_diff(const Vector3<T>& a, const Vector3<T>& b) {
    T dx = std::abs(a.x.value - b.x.value);
    T dy = std::abs(a.y.value - b.y.value);
    T dz = std::abs(a.z.value - b.z.value);
    return std::max({dx, dy, dz});
}

T q_pair_bound(const Quaternion<T>& a, const Quaternion<T>& b) {
    T bw = a.w.total_error() + b.w.total_error();
    T bx = a.x.total_error() + b.x.total_error();
    T by = a.y.total_error() + b.y.total_error();
    T bz = a.z.total_error() + b.z.total_error();
    return std::max({bw, bx, by, bz});
}

T v_pair_bound(const Vector3<T>& a, const Vector3<T>& b) {
    T bx = a.x.total_error() + b.x.total_error();
    T by = a.y.total_error() + b.y.total_error();
    T bz = a.z.total_error() + b.z.total_error();
    return std::max({bx, by, bz});
}

void check_bound(const char* name, T measured, T bound) {
    if (measured <= bound) {
        ++tests_passed;
    } else {
        ++tests_failed;
        std::cerr << "  FAIL: " << name
                  << "  measured=" << std::setprecision(17) << measured
                  << "  bound=" << bound << "\n";
    }
}

// --- Tests ---

void test_associativity() {
    std::cout << "--- AUD-MC-4 (REQ-DQ-4): Hamilton-product associativity ---\n";
    std::mt19937_64 rng(kSeed);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> a = random_unit_q(rng);
        Quaternion<T> b = random_unit_q(rng);
        Quaternion<T> c = random_unit_q(rng);
        Quaternion<T> lhs = (a * b) * c;
        Quaternion<T> rhs = a * (b * c);
        check_bound("(ab)c == a(bc)", q_diff(lhs, rhs), q_pair_bound(lhs, rhs));
    }
}

void test_identity() {
    std::cout << "--- AUD-MC-5 (REQ-DQ-5): identity laws ---\n";
    std::mt19937_64 rng(kSeed + 1);
    Quaternion<T> id = Quaternion<T>::identity();
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q = random_unit_q(rng);
        Quaternion<T> r = q * id;
        Quaternion<T> l = id * q;
        check_bound("q * identity == q", q_diff(r, q), q_pair_bound(r, q));
        check_bound("identity * q == q", q_diff(l, q), q_pair_bound(l, q));
    }
}

void test_conjugate_involutive() {
    std::cout << "--- AUD-MC-6 (REQ-DQ-6): conjugate is involutive ---\n";
    std::mt19937_64 rng(kSeed + 2);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q = random_unit_q(rng);
        Quaternion<T> qcc = q.conjugate().conjugate();
        check_bound("(q*)* == q", q_diff(qcc, q), q_pair_bound(qcc, q));
    }
}

void test_conjugate_of_product() {
    std::cout << "--- AUD-MC-7 (REQ-DQ-7): (ab)* == b* a* ---\n";
    std::mt19937_64 rng(kSeed + 3);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> a = random_unit_q(rng);
        Quaternion<T> b = random_unit_q(rng);
        Quaternion<T> lhs = (a * b).conjugate();
        Quaternion<T> rhs = b.conjugate() * a.conjugate();
        check_bound("(ab)* == b* a*", q_diff(lhs, rhs), q_pair_bound(lhs, rhs));
    }
}

void test_magnitude_multiplicative() {
    std::cout << "--- AUD-MC-8 (REQ-DQ-8): |ab| == |a| |b| ---\n";
    std::mt19937_64 rng(kSeed + 4);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> a = random_unit_q(rng);
        Quaternion<T> b = random_unit_q(rng);
        TV m_ab = (a * b).magnitude();
        TV m_a = a.magnitude();
        TV m_b = b.magnitude();
        TV product = m_a * m_b;
        T diff = std::abs(m_ab.value - product.value);
        T bound = m_ab.total_error() + product.total_error();
        check_bound("|ab| == |a||b|", diff, bound);
    }
}

void test_inverse() {
    std::cout << "--- AUD-MC-9 (REQ-DQ-9): inverse identity ---\n";
    std::mt19937_64 rng(kSeed + 5);
    Quaternion<T> id = Quaternion<T>::identity();
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q = random_unit_q(rng);
        Quaternion<T> inv = q.inverse();
        Quaternion<T> r = q * inv;
        Quaternion<T> l = inv * q;
        check_bound("q * inverse == identity", q_diff(r, id), q_pair_bound(r, id));
        check_bound("inverse * q == identity", q_diff(l, id), q_pair_bound(l, id));
    }
}

void test_rotation_preserves_length() {
    std::cout << "--- AUD-MC-10 (REQ-DQ-10): rotation preserves length ---\n";
    std::mt19937_64 rng(kSeed + 6);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q = random_unit_q(rng);
        Vector3<T> v = random_v(rng);
        Vector3<T> v_rot = q.rotate(v);
        TV m_orig = v.magnitude();
        TV m_rot = v_rot.magnitude();
        T diff = std::abs(m_orig.value - m_rot.value);
        T bound = m_orig.total_error() + m_rot.total_error();
        check_bound("|q.rotate(v)| == |v|", diff, bound);
    }
}

void test_rotation_composition() {
    std::cout << "--- AUD-MC-11 (REQ-DQ-11): rotation composition ---\n";
    std::mt19937_64 rng(kSeed + 7);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q1 = random_unit_q(rng);
        Quaternion<T> q2 = random_unit_q(rng);
        Vector3<T> v = random_v(rng);
        Vector3<T> lhs = (q1 * q2).rotate(v);
        Vector3<T> rhs = q1.rotate(q2.rotate(v));
        check_bound("(q1 q2)(v) == q1(q2(v))",
                    v_diff(lhs, rhs), v_pair_bound(lhs, rhs));
    }
}

void test_exp_log_roundtrip() {
    std::cout << "--- AUD-MC-12 (REQ-DQ-12): exp/log round-trip ---\n";
    std::mt19937_64 rng(kSeed + 8);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q = random_unit_q(rng);
        // Shortest-path: force w >= 0 so we land in log_unit's
        // recoverable branch.
        if (q.w.value < T(0)) {
            q = Quaternion<T>(-q.w, -q.x, -q.y, -q.z);
        }
        Vector3<T> v = q.log_unit();
        Quaternion<T> q_rt = Quaternion<T>::exp_pure(v);
        check_bound("exp(log(q)) == q", q_diff(q_rt, q), q_pair_bound(q_rt, q));
    }
}

void test_axis_angle_roundtrip() {
    std::cout << "--- AUD-MC-13 (REQ-DQ-13): axis-angle round-trip ---\n";
    std::mt19937_64 rng(kSeed + 9);
    std::uniform_real_distribution<T> angle_dist(T(0.1), T(3.0));
    for (int i = 0; i < kSamples; ++i) {
        Vector3<T> axis = random_v(rng);
        TV am = axis.magnitude();
        Vector3<T> axis_unit(axis.x / am, axis.y / am, axis.z / am);
        TV theta = tv(angle_dist(rng));
        Quaternion<T> q = Quaternion<T>::from_axis_angle(axis_unit, theta);
        Vector3<T> rv = q.log_unit();
        TV half_theta = theta / math::exact<T>(2);
        Vector3<T> expected(half_theta * axis_unit.x,
                            half_theta * axis_unit.y,
                            half_theta * axis_unit.z);
        check_bound("log(from_axis_angle(n, theta)) == (theta/2) n",
                    v_diff(rv, expected), v_pair_bound(rv, expected));
    }
}

void test_normalize_idempotent() {
    std::cout << "--- AUD-MC-14 (REQ-DQ-14): normalize idempotent on unit q ---\n";
    std::mt19937_64 rng(kSeed + 10);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q = random_unit_q(rng);
        Quaternion<T> qn = q.normalized();
        check_bound("normalize(unit q) == q", q_diff(qn, q), q_pair_bound(qn, q));
    }
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(17);
    std::cout << "test_quaternion: REQ-DQ-4..14 / AUD-MC-4..14\n\n";

    test_associativity();
    test_identity();
    test_conjugate_involutive();
    test_conjugate_of_product();
    test_magnitude_multiplicative();
    test_inverse();
    test_rotation_preserves_length();
    test_rotation_composition();
    test_exp_log_roundtrip();
    test_axis_angle_roundtrip();
    test_normalize_idempotent();

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << tests_passed
              << "  Failed: " << tests_failed << "\n";
    std::cout << "========================================\n";
    return tests_failed > 0 ? 1 : 0;
}
