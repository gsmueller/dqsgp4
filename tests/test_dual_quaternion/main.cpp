/// test_dual_quaternion — Algebraic correctness for math::DualQuaternion<T>.
///
/// Implements design/audit/mathematical_correctness.md AUD-MC-15..18,
/// verifying REQ-DQ-15..18 from
/// design/specifications/dual_quaternion_algebra.md:
///   AUD-MC-15: DQ composition associativity
///   AUD-MC-16: three conjugate dualities and the unit-DQ inverse identity
///   AUD-MC-17: pose action composition: (M1 M2)(v) == M1(M2(v))
///   AUD-MC-18: screw exp/log round-trip on the half-angle ball
///
/// Plus a `from_pose / translation` accessor round-trip as a separate
/// supporting check.

#include "math/dual_quaternion.h"
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
using math::DualQuaternion;
using math::Quaternion;
using math::Vector3;

namespace {

int tests_passed = 0;
int tests_failed = 0;

constexpr std::uint64_t kSeed = 0xC0FFEEull;
constexpr int kSamples = 8;

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

Vector3<T> random_small_v(std::mt19937_64& rng, T scale) {
    std::normal_distribution<T> dist(T(0), scale);
    return Vector3<T>(tv(dist(rng)), tv(dist(rng)), tv(dist(rng)));
}

DualQuaternion<T> random_unit_dq(std::mt19937_64& rng) {
    Quaternion<T> q_r = random_unit_q(rng);
    Vector3<T> t = random_v(rng);
    return DualQuaternion<T>::from_pose(q_r, t);
}

T q_diff(const Quaternion<T>& a, const Quaternion<T>& b) {
    return std::max({std::abs(a.w.value - b.w.value),
                     std::abs(a.x.value - b.x.value),
                     std::abs(a.y.value - b.y.value),
                     std::abs(a.z.value - b.z.value)});
}

T v_diff(const Vector3<T>& a, const Vector3<T>& b) {
    return std::max({std::abs(a.x.value - b.x.value),
                     std::abs(a.y.value - b.y.value),
                     std::abs(a.z.value - b.z.value)});
}

T dq_diff(const DualQuaternion<T>& a, const DualQuaternion<T>& b) {
    return std::max(q_diff(a.real, b.real), q_diff(a.dual, b.dual));
}

T q_pair_bound(const Quaternion<T>& a, const Quaternion<T>& b) {
    return std::max({a.w.total_error() + b.w.total_error(),
                     a.x.total_error() + b.x.total_error(),
                     a.y.total_error() + b.y.total_error(),
                     a.z.total_error() + b.z.total_error()});
}

T v_pair_bound(const Vector3<T>& a, const Vector3<T>& b) {
    return std::max({a.x.total_error() + b.x.total_error(),
                     a.y.total_error() + b.y.total_error(),
                     a.z.total_error() + b.z.total_error()});
}

T dq_pair_bound(const DualQuaternion<T>& a, const DualQuaternion<T>& b) {
    return std::max(q_pair_bound(a.real, b.real),
                    q_pair_bound(a.dual, b.dual));
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

void test_dq_associativity() {
    std::cout << "--- AUD-MC-15 (REQ-DQ-15): DQ composition associativity ---\n";
    std::mt19937_64 rng(kSeed);
    for (int i = 0; i < kSamples; ++i) {
        DualQuaternion<T> A = random_unit_dq(rng);
        DualQuaternion<T> B = random_unit_dq(rng);
        DualQuaternion<T> C = random_unit_dq(rng);
        DualQuaternion<T> lhs = (A * B) * C;
        DualQuaternion<T> rhs = A * (B * C);
        check_bound("(M1 M2) M3 == M1 (M2 M3)",
                    dq_diff(lhs, rhs), dq_pair_bound(lhs, rhs));
    }
}

void test_conjugate_dualities() {
    std::cout << "--- AUD-MC-16 (REQ-DQ-16): conjugate dualities ---\n";
    std::mt19937_64 rng(kSeed + 1);
    DualQuaternion<T> id = DualQuaternion<T>::identity();
    for (int i = 0; i < kSamples; ++i) {
        DualQuaternion<T> M = random_unit_dq(rng);

        // Quaternion conjugate is involutive.
        DualQuaternion<T> Mcc = M.conjugate().conjugate();
        check_bound("(M*)* == M",
                    dq_diff(Mcc, M), dq_pair_bound(Mcc, M));

        // Dual-number conjugate is involutive.
        DualQuaternion<T> Mdd = M.dual_conjugate().dual_conjugate();
        check_bound("(M_eps)_eps == M",
                    dq_diff(Mdd, M), dq_pair_bound(Mdd, M));

        // Combined conjugate equals composition of the other two.
        DualQuaternion<T> Msharp = M.combined_conjugate();
        DualQuaternion<T> Msharp_check = M.conjugate().dual_conjugate();
        check_bound("M♯ == (M*)_eps",
                    dq_diff(Msharp, Msharp_check),
                    dq_pair_bound(Msharp, Msharp_check));

        // M · M* == identity for unit DQ (REQ-DQ-15 corollary).
        DualQuaternion<T> MMstar = M * M.conjugate();
        check_bound("M * M* == identity (unit pose)",
                    dq_diff(MMstar, id),
                    dq_pair_bound(MMstar, id));
    }
}

void test_pose_action_composition() {
    std::cout << "--- AUD-MC-17 (REQ-DQ-17): pose action composition ---\n";
    std::mt19937_64 rng(kSeed + 2);
    for (int i = 0; i < kSamples; ++i) {
        DualQuaternion<T> M1 = random_unit_dq(rng);
        DualQuaternion<T> M2 = random_unit_dq(rng);
        Vector3<T> v  = random_v(rng);
        Vector3<T> lhs = (M1 * M2).apply(v);
        Vector3<T> rhs = M1.apply(M2.apply(v));
        check_bound("(M1 M2)(v) == M1(M2(v))",
                    v_diff(lhs, rhs), v_pair_bound(lhs, rhs));
    }
}

void test_screw_exp_log_roundtrip() {
    std::cout << "--- AUD-MC-18 (REQ-DQ-18): screw exp/log round-trip ---\n";
    std::mt19937_64 rng(kSeed + 3);
    // |u| must satisfy |u| ≤ π/2 (rotation angle ≤ π) to stay inside the
    // shortest-path recoverable branch of `log_unit` (Quaternion). A
    // per-component normal stddev of 0.2 gives E[|u|] ≈ 0.35 with negligible
    // probability of exceeding π/2 ≈ 1.57. v has no such constraint, so a
    // wider stddev is fine.
    for (int i = 0; i < kSamples; ++i) {
        Vector3<T> u = random_small_v(rng, T(0.2));
        Vector3<T> v = random_small_v(rng, T(0.5));

        DualQuaternion<T> M = DualQuaternion<T>::exp_screw(u, v);
        DualQuaternion<T> logM = M.log_screw();
        Vector3<T> u_rt = logM.angular();
        Vector3<T> v_rt = logM.linear();

        check_bound("log_screw(exp_screw(u, v)).angular == u",
                    v_diff(u_rt, u), v_pair_bound(u_rt, u));
        check_bound("log_screw(exp_screw(u, v)).linear == v",
                    v_diff(v_rt, v), v_pair_bound(v_rt, v));
    }
}

void test_from_pose_translation_roundtrip() {
    std::cout << "--- from_pose / translation accessor round-trip ---\n";
    std::mt19937_64 rng(kSeed + 4);
    for (int i = 0; i < kSamples; ++i) {
        Quaternion<T> q_r = random_unit_q(rng);
        Vector3<T> t   = random_v(rng);
        DualQuaternion<T> M   = DualQuaternion<T>::from_pose(q_r, t);
        Vector3<T> t_rt = M.translation();
        check_bound("translation(from_pose(q, t)) == t",
                    v_diff(t_rt, t), v_pair_bound(t_rt, t));
    }
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(17);
    std::cout << "test_dual_quaternion: REQ-DQ-15..18 / AUD-MC-15..18\n\n";

    test_dq_associativity();
    test_conjugate_dualities();
    test_pose_action_composition();
    test_screw_exp_log_roundtrip();
    test_from_pose_translation_roundtrip();

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << tests_passed
              << "  Failed: " << tests_failed << "\n";
    std::cout << "========================================\n";
    return tests_failed > 0 ? 1 : 0;
}
