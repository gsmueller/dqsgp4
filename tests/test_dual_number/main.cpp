/// test_dual_number — Algebraic correctness for math::DualNumber<T>.
///
/// Implements design/audit/mathematical_correctness.md AUD-MC-1..3,
/// verifying REQ-DQ-1..3 from
/// design/specifications/dual_quaternion_algebra.md:
///   AUD-MC-1: epsilon^2 = 0
///   AUD-MC-2: commutative-ring axioms
///   AUD-MC-3: forward-mode automatic differentiation
///
/// As in `test_quaternion`, each identity is also a check that the
/// reported `total_error()` bound dominates the actual deviation
/// (REQ-EF-2).

#include "math/dual_number.h"
#include "math/tracked_value.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>

using T = double;
using TV = math::TrackedValue<T>;
using math::DualNumber;

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

DualNumber<T> random_dn(std::mt19937_64& rng) {
    std::normal_distribution<T> dist(T(0), T(1));
    return DualNumber<T>(tv(dist(rng)), tv(dist(rng)));
}

T dn_diff(const DualNumber<T>& a, const DualNumber<T>& b) {
    return std::max(std::abs(a.real.value - b.real.value),
                    std::abs(a.dual.value - b.dual.value));
}

T dn_pair_bound(const DualNumber<T>& a, const DualNumber<T>& b) {
    return std::max(a.real.total_error() + b.real.total_error(),
                    a.dual.total_error() + b.dual.total_error());
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

void test_epsilon_squared() {
    std::cout << "--- AUD-MC-1 (REQ-DQ-1): epsilon^2 = 0 ---\n";
    DualNumber<T> eps = DualNumber<T>::epsilon();
    DualNumber<T> eps_sq = eps * eps;
    DualNumber<T> zero = DualNumber<T>::zero();
    check_bound("epsilon * epsilon == zero",
                dn_diff(eps_sq, zero), dn_pair_bound(eps_sq, zero));
}

void test_ring_axioms() {
    std::cout << "--- AUD-MC-2 (REQ-DQ-2): commutative ring axioms ---\n";
    std::mt19937_64 rng(kSeed);
    DualNumber<T> zero = DualNumber<T>::zero();
    DualNumber<T> id = DualNumber<T>::identity();
    for (int i = 0; i < kSamples; ++i) {
        DualNumber<T> a = random_dn(rng);
        DualNumber<T> b = random_dn(rng);
        DualNumber<T> c = random_dn(rng);

        // Additive associativity, commutativity, identity, inverse.
        DualNumber<T> l1 = (a + b) + c;
        DualNumber<T> r1 = a + (b + c);
        check_bound("(a+b)+c == a+(b+c)",
                    dn_diff(l1, r1), dn_pair_bound(l1, r1));

        DualNumber<T> l2 = a + b;
        DualNumber<T> r2 = b + a;
        check_bound("a + b == b + a",
                    dn_diff(l2, r2), dn_pair_bound(l2, r2));

        DualNumber<T> l3 = a + zero;
        check_bound("a + zero == a",
                    dn_diff(l3, a), dn_pair_bound(l3, a));

        DualNumber<T> l4 = a + (-a);
        check_bound("a + (-a) == zero",
                    dn_diff(l4, zero), dn_pair_bound(l4, zero));

        // Multiplicative associativity, commutativity, identity.
        DualNumber<T> l5 = (a * b) * c;
        DualNumber<T> r5 = a * (b * c);
        check_bound("(ab)c == a(bc)",
                    dn_diff(l5, r5), dn_pair_bound(l5, r5));

        DualNumber<T> l6 = a * b;
        DualNumber<T> r6 = b * a;
        check_bound("ab == ba",
                    dn_diff(l6, r6), dn_pair_bound(l6, r6));

        DualNumber<T> l7 = a * id;
        check_bound("a * identity == a",
                    dn_diff(l7, a), dn_pair_bound(l7, a));

        // Distributivity.
        DualNumber<T> l8 = a * (b + c);
        DualNumber<T> r8 = a * b + a * c;
        check_bound("a(b+c) == ab + ac",
                    dn_diff(l8, r8), dn_pair_bound(l8, r8));
    }
}

void test_forward_ad() {
    std::cout << "--- AUD-MC-3 (REQ-DQ-3): forward-mode automatic differentiation ---\n";
    std::mt19937_64 rng(kSeed + 100);
    std::uniform_real_distribution<T> dist(T(0.1), T(5.0));
    for (int i = 0; i < kSamples; ++i) {
        T a_val = dist(rng);
        // a = a_val + ε · 1 → f(a).dual is f'(a_val).
        DualNumber<T> a = DualNumber<T>(tv(a_val), tv(1.0));

        // sqrt: f'(x) = 1 / (2 √x).
        DualNumber<T> s = sqrt(a);
        T sqrt_real_expected = std::sqrt(a_val);
        T sqrt_dual_expected = T(1) / (T(2) * std::sqrt(a_val));
        check_bound(
            "sqrt(a).real == sqrt(a_val)",
            std::abs(s.real.value - sqrt_real_expected),
            s.real.total_error() + TV::representation_bound(sqrt_real_expected));
        check_bound(
            "sqrt(a).dual == 1/(2 sqrt(a_val))",
            std::abs(s.dual.value - sqrt_dual_expected),
            s.dual.total_error() + TV::representation_bound(sqrt_dual_expected));

        // sin: f'(x) = cos x.
        DualNumber<T> si = sin(a);
        T sin_real_expected = std::sin(a_val);
        T sin_dual_expected = std::cos(a_val);
        check_bound(
            "sin(a).real == sin(a_val)",
            std::abs(si.real.value - sin_real_expected),
            si.real.total_error() + TV::representation_bound(sin_real_expected));
        check_bound(
            "sin(a).dual == cos(a_val)",
            std::abs(si.dual.value - sin_dual_expected),
            si.dual.total_error() + TV::representation_bound(sin_dual_expected));

        // cos: f'(x) = -sin x.
        DualNumber<T> co = cos(a);
        T cos_real_expected = std::cos(a_val);
        T cos_dual_expected = -std::sin(a_val);
        check_bound(
            "cos(a).real == cos(a_val)",
            std::abs(co.real.value - cos_real_expected),
            co.real.total_error() + TV::representation_bound(cos_real_expected));
        check_bound(
            "cos(a).dual == -sin(a_val)",
            std::abs(co.dual.value - cos_dual_expected),
            co.dual.total_error() + TV::representation_bound(cos_dual_expected));
    }
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(17);
    std::cout << "test_dual_number: REQ-DQ-1..3 / AUD-MC-1..3\n\n";

    test_epsilon_squared();
    test_ring_axioms();
    test_forward_ad();

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << tests_passed
              << "  Failed: " << tests_failed << "\n";
    std::cout << "========================================\n";
    return tests_failed > 0 ? 1 : 0;
}
