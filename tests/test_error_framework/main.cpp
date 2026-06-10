/// test_error_framework — property-based verification of the three-error
/// framework (AUD-EF-1, -2, -3, -4, -5, -7, -9, -10 / REQ-EF-1..9, -12).
///
/// The load-bearing audit item is AUD-EF-10 / REQ-EF-2: the reported
/// `total_error()` must be a RIGOROUS UPPER BOUND on the true deviation. This
/// test verifies that empirically for the closed-form scalar operators and for
/// the composite types, against a `cpp_bin_float_50` reference.
///
/// Method (per AUD-EF-10). For each operation f and many random trials:
///   1. sample a double input x at a known measurement error δ;
///   2. compute the tracked result y = f(TrackedValue(x, δ)) in double;
///   3. place the "true" input within the error interval, x_true = x + ε with
///      |ε| ≤ 0.9 δ, and compute f(x_true) in 50-digit precision (the model's
///      infinite-precision value, per REQ-EF-2);
///   4. assert |y.value − f(x_true)| ≤ y.total_error().
///
/// We validate error PROPAGATION (δ chosen ≫ machine epsilon) rather than the
/// bare last-ulp rounding of libm transcendentals, which the framework's bounds
/// do not separately model; δ dominates, so the test is rigorous for the
/// framework's actual guarantee and free of floating-point-boundary flakiness.
/// ε is drawn strictly inside the guaranteed domain (|ε| ≤ 0.9 δ) so a correct
/// bound always holds with margin while a too-small bound still fails.
///
/// AUD-EF-7: the same check is applied to composite types (Vector3, DualNumber),
/// confirming they compose the scalar `TrackedValue` propagation rather than
/// duplicating or dropping it. AUD-EF-9: the catastrophic branches (sqrt of a
/// value below its own error, division by a near-zero quantity) are checked to
/// return a bound that dominates the result. AUD-EF-1: every operation returns
/// a tracked value carrying a populated, nonzero error budget.
///
/// Deterministic (fixed RNG seed, per OBJ-8). Exit 0 iff every trial holds.

#include "math/dual_number.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using T = double;
using TV = math::TrackedValue<T>;
using HP = boost::multiprecision::cpp_bin_float_50;
using math::DualNumber;
using math::Vector3;

namespace {

int checks_passed = 0;
int checks_failed = 0;
std::mt19937_64 rng(0xE2202624ull);  // fixed seed: deterministic (OBJ-8)

double runif(double lo, double hi) {
    std::uniform_real_distribution<double> d(lo, hi);
    return d(rng);
}

/// Tracked input at value x with measurement error δ (relative, with an
/// absolute floor) plus the unavoidable representation precision.
TV tracked_input(double x, double delta) {
    return TV(x, delta, TV::representation_bound(x), T(0));
}
double input_delta(double x) { return std::max(1e-7 * std::abs(x), 1e-12); }

/// Report a property block: pass iff zero violations; prints the tightest
/// margin seen (max actual/bound ratio) so the bounds' slack is visible.
void report(const std::string& name, int violations, double worst_ratio) {
    if (violations == 0) {
        ++checks_passed;
        std::cout << "  PASS: " << name << "  (worst actual/bound = "
                  << std::setprecision(3) << worst_ratio << ")\n";
    } else {
        ++checks_failed;
        std::cerr << "  FAIL: " << name << "  (" << violations
                  << " trials with actual > bound)\n";
    }
}

constexpr int kTrials = 300;

/// AUD-EF-10 property check for a unary operator.
void check_unary(const std::string& name,
                 std::function<TV(const TV&)> fd,
                 std::function<HP(const HP&)> fh,
                 double lo, double hi) {
    int viol = 0;
    double worst = 0.0;
    for (int t = 0; t < kTrials; ++t) {
        double x = runif(lo, hi);
        double dx = input_delta(x);
        TV y = fd(tracked_input(x, dx));
        double ex = runif(-0.9 * dx, 0.9 * dx);
        HP truth = fh(HP(x) + HP(ex));
        HP actual = abs(HP(y.value) - truth);
        HP bound = HP(y.total_error());
        if (bound > 0) worst = std::max(worst, static_cast<double>(actual / bound));
        if (actual > bound) ++viol;
    }
    report(name, viol, worst);
}

/// AUD-EF-10 property check for a binary operator.
void check_binary(const std::string& name,
                  std::function<TV(const TV&, const TV&)> fd,
                  std::function<HP(const HP&, const HP&)> fh,
                  double loa, double hia, double lob, double hib) {
    int viol = 0;
    double worst = 0.0;
    for (int t = 0; t < kTrials; ++t) {
        double a = runif(loa, hia), b = runif(lob, hib);
        double da = input_delta(a), db = input_delta(b);
        TV y = fd(tracked_input(a, da), tracked_input(b, db));
        double ea = runif(-0.9 * da, 0.9 * da), eb = runif(-0.9 * db, 0.9 * db);
        HP truth = fh(HP(a) + HP(ea), HP(b) + HP(eb));
        HP actual = abs(HP(y.value) - truth);
        HP bound = HP(y.total_error());
        if (bound > 0) worst = std::max(worst, static_cast<double>(actual / bound));
        if (actual > bound) ++viol;
    }
    report(name, viol, worst);
}

} // anonymous namespace

int main() {
    using std::sqrt; using std::sin; using std::cos; using std::atan;
    using std::atan2; using std::abs;
    std::cout << "test_error_framework: AUD-EF property tests vs cpp_bin_float_50\n\n";

    std::cout << "=== AUD-EF-10 / REQ-EF-2: total_error() bounds the true error ===\n";
    // Closed-form scalar operators (AUD-EF-2: each carries a per-category bound).
    check_binary("REQ-EF-3 operator+  bound >= actual",
                 [](const TV& a, const TV& b) { return a + b; },
                 [](const HP& a, const HP& b) { return a + b; }, -10, 10, -10, 10);
    check_binary("REQ-EF-3 operator-  bound >= actual",
                 [](const TV& a, const TV& b) { return a - b; },
                 [](const HP& a, const HP& b) { return a - b; }, -10, 10, -10, 10);
    check_binary("REQ-EF-3 operator*  bound >= actual",
                 [](const TV& a, const TV& b) { return a * b; },
                 [](const HP& a, const HP& b) { return a * b; }, -10, 10, -10, 10);
    check_binary("REQ-EF-3 operator/  bound >= actual",
                 [](const TV& a, const TV& b) { return a / b; },
                 [](const HP& a, const HP& b) { return a / b; }, -10, 10, 0.5, 10);
    check_unary("REQ-EF-3 sqrt       bound >= actual",
                [](const TV& a) { return sqrt(a); },
                [](const HP& a) { using std::sqrt; return sqrt(a); }, 0.1, 10);
    check_unary("REQ-EF-3 sin        bound >= actual",
                [](const TV& a) { return sin(a); },
                [](const HP& a) { using std::sin; return sin(a); }, -3, 3);
    check_unary("REQ-EF-3 cos        bound >= actual",
                [](const TV& a) { return cos(a); },
                [](const HP& a) { using std::cos; return cos(a); }, -3, 3);
    check_unary("REQ-EF-3 atan       bound >= actual",
                [](const TV& a) { return atan(a); },
                [](const HP& a) { using std::atan; return atan(a); }, -5, 5);
    check_binary("REQ-EF-3 atan2      bound >= actual",
                 [](const TV& y, const TV& x) { return atan2(y, x); },
                 [](const HP& y, const HP& x) { using std::atan2; return atan2(y, x); },
                 1, 5, 1, 5);

    std::cout << "\n=== AUD-EF-7 / REQ-EF-12: composites compose the scalar budget ===\n";
    // Vector3::magnitude — √(x²+y²+z²) through TrackedValue arithmetic.
    {
        int viol = 0; double worst = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            double x = runif(-10, 10), y = runif(-10, 10), z = runif(-10, 10);
            double dx = input_delta(x), dy = input_delta(y), dz = input_delta(z);
            Vector3<T> v(tracked_input(x, dx), tracked_input(y, dy), tracked_input(z, dz));
            TV m = v.magnitude();
            double ex = runif(-0.9*dx, 0.9*dx), ey = runif(-0.9*dy, 0.9*dy),
                   ez = runif(-0.9*dz, 0.9*dz);
            HP X = HP(x)+HP(ex), Y = HP(y)+HP(ey), Z = HP(z)+HP(ez);
            HP truth = sqrt(X*X + Y*Y + Z*Z);
            HP actual = abs(HP(m.value) - truth);
            HP bound = HP(m.total_error());
            if (bound > 0) worst = std::max(worst, static_cast<double>(actual / bound));
            if (actual > bound) ++viol;
        }
        report("REQ-EF-12 Vector3::magnitude bound >= actual", viol, worst);
    }
    // DualNumber multiplication — exercises the dual-arithmetic composition.
    {
        int viol = 0; double worst = 0.0;
        for (int t = 0; t < kTrials; ++t) {
            double ar = runif(-5,5), ad = runif(-5,5), br = runif(-5,5), bd = runif(-5,5);
            double dar=input_delta(ar), dad=input_delta(ad), dbr=input_delta(br), dbd=input_delta(bd);
            DualNumber<T> a(tracked_input(ar,dar), tracked_input(ad,dad));
            DualNumber<T> b(tracked_input(br,dbr), tracked_input(bd,dbd));
            DualNumber<T> p = a * b;
            double er=runif(-0.9*dar,0.9*dar), ed=runif(-0.9*dad,0.9*dad),
                   fr=runif(-0.9*dbr,0.9*dbr), fd2=runif(-0.9*dbd,0.9*dbd);
            HP AR=HP(ar)+HP(er), AD=HP(ad)+HP(ed), BR=HP(br)+HP(fr), BD=HP(bd)+HP(fd2);
            HP truth_real = AR*BR;                 // (real part of a*b)
            HP truth_dual = AR*BD + AD*BR;         // (dual part of a*b)
            HP actual = abs(HP(p.real.value) - truth_real) + abs(HP(p.dual.value) - truth_dual);
            HP bound = HP(p.real.total_error()) + HP(p.dual.total_error());
            if (bound > 0) worst = std::max(worst, static_cast<double>(actual / bound));
            if (actual > bound) ++viol;
        }
        report("REQ-EF-12 DualNumber operator* bound >= actual", viol, worst);
    }

    std::cout << "\n=== AUD-EF-9 / REQ-EF-9: catastrophic regions return a dominating bound ===\n";
    {
        // sqrt of a value smaller than its own error: result is fully uncertain.
        TV tiny(1e-20, T(0), T(0), T(0));
        tiny.errors.precision = 1e-10;   // error >> value
        TV r = sqrt(tiny);
        bool ok = r.total_error() >= std::abs(r.value) && r.total_error() > 0;
        if (ok) { ++checks_passed; std::cout << "  PASS: AUD-EF-9 sqrt(value<error) error dominates result"
                  << "  (total_error=" << r.total_error() << ")\n"; }
        else { ++checks_failed; std::cerr << "  FAIL: AUD-EF-9 sqrt(value<error)\n"; }
    }
    {
        // Division by a denominator smaller than its error: bound blows up.
        TV num(1.0, T(0), 1e-16, T(0));
        TV den(1e-20, T(0), T(0), T(0));
        den.errors.precision = 1e-10;    // |den| < error: catastrophic
        TV q = num / den;
        bool ok = q.total_error() >= std::abs(q.value);
        if (ok) { ++checks_passed; std::cout << "  PASS: AUD-EF-9 div by near-zero error dominates result\n"; }
        else { ++checks_failed; std::cerr << "  FAIL: AUD-EF-9 div by near-zero\n"; }
    }

    std::cout << "\n=== AUD-EF-1 / REQ-EF-1: operations return a populated error budget ===\n";
    {
        TV a = tracked_input(2.0, 1e-7), b = tracked_input(3.0, 1e-7);
        TV s = a + b, p = a * b, q = a / b, r = sqrt(a);
        bool alive = s.total_error() > 0 && p.total_error() > 0 &&
                     q.total_error() > 0 && r.total_error() > 0;
        if (alive) { ++checks_passed; std::cout << "  PASS: AUD-EF-1 +,*,/,sqrt all return nonzero total_error\n"; }
        else { ++checks_failed; std::cerr << "  FAIL: AUD-EF-1 an operation returned zero total_error\n"; }
    }

    std::cout << "\n=== add_bound: the error-channel chokepoint is value-neutral, channel-explicit ===\n";
    {
        // The single sanctioned splice for widening a budget after the fact must
        // (a) leave .value and the two non-chosen channels bit-identical, and
        // (b) move EXACTLY the chosen channel by EXACTLY the bound — so it is a
        // pure rewrite of the hand-typed `v.errors.<field> = ... + bound` lines.
        TV v(1.25, 0.1, 2e-3, 5e-4);   // value, measurement, precision, accuracy
        const T b = 7e-4;
        TV a = math::add_bound(v, b, math::ErrorChannel::accuracy);
        TV p = math::add_bound(v, b, math::ErrorChannel::precision);
        bool ok =
            a.value == v.value && p.value == v.value &&
            a.errors.measurement == v.errors.measurement &&
            p.errors.measurement == v.errors.measurement &&
            a.errors.precision == v.errors.precision &&                 // accuracy path: precision untouched
            a.errors.accuracy  == v.errors.accuracy + b &&              // accuracy moved by exactly the bound
            p.errors.accuracy  == v.errors.accuracy &&                  // precision path: accuracy untouched
            p.errors.precision == v.errors.precision + b;               // precision moved by exactly the bound
        if (ok) { ++checks_passed; std::cout << "  PASS: add_bound moves exactly the chosen channel, value bit-identical\n"; }
        else { ++checks_failed; std::cerr << "  FAIL: add_bound channel/value semantics\n"; }
    }

    std::cout << "\n========================================\n";
    std::cout << "Checks passed: " << checks_passed
              << "  failed: " << checks_failed << "\n";
    std::cout << "========================================\n";
    return checks_failed > 0 ? 1 : 0;
}
