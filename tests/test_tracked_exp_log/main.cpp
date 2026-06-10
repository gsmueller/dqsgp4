/// test_tracked_exp_log — B1 (DQSGP4 Completion Roadmap, issue register B1).
///
/// Verifies the new TrackedValue<T> exp() and log() primitives:
///   - compute the correct value (matching std::exp / std::log);
///   - propagate the three-error budget through the rigorous analytic bound
///       exp:  |exp(x+d) - exp(x)| = exp(x)*(exp(d) - 1)   (convex, +d worst)
///       log:  |log(x) - log(x-d)| <= d / (x - d)          (concave, -d worst)
///   - round-trip log(exp(x)) == x in value;
///   - degrade gracefully when log's error swamps its argument (d >= x);
///   - and, the point of B1-c: the precision budget TIGHTENS with a wider
///     numeric type T, so a higher-precision mode actually improves the tracked
///     drag-density precision that was previously dropped (drag.h raw-`exp`).
///
/// ExeGate B1: returns a nonzero exit code on any failed check.

#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool ok) {
    if (ok) {
        ++passed;
        std::cout << "  PASS: " << name << "\n";
    } else {
        ++failed;
        std::cout << "  FAIL: " << name << "\n";
    }
}

bool close_d(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

// --- exp: value + rigorous error propagation (double) ---------------------
void test_exp_double() {
    using T = double;
    using TV = math::TrackedValue<T>;

    const T x = 0.75;
    const T meas = 1e-3;
    TV a(x, meas, TV::representation_bound(x), T(0));
    TV r = exp(a);  // ADL -> math::exp

    check("exp value", close_d(r.value, std::exp(x), 1e-14));
    // measurement bound = exp(x)*(exp(meas) - 1)
    const T expect = std::exp(x) * (std::exp(meas) - T(1));
    check("exp measurement bound", close_d(r.errors.measurement, expect, 1e-12));
    check("exp precision > 0 (propagated from repr bound)", r.errors.precision > T(0));
    check("exp accuracy == 0 (input had none)", r.errors.accuracy == T(0));
}

// --- log: value + rigorous error propagation (double) ---------------------
void test_log_double() {
    using T = double;
    using TV = math::TrackedValue<T>;

    const T x = 2.5;
    const T meas = 1e-3;
    TV a(x, meas, TV::representation_bound(x), T(0));
    TV r = log(a);

    check("log value", close_d(r.value, std::log(x), 1e-14));
    // measurement bound = meas / (x - meas)
    const T expect = meas / (x - meas);
    check("log measurement bound", close_d(r.errors.measurement, expect, 1e-12));
    check("log precision > 0", r.errors.precision > T(0));
}

// --- round-trip: log(exp(x)) == x in value --------------------------------
void test_round_trip() {
    using T = double;
    using TV = math::TrackedValue<T>;
    const T x = 1.3;
    TV a(x, T(0), TV::representation_bound(x), T(0));
    TV r = log(exp(a));
    check("round-trip log(exp(x)) value", close_d(r.value, x, 1e-12));
}

// --- log degenerate guard: error swamps the argument (d >= x) -------------
void test_log_degenerate() {
    using T = double;
    using TV = math::TrackedValue<T>;
    // value 1.0 with a measurement error of 2.0 -> argument could be <= 0.
    TV a(T(1), T(2), TV::representation_bound(T(1)), T(0));
    TV r = log(a);
    check("log value finite", std::isfinite(r.value));
    check("log degenerate bound large & finite",
          r.errors.measurement > T(0) && std::isfinite(r.errors.measurement));
}

// --- B1-c: precision tightens with a wider T ------------------------------
template<typename T>
double exp_precision_at(double x_d) {
    using TV = math::TrackedValue<T>;
    T x = T(x_d);
    TV a(x, T(0), TV::representation_bound(x), T(0));
    return static_cast<double>(exp(a).errors.precision);
}

void test_precision_scales_with_T() {
    using boost::multiprecision::cpp_bin_float_50;
    const double p_double = exp_precision_at<double>(0.5);
    const double p_bf50   = exp_precision_at<cpp_bin_float_50>(0.5);

    // double's exp precision ~ exp(0.5)*2^-53 ~ 1.8e-16; cpp_bin_float_50's is
    // ~1e-50. The wider type must shrink the tracked precision by many orders.
    check("exp precision > 0 (double)", p_double > 0.0);
    check("exp precision > 0 (cpp_bin_float_50)", p_bf50 > 0.0);
    check("exp precision tightens with wider T (B1-c)", p_bf50 < p_double * 1e-20);
}

} // namespace

int main() {
    test_exp_double();
    test_log_double();
    test_round_trip();
    test_log_degenerate();
    test_precision_scales_with_T();

    std::cout << "\n  tracked exp/log: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
