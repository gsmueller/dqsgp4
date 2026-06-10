/// test_tracked_pow_cbrt — B2 (DQSGP4 Completion Roadmap, issue register B2).
///
/// Verifies the new TrackedValue<T> cbrt() and pow(x, const T& p) primitives:
///   - correct values (matching std::cbrt / std::pow, including negative cbrt);
///   - rigorous per-category error bounds in cancellation-safe form
///       cbrt: d / (3*cbrt(|x|-d)^2)        (derivative at the worst base)
///       pow : |p| * worst^(p-1) * d         (worst toward/away from 0 by sign(p-1))
///   - precision budgets nonzero at representation scale that tighten with a
///     wider T (the B-phase precision goal);
///   - graceful degeneration when the error swamps the base.
///
/// ExeGate B2: returns a nonzero exit code on any failed check.

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

void test_cbrt() {
    using T = double;
    using TV = math::TrackedValue<T>;

    const T x = 27.0;
    const T meas = 1e-3;
    TV a(x, meas, TV::representation_bound(x), T(0));
    TV r = cbrt(a);

    check("cbrt value", close_d(r.value, 3.0, 1e-13));
    const T c = std::cbrt(x - meas);
    const T expect = meas / (T(3) * c * c);
    check("cbrt measurement bound", close_d(r.errors.measurement, expect, 1e-12));
    check("cbrt precision > 0", r.errors.precision > T(0));

    // cbrt is total over R: negative arguments are fine.
    TV n(T(-8), T(0), TV::representation_bound(T(8)), T(0));
    check("cbrt negative value", close_d(cbrt(n).value, -2.0, 1e-13));
}

void test_pow() {
    using T = double;
    using TV = math::TrackedValue<T>;

    const T x = 8.0;
    const T meas = 1e-3;
    const T p = T(2) / T(3);
    TV a(x, meas, TV::representation_bound(x), T(0));
    TV r = pow(a, p);

    check("pow value (8^(2/3) = 4)", close_d(r.value, 4.0, 1e-12));
    // p < 1 so the worst base is x - meas; bound = p * (x-meas)^(p-1) * meas.
    const T expect = p * std::pow(x - meas, p - T(1)) * meas;
    check("pow measurement bound", close_d(r.errors.measurement, expect, 1e-12));
    check("pow precision > 0", r.errors.precision > T(0));

    // degenerate: the error swamps the base (worst base <= 0).
    TV d(T(1), T(2), TV::representation_bound(T(1)), T(0));
    TV rd = pow(d, p);
    check("pow degenerate finite & large",
          std::isfinite(rd.errors.measurement) && rd.errors.measurement > T(0));
}

void test_pow_tracked_exponent() {
    using T = double;
    using TV = math::TrackedValue<T>;

    // x^y with BOTH base and exponent tracked: x=2, y=10 -> 1024.
    const T xv = 2.0, yv = 10.0;
    const T d = 1e-4;
    TV x(xv, d, TV::representation_bound(xv), T(0));
    TV y(yv, d, TV::representation_bound(yv), T(0));
    TV r = pow(x, y);

    check("pow(x,y) value (2^10 = 1024)", close_d(r.value, 1024.0, 1e-9));
    check("pow(x,y) precision > 0", r.errors.precision > T(0));
    check("pow(x,y) measurement > 0 (propagates both args)", r.errors.measurement > T(0));

    // The composed bound (exp o (* ) o log, each a verified worst-endpoint bound)
    // must be a TRUE upper bound over [x +/- d] x [y +/- d]. x^y is monotone in
    // each arg for x>1, so the extreme deviations are at the corners.
    T maxdev = T(0);
    for (int sx = -1; sx <= 1; sx += 2)
        for (int sy = -1; sy <= 1; sy += 2) {
            T dev = std::abs(std::pow(xv + sx * d, yv + sy * d) - r.value);
            if (dev > maxdev) maxdev = dev;
        }
    check("pow(x,y) measurement bound is a true upper bound", r.errors.measurement >= maxdev);

    // Consistency: with a zero-error exponent, the tracked-exponent overload
    // agrees on value with the constant-exponent overload.
    TV y0(yv, T(0), TV::representation_bound(yv), T(0));
    check("pow(x, tracked y) value == pow(x, const y)",
          close_d(pow(x, y0).value, pow(x, yv).value, 1e-9));

    // Degenerate base (error swamps it): log()'s guard keeps the result finite.
    TV xb(T(1), T(2), TV::representation_bound(T(1)), T(0));
    TV rb = pow(xb, y);
    check("pow(x,y) degenerate base stays finite",
          std::isfinite(rb.value) && std::isfinite(static_cast<double>(rb.errors.measurement)));
}

template<typename T>
double powyy_precision_at(double x_d, double y_d) {
    using TV = math::TrackedValue<T>;
    TV x(T(x_d), T(0), TV::representation_bound(T(x_d)), T(0));
    TV y(T(y_d), T(0), TV::representation_bound(T(y_d)), T(0));
    return static_cast<double>(pow(x, y).errors.precision);
}

template<typename T>
double cbrt_precision_at(double x_d) {
    using TV = math::TrackedValue<T>;
    T x = T(x_d);
    TV a(x, T(0), TV::representation_bound(x), T(0));
    return static_cast<double>(cbrt(a).errors.precision);
}

void test_precision_scales_with_T() {
    using boost::multiprecision::cpp_bin_float_50;
    const double p_double = cbrt_precision_at<double>(27.0);
    const double p_bf50   = cbrt_precision_at<cpp_bin_float_50>(27.0);
    check("cbrt precision > 0 (double)", p_double > 0.0);
    check("cbrt precision > 0 (cpp_bin_float_50)", p_bf50 > 0.0);
    check("cbrt precision tightens with wider T", p_bf50 < p_double * 1e-20);

    // Same for the tracked-exponent pow(x, y) = exp(y*log x).
    const double pp_double = powyy_precision_at<double>(2.0, 10.0);
    const double pp_bf50   = powyy_precision_at<cpp_bin_float_50>(2.0, 10.0);
    check("pow(x,y) precision > 0 (double)", pp_double > 0.0);
    check("pow(x,y) precision tightens with wider T", pp_bf50 < pp_double * 1e-20);
}

} // namespace

int main() {
    test_cbrt();
    test_pow();
    test_pow_tracked_exponent();
    test_precision_scales_with_T();

    std::cout << "\n  tracked pow/cbrt: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
