/// test_constant_representation — CR1 (DQSGP4 Completion Roadmap, issue CR1).
///
/// A physical/astronomical constant is not its truncated decimal. This pins the
/// encoding discipline (feedback_constants_generative_or_bounded):
///   - defined("...") stamps the BINARY representation bound, so a 6-digit source
///     appears good to ~50 digits in cpp_bin_float_50 — Phase B then propagates
///     that over-claim faithfully (the bug);
///   - from_truncated_decimal("...") floors precision at the source's decimal
///     truncation 0.5*ulp10, which is FIXED across T (the honest representation
///     bound — no further digits exist to recover);
///   - a generative form (2*pi to full T, exact integers, arcsec*pi/648000)
///     genuinely scales with T.
///
/// ExeGate CR1: nonzero exit on any failed check.

#include "math/angles.h"
#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>

namespace {

using boost::multiprecision::cpp_bin_float_50;

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

template<typename T>
double prec_defined(const char* s) {
    return static_cast<double>(math::TrackedValue<T>::defined(s).errors.precision);
}

template<typename T>
double prec_truncated(const char* s) {
    return static_cast<double>(math::TrackedValue<T>::from_truncated_decimal(s).errors.precision);
}

template<typename T>
double prec_two_pi() {
    return static_cast<double>(math::two_pi<T>().errors.precision);
}

// defined() over-claims: its precision is the binary bound, so it shrinks with
// a wider T even though "6.283185" specified only 6 fractional digits.
void test_defined_overclaims() {
    double d = prec_defined<double>("6.283185");
    double b = prec_defined<cpp_bin_float_50>("6.283185");
    check("defined() precision shrinks with wider T (documents the over-claim)", b < d * 1e-20);
}

// from_truncated_decimal floors precision at 0.5*1e-6 = 5e-7, fixed across T.
void test_truncated_is_honest() {
    double d = prec_truncated<double>("6.283185");
    double b = prec_truncated<cpp_bin_float_50>("6.283185");
    check("truncated precision == 0.5e-6 (double)", close_d(d, 0.5e-6, 1e-12));
    check("truncated precision fixed across T (honest)", close_d(b, d, 1e-12));
    int rd = math::TrackedValue<cpp_bin_float_50>::from_truncated_decimal("6.283185").reliable_digits();
    check("truncated reliable_digits <= 8, not ~50", rd <= 8);
}

// A genuinely generative constant (2*pi, pi to full T) scales with T.
void test_generative_scales() {
    double d = prec_two_pi<double>();
    double b = prec_two_pi<cpp_bin_float_50>();
    check("generative 2*pi precision shrinks with wider T", b < d * 1e-20);
}

// scientific-notation and integer-place decimal_ulp.
void test_decimal_ulp() {
    using TV = math::TrackedValue<double>;
    // "3.9860e14": last digit place = 14 - 4 = 1e10.
    check("decimal_ulp scientific", close_d(TV::decimal_ulp("3.9860e14"), 1e10, 1e-3));
    // "84381.406": 3 fractional digits -> 1e-3.
    check("decimal_ulp fractional", close_d(TV::decimal_ulp("84381.406"), 1e-3, 1e-15));
}

// Obliquity exemplar: IAU 2006 epsilon_0 = 84381.406" (adopted to the mas). The
// honest encoding is the arcsecond form * (pi/648000) — materially tighter than
// the 6-figure degree truncation 23.4393 deg the C1 stub used.
void test_obliquity_exemplar() {
    using T = cpp_bin_float_50;
    using TV = math::TrackedValue<T>;
    TV arcsec = TV::from_truncated_decimal("84381.406");
    TV eps_arcsec = arcsec * (math::pi<T>() / TV::exact_integer(648000));
    TV eps_degree = math::degrees_to_radians(TV::from_truncated_decimal("23.4393"));
    check("obliquity arcsec value ~ degree value",
          close_d(static_cast<double>(eps_arcsec.value),
                  static_cast<double>(eps_degree.value), 1e-4));
    check("arcsec obliquity precision tighter than degree obliquity",
          eps_arcsec.errors.precision < eps_degree.errors.precision);
}

} // namespace

int main() {
    test_defined_overclaims();
    test_truncated_is_honest();
    test_generative_scales();
    test_decimal_ulp();
    test_obliquity_exemplar();

    std::cout << "\n  constant representation: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
