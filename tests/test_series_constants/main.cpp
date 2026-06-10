/// test_series_constants — SC1 / SC2 (user directives 2026-06-05, 2026-06-06).
///
/// SC1: "Constants that are not true by definition need to have their series-based
/// precision and accuracy tracked." The obliquity of the ecliptic is the
/// exemplar: it is NOT a defined constant — it is the leading terms of the IAU
/// 2006 time polynomial ε_A(t). astronomy::obliquity_iau2006 generates it from
/// that series and tracks BOTH error channels through the evaluation:
///   - PRECISION  = representation/rounding -> TIGHTENS with a wider type T;
///   - ACCURACY   = the series-TRUNCATION bound (omitted higher-order terms)
///                  -> TIGHTENS as more terms n are kept.
///
/// SC2: "All constants need accuracy and precision tracked; each J_k value is the
/// result of a series truncation." The zonal harmonics are the case in point —
/// the even J_{2n} ARE generated from the equipotential-ellipsoid series
/// (Heiskanen & Moritz; geodesy::EquipotentialEllipsoid::J2n), so the stored J_2
/// is a truncated-series result that carries accuracy AND precision, never a
/// defined()=0 "exact-by-convention" value. Section J proves the J_2 generated
/// from the series recovers the adopted WGS72 value, carries both error channels,
/// and that its precision tightens with a wider T.
///
/// ExeGate SC1: nonzero exit code on any failed check.

#include "astronomy/obliquity.h"
#include "astronomy/earth_orbit.h"
#include "constants/zonal_harmonics.h"
#include "math/tracked_polynomial.h"
#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

using boost::multiprecision::cpp_bin_float_50;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

} // namespace

int main() {
    using T = double;
    const double arcsec_to_rad = boost::math::constants::pi<double>() / 648000.0;

    // TrackedPolynomial (the new math/ primitive obliquity/earth_orbit now consume,
    // Stage 1) — eval() must be BIT-IDENTICAL to a hand-rolled Horner + Σ|c_k||t|^k
    // tail in ALL THREE error channels (the behavior-preservation contract).
    {
        using TV = math::TrackedValue<T>;
        std::vector<TV> c{TV::model_coefficient("1.5"), TV::model_coefficient("-0.25"),
                          TV::model_coefficient("0.0007")};
        std::vector<T> mag{T(1.5), T(0.25), T(0.0007)};
        math::TrackedPolynomial<T> poly(c, mag);
        TV t = tv<T>(0.37);
        TV got = poly.eval(t, 2, math::ErrorChannel::accuracy, 1);  // keep 2 of 3 + one extra tail
        // Hand reference: Horner over the first 2 terms, then the omitted k=2 term
        // and one extra conservative |c_2|·|t|^3 term, deposited in accuracy.
        TV ref = c[1] * t + c[0];
        using std::abs;
        T ta = abs(t.value), tp = ta * ta;            // |t|^2
        T tail = mag[2] * tp; tp = tp * ta;           // k=2 omitted term, then |t|^3
        tail = tail + mag[2] * tp;                    // one extra conservative term
        ref = math::add_bound(ref, tail, math::ErrorChannel::accuracy);
        check("TrackedPolynomial::eval == hand Horner+tail (value + all 3 channels)",
              got.value == ref.value &&
              got.errors.measurement == ref.errors.measurement &&
              got.errors.precision == ref.errors.precision &&
              got.errors.accuracy == ref.errors.accuracy);
    }

    // Value at the J2000 reference (t = 0): the full polynomial collapses to the
    // leading term 84381.406" -> 0.40909280 rad.
    math::TrackedValue<T> eps0 = astronomy::obliquity_iau2006<T>(tv<T>(0.0));
    check("obliquity(t=0) == 84381.406 arcsec (0.40909280 rad)",
          std::abs(eps0.value - 84381.406 * arcsec_to_rad) < 1e-12);

    // Secular drift: at t = 0.21 cen (~2021) the obliquity has DECREASED by the
    // -46.836769"·t term (~9.84" ~ 4.77e-5 rad).
    math::TrackedValue<T> eps21 = astronomy::obliquity_iau2006<T>(tv<T>(0.21));
    check("obliquity decreases from J2000 (secular -46.84\"/cen drift)",
          eps21.value < eps0.value);
    double drift = (eps0.value - eps21.value);
    check("drift over 0.21 cen ~ 46.84\"*0.21 ~ 4.77e-5 rad",
          std::abs(drift - 46.836769 * 0.21 * arcsec_to_rad) < 1e-7);

    // ACCURACY is the SERIES-TRUNCATION bound and TIGHTENS as more terms are
    // kept. n=1 keeps only 84381.406" and omits the -46.84"·t secular term, so
    // its accuracy is ~4.77e-5 rad; n=2 includes the secular term and the
    // omitted tail collapses to ~1e-11, leaving only the leading-coefficient
    // adoption floor (~2.4e-9 rad) -- orders of magnitude smaller.
    double acc_n1 = astronomy::obliquity_iau2006<T>(tv<T>(0.21), 1).errors.accuracy;
    double acc_n2 = astronomy::obliquity_iau2006<T>(tv<T>(0.21), 2).errors.accuracy;
    double acc_n6 = astronomy::obliquity_iau2006<T>(tv<T>(0.21), 6).errors.accuracy;
    check("accuracy(n=1) ~ omitted secular term (~4.77e-5 rad)",
          acc_n1 > 4e-5 && acc_n1 < 5e-5);
    check("accuracy tightens as more series terms are kept (n=1 >> n=2)",
          acc_n1 > acc_n2 * 100.0);
    check("accuracy(n>=2) bounded by the leading-coefficient adoption floor",
          acc_n2 < 1e-7 && acc_n6 <= acc_n2 * (1.0 + 1e-9));
    check("series accuracy is tracked (> 0) at every order",
          acc_n1 > 0.0 && acc_n2 > 0.0 && acc_n6 > 0.0);

    // PRECISION (representation) TIGHTENS with a wider numeric type T — the
    // calling card, now reaching a series-defined astronomical constant.
    double pd = astronomy::obliquity_iau2006<double>(tv<double>(0.21)).errors.precision;
    double pb = static_cast<double>(
        astronomy::obliquity_iau2006<cpp_bin_float_50>(tv<cpp_bin_float_50>(0.21)).errors.precision);
    check("obliquity precision > 0", pd > 0.0 && pb > 0.0);
    check("obliquity precision tightens with wider T", pb < pd * 1e-20);

    // --- Section E (SC2): Earth's orbital eccentricity as a secular series ----
    // The eccentricity twin of the obliquity: e(t) = Σ c_k t^k (Meeus 25.4 /
    // VSOP87). C1 pins its t=0 value; here we gate the SERIES behaviour the
    // directive demands — secular drift, truncation accuracy, T-scaling.
    math::TrackedValue<T> e0 = astronomy::earth_eccentricity<T>(tv<T>(0.0));
    check("eccentricity(t=0) == 0.016708634 (VSOP leading term)",
          std::abs(e0.value - 0.016708634) < 1e-12);
    math::TrackedValue<T> e21 = astronomy::earth_eccentricity<T>(tv<T>(0.21));
    check("eccentricity decreases from J2000 (secular -4.2e-5/cen drift)",
          e21.value < e0.value);
    double edrift = (e0.value - e21.value);
    check("ecc drift over 0.21 cen ~ 4.2037e-5*0.21 ~ 8.83e-6",
          std::abs(edrift - 0.000042037 * 0.21) < 1e-7);
    double eacc_n1 = astronomy::earth_eccentricity<T>(tv<T>(0.21), 1).errors.accuracy;
    double eacc_n2 = astronomy::earth_eccentricity<T>(tv<T>(0.21), 2).errors.accuracy;
    check("ecc accuracy(n=1) ~ omitted secular term (~8.8e-6)",
          eacc_n1 > 8e-6 && eacc_n1 < 1e-5);
    check("ecc accuracy tightens as more series terms are kept (n=1 >> n=2)",
          eacc_n1 > eacc_n2 * 100.0);
    double epd = astronomy::earth_eccentricity<double>(tv<double>(0.21)).errors.precision;
    double epb = static_cast<double>(
        astronomy::earth_eccentricity<cpp_bin_float_50>(tv<cpp_bin_float_50>(0.21)).errors.precision);
    check("ecc precision > 0", epd > 0.0 && epb > 0.0);
    check("ecc precision tightens with wider T", epb < epd * 1e-20);

    // --- Section J (SC2): each J_k is the result of a series truncation -------
    // The even zonal harmonics are GENERATED from the equipotential-ellipsoid
    // series J_{2n} = (-1)^{n+1} 3e^{2n}/((2n+1)(2n+3)) (1 - n + 5n J2/e^2), with
    // J2 = (e^2/3)(1 - 2 m e'/(15 q0)) and q0 itself an alternating series.
    // ZonalHarmonics realizes this: Jn(2) = earth.J2n(1).
    {
        constants::ZonalHarmonics<T> zh = constants::ZonalHarmonics<T>::wgs72(T(1e-15));
        math::TrackedValue<T> J2 = zh.Jn(2);   // generated from the ellipsoid series

        // Round-trip: from_J2 inverts the SAME Heiskanen-Moritz relation that
        // J2n(1) evaluates (b = a(1-f) = a*sqrt(1-e^2)), so the generated J_2
        // recovers the adopted WGS72 input 0.001082616.
        check("generated J_2 (ellipsoid series) recovers WGS72 0.001082616",
              std::abs(J2.value - 0.001082616) < 1e-9);

        // It is a TRUNCATED SERIES, not exact-by-definition: BOTH error channels
        // are tracked and nonzero (this is what defined()=0 accuracy denied).
        check("generated J_2 carries series accuracy (> 0)", J2.errors.accuracy > T(0));
        check("generated J_2 carries precision (> 0)", J2.errors.precision > T(0));

        // The calling card, now on a zonal harmonic: precision tightens with a
        // wider T (more series terms + finer representation).
        double pj_d = static_cast<double>(
            constants::ZonalHarmonics<double>::wgs72(1e-15).Jn(2).errors.precision);
        double pj_b = static_cast<double>(
            constants::ZonalHarmonics<cpp_bin_float_50>::wgs72(cpp_bin_float_50("1e-30"))
                .Jn(2).errors.precision);
        check("J_2 precision > 0 at double", pj_d > 0.0);
        check("J_2 precision tightens with wider T", pj_b < pj_d * 1e-6);

        // J_4 likewise generated from the series (the NORMAL field value, which
        // differs from the real Earth's empirical J_4 by the anomalous part — a
        // documented physical fact; here it carries its truncation accuracy).
        math::TrackedValue<T> J4 = zh.Jn(4);
        check("generated J_4 (normal field) is negative and tracked",
              J4.value < T(0) && J4.errors.accuracy > T(0));
    }

    std::cout << "\n  series constants: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
