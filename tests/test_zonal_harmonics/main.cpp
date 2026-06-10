/// test_zonal_harmonics — first slice of the common generative astrophysics
/// library (the home for the Jₖ; resolves the old D1 "where do odd zonals live"
/// question and the R04 provenance concern).
///
/// Verifies constants::ZonalHarmonics:
///   - precomputed once at construction, standard-tagged, cheap lookup;
///   - even zonals derived from the level ellipsoid;
///   - odd zonals generated from the born-digital normalized coefficient via
///     Jₙ = -√(2n+1)·C̄ₙ₀ (computed through the tracked sqrt);
///   - R5/O2 (closes R04): the C̄ₙ₀ are the FULL-PRECISION in-repo EGM2008
///     file values (15 digits, gen_egm_zonals.py), encoded `measured` with the
///     IERS Table 6.2 formal-error grade — the identity Jₙ = -√(2n+1)·C̄ₙ₀ is
///     re-checked against the file values, the budget carries the genuine σ
///     (measurement channel) instead of a digit floor, and the two stores
///     (constants::ZonalHarmonics::egm2008 and the ModelSelector
///     wgs84_precise preset) agree cross-site.
///
/// ExeGate: nonzero exit on any failed check.

#include "constants/zonal_harmonics.h"
#include "sgp4/model_selector.h"

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
double j3_precision() {
    return static_cast<double>(
        constants::ZonalHarmonics<T>::egm2008(T(1e-12)).Jn(3).errors.precision);
}

template<typename T>
double j3_measurement() {
    return static_cast<double>(
        constants::ZonalHarmonics<T>::egm2008(T(1e-12)).Jn(3).errors.measurement);
}

} // namespace

int main() {
    using T = double;
    constants::ZonalHarmonics<T> egm = constants::ZonalHarmonics<T>::egm2008(T(1e-12));

    check("standard tag EGM2008/IERS2010", egm.standard() == "EGM2008/IERS2010");
    check("max degree 9", egm.max_degree() == 9);

    // J2: even, derived from the ellipsoid (~1.0826e-3).
    check("J2 from ellipsoid (~1.0826e-3)", close_d(egm.Jn(2).value, 1.0826e-3, 1e-6));

    // O2 identity (closes R04): Jₙ == -√(2n+1)·C̄ₙ₀ with the C̄ₙ₀ VERBATIM from
    // the in-repo EGM2008 file (datalib/EGM-08norm100.txt, full precision).
    const double kCbarFile[4][2] = {{3, 0.957161207093473e-6},
                                    {5, 0.686702913736681e-7},
                                    {7, 0.905120844521618e-7},
                                    {9, 0.280180753216300e-7}};
    for (const auto& row : kCbarFile) {
        const int n = static_cast<int>(row[0]);
        const double expect = -std::sqrt(2.0 * n + 1.0) * row[1];
        check("Jn == -sqrt(2n+1)*Cbar_n0(file) (O2 identity)",
              close_d(egm.Jn(n).value, expect, std::abs(expect) * 1e-14));
    }

    // O2 cross-site: the ModelSelector modern preset (gravity "wgs84_precise",
    // which stores the generator-derived Jₙ directly) must agree with the
    // generative store — double-entry bookkeeping on the same file source.
    {
        sgp4::ModelConfiguration<T> cfg =
            sgp4::ModelSelector<T>::select("modern_2020", T(1e-12));
        bool cross_ok = true;
        for (int n = 3; n <= 9; n += 2) {
            const double a = cfg.Jn(n).value;
            const double b = egm.Jn(n).value;
            if (std::abs(a - b) > std::abs(b) * 1e-12) cross_ok = false;
        }
        check("ModelSelector wgs84_precise Jn == ZonalHarmonics egm2008 Jn (cross-site)",
              cross_ok);
    }

    // J4: even, derived (nonzero).
    check("J4 even derived, nonzero", std::abs(egm.Jn(4).value) > 0.0);

    // Outside the stored range -> exact zero.
    check("J11 absent -> 0", egm.Jn(11).value == 0.0);

    // WGS72: J3 given directly by SR3 (COMMON/C1).
    constants::ZonalHarmonics<T> wgs72 = constants::ZonalHarmonics<T>::wgs72(T(1e-12));
    check("WGS72 standard tag", wgs72.standard() == "WGS72/SR3");
    check("WGS72 J3 == -2.53881e-6", close_d(wgs72.Jn(3).value, -0.253881e-5, 1e-12));

    // R5/O2 re-encode: C̄ₙ₀ is `measured` with the IERS Table 6.2 formal-error
    // grade (±0.49e-11), so the genuine field-model σ lives in the MEASUREMENT
    // channel (propagated through √(2n+1): ~1.3e-11 on J3) — replacing the prior
    // digit-floor accuracy, which would UNDER-claim at the file's 15 digits.
    // Binary storage remains a computational precision that tightens with T.
    const double pd = j3_precision<double>();
    const double p50 = j3_precision<cpp_bin_float_50>();
    check("J3 precision > 0", pd > 0.0 && p50 > 0.0);
    check("J3 precision tightens with wider T (computational)", p50 < pd * 1e-20);
    const double md = j3_measurement<double>();
    const double m50 = j3_measurement<cpp_bin_float_50>();
    check("J3 measurement carries the formal-error grade (~1.3e-11)",
          md > 1e-12 && md < 1e-10);
    check("J3 measurement T-independent (a physical sigma, not arithmetic)",
          std::abs(m50 - md) <= md * 1e-9);

    std::cout << "\n  zonal harmonics: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
