/// test_ephemeris — EPH (DQSGP4 Completion Roadmap, register EPH).
///
/// Verifies the ephemeris fidelity work:
///   - SOLAR equation-of-center 2-term truncation is FOLDED into errors.accuracy
///     ((13/12)|e|³), and the distance approximation carries its O(e²) bound;
///   - LUNAR longitude/latitude now carry honest model-fidelity accuracy bounds
///     for the omitted Meeus periodics (~2.5° / ~0.9°) — the register's weakest
///     link — and the equation of center gained the (5/4)e² sin(2l) term;
///   - GMST gains a UT1−UTC correction: ΔUT1 = 0 reproduces compute_gmst EXACTLY
///     (the authentic SGP4 path is unaffected) and a nonzero ΔUT1 shifts GMST by
///     the expected Earth rotation;
///   - born-digital reference: Earth's max equation of center ≈ 1.915°;
///   - L3 MODERN solar instance (sun_meeus.h, Epoch-typed) reproduces an
///     INDEPENDENT astropy/ERFA oracle to ≤ 0.01° in λ and ≤ 1e-4 AU in R, and
///     its tracked accuracy MAJORIZES the measured residual (no perceived
///     fidelity); reference from tools/gen_solar_oracle.py;
///   - precision tightens with a wider numeric type T.
///
/// ExeGate EPH: nonzero exit code on any failed check.

#include "ephemeris/solar_ephemeris.h"
#include "ephemeris/lunar_ephemeris.h"
#include "ephemeris/celestial_body.h"
#include "ephemeris/sun_meeus.h"
#include "ephemeris/moon_meeus.h"
#include "orbit/kepler_series.h"
#include "astronomy/epoch.h"
#include "astronomy/sidereal_time.h"
#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_bin_float_50;

const double kPI = boost::math::constants::pi<double>();

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

// Clean input: measurement = accuracy = 0, only the binary-representation
// precision (which scales with T). Any accuracy in the output is therefore
// attributable to the ephemeris model; precision exercises the T-scaling.
template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

// A Sun with mean anomaly = M0 at delta_t=0 (zero rates), clean inputs.
template<typename T>
ephemeris::CelestialBody<T> make_test_sun(double e, double M0, double argp) {
    ephemeris::CelestialBody<T> s;
    s.name = "Sun";
    s.mean_anomaly_rate_rad_min = tv<T>(0.0);
    s.mean_anomaly_epoch = tv<T>(M0);
    s.longitude_rate_rad_day = tv<T>(0.0);
    s.longitude_epoch = tv<T>(0.0);
    s.orbit_eccentricity = tv<T>(e);
    s.orbit_arg_perigee = tv<T>(argp);
    s.orbit_inclination = math::exact<T>(0);
    s.node_epoch = math::exact<T>(0);
    s.node_rate_rad_day = math::exact<T>(0);
    return s;
}

template<typename T>
ephemeris::CelestialBody<T> make_test_moon(double e, double l0) {
    ephemeris::CelestialBody<T> m;
    m.name = "Moon";
    m.mean_anomaly_rate_rad_min = tv<T>(0.0);
    m.mean_anomaly_epoch = tv<T>(l0);
    m.longitude_rate_rad_day = tv<T>(0.0);
    m.longitude_epoch = tv<T>(0.0);
    m.orbit_eccentricity = tv<T>(e);
    m.orbit_arg_perigee = tv<T>(0.0);
    m.orbit_inclination = tv<T>(5.145 * boost::math::constants::pi<T>() / T(180));  // 5.145°
    m.node_epoch = tv<T>(0.0);
    m.node_rate_rad_day = tv<T>(0.0);
    return m;
}

void test_solar() {
    using T = double;
    const double e = 0.0167086;
    auto sun = make_test_sun<T>(e, boost::math::constants::pi<T>() / 2.0, 0.0);  // M=90°
    auto pos = ephemeris::compute_solar_position<T>(sun, tv<T>(0.0));

    // The 2-term EoC truncation (13/12)e³ is folded into the longitude accuracy.
    double eoc_trunc = (13.0 / 12.0) * std::pow(e, 3);
    check("solar longitude carries EoC truncation accuracy",
          std::abs(pos.ecliptic_longitude.errors.accuracy - eoc_trunc) < 1e-12);
    // Distance carries its O(e²) model bound.
    check("solar distance carries O(e^2) accuracy",
          std::abs(pos.distance_au.errors.accuracy - e * e) < 1e-12);

    // Born-digital reference: Earth's maximum equation of center is ~1.915°
    // (at M ≈ 90°, EoC ≈ 2e = 1.914° for e=0.0167). The true longitude minus
    // the mean anomaly is the equation of center here (argp = 0).
    double eoc_deg = (pos.ecliptic_longitude.value - boost::math::constants::pi<T>() / 2.0)
                     * 180.0 / boost::math::constants::pi<T>();
    check("solar EoC at quadrature ~ 1.915 deg (born-digital ref)",
          eoc_deg > 1.90 && eoc_deg < 1.93);
}

void test_lunar() {
    using T = double;
    const double e = 0.0549006;
    const double deg = boost::math::constants::pi<T>() / 180.0;

    // EoC fidelity: at l=45° the (5/4)e² sin(2l) term is at its max and must be
    // present (distinguishes the 2-term form from the old single-term 2e sin l).
    auto moon45 = make_test_moon<T>(e, boost::math::constants::pi<T>() / 4.0);
    auto p45 = ephemeris::compute_lunar_position<T>(moon45, tv<T>(0.0));
    double eoc = 2.0 * e * std::sin(kPI / 4) + 1.25 * e * e * std::sin(kPI / 2);
    check("lunar 2-term EoC value (includes 5/4 e^2 sin 2l)",
          std::abs(p45.ecliptic_longitude.value - eoc) < 1e-12);

    // Longitude carries the ~2.5° omitted-periodics model bound (weakest link).
    check("lunar longitude carries ~2.5deg model accuracy",
          std::abs(p45.ecliptic_longitude.errors.accuracy - 2.5 * deg) < 1e-9);
    // Latitude carries the ~0.9° omitted-periodics model bound.
    check("lunar latitude carries ~0.9deg model accuracy",
          std::abs(p45.ecliptic_latitude.errors.accuracy - 0.9 * deg) < 1e-9);
    // Sanity: the trig values inherit a nonzero accuracy from the angle bound.
    check("lunar sin/cos longitude carry the model accuracy",
          p45.sin_longitude.errors.accuracy > 0.0 && p45.cos_longitude.errors.accuracy > 0.0);
}

void test_ut1_utc() {
    using T = double;
    auto jd = tv<T>(2451545.3);  // J2000 + 0.3 day (away from midnight)

    // delta_ut1 = 0 reproduces compute_gmst exactly (SGP4 path unaffected).
    T g0 = astronomy::compute_gmst<T>(jd).value;
    T g0u = astronomy::compute_gmst_ut1<T>(jd, tv<T>(0.0)).value;
    check("UT1 correction with dUT1=0 == compute_gmst (bit-exact)", g0 == g0u);

    // A 0.5 s UT1−UTC shifts GMST by sidereal_ratio·(0.5/86400)·2π.
    T gShift = astronomy::compute_gmst_ut1<T>(jd, tv<T>(0.5)).value;
    double expect = 1.00273790935 * (0.5 / 86400.0) * 2.0 * kPI;
    double diff = gShift - g0;
    // reduce to (-pi, pi] in case of a 2pi wrap boundary
    while (diff > kPI) diff -= 2.0 * kPI;
    while (diff < -kPI) diff += 2.0 * kPI;
    // 1e-8 rad tolerance: above the ~1e-9 double-rounding noise from the large
    // JD (~2.45e6) in the frac_day -> ut1_seconds chain, yet tight enough to
    // catch a missing sidereal_ratio (~1e-7 off). The shift's own precision
    // budget tracks the JD-magnitude rounding separately.
    check("UT1 correction shifts GMST by the expected rotation",
          std::abs(diff - expect) < 1e-8);
}

void test_precision_scales() {
    const double e = 0.0167086;
    double pd = ephemeris::compute_solar_position<double>(
        make_test_sun<double>(e, 1.0, 0.0), tv<double>(0.0)).ecliptic_longitude.errors.precision;
    double pb = static_cast<double>(ephemeris::compute_solar_position<cpp_bin_float_50>(
        make_test_sun<cpp_bin_float_50>(e, 1.0, 0.0),
        tv<cpp_bin_float_50>(0.0)).ecliptic_longitude.errors.precision);
    check("solar longitude precision > 0", pd > 0.0 && pb > 0.0);
    check("solar longitude precision tightens with wider T", pb < pd * 1e-20);
}

// L3: the modern Epoch-typed solar instance (sun_meeus.h) vs the astropy/ERFA
// oracle — the "no perceived fidelity" gate that the SR3 instance lacked.
void test_solar_astropy_oracle() {
    using T = double;
    // Born-digital astropy/ERFA reference — Sun GEOMETRIC geocentric ecliptic-of-
    // date longitude [deg] and distance [AU] — from tools/gen_solar_oracle.py
    // (astropy 7.1.1, builtin ephemeris; geometric = apparent + κ/R aberration).
    struct Row { double jd, lon_deg, R_au; };
    const Row oracle[] = {
        {2451545.0, 280.377823251, 0.983327625},
        {2455000.5,  86.907575467, 1.016029635},
        {2459000.5,  70.023147400, 1.013903983},
        {2460310.5, 280.045456475, 0.983318284},
        {2462000.5, 144.623306677, 1.012437483},
    };
    const double rad2deg = 180.0 / kPI;
    double max_dlon = 0.0, max_dR = 0.0;
    bool majorized = true;
    for (const auto& r : oracle) {
        auto epoch = astronomy::Epoch<T>::from_jd(tv<T>(r.jd), astronomy::TimeScale::TT);
        auto s = ephemeris::sun_meeus_of_date<T>(epoch);
        double dlon = s.lon.value * rad2deg - r.lon_deg;
        while (dlon > 180.0) dlon -= 360.0;
        while (dlon < -180.0) dlon += 360.0;
        double dR = s.radius.value - r.R_au;
        max_dlon = std::max(max_dlon, std::abs(dlon));
        max_dR = std::max(max_dR, std::abs(dR));
        // The tracked accuracy must MAJORIZE the measured residual (the oracle
        // bound is honest, not asserted): |Δλ| ≤ acc(λ), |ΔR| ≤ acc(R).
        if (std::abs(dlon) > s.lon.errors.accuracy * rad2deg) majorized = false;
        if (std::abs(dR) > s.radius.errors.accuracy) majorized = false;
    }
    std::cout << "    sun_meeus vs astropy: max|Δλ| = " << max_dlon * 3600.0
              << "\"  max|ΔR| = " << max_dR << " AU\n";
    check("sun_meeus λ reproduces astropy oracle ≤ 0.01° (Meeus grade)", max_dlon < 0.01);
    check("sun_meeus R reproduces astropy oracle ≤ 1e-4 AU", max_dR < 1e-4);
    check("sun_meeus tracked accuracy MAJORIZES measured astropy residual", majorized);
}

void test_solar_meeus_precision() {
    auto mk = [](auto dummy) {
        using U = decltype(dummy);
        auto ep = astronomy::Epoch<U>::from_jd(tv<U>(2459000.5), astronomy::TimeScale::TT);
        return ephemeris::sun_meeus_of_date<U>(ep).lon.errors.precision;
    };
    double pd = static_cast<double>(mk(double(0)));
    double pb = static_cast<double>(mk(cpp_bin_float_50(0)));
    check("sun_meeus λ precision > 0 (double)", pd > 0.0);
    check("sun_meeus λ precision tightens with wider T", pb < pd * 1e-20);
}

// The GENERATIVE Fourier-Bessel equation-of-center (orbit/kepler_series.h) — the
// dial-up, tracked, extensible series the user directed (2026-06-07).
void test_kepler_series_generator() {
    using T = double;
    auto exact_C = [](double M, double e) {
        double E = M;
        for (int i = 0; i < 200; ++i) E = E - (E - e * std::sin(E) - M) / (1 - e * std::cos(E));
        double nu = std::atan2(std::sqrt(1 - e * e) * std::sin(E), std::cos(E) - e);
        double C = nu - M;
        while (C > kPI) C -= 2 * kPI;
        while (C < -kPI) C += 2 * kPI;
        return C;
    };
    auto maxerr = [&](double e, int K) {
        double m = 0.0;
        for (double M = 0.1; M < 6.28; M += 0.1) {
            auto C = orbit::equation_of_center<T>(tv<T>(M), tv<T>(e), K);
            m = std::max(m, std::abs(C.value - exact_C(M, e)));
        }
        return m;
    };
    // Reproduces the EXACT Kepler solution to machine precision at K=8 (ephemeris e).
    check("EoC series = exact Kepler @K8 (Sun)", maxerr(0.016709, 8) < 1e-12);
    check("EoC series = exact Kepler @K8 (Moon)", maxerr(0.054900, 8) < 1e-10);
    // Dial-up: more harmonics → strictly smaller error.
    check("EoC series dials up (e=0.3: K8<K4<K2)",
          maxerr(0.3, 8) < maxerr(0.3, 4) && maxerr(0.3, 4) < maxerr(0.3, 2));
    // Extensible PAST the Laplace limit e=0.6627 (where any e-power series diverges).
    check("EoC series converges past Laplace limit (e=0.75: K32 << K8)",
          maxerr(0.75, 32) < maxerr(0.75, 8) * 0.2);
    // The tracked total error budget majorizes the actual truncation residual.
    bool maj = true;
    for (double e : {0.016709, 0.054900, 0.30})
        for (int K : {2, 4, 8})
            for (double M = 0.1; M < 6.28; M += 0.2) {
                auto C = orbit::equation_of_center<T>(tv<T>(M), tv<T>(e), K);
                double tot = C.errors.precision + C.errors.accuracy;
                if (tot + 1e-18 < std::abs(C.value - exact_C(M, e))) maj = false;
            }
    check("EoC series tracked error majorizes the residual", maj);
    // Precision tightens with a wider numeric type T (the calling card).
    double pd = orbit::equation_of_center<double>(tv<double>(1.0), tv<double>(0.0549), 8).errors.precision;
    double pb = (double)orbit::equation_of_center<cpp_bin_float_50>(
        tv<cpp_bin_float_50>(1.0), tv<cpp_bin_float_50>(0.0549), 8).errors.precision;
    check("EoC series precision tightens with wider T", pd > 0.0 && pb < pd * 1e-20);
}

// L3: the lunar Poisson-series instance (moon_meeus.h) vs the INDEPENDENT JPL
// DE430 oracle (numerical integration — NOT astropy's builtin Meeus/ELP Moon).
void test_moon_de430_oracle() {
    using T = double;
    // Born-digital JPL DE430 (in-repo) geocentric ecliptic-of-date λ[deg], β[deg],
    // Δ[km], from tools/gen_lunar_oracle.py (DE430 table + erfa.eqec06).
    struct Row { double jd, lon_deg, lat_deg, dist_km; };
    const Row oracle[] = {
        {2451545.0, 223.3189269,  5.1708693, 402448.640},
        {2455455.0, 265.5408049, -1.1805569, 390817.018},
        {2458849.5, 346.1336199, -4.8937547, 403859.527},
        {2460462.5, 357.9329146, -1.4223734, 368744.527},
        {2449796.5, 217.8519101,  0.1689323, 367127.677},
    };
    const double rad2deg = 180.0 / kPI;
    double max_dlon = 0, max_dlat = 0, max_ddist = 0; bool maj = true;
    for (const auto& r : oracle) {
        auto ep = astronomy::Epoch<T>::from_jd(tv<T>(r.jd), astronomy::TimeScale::TT);
        auto s = ephemeris::moon_meeus_of_date<T>(ep, 60);
        double dlon = s.lon.value * rad2deg - r.lon_deg;
        while (dlon > 180) dlon -= 360;
        while (dlon < -180) dlon += 360;
        double dlat = s.lat.value * rad2deg - r.lat_deg;
        double ddist = s.radius.value - r.dist_km;
        max_dlon = std::max(max_dlon, std::abs(dlon));
        max_dlat = std::max(max_dlat, std::abs(dlat));
        max_ddist = std::max(max_ddist, std::abs(ddist));
        if (std::abs(dlon) > s.lon.errors.accuracy * rad2deg) maj = false;
        if (std::abs(dlat) > s.lat.errors.accuracy * rad2deg) maj = false;
        if (std::abs(ddist) > s.radius.errors.accuracy) maj = false;
    }
    std::cout << "    moon_meeus vs JPL DE430: max|Δλ|=" << max_dlon * 3600.0 << "\"  |Δβ|="
              << max_dlat * 3600.0 << "\"  |ΔΔ|=" << max_ddist << " km\n";
    check("moon_meeus λ reproduces JPL DE430 ≤ 10\" (Meeus §47 grade)", max_dlon * 3600.0 < 10.0);
    check("moon_meeus β reproduces JPL DE430 ≤ 8\"", max_dlat * 3600.0 < 8.0);
    check("moon_meeus Δ reproduces JPL DE430 ≤ 30 km", max_ddist < 30.0);
    check("moon_meeus tracked accuracy MAJORIZES the DE430 residual", maj);
    // Dial-up: fewer Poisson terms → markedly coarser.
    auto ep0 = astronomy::Epoch<T>::from_jd(tv<T>(2455455.0), astronomy::TimeScale::TT);
    double e12 = std::abs(ephemeris::moon_meeus_of_date<T>(ep0, 12).lon.value * rad2deg - 265.5408049) * 3600.0;
    check("moon_meeus dials up (12-term ≫ coarser than 60-term)", e12 > max_dlon * 3600.0 * 3.0);
    // Precision tightens with a wider numeric type T.
    double pd = ephemeris::moon_meeus_of_date<double>(ep0, 60).lon.errors.precision;
    double pb = (double)ephemeris::moon_meeus_of_date<cpp_bin_float_50>(
        astronomy::Epoch<cpp_bin_float_50>::from_jd(tv<cpp_bin_float_50>(2455455.0), astronomy::TimeScale::TT), 60)
        .lon.errors.precision;
    check("moon_meeus λ precision tightens with wider T", pd > 0.0 && pb < pd * 1e-20);
}

} // namespace

int main() {
    test_solar();
    test_lunar();
    test_ut1_utc();
    test_precision_scales();
    test_solar_astropy_oracle();
    test_solar_meeus_precision();
    test_kepler_series_generator();
    test_moon_de430_oracle();

    std::cout << "\n  ephemeris: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
