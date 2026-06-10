#pragma once

/**
 * @file sun_meeus.h
 * @brief Modern (Meeus §25 / VSOP-truncation) solar ephemeris instance — L3.
 *
 * The Sun's GEOMETRIC geocentric position in the MEAN ecliptic & equinox OF DATE
 * (λ☉, r☉), as a function of an `Epoch<T>` (L1). This is the "modern instance"
 * of the L3 generic-ephemeris theory (design/derivations/ephemeris.md §3/§5),
 * distinct from — and independently verifiable against — the SR3 1970s instance
 * in solar_ephemeris.h. Unlike that one it is absolute-time (typed by `Epoch`,
 * not "minutes since an unstated epoch") and **astropy-verified**, curing the
 * "accuracy asserted, never measured" gap (PROFESSIONAL_LIBRARY_PLAN §2).
 *
 * @par Conventions (theory note §1, pinned)
 * - Time: TT Julian centuries since J2000 (L1 `centuries_since_j2000`).
 * - Frame: mean ecliptic & equinox OF DATE; reaching the equatorial/GCRF frame
 *   is the L2 obliquity/precession/nutation chain (a later L3 increment).
 * - GEOMETRIC (instantaneous, no light-time/aberration) — the quantity the
 *   third-body perturbation needs (gravity acts from where the body IS). The
 *   light-time/aberration offset to an apparent oracle (~20.5″) is handled in
 *   the gate, not the model.
 *
 * @par Fidelity (no perceived fidelity)
 * Meeus §25 low-precision; coefficients born-digital (pymeeus realises them
 * verbatim), each a finite-digit fit → `model_coefficient`. Verified against the
 * astropy/ERFA oracle to ≤ ~0.01° (the measured residual is ≤ ~0.007″·10³ =
 * 0.007°; the tracked accuracy carries Meeus's documented 0.01° grade, which
 * MAJORIZES it). Generator: tools/gen_solar_oracle.py; gate: test_ephemeris (EPH).
 *
 * @par References (born-digital)
 * - Meeus, Astronomical Algorithms (1998), Ch. 25 (Sun) — Eqs. 25.x.
 * - astropy `coordinates` + ERFA — the verification oracle.
 */

#include "../astronomy/earth_orbit.h"   // earth_eccentricity generator (Meeus 25.4 / VSOP)
#include "../astronomy/epoch.h"         // Epoch, centuries_since_j2000 (L1)
#include "../orbit/kepler_series.h"     // equation_of_center (Fourier-Bessel, generative)
#include "../math/angles.h"
#include "../math/tracked_value.h"
#include "ecliptic_state.h"             // shared EclipticState (also used by moon_meeus)

#include <boost/math/constants/constants.hpp>

namespace ephemeris {

/// Sun's geometric geocentric position in the mean ecliptic & equinox of date,
/// Meeus §25 low-precision. `epoch` is interpreted on the TT scale (the J2000
/// century view). Every coefficient is a born-digital finite-digit fit value.
template<typename T>
EclipticState<T> sun_meeus_of_date(const astronomy::Epoch<T>& epoch) {
    using TV = math::TrackedValue<T>;
    using math::exact;
    using math::degrees_to_radians;
    using math::wrap_two_pi;

    // TT Julian centuries since J2000 (L1 view).
    TV t = astronomy::centuries_since_j2000(epoch.jd());

    // Geometric mean longitude L0 and mean anomaly M [deg] (Meeus 25.x). Each
    // coefficient is a finite-digit VSOP fit → model_coefficient (digits →
    // accuracy, binary storage → T-scaling precision).
    TV L0_deg = TV::model_coefficient("280.46646")
              + TV::model_coefficient("36000.76983") * t
              + TV::model_coefficient("0.0003032") * t * t;
    TV M_deg = TV::model_coefficient("357.52911")
             + TV::model_coefficient("35999.05029") * t
             + TV::model_coefficient("-0.0001537") * t * t;
    TV M = degrees_to_radians(M_deg);

    // Earth's orbital (apparent-solar) eccentricity from the shared generator
    // (Meeus 25.4 / VSOP); carries its own series-truncation accuracy.
    TV e = astronomy::earth_eccentricity(t);

    // Equation of centre C [rad] from the GENERATIVE Fourier–Bessel series
    // (orbit/kepler_series.h; theory design/derivations/ephemeris_series.md) for
    // the time-varying eccentricity e(t): EXACT two-body, dial-up K, tracked
    // truncation — replacing the fixed 3-term Meeus form. K = 8 reaches machine
    // precision at the Sun's eccentricity, and the T-dependence Meeus folded into
    // the 3-term amplitudes (e.g. −0.004817·t·sinM = d(2e)/dt) is exactly the
    // e(t) dependence the generator already carries.
    TV C = orbit::equation_of_center(M, e, 8);

    EclipticState<T> s;
    // True (geometric) longitude = L0 + C, in the mean equinox of date.
    s.lon = wrap_two_pi(degrees_to_radians(L0_deg) + C);
    // Sun's ecliptic latitude is 0 by definition (the ecliptic IS Earth's
    // orbital plane); the ~0.7″ geocentric offset (Earth–Moon barycentre) is
    // carried as an accuracy bound, below this model's grade.
    s.lat = exact<T>(0);

    // Radius vector R = 1.000001018 · (1 − e²)/(1 + e cos ν) [AU], ν = M + C.
    TV nu = M + C;
    TV one = exact<T>(1);
    s.radius = TV::model_coefficient("1.000001018") * (one - e * e) / (one + e * cos(nu));

    // Model-fidelity ACCURACY (born-digital). Longitude: Meeus §25's documented
    // ≈ 0.01° low-precision grade (the omitted planetary periodics), verified to
    // MAJORIZE the measured astropy residual (≤ ~0.007°). Latitude: ~2″ (the
    // Earth–Moon barycentre offset). Radius: ~1e-4 AU (majorizes the ≤ 3.4e-5 AU
    // measured residual). These are the model floors, NOT typographic.
    const T arcsec = boost::math::constants::pi<T>() / T(180) / T(3600);
    s.lon = math::add_bound(s.lon, T(36) * arcsec, math::ErrorChannel::accuracy);   // 0.01° = 36″
    s.lat = math::add_bound(s.lat, T(2) * arcsec, math::ErrorChannel::accuracy);
    s.radius = math::add_bound(s.radius, T(1) / T(10000), math::ErrorChannel::accuracy);  // 1e-4 AU

    return s;
}

} // namespace ephemeris
