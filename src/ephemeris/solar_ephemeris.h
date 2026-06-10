#pragma once

/**
 * @file solar_ephemeris.h
 * @brief Low-precision solar ephemeris for SGP4 deep-space perturbations.
 *
 * ★ DISPOSITION (R4a, the professional-library re-architecture): this is the
 * SR3-HISTORICAL instance — 1970s-style mean elements on a relative clock.
 * SUPERSEDED for fidelity by `sun_meeus.h` (Epoch-typed, generative series,
 * DE430-gated — the instance the third-body/SRP forces consume). RETAINED,
 * unwired from any propagator, as the documented SR3 element set with its
 * gates (test_ephemeris) — do not extend it; extend the Meeus instance.
 *
 * Computes the Sun's geocentric ecliptic longitude and distance.
 * Accuracy ~0.01° — sufficient for the third-body perturbation.
 *
 * @par References
 * - Meeus, J. (1998), "Astronomical Algorithms", Chapter 25
 * - design/derivations/015_solar_ephemeris.md
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "celestial_body.h"

namespace ephemeris {

/**
 * @brief Solar position in ecliptic coordinates.
 */
template<typename T>
struct SolarPosition {
    math::TrackedValue<T> ecliptic_longitude;  ///< λ☉ [rad]
    math::TrackedValue<T> sin_longitude;       ///< sin(λ☉)
    math::TrackedValue<T> cos_longitude;       ///< cos(λ☉)
    math::TrackedValue<T> distance_au;         ///< r☉ [AU] (approximate)
    math::TrackedValue<T> mean_anomaly;        ///< M☉ [rad]
};

/**
 * @brief Compute the Sun's geocentric position.
 *
 * @param sun          Solar orbital elements (from CelestialBody::make_sun)
 * @param delta_t_min  Time since reference epoch [minutes]
 * @return SolarPosition with ecliptic longitude and distance
 */
template<typename T>
SolarPosition<T> compute_solar_position(
    const CelestialBody<T>& sun,
    const math::TrackedValue<T>& delta_t_min)
{
    using math::exact;
    using math::ratio;

    SolarPosition<T> pos;

    // Mean anomaly at time t
    pos.mean_anomaly = sun.mean_anomaly_epoch + sun.mean_anomaly_rate_rad_min * delta_t_min;

    // Equation of center (Lagrange-series, 2-term truncation):
    //   ν − M = 2e sin(M) + (5/4) e² sin(2M) + O(e³)
    //
    // TRUNCATION ANNOTATION (R14 remediation, 2026-05-13):
    //   The next omitted term in the Lagrange series for the
    //   equation of center is
    //       (13/12) e³ sin(3M),                         (Meeus §25.3)
    //   bounded in magnitude by (13/12)|e|³.
    //
    //   For the Sun, e ≈ 0.01671 (J2000 epoch), so the dropped term
    //   is bounded by
    //       |drop| ≤ (13/12)(0.01671)³ ≈ 5.0 × 10⁻⁸ rad
    //             ≈ 0.01 arcsec
    //       ≪ stated ephemeris accuracy target ~0.01° = 36 arcsec.
    //
    //   Convergence holds for e < e_Laplace ≈ 0.6627; the Sun's
    //   eccentricity is two orders of magnitude below this, so
    //   convergence is rapid and the O(e³) bound is uniform.
    //
    //   EPH: this omitted-term magnitude is now FOLDED into the returned
    //   value's errors.accuracy (model-fidelity), so the third-error budget
    //   honestly reflects the 2-term truncation rather than silently dropping
    //   it. It remains far below the deep-space model floor, but it is tracked.
    //
    //   References:
    //     - Meeus (1998) "Astronomical Algorithms" §25.3 — Lagrange
    //       series for equation of center.
    //     - Vallado, Crawford & Hujsa (2006) Revisiting SR3 §B5.
    //     - Brouwer & Clemence (1961) "Methods of Celestial Mechanics"
    //       Ch. 10 — Lagrange-inversion convergence and Laplace radius.
    using std::abs;
    math::TrackedValue<T> sin_M = sin(pos.mean_anomaly);
    math::TrackedValue<T> sin_2M = sin(exact<T>(2) * pos.mean_anomaly);
    math::TrackedValue<T> e = sun.orbit_eccentricity;
    math::TrackedValue<T> equation_of_center = exact<T>(2) * e * sin_M + ratio<T>(5, 4) * e * e * sin_2M;

    // Fold the dropped (13/12) e³ sin(3M) Lagrange term as a model-truncation
    // accuracy bound (|next term| ≤ (13/12)|e|³). This flows into the longitude.
    T abs_e = abs(e.value);
    T eoc_trunc = (T(13) / T(12)) * abs_e * abs_e * abs_e;
    equation_of_center.errors.accuracy = equation_of_center.errors.accuracy + eoc_trunc;

    // True anomaly
    math::TrackedValue<T> true_anomaly = pos.mean_anomaly + equation_of_center;

    // Ecliptic longitude = true anomaly + argument of perigee
    pos.ecliptic_longitude = true_anomaly + sun.orbit_arg_perigee;
    pos.sin_longitude = sin(pos.ecliptic_longitude);
    pos.cos_longitude = cos(pos.ecliptic_longitude);

    // Distance [AU]: r/a = 1 − e·cos(M) drops the O(e²) terms (+e²/2 −
    // (e²/2)cos 2M, Meeus §25). Carry that as a model-fidelity bound (≤ e²).
    pos.distance_au = exact<T>(1) - e * cos(pos.mean_anomaly);
    pos.distance_au.errors.accuracy = pos.distance_au.errors.accuracy + abs_e * abs_e;

    return pos;
}

} // namespace ephemeris
