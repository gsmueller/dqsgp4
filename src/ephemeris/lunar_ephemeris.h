#pragma once

/**
 * @file lunar_ephemeris.h
 * @brief Low-precision lunar ephemeris for SGP4 deep-space perturbations.
 *
 * ★ DISPOSITION (R4a): the SR3-HISTORICAL instance — SUPERSEDED for fidelity
 * by `moon_meeus.h` (the 60-term Meeus §47 Poisson series, DE430-gated to
 * ≤3.7″ vs this file's ~2.5° lumped budget). RETAINED, unwired from any
 * propagator, as the documented SR3 element set with its gates — do not
 * extend it; extend the Meeus instance.
 *
 * Computes the Moon's geocentric ecliptic longitude and latitude.
 * Uses the anomalistic rate for mean anomaly and the sidereal rate for
 * mean longitude — these are physically distinct periods.
 *
 * @par References
 * - Meeus, J. (1998), "Astronomical Algorithms", Chapter 47
 * - design/derivations/016_lunar_ephemeris.md
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "celestial_body.h"

namespace ephemeris {

/**
 * @brief Lunar position in ecliptic coordinates.
 */
template<typename T>
struct LunarPosition {
    math::TrackedValue<T> ecliptic_longitude;  ///< λ [rad]
    math::TrackedValue<T> ecliptic_latitude;   ///< β [rad]
    math::TrackedValue<T> sin_longitude;       ///< sin(λ)
    math::TrackedValue<T> cos_longitude;       ///< cos(λ)
    math::TrackedValue<T> sin_latitude;        ///< sin(β)
    math::TrackedValue<T> cos_latitude;        ///< cos(β)
    math::TrackedValue<T> mean_anomaly;        ///< l [rad]
    math::TrackedValue<T> mean_longitude;      ///< L [rad]
    math::TrackedValue<T> node_longitude;      ///< Ω [rad]
};

/**
 * @brief Compute the Moon's geocentric position.
 *
 * @param moon         Lunar orbital elements (from CelestialBody::make_moon)
 * @param delta_t_min  Time since reference epoch [minutes]
 * @return LunarPosition with ecliptic longitude, latitude, and trig values
 */
template<typename T>
LunarPosition<T> compute_lunar_position(
    const CelestialBody<T>& moon,
    const math::TrackedValue<T>& delta_t_min)
{
    using math::exact;
    using math::ratio;
    using std::abs;

    LunarPosition<T> pos;

    // Time in days for the daily-rate quantities
    math::TrackedValue<T> delta_t_day = delta_t_min / exact<T>(1440);

    // Mean anomaly: anomalistic rate (perigee-to-perigee, ZNL)
    pos.mean_anomaly = moon.mean_anomaly_epoch + moon.mean_anomaly_rate_rad_min * delta_t_min;

    // Mean longitude: sidereal rate (star-to-star, distinct from anomalistic)
    pos.mean_longitude = moon.longitude_epoch + moon.longitude_rate_rad_day * delta_t_day;

    // Node longitude: regresses at node_rate
    pos.node_longitude = moon.node_epoch + moon.node_rate_rad_day * delta_t_day;

    // Equation of center, 2-term Lagrange series (Meeus §47):
    //   ν − l = 2e sin(l) + (5/4) e² sin(2l) + O(e³)
    // The main term is ~6.29° for e=0.0549; the (5/4)e² sin(2l) term ~0.21° is
    // Meeus's sin(2l) coefficient (a real fidelity gain over the 1-term form).
    math::TrackedValue<T> sin_l = sin(pos.mean_anomaly);
    math::TrackedValue<T> sin_2l = sin(exact<T>(2) * pos.mean_anomaly);
    math::TrackedValue<T> e = moon.orbit_eccentricity;
    math::TrackedValue<T> equation_of_center =
        exact<T>(2) * e * sin_l + ratio<T>(5, 4) * e * e * sin_2l;

    const T deg = boost::math::constants::pi<T>() / T(180);

    // Ecliptic longitude.
    pos.ecliptic_longitude = pos.mean_longitude + equation_of_center;
    // Model-fidelity (accuracy): this low-precision model omits the major lunar
    // periodics that require the solar elongation D and solar anomaly M' —
    // evection (1.27°), variation (0.66°), annual equation (0.19°), and the
    // smaller Meeus §47 terms. Their summed magnitude bounds the longitude
    // residual at ~2.5° — carried explicitly (the register's "weakest link"),
    // so the third-error budget honestly reflects this is a degree-grade model.
    pos.ecliptic_longitude.errors.accuracy =
        pos.ecliptic_longitude.errors.accuracy + (T(5) / T(2)) * deg;  // 2.5°
    pos.sin_longitude = sin(pos.ecliptic_longitude);
    pos.cos_longitude = cos(pos.ecliptic_longitude);

    // Ecliptic latitude: β ≈ i·sin(F), F = L − Ω (main term). Omitted latitude
    // periodics (Meeus §47: 0.28° sin(F±M'), 0.17° sin(F−2D), …) sum to ~0.9°.
    math::TrackedValue<T> arg_latitude = pos.mean_longitude - pos.node_longitude;
    pos.ecliptic_latitude = moon.orbit_inclination * sin(arg_latitude);
    pos.ecliptic_latitude.errors.accuracy =
        pos.ecliptic_latitude.errors.accuracy + (T(9) / T(10)) * deg;  // 0.9°
    pos.sin_latitude = sin(pos.ecliptic_latitude);
    pos.cos_latitude = cos(pos.ecliptic_latitude);

    return pos;
}

} // namespace ephemeris
