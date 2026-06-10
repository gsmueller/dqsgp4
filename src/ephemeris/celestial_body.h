#pragma once

/**
 * @file celestial_body.h
 * @brief Orbital elements for a perturbing celestial body (Sun or Moon).
 *
 * ★ DISPOSITION (R4a): part of the SR3-HISTORICAL ephemeris family
 * (solar_ephemeris / lunar_ephemeris) — SUPERSEDED for fidelity by the
 * Epoch-typed Meeus instances (`sun_meeus.h` / `moon_meeus.h`, DE430-gated)
 * that the third-body/SRP forces consume. RETAINED, unwired from any
 * propagator, with its gate (test_celestial_body_fields) as the documented
 * SR3 element-set record — do not extend.
 *
 * Bundles all the orbital parameters needed to compute a body's
 * ephemeris and its gravitational perturbation on a satellite.
 * The body is explicitly named so there is no ambiguity about
 * which eccentricity, which mean anomaly rate, etc.
 */

#include "../math/tracked_value.h"
#include <string>

namespace ephemeris {

/**
 * @brief Geocentric orbital elements for a perturbing body (Sun or Moon).
 *
 * Every element describes the body's **apparent orbit as seen from Earth**.
 * For the Sun, this is Earth's heliocentric orbit reflected to the geocentric
 * frame. For the Moon, it's the physical geocentric lunar orbit.
 *
 * The distinction matters: `orbit_eccentricity` is the eccentricity of
 * the Earth-Sun system (0.0167) or the Earth-Moon system (0.0549), NOT
 * an intrinsic property of the body alone.
 *
 * All angular quantities in radians. Rates in rad/min or rad/day as noted.
 *
 * @tparam T  The underlying numeric type
 */
template<typename T>
struct CelestialBody {
    std::string name;  ///< "Sun" or "Moon" — identifies the perturber

    // --- Geocentric orbital rates ---
    math::TrackedValue<T> mean_anomaly_rate_rad_min;  ///< Geocentric mean anomaly advance [rad/min]
    math::TrackedValue<T> mean_anomaly_epoch;         ///< Geocentric mean anomaly at reference epoch [rad]

    // --- For bodies with distinct longitude rate (Moon) ---
    math::TrackedValue<T> longitude_rate_rad_day;     ///< Geocentric mean longitude rate [rad/day]
    math::TrackedValue<T> longitude_epoch;            ///< Geocentric mean longitude at epoch [rad]

    // --- Geocentric orbital shape ---
    math::TrackedValue<T> orbit_eccentricity;         ///< Eccentricity of the geocentric orbit
    math::TrackedValue<T> orbit_arg_perigee;          ///< Argument of perigee of geocentric orbit [rad]

    // --- Geocentric orbital plane ---
    math::TrackedValue<T> orbit_inclination;          ///< Inclination of geocentric orbit to ecliptic [rad]
    math::TrackedValue<T> node_epoch;                 ///< Ascending node at epoch [rad]
    math::TrackedValue<T> node_rate_rad_day;          ///< Node regression rate [rad/day]

    // (The third-body perturbation strength μ_body/r³ is supplied to
    //  perturbation::compute_third_body_rates as its `perturbation_coef`
    //  argument at evaluation time; it is not a stored CelestialBody field —
    //  a never-assigned `mu_over_r3` member was removed here, issue BUG1.)

    /// Construct the Sun from FundamentalConstants.
    template<typename FC>
    static CelestialBody make_sun(const FC& fc) {
        CelestialBody sun;
        sun.name = "Sun";

        // Solar mean anomaly rate: ZNS
        math::TrackedValue<T> two_pi = math::two_pi<T>();
        math::TrackedValue<T> min_per_day = fc.minutes_per_day;
        sun.mean_anomaly_rate_rad_min = two_pi / (fc.solar_anomaly_period_days * min_per_day);
        sun.mean_anomaly_epoch = fc.solar_anomaly_epoch;

        // Sun has no separate longitude rate — longitude = M + ω̃
        sun.longitude_rate_rad_day = two_pi / fc.solar_anomaly_period_days;
        sun.longitude_epoch = fc.solar_anomaly_epoch; // simplified

        sun.orbit_eccentricity = fc.solar_eccentricity;
        sun.orbit_arg_perigee = fc.solar_arg_perigee;

        // Sun's ecliptic latitude is 0 by definition
        sun.orbit_inclination = math::exact<T>(0);
        sun.node_epoch = math::exact<T>(0);
        sun.node_rate_rad_day = math::exact<T>(0);

        return sun;
    }

    /// Construct the Moon from FundamentalConstants.
    template<typename FC>
    static CelestialBody make_moon(const FC& fc) {
        CelestialBody moon;
        moon.name = "Moon";

        // Lunar mean anomaly rate: ZNL (anomalistic)
        math::TrackedValue<T> two_pi = math::two_pi<T>();
        math::TrackedValue<T> min_per_day = fc.minutes_per_day;
        moon.mean_anomaly_rate_rad_min = two_pi / (fc.lunar_anomalistic_period_days * min_per_day);
        moon.mean_anomaly_epoch = math::exact<T>(0); // from SR3 epoch handling

        // Lunar mean longitude rate: sidereal rate
        moon.longitude_rate_rad_day = two_pi / fc.lunar_sidereal_period_days;
        moon.longitude_epoch = fc.lunar_longitude_epoch;

        moon.orbit_eccentricity = fc.lunar_eccentricity;
        moon.orbit_arg_perigee = fc.lunar_perigee_epoch; // perigee longitude at epoch

        // Moon's ecliptic inclination ~5.145°
        moon.orbit_inclination = fc.lunar_inclination;
        moon.node_epoch = fc.lunar_node_epoch;
        // Node regresses: negative rate
        moon.node_rate_rad_day = -two_pi / fc.lunar_node_period_days;

        return moon;
    }
};

} // namespace ephemeris
