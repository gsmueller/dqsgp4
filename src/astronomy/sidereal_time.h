#pragma once

/**
 * @file sidereal_time.h
 * @brief Greenwich Mean Sidereal Time (GMST) computation.
 *
 * GMST is the hour angle of the mean vernal equinox at Greenwich.
 * It determines Earth's rotational orientation and is needed for the
 * TEME→ECEF frame transformation.
 *
 * The polynomial coefficients come from the Aoki et al. (1982) expression,
 * adopted by the IERS and used in the SGP4 THETAG function.
 *
 * @par References
 * - Aoki, S. et al. (1982), "The new definition of Universal Time",
 *   Astronomy & Astrophysics 105, pp 359-361
 * - IERS Conventions (2010), Chapter 5.5
 * - design/derivations/014_sidereal_time.md
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "epoch.h"   // L1: centuries_since_j2000 view (single-sourced J2000 base)

namespace astronomy {

/**
 * @brief Compute GMST at a given Julian date.
 *
 * Uses the Aoki et al. (1982) polynomial:
 *
 *   GMST(0h) = 24110.54841 + 8640184.812866·T + 0.093104·T² − 6.2e-6·T³
 *
 * in seconds of sidereal time, where T is Julian centuries of UT1 from J2000.0.
 * For times not at 0h UT1, adds the Earth rotation since midnight.
 *
 * The coefficients are derived from:
 * - 24110.54841: GMST at J2000.0 midnight (VLBI measurement)
 * - 8640184.812866: sidereal seconds per Julian century (from the
 *   ratio of sidereal to solar day: 366.2422/365.2422)
 * - 0.093104: quadratic precession correction
 * - −6.2e-6: cubic precession correction
 *
 * @param jd_ut1  Julian date (UT1 time scale)
 * @return GMST in radians [0, 2π), as TrackedValue with precision error
 *         from the polynomial evaluation and measurement error from
 *         the coefficient uncertainties (~0.1 ms in Earth orientation).
 */
template<typename T>
math::TrackedValue<T> compute_gmst(const math::TrackedValue<T>& jd_ut1) {
    using math::exact;
    using math::TrackedValue;

    // Julian centuries from J2000.0 via the single-sourced view (L1). Bit-
    // identical .value to the prior (jd - 2451545)/36525; the J2000 base now
    // lives once in astronomy/epoch.h as the exact-by-convention integer.
    math::TrackedValue<T> t_ut1 = centuries_since_j2000(jd_ut1);

    // Separate the Julian date into 0h UT1 and fractional day.
    // JD at 0h UT1 = floor(JD - 0.5) + 0.5
    // Fractional day since 0h = JD - JD_0h
    using std::floor;
    T jd_0h_val = floor(jd_ut1.value - T(0.5)) + T(0.5);
    math::TrackedValue<T> jd_0h = TrackedValue<T>(jd_0h_val, jd_ut1.errors);
    math::TrackedValue<T> t_0h = centuries_since_j2000(jd_0h);
    math::TrackedValue<T> frac_day = jd_ut1 - jd_0h;  // fraction of day since 0h UT1 [days]

    // GMST at 0h UT1 [seconds of sidereal time]
    // Polynomial: c₀ + c₁·T + c₂·T² + c₃·T³ — the Aoki et al. (1982) fit.
    // c₀/c₁ carry GENUINE adoption uncertainties (σ = 100× their decimal-ULP:
    // ±1 ms GMST-at-epoch from VLBI, ±0.1 ms/century rate) → measured. c₂/c₃
    // had σ EXACTLY equal to their decimal-ULP — a written-digit truncation
    // mis-filed as measurement noise (the CR1-b error), with no separately
    // published σ → model_coefficient (finite-digit fit coefficient: digits →
    // accuracy, binary storage → T-scaling precision). Value-preserving.
    math::TrackedValue<T> c0 = TrackedValue<T>::measured("24110.54841", "0.001");       // ±1 ms in GMST
    math::TrackedValue<T> c1 = TrackedValue<T>::measured("8640184.812866", "0.0001");   // ±0.1 ms/century
    math::TrackedValue<T> c2 = TrackedValue<T>::model_coefficient("0.093104");          // quadratic fit term
    math::TrackedValue<T> c3 = TrackedValue<T>::model_coefficient("-0.0000062");        // cubic fit term

    math::TrackedValue<T> gmst_0h_sec = c0 + t_0h * (c1 + t_0h * (c2 + t_0h * c3));

    // Add Earth rotation since 0h UT1.
    // Sidereal/solar day ratio 1.00273790935: a DERIVED/adopted ratio
    // (≈ 366.2422/365.2422, with the IERS precession correction), NOT a
    // convention — so it is not defined(). No separately published σ →
    // model_coefficient (digits → accuracy, storage → T-scaling precision).
    // ★ PROVENANCE CLOSED (O4, ratified in test_nutation): the adopted value
    // IS this polynomial's own rate, exactly — 1 + c1/(86400·36525) =
    // 1.0027379093508 reproduces all written digits. The value stays the
    // adopted model_coefficient (bit-frozen); its origin is now gate-verified.
    // GMST = GMST(0h) + 1.00273790935 × UT1_seconds
    math::TrackedValue<T> sidereal_ratio = TrackedValue<T>::model_coefficient("1.00273790935");
    math::TrackedValue<T> seconds_per_day = TrackedValue<T>::defined("86400.0");
    math::TrackedValue<T> ut1_seconds = frac_day * seconds_per_day;
    math::TrackedValue<T> gmst_sec = gmst_0h_sec + sidereal_ratio * ut1_seconds;

    // Convert seconds of sidereal time to radians: θ = GMST × 2π/86400
    math::TrackedValue<T> two_pi = math::two_pi<T>();
    math::TrackedValue<T> theta = gmst_sec * two_pi / seconds_per_day;

    // Wrap to [0, 2π)
    return math::wrap_two_pi(theta);
}

/// GMST from a UTC Julian date with an explicit UT1−UTC correction (EPH).
///
/// GMST is defined on the UT1 timescale, but a TLE/OMM epoch is nominally UTC.
/// Feeding UTC directly (ΔUT1 = 0) is the SGP4 convention and is what
/// compute_gmst(jd) does; the omission costs up to |ΔUT1| ≤ 0.9 s of rotation
/// (IERS bounds it by leap seconds), i.e. ≤ ~3.8 milli-degrees. Pass the IERS
/// Bulletin-A/B ΔUT1 [seconds] to apply it. ΔUT1 = 0 reproduces compute_gmst
/// EXACTLY (UT1 = UTC + ΔUT1/86400 days), so the authentic SGP4 path — which
/// calls compute_gmst directly — is unaffected.
template<typename T>
math::TrackedValue<T> compute_gmst_ut1(
    const math::TrackedValue<T>& jd_utc,
    const math::TrackedValue<T>& delta_ut1_seconds)
{
    math::TrackedValue<T> seconds_per_day = math::TrackedValue<T>::defined("86400.0");
    math::TrackedValue<T> jd_ut1 = jd_utc + delta_ut1_seconds / seconds_per_day;
    return compute_gmst(jd_ut1);
}

} // namespace astronomy
