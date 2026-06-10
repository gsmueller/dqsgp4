#pragma once

/**
 * @file element_recovery.h
 * @brief Recover mean orbital elements from TLE osculating elements.
 *
 * The Two-Line Element set provides osculating mean motion n₀ that includes
 * secular perturbation effects. To use it in the Brouwer analytical theory,
 * we must "recover" the mean elements (a₀″, n₀″) by iteratively removing
 * the J₂ secular effect.
 *
 * The iteration uses the Brouwer-Lyddane correction:
 *   1. Compute a₁ = (kₑ/n₀)^(2/3) — Keplerian semi-major axis
 *   2. Compute δ₁ = (3/2)k₂(3cos²i−1)/(a₁²β₀³) — first J₂ correction
 *   3. Compute a₀ = a₁(1 − δ₁/3 − δ₁² − 134/81·δ₁³) — corrected a
 *   4. Compute δ₀ using a₀ — refined correction
 *   5. n₀″ = n₀/(1+δ₀), a₀″ = a₀/(1−δ₀) — recovered mean elements
 *
 * The 134/81 coefficient comes from the third-order expansion of the
 * Brouwer generating function applied to the semi-major axis.
 *
 * @par References
 * - Brouwer (1959), Solution of the Problem of Artificial Satellite Theory
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, page 10
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"

namespace orbit {

/**
 * @brief Result of the mean element recovery.
 *
 * All quantities in SGP4 normalized units: Earth radii for distance,
 * radians/minute for angular rates.
 */
template<typename T>
struct RecoveredElements {
    math::TrackedValue<T> a0;          ///< Recovered semi-major axis [Earth radii]
    math::TrackedValue<T> n0;          ///< Recovered mean motion [rad/min]
    math::TrackedValue<T> perigee_km;  ///< Perigee altitude [km]
    math::TrackedValue<T> period_min;  ///< Orbital period [minutes]
    math::TrackedValue<T> beta0;       ///< √(1 − e₀²)
    bool is_deep_space;                ///< True if period ≥ 225 min
    bool use_simple_model;             ///< True if perigee < 220 km
};

/**
 * @brief Recover mean elements from TLE osculating elements.
 *
 * @param n_input     TLE mean motion [rad/min]
 * @param e0          TLE eccentricity
 * @param cos_i0      Cosine of TLE inclination
 * @param half_J2     J₂/2 (= CK2 in SGP4 notation)
 * @param ke          √(GM/aₑ³) × 60 [rad/min] (Gaussian gravitational constant)
 * @param re_km       Earth equatorial radius [km]
 * @param tolerance   Reserved (the cube roots are now closed-form via the
 *                    tracked pow primitive); retained for API stability.
 * @return RecoveredElements with all three errors propagated
 */
template<typename T>
RecoveredElements<T> recover_mean_elements(
    const math::TrackedValue<T>& n_input,
    const math::TrackedValue<T>& e0,
    const math::TrackedValue<T>& cos_i0,
    const math::TrackedValue<T>& half_J2,
    const math::TrackedValue<T>& ke,
    const math::TrackedValue<T>& re_km,
    const T& tolerance)
{
    using math::exact;
    using math::ratio;

    // The semi-major-axis cube roots are now closed-form (tracked pow, B2); the
    // tolerance parameter is retained for API stability but drives no iteration.
    (void)tolerance;

    RecoveredElements<T> result;

    // Trig combination: 3cos²i − 1 (= 2P₂(cos i))
    math::TrackedValue<T> cos2_i0 = cos_i0 * cos_i0;
    math::TrackedValue<T> x3thm1 = exact<T>(3) * cos2_i0 - exact<T>(1);

    // --- Step 1: Keplerian semi-major axis ---
    // a₁ = (kₑ/n₀)^(2/3), via the tracked pow primitive (B2) — matches the
    // reference SGP4 pow(xke/no, x2o3) directly and retires the open-coded Newton.
    math::TrackedValue<T> a1 = pow(ke / n_input, T(2) / T(3));

    // Eccentricity-derived quantities
    math::TrackedValue<T> eosq = e0 * e0;
    math::TrackedValue<T> betao2 = exact<T>(1) - eosq;
    math::TrackedValue<T> betao = sqrt(betao2);

    // --- Step 2: First J₂ correction ---
    // δ₁ = (3/2)·k₂·(3cos²i − 1) / (a₁²·β₀³)
    math::TrackedValue<T> del1 = ratio<T>(3, 2) * half_J2 * x3thm1 / (a1 * a1 * betao * betao2);

    // --- Step 3: Corrected semi-major axis ---
    // a₀ = a₁·(1 − δ₁/3 − δ₁² − 134/81·δ₁³)
    // The 134/81 comes from the third-order Brouwer generating function.
    math::TrackedValue<T> a0 = a1 * (exact<T>(1) - del1 / exact<T>(3) - del1 * del1
                     - ratio<T>(134, 81) * del1 * del1 * del1);

    // --- Step 4: Refined correction ---
    // δ₀ = (3/2)·k₂·(3cos²i − 1) / (a₀²·β₀³)
    math::TrackedValue<T> del0 = ratio<T>(3, 2) * half_J2 * x3thm1 / (a0 * a0 * betao * betao2);

    // --- Step 5: Recovered mean elements ---
    result.n0 = n_input / (exact<T>(1) + del0);

    // Compute a₀ from n₀ via Kepler's third law: a = (kₑ/n)^(2/3)
    // This matches the reference SGP4 (Vallado SGP4.c line 1097):
    //   ao = pow(xke/no_unkozai, x2o3)
    // NOT the approximation a₀ = a₀_corrected/(1−δ₀), which differs
    // at O(δ₀²) and produces systematic position errors at t>0.
    // a₀ = (kₑ/n₀)^(2/3), via the tracked pow primitive (B2) — matches the
    // reference SGP4 pow(xke/no_unkozai, x2o3) directly.
    result.a0 = pow(ke / result.n0, T(2) / T(3));

    result.beta0 = betao;

    // --- REQ-EF-7: Model-truncation bound on the dropped O(δ₁⁴) series tail ---
    // The Brouwer-Lyddane series reversion a₀/a₁ = 1 − δ₁/3 − δ₁² − (134/81)δ₁³ + O(δ₁⁴)
    // truncates at the cubic term. The next coefficient has magnitude O(1), so the
    // absolute tail bound on a₀ is |δ₁|⁴·a₁. Propagate to n₀″ via the Kepler-3rd-law
    // derivative |dn/da| = (3/2)·n/a (from n² = kₑ²/a³ ⇒ n = kₑ/a^(3/2)).
    using std::abs;
    using std::pow;
    T tail_bound = pow(abs(del1.value), T(4)) * a1.value;
    result.a0.errors.accuracy = result.a0.errors.accuracy + tail_bound;
    result.n0.errors.accuracy = result.n0.errors.accuracy
        + (T(3) / T(2)) * (abs(result.n0.value) / abs(result.a0.value)) * tail_bound;

    // --- Derived quantities ---
    result.perigee_km = (result.a0 * (exact<T>(1) - e0) - exact<T>(1)) * re_km;
    result.period_min = math::two_pi<T>() / result.n0;
    result.is_deep_space = (result.period_min.value >= T(225));
    result.use_simple_model = (result.perigee_km.value < T(220));

    return result;
}

} // namespace orbit
