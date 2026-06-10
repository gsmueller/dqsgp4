#pragma once

/**
 * @file short_period.h
 * @brief J₂ short-period corrections to osculating elements.
 *
 * The Brouwer theory separates perturbations into secular (growing with time),
 * long-period (varying with argument of perigee, period ~months), and
 * short-period (varying with argument of latitude, period ~90 min).
 *
 * This module applies the short-period corrections from the J₂ first-order
 * generating function S₁. These correct the osculating radius, argument of
 * latitude, inclination, RAAN, and velocities for the periodic J₂ effect
 * within a single orbit.
 *
 * These corrections are shared between near-space and deep-space propagation.
 *
 * @par References
 * - Brouwer (1959), Section on short-period generating function
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, page 14
 */

#include "../math/tracked_value.h"

namespace perturbation {

/**
 * @brief Corrected orbital elements after short-period J₂ corrections.
 */
template<typename T>
struct CorrectedElements {
    math::TrackedValue<T> r;       ///< Corrected radius [Earth radii]
    math::TrackedValue<T> u;       ///< Corrected argument of latitude [rad]
    math::TrackedValue<T> i;       ///< Corrected inclination [rad]
    math::TrackedValue<T> Omega;   ///< Corrected RAAN [rad]
    math::TrackedValue<T> rdot;    ///< Corrected radial velocity [ER/min]
    math::TrackedValue<T> rfdot;   ///< Corrected transverse velocity [ER/min]
};

/**
 * @brief Apply J₂ short-period corrections to osculating elements.
 *
 * @param r_osc      Osculating radius [Earth radii]
 * @param u_osc      Osculating argument of latitude [rad]
 * @param i0         Inclination [rad] (from epoch — assumed constant for short-period)
 * @param Omega_sec  Secularly-advanced RAAN [rad]
 * @param rdot_osc   Osculating radial velocity [ER/min]
 * @param rfdot_osc  Osculating transverse velocity [ER/min]
 * @param n          Mean motion at this time [rad/min]
 * @param pl         Semi-latus rectum p = a(1−e²) [Earth radii]
 * @param beta_l     √(1 − e²) at this time
 * @param half_J2    J₂/2
 * @param cos2_i0    cos²(i₀)
 * @param sin2_i0    sin²(i₀) = 1 − cos²(i₀)
 * @param cos_i0     cos(i₀)
 * @param sin_i0     sin(i₀)
 * @param three_cos2i_minus_1  3cos²i₀ − 1
 * @param seven_cos2i_minus_1  7cos²i₀ − 1
 * @param sin_2u     sin(2u)
 * @param cos_2u     cos(2u)
 * @return CorrectedElements with short-period J₂ corrections applied.
 *
 * @par Accuracy limitation
 * These corrections are first-order in J₂ only. The omitted J₂² and J₄
 * short-period terms contribute ~O(J₂² × r) ≈ 1e-6 × 7000 km ≈ 7 mm.
 * This is not tracked as an explicit accuracy error because it is below
 * the model's inherent accuracy from the secular rate truncation (~100 m).
 */
template<typename T>
CorrectedElements<T> apply_short_period(
    const math::TrackedValue<T>& r_osc,
    const math::TrackedValue<T>& u_osc,
    const math::TrackedValue<T>& i0,
    const math::TrackedValue<T>& Omega_sec,
    const math::TrackedValue<T>& rdot_osc,
    const math::TrackedValue<T>& rfdot_osc,
    const math::TrackedValue<T>& n,
    const math::TrackedValue<T>& pl,
    const math::TrackedValue<T>& beta_l,
    const math::TrackedValue<T>& half_J2,
    const math::TrackedValue<T>& cos2_i0,
    const math::TrackedValue<T>& sin2_i0,
    const math::TrackedValue<T>& cos_i0,
    const math::TrackedValue<T>& sin_i0,
    const math::TrackedValue<T>& three_cos2i_minus_1,
    const math::TrackedValue<T>& seven_cos2i_minus_1,
    const math::TrackedValue<T>& sin_2u,
    const math::TrackedValue<T>& cos_2u)
{
    using math::exact;
    using math::ratio;

    CorrectedElements<T> corr;

    math::TrackedValue<T> inv_pl = exact<T>(1) / pl;
    math::TrackedValue<T> sp2 = half_J2 * inv_pl;          // CK2/p
    math::TrackedValue<T> sp3 = sp2 * inv_pl;              // CK2/p²

    // Radius correction: periodic J₂ in r
    corr.r = r_osc * (exact<T>(1) - ratio<T>(3, 2) * sp3 * beta_l * three_cos2i_minus_1)
             + ratio<T>(1, 2) * sp2 * sin2_i0 * cos_2u;

    // Argument of latitude correction
    corr.u = u_osc - ratio<T>(1, 4) * sp3 * seven_cos2i_minus_1 * sin_2u;

    // Inclination correction
    corr.i = i0 + ratio<T>(3, 2) * sp3 * cos_i0 * sin_i0 * cos_2u;

    // RAAN correction
    corr.Omega = Omega_sec + ratio<T>(3, 2) * sp3 * cos_i0 * sin_2u;

    // Velocity corrections
    corr.rdot = rdot_osc - n * sp2 * sin2_i0 * sin_2u;
    corr.rfdot = rfdot_osc + n * sp2 * (sin2_i0 * cos_2u + ratio<T>(3, 2) * three_cos2i_minus_1);

    return corr;
}

} // namespace perturbation
