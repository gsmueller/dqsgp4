#pragma once

/**
 * @file osculating_elements.h
 * @brief Convert Kepler solution to osculating orbital quantities.
 *
 * After solving the modified Kepler equation for (E+ω), this module
 * extracts the osculating orbital quantities: radius, radial velocity,
 * transverse velocity, and argument of latitude.
 *
 * @par References
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, page 13-14
 */

#include "../math/tracked_value.h"

namespace orbit {

/**
 * @brief Osculating orbital quantities from the Kepler solution.
 */
template<typename T>
struct OsculatingState {
    math::TrackedValue<T> r;       ///< Radius [Earth radii]
    math::TrackedValue<T> rdot;    ///< Radial velocity dr/dt [ER/min]
    math::TrackedValue<T> rfdot;   ///< Transverse velocity r·dφ/dt [ER/min]
    math::TrackedValue<T> sin_u;   ///< sin(argument of latitude)
    math::TrackedValue<T> cos_u;   ///< cos(argument of latitude)
    math::TrackedValue<T> u;       ///< Argument of latitude [rad]
    math::TrackedValue<T> sin_2u;  ///< sin(2u) for short-period corrections
    math::TrackedValue<T> cos_2u;  ///< cos(2u) for short-period corrections
    math::TrackedValue<T> e_sq;    ///< axN² + ayn² (eccentricity vector magnitude²)
    math::TrackedValue<T> pl;      ///< Semi-latus rectum p = a(1 − e²) [ER]
    math::TrackedValue<T> beta_l;  ///< √(1 − e²) at this instant
};

/**
 * @brief Extract osculating quantities from the modified Kepler solution.
 *
 * @param E_plus_w   (E+ω) from Kepler solver [rad]
 * @param axN        e·cos(ω) — eccentricity vector x-component
 * @param ayn        e·sin(ω) + long-period correction — eccentricity vector y-component
 * @param a          Semi-major axis at this time [Earth radii]
 * @param ke         √(GM/aₑ³)×60 [rad/min]
 * @return OsculatingState with all quantities and errors
 */
template<typename T>
OsculatingState<T> compute_osculating(
    const math::TrackedValue<T>& E_plus_w,
    const math::TrackedValue<T>& axN,
    const math::TrackedValue<T>& ayn,
    const math::TrackedValue<T>& a,
    const math::TrackedValue<T>& ke)
{
    using math::exact;

    OsculatingState<T> osc;

    math::TrackedValue<T> sinepw = sin(E_plus_w);
    math::TrackedValue<T> cosepw = cos(E_plus_w);

    // e·cos(E) and e·sin(E) from the eccentricity vector
    math::TrackedValue<T> ecose = axN * cosepw + ayn * sinepw;
    math::TrackedValue<T> esine = axN * sinepw - ayn * cosepw;

    osc.e_sq = axN * axN + ayn * ayn;
    osc.pl = a * (exact<T>(1) - osc.e_sq);
    osc.r = a * (exact<T>(1) - ecose);
    osc.rdot = ke * sqrt(a) * esine / osc.r;
    osc.rfdot = ke * sqrt(osc.pl) / osc.r;

    osc.beta_l = sqrt(exact<T>(1) - osc.e_sq);
    math::TrackedValue<T> temp_kep = exact<T>(1) / (exact<T>(1) + osc.beta_l);

    osc.cos_u = (a / osc.r) * (cosepw - axN + ayn * esine * temp_kep);
    osc.sin_u = (a / osc.r) * (sinepw - ayn - axN * esine * temp_kep);

    osc.u = atan2(osc.sin_u, osc.cos_u);
    osc.sin_2u = exact<T>(2) * osc.sin_u * osc.cos_u;
    osc.cos_2u = exact<T>(2) * osc.cos_u * osc.cos_u - exact<T>(1);

    return osc;
}

} // namespace orbit
