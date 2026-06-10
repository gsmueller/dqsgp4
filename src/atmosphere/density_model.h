#pragma once

/**
 * @file density_model.h
 * @brief Power-law atmospheric density model for SGP4 drag computation.
 *
 * The Lane-Hoots atmospheric model uses a power-law density profile:
 *
 *   ρ(r) = ρ₀ × ((q₀ − s) / (r − s))⁴
 *
 * where:
 *   q₀ = 1 + 120/XKMPER — reference altitude (120 km above Earth surface)
 *   s  = 1 + s₀/XKMPER  — atmospheric fitting parameter (~78 km)
 *
 * The exponent 4 comes from an empirical fit to the actual atmospheric
 * density profile in the 150-500 km altitude range where most LEO
 * satellites orbit.
 *
 * For low-perigee orbits (perigee < 156 km), the atmospheric parameter s
 * is adjusted to track the actual perigee altitude, since the power-law
 * approximation breaks down when the satellite dips into denser atmosphere.
 *
 * @par References
 * - Lane, M.H. and Hoots, F.R. (1979), "General Perturbations Theories
 *   Derived from the 1965 Lane Drag Theory"
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, page 10
 */

#include "../math/tracked_value.h"

namespace atmosphere {

/**
 * @brief Atmospheric density parameters for the SGP4 drag model.
 */
template<typename T>
struct DensityParameters {
    math::TrackedValue<T> s;       ///< Atmospheric fitting parameter s [Earth radii]
    math::TrackedValue<T> qoms4;   ///< (q₀ − s)⁴ — density scaling factor
};

/**
 * @brief Compute atmospheric density parameters from perigee altitude.
 *
 * Three regimes:
 * - Perigee ≥ 156 km: standard parameters (s₀ = 78 km)
 * - 98 km ≤ perigee < 156 km: s adjusted to perigee altitude
 * - Perigee < 98 km: s set to 20 km (deep atmosphere)
 *
 * @param perigee_km   Perigee altitude [km]
 * @param a0           Recovered semi-major axis [Earth radii]
 * @param e0           Eccentricity
 * @param re_km        Earth equatorial radius [km]
 * @return DensityParameters with adjusted s and (q₀−s)⁴
 */
template<typename T>
DensityParameters<T> compute_density_parameters(
    const math::TrackedValue<T>& perigee_km,
    const math::TrackedValue<T>& a0,
    const math::TrackedValue<T>& e0,
    const math::TrackedValue<T>& re_km)
{
    using math::exact;

    DensityParameters<T> params;

    // Default: s = 1 + 78/XKMPER (78 km above surface, in Earth radii)
    params.s = exact<T>(1) + exact<T>(78) / re_km;

    // q₀ = 1 + 120/XKMPER (120 km reference altitude)
    math::TrackedValue<T> q0 = exact<T>(1) + exact<T>(120) / re_km;

    // (q₀ − s)⁴ = ((120 − 78)/XKMPER)⁴ = (42/XKMPER)⁴
    math::TrackedValue<T> diff = q0 - params.s;
    params.qoms4 = diff * diff * diff * diff;

    // --- Perigee-dependent adjustment ---
    //
    // The transition at perigee = 156 km is C⁰-continuous (s and (q₀-s)⁴
    // values match across the boundary by construction) but NOT C¹: the
    // derivatives ds/d(perigee) and d(qoms4)/d(perigee) jump across the
    // boundary. This introduces a ~1% drag-rate discontinuity for orbits
    // near the boundary. The discontinuity is a property of the Lane
    // piecewise fitting recipe (Def 2.3) and is documented as ERROR
    // SOURCE A-D4 in design/derivations/sgp4_near_earth_drag_theoretical_basis.md §2.
    // The branch selection consults perigee_km.value directly (not its
    // measurement uncertainty); callers whose perigee sits within its
    // uncertainty of 156 km should treat the resulting (s, qoms4) as
    // carrying an additional ~1% accuracy bound from this C¹ kink.
    if (perigee_km.value < T(156) && perigee_km.value >= T(98)) {
        // s* = perigee-altitude part (ER): s* = a₀(1−e₀) − s_default = (perige−78)/re,
        // so the adjusted s_new = 1 + s* = 1 + (perige−78)/re (reference SGP4 sfour).
        // FIX (2026-06-03): removed a spurious +1 here that added a full Earth radius
        // (~6378 km) to s for 98 ≤ perigee < 156 km — it made s_new = 1 + (correct s),
        // catastrophically corrupting the drag coefficients (e.g. sat 28350, perige 135.75 km).
        // Now consistent with the <98 km branch (s* = 20/re, the altitude part).
        math::TrackedValue<T> s_star = a0 * (exact<T>(1) - e0) - params.s;
        // Adjust (q₀−s)⁴: take 4th root, add (s − s_new), raise to 4th
        math::TrackedValue<T> qoms_root4 = sqrt(sqrt(params.qoms4));
        math::TrackedValue<T> new_qoms = qoms_root4 + params.s - (exact<T>(1) + s_star);
        params.qoms4 = new_qoms * new_qoms * new_qoms * new_qoms;
        params.s = exact<T>(1) + s_star;
    } else if (perigee_km.value < T(98)) {
        // Deep atmosphere: s₀ = 20 km
        math::TrackedValue<T> s_star = exact<T>(20) / re_km;
        math::TrackedValue<T> qoms_root4 = sqrt(sqrt(params.qoms4));
        math::TrackedValue<T> new_qoms = qoms_root4 + params.s - (exact<T>(1) + s_star);
        params.qoms4 = new_qoms * new_qoms * new_qoms * new_qoms;
        params.s = exact<T>(1) + s_star;
    }

    return params;
}

} // namespace atmosphere
