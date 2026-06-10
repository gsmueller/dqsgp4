#pragma once

/**
 * @file secular_update.h
 * @brief Advance mean orbital elements forward in time under secular perturbations and drag.
 *
 * Given the mean elements at epoch and the secular rates (from Brouwer theory)
 * plus drag coefficients (from Lane-Hoots theory), this module computes the
 * mean elements at an arbitrary future time.
 *
 * The time dependence is a polynomial expansion:
 *   a(t) = a₀ × (1 − C₁t − D₂t² − D₃t³ − D₄t⁴)²
 *   e(t) = e₀ − B*C₄t − B*C₅(sin M(t) − sin M₀)
 *   M(t) = M₀ + Ṁt + drag corrections
 *   ω(t) = ω₀ + ω̇t − drag corrections
 *   Ω(t) = Ω₀ + Ω̇t + Ω̇_drag·t²
 *
 * @par References
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, pages 11-12
 * - design/derivations/016_secular_time_advance.md
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "../atmosphere/drag_coefficients.h"
#include "../perturbation/brouwer.h"

namespace orbit {

/**
 * @brief Mean elements at a propagated time.
 */
template<typename T>
struct SecularState {
    math::TrackedValue<T> a;       ///< Semi-major axis at time t [Earth radii]
    math::TrackedValue<T> e;       ///< Eccentricity at time t
    math::TrackedValue<T> M;       ///< Mean anomaly at time t [rad]
    math::TrackedValue<T> omega;   ///< Argument of perigee at time t [rad]
    math::TrackedValue<T> Omega;   ///< RAAN at time t [rad]
    math::TrackedValue<T> n;       ///< Mean motion at time t [rad/min]
    math::TrackedValue<T> mean_longitude; ///< M + ω + Ω + drag correction [rad]
};

/**
 * @brief Advance mean elements from epoch by time tsince.
 *
 * @param a0, e0, M0, omega0, Omega0  Epoch mean elements
 * @param n0           Recovered mean motion [rad/min] (= no_unkozai in Vallado code)
 * @param rates        Brouwer secular rates (M_dot, omega_dot, Omega_dot)
 * @param drag         Drag coefficients from Lane-Hoots model
 * @param bstar        B* drag term from TLE
 * @param ke           √(GM/aₑ³)×60 [rad/min]
 * @param tsince       Time since epoch [minutes]
 * @param simple_model True if perigee < 220 km (omit higher-order drag)
 * @return SecularState at time tsince
 */
template<typename T>
SecularState<T> secular_advance(
    const math::TrackedValue<T>& a0,
    const math::TrackedValue<T>& e0,
    const math::TrackedValue<T>& M0,
    const math::TrackedValue<T>& omega0,
    const math::TrackedValue<T>& Omega0,
    const math::TrackedValue<T>& n0,
    const perturbation::BrouwerSecularRates<T>& rates,
    const atmosphere::DragCoefficients<T>& drag,
    const math::TrackedValue<T>& bstar,
    const math::TrackedValue<T>& ke,
    const math::TrackedValue<T>& tsince,
    bool simple_model)
{
    using math::exact;
    using math::ratio;

    SecularState<T> state;

    math::TrackedValue<T> tsq = tsince * tsince;
    math::TrackedValue<T> tcube = tsq * tsince;
    math::TrackedValue<T> tfour = tcube * tsince;

    // --- Secular advances ---
    math::TrackedValue<T> M_secular = M0 + rates.M_dot * tsince;
    math::TrackedValue<T> omega_secular = omega0 + rates.omega_dot * tsince;
    math::TrackedValue<T> Omega_secular = Omega0 + rates.Omega_dot * tsince;

    // RAAN quadratic correction from drag-gravity coupling
    state.Omega = Omega_secular + drag.Omega_dot_drag * tsq;

    // --- Drag accumulation ---
    math::TrackedValue<T> a_drag_factor = exact<T>(1) - drag.C1 * tsince;
    math::TrackedValue<T> e_drag = bstar * drag.C4 * tsince;
    math::TrackedValue<T> L_drag = drag.t2cof * tsq;

    // --- Non-simple model: higher-order drag corrections ---
    if (!simple_model) {
        math::TrackedValue<T> delomg = drag.omega_drag_coef * tsince;
        math::TrackedValue<T> delm = drag.M_drag_coef
            * ((exact<T>(1) + drag.eta * cos(M_secular))
               * (exact<T>(1) + drag.eta * cos(M_secular))
               * (exact<T>(1) + drag.eta * cos(M_secular))
               - drag.eta_cos_M0_cubed);
        math::TrackedValue<T> temp_drag = delomg + delm;

        state.M = M_secular + temp_drag;
        state.omega = omega_secular - temp_drag;

        a_drag_factor = a_drag_factor - drag.D2 * tsq - drag.D3 * tcube - drag.D4 * tfour;
        e_drag = e_drag + bstar * drag.C5 * (sin(state.M) - drag.sin_M0);
        L_drag = L_drag + drag.t3cof * tcube + tfour * (drag.t4cof + tsince * drag.t5cof);
    } else {
        state.M = M_secular;
        state.omega = omega_secular;
    }

    // --- Updated elements ---
    state.a = a0 * a_drag_factor * a_drag_factor;
    state.e = e0 - e_drag;

    // Eccentricity floor
    // The floor prevents downstream divisions by tiny e (long-period
    // corrections, eccentricity-vector decomposition). When it fires, the
    // true eccentricity lies somewhere in [0, 1e-6] but we report 1e-6 —
    // a model-approximation (not numerical) systematic offset. Per REQ-EF-9
    // ("catastrophic loss is signaled, not silenced"), add the full range
    // of the offset (1e-6) to errors.accuracy so total_error() bounds truth.
    if (state.e.value < T(1e-6)) {
        state.e = math::TrackedValue<T>(T(1e-6), state.e.errors);
        state.e.errors.accuracy = state.e.errors.accuracy + T(1e-6);
    }

    // Mean longitude computation matching the reference SGP4 arithmetic sequence:
    // The reference code wraps individual angles before computing mean longitude,
    // then recomputes M from the wrapped sum. This affects floating-point
    // cancellation and must be matched exactly for bit-identical results.
    //
    // Reference sequence (Vallado SGP4.c lines 1651-1659):
    //   mm = mm + no_unkozai * templ
    //   xlm = mm + argpm + nodem
    //   nodem = fmod(nodem, 2π)
    //   argpm = fmod(argpm, 2π)
    //   xlm = fmod(xlm, 2π)
    //   mm = fmod(xlm - argpm - nodem, 2π)
    state.M = state.M + n0 * L_drag;
    math::TrackedValue<T> xlm = state.M + state.omega + state.Omega;
    state.Omega = math::wrap_two_pi(state.Omega);
    state.omega = math::wrap_two_pi(state.omega);
    xlm = math::wrap_two_pi(xlm);
    state.M = math::wrap_two_pi(xlm - state.omega - state.Omega);
    state.mean_longitude = xlm;

    // Mean motion from updated semi-major axis
    state.n = ke / (state.a * sqrt(state.a));

    return state;
}

} // namespace orbit
