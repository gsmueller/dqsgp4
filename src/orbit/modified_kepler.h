#pragma once

/**
 * @file modified_kepler.h
 * @brief SGP4 modified Kepler equation solvers.
 *
 * SGP4 solves Kepler's equation in a rotated form using x = (E + ω):
 *
 *   x + ayn·cos(x) − axN·sin(x) = U
 *
 * where axN = e·cos(ω), ayn = e·sin(ω) + long-period correction,
 * and U is the mean argument of latitude.
 *
 * This formulation avoids the singularity at e=0 (where ω is undefined)
 * and naturally incorporates the Lyddane eccentricity vector decomposition.
 *
 * Two solver algorithms are provided:
 * - Newton-Raphson: quadratic convergence, simple, ~6 iterations for double
 * - Halley: cubic convergence, ~4 iterations for double
 *
 * These functions are the single source for the SGP4 modified-Kepler solve;
 * the model_functions / model_selector factories delegate to them rather than
 * inlining their own copies.
 *
 * @par References
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, page 13.
 */

#include "../math/tracked_value.h"

namespace orbit {

/**
 * @brief Newton-Raphson solver for the SGP4 modified Kepler equation.
 *
 * f(x) = U − ayn·cos(x) + axN·sin(x) − x
 * f'(x) = ayn·sin(x) + axN·cos(x) − 1
 * Δx = f/f'
 *
 * Convergence: quadratic (doubles correct digits per iteration).
 */
template<typename T>
math::TrackedValue<T> solve_kepler_newton(
    const math::TrackedValue<T>& axN,
    const math::TrackedValue<T>& ayn,
    const math::TrackedValue<T>& U,
    const T& tolerance)
{
    math::TrackedValue<T> x = U;  // starting value
    math::TrackedValue<T> last_correction = math::exact<T>(0);
    for (int iter = 0; iter < 30; ++iter) {
        math::TrackedValue<T> sx = sin(x); math::TrackedValue<T> cx = cos(x);
        math::TrackedValue<T> f = U - ayn * cx + axN * sx - x;
        math::TrackedValue<T> fp = -ayn * sx - axN * cx + math::exact<T>(1);
        math::TrackedValue<T> delta = f / fp;
        x = x + delta;
        last_correction = delta;
        using std::abs;
        if (abs(delta.value) < tolerance) {
            x.errors.precision = x.errors.precision + abs(delta.value);
            return x;
        }
    }
    // Cap-hit fallback: loop exited without converging. Record the residual
    // of the last accepted correction as the precision contribution so the
    // caller's error budget reflects the non-converged state (REQ-EF-5).
    using std::abs;
    x.errors.precision = x.errors.precision + abs(last_correction.value);
    return x;
}

/**
 * @brief Halley solver for the SGP4 modified Kepler equation.
 *
 * Uses the second derivative for cubic convergence:
 * f''(x) = ayn·cos(x) − axN·sin(x)
 * Δx = 2·f·f' / (2·f'² − f·f'')
 *
 * Convergence: cubic (triples correct digits per iteration).
 */
template<typename T>
math::TrackedValue<T> solve_kepler_halley(
    const math::TrackedValue<T>& axN,
    const math::TrackedValue<T>& ayn,
    const math::TrackedValue<T>& U,
    const T& tolerance)
{
    math::TrackedValue<T> x = U;
    math::TrackedValue<T> last_correction = math::exact<T>(0);
    for (int iter = 0; iter < 15; ++iter) {
        math::TrackedValue<T> sx = sin(x); math::TrackedValue<T> cx = cos(x);
        math::TrackedValue<T> f = U - ayn * cx + axN * sx - x;
        math::TrackedValue<T> fp = -ayn * sx - axN * cx + math::exact<T>(1);
        math::TrackedValue<T> fpp = ayn * cx - axN * sx;
        math::TrackedValue<T> num = math::exact<T>(2) * f * fp;
        math::TrackedValue<T> den = math::exact<T>(2) * fp * fp - f * fpp;
        math::TrackedValue<T> delta = num / den;
        x = x + delta;
        last_correction = delta;
        using std::abs;
        if (abs(delta.value) < tolerance) {
            x.errors.precision = x.errors.precision + abs(delta.value);
            return x;
        }
    }
    // Cap-hit fallback: loop exited without converging. Record the residual
    // of the last accepted correction as the precision contribution so the
    // caller's error budget reflects the non-converged state (REQ-EF-5).
    using std::abs;
    x.errors.precision = x.errors.precision + abs(last_correction.value);
    return x;
}

} // namespace orbit
