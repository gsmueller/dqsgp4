#pragma once

/**
 * @file model_functions.h
 * @brief Lambda function types and ModelFunctions bundle for injectable computation.
 *
 * The per-step propagator computations — secular rates, the modified-Kepler
 * solve, and sidereal time — are provided as injectable lambdas in the
 * ModelFunctions struct; the propagator calls these lambdas, so swapping a
 * perturbation theory, solver, or sidereal-time model needs no propagator code
 * change (the test_injection gate proves a swapped kepler_solver changes the
 * result). To change the geodetic model, provide a different ModelFunctions
 * instance.
 *
 * Exception (INJ1, honest scope): the deep-space resonance inclination functions
 * are NOT injected — perturbation/resonance.h uses the bit-tuned SR3 closed
 * forms directly, because they are not all bit-identical to the general Kaula
 * F_lmp and routing them through one lambda would drift the regression-frozen
 * SGP4 path. The former dead `inclination_function` slot was removed rather than
 * left as a misleading injection point; the Kaula function stays reusable as
 * perturbation::inclination_function.
 *
 * @par References
 * - design/lambda_injection.md
 * - design/library_architecture.md
 */

#include "../math/tracked_value.h"
#include "../perturbation/brouwer.h"
#include "../perturbation/kaula.h"
#include "../math/kepler.h"
#include "../orbit/modified_kepler.h"
#include "../astronomy/sidereal_time.h"
#include <functional>
#include <stdexcept>
#include <string>

namespace sgp4 {

// ============================================================================
// Lambda function type definitions
// ============================================================================

/**
 * @brief Compute Brouwer secular rates from orbital parameters.
 * @param n     Recovered mean motion [rad/min]
 * @param a     Recovered semi-major axis [Earth radii]
 * @param e2    Eccentricity squared
 * @param cos_i Cosine of inclination
 * @param J2    Second zonal harmonic
 * @param J4    Fourth zonal harmonic
 * @return Secular rates (M_dot, omega_dot, Omega_dot) with three errors
 */
template<typename T>
using SecularRatesFn = std::function<
    perturbation::BrouwerSecularRates<T>(
        const math::TrackedValue<T>& n,
        const math::TrackedValue<T>& a,
        const math::TrackedValue<T>& e2,
        const math::TrackedValue<T>& cos_i,
        const math::TrackedValue<T>& J2,
        const math::TrackedValue<T>& J4)>;

/**
 * @brief Solve the SGP4 modified Kepler equation for (E + ω).
 *
 * The SGP4 Kepler equation uses the eccentricity vector components:
 *   x + ayn·cos(x) - axN·sin(x) = U
 * where x = E + ω, and axN = e·cos(ω), ayn = e·sin(ω) + long-period correction.
 *
 * This is algebraically equivalent to the standard Kepler equation
 * E - e·sin(E) = M but parameterized to avoid the ω singularity at e=0
 * and to incorporate Lyddane long-period corrections naturally.
 *
 * The standard form is a special case: set axN=e, ayn=0, U=M.
 *
 * Reference: [SR3] page 13, equations for U, (E+ω), and Δ(E+ω).
 *
 * All three error types (measurement, precision, accuracy) propagate
 * through the iteration. The convergence residual is added to the
 * precision error of the result.
 *
 * @param axN       x-component of eccentricity vector (e·cos ω)
 * @param ayn       y-component of eccentricity vector (e·sin ω + long-period correction)
 * @param U         Mean argument of latitude (IL_T - Ω) [rad]
 * @param tolerance Convergence tolerance
 * @return (E + ω) [rad] with precision error from convergence
 */
template<typename T>
using KeplerFn = std::function<
    math::TrackedValue<T>(
        const math::TrackedValue<T>& axN,
        const math::TrackedValue<T>& ayn,
        const math::TrackedValue<T>& U,
        const T& tolerance)>;

/**
 * @brief Compute Greenwich sidereal time from Julian date.
 * @param epoch_jd Julian date
 * @return Greenwich sidereal time [rad]
 */
template<typename T>
using SiderealTimeFn = std::function<
    math::TrackedValue<T>(
        const math::TrackedValue<T>& epoch_jd)>;

// ============================================================================
// ModelFunctions — bundles all injectable lambdas
// ============================================================================

/**
 * @brief Bundle of all computational lambdas used by the propagator.
 *
 * The propagator calls these functions instead of calling domain modules
 * directly. This enables swapping perturbation theories, ephemeris models,
 * solver algorithms, and drag models without changing propagator code.
 *
 * Use ModelSelector to obtain pre-configured instances, or construct
 * manually for custom configurations.
 */
template<typename T>
struct ModelFunctions {
    // ==================================================================
    // COMPUTE-ONCE functions: called during initialization only.
    // Their results are cached in NearSpaceInit / DeepSpaceInit.
    // These could be replaced by precomputed structs in a future refactor.
    // ==================================================================

    /// Brouwer (or other) secular rate computation.
    /// Called once per satellite at init. Result stored in NearSpaceInit.
    SecularRatesFn<T> secular_rates;

    // (INJ1: the former `inclination_function` slot was removed — it was
    //  assigned but never called by the engine. The deep-space resonance uses
    //  the bit-tuned SR3 closed forms in perturbation/resonance.h directly, which
    //  are NOT all bit-identical to the general Kaula F_lmp (e.g. the half-day
    //  f220 = (3/4)(1+2c+c²) vs the synchronous (3/4)(1+c)²), so routing them
    //  through one injected F_lmp would drift the frozen SGP4 path. The Kaula
    //  inclination function remains reusable as perturbation::inclination_function.)

    // ==================================================================
    // PER-CALL functions: called every propagation step with new inputs.
    // These MUST remain as lambdas — they have time-varying arguments.
    // ==================================================================

    /// SGP4 modified Kepler equation solver (Newton, Halley, etc.)
    /// Solves: x + ayn·cos(x) - axN·sin(x) = U for x = (E + ω)
    /// Called every propagate() with different U each time.
    KeplerFn<T> kepler_solver;

    /// Greenwich sidereal time computation (Aoki 1982, IAU 2006 ERA, etc.)
    /// Called every propagate() with different Julian date.
    SiderealTimeFn<T> sidereal_time;

    /// Human-readable description of this configuration
    std::string description;

    // ------------------------------------------------------------------
    // Factory: Standard SGP4 (wraps existing module implementations)
    // ------------------------------------------------------------------

    /**
     * @brief Standard SGP4 model functions.
     *
     * Wraps the existing perturbation/brouwer.h, perturbation/kaula.h,
     * and math/kepler.h implementations as lambdas. These are the
     * computations used in Spacetrack Report #3.
     *
     * @par Secular Rates
     * Brouwer (1959) truncated at J₂² and J₄.
     * Polynomial coefficients: 13, 78, 137 for M_dot; 7, 114, 395 for omega_dot.
     *
     * @par Inclination Functions
     * Kaula (1966) for l=2..5 with exact rational coefficients.
     *
     * @par Kepler Solver
     * Newton's method (SR3 page 13); quadratic convergence in a neighbourhood
     * of the root for e < 1.
     *
     * @par Sidereal Time
     * Aoki et al. (1982) GMST polynomial — wired below to
     * astronomy::compute_gmst (used by the deep-space resonance integrator).
     */
    static ModelFunctions standard_sgp4() {
        ModelFunctions mf;

        mf.secular_rates = [](
            const math::TrackedValue<T>& n,
            const math::TrackedValue<T>& a,
            const math::TrackedValue<T>& e2,
            const math::TrackedValue<T>& cos_i,
            const math::TrackedValue<T>& J2,
            const math::TrackedValue<T>& J4)
        {
            return perturbation::compute_secular_rates(n, a, e2, cos_i, J2, J4);
        };

        /// Default SGP4 Kepler solver: Newton iteration on the modified form.
        /// Delegates to the single-source orbit::solve_kepler_newton so the
        /// solver body is not duplicated across model factories.
        /// Reference: [SR3] page 13.
        mf.kepler_solver = [](
            const math::TrackedValue<T>& axN,
            const math::TrackedValue<T>& ayn,
            const math::TrackedValue<T>& U,
            const T& tolerance) -> math::TrackedValue<T>
        {
            return orbit::solve_kepler_newton<T>(axN, ayn, U, tolerance);
        };

        // Sidereal time: Aoki et al. (1982) GMST polynomial.
        // Computes Greenwich Mean Sidereal Time from Julian date.
        // Required for TEME→ECEF frame transformation.
        mf.sidereal_time = [](const math::TrackedValue<T>& jd) -> math::TrackedValue<T> {
            return astronomy::compute_gmst(jd);
        };

        mf.description = "Standard SGP4: Brouwer J2^2+J4, Kaula l<=5, "
                          "Newton Kepler solver, Aoki 1982 GMST";
        return mf;
    }
};

} // namespace sgp4
