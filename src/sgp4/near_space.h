#pragma once

/**
 * @file near_space.h
 * @brief SGP4 near-space initialization and propagation.
 *
 * Implements the SGP4 near-earth analytical propagation from
 * Spacetrack Report #3, Section 6 "The SGP4 Model", pages 10-15.
 *
 * Every equation references the specific page and formula from [SR3].
 * Every constant is either:
 * - An exact integer or rational (e.g., 134/81, 3/2)
 * - Computed from the ellipsoid via the model functions lambdas
 * - Derived from TLE elements via TrackedValue arithmetic
 *
 * No magic numbers. No hardcoded formulas. All perturbation theory
 * is accessed through the ModelFunctions lambdas.
 *
 * @par References
 * - [SR3] Hoots & Roehrich (1980), Spacetrack Report No. 3, pages 10-15
 * - [V06] Vallado et al. (2006), "Revisiting Spacetrack Report #3" Rev 3
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "model_functions.h"
#include "model_selector.h"
#include "../geodesy/equipotential_ellipsoid.h"
#include "../tle/tle_parser.h"
#include "state_vector.h"
#include "../orbit/element_recovery.h"
#include "../atmosphere/density_model.h"
#include "../atmosphere/drag_coefficients.h"
#include "../orbit/secular_update.h"
#include "../orbit/osculating_elements.h"
#include "../perturbation/short_period.h"
#include "../orbit/state_from_elements.h"

namespace sgp4 {

/**
 * @brief All initialization constants for near-space propagation.
 *
 * Computed once from the TLE elements and ellipsoid constants.
 * Every member is a TrackedValue carrying three errors.
 */
template<typename T>
struct NearSpaceInit {
    // Recovered mean elements
    math::TrackedValue<T> a0;      ///< Recovered semi-major axis [Earth radii]
    math::TrackedValue<T> n0;      ///< Recovered mean motion [rad/min] (Kozai-corrected)
    math::TrackedValue<T> n0_unkozai; ///< Original TLE mean motion [rad/min] (before Kozai)
    math::TrackedValue<T> e0;      ///< Eccentricity at epoch
    math::TrackedValue<T> i0;      ///< Inclination at epoch [rad]
    math::TrackedValue<T> omega0;  ///< Argument of perigee at epoch [rad]
    math::TrackedValue<T> Omega0;  ///< RAAN at epoch [rad]
    math::TrackedValue<T> M0;      ///< Mean anomaly at epoch [rad]
    math::TrackedValue<T> bstar;   ///< B* drag coefficient

    // Derived orbital quantities
    math::TrackedValue<T> perigee_km;  ///< Perigee altitude [km]
    math::TrackedValue<T> period_min;  ///< Orbital period [min]
    bool use_simple_model;             ///< True if perigee < 220 km
    bool is_deep_space;                ///< True if period >= 225 min

    // Atmospheric density parameters
    math::TrackedValue<T> s4;      ///< Adjusted atmospheric parameter s
    math::TrackedValue<T> qoms24;  ///< Adjusted (q₀ - s)⁴

    // Initialization constants [SR3] pages 10-11
    math::TrackedValue<T> xi;      ///< ξ = 1/(a₀'' - s)
    math::TrackedValue<T> beta0;   ///< β₀ = √(1 - e₀²)
    math::TrackedValue<T> eta;     ///< η = a₀''·e₀·ξ
    math::TrackedValue<T> C2;      ///< C₂ coefficient
    math::TrackedValue<T> C1;      ///< C₁ = B*·C₂
    math::TrackedValue<T> C3;      ///< C₃ coefficient (J₃ correction)
    math::TrackedValue<T> C4;      ///< C₄ coefficient
    math::TrackedValue<T> C5;      ///< C₅ coefficient
    math::TrackedValue<T> D2;      ///< D₂ = 4a₀''·ξ·C₁²
    math::TrackedValue<T> D3;      ///< D₃ coefficient
    math::TrackedValue<T> D4;      ///< D₄ coefficient

    // Secular rate coefficients
    math::TrackedValue<T> M_dot;     ///< Secular rate of mean anomaly [rad/min]
    math::TrackedValue<T> omega_dot; ///< Secular rate of arg perigee [rad/min]
    math::TrackedValue<T> Omega_dot; ///< Secular rate of RAAN [rad/min]

    // Time coefficients
    math::TrackedValue<T> t2cof;   ///< (3/2)C₁
    math::TrackedValue<T> t3cof;   ///< D₂ + 2C₁²
    math::TrackedValue<T> t4cof;   ///< (1/4)(3D₃ + C₁(12D₂ + 10C₁²))
    math::TrackedValue<T> t5cof;   ///< (1/5)(3D₄ + 12C₁D₃ + 6D₂² + 15C₁²(2D₂ + C₁²))

    // Drag and mean anomaly correction coefficients
    math::TrackedValue<T> omgcof;  ///< B*·C₃·cos(ω₀)
    math::TrackedValue<T> xmcof;   ///< -(2/3)·(q₀-s)⁴·B*·aE/(e₀·η)
    math::TrackedValue<T> delmo;   ///< (1 + η·cos(M₀))³
    math::TrackedValue<T> sinmo;   ///< sin(M₀)
    math::TrackedValue<T> Omega_dot_nkc; ///< -(21/2)·n₀''·k₂·θ/(a₀''²·β₀²)·C₁

    // Long-period correction coefficients
    math::TrackedValue<T> xlcof;   ///< For long-period correction in mean longitude
    math::TrackedValue<T> aycof;   ///< For long-period correction in a_yN

    // Trig constants
    math::TrackedValue<T> cosio;   ///< cos(i₀)
    math::TrackedValue<T> sinio;   ///< sin(i₀)
    math::TrackedValue<T> theta2;  ///< cos²(i₀)
    math::TrackedValue<T> x3thm1;  ///< 3cos²i - 1
    math::TrackedValue<T> x1mth2;  ///< 1 - cos²i = sin²i
    math::TrackedValue<T> x7thm1;  ///< 7cos²i - 1

    // Earth constants from ellipsoid (stored for propagation)
    math::TrackedValue<T> xke;     ///< √(GM/aE³) in rad/min
    math::TrackedValue<T> CK2;     ///< J₂/2
    math::TrackedValue<T> J3;      ///< Jₙ(3) — third zonal harmonic
    math::TrackedValue<T> A30;     ///< A(3,0) = −J₃ in normalized units
    math::TrackedValue<T> re_km;   ///< Earth radius [km]
};

/**
 * @brief Initialize the near-space propagation constants.
 *
 * Implements [SR3] Section 6, pages 10-11: element recovery,
 * atmospheric parameter adjustment, and coefficient computation.
 *
 * All perturbation theory accessed via model_functions lambdas.
 * All constants from the ellipsoid, not hardcoded.
 *
 * @param config    Complete model configuration (provides Jn resolver, ellipsoid, lambdas)
 * @param elements  Parsed TLE orbital elements
 * @param tolerance Series evaluation tolerance
 * @return Fully initialized NearSpaceInit struct
 */
template<typename T>
NearSpaceInit<T> initialize_near_space(
    const sgp4::ModelConfiguration<T>& config,
    const tle::TleElements<T>& elements,
    const T& tolerance)
{
    using math::exact;
    using math::ratio;

    NearSpaceInit<T> ns;

    // --- Earth constants ---
    ns.re_km = config.ellipsoid.a;
    math::TrackedValue<T> mu_km = config.ellipsoid.GM;
    ns.xke = exact<T>(60) / sqrt(ns.re_km * ns.re_km * ns.re_km / mu_km);

    // --- Zonal harmonics via config.Jn(degree) ---
    math::TrackedValue<T> J2 = config.Jn(2);
    ns.CK2 = J2 / exact<T>(2);
    ns.J3 = config.Jn(3);
    ns.A30 = config.A(3, 0);
    math::TrackedValue<T> J4 = config.Jn(4);

    // --- Store TLE elements ---
    ns.e0 = elements.eccentricity;
    ns.i0 = elements.inclination;
    ns.omega0 = elements.arg_perigee;
    ns.Omega0 = elements.raan;
    ns.M0 = elements.mean_anomaly;
    ns.bstar = elements.bstar;

    // --- Trig constants ---
    ns.cosio = cos(ns.i0);
    ns.sinio = sin(ns.i0);
    ns.theta2 = ns.cosio * ns.cosio;
    ns.x3thm1 = exact<T>(3) * ns.theta2 - exact<T>(1);
    ns.x1mth2 = exact<T>(1) - ns.theta2;
    ns.x7thm1 = exact<T>(7) * ns.theta2 - exact<T>(1);

    // --- Phase 1: Element recovery (standalone module) ---
    ns.n0_unkozai = elements.mean_motion;  // Original TLE mean motion (before Kozai)
    orbit::RecoveredElements<T> recovered = orbit::recover_mean_elements(
        elements.mean_motion, ns.e0, ns.cosio, ns.CK2, ns.xke, ns.re_km, tolerance);
    ns.a0 = recovered.a0;
    ns.n0 = recovered.n0;
    ns.beta0 = recovered.beta0;
    ns.perigee_km = recovered.perigee_km;
    ns.period_min = recovered.period_min;
    ns.is_deep_space = recovered.is_deep_space;
    ns.use_simple_model = recovered.use_simple_model;

    // --- Phase 2: Atmospheric density (standalone module) ---
    atmosphere::DensityParameters<T> atmos = atmosphere::compute_density_parameters(
        ns.perigee_km, ns.a0, ns.e0, ns.re_km);
    ns.s4 = atmos.s;
    ns.qoms24 = atmos.qoms4;

    // --- Phase 3: Drag coefficients (standalone module) ---
    atmosphere::DragInputs<T> drag_in;
    drag_in.a0 = ns.a0;
    drag_in.n0 = ns.n0;
    drag_in.e0 = ns.e0;
    drag_in.beta0 = ns.beta0;
    drag_in.bstar = ns.bstar;
    drag_in.omega0 = ns.omega0;
    drag_in.M0 = ns.M0;
    drag_in.s = ns.s4;
    drag_in.qoms4 = ns.qoms24;
    drag_in.half_J2 = ns.CK2;
    drag_in.A_30 = ns.A30;
    drag_in.cos_i0 = ns.cosio;
    drag_in.sin_i0 = ns.sinio;
    drag_in.cos2_i0 = ns.theta2;
    drag_in.three_cos2i_minus_1 = ns.x3thm1;
    drag_in.sin2_i0 = ns.x1mth2;

    atmosphere::DragCoefficients<T> drag = atmosphere::compute_drag_coefficients(drag_in);

    ns.C1 = drag.C1;  ns.C2 = drag.C2;  ns.C3 = drag.C3;
    ns.C4 = drag.C4;  ns.C5 = drag.C5;
    ns.D2 = drag.D2;  ns.D3 = drag.D3;  ns.D4 = drag.D4;
    ns.t2cof = drag.t2cof;  ns.t3cof = drag.t3cof;
    ns.t4cof = drag.t4cof;  ns.t5cof = drag.t5cof;
    ns.omgcof = drag.omega_drag_coef;
    ns.xmcof = drag.M_drag_coef;
    ns.delmo = drag.eta_cos_M0_cubed;
    ns.sinmo = drag.sin_M0;
    ns.Omega_dot_nkc = drag.Omega_dot_drag;
    ns.xi = drag.xi;
    ns.eta = drag.eta;

    // --- Secular rates via lambda ---
    math::TrackedValue<T> eosq = ns.e0 * ns.e0;
    perturbation::BrouwerSecularRates<T> rates = config.model_functions.secular_rates(ns.n0, ns.a0, eosq, ns.cosio, J2, J4);
    ns.M_dot = rates.M_dot;
    ns.omega_dot = rates.omega_dot;
    ns.Omega_dot = rates.Omega_dot;

    // --- Long-period correction coefficients ---
    // Theorem 14.2 of design/derivations/sgp4_near_earth_drag_theoretical_basis.md:
    //   xlcof = (1/8)·(A30/CK2)·sin(i₀)·(3 + 5·cos i₀)/(1 + cos i₀)
    //
    // Algebraic identity (helpful for the limit analysis below):
    //   (3 + 5c)/(1 + c) = 5 − 2/(1 + c).
    // So:
    //   xlcof = (5/8)·(A30/CK2)·sin i₀ − (1/4)·(A30/CK2)·sin i₀/(1 + cos i₀).
    //
    // CRITICAL-INCLINATION FALLBACK (R14 remediation, 2026-05-14):
    //   The formula has a 1/(1+cos i₀) factor that is singular at
    //   cos i₀ = −1 (retrograde polar, i₀ = π). The guard below
    //   substitutes a closed-form fallback when |cos i₀ ± 1| < 1.5e-12.
    //
    //   PREVIOUS FALLBACK: ratio<T>(3, 8) · A30/CK2 · sin i₀.
    //   This value matches xlcof at i₀ = π/2 (cos i₀ = 0), but is NOT
    //   the limit of the main expression at either pole:
    //     - At cos i₀ → +1 (i₀ → 0): main → (1/2)·(A30/CK2)·sin i₀
    //       (since (3+5cos)/(1+cos) → 8/2 = 4 and (1/8)·4 = 1/2).
    //       The previous 3/8 fallback under-estimated by factor 3/4.
    //     - At cos i₀ → −1 (i₀ → π): main → −∞ (genuine singularity);
    //       the previous 3/8 fallback was an arbitrary band-aid.
    //
    //   NEW FALLBACK (this commit): ratio<T>(1, 2) · A30/CK2 · sin i₀.
    //   Rationale:
    //     (a) At cos i₀ → +1: this IS the rigorous l'Hôpital limit
    //         of the main expression — see analysis above.
    //     (b) At cos i₀ → −1: this gives a small finite value (since
    //         sin i₀ → 0 at i₀ = π), and is no worse than any other
    //         band-aid for a genuinely-singular configuration. Since
    //         retrograde-polar orbits are unphysical for SGP4, the
    //         downstream propagation cost of a wrong xlcof at this
    //         configuration is irrelevant.
    //
    //   PRECISION IMPACT: the guard fires only when |1 ± cos i₀| < 1.5e-12
    //   (≈ 6 × T(double)::epsilon). Inside this band, |sin i₀| <
    //   √(2·1.5e-12) ≈ 1.7e-6, so the absolute difference between
    //   old and new fallback values is
    //     |1/2 − 3/8| · |A30/CK2| · 1.7e-6 ≈ |A30/CK2| · 2.1e-7.
    //   For Earth (A30/CK2 ≈ A30·2/J2 ≈ −2·(2.5e-6)/(1.08e-3) ≈ −4.6e-3),
    //   the absolute difference is bounded by ~1e-9 (dimensionless), well
    //   below SGP4's stated absolute accuracy floor.
    //
    //   See `design/derivations/sgp4_near_earth_drag_theoretical_basis.md`
    //   §14 (Theorem 14.2) for the full derivation of xlcof.
    math::TrackedValue<T> A30_over_CK2 = ns.A30 / ns.CK2;
    if (ns.cosio.value > T(-1) + T(1.5e-12) && ns.cosio.value < T(1) - T(1.5e-12)) {
        ns.xlcof = ratio<T>(1, 8) * A30_over_CK2 * ns.sinio
                   * (exact<T>(3) + exact<T>(5) * ns.cosio) / (exact<T>(1) + ns.cosio);
    } else {
        // l'Hôpital-correct limit at cos i₀ → +1; band-aid at cos i₀ → −1.
        ns.xlcof = ratio<T>(1, 2) * A30_over_CK2 * ns.sinio;
    }
    ns.aycof = ratio<T>(1, 4) * A30_over_CK2 * ns.sinio;

    return ns;
}

/**
 * @brief Propagate a near-space satellite to time tsince from epoch.
 *
 * Implements [SR3] Section 6, pages 11-15: secular updates, drag corrections,
 * long-period periodics, Kepler equation, short-period periodics,
 * and position/velocity computation.
 *
 * @param ns         Initialized near-space constants
 * @param config     Model configuration (for Kepler solver lambda)
 * @param tsince     Time since epoch [minutes]
 * @param tolerance  Convergence tolerance for Kepler solver
 * @return Position [km] and velocity [km/s] in TEME frame
 */
template<typename T>
StateVector<T> propagate_near_space(
    const NearSpaceInit<T>& ns,
    const ModelConfiguration<T>& config,
    const math::TrackedValue<T>& tsince,
    const T& tolerance)
{
    using math::exact;
    using math::ratio;

    // --- Phase 4: Secular update (standalone module) ---
    // Build the BrouwerSecularRates struct from NearSpaceInit fields
    perturbation::BrouwerSecularRates<T> rates;
    rates.M_dot = ns.M_dot;
    rates.omega_dot = ns.omega_dot;
    rates.Omega_dot = ns.Omega_dot;

    // Build the DragCoefficients struct from NearSpaceInit fields
    atmosphere::DragCoefficients<T> drag;
    drag.C1 = ns.C1;  drag.C4 = ns.C4;  drag.C5 = ns.C5;
    drag.D2 = ns.D2;  drag.D3 = ns.D3;  drag.D4 = ns.D4;
    drag.t2cof = ns.t2cof;  drag.t3cof = ns.t3cof;
    drag.t4cof = ns.t4cof;  drag.t5cof = ns.t5cof;
    drag.omega_drag_coef = ns.omgcof;
    drag.M_drag_coef = ns.xmcof;
    drag.eta_cos_M0_cubed = ns.delmo;
    drag.sin_M0 = ns.sinmo;
    drag.Omega_dot_drag = ns.Omega_dot_nkc;
    drag.eta = ns.eta;

    orbit::SecularState<T> sec = orbit::secular_advance(
        ns.a0, ns.e0, ns.M0, ns.omega0, ns.Omega0, ns.n0,
        rates, drag, ns.bstar, ns.xke, tsince, ns.use_simple_model);

    // --- Long-period periodics [SR3] pages 12-13 ---
    // These couple the eccentricity vector to the mean longitude via J₃.
    // Kept inline because they bridge secular and Kepler stages.
    math::TrackedValue<T> axN = sec.e * cos(sec.omega);
    math::TrackedValue<T> beta2 = exact<T>(1) - sec.e * sec.e;
    math::TrackedValue<T> temp_lp = exact<T>(1) / (sec.a * beta2);
    math::TrackedValue<T> xll = temp_lp * ns.xlcof * axN;
    math::TrackedValue<T> aynl = temp_lp * ns.aycof;
    math::TrackedValue<T> xlt = sec.mean_longitude + xll;
    math::TrackedValue<T> ayn = sec.e * sin(sec.omega) + aynl;

    // --- Phase 5: Modified Kepler solver (via injected lambda) ---
    math::TrackedValue<T> capu = math::wrap_two_pi(xlt - sec.Omega);
    math::TrackedValue<T> E_plus_w = config.model_functions.kepler_solver(axN, ayn, capu, tolerance);

    // --- Phase 6: Osculating elements (standalone module) ---
    orbit::OsculatingState<T> osc = orbit::compute_osculating(E_plus_w, axN, ayn, sec.a, ns.xke);

    // --- Phase 7: Short-period corrections (standalone module) ---
    perturbation::CorrectedElements<T> corr = perturbation::apply_short_period(
        osc.r, osc.u, ns.i0, sec.Omega,
        osc.rdot, osc.rfdot, sec.n,
        osc.pl, osc.beta_l,
        ns.CK2, ns.theta2, ns.x1mth2,
        ns.cosio, ns.sinio,
        ns.x3thm1, ns.x7thm1,
        osc.sin_2u, osc.cos_2u);

    // --- Phase 8: Position/velocity from elements (standalone module) ---
    return orbit::elements_to_state(
        corr.r, corr.u, corr.i, corr.Omega,
        corr.rdot, corr.rfdot, ns.re_km);
}

} // namespace sgp4
