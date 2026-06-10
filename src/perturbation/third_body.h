#pragma once

/**
 * @file third_body.h
 * @brief Third-body gravitational perturbation (Sun or Moon).
 *
 * Computes the secular rate corrections to a satellite's orbital elements
 * from the gravitational attraction of a third body. The same function
 * works for both the Sun and Moon — just pass different ephemeris inputs.
 *
 * The disturbing function is expanded to second order in (r/r₃) using
 * the Legendre polynomial P₂(cos S), where S is the angle between the
 * satellite and the third body as seen from Earth. This is orbit-averaged
 * over one satellite revolution, then differentiated via Lagrange's
 * planetary equations to extract secular rates.
 *
 * The perturbation strength scales as μ₃/r₃³. The Moon (~8.6e-14 /s²)
 * is ~2.2× stronger than the Sun (~3.9e-14 /s²) because proximity
 * wins over mass for the tidal (∝ 1/r³) perturbation.
 *
 * @par References
 * - Brouwer, D. and Clemence, G. (1961), "Methods of Celestial Mechanics"
 * - design/derivations/017_third_body_perturbation.md
 */

#include "../math/tracked_value.h"

namespace perturbation {

/**
 * @brief Secular rate corrections from one third body.
 *
 * Five rates: three angular (Ω, ω, M) plus eccentricity and inclination
 * changes that are zero in the near-earth Brouwer-only model but nonzero
 * when lunar/solar perturbations act on the orbit.
 */
template<typename T>
struct ThirdBodyRates {
    math::TrackedValue<T> dM_dt;      ///< Mean anomaly rate correction [rad/min]
    math::TrackedValue<T> domega_dt;  ///< Argument of perigee rate correction [rad/min]
    math::TrackedValue<T> dOmega_dt;  ///< RAAN rate correction [rad/min]
    math::TrackedValue<T> de_dt;      ///< Eccentricity rate [1/min]
    math::TrackedValue<T> di_dt;      ///< Inclination rate [rad/min]
};

/**
 * @brief Compute third-body secular rate corrections.
 *
 * This implements the orbit-averaged P₂(cos S) perturbation from a
 * third body at known ecliptic position. The function is general —
 * it works for any third body given its position and gravitational
 * parameter.
 *
 * @param n             Satellite mean motion [rad/min]
 * @param a             Satellite semi-major axis [Earth radii]
 * @param e             Satellite eccentricity
 * @param i             Satellite inclination [rad]
 * @param omega         Satellite argument of perigee [rad]
 * @param Omega         Satellite RAAN [rad]
 * @param body_sin_lon  sin(ecliptic longitude of third body)
 * @param body_cos_lon  cos(ecliptic longitude of third body)
 * @param body_sin_lat  sin(ecliptic latitude of third body) (0 for Sun)
 * @param body_cos_lat  cos(ecliptic latitude of third body) (1 for Sun)
 * @param sin_obliquity sin(ε) — obliquity of the ecliptic
 * @param cos_obliquity cos(ε) — obliquity of the ecliptic
 * @param perturbation_coef  Perturbation coefficient incorporating μ₃/r₃³
 *                           and the orbit-averaging factor
 * @return ThirdBodyRates with all five secular rate corrections
 */
template<typename T>
ThirdBodyRates<T> compute_third_body_rates(
    const math::TrackedValue<T>& n,
    const math::TrackedValue<T>& a,
    const math::TrackedValue<T>& e,
    const math::TrackedValue<T>& i,
    const math::TrackedValue<T>& omega,
    const math::TrackedValue<T>& Omega,
    const math::TrackedValue<T>& body_sin_lon,
    const math::TrackedValue<T>& body_cos_lon,
    const math::TrackedValue<T>& body_sin_lat,
    const math::TrackedValue<T>& body_cos_lat,
    const math::TrackedValue<T>& sin_obliquity,
    const math::TrackedValue<T>& cos_obliquity,
    const math::TrackedValue<T>& perturbation_coef)
{
    using math::exact;
    using math::ratio;

    ThirdBodyRates<T> rates;

    // --- Third body direction in equatorial coordinates ---
    // Transform from ecliptic (λ, β) to equatorial (α, δ) using obliquity ε.
    //
    // The unit vector toward the third body in equatorial coords:
    //   x₃ = cos(β)cos(λ)
    //   y₃ = cos(β)sin(λ)cos(ε) - sin(β)sin(ε)
    //   z₃ = cos(β)sin(λ)sin(ε) + sin(β)cos(ε)
    math::TrackedValue<T> x3 = body_cos_lat * body_cos_lon;
    math::TrackedValue<T> y3 = body_cos_lat * body_sin_lon * cos_obliquity - body_sin_lat * sin_obliquity;
    math::TrackedValue<T> z3 = body_cos_lat * body_sin_lon * sin_obliquity + body_sin_lat * cos_obliquity;

    // --- Project third body onto satellite's orbital frame ---
    // The satellite's node direction: (cos Ω, sin Ω, 0)
    // The satellite's ascending perpendicular in the orbital plane
    math::TrackedValue<T> sin_Omega = sin(Omega);
    math::TrackedValue<T> cos_Omega = cos(Omega);
    math::TrackedValue<T> sin_i = sin(i);
    math::TrackedValue<T> cos_i = cos(i);

    // Components of the third body in the satellite's orbital reference:
    // a₁ = x₃·cos(Ω) + y₃·sin(Ω)
    // a₃ = z₃/sin(i)  [perpendicular to orbital plane]
    // a₂ depends on both i and Ω
    math::TrackedValue<T> a1 = x3 * cos_Omega + y3 * sin_Omega;
    math::TrackedValue<T> a3_sin_i = z3;  // a₃·sin(i) = z₃ (avoids division by sin(i) at equatorial)
    math::TrackedValue<T> a2 = -x3 * sin_Omega * cos_i + y3 * cos_Omega * cos_i + z3 * sin_i;

    // --- Orbit-averaged P₂(cos S) perturbation ---
    // After averaging over one satellite revolution, the disturbing function
    // depends on a₁, a₂, a₃ and the satellite's eccentricity and argument
    // of perigee. The algebra produces terms involving:
    //   a₁² + a₂² (in-plane squared), a₃² (out-of-plane squared)
    //   a₁·cos(ω), a₂·sin(ω), etc.

    // For the secular rates, the orbit-averaged contributions are:
    // (These are the standard Kaula-type expansion terms for the P₂ case)

    math::TrackedValue<T> sin_omega = sin(omega);
    math::TrackedValue<T> cos_omega = cos(omega);
    math::TrackedValue<T> e2 = e * e;
    math::TrackedValue<T> beta = sqrt(exact<T>(1) - e2);

    // In-plane and out-of-plane projections
    math::TrackedValue<T> a1_sq = a1 * a1;
    math::TrackedValue<T> a2_sq = a2 * a2;
    math::TrackedValue<T> a12_sq = a1_sq + a2_sq;  // = 1 - a3² for a unit vector

    // Cross terms with argument of perigee
    math::TrackedValue<T> a1_cos_w = a1 * cos_omega;
    math::TrackedValue<T> a2_sin_w = a2 * sin_omega;
    math::TrackedValue<T> a1_sin_w = a1 * sin_omega;
    math::TrackedValue<T> a2_cos_w = a2 * cos_omega;

    // The orbit-averaged secular perturbation (from Brouwer & Clemence):
    // These formulas give the time derivatives of the orbital elements.
    // The perturbation_coef absorbs μ₃/(r₃³ · n · a²) and averaging factors.

    math::TrackedValue<T> coef = perturbation_coef;

    // dΩ/dt: RAAN precession from third body
    // Proportional to the out-of-plane component
    rates.dOmega_dt = coef * (a1_cos_w + a2_sin_w) * a3_sin_i
                      / (n * a * a * beta);

    // dω/dt: apsidal rate from third body
    rates.domega_dt = coef * (
        -a1_sin_w + a2_cos_w
        - cos_i * (a1_cos_w + a2_sin_w) * a3_sin_i / (beta * sin_i * sin_i)
    ) / (n * a * a * beta);

    // de/dt: eccentricity change
    rates.de_dt = coef * beta * (a1_sin_w - a2_cos_w)
                  / (n * a * a);

    // di/dt: inclination change
    rates.di_dt = coef * (a1_cos_w + a2_sin_w)
                  * (cos_i * a3_sin_i / sin_i - exact<T>(1))
                  / (n * a * a * beta);

    // dM/dt: mean anomaly rate correction (from da/dt and direct contribution)
    rates.dM_dt = -coef * (exact<T>(3) * (a12_sq - ratio<T>(1, 3))
                  + exact<T>(3) * e * (a1_cos_w + a2_sin_w))
                  / (exact<T>(2) * n * a);

    return rates;
}

} // namespace perturbation
