#pragma once

/**
 * @file drag_coefficients.h
 * @brief Lane-Hoots atmospheric drag coefficients for SGP4.
 *
 * Computes the drag perturbation coefficients C₁–C₅, D₂–D₄,
 * time polynomial coefficients, and correction terms from the
 * Lane-Hoots (1979) analytical drag theory.
 *
 * These coefficients determine how atmospheric drag modifies the
 * satellite's semi-major axis, eccentricity, and mean longitude
 * as functions of time:
 *   a(t) = a₀ × tempa²     where tempa = 1 − C₁t − D₂t² − D₃t³ − D₄t⁴
 *   e(t) = e₀ − tempe       where tempe = B*C₄t + B*C₅(sin M − sin M₀)
 *   ℓ(t) = ℓ₀ + n₀·templ   where templ = t₂t² + t₃t³ + t⁴(t₄ + t·t₅)
 *
 * @par References
 * - Lane, M.H. and Hoots, F.R. (1979), "General Perturbations Theories
 *   Derived from the 1965 Lane Drag Theory"
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, pages 10-12
 * - design/derivations/012_lane_hoots_drag_derivation.md
 */

#include "../math/tracked_value.h"

namespace atmosphere {

/**
 * @brief All inputs needed for drag coefficient computation.
 */
template<typename T>
struct DragInputs {
    // Recovered orbital elements
    math::TrackedValue<T> a0;      ///< Recovered semi-major axis [Earth radii]
    math::TrackedValue<T> n0;      ///< Recovered mean motion [rad/min]
    math::TrackedValue<T> e0;      ///< Eccentricity at epoch
    math::TrackedValue<T> beta0;   ///< √(1 − e₀²)
    math::TrackedValue<T> bstar;   ///< B* drag coefficient from TLE
    math::TrackedValue<T> omega0;  ///< Argument of perigee at epoch [rad]
    math::TrackedValue<T> M0;      ///< Mean anomaly at epoch [rad]

    // Atmospheric parameters
    math::TrackedValue<T> s;       ///< Atmospheric fitting parameter [Earth radii]
    math::TrackedValue<T> qoms4;   ///< (q₀ − s)⁴

    // Zonal harmonics
    math::TrackedValue<T> half_J2; ///< J₂/2 (= CK2)
    math::TrackedValue<T> A_30;    ///< A(3,0) = −J₃ harmonic coefficient

    // Trig constants
    math::TrackedValue<T> cos_i0;  ///< cos(i₀)
    math::TrackedValue<T> sin_i0;  ///< sin(i₀)
    math::TrackedValue<T> cos2_i0; ///< cos²(i₀)
    math::TrackedValue<T> three_cos2i_minus_1;  ///< 3cos²i₀ − 1
    math::TrackedValue<T> sin2_i0; ///< 1 − cos²i₀ = sin²i₀
};

/**
 * @brief Complete set of drag coefficients and correction terms.
 */
template<typename T>
struct DragCoefficients {
    // Primary drag coefficients [SR3] page 11
    math::TrackedValue<T> C1;  ///< B* × C₂
    math::TrackedValue<T> C2;  ///< Fundamental drag integral
    math::TrackedValue<T> C3;  ///< J₃ eccentricity-drag coupling
    math::TrackedValue<T> C4;  ///< Eccentricity decay rate
    math::TrackedValue<T> C5;  ///< Second-order eccentricity correction

    // Higher-order secular drag polynomial
    math::TrackedValue<T> D2;  ///< 4a₀ξC₁²
    math::TrackedValue<T> D3;  ///< (17a₀+s)·D₂ξC₁/3
    math::TrackedValue<T> D4;  ///< Higher-order term

    // Time polynomial coefficients for mean longitude drag
    math::TrackedValue<T> t2cof;  ///< (3/2)C₁
    math::TrackedValue<T> t3cof;  ///< D₂ + 2C₁²
    math::TrackedValue<T> t4cof;  ///< (1/4)(3D₃ + C₁(12D₂ + 10C₁²))
    math::TrackedValue<T> t5cof;  ///< (1/5)(3D₄ + 12C₁D₃ + 6D₂² + 15C₁²(2D₂+C₁²))

    // Drag correction terms
    math::TrackedValue<T> omega_drag_coef;  ///< B*C₃cos(ω₀) — perigee drag correction
    math::TrackedValue<T> M_drag_coef;      ///< −(2/3)(q₀−s)⁴B*/(e₀η) — mean anomaly correction
    math::TrackedValue<T> eta_cos_M0_cubed; ///< (1 + η cos M₀)³ — reference for M correction
    math::TrackedValue<T> sin_M0;           ///< sin(M₀) — reference for eccentricity correction

    // RAAN drag-gravity coupling
    math::TrackedValue<T> Omega_dot_drag;   ///< −(21/2)n₀k₂cos(i)C₁/(a₀²β₀²)

    // Intermediates needed by propagation
    math::TrackedValue<T> xi;   ///< 1/(a₀ − s)
    math::TrackedValue<T> eta;  ///< a₀·e₀·ξ
};

/**
 * @brief Compute the Lane-Hoots drag coefficients.
 *
 * Every formula is from the Lane-Hoots (1979) analytical drag theory.
 * See design/derivations/012_lane_hoots_drag_derivation.md for the
 * physical derivation of each term.
 *
 * @param in  Drag computation inputs (orbital elements + atmospheric params)
 * @return Complete DragCoefficients with all three errors propagated
 */
template<typename T>
DragCoefficients<T> compute_drag_coefficients(const DragInputs<T>& in) {
    using math::exact;
    using math::ratio;

    DragCoefficients<T> dc;

    // --- Intermediates ---
    dc.xi = exact<T>(1) / (in.a0 - in.s);
    dc.eta = in.a0 * in.e0 * dc.xi;
    math::TrackedValue<T> eta2 = dc.eta * dc.eta;
    math::TrackedValue<T> eeta = in.e0 * dc.eta;
    math::TrackedValue<T> psisq = abs(exact<T>(1) - eta2);  // |1 − η²|
    math::TrackedValue<T> betao2 = exact<T>(1) - in.e0 * in.e0;

    math::TrackedValue<T> coef = in.qoms4 * dc.xi * dc.xi * dc.xi * dc.xi;  // (q₀−s)⁴·ξ⁴
    math::TrackedValue<T> coef1 = coef / (psisq * psisq * psisq * sqrt(psisq)); // ÷ |1−η²|^(7/2)

    // --- C₂: fundamental drag integral ---
    //
    // DERIVATION STATUS: VERIFIED (see design/derivations/020)
    //
    // C₂ = COEF1 × n₀ × [Part A + Part B]
    //
    // Part A (Keplerian drag): a₀(1 + 3/2 η² + eη(4+η²))
    //   VERIFIED: coefficients from Lane (1965) orbit-averaged density integral.
    //   SGP4 simplification drops O(e²) terms (3/4 e² + 3e²η²) from AFGP4.
    //   See 020.Eq.7 (AFGP4 full form) and 020.Eq.12 (SGP4 simplified).
    //   Source: [LH79] p. 25, [SR3] p. 11.
    //
    // Part B (J₂ correction): (3/8)J₂ × ξ/ψ² × (3cos²i−1) × (8+24η²+3η⁴)
    //   VERIFIED: (3/8)J₂ = (3/4)×(J₂/2) = (3/4)×half_J2.
    //   This is (3/2)k₂ × (-1/2+3θ²/2) = (3/2)k₂×A where A is Brouwer's
    //   secular inclination factor. See 020.Eq.13 and 020.Eq.16-17.
    //   SGP4 simplification drops -5eη(4+3η²) cross-coupling from AFGP4.
    //   Source: [LH79] p. 25, [SR3] p. 11.
    //
    // The (8+24η²+3η⁴) polynomial: VERIFIED from orbit-averaged J₂-density
    // coupling integral. In Horner form: 8+3η²(8+η²). See 020.Eq.13.
    //
    dc.C2 = coef1 * in.n0 * (in.a0 * (exact<T>(1) + ratio<T>(3, 2) * eta2
        + eeta * (exact<T>(4) + eta2))
        + ratio<T>(3, 4) * in.half_J2 * dc.xi / psisq * in.three_cos2i_minus_1
          * (exact<T>(8) + exact<T>(3) * eta2 * (exact<T>(8) + eta2)));

    // --- C₁ = B* × C₂ ---
    dc.C1 = in.bstar * dc.C2;

    // --- C₃: J₃ eccentricity-drag coupling ---
    if (in.e0.value > T(1e-4)) {
        dc.C3 = coef * dc.xi * in.A_30 * in.n0 * in.sin_i0 / (in.half_J2 * in.e0);
    } else {
        dc.C3 = exact<T>(0);
    }

    // --- C₄: eccentricity decay rate ---
    dc.C4 = exact<T>(2) * in.n0 * coef1 * in.a0 * betao2
        * (dc.eta * (exact<T>(2) + ratio<T>(1, 2) * eta2)
           + in.e0 * (ratio<T>(1, 2) + exact<T>(2) * eta2)
           - ratio<T>(2, 1) * in.half_J2 * dc.xi / (in.a0 * psisq)
             * (exact<T>(-3) * in.three_cos2i_minus_1
                * (exact<T>(1) - exact<T>(2) * eeta
                   + eta2 * (ratio<T>(3, 2) - ratio<T>(1, 2) * eeta))
                + ratio<T>(3, 4) * in.sin2_i0
                  * (exact<T>(2) * eta2 - eeta * (exact<T>(1) + eta2))
                  * cos(exact<T>(2) * in.omega0)));

    // --- C₅: second-order eccentricity correction ---
    dc.C5 = exact<T>(2) * coef1 * in.a0 * betao2
        * (exact<T>(1) + ratio<T>(11, 4) * (eta2 + eeta) + eeta * eta2);

    // --- D₂, D₃, D₄: higher-order secular drag ---
    dc.D2 = exact<T>(4) * in.a0 * dc.xi * dc.C1 * dc.C1;
    math::TrackedValue<T> temp_d = dc.D2 * dc.xi * dc.C1 / exact<T>(3);
    dc.D3 = (exact<T>(17) * in.a0 + in.s) * temp_d;
    dc.D4 = ratio<T>(1, 2) * temp_d * in.a0 * dc.xi
            * (exact<T>(221) * in.a0 + exact<T>(31) * in.s) * dc.C1;

    // --- Time polynomial coefficients ---
    dc.t2cof = ratio<T>(3, 2) * dc.C1;
    dc.t3cof = dc.D2 + exact<T>(2) * dc.C1 * dc.C1;
    dc.t4cof = ratio<T>(1, 4) * (exact<T>(3) * dc.D3
        + dc.C1 * (exact<T>(12) * dc.D2 + exact<T>(10) * dc.C1 * dc.C1));
    dc.t5cof = ratio<T>(1, 5) * (exact<T>(3) * dc.D4
        + exact<T>(12) * dc.C1 * dc.D3 + exact<T>(6) * dc.D2 * dc.D2
        + exact<T>(15) * dc.C1 * dc.C1 * (exact<T>(2) * dc.D2 + dc.C1 * dc.C1));

    // --- Drag correction terms ---
    dc.omega_drag_coef = in.bstar * dc.C3 * cos(in.omega0);

    if (in.e0.value > T(1e-4)) {
        dc.M_drag_coef = -ratio<T>(2, 3) * coef * in.bstar / eeta;
    } else {
        dc.M_drag_coef = exact<T>(0);
    }

    dc.eta_cos_M0_cubed = (exact<T>(1) + dc.eta * cos(in.M0));
    dc.eta_cos_M0_cubed = dc.eta_cos_M0_cubed * dc.eta_cos_M0_cubed * dc.eta_cos_M0_cubed;

    dc.sin_M0 = sin(in.M0);

    // --- RAAN drag-gravity coupling ---
    dc.Omega_dot_drag = ratio<T>(-21, 2) * in.n0 * in.half_J2 * in.cos_i0
                        / (in.a0 * in.a0 * in.beta0 * in.beta0) * dc.C1;

    return dc;
}

} // namespace atmosphere
