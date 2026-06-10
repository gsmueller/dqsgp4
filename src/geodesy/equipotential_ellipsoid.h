#pragma once

/**
 * @file equipotential_ellipsoid.h
 * @brief Equipotential (level) ellipsoid and all derived constants.
 *
 * Constructed from four defining parameters. Supports TWO initialization paths:
 *
 * **WGS84 path:** \f$(a, 1/f, GM, \omega)\f$ — flattening is defining, \f$J_2\f$ is derived.
 * Used when \f$1/f\f$ is the fourth defining parameter.
 * \f$e^2 = 2f - f^2\f$ (direct, no iteration).
 *
 * **GRS80/WGS72 path:** \f$(a, J_2, GM, \omega)\f$ — \f$J_2\f$ is defining, \f$f\f$ is derived.
 * Used when \f$J_2\f$ is the fourth defining parameter (as in GRS80 and the SGP4 WGS72 model).
 * \f$e^2\f$ is solved iteratively from \f$J_2\f$ via the Brouwer formula.
 *
 * Each input is a TrackedValue<T> carrying its three errors.
 * All derived constants are computed at construction time and stored.
 * No magic numbers. The theory note is
 * design/derivations/ch14_equipotential_ellipsoid.md (geometric chain
 * §14.3; the q0/q0' series §14.4.2-3; Brouwer's J2<->e² §14.5; even
 * zonals §14.6; Somigliana §14.7). The "Derivation 001-006" labels in
 * the member comments are this file's internal numbering of that chain.
 *
 * @par References
 * - NGA.STND.0036 (2014), Chapter 3 and Appendix B
 * - Moritz, H. (1980), "Geodetic Reference System 1980", Bulletin Geodesique 54(3)
 * - Heiskanen, W.A. and Moritz, H. (1967), "Physical Geodesy", Sections 2-7 to 2-10
 *
 * @par Key Finding: WGS72 J₂
 * The SGP4 WGS72 gravity model specifies J₂ = 0.001082616 as a direct input,
 * NOT derived from 1/f. Using the WGS84-style initialization (from 1/f = 298.26)
 * produces a slightly different J₂ (0.001082611) due to the approximations
 * in the derivation chain. For SGP4 compatibility, the GRS80 initialization
 * path must be used when J₂ is the specified defining parameter.
 *
 * The series tolerance is a bare T — it controls how many series terms
 * are evaluated. It affects only the precision error of the results.
 */

#include "../math/tracked_value.h"
#include "../math/series.h"
#include "../math/factorial.h"
#include <boost/math/constants/constants.hpp>

namespace geodesy {

/// The level (equipotential) ellipsoid: four defining parameters in, every
/// derived geometric/physical constant out (computed once at construction,
/// each a TrackedValue with the series truncation in its budget).
template<typename T>
class EquipotentialEllipsoid {
public:
    // --- The four defining parameters (stored as given) ---
    math::TrackedValue<T> a;       // semi-major axis [m]
    math::TrackedValue<T> inv_f;   // reciprocal flattening
    math::TrackedValue<T> GM;      // geocentric gravitational constant [m³/s²]
    math::TrackedValue<T> omega;   // angular velocity [rad/s]

    // --- Geometric constants (Derivation 003) ---
    math::TrackedValue<T> f;       // flattening
    math::TrackedValue<T> e2;      // first eccentricity squared
    math::TrackedValue<T> e;       // first eccentricity
    math::TrackedValue<T> e_prime; // second eccentricity
    math::TrackedValue<T> ep2;     // second eccentricity squared
    math::TrackedValue<T> b;       // semi-minor axis
    math::TrackedValue<T> E_lin;   // linear eccentricity
    math::TrackedValue<T> c;       // polar radius of curvature

    // --- Series-evaluated intermediate terms (Derivations 001-002) ---
    math::TrackedValue<T> q0;      // from q0 series
    math::TrackedValue<T> q0p;     // q0' from q0' series

    // --- Physical constants (Derivation 004) ---
    math::TrackedValue<T> m_const; // normal gravity formula constant (ω²a²b/GM)
    math::TrackedValue<T> gamma_e; // normal gravity at equator
    math::TrackedValue<T> gamma_p; // normal gravity at pole
    math::TrackedValue<T> k_som;   // Somigliana constant
    math::TrackedValue<T> U0;      // normal potential on ellipsoid

    // --- Radii (Derivation 005) ---
    math::TrackedValue<T> R1;      // mean radius of three semi-axes
    math::TrackedValue<T> R3;      // radius of sphere of equal volume

    // --- Constructors ---

    /**
     * @brief Construct from (a, 1/f, GM, ω) — the WGS84 initialization path.
     * @param a_in      Semi-major axis
     * @param inv_f_in  Reciprocal flattening (1/f)
     * @param GM_in     Geocentric gravitational constant
     * @param omega_in  Angular velocity
     * @param series_tolerance  Tolerance for series evaluations
     */
    EquipotentialEllipsoid(
        const math::TrackedValue<T>& a_in,
        const math::TrackedValue<T>& inv_f_in,
        const math::TrackedValue<T>& GM_in,
        const math::TrackedValue<T>& omega_in,
        const T& series_tolerance)
        : a(a_in), inv_f(inv_f_in), GM(GM_in), omega(omega_in)
    {
        using math::exact;
        using math::ratio;

        // --- Geometric chain (Derivation 003) ---
        f = exact<T>(1) / inv_f;
        e2 = exact<T>(2) * f - f * f;
        e = sqrt(e2);
        b = a * (exact<T>(1) - f);
        E_lin = a * e;
        c = a * a / b;
        ep2 = e2 / (exact<T>(1) - e2);
        e_prime = sqrt(ep2);

        // --- q0 series (Derivation 001) ---
        // 2*q0 = sum_{n=1}^{inf} (-1)^{n+1} * 4n/((2n+1)(2n+3)) * e'^{2n+1}
        std::function<math::TrackedValue<T>(int)> q0_term = [&](int n) -> math::TrackedValue<T> {
            int sign = (n % 2 == 1) ? 1 : -1; // (-1)^{n+1}
            math::TrackedValue<T> coeff = ratio<T>(sign * 4 * n, (2*n+1) * (2*n+3));
            // e'^{2n+1} computed by repeated multiplication (no pow)
            math::TrackedValue<T> ep_power = e_prime;
            for (int i = 1; i < 2*n+1; ++i) {
                ep_power = ep_power * e_prime;
            }
            return coeff * ep_power;
        };

        math::TrackedValue<T> two_q0 = math::alternating_series<T>(1, q0_term, series_tolerance);
        q0 = two_q0 / exact<T>(2);

        // --- q0' series (Derivation 002) ---
        // q0' = sum_{n=1}^{inf} (-1)^{n+1} * 6/((2n+1)(2n+3)) * e'^{2n}
        std::function<math::TrackedValue<T>(int)> q0p_term = [&](int n) -> math::TrackedValue<T> {
            int sign = (n % 2 == 1) ? 1 : -1;
            math::TrackedValue<T> coeff = ratio<T>(sign * 6, (2*n+1) * (2*n+3));
            math::TrackedValue<T> ep_power = exact<T>(1);
            for (int i = 0; i < 2*n; ++i) {
                ep_power = ep_power * e_prime;
            }
            return coeff * ep_power;
        };

        q0p = math::alternating_series<T>(1, q0p_term, series_tolerance);

        // --- Physical chain (Derivation 004) ---
        m_const = omega * omega * a * a * b / GM;

        // gamma_e = (GM/(a*b)) * (1 - m - m*e'*q0'/(6*q0))
        math::TrackedValue<T> me_q0p_over_6q0 = m_const * e_prime * q0p / (exact<T>(6) * q0);
        gamma_e = (GM / (a * b)) * (exact<T>(1) - m_const - me_q0p_over_6q0);

        // gamma_p = (GM/a²) * (1 + m*e'*q0'/(3*q0))
        math::TrackedValue<T> me_q0p_over_3q0 = m_const * e_prime * q0p / (exact<T>(3) * q0);
        gamma_p = (GM / (a * a)) * (exact<T>(1) + me_q0p_over_3q0);

        // k = b*gamma_p/(a*gamma_e) - 1
        k_som = b * gamma_p / (a * gamma_e) - exact<T>(1);

        // U0 via series: (GM/b)*(1 - e'^2/3 + e'^4/5 - ...) + omega^2*a^2/3
        // (Derivation 004, Step 5)
        std::function<math::TrackedValue<T>(int)> u0_term = [&](int n) -> math::TrackedValue<T> {
            int sign = (n % 2 == 0) ? 1 : -1;
            math::TrackedValue<T> coeff = ratio<T>(sign, 2*n+1);
            math::TrackedValue<T> ep_power = exact<T>(1);
            for (int i = 0; i < 2*n; ++i) {
                ep_power = ep_power * e_prime;
            }
            return coeff * ep_power;
        };

        math::TrackedValue<T> u0_series = math::alternating_series<T>(0, u0_term, series_tolerance);
        U0 = GM / b * u0_series + omega * omega * a * a / exact<T>(3);

        // --- Radii ---
        R1 = a * (exact<T>(1) - f / exact<T>(3));

        // R3 = cbrt(a²·b) = a·cbrt(1-f). Single-sourced via the tracked cbrt
        // primitive (B2), retiring the open-coded Newton cube-root.
        R3 = cbrt(a * a * b);
    }

    // --- Named earth-model presets (C2) ---
    // Canonical catalog of standard reference ellipsoids at the geodesy layer.
    // The GEOMETRIC/KINEMATIC defining parameters (a, 1/f, omega) are exact by
    // convention -> TrackedValue::defined (a metre and a flattening are defined,
    // not truncated). The GRAVITATIONAL FIELD quantities are NOT conventions:
    //   - J2 is the coefficient of a truncated spherical-harmonic series (and,
    //     for the normal field, the truncated ellipsoid J2n series itself), so it
    //     carries a model-fidelity accuracy -> model_coefficient (no published σ)
    //     or measured (with σ). It is NEVER defined(): user directive 2026-06-06
    //     "each J_k is the result of a series truncation".
    //   - GM is a measured physical quantity -> measured (with σ) or
    //     model_coefficient (adopted, no σ). Never defined().
    // These are the reusable single source; consumers should delegate here
    // rather than re-typing the numbers (the consumer de-duplication is item F3).

    /// WGS 84 (NGA.STND.0036): a and 1/f definitional, GM measured (±8e5).
    static EquipotentialEllipsoid wgs84(const T& series_tolerance) {
        using TV = math::TrackedValue<T>;
        return EquipotentialEllipsoid(
            TV::defined("6378137.0"), TV::defined("298.257223563"),
            TV::measured("3.986004418e14", "8e5"),
            TV::defined("7.2921151467e-5"), series_tolerance);
    }

    /// GRS 80 (Moritz 1980, Bull. Géodésique 54). a, omega defining (convention);
    /// J₂, GM are gravity-field quantities carrying accuracy -> model_coefficient.
    static EquipotentialEllipsoid grs80(const T& series_tolerance) {
        using TV = math::TrackedValue<T>;
        return from_J2(
            TV::defined("6378137.0"), TV::model_coefficient("1.08263e-3"),
            TV::model_coefficient("3.986005e14"), TV::defined("7.292115e-5"), series_tolerance);
    }

    /// WGS 72 (Spacetrack Report #3). Uses the from_J2 path so J₂ matches the SGP4
    /// model exactly (see "Key Finding: WGS72 J₂" above). a, omega defining
    /// (convention); J₂, GM are gravity-field quantities -> model_coefficient
    /// (each a truncated-series/fit result that must carry accuracy).
    static EquipotentialEllipsoid wgs72(const T& series_tolerance) {
        using TV = math::TrackedValue<T>;
        return from_J2(
            TV::defined("6378135.0"), TV::model_coefficient("1.082616e-3"),
            TV::model_coefficient("3.986008e14"), TV::defined("7.2921151467e-5"), series_tolerance);
    }

    /// IERS 2010 / IAU 2009 NSFA current best estimates (TT-compatible): the
    /// equatorial radius, GM, and J₂ are measured quantities with published
    /// uncertainties. Source: IAU 2009 NSFA (iau-a3.gitlab.io/NSFA).
    static EquipotentialEllipsoid iers2010(const T& series_tolerance) {
        using TV = math::TrackedValue<T>;
        return from_J2(
            TV::measured("6378136.6", "0.1"),
            TV::measured("1.0826359e-3", "1e-10"),
            TV::measured("3.986004418e14", "8e5"),
            TV::defined("7.292115e-5"), series_tolerance);
    }

    /**
     * @brief Construct from (a, J₂, GM, ω) — the GRS80/WGS72 initialization path.
     *
     * When J₂ is the defining parameter (as in GRS80 and the SGP4 WGS72 model),
     * e² must be solved iteratively from J₂ via the Brouwer formula:
     * \f[
     *   e^2 = 3J_2 + \frac{4}{15}\frac{\omega^2 a^3}{GM}\frac{a^3}{2q_0}
     * \f]
     *
     * This factory method performs the iteration and constructs the ellipsoid.
     *
     * @par Key Finding
     * The SGP4 WGS72 model specifies J₂ = 0.001082616 directly. Using the
     * WGS84-style constructor with 1/f = 298.26 gives J₂ = 0.001082611 —
     * a discrepancy of 4.8E-9. This factory avoids that discrepancy by
     * using J₂ as the defining parameter and deriving f from it.
     *
     * @param a_in      Semi-major axis
     * @param J2_in     Second zonal harmonic (defining parameter)
     * @param GM_in     Geocentric gravitational constant
     * @param omega_in  Angular velocity
     * @param series_tolerance  Tolerance for series evaluations
     */
    static EquipotentialEllipsoid from_J2(
        const math::TrackedValue<T>& a_in,
        const math::TrackedValue<T>& J2_in,
        const math::TrackedValue<T>& GM_in,
        const math::TrackedValue<T>& omega_in,
        const T& series_tolerance)
    {
        using math::exact;

        // Iterative solution for e² from J₂.
        // Starting guess: e² ≈ 3*J₂ (first-order approximation)
        math::TrackedValue<T> e2_guess = exact<T>(3) * J2_in;

        // Iterate: e² = 3*J₂ + (4/15)*(ω²a³/GM)*(a³/(2*q₀))
        // where q₀ depends on e' which depends on e²
        for (int iter = 0; iter < 20; ++iter) {
            // Compute e' from current e² guess
            math::TrackedValue<T> e_prime_iter = sqrt(e2_guess / (exact<T>(1) - e2_guess));

            // Compute q₀ from e' (series form from Derivation 001)
            std::function<math::TrackedValue<T>(int)> q0_term = [&](int n) -> math::TrackedValue<T> {
                int sign = (n % 2 == 1) ? 1 : -1;
                math::TrackedValue<T> coeff = math::ratio<T>(sign * 4 * n, (2*n+1) * (2*n+3));
                math::TrackedValue<T> ep_power = e_prime_iter;
                for (int i = 1; i < 2*n+1; ++i) ep_power = ep_power * e_prime_iter;
                return coeff * ep_power;
            };
            math::TrackedValue<T> two_q0 = math::alternating_series<T>(1, q0_term, series_tolerance);
            math::TrackedValue<T> q0_iter = two_q0 / exact<T>(2);

            // m = ω²a²b/GM, but b depends on e². Use b ≈ a*sqrt(1-e²)
            math::TrackedValue<T> b_iter = a_in * sqrt(exact<T>(1) - e2_guess);
            math::TrackedValue<T> m_iter = omega_in * omega_in * a_in * a_in * b_iter / GM_in;

            // New e² estimate
            math::TrackedValue<T> e2_new = exact<T>(3) * J2_in
                + exact<T>(2) * m_iter * e_prime_iter * e2_guess / (exact<T>(15) * q0_iter);

            using std::abs;
            T correction = abs((e2_new - e2_guess).value);
            e2_guess = e2_new;

            if (correction < series_tolerance) {
                // REQ-EF-5: convergence residual added to precision (R02).
                // The Banach fixed-point iteration on 1/f via J2 converges; the
                // final correction magnitude is the rigorous residual bound.
                // Pattern matches kepler.h:79 and series.h:147.
                e2_guess.errors.precision = e2_guess.errors.precision + correction;
                break;
            }
        }

        // Derive 1/f from e²: f = 1 - sqrt(1-e²), so 1/f = 1/f
        math::TrackedValue<T> f_derived = exact<T>(1) - sqrt(exact<T>(1) - e2_guess);
        math::TrackedValue<T> inv_f_derived = exact<T>(1) / f_derived;

        // Construct using the standard (1/f) constructor with the derived value
        return EquipotentialEllipsoid(a_in, inv_f_derived, GM_in, omega_in, series_tolerance);
    }

    // --- Methods ---

    /// Normal gravity at geodetic latitude phi (Somigliana formula)
    /// gamma = gamma_e * (1 + k*sin²φ) / sqrt(1 - e²*sin²φ)
    math::TrackedValue<T> normal_gravity(const math::TrackedValue<T>& phi) const {
        using math::exact;
        math::TrackedValue<T> sin_phi = sin(phi);
        math::TrackedValue<T> sin2_phi = sin_phi * sin_phi;
        math::TrackedValue<T> numerator = exact<T>(1) + k_som * sin2_phi;
        math::TrackedValue<T> denominator = sqrt(exact<T>(1) - e2 * sin2_phi);
        return gamma_e * numerator / denominator;
    }

    /// Even zonal harmonic J_{2n} (Derivation 006)
    /// J_{2n} = (-1)^{n+1} * 3*e^{2n} / ((2n+1)(2n+3)) * (1 - n + 5n*J2/e²)
    math::TrackedValue<T> J2n(int n) const {
        using math::exact; using math::ratio;

        // First compute J2 (n=1 case requires special handling to avoid circular ref)
        // J2 = (e²/3)*(1 - 2*m*e'/(15*q0))
        math::TrackedValue<T> J2 = (e2 / exact<T>(3))
            * (exact<T>(1) - exact<T>(2) * m_const * e_prime / (exact<T>(15) * q0));

        if (n == 1) return J2;

        int sign = (n % 2 == 1) ? 1 : -1; // (-1)^{n+1}

        // e^{2n} by repeated multiplication
        math::TrackedValue<T> e2n = exact<T>(1);
        for (int i = 0; i < n; ++i) {
            e2n = e2n * e2;
        }

        math::TrackedValue<T> prefactor = ratio<T>(sign * 3, (2*n+1) * (2*n+3));
        math::TrackedValue<T> inner = exact<T>(1) - exact<T>(n) + exact<T>(5*n) * J2 / e2;

        return prefactor * e2n * inner;
    }
};

} // namespace geodesy
