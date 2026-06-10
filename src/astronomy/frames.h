#pragma once

/**
 * @file frames.h
 * @brief Reference-frame rotation primitives and named frame transforms (L2).
 *
 * A reference frame is an orthonormal triad; relating two is a proper rotation
 * (∈ SO(3), det = +1). This module provides the elementary body-axis rotations
 * as Matrix3<T> (math/matrix3.h, from H1) and the named transforms that compose
 * out of them. Every rotation is exact-orthonormal by construction, composable
 * (matrix product), and invertible by transpose; it carries its angle's
 * TrackedValue three-error budget through cos/sin.
 *
 * Convention (right-handed, active about a body axis), as fixed in the theory
 * note §2:
 *
 *         ⎡1    0      0  ⎤        ⎡ cosθ 0 −sinθ⎤        ⎡ cosθ  sinθ 0⎤
 *  Rx(θ)= ⎢0  cosθ  sinθ ⎥  Ry(θ)=⎢ 0    1   0  ⎥  Rz(θ)=⎢−sinθ  cosθ 0⎥
 *         ⎣0 −sinθ  cosθ ⎦        ⎣ sinθ 0  cosθ⎦        ⎣  0     0   1⎦
 *
 * Each is orthonormal (RᵀR = I), det = +1, R(θ)⁻¹ = R(−θ) = R(θ)ᵀ. Frame chains
 * are products; the inverse of a chain is the reversed product of transposes.
 *
 * @par Sidereal rotation (the one transform consumed today)
 * The SGP4 output is TEME (True Equator, Mean Equinox of date); the TEME→ECEF
 * rotation is taken as Rz(GMST) — the standard SGP4 convention (no equation of
 * equinoxes, no polar motion). GMST is supplied by L1 (sidereal_time.h). This
 * is exactly the rotation forces/gravity_tesseral.h hand-coded inline; it is now
 * the named, identity-tested sidereal_rotation primitive, value-preserving.
 *
 * @par L3 frame chain (landed here) and the R4b Earth-fixed chain (landed elsewhere)
 * The IAU 2006 bias-precession P(t) (precession_iau2006), the obliquity rotation
 * Rx(−ε_A) (obliquity_rotation, which wires the obliquity_iau2006 generator), and
 * the composite ecliptic-of-date → GCRS map (ecliptic_to_gcrs = Pᵀ·Rx(−ε_A)) live
 * here with their ERFA oracle — gated element-wise vs erfa.pmat06/ecm06/obl06
 * (test_frame_chain). See design/derivations/frame_chain.md. The once-deferred
 * pieces landed with their consumer (R4b, gate NUT1): nutation N(t) is
 * astronomy/nutation.h (IAU 2000A, the in-repo SOFA table) and polar motion W +
 * ERA/GMST06/GAST + the full GCRS→ITRS map are astronomy/earth_rotation.h —
 * erfa-arbitrated (nut06a/gst06a/pom00/c2t06a), theory in
 * design/derivations/nutation_itrs.md. Only a Frame enum remains unbuilt (no
 * consumer has needed a tagged-frame type).
 *
 * @par Verification (no perceived fidelity)
 * The primitives are proved against analytic identities (orthonormality,
 * det = +1, round-trip R(θ)R(−θ) = I, R(0) = I, axis images at π/2). The
 * sidereal rotation reproduces the prior gravity_tesseral.h Rz arithmetic
 * bit-for-bit (the value-preserving migration). The GMST value itself was
 * oracle-verified at L1. So this module introduces no un-oracled fidelity.
 * Gate: tests/test_reference_frames (ExeGate FRAME1).
 *
 * @par References (born-digital) — see design/derivations/reference_frames.md
 * - IERS Conventions (2010), Ch. 5 — the rotation sequence W·R·N·P.
 * - Vallado, Fundamentals of Astrodynamics, Ch. 3 — TEME, the SGP4 convention.
 * - ERFA (IAU SOFA) / astropy coordinates — the L3 verification oracle.
 */

#include "obliquity.h"                   // obliquity_iau2006 (ε_A) — wired here (L3)
#include "../math/angles.h"              // pi (arcsec → rad)
#include "../math/matrix3.h"
#include "../math/tracked_polynomial.h"  // FW precession angle polynomials
#include "../math/tracked_value.h"

#include <vector>

namespace astronomy {

/// Elementary rotation about the body x-axis (theory note §2). Exact
/// orthonormal; the off-diagonal/zero entries default to exact zero and the
/// (0,0) entry to exact one, so only cos/sin slots carry an error budget.
template<typename T>
math::Matrix3<T> rot_x(const math::TrackedValue<T>& a) {
    math::Matrix3<T> R;                 // all entries default to exact zero
    math::TrackedValue<T> c = cos(a);   // ADL: math::cos(TrackedValue)
    math::TrackedValue<T> s = sin(a);
    R.m[0][0] = math::exact<T>(1);
    R.m[1][1] = c;   R.m[1][2] = s;
    R.m[2][1] = -s;  R.m[2][2] = c;
    return R;
}

/// Elementary rotation about the body y-axis (theory note §2).
template<typename T>
math::Matrix3<T> rot_y(const math::TrackedValue<T>& a) {
    math::Matrix3<T> R;
    math::TrackedValue<T> c = cos(a);
    math::TrackedValue<T> s = sin(a);
    R.m[0][0] = c;   R.m[0][2] = -s;
    R.m[1][1] = math::exact<T>(1);
    R.m[2][0] = s;   R.m[2][2] = c;
    return R;
}

/// Elementary rotation about the body z-axis (theory note §2). This is the
/// primitive the sidereal rotation is built from; its [0][1] = +sinθ /
/// [1][0] = −sinθ sign placement is what reproduces the inline tesseral Rz
/// arithmetic bit-for-bit (see sidereal_rotation).
template<typename T>
math::Matrix3<T> rot_z(const math::TrackedValue<T>& a) {
    math::Matrix3<T> R;
    math::TrackedValue<T> c = cos(a);
    math::TrackedValue<T> s = sin(a);
    R.m[0][0] = c;   R.m[0][1] = s;
    R.m[1][0] = -s;  R.m[1][1] = c;
    R.m[2][2] = math::exact<T>(1);
    return R;
}

/// TEME→ECEF (PEF) sidereal rotation R = Rz(GMST), the SGP4 convention (no
/// equation of equinoxes, no polar motion). The inverse (ECEF→TEME) is its
/// transpose. GMST [rad] comes from astronomy::compute_gmst (L1, sidereal_time.h).
/// Applied to a position vector this is BIT-IDENTICAL to the prior hand-coded
/// componentwise rotation in forces/gravity_tesseral.h (the +0·z terms are exact
/// no-ops, 1·z and a+(-b) are IEEE-exact), so the migration is value-preserving.
template<typename T>
math::Matrix3<T> sidereal_rotation(const math::TrackedValue<T>& gmst) {
    return rot_z(gmst);
}

namespace detail {

/// Build an IAU-2006 secular angle polynomial [arcsec] from its born-digital
/// coefficient strings and conservative truncation magnitudes (frame_chain.md §4;
/// mirrors obliquity.h). Each c_k is a model_coefficient (digits → accuracy,
/// binary storage → T-scaling precision).
template<typename T>
math::TrackedPolynomial<T> fw_poly(const char* const c[6], const double m[6]) {
    using TV = math::TrackedValue<T>;
    return math::TrackedPolynomial<T>(
        std::vector<TV>{TV::model_coefficient(c[0]), TV::model_coefficient(c[1]),
                        TV::model_coefficient(c[2]), TV::model_coefficient(c[3]),
                        TV::model_coefficient(c[4]), TV::model_coefficient(c[5])},
        std::vector<T>{T(m[0]), T(m[1]), T(m[2]), T(m[3]), T(m[4]), T(m[5])});
}

} // namespace detail

/// IAU 2006 bias-precession matrix P(t) (GCRS → mean equatorial of date), t = TT
/// Julian centuries from J2000.0. Built from the Fukushima-Williams angles
/// (γ̄, φ̄, ψ̄, ε_A) via fw2m = R1(−ε_A)·R3(−ψ̄)·R1(φ̄)·R3(γ̄) (frame_chain.md §4;
/// born_digital_sources.md §A, VERIFIED bit-exact vs erfa.pfw06/pmat06). The
/// inverse (mean-of-date → GCRS) is the transpose. Reuses obliquity_iau2006 for
/// ε_A, so the orphaned obliquity generator gains its consumer here.
template<typename T>
math::Matrix3<T> precession_iau2006(const math::TrackedValue<T>& t_jcen,
                                    int n_terms = 6) {
    using TV = math::TrackedValue<T>;
    static const char* kGam[6] = {"-0.052928", "10.556378", "0.4932044",
                                  "-0.00031238", "-0.000002788", "0.0000000260"};
    static const char* kPhi[6] = {"84381.412819", "-46.811016", "0.0511268",
                                  "0.00053289", "-0.000000440", "-0.0000000176"};
    static const char* kPsi[6] = {"-0.041775", "5038.481484", "1.5584175",
                                  "-0.00018522", "-0.000026452", "-0.0000000148"};
    static const double kGamMag[6] = {0.052928, 10.556378, 0.4932044,
                                      0.00031238, 0.000002788, 0.0000000260};
    static const double kPhiMag[6] = {84381.412819, 46.811016, 0.0511268,
                                      0.00053289, 0.000000440, 0.0000000176};
    static const double kPsiMag[6] = {0.041775, 5038.481484, 1.5584175,
                                      0.00018522, 0.000026452, 0.0000000148};
    static const math::TrackedPolynomial<T> kGamP = detail::fw_poly<T>(kGam, kGamMag);
    static const math::TrackedPolynomial<T> kPhiP = detail::fw_poly<T>(kPhi, kPhiMag);
    static const math::TrackedPolynomial<T> kPsiP = detail::fw_poly<T>(kPsi, kPsiMag);

    const math::ErrorChannel acc = math::ErrorChannel::accuracy;
    TV a2r = math::pi<T>() / math::exact<T>(648000);   // arcsec → rad
    TV gam = kGamP.eval(t_jcen, n_terms, acc, 1) * a2r;
    TV phi = kPhiP.eval(t_jcen, n_terms, acc, 1) * a2r;
    TV psi = kPsiP.eval(t_jcen, n_terms, acc, 1) * a2r;
    TV eps = obliquity_iau2006<T>(t_jcen, n_terms);
    return rot_x(-eps) * rot_z(-psi) * rot_x(phi) * rot_z(gam);
}

/// Obliquity rotation ecliptic-of-date → equatorial-of-date, Rx(−ε)
/// (frame_chain.md §3). ε is the mean obliquity ε_A (obliquity_iau2006).
template<typename T>
math::Matrix3<T> obliquity_rotation(const math::TrackedValue<T>& eps) {
    return rot_x(-eps);
}

/// Composite ecliptic-of-date → GCRS, M(t) = P(t)ᵀ · Rx(−ε_A) (frame_chain.md §5;
/// = erfa.ecm06ᵀ). Applied to an ecliptic-of-date Cartesian position it yields the
/// GCRS (mean, NO nutation — matching the geometric Meeus ephemeris) position.
template<typename T>
math::Matrix3<T> ecliptic_to_gcrs(const math::TrackedValue<T>& t_jcen,
                                  int n_terms = 6) {
    math::TrackedValue<T> eps = obliquity_iau2006<T>(t_jcen, n_terms);
    return precession_iau2006<T>(t_jcen, n_terms).transpose()
         * obliquity_rotation(eps);
}

} // namespace astronomy
