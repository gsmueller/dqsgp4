#pragma once

/**
 * @file earth_rotation.h
 * @brief The IAU-2006-era Earth rotation angles (ERA, GMST06, GAST06), polar
 *        motion W(t), and the full GCRS→ITRS chain (replan R4b).
 *
 * The Earth-fixed output chain — the honest consumer nutation was waiting for:
 *
 *   r_ITRS = W(t) · R3(GAST) · N(t) · P(t) · r_GCRS
 *
 * with P the FRAME2-gated IAU2006 precession, N the NUT1 nutation matrix, GAST
 * the apparent sidereal angle, and W polar motion from caller-supplied IERS
 * x_p/y_p (measured DATA, not a model — the erfa pom00 signature pattern; the
 * in-repo CSSI EOP file feeds the gate). Distinct from the frozen Aoki-82 GMST
 * (sidereal_time.h), which serves the SGP4 TEME convention and is untouched.
 *
 * Composition order follows erfa exactly (c2tcio/pom00 source sequences); the
 * NUT1 gate verifies every piece AND the full chain element-wise against
 * erfa.gmst06/gst06a/pom00/c2t06a at the chain's stated mas grade.
 * Theory: design/derivations/nutation_itrs.md §3–§5.
 */

#include "epoch.h"
#include "frames.h"
#include "nutation.h"
#include "obliquity.h"
#include "../math/angles.h"
#include "../math/tracked_polynomial.h"
#include "../math/tracked_value.h"

#include <array>
#include <vector>

namespace astronomy {

/// Earth rotation angle ERA(UT1) [rad] — the IAU 2000 Res. B1.8 DEFINING
/// relation: ERA = 2π(0.7790572732640 + 1.00273781191135448·(JD_UT1 − J2000)).
/// Both constants are exact-by-convention (`defined`, CR1B-allowlisted).
template<typename T>
math::TrackedValue<T> earth_rotation_angle(const math::TrackedValue<T>& jd_ut1) {
    using TV = math::TrackedValue<T>;
    TV tu = jd_ut1 - jd_epoch_j2000<T>();
    TV turns = TV::defined("0.7790572732640") + TV::defined("1.00273781191135448") * tu;
    return math::wrap_two_pi(math::two_pi<T>() * turns);
}

/// GMST (IAU 2006 era) [rad] = ERA(UT1) + the gstime00 polynomial in TT
/// centuries (in-repo AstroLib verbatim; measured ≤ 2.4 μas vs erfa.gmst06 over
/// 1995–2024 — the transcription-arbitration bound, deposited in accuracy).
template<typename T>
math::TrackedValue<T> gmst_iau2006(const math::TrackedValue<T>& jd_ut1,
                                   const math::TrackedValue<T>& t_tt) {
    using TV = math::TrackedValue<T>;
    using std::abs;
    static const math::TrackedPolynomial<T> kPoly(
        std::vector<TV>{TV::model_coefficient("0.014506"),
                        TV::model_coefficient("4612.156534"),
                        TV::model_coefficient("1.3915817"),
                        TV::model_coefficient("-0.00000044"),
                        TV::model_coefficient("0.000029956"),
                        TV::model_coefficient("0.0000000368")},
        std::vector<T>{T(0.014506), T(4612.156534), T(1.3915817),
                       T(0.00000044), T(0.000029956), T(0.0000000368)});
    TV a2r = math::pi<T>() / math::exact<T>(648000);
    TV poly = kPoly.eval(t_tt, 6, math::ErrorChannel::accuracy, 1) * a2r;
    TV gmst = earth_rotation_angle<T>(jd_ut1) + poly;
    // Transcription arbitration vs erfa.gmst06 (measured 2.4 μas; note §3).
    T arb = abs((TV::model_coefficient("0.0000025") * a2r).value);
    return math::wrap_two_pi(math::add_bound(gmst, arb, math::ErrorChannel::accuracy));
}

/// GAST (IAU 2006/2000A era) [rad] = GMST06 + the equation of the equinoxes.
/// Convenience form computing everything from scratch at `n_terms` fidelity;
/// the chain composes the pieces itself to share the series evaluation.
template<typename T>
math::TrackedValue<T> gast_iau2006(const math::TrackedValue<T>& jd_ut1,
                                   const math::TrackedValue<T>& t_tt,
                                   int n_terms = kNutationTermCount) {
    NutationAngles<T> ang = nutation_angles<T>(t_tt, n_terms);
    math::TrackedValue<T> eps = obliquity_iau2006<T>(t_tt, 6);
    std::array<math::TrackedValue<T>, 5> args = delaunay_arguments<T>(t_tt);
    math::TrackedValue<T> ee = equation_of_equinoxes<T>(ang, eps, args[4]);
    return math::wrap_two_pi(gmst_iau2006<T>(jd_ut1, t_tt) + ee);
}

/// Polar motion W(t) from caller-supplied IERS pole coordinates x_p, y_p [rad]
/// and t [TT centuries] for the TIO locator s′ = −47 μas·t (IERS 2010 §5.3).
/// Composition follows the erfa pom00 source sequence exactly:
/// W = R1(−y_p) · R2(−x_p) · R3(s′) in SOFA rotations (our rot_*).
template<typename T>
math::Matrix3<T> polar_motion(const math::TrackedValue<T>& xp,
                              const math::TrackedValue<T>& yp,
                              const math::TrackedValue<T>& t_tt) {
    using TV = math::TrackedValue<T>;
    TV a2r = math::pi<T>() / math::exact<T>(648000);
    TV sp = TV::model_coefficient("-0.000047") * t_tt * a2r;
    return rot_x(-yp) * rot_y(-xp) * rot_z(sp);
}

/// The full GCRS → ITRS rotation at TT centuries `t_tt`, UT1 Julian date
/// `jd_ut1`, and IERS pole coordinates x_p, y_p [rad]:
///   M = W(t) · R3(GAST) · N(t) · P(t)            (equinox-based; ≡ erfa c2t06a
/// to the chain's mas grade — NUT1 measures it element-wise). `n_terms` dials
/// the nutation series; every piece's tracked budget composes through.
template<typename T>
math::Matrix3<T> gcrs_to_itrs(const math::TrackedValue<T>& t_tt,
                              const math::TrackedValue<T>& jd_ut1,
                              const math::TrackedValue<T>& xp,
                              const math::TrackedValue<T>& yp,
                              int n_terms = kNutationTermCount) {
    NutationAngles<T> ang = nutation_angles<T>(t_tt, n_terms);
    math::TrackedValue<T> eps = obliquity_iau2006<T>(t_tt, 6);
    std::array<math::TrackedValue<T>, 5> args = delaunay_arguments<T>(t_tt);
    math::TrackedValue<T> ee = equation_of_equinoxes<T>(ang, eps, args[4]);
    math::TrackedValue<T> gast =
        math::wrap_two_pi(gmst_iau2006<T>(jd_ut1, t_tt) + ee);

    math::Matrix3<T> P = precession_iau2006<T>(t_tt);
    math::Matrix3<T> N = nutation_matrix<T>(eps, ang);
    math::Matrix3<T> W = polar_motion<T>(xp, yp, t_tt);
    return W * rot_z(gast) * (N * P);
}

} // namespace astronomy
