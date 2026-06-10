#pragma once

/// @file gravity_tesseral.h
/// Tesseral / sectoral (longitude-dependent) geopotential acceleration. Issue
/// D2 — complements the longitude-independent zonals in gravity_zonal.h.
///
/// Uses the Cunningham (1970) / Montenbruck & Gill ("Satellite Orbits" §3.2.4)
/// real-valued Cartesian recursion for the harmonic functions
///
///   V_nm = (R/r)^{n+1} P_nm(sin φ) cos(mλ),
///   W_nm = (R/r)^{n+1} P_nm(sin φ) sin(mλ),
///
/// evaluated directly in Earth-fixed Cartesian coordinates — singularity-free
/// (no division by cos φ at the poles) and extensible to any degree/order. The
/// acceleration follows Montenbruck & Gill eqn (3.33) using the UNnormalized
/// C_nm, S_nm from constants/tesseral_harmonics.h.
///
/// Frame. The tesseral field is fixed to the rotating Earth, so the satellite's
/// ECI (TEME) position is rotated into ECEF by the L2 sidereal rotation
/// R = Rz(GMST) (astronomy/frames.h), the acceleration is evaluated there, and
/// rotated back to ECI by Rᵀ. GMST is supplied by the caller (the force lambda
/// has the state's elapsed time but not the absolute epoch). Same body-frame
/// return convention as gravity_central / gravity_zonal.

#include "../astronomy/frames.h"   // L2: sidereal_rotation (TEME→ECEF, Rz(GMST))
#include "../constants/constants_provider.h"
#include "../constants/tesseral_harmonics.h"
#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/matrix3.h"
#include "../math/quaternion.h"
#include "../math/spherical_harmonics.h"   // shared Cunningham V/W recursion
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <vector>

namespace forces {

/// Tesseral/sectoral acceleration in the EARTH-FIXED frame at r_ecef [same units
/// as Re]. Sums the m ≥ 1 contributions of the harmonics up to degree
/// min(max_n, table degree). Returns acceleration in the same length/time units
/// implied by GM and Re.
template<typename T>
math::Vector3<T> gravity_tesseral_accel_ecef(
    const math::Vector3<T>& r_ecef,
    const math::TrackedValue<T>& mu,
    const math::TrackedValue<T>& Re,
    const constants::TesseralHarmonics<T>& th,
    int max_n) {
    using math::TrackedValue;
    using math::exact;
    using TV = TrackedValue<T>;

    const int N = (max_n < th.max_degree()) ? max_n : th.max_degree();
    const int top = N + 1;  // V/W needed up to degree+order N+1

    // V[n][m], W[n][m] via the shared Cunningham recursion (math/spherical_harmonics.h).
    // Bit-identical to the prior inline recursion (same operations, same order).
    math::VWTable<T> vw = math::cunningham_vw(r_ecef.x, r_ecef.y, r_ecef.z, Re, top);
    std::vector<std::vector<TV>>& V = vw.V;
    std::vector<std::vector<TV>>& W = vw.W;

    // Acceleration (Montenbruck-Gill 3.33), tesseral terms (m ≥ 1) only.
    TV scale = mu / (Re * Re);
    TV ax = exact<T>(0), ay = exact<T>(0), az = exact<T>(0);
    for (const auto& t : th.terms) {
        const int n = t.n, m = t.m;
        if (n > N || m < 1) continue;
        const TV& C = t.C;
        const TV& S = t.S;
        TV f = exact<T>((n - m + 2) * (n - m + 1));  // (n−m+2)!/(n−m)!
        ax = ax + math::ratio<T>(1, 2) * (
            (-C * V[n + 1][m + 1] - S * W[n + 1][m + 1])
            + f * (C * V[n + 1][m - 1] + S * W[n + 1][m - 1]));
        ay = ay + math::ratio<T>(1, 2) * (
            (-C * W[n + 1][m + 1] + S * V[n + 1][m + 1])
            + f * (-C * W[n + 1][m - 1] + S * V[n + 1][m - 1]));
        az = az + exact<T>(n - m + 1) * (-C * V[n + 1][m] - S * W[n + 1][m]);
    }
    return math::Vector3<T>(scale * ax, scale * ay, scale * az);
}

/// Body-frame tesseral wrench for the propagator force list. Rotates the ECI
/// position to ECEF by the sidereal rotation R = Rz(GMST), evaluates the tesseral
/// acceleration, rotates it back to ECI by Rᵀ, then inverse-rotates to the body
/// frame.
template<typename T>
dynamics::Wrench<T> gravity_tesseral(
    const dynamics::State<T>& state,
    const constants::ConstantsProvider<T>& K,
    const constants::TesseralHarmonics<T>& th,
    int max_n,
    const math::TrackedValue<T>& gmst) {
    math::Vector3<T> r_eci = state.position();

    // ECI(TEME) → ECEF and back via the L2 sidereal rotation R = Rz(GMST)
    // (astronomy/frames.h). R · r_eci reproduces the prior hand-coded
    // componentwise Rz arithmetic BIT-FOR-BIT — value AND error channels: the
    // +0·z matrix terms are exact no-ops, 1·z is IEEE-exact, and the transpose's
    // a + (−b) equals the prior a − b exactly. So this is value-preserving
    // (OR1/33/33 unaffected); verified by test_reference_frames (FRAME1) and
    // test_tesseral (D2). cos/sin computed once inside rot_z, as before.
    math::Matrix3<T> R = astronomy::sidereal_rotation(gmst);
    math::Vector3<T> r_ecef = R * r_eci;

    math::Vector3<T> a_ecef =
        gravity_tesseral_accel_ecef<T>(r_ecef, K.earth.GM, K.earth.a, th, max_n);

    math::Vector3<T> a_eci = R.transpose() * a_ecef;

    math::Vector3<T> a_body = state.pose.rotation().conjugate().rotate(a_eci);
    return dynamics::Wrench<T>(math::Vector3<T>(), a_body);
}

} // namespace forces
