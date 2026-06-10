#pragma once

/// @file body_position_gcrs.h
/// Compose an ephemeris instance's ecliptic-of-date position (λ, β, Δ) into a GCRS
/// Cartesian position via the L2/L3 frame chain (design/derivations/frame_chain.md
/// §2 + §5). This is the CONSUMER the L2-deferred precession / obliquity rotations
/// were held back for (anti-dead-code): they land in astronomy/frames.h, this wires
/// them to a real ephemeris.
///
///   r_ecl  = Δ · (cos β cos λ, cos β sin λ, sin β)        (ecliptic Cartesian, §2)
///   r_gcrs = ecliptic_to_gcrs(t) · r_ecl                  (= Pᵀ·Rx(−ε_A), §5)
///
/// Generic over the body — the caller supplies the EclipticState (sun_meeus /
/// moon_meeus), so one function serves every analytical ephemeris instance.
///
/// Fidelity (no perceived fidelity): gated END-TO-END against the in-repo JPL
/// DE430 ephemeris (the independent NUMERICAL reference, not the analytical theory
/// the model realizes) in test_frame_chain — the angular residual is the body's
/// Meeus truncation, the orthonormal rotation contributing ≤1e-12.

#include "ecliptic_state.h"

#include "../astronomy/epoch.h"
#include "../astronomy/frames.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

namespace ephemeris {

/// GCRS Cartesian position of a body from its ecliptic-of-date state at `epoch`.
/// `n_terms` is the FRAME (precession/obliquity) series depth; the ephemeris depth
/// is fixed when `ecl` is computed. Units follow `ecl.radius` (AU Sun / km Moon).
template<typename T>
math::Vector3<T> body_position_gcrs(const EclipticState<T>& ecl,
                                    const astronomy::Epoch<T>& epoch,
                                    int n_terms = 6) {
    using TV = math::TrackedValue<T>;
    TV cos_b = cos(ecl.lat);
    math::Vector3<T> r_ecl(ecl.radius * cos_b * cos(ecl.lon),
                           ecl.radius * cos_b * sin(ecl.lon),
                           ecl.radius * sin(ecl.lat));
    TV t = astronomy::centuries_since_j2000<T>(epoch);
    return astronomy::ecliptic_to_gcrs<T>(t, n_terms) * r_ecl;
}

} // namespace ephemeris
