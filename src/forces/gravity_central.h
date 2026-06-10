#pragma once

/// @file gravity_central.h
/// Central (point-mass) gravity wrench lambda.
///
/// Computes the body-frame wrench Ŵ = (τ, F) due to a point-mass
/// gravitational attractor at the world-frame origin:
///
///   F_world = −(GM / |r|³) · r
///
/// where r = state.pose.translation() is the body's position in the world
/// frame. The force is then transformed to the body frame via the inverse
/// rotation:
///
///   F_body = q_r* · F_world · q_r.
///
/// No torque contribution: a point-mass attractor exerts no moment.
///
/// Model truncation. Central gravity captures only the spherically
/// symmetric part of Earth's field. The longitude-independent zonals are
/// added by `gravity_zonal` (Jₙ, issue D1) and the longitude-dependent
/// tesseral/sectoral harmonics by `gravity_tesseral` (C_nm/S_nm, issue D2);
/// each adds its model-truncation contribution to `errors.accuracy` per
/// REQ-EF-7.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../constants/constants_provider.h"
#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <functional>

namespace forces {

/// Compute the central-gravity body-frame wrench at a given state.
///
/// @param state  Propagator state (uses state.pose only).
/// @param K      Constants provider (uses K.earth.GM).
/// @return       Body-frame Wrench<T> with zero torque and gravity force.
///
/// @tparam T  Underlying numeric type.
template<typename T>
dynamics::Wrench<T> gravity_central(const dynamics::State<T>& state,
                                    const constants::ConstantsProvider<T>& K) {
    math::Vector3<T> r_world = state.position();
    math::TrackedValue<T> r_sq =
        r_world.x * r_world.x + r_world.y * r_world.y + r_world.z * r_world.z;
    math::TrackedValue<T> r_mag = sqrt(r_sq);
    math::TrackedValue<T> r_cubed = r_mag * r_sq;

    math::TrackedValue<T> neg_gm_over_r3 = -(K.earth.GM / r_cubed);
    math::Vector3<T> F_world(
        neg_gm_over_r3 * r_world.x,
        neg_gm_over_r3 * r_world.y,
        neg_gm_over_r3 * r_world.z
    );

    // Inverse rotation: q_r* · F_world · q_r. `rotate` applies q v q*,
    // so applying the conjugate quaternion gives the inverse rotation.
    math::Vector3<T> F_body =
        state.pose.rotation().conjugate().rotate(F_world);

    // REQ-EF-6 / AUD-EF-6: the spherically-symmetric monopole −GM·r/|r|³ is
    // EXACT (Newton's shell theorem), so its model-truncation contribution to
    // errors.accuracy is identically zero. Record it explicitly so the zero is
    // auditable (every force states its residual) rather than a silent omission;
    // the neglected Jₙ harmonics are accounted for by their own force lambdas
    // (gravity_zonal, …), which add their nonzero residuals.
    const T model_truncation = T(0);
    F_body.x.errors.accuracy = F_body.x.errors.accuracy + model_truncation;
    F_body.y.errors.accuracy = F_body.y.errors.accuracy + model_truncation;
    F_body.z.errors.accuracy = F_body.z.errors.accuracy + model_truncation;
    return dynamics::Wrench<T>(math::Vector3<T>(), F_body);
}

/// Build a wrench-fn lambda suitable for the propagator's force list.
/// Captures the constants provider by const reference; the caller must
/// ensure the provider outlives the propagator.
template<typename T>
std::function<dynamics::Wrench<T>(const dynamics::State<T>&)>
make_gravity_central(const constants::ConstantsProvider<T>& K) {
    return [&K](const dynamics::State<T>& state) -> dynamics::Wrench<T> {
        return gravity_central(state, K);
    };
}

} // namespace forces
