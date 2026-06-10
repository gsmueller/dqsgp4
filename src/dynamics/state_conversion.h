#pragma once

/// @file state_conversion.h
/// Converters between the two propagators' state types, giving the SGP4 and
/// DQSGP4 APIs a shared vocabulary (item F1):
///
///   to_state_vector(State)         dynamics::State<T> (metres, inertial) ->
///                                  sgp4::StateVector<T> (km, TEME)
///   to_state(StateVector[, time])  the inverse (km -> m, identity attitude)
///
/// The DQ twist linear velocity is body-frame; it is rotated into the world
/// frame by the pose orientation before scaling (a no-op for the identity-
/// attitude orbital case, correct for the general attitude case). Units
/// convert exactly (km <-> m is x1000). Including this header opts into both
/// the dynamics and sgp4 state types; it does NOT pull in a propagator.

#include "state.h"
#include "../sgp4/state_vector.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

namespace dynamics {

/// dynamics::State<T> (metres, inertial) -> sgp4::StateVector<T> (km, TEME).
/// The body-frame linear velocity is rotated into the world frame by the pose
/// orientation, then position and velocity are scaled metres -> km.
template<typename T>
sgp4::StateVector<T> to_state_vector(const State<T>& s) {
    const math::TrackedValue<T> m_to_km = math::ratio<T>(1, 1000);
    math::Vector3<T> v_world = s.orientation().rotate(s.linear_velocity());
    sgp4::StateVector<T> sv;
    sv.position_km = s.position() * m_to_km;
    sv.velocity_km_s = v_world * m_to_km;
    return sv;
}

/// sgp4::StateVector<T> (km, TEME) -> dynamics::State<T> (metres, inertial).
/// Scales km -> m and packs into a State with identity attitude and zero body
/// rate (a StateVector carries no attitude), at the given time (default 0).
template<typename T>
State<T> to_state(const sgp4::StateVector<T>& sv,
                  const math::TrackedValue<T>& time = math::exact<T>(0)) {
    const math::TrackedValue<T> km_to_m = math::exact<T>(1000);
    math::Vector3<T> position_m = sv.position_km * km_to_m;
    math::Vector3<T> velocity_m = sv.velocity_km_s * km_to_m;
    math::Vector3<T> zero_rate;
    return State<T>::from_kinematics(math::Quaternion<T>::identity(),
                                     position_m, zero_rate, velocity_m, time);
}

}  // namespace dynamics
