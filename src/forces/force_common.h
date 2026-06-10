#pragma once

/// @file force_common.h
/// Shared epilogue for body-frame force wrenches. Every gravity/drag force lambda
/// ends with the SAME step: rotate a world-frame (ECI) linear acceleration into
/// the body frame by the pose's inverse rotation and package it as a force-only
/// Wrench (zero torque). Lifting it here removes the 4-copy duplication across
/// gravity_central / gravity_zonal / gravity_tesseral / drag
/// (design/derivations/geopotential.md §6).

#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/quaternion.h"
#include "../math/vector3.h"

namespace forces {

/// Body-frame force-only wrench from a world-frame (ECI) linear acceleration.
template<typename T>
dynamics::Wrench<T> force_wrench_from_world_accel(
    const dynamics::State<T>& state, const math::Vector3<T>& a_world) {
    math::Vector3<T> a_body = state.pose.rotation().conjugate().rotate(a_world);
    return dynamics::Wrench<T>(math::Vector3<T>(), a_body);
}

}  // namespace forces
