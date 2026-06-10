#pragma once

/// @file attitude.h
/// Coupled-axis rigid-body attitude dynamics (issue H1).
///
/// The diagonal dynamics/inertia.h computes ω̇_i = τ_i / I_ii — it OMITS the
/// gyroscopic coupling, so the body axes do not exchange angular momentum and a
/// torque-free asymmetric body does not precess. This file adds the full Euler
/// rigid-body equation, valid for a general (off-diagonal) inertia tensor:
///
///   I·ω̇ = τ − ω × (I·ω),   hence   ω̇ = I⁻¹ ( τ − ω × (I·ω) ),
///
/// where I is a math::Matrix3<T> (products of inertia allowed). It pairs this
/// with the attitude kinematics q̇ = ½ q ⊗ (0, ω) and an RK4 step that
/// integrates (q, ω) and retracts q to unit norm. Unlike the diagonal
/// Inertia::acceleration_from_wrench, a nonzero torque on a full inertia is
/// well defined (no zero-moment catastrophe), so the point-mass torque
/// restriction is lifted for bodies that carry a real inertia tensor.

#include "../math/matrix3.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

namespace dynamics {

/// Rigid-body rotational state: orientation quaternion + body angular velocity.
template<typename T>
struct AttitudeState {
    math::Quaternion<T> q;       ///< Body orientation (unit quaternion).
    math::Vector3<T> omega;      ///< Body angular velocity ω [rad/s].
};

/// Body angular acceleration via the full Euler equation (gyroscopic coupling).
/// @param I      Inertia tensor (body frame).
/// @param I_inv  Its inverse (precomputed once).
/// @param omega  Body angular velocity.
/// @param torque Applied body torque (zero for torque-free motion).
template<typename T>
math::Vector3<T> euler_angular_acceleration(
    const math::Matrix3<T>& I, const math::Matrix3<T>& I_inv,
    const math::Vector3<T>& omega, const math::Vector3<T>& torque) {
    math::Vector3<T> Iw = I * omega;
    math::Vector3<T> gyro = omega.cross(Iw);   // ω × (I·ω)
    math::Vector3<T> net = torque - gyro;
    return I_inv * net;
}

/// Attitude kinematics: q̇ = ½ q ⊗ (0, ω).
template<typename T>
math::Quaternion<T> quaternion_rate(const math::Quaternion<T>& q,
                                    const math::Vector3<T>& omega) {
    return q * math::Quaternion<T>::pure(omega) * math::ratio<T>(1, 2);
}

/// Rotational kinetic energy ½ ω·(I·ω). Conserved for torque-free motion.
template<typename T>
math::TrackedValue<T> rotational_energy(const math::Matrix3<T>& I,
                                        const math::Vector3<T>& omega) {
    return math::ratio<T>(1, 2) * omega.dot(I * omega);
}

/// Body-frame angular momentum L = I·ω. Its magnitude equals the conserved
/// inertial |L| for torque-free motion.
template<typename T>
math::Vector3<T> angular_momentum(const math::Matrix3<T>& I,
                                  const math::Vector3<T>& omega) {
    return I * omega;
}

/// One RK4 step of the coupled (q, ω) attitude dynamics under a constant body
/// torque. The quaternion is renormalized after the step (Lie-group retraction
/// back to the unit sphere — RK4 does not preserve the norm exactly).
template<typename T>
AttitudeState<T> attitude_rk4_step(
    const AttitudeState<T>& s, const math::Matrix3<T>& I,
    const math::Matrix3<T>& I_inv, const math::Vector3<T>& torque,
    const math::TrackedValue<T>& dt) {
    using math::exact;
    using math::ratio;

    auto deriv = [&](const AttitudeState<T>& st) -> AttitudeState<T> {
        return AttitudeState<T>{
            quaternion_rate(st.q, st.omega),
            euler_angular_acceleration(I, I_inv, st.omega, torque)};
    };
    auto add = [](const AttitudeState<T>& a, const AttitudeState<T>& b) -> AttitudeState<T> {
        return AttitudeState<T>{a.q + b.q, a.omega + b.omega};
    };
    auto scale = [](const math::TrackedValue<T>& c, const AttitudeState<T>& a) -> AttitudeState<T> {
        return AttitudeState<T>{c * a.q, c * a.omega};
    };

    math::TrackedValue<T> half = dt * ratio<T>(1, 2);
    AttitudeState<T> k1 = deriv(s);
    AttitudeState<T> k2 = deriv(add(s, scale(half, k1)));
    AttitudeState<T> k3 = deriv(add(s, scale(half, k2)));
    AttitudeState<T> k4 = deriv(add(s, scale(dt, k3)));

    // s + (dt/6)(k1 + 2 k2 + 2 k3 + k4)
    AttitudeState<T> sum =
        add(add(k1, scale(exact<T>(2), k2)), add(scale(exact<T>(2), k3), k4));
    AttitudeState<T> next = add(s, scale(dt / exact<T>(6), sum));
    next.q = next.q.normalized();
    return next;
}

} // namespace dynamics
