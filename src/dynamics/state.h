#pragma once

/// @file state.h
/// State: the propagator's integrated quantity at a single instant.
///
/// A `State<T>` bundles
///
///   pose   — `Pose<T>`, an SE(3) element (position + orientation)
///   twist  — `Twist<T>`, body-frame angular + linear velocity
///   time   — `TrackedValue<T>`, elapsed time since epoch
///
/// Together these are the variables the integrator advances. The
/// derivative of each is described in `dynamics/derivative.h`.
///
/// All sub-objects are error-tracked composites. Per-component
/// `errors.measurement / precision / accuracy` accessors propagate
/// uniformly from the leaf `TrackedValue<T>` storage (REQ-EF-12,
/// REQ-EF-13).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-9,
///   AUD-CC-12, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"
#include "pose.h"
#include "twist.h"

namespace dynamics {

/// Propagator state at a single time.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct State {
    Pose<T> pose;                        ///< SE(3) element.
    Twist<T> twist;                      ///< Body-frame velocity.
    math::TrackedValue<T> time;          ///< Elapsed time since epoch.

    // --- Constructors ---

    /// Identity pose, zero twist, t = 0.
    State() : pose(), twist(), time() {}

    /// From explicit pose, twist, and elapsed time.
    State(const Pose<T>& p, const Twist<T>& w, const math::TrackedValue<T>& t)
        : pose(p), twist(w), time(t) {}

    // --- Named constructors ---

    /// Identity state at t = 0: origin, identity rotation, zero twist.
    static State identity_at_epoch() {
        return State(Pose<T>::identity(), Twist<T>::zero(), math::exact<T>(0));
    }

    /// Build from explicit kinematic ingredients: a unit rotation
    /// quaternion, a translation 3-vector, a body angular velocity, a
    /// body linear velocity, and an epoch time.
    static State from_kinematics(const math::Quaternion<T>& q_r,
                                 const math::Vector3<T>& position,
                                 const math::Vector3<T>& angular_velocity,
                                 const math::Vector3<T>& linear_velocity,
                                 const math::TrackedValue<T>& t) {
        return State(
            Pose<T>::from_rotation_translation(q_r, position),
            Twist<T>(angular_velocity, linear_velocity),
            t
        );
    }

    // --- Convenience accessors ---

    /// Position 3-vector recovered from the pose.
    math::Vector3<T> position() const { return pose.translation(); }

    /// Orientation as a unit rotation quaternion.
    math::Quaternion<T> orientation() const { return pose.rotation(); }

    /// Body angular velocity 3-vector.
    math::Vector3<T> angular_velocity() const { return twist.angular; }

    /// Body linear velocity 3-vector.
    math::Vector3<T> linear_velocity() const { return twist.linear; }

    // --- Error book-keeping ---

    /// Deposit the integrators' per-step LTE envelope PER SLOT FAMILY with
    /// DIMENSIONAL magnitudes (runge_kutta_lie_group.md §6.2): `pos` [m] and
    /// `rot` [rad] are the pose-error magnitudes P and Θ; `vel` [m/s] and
    /// `angvel` [rad/s] the twist-error magnitudes V and Ω (the rotational
    /// pair is exactly 0 for torque-free motion). The unit-dual-quaternion
    /// slot map:
    ///   real slots    += Θ/2                       (small-angle |δq| ≤ Θ/2)
    ///   dual slots    += P/2 + ‖q_d‖₁·Θ/2          (q_d = ½ t⊗q_r product rule)
    ///   twist angular += Ω,  twist linear += V     (their own parts)
    /// This replaces the former uniform 14-slot stamp, whose translational-
    /// scale, mis-united magnitude on the O(1) real slots was amplified by
    /// position scale at extraction (r = 2·q_d⊗q_r*) — measured 4.3e13 m
    /// after one 30 s RK4 step, inf over an hour.
    void add_step_lte(const T& pos, const T& rot, const T& vel, const T& angvel) {
        using std::abs;
        const T half_rot = rot / T(2);
        const T qd_l1 = abs(pose.M.dual.w.value) + abs(pose.M.dual.x.value)
                      + abs(pose.M.dual.y.value) + abs(pose.M.dual.z.value);
        const T dual_bound = pos / T(2) + qd_l1 * half_rot;
        pose.M.real.w.errors.accuracy += half_rot;
        pose.M.real.x.errors.accuracy += half_rot;
        pose.M.real.y.errors.accuracy += half_rot;
        pose.M.real.z.errors.accuracy += half_rot;
        pose.M.dual.w.errors.accuracy += dual_bound;
        pose.M.dual.x.errors.accuracy += dual_bound;
        pose.M.dual.y.errors.accuracy += dual_bound;
        pose.M.dual.z.errors.accuracy += dual_bound;
        twist.angular.x.errors.accuracy += angvel;
        twist.angular.y.errors.accuracy += angvel;
        twist.angular.z.errors.accuracy += angvel;
        twist.linear.x.errors.accuracy += vel;
        twist.linear.y.errors.accuracy += vel;
        twist.linear.z.errors.accuracy += vel;
    }
};

} // namespace dynamics
