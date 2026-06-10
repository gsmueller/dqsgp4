#pragma once

/// @file symplectic_leapfrog.h
/// Symplectic leapfrog (kick–drift–kick velocity Verlet) integrator for
/// time-independent conservative problems (REQ-IN-9).
///
/// For the split system  dtwist/dt = a(state),  dpose/dt = lie(twist):
///
///   twist_half = twist + (dt/2) a(y0)            (half kick)
///   pose_new   = lie_advance_pose(pose, dt, twist_half)   (drift)
///   twist_new  = twist_half + (dt/2) a(pose_new)  (half kick)
///
/// The scheme is second-order and SYMPLECTIC: for a time-independent
/// conservative force it has NO secular energy drift (the energy error stays
/// bounded for arbitrarily many steps), and for a central force it conserves
/// angular momentum to round-off — unlike a non-symplectic Runge–Kutta method,
/// whose energy drifts secularly. Use it for long conservative propagations;
/// use `runge_kutta_4` / `rkf78` when per-step accuracy matters more than
/// long-horizon invariant preservation.
///
/// Stage / complexity (REQ-IN-12): 2 acceleration-callback evaluations, 1
/// pose advance (one drift) + 1 retraction per step.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-9,
///   AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../dynamics/state.h"
#include "../dynamics/twist.h"
#include "../math/tracked_value.h"
#include "runge_kutta.h"

namespace integrators {

/// One symplectic leapfrog (velocity-Verlet KDK) step.
template<typename T>
dynamics::State<T> symplectic_leapfrog(const dynamics::State<T>& y0,
                                       const math::TrackedValue<T>& dt,
                                       const AccelFn<T>& accel_fn) {
    math::TrackedValue<T> half_dt = dt * math::ratio<T>(1, 2);

    // Half kick at the current state.
    dynamics::Twist<T> a0 = accel_fn(y0);
    dynamics::Twist<T> twist_half = y0.twist + half_dt * a0;

    // Drift the pose over the full step using the half-step velocity.
    dynamics::Pose<T> pose_new = lie_advance_pose(y0.pose, dt, twist_half);

    // Half kick at the drifted state (force depends on the new pose).
    dynamics::State<T> drifted(pose_new, twist_half, y0.time + dt);
    dynamics::Twist<T> a1 = accel_fn(drifted);
    dynamics::Twist<T> twist_new = twist_half + half_dt * a1;

    dynamics::State<T> next(pose_new, twist_new, y0.time + dt);

    // R07 / REQ-EF-7 / REQ-IN-8: the per-step LTE envelope
    // (runge_kutta_lie_group.md §6.2/§6.3, identical to euler/rk4):
    // position/rotation (h²/2)·A, velocity/angular-rate h·A, with
    // max(A0, A1) as the derivative-cascade proxy, split by part, per slot
    // family — the rotational pair is exactly 0 torque-free. (The GLOBAL
    // energy error is bounded by symplecticity, but the per-step truncation
    // is still recorded.)
    {
        using std::max;
        T a_lin = max(a0.l1_norm_linear(), a1.l1_norm_linear());
        T a_ang = max(a0.l1_norm_angular(), a1.l1_norm_angular());
        T h2_half = dt.value * dt.value / T(2);
        next.add_step_lte(h2_half * a_lin, h2_half * a_ang,
                          dt.value * a_lin, dt.value * a_ang);
    }
    return next;
}

} // namespace integrators
