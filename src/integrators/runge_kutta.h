#pragma once

/// @file runge_kutta.h
/// Lie-group Runge–Kutta integrators for the dual quaternion propagator.
///
/// The state ODE is
///
///   dM̂/dt = (1/2) M̂ · Ω̂_pure          (Lie-group, evolves on SE(3))
///   dΩ̂/dt = a_body                       (linear in body twist space)
///   dt/dt = 1.
///
/// Standard Runge–Kutta cannot be applied directly to the pose component
/// because the unit-DQ manifold is not a vector space. The
/// implementations here use a Munthe-Kaas style update: at each stage,
/// the pose is advanced from y0 by the SE(3) screw exponential of the
/// stage-derived body twist; the twist is advanced linearly in its
/// vector-space tangent. After the final stage, the pose is renormalized
/// (REQ-EF-15 retraction) to project off any sub-ULP drift.
///
/// The body-acceleration callback is supplied by the propagator and
/// captures the constants provider, force lambdas, and inertia. Its
/// signature is intentionally narrow: takes a stage `State<T>`, returns
/// a body-acceleration `Twist<T>`.
///
/// All three explicit members (euler, runge_kutta_4, and rkf78 in
/// runge_kutta_fehlberg.h) are ONE driver — `rk_step` — applied to three
/// `ButcherTableau`s. Each entry point selects its tableau and applies its
/// own truncation-error policy; `rk_step` itself does pure integration.
/// Theory: design/derivations/runge_kutta_lie_group.md.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/twist.h"
#include "../math/dual_quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <algorithm>
#include <array>
#include <functional>

namespace integrators {

/// Body-acceleration callback signature. Takes a state, returns the body
/// acceleration twist (ω_dot, v_dot) for that state.
template<typename T>
using AccelFn =
    std::function<dynamics::Twist<T>(const dynamics::State<T>&)>;

/// Lie-group pose advance: pose ← pose · exp_screw((δt/2) · twist).
///
/// The half-factor in the screw exponential argument is the half-angle
/// convention of `DualQuaternion<T>::exp_screw` (REQ-DQ-18). This routine
/// applies the rigid-body screw motion to `pose` over time `delta_t`
/// using a constant body twist, then projects back to the SE(3) manifold
/// via `Pose::normalized()`.
template<typename T>
dynamics::Pose<T> lie_advance_pose(const dynamics::Pose<T>& pose,
                                   const math::TrackedValue<T>& delta_t,
                                   const dynamics::Twist<T>& twist) {
    math::TrackedValue<T> half_delta = delta_t * math::ratio<T>(1, 2);
    math::Vector3<T> u = half_delta * twist.angular;
    math::Vector3<T> v = half_delta * twist.linear;
    math::DualQuaternion<T> step = math::DualQuaternion<T>::exp_screw(u, v);
    return dynamics::Pose<T>(pose.M * step).normalized();
}

/// Butcher tableau for an explicit S-stage Runge–Kutta method: nodes `c`,
/// 8th/final-order weights `b`, and a strictly-lower-triangular stage matrix
/// `a` (a[i][j] used only for j < i). Coefficients are exact rationals
/// (`ratio<T>`) / exact integers (`exact<T>`); unset entries are exact zero.
/// The consistency identities — row sum Σⱼ a[i][j] = c[i], Σᵢ b[i] = 1, and
/// the quadrature order conditions Σᵢ b[i] c[i]ᵏ = 1/(k+1) — are what
/// `test_butcher_tableau` checks exactly from the rationals (no numerical
/// reference). See design/derivations/runge_kutta_lie_group.md §2.
template<typename T, int S>
struct ButcherTableau {
    std::array<math::TrackedValue<T>, S> c;                   ///< stage nodes
    std::array<math::TrackedValue<T>, S> b;                   ///< combine weights
    std::array<std::array<math::TrackedValue<T>, S>, S> a;    ///< stage matrix

    /// All entries exact zero; the named factories fill the method's rationals.
    ButcherTableau() {
        for (auto& ci : c) ci = math::exact<T>(0);
        for (auto& bi : b) bi = math::exact<T>(0);
        for (auto& row : a) for (auto& aij : row) aij = math::exact<T>(0);
    }
};

/// Output of one `rk_step`: the integrated state and the per-stage body
/// accelerations (the entry points need these for their truncation-error
/// policies; `rk_step` itself deposits no error budget).
template<typename T, int S>
struct RkStepResult {
    dynamics::State<T> state;
    std::array<dynamics::Twist<T>, S> stage_accel;
};

/// Generic Lie-group explicit-RK step (REQ-IN-1/4). For each stage i the twist
/// advances in the flat se(3) tangent, Wᵢ = twist₀ + dt·Σⱼ a[i][j] Aⱼ, and the
/// pose advances from y0 by the SINGLE screw exponential of the RK twist
/// combination, lie_advance_pose(pose₀, dt, Σⱼ a[i][j] Wⱼ); Aᵢ is the callback
/// at that stage state. The output advances the pose by the b-weighted stage
/// twists and the twist by the b-weighted accelerations. This is the one loop
/// euler (S=1), runge_kutta_4 (S=4), and rkf78 (S=13) all run.
/// Derivation: design/derivations/runge_kutta_lie_group.md §3–§4.
template<typename T, int S>
RkStepResult<T, S> rk_step(const ButcherTableau<T, S>& tab,
                           const dynamics::State<T>& y0,
                           const math::TrackedValue<T>& dt,
                           const AccelFn<T>& accel_fn) {
    std::array<dynamics::Twist<T>, S> W;  // stage twists
    std::array<dynamics::Twist<T>, S> A;  // stage accelerations

    for (int i = 0; i < S; ++i) {
        dynamics::Twist<T> twist_inc = dynamics::Twist<T>::zero();   // dt sum a_ij A_j
        dynamics::Twist<T> pose_twist = dynamics::Twist<T>::zero();  // sum a_ij W_j
        for (int j = 0; j < i; ++j) {
            twist_inc = twist_inc + tab.a[i][j] * A[j];
            pose_twist = pose_twist + tab.a[i][j] * W[j];
        }
        W[i] = y0.twist + dt * twist_inc;
        dynamics::State<T> stage(
            lie_advance_pose(y0.pose, dt, pose_twist),
            W[i],
            y0.time + tab.c[i] * dt);
        A[i] = accel_fn(stage);
    }

    dynamics::Twist<T> twist_out = dynamics::Twist<T>::zero();  // sum b_i A_i
    dynamics::Twist<T> pose_out = dynamics::Twist<T>::zero();   // sum b_i W_i
    for (int i = 0; i < S; ++i) {
        twist_out = twist_out + tab.b[i] * A[i];
        pose_out = pose_out + tab.b[i] * W[i];
    }

    dynamics::State<T> next(
        lie_advance_pose(y0.pose, dt, pose_out),
        y0.twist + dt * twist_out,
        y0.time + dt);

    return RkStepResult<T, S>{next, A};
}

/// One-stage tableau: c=(0), b=(1), a=((0)) — forward Euler.
template<typename T>
ButcherTableau<T, 1> euler_tableau() {
    ButcherTableau<T, 1> tab;
    tab.b[0] = math::exact<T>(1);
    return tab;
}

/// Classical RK4 tableau: c=(0,½,½,1), b=(⅙,⅓,⅓,⅙), a₁₀=a₂₁=½, a₃₂=1.
template<typename T>
ButcherTableau<T, 4> rk4_tableau() {
    using math::exact;
    using math::ratio;
    ButcherTableau<T, 4> tab;
    tab.c[1] = ratio<T>(1, 2); tab.c[2] = ratio<T>(1, 2); tab.c[3] = exact<T>(1);
    tab.b[0] = ratio<T>(1, 6); tab.b[1] = ratio<T>(1, 3);
    tab.b[2] = ratio<T>(1, 3); tab.b[3] = ratio<T>(1, 6);
    tab.a[1][0] = ratio<T>(1, 2);
    tab.a[2][1] = ratio<T>(1, 2);
    tab.a[3][2] = exact<T>(1);
    return tab;
}

/// Forward Euler (one-stage). Cheap and simple; first-order accurate in
/// δt. Useful for unit tests and short-step diagnostic runs; not
/// recommended for long propagations. = rk_step(euler_tableau) + LTE.
template<typename T>
dynamics::State<T> euler(const dynamics::State<T>& y0,
                         const math::TrackedValue<T>& dt,
                         const AccelFn<T>& accel_fn) {
    RkStepResult<T, 1> r = rk_step(euler_tableau<T>(), y0, dt, accel_fn);
    dynamics::State<T> next = r.state;

    // R07 / REQ-EF-7: Euler local truncation error O(dt²).
    // The per-step LTE envelope (runge_kutta_lie_group.md §6.2/§6.3):
    // position/rotation (h²/2)·A, velocity/angular-rate h·A, with the stage
    // acceleration magnitudes as the derivative-cascade proxy, split by part
    // and deposited per slot family. Deliberately method-order-blind — with
    // |a| as the only datum, an h^{p+1} claim would be perceived accuracy.
    {
        T a_lin = r.stage_accel[0].l1_norm_linear();
        T a_ang = r.stage_accel[0].l1_norm_angular();
        T h2_half = T(1) / T(2) * dt.value * dt.value;
        next.add_step_lte(h2_half * a_lin, h2_half * a_ang,
                          dt.value * a_lin, dt.value * a_ang);
    }
    return next;
}

/// Classical Runge–Kutta 4 (four-stage). Fourth-order accurate in δt.
/// = rk_step(rk4_tableau) + LTE. The generic combine Σ bᵢ Wᵢ reorders the
/// hand-factored ⅙(W₁+2W₂+2W₃+W₄), so the result is round-off-equal (not
/// bit) to the prior form — see runge_kutta_lie_group.md §4.2.
template<typename T>
dynamics::State<T> runge_kutta_4(const dynamics::State<T>& y0,
                                 const math::TrackedValue<T>& dt,
                                 const AccelFn<T>& accel_fn) {
    RkStepResult<T, 4> r = rk_step(rk4_tableau<T>(), y0, dt, accel_fn);
    dynamics::State<T> next = r.state;

    // R07 / REQ-EF-7: the per-step LTE envelope (runge_kutta_lie_group.md
    // §6.2/§6.3): position/rotation (h²/2)·A, velocity/angular-rate h·A, the
    // max-over-stages acceleration magnitudes as the derivative-cascade
    // proxy, split by part, per slot family. Method-order-blind by design —
    // RK4's true O(h⁵) superiority is measured in the VALUES (gate RK1
    // empirical order 3.998); claiming it in the bound from |a| alone would
    // be perceived accuracy. The order-true estimate is the adaptive path's
    // embedded difference.
    {
        using std::max;
        T a_lin = max(
            max(r.stage_accel[0].l1_norm_linear(), r.stage_accel[1].l1_norm_linear()),
            max(r.stage_accel[2].l1_norm_linear(), r.stage_accel[3].l1_norm_linear()));
        T a_ang = max(
            max(r.stage_accel[0].l1_norm_angular(), r.stage_accel[1].l1_norm_angular()),
            max(r.stage_accel[2].l1_norm_angular(), r.stage_accel[3].l1_norm_angular()));
        T h2_half = T(1) / T(2) * dt.value * dt.value;
        next.add_step_lte(h2_half * a_lin, h2_half * a_ang,
                          dt.value * a_lin, dt.value * a_ang);
    }
    return next;
}

} // namespace integrators
