# Specification: Integrator

The `integrators::` namespace contains time-stepping integrators consumed
by `dynamics::Propagator<T>`. Each integrator is a function matching the
signature

    State<T>(const State<T>& y0,
             const TrackedValue<T>& dt,
             const AccelFn<T>& acc_fn)

where `acc_fn(state)` returns the body acceleration `Twist<T>` at the
given state. Integrators are injected into the propagator at construction
(OBJ-4, CON-3, REQ-PR-12).

Each requirement has a stable identifier `REQ-IN-N`. Tests live in
`tests/test_propagator.cpp` (Lie-group RK4 on a closed-form Kepler orbit)
and forthcoming test programs for new integrators.


## REQ-IN-1 — Common signature

Verifies: OBJ-4, REQ-PR-12.

Every integrator matches `integrators::IntegratorFn<T>`:
`(State<T>, TrackedValue<T>, AccelFn<T>) → State<T>`. Adapter helpers may
expose alternate spellings, but the canonical form is the only one the
propagator calls.


## REQ-IN-2 — Pose advances on the Lie group

Verifies: REQ-EF-15.

Pose advance uses the SE(3) screw exponential, never Euclidean addition.
The integrator's helper `lie_advance_pose(pose, dt, twist)` returns

    pose · exp_screw((dt / 2) · twist)

followed by `Pose::normalized()`. This keeps the pose on the unit-DQ
manifold to floating-point precision at every step.


## REQ-IN-3 — Twist advances linearly

Verifies: OBJ-1.

Twist lives in the Lie algebra (a vector space) and advances via vector
arithmetic: `twist_new = twist + dt · acceleration`. Stage combinations
(RK4, RKF7(8), etc.) combine accelerations linearly via
`Twist<T>::operator+` and scalar multiplication.


## REQ-IN-4 — Order claimed by name

Verifies: OBJ-2.

| Name | Order | Stages | Notes |
|---|---|---|---|
| `euler` | 1 | 1 | Diagnostic / unit-test use only |
| `runge_kutta_4` | 4 | 4 | Standard non-adaptive RK4 |
| `rkf78` (forthcoming) | 7/8 | 13 | Embedded error estimate |
| `dormand_prince_8` (future) | 8/7 | 13 | Higher-precision DP variant |
| `symplectic_leapfrog` (future) | 2 | 1 | Conservative systems only |

The integrator's local truncation error is `O(dt^p)` where p is the
claimed order. Global error after N steps is bounded by `O(N · dt^p)` for
non-stiff problems.


## REQ-IN-5 — No global state

Verifies: OBJ-8, CON-2.

An integrator is a pure function. Two calls with equal arguments produce
equal outputs. No internal mutable state, no static variables, no
randomness, no environment reads. Integrators called concurrently on
different states do not interfere.


## REQ-IN-6 — Adaptive integrators expose error feedback

Verifies: OBJ-2.

Adaptive integrators (`rkf78` and later) implement an extended interface
that returns both the advanced state and an embedded-method-derived error
estimate. The propagator uses the estimate to choose the next step size.
Non-adaptive integrators do not implement this interface; the propagator
passes the configured `dt_max` as-is.


## REQ-IN-7 — Pose retraction every step

Verifies: REQ-EF-15, REQ-PR-9.

`Pose::normalized()` is called after every pose advance. The integrator
does not skip retraction even when the underlying drift is sub-ULP —
consistency across step counts and integrator schemes is preferred over
saving one square root and four divisions per step.


## REQ-IN-8 — Three-error propagation

Verifies: OBJ-3, CON-5, REQ-EF-12.

Every `TrackedValue<T>` arithmetic inside the integrator — stage
combinations, pose advance via `exp_screw`, twist update, retraction
projection — propagates the three-error budget through the underlying
operator overloads. The output state's `errors.measurement`,
`errors.precision`, and `errors.accuracy` budgets contain every
contribution accumulated through the integrator with no category dropped.


## REQ-IN-9 — Symplectic schemes for conservative problems

Verifies: OBJ-1.

A future `symplectic_leapfrog<T>` integrator preserves a discrete energy
and angular momentum exactly (up to representational error) for
time-independent conservative problems. Its signature and semantic
obligations match the non-symplectic integrators above; users select
based on the physical problem.


## REQ-IN-10 — Composable with all force lambdas

Verifies: OBJ-4, CON-3, REQ-PR-11.

The integrator treats the acceleration callback as opaque. Any force list
supported by the propagator (gravity_central, gravity_J2, drag,
third_body, …) composes with any integrator. No integrator inspects
which force model is active.


## REQ-IN-11 — Lie-group integrators on a unit-DQ pose

Verifies: REQ-DQ-15..18.

Lie-group integrators (`runge_kutta_4`, the planned `rkf78`) use
`DualQuaternion<T>::exp_screw` and `log_screw` from
`design/specifications/dual_quaternion_algebra.md` REQ-DQ-18. The
retraction is `DualQuaternion<T>::normalized()` (REQ-EF-15). Both are
exercised by `tests/test_dual_quaternion.cpp` (AUD-MC-18) and again by
`tests/test_propagator.cpp` (LEO orbit closure).


## REQ-IN-12 — Stage-count and complexity declared

Verifies: OBJ-7.

Each integrator's file documentation declares its stage count and
per-step complexity (number of acceleration-callback evaluations, number
of pose advances, number of normalizations). Callers can reason about
propagation cost without reading the implementation.


## Cross-reference

Audited by `design/audit/code_consistency.md` (style) and the tests:

- `tests/test_propagator.cpp` — exercises `runge_kutta_4` on a circular
  LEO orbit (REQ-IN-1..8, REQ-IN-11).
- Forthcoming `tests/test_integrator_convergence.cpp` — measures
  empirical order-of-accuracy by halving dt and asserting the expected
  Richardson ratio (REQ-IN-4).
- Forthcoming `tests/test_rkf78_adaptive.cpp` — verifies adaptive
  step-size control and embedded error feedback (REQ-IN-6).
