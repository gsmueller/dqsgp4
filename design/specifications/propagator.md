# Specification: Propagator

The `dynamics::Propagator<T>` class orchestrates state advancement: given
an initial state, a body inertia, a list of force lambdas, an integrator,
and a constants provider, it produces a state at a later time.

Each requirement has a stable identifier (`REQ-PR-N`). Tests live in
`tests/test_propagator.cpp`.


## REQ-PR-1 — Pure orchestrator, no embedded physics

Verifies: OBJ-1, OBJ-4, CON-3.

The propagator class contains NO force model, NO hardcoded numerical
constants, NO time-stepping formula. Every dependency arrives through the
constructor:

    Propagator<T>(K, inertia, forces, integrator)

Algorithmic code in `Propagator` reads `K`, dispatches `forces`, applies
`inertia.acceleration_from_wrench`, and delegates timestepping to
`integrator`. The body of `Propagator` does no physics of its own.


## REQ-PR-2 — Body acceleration = sum of forces, divided by inertia

Verifies: OBJ-1.

For state y, the body acceleration twist is

    a(y) = inertia⁻¹ ( Σ force_i(y, K) )

with the sum taken over every force lambda in `forces`. The inertia
applies its inverse via `acceleration_from_wrench`; the propagator does
not implement physics on its own.


## REQ-PR-3 — Single step delegates to the injected integrator

Verifies: OBJ-1, OBJ-4.

The `step(y0, dt)` method builds a body-acceleration callback closing
over the propagator and passes it to `integrator(y0, dt, accel_fn)`. The
return value is the integrator's output verbatim. No post-processing.


## REQ-PR-4 — Multi-step propagation to a target time

Verifies: OBJ-1.

`propagate_to(y0, t_target, dt_max)` advances the state in steps of size
at most `dt_max` until the elapsed time reaches `t_target`. The final
step is shortened so the time lands at `t_target` (subject to
floating-point representation error). The function returns the final
state.

Precondition: `t_target.value > y0.time.value`.


## REQ-PR-5 — Force-lambda signature is (State, ConstantsProvider) → Wrench

Verifies: OBJ-4, OBJ-5.

Force lambdas in `forces` have signature

    Wrench<T> (const State<T>&, const ConstantsProvider<T>&)

The propagator passes its own `K_` member to each lambda at evaluation
time. Lambdas do not capture `K` themselves; this avoids lifetime
hazards when the caller's `ConstantsProvider` is constructed locally and
goes out of scope after the propagator is built.


## REQ-PR-6 — Once-only construction, repeated use

Verifies: feedback_compute_once, CON-2.

The propagator stores `K`, `inertia`, `forces`, and `integrator` at
construction (`std::move` into members). On every `step` and
`propagate_to` call, the stored objects are read; no re-derivation
occurs.


## REQ-PR-7 — Read-only access to constants and inertia

Verifies: REQ-CP-1.

`constants()`, `inertia()`, and `forces()` return const references to the
stored objects. Callers may inspect Earth parameters, the body inertia,
or the configured force list, but cannot mutate them.


## REQ-PR-8 — Three-error propagation through every advancement

Verifies: OBJ-3, CON-5.

Every state advancement composes `TrackedValue<T>` arithmetic from the
math layer. The output state's `errors.measurement`, `errors.precision`,
and `errors.accuracy` budgets reflect all sources accumulated to that
point. No category is dropped or silently zeroed.

Model truncation entering through a force lambda (e.g., J₂-only gravity
representing a higher-order field) appears in `errors.accuracy`
(REQ-EF-7); numerical truncation (Taylor branches, series, retraction
projections) appears in `errors.precision` (REQ-EF-6, REQ-EF-15).


## REQ-PR-9 — Pose retraction every step

Verifies: REQ-EF-15.

The integrator's `lie_advance_pose` applies `Pose::normalized()` after
each pose update, projecting any sub-ULP drift back onto the SE(3)
manifold. The propagator does not perform a second retraction; doing so
would double-count the retraction precision.


## REQ-PR-10 — Deterministic, no global state

Verifies: OBJ-8, CON-2.

`step` and `propagate_to` are pure functions of their inputs and the
propagator's stored members. No randomness, no environment reads, no
global mutable state, no platform-dependent behavior. Two propagators
constructed with equal arguments produce equal outputs from equal
inputs.


## REQ-PR-11 — No force lambda is special-cased

Verifies: OBJ-4, CON-3.

The propagator iterates `forces` uniformly with `force(state, K)` and
sums the resulting wrenches. The first lambda, the last lambda, and a
J₂ lambda inserted between them all receive identical treatment. There
is no internal switch on "gravity vs drag vs third-body".


## REQ-PR-12 — Composable integrator

Verifies: OBJ-4.

`integrator_` is an `IntegratorFn<T>` — a `std::function` matching the
shared integrator signature. Callers may supply Euler, RK4, RKF7(8), a
symplectic scheme, or a custom test stub. The propagator does not
inspect which integrator is in use.


## Cross-reference

Audited by `design/audit/code_consistency.md` (style) and the planned
`tests/test_propagator.cpp`:

- LEO smoke test against the analytical Kepler orbit (REQ-PR-2..4,
  REQ-PR-9).
- Energy- and angular-momentum-conservation tests (REQ-PR-8).
- Force-list permutation invariance (REQ-PR-11).
- Integrator swap invariance (REQ-PR-12).
- Three-error bound dominance over actual numerical error
  (REQ-EF-2, REQ-PR-8).
