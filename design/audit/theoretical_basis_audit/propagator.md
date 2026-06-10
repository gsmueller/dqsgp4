# Theoretical Basis Audit — `src/dynamics/propagator.h`

**File**: `src/dynamics/propagator.h`  
**Lines**: 150  
**Expected formula count**: 7 (ctor + 3 accessors + compute_acceleration + step + propagate_to)  
**Audit date**: 2026-05-13

---

## FORMULA 1: Constructor initialization

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::Propagator::ctor
Location:               src/dynamics/propagator.h:80-87
Mathematical statement: Construct a Propagator by binding the constants provider K,
                        body inertia, force list, and integrator as members via move
                        semantics (pure storage, no computation).

THEORY
  Underlying theorem:   Object construction is initialization; no numeric theorem
                        applies — only object lifetime and move semantics from C++
                        language rules (REQ-PR-6, REQ-PR-1).
  Primary reference:    C++ standard (ISO/IEC 14882) §9.4.1 constructors, §20.2.1
                        move semantics.
  Domain of validity:   All input types matching the signatures
                        (ConstantsProvider<T>, Inertia<T>, vector<ForceFn<T>>,
                        IntegratorFn<T>).

METHOD
  Method declared:      Memberwise initialization via std::move (no physics).
  Method implemented:   Initializer list binds K_, inertia_, forces_, integrator_
                        from the constructor arguments (lines 84-87). No
                        computation.
  Match verdict:        ✓ matched — pure storage operation, no theory required.

ERROR BOUND
  Bound category:       n/a (non-numeric initialization)
  Bound formula:        Not applicable.
  Bound implemented:    Not applicable.
  Bound verdict:        ✓ N/A — constructor stores references; error propagation
                        occurs downstream when stored members are used.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (TrackedValue composability); errors in K, inertia
                        propagate downstream on use.
  AUD-EF applies:       n/a (no numeric operation)
  AUD-MC applies:       n/a
  Verification test:    REQ-PR-6 (once-only construction), REQ-PR-1 (pure
                        orchestrator with no embedded physics).

NOTES
  - Constructor is responsible for storing all dependencies at once; once
    stored, they are immutable (const references via accessors). This satisfies
    REQ-PR-6 (once-only construction, repeated use).
  - K_ is moved into the propagator so force lambdas receive a stable reference
    (REQ-PR-5) and do not capture K themselves, avoiding lifetime hazards.
`

---

## FORMULA 2: `constants()` accessor

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::constants
Location:               src/dynamics/propagator.h:92
Mathematical statement: Return a const reference to the stored constants provider.

THEORY
  Underlying theorem:   Object access and const-correctness; no numeric theorem.
  Primary reference:    C++ const-reference semantics (ISO/IEC 14882) §9.3.2.
  Domain of validity:   All program points after construction.

METHOD
  Method declared:      Return const reference; no modification.
  Method implemented:   `return K_;` (line 92).
  Match verdict:        ✓ matched — pure read operation.

ERROR BOUND
  Bound category:       n/a (non-numeric read)
  Bound formula:        Not applicable.
  Bound implemented:    Not applicable.
  Bound verdict:        ✓ N/A — accessor does not modify or compute; error state
                        of K_ is unaffected.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (errors in K propagate through downstream reads).
  AUD-EF applies:       n/a
  AUD-MC applies:       n/a
  Verification test:    REQ-PR-7 (read-only access to constants).

NOTES
  - This accessor satisfies REQ-PR-7 (read-only access) and supports callers
    inspecting Earth parameters without mutating the propagator state.
`

---

## FORMULA 3: `inertia()` accessor

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::inertia
Location:               src/dynamics/propagator.h:95
Mathematical statement: Return a const reference to the stored body inertia.

THEORY
  Underlying theorem:   Object access; no numeric theorem.
  Primary reference:    C++ const-reference semantics.
  Domain of validity:   All program points after construction.

METHOD
  Method declared:      Return const reference; no modification.
  Method implemented:   `return inertia_;` (line 95).
  Match verdict:        ✓ matched — pure read operation.

ERROR BOUND
  Bound category:       n/a (non-numeric read)
  Bound formula:        Not applicable.
  Bound implemented:    Not applicable.
  Bound verdict:        ✓ N/A — accessor does not modify or compute.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (errors in inertia propagate downstream).
  AUD-EF applies:       n/a
  AUD-MC applies:       n/a
  Verification test:    REQ-PR-7.

NOTES
  - Satisfies REQ-PR-7 (read-only access to inertia). Used by callers to inspect
    the configured body properties.
`

---

## FORMULA 4: `forces()` accessor

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::forces
Location:               src/dynamics/propagator.h:98
Mathematical statement: Return a const reference to the stored force-lambda list.

THEORY
  Underlying theorem:   Object access; no numeric theorem.
  Primary reference:    C++ const-reference semantics.
  Domain of validity:   All program points after construction.

METHOD
  Method declared:      Return const reference to the force list; no modification.
  Method implemented:   `return forces_;` (line 98).
  Match verdict:        ✓ matched — pure read operation.

ERROR BOUND
  Bound category:       n/a (non-numeric read)
  Bound formula:        Not applicable.
  Bound implemented:    Not applicable.
  Bound verdict:        ✓ N/A — accessor does not modify or compute.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (errors in force outputs propagate from each
                        lambda's return wrench).
  AUD-EF applies:       n/a
  AUD-MC applies:       n/a
  Verification test:    REQ-PR-7, REQ-PR-11 (force-list permutation invariance).

NOTES
  - Satisfies REQ-PR-7 (read-only access) and supports REQ-PR-11 verification
    that no force lambda is special-cased.
`

---

## FORMULA 5: `compute_acceleration(state)`

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::compute_acceleration
Location:               src/dynamics/propagator.h:104-110
Mathematical statement: Body acceleration = (Σ_i force_i(state, K)) / inertia
                        = inertia⁻¹ · (Σ force),
                        where sum is over all force lambdas in the propagator.

THEORY
  Underlying theorem:   Newton's second law: F_total = m · a, hence
                        a = F_total / m. For a rigid body (6-DOF):
                        a_twist = inertia⁻¹ · total_wrench
                        (spatial momentum / mass → body acceleration).
  Primary reference:    Classical mechanics textbook; relevant to propagator
                        specs: REQ-PR-2 (body acceleration = sum of forces,
                        divided by inertia).
  Domain of validity:   All states reachable by the propagator; each force
                        lambda must return a well-defined Wrench<T> at the given
                        state.

METHOD
  Method declared:      Sum all force lambdas (via iteration over forces_),
                        divide the total wrench by inertia via the method
                        `inertia_.acceleration_from_wrench(total)`.
  Method implemented:   Lines 105-109:
                          (1) Initialize Wrench<T> total = Wrench<T>::zero()
                          (2) For each force lambda f in forces_:
                              total = total + f(state, K_)
                          (3) Return inertia_.acceleration_from_wrench(total)
  Match verdict:        ✓ matched — Newton's second law applied as F/m, with
                        iteration summing each force's contribution.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (composed from three
                        sources per REQ-EF-3 triangle inequality)
  Bound formula:        Wrench composition via TrackedValue<T> arithmetic
                        (line 107: total = total + f(...)) propagates errors
                        from each force wrench. Inertia division
                        (line 109: acceleration_from_wrench) applies the
                        matrix-division bound from Inertia<T>.
                        Total error = Σ (errors from force_i) + error from
                        inertia⁻¹ application, per REQ-EF-3 (closure under
                        addition and division).
  Bound implemented:    Implicit in TrackedValue<T> operator overloads:
                          - Wrench<T>::operator+ adds errors from both operands
                            (measurement, precision, accuracy budgets).
                          - Inertia<T>::acceleration_from_wrench applies the
                            inverse and propagates errors.
                        The final Twist<T> contains the composed bounds (lines
                        107, 109 use overloaded +, no explicit bound added here).
  Bound verdict:        ✓ matched — error composition is implicit in the
                        TrackedValue<T> operator overloads; all three error
                        categories propagate per REQ-PR-8 (three-error
                        propagation through every advancement).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (error closure under vector sum and division).
  AUD-EF applies:       AUD-EF-1 (every public op returns TrackedValue<T>);
                        AUD-EF-6 (wrench arithmetic and inertia division add
                        bounds via operator overloads).
  AUD-MC applies:       No algebra axioms tested here (composition of known
                        correct ops).
  Verification test:    REQ-PR-2 (body acceleration = sum of forces / inertia).
                        Energy and angular-momentum conservation tests verify
                        the computed acceleration drives correct dynamics
                        (tests/test_propagator.cpp).

NOTES
  - This is an orchestration point: the propagator does not implement force
    summation or inertia division itself. It relies on Wrench<T>::operator+ and
    Inertia<T>::acceleration_from_wrench to handle the mathematics and error
    propagation.
  - Force lambdas are evaluated at every stage of the integrator's RK4 tableau;
    each evaluation contributes to the error budget independently.
  - The iteration over forces_ is uniform (REQ-PR-11): no force is special-cased,
    and the order of summation does not affect the result (commutativity of
    wrench addition).
`

---

## FORMULA 6: `step(y0, dt)`

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::step
Location:               src/dynamics/propagator.h:116-123
Mathematical statement: Single integrator step: y₁ = Integrator(y₀, dt, accel_fn).

THEORY
  Underlying theorem:   Lie-group numerical integration. The integrator (RK4 or
                        other) is a Runge–Kutta scheme applied to the state ODE:
                          dy/dt = g(t, y)
                        where g(t, y) = (pose, twist, accel(y)).
                        The integrator's local truncation error is O(dt^p)
                        (p = 4 for RK4 per REQ-IN-4). The integrator's
                        lie_advance_pose applies the SE(3) screw exponential
                        and retraction (Pose::normalized) per REQ-IN-2.
  Primary reference:    REQ-IN-1..8 (integrator specifications); Munthe-Kaas
                        (1999) "High order Lie group methods" for the LG-RK4
                        theory. Butcher tableaux for RK4: Hairer et al. (1987)
                        "Solving Ordinary Differential Equations I".
  Domain of validity:   All states y0 reachable during propagation; dt must be
                        positive and small enough for the integrator's
                        convergence region.

METHOD
  Method declared:      Create a closure accel_fn over the propagator's
                        compute_acceleration method, then invoke
                        integrator_(y0, dt, accel_fn) and return the result
                        verbatim (REQ-PR-3).
  Method implemented:   Lines 118-122:
                          (1) Define accel_fn as a lambda capturing `this`,
                              which calls compute_acceleration(state).
                          (2) Return integrator_(y0, dt, accel_fn).
                        The integrator is an IntegratorFn<T> injected at
                        construction (line 87: integrator_).
  Match verdict:        ✓ matched — step creates the acceleration callback and
                        delegates to the integrator function, applying whichever
                        RK4 or other scheme is plugged in (REQ-PR-3,
                        REQ-PR-12).

ERROR BOUND
  Bound category:       precision, accuracy, measurement (composed)
  Bound formula:        For RK4 (Runge–Kutta order 4), the local truncation
                        error per step is O(dt^5). The integrator adds this
                        bound to the output state's precision budget (REQ-IN-8).
                        The acceleration callback's error (from compute_acceleration)
                        is evaluated at up to 4 RK4 stages; each stage's errors
                        compose. The final bound is the sum of errors from all
                        stages plus the RK4 local truncation bound (added by the
                        integrator; the propagator does not add extra bounds here).
                        Pose retraction (Pose::normalized) adds a separate
                        precision bound for SE(3) manifold projection (REQ-EF-15,
                        REQ-IN-7).
  Bound implemented:    Implicit in the integrator's return: the integrator is
                        responsible for returning a State<T> whose errors
                        accumulate the RK4 truncation bound and all
                        acceleration-callback errors. The propagator does not
                        post-process; it returns the integrator's output verbatim
                        (line 122: return integrator_(...)).
  Bound verdict:        ✓ matched — error propagation is delegated to the
                        integrator function, which is responsible for REQ-IN-8
                        (three-error propagation through every advancement).
                        The propagator's role is orchestration (no extra bound
                        added).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-6 (truncation bounds from RK4), REQ-EF-15
                        (retraction projection).
  AUD-EF applies:       AUD-EF-1 (step returns TrackedValue-wrapped state),
                        AUD-EF-7 (integrator adds appropriate bounds).
  AUD-MC applies:       Tested by LEO orbit closure (tests/test_propagator.cpp),
                        which verifies RK4's O(dt^4) convergence rate.
  Verification test:    REQ-PR-3 (single step delegates to integrator); REQ-IN-4
                        (order-of-accuracy via Richardson extrapolation on
                        step size).

NOTES
  - The step function is a thin wrapper. All numerical work is done by the
    injected integrator. This supports REQ-PR-1 (pure orchestrator, no embedded
    physics) and REQ-PR-12 (composable integrator).
  - The closure accel_fn captures `this`, allowing the lambda to call
    compute_acceleration; this is a standard C++ pattern and imposes no
    computational overhead (inline by the compiler).
  - No second retraction is applied here (per REQ-PR-9: the integrator already
    calls Pose::normalized after each pose advance; calling it again would
    double-count the retraction error).
`

---

## FORMULA 7: `propagate_to(y0, t_target, dt_max)`

`
=== FORMULA AUDIT CARD ===
ID:                     propagator::propagate_to
Location:               src/dynamics/propagator.h:130-141
Mathematical statement: Multi-step propagation from time y0.time to t_target,
                        taking steps of size at most dt_max, with the final step
                        shortened to land exactly at t_target.

THEORY
  Underlying theorem:   Composition of single-step integrations. If each step
                        solves the local ODE to order p (REQ-IN-4), then N steps
                        with step sizes summing to T compose to global error
                        O(N · dt^p) ≈ O(T · dt^(p-1)) for fixed step sizes.
                        Variable step sizes require a tighter analysis:
                        global error remains O(T · dt^(p-1)) if all steps are
                        comparable in size (no extreme variation).
                        The final shortened step (line 137) introduces no extra
                        error: it is still integrated to order p by the injected
                        integrator.
  Primary reference:    Hairer et al. (1987) §III.1 (global error of RK methods).
                        Practical propagation to a target time: REQ-PR-4.
  Domain of validity:   t_target > y0.time (precondition line 129);
                        all intermediate states reachable via single steps.

METHOD
  Method declared:      Loop while y.time < t_target:
                          (1) Compute remaining time: remaining = t_target - y.time
                          (2) Choose step size: min(remaining, dt_max)
                          (3) Advance one step: y = step(y, step_size)
                          (4) Repeat until y.time ≥ t_target
                        Return the final state.
  Method implemented:   Lines 133-140:
                          - y = y0 (initialize to the input state)
                          - while (y.time.value < t_target.value) loop:
                              - remaining = t_target - y.time (TrackedValue arithmetic)
                              - step_size = (remaining.value < dt_max.value) ?
                                remaining : dt_max (ternary; applies TrackedValue
                                comparison and selection)
                              - y = step(y, step_size) (invoke Formula 6)
                          - return y (line 140)
  Match verdict:        ✓ matched — loop-and-step composition of single
                        integrator steps (REQ-PR-4).

ERROR BOUND
  Bound category:       precision, accuracy, measurement (composed)
  Bound formula:        Global error from N steps to time T ≈ t_target - y0.time:
                          |error| ≤ C · N · dt^(p-1) = C · (T / dt_avg) · dt^(p-1)
                                  = C · T · dt^(p-2)
                        For RK4 (p=4): |error| ≈ C · T · dt².
                        Each step's errors are composed via TrackedValue<T>
                        arithmetic (line 138: y = step(...)). The output state's
                        error budgets accumulate: errors.precision contains all
                        truncation contributions, errors.accuracy contains all
                        model-truncation contributions, errors.measurement
                        contains all measurement contributions.
                        The output state's error.precision is the sum of all N
                        step errors (each step adds O(dt^5) for RK4, summing to
                        O(N · dt^5) = O(T · dt^4) for constant dt).
  Bound implemented:    Implicit in the loop: each step(y, step_size) returns a
                        state with accumulated errors. The TrackedValue arithmetic
                        in the loop condition and step-size computation (lines 135-137)
                        preserves and propagates the error budgets. The final y
                        returned has error budgets reflecting all N accumulated steps.
                        No explicit bound is added by propagate_to itself; all
                        bounds flow from step() (Formula 6) → integrator →
                        acceleration callback (Formula 5).
  Bound verdict:        ✓ matched — error accumulation follows from composing
                        N single-step bounds. REQ-PR-8 guarantees all three error
                        categories propagate through the loop (TrackedValue
                        arithmetic in lines 135-138 uses overloaded operators).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (error closure under TrackedValue subtraction
                        and comparison), REQ-EF-6 (truncation from RK4 summed
                        over steps).
  AUD-EF applies:       AUD-EF-1 (propagate_to returns TrackedValue-wrapped
                        state), AUD-EF-6 (each step's bounds are accumulated).
  AUD-MC applies:       Tested by LEO orbit closure (tests/test_propagator.cpp):
                        propagate_to over 1 orbital period should reproduce the
                        initial state to within the predicted error bound.
  Verification test:    REQ-PR-4 (multi-step to target time, final step shortened);
                        energy and angular-momentum conservation (tests/test_propagator.cpp);
                        precondition check: t_target > y0.time.

NOTES
  - The loop is simple and terminating: on each iteration, y.time increases by
    step_size > 0 (assuming the integrator advances time correctly), so the loop
    condition eventually becomes false.
  - The final step is shortened (line 137) so the elapsed time lands exactly at
    t_target (subject to floating-point representation error). No "overshooting"
    occurs.
  - Variable step sizes (remaining < dt_max on the final iteration) do not change
    the order of accuracy; the integrator still returns a state with O(dt^p)
    local truncation error. Composition to global error remains O(T · dt^(p-1)).
  - All error composition occurs via TrackedValue<T> operator overloads; no
    explicit summation formula is needed in the propagate_to code.
`

---

## File-level verdict

| Dimension | Result | Summary |
|---|---|---|
| **A. Error wiring** | ✓ PASS | All functions that produce TrackedValue results compose errors implicitly via operator overloads. compute_acceleration (Formula 5), step (Formula 6), and propagate_to (Formula 7) all return states with error budgets set by the integrator and acceleration callback. |
| **B. Algebra axioms** | ✓ PASS (by composition) | No algebraic axioms are verified here; this is an orchestration class. Correctness is delegated to Wrench, Inertia, Twist, Pose, and State classes. Tested by LEO orbit closure (tests/test_propagator.cpp). |
| **C. Theoretical basis** | ✓ PASS | All formulas are identified and their theory-method-bound triples match. (1) Constructor: pure storage. (2–4) Accessors: pure reads. (5) compute_acceleration: Newton's F/m law. (6) step: RK4 Lie-group integration. (7) propagate_to: composition of steps with final step shortened to target time. |

**Overall file verdict: PASS** — All 7 functions match their cited theory and error-propagation model. The propagator is a pure orchestrator with no embedded physics, and all numeric errors are composed via delegated classes (integrator, forces, inertia) per REQ-PR-1 and REQ-EF-3.
