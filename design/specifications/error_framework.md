# Specification: Error Framework

The binding contract for OBJ-3 (three-error tracking) and CON-4, CON-5, CON-6
(rigorous bounds; mandatory tracking; truncations accounted). Every operation in
the library is governed by these requirements.

Each requirement has a stable identifier (`REQ-EF-N`). Tests in
`tests/audit/test_error_framework.cpp` cite these identifiers in their names.


## REQ-EF-1 — Every tracked value carries three error budgets

Verifies: OBJ-3, CON-5.

Every `TrackedValue<T>` carries a `ThreeErrors<T>` aggregate with non-negative
components:

    measurement — propagated from physical input uncertainty
    precision   — accumulated numerical roundoff, truncation, and residual
    accuracy    — accumulated model truncation and simplification

The three budgets are independent. Operations may not merge them and may not drop
any one.


## REQ-EF-2 — `total_error()` is a rigorous upper bound

Verifies: CON-4.

For any `TrackedValue<T> x`,

    |x.value − x_true| ≤ x.total_error()

where x_true is the mathematical value computed in infinite precision under the
same model. This is an upper bound, not an estimate. The triangle-inequality
combination (sum of categories) is correct: the three error sources are not
assumed independent.


## REQ-EF-3 — Closed-form arithmetic propagates errors per category

Verifies: OBJ-3, CON-5.

For binary operators (+, −, ×, ÷) and unary functions
(sqrt, sin, cos, atan, atan2, abs, unary −, fmod), the output bound is derived
analytically from the operation's properties and applied per category. The
per-category formulas are documented in `src/math/tracked_value.h`:

    add/sub  : bound = a_bound + b_bound
    mul      : bound = |a_val|·b_bound + |b_val|·a_bound + a_bound·b_bound
    div      : bound = (|a_val|·b_bound + |b_val|·a_bound) /
                       (|b_val|·(|b_val| − b_bound))
    sqrt     : bound = err / (2·sqrt(val − err))            when err < val
    sin/cos  : bound = |deriv(val)|·err + err²/2,    capped at 2
    atan     : bound = err / (1 + val²)
    atan2    : bound = (|x_val|·y_err + |y_val|·x_err) / (x_val² + y_val²)

Each formula has a derivation comment at point of use. New closed-form operations
in the library must derive and document their per-category bound in the same way.


## REQ-EF-4 — Truncation of a series is added to precision

Verifies: OBJ-3, CON-6.

When a series is evaluated to N terms with the tail bounded by T_N (Leibniz bound,
geometric tail bound, or another rigorous bound), the result's `errors.precision`
is incremented by T_N:

    result.errors.precision += T_N;

Pattern reference: `alternating_series`, `geometric_series`, `series_sqrt` in
`src/math/series.h`; `make_binomial_evaluator` in `src/math/binomial_series.h`.


## REQ-EF-5 — Convergence residual is added to precision

Verifies: OBJ-3, CON-6.

When an iterative algorithm (Newton, Halley, fixed-point) terminates at iteration
k with correction Δ_k below the caller's tolerance, the result's
`errors.precision` is incremented by |Δ_k|:

    result.errors.precision += abs(delta_k);

Pattern reference: `solve_kepler` in `src/math/kepler.h`.


## REQ-EF-6 — Fixed-order Taylor truncation is added to precision

Verifies: OBJ-3, CON-6.

When a function evaluates a fixed-order Taylor expansion in a small-argument
branch (e.g., sinc near zero, atan2/|v| near zero), the truncation bound — the
magnitude of the first omitted term, or a larger rigorous bound — is added to
`errors.precision`:

    result.errors.precision += truncation_bound;

This applies to every fixed-order Taylor branch added in `src/math/quaternion.h`,
`src/math/dual_quaternion.h`, and future modules.


## REQ-EF-7 — Model truncation is added to accuracy

Verifies: CON-6.

When a model is deliberately simplified (J₂-only gravity vs full spherical
harmonics, exponential drag vs an empirical density model, two-body astronomy vs
JPL DE ephemeris, etc.), the model-truncation bound is added to `errors.accuracy`.
Force lambdas in `src/forces/` declare their model-truncation contribution; the
propagator accumulates it into the state at each step.


## REQ-EF-8 — Errors never silently decrease

Verifies: CON-4.

For any non-identity operation f, the output's `total_error()` is at least the
maximum input `total_error()` scaled by the operation's Lipschitz constant.
Operations may not reduce reported error below this lower bound.

Operations may preserve error exactly when the Lipschitz constant is 1
(conjugate, unary negation, `abs`). They may not preserve a smaller error.


## REQ-EF-9 — Catastrophic loss is signaled, not silenced

Verifies: CON-4.

When an analytical error formula fails — division by a quantity smaller than its
error, sqrt of a value smaller than its error, atan2 inside the singular disc —
the result's error is set to a conservative bound (the full magnitude of the
result, or `std::numeric_limits<T>::max()`). Per-operator behavior is documented
in `src/math/tracked_value.h`. The library does not silently return tight bounds
in regions where the analytical formula is invalid.


## REQ-EF-10 — `value` is a plain T

Verifies: CON-7.

`TrackedValue<T>::value` is the unmodified scalar of type T. It does not itself
carry an error budget; this permits comparisons (`x.value < 0`) and ordinary
arithmetic on the raw value where a tracked path is not needed.


## REQ-EF-11 — Reliable-digits query

Verifies: OBJ-3.

`TrackedValue<T>::reliable_digits()` returns the largest integer k such that
10^(−k) ≥ total_error() / |value|. Callers use this to decide whether a quantity
is reliable enough for downstream use without inspecting the three categories
individually. Returns `INT_MAX` for an exactly represented value with zero total
error; returns 0 when error dominates the value.


## REQ-EF-12 — Composite types inherit the framework by composition

Verifies: OBJ-3, CON-5, CON-13.

Composite types — `Vector3<T>`, `Quaternion<T>`, `DualNumber<T>`,
`DualQuaternion<T>`, `Pose<T>`, `Twist<T>`, `Wrench<T>`, `Inertia<T>`,
`State<T>` — store `TrackedValue<T>` components and inherit error propagation
by construction. They must not duplicate the propagation logic; they compose it
via the existing `TrackedValue` operator overloads.


## REQ-EF-13 — Composite error queries

Verifies: OBJ-3.

A composite type exposes per-component error access at its public surface (e.g.,
`pose.position().x.errors.precision`). It may also expose an aggregate
`total_error()` defined as the maximum per-component `total_error()`. Aggregates
do not replace per-component access.


## REQ-EF-14 — Identity operations preserve error exactly

Verifies: CON-4.

Operations whose mathematical effect is identity on a subspace — e.g., conjugate
on a unit quaternion's real part, normalization of an already-unit quaternion —
preserve the error budget exactly. They do not introduce additional precision
error beyond the unavoidable floating-point representation cost (which is itself
recorded if the operation rebuilds the value).


## REQ-EF-15 — Retraction onto a constraint surface accounts for retraction error

Verifies: CON-6.

When a value is projected onto a constraint surface (unit-norm renormalization, an
SE(3) retraction, etc.), the discarded normal-direction component is added to
`errors.precision`:

    result.errors.precision += |projected_off_amount|;

Pattern reference: post-step unit-DQ retraction in the propagator's integration
loop.


## Audit cross-reference

This specification is audited by:

- `design/audit/error_framework.md` — checklist form, applied during code review
  and as a CI gate.
- `tests/audit/test_error_framework.cpp` — automated property-based checks: for
  each operator on each composite type, sample random inputs at known precision,
  compute the reported `total_error()`, compute the actual numerical error
  against a higher-precision reference, and assert dominance in 100 % of trials.
