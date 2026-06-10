# Audit: Error Framework

Procedural checklist verifying conformance to
`design/specifications/error_framework.md`. Every file in `src/math/`,
`src/dynamics/`, `src/constants/`, `src/forces/`, `src/integrators/` must pass.

Each item has a stable identifier (`AUD-EF-N`) and cites the requirement it
verifies. Enforced by `tests/audit/test_error_framework.cpp` and by manual
review at PR time.


## AUD-EF-1 — Public operations return tracked values

Verifies: REQ-EF-1.

Every public function or operator that produces a numeric output returns a
`TrackedValue<T>` or a composite of `TrackedValue<T>`. No public function
returns a raw `T` except for the documented `value` accessors (REQ-EF-10).


## AUD-EF-2 — Closed-form operators carry per-category derivations

Verifies: REQ-EF-3.

Every closed-form operator added to the library has a derivation comment
documenting the per-category error bound formula. The bound is implemented in
code (not deferred to a runtime cross-check).


## AUD-EF-3 — Series functions add the tail bound to precision

Verifies: REQ-EF-4.

Every function that evaluates a series ends with

    result.errors.precision = result.errors.precision + T_N;

(or an equivalent in-place form) before returning. The function signature
accepts a tolerance parameter; iteration terminates when the tail bound drops
below tolerance.


## AUD-EF-4 — Iterative algorithms add the residual to precision

Verifies: REQ-EF-5.

Every iterative function (Newton, Halley, fixed-point) records its final
correction magnitude in the result's precision error:

    result.errors.precision = result.errors.precision + abs(delta_k);


## AUD-EF-5 — Taylor branches add the truncation bound to precision

Verifies: REQ-EF-6.

Every fixed-order Taylor branch ends with

    result.errors.precision = result.errors.precision + truncation_bound;

A comment on the surrounding lines identifies the first omitted term (or a
larger rigorous bound) and explains why the chosen threshold is below the
target type's representational precision.


## AUD-EF-6 — Force models add the model bound to accuracy

Verifies: REQ-EF-7.

Force lambdas in `src/forces/` declare their model-truncation contribution; the
propagator accumulates it into the state's `errors.accuracy` at each integrator
step.


## AUD-EF-7 — Composite types compose, never duplicate

Verifies: REQ-EF-12.

Composite types (`Vector3`, `Quaternion`, `DualNumber`, `DualQuaternion`,
`Pose`, `Twist`, `Wrench`, `Inertia`, `State`) store `TrackedValue<T>` fields.
They do not duplicate error-propagation logic; they compose it via the existing
`TrackedValue` operator overloads. No composite stores a raw `T` for any
field that participates in arithmetic.


## AUD-EF-8 — No silent error reduction

Verifies: REQ-EF-8.

A grep audit identifies every assignment to `errors.measurement`,
`errors.precision`, or `errors.accuracy`. Each such assignment is either:

- an addition (`errors.X = errors.X + bound`) — recording new error sources, or
- a per-category propagation through a documented formula.

Direct overwrites that reduce error (`errors.X = T(0)` and the like) are
prohibited unless guarded by a comment justifying why the operation is
mathematically an identity in that category.


## AUD-EF-9 — Catastrophic regions documented

Verifies: REQ-EF-9.

Operations with conditional catastrophic bounds (sqrt of a value smaller than
its error, division by a quantity smaller than its error, atan2 in the singular
disc, etc.) have inline comments naming the singular condition and stating the
conservative bound returned. The bound is documented to dominate the true error
in the singular region.


## AUD-EF-10 — Property-based tests exist

Verifies: every REQ-EF-N.

`tests/audit/test_error_framework.cpp` exists, compiles, and runs without
external dependencies. For each public operation on each composite type it:

1. Samples random inputs at known precision.
2. Computes the reported `total_error()`.
3. Computes the actual error against a higher-precision reference.
4. Asserts `total_error() ≥ actual_error` in 100 % of trials.

The test exits non-zero on any violation.


## Cross-reference

| AUD-EF-N | Verifies | Verified by |
|---|---|---|
| AUD-EF-1 | REQ-EF-1 | `test_error_framework.cpp` |
| AUD-EF-2 | REQ-EF-3 | `test_error_framework.cpp` |
| AUD-EF-3 | REQ-EF-4 | `test_error_framework.cpp` |
| AUD-EF-4 | REQ-EF-5 | `test_error_framework.cpp` |
| AUD-EF-5 | REQ-EF-6 | `test_error_framework.cpp` |
| AUD-EF-6 | REQ-EF-7 | `test_force_models.cpp` |
| AUD-EF-7 | REQ-EF-12 | `test_error_framework.cpp` |
| AUD-EF-8 | REQ-EF-8 | `test_code_consistency.cpp` |
| AUD-EF-9 | REQ-EF-9 | `test_error_framework.cpp` |
| AUD-EF-10 | every REQ-EF-N | `test_error_framework.cpp` (self) |
