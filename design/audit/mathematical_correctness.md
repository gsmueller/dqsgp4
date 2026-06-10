# Audit: Mathematical Correctness

Verifies that the implementations of `DualNumber<T>`, `Quaternion<T>`,
`DualQuaternion<T>`, and (later) `Pose<T>`, `Twist<T>`, `Wrench<T>`,
`Inertia<T>` satisfy the algebraic laws in
`design/specifications/dual_quaternion_algebra.md` and related specs.

Each item has a stable identifier (`AUD-MC-N`) and cites the requirement it
verifies. Enforced by `tests/test_dual_number.cpp`,
`tests/test_quaternion.cpp`, `tests/test_dual_quaternion.cpp`, and related
test files at CI time.

Throughout, "A == B" is read as `(A - B).max_total_error() ≤ tolerance`, where
tolerance is the sum of A's and B's reported `total_error()` plus a
double-precision-scaled slack term for the comparison itself. Tests fail if
the inequality fails — not if a fixed numerical threshold is exceeded.


## DualNumber audits

### AUD-MC-1 — ε² = 0

Verifies REQ-DQ-1.

Sample `DualNumber` instances. Assert `epsilon() * epsilon() == zero()`.


### AUD-MC-2 — Commutative ring axioms

Verifies REQ-DQ-2.

Sample triples of `DualNumber`. Assert each of the eight axiom identities in
REQ-DQ-2.


### AUD-MC-3 — Forward-mode differentiation

Verifies REQ-DQ-3.

Sample `DualNumber` a = (a.real, 1). For f in {sqrt, sin, cos}, assert

    f(a).real == f(a.real)
    f(a).dual == f'(a.real)

using closed-form derivatives.


## Quaternion audits

### AUD-MC-4 — Hamilton-product associativity

Verifies REQ-DQ-4.

Sample triples (q1, q2, q3) including unit and non-unit instances. Assert
(q1 q2) q3 == q1 (q2 q3).


### AUD-MC-5 — Identity laws

Verifies REQ-DQ-5.

Sample q. Assert q · identity() == q, identity() · q == q, q + zero() == q.


### AUD-MC-6 — Conjugate is involutive

Verifies REQ-DQ-6.

Sample q. Assert q.conjugate().conjugate() == q.


### AUD-MC-7 — Conjugate of product

Verifies REQ-DQ-7.

Sample (q1, q2). Assert (q1 q2)* == q2* q1*.


### AUD-MC-8 — Magnitude is multiplicative

Verifies REQ-DQ-8.

Sample (q1, q2). Assert |q1 q2| == |q1| |q2| and |q1 q2|² == |q1|² |q2|².


### AUD-MC-9 — Inverse identity

Verifies REQ-DQ-9.

Sample non-zero q. Assert q · q⁻¹ == identity() and q⁻¹ · q == identity().


### AUD-MC-10 — Rotation preserves length

Verifies REQ-DQ-10.

Sample unit q and arbitrary v. Assert |q.rotate(v)| == |v|.


### AUD-MC-11 — Rotation composition

Verifies REQ-DQ-11.

Sample unit (q1, q2) and arbitrary v. Assert (q1 q2).rotate(v) ==
q1.rotate(q2.rotate(v)).


### AUD-MC-12 — exp / log round-trip on the half-angle ball

Verifies REQ-DQ-12.

(a) Sample unit q with w.value ≥ 0. Assert exp_pure(q.log_unit()) == q.

(b) Sample v with |v| ≤ π/2. Assert log_unit(exp_pure(v)) == v.

(c) Sample v with |v| > π/2. Assert log_unit(exp_pure(v)) equals the
shortest-path equivalent of v (rescaled so |·| ≤ π/2).


### AUD-MC-13 — Axis-angle round-trip

Verifies REQ-DQ-13.

Sample unit axis and theta ∈ (0, π). Construct q via from_axis_angle. Assert
q.log_unit() == (theta/2) · axis.


### AUD-MC-14 — Normalize is idempotent on unit q

Verifies REQ-DQ-14.

Sample unit q. Assert q.normalized() == q. (Tests representational-error
sensitivity, not algebraic equivalence — slack term should account for one
sqrt and four divisions.)


## DualQuaternion audits

(Detailed once `src/math/dual_quaternion.h` lands.)

### AUD-MC-15 — Dual quaternion composition associativity

Verifies REQ-DQ-15.


### AUD-MC-16 — Conjugate dualities

Verifies REQ-DQ-16.


### AUD-MC-17 — Pose action composition

Verifies REQ-DQ-17.


### AUD-MC-18 — Screw exp / log round-trip

Verifies REQ-DQ-18.


## Sampling discipline

Each audit test uses a deterministic seeded PRNG (`std::mt19937` with a fixed
seed declared at the top of the test file). Samples include:

- 8 random unit quaternions covering low / mid / high angle, all octants of
  the unit ball.
- 8 random non-unit quaternions at various magnitudes (10⁻³ to 10³).
- Edge cases: identity, antipode of identity (q = −1), near-singular
  configurations (w very small, w very near 1).

Random sampling is for coverage breadth; edge cases are mandatory.


## Cross-reference

| AUD-MC-N | Verifies | Tested by |
|---|---|---|
| AUD-MC-1..3 | REQ-DQ-1..3 | `tests/test_dual_number.cpp` |
| AUD-MC-4..14 | REQ-DQ-4..14 | `tests/test_quaternion.cpp` |
| AUD-MC-15..18 | REQ-DQ-15..18 | `tests/test_dual_quaternion.cpp` |
