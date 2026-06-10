# Specification: Dual Quaternion Algebra

Mathematical laws that the implementations of `DualNumber<T>`, `Quaternion<T>`,
and `DualQuaternion<T>` must satisfy. Each law has a stable identifier
(`REQ-DQ-N`) and is verified by a test in `tests/test_dual_quaternion.cpp`
(naming convention: `test_REQ_DQ_N_<short_name>`).

These are correctness specs, not error-budget specs (those live in
`design/specifications/error_framework.md`). "A == B" in this document means the
identity holds to within the reported `total_error()` of both sides; tests
assert this inequality.

Every law here serves OBJ-1: the dual quaternion six-DOF state `(M̂, Ω̂)` is
only well-defined if the underlying dual-number, quaternion, and dual-quaternion
algebra is correct. The DualNumber laws underpin the dual part; the Quaternion
laws underpin attitude; the DualQuaternion laws underpin the full SE(3) pose.


## DualNumber

### REQ-DQ-1 — ε² = 0

Verifies: OBJ-1 — the defining property of dual numbers, which underpins every
closed-form formula in `dual_number.h`.

    epsilon() * epsilon() == zero()


### REQ-DQ-2 — Commutative ring axioms

Verifies: OBJ-1.

For sampled `DualNumber` triples (a, b, c):

    (a + b) + c == a + (b + c)            (additive associativity)
    a + b == b + a                        (additive commutativity)
    a + zero() == a                       (additive identity)
    a + (-a) == zero()                    (additive inverse)
    (a * b) * c == a * (b * c)            (multiplicative associativity)
    a * b == b * a                        (multiplicative commutativity)
    a * identity() == a                   (multiplicative identity)
    a * (b + c) == a*b + a*c              (distributivity)


### REQ-DQ-3 — Forward-mode automatic differentiation

Verifies: OBJ-1.

For sampled `DualNumber` of the form a = (a.real, 1) and any of
{`sqrt`, `sin`, `cos`}:

    f(a).real == f(a.real)
    f(a).dual == derivative_of_f(a.real)

(i.e., dual-part transports the first derivative through `f`).


## Quaternion

### REQ-DQ-4 — Hamilton-product associativity

Verifies: OBJ-1.

For sampled `Quaternion` triples (q1, q2, q3):

    (q1 * q2) * q3 == q1 * (q2 * q3)


### REQ-DQ-5 — Identity laws

Verifies: OBJ-1.

For any q:

    q * identity() == q
    identity() * q == q
    q + zero() == q


### REQ-DQ-6 — Conjugate is involutive

Verifies: OBJ-1.

For any q:

    q.conjugate().conjugate() == q


### REQ-DQ-7 — Conjugate of product

Verifies: OBJ-1.

For any (q1, q2):

    (q1 * q2).conjugate() == q2.conjugate() * q1.conjugate()


### REQ-DQ-8 — Magnitude is multiplicative

Verifies: OBJ-1.

For any (q1, q2):

    (q1 * q2).magnitude() == q1.magnitude() * q2.magnitude()
    (q1 * q2).magnitude_squared() == q1.magnitude_squared() * q2.magnitude_squared()


### REQ-DQ-9 — Inverse identity

Verifies: OBJ-1.

For non-zero q:

    q * q.inverse() == identity()
    q.inverse() * q == identity()


### REQ-DQ-10 — Unit quaternion rotation preserves vector length

Verifies: OBJ-1.

For unit q (|q| = 1) and any 3-vector v:

    q.rotate(v).magnitude() == v.magnitude()


### REQ-DQ-11 — Rotation composition

Verifies: OBJ-1.

For unit (q1, q2) and any v:

    (q1 * q2).rotate(v) == q1.rotate(q2.rotate(v))


### REQ-DQ-12 — exp_pure and log_unit are inverses on the half-angle ball

Verifies: OBJ-1.

For sampled unit q with w.value ≥ 0 (rotation angle ≤ π):

    exp_pure(q.log_unit()) == q

For sampled v with |v| ≤ π/2:

    log_unit(exp_pure(v)) == v

For |v| > π/2, the shortest-path convention in `log_unit` returns the
equivalent shorter rotation, not v itself; the round-trip identity is
*expected* to fail there, and tests assert the shortest-path identity instead:

    log_unit(exp_pure(v)) == shortest_path_equivalent(v)

where `shortest_path_equivalent(v)` is defined as `v` rescaled so |·| ≤ π/2.


### REQ-DQ-13 — Axis-angle round-trip

Verifies: OBJ-1.

For unit axis n̂ and angle theta ∈ (0, π):

    Quaternion::from_axis_angle(n_hat, theta).log_unit() == (theta/2) * n_hat


### REQ-DQ-14 — Normalize is idempotent on unit q

Verifies: OBJ-1.

For sampled unit q:

    q.normalized() == q   (to within representation error)


## DualQuaternion

(Detailed once `src/math/dual_quaternion.h` lands. Anticipated structure:)


### REQ-DQ-15 — Dual quaternion composition

Verifies: OBJ-1.

For sampled DualQuaternion (M1, M2, M3):

    (M1 * M2) * M3 == M1 * (M2 * M3)


### REQ-DQ-16 — Conjugate dualities

Verifies: OBJ-1.

The dual quaternion algebra admits three conjugates; each is involutive and
each has its own product rule. Specified at implementation time.


### REQ-DQ-17 — Pose composition on Vector3

Verifies: OBJ-1.

For unit DualQuaternion (M1, M2) and any v ∈ Vector3:

    (M1 * M2).apply(v) == M1.apply(M2.apply(v))


### REQ-DQ-18 — Screw exp / log are inverses on the half-screw ball

Verifies: OBJ-1.

For sampled unit dual quaternion M with rotation angle ≤ π:

    DualQuaternion::exp_screw(M.log_screw()) == M


## Cross-reference

| REQ-DQ-N | Type | Verified by |
|---|---|---|
| REQ-DQ-1..3 | DualNumber | `tests/test_dual_number.cpp` |
| REQ-DQ-4..14 | Quaternion | `tests/test_quaternion.cpp` |
| REQ-DQ-15..18 | DualQuaternion | `tests/test_dual_quaternion.cpp` |

Audited by `design/audit/mathematical_correctness.md`.
