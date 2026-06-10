# Theoretical Basis Audit — `src/math/dual_quaternion.h`

**File**: `src/math/dual_quaternion.h` (356 lines)
**Scope**: `DualQuaternion<T>` struct, 28 distinct formulas
**Theory anchors**: Study (1903) Geometrie der Dynamen; Clifford (1873) Preliminary Sketch of Biquaternions; Murray, Li & Sastry (1994) §2.5, §3.2; Selig (2005) Geometric Fundamentals of Robotics §3.4, §9.1–9.3; Kavan, Collins, Žára & O'Sullivan (2008) "Geometric Skinning with Approximate Dual Quaternion Blending" §5 (sclerp).
**Headline question** (per audit charter): on `exp_screw`, `log_screw`, `sclerp` — is the implementation a **screw exponential** (the Lie-group geodesic on SE(3)) or a linear interpolation in components? **Answer**: screw exponential. See cards 25, 27, 28.

---

## Formula 1: Default constructor

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::ctor_default
Location:               src/math/dual_quaternion.h:64
Mathematical statement: DualQuaternion() = (0 + 0 i + 0 j + 0 k) + ε (0 + 0 i + 0 j + 0 k)

THEORY
  Underlying theorem:   None — default initialization of both Quaternion<T> members.
  Primary reference:    C++ struct default construction.
  Domain of validity:   All T.

METHOD
  Method declared:      Delegate to Quaternion<T>::Quaternion() for each member.
  Method implemented:   `real(), dual()` — both default-constructed (each gives w=x=y=z=0).
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — no computation performed; error state inherited from
                        default-constructed TrackedValue<T> components.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    tests/test_dual_quaternion/ — default-construct, verify all
                        eight TrackedValue<T> components equal zero.

NOTES
  - No arithmetic, no error budget impact.
```

---

## Formula 2: Constructor from real & dual parts

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::ctor_bipart
Location:               src/math/dual_quaternion.h:66-67
Mathematical statement: DualQuaternion(r, d) = r + ε d  for Quaternion<T> r, d.

THEORY
  Underlying theorem:   Definition of M̂ ∈ H ⊕ εH (Study 1903 §2; Selig 2005 §9.1).
  Primary reference:    Study (1903) Geometrie der Dynamen, Ch. 1.
  Domain of validity:   Any pair of quaternions (no unit-norm constraint at this
                        layer; the SE(3)-unit constraint is enforced separately).

METHOD
  Method declared:      Direct initialization of `real` and `dual` members.
  Method implemented:   `real(r), dual(d)` — direct member-init. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — no arithmetic. Error state inherited from inputs.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    Constructor round-trip — components match inputs.

NOTES
  - Used as the algebraic constructor throughout the file.
```

---

## Formula 3: Factory `zero()`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::zero
Location:               src/math/dual_quaternion.h:72
Mathematical statement: zero() returns the additive identity 0̂ = 0 + ε·0.

THEORY
  Underlying theorem:   Additive identity in the dual-quaternion ring H[ε]/(ε²).
  Primary reference:    Study (1903); Selig (2005) §9.1.
  Domain of validity:   All T.

METHOD
  Method declared:      Return default-constructed DualQuaternion.
  Method implemented:   `return DualQuaternion();` — trivial.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — no computation.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       AUD-MC-1 (additive identity in the dual-quaternion ring)
  Verification test:    M̂ + zero() == M̂ for arbitrary M̂.

NOTES
  - Trivial; ring axiom.
```

---

## Formula 4: Factory `identity()`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::identity
Location:               src/math/dual_quaternion.h:75-77
Mathematical statement: identity() returns the multiplicative identity
                        1̂ = (1 + 0i + 0j + 0k) + ε · 0.

THEORY
  Underlying theorem:   Multiplicative identity in the dual-quaternion ring;
                        also the identity element of SE(3) under the unit-DQ
                        embedding (Murray-Li-Sastry §3.2; Selig §9.2).
  Primary reference:    Study (1903); Murray, Li & Sastry (1994) §2.5.
  Domain of validity:   All T.

METHOD
  Method declared:      Compose Quaternion<T>::identity() (real) with
                        Quaternion<T>::zero() (dual).
  Method implemented:   `DualQuaternion(Quaternion<T>::identity(),
                        Quaternion<T>::zero())` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — no arithmetic; both inputs are exact constants.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       AUD-MC-2 (multiplicative identity in the dual-quaternion ring)
  Verification test:    M̂ · identity() == M̂ and identity() · M̂ == M̂.

NOTES
  - Encodes the SE(3) identity (no rotation, no translation).
```

---

## Formula 5: `from_pose(q_r, t)`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::from_pose
Location:               src/math/dual_quaternion.h:86-91
Mathematical statement: M̂ = q_r + ε · (1/2) · t_pure · q_r,
                        where t_pure is the pure quaternion (0, t).

THEORY
  Underlying theorem:   Standard pose-to-DQ embedding of SE(3): given unit q_r ∈ S³
                        (a rotation) and body-frame translation t ∈ ℝ³, the unit
                        dual quaternion encoding the rigid motion is
                          M̂ = q_r + ε · (1/2) · t · q_r          (1)
                        (with t a pure quaternion). Equivalently
                        q_d = (1/2) t q_r, giving t = 2 q_d q_r* on inversion
                        (cf. Formula 7 below).
  Primary reference:    Selig (2005) §9.3 "Rigid-body motions and dual quaternions",
                        Eq (9.10); Murray-Li-Sastry §3.2 (exponential coordinates,
                        which yield the same expression after specialization).
  Domain of validity:   Caller-supplied q_r unit (no internal renormalization);
                        t arbitrary in ℝ³.

METHOD
  Method declared:      Closed-form algebraic identity (1); no series.
  Method implemented:   Lines 88-90:
                          t_pure = Quaternion<T>::pure(t)
                          q_d = (t_pure * q_r) * ratio<T>(1, 2)
                          return DualQuaternion(q_r, q_d)
                        Matches (1). ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        The dual part is the Hamilton product t_pure · q_r times
                        the exact rational 1/2; errors propagate through the
                        TrackedValue<T> operators in each component.
  Bound implemented:    Delegated to Quaternion<T>::operator*() and scalar mul. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-15 (unit-DQ constraint q_r·q_d^* + q_d·q_r^* = 0;
                        this construction satisfies it identically since
                        2·Re(q_d q_r*) = Re(t_pure · q_r · q_r*) = Re(t_pure) · |q_r|² = 0)
  Verification test:    tests/test_dual_quaternion/ — from_pose → translation() round-trip;
                        magnitude_squared() == DualNumber(|q_r|², 0).

NOTES
  - Caller must ensure |q_r| = 1; otherwise the magnitude has a non-zero dual part
    (Formula 20) and the SE(3) embedding fails (REQ-DQ-15).
  - This is one of two named constructors; pairs with `from_screw` (Formula 6).
```

---

## Formula 6: `from_screw(omega, v)`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::from_screw
Location:               src/math/dual_quaternion.h:97-101
Mathematical statement: Ŝ = (0, ω) + ε (0, v) — a pure dual quaternion encoding
                        a body-frame screw with angular part ω and linear part v.

THEORY
  Underlying theorem:   Pure dual quaternions span the Lie algebra se(3) of SE(3)
                        (Selig 2005 §9.4; Murray-Li-Sastry §2.5 Theorem 2.4). A
                        screw with angular velocity ω and linear velocity v at the
                        origin is the dual vector ω + ε v, lifted to a pure DQ.
  Primary reference:    Selig (2005) §9.4 "Screws and dual quaternions";
                        Murray, Li & Sastry (1994) §2.5 "Twists".
  Domain of validity:   Any ω, v ∈ ℝ³.

METHOD
  Method declared:      Lift each 3-vector to a pure quaternion (zero scalar part).
  Method implemented:   `DualQuaternion(Quaternion<T>::pure(omega),
                        Quaternion<T>::pure(v))` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (no arithmetic; component copies only)
  Bound formula:        Component error state inherited from inputs.
  Bound implemented:    Delegated to Quaternion<T>::pure(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       N/A (no computation)
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a (this is a constructor; algebraic checks live in
                        twist.h / wrench.h tests)
  Verification test:    tests/test_dual_quaternion/ — angular() and linear() round-trip.

NOTES
  - Used by `exp_screw(twist)` overload (Formula 26) and dynamic propagators
    (twist.h, wrench.h).
```

---

## Formula 7: `rotation()` accessor

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::rotation
Location:               src/math/dual_quaternion.h:106
Mathematical statement: q_r ← real (the rotation quaternion of a pose).

THEORY
  Underlying theorem:   Definition — the real part of a unit pose-DQ is the SO(3)
                        rotation quaternion (Selig §9.3).
  Primary reference:    Selig (2005) §9.3.
  Domain of validity:   Meaningful when M̂ is a pose (unit DQ).

METHOD
  Method declared:      Return `real` member by value.
  Method implemented:   `return real;` — trivial.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — no arithmetic; copy of input.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    from_pose(q, t).rotation() == q.

NOTES
  - Pure accessor; no algebraic content.
```

---

## Formula 8: `translation()` accessor

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::translation
Location:               src/math/dual_quaternion.h:111-116
Mathematical statement: t = 2 · (q_d · q_r*).vector()  for a unit pose.

THEORY
  Underlying theorem:   Inversion of from_pose (Formula 5). Starting from
                          q_d = (1/2) t_pure · q_r,
                        right-multiply both sides by q_r*:
                          q_d · q_r* = (1/2) t_pure · q_r · q_r* = (1/2) t_pure
                        (since |q_r|² = 1 for unit q_r). Hence
                          t = 2 · vector_part(q_d · q_r*).
  Primary reference:    Selig (2005) §9.3 Eq (9.11); Kavan et al. (2008) §3.
  Domain of validity:   Unit M̂. For non-unit M̂ the scalar part of q_d · q_r*
                        is non-zero and the result is meaningless (file comment).

METHOD
  Method declared:      Closed-form: form q_d · q_r*, take vector part, scale by 2.
  Method implemented:   Lines 112-115:
                          td_over_2 = dual * real.conjugate()
                          v = td_over_2.vector()
                          two = exact<T>(2)
                          return Vector3<T>(two*v.x, two*v.y, two*v.z)
                        Matches the derivation. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Hamilton product q_d · q_r* propagates errors per
                        Quaternion<T>::operator*; scaling by 2 is exact.
  Bound implemented:    Delegated to component operators. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-16 (point action; consistent with apply() at origin)
  Verification test:    tests/test_dual_quaternion/ — from_pose(q, t).translation() ≈ t.

NOTES
  - The factor of 2 is the Study factor in the DQ pose embedding (1/2 in encoding,
    2 in decoding).
```

---

## Formula 9: `angular()` accessor

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::angular
Location:               src/math/dual_quaternion.h:119
Mathematical statement: ω = real.vector()  for a pure (twist/wrench) DQ.

THEORY
  Underlying theorem:   Definition — the vector part of the real component of a
                        pure DQ is the angular component of the screw
                        (Murray-Li-Sastry §2.5; Selig §9.4).
  Primary reference:    Selig (2005) §9.4.
  Domain of validity:   Meaningful when M̂ is pure (i.e. real.scalar() = 0).

METHOD
  Method declared:      Return real.vector() — extract (x, y, z) from the real
                        quaternion.
  Method implemented:   `return real.vector();` — trivial.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — component copy.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    from_screw(ω, v).angular() == ω.

NOTES
  - Pure accessor. Used by `exp_screw(twist)` overload (Formula 26).
```

---

## Formula 10: `linear()` accessor

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::DualQuaternion::linear
Location:               src/math/dual_quaternion.h:122
Mathematical statement: v = dual.vector()  for a pure (twist/wrench) DQ.

THEORY
  Underlying theorem:   Definition — the vector part of the dual component of a
                        pure DQ is the linear component of the screw
                        (Murray-Li-Sastry §2.5; Selig §9.4).
  Primary reference:    Selig (2005) §9.4.
  Domain of validity:   Meaningful when M̂ is pure.

METHOD
  Method declared:      Return dual.vector().
  Method implemented:   `return dual.vector();` — trivial.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A — component copy.
  Bound implemented:    N/A.
  Bound verdict:        ✓ n/a.

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    from_screw(ω, v).linear() == v.

NOTES
  - Pure accessor; pairs with `angular()` (Formula 9).
```

---

## Formula 11: Component-wise addition

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::operator_add
Location:               src/math/dual_quaternion.h:127-130
Mathematical statement: (r₁ + εd₁) + (r₂ + εd₂) = (r₁ + r₂) + ε(d₁ + d₂).

THEORY
  Underlying theorem:   Ring axiom: component-wise addition in H ⊕ εH ≅ H[ε]/(ε²).
  Primary reference:    Study (1903); Clifford (1873); Selig (2005) §9.1.
  Domain of validity:   All DualQuaternion<T> inputs.

METHOD
  Method declared:      Component-wise addition; delegate to Quaternion<T> + .
  Method implemented:   `DualQuaternion(a.real + b.real, a.dual + b.dual)` —
                        correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Componentwise quaternion addition; each TrackedValue<T>
                        addition obeys the per-category triangle inequality
                        (REQ-EF-3).
  Bound implemented:    Delegated to Quaternion<T>::operator+. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-3 (additive group on H ⊕ εH: associative,
                        commutative, identity, inverse)
  Verification test:    tests/test_dual_quaternion/ — (A+B)+C == A+(B+C),
                        A+B == B+A, A+zero == A.

NOTES
  - Exact dual-algebra operation; no truncation introduced.
```

---

## Formula 12: Component-wise subtraction

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::operator_sub
Location:               src/math/dual_quaternion.h:133-136
Mathematical statement: (r₁ + εd₁) − (r₂ + εd₂) = (r₁ − r₂) + ε(d₁ − d₂).

THEORY
  Underlying theorem:   Ring axiom — additive inverse + componentwise sum.
  Primary reference:    Study (1903); Selig (2005) §9.1.
  Domain of validity:   All DualQuaternion<T> inputs.

METHOD
  Method declared:      Component-wise subtraction.
  Method implemented:   `DualQuaternion(a.real - b.real, a.dual - b.dual)` —
                        correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Componentwise via TrackedValue<T> subtraction.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-3
  Verification test:    A − A == zero; (A+B) − B == A.

NOTES
  - Exact dual-algebra operation.
```

---

## Formula 13: Unary negation

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::operator_neg
Location:               src/math/dual_quaternion.h:139-141
Mathematical statement: −(r + εd) = (−r) + ε(−d).

THEORY
  Underlying theorem:   Ring axiom: additive inverse exists componentwise.
  Primary reference:    Study (1903); Selig (2005) §9.1.
  Domain of validity:   All DualQuaternion<T>.

METHOD
  Method declared:      Component-wise negation.
  Method implemented:   `DualQuaternion(-real, -dual)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Negation preserves the per-category error magnitudes.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-3
  Verification test:    −(−A) == A; A + (−A) == zero.

NOTES
  - Exact.
```

---

## Formula 14: Dual-quaternion product

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::operator_mul
Location:               src/math/dual_quaternion.h:148-154
Mathematical statement: (r₁ + εd₁)(r₂ + εd₂) = r₁r₂ + ε(r₁d₂ + d₁r₂),
                        since ε² = 0.

THEORY
  Underlying theorem:   Algebra structure on H[ε]/(ε²): the full expansion
                        (r₁ + εd₁)(r₂ + εd₂) = r₁r₂ + ε(r₁d₂ + d₁r₂)
                                              + ε²(d₁d₂)
                        is truncated by the defining relation ε² = 0. Each pair
                        of quaternions is multiplied with Hamilton's product
                        (Quaternion<T>::operator*), which is non-commutative.
  Primary reference:    Study (1903) Ch. 1; Clifford (1873); Selig (2005) §9.1
                        Eq (9.4); Murray, Li & Sastry (1994) §2.5.
  Domain of validity:   All DualQuaternion<T>.

METHOD
  Method declared:      Closed-form DQ multiplication using ε² = 0.
  Method implemented:   Lines 150-153:
                          real: a.real * b.real
                          dual: a.real * b.dual + a.dual * b.real
                        Matches (r₁r₂, r₁d₂ + d₁r₂). ✓
  Match verdict:        ✓ matched. The ε² = 0 truncation is exact (definitional),
                        not approximate.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Composed from four Hamilton products and one
                        quaternion addition; each propagates per REQ-EF-3.
                        No algebra-truncation term is needed because ε² = 0
                        is exact.
  Bound implemented:    Delegated to Quaternion<T>::operator* and +. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-4 (associativity (AB)C == A(BC)),
                        AUD-MC-5 (left/right distributivity),
                        AUD-MC-6 (identity 1̂ · A == A · 1̂ == A)
  Verification test:    tests/test_dual_quaternion/ —
                        - associativity over random DQ triples
                        - identity neutrality on left and right
                        - non-commutativity (typically AB ≠ BA)
                        - SE(3) composition matches 4x4 homogeneous product
                          on independent pose representations

NOTES
  - The truncation ε² = 0 is the defining relation of the dual ring; not a
    method-theory mismatch.
  - Non-commutative on both factors (rotations don't commute).
```

---

## Formula 15: Scalar–DQ multiplication (left)

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::operator_mul_scalar_left
Location:               src/math/dual_quaternion.h:157-160
Mathematical statement: s(r + εd) = (sr) + ε(sd).

THEORY
  Underlying theorem:   Module axiom: scalar action distributes over real and
                        dual parts of H[ε]/(ε²).
  Primary reference:    Standard module theory; Study (1903).
  Domain of validity:   All s ∈ TrackedValue<T>, M̂ ∈ DualQuaternion<T>.

METHOD
  Method declared:      Scale real and dual components by s.
  Method implemented:   `DualQuaternion(s * m.real, s * m.dual)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Each component multiplication propagates per REQ-EF-3.
  Bound implemented:    Delegated to Quaternion<T>::operator*(scalar, q). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-7 (module/vector-space scalar law)
  Verification test:    s(A+B) == sA + sB; (s+t)A == sA + tA; 1·A == A.

NOTES
  - Exact module operation.
```

---

## Formula 16: Scalar–DQ multiplication (right)

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::operator_mul_scalar_right
Location:               src/math/dual_quaternion.h:163-166
Mathematical statement: (r + εd)s = (rs) + ε(ds).

THEORY
  Underlying theorem:   Module axiom (symmetric to Formula 15).
  Primary reference:    Standard module theory.
  Domain of validity:   All s, M̂.

METHOD
  Method declared:      Scale real and dual by s on the right.
  Method implemented:   `DualQuaternion(m.real * s, m.dual * s)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Componentwise via TrackedValue<T>.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-7
  Verification test:    A·s == s·A (commutative because TrackedValue<T> scalars
                        commute with quaternion components).

NOTES
  - Exact module operation.
```

---

## Formula 17: Quaternion-style conjugate `M̂*`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::conjugate
Location:               src/math/dual_quaternion.h:173-175
Mathematical statement: M̂* = q_r* + ε q_d*  (REQ-DQ-16).

THEORY
  Underlying theorem:   First of the three DQ conjugates. Equals the multiplicative
                        inverse for a unit DQ, since
                          (q_r + εq_d)(q_r* + εq_d*)
                            = |q_r|² + ε(q_r q_d* + q_d q_r*)
                            = 1 + 0 = 1̂
                        when |q_r| = 1 and q_r q_d* + q_d q_r* = 0 (the SE(3)
                        unit-DQ constraint, REQ-DQ-15).
  Primary reference:    Selig (2005) §9.1, §9.3; Murray, Li & Sastry (1994) §2.5
                        (Theorem 2.5).
  Domain of validity:   All DualQuaternion<T>; M̂ M̂* = 1̂ only for unit pose.

METHOD
  Method declared:      Closed-form: conjugate each quaternion factor.
  Method implemented:   `DualQuaternion(real.conjugate(), dual.conjugate())` —
                        correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Conjugation is exact in each component (sign flips on
                        x,y,z); no new error.
  Bound implemented:    Delegated to Quaternion<T>::conjugate(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-8 (conjugate is anti-involution: (AB)* = B*A*),
                        AUD-MC-15 (unit-pose inverse via conjugate)
  Verification test:    tests/test_dual_quaternion/ — (AB)* == B*A*;
                        from_pose(q, t).conjugate() · from_pose(q, t) == identity().

NOTES
  - One of three conjugates documented in the file header (lines 23-28).
```

---

## Formula 18: Dual-number conjugate `M̂_ε`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::dual_conjugate
Location:               src/math/dual_quaternion.h:180-182
Mathematical statement: M̂_ε = q_r − ε q_d.

THEORY
  Underlying theorem:   Second of the three DQ conjugates — flips ε → −ε at the
                        dual-ring level (no quaternion conjugation).
  Primary reference:    Selig (2005) §9.1 (three conjugates).
  Domain of validity:   All DualQuaternion<T>.

METHOD
  Method declared:      Negate the dual part; leave the real part alone.
  Method implemented:   `DualQuaternion(real, -dual)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Negation preserves per-category magnitudes.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a (involution: (M̂_ε)_ε = M̂; tested as part of
                        Formula 19 combined conjugate consistency)
  Verification test:    tests/test_dual_quaternion/ — double dual-conjugate is
                        identity.

NOTES
  - Distinct from the quaternion conjugate (Formula 17); they commute, and the
    product of the two is the combined conjugate (Formula 19).
```

---

## Formula 19: Combined conjugate `M̂♯` (REQ-DQ-16)

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::combined_conjugate
Location:               src/math/dual_quaternion.h:188-190
Mathematical statement: M̂♯ = q_r* − ε q_d*  (the composition of `conjugate()`
                        and `dual_conjugate()` in either order).

THEORY
  Underlying theorem:   Third DQ conjugate. Used in the point action
                          P̂' = M̂ P̂ M̂♯           (REQ-DQ-17)
                        which gives p' = Rp + t for a unit pose M̂ acting on the
                        lifted point P̂ = 1 + ε p_pure (cf. Formula 23).
                        Equivalently M̂♯ = (M̂*)_ε = (M̂_ε)*.
  Primary reference:    Selig (2005) §9.3, Eq (9.13); Kavan et al. (2008) §3.
  Domain of validity:   All DualQuaternion<T>.

METHOD
  Method declared:      Closed-form: conjugate the real part, conjugate and negate
                        the dual part.
  Method implemented:   `DualQuaternion(real.conjugate(), -dual.conjugate())` —
                        correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Two component conjugations and one negation, each exact;
                        errors inherited per REQ-EF-3.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-16 (consistency of the point action M̂ P̂ M̂♯,
                        cross-checked at Formula 23)
  Verification test:    tests/test_dual_quaternion/ — combined_conjugate() ==
                        conjugate().dual_conjugate() == dual_conjugate().conjugate();
                        point action gives Rp + t.

NOTES
  - Choice of which conjugate is "the" inverse depends on context: M̂* inverts unit
    poses (Formula 17), but M̂♯ is the partner in the M̂ P̂ M̂♯ point action.
```

---

## Formula 20: Dual-valued squared magnitude

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::magnitude_squared
Location:               src/math/dual_quaternion.h:197-202
Mathematical statement: ||M̂||²_dual = |q_r|² + ε · 2 · <q_r, q_d>_ℝ⁴
                        where <q_r, q_d> = q_r.w·q_d.w + q_r.x·q_d.x +
                                            q_r.y·q_d.y + q_r.z·q_d.z.
                        Equivalent form: ||M̂||²_dual = M̂ M̂* with the result
                        returned as a DualNumber<T>.

THEORY
  Underlying theorem:   The DQ norm is a dual-number-valued quantity. Direct
                        expansion of M̂ M̂*:
                          (q_r + εq_d)(q_r* + εq_d*)
                            = q_r q_r* + ε(q_r q_d* + q_d q_r*) + ε²(q_d q_d*)
                            = |q_r|² + ε · 2 Re(q_r q_d*) + 0
                        (since ε² = 0 and Re(q_r q_d*) = q_r·q_d as a 4-vector
                        inner product). For a unit pose, |q_r|² = 1 and the
                        4-vector inner product vanishes (REQ-DQ-15).
  Primary reference:    Selig (2005) §9.1 Eq (9.6); Murray, Li & Sastry (1994)
                        §2.5 Eq (2.39).
  Domain of validity:   All DualQuaternion<T>.

METHOD
  Method declared:      Closed-form: compute the two components above.
  Method implemented:   Lines 198-201:
                          real_sq = real.magnitude_squared()
                          qr_qd_dot = real.w*dual.w + real.x*dual.x +
                                      real.y*dual.y + real.z*dual.z
                          return DualNumber<T>(real_sq, exact<T>(2) * qr_qd_dot)
                        Matches the derivation. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Real part: TrackedValue<T> squared-magnitude bound
                        (Pythagorean, four-term sum of squares). Dual part:
                        four-term inner product times exact 2.
  Bound implemented:    Delegated to TrackedValue<T> ops. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-15 (||M̂||²_dual == 1̂ characterizes unit pose);
                        AUD-MC-17 (||AB||²_dual == ||A||²_dual · ||B||²_dual
                                   multiplicativity of the dual norm)
  Verification test:    tests/test_dual_quaternion/ — from_pose(q, t).magnitude_squared()
                        == DualNumber(|q|², 0); multiplicativity over random pairs.

NOTES
  - Result is DualNumber<T>, not TrackedValue<T>; downstream consumers must
    extract real/dual parts as appropriate.
```

---

## Formula 21: `inverse()` for unit DQ

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::inverse
Location:               src/math/dual_quaternion.h:211-213
Mathematical statement: M̂⁻¹ = M̂*  when |M̂|_dual = 1̂ (i.e. M̂ is a unit pose).

THEORY
  Underlying theorem:   For a unit pose, M̂ M̂* = 1̂ (cf. Formula 20 derivation),
                        so M̂* is the two-sided multiplicative inverse.
                        General non-unit DQ inverse requires dividing by the
                        dual-valued ||M̂||²_dual; that case is NOT handled by
                        this method (documented precondition).
  Primary reference:    Selig (2005) §9.3; Murray-Li-Sastry §2.5 Theorem 2.5.
  Domain of validity:   Unit pose only (preconditioned). For non-unit M̂ the
                        returned value is the quaternion conjugate, NOT the
                        inverse.

METHOD
  Method declared:      Closed-form (return the conjugate).
  Method implemented:   `return conjugate();` — correct under precondition.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Same as Formula 17 (no new error term).
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-14 (M̂ · M̂⁻¹ == 1̂ for unit pose)
  Verification test:    tests/test_dual_quaternion/ — from_pose(q, t).inverse() ·
                        from_pose(q, t) == identity() (to within REQ-EF-3 bound).

NOTES
  - For non-unit DQ the caller must instead compute M̂*/||M̂||²_dual using
    dual-number division (DualNumber<T>::operator/, Formula 13 in dual_number.md).
  - Used by `sclerp` (Formula 28) on unit-pose inputs.
```

---

## Formula 22: `normalized()` — SE(3) retraction

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::normalized
Location:               src/math/dual_quaternion.h:219-226
Mathematical statement: For an approximately-unit M̂ = (q_r, q_d) drifting off
                        SE(3), return (q̂_r, q̂_d) where
                          q̂_r = q_r / |q_r|             (rotation: project to S³)
                          q̂_d = q_d − <q̂_r, q_d>_ℝ⁴ · q̂_r
                                                          (dual: orthogonalize)
                        so that |q̂_r| = 1 and the SE(3) constraint
                        <q̂_r, q̂_d>_ℝ⁴ = 0 holds.

THEORY
  Underlying theorem:   Lie-group retraction back to SE(3) after a numerical
                        integration step (REQ-EF-15). The unit-pose manifold is
                        {(q_r, q_d) ∈ H × H : |q_r| = 1 and <q_r, q_d> = 0}; the
                        retraction is the orthogonal projection onto this
                        submanifold (Munthe-Kaas 1999 for the principle;
                        Kavan et al. 2008 §5.2 for the explicit DQ form).
  Primary reference:    Munthe-Kaas (1999) "High order Runge-Kutta methods on
                        manifolds"; Selig (2005) §9.3 (unit-DQ constraint);
                        Kavan et al. (2008) §5.2.
  Domain of validity:   q_r close to unit (|q_r| ≠ 0). Division by |q_r| is
                        safe for poses near the manifold.

METHOD
  Method declared:      Closed-form: normalize the real part, then subtract the
                        component of the dual along the (normalized) real to
                        enforce orthogonality.
  Method implemented:   Lines 220-225:
                          mr = real.magnitude()
                          qr_n = (real.w/mr, real.x/mr, real.y/mr, real.z/mr)
                          dual_dot = qr_n.w*dual.w + qr_n.x*dual.x +
                                     qr_n.y*dual.y + qr_n.z*dual.z
                          qd_n = dual − dual_dot * qr_n
                          return DualQuaternion(qr_n, qd_n)
                        This is the Gram–Schmidt step on the 4-vector pair
                        (q_r, q_d), normalizing q_r and orthogonalizing q_d
                        against q̂_r. ✓
  Match verdict:        ✓ matched — projection / Gram–Schmidt, not exp/log
                        retraction.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Errors propagate through TrackedValue<T> sqrt
                        (for magnitude), division, 4-term inner product, and
                        scalar–quaternion subtraction. No truncation introduced.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-15 (manifold retraction after step)
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-15 (post-retraction |M̂|_dual = 1̂)
  Verification test:    tests/test_dual_quaternion/ —
                        - normalized() output has unit real magnitude
                        - normalized() output satisfies <q̂_r, q̂_d> = 0
                        - on a unit input, normalized() is approximately the
                          identity (within REQ-EF-3 bound)

NOTES
  - This is the **first-order retraction** (projection); the exact exponential
    retraction would be exp(log(M̂)). Both serve the same purpose for small
    drift, and projection is cheaper.
  - The orthogonality step <q̂_r, q̂_d>=0 is the dual-part Gram–Schmidt; it
    eliminates the "twisted" component along the rotation, leaving a clean
    translation encoding.
```

---

## Formula 23: `apply(p)` — SE(3) point action

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::apply
Location:               src/math/dual_quaternion.h:240-247
Mathematical statement: p' = R(p) + t  via the algebraic identity
                          p' = vector_part(q_r p q_r* + (q_d q_r* − q_r q_d*))
                        derived from the lifted point P̂ = 1 + ε p_pure and
                        the sandwich M̂ P̂ M̂♯ (REQ-DQ-17). The two dual-part
                        terms combine to yield R(p) + t.

THEORY
  Underlying theorem:   SE(3) point action via the combined-conjugate sandwich
                        (Selig 2005 §9.3 Eq (9.13–9.15); Kavan et al. 2008 §3).
                        Expanding M̂ P̂ M̂♯ for P̂ = 1 + ε p_pure with M̂ unit:
                          M̂ P̂ M̂♯ = (q_r + εq_d)(1 + εp_pure)(q_r* − εq_d*)
                                  = (q_r + ε(q_d + q_r p_pure))(q_r* − εq_d*)
                                  = q_r q_r* + ε(q_d q_r* + q_r p_pure q_r*
                                                  − q_r q_d*)
                                  = 1 + ε(q_r p_pure q_r* + q_d q_r* − q_r q_d*)
                        The vector part of the dual coefficient is q_r p q_r* +
                        2 · vec(q_d q_r*) = R(p) + t (since q_d q_r* −
                        q_r q_d* = 2i · vec(q_d q_r*) on a pure-quaternion
                        difference and 2 · vec(q_d q_r*) = t by Formula 8).
  Primary reference:    Selig (2005) §9.3; Murray, Li & Sastry (1994) §2.5.
  Domain of validity:   Unit pose M̂ only (precondition).

METHOD
  Method declared:      Closed-form sandwich product computed in three Hamilton
                        products and one quaternion subtraction.
  Method implemented:   Lines 240-247:
                          p_pure = Quaternion<T>::pure(p)
                          qr_conj = real.conjugate()
                          qd_conj = dual.conjugate()
                          dual_out = real * p_pure * qr_conj
                                     + (dual * qr_conj − real * qd_conj)
                          return dual_out.vector()
                        Matches the derivation. ✓
  Match verdict:        ✓ matched — this is the **algebraic SE(3) action**, not
                        a 4×4 matrix multiply, not a Lie-bracket reformulation.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Three Hamilton products + one quaternion subtraction +
                        vector extraction; each component-level bound is
                        REQ-EF-3.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-DQ-17
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-16 (point action matches Rp + t from independent
                        rotation/translation extraction);
                        AUD-MC-18 (apply ∘ from_pose round-trip)
  Verification test:    tests/test_dual_quaternion/ —
                        - apply(p) matches q_r.rotate(p) + t on random unit poses
                        - apply(p) on identity == p
                        - composition: (M̂₁ · M̂₂).apply(p) == M̂₁.apply(M̂₂.apply(p))

NOTES
  - The "(q_d q_r* − q_r q_d*)" term contributes the translation; the
    "q_r p q_r*" term contributes the rotation.
  - For non-unit M̂ the result is meaningless (documented precondition); norm
    drift in time integration is corrected by `normalized()` (Formula 22).
```

---

## Formula 24: `apply_direction(v)`

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::apply_direction
Location:               src/math/dual_quaternion.h:252-254
Mathematical statement: v' = R(v) = q_r v q_r*  (rotation only, no translation).

THEORY
  Underlying theorem:   For a direction (free vector) there is no translation
                        component; the SE(3) action reduces to the SO(3) action
                        of the rotation quaternion.
  Primary reference:    Selig (2005) §9.3 (direction vs. point vectors);
                        Murray, Li & Sastry (1994) §2.4.
  Domain of validity:   Unit pose M̂ (so q_r is a unit rotation); v ∈ ℝ³.

METHOD
  Method declared:      Delegate to Quaternion<T>::rotate.
  Method implemented:   `return real.rotate(v);` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per REQ-EF-3)
  Bound formula:        Same as Quaternion<T>::rotate (v + 2 w (q_v × v) +
                        2 q_v × (q_v × v) per quaternion.h:178–184).
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-2
  AUD-MC applies:       AUD-MC-9 (rotation algebra closure / SO(3) action)
  Verification test:    tests/test_dual_quaternion/ — apply_direction(v) ==
                        rotation().rotate(v); composition consistency.

NOTES
  - Distinct from `apply(p)` (Formula 23) in that translation is not added.
  - Used by twist.h to rotate a body-frame angular velocity into world frame.
```

---

## Formula 25: `exp_screw(u, v)` — screw exponential

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::exp_screw_uv
Location:               src/math/dual_quaternion.h:269-292
Mathematical statement: For pure DQ Û = u + ε v (u, v lifted to pure
                        quaternions), the **screw exponential** in component form
                        (REQ-DQ-18):
                          q_r = cos|u| + sinc|u| · u    [= Quaternion::exp_pure(u)]
                          q_d.scalar = −(u·v) sinc|u|
                          q_d.vector = sinc|u| · v + (u·v) · β · u
                        with β = (cos|u| − sinc|u|) / |u|².

THEORY
  Underlying theorem:   The Lie-algebra-to-Lie-group exponential
                          exp : se(3) → SE(3),
                          Û ↦ Σ_{k≥0} Û^k / k!
                        closes to a finite formula on the unit-DQ representation
                        because (q_r + εq_d)² has all the structure needed and
                        ε² = 0 truncates the series exactly. Writing Û = u + εv
                        and using the dual-angle θ̂ = θ + ε d with θ = |u|,
                        d = (u·v)/θ:
                          exp(Û) = cos(θ̂) + sinc(θ̂) · Û.
                        Expanding cos(θ̂) = cos θ − ε d sin θ and
                        sinc(θ̂) = sinc θ + ε d · (cos θ − sinc θ)/θ²:
                          q_r = cos θ + sinc θ · u
                          q_d = −d sin θ + (sinc θ) · v
                                + d · (cos θ − sinc θ)/θ² · u  (vector part)
                          With d θ = u·v, the q_d vector becomes
                          (sinc θ) v + (u·v) β u, with β = (cos θ − sinc θ)/θ².
                          The q_d scalar is −d sin θ = −(u·v) sinc θ.
  Primary reference:    Murray, Li & Sastry (1994) §3.2 "Exponential coordinates
                        for rigid-body motion", Theorem 2.8;
                        Selig (2005) §9.4 "Screws and the exponential map";
                        Brockett (1984) "Robotic manipulators and the product of
                        exponentials formula".
                        Small-angle Taylor branches: see helpers in
                        `small_angle_series.h` (PASS per audit doc §5).
  Domain of validity:   All u, v ∈ ℝ³. The exponential is surjective onto
                        the connected component of SE(3); the inverse `log_screw`
                        (Formula 27) is single-valued on the half-angle ball
                        |u| ≤ π/2.

METHOD
  Method declared:      Screw exponential via the closed-form REQ-DQ-18, using
                        Taylor branches in `small_angle_series.h` for
                        `taylor_sinc(θ)` and `taylor_cos_minus_sinc_over_theta_sq(θ)`
                        when |θ| < 1e-4 to avoid 0/0.
  Method implemented:   Lines 270-291:
                          theta_sq = u.x² + u.y² + u.z²
                          theta = sqrt(theta_sq)
                          uv_dot = u·v
                          c = cos(theta)
                          sinc_t = taylor_sinc(theta, theta_sq)
                          beta = taylor_cos_minus_sinc_over_theta_sq(theta, theta_sq)
                          q_r = Quaternion(c, sinc_t·u.x, sinc_t·u.y, sinc_t·u.z)
                          qd_w = −(uv_dot · sinc_t)
                          uv_beta = uv_dot · beta
                          q_d = Quaternion(qd_w,
                                           sinc_t·v.x + uv_beta·u.x,
                                           sinc_t·v.y + uv_beta·u.y,
                                           sinc_t·v.z + uv_beta·u.z)
                          return DualQuaternion(q_r, q_d)
                        Matches the derivation. ✓
  Match verdict:        ✓ matched — this is the **screw exponential** (Lie-group
                        exp on SE(3)), NOT a linear blend, NOT a Padé approximant,
                        NOT a continued fraction. Taylor branches are exactly the
                        Taylor-method helpers per `small_angle_series.h` audit.

ERROR BOUND
  Bound category:       precision (Taylor truncation in the small-angle branches),
                        plus per-category propagation per REQ-EF-3 throughout.
  Bound formula:        Two truncation contributions when small-angle branches
                        fire:
                          - taylor_sinc: |θ|⁶ / 5040       (REQ-EF-6)
                          - taylor_cos_minus_sinc_over_θ²: |θ|⁴ / 840 (REQ-EF-6)
                        Both are added to result.errors.precision inside the
                        helpers (AUD-EF-5). For |θ| ≥ 1e-4 the closed-form
                        cos / sin / sqrt expressions propagate per REQ-EF-3.
  Bound implemented:    Helpers add `trunc_bound` to precision; the DQ-level
                        composition delegates everything else to TrackedValue<T>
                        ops. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-6 (Taylor branch precision bound),
                        REQ-DQ-18 (screw exp component form)
  AUD-EF applies:       AUD-EF-1, AUD-EF-2, AUD-EF-5 (Taylor branch wiring)
  AUD-MC applies:       AUD-MC-18 (exp ∘ log round-trip on the half-angle ball);
                        AUD-MC-15 (exp output is unit pose)
  Verification test:    tests/test_dual_quaternion/ —
                        - exp_screw(0, 0) == identity()
                        - magnitude_squared of exp_screw output is (1, 0)
                          (unit pose property)
                        - exp_screw round-trips with log_screw at small, medium,
                          large |u|
                        - matches independent matrix-exp(SE(3) hat) on test screws

NOTES
  - Headline match-verdict: ✓ — this is the screw exponential (Lie-group
    geodesic), not a linear interpolation in components. The cos, sinc, β
    factors are the closed-form result of the infinite series Σ Û^k/k!
    truncated only by the exact ε² = 0 relation, not by approximation.
  - Small-angle branches in `taylor_sinc` and `taylor_cos_minus_sinc_over_θ²` are
    Taylor (per §5 of the audit framework), with rigorous Leibniz / next-term
    bounds added to precision. No method-theory mismatch.
  - Dependence on `small_angle_series.h` (PASS) is documented in the file
    header (line 48) and inherited.
```

---

## Formula 26: `exp_screw(twist)` — overload

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::exp_screw_twist
Location:               src/math/dual_quaternion.h:297-299
Mathematical statement: exp_screw(Ŝ) for a pure DQ Ŝ; extract u = Ŝ.angular(),
                        v = Ŝ.linear(), delegate to Formula 25.

THEORY
  Underlying theorem:   Same screw exponential as Formula 25; this is a
                        convenience overload that accepts a DualQuaternion
                        directly. Scalar parts of the input are ignored as
                        documented (line 295–296).
  Primary reference:    Same as Formula 25.
  Domain of validity:   Pure DQ inputs (scalar parts ignored). For non-pure
                        inputs the scalar information is silently dropped,
                        consistent with the screw exponential acting on se(3).

METHOD
  Method declared:      Extract angular and linear vectors; delegate to the
                        (u, v) overload.
  Method implemented:   `return exp_screw(twist.angular(), twist.linear());` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       Same as Formula 25.
  Bound formula:        Inherited from Formula 25 via delegation.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       Same as Formula 25.
  AUD-EF applies:       Same as Formula 25.
  AUD-MC applies:       Same as Formula 25.
  Verification test:    tests/test_dual_quaternion/ — exp_screw(from_screw(u, v))
                        == exp_screw(u, v).

NOTES
  - Convenience wrapper; no new theory.
```

---

## Formula 27: `log_screw()` — screw logarithm

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::log_screw
Location:               src/math/dual_quaternion.h:314-333
Mathematical statement: Given unit pose M̂ = q_r + εq_d, return the pure DQ
                        Û = u + εv such that exp_screw(Û) = M̂. Component form:
                          u = q_r.log_unit()          (Rodrigues half-angle)
                          s := u·v = −q_d.scalar / sinc|u|
                          v = (q_d.vector − s β u) / sinc|u|
                        with β = (cos|u| − sinc|u|) / |u|².

THEORY
  Underlying theorem:   Inverse of the screw exponential (Formula 25) on the
                        half-angle ball |u| ≤ π/2. Solve the three component
                        equations from REQ-DQ-18:
                          (i)   q_r = exp_pure(u)             → u = log_unit(q_r)
                          (ii)  q_d.scalar = −(u·v) sinc|u|   → solve for u·v
                          (iii) q_d.vector = sinc|u| · v + (u·v) β u
                                                              → solve for v
                        Equations (ii) and (iii) are linear in v given u, with
                        coefficient sinc|u| bounded below by sinc(π/2) = 2/π
                        on the half-angle ball, so the divisions are bounded.
  Primary reference:    Murray, Li & Sastry (1994) §3.2 (logarithm and exponential
                        coordinates); Selig (2005) §9.4; the shortest-path
                        half-angle convention is in Quaternion<T>::log_unit
                        (quaternion.h:218-230).
  Domain of validity:   Unit M̂ within the half-angle ball |u| ≤ π/2. For
                        |u| → π the q_r.log_unit() singularity dominates;
                        for u → 0 the small-angle Taylor branches in helpers
                        keep divisions safe.

METHOD
  Method declared:      Closed-form inversion of REQ-DQ-18, using the Taylor
                        helpers (`taylor_sinc`, `taylor_cos_minus_sinc_over_θ²`)
                        for the small-angle branch.
  Method implemented:   Lines 314-332:
                          u = real.log_unit()
                          theta_sq = u·u
                          theta = sqrt(theta_sq)
                          sinc_t = taylor_sinc(theta, theta_sq)
                          beta = taylor_cos_minus_sinc_over_theta_sq(theta, theta_sq)
                          s = −dual.w / sinc_t          (inverts q_d.scalar)
                          sb = s · beta
                          v.x = (dual.x − sb · u.x) / sinc_t
                          v.y = (dual.y − sb · u.y) / sinc_t
                          v.z = (dual.z − sb · u.z) / sinc_t
                          return DualQuaternion(pure(u), pure(v))
                        Matches the inversion above. ✓
  Match verdict:        ✓ matched — direct inversion of the screw exp formulae,
                        NOT a series in M̂ − 1̂, NOT a Newton iteration.

ERROR BOUND
  Bound category:       precision (Taylor truncation when small-angle branches
                        fire), plus per-category propagation per REQ-EF-3.
  Bound formula:        Truncation contributions (when |θ| < 1e-4):
                          - taylor_sinc:                |θ|⁶ / 5040
                          - taylor_cos_minus_sinc_over_θ²: |θ|⁴ / 840
                        Plus the precision bound from Quaternion<T>::log_unit
                        on q_r (uses `taylor_half_angle_scale` from
                        `small_angle_series.h`, which adds 5|s|⁶/112 to
                        precision in its small-arg branch per the audit doc §5.2).
                        All errors propagated through TrackedValue<T> divisions
                        per REQ-EF-3.
  Bound implemented:    Helpers add their `trunc_bound`; divisions and products
                        propagate via TrackedValue<T> operators. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-6, REQ-DQ-18
  AUD-EF applies:       AUD-EF-1, AUD-EF-2, AUD-EF-5
  AUD-MC applies:       AUD-MC-18 (exp ∘ log == id on the half-angle ball;
                        log ∘ exp == id for pure screws with |u| ≤ π/2)
  Verification test:    tests/test_dual_quaternion/ —
                        - log_screw(identity()) ≈ (0, 0)
                        - log_screw(exp_screw(u, v)) ≈ (u, v) on the half-angle ball
                        - exp_screw(log_screw(M̂)) ≈ M̂ for random unit poses

NOTES
  - The half-angle ball constraint comes from Quaternion<T>::log_unit's
    shortest-path convention (flip q → −q when w < 0); this keeps |u| ≤ π/2
    and sinc|u| ≥ 2/π, so all divisions are bounded away from zero.
  - Dependence on small_angle_series.h (PASS) and quaternion.h::log_unit
    (pending audit) is documented in the file header.
```

---

## Formula 28: `sclerp(m0, m1, t)` — screw geodesic interpolation

```
=== FORMULA AUDIT CARD ===
ID:                     dual_quaternion::sclerp
Location:               src/math/dual_quaternion.h:345-353
Mathematical statement: M̂(t) = m0 · exp_screw(t · log_screw(m0⁻¹ · m1))
                        for t ∈ [0, 1]; gives the SE(3) Lie-group geodesic
                        from m0 (t = 0) to m1 (t = 1).

THEORY
  Underlying theorem:   Lie-group geodesic on SE(3) (left-invariant
                        connection): given two unit poses m0 and m1, the
                        constant-screw-rate geodesic is
                          M̂(t) = m0 · exp(t · log(m0⁻¹ m1)),
                        which interpolates uniformly along a screw axis. This
                        is the SE(3) analogue of slerp on SO(3) (Shoemake 1985)
                        and is the "ScLERP" of Kavan, Collins, Žára & O'Sullivan
                        (2008).
                        **NOT** a linear interpolation in component space —
                        that would give an off-manifold result with non-unit
                        magnitude that doesn't even encode a rigid motion at
                        intermediate t.
  Primary reference:    Kavan, Collins, Žára & O'Sullivan (2008) §5 "ScLERP";
                        Murray, Li & Sastry (1994) §3.2 (exponential
                        coordinates); Munthe-Kaas (1999) (Lie-group methods
                        more broadly).
  Domain of validity:   m0, m1 both unit (preconditioned). The delta
                        m0⁻¹ m1 should be within the half-angle ball for
                        single-valued logarithm; if its rotation is > π
                        the result silently takes the shortest-path inverse
                        via Quaternion::log_unit's q → −q flip.

METHOD
  Method declared:      Lie-group geodesic via exp ∘ scale ∘ log of the relative
                        pose: delta = m0⁻¹ m1; step = exp_screw(t · log_screw(delta));
                        result = m0 · step.
  Method implemented:   Lines 345-353:
                          delta = m0.inverse() · m1
                          log_delta = delta.log_screw()
                          u_scaled = t · log_delta.angular()
                          v_scaled = t · log_delta.linear()
                          return m0 · exp_screw(u_scaled, v_scaled)
                        Matches the algorithm. ✓
  Match verdict:        ✓ matched — screw geodesic, NOT linear interpolation
                        in components. The exp_screw / log_screw cycle is the
                        SE(3) Lie-group machinery, and m0 · step is left-
                        translation back to the m0-anchored frame.

ERROR BOUND
  Bound category:       precision (Taylor truncation from inner exp/log helpers),
                        plus per-category propagation per REQ-EF-3.
  Bound formula:        Composed:
                          - inverse(): no new error beyond REQ-EF-3
                          - operator*: REQ-EF-3 through Hamilton products
                          - log_screw: Taylor branches (helpers add to precision)
                          - scalar multiplication t · (...): REQ-EF-3
                          - exp_screw: Taylor branches (helpers add to precision)
                          - final m0 · step: REQ-EF-3
                        Final bound = sum of inherited bounds along this chain.
  Bound implemented:    Each step delegates; sclerp itself adds no new term. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-6
  AUD-EF applies:       AUD-EF-1, AUD-EF-2, AUD-EF-5 (Taylor branches inside
                        the chain)
  AUD-MC applies:       AUD-MC-17 (sclerp(m0, m1, 0) == m0, sclerp(m0, m1, 1) == m1);
                        AUD-MC-18 (constant screw rate: dM̂/dt at parameter τ
                        gives the same body-frame screw for all τ)
  Verification test:    tests/test_dual_quaternion/ —
                        - sclerp(m0, m1, 0) ≈ m0
                        - sclerp(m0, m1, 1) ≈ m1
                        - sclerp(m0, m1, 0.5) on commuting screws matches
                          half-axis-half-translation expectation
                        - midpoint independence: sclerp(m0, m1, 0.5) ==
                          sclerp(m1, m0, 0.5) (with appropriate orientation
                          handling of the shortest-path log)

NOTES
  - Headline match-verdict (confirmed): this is the **screw geodesic** on SE(3),
    not a componentwise linear blend. Kavan et al. (2008) introduced ScLERP
    explicitly to avoid the off-manifold drift that componentwise blending
    produces.
  - Inherits Taylor-branch correctness from exp_screw / log_screw (Formulas
    25, 27), which in turn inherit from small_angle_series.h (PASS).
  - The Lie-group interpretation: m0⁻¹ m1 expresses m1 in the m0-anchored
    body frame; t·log scales the screw uniformly; exp recovers a partial
    pose on the geodesic; left-multiply by m0 returns the result to the
    world frame.
```

---

## File-level verdict

**A. Error wiring** (AUD-EF-1, AUD-EF-2, AUD-EF-5, AUD-EF-7):
- All 28 ops return either `DualQuaternion<T>` (a pair of `Quaternion<T>`, each four `TrackedValue<T>`), `Vector3<T>` (three `TrackedValue<T>`), `Quaternion<T>` (four `TrackedValue<T>`), `DualNumber<T>` (two `TrackedValue<T>`), or trivial accessors of the above. No bare `T` slips through. ✓
- The two Lie-algebra ops with truncation (`exp_screw`, `log_screw`) inherit Taylor-bound additions via `small_angle_series.h`. ✓
- `sclerp` inherits through its chain (no own truncation). ✓

**B. Algebra axioms** (AUD-MC-3..18 as cross-referenced per card):
- Additive group on H ⊕ εH (Formulas 11–13): ✓
- Dual-quaternion product (Formula 14): associative, distributive, identity-neutral, non-commutative on factors; truncation ε² = 0 is exact. ✓
- Conjugates (Formulas 17–19): involutions; combined conjugate is the composition of the other two. ✓
- Unit-pose inverse via conjugate (Formula 21): ✓
- SE(3) point action via the M̂ P̂ M̂♯ sandwich (Formula 23): ✓
- Lie-algebra exp / log: closed-form on the half-angle ball, exact round-trip up to TrackedValue<T> propagation + Taylor-truncation bound (Formulas 25–28). ✓

**C. Theoretical basis** (this document):
- All 28 formulas cite a primary source (Study 1903; Clifford 1873; Murray-Li-Sastry 1994; Selig 2005; Kavan et al. 2008; Munthe-Kaas 1999 for retraction; small_angle_series helpers per §5 of the framework). ✓
- All 28 method declarations match the implementation. ✓
- All 28 bound declarations match the implemented per-category error accumulation. ✓
- **Headline match-verdict on `exp_screw` / `log_screw` / `sclerp`**: ✓ **screw exponential** (Lie-group exp on SE(3)) and screw geodesic, NOT linear interpolation in components, NOT Padé approximant, NOT continued fraction. The closed-form factors (cos, sinc, β) follow directly from Σ Û^k / k! truncated only by the exact ε² = 0 relation; small-angle Taylor branches are bounded by Leibniz / next-term formulas in `small_angle_series.h` (PASS per audit doc §5).

**Overall file verdict: PASS**

- Method-theory matches: 28 / 28 ✓
- Bound-method matches: 28 / 28 ✓
- No Padé / continued-fraction / linear-blend masquerading as Taylor or screw exponential.

**Dependencies**:
- Formulas 25, 26, 27 depend on `small_angle_series.h::taylor_sinc` and
  `taylor_cos_minus_sinc_over_theta_sq` (both PASS per audit framework §5).
- Formula 27 also depends on `quaternion.h::log_unit`, which in turn uses
  `small_angle_series.h::taylor_half_angle_scale` (PASS with note for
  wide-T tightening).
- Formula 28 (sclerp) depends on Formulas 14, 21, 25, 27, so its bound is
  inherited along that chain.

---

## Verification checklist

- [ ] tests/test_math/test_dual_quaternion.cc or equivalent:
  - [ ] Default-constructor and named factories (zero, identity)
  - [ ] from_pose / translation round-trip (Formulas 5, 8)
  - [ ] from_screw / angular / linear round-trip (Formulas 6, 9, 10)
  - [ ] Additive group axioms (Formulas 11–13)
  - [ ] Product associativity, identity, distributivity (Formula 14)
  - [ ] Scalar–DQ commutativity left/right (Formulas 15, 16)
  - [ ] All three conjugates: involution, combined = comp(real conj, dual conj)
        (Formulas 17–19)
  - [ ] magnitude_squared multiplicativity and unit-pose value (Formula 20)
  - [ ] inverse * M̂ == identity on unit poses (Formula 21)
  - [ ] normalized: post-normalization |q_r| = 1 and <q̂_r, q̂_d> = 0
        (Formula 22)
  - [ ] apply / apply_direction: composition consistency, identity action
        (Formulas 23, 24)
  - [ ] exp_screw: matches Quaternion::exp_pure on the rotation alone,
        agrees with matrix exp(SE(3) hat) on test screws (Formula 25)
  - [ ] log_screw inverts exp_screw on the half-angle ball (Formulas 25 ↔ 27)
  - [ ] sclerp: endpoint conditions, constant-screw-rate property
        (Formula 28)
  - [ ] Cross-type tests (double, cpp_bin_float_50)

---

**Document**: `design/audit/theoretical_basis_audit/dual_quaternion.md`
**Status**: **PASS**
**Date**: 2026-05-13
**Analyst**: internal review
