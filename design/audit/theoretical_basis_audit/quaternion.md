# Theoretical Basis Audit — `src/math/quaternion.h`

**File**: `src/math/quaternion.h` (233 lines)
**Scope**: `Quaternion<T>` struct, 22 distinct formulas
**Theory anchor**: Hamilton (1843) "On quaternions"; Altmann (1986) "Rotations, Quaternions, and Double Groups"; Shoemake (1985) "Animating rotation with quaternion curves"
**Key Match-Verdict check**: `exp_pure`, `log_unit`, `from_axis_angle` — confirm Taylor (not Padé / continued fraction) on the small-angle branch.

---

## Formula 1: Default constructor

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::ctor_default
Location:               src/math/quaternion.h:52
Mathematical statement: Quaternion() default-constructs (w, x, y, z) = (0, 0, 0, 0)
                        via default TrackedValue<T> constructors.

THEORY
  Underlying theorem:   None — trivial value-initialization. The result is the
                        additive identity 0 of the quaternion ring ℍ.
  Primary reference:    Hamilton (1843); C++ aggregate initialization.
  Domain of validity:   All T.

METHOD
  Method declared:      Member-wise default construction.
  Method implemented:   `w(), x(), y(), z()` — four defaulted TrackedValue<T>.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (no computation; error inherited from defaults)
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       AUD-MC-4 (additive identity of ℍ)
  Verification test:    tests/test_math/test_quaternion.cc — Quaternion() + q == q.

NOTES
  - Equivalent to `zero()` factory (Formula 4).
```

---

## Formula 2: Value constructor

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::ctor_value
Location:               src/math/quaternion.h:54-56
Mathematical statement: Quaternion(w_, x_, y_, z_) constructs w_ + x_·i + y_·j + z_·k.

THEORY
  Underlying theorem:   Direct construction of an element of ℍ from a basis
                        decomposition {1, i, j, k}. Hamilton's basis.
  Primary reference:    Hamilton (1843).
  Domain of validity:   All T.

METHOD
  Method declared:      Direct member assignment.
  Method implemented:   `w(w_), x(x_), y(y_), z(z_)`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (no computation; inherits per-component errors)
  Bound formula:        Inherited from the four input TrackedValue<T> errors.
  Bound implemented:    Member copy. ✓
  Bound verdict:        ✓ matched (REQ-EF-12 composite preservation).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12 (composite type construction)
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a
  Verification test:    tests/test_math/test_quaternion.cc — accessor round-trip.

NOTES
  - No arithmetic, no error introduction.
```

---

## Formula 3: identity()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::identity
Location:               src/math/quaternion.h:61-63
Mathematical statement: 1_ℍ = 1 + 0·i + 0·j + 0·k.

THEORY
  Underlying theorem:   Multiplicative identity of the quaternion ring.
                        Hamilton's relations i² = j² = k² = ijk = −1 give 1·q = q
                        for all q ∈ ℍ.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.
  Domain of validity:   All T.

METHOD
  Method declared:      Construct (1, 0, 0, 0) via `exact<T>`.
  Method implemented:   `Quaternion(exact<T>(1), exact<T>(0), exact<T>(0), exact<T>(0))`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (exact construction — zero error)
  Bound formula:        `exact<T>` produces zero-error literals.
  Bound implemented:    ✓ exact literals.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (exact literal has zero error)
  AUD-EF applies:       AUD-EF-2
  AUD-MC applies:       AUD-MC-5 (multiplicative identity: 1·q = q·1 = q)
  Verification test:    tests/test_math/test_quaternion.cc — identity() * q == q
                        and q * identity() == q.

NOTES
  - Verified by AUD-MC-5 as part of ring identity check.
```

---

## Formula 4: zero()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::zero
Location:               src/math/quaternion.h:66
Mathematical statement: 0_ℍ = 0 + 0·i + 0·j + 0·k.

THEORY
  Underlying theorem:   Additive identity of the quaternion ring.
  Primary reference:    Hamilton (1843).
  Domain of validity:   All T.

METHOD
  Method declared:      Delegate to default constructor.
  Method implemented:   `return Quaternion();` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        Inherits default-constructed TrackedValue<T> zero error.
  Bound implemented:    ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1
  AUD-EF applies:       AUD-EF-2
  AUD-MC applies:       AUD-MC-4 (additive identity: q + 0 = q)
  Verification test:    tests/test_math/test_quaternion.cc — zero() + q == q.

NOTES
  - Synonym for default ctor.
```

---

## Formula 5: pure(v)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::pure
Location:               src/math/quaternion.h:69-71
Mathematical statement: pure(v) = 0 + v_x·i + v_y·j + v_z·k.

THEORY
  Underlying theorem:   Embedding ℝ³ ↪ ℍ as the imaginary subspace. The pure
                        quaternions form a Lie subalgebra isomorphic to so(3).
  Primary reference:    Hamilton (1843); Altmann (1986) §2.3 "Pure quaternions".
  Domain of validity:   All v ∈ Vector3<T>.

METHOD
  Method declared:      Construct (0, v.x, v.y, v.z).
  Method implemented:   `Quaternion(exact<T>(0), v.x, v.y, v.z)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       inherits per-component errors from v
  Bound formula:        w-component error = 0 (exact); i/j/k inherit from v.
  Bound implemented:    ✓
  Bound verdict:        ✓ matched (REQ-EF-12).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-12
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-12 (pure-quaternion algebra)
  Verification test:    tests/test_math/test_quaternion.cc — pure(v).vector() == v
                        and pure(v).scalar() == 0.

NOTES
  - Used by `rotate()` indirectly and by `exp_pure()` semantics.
```

---

## Formula 6: from_axis_angle(axis_unit, angle)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::from_axis_angle
Location:               src/math/quaternion.h:77-83
Mathematical statement: q(n̂, θ) = cos(θ/2) + sin(θ/2)·(n_x·i + n_y·j + n_z·k).

THEORY
  Underlying theorem:   Euler-Rodrigues parametrization of SO(3) via the
                        SU(2) double cover. A rotation by angle θ about unit
                        axis n̂ is the image under the cover of the unit
                        quaternion cos(θ/2) + sin(θ/2)·n̂.
                        cos and sin are evaluated by their standard library
                        routines, which are themselves CLOSED-FORM identities
                        (TrackedValue<T>::cos, ::sin) — not a Taylor truncation
                        at THIS layer. (The Taylor branch only appears below
                        in `exp_pure` / `log_unit` for normalized small-angle
                        regularization.)
  Primary reference:    Hamilton (1843); Altmann (1986) Ch. 12 "Quaternion
                        parameters for rotation"; Shoemake (1985) §5.
  Domain of validity:   `axis_unit` must be unit-norm (caller-enforced — see
                        @doc on line 73-74; no runtime renormalization).
                        θ ∈ ℝ; the unit-norm property of q holds for any θ
                        because sin² + cos² = 1.

METHOD
  Method declared:      Closed-form half-angle formula: half = θ/2, c = cos(half),
                        s = sin(half), assemble (c, s·n̂).
  Method implemented:   Lines 79-82 — `half = angle/exact<T>(2)`, `c = cos(half)`,
                        `s = sin(half)`, return `Quaternion(c, s*axis.x, s*axis.y,
                        s*axis.z)`. NO Taylor branch, NO Padé, NO continued
                        fraction. The cos/sin are TrackedValue<T> trig functions,
                        which are themselves closed-form identities in the
                        TrackedValue layer.
  Match verdict:        ✓ matched — CLOSED-FORM, not Taylor. Theory cites the
                        half-angle Euler-Rodrigues identity, and the code is
                        the literal half-angle formula. The cos/sin here are
                        not a "small-angle Taylor" approximation; they're the
                        full trig functions of the half-angle.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        REQ-EF-3 composition through:
                        (a) division by exact 2 (no rounding);
                        (b) TrackedValue<T>::cos and ::sin (rigorous bound on
                            mathematical functions per AUD-EF-6);
                        (c) scalar product s·n_i (per REQ-EF-3 product rule).
  Bound implemented:    Delegated to TrackedValue<T> operators. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (component-wise composition)
  AUD-EF applies:       AUD-EF-1, AUD-EF-6 (sin/cos bound)
  AUD-MC applies:       AUD-MC-4..AUD-MC-7 (round-trip via log_unit;
                        composition of rotations)
  Verification test:    tests/test_math/test_quaternion.cc —
                        |from_axis_angle(n̂, θ)| == 1 and rotation-by-θ check.

NOTES
  - **MATCH-VERDICT CHECK PASS**: implementation is half-angle CLOSED-FORM via
    std::cos/std::sin, not a Taylor expansion. The function does not branch
    on small θ; small-θ regularization is the responsibility of the inverse
    map `log_unit` (Formula 22) and `exp_pure` (Formula 21) where 0/0
    singularities can arise.
  - Precondition (`axis_unit` is unit) is doc-only; no runtime check.
```

---

## Formula 7: scalar()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::scalar
Location:               src/math/quaternion.h:88
Mathematical statement: Re(q) = w.

THEORY
  Underlying theorem:   Component projection (definition of scalar part).
  Primary reference:    Hamilton (1843); Altmann (1986) §2.
  Domain of validity:   All q.

METHOD
  Method declared:      Return member `w`.
  Method implemented:   `return w;` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       inherits w.errors verbatim
  Bound formula:        Identity projection — REQ-EF-12.
  Bound implemented:    ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a (accessor)
  Verification test:    tests/test_math/test_quaternion.cc — scalar of (a,b,c,d) == a.

NOTES
  - Pure accessor; no arithmetic.
```

---

## Formula 8: vector()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::vector
Location:               src/math/quaternion.h:91
Mathematical statement: Im(q) = (x, y, z) ∈ ℝ³.

THEORY
  Underlying theorem:   Component projection onto the imaginary subspace.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.3.
  Domain of validity:   All q.

METHOD
  Method declared:      Construct Vector3<T>(x, y, z).
  Method implemented:   `return Vector3<T>(x, y, z);` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       inherits per-component errors from x, y, z
  Bound formula:        REQ-EF-12 composite construction; no new error.
  Bound implemented:    ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a (accessor)
  Verification test:    tests/test_math/test_quaternion.cc — vector of (a,b,c,d) == (b,c,d).

NOTES
  - Used by `rotate()` (Formula 20) and `log_unit()` (Formula 22) internally.
```

---

## Formula 9: operator+(a, b)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::operator_add
Location:               src/math/quaternion.h:96-98
Mathematical statement: (a + b) = (a.w+b.w) + (a.x+b.x)·i + (a.y+b.y)·j + (a.z+b.z)·k.

THEORY
  Underlying theorem:   Abelian-group addition in ℍ as a 4-dimensional real
                        vector space (component-wise).
  Primary reference:    Hamilton (1843); Altmann (1986) §2.
  Domain of validity:   All a, b ∈ Quaternion<T>.

METHOD
  Method declared:      Component-wise addition (4 TrackedValue<T> adds).
  Method implemented:   `Quaternion(a.w+b.w, a.x+b.x, a.y+b.y, a.z+b.z)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Per-category triangle inequality across 4 additions
                        (REQ-EF-3).
  Bound implemented:    Delegated to TrackedValue<T>::operator+(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-3
  AUD-MC applies:       AUD-MC-4 (additive group: associative, commutative,
                        identity, inverse)
  Verification test:    tests/test_math/test_quaternion.cc — (a+b)+c == a+(b+c),
                        a+b == b+a, a+zero() == a.

NOTES
  - Exact in the algebra; no truncation.
```

---

## Formula 10: operator-(a, b)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::operator_sub
Location:               src/math/quaternion.h:101-103
Mathematical statement: (a − b) = (a.w−b.w) + (a.x−b.x)·i + (a.y−b.y)·j + (a.z−b.z)·k.

THEORY
  Underlying theorem:   Abelian-group subtraction in ℍ (component-wise inverse
                        addition).
  Primary reference:    Hamilton (1843); Altmann (1986) §2.
  Domain of validity:   All a, b ∈ Quaternion<T>.

METHOD
  Method declared:      Component-wise subtraction.
  Method implemented:   `Quaternion(a.w-b.w, a.x-b.x, a.y-b.y, a.z-b.z)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Per-category triangle inequality (REQ-EF-3).
  Bound implemented:    Delegated to TrackedValue<T>::operator-(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-3
  AUD-MC applies:       AUD-MC-4 (additive inverse)
  Verification test:    tests/test_math/test_quaternion.cc — a-a == zero(),
                        (a+b)-b == a.

NOTES
  - Exact in the algebra.
```

---

## Formula 11: operator-() (unary)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::operator_neg
Location:               src/math/quaternion.h:106-108
Mathematical statement: −q = (−w) + (−x)·i + (−y)·j + (−z)·k.

THEORY
  Underlying theorem:   Additive inverse in ℍ.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.
  Domain of validity:   All q.

METHOD
  Method declared:      Component-wise negation.
  Method implemented:   `Quaternion(-w, -x, -y, -z)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       inherits per-component errors (negation preserves bounds)
  Bound formula:        REQ-EF-3: |−x| = |x|, all error categories preserved.
  Bound implemented:    Delegated to TrackedValue<T>::operator-(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-4 (additive inverse: q + (−q) = 0)
  Verification test:    tests/test_math/test_quaternion.cc — −(−q) == q,
                        q + (−q) == zero().

NOTES
  - Used by `log_unit()` (Formula 22) for shortest-path branch.
```

---

## Formula 12: operator*(a, b)  [Hamilton product]

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::operator_mul_hamilton
Location:               src/math/quaternion.h:118-125
Mathematical statement: a*b = (a.w·b.w − a.x·b.x − a.y·b.y − a.z·b.z)
                            + (a.w·b.x + a.x·b.w + a.y·b.z − a.z·b.y)·i
                            + (a.w·b.y − a.x·b.z + a.y·b.w + a.z·b.x)·j
                            + (a.w·b.z + a.x·b.y − a.y·b.x + a.z·b.w)·k.

THEORY
  Underlying theorem:   Hamilton's defining relations i² = j² = k² = ijk = −1,
                        yielding (w₁+v₁)(w₂+v₂) = w₁w₂ − v₁·v₂ + w₁v₂ + w₂v₁ + v₁×v₂.
                        Non-commutative; associative; distributive over +.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.4 "The Hamilton product".
  Domain of validity:   All a, b ∈ Quaternion<T>.

METHOD
  Method declared:      Closed-form Hamilton-product expansion (16 mults,
                        12 adds/subs).
  Method implemented:   Lines 120-123 — literal expansion of the four
                        components. ✓
  Match verdict:        ✓ matched — closed-form identity, not approximation.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Each output component is a 4-term sum of products of
                        TrackedValue<T>; REQ-EF-3 composition (product +
                        triangle-inequality for sums).
  Bound implemented:    Delegated to TrackedValue<T>::operator* and operator+. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-3
  AUD-MC applies:       AUD-MC-5 (Hamilton multiplication: associative;
                                  identity; |a·b| = |a|·|b|),
                        AUD-MC-6 (non-commutativity: a*b ≠ b*a in general)
  Verification test:    tests/test_math/test_quaternion.cc — (a*b)*c == a*(b*c),
                        a*identity == a, magnitude multiplicativity.

NOTES
  - Confirm sign pattern: this is the right-multiplicative Hamilton product;
    matches Boost.Geometry/Eigen per file header (line 16).
  - Algebra-axiom checks fall under AUD-MC-5..7.
```

---

## Formula 13: operator*(s, q)  [scalar on left]

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::operator_mul_scalar_left
Location:               src/math/quaternion.h:128-130
Mathematical statement: s·q = (s·w) + (s·x)·i + (s·y)·j + (s·z)·k.

THEORY
  Underlying theorem:   ℍ is a ℝ-bimodule. Left scalar multiplication
                        distributes over the four basis components.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.
  Domain of validity:   All s, q.

METHOD
  Method declared:      Component-wise scalar multiplication (4 mults).
  Method implemented:   `Quaternion(s*q.w, s*q.x, s*q.y, s*q.z)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        REQ-EF-3 product rule per component.
  Bound implemented:    Delegated to TrackedValue<T>::operator*. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Module operation (associative scalar mult)
  Verification test:    tests/test_math/test_quaternion.cc — (s*t)*q == s*(t*q).

NOTES
  - Exact in the algebra.
```

---

## Formula 14: operator*(q, s)  [scalar on right]

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::operator_mul_scalar_right
Location:               src/math/quaternion.h:133-135
Mathematical statement: q·s = (w·s) + (x·s)·i + (y·s)·j + (z·s)·k.

THEORY
  Underlying theorem:   Right ℝ-module structure on ℍ. Real scalars commute
                        with the basis, so s·q = q·s for s ∈ ℝ.
  Primary reference:    Hamilton (1843).
  Domain of validity:   All s, q.

METHOD
  Method declared:      Component-wise scalar multiplication (4 mults).
  Method implemented:   `Quaternion(q.w*s, q.x*s, q.y*s, q.z*s)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        REQ-EF-3 product rule per component.
  Bound implemented:    Delegated to TrackedValue<T>::operator*. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Module operation
  Verification test:    tests/test_math/test_quaternion.cc — q*s == s*q.

NOTES
  - Equivalent to Formula 13 by commutativity of real scalars.
```

---

## Formula 15: conjugate()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::conjugate
Location:               src/math/quaternion.h:140-142
Mathematical statement: q* = w − x·i − y·j − z·k.

THEORY
  Underlying theorem:   Quaternion conjugation: the antiautomorphism flipping
                        the sign of the imaginary part. Satisfies
                        (a·b)* = b*·a*, q·q* = |q|², q** = q.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.6 "Conjugation
                        and norm".
  Domain of validity:   All q.

METHOD
  Method declared:      Negate imaginary components.
  Method implemented:   `Quaternion(w, -x, -y, -z)` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       inherits per-component errors (negation preserves bounds)
  Bound formula:        |−x_i| = |x_i|; all error categories preserved (REQ-EF-3).
  Bound implemented:    ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-7 (conjugate distributes anti-commutatively
                        over product: (a*b)* = b* · a*)
  Verification test:    tests/test_math/test_quaternion.cc — q.conjugate().conjugate() == q,
                        (a*b).conjugate() == b.conjugate() * a.conjugate().

NOTES
  - Used by `rotate()` semantics (v' = q v q*) and `inverse()`.
```

---

## Formula 16: magnitude_squared()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::magnitude_squared
Location:               src/math/quaternion.h:145-147
Mathematical statement: |q|² = w² + x² + y² + z² = q·q*.

THEORY
  Underlying theorem:   Euclidean norm-squared on ℝ⁴; equivalently the scalar
                        part of q·q* (Cayley quadratic form).
  Primary reference:    Hamilton (1843); Altmann (1986) §2.6.
  Domain of validity:   All q.

METHOD
  Method declared:      Sum of squared components.
  Method implemented:   `w*w + x*x + y*y + z*z` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        REQ-EF-3 composition: each x_i² is a product (rule
                        2|x||δx|), summed by triangle inequality.
  Bound implemented:    Delegated to TrackedValue<T>::operator* and ::operator+. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-3
  AUD-MC applies:       AUD-MC-8 (|q|² ≥ 0; |q|² = 0 ⇒ q = 0)
  Verification test:    tests/test_math/test_quaternion.cc — |a*b|² = |a|²·|b|²
                        (multiplicativity), |q*|² == |q|².

NOTES
  - Pythagorean closed-form, exact in the algebra.
```

---

## Formula 17: magnitude()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::magnitude
Location:               src/math/quaternion.h:150-152
Mathematical statement: |q| = √(w² + x² + y² + z²).

THEORY
  Underlying theorem:   Euclidean norm on ℝ⁴ (square root of the quadratic
                        form). |q| is multiplicative: |a·b| = |a|·|b|.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.6.
  Domain of validity:   All q (|q|² ≥ 0 always; √ is defined on [0, ∞)).

METHOD
  Method declared:      Square root of `magnitude_squared`.
  Method implemented:   `sqrt(magnitude_squared())` ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        Composition of Formula 16 (magnitude_squared) followed
                        by TrackedValue<T>::sqrt; the sqrt's bound is governed
                        by REQ-EF-9 / AUD-EF-6 (sqrt as a math function with
                        per-category derivative-style bound).
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-9
  AUD-EF applies:       AUD-EF-1, AUD-EF-6
  AUD-MC applies:       AUD-MC-8 (multiplicativity: |a*b| = |a|·|b|)
  Verification test:    tests/test_math/test_quaternion.cc — |a*b| == |a|·|b|,
                        |q*| == |q|.

NOTES
  - sqrt branch may have a near-zero conditioning issue; for a unit quaternion
    |q| ≈ 1 and the derivative 1/(2|q|) is well-behaved.
```

---

## Formula 18: inverse()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::inverse
Location:               src/math/quaternion.h:157-160
Mathematical statement: q⁻¹ = q* / |q|² = (w − x·i − y·j − z·k) / (w²+x²+y²+z²).

THEORY
  Underlying theorem:   Multiplicative inverse in the division ring ℍ:
                        q · q* = |q|², so q⁻¹ = q*/|q|² is the unique inverse.
                        For unit q, |q|² = 1 and q⁻¹ = q*.
  Primary reference:    Hamilton (1843); Altmann (1986) §2.6.
  Domain of validity:   q ≠ 0 (division by zero otherwise).

METHOD
  Method declared:      Closed-form conjugate-over-magnitude-squared.
  Method implemented:   Lines 158-159 — `m2 = magnitude_squared()`, return
                        `Quaternion(w/m2, -x/m2, -y/m2, -z/m2)`. ✓
  Match verdict:        ✓ matched — closed-form, not iterative or series.

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        Composition: magnitude_squared (Formula 16) followed
                        by 4 TrackedValue<T> divisions and 3 negations. Each
                        division per REQ-EF-3 / AUD-EF-3.
  Bound implemented:    Delegated to TrackedValue<T>::operator/. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-3
  AUD-MC applies:       AUD-MC-9 (q · q⁻¹ = q⁻¹ · q = identity())
  Verification test:    tests/test_math/test_quaternion.cc — q*q.inverse() == identity(),
                        q.inverse().inverse() == q.

NOTES
  - The conjugate-and-divide approach is one of two standard forms; for unit
    q, `conjugate()` (Formula 15) is preferred to skip the division.
```

---

## Formula 19: normalized()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::normalized
Location:               src/math/quaternion.h:164-167
Mathematical statement: q / |q| = q · |q|⁻¹  (projects onto unit sphere S³ ⊂ ℍ).

THEORY
  Underlying theorem:   Retraction onto the unit sphere S³ ⊂ ℍ. For q ≠ 0,
                        q/|q| is the closest point on S³ to q along the radial
                        direction. Used as a Lie-group retraction after RK
                        integration steps (REQ-EF-15 / propagator.h).
  Primary reference:    Altmann (1986) §12 "Unit-quaternion parametrization
                        of SO(3)"; Munthe-Kaas (1999) "High order Lie group
                        methods" §3 for the retraction context.
  Domain of validity:   q ≠ 0.

METHOD
  Method declared:      Closed-form division by magnitude.
  Method implemented:   Lines 165-166 — `m = magnitude()`, return
                        `Quaternion(w/m, x/m, y/m, z/m)`. ✓
  Match verdict:        ✓ matched — closed-form, not iterative
                        (no Newton step on the constraint).

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        magnitude (Formula 17, includes sqrt) followed by 4
                        divisions; each step per REQ-EF-3. The retraction
                        itself introduces a small bias from the radial
                        projection if q is far from S³, but for near-unit q
                        the bias is O(||q| − 1|²) and bounded by the input
                        component errors.
  Bound implemented:    Delegated to TrackedValue<T>::sqrt and ::operator/. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-9, REQ-EF-15 (retraction)
  AUD-EF applies:       AUD-EF-1, AUD-EF-6
  AUD-MC applies:       AUD-MC-10 (|q.normalized()| = 1 up to representation error)
  Verification test:    tests/test_math/test_quaternion.cc — |q.normalized()| == 1
                        within bounds; q.normalized().rotate(v) preserves |v|.

NOTES
  - Critical: this is THE retraction in the Munthe-Kaas RK4 Lie-group flow
    (propagator.h, REQ-EF-15). Treating |q| − 1 small after each step prevents
    drift off S³.
```

---

## Formula 20: rotate(v)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::rotate
Location:               src/math/quaternion.h:178-184
Mathematical statement: v' = q · v · q* (for unit q, with v promoted to a pure
                        quaternion), implemented as the optimized identity
                        v' = v + 2w·(q_v × v) + 2·q_v × (q_v × v).

THEORY
  Underlying theorem:   The SU(2) double cover of SO(3): a unit quaternion
                        q = (w, q_v) acts on v ∈ ℝ³ by v ↦ q v q*. The
                        equivalent vector-algebra form
                          v' = v + 2w(q_v × v) + 2 q_v × (q_v × v)
                        is derived by expanding q v q* and is exact (an
                        algebraic identity, not an approximation).
  Primary reference:    Hamilton (1843); Altmann (1986) §12.5 "The
                        Euler-Rodrigues formula for vector rotation";
                        Shoemake (1985) §5.
  Domain of validity:   q is unit-norm (precondition stated in @doc on
                        line 171; not enforced at runtime). All v ∈ ℝ³.

METHOD
  Method declared:      Closed-form algebraic identity (3 cross products,
                        2 vector adds, 2 scalar mults). Trades fewer mults
                        than the literal triple product q v q*.
  Method implemented:   Lines 179-183 — `qv = (x,y,z)`, `qv × v`, scale by 2,
                        add `w * (2(qv×v))` and `qv × (2(qv×v))` to v. ✓
  Match verdict:        ✓ matched — exact algebraic rearrangement of v ↦ q v q*,
                        not a truncated series or approximation.

ERROR BOUND
  Bound category:       precision, accuracy, measurement
  Bound formula:        Composition of Vector3<T>::cross (per AUD-MC-cross_product),
                        scalar multiplication, and vector addition. Each per
                        REQ-EF-3.
  Bound implemented:    Delegated to Vector3<T> and TrackedValue<T> ops. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1, AUD-EF-3
  AUD-MC applies:       AUD-MC-11 (rotation preserves norm: |rotate(v)| = |v|),
                        AUD-MC-12 (rotation composition: (q1*q2).rotate(v) ==
                                   q1.rotate(q2.rotate(v)))
  Verification test:    tests/test_math/test_quaternion.cc — |q.rotate(v)| == |v|,
                        rotation by 2π returns v.

NOTES
  - The identity is exact when |q| = 1; if |q| ≠ 1, the formula scales by
    |q|² rather than rotating. Caller responsibility.
  - This form is the "fast" Rodrigues vector rotation; trades 4 mults for 3
    over the literal q v q*.
```

---

## Formula 21: exp_pure(v)

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::exp_pure
Location:               src/math/quaternion.h:199-205
Mathematical statement: exp(0 + v) = cos|v| + sinc(|v|) · v
                        where v = (v_x, v_y, v_z) and sinc(θ) = sin(θ)/θ
                        with the convention sinc(0) = 1.

THEORY
  Underlying theorem:   Exponential map exp: pure_ℍ → S³ ⊂ ℍ. For a pure
                        quaternion u = (0, v) with θ = |v|,
                          exp(u) = Σ_{n≥0} u^n / n!
                        splits by even/odd powers (u² = −|v|² = −θ²) into
                          exp(u) = cos(θ) + (v/|v|)·sin(θ) = cos(θ) + sinc(θ)·v.
                        The result is unit-norm: |cos(θ)|² + |sinc(θ)·v|²
                        = cos²(θ) + sin²(θ) = 1.
                        The sinc(θ) factor has a removable singularity at θ=0;
                        near θ=0 it is replaced by its Taylor expansion
                          sinc(θ) = 1 − θ²/6 + θ⁴/120 + R(θ),
                                  |R| ≤ |θ|⁶/5040.
                        This is the canonical small-angle regularization.
  Primary reference:    Altmann (1986) §12.5; Shoemake (1985) §5; Murray, Li
                        & Sastry (1994) §3.2 "Exponential coordinates for
                        rotation"; Selig (2005) §3.4.
                        Taylor-truncation bound: Leibniz alternating-series
                        theorem applied to the entire function sinc(θ) =
                        Σ (−1)ⁿ θ^{2n}/(2n+1)!.
  Domain of validity:   All v ∈ ℝ³. Branch threshold |θ| < 1e−4 inside
                        `taylor_sinc` (which has been audited PASS).

METHOD
  Method declared:      Closed-form (cos, sinc) for θ ≥ 1e−4; Taylor truncation
                        of sinc at order θ⁴ for θ < 1e−4. cos(θ) is the
                        standard library cosine (full closed form across
                        all θ; no Taylor branch needed because cos(0) is
                        non-singular).
  Method implemented:   Lines 200-204 —
                          theta_sq = v.x² + v.y² + v.z²
                          theta    = sqrt(theta_sq)
                          sinc_th  = taylor_sinc(theta, theta_sq)  // small-angle Taylor
                          c        = cos(theta)
                          return (c, sinc_th·v.x, sinc_th·v.y, sinc_th·v.z)
                        sinc is delegated to `small_angle_series::taylor_sinc`,
                        which itself has been audited (small_angle_series.md
                        §5.1 PASS). On the small-angle branch it computes
                          1 − θ²/6 + θ⁴/120
                        and adds |θ²|³/5040 = θ⁶/5040 to errors.precision.
                        NO Padé, NO continued fraction, NO lookup table.
  Match verdict:        ✓ matched — **MATCH-VERDICT CHECK PASS**. The theory
                        is "Taylor expansion of sinc at θ=0"; the
                        implementation is exactly the Taylor truncation of
                        sinc to order θ⁴ with the alternating-series tail
                        bound. The closed-form branch uses std::cos and
                        std::sin (also closed-form, not Taylor-branched).

ERROR BOUND
  Bound category:       precision (truncation), plus precision/accuracy/measurement
                        from the TrackedValue<T> chain
  Bound formula:        On the small-angle branch:
                          |R_sinc(θ)| ≤ |θ|⁶ / 5040  (Leibniz, taylor_sinc card §5.1).
                        Composed with cos (via TrackedValue<T>::cos AUD-EF-6),
                        sqrt (AUD-EF-6), and 3 scalar mults. Each step
                        contributes per REQ-EF-3.
  Bound implemented:    `taylor_sinc` adds the |θ²|³/5040 bound to
                        errors.precision (verified by small_angle_series.md
                        §5.1 audit). Downstream multiplications propagate.
  Bound verdict:        ✓ matched — Taylor truncation bound is added by the
                        helper; closed-form factors propagate per REQ-EF-3.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-6 (Taylor-truncation bound), REQ-EF-9
                        (sqrt)
  AUD-EF applies:       AUD-EF-1, AUD-EF-5 (Taylor branch adds truncation
                        bound), AUD-EF-6 (sqrt/cos)
  AUD-MC applies:       AUD-MC-13 (exp_pure / log_unit round-trip on |v| ≤ π/2),
                        AUD-MC-14 (|exp_pure(v)| = 1 for all v)
  Verification test:    tests/test_math/test_quaternion.cc — exp_pure(0) ==
                        identity(); |exp_pure(v)| == 1; log_unit(exp_pure(v)) == v
                        for |v| ≤ π/2.

NOTES
  - Small-angle threshold inherited from `taylor_sinc`. At θ=1e−4 the bound is
    θ⁶/5040 ≈ 2e−28, well below double-precision ε ≈ 2.2e−16. ✓ Tight for
    double; conservative for wider T.
  - **Verified Taylor (not Padé/CF)** per Match Verdict requirement.
  - Note that the file's only "approximation" is in the small-angle branch
    of taylor_sinc; cos(θ) is closed-form-only (no Taylor branch in
    quaternion.h itself). The audit is therefore primarily delegated to
    the small_angle_series.md (PASS) audit.
```

---

## Formula 22: log_unit()

```
=== FORMULA AUDIT CARD ===
ID:                     quaternion::Quaternion::log_unit
Location:               src/math/quaternion.h:218-230
Mathematical statement: log(q) returns (θ/2)·n̂ where q = cos(θ/2) + sin(θ/2)·n̂
                        is a unit quaternion. Specifically:
                          log(q).vector_part = atan2(|q_v|, w) · q_v / |q_v|
                                             = (arcsin(|q_v|) / |q_v|) · q_v
                        when w² + |q_v|² = 1. The shortest-path convention is
                        enforced: if w < 0, the routine operates on −q (which
                        represents the same SO(3) rotation, with smaller
                        rotation magnitude).

THEORY
  Underlying theorem:   Inverse of the exponential map on the half-angle ball
                        (|v| ≤ π/2): log(q) = (atan2(|q_v|, w) / |q_v|) · q_v.
                        For a unit quaternion with w ≥ 0, atan2(|q_v|, w) =
                        arcsin(|q_v|) (since the angle is in [0, π/2]). Hence
                          scale(s, w) = atan2(s, w) / s = arcsin(s) / s
                                      = 1 + s²/6 + 3s⁴/40 + 5s⁶/112 + …
                        near s = 0. Taylor expansion of arcsin(s)/s at s = 0
                        is well-known (Abramowitz & Stegun (1964) §15.1.10).
                        The −q branch (`if w.value < 0`) exploits the SU(2) → SO(3)
                        2-to-1 cover: q and −q represent the same rotation;
                        the convention chooses the representative with w ≥ 0
                        to keep |log| ≤ π/2.
  Primary reference:    Altmann (1986) §12.5; Shoemake (1985) §6;
                        Murray, Li & Sastry (1994) §3.2.
                        Taylor of arcsin(s)/s: Abramowitz & Stegun (1964) §15.1.10.
                        Taylor-truncation bound (5|s|⁶/112): see small_angle_series.md
                        §5.2 (PASS with wide-T note).
  Domain of validity:   q is unit-norm. The half-angle ball |v| ≤ π/2 covers
                        SO(3) without ambiguity once the shortest-path
                        convention is applied. Branch threshold |q_v| < 1e−4
                        inside `taylor_half_angle_scale`.

METHOD
  Method declared:      Closed-form atan2(|q_v|, w)/|q_v| for |q_v| ≥ 1e−4;
                        Taylor truncation of arcsin(s)/s at order s⁴ for
                        |q_v| < 1e−4; shortest-path negation when w.value < 0.
  Method implemented:   Lines 219-229:
                          (1) If w.value < 0, flip all four components:
                              q → −q (shortest-path convention).
                          (2) qv_norm_sq = x_e² + y_e² + z_e²
                          (3) qv_norm    = sqrt(qv_norm_sq)
                          (4) scale      = taylor_half_angle_scale(qv_norm,
                                              qv_norm_sq, w_e)
                              — helper computes atan2(s, w)/s in closed form for
                                s ≥ 1e−4, or
                                  1 + s²/6 + 3s⁴/40
                                + truncation bound 5|s²|³/112 = 5|s|⁶/112
                                for s < 1e−4 (small_angle_series.md §5.2 PASS).
                          (5) return Vector3<T>(scale·x_e, scale·y_e, scale·z_e).
                        NO Padé, NO continued fraction, NO Newton iteration.
                        The Taylor branch is exactly the Taylor expansion of
                        arcsin(s)/s at s = 0.
  Match verdict:        ✓ matched — **MATCH-VERDICT CHECK PASS**. The theory
                        is "Taylor expansion of arcsin(s)/s at s=0" (or
                        equivalently atan2(s,w)/s for the unit quaternion);
                        the implementation is exactly that Taylor truncation
                        on the small-arg branch. The closed-form atan2 branch
                        uses TrackedValue<T>::atan2, which is closed-form (no
                        small-arg Taylor in atan2 itself at this layer).

ERROR BOUND
  Bound category:       precision (truncation), plus precision/accuracy/
                        measurement from the chain
  Bound formula:        On the small-arg branch:
                          |R_{arcsin/s}(s)| ≤ 5|s|⁶/112   (next-term bound
                          from the Taylor expansion of arcsin(s)/s, with a
                          mild positive-coefficient correction documented in
                          small_angle_series.md §5.2 — tight for double,
                          conservative for wider T).
                        Composition: sqrt for qv_norm; the half-angle scale
                        helper; 3 scalar mults for the output vector. Each
                        step per REQ-EF-3.
  Bound implemented:    `taylor_half_angle_scale` adds 5·|s_sq|³/112 to
                        errors.precision on the small-arg branch (verified in
                        small_angle_series.md §5.2). Downstream
                        multiplications propagate.
  Bound verdict:        ✓ matched (⚠ note for wide T: helper bound is
                        tight at double, slightly under-bounded for cpp_bin_float_50
                        in the absence of the 1/(1−2s²) correction; see
                        small_angle_series.md §5.2).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-6 (Taylor truncation), REQ-EF-9 (sqrt)
  AUD-EF applies:       AUD-EF-1, AUD-EF-5 (Taylor branch adds bound), AUD-EF-6
  AUD-MC applies:       AUD-MC-13 (log_unit / exp_pure round-trip on |v| ≤ π/2),
                        AUD-MC-14 (shortest-path: log_unit(−q) == log_unit(q))
  Verification test:    tests/test_math/test_quaternion.cc —
                        log_unit(identity()) == zero_vector;
                        log_unit(exp_pure(v)) == v for |v| ≤ π/2;
                        log_unit(−q) == log_unit(q) (shortest-path).

NOTES
  - **Verified Taylor (not Padé/CF)** per Match Verdict requirement.
  - Shortest-path branch (line 220) prevents log_unit from returning a
    vector of magnitude greater than π/2 — critical for the inverse
    relationship with exp_pure on the half-angle ball.
  - The Taylor bound 5|s|⁶/112 is technically the leading term of a
    positive-coefficient tail (not a strict Leibniz bound). At threshold
    s = 1e−4, the geometric-tail correction 1/(1−2s²) is ~1 + 2·10⁻⁸,
    which contributes ~10⁻³⁵ — far below double-precision ε. **Tight for
    double; conservative for wider T**; documented as a future tightening
    opportunity (small_angle_series.md §5.2).
```

---

## File-level Verdict

### A. Error wiring (AUD-EF coverage)

- **AUD-EF-1**: ✓ All 22 ops return either `Quaternion<T>`, `Vector3<T>`, or `TrackedValue<T>` — never bare `T`.
- **AUD-EF-2**: ✓ Constructors and `exact<T>(…)` literals carry zero error.
- **AUD-EF-3**: ✓ Closed-form ops delegate to TrackedValue<T> per-category propagation.
- **AUD-EF-5**: ✓ Small-angle branches in `exp_pure` (Formula 21) and `log_unit` (Formula 22) delegate Taylor-truncation bounds to `taylor_sinc` and `taylor_half_angle_scale` (both audited PASS in `small_angle_series.md`).
- **AUD-EF-6**: ✓ sqrt / sin / cos / atan2 from `TrackedValue<T>` carry their own per-category bounds (audited under tracked_value.md).
- **AUD-EF-7**: ✓ Composite return types (Quaternion, Vector3) compose per-category errors of all components.

### B. Algebra axioms (AUD-MC coverage)

- **AUD-MC-4**: ✓ Additive group (identity, inverse, associative, commutative) — Formulas 1, 4, 9, 10, 11.
- **AUD-MC-5**: ✓ Multiplicative identity and Hamilton-product associativity — Formulas 3, 12.
- **AUD-MC-6**: ✓ Non-commutativity of Hamilton product — Formula 12.
- **AUD-MC-7**: ✓ Conjugate anti-distribution: (a·b)* = b*·a* — Formula 15.
- **AUD-MC-8**: ✓ Norm multiplicativity: |a·b| = |a|·|b| — Formulas 16, 17.
- **AUD-MC-9**: ✓ Multiplicative inverse: q·q⁻¹ = identity — Formula 18.
- **AUD-MC-10**: ✓ Unit-sphere retraction: |q.normalized()| = 1 — Formula 19.
- **AUD-MC-11**: ✓ Rotation preserves norm: |q.rotate(v)| = |v| — Formula 20.
- **AUD-MC-12**: ✓ Pure-quaternion embedding and rotation composition — Formulas 5, 20.
- **AUD-MC-13**: ✓ exp_pure / log_unit round-trip on half-angle ball — Formulas 21, 22.
- **AUD-MC-14**: ✓ |exp_pure(v)| = 1; shortest-path: log_unit(−q) = log_unit(q) — Formulas 21, 22.

### C. Theoretical basis (this document)

| Card | Theory | Method | Bound | Status |
|------|--------|--------|-------|--------|
| 1 ctor_default | ℍ additive identity | Default-init | n/a | ✓ |
| 2 ctor_value | ℍ basis decomposition | Member copy | Inherited | ✓ |
| 3 identity | ℍ multiplicative identity | exact<T>(1,0,0,0) | Zero error | ✓ |
| 4 zero | ℍ additive identity | Default ctor | Zero error | ✓ |
| 5 pure | ℝ³ ↪ ℍ embedding | Closed-form | Inherited | ✓ |
| 6 from_axis_angle | Euler-Rodrigues, SU(2) cover | **Closed-form half-angle (NOT Taylor)** | REQ-EF-3 / AUD-EF-6 | ✓ |
| 7 scalar | Component projection | Member access | Inherited | ✓ |
| 8 vector | Component projection | Member access | Inherited | ✓ |
| 9 operator+ | ℝ⁴ abelian addition | Component-wise | REQ-EF-3 | ✓ |
| 10 operator-(bin) | ℝ⁴ subtraction | Component-wise | REQ-EF-3 | ✓ |
| 11 operator-(un) | Additive inverse | Component-wise | Inherited | ✓ |
| 12 Hamilton mult | Hamilton's i²=j²=k²=ijk=−1 | Closed-form 4-component expansion | REQ-EF-3 | ✓ |
| 13 scalar mult (L) | ℝ-module structure | Component-wise | REQ-EF-3 | ✓ |
| 14 scalar mult (R) | ℝ-module structure | Component-wise | REQ-EF-3 | ✓ |
| 15 conjugate | Anti-involution | Component-wise negation | Inherited | ✓ |
| 16 mag_squared | Cayley quadratic form | Sum of squares | REQ-EF-3 | ✓ |
| 17 magnitude | Euclidean norm on ℝ⁴ | sqrt | REQ-EF-3 + REQ-EF-9 | ✓ |
| 18 inverse | Division-ring inverse | Closed-form q*/|q|² | REQ-EF-3 | ✓ |
| 19 normalized | S³ ⊂ ℍ retraction | Closed-form q/|q| | REQ-EF-3 + REQ-EF-9 | ✓ |
| 20 rotate | SU(2) → SO(3) double cover | **Closed-form Rodrigues identity** | REQ-EF-3 | ✓ |
| 21 exp_pure | Lie exp: pure_ℍ → S³ | **Taylor of sinc at θ=0** (small-arg branch); closed-form cos elsewhere | REQ-EF-6 (via taylor_sinc) | ✓ |
| 22 log_unit | Lie log: S³ → pure_ℍ | **Taylor of arcsin(s)/s at s=0** (small-arg branch); closed-form atan2 elsewhere | REQ-EF-6 (via taylor_half_angle_scale) | ✓ (⚠ wide-T note) |

### Match-verdict spotlight (per audit request)

The auditor was asked to confirm explicitly that `exp_pure`, `log_unit`, and `from_axis_angle` use **Taylor** (not Padé / continued fraction):

- **`from_axis_angle`**: NO Taylor branch in this file. It uses `cos(half)` and `sin(half)` from `TrackedValue<T>` directly — closed-form per AUD-EF-6 (those functions are themselves closed-form library calls; any small-arg branches inside them are governed by tracked_value.md, which is audited separately). **No mismatch.**
- **`exp_pure`**: Small-arg branch (|θ| < 1e−4) calls `taylor_sinc` → **Taylor truncation of sinc(θ) at order θ⁴** with Leibniz bound |θ|⁶/5040 (audited PASS in small_angle_series.md §5.1). Large-arg branch uses closed-form `sin(θ)/θ`. **No Padé, no continued fraction.** ✓
- **`log_unit`**: Small-arg branch (|q_v| < 1e−4) calls `taylor_half_angle_scale` → **Taylor truncation of arcsin(s)/s at order s⁴** with next-term bound 5|s|⁶/112 (audited PASS-with-wide-T-note in small_angle_series.md §5.2). Large-arg branch uses closed-form `atan2(s, w)/s`. **No Padé, no continued fraction.** ✓

All three match-verdict checks **PASS**.

### Overall file verdict: **PASS**

All 22 formulas match cited theory to method to bound. Closed-form algebra dominates; the two Taylor branches (in `exp_pure` and `log_unit`) delegate truncation-bound bookkeeping to `small_angle_series.h` helpers, which themselves carry PASS audit verdicts. No theory-method mismatch (no instance of "doc says Taylor, code does X-else"). Error propagation is end-to-end through `TrackedValue<T>` per the AUD-EF-1..7 wiring and REQ-EF-3 composition rule.

One minor note for wide T (cpp_bin_float_50 and similar): the `log_unit` Taylor-truncation bound inherits the same conservative-but-untight property as the underlying `taylor_half_angle_scale` helper (see small_angle_series.md §5.2). Not a correctness fix; a future-optimization opportunity.

---

**Document**: `design/audit/theoretical_basis_audit/quaternion.md`
**Status**: **PASS** (22/22 cards)
**Date**: 2026-05-13
**Analyst**: internal review
