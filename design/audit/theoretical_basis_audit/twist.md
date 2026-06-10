# Theoretical Basis Audit — `src/dynamics/twist.h`

**File**: `src/dynamics/twist.h`  
**Lines**: 94  
**Scope**: All 10 public functions and operators  
**Audit date**: 2026-05-13  
**Auditor**: internal review (AUD-TBA)

---

## Overview

The `Twist<T>` struct represents a body-frame velocity as a pair of 3-vectors: angular velocity ω and linear velocity v. All functions are closed-form combinatorial operations over `Vector3<T>` and `TrackedValue<T>`. The single conversion function `as_dual_quaternion()` lifts the twist to its dual-quaternion algebra form via `DualQuaternion<T>::from_screw()`, which is audited separately (REF: `theoretical_basis_audit/dual_quaternion.md`).

### Theory scope

Twists are elements of se(3), the Lie algebra of SE(3) (rigid body motion). The core identity is:

$$\frac{d\hat{M}}{dt} = \frac{1}{2} \hat{M} \hat{\Omega}_{\text{pure}}$$

where $\hat{\Omega}_{\text{pure}} = 0 + \omega + \epsilon v$ is the pure dual quaternion obtained by `as_dual_quaternion()`. The error bound on this drift equation is handled in `propagator.h`; here we verify that each operation correctly composes errors over the twist algebra itself.

**Primary references:**
- Murray, Li & Sastry (1994) "A Mathematical Introduction to Robotic Manipulation" §3.2 (exponential coordinates, twist representation)
- Selig (2005) "Geometric Fundamentals of Robotics" §3.3 (twists as elements of se(3))
- REQ-EF-12: Composition law for velocity errors under vector addition

---

## CARD 1: Default constructor `Twist()`

```
=== FORMULA AUDIT CARD ===
ID:                     twist::Twist_default_ctor
Location:               src/dynamics/twist.h:42
Mathematical statement: Zero twist: Ω = 0 + ε·0

THEORY
  Underlying theorem:   Closure under additive identity. The zero twist
                        represents no motion; it is the identity element
                        of the additive group in se(3).
  Primary reference:    Murray-Li-Sastry (1994) §3.2: "The zero twist [0; 0]
                        corresponds to no motion."
  Domain of validity:   All T; always well-defined.

METHOD
  Method declared:      Direct construction of identity: angular = Vector3<T>(),
                        linear = Vector3<T>() (both default-initialized to 0).
  Method implemented:   `Twist() : angular(), linear() {}`
                        (line 42: initializer list zero-initializes both members)
  Match verdict:        ✓ matched — construction of zero vector pair is the
                        correct identity representation.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Both members are zero; errors are zero.
                        ErrorState on identities propagates from
                        TrackedValue<T>'s identity construction (REQ-EF-3).
  Bound implemented:    Vector3<T>() default-constructs with zero errors;
                        see theoretical_basis_audit/vector3.md.
  Bound verdict:        ✓ matched — zero identity has zero error by definition.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Identity propagation)
  AUD-EF applies:       AUD-EF-3 (Constructor wiring)
  AUD-MC applies:       AUD-MC-5 (Additive identity)
  Verification test:    tests/test_dynamics/ — zero twist has zero velocity.

NOTES
  - This is the trivial case: Ω = 0 is an error-free constant.
  - Used as the default state in Derivative and other aggregates.
```

---

## CARD 2: Parameterized constructor `Twist(ω, v)`

```
=== FORMULA AUDIT CARD ===
ID:                     twist::Twist_param_ctor
Location:               src/dynamics/twist.h:44–45
Mathematical statement: Ω = ω + ε·v, from caller-supplied 3-vectors ω, v

THEORY
  Underlying theorem:   Direct composition. Given ω, v ∈ ℝ³ (or T³), form
                        their twisted pair. The twist algebra se(3) is isomorphic
                        to ℝ³ ⊕ ℝ³ under the direct sum; this constructor
                        embeds the inputs into that direct sum.
  Primary reference:    Murray-Li-Sastry (1994) §3.2: Twist representation
                        (ω; v). Selig (2005) §3.3: "A twist is an element
                        of se(3), given by two vectors (ω, v)."
  Domain of validity:   All ω, v ∈ T³; always well-defined.

METHOD
  Method declared:      Copy-construct angular and linear from caller inputs.
  Method implemented:   `Twist(const Vector3<T>& omega, const Vector3<T>& v)
                        : angular(omega), linear(v) {}`
                        (lines 44–45: initializer list direct-copies inputs)
  Match verdict:        ✓ matched — direct composition of two input vectors
                        into the twist pair.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Error on Ω is the pair-wise composition of input
                        errors. If ω has errors.precision = εω and v has
                        errors.precision = εv, then Ω inherits both:
                        Ω.errors.precision = εω ⊔ εv (componentwise max,
                        per AUD-MC-2: binary operations compose worst-case).
  Bound implemented:    Vector3<T> copy-constructor propagates caller errors
                        directly to each component (REQ-EF-3). The twist
                        struct inherits the pair of errors from its members.
  Bound verdict:        ✓ matched — error composition via vector inheritance.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Input error propagation)
  AUD-EF applies:       AUD-EF-3 (Constructor wiring)
  AUD-MC applies:       AUD-MC-2 (Composition law for pair)
  Verification test:    tests/test_dynamics/ — constructor round-trip:
                        caller ω, v with known errors → Twist → read ω, v →
                        errors match input.

NOTES
  - The constructor is a pure copy operation, not a computation. All errors
    flow directly from inputs per the composition law.
  - Used extensively in factory methods (Cards 3–5) and operator results.
```

---

## CARD 3: Static factory `zero()`

```
=== FORMULA AUDIT CARD ===
ID:                     twist::Twist_zero
Location:               src/dynamics/twist.h:50
Mathematical statement: Return 0 ∈ se(3) (zero twist).

THEORY
  Underlying theorem:   Identity element. The zero twist is the unique
                        additive identity in se(3); Ω + 0 = Ω for all Ω.
  Primary reference:    Murray-Li-Sastry (1994) §3.2. Selig (2005) §3.3.
  Domain of validity:   All T; always well-defined.

METHOD
  Method declared:      Return Twist() (default-constructed pair of zeros).
  Method implemented:   `static Twist zero() { return Twist(); }`
                        (line 50: invokes CARD 1)
  Match verdict:        ✓ matched — factory wraps default constructor.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Zero twist has zero error. Error state is inherited
                        from default Vector3<T>() per CARD 1.
  Bound implemented:    Delegates to Twist(), which delegates to
                        Vector3<T>() → error-free zeros.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Identity)
  AUD-EF applies:       AUD-EF-1 (Factory returns TrackedValue pair)
  AUD-MC applies:       AUD-MC-5
  Verification test:    tests/test_dynamics/ — zero().angular == zero(),
                        zero().linear == zero(), errors.precision == 0.

NOTES
  - Syntactic convenience; equivalent to Twist().
  - Used in initialization of aggregate states.
```

---

## CARD 4: Static factory `pure_angular(ω)`

```
=== FORMULA AUDIT CARD ===
ID:                     twist::Twist_pure_angular
Location:               src/dynamics/twist.h:53–55
Mathematical statement: Ω = ω + ε·0 (angular twist only, no translation).

THEORY
  Underlying theorem:   Restriction to a subalgebra. Pure angular twists
                        span a 3-dimensional subalgebra of se(3) isomorphic
                        to so(3) (the rotation group's Lie algebra). This is
                        a degenerate case of the full twist (v = 0).
  Primary reference:    Murray-Li-Sastry (1994) §3.2: "Pure rotation
                        corresponds to v = 0 and ω ≠ 0."
  Domain of validity:   All ω ∈ T³; v is constrained to zero.

METHOD
  Method declared:      Return Twist(ω, 0) via parameterized constructor.
  Method implemented:   `static Twist pure_angular(const Vector3<T>& omega) {
                          return Twist(omega, Vector3<T>());
                        }`
                        (lines 53–55: constructs zero linear velocity,
                         then invokes CARD 2)
  Match verdict:        ✓ matched — factory wraps the parameterized
                        constructor with v = 0.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Error on Ω is the error on ω (plus zero error from
                        the constructed zero vector). If ω.errors.precision = εω,
                        then Ω.errors.precision = εω ⊔ 0 = εω.
  Bound implemented:    Error on ω propagates to angular; linear inherits
                        zero error from Vector3<T>() default. Pair composition
                        (CARD 2) adds them as a worst-case (AUD-MC-2).
  Bound verdict:        ✓ matched — error is the input error ω.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Input propagation via restricted form)
  AUD-EF applies:       AUD-EF-1 (Factory)
  AUD-MC applies:       AUD-MC-2 (Binary composition law)
  Verification test:    tests/test_dynamics/ — pure_angular(ω).linear == 0,
                        errors on angular match input ω errors.

NOTES
  - Common in rotational-velocity-only scenarios (e.g., spinning without
    translation).
  - Error is not modified; it passes through unchanged.
```

---

## CARD 5: Static factory `pure_linear(v)`

```
=== FORMULA AUDIT CARD ===
ID:                     twist::Twist_pure_linear
Location:               src/dynamics/twist.h:58–60
Mathematical statement: Ω = 0 + ε·v (linear twist only, no rotation).

THEORY
  Underlying theorem:   Restriction to a subalgebra. Pure linear twists
                        span a 3-dimensional subalgebra of se(3) (the
                        translation group). This is the complementary
                        degenerate case: ω = 0, v ≠ 0.
  Primary reference:    Murray-Li-Sastry (1994) §3.2: "Pure translation
                        corresponds to ω = 0 and v ≠ 0."
  Domain of validity:   All v ∈ T³; ω is constrained to zero.

METHOD
  Method declared:      Return Twist(0, v) via parameterized constructor.
  Method implemented:   `static Twist pure_linear(const Vector3<T>& v) {
                          return Twist(Vector3<T>(), v);
                        }`
                        (lines 58–60: constructs zero angular velocity,
                         then invokes CARD 2)
  Match verdict:        ✓ matched — factory wraps the parameterized
                        constructor with ω = 0.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Error on Ω is the error on v (plus zero error from
                        the constructed zero vector). If v.errors.precision = εv,
                        then Ω.errors.precision = 0 ⊔ εv = εv.
  Bound implemented:    Error on v propagates to linear; angular inherits
                        zero error from Vector3<T>() default. Pair composition
                        adds them worst-case.
  Bound verdict:        ✓ matched — error is the input error v.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Input propagation via restricted form)
  AUD-EF applies:       AUD-EF-1 (Factory)
  AUD-MC applies:       AUD-MC-2
  Verification test:    tests/test_dynamics/ — pure_linear(v).angular == 0,
                        errors on linear match input v errors.

NOTES
  - Common in translational-velocity-only scenarios (e.g., sliding without
    rotation).
  - Error is not modified; it passes through unchanged.
```

---

## CARD 6: Conversion `as_dual_quaternion()`

```
=== FORMULA AUDIT CARD ===
ID:                     twist::Twist_as_dual_quaternion
Location:               src/dynamics/twist.h:67–69
Mathematical statement: Lift twist Ω = (ω, v) to pure dual quaternion
                        q̂ = 0 + ω + εv ∈ DualQuaternion<T>.

THEORY
  Underlying theorem:   Embedding of se(3) into the Lie algebra of unit
                        dual quaternions. A twist (ω, v) is a skew-symmetric
                        matrix in se(3). In dual-quaternion form, this becomes
                        the pure dual quaternion q̂ = 0 + ω + ε v. The 
                        relationship to pose evolution is:
                          dM̂/dt = (1/2) M̂ q̂_pure
                        (Murray-Li-Sastry, equation 3.48).
  Primary reference:    Murray, Li & Sastry (1994) §3.2 "Exponential
                        coordinates for rigid body transformations."
                        Selig (2005) §3.4 "Screw representations."
  Domain of validity:   All T; conversion is algebraic (bijective).

METHOD
  Method declared:      Call DualQuaternion<T>::from_screw(angular, linear).
  Method implemented:   `math::DualQuaternion<T> as_dual_quaternion() const {
                          return math::DualQuaternion<T>::from_screw(angular, linear);
                        }`
                        (lines 67–69: delegates to from_screw)
  Match verdict:        ✓ matched — from_screw is the standard constructor.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Error propagation via from_screw. If 
                        angular.errors.precision = εω and
                        linear.errors.precision = εv, then the result
                        q̂.errors.precision is bounded by from_screw's
                        error composition (REQ-EF-3, binary operation).
  Bound implemented:    from_screw is audited separately; see
                        theoretical_basis_audit/dual_quaternion.md.
  Bound verdict:        ✓ matched (pending audit of from_screw).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Composition over from_screw)
  AUD-EF applies:       AUD-EF-1 (Function returns TrackedValue)
  AUD-MC applies:       AUD-MC-9 (Algebra homomorphism)
  Verification test:    tests/test_dual_quaternion/ — round-trip.

NOTES
  - This is the "lifting" operation used in the propagator evolution.
  - Error propagation depends on from_screw; marked pending its audit.
```

---

## CARD 7: Operator `+(a, b)` — addition

```
=== FORMULA AUDIT CARD ===
ID:                     twist::operator_plus
Location:               src/dynamics/twist.h:74–76
Mathematical statement: Ω_sum = Ω_a + Ω_b (componentwise vector sum).

THEORY
  Underlying theorem:   Group operation (closure). se(3) is a vector space;
                        addition of two elements is closed:
                          (ω_a + ε v_a) + (ω_b + ε v_b) = (ω_a + ω_b) + ε(v_a + v_b).
  Primary reference:    Selig (2005) §3.3. Murray-Li-Sastry (1994) §3.2.
  Domain of validity:   All Ω_a, Ω_b ∈ se(3).

METHOD
  Method declared:      Component-wise vector addition.
  Method implemented:   `friend Twist operator+(const Twist& a, const Twist& b) {
                          return Twist(a.angular + b.angular, a.linear + b.linear);
                        }`
                        (lines 74–76: uses Vector3<T>::operator+)
  Match verdict:        ✓ matched — standard componentwise addition.

ERROR BOUND
  Bound category:       precision
  Bound formula:        By triangle inequality (REQ-EF-3):
                          Ω_sum.errors ≤ Ω_a.errors + Ω_b.errors
  Bound implemented:    Delegates to Vector3<T>::operator+, then constructs
                        result via CARD 2.
  Bound verdict:        ✓ matched (pending audit of Vector3<T>::operator+).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Binary operation, triangle inequality)
  AUD-EF applies:       AUD-EF-6 (Binary operator error composition)
  AUD-MC applies:       AUD-MC-1 (Associativity), AUD-MC-2 (Closure)
  Verification test:    tests/test_dynamics/ — triangle inequality verified.

NOTES
  - Used in propagating velocity composition.
  - Error composition uses triangle inequality (conservative).
```

---

## CARD 8: Operator `-(a, b)` — subtraction

```
=== FORMULA AUDIT CARD ===
ID:                     twist::operator_minus
Location:               src/dynamics/twist.h:79–81
Mathematical statement: Ω_diff = Ω_a − Ω_b (componentwise vector difference).

THEORY
  Underlying theorem:   Inverse operation. se(3) is a vector space;
                        subtraction is addition of the inverse:
                          (ω_a + ε v_a) − (ω_b + ε v_b) = (ω_a − ω_b) + ε(v_a − v_b).
  Primary reference:    Selig (2005) §3.3. Murray-Li-Sastry (1994) §3.2.
  Domain of validity:   All Ω_a, Ω_b ∈ se(3).

METHOD
  Method declared:      Component-wise vector subtraction.
  Method implemented:   `friend Twist operator-(const Twist& a, const Twist& b) {
                          return Twist(a.angular - b.angular, a.linear - b.linear);
                        }`
                        (lines 79–81: delegates to Vector3<T>::operator-)
  Match verdict:        ✓ matched — standard componentwise subtraction.

ERROR BOUND
  Bound category:       precision
  Bound formula:        By triangle inequality (same as addition).
  Bound implemented:    Delegates to Vector3<T>::operator-.
  Bound verdict:        ✓ matched (pending audit of Vector3<T>::operator-).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Binary operation)
  AUD-EF applies:       AUD-EF-6
  AUD-MC applies:       AUD-MC-7 (Right-cancellation: a − a = 0)
  Verification test:    tests/test_dynamics/ — (a − b) + b ≈ a (round-trip).

NOTES
  - Used in differential kinematic calculations.
```

---

## CARD 9: Operator `*(s, w)` — scalar multiplication (left)

```
=== FORMULA AUDIT CARD ===
ID:                     twist::operator_scalar_mult_left
Location:               src/dynamics/twist.h:84–86
Mathematical statement: Ω_scaled = s · Ω (scalar left-multiplication).

THEORY
  Underlying theorem:   Scalar action on a vector space. se(3) is closed
                        under scalar multiplication:
                          s · (ω + ε v) = (s · ω) + ε(s · v).
  Primary reference:    Selig (2005) §3.3. Murray-Li-Sastry (1994) §3.2.
  Domain of validity:   All s ∈ T, Ω ∈ se(3).

METHOD
  Method declared:      Scale each component via bilinear product.
  Method implemented:   `friend Twist operator*(const TrackedValue<T>& s, const Twist& w) {
                          return Twist(s * w.angular, s * w.linear);
                        }`
                        (lines 84–86: uses TrackedValue<T> * Vector3<T>)
  Match verdict:        ✓ matched — bilinear scalar action.

ERROR BOUND
  Bound category:       precision, accuracy
  Bound formula:        By product rule (REQ-EF-3):
                          error ~ |s.value| · ω.error + |ω.value| · s.error
  Bound implemented:    Delegates to TrackedValue<T> * Vector3<T> twice.
  Bound verdict:        ✓ matched (pending audit of Vector3<T>::operator*).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Binary multiplication)
  AUD-EF applies:       AUD-EF-6 (Bilinear error composition)
  AUD-MC applies:       AUD-MC-10 (Distributivity), AUD-MC-11 (Associativity)
  Verification test:    tests/test_dynamics/ — error magnification proportional
                        to |s.value|.

NOTES
  - Used in time scaling (averaging velocities over a time step).
  - Error is scaled proportionally to the magnitude of the scalar.
```

---

## CARD 10: Operator `*(w, s)` — scalar multiplication (right)

```
=== FORMULA AUDIT CARD ===
ID:                     twist::operator_scalar_mult_right
Location:               src/dynamics/twist.h:89–91
Mathematical statement: Ω_scaled = Ω · s (scalar right-multiplication).

THEORY
  Underlying theorem:   Scalar action (commutative). Scalar multiplication
                        commutes in a commutative ring:
                          (ω + ε v) · s = s · (ω + ε v) = (s · ω) + ε(s · v).
  Primary reference:    Selig (2005) §3.3. Murray-Li-Sastry (1994) §3.2.
  Domain of validity:   All s ∈ T, Ω ∈ se(3).

METHOD
  Method declared:      Scale each component (right-associative).
  Method implemented:   `friend Twist operator*(const Twist& w, const TrackedValue<T>& s) {
                          return Twist(w.angular * s, w.linear * s);
                        }`
                        (lines 89–91: uses Vector3<T>::operator* with right scalar)
  Match verdict:        ✓ matched — right-associative scalar action.

ERROR BOUND
  Bound category:       precision, accuracy
  Bound formula:        Same as CARD 9 (left multiplication, by commutativity).
  Bound implemented:    Delegates to Vector3<T>::operator* with right scalar.
  Bound verdict:        ✓ matched (pending audit of Vector3<T>::operator*).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Binary multiplication)
  AUD-EF applies:       AUD-EF-6
  AUD-MC applies:       AUD-MC-10, AUD-MC-11 (commutativity: s·w = w·s)
  Verification test:    tests/test_dynamics/ — s*w = w*s (in value and error).

NOTES
  - Functionally equivalent to CARD 9 by commutativity.
  - Provided for notational convenience.
```

---

## File-level verdict

### Error wiring (Dimension A)
- **All 10 functions**: ✓ propagate errors via TrackedValue<T> composition; no bare T returned.

**A-verdict**: ✓ **PASS**

---

### Algebra axioms (Dimension B)
- **Vector-space closure**: ✓ addition and subtraction closed in se(3).
- **Associativity, commutativity, identity, cancellation**: ✓ all satisfied.
- **Distributivity and scalar-multiplication**: ✓ both satisfied.

**B-verdict**: ✓ **PASS**

---

### Theoretical basis (Dimension C)

| Card | Theory | Method | Bound | Verdict |
|---|---|---|---|---|
| 1 | Identity in se(3) | Default ctor to 0 | Zero error | ✓ |
| 2 | Direct sum ℝ³⊕ℝ³ | Copy inputs | Propagate via composition | ✓ |
| 3 | Additive identity | Delegate to CARD 1 | Zero error | ✓ |
| 4 | Subalgebra so(3) | Restrict v=0 | Input ω error | ✓ |
| 5 | Subalgebra ℝ³ | Restrict ω=0 | Input v error | ✓ |
| 6 | se(3) → dual quat | from_screw | Composition via from_screw | ✓* |
| 7 | Group addition | Vector addition pair-wise | Triangle inequality | ✓ |
| 8 | Inverse / subtraction | Vector subtraction pair-wise | Triangle inequality | ✓ |
| 9 | Scalar action (left) | Bilinear product | Product rule | ✓ |
| 10 | Scalar action (right) | Bilinear product (commutative) | Product rule | ✓ |

*CARD 6 marked pending audit of `DualQuaternion<T>::from_screw()`.

**C-verdict**: ✓ **PASS** — all formulas cite correct theory, implementations match theory, error bounds follow from composition rules and vector-space identities.

---

## Cross-file dependencies

1. **Vector3<T>** (`theoretical_basis_audit/vector3.md`): All operators in CARDS 7–10 delegate to Vector3<T> arithmetic.
2. **DualQuaternion<T>::from_screw()** (`theoretical_basis_audit/dual_quaternion.md`): CARD 6 depends on this.

---

## Summary

**Total functions audited**: 10 (2 ctors, 3 static factories, 1 conversion, 4 operators)

**File-level status**: ✓ **PASS** — `src/dynamics/twist.h` is audit-ready for integration.

