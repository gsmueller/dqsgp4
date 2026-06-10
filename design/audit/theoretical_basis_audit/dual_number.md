# Theoretical Basis Audit — `src/math/dual_number.h`

**File**: `src/math/dual_number.h` (141 lines)  
**Scope**: DualNumber struct, 16 distinct formulas  
**Theory anchor**: Clifford (1873), forward-mode automatic differentiation (Wengert 1964), SE(3) screw calculus  
**Key question**: Line 113 — is "Taylor expansion" the actual method, or misdeclaration?

---

## Formula 1: Default constructor

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::DualNumber::ctor_default
Location:               src/math/dual_number.h:44
Mathematical statement: DualNumber() constructs (0 + ε·0) with both parts uninitialized (defaulted TrackedValue<T>)

THEORY
  Underlying theorem:   None — this is a trivial initialization, equivalent to the
                        default-constructed values of two TrackedValue<T> objects.
  Primary reference:    C++ struct default construction semantics.
  Domain of validity:   All T.

METHOD
  Method declared:      Default constructor delegates to default constructors of
                        real and dual members.
  Method implemented:   Trivial delegation. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (no computation; error state inherited from member defaults)
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A (no numeric operation)
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    tests/test_math/ — constructor can be called.

NOTES
  - No arithmetic is performed; no error propagation.
```

---

## Formula 2: Constructor from TrackedValue<T>

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::DualNumber::ctor_real_only
Location:               src/math/dual_number.h:47
Mathematical statement: DualNumber(const TrackedValue<T>& r) promotes a real number to a dual with zero dual part

THEORY
  Underlying theorem:   None — promotion is identity in the dual extension.
  Primary reference:    N/A
  Domain of validity:   All T.

METHOD
  Method declared:      Direct assignment of r to real; default construction (zero) for dual.
  Method implemented:   `real(r), dual()` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    tests/test_math/ — promotion and round-trip access.

NOTES
  - Error state inherited from input r.errors and default-constructed dual.errors.
```

---

## Formula 3: Constructor from real and dual parts

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::DualNumber::ctor_bipart
Location:               src/math/dual_number.h:50–51
Mathematical statement: DualNumber(const TrackedValue<T>& r, const TrackedValue<T>& d) constructs a + ε b directly

THEORY
  Underlying theorem:   None — direct construction.
  Primary reference:    N/A
  Domain of validity:   All T.

METHOD
  Method declared:      Direct assignment of r to real, d to dual.
  Method implemented:   `real(r), dual(d)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       N/A
  Verification test:    tests/test_math/ — constructor and accessors.

NOTES
  - Error state inherited from inputs.
```

---

## Formula 4: Factory zero()

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::DualNumber::zero
Location:               src/math/dual_number.h:56
Mathematical statement: zero() returns the additive identity (0 + ε·0)

THEORY
  Underlying theorem:   Definition of additive identity in a ring.
  Primary reference:    N/A
  Domain of validity:   All T.

METHOD
  Method declared:      Return DualNumber() — default-constructed.
  Method implemented:   Trivial delegation.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       AUD-MC-1 (additive identity in ring)
  Verification test:    tests/test_math/ — zero() == zero() and a + zero() == a

NOTES
  - Satisfies AUD-MC-1 requirement.
```

---

## Formula 5: Factory identity()

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::DualNumber::identity
Location:               src/math/dual_number.h:59–61
Mathematical statement: identity() returns the multiplicative identity (1 + ε·0)

THEORY
  Underlying theorem:   Definition of multiplicative identity in a ring.
  Primary reference:    N/A
  Domain of validity:   All T.

METHOD
  Method declared:      Return DualNumber(exact<T>(1), exact<T>(0)).
  Method implemented:   Correct — constructs TrackedValue<T>(1, 0) in real, zero in dual.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       AUD-MC-2 (multiplicative identity)
  Verification test:    tests/test_math/ — identity() * a == a for all a

NOTES
  - Satisfies AUD-MC-2 requirement.
```

---

## Formula 6: Factory epsilon()

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::DualNumber::epsilon
Location:               src/math/dual_number.h:64–66
Mathematical statement: epsilon() returns the nilpotent unit (0 + ε·1)

THEORY
  Underlying theorem:   Definition of the free generator ε in the dual-number algebra,
                        satisfying ε² = 0.
  Primary reference:    Clifford (1873) "Preliminary Sketch of Biquaternions"; modern
                        treatment in Yaglom (1968) "Complex Numbers in Geometry".
  Domain of validity:   All T.

METHOD
  Method declared:      Return DualNumber(exact<T>(0), exact<T>(1)).
  Method implemented:   Correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a
  Bound formula:        N/A
  Bound implemented:    N/A
  Bound verdict:        ✓ n/a

CROSS-AUDIT
  REQ-EF applies:       N/A
  AUD-EF applies:       N/A
  AUD-MC applies:       AUD-MC-3 (ε² = 0 identity)
  Verification test:    tests/test_math/ — epsilon() * epsilon() == zero()

NOTES
  - Core definition of the dual-number structure. ε² = 0 is verified in AUD-MC-3.
```

---

## Formula 7: Addition

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_add
Location:               src/math/dual_number.h:71–73
Mathematical statement: (a + εb) + (c + εd) = (a + c) + ε(b + d)

THEORY
  Underlying theorem:   Ring axiom: component-wise addition in the direct sum
                        ℝ ⊕ ε·ℝ ≅ ℝ[ε]/(ε²).
  Primary reference:    Clifford (1873); Yaglom (1968).
  Domain of validity:   All a, c ∈ TrackedValue<T>.

METHOD
  Method declared:      Component-wise addition: add reals, add duals.
  Method implemented:   `DualNumber(a.real + b.real, a.dual + b.dual)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T> rules)
  Bound formula:        Each output component (real, dual) is a TrackedValue<T>
                        addition, governed by REQ-EF-3 (triangle inequality per
                        error category). No term is introduced by the dual-number
                        algebra itself (ε² = 0 is exact).
  Bound implemented:    Delegated to TrackedValue<T>::operator+(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (error propagation via component ops)
  AUD-EF applies:       AUD-EF-1 (all ops return TrackedValue)
  AUD-MC applies:       AUD-MC-1 (addition is associative, commutative)
  Verification test:    tests/test_dual_number/ — (a+b)+c == a+(b+c),
                        a+b == b+a, a+zero == a

NOTES
  - No algebra-specific truncation; all error propagation is exact.
```

---

## Formula 8: Subtraction

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_sub
Location:               src/math/dual_number.h:76–78
Mathematical statement: (a + εb) − (c + εd) = (a − c) + ε(b − d)

THEORY
  Underlying theorem:   Ring axiom: component-wise subtraction.
  Primary reference:    Clifford (1873); Yaglom (1968).
  Domain of validity:   All a, c ∈ TrackedValue<T>.

METHOD
  Method declared:      Component-wise subtraction: subtract reals, subtract duals.
  Method implemented:   `DualNumber(a.real - b.real, a.dual - b.dual)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Each component is a TrackedValue<T> subtraction, governed
                        by REQ-EF-3.
  Bound implemented:    Delegated to TrackedValue<T>::operator-(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-1 (inverse under addition)
  Verification test:    tests/test_dual_number/ — a - a == zero, (a+b)-b == a

NOTES
  - Exact in the dual algebra.
```

---

## Formula 9: Unary negation

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_neg
Location:               src/math/dual_number.h:81–83
Mathematical statement: −(a + εb) = −a + ε(−b)

THEORY
  Underlying theorem:   Ring axiom: component-wise negation (additive inverse).
  Primary reference:    Clifford (1873); Yaglom (1968).
  Domain of validity:   All a ∈ TrackedValue<T>.

METHOD
  Method declared:      Component-wise negation.
  Method implemented:   `DualNumber(-real, -dual)` — correct.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Each component inherits the error state of its input.
                        Negation does not introduce a new error term (REQ-EF-3).
  Bound implemented:    Delegated to TrackedValue<T>::operator-(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       AUD-MC-1 (additive inverse)
  Verification test:    tests/test_dual_number/ — −(−a) == a, a + (−a) == zero

NOTES
  - Exact in the dual algebra.
```

---

## Formula 10: Multiplication (dual × dual)

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_mul_dual
Location:               src/math/dual_number.h:86–91
Mathematical statement: (a + εb)(c + εd) = ac + ε(ad + bc), since ε² = 0

THEORY
  Underlying theorem:   Ring axiom: component-wise multiplication under ε² = 0.
                        The full product would be ac + ε(ad + bc) + ε²(bd) = ac + ε(ad + bc).
  Primary reference:    Clifford (1873); Yaglom (1968). The truncation ε² = 0 is
                        the defining property of the dual-number algebra.
  Domain of validity:   All a, b, c, d ∈ TrackedValue<T>.

METHOD
  Method declared:      Dual-number algebra: compute (a + εb)(c + εd) using ε² = 0.
  Method implemented:   Lines 87–90:
                        - real part: a.real * b.real
                        - dual part: a.real * b.dual + a.dual * b.real
                        Correct. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Each output component (real, dual) is a TrackedValue<T>
                        arithmetic operation (multiplication and addition),
                        governed by REQ-EF-3. The ε² = 0 truncation is exact —
                        no term is dropped (it does not exist in the algebra).
  Bound implemented:    Delegated to TrackedValue<T> operators. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1 (all ops return TrackedValue)
  AUD-MC applies:       AUD-MC-3 (ε² = 0 is used here; multiplication is associative)
  Verification test:    tests/test_dual_number/ — (a*b)*c == a*(b*c),
                        a*identity() == a, ε*ε == zero()

NOTES
  - The ε² = 0 rule is exact (Clifford's definition), not an approximation.
  - Multiplication distributes correctly over addition (AUD-MC-1).
```

---

## Formula 11: Scalar multiplication (left)

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_mul_scalar_left
Location:               src/math/dual_number.h:94–96
Mathematical statement: s(a + εb) = sa + ε(sb)

THEORY
  Underlying theorem:   Module axiom: scalar multiplication distributes over
                        the real and dual parts.
  Primary reference:    Clifford (1873); standard module/vector-space theory.
  Domain of validity:   All s, a, b ∈ TrackedValue<T>.

METHOD
  Method declared:      Component-wise scalar multiplication.
  Method implemented:   `DualNumber(s * a.real, s * a.dual)` — correct. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Each component is a TrackedValue<T> multiplication,
                        governed by REQ-EF-3.
  Bound implemented:    Delegated to TrackedValue<T>::operator*(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Module operations (associativity of scalar multiplication)
  Verification test:    tests/test_dual_number/ — (s*t)*a == s*(t*a), s*(a+b) == s*a + s*b

NOTES
  - Exact in the dual algebra.
```

---

## Formula 12: Scalar multiplication (right)

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_mul_scalar_right
Location:               src/math/dual_number.h:99–101
Mathematical statement: (a + εb)s = as + ε(bs)

THEORY
  Underlying theorem:   Module axiom: scalar multiplication (right) distributes.
  Primary reference:    Clifford (1873); standard algebra.
  Domain of validity:   All a, b, s ∈ TrackedValue<T>.

METHOD
  Method declared:      Component-wise scalar multiplication.
  Method implemented:   `DualNumber(a.real * s, a.dual * s)` — correct. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        Each component is a TrackedValue<T> multiplication,
                        governed by REQ-EF-3.
  Bound implemented:    Delegated to TrackedValue<T>::operator*(). ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Module operations
  Verification test:    tests/test_dual_number/ — a*s == s*a, (a*b)*s == a*(b*s)

NOTES
  - Exact in the dual algebra.
```

---

## Formula 13: Division

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::operator_div
Location:               src/math/dual_number.h:113–119
Mathematical statement: (a + εb) / (c + εd) = a/c + ε(bc − ad)/c²

THEORY
  Underlying theorem:   Taylor's theorem applied to 1/x at x = c:
                        1/(c + εd) = 1/c + (−d/c²)·ε + O(ε²)
                        = 1/c − εd/c² (since ε² = 0 truncates higher orders).
                        Then multiply (a + εb) · [1/c − εd/c²]:
                        = a/c − aεd/c² + εb/c − ε²bd/c²
                        = a/c + ε(b/c − ad/c²)    (dropping ε²)
                        = a/c + ε(bc − ad)/c².
  Primary reference:    Taylor's theorem; the derivation is on lines 105–110 of
                        the source file (exact match).
                        Griewank & Walther (2008) "Evaluating Derivatives",
                        §2.2: forward-mode AD via dual numbers; the division
                        rule is standard.
  Domain of validity:   All c ≠ 0. The result is exact (no truncation error in
                        the ε algebra).

METHOD
  Method declared:      Taylor expansion (line 105: "Taylor expansion of 1/x at x = c").
  Method implemented:   Lines 114–118:
                        TrackedValue<T> c2 = b.real * b.real;
                        return DualNumber(
                            a.real / b.real,
                            (a.dual * b.real - a.real * b.dual) / c2
                        );
                        This computes:
                        - real: a/c (correct)
                        - dual: (bd − ad)/c² = (b.dual * b.real − a.real * b.dual) / c²
                               But the formula says (bc − ad)/c².
                               CRITICAL REVIEW:
                               b.real = c, b.dual = d, a.real = a, a.dual = b.
                               So the code computes (a.dual * b.real − a.real * b.dual) / c²
                               = (b·c − a·d) / c²  ✓
                        Correct. ✓
  Match verdict:        ✓ matched — the cited method is "Taylor expansion of 1/x"
                        and the implementation is exactly that: the first-order
                        Taylor of 1/(c+εd) at c, evaluated using ε² = 0.
                        No continued fraction or rational approximant is used.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        The Taylor truncation is exact in the dual algebra
                        (ε² = 0 is exact, not approximate). The real and dual
                        outputs are TrackedValue<T> arithmetic (divisions and
                        multiplications), so errors propagate per REQ-EF-3.
                        No truncation error is introduced by the dual-number
                        operation itself.
  Bound implemented:    Delegated to TrackedValue<T> division and multiplication. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (error propagation through component ops)
  AUD-EF applies:       AUD-EF-1 (returns TrackedValue)
  AUD-MC applies:       n/a (algebra identity: (a/c)·(c+εd) == (a+εb) would be tested
                            separately if needed, but division is not a primary algebra op)
  Verification test:    tests/test_dual_number/ — (a*c) / c == a, division with known
                        forward derivatives (f(x+εh) / g(x+εh) should give f/g + ε(f'g−fg')/g²)

NOTES
  - **CRITICAL FINDING**: The file comment (line 105) says "Taylor expansion" and the
    implementation is *literally* the first-order Taylor of 1/(c+εd) = 1/c − εd/c².
    This is NOT a continued fraction, NOT a Padé approximant. The comment is correct.
    The user's stated concern ("theory says Taylor, code does continued fraction") does
    not apply here. ✓
  - The exact truncation ε² = 0 eliminates any error term beyond the linear ε term.
    This is the defining feature of the dual-number algebra and is correct.
  - For wide types (cpp_bin_float_50, etc.), the method remains exact in ε and
    error-bounded through TrackedValue<T>.
```

---

## Formula 14: sqrt

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::sqrt
Location:               src/math/dual_number.h:125–128
Mathematical statement: √(a + εb) = √a + ε b/(2√a)

THEORY
  Underlying theorem:   Taylor's theorem applied to √x at x = a:
                        f(x) = √x has f'(x) = 1/(2√x).
                        f(a + εb) = f(a) + ε b f'(a) + O(ε²)
                                  = √a + ε b/(2√a)  (since ε² = 0).
  Primary reference:    Griewank & Walther (2008) "Evaluating Derivatives", §2.2
                        (forward-mode AD); also implicit in any dual-number treatment.
  Domain of validity:   a.value > 0 (sqrt is defined on positive reals).

METHOD
  Method declared:      Taylor series of √x at a, first-order in ε.
  Method implemented:   Lines 126–127:
                        TrackedValue<T> sr = sqrt(a.real);
                        return DualNumber(sr, a.dual / (exact<T>(2) * sr));
                        Computes √a and b/(2√a). Correct. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        The Taylor truncation at order ε is exact (ε² = 0).
                        The real component is TrackedValue<T>::sqrt() (governed
                        by error bounds on sqrt). The dual component is division
                        by the square root, also governed by TrackedValue<T>.
                        No algebra-specific truncation error.
  Bound implemented:    Delegated to TrackedValue<T> operations. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Algebra closure (√(a+εb) ∈ DualNumber if a+εb ∈ DualNumber
                        and a > 0)
  Verification test:    tests/test_dual_number/ — forward derivative check:
                        √(a+εh) should have dual part matching h/(2√a)

NOTES
  - Division by 2√a is safe if a > 0, which is the domain of real sqrt.
  - The derivative 1/(2√x) is correct and well-known.
```

---

## Formula 15: sin

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::sin
Location:               src/math/dual_number.h:131–133
Mathematical statement: sin(a + εb) = sin a + ε b cos a

THEORY
  Underlying theorem:   Taylor's theorem applied to sin x at x = a:
                        f(x) = sin x has f'(x) = cos x.
                        f(a + εb) = f(a) + ε b f'(a) + O(ε²)
                                  = sin a + ε b cos a  (since ε² = 0).
                        Also: sin(θ + εd) = sin θ + ε d cos θ is the
                        standard closed-form rule for dual angles in
                        screw-theory SE(3) exponentials (file header, line 17).
  Primary reference:    Griewank & Walther (2008) §2.2 (forward-mode AD);
                        Selig (2005) "Geometric Fundamentals of Robotics" §3.4
                        (screw theory);
                        Murray, Li & Sastry (1994) §3.2 (exponential coordinates).
  Domain of validity:   All a ∈ ℝ (or ℂ).

METHOD
  Method declared:      Taylor series of sin x at a, first-order in ε.
  Method implemented:   Line 132: `DualNumber(sin(a.real), a.dual * cos(a.real))`
                        Correct. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        The Taylor truncation is exact (ε² = 0).
                        Both components delegate to TrackedValue<T> operations
                        (sin, cos, multiplication), governed by REQ-EF-3.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Algebra closure; interplay with dual_quaternion::exp_screw
                        (AUD-MC-18)
  Verification test:    tests/test_dual_number/ — forward derivative check:
                        sin(a+εh) dual part should match h·cos(a)

NOTES
  - The dual-angle rule sin(θ+εd) = sin θ + ε d cos θ is foundational to
    the library's SE(3) screw-motion calculus (file header, line 17).
  - Derivative of sin is cos — standard.
```

---

## Formula 16: cos

```
=== FORMULA AUDIT CARD ===
ID:                     dual_number::cos
Location:               src/math/dual_number.h:136–138
Mathematical statement: cos(a + εb) = cos a − ε b sin a

THEORY
  Underlying theorem:   Taylor's theorem applied to cos x at x = a:
                        f(x) = cos x has f'(x) = −sin x.
                        f(a + εb) = f(a) + ε b f'(a) + O(ε²)
                                  = cos a − ε b sin a  (since ε² = 0).
                        Also: cos(θ + εd) = cos θ − ε d sin θ is the
                        standard closed-form rule for dual angles in
                        SE(3) screw-theory exponentials (file header, line 18).
  Primary reference:    Griewank & Walther (2008) §2.2 (forward-mode AD);
                        Selig (2005) "Geometric Fundamentals of Robotics" §3.4
                        (screw theory);
                        Murray, Li & Sastry (1994) §3.2 (exponential coordinates).
  Domain of validity:   All a ∈ ℝ (or ℂ).

METHOD
  Method declared:      Taylor series of cos x at a, first-order in ε.
  Method implemented:   Line 137: `DualNumber(cos(a.real), -(a.dual * sin(a.real)))`
                        Correct: computes cos a in the real part and −b sin a
                        in the dual part. ✓
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision, accuracy, measurement (per TrackedValue<T>)
  Bound formula:        The Taylor truncation is exact (ε² = 0).
                        Both components delegate to TrackedValue<T> operations
                        (cos, sin, multiplication, negation), governed by REQ-EF-3.
  Bound implemented:    Delegated. ✓
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       Algebra closure; interplay with dual_quaternion::exp_screw
                        (AUD-MC-18)
  Verification test:    tests/test_dual_number/ — forward derivative check:
                        cos(a+εh) dual part should match −h·sin(a)

NOTES
  - The dual-angle rule cos(θ+εd) = cos θ − ε d sin θ is foundational to
    SE(3) screw-motion calculus (file header, line 18).
  - Derivative of cos is −sin — standard.
  - Note the negation (line 137): `−(a.dual * sin(a.real))` correctly implements
    the minus sign from f'(x) = −sin x.
```

---

## File-level verdict

**A. Error wiring** (AUD-EF-1, AUD-EF-2, AUD-EF-7):
- All constructors and factories return either trivial values or delegates to existing constructors.
- All operators return `DualNumber` with components that are `TrackedValue<T>`. ✓
- All arithmetic and elementary functions propagate error through `TrackedValue<T>` operations (REQ-EF-3). ✓

**B. Algebra axioms** (AUD-MC-1, AUD-MC-3):
- Addition: associative, commutative, identity, inverse. ✓
- Multiplication: associative, distributive, identity (ε² = 0 exact). ✓
- Elementary functions: computed via forward-mode AD rule (f(a+εb) = f(a) + εbf'(a)) which is correct by differentiation. ✓

**C. Theoretical basis** (this document):
- **Constructors** (1–3): Trivial; no theory required. ✓
- **Factories** (4–6): Identity/zero/epsilon; ring definitions. ✓
- **Arithmetic** (7–12): Ring/module axioms. ✓
- **Division** (13): **CRITICAL** — Comment cites "Taylor expansion". Implementation is the first-order Taylor of 1/(c+εd) = 1/c − εd/c², exactly as stated. No continued fraction, no Padé. Verdict: ✓ **matched**.
- **Elementary functions** (14–16): Forward-mode AD via f(a+εb) = f(a) + εbf'(a), using ε² = 0 to truncate. All derivatives are correct (√'=1/(2√), sin'=cos, cos'=−sin). ✓

**Overall file verdict: PASS**

- All 16 formulas match their cited theory to method to bound.
- No Taylor/continued-fraction misdeclaration (the division formula is correctly labeled Taylor and correctly implemented).
- Error propagation is exact through `TrackedValue<T>` (no algebra-specific truncation).
- Algebra axioms satisfied (closure, identities, associativity, distributivity).
- Forward-mode AD is sound by construction (ε² = 0 + linear truncation).

---

## Verification checklist

- [ ] tests/test_math/test_dual_number.cc:
  - [ ] Constructor and factory tests (zero, identity, epsilon)
  - [ ] Addition, subtraction, negation, commutativity, associativity
  - [ ] Multiplication: closure, associativity, distributivity, identity
  - [ ] ε² = 0 verified explicitly: `epsilon() * epsilon() == zero()`
  - [ ] Scalar multiplication: left and right
  - [ ] Division: (a/c) * c ≈ a (within error bounds); forward-derivative check
  - [ ] sqrt: forward derivative h/(2√a)
  - [ ] sin: forward derivative h·cos(a)
  - [ ] cos: forward derivative −h·sin(a)
  - [ ] Cross-type tests (double, cpp_bin_float_50, etc.)

---

**Document**: `design/audit/theoretical_basis_audit/dual_number.md`  
**Status**: **PASS**  
**Date**: 2026-05-13  
**Analyst**: internal review
