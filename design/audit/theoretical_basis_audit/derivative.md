# Theoretical Basis Audit — `src/dynamics/derivative.h`

**File**: `src/dynamics/derivative.h` (93 lines)
**Scope**: 8 functions — 2 constructors, 2 named constructors, 4 arithmetic operators.
**Theory anchor**: REQ-EF-12 (composite types inherit the error framework by composition); Murray, Li & Sastry (1994) *A Mathematical Introduction to Robotic Manipulation* §3, which establishes that the derivative of an SE(3) trajectory at a point lies in the Lie algebra se(3) ≅ ℝ⁶ — a finite-dimensional real vector space. The arithmetic on `Derivative<T>` is therefore exactly the vector-space arithmetic on se(3) × ℝ (the extra ℝ tracks the autonomous-time component dt/dt).

**Audit posture**. `Derivative<T>` carries no novel formula. Its state is `(Twist<T> acceleration, TrackedValue<T> time_rate)`, where `Twist<T>` is itself a pair of `Vector3<T>` (which is a triple of `TrackedValue<T>`). Every numerical operation is a pre-existing scalar/vector op on `TrackedValue<T>`; the audit reduces to confirming that each operator's implementation is exactly that composition with no hidden truncation or method choice. The underlying numeric work, and its error bound, is owned by `tracked_value.h` / `vector3.h` / `twist.h`.

---

## Card 1 — `Derivative()` (default constructor)

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::Derivative()
Location:               src/dynamics/derivative.h:44
Mathematical statement: D₀ := (Twist::zero(), 1).
                        The default Derivative represents the zero
                        acceleration paired with the standard
                        time-parameterization rate dt/dt = 1.

THEORY
  Underlying theorem:   Vector-space identity. The additive identity of
                        se(3) × ℝ is (0_se(3), 0). The constructor does NOT
                        produce the additive identity — it produces
                        (0_se(3), 1), which is the natural-time-rate
                        baseline used by Munthe-Kaas / RK4 stages. (The
                        true additive identity is not constructed by any
                        function in this file; see Card 4.)
  Primary reference:    Murray, Li & Sastry (1994) §3 (se(3) as a real
                        vector space); standard ODE parameterization
                        s ↦ (state(s), t(s)) with dt/ds = 1 used by any
                        autonomous-time RK integrator.
  Domain of validity:   All T for which `math::exact<T>(1)` is exact.

METHOD
  Method declared:      Default-construct both members, then overwrite
                        `time_rate` with `exact<T>(1)`.
  Method implemented:   `acceleration()` invokes Twist's default ctor
                        (zero pair of Vector3 ≡ six exact zeros);
                        `time_rate(math::exact<T>(1))` invokes the
                        exact-1 constructor of TrackedValue<T>.
  Match verdict:        ✓ matched — pure composition, no method choice.

ERROR BOUND
  Bound category:       n/a (no truncation, no iteration, no approximation)
  Bound formula:        All six acceleration components carry the
                        zero-error TrackedValue produced by default
                        Vector3 construction. `time_rate` carries the
                        zero-error TrackedValue produced by `exact<T>(1)`
                        (per REQ-EF-3, `exact<T>` constructs a value with
                        precision = accuracy = measurement = 0).
  Bound implemented:    Inherited from `Twist<T>()` and
                        `math::exact<T>(1)`; nothing added here.
  Bound verdict:        ✓ matched — REQ-EF-12 composition only.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12 (composite types inherit by composition);
                        REQ-EF-3 (exact constants carry zero error).
  AUD-EF applies:       AUD-EF-1 (every public op returns/produces
                        TrackedValue-bearing types — here the members
                        are TrackedValue-bearing aggregates).
  AUD-MC applies:       n/a (no algebra operation).
  Verification test:    tests covering `Derivative()` default state should
                        assert `acceleration == Twist::zero()` and
                        `time_rate.value == 1` with `total_error() == 0`.

NOTES
  - This is the most common starting point for force-model code that
    fills in acceleration and leaves time_rate untouched. It is NOT the
    additive identity (which is `zero()` — see Card 4 — but note that
    even `zero()` here has time_rate = 1, not 0; see Card 4's notes).
```

---

## Card 2 — `Derivative(const Twist<T>&, const TrackedValue<T>&)` (parameterized constructor)

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::Derivative(Twist, TrackedValue)
Location:               src/dynamics/derivative.h:46-47
Mathematical statement: D := (a, dt/dt). Copy-construct both fields.

THEORY
  Underlying theorem:   Direct-sum construction. The product space
                        V := se(3) × ℝ is a vector space; its elements
                        are uniquely identified by their two components.
                        This constructor packages two given components
                        into a single element.
  Primary reference:    Murray, Li & Sastry (1994) §3; any linear algebra
                        text on direct-sum spaces.
  Domain of validity:   Any T, any well-formed Twist<T>, any well-formed
                        TrackedValue<T>.

METHOD
  Method declared:      Member-wise copy.
  Method implemented:   `acceleration(a)` and `time_rate(dt_dt)`, i.e.
                        Twist's copy ctor and TrackedValue's copy ctor.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (copy only, no new arithmetic).
  Bound formula:        TrackedValue copy preserves errors per REQ-EF-3
                        (copy ctor is exact: value, precision, accuracy,
                        measurement all preserved bit-for-bit).
  Bound implemented:    Inherited from member copy constructors.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12.
  AUD-EF applies:       AUD-EF-1.
  AUD-MC applies:       n/a.
  Verification test:    Construct with non-default arguments and assert
                        component-wise equality including errors.

NOTES
  - No new arithmetic. The audit is closed by REQ-EF-12 composition.
```

---

## Card 3 — `from_acceleration(const Twist<T>&)`

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::from_acceleration
Location:               src/dynamics/derivative.h:53-55
Mathematical statement: from_acceleration(a) := (a, 1).
                        Lift a body acceleration into a Derivative under
                        the standard parameterization dt/dt = 1.

THEORY
  Underlying theorem:   Standard autonomous-time ODE parameterization:
                        for an autonomous system ẋ = f(x), augmenting
                        the state with t and treating it as another
                        variable gives the autonomous form
                        (ẋ, ṫ) = (f(x), 1).
  Primary reference:    Hairer, Nørsett & Wanner (1993) *Solving Ordinary
                        Differential Equations I* §I.2 (autonomous form);
                        Murray, Li & Sastry (1994) §3 (acceleration as an
                        element of se(3)).
  Domain of validity:   Any T for which `math::exact<T>(1)` is exact.

METHOD
  Method declared:      Pair the given acceleration with the exact
                        constant 1 in the time-rate slot.
  Method implemented:   `return Derivative(a, math::exact<T>(1));`
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (no arithmetic on `a`; `exact<T>(1)` carries
                        zero error per REQ-EF-3).
  Bound formula:        Errors of `a` are preserved bit-for-bit; the
                        time_rate carries exact-zero error.
  Bound implemented:    Inherited from `exact<T>(1)` and the copy of `a`.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12; REQ-EF-3 (exact constants).
  AUD-EF applies:       AUD-EF-1.
  AUD-MC applies:       n/a.
  Verification test:    Build with a non-zero Twist, assert acceleration
                        round-trips and time_rate.value == 1 with zero
                        error.

NOTES
  - This is the canonical force-model output: the force model computes
    body acceleration; this function packages it for the integrator.
```

---

## Card 4 — `zero()`

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::zero
Location:               src/dynamics/derivative.h:58-60
Mathematical statement: zero() := (0_Twist, 1).
                        "Zero acceleration, time still advances at rate 1."

THEORY
  Underlying theorem:   Convention (NOT a Lie-algebra zero). The
                        Lie-algebra additive identity of se(3) × ℝ would
                        be (0_se(3), 0). This function deliberately
                        returns (0_se(3), 1) to keep time advancing at
                        the standard parameterization rate when no
                        acceleration is applied.
  Primary reference:    same parameterization choice as Card 3.
  Domain of validity:   Any T for which `math::exact<T>(1)` is exact.

METHOD
  Method declared:      Pair Twist::zero() with exact<T>(1).
  Method implemented:   `return Derivative(Twist<T>::zero(), math::exact<T>(1));`
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a.
  Bound formula:        Both members carry zero error (Twist::zero() by
                        construction; exact<T>(1) per REQ-EF-3).
  Bound implemented:    As declared.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12; REQ-EF-3.
  AUD-EF applies:       AUD-EF-1.
  AUD-MC applies:       n/a — note that the four arithmetic operators
                        below do NOT form an abelian group with this
                        `zero()`: `D - D` produces (0, 0), not zero() (
                        which has time_rate = 1). See NOTES below.
  Verification test:    Assert `zero().acceleration == Twist::zero()` and
                        `zero().time_rate.value == 1` with zero error.

NOTES
  - ⚠ Naming caveat. `Derivative::zero()` does NOT satisfy
    `D + Derivative::zero() == D` for arbitrary D: it adds 1 to D's
    time_rate. Strictly, the additive identity is `Derivative()` with
    `time_rate` overwritten to 0, which no named constructor produces.
    This is acceptable for RK4-style integrators that combine
    derivatives but never substitute `zero()` for an additive identity
    — the integrator builds new Derivatives via `from_acceleration` and
    weights them with TrackedValue scalars (Cards 7-8). The audit flag
    `⚠` is recorded here to surface a name-vs-mathematical-property
    mismatch that could mislead future callers; it is not a C-fail in
    the current usage pattern.
  - For the file-level verdict, this is a documentation / API issue,
    not a method-vs-theory mismatch. The implementation matches what
    the comment "(no acceleration, time still advances at rate 1)"
    declares.
```

---

## Card 5 — `operator+(Derivative, Derivative)`

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::operator+
Location:               src/dynamics/derivative.h:65-70
Mathematical statement: a + b := (a.acceleration + b.acceleration,
                                  a.time_rate + b.time_rate).
                        Component-wise vector-space addition on
                        se(3) × ℝ.

THEORY
  Underlying theorem:   Vector-space addition on the direct sum
                        se(3) × ℝ. Component-wise definition is the
                        unique extension that makes both projections
                        linear.
  Primary reference:    Murray, Li & Sastry (1994) §3 (se(3) is a
                        real vector space); standard linear algebra
                        on direct-sum spaces.
  Domain of validity:   All T for which `TrackedValue<T>` addition is
                        defined.

METHOD
  Method declared:      Member-wise call to Twist's `operator+` and
                        TrackedValue's `operator+`.
  Method implemented:   `Derivative(a.acceleration + b.acceleration,
                                    a.time_rate + b.time_rate)`.
  Match verdict:        ✓ matched — pure composition.

ERROR BOUND
  Bound category:       Inherited (precision / accuracy / measurement
                        propagated by the inner adds).
  Bound formula:        TrackedValue addition adds the absolute errors
                        per category (REQ-EF-3): for x + y, the
                        precision/accuracy/measurement of the result are
                        the sums of those of x and y. Vector3 + and
                        Twist + propagate this component-wise. This file
                        adds no new error.
  Bound implemented:    Delegated to `Twist<T>::operator+` (which
                        delegates to `Vector3<T>::operator+`, which
                        delegates to `TrackedValue<T>::operator+`) and
                        `TrackedValue<T>::operator+`.
  Bound verdict:        ✓ matched — REQ-EF-12 composition only.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12; REQ-EF-3 (closed-form propagation).
  AUD-EF applies:       AUD-EF-1; AUD-EF-7 (composite arithmetic
                        propagates errors).
  AUD-MC applies:       Vector-space axioms (commutativity / associativity
                        of `+`) inherited from Twist's and TrackedValue's
                        algebra audits; no new axiom proof needed here.
  Verification test:    Symmetric/associativity sanity test on
                        Derivative: `(a + b) + c == a + (b + c)` and
                        `a + b == b + a` to within composed error.

NOTES
  - No method choice (no Taylor, Padé, etc.). REQ-EF-12 carries the
    audit.
```

---

## Card 6 — `operator-(Derivative, Derivative)`

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::operator-
Location:               src/dynamics/derivative.h:73-78
Mathematical statement: a - b := (a.acceleration - b.acceleration,
                                  a.time_rate - b.time_rate).
                        Component-wise vector-space subtraction on
                        se(3) × ℝ.

THEORY
  Underlying theorem:   Vector-space subtraction = addition with the
                        additive inverse. Same direct-sum reasoning as
                        Card 5.
  Primary reference:    Murray, Li & Sastry (1994) §3.
  Domain of validity:   All T for which `TrackedValue<T>` subtraction
                        is defined.

METHOD
  Method declared:      Member-wise call to Twist's `operator-` and
                        TrackedValue's `operator-`.
  Method implemented:   `Derivative(a.acceleration - b.acceleration,
                                    a.time_rate - b.time_rate)`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       Inherited (precision / accuracy / measurement).
  Bound formula:        TrackedValue subtraction adds absolute errors
                        per category (errors do not cancel — REQ-EF-3
                        treats subtraction as |x| + |y| in the error
                        components, with potential catastrophic-
                        cancellation amplification when value(a-b) is
                        small relative to |a|+|b|; the bound itself is
                        still rigorous, the relative degradation is a
                        property of the inputs, not the method).
  Bound implemented:    Delegated to `Twist<T>::operator-`, ultimately
                        to `TrackedValue<T>::operator-`.
  Bound verdict:        ✓ matched — REQ-EF-12 composition only.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12; REQ-EF-3.
  AUD-EF applies:       AUD-EF-1; AUD-EF-7.
  AUD-MC applies:       Inherited; `a - a == zero-pair` (NOT
                        `Derivative::zero()` — see Card 4 NOTES).
  Verification test:    `a - a` should produce zero acceleration AND
                        zero time_rate (the additive identity of the
                        product vector space, not the named `zero()`).

NOTES
  - Subtraction is symmetric in structure with addition; the audit is
    the same modulo sign.
  - The "additive identity is not equal to `Derivative::zero()`" point
    flagged in Card 4 surfaces here: `a - a` yields (0, 0), whereas
    `Derivative::zero()` yields (0, 1). Currently no integrator code
    in the repo relies on this distinction, but a unit test that
    asserted `a - a == Derivative::zero()` would (correctly) fail.
```

---

## Card 7 — `operator*(TrackedValue<T> scalar, Derivative)` (scalar on left)

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::operator*(scalar, derivative)
Location:               src/dynamics/derivative.h:81-84
Mathematical statement: s · d := (s · d.acceleration, s · d.time_rate).
                        Scalar-vector multiplication on se(3) × ℝ.

THEORY
  Underlying theorem:   Scalar multiplication is the bilinear map
                        ℝ × V → V (here V = se(3) × ℝ), distributed
                        component-wise on the direct sum.
  Primary reference:    Murray, Li & Sastry (1994) §3 (se(3) as a real
                        vector space — scalar multiplication is defined
                        component-wise on ω, v).
  Domain of validity:   All T; any TrackedValue scalar `s`.

METHOD
  Method declared:      Member-wise scalar multiplication via Twist's
                        `operator*(scalar, twist)` and TrackedValue's
                        `operator*`.
  Method implemented:   `Derivative(s * d.acceleration, s * d.time_rate)`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       Inherited (precision / accuracy / measurement)
                        via the multiplicative error-propagation rule.
  Bound formula:        TrackedValue scalar product per REQ-EF-3 uses
                        |s|·err(d) + |d|·err(s) per category (the
                        Leibniz-style first-order propagation). Vector3
                        and Twist multiplications apply this component-
                        wise. This file adds no new error.
  Bound implemented:    Delegated to `Twist<T>::operator*(scalar, twist)`
                        and `TrackedValue<T>::operator*`.
  Bound verdict:        ✓ matched — REQ-EF-12 composition only.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12; REQ-EF-3.
  AUD-EF applies:       AUD-EF-1; AUD-EF-7.
  AUD-MC applies:       Vector-space axiom: scalar associativity
                        `(s·t)·d == s·(t·d)` and distributivity
                        `s·(a + b) == s·a + s·b` — inherited from
                        Twist's and TrackedValue's algebra audits.
  Verification test:    RK4-stage weighting uses this operator;
                        coverage via tests/test_dynamics integrators.

NOTES
  - Scalar on left is the convention used by RK4 / Munthe-Kaas stage
    combiners (`(h/6) * (k1 + 2*k2 + 2*k3 + k4)` form).
```

---

## Card 8 — `operator*(Derivative, TrackedValue<T> scalar)` (scalar on right)

```
=== FORMULA AUDIT CARD ===
ID:                     dynamics::Derivative::operator*(derivative, scalar)
Location:               src/dynamics/derivative.h:87-90
Mathematical statement: d · s := (d.acceleration · s, d.time_rate · s).
                        Scalar-vector multiplication, scalar on the right.

THEORY
  Underlying theorem:   Same as Card 7. Real scalar multiplication is
                        commutative on the ℝ-vector space se(3) × ℝ,
                        so left- and right-multiplication agree.
                        (Caveat: if T is non-commutative — which it is
                        not for any T used in this codebase: T ∈ {double,
                        cpp_bin_float_50, …} — left and right forms
                        would differ. The implementation respects the
                        left/right order in the inner calls, so a future
                        non-commutative T would be handled correctly at
                        the composition layer; the algebra-axiom audit
                        would then be the determinative check.)
  Primary reference:    Murray, Li & Sastry (1994) §3.
  Domain of validity:   All T; any TrackedValue scalar `s`.

METHOD
  Method declared:      Member-wise scalar multiplication with scalar
                        on the right.
  Method implemented:   `Derivative(d.acceleration * s, d.time_rate * s)`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       Inherited (precision / accuracy / measurement).
  Bound formula:        Same as Card 7 — Leibniz-style first-order
                        propagation, delegated to inner ops.
  Bound implemented:    Delegated to `Twist<T>::operator*(twist, scalar)`
                        and `TrackedValue<T>::operator*`.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12; REQ-EF-3.
  AUD-EF applies:       AUD-EF-1; AUD-EF-7.
  AUD-MC applies:       Commutativity of scalar multiplication
                        `s·d == d·s` for commutative T — verified by the
                        algebra audit on TrackedValue, inherited here.
  Verification test:    `s * d == d * s` for all-commutative T should be
                        a unit assertion at the Derivative level.

NOTES
  - Provided for symmetry / call-site ergonomics; mathematically
    identical to Card 7 under any commutative T.
```

---

## File-level verdict — `derivative.h`

- **A. Error wiring**: ✓ All eight functions are pure compositions of `Twist<T>` and `TrackedValue<T>` operations; the error budget is propagated entirely by the inner ops (REQ-EF-12). No bare `T` slips through; no error category is mishandled at this layer.
- **B. Algebra axioms**: vector-space axioms on `Derivative<T>` (commutativity, associativity of `+`; distributivity and associativity of scalar `·`) are inherited from `Twist<T>` and `TrackedValue<T>`. No new axiom obligation introduced at this layer. ✓
- **C. Theoretical basis**:
  - Cards 1, 2, 3: ✓ matched (constructors / copy / lift).
  - Card 4 `zero()`: ⚠ name-vs-mathematical-property caveat — `zero()` is not the additive identity of the underlying direct-sum vector space (it sets `time_rate = 1`, not `0`). Implementation matches what the comment declares; the warning is a documentation / API smell, not a method-theory mismatch. No C-fail. Flag for future API tightening: either rename to `natural_time_rate()` or introduce an explicit `additive_identity()` that returns `(0_Twist, 0)`.
  - Cards 5, 6, 7, 8: ✓ matched (vector-space addition, subtraction, scalar multiplication).

**File verdict: PASS** — all formulas are vector-space-composition per cited theory; bounds are inherited from `tracked_value.h` / `vector3.h` / `twist.h`. One ⚠ note on `zero()` API naming; no C-fails; no method-theory mismatches; no novel arithmetic at this layer.

**Cross-references**:
- Composition chain audited end-to-end requires the corresponding PASS in `theoretical_basis_audit/tracked_value.md`, `theoretical_basis_audit/vector3.md`, and `theoretical_basis_audit/twist.md`. Until those land, this file's PASS is conditional on REQ-EF-12 propagation through those (currently `NEEDED`) cards.
