# Theoretical Basis Audit — `src/dynamics/inertia.h`

**Source file**: `src/dynamics/inertia.h` (150 lines)
**Functions audited**: 8
**Companion framework**: `design/audit/theoretical_basis_audit.md` §1 (per-formula card schema), §5 (worked example).

## Scope

`Inertia<T>` is the rigid-body **dual inertia operator** **H** in the dual-quaternion formulation. Per Selig (2005) §5 and Murray-Li-Sastry (1994) Ch.4, the dual inertia maps a body twist Ω̂ = ω + ε v to body generalized momentum H Ω̂ = (I_body ω) + ε (m v). The file implements **H** (`momentum_of`), **H⁻¹** (`acceleration_from_wrench`), the point-mass restriction of **H⁻¹** (`linear_acceleration_from_force`), plus three constructors that build elements of the **diagonal subclass** of inertia tensors.

Theory anchors used throughout this document:
- **Murray, Li & Sastry (1994)** *A Mathematical Introduction to Robotic Manipulation*, Ch.4 "Rigid Body Velocity and Acceleration", §4.2 generalized inertia tensor for SE(3).
- **Featherstone (2008)** *Rigid Body Dynamics Algorithms*, Ch.2 "Spatial Vector Algebra", §2.13 spatial inertia matrix.
- **Selig (2005)** *Geometric Fundamentals of Robotics*, §5 "Robot Dynamics" — the dual quaternion / Clifford algebra inertia.
- **REQ-EF-9 / AUD-EF-9**: safe-div semantics for the documented point-mass zero/zero degenerate case.

---

## 1 `Inertia()` — default constructor

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::default_ctor
Location:               src/dynamics/inertia.h:47
Mathematical statement: H_default = (m=0, I_body = diag(0,0,0)) ∈ S²(se(3)*)
                        The zero element of the cone of dual inertias.

THEORY
  Underlying theorem:   Identity / default initialization of a value type.
                        The diagonal subclass of SE(3) generalized inertia
                        tensors S²_diag ⊂ S²(se(3)*) is a closed convex cone
                        containing 0 as its origin (Murray-Li-Sastry §4.2).
                        Zero is a legal (degenerate) point in the cone.
  Primary reference:    Murray-Li-Sastry (1994) Ch.4 eq.(4.27) — generalized
                        inertia is a symmetric positive semidefinite operator;
                        Selig (2005) §5.1 — dual inertia is the SE(3) analog.
  Domain of validity:   Always valid as a value-initialized object; the
                        resulting Inertia is degenerate (both `momentum_of`
                        and `acceleration_from_wrench` then have all-zero or
                        all-zero/zero responses governed by REQ-EF-9).

METHOD
  Method declared:      Closed-form: m ← 0, I_body ← (0, 0, 0).
  Method implemented:   `mass()` and `principal_moments()` value-initialize
                        their `TrackedValue<T>` / `Vector3<T>` members to
                        zero with `errors = (0,0,0,0)`.
  Match verdict:        ✓ matched — value-initialization yields the zero
                        element with zero error.

ERROR BOUND
  Bound category:       n/a (no operation; pure initialization)
  Bound formula:        N/A — a defaulted value introduces no error.
                        Per REQ-EF-3, `exact` constants carry zero error.
  Bound implemented:    All four error categories are 0.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form construction; zero error for
                        exact constants).
  AUD-EF applies:       AUD-EF-1 (TrackedValue<T> default carries the
                        empty error budget).
  AUD-MC applies:       n/a (no algebra operation invoked).
  Verification test:    Implicit in any test that uses `Inertia<T>{}` and
                        checks `total_error() == 0`.

NOTES
  - The default-constructed `Inertia` is **degenerate** (zero moment +
    zero mass). Using it as input to `acceleration_from_wrench` triggers
    REQ-EF-9 safe-div pathology on both linear and angular sides.
    This is intentional; documented in the file header.
```

---

## 2 `Inertia(const TrackedValue<T>& m, const Vector3<T>& I)` — parameterized constructor

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::param_ctor
Location:               src/dynamics/inertia.h:49-50
Mathematical statement: H = diag(I_xx, I_yy, I_zz, m, m, m) — the body-frame
                        dual inertia operator with explicit principal moments
                        and mass. Element of the diagonal subclass
                        S²_diag(se(3)*).

THEORY
  Underlying theorem:   Definition. The dual quaternion inertia operator is
                        the block-diagonal map (ω, v) ↦ (I_body ω, m v); for
                        a diagonal I_body this collapses to six independent
                        scalar multiplications (Featherstone (2008) §2.13;
                        Selig (2005) §5.2 for the dual quaternion form).
  Primary reference:    Featherstone (2008) eq.(2.62) spatial inertia matrix;
                        Murray-Li-Sastry (1994) eq.(4.30) (block-diagonal in
                        the principal axes of I_body).
  Domain of validity:   m ≥ 0, I_ii ≥ 0 for physical realizability (positive
                        semidefiniteness). The constructor does NOT enforce
                        this; callers / named constructors do.

METHOD
  Method declared:      Closed-form copy: store (m, I) by value.
  Method implemented:   Member-init list `mass(m), principal_moments(I)`.
  Match verdict:        ✓ matched — exact copy, no arithmetic.

ERROR BOUND
  Bound category:       n/a (no new error is introduced; inputs are stored
                        verbatim).
  Bound formula:        Per REQ-EF-3, copy/assignment is closed-form: output
                        carries input's error vector unchanged.
  Bound implemented:    Errors flow through the `TrackedValue` and `Vector3`
                        copy constructors unmodified.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form copy).
  AUD-EF applies:       AUD-EF-1 (no error category is wrongly cleared).
  AUD-MC applies:       n/a.
  Verification test:    Implicit via any test that constructs an Inertia
                        from a TrackedValue with known error budget and
                        verifies the budget is preserved.

NOTES
  - The struct is the carrier of "diagonal subclass" inertia. Positivity
    of m, I_ii is **caller's responsibility**; the audit does not flag this
    as a C-fail because it is not a method/bound mismatch — it is a
    domain enforcement decision documented in the file header.
```

---

## 3 `point_mass(m)` — named constructor

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::point_mass
Location:               src/dynamics/inertia.h:56-58
Mathematical statement: H_pt(m) = diag(0, 0, 0, m, m, m).
                        A degenerate inertia with zero rotational moment;
                        applicable to point-mass orbital propagation where
                        attitude dynamics are not integrated.

THEORY
  Underlying theorem:   Definition — the point-mass model is the m → m,
                        I_body → 0 limit of the diagonal inertia. In the
                        Newtonian formulation, a particle of mass m has
                        zero second-moment-of-mass and so zero I_body
                        (Murray-Li-Sastry §4.2.1 distinguishes particle
                        from rigid body precisely on the I_body ≠ 0 property).
  Primary reference:    Murray-Li-Sastry (1994) §4.2.1; Featherstone (2008)
                        §2.13 (the m=0 → invalid case is the symmetric
                        opposite; point-mass with m>0 and I=0 is the valid
                        partner).
  Domain of validity:   m > 0 (caller's responsibility). The propagator
                        must be configured to skip attitude integration
                        (documented file-header precondition).

METHOD
  Method declared:      Closed-form: I_body ← (0, 0, 0), mass ← m.
  Method implemented:   `Inertia(m, math::Vector3<T>())` — default Vector3<T>
                        is zero per Vector3<T>'s value-initialization.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (closed-form copy).
  Bound formula:        Per REQ-EF-3, the constructed Inertia carries m's
                        full error vector and zero error on I_body.
  Bound implemented:    Identical to declared.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form), REQ-EF-9 (safe-div for the
                        degenerate angular slot when this Inertia is later
                        passed to `acceleration_from_wrench`).
  AUD-EF applies:       AUD-EF-1, AUD-EF-9.
  AUD-MC applies:       n/a (this constructor doesn't do algebra; the
                        consumer `acceleration_from_wrench` does).
  Verification test:    tests/test_dynamics — point-mass orbital propagation
                        should round-trip an `Inertia::point_mass(m)` through
                        `linear_acceleration_from_force` with zero attitude
                        contribution; the test should also verify that
                        `acceleration_from_wrench` with a torque-free Wrench
                        returns ω_dot = 0 (the documented zero/zero branch
                        of REQ-EF-9).

NOTES
  - This is the **intended use-mode** for the orbital propagator (single
    point mass in central + zonal gravity + drag; no attitude). The named
    constructor exists to make this intent explicit at call sites and to
    pair semantically with the REQ-EF-9 safe-div in
    `acceleration_from_wrench`.
```

---

## 4 `uniform_sphere(m, r)` — named constructor

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::uniform_sphere
Location:               src/dynamics/inertia.h:62-69
Mathematical statement: I_sphere(m, r) = (2/5) m r²  (per principal axis).
                        H_sphere = diag(2/5 mr², 2/5 mr², 2/5 mr², m, m, m).

THEORY
  Underlying theorem:   The principal moment of inertia of a uniform solid
                        sphere of mass m and radius r about any axis through
                        its center is I = (2/5) m r². This is the integral
                        ∫ ρ(x² + y²) dV over the ball B(r), evaluated in
                        spherical coordinates, with the ball's symmetry
                        making every axis principal.
  Primary reference:    Goldstein, Poole & Safko (2002) *Classical Mechanics*,
                        Ch.5 "The Rigid Body Equations of Motion", §5.3
                        Table 5.1 (or any introductory mechanics text);
                        ProofWiki — *Moment of Inertia of Uniform Solid Ball*.
  Domain of validity:   m > 0, r > 0; density assumed uniform.

METHOD
  Method declared:      Closed-form scalar multiplication: I = (2/5) · m · r²
                        replicated on all three principal axes.
  Method implemented:   `math::ratio<T>(2, 5) * m * r * r` (one
                        TrackedValue<T> value), then broadcast into a
                        Vector3<T> with all three components equal.
  Match verdict:        ✓ matched — explicit (2/5) · m · r · r closed-form;
                        no series, no iteration.

ERROR BOUND
  Bound category:       precision (and accuracy on m, r inputs).
  Bound formula:        Per REQ-EF-3, the resulting TrackedValue error is
                        the closed-form propagated bound for the product
                        (2/5) · m · r · r:
                          δI ≈ |I| · (δm/|m| + 2·δr/|r|),
                        with the ratio (2/5) introducing zero error
                        (it is an `exact` rational).
                        Each principal-moment component carries this same
                        scalar error.
  Bound implemented:    TrackedValue<T>'s overloaded `*` propagates the
                        product-rule error per REQ-EF-3 / AUD-EF-1. The
                        Vector3<T> constructor copies the scalar
                        TrackedValue into all three components.
  Bound verdict:        ✓ matched — closed-form product-rule propagation
                        as specified by REQ-EF-3.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form product error propagation).
  AUD-EF applies:       AUD-EF-1 (TrackedValue<T> arithmetic preserves
                        error categories).
  AUD-MC applies:       n/a (this is a parameter-builder, not an algebra op).
  Verification test:    tests/test_dynamics — verify uniform_sphere(m, r)
                        agrees with diagonal(m, 2mr²/5, 2mr²/5, 2mr²/5)
                        and that the error budget matches a manual
                        computation of δI = |I|·(δm/|m| + 2δr/|r|).

NOTES
  - The factor (2/5) is `math::ratio<T>(2,5)`, which per the constants
    convention carries zero error (it is a closed-form exact rational).
  - All three principal-moment components share the **same**
    TrackedValue<T> instance (broadcast assignment). This is correct —
    the components are deterministically equal by spherical symmetry;
    no statistical independence is implied that would require uncorrelated
    error budgets.
```

---

## 5 `diagonal(m, I_xx, I_yy, I_zz)` — named constructor

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::diagonal
Location:               src/dynamics/inertia.h:72-77
Mathematical statement: H_diag = diag(I_xx, I_yy, I_zz, m, m, m). The
                        general element of the diagonal subclass
                        S²_diag(se(3)*) with explicit principal moments.

THEORY
  Underlying theorem:   Definition / packaging. In the principal-axes frame
                        of a rigid body, I_body is diagonal by the spectral
                        theorem applied to the symmetric inertia tensor
                        (Murray-Li-Sastry §4.2 prop.4.5; Goldstein §5.4).
                        Any I_body in body-frame coordinates aligned with
                        principal axes has only three independent entries.
  Primary reference:    Murray-Li-Sastry (1994) §4.2; Featherstone (2008)
                        §2.13; the diagonal form is the canonical
                        parameterization in body coordinates.
  Domain of validity:   Triangle inequalities I_xx + I_yy ≥ I_zz (and
                        cyclic permutations) for physical realizability;
                        I_ii > 0 for non-degenerate rotational dynamics.
                        The constructor does NOT enforce these; callers do.

METHOD
  Method declared:      Closed-form copy: package three scalars into a
                        Vector3<T> and forward to the parameterized
                        constructor.
  Method implemented:   `Inertia(m, math::Vector3<T>(I_xx, I_yy, I_zz))`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       n/a (closed-form copy; no arithmetic).
  Bound formula:        Per REQ-EF-3, each input TrackedValue's error
                        vector flows verbatim into the corresponding
                        component of `principal_moments`.
  Bound implemented:    Identical to declared via Vector3<T>'s constructor.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3.
  AUD-EF applies:       AUD-EF-1.
  AUD-MC applies:       n/a.
  Verification test:    tests/test_dynamics — verify diagonal(m, I_xx, I_yy,
                        I_zz) yields the expected diagonal H by inspecting
                        the resulting `principal_moments` and `mass`.

NOTES
  - Triangle-inequality enforcement is the caller's responsibility (same
    convention as for the parameterized constructor in §2). This is a
    domain-of-input precondition, not a method/bound mismatch — not a
    C-fail.
```

---

## 6 `momentum_of(twist)` — body generalized momentum

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::momentum_of
Location:               src/dynamics/inertia.h:86-98
Mathematical statement: Given a body twist Ω̂ = (ω, v), compute the body
                        generalized momentum (L, p) = H Ω̂ where:
                          L_i = I_ii · ω_i    (angular momentum, i ∈ {x,y,z})
                          p_i = m · v_i       (linear momentum,  i ∈ {x,y,z})
                        Packed in `Wrench<T>` shape: (torque=L, force=p).

THEORY
  Underlying theorem:   The dual inertia operator H : se(3) → se(3)* acting
                        on a twist gives the body generalized momentum.
                        For a diagonal I_body, the angular leg is the
                        principal-moments diagonal product L_i = I_ii ω_i;
                        the linear leg is the scalar mass product p = m v.
                        This is the diagonal restriction of Murray-Li-Sastry
                        eq.(4.32) and Featherstone eq.(2.63), and the
                        Clifford-algebra dual-quaternion form in
                        Selig (2005) §5.2.
  Primary reference:    Murray, Li & Sastry (1994) §4.2 eq.(4.32);
                        Featherstone (2008) §2.13 eq.(2.63);
                        Selig (2005) §5.2.
  Domain of validity:   All (ω, v) ∈ se(3). No singularities in this
                        direction (H is bounded above by max(m, max I_ii));
                        only the inverse H⁻¹ has the m=0 / I_ii=0 issues.

METHOD
  Method declared:      Closed-form: six independent scalar products,
                        component-by-component. No iteration, no series.
  Method implemented:   Six explicit scalar multiplies via TrackedValue<T>
                        `*` operator, packed into two Vector3<T> and then
                        a Wrench<T>:
                          L = (I_xx*ω_x, I_yy*ω_y, I_zz*ω_z)
                          p = (m*v_x,    m*v_y,    m*v_z)
                        Match verdict: ✓ matched.

ERROR BOUND
  Bound category:       precision and accuracy propagated from inputs;
                        no new bound introduced (closed-form arithmetic).
  Bound formula:        Per REQ-EF-3, scalar product `a * b` propagates
                        error per the product rule:
                          δ(ab) ≈ |a|·δb + |b|·δa,
                        across all four error categories independently
                        (precision, accuracy, measurement, model).
                        Per AUD-EF-1, the bound is what TrackedValue<T>'s
                        operator* implements.
  Bound implemented:    TrackedValue<T>::operator* propagates per REQ-EF-3.
                        The Vector3<T> and Wrench<T> packaging carry the
                        per-component errors through verbatim.
  Bound verdict:        ✓ matched (delegates to AUD-EF-1's product-rule
                        verification for TrackedValue<T>).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form product).
  AUD-EF applies:       AUD-EF-1 (TrackedValue<T> product-rule wiring).
  AUD-MC applies:       AUD-MC-X (SE(3) inertia algebra — the symmetry
                        property H = H^T is checked via the pairing
                        ⟨H Ω̂_1, Ω̂_2⟩ = ⟨Ω̂_1, H Ω̂_2⟩; for a diagonal
                        H this is trivially true). Use Wrench-Twist
                        pairing tests if available.
  Verification test:    tests/test_dynamics — exercise `momentum_of`
                        with a known (ω, v) and a diagonal I_body, verify
                        L_i = I_ii ω_i numerically and that the error
                        budget matches the product-rule propagation.

NOTES
  - The packing of (L, p) into `Wrench<T>` is a **type reuse** for shape,
    not a semantic claim that body momentum **is** a wrench. The class-
    docstring (lines 84-85) explicitly notes this. From the TBA viewpoint,
    this does not introduce any method/bound mismatch — the math is
    identical regardless of which carrier struct holds the six numbers.
  - The off-diagonal I_body case (full Matrix3<T>) is deferred; the
    file-header docstring notes this. The audit therefore covers only
    the diagonal subclass.
```

---

## 7 `acceleration_from_wrench(wrench)` — inverse dual inertia H⁻¹ with safe-div

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::acceleration_from_wrench
Location:               src/dynamics/inertia.h:114-140
Mathematical statement: Given a body wrench (τ, F), compute the body
                        acceleration (ω_dot, v_dot) = H⁻¹ (τ, F) where:
                          ω_dot_i = τ_i / I_ii          (i ∈ {x,y,z})
                          v_dot_i = F_i / m
                        With the degenerate-axis handling per REQ-EF-9:
                          - I_ii = 0 ∧ τ_i = 0   →   ω_dot_i := 0
                          - I_ii = 0 ∧ τ_i ≠ 0   →   ω_dot_i := MAX (catastrophe)
                          - I_ii ≠ 0             →   ω_dot_i := τ_i / I_ii

THEORY
  Underlying theorem:   H is positive-definite on its rank-6 image (Murray-
                        Li-Sastry §4.2 prop.4.5) when m > 0 ∧ I_body ≻ 0;
                        H⁻¹ is then the component-wise reciprocal product
                        on the diagonal subclass: (τ_i, F_i) ↦ (τ_i/I_ii,
                        F_i/m). The zero-moment limit is a rank-defective
                        boundary point of the cone S²_diag(se(3)*); it is
                        outside the domain of H⁻¹ as an operator, but is
                        handled here as a documented degenerate case per
                        the file-header point-mass convention.
  Primary reference:    Murray-Li-Sastry (1994) §4.2; Selig (2005) §5.2;
                        for the REQ-EF-9 safe-div semantics:
                        design/specifications/error_framework.md REQ-EF-9
                        and design/audit/error_framework.md AUD-EF-9.
  Domain of validity:   m > 0 (precondition; documented at line 113); I_ii
                        either > 0 or exactly 0 with τ_i = 0 (point-mass
                        propagator path).

METHOD
  Method declared:      Closed-form: six independent scalar quotients with
                        an inline `safe_div` lambda guarding the angular
                        leg. The linear leg unconditionally divides by m
                        (per precondition mass > 0).
  Method implemented:   Lines 115-128: lambda `safe_div(num, den)` with
                        branching on `den.value == T(0)`:
                          - num.value == 0  → return `exact<T>(0)`
                          - num.value ≠ 0  → return TrackedValue(MAX, 0,
                                              MAX, 0)
                          - else            → return `num / den`
                        Lines 129-133: `alpha = (safe_div(τ_i, I_ii))_i`.
                        Lines 134-138: `a = (F_i / m)_i`, no guard
                        (m > 0 precondition).
  Match verdict:        ✓ matched — declared and implemented match
                        line-for-line. The zero/zero return value
                        `exact<T>(0)` is per REQ-EF-9. The nonzero/zero
                        return `TrackedValue(MAX, 0, MAX, 0)` is the
                        catastrophe-signal value per AUD-EF-9.

ERROR BOUND
  Bound category:       precision / accuracy on the safe path; the
                        catastrophe path overrides all error categories
                        with MAX (signal).
  Bound formula:        Safe path (den.value ≠ 0): per REQ-EF-3, the
                        quotient `num / den` propagates per the quotient
                        rule:
                          δ(a/b) ≈ |1/b|·δa + |a/b²|·δb
                        across all four categories.
                        Zero/zero path (point mass / no-torque): the
                        result is exactly 0 with zero error budget (the
                        result is **definitional**, not an approximation;
                        zero rotational inertia plus zero torque has zero
                        angular acceleration as a *boundary value* per the
                        point-mass model documented in the file header).
                        Nonzero/zero path (catastrophe): the bound is
                        irrelevant — the value is set to MAX with
                        precision = MAX, accuracy = 0, measurement = MAX,
                        model = 0 to signal a physically undefined regime
                        per REQ-EF-9 / AUD-EF-9.
  Bound implemented:    Safe path: TrackedValue<T>::operator/ implements
                        the quotient rule per AUD-EF-1.
                        Zero/zero path: `math::exact<T>(0)` returns a
                        zero-error TrackedValue.
                        Nonzero/zero path: explicit construction
                        `TrackedValue<T>(MAX, 0, MAX, 0)` matches the
                        REQ-EF-9 catastrophe-signal contract.
  Bound verdict:        ✓ matched on all three paths — declared and
                        implemented behaviors agree.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form quotient), REQ-EF-9 (safe-div
                        for point-mass / catastrophe).
  AUD-EF applies:       AUD-EF-1 (quotient-rule wiring), AUD-EF-9
                        (catastrophe-signal value verification).
  AUD-MC applies:       AUD-MC-X (H⁻¹H = I check — round-trip
                        `momentum_of` then `acceleration_from_wrench`
                        should recover the input twist when m > 0,
                        I_body ≻ 0).
  Verification test:    tests/test_dynamics — exercise (a) the safe path
                        with m > 0, I_body ≻ 0 and verify H⁻¹H = I to
                        machine precision; (b) the zero/zero path with
                        Inertia::point_mass(m) and a torque-free wrench,
                        verify α = 0; (c) the catastrophe path with
                        Inertia::point_mass(m) and a torqued wrench,
                        verify the result's precision/measurement are
                        MAX.

NOTES
  - The catastrophe-signal pattern `TrackedValue(MAX, 0, MAX, 0)` is the
    REQ-EF-9-defined sentinel: precision-MAX and measurement-MAX dominate
    any downstream `total_error()` so that any consumer immediately
    fails go/no-go; accuracy and model remain 0 because the failure is
    not a model error per se but an out-of-domain operation.
  - The linear leg `F/m` has **no guard**; this is explicitly relying
    on the file-header precondition `mass.value > 0` (line 113). For
    `Inertia::point_mass(m=0)` the linear leg would itself divide by
    zero — but this is forbidden by the precondition, and would surface
    as IEEE inf/NaN at runtime rather than via the safe-div sentinel.
    **Flag (?): consider symmetric safe-div on the linear leg, or
    explicit assertion at function entry.** Not a C-fail because the
    cited theory (Murray-Li-Sastry §4.2 prop.4.5) requires m > 0 for
    H⁻¹ to exist; the precondition is well-defined.
  - The angular safe-div uses `den.value == T(0)` for the bypass test.
    This is **exact equality** on `T`. For double, this means the test
    fires only when I_ii is bit-for-bit zero — the typical point-mass
    case. A near-zero but not-exactly-zero I_ii (e.g. from rounding in
    a derived computation) would go through the normal `num/den` path
    and report large errors from REQ-EF-3 quotient propagation — this
    is the correct behavior (the user **declared** zero by passing 0;
    a near-zero value indicates a real, very large angular
    acceleration).
```

---

## 8 `linear_acceleration_from_force(F)` — point-mass-only inverse mass

```
=== FORMULA AUDIT CARD ===
ID:                     inertia::Inertia::linear_acceleration_from_force
Location:               src/dynamics/inertia.h:144-147
Mathematical statement: Given a body-frame force F = (F_x, F_y, F_z),
                        compute the linear acceleration
                          a_i = F_i / m,   i ∈ {x, y, z}.
                        Newton's second law restricted to translation
                        (the linear leg of H⁻¹).

THEORY
  Underlying theorem:   Newton's second law for a point particle: F = m a,
                        equivalently a = F / m. This is the translational
                        component of the SE(3) inverse inertia operator
                        H⁻¹ for a body whose rotational inertia is being
                        skipped by the propagator (point-mass mode).
  Primary reference:    Newton (1687) *Principia*, Lex II; Goldstein, Poole
                        & Safko (2002) §1.1; for the SE(3) embedding,
                        Murray-Li-Sastry (1994) §4.2.
  Domain of validity:   m > 0. The same precondition as
                        `acceleration_from_wrench`'s linear leg (line 113
                        of the header); applies a fortiori here since
                        only the linear leg is invoked.

METHOD
  Method declared:      Closed-form: three independent scalar quotients
                        F_i / m, packaged in Vector3<T>.
  Method implemented:   `math::Vector3<T>(F.x / mass, F.y / mass, F.z / mass)`
                        — three TrackedValue<T>::operator/ invocations.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision / accuracy propagated from inputs.
  Bound formula:        Per REQ-EF-3, the quotient `F_i / m` propagates per
                        the quotient rule:
                          δ(F_i / m) ≈ |1/m|·δF_i + |F_i/m²|·δm,
                        across all four error categories independently.
  Bound implemented:    TrackedValue<T>::operator/ implements the quotient
                        rule per AUD-EF-1; the Vector3<T> constructor
                        copies the three resulting TrackedValues.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form quotient).
  AUD-EF applies:       AUD-EF-1 (quotient-rule wiring).
  AUD-MC applies:       n/a (this is the linear restriction; the SE(3)
                        algebra audits live with the full
                        `acceleration_from_wrench`).
  Verification test:    tests/test_dynamics — exercise with a known F
                        and a known m, verify a = F/m to machine precision
                        and that the error budget matches the quotient
                        rule.

NOTES
  - Like the linear leg of `acceleration_from_wrench` (§7), this function
    has **no guard** against m = 0; it relies on the file-header
    precondition `mass.value > 0`. For `Inertia::point_mass(m=0)`,
    this would emit IEEE inf/NaN rather than the REQ-EF-9 sentinel.
    **Same (?) flag as §7's NOTE — consider an entry assertion or a
    symmetric safe-div.** Not a C-fail.
  - There is no torque input here, by construction. The function is
    explicitly the point-mass mode partner of
    `acceleration_from_wrench`; callers using
    `Inertia::point_mass(m)` should prefer this function for
    translation-only kinematics so that no angular safe-div fires.
```

---

## File-level verdict

| Function | Theory | Method | Bound | Verdict |
|---|---|---|---|---|
| 1 `Inertia()` default ctor | Murray-Li-Sastry §4.2 (zero ∈ cone) | Closed-form value-init | ✓ zero error | ✓ PASS |
| 2 `Inertia(m, I)` param ctor | Featherstone (2008) §2.13; Murray-Li-Sastry §4.2 | Closed-form copy | ✓ pass-through | ✓ PASS |
| 3 `point_mass(m)` | Murray-Li-Sastry §4.2.1 | Closed-form | ✓ pass-through | ✓ PASS |
| 4 `uniform_sphere(m, r)` | Goldstein §5.3 Table 5.1; ProofWiki | Closed-form `(2/5)·m·r²` | ✓ product rule | ✓ PASS |
| 5 `diagonal(m, I_xx, I_yy, I_zz)` | Murray-Li-Sastry §4.2 (principal axes) | Closed-form copy | ✓ pass-through | ✓ PASS |
| 6 `momentum_of(twist)` | Murray-Li-Sastry §4.2 eq.(4.32); Featherstone §2.13; Selig §5.2 | Closed-form six scalar products | ✓ product rule | ✓ PASS |
| 7 `acceleration_from_wrench(wrench)` | Murray-Li-Sastry §4.2 prop.4.5; REQ-EF-9 | Closed-form quotients + safe-div | ✓ quotient rule + REQ-EF-9 sentinel | ✓ PASS, ? linear-leg unguarded |
| 8 `linear_acceleration_from_force(F)` | Newton Lex II; Murray-Li-Sastry §4.2 | Closed-form three quotients | ✓ quotient rule | ✓ PASS, ? unguarded |

**File verdict: PASS**

- All eight functions are **closed-form** implementations of their cited theorems (no Taylor / Padé / iterative methods are used, so the §5 worked-example category of theory-method mismatch cannot arise here).
- The only nontrivial pattern is **REQ-EF-9 safe-div** in `acceleration_from_wrench`, which is correctly wired to the catastrophe-signal value per AUD-EF-9.
- One **(?) note**: the linear-leg `F/m` in both §7 and §8 is **unguarded** and relies on the documented `mass.value > 0` precondition. This is not a C-fail (the cited theory itself requires m > 0 for H⁻¹ to exist), but a hardening opportunity if symmetric guarding is desired.
- No (⚠) bound issues found.
- No (✗) method-theory mismatches.

## Cross-audit references

- **REQ-EF**: REQ-EF-3 (closed-form propagation), REQ-EF-9 (safe-div semantics).
- **AUD-EF**: AUD-EF-1 (`TrackedValue<T>` arithmetic wiring), AUD-EF-9 (catastrophe-signal verification).
- **AUD-MC**: The H⁻¹H = I round-trip and the symmetry of H are the relevant algebra-axiom tests; AUD-MC-X (the SE(3)-inertia-specific algebra audit) should be authored if not yet in `design/audit/mathematical_correctness.md`.
- **AUD-CC**: The file header at lines 24-27 self-declares conformance to AUD-CC-1, -2, -3, -5, -6, -7, -8, -9, -10, -12, -17, -18 and AUD-EF-1, -7. These are independently audited in their respective documents.

## Open items / future work

1. **(?) Linear-leg guard symmetry** — Consider adding either:
   - An `assert(mass.value > T(0))` at the top of `acceleration_from_wrench` and `linear_acceleration_from_force`, **or**
   - A symmetric `safe_div` on the linear leg matching the angular leg's REQ-EF-9 behavior.
   The current design is consistent with the file-header precondition, but a runtime guard would make the failure mode explicit rather than implicit (IEEE inf/NaN).
2. **Off-diagonal `I_body`** — File header notes the Matrix3<T> generalization is deferred. When added, the audit will gain three or more cards (for the full I_body product, the inverse, and the eigen-decomposition if any). Theory anchor for the full case: Murray-Li-Sastry §4.2 eq.(4.30), Featherstone §2.13.
3. **AUD-MC-X authorship** — If not already in `mathematical_correctness.md`, add an audit for SE(3)-inertia symmetry and `momentum_of` ∘ `acceleration_from_wrench` = I (within the non-degenerate cone).
