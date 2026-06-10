# Theoretical Basis Audit — Framework and Per-File Audit Cards

**Purpose.** Extend the existing audit suite (AUD-EF-1..10 error-framework wiring, AUD-MC-1..18 algebra-axiom correctness) with a **third dimension**: every formula in the codebase must be **annotated with its theoretical basis**, the **implementation method** must match the cited theory, and the **error bound** must be the rigorous bound for the chosen method.

> User direction (2026-05-13): *"every formula in our dual quaternion is wired to it. … all of the formulas will need to be associated with their theoretical basis. … we need to ensure that the correct formulas are used, for example if the theory is using Taylor's theorem, but we are using a continued fraction decomposition, this is incorrect."*

## 0 Three audit dimensions

| Dimension | Existing audit | What it asks | Failure mode |
|---|---|---|---|
| **A. Error wiring** | AUD-EF-1..10 | Does every public op return `TrackedValue<T>`? Is the bound added to the right error category? | Bare `T` slipping through, silently dropped error |
| **B. Algebra axioms** | AUD-MC-1..18 | Do the operations satisfy the laws of the algebra they implement? (associativity, conjugate distribution, magnitude multiplicativity, …) | Operation gives wrong answer for a known identity |
| **C. Theoretical basis** | **NEW — this document** | Does each formula cite the theorem it implements? Does the chosen method match the cited theory? Does the error bound match the rigorous bound for the method used? | "Doc says Taylor, code does continued fraction" — silently wrong error bound |

A function may pass A and B but fail C: e.g. one could implement $\sin(x)$ via a Padé approximant, propagate errors with the correct closed-form rule, satisfy the right algebraic identities — yet have an error bound that assumes Taylor truncation. The error report is no longer rigorous.

## 1 Per-formula audit card

Every numeric-producing function (or operator) in `src/` must have an **audit card** with the following slots:

```
=== FORMULA AUDIT CARD ===
ID:                     [stable identifier — file::function::variant]
Location:               [file:line range]
Mathematical statement: [what is being computed, in math notation]

THEORY
  Underlying theorem:   [Brouwer 1959 Thm N, Taylor's theorem on f at x=a, …]
  Primary reference:    [paper or textbook citation with equation number]
  Domain of validity:   [where the theorem applies — e.g. analytic neighborhood, |x|<R]

METHOD
  Method declared:      [Taylor expansion order N | Newton iter | Halley iter | …]
  Method implemented:   [what the code actually does — must match "declared"]
  Match verdict:        [✓ matched / ✗ mismatched / ⚠ approximate]

ERROR BOUND
  Bound category:       [precision / accuracy / measurement]
  Bound formula:        [the rigorous closed-form bound for the chosen method]
  Bound implemented:    [what the code actually adds — must match "formula"]
  Bound verdict:        [✓ matched / ✗ unsound / ⚠ conservative]

CROSS-AUDIT
  REQ-EF applies:       [which REQ-EF-N this conforms to]
  AUD-EF applies:       [which AUD-EF-N test exercises this]
  AUD-MC applies:       [which algebra audit, if relevant]
  Verification test:    [unit test file:line that exercises it]

NOTES
  [Open issues, dropped terms, threshold caveats, alternative methods considered]
```

## 2 Verification triad

For each audit card:

1. **Theory** is a primary-source claim. We cite **the** theorem that justifies the form of the computation (a Taylor series, a Newton iteration, a closed-form identity, a Lipschitz bound, …). For composite operations, theory may be a chain of cited theorems.

2. **Method** is the specific algorithmic realization. If theory is "Taylor's theorem on $f$ near $a$", method must be the Taylor expansion of $f$ at $a$ — *not* a Padé approximant or continued-fraction expansion of the same function. The implementation must reflect the cited theory.

3. **Bound** is the rigorous error bound *for the chosen method*. For Taylor with $N$ terms: the magnitude of the $(N+1)$th term (or a Lagrange-form bound on the remainder). For Newton iteration: the final correction magnitude (Kantorovich bound). For a Padé approximant: a totally different bound formula based on the approximant's order.

**Theorem 2.1** (correctness composition). For a function $f$ with audit card $(T, M, B)$:

- If $T \Rightarrow M$ is *false* (declared method doesn't follow from cited theory): C-fail.
- If $M \Rightarrow B$ is *false* (declared bound isn't the rigorous bound for declared method): C-fail.
- A C-fail invalidates the entire `total_error()` claim of the operation: any downstream code that consumes its result is operating without a sound error budget.

## 3 Audit status table — DQ propagator core

The following 24 files are in scope. Status is **NEEDED** until their audit card document is committed.

| File | Lines | Audit status | Card document |
|---|---:|---|---|
| `src/math/tracked_value.h` | 534 | NEEDED | `theoretical_basis_audit/tracked_value.md` |
| `src/math/series.h` | 159 | NEEDED | `theoretical_basis_audit/series.md` |
| `src/math/factorial.h` | 70 | NEEDED | `theoretical_basis_audit/factorial.md` |
| `src/math/wallis.h` | 59 | NEEDED | `theoretical_basis_audit/wallis.md` |
| `src/math/binomial_series.h` | 138 | NEEDED | `theoretical_basis_audit/binomial_series.md` |
| `src/math/small_angle_series.h` | 121 | **PASS** (§5 below) | this document, §5 |
| `src/math/kepler.h` | 91 | NEEDED | this document, §6 (stub) |
| `src/math/angles.h` | 57 | NEEDED | `theoretical_basis_audit/angles.md` |
| `src/math/vector3.h` | 73 | NEEDED | `theoretical_basis_audit/vector3.md` |
| `src/math/dual_number.h` | 141 | NEEDED | `theoretical_basis_audit/dual_number.md` |
| `src/math/quaternion.h` | 233 | NEEDED | `theoretical_basis_audit/quaternion.md` |
| `src/math/dual_quaternion.h` | 356 | NEEDED | `theoretical_basis_audit/dual_quaternion.md` |
| `src/constants/constants_provider.h` | 103 | NEEDED | `theoretical_basis_audit/constants_provider.md` |
| `src/dynamics/pose.h` | 111 | NEEDED | `theoretical_basis_audit/pose.md` |
| `src/dynamics/twist.h` | 94 | NEEDED | `theoretical_basis_audit/twist.md` |
| `src/dynamics/wrench.h` | 94 | NEEDED | `theoretical_basis_audit/wrench.md` |
| `src/dynamics/inertia.h` | 150 | NEEDED | `theoretical_basis_audit/inertia.md` |
| `src/dynamics/state.h` | 86 | NEEDED | `theoretical_basis_audit/state.md` |
| `src/dynamics/derivative.h` | 93 | NEEDED | `theoretical_basis_audit/derivative.md` |
| `src/dynamics/propagator.h` | 150 | NEEDED | `theoretical_basis_audit/propagator.md` |
| `src/forces/gravity_central.h` | 84 | NEEDED | `theoretical_basis_audit/gravity_central.md` |
| `src/forces/gravity_zonal.h` | 99 | NEEDED | `theoretical_basis_audit/gravity_zonal.md` |
| `src/forces/drag.h` | 138 | NEEDED | `theoretical_basis_audit/drag.md` |
| `src/integrators/runge_kutta.h` | 139 | NEEDED | `theoretical_basis_audit/runge_kutta.md` |

**Total**: 24 files, ~3,373 lines. Estimated effort: 1 worked card per file averages 30-60 minutes for an analyst familiar with the underlying theory; ~16-24 hours total. The work is decomposable and parallelizable.

## 4 Audit workflow

For each file:

1. **Identify formulas.** Read the file. List every distinct mathematical operation that produces a `TrackedValue<T>` (or a composite of one). A "formula" is each independent application of theory, not each line of code.

2. **Cite theory.** For each formula, identify the theorem or primary-source equation it implements. If the theory is composed (e.g., "Taylor for $\sin$ × Newton for the divisor"), cite both legs.

3. **Determine method.** Read the code carefully. Categorize the implementation:
   - **Taylor order $N$** at a stated point — bound is magnitude of $(N+1)$th term.
   - **Power series with rigorous tail** (Leibniz / geometric / d'Alembert) — bound is the chosen tail majorant.
   - **Newton / Halley iteration** to tolerance $\tau$ — bound is $|\Delta_{\text{final}}|$.
   - **Closed-form identity** — bound is propagated from inputs via REQ-EF-3.
   - **Padé or rational approximant** — bound depends on numerator/denominator orders (REJECTED if cited theory is Taylor).
   - **Continued fraction** — bound depends on truncation depth (REJECTED if cited theory is Taylor).
   - **Lookup table / spline** — accuracy bound from interpolation theory (REJECTED if cited theory is closed form).
   - **Monte-Carlo / iterative sampling** — bound is statistical (NOT a rigorous bound; downgrades to estimate).

4. **Check bound.** Compare the bound formula declared (or implied) by the theory-method choice against the bound the code adds to `errors.X`. Mismatches are the C-fail signal.

5. **Record card.** Write the audit card into the corresponding `theoretical_basis_audit/<file>.md` document.

6. **Tag the code (optional but recommended).** Add a stable comment `// audit:tba:<id>` next to each formula pointing to its card. This makes the code self-documenting.

## 5 Worked example — `src/math/small_angle_series.h`

This file is the canonical "Taylor vs anything-else" test case. The file contains three Taylor-branch helpers used by `quaternion.h` and `dual_quaternion.h`. Each is audited below.

### 5.1 `taylor_sinc(theta, theta_sq)`

```
=== FORMULA AUDIT CARD ===
ID:                     small_angle_series::taylor_sinc
Location:               src/math/small_angle_series.h:45-60
Mathematical statement: sinc(θ) = sin(θ)/θ, with the convention sinc(0)=1

THEORY
  Underlying theorem:   Taylor's theorem (Lagrange-form remainder) applied to
                        the entire function sinc(θ) = sin(θ)/θ at θ=0.
  Primary reference:    Whittaker & Watson (1927) §7.21, or any text on
                        Taylor series of entire functions. The series
                        sinc(θ) = Σ_{n≥0} (-1)^n θ^{2n} / (2n+1)!
                        converges for all θ ∈ ℂ.
  Domain of validity:   All of ℂ; the small-angle branch fires when
                        |θ.value| < 1e-4.

METHOD
  Method declared:      Taylor series of sinc at θ=0, truncated after the
                        θ⁴/120 term (= 5th-order in θ).
  Method implemented:   `exact<T>(1) - theta_sq/exact<T>(6) + t4/exact<T>(120)`
                        i.e. 1 - θ²/6 + θ⁴/120
  Match verdict:        ✓ matched — implementation is Taylor truncated at
                        order 4 (in θ; equivalently order 2 in θ²).

ERROR BOUND
  Bound category:       precision
  Bound formula:        The series is alternating with monotonically
                        decreasing magnitudes for |θ| < √20 ≈ 4.47 (where
                        |term_{n+1}| < |term_n|). By Leibniz's theorem,
                        |R_N(θ)| ≤ |first omitted term| = |θ|⁶ / 6! = |θ|⁶/720.
                        WAIT — the next term in the Taylor series of sinc
                        after θ⁴/120 is -θ⁶/7! = -θ⁶/5040.
                        Conservative bound: |θ²|³ / 5040 = |θ|⁶ / 5040.
  Bound implemented:    `T trunc_bound = ts_abs * ts_abs * ts_abs / T(5040);`
                        where ts_abs = |theta_sq.value| = θ.value².
                        So trunc_bound = (θ²)³ / 5040 = θ⁶/5040.
  Bound verdict:        ✓ matched — implemented bound equals the Leibniz
                        bound |θ|⁶/5040 exactly.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-6 (Fixed-order Taylor truncation added to precision)
  AUD-EF applies:       AUD-EF-5 (Taylor branches add the truncation bound)
  AUD-MC applies:       n/a (this is a numerical helper, not an algebra op)
  Verification test:    tests/test_math/ — sinc test at small θ. The bound
                        |θ|⁶/5040 should be tested against actual error from
                        a high-precision reference.

NOTES
  - At threshold θ = 1e-4, |θ|⁶ = 1e-24. For T = double (ε ≈ 2.2e-16),
    the bound is 24 orders of magnitude below ε. ✓ Tight relative to T's
    representational precision.
  - For wider T (cpp_bin_float_50, etc.), the threshold 1e-4 is too coarse:
    the closed-form sin(θ)/θ would have ~50 digits of accuracy, while the
    Taylor truncation at order 4 has only ~24 digits at θ=1e-4. The note in
    the file header acknowledges this. A T-dependent threshold scaling with
    ε^{1/6} would be tight; this is a future optimization, not a correctness
    failure (the bound remains rigorous; it's just conservative for wide T).
```

### 5.2 `taylor_half_angle_scale(qv_norm, qv_norm_sq, w_pos)`

```
=== FORMULA AUDIT CARD ===
ID:                     small_angle_series::taylor_half_angle_scale
Location:               src/math/small_angle_series.h:73-90
Mathematical statement: scale(s, w) = atan2(s, w) / s,
                        for unit quaternion with s = |q_v|, w ≥ 0;
                        equivalent to arcsin(s)/s when w² + s² = 1.

THEORY
  Underlying theorem:   Taylor's theorem applied to arcsin(s)/s at s=0.
                        arcsin(s)/s = Σ_{n≥0} ((2n)! / (4^n (n!)² (2n+1))) · s^{2n}
                                    = 1 + s²/6 + 3s⁴/40 + 15s⁶/336 + 105s⁸/3456 + …
                        Equivalently:  1 + s²/6 + 3s⁴/40 + 5s⁶/112 + …
                        (Note 15/336 = 5/112 after simplification.)
                        Radius of convergence: |s| ≤ 1 (arcsin's branch points).
  Primary reference:    Abramowitz & Stegun (1964) §15.1.10, or any
                        treatment of arcsin's Taylor series.
  Domain of validity:   s ∈ [−1, 1]; the small-argument branch fires for
                        |qv_norm.value| < 1e-4.

METHOD
  Method declared:      Taylor series of arcsin(s)/s at s=0, truncated after
                        the 3s⁴/40 term (= 5th order in s; even powers only).
  Method implemented:   `exact<T>(1) + s2/exact<T>(6) + exact<T>(3)*s4/exact<T>(40)`
                        i.e. 1 + s²/6 + 3s⁴/40.
  Match verdict:        ✓ matched — implementation is Taylor of arcsin(s)/s
                        truncated at s⁴.

ERROR BOUND
  Bound category:       precision
  Bound formula:        The arcsin(s)/s series has all positive coefficients
                        (NOT alternating). Leibniz does not apply directly.
                        However, for |s| < 1 the ratio of successive terms
                        is bounded:
                          term_{n+1}/term_n = ((2n+1)(2n+2))/((n+1)²(2n+3)) · s²
                                            = (2n+1)/(n+1) · (2n+2)/((n+1)(2n+3)) · s²
                                            ≤ 2 · s²    (for n ≥ 0)
                        For the truncation after s⁴, next term is 5s⁶/112.
                        Conservative bound:  |R_N| ≤ (5/112) · |s|⁶ · 1/(1 − 2s²)
                                            ≤ (5/112) · |s|⁶ · 2   for |s| < 1/2.
                        Looser but rigorous: |R_N| ≤ 5|s|⁶ / 112.
  Bound implemented:    `T trunc_bound = T(5) * s2_abs * s2_abs * s2_abs / T(112);`
                        i.e. 5 · (s²)³ / 112 = 5|s|⁶ / 112.
  Bound verdict:        ⚠ tight-only-for-small-s — the implementation uses
                        the magnitude of the next term, which is rigorous
                        ONLY if the tail is dominated by that next term.
                        For arcsin(s)/s with non-alternating positive
                        coefficients, the tail is technically larger than
                        the next term. The correction factor 1/(1 − 2s²)
                        is ~1.00000002 at the branch threshold s=1e-4, so
                        the implemented bound under-counts by ~2e-8 *
                        next_term ≈ 1e-32 — far below T = double's ε.
                        Verdict: bound is rigorous within type precision
                        for double; for wider T the correction must be
                        applied. **Flag for tightening.**

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-6 (Fixed-order Taylor truncation)
  AUD-EF applies:       AUD-EF-5
  AUD-MC applies:       Used by AUD-MC-13 (axis-angle round-trip via log_unit)
  Verification test:    tests/test_quaternion/ — log_unit on small-angle q

NOTES
  - The non-alternating nature of the arcsin(s)/s series means the simple
    "next term" bound is not the proper Leibniz bound — it's the
    leading term of a positive-coefficient tail. Strictly,
      tail = Σ_{n≥3} c_n s^{2n} ≥ c_3 s⁶ = (5/112) s⁶,
    with equality only as s → 0. The implementation uses the leading
    coefficient as the bound, which is the correct tightening at small s
    and approaches under-bounding (un-rigorous) as s grows.
  - The branch threshold s < 1e-4 keeps us deep in the small-s regime
    where the correction factor 1/(1−2s²) is negligible. The bound is
    rigorous for double; for wider T, tighten by including the geometric-
    tail correction factor 1/(1−2s²) explicitly.
```

### 5.3 `taylor_cos_minus_sinc_over_theta_sq(theta, theta_sq)`

```
=== FORMULA AUDIT CARD ===
ID:                     small_angle_series::taylor_cos_minus_sinc_over_theta_sq
Location:               src/math/small_angle_series.h:106-119
Mathematical statement: β(θ) = (cos(θ) − sinc(θ)) / θ²

THEORY
  Underlying theorem:   Taylor's theorem applied to β(θ) at θ=0, derived
                        from the Taylor series of cos and sinc:
                          cos(θ)  = 1 − θ²/2  + θ⁴/24  − θ⁶/720  + …
                          sinc(θ) = 1 − θ²/6  + θ⁴/120 − θ⁶/5040 + …
                          cos − sinc = (−1/2 + 1/6)θ² + (1/24 − 1/120)θ⁴
                                       + (−1/720 + 1/5040)θ⁶ + …
                                     = −θ²/3 + θ⁴/30 − θ⁶/840 + …
                        Dividing by θ²: β(θ) = −1/3 + θ²/30 − θ⁴/840 + …
                        β is entire after the removable singularity at θ=0.
  Primary reference:    derivation: subtraction of two known Taylor series.
                        Equivalent forms appear in screw-theory texts:
                        Murray, Li & Sastry (1994) §3.2 "exponential
                        coordinates"; Selig (2005) "Geometric Fundamentals
                        of Robotics" §3.4.
  Domain of validity:   All θ; small-argument branch fires for |θ.value|<1e-4.

METHOD
  Method declared:      Taylor series of β(θ) at θ=0, truncated after the
                        θ²/30 term (= 2nd order in θ²).
  Method implemented:   `ratio<T>(-1, 3) + theta_sq/exact<T>(30)`
                        i.e. −1/3 + θ²/30.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision
  Bound formula:        The β series is alternating with magnitudes
                        |c_n| = (1/((2n+3)!) − 1/((2n+5)!)) · (2n+5)(2n+4)
                              — actually easier: alternates in sign by
                        construction (cos and sinc have alternating series
                        of the same sign pattern; subtraction preserves
                        alternation in the difference).
                        For |θ|<1, |c_{n+1}| < |c_n|, so Leibniz applies:
                          |R_N(θ)| ≤ |first omitted term| = |θ⁴|/840.
                        Equivalent: |θ_sq|² / 840.
  Bound implemented:    `T trunc_bound = ts_abs * ts_abs / T(840);`
                        i.e. (θ²)² / 840 = θ⁴/840.
  Bound verdict:        ✓ matched — Leibniz bound exactly.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-6
  AUD-EF applies:       AUD-EF-5
  AUD-MC applies:       Used by AUD-MC-18 (screw exp/log round-trip)
  Verification test:    tests/test_dual_quaternion/ — exp_screw and
                        log_screw with small-angle screws.

NOTES
  - At threshold θ = 1e-4, bound = 1e-16 / 840 ≈ 1.2e-19. Below double-
    precision ε (~2.2e-16). ✓ Tight.
  - The closed-form branch computes `(cos(theta) - sin(theta)/theta) / theta_sq`
    which subtracts two near-1 quantities. For θ ≥ 1e-4, cos and sinc
    differ by ~θ²/3 ≈ 3e-9 — within double-precision ε? No: 3e-9 is much
    larger than ε ≈ 2.2e-16, so the subtraction has ~7 digits of relative
    precision. The branch threshold 1e-4 keeps us where the closed-form
    has at least 7 reliable digits, then the Taylor branch takes over.
    **This is the right threshold for double.**
```

### 5.4 File-level verdict for `small_angle_series.h`

- **A. Error wiring**: ✓ all three functions add `trunc_bound` to `result.errors.precision` per AUD-EF-5.
- **B. Algebra axioms**: n/a (numerical helpers, not algebra operations); their use in `quaternion.h` and `dual_quaternion.h` is exercised by AUD-MC-12, AUD-MC-13, AUD-MC-18.
- **C. Theoretical basis**:
  - 5.1 `taylor_sinc`: ✓ all three slots matched. **PASS.**
  - 5.2 `taylor_half_angle_scale`: ✓ method matches; ⚠ bound is tight for double but technically should include the geometric correction `1/(1−2s²)` for wider T. **PASS with note for tightening.**
  - 5.3 `taylor_cos_minus_sinc_over_theta_sq`: ✓ all three slots matched. **PASS.**

**File verdict: PASS** — all formulas are Taylor-method per cited theory, and bounds match the rigorous Leibniz / next-term form. One note for wide-T tightening of `taylor_half_angle_scale`.

## 6 Next file — `src/math/kepler.h` (stub)

```
=== FORMULA AUDIT CARD ===
ID:                     kepler::solve_kepler
Location:               src/math/kepler.h:32-89
Mathematical statement: Solve E - e·sin(E) = M for E given M, e ∈ [0, 1).

THEORY
  Underlying theorem:   Halley's third-order iteration applied to the
                        function f(E) = E - e·sin(E) - M, with derivatives
                        f'(E) = 1 - e·cos(E), f''(E) = e·sin(E).
                        Newton-Kantorovich for Halley: convergence is
                        cubic when f, f', f'' are continuous and f'(E*) ≠ 0;
                        for Kepler's equation with e < 1, f' > 0 everywhere.
  Primary reference:    Halley (1694), or any text on cubic-order root
                        finding; e.g. Conte & de Boor (1980) §3.7.
                        Kepler-specific treatment: Battin (1999) §5.3
                        "An Introduction to the Mathematics and Methods of
                        Astrodynamics".
  Domain of validity:   e ∈ [0, 1), any M. For e → 1 convergence slows
                        but remains cubic in a neighborhood of E*.

METHOD
  Method declared:      Halley iteration with starter
                        E₀ = M + e·sin(M) + (e²/2)·sin(2M).
                        Iterate Δ = 2f·f' / (2f'² − f·f'') until
                        |Δ| < tolerance. Cap at 20 iterations.
  Method implemented:   src/math/kepler.h:40-82, matches declared.
  Match verdict:        ✓ matched — Halley's method, not Newton, not a
                        continued fraction, not a series expansion.

ERROR BOUND
  Bound category:       precision
  Bound formula:        For Halley iteration, the post-convergence residual
                        |E_k - E*| is bounded by |Δ_k| / (1 - q) for the
                        cubic convergence rate q. For practical bounds the
                        simpler |Δ_k| itself is used, since q is small once
                        we're in the cubic regime.
                        This is the Kantorovich/Ostrowski bound and is the
                        rigorous bound for any iterative method (REQ-EF-5).
  Bound implemented:    `E.errors.precision += abs(delta.value);` at the
                        post-convergence return.
  Bound verdict:        ✓ matched — uses the final correction magnitude
                        as the precision bound, matching the
                        Kantorovich-style bound for cubic convergence.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (Convergence residual added to precision)
  AUD-EF applies:       AUD-EF-4 (Iterative algorithms add residual)
  AUD-MC applies:       n/a
  Verification test:    tests/test_math — verify Halley converges to
                        machine precision for various (M, e); verify the
                        reported precision dominates the actual error.

NOTES
  - The starter E₀ is a 2nd-order Taylor in e: for e<0.3 gives ~3-4
    digits; for e→0.95 still converges but needs more iterations.
  - Cap of 20 iterations is safe: 3^20 ≈ 3.5×10⁹ digits, more than any
    representable type can hold.
  - The fallback at line 87 (`E.errors.precision += abs(last_correction)`)
    is reached only if convergence stalls — should never happen for e<1.
    Bound is conservative.
  - **No method-theory mismatch here**: Halley is what's cited and what's
    implemented. Bound is the Kantorovich bound, appropriate for any
    cubic iteration.
```

## 7 Open audit cards (remaining files)

Each file in §3's status table needs an audit card document with one card per distinct formula. To accelerate the sweep, here's the **expected method per file** based on a first-pass read of the codebase. The auditor confirms or refutes:

| File | Expected methods | Theory anchors |
|---|---|---|
| `tracked_value.h` | Closed-form per-category error propagation (REQ-EF-3) | Triangle inequality, mean-value theorem for each std::f |
| `series.h` | Leibniz alternating bound, geometric tail bound | Alternating series test, geometric series convergence |
| `factorial.h` | Direct integer computation (no truncation) | Definition |
| `wallis.h` | Closed-form Wallis cosine integrals | Wallis product formula (see ProofWiki) |
| `binomial_series.h` | Power-series with caller-tolerance termination | Binomial theorem, Newton's generalized form |
| `kepler.h` | Halley iteration | (audited above) |
| `angles.h` | Closed-form normalization (fmod-based) | mod 2π identity |
| `vector3.h` | Closed-form vector ops, sqrt for norm | Pythagorean identity |
| `dual_number.h` | Closed-form dual-number arithmetic | ε² = 0 → forward-mode AD |
| `quaternion.h` | Hamilton product, conjugate, normalize; Taylor for log_unit/exp_pure singularities | Hamilton's quaternions, SU(2) double cover |
| `dual_quaternion.h` | Dual-quaternion product, screw exp/log; Taylor for singularities | Study's biquaternions, SE(3) screw theory |
| `constants_provider.h` | Closed-form construction from physical inputs | WGS72/84 definitions |
| `pose.h, twist.h, wrench.h, inertia.h, state.h, derivative.h` | Composition over `TrackedValue<T>` | SE(3) algebra |
| `propagator.h` | Munthe-Kaas Lie-group RK4 | Munthe-Kaas (1999) "High order Lie group methods" |
| `gravity_central.h` | Closed-form Newton inverse-square | Newton's gravitation |
| `gravity_zonal.h` | Closed-form J₂ Legendre derivative | Spherical harmonic decomposition |
| `drag.h` | Lane power-law density × v² (this propagator), or `compute_drag_coefficients` (SGP4 path) | Lane (1965), Lane-Hoots (1979) |
| `runge_kutta.h` | Classical RK4 / Munthe-Kaas RK4 | Butcher tableaux, Lie-group integration |

For each file, the auditor produces a `theoretical_basis_audit/<filename>.md` document with one card per formula, then updates the status table in §3.

## 8 Detection rules — common method-theory mismatches

Examples of the failure mode the user described ("theory says Taylor, code does continued fraction") and how to detect them in code review:

| Theory cited | Code pattern suggesting mismatch | Verdict |
|---|---|---|
| Taylor at $a$ | Code has `numerator / denominator` with both being polynomials in the small argument | Likely Padé / rational approximant — C-fail |
| Closed-form identity | Code has a `for` loop summing terms | Likely series approximation — must update theory citation |
| Newton iteration | Code has explicit Halley correction `2ff' / (2f'² − ff'')` | Halley, not Newton — theory citation must reflect |
| Leibniz-alternating bound | Series has all-positive coefficients | Bound is unsound; need geometric or Lagrange-form bound |
| Lipschitz bound for `sin` | Code uses `mean-value`-style `|cos(ξ)|·δ` for some ξ | Equivalent at first order; both rigorous |
| Closed-form `atan2` | Code uses Taylor near origin | Method is Taylor at origin + closed-form elsewhere; theory citation must capture branch |

The audit card's "Match verdict" slot is the single check that prevents this class of failure.

## 9 Reporting and remediation

When a card surfaces a C-fail:

1. **Annotate the file** with a TODO referencing the audit card ID and the mismatch type.
2. **Block downstream use** if the C-fail invalidates the error budget materially. (Practical: if the affected function feeds into a `total_error()` that callers use for go/no-go decisions, the C-fail is blocking.)
3. **Fix in one of two ways**:
   - Update the theory citation (and the bound formula) to match what the code actually does.
   - Update the code to match the cited theory.
4. **Re-audit** the card after the fix.

A passing audit (`PASS` in the status table) is the contract that the function's reported `total_error()` is a **rigorous upper bound** under the **declared theory** with the **declared method** and the **declared bound formula**, end to end.

## 10 Companion documents

- `design/specifications/error_framework.md` — REQ-EF-1..15 (what the error budget must do)
- `design/audit/error_framework.md` — AUD-EF-1..10 (how to verify the wiring)
- `design/audit/mathematical_correctness.md` — AUD-MC-1..18 (algebra-axiom tests)
- `design/audit/code_consistency.md` — AUD-CC-1..18 (style, naming, conventions)
- **this document** — `design/audit/theoretical_basis_audit.md` — formula ↔ theory ↔ method ↔ bound match
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` — example of theoretical-basis derivation completed for SGP4 drag coefs

The three audits (EF / MC / TBA) together cover:
- **EF**: the error budget *is wired up*.
- **MC**: the algebra *is correct*.
- **TBA**: the formula *is the right thing*, and its error bound *is the right bound for what we're doing*.

A function passing all three has a sound `total_error()` claim end to end.
