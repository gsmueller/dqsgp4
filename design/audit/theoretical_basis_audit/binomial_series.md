# Theoretical Basis Audit — `src/math/binomial_series.h`

**File**: `src/math/binomial_series.h` (138 lines)  
**Expected functions**: 2 (`make_binomial_evaluator`, `geodetic_binomial_coefficient`)  
**Audit date**: 2026-05-13

---

## Card 1: `make_binomial_evaluator`

```
=== FORMULA AUDIT CARD ===
ID:                     binomial_series::make_binomial_evaluator
Location:               src/math/binomial_series.h:47-96
Mathematical statement: f(x) = (1+x)^α  for |x| < 1, α ∈ ℝ,
                        evaluated as
                          f(x) = Σ_{k=0}^{N} C(α, k) x^k + R_N(x),
                        with N chosen at call time so that the tail bound
                        |R_N(x)| < tolerance.
                        C(α, k) is the generalized binomial coefficient
                        C(α, 0) = 1,  C(α, k) = C(α, k−1)·(α − k + 1)/k.

THEORY
  Underlying theorem:   Newton's generalized binomial theorem
                        (Newton, 1665; published 1736 in *De Methodis*):
                        for α ∈ ℝ and |x| < 1,
                          (1+x)^α = Σ_{k=0}^{∞} C(α, k) x^k,
                        converges absolutely on the open disk |x| < 1.
                        Coefficient recurrence
                        C(α, k) = C(α, k−1)·(α−k+1)/k is the definition
                        via the falling factorial / k! (Ch 4, Thm 4.10.2).
  Primary reference:    Newton (1665) generalized binomial theorem;
                        modern treatment: Whittaker & Watson (1927)
                        *A Course of Modern Analysis* §4.2 "Binomial
                        theorem". Convergence + remainder:
                        Knopp (1928) *Theory and Application of Infinite
                        Series* §57 ("The binomial series for arbitrary
                        exponent").
                        Project-internal anchor: Ch 4, Theorem 4.10.1
                        (binomial series convergence and remainder);
                        Ch 4, Theorem 4.10.2 (exact rational coefficients
                        for α = p/q).
  Domain of validity:   x ∈ ℝ (or ℂ) with |x| < 1; α ∈ ℝ unrestricted.
                        The series fails to converge for |x| ≥ 1 unless
                        α ∈ ℤ_{≥0} (which terminates).

METHOD
  Method declared:      Power-series evaluation of the Newton binomial
                        series, with caller-tolerance termination via a
                        geometric tail bound (Ch 1, Thm 1.5.3 cited in
                        the source header):
                          |R_N| ≤ |C(α, N+1) x^{N+1}| / (1 − |x|)
                                = |a_{N+1}| / (1 − |x|).
                        Coefficients computed incrementally by the
                        recurrence C(α, k) = C(α, k−1)·(α−k+1)/k —
                        exact rational arithmetic at the algebra level
                        (Ch 4, Thm 4.10.2).
  Method implemented:   Lines 72-90: initialises sum = 1, x_power = 1,
                        coeff = 1; for k = 1..9999:
                          coeff   *= (alpha − (k−1)) / k        (line 78)
                          x_power *= x                          (line 79)
                          term     = coeff * x_power            (line 81)
                          sum     += term                       (line 82)
                          tail_bound = |term| · |x| / (1 − |x|) (line 85)
                          if tail_bound < tolerance: stop,
                            adding tail_bound to sum.errors.precision.
                        For |x| ≥ 1 (line 63) the function returns
                        TrackedValue(0, 0, |x|, 0) — i.e., a precision
                        error equal to |x| signalling divergence.
  Match verdict:        ✓ matched — implementation is the Newton binomial
                        power series with d'Alembert / geometric ratio
                        tail bound. NOT a Padé approximant, NOT a
                        continued fraction, NOT a Chebyshev expansion.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Geometric / d'Alembert tail bound. For the
                        Newton series the ratio of successive terms is
                          |a_{k+1}/a_k|
                            = |(α − k)/(k+1)| · |x|
                            → |x| as k → ∞.
                        Hence for sufficiently large k the tail is
                        dominated by a geometric series with ratio |x|:
                          |R_N| = |Σ_{k=N+1}^{∞} a_k|
                                ≤ |a_{N+1}| · Σ_{j=0}^{∞} |x|^j
                                = |a_{N+1}| / (1 − |x|).
                        This is the standard ratio-test tail bound
                        (Knopp §57; Rudin (1976) §3.23 for geometric
                        majorant). Equivalent to the form printed in
                        the header comment line 43:
                          |R_N| ≤ |C(α, N+1)| · |x|^{N+1} / (1 − |x|).
  Bound implemented:    Line 85: `tail_bound = |term| · |x| / (1 − |x|)`
                        where `term = C(α, k) · x^k` at the current k.
                        So tail_bound = |C(α, k) x^k| · |x| / (1 − |x|)
                                      = |C(α, k) x^{k+1}| / (1 − |x|)
                                      = |a_{k+1}| / (1 − |x|).
                        Termination check (line 86) uses this value;
                        on accept, line 87 adds it to
                        `sum.errors.precision`.
  Bound verdict:        ⚠ tight-eventually — the bound |a_{N+1}|/(1−|x|)
                        rigorously dominates the geometric tail ONCE
                        |a_{k+1}/a_k| ≤ |x| for all k ≥ N. From the
                        ratio |(α−k)/(k+1)|·|x|, this holds when
                        |α − k| ≤ k + 1, i.e. k ≥ (α − 1)/2 (for α > 0)
                        or k ≥ −(α + 1)/2 (for α < 0). For small
                        |α| ≤ 2 (the typical geodetic case, α = −3/2 or
                        −1/2), this is satisfied from k = 1.
                        For very large |α| the early terms may exceed
                        the geometric bound until k passes (|α|−1)/2;
                        in that regime the bound is NOT yet rigorous
                        and the caller should provide a tolerance
                        commensurate with the early-k behaviour.
                        **For all uses inside this codebase** (geodetic
                        meridian-arc series, J_{2n} ellipsoidal series
                        — α ∈ {−3/2, −1/2, 1/2, 3/2}, |x| ≈ e² ≤ 0.007)
                        the bound is rigorous from k = 1 onward.
                        Flag for documentation: state the rigorous
                        regime |α| ≤ 1 (or k ≥ (|α|+1)/2) in the header.

                        The divergence-signalling path on line 63-66
                        (|x| ≥ 1) returns a non-rigorous "value = 0
                        with precision = |x|" sentinel. This is not a
                        true error bound; it is a failure indicator.
                        Documented behaviour, but flag for the caller:
                        downstream consumers must check the magnitude
                        of precision against value before trusting the
                        result.

                        Non-convergence fallback on line 92-94 (loop
                        exhausted after 9999 terms) adds |x_power| =
                        |x|^{9999} as a conservative precision bound.
                        For |x| < 1 this is below any practical
                        tolerance; for |x| ≥ 1 the loop is never
                        entered because of the line-63 early return.
                        ✓ Conservative but never reached in practice.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-4 (Series tail bound → precision)
                        REQ-EF-6 (Caller-tolerance termination)
  AUD-EF applies:       AUD-EF-5 (Series bounds added to precision)
  AUD-MC applies:       n/a (numerical helper; not an algebra op)
  Verification test:    tests/test_math — verify the bound exceeds
                        actual error for representative α ∈ {−3/2, −1/2,
                        1/2, 3/2} and |x| ∈ {0.001, 0.007, 0.5, 0.9}.
                        Convergence test: integer α (e.g. α = 3) must
                        terminate at k = α + 1 with C(α, α+1) = 0 — the
                        bound becomes exactly 0 and accuracy is exact.

NOTES
  - The closure-pattern documented in the header (lines 36-46) is the
    init-time capture pattern from Ch 3, §3.4: α and tolerance are
    fixed at creation; x varies per call. Coefficients can NOT be
    precomputed because N depends on |x| — confirmed by line 76's
    per-call loop.
  - The use of `exact<T>(1)` for sum/x_power/coeff initialisation
    means the recurrence starts with TrackedValues whose precision =
    accuracy = measurement = 0; any error in the final result is
    fully attributable to (a) the truncation bound on line 87 and
    (b) propagated errors in `alpha` and `x`. ✓ Sound wiring.
  - The hard cap `k < 10000` on line 76 is a safety guard; for
    |x| < 1 the geometric bound forces termination in
    N ≈ log(tolerance) / log(|x|) terms. For |x| = 0.007 (e² for
    WGS84), tolerance = 1e-40 requires N ≈ 19 terms; the cap of
    10000 is generous by ~500×.
  - **Header comment** at lines 36-46 already documents the geometric
    tail bound formula. ✓ Self-documenting; no AUD-CC mismatch.
```

---

## Card 2: `geodetic_binomial_coefficient`

```
=== FORMULA AUDIT CARD ===
ID:                     binomial_series::geodetic_binomial_coefficient
Location:               src/math/binomial_series.h:113-136
Mathematical statement: For k ∈ ℤ_{≥0}, α ∈ ℝ, return the k-th
                        coefficient of one of two integrated geodetic
                        binomial series:

                        Case A (cosine_weight = false, w = 1):
                          a_k = C(α, k) · I_k^{(0)}
                          where I_k^{(0)} = ∫_0^{π/2} sin^{2k}(φ) dφ
                                          ÷ (π/2)
                                         = (2k−1)!! / (2k)!!     (k ≥ 1)
                          and I_0^{(0)} = 1 (Wallis even-power identity,
                          normalised by π/2).
                        Case B (cosine_weight = true, w = cos φ):
                          a_k = C(α, k) · I_k^{(c)}
                          where I_k^{(c)} = ∫_0^{π/2} sin^{2k}(φ) cos(φ) dφ
                                          = 1 / (2k + 1).
                        Sign factors (−1)^k from the argument
                        (−e^2 sin^2 φ)^k are deliberately NOT included —
                        the caller multiplies by (−e²)^k or (e²)^k as
                        appropriate (header lines 119-122).

THEORY
  Underlying theorem:   Composition of two theorems (the "three-stage"
                        derivation pattern of Ch 5, Theorem 5.6.2):
                          (Stage 1) Newton's generalized binomial
                          theorem (Card 1 above) supplies C(α, k) as
                          the coefficient of x^k in (1 + x)^α with
                          x = −e² sin² φ.
                          (Stage 2) Term-by-term integration of an
                          absolutely-convergent power series over
                          φ ∈ [0, π/2], with weight w(φ) ∈ {1, cos φ}:
                            ∫_0^{π/2} (1 − e² sin² φ)^α w(φ) dφ
                              = Σ_{k=0}^{∞} (−e²)^k C(α, k) ·
                                ∫_0^{π/2} sin^{2k}(φ) w(φ) dφ.
                          For w = 1: the φ-integral is the even-power
                          Wallis integral W_{2k} = (2k−1)!!/(2k)!! · π/2
                          (here normalised by π/2; "integrated form" of
                          Heiskanen & Moritz (1967) §2-14 meridian-arc
                          identity).
                          For w = cos φ: the φ-integral is the
                          elementary identity
                          ∫_0^{π/2} sin^{2k}(φ) cos(φ) dφ = 1/(2k+1)
                          (substitution u = sin φ).
                        The right of Fubini / dominated convergence to
                        interchange sum and integral is the dominated-
                        convergence theorem applied to the geometric
                        majorant of the binomial series (Card 1's
                        domain |x| < 1, with x = −e² sin² φ and
                        |x| ≤ e² < 1).
  Primary reference:    Stage 1: Newton (1665) / Whittaker & Watson
                        §4.2 (as Card 1).
                        Stage 2 (Wallis even-power identity):
                        Wallis (1656) *Arithmetica Infinitorum*;
                        modern statement Whittaker & Watson §12.31,
                        or ProofWiki "Wallis's product / integral
                        form". Project-internal anchor:
                        Ch 5, Thm 5.6.2 "Geodetic three-stage
                        coefficient" (header line 24).
                        Stage 2 (cos-weighted): elementary integration
                        by substitution u = sin φ; appears as
                        Heiskanen-Moritz (1967) *Physical Geodesy*
                        §2-13/2-14 in the integrated form for
                        meridian-arc and gravity-anomaly series.
                        Heiskanen-Moritz § identification: "?" —
                        the precise §-number depends on edition; the
                        derivation is canonical to ellipsoidal
                        geodesy and not specific to one author. Flag
                        for verification when Ch 5 Theorem 5.6.2 is
                        finalised in `design/derivations/ch05_series.md`.
  Domain of validity:   k ∈ ℤ_{≥0}; α ∈ ℝ (typically α ∈ {−3/2, −1/2,
                        1/2, 3/2} for ellipsoidal series). The
                        underlying series converges for |e²| < 1
                        (always true for an ellipsoid; for WGS84,
                        e² ≈ 0.00669).

METHOD
  Method declared:      Closed-form rational construction:
                          Stage 1: C(α, k) via the generalized binomial
                          coefficient (factorial.h::generalized_binomial,
                          which uses falling_factorial / k!).
                          Stage 2A (cosine_weight = true):
                            return binom * sin_power_cos_integral<T>(k)
                            = C(α, k) / (2k + 1).
                          Stage 2B (cosine_weight = false):
                            k = 0  → return C(α, 0) = 1.
                            k ≥ 1  → return
                              C(α, k) · (2k − 1)!! / (2k)!!.
                        No truncation; no iteration; pure algebra over
                        TrackedValue<T>.
  Method implemented:   Line 124: `binom = generalized_binomial(alpha, k)`
                                  — Stage 1 closed form.
                        Line 127: `if (cosine_weight)`
                                  branches the Stage 2 multiplier.
                        Line 129: cosine-weighted branch returns
                                  `binom * sin_power_cos_integral<T>(k)`
                                  = `binom * 1/(2k+1)` (wallis.h:55-57).
                        Line 133: even-Wallis branch, k = 0 returns
                                  `binom` (since W_0/(π/2) = 1, header
                                  line 132).
                        Line 134: even-Wallis branch, k ≥ 1 returns
                                  `binom * double_factorial(2k−1)
                                          / double_factorial(2k)`
                                  = `binom · (2k−1)!! / (2k)!!`.
  Match verdict:        ✓ matched — implementation is the literal
                        closed-form rational composition of Stage 1
                        and Stage 2. NOT an approximation; NOT an
                        iterative method.

ERROR BOUND
  Bound category:       precision (only)
  Bound formula:        Closed-form propagation only (REQ-EF-3): the
                        result inherits the precision/accuracy/
                        measurement components of `alpha` through the
                        falling-factorial multiplications inside
                        `generalized_binomial`, plus zero new bound
                        contribution from the integer-valued double
                        factorials and from `1/(2k+1)`.
                        Specifically:
                          - C(α, k) is computed via k multiplications
                            and one division, each of which propagates
                            inputs' errors via REQ-EF-3.
                          - (2k−1)!! and (2k)!! are exact integers in
                            T; representational precision error
                            (round-to-T) is the only contribution and
                            is captured at the `*` operator level.
                          - 1/(2k+1) is `ratio<T>(1, 2k+1)`, which is
                            an exact rational converted to T; the
                            representational precision is captured by
                            `ratio<T>` (factorial.h / tracked_value.h).
                        There is NO truncation here — Stage 2 is a
                        closed-form identity, not a series.
                        The "error bound for this method" is therefore
                        just the propagated bound from inputs, no
                        additional truncation term is required.
  Bound implemented:    The function performs the multiplications
                        directly on TrackedValue<T> objects (line 124,
                        129, 134). Errors propagate via the overloaded
                        operators in tracked_value.h per REQ-EF-3.
                        No extra `errors.precision +=` term is added,
                        because no method-specific truncation occurs.
  Bound verdict:        ✓ matched — closed-form identity, errors
                        propagate from inputs via the standard
                        TrackedValue arithmetic. No truncation bound
                        applies because there is no truncation.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form propagation)
  AUD-EF applies:       AUD-EF-3 (Operator-level propagation correctness)
  AUD-MC applies:       n/a (numerical helper)
  Verification test:    tests/test_math — verify the returned
                        coefficient against hand-computed values for
                        small k. Example: α = −3/2, k = 1:
                          C(−3/2, 1) = −3/2;
                          Stage 2A (cos-weighted): −3/2 · 1/3 = −1/2.
                          Stage 2B (Wallis): −3/2 · 1!!/2!!
                                           = −3/2 · 1/2 = −3/4.
                        Verify exact rational result holds when T is
                        cpp_rational; verify ε-level rounding when T
                        is double.

NOTES
  - Header lines 8-16 correctly describe the function's role as the
    Stage-1 × Stage-2 building block of Ch 5 §5.6's three-stage
    geodetic series pattern. The (-1)^k sign factor from the
    argument is deliberately NOT included in this function (header
    line 119-122), consistent with the "Stage 1 × Stage 2" splitting.
  - The Ch 5 Theorem 5.6.2 citation in the header (line 23) is
    project-internal; the underlying mathematics is the
    composition of Newton's binomial theorem + Wallis even-power
    integral + dominated convergence (or the elementary
    cos-weighted substitution). Heiskanen-Moritz (1967) is a
    plausible textbook anchor for the integrated form but the
    precise §-number is "?" pending finalisation of Ch 5.
  - The function takes `alpha` as a `TrackedValue<T>` (not just T),
    so non-integer α with its own error budget (e.g., α computed
    from f = (a − b)/a with measurement uncertainty in a, b)
    propagates correctly. ✓
  - **No method-theory mismatch**: this is a closed-form rational
    construction, declared and implemented as such. The only
    "approximation" present is the representational rounding of
    the rational result into type T, which is the universal
    REQ-EF-3 closed-form bound.
```

---

## File-level verdict — `binomial_series.h`

- **A. Error wiring**: ✓ both functions add precision bounds (or propagate them) per REQ-EF-3/4/6 and AUD-EF-3/5.
- **B. Algebra axioms**: n/a (numerical helpers; not algebra operations).
- **C. Theoretical basis**:
  - Card 1 `make_binomial_evaluator`: ✓ theory (Newton's binomial theorem) matches method (power series with d'Alembert/geometric tail bound). ⚠ bound is rigorous only once the recurrence ratio falls below |x|; for |α| ≤ 2 (all known users), this holds from k = 1. Flag header documentation to state the rigorous regime explicitly. **PASS with note.**
  - Card 2 `geodetic_binomial_coefficient`: ✓ all three slots matched. Closed-form rational; no truncation. Heiskanen-Moritz §-anchor marked "?" pending Ch 5 Theorem 5.6.2 finalisation in `design/derivations/ch05_series.md`. **PASS.**

**File verdict: PASS** — both formulas correctly match Newton's generalized binomial theorem (Card 1) and the three-stage geodetic coefficient construction (Card 2). Two flags: (i) document the rigorous-regime condition on the geometric tail bound in `make_binomial_evaluator`'s header; (ii) verify the Heiskanen-Moritz §-number once Ch 5 Theorem 5.6.2 is finalised.
