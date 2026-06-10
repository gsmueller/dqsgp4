# Theoretical Basis Audit — `src/math/series.h`

**File**: `src/math/series.h` (159 lines)  
**Expected functions**: 4 (alternating_series, geometric_series, horner, series_sqrt)  
**Audit date**: 2026-05-13

---

## Card 1: `alternating_series`

```
=== FORMULA AUDIT CARD ===
ID:                     series::alternating_series
Location:               src/math/series.h:27-69
Mathematical statement: S = Σ_{k=start}^{∞} term(k) where terms alternate in sign
                        and decrease in absolute value, truncated when
                        truncation_bound < tolerance.

THEORY
  Underlying theorem:   Leibniz alternating series test: if a_k is
                        monotonically decreasing in |a_k| and a_k → 0,
                        then |remainder after N terms| ≤ |a_{N+1}|.
                        (Knopp §15)
  Primary reference:    Knopp (1928) Theory and Application of Infinite Series §15.
                        Geometric refinement: tail bound |tail| ≤ |a_N| * r/(1-r)
                        when |a_{k+1}/a_k| ≤ r.
  Domain of validity:   Alternating series with monotonic |term| decrease.
                        Geometric ratio r < 1.

METHOD
  Method declared:      Two modes:
                        (1) Leibniz (convergence_ratio < 0): bound = |term(N)|.
                        (2) Geometric (convergence_ratio ≥ 0): bound = |term(N)| * r/(1-r).
  Method implemented:   Lines 40-69: accumulates sum, computes term_mag = |last_term|.
                        If use_geometric, truncation = term_mag * r/(1-r); else
                        truncation = term_mag. Returns with bound added.
  Match verdict:        ✓ matched — Leibniz when convergence_ratio < 0, geometric
                        tail when convergence_ratio ≥ 0.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Leibniz: |remainder| ≤ |term(N)|.
                        Geometric: |remainder| ≤ |term(N)| * r/(1-r).
  Bound implemented:    Line 54: sum.errors.precision += truncation where
                        truncation = term_mag (Leibniz) or
                        term_mag * r/(1-r) (geometric).
  Bound verdict:        ✓ matched — truncation = Leibniz or geometric tail bound.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-4 (Series tail bound → precision)
  AUD-EF applies:       AUD-EF-5 (Series bounds to precision)
  AUD-MC applies:       n/a (numeric helper)
  Verification test:    tests/test_math — verify bound dominates actual error.

NOTES
  - Leibniz bound is the next omitted term, not a coarser bound (requirement).
  - Geometric mode tightens ~150x when r = 0.007.
```

---

## Card 2: `geometric_series`

```
=== FORMULA AUDIT CARD ===
ID:                     series::geometric_series
Location:               src/math/series.h:77-102
Mathematical statement: S = Σ_{k=start}^{∞} term(k) where |term(k+1)/term(k)| ≤ r < 1,
                        truncated when tail_bound < tolerance.

THEORY
  Underlying theorem:   Geometric series bound: |term(k+1)/term(k)| ≤ r < 1
                        ⟹ |tail| = |Σ_{j=N+1}^{∞} term(j)| ≤ |term(N)| / (1-r).
  Primary reference:    Rudin (1976) Principles of Mathematical Analysis §3.23.
  Domain of validity:   Any series where successive terms satisfy
                        |term(k+1)/term(k)| ≤ r < 1.

METHOD
  Method declared:      Accumulate sum. At each step compute
                        tail_bound = |last_term| * r / (1-r).
                        Truncate when tail_bound < tolerance.
  Method implemented:   Lines 87-97: loop accumulates sum. Line 92 computes
                        tail_bound = abs(last_term) * convergence_ratio / (1 - convergence_ratio).
                        Returns when bound drops below tolerance.
  Match verdict:        ✓ matched — geometric-ratio tail bound, not Leibniz,
                        not a coarser approximation.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Rigorous bound: |tail| ≤ |term(N)| * r / (1-r).
  Bound implemented:    Line 94: sum.errors.precision += tail_bound where
                        tail_bound = abs(last_term) * r / (1-r).
  Bound verdict:        ✓ matched — bound = a_(N+1)/(1-r), NOT a smaller approx.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-4 (Series tail bound → precision)
  AUD-EF applies:       AUD-EF-5 (Series bounds to precision)
  AUD-MC applies:       n/a (numeric helper)
  Verification test:    tests/test_math — verify with varying ratios (0.01, 0.5, 0.9).

NOTES
  - Caller must supply convergence_ratio ≥ |term(k+1)/term(k)| for all k.
  - Non-alternating series valid; geometric bound applies to ratio of magnitudes.
```

---

## Card 3: `horner`

```
=== FORMULA AUDIT CARD ===
ID:                     series::horner
Location:               src/math/series.h:108-118
Mathematical statement: p(x) = c₀ + x·(c₁ + x·(c₂ + ... + x·c_n))
                        (Horner nested multiplication scheme)

THEORY
  Underlying theorem:   Horner (1819): closed-form polynomial evaluation
                        via nested multiplication. Mathematically equivalent
                        to Σ c_i · x^i; reduced operations, better stability.
  Primary reference:    Horner (1819). Golub & van Loan (1996)
                        Matrix Computations §4.6.
  Domain of validity:   All polynomials; any numeric type.

METHOD
  Method declared:      Horner nested form: result ← result * x + c_i,
                        i = n-1 down to 0, starting from c_n.
  Method implemented:   Lines 113-116: loop from i = n_coeffs-2 down to 0,
                        executing result = result * x + coeffs[i].
                        Initial result = coeffs[n_coeffs-1].
  Match verdict:        ✓ matched — Horner nested multiplication exactly.

ERROR BOUND
  Bound category:       precision, accuracy, accuracy_measurement (composed)
  Bound formula:        Error propagates through chain via REQ-EF-3 (closed-form
                        composition). result.errors inherits from inputs.
  Bound implemented:    Line 115: result = result * x + coeffs[i] performs
                        TrackedValue arithmetic, composing error per REQ-EF-3.
  Bound verdict:        ✓ matched — error propagates via REQ-EF-3; no explicit
                        truncation bound. Bound inherited from inputs.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form operations compose error)
  AUD-EF applies:       AUD-EF-3 (Arithmetic ops add to error)
  AUD-MC applies:       n/a (algebraically isomorphic to expansion)
  Verification test:    tests/test_math — verify horner() error matches
                        composed error of explicit summation.

NOTES
  - Horner is closed-form, not series truncation. All error is compositional.
  - Reduced operations (n mult + n add vs. more for direct form).
  - No approximation; polynomial computed exactly (to type precision).
```

---

## Card 4: `series_sqrt`

```
=== FORMULA AUDIT CARD ===
ID:                     series::series_sqrt
Location:               src/math/series.h:132-157
Mathematical statement: Compute s where s² = I via Newton iteration:
                        s_{j+1} = (s_j + I/s_j) / 2, starting from s₀ = 1.

THEORY
  Underlying theorem:   Newton iteration for square root. f(s) = s² - I;
                        s_{j+1} = s_j - f(s_j)/f'(s_j) = (s_j + I/s_j)/2.
                        Quadratic convergence: |e_{j+1}| ≤ C·|e_j|² (Kantorovich).
  Primary reference:    Knopp (1928). Codebase Ch 5, Lemma 5.6.1
                        (geodetic: |U| ~ 0.007 → 4 iterations to δ_p < 1e-40).
  Domain of validity:   I close to 1. |I.value - 1| < 1 guarantees convergence.

METHOD
  Method declared:      Newton iteration (NOT Taylor expansion, NOT quadratic-
                        bound refinement). Start s = 1. Iterate s ← (s + I/s) / 2
                        until |s_new - s| < tolerance.
  Method implemented:   Lines 139-150: loop iterates. Line 140: s_new = (s + I/s)/2.
                        Line 143: correction = |s_new - s|. Returns when
                        correction < tolerance.
  Match verdict:        ✓ matched — Newton iteration exactly. Not Taylor,
                        not quadratic-bound refinement.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Newton iteration: post-convergence error bounded by
                        final correction magnitude |e_final| ≤ |Δ_final|.
                        Kantorovich bound, rigorous for Newton reaching
                        quadratic-convergence basin.
  Bound implemented:    Line 147: s.errors.precision += correction where
                        correction = |s_new - s| at convergence.
  Bound verdict:        ✓ matched — final correction = Kantorovich bound.
                        Matches REQ-EF-5 (iteration residual → precision).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (Iteration residual → precision)
  AUD-EF applies:       AUD-EF-4 (Iterative algorithms add residual)
  AUD-MC applies:       n/a (numeric helper)
  Verification test:    tests/test_math — verify convergence for I in [0.99,1.01];
                        verify reported bound dominates |s_computed² - I|.

NOTES
  - Quadratic: |e_j| ≤ |U|^{2^j} / 2^{j+1} where U = I.value - 1.
  - s₀ = 1 is exact, no initial approximation error.
  - Fallback (line 154) uses equation residual if stalled; conservative, rigorous.
```

---

## File-level verdict for `series.h`

- **A. Error wiring**: ✓ all functions add to precision or propagate via REQ-EF-3.
- **B. Algebra axioms**: n/a (numeric helpers).
- **C. Theoretical basis**:
  - alternating_series: ✓ Leibniz bound = next omitted term.
  - geometric_series: ✓ Bound = a_(N+1)/(1-r).
  - horner: ✓ Closed-form nested multiplication.
  - series_sqrt: ✓ Newton iteration with Kantorovich bound.

**File verdict: PASS** — all functions match cited theory, method, and bound.