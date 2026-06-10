# Theoretical Basis Audit — `src/orbit/element_recovery.h`

## 1. File overview

**File**: `src/orbit/element_recovery.h`
**Lines**: 61–147 (`recover_mean_elements` template function)
**Status**: **PASS-with-notes**
**Expected method**: Brouwer–Lyddane "un-Kozai" series reversion + Newton cube-root inversion of Kepler's 3rd law
**Function count**: 1

The file contains one public function, `recover_mean_elements`, but it composes **three distinct theoretical sub-formulas**:

1. **(F1)** Newton cube-root iteration to compute $a_1 = (k_e/n_o)^{2/3}$ (lines 83–96).
2. **(F2)** Brouwer–Lyddane third-order series-reversion polynomial $a_o = a_1 (1 - \delta_1/3 - \delta_1^2 - (134/81)\delta_1^3)$ (lines 105–115).
3. **(F3)** Recovered Brouwer mean motion and SMA via Kepler's 3rd law on $n_o''$ (lines 118–136).

Below: one composite audit card for the function, with explicit sub-slots for each of (F1)/(F2)/(F3).

---

## recover_mean_elements — TLE → Brouwer mean elements

```
=== FORMULA AUDIT CARD ===
ID:                     element_recovery::recover_mean_elements
Location:               src/orbit/element_recovery.h:61-147
Mathematical statement: Given TLE osculating mean motion n₀ (Kozai-mean, with
                        secular J₂ included) and eccentricity e₀, inclination
                        i₀, recover the Brouwer double-primed mean elements
                        (a₀″, n₀″) used by the analytical SGP4 theory:
                          a₁  = (kₑ/n₀)^(2/3)                  (Kepler's 3rd law)
                          δ₁  = (3/2) k₂ (3cos²i₀ − 1) / (a₁² β₀³)
                          a₀  = a₁ (1 − δ₁/3 − δ₁² − (134/81)δ₁³)
                          δ₀  = (3/2) k₂ (3cos²i₀ − 1) / (a₀² β₀³)
                          n₀″ = n₀ / (1 + δ₀)
                          a₀″ = (kₑ/n₀″)^(2/3)
                        with β₀ := √(1 − e₀²).

THEORY
  Underlying theorem:
    (F1) Cube-root inversion of Kepler's third law n²a³ = kₑ² via the
         Newton-Raphson iteration x_{k+1} = (2x_k + c/x_k²)/3 for the
         function f(x) = x³ − c with c = (kₑ/n)². Quadratic convergence
         on the contraction region x > 0 (Conte & de Boor §3.4).
    (F2) Series reversion of the Kozai-to-Brouwer mean-motion relation
         n_o = n_Brouwer (1 + δ(a_Brouwer)) with δ(a) ∝ 1/a², yielding
         the third-order asymptotic series in δ₁:
           a_o/a₁ = 1 − δ₁/3 − δ₁² − (134/81)δ₁³ + O(δ₁⁴).
         This is the Lyddane reorganization of Brouwer's secular Ṁ
         correction (Brouwer 1959 Eq. (39), p. 393), one-step refined
         by re-evaluating δ at a_o (the "δ₀ refinement").
    (F3) Closed-form application of Kepler's 3rd law to convert the
         recovered n₀″ to a₀″ via the same cube-root inversion as (F1).
  Primary reference:
    - Brouwer (1959), "Solution of the Problem of Artificial Satellite
      Theory Without Drag," Astronomical Journal 64, pp. 378–397.
      Specifically §6 (Eq. 39) for the secular Ṁ rate that defines δ.
    - Hoots & Roehrich (1980), Spacetrack Report No. 3, §6 "SGP4 Model,"
      p. 10, Equations (sequence un-numbered in SR3 but reproducible as):
        a₁  = (kₑ/n₀)^(2/3)
        δ₁  = (3/2)·k₂·(3cos²i₀ − 1)/(a₁²·(1 − e₀²)^(3/2))
        a₀  = a₁(1 − (1/3)δ₁ − δ₁² − (134/81)δ₁³)
        δ₀  = (3/2)·k₂·(3cos²i₀ − 1)/(a₀²·(1 − e₀²)^(3/2))
        n₀″ = n₀/(1 + δ₀)
        a₀″ = (kₑ/n₀″)^(2/3)        ← Vallado SGP4.cpp initl() L1097
    - Newton's method for cube root: Conte & de Boor (1980),
      "Elementary Numerical Analysis," §3.4 ("Newton's method").
    - Companion project derivation:
      design/derivations/deprecated/013_element_recovery.md
      (TEMPLATE — symbol table & sources captured; series-reversion
       coefficients 1/3, 1, 134/81 still flagged TODO-DERIVE).
  Domain of validity:
    (F1) c = (kₑ/n)² > 0, which holds for any TLE n₀ > 0. Newton
         iteration converges quadratically from any x₀ > 0 to the
         positive real cube root.
    (F2) Asymptotic series in δ₁. For LEO/MEO/GEO, |δ₁| ≲ 1e-3, so
         truncation residual O(δ₁⁴) ≲ 1e-12. The series is the
         third-order truncation of an infinite reversion; the bound
         on the dropped tail is geometric in δ₁ provided |δ₁| < 1.
         Tightens near equator (3cos²i₀ − 1 → 2); loosens near the
         critical inclination cos²i₀ = 1/3 where δ₁ ≈ 0 trivially.
    (F3) Same as (F1).

METHOD
  Method declared:
    (F1) Newton iteration with starter x₀ = c, cap 30 iterations,
         convergence test |x_{k+1} − x_k| < tolerance. Cube-root
         Newton step: x_new = (2x + c/x²)/3.
    (F2) Direct closed-form evaluation of the third-order polynomial
         in δ₁ — no series sum or iteration, just multiply-add chain.
    (F3) Newton cube-root iteration identical to (F1) but driven by
         the just-computed n₀″ instead of the input n₀.
  Method implemented:
    (F1) Lines 83-96: ke_over_n_sq = (ke/n_input)²; a1 = ke_over_n_sq;
         loop i ∈ [0,30): a1_new = (2·a1 + ke_over_n_sq/a1²)/3;
         residual = |a1_new − a1|; a1 = a1_new; break if residual <
         tolerance and add residual to a1.errors.precision.
         ✓ This is the Newton step for f(x) = x³ − c written
         (using f'(x) = 3x²) as x − f/f' = x − (x³−c)/(3x²) =
         (2x³ + c)/(3x²) = (2x + c/x²)/3.
    (F2) Lines 110-111: a0 = a1 · (1 − del1/3 − del1² − (134/81)·del1³).
         Direct algebraic evaluation; no loop, no series sum.
    (F3) Lines 125-136: identical Newton cube-root structure driven by
         result.n0 (= n₀″), result stored in result.a0.
  Match verdict:
    (F1) ✓ matched — implementation is Newton on f(x) = x³ − c.
    (F2) ✓ matched — implementation is the SR3 §6 polynomial verbatim.
         Note: companion derivation 013_element_recovery.md still has
         the coefficients 1/3, 1, 134/81 flagged TODO-DERIVE; the code
         imports them from SR3 §6 / Vallado.
    (F3) ✓ matched.

ERROR BOUND
  Bound category:
    (F1) precision (iterative residual, REQ-EF-5).
    (F2) accuracy (closed-form propagation, REQ-EF-3) + model
         truncation O(δ₁⁴) for the dropped series tail (REQ-EF-7).
    (F3) precision (same as F1).
  Bound formula:
    (F1) Newton-Kantorovich post-convergence residual:
           |x_k − x*| ≤ |x_{k+1} − x_k| / (1 − q)
         with q = local contraction ratio. For Newton with
         quadratic convergence in the basin of attraction (any
         x₀ > 0 here), q ≪ 1 once converged, so the simpler
         |x_{k+1} − x_k| is a rigorous upper bound at machine ε.
    (F2) Closed-form: ε(a₀) ≤ Σᵢ |∂a₀/∂xᵢ|·ε(xᵢ) for xᵢ ∈
         {a₁, δ₁}. Plus geometric tail bound on the dropped
         series: |R₄| ≤ C·|δ₁|⁴ / (1 − |δ₁|) for the next-order
         term coefficient C (not enumerated here; numerically
         negligible for |δ₁| ≲ 1e-3, giving |R₄| ≲ 1e-12).
    (F3) Same Kantorovich form as (F1).
  Bound implemented:
    (F1) Line 93: `a1.errors.precision = a1.errors.precision + residual;`
         where residual = |a1_new − a1|.value at convergence.
         ✓ matches Kantorovich / final-correction-magnitude form.
    (F2) Lines 110-111 use TrackedValue arithmetic, so accuracy is
         propagated via REQ-EF-3 operator overloads (no explicit
         additive bound from the series tail).
         ⚠ The dropped-tail residual O(δ₁⁴) from truncating the
         series at the cubic term is NOT explicitly added.
    (F3) Line 133: `result.a0.errors.precision = result.a0.errors.precision + res;`
         ✓ matches Kantorovich form, identical to (F1).
  Bound verdict:
    (F1) ✓ matched.
    (F2) ⚠ PASS-with-note — the closed-form Jacobian propagation
         carries through correctly via TrackedValue, but the series
         truncation residual O(δ₁⁴) at the 4th-order coefficient is
         not added to accuracy or model-truncation. For TLE inputs
         this is numerically ≲ 1e-12 (well below the secular drift
         floor) but is formally a REQ-EF-7 gap.
    (F3) ✓ matched.

CROSS-AUDIT
  REQ-EF applies:
    - REQ-EF-3 (closed-form Jacobian propagation) — F2 algebraic chain.
    - REQ-EF-5 (iterative convergence residual to precision) — F1, F3.
    - REQ-EF-7 (model-truncation bound) — F2 dropped O(δ₁⁴) tail
      (currently uncounted, flagged ⚠).
  AUD-EF applies:
    - AUD-EF-4 (iterative algorithms add residual to precision) —
      F1, F3 Newton cube-root loops.
    - AUD-EF-7 (model-truncation bound is added) — F2 (currently
      a partial-fail: bound exists in theory, not in code).
  AUD-MC applies:
    n/a (this function consumes TrackedValue scalars and produces
    a RecoveredElements aggregate of TrackedValue scalars; no
    algebra-axiom obligations beyond REQ-EF wiring).
  Verification test:
    tests/test_orbit/ — recover_mean_elements should be exercised
    against a TLE with known Vallado-reference (a₀″, n₀″) pair.
    Convergence-residual claim from F1/F3 should be tested by
    running tolerance = 1e-16 and verifying the reported precision
    bound brackets the true error against an arbitrary-precision
    reference. Series-truncation claim from F2 should be tested by
    sweeping |δ₁| up to ~0.1 and verifying the implicit residual
    stays below 10·|δ₁|⁴.

NOTES
  1. METHOD/THEORY MATCH (CRITICAL — was a project-historical bug
     site): At lines 120-124 there is an explicit comment block
     warning that the FINAL a₀″ MUST be recomputed via Kepler's
     3rd law (a₀″ = (kₑ/n₀″)^(2/3)) rather than via the algebraic
     shortcut a₀″ = a_o/(1 − δ₀). The shortcut differs at O(δ₀²)
     and was previously found to produce systematic t > 0 position
     errors. Current implementation uses the Kepler-law route at
     line 125, matching Vallado SGP4.cpp initl() line 1097. ✓ This
     is a correct method-theory match and resolves a prior C-fail.
  2. THE 134/81 COEFFICIENT: The SR3 §6 cubic coefficient 134/81
     comes from the third-order Brouwer generating-function
     reversion. Companion derivation 013_element_recovery.md
     captures the symbol table and sources but still has the
     coefficients 1/3, 1, 134/81 marked TODO-DERIVE. The code
     correctly imports the SR3 values; the missing derivation is
     a documentation gap, not a code-correctness issue. **Flag:
     close the loop in 013_element_recovery.md.**
  3. NEWTON STARTER ROBUSTNESS: x₀ = c = (kₑ/n₀)² is roughly
     a^(3) for typical Earth orbits, so the starter is several
     orders of magnitude larger than the true cube root. Newton
     on x³ − c with x₀ > x* converges monotonically downward at
     quadratic rate; the 30-iteration cap is conservative — for
     T = double, convergence to machine ε is reached in ≲ 10
     iterations from this starter. For wide T (cpp_bin_float_50)
     more iterations are needed, but 30 still suffices for ~50
     decimal digits.
  4. DEEP-SPACE / SIMPLE-MODEL FLAGS (lines 143-144): The boolean
     flags `is_deep_space` (period ≥ 225 min) and `use_simple_model`
     (perigee < 220 km) are SR3 §6 / SGP4 conventional thresholds.
     They are not formulas with error bounds; they are dispatch
     decisions. NO audit slot applies to them beyond the precondition
     that period_min and perigee_km themselves carry valid
     TrackedValue error state, which they do (lines 141-142).
  5. THE δ₁ → δ₀ REFINEMENT (lines 105 → 115): The single Picard
     iteration on δ — compute δ₁ at a₁, get a₀ via the polynomial,
     then re-compute δ at a₀ as δ₀ — is itself a half-step series
     reversion (the full reversion would iterate to fixed point).
     SR3 §6 prescribes exactly this two-step pattern. The residual
     of stopping at one refinement is O(δ²·(a_o − a₁)/a₁) ≈ O(δ³),
     which is absorbed into the same series-truncation bound noted
     in F2 ⚠.
  6. DIMENSIONAL CONSISTENCY: All quantities are in SGP4-normalized
     units (Earth radii, rad/min). re_km is multiplied in only at
     the perigee-altitude conversion (line 141). half_J2 is the
     dimensionless k₂ = J₂/2 (= CK2). No unit-conversion errors
     possible inside the recovery iteration itself.
```

---

## 2. File-level verdict

**A. Error wiring**: ✓ all three sub-formulas (F1, F2, F3) return `RecoveredElements<T>` composed of `TrackedValue<T>` scalars. F1/F3 Newton residuals correctly added to `.errors.precision` per REQ-EF-5/AUD-EF-4. F2 closed-form propagates accuracy via overloaded `*` / `+` / `-` operators per REQ-EF-3.

**B. Algebra axioms**: n/a — this is an orbit-mechanics composite, not an algebra operator. The `TrackedValue` operators it composes are exercised by AUD-MC-1..18 at the math-library level.

**C. Theoretical basis**:
- **(F1) Newton cube root**: ✓ all slots matched. Theorem (Newton on $f(x) = x^3 - c$), method (`x_new = (2x + c/x^2)/3`), and bound (final correction magnitude as Kantorovich proxy) all match. **PASS.**
- **(F2) Brouwer–Lyddane series reversion** ($a_o = a_1(1 - \delta_1/3 - \delta_1^2 - (134/81)\delta_1^3)$): ✓ method matches SR3 §6 verbatim; ⚠ the dropped $O(\delta_1^4)$ series-truncation residual is not explicitly added to a REQ-EF-7 model-truncation slot, though it is numerically ≲ 1e-12 for any realistic TLE. **PASS-with-note** (file a small TBA-followup to add the geometric tail bound, and a separate followup to close out the 134/81 derivation in `013_element_recovery.md`).
- **(F3) Kepler-3rd-law recovery of $a_o''$**: ✓ all slots matched. Critically, the code uses Kepler-3rd-law cube-root recomputation (`(ke/n_unkozai)^(2/3)`) instead of the algebraic shortcut $a_o/(1 - \delta_0)$ — this is an *explicit historical bug-fix*, documented in-line at lines 120-124. **PASS.**

**File verdict: PASS-with-notes.** The implementation matches the cited theory (Brouwer 1959 / SR3 §6) and reference (Vallado SGP4.cpp `initl`) verbatim. The Newton iterations carry proper Kantorovich-style precision bounds. The one substantive note is the un-added $O(\delta_1^4)$ series-truncation residual in F2 — numerically negligible for TLEs but formally a REQ-EF-7 gap. A second documentation note: the 134/81 coefficient is still flagged TODO-DERIVE in companion `013_element_recovery.md`.
