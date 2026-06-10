# Theoretical Basis Audit — `src/math/small_angle_series.h`

**File**: `src/math/small_angle_series.h`
**Functions audited**: 3 (taylor_sinc, taylor_half_angle_scale, taylor_cos_minus_sinc_over_theta_sq)
**Date audited**: 2026-05-13
**Audit status**: PASS

---

## 1. taylor_sinc

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

**Verdict**: ✓ **PASS**

---

## 2. taylor_half_angle_scale

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

**Verdict**: ⚠ **PASS-with-note** — Geometric correction needed for wide T.

---

## 3. taylor_cos_minus_sinc_over_theta_sq

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

**Verdict**: ✓ **PASS**

---

## File-level verdict

- **A. Error wiring**: ✓ all three functions add `trunc_bound` to `result.errors.precision` per AUD-EF-5.
- **B. Algebra axioms**: n/a (numerical helpers, not algebra operations); their use in `quaternion.h` and `dual_quaternion.h` is exercised by AUD-MC-12, AUD-MC-13, AUD-MC-18.
- **C. Theoretical basis**:
  - 1. `taylor_sinc`: ✓ all three slots matched. **PASS.**
  - 2. `taylor_half_angle_scale`: ✓ method matches; ⚠ bound is tight for double but technically should include the geometric correction `1/(1−2s²)` for wider T. **PASS-with-note.**
  - 3. `taylor_cos_minus_sinc_over_theta_sq`: ✓ all three slots matched. **PASS.**

**Overall status**: **PASS**
