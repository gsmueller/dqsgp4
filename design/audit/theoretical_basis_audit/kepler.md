# Theoretical Basis Audit — `src/math/kepler.h`

**Source**: `src/math/kepler.h` (91 lines)
**Functions audited**: 1 (`solve_kepler`)
**Framework**: `design/audit/theoretical_basis_audit.md` (§1 template, §6 stub elevated here)
**Date audited**: 2026-05-13

---

## Overview

`solve_kepler` solves the **classical** Kepler equation

  **E − e·sin(E) = M**

for the eccentric anomaly `E` given the mean anomaly `M ∈ ℝ` and the
eccentricity `e ∈ [0, 1)`. The implementation is a textbook Halley cubic
iteration on the residual

  **f(E) := E − e·sin(E) − M,**     `f'(E) = 1 − e·cos(E),`   `f''(E) = e·sin(E),`

starting from a 2nd-order Taylor-in-e expansion

  **E₀ = M + e·sin(M) + (e²/2)·sin(2M).**

The cubic Halley correction `Δ = 2·f·f' / (2·f'² − f·f'')` is applied with
the standard `E ← E − Δ` update (no sign games: the code's `f` *is* the
residual `R`, unlike the modified-Kepler companion file which uses
`f = −R`). At convergence, `|Δ_final|` is added to `E.errors.precision`
as the Kantorovich/Ostrowski iteration-residual bound per REQ-EF-5.

**Companion file**: `src/orbit/modified_kepler.h` solves the SGP4 *modified*
Kepler equation `x − a_yn·cos(x) + a_xN·sin(x) = U`; see
`theoretical_basis_audit/modified_kepler.md` for that audit. The classical
form audited here is the textbook two-body equation; SGP4 itself does not
call this function directly but uses the modified-Kepler variant.

---

## Card 1 — `solve_kepler`

```
=== FORMULA AUDIT CARD ===
ID:                     kepler::solve_kepler
Location:               src/math/kepler.h:32-89
Mathematical statement: Solve f(E) := E − e·sin(E) − M = 0 for E given
                        M ∈ ℝ and e ∈ [0, 1). Output is the unique
                        eccentric anomaly E satisfying Kepler's equation
                        modulo 2π.

THEORY
  Underlying theorem:   Halley's third-order iteration applied to the
                        smooth scalar residual f(E) = E − e·sin(E) − M.
                        Derivatives:
                          f'(E)  = 1 − e·cos(E),
                          f''(E) = e·sin(E).
                        Cubic-convergence theorem (Halley):
                        For f ∈ C³ with f'(E*) ≠ 0, the iteration
                          E_{k+1} = E_k − 2·f(E_k)·f'(E_k) /
                                          (2·f'(E_k)² − f(E_k)·f''(E_k))
                        converges cubically to the simple root E* in
                        a neighborhood:
                          |E_{k+1} − E*| ≤ C · |E_k − E*|³.
                        For Kepler with e < 1:
                          f'(E) = 1 − e·cos(E) ∈ [1 − e, 1 + e],
                        so f' > 0 globally (bounded away from 0 by 1−e),
                        the simple-root condition holds everywhere, and
                        cubic convergence is global (not just local) once
                        the iterate enters any bounded interval.
                        Existence-uniqueness of the root mod 2π is a
                        classical result (Kepler's equation has a unique
                        solution E ∈ ℝ for each M ∈ ℝ when e < 1; see
                        Battin (1999) §3.2).
  Primary reference:    Halley, E. (1694) "Methodus Nova Accurata et
                        Facilis Inveniendi Radices Aequationum
                        quarumcumque generaliter, sine praevia
                        Reductione," Phil. Trans. Royal Soc. 18,
                        pp. 136-148 — the original cubic-iteration
                        scheme.
                        Modern textbook treatment: Conte & de Boor
                        (1980), "Elementary Numerical Analysis,"
                        §3.7 (Higher-order methods, including Halley's
                        formula Eq. (3.7.5)).
                        Kepler-specific treatment: Battin (1999),
                        "An Introduction to the Mathematics and Methods
                        of Astrodynamics," §5.3 (Solution of Kepler's
                        equation, including Halley's method as a
                        cubic-convergence solver).
                        Starter formula E₀ = M + e·sin(M) + (e²/2)·sin(2M):
                        2nd-order Taylor-in-e expansion of E about e=0;
                        Battin (1999) §5.3 Eq. (5.41) gives the same
                        truncation as a Bessel-function series.
                        Header references `design/derivations/
                        equation_convergence_tricks.md` §5.
  Domain of validity:   e ∈ [0, 1), any M ∈ ℝ. As e → 1⁻ the
                        convergence constant C grows (cusps in f near
                        E = π for e = 1), but for any fixed e < 1 the
                        iteration retains cubic convergence in a
                        neighborhood of E*. Cap of 20 iterations is
                        adequate for any e < 1 and any representable
                        precision.

METHOD
  Method declared:      Halley's third-order iteration with starter
                          E₀ = M + e·sin(M) + (e²/2)·sin(2M),
                        update rule
                          Δ = 2·f·f' / (2·f'² − f·f''),
                          E ← E − Δ,
                        termination at |Δ.value| < tolerance, safety cap
                        of 20 iterations.
                        Sign convention: code's `f` IS the residual
                        R(E) = E − e·sin(E) − M (unlike the modified-
                        Kepler companion which uses f = −R). Therefore
                        the standard Halley step is the textbook form
                        Δ = 2·R·R' / (2·R'² − R·R''),
                        with E ← E − Δ. No sign bookkeeping required.
  Method implemented:   Lines 40-82:
                          line 40-41: starter
                            E = M + e·sin(M) + e²·sin(M+M)/2
                              = M + e·sin(M) + (e²/2)·sin(2M).
                          line 55-82: iteration loop, ≤ 20 iterations:
                            line 56-57: sinE = sin(E), cosE = cos(E)
                            line 60: f   = E − e·sinE − M
                            line 62: fp  = 1 − e·cosE
                            line 64: fpp = e·sinE
                            line 67: num = 2·f·fp
                            line 68: den = 2·fp·fp − f·fpp
                            line 69: Δ   = num/den
                            line 71: E ← E − Δ
                            line 75: if |Δ.value| < tolerance → exit
                            line 79: errors.precision += |Δ.value|
                          line 87 (fallback): same precision bookkeeping
                            on cap exit (with |last_correction.value|).
                        Per iteration: 1 sin, 1 cos, 1 division. Explicit
                        f'' = e·sin(E) term — the cubic-correction
                        signature.
  Match verdict:        ✓ matched — Halley's third-order method on the
                        classical Kepler residual, with starter
                        E₀ = M + e·sin(M) + (e²/2)·sin(2M). NOT Newton
                        (Newton would lack the f·f'' term in the
                        denominator), NOT a Padé approximant of sin or
                        of f⁻¹, NOT a continued-fraction expansion of
                        E in M and e, NOT a series in e truncated at
                        fixed order. The iteration is genuinely cubic,
                        using f, f', f'' explicitly per step.

ERROR BOUND
  Bound category:       precision (REQ-EF-5: iteration-residual
                        contributes to precision)
  Bound formula:        Kantorovich/Ostrowski final-correction bound
                        for any iterative root-finder once in its
                        contraction basin. For Halley in the cubic
                        regime,
                          |E_{k+1} − E*| ≤ C·|E_k − E*|³,
                        and writing q := |E_{k+1} − E*| / |E_k − E*|,
                        the post-step bound is
                          |E_k − E*| ≤ |Δ_k| / (1 − q²)
                        with q → 0 cubically. Conservative leading-order
                        form: |E_k − E*| ≈ |Δ_k|. This is the same
                        Kantorovich-style termination bound used by all
                        iterative methods (REQ-EF-5), tighter for Halley
                        than for Newton at the same iteration count
                        because q² shrinks faster than q.
                        Sharper bound (rarely used): for f' bounded below
                        by 1 − e on the iteration interval, one has
                          |E_k − E*| ≤ |f(E_k)| / (1 − e),
                        which is the *residual-form* bound; the
                        correction-form bound used in code is its first-
                        order proxy and is tight in the cubic regime.
  Bound implemented:    Line 79 (convergence path):
                          E.errors.precision = E.errors.precision
                                              + abs(delta.value);
                        Line 87 (safety-cap fallback path):
                          E.errors.precision = E.errors.precision
                                              + abs(last_correction.value);
                        Both paths add the magnitude of the final
                        Halley correction to the precision budget.
  Bound verdict:        ✓ matched — uses |Δ_final| as the precision
                        contribution, matching the Kantorovich-style
                        bound for cubic convergence. Conservative by a
                        factor of (1 − q²)⁻¹ ≈ 1 (so very tight) once
                        q is small. The same bound used in Card 2 of
                        `modified_kepler.md`.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (Convergence residual added to
                        precision)
  AUD-EF applies:       AUD-EF-4 (Iterative algorithms add residual)
  AUD-MC applies:       n/a (root-finder, not an algebra operation)
  Verification test:    tests/test_math/ — verify Halley converges for
                        sweeps over M and e (e ∈ [0, 0.99], M sampled
                        across [0, 2π]). Verify that
                          - reported precision dominates the true
                            residual |E_returned − E_reference|;
                          - iteration count is ≤ ~6 for double, ≤ ~7
                            for cpp_bin_float_50, well under the
                            20-iteration cap;
                          - the result reduces to E = M when e = 0
                            (within one iteration).

NOTES
  - The starter E₀ = M + e·sin(M) + (e²/2)·sin(2M) is the 2nd-order
    Taylor-in-e truncation of the Bessel-function expansion of E in
    Fourier modes of M. Battin (1999) §5.3 / Eq. (5.41) and Brouwer
    & Clemence (1961) §6.5 both give this expansion. For e ≤ 0.3
    the starter is accurate to ~3-4 digits; for e → 0.95 it provides
    only ~1 digit but Halley's cubic convergence still drives the
    iterate to machine precision in 6-7 steps.
  - Cap of 20 iterations is *very* safe: at cubic convergence, 3²⁰
    ≈ 3.5 × 10⁹ digits of headroom — more than any representable
    type can hold by ~7 orders of magnitude. For double precision
    (16 digits) Halley typically terminates in 3-4 iterations.
  - **Fallback path** (line 87): unlike the modified-Kepler companion
    file (whose cap-exit path silently returns x with no residual
    added), this `solve_kepler` adds `|last_correction.value|` to
    `errors.precision` in the cap-exit branch as well. ✓ The error
    budget is sound even if convergence stalls. The header note at
    line 47 ("we don't precompute — we just iterate until convergence")
    is supported by this fallback bookkeeping.
  - **No method-theory mismatch**: Halley is what's cited (Halley 1694
    via Battin §5.3) and what's implemented (Δ = 2·f·f' / (2·f'² − f·f''),
    cubic correction with explicit f''). Bound is the Kantorovich
    bound, appropriate for any cubic-order iterative method.
  - Sign convention is the textbook one (f = R, E ← E − Δ). This is
    the simpler convention; the companion modified-Kepler file uses
    f = −R, x ← x + Δ which is equivalent but requires sign-tracking
    care (see Card 2 of `modified_kepler.md` for the sign audit).
  - Halley is the m = 2 case of Householder's iteration family
    (m = 1 → Newton, m = 2 → Halley, m = 3 → Householder quartic).
    The modified-Kepler companion header advertises Householder as
    an option; this classical file commits to Halley as the only
    solver, matching the docstring exactly.
```

**Verdict**: ✓ **PASS**

---

## File-level verdict

| Dimension | Verdict | Notes |
|---|---|---|
| **A. Error wiring** | ✓ | `solve_kepler` adds `|Δ_final|` (convergence path) and `|last_correction|` (fallback path) to `errors.precision` per REQ-EF-5 / AUD-EF-4. |
| **B. Algebra axioms** | n/a | Numerical root-finder, not an algebra operation. |
| **C. Theoretical basis** | ✓ PASS | Halley on f(E) = E − e·sin(E) − M; uses f, f', f''; cubic-convergence bound = Kantorovich \|Δ\|. Match verdict: ✓ across THEORY → METHOD → BOUND. |

**Overall status**: ✓ **PASS** — Halley is what's declared in the header,
what the math derivation cites (Halley 1694 via Battin §5.3 and
Conte & de Boor §3.7), and what the code implements (Δ = 2·f·f' /
(2·f'² − f·f''), cubic correction with explicit f''-term). The
precision bound `|Δ_final|` is the rigorous Kantorovich/Ostrowski
termination bound for any iterative method (REQ-EF-5), and is tight
in the cubic regime.

## Constraints summary

| Constraint | Status |
|---|---|
| Theory cited matches Halley (not Newton, not Padé, not continued fraction, not series) | ✓ |
| Method implemented matches Halley's cubic Δ = 2·f·f'/(2·f'² − f·f'') | ✓ |
| Precision bound is Kantorovich/final-correction `|Δ_final|` (REQ-EF-5) | ✓ |
| Bound added on both convergence and safety-cap paths | ✓ |
| Sign convention (f = R, E ← E − Δ) is the textbook form | ✓ |
| Starter formula (2nd-order Taylor-in-e) traces to Battin §5.3 Eq. (5.41) | ✓ |
| Domain of validity (e ∈ [0,1)) matches Kepler's equation existence-uniqueness | ✓ |
| Safety cap (20 iterations) is far above any representable precision | ✓ |

No `?`, `⚠`, or `✗` constraints.

## Cross-references

- Framework: `design/audit/theoretical_basis_audit.md` (§1 template, §6
  stub — this document elevates that stub to a full audit per §3
  status table).
- Companion file (SGP4 modified Kepler): `src/orbit/modified_kepler.h`
  → audit `theoretical_basis_audit/modified_kepler.md`. Same Halley
  theorem applied to a different (rotated, eccentricity-vector)
  residual; sign convention differs.
- Derivation reference (per header line 13):
  `design/derivations/equation_convergence_tricks.md` §5.
- REQ-EF: `design/specifications/error_framework.md` REQ-EF-5
  (iteration residual → precision).
- AUD-EF: `design/audit/error_framework.md` AUD-EF-4 (iterative
  algorithms add residual).
- Primary sources:
  - Halley, E. (1694) — original cubic-iteration paper.
  - Battin (1999), *An Introduction to the Mathematics and Methods
    of Astrodynamics*, §5.3 — Kepler-specific Halley treatment.
  - Conte & de Boor (1980), *Elementary Numerical Analysis*, §3.7,
    Eq. (3.7.5) — general Halley formula.
