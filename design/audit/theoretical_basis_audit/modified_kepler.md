# Theoretical Basis Audit — `src/orbit/modified_kepler.h`

**Source**: `src/orbit/modified_kepler.h` (99 lines)
**Functions audited**: 2 (`solve_kepler_newton`, `solve_kepler_halley`)
**Framework**: `design/audit/theoretical_basis_audit.md`

---

## Overview

The SGP4 *modified* Kepler equation differs from the classical Kepler form
`E − e·sin(E) = M`. SGP4 instead solves, in the rotated variable
`x = E + ω`:

  **x − a_yn·cos(x) + a_xN·sin(x) = U**

with `a_xN = e·cos(ω)`, `a_yn = e·sin(ω) + (long-period correction)`, and
`U` the mean argument of latitude. This formulation is the Lyddane
eccentricity-vector form: it avoids the singularity at `e = 0` (where `ω`
is undefined) by working with the components `(a_xN, a_yn)` of the
eccentricity vector rather than `(e, ω)` separately.

Defining the residual
  **R(x) := x − a_yn·cos(x) + a_xN·sin(x) − U,**
both solvers seek the root `R(x*) = 0`. In code the file uses
`f := U − a_yn·cos(x) + a_xN·sin(x) − x = −R(x)`; signs cascade through
`f'` and the Newton/Halley updates so that `x ← x + Δ` (rather than
`x ← x − Δ`) is the correct step. This sign convention is verified
explicitly below.

---

## Card 1 — `solve_kepler_newton`

```
=== FORMULA AUDIT CARD ===
ID:                     modified_kepler::solve_kepler_newton
Location:               src/orbit/modified_kepler.h:40-61
Mathematical statement: Solve R(x) := x − a_yn·cos(x) + a_xN·sin(x) − U = 0
                        for x given (a_xN, a_yn, U). The unknown x = E + ω
                        is a rotated eccentric anomaly; the equation is the
                        SGP4 modified Kepler equation (Lyddane form).

THEORY
  Underlying theorem:   Newton–Raphson iteration on R(x). For R ∈ C² with
                        R'(x*) ≠ 0, the iteration
                          x_{k+1} = x_k − R(x_k) / R'(x_k)
                        converges quadratically to the simple root x* in a
                        neighborhood (Kantorovich; standard).
                        For SGP4 e ≤ ~0.05 (near-Earth regime) and any e<1
                        more generally, R'(x) = 1 + a_yn·sin(x) + a_xN·cos(x)
                        satisfies |R'(x) − 1| ≤ |a_yn| + |a_xN| ≤ e < 1, so
                        R' > 0 everywhere and the simple-root condition
                        holds globally.
  Primary reference:    Hoots & Roehrich (1980), Spacetrack Report No. 3,
                        §6 (page 13) — modified Kepler equation in form
                        x − a_yn·cos(x) + a_xN·sin(x) = U with Newton
                        iteration. Standard treatment in Battin (1999)
                        §5.3 and Conte & de Boor (1980) §3.5 for
                        Newton on a smooth nonlinear residual.
                        Header references `design/derivations/
                        009_sgp4_modified_kepler.md`.
  Domain of validity:   e ∈ [0, 1) via |(a_xN, a_yn)| < 1 (eccentricity
                        vector norm < 1). Any U ∈ ℝ; periodic in U mod 2π.
                        Starter x₀ := U gives correct ~e-leading order.

METHOD
  Method declared:      Newton–Raphson. Header docstring states
                          f(x) = U − a_yn·cos(x) + a_xN·sin(x) − x,
                          f'(x) = a_yn·sin(x) + a_xN·cos(x) − 1,
                          Δx = f / f',
                        update x ← x + Δx, terminate at |Δx| < tolerance,
                        cap 30 iterations.
                        Note: code uses f = −R and f' = −R', so
                        Δ = f / f' = (−R)/(−R') = R/R', and the update
                        x ← x + Δ equals the standard x ← x − R/R'·(−1)
                        ... let us verify carefully:
                          R = −f  ⇒  Newton:  x ← x − R/R' = x − (−f)/(−fp)
                                                            = x − f/fp.
                        But code uses x ← x + f/fp. Re-check signs:
                          f := U − a_yn·cos(x) + a_xN·sin(x) − x = −R(x).
                          df/dx = a_yn·sin(x) + a_xN·cos(x) − 1.
                          fp (code) = −a_yn·sin(x) − a_xN·cos(x) + 1
                                    = −df/dx.
                        So fp_code = −f'_math. Then
                          Δ_code = f / fp_code = (−R) / (−R') · (−1) = ...
                        Cleaner: x ← x − R/R' is standard Newton. Code's
                        Δ = f/fp = (−R)/(1 + a_yn·sin(x) + a_xN·cos(x))
                                 = (−R) / R' = − R/R'.
                        Then x ← x + Δ = x − R/R'. ✓ matches Newton on R.
  Method implemented:   Lines 47-58: starter x₀ ← U; loop ≤ 30 iterations;
                        evaluate sx, cx, f, fp; Δ ← f/fp;
                        x ← x + Δ; exit on |Δ.value| < tolerance.
                        Each iteration uses one division (the Newton step).
                        Strictly Newton — no Halley term, no continued
                        fraction, no fixed-point unrolling.
  Match verdict:        ✓ matched — single-derivative Newton–Raphson on
                        the SGP4 modified-Kepler residual R(x). NOT a
                        Halley step, NOT a continued-fraction iteration,
                        NOT a series in e.

ERROR BOUND
  Bound category:       precision (REQ-EF-5: iteration residual contributes
                        to precision)
  Bound formula:        For Newton–Raphson with quadratic convergence,
                        once the iterates enter the basin where
                          |x_{k+1} − x*| ≤ C · |x_k − x*|²
                        with C·|x_0 − x*| < 1 (Kantorovich condition),
                        the final correction Δ_k = x_{k+1} − x_k satisfies
                          |x_k − x*| ≤ |Δ_k| / (1 − q)
                        where q := |x_{k+1} − x*| / |x_k − x*| → 0 in
                        the quadratic regime. The simpler conservative
                        bound |x_k − x*| ≈ |Δ_k| is rigorous to leading
                        order once q is small. This is the standard
                        Ostrowski/Kantorovich-style termination bound
                        used for any iterative method (REQ-EF-5).
  Bound implemented:    Line 56:
                          x.errors.precision = x.errors.precision
                                              + abs(delta.value);
                        i.e. the magnitude of the final accepted
                        correction is added to the precision budget.
  Bound verdict:        ✓ matched — the implementation adds |Δ_final|
                        as the precision-bound contribution, which is
                        the standard Kantorovich-style bound for Newton
                        in the quadratic regime.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (iteration residual → precision)
  AUD-EF applies:       AUD-EF-4 (iterative algorithms add residual)
  AUD-MC applies:       n/a (root-finder, not an algebra operation)
  Verification test:    tests/test_orbit/ — Newton convergence vs Halley
                        for a range of (a_xN, a_yn, U); compare reported
                        precision to actual residual after termination.
                        Cross-check x against SR3 reference values from
                        Hoots-Roehrich Appendix.

NOTES
  - Starter x₀ = U is "zeroth-order" in e (correct when a_xN = a_yn = 0,
    where x* = U exactly). For SGP4 near-Earth e ≤ 0.05 this gives ~1.5
    digits at iteration 0; quadratic Newton then doubles digits per step,
    so 6 iterations easily exhaust double precision and 30 iterations
    exhaust any practical multiprecision type.
  - The 30-iteration cap is a safety guard for pathological e → 1 cases.
    Fallback (no return inside loop) at line 60 returns x without adding
    a precision bound — a non-converged result is silently passed
    upstream. ⚠ Open: when the cap is hit, the returned value has no
    iteration-residual contribution to errors.precision. Caller cannot
    distinguish converged vs non-converged from the error budget alone.
    Recommend: in the fallback path, add |Δ_last| to precision OR set a
    convergence-failure flag. For SGP4 e ≤ 0.05 this never trips, but
    the file is templated and could be used at larger e.
  - The starter is *not* the Hoots-Roehrich SR3 starter (which uses a
    short-period correction). Header in `solve_kepler` (the classical
    Kepler counterpart in `src/math/kepler.h`) uses a 2nd-order Taylor
    starter; this modified-Kepler version uses the simpler x₀ = U.
    Acceptable because R'(x) is bounded away from 0 by 1 − e.
```

---

## Card 2 — `solve_kepler_halley`

```
=== FORMULA AUDIT CARD ===
ID:                     modified_kepler::solve_kepler_halley
Location:               src/orbit/modified_kepler.h:72-96
Mathematical statement: Solve R(x) := x − a_yn·cos(x) + a_xN·sin(x) − U = 0
                        for x, same residual as Card 1, using Halley's
                        cubic iteration.

THEORY
  Underlying theorem:   Halley's third-order iteration. For R ∈ C³ with
                        R'(x*) ≠ 0, the iteration
                          x_{k+1} = x_k − 2·R(x_k)·R'(x_k) /
                                          (2·R'(x_k)² − R(x_k)·R''(x_k))
                        converges cubically to x* in a neighborhood:
                          |x_{k+1} − x*| ≤ C · |x_k − x*|³.
                        For the SGP4 modified Kepler residual, R, R', R''
                        are entire (trig polynomials), R' > 0 globally
                        for e < 1 (see Card 1), so the simple-root
                        condition holds and cubic convergence applies
                        globally.
  Primary reference:    Halley (1694), "Methodus Nova Accurata et Facilis
                        Inveniendi Radices Aequationum quarumcumque
                        generaliter, sine praevia Reductione". Modern
                        treatment: Conte & de Boor (1980) §3.7, or any
                        text on higher-order root-finding.
                        Application to Kepler (classical form):
                        Battin (1999) §5.3. The SGP4 modified-Kepler
                        residual is structurally identical (smooth trig
                        polynomial residual with bounded derivative), so
                        the cubic-convergence theorem transfers without
                        modification.
                        Reference: `design/derivations/
                        009_sgp4_modified_kepler.md` (per file header).
  Domain of validity:   e ∈ [0, 1), any U. Same domain as Card 1.

METHOD
  Method declared:      Halley iteration. Header docstring states
                          f''(x) = a_yn·cos(x) − a_xN·sin(x)
                          Δx = 2·f·f' / (2·f'² − f·f''),
                        update x ← x + Δx, tolerance termination, cap 15.
                        Sign check (parallel to Card 1):
                          R'  =  1 + a_yn·sin(x) + a_xN·cos(x) = −fp_code.
                          R'' =  a_yn·cos(x) − a_xN·sin(x) = +fpp_code.
                        Halley step on R:
                          ΔR = −2·R·R' / (2·R'² − R·R'')
                             = −2·(−f)·(−fp_code) / (2·fp_code² − (−f)·fpp_code)
                             = −2·f·fp_code / (2·fp_code² + f·fpp_code).
                        Hmm — code has − f·fpp in denominator, not +.
                        Let me redo with full care. With R = −f, R' = −fp,
                        R'' = fpp (where fp_code, fpp_code as in source):
                          numerator (Halley on R):    −2 R R'
                                                    = −2 (−f)(−fp) = −2 f fp.
                          denominator (Halley on R):  2 R'² − R R''
                                                    = 2 fp² − (−f)(fpp)
                                                    = 2 fp² + f fpp.
                          ΔR = −2 f fp / (2 fp² + f fpp).
                        x_new = x + ΔR = x − 2 f fp / (2 fp² + f fpp).
                        But code computes
                          delta = 2 f fp / (2 fp² − f fpp),
                          x_new = x + delta
                                = x + 2 f fp / (2 fp² − f fpp).
                        These differ in the sign of the f·fpp term in the
                        denominator AND in the overall sign of the
                        numerator. Re-check by working directly on f
                        (since the iteration is invariant under sign of R):
                        Halley on f (find root of f):
                          Δ_f = −2 f f' / (2 f'² − f f''),
                          f' = df/dx = a_yn·sin(x) + a_xN·cos(x) − 1
                                     = −fp_code,
                          f'' = d²f/dx² = a_yn·cos(x) − a_xN·sin(x)
                                        = fpp_code.
                          Δ_f = −2 f (−fp_code) /
                                 (2 (−fp_code)² − f·fpp_code)
                              = 2 f fp_code / (2 fp_code² − f fpp_code).
                          x_new = x + Δ_f (since Halley on f gives the
                                            same root as on R = −f).
                        This exactly matches the code. ✓
                        So the code applies Halley to f := −R with the
                        signs of f', f'' correctly carried through. The
                        result is identical to Halley on R; the sign
                        gymnastics are bookkeeping.
  Method implemented:   Lines 79-94: starter x₀ ← U; loop ≤ 15 iterations;
                        evaluate sx, cx, f, fp, fpp; build
                          num = 2·f·fp, den = 2·fp² − f·fpp;
                        Δ ← num/den; x ← x + Δ; exit on |Δ| < tol.
                        Three sin/cos uses per iteration (one each + reuses);
                        explicit second-derivative term — cubic-correction
                        signature.
  Match verdict:        ✓ matched — Halley's method (uses f, f', f''),
                        cubic convergence by construction. NOT Newton (no
                        Halley term would be present), NOT a continued
                        fraction, NOT a series expansion. The sign
                        bookkeeping (working in f = −R rather than R)
                        produces the standard cubic correction with
                        denominator (2·f'² − f·f'') and numerator (2·f·f'),
                        matching e.g. Conte & de Boor (1980) §3.7
                        Eq. (3.7.5).

ERROR BOUND
  Bound category:       precision (REQ-EF-5)
  Bound formula:        For Halley iteration in the cubic regime,
                          |x_{k+1} − x*| ≤ C · |x_k − x*|³,
                        and the final-correction bound
                          |x_k − x*| ≤ |Δ_k| / (1 − q²)
                        with q := |x_{k+1} − x*| / |x_k − x*| → 0
                        cubically. Conservative leading-order bound:
                          |x_k − x*| ≈ |Δ_k|.
                        This is the same Kantorovich-style bound as
                        Newton, with the additional safety that q
                        contracts faster (cubic vs quadratic). Used by
                        REQ-EF-5 for any iterative method.
  Bound implemented:    Line 91:
                          x.errors.precision = x.errors.precision
                                              + abs(delta.value);
                        Final-correction magnitude added to precision,
                        identical to Card 1.
  Bound verdict:        ✓ matched — Kantorovich/final-correction bound
                        is the standard precision contribution for any
                        iterative root-finder and is tighter for Halley
                        than Newton at the same iteration count.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (iteration residual → precision)
  AUD-EF applies:       AUD-EF-4 (iterative algorithms add residual)
  AUD-MC applies:       n/a
  Verification test:    tests/test_orbit/ — Halley should converge in
                        ~4 iterations for double precision, ~5-6 for
                        cpp_bin_float_50, well under the 15-iteration
                        cap. Compare digits-per-iter against Newton (~6
                        iterations for double): Halley should win at
                        ~1.5× speed per iteration count, ~equal cost
                        per CPU time due to the extra sin/cos and the
                        second-derivative arithmetic.

NOTES
  - The 15-iteration cap is half Newton's 30 (rationale: 3^15 ≈ 1.4×10⁷
    digits of headroom, more than any representable type holds).
  - Same fallback caveat as Card 1: if the loop exits without converging,
    line 95 returns x with no residual added to errors.precision. ⚠ Open:
    cap-hit path silently returns under-bounded x.
  - Header line 19 advertises "Householder" as a third solver; this is
    NOT implemented in this file. Only Newton and Halley are present.
    The Householder method (quartic convergence using f, f', f'', f''')
    would be a natural extension but is out of scope for this audit.
  - Halley's method is the m = 2 case of Householder's family, where
    m = 1 is Newton. The header's "Householder quartic" would be m = 3.
  - **No theory-method mismatch**: the code implements the standard
    cubic Halley correction Δ = 2·f·f' / (2·f'² − f·f'') verbatim,
    matching the theory citation.
```

---

## File-level verdict

| Dimension | Verdict | Notes |
|---|---|---|
| **A. Error wiring** | ✓ | Both functions add `|Δ_final|` to `errors.precision` per REQ-EF-5 / AUD-EF-4. |
| **B. Algebra axioms** | n/a | Numerical root-finders, not algebra operations. |
| **C. Theoretical basis (Card 1)** | ✓ PASS | Newton on R(x); single-derivative correction; bound = Kantorovich \|Δ\|. |
| **C. Theoretical basis (Card 2)** | ✓ PASS | Halley on R(x); uses f, f', f''; cubic-convergence bound = Kantorovich \|Δ\|. |

**File verdict: PASS.**

## Open items

1. **⚠ Cap-hit fallback** (both functions): when the iteration cap is
   exhausted without `|Δ| < tolerance`, the function returns `x` with
   *no* contribution to `errors.precision`. A non-converged result is
   indistinguishable from a converged one in the error budget. Either
   add `|Δ_last|` in the fallback path or set a convergence-failure
   sentinel. For SGP4 near-Earth (e ≤ 0.05) this never trips, but the
   templated code admits arbitrary eccentricity.
2. **Header advertises Householder** (line 20) but the file implements
   only Newton and Halley. Either implement Householder or amend the
   header to drop the claim. Out of scope for this audit card.
3. **Starter x₀ = U** is "zeroth-order in e" — accurate to ~e for
   small eccentricity. For SGP4 near-Earth this is fine; for larger e
   a Taylor-in-e starter (e.g. `x₀ = U + a_yn·cos(U) − a_xN·sin(U)`
   for first-order) would reduce iteration count. Not a correctness
   issue; a performance tuning.

## Cross-references

- Companion file `src/math/kepler.h` (classical Kepler equation
  `E − e·sin(E) = M`) — audit card stubbed in
  `theoretical_basis_audit.md` §6.
- Framework: `design/audit/theoretical_basis_audit.md` §1 (template),
  §5 (worked example), §6 (kepler.h stub).
- Derivation reference: `design/derivations/009_sgp4_modified_kepler.md`
  (per source header).
- Primary source: Hoots & Roehrich (1980), SR3 §6 page 13.
