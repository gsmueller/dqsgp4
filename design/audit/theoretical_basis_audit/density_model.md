# Theoretical Basis Audit — `src/atmosphere/density_model.h`

Audit framework: `design/audit/theoretical_basis_audit.md` §§1, 5.
Source: `src/atmosphere/density_model.h` (99 lines, 1 function).

This file contains a single distinct mathematical operation: `compute_density_parameters`. It is a piecewise closed-form construction of the Lane power-law atmosphere fitting parameters $(s, (q_0-s)^4)$ as functions of perigee altitude. There is no truncation, no iteration, no transcendental series — every step is a finite arithmetic composition over `TrackedValue<T>` plus three `sqrt` calls in the adjusted branches.

---

## Card 1 — `compute_density_parameters`

```
=== FORMULA AUDIT CARD ===
ID:                     density_model::compute_density_parameters
Location:               src/atmosphere/density_model.h:56-96
Mathematical statement: Given perigee altitude h_p (km), semi-major axis a₀ (ER),
                        eccentricity e₀, and Earth radius R_E (km), construct
                        the Lane atmosphere fitting parameters

                          s ∈ ER,  (q₀ − s)⁴ ∈ ER⁴

                        via the three-regime piecewise rule (h_p in km):

                          Case A  h_p ≥ 156:
                            s        = 1 + 78/R_E
                            (q₀−s)⁴  = (42/R_E)⁴

                          Case B  98 ≤ h_p < 156:
                            s*       = a₀(1 − e₀) − s_nom + 1
                            s        = 1 + s*
                            (q₀−s)⁴  = [ (q₀−s_nom)        ← nominal 42/R_E
                                        + s_nom − (1 + s*) ]⁴

                          Case C  h_p < 98:
                            s*       = 20/R_E
                            s        = 1 + s*
                            (q₀−s)⁴  = [ (q₀−s_nom) + s_nom − (1 + s*) ]⁴

                        with constants q₀ = 1 + 120/R_E, s_nom = 1 + 78/R_E.
                        (Code reads (q₀−s_nom) as `qoms_root4 := √√((q₀−s_nom)⁴)`
                        rather than re-deriving 42/R_E; algebraically identical
                        for non-negative arguments.)

THEORY
  Underlying theorem:   None — this is a *definition* (Lane 1965 / Hoots-Roehrich
                        1980 §6, the perigee-dependent piecewise fitting recipe
                        for the SGP4 power-law atmosphere). It is **not** the
                        application of a convergence theorem (Taylor, Newton,
                        Leibniz, Halley, …); it is the construction of the
                        density model's parameters.

                        Supporting definitions:
                          • Definition 2.1 of `sgp4_near_earth_drag_theoretical_basis.md`
                            (Lane power-law form ρ(r) = ρ₀ ((q₀−s)/(r−s))⁴).
                          • Definition 2.3 (perigee-adjusted parameters,
                            three cases at 156 km and 98 km).

                        Underlying numerical theorem (for the bound on each
                        intermediate `TrackedValue<T>` arithmetic step):
                        triangle / mean-value bounds for `+`, `−`, `*`, `/`
                        as wired in REQ-EF-3, and the Lipschitz bound
                        |√(x+δ) − √x| ≤ δ/(2√x) for the two `sqrt` calls in
                        the adjusted branches (Cases B and C).

  Primary reference:    1. Lane, M.H. (1965), "The development of an artificial
                           satellite theory using a power-law atmospheric
                           density representation", AIAA paper 65-35.
                        2. Hoots, F.R. and Roehrich, R.L. (1980), "Models for
                           Propagation of NORAD Element Sets" (Spacetrack Report
                           No. 3), §6 (page 10 — perigee adjustment block).
                        3. Lane, M.H. and Hoots, F.R. (1979), "General Perturbations
                           Theories Derived from the 1965 Lane Drag Theory",
                           the s₀=78 km, q₀=120 km fitting constants.
                        4. `design/derivations/sgp4_near_earth_drag_theoretical_basis.md`
                           §2.3 — Definition 2.3 "Perigee-adjusted density parameters".

  Domain of validity:   • All cases require a₀ ≥ 1 (orbit above Earth surface) and
                          e₀ ∈ [0, 1) (bound orbit) for the perigee to be a real,
                          non-degenerate altitude.
                        • Case A: nominal physical regime (perigee in the
                          150-500 km band the τ=4 fit was tuned for).
                        • Case B: degraded-accuracy regime; the adjustment is a
                          documented engineering workaround, not a theorem (the
                          power-law form was never fit at h_p < 156 km).
                        • Case C: h_p < 98 km; orbits at this altitude are
                          decaying within hours; SGP4 is being used out of its
                          intended regime.
                        • All branches require the inner factor
                          (q₀−s_nom) + s_nom − (1+s*) ≥ 0 for the fourth power
                          to be meaningful (the model assumes the perigee lies
                          above the chosen s* surface, so the inner factor is
                          (q₀ − (1+s*)) = (120/R_E − s*) which is positive for
                          s* < 120/R_E ⇔ adjusted altitude < 120 km — satisfied
                          for all three cases by construction).

METHOD
  Method declared:      Closed-form piecewise evaluation of Definition 2.3,
                        following the Spacetrack Report No. 3 §6 recipe verbatim.
                        Three regimes selected by raw-comparison of
                        `perigee_km.value` to T-constants 156 and 98. No series,
                        no iteration.

  Method implemented:   src/atmosphere/density_model.h:62-95.
                        Step-by-step:
                          L67  params.s_nom  = 1 + 78/R_E          (+,/ on TV<T>)
                          L70  q₀            = 1 + 120/R_E         (+,/ on TV<T>)
                          L73  diff          = q₀ − s_nom          (− on TV<T>)
                          L74  qoms4_nom     = diff⁴              (3 × * on TV<T>)
                          L77  branch on perigee_km.value           (raw T cmp)
                          L79  s*            = a₀(1−e₀) − s_nom + 1 (∗,−,−,+)
                          L81  qoms_root4    = √√qoms4_nom          (2 × sqrt)
                          L82  new_qoms      = qoms_root4 + s_nom − (1+s*)
                                                                    (+,−,+)
                          L83  params.qoms4  = new_qoms⁴            (3 × *)
                          L84  params.s      = 1 + s*               (+)
                        Case C is the same as Case B with s* = 20/R_E.

  Match verdict:        ✓ matched. The implementation is the literal evaluation
                        of Definition 2.3 in the cited derivation document. No
                        approximation is introduced beyond the framework choice
                        of the Lane model itself (which is an empirical fit, not
                        a derivation — see U-D2 in the §13 error budget of the
                        drag basis document).

                        Note on Case B's algebraic structure: the recipe
                        "take 4th-root, add (s_nom − s_new), raise to 4th"
                        is mathematically equivalent to evaluating
                        (q₀ − s_new)⁴ directly with q₀ = 1 + 120/R_E. The code
                        uses the √√/⁴ route, which matches the SR3 page-10
                        printed recipe. Both forms agree exactly in real
                        arithmetic and within `TrackedValue` propagation;
                        the √√/⁴ route accumulates two extra Lipschitz bounds
                        (one per sqrt) — see ERROR BOUND below.

ERROR BOUND
  Bound category:       Multiple, by per-operation routing through
                        `TrackedValue<T>` arithmetic:

                          • measurement     ← propagates from a₀, e₀, perigee_km,
                                              re_km input categories (these are
                                              recovered from TLE measurements
                                              by `element_recovery`).
                          • precision       ← rounding of each TV<T> op
                                              (REQ-EF-3 per-category bounds);
                                              plus the two sqrt Lipschitz
                                              residuals on Cases B, C.
                          • accuracy        ← the *model* accuracy carried in
                                              the inputs; this routine does
                                              not add a model-error term of
                                              its own (the Lane-model 30%
                                              fit error documented in
                                              `sgp4_near_earth_drag_theoretical_basis.md`
                                              §2 is a property of the *density
                                              evaluator*, not of the parameter
                                              construction).

  Bound formula:        Each elementary step (`+`, `−`, `*`, `/`, `sqrt`)
                        contributes its REQ-EF-3 closed-form bound per
                        category. Specifically:

                          δ(x + y)  =  δx + δy            per category
                          δ(x · y)  ≤  |y|·δx + |x|·δy + δx·δy  per category
                                                          (mul_bound in tracked_value.h:227)
                          δ(x / y)  ≤  (|y|·δx + |x|·δy)/y² + …  (div_bound, l.266-268)
                          δ(√x)     ≤  δx / (2√x)         per category (Lipschitz)

                        For Case A the inner factor (q₀−s_nom) involves only
                        +, −, /, *; no sqrt. The fourth power is three
                        multiplications. Bound is the standard product
                        propagation of three `mul_bound` calls.

                        For Cases B and C the bound additionally includes
                        two sqrt-Lipschitz residuals on `qoms_root4`. Since
                        the argument is (42/R_E)⁴ ≈ 4.4×10⁻¹⁰ ER⁴, the
                        inner √ produces 42²/R_E² and the outer √ produces
                        42/R_E ≈ 6.6×10⁻⁶. Both sqrt evaluations are
                        deep into the positive interior, so the Lipschitz
                        denominator 2√x is well-bounded away from zero, and
                        the residual is the standard δx / (2√x) per category.

                        No truncation bound is added at the
                        `compute_density_parameters` level (REQ-EF-6 does
                        not apply — there is no Taylor truncation here);
                        no iteration residual is added (REQ-EF-5 does not
                        apply — there is no Newton/Halley loop).

                        The *modeling* error (the 30% Lane-vs-true-atmosphere
                        gap; the C⁰-only continuity across h_p = 156 km
                        producing ~1% drag-rate jumps) is a separate accuracy
                        term that lives on the density *evaluator* ρ(r),
                        not on the parameters themselves. The parameter
                        construction is exact in the model; the model is
                        wrong by 30%, but that is a downstream consumer's
                        concern.

  Bound implemented:    The code does NOT explicitly add a bound at
                        `compute_density_parameters` scope. All bounds
                        propagate implicitly via the `TrackedValue<T>`
                        operator overloads called on lines 67-91:
                          • +, −, *, /  call the REQ-EF-3-conformant
                            overloads in src/math/tracked_value.h.
                          • sqrt calls the sqrt overload (Lipschitz bound).
                          • Branch selection on `perigee_km.value` uses the
                            raw T component, NOT the tracked uncertainty
                            (see NOTES below — flagged caveat).

  Bound verdict:        ✓ matched, with one ⚠ note on branch-selection
                        sensitivity. The construction is fully decomposed
                        into `TrackedValue` ops, each of which is audited
                        elsewhere (REQ-EF-3 and the per-op cards in
                        tracked_value.md when those land). No bounds are
                        dropped at this level. The branch boundary is the
                        only non-smooth feature, and that is a
                        modeling/discontinuity issue covered by the
                        cited derivation document's error budget rather
                        than by this card's method-bound contract.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form per-category propagation, applies
                                  to every +/−/*//* and sqrt in the body).
                        REQ-EF-5/6  N/A (no iteration, no Taylor truncation).
  AUD-EF applies:       AUD-EF-3 (per-category bounds wired through
                                  composite operations). The composite is
                                  exercised whenever `compute_density_parameters`
                                  is called inside `drag.h`/`drag_coefficients.h`.
  AUD-MC applies:       N/A — this is a parameter-construction utility, not
                        an algebra operation; no algebraic identity is
                        being implemented.
  Verification test:    Expected: `tests/test_atmosphere/` exercising
                        compute_density_parameters against Hoots-Roehrich
                        1980 SR3 page-10 worked numbers across all three
                        regimes (perigee at 200 km, 120 km, 50 km). Also
                        a continuity test at h_p = 156 km⁻ vs 156 km⁺
                        (C⁰ continuity should hold to machine ε; C¹ should
                        NOT — the documented kink).
                        Status: ?  — test presence not confirmed by this
                        audit; flag for follow-up.

NOTES
  - **Branch selection uses bare T comparisons** on lines 78 and 86:
    `perigee_km.value < T(156)` and `perigee_km.value < T(98)`. The
    measurement uncertainty on `perigee_km` (carried in `perigee_km.errors`)
    is **NOT** consulted. This is the correct SR3-faithful behavior — the
    perigee threshold check is part of the model definition — but it means
    a TLE whose perigee sits within its measurement uncertainty of the
    156 km boundary will deterministically pick one branch even though
    "the other branch" is statistically possible. The C⁰-only continuity
    cited in Definition 2.3 of the drag basis document means the choice
    of branch matters: a 1% drag-rate jump across the boundary becomes a
    deterministic-vs-actual selection sensitivity. **Flag for downstream
    review** (whether the propagator should re-evaluate near boundaries,
    or carry a "near-boundary" accuracy hit).

  - **The §2.3 30% atmosphere fit error is a Definition-2.1 property**,
    not a Definition-2.3 property. It belongs on the density *evaluator*
    ρ(r) wherever that is implemented, not on `compute_density_parameters`.
    Confirm `drag.h`'s ρ(r) call adds this accuracy term — the search is
    a downstream audit item.

  - **The √√/⁴ identity is benign for double**. Case A bypasses both
    sqrts; Cases B and C take √√((42/R_E)⁴) = 42/R_E with the only
    extra cost being two Lipschitz residuals, each ~ε|42/R_E|/2 ≈
    7×10⁻²² ER — 6 orders of magnitude below the input measurement
    precision on R_E itself. The detour costs nothing in practice but
    is **un-needed**: the equivalent expression
    `(q₀ − (1 + s_star))⁴` evaluates the same thing without the
    √√/⁴ round-trip and saves the two Lipschitz residuals. **Flag for
    cleanup** (a future tightening; correctness is unaffected).

  - **No theory-method mismatch**. This is a definition, not a theorem
    application. The §1 framework's "method declared vs method implemented"
    triad collapses to "Definition 2.3 declared, Definition 2.3
    implemented." The Match Verdict accordingly checks only the
    typographic-equivalence of the code to the printed recipe — which
    holds.

  - **Lane drag is a model used by the SGP4 path, not by the dual-quaternion
    propagator's drag** (see `MEMORY.md` `project_validation_status.md` —
    "SGP4 uses Lane drag theory not BH61"). The audit framework's status
    table in §3 of the framework document lists this file in the DQ-propagator
    scope; this should be cross-checked. If `compute_density_parameters` is
    only consumed by the SGP4 path, it is out of scope for the DQ-propagator
    audit, and the §3 table entry should be removed from that table or
    routed to a parallel "SGP4-path TBA" sweep. **Scope-flag for the
    parent audit framework owner**.
```

---

## File-level verdict for `density_model.h`

- **A. Error wiring**: ✓ propagation happens implicitly via the `TrackedValue<T>` operator overloads called on every line of the body. No bare-`T` slipthrough. (Verifier: AUD-EF-3.)
- **B. Algebra axioms**: n/a — this is a parameter-construction utility, not an algebra operation.
- **C. Theoretical basis**:
  - Card 1 `compute_density_parameters`:
    - Theory ↔ Method: ✓ (Definition 2.3 declared, Definition 2.3 implemented; this is a definition not a theorem so no Taylor-vs-Padé class of mismatch is available).
    - Method ↔ Bound: ✓ for the per-operation level (REQ-EF-3 propagation through every `+`, `−`, `*`, `/`, `sqrt`); ⚠ for the model-discontinuity term at h_p = 156 km, which is a known C⁰-only feature documented in the derivation document but not surfaced by this routine into the `accuracy` category.
  - Three open flags:
    - ⚠ branch-selection ignores the measurement uncertainty on `perigee_km` (correct per SR3; flagged for downstream review of near-boundary handling).
    - ⚠ the √√/⁴ detour in Cases B and C is algebraically equivalent to the direct `(q₀−(1+s*))⁴` form and adds two Lipschitz residuals; a tightening candidate.
    - ? test coverage at `tests/test_atmosphere/` not confirmed by this audit.
    - ? scope: this file is on the SGP4 drag path; whether it belongs in the DQ-propagator §3 audit table is a parent-framework question.

**File verdict: PASS-with-flags** — the formula matches the cited definition, the bound is the per-operation `TrackedValue<T>` propagation it should be, and the modeling-error term (Lane 30%, C⁰ kink) is correctly externalized to the density evaluator rather than the parameter constructor. The four flags are clarifications, not C-fails.
