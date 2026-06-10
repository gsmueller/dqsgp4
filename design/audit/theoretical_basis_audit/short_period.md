# Theoretical Basis Audit — `src/perturbation/short_period.h`

**File**: `src/perturbation/short_period.h` (119 lines)
**Functions audited**: 1 (`apply_short_period`)
**Expected method per §7**: Closed-form algebraic application of SR3 §6 short-period J₂ corrections.

---

## Card 1 — `apply_short_period`

```
=== FORMULA AUDIT CARD ===
ID:                     short_period::apply_short_period
Location:               src/perturbation/short_period.h:69-117
Mathematical statement: Six closed-form first-order J₂ short-period
                        corrections {rk, uk, xinck, xnodek, rdotk, rfdotk}
                        from the SR3 §6 expressions:
                          rk      = r·(1 − (3/2) sp3 β (3cos²i−1)) + (1/2) sp2 sin²i cos2u
                          uk      = u − (1/4) sp3 (7cos²i−1) sin2u
                          xinck   = i₀ + (3/2) sp3 cosi sini cos2u
                          xnodek  = Ω + (3/2) sp3 cosi sin2u
                          rdotk   = ṙ − n sp2 sin²i sin2u
                          rfdotk  = rḟ + n sp2 (sin²i cos2u + (3/2)(3cos²i−1))
                        with sp2 = (J₂/2)/p, sp3 = (J₂/2)/p².

THEORY
  Underlying theorem:   Brouwer (1959) first-order von-Zeipel short-period
                        transformation: the periodic part of the J₂
                        generating function S₁(L,G,H,l,g,h) produces closed-
                        form (l-periodic, equivalently u = f+g periodic)
                        corrections to the osculating elements via
                          δX = {X, S₁}  (Poisson bracket).
                        After projection to (r, u, i, Ω, ṙ, rḟ) the
                        corrections are *closed-form algebraic identities*
                        in {sin2u, cos2u, sin²i, cos²i, β=√(1−e²), p, n} —
                        no truncation, no iteration, no series.
  Primary reference:    Hoots & Roehrich (1980), Spacetrack Report No. 3,
                        p. 14 — equations for (rk, uk, xinck, xnodek,
                        rdotk, rfdotk). Brouwer (1959) §3 (S₁ definition).
                        Doc comment cites both correctly.
  Domain of validity:   e < 1 (so β=√(1−e²) real and p>0); i₀ avoidance
                        of the critical inclination is not required for
                        these *short-period* terms (the critical-inclination
                        singularity sits in the long-period generator
                        S₁* via the factor 1/(5cos²i−1), which does NOT
                        appear here). All trigonometric inputs assumed
                        in standard ranges. Function applies to both
                        near-space and deep-space paths per doc.

METHOD
  Method declared:      Closed-form algebraic substitution of SR3 §6
                        formulas — no truncation, no series, no iteration.
  Method implemented:   Lines 95-114: direct arithmetic on TrackedValue
                        inputs using `exact<T>` (1) and `ratio<T>` (3/2,
                        1/2, 1/4) prefactors. Each of the six output
                        components is one straight algebraic expression
                        in the inputs. No loops, no convergence test,
                        no truncated series.
  Match verdict:        ✓ matched — closed-form identity, exactly as SR3
                        §6 (and Brouwer 1959) prescribe.

ERROR BOUND
  Bound category:       precision (per-op TrackedValue propagation)
                        + accuracy (model-level via doc note)
  Bound formula:        Per REQ-EF-3, every TrackedValue arithmetic op
                        (* + − /) composes the operand precision bounds
                        via the closed-form rules in tracked_value.h.
                        For this function specifically there is NO
                        method-level truncation bound: the corrections
                        are *exact* algebraic transcriptions of S₁'s
                        first-order short-period contribution. Therefore
                        the precision bound is the cumulative composition
                        of input bounds through the ~30 arithmetic ops
                        per output component.
  Bound implemented:    Implicitly: every `*`, `+`, `−`, `/` on
                        TrackedValue<T> triggers the closed-form bound
                        composition defined in tracked_value.h. No
                        explicit `errors.precision += …` adjustment is
                        made by this function (and none is required —
                        closed-form ops carry their own bounds).
  Bound verdict:        ✓ matched — closed-form composition is the
                        rigorous bound for a closed-form identity
                        (REQ-EF-3). No additional truncation term is
                        owed because there is no truncation.

                        ⚠ One *accuracy* item, properly handled: the
                        doc comment (lines 63-68) flags the omitted
                        higher-order J₂² and J₄ short-period terms
                        (~O(J₂² · r) ≈ 7 mm) as below the secular-rate
                        truncation accuracy floor (~100 m), so they
                        are deliberately NOT added to `errors.accuracy`.
                        This is a *model accuracy* choice, not a
                        formula-method mismatch. Per §1 the audit card
                        records this as a note, not a fail.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation: every
                        TrackedValue arithmetic op composes precision
                        bounds in closed form).
  AUD-EF applies:       AUD-EF-3 (closed-form ops carry input bounds
                        without explicit additions).
  AUD-MC applies:       n/a (this function is a physical-model
                        application, not an algebra operation).
  Verification test:    tests/ — short-period corrections should be
                        exercised against a Vallado / SR3 reference
                        propagation at chosen (e, i, u). Cross-check
                        of `rk, uk, xinck, xnodek, rdotk, rfdotk`
                        formulas against SR3 page 14 verbatim.

NOTES
  - **No method-theory mismatch.** Closed-form identity cited, closed-form
    arithmetic implemented, closed-form bound composition via
    TrackedValue.
  - **Critical-inclination clean.** No 1/(5cos²i−1) factor anywhere —
    that singularity is confined to S₁* / long-period and does not
    contaminate this short-period path.
  - **Conventions match SR3 page 14:**
      sp2 = CK2/p  (with CK2 = J₂/2 in SR3 notation; `half_J2` here)
      sp3 = CK2/p² (= sp2/p)
    Both built once per call (lines 96-97), reused across all six
    outputs — REQ-EF-style reuse, no recomputation.
  - **Input dependencies are pre-derived.** The auxiliaries
    `three_cos2i_minus_1`, `seven_cos2i_minus_1`, `cos2_i0`,
    `sin2_i0`, `sin_2u`, `cos_2u`, `beta_l` are passed in as
    TrackedValues already carrying their precision bounds. This
    means the entire correction pipeline (inputs → corrections)
    is a single chain of closed-form ops, fully audited under
    REQ-EF-3.
  - **i₀ assumption.** Doc states `i0` is held at-epoch for the
    short-period (i.e., short-period correction uses the unperturbed
    inclination, not the secularly-advanced one). This matches SR3
    §6, where the short-period correction to inclination is a
    *correction* to i₀, not to a secular i(t).
  - **Sign conventions.** rk, uk, xinck, xnodek match SR3 signs
    on first inspection (rk has −sp3 β (3cos²i−1) leading;
    uk has −(1/4) sp3 (7cos²i−1) sin2u; xinck has +(3/2) sp3
    cosi sini cos2u; xnodek has +(3/2) sp3 cosi sin2u). Velocity
    corrections rdotk, rfdotk match: rdot loses an n sp2 sin²i
    sin2u term; rfdot gains an n sp2 (sin²i cos2u + (3/2)(3cos²i−1))
    term. ✓ consistent with SR3 page 14 verbatim sign pattern.
  - **Accuracy floor.** Doc comment correctly identifies the J₂² /
    J₄ omission as ~7 mm vs the ~100 m secular-rate-truncation floor;
    no `errors.accuracy +=` is owed for the omitted terms because
    they are below the model's stated accuracy floor.
```

---

## File-level verdict for `short_period.h`

- **A. Error wiring**: ✓ Output components are built entirely from
  TrackedValue arithmetic; precision bounds propagate via REQ-EF-3
  with no explicit per-formula additions required (correct for
  closed-form identities).
- **B. Algebra axioms**: n/a — not an algebra operation; physical-model
  application.
- **C. Theoretical basis**:
  - 1× formula (`apply_short_period`): ✓ closed-form SR3 §6 identity
    cited, closed-form arithmetic implemented, closed-form bound
    composition. **PASS.**

**File verdict: PASS** — one card, all slots matched. Accuracy note on
omitted J₂² / J₄ short-period terms (~7 mm) is correctly documented as
below the model accuracy floor and not added to `errors.accuracy`.
