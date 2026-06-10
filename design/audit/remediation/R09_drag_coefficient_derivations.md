# R09 — Drag coefficient literal derivations

**Status:** DEFERRED — agent's 625-line derivation extension did not persist (file-revert pattern during parallel dispatch). Needs re-execution in follow-on session.
**Severity:** P3 (theory-rigor gap; not correctness-blocking)
**Estimated scope:** ~4-6 hours of symbolic derivation work

---

## Files

**Write:**
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` — extend Part III §11 (D₃/D₄), §12 (t-cofs); add new section formalizing Lane f† substitution

**Read-and-annotate-cite-only (no formula changes):**
- `src/atmosphere/drag_coefficients.h`

**Read:**
- `sgp4_references/vallado_celestrak/documentation/SGP4/Lane_Hoots_1979_General_Perturbations_Lane_Drag.pdf` [LH79] — primary source for D₃/D₄ literals
- `sgp4_references/Hoots_Roehrich_1980_Spacetrack_Report_No3.pdf` §6

**Audit card:**
- `design/audit/theoretical_basis_audit/drag_coefficients.md`

---

## Primary issues (theory-rigor gaps)

Five gaps flagged in `AUD_TBA_results.md` Part VII and the audit card for drag_coefficients.h:

### 1. D₃ literal coefficients `17·a₀ + s`

Currently transcribed from [LH79 p. 26], not symbolically re-derived. Need to show step-by-step that the Taylor expansion of the orbit-averaged drag integrand at order 3 gives exactly `(17a₀ + s)/3` for the `a²ξ²C₁³` coefficient.

### 2. D₄ literal coefficients `221·a₀ + 31·s`

Same — Taylor at order 4 should give exactly `(221a₀ + 31s)/2`. The literals 221 and 31 come from a 4th-derivative chain in `(a−s)⁻ᵏ`; need explicit binomial accounting.

### 3. C₅ literal `11/4`

Sketched in §9 as `11/4 = 3/2 + 5/4`. The `3/2` comes from the η² coefficient in C₂ Part A (already derived). The `5/4` is claimed to come from the `e₀·η` cross-coupling in the orbit average of `∂(dM)/∂t`. Needs rigorous evaluation of the orbit average:

```
⟨∂(dM)/∂t · cos M⟩_M = ?
```

with the next-order eccentricity expansion.

### 4. Lane's f† harmonic substitution

§4 Note 4.3 currently says the substitution `r - s = (1/ξ)(1 - η·cos f†)` is exact only at first order in e, with the dropped corrections O(eη). Should be:

- Define f† explicitly as a series in f, e
- Bound the truncation error in the orbit-averaged density integral
- Express the bound as a contribution to `errors.accuracy` (will feed into R06)

### 5. t₃cof, t₄cof, t₅cof term-by-term assembly

Part III §12 sketches the derivation via binomial expansion of `(1-x)^{-3/2}` for `g(τ)^{-3/2}`. The full term-by-term matching at τ⁴ and τ⁵ is straightforward but tedious; should be written out as a numbered table.

---

## Theory anchor

Lane (1965), Lane-Hoots (1979) [LH79], Hoots-Roehrich 1980 SR3 §6.

The derivation document `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` is the canonical artefact. This R-item extends it.

---

## Fix scope

~4-6 hours of focused symbolic-derivation work. No source code changes (the numerical values are already verified against 7 independent reference SGP4 implementations and the original 1980 FORTRAN; this is only about closing theoretical-rigor gaps in the derivation document).

Approach:

1. Open `design/derivations/sgp4_near_earth_drag_theoretical_basis.md`.
2. Extend §11 with full Taylor-derivative chain for D₃, D₄. Show literal 17, 221, 31 as outputs of explicit binomial coefficients.
3. Extend §9 with the `5/4` orbit-average evaluation for C₅.
4. Add new §4a (or extend §4) formalizing the f† substitution with explicit truncation bound.
5. Extend §12 with the term-by-term table for t-cofs.
6. Cross-reference the corresponding `drag_coefficients.h` lines in each new theorem.

---

## Verification

No new code; verification is documentary:

- Each previously-transcribed literal is now traceable to a numbered theorem in the derivation document.
- Optional: SymPy or Mathematica verifier in `design/derivations/` that mechanically checks the closed forms (matches existing style in BH61 work).
- Audit card `drag_coefficients.md` updated to reflect that the rigor gaps are closed.

---

## References

- Audit card: `design/audit/theoretical_basis_audit/drag_coefficients.md`
- Existing derivation: `design/derivations/sgp4_near_earth_drag_theoretical_basis.md`
- Primary sources: Lane-Hoots (1979) p. 16 (orbit-average integrals), p. 25-26 (D₂..D₄ formulas)

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
