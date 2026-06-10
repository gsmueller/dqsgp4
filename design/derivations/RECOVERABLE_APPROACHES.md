# Recoverable Aspects of Abandoned / Superseded Approaches

> Created 2026-06-03 (maintenance/accounting pass). The project has iterated through several
> superseded approaches (early `deprecated/` derivations, alternative `drafts/`, the retired 38-chapter
> "Mathematical Foundations" outline, the legacy unified drag doc). The *current* approaches are
> generally more correct and more complete — but a few abandoned approaches had a genuinely **superior
> aspect** that was lost in the transition and is worth recovering. This note records only those, vetted
> for real value (an Explore agent surveyed structure/pedagogy — no derivation equations transcribed —
> and these were filtered down by hand). "Recover" = cheap doc/index work, **not** rewriting the math.

## Worth recovering (high value)

### 1. The legacy drag doc's embedded audit trail (D-1 … D-11)
- **What was lost:** `sgp4_near_earth_drag_theoretical_basis.md` (now "legacy/synopsis") carried the
  error-source annotations **D-1…D-11 inline with the derivation** — *where* an approximation/error was,
  and how it was corrected. The per-coefficient trace docs (`sgp4_drag_*_trace.md`) have the correct
  math but **scattered the audit trail**; there is no single index of "which D-finding was resolved where."
- **Recover:** a one-page **D-finding resolution index** (D-1…D-11 → the trace doc/section + verifier that
  closed it). Pure cross-reference; high pedagogical + audit value; complements `AUDIT_BACKLOG.md`.

### 2. A written multi-agent / cleanroom verification SOP
- **What was lost:** the `phase_b_*` workbooks + the 3-agent symbolic cleanroom (solver A/B/C pairwise
  `simplify(diff)=0`) embody a *reusable process* for catching derivation errors at scale (it caught a
  real solver error and validated results 81/81). The process lives in memory + `BIAS_PREVENTION_PROTOCOL.md`
  but there is no concise, in-repo SOP a newcomer can follow.
- **Recover:** a short **Derivation-Verification SOP** (outline → N independent solvers → pairwise
  symbolic cross-validation → synthesis → commit; with the caught-error case as the worked example).
  High value for continuity; the discipline already exists, just isn't packaged.

## Worth a lightweight note (medium value)

### 3. The 38-chapter outline's prerequisite structure
- **What was lost:** `00_table_of_contents.md`'s 9-part outline made chapter *dependencies* visible
  (e.g. ch14 equipotential-ellipsoid → ch15 Kaula). The current organization (scattered trace docs +
  a separate submodule `MASTER_INDEX.md`) is not reader-discoverable as a dependency graph.
- **Recover (optional):** a dependency map of **what actually exists** (ch01-05, ch14, the trace docs)
  overlaid on the original part structure — without committing to writing the ~30 missing chapters.

### 4. Pedagogical clarity in the `drafts/` ch02 variants
- **Observation:** the live `ch02_state_matrix.md` is rigorous/complete but interleaves SU(2) algebra
  with its matrix expansions. The abandoned drafts (`ch02_dual_quaternion_draft.md`,
  `ch02_alt_two_layer.md`, `ch02_reference_architecture.md`) led with a *clearer conceptual frame*:
  the entire rigid-body state as one composable object; explicit "construct in SU(2) (singularity-free)
  / act via matrices" layering; cross-product = commutator, dot = trace, norm = determinant as native ops.
- **Recover (optional):** a 1-paragraph "conceptual summary / reading guide" at the top of §2 pointing
  at that framing. Cheap; improves onboarding. Lower priority than the math/audit items above.

## Reviewed and deliberately NOT recovered
- The deprecated `007_matched_pair_principle.md` was *shorter* than the live `ch03`, but the live
  chapter's rigor is the point; at most quote its plain-English precis as an intro (cosmetic).
- The `plan/` aspirational appendix/chapter stubs and the full 38-chapter book: archival is correct;
  recover only if that book is revived.
- `deprecated/` derivations: pedagogy-deprecated and **not citable** as sources (binding); kept for
  history only.

## Provenance
Agent-surveyed (structure/pedagogy only), hand-vetted. The agent over-valued some items (it tagged
4 of 9 as "recover"); items 1–2 here are the ones with real, low-cost, high-return recovery; 3–4 are
optional onboarding improvements.
