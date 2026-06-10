# SGP4 Near-Earth Drag Derivation — Full-Session Multi-Pass Audit

**Session:** 2026-05-15 (R09 resumption)
**Auditor:** main session (dedicated audit worktree @ c25ff56)
**Subject:** `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` (832 lines, §0 … §20) and the
companion code in `src/atmosphere/drag_coefficients.h`, `src/sgp4/near_space.h`, `src/forces/drag.h`.

## Standard

Per user directive 2026-05-15 (session-start):

> *"All math work that has been derived has errors in it. Errors in derivation, errors
> in it being applied to our specific use as applied to the original SGP4, and now the
> different application DQSGP4, errors in sufficiency for the immediate task. You will
> be required to do several accuracy passes over each solution. These are not spot checks,
> these are deliberate tracked full-session audits of all code. Every line of proof must
> be viewed."*

This audit log is the **tracked artifact** for those passes. Each pass produces additional
finding rows and updated verdicts. No claim is "spot-checked PASS" — every line of every
proof must be traced.

## STANDING BINDING RULE (Standard 10, user directive 2026-05-15)

> *"The solver shall not use theoretical results that can not be matched to the code.
> Numerical testing at points is not acceptable."*

This is a **standing shall statement** — it applies to all current and future work in
this corpus. The rule has the following operational consequences:

1. **Theoretical result usage requires symbolic code-match.** Any theoretical
   result that is intended for use by the SGP4 implementation must be
   symbolically transformable, by a chain of named theorem applications and
   definition-based substitutions, to a closed-form expression algebraically
   equal to the corresponding code expression. The chain must be demonstrable
   step-by-step in the cleanroom proof format.

2. **Numerical testing at sample points is INSUFFICIENT.** Verifying that a
   derived expression evaluates to the same value as a code expression at
   sample points (e.g. η=0, η=0.5, e=0.01) is **not acceptable as proof**.
   Numerical agreement is necessary but not sufficient — it does not rule out
   coincidental agreement at sampled points while disagreeing at unsampled
   points, nor does it rule out a theoretical result computing a slightly
   different quantity that happens to agree numerically in the tested regime.

3. **Verifier scope (what counts as symbolic equality).** Octave + symbolic-pkg
   / SymPy `simplify(lhs − rhs) == 0` IS symbolic equality verification
   (acceptable). Octave / Python evaluating `lhs(η=0.5) − rhs(η=0.5)`
   numerically and comparing to a tolerance is NOT acceptable as proof of
   equality, only as a sanity check that the symbolic equality is plausible.

4. **Status of audit findings under this rule.** Findings D-1 through D-9
   (and any future D-N) are observations of inconsistency between the
   existing main doc and *something* (the code, the mathematical definition,
   the SGP4 application, or the re-derivation). A finding is **CLOSED** only
   when:
   - The discrepancy is resolved by a symbolic code-match across the
     affected derivation chain, OR
   - The finding is escalated to a documented model-choice or accuracy-bound
     change in the SGP4 specification with explicit per-step error bookkeeping.

5. **Status of intermediate theoretical results under this rule.** Phase 0
   (Gauss VE) and Phase 1 (Lane integrals) produce intermediate theoretical
   primitives that do not individually have direct code matches in
   `drag_coefficients.h`. They are PREREQUISITES for downstream Phases that
   DO produce code-matched results (the C-coefficients, D-coefficients,
   t-cofs, etc.). The validity of Phase 0 + Phase 1 is established **only
   retrospectively** when a downstream Phase symbolically derives a code
   expression and the chain lands cleanly.

6. **Phase 2.A is the critical gating task** before any further Phase work
   proceeds. Phase 2.A: symbolically derive the SGP4 code's C₂ expression
   from Phase 0 + Phase 1 inputs, with explicit step-by-step symbolic
   equality at each algebraic transformation. The outcome decides:
   - If the chain lands cleanly → Phases 0+1 are retrospectively validated;
     D-9 is closed as "main doc Eq (4.4) is a transcription error in a
     consistent-with-code direction".
   - If the chain fails to land → some step in Phases 0+1 is wrong or
     mis-applied to SGP4; the failure point identifies the actual error.

**No theoretical result is ACCEPTED into the corpus until it has either been
directly code-matched OR is an intermediate in a chain whose final step has
been code-matched.** The phrase "the corpus consumes Theorem X" downstream of
a Phase N is *gated* on a Phase ≥ N+1 producing a code-matched chain that
uses Theorem X.

## Five-pass schedule

| Pass | Scope | Output |
|---:|---|---|
| **1** | Internal-consistency / algebra: every line of every proof in §0…§20 read; algebraic chains re-derived; equation labels and citations cross-checked within the document. | List of D-#  /  N-#  /  Q-# findings. |
| **2** | Primary-source fidelity: every [B59] / [BH61] / [L65] / [LH79] / [SR3] citation cross-checked against the cited primary source page/equation. | Per-citation PASS / NOT-FOUND / MISMATCH. |
| **3** | Code fidelity: every "Code: <file>:<lines>" cross-reference cross-checked against the actual source file. | Per-cross-reference PASS / MISMATCH. |
| **4** | BH61 cleanroom theory-corpus consistency: the SGP4-drag derivation cross-checked against the BH61 cleanroom derivations and the SR3 transcript at `sgp4_references/.../derivation/` and `sgp4_references/hoots_roehrich_1980/hoots_roehrich_1980_math_derivation.md`. | Cross-corpus inconsistency list. |
| **5** | DQSGP4 application: how does each SGP4-near-earth derivation apply (or not) to the dual-quaternion / SE(3) path? Where does the gravity-acceleration convention (Wrench linear-slot = acceleration, not force) interact with the drag-rate-of-element-evolution? | DQ-application fitness verdict per coefficient. |

After Passes 1-5 produce findings, the FINDING triage step assigns each finding
one of:

- **DERIVATION-ERROR** — the proof / algebra is wrong; must be re-derived
- **APPLICATION-ERROR** — proof is correct in isolation, but the application to SGP4 (or DQSGP4) misuses it
- **SUFFICIENCY-GAP** — proof needs additional content (e.g. R09's §11a/§9a/§4a/§12a) but what's there is not wrong
- **TRANSCRIPTION-ERROR** — typo / formatting in the document; math is right
- **NOTATIONAL** — confusing but not incorrect

## Pre-contamination checkpoint (per feedback_pre_contamination_checkpoint.md)

- **HEAD before audit:** `c25ff56` (dedicated audit worktree)
- **Derivation doc HEAD:** 832 lines, last touched in commit `94b7584` ("Add theoretical basis document for SGP4 near-earth drag coefficients")
- **Code file HEAD:** `src/atmosphere/drag_coefficients.h` last touched in commit `03b620a` (phantom-C₂ fix). Per the AUD-TBA audit card the file is PASS-with-notes; this audit re-opens that verdict.

If any finding rises to DERIVATION-ERROR severity, that closes the prior PASS-with-notes verdict and the file is re-classed as **AUDIT-OPEN**.

---

