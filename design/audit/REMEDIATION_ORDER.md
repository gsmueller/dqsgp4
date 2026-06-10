# Remediation Order — SGP4 Drag Re-derivation Work Plan

> **SUPERSEDED 2026-05-16** by a successor plan, itself since superseded; see `design/DQSGP4_EXECUTION_PLAN.md`.
> This document was the initial flat 45-item plan; the user rejected it as "pretty
> weak" and the proper plan-mode plan replaces it. The canonical work order and
> execution pointer now live in the plan file. This file is retained as a
> historical record only — no further updates.

**Created:** 2026-05-15
**Authority:** User directive 2026-05-15 — "go through a proper planning cycle. Create an order for each item that needs to be done, and then only do one at a time."

This is the canonical **ordered work plan** for the SGP4 drag re-derivation
project. Items must be executed **one at a time, in order**. Each item is
small enough to be a single focused work cycle (1–3 tool calls). After each
item, I yield to the user for review before proceeding to the next.

If an item surfaces a finding that re-orders subsequent items, the plan is
**revised** before continuation — the revision itself becomes an item.

## Status legend

- ☐ pending
- ⧖ in progress (only ONE item may be ⧖ at a time)
- ✓ complete (with commit hash)
- ⊘ blocked (waiting on a prerequisite)
- ✗ failed (must be re-planned)

---

## Block A: Phase 0-rev1 (fix D-10 missing-|v| error)

D-10 is the critical Phase 0 §0.4 error: the drag specialization wrote
`F_drag = -B*·ρ·v` (linear in v) instead of the correct quadratic Newtonian
drag `F_drag = -B*·ρ·|v|·v`. All Phase 0 §0.4 closed forms are off by a `|v|`
factor. Must be fixed before any Phase ≥ 2 work can proceed.

| # | Item | File(s) | Success criterion | Status |
|---:|---|---|---|---|
| 1 | Write this plan | `design/audit/REMEDIATION_ORDER.md` | Plan file exists, ordered, approved | ✓ (this turn) |
| 2 | Theorem 0.4.1 rewrite — quadratic drag postulate + `|v|` in R/T components + corrected proof | `sgp4_drag_phase0_foundations.md` §0.4 Theorem 0.4.1 | Statement reflects `F_drag = -B*·ρ·|v|·v`; proof derives `R_drag = -B*·ρ·|v|·ṙ`, `T_drag = -B*·ρ·|v|·r ḟ`, `N_drag = 0` with `|v| = √(ṙ² + (rḟ)²)` | ☐ |
| 3 | Theorem 0.4.2 (0.4.2.1) — `ȧ_drag` recompute with corrected drag | `sgp4_drag_phase0_foundations.md` §0.4 Theorem 0.4.2 (0.4.2.1) | Closed form `ȧ_drag = -(2 n B* ρ a²/β³)·(1+e²+2e cos f)^{3/2}`; derivation shows where `|v|·v = |v|² · v̂ = v² · v̂` (or equivalently `(1+e²+2e cos f)^{3/2}`) emerges | ☐ |
| 4 | Theorem 0.4.2 (0.4.2.2) — `ė_drag` recompute | same file, §0.4 (0.4.2.2) | Closed form derived consistently with corrected `R_drag`, `T_drag` | ☐ |
| 5 | Theorem 0.4.2 (0.4.2.3) — `ω̇_drag` recompute | same file, §0.4 (0.4.2.3) | Closed form derived; `−Ω̇ cos i` correction still vanishes (N=0) | ☐ |
| 6 | Theorem 0.4.2 (0.4.2.4) — confirm `Ω̇_drag = i̇_drag = 0` still | same file, §0.4 (0.4.2.4) | One-sentence proof that N=0 ⇒ both zero (unchanged from current) | ☐ |
| 7 | Theorem 0.4.2 (0.4.2.5) — `Ṁ_drag` recompute | same file, §0.4 (0.4.2.5) | Closed form derived consistently | ☐ |
| 8 | Alignment remarks revision — §0.3.1 (drag preview), §0.4 alignment-to-SGP4 / alignment-to-implementation remarks for both theorems | same file, §0.3.1, §0.4 alignment remarks | Remarks reference the corrected forms; explicit acknowledgement that D-10 was found by Standard 10 | ☐ |
| 9 | §0.5 Theorem 0.5.3 — orbit-averaged drag rate setup with new integrand | same file, §0.5 Theorem 0.5.3 | Integrand updated to `(1+e²+2e cos f)^{3/2} · ρ(r(f))`; Phase 2.A handoff annotated | ☐ |
| 10 | §0.7 closed-status summary update — add 0.4.1/0.4.2 corrections; remove "D-2 closed" from §0.2.5 if invalidated, etc. | same file, §0.7 | Summary reflects the rev1 state; Phase 0 marked as "rev1, with D-10 fixed" | ☐ |
| 11 | `verify_phase0.m` — add F_drag-vs-Postulate-1.2 check (the check that should have caught D-10) | `verify_phase0.m` | New check: `R_drag·r̂ + T_drag·t̂ + N_drag·n̂` ≡ `-B*·ρ·|v|·v` symbolically | ☐ |
| 12 | `verify_phase0.m` — update existing Theorem 0.4.1 ratio check to reflect that `|v|` factor cancels in ratio (unchanged), and add new check for the closed forms with `(...)^{3/2}` factor | same file | Existing 0.4.1 ratio check still passes; new closed-form check for ȧ_drag passes | ☐ |
| 13 | Run `verify_phase0.m` and confirm all PASS (now ~30 checks expected) | (run only) | Octave: "X / X PASS" with X ≥ 28 | ☐ |
| 14 | Commit Phase 0-rev1 | (git) | Single commit with message describing D-10 fix; commit hash recorded | ☐ |

## Block B: Phase 2.A (symbolic trace of SGP4 C₂ from corrected Phase 0)

Phase 2.A is the **critical gating Phase** that determines (a) whether the
Phase 0 + Phase 1 framework actually produces the SGP4 code's C₂, and
(b) the classification of D-9 (preliminary).

| # | Item | File(s) | Success criterion | Status |
|---:|---|---|---|---|
| 15 | Create Phase 2.A document skeleton (header, scope, sections list) | `design/derivations/sgp4_drag_phase2a_C2_trace.md` (new) | File exists with §2A.0 scope + §§2A.1-2A.9 stubs | ⊘ on Block A |
| 16 | Phase 2.A §2A.1 — Setup orbit-averaged ȧ_drag integral from corrected Phase 0 | same file, §2A.1 | Integral expressed in `(f, e, a, β, n, ρ(r), |v|)` with the `(1+e²+2e cos f)^{3/2}` factor explicit | ⊘ |
| 17 | Phase 2.A §2A.2 — Substitute Lane density model ρ(r) = ρ_0 ((q-s)/(r-s))^4 | same file, §2A.2 | `(r-s)^{-4}` substituted in terms of `(f, e, a, s)` | ⊘ |
| 18 | Phase 2.A §2A.3 — Apply Lane f† substitution; state the O(eη) approximation explicitly | same file, §2A.3 | `r-s = (1/ξ)(1-η cos f†)`; approximation error bounded; integration variable change | ⊘ |
| 19 | Phase 2.A §2A.4 — Reduce to Lane integrals; identify which `I^{(p,m)}` appear | same file, §2A.4 | Specific `(p, m)` values identified, mapped to Phase 1 Theorem 1.3.2 closed forms | ⊘ |
| 20 | Phase 2.A §2A.5 — Substitute Phase 1 Theorem 1.3.2 closed forms; assemble Part A (Keplerian drag) closed form | same file, §2A.5 | Part A: `a₀·(1 + ... )` polynomial derived | ⊘ |
| 21 | Phase 2.A §2A.6 — Add J₂ secular density perturbation Part B (from Phase 0 §3.2) | same file, §2A.6 | Part B: `(3/8)J₂·ξ·ψ^{-2}·(3cos²i-1)·(...)` polynomial derived | ⊘ |
| 22 | Phase 2.A §2A.7 — AFGP4 → SGP4 simplification (drop O(e²) Part A, drop O(e) Part B); explicit accuracy bounds | same file, §2A.7 | SGP4 Part A: `(1 + (3/2)η² + 4eη + eη³)`; SGP4 Part B: `(8 + 24η² + 3η⁴)`; dropped terms catalogued with order | ⊘ |
| 23 | Phase 2.A §2A.8 — Code-match witness: step-by-step symbolic equality with `drag_coefficients.h:146-149` | same file, §2A.8 | Each algebraic transformation from derivation → code expression is shown explicitly | ⊘ |
| 24 | `verify_phase2a.m` — Octave + symbolic-pkg verifier that runs the derivation chain and asserts equality with code expression | `verify_phase2a.m` (new) | Symbolic equality between derivation and code C₂ formula (parsed from header file) | ⊘ |
| 25 | Run `verify_phase2a.m` and confirm symbolic equality | (run only) | All checks PASS | ⊘ |
| 26 | Audit log update — D-9 classified (closed/escalated based on Phase 2.A outcome) | `2026_05_15_sgp4_drag_derivation_full_audit.md` D-9 section | D-9 status changed to CLOSED with commit hash citation OR ESCALATED with explicit re-derivation task | ⊘ |
| 27 | Commit Phase 2.A | (git) | Single commit | ⊘ |

## Block C: Phase 1 §1.4 (Lane integrals for p ≥ 1)

Gated on Phase 2.A outcome. If Phase 2.A only needs `I^{(0,m)}` already
derived in Phase 1 §1.3, then §1.4 may not be needed at all. If §2.A
needs `I^{(p,m)}` for `p ≥ 1`, then §1.4 derives them.

| # | Item | Status |
|---:|---|---|
| 28 | Determine `I^{(p,m)}` needs from Phase 2.A | ⊘ on Block B |
| 29 | Phase 1 §1.4 — `I^{(p,m)}` for `p ≥ 1` (if needed) | ⊘ |
| 30 | Phase 1 §1.5 — Residue calculus alternative path (second verification) | ⊘ |
| 31 | `verify_phase1.m` mechanical verifier | ⊘ |
| 32 | Commit Phase 1 finalization | ⊘ |

## Block D: Remaining phases

Each subsequent phase has its own ordered item-list — to be added when the
prior phase is complete and the actual subsequent-phase requirements are
known (since they may depend on the resolution of upstream phases).

| # | Item | Status |
|---:|---|---|
| 33 | Plan Phase 3 (C₃, J₃ coupling) — items TBD when Phase 2 complete | ⊘ |
| 34 | Plan Phase 4 (C₄) — items TBD | ⊘ |
| 35 | Plan Phase 5 (C₅) — items TBD | ⊘ |
| 36 | Plan Phase 6 (D₂) — items TBD | ⊘ |
| 37 | Plan Phase 7 (D₃/D₄) — items TBD | ⊘ |
| 38 | Plan Phase 8 (t-cofs) — items TBD | ⊘ |
| 39 | Plan Phase 9 (xnodcf) — items TBD | ⊘ |
| 40 | Plan Phase 10 (omgcof / xmcof / delmo / sinmo) — items TBD | ⊘ |
| 41 | Plan Phase 11 (Lane f† formalization) — items TBD | ⊘ |
| 42 | BH61 cleanroom dispatch — required by Phase 3 + Phase 9 | ⊘ |

## Block E: End-of-session

| # | Item | Status |
|---:|---|---|
| 43 | Memory writeback — record D-10 finding, Standard 10, Phase 0-rev1 outcome | ⊘ |
| 44 | Worktree management — merge the audit worktree to `session/2026-04-23` per policy | ⊘ |
| 45 | Final regression test — SGP4-VER 23/33 satellites (unchanged baseline; documentation-only changes should not alter test results) | ⊘ |

---

## Execution discipline

- **One item at a time.** I do not move to item N+1 until item N is complete AND user-approved.
- **After each item:** brief status surface to user with what was done and what comes next. User can approve, redirect, or revise the plan.
- **No batch optimization.** Even when two items look like they could be combined, I keep them separate. Standard 10 demands that each algebraic transformation be visible and reviewable.
- **Plan revision rule.** If an item surfaces a finding that invalidates a later item, I STOP and re-plan before proceeding. The plan revision is itself an item.
- **Verification rule.** Per Standard 10: every item that produces a theoretical claim must have a symbolic-equality check before the item is marked complete. Numerical agreement at sample points is sanity-check-only.

---

## Current execution pointer

- **Item 1 (this plan):** ✓ complete
- **Item 2 (Theorem 0.4.1 rewrite):** ☐ — next, awaiting approval to proceed
- **Items 3..45:** ☐ / ⊘ — gated on Item 2 and downstream
