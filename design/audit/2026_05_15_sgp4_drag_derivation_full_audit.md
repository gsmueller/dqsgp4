# SGP4 Near-Earth Drag Derivation — Full-Session Multi-Pass Audit — INDEX

> Reconstructed 2026-06-01 into single-pass-readable parts (directive: documents
> required to be read in full must be split to manageable size). Full content
> (799 lines) now lives in `2026_05_15_sgp4_drag_derivation_full_audit/`. The
> pre-split monolith is preserved in git history at commit `f790297`.

Status @ 915b9c9: CLOSED D-2 / D-3 / D-7 / D-9 / D-10 · OPEN D-1 / D-4 / D-5 / D-6 / D-8 · NEW D-11 (main-doc §3 Thm 3.1 J₂-radial secular multiplier 3× too small — found by Phase 2.A §A.6; code-correct form `δr/r = −(3/2)(k₂/p²)β(3cos²i−1)`).
Phase 2.A (C₂ code-match) **COMPLETE** — `verify_phase2a.m` 11/11, `simplify(C₂_derived − C₂_code)=0` vs `drag_coefficients.h:146-149`; Phases 0-rev1 + 1 retrospectively validated. **D-9 CLOSED** (`I^(0,4)=(2+3η²)/(2(1−η²)^(7/2))` code-matched; main-doc Eq (4.4) = transcription error). **D-3 resolved** (correct β³ master orbit-average, §A.2). Trace: `design/derivations/sgp4_drag_phase2a_C2_trace.md`.

| Part | Scope |
|---|---|
| [01_standards_and_schedule.md](2026_05_15_sgp4_drag_derivation_full_audit/01_standards_and_schedule.md) | Standard; STANDING Standard 10 (code-match shall; numerical-at-points insufficient); five-pass schedule; pre-contamination checkpoint |
| [02_pass1_findings_sec0-sec4.md](2026_05_15_sgp4_drag_derivation_full_audit/02_pass1_findings_sec0-sec4.md) | PASS 1 §0-§4: D-1 (B* missing rho0); §4 orbit-averaging incl D-9 (I^(0,4)), D-10 (missing |v|, CRITICAL), D-2, D-3 (master orbit-avg Eq 4.2) |
| [03_pass1_findings_sec5-sec20.md](2026_05_15_sgp4_drag_derivation_full_audit/03_pass1_findings_sec5-sec20.md) | PASS 1 §5-§20: C2, C1, C3, C4, C5, D2-D4, t-cofs, RAAN coupling, long-period, drag corrections, error catalog; Pass-1 summary |
| [04_pass1.5_ocr_matching.md](2026_05_15_sgp4_drag_derivation_full_audit/04_pass1.5_ocr_matching.md) | PASS 1.5 OCR-matching audit; boxed-equation catalog + tally; re-derivation queue |
| [05_sanity_logging_deprecated_audit.md](2026_05_15_sgp4_drag_derivation_full_audit/05_sanity_logging_deprecated_audit.md) | Pre-Phase-0 sanity checks; logging convention; PASS A.6 deprecated-folder audit |
