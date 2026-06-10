# Draft Plan: Appendix C — Code-to-Theorem Mapping

## Objectives

1. For every source module (`src/*.h`, `src/*.cpp`), list the theorems, definitions, and lemmas it implements.
2. For every theorem in the textbook that has a computational realization, identify the implementing code location.
3. List all error tags (by [M/P/A.ch.n]) that apply to each module.
4. Serve as the primary audit tool: a reviewer can check that every code path corresponds to a derived result.

## Section Structure

### §C.1 Organization

This section explains the two-way cross-reference structure and the audit workflow for a reviewer checking code against derivations.

Stub: Two-way cross-reference. Part C.2 is organized by source file. Part C.3 is organized by theorem number. Each entry gives: source file or theorem number, counterpart reference, and a one-line description of the correspondence.

### §C.2 Source File Index

This section provides a table indexed by source file, listing the theorems implemented and the error tags applicable to each.

Stub: Table with one row per source file. Columns: file, description, implements (theorem list), error tags applied. Files: `angles.h`, `binomial_series.h`, `density_model.h`, `drag_coefficients.h`, `element_recovery.h`, `kepler.h`, `modified_kepler.h`, `near_space.h`, `deep_space.h`, `precomputed.h`, `resonance.h`, `secular_update.h`, `series.h`, `sgp4_propagator.h`, `sidereal_time.h`, `state_from_elements.h`, `state_vector.h`, `tle_parser.h`, and the math library headers (`factorial.h`, `wallis.h`, etc.).

### §C.3 Theorem Index

This section provides a table indexed by theorem identifier, listing the implementing source file(s) and code symbol names.

Stub: Table with one row per theorem. Columns: theorem identifier (e.g., Thm 16.4.1), description, implementing file(s), code symbol name if applicable, error tags. Organized by chapter.

### §C.4 Error Tag Summary

This section provides a complete alphabetically-ordered list of all error tags with their inline meanings and whether each bound is quantified or qualitative.

Stub: Complete list of all [M.ch.n], [P.ch.n], [A.ch.n] tags in the textbook, organized by chapter. For each tag: the inline meaning, the footnote bound, and whether the bound is quantified or qualitative.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1–5 (foundation) | Error framework, matched pair, series | Core theorems mapped to math library |
| Ch 8–10 (Keplerian) | Kepler equation, orbital elements | Mapped to `kepler.h`, `angles.h` |
| Ch 14–22 (perturbations) | Gravity, drag, Brouwer corrections | Mapped to perturbation source files |
| Ch 25–35 (pipeline) | Ephemerides through propagation | Mapped to pipeline source files |
| Ch 36–38 (architecture) | Propagator, constants, output | Mapped to `sgp4_propagator.h`, `precomputed.h`, `state_vector.h` |
| App B | Symbol names | Symbol-to-code variable correspondence |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|

(None -- error tags are catalogued in §C.4 but originate in their respective chapters.)

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 1 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 0 |
| Examples | 2 |
| Error Notes | 0 |
| Equations | 0 |
| Sections | 4 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §C.1 Organization | Draft | |
| §C.2 Source File Index | Draft | |
| §C.3 Theorem Index | Draft | |
| §C.4 Error Tag Summary | Draft | |
