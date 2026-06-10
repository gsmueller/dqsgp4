# Draft Plan: Appendix D — Source Document Index

## Objectives

1. Provide a complete, annotated bibliography of every source document used in the textbook.
2. Classify each source as: (a) clean digital (post-LaTeX, directly typeset); (b) scanned historical (OCR-prone; all values cross-verified); (c) reference implementation (SR3 code, Vallado code — treated as hypothesis, not proof).
3. For each source: give the notation map to textbook notation, the specific sections used, and any known transcription errors found during cross-verification.
4. Establish provenance for every derived coefficient.

## Section Structure

### §D.1 Clean Digital Sources

This section provides annotated entries for post-LaTeX, directly-typeset sources with their verification status and notation differences from textbook notation.

Stub: Annotated entries for: Lara (2021) — "Brouwer's satellite problem revisited" (Journal of the Astronautical Sciences); Sneeuw (2022) — geodesy lecture notes (direct-typeset PDF); Na et al. (2012) — drag coefficient derivation. For each: full citation, availability, sections used in textbook, notation differences from textbook notation, status (verified, partially verified, not yet cross-checked).

### §D.2 Historical Sources: Verified

This section provides annotated entries for scanned historical sources, documenting known OCR issues and the cross-verification sources used to validate each value.

Stub: Annotated entries for: Brouwer (1959) — "Solution of the Problem of Artificial Satellite Theory Without Drag" (Astronomical Journal); Lane and Hoots (1979) — "General Perturbation Theories Derived from the 1965 Lane Drag Theory"; Hoots and Roehrich (1980) — "Models for Propagation of NORAD Element Sets" (SR3); Kaula (1966) — "Theory of Satellite Geodesy"; Heiskanen and Moritz (1967) — "Physical Geodesy"; Aoki et al. (1982) — "The new definition of Universal Time." For each: full citation, ADS/DOI location, sections used, known OCR issues (if scanned), notation map to textbook. Explicit flag: values from scanned sources are cross-verified against at least one clean digital source or independent derivation.

### §D.3 Reference Implementations: Not Cited as Proof

This section documents reference implementations (SR3, Vallado code, Celestrak tools) with an explicit statement that these are treated as hypotheses to be tested against independent derivations, not as proofs.

Stub: Entries for: SR3 (Hoots, Schumacher, Glover 2004 FORTRAN); Vallado MATLAB/C++ implementation; Celestrak online tools. Explicit statement: these implementations are treated as hypotheses to be tested, not as proofs. Coefficients appearing in these codes are derived from first principles in the textbook chapters; matching is verified but not used as the primary justification. Cross-reference feedback_derivation_first and feedback_vet_all_claims.

### §D.4 Notation Map

This section provides the global notation correspondence table between the textbook and the four major primary sources, enabling readers to cross-reference equations directly.

Stub: Table D.4.1: global notation correspondences between the textbook and the four major sources (Brouwer 1959, Hoots and Roehrich 1980, Lara 2021, Kaula 1966). Rows: each symbol with its meaning and its name in each source. This table is the key to reading the source documents alongside the textbook.

### §D.5 Known Issues Log

This section maintains a running log of documented transcription errors, OCR artifacts, and notation inconsistencies discovered in source documents during the textbook's development.

Stub: Table D.5.1: documented transcription errors, OCR artifacts, and notation inconsistencies found in source documents during the textbook's development. For each: source, location (equation or table number), issue description, textbook resolution. Updated as new issues are found.

**Pre-populated entries (known at Draft stage):**

| Source | Location | Issue | Textbook Resolution | Ref Tag |
|--------|----------|-------|--------------------|----|
| Hoots & Roehrich (1980) SR3 | Multiple formulas throughout | Known transcription and coding errors in the SGP4 equations, documented by Vallado et al. (2006) AIAA 2006-6753 | All SR3 formulas re-derived from first principles; Vallado et al. (2006) used as secondary check (not primary source); each discrepancy documented in-chapter | [A.27.4], [A.35.x] |
| Brouwer (1959) | Long-period inclination rate formula | Sign error identified by Lara (2021, J. Astronaut. Sci. 68) in the long-period inclination rate; original paper gives wrong sign in one coefficient | Use Lara (2021) corrected formula as cross-check; the textbook derives from the Hamiltonian and verifies both signs | [A.16.x], [A.19.x] |
| Chapront-Touzé & Chapront (1988) ELP 2000-82 | Longitude series coefficients | Published corrections by Chapront-Touzé & Chapront (2003, A&A 404) affect several series coefficients including terms used in the SGP4 lunar ephemeris | Ch 26 coefficients verified against the 2003 corrected version; original 1988 paper not used as primary source | [A.26.4] |
| Aoki et al. (1982) | Eq. (14) GMST polynomial | Superseded for precise applications by IAU 2006 Earth Orientation Model (Capitaine et al. 2003, A&A 412); the 1982 formula is retained in SGP4 as part of the matched pair | SGP4 matched-pair implementation uses 1982 formula; §29.3 notes the IAU 2006 alternative | [A.29.2] |
| Lane (1965) drag theory | Drag coefficient formulas | Errors in the original 1965 formulas corrected by Lane & Hoots (1979) | Ch 22 uses the 1979 corrected formulas; the 1965 paper is not cited as primary | [A.22.4] |
| Kaula (1966) | Table 1 ($F_{lmp}$ inclination functions) | $(l,m,p)=(4,2,0)$: $\sin i$ should be $\sin^2 i$; $(l,m,p)=(4,2,2)$: sign error; confirmed in multiple implementations | Derive $F_{lmp}$ analytically from Kaula's closed-form expression; do not transcribe from Table 1 | [A.15.4] |
| SR3 (1980) | DPPER Lyddane switching logic | Incorrect conditional for direct/Lyddane choice causes inclination to go negative during propagation | Use Vallado et al. (2006) AFSPC-mode conditional: if `nodep < 0` and `opsmode == 'a'`, add $2\pi$ | [A.35.5]–[A.35.7] |
| SR3 (1980) | DSPACE integrator step direction | Sign error in step logic for backward propagation causes jump at $t = -720$ min for GEO/Molniya | Use Vallado et al. (2006) corrected integrator | [A.35.5] |
| SR3 (1980) | SAVTSN 30-min shortcut in DPPER | Skipping lunar-solar terms when $|\Delta t| < 30$ min causes discontinuous ephemeris | Shortcut removed; recompute on every call | [A.35.6] |
| Jacchia (1970) | Semi-annual density variation | Modeled as temperature-driven; Jacchia (1971) showed this is incorrect | Use J71 or J77; do not use J70 for semi-annual variation | [A.21.5] |
| Lieske et al. (1977) | IAU 1976 precession constant | Approximately 3 mas/yr too high; superseded by IAU 2000 (Capitaine et al. 2003) | For TEME-to-ICRF better than 0.1 arcsec: use IERS Conventions 2003 | [A.30.3] |

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| All chapters | Specific citations | Each chapter's source references collected here |
| App C | Error tags | Tags associated with values from scanned sources |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|

(None)

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 1 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 0 |
| Examples | 0 |
| Error Notes | 0 |
| Equations | 0 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §D.1 Clean Digital Sources | Draft | |
| §D.2 Historical Sources: Verified | Draft | |
| §D.3 Reference Implementations: Not Cited as Proof | Draft | |
| §D.4 Notation Map | Draft | |
| §D.5 Known Issues Log | Draft | |
