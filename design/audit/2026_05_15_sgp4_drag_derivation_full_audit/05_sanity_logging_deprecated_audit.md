## Pre-Phase-0 sanity checks

Before any re-derivation:
1. **Verify BH61 cleanroom paths.** Confirm the BH61 cleanroom files cited in §20 (`ch10c_secular_average.md`, `ch11b_S1star_partials.md`, `ch11d_secular_rates.md`) exist and contain what (3.1), (7.1), (14.2), (14.3) need.
2. **Verify deprecated/020 status.** §5 currently delegates the full C₂ derivation to `deprecated/020_c2_drag_integral_derivation.md`. If 020 has been retired (per the `deprecated/` location), Phase 2 must absorb its content; otherwise, Phase 2 may cite it.
3. **Source policy.** Per `feedback_born_digital_latex.md`: where Wikipedia / ProofWiki / nLab covers a result (Kepler geometry, residue calculus, binomial series), cite that. Where born-digital astrodynamics (Lara 2021, Sneeuw 2022) covers a result, prefer that. Only as a last resort consult the OCR'd primary sources [B59], [LH79], [SR3].

## Logging convention

When a finding is **fixed**, the row above gets a `Status: FIXED in commit <hash>` line. When a finding is **DECLINED** (e.g. determined to be a primary-source quirk, not a real error), the row gets `Status: DECLINED — <rationale>`. No finding is removed.

When a re-derivation Phase is **completed**, the Phase row gets a `Status: COMPLETE in commit <hash> + verifier <path>` line. Each Phase produces a markdown file with the symbolic derivation and (per BH61 cleanroom style) an Octave or SymPy verifier that mechanically checks the result.

---

## PASS A.6 — Deprecated-Document Audit (user directive 2026-05-15)

> *"Files in the deprecated folders have various levels of quality. There should be a pass over the top of each document to ensure that it is symbolically verified and traced to an accurate proof."*

`design/derivations/deprecated/` contains **28 files** — 25 numbered (000…025) plus 3 ch02 variants (state-matrix DCM draft iterations). Several of these are cited by the main `sgp4_near_earth_drag_theoretical_basis.md` (e.g. §5 delegates the full C₂ derivation to `deprecated/020`). Per the directive, each deprecated file must be audited the same way the main derivation doc was: every line of proof read; every "matched to OCR" replaced with a symbolic re-derivation traced to a theoretical source.

### Per-file audit catalog

Priority A = on the critical path for current R09 / SGP4 drag work; must be audited **before** the corresponding Phase in Pass 1.5 can complete.
Priority B = related to SGP4 broader work; audit after Phase 0-2 of the drag re-derivation is settled.
Priority C = other modules (astronomy, coord transform, state matrix); audit independently.

| File | Priority | Cited by main doc | Subject |
|---|---|---|---|
| `003_geometric_constants.md` | **A** | (constants table; §0.2) | Geometric primitives used throughout. |
| `004_physical_constants.md` | **A** | (constants table; §0.2) | Physical constants used throughout. |
| `010_brouwer_averaged_hamiltonian.md` | **A** | implicit (Brouwer foundation for §3, §13) | The Brouwer 1959 secular Hamiltonian F₂*. |
| `011_brouwer_secular_rate_polynomials.md` | **A** | §13 (drag-RAAN coupling) | Secular rates Ṁ, ω̇, Ω̇. |
| `012_lane_hoots_drag_derivation.md` | **A** | §0.4 (overview cite) | Top-level Lane-Hoots drag model overview. |
| `019_short_period_corrections.md` | **A** | §14 (xlcof algebra) | Short-period J₂ corrections. |
| **`020_c2_drag_integral_derivation.md`** | **A (top)** | **§5 delegates entire C₂ derivation here** | Full C₂ from Lane integrals. The most critical deprecated file for current work. |
| `021_long_period_corrections.md` | **A** | §14 (long-period J₃) | Long-period J₃ — informs aycof/xlcof. |
| `022_density_model.md` | **A** | §2 (referenced) | Lane density model. |
| `000_three_errors.md` | B | — | Three-error framework. |
| `001_q0_series.md` | B | possibly §2 (atmosphere) | q₀ series expansion. |
| `002_q0_prime_series.md` | B | possibly §2 (atmosphere) | q₀' series expansion. |
| `005_R2_series.md` | B | — | R² series. |
| `006_J2n_from_ellipsoid.md` | B | (§2 perigee adjustment indirectly) | J₂n from ellipsoid parameters. |
| `007_matched_pair_principle.md` | B | — | Matched-pair principle (likely numerical-method). |
| `008_deep_space_constants.md` | B | — | Deep-space constants. |
| `009_sgp4_modified_kepler.md` | B | (R12 audit covers code) | Modified Kepler iteration. |
| `013_element_recovery.md` | B | (R08 audit covers code) | Brouwer-Lyddane element recovery. |
| `018_resonance.md` | B | (R01 audit covers code) | Deep-space tesseral resonance. |
| `023_coordinate_transformation.md` | B | — | Coordinate transformation. |
| `024_numerical_guards.md` | B | — | Numerical guards / cutoffs. |
| `014_sidereal_time.md` | C | — | IAU 1982 sidereal time. |
| `015_solar_ephemeris.md` | C | — | Solar ephemeris (Meeus). |
| `016_lunar_ephemeris.md` | C | — | Lunar ephemeris (Meeus). |
| `017_third_body_perturbation.md` | C | — | Third-body perturbation (Kozai). |
| `025_astronomical_constants.md` | C | — | Astronomical constants. |
| `ch02_state_matrix_dcm_complete.md` | C | — | State-matrix DCM (DRAFT-3). |
| `ch02_state_matrix_v1_dcm.md` | C | — | State-matrix DCM (DRAFT-1). |
| `ch02_state_matrix_v2_pre_dual_velocity.md` | C | — | State-matrix DCM (DRAFT-2). |

**Tally:** 9 Priority-A, 12 Priority-B, 7 Priority-C.

### Per-file audit method (each file gets its own 4-pass mini-audit)

For each deprecated file, in the order of Priority A → B → C:

- **A.6.P1** Read every line. Identify boxed equations / theorems. Catalog status as DEF / POST / PRIM / BH61-CR / OCR-MATCH.
- **A.6.P2** Internal consistency: for each OCR-MATCH equation, does the file include a symbolic derivation or is it just transcription?
- **A.6.P3** Code fidelity: for each numerical claim, cross-check against the corresponding source file (where applicable).
- **A.6.P4** Triage: if symbolically derived → KEEP-CITE; if not derived → mark as needs-re-derivation; if obsolete/redundant → mark for deletion.

The output per file is an audit row in this table below (table is populated lazily as each Priority-A file is audited):

| File | Boxed equations | OCR-MATCH count | Symbolic-derivation count | Verdict | Audit notes |
|---|---:|---:|---:|---|---|
| `020_c2_drag_integral_derivation.md` | 22 (020.Eq.1 … 020.Eq.22) | **9 (self-tagged `[UNVERIFIED]`)** | 1 (only the (3/8) J₂ rearrangement at Eq.16-17) | **NEEDS-RE-DERIVATION** | Honest disclosure — the author explicitly inserts "`[UNVERIFIED — transcribed from scanned PDF: Lane_Hoots_1979 p. X]`" warnings ABOVE every transcribed equation. The "BUILD" section interprets the structure but does **not** evaluate the orbit-averaged integrals. The single piece of actual algebra is the trivial `(3/2)·k₂·A = (3/8)·J₂·(3cos²i-1)` rearrangement. Cross-check against main doc §5/§8: 020 Eq.18 (C₄) is algebraically identical to main §8 modulo equivalent factoring; 020 Eq.21 (D₃) matches main §11.1 expansion. Status: 020 is **transcription-with-cross-check, not derivation**. |
| `012_lane_hoots_drag_derivation.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `022_density_model.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `010_brouwer_averaged_hamiltonian.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `011_brouwer_secular_rate_polynomials.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `003_geometric_constants.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `004_physical_constants.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `019_short_period_corrections.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| `021_long_period_corrections.md` | (TBD) | (TBD) | (TBD) | (TBD) | (TBD) |
| ... (Priority B / C lazily filled) | | | | | |

### Key 020 observations (impact on Phase 2 of the re-derivation queue)

1. **020 is honest and explicitly OCR-transcription-only.** Every transcribed equation carries an `[UNVERIFIED — transcribed from scanned PDF]` warning. The author did not claim symbolic derivation.
2. **Phase 2 of the re-derivation queue cannot delegate to 020.** Phase 2 must produce a from-scratch symbolic derivation of C₂ from Phase 0 (Kepler primitives) + Phase 1 (Lane integrals via residue calculus). The 020 transcript may be referenced as a **structure guide** (Part A vs Part B decomposition is useful), but every equation must be independently derived.
3. **020 contains the AFGP4 → SGP4 simplification record** (020.Eq.7 full form → 020.Eq.12-13 simplified). This record IS useful for the "error sources" inventory but is itself transcribed from LH79.
4. **Other coefficient transcriptions** (C₃ at 020.Eq.17, C₄ at 020.Eq.18, C₅ at 020.Eq.19, D₂-D₄ at 020.Eq.20-22) match the main doc's §7-§11 — same OCR-MATCH inheritance, same NEEDS-RE-DERIVATION verdict.
5. **There is no Brouwer-side derivation in 020.** It treats Brouwer J₂ secular results as given. The Brouwer foundation is in 010 / 011 (still to be audited).

---

### User triage 2026-05-15 — entire deprecated/ folder is pedagogy-deprecated

> *"the truth is they are deprecated because they all share the same pedagogy."*

The deprecation is itself the verdict. All 28 files in `design/derivations/deprecated/`
share the **same flawed pedagogy** that 020 made visible: OCR-transcription-with-cross-check
masquerading as a derivation document, with the result that no equation is actually
symbolically derived from first-principles theoretical sources. Per-file audit of the
remaining 27 would only re-confirm this pattern; the user has saved that work by
making the folder-level classification explicit.

#### Implications

- **No deprecated file may be cited as a "derivation source"** under the new standard.
- The deprecated/ folder is retained as a **historical OCR-transcription record** and as a **structure / sketch reference** (the Part-A/Part-B decomposition in 020, the AFGP4→SGP4 simplification record, etc. are useful organizational aids).
- All re-derivations from this audit onward must use the **correct pedagogy**: BH61-cleanroom style with named theorems / lemmas, (T)/(D)-labeled steps, full proof chains (every step shown, none skipped), born-digital theoretical sources (Wikipedia / ProofWiki / nLab / Battin / Roy / Vallado-LaTeX / Lara 2021 / Sneeuw 2022), and a mechanical (Octave / SymPy) verifier alongside each derivation file.
- The Brouwer-side foundations (010, 011, the BH61 J₂ secular and S₁ long-period generator results) — when needed — are pulled from the **BH61 cleanroom** at `sgp4_references/vallado_celestrak/documentation/SGP4/Brouwer_Hori 1961/derivation/` via sub-agent dispatch per `feedback_consolidator_no_BH61_reads.md`, not from the deprecated/ folder.
- Phase 0 of the re-derivation queue starts now. It draws only on born-digital primitives (Kepler geometry, Gauss-Lagrange VE) — no OCR, no deprecated, no BH61 dispatch required.

#### Single-row audit verdict for the folder

| File / folder | Verdict | Reason |
|---|---|---|
| `design/derivations/deprecated/` (28 files) | **PEDAGOGY-DEPRECATED — do not cite as derivation source** | Per user 2026-05-15. All files share the OCR-transcription-with-cross-check pedagogy that 020 made visible. Retained as historical record. |

### Audit-ordering rationale

`020_c2_drag_integral_derivation.md` is **audited first** because:
1. The main `sgp4_near_earth_drag_theoretical_basis.md` §5 delegates the entire C₂ derivation to it.
2. If 020 is symbolically rigorous, Phase 2 of the main re-derivation queue can cite it.
3. If 020 is OCR-matching-only, Phase 2 must absorb it and re-derive from scratch.

This audit is on the critical path; no Phase-2 work can begin until 020's status is decided.

After 020, the order for Priority-A:
- 022 density model (Phase 0 needs it)
- 012 Lane-Hoots overview (orientation)
- 010 / 011 Brouwer (Phase 0 + Phase 6/9 need Brouwer secular & long-period generator)
- 003 / 004 constants (used throughout)
- 019 / 021 short / long period (Phase 7 needs these for xlcof/aycof)

### Caveat — BH61 read restriction

Per `feedback_consolidator_no_BH61_reads.md`, the main session cannot directly read BH61-* / Brouwer* / Hoots* / Lara* / ch05+ files / primary PDFs. **But files in `design/derivations/deprecated/` are project files**, not BH61 cleanroom files. They appear to be in-house derivations and may be read directly. Per content, however, any `ch10+` filenames or BH61-result quotes would still be off-limits. If a deprecated file turns out to quote BH61 wholesale, the audit pivots to sub-agent dispatch for that file's audit.
