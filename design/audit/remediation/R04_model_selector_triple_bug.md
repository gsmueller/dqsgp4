# R04 — model_selector triple bug

**Status:** PARTIAL — Halley naming + GM σ fixed; J5..J9 EGM2008 re-derivation TODO'd (commit 3cc618f)
**Severity:** P1
**Estimated scope:** ~30 minutes

---

## Files

**Write:**
- `src/sgp4/model_selector.h`

**Read:**
- EGM2008 published values (e.g. via `sgp4_references/hoots_roehrich_1980/hoots_roehrich_1980_spacetrack_report_3.md` §B or web source — Pavlis et al. 2012)
- `src/constants/constants_provider.h` (for GM uncertainty convention)
- `tests/test_sgp4/main.cpp` (model-selector tests)

**Audit card:**
- `design/audit/theoretical_basis_audit/model_selector.md`

---

## Primary issues

Three independent bugs in one file:

### 1. Preset names mislabel Halley as cubic/quartic

- Preset name `"householder"` is a misnomer for what is implemented as Halley (which is cubic, not quartic order)
- The `"eleven"` preset docstring further mislabels Halley as "quartic"

### 2. GM uncertainty inconsistency in WGS84 preset

- Docstring claims `±8e6 m³/s²` measurement uncertainty
- Code uses `±0.008 km³/s² = ±8e3 m³/s²`
- **Discrepancy: 10³× too small.**

NIMA TR 8350.2 published WGS84 GM = (3.986004418 ± 0.000000008) × 10¹⁴ m³/s² = (398600.4418 ± 0.0008) km³/s², so the docstring `±8e6 m³/s²` is wrong (decimal-shift error) — should be `±8e5 m³/s²`. Code should be `±8e-4 km³/s²`.

### 3. Unnormalized Jₙ values for n=5..9 disagree with EGM2008

The unnormalized Jₙ values for n=5, 6, 7, 8, 9 in the gravity-model factory disagree with the comment-derived (C̄ₙ₀ → Jₙ) values at the 4th-5th significant digit. The Jₙ values were either:
- Truncated incorrectly during transcription, or
- Computed from a pre-EGM2008 source (EGM96?) without source citation

Need re-derivation from full-precision EGM2008 published normalized C̄ₙ₀ values via `Jₙ = -C̄ₙ₀ · √(2n+1)`.

---

## Theory anchor

- EGM2008: Pavlis, Holmes, Kenyon, Factor (2012) "The development and evaluation of the Earth Gravitational Model 2008 (EGM2008)" JGR Solid Earth 117 B04406.
- WGS84 GM uncertainty: NIMA TR 8350.2 Table 3.4.1.
- Halley method: Halley (1694), Conte & de Boor (1980) §3.7 (cubic convergence).

---

## Fix scope

~30 minutes:

1. **Rename preset.** Either rename `"householder"` to `"halley"` (correct name) or document that the user-facing string is a stable API contract and only fix the docstring/comments.
2. **Fix GM σ.** Either docstring or code — make them agree. Recommend docstring fix to `±8e5 m³/s²` (matching NIMA TR 8350.2) and code fix to `±8e-4 km³/s²`.
3. **Recompute Jₙ for n=5..9.** From full-precision EGM2008 C̄ₙ₀ values, compute unnormalized Jₙ to ~12 significant digits. Update the table.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_sgp4.exe`:
   - All near-earth cases PASS (regression — J₅..J₉ don't affect SGP4 near-earth, only deep-space sensitivity)
3. Add a regression test (in `tests/test_sgp4` or `tests/test_perturbation`):
   - Verify `J5..J9` values match EGM2008 to 10 sig figs
   - Verify preset name returns expected solver type
4. Manual: print `J5..J9` after model selection; cross-check against EGM2008 reference values.

---

## References

- Audit card: `design/audit/theoretical_basis_audit/model_selector.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md`

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
