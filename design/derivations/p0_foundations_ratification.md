# P0 — Foundations ratification (math / constants / geodesy)

**The R5 ratification pass** (`design/PROFESSIONAL_LIBRARY_PLAN.md` §5/§6): the roadmap's P0 said the
foundations are "largely already professional — write/complete their theory notes, confirm each callable is
isolatable + oracle-gated; rebuild only where a theory note exposes a gap. Mostly ratification." This note IS
that ratification: the module-by-module audit of theory home, gate, and oracle, with the verdicts and the
honestly-flagged thin spots. No rebuilds were warranted; one provenance defect found en route was fixed (O2,
below).

## 1. The audit table

| module | theory home | gate(s) | oracle / basis | verdict |
|---|---|---|---|---|
| math/tracked_value.h | ch01_three_errors.md + the Constants Initiative (8 phases) | W12, B1, B2, B3, CR1, CR1B | per-function analytic bounds; honesty scans | ✅ ratified (the substrate; do-not-refactor per the reuse survey) |
| math/vector3, quaternion, dual_number, dual_quaternion | dual_quaternion_propagator.md + ch02/ch03 | I0.test_math, test_quaternion, test_dual_number, test_dual_quaternion | algebraic identities | ✅ ratified |
| math/matrix3.h | frame_chain.md (product), H1 register | H1, FRAME1/FRAME2 (R·Rᵀ=I, dets) | analytic identities + erfa (chain) | ✅ ratified |
| math/series.h | ch05_series.md | I0.test_math, geodesy consumers | alternating-series bound (analytic) | ✅ ratified (dead `horner` already removed, Stage 1) |
| math/tracked_polynomial.h | the corrected formula plan (Stage 1) | SC1 + consumers (obliquity/earth_orbit/FW/nutation args) | bit-identity vs prior sites; erfa downstream | ✅ ratified |
| math/bessel.h, orbit/kepler_series.h | ephemeris_series.md | EPH (dial-up, exact-Kepler cross-check, Laplace-limit) | exact Newton solve; convergence theory | ✅ ratified |
| math/spherical_harmonics.h | geopotential.md §3 | GEOPOT, D2 (bit) | summed-legacy + closed-form J₂ + monopole | ✅ ratified |
| math/angles, factorial, binomial_series, small_angle_series, wallis, model_value, kepler | ch04/ch05 + per-header derivations | I0.test_math, W16, orbit/SGP4 consumers | analytic identities | ✅ ratified — *thin spot flagged*: these carry header-level derivations rather than standalone notes; consumers' gates (SGP4 33/33, EPH) exercise them end-to-end. Acceptable for utility mathematics; a standalone note is warranted only if one becomes load-bearing in a new derivation. |
| constants/constants_provider.h | CONSTANTS_INITIATIVE_PLAN.md (8 phases, complete) | C1, C2, F3, CR1B | IAU/IERS/NGA born-digital values + provenance scans | ✅ ratified |
| constants/zonal_harmonics.h + sgp4/model_selector.h zonals | this pass (O2) | GAL1 (extended), D1 | **the in-repo EGM2008 file** | ✅ **fixed + ratified** (§2) |
| constants/tesseral_harmonics.h, gravity_field.h | geopotential.md | GEOPOT, D2, CR1B | EGM2008 coefficients verbatim; honesty tags | ✅ ratified |
| geodesy/equipotential_ellipsoid.h | ch14_equipotential_ellipsoid.md | I0.test_geodesy, test_wgs84, C2 | Somigliana/level-ellipsoid closed forms; published WGS84/GRS80 values | ✅ ratified |

**Summary verdict: the P0 prediction held.** The foundations were already theory-homed and gate-covered —
the Constants Initiative and the formula-layer stages had done P0's heavy lifting for constants/ and the
math series machinery before this pass. The one real defect found was provenance, not structure:

## 2. O2 — the EGM2008 Jₙ provenance defect (issue R04), found and fixed

The audit's one FAIL-grade finding, predicted by audit card #17 and issue R04, confirmed by measurement
(`tools/gen_egm_zonals.py`, which reads the in-repo born-digital NGA EGM2008 file
`datalib/EGM-08norm100.txt` — the source the old TODO wished for was in-repo all along):

- The stored modern-preset C̄ₙ₀ (n = 5, 7, 9) were transcription artifacts at the **5th significant figure**
  (n = 5: stored 0.0686729e-6 vs file 0.0686703e-6, rel 3.8e-5; n = 9: rel 1.9e-4).
- Worse, the stored **J₃ = −2.53215306e-6 never equaled −√7·C̄₃₀ at all** (the file-derived value is
  −2.53241052e-6, rel 1.0e-4) — the original arithmetic itself had been inconsistent, exactly the
  "comment reproduces the stored Jₙ only to 4–5 sig figs" smell the audit flagged.
- All J₅…J₉ were off by rel 0.6–3.0e-4.

**The fix** (modern presets only; the frozen `sgp4_standard`/WGS72 path is untouched):
- Both stores re-derived from the in-repo file at full 15-digit precision by the generator
  (`constants::ZonalHarmonics::egm2008` C̄ strings; the ModelSelector `wgs84_precise` Jₙ table, its comment
  rewritten and the R04 TODO closed).
- **The honest re-encode**: C̄ₙ₀ now carries the published formal-error grade (IERS Conventions 2010
  Table 6.2, ±0.49e-11-class) in the **measurement** channel via `measured` — the prior digit-floor
  `model_coefficient` convention was honest at 7 written digits but would *under-claim* at the file's 15
  (a 1e-21 floor against a genuine 5e-12 σ). CR1B-clean (σ ≠ decimal ULP).
- **The gate** (GAL1, 16 checks): the identity Jₙ = −√(2n+1)·C̄ₙ₀ re-checked against file values embedded
  in the test; **cross-site double-entry** — the generative store and the ModelSelector preset must agree
  to 1e-12 (the old state fails this at 1e-4); the measurement channel carries the σ grade, T-independent;
  precision still tightens with T.

Effect size: the boosted/modern-preset J₅–J₉ shift at the 4th–5th significant figure (a ~1e-4 relative
change on perturbations that are themselves ~1e-6 of gravity). No gate asserted the old absolute values
except GAL1's own expected constant (updated with the fix).

## 3. Dispositions carried forward

- The SR3-historical ephemeris family (`solar/lunar_ephemeris`, `celestial_body`): dispositioned in R4a
  (superseded-by-Meeus banners, retained with gates).
- `perturbation/third_body.h` (the orbit-averaged secular-rate form): remains unwired by design — it
  belongs to an averaged-element consumer; the Cartesian force (L4) serves the DQ propagator.
- `kaula::inclination_function`: still a flagged orphan (the resonance F_lmp home); wiring it is a
  documented future unification, not a P0 gap.
- O4 (sidereal ratio): closed in R4a as a generative ratification (NUT1).

With this pass, **P0 is complete** — every foundations module is theory-homed, gate-covered, and
oracle-anchored, and the one provenance defect the foundations carried (R04) is fixed and double-entry
gated. The replan's remaining item is R6 (NRLMSISE-00, the flagged opt-in).
