# DQSGP4 — Remaining-Actions Execution Plan

**What this is.** The complete, ordered backlog of *every* remaining action to reach a genuinely-complete
DQSGP4. The issue register (`DQSGP4_ISSUE_REGISTER.md`) is the victory *checklist*; this is the execution
*order*. **"Continue" means: work straight down this plan — build, gate, validate, commit each item — and do
NOT stop after one.** I only pause for (a) a scope decision explicitly marked ⚠ below that you must make, or
(b) a failure I cannot resolve. Everything else is autonomous.

**Status at plan creation:** 8 items shipped this session (OR1, F3-a, D1, F2, E1, E2, G1, AD1); gate suite
**57/58**; the only open gated item is DS1. The plan below covers DS1 **plus all seven NEEDS-RULING items plus
final consolidation** — i.e. the "lot of actions" beyond the harness.

**Standing invariants (every single commit).** `test_sgp4` 33/33 623/623 · OR1 bit-exact · W11 AUD-CC ·
no NEW red in `run_acceptance`. Born-digital constants only (IAU NSFA / IERS / Wikipedia-cited primaries).
Validate each item with the full `-SkipVerifiers` suite (W11 + all ExeGates), then the full verifier sweep at
phase boundaries. Per-item: a named harness gate is the Definition-of-Done.

---

## ⚠ Scope decisions to set once (so I don't block mid-run)

These four items have a genuine "how far" that changes the work. My **recommended default** is in brackets;
tell me if you want different and I'll proceed on the default otherwise:

1. **D2 tesseral degree/order** — [cap at 4×4, EGM2008 born-digital, with a clear extension path to full degree].
2. **EPH ephemeris fidelity** — [Meeus-grade lunar/solar periodics with carried error terms + UT1−UTC; not full JPL DE].
3. **DRAG1 density depth** — [pluggable density-callback interface + the exponential model + a documented
   NRLMSISE-00 hook stub; not a full MSIS port].
4. **H1 attitude dynamics** — [BUILD it: `Matrix3<T>` + coupled-axis rigid-body attitude. It's large and
   orbital propagation doesn't need it, but it completes the dual-quaternion promise. Say "skip H1" to defer].

---

## Phase 0 — Tracked transcendentals (B3) · foundational ✅ DONE
Why first: DS1's DPPER and the ephemeris need tracked `atan2`/`asin`; doing them first keeps every downstream
conversion fully budget-propagating. **Landed: gate `B3` = `test_tracked_transcendentals` (55 checks), 33/33 +
OR1 + W11 green, no regression in the quaternion/precision modules that consume the new Taylor threshold.**
- **B3.1** `math/tracked_value.h` — add tracked `asin, acos, atan2, tan, sinh, cosh, tanh, expm1, log1p,
  log10, log2, hypot, min, max, clamp`, each with value + forward error bound + degenerate guards, in the
  style of the existing `exp/log/pow/cbrt`.
- **B3.2** `math/small_angle_series.h` — replace the hardcoded `1e-4` Taylor threshold with a `T`-dependent bound.
- **B3.3** `math/vector3.h` — add `normalize()` and scalar division.
- **Gate** `B3` = `test_tracked_transcendentals` — each function's precision tightens with wider `T`; bounds verified.

## Phase 1 — DS1: deep-space precision-tracking · the calling card ✅ DONE
Already scoped in register §DS1. Five sub-phases, **each bit-exact vs OR1 + 33/33 before the next**.
**Landed** (commits `0a35f40`, `57bc2ed`, `b04f0c1`, + gate): the whole deep-space evolution runs in
`TrackedValue<T>`; OR1 max drift 0; gate `DS1` = `test_sdp4_precision` (28 checks, 413× sharpening on a clean
case). Introduced `TrackedValue::model_coefficient` (finite-digit model floor → accuracy, binary storage →
precision-that-scales). Residual J₃/earth-model precision floor on high-ecc sats handed to CR1-b.
- **DS1.1** struct coefficient fields + `build_dpper_coefficients` (the dense 2-pass DSCOM loop) →
  `TrackedValue<T>`; DSCOM magic constants (ZNS, C1SS, ZES, ZNL, C1L, ZEL, …) → `from_truncated_decimal`.
- **DS1.2** `apply_deep_periodics` (DPPER) → `TrackedValue<T>` (uses tracked `atan2` from B3).
- **DS1.3** `propagate_deep_space` steps 1–9 (secular advance, sun/moon linear terms, drag scalings, J₃
  long-period) → `TrackedValue<T>`; **drop the epoch-budget re-injection** at the Kepler stage.
- **DS1.4** `perturbation::step_resonance` (the irez leapfrog integrator) → `TrackedValue<T>` — the hardest piece.
- **DS1.5** **Gate** `DS1` = `test_sdp4_precision` — a deep-space sat's precision budget tightens at
  `cpp_bin_float_50` vs `double`. Closes DS1 → **58/58**.

## Phase 2 — CR1-b: precision-honest constant sweep ✅ DONE
- Audited `defined(...)`/`measured(...)`. Fixed the real over-claim: six `model_selector.h` gravity-zonal
  `measured("Jₙ","σ")` whose σ was the decimal-ULP → `from_truncated_decimal` (value-preserving; 33/33 + OR1
  bit-exact). Genuine σ (EGM formal errors, resonance σ, geodesy GM/a⊕) kept as `measured`.
- **Gate** `CR1B` = `check_constant_honesty.ps1` (scans model files for σ == decimal-ULP) — PASS; all
  precision tests green. Categorization decision (model_coefficient in deep-space, from_truncated_decimal for
  static constants) documented in register §CR1-b for user review.

## Phase 3 — Force-model completeness
- **DRAG1** ✅ DONE — `forces/drag.h` pluggable `DensityModel<T>` callback + `exponential_density_model` +
  `nrlmsise00_density_model_stub` (documented MSIS seam) + generic `make_drag`; the flat 30% is gone, the
  density model's declared accuracy now propagates (model-derived). **Gate** `DRAG1` = `test_drag_density`.
- **EPH** ✅ DONE — `solar_ephemeris.h` folds the 2-term EoC `(13/12)e³` + distance `O(e²)` bounds into
  accuracy; `lunar_ephemeris.h` adds the `(5/4)e² sin2l` term + carries the ~2.5°/~0.9° omitted-periodic
  bounds (weakest-link fix); `sidereal_time.h` gains `compute_gmst_ut1` (ΔUT1=0 bit-exact). **Gate** `EPH` =
  `test_ephemeris` (11 checks, born-digital Earth-EoC ref).
- **D2** ✅ DONE — `constants/tesseral_harmonics.h` (C_nm/S_nm, EGM2008 C̄_22/S̄_22 born-digital + denorm) +
  `forces/gravity_tesseral.h` (Cunningham V_nm/W_nm Cartesian recursion, extensible; ECI↔ECEF GMST rotation).
  **Gate** `D2` = `test_tesseral` (9 checks: longitude-dependent sectoral accel, GMST coupling). General
  `F_lmp` (orbital-element Kaula path, distinct from this Cartesian gravity) left as-is, noted in register.

## Phase 4 — Injection + format completeness
- **INJ1** ✅ DONE — chose simplify (full wiring of inclination_function is OR1-impossible: the two resonance
  f220 forms are bit-different). Removed the dead `inclination_function` slot + the label-only `.drag()`
  selector; the remaining slots (secular_rates/kepler_solver/sidereal_time) are genuinely engine-dispatched.
  **Gate** `INJ1` = `test_injection` (5 checks: swapped kepler_solver/secular_rates change the result).
- **E3** ✅ DONE — `tle/omm_xml_parser.h` `parse_omm_xml` extracts the CCSDS OMM XML into the same `TleData`
  via the shared `populate_from_kv` (refactored from the KVN parser); boundary-aware element scan
  (MEAN_MOTION ≠ MEAN_MOTION_DOT), attribute-skipping. **Gate** `E3` = `test_omm_xml` (28 checks, XML==KVN).

## Phase 5 — Attitude dynamics (H1) ✅ DONE
- Built `math/matrix3.h` (Matrix3<T>: mat·vec, det, inverse, symmetric/diagonal factories) +
  `dynamics/attitude.h` (full Euler ω̇ = I⁻¹(τ − ω×Iω) with the gyroscopic coupling the diagonal path
  omitted, quaternion kinematics, RK4 step; point-mass torque restriction lifted). **Gate** `H1` =
  `test_attitude_dynamics` (11 checks: analytic free symmetric-top precession Ω=ω_z(C−A)/A, energy + |L|
  conserved, gyroscopic cross-feed, Matrix3 inverse on an off-diagonal tensor).

## Phase 6 — Final consolidation ✅ DONE
- **Completeness audit** ✅ — 3 parallel Explore agents swept sgp4/perturbation/orbit, forces/dynamics/
  integrators, and math/constants/geodesy/astro/ephemeris/tle. Findings: a `gravity_zonal.h` stale doc-lie
  (residual IS tracked — fixed), dead `ds.xlcof/ds.aycof` in `deep_space.h` (removed; retired a moot
  ratio(3,8) fallback discrepancy), a stale `gravity_central.h` "forthcoming" pointer (fixed). All else was
  already-tracked (C2-b R04, deferred pow, EPH carried bounds) or BY-DESIGN. Committed `b97374b`.
- **Full `run_acceptance` → `COMPLETE`** ✅ — **66/66 gates pass, VERDICT: COMPLETE** (85 Octave + 2 Python
  verifiers, all W-gates, OR1, 33/33 623/623, and B3/DS1/CR1B/DRAG1/EPH/D2/INJ1/E3/H1).
  **Then 67/67** with SC1 (series-based constants, Phase 7 below).
- **BY-DESIGN dispositions** ✅ — L1-L4 dispositioned ACCEPTED (inherent model limitations / physical facts;
  upgrade paths noted), register §BY-DESIGN; flagged for user re-open.
- **MEMORY consolidation** + **Final review** — in progress.

## Phase 7 — Series-based non-defined constants ✅ DONE (directive extension)
- **Directive (user, 2026-06-05, post-COMPLETE):** *"Constants that are not true by definition need to have
  their series-based precision and accuracy tracked."* Sharpens the panel's stamped-decimal encoding for the
  subclass of non-defined constants that are leading terms of a known series: generate from the series, carry
  PRECISION = representation (scales with `T`) + ACCURACY = series-truncation Σ_{k≥n}|c_k||t|^k (tightens with
  kept terms `n`) — not a typographic digit floor.
- **SC1 — obliquity exemplar** ✅ — `src/astronomy/obliquity.h` `obliquity_iau2006<T>(t, n_terms)` (IAU 2006
  ε_A(t), IERS Conventions 2010 Eq. 5.40, born-digital coefficients as `model_coefficient`); `model_selector.h`
  wires the modern preset's `fc.obliquity` to it at t=0 (value-preserving → C1/OR1 bit-exact). Gate
  `test_series_constants` (SC1): value matches IAU at t=0, secular drift + truncation accuracy at t≠0, accuracy
  tightens with `n_terms`, precision tightens ~20 orders double→`cpp_bin_float_50`. See register §SC1.
- **SC2 — all constants tracked; each J_k is a series truncation** ✅ (directive 2, user 2026-06-06). (a) J₂ and
  the SGP4-adopted GM moved `defined()`→`model_coefficient` in `model_selector.h` + `equipotential_ellipsoid.h`
  (value-preserving → OR1/33/33 bit-exact); datum geometry a/1f/ω stays `defined()`. The even J_{2n} are already
  generated from the ellipsoid series (`EquipotentialEllipsoid::J2n`); `test_series_constants` §J proves the
  round-trip. `check_constant_honesty.ps1` now forbids any J_k/GM on `defined()`. (b) `solar_eccentricity`
  migrated to the VSOP secular series `astronomy::earth_eccentricity` (Meeus 25.4; born-digital; 0.016708634 at
  J2000, sharper than the old stamp; not on the SGP4 path so OR1 unaffected); C1 updated, §E gates it. Lunar
  mean eccentricity stays `measured` (periodic, not secular — the honest encoding). Sidereal-year series is the
  one remaining minor forward item. Register §SC1/SC2 details.

---

## Execution contract
When you say **continue**, I start at the top of the first unfinished phase and work down — each item built,
gated, OR1+33/33-validated, and committed on `session/2026-04-23` — surfacing only the ⚠ scope decisions above
(once, up front, with defaults) or a genuine blocker. I will not interpret "continue" as "do one thing." This
file is updated as items land so the plan always reflects true remaining work.
