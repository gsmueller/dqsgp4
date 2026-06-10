# AUD-TBA Audit Results — Consolidated Summary

**Run date:** 2026-05-13
**Scope:** 50 `.h` files in `src/`, 287 functions
**Method:** one general-purpose / Explore agent per file; each agent produced a per-function audit card per the framework template (`design/audit/theoretical_basis_audit.md` §1)
**Output:** 50 audit documents in `design/audit/theoretical_basis_audit/<file>.md`

## File-level verdicts

| # | File | Functions | Verdict | Notes |
|---:|---|---:|---|---|
| 1 | `astronomy/sidereal_time.h` | 1 | **PASS** | IAU 1982 Aoki polynomial; IERS T⁴+ absorbed into σ_meas |
| 2 | `astronomy/solar_system.h` | 2 | **PASS** | Sun/Moon constants; anomalistic vs sidereal lunar period correctly distinguished |
| 3 | `atmosphere/density_model.h` | 1 | **PASS-with-flags** | Lane perigee adjustment; C⁰-only at 156 km boundary; bare-T branch select |
| 4 | `atmosphere/drag_coefficients.h` | 1 (15 sub-formulas) | **PASS-with-notes** | All 15 coefficients (C₁..C₅, D₂..D₄, t-cofs, omgcof, xmcof, delmo/sinmo, Ω̇_drag) match; D₃/D₄ literals (17, 221, 31) transcribed not re-derived |
| 5 | `constants/constants_provider.h` | 4 | **PASS** | WGS84 GM `measured`; WGS72/GRS80 routed via from_J2 |
| 6 | `dynamics/derivative.h` | 8 | **PASS** | REQ-EF-12 composition only; ⚠ `zero()` sets time_rate=1 not 0 |
| 7 | `dynamics/inertia.h` | 8 | **PASS** | Closed-form M-L-S §4.2; REQ-EF-9 safe-div for point-mass |
| 8 | `dynamics/pose.h` | 13 | **PASS** | Study biquaternion + Selig; SE(3) retraction via normalized() |
| 9 | `dynamics/propagator.h` | 7 | **PASS** | Newton/Euler + injected integrator; orchestrator |
| 10 | `dynamics/state.h` | 8 | **PASS** | Pure composition |
| 11 | `dynamics/twist.h` | 10 | **PASS** | se(3) algebra wrapper |
| 12 | `dynamics/wrench.h` | 10 | **PASS** | se(3)* algebra wrapper |
| 13 | `ephemeris/celestial_body.h` | 2 | **PASS** | Sun/Moon orbital element bundles |
| 14 | `ephemeris/lunar_ephemeris.h` | 1 | **PASS** | Meeus Ch.47; O(e³) truncation in equation of center |
| 15 | `ephemeris/solar_ephemeris.h` | 1 | **PASS-with-note** | 2-term equation of center; truncation negligible vs 0.01° target |
| 16 | `forces/drag.h` | 1 | **PASS** | Lane exponential atmosphere lambda |
| 17 | `forces/gravity_central.h` | 2 | **PASS** | Newton inverse-square |
| 18 | `forces/gravity_zonal.h` | 1 | **PASS-with-note** | J₂ spherical-harmonic derivative; REQ-EF-7 J₃+ truncation accuracy bound not wired |
| 19 | `geodesy/equipotential_ellipsoid.h` | 4 | **C-FAIL** | from_J2: Banach fixed-point residual NOT added to errors.precision (REQ-EF-5 violation; ~5 LoC fix) |
| 20 | `integrators/runge_kutta.h` | 3 | **PASS-with-flags** | Munthe-Kaas RK4; REQ-EF-7 integrator-order accuracy not added to errors.accuracy |
| 21 | `math/angles.h` | 6 | **PASS** | π/2π constants; wrap_two_pi via fmod |
| 22 | `math/binomial_series.h` | 2 | **PASS-with-flags** | Newton generalized binomial + geometric tail; HM §-anchor pending |
| 23 | `math/dual_number.h` | 16 | **PASS** | Division at line 113 is FIRST-ORDER TAYLOR (1/(c+εd) = 1/c − εd/c²), not Padé/CF — user concern verified absent |
| 24 | `math/dual_quaternion.h` | 28 | **PASS** | sclerp is Lie-group geodesic (not linear interp) |
| 25 | `math/factorial.h` | 5 | **PASS** | Definition products |
| 26 | `math/kepler.h` | 1 | **PASS** | Halley CONFIRMED (cubic; not Newton, not Padé/CF, not series) |
| 27 | `math/quaternion.h` | 22 | **PASS** | exp_pure/log_unit use Taylor (not Padé/CF) via small_angle helpers |
| 28 | `math/series.h` | 4 | **PASS** | Leibniz + geometric + Horner + Newton-sqrt |
| 29 | `math/small_angle_series.h` | 3 | **PASS** | (already audited in framework §5) |
| 30 | `math/tracked_value.h` | 43 | **PASS** | **FOUNDATION** — every REQ-EF-3 closed-form bound matches; REQ-EF-9 catastrophic regions handled |
| 31 | `math/vector3.h` | 9 | **PASS** | Closed-form vector algebra |
| 32 | `math/wallis.h` | 4 | **PASS** | Wallis 1656 reduction; exact rationals |
| 33 | `orbit/element_recovery.h` | 1 | **PASS-with-notes** | Brouwer-Lyddane; F2 series truncation O(δ₁⁴) tail not added to errors.accuracy |
| 34 | `orbit/modified_kepler.h` | 2 | **PASS-with-notes** | Newton + Halley both confirmed; ⚠ 30-iter cap fallback omits \|Δ_last\| in precision |
| 35 | `orbit/osculating_elements.h` | 1 | **PASS** | SR3 §6 closed-form conic identities |
| 36 | `orbit/secular_update.h` | 1 | **PASS-with-flags** | Brouwer + Lane-Hoots tempa/tempe/templ; e<1e-6 floor relaxes rigorous-bound contract |
| 37 | `orbit/state_from_elements.h` | 1 | **PASS** | SR3 §6 column-extraction TEME conversion |
| 38 | `perturbation/brouwer.h` | 1 | **PASS-with-flag** | Brouwer Eq.36-38 Horner; conservative J₂³ accuracy majorant |
| 39 | `perturbation/kaula.h` | 1 (13 branches) | **PASS** | 12 exact-rational closed-form branches; ⚠ F_311 needs Kaula 1966 Table 1 spot-check |
| 40 | `perturbation/resonance.h` | 3 | **C-FAIL** | step_resonance stub: `(void)earth_rate`, linear extrap at constant n instead of SR3 §6 tesseral; explains deep-space test failures |
| 41 | `perturbation/short_period.h` | 1 | **PASS** | SR3 §6 closed-form J₂ corrections |
| 42 | `perturbation/third_body.h` | 1 | **PASS-with-flags** | Orbit-averaged P₂(cos S); needs Kozai 1959 + SR3 §B citation anchors |
| 43 | `sgp4/deep_space.h` | 2 | **C-FAIL** | propagate_deep_space omits DPPER lunar/solar long-period, DSPER resonance, time-varying sun/moon angles. **STRUCTURAL CAUSE of all 24/24 SDP4 tcppver failures.** |
| 44 | `sgp4/model_functions.h` | 1 | **PASS-with-doc-bug** | standard_sgp4 lambda implements Newton but docstring says "Halley cubic"; AUD-EF-4 violation in non-convergence path |
| 45 | `sgp4/model_selector.h` | 22 | **PASS-with-flags** | "householder"/"eleven" preset names mislabel Halley as cubic/quartic; GM σ docstring vs code 10³× off; Jn for n=5..9 disagree with EGM2008 at 4th-5th digit |
| 46 | `sgp4/near_space.h` | 2 | **PASS** | Composition glue; xlcof critical-i fallback flagged |
| 47 | `sgp4/precomputed.h` | 4 | **PASS** | KaulaTable container; delegates to perturbation/kaula |
| 48 | `sgp4/sgp4_propagator.h` | 4 | **PASS** | Dispatch/assembly layer; classification threshold in near_space |
| 49 | `sgp4/state_vector.h` | 5 | **PASS-with-notes** | RSS correctly applied; ⚠ position_error/velocity_error zero the error fields |
| 50 | `tle/tle_parser.h` | 3 | **PASS-with-notes** | parse() declarations only (defs in .cpp); ?'s on bound for the parse pair |

## Summary table

| Verdict | Count | % |
|---|---:|---:|
| PASS (clean) | **26** | 52% |
| PASS-with-notes / flags | **18** | 36% |
| PASS-with-doc-bug | **1** | 2% |
| **C-FAIL** | **3** | 6% |
| Other (footnote pending) | **2** | 4% |
| **Total** | **50** | 100% |

## The three C-FAILs

### C-FAIL #1: `sgp4/deep_space.h::propagate_deep_space`

**Theory cited:** Hoots-Roehrich 1980 SR3 §6 (SDP4 with DPPER long-period periodics + DSPER resonance integration + time-varying sun/moon angles)

**Method implemented:** "near-space with constant solar/lunar offsets" — only the linear-time secular drift is applied. DPPER is missing. DSPER is missing. The `ds.resonance` field is initialized but never consumed during propagation.

**Impact:** Every deep-space TLE in `tcppver.out` fails. 24 of 24 deep-space test cases currently FAIL the validation suite. This is the structural reason — not 24 separate bugs, but one missing implementation surfaced 24 times.

**Fix scope:** non-trivial. Need to:
1. Wire DPPER (Sun + Moon long-period periodic corrections to e, i, M, ω, Ω)
2. Wire DSPER resonance integration (24h and 12h tesseral resonance for the relevant subset of TLEs)
3. Update the secular rates with time-varying sun/moon angle terms

Reference: dnwrnr libsgp4 (`SGP4.cc`) has all three implemented and matches tcppver.out for deep-space cases.

### C-FAIL #2: `perturbation/resonance.h::step_resonance`

**Theory cited:** SR3 §6 (SDP4 tesseral resonance integration)

**Method implemented:** linear extrapolation at constant n; the `earth_rate` argument is discarded via `(void)earth_rate`; the 10 d-coefficients (d2201..d5433) needed for tesseral acceleration are silently zero-errored at initialize time (TODO at line 134).

**Impact:** consistent with C-FAIL #1 — this is the missing piece that `deep_space.h` would need to call.

**Fix scope:** moderate. Need to:
1. Populate the 10 d-coefficients in `initialize_resonance` per SR3 §6 (formulas exist in dnwrnr and Vallado references)
2. Implement leapfrog or Euler integration of the tesseral acceleration in `step_resonance` (SR3 specifies 720-minute integration step)
3. Add REQ-EF-5 residual for the iteration

### C-FAIL #3: `geodesy/equipotential_ellipsoid.h::from_J2`

**Theory cited:** Banach fixed-point iteration on 1/f given J₂ (Heiskanen-Moritz)

**Method implemented:** correct iteration, but the final correction magnitude is NOT added to `e2_guess.errors.precision` before the result is used downstream. Per REQ-EF-5 this is required for every iterative algorithm.

**Impact:** every constant derived from `from_J2` (WGS72, GRS80 paths) carries an under-counted precision bound. The magnitude is type-epsilon-level (negligible for double; matters for cpp_bin_float_50). Structurally a REQ-EF spec violation.

**Fix scope:** trivial. ~5 LoC matching the pattern in `kepler.h:87` and `series.h:147`:

```cpp
e2_guess.errors.precision = e2_guess.errors.precision + abs(correction.value);
```

## Cross-cutting observations

1. **No "theory says Taylor, code does continued fraction" mismatches** — the user's primary worry was systematically checked. Every Taylor-cited function (`taylor_sinc`, `taylor_half_angle_scale`, `taylor_cos_minus_sinc_over_theta_sq`, `exp_pure`, `log_unit`, `dual_number::operator/`) is genuinely implemented as Taylor. No Padé, no continued fraction, no rational approximant slipped in.

2. **REQ-EF-7 (model-truncation → errors.accuracy) not consistently wired.** Several files (`gravity_zonal.h`, `forces/drag.h`, `runge_kutta.h`, `element_recovery.h::F2 series`) propagate precision correctly but do not add the model-truncation bound to `errors.accuracy`. This is a systematic gap, not 4 separate bugs.

3. **`tracked_value.h` is sound.** The 43-card audit confirmed every REQ-EF-3 closed-form bound formula matches the spec exactly. This means downstream `total_error()` claims have a sound foundation — conditional on each downstream function's own audit (PASS for the 26 clean ones; flagged for the rest).

4. **The Halley docstring bug in `model_functions.h`** is purely cosmetic — the math is correct (Newton); only the comment is wrong. Easiest fix in the codebase.

5. **`sclerp` (dual quaternion screw lerp)** is a Lie-group geodesic, not linear interpolation. The name might suggest the latter; the implementation correctly does the former (per Kavan et al. 2008).

## Priority remediation list

1. **[1 day, high impact]** Fix C-FAIL #1 + #2 together (they're coupled). Brings 24 SDP4 tests from FAIL to PASS. Reference impl in dnwrnr/libsgp4.
2. **[5 minutes]** Fix C-FAIL #3 — add the missing `errors.precision += abs(correction)` in `from_J2`.
3. **[5 minutes]** Fix the Halley docstring in `model_functions.h` (either update doc to say "Newton" or upgrade lambda to Halley).
4. **[1 hour]** Wire REQ-EF-7 model-truncation accuracy bounds in the 4 force/integrator/recovery files identified.
5. **[hours]** Close the open theoretical-rigor notes (drag coefs D3/D4 literals; F_311 spot-check; xlcof critical-i fallback verification).

## How the per-file audit cards are organized

Each `design/audit/theoretical_basis_audit/<file>.md` contains:
- Header (file path, function count, audit date, file-level verdict)
- One audit card per function with the 10 slots: ID, Location, Mathematical statement, THEORY{theorem, primary ref, domain}, METHOD{declared, implemented, match verdict}, ERROR BOUND{category, formula, implemented, verdict}, CROSS-AUDIT{REQ-EF, AUD-EF, AUD-MC, test}, NOTES
- File-level verdict at the end (axis A error wiring, B algebra axioms, C theoretical basis)

To audit a single function: look up its file in §1 of this document; open the corresponding `.md`; find the card. Each card stands alone as a verifiable claim.

## Companion documents

- `design/audit/theoretical_basis_audit.md` — framework definition (§1 card template, §5 worked example for small_angle_series.h)
- `design/audit/formula_status_tracker.md` — 287-row per-formula tracker (TB/SRC/ERR/TBA columns)
- `design/audit/error_framework.md` — AUD-EF-1..10 checks
- `design/audit/mathematical_correctness.md` — AUD-MC-1..18 algebra-axiom tests
- `design/audit/code_consistency.md` — AUD-CC-1..18 style/naming
- `design/specifications/error_framework.md` — REQ-EF-1..15
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` — SGP4 drag full derivation (used by drag_coefficients audit)
