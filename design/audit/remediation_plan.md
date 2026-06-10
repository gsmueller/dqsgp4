# AUD-TBA Remediation Plan — Index

**Source:** `design/audit/AUD_TBA_results.md` (50-file audit, 287 functions, 2026-05-13)
**Companion tracker:** `design/audit/formula_status_tracker.md` (TBA verdict directory)

This document is the **index** to the 14 independently executable remediation items produced by the AUD-TBA audit. Each item lives in its own file under `design/audit/remediation/` so that separate agents can pick up an item in a follow-on plan-mode session without needing to read the others.

## Index

| ID | Title | Tier | File |
|---:|---|---|---|
| R01 | SDP4 deep-space implementation | P0 | [`remediation/R01_sdp4_deep_space.md`](remediation/R01_sdp4_deep_space.md) |
| R02 | from_J2 Banach fixed-point residual | P0 | [`remediation/R02_from_j2_residual.md`](remediation/R02_from_j2_residual.md) |
| R03 | model_functions docstring + AUD-EF-4 | P1 | [`remediation/R03_model_functions_docstring.md`](remediation/R03_model_functions_docstring.md) |
| R04 | model_selector triple bug | P1 | [`remediation/R04_model_selector_triple_bug.md`](remediation/R04_model_selector_triple_bug.md) |
| R05 | gravity_zonal.h REQ-EF-7 J₃+ truncation | P2 | [`remediation/R05_gravity_zonal_req_ef_7.md`](remediation/R05_gravity_zonal_req_ef_7.md) |
| R06 | forces/drag.h REQ-EF-7 atmospheric-model truncation | P2 | [`remediation/R06_forces_drag_req_ef_7.md`](remediation/R06_forces_drag_req_ef_7.md) |
| R07 | runge_kutta.h integrator-order accuracy | P2 | [`remediation/R07_runge_kutta_req_ef_7.md`](remediation/R07_runge_kutta_req_ef_7.md) |
| R08 | element_recovery F2 series tail accuracy | P2 | [`remediation/R08_element_recovery_req_ef_7.md`](remediation/R08_element_recovery_req_ef_7.md) |
| R09 | Drag coefficient literal derivations | P3 | [`remediation/R09_drag_coefficient_derivations.md`](remediation/R09_drag_coefficient_derivations.md) |
| R10 | tracked_value representation_bound 0.5 ULP | P3 | [`remediation/R10_tracked_value_ulp.md`](remediation/R10_tracked_value_ulp.md) |
| R11 | taylor_half_angle_scale wide-T correction | P3 | [`remediation/R11_taylor_half_angle_wide_t.md`](remediation/R11_taylor_half_angle_wide_t.md) |
| R12 | Modified-Kepler + secular-update boundary handling | P3 | [`remediation/R12_kepler_boundary_handling.md`](remediation/R12_kepler_boundary_handling.md) |
| R13 | Composite-API hygiene sweep | P4 | [`remediation/R13_composite_api_hygiene.md`](remediation/R13_composite_api_hygiene.md) |
| R14 | Citation / spot-check / decl-vs-def completeness | P4 | [`remediation/R14_citation_completeness.md`](remediation/R14_citation_completeness.md) |

## Severity tiers

| Tier | Meaning | Count |
|---|---|---:|
| **P0** | Blocking C-FAIL (theory ↔ implementation mismatch or REQ-EF spec violation) | 2 |
| **P1** | High-value quick win (doc bug or small triple-fix) | 2 |
| **P2** | Systematic REQ-EF-7 wiring (model-truncation → `errors.accuracy`) | 4 |
| **P3** | Theory-rigor gap (transcribed-but-not-derived; conservative-but-rigorous bound) | 4 |
| **P4** | API hygiene / citation completeness | 2 |

## File-collision matrix (verification)

Each item touches a disjoint set of source files. No file appears in more than one item; parallel agents may execute in any order without conflict.

| Item | Files (write only) |
|---:|---|
| R01 | `deep_space.h`, `resonance.h` |
| R02 | `equipotential_ellipsoid.h` |
| R03 | `model_functions.h` |
| R04 | `model_selector.h` |
| R05 | `gravity_zonal.h` |
| R06 | `forces/drag.h` |
| R07 | `runge_kutta.h` |
| R08 | `element_recovery.h` |
| R09 | `drag_coefficients.h`, `sgp4_near_earth_drag_theoretical_basis.md` |
| R10 | `tracked_value.h` |
| R11 | `small_angle_series.h` |
| R12 | `modified_kepler.h`, `secular_update.h` |
| R13 | `state_vector.h`, `derivative.h`, `density_model.h` |
| R14 | `kaula.h`, `brouwer.h`, `tle_parser.h`, `solar_ephemeris.h`, `near_space.h` |

**No file appears in more than one item.** ✓

Theoretical anchors with read-only overlap (no write conflict):
- R01 and R04 both reference SR3 §B (different formulas)
- R05-R08 all touch REQ-EF-7 (each a different file with a different bound formula)
- R09 *writes* `sgp4_near_earth_drag_theoretical_basis.md`; R14 sub-item 5 *reads* it

## Recommended dispatch sequencing

| Batch | Items | Rationale |
|---:|---|---|
| 1 | R02, R03, R10, R11 | Trivial (≤5 LoC each); warm up the audit→code feedback loop |
| 2 | R05, R06, R07, R08 | All REQ-EF-7 wiring (same pattern, different files) |
| 3 | R04, R12, R13, R14 | Medium-sized cleanup items |
| 4 | R09 | Theory-rigor extension to derivation doc |
| 5 | R01 | SDP4 deep-space — biggest item, gates 24/24 SDP4 tests |

Two or three items from the same batch may be dispatched as parallel agents in a single message.

## Project-level verification (after all 14 land)

- ✓ All near-earth `tcppver.out` cases still PASS at machine precision (regression — currently 8/8)
- ✓ All 24 deep-space `tcppver.out` cases PASS (currently 0/24; R01 closes this)
- ✓ `test_propagator.exe` Phase 1-3 LEO smoke: all green (regression — currently 12/12)
- ✓ `test_math.exe`, `test_quaternion.exe`, `test_dual_quaternion.exe`: all green (regression)
- ✓ `total_error()` on a 1-day RK4 LEO state has a non-trivial `errors.accuracy` component (R05, R06, R07, R08 close this)
- ✓ All 287 rows in `formula_status_tracker.md` have non-`?` TBA values via the directory mapping

## Workflow for picking up an R-item

When an agent (or human) is ready to execute an R-item:

1. Open `design/audit/remediation/R##_<slug>.md` for the item.
2. The file is **self-contained**: it lists the source files to read/write, the audit cards that document the issue, the theory anchor, the fix scope, and the verification plan.
3. Either start work directly, or enter plan mode and use the file as the basis for a finer-grained plan.
4. When complete, update the item file's **Status** field from `OPEN` → `DONE` and update `formula_status_tracker.md`'s TBA directory accordingly.

## Out of scope for these 14 items

After all 14 land, the following matters from `AUD_TBA_results.md` are tracked separately:

- High-precision test at `T = cpp_bin_float_50` (REQ-SY-7 demonstration) — separate test-coverage effort
- Constants-swap test (wgs84 vs wgs72 vs grs80 round-trip — REQ-SY-9) — separate test-coverage effort
- Automated AUD-CC-* checks (`tests/audit/test_code_consistency.cpp`) — separate test-coverage effort
- Force-lambda mass convention ambiguity (linear slot stores acceleration, not actual force) — separate design decision

These can be tracked in a future `tests/coverage_plan.md` document or as a new R15..R18.

## History

- **2026-05-13** — Document created from approved plan `peppy-lobster`. All 14 items split into individual files under `remediation/`. All items OPEN.
