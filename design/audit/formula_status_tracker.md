# Formula Status Tracker

**Purpose.** Per-formula status tracker for every function in `src/*.h`. Tracks three states for each formula:

1. **TB** — *Theoretical basis*: is the underlying theorem/method identified?
2. **SRC** — *Source citation*: what document records the theoretical basis?
3. **ERR** — *Error integration*: does the implementation wire into the three-error framework (REQ-EF) properly?

Plus a derived state:

4. **TBA** — *Theoretical Basis Audit verdict*: has a per-formula audit card been written confirming theory ↔ method ↔ bound match?

Status values are explicit and human-updatable. **Do not pre-fill claims that have not been verified**: an empty `?` is honest; a fake `✓` is harmful.

---

## TBA Verdict Directory — POST-AUDIT (2026-05-13)

The AUD-TBA audit ran on 2026-05-13 and produced one audit document per `.h` file in `design/audit/theoretical_basis_audit/<file>.md`. Each per-row `TBA` cell in the tables below now inherits the **file-level verdict** from this directory. The `R##` column maps each file to its remediation item in `design/audit/remediation_plan.md`.

**Reading rule**: a row's `TBA` cell shows the file-level verdict from this directory unless the row's own audit card flags a per-function override. The few overrides are noted in the Notes column of the row.

| File | Verdict | R-item | Audit card |
|---|---|---|---|
| `astronomy/sidereal_time.h` | **PASS** | — | `sidereal_time.md` |
| `astronomy/solar_system.h` | **PASS** | — | `solar_system.md` |
| `atmosphere/density_model.h` | **PASS-with-flags** | R13 | `density_model.md` |
| `atmosphere/drag_coefficients.h` | **PASS-with-notes** | R09 | `drag_coefficients.md` |
| `constants/constants_provider.h` | **PASS** | — | `constants_provider.md` |
| `dynamics/derivative.h` | **PASS-with-note** | R13 | `derivative.md` |
| `dynamics/inertia.h` | **PASS** | — | `inertia.md` |
| `dynamics/pose.h` | **PASS** | — | `pose.md` |
| `dynamics/propagator.h` | **PASS** | — | `propagator.md` |
| `dynamics/state.h` | **PASS** | — | `state.md` |
| `dynamics/twist.h` | **PASS** | — | `twist.md` |
| `dynamics/wrench.h` | **PASS** | — | `wrench.md` |
| `ephemeris/celestial_body.h` | **PASS** | — | `celestial_body.md` |
| `ephemeris/lunar_ephemeris.h` | **PASS** | — | `lunar_ephemeris.md` |
| `ephemeris/solar_ephemeris.h` | **PASS-with-note** | R14 | `solar_ephemeris.md` |
| `forces/drag.h` | **PASS** | R06 | `drag.md` |
| `forces/gravity_central.h` | **PASS** | — | `gravity_central.md` |
| `forces/gravity_zonal.h` | **PASS-with-note** | R05 | `gravity_zonal.md` |
| `geodesy/equipotential_ellipsoid.h` | **C-FAIL** | **R02** | `equipotential_ellipsoid.md` |
| `integrators/runge_kutta.h` | **PASS-with-flags** | R07 | `runge_kutta.md` |
| `math/angles.h` | **PASS** | — | `angles.md` |
| `math/binomial_series.h` | **PASS-with-flags** | — | `binomial_series.md` |
| `math/dual_number.h` | **PASS** | — | `dual_number.md` |
| `math/dual_quaternion.h` | **PASS** | — | `dual_quaternion.md` |
| `math/factorial.h` | **PASS** | — | `factorial.md` |
| `math/kepler.h` | **PASS** | — | `kepler.md` |
| `math/quaternion.h` | **PASS** | — | `quaternion.md` |
| `math/series.h` | **PASS** | — | `series.md` |
| `math/small_angle_series.h` | **PASS-with-note** | R11 | `small_angle_series.md` |
| `math/tracked_value.h` | **PASS-with-note** | R10 | `tracked_value.md` |
| `math/vector3.h` | **PASS** | — | `vector3.md` |
| `math/wallis.h` | **PASS** | — | `wallis.md` |
| `orbit/element_recovery.h` | **PASS-with-notes** | R08 | `element_recovery.md` |
| `orbit/modified_kepler.h` | **PASS-with-notes** | R12 | `modified_kepler.md` |
| `orbit/osculating_elements.h` | **PASS** | — | `osculating_elements.md` |
| `orbit/secular_update.h` | **PASS-with-flags** | R12 | `secular_update.md` |
| `orbit/state_from_elements.h` | **PASS** | — | `state_from_elements.md` |
| `perturbation/brouwer.h` | **PASS-with-flag** | R14 | `brouwer.md` |
| `perturbation/kaula.h` | **PASS** | R14 | `kaula.md` |
| `perturbation/resonance.h` | **C-FAIL** | **R01** | `resonance.md` |
| `perturbation/short_period.h` | **PASS** | — | `short_period.md` |
| `perturbation/third_body.h` | **PASS-with-flags** | — | `third_body.md` |
| `sgp4/deep_space.h` | **C-FAIL** | **R01** | `deep_space.md` |
| `sgp4/model_functions.h` | **PASS-with-doc-bug** | R03 | `model_functions.md` |
| `sgp4/model_selector.h` | **PASS-with-flags** | R04 | `model_selector.md` |
| `sgp4/near_space.h` | **PASS** | R14 | `near_space.md` |
| `sgp4/precomputed.h` | **PASS** | — | `precomputed.md` |
| `sgp4/sgp4_propagator.h` | **PASS** | — | `sgp4_propagator.md` |
| `sgp4/state_vector.h` | **PASS-with-notes** | R13 | `state_vector.md` |
| `tle/tle_parser.h` | **PASS-with-notes** | R14 | `tle_parser.md` |

**Verdict aggregates:**
- **PASS** (clean): 26 files
- **PASS-with-notes / flags / note**: 20 files
- **PASS-with-doc-bug**: 1 file (model_functions.h)
- **C-FAIL**: 3 files (deep_space.h, resonance.h, equipotential_ellipsoid.h)

**R-items in scope:** R01..R14 (see `design/audit/remediation_plan.md`).

**Coverage:** every file is mapped. The `R-item` column is `—` when the file is fully clean PASS (no remediation needed). It is also `—` when a PASS-with-notes file's notes are minor enough to be tracked indirectly via another item's audit card cross-reference.

### R-item resolution status (post-2026-05-14)

| R-item | Status | Note |
|---:|---|---|
| R01 | **PARTIAL** | 14/24 deep-space PASS (was 0/24); 10 high-B*/Lyddane edges remain |
| R02 | **DONE** | commit 3cc618f |
| R03 | **DONE** | commit 3cc618f |
| R04 | **PARTIAL** | Halley naming + GM σ done; J5..J9 EGM2008 re-derivation TODO'd |
| R05 | **DONE** | commit 3cc618f |
| R06 | **DONE** | commit 3cc618f |
| R07 | **DONE** | commit 3cc618f |
| R08 | **DONE** | commit 3cc618f |
| R09 | **DEFERRED** | 625-line derivation extension didn't persist mid-dispatch; re-execute |
| R10 | **ALREADY-CORRECT** | file matches HEAD with 0.5 ULP; no change required |
| R11 | **DONE** | commit 3cc618f |
| R12 | **DONE** | commit 3cc618f |
| R13 | **DONE** | commit 3cc618f |
| R14 | **DONE** | commit 3cc618f + bonus xlcof critical-i bug fix |

**Project-level SGP4-VER validation (post-remediation):**

- **23/33 satellites PASS, 506/623 data points** (was 8/33, 146/623 baseline)
- All 8 near-earth at machine precision (~7e-9 km); sat 28350 pre-existing drag bug
- 15/24 deep-space PASS via R01's SDP4 implementation
- 10 remaining FAILs: 1 near-earth (28350) + 9 deep-space (drag-dominated or Lyddane edges)

---

## Status legend

| Symbol | Meaning |
|---|---|
| `✓` | confirmed by primary-source citation or written audit card |
| `?` | not yet verified — needs audit |
| `⚠` | partial — see Notes |
| `✗` | confirmed wrong / known mismatch |
| `n/a` | not applicable (passthrough, accessor, type construction with no formula) |

## Column semantics

| Column | What it tracks |
|---|---|
| **Line** | starting line of the function in its file |
| **Function** | bare function or operator name |
| **Computes** | brief description of the math (from the function inventory) |
| **TB** | is the theoretical basis identified and documented? |
| **SRC** | which document records the theory? (or "?" if none yet) |
| **ERR** | does the code add the correct bound to the correct error category? (REQ-EF check) |
| **TBA** | has a formal audit card been written? (writes are committed in `design/audit/theoretical_basis_audit.md` or per-file expansions) |
| **Notes** | open issues, dropped terms, related audits, known mismatches |

## Inventory summary — EXACT COUNT

**50 files, 287 functions total.** Verified by full-Read pass of every `.h` file (2026-05-13); per-file counts summed and arithmetic checked.

| Directory | Files | Functions |
|---|---:|---:|
| `astronomy/` | 2 | 3 |
| `atmosphere/` | 2 | 2 |
| `constants/` | 1 | 4 |
| `dynamics/` | 7 | 64 |
| `ephemeris/` | 3 | 4 |
| `forces/` | 3 | 4 |
| `geodesy/` | 1 | 4 |
| `integrators/` | 1 | 3 |
| `math/` | 12 | 143 |
| `orbit/` | 5 | 6 |
| `perturbation/` | 5 | 7 |
| `sgp4/` | 7 | 40 |
| `tle/` | 1 | 3 |
| **Total** | **50** | **287** |

**Per-file function counts** (verified):

| File | Functions |
|---|---:|
| `src/astronomy/sidereal_time.h` | 1 |
| `src/astronomy/solar_system.h` | 2 |
| `src/atmosphere/density_model.h` | 1 |
| `src/atmosphere/drag_coefficients.h` | 1 |
| `src/constants/constants_provider.h` | 4 |
| `src/dynamics/derivative.h` | 8 |
| `src/dynamics/inertia.h` | 8 |
| `src/dynamics/pose.h` | 13 |
| `src/dynamics/propagator.h` | 7 |
| `src/dynamics/state.h` | 8 |
| `src/dynamics/twist.h` | 10 |
| `src/dynamics/wrench.h` | 10 |
| `src/ephemeris/celestial_body.h` | 2 |
| `src/ephemeris/lunar_ephemeris.h` | 1 |
| `src/ephemeris/solar_ephemeris.h` | 1 |
| `src/forces/drag.h` | 1 |
| `src/forces/gravity_central.h` | 2 |
| `src/forces/gravity_zonal.h` | 1 |
| `src/geodesy/equipotential_ellipsoid.h` | 4 |
| `src/integrators/runge_kutta.h` | 3 |
| `src/math/angles.h` | 6 |
| `src/math/binomial_series.h` | 2 |
| `src/math/dual_number.h` | 16 |
| `src/math/dual_quaternion.h` | 28 |
| `src/math/factorial.h` | 5 |
| `src/math/kepler.h` | 1 |
| `src/math/quaternion.h` | 22 |
| `src/math/series.h` | 4 |
| `src/math/small_angle_series.h` | 3 |
| `src/math/tracked_value.h` | 43 |
| `src/math/vector3.h` | 9 |
| `src/math/wallis.h` | 4 |
| `src/orbit/element_recovery.h` | 1 |
| `src/orbit/modified_kepler.h` | 2 |
| `src/orbit/osculating_elements.h` | 1 |
| `src/orbit/secular_update.h` | 1 |
| `src/orbit/state_from_elements.h` | 1 |
| `src/perturbation/brouwer.h` | 1 |
| `src/perturbation/kaula.h` | 1 |
| `src/perturbation/resonance.h` | 3 |
| `src/perturbation/short_period.h` | 1 |
| `src/perturbation/third_body.h` | 1 |
| `src/sgp4/deep_space.h` | 2 |
| `src/sgp4/model_functions.h` | 1 |
| `src/sgp4/model_selector.h` | 22 |
| `src/sgp4/near_space.h` | 2 |
| `src/sgp4/precomputed.h` | 4 |
| `src/sgp4/sgp4_propagator.h` | 4 |
| `src/sgp4/state_vector.h` | 5 |
| `src/tle/tle_parser.h` | 3 |
| **TOTAL** | **287** |

**Audit progress (4 of 287 = 1.4%):**

- **3 functions PASS** (`taylor_sinc`, `taylor_half_angle_scale`, `taylor_cos_minus_sinc_over_theta_sq`) per `design/audit/theoretical_basis_audit.md` §5
- **1 function stub-PASS** (`solve_kepler` per same document §6)
- **283 functions remaining** need audit cards

**Files missing from earlier draft of this tracker** (now added below in their respective sections):
- `src/astronomy/solar_system.h` (2 functions: `sgp4_standard()`, `compute()`)
- `src/perturbation/kaula.h` (1 function: `inclination_function`)
- `src/perturbation/resonance.h` (3 functions: `detect_resonance`, `initialize_resonance`, `step_resonance`)
- `src/perturbation/short_period.h` (1 function: `apply_short_period`)
- `src/perturbation/third_body.h` (1 function: `compute_third_body_rates`)

---

# Math Layer

## src/math/tracked_value.h

Core numeric type. Every error formula is the **bound formula** that downstream code inherits — auditing this file is **prerequisite** to every other audit.

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 28 | `ThreeErrors()` | default init | n/a | n/a | n/a | n/a | type ctor |
| 29 | `ThreeErrors(T,T,T)` | component init | n/a | n/a | n/a | n/a | type ctor |
| 33 | `total()` | triangle-inequality sum | ✓ | spec REQ-EF-2 | ✓ | ? | bound is m+p+a (correct upper bound; no independence assumed) |
| 37 | `rss()` | statistical RSS | ✓ | (statistics) | ✓ | ? | only valid if independent — not used for upper-bound paths |
| 48 | `operator+(ThreeErrors)` | per-category add | ✓ | REQ-EF-3 | ✓ | ? | add: bound = a + b, applied per category |
| 58 | `operator-(ThreeErrors)` | per-category sub | ⚠ | REQ-EF-3 | ⚠ | ? | clamps to ≥0; subtraction of errors is unusual — audit when this gets called |
| 69 | `operator*(T,ThreeErrors)` | scale errors | ✓ | REQ-EF-3 | ✓ | ? | uses \|s\|; correct for any sign of s |
| 75 | `operator*(ThreeErrors,T)` | scale errors | ✓ | REQ-EF-3 | ✓ | ? | symmetric to line 69 |
| 82 | `operator/(ThreeErrors,T)` | divide errors | ⚠ | REQ-EF-3 | ⚠ | ? | handles s=0 → max; check vs division formula in line 238 |
| 93 | `max_per_category` | per-cat max | ✓ | (utility) | ✓ | ? | used when result inherits the larger of two error sources |
| 104 | `apply<F>` | function on each | ⚠ | (utility) | ? | ? | accepts F as black box — caller must ensure F is non-negative-preserving |
| 118 | `TrackedValue()` | default init | n/a | n/a | n/a | n/a | type ctor |
| 120 | `TrackedValue(T, ThreeErrors)` | full init | n/a | n/a | n/a | n/a | type ctor |
| 123 | `TrackedValue(T,T,T,T)` | full init flat | n/a | n/a | n/a | n/a | type ctor |
| 129 | `exact_integer(int)` | zero error | ✓ | (def) | ✓ | ? | exact integers have zero error in all categories |
| 135 | `defined(str)` | from decimal str | ✓ | REQ-EF-1 | ✓ | ? | precision = `representation_bound(val)` |
| 143 | `defined_with_physical_uncertainty` | measured | ✓ | REQ-EF-1 | ✓ | ? | measurement = parsed sigma |
| 153 | `measured(val, sigma)` | measured | ✓ | REQ-EF-1 | ✓ | ? | same shape as 143 |
| 162 | `total_error()` | sum of categories | ✓ | REQ-EF-2 | ✓ | ? | passthrough to ThreeErrors::total |
| 164 | `reliable_digits()` | log10 ratio | ✓ | REQ-EF-11 | ✓ | ? | formula: ⌊−log10(total/|val|)⌋ |
| 176 | `operator+(TV,TV)` | add | ✓ | REQ-EF-3 | ✓ | ? | bound: a_bound + b_bound (exact for addition) |
| 185 | `operator-(TV,TV)` | sub | ✓ | REQ-EF-3 | ✓ | ? | same bound as add; cancellation handled by reliable_digits() |
| 192 | `operator-()` | unary neg | ✓ | REQ-EF-14 | ✓ | ? | preserves error exactly |
| 204 | `operator*(TV,TV)` | mul | ✓ | REQ-EF-3 | ✓ | ? | bound: \|a\|·b_err + \|b\|·a_err + a_err·b_err |
| 238 | `operator/(TV,TV)` | div | ✓ | REQ-EF-3 | ⚠ | ? | bound formula in REQ-EF-3; verify catastrophic-region handling (REQ-EF-9) |
| 280 | `sqrt(TV)` | square root | ✓ | REQ-EF-3 | ⚠ | ? | bound: err/(2·sqrt(val−err)); REQ-EF-9 when err≥val |
| 308 | `sin(TV)` | sine | ✓ | REQ-EF-3 | ✓ | ? | bound: \|cos(val)\|·err + err²/2, capped at 2 |
| 333 | `cos(TV)` | cosine | ✓ | REQ-EF-3 | ✓ | ? | bound: \|sin(val)\|·err + err²/2, capped at 2 |
| 349 | `atan(TV)` | arctan | ✓ | REQ-EF-3 | ✓ | ? | bound: err/(1+val²) |
| 365 | `atan2(TV,TV)` | two-arg arctan | ✓ | REQ-EF-3 | ⚠ | ? | bound: (\|x\|·y_err + \|y\|·x_err)/(x²+y²); REQ-EF-9 singular disc |
| 396 | `abs(TV)` | absolute value | ✓ | REQ-EF-14 | ✓ | ? | Lipschitz-1; bound passes through |
| 404 | `fmod(TV,TV)` | floating remainder | ? | ? | ? | ? | error bound on fmod is tricky — needs audit |
| 414 | `operator<` | bool result | n/a | n/a | n/a | n/a | passthrough |
| 415 | `operator>` | bool result | n/a | n/a | n/a | n/a | passthrough |
| 416 | `operator<=` | bool result | n/a | n/a | n/a | n/a | passthrough |
| 417 | `operator>=` | bool result | n/a | n/a | n/a | n/a | passthrough |
| 418 | `operator==` | bool result | n/a | n/a | n/a | n/a | passthrough |
| 419 | `operator!=` | bool result | n/a | n/a | n/a | n/a | passthrough |
| 423 | `TrackedValue(int)` | from int | ✓ | (def) | ✓ | ? | implicit conv from int = exact |
| 441 | `representation_bound(T)` | 0.5 ULP | ✓ | IEEE-754 | ✓ | ? | float-type rounding bound |
| 498 | `from_string(const char*)` | string→T | ⚠ | (decimal parse) | ⚠ | ? | private static; precision of string-to-T parse needs audit |
| 513 | `exact<T>(int)` | wrapper | ✓ | (utility) | ✓ | ? | thin wrapper around exact_integer |
| 528 | `ratio<T>(num,den)` | exact ratio | ✓ | (utility) | ✓ | ? | thin wrapper for exact rationals |

**File summary**: 36 entries; foundation file; **AUD-EF prerequisite**. The per-operation bound formulas at lines 176, 185, 204, 238, 280, 308, 333, 349, 365 are the **closed-form bounds that the entire rest of the codebase inherits** — they all need formal audit cards under AUD-TBA before downstream audits can claim soundness.

## src/math/series.h

Series evaluation with rigorous tail bounds (Leibniz / geometric).

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 27 | `alternating_series` | alt. series eval | ✓ | Leibniz thm | ✓ | ? | bound = magnitude of next term; cited in REQ-EF-4 |
| 77 | `geometric_series` | geom. series eval | ✓ | geom. convergence | ✓ | ? | bound = a_(N+1)/(1−r); REQ-EF-4 |
| 108 | `horner` | polynomial eval | ✓ | Horner's scheme | ✓ | ? | bound propagates via REQ-EF-3 per multiplication |
| 132 | `series_sqrt` | sqrt near 1 | ✓ | Newton iteration | ✓ | ? | REQ-EF-5 (iteration residual) |

## src/math/factorial.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 19 | `factorial(n)` | n! | ✓ | def | ✓ | ? | exact integer, no error |
| 30 | `double_factorial(n)` | n!! | ✓ | def | ✓ | ? | exact integer |
| 43 | `falling_factorial` | α(α−1)…(α−n+1) | ✓ | def | ⚠ | ? | for non-integer α: products propagate via REQ-EF-3 |
| 54 | `rising_factorial` | Pochhammer | ✓ | def | ⚠ | ? | same as falling, see Notes |
| 65 | `generalized_binomial` | Γ-form binom | ✓ | def (Newton) | ⚠ | ? | for non-integer α: needs audit of how the closed form is implemented |

## src/math/wallis.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 22 | `wallis_odd(n)` | ∫cos^n dφ (n odd) | ✓ | Wallis formula | ✓ | ? | closed-form rational |
| 29 | `wallis_even(n)` | ∫cos^n dφ (n even) | ✓ | Wallis formula | ✓ | ? | closed-form rational with π |
| 43 | `wallis(n)` | dispatch | ✓ | Wallis formula | ✓ | ? | calls odd/even based on n parity |
| 55 | `sin_power_cos_integral` | ∫sin^(2k)cos dφ | ✓ | direct integration | ✓ | ? | = 1/(2k+1) |

## src/math/binomial_series.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 48 | `make_binomial_evaluator` | (1+x)^α closure | ✓ | Newton gen. binom | ⚠ | ? | series convergence depends on |x|<1; bound formula needs audit |
| 114 | `geodetic_binomial_coefficient` | integrated binom | ? | ? | ? | ? | "geodetic integrated" — what theorem? needs source citation |

## src/math/small_angle_series.h

**AUDITED** — see `design/audit/theoretical_basis_audit.md` §5.

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 46 | `taylor_sinc` | sinc Taylor | ✓ | TBA §5.1 | ✓ | **PASS** | Leibniz bound \|θ\|⁶/5040 |
| 74 | `taylor_half_angle_scale` | arcsin(s)/s Taylor | ✓ | TBA §5.2 | ⚠ | **PASS-with-note** | non-alternating series; bound tight only for double; needs `1/(1−2s²)` correction for wide T |
| 106 | `taylor_cos_minus_sinc_over_theta_sq` | β(θ) Taylor | ✓ | TBA §5.3 | ✓ | **PASS** | Leibniz bound \|θ\|⁴/840 |

## src/math/kepler.h

**STUB-AUDITED** — see `design/audit/theoretical_basis_audit.md` §6.

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 32 | `solve_kepler` | Halley E−e·sinE=M | ✓ | Battin (1999) §5.3 | ✓ | **PASS-stub** | Kantorovich bound `|Δ_final|` |

## src/math/angles.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 14 | `pi<T>()` | π constant | ✓ | def | ✓ | ? | from boost::math constants; precision = representation_bound |
| 21 | `two_pi<T>()` | 2π | ✓ | def | ✓ | ? | same |
| 27 | `degrees_to_radians` | × π/180 | ✓ | def | ✓ | ? | propagates via REQ-EF-3 |
| 33 | `radians_to_degrees` | × 180/π | ✓ | def | ✓ | ? | same |
| 39 | `wrap_two_pi` | wrap to [0,2π) | ⚠ | def | ? | ? | uses fmod; tracked_value::fmod bound (line 404 tracked_value.h) is itself "?" |
| 51 | `wrap_neg_pos_pi` | wrap to [−π,π) | ⚠ | def | ? | ? | same; depends on fmod audit |

## src/math/vector3.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 29-30 | constructors | type init | n/a | n/a | n/a | n/a | |
| 34 | `operator+` | v+w | ✓ | REQ-EF-12 | ✓ | ? | per-component add |
| 39 | `operator-` | v−w | ✓ | REQ-EF-12 | ✓ | ? | per-component sub |
| 44 | `operator*(T,v)` | s·v | ✓ | REQ-EF-12 | ✓ | ? | per-component scale |
| 49 | `operator*(v,T)` | v·s | ✓ | REQ-EF-12 | ✓ | ? | symmetric |
| 54 | `dot(v,w)` | Σ v_i w_i | ✓ | def | ✓ | ? | propagates via REQ-EF-3 mul + add |
| 59 | `cross(v,w)` | det formula | ✓ | def | ✓ | ? | propagates via REQ-EF-3 |
| 68 | `magnitude()` | √(Σv²) | ✓ | Pythagorean | ✓ | ? | uses sqrt — REQ-EF-9 catastrophic-region inheritance |

## src/math/dual_number.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 44-50 | constructors | type init | n/a | n/a | n/a | n/a | |
| 56 | `zero()` | 0+ε·0 | ✓ | def | ✓ | ? | named ctor |
| 59 | `identity()` | 1+ε·0 | ✓ | def | ✓ | ? | named ctor |
| 64 | `epsilon()` | 0+ε·1 | ✓ | def | ✓ | ? | named ctor |
| 71 | `operator+` | per-cmp add | ✓ | def | ✓ | ? | REQ-EF-12 |
| 76 | `operator-` | per-cmp sub | ✓ | def | ✓ | ? | REQ-EF-12 |
| 81 | `operator-()` | unary neg | ✓ | def | ✓ | ? | REQ-EF-14 |
| 86 | `operator*` | (a+εb)(c+εd) = ac + ε(ad+bc) | ✓ | ε²=0 axiom | ✓ | ? | foundational rule; verifies AUD-MC-1 |
| 94 | `operator*(T,DN)` | scalar left | ✓ | def | ✓ | ? | per-cmp scale |
| 99 | `operator*(DN,T)` | scalar right | ✓ | def | ✓ | ? | symmetric |
| 113 | `operator/` | inverse via Taylor | ⚠ | def | ⚠ | ? | "Taylor expansion" comment — verify form & error bound |
| 125 | `sqrt(DN)` | sqrt derivative form | ✓ | forward-mode AD | ✓ | ? | AUD-MC-3 |
| 131 | `sin(DN)` | sin derivative form | ✓ | forward-mode AD | ✓ | ? | AUD-MC-3 |
| 136 | `cos(DN)` | cos derivative form | ✓ | forward-mode AD | ✓ | ? | AUD-MC-3 |

## src/math/quaternion.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 52-54 | constructors | type init | n/a | n/a | n/a | n/a | |
| 61 | `identity()` | 1+0i+0j+0k | ✓ | def | ✓ | ? | |
| 66 | `zero()` | 0 | ✓ | def | ✓ | ? | |
| 69 | `pure(v)` | from vector | ✓ | def | ✓ | ? | |
| 77 | `from_axis_angle` | exp(½θ·n̂) | ✓ | SU(2) parameterization | ✓ | ? | uses sin/cos of θ/2 |
| 88-91 | `scalar()`/`vector()` | accessors | n/a | n/a | n/a | n/a | |
| 96 | `operator+` | per-cmp add | ✓ | REQ-EF-12 | ✓ | ? | |
| 101 | `operator-` | per-cmp sub | ✓ | REQ-EF-12 | ✓ | ? | |
| 106 | `operator-()` | unary neg | ✓ | REQ-EF-14 | ✓ | ? | |
| 118 | `operator*(q1,q2)` | Hamilton product | ✓ | Hamilton (1843) | ✓ | ? | foundational; AUD-MC-4, AUD-MC-7 |
| 128 | `operator*(T,q)` | scalar left | ✓ | def | ✓ | ? | |
| 133 | `operator*(q,T)` | scalar right | ✓ | def | ✓ | ? | |
| 140 | `conjugate()` | (w, −v) | ✓ | def | ✓ | ? | AUD-MC-6, AUD-MC-7 |
| 145 | `magnitude_squared()` | q·q* | ✓ | def | ✓ | ? | |
| 150 | `magnitude()` | sqrt(\|q\|²) | ✓ | Pythagorean | ✓ | ? | uses sqrt — REQ-EF-9 |
| 157 | `inverse()` | q*/\|q\|² | ✓ | def | ✓ | ? | AUD-MC-9; needs guard on \|q\|=0 |
| 164 | `normalized()` | q/\|q\| | ✓ | def | ✓ | ? | REQ-EF-15 (retraction); REQ-EF-14 if already unit |
| 178 | `rotate(v)` | q·pure(v)·q* | ✓ | SU(2) double cover | ✓ | ? | AUD-MC-10, AUD-MC-11 |
| 199 | `exp_pure(v)` | cos(\|v\|) + sinc(\|v\|)·v | ✓ | exp on quaternion algebra | ✓ | ? | uses taylor_sinc when \|v\|<1e−4; PASS by composition |
| 218 | `log_unit()` | half-angle log | ✓ | exp inverse | ✓ | ? | uses taylor_half_angle_scale; AUD-MC-12, AUD-MC-13 |

## src/math/dual_quaternion.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 64-66 | constructors | type init | n/a | n/a | n/a | n/a | |
| 72 | `zero()` | (0, 0) | ✓ | def | ✓ | ? | |
| 75 | `identity()` | (1, 0) | ✓ | def | ✓ | ? | |
| 86 | `from_pose(q, t)` | (q, ½ t·q) | ✓ | Study biquaternion | ✓ | ? | |
| 97 | `from_screw(ω, v)` | pure DQ | ✓ | screw theory | ✓ | ? | |
| 106-122 | rotation/translation/angular/linear | accessors | n/a | n/a | n/a | n/a | |
| 127 | `operator+` | per-cmp add | ✓ | REQ-EF-12 | ✓ | ? | |
| 133 | `operator-` | per-cmp sub | ✓ | REQ-EF-12 | ✓ | ? | |
| 139 | `operator-()` | unary neg | ✓ | REQ-EF-14 | ✓ | ? | |
| 148 | `operator*(M̂,N̂)` | DQ product | ✓ | Study (1903) | ✓ | ? | AUD-MC-15 |
| 157-163 | `operator*(T,M̂)` etc | scalar mul | ✓ | def | ✓ | ? | |
| 173 | `conjugate()` | (q*, q'*) | ✓ | Study | ✓ | ? | AUD-MC-16 |
| 180 | `dual_conjugate()` | (q, −q') | ✓ | Study | ✓ | ? | |
| 188 | `combined_conjugate()` | (q*, −q'*) | ✓ | Study | ✓ | ? | |
| 197 | `magnitude_squared()` | M̂·M̂* | ✓ | def | ✓ | ? | |
| 211 | `inverse()` | M̂*/\|M̂\|² | ✓ | def | ✓ | ? | needs guard |
| 219 | `normalized()` | SE(3) retraction | ✓ | retraction map | ✓ | ? | REQ-EF-15 |
| 240 | `apply(p)` | rigid-body action | ✓ | SE(3) action | ✓ | ? | AUD-MC-17 |
| 252 | `apply_direction(d)` | rotation only | ✓ | SO(3) action | ✓ | ? | |
| 269 | `exp_screw(ω, v)` | screw exp | ✓ | Lie exp on se(3) | ✓ | ? | uses taylor helpers; AUD-MC-18 |
| 297 | `exp_screw(M̂)` | from pure DQ | ✓ | same | ✓ | ? | |
| 314 | `log_screw()` | screw log | ✓ | screw log inverse | ✓ | ? | AUD-MC-18 |
| 345 | `sclerp(M̂₁, M̂₂, t)` | screw linear interp | ✓ | Kavan et al. (2008) | ✓ | ? | |

---

# Dynamics Layer

## src/dynamics/state.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 44-46 | constructors | type init | n/a | n/a | n/a | n/a | |
| 52 | `identity_at_epoch()` | identity state | ✓ | def | ✓ | ? | named ctor |
| 59 | `from_kinematics(...)` | build from prim | ✓ | def | ✓ | ? | named ctor |
| 74-83 | accessors | passthrough | n/a | n/a | n/a | n/a | |

## src/dynamics/pose.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 41-45 | constructors | type init | n/a | n/a | n/a | n/a | |
| 50 | `identity()` | SE(3) identity | ✓ | def | ✓ | ? | |
| 56 | `from_rotation_translation` | (R, t) → DQ | ✓ | Study biquaternion | ✓ | ? | inherits from DQ::from_pose |
| 64-70 | accessors | passthrough | n/a | n/a | n/a | n/a | |
| 75 | `operator*` | SE(3) compose | ✓ | SE(3) group law | ✓ | ? | inherits DQ::operator* |
| 80 | `inverse()` | SE(3) inv | ✓ | SE(3) group law | ✓ | ? | inherits DQ::inverse |
| 84 | `normalized()` | retraction | ✓ | SE(3) retraction | ✓ | ? | REQ-EF-15 |
| 90 | `apply(p)` | rigid-body | ✓ | SE(3) action | ✓ | ? | |
| 96 | `apply_direction(d)` | rotation | ✓ | SO(3) action | ✓ | ? | |
| 105 | `sclerp(P₁, P₂, t)` | screw interp | ✓ | Kavan et al. | ✓ | ? | inherits |

## src/dynamics/twist.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 42-44 | constructors | type init | n/a | n/a | n/a | n/a | |
| 50 | `zero()` | zero twist | ✓ | def | ✓ | ? | |
| 53 | `pure_angular(ω)` | (ω, 0) | ✓ | def | ✓ | ? | |
| 58 | `pure_linear(v)` | (0, v) | ✓ | def | ✓ | ? | |
| 67 | `as_dual_quaternion()` | lift to pure DQ | ✓ | se(3) ↔ pure DQ | ✓ | ? | |
| 74-89 | linear ops | ✓ | REQ-EF-12 | ✓ | ? | per-component, no theory beyond linearity |

## src/dynamics/wrench.h

Mirror of twist.h — same audit pattern.

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 43-45 | constructors | type init | n/a | n/a | n/a | n/a | |
| 51 | `zero()` | zero wrench | ✓ | def | ✓ | ? | |
| 54 | `pure_torque(τ)` | (τ, 0) | ✓ | def | ✓ | ? | |
| 59 | `pure_force(F)` | (0, F) | ✓ | def | ✓ | ? | |
| 67 | `as_dual_quaternion()` | lift | ✓ | se(3)* ↔ pure DQ | ✓ | ? | wrench lives on dual of se(3) |
| 74-89 | linear ops | ✓ | REQ-EF-12 | ✓ | ? | per-component |

## src/dynamics/inertia.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 47-49 | constructors | type init | n/a | n/a | n/a | n/a | |
| 56 | `point_mass(m)` | I=0, diag=0 | ✓ | def | ✓ | ? | named ctor |
| 62 | `uniform_sphere(m,r)` | I=⅖mr² | ✓ | mech. textbook | ✓ | ? | named ctor |
| 72 | `diagonal(m, Ix, Iy, Iz)` | diag inertia | ✓ | def | ✓ | ? | named ctor |
| 86 | `momentum_of(twist)` | p = I·ξ | ✓ | rigid-body dyn | ✓ | ? | matrix-vector form for diag I |
| 114 | `acceleration_from_wrench(W)` | ξ̇ = I⁻¹·W | ✓ | Newton-Euler | ⚠ | ? | safe-div for point-mass; REQ-EF-9 |
| 144 | `linear_acceleration_from_force` | F/m | ✓ | F=ma | ✓ | ? | |

## src/dynamics/derivative.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 44-46 | constructors | type init | n/a | n/a | n/a | n/a | |
| 53 | `from_acceleration(twist)` | std time param | ✓ | def | ✓ | ? | |
| 58 | `zero()` | zero | ✓ | def | ✓ | ? | |
| 65-89 | linear ops | ✓ | REQ-EF-12 | ✓ | ? | per-component |

## src/dynamics/propagator.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 80 | ctor | type init | n/a | n/a | n/a | n/a | |
| 92-98 | accessors | passthrough | n/a | n/a | n/a | n/a | |
| 104 | `compute_acceleration(state)` | I⁻¹ Σ Wᵢ | ✓ | Newton-Euler | ✓ | ? | sums force lambdas, divides by I |
| 116 | `step(state, dt)` | integrator step | ✓ | (delegated) | ✓ | ? | dispatches to injected integrator |
| 130 | `propagate_to(state, t₀, tₑ)` | multi-step | ✓ | (delegated) | ✓ | ? | loop of step calls |

---

# Forces Layer

## src/forces/gravity_central.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 51 | `gravity_central(state, K)` | −GM·r̂/r² | ✓ | Newton (1687) | ✓ | ? | inverse-square law |
| 78 | `make_gravity_central(K)` | lambda closure | ✓ | (utility) | ✓ | ? | captures K |

## src/forces/gravity_zonal.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 62 | `gravity_J2(state, K)` | J₂ perturbation acc | ✓ | spherical harmonic | ⚠ | ? | accuracy bound (model truncation) needs to be added to errors.accuracy per REQ-EF-7 |

## src/forces/drag.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 77 | `make_drag_exponential(...)` | drag wrench lambda | ✓ | exp. atmosphere | ⚠ | ? | accuracy bound for exp. atmosphere model truncation per REQ-EF-7 |

---

# Integrators

## src/integrators/runge_kutta.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 55 | `lie_advance_pose(P, dt, ξ)` | SE(3) exp advance | ✓ | Munthe-Kaas (1999) | ✓ | ? | Lie-group integration; exp_screw inside |
| 69 | `euler(state, dt, accel)` | 1st-order Euler | ✓ | Euler (1768) | ✓ | ? | O(dt) local error → REQ-EF-7 |
| 88 | `runge_kutta_4(state, dt, accel)` | classical RK4 | ✓ | Runge-Kutta | ⚠ | ? | O(dt⁵) local error; theoretical bound for RK4 truncation needs to be added per REQ-EF-7 |

---

# Constants

## src/constants/constants_provider.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 39 | ctor | from ellipsoid | n/a | n/a | n/a | n/a | |
| 51 | `wgs84(tol)` | WGS84 const | ✓ | NIMA TR 8350.2 | ✓ | ? | uses EquipotentialEllipsoid::ctor |
| 73 | `wgs72(tol)` | WGS72 const | ✓ | SR3 | ✓ | ? | uses EquipotentialEllipsoid::from_J2 |
| 91 | `grs80(tol)` | GRS80 const | ✓ | Moritz (1980) | ✓ | ? | uses EquipotentialEllipsoid::from_J2 |

---

# Atmosphere

## src/atmosphere/density_model.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 57 | `compute_density_parameters(...)` | Lane perigee-adjust | ✓ | Lane (1965) | ⚠ | ? | derivation §2.3 of sgp4_near_earth_drag_theoretical_basis.md; C⁰-only continuity flagged A-D4 |

## src/atmosphere/drag_coefficients.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 107 | `compute_drag_coefficients(in)` | Lane-Hoots C₁..C₅, D₂..D₄, t-cofs, etc. | ✓ | sgp4_near_earth_drag_theoretical_basis.md §§5-15 | ⚠ | **partial-PASS** | per-coefficient theory done; model-truncation `errors.accuracy` not yet wired in per REQ-EF-7 |

---

# Orbit

## src/orbit/element_recovery.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 62 | `recover_mean_elements(...)` | Brouwer-Lyddane recovery | ✓ | Brouwer 1959 | ⚠ | ? | the iteration is documented in deprecated/013; status of error wiring needs audit |

## src/orbit/osculating_elements.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 48 | `compute_osculating(...)` | osculating from (E+ω) | ✓ | SR3 §6 | ⚠ | ? | derivation NEEDED — currently only cited |

## src/orbit/state_from_elements.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 40 | `elements_to_state(...)` | TEME r,v from elements | ✓ | SR3 §6 | ✓ | ? | classical Keplerian transform |

## src/orbit/modified_kepler.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 41 | `solve_kepler_newton(...)` | SGP4 mod Kepler (Newton) | ✓ | SR3 §6 | ✓ | ? | Newton; REQ-EF-5; verify same bound as math/kepler.h |
| 73 | `solve_kepler_halley(...)` | SGP4 mod Kepler (Halley) | ✓ | SR3 §6 | ✓ | ? | Halley; REQ-EF-5; same as math/kepler.h |

## src/orbit/secular_update.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 58 | `secular_advance(...)` | mean elements at t | ✓ | SR3 §6 + sgp4 derivation | ⚠ | ? | uses C₁ + D₂..D₄ + t-cofs; relies on drag_coefficients audit |

---

# Geodesy

## src/geodesy/equipotential_ellipsoid.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 89 | ctor (a, 1/f, GM, ω) | WGS84 path | ✓ | NIMA TR 8350.2 | ✓ | ? | direct from defining parameters |
| 214 | `from_J2(...)` | (a, J₂, GM, ω) | ✓ | GRS80 / WGS72 | ⚠ | ? | iteration on 1/f from J₂; REQ-EF-5 if convergence loop |
| 271 | `normal_gravity(lat)` | Somigliana formula | ✓ | Heiskanen-Moritz | ✓ | ? | |
| 282 | `Jn(2n)` | even zonal | ✓ | Moritz / WGS | ✓ | ? | from MacCullagh's formula |

---

# SGP4

## src/sgp4/state_vector.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 31 | ctor | type init | n/a | n/a | n/a | n/a | |
| 37 | `position_error()` | RSS pos err | ✓ | RSS | ✓ | ? | uses ThreeErrors::rss |
| 48 | `velocity_error()` | RSS vel err | ✓ | RSS | ✓ | ? | same |
| 60-80 | per-category accessors | passthrough | n/a | n/a | n/a | n/a | |

## src/sgp4/sgp4_propagator.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 104 | ctor | dispatches to NS / DS init | n/a | n/a | n/a | n/a | |
| 121 | `propagate(t)` | dispatch | ✓ | (delegate) | ✓ | ? | calls near_space::propagate or deep_space::propagate |
| 130-133 | accessors | passthrough | n/a | n/a | n/a | n/a | |

## src/sgp4/precomputed.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 61 | `set(l,m,p, val)` | cache | n/a | n/a | n/a | n/a | utility |
| 69 | `F(l,m,p)` | lookup | n/a | n/a | n/a | n/a | utility |
| 79 | `has(l,m,p)` | lookup | n/a | n/a | n/a | n/a | utility |
| 90 | `precompute(lat, ratio)` | Kaula inc fn | ✓ | Kaula (1966) | ⚠ | ? | inclination function — needs SRC for the closed forms |

## src/sgp4/model_functions.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 209 | `standard_sgp4()` | model fn bundle | ✓ | SR3 | ✓ | ? | wraps lambdas; trivial |

## src/sgp4/near_space.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 136 | `initialize_near_space(...)` | NS const setup | ✓ | SR3 §6 | ⚠ | ? | composes element recovery + drag coefs + secular rates + long-period; transitively audited via parts |
| 262 | `propagate_near_space(...)` | propagate | ✓ | SR3 §6 | ⚠ | ? | composes secular + long-period + Kepler + osculating + short-period |

## src/sgp4/deep_space.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 108 | `initialize_deep_space(...)` | DS const setup | ? | SR3 §6 | ? | ? | SDP4 init — Sun/Moon, resonance, secular; complex |
| 261 | `propagate_deep_space(...)` | DS propagate | ? | SR3 §6 | ? | ? | currently failing tcppver for all 24 deep-space cases |

## src/sgp4/model_selector.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 83 | `ZonalHarmonics::add` | cache | n/a | n/a | n/a | n/a | container util |
| 89 | `ZonalHarmonics::Jn` | lookup | n/a | n/a | n/a | n/a | |
| 97-100 | `J3`/`J5`/`J7`/`J9` | shortcut | n/a | n/a | n/a | n/a | |
| 103 | `max_degree` | container query | n/a | n/a | n/a | n/a | |
| 112 | `has` | container query | n/a | n/a | n/a | n/a | |
| 138 | `ModelConfiguration::Jn` | resolver | ✓ | (composite) | ✓ | ? | per-degree lookup with fallback |
| 167 | `ModelConfiguration::A(n,m)` | A_(n,m) | ✓ | Kaula | ✓ | ? | A_(n,m) = −J_n · a_E^n (when m=0) |
| 191-197 | CustomBuilder | builder pattern | n/a | n/a | n/a | n/a | |
| 510 | `select(preset, tol)` | preset config | n/a | n/a | n/a | n/a | utility |
| 625 | `custom()` | builder | n/a | n/a | n/a | n/a | utility |

---

# TLE

## src/tle/tle_parser.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 48-49 | `parse(line1, line2[, name])` | TLE format parse | ✓ | SR3 §3 | ✓ | ? | fixed-column format |
| 69 | `from_tle_data(td)` | to TrackedValue | ✓ | def | ✓ | ? | bstar / mean motion / inclination / eccentricity / etc. |

---

# Ephemeris

## src/ephemeris/celestial_body.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 58 | `make_sun(fc)` | Sun orbital elems | ✓ | Almanac | ✓ | ? | constants from FundamentalConstants |
| 84 | `make_moon(fc)` | Moon orbital elems | ✓ | Almanac | ✓ | ? | constants from FundamentalConstants |

## src/ephemeris/solar_ephemeris.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 41 | `compute_solar_position(body, t)` | ecliptic Sun pos | ⚠ | SR3 §B5 / almanac | ⚠ | ? | low-precision; truncation needs accuracy bound (REQ-EF-7) |

## src/ephemeris/lunar_ephemeris.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 46 | `compute_lunar_position(body, t)` | ecliptic Moon pos | ⚠ | SR3 §B5 / almanac | ⚠ | ? | low-precision; truncation needs accuracy bound (REQ-EF-7) |

---

# Perturbation

## src/perturbation/brouwer.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 83 | `compute_secular_rates(...)` | Brouwer J₂²+J₄ rates | ✓ | Brouwer (1959) | ✓ | ? | full derivation in BH61 cleanroom work (separate silo); ch10c, ch11d |

## src/perturbation/kaula.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 58 | `inclination_function(l,m,p, sin_i, cos_i)` | Kaula F_(l,m,p) | ✓ | Kaula (1966) ch.3 | ⚠ | ? | exact rational coefficients per AUD-CC; recursion structure needs audit |

## src/perturbation/resonance.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 63 | `detect_resonance(period_min)` | classify 24h/12h/other | ✓ | SR3 §6 | n/a | ? | inline; bool/enum result based on period threshold |
| 94 | `initialize_resonance(...)` | resonance state setup | ⚠ | SR3 §6 (SDP4) | ? | ? | for SDP4; currently failing all 24 deep-space tcppver tests — audit priority |
| 151 | `step_resonance(state, t, ω⊕)` | step the resonance | ⚠ | SR3 §6 (SDP4) | ? | ? | mutates state in place — REQ-EF-8 check |

## src/perturbation/short_period.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 70 | `apply_short_period(...)` | J₂ short-period corrections | ✓ | SR3 §6 | ⚠ | ? | derivation in deprecated/019; explicit closed form; many trig args propagate |

## src/perturbation/third_body.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 71 | `compute_third_body_rates(...)` | Sun/Moon secular rates | ✓ | SR3 §B (SDP4) | ⚠ | ? | used in deep-space; relies on ephemeris ⚠ |

---

# Astronomy

## src/astronomy/sidereal_time.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 49 | `compute_gmst(JD)` | GMST from Julian Date | ✓ | IAU 1982 / Vallado | ⚠ | ? | series formula; truncation needs accuracy bound |

## src/astronomy/solar_system.h

| Line | Function | Computes | TB | SRC | ERR | TBA | Notes |
|---:|---|---|:-:|---|:-:|:-:|---|
| 113 | `FundamentalConstants::sgp4_standard()` | SGP4 std astronomical const | ✓ | SR3 §B + Almanac | ✓ | ? | named ctor — ZNS, ZES, ZSINIS, ZCOSIS, ZNL, ZEL etc. |
| 227 | `DerivedOrbitalElements::compute(fc)` | derived orbital elements | ✓ | SR3 §B | ⚠ | ? | computes Sun/Moon secular params from FundamentalConstants |

---

# Statistics (exact)

| Layer | Files | Functions | TBA=PASS | Remaining audit |
|---|---:|---:|---:|---:|
| math         | 12 | 143 | 4 | 139 |
| dynamics     |  7 |  64 | 0 |  64 |
| forces       |  3 |   4 | 0 |   4 |
| integrators  |  1 |   3 | 0 |   3 |
| constants    |  1 |   4 | 0 |   4 |
| atmosphere   |  2 |   2 | 0 (partial drag derivation) | 2 |
| orbit        |  5 |   6 | 0 |   6 |
| geodesy      |  1 |   4 | 0 |   4 |
| sgp4         |  7 |  40 | 0 |  40 |
| tle          |  1 |   3 | 0 |   3 |
| ephemeris    |  3 |   4 | 0 |   4 |
| perturbation |  5 |   7 | 0 |   7 |
| astronomy    |  2 |   3 | 0 |   3 |
| **Total**    | **50** | **287** | **4** | **283** |

**Notes:**
- All 287 functions counted, including constructors, accessors, comparison operators, and container utilities.
- The 4 PASS functions are: `taylor_sinc`, `taylor_half_angle_scale`, `taylor_cos_minus_sinc_over_theta_sq` (all in `small_angle_series.h`), and `solve_kepler` (in `kepler.h`, stub PASS).
- TB/SRC/ERR status per-row is filled in across the file sections above; this table aggregates by directory.
- A row marked `n/a` for a category (e.g., default constructor with no math) is **not** counted as an audit candidate but **is** in the function total.

---

# Process for filling in the tracker

For each row with **TB=?** or **ERR=?** or **TBA=?**:

1. **Read the function in the source file**.
2. **Identify the underlying theorem**. Cite primary source (paper, textbook, definition).
3. **Update `TB` and `SRC`** columns.
4. **Trace the error-propagation path**. Confirm the bound formula matches REQ-EF-3..7 for the method used.
5. **Update `ERR`** column.
6. **Write a per-formula audit card** in `design/audit/theoretical_basis_audit.md` (or a per-file expansion under `design/audit/theoretical_basis_audit/<file>.md`).
7. **Update `TBA`** column with PASS / PASS-with-note / FAIL.

**Audit priority** (foundation up, with exact function counts):

1. `tracked_value.h` (**43 entries**) — every closed-form bound formula. Until this is PASS, every downstream bound is conditionally suspect.
2. `series.h` (4 entries) — Leibniz / geometric tail bounds.
3. `factorial.h` (5), `wallis.h` (4), `binomial_series.h` (2) — small math primitives = **11 entries**.
4. `angles.h` (6), `vector3.h` (9), `small_angle_series.h` (3 PASS already) — **15 entries**, of which 3 done.
5. `dual_number.h` (16), `quaternion.h` (22), `dual_quaternion.h` (28) — DQ library proper = **66 entries**.
6. `dynamics/*.h` (64 entries) — mostly inherit by composition.
7. `forces/*.h` (4), `integrators/*.h` (3) — model truncation accuracy bounds (REQ-EF-7) = **7 entries**.
8. `constants/` (4), `geodesy/` (4), `atmosphere/` (2), `tle/` (3), `ephemeris/` (4), `astronomy/` (3) — applied math = **20 entries**.
9. `orbit/` (6), `perturbation/` (7) — mid-level applied math = **13 entries**.
10. `sgp4/` (40 entries) — composite of everything above; partially covered by SGP4 derivation.

Running total: 43 + 4 + 11 + 15 + 66 + 64 + 7 + 20 + 13 + 40 = **283 remaining + 4 PASS = 287 total.** ✓

---

# Companion documents

- `design/specifications/error_framework.md` — REQ-EF-1..15 (what the budget must do)
- `design/audit/error_framework.md` — AUD-EF-1..10 (how to verify wiring)
- `design/audit/mathematical_correctness.md` — AUD-MC-1..18 (algebra-axiom tests)
- `design/audit/theoretical_basis_audit.md` — AUD-TBA framework + worked examples
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` — SGP4 drag full derivation
- `sgp4_references/hoots_roehrich_1980/hoots_roehrich_1980_math_derivation.md` — SR3 textbook-style
- BH61 cleanroom derivations in `sgp4_references/.../derivation/` (separate silo for Brouwer rates)

---

# Open structural questions

1. **Is the `?`-marked SRC for ERR=? entries actually a code defect**, or just an undocumented but correct bound?
2. ~~**Some files were missed by the first Explore agent**~~ — RESOLVED: second-pass enumeration covered all 50 files including `astronomy/solar_system.h` (2 functions) and `perturbation/{kaula,resonance,short_period,third_body}.h` (6 functions across 4 files). Total now exact at 287.
3. **TB=? and ERR=? entries overlap heavily** — same functions likely. Track jointly for efficiency.
4. **Two functions are concentrated in 1.4% of code but underpin everything**: `tracked_value.h`'s 43 entries (operator+, -, *, /, sqrt, sin, cos, atan, atan2 etc.) define the closed-form bound formulas for **every** other tracked computation in the codebase. An audit pass over `tracked_value.h` alone unlocks confidence in every downstream `total_error()` claim, even before per-file audits below it complete.

# How counts were verified

- **Method**: Glob enumeration of every `*.h` file in `src/` → 50 files.
- **Per-file**: an Explore agent Read every file in full and produced a numbered list of every function definition, then summed per-file counts and per-directory counts.
- **Cross-check**: I verified the per-file sum manually:
  - per-directory totals (3+2+4+64+4+4+4+3+143+6+7+40+3 = 287); the agent's grand total claimed 289 but its per-file sums add to 287 (a 2-off arithmetic slip in their summary, with "orbit: 7" instead of 6).
  - Authoritative number: **287 functions across 50 files** as of 2026-05-13.
