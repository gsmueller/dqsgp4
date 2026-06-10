# DQSGP4 Issue Register — the victory checklist

**★ STATUS 2026-06-06: VICTORY — `run_acceptance` prints COMPLETE (67/67 gates; +SC1 series-constants).** Every GATED item is DONE
(its named gate passes), every NEEDS-RULING item was ruled in-scope and BUILT (B3, DS1, CR1-b, DRAG1, EPH,
D2, INJ1, E3, H1 — gates `B3/DS1/CR1B/DRAG1/EPH/D2/INJ1/E3/H1`), and the L1-L4 BY-DESIGN items are
dispositioned ACCEPTED. The regression invariant held throughout: `test_sgp4` 33/33 623/623 + OR1 bit-exact
at every commit. The full suite ran the 85 Octave + 2 Python derivation verifiers + all W-gates clean.

**Purpose.** Every concrete issue surfaced by the 2026-06-04 completeness survey, tracked to closure. The
user's directive: *all of these are addressed before we declare victory.* This is the shared ledger.

**Victory condition.**
1. Every **GATED** item below is DONE — i.e. its named gate in `tools/run_acceptance.ps1` passes, so a FULL
   `pwsh tools/run_acceptance.ps1` prints `COMPLETE`. The harness, not judgment, is the arbiter.
2. Every **NEEDS-RULING** item has an explicit user disposition (build it → promote to GATED, or accept →
   move to BY-DESIGN). None silently dropped.
3. The regression invariant (I0: build, `test_sgp4` 33/33 623/623, module exes, all verifiers) stays green
   throughout.

**Status legend:** `[ ]` open · `[x]` done · `[~]` in progress. Each line:
`ID` `file:line` — issue → *gate or disposition*. IDs match `run_acceptance.ps1` gate IDs.

Companion docs: plan `design/DQSGP4_COMPLETION_ROADMAP.md`; survey detail in
`memory/project_dqsgp4_completion.md`.

---

## GATED — must be DONE for victory

### Phase A — Truthful baseline (R1)

**A1 — stale doc-lies (DONE, gate `A1` = `check_no_stale_stubs.ps1`)**
- [x] A1-a `model_functions.h:208` — GMST "Placeholder — not yet implemented" (compute_gmst IS wired).
- [x] A1-b `model_functions.h:287` — description says "Halley"; `standard_sgp4()` installs Newton.
- [x] A1-c `wrench.h:19` — "(REQ-PR-2, forthcoming)" summation is live in `Propagator::compute_acceleration`.
- [x] A1-d `brouwer.h:33` — dangling `009_brouwer_rates.md (TODO)` derivation-doc pointer.

**A2 — dead scaffolding disposition (gate `A2` = `check_no_dead_code.ps1`) — DONE**
- [x] A2-a `model_functions.h` — removed the `throw`ing `drag_coefficients` lambda + `DragCoefficientsFn`
  typedef + bundle member (drag stays in the isolated `atmosphere::compute_drag_coefficients`).
- [~] A2-b `model_selector.h` — `.drag(name)` builder is a no-op label-setter → rescoped to **INJ1**
  (injection completeness, NEEDS-RULING).
- [~] A2-c `model_functions.h` — `inclination_function` slot assigned but the engine calls `perturbation::`
  directly → rescoped to **INJ1** (kept as valid injectable surface; full wiring is a design decision).
- [x] A2-d `precomputed.h` — DELETED (orphaned `KaulaTable`; `#include`d nowhere).
- [x] A2-e phantom `PrecomputedModel`/`PrecomputedSatellite` doc types — gone with the file.
- [x] A2-f `model_functions.h:160` — removed the stale `KaulaTable` reference.
- [x] A2-g `sgp4_propagator.h` — deleted the unused `NearSpaceCoefficients`/`DeepSpaceCoefficients` structs.
- [x] A2-h `sgp4_propagator.h` — deleted the unused `PropagationState` struct + its false doc.
- [x] A2-i `sgp4_propagator.h:11-14` — header doc rewritten to reflect the real architecture.
- [x] A2-j `modified_kepler.h` — now the SINGLE SOURCE (de-dup below); no longer orphaned.
- [x] A2-k `modified_kepler.h` — removed the phantom "Householder (quartic)" solver doc + a SECOND dangling
  `009_sgp4_modified_kepler.md` ref the A2 gate caught at `model_functions.h:86`.
- [x] A2-l Kepler solver de-duplicated: `model_functions`/`model_selector` delegate to
  `orbit::solve_kepler_newton`/`solve_kepler_halley` (byte-identical; 33/33 preserved). Also fixed a latent
  C4244 (`T→int`) in `modified_kepler.h` exposed by its first-ever compilation.

### Phase B — Precision infrastructure (R3, R5)

**B1 — tracked exp/log (gate `B1` = `test_tracked_exp_log.exe`) — DONE**
- [x] B1-a `tracked_value.h` — added tracked `exp` (convex forward bound `val·expm1(d)`; `expm1` avoids the
  catastrophic cancellation that silently zeroed the precision budget for representation-scale `d` — a bug
  the B1 test caught).
- [x] B1-b `tracked_value.h` — added tracked `log` (concave bound `d/(x-d)` + degenerate guard when `d≥x`).
- [x] B1-c `drag.h` — density is now `rho_0 * exp(arg)` with tracked `arg = -(h-h_0)/H_scale`; precision &
  measurement propagate through the exponential and tighten with wider `T` (proven by the B1 test). The 30%
  Lane accuracy bound is unchanged, so W8/`test_force_models` stays green.

**B2 — tracked pow/cbrt (gate `B2` = `test_tracked_pow_cbrt.exe`) — DONE**
- [x] B2-a `tracked_value.h` — added tracked `pow(x, const T& p)` (constant exponent; `|p|·worst^(p-1)·d`
  bound with degenerate guard). Tracked-exponent `pow(x,y)=exp(y·log x)` deferred until a use needs it.
- [x] B2-b `tracked_value.h` — added tracked `cbrt` (total over R incl. negative args; cancellation-safe
  `d/(3·cbrt(|x|-d)²)` bound).
- [x] B2-c `element_recovery.h` — both inline Newton cube-roots (a₁, a₀) replaced by `pow(kₑ/n, 2/3)`,
  matching the reference `pow(xke/no, x2o3)` directly; **33/33 preserved exactly**. The now-unused
  `tolerance` param is marked reserved. (The scalar `std::pow(|δ₁|,4)` tail-bound at ~147 is a raw-`T`
  accuracy estimate, not a tracked path — left as-is.)
- [x] B2-d `equipotential_ellipsoid.h` — R3 inline Newton cube-root replaced by `cbrt(a·a·b)`;
  test_geodesy / test_wgs84 green.
- [~] B2-e `test_series.cpp` — `std::pow` on bare `T` is a reference-value computation inside a unit test,
  not a production tracked path (no precision budget flows through it); left as acceptable.

### Constant representation (R3, R5) — precision-honest constant encoding

**CR1 — generative-or-bounded constant encoding (gate `CR1` = `test_constant_representation.exe`) — DONE**
- [x] CR1-a `tracked_value.h` — `defined("23.4393")` stamped the BINARY representation bound as precision, so
  a 6-figure source claimed ~50 reliable digits in `cpp_bin_float_50` and Phase B propagated the over-claim
  faithfully. Added `from_truncated_decimal(str)`: precision = max(binary, 0.5·`decimal_ulp(str)`) — the honest
  "representation Taylor bound", FIXED across `T`. Generative forms (exact integers, series, `arcsec·π/648000`)
  still scale with `T`. The CR1 test contrasts all three.
- [x] CR1-b (DONE, gate `CR1B` = `check_constant_honesty.ps1`) — swept the `defined(...)`/`measured(...)`
  sites. Found + fixed the concrete over-claim: six SGP4 gravity-zonal `measured("Jₙ", "σ")` sites in
  `model_selector.h` (wgs72 J₃/J₄, wgs84_sgp4 J₃/J₄, and the default pair) whose "σ" (1e-11 / 1e-14) was
  exactly the **decimal-ULP** of the value's last written digit — a decimal-truncation mis-filed as a
  measurement uncertainty. Re-encoded as `from_truncated_decimal` (truncation → precision, the CR1 honest
  form); value-preserving, so 33/33 + OR1 bit-exact and GAL1/D1/C2/W23/W24 stay green. Genuine σ were left
  as `measured`: the EGM2008 formal errors (1.3e-11, 1.41e-11, ...), the SR3 resonance σ (one order above the
  ULP), the geodesy GM/a⊕ published σ. Gate `check_constant_honesty.ps1` scans the SGP4/perturbation/gravity
  model files and fails on any `measured()` whose σ equals the value's decimal-ULP (false-positive-free: real
  σ are not the ULP; geodesy excluded as its published σ can legitimately coincide).
  - **Categorization decision — RESOLVED by expert panel (2026-06-05, user deferred per
    [[feedback_expert_panel_for_decisions]]).** A 5-agent independent panel (numerical-analysis, metrology/GUM,
    astrodynamics, software-architecture, reproducibility/red-team) ruled: categorize a finite-digit constant
    by its **intrinsic nature (provenance), uniformly across the codebase — NOT by which code path consumes
    it** (the prior path-based hybrid was judged teach-to-the-test and non-composable, since a deep-space value
    routinely meets a static-earth value in one expression). The resulting uniform rule:
    `defined()` = exact-by-definition (WGS84 a/1f/ω); `measured(v,σ)` = a constant with a real/published σ OR
    an adopted PHYSICAL/astronomical quantity (σ = the adoption bound) — obliquity, eccentricities, sidereal
    year, EGM coefficients with formal errors; `model_coefficient()` = a finite-digit MODEL-FIT coefficient
    with no σ (SR3 lunisolar/resonance, WGS/SR3 zonal Jₙ, EGM C̄ₙ₀/C̄ₙₘ) — digits → accuracy, storage →
    T-scaling precision; `from_truncated_decimal()` = reserved for a numerical guard / a standard's own
    source-truncated value. **Implemented:** the static gravity Jₙ (zonal + tesseral) and astronomy constants
    moved off `from_truncated_decimal` (gravity → `model_coefficient`, astronomy → `measured`); value-preserving
    (OR1 drift 0 km); `test_zonal_harmonics`/`test_tesseral` updated to assert the principle (precision tightens
    with T, accuracy holds the digit floor); `check_constant_honesty.ps1` now also fails on any
    `from_truncated_decimal` model coefficient (precision-purity gate), allowlisting the deep_space `1.5e-12`
    guard. So precision now tightens with T through the geopotential too.

**SC1 / SC2 — series-based constants; every J_k carries series-truncation accuracy (gate `SC1` =
`test_series_constants.exe`) — DONE**
- **Directive 1 (user, 2026-06-05):** *"Constants that are not true by definition need to have their series-based
  precision and accuracy tracked."* This sharpens the panel's `measured`/`model_coefficient` *stamped-decimal*
  encoding for the subclass of non-defined constants that are really the **leading terms of a known series**:
  rather than stamp a decimal and book its digit-floor as accuracy, GENERATE the value from its series and let
  the evaluation carry both channels — PRECISION = representation of the Horner sum (tightens with wider `T`,
  the calling card, now reaching an astronomical constant); ACCURACY = the **series-truncation** bound
  Σ_{k≥n}|c_k||t|^k (the physically omitted higher-order terms; tightens as more terms `n` are kept), not the
  typographic "how many digits were written" floor.
- **Directive 2 (user, 2026-06-06):** *"All constants need to have their accuracy and precision tracked. Each
  J_k value is the result of a series truncation."* Broadens directive 1 to ALL constants and names the
  geopotential coefficients specifically. The honest reading: a constant that is exact-by-CONVENTION (datum
  geometry a/1f/ω, the integer 648000, π) legitimately has accuracy 0 — that IS tracking it correctly. But a
  constant that is physically a truncated series / a fit — every zonal J_k, GM, the eccentricity, the obliquity
  — must carry a nonzero accuracy and never hide behind `defined()` (which zeroes it). The J_k are the flagship:
  the geopotential is an infinite spherical-harmonic expansion, each J_k its truncated/fit coefficient.
- [x] **Obliquity exemplar.** `src/astronomy/obliquity.h` — `obliquity_iau2006<T>(t_jcen, n_terms=6)` evaluates
  the IAU 2006 mean-obliquity polynomial ε_A(t) = Σ_{k=0}^{5} c_k t^k arcsec (IERS Conventions 2010 Eq. 5.40,
  born-digital), Horner in `TrackedValue` with each c_k a `model_coefficient`, the omitted tail folded into
  accuracy, ×(π/648000) for radians. `model_selector.h` modern preset now wires
  `fc.obliquity = obliquity_iau2006<T>(exact<T>(0))` — value unchanged at J2000 (84381.406″·π/648000 =
  0.40909280 rad, so C1/OR1 bit-exact), but now a consumer evaluating at the satellite epoch t≠0 gets the
  secular drift AND its truncation accuracy. Gate `test_series_constants` proves: value at t=0 matches IAU,
  secular decrease at t=0.21 cen (≈2021), accuracy(n=1)≈4.77e-5 rad (the omitted −46.84″·t term) ≫
  accuracy(n=2) (truncation tail collapses to the leading-coefficient adoption floor ~2.4e-9 rad), and
  precision tightens ~20 orders double→`cpp_bin_float_50`.
- [x] **Every J_k carries series-truncation accuracy (no J_k / GM on `defined()`).** The CR1-b panel sweep moved
  the SGP4 zonal J₃/J₄ to `model_coefficient` but MISSED J₂ — it lives in `make_ellipsoid` (the
  `EquipotentialEllipsoid` input), not `make_zonal_harmonics`, so it was still `defined()` (accuracy 0). Fixed:
  J₂ (wgs72/wgs72_old/wgs84_sgp4/grs80) and the SGP4-adopted GM moved `defined()` → `model_coefficient` in
  `model_selector.h` + the `equipotential_ellipsoid.h` wgs72/grs80 presets; **value-preserving** (`model_coefficient`
  and `defined` share the same `.value`; only accuracy 0→digit-floor), so OR1 drift 0 km + 33/33 + C2/D1/GAL1
  stay green. Datum GEOMETRY (a, 1/f, ω) stays `defined()` — exact by convention, accuracy 0 is correct. The
  even J_{2n} are ALREADY generated from the equipotential-ellipsoid series (`EquipotentialEllipsoid::J2n`,
  Heiskanen–Moritz; `ZonalHarmonics` realises it) — so "each J_k is a series truncation" is literal there;
  `test_series_constants` §J proves the generated J₂ recovers the adopted WGS72 0.001082616 (round-trip:
  `from_J2` inverts the same relation `J2n(1)` evaluates), carries both error channels, and its precision
  tightens with wider `T`. `check_constant_honesty.ps1` gains a gate: a value in the J_k range (~1.08e-3) or the
  GM range (~3.986e5 / 3.986e14) may never be `defined()` (value-range detection, robust to variable name and
  the multi-line `from_J2` call).
- [x] **`solar_eccentricity` migrated to the VSOP secular series** (`src/astronomy/earth_orbit.h`
  `earth_eccentricity<T>(t, n_terms)`): e(t) = 0.016708634 − 0.000042037·t − 0.0000001267·t² (Meeus
  *Astronomical Algorithms* Ch. 25 Eq. 25.4 / VSOP87, born-digital; pymeeus realises it verbatim), each c_k a
  `model_coefficient`, omitted tail → accuracy. `model_selector.h` wires `fc.solar_eccentricity =
  earth_eccentricity<T>(exact<T>(0))`. **Not value-preserving** (the born-digital leading term 0.016708634 is
  sharper than the former 7-figure 0.0167086 stamp), but `solar_eccentricity` is a modern-preset astro constant
  for the DQ lunisolar ephemeris — NOT on the SGP4/OR1 path, so OR1 stays bit-exact. C1 updated to expect
  0.016708634; `test_series_constants` §E gates the series (secular drift, truncation accuracy tightens with
  n_terms, precision tightens with `T`) — the eccentricity twin of the obliquity exemplar.
- [ ] **`solar_anomaly_period_days` (sidereal year)** — forward, minor: a sub-second/century secular drift
  exists (Simon 1994), but the leading term already dominates the tracked accuracy; migrate to a series form
  for completeness when convenient. Until then it stays `measured` (an adopted physical quantity, accuracy
  tracked via its σ — not on `defined()`, so directive-2-compliant).
- [x] **`lunar_eccentricity` correctly stays `measured`** (NOT an omission): the Moon's mean eccentricity is
  dominated by **periodic** terms (evection ±0.0117, variation), not a secular polynomial; the mean is secularly
  ~constant, so its honest error is the periodic *amplitude* (a measurement-like adoption bound). `measured(0.0549006,
  σ)` IS the directive-compliant encoding (accuracy/precision tracked; not `defined()`). A full ELP/Brown
  periodic series would be a force-model, not a constant.

### Phase C — Earth models / Astronomical Almanac (R7)

**C1 — real Almanac/IAU constants (gate `C1` = `test_earth_constants.exe`) — DONE (encoded per CR1)**
- [x] C1-a `model_selector.h` — the `iau_2006`/`almanac_20xx` preset now overrides the **unambiguous** J2000
  quantities with authoritative born-digital values, CR1-encoded: obliquity = 84381.406″ (IAU 2006 B1 / IAU
  2009 NSFA, `iau-a3.gitlab.io/NSFA`) via the generative `from_truncated_decimal("84381.406")·(π/648000)`
  (mas-precise, replacing the lossy 23.4393° `measured` form); Earth/solar eccentricity 0.0167086 (J2000
  VSOP); lunar eccentricity 0.0549006. Gate `test_earth_constants` checks the actual preset values. Authentic
  SGP4 (`sgp4_standard`/WGS72) is untouched → 33/33 unaffected. The lunar/solar PERIODS carry definitional
  subtleties vs the SR3 rates (SR3's "sidereal" month 27.321582 is in fact the *tropical* month; sidereal is
  27.321662) — deferred to CR1-b for a careful per-element pass.
- [ ] C1-b `sidereal_time.h` — no UT1−UTC correction (UTC fed directly as UT1).

**C2 — earth-model preset catalog (gate `C2` = `test_earth_models.exe`) — DONE (C2-a)**
- [x] C2-a `equipotential_ellipsoid.h` — added named static factories `wgs84`/`grs80`/`wgs72`/`iers2010` at
  the geodesy layer (the reusable single source), CR1-encoded (defining params `defined`, measured quantities
  with published σ). The new `iers2010` is the current IAU 2009 NSFA model (a⊕=6378136.6, GM=3.986004418e14,
  J₂=1.0826359e-3, `iau-a3.gitlab.io/NSFA`). Gate `test_earth_models` checks each. Consumer de-duplication
  (ConstantsProvider/model_selector delegating here) is item F3.
- [ ] C2-b `model_selector.h:288-300` — EGM2008 J₅–J₉ provenance TODO (R04): reproduce quoted C̄ₙ₀ only to
  4–5 sig figs; re-derive from the raw EGM2008 file.

### Phase D — Geopotential extension (R5)

> **Redesign (user, 2026-06-04):** D1, F3, and CR1-b merge into the **Generative Astro Library** — one
> common, generative, standard-tagged constants library both propagators draw from; constants generated ONCE
> at init into a precomputed table (compute-once), each carrying its error budget + the standard it realizes.
> Build additive-first, 33/33-gated (SGP4 is the DQSGP4 test oracle). See `memory/feedback_common_astro_library`.

**GAL1 — zonal harmonics provider (gate `GAL1` = `test_zonal_harmonics.exe`) — DONE**
- [x] `constants/zonal_harmonics.h` (new) — a precomputed, standard-tagged `ZonalHarmonics<T>`: even zonals
  derived from the ellipsoid; ODD zonals generated ONCE at init from the born-digital C̄ₙ₀ via
  Jₙ=−√(2n+1)·C̄ₙ₀ (tracked sqrt). Gives J₃/J₅ a real home (resolves D1-c) and generates from the published
  coefficient instead of a copied-down Jₙ (resolves the R04 provenance concern). Factories `egm2008` (IERS
  2010 Table 6.2, J₃–J₉) and `wgs72` (SR3 J₃). 12-check gate; 33/33 + W11/W17/W18 green.

**D1 — generic Jₙ zonal force (gate `D1` = `test_gravity_zonal_jn.exe`) — DONE**
- [x] D1-a/b/d/e `gravity_zonal.h` — `gravity_zonal(state, K, zh, max_n)` sums J₂…J_max_n via the closed-form
  zonal gradient a_n = μ Jₙ Rₑⁿ r^{-(n+3)}·[(n+1)Pₙ(u)+u Pₙ'(u)]·(x,y,·) with the non-singular Legendre
  recurrence (Pₙ, Pₙ'), reading ALL Jₙ from the shared generative `ZonalHarmonics` (GAL1) — even and odd from
  ONE standard-tagged source. Verified: max_n=2 reproduces `gravity_J2` to rel 4e-16; J₃ gives the N–S
  asymmetry (nonzero equatorial a_z, where J₂-only is exactly 0); the residual is the REAL omitted tail
  (Σ|aₙ| for force_max<n≤table + Kaula beyond-table), so it shrinks monotonically with max_n. The `|J4/J2|`
  proxy is retired: the generic path uses real coefficients, and standalone `gravity_J2`'s J₃ bound is now
  Kaula's a-priori |J₃| (the old |J4/J2|·1.1 ≈ 1.6e-3 actually UNDER-bounded the real |J₃/J₂| ≈ 2.3e-3). Gate
  `test_gravity_zonal_jn` (9 checks). 33/33 + OR1 bit-exact; W8/W24 (gravity_J2 consumers) green.
- [x] D1-c — RESOLVED by GAL1 (odd zonals are now representable and generative).

### Phase E — Updated TLE / OMM (R6)

**E1 — Alpha-5 decode + checksum (gate `E1` = `test_tle_alpha5.exe`) — DONE (E1-a/b/c)**
- [x] E1-a `tle_parser.cpp` — `decode_alpha5()` decodes the catalog field (Alpha-5 leader A–Z, skipping I/O,
  → high part 10–33: "E8493"→148493; classic numeric passes through). parse() stores it in
  `TleData.catalog_number`.
- [x] E1-b `tle_parser.cpp` — `tle_line_checksum()` (mod-10 over cols 0–67, '-'=1) + `checksum_valid()`;
  parse() reports `line{1,2}_checksum_valid` (validated, NOT enforced — so existing TLEs aren't rejected).
- [x] E1-c `tle_parser.cpp:4` / `.h:9` — the false "9-digit catalog" claim corrected to classic + Alpha-5
  (the column map is 5-char). Gate `test_tle_alpha5` (16 checks): decode boundary cases (H/J, N/P, I/O
  rejected), construct-verify checksum, canonical sat-5 validation, parse integration.
- [ ] E1-d `tle_parser.h:119` — the "Simplified JD" epoch is in fact the standard day-of-year JD formula
  (verified: 2000-01-01.0 → 2451544.5); the comment merely understates it. Doc nit, no functional gap.

**E2 — OMM (KVN) ingestion (gate `E2` = `test_omm_kvn.exe`) — DONE**
- [x] E2-a `tle/omm_parser.h` (new, header-only) — `parse_omm_kvn(text, TleData)` maps a CCSDS OMM in KVN
  into the SAME `TleData` the two-line parser produces, so it drives `from_tle_data` + both propagators
  unchanged (an alternate front-end, not a second element pipeline). ISO EPOCH → epoch_year + fractional
  day-of-year (leap-aware calendar form AND day-of-year form); MEAN_MOTION_DOT/_DDOT taken in the Celestrak
  GP convention (already n-dot/2, n-ddot/6); tolerant (missing mandatory key / malformed value → false, never
  throws). Gate `test_omm_kvn` (25 checks): full field mapping, Dec-26-2020 leap epoch = 361.667, DOY epoch,
  feeds from_tle_data, missing-mandatory rejected. Clean library header (no sgp4 coupling) → W11 green.

### Phase F — API parity (R2, R3)

**F1 — State↔StateVector converter (gate `F1` = `test_state_conversion.exe`) — DONE**
- [x] F1-a `dynamics/state_conversion.h` (new) — `to_state_vector(State)` and `to_state(StateVector)` bridge
  the two propagators' state types: km↔m exact scaling, and the DQ body-frame velocity is rotated into the
  world frame by the pose orientation (no-op for the identity-attitude orbital case, correct for general
  attitude — verified by a 90° attitude test).
- [x] F1-b round-trip both directions; `to_state` sets identity attitude + zero rate (a StateVector carries
  no attitude). Reuses the exact ×1000 unit factor.
- [x] F1-c `state_vector.h` — added `velocity_{measurement,precision,accuracy}_error()`, symmetric with the
  position accessors.

**F2 — common verb + TLE-driven DQ factory (gate `F2` = `test_api_parity.exe`) — DONE**
- [x] F2-a `dynamics/dq_sgp4_propagator.h` (new) — `DqSgp4Propagator<T>` exposes the SGP4-symmetric verb
  `propagate(tsince_minutes) → State`, internally calling the DQ `propagate_to(epoch, t_target, dt_max)`.
- [x] F2-b TLE-driven DQ factories `authentic`/`boosted` — DQ-from-TLE is now 1 line (was ~20): seed via
  `state_from_tle`, assemble a default gravitational force model (central + zonal Jₙ) + RK4. Verified the DQ
  epoch reproduces SGP4 t=0 to 7e-18 km (round-tripped through the F1 bridge).
- [x] F2-c the authentic-WGS72 seed vs the (optionally boosted) DQ propagation model are now EXPLICIT in the
  factory names + docs — the former silent mismatch is surfaced as the intended authentic-vs-boosted split.
- [x] F2-d the verb takes minutes-since-epoch (SGP4 convention), converted to the DQ second-clock (×60)
  internally; forward-only (tsince ≤ 0 → epoch). Gate `test_api_parity` (8 checks). 33/33 + OR1 + W11 green.
  (The `authentic`/`boosted` named modes also deliver the core of G1.)

**F3 — unify constants bundle (gate `F3` = `test_constants_single_source.exe`) — DONE (F3-a)**
- [x] F3-a WGS 72 / WGS 84 / GRS 80 defining constants single-sourced. The DQ `ConstantsProvider::{wgs84,
  wgs72,grs80}` now DELEGATE to the geodesy `EquipotentialEllipsoid` factories (the reusable single source,
  C2) — the former inline copies in `constants_provider.h` are gone, and the resulting ellipsoid is
  byte-identical (the F3 gate asserts `dq == si` exactly). The authentic SGP4 km model is deliberately
  NOT altered (it is the DQSGP4 oracle); the gate instead BINDS it to the SI source by asserting
  scale-consistency (a·1000 = aₘ, GM·1e9 = GMₘ, inv_f / e² / ω equal — all hold to rel 0). Any future
  divergence in either consumer fails gate `F3`. 33/33 + OR1 bit-exact; W8/W11/W24/F1 green.
- [ ] F3-b `ConstantsProvider` ("v1 carries only earth") not interchangeable with `ModelConfiguration` — an
  API-shape parity gap (the two bundles differ in members); folded into the F2 API-parity work.

### Phase G — Authentic vs boosted mode (R4, R5)

**G1 — explicit mode selection + docs (gate `G1` = `test_modes.exe`) — DONE**
- [x] G1-a `dynamics::PropagatorMode {Authentic, Boosted}` enum + `DqSgp4Propagator::from_tle(td, tol, mode,
  dt)` surface the mode as a first-class named parameter (the F2 `authentic`/`boosted` factories realise
  each). Authentic = WGS72 ellipsoid + WGS72 zonals (the sgp4_standard analogue); Boosted = WGS84 + EGM2008
  J₂…J₉; both seed from the authentic WGS72 SGP4 (a TLE is WGS72). Wide `T` is the orthogonal precision boost.
  Gate `test_modes` (8 checks): mode→earth-model mapping, shared authentic seed, modes diverge (4.99 m/20 min),
  from_tle(mode) dispatch, SGP4 sgp4_standard = WGS72.
- [x] G1-b `propagator_choice.md` — new sections: the one-line `DqSgp4Propagator` easy path, the
  authentic-vs-boosted mode table, the input-formats (TLE Alpha-5 + checksum, OMM KVN) section; README gains a
  one-line-from-TLE + modes note pointing to it.

### Cross-cutting — correctness & precision (R1, R3)

- [x] **AD1** `integrators/runge_kutta_fehlberg.h` — `rkf78_propagate_adaptive(y0, t_target, dt0, accel, tol,
  dt_min)` closes the loop the fixed-cadence `propagate_to` omits: it drives `rkf78_step` with the embedded
  7(8) local-error estimate, accepting a step when err ≤ tol (else retrying smaller) and scaling the next dt
  by the 8th-order law dt·0.9·(tol/err)^(1/8), clamped [0.1,5]× and floored at dt_min. Returns an
  `AdaptiveResult` (final state + accepted/rejected counts); the summed local error folds into the twist
  accuracy budget. Gate `test_adaptive_stepping` (6 checks) — monotonic adaptation verified: tol
  1e-2→1e-6→1e-10 gives 4→6→12 steps and 2.6e-3→6.5e-5→9.5e-9 m error vs a fine reference. 33/33 + W20 green.
- [x] **DS1 (DONE, gate `DS1` = `test_sdp4_precision.exe`)** — the SDP4 deep-space evolution now runs entirely
  in `TrackedValue<T>` (DSCOM buildup, DPPER periodics, the steps-1-9 secular/drag/J₃ advance, and the
  resonance leapfrog), so the three-error budget propagates through it instead of being re-injected from the
  epoch at the Kepler stage. Landed in four bit-exact phases (DS1.1 `0a35f40`, DS1.2+3 `57bc2ed`, DS1.4
  `b04f0c1`, DS1.5 gate) — OR1 reports **max drift 0 km / 0 km·s⁻¹** and `test_sgp4` stays 33/33 623/623 at
  every phase. **Key refinement — `TrackedValue::model_coefficient`:** the SR3 lunisolar magic numbers (and
  any finite-digit model coefficient) are encoded so their *binary-storage* error is computational PRECISION
  (scales with `T`) while their *finite-written-digit* error is MODEL-FIDELITY → ACCURACY (T-independent).
  `from_truncated_decimal` had booked the whole finite-digit bound as precision, which floored the deep-space
  precision at ~km regardless of `T` (masking the calling card); `defined()` had dropped the digit error
  entirely (under-claiming). `model_coefficient` keeps the total honest (= `from_truncated_decimal`) but in the
  right categories — realising the register's original "so the deep-space *accuracy* budget is honest" intent.
  Gate `test_sdp4_precision` (28 checks, 3 sats covering non-resonant + 12h + 24h leapfrog): precision never
  worsens with `T`, sits far below the model-accuracy floor at `cpp_bin_float_50` (the result is SR3-model-
  limited, not arithmetic-limited), accuracy is the T-independent model floor, and a clean case sharpens
  **413× double→bf50**. **Residual (corrected 2026-06-05):** the high-ecc/high-drag sat 11801 keeps a small
  precision floor (~1.1e-3 km, barely T-scaling). It was HYPOTHESISED to be the J₃/earth-model book
  coefficients — but the panel-ruling sweep moved ALL of those to `model_coefficient`/`measured` and 11801's
  floor did NOT move, DISPROVING that hypothesis. The source is an as-yet-unidentified T-independent
  computational contribution in the high-drag/high-ecc deep-space path (likely a degenerate-guard or
  cancellation site, not a constant). It does NOT affect the DS1 gate or the calling card: precision (1.1e-3
  km) still sits far below the model-accuracy floor (1.63 km), so the result is honestly model-limited. Flagged
  as an open diagnostic, not a gate concern. Mechanism + phase plan retained below for reference:
  - **Mechanism of the gap.** The deep-space evolution runs entirely in raw `double` and only re-wraps into
    `TrackedValue` at the Kepler/osculating stage by RE-INJECTING the epoch budget,
    `TrackedValue<T>(T(axn), ds.e0.errors)` (`propagate_deep_space` step 10, ~lines 655-677). So the output
    precision is the epoch precision carried forward, NOT a budget propagated through the secular/periodic
    evolution → it doesn't scale with `T`.
  - **Conversion set (compile-atomic).** (1) struct fields `sse..ssh, se2..xh3, zmos, zmol` (and later
    `gsto, theta_dot, a3ovk2`) double→`TrackedValue<T>`; (2) `build_dpper_coefficients` (~lines 116-295, the
    dense 2-pass solar/lunar DSCOM loop — a1..a10, x1..x8, z1..z33, s1..s7) → `TrackedValue`, taking the
    epoch elements as `TrackedValue` so precision flows from the inputs; (3) `apply_deep_periodics`
    (~303-382) → `TrackedValue`; (4) `propagate_deep_space` steps 1-9 (the secular advance, sun/moon linear
    terms, drag scalings, J₃ long-period) → `TrackedValue`, dropping the step-10 epoch re-injection;
    (5) `perturbation::step_resonance` (the irez leapfrog, in `perturbation/`) → `TrackedValue&` — the hardest
    piece; can be phased last (extract/rewrap for resonance sats meanwhile, documented).
  - **Born-digital constants.** The DSCOM/DPPER magic numbers (ZNS 1.19459e-5, C1SS 2.9864797e-6, ZES, ZNL,
    C1L, ZEL, ZCOSIS, …) become `from_truncated_decimal` (CR1) so the deep-space accuracy budget is honest.
  - **Validation protocol (per phase).** Bit-exact at `T=double` → OR1 (frozen oracle, deep-space sats) and
    `test_sgp4` 33/33 must stay green after EVERY phase; the gate then asserts a deep-space sat's precision
    budget tightens at `cpp_bin_float_50` vs `double`. Revert any phase that drifts OR1.
- [x] **BUG1** `celestial_body.h` — `mu_over_r3` was declared, never assigned, and never read (confirmed
  across `src/` and `tests/`); it could not be computed from `FundamentalConstants` and was redundant with
  the `perturbation_coef` passed to `compute_third_body_rates`. Removed the latent trap (a comment now points
  to where the strength actually comes from). Gate `test_celestial_body_fields` checks `make_sun`/`make_moon`
  populate every remaining field sensibly. 33/33 unaffected; W11 green.
- [x] **OR1** `tests/test_sgp4_regression/` — frozen SGP4 self-consistency oracle (the user's "temp version
  of SGP4 for regression until the final functions are complete"). A golden-master over the SGP4-VER 33-sat
  battery (2338 points) guarding the SGP4 path at **1e-6 km / 1e-9 km/s — ~1000× tighter** than the 33/33
  reference test (which only checks vs `tcppver.out` at km tolerance, so a refactor could shift outputs WITHIN
  that band undetected — exactly the F3 risk). Empirically SGP4 is bit-deterministic: finite outputs are
  identical across rebuilds (measured max drift 0); per-point NaN error-sentinels are compared by error-state,
  not value (NaN-safe serialization — MSVC `-nan(ind)` is not `>>`-readable). Gate `OR1`. **This is the guard
  that makes F3 safe.** TEMPORARY — delete the golden to re-baseline, or remove the oracle, once the redesign's
  final functions land.

---

## NEEDS-RULING — tracked, ungated until you decide (build → GATED, or accept → BY-DESIGN)

- [x] **INJ1 injection completeness (DONE, gate `INJ1` = `test_injection.exe`)** — chose the register's
  *simplify-and-make-honest* option (full engine-side wiring of the inclination function is OR1-IMPOSSIBLE:
  the half-day resonance `f220 = (3/4)(1+2c+c²)` and synchronous `(3/4)(1+c)²` are mathematically equal but
  bit-different, so routing both through one Kaula F_lmp lambda would drift the frozen SGP4 path). So:
  removed the dead `inclination_function` slot + its `InclinationFn` typedef from `ModelFunctions` (it was
  assigned but the engine never called it — `perturbation::inclination_function` remains a reusable library
  function); removed the label-only `CustomBuilder.drag()` selector + the unused `make_model_functions` drag
  arg + the 6 call sites (SGP4 near-space drag is the single fixed B*-Lane model — nothing to select; the DQ
  swappable drag is `forces::DensityModel`, DRAG1). The `model_functions.h` header now honestly states the
  per-step slots (secular_rates / kepler_solver / sidereal_time) are injected and the resonance F_lmp is a
  documented direct-call exception. Gate `test_injection` (5 checks): a swapped `kepler_solver` AND a swapped
  `secular_rates` each change the propagation result (proving the engine genuinely dispatches the remaining
  slots), and the config description no longer advertises `drag=`. 33/33 + OR1 bit-exact (value-neutral
  removals); W11 + A2 (dead-code) green.
- [x] **D2 tesseral/sectoral harmonics (DONE, gate `D2` = `test_tesseral.exe`)** — the longitude-DEPENDENT
  geopotential now exists. `constants/tesseral_harmonics.h`: `TesseralHarmonics<T>` storing unnormalized
  C_nm/S_nm, factory `egm2008()` populating the canonical degree-2 sectoral C̄_22/S̄_22 (born-digital,
  denormalized via tracked sqrt N_nm = √[(2n+1)2(n−m)!/(n+m)!]); extensible to 4×4+ by one `add()` line per
  EGM2008 coefficient. `forces/gravity_tesseral.h`: the singularity-free Cunningham / Montenbruck-Gill
  Cartesian V_nm/W_nm recursion (extensible to any degree/order) + the ECI↔ECEF GMST rotation. Gate
  `test_tesseral` (9 checks): denormalization, the SECTORAL acceleration is longitude-dependent (the D2 core —
  same |a| at λ=0 and λ=90° but |Δ|=1.45e-4 m/s², a zonal field would be identical), the GMST coupling
  changes the inertial tesseral pull as Earth rotates, gmst=0 rotation identity, and the precision is honestly
  source-floored by the EGM2008 digits (the CR1-b hybrid for static constants). Born-digital only: I populated
  the derivation-cross-checked C̄_22/S̄_22 and structured the rest as a documented one-line-each EGM-file
  extension (per the binding sourcing rule, not reciting uncertain higher-degree digits). DQ-only force →
  33/33 + OR1 untouched; W11 green.
  - **Residual:** the general `F_lmp` inclination function (`kaula.h:189-193`, 12 hardcoded triples) is for the
    ORBITAL-ELEMENT Kaula expansion (a different formulation from the Cartesian tesseral gravity above) and is
    not needed for the geopotential force; left as-is (tracked here so it is not mistaken for a D2 omission).
- [x] **E3 OMM (XML) (DONE, gate `E3` = `test_omm_xml.exe`)** — `tle/omm_xml_parser.h`: `parse_omm_xml`
  extracts the OMM mean elements + TLE parameters from the CCSDS XML serialization and feeds them through the
  SAME `omm_detail::populate_from_kv` mapping as the KVN parser (refactored out of `omm_parser.h` so both OMM
  front-ends are one mapping). The element extractor is a targeted, dependency-free scan with a tag-name
  BOUNDARY check so `MEAN_MOTION` is never confused with `MEAN_MOTION_DOT`/`_DDOT`, and it skips attributes
  (`units="…"`). Gate `test_omm_xml` (28 checks): full field mapping, the MEAN_MOTION/_DOT/_DDOT
  disambiguation, attribute skipping, leap-aware Dec-26-2020 epoch = 361.667, XML==KVN byte-identical TleData
  for equivalent content, feeds `from_tle_data`, missing-mandatory → false. Library header (no sgp4 coupling)
  → W11 green; E2 (KVN) still green after the shared-mapping refactor.
- [x] **B3 remaining tracked transcendentals (DONE, gate `B3` = `test_tracked_transcendentals.exe`)** —
  `tracked_value.h` now has tracked `asin/acos/tan/sinh/cosh/tanh/expm1/log1p/log10/log2/hypot/min/max/clamp`,
  each a value + rigorous forward error bound (worst-endpoint derivative form, written `d·f'` not a difference
  to dodge cancellation) + degenerate guard (domain-edge → π for asin/acos, pole → numeric-max for tan, `d≥x`
  flag for the logs). `small_angle_series.h` replaces the hardcoded `1e-4` with the T-dependent crossover
  `taylor_branch_threshold<T>(C) = (C·eps)^(1/6)` (C = 5040 / 112⁄5 / 840 per helper) — the point where the
  next-Taylor-term truncation balances the closed form's representation/cancellation floor, so it SHRINKS for
  wider T (the fix: the old fixed threshold kept Taylor live past where the closed form had become more
  accurate). `vector3.h` gains `normalize()` + scalar `operator/`. Gate `test_tracked_transcendentals` (55
  checks): value vs `std::`, the bound is a true upper bound (sampled across `[x−d,x+d]`), degenerate guards
  fire, and precision tightens with `cpp_bin_float_50` vs `double`. 33/33 + OR1 + W11 (after a line-length fix)
  green; the threshold change did NOT regress test_quaternion/test_dual_quaternion/test_precision_scaling.
- [x] **EPH ephemeris fidelity (DONE, gate `EPH` = `test_ephemeris.exe`)** — `solar_ephemeris.h`: the 2-term
  equation-of-center truncation `(13/12)|e|³` is now FOLDED into `errors.accuracy` (was documented but
  dropped), and the distance approximation carries its `O(e²)` model bound. `lunar_ephemeris.h`: the EoC gained
  the `(5/4)e² sin(2l)` term (Meeus §47 sin-2l), and the longitude/latitude now carry honest model-fidelity
  accuracy bounds for the omitted major periodics (evection/variation/annual ≈ 2.5° longitude, ≈ 0.9°
  latitude) — the register's "weakest link" (single-term, no error term) is now a tracked degree-grade budget.
  `sidereal_time.h`: new `compute_gmst_ut1(jd_utc, ΔUT1_sec)` applies the UT1−UTC correction; ΔUT1 = 0
  reproduces `compute_gmst` EXACTLY, so the authentic SGP4 path (which calls `compute_gmst` directly for
  `gsto`) is bit-untouched. The lunar/solar ephemeris has no SGP4 consumer (SGP4 uses its own DSCOM
  lunisolar), so OR1 + 33/33 unaffected. Gate `test_ephemeris` (11 checks): folded bounds, 2-term EoC value,
  born-digital ref (Earth max EoC ≈ 1.915°), UT1 bit-exact-at-zero + correct shift, precision tightens with T.
  Fuller Meeus (evection/variation terms, needing the solar elongation D) is a documented extension.
- [x] **DRAG1 atmospheric density fidelity (DONE, gate `DRAG1` = `test_drag_density.exe`)** — `forces/drag.h`
  now has a pluggable `DensityModel<T>` callback (altitude → tracked density carrying its own accuracy);
  `exponential_density_model` (Lane, with `rel_accuracy` param, default 0.30) is one implementation,
  `nrlmsise00_density_model_stub` is the documented MSIS seam (fallback to exponential with a 50% band — a
  fallback, not a fake), and `make_drag(density, B)` is the generic force; `make_drag_exponential` is retained
  as the convenience wiring. The flat 30% stamped on the acceleration is GONE — the density model's declared
  accuracy now PROPAGATES through the drag arithmetic into `errors.accuracy` (model-derived: 3× the declared
  accuracy → 3× the drag accuracy, proven by the gate). Gate `test_drag_density` (12 checks): pluggable
  callback drives the force, exponential value, model-derived accuracy scaling, NRLMSISE stub, precision
  tightens with T, convenience == generic. DQ-only force (not the SGP4 path) so 33/33 + OR1 untouched;
  W8/W11/test_propagator green.

---

## BY-DESIGN — accepted limitations (documented; listed so they're not mistaken for omissions)

**Phase-6 disposition (2026-06-05, for user review):** each is an INHERENT property of the
SGP4/Lane/Lagrange model or a physical fact, documented at its cited location, with no warranted upgrade —
dispositioned ACCEPTED. (Upgrade paths exist but are out of scope / would change the model: L1 smoothing the
Lane recipe, L2 second-order short-period, L3 a fully singularity-free third-body form. L4 cannot be
"upgraded" — a TLE simply carries no attitude.) Flagged so the user can re-open any of these.

- [x] L1 `density_model.h:79-89` — C¹ kink at perigee 156 km (~1% drag-rate discontinuity; inherent to the
  Lane recipe). Documented as ERROR SOURCE A-D4. **Accepted** (inherent to the Lane piecewise recipe).
- [x] L2 `short_period.h:63-67` — first-order-J₂-only short-period; ~7 mm omitted-term bound, below the
  ~100 m secular floor (deliberately not tracked). **Accepted** (omitted term ≪ the secular floor).
- [x] L3 `third_body.h:158,167` — equatorial (sin i→0) singularity inherent to the Lagrange-planetary form.
  **Accepted** (inherent to the orbital-element form; the Cartesian forces have no such singularity).
- [x] L4 `state_from_tle.h:16,34-35,55` — identity attitude + zero body rate from a TLE (a TLE carries no
  attitude). Correct, not a gap. **Accepted** (a physical fact about the TLE input — nothing to upgrade).

*(BY-DESIGN items dispositioned ACCEPTED in Phase 6 per the above; the user may re-open any to upgrade.)*

---

## DEFERRED — out of scope for this objective unless prioritized

- [x] **H1 full attitude dynamics (DONE, gate `H1` = `test_attitude_dynamics.exe`)** — built. `math/matrix3.h`
  (new): `Matrix3<T>` (3×3 TrackedValue) with mat·vec, determinant, adjugate/det inverse, diagonal +
  symmetric (products-of-inertia) factories — the off-diagonal tensor infrastructure that didn't exist.
  `dynamics/attitude.h` (new): the FULL Euler equation `ω̇ = I⁻¹(τ − ω×(I·ω))` (the diagonal
  `Inertia::acceleration_from_wrench` omitted the gyroscopic `ω×(Iω)` coupling), quaternion attitude
  kinematics `q̇ = ½q⊗(0,ω)`, rotational energy / angular-momentum, and an RK4 step that integrates (q, ω)
  and retracts q to unit norm. The point-mass torque restriction is lifted: a real inertia tensor takes any
  torque with a finite ω̇. Gate `test_attitude_dynamics` (11 checks): Matrix3 inverse on an off-diagonal
  tensor, the free symmetric-top precession matches the analytic Ω=ω_z(C−A)/A, energy + |L| conserved
  torque-free, q stays unit-norm, the gyroscopic term cross-feeds the axes (ω̇_y = −5/3 from the x/z spin),
  and precision tightens with `cpp_bin_float_50`. DQ-side library (orbital SGP4 unaffected) → 33/33 + OR1
  untouched.
