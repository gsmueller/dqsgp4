# DQSGP4 Completion Roadmap — authentic + boosted dual-mode propagator

**Status:** ACTIVE (opened 2026-06-04). Supersedes the "ROADMAP COMPLETE" resting state of
`recursive-forging-rain.md` (W0–W24), which remains done and green. This is a NEW initiative on a
NEW user objective.

**Branch:** `session/2026-04-23`, MAIN repo. **Completion = `pwsh tools/run_acceptance.ps1` prints
`COMPLETE`** with the new gates below added. Regression invariant unchanged: `test_sgp4` 33/33, 623/623
on every code-touching commit; one session = one branch; commit per item.

---

## 1. The objective (user, 2026-06-04, verbatim intent)

> Survey the codebase and ensure all logic is complete for the DQSGP4. The API for SGP4 and DQSGP4 should
> be similar in nature. Functions should be reusable, isolated, and capable of supporting multiple levels
> of increased accuracy and precision. There should be a mode that directly models SGP4 as advertised, but
> an ability to increase the precision/accuracy of the constants, which potentially extends things such as
> the modeled geopotential constants. It should be capable of using the updated TLE formats, and updated
> earth models present in the Astronomical Almanac.

Decomposed requirements:

- **R1 — Logic complete.** No stubs, dead injection slots, orphaned modules, or comments that lie about
  the code, in the DQSGP4 path or its dependencies.
- **R2 — API parity.** `sgp4::Propagator` and `dynamics::Propagator` should be "similar in nature":
  a common propagate-to-time verb, a shared State↔StateVector conversion, a TLE-driven DQ entry point.
- **R3 — Reusable / isolated / multi-precision.** Functions templated on `T`, isolated, no duplicated
  constant definitions; precision/accuracy tracked through *every* computation.
- **R4 — Authentic mode.** A mode that reproduces SGP4 bit-faithfully "as advertised" (exists today as
  `sgp4::Propagator` + WGS72; must be named/selectable and protected by the 33/33 gate).
- **R5 — Boosted constants / geopotential.** Selectable higher-accuracy constants that *extend* the
  modelled field (odd + higher even zonals wired into a force; eventually tesseral/sectoral).
- **R6 — Updated TLE formats.** Alpha-5 actually decoded (not stored raw), checksum validation, and the
  modern OMM (CCSDS KVN/XML) mean-element message Celestrak/Space-Track now publish.
- **R7 — Updated earth models.** Real Astronomical-Almanac / IAU-IERS earth-model + celestial constants,
  replacing the `model_selector.h` "TODO: populate from specific Almanac editions" stub.

---

## 2. Survey findings (2026-06-04) — evidence-based, mapped to requirements

Six read-only subsystem surveys + three lead reads. The codebase is **green but not complete**: the
incompleteness splits cleanly into (a) *stale scaffolding / doc-lies* that make working code look unfinished
and (b) *genuine functional gaps* — and the genuine gaps map almost 1:1 onto R2–R7. The prior roadmap's
33/33 fidelity gate never exercised the half-built seams, which is exactly why it stayed green.

### 2a. Dead/stale scaffolding & doc-lies (serve R1)

- `src/sgp4/model_functions.h:208` — "Aoki et al. (1982) GMST polynomial. (Placeholder — not yet
  implemented.)" — **FALSE**; `compute_gmst` is fully implemented (`astronomy/sidereal_time.h`) and wired at
  `model_functions.h:264-269`, used in `deep_space.h:448`.
- `src/sgp4/model_functions.h:271-285` — `drag_coefficients` lambda is a hard `throw` ("not yet extracted
  as injectable lambda"); the comment says drag is "inline in near_space.h" — also stale (it lives in
  `atmosphere/drag_coefficients.h:107`, which `near_space.h:211` calls directly, bypassing the slot). The
  `.drag(name)` builder (`model_selector.h:195`) is a silent no-op (`make_model_functions` ignores it).
- `src/sgp4/model_functions.h:224-230` — `inclination_function` lambda assigned but never invoked.
- `src/sgp4/model_functions.h:287` — description says "Halley Kepler solver"; `standard_sgp4()` installs
  **Newton** (the Halley body lives only in `model_selector.h`).
- `src/sgp4/precomputed.h` — `KaulaTable` never instantiated (orphaned); doc-table references phantom types
  `PrecomputedModel` / `PrecomputedSatellite` that do not exist.
- `src/sgp4/sgp4_propagator.h:44-80` — `NearSpaceCoefficients`, `DeepSpaceCoefficients`, `PropagationState`
  declared but never used (real state is `NearSpaceInit`/`DeepSpaceInit`/`ResonanceState`); the
  `PropagationState` doc-comment ("the ONLY mutable state") is wrong.
- `src/orbit/modified_kepler.h:17-20,41,80` — standalone Newton/Halley solvers never called (live solvers
  inline-duplicated in `model_functions.h:235` and `model_selector.h:433,471` — triplicated); docstring
  advertises a "Householder (quartic)" solver that does not exist.
- `src/dynamics/wrench.h:19` — stale "(REQ-PR-2, forthcoming)"; the summation is in
  `propagator.h:126-132`.
- `src/perturbation/brouwer.h:33` — cites `design/derivations/009_brouwer_rates.md (TODO)`; file absent.

### 2b. Precision infrastructure — missing tracked transcendentals (serve R3, R5)

- `src/math/tracked_value.h` tracks exactly: `+ − * /`, `sqrt sin cos atan atan2 abs fmod`. **Not tracked:**
  `exp`, `log`, `pow`, `cbrt`, `tan`, `asin`, `acos`, hyperbolics, and `min/max/clamp` on `TrackedValue`.
- Consequence: `src/forces/drag.h:104-110` computes density `ρ = ρ₀·exp(arg)` in raw `T` ("`std::exp` is
  not yet a `TrackedValue<T>` operation"); precision error of the exponential is silently dropped — so
  dialing up `T` does **not** improve drag's tracked precision. This is precisely where "track precision AND
  accuracy through every computation" fails.
- `src/orbit/element_recovery.h:122-147` and `src/geodesy/equipotential_ellipsoid.h:172-188` open-code
  tracked-Newton cube-roots because there is no reusable tracked `pow`/`cbrt`.

### 2c. Earth models / Astronomical Almanac (serve R7)

- `src/sgp4/model_selector.h:397-414` — `iau_2006` / `almanac_2010..2025` presets are stubs: they return
  SR3 constants with **only** obliquity + solar period overridden ("TODO: populate from specific Almanac
  editions"). The headline stub for R7.
- `src/geodesy/equipotential_ellipsoid.h` — a complete level-ellipsoid *engine* but **no named presets**;
  `constants_provider.h` offers only `wgs84/wgs72/grs80`. No current IERS/Almanac GM, obliquity, rotation,
  precession catalog.
- `src/sgp4/model_selector.h:288-300` — EGM2008 J₅–J₉ carry a provenance TODO (R04): stored values
  reproduce the quoted C̄ₙ₀ only to 4–5 sig figs.

### 2d. Geopotential fidelity (serve R5)

- DQ force path implements **J₂ only** (`forces/gravity_zonal.h:74`); the generic Jₙ (n≥2) force is a
  doc-only stub (`gravity_zonal.h:111-116`). Higher even zonals (J₄, J₆) are available from the constants
  but never added to a force; J₃ is faked by a `|J4/J2|` proxy (`gravity_zonal.h:88`).
- `src/geodesy/equipotential_ellipsoid.h:287-289` — `J2n(n)` produces **even zonals only**; odd zonals
  (J₃, J₅) are not representable in the geodesy layer at all.
- **No tesseral/sectoral (C_nm, S_nm) field anywhere** — `model_selector.h:168-172` `A(n,m)` returns 0 for
  m≠0; the geopotential is strictly axisymmetric.

### 2e. TLE / OMM (serve R6)

- `src/tle/tle_parser.cpp:134-136` — Alpha-5 is **stored raw, not decoded** (`parse_catalog_number` just
  returns the field); the A-Z→10-33/skip-I,O decode in the comment is unimplemented. No checksum
  validation (col 68 parsed in the map but never checked). No 9-digit catalog despite the `.cpp` claim.
- No OMM (CCSDS Orbit Mean-elements Message, KVN or XML) ingestion — the modern format. Epoch→JD is a
  "Simplified JD computation" (`tle_parser.h:119`).

### 2f. API parity & reuse (serve R2, R3)

- No common propagate verb: SGP4 `propagate(tsince_minutes) → StateVector` (km/TEME, stateful) vs DQ
  `propagate_to(state, t_target, dt_max) → State` (m/inertial, functional). No `State↔StateVector`
  converter exists. No TLE-driven DQ constructor (the bridge seeds a state but the user still hand-builds
  4 objects + lambdas; ~20 lines vs SGP4's 3).
- WGS72 constants are typed in **two** places (`constants_provider.h:73-82` km-units vs
  `model_selector.h:326-345` — duplication / drift risk). The DQ `ConstantsProvider` ("v1 carries only
  earth") and SGP4 `ModelConfiguration` are not interchangeable.
- The TLE→DQ bridge hard-codes `"sgp4_standard"` (WGS72) for the seed while the user typically injects
  `ConstantsProvider::wgs84` — seed and propagation ellipsoids can silently disagree
  (`state_from_tle.h:39`).

### 2g. DQ core — genuinely complete, one bounded deferral (serve R1)

- The DQ algebra (`dual_number/quaternion/dual_quaternion`), `state/pose/twist/wrench`, the propagator, and
  `lie_advance_pose` SE(3) retraction across all four integrators are **complete** for translational
  (orbital, point-mass) propagation and diagonal-inertia attitude. `state_from_tle.h` fully bridges a TLE
  to a DQ `State`.
- Sole substantive deferral: `src/dynamics/inertia.h:13-15` — general off-diagonal inertia tensor needs a
  `Matrix3<T>` type that doesn't exist; only diagonal/principal-axis inertia inverts. Full coupled-axis
  attitude dynamics is out of scope for orbital DQSGP4 (deferred phase).
- RKF7(8) computes an embedded per-step error but `propagate_to` runs fixed-cadence and ignores it
  (`propagator.h:152-163`) — adaptivity exists at the integrator boundary but isn't closed-loop.

---

## 3. The Work-Queue

Selection rule: lowest-numbered not-Done item. Each item has an objective **Definition of Done** that is a
mechanical `run_acceptance` gate (a new `ExeGate`/`VerifierGate`/`ScriptGate`, or a `git grep` assertion).
Phases are ordered so a *truthful baseline* precedes feature work and *enablers* precede consumers.

### Phase A — Truthful baseline (R1). Make "logic complete" literally true.

- **A1 — Stale-comment truth pass.** Correct/delete every comment that lies about implemented code: the
  GMST "not yet implemented" (`model_functions.h:208`), the drag "inline in near_space.h" note, the
  Newton/Halley description mislabel (`:287`), the `wrench.h:19` "forthcoming", the `brouwer.h:33` TODO doc
  cite. **DoD:** new `tools/check_no_stale_stubs.ps1` ScriptGate greps the specific dead phrases at the
  specific sites and exits nonzero if any remain; `test_sgp4` 33/33.
- **A2 — Dead-scaffolding disposition.** For each orphan, decide wire-vs-delete and act: delete the phantom
  `PrecomputedModel`/`PrecomputedSatellite` doc-types and unused `KaulaTable` (or wire it into deep-space);
  delete the unused `NearSpaceCoefficients`/`DeepSpaceCoefficients`/`PropagationState` structs; remove the
  non-existent "Householder" doc + de-duplicate the triplicated Kepler bodies into the
  `modified_kepler.h` functions (or delete the unused ones). **DoD:** a verifier asserts no
  declared-but-unreferenced type in the named files; `test_sgp4` 33/33; full module tests green.

### Phase B — Precision infrastructure (R3, R5). Enablers.

- **B1 — Tracked `exp` / `log`.** Add `TrackedValue<T>` overloads with rigorous error propagation; route
  `drag.h` density through tracked `exp`. **DoD:** `tests/test_tracked_transcendental` exercises exp/log
  error budgets + a verifier checks drag precision now scales with `T`; `test_sgp4` 33/33.
- **B2 — Tracked `pow` / `cbrt`.** Reusable tracked `pow`/`cbrt`; retire the open-coded Newton cube-roots
  in `element_recovery.h` and `equipotential_ellipsoid.h`. **DoD:** same test exe extended; `test_sgp4`
  33/33 (recovery path unchanged to tolerance).

### Phase C — Earth models from the Astronomical Almanac (R7).

- **C1 — Real Almanac/IAU-IERS constants.** Populate `FundamentalConstants` for `iau_2006`/`almanac_20xx`
  with actual published values (GM⊕, obliquity, precession, Earth rotation, lunar/solar rates) from the
  Astronomical Almanac / IERS Conventions; remove the "TODO: populate" stub. **DoD:**
  `verify_almanac_constants.m` checks each against its published value to stated σ; no "TODO: populate"
  remains; `test_sgp4` 33/33.
- **C2 — Earth-model preset catalog.** Named `EquipotentialEllipsoid`/`ConstantsProvider` presets for the
  current models (WGS84(G2296), GRS80, IERS/EGM2008 current, Almanac) with provenance. **DoD:**
  `test_earth_models` constructs each + checks defining params; `test_sgp4` 33/33.

### Phase D — Geopotential extension (R5).

- **D1 — Generic Jₙ zonal force.** Implement `gravity_Jn` (J₃, J₄, J₅, J₆…) in the DQ path, wired to the
  constants; add a real odd-zonal source so J₃ is not a `J4/J2` proxy. **DoD:** `test_force_models`
  extended with Jₙ acceleration checks vs an independent reference; `test_sgp4` 33/33.
- **D2 — Tesseral/sectoral harmonics (LARGE; depth TBD).** `C_nm/S_nm` field with Earth-fixed→inertial
  rotation; `A(n,m)`, m≠0 path; EGM coefficient ingestion. **DoD:** `verify_tesseral.m` + a force test;
  `test_sgp4` 33/33. *Scope to be confirmed before starting — may be capped at low degree/order.*

### Phase E — Updated TLE / OMM ingestion (R6).

- **E1 — Alpha-5 decode + checksum.** Decode the Alpha-5 leading letter to the numeric catalog id;
  validate the line-68 checksum; reject malformed sets with a reason. **DoD:** `test_tle` extended with
  Alpha-5 vectors (e.g. `T00001`→`270001`) + checksum pass/fail vectors; `test_sgp4` 33/33.
- **E2 — OMM (KVN) ingestion.** Parse CCSDS OMM Key-Value-Notation into the same `TleData`/`TleElements`.
  **DoD:** `test_omm` parses a reference OMM and matches the equivalent TLE elements to tolerance.
- **E3 — OMM (XML) ingestion (optional).** Same, for the XML serialization. **DoD:** `test_omm` extended.

### Phase F — API parity (R2, R3).

- **F1 — `State`↔`StateVector` converter.** Bidirectional, with the km↔m factor lifted into a reusable
  unit helper; attitude dropped/identity per direction. **DoD:** `test_api_parity` round-trips to ULP;
  `test_sgp4` 33/33.
- **F2 — Common propagate verb + TLE-driven DQ factory.** A shared "propagate to time t" entry point with
  one agreed time unit; `dynamics::Propagator::from_tle(...)` that reuses the seed's gravity preset (closing
  the WGS72/WGS84 mismatch). **DoD:** `test_api_parity` drives both propagators through the common verb;
  quickstart DQ setup collapses toward the SGP4 line-count.
- **F3 — Unify the constants bundle.** Build `ConstantsProvider` from a `ModelConfiguration` (or share one
  geodesy bundle); remove the duplicated WGS72 literals. **DoD:** `check_magic_numbers`/a verifier asserts
  the WGS72 defining set appears once; `test_sgp4` 33/33.

### Phase G — Authentic vs boosted mode (R4, R5).

- **G1 — Explicit mode selection + docs.** Name the modes: **authentic** (WGS72 + SR3, bit-faithful, 33/33)
  vs **boosted** (Almanac constants + extended geopotential + wide `T`). A single selector surfaces both on
  each propagator; `design/propagator_choice.md` + README updated. **DoD:** `test_modes` asserts authentic
  mode reproduces the 33/33 reference and boosted mode engages the extended field; full `run_acceptance`
  COMPLETE.

### Phase H — Attitude dynamics (R1 completeness; LARGE; deferred).

- **H1 — `Matrix3<T>` + full inertia tensor.** Off-diagonal inertia inversion for coupled-axis attitude.
  **DoD:** `test_attitude` integrates a torque-free asymmetric top (Euler) and conserves the polhode;
  `test_sgp4` 33/33. *Deferred unless the user prioritizes attitude — orbital DQSGP4 does not need it.*

---

## 4. Acceptance protocol

Each completed item adds its gate to `tools/run_acceptance.ps1` (mirroring the W1–W24 wiring). A milestone
is declared only after a FULL `run_acceptance` (no `-SkipVerifiers`) prints `COMPLETE`. The standing I0
invariant (sln builds; `test_sgp4` 33/33, 623/623; all module exes; all Octave/Python verifiers) gates
every commit. New test exes mirror `tests/test_force_models/` (fresh-GUID `.vcxproj` registered in
`sgp4.sln`); `src/` stays AUD-CC + magic-number clean.

**Ordering note:** A → B → C → D1 → E1 → F → G is the high-value / lower-risk spine; D2 (tesseral), E2/E3
(OMM), and H1 (attitude) are larger and their scope is confirmed before starting.
