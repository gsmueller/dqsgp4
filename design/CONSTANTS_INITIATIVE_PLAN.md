# Constants Honesty & Generative Initiative — Governing Plan

**★ This document is the MAIN PROJECT THREAD for the constant-encoding initiative.**
This session owns it end-to-end (user directive, 2026-06-06: *"Create a plan that
captures you being responsible for the entire project."*). A prior parallel
session ("DQSGP4 execution plan completion") drove the gravity-field slice and
then stopped mid-flight; its work is merged, coherent, and green (see Baseline).
All further work on constant encoding is tracked here.

## The directive

Two user directives, one principle:
- 2026-06-05: *"Constants that are not true by definition need to have their
  series-based precision and accuracy tracked."* → a non-defined constant that is
  the leading terms of a known series is **generated** from that series, carrying
  series-truncation **accuracy** + representation **precision**, not stamped as a
  decimal.
- 2026-06-06: *"All constants need accuracy and precision tracked; each Jₖ value
  is the result of a series truncation."* → no physical/model quantity may hide
  behind `defined()` (which zeroes accuracy and claims exact-by-convention).

Restated: **every constant carries an honest three-error budget, and every
constant that is *derived from a formula* has a lambda generator that can be
cranked up (more series terms, wider type T) to sharpen it.**

## The encoding decision tree (binding — the 2026-06-05 panel ruling, by PROVENANCE)

| Provenance | Constructor | precision | accuracy |
|---|---|---|---|
| Exact by convention — datum geometry (a, 1/f, ω), epoch JD, integer day/sec counts | `defined()` | 0.5 ULP (scales w/ T) | 0 |
| Adopted physical/astronomical quantity (σ = adoption bound) **or** real published σ | `measured(v,σ)` | 0.5 ULP | 0 (σ in measurement) |
| Finite-digit model-fit coefficient, no σ — SR3 lunisolar, gravity Jₖ/C̄ₙₘ, GM | `model_coefficient(s)` | 0.5 ULP (scales w/ T) | 0.5·decimal_ulp (T-independent) |
| **Leading terms of a known SERIES** — obliquity, eccentricity, … | **generator** `f<T>(t, n_terms)` | representation of the Horner sum (scales w/ T) | series-truncation Σ_{k≥n}|cₖ||t|ᵏ (tightens w/ n) |
| Numerical guard / source-truncated `defined` value | `from_truncated_decimal(s)` | max(0.5 ULP, 0.5·decimal_ulp) | 0 |

The generator row is the calling card: precision tightens with T *and* accuracy
tightens with kept terms — the only encoding where BOTH error channels improve on
demand. `defined()` is reserved for things that are exact *because we declared
them so* (a metre, a flattening, the J2000 epoch), never for a measured or
fitted physical number.

## Completed (verified green this session)

- **Obliquity** εA(t) generator — `src/astronomy/obliquity.h` (SC1, committed `4100fe2`).
- **Earth eccentricity** e(t) generator — `src/astronomy/earth_orbit.h` (SC2-astro, THIS session).
  Born-digital: Meeus *Astronomical Algorithms* Ch. 25 Eq. 25.4 (VSOP87), per
  Julian century; confirmed verbatim vs the pymeeus source. Wired into the
  `iau_2006` preset; J2000 value sharpened 0.0167086 → 0.016708634.
- **Gravity J₂/GM** `defined()`→`model_coefficient()` in `model_selector.h`
  (make_ellipsoid) and `equipotential_ellipsoid.h` (grs80/wgs72 presets), +
  zonal SC2 test section, + a J_k/GM `defined()` honesty scan (prior session, merged).
- **Zonal/tesseral C̄** already generative (−√(2n+1)·C̄ₙ₀, tracked √) + `model_coefficient`.

## Baseline (this session, `run_acceptance -SkipVerifiers`)

`I0.build` PASS · `I0.sgp4` 33/33 623/623 PASS · `OR1` bit-exact PASS · all module
tests PASS · every ExeGate incl. `C1`/`SC1`/`CR1B` PASS · `W18` zero unclassified
literals PASS. **Full sweep (verifiers included): 67/67 gates, VERDICT COMPLETE;
OR1 drift 0 km; 33/33 623/623.** All 8 phases done — see Status: COMPLETE below.

## Phases

Each phase: value-preserving on the SR3/OR1 path (bit-exact) unless a born-digital
source legitimately sharpens a modern-preset value (then update its gating test);
verify with `.sln` build + 33/33 + OR1 + the relevant gate; commit per phase.

- **Phase 1 — eccentricity gate.** ✅ `test_series_constants` §E: value at t=0
  (0.016708634), secular drift, accuracy tightens with n_terms, precision tightens
  with T — parity with the obliquity SC1 block. (Committed with SC2, `dc54846`.)
- **Phase 2 — resonance/deep-space model coefficients.** ✅ `resonance.h` g300
  `6.60937` `defined()`→`model_coefficient` (SR3 synchronous-resonance fit
  coefficient). It was the ONLY `defined()` in `resonance.h`/`deep_space.h`.
  Value-preserving: OR1 0 km, 33/33, DS1. Committed `e4cd4b9`.
- **Phase 3 — sidereal-time honesty.** ✅ `sidereal_time.h`: c₂/c₃ (σ == decimal
  ULP, mis-filed) → `model_coefficient`; `sidereal_ratio` `defined()`→
  `model_coefficient` (a derived ratio incl. precession, so a 1+1/tropical-year
  generation would shift the value ~6e-8 — not value-preserving, deferred). c₀/c₁
  kept `measured` (σ = 100× ULP, genuine). Added to the honesty scan. OR1 0 km,
  EPH, CR1B. Committed `21b8457`.
- **Phase 4 — consumer delegation.** ✅ `src/test_series.cpp` (3× inline WGS84) →
  `EquipotentialEllipsoid::wgs84()` factory (one honest source; GM now measured);
  verified compile+link+run 29/29. `src/test_sgp4_ver.cpp` J₂/GM → `model_coefficient`
  (kept the inline `from_J2` km demo; verified compiles). `src/main.cpp` left as-is
  — it does NOT re-type J₂/GM (GM already `measured`, a/1f/ω correctly `defined()`);
  it is a pedagogical 4-param-construction demo. All three are orphaned (not in the
  sln/build). Committed `532a03f`.
- **Phase 5 — honesty-check completeness.** ✅ Broadened the J_k/GM `defined()`
  scan from the two gravity files to the ENTIRE `src/` tree (68 files); clean +
  false-positive-free. Documented exclusions (SR3 `sgp4_standard` uses `measured`
  with adoption-bound σ; model_coefficient encodings carry accuracy). Committed `8211310`.
- **Phase 6 — secular-series forward (sourcing-gated).** ✅ Sidereal year
  (`solar_anomaly_period_days`) KEPT `measured` (already directive-compliant:
  adopted astronomical quantity, σ = adoption bound, never `defined()`).
  Documented why no secular series: its drift is far below the σ≈43 µs adoption
  bound over any practical span (mirrors `lunar_eccentricity`). Doc-only; W11/CR1B.
  Committed `5e2943c`.
- **Phase 7 — full sweep + commit + memory.** ✅ Full `run_acceptance` (85 Octave +
  2 Python verifiers + all gates) → **67/67 COMPLETE**; OR1 0 km; 33/33 623/623.
  MEMORY.md consolidated.

- **Phase 8 — constant-audit completeness (the strongest reading of directive 2).**
  ✅ Audited all 134 encoding calls (27 `defined()`, 40 `measured()`, 66
  `model_coefficient()`, 1 `from_truncated_decimal()`). (a) Every `defined()` is a
  genuine convention EXCEPT `src/main.cpp` π/4 frozen as `defined("0.7853981633974483")`
  → fixed to generative `quarter_pi<T>()`. (b) Every `measured()` σ accounted for:
  model-file σ genuine (≠ ULP); the σ==ULP cases are all SR3 `sgp4_standard`
  (`solar_system.h`) + IERS/geodesy adopted astronomical measurements (σ = adoption
  bound, genuine, bit-exact) — documented exclusions. (c) Replaced the J_k/GM
  value-range scan with a **convention-allowlist** scan over all 68 `src/` files:
  every `defined()` must match a registered exact-by-convention value or it is
  flagged (proven to flag J2/GM/ecc/obliquity/π-decimal). Committed `d5f21fb`.
  run_acceptance **67/67 COMPLETE**; OR1 0 km; 33/33 623/623.

## Status: COMPLETE

All 8 phases done. Every constant in `src/` is now one of: a registered
exact-by-convention `defined()` (accuracy 0, precision scales w/T); a `measured()`
adopted/observed quantity with a genuine σ; a `model_coefficient()` finite-digit
fit (digit-floor → accuracy, storage → T-precision); or a **generator** whose
precision (wider T) and accuracy (more terms) both sharpen on demand. Each
category is gate-enforced (`CR1B` convention-allowlist + σ-honesty +
precision-purity; `SC1` series generators). The directive — *every* constant an
honest three-error budget, *each J_k* a series truncation — is fully realized and
guarded. Standing guardrail remains the DQSGP4 suite (33/33 + OR1 + 67/67).

## Out of scope / guardrails

- The SR3 `sgp4_standard()` frozen 1970s constants stay `measured` — they are
  adopted model inputs whose σ is an adoption bound (panel ruling), and they must
  reproduce SR3 bit-for-bit (OR1). Do NOT "modernize" them.
- DQSGP4 itself is 66/66 COMPLETE; this initiative extends constant *honesty*, and
  the DQSGP4 regression suite (33/33 + OR1 + acceptance) is the standing guardrail.
- Born-digital sourcing for every series coefficient ([[feedback_born_digital_latex]]).
