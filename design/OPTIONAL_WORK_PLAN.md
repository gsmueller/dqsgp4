# Optional Work Plan — DQSGP4 fidelity extensions

**★ Governing doc for the OPTIONAL work** surfaced by the 2026-06-06 code sweep. The DQSGP4 propagator and
the Constants Honesty initiative are COMPLETE (`run_acceptance` 67/67; OR1 0 km; 33/33 623/623; Debug+Release
both build; sln folder-organized). Nothing here is required for correctness — each item is a real fidelity
upgrade that was deferred because it needs external resources, a large port, or is a documented extension.

This doc is the source of truth; like the constants plan, it lists scope DECISIONS for the user up front,
then phases. Pick the scope dials and I execute top-to-bottom, each item gated + committed.

## Standing invariants (every item)

1. **★ NO PERCEIVED FIDELITY — verify against an INDEPENDENT reference oracle (user directive 2026-06-06).**
   *"Whatever work we take on needs an adequate methodology to test it to verify that all code and results are
   believable. We can easily create a seam by adding perceived fidelity."* A model that merely *looks* more
   accurate — or is only self-consistency-checked — is **perceived fidelity, a seam, and is REJECTED**. Each
   phase must reproduce the output of an **independent, trusted oracle** (in-repo Vallado `AstroLib` sun/moon,
   `MSIS_Vers` density; published born-digital tables; the source coefficient itself for a provenance check)
   to a **stated tolerance**, gate-enforced. The fidelity claim (and its declared `errors.accuracy`) is only
   believable once it MATCHES the oracle — never asserted from the model's own internals. This generalises the
   OR1 pattern (frozen SGP4 oracle) to every phase below; the oracle is named per-phase and the gate fails if
   the port drifts from it. Establishing the oracle comes FIRST, before the fidelity is claimed.
2. **OR1 bit-exact + 33/33.** Every item below is OFF the authentic SGP4 path (DQ-side force / modern-preset
   ephemeris / provenance-only). None may perturb the frozen SGP4 oracle — verified per item.
3. **Born-digital sourcing only.** Coefficients come from a citable born-digital source — the in-repo Vallado
   `MSIS_Vers.cpp`/`EopSpw.cpp`, Meeus *Astronomical Algorithms* (pymeeus realises it verbatim), IERS/NGA
   tables — **never recited from memory**. Large coefficient tables are transcribed from the in-tree source.
4. **Each item → a named `run_acceptance` gate** (extend an existing gate or add one), and the three-error
   budget is tracked (the model DECLARES its accuracy; it propagates).
5. **Build via the `.sln`** (Debug+Release); validate with the full `run_acceptance`.

## Reference oracles (the methodology, per phase — invariant #1)

| Phase | Independent trusted oracle | Believability check |
|---|---|---|
| O1 ephemeris | Vallado `AstroLib` sun/moon (in-repo, born-digital) and/or JPL Horizons values | port's λ, β, r reproduce the oracle to the model's stated grade (e.g. ≤ arcmin after main terms) |
| O2 EGM Jₙ | the born-digital C̄ₙ₀ table itself (NGA/ICGEM/IERS) | Jₙ == −√(2n+1)·C̄ₙ₀ to the stored digits; mismatch ⇒ transcription error to fix |
| O3 density | Vallado `MSIS_Vers` (NRLMSISE-00); the published USSA76 table | port's ρ reproduces the oracle (≤ ~0.1 % MSIS / table-exact USSA76) across sample inputs |
| O4 sidereal ratio | IERS published ERA rate / a known GMST | generated ratio matches the published value to its digits |

**Decisions (user, 2026-06-06):** O2 source = *I source a born-digital table* (D3-b). Execute the *full plan
O1→O4*. Atmosphere (D1): the binding answer is the methodology above — fidelity is only "added" once it is
verified against `MSIS_Vers` / the USSA76 table; staged USSA76→MSIS lets each increment be proven believable.

## Scope decisions (please choose; my recommendation in brackets)

- **D1 — Atmospheric-density scope (Phase O3).** (a) Full **NRLMSISE-00** port (space-weather-responsive, the
  real driver of drag uncertainty; large); (b) **US Standard Atmosphere 1976** tabulated model (static,
  born-digital, far smaller — improves the Lane ~30 % *mean* profile to ~few %, but has NO solar/geomagnetic
  response); (c) **both, staged** — USSA76 now, NRLMSISE-00 later. *[Recommend (c): land USSA76 as a quick
  honest win, then decide on the full MSIS port — the space-weather response is where the real value is, but
  it carries the F10.7/Ap-input complexity below.]*
- **D2 — NRLMSISE-00 space-weather inputs (if O3 includes MSIS).** (a) **Caller-supplied** F10.7 / F10.7avg /
  Ap as parameters on the model (no data-file dependency; tests pass representative values); (b) port `EopSpw`
  and ship a CelesTrak `SpaceWeather-All.txt` data file (operational, but needs a periodically-updated file).
  *[Recommend (a): keeps the port pure and gateable; a data-file reader is a separable follow-on.]*
- **D3 — EGM2008 Jₙ source (Phase O2).** The n=5–9 cross-check needs full-precision normalized C̄ₙ₀.
  (a) You provide the 104 MB EGM2008 spherical-harmonic file (`earth-info.nga.mil`); (b) I source a smaller
  born-digital low-degree table (ICGEM / IERS) for n=5–9. *[Recommend (b) if a citable born-digital low-degree
  table is locatable; else (a).]*
- **D4 — Lunar periodic-term count (Phase O1).** (a) The **main terms** (evection, variation, annual equation,
  + the largest few — ~arcminute longitude); (b) the **fuller Meeus §47** set (~dozens of terms — ~arcsecond).
  *[Recommend (a): the top 3–4 terms remove ~2.1° of the current ~2.5° budget; diminishing returns after.]*

## Phases (recommended sequence — self-contained first)

### Phase O1 — Lunar/solar ephemeris periodics (EPH extension)
- **Now:** `lunar_ephemeris.h` uses a 2-term equation-of-center and carries a *lumped* ~2.5° longitude /
  ~0.9° latitude model-accuracy budget for the omitted periodics. `solar_ephemeris.h` is already ~0.01° (its
  omitted terms are tiny) — the win is almost entirely lunar.
- **Do:** add the Meeus §47 major periodics — evection 1.274°·sin(2D−l), variation 0.658°·sin(2D), annual
  equation −0.186°·sin(M), (+ the next few per D4), and the latitude terms (0.280°·sin(F+l), 0.173°·sin(F−2D),
  …). These need the **solar elongation D = L_moon − L_sun** and solar anomaly M — i.e. a sun↔moon coupling:
  compute the solar mean longitude/anomaly and pass D, M into `compute_lunar_position` (or compute jointly).
  Replace the lumped 2.5°/0.9° accuracy with the residual of the *still*-omitted terms (tracked, now ~arcmin).
- **Source:** Meeus *Astronomical Algorithms* Ch. 47 (born-digital; cross-check vs pymeeus). NOT from memory.
- **Gate:** extend `EPH` (`test_ephemeris`) — added-term values vs Meeus reference, accuracy budget shrinks
  (longitude accuracy < ~0.1° after the main terms), precision tightens with T. Off-SGP4-path → OR1 untouched.
- **Effort:** moderate (~½ session). **Risk:** low (born-digital, DQ-side).

### Phase O2 — EGM2008 Jₙ (n=5–9) provenance cross-check (closes R04)
- **Now:** `model_selector.h` (`wgs84_precise`/`grs80` zonals) stores Jₙ for n=5–9 as `measured(...)`; the
  comment-quoted C̄ₙ₀ reproduce the stored Jₙ only to 4–5 sig figs (R04: a transcription artifact, or the
  comment lost decimals). C̄₃₀/C̄₄₀ (IERS Table 6.2) are solid; n=5–9 are the open ones.
- **Do:** obtain full-precision C̄ₙ₀ (n=5–9), verify Jₙ = −√(2n+1)·C̄ₙ₀ to full precision, fix any
  transcription error, and re-encode the zonals as `model_coefficient` (per the constants initiative —
  digit-floor → accuracy, storage → T-precision) rather than `measured` with a coincidental σ.
- **Source:** per D3. **Gate:** extend `D1`/`GAL1` — Jₙ == −√(2n+1)·C̄ₙ₀ to the stored digits; honesty scan
  already enforces the encoding. Modern-preset only (the SGP4 WGS72 zonals are frozen) → OR1 untouched.
- **Effort:** small (~¼ session once sourced). **Risk:** sourcing-gated (D3).

### Phase O3 — Atmospheric-density fidelity (DRAG1 stub → real model)
- **Now:** `DensityModel<T>` is altitude→density; `nrlmsise00_density_model_stub` is the Lane exponential
  with a 50 % band. The seam exists; the model does not.
- **Do (per D1):**
  - **O3a — input model.** Define `AtmosphericState<T>` (geodetic altitude, latitude, longitude, UT/day-of-year,
    and — for MSIS — F10.7, F10.7avg, Ap) and a richer model type, OR extend `DensityModel`. `make_drag`
    derives alt/lat/lon/time from the `State` + epoch; space-weather per D2.
  - **O3b — the model.** USSA76 (piecewise-layer, born-digital, ~hundreds of lines) and/or NRLMSISE-00
    (templatized + tracked from the in-repo `MSIS_Vers.cpp` NRLMSISE-00 routine — ~2–3 k lines + large
    coefficient tables, transcribed born-digital). Declared model accuracy: USSA76 ~few %, NRLMSISE-00 ~15 %
    — propagates into the drag `errors.accuracy` (replacing the 30–50 % band).
  - **O3c — space weather (MSIS only, per D2).**
  - **O3d — gate.** New `test_atmosphere`: density vs reference (USSA76 published table values; NRLMSISE-00 vs
    Vallado `MSIS_Vers` output at sample (alt,lat,lon,F10.7,Ap) to ~0.1 %), accuracy budget, precision with T.
- **Source:** in-repo `MSIS_Vers.cpp` / `EopSpw.cpp` (born-digital); USSA76 from the 1976 standard's
  born-digital layer table. DQ drag force only (SGP4 uses its own Lane drag) → OR1 untouched.
- **Effort:** USSA76 ~½ session; full NRLMSISE-00 multi-session (transcription volume + templatizing +
  tracked error through a large model). **Risk:** medium (volume; the tracked budget through MSIS).

### Phase O4 — `sidereal_ratio` generative (minor, sourcing-gated)
- `sidereal_time.h` `sidereal_ratio = model_coefficient("1.00273790935")`. Generate as 1 + 1/(tropical-year
  in days) **iff** a born-digital IERS Earth-Rotation-Angle rate is sourced. The value would shift ~6e-8 (the
  adopted ratio includes precession), so it is NOT value-preserving — gate update required. Tiny; low priority.

## Effort / risk summary

| Phase | Effort | Risk | Blocking input |
|---|---|---|---|
| O1 lunar periodics | ~½ session | low | none (Meeus born-digital) |
| O2 EGM2008 Jₙ check | ~¼ session | low | D3 source |
| O3 USSA76 | ~½ session | low | none (born-digital table) |
| O3 NRLMSISE-00 | multi-session | medium | D1/D2 |
| O4 sidereal ratio | ~minutes | low | born-digital IERS rate |

**Recommended order:** O1 → O2 → O3 (USSA76 first, then MSIS if chosen) → O4. O1 and the USSA76 slice of O3
are fully self-contained (no external input) and deliver real, gateable fidelity now; O2 and full MSIS are
gated on the sourcing decisions above.
