# Static-atmosphere density: the piecewise-exponential table model (Vallado 8-4)

**Theory note (FIRST, before code) for replan item R1** (`design/PROFESSIONAL_LIBRARY_PLAN.md` §6) — the first
stage of the user-decided O3 atmosphere ladder (USSA76-grade static table now, NRLMSISE-00 later/R6). Governed
by [[feedback_theory_first_library]] and [[feedback_no_perceived_fidelity]]. The companion code is
`src/atmosphere/exponential_table.h` (+ the transcription generator `tools/gen_atmosphere_table.py`); the gate
is `tests/test_atmosphere` (ATM1).

---

## 1. What this replaces, and what it honestly claims

The DQ drag path (`forces/drag.h`, REQ-EF-7) evaluates a pluggable `DensityModel<T>`. The only real
implementation is the **single-band Lane exponential** `ρ₀·exp(−(h−h₀)/H)` with caller-chosen (ρ₀, h₀, H) and a
declared **30 % accuracy band** ("reproduces the U.S. Standard Atmosphere within ~30 %, up to ~100 % at solar
max" — `sgp4_near_earth_drag_theoretical_basis.md` §2 Lemma 2.2, §15 ERROR SOURCE A-D2); the
`nrlmsise00_density_model_stub` is the same with a 50 % band. One scale height cannot follow the thermosphere's
strongly altitude-dependent temperature, hence the 30 %.

This increment replaces that with the **published 27-band piecewise-exponential static atmosphere** — Vallado,
*Fundamentals of Astrodynamics*, Table 8-4 ("exponential atmospheric model", densities from USSA76 at low
altitude and CIRA-72 above) — whose born-digital realization ships **in-repo**:
`sgp4_references/vallado_celestrak/software/misc/pascal/ATMOSEXP.DAT` (27 rows: base altitude h₀ [km], nominal
density ρ₀ [kg/m³], scale height H [km], from 0.01 km to 900 km).

**The honest claim (two distinct layers, never conflated):**

1. *Code vs the published model*: the implementation **is** the model — within a band, ρ(h) = ρ₀·exp(−(h−h₀)/H)
   is the model's *definition*, not an approximation of it. So code-vs-model accuracy is the coefficients'
   **digit floor** only (`model_coefficient`: 4 significant digits → ~1e-4 relative) plus the band-edge
   chaining residual (§4, measured max 9.6e-5). This is what the gate verifies, table-exactly.
2. *The static standard vs the real thermosphere*: a static model **cannot represent solar/geomagnetic
   variability** — the real LEO density varies by factors of ~2–5 over the solar cycle (ERROR SOURCE A-D2,
   already documented in-repo). That layer is **NOT claimed and NOT silently buried in the accuracy channel**:
   it is documented here, the factory exposes a caller-widenable `representativeness` band for it, and it is
   the precise gap NRLMSISE-00 (R6) exists to close. Claiming "few-% real-atmosphere density" from a static
   table would be perceived fidelity; we claim few-% (actually digit-floor) fidelity *to the named published
   static model*, gate-verified.

This is the same convention the Lane band used (its 30 % was Lane-vs-USSA76, not Lane-vs-reality), with layer
(1) improved from 30 % to digit-floor.

---

## 2. The physics — why piecewise exponential is the natural static encoding

Hydrostatic equilibrium with the ideal gas law:

```
  dp/dh = −ρ g ,    p = ρ R* T / M     ⇒     dρ/ρ ≈ −(M g / R* T) dh ≡ −dh / H(T)
```

with `H = R*T/(Mg)` the **scale height**. Where T (and the mean molecular mass M, in the heterosphere) is
locally constant, density decays exactly exponentially with scale H. The real atmosphere's T(h) profile is
piecewise-smooth (troposphere lapse, stratospheric inversion, thermospheric rise toward T∞ ~ 1000 K), so a
**piecewise-exponential fit with per-band H** is the natural compact encoding of a static standard: each band's
H is the *effective* (fitted) scale height absorbing the in-band T and M variation. That is exactly what Table
8-4 is: 27 bands whose (ρ₀, H) were fitted to USSA76/CIRA-72 so consecutive bands **chain** — row i evaluated at
the next base altitude reproduces row i+1's ρ₀ (verified to the displayed digits, §4).

Evaluation, for geometric altitude h above the reference ellipsoid's equatorial radius (the convention
`forces/drag.h` already uses, `h = |r| − a`):

```
  band(h) = the row with the largest h₀ ≤ h          (h < h₀(first) clamps to the first band;
  ρ(h)    = ρ₀(band) · exp( −(h − h₀(band)) / H(band) )    h > 900 km extrapolates on the last band)
```

The branch on h is a **value comparison** (like the Lane perigee regimes in `atmosphere/density_model.h`); the
arithmetic inside the band is fully tracked (`TrackedValue` exp, gate B1), so precision and the coefficient
digit floors propagate into ρ automatically.

---

## 3. The data — born-digital, transcription mechanized

The 27 rows are transcribed from the in-repo `ATMOSEXP.DAT` by `tools/gen_atmosphere_table.py` (the
`gen_lunar_terms.py` pattern: the generator reads the source artifact and EMITS the C++ rows, so transcription
is mechanical and reproducible — never typed from memory). Encoding per the constants mandate
([[feedback_constants_generative_or_bounded]]):

- **h₀ in metres as exact integers** — the file's km values (0.010, 25.000, …, 900.000) are all exact at the
  metre: 10, 25 000, …, 900 000 m → `exact<T>`. (The file's first base is 0.010 km = 10 m, not 0 — transcribed
  verbatim; the source artifact is the authority. ρ(0) then evaluates on the first band, 0.14 % above sea-level
  1.225.)
- **ρ₀ as `model_coefficient` strings** (e.g. "3.899e-2") — 4 significant digits → the decimal ULP digit floor
  lands in accuracy; binary storage in precision (tightens with T).
- **H as `model_coefficient` strings in km, × exact(1000) to metres** — the digit floor scales exactly.

Units: altitude in **metres** in, density in **kg/m³** out — the `DensityModel<T>` contract of `make_drag`.

---

## 4. The oracle & gate (ATM1) — no perceived fidelity

The published table itself is the oracle (the "source coefficient itself" clause), with four independent
verification layers:

1. **Node exactness** — at every base altitude, ρ(h₀) must equal ρ₀ **bit-for-bit** (Δh = 0 ⇒ exp(0) = 1
   exactly ⇒ ρ = ρ₀·1). Any transcription slip in ρ₀ fails its row.
2. **Neighbor chaining** — row i evaluated at h₀(i+1) must reproduce ρ₀(i+1) to the table's construction grade.
   MEASURED (gen_atmosphere_table.py, all 26 edges): max relative residual **9.6e-5** — confirming the table
   was built by exact chaining with H rounded to 3 decimals (the residual is the H rounding). Gate tolerance
   2e-4 (2× the measured majorant). A wrong digit in ρ₀ *or H* breaks the chain at TWO rows at the ≥1e-2
   level — the strongest transcription catcher.
3. **Cross-repo consistency** — the 200-km row's ρ₀ is 2.789e-10 kg/m³, the very constant `test_propagator`'s
   drag phase has used all along (that value's provenance is hereby identified: it *is* this table's row).
4. **Lane-divergence demonstration** — the prior single-band Lane model (ρ₀ = 2.789e-10 @ 200 km, H = 50 km,
   the test_propagator parameters) agrees with the table near its anchor (~15 % at 220 km, inside its 30 %
   band) but drifts BEYOND the band away from it: MEASURED 41 % at 250 km, 56 % at 300 km. (The original
   "reproduces USSA76 within ~30 %" claim is per-band-fit; a single band anchored at one altitude cannot hold
   it globally — corrected here by measurement, and exactly why the per-band table replaces it.) The gate
   asserts both sides: within-band near the anchor, broken beyond it.

Plus the standing framework checks: the default declared band majorizes the measured chaining residual; the
digit-floor accuracy is nonzero; precision tightens with a wider T; and a `make_drag` wiring check (the model
plugged into the existing seam produces the hand-computed ½ρv²B at 400 km).

OR1-untouched: the SGP4 path keeps its own Lane power-law drag (`atmosphere/density_model.h`, frozen); this
model is DQ-side only, reachable through the `DensityModel` seam (and the R3 presets later).

---

## 5. Error budget summary

| channel | content |
|---|---|
| precision | binary storage of ρ₀, H + the tracked exp arithmetic — tightens with T |
| accuracy (model) | digit floors of ρ₀, H (`model_coefficient`, ~1e-4 rel) + the default band (1e-3, a 10× majorant of the measured 9.6e-5 chaining residual) |
| representativeness (documented, NOT defaulted-in) | static-standard vs real thermosphere: solar-cycle factors ~2–5 at LEO (A-D2); caller-widenable parameter; closed by NRLMSISE-00 (R6) |

The model: 27 born-digital rows, one tracked exponential, node-exact gates — replacing a 30–50 % band with a
digit-floor claim *about a named published model*, while saying out loud what a static atmosphere cannot know.
