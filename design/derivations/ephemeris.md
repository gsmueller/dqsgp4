# Generic Ephemeris — Theory Note (L3)

**Module L3 of the professional-library re-architecture** (`design/PROFESSIONAL_LIBRARY_PLAN.md`). THEORY
FIRST: this note develops the analytical-ephemeris theory exhaustively; the `OrbitalModel` /
`Ephemeris` types and their functions (§5) are its direct consequence. No ephemeris code is rebuilt before
this note stands ([[feedback_theory_first_library]]). Fidelity is verified against an INDEPENDENT oracle —
astropy / ERFA + Vallado `AstroLib` + JPL DE430 — established FIRST, never asserted
([[feedback_no_perceived_fidelity]]). Builds on L1 (`Epoch`/`TimeScale`, [[time_scales_and_epochs]]) and L2
(`Frame` transforms, [[reference_frames]]). **Subsumes optional O1** (`design/OPTIONAL_WORK_PLAN.md`).

## 0. Why this module exists (grounding: the existing state)

The library's lunar/solar ephemeris is the textbook "modify-first" symptom (`PROFESSIONAL_LIBRARY_PLAN` §2):

- `ephemeris/{celestial_body,solar_ephemeris,lunar_ephemeris}.h`: `CelestialBody<T>` bundles geocentric mean
  elements (built from the SR3 1970s `FundamentalConstants` in `astronomy/solar_system.h`);
  `compute_solar_position` / `compute_lunar_position` return ecliptic (λ, β, distance) from a 2-term equation
  of center plus, for the Moon, `β ≈ i·sin(L − Ω)`.
- It is **UNWIRED**: no `src/` file consumes `compute_*_position` (only `test_ephemeris`). The authentic SDP4
  deep-space path uses the *separate* SR3 `DerivedOrbitalElements` inline in DPINIT/DPPER — **not** this
  module. So the ephemeris is entirely **off the OR1 path** (verified: grep finds no consumer).
- It has **no frame or time-scale type**: `delta_t_min` is "minutes since reference epoch" with the epoch
  unstated; "ecliptic longitude" never says of-date vs J2000; no `Epoch<T>` (L1) appears.
- Its accuracy is **ASSERTED, not measured**: the lumped ~2.5° longitude / ~0.9° latitude lunar bounds and the
  ~0.01° solar bound are folded into `errors.accuracy`, but nothing ever compares an actual Sun/Moon position
  to truth. "Nobody could say which frame/epoch `compute_lunar_position` returns" (§2).

That is a *cause* (missing theory), not a bug to patch. Writing the ephemeris theory down yields a typed,
frame/scale-explicit, body-generic, source-pluggable, **astropy-verified** `Ephemeris` — and it is the natural
**consumer** of exactly the precession / nutation / obliquity_rotation that L2 deferred to "their L3 consumer".

## 1. The model and the conventions it pins

A low-precision analytical ephemeris of a body B is: **mean elements advancing linearly in a named time
scale**, plus a **truncated periodic series**, yielding the body's direction (λ, β) and distance r in a named
frame. A claim missing any of the following conventions is a defect:

- **Time scale & epoch.** The independent variable is an `Epoch<T>` (L1), carrying its `TimeScale`.
  `mean_element(t) = element(epoch) + rate·(t − epoch)`, both sides in the SAME named scale. The **SR3
  instance** uses the 1900-based day (days since 1899-12-31 12:00 = J1900 = L1 `days_since_1900`) on the SGP4
  UT1≈UTC convention; a **modern instance** uses TT Julian centuries since J2000 (L1
  `centuries_since_j2000`). The epoch is never a bare number ([[time_scales_and_epochs]]).
- **Frame of (λ, β).** Mean **ecliptic and equinox OF DATE** — the mean elements and periodics are referred to
  the moving ecliptic/equinox. Reaching the equatorial propagation frame is `Rx(+ε(t))` (obliquity, SC1
  generator) to equator-of-date, then precession/nutation to TOD/MOD/GCRF — the L2 transforms (§3).
- **Units & ranges.** λ, β and all angles in radians, wrapped to a stated range; distance per body (Sun in AU;
  Moon in Earth-radii or km — pinned at the instance); rates in rad per scale-unit.
- **Origin = geocentric** (the body's *apparent* orbit seen from Earth's centre). The Sun's "orbit" is Earth's
  heliocentric orbit reflected (e ≈ 0.0167, β☉ ≡ 0 by definition); the Moon's is the physical geocentric orbit
  (e ≈ 0.0549, i ≈ 5.145° to the ecliptic). `orbit_eccentricity` is the Earth–body system's, not the body's.

## 2. The orbital-model algebra (the structure)

Mean angles advance linearly in the scale (n = 2π/period):

```
 mean anomaly      l(t) = l₀ + n·(t − epoch)
 mean longitude    L(t) = L₀ + ṅ_L·(t − epoch)        (Moon: SIDEREAL rate ≠ anomalistic n — distinct periods)
 node              Ω(t) = Ω₀ + Ω̇·(t − epoch)          (regresses: Ω̇ < 0)
```

**Equation of center** — the Kepler/Lagrange inversion of `M = E − e sin E`, expanded as the Fourier–Bessel /
Lagrange series in the mean anomaly:

```
 ν − M = 2e sin M + (5/4)e² sin 2M + (13/12)e³ sin 3M + O(e⁴)            (Meeus §25/§47)
```

It converges for e < e_Laplace ≈ 0.6627 (both bodies are ~2 orders below this → rapid, uniform convergence).
Truncating after K terms costs **accuracy = Σ_{k>K} |c_k e^k|**, dominated by the first omitted term — the
same generator pattern as obliquity/eccentricity ([[feedback_constants_generative_or_bounded]]): the bound
tightens as K grows. Then

```
 λ = (L  or  M + ϖ) + (ν − M) + Σ longitude periodics
 β = i·sin F + Σ latitude periodics            (Moon; F = L − Ω = argument of latitude)
 r/a = 1 − e cos M + (e²/2)(1 − cos 2M) + O(e³)                          (Meeus §25)
```

**Lunar periodics** (the major Meeus §47 terms; physical origin = the **Sun's perturbation of the lunar
orbit**, so the arguments are Delaunay combinations that COUPLE to the Sun):

```
 evection         ≈ 1.274° · sin(2D − l)      solar tide modulating lunar eccentricity
 variation        ≈ 0.658° · sin(2D)          solar tangential perturbation over a synodic cycle
 annual equation  ≈ 0.186° · sin(l′)          Earth's orbital eccentricity modulating the solar tide
   D  = mean elongation (Moon − Sun)    l = lunar mean anomaly
   l′ = solar mean anomaly              F = lunar argument of latitude
```

D and l′ require the **solar** mean longitude/anomaly, so the lunar model consumes the solar model — the
`Ephemeris` evaluates both. (The amplitudes are **born-digital from Meeus Ch. 47**, to be VERIFIED term-by-term
against astropy/pymeeus at implementation — NEVER asserted from memory here, [[feedback_no_perceived_fidelity]].)

## 3. The component transforms (theory + provenance)

- **Mean elements & rates** — each a `{value @ epoch, rate}` pair in a named scale. The **SR3 instance** =
  the 1970s `FundamentalConstants` (measured σ). A **modern instance** = VSOP87 (Sun) / ELP2000 (Moon) or
  Meeus mean-element polynomials @ J2000 (finite-digit fit ⇒ `model_coefficient`: digits → accuracy, storage →
  T-precision). Two INSTANCES of one structure — the element **source is pluggable**, exactly like
  `ModelSelector` (this is what makes O1 principled: SR3 and a modern set are independently verifiable, never
  silently mixed).
- **Periodic amplitudes A_k** — born-digital (Meeus/VSOP/ELP); each term is `A_k·sin(c_D D + c_l l + c_{l′} l′
  + c_F F)` with integer arg-coefficients. Keeping the top-K by amplitude is a truncation; **accuracy = Σ
  omitted |A_k|** (tightens with K).
- **Obliquity ε(t)** — the SC1 generator `obliquity_iau2006` (`astronomy/obliquity.h`), oracle-verified;
  ecliptic↔equator is `Rx(±ε)` (L2's `obliquity_rotation`, consumed here).
- **To GCRF/J2000** — precession `P(t)` + nutation `N(t)` (L2, consumed here), each ERFA-verified. For
  SGP4-grade third-body work nutation may be dropped with a tracked ~arcsec accuracy.

## 4. The three-error budget

- **measurement** — element σ (SR3 measured values) or adopted-constant σ.
- **precision** — T-representation, tightening with a wider T (the calling card); it must flow through the
  whole mean-element + series + frame chain.
- **accuracy** — series truncation (Σ omitted periodic amplitudes + the EoC next-term), tightening with kept
  terms K; PLUS the model-instance fidelity floor (1970s SR3 elements vs a modern set). Today's *lumped* ~2.5°
  lunar longitude bound becomes the **explicit residual after the named periodics** (→ ~arcmin), and that
  residual is **VERIFIED to majorize the measured astropy difference** — not asserted (§6).

## 5. The abstraction that falls out (the functions, as consequences)

```cpp
// One element as a linear model in a named scale; one periodic row on Delaunay args.
template<typename T> struct MeanElement { math::TrackedValue<T> at_epoch, rate; };
template<typename T> struct PeriodicTerm { math::TrackedValue<T> amplitude; int cD, cl, clp, cF; };

// A body's analytical model: mean elements @ epoch (a SCALE) + a truncated series.
// Generic over the BODY and over the element SOURCE (SR3 vs Meeus/VSOP/ELP).
template<typename T>
struct OrbitalModel {
    std::string body;  TimeScale scale;  Epoch<T> epoch;
    MeanElement<T> mean_anomaly, mean_longitude, node;
    math::TrackedValue<T> eccentricity, inclination, arg_perigee;
    std::vector<PeriodicTerm<T>> longitude_terms, latitude_terms;   // truncation → accuracy
    EclipticState<T> evaluate(const Epoch<T>& t, const OrbitalModel<T>* sun) const;  // (λ,β,r) ecliptic-of-date
};

// Body × source instances (pluggable like ModelSelector):
template<typename T> OrbitalModel<T> sun_sr3();   template<typename T> OrbitalModel<T> sun_meeus();
template<typename T> OrbitalModel<T> moon_sr3();  template<typename T> OrbitalModel<T> moon_meeus();

// Ephemeris = OrbitalModel + Frame(L2) + TimeScale(L1): a body's position in ANY frame.
//   = frame_transform(EclipticOfDate → target, t) · [ r(t) · direction(λ(t), β(t)) ]
template<typename T>
math::Vector3<T> body_position(const OrbitalModel<T>& m, const Epoch<T>& t, Frame target,
                               const OrbitalModel<T>* sun = nullptr);
```

The existing `CelestialBody` + `compute_*_position` are precisely the **SR3 instance** of `OrbitalModel` plus
its ecliptic `evaluate()`; they become one instance — now typed by `Frame`/`Epoch`, astropy-verified, with the
Cartesian `body_position` and the frame chain added, and the major lunar periodics filled in (O1).

## 6. Oracle & verification (independent truth, never self-consistency)

Established FIRST (invariant #1 — the oracle precedes the fidelity claim):

| Claim | Independent oracle / check |
|---|---|
| Sun/Moon geometric ecliptic λ, β, r at sample epochs | **astropy** `get_body` geometric ecliptic-of-date; cross-check **JPL DE430** (in-repo `sgp4_references/.../sunmooneph_430t12.txt`) |
| equatorial / GCRF Cartesian after the frame chain | **astropy** GCRS/ICRS position to a stated tolerance |
| each periodic term | **Meeus Ch. 47 / pymeeus**, term-by-term |
| truncation accuracy bound | residual vs astropy ≤ Σ omitted amplitudes — the bound MUST **majorize** the measured residual |
| precession / nutation / obliquity rotations | **ERFA (IAU SOFA)** (via L2) |
| precision scaling | double vs `cpp_bin_float_50` |
| (unchanged) SGP4 path | **OR1 0 km** — ephemeris is UNWIRED/off-path |

The stated grade is **gate-enforced**: solar λ ≤ ~arcmin, lunar λ ≤ a few arcmin after the named periodics,
vs astropy. A claimed bound that does NOT majorize the measured astropy residual is **perceived fidelity →
REJECTED**. astropy/ERFA and JPL DE are independent of any in-repo code, so a match is proof.

## 7. Migration (off the OR1 path; freely rebuilt, oracle-gated)

The ephemeris module is unwired (no `src/` consumer; SDP4 uses the separate SR3 `DerivedOrbitalElements`), so
L3 is **not** value-preserving-constrained: it is rebuilt into `OrbitalModel`/`Ephemeris` with `Frame`/`Epoch`
typing, the SR3 element set preserved as one instance (its numbers unchanged), a modern Meeus/J2000 instance
added, the major lunar periodics added (O1), and the whole verified against astropy — **OR1 untouched
(verified)**. Implementation is bottom-up, each unit isolated + oracle-gated: `OrbitalModel` structure → SR3
instance (reproduces today's λ,β) → the astropy oracle harness → periodic terms (accuracy shrinks,
astropy-verified) → frame chain (consume L2 precession/nutation/obliquity) → `body_position` vector. Gate:
extend `EPH` (`test_ephemeris`) + a new astropy-oracle comparison gate.

## References (born-digital)

- Meeus, *Astronomical Algorithms* (1998), Ch. 25 (solar), Ch. 47 (lunar) — the periodic series; `pymeeus`
  realizes them verbatim ([[feedback_born_digital_latex]]).
- VSOP87 (Bretagnon & Francou) / ELP2000 (Chapront-Touzé & Chapront) — the full analytical theories the
  truncations come from.
- IERS Conventions (2010) / IAU 2006 — precession / nutation / obliquity (consumed via L2).
- astropy `coordinates` + ERFA (IAU SOFA); JPL DE430 (`sgp4_references/vallado_celestrak/datalib/
  sunmooneph_430t12.txt`); Vallado `AstroLib` sun/moon — the verification oracle.
