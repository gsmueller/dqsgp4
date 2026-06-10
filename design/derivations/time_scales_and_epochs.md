# Time Scales & Epochs — Theory Note (L1)

**Module L1 of the professional-library re-architecture** (`design/PROFESSIONAL_LIBRARY_PLAN.md`). Theory
FIRST: this note develops the time/epoch theory exhaustively; the `Epoch`/`TimeScale` types (§7) are the
direct consequence. **No time/epoch code is written or changed before this note stands on its own**
([[feedback_theory_first_library]]). Every fidelity claim is verified against an independent oracle
([[feedback_no_perceived_fidelity]]) — here astropy `Time` / ERFA (IAU SOFA) and the in-repo conventions.

## 0. Why this module exists (the cause, not the symptom)

The codebase carries **three coexisting epoch conventions with no unifying abstraction**:

1. `elements.epoch_jd` — a bare Julian Date, the canonical satellite epoch (from the TLE/OMM).
2. `deep_space.h:472` — `epoch_jd − 2415020`, a **day count from the 1900 base** (JD 2415020.0), feeding the
   SR3 lunisolar mean elements (`ZMOS = 6.2565837 + 0.017201977·DAY`, etc., `solar_system.h`).
3. `sidereal_time.h` — a **Julian century from J2000** (`(jd − 2451545.0)/36525`) for GMST.

No type records *which time scale* a JD is in, and the TLE epoch (nominally UTC) is fed directly to GMST
(which is defined on UT1) — i.e. `ΔUT1 = 0` is assumed silently. Nobody could state, from the types, the
scale or epoch base of a given instant. That ambiguity is the *cause* behind the ephemeris confusion (an
"epoch" that was neither clearly J2000 nor clearly 1900). The cure is to write the theory down; a correct,
reusable `Epoch`/`TimeScale` abstraction then falls out, and the three views above become derived functions of
one canonical instant.

## 1. The physical time scales

An *instant* is a point on a time line; a *time scale* is a named, physically-defined parameterization of that
line. The scales relevant to satellite dynamics, with their defining relationships:

| Scale | Definition | Relationship | Status |
|---|---|---|---|
| **TAI** | International Atomic Time; the SI second on the rotating geoid | reference | realized, continuous |
| **TT** | Terrestrial Time; the independent argument of geocentric ephemerides | **TT = TAI + 32.184 s** | offset EXACT by convention (IAU 1991/2006) |
| **TDB** | Barycentric Dynamical Time; argument of barycentric ephemerides | **TDB = TT + Σ periodic** (amplitude ≈ 1.7 ms, dominant 1-yr term) | periodic series (IAU 2006 / Fairhead–Metris); convention-defined, series-truncated |
| **UT1** | Universal Time; proportional to Earth's rotation angle (ERA) | **UT1 = UTC + ΔUT1**, \|ΔUT1\| < 0.9 s | ΔUT1 OBSERVED (IERS Bulletin A/B), not predictable |
| **UTC** | Coordinated Universal Time; atomic rate, integer-second steps | **UTC = TAI − (leap seconds)** | leap seconds OBSERVED/decreed (IERS Bulletin C) |

Provenance taxonomy (cf. the constants initiative, [[feedback_constants_generative_or_bounded]]):
- **Exact by convention** — `TT − TAI = 32.184 s`; J2000 = JD 2451545.0; MJD = JD − 2400000.5; Julian century
  = 36525 d; 1 d = 86400 SI s. These are `defined()`-grade: accuracy 0, precision scales with `T`.
- **Series (generator)** — `TDB − TT`: leading periodic terms of a known series; truncation → accuracy,
  representation → precision. The same generator pattern as obliquity/eccentricity (SC1/SC2).
- **Observed** — `ΔUT1` (a measured quantity with a real σ, bounded \|·\| < 0.9 s by the UTC leap-second
  design), and the **leap-second table** (exact integers, but their *occurrence* is decreed — an accuracy
  bound if the table is stale). Both are `measured`-grade inputs the caller supplies (IERS), never invented.

## 2. Epoch representations

- **Julian Date (JD):** continuous count of days from −4712-01-01 12:00. The astronomical day starts at noon.
- **Modified Julian Date:** `MJD = JD − 2400000.5` (exact; starts at midnight).
- **J2000.0:** `JD 2451545.0 TT` = 2000-01-01 12:00:00 TT. The modern reference epoch. (Exact by definition.)
- **The "1900" base:** `JD 2415020.0` = 1899-12-31 12:00 (the "1900 January 0.5" of the older day-numbering).
  Exactly **36525 days = 100 Julian years** before J2000 (`2451545 − 2415020 = 36525`), which is precisely the
  Brouwer/SR3 lunisolar epoch base used in §0 item 2. Stating it as `J2000 − 100·(Julian century)` is the
  exact, convention-grounded relationship; the calendar date is its consequence.
- **Julian century from J2000:** `T = (JD_TT − 2451545.0) / 36525`. The argument of the IAU polynomials
  (obliquity SC1, GMST, precession L2).
- **Calendar / day-of-year → JD:** the Gregorian algorithm (Meeus Ch. 7; in-repo `tle_parser.h` specialised to
  month = 1 + day-of-year, verified correct in the Phase-8 audit).

## 3. Precision of a Julian Date (why a single double is not enough)

A JD near the present is ≈ 2.46×10⁶. In IEEE-754 double (53-bit mantissa, ulp ≈ 2.2×10⁻¹⁶ relative) the
representable resolution is `2.46e6 × 2.2e-16 ≈ 5.4e-10 d ≈ 47 µs`. Satellite work wanting sub-µs epoch
fidelity (and the wider-`T` calling card) therefore requires the **two-part split** `JD = jd_day + jd_frac`
(Vallado's `jdtdb`, `jdtdbF`; astropy's `jd1`, `jd2`), where `jd_day` carries the integer/half-integer part and
`jd_frac ∈ [−0.5, 0.5)` carries the sub-day part at full mantissa precision. Under `TrackedValue<T>` the
representation precision then scales with `T`, and the two-part form keeps `double` honest to ≈ 1e-11 d ≈ 1 µs.
This is a *precision* matter (the value is the same instant); it is independent of the *time-scale* (accuracy)
question of §1.

## 4. The three-error budget for an epoch and its conversions

- **TT ↔ TAI:** exact offset → conversion adds only representation *precision* (zero accuracy).
- **TDB ↔ TT:** the periodic series → *accuracy* = the truncation bound (Σ omitted terms; ≪ 1.7 ms), *precision*
  scales with `T`. Negligible for SGP4 (which uses neither), tracked for completeness.
- **UTC → TAI (leap seconds):** exact when the table covers the date; *accuracy* = a flag if the date is beyond
  the known table (cannot be predicted).
- **UT1 ← UTC (ΔUT1):** *measurement* σ from IERS; with the SGP4 convention `ΔUT1 = 0`, the unmodelled rotation
  is bounded by \|ΔUT1\| < 0.9 s ⇒ ≤ 0.9·(ω_E) ≈ **3.8 milli-degrees** of GMST — carried as an *accuracy* on any
  UT1-derived quantity (this is exactly the EPH bound already documented for `compute_gmst_ut1`).
- **JD representation:** §3 — *precision* only, scales with `T` (two-part form keeps double ≈ 1 µs).

## 5. The conventions in THIS library, reconciled (grounded in the code)

One canonical instant — the satellite epoch — viewed three ways:

```
            epoch  ──(JD, time scale = UTC≈UT1 per SGP4 convention, ΔUT1=0)
              │
   days_since_1900(epoch)  = JD − 2415020.0            → SR3 lunisolar mean elements   (deep_space.h:472)
   centuries_since_j2000(epoch) = (JD − 2451545.0)/36525 → GMST, IAU polynomials        (sidereal_time.h)
   (calendar/day-of-year → JD)                         → TLE/OMM epoch ingest           (tle_parser.h)
```

These are not three "epochs" — they are pure functions of one `Epoch`. The SGP4 path treats the TLE epoch JD as
UT1 = UTC (`ΔUT1 = 0`); that assumption becomes *explicit and overridable* under the abstraction, without
changing the frozen path (`ΔUT1 = 0` reproduces it bit-for-bit).

## 6. Oracle & verification (no perceived fidelity, invariant #1)

| Claim | Independent oracle / check |
|---|---|
| `TT−TAI = 32.184 s`, `J2000 = 2451545.0`, `MJD offset`, `century = 36525`, `1 d = 86400 s` | exact identities (convention) |
| calendar/day-of-year ↔ JD round-trip; J2000 ↔ 2000-01-01 12:00 | self-inverse + known date; cross-check vs **astropy `Time`** / ERFA |
| `days_since_1900`, `centuries_since_j2000` views | algebraic identity vs the bare JD; reproduce `deep_space.h`/`sidereal_time.h` **bit-exact** (value-preserving) |
| scale conversions (UTC→TAI→TT, ΔUT1) | **astropy `Time` scale transforms** (ERFA); leap-second table vs IERS |
| two-part JD precision | vs single-JD in `cpp_bin_float_50` (precision tightens with `T`) |
| **whole-library regression** | OR1 0 km + 33/33 + 67/67 — the migration (§8) is value-preserving on the SGP4 path |

astropy/ERFA is the IAU-standard implementation, independent of any in-repo code, so a match is *proof*, not
self-consistency. `ΔUT1 = 0` reproduces today's behavior exactly — the EPH invariant.

## 7. The abstraction that falls out (the functions, as consequences)

The theory above forces exactly this shape — nothing invented:

```cpp
enum class TimeScale { TAI, TT, TDB, UT1, UTC };

template<typename T>
struct Epoch {                       // one canonical instant
    math::TrackedValue<T> jd_day;    // integer/half-integer part  (two-part split, §3)
    math::TrackedValue<T> jd_frac;   // sub-day part, |·| < 0.5     (full precision)
    TimeScale scale;                 // which line we are on        (§1)
    // jd() = jd_day + jd_frac  (combine only when the magnitude loss is acceptable)
};

// Conversions — each a documented, oracle-verifiable function (§1 relationships):
Epoch<T> to_TT (const Epoch<T>&);                 // TAI+32.184 (exact); UTC via leap table; UT1 via ΔUT1
Epoch<T> to_UT1(const Epoch<T>&, TrackedValue<T> dUT1);   // ΔUT1=0 default → today's behavior
// TDB via the periodic generator (truncation→accuracy), mirroring obliquity_iau2006.

// Derived VIEWS (pure functions — the §5 reconciliation):
TrackedValue<T> centuries_since_j2000(const Epoch<T>&);   // (jd_TT − 2451545)/36525   → GMST, IAU polys
TrackedValue<T> days_since_1900     (const Epoch<T>&);    //  jd     − 2415020          → SR3 lunisolar
```

Generic over `T`; every value `TrackedValue`-tracked; conventions are `exact`/`defined`, ΔUT1/leap-seconds are
`measured` inputs, TDB is a generator. Consumers stop juggling bare JDs and named offsets; they ask the epoch
for the view they need.

## 8. Migration (value-preserving; SGP4 path stays OR1-frozen)

`elements.epoch_jd` (bare JD) → `Epoch` (scale = UTC, the SGP4 convention). `deep_space.h`'s
`epoch_jd − 2415020` → `days_since_1900(epoch)`; `sidereal_time.h`'s century → `centuries_since_j2000(epoch)`.
Each substitution is the *same number* → OR1 0 km, 33/33, 67/67 bit-exact. The new capability — explicit time
scales, an overridable ΔUT1, a TDB generator — becomes available **without touching the frozen authentic
path**. Implementation proceeds only after this note is accepted, bottom-up, each unit isolated and
oracle-gated (a new `test_time_scales` gate).

## References (born-digital)

- IERS Conventions (2010), Ch. 1 (time), Ch. 5 (rotation) — scale relationships, ΔUT1, leap seconds.
- IAU 2006 Resolutions — TT, TDB definitions; the J2000 epoch.
- Vallado, *Fundamentals of Astrodynamics and Applications* — time systems; the two-part JD; the in-repo
  `tle_parser` calendar→JD formula.
- Meeus, *Astronomical Algorithms*, Ch. 7 — Julian Date algorithm.
- astropy `astropy.time.Time` / ERFA (IAU SOFA) — the verification oracle.
