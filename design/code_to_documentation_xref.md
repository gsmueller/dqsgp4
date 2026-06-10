# Code-to-Documentation Cross-Reference

Maps every constant and formula in the reference SGP4 implementations to the source documents.

## 1. Gravity Model Constants

### 1.1 Earth Radius

| Code | Value | Document | Eq/Table | Notes |
|---|---|---|---|---|
| dnwrnr `kXKMPER` (Globals.h:29) | 6378.135 | NOT in [NGA] WGS84 | — | This is the WGS72 value |
| Vallado wgs72old (SGP4.cpp:2101) | 6378.135 | — | — | WGS72 old |
| Vallado wgs72 (SGP4.cpp:2112) | 6378.135 | — | — | WGS72 |
| Vallado wgs84 (SGP4.cpp:2123) | 6378.137 | [NGA] Table 3.1, Eq 3-2 | $a = 6378137.0$ m | **Match** (km vs m) |
| [NGA] WGS84 | 6378.137 km | Table 3.1 | Defining parameter | |
| [M80] GRS80 | 6378.137 km | Sec. 4 | Defining parameter | Same as WGS84 |

**Issue:** dnwrnr hardcodes the WGS72 value (6378.135) and has no WGS84 option. Vallado supports all three.

### 1.2 Gravitational Parameter (GM / mu)

| Code | Value (km³/s²) | Document | Eq/Table | Notes |
|---|---|---|---|---|
| dnwrnr `kMU` (Globals.h:28) | 398600.8 | — | — | WGS72 value |
| Vallado wgs72old (SGP4.cpp:2100) | 398600.79964 | — | — | WGS72 old (different from dnwrnr!) |
| Vallado wgs72 (SGP4.cpp:2111) | 398600.8 | — | — | WGS72 |
| Vallado wgs84 (SGP4.cpp:2122) | 398600.5 | [NGA] Table 3.4 | $GM_\text{GPSNAV}$ Eq 3-12 | NOT the refined WGS84 GM |
| [NGA] WGS84 refined GM | 398600.4418 | Table 3.1, Eq 3-4 | Defining parameter | |
| [NGA] WGS84 GPS nav GM | 398600.50 | Table 3.4, Eq 3-12 | $GM_\text{GPSNAV}$ | Original rounded value |
| [M80] GRS80 | 398600.5 | Sec. 4 | Defining parameter | Same as WGS84 GPSNAV |

**Issue:** SGP4's "wgs84" uses $GM_\text{GPSNAV} = 398600.5$, NOT the refined $GM = 398600.4418$. This is documented in [NGA] Sec 3.7.1 — the original value is retained for GPS receiver compatibility. The difference is $0.0582$ km³/s² which caused a 1.3 m radial bias (fixed in OCS orbit estimation but retained in broadcast).

### 1.3 J2, J3, J4 (Zonal Harmonics)

| Code | J2 | J3 | J4 | Model |
|---|---|---|---|---|
| dnwrnr (Globals.h:30-32) | 1.082616e-3 | -2.53881e-6 | -1.65597e-6 | WGS72 |
| Vallado wgs72old (2104-2106) | 0.001082616 | -0.00000253881 | -0.00000165597 | WGS72old |
| Vallado wgs72 (2115-2117) | 0.001082616 | -0.00000253881 | -0.00000165597 | WGS72 |
| Vallado wgs84 (2126-2128) | 0.00108262998905 | -0.00000253215306 | -0.00000161098761 | WGS84 |
| [NGA] Table 3.5 | 1.082629821313e-3 | — | — | $J_{2\text{geo}}$ (derived from $f$) |
| [NGA] Appendix B, Eq B-21 | 1.082629821313e-3 | — | — | Computed from $e^2$, $m$, $q_0$ |

**Issue:** The Vallado WGS84 J2 value (0.00108262998905) has 12 significant figures. The [NGA] geometric J2 (1.082629821313e-3) has 13 figures. These are DIFFERENT values because:
- Vallado's J2 comes from the SGP4 gravity model (historical)
- [NGA]'s $J_{2\text{geo}}$ is derived from $f$ via Appendix B formula B-21

### 1.4 Derived Constants (xke, ck2, ck4, etc.)

| Code | Formula | Document |
|---|---|---|
| dnwrnr `kXKE` (Globals.h:42) | $60 / \sqrt{r_e^3 / \mu}$ | Not in [NGA] — SGP4-specific |
| dnwrnr `kCK2` (Globals.h:43) | $0.5 \times J_2$ | Not in [NGA] — SGP4-specific |
| dnwrnr `kCK4` (Globals.h:44) | $-0.375 \times J_4$ | Not in [NGA] — SGP4-specific |
| dnwrnr `kQOMS2T` (Globals.h:54) | $((120 - 78)/r_e)^4$ | Not in [NGA] — SGP4-specific |
| dnwrnr `kS` (Globals.h:56) | $1 + 78/r_e$ | Not in [NGA] — SGP4-specific |
| dnwrnr `kA3OVK2` (Globals.h:81) | $-J_3 / (CK2)$ | Not in [NGA] — SGP4-specific |
| Vallado `j3oj2` (SGP4.cpp:2107) | $J_3 / J_2$ | Not in [NGA] — SGP4-specific |
| Vallado `tumin` | $1 / xke$ | Not in [NGA] — SGP4-specific |

**These SGP4-specific constants have NO documentation in [NGA], [M80], or [HM67].** They originate from the Spacetrack Report #3 (1980) and the AIAA 2006-6573 paper.

### 1.5 Pi

| Code | Value | Digits | Notes |
|---|---|---|---|
| dnwrnr `kPI` (Globals.h:57) | 3.14159265358979323846264338327950288419716939937510582 | 53 | High precision |
| Vallado `#define pi` (SGP4.cpp:66) | 3.14159265358979323846 | 21 | Standard double precision |
| [NGA] — | Not specified | — | Uses $\pi$ symbol |
| [M80] — | Not specified | — | Uses $\pi$ symbol |

**Issue:** Vallado uses a `#define` macro for pi — not type-safe, not arbitrary precision. dnwrnr uses a `const double` with 53 digits (but `double` only uses ~16). Our implementation uses `boost::math::constants::pi<T>()` which provides pi to the full precision of `T`.

## 2. Deep Space Constants

### 2.1 Solar/Lunar Orbital Parameters

| Code Variable | dnwrnr (SGP4.cc) | Vallado (SGP4.cpp) | Document | Notes |
|---|---|---|---|---|
| ZNS (solar node rate) | Line 695: 1.19459e-5 | Line 254/450: 1.19459e-5 | Not in [NGA] | Spacetrack Report #3 |
| ZES (solar eccentricity) | Line 697: 0.01675 | Line 255/450: 0.01675 | Not in [NGA] | |
| ZNL (lunar node rate) | Line 698: 1.5835218e-4 | Line 256/451: 1.5835218e-4 | Not in [NGA] | |
| ZEL (lunar eccentricity) | Line 700: 0.05490 | Line 257/451: 0.05490 | Not in [NGA] | |
| ZSINIS | Line 702: 0.39785416 | Line 454: 0.39785416 | Not in [NGA] | sin(23.4°) ≈ obliquity |
| ZCOSIS | Line 701: 0.91744867 | Line 455: 0.91744867 | Not in [NGA] | cos(23.4°) ≈ obliquity |
| ZSINGS | Line 703: -0.98088458 | Line 457: -0.98088458 | Not in [NGA] | |
| ZCOSGS | Line 704: 0.1945905 | Line 456: 0.1945905 | Not in [NGA] | |
| C1SS | Line 696: 2.9864797e-6 | Line 452: 2.9864797e-6 | Not in [NGA] | |
| C1L | Line 699: 4.7968065e-7 | Line 453: 4.7968065e-7 | Not in [NGA] | |

**All solar/lunar constants match between dnwrnr and Vallado.** None appear in the WGS84 standard or GRS80 — they are SGP4-specific from Spacetrack Report #3.

### 2.2 Resonance Constants

| Constant | Value | Both implementations match? | Document |
|---|---|---|---|
| Q22 | 1.7891679e-6 | Yes | Spacetrack #3 |
| Q31 | 2.1460748e-6 | Yes | Spacetrack #3 |
| Q33 | 2.2123015e-7 | Yes | Spacetrack #3 |
| ROOT22 | 1.7891679e-6 | Yes | Spacetrack #3 |
| ROOT32 | 3.7393792e-7 | Yes | Spacetrack #3 |
| ROOT44 | 7.3636953e-9 | Yes | Spacetrack #3 |
| ROOT52 | 1.1428639e-7 | Yes | Spacetrack #3 |
| ROOT54 | 2.1765803e-9 | Yes | Spacetrack #3 |

### 2.3 Integration and Resonance Detection

| Constant | Value | Both match? | Meaning |
|---|---|---|---|
| STEP | 720.0 min | Yes | Deep space integration step |
| STEP2 | 259200.0 | Yes | STEP² / 2 |
| rptim / kTHDT | 4.37526908801129966e-3 | Yes | Earth rotation rate (rad/min) = $\omega \times 60$ |
| Synch lower | 0.0034906585 rad/min | Yes | ~225 min period threshold |
| Synch upper | 0.0052359877 rad/min | Yes | ~150 min period threshold (note: 0.0052359877 ≈ π/600) |
| Half-day lower | 8.26e-3 rad/min | Yes | |
| Half-day upper | 9.24e-3 rad/min | Yes | |

### 2.4 Lunar Geometry Initialization

| Code (both implementations) | Value | Meaning |
|---|---|---|
| 4.5236020 | Initial lunar node | Epoch 1900 reference |
| 9.2422029e-4 | Lunar node rate | rad/day |
| 0.91375164 | cos(obliquity of lunar orbit) | |
| 0.03568096 | Perturbation of lunar cos(i) | |
| 0.089683511 | Lunar ascending node factor | |
| 5.8351514 | Solar mean anomaly epoch | |
| 0.0019443680 | Solar mean anomaly rate | rad/day |
| 4.7199672 | Lunar argument epoch | |
| 0.22997150 | Lunar argument rate | rad/day |
| 6.2565837 | Solar mean longitude epoch | |
| 0.017201977 | Solar mean longitude rate | rad/day |

**None of these appear in [NGA] or [M80].** They originate from simplified solar/lunar ephemerides in the SGP4 analytical theory.

## 3. Sidereal Time Constants

| Code | Value | Document |
|---|---|---|
| dnwrnr: computed via `DateTime::ToGreenwichSiderealTime()` | Complex formula | [NGA] App. A, Eq A-16 (ERA) |
| Vallado `thgr70` (initl:1271) | 1.7321343856509374 | Not in [NGA] — epoch 1970 GMST |
| Vallado `c1` (initl:1270) | 1.72027916940703639e-2 | Not in [NGA] — GMST rate |
| Vallado `fk5r` (initl:1272) | 5.07551419432269442e-15 | Not in [NGA] — FK5 correction |

**Issue:** The SGP4 sidereal time computation uses a simplified formula from the 1970 epoch, NOT the modern ERA formula from [NGA] Appendix A (Eq A-16). This is intentional — SGP4 must use the same sidereal time model that generated the TLE elements.

## 4. Element Recovery Formulas

| Formula | dnwrnr (OrbitalElements.cc) | Vallado (sgp4init) | Document |
|---|---|---|---|
| $a_1 = (k_e / n_0)^{2/3}$ | Line 43 | Line ~1430 | Spacetrack #3 |
| $\delta_1 = \frac{3}{2}\frac{CK2}{a_1^2}\frac{3\cos^2 i - 1}{\beta_0^3}$ | Line 50-51 | Line ~1440 | Spacetrack #3 |
| $a_0 = a_1(1 - \frac{1}{3}\delta_1 - \delta_1^2 - \frac{134}{81}\delta_1^3)$ | Line 52 | Line 1246-1247 | Spacetrack #3 |
| $n_0'' = n_0 / (1 + \delta_0)$ | Line 55 | Line ~1250 | Spacetrack #3 |

**The $\frac{134}{81}$ coefficient** appears in both implementations. This is a third-order correction to the Brouwer mean motion recovery. It does NOT appear in [NGA], [M80], or [HM67] — it's specific to the SGP4 analytical theory.

## 5. Kepler Equation Solver

| Parameter | dnwrnr | Vallado | Notes |
|---|---|---|---|
| Max iterations | 10 | 10 | Both same |
| Convergence tolerance | 1.0e-12 | 1.0e-12 | Both same |
| Damping threshold | 0.95 | 0.95 | Both same — limits Newton step to ±0.95 |
| Method | Newton-Raphson | Newton-Raphson | Standard |

## 6. Protection Thresholds

| Threshold | dnwrnr | Vallado | Purpose |
|---|---|---|---|
| Min eccentricity | 0.0 (exception) | 1.0e-6 (clamp) | **Different!** dnwrnr throws; Vallado clamps |
| Max eccentricity | 1.0 (exception) | 1 - 1.0e-6 (clamp) | |
| Negative eccentricity | < 0 (exception) | < -0.001 (error code) | **Different threshold!** |
| Inclination singularity | Not explicit | 1.5e-12 (sgp4fix) | Vallado has explicit fix |
| 180° inclination | Implicit | 5.2359877e-2 rad ≈ 3° (sgp4fix) | Vallado has explicit fix |
| Lyddane threshold | kPI (≈ 0.2 rad) | 0.2 rad | Both same |
| Decay detection | r < 1.0 Earth radii | mrt < 1.0 | Both same |
| Deep space period | ≥ 225 min | ≥ 225 min | Both same |

## 7. Constants NOT in Any Standard Document

The following constants appear in the code but have NO traceability to [NGA], [M80], [HM67], or [IERS10]. They originate from the SGP4 analytical theory (Spacetrack Reports #3 and #6, AIAA 2006-6573):

- All `Q22, Q31, Q33, ROOT22...ROOT54` resonance coefficients
- All solar/lunar orbital parameters (`ZNS, ZES, ZNL, ZEL, ZSINIS, ZCOSIS...`)
- All lunar geometry initialization constants (`4.5236020, 9.2422029e-4, 0.91375164...`)
- The `134/81` element recovery coefficient
- Deep space integration step (720 min)
- Resonance detection thresholds
- The `g` polynomial coefficients for half-day resonance (3.616, -13.247, 16.290, etc.)
- The `f` coefficients for resonance (0.75, 1.5, 1.875, 35.0, 39.3750, 9.84375, etc.)
- Sidereal time epoch constants (thgr70, c1, fk5r)

**These need to be documented with `MeasuredValue` or a new `ModelValue` category** — they come from analytical theory, not measurement or definition. Their "accuracy" is the accuracy of the SGP4 analytical model itself.

## 8. Discrepancies Between Implementations

| Item | dnwrnr | Vallado | Impact |
|---|---|---|---|
| Gravity model | WGS72 only (hardcoded) | WGS72old, WGS72, WGS84 (selectable) | dnwrnr can't do WGS84 |
| mu (wgs72old) | 398600.8 | 398600.79964 | 0.00036 km³/s² difference |
| xke (wgs72old) | computed | 0.0743669161 (hardcoded) | Vallado hardcodes; dnwrnr computes |
| Error handling | C++ exceptions | Error codes (satrec.error) | Architectural difference |
| AFSPC mode | Not supported | Supported ('a' mode) | dnwrnr is "improved" mode only |
| 180° inclination | No explicit handling | sgp4fix at 5.2359877e-2 rad | Vallado more robust |
| Min eccentricity | 0 (throws exception) | 1.0e-6 (clamped) | Vallado more permissive |
| Alpha-5 TLE | Not supported | Supported | Vallado handles modern TLEs |
