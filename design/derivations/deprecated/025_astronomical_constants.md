# Derivation 025: Astronomical Constants and Epoch Values

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Trace every astronomical constant used in the deep-space model to its physical origin.
For each constant in the SR3 FORTRAN DATA statements (ZNS, ZES, ZNL, ZEL, ZCOSIS, ZSINIS,
ZCOSGS, ZSINGS, C1SS, C1L, and the lunar geometry initialization values), determine:

1. What physical quantity it represents
2. What epoch it was evaluated at
3. What its modern value would be
4. Whether it should be a computed constant or a frozen `ModelValue` for SGP4 compatibility

**Code:** `src/astronomy/solar_system.h`, `src/ephemeris/celestial_body.h`

---

## Source Documents Required

| Source | Location | Availability | What it provides |
|--------|----------|:------------:|-----------------|
| Spacetrack Report No. 3 pp. 58-59 | `Spacetrack_Report_No3...pdf` | ✓ | FORTRAN DATA statements with all constants |
| IERS Conventions (2010) | — | NOT in repo | Modern values for obliquity, eccentricity, periods — needed for comparison |
| Meeus (1998) "Astronomical Algorithms" | — | NOT in repo | Fundamental astronomical periods, lunar/solar orbital elements |
| Explanatory Supplement to the Astronomical Almanac | — | NOT in repo | Definitive reference for astronomical constants and their epochs |

---

## Items to Resolve

| Constant | Value | Status | What must be determined |
|----------|-------|--------|------------------------|
| ZNS | $1.19459\times10^{-5}$ | TEMPLATE | Is this $2\pi/(365.25\times1440)$? What year definition — tropical, sidereal, or anomalistic? Derive from fundamental period |
| ZES | $0.01675$ | TEMPLATE | Earth's orbital eccentricity — at what epoch? Compare with modern value ($\approx 0.0167086$ at J2000). Rate of change $\approx -4\times10^{-5}$/century |
| ZNL | $1.5835218\times10^{-4}$ | TEMPLATE | Anomalistic month $\approx 27.554551$ d — derive from source. Compare with modern value. Why anomalistic (not sidereal or synodic)? |
| ZEL | $0.05490$ | TEMPLATE | Lunar eccentricity — constant or epoch-dependent? Modern value $\approx 0.0549$. Secular variation is small but nonzero |
| ZCOSIS | $0.91744867$ | TEMPLATE | $\cos\epsilon$ where $\epsilon$ is obliquity — what epoch? IAU 1976 ($\epsilon_0 = 23°26'21.448''$) or IAU 2006? Compute $\cos(23.4393°)$ and compare |
| ZSINIS | $0.39785416$ | TEMPLATE | $\sin\epsilon$ — must be consistent with ZCOSIS. Verify $\sin^2 + \cos^2 = 1$ to full precision |
| ZCOSGS | value from DATA | TEMPLATE | $\cos\bar{g}_\odot$ — solar argument of perigee at what epoch? Is this precessing or frozen? |
| ZSINGS | value from DATA | TEMPLATE | $\sin\bar{g}_\odot$ — must be consistent with ZCOSGS. Verify Pythagorean identity |
| C1SS | $2.9864797\times10^{-6}$ | TEMPLATE | Third-body secular coefficient for Sun — derive from $GM_\odot$, $a_\odot$, and the disturbing function expansion. Show the physical formula |
| C1L | $4.7968065\times10^{-7}$ | TEMPLATE | Third-body secular coefficient for Moon — derive from $GM_☾$, $a_☾$, and the disturbing function expansion. Show ratio C1L/C1SS and relate to mass/distance ratios |
| Lunar geometry: $4.5236020$ | value from DATA | TEMPLATE | J1900 epoch lunar mean longitude or argument — derive from fundamental periods and J1900.0 reference |
| Lunar geometry: $9.2422029\times10^{-4}$ | value from DATA | TEMPLATE | Lunar rate in rad/min — derive from sidereal month and verify against ZNL |
| $0.22997150$ rad/day | value from DATA | TEMPLATE | Lunar sidereal rate vs. ZNL anomalistic rate — why two different lunar rates? (Referenced in 008 but needs formal derivation showing sidereal vs. anomalistic distinction) |

---

## Assessment

The FORTRAN DATA constants are frozen values from approximately the 1960s-1970s. Most can be
traced to well-known astronomical quantities, but their exact epoch and source convention
(IAU 1976, pre-IAU, etc.) must be determined. Three of the four required references are not
in the local repository, but the [SR3] DATA statements provide the values to be analyzed,
and many can be reverse-engineered from known astronomical periods and constants.

The key question for each constant is: **compute or freeze?** For SGP4 TLE compatibility,
all must be frozen at their original values. But documenting the modern equivalents enables
future high-precision implementations to use updated values when TLE compatibility is not required.
