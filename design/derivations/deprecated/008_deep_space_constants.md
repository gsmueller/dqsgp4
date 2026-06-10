# Derivation 008: Deep Space Constants from Spacetrack Report #3

## References

- [SR3] Hoots, F.R. and Roehrich, R.L. (1980), "Spacetrack Report No. 3: Models for Propagation of NORAD Element Sets", Aerospace Defense Command, Peterson AFB. Available: https://celestrak.org/NORAD/documentation/spacetrk.pdf
- [Hujsak79] Hujsak, R.S. (1979), "A Restricted Four Body Solution for Resonating Satellites Without Drag", AIAA Paper No. 79-136, Aerospace Defense Command.
- [NASA-Eclipse] NASA Goddard, "Eclipses and the Moon's Orbit", https://eclipse.gsfc.nasa.gov/SEhelp/moonorbit.html
- [Wikipedia-LunarMonth] "Lunar month", Wikipedia, https://en.wikipedia.org/wiki/Lunar_month
- [NASA-ECI-SGP4] NASA ECI SGP4 implementation, https://github.com/nasa/ECI/blob/master/examples/sgp4Prop/fsw/src/SGP4.c

## Source

[SR3] Section 10 "The Deep-Space Subroutine", pages 58-68. The FORTRAN DATA statements on page 59 are the authoritative source of all deep-space constants.

## Solar/Lunar Mean Motions

### ZNS = 1.19459E-5 rad/min — Solar Mean Anomaly Rate

**Derivation:**
$$ZNS = \frac{2\pi}{T_\text{sidereal year} \times 1440 \text{ min/day}}$$

Computing from the sidereal year $T = 365.25636$ days:
$$ZNS = \frac{2\pi}{365.25636 \times 1440} = \frac{6.28318}{525969.2} = 1.19443 \times 10^{-5}$$

This doesn't match exactly. The period that gives ZNS = 1.19459E-5 exactly is:
$$T = \frac{2\pi}{1.19459 \times 10^{-5} \times 1440} = 365.257 \text{ days}$$

This is between the sidereal year (365.256 days) and the anomalistic year (365.260 days). The SGP4 value is likely a rounded value from circa 1970s astronomical constants. **Exact provenance unclear — treat as ModelValue.**

### ZNL = 1.5835218E-4 rad/min — Lunar Mean Anomaly Rate

**Derivation:**
$$ZNL = \frac{2\pi}{T_\text{anomalistic month} \times 1440 \text{ min/day}}$$

Computing from the anomalistic month $T = 27.554551$ days:
$$ZNL = \frac{2\pi}{27.554551 \times 1440} = \frac{6.28318}{39678.6} = 1.5835218 \times 10^{-4}$$

**Confirmed:** Matches to 12 significant figures ($3.3 \times 10^{-12}$ difference). ZNL advances the Moon's **mean anomaly**, which requires the **anomalistic** period (perigee-to-perigee), not the sidereal period.

### Lunar Longitude Rate: 0.22997150 rad/day

This appears in DPINIT (page 59): `C=4.7199672+.22997150*DAY`

Computing from the sidereal month $T = 27.321582$ days:
$$\frac{2\pi}{27.321582} = 0.229971 \text{ rad/day} = 13.1764°/\text{day}$$

**Confirmed:** This is the Moon's **sidereal** mean longitude rate (not anomalistic). Used for the Moon's mean longitude, which advances at the sidereal rate.

### Solar Mean Anomaly Rate: 0.017201977 rad/day

This appears in DPINIT: `ZMOS=6.2565837D0+.017201977D0*DAY`

$$0.017201977 \times \frac{180}{\pi} = 0.98560°/\text{day}$$

The tropical year gives $360°/365.2422 = 0.98565°$/day. The anomalistic year gives $360°/365.2596 = 0.98555°$/day. The SGP4 value 0.98560° is between these — it corresponds to a year of approximately 365.25 days (the Julian year). **The SGP4 solar mean anomaly rate uses the Julian year convention.**

## Obliquity Constants

### ZCOSIS = 0.91744867, ZSINIS = 0.39785416

These are $\cos(\varepsilon)$ and $\sin(\varepsilon)$ where $\varepsilon$ is the obliquity of the ecliptic.

$$\varepsilon = \arcsin(0.39785416) = 23.4441°$$

The IAU 1976 obliquity at J2000.0 is 23.4393°. The SGP4 value 23.4441° is the obliquity at approximately epoch 1970, which is correct for the 1970s-era constants used in Spacetrack Report #3.

## Lunar Inclination Constants

### ZCOSIL = 0.91375164 - 0.03568096*CTEM

The base value 0.91375164 corresponds to $\cos(23.96°)$. The Moon's orbital inclination to the ecliptic is 5.145°, so the Moon's inclination to the equator oscillates between $23.44° - 5.15° = 18.29°$ and $23.44° + 5.15° = 28.59°$. The value $\arccos(0.91375164) = 23.96°$ is approximately the mean of these extremes. The coefficient 0.03568096 modulates this with the 18.6-year nodal regression cycle (via CTEM = cos(XNODCE)).

## Resonance Constants

### Q22 = 1.7891679E-6, Q31 = 2.1460748E-6, Q33 = 2.2123015E-7

These are tesseral harmonic coupling coefficients from the Earth's gravity field. They are related to the sectoral and tesseral harmonics $J_{22}$, $J_{31}$, $J_{33}$ of the geopotential model used in the original SDP4 derivation. Their exact provenance traces to Hujsak (1979).

### ROOT22 = 1.7891679E-6, ROOT32 = 3.7393792E-7, etc.

These are square roots of products of the Q coefficients and are used in computing resonance perturbation strengths.

## Key Finding: Two Different Lunar Rates

The SGP4 deep-space model uses TWO different lunar angular rates:

| Rate | Value | Period | Used For |
|---|---|---|---|
| ZNL | 1.5835218E-4 rad/min | 27.5546 days (anomalistic) | Advancing lunar **mean anomaly** in DPPER |
| 0.22997150 | rad/day | 27.3216 days (sidereal) | Advancing lunar **mean longitude** in DPINIT |

This is physically correct: the mean anomaly advances at the anomalistic rate, while the mean longitude advances at the sidereal rate. Our `astronomy/` module must provide BOTH rates from the appropriate fundamental periods.

## Impact on Our Design

The `FundamentalConstants<T>` struct in `astronomy/solar_system.h` needs:
- `lunar_anomalistic_period_days` (for ZNL = mean anomaly rate)
- `lunar_sidereal_period_days` (for mean longitude rate)
- `solar_anomalistic_year_days` or `julian_year_days` (for ZMOS rate)
- `sidereal_year_days` (for ZNS = mean anomaly rate, though the SGP4 value doesn't exactly match either year definition)

The computation functions should be injectable (lambda or template parameter) so that different Almanac editions can provide different formulas for sidereal time, ephemeris rates, etc.

## Equations to Resolve

| Equation | Status | What must be shown |
|---|---|---|
| ZNS solar mean motion derivation | CITED | = 2π/(year×1440)? Which year definition? |
| ZNL lunar anomalistic rate | CITED | = 2π/(27.554551×1440) — source of 27.554551d? |
| Obliquity sin/cos (ZCOSIS, ZSINIS) | CITED | Epoch? IAU 1976 or 2006? |
| Solar perigee (ZCOSGS, ZSINGS) | CITED | Epoch-frozen value — what epoch? |
| Q₂₂ = 1.7891679e-6 | NOT DERIVED | Tesseral harmonic from J₂₂ gravity coefficient |
| Q₃₁ = 2.1460748e-6 | NOT DERIVED | From J₃₁ |
| Q₃₃ = 2.2123015e-7 | NOT DERIVED | From J₃₃ |
| ROOT22, ROOT32, ROOT44, ROOT52, ROOT54 | NOT DERIVED | Square roots of Q products |
| C1SS, C1L third-body secular coefficients | NOT DERIVED | From GM_sun, GM_moon, distances |
| Lunar geometry epoch values (4.5236020, 9.2422029e-4, etc.) | CITED | J1900 epoch — derive from fundamental periods |

## Source Documents Required

| Source | Location | Availability | What it provides |
|---|---|---|---|
| Spacetrack Report No. 3 pp. 58-59 | In repo | ✓ | FORTRAN DATA statements with all deep-space constants |
| Hujsak (1979) AIAA 79-136 | External | NOT in repo | Q₂₂/Q₃₁/Q₃₃ provenance |
| EGM2008 coefficients | wgs84_markdown_parts/ | ✓ | For recomputing Qₙₘ from modern gravity field |
| IERS Conventions (2010) | External | NOT in repo | Modern obliquity, precession |
