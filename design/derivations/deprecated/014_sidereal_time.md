# Derivation 014: Greenwich Mean Sidereal Time (GMST)

## Purpose

Derive the GMST polynomial needed for the TEME→ECEF frame transformation.
GMST tells us the rotational orientation of the Earth at a given instant.

## Physical Definition

Greenwich Mean Sidereal Time is the hour angle of the **mean vernal equinox**
measured westward from the Greenwich meridian. It accumulates due to two effects:

1. **Earth's diurnal rotation** — one sidereal day ≈ 23h 56m 4.09s
2. **Precession of the equinox** — the equinox drifts westward at ~50.3"/year
   due to lunisolar torques on Earth's equatorial bulge (the same J₂ oblateness
   that causes satellite nodal regression)

## Step 1: Sidereal Day vs Solar Day

Earth rotates 360° in one sidereal day. But in that time, Earth also moves
~0.986° along its orbit, so it needs an extra ~3m 56s to face the Sun again
(one solar day).

- Solar day: 86400 s (by definition of SI second)
- Sidereal day: 86164.09054 s (measured)
- Earth rotation rate: ωₑ = 2π / 86164.09054 = 7.292115×10⁻⁵ rad/s

This ωₑ is the same value used in the EquipotentialEllipsoid.

## Step 2: GMST at the J2000.0 Epoch

The IAU defines GMST at the J2000.0 epoch (2000 January 1, 12h TT = JD 2451545.0)
from the Earth Rotation Angle (ERA):

The GMST at 0h UT1 on a given Julian date is computed via a polynomial in Julian
centuries T_UT1 from the J2000.0 epoch:

$$\theta_{GMST} = 67310.54841 + (876600h + 8640184.812866)T + 0.093104 T^2 - 6.2 \times 10^{-6} T^3 \tag{014.Eq.1}$$

where:
- T = (JD_UT1 - 2451545.0) / 36525.0 (Julian centuries from J2000.0)
- The result is in seconds of time
- Convert to radians: multiply by 2π/86400

**Origin of the coefficients:**

The constant 67310.54841 seconds is GMST at J2000.0 midnight (0h UT1 on
2000 Jan 1). This is a measured quantity from VLBI observations.

The linear term 8640184.812866 T represents the number of sidereal time
seconds that accumulate per Julian century. Since there are 36525 solar days
per Julian century, and each sidereal day is 86400 sidereal seconds:

86400 × 36525 × (sidereal/solar day ratio - 1) + 86400×36525
= 86400 × 36525 × (366.2422/365.2422)
≈ 8640184.8 seconds per century

The quadratic term 0.093104 T² accounts for the secular change in Earth's
rotation rate (precession of the equinox). Precession causes the equinox to
move, which changes how we measure sidereal time relative to the stars.

The cubic term −6.2×10⁻⁶ T³ accounts for the higher-order precession rate
change, negligible for most applications.

## Step 3: GMST at an Arbitrary UT1

For a Julian date that is NOT at 0h UT1, we add the Earth rotation since
midnight:

$$\theta_{GMST}(JD) = \theta_{GMST}(JD_0) + \omega_E \times (JD - JD_0) \times 86400 \tag{014.Eq.2}$$

where JD₀ is the Julian date at the preceding 0h UT1.

Alternatively, using the full polynomial with T computed from the actual JD:

$$T_{UT1} = (JD_{UT1} - 2451545.0) / 36525.0 \tag{014.Eq.3}$$

$$\theta = 67310.54841 + (876600 \times 3600 + 8640184.812866) \times T_{UT1} + 0.093104 \times T_{UT1}^2 - 6.2 \times 10^{-6} \times T_{UT1}^3 \tag{014.Eq.4}$$

The 876600h term represents ~100 years of sidereal rotations and is absorbed differently depending on the formulation. The expression from Aoki et al. (1982) and the IERS Conventions is:

$$GMST(0h UT1) = 24110.54841 + 8640184.812866 T_{UT1} + 0.093104 T_{UT1}^2 - 6.2 \times 10^{-6} T_{UT1}^3 \tag{014.Eq.5}$$

in seconds of sidereal time at 0h UT1, where T_UT1 is in Julian centuries.

For an arbitrary time within the day:

$$GMST = GMST(0h) + 1.00273790935 \times UT1 \tag{014.Eq.6}$$

where UT1 is in seconds since 0h UT1, and the factor 1.00273790935 = 366.2422/365.2422
converts solar time intervals to sidereal time intervals.

## Step 4: Conversion to Radians

To convert from seconds of sidereal time to radians:

$$\theta_{rad} = GMST_{sec} \times \frac{2\pi}{86400} \tag{014.Eq.7}$$

Then wrap to [0, 2π).

## Step 5: SGP4 Implementation

The SGP4 THETAG function computes GMST at a given epoch. The TLE epoch gives
a Julian date, from which T_UT1 is computed, and the polynomial evaluated.

For the SGP4 standard model (SR3 era), the polynomial coefficients from
Aoki et al. (1982) are used. For modern configurations (IAU 2006), the
ERA-based formula is available but gives the same result to ~0.01" over
the SGP4 prediction horizon.

## Numerical Verification

At J2000.0 epoch (JD 2451545.0, T=0):
- GMST = 24110.54841 seconds = 6.697374558 hours
- In radians: 24110.54841 × 2π/86400 = 1.7535... rad ≈ 100.46°

This is the sidereal time at Greenwich at noon on 2000 Jan 1 — the Sun
was approximately overhead at ~80°W longitude, consistent with GMST ≈ 100°
(the equinox was ~100° west of Greenwich at that moment).
