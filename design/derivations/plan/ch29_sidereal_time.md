# Draft Plan: Chapter 29 — Sidereal Time and Earth Rotation

**Part VIII: Reference Frames and Time** | Implementation file: `sidereal_time.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathrm{JD}$ | Julian Date | §29.1 |
| $\mathrm{JD}_0$ | Julian Date of J2000.0: 2451545.0 | §29.1 |
| $T_{\mathrm{UT1}}$ | Julian centuries of UT1 from J2000.0 | §29.2 |
| $\theta_{\mathrm{GMST}}$ | Greenwich Mean Sidereal Time (radians or degrees) | §29.3 |
| $\omega_E$ | Earth rotation rate: $2\pi / 86164.0905$ rad/s | §29.5 |
| $\Delta\mathrm{UT1}$ | UT1 − UTC correction (from IERS bulletins) | §29.4 |
| $\Delta\mathrm{AT}$ | TAI − UTC (leap seconds) | §29.4 |
| $r_{\mathrm{sid/sol}}$ | Sidereal-to-solar day ratio: $366.2422/365.2422$ | §29.5 |

---

## Objectives

1. Define Julian Date and Julian centuries from first principles; derive the conversion from calendar date.
2. Derive the GMST polynomial coefficients (Aoki et al. 1982) and give the physical meaning of each term.
3. Establish the sidereal/solar day ratio from the definitions of tropical year and sidereal day.
4. Develop the UTC ↔ UT1 ↔ GMST conversion chain with error classification at each step.
5. Provide a precise definition of Earth rotation rate $\omega_E$ and its role in coordinate transforms (Ch 30) and resonance (Ch 28).

## Section Structure

### §29.1 Julian Date

This section defines Julian Date from first principles and derives an exact integer-arithmetic algorithm for conversion from calendar date.

Stub:

1. **Definition 29.1.1** (Julian Date)**.** *A continuous count of days and fractional days elapsed since noon Universal Time on 1 January 4713 BC (Julian proleptic calendar), denoted $\mathrm{JD}$.*
2. **Definition 29.1.2** (J2000.0 epoch)**.** *The standard astronomical epoch J2000.0 corresponds to $\mathrm{JD} = 2\,451\,545.0$, equal to 2000 January 1, 12:00:00 TT.*
3. **Theorem 29.1.1** (calendar-to-JD exactness)**.** *The integer-arithmetic algorithm for $\mathrm{JD}$ from $(Y, M, D, h, m, s)$ via the Modified Julian Day Number is exact (zero rounding error) for all Gregorian calendar dates after 15 October 1582.*
   — *Proof approach: modular arithmetic analysis — verify the algorithm reduces to floor-division identities that map each Gregorian date to a unique integer day number; the fractional day part $h/24 + m/1440 + s/86400$ is representable exactly in IEEE 754 double for all $h, m, s$ in their standard ranges.*

[M.29.1] UTC time tags on TLEs have finite resolution; $\sigma_m \sim 0.001$ s typical.

**Example 29.1.1** (ISS TLE epoch to Julian Date)**.** *Date: 2024 January 15, 18:30:00 UTC. Input: $Y=2024$, $M=1$, $D=15$, $h=18$, $m=30$, $s=0$. Expected output: $\mathrm{JD} = 2\,460\,324.270\,833\,33$. Source: direct integer arithmetic per the stated algorithm; cross-checked against USNO Julian Date converter.*

### §29.2 Julian Centuries

This section defines Julian centuries of UT1 from J2000.0 and derives the scaling factor from the Julian year definition.

Stub:

1. **Definition 29.2.1** (Julian centuries of UT1)**.** *$T_{\mathrm{UT1}} = (\mathrm{JD}_{\mathrm{UT1}} - 2\,451\,545.0)/36525$, where the divisor $36525 = 365.25 \times 100$ is exact by the definition of the Julian year.*

The argument uses UT1, not UTC or TT; this distinction is essential because the GMST polynomial (§29.3) is defined as a function of UT1 only.

**Example 29.2.1** (Julian centuries for J2000.0)**.** *At the J2000.0 epoch, $\mathrm{JD}_{\mathrm{UT1}} = 2\,451\,545.0$, so $T_{\mathrm{UT1}} = 0.0$ exactly. Source: definition of J2000.0 epoch [Vallado 2013, §3.5].*

### §29.3 GMST Polynomial

This section derives the IAU 1982 GMST polynomial, establishes the physical meaning of each coefficient, and quantifies the precision and accuracy of Horner evaluation.

Stub:

1. **Definition 29.3.1** (GMST polynomial)**.** *The Greenwich Mean Sidereal Time is given by*

$$\theta_{\mathrm{GMST}} = 67310.54841 + \bigl(8\,640\,184.812\,866 + 0.093\,104\,T_{\mathrm{UT1}} - 6.2\times10^{-6}T_{\mathrm{UT1}}^2\bigr)T_{\mathrm{UT1}} \tag{29.1}$$

*(seconds of time; convert to radians via $2\pi/86400$). [Aoki et al. 1982, Eq. 14]*

   — *Proof approach: polynomial fit to IERS tabulated data — derive the physical meaning of each coefficient: $67310.548\,41$ s is the GMST at J2000.0 at 0h UT1; $8\,640\,184.812\,866$ s/cy is Earth's mean sidereal rotation rate in excess of one sidereal day per century; $0.093\,104$ s/cy² captures secular variation in rotation rate from tidal deceleration and precession; $-6.2\times10^{-6}$ s/cy³ is a higher-order empirical correction. Each coefficient must be verified against Aoki et al. (1982) — do not accept from secondary sources without cross-check.*

2. **Theorem 29.3.1** (Horner precision for GMST polynomial)**.** *Horner evaluation of Eq. (29.1) requires three fused multiply-add operations and accumulates a relative precision error no larger than $3\,\mathbf{u}$ at double precision, where $\mathbf{u} = 2^{-53}$ is the unit roundoff.*
   — *Proof approach: argument reduction via Cody-Waite extended precision — for $|T_{\mathrm{UT1}}| \leq 2$ (year range 1900–2100), the dominant term $8\,640\,184.8\, T_{\mathrm{UT1}}$ is at most $\sim 1.7 \times 10^7$; the constant term $67310.5$ loses no significance relative to the polynomial range. Cody-Waite two-part argument splitting is recommended when $|T_{\mathrm{UT1}}| > 2$. Apply Ch 4, Theorem 4.8.1 for the general Horner error bound.*

[A.29.1] The polynomial is valid to $\sim 0.1$ s over the period 1900–2100; outside this range extrapolation error grows. [P.29.1] Horner evaluation recommended to minimize rounding (Ch 4, §4.8).

**Example 29.3.1** (GMST at J2000.0)**.** *At J2000.0 ($T_{\mathrm{UT1}} = 0$): from Eq. (29.1) the constant term gives $\theta_{\mathrm{GMST}} = 67310.548\,41$ s $= 280.460\,618\,37°$. Expected output: $\theta_{\mathrm{GMST}} = 280.460\,618\,37°$. Source: IAU 1982 GMST formula [Aoki et al. 1982, Eq. 14]; cross-check: [Vallado 2013, Example 3-5].*

### §29.4 UTC, UT1, and GMST Conversion Chain

This section defines the UTC, UT1, TAI, and TT time scales and establishes the complete conversion chain required to compute GMST from a TLE epoch tag.

Stub:

1. **Definition 29.4.1** (UTC)**.** *Coordinated Universal Time is an atomic time scale maintained within $\pm 0.9$ s of UT1 via integer leap-second insertions; $\mathrm{UTC} = \mathrm{TAI} - \Delta\mathrm{AT}$, where $\Delta\mathrm{AT}$ [s] is the accumulated leap-second count.*
2. **Definition 29.4.2** (UT1)**.** *Universal Time UT1 is tied to Earth's instantaneous rotation angle; $\mathrm{UT1} = \mathrm{UTC} + \Delta\mathrm{UT1}$ where $\Delta\mathrm{UT1}$ is the IERS Bulletin A correction, $|\Delta\mathrm{UT1}| < 0.9$ s by construction.*
3. **Definition 29.4.3** (TT)**.** *Terrestrial Time satisfies $\mathrm{TT} = \mathrm{TAI} + 32.184$ s (exact by convention).*

The SGP4 use case: TLE epoch is in UTC; propagation requires $T_{\mathrm{UT1}}$; the correction $\Delta\mathrm{UT1}$ is at most $\pm 0.9$ s. [M.29.2] $\Delta\mathrm{UT1}$ is measured and published by IERS with uncertainty $\sim 0.02$ ms; Tier I.

**Example 29.4.1** (UTC to UT1 conversion)**.** *Given UTC epoch 2024 January 15, 18:30:00 UTC with $\Delta\mathrm{UT1} = -0.1082$ s (IERS Bulletin A, week of 2024-01-11): $\mathrm{UT1} = 18$:29:59.8918. Propagated $\mathrm{JD}_{\mathrm{UT1}} = 2\,460\,324.270\,831\,08$, giving $T_{\mathrm{UT1}} = 0.240\,055\,8$ cy. Source: IERS Bulletin A values are illustrative; the correction procedure is per [Vallado 2013, §3.5].*

### §29.5 Earth Rotation Rate

This section derives the Earth rotation rate $\omega_E$ from the sidereal/solar day ratio and establishes its role in coordinate transformations and resonance computations.

Stub:

1. **Definition 29.5.1** (Earth rotation rate)**.** *$\omega_E = 2\pi / T_{\mathrm{sidereal}}$ where $T_{\mathrm{sidereal}} = 86\,164.0905$ s is the mean sidereal day. Numerical value: $\omega_E = 7.292\,115 \times 10^{-5}$ rad/s.*
2. **Theorem 29.5.1** (sidereal/solar day ratio)**.** *The sidereal day and solar day satisfy $T_{\mathrm{solar}}/T_{\mathrm{sidereal}} = n_{\mathrm{rev}}/(n_{\mathrm{rev}}-1)$ where $n_{\mathrm{rev}} = 365.2422$ is the number of tropical revolutions per year.*
   — *Proof approach: modular arithmetic analysis — the Earth completes exactly one more rotation relative to the stars than relative to the Sun per tropical year; this geometric identity fixes the ratio without free parameters.*

Role: $\omega_E$ appears in the resonance condition $\dot{\psi}_{lmpq}$ (Ch 28, §ch28-tesseral-resonance) and the TEME → PEF Coriolis correction (Ch 30, §30.3). [M.29.3] $T_{\mathrm{sidereal}}$ is derived from the tropical year via the ratio; Tier I.

**Example 29.5.1** (Earth rotation rate)**.** *From $T_{\mathrm{sidereal}} = 86\,164.0905$ s: $\omega_E = 2\pi / 86\,164.0905 = 7.292\,115\,486\,7 \times 10^{-5}$ rad/s. Source: [NGA.STND.0036, 2014, Table 3.1]; matched-pair value for SGP4 is $\omega_E = 7.292\,115 \times 10^{-5}$ rad/s **[MP]** from WGS72 [DMA 1974].*

### §29.6 Implementation Notes

This section documents the modular reduction of $\theta_{\mathrm{GMST}}$, the Cody-Waite large-argument case, and the initialization boundary between precomputed constants and per-step evaluations.

Stub: Modular reduction of $\theta_{\mathrm{GMST}}$ to $[0, 2\pi)$ via the angle-arithmetic techniques of Ch 7 (§ch07-angle-arithmetic). The reduction is Lipschitz-1 so carries no additional precision error beyond one addition. Constants ($\mathrm{JD}_0$, polynomial coefficients, $\omega_E$) computed once at initialization per the precompute protocol (Ch 37, §ch37-precomputed-constants). GMST is a function of propagation time $t$, not a precomputed constant; it must be evaluated per step.

1. **Proposition 29.6.1** (per-step GMST evaluation)**.** *The sidereal angle $\theta_{\mathrm{GMST}}(t)$ must be evaluated at every propagation step; it is not a precomputable constant because it depends on $T_{\mathrm{UT1}}(t)$, which varies linearly with propagation time $t$.*

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 4, Theorem 4.8.1 | §29.3 | Horner evaluation for polynomial precision |
| Ch 7, angle arithmetic | §29.6 | Modular reduction of $\theta_{\mathrm{GMST}}$ to $[0, 2\pi)$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 28, tesseral resonance | §29.5 | $\omega_E$ in commensurability condition |
| Ch 30, coordinate transformations | §29.3, §29.5 | $\theta_{\mathrm{GMST}}$ and $\omega_E$ for TEME to PEF rotation |
| Ch 37, precomputed constants | §29.6 | Initialization boundary for constants |
| Appendix A | §29.5 | WGS72/WGS84 $\omega_E$ values |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.29.1] | M | §29.1 | UTC timestamp on TLE has finite resolution; propagated to $T_{\mathrm{UT1}}$ |
| [M.29.2] | M | §29.4 | $\Delta$UT1 from IERS; uncertainty ~0.02 ms |
| [M.29.3] | M | §29.5 | $\omega_E$ derived from measured tropical year; Tier I |
| [P.29.1] | P | §29.3 | GMST polynomial evaluated via Horner; cubic in $T_{\mathrm{UT1}}$ |
| [A.29.1] | A | §29.3 | GMST polynomial valid to ~0.1 s for 1900–2100; model error outside this range |
| [A.29.2] | A | §29.3 | The Aoki et al. (1982) IAU 1982 GMST polynomial is the SGP4 matched-pair formula; the IAU 2006 Earth Orientation Model (Capitaine et al. 2003, A\&A 412) supersedes it with a correction of up to 0.01 s for dates outside 1980–2020. For SGP4 matched-pair use, the 1982 formula must be used as-is; for standalone time conversion, the IAU 2006 formula should be preferred. |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 7 |
| Theorems | 3 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 1 |
| Examples | 6 |
| Error Notes | 6 |
| Equations | ~10 |
| Sections | 6 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §29.1 | Draft | Julian Date |
| §29.2 | Draft | Julian centuries |
| §29.3 | Draft | GMST polynomial |
| §29.4 | Draft | UTC, UT1, and GMST conversion chain |
| §29.5 | Draft | Earth rotation rate |
| §29.6 | Draft | Implementation notes |
