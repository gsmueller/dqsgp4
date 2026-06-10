# Cross-Reference: Geodetic Standards, Naming Conventions, and Class Design

## 1. Source Documents

| ID | Document | Citation | Scope |
|---|---|---|---|
| **[NGA]** | NGA.STND.0036\_1.0.0\_WGS84 (2014) | NGA Standard, 2014-07-08 | WGS 84 definition, all constants, datum transforms, EGM2008, WMM, ECI-ECEF |
| **[M80]** | Geodetic Reference System 1980 (Moritz) | Bulletin Geodesique, Vol. 54 No. 3, 1980 | GRS 80 defining parameters, closed-form derivation formulas, numerical values |
| **[HM67]** | Physical Geodesy (Heiskanen & Moritz) | W.H. Freeman, 1967 | Theoretical foundation — equipotential ellipsoid theory, sections 2-7 through 2-10 |
| **[DMA87]** | DMA TR 8350.2-A Supplement Part I | DMA, 1 December 1987 | Supplementary derivations for WGS 84 TR 8350.2 |
| **[IERS10]** | IERS Conventions 2010 (TN 36) | Petit & Luzum, 2010 | EOP, coordinate transforms, precession-nutation, time systems |

## 2. Defining Parameters Cross-Reference

### 2.1 Four Defining Parameters

| Parameter | Symbol | [NGA] WGS 84 | [M80] GRS 80 | Units | Notes |
|---|---|---|---|---|---|
| Semi-major axis | $a$ | 6378137.0 | 6378137 | m | Identical |
| Reciprocal flattening | $1/f$ | 298.257223563 | 298.257222101 (derived) | — | WGS84: defining. GRS80: derived from $J_2$ |
| Geocentric gravitational constant | $GM$ | $3.986004418 \times 10^{14}$ | $3.986005 \times 10^{14}$ | m$^3$/s$^2$ | WGS84 refined in 1994 |
| Dynamical form factor | $J_2$ | derived: $1.082629821313 \times 10^{-3}$ | $108263 \times 10^{-8}$ (defining) | — | WGS84: derived from $f$. GRS80: defining |
| Angular velocity | $\omega$ | $7.292115 \times 10^{-5}$ | $7.292115 \times 10^{-11} \times 10^6$ | rad/s | Identical |

**Key design note:** The `EquipotentialEllipsoid<T>` class must support two initialization paths:
- **GRS80 path:** ($a$, $GM$, $J_2$, $\omega$) → solve $e^2$ iteratively → derive $f$
- **WGS84 path:** ($a$, $1/f$, $GM$, $\omega$) → $e^2 = 2f - f^2$ directly → derive $J_2$

### 2.2 Special WGS 84 Parameters (no GRS80 equivalent)

| Parameter | Symbol | [NGA] Value | Source |
|---|---|---|---|
| GPS navigation GM | $GM_\text{GPSNAV}$ | $3.9860050 \times 10^{14}$ m$^3$/s$^2$ | [NGA] Eq. 3-12 |
| IAU angular velocity | $\omega'$ | $7.2921151467 \times 10^{-5}$ rad/s | [NGA] Eq. 3-13 |
| Precessing frame angular velocity | $\omega^*$ | $7.2921158553 \times 10^{-5} + 4.3 \times 10^{-15} T_U$ | [NGA] Eq. 3-16 |
| EGM2008 dynamic $\bar{C}_{2,0}$ | $\bar{C}_{2,0\text{dyn}}$ | $-4.84165143790815 \times 10^{-4}$ | [NGA] Eq. 3-6 |
| EGM2008 dynamic $\bar{C}_{2,2}$ | $\bar{C}_{2,2\text{dyn}}$ | $2.43938357328313 \times 10^{-6}$ | [NGA] Eq. 3-7 |

## 3. Derived Constants Cross-Reference

### 3.1 Geometric Constants

| Constant | Symbol | [NGA] Table 3.5 | [M80] Sec. 4 | Formula | Eq. Refs |
|---|---|---|---|---|---|
| Flattening | $f$ | $3.3528106647475 \times 10^{-3}$ | $0.00335281068118$ | $f = (a-b)/a$ | [NGA] B-1, [M80] Sec. 3 |
| Semi-minor axis | $b$ | $6356752.3142$ m | $6356752.3141$ m | $b = a\sqrt{1-e^2}$ | [NGA] B-2, [M80] Sec. 3, [HM67] p.66 |
| First eccentricity squared | $e^2$ | $6.694379990141 \times 10^{-3}$ | $0.00669438002290$ | $e^2 = 2f - f^2$ | [NGA] B-4, [M80] Sec. 3, [HM67] (2-90) |
| Second eccentricity squared | $e'^2$ | $6.739496742276 \times 10^{-3}$ | $0.00673949677548$ | $e'^2 = e^2/(1-e^2)$ | [NGA] B-7, [M80] Sec. 3, [HM67] (2-94) |
| Linear eccentricity | $E$ | $5.2185400842339 \times 10^5$ m | $521854.0097$ m | $E = ae$ | [NGA] B-10, [M80] Sec. 3 |
| Polar radius of curvature | $c$ or $R_p$ | $6399593.6258$ m | $6399593.6259$ m | $c = a^2/b = a/(1-f)$ | [NGA] B-13, [M80] Sec. 3 |
| Meridian quadrant | $Q$ | — | $10001965.7293$ m | Series in $e'^2$ | [M80] Sec. 3 |
| Mean radius (arithmetic) | $R_1$ | $6371008.7714$ m | $6371008.7714$ m | $R_1 = a(1-f/3)$ | [NGA] B-14, [M80] Sec. 3 |
| Radius equal area | $R_2$ | $6371007.1810$ m | $6371007.1810$ m | Wallis integral series in $e'^2$ | [NGA] B-15, [M80] Sec. 3 |
| Radius equal volume | $R_3$ | $6371000.7900$ m | $6371000.7900$ m | $R_3 = \sqrt[3]{a^2 b}$ | [NGA] B-17, [M80] Sec. 3 |
| Axis ratio | $AR$ | $9.96647189335 \times 10^{-1}$ | — | $AR = b/a$ | [NGA] B-11 |

### 3.2 Physical Constants

| Constant | Symbol | [NGA] Table 3.6 | [M80] Sec. 4 | Formula | Eq. Refs |
|---|---|---|---|---|---|
| Normal potential | $U_0$ | $6.26368517146 \times 10^{7}$ | $62636860.850$ | $(GM/E)\arctan e' + \omega^2 a^2/3$ | [NGA] B-23, [M80] Sec. 3, [HM67] (2-61) |
| Gravity formula const | $m$ | $3.449786506841 \times 10^{-3}$ | $0.00344978600308$ | $m = \omega^2 a^2 b/GM$ | [NGA] B-20, [M80] Sec. 3, [HM67] (2-70) |
| Intermediate $q_0$ | $q_0$ | $7.334625787083 \times 10^{-5}$ | — | $(1/2)[(1+3/e'^2)\arctan e' - 3/e']$ | [NGA] B-18, [M80] Sec. 3, [HM67] (2-58),(2-71) |
| Intermediate $q_0'$ | $q_0'$ | $2.688041300461 \times 10^{-3}$ | — | $3(1+1/e'^2)(1-(1/e')\arctan e') - 1$ | [NGA] B-19, [M80] Sec. 3, [HM67] (2-67) |
| Normal gravity (equator) | $\gamma_e$ | $9.7803253359$ | $9.7803267715$ | $(GM/ab)(1-m-me'q_0'/(6q_0))$ | [NGA] B-24, [M80] Sec. 3, [HM67] (2-73) |
| Normal gravity (pole) | $\gamma_p$ | $9.8321849379$ | $9.8321863685$ | $(GM/a^2)(1+me'q_0'/(3q_0))$ | [NGA] B-25, [M80] Sec. 3, [HM67] (2-74) |
| Somigliana constant | $k$ | $1.931852652458 \times 10^{-3}$ | $0.001931851353$ | $k = b\gamma_p/(a\gamma_e) - 1$ | [NGA] B-26, [M80] Sec. 3 |
| Gravity flattening | $f^*$ | — | $0.005302440112$ | $f^* = (\gamma_p - \gamma_e)/\gamma_e$ | [M80] Sec. 3, [HM67] (2-75) |
| Mean gravity | $\bar{\gamma}$ | $9.7976432223$ | $9.797644656$ | Series in $e^2, k$ | [NGA] B-27, [M80] Sec. 3 |
| Dynamical form factor | $J_2$ | $1.082629821313 \times 10^{-3}$ | $1.08263 \times 10^{-3}$ (defining) | $(e^2/3)(1-2me'/(15q_0))$ | [NGA] B-21, [M80] Sec. 3, [HM67] (2-90) |
| Normalized zonal | $\bar{C}_{2,0\text{geo}}$ | $-4.84166774985 \times 10^{-4}$ | — | $-J_2/\sqrt{5}$ | [NGA] B-22, [NGA] Eq. 3-1 |
| Mass of Earth | $M$ | $5.9721864 \times 10^{24}$ kg | — | $M = GM/G$ | [NGA] B-28 |

### 3.3 Even Zonal Harmonics ($J_{2n}$)

| Harmonic | [NGA] Table 3.5 $J_{2\text{geo}}$ | [M80] Sec. 4 | Formula |
|---|---|---|---|
| $J_2$ | $1.082629821313 \times 10^{-3}$ | $108263 \times 10^{-8}$ (defining) | See above |
| $J_4$ | — | $-0.00000237091222$ | $J_{2n} = (-1)^{n+1}\frac{3e^{2n}}{(2n+1)(2n+3)}(1-n+5nJ_2/e^2)$ |
| $J_6$ | — | $0.00000000608347$ | Same formula |
| $J_8$ | — | $-0.00000000001427$ | Same formula |

**Formula source:** [M80] Sec. 3, [HM67] (2-92), (2-92')

## 4. Formula Cross-Reference by Topic

### 4.1 Somigliana's Formula (Normal Gravity on Ellipsoid)

| Variant | Formula | Sources |
|---|---|---|
| Closed (original) | $\gamma = \frac{a\gamma_e\cos^2\Phi + b\gamma_p\sin^2\Phi}{\sqrt{a^2\cos^2\Phi + b^2\sin^2\Phi}}$ | [M80] Sec. 2, [HM67] p.70 |
| Closed (convenient) | $\gamma = \gamma_e\frac{1+k\sin^2\Phi}{\sqrt{1-e^2\sin^2\Phi}}$ | [M80] Sec. 3, [NGA] Eq. 4-1 |
| Series (abbreviated) | $\gamma = \gamma_e(1 + f^*\sin^2\Phi - \frac{1}{4}f_4\sin^2 2\Phi)$ | [M80] Sec. 3, [HM67] (2-115),(2-116) |
| Series (full) | $\gamma = \gamma_e(1 + \sum a_{2n}\sin^{2n}\Phi)$ | [M80] Sec. 3 |

### 4.2 Mean Gravity (Wallis-type Integrals)

The mean value of normal gravity over the ellipsoid surface is defined as the surface integral of $\gamma$ divided by the surface area. Using geodetic latitude $\Phi$ and the Moritz (1980) formulation, the formula from [M80] Sec. 3 is:

$$\bar{\gamma} = \gamma_e \left(1 + \frac{1}{6}e^2 + \frac{1}{3}k + \frac{59}{360}e^4 + \frac{5}{18}e^2k + \cdots \right)$$

> **Note:** The integral formulation for $\bar{\gamma}$ involves the ratio of two elliptic-type integrals whose exact form depends on the choice of integration variable and surface element. The series coefficients above are verified against both [M80] Sec. 3 and [NGA] Eq. B-27 and are the authoritative values. The intermediate integral derivation in `series_coefficient_derivation.md` Section 3.3 contains errors in the exponents and intermediate coefficients that do not affect the final verified result. A rigorous re-derivation from Heiskanen & Moritz (1967) Section 2-13 is needed before using the integral form directly.

**Series coefficients** (verified across [NGA] B-27 and [M80] Sec. 3):

| Term | Coefficient | Type |
|---|---|---|
| $e^2$ | $1/6$ | Geometric |
| $k$ | $1/3$ | Physical |
| $e^4$ | $59/360$ | Geometric |
| $e^2 k$ | $5/18$ | Mixed |
| $e^6$ | $2371/15120$ | Geometric |
| $e^4 k$ | $259/1080$ | Mixed |
| $e^8$ | $270229/1814400$ | Geometric |
| $e^6 k$ | $9623/45360$ | Mixed |

These coefficients arise from generalized binomial expansion of $(1-e^2\sin^2\Phi)^{-n}$ integrated with $\cos\Phi$ weighting (Wallis-type).

### 4.3 Sphere of Equal Area ($R_2$) Series

$$R_2 = c\left(1 - \frac{2}{3}e'^2 + \frac{26}{45}e'^4 - \frac{100}{189}e'^6 + \frac{7034}{14175}e'^8\right)$$

**Source:** [M80] Sec. 3, [NGA] B-15

Coefficients: $2/3$, $26/45$, $100/189$, $7034/14175$ (alternating sign series in $e'^{2n}$)

### 4.4 Meridian Quadrant ($Q$) Series

$$Q = c\frac{\pi}{2}\left(1 - \frac{3}{4}e'^2 + \frac{45}{64}e'^4 - \frac{175}{256}e'^6 + \frac{11025}{16384}e'^8\right)$$

**Source:** [M80] Sec. 3 only (not in [NGA])

Coefficients: $3/4$, $45/64$, $175/256$, $11025/16384$ — these are generalized binomial coefficients of $(1+x)^{-3/2}$.

### 4.5 Normal Gravity Above Ellipsoid

| Approach | Formula | Accuracy | Source |
|---|---|---|---|
| Taylor series | $\gamma_h = \gamma[1-(2/a)(1+f+m-2f\sin^2\varphi)h + (3/a^2)h^2]$ | Near-surface | [NGA] Eq. 4-3, [HM67] |
| Ellipsoidal coords (approx) | $\gamma_h \approx \sqrt{\gamma_\beta^2 + \gamma_u^2}$ | Sub-$\mu$gal to 20km | [NGA] Eq. 4-4 |
| Ellipsoidal coords (exact) | $\gamma_h = -\gamma_r\cos\alpha - \gamma_\psi\sin\alpha$ | Exact | [NGA] Eq. 4-16 |

### 4.6 Iterative $e^2$ from $J_2$ (GRS80 initialization path)

$$e^2 = 3J_2 + \frac{4}{15}\frac{\omega^2 a^3}{GM}\frac{a^3}{2q_0}$$

with $2q_0 = (1+3/e'^2)\arctan e' - 3/e'$

**Source:** [M80] Sec. 3 — This is the equation solved iteratively when $J_2$ is the defining parameter instead of $f$.

### 4.7 Clairaut's Theorem (Validation Check)

$$f + f^* = \frac{\omega^2 b}{\gamma_e}\left(1 + \frac{e'q_0}{2q_0}\right)$$

**Source:** [M80] Sec. 3, [HM67] (2-75) — Provides a closed-form validation that the computed $f$, $f^*$, $\gamma_e$ are self-consistent.

## 5. Proposed Class Hierarchy

```
math/
  constants.h          — pi<T>(), two_pi<T>(), etc.
  functions.h          — Mod<T>, WrapTwoPI<T>, DegreesToRadians<T>
  vector3.h            — Vector3<T>
  binomial.h           — GeneralizedBinomialCoefficient<T>, BinomialExpansion<T>
  wallis_integrals.h   — WallisCosineIntegral<T>

geodesy/
  equipotential_ellipsoid.h  — EquipotentialEllipsoid<T>
                                (central class: 4 defining params → all derived constants)
                                Methods: NormalGravity(), NormalGravityAbove()
                                Validation: Clairaut's theorem check
  wgs84.h                    — WGS84<T> : static instance of EquipotentialEllipsoid<T>
  grs80.h                    — GRS80<T> : static instance of EquipotentialEllipsoid<T>
  reference_ellipsoids.h     — Catalog of named ellipsoids (a, 1/f pairs from App. C)
  datum_transformation.h     — 3/7/10/14 param transforms

sgp4/
  gravity_model.h            — SGP4GravityModel<T> (wgs72old/wgs72/wgs84 variants)
  tle_parser.h/.cpp          — TLE parsing
  orbital_elements.h         — OrbitalElements<T>
  epoch.h                    — JulianDate, GMST
  state_vector.h             — StateVector<T>
  exceptions.h               — Exceptions
  sgp4_propagator.h          — SGP4Propagator<T>
```

## 6. Initialization Sequence (Compute Once)

```
1. Construct EquipotentialEllipsoid<T> with defining parameters
2. Compute e² (directly from f, or iteratively from J₂)
3. Compute b, e, e', E, c from e²
4. Compute q₀, q₀' (arctan-based, numerically sensitive)
5. Compute m = ω²a²b/GM
6. Compute γₑ, γₚ from closed formulas
7. Compute k = bγₚ/(aγₑ) - 1
8. Compute f* = (γₚ - γₑ)/γₑ
9. Validate: Clairaut's theorem
10. Compute U₀, γ̄, J₂ₙ (optional, on demand)
11. Compute R₁, R₂, R₃, Q (optional, on demand)
```

All values stored as members. No recomputation per call.
