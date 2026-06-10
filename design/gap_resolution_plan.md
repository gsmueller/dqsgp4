# Gap Resolution Plan: Addressing Every Undocumented Constant

## Gaps Identified in Code Review (from `documentation_gaps.md`)

### Category A: Resolvable from First Principles

These constants have known mathematical derivations. We will re-derive them in arbitrary precision rather than using the hardcoded values.

#### A.1 Brouwer Secular Rate Polynomials

| Constants | Current Code | Resolution |
|---|---|---|
| 13, 78, 137 (xmdot J₂² term) | Hardcoded in SGP4 | Re-derive from Brouwer (1959) Hamiltonian perturbation theory |
| 7, 114, 395 (omgdot J₂² term) | Hardcoded | Same |
| 3, 36, 49 (omgdot J₄ term) | Hardcoded | Same |
| 4, 19 (xnodot J₂² term) | Hardcoded | Same |
| 3, 7 (xnodot J₄ term) | Hardcoded | Same |
| 134/81 (element recovery) | Hardcoded | Series reversion of Brouwer mean-to-osculating |

**Derivation approach:** The Brouwer theory expands the satellite Hamiltonian in powers of $J_2$ and solves via canonical transformations (von Zeipel method). The secular rates are the averaged Hamiltonian's partial derivatives with respect to the action variables. The polynomial coefficients in $\cos^2 i$ come from expanding Legendre polynomials $P_2(\cos\theta)$ and their powers.

**We have:** Brouwer (1959) paper downloaded. The equations are in Sections 3-6.

#### A.2 Kaula Inclination Functions $F_{lmp}(i)$

| Constants | Current Code Values | Resolution |
|---|---|---|
| 3/4, 3/2, 15/8, 15/16 | f220, f221, f321, f322 (exact rationals) | Re-derive from Kaula recursion |
| 35, 315/8, 315/32, 315/64, 105/16 | f441, f442, f522, f523 | Same |
| 945/32 | f542, f543 | Same |

**Derivation approach:** The Kaula inclination function is a double sum over indices $t$ and $s$:

$$F_{lmp}(i) = \sum_{t} \sum_{s} (-1)^t \frac{(2l-2t)!}{2^{2l} \, t! \, (l-t)! \, (l-m-2t)!} \binom{l-m-2t}{p-t-s} \binom{m+2t-2s}{s} (\sin i)^{l-m-2t+2s} (-\cos i)^{m+2t-2s-s'}$$

> **Note:** The exact form of the Kaula F-function varies by reference (Kaula 1966 vs. Allan 1965 vs. Giacaglia 1977). The key property is that for each specific $(l,m,p)$, the double sum collapses to a finite polynomial in $\sin i$ and $\cos i$ with exact rational coefficients. The implementation should use a verified recursion (e.g., the Gooding & King 1989 recursion) rather than the explicit double-sum formula, which is error-prone for hand computation. For each specific $(l,m,p)$, this reduces to an explicit polynomial in $\sin i$ and $\cos i$ with exact rational coefficients.

These are computable to arbitrary precision from the formula. The FORTRAN code has truncated decimal approximations of exact rationals (e.g., 4.92187512 instead of 315/64).

**Our approach:** Implement the Kaula recursion as a function `KaulaF<T>(l, m, p, sin_i, cos_i)` that computes the exact rational coefficients at initialization time.

#### A.3 Even Zonal Harmonics $J_{2n}$

Already have the closed formula from [M80]:

$$J_{2n} = (-1)^{n+1}\frac{3e^{2n}}{(2n+1)(2n+3)}\left(1 - n + 5n\frac{J_2}{e^2}\right)$$

Computable to arbitrary precision. No gap.

### Category B: Require Specific Physical Model Knowledge

These constants come from simplified physical models. We use them as `ModelValue` with provenance tracking.

#### B.1 Solar/Lunar Orbital Parameters

| Constant | Value | Physical Meaning | Resolution |
|---|---|---|---|
| ZNS = 1.19459e-5 | Solar mean motion | = $2\pi / (365.25 \times 1440)$ rad/min | **Derivable:** exact from the tropical year |
| ZES = 0.01675 | Solar eccentricity | Earth's orbital eccentricity ~1970 | **ModelValue:** epoch-dependent, ~0.0167 currently |
| ZNL = 1.5835218e-4 | Lunar mean motion | = $2\pi / (27.5546 \times 1440)$ rad/min | **Derivable:** from the sidereal month |
| ZEL = 0.05490 | Lunar eccentricity | Moon's orbital eccentricity | **ModelValue:** approximately constant |
| ZSINIS = 0.39785416 | $\sin(23.4441°)$ | Obliquity of the ecliptic | **Derivable:** $\sin(\varepsilon)$ from IAU value |
| ZCOSIS = 0.91744867 | $\cos(23.4441°)$ | Obliquity of the ecliptic | **Derivable:** $\cos(\varepsilon)$ from IAU value |
| ZSINGS = -0.98088458 | Solar argument of perigee | $\sin(\omega_\odot)$ at epoch | **ModelValue:** epoch-dependent |
| ZCOSGS = 0.1945905 | Solar argument of perigee | $\cos(\omega_\odot)$ at epoch | **ModelValue:** epoch-dependent |
| C1SS = 2.9864797e-6 | Solar secular coefficient | Derived from solar GM and distance | **Derivable** from fundamental constants |
| C1L = 4.7968065e-7 | Lunar secular coefficient | Derived from lunar GM and distance | **Derivable** from fundamental constants |

**Resolution strategy:** For the derivable ones, compute from fundamental astronomical constants (obliquity, orbital periods, eccentricities) at arbitrary precision. For epoch-dependent ones, either:
- Use the SGP4 frozen values (for compatibility), or
- Compute from modern ephemeris at the TLE epoch (for enhanced mode)

#### B.2 Lunar Geometry Initialization

| Constant | Value | Resolution |
|---|---|---|
| 4.5236020 | Lunar ascending node at J1900 | **Derivable** from fundamental lunar theory |
| 9.2422029e-4 | Lunar node regression rate | **Derivable** from lunar nodal period (18.613 years) |
| 0.91375164 | $\cos(i_\text{moon})$ base | **Derivable** from lunar inclination to ecliptic (5.145°) and obliquity |
| 0.03568096 | Nutation of lunar $\cos(i)$ | **Derivable** from lunar theory |
| 5.8351514, 0.0019443680 | Lunar perigee longitude | **Derivable** from lunar apsidal period (8.85 years) |
| 4.7199672, 0.22997150 | Lunar mean longitude | **Derivable** from sidereal month |
| 6.2565837, 0.017201977 | Solar mean anomaly | **Derivable** from tropical year |

**These are all derivable from fundamental astronomical periods.** The hardcoded values are simply these fundamental quantities evaluated at J1900 epoch. We can re-derive them from the defining orbital elements of the Moon and Sun.

### Category C: Opaque Polynomial Curve Fits (Hardest to Resolve)

#### C.1 Hansen/Kaula G-Function Polynomials (~100 coefficients)

| Function | Eccentricity Range | # Coefficients | Current Source |
|---|---|---|---|
| g201 | All | 2 (linear) | Hough curve fit |
| g211 | e ≤ 0.65 / e > 0.65 | 4 + 4 | Hough curve fit |
| g310 | e ≤ 0.65 / e > 0.65 | 4 + 4 | Hough curve fit |
| g322 | e ≤ 0.65 / e > 0.65 | 4 + 4 | Hough curve fit |
| g410 | e ≤ 0.65 / e > 0.65 | 4 + 4 | Hough curve fit |
| g422 | e ≤ 0.65 / e > 0.65 | 4 + 4 | Hough curve fit |
| g520 | e ≤ 0.65 / e ≤ 0.715 / e > 0.715 | 4 + 4 + 4 | Hough curve fit |
| g533, g521, g532 | e < 0.7 / e ≥ 0.7 | 4 + 4 each × 3 | Hough curve fit |
| g200, g310, g300 | All (sync) | 3 + 2 + 3 | Hough curve fit |

**Total: ~100 polynomial coefficients that are curve fits to exact integrals.**

**The exact Hansen coefficient $G_{lpq}(e)$ is:**

$$G_{lpq}(e) = \frac{1}{2\pi}\int_0^{2\pi}\left(\frac{r}{a}\right)^{-(l+1)} e^{i(l-2p+q)v} e^{-i(l-2p)M}\,dM$$

where $v$ is the true anomaly as a function of $M$ (mean anomaly) and $e$.

This integral can be evaluated as a series in $e$ using Bessel function expansions:

$$G_{lpq}(e) = \sum_{s} \text{(rational coefficient)} \times e^{|s|}$$

The series converges for $e < 1$ and can be computed to arbitrary precision.

**Resolution:** Implement the Hansen coefficient as a function `HansenG<T>(l, p, q, e)` that evaluates the Bessel-function series to the required precision. This replaces all ~100 polynomial coefficients with a single function that computes the exact values. The polynomial fits become validation targets.

#### C.2 Tesseral Resonance Constants

| Constant | Value | Resolution |
|---|---|---|
| Q22 = 1.7891679e-6 | Tesseral harmonic coupling | **Derivable** from $J_{2,2}$ (sectoral harmonic of Earth's gravity field) |
| Q31 = 2.1460748e-6 | | **Derivable** from $J_{3,1}$ |
| Q33 = 2.2123015e-7 | | **Derivable** from $J_{3,3}$ |
| ROOT22...ROOT54 | Square roots of products | **Derivable** from the Q values |

These relate to specific tesseral harmonics of the Earth's gravity field. They can be recomputed from EGM2008 coefficients (which we have in Table 5.1 / Table 6.1 of the NGA document).

$$Q_{nm} \propto \sqrt{C_{nm}^2 + S_{nm}^2}$$

**Resolution:** Compute from the EGM2008 $\bar{C}_{nm}$ and $\bar{S}_{nm}$ coefficients with proper denormalization.

#### C.3 Resonance Phase Angles

| Constant | Value | Resolution |
|---|---|---|
| G22 = 5.7686396 | Phase of J₂₂ tesseral | **Derivable:** $\lambda_{22} = \arctan(S_{22}/C_{22})$ from EGM |
| G32 = 0.95240898 | Phase of J₃₂ | Same approach |
| G44 = 1.8014998 | Phase of J₄₄ | Same |
| G52 = 1.0508330 | Phase of J₅₂ | Same |
| G54 = 4.4108898 | Phase of J₅₄ | Same |
| FASX2 = 0.13130908 | Synchronous phase 2 | Related to above |
| FASX4 = 2.8843198 | Synchronous phase 4 | Related to above |
| FASX6 = 0.37448087 | Synchronous phase 6 | Related to above |

**Resolution:** Compute from EGM coefficients. The phase angle of tesseral harmonic $(n,m)$ is:

$$\lambda_{nm} = \frac{1}{m}\arctan\frac{S_{nm}}{C_{nm}}$$

The FASX values are related to these phases at specific resonant combinations.

## Implementation Priority

| Priority | Category | Action | Effort |
|---|---|---|---|
| **P0** | A.3 | $J_{2n}$ from closed formula | Already designed |
| **P0** | A.1 | Brouwer polynomials — verify the exact rational coefficients | Read Brouwer (1959) |
| **P1** | A.2 | Kaula $F_{lmp}$ — implement recursion formula | Moderate |
| **P1** | B.1 | Solar/lunar constants — derive from fundamental astronomy | Moderate |
| **P1** | C.2, C.3 | Tesseral resonance — compute from EGM2008 | Moderate |
| **P2** | B.2 | Lunar geometry — derive from fundamental lunar theory | Significant |
| **P3** | C.1 | Hansen $G_{lpq}$ — implement Bessel series evaluation | Significant (replaces ~100 magic numbers) |

P0 is needed before any implementation. P1 can be done during implementation. P2 and P3 are enhancement phases that can use the hardcoded values initially with `ModelValue` tracking.
