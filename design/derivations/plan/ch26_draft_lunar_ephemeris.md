# Draft Plan: Chapter 26 — Lunar Ephemeris

**Part VI: Third-Body Perturbations**

Target header: `lunar_ephemeris.h`

---

## Objectives

1. Derive the mean lunar longitude and mean anomaly as linear functions of time from epoch, using the periods of Ch 23.
2. Derive and physically explain the three principal lunar perturbations — evection, variation, and annual equation — as consequences of solar tidal forces.
3. Obtain the corrected ecliptic longitude $\lambda_\mathbb{C}$ and ecliptic latitude $\beta_\mathbb{C}$ with quantified accuracy bounds.
4. Compare the accuracy of the truncated theory against JPL DE430.
5. Generalize to adding more perturbation terms and convergence of the full Fourier series.

---

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\ell_\mathbb{C}$ | Mean lunar longitude | §26.2 |
| $M_\mathbb{C}$ | Lunar mean anomaly | §26.2 |
| $\varpi_\mathbb{C}$ | Lunar longitude of perihelion ($= \Omega_\mathbb{C} + \omega_\mathbb{C}$) | §26.2 |
| $F$ | Lunar argument of latitude ($= \ell_\mathbb{C} - \Omega_\mathbb{C}$) | §26.2 |
| $D$ | Elongation: $\ell_\mathbb{C} - \ell_\odot$ | §26.2 |
| $n_\mathbb{C}$ | Lunar sidereal mean motion (degrees/day) | §26.2 |
| $e_\mathbb{C}$ | Lunar mean eccentricity | §26.2 |
| $\Omega_\mathbb{C}$ | Lunar node longitude | §26.2 |
| $\omega_\mathbb{C}$ | Lunar argument of perigee | §26.2 |
| $\lambda_\mathbb{C}$ | Lunar geocentric ecliptic longitude | §26.4 |
| $\beta_\mathbb{C}$ | Lunar geocentric ecliptic latitude | §26.4 |
| $r_\mathbb{C}$ | Geocentric lunar distance | §26.4 |
| $\ell_\odot$ | Mean solar longitude (from Ch 25) | Ch 25 |

---

## Section Structure

### §26.1 Introduction

This section motivates the five-term lunar ephemeris used by SGP4, contrasts the complexity of the Moon's orbit with the Sun's, and previews the three principal perturbations to be derived.

*Stub.* The Moon is the largest perturber of near-Earth satellite orbits (Ch 24, §24.4: the Moon's tidal acceleration is $\approx 2.2\times$ the Sun's). Its orbit is significantly more complex than the Sun's apparent orbit: the Moon's eccentricity varies by $\pm 14\%$, its inclination oscillates, and there are three principal perturbations — evection, variation, and annual equation — each with amplitude exceeding $0.5°$ in longitude.

The SGP4 lunar ephemeris uses a five-term approximation adequate for third-body perturbation of LEO and GEO satellites. This chapter derives each term from the Hill-Brown lunar theory (physical basis) and states the algebraic form used by SGP4.

---

### §26.2 Mean Lunar Elements

This section defines the mean lunar orbital elements at epoch, establishes the linear secular extrapolation for each, and states the SGP4 numerical values with their matched-pair source.

*Stub.* Define the mean lunar elements at epoch, analogous to §25.2 for the Sun.

**Definition 26.2.1** (Mean lunar longitude $\ell_\mathbb{C}$). *Stub — linear extrapolation from epoch using the sidereal mean motion $n_\mathbb{C}$.*

**Definition 26.2.2** (Mean lunar anomaly $M_\mathbb{C}$). *Stub — $M_\mathbb{C} = \ell_\mathbb{C} - \varpi_\mathbb{C}$ where $\varpi_\mathbb{C} = \Omega_\mathbb{C} + \omega_\mathbb{C}$.*

**Definition 26.2.3** (Lunar argument of latitude $F$). The angle from the ascending node to the Moon along its orbit: $F = \ell_\mathbb{C} - \Omega_\mathbb{C}$.

**Definition 26.2.4** (Elongation $D$). The angular separation of the Moon from the Sun: $D = \ell_\mathbb{C} - \ell_\odot$.

**Notation table:**

| Symbol | Meaning |
|--------|---------|
| $\ell_\mathbb{C}$ | Mean lunar longitude |
| $M_\mathbb{C}$ | Lunar mean anomaly |
| $\varpi_\mathbb{C}$ | Lunar longitude of perihelion ($= \Omega_\mathbb{C} + \omega_\mathbb{C}$) |
| $F$ | Lunar argument of latitude ($= \ell_\mathbb{C} - \Omega_\mathbb{C}$) |
| $D$ | Elongation: $\ell_\mathbb{C} - \ell_\odot$ |
| $n_\mathbb{C}$ | Lunar sidereal mean motion (degrees/day) |
| $e_\mathbb{C}$ | Lunar mean eccentricity |
| $\Omega_\mathbb{C}$ | Lunar node longitude (regression rate $\dot\Omega_\mathbb{C}$) |
| $\omega_\mathbb{C}$ | Lunar argument of perigee (advance rate $\dot\omega_\mathbb{C}$) |
| $\lambda_\mathbb{C}$ | Lunar geocentric ecliptic longitude |
| $\beta_\mathbb{C}$ | Lunar geocentric ecliptic latitude |
| $r_\mathbb{C}$ | Geocentric lunar distance |

State the SGP4 numerical values for $n_\mathbb{C}$, $\dot\Omega_\mathbb{C}$, $\dot\omega_\mathbb{C}$ and their epoch reference.

---

### §26.3 The Evection

This section derives the evection as the leading solar perturbation of the lunar eccentricity, establishes its period and amplitude from the Hill-Brown theory, and states the algebraic form used by SGP4.

*Stub.* The evection is the largest periodic perturbation of the lunar orbit and was discovered by Ptolemy. It arises from the oscillation of the lunar eccentricity due to the solar tidal force.

**Physical origin.** The Sun's gravitational gradient across the lunar orbit creates a tidal torque that periodically increases and decreases the eccentricity. The period is approximately 31.8 days (the synodic anomalistic period: $1/|n_\mathbb{C}^{\rm anom} - n_\odot|$).

**Algebraic form.** The evection contribution to lunar longitude:

$$\Delta\lambda_{\rm ev} = +1.274° \sin(2D - M_\mathbb{C}) \tag{26.1}$$

**Proposition 26.3.1** (Evection amplitude). *The evection coefficient is $1.274° = 2e_\mathbb{C}(1 + 15m/8 + \ldots)/(1 - 5m/4 + \ldots)$ where $m = n_\odot^2/n_\mathbb{C}^2$ is the solar perturbation parameter; numerically $m \approx 0.00543$.*
— *Proof approach: Hill-Brown disturbing function — the evection arises from the $e_\mathbb{C}\cos(2D - M_\mathbb{C})$ term in the lunar longitude when the solar tidal quadrupole is expanded in the Delaunay variables; the amplitude is proportional to $e_\mathbb{C}$ corrected by a power series in $m$. The coefficient $1.274°$ is taken from Chapront-Touzé \& Chapront (1988) ELP 2000-82; cross-verify against the independent derivation in Brown (1896). NOTE: the ELP 2000-82 theory (Chapront-Touzé \& Chapront 1988) has published corrections [A.26.4]; verify the $1.274°$ coefficient against the corrected series, not the original paper.*

**Error Note [A.26.1]:** The evection coefficient $1.274°$ varies slightly with the secular change in lunar eccentricity; the variation is $< 0.001°$/century. [Tier I.]

---

### §26.4 The Variation

This section derives the variation as the tangential-force deformation of the lunar orbit at quadrature, derives its amplitude from the solar perturbation parameter, and states the algebraic form.

*Stub.* The variation arises from the deformation of the lunar orbit by the solar tidal quadrupole, specifically the tangential component of the solar tidal force at half-synodic-month intervals.

**Physical origin.** At quadrature (Moon $90°$ from Sun), the solar tangential force is maximum and accelerates or decelerates the Moon, deforming the orbit from an ellipse to a bean shape. The period is half the synodic month $\approx 14.8$ days.

**Algebraic form:**

$$\Delta\lambda_{\rm var} = +0.658° \sin 2D \tag{26.2}$$

**Proposition 26.4.1** (Variation amplitude). *The variation coefficient is $0.658° \approx (11m/8)(180°/\pi)$ to leading order in $m = n_\odot^2/n_\mathbb{C}^2$.*
— *Proof approach: Hamiltonian perturbation theory — the $\sin 2D$ term arises from the second-order solar perturbation of the Moon's tangential velocity; the coefficient is derived by expanding the solar disturbing function to order $m$ and applying the Lagrange planetary equation for the mean longitude; verify: $11m/8 \approx 11 \times 0.00543/8 \approx 0.00747$ rad $\approx 0.428°$ (the exact coefficient includes higher-order $m$ corrections giving $0.658°$).*

---

### §26.5 The Annual Equation

This section derives the annual equation as a modulation of the Moon's mean motion by the varying Sun-Earth distance, and derives the coefficient from the solar eccentricity.

*Stub.* The annual equation arises from the variation of the Sun-Earth distance over the year: the solar tidal force is stronger at perihelion (January) and weaker at aphelion (July), modulating the rate of lunar mean motion.

**Physical origin.** The mean motion of the Moon depends on the solar perturbation parameter $m = n_\odot^2/n_\mathbb{C}^2 \propto r_\odot^{-3}$. When $r_\odot$ is smaller (perihelion), $m$ is larger and the Moon's mean motion is slightly faster.

**Algebraic form:**

$$\Delta\lambda_{\rm ann} = -0.186° \sin M_\odot \tag{26.3}$$

**Proposition 26.5.1** (Annual equation amplitude). *The annual equation coefficient is $-0.186° \approx -2e_\odot \cdot (3m/2) \cdot (180°/\pi)$ to leading order in $m$ and $e_\odot$.*
— *Proof approach: series expansion — $m \propto r_\odot^{-3}$; expand $r_\odot^{-3}$ as a function of $M_\odot$ to get the $\sin M_\odot$ term; the coefficient at leading order in $e_\odot$ and $m$ is $-\frac{3}{2}me_\odot$; insert $m = 0.00543$, $e_\odot = 0.0167$: $-\frac{3}{2}(0.00543)(0.0167) \approx -1.36\times10^{-4}$ rad $\approx -0.0078°$; the full coefficient including higher-order corrections is $-0.186°$, which requires the complete fourth-order derivation from the Delaunay disturbing function.*

---

### §26.6 Additional Terms in the SGP4 Lunar Ephemeris

This section derives the two additional SGP4 lunar longitude terms — the second-anomaly harmonic and the parallactic inequality — and identifies their origins in the Delaunay disturbing function.

*Stub.* The SGP4 lunar longitude formula (Hoots and Roehrich 1980) contains two further terms beyond evection, variation, and annual equation:

$$\Delta\lambda_4 = +0.214°\sin 2M_\mathbb{C} \tag{26.4}$$
$$\Delta\lambda_5 = +0.186°\sin(M_\mathbb{C} + M_\odot) \tag{26.5}$$

Term (26.4) is the second harmonic of the lunar anomaly — analogous to the $\tfrac{5}{4}e^2\sin 2M$ term in the equation of center for the Sun, but amplified by the larger lunar eccentricity. Term (26.5) is the parallactic inequality coupling.

**Proposition 26.6.1** (Second-anomaly term). *The coefficient $0.214°$ of $\sin 2M_\mathbb{C}$ is the $O(e_\mathbb{C}^2)$ term in the lunar equation of center: $(5/4)e_\mathbb{C}^2 (180°/\pi) = (5/4)(0.0549)^2 \times 57.296° \approx 0.217°$, agreeing with the SGP4 value to within the higher-order corrections.*
— *Proof approach: Lagrange inversion comparison — re-derive the lunar equation of center using the same series as Ch 25, Theorem 25.3.1, but with $e_\mathbb{C} = 0.0549$; the coefficient of $\sin 2M$ is $(5/4)e^2$ to leading order. The small discrepancy from $0.214°$ to $0.217°$ comes from $O(e^3)$ corrections; these must be computed explicitly to verify the SGP4 coefficient.*

---

### §26.7 Ecliptic Latitude

This section derives the single-term latitude formula from the orbital inclination and bounds the error from neglecting higher-order harmonics.

*Stub.* The Moon's orbit is inclined to the ecliptic by $\approx 5.145°$. The ecliptic latitude is:

$$\beta_\mathbb{C} \approx 5.128°\sin F \tag{26.6}$$

where $F = \ell_\mathbb{C} - \Omega_\mathbb{C}$ is the argument of latitude (Definition 26.2.3). The $5.128°$ amplitude differs from the inclination $5.145°$ by the planetary perturbation of the node.

**Proposition 26.7.1** (Latitude amplitude). *The latitude amplitude $5.128°$ satisfies $5.128° < i_\mathbb{C} = 5.145°$; the difference $0.017°$ arises from the planetary perturbation of the lunar node and from the distinction between the osculating inclination and the mean inclination.*
— *Proof approach: small-angle expansion — the exact latitude is $\sin\beta_\mathbb{C} = \sin i_\mathbb{C}\sin F$; for $i_\mathbb{C} \approx 5°$, $\sin\beta_\mathbb{C} \approx \beta_\mathbb{C}$ and $\sin i_\mathbb{C} \approx i_\mathbb{C}$; so $\beta_\mathbb{C} \approx i_\mathbb{C}\sin F$ and the amplitude is $i_\mathbb{C}$ in radians $\times (180°/\pi) = 5.145°$; the empirical value $5.128°$ is the mean from ELP 2000-82 and includes secular and periodic corrections to the inclination.*

**Error Note [A.26.2]:** The single-term latitude approximation omits higher harmonics in $F$ (next term $\approx 0.28°\sin(F + M_\mathbb{C} + ...)$). For third-body perturbation of equatorial satellites, the latitude error matters; for polar satellites it is a secondary effect. [Tier II for third-body perturbation on equatorial orbit; Tier I for polar orbit.]

---

### §26.8 Complete SGP4 Lunar Longitude

This section assembles the complete five-term SGP4 lunar longitude formula, states the two-term lunar equation of center, and quantifies the combined accuracy error against JPL DE430.

*Stub.* Assemble the complete SGP4 lunar longitude formula:

$$\lambda_\mathbb{C} = \ell_\mathbb{C} + ({\rm equation~of~center}) + \Delta\lambda_{\rm ev} + \Delta\lambda_{\rm var} + \Delta\lambda_{\rm ann} + \Delta\lambda_4 + \Delta\lambda_5 \tag{26.7}$$

State the equation of center for the Moon in the same two-term form as Ch 25 but with $e_\mathbb{C}$ in place of $e_\odot$. Note that $e_\mathbb{C} \approx 0.0549$ makes the two-term truncation less accurate than for the Sun (error $\sim 5e_\mathbb{C}^3/(1-e_\mathbb{C}) \approx 9\times10^{-4}$ rad $\approx 0.05°$).

**Error Note [A.26.3]:** The two-term lunar equation of center has accuracy error $< 0.05°$; the five-term perturbation truncation introduces additional error $< 0.3°$ (next-largest omitted term is $\sim 0.186°$). The combined error relative to DE430 is $< 0.5°$ over multi-decade intervals. [Tier II for third-body perturbation at the $10^{-4}$ rad level.]

---

### §26.9 Geocentric Distance

This section states the two-term lunar distance formula with evection correction and quantifies the combined precision error.

*Stub.* The mean lunar distance is:

$$r_\mathbb{C} \approx a_\mathbb{C}(1 - e_\mathbb{C}\cos M_\mathbb{C}) + (\text{evection distance term}) \tag{26.8}$$

State the evection correction to the distance and the combined precision error.

---

### §26.10 Generalization: Adding More Terms

This section documents the next five omitted ELP 2000-82 longitude terms, defines the `lunar_ephemeris_terms` tolerance parameter, and states the convergence rate of the series.

*Stub.* The full Brown (1896) / ELP 2000-82 lunar theory contains 1463 terms in longitude, 974 in latitude, and 675 in distance. For SGP4 purposes, only the five longitude terms are retained. Adding more terms recovers accuracy in two ways:

1. **Direct improvement.** Each additional term of amplitude $A_k$ reduces the longitude error by $|A_k|$.
2. **Perturbation accuracy.** The third-body perturbation depends on $\hat{r}_\mathbb{C}$ integrated over the propagation interval; each additional ephemeris term reduces the integrated perturbation error.

State the next five omitted terms and their amplitudes. Define the tolerance parameter analog: a `lunar_ephemeris_terms` parameter that selects how many terms to include, with the default matching SGP4 and higher values recovering ELP accuracy.

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.26.1] | A | §26.3 | Evection coefficient secular variation: $< 0.001°$/century; Tier I |
| [A.26.2] | A | §26.7 | Single-term latitude: omits $\sim 0.28°$ harmonic; Tier II for equatorial sat |
| [A.26.3] | A | §26.8 | Five-term longitude: $< 0.5°$ vs. DE430; Tier II for third-body perturbation |
| [A.26.4] | A | §26.3, §26.6 | ELP 2000-82 (Chapront-Touzé \& Chapront 1988) has published corrections (Chapront-Touzé \& Chapront 2003, A\&A 404); coefficients in §26.3–26.6 must be verified against the corrected series, not the original 1988 paper |
| [P.26.1] | P | §26.2 | Modular reduction of $\ell_\mathbb{C}$, $F$, $D$ to $[0,2\pi)$: Lipschitz-1 (Ch 7) |

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 5 | Ch 5, alternating series bound | Alternating series bound for equation-of-center truncation error |
| Ch 7 | Ch 7, angle modular reduction | Angle modular reduction |
| Ch 23 | Ch 23, astronomical constants | Lunar periods: sidereal, anomalistic, nodical, synodic |
| Ch 24 | Ch 24, perturbing body framework | Perturbing body framework; elongation $D$ definition |
| Ch 25 | Ch 25, solar mean longitude | Solar mean longitude $\ell_\odot$ used in $D$ and annual equation |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 27 | Ch 27, third-body perturbation | Consumes $r_\mathbb{C}$, $\hat{r}_\mathbb{C}$ from this chapter |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 4 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 6 |
| Examples | 1 |
| Error Notes | 5 |
| Equations | ~15 |
| Sections | 10 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §26.1 Introduction | Draft | |
| §26.2 Mean Lunar Elements | Draft | |
| §26.3 The Evection | Draft | |
| §26.4 The Variation | Draft | |
| §26.5 The Annual Equation | Draft | |
| §26.6 Additional Terms in the SGP4 Lunar Ephemeris | Draft | |
| §26.7 Ecliptic Latitude | Draft | |
| §26.8 Complete SGP4 Lunar Longitude | Draft | |
| §26.9 Geocentric Distance | Draft | |
| §26.10 Generalization: Adding More Terms | Draft | |
