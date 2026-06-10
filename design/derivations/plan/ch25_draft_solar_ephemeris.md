# Draft Plan: Chapter 25 — Solar Ephemeris

**Part VI: Third-Body Perturbations**

Target header: `solar_ephemeris.h`

---

## Objectives

1. Derive the mean solar longitude and mean anomaly as linear functions of time from epoch.
2. Derive the equation of center — the difference between true and mean anomaly — from Kepler's equation inversion.
3. Obtain the two-term truncation $\lambda_\odot \approx \ell_\odot + 2e_\odot\sin M_\odot + \tfrac{5}{4}e_\odot^2\sin 2M_\odot$ and quantify its accuracy.
4. Derive the geocentric ecliptic longitude $\lambda_\odot$ and distance $r_\odot$.
5. Generalize to the continued-fraction inversion of Kepler's equation for higher accuracy.

---

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\ell_\odot$ | Mean solar longitude (linear in time) | §25.2 |
| $M_\odot$ | Solar mean anomaly | §25.2 |
| $\nu_\odot$ | Solar true anomaly | §25.3 |
| $\lambda_\odot$ | Geocentric ecliptic longitude of the Sun | §25.4 |
| $r_\odot$ | Geocentric distance to the Sun (AU) | §25.4 |
| $e_\odot$ | Solar orbital eccentricity ($\approx 0.01671$) | §25.2 |
| $n_\odot$ | Solar mean motion (rad/day or rad/min) | §25.2 |
| $T_0$ | Julian centuries from J2000.0 | §25.2 |

---

## Section Structure

### §25.1 Introduction

This section establishes the scope and accuracy requirements of the SGP4 solar ephemeris and motivates the two-term equation-of-center approximation.

*Stub.* The solar ephemeris provides the instantaneous Sun position needed for the solar third-body perturbation (Ch 27) and for the shadow/eclipse model (Ch 31). The SGP4 model uses a low-precision solar position adequate for third-body perturbation, not for precise solar radiation pressure. This chapter derives the SGP4-compatible solar position together with its accuracy bounds.

The solar orbit is treated as a fixed ellipse with constant eccentricity and period (the secular variation of $e_\odot$ is negligible over the SGP4 validity period — see [A.23.1]). The equation of center converts mean anomaly to true anomaly.

---

### §25.2 Mean Solar Motion

This section defines the mean solar longitude and mean anomaly as linear functions of time from epoch, identifies the SGP4 numerical constants, and traces their source to a specific Almanac edition.

*Stub.* The mean solar longitude $\ell_\odot$ increases linearly with time from a reference epoch:

$$\ell_\odot(t) = \ell_{\odot,0} + n_\odot(t - t_0) \tag{25.1}$$

where $n_\odot = 360°/T_{\rm trop}$ is the mean solar motion in degrees/day and $t_0$ is the epoch. The mean solar anomaly is:

$$M_\odot(t) = \ell_\odot(t) - \varpi_\odot \tag{25.2}$$

where $\varpi_\odot = \Omega_\odot + \omega_\odot$ is the solar longitude of perihelion.

**Definition 25.2.1** (Mean solar longitude). *Stub — formal definition with epoch reference.*

**Definition 25.2.2** (Mean solar anomaly $M_\odot$). *Stub — relation to mean longitude and perihelion longitude.*

**Notation table:**

| Symbol | Meaning |
|--------|---------|
| $\ell_\odot$ | Mean solar longitude |
| $M_\odot$ | Solar mean anomaly |
| $\varpi_\odot$ | Solar longitude of perihelion: $\Omega_\odot + \omega_\odot$ |
| $n_\odot$ | Solar mean motion (degrees/day) |
| $e_\odot$ | Solar orbital eccentricity |
| $\nu_\odot$ | Solar true anomaly |
| $\lambda_\odot$ | Solar ecliptic longitude (geocentric) |
| $r_\odot$ | Geocentric solar distance |
| $a_\odot$ | Solar semi-major axis (= 1 AU) |
| $t_0$ | Reference epoch (SGP4: 2000 January 0.5 TT) |
| $\ell_{\odot,0}$ | Mean solar longitude at epoch $t_0$ |

State the SGP4 numerical values of $n_\odot$, $\ell_{\odot,0}$, and $\varpi_\odot$. Identify which Almanac edition they derive from (matched-pair principle, Ch 3).

---

### §25.3 The Equation of Center

This section derives the equation-of-center power series by Lagrange inversion of Kepler's equation and establishes the truncation error bound for the two-term approximation.

*Stub.* The equation of center $q_\odot = \nu_\odot - M_\odot$ is the angular difference between the true anomaly and mean anomaly; it is the inverse of Kepler's equation in the near-circular limit.

Derive the equation-of-center series by inverting $M = E - e\sin E$ and converting to true anomaly via the half-angle formula. The series in powers of $e$:

$$\nu = M + 2e\sin M + \frac{5}{4}e^2\sin 2M + \frac{e^3}{12}(13\sin 3M - 3\sin M) + \cdots \tag{25.3}$$

**Theorem 25.3.1** (Equation of center series). *Stub — derive the series (25.3) from the Lagrange inversion of Kepler's equation, showing the origin of each coefficient.* — *Proof approach: use Lagrange's inversion theorem (Ch 4) on $E = M + e\sin E$ to obtain $E$ as a power series in $e$, then convert $E$ to $\nu$ via the half-angle formulae $\tan(\nu/2) = \sqrt{(1+e)/(1-e)}\tan(E/2)$.*

**Corollary 25.3.1** (Two-term truncation accuracy). *For $e_\odot \approx 0.0167$, the truncation error after two terms is $O(e_\odot^3) \approx 4.7\times10^{-6}$ rad $\approx 0.001°$.* Stub — derive the bound from the Leibniz alternating series argument applied to the higher-order terms.

---

### §25.4 Two-Term Truncation

This section states and bounds the two-term solar longitude approximation used by SGP4 and quantifies the induced error in the third-body perturbation rate.

*Stub.* The SGP4 solar ephemeris uses the two-term equation of center:

$$\lambda_\odot \approx \ell_\odot + 2e_\odot\sin M_\odot + \tfrac{5}{4}e_\odot^2\sin 2M_\odot \tag{25.4}$$

**Proposition 25.4.1** (Two-term truncation). *The two-term approximation of (25.4) satisfies $|\lambda_\odot^{\rm exact} - \lambda_\odot^{\rm 2\text{-}term}| < 5e_\odot^3/(1 - e_\odot) \approx 2.4\times10^{-5}$ rad $\approx 0.0014°$.*
— *Proof approach: series remainder bound — the equation-of-center series is an alternating power series in $e$ (from Theorem 25.3.1); the Leibniz bound (Ch 5, alternating series bound) gives the remainder after $n$ terms as $O(e_\odot^{n+1})$; insert $n=2$, $e_\odot = 0.01671$.*

**Error Note [A.25.1]:** The two-term truncation introduces accuracy error $< 0.0014°$ in solar ecliptic longitude. For the third-body perturbation of a LEO satellite (Ch 27), the solar perturbation amplitude is $\sim 0.003°$/day in the secular node rate; the ephemeris error propagates at the $\sim 0.5\%$ level. [Tier I for third-body perturbation.]

---

### §25.5 Geocentric Distance

This section derives the two-term approximation for the geocentric solar distance and bounds the error from neglecting the $O(e_\odot^2)$ remainder.

*Stub.* The geocentric solar distance is:

$$r_\odot = a_\odot(1 - e_\odot^2)/(1 + e_\odot\cos\nu_\odot) = a_\odot(1 - e_\odot\cos E_\odot) \tag{25.5}$$

For the two-term approximation, express $r_\odot$ directly in terms of $M_\odot$:

$$r_\odot \approx a_\odot(1 - e_\odot\cos M_\odot) + O(e_\odot^2) \tag{25.6}$$

State the precision error from using (25.6) vs. the exact (25.5).

**Proposition 25.5.1** (Solar distance two-term approximation). *The approximation $r_\odot \approx a_\odot(1 - e_\odot\cos M_\odot)$ has error $|r_\odot^{\rm exact} - r_\odot^{\rm approx}| \leq a_\odot e_\odot^2 / (1 - e_\odot) \approx 2.8\times10^{-4}$ AU.*
— *Proof approach: Taylor expansion — expand $r_\odot = a_\odot(1 - e_\odot\cos E)$ in powers of $e_\odot$ using $E = M + e_\odot\sin M + O(e_\odot^2)$; the leading correction is $-a_\odot e_\odot^2(\cos^2 M - 1/2)/2 + O(e_\odot^3)$; bound the maximum of $|{-}a_\odot e_\odot^2\cos^2 M/2|$ over all $M$.*

---

### §25.6 Ecliptic-to-Equatorial Conversion

This section applies the obliquity rotation to convert the solar ecliptic longitude and distance to equatorial Cartesian coordinates, and identifies the SGP4 obliquity value.

*Stub.* The solar ephemeris yields an ecliptic longitude and distance. Convert to equatorial coordinates using the obliquity of the ecliptic $\epsilon$:

$$x_\odot = r_\odot\cos\lambda_\odot, \quad y_\odot = r_\odot\sin\lambda_\odot\cos\epsilon, \quad z_\odot = r_\odot\sin\lambda_\odot\sin\epsilon \tag{25.7}$$

(approximating solar latitude $\beta_\odot = 0$, which is exact in this model). State the obliquity value used by SGP4 and its source (matched-pair principle).

**Error Note [A.25.2]:** The obliquity value introduces accuracy error; see [A.23.4]. The $\beta_\odot = 0$ approximation is exact for the two-body Sun-Earth system; the actual solar latitude is $< 1''$ due to planetary perturbations. [Tier I.]

---

### §25.7 Generalization: Continued Fraction Inversion

This section references the continued-fraction form of Kepler's equation from Ch 9 and quantifies the convergence improvement over the power series for the solar case.

*Stub.* The two-term truncation of §25.4 is adequate for third-body perturbation but not for precise solar radiation pressure or eclipse modeling. The Lagrange inversion series (25.3) converges slowly for larger $e$. The continued fraction representation of Kepler's equation (Ch 4, §4.3; Ch 9, §9.7) provides faster convergence:

$$\nu(M, e) = M + e\sin M / \left(1 - \frac{e\cos M}{1 + \cdots}\right) \tag{25.8}$$

*Stub — reference the general continued fraction from Ch 9 and state the convergence rate improvement over (25.3) for $e_\odot$.*

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.25.1] | A | §25.4 | Two-term equation-of-center truncation: $< 0.0014°$; Tier I for third-body use |
| [A.25.2] | A | §25.6 | Obliquity value and $\beta_\odot = 0$: Tier I |
| [A.25.3] | A | §25.7 | Solar ephemeris vs. DE430: $< 0.01°$ over 200-year window; Tier I for SGP4 purposes |
| [P.25.1] | P | §25.2 | Modular reduction of $\ell_\odot$ to $[0,2\pi)$: no cancellation; Lipschitz-1 (Ch 7) |

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 4 | Ch 4, Lagrange inversion theorem | Lagrange inversion for equation-of-center series; continued fractions |
| Ch 5 | Ch 5, alternating series bound | Alternating series bound for truncation error |
| Ch 7 | Ch 7, angle modular reduction | Angle modular reduction |
| Ch 9 | Ch 9, continued fraction Kepler | Continued fraction form of Kepler's equation |
| Ch 23 | Ch 23, astronomical constants | Solar periods and epoch conventions |
| Ch 24 | Ch 24, perturbing body framework | Perturbing body framework; obliquity rotation |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 27 | Ch 27, third-body perturbation | Consumes $r_\odot$, $\hat{r}_\odot$ from this chapter |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 1 |
| Lemmas | 0 |
| Corollaries | 1 |
| Propositions | 3 |
| Examples | 1 |
| Error Notes | 4 |
| Equations | ~12 |
| Sections | 7 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §25.1 Introduction | Draft | |
| §25.2 Mean Solar Motion | Draft | |
| §25.3 The Equation of Center | Draft | |
| §25.4 Two-Term Truncation | Draft | |
| §25.5 Geocentric Distance | Draft | |
| §25.6 Ecliptic-to-Equatorial Conversion | Draft | |
| §25.7 Generalization: Continued Fraction Inversion | Draft | |
