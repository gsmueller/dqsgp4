# Draft Plan: Chapter 24 — The Celestial Body Framework

**Part VI: Third-Body Perturbations**

Target header: `celestial_body.h`

---

## Objectives

1. Define the data structure that represents a perturbing celestial body (Sun or Moon) for use in the third-body perturbation theory of Ch 27.
2. Establish the geocentric orbital element representation for each body: mean longitude, mean anomaly, node, inclination, and their secular rates.
3. Derive the geocentric rates (mean anomaly, node, perigee) from the periods of Ch 23.
4. Explain why the Moon requires significantly more perturbation terms than the Sun, and provide a precise quantitative justification.
5. Classify all quantities by tier.

---

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $a_*$ | Semi-major axis of perturbing body orbit | §24.2 |
| $e_*$ | Eccentricity of perturbing body orbit | §24.2 |
| $i_*$ | Inclination of perturbing body orbit to ecliptic | §24.2 |
| $\Omega_*$ | Longitude of ascending node of perturbing body | §24.2 |
| $\omega_*$ | Argument of perigee of perturbing body | §24.2 |
| $M_*$ | Mean anomaly of perturbing body | §24.2 |
| $n_*$ | Mean motion of perturbing body: $n_* = 2\pi/T_*$ | §24.3 |
| $r_*$ | Geocentric distance to perturbing body | §24.5 |
| $\mu_*$ | Gravitational parameter of perturbing body | §24.4 |
| $\hat{r}_*$ | Unit vector from Earth to perturbing body | §24.5 |
| $T_{\rm trop}$ | Tropical year (from Ch 23) | Ch 23 |

---

## Section Structure

### §24.1 Introduction

This section introduces the abstract celestial body framework shared by both the solar and lunar ephemerides and previews the common data structure consumed by the third-body perturbation theory of Ch 27.

*Stub.* Chapters 25 and 26 derive specific ephemerides for the Sun and Moon. Chapter 27 uses those positions in the third-body perturbation calculation. This chapter provides the abstract framework shared by both perturbing bodies: a common data structure and a common set of geometric quantities. The framework is general enough to admit additional perturbing bodies (Jupiter, Venus) with only new constant tables.

---

### §24.2 The Perturbing Body Data Structure

This section defines the orbital element data structure for a perturbing body and specifies the linear secular extrapolation used to compute its state at an arbitrary epoch.

*Stub.* Define the quantities that characterize a perturbing body for use in the perturbation calculation:

**Definition 24.2.1** (Geocentric orbital elements for a perturbing body). The set $(a_*, e_*, i_*, \Omega_*, \omega_*, M_*)$ at a reference epoch, together with secular rates $(\dot{a}_*, \dot{e}_*, \dot{i}_*, \dot{\Omega}_*, \dot{\omega}_*, \dot{M}_*)$. For the Sun and Moon these rates are tabulated from the ephemeris models of Ch 25–26.

**Definition 24.2.2** (Perturbing body state at epoch $t$). *Stub — linear extrapolation from reference epoch using secular rates.*

Note: the perturbing body orbital elements are heliocentric for the Sun (degenerate case: the Sun is at the origin) and geocentric for the Moon. For the Sun, the geocentric position is the negative of the heliocentric Earth position.

**Notation table:**

| Symbol | Meaning |
|--------|---------|
| $a_*$ | Semi-major axis of perturbing body orbit |
| $e_*$ | Eccentricity of perturbing body orbit |
| $i_*$ | Inclination of perturbing body orbit to ecliptic (or equator) |
| $\Omega_*$ | Longitude of ascending node of perturbing body |
| $\omega_*$ | Argument of perigee (perihelion) of perturbing body |
| $M_*$ | Mean anomaly of perturbing body |
| $n_*$ | Mean motion of perturbing body: $n_* = 2\pi / T_*$ |
| $r_*$ | Geocentric distance to perturbing body |
| $\mu_*$ | Gravitational parameter of perturbing body |
| $\hat{r}_*$ | Unit vector from Earth to perturbing body |

---

### §24.3 Geocentric Rates

This section derives the geocentric mean motion and precession rates for the Sun and Moon directly from the period definitions of Ch 23.

*Stub.* Derive the secular rates from the periods of Ch 23.

**Proposition 24.3.1** (Solar geocentric mean motion). $n_\odot = 2\pi / T_{\rm trop}$ rad/day (using the tropical year since the reference is the vernal equinox). State the numerical value.
— *Proof approach: direct substitution — insert $T_{\rm trop} = 365.2422$ days; verify agreement with the SGP4 matched-pair value from Ch 23, §23.6.*

**Proposition 24.3.2** (Lunar geocentric mean motion). $n_\mathbb{C} = 2\pi / T_{\rm sid}$ rad/day (using the sidereal month). State the numerical value.
— *Proof approach: direct substitution — insert $T_{\rm sid} = 27.3217$ days; verify against the SGP4 constant from Ch 23, §23.7.*

**Proposition 24.3.3** (Perigee and node rates). Derive $\dot\omega_\mathbb{C}$ and $\dot\Omega_\mathbb{C}$ from the anomalistic and nodical months (Ch 23): $\dot\omega_\mathbb{C} = n_\mathbb{C}^{\rm anom} - n_\mathbb{C}^{\rm sid}$ and $\dot\Omega_\mathbb{C} = n_\mathbb{C}^{\rm sid} - n_\mathbb{C}^{\rm node}$.
— *Proof approach: phase counting argument — the argument is the same as Ch 23, Proposition 23.3.1; rearrange the period relations to express rates as differences of mean motions.*

---

### §24.4 Why the Moon Requires More Terms than the Sun

This section derives the quantitative criterion for ephemeris truncation order and proves that equivalent perturbation accuracy requires more lunar terms than solar terms.

*Stub.* The truncation order of the ephemeris expansion is set by the required accuracy of the third-body perturbation calculation. The perturbation amplitude scales as $(\mu_*/r_*^3) \cdot a_{\rm sat}^2$ (quadrupole term from Ch 27). The ratio of Moon-to-Sun perturbation amplitudes is:

$$\frac{\mu_\mathbb{C}/r_\mathbb{C}^3}{\mu_\odot/r_\odot^3} \approx \frac{4.90 \times 10^{12} / (3.84 \times 10^5)^3}{1.33 \times 10^{20} / (1.50 \times 10^8)^3} \approx 2.2 \tag{24.1}$$

So the Moon produces $\approx 2.2 \times$ larger secular perturbations than the Sun on a satellite orbit.

The additional factor driving more lunar terms is the large fractional variation in the Moon's distance ($r_\mathbb{C}$ varies by $\pm 5.5\%$ due to eccentricity) and the large amplitude of the evection ($\pm 1.27°$ in longitude) and variation ($\pm 0.66°$) perturbations. For the Sun, the orbit is nearly circular ($e_\odot \approx 0.0167$) and the equation of center is at most $\pm 1.9°$. The Sun position is adequately described by a two-term approximation (Ch 25); the Moon requires at least five terms for equivalent fractional accuracy (Ch 26).

**Theorem 24.4.1** (Required ephemeris terms for equivalent perturbation accuracy). *The number of ephemeris terms required for a fractional third-body perturbation rate error $\varepsilon$ scales as $k \geq \lceil \log(C\varepsilon) / \log(e_*) \rceil$ for a body with orbital eccentricity $e_*$, where $C$ is a constant depending on the longitude amplitude of the first omitted term.*
— *Proof approach: series truncation bound — use the Leibniz alternating series bound on the equation-of-center expansion: the $k$-th omitted term has amplitude $O(e_*^k)$; the resulting longitude error $\delta\lambda_* \sim A_k$ induces a perturbation rate error $\delta\dot\Omega \sim (\mu_*/r_*^3)\,a^2\,\delta\hat{r}_*$ where $\delta\hat{r}_* \sim \delta\lambda_*$; set this less than the target tolerance and solve for $k$.*

**Error Note [A.24.1]:** Using a one-term solar ephemeris (mean longitude only) introduces position error of order $2e_\odot \approx 0.033$ rad, which at $r_\odot = 1.5\times10^8$ km gives $\sim 5000$ km absolute error. The resulting error in $\hat{r}_\odot$ is $\sim 3\times10^{-5}$ rad, small for perturbation purposes. [Tier I for solar ephemeris in third-body computation.]

---

### §24.5 The Common Perturbation Geometry

This section defines the geometric interface between the ephemeris models and the third-body perturbation theory, specifying exactly what Ch 27 requires from each perturbing body.

*Stub.* For both Sun and Moon, the perturbation theory of Ch 27 requires the geocentric unit vector $\hat{r}_*$ and the distance $r_*$ at each epoch. Define the common geometric interface:

**Definition 24.5.1** (Perturbing body geometry). *Given a perturbing body state, the perturbation geometry consists of: (a) the geocentric distance $r_*$; (b) the ecliptic longitude $\lambda_*$ and latitude $\beta_*$ (for the Moon: nonzero); (c) the equatorial unit vector $\hat{r}_{*,{\rm eq}}$ after obliquity rotation.*

State the obliquity rotation (ecliptic to equatorial) used by SGP4. Note that Ch 29 provides the full frame transformation; this section uses a simplified two-rotation model.

---

### §24.6 Tier Classification

This section classifies all perturbing body constants and their derived rates by accuracy tier, following the Ch 1 framework.

*Stub.* Classify all perturbing body constants and rates by tier.

| Quantity | Body | Tier | Error Type |
|---------|------|------|-----------|
| $\mu_\odot$ | Sun | I | A |
| $\mu_\mathbb{C}$ | Moon | I | A |
| $r_\odot$ (mean) | Sun | I | A |
| $r_\mathbb{C}$ (mean) | Moon | II | A — varies $\pm 5.5\%$ |
| $n_\odot$ | Sun | I | A |
| $n_\mathbb{C}$ | Moon | I | A |
| $\dot\Omega_\mathbb{C}$ | Moon | I | A |
| $\dot\omega_\mathbb{C}$ | Moon | I | A |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.24.1] | A | §24.4 | One-term solar ephemeris: 5000 km position error; $3\times10^{-5}$ rad unit vector error; Tier I for perturbation |
| [A.24.2] | A | §24.3 | Linear secular rate extrapolation: valid to $\sim 0.01°$/century; Tier I for $|t - t_0| < 100$ years |

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 23 | Ch 23, astronomical constants | Provides all fundamental periods and epochs |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 25 | Ch 25, solar ephemeris | Implements the solar ephemeris using this framework |
| Ch 26 | Ch 26, lunar ephemeris | Implements the lunar ephemeris using this framework |
| Ch 27 | Ch 27, third-body perturbation | Consumes $r_*$, $\hat{r}_*$ from this framework |
| Ch 29 | Ch 29, frame definitions | Full ecliptic-to-equatorial frame transformation |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 5 |
| Theorems | 1 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 3 |
| Examples | 1 |
| Error Notes | 2 |
| Equations | ~6 |
| Sections | 6 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §24.1 Introduction | Draft | |
| §24.2 The Perturbing Body Data Structure | Draft | |
| §24.3 Geocentric Rates | Draft | |
| §24.4 Why the Moon Requires More Terms than the Sun | Draft | |
| §24.5 The Common Perturbation Geometry | Draft | |
| §24.6 Tier Classification | Draft | |
