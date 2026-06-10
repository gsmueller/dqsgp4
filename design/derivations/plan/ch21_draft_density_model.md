# Draft Plan: Chapter 21 — The Power-Law Density Model

**Part V: Atmospheric Drag**

Target header: `density_model.h`

---

## Objectives

1. Derive the spherically symmetric power-law density model $\rho(r) = \rho_0\bigl((q_0 - s)/(r - s)\bigr)^\tau$ from first principles (dimensional analysis and the motivation for a geometric-tail profile).
2. Compare this model against NRLMSISE-00 and Jacchia reference models to bound the accuracy error it introduces.
3. Characterize the fitting parameters $s$ (reference altitude offset), $q_0$ (reference geocentric distance), and the perigee-dependent adjustments.
4. Quantify the non-rotating atmosphere assumption and its ~5% error for prograde LEO.
5. Develop generalizations: rotating atmosphere, variable $\tau$ across altitude bands, solar-flux scaling.

---

## Section Structure

### §21.1 Introduction

This section motivates the power-law density model as the analytically tractable surrogate required by orbit-averaged drag integrals.

*Stub.* Part V concerns the secular decay of satellite orbits due to aerodynamic drag. This chapter establishes the atmospheric density model that is the input to the drag-averaged integrals of Ch 22. The model must be analytically integrable over an elliptic orbit — this requirement, not physical fidelity, drives the choice of a power-law profile.

Contrast the philosophy with Ch 14 (equipotential ellipsoid): there, physical fidelity (equipotential surface theory) was paramount. Here, analytic tractability is the binding constraint, with accuracy error [A.21.x] as the price paid.

---

### §21.2 Physical Atmosphere Background

This section surveys thermospheric density structure and identifies why empirical models are not analytically integrable over an elliptic orbit.

*Stub.* Brief survey of the structure of Earth's atmosphere in the altitude range relevant to LEO satellites (~150–2000 km): the thermosphere, the exosphere, and the principal density drivers (solar EUV flux, geomagnetic activity, diurnal bulge). Reference NRLMSISE-00 (Picone et al. 2002) and Jacchia (1970/1977) as the two principal empirical models.

- Definition 21.2.1 (Thermospheric density regime). Formal statement of the altitude band 150--2000 km and the dominant physical processes (EUV heating, molecular diffusion).
- Definition 21.2.2 (Scale height $H$). The local exponential decay length $H = -\rho/(d\rho/dr)$ and its altitude dependence.

State the key fact that drives model choice: neither NRLMSISE-00 nor Jacchia is analytically integrable over an elliptic orbit in closed form. A surrogate model is needed for orbit propagation.

**Error Note [A.21.1]:** Replacing the physical atmosphere with a power-law surrogate introduces accuracy error. Typical density errors are 10–20% at quiet solar conditions, rising to factors of 2–3 during geomagnetic storms. This is the dominant accuracy error in the drag model. [Tier III: model-dependent, environment-variable.]

---

### §21.3 The Power-Law Density Profile

This section derives the power-law density model from the requirement of analytic integrability over a Keplerian orbit.

*Stub.* Motivate the power-law model $\rho(r) = \rho_0\bigl((q_0-s)/(r-s)\bigr)^\tau$ from two directions:

1. **Exponential approximation:** In the thermosphere, $\rho \approx \rho_{\rm ref}\exp(-(r - r_{\rm ref})/H)$ where $H$ is the scale height. A power-law with large $\tau$ approximates an exponential over a limited altitude range.
2. **Analytic tractability:** The orbit-averaged drag integral (Ch 22) becomes a combinatorial sum (Lane 1965) only when the density profile is a power law. For other profiles the integral has no closed form.

**Definition 21.3.1** (Power-law density model). *Formal definition of $\rho(r) = \rho_0[(q_0-s)/(r-s)]^\tau$ with domain $r > s$.*

**Proposition 21.3.1** (Limit to exponential model). *Stub — as $\tau \to \infty$ with $\rho_0(q_0-s)^\tau$ held fixed and the profile normalized at $r_0$, the power-law converges to the exponential $\rho_0\exp(-(r-r_0)/H)$.* — *Proof approach: Taylor expansion of $\ln(1 - x/\tau)$ in the limit $\tau \to \infty$ with $x = \tau(r - r_0)/(q_0 - s)$ fixed.*

**Notation table for this section:**

| Symbol | Meaning |
|--------|---------|
| $\rho(r)$ | Atmospheric mass density at geocentric distance $r$ |
| $\rho_0$ | Reference density (at perigee or reference altitude) |
| $s$ | Altitude offset (reference surface shift, $\approx r_E$) |
| $q_0$ | Reference geocentric distance; $q_0 - s$ is the reference altitude above the shifted surface |
| $\tau$ | Power-law exponent; SGP4 uses $\tau = 4$ |
| $r_p$ | Perigee geocentric distance |
| $H$ | Atmospheric scale height (for exponential comparison) |

---

### §21.4 Parameter Fitting

This section derives the piecewise altitude-band fitting of $(\rho_0, q_0, s)$ to the Jacchia-70 reference model.

*Stub.* The SGP4 fitting procedure: for a given perigee altitude, a reference density $\rho_0$ and reference parameters $(q_0, s)$ are chosen from a lookup table (Hoots and Roehrich 1980, Table 1). Characterize the structure of this table: altitude-band discretization, how the bands were determined from Jacchia model fits.

**Definition 21.4.1** (Perigee density parameters). *Stub — formal definition of the altitude-band piecewise model used by SGP4.*

**Proposition 21.4.1** (Continuity across altitude bands). *Stub — density is continuous at band boundaries by construction; the exponent $\tau$ is constant within a band.* — *Proof approach: piecewise matching conditions at band boundaries with continuity constraint on $\rho(r)$.*

**Error Note [A.21.2]:** The SGP4 density table was fitted to the Jacchia-70 model, not NRLMSISE-00. Post-1980 density models differ from Jacchia-70 by up to 15% in the thermosphere under nominal solar conditions. [Tier III.]

---

### §21.5 Density at Perigee

This section derives the perigee density $\rho_p$ and its sensitivity to perigee altitude errors.

*Stub.* In the orbit-averaged drag integral (Ch 22), the reference density $\rho_0$ is evaluated at perigee: $\rho_p = \rho(r_p)$. Express:

$$\rho_p = \rho_0\left(\frac{q_0 - s}{r_p - s}\right)^\tau \tag{21.1}$$

Derive the perigee sensitivity: $\partial\log\rho_p/\partial r_p = -\tau/(r_p - s)$. For $\tau = 4$ and LEO altitudes (~300 km, $r_p - s \approx 300$ km), this gives approximately $-13.3 \times 10^{-3}$ km$^{-1}$, so a 1 km error in perigee altitude produces roughly 1.3% density error.

- Theorem 21.5.1 (Perigee density sensitivity). $\partial\log\rho_p/\partial r_p = -\tau/(r_p - s)$. — *Proof approach: logarithmic differentiation of the power-law density profile.*

**Example 21.5.1** (ISS perigee density). Compute $\rho_p$ for the ISS orbit: $a = 6780$ km, $e = 0.001$, perigee altitude $\approx 408$ km, with SGP4 parameters $\tau = 4$, $q_0 = 120$ km + $r_E$, $s \approx r_E$. Evaluate the sensitivity $\partial\log\rho_p/\partial r_p \approx -13.3 \times 10^{-3}$ km$^{-1}$ and compare against NRLMSISE-00 tabulated density at 408 km altitude under mean solar conditions ($F_{10.7} = 150$).

**Error Note [P.21.1]:** Precision error in $r_p$ (from the Kepler solution) propagates into $\rho_p$ with amplification factor $\tau/(r_p - s)$.

---

### §21.6 The Non-Rotating Atmosphere Assumption

This section derives the velocity error from neglecting atmospheric co-rotation and bounds the resulting drag force error.

*Stub.* The SGP4 model assumes the atmosphere rotates with the Earth (so the relative velocity of the satellite with respect to the atmosphere is $\vec{v}_{\rm rel} = \vec{v}_{\rm ECI} - \vec{\omega}_E \times \vec{r}$) but the standard derivation of the orbit-averaged drag integral uses the geocentric inertial velocity $v$ rather than the relative velocity $v_{\rm rel}$.

Estimate the error: $|\vec{\omega}_E \times \vec{r}| \approx 0.465$ km/s at the equator; $v_{\rm orbital} \approx 7.7$ km/s for LEO. The fractional velocity error is $\approx 6\%$, and the drag force scales as $v^2$, so the resulting force error is $\approx 12\%$ for equatorial prograde orbits, ~5% when orbit-averaged.

**Definition 21.6.1** (Non-rotating atmosphere approximation). *Stub — formal statement of the approximation.*

**Error Note [A.21.3]:** The non-rotating atmosphere approximation introduces accuracy error of approximately 5% in orbit-averaged drag for prograde LEO. [Tier II for LEO; Tier I for polar orbits.]

---

### §21.7 Comparison with Physical Models

This section quantifies the accuracy error of the power-law surrogate by tabulating density ratios against NRLMSISE-00 at representative altitudes.

*Stub.* Tabulate representative density comparisons between the SGP4 power-law model (fitted to Jacchia-70), the Jacchia-70 model itself, and NRLMSISE-00, at selected altitudes and solar flux levels. Identify the altitude regime where the power-law approximation is worst (near the thermopause, ~500 km).

**Example 21.7.1** (Density comparison at selected altitudes). Compute $\rho_{\rm PL}(r)$ from the SGP4 power-law at altitudes 200, 300, 400, 500, and 800 km. Compare against NRLMSISE-00 densities at $F_{10.7} = 100$ (solar minimum) and $F_{10.7} = 200$ (solar maximum). Tabulate absolute density values (kg/m$^3$) and percentage deviations.

**Error Note [A.21.4]:** The combined accuracy error of the density model relative to NRLMSISE-00 is 10–30% under typical conditions, with no single correction that recovers this without refitting the drag coefficient $C_D(A/m)B^*$ to a new model. [Tier III.]

---

### §21.8 Generalization

This section extends the power-law model to rotating atmospheres, variable exponents, and solar-flux scaling.

*Stub.*

**21.8.1 Rotating atmosphere.** Replace $v^2$ with $v_{\rm rel}^2 = |\vec{v} - \vec{\omega}_E \times \vec{r}|^2$ in the drag force. Expand and orbit-average to first order in $\omega_E/n$. Identify which terms enter the mean-motion equation and which are negligible.

- Proposition 21.8.1 (First-order co-rotation correction). The orbit-averaged drag correction from atmospheric co-rotation is proportional to $\omega_E \cos i / n$. — *Proof approach: binomial expansion of $v_{\rm rel}^2$ and orbit averaging of the cross term $\vec{v}\cdot(\vec{\omega}_E \times \vec{r})$.*

**21.8.2 Variable $\tau$ across altitude bands.** The SGP4 lookup table implicitly uses different effective $\tau$ values in different altitude regimes. A continuous-$\tau(r)$ model can represent this; however, the orbit-averaged integral (Ch 22) then requires numerical evaluation.

- Definition 21.8.1 (Variable-exponent density model). Piecewise power-law model $\rho(r) = \rho_k[(q_k - s_k)/(r - s_k)]^{\tau_k}$ on altitude band $[r_k, r_{k+1})$.

**21.8.3 Bessel function expansion for exponential density (Fitzgibbon 1982).** When the density model is exponential ($\rho \propto \exp(-\alpha r)$) rather than power-law, the orbit average requires expanding $\exp(-\alpha r)$ with $r = a(1-e\cos E)$ over an orbit. Two approaches exist:

- **Taylor series in $\alpha a e$** (BH61): $\exp(-\alpha a e \cos E) = \sum_k (-\alpha a e)^k \cos^k E / k!$. Diverges for $\alpha a e > 1$ (moderate eccentricity at low altitude). Requires BH61's Appendix transformation [7] for convergence.
- **Bessel function expansion** (Fitzgibbon 1982, INPE-2746 §4): $\exp[-\alpha(r-r_p)] = \sum_{n=0}^\infty f_{2n} \cos nE$ where $f_{2n} = \exp[-\alpha(a-r_p)] \cdot F_n(\alpha a e)$ and $F_n$ are modified Bessel functions $I_n(\alpha a e)$. Converges uniformly for all $\alpha a e$.

The Bessel expansion is mathematically superior but not used by SGP4 (which uses the power-law model instead). It is relevant if extending the theory to exponential density models or validating BH61's orbit-averaged results.

**21.8.4 Solar flux scaling.** The reference density $\rho_0$ varies with the 10.7 cm solar flux index $F_{10.7}$. A simple linear correction $\rho_0 \to \rho_0(1 + \alpha(F_{10.7} - 150))$ can reduce the accuracy error from Tier III to Tier II for moderate solar conditions.

- Proposition 21.8.2 (Linear solar flux correction). The first-order correction $\rho_0(F_{10.7}) \approx \rho_0(150)(1 + \alpha(F_{10.7} - 150))$ with $\alpha$ fitted to NRLMSISE-00 reduces density error from Tier III to Tier II for $80 < F_{10.7} < 250$. — *Proof approach: linear regression of $\ln\rho$ against $F_{10.7}$ from NRLMSISE-00 tabulations at fixed altitude.*

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | accuracy error tier classification | Tier classification and accuracy error |
| Ch 3 | matched-pair principle | Matched-pair principle: SGP4 $B^*$ absorbs density model error |
| Ch 5 | geometric-tail series bound | Geometric-tail series evaluation (power-law profile analogy) |
| Ch 9 | Kepler equation solution | Kepler solution for $r_p$ (perigee distance) |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 22 | orbit-averaged drag integral | Power-law $\rho(r)$ is the density input to the Lane expansion |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.21.1] | A | §21.2 | Power-law vs. physical atmosphere: 10–20% typical; Tier III |
| [A.21.2] | A | §21.4 | SGP4 table fitted to Jacchia-70; post-1980 model differences up to 15%; Tier III |
| [A.21.3] | A | §21.6 | Non-rotating atmosphere approximation: ~5% orbit-averaged drag error for prograde LEO; Tier II |
| [A.21.4] | A | §21.7 | Combined density model vs. NRLMSISE-00: 10–30%; Tier III |
| [P.21.1] | P | §21.5 | Perigee altitude error amplified by $\tau/(r_p - s)$ into density error |
| [A.21.5] | A | §21.2 | Jacchia (1970) modeled the semi-annual density variation as driven by temperature change; Jacchia (1971, SAO Spec. Rep. 332) showed this is incorrect — the apparent solar-cycle dependence was an artifact — and decoupled the semi-annual variation from temperature; J71 is a direct correction to J70 and should be used in preference to it |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 7 |
| Theorems | 1 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 4 |
| Examples | 2 |
| Error Notes | 6 |
| Equations | ~6 |
| Sections | 9 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §21.1 | Draft | Stub |
| §21.2 | Draft | Stub |
| §21.3 | Draft | Stub |
| §21.4 | Draft | Stub |
| §21.5 | Draft | Stub |
| §21.6 | Draft | Stub |
| §21.7 | Draft | Stub |
| §21.8 | Draft | Stub |
| §21.9 | Draft | Error summary table |
