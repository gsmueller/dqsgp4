# Draft Plan: Chapter 28 — Tesseral Resonance

**Part VII: Orbital Resonance** | Implementation file: `resonance.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\psi_{lmpq}$ | Angular argument of a Kaula term | Ch 15 |
| $\dot{\psi}_{lmpq}$ | Time derivative of angular argument | Ch 15 |
| $\omega_E$ | Earth rotation rate (rad/s) | Ch 29, Def. 29.5.1 |
| $n$ | Satellite mean motion (rad/s) | Ch 8 |
| $\alpha, \beta$ | Commensurability integers, $\gcd(\alpha,\beta)=1$ | §28.2 |
| $\lambda$ | Resonant phase angle (libration coordinate) | §28.4 |
| $K$ | Pendulum equation driving coefficient | §28.4 |
| $\lambda_{nm}$ | Phase angle from $C_{nm}, S_{nm}$: $\frac{1}{m}\arctan(S_{nm}/C_{nm})$ | §28.5 |
| $F_{lmp}(I)$ | Inclination function (Ch 15, Def. 15.3.1) | Ch 15 |
| $G_{lpq}(e)$ | Eccentricity function / Hansen coefficient (Ch 15, Def. 15.4.1) | Ch 15 |
| $d_{lmpq}$ | Resonant forcing coefficient (specific Kaula-mode amplitude) | §28.6 |
| $h$ | Integration step size for secular resonance integrator | §28.7 |
| $t_{\mathrm{epoch}}$ | Epoch restart time for integrator | §28.7 |

---

## Objectives

1. Derive the commensurability condition for tesseral resonance from the Kaula frequency spectrum.
2. Show how repeat-ground-track orbits satisfy this condition for integer $\alpha:\beta$ ratios.
3. Derive the pendulum equation for the resonant phase angle from the disturbing function, identifying equilibrium points and their stability.
4. Compute the driving coefficient $K$ as a sum of $F_{lmp} \times G_{lpq}$ products over the resonant harmonic set.
5. Relate phase angles $\lambda_{nm}$ to the spherical harmonic coefficients $C_{nm}, S_{nm}$.
6. Develop the secular integrator for near-resonant and resonant orbits: step size selection, epoch restart, energy conservation check.
7. Derive the specific coefficients d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 from first principles (no matching to SR3 code values).
8. State Framework: resonant perturbations in both dual-quaternion and 7×7 form.
9. Generalization: higher-order tesseral harmonics.

## Section Structure

### §28.1 Introduction

This section motivates tesseral resonance physically and places it within the perturbation hierarchy, providing a road map of the chapter and forward references.

Stub: Context — why tesseral harmonics are negligible for most orbits but dominant near resonance. Physical picture: the satellite repeatedly passes over the same longitude where the gravity anomaly is greatest, so small perturbations accumulate coherently. Forward-reference table:

| Section | Feeds | Role |
|---------|-------|------|
| §28.2 | Ch 15, Ch 29 | Kaula frequency spectrum and Earth rotation rate |
| §28.4–§28.6 | Ch 35 | Pendulum equation and d-coefficients used by deep-space integrator |
| §28.8 | Ch 2, Ch 27 | State Framework pattern for perturbation forces |

### §28.2 The Perturbation Frequency Spectrum

This section derives the time derivative of each Kaula angular argument and identifies the commensurability condition that makes a term resonant.

Stub: From the Kaula expansion (Ch 15, §15.6), each term has angular argument $\psi_{lmpq} = (l-2p)\omega + (l-2p+q)M + m(\Omega - \theta_{\mathrm{GMST}})$. Compute $\dot{\psi}_{lmpq}$ using secular rates from Brouwer (Ch 16–17). State the commensurability condition: $\dot{\psi}_{lmpq} \approx 0$ when $n/\omega_E = \beta/\alpha$ (near-integer ratio).

1. **Definition 28.2.1** (Kaula angular argument): $\psi_{lmpq} = (l-2p)\omega + (l-2p+q)M + m(\Omega - \theta_{\mathrm{GMST}})$.
2. **Definition 28.2.2** (resonant term): a Kaula term $(l,m,p,q)$ is resonant when $|\dot{\psi}_{lmpq}| < \epsilon$ for a specified threshold $\epsilon$.
3. **Lemma 28.2.1** ($\dot{\psi}$ zero-crossing): For fixed $(l,m,p,q)$, $\dot{\psi}_{lmpq}(a) = 0$ defines a resonant semi-major axis $a_{\mathrm{res}}$, and the zero-crossing set forms a continuous shell in $(a, e, I)$ space. — *Proof approach: express $\dot{\psi}$ as a smooth function of $a$ via Kepler's third law, show $\partial\dot{\psi}/\partial a \neq 0$ at the zero, and apply the implicit function theorem.*

### §28.3 Repeat Ground Tracks

This section connects the commensurability condition to the physical geometry of repeat ground tracks and computes the resonant semi-major axis for standard satellite classes.

Stub:

1. **Definition 28.3.1** (repeat ground track orbit): an orbit for which $n \cdot T_{\mathrm{repeat}} = 2\pi\beta$ and $\omega_E \cdot T_{\mathrm{repeat}} = 2\pi\alpha$, with $\gcd(\alpha,\beta)=1$.
2. **Proposition 28.3.1** (resonant semi-major axis): $a_{\mathrm{res}} = (\mu/\omega_E^2)^{1/3}(\alpha/\beta)^{2/3}$. — *Proof approach: substitute $n = \omega_E\beta/\alpha$ into Kepler's third law $n^2 a^3 = \mu$ and solve for $a$.*

Table of $\alpha:\beta$ ratios: GPS ($\alpha=2, \beta=1$, $a_{\mathrm{res}} \approx 26\,560$ km), GEO ($\alpha=1,\beta=1$, $a_{\mathrm{res}} \approx 42\,164$ km), LEO ($\alpha \approx 16, \beta=1$, $a_{\mathrm{res}} \approx 6\,570$ km). [A.28.1] The Brouwer secular node regression shifts the effective $n$ slightly; first-order correction to $a_{\mathrm{res}}$.

**Example 28.3.1** (GPS resonance): For GPS with $\alpha=2$, $\beta=1$: $n = \omega_E/2 = 3.646\times10^{-5}$ rad/s. From Kepler's third law with $\mu = 398600.8$ km$^3$/s$^2$: $a_{\mathrm{res}} = (\mu/n^2)^{1/3} = 26\,559.7$ km. Source: direct computation from WGS72 constants (Appendix A).

### §28.4 The Pendulum Equation

This section derives the pendulum equation governing resonant libration from the truncated disturbing function and classifies the stability of the equilibrium points.

Stub: Near $\dot{\psi}=0$, expand $\ddot{\psi}$ to linear order in $\Delta a = a - a_{\mathrm{res}}$ and to zeroth order in $\lambda = \psi \mod 2\pi$. Derive pendulum equation $\ddot{\lambda} + K\sin\lambda = 0$ where $K > 0$ for the dominant harmonic.

1. **Definition 28.4.1** (libration): bounded oscillation of $\lambda$ about a stable equilibrium; circulation: unbounded monotone evolution of $\lambda$.
2. **Theorem 28.4.1** (stability of equilibria): For $K > 0$, the equilibrium $\lambda^* = 0$ is stable (center); $\lambda^* = \pi$ is unstable (saddle). The separatrix has energy $E_{\mathrm{sep}} = K$, and the libration half-width in $\Delta a$ is $\Delta a_{\mathrm{max}} = 2\sqrt{K/|(\partial\dot{\psi}/\partial a)|}$. — *Proof approach: linearize the pendulum equation about each equilibrium and compute eigenvalues of the linearized system ($\pm i\sqrt{K}$ for center, $\pm\sqrt{K}$ for saddle); separatrix energy from the conserved Hamiltonian $H = \tfrac{1}{2}\dot{\lambda}^2 - K\cos\lambda$ evaluated at the saddle.*

[A.28.2] Truncation to single harmonic; coupling between multiple resonant terms gives chaos near the separatrix.

**Example 28.4.1** (GEO libration period): For a geostationary orbit, $K \approx 1.5\times10^{-15}$ rad$^2$/s$^2$ from the dominant $J_{22}$ term. Libration period $T_{\mathrm{lib}} = 2\pi/\sqrt{K} \approx 840$ days. Source: derived from $F_{220}(0^\circ) \cdot G_{200}(0) \cdot R_{22}$ product with WGS72 $C_{22}, S_{22}$.

### §28.5 Phase Angles from Spherical Harmonic Coefficients

This section expresses the Kaula phase angles in terms of the measured spherical harmonic coefficients $C_{nm}, S_{nm}$ and classifies the inherited measurement error.

Stub:

1. **Definition 28.5.1** (amplitude and phase): $R_{nm} = \sqrt{C_{nm}^2 + S_{nm}^2}$, $\lambda_{nm} = \frac{1}{m}\operatorname{atan2}(S_{nm}, C_{nm})$.
2. **Lemma 28.5.1** (harmonic decomposition): $C_{nm}\cos(m\lambda) + S_{nm}\sin(m\lambda) = R_{nm}\cos(m\lambda - m\lambda_{nm})$. — *Proof approach: apply the trigonometric addition formula to expand $R_{nm}\cos(m\lambda - m\lambda_{nm})$, then match coefficients of $\cos(m\lambda)$ and $\sin(m\lambda)$ and verify consistency with Definition 28.5.1.*

Note atan2 branch selection for correct quadrant. [M.28.1] $C_{nm}, S_{nm}$ are measured quantities (Tier II); $\lambda_{nm}$ inherits $\sigma_m$.

**Example 28.5.1** (WGS72 $J_{22}$ phase angle): $C_{22} = 1.57447\times10^{-6}$, $S_{22} = -9.03812\times10^{-7}$. Then $R_{22} = 1.815\times10^{-6}$, $\lambda_{22} = \frac{1}{2}\operatorname{atan2}(-9.038\times10^{-7}, 1.574\times10^{-6}) \approx -15.0^\circ$. Source: WGS72 gravity model coefficients (Appendix A).

### §28.6 The Driving Coefficient K

This section computes the driving coefficient $K$ as an explicit sum over the resonant harmonic set and derives each d-coefficient from first principles using the inclination and eccentricity function recurrences.

Stub: $K$ is the sum over all resonant $(l,m,p,q)$ of $\frac{\mu}{a^{l+1}} \left(\frac{a_E}{a}\right)^l F_{lmp}(I)\, G_{lpq}(e)\, R_{lm}$ evaluated at $a_{\mathrm{res}}$. Derive by substituting the resonant term of the Kaula expansion into the pendulum equation.

1. **Definition 28.6.1** (driving coefficient): $K = \sum_{\mathrm{res}} d_{lmpq}$ where $d_{lmpq} = \frac{\mu}{a_{\mathrm{res}}^{l+1}} (a_E/a_{\mathrm{res}})^l F_{lmp}(I) G_{lpq}(e) R_{lm}$.
2. **Theorem 28.6.1** (d-coefficient derivation): Each $d_{lmpq}$ equals the product $F_{lmp}(I) \cdot G_{lpq}(e) \cdot R_{lm} \cdot \mu a_E^l / a_{\mathrm{res}}^{2l+1}$, where $F_{lmp}$ and $G_{lpq}$ are computed from the recurrence relations of Ch 15 (inclination function definition, eccentricity function definition). — *Proof approach: substitute the single resonant Kaula term into $\ddot{\lambda} = -\partial \mathcal{R}/\partial\lambda$, identify the coefficient of $\sin\lambda$ after differentiation, and collect all constant prefactors into the d-coefficient definition.*

Table of dominant contributors for synchronous and half-day resonance. Derivation of each specific coefficient: **d2201** $(l=2,m=2,p=0,q=1)$, **d2211** $(l=2,m=2,p=1,q=1)$, **d3210** $(l=3,m=2,p=1,q=0)$, **d3222** $(l=3,m=2,p=2,q=2)$, **d4410** $(l=4,m=4,p=1,q=0)$, **d4422** $(l=4,m=4,p=2,q=2)$, **d5220** $(l=5,m=2,p=2,q=0)$, **d5232** $(l=5,m=2,p=3,q=2)$, **d5421** $(l=5,m=4,p=2,q=1)$, **d5433** $(l=5,m=4,p=3,q=3)$. Each derived from the $F_{lmp} \times G_{lpq}$ product at nominal inclination and eccentricity. [A.28.3] Coefficients depend on inclination and eccentricity; the SR3 implementation uses fixed values valid near reference orbit.

**Example 28.6.1** (d2201 for half-day resonance): For $I = 63.4^\circ$, $e = 0.01$: $F_{220}(63.4^\circ) = 1.781$, $G_{201}(0.01) \approx -0.010$ (Hansen coefficient), $R_{22} = 1.815\times10^{-6}$, $a_{\mathrm{res}} = 26\,560$ km. Product yields $d_{2201} \approx -2.4\times10^{-16}$ rad$^2$/s$^2$. Source: $F_{lmp}$, $G_{lpq}$ from Ch 15 recurrences evaluated at stated $(I,e)$; $R_{22}$ from WGS72.

### §28.7 The Secular Resonance Integrator

This section specifies the fixed-step numerical integrator for resonant orbits, derives the step size selection criterion from the libration period, and defines the epoch restart protocol.

Stub: Near resonance the long-period terms cannot be averaged; a numerical integrator is required.

1. **Definition 28.7.1** (secular integrator): advances $\lambda, \dot{\lambda}$ using only the slowly varying (resonant) terms of the disturbing function; step size $h = 720$ min (standard).
2. **Definition 28.7.2** (epoch restart): re-initialization of the integrator state from the current osculating elements when accumulated error exceeds a threshold or when the propagation direction reverses.
3. **Proposition 28.7.1** (step size requirement): $h \ll 2\pi/K^{1/2}$ (libration period). — *Proof approach: the Nyquist condition requires sampling at least twice per libration period; apply the local truncation error bound $O(h^3)$ for the fixed-step scheme, then solve $h < T_{\mathrm{lib}}/10$ as the practical stability criterion.*

[P.28.1] Fixed-step integrator; no adaptive control; accumulates precision error $O(h^2)$ per step. Energy conservation check as a diagnostic.

**Example 28.7.1** (GEO step size): For GEO, $T_{\mathrm{lib}} \approx 840$ days $= 1.21\times10^6$ min. Step size $h = 720$ min gives $h/T_{\mathrm{lib}} \approx 6\times10^{-4}$, well within the stability bound. After 30 days (60 steps), accumulated precision error $\sim 60 \cdot O(h^2) \approx 3\times10^7$ min$^2$ in phase. Source: direct computation from $K$ in Example 28.4.1.

### §28.8 State Framework: Resonant Corrections in Both Forms

This section expresses the resonant perturbation force in both the dual-quaternion and 7×7 matrix representations and places it within the parallel-row error budget.

Stub: The resonant disturbing function adds a force perturbation to the satellite dynamics.

1. **Definition 28.8.1** (resonant force perturbation): the gradient $-\nabla_{\mathbf{r}} \mathcal{R}_{\mathrm{res}}$ of the resonant part of the Kaula disturbing function.
2. **Proposition 28.8.1** (dual-quaternion resonant correction): The perturbation $\delta\hat{\Omega}$ added to $\hat{\Omega}_b$ follows the same structure as the third-body force (Ch 27, third-body state framework section), with the Kaula gradient replacing the lunisolar gradient. — *Proof approach: identify the resonant force as a conservative gradient field $-\nabla_{\mathbf{r}}\mathcal{R}_{\mathrm{res}}$, then apply the Ch 2 Principle 2 (force/perturbation in dual-quaternion form) by direct substitution of the gradient into the velocity dual-quaternion update equation.*

In 7x7 form: force vector in the translation subblock, assembled per Ch 2, §2.8. Both forms carry the parallel-row error budget (Ch 2, §2.10). Cross-reference Ch 27, §27.8 for the pattern.

### §28.9 Generalization: Higher-Order Tesseral Harmonics

This section extends the resonance framework to higher-degree tesseral harmonics and states the Chirikov overlap criterion for the onset of chaotic motion.

Stub: The same framework applies to $(l,m)=(6,\ldots)$ and higher.

1. **Definition 28.9.1** (resonance overlap): two resonant terms with adjacent commensurabilities whose separatrices intersect in phase space.
2. **Proposition 28.9.1** (Chirikov overlap criterion): chaotic motion onset when the sum of libration half-widths exceeds the frequency separation between adjacent resonances. — *Proof approach: compute the separatrix half-widths in semi-major axis from Theorem 28.4.1 applied to each resonance, then compare to the frequency gap $|\dot{\psi}_1 - \dot{\psi}_2|$ evaluated at the midpoint semi-major axis; overlap ($\Delta a_1 + \Delta a_2 > |a_{\mathrm{res},1} - a_{\mathrm{res},2}|$) signals chaos.*

List the next terms in order of magnitude for GEO, GPS, and LEO resonances. Note the condition for multiple simultaneous resonances.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 15, Kaula expansion | §28.2, §28.6 | $F_{lmp}$, $G_{lpq}$, angular arguments $\psi_{lmpq}$ |
| Ch 16–17, Brouwer secular rates | §28.2 | Secular rates entering $\dot{\psi}_{lmpq}$ |
| Ch 27, third-body force architecture | §28.8 | State Framework pattern for perturbation forces |
| Ch 29, sidereal time | §28.2, §28.5 | Earth rotation rate $\omega_E$ in commensurability condition |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 35, deep-space pipeline | §28.7 | Resonance integrator with 720-min step |
| Appendix A | §28.5 | WGS72/WGS84 $C_{nm}, S_{nm}$ values |
| Appendix C | §28.6 | d-coefficient code locations |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.28.1] | M | §28.5 | Phase angles $\lambda_{nm}$ inherit Tier II measurement error from $C_{nm}, S_{nm}$ |
| [P.28.1] | P | §28.7 | Fixed-step secular integrator accumulates $O(h^2)$ precision error per step; no adaptive control |
| [A.28.1] | A | §28.3 | Brouwer secular node regression shifts effective $n$; first-order correction to $a_{\mathrm{res}}$ needed |
| [A.28.2] | A | §28.4 | Single-harmonic pendulum approximation; multi-harmonic coupling produces chaos near separatrix |
| [A.28.3] | A | §28.6 | Specific d-coefficients derived at nominal $(I, e)$; errors grow away from reference orbit |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 9 |
| Theorems | 2 |
| Lemmas | 2 |
| Corollaries | 0 |
| Propositions | 4 |
| Examples | 5 |
| Error Notes | 5 |
| Equations | ~20 |
| Sections | 9 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §28.1 | Draft | Introduction and context |
| §28.2 | Draft | Perturbation frequency spectrum |
| §28.3 | Draft | Repeat ground tracks |
| §28.4 | Draft | Pendulum equation |
| §28.5 | Draft | Phase angles from spherical harmonic coefficients |
| §28.6 | Draft | Driving coefficient K |
| §28.7 | Draft | Secular resonance integrator |
| §28.8 | Draft | State Framework: resonant corrections |
| §28.9 | Draft | Generalization: higher-order tesseral harmonics |
