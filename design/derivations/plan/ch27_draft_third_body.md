# Draft Plan: Chapter 27 — Third-Body Perturbation Theory

**Part VI: Third-Body Perturbations**

Target header: `third_body.h`

---

## Objectives

1. Derive the third-body disturbing potential from the two-body gravitational force of an external mass, via the Legendre polynomial expansion.
2. Show that the $n = 0$ term is a constant (irrelevant to orbital motion) and the $n = 1$ term cancels identically (indirect acceleration, Jacobi's theorem).
3. Derive the quadrupole ($n = 2$) disturbing function and show it dominates for solar/lunar perturbations.
4. Express $\cos S$ — the cosine of the angle between the satellite position and the perturbing body direction — in terms of orbital elements via spherical trigonometry.
5. Derive the orbit-averaged $\langle r^2 P_2(\cos S)\rangle$ and decompose into secular, long-period, and short-period parts.
6. Derive the five secular rate equations $(\dot\Omega, \dot\omega, \dot{M}, \dot{e}, \dot{i})$ from the von Zeipel averaging of the quadrupole disturbing function.
7. State the relative perturbation strength $\mu_\mathbb{C}/r_\mathbb{C}^3 \approx 2.2\,\mu_\odot/r_\odot^3$.
8. State both the dual quaternion and $7\times7$ matrix form of the third-body perturbation.

---

## Section Structure

### §27.1 Introduction

This section previews the structure of the third-body perturbation derivation, places it in the context of the SGP4 propagator, and states the derivation-first methodology applied throughout.

*Stub.* This chapter is the theoretical culmination of Part VI. The solar and lunar position models of Ch 25–26 are inputs; the outputs are the secular orbital element rates driven by third-body gravity. These rates are added to the Brouwer secular rates of Ch 16–17 to produce the complete secular model of the SGP4 propagator.

The derivation is entirely ab initio: every coefficient is derived from the Legendre expansion, not accepted from a reference implementation.

---

### §27.2 The Third-Body Disturbing Potential

This section derives the disturbing potential of an external mass from Newton's law, expands it as a Legendre series, and establishes the convergence domain $r < r_*$.

*Stub.* For a perturbing body of mass $m_*$ at position $\vec{r}_*$ (geocentric), the potential at the satellite position $\vec{r}$ is:

$$U_* = -\frac{\mu_*}{|\vec{r}_* - \vec{r}|} + \frac{\mu_*}{r_*^2}\vec{r}\cdot\hat{r}_* \tag{27.1}$$

The first term is the direct attraction of the perturber; the second term is the indirect (tidal) term that accounts for the acceleration of the Earth by the perturber. Together they form the disturbing function $R_* = U_*$.

Expand $1/|\vec{r}_* - \vec{r}|$ in Legendre polynomials using the generating function:

$$\frac{1}{|\vec{r}_* - \vec{r}|} = \frac{1}{r_*}\sum_{n=0}^\infty\left(\frac{r}{r_*}\right)^n P_n(\cos S) \tag{27.2}$$

where $S$ is the angle between $\vec{r}$ and $\vec{r}_*$ (valid for $r < r_*$).

**Notation table:**

| Symbol | Meaning |
|--------|---------|
| $\mu_*$ | Gravitational parameter of the perturbing body |
| $\vec{r}_*$ | Geocentric position of perturbing body |
| $r_*$ | $|\vec{r}_*|$ |
| $\hat{r}_*$ | Unit vector from Earth to perturbing body |
| $S$ | Angle between satellite position and perturbing body direction |
| $\cos S$ | $\hat{r}\cdot\hat{r}_*$ |
| $P_n$ | Legendre polynomial of degree $n$ |
| $R_*$ | Third-body disturbing function |
| $n_*$ | Mean motion of perturbing body |
| $\lambda_*$ | Ecliptic longitude of perturbing body |
| $\delta_*$ | Declination of perturbing body |
| $\alpha_*$ | Right ascension of perturbing body |

---

### §27.3 Cancellation of $n=0$ and $n=1$ Terms

This section proves that the $n=0$ Legendre term produces no orbital force and that the $n=1$ direct term exactly cancels the indirect (tidal) term, so the disturbing function begins at $n=2$.

*Stub.*

**$n = 0$ term.** $P_0(\cos S) = 1$, so the $n=0$ contribution to $R_*$ is $-\mu_*/r_*$ — a constant independent of the satellite position. Constants in the disturbing function produce no force.

**$n = 1$ term.** $P_1(\cos S) = \cos S = \hat{r}\cdot\hat{r}_*$. The $n=1$ contribution to $R_*$ is $-(\mu_*/r_*^2)\vec{r}\cdot\hat{r}_*$, which exactly cancels the indirect term $+(\mu_*/r_*^2)\vec{r}\cdot\hat{r}_*$ in $U_*$. This cancellation is the statement that the Earth's frame is non-inertial (it accelerates toward the perturber) and the indirect term corrects for this.

**Theorem 27.3.1** (Indirect acceleration cancellation). *In the geocentric frame, the $n=1$ term of the direct expansion $-(\mu_*/r_*^2)(r/r_*)\cos S$ and the indirect term $+(\mu_*/r_*^2)\vec{r}\cdot\hat{r}_*$ cancel exactly, leaving $R_* = \sum_{n\geq 2}(-\mu_*/r_*)(r/r_*)^n P_n(\cos S)$.*
— *Proof approach: direct substitution — write $P_1(\cos S) = \cos S = \hat{r}\cdot\hat{r}_*$; the $n=1$ direct term is $-(\mu_*/r_*^2)\vec{r}\cdot\hat{r}_*$ (with sign from the disturbing function definition); the indirect term is $+(\mu_*/r_*^2)\vec{r}\cdot\hat{r}_*$ (the tidal correction for the Earth's acceleration toward the perturber); sum is exactly zero. The cancellation is an identity requiring no approximation.*

---

### §27.4 The Quadrupole Approximation

This section retains only the $n=2$ term of the Legendre series, quantifies the octupole correction factor $r/r_*$, and establishes the error tier for the quadrupole approximation at LEO and GEO altitudes.

*Stub.* The leading non-zero term in the disturbing function is the $n=2$ quadrupole:

$$R_*^{(2)} = -\frac{\mu_*}{r_*^3}r^2 P_2(\cos S) = -\frac{\mu_*}{r_*^3}\frac{r^2}{2}(3\cos^2 S - 1) \tag{27.3}$$

This is the tidal potential — it depends on the satellite's position $\vec{r}$ through both $r^2$ and $\cos S$.

**Proposition 27.4.1** (Validity of quadrupole approximation). The neglected $n=3$ term is smaller than the $n=2$ term by a factor of $r/r_* \approx 7000/385000 \approx 0.018$ for LEO vs. Moon, and $7000/(1.5\times10^8) \approx 4.7\times10^{-5}$ for LEO vs. Sun. [Tier I for Sun; Tier I for LEO Moon; Tier II for high-altitude Moon perturbation.]
— *Proof approach: Legendre ratio bound — the ratio of consecutive terms in the Legendre expansion is $|P_{n+1}/P_n| \leq r/r_*$ for all argument values $|\cos S| \leq 1$; the formal bound follows directly from the generating function.*

**Error Note [A.27.1]:** The quadrupole truncation omits the octupole ($n=3$) and higher terms. For the Moon, the octupole introduces an additional $\sim 1.8\%$ correction to the secular perturbation rate at LEO altitudes. [Tier I for LEO; Tier II for altitudes above 20000 km.]

---

### §27.5 Decomposition of $\cos S$ in Orbital Elements

This section expands $\cos S = \hat{r}\cdot\hat{r}_*$ as an explicit function of the satellite orbital angles $(u, \Omega, i)$ and the perturbing body coordinates $(\alpha_*, \delta_*)$ via spherical trigonometry.

*Stub.* The satellite's geocentric unit vector $\hat{r}$ and the perturbing body's unit vector $\hat{r}_*$ can each be expressed in equatorial coordinates. The inner product $\cos S = \hat{r}\cdot\hat{r}_*$ then expands in terms of the satellite's orbital angles $(u, \Omega, i)$ and the perturbing body's equatorial coordinates $(\alpha_*, \delta_*)$.

Using spherical trigonometry (spherical cosine formula):

$$\cos S = \sin\delta_*\sin\phi + \cos\delta_*\cos\phi\cos(\lambda - \alpha_*) \tag{27.4}$$

where $\phi$ is the satellite latitude and $\lambda$ is the satellite longitude. Express $\phi$ and $\lambda$ in orbital elements $(u, \Omega, i)$:

$$\sin\phi = \sin i\sin u, \qquad \tan(\lambda - \Omega) = \cos i\tan u \tag{27.5}$$

**Theorem 27.5.1** ($\cos S$ in orbital elements). *The angle $\cos S$ expands as $\cos S = A_1\cos u + A_2\sin u + A_3$ where $A_1 = \cos\delta_*(\cos(\Omega - \alpha_*) - \sin i\sin(\Omega - \alpha_*)\cot i_*)$, with analogous expressions for $A_2$ and $A_3$ involving $\sin\delta_*$, $\cos i$, $\sin i$.*
— *Proof approach: spherical cosine formula — apply the spherical triangle formula $\cos S = \sin\delta_*\sin\phi + \cos\delta_*\cos\phi\cos(\lambda - \alpha_*)$; substitute the satellite latitude and longitude in orbital elements from equations (27.5); expand using sum-to-product identities and collect terms by their dependence on $u$. The result is a linear combination of $\cos u$, $\sin u$, and constant terms — verify by direct substitution of the equatorial limit $i=0$.*

---

### §27.6 Orbit Average of $\langle r^2 P_2(\cos S)\rangle$

This section performs the orbit average of the quadrupole disturbing function over the satellite mean anomaly, identifies which terms survive as secular, long-period, and short-period, and states the key Keplerian integrals.

*Stub.* The orbit average is taken over the mean anomaly $M$ (equivalently the eccentric anomaly $E$ or true anomaly $\nu$) at fixed $(\Omega, \omega, i)$. The key Keplerian averages are:

$$\langle r^2 \rangle = a^2\left(1 + \tfrac{3}{2}e^2\right), \qquad \langle r^2\cos 2u \rangle = -\tfrac{5}{2}a^2 e^2\cos 2\omega \tag{27.6}$$

Substitute the decomposition from §27.5. The $\cos 2u$ and $\sin 2u$ terms survive the orbit average only in combination with $\omega$ (they produce the long-period terms). The purely $u$-dependent terms average to zero (short-period). The $u$-independent terms survive directly (secular).

**Theorem 27.6.1** (Orbit-averaged quadrupole). *The secular part is $\langle r^2 P_2(\cos S)\rangle_{\rm sec} = \frac{1}{2}a^2(1 + \frac{3}{2}e^2)(A_3^2 - \frac{1}{3}) + \text{terms involving } A_1, A_2, \omega$; the long-period part contains $\sin 2\omega$ and $\cos 2\omega$ multiplied by $e^2$; the short-period part vanishes upon averaging.*
— *Proof approach: Keplerian average integrals — substitute $r^2 = a^2(1 - e\cos E)^2$ and change variable to eccentric anomaly $E$; use the standard Keplerian orbit averages $\langle r^2\rangle = a^2(1 + \frac{3}{2}e^2)$, $\langle r^2\cos 2u\rangle = -\frac{5}{2}a^2 e^2\cos 2\omega$, $\langle r^2\sin 2u\rangle = -\frac{5}{2}a^2 e^2\sin 2\omega$ (Ch 12, orbit-averaging lemmas); insert the expansion from Theorem 27.5.1. Cross-reference Ch 12 for the derivation of these Keplerian integrals.*

---

### §27.7 The Five Secular Rate Corrections

This section applies the Lagrange planetary equations to the orbit-averaged quadrupole to derive explicit secular rate formulas for all five orbital elements affected by the third-body perturbation.

*Stub.* From the orbit-averaged quadrupole disturbing function $\langle R_*^{(2)}\rangle$, apply the Lagrange planetary equations (Ch 12) to obtain secular rates:

$$\dot\Omega_* = \frac{\partial\langle R_*^{(2)}\rangle}{\partial H} \cdot \frac{1}{Gn_0 a^2 \sqrt{1-e^2}\sin i} \tag{27.7}$$

and similarly for $\dot\omega_*$, $\dot{M}_*$, $\dot{e}_*$, $\dot{i}_*$.

**Theorem 27.7.1** (Third-body secular node rate). *$\dot\Omega_* = -\frac{3n_*^4}{2n}\frac{\mu_*}{\mu}\frac{a^3}{(1-e^2)^2}\cos i_* f_1(i, i_*, \Omega - \Omega_*)$ where $f_1$ is an explicit function of the satellite and perturber inclinations and node difference.*
— *Proof approach: Lagrange planetary equations — differentiate the secular part of $\langle R_*^{(2)}\rangle$ with respect to the Delaunay momentum $H = \mu a(1-e^2)\cos i$ using the Delaunay form of the Lagrange equations (Ch 11); the chain rule introduces the $1/\sin i$ factor at the denominator, but the $\partial\langle R\rangle/\partial H$ numerator has a compensating $\sin i$ from the $A_3$ term so the singularity cancels for $i = 0$ (verify explicitly).*

**Theorem 27.7.2** (Third-body secular perigee rate). *$\dot\omega_* = \frac{3n_*^4}{2n}\frac{\mu_*}{\mu}\frac{a^3}{(1-e^2)^2} f_2(i, i_*, \Omega - \Omega_*)$ where $f_2$ is an explicit trigonometric function.*
— *Proof approach: Lagrange planetary equations — differentiate $\langle R_*^{(2)}\rangle_{\rm sec}$ with respect to the Delaunay momentum $G = \mu a(1-e^2)$; identify the $\cos i$ terms that give the node contribution and the $e^2$ terms that give the eccentricity contribution.*

**Theorem 27.7.3** (Third-body secular mean motion rate). *$\dot{M}_{*,\rm sec} = -\frac{3n_*^4}{n}\frac{\mu_*}{\mu}\frac{a^3}{(1-e^2)^{3/2}} f_3(i, i_*, \Omega - \Omega_*)$.*
— *Proof approach: Lagrange planetary equations — differentiate $\langle R_*^{(2)}\rangle_{\rm sec}$ with respect to the Delaunay action $L = \mu a$; the result is an overall shift in the mean motion that appears as a secular drift in $M$.*

**Theorem 27.7.4** (Third-body secular eccentricity rate). *The secular eccentricity rate $\dot{e}_{*,\rm sec} = 0$ for the pure quadrupole when the perturber orbit is coplanar with the satellite orbit ($\Omega_* = \Omega$, $i_* = i$); when the orbits are inclined to each other, $\dot{e}_{*,\rm sec}$ is generally non-zero and driven by the $\sin 2\omega$ long-period term.*
— *Proof approach: Lagrange planetary equations — differentiate $\langle R_*^{(2)}\rangle_{\rm sec}$ with respect to $\omega$ (the argument of perigee); the secular part is $\omega$-independent (by construction of the secular average), so $\dot{e} = -\partial\langle R\rangle/\partial M \cdot (1/na^2 e)\sqrt{1-e^2} = 0$ for the secular average. The non-zero $\dot{e}$ arises only from the long-period terms.*

**Theorem 27.7.5** (Third-body secular inclination rate). *$\dot{i}_{*,\rm sec} = \frac{3n_*^4}{2n}\frac{\mu_*}{\mu}\frac{a^3}{(1-e^2)^2\sin i} f_5(i, i_*, \Omega - \Omega_*)$.*
— *Proof approach: Lagrange planetary equations — differentiate $\langle R_*^{(2)}\rangle_{\rm sec}$ with respect to $\Omega$ (the longitude of the ascending node); the chain rule through the $\cos i$ dependence introduces the $1/\sin i$ factor; verify the result is consistent with the secular Hamiltonian structure (the secular inclination rate must equal $-\partial\langle R\rangle/\partial H \cdot \cot i / (n a^2\sqrt{1-e^2})$ by the Delaunay canonical form).*

All five theorems require complete algebraic derivation in the Develop phase; here only the structure is stated.

---

### §27.8 Relative Perturbation Strength

This section computes the numerical ratio of lunar to solar tidal perturbation amplitudes and states the matched-pair values of $\mu_\mathbb{C}$, $r_\mathbb{C}$, $\mu_\odot$, $r_\odot$ used by SGP4.

*Stub.* The ratio of the lunar to solar tidal perturbation strength is:

$$\frac{\mu_\mathbb{C}/r_\mathbb{C}^3}{\mu_\odot/r_\odot^3} = \frac{4.903\times10^{12}\ {\rm m^3/s^2}}{(3.844\times10^8\ {\rm m})^3} \cdot \frac{(1.496\times10^{11}\ {\rm m})^3}{1.327\times10^{20}\ {\rm m^3/s^2}} \approx 2.19 \tag{27.8}$$

**Proposition 27.8.1** (Moon-to-Sun perturbation ratio). *The Moon's secular perturbation rate is $\approx 2.19\times$ the Sun's for the same satellite orbit.*
— *Proof approach: direct ratio computation — insert $\mu_\mathbb{C} = 4.903\times10^{12}$ m³/s², $r_\mathbb{C} = 3.844\times10^8$ m, $\mu_\odot = 1.327\times10^{20}$ m³/s², $r_\odot = 1.496\times10^{11}$ m into equation (27.8); verify each constant against the SGP4 matched-pair sources from App A.*

State the numerical values of $\mu_\mathbb{C}$, $r_\mathbb{C}$, $\mu_\odot$, $r_\odot$ used by SGP4 and identify the matched-pair sources.

**Error Note [A.27.2]:** The ratio 2.19 is computed at mean distances; the instantaneous ratio varies with the lunar orbit eccentricity and phase, ranging from approximately 1.9 to 2.5. The secular perturbation uses the mean value; the variation averages out over the lunar anomalistic period. [Tier I for secular rates; Tier II for short-period rates.]

---

### §27.9 Generalization: Octupole ($n = 3$)

This section states the structure of the octupole disturbing function, identifies which octupole terms are secular vs. long-period, and states when the octupole correction must be included.

*Stub.* The octupole term is:

$$R_*^{(3)} = -\frac{\mu_*}{r_*^4}r^3 P_3(\cos S) \tag{27.9}$$

For the Moon at LEO, the amplitude ratio $R_*^{(3)}/R_*^{(2)} \sim r_{\rm sat}/r_\mathbb{C} \approx 0.018$, so the octupole contributes at the $\sim 1.8\%$ level.

The octupole introduces new angular dependencies in $\langle r^3 P_3(\cos S)\rangle$ that do not simplify as cleanly as the quadrupole. State the structure of the octupole orbit average and identify which terms are secular vs. long-period. Note that the octupole secular terms involve $\langle r^3\rangle$ (an odd moment that vanishes when averaged over the full Kepler orbit in the Keplerian mean motion), so the octupole contributes only to long-period and short-period effects for circular orbits.

---

### §27.10 State Framework

This section expresses the third-body tidal acceleration in both the dual-quaternion and $7\times7$ matrix representations of Ch 2 and identifies the tidal tensor block.

*Stub.*

**Dual quaternion form.** The third-body tidal acceleration $\vec{a}_* = \nabla R_*^{(2)}$ is a force in the inertial frame. It enters the velocity dual quaternion equation of motion as a perturbation to the dual part of $\dot{\hat{\Omega}}_b$:

$$\dot{\hat{\Omega}}_b = [\hat\Omega_b, \hat\Omega_b]/2 + \hat{F}_*/2 \tag{27.10}$$

(symbolic form — the precise commutator structure follows from Ch 2, §2.4). The tidal force $\vec{a}_*$ is encoded as a traceless Hermitian perturbation in the dual part of $\hat{F}_*$ (Principle 2 from Ch 2).

**7×7 matrix form.** The tidal acceleration adds a $3\times 3$ symmetric matrix (the tidal tensor $\partial^2 U_*/\partial r_i\partial r_j$) to the force block of the $7\times7$ state equation. The secular perturbation is represented as a near-identity $7\times7$ matrix applied at each time step.

**Transport coupling.** The velocity perturbation from the tidal force produces a position change through the transport coupling term in the $7\times7$ block structure (Ch 2, §2.8). Over one orbit period, this coupling integrates to the secular rate expressions of §27.7.

**Proposition 27.10.1** (Tidal tensor encoding). *The tidal tensor $T_{ij} = \partial^2 R_*^{(2)}/\partial r_i\partial r_j = (\mu_*/r_*^3)(3\hat{r}_{*i}\hat{r}_{*j} - \delta_{ij})$ is traceless and symmetric; it embeds in the velocity dual quaternion as the imaginary dual part of the perturbation $\hat{F}_*/2$.*
— *Proof approach: direct differentiation — compute $\partial^2/\partial r_i\partial r_j[-(\mu_*/r_*^3)(3\cos^2 S - 1)r^2/2]$ at constant $r_*$, $\hat{r}_*$; collect into matrix form; verify tracelessness: $\mathrm{tr}(T) = (\mu_*/r_*^3)(3 - 3) = 0$. The embedding in the dual quaternion follows from Ch 2, Principle 2 (force encoding).*

**Proposition 27.10.2** (7×7 tidal block). *In the $7\times7$ matrix representation (Ch 2, §2.8), the tidal perturbation adds the $3\times3$ tidal tensor $T$ to the $(4{:}6, 1{:}3)$ force block; the $(1{:}3, 1{:}3)$ transport block is unchanged.*
— *Proof approach: block matrix expansion — write the $7\times7$ equation $\dot{\mathbf{x}} = A\mathbf{x}$ where $\mathbf{x} = (r, v, 1)^T$; the tidal force $\vec{a}_* = T\vec{r}$ contributes $T$ to the $(4{:}6, 1{:}3)$ block; verify that the remaining blocks are not modified by this force.*

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.27.1] | A | §27.4 | Quadrupole truncation: octupole $\sim 1.8\%$ for LEO Moon; Tier I LEO; Tier II high altitude |
| [A.27.2] | A | §27.8 | Mean vs. instantaneous Moon-to-Sun ratio: factor 1.9–2.5 range; Tier I for secular |
| [A.27.3] | A | §27.7 | Solar/lunar ephemeris errors (Ch 25–26) propagate into secular rates; bounded by [A.25.1], [A.26.3] |
| [A.27.4] | A | §27.7 | Third-body secular rate formulas in Hoots \& Roehrich (1980) SR3 have known transcription errors corrected by Vallado et al. (2006); all five secular rate theorems must be derived from first principles and cross-checked against the corrected Vallado source, not the original SR3 |
| [P.27.1] | P | §27.6 | Orbit average commutes with Legendre expansion only for $r < r_*$; always satisfied for LEO vs. Moon |

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 2 | Thm 2.4.x | Velocity dual quaternion equation of motion |
| Ch 2 | Thm 2.7.x–2.8.x | Adjoint bridge; $7\times7$ matrix; tidal tensor block |
| Ch 11 | Ch 11, Delaunay variables | Delaunay variables; Lagrange planetary equations |
| Ch 12 | Ch 12, orbit averaging | Orbit averaging; disturbing function; Gauss/Lagrange variational equations |
| Ch 13 | Ch 13, Legendre polynomials | Legendre polynomials and spherical harmonic expansion |
| Ch 16–17 | Ch 16–17, Brouwer secular rates | Secular rate structure; third-body rates add to Brouwer rates |
| Ch 24 | Ch 24, perturbing body framework | Perturbing body framework; $r_*$, $\hat{r}_*$ |
| Ch 25 | Ch 25, solar ephemeris | Solar position and distance |
| Ch 26 | Ch 26, lunar ephemeris | Lunar position and distance |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 33 | Ch 33, secular update | Third-body secular rates $(\dot\Omega_*, \dot\omega_*, \dot{M}_*)$ added to Brouwer secular rates |
| Ch 35 | Ch 35, DPINIT/DPSEC/DPPER | Deep-space initialization and secular accumulation use third-body rates from this chapter |
| Ch 38 | Ch 38, output error budget | Third-body perturbation errors contribute to the overall error budget |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 1 |
| Theorems | 7 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 4 |
| Examples | 1 |
| Error Notes | 5 |
| Equations | ~20 |
| Sections | 10 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §27.1 Introduction | Draft | |
| §27.2 The Third-Body Disturbing Potential | Draft | |
| §27.3 Cancellation of $n=0$ and $n=1$ Terms | Draft | |
| §27.4 The Quadrupole Approximation | Draft | |
| §27.5 Decomposition of $\cos S$ in Orbital Elements | Draft | |
| §27.6 Orbit Average of $\langle r^2 P_2(\cos S)\rangle$ | Draft | |
| §27.7 The Five Secular Rate Corrections | Draft | |
| §27.8 Relative Perturbation Strength | Draft | |
| §27.9 Generalization: Octupole ($n=3$) | Draft | |
| §27.10 State Framework | Draft | |
