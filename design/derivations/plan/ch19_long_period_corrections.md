# Draft Plan: Chapter 19 — Long-Period Corrections

**Part IV: Brouwer's Gravitational Perturbation Theory**
**Phase:** Draft
*(No dedicated header — long-period corrections are applied within `short_period.h` or `brouwer.h`)*

---

## Objectives

1. Derive the J₃ long-period perturbation of eccentricity and its $\sin\omega$ dependence.
2. Prove that J₃ (and all odd zonals) contribute no secular terms, from first principles.
3. Derive long-period corrections to argument of perigee $\omega$ and mean longitude.
4. Identify the $\cos(2\omega)$ terms arising from the second-order generating function $S_2$ (from Ch 17, §17.5.3).
5. Derive the long-period mean longitude correction $IL_L$.
6. Generalize to J₅ long-period contributions and their structure.

---

## Section Structure

### §19.1 Introduction

Narrative: the secular theory of Chs. 16–17 captures constant rates of change. The short-period theory of Ch 18 captures oscillations completing one cycle per orbit. Between these two timescales lies a third class: long-period perturbations, completing one cycle per precession period of the argument of perigee (months to years for LEO). Long-period perturbations are not negligible over propagation arcs of days or weeks.

Two distinct sources drive long-period perturbations:
(a) The odd zonal harmonic J₃, which produces a perturbation with $\sin\omega$ dependence and characteristic magnitude $\sim J_3/(J_2 e)$.
(b) The $\cos(2\omega)$ term in the second-order Hamiltonian $F_2^*$ (Ch 17, Theorem 17.5.3), which produces a $J_2^2$ long-period variation.

The J₃ contribution is of particular physical significance: it breaks the symmetry between the ascending and descending halves of the orbit, and it couples with atmospheric drag (Ch 21) in a way that can produce a secular-looking evolution of eccentricity over many orbits.

Roadmap through the chapter. Forward reference to Ch 20 (all corrections combined for osculating elements), Ch 21 (J₃–drag coupling).

**Error orientation.** Long-period corrections are $O(J_3/J_2)$ relative to the mean elements for the J₃ term, and $O(J_2)$ for the $\cos(2\omega)$ term. The precision error [P.19.1] arises from small denominators proportional to $e\dot{\omega}$ in the J₃ formula.

---

### §19.2 Why Odd Zonals Produce No Secular Terms

This section proves rigorously that the orbit average of odd-degree zonal harmonics vanishes, establishing that J₃ (and J₅) contribute only long-period and short-period perturbations.

**Theorem 19.2.1** (Vanishing orbit average of odd zonals)**.** *Stub: for any odd zonal $J_{2k+1}$ ($k \geq 1$), the orbit average of the disturbing function $\langle R_{2k+1} \rangle = 0$. Proof: $P_{2k+1}(\sin\phi) = P_{2k+1}(\sin i \sin u)$ is an odd function of $\sin u$ for $u \in [0, 2\pi)$; the integral over one complete orbit of an odd function of $u$ vanishes by symmetry of the Keplerian orbit under $u \to \pi - u$.* — *Proof approach: parity argument; decompose $P_{2k+1}(\sin i \sin u)$ into Fourier sine series in $u$ (only odd harmonics appear by Theorem 13.6.3), then orbit-average each $\sin(mu)/r^{2k+2}$ term, showing it vanishes by the symmetry $E \to -E$ of the Keplerian orbit*

**Corollary 19.2.1** (No secular rates from J₃)**.** *Stub: $\langle R_3 \rangle = 0$ implies $\partial \langle R_3 \rangle / \partial L = \partial \langle R_3 \rangle / \partial G = \partial \langle R_3 \rangle / \partial H = 0$; therefore J₃ contributes no secular rates at first order. Its non-zero contribution is entirely in the long-period terms.*

**Remark.** *Stub: this result is exact for a Keplerian reference orbit. At second order, the coupling $J_2 \times J_3$ can produce a small secular-like drift; this is beyond the scope of the first-order theory here.*

---

### §19.3 J₃ Long-Period Perturbation of Eccentricity

This section derives the J₃-induced long-period eccentricity correction with its characteristic $\sin\omega$ dependence and $1/e$ singularity structure.

**Definition 19.3.1** (Long-period part of J₃ disturbing function)**.** *Stub: after removing the short-period terms from $R_3$ via the short-period generating function $S_1$, the residual contains terms that depend on $g = \omega$ (argument of perigee) but not on $l = M$. These are the long-period terms $\hat{R}_3^{\text{lp}}$.*

**Theorem 19.3.1** (Long-period J₃ eccentricity correction)**.** *Stub: the long-period change in eccentricity due to J₃ is*

$$\Delta e_{\text{lp}} = \frac{R_E^3 J_3}{2p^3}\frac{\sin i}{e}\left(1 - \tfrac{5}{4}\sin^2 i\right)\sin\omega$$

*Proof via the Lagrange planetary equation $\dot{e} = \ldots \partial \hat{R}^{\text{lp}}/\partial \omega$, combined with integration over the long period. Explicit derivation at Develop.* — *Proof approach: apply the Lagrange planetary equation for $\dot{e}$ to the J₃ long-period disturbing function, integrate over the long-period cycle $2\pi/\dot{\omega}$, and simplify using $p = a(1-e^2)$*

**Corollary 19.3.1** ($\sin\omega$ dependence)**.** *Stub: the factor $\sin\omega$ means the correction is maximum at $\omega = 90°$ (perigee above north pole) and zero at $\omega = 0°$ or $180°$. Physical interpretation: J₃ represents the pear-shape asymmetry of the Earth; it pushes the orbit away from circular symmetry in the direction of $\omega$.*

**[P.19.1]** The formula contains the factor $1/e$ — it diverges as $e \to 0$. This is a genuine singularity of the Brouwer theory. In the Lyddane variables (Ch 18, Definition 18.4.1), the product $\Delta e \cdot \sin\omega = \Delta\beta$ and $\Delta e \cdot \cos\omega = \Delta\alpha$ are both regular as $e \to 0$; the Lyddane form must be used for near-circular orbits.

**[A.19.1]** For LEO with small $e$ ($e \sim 0.001$), $\Delta e_{\text{lp}} / e$ can be $O(1)$ — i.e., the J₃ long-period correction can be comparable in magnitude to the mean eccentricity itself. This is the case that makes J₃ important despite being a small zonal: $J_3/J_2 \approx 10^{-3}$, but $J_3/(J_2 e) \sim O(1)$ for near-circular orbits.

---

### §19.4 Long-Period Corrections to $\omega$ and Mean Longitude

This section derives the J₃ long-period corrections to argument of perigee and mean longitude, including the Lyddane-regular forms.

**Theorem 19.4.1** (Long-period J₃ correction to $\omega$)**.** *Stub:*

$$\Delta\omega_{\text{lp}} = -\frac{R_E^3 J_3}{2p^3}\frac{\cos i\sin i}{e}\left(1 - \tfrac{5}{4}\sin^2 i\right)\cos\omega + \ldots$$

*Derivation from Lagrange equation for $\dot{\omega}$. Note the $\cos\omega$ dependence -- 90° out of phase with $\Delta e_{\text{lp}}$.* — *Proof approach: apply the Lagrange planetary equation for $\dot{\omega}$ to the J₃ long-period disturbing function and integrate; the $\cos\omega$ phase arises from $\partial/\partial e$ acting on $\sin\omega$ terms*

**Theorem 19.4.2** (Long-period J₃ correction to mean longitude)**.** *Stub: the correction $\Delta M_{\text{lp}}$ or equivalently $\Delta\lambda_{\text{lp}} = \Delta M + \Delta\omega + \Delta\Omega$ -- the mean longitude correction that accounts for the long-period drift in the argument of perigee.* — *Proof approach: combine $\Delta M$, $\Delta\omega$, $\Delta\Omega$ from the respective Lagrange equations, simplifying the sum to a single expression in $(a, e, i, \omega)$*

**Corollary 19.4.1** (Lyddane form)**.** *Stub: write $\Delta\alpha_{\text{lp}} = \Delta(e\cos\omega)$ and $\Delta\beta_{\text{lp}} = \Delta(e\sin\omega)$ as the regular forms.*

---

### §19.5 Long-Period Terms from $F_2^*$: the $\cos(2\omega)$ Part

This section derives the $J_2^2$ long-period eccentricity and inclination variations from the $\cos(2g')$ term of the second-order Hamiltonian.

**Theorem 19.5.1** ($J_2^2$ long-period eccentricity variation)**.** *Stub: the $\cos(2g')$ term of $F_2^*$ (Ch 17, Theorem 17.5.3) generates a long-period variation in eccentricity and inclination with frequency $2\dot{\omega}$. From Hamilton's equations applied to $F_2^{*\text{lp}}\cos(2g')$:*

$$\delta e_{J_2^2} = C_{J_2^2}(\eta, \theta) \sin(2\omega), \quad \delta i_{J_2^2} = D_{J_2^2}(\eta, \theta) \sin(2\omega)$$

*where $C$, $D$ are functions of $\eta$ and $\theta = \cos i$ to be derived at Develop.* — *Proof approach: apply Hamilton's equations $\dot{G} = -\partial F_2^{*\text{lp}}\cos(2g')/\partial g'$ and $\dot{H} = -\partial F_2^{*\text{lp}}\cos(2g')/\partial h'$, then integrate the resulting $\sin(2\omega)$ terms over the long-period timescale $\pi/\dot{\omega}$*

**Corollary 19.5.1** (Magnitude comparison)**.** *Stub: $\delta e_{J_2^2}$ is $O(J_2)$ while $\Delta e_{J_3}$ is $O(J_3/e)$. For LEO with $e \sim 0.001$, $J_3/e \sim J_3 \times 10^3 \sim 10^{-4}$, while $J_2 \sim 10^{-3}$. The J₃ term dominates for small $e$; the $J_2^2$ term dominates for larger $e$.*

---

### §19.6 The Long-Period Mean Longitude Correction $IL_L$

This section defines and derives the SGP4 long-period mean longitude correction $IL_L$, connecting the Brouwer theory to the SR3 implementation.

**Definition 19.6.1** (Long-period mean longitude correction)**.** *Stub: $IL_L$ denotes the accumulated long-period correction to the mean longitude $\lambda = M + \omega + \Omega$. In the SGP4 implementation (SR3, Hoots and Roehrich 1980), $IL_L$ is the notation for the long-period correction applied before the mean-to-osculating step.*

**Theorem 19.6.1** (Formula for $IL_L$)**.** *Stub: $IL_L = \Delta M_{\text{lp}} + \Delta\omega_{\text{lp}} + \Delta\Omega_{\text{lp}}$ from Theorems 19.4.1--19.4.2 and the nodal long-period term. Express in terms of $(a, e, i, \omega)$.* — *Proof approach: algebraic summation of the three long-period angle corrections from Theorems 19.4.1--19.4.2, simplifying common factors*

**Example 19.6.1.** *Stub: compute $IL_L$ for a LEO orbit with $a = 7000$ km, $e = 0.01$, $i = 51.6°$, $\omega = 45°$ using WGS84 ($J_2 = 1.08263 \times 10^{-3}$, $J_3 = -2.53881 \times 10^{-6}$, $a_E = 6378.137$ km). Compute: $p = a(1-e^2) = 6999.3$ km, $\Delta e_{\text{lp}} \approx 2.6 \times 10^{-5}$, $\Delta\omega_{\text{lp}} \approx -1.2 \times 10^{-4}$ rad, $IL_L \approx \Delta M_{\text{lp}} + \Delta\omega_{\text{lp}} + \Delta\Omega_{\text{lp}} \approx -8 \times 10^{-5}$ rad. Verify against SR3 (Hoots and Roehrich 1980) reference values. Source: WGS84.*

**[P.19.2]** The long-period corrections in $IL_L$ add an $O(J_3/(J_2 e))$ term to the mean longitude. For small $e$ this can be $O(\text{arcseconds})$ — larger than the expected position accuracy — confirming that $IL_L$ is not negligible even for the reduced SGP4 model.

---

### §19.7 Generalization: J₅ Long-Period Contributions

This section derives the J₅ long-period correction structure and quantifies the model error from its omission in SGP4.

*Stub: the J₅ zonal harmonic is odd, so $\langle R_5 \rangle = 0$ by Theorem 19.2.1. Its long-period perturbation has the same structural form as J₃: corrections proportional to $\sin\omega$, $\cos\omega$, with factors $(J_5 R_E^5)/(p^5 e)$ and polynomials in $\sin^2 i$. The leading coefficient is $J_5/(J_3) \approx (5.8 \times 10^{-7})/(2.5 \times 10^{-6}) \approx 0.23$, so J₅ contributes roughly 23% of the J₃ long-period effect.*

**Theorem 19.7.1** (J₅ long-period eccentricity correction)**.** *Stub: $\Delta e_{\text{lp}}^{(5)}$ -- analogous to Theorem 19.3.1 with $J_5$, $R_E^5$, and a different inclination polynomial.* — *Proof approach: same Lagrange equation method as Theorem 19.3.1, applied to $P_5(\sin\phi)$ expanded in orbital elements (Ch 13, Theorem 13.6.3 extended to degree 5)*

**Corollary 19.7.1** (Combined J₃ + J₅ correction)**.** *Stub: for high-accuracy propagation, the J₃ and J₅ terms add directly (both $\sin\omega$ dependence). The combined coefficient is $J_3 + (R_E^2/p^2) J_5 P(\sin^2 i)$ where $P$ is a polynomial.*

**[A.19.2]** SGP4 includes J₃ but not J₅ in the long-period corrections. The omitted J₅ term is a model error $\delta_a \sim 23\%$ of $\Delta e_{\text{lp}}^{(3)}$, giving a position error of order $1$–$5$ m for typical LEO.

---

## Cross-References

**Uses (backward)**

| Source | Section | Role |
|---|---|---|
| Ch 12, Lagrange planetary equations | §19.3–19.4 | Mechanism for eccentricity and $\omega$ corrections |
| Ch 13, J₃ disturbing function | §19.3 | Starting point for long-period eccentricity correction |
| Ch 16, Odd zonal orbit average vanishing | §19.2 | Theorem 19.2.1 foundation |
| Ch 16, First-order secular rate $\dot{\omega}$ | §19.3 | Long-period period = $2\pi/\dot{\omega}$ |
| Ch 17, $\cos(2g')$ term of $F_2^*$ | §19.5 | Source of $J_2^2$ long-period terms |
| Ch 18, Lyddane variables | §19.3–19.4 | Regularity at $e = 0$ |

**Feeds (forward)**

| Target | Section | Role |
|---|---|---|
| Ch 20, Osculating elements | §19.3–19.6 | $\Delta e_{\text{lp}}$, $IL_L$, and all long-period corrections applied in final osculating computation |
| Ch 21, Atmospheric drag coupling | §19.3 | J₃ $\sin\omega$ structure interacts with drag to produce slow secular eccentricity evolution |

---

## Error Notes

| Tag | Type | § | Description |
|---|---|---|---|
| [A.19.1] | A | §19.3 | J₃ long-period eccentricity correction $\propto 1/e$; for $e \sim 0.001$ LEO orbits, correction is comparable to mean eccentricity |
| [A.19.2] | A | §19.7 | Omitting J₅ long-period correction (SGP4 standard) introduces $\delta_a \sim 1$–$5$ m position error |
| [P.19.1] | P | §19.3 | Singularity at $e = 0$ in standard form of $\Delta e_{\text{lp}}$; resolved by Lyddane variables |
| [P.19.2] | P | §19.6 | Long-period correction $IL_L$ can reach arcsecond magnitude; not negligible in mean longitude |
| [A.19.3] | A | §19.3 | Brouwer (1959) long-period correction formulas contain a sign error in the $\Delta i_{\text{lp}}$ term; corrected by Lara (2021, J. Astronaut. Sci. 68, 149–180); verify $\Delta e_{\text{lp}}$ and $\Delta i_{\text{lp}}$ signs against corrected source |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 7 |
| Lemmas | 0 |
| Corollaries | 4 |
| Propositions | 0 |
| Examples | 1 |
| Error Notes | 5 |
| Equations | ~12 |
| Sections | 7 |

**Tier classification:** Vanishing of odd-zonal orbit average — Tier I (exact symmetry argument). J₃ long-period eccentricity formula — Tier II (first-order long-period perturbation, derivation from Lagrange equations). $J_2^2$ $\cos(2\omega)$ term — Tier II (second-order, Brouwer Eq. 29 derivation chain). SGP4 J₃ implementation ($IL_L$) — Tier III (matched-pair constant, Chapter 3 applies).

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §19.1 Introduction | Draft | |
| §19.2 Why Odd Zonals Produce No Secular Terms | Draft | |
| §19.3 J₃ Long-Period Perturbation of Eccentricity | Draft | |
| §19.4 Long-Period Corrections to $\omega$ and Mean Longitude | Draft | |
| §19.5 Long-Period Terms from $F_2^*$: the $\cos(2\omega)$ Part | Draft | |
| §19.6 The Long-Period Mean Longitude Correction $IL_L$ | Draft | |
| §19.7 Generalization: J₅ Long-Period Contributions | Draft | |
