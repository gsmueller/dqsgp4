# Draft Plan: Chapter 17 — Second-Order (J₂² + J₄) Secular Perturbations

**Part IV: Brouwer's Gravitational Perturbation Theory**
**Target header:** `brouwer.h` (second-order secular portion)
**Phase:** Draft

---

## Objectives

1. Construct the first-order generating function $S_1$ by integrating the short-period part of Brouwer's $F_1$.
2. Write the von Zeipel order-2 equation for the second-order averaged Hamiltonian $F_2^*$.
3. Evaluate all orbit-averaged Poisson-bracket products that comprise $F_2^*$.
4. State the complete $F_2^*$ (Brouwer Eq. 29 / Lara 2021 Eq. 14): secular part and $\cos(2g')$ part.
5. Derive the second-order secular rates in Delaunay variables and convert to Keplerian form.
6. Identify the algebraic origin of the polynomial coefficients $(13, 78, 137)$, $(7, 114, 395)$, $(4, 19)$.
7. Add the $J_4$ contribution to secular rates.
8. Present both the full $\eta$-dependent formulas and the $\eta = 1$ specialization used in SGP4.
9. Provide factored forms avoiding subtractive cancellation.
10. Develop the State Framework. Discuss generalization to third-order $J_2^3$.

---

## Section Structure

### §17.1 Introduction

Narrative: the first-order secular theory of Ch 16 predicts $\dot{\Omega}$, $\dot{\omega}$, $\dot{M}$ to $O(J_2)$. Brouwer's second-order theory extends this to $O(J_2^2)$ by applying the von Zeipel method a second time. The key inputs are the first-order generating function $S_1$ (constructed from the short-period part of $F_1$) and orbit-averaged Poisson brackets of $S_1$ with itself and with the geopotential.

The second-order Hamiltonian $F_2^*$ has two parts: a secular part (function of $a, e, i$ only) and a $\cos(2\omega)$ part (long-period; treated in Ch 19). The secular part yields corrections to $\dot{\Omega}$, $\dot{\omega}$, $\dot{M}$ at $O(J_2^2)$.

The integers $(13, 78, 137)$, $(7, 114, 395)$, $(4, 19)$ appearing in the polynomial coefficients are not arbitrary — they arise from specific combinations of the orbit-averaged Poisson brackets. Each is derived algebraically here.

Roadmap through the chapter. Forward references to Ch 18 (short-period corrections use $S_1$ derived in §17.2), Ch 19 (the $\cos(2g')$ part of $F_2^*$ drives long-period corrections).

**Error orientation.** The second-order theory is exact within its assumptions (Keplerian reference orbit, expansion in $J_2$). Model error [A.17.1] is $O(J_2^3)$, negligible for most applications. Precision errors arise from evaluating the polynomial coefficients and the $\eta$-dependent factors [P.17.1].

---

### §17.2 First-Order Generating Function $S_1$

This section constructs the first-order generating function $S_1$ by integrating the short-period part of Brouwer's $F_1$, establishing the input to the second-order von Zeipel equation.

**Definition 17.2.1** (Short-period part of $F_1$)**.** *Stub: split Brouwer's $F_1$ into secular part $F_1^* = \langle F_1 \rangle$ (Ch 16, Theorem 16.3.1) and short-period part $F_{1p} = F_1 - F_1^*$. Express $F_{1p}$ as a function of $(l, g, h, L, G, H)$.*

**Theorem 17.2.1** (Generating function $S_1$)**.** *Stub: the first-order generating function is defined by the von Zeipel equation $n_0 \partial S_1/\partial l = F_{1p}$, giving*

$$S_1 = \frac{1}{n_0}\int F_{1p}\,dl$$

*where the integral is taken holding $g, h, L, G, H$ fixed and $n_0 = \mu^2/L^3$ is the mean motion. Explicit evaluation via the eccentric anomaly substitution. Cross-reference Brouwer (1959) Eqs. 14--15.* — *Proof approach: indefinite integration of $F_{1p}$ with respect to $l = M$ using the eccentric anomaly substitution $dl = (1 - e\cos E)\,dE$; the integration constant is fixed by the zero-average condition $\langle S_1 \rangle = 0$*

**Corollary 17.2.1** (Structure of $S_1$)**.** *Stub: $S_1$ is a sum of terms involving $\sin(l + 2g)$, $\sin(2l + 2g)$, $\sin(l)$, $\sin(2g)$ — i.e., mixed mean-anomaly and argument-of-perigee harmonics. The coefficient of each term involves $L$, $G$, $H$ through powers of $\eta = G/L$ and $\theta = H/G$.*

**[P.17.1]** $S_1$ involves division by $n_0$ and by small combinations of $\eta$; near-circular orbits require Lyddane substitution (Ch 18) to avoid singularity in $e \to 0$. Near critical inclination, denominator $(5\theta^2 - 1)^{-1}$ diverges; treated in Ch 18, §18.6.

---

### §17.3 The von Zeipel Order-2 Equation

This section states the von Zeipel equation at second order and identifies the Poisson-bracket products that must be orbit-averaged.

**Theorem 17.3.1** (von Zeipel order-2 equation)**.** *Stub: the second-order averaged Hamiltonian $F_2^*$ satisfies*

$$F_2^* = \langle F_2 \rangle + \frac{1}{2}\langle \{S_1, F_1\} \rangle + \frac{1}{6}\langle \{S_1, \{S_1, F_0\}\} \rangle$$

*where $\{\cdot,\cdot\}$ denotes the Poisson bracket in Delaunay variables, $F_0 = -\mu^2/(2L^2)$ is the Keplerian Hamiltonian, and $F_2$ is the second-order perturbation (including $J_4$ -- see §17.6). Derivation via the Deprit triangle (Ch 12, §12.8).* — *Proof approach: apply the Lie-Deprit perturbation method (Ch 12) to order 2, collecting terms $\sim J_2^2$ from the Deprit triangle; the three terms correspond to the direct $F_2$, the Poisson bracket $\{S_1, F_1\}$, and the iterated bracket $\{S_1, \{S_1, F_0\}\}$*

**Remark.** *Stub: the orbit average $\langle \cdot \rangle$ in the von Zeipel equation is taken over the mean anomaly $l$, holding $g, h, L, G, H$ fixed. Because $\{S_1, F_0\} = n_0 \partial S_1/\partial l = F_{1p}$, the last term above reduces to $\frac{1}{2n_0}\langle F_{1p}\, \partial S_1/\partial l \rangle$.*

---

### §17.4 Orbit-Averaged Poisson-Bracket Products

This section evaluates all orbit-averaged Poisson-bracket products that appear in the von Zeipel order-2 equation, reducing $F_2^*$ to explicit polynomials in $\theta = \cos i$ with $\eta$-dependent coefficients.

*Stub: this section evaluates the orbit averages of the products $\langle \sigma_i^2 \rangle$, $\langle \sigma_i\sigma_j \rangle$, $\langle \sigma_i\tau_j \rangle$, $\langle \sigma_i\rho_2 \rangle$ that arise in expanding $\{S_1, F_1\}$, where $\sigma_i$, $\tau_j$, $\rho_2$ denote the trigonometric components of $S_1$ and $F_{1p}$ respectively.*

**Lemma 17.4.1** (Average of $\sin^2$ and $\cos^2$ terms)**.** *Stub: $\langle \sin^2(l + 2g) \rangle = \langle \cos^2(l + 2g) \rangle = 1/2$; $\langle \sin(l+2g)\cos(l+2g) \rangle = 0$.* — *Proof approach: trigonometric orthogonality over the mean anomaly period*

**Lemma 17.4.2** (Coupling terms)**.** *Stub: $\langle \sigma_1 \tau_2 \rangle$ for mixed $(l, 2l+2g)$ products. These vanish for many index combinations by orthogonality.* — *Proof approach: expand the product of trigonometric terms using product-to-sum identities, then apply the orbit average operator; non-zero contributions arise only from matching frequency combinations*

**Lemma 17.4.3** (Secular part of $\langle \{S_1, F_1\} \rangle$)**.** *Stub: after collecting all surviving averages, the secular part of $\langle \{S_1, F_1\} \rangle$ is a polynomial in $\theta^2 = \cos^2 i$ with $\eta$-dependent coefficients.* — *Proof approach: evaluate the Poisson bracket $\{S_1, F_1\} = \sum_k (\partial S_1/\partial q_k \cdot \partial F_1/\partial p_k - \partial S_1/\partial p_k \cdot \partial F_1/\partial q_k)$ in Delaunay variables, then orbit-average using Lemmas 17.4.1--17.4.2*

**[P.17.2]** The intermediate products involve differences of large terms; subtractive cancellation can occur. Factored forms (§17.8) are essential for numerical precision.

---

### §17.5 The Complete $F_2^*$

This section states the complete second-order averaged Hamiltonian $F_2^*$, separating it into secular and long-period ($\cos 2g'$) parts.

**Theorem 17.5.1** (Brouwer's second-order Hamiltonian)**.** *Stub: the complete $F_2^*$ is*

$$F_2^* = F_2^{*\text{sec}} + F_2^{*\text{lp}}\cos(2g')$$

*where the secular part $F_2^{*\text{sec}}$ depends only on $L, G, H$ and the long-period part $F_2^{*\text{lp}}$ is the coefficient of $\cos(2g')$. This is Brouwer (1959) Eq. 29, cross-verified against Lara (2021) Eq. 14.*

**Theorem 17.5.2** (Secular part of $F_2^*$)**.** *Stub:*

$$F_2^{*\text{sec}} = \frac{\mu^7 R_E^4 J_2^2}{L^7} \cdot \frac{1}{G^7}\left[A(\theta)\,G^2 + B(\theta)\right]$$

*where $A(\theta)$ and $B(\theta)$ are polynomials in $\theta^2 = H^2/G^2$ to be derived. Placeholder form — exact coefficients derived at Develop phase from the Poisson-bracket products of §17.4.*

**Theorem 17.5.3** ($\cos(2g')$ part of $F_2^*$)**.** *Stub: the $\cos(2g')$ coefficient is*

$$F_2^{*\text{lp}} = \frac{\mu^7 R_E^4 J_2^2}{L^7 G^7} C(\theta, \eta)$$

*Drives the long-period eccentricity and inclination variations of Ch 19.*

**Remark (Sign verification).** *Stub: Brouwer (1959) uses positive kinetic energy convention $F = -H_{\text{total}}$; Lara (2021) uses negative total energy. The sign of $F_2^*$ must be verified at Build by tracing both derivations. Cross-reference session note (2026-03-28): the $L'^5/G'^5$ term sign was +15/32, not −15/32.*

---

### §17.6 Second-Order Secular Rates in Delaunay Variables

This section derives the explicit second-order corrections to $\dot{\Omega}$, $\dot{\omega}$, $\dot{M}$ by differentiating the secular part of $F_2^*$ with respect to the Delaunay momenta.

**Theorem 17.6.1** (Second-order secular rate of $h = \Omega$)**.** *Stub:*

$$\dot{h}^{(2)} = \frac{\partial F_2^{*\text{sec}}}{\partial H}$$

*Explicit polynomial in $\theta$ after computing $\partial/\partial H$ via chain rule. Combine with first-order result of Ch 16.* — *Proof approach: partial differentiation of $F_2^{*\text{sec}}(L, G, H)$ with respect to $H$, applying the chain rule through the $\theta = H/G$ dependence*

**Theorem 17.6.2** (Second-order secular rate of $g = \omega$)**.** *Stub:*

$$\dot{g}^{(2)} = \frac{\partial F_2^{*\text{sec}}}{\partial G}$$

*Explicit polynomial in $\theta$ and $\eta$. Identifies the polynomial coefficients $(13, 78, 137)$.* — *Proof approach: partial differentiation of $F_2^{*\text{sec}}$ with respect to $G$, which involves both the explicit $G$-dependence and the implicit dependence through $\eta = G/L$ and $\theta = H/G$; careful accounting of each contribution produces the integer coefficients*

**Theorem 17.6.3** (Second-order secular rate of $l = M$)**.** *Stub:*

$$\dot{l}^{(2)} = \frac{\partial F_2^{*\text{sec}}}{\partial L}$$

*Identifies the polynomial coefficients $(7, 114, 395)$ and $(4, 19)$.* — *Proof approach: partial differentiation of $(F_0 + F_1^* + F_2^{*\text{sec}})$ with respect to $L$, where the dominant contribution comes from the $L^{-7}$ factor in $F_2^{*\text{sec}}$*

**Corollary 17.6.1** (Keplerian form of second-order rates)**.** *Stub: convert Theorems 17.6.1–17.6.3 from Delaunay to Keplerian variables $(a, e, i)$ using Ch 11. Express in terms of $n$, $p = a\eta^2$, $\theta = \cos i$, $\eta$.*

---

### §17.7 Algebraic Origin of Polynomial Coefficients

This section traces each integer coefficient appearing in the second-order secular rates to specific Poisson-bracket products, ensuring no coefficient is accepted without derivation.

**Theorem 17.7.1** (Origin of $(13, 78, 137)$)**.** *Stub: the coefficients in the $\dot{\omega}^{(2)}$ polynomial arise from three groups of Poisson-bracket averages. Show that combining specific lemma results from §17.4 yields $13$, $78$, $137$ as integer-linear combinations. Each coefficient is a sum of at most four integers.* — *Proof approach: collect terms from the $\partial/\partial G$ differentiation of each Poisson-bracket average in Lemma 17.4.3, grouping by powers of $\theta^2$; the integer arithmetic is verified by symbolic computation*

**Theorem 17.7.2** (Origin of $(7, 114, 395)$)**.** *Stub: the $\dot{M}^{(2)}$ polynomial coefficients. Per the session note (2026-03-28): $-35 + 24 + 25 = 14 = 2 \times 7$; $90 - 192 - 126 = -228 = 2 \times (-114)$; $385 + 360 + 45 = 790 = 2 \times 395$; factor of 2 absorbed by $3/32 \to 3/16$. Full derivation at Develop.* — *Proof approach: same method as Theorem 17.7.1 applied to $\partial/\partial L$; the three addends in each coefficient correspond to contributions from $\langle F_2 \rangle$, $\langle \{S_1, F_1\} \rangle$, and $\langle \{S_1, \{S_1, F_0\}\} \rangle$ respectively*

**Theorem 17.7.3** (Origin of $(4, 19)$)**.** *Stub: the small polynomial in the $\dot{\Omega}^{(2)}$ correction. Identify the specific Poisson-bracket product that produces each integer.* — *Proof approach: $\partial/\partial H$ differentiation of $F_2^{*\text{sec}}$; the polynomial $(4 + 19\theta^2)$ arises because $H$ enters only through $\theta = H/G$, reducing the number of contributing terms*

**[A.17.1]** These integer coefficients are exact — they are rational numbers appearing in the perturbation expansion, not approximations. Any numerical discrepancy with a reference is a sign error or transcription error, not a rounding issue. Brouwer (1959) Eq. 29 and Lara (2021) Eq. 14 are the two verification sources.

---

### §17.8 J₄ Contribution to Secular Rates

This section derives the first-order $J_4$ secular rate corrections, which enter at the same order of magnitude as the $J_2^2$ terms.

**Theorem 17.8.1** (First-order J₄ secular rates)**.** *Stub: $J_4$ contributes a first-order secular correction (independent of $J_2$) through the $F_2$ term in the von Zeipel equation. The orbit average $\langle R_4 \rangle$ is evaluated analogously to Ch 16, §16.2. Result:* — *Proof approach: orbit averaging of $J_4(a_E/r)^4 P_4(\sin\phi)/r$ over mean anomaly, using $P_4(\sin\phi)$ expressed in orbital elements (Ch 13, Theorem 13.4.1) and the radial averages $\langle r^{-5} \rangle$, $\langle r^{-5}\cos 2u \rangle$, $\langle r^{-5}\cos 4u \rangle$*

$$\dot{\omega}_{J_4} = \frac{15n R_E^4 J_4}{16p^4}\eta^2\left[\ldots\right]$$

*Placeholder — explicit polynomial in $\theta$ to be derived at Develop.*

**Remark.** *Stub: the $J_4$ contribution and the $J_2^2$ contribution both enter at the same order of magnitude for low Earth orbit; neither can be neglected relative to the other.*

---

### §17.9 $\eta$-Dependent Formulas and the $\eta = 1$ Specialization

This section presents the complete second-order secular rates retaining all $\eta$-dependence, and quantifies the model error introduced by the SGP4 simplification $\eta = 1$.

**Theorem 17.9.1** (Full $\eta$-dependent secular rates)**.** *Stub: complete formulas for $\dot{\Omega}^{(1+2)}$, $\dot{\omega}^{(1+2)}$, $\dot{M}^{(1+2)}$ with all $\eta = \sqrt{1-e^2}$ factors retained. Valid for any eccentricity $e \in [0, 1)$.*

**Corollary 17.9.1** ($\eta = 1$ specialization)**.** *Stub: setting $\eta = 1$ (i.e., $e = 0$) in Theorem 17.9.1 recovers the SGP4 standard formulas. Identify exactly which $\eta$-dependent terms are dropped; quantify the resulting model error [A.17.2] as a function of eccentricity.*

**Example 17.9.1.** *Stub: Molniya orbit with $a = 26560$ km, $e = 0.7$, $i = 63.4°$ using WGS84 ($\mu = 398600.4418$ km$^3$/s$^2$, $a_E = 6378.137$ km, $J_2 = 1.08263 \times 10^{-3}$, $J_4 = -1.62 \times 10^{-6}$). Compute: $\eta = \sqrt{1 - 0.49} = 0.7141$, $p = a\eta^2 = 13539$ km. Full $\eta$-dependent $\dot{\omega}^{(2)}$: evaluate numerically. $\eta = 1$ approximation: same formula with $\eta = 1$. Demonstrate that the discrepancy is $O(e^2) \approx 49\%$ of the full second-order correction. Source: WGS84.*

**[A.17.2]** The SGP4 $\eta = 1$ approximation introduces a model error proportional to $e^2$. For circular orbits this error is negligible; for Molniya-type orbits ($e \approx 0.7$) it is the dominant source of secular rate error. Quantify at Build.

---

### §17.10 Factored Forms Avoiding Cancellation

This section derives algebraically factored forms of the second-order polynomials that avoid subtractive cancellation in floating-point evaluation.

**Theorem 17.10.1** (Factored $Q_\omega$)**.** *Stub: the argument-of-perigee rate polynomial factors as*

$$Q_\omega = (5\theta^2 - 1)(A\theta^2 + B)$$

*or a similar product form that avoids computing a large positive term minus a large positive term. Derive the factoring algebraically.* — *Proof approach: polynomial division and root-finding to express $Q_\omega(\theta^2)$ as a product of irreducible factors, verified by expanding the product back to the original polynomial*

**Theorem 17.10.2** (Factored $Q_\Omega$)**.** *Stub: analogous factored form for the nodal rate polynomial.* — *Proof approach: same polynomial factoring technique as Theorem 17.10.1 applied to $Q_\Omega(\theta^2)$*

**[P.17.3]** The unfactored polynomial $a_0 + a_2\theta^2 + a_4\theta^4$ with mixed-sign coefficients $a_i$ loses digits when $|a_0| \approx |a_2|\theta^2 + |a_4|\theta^4$. The factored form eliminates this cancellation by ensuring all factors are evaluated with the same sign at the current inclination.

---

### §17.11 State Framework

This section expresses the second-order secular corrections as additive perturbations in both state representations, composing with the first-order update of Ch 16.

**§17.11.1 Additive Correction to $({\hat{M}}, {\hat{\Omega}})$ Form**

*Stub: the second-order secular correction is an additive increment to the first-order secular rates, implemented as a near-identity dual quaternion perturbation $\delta\hat{M}^{(2)}$, $\delta\hat{\Omega}^{(2)}$ applied after the first-order update of Ch 16.*

**§17.11.2 Additive Correction in $7\times7$ Form**

*Stub: the second-order rates enter as additional terms in the $\delta W$ matrix of Ch 16, §16.8.2. Composition of first-order and second-order secular updates.*

**§17.11.3 Error Propagation (Principle 1)**

*Stub: the dominant error at second order is model error $\delta_a = O(J_2^3)$ [A.17.1]. Precision error from evaluating $F_2^*$ coefficients [P.17.1]–[P.17.3] is $O(e^2 \epsilon_{\text{mach}})$ in the $\eta = 1$ path.*

---

### §17.12 Generalization: Third-Order $J_2^3$ Terms

This section outlines the structure and magnitude of third-order secular perturbations without full derivation, establishing the error floor of the second-order theory.

*Stub:*

1. *Structure: the von Zeipel order-3 equation involves triple Poisson brackets $\{S_1, \{S_1, \{S_1, F_0\}\}\}$ and cross-terms $\{S_2, F_1\}$, $\{S_1, F_2\}$*
2. *Identify which bracket products enter at order $J_2^3$: at least 10 distinct averaged products*
3. *Magnitude estimate: for LEO ($p \approx 1.1 R_E$, $J_2 \approx 10^{-3}$), the third-order rate is $\sim (J_2/p^2)^3 n \sim 10^{-12} n$*
4. *Satellite classes requiring third-order terms: high-precision orbit determination over multi-week arcs, geodetic satellites (LAGEOS, GRACE-FO)*

**[A.17.3]** Third-order $J_2^3$ secular rates are $\sim (J_2/p^2)^3$ smaller than the Keplerian mean motion. For LEO ($p \approx 1.1 R_E$, $J_2 \approx 10^{-3}$), this is $\sim 10^{-12} n$, negligible for SGP4 but potentially relevant for high-precision orbit determination over weeks.

---

## Cross-References

**Uses (backward)**

| Source | Section | Role |
|---|---|---|
| Ch 2, Theorem 2.4.1 (Velocity dual quaternion $\hat{\Omega}_b$) | §17.11.1 | State Framework dual quaternion form |
| Ch 2, Theorem 2.8.1 ($7\times7$ matrix equivalence) | §17.11.2 | State Framework matrix form |
| Ch 4, Theorem 4.9.1 (Factored polynomial evaluation) | §17.10 | Factored forms avoiding cancellation |
| Ch 11, Delaunay partial derivatives | §17.6 | Rate conversions |
| Ch 12, von Zeipel method and Deprit triangle | §17.3 | Theoretical framework |
| Ch 16, Averaged disturbing function $F_1$ | §17.2–17.3 | Input to von Zeipel order-2 equation |
| Brouwer (1959) Eqs. 28–29, 39–41 | §17.5 | Cross-verification of $F_2^*$ |
| Lara (2021) Eq. 14 | §17.5 | Cross-verification of $F_2^*$ |

**Feeds (forward)**

| Target | Section | Role |
|---|---|---|
| Ch 18, Short-period corrections | §17.2 | $S_1$ constructed here used in short-period generating function |
| Ch 19, Long-period corrections | §17.5 | $\cos(2g')$ part of $F_2^*$ drives long-period corrections |
| Ch 20, Osculating elements | §17.6 | Second-order secular rates enter mean-to-osculating transformation |

---

## Error Notes

| Tag | Type | § | Description |
|---|---|---|---|
| [A.17.1] | A | §17.12 | Second-order theory is exact within $O(J_2^2)$; next error is $O(J_2^3)$ |
| [A.17.2] | A | §17.9 | SGP4 $\eta = 1$ approximation; error $\propto e^2$, significant for $e > 0.3$ |
| [A.17.3] | A | §17.12 | Third-order $J_2^3$ terms; negligible for SGP4 applications |
| [P.17.1] | P | §17.2 | Precision error in $S_1$ near $e = 0$ (eccentricity singularity) and near $i = i_c$ (critical inclination denominator) |
| [P.17.2] | P | §17.4 | Subtractive cancellation in Poisson-bracket products; magnitude quantified for each product |
| [P.17.3] | P | §17.10 | Unfactored polynomial evaluation; factored forms required for reliable computation |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 1 |
| Theorems | 15 |
| Lemmas | 3 |
| Corollaries | 3 |
| Propositions | 0 |
| Examples | 1 |
| Error Notes | 6 |
| Equations | ~20 |
| Sections | 12 |

**Tier classification:** $F_2^*$ secular part — Tier I (exact within perturbation order, integer coefficients). $\eta = 1$ simplification — Tier III (matched-pair SGP4 approximation, Chapter 3 applies). $J_4$ contribution — Tier II (orbit-averaged, same rigour as first-order J₂).

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §17.1 Introduction | Draft | |
| §17.2 First-Order Generating Function $S_1$ | Draft | |
| §17.3 Von Zeipel Order-2 Equation | Draft | |
| §17.4 Orbit-Averaged Poisson-Bracket Products | Draft | |
| §17.5 The Complete $F_2^*$ | Draft | |
| §17.6 Second-Order Secular Rates in Delaunay Variables | Draft | |
| §17.7 Algebraic Origin of Polynomial Coefficients | Draft | |
| §17.8 J₄ Contribution to Secular Rates | Draft | |
| §17.9 $\eta$-Dependent Formulas and $\eta = 1$ Specialization | Draft | |
| §17.10 Factored Forms Avoiding Cancellation | Draft | |
| §17.11 State Framework | Draft | |
| §17.12 Generalization: Third-Order $J_2^3$ Terms | Draft | |
