# Analytic Remediation of the Hansen Library for Second-Order Celestial Mechanics

**Author:** Critical Mathematical Consultant  
**Subject:** Formal Rectification of Step 05a (Hansen Extension)  
**Constraint:** False Confidence Protocol (Zero-Correctness-Assertion)

## I. Symplectic Context and Global Assumptions

We consider a Hamiltonian system on a six-dimensional symplectic manifold $(\mathcal{M}, \omega)$, where $\mathcal{M} = T^* \mathbb{R}^3 \cong \mathbb{R}^3 \times \mathbb{T}^3$ is parameterized by the Delaunay action-angle coordinates $(L, G, H, l, g, h)$. We define the regular domain $\mathcal{D} \subset \mathcal{M}$ as the set of points where the eccentricity $e := \sqrt{1 - (G/L)^2}$ satisfies $e \in [0, 1)$ and the inclination $I := \arccos(H/G)$ is well-defined.

The perturbation theory developed by Brouwer (1959) and refined by Hori (1961) requires the evaluation of the second-order secular Hamiltonian $F_2^*$, which depends on the harmonic average $\mathcal{I}_g := \langle \{F_1, S_1\} \rangle_l$. This average necessitates the construction of an analytic extension of the Hansen coefficients $X_k^{n,m}(e)$ to include kernels $n \in \{-3, -4, -5\}$. 

We assume the following project primitives as established truths:
1.  **Keplerian measure:** The volume form on the torus $\mathbb{T}^1$ transforms as $dl = \kappa dE$, where $\kappa = 1 - e \cos E$.
2.  **Conformal parameterization:** The mapping from the eccentric anomaly $E$ to the true anomaly $f$ is given by the Möbius transformation $e^{if} = (u-\alpha)/(1-\alpha u)$ on the complex unit circle $\Gamma = \{u \in \mathbb{C} : |u|=1\}$, where $\alpha = e/(1+\eta)$ and $\eta = \sqrt{1-e^2}$.
3.  **Homological Identity:** The generator $S_1$ satisfies $n \frac{\partial S_1}{\partial l} = F_1 - F_1^*$, where $n = \mu^2/L^3$.

---

## II. Measure Transformations and the Residue Prefactor (Lemma 6.4.1)

To evaluate the Hansen integral $X_k^{n,m} = \frac{1}{2\pi} \int \kappa^n e^{i(mf-kl)} dl$, we must lift the integration to the complex plane. We demonstrate that for the kernel $n=-3$, the transformation introduces a specific scalar weight.

**Proposition 1.** Let $n=-3$. The transformation of the Hansen integral to a contour integral over $\Gamma$ yields a scalar prefactor $\zeta = 4/(1+\eta)^2$.

**Proof.** 
By the Keplerian relation, $dl = \kappa dE$. Substituting this into the integral definition:
$$ X_k^{n,m} = \frac{1}{2\pi} \int_{0}^{2\pi} \kappa^{n+1} e^{i(mf-kl)} dE. $$
For $n=-3$, the kernel becomes $\kappa^{-2}$. We parameterize the circle by $u = e^{iE}$, so $dE = du/(iu)$. The radius $\kappa(u)$ is expressed as:
$$ \kappa(u) = 1 - \frac{e}{2}(u + u^{-1}) = \frac{1+\eta}{2u} (u-\alpha)(1-\alpha u). $$
Taking the reciprocal square:
$$ \kappa^{-2}(u) = \frac{4u^2}{(1+\eta)^2} (u-\alpha)^{-2} (1-\alpha u)^{-2}. $$
The integral on $\Gamma$ is:
$$ X_k^{-3,m} = \frac{1}{2\pi} \oint_{\Gamma} \left[ \frac{4u^2}{(1+\eta)^2} (u-\alpha)^{-2} (1-\alpha u)^{-2} \right] e^{i(mf-kl)} \frac{du}{iu}. $$
Rearranging constants and $u$ powers:
$$ X_k^{-3,m} = \frac{4}{(1+\eta)^2} \left[ \frac{1}{2\pi i} \oint_{\Gamma} \frac{u e^{i(mf-kl)}}{(u-\alpha)^2 (1-\alpha u)^2} du \right]. $$
The bracketed expression is the residue $\Phi_{n,k,m}$. The identity $\zeta = 4/(1+\eta)^2$ is thus established. $\square$

---

## III. Mapping the $\mathcal{L}$-residual onto the Hansen Basis (§6.0)

The failed derivation erroneously used the Bessel expansion for the equation of center. We provide the correct mapping based on the analytic properties of the $\mathcal{L}$-residual $\mathcal{A} := f - l + e \sin f$.

**Proposition 2.** The $\mathcal{L}$-residual interaction in Route $\alpha$ is governed by the $X_k^{-3, 0}$ Hansen basis.

**Proof.** 
From `ch07c` Proposition C.4.8.3, we have the exact Fourier expansion:
$$ \mathcal{A}(l, e) = 2 \eta^3 \sum_{k=1}^{\infty} \frac{X_k^{-3, 0}(e)}{k} \sin(kl). $$
The contribution to $\mathcal{I}_g$ involves the partial derivative $\frac{\partial S_1}{\partial G}$. Let $S_{1, \mathcal{A}} = C_{\mathcal{A}}(\theta, \eta) \mathcal{A}(l, e)$, where $C_{\mathcal{A}} = \frac{\mu k_2 A(\theta)}{n a^3 \eta^3}$. 
We differentiate $C_{\mathcal{A}}$ with respect to the action $G$. Given $\theta = H/G$ and $\eta = G/L$:
$$ \frac{\partial \theta}{\partial G} = -\frac{\theta}{G}, \quad \frac{\partial \eta^{-3}}{\partial G} = -\frac{3}{G \eta^3}. $$
The derivative of the product $A(\theta) \eta^{-3}$ is:
$$ \frac{\partial (A \eta^{-3})}{\partial G} = \left( 3\theta \frac{\partial \theta}{\partial G} \right) \eta^{-3} + A \left( -\frac{3}{G \eta^3} \right) = \frac{-3\theta^2}{G \eta^3} - \frac{3(3/2\theta^2 - 1/2)}{G \eta^3} = -\frac{3(5\theta^2 - 1)}{2 G \eta^3}. $$
The resulting secular weight for the interaction is:
$$ \frac{\partial C_{\mathcal{A}}}{\partial G} = -\frac{3 \mu k_2 (5\theta^2 - 1)}{2 G \eta^3 n a^3}. $$
This confirms that the interaction with the Hamiltonian's $\sin(2f+2g)$ harmonics must be expanded in the $X_k^{-3, 0}$ basis to maintain structural integrity. $\square$

---

## IV. Trigonometric Closure and the $\cos(2g)$ Mode (Theorem 6.4.5)

The previous agent provided a "sketch" of Theorem 6.4.5. We provide the logical reduction for the harmonic collapse.

**Proposition 3.** The $l$-average of the interaction between the $\mathcal{L}$-residual and the Hamiltonian's primary harmonic reduces to a pure $\cos(2g)$ mode.

**Proof.** 
The relevant integrand is proportional to $\mathcal{A}(l, e) \sin(2f+2g)$. Substituting the Hansen expansion for $\mathcal{A}$ and the complex exponential form for the sine:
$$ \langle \sin(kl) \sin(2f+2g) \rangle_l = \frac{1}{2} \langle \cos(kl - 2f - 2g) - \cos(kl + 2f + 2g) \rangle_l. $$
Expanding the arguments:
$$ \cos(kl - 2f - 2g) = \cos(kl - 2f)\cos(2g) + \sin(kl - 2f)\sin(2g). $$
Under the $l$-average, the sine terms vanish by parity. The cosine terms yield:
$$ \langle \cos(kl - 2f) \rangle_l = X_k^{-p, 2} \quad \text{(Hansen definition)}. $$
Summing over $k$, the interaction collapses into the form $\Lambda_{\mathcal{A}} \cos(2g)$, where $\Lambda_{\mathcal{A}}$ is a series in $X_k^{-3, 0} X_k^{-3, \pm 2}$. $\square$

---

## V. Remediation Directives for the Principal Author

1.  **Residue Convention:** Fix $\Phi_{n,k,m}$ such that the prefactor is exactly $4/(1+\eta)^2$. Eliminate any drift to $2/\pi i$.
2.  **Basis Correction:** Replace all $J_k$ Bessel coefficients for $\mathcal{A}$ with the correct $X_k^{-3, 0}$ Hansen coefficients.
3.  **Algebraic Brute-Force:** In the Section 7 extensions, use the identity $1 - 3\eta^2 + 2\eta^3 = (1-\eta)^2(1+2\eta)$ to prove the equivalence of the unified $X_0^{0,2}$ form.
4.  **Parity Audit:** Correct the citation in §9. The $\eta$-denominator power is $\eta^{2p-3}$ (derived from the Wallis integral $I_{p-1}$), not $\eta^{2p+1}$.

**Path to Analysis:** `design/cleanroom/consultant_recommendations/step_05ar_critical_analysis.md`
