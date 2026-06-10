# Theorem: δp₁ = D(∂S₁/∂l) = 6ΓB₀(a/r)

## Definitions

- (ξ_j, η_j) = (dx_j/dt, x_j): Cartesian velocities and positions
- (L_j, l_j): Delaunay canonical variables
- D = -Σ_i ξ_i ∂/∂ξ_i: Euler velocity-homogeneity operator (BH61 Eq 10)
- S₁: Brouwer's short-period generating function (verified Phase A, 72/72)
- Γ = μ²/(a³η³), B₀ = -½ + ³⁄₂θ², B₁' = ³⁄₂(1-θ²)
- p_j = Σ_k ξ_k ∂η_k/∂l_j, q_j = Σ_k ξ_k ∂η_k/∂L_j (BH61 Eq 4)

---

## Lemma 1: D annihilates position functions

**Statement:** D(η_k) = 0 for all k. Consequently D(f(η)) = 0 for any function f of positions only.

**Proof:** D(η_k) = -Σ_j ξ_j ∂η_k/∂ξ_j = 0, because η_k = x_k does not depend on ξ_j = dx_j/dt. The variables (ξ, η) are independent coordinates in phase space. For f(η), the chain rule gives D(f) = Σ_k (∂f/∂η_k) D(η_k) = 0. □

**Corollary:** Dr = 0, since r = |η| = √(Σ η_k²) is a pure position function.

---

## Lemma 2: D-action on Delaunay variables

**Statement (BH61 Eq 12):** DL_j = -(p_j)'', Dl_j = +(q_j)''

**Proof:** The canonical transformation (ξ,η) → (L,l) satisfies ∂L_j/∂ξ_k = ∂η_k/∂l_j (a property of canonical transforms — the generating function identity). Therefore:

DL_j = -Σ_k ξ_k ∂L_j/∂ξ_k = -Σ_k ξ_k ∂η_k/∂l_j = -p_j

Similarly, using ∂l_j/∂ξ_k = -∂η_k/∂L_j:

Dl_j = -Σ_k ξ_k ∂l_j/∂ξ_k = +Σ_k ξ_k ∂η_k/∂L_j = +q_j □

---

## Lemma 3: Dθ = 0 (inclination invariance)

**Statement:** D(θ) = 0 where θ = H/G = cos I.

**Proof:** DH = -p₃ = -H and DG = -p₂ = -G (Lemma 2 with p₃ = H, p₂ = G from BH61 Eq 5).

D(H/G) = (G·DH - H·DG)/G² = (G(-H) - H(-G))/G² = (-GH + GH)/G² = 0. □

**Corollary:** D(B₀) = D(B₁') = 0, since B₀ = -½+³⁄₂θ² and B₁' = ³⁄₂(1-θ²) depend only on θ.

---

## Lemma 4: Df = 2sinf/e

**Statement:** The D-action on the true anomaly is Df = 2sinf/e.

**Proof:** Since f = f(L, G, l) (depends on Delaunay momenta and mean anomaly through Kepler's equation, but not on g, h, H):

Df = -p₁ ∂f/∂L - p₂ ∂f/∂G + q₁ ∂f/∂l

The partial derivatives of f:
- ∂f/∂l|_{L,G} = a²η/r² (standard Kepler relation)
- ∂f/∂L|_{G,l} = sinf(2+ecosf)/(eL) (implicit differentiation through Kepler; derived in Phase A2, verified numerically)
- ∂f/∂G|_{L,l} = -sinf(2+ecosf)/(eG) (similarly derived)

Substituting p₁ = L(2a/r-1), p₂ = G, q₁ = 2esinE + 2ηsinf/e:

Df = -L(2a/r-1)·sinf(2+ecosf)/(eL) - G·(-sinf(2+ecosf)/(eG)) + (2esinE + 2ηsinf/e)·a²η/r²

Simplifying the first two terms, with ρ = a/r:
= sinf(2+ecosf)/e · [-(2ρ-1) + 1] + (2esinE + 2ηsinf/e)·a²η/r²
= sinf(2+ecosf)/e · (2-2ρ) + (2esinE + 2ηsinf/e)·ηρ²

Using 2+ecosf = 1+η²ρ (since ecosf = η²ρ-1):
= -2sinf(1+η²ρ)(ρ-1)/e + 2esinE·ηρ²/1 + 2η²ρ²sinf/e  ... [*]

For the sinE term, use sinE = rsinf/(aη) (from the relation aηsinE = rsinf):
2e·a²η·sinE/r² = 2e·a²η·rsinf/(aη·r²) = 2easinf/r = 2eρsinf

So [*] becomes:
= -2sinf(1+η²ρ)(ρ-1)/e + 2eρsinf + 2η²ρ²sinf/e

= (2sinf/e)·[-(1+η²ρ)(ρ-1) + e²ρ + η²ρ²]

Expanding -(1+η²ρ)(ρ-1):
= -(ρ-1) - η²ρ(ρ-1) = -ρ+1-η²ρ²+η²ρ

Adding e²ρ + η²ρ²:
= -ρ+1-η²ρ²+η²ρ + e²ρ + η²ρ²
= 1 + ρ(-1+η²+e²)

**Since η² + e² = 1:**
= 1 + ρ·0 = 1

Therefore: **Df = (2sinf/e)·1 = 2sinf/e.** □

---

## Lemma 5: D(f+g) = 0

**Proof:** Dg = q₂ = -2sinf/e (BH61 Eq 5, definition).
D(f+g) = Df + Dg = 2sinf/e + (-2sinf/e) = 0. □

**Corollary:** D(cos(2f+2g)) = -sin(2f+2g)·2·D(f+g) = 0.

---

## Lemma 6: D(Γ) = 6Γ(a/r)

**Proof:** Γ = μ²·a⁻³·η⁻³.

Da = -2a(2a/r-1) (from Lemma 2: Da = (2a/L)DL = (2a/L)(-p₁) = -2a(2a/r-1)).
Dη = 2η(a/r-1) (from Dη = (1/L)DG + (-G/L²)DL = (-G/L) + (G/L)(2a/r-1) = η(2a/r-2) = 2η(a/r-1)). [Note: Dη uses DG = -G and DL = -p₁ = -L(2a/r-1), via η = G/L and the quotient rule.]

D(a⁻³) = -3a⁻⁴·Da = -3a⁻⁴·(-2a)(2a/r-1) = 6a⁻³(2a/r-1)

D(η⁻³) = -3η⁻⁴·Dη = -3η⁻⁴·2η(a/r-1) = -6η⁻³(a/r-1)

D(a⁻³η⁻³) = η⁻³·D(a⁻³) + a⁻³·D(η⁻³)
= 6a⁻³η⁻³(2a/r-1) - 6a⁻³η⁻³(a/r-1)
= 6a⁻³η⁻³[(2a/r-1)-(a/r-1)]
= 6a⁻³η⁻³·(a/r)

D(Γ) = μ²·D(a⁻³η⁻³) = 6μ²a⁻³η⁻³·(a/r) = 6Γ(a/r). □

---

## Main Theorem

**Statement:** δp₁ = D(∂S₁/∂l) = 6Γ B₀ (a/r)

**Proof:**

From the homological equation (verified Phase A, 72/72 quadrature tests):

$$\frac{\partial S_1}{\partial l} = \underbrace{\Gamma B_0}_{\text{Term 1}} - \underbrace{\frac{\mu^2}{r^3}\left[B_0 + B_1'\cos(2f+2g)\right]}_{\text{Term 2}}$$

Apply D to each term using the product rule (D is a derivation):

**D(Term 1) = D(ΓB₀):**
= D(Γ)·B₀ + Γ·D(B₀)
= 6Γ(a/r)·B₀ + 0    [Lemmas 6 and 3]
= **6ΓB₀(a/r)**

**D(Term 2) = D(μ²r⁻³[B₀ + B₁'cos(2f+2g)]):**
= μ²·{D(r⁻³)·[B₀+B₁'cos(2f+2g)] + r⁻³·D(B₀) + r⁻³·D(B₁')·cos(2f+2g) + r⁻³·B₁'·D(cos(2f+2g))}

Each factor:
- D(r⁻³) = -3r⁻⁴·Dr = 0 [Lemma 1 corollary]
- D(B₀) = 0 [Lemma 3 corollary]
- D(B₁') = 0 [Lemma 3 corollary]
- D(cos(2f+2g)) = 0 [Lemma 5 corollary]

Therefore D(Term 2) = **0**.

**Result:**

$$\delta p_1 = D\!\left(\frac{\partial S_1}{\partial l}\right) = 6\,\Gamma\, B_0\,\frac{a}{r} - 0 = \frac{6\mu^2 B_0}{a^2 r\,\eta^3}$$

□

---

## Discussion

### Scope: This theorem applies to OUR normalization only

**IMPORTANT (2026-04-05 resolution):** This theorem computes D of OUR form of ∂S₁/∂l (using Γ = μ²/(a³η³)). The result 6ΓB₀(a/r) is correct FOR THIS FUNCTION.

However, BH61 Eq(14) computes D of BROUWER'S form of ∂S₁/∂l (using n = μ²/L³). These are DIFFERENT functions — they differ by a velocity-dependent factor -(μ/a)^(3/2), which equals -1 only when a=μ (the normalized case used in all prior tests). See `verify_normalization_ratio.py`.

### BH61 Eq(14) structural form verified — the cos(2f+2g) terms are present in D(Brouwer's form)

BH61 Eq(14) contains cos(2f+2g) terms that do not appear in our result. These are NOT errors. They arise from D(ρ³) in the Leibniz expansion of D(n·(content)):

- D(cos(2f+2g)) = 0 (Lemma 5) — everyone agrees
- But D((μ²/L³)·ρ³·cos(2f+2g)) ≠ 0, because D(μ��/L³) ≠ 0 and D(ρ³) ≠ 0
- The cos(2f+2g) harmonic survives through multiplication by velocity-dependent coefficients

In our form, the same content is grouped as μ²/r³·[...], and D(μ²/r³) = 0 because r is a position function. The grouping matters because D is a derivation: D(fg) = D(f)g + fD(g), and different factorizations distribute the velocity-dependence differently.

**Numerically confirmed:** D(Brouwer's form) = BH61 Eq(14) at 81/81 test points (max err 1.85e-05). See `verify_D_brouwer_form.py`, `verify_breakthrough_algebra.py`.

### INPE-2746 "spurious Poisson terms" — not about Eq(14)

INPE-2746-PRE/322 (Fitzgibbon, De Moraes, Lobão, 1983) claims BH61 contains "spurious Poisson terms." Reading the paper (OCR translation from Portuguese) reveals:

1. The paper never mentions Eq(14), the D operator, δp₁, or cos(2f+2g) terms specifically.
2. Section 3 states: "This theory contains spurious Poisson terms [5]. To use it, all these terms must be neglected [6]" — citing Fitzgibbon's thesis [5] and Vilhena de Moraes (1981) [6]. The claim is not derived in the paper itself.
3. Section 4 ("Elimination of Spurious Terms") describes their alternative: Variation of Parameters with Bessel function expansions. Their solution "does not contain Poisson terms."
4. Section 5 separately mentions "spurious **secular** terms" alongside the Taylor series convergence problem.

The paper does not explain what the spurious terms are, where they appear, or why they arise. The underlying analysis is in Fitzgibbon's thesis [5] and Vilhena de Moraes (1981) [6], which we do not have.

What we can say: INPE-2746 does not support or contradict our normalization resolution regarding Eq(14). It is about the full BH61 theory at a level of detail the conference paper does not specify.

### Verification status

The theorem has been verified numerically: 81/81 test cases pass with max relative error 5.09×10⁻⁶ (finite-difference limited). See `cleanroom_phase_a3_analysis.py`. This confirms D(our ∂S₁/∂l) = 6ΓB₀(a/r), which is correct but distinct from BH61 Eq(14).
