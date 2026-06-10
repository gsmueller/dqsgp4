# Phase B Proof Workbook — Deriving BH61 Eq(14) as a Theorem

## Mission

Prove BH61 Eq(14) for δp₁ algebraically. No numerical verification until the proof is complete. The proof should be readable as a textbook derivation.

## What we need to prove

BH61 Eq(14) claims:

$$\delta p_1 = 3\frac{\mu^2 k_2}{L^3}\left\{B_0\left[-\eta^{-3} + \frac{a^3}{r^3}\left(1 - 2\frac{a}{r}\right)\right] + B_1\frac{a^3}{r^3}\left(1 - 2\frac{a}{r}\right)\cos(2g+2f)\right\}$$

where δp₁ = D(∂S₁/∂l) and D = -Σ ξ_i ∂/∂ξ_i is the velocity-homogeneity operator.

## Strategy

### Part I: Establish the D operator framework
- Define D from the Cartesian velocity-homogeneity
- Derive its action on Delaunay variables (get p_j, q_j)
- Derive its action on orbital elements (a, e, η, θ, r, f, g)
- Establish key lemma: what is D's action on products and compositions?

### Part II: Understand ∂S₁/∂l
- Write ∂S₁/∂l in the most useful form for applying D
- Identify what D acts on: which factors depend on velocities?

### Part III: The core proof
- Apply D to ∂S₁/∂l term by term
- Use the homogeneity properties to simplify
- Arrive at Eq(14) or identify where BH61's form comes from

### Part IV: Numerical verification
- Only AFTER the proof is complete, verify numerically

## Key insight I'm missing

The D operator is the VELOCITY Euler operator. For a function homogeneous of degree n in velocities: D(f) = -nf. This is the homogeneity theorem.

BH61 Eq(10) states: "If f(L_j, l_j) is a homogeneous function of degree n in the ξ's only, Df = -nf."

So the question becomes: **what is the velocity-homogeneity degree of each factor in ∂S₁/∂l?**

If ∂S₁/∂l = Σ (terms), and each term has a definite homogeneity degree in the velocities, then D acts by simply multiplying by -n.

This is probably how BH61 gets their result — not by chain-ruling through all the orbital elements, but by identifying the homogeneity degree of each piece.

## Homogeneity analysis needed

In Keplerian motion:
- ξ_i = dx_i/dt (velocities) scale like v
- η_i = x_i (positions) scale like... well, they don't scale with v

If we scale all velocities by λ (ξ → λξ) while keeping positions fixed:
- v → λv
- What happens to a, e, r, f, l, g, L, G, H?

This is the key question. The vis-viva relation v² = μ(2/r - 1/a) determines a from v and r:
- If v → λv: v² → λ²v², so μ(2/r - 1/a) → λ²μ(2/r - 1/a)
- Since r is fixed (position): 2/r - 1/a → λ²(2/r - 1/a)
- So 1/a → 2/r - λ²(2/r - 1/a) = 2(1-λ²)/r + λ²/a
- Hmm, this doesn't give a clean scaling for a.

Wait — I think the homogeneity is about scaling velocities at FIXED positions in the KEPLER problem. Let me think about this more carefully.

Actually, BH61 defines (Eq 1): ξ_j = dx_j/dt (velocities), η_j = x_j (positions), F = -(1/2)Σξ_j² + U.

The canonical momenta in this (ξ,η) system ARE the velocities. So L_j are functions of (ξ, η) through the canonical transform. The D operator scales ξ while holding η fixed.

For Kepler motion at a given position r:
- L = √(μa), and a = μr/(2μ - rv²), so a depends on v² and r
- G = |r × v|, so G depends on |v| and the angle
- f is determined by the position r in the orbit, which depends on (a, e) and the direction

This is getting complicated. I need to work through the homogeneity degrees systematically.

## TODO for the proof

1. [ ] Determine homogeneity degree of each Delaunay variable in ξ
2. [ ] Determine homogeneity degree of each orbital element in ξ
3. [ ] Determine homogeneity degree of Γ = μ²/(a³η³) in ξ
4. [ ] Determine homogeneity degree of 1/r³ in ξ (should be 0 since r is position)
5. [ ] Determine homogeneity degree of cos(2f+2g) in ξ
6. [ ] Apply the homogeneity theorem to each factor in ∂S₁/∂l
7. [ ] Assemble the result
8. [ ] Compare to BH61 Eq(14)

## Session log

(Will be updated as work progresses)
