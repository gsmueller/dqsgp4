# Phase B Synthesis: Resolving Agent Results

## The Critical Error Found

**Agent 3 used incorrect D-actions on Delaunay variables.**

Agent 3 assumed:
- DL = -L, Dl = l (simple homogeneity scaling)
- This gives Da = -2a, De = 0, Dη = 0

But BH61 Eq(12) actually says:
- DL = -(p₁)'' = -L(2a/r - 1)
- Dl = +(q₁)'' = 2esinE + 2ηsinf/e

The p_j, q_j are NOT the Delaunay variables themselves. They are specific functions defined by Eq(5).

Specifically:
- p₁ = L(2a/r - 1) ≠ L  (differs by the factor 2a/r - 1)
- p₂ = G = L₂  (happens to equal the momentum)
- p₃ = H = L₃  (happens to equal the momentum)
- q₁ = 2esinE + 2ηsinf/e  (NOT equal to l)
- q₂ = -2sinf/e  (NOT equal to g)
- q₃ = 0  (NOT equal to h)

This means Agent 3's entire derivation is invalid because every D-action it computed (Da, De, Dη, Df, Dr, etc.) used the wrong starting identities.

## The Correct D-Actions (from my earlier numerical work)

Using DL = -L(2a/r-1), DG = -G, DH = -H, Dl = 2esinE+2ηsinf/e, Dg = -2sinf/e, Dh = 0:

| Quantity | D action | How derived |
|----------|----------|-------------|
| Da | -2a(2a/r - 1) | Da = (2a/L)DL = (2a/L)(-L(2a/r-1)) |
| De | -2(e + cosf) | Chain rule from De² = -D(G²/L²), using DG=-G, DL=-L(2a/r-1) |
| Dη | 2η(a/r - 1) | From DG=-G, DL=-L(2a/r-1), quotient rule on η=G/L |
| Dθ | 0 | From DH=-H, DG=-G, quotient rule on θ=H/G |
| Dr | 0 | Numerically verified to 1e-7; follows from Lemma 1 (r is position) |
| Df | +2sinf/e | Numerically verified to 1e-9; derived via chain through Kepler |
| Dg | -2sinf/e | Definition: Dg = q₂ = -2sinf/e |
| D(f+g) | 0 | Df + Dg = 0 |

These were verified numerically in Phase A3 (81/81 PASS for the resulting D(∂S₁/∂l) = 6ΓB₀(a/r)).

## What Agent 1 Tells Us

Agent 1 correctly proved that most Keplerian quantities are NOT velocity-homogeneous:
- a, e, L, l, f, g, E, Γ: all non-homogeneous
- G, H, cosI: homogeneous (degrees 1, 1, 0)

This means the shortcut D(f) = -nf cannot be used for these quantities. The D-actions must be computed via the chain rule from the fundamental DL_j = -p_j, Dl_j = q_j.

This is CONSISTENT with our numerically verified results: Da = -2a(2a/r-1) is NOT equal to any constant times a.

## What Agent 2 Tells Us

Agent 2 correctly proved the D-operator reduction:
- Starting from the Jacobian definition of δp_j (Eq 8')
- Through substitution of the generating function relations (Eq 9)
- Through conversion to Cartesian partials
- Through the fundamental Poisson bracket identity {η_i, ξ_p} = -δ_ip
- Arriving at δp_j = D(∂(S₁+S₁*)/∂l_j'')

This proof is valid and complete. The theorem δp_j = D(∂S₁/∂l_j) is established.

## The Remaining Question

The theorem says δp₁ = D(∂S₁/∂l). Our proof in theorem_dp1.md derives this as 6ΓB₀(a/r), which was verified 81/81. But BH61 Eq(14) claims a different (more complex) result.

The proof in theorem_dp1.md uses:
1. D(ΓB₀) = 6Γ(a/r)B₀ (Lemma 6 — chain rule through Da and Dη)
2. D(μ²[B₀+B₁'cos(2f+2g)]/r³) = 0 (from Dr=0, DB₀=0, DB₁'=0, D[cos(2f+2g)]=0)

The key step (2) depends on D[cos(2f+2g)] = 0, which follows from D(f+g) = 0.

**If any of the D-action identities are wrong, step (2) could be wrong.**

The D-actions were verified NUMERICALLY. But Agent 3's algebraic attempt (with wrong assumptions) got completely different results — showing that algebraic errors in the D-actions propagate catastrophically.

## Path Forward

The proof stands or falls on the D-action identities. These were verified numerically but not proven algebraically in a single complete derivation from first principles. The Agent 3 attempt shows how easy it is to get them wrong algebraically.

What I need:
1. A complete algebraic proof of Df = 2sinf/e from DL=-p₁, Dl=q₁ (my theorem_dp1.md Lemma 4 has this but should be checked)
2. Numerical verification that specifically tests D[cos(2f+2g)] = 0 at a grid of points
3. A reconciliation: WHY does BH61 get a different result? Did they use incorrect D-actions? Did they compute D(∂S₁/∂l) differently?
