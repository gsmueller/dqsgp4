# Phase B Proof Synthesis Plan

## What the agents are computing

### Agent 0 (Phase B spec execution)
- Full cleanroom derivation following BH61's assumptions
- Numerical verification of the result
- Output: `cleanroom_phase_b_results.md`

### Agent 1 (Homogeneity analysis)
- Key question: are Keplerian quantities homogeneous in velocity?
- If YES: D(f) = -nf gives the answer directly
- If NO: must use the chain-rule approach through Delaunay variables
- This determines which proof technique is correct

### Agent 2 (BH61 proof structure)
- Complete proof that δp_j = D(∂S₁/∂l_j) from the Jacobian definition
- Analysis of BH61 Eq(10)-(13) properties
- Understanding of WHY Dr = 0, Dβ = 0

### Agent 3 (Algebraic D(∂S₁/∂l))
- Step-by-step algebra applying D to each term
- Identification of which terms survive
- Comparison to BH61 Eq(14) form

## The critical question for synthesis

The homogeneity analysis (Agent 1) will determine which of two proof paths is correct:

### Path A: Direct homogeneity
If ∂S₁/∂l or its components are homogeneous of definite degree n in velocities,
then D(∂S₁/∂l) = -n · ∂S₁/∂l (or applied term by term).
This would immediately give the structure of Eq(14).

### Path B: Chain rule through Delaunay
If the components are NOT homogeneous (likely, since 1/a = 2/r - v²/μ is affine, not power-law),
then we must use D(F) = Σ[-p_j ∂F/∂L_j + q_j ∂F/∂l_j] and work through the algebra.
This is what we did in Phase A3 (getting 6ΓB₀(a/r)), but possibly with errors.

### Path C: Hybrid
Perhaps ∂S₁/∂l can be decomposed into pieces that ARE individually homogeneous,
even though the orbital elements themselves are not.

For example: μ²/r³ is degree 0 in velocities (pure position function).
And Γ = μ²/(a³η³) has some definite degree.
If Γ has degree n, then D[Γ] = -nΓ, and D[Γ B₀] = -nΓ B₀ (since B₀ depends on θ = cosI, and the homogeneity of θ needs to be determined).

## What I expect to find

The key identity is probably:
- Γ = μ²/(a³η³) is homogeneous of degree -6 in velocities (since a has degree -2 from vis-viva, η has degree... needs analysis)
- μ²/r³ is degree 0 (position only)
- So D[Γ] = 6Γ, D[μ²/r³] = 0
- But then D[∂S₁/∂l] = D[ΓB₀] - D[μ²(B₀+B₁'cos(2f+2g))/r³] = 6ΓB₀ - 0 = 6ΓB₀

Wait — that gives a CONSTANT (no a/r factor), not 6ΓB₀(a/r).

Unless Γ is NOT homogeneous of degree -6. Let me think...

From vis-viva: v² = μ(2/r - 1/a), so 1/a = 2/r - v²/μ.
When v → λv at fixed r: 1/a → 2/r - λ²v²/μ = 2/r(1) - λ²·(2/r - 1/a) + ...

This shows a is NOT a homogeneous function of v. So Γ is NOT homogeneous either.

Therefore Path A fails. We must use Path B or Path C.

But BH61 claims D(f) = -nf for homogeneous functions, and then evaluates D(∂S₁/∂l) to get Eq(14). They must be decomposing ∂S₁/∂l into pieces that ARE homogeneous.

What pieces could be homogeneous?
- p₁ = L(2a/r - 1): from BH61's analysis, p_j = Σ ξ_k ∂η_k/∂l_j. Since ξ_k has degree 1 and η_k has degree 0, and ∂η_k/∂l_j depends on the canonical transform (which mixes positions and velocities)... this needs the Agent 1 analysis.

## After agents complete

1. Read all agent outputs
2. Determine which proof path is correct
3. Assemble the proof
4. Write numerical verification
5. Update `cleanroom_phase_b_results.md`
6. Update memory notes
