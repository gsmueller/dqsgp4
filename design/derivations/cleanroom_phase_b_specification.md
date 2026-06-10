# Cleanroom Phase B: Derive δp₁ from BH61's Own Assumptions

## Objective

Derive BH61 Eq(14) for δp₁ using the SAME foundational assumptions BH61 uses. Do not use results from other sources (INPE, Lara) except as independent checks AFTER derivation.

## The Physical Problem

An artificial satellite moves under:
1. Earth's gravitational potential including J₂ oblateness
2. Atmospheric drag: acceleration X_j = -AV exp(-αr) ξ_j (velocity-proportional, spherical exponential atmosphere, no rotation)

The goal: express the equations of motion in terms of Brouwer's mean elements (L_j'', l_j''), which are the solution of the drag-free problem.

## The Transformation Chain

### Stage 1: Cartesian → Delaunay

The drag force X_j = -AV exp(-αr) ξ_j in Cartesian velocities transforms to Delaunay via:

P_j = Σ_k X_k ∂x_k/∂l_j = -AV exp(-αr) · p_j

where p_j = Σ_k ξ_k ∂η_k/∂l_j are the geometric projections (BH61 Eq 4-5):

p₁ = L₁(2a/r - 1),  q₁ = 2e sinE + (2/e)(L₂/L₁) sinf
p₂ = L₂,             q₂ = -(2/e) sinf
p₃ = L₃,             q₃ = 0

These are GIVEN. They follow from Keplerian orbital mechanics.

### Stage 2: Delaunay → Mean elements

The Brouwer canonical transformation (L_j, l_j) → (L_j'', l_j'') eliminates short-period and long-period terms from the gravitational Hamiltonian. The generating function is S₁ + S₁*.

Under this canonical transformation, the non-conservative force transforms as (BH61 Eq 7):

P_j'' = Σ_k P_k ∂l_k/∂l_j'' + Σ_k Q_k ∂L_k/∂l_j''

This is the COVARIANT transformation of a force 1-form — it mixes all components P_k and Q_k through the Jacobian of the coordinate change.

### Stage 3: Factor out AV exp(-αr)

Since P_k = -AV exp(-αr) p_k and Q_k = -AV exp(-αr) q_k, and AV exp(-αr) is the SAME scalar for all k at a given instant:

P_j'' = -AV exp(-αr) · [Σ_k p_k ∂l_k/∂l_j'' + Σ_k q_k ∂L_k/∂l_j'']
      = -AV exp(-αr) · p_j''

where p_j'' = p_j + δp_j (BH61 Eq 8).

### Stage 4: BH61's proof that δp_j = D(∂S₁/∂l_j)

Starting from:
δp_j = Σ_k p_k ∂(l_k - l_k'')/∂l_j'' + Σ_k q_k ∂(L_k - L_k'')/∂l_j''

Using the generating function (Eq 9):
l_k - l_k'' = -∂(S₁+S₁*)/∂L_k''
L_k - L_k'' = +∂(S₁+S₁*)/∂l_k''

Substituting, then using the chain rule to convert ∂/∂L_k and ∂/∂l_k to ∂/∂ξ_i and ∂/∂η_i (Cartesian), and applying the fundamental Poisson bracket identity {η_i, ξ_p} = -δ_ip:

δp_j = -Σ_i ξ_i · ∂[∂(S₁+S₁*)/∂l_j'']/∂ξ_i = D(∂(S₁+S₁*)/∂l_j'')

where D = -Σ_i ξ_i ∂/∂ξ_i is the Euler velocity-homogeneity operator.

## YOUR TASK

### Step 1: Verify the D operator action on Delaunay variables

From D = -Σ ξ_i ∂/∂ξ_i acting through the canonical transform (ξ,η) → (L,l), derive:

DL_j = -p_j,  Dl_j = q_j

Then derive the chain-rule actions on orbital elements:
Da, De, Dη, Dθ, Dr, Df, Dg

For EACH, show your work. Do NOT assume Df=0 — compute it explicitly from the chain rule through Kepler's equation.

### Step 2: Compute ∂S₁/∂l

This is the homological equation result (already verified in Phase A):

∂S₁/∂l = (1/n)(H₁ - ⟨H₁⟩) = Γ B₀ - μ²[B₀ + B₁'cos(2f+2g)]/r³

where Γ = μ²/(a³η³), B₀ = -1/2 + 3θ²/2, B₁' = 3(1-θ²)/2.

### Step 3: Apply D to ∂S₁/∂l

Express ∂S₁/∂l as a function of (a, η, θ, r, f, g).

Apply D using the product rule and the chain-rule identities from Step 1.

**Show every term.** The key question is which terms survive and which vanish.

### Step 4: Also compute ∂S₁*/∂l

S₁* depends on (L, G, H, g) but NOT on l. Therefore ∂S₁*/∂l = 0, and D(∂S₁*/∂l) = 0.

State this explicitly.

### Step 5: Write a numerical verification script

Compute D(∂S₁/∂l) two ways:
(a) Your closed-form result from Step 3
(b) Numerical finite differences using the Delaunay D-actions: D(F) = Σ [-p_j ∂F/∂L_j + q_j ∂F/∂l_j]

Run for: e ∈ {0.01, 0.1, 0.3}, I ∈ {30°, 60°, 85°}, g ∈ {0°, 45°, 90°}, l ∈ {0.5, 1.5, 3.0}

### Step 6: Compare against BH61 Eq(14)

BH61 Eq(14) claims (with k₂ = 1):

δp₁ = 3(μ²/L³){(-1/2+3θ²/2)[-η⁻³ + (a³/r³)(1-2a/r)] + (3/2-3θ²/2)(a³/r³)(1-2a/r)cos(2g+2f)}

Compare your derived result against this formula numerically.

If they AGREE: verify the coefficients.
If they DISAGREE: identify which terms differ and trace the discrepancy to a specific step.

### Step 7: Document

Write all results to `design/derivations/cleanroom_phase_b_results.md`.

## Constraints

- You may read: this specification, `cleanroom_phase_a_results.md` (verified S₁), `lara_2021_study_notes.md` (for cross-checking AFTER derivation only)
- You may NOT read any file containing "Brouwer_Hori" or "VERIFICATION" until Step 6
- You may NOT search the web
- Show ALL intermediate steps — do not skip any algebra

## Key Identities You Will Need

- r = a(1 - e cosE) = aη²/(1 + e cosf)
- l = E - e sinE (Kepler's equation)
- ∂f/∂l|_{a,e} = a²η/r²
- sinE = ηa sinf/r... derive from r sinf = aη sinE
- Product rule: D(FG) = (DF)G + F(DG)
- D acts as a derivation (linear, Leibniz rule)
