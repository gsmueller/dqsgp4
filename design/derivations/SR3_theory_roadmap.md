# SR3/SGP4 Theory Roadmap

## The actual source: Lane & Hoots (1979) "General Perturbations Theories Derived from the 1965 Lane Drag Theory"

Located at: `sgp4_references/vallado_celestrak/documentation/SGP4/Lane_Hoots_1979_General_Perturbations_Lane_Drag.pdf`

## Document structure

| Section | Pages | Content |
|---------|-------|---------|
| 1. Introduction | 1 | Lineage: Lane (1965) → Lane & Cranford (1969) → AFGP4 → IGP4 → SGP4 |
| 2. AFGP4 | 1-24 | **THE FULL THEORY** — every equation derived |
| 3. Drag Simplification | 25-33 | AFGP4 → IGP4 simplification |
| 4. IGP4 | 34-38 | IGP4 equations |
| 5. Geopotential Simplification | 39-46 | IGP4 → SGP4 simplification |
| 6. SGP4 | 47-51 | SGP4 equations (same as SR3 Section 6) |
| 7. References | 52 | |
| 8. Appendix | 53+ | Symbol definitions |

## The AFGP4 theory structure (Section 2)

### Assumptions
- Power-law density: ρ = ρ₀((q₀-s)/(r-s))^τ
- Brouwer (1959) gravitational model (J₂ through J₅)
- Non-rotating, spherically symmetric atmosphere
- Method of successive approximations on Gauss variational equations

### Key quantities (page 2)
- k = Bμ(1-η²)^(-7/2) {drag coefficient with J₂ coupling}
- The k formula contains the (-1/2 + 3/2 θ²) term — this is the drag-oblateness coupling
- A₁ through A₄: auxiliary constants for the time polynomial
- Q₁, Q₂: time-dependent drag integrals

### Secular rates (page 3)
- l'' (mean anomaly rate) — Brouwer secular with J₂, J₂², J₄
- g'' (arg perigee rate) — Brouwer secular
- h'' (RAAN rate) — Brouwer secular

### The core integrals (pages 4-7)
- I₁ = ∫∫ (density integral over true anomaly and time) — THE MAIN DRAG INTEGRAL
- I₂ = ∫∫ (similar structure, different powers)
- I₃, I₄ = single integrals over true anomaly
- These use γ₁, γ₂ = functions of e, η, and the density model
- J_{τ-k} = Bessel-like integrals evaluated between λ₀'' and λ_s''

### The variational equations (pages 5-8)
- ∫ΔL_i dt = main drag integral for semimajor axis decay
- ∫ΔG_i dt = angular momentum change
- ∫ΔH_i dt = nodal angular momentum change
- These contain the J₂ coupling terms (the 3k₂/2a² factors)
- F₂ expression (page 10) — the secular mean longitude update

### The iterative scheme (pages 8-9)
- l''_{i+1}, g''_{i+1}, h''_{i+1} computed from ΔL_i, ΔG_i, ΔH_i
- Number of iterations depends on |n∫ΔL_i dt|
- This is the Method of Successive Approximations

### From AFGP4 to SGP4 (Sections 3-6)
- Section 3: Drag simplification — drops higher-order eccentricity terms
- Section 4: IGP4 — intermediate simplified theory
- Section 5: Geopotential simplification — drops J₂² cos(2g) terms
- Section 6: SGP4 — final simplified equations (= SR3 Section 6)

## What needs to be done

1. **Transcribe Section 2 (AFGP4) into markdown** — this is the theory
2. **For each equation, identify the physical origin** — which term is Brouwer, which is drag, which is coupling
3. **Trace the simplification chain** — what AFGP4 terms survive in SGP4, what gets dropped
4. **Write Ch 21-22 of the textbook** from this understanding

## Key insight from page 2

The drag coefficient k contains:
```
k = Bμ(1-η²)^(-7/2) { (1 + 3/4 e² + 3/2 η² + 3e²η² + 4eη + eη³)
    + (3k₂/2a)(-1/2 + 3/2 θ²) ξ(1-η²)^(-1) [(8+24η²+3η⁴) - 5eη(4+3η²)] }
```

The second line is the J₂-drag coupling. It enters through the DENSITY VARIATION — when the orbit precesses under J₂, perigee samples different atmospheric densities. This is NOT the same as BH61's δp₁ (which is the force-direction correction). Lane's theory handles the coupling through the density integral, not through the D operator.

This is why BH61 was the wrong map. The coupling mechanism in SGP4 is fundamentally different.
