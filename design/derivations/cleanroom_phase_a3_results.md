# Cleanroom Phase A3: D Operator vs Poisson Bracket Analysis

**Date:** 2026-04-05
**Script:** `design/derivations/cleanroom_phase_a3_analysis.py`
**Verification:** All claims backed by 81-point numerical grids (e ∈ {0.01, 0.1, 0.3}, I ∈ {30°, 60°, 85°}, g ∈ {0°, 45°, 90°}, l ∈ {0.5, 1.5, 3.0})

---

## The Question

BH61 defines the oblateness correction to the drag force as:

$$\delta p_j = D\left(\frac{\partial(S_1 + S_1^*)}{\partial l_j''}\right)$$

where D is the Euler velocity-homogeneity operator. An earlier phase (A2) computed the Poisson bracket {p₁, S₁} instead. This analysis determines:

1. Whether these are the same quantity
2. What each equals
3. Which matches BH61 Eq(14)

---

## Phase 1: D-Action Identities (Verified)

The D operator is defined on Delaunay variables via BH61 Eq(5) and Eq(10):

| Variable | D action | Verification |
|----------|----------|-------------|
| DL = -p₁ = -L(2a/r-1) | Definition | — |
| DG = -p₂ = -G | Definition | — |
| DH = -p₃ = -H | Definition | — |
| Dl = q₁ = 2esinE + 2ηsinf/e | Definition | — |
| Dg = q₂ = -2sinf/e | Definition | — |
| Dh = q₃ = 0 | Definition | — |
| Da = -2a(2a/r-1) | Chain rule from DL | ✓ matches to ~1e-9 |
| De = -2(e+cosf) | Chain rule | ✓ matches to ~1e-9 |
| Dη = 2η(a/r-1) | Chain rule | ✓ |
| Dθ = 0 | D(H/G) = (G·DH-H·DG)/G² = 0 | ✓ ~1e-10 |
| Dr = 0 | r is pure position function | ✓ ~1e-7 |
| **Df = +2sinf/e** | Chain rule through Kepler | ✓ matches to ~1e-9 |
| **D(f+g) = 0** | Df = -Dg | ✓ ~1e-7 |

**Critical finding:** Df ≠ 0. Phase A §3.1 incorrectly claimed Df = 0. The correct value is Df = 2sinf/e = -Dg, giving D(f+g) = 0.

---

## Phase 2: Three Computations of δp₁ (Compared)

| Computation | Definition | Result |
|-------------|-----------|--------|
| **(A) -{p₁, S₁}** | Negative Poisson bracket of p₁ with S₁ | Rich harmonic structure: cos(2g) amplitude ~4.6 at e=0.1 |
| **(B) D(∂S₁/∂l), Delaunay FD** | D operator via DL=-p₁, Dl=q₁, etc. | **6ΓB₀(a/r)**, no g-dependence |
| **(C) D(∂S₁/∂l), Cartesian FD** | D = -Σ vᵢ ∂/∂vᵢ applied directly | Agrees with (B) at g=45°, 90° to ~1e-9; deviates at g=0° by ~0.2-3.5% due to coordinate conversion numerics |
| **(D) BH61 Eq(14)** | 3(μ²/L³){B₀[-η⁻³+ρ³(1-2ρ)] + B₁ρ³(1-2ρ)cos(2f+2g)} | Rich harmonic structure, different from all of A, B, C |

### Key ratios from the data:

**A/B (Poisson bracket vs D operator):** Varies wildly from -14.8 to +5.5. These are **NOT the same quantity.**

**B/C (Delaunay D vs Cartesian D):** 1.0000 at g=45° and g=90°. Deviates at g=0° due to coordinate conversion singularity (arccos near 0). **The D operator definition is correct.**

**B/D (D operator vs BH61 Eq 14):** Varies wildly. **BH61 Eq(14) does not match D(∂S₁/∂l).**

**A/D (Poisson bracket vs BH61):** Also varies. **BH61 Eq(14) does not match -{p₁, S₁} either.**

### Conclusion: All four quantities are different.

---

## Phase 3: Closed-Form for D(∂S₁/∂l) (Verified)

**Derivation:**

$$\frac{\partial S_1}{\partial l} = \Gamma B_0 - \frac{\mu^2}{r^3}\left[B_0 + B_1'\cos(2f+2g)\right]$$

Applying D:
- D[ΓB₀] = B₀·D[Γ] = B₀·6Γ(a/r) ← from D[a⁻³η⁻³] = 6(a/r)/(a³η³)
- D[μ²B₀/r³] = 0 ← because Dr = 0, so D[r⁻³] = 0, and DB₀ = 0
- D[μ²B₁'cos(2f+2g)/r³] = 0 ← because D[r⁻³] = 0 AND D[cos(2f+2g)] = -sin(2f+2g)·2·D(f+g) = 0

$$\boxed{D\left(\frac{\partial S_1}{\partial l}\right) = 6\,\Gamma\,B_0\,\frac{a}{r} = \frac{6\mu^2 B_0}{a^2 r \eta^3}}$$

**Verification:** 81/81 PASS, max relative error 5.09e-6 (finite-difference limited).

**Properties:**
- Proportional to B₀ = -1/2 + 3cos²I/2 only
- No B₁' = 3sin²I/2 term
- No cos(2f+2g) harmonic
- No constant offset
- Purely proportional to a/r

---

## Phase 4: Harmonic Decomposition

| Quantity | Mean | |cos(2g)| | |sin(2g)| | Has long-period terms? |
|----------|------|-----------|-----------|----------------------|
| -{p₁, S₁} (e=0.1) | +0.718 | 4.581 | 1.555 | **YES** |
| D(∂S₁/∂l) (e=0.1) | -0.759 | 0.000 | 0.000 | **NO** |
| BH61 Eq(14) (e=0.1) | +0.750 | 3.216 | 0.850 | **YES** |

D(∂S₁/∂l) has **zero** g-dependence. The Poisson bracket and BH61 both have strong cos(2g) components, but with different amplitudes. None of the three quantities agree with each other.

---

## External References (CORRECTED 2026-04-05)

### INPE-2746 — misattributed corroboration

INPE-2746-PRE/322 (Fitzgibbon, De Moraes, Lobão, 1983) states:

> "Esta teoria contém termos de Poisson espúrios [5]. Para utilizá-la devem-se negligenciar todos estes termos [6]"

("This theory contains spurious Poisson terms. To use it, all these terms must be neglected.")

**CORRECTION:** This was originally cited as corroboration that BH61 Eq(14) has spurious terms. Reading the paper reveals INPE-2746 never mentions Eq(14), the D operator, or cos(2f+2g). The paper does not explain what the spurious terms are or where they appear — the underlying analysis is in Fitzgibbon's thesis [5] and Vilhena de Moraes (1981) [6], which we do not have.

De Moraes (1981), *Celestial Mechanics* **25**, 281-292, uses Lagrange's Variation of Parameters to avoid the integration artifacts. This is about the integration method, not about Eq(14).

---

## Open Questions (UPDATED 2026-04-05)

1. **The Poisson bracket {p₁, S₁} and the D operator D(∂S₁/∂l) are different quantities.** Confirmed numerically (0/81 match). The physical interpretation of each needs clarification: which quantity enters the actual equations of motion for the mean elements under drag+oblateness?

2. ~~**BH61's Eq(14) matches neither.**~~ **RESOLVED:** BH61 Eq(14) matches D(Brouwer's ∂S₁/∂l) — a different function from D(our ∂S₁/∂l). The two forms of ∂S₁/∂l differ by a velocity-dependent factor -(μ/a)^(3/2) that equals -1 only when a=μ (the normalized test case). Both D computations are correct for their respective functions. See `verify_normalization_ratio.py`, `verify_D_brouwer_form.py`.

3. **The Lagrange approach (INPE/De Moraes)** computes the total rate dL₁''/dt directly via Variation of Parameters. INPE-2746 claims this avoids spurious Poisson terms present in BH61, but does not detail what those terms are. The full analysis is in Fitzgibbon (1982) and Vilhena de Moraes (1981), which we do not have.

---

## Files

- `design/derivations/cleanroom_phase_a3_analysis.py` — Single script reproducing all results
- `design/derivations/cleanroom_phase_a_results.md` — S₁ derivation (72/72 verified)
- `design/derivations/cleanroom_phase_a2_results.md` — Poisson bracket framework (81/81 verified)
