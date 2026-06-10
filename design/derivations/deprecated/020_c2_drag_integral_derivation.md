# Derivation 020: The C₂ Drag Coefficient

## Purpose

Derive the C₂ coefficient that governs the orbit-averaged effect of atmospheric
drag on the semi-major axis (via the Delaunay momentum $L$). Show how each
numerical coefficient in the SGP4 C₂ formula arises from the orbit-averaged
drag integral, and identify the terms dropped during the AFGP4→SGP4 simplification.

**Sources:**
- Lane, M.H. (1965), "The Development of an Artificial Satellite Theory Using a
  Power-Law Atmospheric Density Model" — original orbit-averaged drag integrals
- Lane, M.H. and Hoots, F.R. (1979), "General Perturbations Theories Derived
  from the 1965 Lane Drag Theory" — Project Space Track Report No. 2.
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Lane_Hoots_1979_General_Perturbations_Lane_Drag.pdf`
  (cited as [LH79])
- Hoots, F.R. and Roehrich, R.L. (1980), Spacetrack Report No. 3.
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Spacetrack_Report_No3_Hoots_Roehrich_1980.pdf`
  (cited as [SR3])

---

## TRANSLATE: Notation for Orbit-Averaged Drag Theory

### Symbol Table

| Symbol | Type | Definition | Units | Code Identifier |
|--------|------|-----------|-------|-----------------|
| $B^*$ | scalar | $C_D A/(2m)$ — TLE ballistic coefficient | m⁻¹ (Earth radii⁻¹) | `bstar` |
| $B$ | scalar | Brouwer notation for $B^*$ [LH79 uses $B$] | ER⁻¹ | `bstar` |
| $\rho(r)$ | function | $\rho_0((q_0-s)/(r-s))^\tau$ — power-law density | kg/m³ | — |
| $\tau$ | integer | Density exponent = 4 in SGP4. [LH79] keeps general $\tau$ | — | — |
| $q_0$ | scalar | 120 km + $a_E$ (base altitude for density model) | ER | `q0` |
| $s$ | scalar | Atmospheric fitting parameter (78-120 km + $a_E$, perigee-dependent) | ER | `s` |
| $\xi$ | scalar | $1/(a_0-s)$ — inverse atmospheric scale length | ER⁻¹ | `xi` |
| $\eta$ | scalar | $a_0 e_0 \xi$ — eccentricity-scale product | — | `eta` |
| $\psi$ | scalar | $\sqrt{1-\eta^2}$ — analogous to $\beta$ but for the density integral | — | `psisq = 1-eta^2` |
| $\gamma_1, \gamma_2$ | functions | $\gamma_1 = 1+e\cos\lambda$, $\gamma_2 = a/r = 1/(1-e\cos E)$ — orbit geometry | — | — |
| $\lambda$ | angle | Mean longitude $= l + g + h$ | rad | — |
| $I_n$ | integrals | Orbit-averaged integrals $\int \gamma_1^p \gamma_2^m d\gamma$ [LH79 notation] | — | — |
| $C_1$ | scalar | $B^* \times C_2$ — the drag rate coefficient | min⁻¹ | `C1` |
| $C_2$ | scalar | Orbit-averaged drag effect on $a$ | ER | `C2` |
| $C_3$ | scalar | $J_3$-eccentricity drag coupling | — | `C3` |
| $C_4$ | scalar | Eccentricity decay rate | — | `C4` |
| $C_5$ | scalar | Mean anomaly drag correction | — | `C5` |
| $D_2, D_3, D_4$ | scalars | Higher-order drag time polynomials | — | `D2`, `D3`, `D4` |
| $A$ | scalar | $(-1/2+3\theta^2/2)$ — secular inclination factor (same as Brouwer $A$, Eq. 012 in [B59]) | — | — |
| $n_0$ | scalar | Recovered mean motion | rad/min | `n0` |
| $a_0$ | scalar | Recovered semi-major axis | ER | `a0` |

**Key notation bridge:** [LH79] uses $B$ where SGP4 uses $B^*$. [LH79] uses $\psi = |1-\eta^2|^{1/2}$. [LH79] keeps general power-law exponent $\tau$; SGP4 fixes $\tau = 4$.

---

## BUILD Step 1: Physical Setup (020.Eq.1-5)

### 1.1: Drag acceleration

The drag force on a satellite of mass $m$, cross-section $A$, drag coefficient $C_D$:

(020.Eq.1) $\mathbf{F}_{drag} = -\frac{1}{2}C_D A \rho |\mathbf{v}_{rel}|^2 \hat{\mathbf{v}}_{rel}$

For a non-rotating atmosphere (the SGP4 assumption), $\mathbf{v}_{rel} = \mathbf{v}$.

### 1.2: Ballistic coefficient

(020.Eq.2) $B^* \equiv \frac{C_D A}{2m}$

This is the quantity stored in the TLE. It encapsulates the satellite's susceptibility to drag.

### 1.3: Power-law atmospheric density model

(020.Eq.3) $\rho(r) = \rho_0\left(\frac{q_0 - s}{r - s}\right)^\tau$

where $\tau = 4$ for SGP4. The parameters $s$ and $q_0$ are fitting constants chosen so that the model approximates the actual atmosphere over the altitude range of interest.

### 1.4: Gauss variational equation for $L$

The drag effect on the Delaunay momentum $L = \sqrt{\mu a}$ is ([LH79], Section 2):

(020.Eq.4) $\frac{dL}{dt}\bigg|_{drag} = -B^* \cdot n \cdot a \cdot \rho(r) \cdot v^2 \cdot \frac{r^2}{2\mu}$

This is derived from the Gauss VE $da/dt = (2a^2/\mu)\mathbf{v}\cdot\mathbf{F}/m$ with $\mathbf{F}/m = -(1/2)C_D(A/m)\rho v^2\hat{v}$ and $L = \sqrt{\mu a}$, so $dL/dt = (\mu/(2L))da/dt$.

### 1.5: Orbit geometry in terms of true anomaly

(020.Eq.5a) $r(f) = \frac{a(1-e^2)}{1+e\cos f} = \frac{p}{1+e\cos f}$

(020.Eq.5b) $v^2(f) = \frac{\mu}{p}(1+2e\cos f + e^2)$ (from vis-viva)

---

## BUILD Step 2: The AFGP4 Orbit-Averaged Integral (020.Eq.6-15)

### 2.1: The orbit average

The secular change in $L$ over one orbit is obtained by averaging (020.Eq.4) over one orbital period. [LH79] performs this in terms of the eccentric anomaly $\lambda$ (their notation for mean longitude) using the integral:

> **[UNVERIFIED — transcribed from scanned PDF: Lane_Hoots_1979, p. 16. No clean post-LaTeX source available for Lane drag theory.]**

(020.Eq.6) $\langle\Delta L\rangle = L_0\left\{-B\alpha a \int_{\lambda_o}^{\lambda_s} (\gamma_2^{\tau-1} + 2e\gamma_1\gamma_2^{\tau-2} + \frac{3}{2}e^2\gamma_1^2\gamma_2^{\tau-3} + \ldots)d\lambda + \text{J}_2\text{ terms}\right\}$

where $\gamma_1 = 1+e\cos\lambda$, $\gamma_2 = a/r$, and $\alpha$ is a constant involving the density model parameters.

[LH79] evaluates these integrals for general $\tau$ using the integral formula on p. 16 (the combinatorial expansion of $\int\gamma_1^p\gamma_2^m d\gamma$). For $\tau = 4$, the integrals reduce to specific combinations.

### 2.2: The AFGP4 result (full form)

After evaluating all integrals, [LH79] p. 25 gives the secular part of $L''$:

> **[UNVERIFIED — transcribed from scanned PDF: Lane_Hoots_1979, p. 25. No clean post-LaTeX source available. Code output validated against Vallado test cases.]**

(020.Eq.7) $L'' = L_0''\left\{1 - B\psi^{-7}\xi^4 a\left(1+\frac{3}{2}\eta^2+4e\eta+e\eta^3+\frac{3}{4}e^2+3e^2\eta^2\right)(\lambda_s''-\lambda_o'')\ \right.$
$\left. - \frac{3}{2}Bk_2\xi^4(-\frac{1}{2}+\frac{3}{2}\theta^2)[\xi\psi^{-9}(8+24\eta^2+3\eta^4) - 5e\eta(4+3\eta^2)](\ell_s''-\ell_o'')\right\} + L_o'' Q_2$

This is the **complete AFGP4 drag integral** before any simplification. The two parts are:

**Part A (Keplerian drag):** $a(1+\frac{3}{2}\eta^2+4e\eta+e\eta^3+\frac{3}{4}e^2+3e^2\eta^2)$

The coefficients arise from the orbit-averaged integrals of $(r/a)^n \cdot v^m \cdot \rho(r)$ over one revolution:
- $1$: the leading term (circular orbit, no eccentricity correction)
- $\frac{3}{2}\eta^2$: from the $e^2$ correction to the orbit-averaged density
- $4e\eta$: from the first-order eccentricity-velocity coupling
- $e\eta^3$: from the third-order eccentricity-density coupling
- $\frac{3}{4}e^2$: from the vis-viva velocity correction at $O(e^2)$ [dropped in SGP4]
- $3e^2\eta^2$: from the cross-term of eccentricity and density at $O(e^2\eta^2)$ [dropped in SGP4]

**Part B (J₂ correction):** $\frac{3}{2}k_2(-\frac{1}{2}+\frac{3}{2}\theta^2)\xi\psi^{-1}[(8+24\eta^2+3\eta^4) - 5e\eta(4+3\eta^2)]$

This comes from the J₂ radial perturbation's effect on the atmospheric density:
- $(8+24\eta^2+3\eta^4)$: from orbit-averaging the density perturbation $\delta\rho = (\partial\rho/\partial r)\delta r_{J_2}$
- $-5e\eta(4+3\eta^2)$: from the eccentricity-dependent part of the J₂ radial correction [dropped in SGP4]

### 2.3: The intermediate quantities

(020.Eq.8) $\xi = \frac{1}{a_0 - s}$

(020.Eq.9) $\eta = a_0 e_0 \xi$

(020.Eq.10) $\psi^2 = |1-\eta^2|$

The COEF1 prefactor combines density model parameters:

(020.Eq.11) $\text{COEF1} = (q_0-s)^4 \xi^4 \psi^{-7}$

This arises from the density function evaluation:
$\rho(r) = \rho_0((q_0-s)/(r-s))^4$. At perigee, $r_{perigee} = a(1-e)$, so $r_{perigee}-s = a-s-ae = 1/\xi - a_0 e_0 = (1-\eta)/\xi$. The $(q_0-s)^4\xi^4$ factor comes from $\rho_0 \times (q_0-s)^4 \times \xi^4$, and $\psi^{-7} = (1-\eta^2)^{-7/2}$ comes from orbit-averaging the density with its $\tau = 4$ exponent: the $(r-s)^{-4}$ factor, when expanded around the mean value and averaged, produces the $(1-\eta^2)^{-7/2}$ factor.

---

## BUILD Step 3: The SGP4 Simplification (020.Eq.16-25)

### 3.1: Terms dropped from Part A

[LH79] Section 3a, p. 26: "dropping certain smaller terms" produces:

(020.Eq.12) Part A (SGP4): $a_0(1+\frac{3}{2}\eta^2+4e_0\eta+e_0\eta^3)$

The dropped terms $\frac{3}{4}e^2$ and $3e^2\eta^2$ are $O(e^2)$ relative to the leading term. For a typical LEO with $e = 0.01$, these are $\sim 10^{-4}$ corrections — negligible compared to the $\eta^2$ and $e\eta$ terms.

### 3.2: Terms dropped from Part B

(020.Eq.13) Part B (SGP4): $\frac{3}{2}k_2\xi\psi^{-2}(-\frac{1}{2}+\frac{3}{2}\theta^2)(8+24\eta^2+3\eta^4)$

The dropped term $-5e\eta(4+3\eta^2)$ is $O(e)$ relative to the $(8+24\eta^2+3\eta^4)$ term. This is a larger simplification than Part A's $O(e^2)$ drops.

### 3.3: Assembly of SGP4 C₂

The SGP4 C₂ is defined so that $C_1 = B^* C_2$ and $L'' = L_0''[1-C_1(t-t_0)+...]$.

From [LH79] p. 26:

> **[UNVERIFIED — transcribed from scanned PDF: Lane_Hoots_1979, p. 26. No clean post-LaTeX source available. Code output validated against Vallado test cases.]**

(020.Eq.14) $C_1 = B\xi^4 n\psi^{-7}\left[a\left(1+\frac{3}{2}\eta^2+4e\eta+e\eta^3\right) + \frac{3}{2}k_2\xi\psi^{-2}\left(-\frac{1}{2}+\frac{3}{2}\theta^2\right)(8+24\eta^2+3\eta^4)\right]$

Therefore:

> **[PARTIALLY VERIFIED — transcribed from scanned PDF: Lane_Hoots_1979, p. 26 / SR3 p. 11. The 3/8 J₂ coefficient is verified by algebra (020.Eq.16-17). Code output validated against Vallado test cases.]**

(020.Eq.15) $C_2 = (q_0-s)^4\xi^4 n_0\psi^{-7}\left[a_0\left(1+\frac{3}{2}\eta^2+4e_0\eta+e_0\eta^3\right) + \frac{3}{2}\frac{k_2\xi}{\psi^2}\left(-\frac{1}{2}+\frac{3}{2}\theta^2\right)(8+24\eta^2+3\eta^4)\right]$

This matches [SR3] p. 11 exactly. ✓

### 3.4: Correspondence to code (drag_coefficients.h lines 141-144)

The [SR3] p. 11 formula for Part B is:

$$+\frac{3}{2}\frac{k_2\xi}{(1-\eta^2)}\left(-\frac{1}{2}+\frac{3}{2}\theta^2\right)(8+24\eta^2+3\eta^4)$$

Substituting $k_2 = J_2/2$ and $(-1/2+3\theta^2/2) = (3\cos^2 i-1)/2$:

$$= \frac{3}{2}\cdot\frac{J_2}{2}\cdot\frac{\xi}{\psi^2}\cdot\frac{(3\cos^2 i-1)}{2}\cdot(8+24\eta^2+3\eta^4)$$

$$= \frac{3J_2}{8}\cdot\frac{\xi}{\psi^2}\cdot(3\cos^2 i-1)\cdot(8+24\eta^2+3\eta^4) \tag{020.Eq.16}$$

The code computes this as:

```
ratio<T>(3, 4) * in.half_J2 * dc.xi / psisq * in.three_cos2i_minus_1 * (...)
```

where `ratio<T>(3,4)` = 3/4, `in.half_J2` = J₂/2, and `in.three_cos2i_minus_1` = 3cos²i - 1.

$(3/4) \times (J_2/2) \times \xi/\psi^2 \times (3\cos^2 i-1) = (3J_2/8) \times \xi/\psi^2 \times (3\cos^2 i-1)$ ✓

(020.Eq.17) **The existing code comment "(3/8)J₂" is CORRECT.** The coefficient of $J_2$ in Part B, when factored with $(3\cos^2 i - 1)$ as the code does, is indeed $3/8$. This coefficient traces to $\frac{3}{2}k_2 \times A = \frac{3}{2}\cdot\frac{J_2}{2}\cdot\frac{(3\cos^2 i-1)}{2} = \frac{3J_2}{8}(3\cos^2 i-1)$. The factor $3/2$ in front of $k_2$ comes from the orbit-averaged J₂-density coupling integral [LH79, p. 25].

---

## BUILD Step 4: Higher-Order Drag Coefficients (020.Eq.17-25)

From [LH79] Section 3, the remaining drag coefficients:

### C₃ (J₃ eccentricity-drag coupling) [LH79, p. 29]

> **[UNVERIFIED — transcribed from scanned PDF: Lane_Hoots_1979, p. 29. No clean post-LaTeX source available.]**

(020.Eq.17) $C_3 = \frac{(q_0-s)^4\xi^5 A_{3,0} n_0 \sin i_0}{k_2 e_0}$

where $A_{3,0} = -J_3 a_E^3/2$. This comes from the $J_3$ term in the radial perturbation, which produces a $\sin g$ (long-period) drag effect.

### C₄ (eccentricity decay rate) [SR3, p. 11]

> **[UNVERIFIED — transcribed from scanned PDF: SR3 p. 11. No clean post-LaTeX source available. Code output validated against Vallado test cases.]**

(020.Eq.18) $C_4 = 2n_0(q_0-s)^4\xi^4 a_0\beta_0^2\psi^{-7}\left\{\left[2\eta(1+e_0\eta)+\frac{1}{2}e_0+\frac{1}{2}\eta^3\right] - \frac{2k_2\xi}{a_0\psi^2}\times\right.$
$\left.[3(1-3\theta^2)(1+\frac{3}{2}\eta^2-2e_0\eta-\frac{1}{2}e_0\eta^3)+\frac{3}{4}(1-\theta^2)(2\eta^2-e_0\eta-e_0\eta^3)\cos 2\omega_0]\right\}$

### C₅ (mean anomaly drag correction) [SR3, p. 11]

> **[UNVERIFIED — transcribed from scanned PDF: SR3 p. 11. No clean post-LaTeX source available. Code output validated against Vallado test cases.]**

(020.Eq.19) $C_5 = 2(q_0-s)^4\xi^4 a_0\beta_0^2\psi^{-7}\left[1+\frac{11}{4}\eta(\eta+e_0)+e_0\eta^3\right]$

### D₂, D₃, D₄ (time polynomial coefficients) [LH79, p. 26]

> **[UNVERIFIED — transcribed from scanned PDF: Lane_Hoots_1979, p. 26. No clean post-LaTeX source available. Code output validated against Vallado test cases.]**

(020.Eq.20) $D_2 = 4a_0\xi C_1^2$

(020.Eq.21) $D_3 = \frac{4}{3}a_0\xi^2(17a_0+s)C_1^3$

(020.Eq.22) $D_4 = \frac{2}{3}a_0^2\xi^3(221a_0+31s)C_1^4$

These arise from the Taylor expansion of the orbit-averaged drag integral to higher powers of $(t-t_0)$, accounting for the changing semi-major axis on the drag rate itself.

---

## PAUSE on C₂ Derivation

### Correctness

1. **Code coefficient verified:** The code computes `ratio<T>(3,4) * in.half_J2` $= (3/4)(J_2/2) = 3J_2/8$, multiplied by `three_cos2i_minus_1` $= (3\cos^2 i - 1) = 2A$. This gives $(3J_2/8)(3\cos^2 i - 1) = (3/2)k_2 A$, matching [SR3] p. 11. The existing code comment "(3/8)J₂" correctly describes the coefficient of $J_2$ when factored with $(3\cos^2 i - 1)$, as shown in 020.Eq.16-17.
2. **Part A polynomial verified:** $1+\frac{3}{2}\eta^2+4e\eta+e\eta^3$ matches [LH79] p. 25 and [SR3] p. 11.
3. **Part B polynomial verified:** $(8+24\eta^2+3\eta^4)$ matches [LH79] p. 25 and [SR3] p. 11.
4. **COEF1 verified:** $(q_0-s)^4\xi^4(1-\eta^2)^{-7/2}$ matches [SR3] p. 11.

### Accuracy-limiting assumptions

1. **Power-law exponent fixed at $\tau=4$:** The AFGP4 theory keeps general $\tau$, but SGP4 fixes $\tau=4$. A generalized exponent would change all the orbit-averaged integrals. For the enhanced preset: parameterize $\tau$ and recompute the integrals from [LH79]'s general formula (p. 16).

2. **Non-rotating atmosphere:** SGP4 assumes $\mathbf{v}_{rel} = \mathbf{v}$, ignoring atmospheric co-rotation. This underestimates drag for prograde LEO by ~5% and overestimates for retrograde. The correction is $v_{rel} = v - \omega_E r\cos i$, which modifies the drag integral by a factor $(1-\omega_E r\cos i/v)$.

3. **Dropped $O(e^2)$ terms in Part A:** $\frac{3}{4}e^2 + 3e^2\eta^2$ — for $e=0.1$, this is ~1% of Part A.

4. **Dropped $-5e\eta(4+3\eta^2)$ in Part B:** For $e=0.1$, $\eta\approx 0.1$, this is $-5\times 0.01\times(4+0.03) \approx -0.2$, compared to $(8+24\times 0.01+3\times 0.0001) \approx 8.24$. So ~2.4% of Part B. Significant for eccentric orbits.

5. **J₂-density cross-coupling truncated at first order:** The density perturbation $\delta\rho = (\partial\rho/\partial r)\delta r_{J_2}$ was linearized. The second-order correction $\frac{1}{2}(\partial^2\rho/\partial r^2)(\delta r)^2$ is $O(J_2^2)$ and likely negligible.

### Precision improvements

1. **The polynomial $(8+24\eta^2+3\eta^4)$:** In Horner form: $8+\eta^2(24+3\eta^2)$. The code already uses this: `exact<T>(8) + exact<T>(3) * eta2 * (exact<T>(8) + eta2)`. ✓

2. **Near-circular orbits ($e\to 0$):** When $e\to 0$, $\eta\to 0$ and Part A $\to a_0$, Part B $\to (3/2)k_2\xi\psi^{-2}(-1/2+3\theta^2/2)\times 8$. No cancellation issues.

3. **COEF1 overflow/underflow:** $(1-\eta^2)^{-7/2}$ can be very large when $\eta\to 1$ (perigee approaching the atmospheric fitting altitude $s$). The code should track this via error bounds on $\eta$.

### Generalizations for Enhanced Preset

1. **Full AFGP4 Part A:** Restore $+\frac{3}{4}e^2+3e^2\eta^2$ terms.
2. **Full AFGP4 Part B:** Restore $-5e\eta(4+3\eta^2)$ term.
3. **General power-law exponent $\tau$:** Replace all $\tau=4$ specializations with general integrals from [LH79] p. 16.
4. **Co-rotating atmosphere correction:** Modify $v^2$ in the drag integral.

These generalizations are injectable via lambdas in the existing architecture — the standard preset uses the SGP4 formulas, while the enhanced preset uses the AFGP4 or further-generalized forms.

---

## Derivation Chain Summary

The AFGP4 constant `k` on [LH79] p. 2 is the orbit-averaged drag integral:

$$k = B\mu(1-\eta^2)^{-7/2}\{(\text{Part A: Keplerian}) + (\text{Part B: J}_2\text{ correction})\}$$

**Lane (1965)** → orbit-averaged integrals with power-law density → AFGP4 `k`

**[LH79] Section 3** → drops $O(e^2)$ from Part A and $O(e)$ from Part B → SGP4 $C_1$

**[SR3] p. 11** → presents the final SGP4 form with $C_2 = C_1/B^*$

**This derivation (020)** → traces each coefficient in the code to the specific integral evaluation in [LH79], identifies the dropped terms, and proposes generalizations.

### Equation Reference for Code Verification

| Code location | Expression | Derivation equation |
|--------------|-----------|-------------------|
| `drag_coefficients.h:121` | `coef1` = $(q_0-s)^4\xi^4\psi^{-7}$ | 020.Eq.11 |
| `drag_coefficients.h:141` | Part A: $a_0(1+\frac{3}{2}\eta^2+4e_0\eta+e_0\eta^3)$ | 020.Eq.12 |
| `drag_coefficients.h:143` | Part B: $\frac{3J_2}{8}\xi\psi^{-2}(3\cos^2 i-1)(8+24\eta^2+3\eta^4)$ | 020.Eq.13, 020.Eq.16-17 |
| `drag_coefficients.h:147` | $C_1 = B^* C_2$ | 020.Eq.14 |
