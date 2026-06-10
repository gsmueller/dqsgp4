# Draft Plan: Chapter 32 — Element Recovery

**Part IX: The SGP4 Propagator** | Implementation file: `element_recovery.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $n_o$ | TLE mean motion (rad/min after conversion) | Ch 31 |
| $n_o''$ | Recovered mean motion after $\delta$ corrections | §32.2 |
| $a_o''$ | Recovered semi-major axis after corrections | §32.2 |
| $\delta_1$ | First correction: $\delta_1 = \frac{3}{2}k_2 \frac{(3\cos^2 i_o - 1)}{a_1^2(1-e_o^2)^{3/2}} \cdot \frac{1}{3}$ (form stated in §32.2) | §32.2 |
| $a_1$ | First estimate: $a_1 = (k_e / n_o)^{2/3}$ | §32.2 |
| $a_o$ | Intermediate estimate: $a_o = a_1(1 - \delta_1/3 - \delta_1^2 - 134\delta_1^3/81)$ | §32.2 |
| $\delta_o$ | Final correction: $\delta_o = \frac{3}{2}k_2\frac{3\cos^2 i_o - 1}{a_o^2(1-e_o^2)^{3/2}}$ | §32.2 |
| $B^*$ | SGP4 drag coefficient from TLE | Ch 31 |
| $C_D A/(2m)$ | Ballistic coefficient | §32.4 |
| $k_e$ | Gaussian gravitational constant (Ch App A) | Appendix A |
| $k_2$ | SGP4 second zonal constant: $k_2 = \tfrac{1}{2}J_2 a_E^2$ | Appendix A |
| $s^*$ | Adjustable drag parameter; perigee-dependent | §32.5 |
| $q_0 - s^*$ | Effective parameter for density model | Ch 21 |
| $r_{perigee}$ | Perigee radius: $a_o''(1 - e_o) - a_E$ in km altitude | §32.5 |

---

## Objectives

1. Derive the iterative sequence $a_1 \to \delta_1 \to a_o \to \delta_o \to a_o'', n_o''$ from the Brouwer mean-to-osculating relationship.
2. Prove convergence of the iteration and bound the residual precision error.
3. Derive the third-order series $a_o = a_1(1 - \delta_1/3 - \delta_1^2 - 134\delta_1^3/81)$ algebraically (not by matching to SR3).
4. Define $B^*$ and its relation to the ballistic coefficient $C_D A/(2m)$.
5. Derive the three perigee altitude thresholds (98, 156, 220 km) and their effect on $s^*$ and $(q_0-s^*)$.
6. Specify the cube root Newton iteration and prove its convergence.
7. Generalization: higher-order series inversion.

## Section Structure

### §32.1 Introduction

This section explains why TLE mean motion is not the Keplerian mean motion and provides a road map of the recovery iteration.

Stub: The TLE mean motion $n_o$ is not the pure Keplerian mean motion; it incorporates the SGP4 mean-element conventions. Element recovery is the process of inverting this relationship to obtain $n_o''$ and $a_o''$ — the fully corrected mean motion and semi-major axis — from which all subsequent quantities are computed. Every coefficient in the iteration must be derived, not matched to code. [A.32.1] The recovery procedure is model-specific: it is valid only within the SGP4 matched pair (Ch 3). Using a different $k_2$ or $J_2$ breaks the pair and introduces error.

### §32.2 The Recovery Iteration

This section derives the five-step sequence $a_1 \to \delta_1 \to a_o \to \delta_o \to (a_o'', n_o'')$ from the first-order Brouwer secular correction and proves convergence.

Stub: Step 1: $a_1 = (k_e/n_o)^{2/3}$ (Kepler's third law; $k_e$ from Appendix A). Step 2: $\delta_1 = \tfrac{3}{2}k_2(3\cos^2 i_o - 1) / (a_1^2(1-e_o^2)^{3/2}) / 3$. Derive $\delta_1$ as the first Brouwer secular correction to $n$ expressed in terms of $a_1$ (cross-reference Ch 16, Theorem 16.4.1). Step 3: expand $a_o = a_1(1 - \delta_1/3 - \delta_1^2 - 134\delta_1^3/81)$ — derive this as the third-order Maclaurin expansion of $(1 + \delta_1)^{-2/3}$ to third order [derive each coefficient algebraically]. Step 4: $\delta_o = $ analogous expression with $a_o$ in place of $a_1$. Step 5: $n_o'' = n_o/(1 + \delta_o)$, $a_o'' = (k_e/n_o'')^{2/3}$. **Theorem 32.2.1** (convergence): $|\delta_1| < 0.01$ for all physical orbits; the third-order series is accurate to $|\delta_1^4|/4 < 10^{-8}$ in $a_o$. — *Proof approach: bound $|\delta_1|$ by substituting the extremal physical values $i = 0$, $e = 0$ and the minimum physical $a_1 > a_E$ into the $\delta_1$ formula; evaluate the fourth-order Taylor remainder using $|\delta_1| \leq 0.01$.* [P.32.1] The cube root $a_1 = (k_e/n_o)^{2/3}$ is evaluated by Newton iteration (§32.6); introduces precision error.

### §32.3 Derivation of the 134/81 Coefficient

This section derives the coefficient 134/81 by explicit term-by-term evaluation of the Maclaurin expansion of $(1+x)^{-2/3}$ at $x = -\delta_1$.

Stub: Expand $(1 + x)^{-2/3}$ as a Maclaurin series: $(1+x)^{-2/3} = 1 - \tfrac{2}{3}x + \tfrac{2}{3}\cdot\tfrac{5}{3}\cdot\tfrac{x^2}{2!} - \tfrac{2}{3}\cdot\tfrac{5}{3}\cdot\tfrac{8}{3}\cdot\tfrac{x^3}{3!} + \cdots$. Substitute $x = -\delta_1$ and collect: coefficient of $\delta_1^3$ is $\tfrac{2\cdot5\cdot8}{3^3\cdot6} = \tfrac{80}{162} = \tfrac{40}{81}$... [work through the algebra explicitly]. Final coefficient $134/81$ — derive it step by step. This section exists because the value 134/81 must be algebraically verified, not read from SR3. [A.32.2] Third-order truncation; the fourth-order term is $O(\delta_1^4) < 10^{-8}$ for physical orbits.

### §32.4 B* and the Ballistic Coefficient

This section defines $B^*$ in the SGP4 convention, derives its relation to the physical ballistic coefficient $C_D A/(2m)$, and classifies the resulting measurement and accuracy errors.

Stub: Definition 32.4.1 ($B^*$): the SGP4 drag parameter satisfying $B^* = C_D A/(2m) \cdot \rho_0$ where $\rho_0$ is the reference density. Historically derived from atmospheric drag observations; carried in the TLE as a dimensionless number in units of $a_E^{-1}$ (Earth radii⁻¹). Conversion formula. [M.32.1] $B^*$ inherits $\sigma_m$ from TLE parsing (Ch 31, §31.2); 5 significant digits. [A.32.3] $B^*$ is a fitted parameter absorbing solar flux, attitude, and area/mass ratio variations; it is not a physical constant.

### §32.5 Perigee Altitude Thresholds and s* Adjustment

This section derives the three perigee altitude thresholds from the atmospheric density model and computes their effect on the fitted density parameters $s^*$ and $(q_0 - s^*)^4$.

**Example 32.5.1** (threshold effect): For perigee altitude 100 km ($r_{\mathrm{perigee}} = 6471$ km), compute $s^* = r_{\mathrm{perigee}} - 78$ km $= 22$ km$/a_E$ and the resulting $(q_0 - s^*)^4$; compare to the standard $s^* = 78$ km case. Source: SGP4 density model from Ch 21, WGS72 constants from Appendix A.

Stub: The drag parameter $s^*$ controls the reference atmospheric density layer. Three regimes: $r_{perigee} < 98$ km: $s^* = 20/a_E$; $98 \leq r_{perigee} < 156$ km: $s^* = r_{perigee} - 78$ km (interpolated); $r_{perigee} \geq 156$ km: $s^* = 78$ km/$a_E$ (standard). Above 220 km: simplified drag mode. Derive each threshold from the atmospheric density model (Ch 21, §21.2). Effect on $q_0 - s^*$ and $(q_0 - s^*)^4$: quantify the change across each threshold. [A.32.4] The 98, 156, 220 km values are empirical thresholds from the Hoots and Roehrich (1980) model, not derived from first principles.

### §32.6 Cube Root via Newton Iteration

This section specifies the Newton iteration for the cube root, proves quadratic convergence, and bounds the precision error of the final iterate.

Stub: $a_1 = (k_e/n_o)^{2/3}$ requires a cube root. Definition 32.6.1 (Newton cube-root iteration): $x_{k+1} = \tfrac{1}{3}(2x_k + c/x_k^2)$ for $x^3 = c$. Theorem 32.6.1 (quadratic convergence): the iteration converges quadratically; from a starting value within factor 2 of the true root, 5 iterations achieve double precision. — *Proof approach: compute the error $e_k = x_k - c^{1/3}$; substitute $x_k = c^{1/3}(1 + \epsilon_k)$ into the Newton formula; Taylor expand to show $\epsilon_{k+1} = \epsilon_k^2 + O(\epsilon_k^3)$, establishing quadratic convergence with rate constant 1.* Starting value: $x_0 = c^{1/3}$ approximated by manipulating the IEEE 754 exponent (bit manipulation for fast starting value). [P.32.2] Starting value bias introduces one extra iteration; precision error from final iterate $< \epsilon_{\mathrm{mach}}$.

### §32.7 Generalization: Higher-Order Series Inversion

This section presents two generalization paths for element recovery: the full generalized binomial series and direct Newton iteration on the mean-motion equation, both achieving arbitrary precision.

Stub: The third-order expansion in §32.2 is the SGP4 standard. The exact inversion is $a_o = a_1(1 + \delta_1)^{-2/3}$, which can be evaluated to arbitrary precision using the generalized binomial series (Ch 5, §5.6) with $\alpha = -2/3$ and $x = -\delta_1$. Alternatively, Newton iteration on $a^{3/2} n_o'' = k_e$ converges to full precision in 3–4 steps. This is the Generalization path for element recovery.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Theorem 1.2.4 | §32.1 | TrackedValue and precision error framework |
| Ch 3, Theorem 3.2.1 | §32.1 | Matched pair: $k_2$ must be the SGP4 value |
| Ch 5, Theorem 5.6.1 | §32.7 | Generalized binomial series for higher-order inversion |
| Ch 16, Brouwer secular rates | §32.2 | Secular $n$-correction that defines $\delta_1$ |
| Ch 21, atmospheric density model | §32.5 | Density model for perigee threshold derivation |
| Ch 31, TLE parsing | §32.2, §32.4 | Parsed $n_o, e_o, i_o, B^*$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 33, secular update | §32.2 | Recovered $a_o'', n_o'', e_o, B^*$ for propagation |
| Ch 34, near-space pipeline | §32.2, §32.5 | Epoch elements and $s^*$ parameter |
| Ch 35, deep-space pipeline | §32.2, §32.5 | Epoch elements and $s^*$ parameter |
| Appendix A | §32.2 | $k_e, k_2, a_E$ constants |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.32.1] | M | §32.4 | $B^*$ precision from TLE; 5 significant digits |
| [P.32.1] | P | §32.6 | Cube root via Newton introduces precision error $< \epsilon_{\mathrm{mach}}$ |
| [P.32.2] | P | §32.6 | Starting value selection for Newton iteration |
| [A.32.1] | A | §32.1 | Recovery procedure matched-pair dependent |
| [A.32.2] | A | §32.3 | Third-order truncation in $a_o$ expansion; bound on fourth-order term |
| [A.32.3] | A | §32.4 | $B^*$ is a fitted parameter, not a physical constant |
| [A.32.4] | A | §32.5 | Altitude thresholds (98, 156, 220 km) are empirical |
| [A.32.5] | A | §32.2 | The Kozai→Brouwer mean motion conversion (solving $a_o''^3 n_o''^2 = \mu$ after removing the $J_2$ secular correction) is a mandatory initialization step; omitting it or applying it incorrectly introduces an $O(J_2) \approx 10^{-3}$ systematic error in $a_o''$ that propagates into all secular rates; documented in Vallado et al. (2006) AIAA-2006-6753 |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 4 |
| Theorems | 2 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 1 |
| Examples | 1 |
| Error Notes | 8 |
| Equations | ~15 |
| Sections | 7 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §32.1 | Draft | Introduction |
| §32.2 | Draft | Recovery iteration |
| §32.3 | Draft | Derivation of the 134/81 coefficient |
| §32.4 | Draft | B* and the ballistic coefficient |
| §32.5 | Draft | Perigee altitude thresholds and s* adjustment |
| §32.6 | Draft | Cube root via Newton iteration |
| §32.7 | Draft | Generalization: higher-order series inversion |
