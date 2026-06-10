# Draft Plan: Chapter 34 — Near-Space Propagation

**Part IX: The SGP4 Propagator** | Implementation file: `near_space.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $T_{\mathrm{period}}$ | Orbital period in minutes | §34.1 |
| $a_t, e_t$ | Updated semi-major axis and eccentricity at time $t$ | Ch 33 |
| $M_p, \omega_p, \Omega_p$ | Mean elements after secular update | Ch 33 |
| $e_L, \omega_L, i_L, M_L, \Omega_L$ | Elements after long-period correction | §34.3 |
| $E_L$ | Eccentric anomaly from modified Kepler equation | Ch 10 |
| $r_k, \dot{r}_k, (r\dot{f})_k$ | Osculating radius, radial and transverse velocity | Ch 20 |
| $\hat{M}_b, \hat{\Omega}_b$ | Configuration and velocity dual quaternions in TEME | Ch 2 |
| $\mathbf{X}_{7\times7}$ | State in homogeneous 7×7 form | Ch 2 |

---

## Objectives

1. Define the near-space classification ($T_{\mathrm{period}} < 225$ min) and its physical basis.
2. Specify the initialization procedure: which constants are precomputed from the TLE, and which are computed per propagation step.
3. Describe the five-stage propagation pipeline as a composition of chapters.
4. Characterize simplified drag mode: which pipeline stages are affected.
5. State Framework: the full pipeline in both dual-quaternion and 7×7 forms, with parallel-row error budget.

## Section Structure

### §34.1 Near-Space Classification

This section formally defines the near-space orbit class and justifies the 225-minute threshold from perturbation magnitude analysis.

**Definition 34.1.1** (Near-space orbit): An orbit is classified as near-space if $T_{\mathrm{period}} = 2\pi/n_o'' < 225$ min ($a \lesssim 12200$ km). Physical basis: orbits below this threshold have insufficient period for third-body (Moon, Sun) perturbations and tesseral resonance to be significant.

**Proposition 34.1.1** (Threshold justification): For $T < 225$ min, the magnitude of the lunar secular perturbation rate $|\dot{\Omega}_{\leftmoon}|$ is less than $10^{-7}$ rad/min, contributing $< 0.1$ m position error per orbit. — *Proof approach: bound the third-body perturbation magnitude using the Ch 27 secular rate formulas evaluated at $a = 12200$ km, showing the ratio to the $J_2$ secular rate is $O(10^{-4})$.*

[A.34.1] The 225-min threshold is empirical; it demarcates where deep-space corrections are numerically small but not rigorously negligible.

### §34.2 Initialization: Precomputed Constants

This section enumerates all constants computed once from the TLE at initialization and establishes the initialization-vs-propagation boundary for the near-space path.

**Definition 34.2.1** (Near-space precomputed set): The set $\mathcal{P}_{\mathrm{near}} \subset \mathcal{P}_{\mathrm{sat}}$ consisting of all quantities computed once from TLE epoch elements for near-space propagation.

Stub: List all constants computed once from the TLE at initialization:

1. Element recovery group: $n_o'', a_o'', s^*, (q_0-s^*)^4$ (Ch 32)
2. Drag coefficient group: $C_1, C_2, C_3, C_4, C_5, D_2, D_3, D_4$ (Ch 22)
3. Brouwer secular rates: $\dot{M}_0, \dot{\omega}_0, \dot{\Omega}_0$ (Ch 16--17)
4. Short-period coefficient group (Ch 18)
5. Long-period coefficient group (Ch 19)
6. Kaula inclination functions $F_{lmp}(i_o)$ for $l \leq 5$ (Ch 15)

The initialization-vs-propagation boundary: all quantities depending only on $a_o'', e_o, i_o$ are initialized; all quantities depending on $t$ (through $a_t, e_t$, or angles) are computed per step. Cross-reference Ch 37 (full precomputed constant taxonomy).

**Example 34.2.1** (ISS initialization): For ISS TLE (NORAD 25544, $i \approx 51.6°$, $e \approx 0.0001$, $n \approx 15.5$ rev/day): $a_o'' \approx 6796$ km, $(q_0 - s^*)^4 \approx 1.19 \times 10^{-12}$, $C_1 \approx 2.1 \times 10^{-11}$, $\dot{\Omega}_0 \approx -5.06°/\text{day}$. Verify these against Ch 32 element recovery and Ch 16 secular rate derivations.

### §34.3 Stage 1: Secular Update

This section specifies the secular update stage as a direct evaluation of Ch 33 formulas, producing time-dependent mean elements $(a_t, e_t, M_p, \omega_p, \Omega_p)$.

Stub: Call Ch 33 secular update. Output: $(a_t, e_t, M_p, \omega_p, \Omega_p)$. This is purely arithmetic; no iteration required.

### §34.4 Stage 2: Long-Period Correction

This section applies Brouwer long-period corrections as near-identity perturbations to the secular mean elements, producing the long-period-corrected element set $(e_L, \omega_L, i_L, M_L, \Omega_L)$.

Stub: Apply Brouwer long-period corrections from Ch 19 (and Ch 18 for the $\omega$ correction). Output: $(e_L, \omega_L, i_L, M_L, \Omega_L)$. Note that in simplified drag mode ($< 220$ km perigee) the long-period corrections are still applied. [A.34.2] Long-period corrections derived only to first order in $J_3/J_2$; higher-order terms dropped.

### §34.5 Stage 3: Modified Kepler Equation

**Objective.** Solve the modified Kepler equation for eccentric anomaly and derive osculating orbital quantities, characterizing Newton iteration convergence.

Stub: Solve the modified Kepler equation (Ch 10) for eccentric anomaly $E_L$ given $M_L, e_L, \omega_L$ in Lyddane variables $(a_{xN} = e_L\cos\omega_L, a_{yN} = e_L\sin\omega_L + $ long-period $e_{LL}$ correction). Newton iteration converges in 3–5 steps to double precision. Output: $E_L$, and the derived osculating $(r_k, \dot{r}_k, (r\dot{f})_k)$ via Ch 20.

**Example 34.5.1** (Kepler convergence for ISS): ISS TLE propagated to $t = 120$ min. With $e_L \approx 0.0001$, $M_L \approx 3.14$ rad: Newton iteration converges in 3 steps to $|E_{k+1} - E_k| < 10^{-12}$. Compare resulting $(r_k, \dot{r}_k, (r\dot{f})_k)$ against Ch 20 formulas; verify position radius $r_k \approx 6796$ km.

### §34.6 Stage 4: Short-Period Correction

This section applies Brouwer short-period corrections as additive near-identity perturbations to the osculating orbital quantities and quantifies the first-order truncation error.

Stub: Apply Brouwer short-period corrections from Ch 18 to the osculating orbital quantities. Output: $(r, u, \Omega_{\mathrm{osc}}, i_{\mathrm{osc}}, \dot{r}, r\dot{f})$. [A.34.3] Short-period corrections to first order in $J_2$; second-order Lara (2021) terms provide ~1 m improvement for LEO.

### §34.7 Stage 5: State Vector Construction

This section constructs the TEME position and velocity state vector from the osculating orbital quantities via the Ch 30 coordinate transformation, completing the near-space pipeline.

Stub: Call Ch 30 to construct TEME position and velocity from $(r, u, \Omega_{\mathrm{osc}}, i_{\mathrm{osc}}, \dot{r}, r\dot{f})$. Output: $\mathbf{r}_b, \mathbf{v}_b$ in TEME.

### §34.8 State Framework: Full Pipeline in Both Forms

This section expresses the complete five-stage near-space pipeline as a composition of near-identity transformations in both the dual-quaternion and 7×7 matrix representations, with parallel-row error propagation.

**Theorem 34.8.1** (Pipeline composition): The near-space propagator is a composition of five stages, each a near-identity perturbation, whose total error decomposes as the sum of per-stage error contributions in each of the three error categories $(\sigma_m, \delta_p, \delta_a)$. — *Proof approach: model each stage as a near-identity perturbation adding $O(\epsilon_k)$ to the state; apply the triangle inequality sequentially through all five stages using the three-error decomposition theorem (Ch 1, Thm 1.3.1) to bound the cumulative error.*

Stub: **Dual-quaternion form.** Each pipeline stage is a transformation on the pair $(\hat{M}_b, \hat{\Omega}_b)$. Stages 1–2 update the mean-element representation (near-identity perturbations to the velocity dual quaternion $\hat{\Omega}_b$). Stage 3 is the Kepler step: a near-identity $\hat{M}_b$ update representing one orbit. Stage 4 adds short-period perturbations to $\hat{M}_b$. Stage 5 is the final frame composition giving the TEME state. The composition is left-to-right multiplication: $(\hat{M}, \hat{\Omega})_{\mathrm{TEME}} = $ (frame rotation) $\circ$ (Kepler) $\circ$ (corrections) $\circ$ (epoch state). **7×7 matrix form.** Each stage corresponds to a 7×7 matrix acting on $(r, v, 1)^T$; corrections are additive perturbation matrices (near-identity). The complete pipeline is a product of 7×7 matrices. **Parallel-row error budget** (Ch 2, §2.10): the error state is a 4-row structure in each form; each stage propagates $(σ_m, δ_p, δ_a)$ per Ch 1 rules.

### §34.9 Simplified Drag Mode: Pipeline Changes

This section documents precisely which terms are dropped in simplified drag mode and quantifies the resulting accuracy cost relative to the standard pipeline.

**Proposition 34.9.1** (Simplified mode accuracy bound): In simplified drag mode ($h_p < 220$ km), the position error attributable to term truncation is bounded by $\delta_a^{\mathrm{simp}} \leq C \cdot B^* \cdot (a_o'' - a_E)^2 \cdot \Delta t$ where $C$ depends on the dropped higher-order drag terms. — *Proof approach: bound the magnitude of each dropped $D_k$ term from Ch 33 §33.7 by evaluating the ratio $D_{k+1}(\Delta t)^{k+2}/D_k(\Delta t)^{k+1} \approx C_1\Delta t$ and using the power-law atmosphere model (Ch 21) to estimate $C_1$ at perigee altitude 220 km.*

**Example 34.9.1** (Simplified vs. standard for low-perigee satellite): Satellite with $h_p = 180$ km, $B^* = 0.01$, propagated 24h. Compare position output from standard pipeline vs. simplified mode; expected difference $\sim 50$--$200$ m depending on drag coefficient, dominated by truncated $D_3, D_4$ terms.

Stub: When perigee < 220 km: Stage 1 uses truncated $a(t)$ (Ch 33, §33.7); Stage 3 uses fewer Kepler iterations (same algorithm, different entry values); all other stages unchanged. Document precisely which terms are dropped and quantify the accuracy cost. [A.34.4] Simplified mode sacrifices accuracy below 220 km in exchange for numerical stability in rapidly changing $a_t$.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 10 | Modified Kepler equation | Eccentric anomaly solution in Stage 3 |
| Ch 18 | Short-period corrections | Stage 4 short-period correction formulas |
| Ch 19 | Long-period corrections | Stage 2 long-period correction formulas |
| Ch 20 | Osculating orbital quantities | Osculating $(r_k, \dot{r}_k, (r\dot{f})_k)$ from Kepler step |
| Ch 30 | State vector from elements | Stage 5 TEME coordinate construction |
| Ch 33 | Secular update | Stage 1 secular element propagation |
| Ch 37 | Precomputed constants taxonomy | Initialization-vs-propagation boundary |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 35 | Deep-space classification | Complementary path for $T \geq 225$ min |
| Ch 38 | Output error budget | Near-space pipeline error contributions |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.34.1] | A | §34.1 | 225-min threshold empirical; demarcation not rigorously negligible |
| [A.34.2] | A | §34.4 | Long-period corrections first-order in $J_3/J_2$ only |
| [A.34.3] | A | §34.6 | Short-period corrections first-order in $J_2$ only |
| [A.34.4] | A | §34.9 | Simplified drag mode sacrifices accuracy below 220 km |
| [A.34.5] | A | §34.5 | SR3 Kepler iteration fails for very high eccentricity orbits (demonstrated on satellite 23333, WIND spacecraft); produces dramatic jumps in computed inclination around 200 min propagation time; corrected by Vallado et al. (2006) AIAA-2006-6753 with an improved convergence strategy |
| [A.34.6] | A | §34.8 | Lara et al. (2023, arXiv:2307.06864) found two discrepancies between SGP4 equations and the analytical theory: (a) slightly different coefficient on the linear secular term in the along-track angle; (b) a missing short-period sinusoidal term with orbital period; the secular discrepancy is suspected to be absorbed into $B^*$ during TLE fitting |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 3 |
| Theorems | 1 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 2 |
| Examples | 3 |
| Error Notes | 6 |
| Equations | ~20 |
| Sections | 9 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §34.1 Near-Space Classification | Draft | |
| §34.2 Initialization: Precomputed Constants | Draft | |
| §34.3 Stage 1: Secular Update | Draft | |
| §34.4 Stage 2: Long-Period Correction | Draft | |
| §34.5 Stage 3: Modified Kepler Equation | Draft | |
| §34.6 Stage 4: Short-Period Correction | Draft | |
| §34.7 Stage 5: State Vector Construction | Draft | |
| §34.8 State Framework: Full Pipeline in Both Forms | Draft | |
| §34.9 Simplified Drag Mode: Pipeline Changes | Draft | |
