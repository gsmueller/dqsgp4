# Draft Plan: Chapter 35 — Deep-Space Propagation

**Part IX: The SGP4 Propagator** | Implementation file: `deep_space.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $T_{\mathrm{period}}$ | Orbital period in minutes | Ch 34 |
| $\dot{M}_{\odot}, \dot{\omega}_{\odot}, \dot{\Omega}_{\odot}$ | Solar secular perturbation rates | Ch 27 |
| $\dot{M}_{\leftmoon}, \dot{\omega}_{\leftmoon}, \dot{\Omega}_{\leftmoon}$ | Lunar secular perturbation rates | Ch 27 |
| $\Delta M_{\mathrm{sol}}, \Delta \omega_{\mathrm{sol}}, \ldots$ | Periodic solar corrections (DPPER outputs) | §35.5 |
| $h_{\mathrm{int}}$ | Deep-space integrator step size: 720 min | §35.6 |
| $t_{\mathrm{epoch}}$ | Epoch restart time for integrator | §35.6 |
| $\lambda_{\mathrm{res}}$ | Resonant phase angle (from Ch 28) | Ch 28 |
| $\mathrm{DPINIT}$ | Deep-space initialization routine | §35.2 |
| $\mathrm{DPSEC}$ | Deep-space secular accumulation routine | §35.3 |
| $\mathrm{DPPER}$ | Deep-space periodic correction routine | §35.5 |
| ZF, SES, SIS, SLS, SGHS, SHS | Intermediate quantities in DPPER | §35.5 |

---

## Objectives

1. Establish the $T_{\mathrm{period}} \geq 225$ min classification and its physical basis.
2. Derive the DPINIT initialization: third-body secular rates and resonance coefficient computation.
3. Derive DPSEC: accumulation of gravity and third-body secular rates.
4. Derive resonance detection condition and the Lyddane modification for near-zero inclination in deep space.
5. Derive DPPER: the periodic solar and lunar corrections, identifying each intermediate quantity.
6. Describe the 720-min secular integrator: derivation of step size, epoch restart, energy conservation check.
7. State Framework in both forms.

## Section Structure

### §35.1 Deep-Space Classification

This section defines the deep-space orbit class and its two subclasses (resonant and non-resonant), justifying the classification from perturbation magnitude analysis.

**Definition 35.1.1** (Deep-space orbit): An orbit is classified as deep-space if $T_{\mathrm{period}} \geq 225$ min ($a \gtrsim 12200$ km).

**Definition 35.1.2** (Resonant subclass): A deep-space orbit is resonant if its mean motion falls within the synchronous band $0.0034906585 \leq n_o'' \leq 0.0052359877$ rad/min (GEO) or the half-day band (GPS/Molniya).

Physical basis: at these altitudes the lunar and solar third-body perturbations accumulate significantly over one period; the secular and periodic corrections are no longer small corrections but become the dominant dynamics between Kepler steps. [A.35.1] The 225-min boundary is the same threshold as Ch 34; together they partition the orbit population. Two subclasses: resonant (GEO, GPS, Molniya) and non-resonant (HEO, highly elliptical).

### §35.2 DPINIT: Initialization

This section derives the deep-space initialization procedure, computing third-body secular rates and resonance coefficients from the epoch elements.

Stub: Compute lunar and solar secular rate contributions $(\dot{M}, \dot{\omega}, \dot{\Omega})_{\odot}$ and $(\dot{M}, \dot{\omega}, \dot{\Omega})_{\leftmoon}$ from Ch 27 formulas using epoch elements. Resonance detection: if $0.0034906585 \leq n_o'' \leq 0.0052359877$ (synchronous resonance band), or half-day resonance band: set resonance flag and compute $K, \lambda_{nm}$ coefficients from Ch 28. Store all DPINIT results as precomputed constants (Ch 37). List each quantity derived and its cross-reference. [A.35.2] Lunar and solar positions evaluated at epoch, not updated during propagation; secular rate error grows with $\Delta t$.

### §35.3 DPSEC: Secular Rate Accumulation

This section specifies the secular accumulation step combining gravity and third-body rates, distinguishing the resonant from the non-resonant paths.

Stub: At each propagation step, accumulate secular contributions: $M, \omega, \Omega$ updated by gravity rates (Ch 33) plus third-body rates (DPINIT). For resonant orbits: also accumulate the resonant drive from Ch 28. DPSEC is an inner loop within the deep-space secular integrator (§35.6). Output: updated $(M, \omega, \Omega, e, a, i)$ at the integration time.

**Example 35.3.1** (GEO secular accumulation): GPS Block IIF TLE (NORAD 37753, $a \approx 26560$ km, $e \approx 0.005$, $i \approx 55°$), propagated 7 days. Gravity secular rate $\dot{\Omega}_0 \approx -0.04°/\text{day}$; lunar secular contribution $\dot{\Omega}_{\leftmoon} \approx -0.003°/\text{day}$; solar $\dot{\Omega}_{\odot} \approx -0.001°/\text{day}$. Verify accumulated $\Omega$ after 7 days against full pipeline output.

### §35.4 Resonance Detection and Lyddane Modification

This section derives the resonance detection conditions and the Lyddane variable coordinate change that removes the $i \to 0$ singularity in deep-space secular accumulation.

**Theorem 35.4.1** (Lyddane regularity): The Lyddane variables $(L_1, L_2, L_3, L_4)$ remove the singularity at $i = 0$ in the deep-space secular equations, and the inverse transformation $L \to (\Omega, \omega, M)$ is well-defined for all $i > 0$. — *Proof approach: compute the Jacobian $\partial(L)/\partial(\Omega, \omega, M)$ explicitly and show it is nonsingular for $i > 0$; at $i = 0$ the original variables are singular but the $L$-variables remain bounded by construction — verify analytically that no $1/\sin i$ terms survive in the Lyddane form.*

Stub: Resonance flag set in DPINIT (§35.2). When active: the secular integrator uses the pendulum-equation drive (Ch 28, §28.4) in addition to the linear secular rates. Lyddane modification for deep-space near-zero inclination: replace $(\Omega, \omega, M)$ with Lyddane variables $(L_1, L_2, L_3, L_4)$ to avoid singularity. Derive the Lyddane variable form of DPSEC. [A.35.3] Lyddane modification introduces a coordinate change; verify the inverse transformation is correct.

### §35.5 DPPER: Periodic Corrections

This section derives each DPPER intermediate quantity from the periodic part of the third-body disturbing function and tables their formulas and derivation cross-references.

Stub: DPPER applies the periodic (non-secular) lunar and solar corrections to the deep-space elements. Intermediate quantities: ZF (lunar eccentricity argument), SES (solar eccentricity sine term), SIS (solar inclination sine term), SLS (solar longitude correction), SGHS (solar-lunar coupling angle), SHS (solar-hemisphere correction). Derive each from the periodic part of the third-body disturbing function (Ch 27, §27.7) specialized to the long-period terms. Table 35.5.1: each intermediate quantity, its formula, and its derivation cross-reference. [A.35.4] Periodic corrections derived at the current lunar/solar positions, which are themselves computed from the ephemeris approximations of Ch 25–26; accuracy limited by those approximations.

### §35.6 The 720-Minute Secular Integrator

This section specifies the fixed-step secular integrator, derives the epoch restart protocol, and proves that the 720-minute step satisfies the libration period requirement for all physical resonant configurations.

**Definition 35.6.1** (epoch restart): integration begins fresh at $t_o$ with the epoch initial conditions, then steps forward to time $t$ using DPSEC with step $h = 720$ min.

**Lemma 35.6.1** (step size requirement): $h = 720$ min $\ll T_{\mathrm{libration}}$ for all physical resonant orbits. — *Proof approach: from the Ch 28 GEO libration period $T_{\mathrm{lib}} \approx 840$ days $\approx 1.21\times10^6$ min; $h/T_{\mathrm{lib}} \approx 6\times10^{-4}$, satisfying the Nyquist condition $h < T_{\mathrm{lib}}/10$ for GPS and GEO resonances, and the GPS libration period is even larger.*

Stub: Deep-space satellites propagate by integrating the secular equation of motion (DPSEC) numerically, using a fixed step of $h = 720$ min (12 hours) starting from epoch. At each target time $t$: integrate forward or backward to the nearest multiple of 720 min, then step to $t$. Epoch restart: when the target $t$ changes sign (crossing epoch), restart integration from epoch. Definition 35.6.1 (epoch restart): integration begins fresh at $t_o$ with the epoch initial conditions, then steps to $t$. Energy conservation check: $|n^2 a^3 - \mu| < \epsilon$ as a diagnostic. Lemma 35.6.1 (step size requirement): $h \ll T_{\mathrm{libration}}$ from Ch 28 for resonant orbits; 720 min satisfies this for all physical resonant configurations. [P.35.1] Fixed-step integrator; no adaptive error control; accumulates $O(h^2)$ precision error per step.

### §35.7 State Framework: Deep-Space Corrections in Both Forms

This section expresses the DPPER periodic corrections and DPSEC secular accumulation in both the dual-quaternion and 7×7 matrix representations, placing them within the parallel-row error budget.

Stub: DPPER corrections are additive perturbations $(\Delta M, \Delta\omega, \Delta\Omega, \Delta e, \Delta i)$. **Dual-quaternion form:** each additive correction to a Keplerian element corresponds to a near-identity perturbation $\delta\hat{\Omega}_b$ to the velocity dual quaternion (Principle 2, Ch 2, §2.6). Composition: the corrected $(\hat{M}_b, \hat{\Omega}_b)$ is formed by left-multiplying $\hat{\Omega}_b$ by $(I + \delta\hat{\Omega}_b)$. **7×7 matrix form:** additive corrections form near-identity 7×7 matrices; composition by matrix multiplication. Parallel-row error budget: each correction adds to $\delta_a$ (model error from the third-body approximation) and $\delta_p$ (from the integration step).

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 2 | State framework, parallel-row architecture | Dual-quaternion and 7x7 forms for corrections |
| Ch 10 | Lyddane variables | Kepler step variable formulation |
| Ch 25–26 | Solar and lunar ephemerides | Epoch positions for third-body rates |
| Ch 27 | Third-body perturbation rates | DPINIT, DPSEC, DPPER source formulas |
| Ch 28 | Resonance coefficients | Pendulum-equation drive for resonant orbits |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 34 | Near-space pipeline | Complementary path for $T < 225$ min |
| Ch 37 | Precomputed constants | DPINIT outputs stored as precomputed |
| Ch 38 | Output error budget | Deep-space pipeline error contributions |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.35.1] | P | §35.6 | Fixed-step integrator accumulates $O(h^2)$ precision error per step |
| [A.35.1] | A | §35.1 | 225-min threshold empirical |
| [A.35.2] | A | §35.2 | Third-body rates from epoch positions; secular drift not updated |
| [A.35.3] | A | §35.4 | Lyddane inverse transformation must be verified |
| [A.35.4] | A | §35.5 | DPPER accuracy limited by ephemeris approximations of Ch 25–26 |
| [A.35.5] | A | §35.6 | SR3 deep-space integrator has a step-direction sign error for backward propagation, producing a jump in semi-major axis near $t = -720$ min for GEO and Molniya orbits; corrected by Vallado et al. (2006) AIAA 2006-6753 (tagged `sgp4fix` in source) |
| [A.35.6] | A | §35.5 | SR3 DPPER skips lunar-solar recomputation if $|\Delta t| < 30$ min (SAVTSN shortcut), producing a discontinuous ephemeris for small time steps; shortcut removed in Vallado et al. (2006) |
| [A.35.7] | A | §35.2 | SR3 DPINIT initializes deep-space resonance terms using unperturbed element values; correct behavior uses element values updated by deep-space perturbations; corrected in Vallado et al. (2006) |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 5 |
| Theorems | 1 |
| Lemmas | 2 |
| Corollaries | 0 |
| Propositions | 1 |
| Examples | 2 |
| Error Notes | 8 |
| Equations | ~30 |
| Sections | 7 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §35.1 Deep-Space Classification | Draft | |
| §35.2 DPINIT: Initialization | Draft | |
| §35.3 DPSEC: Secular Rate Accumulation | Draft | |
| §35.4 Resonance Detection and Lyddane Modification | Draft | |
| §35.5 DPPER: Periodic Corrections | Draft | |
| §35.6 The 720-Minute Secular Integrator | Draft | |
| §35.7 State Framework: Deep-Space Corrections in Both Forms | Draft | |
