# Draft Plan: Chapter 20 — Osculating Orbital Quantities

**Part IV: Brouwer's Gravitational Perturbation Theory**

Target header: `osculating_elements.h`

---

## Objectives

1. Convert the mean elements (after all secular and periodic corrections from Ch 16–19) into osculating orbital quantities: radius $r$, radial velocity $\dot{r}$, transverse velocity $r\dot{f}$, and argument of latitude $u$.
2. Derive each quantity from the modified Kepler solution of Ch 10, treating $E$ as a function of the corrected mean anomaly and the Lyddane variables.
3. State the complete mean-to-osculating pipeline in both representations: dual quaternion composition on $(\hat{M}, \hat{\Omega})$ and matrix product on the $7\times7$ state.
4. Propagate the parallel-row error structure through the final assembly step.

---

## Section Structure

### §20.1 Introduction

This section summarizes the complete mean-to-osculating pipeline that assembles corrections from Chs 16--19 into physical position-velocity state vectors.

*Stub.* Summarize the role of this chapter as the terminal assembly step of Part IV. The secular corrections of Ch 16–17, the short-period corrections of Ch 18, and the long-period corrections of Ch 19 have each modified the mean element set. This chapter applies the corrected elements to obtain quantities that can be assembled into a physical position-velocity state.

Forward references: Ch 8 (Keplerian geometry), Ch 9 (Kepler's equation), Ch 10 (modified Kepler equation and Lyddane variables), Ch 18 (short-period corrections), Ch 19 (long-period corrections).

---

### §20.2 The Modified Kepler Solution

This section sets up the Kepler equation in Lyddane variables and defines the intermediate angle $U = E + \omega$ used in all subsequent osculating quantity formulas.

*Stub.* Recall from Ch 10 the formulation of Kepler's equation in Lyddane variables $(a_{xN}, a_{yN}) = (e\cos\omega, e\sin\omega)$ that avoids the singularity at small $e$ and small $i$. The eccentric anomaly $E$ is the solution variable; define the intermediate angle $U = E + \omega$ used throughout SGP4.

**Notation table for this section:**

| Symbol | Meaning |
|--------|---------|
| $M'$ | Corrected mean anomaly (after Ch 16–19) |
| $e$ | Corrected eccentricity (after long-period corrections) |
| $\omega$ | Corrected argument of perigee (after long-period corrections) |
| $a_{xN}$ | $e\cos\omega$ (Lyddane $x$-component) |
| $a_{yN}$ | $e\sin\omega$ (Lyddane $y$-component) |
| $E + \omega$ | Combined eccentric anomaly plus argument of perigee |
| $\ell$ | $M' + \omega + \Omega$ (mean longitude, Lyddane variable) |

**Definition 20.2.1** (Lyddane mean longitude). *Stub — formal definition.*

**Proposition 20.2.1** (Newton iteration in Lyddane variables). *Stub -- convergence in Lyddane variables from Ch 10.* — *Proof approach: apply the Newton-Raphson convergence theorem from Ch 10 to the modified Kepler equation $U - a_{xN}\sin U + a_{yN}\cos U = \ell$, with convergence guaranteed by the contraction condition $|a_{xN}| + |a_{yN}| < 1$*

---

### §20.3 Radius and Radial Velocity

This section derives the osculating radius and radial velocity from the Lyddane-variable Kepler solution, with error analysis for the near-circular case.

*Stub.* From the solution $U = E + \omega$ of §20.2, derive:

$$r = a(1 - e\cos E) \tag{20.1}$$

Expand $\cos E = \cos(U - \omega) = \cos U \cos\omega + \sin U \sin\omega$ and express in terms of $a_{xN}, a_{yN}$. Derive:

$$\dot{r} = \frac{\sqrt{\mu a}}{r}e\sin E \tag{20.2}$$

in terms of the Lyddane variables. State the error propagation: $r$ accumulates precision error from the Kepler solution and from $a$, $e$ (which carry corrections from Ch 16–19). Identify [P.20.1]: subtractive cancellation in $1 - e\cos E$ for near-circular orbits.

**Theorem 20.3.1** (Radius from Lyddane variables). *Stub -- express $r$ directly in $(a_{xN}, a_{yN}, U)$: $r = a(1 - a_{xN}\cos U - a_{yN}\sin U)$.* — *Proof approach: substitute $e\cos E = a_{xN}\cos(U-\omega) + a_{yN}\sin(U-\omega) = a_{xN}\cos U + a_{yN}\sin U$ (using $U = E + \omega$) into $r = a(1 - e\cos E)$; the key identity is verified by expanding $\cos(U - \omega)$ and $\sin(U - \omega)$*

---

### §20.4 Transverse Velocity

This section expresses the transverse velocity component in Lyddane variables via the angular momentum magnitude.

*Stub.* The transverse velocity (sometimes called the tangential component in the radial-transverse-normal frame) is:

$$r\dot{f} = \frac{\sqrt{\mu p}}{r} \tag{20.3}$$

where $p = a(1 - e^2)$ is the semi-latus rectum. Express in Lyddane variables. Derive $\dot{f}$ as a function of $E$. Note the relation to the angular momentum magnitude $h = \sqrt{\mu p}$.

**Proposition 20.4.1** (Transverse velocity from angular momentum). *Stub: $r\dot{f} = h/r = \sqrt{\mu p}/r$ where $p = a(1 - a_{xN}^2 - a_{yN}^2)$.* — *Proof approach: use the vis-viva relation and the angular momentum conservation $h = \sqrt{\mu p}$, substituting $e^2 = a_{xN}^2 + a_{yN}^2$*

---

### §20.5 Argument of Latitude

This section recovers the argument of latitude $u = \omega + \nu$ from the Lyddane Kepler solution, using atan2 for quadrant-correct evaluation.

*Stub.* The argument of latitude $u = \omega + \nu$ where $\nu$ is the true anomaly. Derive $\sin u$ and $\cos u$ from the eccentric anomaly solution via:

$$\cos\nu = \frac{\cos E - e}{1 - e\cos E}, \qquad \sin\nu = \frac{\sqrt{1-e^2}\sin E}{1 - e\cos E} \tag{20.4}$$

Then $u$ is recovered by `atan2(sin u, cos u)`. In the Lyddane formulation, express $\sin u$ and $\cos u$ directly in terms of $(a_{xN}, a_{yN}, U)$ to avoid separate computation of $\nu$.

**Proposition 20.5.1** (Argument of latitude in Lyddane variables). *Stub -- algebraic form: $\sin u = (sin U - a_{yN} - a_{xN}\sin U \cdot a_{yN}\cos U / (1 + \sqrt{1-e^2}))\cdot a/r$ and analogous $\cos u$ expression in $(a_{xN}, a_{yN}, U)$; then $u = \text{atan2}(\sin u, \cos u)$.* — *Proof approach: express $\sin\nu$ and $\cos\nu$ via the eccentric anomaly relations (Ch 8) in Lyddane form, then apply the addition formula $u = \omega + \nu$ absorbed into the $U$ variable*

**Error Note [A.20.1]:** The `atan2` step is Tier I (full double precision, Lipschitz-1), so no additional accuracy error is introduced here. Precision error propagates continuously.

---

### §20.6 Application of All Corrections

This section assembles the complete mean-to-osculating pipeline, specifying the correct order of correction application and the resulting osculating state.

*Stub.* Assemble the complete pipeline:

1. Start from mean Keplerian elements $(a_0, e_0, i_0, \Omega_0, \omega_0, M_0)$.
2. Apply secular corrections (Ch 16–17): obtain secular mean elements at time $t$.
3. Apply long-period corrections (Ch 19): obtain long-period-corrected mean elements.
4. Solve the modified Kepler equation (§20.2): obtain $U$, $r$, $\dot{r}$, $r\dot{f}$, $u$.
5. Apply short-period corrections (Ch 18): add the Brouwer short-period $\Delta r$, $\Delta u$, $\Delta\Omega$, $\Delta i$, $\Delta\dot{r}$, $\Delta(r\dot{f})$ to the quantities from step 4.
6. Compute position and velocity unit vectors from $u$ and $\Omega$, $i$.

Note the order dependence: short-period corrections are applied to the post-Kepler osculating quantities, not to mean elements.

**Theorem 20.6.1** (Brouwer osculating state). *Stub -- complete expression for $(r, \dot{r}, r\dot{f}, u, \Omega, i)$ as osculating elements after all corrections:*

1. *$r_{\text{osc}} = r_{\text{Kepler}} + \Delta r^{\text{sp}}$*
2. *$u_{\text{osc}} = u_{\text{Kepler}} + \Delta u^{\text{sp}}$*
3. *$\Omega_{\text{osc}} = \bar{\Omega}_{\text{sec}} + \Delta\Omega^{\text{sp}}$*
4. *$i_{\text{osc}} = \bar{i} + \Delta i^{\text{sp}}$*
5. *$\dot{r}_{\text{osc}} = \dot{r}_{\text{Kepler}} + \Delta\dot{r}^{\text{sp}}$*
6. *$(r\dot{f})_{\text{osc}} = (r\dot{f})_{\text{Kepler}} + \Delta(r\dot{f})^{\text{sp}}$*

— *Proof approach: direct substitution of the correction formulas from Ch 18 (short-period) and Ch 19 (long-period, already applied in the mean elements) into the Kepler solution quantities from §§20.2--20.5*

---

### §20.7 Position and Velocity Vectors

This section computes the ECI position and velocity vectors from the osculating orbital quantities using the perifocal-to-inertial frame rotation.

*Stub.* From the osculating quantities, compute the ECI position and velocity vectors. Define the unit vectors in the orbital plane:

$$\hat{r} = \cos u\, \hat{q} + \sin u\, \hat{p} \tag{20.5}$$

where $\hat{p}$, $\hat{q}$ are the perifocal frame unit vectors (functions of $\Omega$ and $i$). Explicitly:

$$\vec{r}_I = r\hat{r}_I, \qquad \dot{\vec{r}}_I = \dot{r}\hat{r}_I + r\dot{f}\hat{\theta}_I \tag{20.6}$$

State which frame each vector is expressed in (inertial $I$ frame throughout this section). Note [A.20.2]: the ECI frame is a model-dependent choice; the frame definition carries accuracy error from the pole and equinox conventions (Ch 29).

**Example 20.7.1** (Complete osculating state for ISS-like orbit)**.** *Stub: for an ISS-like orbit with mean elements $a = 6780$ km, $e = 0.001$, $i = 51.6°$, $\Omega = 120°$, $\omega = 45°$, $M = 90°$ at epoch, propagated $\Delta t = 1$ hour using WGS84 ($\mu = 398600.4418$ km$^3$/s$^2$, $a_E = 6378.137$ km, $J_2 = 1.08263 \times 10^{-3}$). Execute the full pipeline:*

1. *Secular update (Ch 16--17): $\Omega_{\text{sec}} = \Omega_0 + \dot{\Omega}\Delta t$, $\omega_{\text{sec}} = \omega_0 + \dot{\omega}\Delta t$, $M_{\text{sec}} = M_0 + \dot{M}_{\text{sec}}\Delta t$*
2. *Long-period corrections (Ch 19): $\Delta e_{\text{lp}} \approx 10^{-4}$, $IL_L \approx 10^{-5}$ rad*
3. *Kepler solution (§20.2): $U$ from Newton iteration, 3--4 iterations to $10^{-15}$*
4. *Osculating quantities (§§20.3--20.5): $r_{\text{osc}} \approx 6780$ km, $\dot{r} \approx 0$ km/s, $r\dot{f} \approx 7.66$ km/s*
5. *Short-period corrections (Ch 18): $\Delta r \approx 2$ km, $\Delta u \approx 0.002$ rad*
6. *Position/velocity vectors (§20.7): $\vec{r}_I$, $\dot{\vec{r}}_I$ in TEME frame*

*Source: WGS84.*

---

### §20.8 State Framework

This section expresses the complete mean-to-osculating transformation in both the dual quaternion and 7x7 matrix forms, with parallel-row error propagation.

*Stub.* State the complete mean-to-osculating transformation in both representations:

**Dual quaternion form.** The configuration dual quaternion $\hat{M}$ encodes both the perifocal-to-inertial rotation (via the SU(2) part) and the position (via the dual part). The velocity dual quaternion $\hat{\Omega}_b$ encodes the angular velocity and translational velocity. The Brouwer corrections (Ch 18) modify $\hat{M}$ by a near-identity dual quaternion perturbation; the velocity corrections modify $\hat{\Omega}_b$.

**7×7 matrix form.** The complete transformation from mean elements to osculating state is a $7\times7$ matrix product: secular + long-period rotation, Kepler step (near-identity block), short-period additive correction assembled as a near-identity $7\times7$.

**Parallel-row error propagation.** The four-row structure (state row + three error rows) passes through each step unchanged in structure; the error rows accumulate contributions from each correction step according to the rules of Ch 1.

**Proposition 20.8.1** (Near-identity composition). *Stub -- dual quaternion near-identity perturbation from the short-period correction: $\hat{M}_{\text{osc}} = \hat{M}_{\text{mean}} \cdot (1 + \epsilon\,\delta\hat{M}^{\text{sp}})$ where $\delta\hat{M}^{\text{sp}}$ encodes the six short-period corrections.* — *Proof approach: linearize the dual quaternion composition (Ch 2, Theorem 2.3.1) to first order in $J_2$, showing that the near-identity perturbation adds the correction terms to both the rotation and translation components*

**Proposition 20.8.2** (7x7 assembly). *Stub -- equivalent matrix form: the total transformation is $(I + \delta W_{\text{sec}})(I + \delta W_{\text{lp}})(I + \delta W_{\text{Kepler}})(I + \delta W_{\text{sp}})$ acting on $(r, v, 1)^T$, where each $\delta W$ is a sparse $7\times7$ near-identity perturbation.* — *Proof approach: direct matrix composition using the algebraic equivalence of Ch 2, Theorem 2.8.1; to first order in $J_2$ the composition is additive*

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.20.1] | P | §20.3 | Subtractive cancellation in $r = a(1-e\cos E)$ for $e \ll 1$: relative error $\sim e\epsilon_{\rm mach}/(1-e)$ |
| [P.20.2] | P | §20.6 | Accumulation of correction terms from Ch 16–19 into $e$, $\omega$ before Kepler solve |
| [A.20.1] | A | §20.6 | Brouwer theory truncated at first-order short-period; second-order contributions (Lara 2021) omitted |
| [A.20.2] | A | §20.7 | ECI frame realization depends on pole/equinox conventions (Ch 29); frame error not tracked until Ch 29 |
| [M.20.1] | M | §20.6 | Osculating elements depend on TLE-fitted mean elements; TLE fit residuals propagate into all derived quantities |

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Ch 1, error propagation through arithmetic | Error propagation rules |
| Ch 2 | Ch 2, dual quaternion composition | State Framework dual quaternion form |
| Ch 2 | Ch 2, velocity dual quaternion | Velocity state representation |
| Ch 2 | Ch 2, 7×7 matrix equivalence | State Framework matrix form |
| Ch 2 | Ch 2, parallel-row error architecture | Four-row error structure |
| Ch 10 | Ch 10, Lyddane variables and modified Kepler equation | Non-singular Kepler solution |
| Ch 10 | Ch 10, Newton-Raphson convergence | Kepler iteration in Lyddane variables |
| Ch 16–17 | Ch 16–17, secular corrections | Applied before Kepler solve |
| Ch 18 | Ch 18, short-period corrections | Applied after Kepler solve |
| Ch 19 | Ch 19, long-period corrections | Applied before Kepler solve |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 29 | Ch 29, ECI frame conventions | Frame definition for position/velocity vectors |
| Ch 30 | Ch 30, coordinate transformations | Osculating quantities $r$, $\dot{r}$, $r\dot{f}$, $u$ consumed by TEME state vector |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 2 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 5 |
| Examples | 1 |
| Error Notes | 5 |
| Equations | ~10 |
| Sections | 8 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §20.1 Introduction | Draft | |
| §20.2 The Modified Kepler Solution | Draft | |
| §20.3 Radius and Radial Velocity | Draft | |
| §20.4 Transverse Velocity | Draft | |
| §20.5 Argument of Latitude | Draft | |
| §20.6 Application of All Corrections | Draft | |
| §20.7 Position and Velocity Vectors | Draft | |
| §20.8 State Framework | Draft | |
