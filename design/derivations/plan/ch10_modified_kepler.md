# Draft Plan: Chapter 10 — The Modified Kepler Equation

**Part II: The Two-Body Problem**

**Phase:** Draft

**Implementation target:** `modified_kepler.h`

---

## Objectives

1. Derive the rotated (shifted) form of Kepler's equation in which the combined angle $(E + \omega)$ is the unknown, eliminating the explicit appearance of $\omega$ in the iteration.
2. Introduce the Lyddane variables $(a_{xN}, a_{yN}) = (e\cos\omega, e\sin\omega)$ and show how they absorb the eccentricity singularity (undefined $\omega$ at $e = 0$).
3. Reformulate Newton-Raphson in Lyddane variables and derive the explicit iteration formula used in SGP4.
4. Establish convergence bounds for the Lyddane-form iteration analogous to those of Ch 9.
5. Express the Kepler propagation step in both the dual quaternion representation $(\hat{M}, \hat{\Omega})$ and the 7×7 matrix form, completing the State Framework for the two-body orbit advance.
6. Tag all precision and accuracy implications of the Lyddane formulation.

---

## Notation Table

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $M$ | Mean anomaly | Ch 8, §8.5 |
| $E$ | Eccentric anomaly | Ch 9, §9.2 |
| $\nu$ | True anomaly | Ch 8, §8.4 |
| $\omega$ | Argument of perigee | Ch 8, §8.5 |
| $e$ | Eccentricity | Ch 8, §8.4 |
| $\eta$ | Eccentricity complement: $\eta = \sqrt{1-e^2}$ | Ch 9, §9.2 |
| $U$ | Combined angle (shifted eccentric anomaly): $U = E + \omega$ | §10.2 |
| $M_L$ | Lyddane mean longitude: $M_L = M + \omega$ | §10.2 |
| $a_{xN}$ | Lyddane $x$-eccentricity variable: $a_{xN} = e\cos\omega$ | §10.3 |
| $a_{yN}$ | Lyddane $y$-eccentricity variable: $a_{yN} = e\sin\omega$ | §10.3 |
| $\mathbf{u}_L$ | Lyddane eccentricity vector: $(a_{xN}, a_{yN})$ | §10.3 |
| $e_L$ | Eccentricity recovered from Lyddane variables: $e_L = \sqrt{a_{xN}^2 + a_{yN}^2}$ | §10.3 |
| $f_L(U)$ | Lyddane Kepler residual: $f_L(U) = U - a_{xN}\sin U + a_{yN}\cos U - M_L$ | §10.4 |
| $\tau$ | Caller-supplied tolerance | Ch 3 |
| $\hat{M}_{\text{step}}$ | Configuration dual quaternion for the Kepler propagation step | §10.6 |
| $T_{\text{step}}$ | 7×7 matrix for the Kepler propagation step | §10.6 |

---

## Section Structure

### §10.1 Introduction

Narrative: Chapter 9 solves Kepler's equation in the form $M = E - e\sin E$, requiring knowledge of both $e$ and $\omega$ as separate quantities. Two difficulties arise in practice. First, when the SGP4 mean elements are updated by secular and long-period corrections (Chapters 16–20), the corrected $\omega$ must be used together with the corrected $e$; the equation can be rewritten directly in terms of the corrected combined angle. Second, near $e = 0$, the argument of perigee $\omega$ is geometrically undefined, yet the product $e\sin\omega$ appears in the equation and can be evaluated stably even as $\omega$ loses meaning individually.

The Lyddane variable approach, used in the SGP4 implementation, addresses both issues simultaneously. This chapter derives the Lyddane formulation from first principles and certifies the Newton-Raphson iteration.

Forward reference table:

| Section | Feeds | Role |
|---------|-------|------|
| §10.3 (Lyddane variables) | Ch 18 | Lyddane modification of Brouwer short-period corrections |
| §10.4 (Newton in Lyddane form) | Ch 20 | Core iteration in osculating element reconstruction |
| §10.6 (State Framework) | Ch 20 | Kepler step in full propagation pipeline |
| §10.5 (convergence bounds) | Ch 9, §9.4 | Parallel treatment; Lyddane bounds reduce to Ch 9 at $a_{yN} = 0$ |

---

### §10.2 The Rotated Form

This section derives the rotated (shifted) Kepler equation in terms of the combined angle $U = E + \omega$, eliminating separate dependence on $\omega$.

- **Definition 10.2.1** (Shifted eccentric anomaly). $U = E + \omega$; the combined angle measured from the ascending node (rather than from periapsis).
- **Theorem 10.2.1** (Rotated Kepler equation). Kepler's equation $M = E - e\sin E$ can be rewritten as $M_L = U - e\sin(U - \omega)$ where $M_L = M + \omega$. — *Proof approach: substitute $E = U - \omega$ into $M = E - e\sin E$; add $\omega$ to both sides to obtain $M_L = U - e\sin(U - \omega)$.*
- **Corollary 10.2.1** (Explicit form). Expanding $\sin(U - \omega)$: $M_L = U - e\cos\omega\sin U + e\sin\omega\cos U$. This is the standard form before introducing Lyddane variables.
- **Remark** (SGP4 usage). In SGP4, the secular rates $\dot{\omega}$ and $\dot{M}$ are accumulated into $M_L = M + \omega$ before solving; the equation is solved directly for $U$ without ever extracting $\omega$ individually at the iteration step.
- **Example 10.2.1** (Rotated form equivalence check). For $a = 7000$ km, $e = 0.15$, $\omega = 60°$, $M = 1.5$ rad: solve both the standard and rotated forms, verify identical $r$ and $\nu$ to machine precision.

---

### §10.3 Lyddane Variables

This section introduces the Lyddane eccentricity variables $(a_{xN}, a_{yN})$ and shows they remove the eccentricity singularity from Kepler's equation.

- **Definition 10.3.1** (Lyddane eccentricity variables). $a_{xN} = e\cos\omega$, $a_{yN} = e\sin\omega$; the Cartesian components of the eccentricity vector in the equatorial plane projection.
- **Theorem 10.3.1** (Lyddane form of Kepler equation). $M_L = U - a_{xN}\sin U + a_{yN}\cos U$. — *Proof approach: substitute $e\cos\omega = a_{xN}$ and $e\sin\omega = a_{yN}$ into the expanded rotated form; apply the angle-addition identity $e\sin(U-\omega) = a_{xN}\sin U - a_{yN}\cos U$.*
- **Corollary 10.3.1** (Singularity removal). At $e = 0$, $a_{xN} = a_{yN} = 0$ and the equation reduces to $M_L = U$, which is nonsingular. The orbit still has a well-defined position; only $\omega$ individually is undefined.
- **Corollary 10.3.2** (Recovery of standard variables). $e_L = \sqrt{a_{xN}^2 + a_{yN}^2}$; $\omega = \text{atan2}(a_{yN}, a_{xN})$ when $e_L > 0$. The atan2 is from Ch 7.
- **Remark** (Relation to Equinoctial elements). The Lyddane pair $(a_{xN}, a_{yN})$ is a subset of the equinoctial element set; the connection is noted as a generalization.
- **Generalization.** Equinoctial elements $(p, f, g, h, k, L)$ where $f = e\cos\omega$, $g = e\sin\omega$, etc., remove both the eccentricity and inclination singularities simultaneously. The Lyddane variables of this chapter are the equinoctial eccentricity components.
- **Example 10.3.1** (Lyddane variables for circular orbit). $e = 0.0001$, $\omega = 123°$: compute $(a_{xN}, a_{yN})$, verify $e_L = \sqrt{a_{xN}^2 + a_{yN}^2}$ recovers $e$, and show the Lyddane equation converges in a single Newton step.

---

### §10.4 Newton-Raphson in Lyddane Form

This section derives the explicit Newton iteration formula in Lyddane variables and verifies its equivalence to the standard form.

- **Definition 10.4.1** (Lyddane Kepler residual). $f_L(U) = U - a_{xN}\sin U + a_{yN}\cos U - M_L$; $f_L'(U) = 1 - a_{xN}\cos U - a_{yN}\sin U$.
- **Proposition 10.4.1** (Second derivative). $f_L''(U) = a_{xN}\sin U - a_{yN}\cos U$. Note that $f_L''(U) = -(f_L(U) - (U - M_L))$, a useful identity for the Halley extension.
- **Definition 10.4.2** (Newton iteration in Lyddane form). $U_{k+1} = U_k - f_L(U_k)/f_L'(U_k)$.
- **Theorem 10.4.1** (Explicit update formula). The Newton step simplifies to: $U_{k+1} = U_k + (M_L - U_k + a_{xN}\sin U_k - a_{yN}\cos U_k)/(1 - a_{xN}\cos U_k - a_{yN}\sin U_k)$. — *Proof approach: substitute the definitions of $f_L$ and $f_L'$ into the Newton formula $U_{k+1} = U_k - f_L/f_L'$; simplify the double negative.*
- **Corollary 10.4.1** (Reduction to standard Newton). When $a_{yN} = 0$ (i.e., $\sin\omega = 0$, equivalently $\omega = 0$ or $\pi$) and $U_k$ is replaced by $E_k + \omega$, Theorem 10.4.1 reduces to the Newton iteration of Ch 9, §9.5.
- **Example 10.4.1** (Convergence table in Lyddane variables). Tabulate $U_k$, $f_L(U_k)$, $|U_{k+1} - U_k|$ for $a_{xN} = 0.45$, $a_{yN} = 0.35$ (corresponding to $e \approx 0.57$, $\omega \approx 37.9°$), $M_L = 3.0$ rad. Show convergence in 4--5 Newton steps.
- **Example 10.4.2** (Near-circular orbit). $e = 0.001$, $\omega = 45°$, $M = 2.0$ rad: compute $(a_{xN}, a_{yN})$, verify that the Lyddane form converges in a single Newton step, and confirm the residual is $O(e^2) \approx 10^{-6}$ after one step.
- *Error Note placeholder:* [P.10.1] the denominator $f_L'(U) = 1 - a_{xN}\cos U - a_{yN}\sin U$ is bounded below by $1 - e$ (since $|a_{xN}\cos U + a_{yN}\sin U| \leq e$); same safe-division bound as Ch 9.

---

### §10.5 Convergence Analysis

This section proves that the Lyddane Newton iteration has identical convergence bounds to the standard Kepler iteration, inheriting all results of Ch 9.

- **Theorem 10.5.1** (Monotonicity). $f_L'(U) = 1 - a_{xN}\cos U - a_{yN}\sin U \geq 1 - \sqrt{a_{xN}^2 + a_{yN}^2} = 1 - e > 0$ for all $U$. Therefore $f_L$ is strictly monotone and has a unique root $U_*$. — *Proof approach: bound $|a_{xN}\cos U + a_{yN}\sin U| \leq \sqrt{a_{xN}^2 + a_{yN}^2}$ by the Cauchy--Schwarz inequality; conclude $f_L' \geq 1 - e > 0$.*
- **Theorem 10.5.2** (Contraction constant). The Newton contraction constant for $f_L$ is $C_L = |f_L''(U_*)| / (2f_L'(U_*)^2)$. Bounding $|f_L''(U)| \leq e$ gives $C_L \leq e / (2(1-e)^2)$, identical to the bound in Ch 9, Theorem 9.5.1. — *Proof approach: compute $f_L''(U) = a_{xN}\sin U - a_{yN}\cos U$; bound $|f_L''| \leq \sqrt{a_{xN}^2 + a_{yN}^2} = e$; substitute into the standard Newton contraction formula.*
- **Corollary 10.5.1** (Starting value). $U_0 = M_L$ satisfies $|f_L(U_0)| = |{-a_{xN}\sin M_L + a_{yN}\cos M_L}| \leq e$; same bound as Ch 9, Proposition 9.9.1.
- **Corollary 10.5.2** (Iteration count bound). The Lyddane Newton iteration starting from $U_0 = M_L$ converges to tolerance $\tau$ in the same number of steps as the standard Newton iteration, because the contraction bound is identical.
- **Proposition 10.5.1** (Halley extension). Halley's method applies to $f_L$ with identical cubic convergence; $f_L'''(U) = -f_L''(U)$, so the Halley denominator evaluates without additional trig calls beyond those already computed for the Newton step.
- *Error Note placeholder:* [P.10.2] precision of $e_L = \sqrt{a_{xN}^2 + a_{yN}^2}$ when $e$ is small; catastrophic cancellation is not possible here since both terms are non-negative.

---

### §10.6 State Framework: Kepler Propagation Step

This section expresses the complete Kepler propagation step — advance mean anomaly by $\Delta M = n\,\Delta t$, solve for the new eccentric anomaly, and compute the new position and velocity — in both the dual quaternion and 7×7 matrix representations.

- **Definition 10.6.1** (Kepler propagation step, dual quaternion form). Given initial state $(\hat{M}_0, \hat{\Omega}_0)$ and time advance $\Delta t$:
  1. Compute $\Delta M = n\,\Delta t$ and update $M_L$.
  2. Solve the Lyddane equation (§10.4) for $U_*$.
  3. Recover $\nu_*$ from $U_*$ and orbital elements via Ch 9, §9.3.
  4. Compute new $\mathbf{r}_F$, $\mathbf{v}_F$ in perifocal frame (Ch 8, §8.6).
  5. Apply perifocal-to-inertial dual quaternion $\hat{M}_{F\to I}$ (Ch 8, §8.7) to obtain updated $(\hat{M}_1, \hat{\Omega}_1)$.
- **Definition 10.6.2** (Kepler propagation step, 7×7 form). The same sequence in the matrix representation: the step matrix $T_{\text{step}}$ maps $(r_0, v_0, 1)^T \to (r_1, v_1, 1)^T$ where the mapping passes through the anomaly solution.
- **Theorem 10.6.1** (Kepler step is near-identity for small $\Delta t$). For $|\Delta t| \ll T$ (orbital period), the configuration dual quaternion $\hat{M}_{\text{step}} = \hat{M}_1 \hat{M}_0^{-1}$ is a near-identity dual quaternion with rotation angle $O(\Delta M)$ and translation $O(v\,\Delta t)$. The 7x7 matrix $T_{\text{step}} = T_1 T_0^{-1}$ is correspondingly near-identity. — *Proof approach: expand the dual quaternion composition for small rotation angle $\Delta\nu$; show the real part is $1 + O((\Delta\nu)^2)$ and the dual part is $O(\Delta\nu)$.*
- **Corollary 10.6.1** (Error propagation for Kepler step). By Ch 2, Principle 2, the precision error in the updated state is bounded by the precision error in the anomaly solution ($O(\tau/f'(U_*))$ from the iteration) composed with the perifocal-frame position errors.
- **Example 10.6.1** (One-step propagation). For an ISS-class orbit ($a = 6778$ km, $e = 0.0007$, $i = 51.6°$, $\omega = 30°$, $\Omega = 200°$, $M_0 = 0°$): advance by one orbital period ($T \approx 5553$ s), verify return to initial elements within rounding ($< 10^{-12}$ rad in angles). Demonstrate in both dual quaternion and 7x7 representations; verify numerical agreement.
- **Remark** (Velocity dual quaternion update). The Kepler step changes the velocity direction (tangent to the orbit) according to the new $\nu$; the velocity dual quaternion $\hat{\Omega}_1$ encodes the new angular velocity and translational velocity of the satellite in the inertial frame.
- **Remark** (Connection to Ch 20). In the full SGP4 pipeline, the Kepler step is embedded inside the mean-to-osculating correction sequence of Ch 20. The notation $U_*$ maps to the SGP4 variable `U` (equation of center argument) in the Hoots–Roehrich 1980 notation.
- *Error Note placeholder:* [A.10.1] the Kepler step here is the unperturbed two-body propagation; it neglects $J_2$, drag, and other perturbations. The error accumulated per orbital period is $O(J_2)$ in relative position, corrected by the Brouwer secular terms in Ch 16–19.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Thm 1.4.1 (convergence error) | Error tracking; convergence error classification |
| Ch 2 | Thm 2.3 composition, Thm 2.8 matrix, Thm 2.10 parallel-row | Dual SU(2) composition; 7x7 matrix; parallel-row architecture; Principles 1 and 2 |
| Ch 3 | Thm 3.4.1 (tolerance parameter) | Tolerance parameter; matched-pair principle; SGP4-standard Lyddane iteration reproduces Hoots-Roehrich at tau_standard |
| Ch 7 | atan2 quadrant resolution | atan2 for omega recovery from (a_xN, a_yN) |
| Ch 8 | Perifocal frame and State Framework | Perifocal frame, position/velocity formulae, State Framework template |
| Ch 9 | Convergence theory and iteration methods | All convergence theory; contraction bounds; Halley extension; starting values carry over to Lyddane form |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 18 | Lyddane short-period corrections | Lyddane modification of Brouwer short-period corrections uses (a_xN, a_yN) |
| Ch 20 | Osculating element reconstruction | The modified Kepler equation is the central computational step |
| Appendix B | modified_kepler.h implementation | Implementation specifications |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.10.1] | P | §10.4 | Denominator bound for Lyddane Newton iteration; $f_L'(U) \geq 1 - e$; division amplification at most $1/(1-e)$, identical to standard Kepler case |
| [P.10.2] | P | §10.5 | Recovery of eccentricity from Lyddane variables; relative error $O(\epsilon_{\rm mach}/e_L)$ for small $e$; acceptable since $a_{xN}$, $a_{yN}$ remain well-defined |
| [A.10.1] | A | §10.6 | Unperturbed Kepler step neglects all perturbations; position error ~10 km/day at 400 km; perturbations reinstated in Parts IV–VI |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 6 |
| Theorems | 6 |
| Lemmas | 0 |
| Corollaries | 7 |
| Propositions | 2 |
| Examples | 5 |
| Error Notes | 3 |
| Equations | ~20 |
| Sections | 6 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §10.1 Introduction | Draft | |
| §10.2 The Rotated Form | Draft | |
| §10.3 Lyddane Variables | Draft | |
| §10.4 Newton-Raphson in Lyddane Form | Draft | |
| §10.5 Convergence Analysis | Draft | |
| §10.6 State Framework: Kepler Propagation Step | Draft | |
