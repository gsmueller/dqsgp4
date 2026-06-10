# Draft Plan: Chapter 30 — Coordinate Transformations

**Part VIII: Reference Frames and Time** | Implementation file: `state_from_elements.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\hat{U}_b, \hat{V}_b$ | Unit vectors in TEME frame for radial and transverse directions | §30.2 |
| $\mathbf{r}_b$ | Position vector expressed in TEME basis $b$ | §30.2 |
| $\mathbf{v}_b$ | Velocity vector expressed in TEME basis $b$ | §30.2 |
| $\theta_{\mathrm{GMST}}$ | Greenwich Mean Sidereal Time (Ch 29, Def. 29.3.1) | Ch 29 |
| $R_3(\theta)$ | Rotation about $z$-axis by angle $\theta$ (active or passive, convention stated) | §30.1 |
| $R_1(\theta)$ | Rotation about $x$-axis by angle $\theta$ | §30.1 |
| TEME | True Equator, Mean Equinox reference frame | §30.1 |
| PEF | Pseudo Earth Fixed frame | §30.3 |
| ECEF | Earth-Centered Earth-Fixed frame (WGS84) | §30.3 |
| $u$ | Argument of latitude: $u = \omega + \nu$ | Ch 20 |
| $r_k, \dot{r}_k, (r\dot{f})_k$ | Radial distance, radial velocity, transverse velocity from osculating solution | Ch 20 |
| $M_b$ | Configuration dual quaternion in basis $b$ | Ch 2 |
| $\hat{\Omega}_b$ | Velocity dual quaternion in basis $b$ | Ch 2 |

---

## Objectives

1. Define the TEME frame precisely and establish its limitations as a practical (not rigorous) approximation.
2. Derive unit vectors $\hat{U}_b, \hat{V}_b$ from the orbital element sequence $R_3(-\Omega)R_1(-i)R_3(-u)$ applied to $\hat{x}_b$.
3. Construct position and velocity vectors in TEME.
4. Derive the TEME → PEF rotation via $R_3(\theta_{\mathrm{GMST}})$; state the PEF → ECEF correction (polar motion).
5. State Framework: all transforms in dual-quaternion form acting on $(\hat{M}, \hat{\Omega}_b)$ and in 7×7 matrix form acting on $(r, v, 1)^T$.

## Section Structure

### §30.1 Reference Frame Definitions

This section defines the TEME, PEF, and ECEF reference frames and establishes the passive rotation convention used throughout the chapter.

Stub:

1. **Definition 30.1.1** (TEME — True Equator, Mean Equinox)**.** *The TEME frame has its $z$-axis aligned with Earth's instantaneous rotation axis (true equator of date) and its $x$-axis along the mean vernal equinox of date. It is the native output frame of the SGP4 propagator by convention.*
2. **Definition 30.1.2** (PEF — Pseudo Earth Fixed)**.** *PEF shares the TEME $z$-axis; its $x$-axis is rotated from the mean equinox to the Greenwich meridian by the passive rotation $R_3(\theta_{\mathrm{GMST}})$.*
3. **Definition 30.1.3** (ECEF — Earth-Centered Earth-Fixed)**.** *The WGS84 ECEF frame differs from PEF by a small polar motion rotation $W$ (IERS $x_p, y_p \lesssim 0.3''$) applied to align the $z$-axis with the Conventional International Origin.*
4. **Proposition 30.1.1** (passive rotation convention)**.** *All rotation matrices $R_k(\theta)$ in this chapter are passive (alias) rotations: they re-express the same vector in a rotated basis, so $\mathbf{r}_{\mathrm{new}} = R(\theta)\,\mathbf{r}_{\mathrm{old}}$ transforms coordinates, not the physical vector.*
   — *Proof approach: orthogonality verification by direct matrix multiplication — verify $R_k(\theta)^T R_k(\theta) = I$ for each elementary rotation.*

[A.30.1] TEME is not a standard IAU frame; it is the output frame of SGP4 by convention. Users requiring ECEF must apply the TEME → ECEF chain; users requiring J2000 ECI must apply an additional nutation rotation.

**Example 30.1.1** (frame chain summary)**.** *The complete chain is: perifocal $\xrightarrow{R_3(-\Omega)R_1(-i)R_3(-u)}$ TEME $\xrightarrow{R_3(\theta_{\mathrm{GMST}})}$ PEF $\xrightarrow{W}$ ECEF. Source: [Vallado 2013, §3.7]; the chain is also documented in [Hoots and Roehrich 1980, §6] with the TEME convention implicit.*

### §30.2 Unit Vectors and State in TEME

This section derives the explicit Cartesian components of the TEME radial and transverse unit vectors from the three-rotation sequence and constructs the satellite state vector.

Stub: Starting from $\hat{x}_b = (1,0,0)^T_b$ in the perifocal frame: apply $R_3(-\Omega)$ (rotate about $z$ by $-\Omega$ to align node), then $R_1(-i)$ (tilt by inclination), then $R_3(-u)$ (rotate by argument of latitude). Result: $\hat{U}_b$ along the satellite radial, $\hat{V}_b = \hat{z}_b \times \hat{U}_b / \|\cdot\|$ along the in-track transverse. Position $\mathbf{r}_b = r_k\,\hat{U}_b$, velocity $\mathbf{v}_b = \dot{r}_k\,\hat{U}_b + (r\dot{f})_k\,\hat{V}_b$.

1. **Theorem 30.2.1** (TEME unit vector components)**.** *The components of $\hat{U}_b$ and $\hat{V}_b$ in the TEME frame are:*

$$\hat{U}_b = \begin{pmatrix} \cos\Omega\cos u - \sin\Omega\sin u\cos i \\ \sin\Omega\cos u + \cos\Omega\sin u\cos i \\ \sin u\sin i \end{pmatrix}, \quad \hat{V}_b = \begin{pmatrix} -\cos\Omega\sin u - \sin\Omega\cos u\cos i \\ -\sin\Omega\sin u + \cos\Omega\cos u\cos i \\ \cos u\sin i \end{pmatrix} \tag{30.1}$$

   — *Proof approach: rotation matrix composition ($R_z \cdot R_x \cdot R_z$ Euler angles) — expand $R_3(-\Omega)R_1(-i)R_3(-u)\hat{x}_b$ by explicit $3\times3$ matrix multiplication; collect trigonometric terms. Orthogonality $\hat{U}_b \cdot \hat{V}_b = 0$ and unit norm $\|\hat{U}_b\| = \|\hat{V}_b\| = 1$ verified by direct computation.*

[P.30.1] Near $i=0$ or $e=0$, $\Omega$ and $\omega$ are individually singular; $u = \omega + \nu$ remains well-defined (Lyddane variables protect $e\cos\omega, e\sin\omega$ from singularity, but $u$ itself does not depend on them separately).

**Example 30.2.1** (TEME state vector at J2000.0)**.** *For a circular orbit with $a = 7000$ km, $e = 0.001$, $i = 51.6°$, $\Omega = 0°$, $\omega = 0°$, $\nu = 0°$ at epoch J2000.0: $u = 0°$, $\hat{U}_b = (1, 0, 0)^T$, $\mathbf{r}_b \approx (7000, 0, 0)^T$ km. Source: trivial case of Eq. (30.1); cross-check against [Vallado 2013, Example 2-5] for a non-trivial case.*

### §30.3 TEME to PEF to ECEF

This section derives the TEME → PEF rotation and velocity Coriolis correction, states the PEF → ECEF polar motion correction, and bounds the error from neglecting polar motion.

Stub:

1. **Theorem 30.3.1** (TEME to PEF position)**.** *The position vector in PEF is $\mathbf{r}_{\mathrm{PEF}} = R_3(\theta_{\mathrm{GMST}})\,\mathbf{r}_b$, where $R_3(\theta_{\mathrm{GMST}})$ is the passive rotation about $z$ by $\theta_{\mathrm{GMST}}$ (Ch 29, Def. 29.3.1).*
   — *Proof approach: rotation matrix composition — $R_3(\theta_{\mathrm{GMST}})$ is the unique $SO(3)$ element mapping the mean equinox direction to the Greenwich meridian; orthogonality verified by direct matrix multiplication.*

2. **Theorem 30.3.2** (TEME to PEF velocity)**.** *The velocity vector in PEF includes a Coriolis term from Earth rotation: $\mathbf{v}_{\mathrm{PEF}} = R_3(\theta_{\mathrm{GMST}})\,\mathbf{v}_b - \omega_E\,\hat{z}_{\mathrm{PEF}} \times \mathbf{r}_{\mathrm{PEF}}$.*
   — *Proof approach: rotation matrix composition — differentiate $\mathbf{r}_{\mathrm{PEF}}(t) = R_3(\theta_{\mathrm{GMST}}(t))\,\mathbf{r}_b(t)$ with respect to time; $\dot{R}_3 = \omega_E\,\hat{z} \times R_3$ yields the transport term.*

3. **Lemma 30.3.1** (polar motion bound)**.** *The PEF → ECEF polar motion rotation $W$ satisfies $\|W - I\| < 1.5\times10^{-6}$ for all tabulated IERS polar motion values $x_p, y_p \lesssim 0.3''$.*
   — *Proof approach: orthogonality verification by direct matrix multiplication — bound $\|W - I\|$ in terms of $\max(|x_p|, |y_p|)$ using the small-angle approximation $\sin\theta \approx \theta$.*

[A.30.2] SGP4 convention ignores polar motion; error $< 15$ m in position.

**Example 30.3.1** (TEME to ECEF at J2000.0)**.** *At J2000.0, $\theta_{\mathrm{GMST}} = 280.460\,618\,37°$ (Ch 29, Example 29.3.1). For $\mathbf{r}_b = (7000, 0, 0)^T$ km: $\mathbf{r}_{\mathrm{PEF}} = R_3(280.460\,618\,37°)\,(7000, 0, 0)^T = (7000\cos(280.46°),\; 7000\sin(280.46°),\; 0)^T \approx (1274.9,\; -6882.3,\; 0)^T$ km. Source: trigonometric evaluation of Eq. (30.1) with GMST from [Aoki et al. 1982]; cross-check against [Vallado 2013, Example 3-5].*

### §30.4 State Framework: Both Representations

This section expresses all coordinate transformation steps in dual-quaternion form and in 7×7 matrix form, following the framework of Ch 2.

Stub: **Dual-quaternion form.** The sequence $R_3(-\Omega)R_1(-i)R_3(-u)$ corresponds to a composition of three configuration dual quaternions (Ch 2, §2.4):

$$\hat{M}_b^{\mathrm{TEME}} = \hat{M}_{R_3(-\Omega)} \cdot \hat{M}_{R_1(-i)} \cdot \hat{M}_{R_3(-u)} \tag{30.2}$$

Each elementary rotation $R_k(\alpha)$ maps to a unit dual quaternion with real part $\cos(\alpha/2)$ and imaginary part determined by the rotation axis. Velocity dual quaternion $\hat{\Omega}_b$ transforms under the adjoint action (Ch 2, §2.7). The TEME → PEF step adds the $\omega_E$-coupling term to the velocity dual part.

1. **Proposition 30.4.1** (dual-quaternion TEME state)**.** *The configuration dual quaternion $\hat{M}_b^{\mathrm{TEME}}$ encodes the full TEME orientation; the TEME position $\mathbf{r}_b$ is recovered as the dual part of $\hat{M}_b^{\mathrm{TEME}}\,\hat{e}_r\,(\hat{M}_b^{\mathrm{TEME}})^*$ where $\hat{e}_r$ is the radial unit pure-quaternion.*
   — *Proof approach: dual-quaternion conjugation action — verify by expanding the sandwich product $\hat{M}\,\hat{p}\,\hat{M}^*$ for a pure-translation dual quaternion $\hat{p}$; the result is the Ch 2 state framework (Ch 2, §2.5).*

**7×7 matrix form.** Each $R_k$ becomes a $4\times4$ block (rotation sub-block, affine column = 0, homogeneous row); the TEME → PEF velocity coupling populates the off-diagonal transport block. The full TEME → ECEF chain is a product of 7×7 matrices. Error propagation via Principle 1 (composition) and Principle 2 (force/perturbation) of Ch 2. Cross-reference Ch 2, §2.8 for the 7×7 assembly.

**Example 30.4.1** (7×7 TEME matrix structure)**.** *For a pure rotation $R_3(\alpha)$ acting on the position-velocity state $(r, v, 1)^T$: the 7×7 matrix has $R_3(\alpha)$ in the $(1{:}3, 1{:}3)$ block, $R_3(\alpha)$ in the $(4{:}6, 4{:}6)$ block, the transport term $-\omega_E[\hat{z}_{\times}]R_3(\alpha)$ in the $(4{:}6, 1{:}3)$ block, and 1 in position $(7,7)$. Source: Ch 2, §2.8 assembly rules; [Vallado 2013, §3.7] for the physics.*

### §30.5 Implementation Notes

This section specifies computational strategies for efficient TEME state vector construction, including Lyddane variable protection and per-step constant management.

Stub: All six Cartesian components of $\hat{U}_b, \hat{V}_b$ computed in one pass; no redundant trig evaluations. The $\omega_E$ factor applied to the PEF velocity correction is the constant from Ch 29, §29.5 — precomputed at initialization. Near-circular or near-equatorial orbits use Lyddane variable protection (Ch 10, §ch10-modified-kepler; Ch 18, §ch18-short-period-corrections) to avoid division by small $e$ or $\sin i$.

1. **Proposition 30.5.1** (one-pass unit vector computation)**.** *The six components of $\hat{U}_b$ and $\hat{V}_b$ from Eq. (30.1) require exactly six trigonometric products, two of which ($\cos u$, $\sin u$, $\cos i$, $\sin i$, $\cos\Omega$, $\sin\Omega$) can be evaluated once and reused. No intermediate rotation matrices need be formed explicitly.*

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 2, Theorem 2.8.1 | §30.4 | 7x7 matrix assembly and adjoint action |
| Ch 2, dual quaternion composition | §30.4 | Configuration and velocity dual quaternion transforms |
| Ch 9–10, Kepler equation outputs | §30.2 | $r_k, \dot{r}_k, (r\dot{f})_k$ from osculating solution |
| Ch 18, Lyddane variables | §30.5 | Protection for near-circular and near-equatorial orbits |
| Ch 20, osculating elements | §30.2 | Argument of latitude $u$ and radial distance $r_k$ |
| Ch 29, sidereal time | §30.3 | $\theta_{\mathrm{GMST}}$ and $\omega_E$ for TEME to PEF rotation |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 34, near-space pipeline | §30.2, §30.3 | TEME state vector and coordinate chain |
| Ch 35, deep-space pipeline | §30.2, §30.3 | TEME state vector and coordinate chain |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.30.1] | P | §30.2 | Near $i=0$: $\Omega, \omega$ singular individually; $u$ well-defined |
| [A.30.1] | A | §30.1 | TEME is a non-standard approximation to ECI; nutation ignored |
| [A.30.2] | A | §30.3 | Polar motion neglected in SGP4 convention; $< 15$ m position error |
| [A.30.3] | A | §30.4 | Lieske et al. (1977) IAU 1976 precession constant is approximately 3 mas/yr too high (Williams et al. 1991); fully superseded by IAU 2000 model (Capitaine et al. 2003, A&A 412, 567); for TEME-to-ICRF conversions better than 0.1 arcsecond accuracy, use IERS Conventions 2003 or later |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 4 |
| Theorems | 3 |
| Lemmas | 1 |
| Corollaries | 0 |
| Propositions | 3 |
| Examples | 5 |
| Error Notes | 4 |
| Equations | ~14 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §30.1 | Draft | Reference frame definitions |
| §30.2 | Draft | Unit vectors and state in TEME |
| §30.3 | Draft | TEME to PEF to ECEF |
| §30.4 | Draft | State Framework: both representations |
| §30.5 | Draft | Implementation notes |
