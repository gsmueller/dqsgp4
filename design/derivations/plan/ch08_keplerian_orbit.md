# Draft Plan: Chapter 8 — The Keplerian Orbit

**Part II: The Two-Body Problem**

**Phase:** Draft

---

## Objectives

1. Derive the inverse-square central force law from Newton's law of gravitation and prove that the resulting orbit is a conic section (Binet equation approach).
2. Establish the two conservation laws — angular momentum and energy — as first integrals of the equations of motion, and derive the vis-viva equation.
3. Derive the orbit equation $r = p/(1 + e\cos\nu)$ from the Binet equation.
4. Define the six Keplerian orbital elements $(a, e, i, \Omega, \omega, M)$ and establish their geometric meaning.
5. Derive position and velocity in the perifocal frame, including the direction cosine matrix relating the perifocal and geocentric equatorial frames.
6. Express the perifocal-to-inertial transform in both the dual quaternion representation $(\hat{M}, \hat{\Omega})$ and the equivalent 7×7 matrix acting on $(r, v, 1)^T$.
7. Establish error tier classifications for the Keplerian elements and the derived position/velocity.

---

## Notation Table

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mu$ | Geocentric gravitational constant $GM$ | §8.2 |
| $\mathbf{r}_I$ | Position vector in inertial frame (basis subscript I) | §8.2 |
| $\mathbf{v}_I$ | Velocity vector in inertial frame | §8.2 |
| $r$ | Scalar distance: $r = \|\mathbf{r}_I\|$ | §8.2 |
| $\mathbf{h}_I$ | Specific angular momentum vector: $\mathbf{h}_I = \mathbf{r}_I \times \mathbf{v}_I$ | §8.3 |
| $h$ | Magnitude of specific angular momentum: $h = \|\mathbf{h}_I\|$ | §8.3 |
| $\mathcal{E}$ | Specific mechanical energy: $\mathcal{E} = v^2/2 - \mu/r$ | §8.3 |
| $\mathbf{e}_I$ | Eccentricity vector (Laplace–Runge–Lenz), in inertial frame | §8.3 |
| $p$ | Semi-latus rectum: $p = h^2/\mu$ | §8.4 |
| $e$ | Eccentricity (scalar): $e = \|\mathbf{e}_I\|$ | §8.4 |
| $\nu$ | True anomaly | §8.4 |
| $a$ | Semi-major axis: $a = -\mu/(2\mathcal{E})$ | §8.5 |
| $i$ | Inclination | §8.5 |
| $\Omega$ | Right ascension of the ascending node (RAAN) | §8.5 |
| $\omega$ | Argument of perigee | §8.5 |
| $M$ | Mean anomaly | §8.5 |
| $n$ | Mean motion: $n = \sqrt{\mu/a^3}$ | §8.5 |
| $T$ | Orbital period: $T = 2\pi/n$ | §8.5 |
| $\hat{P}$ | Perifocal unit vector in the direction of periapsis | §8.6 |
| $\hat{Q}$ | Perifocal unit vector 90° ahead of periapsis in the orbit plane | §8.6 |
| $\hat{W}$ | Orbit normal unit vector: $\hat{W} = \hat{P} \times \hat{Q}$ | §8.6 |
| $\mathbf{r}_F$ | Position in perifocal frame (basis subscript F) | §8.6 |
| $\mathbf{v}_F$ | Velocity in perifocal frame | §8.6 |
| $R(\theta, \hat{k})$ | Rotation by angle $\theta$ about axis $\hat{k}$ (dual quaternion form) | §8.7 |
| $\hat{M}_{F\to I}$ | Configuration dual quaternion for perifocal-to-inertial transform | §8.7 |
| $T_{F\to I}$ | 7×7 matrix for perifocal-to-inertial transform | §8.7 |

---

## Section Structure

### §8.1 Introduction

This section motivates the two-body problem as the zeroth-order model and maps forward dependencies to subsequent chapters.

Narrative placing the two-body problem as the zeroth-order model for all orbit propagation. Motivation: the exact solution is a conic section; all subsequent chapters add perturbations to this baseline.

Forward reference table:

| Section | Feeds | Role |
|---------|-------|------|
| §8.5 (elements) | Ch 9 | Mean anomaly defines the Kepler equation argument |
| §8.5 (elements) | Ch 10 | Modified Kepler uses $\omega$ and eccentricity vector |
| §8.6 (perifocal state) | Ch 11 | Delaunay variables derived from orbital elements |
| §8.7 (dual quaternion transform) | Ch 2 | Application of §2.3 composition to orbit geometry |
| §8.7 (7×7 form) | Ch 20 | Osculating element reconstruction uses inverse transform |

---

### §8.2 Newton's Law and the Equations of Motion

This section derives the inverse-square equations of motion and proves the orbit is a conic section via the Binet equation.

- **Definition 8.2.1** (Two-body problem). Equations of motion $\ddot{\mathbf{r}}_I = -(\mu/r^3)\mathbf{r}_I$; reduction from two-body to one-body; reduced mass; identification of $\mu = G(M_1 + M_2) \approx GM_\oplus$.
- **Theorem 8.2.1** (Equations of motion in polar form). Radial and transverse components of $\ddot{\mathbf{r}}_I$ in the orbit plane. — *Proof approach: decompose acceleration into radial/transverse components using rotating frame identities.*
- **Definition 8.2.2** (Binet variable). $u = 1/r$; the Binet substitution $\dot{r} = -h\,du/d\nu$.
- **Theorem 8.2.2** (Binet equation). $d^2u/d\nu^2 + u = \mu/h^2$. — *Proof approach: substitute Binet variable into radial equation; use $\dot{r} = -h\,du/d\nu$ and chain rule through $\nu$.*
- **Theorem 8.2.3** (General solution). $u = \mu/h^2 + C\cos(\nu - \nu_0)$; geometric interpretation as a conic section. — *Proof approach: solve the linear second-order ODE with constant coefficients; identify the homogeneous and particular solutions.*
- **Example 8.2.1** (Circular orbit check). $e = 0$: verify $r = p = h^2/\mu$ is constant; $v^2 = \mu/r$. Values: $\mu = 398600.4418$ km$^3$/s$^2$, $r = 6778$ km (ISS-class); compute $v_{\text{circ}}$ and verify orbit equation reduces to constant radius.
- *Error Note placeholder:* [A.8.1] point mass approximation; [A.8.2] neglect of relativistic precession.

---

### §8.3 Conservation Laws

This section establishes the three first integrals (angular momentum vector, eccentricity vector, energy) and derives the vis-viva equation.

- **Theorem 8.3.1** (Conservation of angular momentum). $d\mathbf{h}_I/dt = \mathbf{0}$; the orbit is planar. — *Proof approach: compute $d\mathbf{h}_I/dt = \mathbf{r}_I \times \ddot{\mathbf{r}}_I$; central force gives zero cross product.*
- **Corollary 8.3.1** (Areal velocity). $dA/dt = h/2$ (Kepler's second law). — *Proof approach: express swept area in polar coordinates; relate $d\nu/dt$ to $h/r^2$.*
- **Theorem 8.3.2** (Eccentricity vector). $\mathbf{e}_I = \mathbf{v}_I \times \mathbf{h}_I/\mu - \hat{r}_I$ is a constant of motion. — *Proof approach: differentiate $\mathbf{e}_I$ with respect to time; substitute $\ddot{\mathbf{r}} = -\mu\mathbf{r}/r^3$ and use vector triple product identity.*
- **Theorem 8.3.3** (Conservation of energy / vis-viva equation). $\mathcal{E} = v^2/2 - \mu/r = -\mu/(2a)$ is constant; derive $v^2 = \mu(2/r - 1/a)$. — *Proof approach: compute $d\mathcal{E}/dt = \dot{\mathbf{v}}\cdot\mathbf{v} + \mu\dot{r}/r^2$; substitute equations of motion to show cancellation.*
- **Corollary 8.3.2** (Energy–semi-major axis relation). $a = -\mu/(2\mathcal{E})$; sign classification for ellipse/parabola/hyperbola.
- **Example 8.3.1** (ISS approximate orbit). Numerical vis-viva: $a \approx 6770$ km, $e \approx 0.0007$, $\mu = 398600.4418$ km$^3$/s$^2$; compute $v_{\text{circ}} \approx 7.67$ km/s, $\mathcal{E}$, and $v$ at perigee and apogee. Source: representative ISS-class elements from Vallado (2013), Table 2-1.
- *Error Note placeholder:* [A.8.3] two-body energy neglects $J_2$, drag, third-body; magnitude estimate.

---

### §8.4 The Orbit Equation

This section derives the orbit equation in polar form and classifies conic sections by eccentricity.

- **Theorem 8.4.1** (Orbit equation). $r = p/(1 + e\cos\nu)$ where $p = h^2/\mu$ and $e = \|\mathbf{e}_I\|$. — *Proof approach: take the dot product $\mathbf{e}_I \cdot \mathbf{r}_I$; expand using the definition of $\mathbf{e}_I$ and the identity $\mathbf{r}\cdot(\mathbf{v}\times\mathbf{h}) = \mathbf{h}\cdot(\mathbf{r}\times\mathbf{v}) = h^2$; solve for $r$.*
- **Definition 8.4.1** (True anomaly $\nu$). Angle from the eccentricity vector to the position vector, in the direction of motion.
- **Corollary 8.4.1** (Conic classification). $e < 1$: ellipse; $e = 1$: parabola; $e > 1$: hyperbola.
- **Corollary 8.4.2** (Apoapsis and periapsis). $r_p = p/(1+e) = a(1-e)$; $r_a = p/(1-e) = a(1+e)$ (elliptic case).
- **Corollary 8.4.3** (Semi-latus rectum–semi-major axis relation). $p = a(1-e^2)$.
- **Example 8.4.1** (Orbit parameters from energy and momentum). Given $\mathcal{E} = -29.5$ km$^2$/s$^2$ and $h = 52000$ km$^2$/s (representative LEO), compute $a$, $e$, $p$ and verify consistency with the orbit equation at $\nu = 0$ and $\nu = \pi$.

---

### §8.5 Keplerian Orbital Elements

This section defines the six Keplerian elements, derives Kepler's third law, and establishes the state vector conversion algorithms.

- **Definition 8.5.1** (The six Keplerian elements). $(a, e, i, \Omega, \omega, M)$: semi-major axis, eccentricity, inclination, RAAN, argument of perigee, mean anomaly. Geometric definition of each.
- **Definition 8.5.2** (Mean motion). $n = \sqrt{\mu/a^3}$ (Kepler's third law).
- **Definition 8.5.3** (Mean anomaly and epoch). $M = M_0 + n(t - t_0)$; the linear propagation valid in the unperturbed two-body problem.
- **Theorem 8.5.1** (Kepler's third law). $T^2 = (4\pi^2/\mu)a^3$. — *Proof approach: integrate the areal velocity $dA/dt = h/2$ over one full orbit; equate total area $\pi ab$ to $hT/2$; substitute $b = a\sqrt{1-e^2}$ and $h = \sqrt{\mu p}$.*
- **Definition 8.5.4** (Ascending node and RAAN). The line of nodes as $\hat{n} = \hat{k}_I \times \hat{h}_I / \|\hat{k}_I \times \hat{h}_I\|$; $\Omega$ measured in the equatorial plane.
- **Remark** (Singularities). $e = 0$: $\omega$ undefined; $i = 0$ or $\pi$: $\Omega$ undefined. Motivation for Chapters 9–10 modifications.
- **Example 8.5.1** (State vector to elements). Given $\mathbf{r}_I$ and $\mathbf{v}_I$, compute all six elements. Values: $\mathbf{r}_I = (6524.834, 6862.875, 6448.296)$ km, $\mathbf{v}_I = (4.901, 5.533, -1.976)$ km/s (Vallado 2013, Example 2-6 class). Compute $a, e, i, \Omega, \omega, M$ step-by-step.
- **Example 8.5.2** (Elements to state vector). Inverse of Example 8.5.1 via perifocal frame (§8.6). Start from the six elements computed in Example 8.5.1, reconstruct $\mathbf{r}_I$ and $\mathbf{v}_I$, verify round-trip agreement to machine precision.
- *Error Note placeholder:* [M.8.1] TLE element measurement errors; [A.8.4] osculating vs. mean elements distinction deferred to Ch 12.

---

### §8.6 Position and Velocity in the Perifocal Frame

This section derives explicit position and velocity expressions in the perifocal frame and constructs the perifocal-to-inertial DCM.

- **Definition 8.6.1** (Perifocal frame). Orthonormal basis $(\hat{P}, \hat{Q}, \hat{W})$ fixed to the orbit plane with $\hat{P}$ toward periapsis.
- **Theorem 8.6.1** (Position in perifocal frame). $\mathbf{r}_F = r\cos\nu\,\hat{P} + r\sin\nu\,\hat{Q}$ with $r = p/(1 + e\cos\nu)$. — *Proof approach: direct projection of the orbit equation onto the perifocal basis vectors.*
- **Theorem 8.6.2** (Velocity in perifocal frame). $\mathbf{v}_F = (\mu/h)(-\sin\nu\,\hat{P} + (e + \cos\nu)\,\hat{Q})$. — *Proof approach: differentiate $\mathbf{r}_F$ with respect to time; substitute $\dot{\nu} = h/r^2$ and simplify using the orbit equation.*
- **Corollary 8.6.1** (Speed from elements). $v^2 = (\mu/h)^2(1 + e^2 + 2e\cos\nu)$; consistency check with vis-viva. — *Proof approach: compute $|\mathbf{v}_F|^2$ from Theorem 8.6.2; substitute $p = h^2/\mu$ and verify equivalence with vis-viva.*
- **Definition 8.6.2** (Perifocal-to-inertial DCM). $R_{F\to I}$: the $3\times 3$ direction cosine matrix as three successive rotations $R_3(-\Omega)\,R_1(-i)\,R_3(-\omega)$ where $R_j$ denotes rotation about axis $j$.
- **Theorem 8.6.3** (DCM explicit form). Expansion of $R_{F\to I}$ in terms of $\cos\Omega$, $\sin\Omega$, $\cos i$, $\sin i$, $\cos\omega$, $\sin\omega$. — *Proof approach: multiply the three rotation matrices $R_3(-\Omega)\,R_1(-i)\,R_3(-\omega)$ and expand each entry.*
- **Example 8.6.1** (Perifocal state for ISS-class orbit). Numerical evaluation of $\mathbf{r}_F$, $\mathbf{v}_F$, and $R_{F\to I}$ for elements from Example 8.5.1. Compute perifocal position and velocity at $\nu = 60°$, then apply $R_{F\to I}$ to recover $\mathbf{r}_I$, $\mathbf{v}_I$.
- *Error Note placeholder:* [P.8.1] finite-precision loss in $R_{F\to I}$ near $i \approx 0$ (RAAN singularity); magnitude estimate.

---

### §8.7 State Framework: Perifocal-to-Inertial Transform

The two-body state $(\mathbf{r}_I, \mathbf{v}_I)$ lives in the inertial frame. The perifocal frame provides a natural coordinate system for the instantaneous orbit. This section expresses the frame transform in both the dual quaternion representation and the 7×7 matrix form, establishing the template used by all subsequent orbit transform steps.

- **Definition 8.7.1** (Configuration dual quaternion for perifocal-to-inertial). $\hat{M}_{F\to I}$: the dual SU(2) element encoding the three successive rotations $R_3(-\Omega)$, $R_1(-i)$, $R_3(-\omega)$ as a single composition. Composition law from Ch 2, §2.3.
- **Theorem 8.7.1** (Dual quaternion perifocal-to-inertial). Given the configuration dual quaternion $\hat{M}_{F\to I}$, the inertial position vector is $\mathbf{r}_I = \hat{M}_{F\to I}\,\mathbf{r}_F\,\hat{M}_{F\to I}^\dagger$ (adjoint action on pure dual quaternion); inertial velocity similarly via velocity dual quaternion $\hat{\Omega}$. — *Proof approach: compose three dual quaternion rotations from Ch 2, §2.3; apply the adjoint action formula and expand.*
- **Theorem 8.7.2** (7x7 matrix perifocal-to-inertial). The transform acts as $\begin{pmatrix}\mathbf{r}_I \\ \mathbf{v}_I \\ 1\end{pmatrix} = T_{F\to I}\begin{pmatrix}\mathbf{r}_F \\ \mathbf{v}_F \\ 1\end{pmatrix}$ where $T_{F\to I}$ has the explicit block structure from Ch 2, §2.8 applied to the rotation-only case (no translation in the perifocal origin). — *Proof approach: specialize the general 7x7 block structure of Ch 2 to a pure rotation (zero translation, zero angular velocity); write out the block entries explicitly.*
- **Corollary 8.7.1** (Equivalence). The dual quaternion action and the 7x7 matrix produce identical numerical results on the 7-vector state. — *Proof approach: expand both representations component-by-component and verify algebraic identity of all 7 output entries.*
- **Example 8.7.1** (Numerical perifocal-to-inertial). From Example 8.6.1, apply both representations (dual quaternion and 7x7 matrix) to the perifocal state and verify agreement to machine precision ($< 10^{-15}$ relative error in each component).
- **Remark** (Velocity dual quaternion). In the orbit transform, the velocity dual quaternion $\hat{\Omega}$ encodes the angular velocity of the perifocal frame relative to inertial; this is zero for the Keplerian orbit (the perifocal frame is inertial relative to the orbit). The $\hat{\Omega}$ contribution enters in Ch 16 when secular precession of $\Omega$ and $\omega$ is added.
- *Error Note placeholder:* [P.8.2] quaternion normalization drift over many compositions; mitigation by renormalization.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Thm 1.2.4 (TrackedValue) | Error categories, tier classification |
| Ch 2 | Thm 2.3 composition, Thm 2.8 matrix, Thm 2.10 parallel-row | Dual SU(2) composition, 7x7 matrix assembly, parallel-row architecture |
| Ch 3 | Matched pair principle, mu selection | Matched pair principle; mu value selection (SGP4 standard vs. physical) |
| Ch 7 | Angle arithmetic definitions | Angle arithmetic for Omega, omega, nu, M in [0, 2pi) |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 9 | Kepler equation inputs | Mean anomaly M and semi-major axis a are the inputs to Kepler's equation |
| Ch 10 | Modified Kepler inputs | Eccentricity e and argument of perigee omega enter the modified Kepler equation |
| Ch 11 | Delaunay variable definitions | Delaunay variables derived from (a, e, i, Omega, omega, M) |
| Ch 12 | Osculating vs. mean element framework | Osculating vs. mean element distinction extends the Keplerian framework |
| Ch 16 | Secular precession rates | Secular RAAN and argument-of-perigee precession via dot-Omega, dot-omega |
| Ch 20 | Osculating element reconstruction | Perifocal-to-inertial transform applied to osculating elements |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.8.1] | M | §8.5 | TLE measurement errors in Keplerian elements; $\sigma_m$ reflects TLE fitting residuals, typically ~1 km in position for fresh TLEs |
| [P.8.1] | P | §8.6 | Precision loss in DCM near $i \approx 0$ or $i \approx \pi$ (RAAN singularity) and near $e \approx 0$; dual quaternion avoids orientation singularity; eccentricity singularity addressed in Ch 10 |
| [P.8.2] | P | §8.7 | Quaternion normalization drift under multiple compositions; unit-norm constraint violated by $O(N\epsilon_{\rm mach})$; renormalization restores unit norm |
| [A.8.1] | A | §8.2 | Point-mass approximation; two-body potential neglects oblateness, tesseral harmonics, and drag; accuracy error $O(J_2) \approx 10^{-3}$ per orbit radially |
| [A.8.2] | A | §8.2 | Relativistic precession (Schwarzschild); ~20 m/day positional error at 400 km; negligible for SGP4 |
| [A.8.3] | A | §8.3 | Two-body energy neglects perturbations; semi-major axis decays under drag (Ch 22) |
| [A.8.4] | A | §8.5 | Osculating vs. mean elements; osculating are instantaneous best-fit Keplerian; mean elements defined in Ch 12 |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 10 |
| Theorems | 8 |
| Lemmas | 0 |
| Corollaries | 7 |
| Propositions | 0 |
| Examples | 6 |
| Error Notes | 7 |
| Equations | ~25 |
| Sections | 7 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §8.1 Introduction | Draft | |
| §8.2 Newton's Law and the Equations of Motion | Draft | |
| §8.3 Conservation Laws | Draft | |
| §8.4 The Orbit Equation | Draft | |
| §8.5 Keplerian Orbital Elements | Draft | |
| §8.6 Position and Velocity in the Perifocal Frame | Draft | |
| §8.7 State Framework: Perifocal-to-Inertial Transform | Draft | |
