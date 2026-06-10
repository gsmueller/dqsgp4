# Derivation 023: State Transformation via 7×7 Homogeneous Matrices

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Derive the **7×7 homogeneous transformation matrix** framework for converting
satellite state vectors between reference frames. This replaces the traditional
approach of separate 3×3 rotations, translations, and transport-term corrections
with a single composable linear operation:

$$\mathbf{x}_{new} = T \cdot \mathbf{x}_{old}$$

where $\mathbf{x} = [r_x, r_y, r_z, v_x, v_y, v_z, 1]^T$ is the extended state vector.

### Matrix Topology

The 7×7 matrix is partitioned as:

$$T = \begin{bmatrix} R & \mathbf{0} & \Delta\vec{r} \\ \dot{R} & R & \Delta\vec{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}$$

| Block | Size | Function |
|-------|------|----------|
| $R$ | 3×3 | Direction Cosine Matrix (primary rotation) |
| $\mathbf{0}$ | 3×3 | Null — position does not depend on velocity |
| $\Delta\vec{r}$ | 3×1 | Origin offset (e.g., geocenter → ground station) |
| $\dot{R}$ | 3×3 | Transport term: $[\vec{\omega}]_\times R$ (skew-symmetric rotation rate) |
| $R$ | 3×3 | Velocity rotation (same DCM as position) |
| $\Delta\vec{v}$ | 3×1 | Linear velocity of new frame origin |
| Bottom row | 1×7 | $[0, 0, 0, 0, 0, 0, 1]$ (homogeneous fill) |

### Application Chain

The framework supports composable reference frame transformations:

1. **SGP4 output → TEME:** Perifocal-to-TEME via 3-1-3 Euler rotation $T_3(-\Omega)\cdot T_1(-i)\cdot T_3(-u)$.
   Pure rotation: $\dot{R} = 0$, $\Delta\vec{r} = 0$, $\Delta\vec{v} = 0$.

2. **TEME → ECEF:** Earth rotation via $T_3(\theta_{GMST})$.
   Rotating frame: $\dot{R} = [\vec{\omega}_E]_\times R \neq 0$, $\Delta\vec{r} = 0$, $\Delta\vec{v} = 0$.

3. **ECEF → SEZ (topocentric):** Ground station orientation and location.
   Rotation + translation: $R$ from geodetic coordinates, $\Delta\vec{r}$ = ground station ECEF position,
   $\Delta\vec{v}$ = ground station velocity (from Earth rotation).

Total transformation: $T_{total} = T_{ECEF \to SEZ} \cdot T_{TEME \to ECEF} \cdot T_{perifocal \to TEME}$

The range rate (Doppler) emerges automatically from the velocity block of the
matrix multiplication — no manual $\vec{\omega} \times \vec{r}$ cross-product terms needed.

**Code:** `src/orbit/state_from_elements.h` (current 3×3 implementation to be extended)

---

## Source Documents Required

| Source | Location | Availability | What it provides |
|--------|----------|:------------:|-----------------|
| Spacetrack Report No. 3 pp. 14-15 | `Spacetrack_Report_No3...pdf` | ✓ | Current M, N, U, V unit vectors and position/velocity formulas |
| Battin (1999) Ch. 2 | — | NOT in repo | Perifocal frame definition, standard rotation sequence |
| Vallado (2013) *Fundamentals of Astrodynamics* Ch. 3 | — | NOT in repo | TEME/ECEF/SEZ frame definitions, transport theorem |
| Craig (2005) *Introduction to Robotics* Ch. 2 | — | NOT in repo | 4×4 homogeneous transformation matrix theory (robotics standard) |
| IERS Conventions (2010) | — | NOT in repo | Precise TEME→ITRF transformation (precession/nutation) |

---

## Equations to Resolve

### Part A: Elementary 7×7 Matrices

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| Elementary rotation $T_1(\alpha)$ | NOT DERIVED | 7×7 matrix for rotation about x-axis: $R_1(\alpha)$ in rotation block, $\dot{R}_1 = 0$, $\Delta\vec{r} = 0$ |
| Elementary rotation $T_3(\alpha)$ | NOT DERIVED | 7×7 matrix for rotation about z-axis |
| Translation matrix $T_{trans}(\vec{d})$ | NOT DERIVED | Identity rotation, $\Delta\vec{r} = \vec{d}$, $\Delta\vec{v} = 0$ |
| Rotating frame matrix | NOT DERIVED | $\dot{R} = [\vec{\omega}]_\times R$ for constant angular velocity $\vec{\omega}$ |
| Composition rule | NOT DERIVED | Show $T_{AB} \cdot T_{BC} = T_{AC}$ with correct block products |

### Part B: Perifocal → TEME (SGP4 Output)

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| 3-1-3 Euler decomposition | NOT DERIVED | $T_{TEME \leftarrow peri} = T_3(-\Omega) \cdot T_1(-i) \cdot T_3(-u)$ |
| Explicit matrix entries | NOT DERIVED | The 9 DCM entries as products of sin/cos of $\Omega$, $i$, $u$ |
| Equivalence to [SR3] M, N, U, V vectors | NOT DERIVED | Show M, N are columns of the DCM, U = M sin u + N cos u |
| Velocity in perifocal frame | NOT DERIVED | $\vec{v}_{peri} = [\dot{r}, r\dot{f}, 0]^T$ — derive from $\vec{r} = r\hat{u}$ |
| Unit conversion | NOT DERIVED | ER/min → km/s factor: $a_{E,km}/60$ |

### Part C: TEME → ECEF (Earth Rotation)

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| GMST rotation | NOT DERIVED | $R = R_3(\theta_{GMST})$, $\dot{R} = [\omega_E \hat{z}]_\times R$ |
| Transport term | NOT DERIVED | $\dot{R}\vec{r}$ gives the $\vec{\omega}_E \times \vec{r}$ velocity correction |
| $\Delta\vec{v}$ = 0 for geocentric | NOT DERIVED | No translational velocity offset for concentric frames |

### Part D: ECEF → SEZ (Topocentric/Ground Station)

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| SEZ rotation from geodetic lat/lon | NOT DERIVED | $R$ from $(\phi, \lambda)$ of ground station |
| Ground station ECEF position | NOT DERIVED | $\Delta\vec{r}$ from ellipsoidal coordinates → Cartesian |
| Ground station velocity | NOT DERIVED | $\Delta\vec{v} = \vec{\omega}_E \times \vec{r}_{station}$ |
| Look angles from SEZ state | NOT DERIVED | Azimuth = atan2(E, S), Elevation = asin(Z/|r|) |
| Range rate (Doppler) | NOT DERIVED | $\dot{\rho} = \vec{v}_{rel} \cdot \hat{\rho}$ — show this emerges from matrix multiplication |

### Part E: Precision and Accuracy

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| Trig evaluation order | NOT DERIVED | Precision loss from computing sin(Ω+δΩ) vs sin Ω cos δΩ + cos Ω sin δΩ |
| Near-polar singularity | NOT DERIVED | Gimbal lock at $i = 0$ or $i = \pi$ — does the 7×7 form mitigate? |
| Error propagation through T | NOT DERIVED | How TrackedValue errors propagate through the 7×7 multiplication |
| TEME vs J2000/GCRF | NOT DERIVED | Accuracy limitation of TEME frame (~1 arcsec from J2000) |

### Part F: Extensions (LVLH / Clohessy-Wiltshire)

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| LVLH frame definition | NOT DERIVED | $\hat{R} = \vec{r}/r$, $\hat{N} = \vec{r}\times\vec{v}/|\vec{r}\times\vec{v}|$, $\hat{T} = \hat{N}\times\hat{R}$ |
| TEME → LVLH matrix | NOT DERIVED | $R$ from target satellite state, $\Delta\vec{r}$ = target position, $\Delta\vec{v}$ = target velocity |
| Non-inertial terms | NOT DERIVED | $\dot{R}$ includes orbital angular velocity of LVLH frame |

---

## Assessment

The 7×7 homogeneous matrix framework requires no obscure sources — it is a
well-established technique from robotics (Craig 2005) and spacecraft attitude
dynamics (Schaub & Junkins 2009). The novel contribution is applying it
systematically to the SGP4 output chain, unifying the traditionally separate
rotation/translation/transport-term steps.

Parts A-B are achievable with [SR3] and standard linear algebra. Parts C-D
require the GMST computation (Derivation 014, complete) and ellipsoidal
geodesy (Derivations 003-004, complete). Part E is precision analysis using
the TrackedValue framework (Derivation 000, complete). Part F is an extension
for relative navigation.

**The primary source gap is the TEME frame definition itself** — TEME is not
precisely defined by any standards body (unlike GCRF/ITRF). Its relationship
to J2000/GCRF involves the equation of the equinoxes and is documented only
informally in Vallado (2013). For SGP4 compatibility this is acceptable, but
for enhanced precision the TEME→GCRF offset (~1 arcsec) must be characterized.
