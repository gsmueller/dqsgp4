# Chapter 2: The State Framework

**Part I: Mathematical Foundations**

**Status: DRAFT — Architectural template. No theorems yet. Next phase: Design.**

---

## §2.1 Why a Unified State Framework Is Needed

### Objective

A satellite's complete physical state includes where it is, which way it's oriented, how fast it's moving, and how fast it's rotating. The SGP4 propagation pipeline transforms this state through a sequence of reference frames — from the orbital plane, through an inertial frame, to the rotating Earth. Each transformation must handle all four aspects simultaneously, and each must propagate the three error categories (Ch 1) through the same operations.

This section establishes why a single algebraic framework is needed, rather than treating position, orientation, velocity, and angular velocity as independent objects transformed by separate mechanisms.

### Problem Statement

Without a unified framework:
- Position and velocity are transformed by separate formulas, and the velocity coupling (transport theorem) must be added as a manual correction
- Orientation is handled by Euler angles with known singularities, requiring special cases (Lyddane modification)
- Error propagation requires separate derivations for each operation, each frame transform, and each combination of state components
- Composing multiple transforms requires tracking the bookkeeping of which formula applies where

With a unified framework:
- A single algebraic product transforms the entire state (position, orientation, velocity, angular velocity)
- Composition of transforms is a product of representations — no manual bookkeeping
- The velocity coupling (transport theorem) emerges automatically from differentiating the representation
- Error propagation uses the same product applied to the error state — no separate derivations for linear steps

### Downstream Dependencies

Every subsequent chapter with a "State Matrix Formulation" section builds on the framework defined here:

| Chapter | What it needs from Ch 2 |
|---------|------------------------|
| Ch 8 (Keplerian orbit) | How to express the perifocal frame orientation |
| Ch 13 (Geopotential) | How force gradients enter the state framework |
| Ch 16–17 (Brouwer theory) | How secular updates compose with the state |
| Ch 18 (Short-period corrections) | Near-identity perturbations in the framework |
| Ch 20 (Osculating elements) | The complete mean-to-osculating state transform |
| Ch 27 (Third-body perturbations) | Third-body force in the state framework |
| Ch 29 (Sidereal time) | GMST angle as input to the Earth-rotation transform |
| Ch 30 (Coordinate transforms) | Every frame rotation and composition |
| Ch 34–35 (Propagation) | The full pipeline as a chain of composed transforms |
| Ch 38 (State vector output) | The error budget from the composed pipeline |

### Connection to Error Tracking

The three-error framework of Ch 1 operates on scalar values. This section motivates extending it to the full state: each of the 13 state components carries its own $(v, \sigma_m, \delta_p, \delta_a)$ tracked value, and the framework's transforms must propagate all three error categories through every operation.

The parallel-row architecture (§2.8) — where the same transform is applied to the physical state and to the error states — is previewed here as the design goal.

---

## §2.2 The Algebra of Rotation Without Singularity

### Objective

Introduce a representation of three-dimensional rotation that:
1. Is singularity-free for all orientations
2. Composes by a single algebraic product
3. Has a simple inverse
4. Decomposes into scalar arithmetic (so Ch 1 error rules apply)
5. Connects to the traditional rotation matrix as a derived consequence

### Required Operations

The following operations must be defined and shown to be closed (the result is in the same representation):

| Operation | Input | Output | Purpose |
|-----------|-------|--------|---------|
| **Apply rotation to vector** | rotation $q$, vector $\mathbf{v}$ | rotated vector $\mathbf{v}'$ | Transforming position/velocity between frames |
| **Compose two rotations** | rotations $q_1$, $q_2$ | combined rotation $q_{12}$ | Chaining frame transforms |
| **Invert a rotation** | rotation $q$ | inverse $q^{-1}$ | Reversing a frame transform |
| **Construct from axis + angle** | axis $\hat{\mathbf{e}}$, angle $\theta$ | rotation $q$ | Building transforms from orbital elements |
| **Construct from 3-1-3 Euler sequence** | angles $\Omega, i, u$ | rotation $q$ | The specific SGP4 frame transform |
| **Extract rotation matrix** | rotation $q$ | $3 \times 3$ matrix $R$ | Connection to existing literature and code |
| **Extract rotation from matrix** | matrix $R$ | rotation $q$ | Converting from existing code representations |

### Constraints

- **Right-handed coordinate systems.** All coordinate systems are right-handed throughout; rotation angles follow the right-hand rule. This is a global assumption inherited by all subsequent chapters.
- **No singularity at any orientation.** The Euler angle representation $(\Omega, i, u)$ is singular when $i = 0$ (the first and third rotations share an axis). The new representation must handle $i = 0$ and all other orientations uniformly.
- **Scalar decomposition.** Every operation above must reduce to a finite sequence of additions and multiplications on the representation's components. This enables direct application of Ch 1 error propagation (Theorems 1.3.1–1.3.2).
- **Closure under composition.** Composing two rotations must produce a result in the same representation, not a different type.

### Dependency on the Cross Product

The transport theorem (§2.4) involves the vector cross product $\boldsymbol{\omega} \times \mathbf{r}$ regardless of rotation representation. The algebraic properties of the cross product — its antisymmetry, its behavior under rotation, and its matrix representation — are needed when deriving the traditional matrix form in §2.7. These properties are a Build-phase concern for §2.7 specifically; they are not needed for the primary rotation framework.

### Error Responsibility

This section must establish:
- Applying an **exact** rotation to a vector preserves the error magnitude (the rotation is an isometry — this is the analog of the DCM result in Theorem 2.6.3 of the archived chapter)
- Applying an **inexact** rotation (where the representation's components carry error from computing sin/cos of input angles) introduces additional error proportional to $\|\mathbf{v}\| \cdot \delta(q)$
- The error in composing two rotations follows from the scalar product rules applied to the composition formula

---

## §2.3 Incorporating Translation into the Rotation Framework

### Objective

Extend the rotation representation to handle both rotation and translation in a single algebraic object. This is the analog of the 4×4 homogeneous matrix, but built on the singularity-free rotation representation of §2.2.

### Required Operations

| Operation | Input | Output | Purpose |
|-----------|-------|--------|---------|
| **Apply rotation + translation to point** | transform $\hat{q}$, point $\mathbf{p}$ | transformed point $\mathbf{p}'$ | Position frame transform |
| **Compose two transforms** | transforms $\hat{q}_1$, $\hat{q}_2$ | combined $\hat{q}_{12}$ | Chaining frame transforms |
| **Invert a transform** | transform $\hat{q}$ | inverse $\hat{q}^{-1}$ | Reversing a frame transform |
| **Extract rotation** | transform $\hat{q}$ | rotation $q$ | Separating rotation from translation |
| **Extract translation** | transform $\hat{q}$ | vector $\mathbf{t}$ | Separating translation from rotation |
| **Construct from rotation + translation** | rotation $q$, vector $\mathbf{t}$ | transform $\hat{q}$ | Building from known components |

### Constraints

- **Composition must be a single product** (not "rotate then add"). The algebraic product of two transforms must yield the correct combined rotation and combined translation.
- **Must be equivalent to the 4×4 homogeneous matrix** — given the same rotation and translation, both representations must produce the same transformed point.
- **Compact.** The representation should use fewer numbers than the 4×4 matrix (16 entries) while encoding the same 6 degrees of freedom.

### Error Responsibility

- Translation by an exact offset adds zero error to the error state (the offset is a known constant)
- Translation by an uncertain offset ($\Delta\mathbf{r}$ with error $\delta(\Delta\mathbf{r})$) adds that error directly to the position error row

---

## §2.4 The Complete Rigid Body State

### Objective

Define the full state of a rigid body in a single fixed-length vector:

| Component | DOF | What it represents |
|-----------|-----|-------------------|
| Position | 3 | Where the body is in space |
| Orientation | 3 DOF (represented with 4 numbers + 1 constraint) | Which way the body is pointing |
| Velocity | 3 | How fast the position is changing |
| Angular velocity | 3 | How fast the orientation is changing |
| **Total** | **12 DOF, 13 components** | |

The dual extension of the rotation algebra (§2.3) absorbs translation into the algebraic product, eliminating the need for a separate homogeneous coordinate. Affine terms (angular velocity offsets) are known constants handled as post-additions.

Derive how this state transforms between frames, with the velocity coupling (transport theorem) emerging naturally.

### The Transport Theorem

This is the central physical result of the section. When frame B rotates relative to frame A with angular velocity $\boldsymbol{\omega}$, the velocity observed in frame A differs from the velocity in frame B by the coupling term $\boldsymbol{\omega} \times \mathbf{r}$.

**Requirement:** The transport theorem must emerge from differentiating the configuration representation with respect to time. It must not be stated as an independent formula — it must be a consequence of the algebra.

This means: if the configuration is represented by some object $\hat{q}(t)$ that changes with time as the frame rotates, then $\frac{d}{dt}(\hat{q} \cdot \mathbf{state})$ must naturally produce the velocity coupling term. The product rule of the algebra must give the transport theorem.

### SGP4 Standard Mode

For SGP4 (point-mass satellite with fixed ballistic coefficient):
- Orientation = identity (no body attitude tracked)
- Angular velocity = zero (no body rotation tracked)
- The 13-component state reduces effectively to 6 active components: position (3) + velocity (3), with orientation = identity and angular velocity = zero carried trivially
- The framework's operations, restricted to this subspace, must reproduce exactly the transformations derived in the archived DCM chapter

This reduction must happen naturally from setting orientation = identity and ω = 0, without special-casing.

### Assumptions

- **Instantaneous angular velocity.** The angular velocity $\boldsymbol{\omega}$ is treated as constant during a single transformation step. Time-varying $\boldsymbol{\omega}$ requires re-evaluation of the transport term at each epoch. This assumption is valid for SGP4 because the Earth's rotation rate varies by $< 10^{-8}$ rad/s over a propagation step.

### Error Responsibility

- The transport term is the **primary error amplification mechanism** in the pipeline: it couples position error into velocity error with amplification factor $|\omega|$
- For Earth rotation: $\omega_E \approx 7.29 \times 10^{-5}$ rad/s, so a 1 km position error produces $\approx 0.073$ m/s velocity error
- This coupling is present in the error rows exactly as in the physical row (the transport term is a linear operation)

---

## §2.5 The Specific Transforms of the SGP4 Pipeline

### Objective

Express each coordinate transformation in the SGP4 pipeline as an instance of the general framework from §2.4. Show that the complete pipeline is a composition of these instances.

### The Pipeline

```
Perifocal frame        →  TEME (inertial)    →  PEF (Earth-fixed)
(orbital plane)           (True Equator,         (Pseudo Earth-Fixed,
                           Mean Equinox)          rotated by GMST)
```

**Transform 1: Perifocal → TEME**
- Pure rotation, no translation, no frame rotation rate
- Rotation angles: node $\Omega$, inclination $i$, argument of latitude $u = \omega + \nu$
- This is the 3-1-3 Euler sequence $R_3(-\Omega) R_1(-i) R_3(-u)$, expressed in the framework representation
- Error sources: $\sigma_m(\Omega, i, u)$ from TLE measurement, $\delta_p$ from sin/cos evaluation

**Transform 2: TEME → PEF**
- Rotation by GMST about the polar axis, WITH nonzero frame rotation rate $\omega_E$
- The transport theorem produces velocity coupling: $\boldsymbol{\omega}_E \times \mathbf{r}$
- Error sources: $\delta_p(\theta_{\mathrm{GMST}})$ from polynomial evaluation (Ch 29), $\delta_a$ from TEME frame limitations [A.2.1]

**Composition: Perifocal → PEF**
- Single framework product of Transform 1 and Transform 2
- The transport term from Transform 2 acts on the already-rotated position from Transform 1

### Near-Identity Perturbations

**Objective:** Many pipeline steps (short-period corrections, long-period corrections, drag adjustments) are small perturbations to the state — the transform is close to "do nothing." The framework must support a near-identity regime where composing multiple small corrections reduces to simple addition of the corrections, rather than the full nonlinear composition.

**What depends on this:**
- Ch 18 (short-period corrections): expressed as a small perturbation on top of the secular state
- Ch 19 (long-period corrections): same structure
- Ch 33 (secular update): linear-in-time rate additions
- §2.6 (error propagation): perturbative corrections are linear, so the error rows pass through them identically

**Constraint:** The framework's composition operation must have a well-defined first-order approximation for small perturbations, and the error introduced by using the first-order approximation instead of the exact composition must be quantifiable.

### Constraints

- Every angle input must be traced to its source chapter: orbital elements from TLE (Ch 31), GMST from sidereal time polynomial (Ch 29)
- The composed result must be identical to the archived DCM chapter's Theorem 2.5.3 when evaluated numerically
- The matched pair principle (Ch 3) constrains which constants may be used

### Error Source Map

| Input | Error type | Source | Magnitude |
|-------|-----------|--------|-----------|
| $\Omega$ | $\sigma_m$ | TLE measurement | ~8 decimal digits |
| $i$ | $\sigma_m$ | TLE measurement | ~8 decimal digits |
| $u$ | $\sigma_m$ | TLE measurement | ~8 decimal digits |
| $\sin(\Omega)$, $\cos(\Omega)$, etc. | $\delta_p$ | Trig evaluation | $\leq |\cos x| \cdot \delta(x)$ (Ch 1, Cor. 1.4.1) |
| $\theta_{\mathrm{GMST}}$ | $\delta_p$ | Polynomial evaluation | Ch 29 |
| TEME frame | $\delta_a$ | Truncated precession/nutation | ~0.1 arcsec |

---

## §2.6 Error Propagation Through the Framework

### Objective

Formalize two principles:

**Principle 1 (Linear operations):** For any linear operation $\mathbf{y} = f(\mathbf{x})$, the error state transforms by the same operation: $\boldsymbol{\delta}(\mathbf{y}) = f(\boldsymbol{\delta}(\mathbf{x}))$. No separate error formula is needed.

**Principle 2 (Nonlinear operations):** For a nonlinear operation $\mathbf{y} = g(\mathbf{x})$, the error state transforms by a rigorous bound on the operation's sensitivity: $\|\boldsymbol{\delta}(\mathbf{y})\| \leq B(g, \mathbf{x}) \cdot \|\boldsymbol{\delta}(\mathbf{x})\|$. The bound $B$ must be a proven upper bound, derived in the chapter where $g$ is developed.

### Pipeline Classification

Every operation in the SGP4 pipeline must be classified as linear or nonlinear for error propagation purposes. This classification determines whether the error row can use the same transform as the physical row, or requires a separate sensitivity bound.

**Linear operations** (error row uses same transform):

| Operation | Why linear | Cancellation risk? |
|-----------|-----------|-------------------|
| Frame rotation (exact angles) | Rotation is a linear map on vectors | No |
| Translation by exact offset | Addition of a known constant | No |
| Transport term ($\boldsymbol{\omega} \times \mathbf{r}$) | Cross product is bilinear in its arguments | No |
| Secular rate updates | Linear-in-time additions | No |
| Short/long-period additive corrections | Additive perturbations | No |
| Dot product | Sum of products | No |
| Cross product | Difference of products | **Yes** — when inputs are nearly parallel |

**Nonlinear operations** (error row requires a rigorous sensitivity bound):

| Operation | Why nonlinear | Bound derived in |
|-----------|--------------|-----------------|
| Frame rotation (inexact angles) | Rotation components are trigonometric functions of the angles | Ch 1 (trig bounds), §2.2 |
| Kepler equation solver | Iterative, transcendental | Ch 9 |
| Trigonometric evaluation | Nonlinear functions | Ch 1 |
| Representation renormalization | Division by a computed norm | §2.6 |
| Cube root (element recovery) | Nonlinear root extraction | Ch 32 |
| Density model (power law) | Nonlinear exponentiation | Ch 21 |
| Euclidean norm | Involves square root | Ch 1 |

### Constraints

- All bounds must be **rigorous** — proven upper bounds, not approximations
- The three error categories propagate independently through all operations (same functional form, different input magnitudes)
- **Sensitivity measure.** For each nonlinear step, a standard measure of its sensitivity to input perturbations must be established and bounded. This quantifies "how much does the output error grow relative to the input error" for that step.
- **Linearization validity.** The parallel-row principle relies on error propagation being approximately linear. There must be a well-defined criterion — connected to the reliable-digits framework of Ch 1, §1.9 — that identifies when the linear approximation breaks down and full nonlinear propagation from Ch 1 must be used instead. This criterion must be stated, not left implicit.

### Error Responsibility

This section establishes what every downstream chapter's "State Matrix Formulation" section must provide:
- For linear operations: a statement that the parallel-row principle applies
- For nonlinear operations: a rigorous bound on the sensitivity, with proof

---

## ��2.7 The Traditional Matrix Form as a Derived Result

### Objective

Derive the familiar 3×3 rotation matrix and 7×7 state matrix as component-wise expansions of the framework's operations. This connects the framework to the existing literature (Vallado, Battin, Brouwer) and existing code (`state_from_elements.h`).

### Required Derivations

1. The 3×3 rotation matrix entries expressed in terms of the framework representation's components
2. The 7×7 state matrix block structure as a consequence of the framework's state transform
3. Numerical equivalence: for any input state, the framework transform and the matrix multiply produce the same output

### Specializations

The matrix derivation must show how the general framework reduces to the specific cases used in the SGP4 pipeline:
- **Pure rotation** (no translation, no frame rotation): the simplification that applies to the Perifocal → TEME transform
- **Rotating frame** (rotation with nonzero angular velocity): the case that applies to the TEME → PEF transform, where the velocity coupling emerges
- **SGP4 reduction**: the full framework with orientation and angular velocity set to trivial values must produce exactly the traditional matrix structure used in existing SGP4 literature

### Why This Section Exists

The archived DCM chapter (626 lines, `deprecated/ch02_state_matrix_dcm_complete.md`) contains complete, proven results for the matrix form. Those results are mathematically correct. This section:
- Proves they are **consequences** of the framework (not independent results)
- Identifies which archived proofs carry over as corollaries
- Explains when the matrix form is preferable (pedagogical clarity, connection to existing code) vs. when the framework form is preferable (singularity avoidance, compact composition)

---

## §2.8 The Parallel-Row Error Architecture

### Objective

Define the complete state representation as a set of parallel rows:

| Row | Content | Length |
|-----|---------|--------|
| 1 | Physical state | 13 components |
| 2 | $\sigma_m$ error state | 13 components (same structure) |
| 3 | $\delta_p$ error state | 13 components (same structure) |
| 4 | $\delta_a$ error state | 13 components (same structure) |

All four rows pass through the same transform pipeline. For linear operations, the transform is identical. For nonlinear operations, the error rows are multiplied by the rigorous sensitivity bound instead.

### Scalar Reduction

The 13 per-component error values in each row are reduced to 4 scalar summaries:

| Summary | Components reduced | Reduction method |
|---------|-------------------|-----------------|
| $\delta(\text{position})$ | $r_x, r_y, r_z$ | Euclidean norm of component errors |
| $\delta(\text{orientation})$ | quaternion components | [To be defined — depends on quaternion error geometry] |
| $\delta(\text{velocity})$ | $v_x, v_y, v_z$ | Euclidean norm of component errors |
| $\delta(\text{angular velocity})$ | $\omega_x, \omega_y, \omega_z$ | Euclidean norm of component errors |

These scalar summaries feed into the reliable-digits criterion (Ch 1, Definition 1.9.1) for formula-switching decisions.

### Merged vs. Separate Categories

The architecture supports two modes:
- **Separate** (4 rows): full diagnostic — identifies whether error comes from measurement, precision, or model
- **Merged** (2 rows): physical state + one total error row with $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$ per component

The merged mode uses 1 extra transform pass. The separate mode uses 3 extra passes. Same transform code in both cases.

### SGP4 Standard Mode

For SGP4 (point-mass, no attitude):
- Orientation components = identity representation (error = 0)
- Angular velocity components = 0 (error = 0)
- The 13-component rows effectively reduce to 6 active components (position + velocity), with orientation and angular velocity trivially zero
- No special case in the algebra — the trivial components are simply carried along

---

## §2.9 Summary: How Subsequent Chapters Use This Framework

### Objective

Reference table for all downstream chapters. Each row identifies what the chapter needs from Ch 2 and which section provides it.

| Chapter | Framework aspect used | Ch 2 section |
|---------|---------------------|-------------|
| Ch 8 | Perifocal frame orientation | §2.5 |
| Ch 13 | Force gradient in state framework | §2.4 |
| Ch 15 | Kaula expansion in framework | §2.4 |
| Ch 16 | Secular rates as state perturbation | §2.4, §2.6 |
| Ch 17 | Second-order secular as composed perturbation | §2.4, §2.6 |
| Ch 18 | Short-period as near-identity transform | §2.5, §2.6 |
| Ch 20 | Mean-to-osculating complete transform | §2.5 |
| Ch 27 | Third-body perturbation in framework | §2.4 |
| Ch 29 | GMST angle as framework input | §2.5 |
| Ch 30 | All frame rotations and compositions | §2.2, §2.3, §2.5 |
| Ch 34 | Near-space pipeline composition | §2.4, §2.5, §2.8 |
| Ch 35 | Deep-space pipeline composition | §2.4, §2.5, §2.8 |
| Ch 38 | Error budget from composed pipeline | §2.6, §2.8 |

### Pipeline Composition Pattern

Every "State Matrix Formulation" section in a downstream chapter follows this pattern:
1. Express the chapter's operation as a framework transform (§2.2–§2.4)
2. Classify it as linear or nonlinear for error propagation (§2.6)
3. For nonlinear: provide the rigorous sensitivity bound
4. Show how it composes with adjacent pipeline steps (§2.4 composition)
5. Identify the error sources and their categories (§2.8)

---

## Error Notes

Error notes will be assigned when specific approximations are introduced in the Design phase. Placeholders from the architectural analysis:

- **[P.2.1]**: Trigonometric evaluation errors when constructing the rotation representation from Euler angles
- **[P.2.2]**: Accumulation of rounding through the representation's algebraic product
- **[P.2.3]**: Renormalization error after chains of composition
- **[A.2.1]**: TEME frame accuracy (~0.1 arcsec from truncated precession/nutation)
- **[A.2.2]**: Linearization in error propagation for nonlinear steps
- **[M.2.1]**: GMST polynomial coefficient measurement uncertainty

---

## Phase Tracking

| Phase | Status | Content |
|-------|--------|---------|
| **Draft** (this document) | **CURRENT** | Architectural template: objectives, dependencies, constraints, error responsibilities |
| **Design** | Next | Reference architecture: specific algebraic structures, operation definitions, dependency graph |
| **Build** | After Design | Full theorems with proofs, equation numbers, cross-references |
