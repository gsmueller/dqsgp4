# Chapter 2 — Reference Architecture (Scratchpad)

**Purpose:** Working document for the Design phase. Shows the concrete algebra, how operations compose, how errors propagate, and what each template section requires. Not part of the final chapter.

---

## Part 1: The Algebra

### 1.1 Quaternion Fundamentals

A quaternion has four components:

$$q = w + xi + yj + zk = (w, x, y, z)$$

where $i^2 = j^2 = k^2 = ijk = -1$.

We write $q = (w, \mathbf{v})$ where $w$ is the scalar part and $\mathbf{v} = (x, y, z)$ is the vector part.

### 1.2 The Hamilton Product

Given $q_1 = (w_1, \mathbf{v}_1)$ and $q_2 = (w_2, \mathbf{v}_2)$:

$$q_1 q_2 = (w_1 w_2 - \mathbf{v}_1 \cdot \mathbf{v}_2,\; w_1 \mathbf{v}_2 + w_2 \mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2)$$

In components:

$$q_1 q_2 = \begin{pmatrix} w_1 w_2 - x_1 x_2 - y_1 y_2 - z_1 z_2 \\ w_1 x_2 + x_1 w_2 + y_1 z_2 - z_1 y_2 \\ w_1 y_2 - x_1 z_2 + y_1 w_2 + z_1 x_2 \\ w_1 z_2 + x_1 y_2 - y_1 x_2 + z_1 w_2 \end{pmatrix}$$

### 1.3 The Hamilton Product as a Matrix Multiplication

**This is the key structural insight.** The product $q_1 q_2$ is linear in $q_2$ for fixed $q_1$. Therefore it can be written as a matrix times a vector:

$$q_1 q_2 = L(q_1) \, q_2$$

where $L(q_1)$ is the **left multiplication matrix**:

$$L(q_1) = \begin{pmatrix} w_1 & -x_1 & -y_1 & -z_1 \\ x_1 & w_1 & -z_1 & y_1 \\ y_1 & z_1 & w_1 & -x_1 \\ z_1 & -y_1 & x_1 & w_1 \end{pmatrix}$$

Similarly, $q_1 q_2$ is linear in $q_1$ for fixed $q_2$:

$$q_1 q_2 = R(q_2) \, q_1$$

where $R(q_2)$ is the **right multiplication matrix**:

$$R(q_2) = \begin{pmatrix} w_2 & -x_2 & -y_2 & -z_2 \\ x_2 & w_2 & z_2 & -y_2 \\ y_2 & -z_2 & w_2 & x_2 \\ z_2 & y_2 & -x_2 & w_2 \end{pmatrix}$$

**Verification:** $L(q_1)$ and $R(q_2)$ commute: $L(q_1) R(q_2) = R(q_2) L(q_1)$, which corresponds to the associativity $(q_1 q_2) q_3 = q_1 (q_2 q_3)$.

**Why this matters:** Since the Hamilton product is a matrix multiplication, it is a **linear operation** on the second (or first) argument. This means the parallel-row error architecture works: applying $L(q_1)$ to the error vector is the correct error propagation, by linearity.

### 1.4 Conjugate, Norm, Inverse

$$q^* = (w, -\mathbf{v}) = (w, -x, -y, -z)$$

$$\|q\|^2 = q q^* = w^2 + x^2 + y^2 + z^2$$

$$q^{-1} = \frac{q^*}{\|q\|^2}$$

For a **unit quaternion** ($\|q\| = 1$): $q^{-1} = q^*$.

### 1.5 Rotation Action

A unit quaternion $q$ rotates a vector $\mathbf{v} \in \mathbb{R}^3$ by embedding $\mathbf{v}$ as a pure quaternion $p = (0, \mathbf{v})$:

$$p' = q \, p \, q^*$$

The result $p'$ is also a pure quaternion $(0, \mathbf{v}')$, and $\mathbf{v}'$ is $\mathbf{v}$ rotated by angle $2\theta$ about axis $\hat{\mathbf{v}}_q / \|\hat{\mathbf{v}}_q\|$, where $q = (\cos\theta, \sin\theta\,\hat{\mathbf{e}})$.

**As a matrix operation:** The sandwich product $q \, p \, q^*$ can be written using the $L$ and $R$ matrices:

$$q \, p \, q^* = L(q) \, R(q^*) \, p$$

The $4 \times 4$ matrix $M_{\mathrm{rot}}(q) = L(q) \, R(q^*)$ acts on the 4-vector $p = (0, v_x, v_y, v_z)$. The scalar component of the output is guaranteed to be zero (the rotation preserves the "pure quaternion" structure), so the effective operation is a $3 \times 3$ map on the vector part. This $3 \times 3$ submatrix IS the rotation matrix $R(q)$.

### 1.6 Elementary Rotations

Rotation by angle $\alpha$ about axis $\hat{\mathbf{e}}_k$:

$$q_k(\alpha) = \left(\cos\frac{\alpha}{2},\; \sin\frac{\alpha}{2}\,\hat{\mathbf{e}}_k\right)$$

Explicitly:

$$q_1(\alpha) = \left(\cos\frac{\alpha}{2},\; \sin\frac{\alpha}{2},\; 0,\; 0\right)$$

$$q_2(\alpha) = \left(\cos\frac{\alpha}{2},\; 0,\; \sin\frac{\alpha}{2},\; 0\right)$$

$$q_3(\alpha) = \left(\cos\frac{\alpha}{2},\; 0,\; 0,\; \sin\frac{\alpha}{2}\right)$$

The 3-1-3 Euler sequence $R_3(-\Omega) R_1(-i) R_3(-u)$ becomes:

$$q_{\mathrm{3-1-3}} = q_3(-\Omega) \cdot q_1(-i) \cdot q_3(-u)$$

This is three quaternion multiplications (or two, since the first two can be composed first).

### 1.7 Double Cover

$q$ and $-q$ produce the same rotation: $q p q^* = (-q) p (-q)^*$. This is a topological fact (SO(3) $\cong$ SU(2)/$\mathbb{Z}_2$). It does not affect computations — we conventionally choose $w \geq 0$.

### 1.7a Quaternion from Rotation Matrix (Shepperd's Method)

Given a 3×3 rotation matrix $R$, extract the quaternion $q = (w, x, y, z)$:

1. Compute the trace and diagonal: $T = R_{11} + R_{22} + R_{33}$
2. Find the largest of $\{T, R_{11}, R_{22}, R_{33}\}$ to determine which component to extract first (avoids dividing by near-zero)
3. Extract the dominant component via a square root, then the remaining three via ratios

The four cases:
- If $T$ is largest: $w = \frac{1}{2}\sqrt{1+T}$, then $x = (R_{32}-R_{23})/(4w)$, etc.
- If $R_{11}$ is largest: $x = \frac{1}{2}\sqrt{1+R_{11}-R_{22}-R_{33}}$, then $w = (R_{32}-R_{23})/(4x)$, etc.
- Similarly for $R_{22}$ and $R_{33}$

This is numerically stable because the branch selection guarantees the denominator is bounded away from zero. The square root is the only nonlinear step (sensitivity bound from Ch 1, Corollary 1.4.4).

### 1.7b Convention: Right-Handed Coordinates

All coordinate systems throughout this project are right-handed. Rotation by angle $\alpha$ about axis $\hat{\mathbf{e}}$ follows the right-hand rule: curl the fingers of the right hand in the direction of positive $\alpha$, and the thumb points along $\hat{\mathbf{e}}$. The quaternion encoding $q = (\cos(\alpha/2), \sin(\alpha/2)\hat{\mathbf{e}})$ is consistent with this convention.

---

### 1.8 Dual Numbers

A dual number is $a + \varepsilon b$ where $\varepsilon^2 = 0$ (the dual unit). Arithmetic:

$$(a_1 + \varepsilon b_1)(a_2 + \varepsilon b_2) = a_1 a_2 + \varepsilon(a_1 b_2 + b_1 a_2)$$

The dual part encodes first-order perturbation information — it is essentially automatic differentiation in algebraic form.

### 1.9 Dual Quaternion

A dual quaternion is:

$$\hat{q} = q_r + \varepsilon\, q_d$$

where $q_r$ and $q_d$ are both quaternions (8 components total).

**Unit dual quaternion constraints:**
1. $\|q_r\| = 1$ (the real part is a unit quaternion)
2. $q_r \cdot q_d = w_r w_d + x_r x_d + y_r y_d + z_r z_d = 0$ (orthogonality)

This gives 8 components − 2 constraints = **6 degrees of freedom** (matching a rigid body configuration: 3 position + 3 orientation).

### 1.10 Dual Quaternion Product

$$\hat{q}_1 \hat{q}_2 = (q_{r1} + \varepsilon q_{d1})(q_{r2} + \varepsilon q_{d2}) = q_{r1} q_{r2} + \varepsilon(q_{r1} q_{d2} + q_{d1} q_{r2})$$

Since $\varepsilon^2 = 0$, the dual-dual cross term vanishes.

In components: this is two quaternion multiplications for the real part and two quaternion multiplications plus an addition for the dual part — total: 3 quaternion multiplications + 1 quaternion addition.

### 1.11 Dual Quaternion Product as an 8×8 Matrix

The product $\hat{q}_1 \hat{q}_2$ is linear in $\hat{q}_2$ for fixed $\hat{q}_1$. Write it as:

$$\hat{q}_1 \hat{q}_2 = \hat{L}(\hat{q}_1)\, \hat{q}_2$$

where $\hat{q}_2$ is treated as an 8-vector $(q_{r2}, q_{d2})$ and:

$$\hat{L}(\hat{q}_1) = \begin{pmatrix} L(q_{r1}) & \mathbf{0}_{4 \times 4} \\ L(q_{d1}) & L(q_{r1}) \end{pmatrix}$$

This is a **block lower-triangular** $8 \times 8$ matrix. The upper-left block rotates the real part. The lower-left block mixes the real part of $\hat{q}_2$ into the dual part of the result (this is where translation enters). The lower-right block rotates the dual part.

**This is the dual quaternion analog of the 4×4 homogeneous matrix:** the translation is absorbed into the matrix entries (the $L(q_{d1})$ block), and the entire rotation+translation is a single linear matrix-vector multiply.

Similarly for right multiplication:

$$\hat{q}_1 \hat{q}_2 = \hat{R}(\hat{q}_2)\, \hat{q}_1$$

$$\hat{R}(\hat{q}_2) = \begin{pmatrix} R(q_{r2}) & \mathbf{0}_{4 \times 4} \\ R(q_{d2}) & R(q_{r2}) \end{pmatrix}$$

### 1.12 Encoding Rotation + Translation

Given rotation quaternion $q_r$ and translation vector $\mathbf{t}$:

$$q_d = \frac{1}{2}(0, \mathbf{t}) \cdot q_r = \frac{1}{2} t \, q_r$$

where $t = (0, t_x, t_y, t_z)$ is the translation embedded as a pure quaternion.

**Extracting translation back:**

$$\mathbf{t} = 2\, q_d \, q_r^*$$

(take the vector part of the result).

**Extracting rotation:** $q_r$ is the real part of $\hat{q}$.

### 1.13 Dual Quaternion Conjugates

Three distinct conjugates:

1. **Quaternion conjugate:** $\hat{q}^* = q_r^* + \varepsilon q_d^*$
2. **Dual conjugate:** $\bar{\hat{q}} = q_r - \varepsilon q_d$
3. **Combined conjugate:** $\hat{q}^\dagger = q_r^* - \varepsilon q_d^*$

The rotation+translation action uses the combined conjugate:

$$p' = \hat{q}\, p\, \hat{q}^\dagger$$

### 1.14 Dual Quaternion Inverse

For a unit dual quaternion:

$$\hat{q}^{-1} = \hat{q}^\dagger = q_r^* - \varepsilon q_d^*$$

---

## Part 2: The State Vector

### 2.1 The 14-Component Layout

| Index | Component | Physical meaning |
|-------|-----------|-----------------|
| 0 | $w_r$ | Rotation scalar (cos half-angle) |
| 1 | $x_r$ | Rotation vector, x |
| 2 | $y_r$ | Rotation vector, y |
| 3 | $z_r$ | Rotation vector, z |
| 4 | $w_d$ | Dual scalar (encodes translation) |
| 5 | $x_d$ | Dual vector, x (encodes translation) |
| 6 | $y_d$ | Dual vector, y (encodes translation) |
| 7 | $z_d$ | Dual vector, z (encodes translation) |
| 8 | $v_x$ | Translational velocity, x |
| 9 | $v_y$ | Translational velocity, y |
| 10 | $v_z$ | Translational velocity, z |
| 11 | $\omega_x$ | Angular velocity, x |
| 12 | $\omega_y$ | Angular velocity, y |
| 13 | $\omega_z$ | Angular velocity, z |

**Constraints:** $w_r^2 + x_r^2 + y_r^2 + z_r^2 = 1$ and $w_r w_d + x_r x_d + y_r y_d + z_r z_d = 0$.

### 2.2 Construction from Orbital Elements

Given SGP4 outputs: orbital angles $(\Omega, i, u)$, position in perifocal frame $(r_{\mathrm{PF}}, 0, 0)$ rotated by argument of latitude, velocities $(\dot{r}, r\dot{f})$:

1. Build rotation quaternion: $q_r = q_3(-\Omega) \cdot q_1(-i) \cdot q_3(-u)$
2. Position in TEME: $\mathbf{r} = q_r (0, r_{\mathrm{PF}}, 0, 0) q_r^*$ (extract vector part)
3. Build dual part: $q_d = \frac{1}{2}(0, \mathbf{r}) \cdot q_r$
4. Velocity in TEME: $\mathbf{v} = q_r (0, \dot{r}, r\dot{f}, 0) q_r^*$ (extract vector part)
5. Angular velocity: $\boldsymbol{\omega} = (0, 0, 0)$ for SGP4 standard mode

### 2.3 Extraction

- **Position:** $\mathbf{r} = 2\, q_d\, q_r^*$ (vector part)
- **Orientation:** $q_r$ directly (the real part of $\hat{q}$)
- **Velocity:** components 8–10 directly
- **Angular velocity:** components 11–13 directly

### 2.4 SGP4 Reduction

For SGP4 standard mode (point-mass satellite):
- $q_r$ encodes only the frame rotation (not body attitude)
- $(w_r, x_r, y_r, z_r)$ = rotation from perifocal to TEME
- $(\omega_x, \omega_y, \omega_z) = (0, 0, 0)$
- Active components: 8 (dual quaternion) + 3 (velocity) = 11 physically meaningful; the angular velocity slots are carried but trivially zero

---

## Part 3: Operations

### 3.1 Transforming the State Between Frames

Given a frame transform represented by dual quaternion $\hat{q}_T$ with angular velocity $\boldsymbol{\omega}_T$:

**Configuration transforms by dual quaternion product:**

$$\hat{q}' = \hat{q}_T \cdot \hat{q}_{\mathrm{state}}$$

This is the $8 \times 8$ matrix $\hat{L}(\hat{q}_T)$ applied to the 8-vector $\hat{q}_{\mathrm{state}}$.

**Velocity transforms with transport term:**

$$\mathbf{v}' = q_{r,T}\, \mathbf{v}\, q_{r,T}^* + \boldsymbol{\omega}_T \times \mathbf{r}'$$

The first term rotates the velocity. The second term is the transport coupling.

In matrix form, using the $3 \times 3$ rotation matrix $R(q_{r,T})$ extracted from the quaternion:

$$\mathbf{v}' = R \, \mathbf{v} + [\boldsymbol{\omega}_T]_\times \, \mathbf{r}'$$

This is a linear operation on $(\mathbf{v}, \mathbf{r}')$.

**Angular velocity transforms by rotation:**

$$\boldsymbol{\omega}' = R \, \boldsymbol{\omega}_{\mathrm{state}} + \boldsymbol{\omega}_T$$

### 3.2 The Full State Transform as a Single Linear Operation

Collecting all components, the full 14-component state transform is:

$$\begin{pmatrix} \hat{q}' \\ \mathbf{v}' \\ \boldsymbol{\omega}' \end{pmatrix} = \begin{pmatrix} \hat{L}(\hat{q}_T) & \mathbf{0}_{8 \times 3} & \mathbf{0}_{8 \times 3} \\ \mathbf{V} & R & \mathbf{0}_{3 \times 3} \\ \mathbf{0}_{3 \times 8} & \mathbf{0}_{3 \times 3} & R \end{pmatrix} \begin{pmatrix} \hat{q} \\ \mathbf{v} \\ \boldsymbol{\omega} \end{pmatrix} + \begin{pmatrix} \mathbf{0}_8 \\ \mathbf{0}_3 \\ \boldsymbol{\omega}_T \end{pmatrix}$$

where $\mathbf{V}$ is the $3 \times 8$ matrix that extracts position from $\hat{q}$ and applies $[\boldsymbol{\omega}_T]_\times$ to it (the transport term). The additive term $\boldsymbol{\omega}_T$ is a constant offset.

**Note:** The additive $\boldsymbol{\omega}_T$ offset makes this affine, not purely linear. To make it a single matrix multiplication, extend to homogeneous form:

$$\begin{pmatrix} \hat{q}' \\ \mathbf{v}' \\ \boldsymbol{\omega}' \\ 1 \end{pmatrix} = \begin{pmatrix} \hat{L}(\hat{q}_T) & \mathbf{0} & \mathbf{0} & \mathbf{0}_8 \\ \mathbf{V} & R & \mathbf{0} & \mathbf{0}_3 \\ \mathbf{0} & \mathbf{0} & R & \boldsymbol{\omega}_T \\ \mathbf{0}^T & \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \begin{pmatrix} \hat{q} \\ \mathbf{v} \\ \boldsymbol{\omega} \\ 1 \end{pmatrix}$$

This is a **15×15 matrix** acting on a 15-vector. The homogeneous coordinate absorbs the angular velocity offset, just as it absorbs translation in the 4×4 case.

**Wait** — but the translation is already absorbed in the dual quaternion (the $\hat{L}$ block). The only remaining affine term is $\boldsymbol{\omega}_T$. If we include the homogeneous coordinate, the total is 15 components.

**Alternatively:** If $\boldsymbol{\omega}_T$ is treated as a known constant applied outside the matrix (not embedded), the state is 14 components and the matrix is 14×14. The angular velocity offset is handled as a post-addition, which is acceptable because it's a known constant with known error — it doesn't need the parallel-row treatment.

**Decision point:** Is the angular velocity offset common enough to justify the 15th dimension? For SGP4, the only angular velocity offset is Earth rotation $\boldsymbol{\omega}_E$ in the TEME→PEF transform. All other transforms have $\boldsymbol{\omega}_T = 0$.

### 3.3 Transport Theorem from Differentiation

The transform $\hat{q}_T(t)$ varies with time as the frame rotates. The time derivative of the transformed configuration is:

$$\frac{d}{dt}\left(\hat{q}_T(t) \cdot \hat{q}_{\mathrm{state}}\right) = \dot{\hat{q}}_T(t) \cdot \hat{q}_{\mathrm{state}} + \hat{q}_T(t) \cdot \dot{\hat{q}}_{\mathrm{state}}$$

The first term ($\dot{\hat{q}}_T \cdot \hat{q}_{\mathrm{state}}$) is the transport contribution: the frame is rotating, which changes the observed position even if the body isn't moving. The angular velocity is extracted via:

$$\boldsymbol{\omega}_T = 2\, q_{r,T}^* \, \dot{q}_{r,T}$$

and the transport velocity is $\boldsymbol{\omega}_T \times \mathbf{r}'$.

The second term ($\hat{q}_T \cdot \dot{\hat{q}}_{\mathrm{state}}$) is the body's own velocity rotated into the new frame.

**This is exactly the transport theorem**, obtained from the product rule of dual quaternion multiplication — not from a separate physical argument.

**Assumption: Instantaneous angular velocity.** The derivation assumes $\boldsymbol{\omega}_T$ is constant during a single transform step. For Earth rotation, $\omega_E \approx 7.292115 \times 10^{-5}$ rad/s varies by $< 10^{-8}$ rad/s over any SGP4 propagation step, so this is an excellent approximation.

**Quantitative transport error example:** For a position error $\delta r = 1$ km at LEO:

$$\delta v_{\mathrm{transport}} = \omega_E \cdot \delta r = 7.29 \times 10^{-5} \times 1 = 7.29 \times 10^{-5} \text{ km/s} \approx 0.073 \text{ m/s}$$

This is the dominant mechanism by which position error couples into velocity error through the pipeline. For the error rows, the same $[\boldsymbol{\omega}_E]_\times$ matrix multiplies $\delta(\mathbf{r})$ to give $\delta(\mathbf{v}_{\mathrm{transport}})$ — no separate formula needed (linearity of the cross product in $\mathbf{r}$).

### 3.4 The SGP4 Pipeline

**Transform 1: Perifocal → TEME**

Pure rotation, no translation, no frame angular velocity:

$$q_T = q_3(-\Omega) \cdot q_1(-i) \cdot q_3(-u)$$
$$\hat{q}_T = q_T + \varepsilon\, \mathbf{0}$$
$$\boldsymbol{\omega}_T = \mathbf{0}$$

The $\hat{L}$ matrix is:

$$\hat{L} = \begin{pmatrix} L(q_T) & \mathbf{0} \\ \mathbf{0} & L(q_T) \end{pmatrix}$$

(block diagonal — no translation mixing).

**Transform 2: TEME → PEF**

Rotation by $\theta_{\mathrm{GMST}}$ about $\hat{\mathbf{e}}_3$, with frame angular velocity $\boldsymbol{\omega}_E = \omega_E \hat{\mathbf{e}}_3$:

$$q_T = q_3(\theta_{\mathrm{GMST}})$$
$$\hat{q}_T = q_T + \varepsilon\, \mathbf{0}$$
$$\boldsymbol{\omega}_T = (0, 0, \omega_E)$$

The velocity transform includes the transport term:

$$\mathbf{v}_{\mathrm{PEF}} = R_3(\theta) \, \mathbf{v}_{\mathrm{TEME}} + \boldsymbol{\omega}_E \times \mathbf{r}_{\mathrm{PEF}}$$

**Composition: Perifocal → PEF**

$$\hat{q}_{\mathrm{PF} \to \mathrm{PEF}} = \hat{q}_{\mathrm{TEME} \to \mathrm{PEF}} \cdot \hat{q}_{\mathrm{PF} \to \mathrm{TEME}}$$

The composed velocity transform accumulates both the rotation and the transport term from the second transform acting on the position produced by the first.

### 3.5 Near-Identity Perturbations

A quaternion near identity: $q = (1, \epsilon_1, \epsilon_2, \epsilon_3)$ where $|\epsilon_i| \ll 1$.

To maintain unit norm: $w = \sqrt{1 - \epsilon_1^2 - \epsilon_2^2 - \epsilon_3^2} \approx 1 - \frac{1}{2}(\epsilon_1^2 + \epsilon_2^2 + \epsilon_3^2)$.

Composition of two near-identity quaternions:

$$q_1 q_2 \approx (1, \epsilon_1 + \epsilon_1', \epsilon_2 + \epsilon_2', \epsilon_3 + \epsilon_3')$$

to first order in the small components. Composition is approximately **additive in the perturbation**.

For dual quaternions, the same applies: $\hat{q}_1 \hat{q}_2 \approx \hat{1} + \delta\hat{q}_1 + \delta\hat{q}_2$ where $\delta\hat{q}$ is the deviation from identity.

This is the mathematical basis for treating short-period and long-period corrections as additive.

---

## Part 4: Error Calculus

### 4.1 The Parallel-Row Structure

$$\begin{pmatrix} \text{Row 1: physical} \\ \text{Row 2: } \sigma_m \\ \text{Row 3: } \delta_p \\ \text{Row 4: } \delta_a \end{pmatrix} = \begin{pmatrix} \hat{q}_1, \hat{q}_2, \ldots, \hat{q}_8, v_x, v_y, v_z, \omega_x, \omega_y, \omega_z \\ \delta_m(\hat{q}_1), \delta_m(\hat{q}_2), \ldots, \delta_m(\omega_z) \\ \delta_p(\hat{q}_1), \delta_p(\hat{q}_2), \ldots, \delta_p(\omega_z) \\ \delta_a(\hat{q}_1), \delta_a(\hat{q}_2), \ldots, \delta_a(\omega_z) \end{pmatrix}$$

Each row is a 14-component vector. The same transform matrix acts on each row.

### 4.2 Why the Same Matrix Works (Linear Operations)

The core of the parallel-row principle:

**For quaternion multiplication:** $q_1 q_2 = L(q_1) q_2$. The matrix $L(q_1)$ depends on $q_1$ (the transform), not on $q_2$ (the thing being transformed). So:

- Physical: $q_{\mathrm{out}} = L(q_T) \, q_{\mathrm{state}}$
- Error: $\delta(q_{\mathrm{out}}) = L(q_T) \, \delta(q_{\mathrm{state}})$

Same matrix $L(q_T)$, different input vector. Valid because $L$ is a linear map.

**For the transport term:** $\boldsymbol{\omega} \times \mathbf{r} = [\boldsymbol{\omega}]_\times \mathbf{r}$. The matrix $[\boldsymbol{\omega}]_\times$ depends on $\boldsymbol{\omega}$ (the frame rotation rate), not on $\mathbf{r}$. So:

- Physical: $\mathbf{v}_{\mathrm{transport}} = [\boldsymbol{\omega}_T]_\times \, \mathbf{r}$
- Error: $\delta(\mathbf{v}_{\mathrm{transport}}) = [\boldsymbol{\omega}_T]_\times \, \delta(\mathbf{r})$

Same matrix, different input. Valid because the cross product is linear in its second argument.

**For the full state matrix:** Each block is a linear map depending on the transform parameters, not on the state. So the entire state matrix applies to the error rows identically.

### 4.3 Why the Same Matrix Fails (Nonlinear Operations)

**Constructing $q$ from angle $\alpha$:**

$$q_3(\alpha) = (\cos(\alpha/2), 0, 0, \sin(\alpha/2))$$

This is nonlinear in $\alpha$. The error is:

$$\delta(q) \approx \frac{dq}{d\alpha} \delta(\alpha) = \frac{1}{2}(-\sin(\alpha/2), 0, 0, \cos(\alpha/2)) \, \delta(\alpha)$$

This is NOT the same operation as constructing $q$ — it requires the derivative, which is a different computation.

**Renormalization $q/\|q\|$:**

This is nonlinear (division by the norm). The sensitivity is:

$$\delta(q/\|q\|) \approx \frac{1}{\|q\|}\left(I - \frac{q q^T}{\|q\|^2}\right) \delta(q)$$

The matrix $(I - qq^T/\|q\|^2)$ projects the error onto the tangent plane of the unit sphere. For a near-unit quaternion, $\|q\| \approx 1$ and the projection is close to identity, so $\delta(q/\|q\|) \approx \delta(q)$ — renormalization barely changes the error.

### 4.4 The Isometry Property

For a unit quaternion $q_T$ with exact components ($\delta(q_T) = 0$):

$$\|L(q_T) \, \delta(\mathbf{p})\| = \|\delta(\mathbf{p})\|$$

because $L(q_T)$ is an orthogonal matrix when $q_T$ is a unit quaternion. (Proof: $L(q_T)^T L(q_T) = \|q_T\|^2 I = I$.)

**An exact rotation does not amplify errors.** This is the quaternion analog of the DCM isometry property.

### 4.5 Error from Inexact Rotation

When the transform quaternion $q_T$ carries error $\delta(q_T)$ (from computing sin/cos of input angles):

$$\delta(q_T \, p \, q_T^*) \leq \underbrace{\|L(q_T) R(q_T^*)\| \cdot \delta(p)}_{\text{isometry: } = \delta(p)} + \underbrace{\|p\| \cdot g(\delta(q_T))}_{\text{rotation error}}$$

where $g(\delta(q_T))$ is a bound depending on the quaternion error magnitude. For small $\delta(q_T)$:

$$g(\delta(q_T)) \approx 2\|\delta(q_T)\|$$

So: $\delta(\mathbf{v}') \leq \delta(\mathbf{v}) + 2\|\mathbf{v}\| \cdot \|\delta(q_T)\|$.

### 4.6 Scalar Reduction

**Position error:** Extract position $\mathbf{r} = 2 q_d q_r^*$. This involves a quaternion multiplication, so the position error is a function of the dual quaternion errors. The scalar summary:

$$\delta(\text{position}) = \|\delta(\mathbf{r})\| = \sqrt{\delta(r_x)^2 + \delta(r_y)^2 + \delta(r_z)^2}$$

**Velocity error:**

$$\delta(\text{velocity}) = \|\delta(\mathbf{v})\| = \sqrt{\delta(v_x)^2 + \delta(v_y)^2 + \delta(v_z)^2}$$

**Orientation error:** For a unit quaternion, the angular error is approximately:

$$\delta(\text{orientation}) \approx 2\|\delta(\mathbf{v}_q)\|$$

where $\mathbf{v}_q = (x_r, y_r, z_r)$ is the vector part. The factor of 2 comes from the half-angle encoding.

**Angular velocity error:**

$$\delta(\text{angular velocity}) = \|\delta(\boldsymbol{\omega})\| = \sqrt{\delta(\omega_x)^2 + \delta(\omega_y)^2 + \delta(\omega_z)^2}$$

### 4.7 Error Through a Composition Chain

The 3-1-3 Euler rotation $q = q_3(-\Omega) \cdot q_1(-i) \cdot q_3(-u)$ is three Hamilton products. Each product uses $L(q_{\mathrm{prev}})$, and each introduces error from:
1. The exactness of $q_{\mathrm{prev}}$ (which itself carries error from angle uncertainty)
2. Rounding in the 16 multiply-adds of the Hamilton product

For the chain $q_{12} = q_1 \cdot q_2$ then $q_{123} = q_{12} \cdot q_3$:

- Error in $q_{12}$: $\delta(q_{12}) = L(q_1)\,\delta(q_2) + R(q_2)\,\delta(q_1)$ (Ch 1, Theorem 1.3.2 applied to the bilinear product, plus rounding)
- Error in $q_{123}$: $\delta(q_{123}) = L(q_{12})\,\delta(q_3) + R(q_3)\,\delta(q_{12})$

The error in $q_{12}$ feeds into the second step. By the isometry property ($\|L(q)\| = \|R(q)\| = 1$ for unit quaternions), the norm of the error does not grow: $\|\delta(q_{123})\| \leq \|\delta(q_1)\| + \|\delta(q_2)\| + \|\delta(q_3)\|$ (plus rounding terms). **Composition of exact-norm rotations accumulates errors additively, not multiplicatively.**

### 4.8 Missing Pipeline Rows

**Secular rate updates** (linear): These are additions $M' = M + \dot{M}\,\Delta t$, $\omega' = \omega + \dot{\omega}\,\Delta t$. Since addition is linear, $\delta(M') = \delta(M) + |\Delta t|\,\delta(\dot{M}) + |\dot{M}|\,\delta(\Delta t)$. The same addition applies to the error row. (Fully covered by Ch 1, Theorem 1.3.1.)

**Dot product** $s = \mathbf{a} \cdot \mathbf{b}$ (linear in each argument separately): $\delta(s) \leq \sum_i |a_i|\,\delta(b_i) + |b_i|\,\delta(a_i)$ (Ch 1, Theorems 1.3.1–1.3.2). No cancellation risk since all terms are same-sign products summed.

**Cross product** $\mathbf{c} = \mathbf{a} \times \mathbf{b}$: Each component is a difference of two products ($c_1 = a_2 b_3 - a_3 b_2$). **Subtractive cancellation risk** when $\mathbf{a} \approx k\mathbf{b}$ (nearly parallel): in this case $|a_2 b_3| \approx |a_3 b_2|$ and the subtraction loses digits. The condition number is $\kappa \approx \|\mathbf{a}\|\|\mathbf{b}\| / \|\mathbf{a} \times \mathbf{b}\|$, which diverges as the vectors become parallel. This is detected by $\mathrm{rd}(\mathbf{c}) \to 0$ (Ch 1, Definition 1.9.1). Occurs in angular momentum $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ when the orbit is nearly rectilinear.

**Euclidean norm** $\rho = \|\mathbf{v}\| = \sqrt{\sum v_i^2}$: Nonlinear (square root). Sensitivity: $\delta(\rho) \leq \sum_i (|v_i|/\rho)\,\delta(v_i) \leq \|\boldsymbol{\delta}(\mathbf{v})\|$ (Cauchy-Schwarz). The norm never amplifies the total error magnitude, but the square root is a nonlinear step requiring Ch 1, Corollary 1.4.4.

### 4.9 SGP4 Error Source Map

Tracing each error source through the quaternion construction:

| Source | Type | Enters at | Propagates via |
|--------|------|-----------|---------------|
| $\Omega$ from TLE | $\sigma_m$ (~8 digits) | $q_3(-\Omega) = (\cos(\Omega/2), 0, 0, -\sin(\Omega/2))$ | Nonlinear: $\delta(q) = \frac{1}{2}\delta(\Omega)$ in the trig components |
| $i$ from TLE | $\sigma_m$ (~8 digits) | $q_1(-i) = (\cos(i/2), -\sin(i/2), 0, 0)$ | Same |
| $u$ from TLE | $\sigma_m$ (~8 digits) | $q_3(-u)$ | Same |
| $\sin(\Omega/2)$, $\cos(\Omega/2)$ | $\delta_p$ (repr. + eval) | Elementary quaternion components | Ch 1, Cors. 1.4.1–1.4.2 |
| Hamilton product rounding | $\delta_p$ ($\sim 4\epsilon_{\mathrm{mach}}$ per product) | Each of 2 compositions in the 3-1-3 chain | §4.7 accumulation |
| $\theta_{\mathrm{GMST}}$ polynomial | $\delta_p$ (polynomial eval) | $q_3(\theta)$ for TEME→PEF | Same as angle-to-quaternion |
| $\omega_E$ (Earth rotation rate) | $\sigma_m$ (measured) | Transport term $[\boldsymbol{\omega}_E]_\times \mathbf{r}$ | Linear: $\delta(\mathbf{v}) = [\delta(\boldsymbol{\omega}_E)]_\times \mathbf{r} + [\boldsymbol{\omega}_E]_\times \delta(\mathbf{r})$ |
| TEME frame definition | $\delta_a$ (~0.1 arcsec) | The frame itself | Not reducible — entire output has this floor |
| GMST coefficients | $\sigma_m$ (IAU 1982 determination) | $\theta_{\mathrm{GMST}}$ value | Propagates through trig into rotation quaternion |

### 4.10 Linearization Validity Criterion

The parallel-row principle uses $\delta(\mathbf{y}) = A\,\delta(\mathbf{x})$ which is the first-order (linear) approximation. This is valid when:

$$\delta(\mathbf{x}) \ll |\mathbf{x}| \quad \text{i.e.,} \quad \mathrm{rd}(\mathbf{x}) \geq d_{\min}$$

where $d_{\min}$ is a threshold (typically 4, meaning the error is $< 10^{-4}$ relative to the value).

When $\mathrm{rd}(v) < d_{\min}$: the linear approximation breaks down, and the full nonlinear propagation from Ch 1 must be used — applying the exact error formulas (Theorems 1.3.1–1.3.3, Corollaries 1.4.1–1.4.10) to each scalar operation individually, rather than the matrix shortcut.

The threshold $d_{\min} = 4$ is conservative. Analysis during the Build phase may refine this for specific operations. The criterion connects to Ch 1, Definition 1.9.1 (reliable digits) and Definition 1.9.3 (comparison reliability).

### 4.11 Merged vs. Separate Error Modes

**Separate mode (4 rows):** Physical state + $\sigma_m$ row + $\delta_p$ row + $\delta_a$ row. The same transform is applied 4 times. At the end, the diagnostic tells you which error category dominates and what remedy applies:
- $\sigma_m$ dominant → need better measurements (not our problem)
- $\delta_p$ dominant → use wider arithmetic or more terms
- $\delta_a$ dominant → use a better model

**Merged mode (2 rows):** Physical state + $\delta_{\mathrm{total}}$ row where $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$ per component. The transform is applied twice. The diagnostic tells you the total error magnitude but not the category breakdown.

Both modes use the same transform code. The choice is a runtime decision — separate for diagnostic analysis, merged for production propagation.

### 4.12 Reliable-Digits Connection

The scalar reductions (§4.6) feed directly into Ch 1, Definition 1.9.1:

$$\mathrm{rd}(\text{position}) = \left\lfloor -\log_{10}\!\left(\frac{\delta(\text{position})}{\|\mathbf{r}\|}\right) \right\rfloor$$

This is the decision criterion for:
- **Formula switching** (Ch 18): when $\mathrm{rd}(\sin i) < 2$, switch to Lyddane variables
- **Convergence adequacy** (Ch 9): when $\mathrm{rd}(E) < d_{\min}$, the Kepler solver needs more iterations
- **Output quality** (Ch 38): the number of reliable digits in each output component

---

## Part 5: Affine Operations via Linear Multiplication

### 5.1 The Core Principle

**The power of homogeneous coordinates**, generalized:

In the 4×4 homogeneous matrix framework:
- Affine: $\mathbf{r}' = R\,\mathbf{r} + \mathbf{t}$ (rotation + translation = two operations)
- Linear: $\bar{\mathbf{p}}' = H\,\bar{\mathbf{p}}$ where $\bar{\mathbf{p}} = (\mathbf{r}, 1)$ and $H$ absorbs $\mathbf{t}$ (one operation)

The translation is not a separate addition — it's encoded in the matrix, so the entire affine transform is a single matrix multiplication.

In the dual quaternion framework:
- Affine: $\mathbf{r}' = q_r\,\mathbf{r}\,q_r^* + \mathbf{t}$ (rotation + translation = two operations)
- Linear: $\hat{p}' = \hat{q}\,\hat{p}\,\hat{q}^\dagger$ where $\hat{p}$ embeds $\mathbf{r}$ as a dual quaternion (one operation)

Again, the translation is not a separate addition — it's encoded in the dual part of $\hat{q}$, so the entire affine transform is a single algebraic product.

**As a matrix multiplication:** Via the 8×8 matrix $\hat{L}(\hat{q}) \cdot \hat{R}(\hat{q}^\dagger)$, the action on an embedded point is:

$$\hat{p}' = \hat{M}(\hat{q}) \, \hat{p}$$

This is an 8×8 matrix times an 8-vector. The matrix encodes both rotation and translation. Composition of two transforms is multiplication of their matrices:

$$\hat{M}(\hat{q}_1 \hat{q}_2) = \hat{M}(\hat{q}_1) \cdot \hat{M}(\hat{q}_2)$$

### 5.2 Why This Enables the Parallel-Row Principle

The parallel-row principle requires that the transform be a **linear map** on the state vector — specifically, a matrix times a vector.

- If the transform is $\mathbf{y} = A\mathbf{x}$ (linear), then $\delta(\mathbf{y}) = A\,\delta(\mathbf{x})$ — same matrix, different input.
- If the transform is $\mathbf{y} = A\mathbf{x} + \mathbf{b}$ (affine), then $\delta(\mathbf{y}) = A\,\delta(\mathbf{x})$ — the constant $\mathbf{b}$ drops out of the error. Still works, but requires knowing that $\mathbf{b}$ is exact.
- If the transform is $\mathbf{y} = f(\mathbf{x})$ (nonlinear), then $\delta(\mathbf{y}) \approx J_f\,\delta(\mathbf{x})$ — requires the Jacobian, which is a DIFFERENT matrix.

By encoding translation in the matrix (homogeneous or dual quaternion), we convert the affine case into the linear case. The error row uses the SAME matrix as the physical row, with no special treatment for the translation.

### 5.3 The Velocity Extension

The 7×7 state matrix extends this principle to velocity:

- Without homogeneous trick: $\mathbf{v}' = R\,\mathbf{v} + \boldsymbol{\omega} \times \mathbf{r}' + \Delta\mathbf{v}$ (affine in $(\mathbf{r}, \mathbf{v})$ jointly)
- With the 7×7 matrix: $(\mathbf{r}', \mathbf{v}', 1)^T = T_{7 \times 7} \, (\mathbf{r}, \mathbf{v}, 1)^T$ (single matrix multiply)

The transport term $\boldsymbol{\omega} \times \mathbf{r}'$ is absorbed into the matrix (the $\dot{R}$ block), and the velocity offset $\Delta\mathbf{v}$ is absorbed into the homogeneous column.

In the dual quaternion framework, the transport theorem falls out of the time derivative (Part 3, §3.3). The velocity transform can be expressed as a linear map on $(\hat{q}, \mathbf{v})$, with the transport coupling in the off-diagonal block $\mathbf{V}$ (Part 3, §3.2).

### 5.4 Summary of Linearity

| Transform | Representation | Linear in state? | Error row treatment |
|-----------|---------------|-----------------|-------------------|
| Rotation | Quaternion $q$ | Yes ($L(q)$ is a matrix) | Same matrix |
| Rotation + translation | Dual quaternion $\hat{q}$ | Yes ($\hat{L}(\hat{q})$ is a matrix) | Same matrix |
| State (pos + vel) | 7×7 matrix or 14-vector | Yes (block matrix) | Same matrix |
| Transport term | $[\boldsymbol{\omega}]_\times$ | Yes (skew-symmetric matrix) | Same matrix |
| Angular velocity offset | additive $\boldsymbol{\omega}_T$ | Affine (constant) | Drops out of error |

Every operation in the pipeline is either linear (same matrix applies to error) or a known constant offset (drops out). The only nonlinear steps are the CONSTRUCTION of the transform (computing sin/cos, solving Kepler's equation, etc.), not the APPLICATION of the transform.

---

## Part 6: Matrix Form Derivation

### 6.1 Rotation Matrix from Quaternion

Expand the sandwich product $q p q^*$ component by component for $q = (w, x, y, z)$ and $p = (0, p_x, p_y, p_z)$:

$$R(q) = \begin{pmatrix} 1 - 2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\ 2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\ 2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2) \end{pmatrix}$$

This IS the 3×3 rotation matrix. Every entry is a quadratic polynomial in the quaternion components.

### 6.2 State Matrix from Dual Quaternion

The dual quaternion transform on the full state produces, when written in components:

- Position block: $R(q_r)$ — the 3×3 rotation matrix (from §6.1)
- Translation: $\mathbf{t} = 2 q_d q_r^*$ — the position offset
- Velocity rotation block: $R(q_r)$ — same rotation
- Transport block: $[\boldsymbol{\omega}]_\times R(q_r)$ — from the time derivative
- Velocity offset: $\Delta\mathbf{v}$ — from the frame velocity

These assemble into:

$$T_{7 \times 7} = \begin{pmatrix} R & \mathbf{0} & \mathbf{t} \\ [\boldsymbol{\omega}]_\times R & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix}$$

which is exactly the 7×7 state matrix from the archived DCM chapter. **The matrix form is a derived consequence of the dual quaternion algebra.**

---

## Part 7: Section-by-Section Requirements

| Template Section | What it needs from this reference architecture |
|-----------------|----------------------------------------------|
| **§2.1** (Motivation) | Part 5 §5.1: the linearity principle that unifies the framework |
| **§2.2** (Rotation) | Part 1 §§1.1–1.7: quaternion algebra, Hamilton product, $L(q)$ matrix, elementary rotations, double cover |
| **§2.3** (Translation) | Part 1 §§1.8–1.14: dual numbers, dual quaternions, $\hat{L}$ matrix, encoding/extraction |
| **§2.4** (Rigid body state) | Part 2: state vector layout; Part 3 §§3.1–3.3: transform formula, transport theorem from differentiation |
| **§2.5** (SGP4 pipeline) | Part 3 §§3.4–3.5: concrete transforms, near-identity perturbations |
| **§2.6** (Error propagation) | Part 4 §§4.1–4.5: parallel rows, linear/nonlinear split, isometry, inexact rotation |
| **§2.7** (Matrix form) | Part 6: $R(q)$ derivation, 7×7 from dual quaternion |
| **§2.8** (Parallel-row architecture) | Part 4 §4.6: scalar reduction; Part 5: linearity justification |
| **§2.9** (Summary) | This table (Part 7) |

---

## Resolved Design Decisions

**State vector size: 14 components.** The angular velocity offset $\boldsymbol{\omega}_T$ (used only for Earth rotation in the TEME→PEF transform) is handled as a known-constant post-addition, not embedded in a homogeneous coordinate. Rationale: only one transform in the pipeline has a nonzero $\boldsymbol{\omega}_T$, so adding a 15th dimension for all transforms is unjustified. The post-addition contributes zero to the error rows (the offset is a known constant with its own tracked uncertainty).

**SGP4 active components: reconciliation.** The template's "7 active components" counts position(3) + velocity(3) + homogeneous(1) in the traditional formulation. The reference architecture's "11 active components" counts dual-quaternion(8) + velocity(3) because position is encoded in the dual quaternion rather than stored as raw Cartesian. Both are correct — they describe the same 6 DOF (3 position + 3 velocity) in different representations. For SGP4 standard mode, the rotation part of the dual quaternion encodes the frame rotation (not body attitude), and the dual part encodes position. The orientation DOF (body attitude) is trivially identity; the angular velocity slots are zero.

## Open Questions for Build Phase

1. **Orientation error reduction.** The factor of 2 in $\delta(\text{orientation}) \approx 2\|\delta(\mathbf{v}_q)\|$ comes from the half-angle encoding. Verify this is the correct geometric interpretation (small-angle approximation to the rotation angle error).

2. **Position extraction error.** Extracting $\mathbf{r} = 2 q_d q_r^*$ involves a quaternion multiplication, which is a nonlinear function of BOTH $q_d$ and $q_r$. The error in the extracted position depends on errors in both the dual and real parts. This needs careful treatment — it's not just "read off components 5–7."

3. **Renormalization frequency.** After how many quaternion compositions does the norm drift require renormalization? This determines whether renormalization is a rare correction or a per-step operation, which affects the error budget.
