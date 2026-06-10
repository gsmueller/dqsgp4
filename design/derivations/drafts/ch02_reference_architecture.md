# Chapter 2 — Reference Architecture v2: SU(2) Framework (Scratchpad)

**Purpose:** Working document for the Design phase. Reworked from v1 to use the 2×2 complex matrix (SU(2)) representation as the foundational algebra. All state components — rotations, positions, velocities, angular velocities — live in the same algebra of 2×2 complex matrices. Not part of the final chapter.

**Supersedes:** `ch02_reference_architecture_v1_real4x4.md` (archived).

---

## Part 1: The Algebra — Everything is a 2×2 Complex Matrix

### 1.1 The Pauli Matrices

The three Pauli matrices are the basis for encoding 3D vectors as 2×2 matrices:

$$\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$$

Together with the identity:

$$I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$$

these form a basis for all 2×2 Hermitian matrices.

Key properties:
- $\sigma_k^2 = I$ for each $k$
- $\sigma_j \sigma_k = \delta_{jk} I + i\,\varepsilon_{jkl}\,\sigma_l$ (the Pauli algebra)
- The anticommutator: $\{\sigma_j, \sigma_k\} = 2\delta_{jk}I$
- The commutator: $[\sigma_j, \sigma_k] = 2i\,\varepsilon_{jkl}\,\sigma_l$

### 1.2 Vectors as 2×2 Hermitian Matrices

Any 3D vector $\mathbf{v} = (v_x, v_y, v_z)$ maps to a 2×2 traceless Hermitian matrix:

$$V = v_x \sigma_1 + v_y \sigma_2 + v_z \sigma_3 = \begin{pmatrix} v_z & v_x - iv_y \\ v_x + iv_y & -v_z \end{pmatrix}$$

**Properties:**
- $V$ is Hermitian: $V^\dagger = V$
- $V$ is traceless: $\mathrm{tr}(V) = 0$
- The determinant encodes the squared norm: $\det(V) = -(v_x^2 + v_y^2 + v_z^2) = -\|\mathbf{v}\|^2$
- Extraction: $v_k = \frac{1}{2}\mathrm{tr}(\sigma_k V)$

This is an **isomorphism**: every 3D vector corresponds to exactly one traceless Hermitian 2×2 matrix, and vice versa.

### 1.3 Rotations as SU(2) Matrices

A rotation is a 2×2 unitary matrix with determinant 1:

$$M \in SU(2): \quad M^\dagger M = I, \quad \det M = 1$$

The general form, parameterized by $q = (w, x, y, z)$ with $w^2 + x^2 + y^2 + z^2 = 1$:

$$M = wI + i(x\sigma_1 + y\sigma_2 + z\sigma_3) = \begin{pmatrix} w + iz & ix - y \\ ix + y & w - iz \end{pmatrix}$$

Equivalently, with $\alpha = w + iz$ and $\beta = ix + y$:

$$M = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix}, \qquad |\alpha|^2 + |\beta|^2 = 1$$

**Correspondence with quaternions:** The quaternion $q = w + xi + yj + zk$ maps directly to $M$. The Hamilton product of quaternions equals the matrix product of their SU(2) matrices. This is not an analogy — it is an identity.

### 1.4 The Rotation Action: Matrix Conjugation

The rotation $M$ acts on a vector $V$ by:

$$V' = M\, V\, M^\dagger$$

This is the **same operation** for every vector — position, velocity, angular velocity, any 3D vector.

**Properties:**
- $V'$ is traceless Hermitian (rotation preserves the vector space)
- $\det(V') = \det(V)$ (rotation preserves the norm: $\|\mathbf{v}'\| = \|\mathbf{v}\|$)
- The map $V \mapsto M V M^\dagger$ is **linear in $V$**: $M(aV_1 + bV_2)M^\dagger = a\,MV_1M^\dagger + b\,MV_2M^\dagger$

The linearity in $V$ is the foundation of the parallel-row error principle.

### 1.5 The Cross Product as a Commutator

For two vectors encoded as Hermitian matrices $A$ and $B$:

$$\mathbf{a} \times \mathbf{b} \;\longleftrightarrow\; \frac{1}{2i}[A, B] = \frac{1}{2i}(AB - BA)$$

**Proof sketch:** From the Pauli algebra $\sigma_j \sigma_k = \delta_{jk}I + i\varepsilon_{jkl}\sigma_l$, the commutator $[\sigma_j, \sigma_k] = 2i\varepsilon_{jkl}\sigma_l$. Expanding $[A,B]$ with $A = \sum a_j\sigma_j$ and $B = \sum b_k\sigma_k$ gives $[A,B] = 2i\sum_{jkl}\varepsilon_{jkl}a_j b_k\sigma_l = 2i(\mathbf{a}\times\mathbf{b})\cdot\boldsymbol{\sigma}$, hence $\frac{1}{2i}[A,B] = (\mathbf{a}\times\mathbf{b})\cdot\boldsymbol{\sigma}$.

**Why this matters:** The transport term $\boldsymbol{\omega} \times \mathbf{r}$ is the commutator $\frac{1}{2i}[\Omega, R_{\mathrm{pos}}]$. It's a native operation in the matrix algebra — not an imported formula.

### 1.6 The Dot Product

$$\mathbf{a} \cdot \mathbf{b} = -\frac{1}{2}\mathrm{tr}(AB)$$

**Proof sketch:** $\mathrm{tr}(\sigma_j\sigma_k) = \mathrm{tr}(\delta_{jk}I + i\varepsilon_{jkl}\sigma_l) = 2\delta_{jk}$ (since $\mathrm{tr}(\sigma_l) = 0$). So $\mathrm{tr}(AB) = \sum_{jk}a_j b_k\,\mathrm{tr}(\sigma_j\sigma_k) = 2\sum_j a_j b_j = 2(\mathbf{a}\cdot\mathbf{b})$.

### 1.7 The Norm

$$\|\mathbf{v}\|^2 = -\det(V) = \frac{1}{2}\mathrm{tr}(V^2)$$

Both expressions give $v_x^2 + v_y^2 + v_z^2$.

### 1.8 Composition of Rotations

$$M_{AC} = M_{AB} \cdot M_{BC}$$

This is ordinary 2×2 complex matrix multiplication. The product of two SU(2) matrices is SU(2) (closure). Inverse: $M^{-1} = M^\dagger$ (adjoint = inverse for unitary matrices).

### 1.9 Elementary Rotations

Rotation by angle $\alpha$ about axis $\hat{\mathbf{e}}_k$:

$$M_k(\alpha) = \cos\frac{\alpha}{2}\,I + i\sin\frac{\alpha}{2}\,\sigma_k$$

Explicitly:

$$M_1(\alpha) = \begin{pmatrix} \cos\frac{\alpha}{2} & i\sin\frac{\alpha}{2} \\ i\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}$$

$$M_2(\alpha) = \begin{pmatrix} \cos\frac{\alpha}{2} & \sin\frac{\alpha}{2} \\ -\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}$$

$$M_3(\alpha) = \begin{pmatrix} e^{i\alpha/2} & 0 \\ 0 & e^{-i\alpha/2} \end{pmatrix}$$

Note: $M_3(\alpha)$ is diagonal — rotation about $\hat{\mathbf{e}}_3$ is just phase multiplication. This is a computational advantage.

### 1.10 Double Cover

$M$ and $-M$ produce the same rotation: $M V M^\dagger = (-M) V (-M)^\dagger$. Convention: choose the representative with $\mathrm{Re}(\alpha) = w \geq 0$.

### 1.11 Extracting Rotation from Matrix (Shepperd analog)

Given a 3×3 real rotation matrix $R$, compute $M \in SU(2)$:

From the relation $R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j M \sigma_k M^\dagger)$, extract $w, x, y, z$ using the standard Shepperd branching (choose the largest component to avoid dividing by near-zero). The mapping $M \to R$ is 2-to-1 ($\pm M$ give the same $R$); the inverse $R \to M$ is single-valued with the $w \geq 0$ convention.

---

## Part 2: Translation — The Dual Extension

### 2.1 Dual Complex Numbers

A dual complex number is $\hat{z} = z + \varepsilon z'$ where $z, z' \in \mathbb{C}$ and $\varepsilon^2 = 0$.

Arithmetic: $\hat{z}_1 \hat{z}_2 = z_1 z_2 + \varepsilon(z_1 z_2' + z_1' z_2)$.

### 2.2 Dual SU(2) Matrices

Extend the 2×2 complex matrix to have dual complex entries:

$$\hat{M} = M + \varepsilon D$$

where $M \in SU(2)$ is the rotation and $D$ is a 2×2 complex matrix encoding the translation.

**Unit dual SU(2) constraints:**
1. $M^\dagger M = I$ (the real part is unitary)
2. $M^\dagger D + D^\dagger M = 0$ (orthogonality: $D$ is in the tangent space of SU(2) at $M$)

### 2.3 Encoding Translation

Given rotation $M$ and translation vector $\mathbf{t}$, encode as:

$$D = \frac{1}{2} T \cdot M$$

where $T = t_x\sigma_1 + t_y\sigma_2 + t_z\sigma_3$ is the translation vector as a traceless Hermitian matrix.

**Extraction:**

$$T = 2\, D\, M^\dagger$$

(then read off the vector components from the Hermitian matrix).

### 2.4 The Rotation+Translation Action

The combined action on a point $P$ (encoded as traceless Hermitian):

$$P' = M\, P\, M^\dagger + T$$

where $T = 2\,D\,M^\dagger$.

In the dual framework, this is encoded as:

$$\hat{P}' = \hat{M}\, \hat{P}\, \hat{M}^\dagger$$

where $\hat{P} = P + \varepsilon\cdot 0$ embeds the point, and the dual part of the result carries the translated position.

### 2.5 Composition

$$\hat{M}_{AC} = \hat{M}_{AB} \cdot \hat{M}_{BC}$$

Expanding:
- Real part: $M_{AC} = M_{AB} M_{BC}$ (composed rotation)
- Dual part: $D_{AC} = M_{AB} D_{BC} + D_{AB} M_{BC}$ (composed translation)

This is 2×2 complex matrix multiplication over the dual numbers — the same algebraic operation, just extended.

### 2.6 Inverse

$$\hat{M}^{-1} = M^\dagger - \varepsilon\, D^\dagger$$

(adjoint in the dual algebra).

---

## Part 3: The Complete State

### 3.1 The Uniform Representation

Every component of the rigid body state is a 2×2 matrix:

| Quantity | Matrix type | Encoding | Real parameters |
|----------|------------|----------|-----------------|
| **Orientation** | $M \in SU(2)$ | Unitary, det=1 | 4 (with 1 constraint → 3 DOF) |
| **Position** | $R_{\mathrm{pos}}$: traceless Hermitian | $r_x\sigma_1 + r_y\sigma_2 + r_z\sigma_3$ | 3 |
| **Velocity** | $V$: traceless Hermitian | $v_x\sigma_1 + v_y\sigma_2 + v_z\sigma_3$ | 3 |
| **Angular velocity** | $\Omega$: traceless Hermitian | $\omega_x\sigma_1 + \omega_y\sigma_2 + \omega_z\sigma_3$ | 3 |

**Total:** 4 + 3 + 3 + 3 = 13 real parameters, 12 DOF.

Or equivalently as complex numbers: 2 (for $M$'s $\alpha, \beta$) + 3 (position reals) + 3 (velocity reals) + 3 (angular velocity reals) = 2 complex + 9 real.

**The combined configuration** (position + orientation) is the dual SU(2) matrix $\hat{M}$ with the translation encoded in the dual part (§2.3). This gives 8 real parameters (4 from $M$, 4 from $D$) with 2 constraints = 6 DOF, matching the rigid body configuration.

### 3.2 The State Transform

Given a frame transform with rotation $M_T$, angular velocity $\Omega_T$, and translation $T_T$:

**Every vector transforms by conjugation with the same $M_T$:**

$$R_{\mathrm{pos}}' = M_T\, R_{\mathrm{pos}}\, M_T^\dagger + T_T$$

$$V' = M_T\, V\, M_T^\dagger + \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$$

$$\Omega' = M_T\, \Omega\, M_T^\dagger + \Omega_T$$

$$M' = M_T \cdot M$$

**The structure:**
- Position: rotate + translate (affine)
- Velocity: rotate + transport term (the commutator $\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$)
- Angular velocity: rotate + add frame rotation rate (affine)
- Orientation: matrix multiply

All four use the **same rotation** $M_T (\cdot) M_T^\dagger$. The additional terms (translation, transport, angular velocity offset) are either known constants or commutators within the algebra.

### 3.3 The Transport Theorem as a Commutator

The velocity coupling in a rotating frame is:

$$\boldsymbol{\omega}_T \times \mathbf{r}' \;\longleftrightarrow\; \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$$

This is a **native algebraic operation**: the commutator of two matrices in the same algebra. No cross product formula needed — it falls out of the matrix multiplication rules.

**Properties of the commutator relevant to error propagation:**
- The commutator is **bilinear**: $[\alpha A_1 + \beta A_2, B] = \alpha[A_1, B] + \beta[A_2, B]$
- For fixed $\Omega_T$, the map $R_{\mathrm{pos}} \mapsto \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ is **linear in $R_{\mathrm{pos}}$**
- Therefore: $\frac{1}{2i}[\Omega_T, \delta(R_{\mathrm{pos}})]$ gives the transport error from position error — same operation, different input

### 3.4 Derivation from Time Derivative

If $M_T(t)$ is a time-varying rotation, the velocity of a rotated position is obtained by the product rule:

$$\frac{d}{dt}\left(M_T R_{\mathrm{pos}} M_T^\dagger\right) = \dot{M}_T R_{\mathrm{pos}} M_T^\dagger + M_T R_{\mathrm{pos}} \dot{M}_T^\dagger$$

For a unitary matrix rotating at angular velocity $\boldsymbol{\omega}_T$, the time derivative satisfies $\dot{M}_T = \frac{i}{2}\Omega_T' M_T$ where $\Omega_T' = M_T \Omega_T M_T^\dagger$ is the angular velocity expressed in the target frame, and correspondingly $\dot{M}_T^\dagger = -M_T^\dagger \frac{i}{2}\Omega_T'$. Substituting into the product rule and writing $R_{\mathrm{pos}}' = M_T R_{\mathrm{pos}} M_T^\dagger$:

$$\frac{d}{dt}\left(M_T R_{\mathrm{pos}} M_T^\dagger\right) = \frac{i}{2}\Omega_T' R_{\mathrm{pos}}' - R_{\mathrm{pos}}' \frac{i}{2}\Omega_T' = \frac{i}{2}[\Omega_T', R_{\mathrm{pos}}']$$

Now connect to the cross product convention. From §1.5, the cross product is $\boldsymbol{\omega} \times \mathbf{r} \leftrightarrow \frac{1}{2i}[\Omega, R_{\mathrm{pos}}]$. Since $\frac{i}{2} = -\frac{1}{2i}$ (multiply numerator and denominator by $i$: $\frac{i}{2} = \frac{i \cdot i}{2 \cdot i} = \frac{-1}{2i}$), we have:

$$\frac{d}{dt}\left(M_T R_{\mathrm{pos}} M_T^\dagger\right) = -\frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}'] = \frac{1}{2i}[R_{\mathrm{pos}}', \Omega_T']$$

The antisymmetry of the commutator ($[A,B] = -[B,A]$) means $\frac{1}{2i}[R_{\mathrm{pos}}', \Omega_T'] = -\frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}']$, which corresponds to the cross product $\mathbf{r}' \times \boldsymbol{\omega}_T = -\boldsymbol{\omega}_T \times \mathbf{r}'$. The transport velocity in the right-hand convention ($\boldsymbol{\omega}_T \times \mathbf{r}'$) is therefore:

$$V_{\mathrm{transport}} = \frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}']$$

which matches §3.3 exactly. **The transport theorem IS the product rule of matrix differentiation.** No separate physical argument needed.

### 3.5 Instantaneous Angular Velocity Assumption

The derivation in §3.4 treats $\Omega_T$ as constant during the differentiation. For time-varying $\Omega_T(t)$, additional terms appear (angular acceleration). For SGP4, $\omega_E$ varies by $< 10^{-8}$ rad/s over any propagation step, so the constant-$\Omega_T$ assumption is excellent.

---

## Part 4: The SGP4 Pipeline

### 4.1 Transform 1: Perifocal → TEME

Pure rotation, no translation, no frame angular velocity:

$$M_T = M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u)$$

Since $M_3(\alpha)$ is diagonal ($e^{\pm i\alpha/2}$ on the diagonal), the first and third multiplications are just phase rotations — computationally cheap.

$$R_{\mathrm{pos,TEME}} = M_T\, R_{\mathrm{pos,PF}}\, M_T^\dagger$$
$$V_{\mathrm{TEME}} = M_T\, V_{\mathrm{PF}}\, M_T^\dagger$$
$$\Omega_T = 0, \quad T_T = 0$$

**Velocity transforms identically to position** — same conjugation, no transport term, no additional terms. In this picture there is no asymmetry whatsoever between position and velocity.

### 4.2 Transform 2: TEME → PEF

Rotation by GMST about $\hat{\mathbf{e}}_3$, with frame angular velocity $\boldsymbol{\omega}_E = \omega_E \hat{\mathbf{e}}_3$:

$$M_T = M_3(\theta_{\mathrm{GMST}}) = \begin{pmatrix} e^{i\theta/2} & 0 \\ 0 & e^{-i\theta/2} \end{pmatrix}$$

$$\Omega_T = \omega_E \sigma_3 = \begin{pmatrix} \omega_E & 0 \\ 0 & -\omega_E \end{pmatrix}$$

The velocity transform includes the transport term:

$$V_{\mathrm{PEF}} = M_T\, V_{\mathrm{TEME}}\, M_T^\dagger + \frac{1}{2i}[\Omega_T, R_{\mathrm{pos,PEF}}]$$

Since both $M_T$ and $\Omega_T$ are diagonal, the transport commutator has a particularly simple form. For $R_{\mathrm{pos,PEF}} = \begin{pmatrix} r_z & r_- \\ r_+ & -r_z \end{pmatrix}$ where $r_\pm = r_x \pm ir_y$:

$$\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}] = \frac{\omega_E}{2i}\begin{pmatrix} 0 & 2r_- \\ -2r_+ & 0 \end{pmatrix} = \omega_E \begin{pmatrix} 0 & -ir_- \\ ir_+ & 0 \end{pmatrix}$$

Extracting the vector: the transport velocity is $(-\omega_E r_y, \omega_E r_x, 0)$, which is indeed $\boldsymbol{\omega}_E \times \mathbf{r} = (-\omega_E r_y, \omega_E r_x, 0)$. ✓

### 4.3 Composition: Perifocal → PEF

$$M_{\mathrm{PF \to PEF}} = M_3(\theta) \cdot M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u)$$

Since $M_3$ is diagonal, the first two compose trivially:

$$M_3(\theta) \cdot M_3(-\Omega) = M_3(\theta - \Omega)$$

This is a computational simplification: two diagonal multiplications reduce to one.

The velocity transform accumulates the transport term from Transform 2:

$$V_{\mathrm{PEF}} = M_{\mathrm{total}}\, V_{\mathrm{PF}}\, M_{\mathrm{total}}^\dagger + \frac{1}{2i}[\Omega_E, R_{\mathrm{pos,PEF}}]$$

The transport term acts on the already-rotated position.

### 4.4 Near-Identity Perturbations

A small rotation: $M = I + \frac{i}{2}\varepsilon_k\sigma_k + O(\varepsilon^2)$ where $\varepsilon_k \ll 1$.

Composition of two small rotations:

$$(I + \frac{i}{2}A)(I + \frac{i}{2}B) = I + \frac{i}{2}(A + B) + O(\varepsilon^2)$$

**Perturbations compose additively to first order.** This is the basis for treating short-period and long-period corrections as additive.

The transport term for a near-identity perturbation $\delta M = \frac{i}{2}\varepsilon_k\sigma_k$ applied to position $R_{\mathrm{pos}}$:

$$\delta V = \frac{1}{2i}[\delta\Omega, R_{\mathrm{pos}}] \approx \frac{1}{2i}\varepsilon_k[\sigma_k, R_{\mathrm{pos}}]$$

This is linear in $\varepsilon_k$ — the perturbation is a linear correction to velocity, proportional to the perturbation magnitude.

---

## Part 5: Error Calculus

### 5.1 The Parallel-Row Principle — Unified

Since all vectors transform by $V' = M V M^\dagger$ (conjugation), and conjugation is **linear in $V$** for fixed $M$:

$$\delta(V') = M\,\delta(V)\, M^\dagger$$

**The error transforms by the same conjugation as the physical vector.** This applies uniformly to:
- Position error: $\delta(R_{\mathrm{pos}}') = M_T\,\delta(R_{\mathrm{pos}})\, M_T^\dagger$
- Velocity error: $\delta(V') = M_T\,\delta(V)\, M_T^\dagger$
- Angular velocity error: $\delta(\Omega') = M_T\,\delta(\Omega)\, M_T^\dagger$

No asymmetry. Every vector error is a traceless Hermitian 2×2 matrix, transformed by the same operation.

For the transport term (a commutator), linearity gives:

$$\delta\!\left(\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']\right) = \frac{1}{2i}[\Omega_T, \delta(R_{\mathrm{pos}}')]$$

Same commutator, different input. The transport error IS the transport operation applied to the position error.

### 5.2 The Four-Row Structure

$$\begin{array}{c|cccc}
 & \text{Orientation} & \text{Position} & \text{Velocity} & \text{Ang. vel.} \\
 & M \in SU(2) & V_r: \text{tl-Herm} & V_v: \text{tl-Herm} & V_\omega: \text{tl-Herm} \\
\hline
\text{Physical} & M & R_{\mathrm{pos}} & V & \Omega \\
\sigma_m & \delta_m(M) & \delta_m(R_{\mathrm{pos}}) & \delta_m(V) & \delta_m(\Omega) \\
\delta_p & \delta_p(M) & \delta_p(R_{\mathrm{pos}}) & \delta_p(V) & \delta_p(\Omega) \\
\delta_a & \delta_a(M) & \delta_a(R_{\mathrm{pos}}) & \delta_a(V) & \delta_a(\Omega) \\
\end{array}$$

Each error entry is a 2×2 traceless Hermitian matrix (for vectors) or a tangent vector to SU(2) at $M$ (for orientation). The same transform applies to all rows.

**Storage per row:** 4 (orientation) + 3 (position) + 3 (velocity) + 3 (angular velocity) = 13 real numbers.

**Total for all four rows:** 4 × 13 = 52 real numbers. (Physical state + 3 error categories.)

### 5.3 The Isometry Property

For exact $M_T$ (no error in the rotation):

$$\|\mathbf{v}'\|^2 = -\det(V') = -\det(M_T V M_T^\dagger) = -\det(M_T)\det(V)\det(M_T^\dagger) = -\det(V) = \|\mathbf{v}\|^2$$

since $\det(M_T) = 1$ for $M_T \in SU(2)$. Exact rotations preserve the norm. Applied to the error rows:

$$\|\delta(\mathbf{v}')\| = \|\delta(\mathbf{v})\|$$

**An exact rotation does not amplify errors.** This is immediate from $\det(M) = 1$.

### 5.4 Error from Inexact Rotation

When $M_T$ carries error (from computing sin/cos of orbital angles):

$$V' = (M_T + \delta M_T)\, V\, (M_T + \delta M_T)^\dagger$$

Expanding to first order in $\delta M_T$:

$$V' \approx M_T V M_T^\dagger + \delta M_T\, V\, M_T^\dagger + M_T\, V\, \delta M_T^\dagger$$

The additional error has magnitude $\approx 2\|\mathbf{v}\|\,\|\delta M_T\|$. The factor of 2 comes from the two cross terms (left and right multiplication).

For angles known to $d$ decimal digits: $\|\delta M_T\| \sim 10^{-d}$, so the rotation-induced position error is $\sim \|\mathbf{r}\| \cdot 10^{-d}$.

### 5.5 Error Through a Composition Chain

The 3-1-3 rotation $M = M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u)$ is three matrix multiplications. Each introduces error from:
1. The input angle uncertainty ($\sigma_m$)
2. The sin/cos evaluation ($\delta_p$)
3. The complex arithmetic rounding ($\delta_p$)

Since each $M_k$ is unitary, $\|M_k\| = 1$. The error accumulates additively through the chain:

$$\|\delta M\| \leq \|\delta M_3(-\Omega)\| + \|\delta M_1(-i)\| + \|\delta M_3(-u)\| + \text{rounding}$$

(Rigorous bound from the submultiplicativity of the operator norm and the unitarity of each factor.)

### 5.6 Nonlinear Operations

**Constructing $M$ from angle $\alpha$:**

$$M_k(\alpha) = \cos\frac{\alpha}{2}\,I + i\sin\frac{\alpha}{2}\,\sigma_k$$

This is nonlinear in $\alpha$. The sensitivity:

$$\frac{dM_k}{d\alpha} = -\frac{1}{2}\sin\frac{\alpha}{2}\,I + \frac{i}{2}\cos\frac{\alpha}{2}\,\sigma_k$$

The norm: $\left\|\frac{dM_k}{d\alpha}\right\| = \frac{1}{2}$ (constant, independent of $\alpha$). So:

$$\|\delta M_k\| \leq \frac{1}{2}|\delta\alpha|$$

This is a **rigorous bound** (not an approximation): the derivative has constant norm $\frac{1}{2}$, so the first-order bound is tight.

**Renormalization $M/\sqrt{\det M}$:** When $M$ drifts from unitarity after accumulated rounding, renormalization divides by $\sqrt{\det M}$. For near-unitary $M$ with $\det M = 1 + \varepsilon$:

$$\frac{M}{\sqrt{\det M}} \approx M(1 - \varepsilon/2) = M - \frac{\varepsilon}{2}M$$

The correction is proportional to the determinant drift $\varepsilon$, which is $O(n\epsilon_{\mathrm{mach}})$ after $n$ multiplications.

**Euclidean norm** $\|\mathbf{v}\| = \sqrt{-\det V}$: Nonlinear (square root). Sensitivity from Ch 1, Corollary 1.4.4. The norm never amplifies total error magnitude (§5.3).

### 5.7 SGP4 Error Source Map

| Source | Error type | Enters as | Propagates via |
|--------|-----------|-----------|---------------|
| $\Omega$ from TLE | $\sigma_m$ | $M_3(-\Omega)$: $\|\delta M\| \leq \frac{1}{2}\delta\Omega$ | Conjugation of all vectors |
| $i$ from TLE | $\sigma_m$ | $M_1(-i)$: $\|\delta M\| \leq \frac{1}{2}\delta i$ | Same |
| $u$ from TLE | $\sigma_m$ | $M_3(-u)$: $\|\delta M\| \leq \frac{1}{2}\delta u$ | Same |
| sin/cos evaluation | $\delta_p$ | $M_k$ matrix entries | $\leq \epsilon_{\mathrm{mach}}$ per entry |
| Complex arithmetic rounding | $\delta_p$ | Each 2×2 multiply | $\sim 4\epsilon_{\mathrm{mach}}$ per multiply |
| $\theta_{\mathrm{GMST}}$ | $\delta_p$ | $M_3(\theta)$ for TEME→PEF | Same as angle sensitivity |
| GMST coefficients | $\sigma_m$ | $\theta_{\mathrm{GMST}}$ value | Through $M_3(\theta)$ |
| $\omega_E$ | $\sigma_m$ | $\Omega_T = \omega_E\sigma_3$ | Commutator $[\Omega_T, R_{\mathrm{pos}}]$ |
| TEME frame | $\delta_a$ | Entire frame definition | Irreducible floor on output |

### 5.8 Linearization Validity

The parallel-row principle is exact for linear operations (conjugation, commutator). For nonlinear operations (angle→matrix, renormalization), it uses the first-order sensitivity. This is valid when:

$$\mathrm{rd}(v) \geq d_{\min}$$

where $d_{\min}$ is a threshold (typically 4) and $\mathrm{rd}$ is the reliable-digits count from Ch 1, Definition 1.9.1. When $\mathrm{rd} < d_{\min}$, the full nonlinear error formulas of Ch 1 must be used.

### 5.9 Scalar Reduction

**Position error:**
$$\delta(\text{position}) = \sqrt{-\det(\delta R_{\mathrm{pos}})} = \sqrt{\delta r_x^2 + \delta r_y^2 + \delta r_z^2}$$

**Velocity error:**
$$\delta(\text{velocity}) = \sqrt{-\det(\delta V)} = \sqrt{\delta v_x^2 + \delta v_y^2 + \delta v_z^2}$$

**Orientation error:** The angular error for a near-identity perturbation $\delta M = \frac{i}{2}\varepsilon_k\sigma_k$ is:
$$\delta(\text{orientation}) = \sqrt{\varepsilon_1^2 + \varepsilon_2^2 + \varepsilon_3^2}$$
This is the magnitude of the rotation perturbation vector. For the half-angle encoding, $\varepsilon_k$ are the small-angle perturbations, so this directly gives the angular error in radians.

**Angular velocity error:**
$$\delta(\text{angular velocity}) = \sqrt{-\det(\delta\Omega)} = \sqrt{\delta\omega_x^2 + \delta\omega_y^2 + \delta\omega_z^2}$$

All four reductions use the same algebraic operation: $\sqrt{-\det(\cdot)}$ of a traceless Hermitian matrix. This is the norm in the algebra.

### 5.10 Merged vs. Separate Error Modes

**Separate (4 rows):** Physical + $\sigma_m$ + $\delta_p$ + $\delta_a$. Same transform applied 4 times. Diagnostic: which category dominates.

**Merged (2 rows):** Physical + $\delta_{\mathrm{total}}$ where $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$ per component. Same transform applied twice. Cheaper but loses category breakdown.

### 5.11 Reliable-Digits Connection

The scalar reductions feed into Ch 1, Definition 1.9.1:

$$\mathrm{rd}(\text{position}) = \left\lfloor -\log_{10}\!\left(\frac{\delta(\text{position})}{\|\mathbf{r}\|}\right) \right\rfloor$$

Decision criteria:
- $\mathrm{rd}(\sin i) < 2$: switch to alternative angle parameterization
- $\mathrm{rd}(E) < d_{\min}$: Kepler solver needs more iterations
- $\mathrm{rd}(\text{position}) < 0$: output is meaningless

---

## Part 6: Affine Operations via Linear Multiplication

### 6.1 The Core Principle in the SU(2) Framework

The rotation action $V' = M V M^\dagger$ is **linear in $V$** for fixed $M$. This is immediately apparent: conjugation is a linear map on the vector space of 2×2 Hermitian matrices.

The **translation** $V' = M V M^\dagger + T$ is **affine** in $V$ — a linear map plus a constant.

The **transport term** $\frac{1}{2i}[\Omega_T, V']$ is **linear in $V'$** for fixed $\Omega_T$ (the commutator is linear in each argument).

### 6.2 Embedding Affine Operations as Linear

To absorb the translation into a linear multiplication, we extend the algebra. The dual SU(2) matrix $\hat{M} = M + \varepsilon D$ achieves this:

$$\hat{V}' = \hat{M}\, \hat{V}\, \hat{M}^\dagger$$

The real part gives the rotation; the dual part gives the translation. A single algebraic product — linear in $\hat{V}$ for fixed $\hat{M}$ — replaces "rotate then add."

### 6.3 Why Linearity Enables the Parallel-Row Principle

For any linear map $\mathcal{L}$:

$$\mathcal{L}(\mathbf{x} + \delta\mathbf{x}) = \mathcal{L}(\mathbf{x}) + \mathcal{L}(\delta\mathbf{x})$$

So:

$$\delta(\mathcal{L}(\mathbf{x})) = \mathcal{L}(\delta\mathbf{x})$$

The same linear map $\mathcal{L}$ applied to the error gives the error of the output. This is exact — not a first-order approximation. No Taylor expansion, no neglected terms.

For **affine** operations $\mathcal{A}(\mathbf{x}) = \mathcal{L}(\mathbf{x}) + \mathbf{c}$:

$$\delta(\mathcal{A}(\mathbf{x})) = \mathcal{L}(\delta\mathbf{x})$$

The constant drops out. The error propagation is the same as for the linear part.

### 6.4 What Linearity Covers in the SGP4 Pipeline

| Operation | Expression | Linear in state? | Consequence |
|-----------|-----------|-----------------|-------------|
| Rotation of any vector | $M V M^\dagger$ | Yes (in $V$) | Error row: same operation |
| Translation | $+ T$ | Affine (const) | Error row: zero contribution |
| Transport term | $\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ | Yes (in $R_{\mathrm{pos}}$) | Error row: same commutator |
| Angular velocity offset | $+ \Omega_T$ | Affine (const) | Error row: zero contribution |
| Composition of rotations | $M_{AB} \cdot M_{BC}$ | Yes (in $M_{BC}$) | Error row: same multiply |
| Near-identity perturbation | $I + \frac{i}{2}\delta\sigma$ | Linear (in $\delta$) | Error row: same addition |
| Secular rate addition | $+ \dot{M}\Delta t$ | Affine (const) | Error row: from $\delta(\dot{M})\Delta t$ |

**The only nonlinear steps are constructing the transform itself** (computing $M$ from angles, renormalization). The APPLICATION of the transform to the state is always linear.

---

## Part 7: The Matrix Form as a Derived Consequence

### 7.1 The 3×3 Rotation Matrix from SU(2)

The map $V \mapsto M V M^\dagger$ is a linear map on the 3D real vector space of traceless Hermitian matrices. In the Pauli basis, this map has a 3×3 real matrix representation:

$$R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j\, M\, \sigma_k\, M^\dagger)$$

For $M = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix}$ with $\alpha = w + iz$, $\beta = ix + y$, this produces exactly the standard rotation matrix:

$$R(q) = \begin{pmatrix} 1 - 2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\ 2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\ 2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2) \end{pmatrix}$$

The matrix IS the conjugation map written in coordinates.

### 7.2 The 7×7 State Matrix from the Full Transform

Writing the full state transform (§3.2) in Cartesian components:

- Position: $\mathbf{r}' = R\,\mathbf{r} + \mathbf{t}$
- Velocity: $\mathbf{v}' = R\,\mathbf{v} + [\boldsymbol{\omega}_T]_\times R\,\mathbf{r} + \Delta\mathbf{v}$

where $[\boldsymbol{\omega}_T]_\times$ is the skew-symmetric matrix corresponding to the commutator $\frac{1}{2i}[\Omega_T, \cdot]$ written in coordinates.

In block form:

$$\begin{pmatrix} \mathbf{r}' \\ \mathbf{v}' \\ 1 \end{pmatrix} = \begin{pmatrix} R & \mathbf{0} & \mathbf{t} \\ [\boldsymbol{\omega}]_\times R & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \begin{pmatrix} \mathbf{r} \\ \mathbf{v} \\ 1 \end{pmatrix}$$

This is exactly the 7×7 state matrix from the archived DCM chapter. **The matrix form is the SU(2) conjugation written in coordinates.**

---

## Part 8: Section-by-Section Requirements Map

| Template Section | Required from this reference architecture |
|-----------------|------------------------------------------|
| **§2.1** (Motivation) | §6.3: linearity as the unifying principle |
| **§2.2** (Rotation) | Part 1: Pauli matrices, SU(2), conjugation action, elementary rotations, Shepperd extraction |
| **§2.3** (Translation) | Part 2: dual extension, encoding/extraction, composition |
| **§2.4** (Rigid body state) | §3.1: uniform representation; §3.2: state transform; §3.3–3.4: transport as commutator/derivative |
| **§2.5** (SGP4 pipeline) | Part 4: concrete transforms, near-identity; §5.7: error source map |
| **§2.6** (Error propagation) | Part 5: parallel-row principle, isometry, inexact rotation, composition chain, linearization validity |
| **§2.7** (Matrix form) | Part 7: $R(q)$ derivation, 7×7 from SU(2) |
| **§2.8** (Parallel-row architecture) | §5.2: four-row structure; §5.9: scalar reduction; §5.10: merged vs separate; §5.11: reliable-digits |
| **§2.9** (Summary) | This table |

---

## Part 9: Source Code Audit — Gap Remediation

The following subsections address gaps identified by auditing 206 operations across 11 source files. Each gap is stated in terms of the **problem** and the **requirements** the framework must satisfy — without committing to specific formulas or numerical values, which may change when expressed in the SU(2) algebra.

### 9.1 Chained Inverse-Power Operations

**Problem:** Perturbation computations frequently involve chained reciprocal powers ($p^{-2}$, $p^{-4}$, $\xi^4$). Each reciprocal amplifies error, and chaining compounds the amplification. These appear throughout the secular rate, drag, and element recovery computations.

**What the framework must provide:**
- A general bound on error amplification through $n$ chained divisions or reciprocal powers, expressed in terms of the operand's magnitude and error
- Identification that these are scalar (non-rotational) operations — they live outside the conjugation framework and require Ch 1 division bounds (Theorem 1.3.3) applied iteratively
- Classification as nonlinear operations in the pipeline — the error row cannot reuse the same computation; it must use the sensitivity of each step

**What the Build phase must determine:**
- Whether the specific amplification factors derived from the traditional formulation still hold when the perturbation coefficients are re-expressed in the SU(2) framework
- Whether alternative formulations (e.g., computing $p^{-4}$ directly rather than as $(p^{-2})^2$) reduce the error amplification
- The orbit regimes where the amplification becomes large enough to degrade the output

**Code locations:** `brouwer.h:102-103`, `element_recovery.h:83,105,115`, `drag_coefficients.h:121`.

### 9.2 Subtractive Cancellation in Scalar Denominators

**Problem:** The drag coefficient computation involves a denominator that can approach zero when orbital parameters reach certain physical regimes (perigee near the density scale height). This is a scalar singularity — it does not involve the rotation representation — but it appears in the pipeline and affects the error budget.

**What the framework must provide:**
- Recognition that certain scalar quantities in the perturbation models are susceptible to subtractive cancellation (Ch 1, §1.7)
- A detection criterion based on reliable digits (Ch 1, Definition 1.9.1) rather than hardcoded thresholds
- A placeholder for alternative formulations that avoid the cancellation (to be developed in the chapters where each perturbation model is derived)

**What the Build phase must determine:**
- The exact conditions under which each scalar singularity is triggered in the SU(2) formulation (these may differ from the traditional formulation if the perturbation expressions are restructured)
- Whether refactoring the expression (e.g., factoring a difference of squares) eliminates the cancellation
- Rigorous bounds on the sensitivity amplification at each singularity

**Code locations:** `drag_coefficients.h:118,122,155,196`.

### 9.3 Large-Angle Accumulation and the SU(2) Alternative

**Problem:** The traditional pipeline accumulates orbital angles over time, producing large values (tens of thousands of radians for long propagations). Modular reduction (wrap to $[0, 2\pi)$) introduces precision loss from the reduction itself and from evaluating trigonometric functions at large arguments.

**What the framework must provide:**
- An explanation of how the SU(2) representation avoids this problem: instead of accumulating angles and then evaluating trig functions, the rotation itself is accumulated by composing small-angle rotations directly
- Analysis of the tradeoffs: accumulating rotations introduces norm drift (requiring periodic renormalization), while accumulating angles introduces large-argument trig loss. The framework should characterize both error modes so the Build phase can compare them
- A clear statement of the assumed operating mode: does the SU(2) propagator accumulate rotations incrementally, or does it compute the rotation from the accumulated angle at each step?

**What the Build phase must determine:**
- The norm drift rate for accumulated SU(2) rotations (how many compositions before renormalization is needed)
- The precision advantage (if any) of incremental rotation accumulation over angle accumulation for typical SGP4 propagation durations
- Whether the matched-pair principle (Ch 3) constrains which accumulation method is permitted in standard mode

**Code locations:** `secular_update.h:138-141`, `near_space.h:307`, `deep_space.h:299`.

### 9.4 Conditional Guards and Singularity Avoidance

**Problem:** The existing code uses hardcoded numerical thresholds to guard against singularities (eccentricity near zero, inclination near zero or $\pi$). These thresholds are not derived from the error analysis — they are empirical values that may not be optimal.

**What the framework must provide:**
- The principle that singularity detection should use the reliable-digits criterion (Ch 1, Definitions 1.9.1 and 1.9.3) rather than hardcoded thresholds
- Identification of which singularities are eliminated by the SU(2) representation (the Euler angle singularities at $i = 0$ and $e = 0$ for angular quantities) and which remain (scalar singularities in the perturbation models)
- For the remaining singularities: a framework for expressing the guard condition in terms of the error state, so that the threshold adapts to the actual precision of the computation

**What the Build phase must determine:**
- Whether the SU(2) representation truly eliminates the eccentricity and inclination guards, or whether those singularities resurface in a different form (e.g., in the perturbation coefficients that feed into the rotation)
- What the reliable-digits thresholds should be for each guard — this requires the specific error bounds from the perturbation chapters (Ch 16–19, 21–22)
- Whether eliminating a guard changes the matched-pair relationship with the standard SGP4 model

**Code locations:** `secular_update.h:120`, `deep_space.h:278`, `near_space.h:237`, `drag_coefficients.h:155,196`.

### 9.5 Time-Polynomial Error Growth

**Problem:** Secular updates advance orbital elements as polynomials in propagation time (up to $t^5$). The error in these polynomials grows with time because each coefficient carries uncertainty that is multiplied by increasing powers of $t$.

**What the framework must provide:**
- Recognition that polynomial-in-time evaluation is a scalar operation happening before the frame transform — it produces updated orbital elements that are then fed into the SU(2) rotation construction
- The observation that polynomial evaluation is linear in the coefficients (for fixed $t$), so the error rows transform by the same polynomial evaluation — the parallel-row principle applies to this stage
- A characterization of the error growth regime: which power of $t$ dominates the error at typical propagation durations

**What the Build phase must determine:**
- Whether the polynomial coefficients themselves change when re-derived in the SU(2) framework (they should not — the secular rates are scalar quantities independent of the rotation representation)
- The reliable-digits count on each orbital element as a function of propagation time — this is the fundamental limit on how far an analytical propagator can predict
- Whether incremental propagation (re-initializing the polynomial at intermediate epochs) reduces the time-polynomial error growth, and whether this is compatible with the matched-pair principle

**Code locations:** `secular_update.h:77-145`.

### 9.6 Unverified Coefficients

**Problem:** At least one numerical coefficient in the existing code (`134/81` in the element recovery third-order series inversion) lacks a derivation.

**What the framework must provide:**
- The principle (already established in the project rules) that every coefficient must be derived, not matched to reference code
- Identification of all such unverified coefficients across the codebase
- A tracking mechanism so the Build phase can confirm each derivation is complete before the coefficient is used

**What the Build phase must determine:**
- Whether the coefficient is correct (derive it from first principles in the appropriate chapter)
- Whether the SU(2) formulation changes the expression in which the coefficient appears, potentially altering or eliminating it

**Code locations:** `element_recovery.h:110`.

---

## Resolved Design Decisions

**State representation:** 13 real parameters per row (4 for SU(2) rotation + 3×3 for position, velocity, angular velocity as traceless Hermitian). Four rows (physical + 3 error categories) = 52 real numbers total.

**14 vs 15 components:** With the SU(2) framework, the state is 13 real parameters (not 14). The "14th component" from v1 was the homogeneous coordinate in the dual quaternion. Here, the dual extension handles translation algebraically, and the angular velocity offset $\Omega_T$ is a known constant handled as a post-addition. No homogeneous coordinate needed.

**Affine terms:** Translation and angular velocity offsets are known constants. They drop out of the error rows. Handled as additions after the linear conjugation.

## Implementation Items — Dependency Order and Resolutions

### Dependency Graph

```
FOUNDATIONAL (no upstream deps):        RESOLVED BELOW
  #4  Complex vs real arithmetic         → §R.1
  #13 Polynomial coefficients unchanged  → §R.2
  #19 Template editorial update          → DONE (committed b793b34)

SECOND TIER:
  #1  Orientation error geometry          → §R.3 (depends on understanding #4)
  #3/#7 Renormalization / norm drift      → §R.4 (depends on #4)

THIRD TIER:
  #2  Position extraction error           → OPEN (needs #1)
  #8  Rotation vs angle accumulation      → OPEN (needs #3/#7)

FOURTH TIER:
  #18 d_min threshold                     → OPEN (needs #1, #2)
  #10 Guard elimination verification      → OPEN (needs #1, #18)
```

### §R.1 Complex vs. Real Arithmetic (Item #4) — RESOLVED

The SU(2) framework uses 2×2 complex matrices for mathematical clarity. An implementation stores each complex entry as two reals. The SU(2) matrix $M = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix}$ has two independent complex numbers $\alpha, \beta$ (the other two entries are determined by them). Storage: 4 real numbers — identical to the quaternion representation $(w, x, y, z)$ with the identification $\alpha = w + iz$, $\beta = ix + y$.

**Resolution:** The SU(2) framework and the quaternion representation use the same 4 real numbers. The choice of notation (complex 2×2 vs. quaternion 4-tuple) is a mathematical convenience, not an implementation choice. Error propagation through "complex multiplication" is identical to error propagation through the "Hamilton product" — they are the same 16 real multiply-adds (Ch 1, Theorems 1.3.1–1.3.2). The TrackedValue machinery operates on real numbers regardless of which mathematical notation is used.

The reference architecture uses SU(2) notation for mathematical clarity (the linearity of conjugation $V \mapsto M V M^\dagger$ is immediate from the matrix framework). The implementation uses 4 real TrackedValues per rotation, as it does today.

### §R.2 Polynomial Coefficients Unchanged (Item #13) — RESOLVED

The secular rate coefficients (Brouwer's polynomials in $\cos^2 i$, the drag coefficients $C_1$–$C_5$, $D_2$–$D_4$) are scalar quantities computed from orbital elements. They do not depend on the rotation representation — they are the same whether the frame transform uses Euler angles, DCMs, quaternions, or SU(2) matrices.

**Resolution:** The secular update is a scalar polynomial evaluated BEFORE the frame transform. The polynomial coefficients are derived from the perturbation theory (Ch 16–22) and from the TLE elements. The SU(2) framework changes only how the resulting updated elements are converted into a frame rotation — it does not change the elements themselves or the perturbation physics that updates them.

Verification: the inputs to `secular_update.h` are orbital elements $(a, e, i, \Omega, \omega, M)$ and the outputs are updated orbital elements. No rotation representation appears in this computation. The SU(2) framework enters only after the secular update, when `state_from_elements.h` constructs the rotation from the updated angles.

### §R.3 Orientation Error Geometry (Item #1) — RESOLVED

A small perturbation of an SU(2) rotation $M$ is:

$$M + \delta M = M(I + \frac{i}{2}\varepsilon_k\sigma_k) = M \cdot M_{\mathrm{pert}}$$

where $M_{\mathrm{pert}} = I + \frac{i}{2}\varepsilon_k\sigma_k + O(\varepsilon^2)$ is a near-identity rotation. The perturbation parameters $\varepsilon_k$ are the components of the rotation error vector.

From §1.9, a rotation by angle $\alpha$ about axis $\hat{\mathbf{e}}$ gives $M = \cos(\alpha/2)\,I + i\sin(\alpha/2)\,\hat{e}_k\sigma_k$. For a small rotation ($\alpha \ll 1$): $M \approx I + i(\alpha/2)\hat{e}_k\sigma_k$. Comparing with the perturbation form: $\varepsilon_k = \alpha\,\hat{e}_k$, so $\|\boldsymbol{\varepsilon}\| = |\alpha|$.

**Resolution:** The scalar orientation error $\delta(\text{orientation}) = \|\boldsymbol{\varepsilon}\| = \sqrt{\varepsilon_1^2 + \varepsilon_2^2 + \varepsilon_3^2}$ gives the rotation angle error **directly in radians** — not half-radians. The factor of $1/2$ in the parameterization $\frac{i}{2}\varepsilon_k\sigma_k$ is absorbed by the half-angle encoding: the perturbation parameters $\varepsilon_k$ already represent the full rotation angle, not the half-angle.

Verification: a rotation by $\alpha = 0.01$ rad about $\hat{\mathbf{e}}_3$ gives $M_3(0.01) = I + i(0.005)\sigma_3 + O(10^{-4})$, so $\varepsilon_3 = 0.01$ (the full angle). ✓

### §R.4 Renormalization and Norm Drift (Items #3, #7) — RESOLVED

For exact arithmetic, the product of two SU(2) matrices is SU(2) ($\det(M_1 M_2) = \det M_1 \cdot \det M_2 = 1$). In finite precision, each multiplication introduces rounding error that causes $\det M$ to drift from 1.

Each 2×2 complex matrix multiply involves 8 complex multiply-adds = 16 real multiply-adds. Each introduces relative error $\leq \epsilon_{\mathrm{mach}}$. The determinant $\det M = |\alpha|^2 + |\beta|^2$ accumulates error at rate:

$$|\det M_n - 1| \leq n \cdot c \cdot \epsilon_{\mathrm{mach}}$$

where $n$ is the number of multiplications and $c$ is a small constant depending on the specific arithmetic (typically $c \leq 8$ for the 2×2 case).

For binary64 ($\epsilon_{\mathrm{mach}} \approx 1.1 \times 10^{-16}$):
- After 100 multiplications: drift $\leq 10^{-13}$ — negligible
- After 10,000 multiplications: drift $\leq 10^{-11}$ — still negligible
- After $10^6$ multiplications: drift $\leq 10^{-9}$ — may warrant renormalization

**Resolution:** For the SGP4 pipeline, the total number of SU(2) multiplications per propagation step is small (3 for the 3-1-3 Euler composition + 1 for the GMST rotation = 4 per step). Even for a year of propagation at 1-minute steps ($\sim 5 \times 10^5$ steps), the accumulated drift is $\leq 4 \times 5 \times 10^5 \times 8 \times 10^{-16} \approx 10^{-9}$. This is far below the measurement error floor ($\sigma_m \sim 10^{-8}$ from TLE precision).

**Renormalization is not needed for SGP4-class propagation durations.** If the incremental rotation accumulation approach (§9.3) is adopted — where rotations are composed across time steps rather than recomputed from accumulated angles — then renormalization should be applied when $\mathrm{rd}(\det M) < d_{\min}$, using the reliable-digits criterion. The renormalization itself ($M \leftarrow M / \sqrt{\det M}$) has sensitivity $\leq 1 + O(\det M - 1)$, which is near-unity when the drift is small.

## Remaining Open Items

### Third Tier (depend on resolved items above)

**#2 Position extraction error.** Extracting $\mathbf{r}$ from the dual part: $T = 2 D M^\dagger$. This is a product of $D$ (carrying its own error) and $M^\dagger$ (carrying orientation error from §R.3). By Ch 1, Theorem 1.3.2 (multiplication): $\delta(T) \leq |D| \cdot \delta(M^\dagger) + |M^\dagger| \cdot \delta(D)$. Since $|M^\dagger| = 1$ (unitary), the second term simplifies: $\delta(T) \leq |D| \cdot \delta(M) + \delta(D)$. The position error is therefore bounded by the dual-part error plus a coupling from the rotation error scaled by the position magnitude. **Status: framework clear, specific bound to be derived in Build.**

**#8 Rotation vs angle accumulation.** From §R.4, the rotation accumulation drift is $\sim 10^{-9}$ for a year of propagation. The angle accumulation loss (§9.3) is $\sim 4$ digits from large-argument trig reduction at $M \sim 34000$ rad. The rotation accumulation approach is therefore more precise by $\sim 7$ orders of magnitude for long propagations. **Status: comparison clear, but whether this is compatible with the matched pair principle (Ch 3) is deferred.**

### Fourth Tier

**#18 d_min threshold.** The reliable-digits threshold that triggers formula switching or nonlinear error propagation. From the resolved items: orientation error is in radians (§R.3), position extraction error has a known structure (#2 above), and the pipeline's linear steps preserve error magnitude (isometry). The threshold should be set where the linearization error (second-order terms neglected in the parallel-row principle) exceeds the first-order bound by a specified fraction. **Status: requires the specific second-order remainder bounds, which are Build-phase work.**

**#10 Guard elimination.** From §R.3 and §9.4: the SU(2) rotation is singularity-free, so the inclination and eccentricity guards for *angular quantities* are eliminated. The scalar guards (eccentricity floor for drag coefficients, density singularities) remain and should use the reliable-digits criterion with threshold from #18. **Status: principle clear, specific threshold values are Build-phase work.**
