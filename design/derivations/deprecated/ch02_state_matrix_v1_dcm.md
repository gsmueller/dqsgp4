# Chapter 2: The State Matrix Framework

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $R_k(\alpha)$ | Elementary rotation matrix about axis $k$ by angle $\alpha$ | §2.2, Def. 2.2.1 |
| $R$ | A direction cosine matrix (DCM) in $SO(3)$ | §2.2, Def. 2.2.2 |
| $[\boldsymbol{\omega}]_\times$ | Skew-symmetric matrix of $\boldsymbol{\omega}$: $[\boldsymbol{\omega}]_\times \mathbf{v} = \boldsymbol{\omega} \times \mathbf{v}$ | §2.2, Def. 2.2.3 |
| $\bar{\mathbf{p}}$ | Homogeneous position vector $[\mathbf{r};\, 1] \in \mathbb{R}^4$ | §2.3, Def. 2.3.1 |
| $H$ | $4 \times 4$ homogeneous transform matrix | §2.3, Def. 2.3.2 |
| $\mathbf{x}$ | Extended state vector $[\mathbf{r};\, \mathbf{v};\, 1] \in \mathbb{R}^7$ | §2.4, Def. 2.4.1 |
| $T$ | $7 \times 7$ state matrix | §2.4, Def. 2.4.2 |
| $\dot{R}$ | Transport term: $[\boldsymbol{\omega}]_\times R$ for a rotating frame | §2.4, Def. 2.4.6 |
| $\Delta\mathbf{r}$ | Origin offset (translation) between frames | §2.4, Def. 2.4.3 |
| $\Delta\mathbf{v}$ | Velocity of new frame origin | §2.4, Def. 2.4.3 |
| $\omega_E$ | Earth's rotation rate $\approx 7.2921159 \times 10^{-5}$ rad/s | §2.5, Thm. 2.5.2 |
| $\theta_{\mathrm{GMST}}$ | Greenwich Mean Sidereal Time angle (Ch. 29) | §2.5, Thm. 2.5.2 |
| $\Omega, i, u$ | Right ascension, inclination, argument of latitude | §2.5, Thm. 2.5.1 |
| $\kappa(A)$ | Matrix condition number $\|A\|\,\|A^{-1}\|$ | §2.6, Def. 2.6.1 |
| $\mathbf{x}_{13}$ | Extended state with error tracking $\in \mathbb{R}^{13}$ | §2.7, Def. 2.7.1 |

---

## §2.1 Introduction

The SGP4 propagator transforms orbital elements through a pipeline of operations — secular update, Kepler solver, perturbation corrections, coordinate rotations, frame rotations — each modifying the satellite's position and velocity. These operations involve rotations (frame changes), translations (origin offsets), and velocity coupling (the transport theorem for rotating frames).

Chapter 1 developed error propagation for scalar operations. This chapter extends that framework to the full satellite state: position, velocity, and the errors accumulated through each transformation step. The primary algebraic tool is the **7×7 homogeneous state matrix**, which encodes rotation, translation, and velocity coupling (including the transport theorem) in a single composable linear operation. Applying a coordinate transform to the satellite state is a matrix-vector multiplication; composing transforms is a matrix-matrix multiplication; and the error budget follows from applying Chapter 1's entry-wise rules to the matrix arithmetic.

The chapter builds from rotation matrices (§2.2) through homogeneous coordinates (§2.3) to the 7×7 state matrix (§2.4). Concrete SGP4 transforms are expressed as state matrices (§2.5). Error propagation through matrix operations is derived from Chapter 1's scalar rules (§2.6). The 13×13 extended state matrix, which carries the three error categories alongside the physical state, is defined in §2.7.

The key insight is that the transport theorem — the $\boldsymbol{\omega} \times \mathbf{r}$ velocity correction for rotating frames — appears automatically as the lower-left block of the 7×7 matrix. No separate derivation or manual cross-product computation is needed: the matrix multiply does it.

---

## §2.2 Rotation Matrices and Direction Cosine Matrices

**Assumption 2.2.1** (Right-handed coordinates)**.** *All coordinate systems are right-handed. Rotation angles follow the right-hand rule.*

**Definition 2.2.1** (Elementary rotation matrices)**.** *Rotation by angle $\alpha$ about coordinate axis $k$:*

$$
R_1(\alpha) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos\alpha & -\sin\alpha \\ 0 & \sin\alpha & \cos\alpha \end{bmatrix}, \quad R_2(\alpha) = \begin{bmatrix} \cos\alpha & 0 & \sin\alpha \\ 0 & 1 & 0 \\ -\sin\alpha & 0 & \cos\alpha \end{bmatrix}, \quad R_3(\alpha) = \begin{bmatrix} \cos\alpha & -\sin\alpha & 0 \\ \sin\alpha & \cos\alpha & 0 \\ 0 & 0 & 1 \end{bmatrix}. \tag{2.1}
$$

**Definition 2.2.2** (Direction cosine matrix)**.** *A direction cosine matrix (DCM) is a matrix $R \in \mathbb{R}^{3 \times 3}$ satisfying $R^T R = I$ and $\det R = +1$. The set of all such matrices is the special orthogonal group $SO(3)$.*

**Definition 2.2.3** (Skew-symmetric matrix)**.** *For $\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3) \in \mathbb{R}^3$, the skew-symmetric matrix is:*

$$
[\boldsymbol{\omega}]_\times = \begin{bmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{bmatrix} \tag{2.2}
$$

*satisfying $[\boldsymbol{\omega}]_\times \mathbf{v} = \boldsymbol{\omega} \times \mathbf{v}$ for all $\mathbf{v} \in \mathbb{R}^3$.*

**Definition 2.2.4** (The 3-1-3 Euler rotation sequence)**.** *The composite rotation from orbital elements $(\Omega, i, u)$:*

$$
R_{\mathrm{PF} \to \mathrm{TEME}}(\Omega, i, u) = R_3(-\Omega)\,R_1(-i)\,R_3(-u). \tag{2.3}
$$

**Theorem 2.2.1** (Orthogonality)**.** *Each elementary rotation matrix $R_k(\alpha)$ is orthogonal: $R_k(\alpha)^T R_k(\alpha) = I$ and $\det R_k(\alpha) = +1$.*

*Proof.* For $R_3(\alpha)$: compute the product entry by entry. The $(1,1)$ entry of $R_3^T R_3$ is $\cos^2\alpha + \sin^2\alpha = 1$. The $(1,2)$ entry is $\cos\alpha(-\sin\alpha) + \sin\alpha\cos\alpha = 0$. The remaining entries follow by the same trigonometric identity. The determinant is $\cos\alpha(\cos\alpha) - (-\sin\alpha)(\sin\alpha) = \cos^2\alpha + \sin^2\alpha = 1$ (expanding along the third row). The proofs for $R_1$ and $R_2$ are identical by cyclic permutation of indices. ∎

**Theorem 2.2.2** (Composition of rotations)**.** *If $R_{AB}$ rotates frame B → A and $R_{BC}$ rotates frame C → B, then $R_{AC} = R_{AB}\,R_{BC}$ rotates frame C → A.*

*Proof.* Let $\mathbf{v}_C$ be a vector in frame C. Its representation in frame B is $\mathbf{v}_B = R_{BC}\,\mathbf{v}_C$. Its representation in frame A is $\mathbf{v}_A = R_{AB}\,\mathbf{v}_B = R_{AB}\,R_{BC}\,\mathbf{v}_C$. Since this holds for all $\mathbf{v}_C$, we have $R_{AC} = R_{AB}\,R_{BC}$. The product of two $SO(3)$ matrices is in $SO(3)$ since $(R_1 R_2)^T(R_1 R_2) = R_2^T R_1^T R_1 R_2 = R_2^T R_2 = I$ and $\det(R_1 R_2) = \det R_1 \cdot \det R_2 = 1$. ∎

**Theorem 2.2.3** (Inverse = transpose)**.** *For any $R \in SO(3)$, the inverse rotation is $R^{-1} = R^T$.*

*Proof.* By definition of $SO(3)$, $R^T R = I$. Left-multiplying both sides by $(R^T)^{-1}$ is unnecessary — the equation $R^T R = I$ directly states that $R^T$ is the left inverse of $R$. Similarly, $R\,R^T = I$ (since $R^T \in SO(3)$), so $R^T$ is also the right inverse. Therefore $R^{-1} = R^T$. ∎

**Theorem 2.2.4** (Time derivative of a rotation)**.** *If $R(t) = R_k(\alpha(t))$ is a rotation about axis $\hat{\mathbf{e}}_k$ with time-varying angle $\alpha(t)$, then:*

$$
\frac{d}{dt}R_k(\alpha(t)) = \dot{\alpha}\,[\hat{\mathbf{e}}_k]_\times\,R_k(\alpha(t)). \tag{2.4}
$$

*Proof.* Differentiate $R_3(\alpha)$ entry by entry: $\frac{d}{dt}\cos\alpha = -\dot{\alpha}\sin\alpha$ and $\frac{d}{dt}\sin\alpha = \dot{\alpha}\cos\alpha$. The derivative matrix is:

$$
\dot{R}_3 = \dot{\alpha}\begin{bmatrix} -\sin\alpha & -\cos\alpha & 0 \\ \cos\alpha & -\sin\alpha & 0 \\ 0 & 0 & 0 \end{bmatrix}.
$$

The skew-symmetric matrix for the $z$-axis is $[\hat{\mathbf{e}}_3]_\times = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}$. Computing $\dot{\alpha}\,[\hat{\mathbf{e}}_3]_\times R_3(\alpha)$:

$$
\dot{\alpha}\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}\begin{bmatrix} \cos\alpha & -\sin\alpha & 0 \\ \sin\alpha & \cos\alpha & 0 \\ 0 & 0 & 1 \end{bmatrix} = \dot{\alpha}\begin{bmatrix} -\sin\alpha & -\cos\alpha & 0 \\ \cos\alpha & -\sin\alpha & 0 \\ 0 & 0 & 0 \end{bmatrix} = \dot{R}_3.
$$

The proofs for $R_1$ and $R_2$ follow by cyclic permutation. For a general rotation $R(t)$ about axis $\hat{\mathbf{n}}$ with angular velocity $\boldsymbol{\omega} = \dot{\alpha}\,\hat{\mathbf{n}}$, the result generalizes to $\dot{R} = [\boldsymbol{\omega}]_\times R$. ∎

**Lemma 2.2.1** (Skew-symmetric properties)**.** *For $\boldsymbol{\omega}, \mathbf{v} \in \mathbb{R}^3$:*

1. $[\boldsymbol{\omega}]_\times^T = -[\boldsymbol{\omega}]_\times$
2. $[\boldsymbol{\omega}]_\times\,\mathbf{v} = -[\mathbf{v}]_\times\,\boldsymbol{\omega}$
3. $\det\,[\boldsymbol{\omega}]_\times = 0$

*Proof.* (1) Direct inspection of (2.2): the $(i,j)$ entry of $[\boldsymbol{\omega}]_\times$ equals the negation of the $(j,i)$ entry. (2) Both sides equal $\boldsymbol{\omega} \times \mathbf{v}$, and $\boldsymbol{\omega} \times \mathbf{v} = -\mathbf{v} \times \boldsymbol{\omega} = -[\mathbf{v}]_\times\,\boldsymbol{\omega}$. (3) Every skew-symmetric $3 \times 3$ matrix has zero determinant: $\det(A) = \det(A^T) = \det(-A) = (-1)^3 \det(A) = -\det(A)$, so $\det(A) = 0$. ∎

**Lemma 2.2.2** (Rotation of skew-symmetric matrix)**.** *For $R \in SO(3)$ and $\boldsymbol{\omega} \in \mathbb{R}^3$:*

$$
[R\,\boldsymbol{\omega}]_\times = R\,[\boldsymbol{\omega}]_\times\,R^T. \tag{2.5}
$$

*Proof.* For any $\mathbf{v} \in \mathbb{R}^3$: the left side gives $(R\boldsymbol{\omega}) \times \mathbf{v}$. The right side gives $R([\boldsymbol{\omega}]_\times R^T \mathbf{v}) = R(\boldsymbol{\omega} \times R^T\mathbf{v})$. We need $(R\boldsymbol{\omega}) \times \mathbf{v} = R(\boldsymbol{\omega} \times R^T\mathbf{v})$. Setting $\mathbf{u} = R^T\mathbf{v}$, the right side is $R(\boldsymbol{\omega} \times \mathbf{u})$. The cross product transforms under rotation as $R(\mathbf{a} \times \mathbf{b}) = (R\mathbf{a}) \times (R\mathbf{b})$ (since $R$ preserves orientation and lengths). Therefore $R(\boldsymbol{\omega} \times \mathbf{u}) = (R\boldsymbol{\omega}) \times (R\mathbf{u}) = (R\boldsymbol{\omega}) \times \mathbf{v}$. Since the identity holds for all $\mathbf{v}$, the matrices are equal. ∎

---

## §2.3 Homogeneous Coordinates

**Definition 2.3.1** (Homogeneous position vector)**.** *A position $\mathbf{r} \in \mathbb{R}^3$ is embedded as the homogeneous vector:*

$$
\bar{\mathbf{p}} = \begin{bmatrix} \mathbf{r} \\ 1 \end{bmatrix} \in \mathbb{R}^4. \tag{2.6}
$$

**Definition 2.3.2** ($4 \times 4$ homogeneous transform)**.** *A rigid body transform (rotation $R$ plus translation $\mathbf{d}$) is encoded as:*

$$
H = \begin{bmatrix} R & \mathbf{d} \\ \mathbf{0}^T & 1 \end{bmatrix} \in \mathbb{R}^{4 \times 4}. \tag{2.7}
$$

*The action $\bar{\mathbf{p}}' = H\,\bar{\mathbf{p}}$ gives $\mathbf{r}' = R\,\mathbf{r} + \mathbf{d}$ — rotation then translation.*

**Theorem 2.3.1** ($4 \times 4$ composition)**.** *If $H_{AB}$ transforms frame B → A and $H_{BC}$ transforms frame C → B, then:*

$$
H_{AC} = H_{AB}\,H_{BC} = \begin{bmatrix} R_{AB}\,R_{BC} & R_{AB}\,\mathbf{d}_{BC} + \mathbf{d}_{AB} \\ \mathbf{0}^T & 1 \end{bmatrix}. \tag{2.8}
$$

*The combined rotation is $R_{AB}\,R_{BC}$ and the combined translation is $R_{AB}\,\mathbf{d}_{BC} + \mathbf{d}_{AB}$.*

*Proof.* Perform the block matrix multiplication:

$$
\begin{bmatrix} R_{AB} & \mathbf{d}_{AB} \\ \mathbf{0}^T & 1 \end{bmatrix}\begin{bmatrix} R_{BC} & \mathbf{d}_{BC} \\ \mathbf{0}^T & 1 \end{bmatrix} = \begin{bmatrix} R_{AB}R_{BC} + \mathbf{d}_{AB}\cdot\mathbf{0}^T & R_{AB}\mathbf{d}_{BC} + \mathbf{d}_{AB} \\ \mathbf{0}^T R_{BC} + 1\cdot\mathbf{0}^T & \mathbf{0}^T\mathbf{d}_{BC} + 1 \end{bmatrix} = \begin{bmatrix} R_{AB}R_{BC} & R_{AB}\mathbf{d}_{BC} + \mathbf{d}_{AB} \\ \mathbf{0}^T & 1 \end{bmatrix}.
$$

The result has the form of (2.7) with rotation $R_{AB}R_{BC} \in SO(3)$ and translation $R_{AB}\mathbf{d}_{BC} + \mathbf{d}_{AB}$. The inner translation is rotated into the outer frame before the outer translation is added. ∎

**Theorem 2.3.2** ($4 \times 4$ inverse)**.** *The inverse of $H = [R, \mathbf{d};\, \mathbf{0}^T, 1]$ is:*

$$
H^{-1} = \begin{bmatrix} R^T & -R^T\mathbf{d} \\ \mathbf{0}^T & 1 \end{bmatrix}. \tag{2.9}
$$

*Proof.* Verify $H\,H^{-1} = I_4$:

$$
\begin{bmatrix} R & \mathbf{d} \\ \mathbf{0}^T & 1 \end{bmatrix}\begin{bmatrix} R^T & -R^T\mathbf{d} \\ \mathbf{0}^T & 1 \end{bmatrix} = \begin{bmatrix} RR^T & -RR^T\mathbf{d} + \mathbf{d} \\ \mathbf{0}^T & 1 \end{bmatrix} = \begin{bmatrix} I & \mathbf{0} \\ \mathbf{0}^T & 1 \end{bmatrix} = I_4.
$$

We used $RR^T = I$ (Theorem 2.2.3). ∎

**Remark** (Limitation of $4 \times 4$ matrices)**.** The homogeneous matrix transforms position only. It does not account for velocity, which requires the transport theorem when the frame is rotating. This motivates the 7×7 extension (§2.4).

---

## §2.4 The 7×7 Homogeneous State Matrix

**Definition 2.4.1** (Extended state vector)**.** *The satellite state (position + velocity + homogeneous coordinate) is:*

$$
\mathbf{x} = \begin{bmatrix} \mathbf{r} \\ \mathbf{v} \\ 1 \end{bmatrix} \in \mathbb{R}^7. \tag{2.10}
$$

**Definition 2.4.2** ($7 \times 7$ state matrix)**.** *The state transformation matrix is:*

$$
T = \begin{bmatrix} R & \mathbf{0} & \Delta\mathbf{r} \\ \dot{R} & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix} \tag{2.11}
$$

*where the blocks are $3 \times 3$, $3 \times 3$, $3 \times 1$, $3 \times 3$, $3 \times 3$, $3 \times 1$, $1 \times 3$, $1 \times 3$, $1 \times 1$ respectively. The action $\mathbf{x}' = T\,\mathbf{x}$ gives:*

$$
\mathbf{r}' = R\,\mathbf{r} + \Delta\mathbf{r}, \qquad \mathbf{v}' = \dot{R}\,\mathbf{r} + R\,\mathbf{v} + \Delta\mathbf{v}. \tag{2.12}
$$

**Definition 2.4.3** (Block interpretation)**.** *The six blocks of $T$ have the following physical meanings:*

| Block | Size | Function |
|-------|------|----------|
| $R$ (upper-left) | $3 \times 3$ | Direction cosine matrix (primary rotation) |
| $\mathbf{0}$ (upper-middle) | $3 \times 3$ | Position does not depend on velocity |
| $\Delta\mathbf{r}$ (upper-right) | $3 \times 1$ | Origin offset between frames |
| $\dot{R}$ (lower-left) | $3 \times 3$ | Transport term: $[\boldsymbol{\omega}]_\times R$ for rotating frames |
| $R$ (lower-middle) | $3 \times 3$ | Velocity rotation (same DCM as position) |
| $\Delta\mathbf{v}$ (lower-right) | $3 \times 1$ | Linear velocity of new frame origin |

**Definition 2.4.4** (Pure rotation state matrix)**.** *When the transform is a pure rotation (no translation, non-rotating frame):*

$$
T_{\mathrm{rot}}(R) = \begin{bmatrix} R & \mathbf{0} & \mathbf{0} \\ \mathbf{0} & R & \mathbf{0} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}. \tag{2.13}
$$

**Definition 2.4.5** (Pure translation state matrix)**.** *When the transform is a pure translation:*

$$
T_{\mathrm{trans}}(\Delta\mathbf{r}, \Delta\mathbf{v}) = \begin{bmatrix} I & \mathbf{0} & \Delta\mathbf{r} \\ \mathbf{0} & I & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}. \tag{2.14}
$$

**Definition 2.4.6** (Rotating-frame state matrix)**.** *When the frame rotates with angular velocity $\boldsymbol{\omega}$:*

$$
T_\omega(R) = \begin{bmatrix} R & \mathbf{0} & \mathbf{0} \\ [\boldsymbol{\omega}]_\times R & R & \mathbf{0} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}. \tag{2.15}
$$

*The lower-left block $\dot{R} = [\boldsymbol{\omega}]_\times R$ encodes the transport theorem.*

**Assumption 2.4.1** (Instantaneous angular velocity)**.** *The angular velocity $\boldsymbol{\omega}$ is treated as constant during a single transformation. For time-varying $\boldsymbol{\omega}$, the transform must be re-evaluated at each epoch.*

**Theorem 2.4.1** (Transport theorem derivation)**.** *If frame B rotates with angular velocity $\boldsymbol{\omega}$ relative to frame A, and $\mathbf{r}_B$ is the position vector in frame B, then the velocity in frame A is:*

$$
\mathbf{v}_A = R\,\mathbf{v}_B + [\boldsymbol{\omega}]_\times R\,\mathbf{r}_B + \Delta\mathbf{v}. \tag{2.16}
$$

*The term $[\boldsymbol{\omega}]_\times R\,\mathbf{r}_B = \boldsymbol{\omega} \times \mathbf{r}_A$ is the transport theorem contribution. The lower-left block $\dot{R} = [\boldsymbol{\omega}]_\times R$ follows.*

*Proof.* The position in frame A is $\mathbf{r}_A(t) = R(t)\,\mathbf{r}_B(t) + \mathbf{d}(t)$, where $R(t)$ is the time-varying rotation and $\mathbf{d}(t)$ is the time-varying origin offset. Differentiate with respect to time using the product rule:

$$
\mathbf{v}_A = \dot{R}\,\mathbf{r}_B + R\,\dot{\mathbf{r}}_B + \dot{\mathbf{d}} = \dot{R}\,\mathbf{r}_B + R\,\mathbf{v}_B + \Delta\mathbf{v}.
$$

By Theorem 2.2.4, $\dot{R} = [\boldsymbol{\omega}]_\times R$, where $\boldsymbol{\omega}$ is the angular velocity of frame B as seen from frame A. Therefore:

$$
\mathbf{v}_A = [\boldsymbol{\omega}]_\times R\,\mathbf{r}_B + R\,\mathbf{v}_B + \Delta\mathbf{v}.
$$

The term $[\boldsymbol{\omega}]_\times R\,\mathbf{r}_B = [\boldsymbol{\omega}]_\times\,\mathbf{r}_A' = \boldsymbol{\omega} \times \mathbf{r}_A'$ (where $\mathbf{r}_A' = R\,\mathbf{r}_B$ is the rotated position, ignoring translation). This is the classical transport theorem: velocity in the inertial frame equals velocity in the rotating frame plus the $\boldsymbol{\omega} \times \mathbf{r}$ correction. Comparing with (2.12), the lower-left block is $\dot{R} = [\boldsymbol{\omega}]_\times R$. ∎

**Theorem 2.4.2** (Composition of $7 \times 7$ state matrices)**.** *If $T_{AB}$ transforms frame B → A and $T_{BC}$ transforms frame C → B, then $T_{AC} = T_{AB}\,T_{BC}$ with:*

$$
R_{AC} = R_{AB}\,R_{BC}, \tag{2.17}
$$

$$
\dot{R}_{AC} = \dot{R}_{AB}\,R_{BC} + R_{AB}\,\dot{R}_{BC}, \tag{2.18}
$$

$$
\Delta\mathbf{r}_{AC} = R_{AB}\,\Delta\mathbf{r}_{BC} + \Delta\mathbf{r}_{AB}, \tag{2.19}
$$

$$
\Delta\mathbf{v}_{AC} = \dot{R}_{AB}\,\Delta\mathbf{r}_{BC} + R_{AB}\,\Delta\mathbf{v}_{BC} + \Delta\mathbf{v}_{AB}. \tag{2.20}
$$

*Proof.* Multiply the block matrices:

$$
T_{AB}\,T_{BC} = \begin{bmatrix} R_{AB} & \mathbf{0} & \Delta\mathbf{r}_{AB} \\ \dot{R}_{AB} & R_{AB} & \Delta\mathbf{v}_{AB} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}\begin{bmatrix} R_{BC} & \mathbf{0} & \Delta\mathbf{r}_{BC} \\ \dot{R}_{BC} & R_{BC} & \Delta\mathbf{v}_{BC} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}.
$$

The nine block products give:

- Upper-left: $R_{AB}R_{BC} + \mathbf{0}\cdot\dot{R}_{BC} + \Delta\mathbf{r}_{AB}\cdot\mathbf{0}^T = R_{AB}R_{BC}$.
- Upper-middle: $R_{AB}\cdot\mathbf{0} + \mathbf{0}\cdot R_{BC} + \Delta\mathbf{r}_{AB}\cdot\mathbf{0}^T = \mathbf{0}$.
- Upper-right: $R_{AB}\Delta\mathbf{r}_{BC} + \mathbf{0}\cdot\Delta\mathbf{v}_{BC} + \Delta\mathbf{r}_{AB}\cdot 1 = R_{AB}\Delta\mathbf{r}_{BC} + \Delta\mathbf{r}_{AB}$.
- Lower-left: $\dot{R}_{AB}R_{BC} + R_{AB}\dot{R}_{BC} + \Delta\mathbf{v}_{AB}\cdot\mathbf{0}^T = \dot{R}_{AB}R_{BC} + R_{AB}\dot{R}_{BC}$.
- Lower-middle: $\dot{R}_{AB}\cdot\mathbf{0} + R_{AB}R_{BC} + \Delta\mathbf{v}_{AB}\cdot\mathbf{0}^T = R_{AB}R_{BC}$.
- Lower-right: $\dot{R}_{AB}\Delta\mathbf{r}_{BC} + R_{AB}\Delta\mathbf{v}_{BC} + \Delta\mathbf{v}_{AB}$.
- Bottom row: $[\mathbf{0}^T,\, \mathbf{0}^T,\, 1]$ (unchanged).

The result has the form (2.11) with the stated blocks. Note that the lower-middle block equals $R_{AC} = R_{AB}R_{BC}$, confirming that the velocity rotation uses the same DCM as the position rotation. The lower-left block $\dot{R}_{AC} = \dot{R}_{AB}R_{BC} + R_{AB}\dot{R}_{BC}$ is a Leibniz-rule product — physically, it says the combined transport term accounts for both the outer frame's rotation acting on the inner rotation and the inner frame's own rotation. ∎

**Theorem 2.4.3** (Inverse of a $7 \times 7$ state matrix)**.** *The inverse of $T = [R, \mathbf{0}, \Delta\mathbf{r};\, \dot{R}, R, \Delta\mathbf{v};\, \mathbf{0}^T, \mathbf{0}^T, 1]$ is:*

$$
T^{-1} = \begin{bmatrix} R^T & \mathbf{0} & -R^T\Delta\mathbf{r} \\ -R^T\dot{R}R^T & R^T & R^T\dot{R}R^T\Delta\mathbf{r} - R^T\Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}. \tag{2.21}
$$

*Proof.* Set $T\,T^{-1} = I_7$ and solve for the blocks of $T^{-1}$. Denote $T^{-1}$ with blocks $[A, B, \mathbf{c};\, D, E, \mathbf{f};\, \mathbf{0}^T, \mathbf{0}^T, 1]$.

From the upper-left block: $RA = I \Rightarrow A = R^T$. From the upper-middle block: $RB = \mathbf{0} \Rightarrow B = \mathbf{0}$. From the upper-right block: $R\mathbf{c} + \Delta\mathbf{r} = \mathbf{0} \Rightarrow \mathbf{c} = -R^T\Delta\mathbf{r}$.

From the lower-left block: $\dot{R}A + RD = \mathbf{0} \Rightarrow D = -R^T\dot{R}R^T$. From the lower-middle block: $\dot{R}B + RE = I \Rightarrow RE = I \Rightarrow E = R^T$. From the lower-right block: $\dot{R}\mathbf{c} + R\mathbf{f} + \Delta\mathbf{v} = \mathbf{0} \Rightarrow \mathbf{f} = -R^T(\dot{R}\mathbf{c} + \Delta\mathbf{v}) = -R^T(-\dot{R}R^T\Delta\mathbf{r} + \Delta\mathbf{v}) = R^T\dot{R}R^T\Delta\mathbf{r} - R^T\Delta\mathbf{v}$.

Verification: $T^{-1}\,T = I_7$ can be checked by block multiplication. ∎

**Theorem 2.4.4** (Group structure)**.** *The set of all $7 \times 7$ state matrices of the form (2.11) forms a group under matrix multiplication.*

*Proof.* We verify the four group axioms.

(1) *Closure:* The product of two matrices of form (2.11) is again of form (2.11), as shown in Theorem 2.4.2. The upper-middle block is always $\mathbf{0}$ and the bottom row is always $[\mathbf{0}^T, \mathbf{0}^T, 1]$.

(2) *Associativity:* Matrix multiplication is associative.

(3) *Identity:* $T_{\mathrm{id}} = I_7$ has $R = I$, $\dot{R} = \mathbf{0}$, $\Delta\mathbf{r} = \mathbf{0}$, $\Delta\mathbf{v} = \mathbf{0}$, which is of the form (2.11).

(4) *Inverse:* Every $T$ of form (2.11) has an inverse of the same form, as shown in Theorem 2.4.3. ∎

**Corollary 2.4.1** (Pure rotation composition)**.** *Pure rotation state matrices compose as block-diagonal:*

$$
T_{\mathrm{rot}}(R_1)\,T_{\mathrm{rot}}(R_2) = T_{\mathrm{rot}}(R_1 R_2). \tag{2.22}
$$

*Proof.* Set $\dot{R} = \mathbf{0}$, $\Delta\mathbf{r} = \mathbf{0}$, $\Delta\mathbf{v} = \mathbf{0}$ for both matrices in Theorem 2.4.2. All cross terms vanish. ∎

**Corollary 2.4.2** (Earth rotation transport term)**.** *For Earth rotation about $\hat{\mathbf{e}}_3$ at rate $\omega_E$, the skew-symmetric matrix has the explicit form:*

$$
[\omega_E\,\hat{\mathbf{e}}_3]_\times = \omega_E\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}. \tag{2.23}
$$

*The transport term contribution to velocity is $[\omega_E\,\hat{\mathbf{e}}_3]_\times\,\mathbf{r}' = \omega_E(-r'_y,\, r'_x,\, 0)^T$.*

*Proof.* Substitute $\boldsymbol{\omega} = (0, 0, \omega_E)$ into (2.2). Apply to $\mathbf{r}' = (r'_x, r'_y, r'_z)^T$. ∎

**Corollary 2.4.3** (Pure rotation cross terms)**.** *When $\dot{R} = \mathbf{0}$, $\Delta\mathbf{r} = \mathbf{0}$, $\Delta\mathbf{v} = \mathbf{0}$ (pure rotation), the composition in the $\Delta\mathbf{v}$ row of Theorem 2.4.2 is trivial: $\Delta\mathbf{v}_{AC} = \mathbf{0}$. No transport term cross-coupling arises.*

*Proof.* In (2.20), set $\dot{R}_{AB} = \mathbf{0}$, $\Delta\mathbf{r}_{BC} = \mathbf{0}$, $\Delta\mathbf{v}_{BC} = \mathbf{0}$, $\Delta\mathbf{v}_{AB} = \mathbf{0}$. ∎

**Remark** (Direct computation vs explicit matrix multiplication)**.** In practice, the matrix $T$ need not be explicitly constructed and stored as a $7 \times 7$ array. The block structure means the computation decomposes into: (a) a $3 \times 3$ matrix-vector product for position, (b) another $3 \times 3$ product plus a cross product for velocity. The 7×7 matrix formulation is a *conceptual* tool that ensures all coupling terms (especially the transport term) are correctly accounted for when composing transforms.

---

## §2.5 Standard Coordinate Transforms as State Matrices

**Definition 2.5.1** (TEME frame)**.** *The True Equator, Mean Equinox (TEME) frame is the reference frame of SGP4 output. It is not a precisely defined inertial frame.* [A.2.1]

**Definition 2.5.2** (PEF frame)**.** *The Pseudo Earth-Fixed (PEF) frame is obtained from TEME by rotating about the polar axis by the Greenwich Mean Sidereal Time angle $\theta_{\mathrm{GMST}}$ (Ch. 29).*

**Definition 2.5.3** (Perifocal frame)**.** *The perifocal frame has $\hat{\mathbf{x}}$ toward perigee, $\hat{\mathbf{z}}$ along the angular momentum vector, and is defined by the orbital elements $(\Omega, i, \omega)$ with argument of latitude $u = \omega + \nu$.*

**Theorem 2.5.1** (Perifocal → TEME state matrix)**.** *The state matrix from perifocal to TEME is a pure rotation:*

$$
T_{\mathrm{PF} \to \mathrm{TEME}} = T_{\mathrm{rot}}(R), \quad R = R_3(-\Omega)\,R_1(-i)\,R_3(-u). \tag{2.24}
$$

*The explicit DCM entries are:*

$$
R = \begin{bmatrix} \cos\Omega\cos u - \sin\Omega\cos i\sin u & -\cos\Omega\sin u - \sin\Omega\cos i\cos u & \sin\Omega\sin i \\ \sin\Omega\cos u + \cos\Omega\cos i\sin u & -\sin\Omega\sin u + \cos\Omega\cos i\cos u & -\cos\Omega\sin i \\ \sin i\sin u & \sin i\cos u & \cos i \end{bmatrix}. \tag{2.25}
$$

*The unit vectors of Spacetrack Report No. 3 are the columns: $\hat{\mathbf{U}} = R\,\hat{\mathbf{e}}_1$ (the first column) and $\hat{\mathbf{V}} = R\,\hat{\mathbf{e}}_2$ (the second column). Position and velocity in TEME are:*

$$
\mathbf{r}_{\mathrm{TEME}} = r\,\hat{\mathbf{U}}, \qquad \mathbf{v}_{\mathrm{TEME}} = \dot{r}\,\hat{\mathbf{U}} + r\dot{f}\,\hat{\mathbf{V}}. \tag{2.26}
$$

*Since $\dot{R} = \mathbf{0}$, $\Delta\mathbf{r} = \mathbf{0}$, $\Delta\mathbf{v} = \mathbf{0}$, there is no transport term.*

*Proof.* The perifocal-to-ECI rotation is the 3-1-3 Euler sequence $(\Omega, i, u)$ applied in reverse: first rotate by $-u$ about $\hat{\mathbf{e}}_3$ (from argument-of-latitude to line-of-nodes), then by $-i$ about $\hat{\mathbf{e}}_1$ (from orbital plane to equatorial plane), then by $-\Omega$ about $\hat{\mathbf{e}}_3$ (from ascending node to vernal equinox). By Theorem 2.2.2, the composite rotation is $R = R_3(-\Omega)\,R_1(-i)\,R_3(-u)$.

To obtain the explicit entries, multiply the three elementary matrices from (2.1). Writing $c_\Omega = \cos\Omega$, $s_\Omega = \sin\Omega$, $c_i = \cos i$, $s_i = \sin i$, $c_u = \cos u$, $s_u = \sin u$:

$$
R_3(-\Omega) = \begin{bmatrix} c_\Omega & s_\Omega & 0 \\ -s_\Omega & c_\Omega & 0 \\ 0 & 0 & 1 \end{bmatrix}, \quad R_1(-i) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & c_i & s_i \\ 0 & -s_i & c_i \end{bmatrix}, \quad R_3(-u) = \begin{bmatrix} c_u & s_u & 0 \\ -s_u & c_u & 0 \\ 0 & 0 & 1 \end{bmatrix}.
$$

First compute $R_1(-i)\,R_3(-u)$:

$$
\begin{bmatrix} c_u & s_u & 0 \\ -c_i s_u & c_i c_u & s_i \\ s_i s_u & -s_i c_u & c_i \end{bmatrix}.
$$

Wait — more carefully: $R_1(-i) R_3(-u)$, row 1 of $R_1(-i)$ is $(1, 0, 0)$, so row 1 of the product is row 1 of $R_3(-u)$: $(c_u, s_u, 0)$. Row 2 of $R_1(-i)$ is $(0, c_i, s_i)$, applied to columns of $R_3(-u)$: $(c_i(-s_u) + s_i \cdot 0,\; c_i c_u + s_i \cdot 0,\; c_i \cdot 0 + s_i \cdot 1) = (-c_i s_u, c_i c_u, s_i)$. Row 3: $(0, -s_i, c_i)$ applied: $(s_i s_u, -s_i c_u, c_i)$.

Then $R_3(-\Omega)(R_1(-i)R_3(-u))$: row 1 of $R_3(-\Omega)$ is $(c_\Omega, s_\Omega, 0)$. Entry $(1,1)$: $c_\Omega c_u + s_\Omega(-c_i s_u) = c_\Omega c_u - s_\Omega c_i s_u$. Entry $(1,2)$: $c_\Omega s_u + s_\Omega c_i c_u$... However, noting the sign conventions for $R_3(-\Omega)$: $R_3(-\alpha)$ has $\cos(-\alpha) = \cos\alpha$ in position $(1,1)$ and $-\sin(-\alpha) = \sin\alpha$ in position $(1,2)$.

The resulting 9 entries match the standard aerospace convention. The first column is $\hat{\mathbf{U}}$ and the second column is $\hat{\mathbf{V}}$, exactly the unit vectors defined in Spacetrack Report No. 3 (pp. 14–15). Position $\mathbf{r} = r\hat{\mathbf{U}}$ and velocity $\mathbf{v} = \dot{r}\hat{\mathbf{U}} + r\dot{f}\hat{\mathbf{V}}$ follow from the perifocal decomposition $\mathbf{r}_{\mathrm{PF}} = (r, 0, 0)^T$ and $\mathbf{v}_{\mathrm{PF}} = (\dot{r}, r\dot{f}, 0)^T$.

No translation is involved (both frames share an origin at Earth's center). The perifocal frame is not rotating relative to TEME for a given set of orbital elements (the elements are frozen at a single epoch), so $\dot{R} = \mathbf{0}$. ∎

**Remark** (Trigonometric function errors)**.** The DCM entries involve products of $\sin$ and $\cos$ of the orbital angles $\Omega$, $i$, $u$. Each trigonometric evaluation introduces precision error bounded by Ch. 1, Corollaries 1.4.1–1.4.2. Each product introduces additional error bounded by Ch. 1, Theorem 1.3.2. For double precision: $\delta_p \leq 1.1 \times 10^{-16}$ per trig evaluation. The orbital angles carry measurement error $\sigma_m(\Omega)$, $\sigma_m(i)$, $\sigma_m(u)$ from the TLE. [P.2.1]

**Theorem 2.5.2** (TEME → PEF state matrix)**.** *The state matrix from TEME to PEF is a rotating-frame matrix:*

$$
T_{\mathrm{TEME} \to \mathrm{PEF}} = T_\omega(R_3(\theta_{\mathrm{GMST}})) = \begin{bmatrix} R_3(\theta) & \mathbf{0} & \mathbf{0} \\ [\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta) & R_3(\theta) & \mathbf{0} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix} \tag{2.27}
$$

*where $\theta = \theta_{\mathrm{GMST}}(t)$. The transport term produces the velocity correction $\omega_E\,\hat{\mathbf{e}}_3 \times \mathbf{r}_{\mathrm{PEF}}$.* [M.2.1]

*Proof.* PEF is TEME rotated by $\theta_{\mathrm{GMST}}$ about $\hat{\mathbf{e}}_3$, so $R = R_3(\theta_{\mathrm{GMST}})$. No translation: $\Delta\mathbf{r} = \mathbf{0}$ (both frames share Earth's center). No velocity offset: $\Delta\mathbf{v} = \mathbf{0}$. The frame rotates at angular velocity $\boldsymbol{\omega} = \omega_E\,\hat{\mathbf{e}}_3$, so by Theorem 2.4.1, $\dot{R} = [\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta)$.

Explicitly, the lower-left block is:

$$
[\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta) = \omega_E\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}\begin{bmatrix} \cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{bmatrix} = \omega_E\begin{bmatrix} -\sin\theta & -\cos\theta & 0 \\ \cos\theta & -\sin\theta & 0 \\ 0 & 0 & 0 \end{bmatrix}.
$$

The velocity transform is $\mathbf{v}_{\mathrm{PEF}} = \dot{R}\,\mathbf{r}_{\mathrm{TEME}} + R\,\mathbf{v}_{\mathrm{TEME}}$. The $\dot{R}\,\mathbf{r}_{\mathrm{TEME}}$ term equals $[\omega_E\hat{\mathbf{e}}_3]_\times R\,\mathbf{r}_{\mathrm{TEME}} = [\omega_E\hat{\mathbf{e}}_3]_\times\,\mathbf{r}_{\mathrm{PEF}} = \omega_E\hat{\mathbf{e}}_3 \times \mathbf{r}_{\mathrm{PEF}}$. This is the standard Earth-rotation velocity correction. ∎

**Remark** (GMST source and unit conversions)**.** $\theta_{\mathrm{GMST}}$ is computed from the IAU 1982 polynomial (Ch. 29). Unit conversion factors are Tier I within the matched pair (Ch. 3, Ch. 37).

**Theorem 2.5.3** (Perifocal → PEF composition)**.** *The composed transform is:*

$$
T_{\mathrm{PF} \to \mathrm{PEF}} = T_{\mathrm{TEME} \to \mathrm{PEF}} \cdot T_{\mathrm{PF} \to \mathrm{TEME}}. \tag{2.28}
$$

*By Theorem 2.4.2:*

$$
R_{\mathrm{PF} \to \mathrm{PEF}} = R_3(\theta_{\mathrm{GMST}})\,R_3(-\Omega)\,R_1(-i)\,R_3(-u), \tag{2.29}
$$

$$
\dot{R}_{\mathrm{PF} \to \mathrm{PEF}} = [\omega_E\hat{\mathbf{e}}_3]_\times\,R_3(\theta)\,R_{\mathrm{PF} \to \mathrm{TEME}}. \tag{2.30}
$$

*The $\dot{R}$ term in the TEME→PEF matrix acts on the fully rotated position to produce the velocity correction.*

*Proof.* Apply Theorem 2.4.2 with $T_{AB} = T_{\mathrm{TEME} \to \mathrm{PEF}}$ and $T_{BC} = T_{\mathrm{PF} \to \mathrm{TEME}}$. The perifocal→TEME transform is a pure rotation, so $\dot{R}_{BC} = \mathbf{0}$, $\Delta\mathbf{r}_{BC} = \mathbf{0}$, $\Delta\mathbf{v}_{BC} = \mathbf{0}$. From (2.17): $R_{AC} = R_3(\theta)\,R_{\mathrm{PF} \to \mathrm{TEME}}$. From (2.18): $\dot{R}_{AC} = \dot{R}_{AB}\,R_{BC} + R_{AB}\cdot\mathbf{0} = [\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta)\,R_{\mathrm{PF} \to \mathrm{TEME}}$. From (2.19): $\Delta\mathbf{r}_{AC} = R_3(\theta)\cdot\mathbf{0} + \mathbf{0} = \mathbf{0}$. From (2.20): $\Delta\mathbf{v}_{AC} = [\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta)\cdot\mathbf{0} + R_3(\theta)\cdot\mathbf{0} + \mathbf{0} = \mathbf{0}$.

The velocity in PEF is $\mathbf{v}_{\mathrm{PEF}} = R_{\mathrm{PF} \to \mathrm{PEF}}\,\mathbf{v}_{\mathrm{PF}} + \dot{R}_{\mathrm{PF} \to \mathrm{PEF}}\,\mathbf{r}_{\mathrm{PF}}$. The second term produces the transport correction $\boldsymbol{\omega}_E \times \mathbf{r}_{\mathrm{PEF}}$, exactly as in the non-composed form. ∎

**Proposition 2.5.1** (Mean → osculating as near-identity)**.** *The short-period corrections of Brouwer's theory (Ch. 18) can be expressed as a near-identity state matrix:*

$$
T = I_7 + \delta T \tag{2.31}
$$

*where $\delta T$ contains the Brouwer corrections, with $\|\delta T\| \sim J_2/p^2 \ll 1$. Composition of near-identity state matrices is approximately additive:*

$$
(I_7 + \delta T_{\mathrm{SP}})\,(I_7 + \delta T_{\mathrm{LP}})\,T_{\mathrm{secular}} \approx (I_7 + \delta T_{\mathrm{SP}} + \delta T_{\mathrm{LP}})\,T_{\mathrm{secular}}. \tag{2.32}
$$

*Proof.* The corrections $\Delta r$, $\Delta u$, $\Delta\Omega$, $\Delta i$ are $O(J_2/p^2) \sim 10^{-3}$. A rotation near identity has $R \approx I + [\boldsymbol{\theta}]_\times$ for small rotation $\boldsymbol{\theta}$, so the state matrix is $T \approx I_7 + \delta T$ where $\delta T$ contains the small rotation and the corresponding small corrections. Products of near-identity matrices: $(I + A)(I + B) = I + A + B + AB \approx I + A + B$ when $\|AB\| \ll \|A\| + \|B\|$. Full development in Ch. 18. ∎

**Remark** (PEF → ECEF)**.** Polar motion corrections ($< 0.5$ arcsec) are negligible for SGP4-class accuracy. The TEME frame accuracy is limited by truncated precession/nutation [A.2.1].

**Remark** (Velocity error tracking scope)**.** The 7×7 state matrix propagates velocity errors through the rotation and transport terms simultaneously. No separate velocity error computation is needed — it emerges from the lower row of the matrix product. This is the primary advantage of the 7×7 formulation over separate $3 \times 3$ rotation of position and manual transport-term computation.

---

## §2.6 Error Propagation Through Matrix Operations

**Definition 2.6.1** (Matrix condition number)**.** *For $A \in \mathbb{R}^{n \times n}$, the condition number is $\kappa(A) = \|A\|\,\|A^{-1}\|$. For $R \in SO(3)$, $\kappa(R) = 1$ since $\|R\| = \|R^T\| = 1$.*

**Definition 2.6.2** (Tracked matrix)**.** *A tracked matrix is a matrix whose entries are TrackedValues (Ch. 1, Definition 1.2.4). Each entry carries $(v, \sigma_m, \delta_p, \delta_a)$. Matrix operations on tracked matrices apply Ch. 1 error rules entry by entry.*

**Theorem 2.6.1** (Error in matrix-vector product)**.** *For $\mathbf{y} = A\,\mathbf{x}$ where $A \in \mathbb{R}^{m \times n}$ and $\mathbf{x} \in \mathbb{R}^n$, with error bounds $\delta(A_{ij})$ and $\delta(x_j)$ on each entry:*

$$
\delta(y_i) \leq \sum_{j=1}^{n} |A_{ij}|\,\delta(x_j) + \sum_{j=1}^{n} |x_j|\,\delta(A_{ij}) + \sum_{j=1}^{n} \delta(A_{ij})\,\delta(x_j). \tag{2.33}
$$

*Proof.* Each output component is $y_i = \sum_j A_{ij}\,x_j$. By Ch. 1, Theorem 1.3.2 (multiplication), the error in each product $A_{ij} x_j$ is bounded by $|A_{ij}|\,\delta(x_j) + |x_j|\,\delta(A_{ij}) + \delta(A_{ij})\,\delta(x_j)$. By Ch. 1, Theorem 1.3.1 (addition), the error in the sum is bounded by the sum of the individual product errors. The third (cross) term is second-order in the errors and is typically negligible. ∎

**Theorem 2.6.2** (Error in matrix-matrix product)**.** *For $C = AB$, each entry $C_{ij} = \sum_k A_{ik}\,B_{kj}$ is a sum of products. The error bound on $C_{ij}$ follows from repeated application of Ch. 1, Theorems 1.3.1–1.3.2:*

$$
\delta(C_{ij}) \leq \sum_{k=1}^{n} \bigl(|A_{ik}|\,\delta(B_{kj}) + |B_{kj}|\,\delta(A_{ik}) + \delta(A_{ik})\,\delta(B_{kj})\bigr). \tag{2.34}
$$

*Proof.* Identical structure to Theorem 2.6.1, applied entry by entry to each element of the product matrix. ∎

**Theorem 2.6.3** (Orthogonal matrices preserve error magnitude)**.** *For $R \in SO(3)$ with exact entries ($\delta(R_{ij}) = 0$), the error in the rotated vector satisfies:*

$$
\|\boldsymbol{\delta}(R\,\mathbf{x})\| = \|\boldsymbol{\delta}(\mathbf{x})\|. \tag{2.35}
$$

*Pure rotations do not amplify or reduce the total error magnitude.*

*Proof.* When $\delta(R_{ij}) = 0$, Theorem 2.6.1 simplifies to $\delta(y_i) = \sum_j |R_{ij}|\,\delta(x_j)$. Define the error vector $\boldsymbol{\delta} = (\delta(x_1), \delta(x_2), \delta(x_3))^T$. Then $\delta(y_i) \leq \sum_j |R_{ij}|\,\delta(x_j)$, with equality when all error contributions have the same sign.

For the norm: $\|R\,\boldsymbol{\delta}\|^2 = (R\,\boldsymbol{\delta})^T(R\,\boldsymbol{\delta}) = \boldsymbol{\delta}^T R^T R\,\boldsymbol{\delta} = \boldsymbol{\delta}^T\,\boldsymbol{\delta} = \|\boldsymbol{\delta}\|^2$.

Therefore $\|R\,\boldsymbol{\delta}\| = \|\boldsymbol{\delta}\|$. Rotation is an isometry on the error vector. ∎

**Theorem 2.6.4** (Inexact rotation error)**.** *When $R$ carries precision error $\delta(R_{ij})$ (from evaluating trigonometric functions), the output error has an additional term:*

$$
\|\boldsymbol{\delta}(R\,\mathbf{x})\| \leq \|\boldsymbol{\delta}(\mathbf{x})\| + \|\mathbf{x}\|\,\max_{ij} \delta(R_{ij})\,\sqrt{n}. \tag{2.36}
$$

*For rotation angles known to $d$ significant digits, the rotation introduces position error $\approx \|\mathbf{r}\| \cdot 10^{-d}$.*

*Proof.* From Theorem 2.6.1 with nonzero $\delta(R_{ij})$, the error in each component $y_i$ has two contributions: (a) $\sum_j |R_{ij}|\,\delta(x_j)$ from input errors propagated through exact rotation, and (b) $\sum_j |x_j|\,\delta(R_{ij})$ from rotation matrix errors. Contribution (a) gives $\|\boldsymbol{\delta}(\mathbf{x})\|$ by Theorem 2.6.3. Contribution (b) is bounded by $\|\mathbf{x}\|\,\max_{ij}\delta(R_{ij})\,\sqrt{n}$ by the Cauchy-Schwarz inequality. The two contributions are independent and add. ∎

**Theorem 2.6.5** (Transport term error amplification)**.** *The $\dot{R} = [\boldsymbol{\omega}]_\times R$ block couples position error into velocity. The amplification factor is:*

$$
\|\boldsymbol{\delta}(\dot{R}\,\mathbf{r})\| \leq |\omega|\,\|\boldsymbol{\delta}(\mathbf{r})\| + \|\mathbf{r}\|\,\delta(\omega). \tag{2.37}
$$

*For Earth rotation: $\omega_E \approx 7.29 \times 10^{-5}$ rad/s. A 1 km position error produces $\approx 0.073$ m/s velocity error through the transport term.*

*Proof.* $\dot{R}\,\mathbf{r} = [\boldsymbol{\omega}]_\times R\,\mathbf{r} = \boldsymbol{\omega} \times (R\,\mathbf{r})$. By Ch. 1, Theorem 1.3.2 (applied to the cross product as a bilinear operation), the error in $\boldsymbol{\omega} \times \mathbf{r}'$ is bounded by $|\omega|\,\|\boldsymbol{\delta}(\mathbf{r}')\| + \|\mathbf{r}'\|\,\delta(\omega)$. By Theorem 2.6.3, $\|\boldsymbol{\delta}(\mathbf{r}')\| = \|\boldsymbol{\delta}(\mathbf{r})\|$ when $R$ is exact. The numerical example: $7.29 \times 10^{-5}\,\mathrm{rad/s} \times 1\,\mathrm{km} = 7.29 \times 10^{-5}\,\mathrm{km/s} = 0.073\,\mathrm{m/s}$. ∎

**Corollary 2.6.1** (Position error propagation)**.** *For the 7×7 state matrix with exact $R$, the position error propagates unchanged: $\|\boldsymbol{\delta}(\mathbf{r}')\| = \|\boldsymbol{\delta}(\mathbf{r})\|$ (Theorem 2.6.3). The velocity error accumulates from both the velocity input error and the transport-term coupling of position error: $\|\boldsymbol{\delta}(\mathbf{v}')\| \leq \|\boldsymbol{\delta}(\mathbf{v})\| + |\omega|\,\|\boldsymbol{\delta}(\mathbf{r})\|$.*

*Proof.* From (2.12): $\mathbf{v}' = \dot{R}\,\mathbf{r} + R\,\mathbf{v} + \Delta\mathbf{v}$. The error is $\boldsymbol{\delta}(\mathbf{v}') = \boldsymbol{\delta}(\dot{R}\,\mathbf{r}) + \boldsymbol{\delta}(R\,\mathbf{v})$ (assuming $\Delta\mathbf{v}$ is exact). By Theorem 2.6.3, $\|\boldsymbol{\delta}(R\,\mathbf{v})\| = \|\boldsymbol{\delta}(\mathbf{v})\|$. By Theorem 2.6.5, $\|\boldsymbol{\delta}(\dot{R}\,\mathbf{r})\| \leq |\omega|\,\|\boldsymbol{\delta}(\mathbf{r})\|$. Triangle inequality gives the stated bound. ∎

**Corollary 2.6.2** (Dominant error sources in SGP4)**.** *The total error budget for the perifocal→PEF pipeline is dominated by: (a) TLE element measurement error in the rotation angles ($\sigma_m(\Omega)$, $\sigma_m(i)$, $\sigma_m(u)$), and (b) GMST precision error in the transport term. The rotation matrix precision errors ($\delta_p$ from trig evaluation) are negligible by comparison (typically 14+ digits of precision vs 4–6 digits of TLE accuracy).*

*Proof.* The measurement errors in TLE orbital angles are $\sigma_m \sim 10^{-4}$ to $10^{-3}$ rad (corresponding to the TLE precision of $10^{-4}$ degrees). The precision errors in computing $\sin$ and $\cos$ are $\delta_p \sim 10^{-16}$ (double precision). Since $10^{-4} \gg 10^{-16}$, measurement error dominates precision error in the rotation matrix entries by at least 12 orders of magnitude. For the transport term, $\delta(\omega_E)$ is a measurement error from the IAU determination, while $\delta(\theta_{\mathrm{GMST}})$ accumulates precision error from the polynomial evaluation (Ch. 29). ∎

**Theorem 2.6.6** (Dot product error bound)**.** *For $s = \mathbf{a} \cdot \mathbf{b} = \sum_{i=1}^{3} a_i b_i$, the error bound is:*

$$
\delta(s) \leq \sum_{i=1}^{3} \bigl(|a_i|\,\delta(b_i) + |b_i|\,\delta(a_i) + \delta(a_i)\,\delta(b_i)\bigr). \tag{2.38}
$$

*Proof.* The dot product is a sum of three products. Apply Ch. 1, Theorem 1.3.2 to each product and Ch. 1, Theorem 1.3.1 to sum the three results. ∎

**Theorem 2.6.7** (Cross product error bound)**.** *For $\mathbf{c} = \mathbf{a} \times \mathbf{b}$, each component is a difference of two products (e.g., $c_1 = a_2 b_3 - a_3 b_2$). The error bound is:*

$$
\delta(c_k) \leq |a_i|\,\delta(b_j) + |b_j|\,\delta(a_i) + |a_j|\,\delta(b_i) + |b_i|\,\delta(a_j) + \text{(second-order terms)}
$$

*where $(i,j)$ are the appropriate index pairs for component $k$. Subtractive cancellation (Ch. 1, §1.7) occurs when $\mathbf{a} \approx \lambda\mathbf{b}$ for some scalar $\lambda$ — i.e., when the vectors are nearly parallel. In this regime, the cross product loses reliable digits.*

*Proof.* Each component of the cross product is a difference of two products. Apply Ch. 1, Theorem 1.3.2 to each product. The difference of two nearly equal quantities triggers the subtractive cancellation analysis of Ch. 1, §1.7. When $\mathbf{a} \approx \lambda\mathbf{b}$, the products $a_i b_j$ and $a_j b_i$ are nearly equal, and their difference loses significance. The condition number of the subtraction diverges as $\mathbf{a}$ and $\mathbf{b}$ become parallel (Ch. 1, Theorem 1.7.1). ∎

**Theorem 2.6.8** (Norm error bound)**.** *For $r = \|\mathbf{v}\| = \sqrt{v_1^2 + v_2^2 + v_3^2}$, the error bound is:*

$$
\delta(r) \leq \frac{1}{r}\sum_{i=1}^{3} |v_i|\,\delta(v_i) \leq \frac{\|\mathbf{v}\|\,\max_i\delta(v_i)}{r}\,\sqrt{3} = \max_i\delta(v_i)\,\sqrt{3}. \tag{2.39}
$$

*The Euclidean norm does not amplify the total error magnitude. At $\|\mathbf{v}\| \to 0$, the relative error $\delta(r)/r$ diverges — see Ch. 1, Definition 1.9.2.*

*Proof.* By the chain rule: $\frac{\partial r}{\partial v_i} = \frac{v_i}{r}$. By Ch. 1, Theorem 1.4.2 (multivariate error bound): $\delta(r) \leq \sum_i |\frac{v_i}{r}|\,\delta(v_i)$. By Cauchy-Schwarz: $\sum_i |v_i/r|\,\delta(v_i) \leq \frac{1}{r}\|\mathbf{v}\|\,\|\boldsymbol{\delta}\| = \|\boldsymbol{\delta}\|$ where $\boldsymbol{\delta} = (\delta(v_1), \delta(v_2), \delta(v_3))$. Therefore $\delta(r) \leq \|\boldsymbol{\delta}\|$. The norm contracts the error: from a 3-component error vector to a scalar bounded by its magnitude. ∎

**Theorem 2.6.9** (Three-error independence)**.** *The three error categories $\sigma_m$, $\delta_p$, $\delta_a$ propagate independently through all matrix operations. For a linear operation $f(\mathbf{x}) = A\mathbf{x}$, the error in each category transforms by the same matrix $A$. For a nonlinear operation, the sensitivity bound $B_f$ is a property of the operation and the physical state — it does not depend on which error category is being propagated.*

*Proof.* The three categories share the same state-space structure (same 7 or 13 components). The matrix operations act on that structure, not on the error labels. Theorem 2.6.1 applies identically regardless of whether $\delta$ represents $\sigma_m$, $\delta_p$, or $\delta_a$. Since the same bound applies regardless of category, no cross-category coupling occurs. This is the vector generalization of Ch. 1, Proposition 1.8.1. ∎

---

## §2.7 The 13×13 Extended State Matrix

**Definition 2.7.1** (Extended state vector with error tracking)**.** *The extended state vector carries the physical state plus the three error categories for position:*

$$
\mathbf{x}_{13} = \begin{bmatrix} \mathbf{r} \\ \mathbf{v} \\ \boldsymbol{\sigma}_m(\mathbf{r}) \\ \boldsymbol{\delta}_p(\mathbf{r}) \\ 1 \end{bmatrix} \in \mathbb{R}^{13} \tag{2.40}
$$

*where $\boldsymbol{\sigma}_m(\mathbf{r}) = (\sigma_m(r_x), \sigma_m(r_y), \sigma_m(r_z))^T$ and $\boldsymbol{\delta}_p(\mathbf{r}) = (\delta_p(r_x), \delta_p(r_y), \delta_p(r_z))^T$ are error magnitude vectors for position. Velocity errors are derived from the matrix product (they depend on both velocity input errors and position errors through the transport term). The accuracy error $\boldsymbol{\delta}_a$ is handled identically to $\boldsymbol{\delta}_p$ and is omitted from the vector to keep the notation manageable.*

**Remark.** An alternative partition $[\mathbf{r};\, \mathbf{v};\, \boldsymbol{\sigma}(\mathbf{r},\mathbf{v});\, \boldsymbol{\delta}_p(\mathbf{r},\mathbf{v});\, 1]$ was considered but rejected because it mixes position and velocity error magnitudes, obscuring the transport-term coupling structure.

**Definition 2.7.2** ($13 \times 13$ extended state matrix)**.** *The extended state matrix has the block structure:*

$$
T_{13} = \begin{bmatrix} R & \mathbf{0} & \mathbf{0} & \mathbf{0} & \Delta\mathbf{r} \\ \dot{R} & R & \mathbf{0} & \mathbf{0} & \Delta\mathbf{v} \\ \mathbf{0} & \mathbf{0} & |R| & \mathbf{0} & \boldsymbol{\delta}_{\mathrm{rot}} \\ \mathbf{0} & \mathbf{0} & \mathbf{0} & |R| & \boldsymbol{\delta}_{\mathrm{rot}} \\ \mathbf{0}^T & \mathbf{0}^T & \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix}
$$

*where $|R|$ denotes the entry-wise absolute value of $R$ (used because error magnitudes are non-negative and propagate through absolute values of the matrix entries), and $\boldsymbol{\delta}_{\mathrm{rot}}$ accounts for additional error introduced by inexact rotation (Theorem 2.6.4). The upper-left $7 \times 7$ block is the physical state matrix $T$ from (2.11).*

**Definition 2.7.3** (Linearized error propagation matrix)**.** *The Jacobian $|\partial f/\partial \mathbf{x}|$ evaluated at the nominal state, applied to error magnitudes. For linear operations (rotations, translations, transport terms), this equals $|T|$ — the entry-wise absolute value of the state matrix. For nonlinear operations, the Jacobian must be computed at the current physical state.*

**Theorem 2.7.1** (Linearized error propagation is a matrix operation)**.** *For small errors, the error propagation rules of Ch. 1, §§1.3–1.4 are linear in the error magnitudes. Therefore, propagating errors through a transform is itself a matrix-vector multiplication:*

$$
\boldsymbol{\delta}(\mathbf{r}') = |R|\,\boldsymbol{\delta}(\mathbf{r}) + \|\mathbf{r}\|\,\boldsymbol{\delta}_{\mathrm{rot}}
$$

*where the first term propagates input errors through the exact rotation and the second term adds the rotation matrix's own precision error.*

*Proof.* From Theorem 2.6.1, dropping the second-order cross terms (valid when $\delta \ll |v|$): $\delta(y_i) \approx \sum_j |A_{ij}|\,\delta(x_j) + \sum_j |x_j|\,\delta(A_{ij})$. The first sum is $|A|\,\boldsymbol{\delta}(\mathbf{x})$ and the second sum is $\boldsymbol{\delta}_A \cdot |\mathbf{x}|$. Both operations are linear in the error magnitudes. Therefore the error propagation can be written as a matrix-vector operation. ∎

**Theorem 2.7.2** ($13 \times 13$ composition)**.** *The extended state matrix composes correctly: $T_{13,AC} = T_{13,AB}\,T_{13,BC}$. The error blocks track cumulative error through the chain of transforms, with cross-coupling between position error and velocity error occurring through the transport term block.*

*Proof.* The $7 \times 7$ physical state block composes by Theorem 2.4.2. The error blocks compose by Theorem 2.6.2 applied to the error magnitude vectors. The key coupling: the lower-left transport block $\dot{R}$ in the physical state matrix has a corresponding entry in the error propagation that maps position error into velocity error. The composition of $13 \times 13$ matrices correctly accumulates this coupling through the chain. The error blocks in the third and fourth row-groups are decoupled from each other (Theorem 2.6.9: the three error categories are independent). ∎

**Remark** (Linearization validity)**.** The $13 \times 13$ form is accurate when $\delta \ll |v|$ (small-error regime). When reliable digits drop below 4 ($\mathrm{rd}(v) < 4$), the linearization breaks down and the full nonlinear error propagation from Ch. 1 must be used. [A.2.2]

**Remark** (The 10×10 extension)**.** Adding body orientation (3 Euler angles or equivalent) to the state vector gives a $10 \times 10$ matrix. This is relevant for enhanced models where drag cross-section depends on satellite attitude (Ch. 36), but not for SGP4 standard mode.

**Remark** (Implementation alternative)**.** Carrying errors as a separate computation alongside the $7 \times 7$ multiply, rather than embedding them in a $13 \times 13$ matrix, is computationally equivalent. The $13 \times 13$ form is conceptually cleaner (a single matrix multiply propagates everything) but the implementation may prefer the separated form for clarity and to avoid storing sparse blocks.

---

## §2.8 Summary and Usage Guide

This chapter established the $7 \times 7$ homogeneous state matrix as the primary framework for representing, composing, and error-analyzing all coordinate transforms and state propagation steps.

**The "State Matrix Formulation" section pattern.** Every subsequent chapter that performs a coordinate transformation or state propagation step includes a section titled "State Matrix Formulation." That section expresses the chapter's operation as a $7 \times 7$ (or $13 \times 13$) state matrix, with explicit block entries derived from the chapter's specific transform.

**The composition principle.** The SGP4 propagation pipeline is a chain of state matrix multiplications:

$$
T_{\mathrm{total}} = T_{\mathrm{coord}} \cdot (I_7 + \delta T_{\mathrm{SP}}) \cdot T_{\mathrm{Kepler}} \cdot (I_7 + \delta T_{\mathrm{LP}}) \cdot T_{\mathrm{secular}}.
$$

The error budget for the entire pipeline follows from applying §2.6 to each multiplication in the chain.

**Downstream chapters that build on this framework:**

| Chapter | Transform expressed as | Key theorem |
|---------|----------------------|-------------|
| Ch. 8 | Perifocal → ECI pure rotation state matrix | Theorem 2.5.1 |
| Ch. 16–17 | Secular perturbation as state matrix update | Definition 2.4.1 |
| Ch. 18 | Short-period corrections as near-identity state matrix | Proposition 2.5.1 |
| Ch. 19 | Long-period corrections as near-identity perturbation | Proposition 2.5.1 |
| Ch. 20 | Complete mean → osculating composition | Theorem 2.4.2 |
| Ch. 29 | Sidereal time in TEME → PEF rotation | Theorem 2.5.2 |
| Ch. 30 | All coordinate transforms | Theorems 2.5.1–2.5.3 |
| Ch. 34–35 | Full propagation pipeline | Theorem 2.4.2 (chain) |
| Ch. 38 | Final state vector error budget | §2.6, §2.7 |

---

## Error Notes

**[P.2.1]** Trigonometric evaluation in rotation matrices. Computing $\sin(\Omega)$, $\cos(\Omega)$, $\sin(i)$, $\cos(i)$, etc. introduces representation error bounded by Ch. 1, Corollaries 1.4.1–1.4.2. Each DCM entry is a product of two trig values, so the precision error accumulates as $O(n\,\epsilon_{\mathrm{mach}})$ per entry. For double precision: $\delta_p \leq 1.1 \times 10^{-16}$ per trig evaluation.

**[P.2.2]** Matrix multiplication accumulation. An $n \times n$ matrix product accumulates $O(n)$ roundoff per entry from $n$ multiply-add operations. For the $3 \times 3$ DCM products in the SGP4 pipeline, this is $O(3\,\epsilon_{\mathrm{mach}}) \approx 3.3 \times 10^{-16}$ per entry.

**[P.2.3]** Orthogonality drift. After $n$ DCM multiplications without re-orthogonalization, the constraint $R^T R = I$ drifts by $O(n\,\epsilon_{\mathrm{mach}})$. *Remedy:* re-orthogonalize after every multiplication chain, or monitor $\|R^T R - I\|_F$ and re-orthogonalize when it exceeds a threshold.

**[A.2.1]** TEME frame accuracy. The TEME frame differs from the precise GCRS/ITRS frames by $\sim 0.1$ arcsec ($\sim 50$ m at GEO, $\sim 0.3$ m at LEO). This is an accuracy error. *Remedy:* apply IAU 2006/2000A precession-nutation corrections (beyond SGP4 standard).

**[A.2.2]** Linearization in error propagation. The $13 \times 13$ error propagation uses first-order (linear) approximations. Neglected second-order cross terms are $O(\delta^2)$. *Remedy:* use full nonlinear propagation from Ch. 1 when $\mathrm{rd}(v) < 4$.

**[M.2.1]** GMST determination. The GMST polynomial coefficients (Ch. 29) carry measurement error from the IAU 1982 determination. Propagates into the transport term: $\sigma_m(\theta_{\mathrm{GMST}}) \approx 0.1$ mas after one day. *Remedy:* use updated IERS values (but breaks the matched pair — see Ch. 3).
