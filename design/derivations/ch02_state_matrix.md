# Chapter 2: The State Framework

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| **SU(2) algebra** | | |
| $\sigma_k$ | Pauli matrices ($k = 1, 2, 3$) | §2.2, Def. 2.2.1 |
| $\boldsymbol{\sigma}$ | The Pauli vector $(\sigma_1, \sigma_2, \sigma_3)$ | §2.2 |
| $I$ | The $2 \times 2$ identity matrix | §2.2 |
| $A^\dagger$ | Conjugate transpose (Definition 2.2.2) | §2.2 |
| $M$ | An element of $SU(2)$: unitary with $\det M = 1$ | §2.2, Def. 2.2.6 |
| $\alpha$, $\beta$ | Complex parameters of $M$ satisfying $\lvert\alpha\rvert^2+\lvert\beta\rvert^2=1$ | §2.2, Def. 2.2.6 |
| $(w, x, y, z)$ | Quaternion components: $\alpha = w + iz$, $\beta = ix + y$ | §2.2, Def. 2.2.6 |
| $V$ | A traceless Hermitian $2 \times 2$ matrix encoding a 3D vector | §2.2, Def. 2.2.5 |
| $\delta_{jk}$ | Kronecker delta: $1$ if $j=k$, else $0$ | §2.2, Def. 2.2.9 |
| $\varepsilon_{jkl}$ | Levi-Civita symbol: fully antisymmetric tensor on $\{1,2,3\}$ | §2.2, Def. 2.2.9 |
| $[A, B]$ | Matrix commutator $AB - BA$ | §2.2 |
| $M_k(\alpha)$ | Elementary rotation about axis $k$ by angle $\alpha$ | §2.2, Def. 2.2.8 |
| $\hat{\mathbf{e}}_k$ | Standard unit basis vector along axis $k$ | §2.2 |
| **Dual extension** | | |
| $\hat{M}$ | A dual $SU(2)$ matrix: $\hat{M} = M + \varepsilon D$ | §2.3, Def. 2.3.2 |
| $\varepsilon$ | Dual unit satisfying $\varepsilon^2 = 0$ (distinct from $\varepsilon_{jkl}$) | §2.3, Def. 2.3.1 |
| $D$ | Dual part of $\hat{M}$, encoding translation via $D = \frac{1}{2}TM$ | §2.3, Def. 2.3.4 |
| $T$ | Translation encoded as traceless Hermitian matrix: $T = t_k \sigma_k$ | §2.3, Def. 2.3.4 |
| **Rigid body state** | | |
| $R_{\mathrm{pos}}$ | Position encoded as traceless Hermitian matrix | §2.4, Def. 2.4.1 |
| $V_{\mathrm{vel}}$ | Velocity encoded as traceless Hermitian matrix | §2.4, Def. 2.4.1 |
| $\Omega$ | Angular velocity encoded as traceless Hermitian matrix | §2.4, Def. 2.4.1 |
| $\Omega_T$ | Angular velocity of the target reference frame | §2.4, Thm. 2.4.1 |
| **Orbital elements** (introduced in §2.5, defined fully in Ch 8) | | |
| $\Omega_{\mathrm{node}}$ | Right ascension of the ascending node | §2.5 |
| $i$ | Inclination (also used as $\sqrt{-1}$; context distinguishes) | §2.5 |
| $u = \omega + \nu$ | Argument of latitude ($\omega$: argument of perigee, $\nu$: true anomaly) | §2.5 |
| $e$ | Orbital eccentricity | §2.5 |
| $r$ | Orbital radius (scalar, distinct from $\mathbf{r}$ position vector) | §2.5 |
| **Reference frames** | | |
| $\omega_E$ | Earth's rotation rate ($\approx 7.292 \times 10^{-5}$ rad/s) | §2.5, Thm. 2.5.2 |
| $\theta_{\mathrm{GMST}}$ | Greenwich Mean Sidereal Time angle (Ch 29) | §2.5, Def. 2.5.3 |
| **Matrix forms** (derived from SU(2) in §§2.7–2.8) | | |
| $R(M)$ | The $3 \times 3$ rotation matrix from $M$ via the adjoint representation | §2.7, Def. 2.7.1 |
| $H$ | A $4 \times 4$ homogeneous matrix | §2.8, Def. 2.8.1 |
| $w_H$ | Homogeneous coordinate: $w_H = 1$ for points, $w_H = 0$ for displacements | §2.8, Def. 2.8.1 |
| $T_7$ | The $7 \times 7$ state matrix | §2.8, Def. 2.8.3 |
| $[\mathbf{a}]_\times$ | Skew-symmetric matrix: $[\mathbf{a}]_\times \mathbf{b} = \mathbf{a} \times \mathbf{b}$ | §2.8, Def. 2.8.2 |
| $\Delta\mathbf{v}$ | Velocity offset in the $7 \times 7$ state matrix (zero for standard transforms) | §2.8, Def. 2.8.3 |
| **Error framework** (defined in Ch 1, used throughout) | | |
| $\delta(\cdot)$ | Error bound on a quantity (generic; applies to each category) | Ch 1 |
| $\sigma_m$, $\delta_p$, $\delta_a$ | Measurement, precision, accuracy error bounds | Ch 1, Defs. 1.2.1–1.2.3 |
| $\epsilon_{\mathrm{mach}}$ | Machine epsilon of the arithmetic type ($2^{-52}$ for binary64) | Ch 1, Def. 1.2.2 |
| $\mathrm{rd}(v)$ | Reliable digits (Ch 1, Def. 1.9.1) | Ch 1 |

---

## §2.1 Introduction

Orbit propagation produces a satellite's position and velocity in one reference frame. Comparing with observations, planning maneuvers, or computing ground tracks requires transforming this state into other frames. Each transformation involves rotation (the frames are oriented differently), translation (the frames may have different origins), and — for rotating frames — velocity coupling (the transport theorem adds $\boldsymbol{\omega}_T \times \mathbf{r}$ to the observed velocity, where $\boldsymbol{\omega}_T$ is the frame's angular velocity).

The traditional approach handles these separately: a $3 \times 3$ rotation matrix rotates vectors, a translation vector is added, and the transport term is computed from a cross product. Composing multiple transforms requires tracking which matrix, which offset, and which angular velocity applies at each stage.

This chapter develops a unified algebraic framework in which a single mathematical object encodes rotation, translation, and velocity coupling simultaneously. The framework has three properties that make it the foundation for all subsequent chapters:

1. **Singularity-free.** The representation of orientation has no coordinate singularity (no gimbal lock), eliminating the need for the Lyddane modification and similar patches at special inclinations or eccentricities.

2. **Compositional.** Two transforms compose by a single algebraic product, producing a third transform of the same type. The entire propagation pipeline — from orbital-plane coordinates through inertial frame to Earth-fixed frame — is a chain of such products.

3. **Error-compatible.** The transform is linear in the state, so the error state transforms by the same operation as the physical state (§2.6, Theorem 2.6.1). The composed error bounds remain valid through chains of operations (Ch 1, Theorem 1.3.4). The three error categories ($\sigma_m$, $\delta_p$, $\delta_a$) propagate independently through the same transform (Ch 1, Proposition 1.8.1).

The algebraic foundation is the group $SU(2)$ of $2 \times 2$ unitary matrices with unit determinant, together with its dual extension for incorporating translation (§2.3). Every 3D vector — position, velocity, angular velocity — is encoded as a $2 \times 2$ traceless Hermitian matrix (Definition 2.2.5), and every rotation acts by conjugation: $V' = MVM^\dagger$ (Definition 2.2.7). The cross product, dot product, and norm are native operations in this algebra (Theorems 2.2.1–2.2.3), so the transport theorem and error propagation emerge from standard matrix calculus rather than vector-specific formulas.

The traditional $3 \times 3$ rotation matrix and $7 \times 7$ state matrix are derived in §§2.7–2.9 as component-wise expansions of the $SU(2)$ operations. These provide the bridge to the existing orbit propagation literature and a linear operator form suitable for SVD and matrix analysis. The conversion between representations is exact and bidirectional (§2.7 forward, Theorem 2.9.1 inverse), so the framework connects to existing DCM-based implementations without loss.

---

## §2.2 The Algebra of Rotation

### Foundations

Three-dimensional rotations can be represented by $2 \times 2$ complex matrices. This section builds the representation from its basic elements: the Pauli basis, the encoding of vectors, the definition of rotation, and the algebraic identities that connect matrix operations to vector operations.

**Definition 2.2.1** (Pauli matrices)**.** *The three Pauli matrices are:*

$$\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} \tag{2.1}$$

*Together with the identity $I$, they form a basis for the real vector space of $2 \times 2$ Hermitian matrices (Definition 2.2.3).*

**Definition 2.2.2** (Conjugate transpose)**.** *For a matrix $A \in \mathbb{C}^{m \times n}$ with entries $a_{jk}$, the conjugate transpose $A^\dagger \in \mathbb{C}^{n \times m}$ is defined by:*

$$(A^\dagger)_{jk} = \overline{a_{kj}} \tag{2.2}$$

*That is: transpose the matrix and take the complex conjugate of every entry. For a $2 \times 2$ matrix:*

$$\begin{pmatrix} a & b \\ c & d \end{pmatrix}^\dagger = \begin{pmatrix} \bar{a} & \bar{c} \\ \bar{b} & \bar{d} \end{pmatrix}$$

*Key properties: $(AB)^\dagger = B^\dagger A^\dagger$ (reversal of order), $(A^\dagger)^\dagger = A$, and $\det(A^\dagger) = \overline{\det(A)}$. For real matrices, $A^\dagger = A^T$ (ordinary transpose).*

**Definition 2.2.3** (Hermitian matrix)**.** *A matrix $A \in \mathbb{C}^{2 \times 2}$ is Hermitian if $A = A^\dagger$. The constraint $A = A^\dagger$ forces the diagonal entries to satisfy $a_{jj} = \overline{a_{jj}}$, which requires $\mathrm{Im}(a_{jj}) = 0$: the diagonal entries are real. The off-diagonal entries satisfy $a_{12} = \overline{a_{21}}$. A general $2 \times 2$ Hermitian matrix therefore has the form:*

$$A = \begin{pmatrix} a & b \\ \bar{b} & c \end{pmatrix} \tag{2.3}$$

*where $a, c \in \mathbb{R}$ and $b \in \mathbb{C}$. This is a 4-dimensional real vector space.*

**Definition 2.2.4** (Traceless Hermitian matrix)**.** *A Hermitian matrix $A$ is traceless if $\mathrm{tr}(A) = a + c = 0$, i.e., $c = -a$. Every traceless Hermitian $2 \times 2$ matrix has the form:*

$$A = \begin{pmatrix} a & b \\ \bar{b} & -a \end{pmatrix} \tag{2.4}$$

*where $a \in \mathbb{R}$ (forced by the Hermitian constraint) and $b \in \mathbb{C}$. This gives 3 real degrees of freedom ($a$, $\mathrm{Re}(b)$, $\mathrm{Im}(b)$), forming a 3-dimensional real vector space with basis $\{\sigma_1, \sigma_2, \sigma_3\}$.*

**Definition 2.2.5** (Vector encoding)**.** *A vector $\mathbf{v} = (v_1, v_2, v_3) \in \mathbb{R}^3$ is encoded as the traceless Hermitian matrix:*

$$V = v_1 \sigma_1 + v_2 \sigma_2 + v_3 \sigma_3 = \begin{pmatrix} v_3 & v_1 - iv_2 \\ v_1 + iv_2 & -v_3 \end{pmatrix} \tag{2.5}$$

*This is an isomorphism: every traceless Hermitian $2 \times 2$ matrix corresponds to exactly one vector in $\mathbb{R}^3$, and vice versa. The components are recovered by $v_k = \frac{1}{2}\mathrm{tr}(\sigma_k V)$.*

**Example 2.2.1** (Vector encoding)**.** *The vector $\mathbf{v} = (3, 4, 0)$ encodes as:*

$$V = 3\sigma_1 + 4\sigma_2 + 0\cdot\sigma_3 = \begin{pmatrix} 0 & 3 - 4i \\ 3 + 4i & 0 \end{pmatrix}$$

*The norm: $\|\mathbf{v}\|^2 = -\det(V) = -(0 \cdot 0 - (3-4i)(3+4i)) = -(-(9+16)) = 25$, so $\|\mathbf{v}\| = 5$. Recovery: $v_1 = \frac{1}{2}\mathrm{tr}(\sigma_1 V) = \frac{1}{2}\mathrm{tr}\left(\begin{smallmatrix} 3+4i & 0 \\ 0 & 3-4i \end{smallmatrix}\right) = \frac{1}{2}(6) = 3$. ✓*

**Definition 2.2.6** (Rotation element)**.** *A rotation is an element $M \in SU(2)$: a $2 \times 2$ unitary matrix with unit determinant:*

$$M^\dagger M = I, \qquad \det M = 1 \tag{2.6}$$

*The general form, parameterized by $\alpha, \beta \in \mathbb{C}$ with $|\alpha|^2 + |\beta|^2 = 1$:*

$$M = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix} \tag{2.7}$$

*This is equivalent to the unit quaternion $q = (w, x, y, z)$ via the identification $\alpha = w + iz$, $\beta = ix + y$.*

**Assumption 2.2.1** (Right-handed coordinates)**.** *All coordinate systems are right-handed throughout. Rotation by angle $\theta$ about axis $\hat{\mathbf{e}}$ follows the right-hand rule.*

**Definition 2.2.7** (Conjugation)**.** *The conjugation of a matrix $V$ by an invertible matrix $M$ is the operation $V \mapsto MVM^\dagger$. When $M \in SU(2)$ and $V$ is a traceless Hermitian matrix encoding a vector $\mathbf{v}$, the conjugation produces another traceless Hermitian matrix encoding the rotated vector (Theorems 2.2.4–2.2.5):*

$$V' = M\,V\,M^\dagger \tag{2.8}$$

**Definition 2.2.8** (Elementary rotations)**.** *Rotation by angle $\alpha$ about coordinate axis $\hat{\mathbf{e}}_k$ is:*

$$M_k(\alpha) = \cos\frac{\alpha}{2}\,I + i\sin\frac{\alpha}{2}\,\sigma_k = \exp\!\left(\frac{i\alpha}{2}\,\sigma_k\right) \tag{2.9}$$

*The exponential form is the exponential map from the Lie algebra $\mathfrak{su}(2)$ (spanned by $\{i\sigma_1/2, i\sigma_2/2, i\sigma_3/2\}$) to the Lie group $SU(2)$. The half-angle $\alpha/2$ arises because $SU(2)$ is a double cover of $SO(3)$ (Lemma 2.2.2). The three elementary rotations are:*

$$M_1(\alpha) = \begin{pmatrix} \cos\frac{\alpha}{2} & i\sin\frac{\alpha}{2} \\ i\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}, \quad M_2(\alpha) = \begin{pmatrix} \cos\frac{\alpha}{2} & \sin\frac{\alpha}{2} \\ -\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}, \quad M_3(\alpha) = \begin{pmatrix} e^{i\alpha/2} & 0 \\ 0 & e^{-i\alpha/2} \end{pmatrix} \tag{2.10}$$

*Note: $M_3(\alpha)$ is diagonal — rotation about $\hat{\mathbf{e}}_3$ reduces to a phase multiplication.*

**Example 2.2.2** (Elementary rotation matrix)**.** *A $90°$ rotation about $\hat{\mathbf{e}}_3$ ($\alpha = \pi/2$, so $\alpha/2 = \pi/4$):*

$$M_3(90°) = \begin{pmatrix} e^{i\pi/4} & 0 \\ 0 & e^{-i\pi/4} \end{pmatrix} = \begin{pmatrix} \frac{1+i}{\sqrt{2}} & 0 \\ 0 & \frac{1-i}{\sqrt{2}} \end{pmatrix}$$

*Verification: $|\alpha|^2 + |\beta|^2 = \frac{1}{2}(1+1) + 0 = 1$. ✓ The determinant is $e^{i\pi/4} \cdot e^{-i\pi/4} = 1$. ✓*

### Algebraic Properties

The Pauli matrices satisfy a set of multiplication identities that connect the $2 \times 2$ matrix algebra to the familiar vector operations of $\mathbb{R}^3$. The cross product, dot product, and norm each correspond to a standard matrix operation. These identities are used throughout the chapter: in §2.3 for the dual extension, in §2.4 for the transport theorem, and in §2.7 for deriving the $3 \times 3$ rotation matrix.

**Definition 2.2.9** (Kronecker delta and Levi-Civita symbol)**.** *The Kronecker delta $\delta_{jk}$ equals $1$ when $j = k$ and $0$ otherwise. The Levi-Civita symbol $\varepsilon_{jkl}$ is the fully antisymmetric tensor on indices $\{1, 2, 3\}$:*

$$\varepsilon_{123} = \varepsilon_{231} = \varepsilon_{312} = +1, \qquad \varepsilon_{132} = \varepsilon_{213} = \varepsilon_{321} = -1 \tag{2.11}$$

*All other values (any repeated index) are zero. The Levi-Civita symbol encodes the cross product: $(\mathbf{a} \times \mathbf{b})_l = \sum_{jk}\varepsilon_{jkl}\,a_j\,b_k$.*

**Lemma 2.2.1** (Pauli algebra)**.** *The Pauli matrices satisfy:*

$$\sigma_j \sigma_k = \delta_{jk}\,I + i\,\varepsilon_{jkl}\,\sigma_l \tag{2.12}$$

*where summation over the repeated index $l$ is implied. Equivalently:*
- *Anticommutator: $\{\sigma_j, \sigma_k\} = \sigma_j\sigma_k + \sigma_k\sigma_j = 2\delta_{jk}\,I$*
- *Commutator: $[\sigma_j, \sigma_k] = \sigma_j\sigma_k - \sigma_k\sigma_j = 2i\,\varepsilon_{jkl}\,\sigma_l$*
- *Trace: $\mathrm{tr}(\sigma_j\sigma_k) = 2\delta_{jk}$*

*Proof.* The complete multiplication table, computed from (2.1):

| | $\sigma_1$ | $\sigma_2$ | $\sigma_3$ |
|---|---|---|---|
| $\sigma_1$ | $I$ | $i\sigma_3$ | $-i\sigma_2$ |
| $\sigma_2$ | $-i\sigma_3$ | $I$ | $i\sigma_1$ |
| $\sigma_3$ | $i\sigma_2$ | $-i\sigma_1$ | $I$ |

Each diagonal entry ($j = k$): $\sigma_k^2 = I$, confirming $\delta_{kk}I = I$ with no $\varepsilon$ term (since $\varepsilon_{kkl} = 0$). Each off-diagonal entry ($j \neq k$): $\sigma_j\sigma_k = i\varepsilon_{jkl}\sigma_l$, verified for all six cases. For example: $\sigma_1\sigma_2 = \left(\begin{smallmatrix} 0 & 1 \\ 1 & 0 \end{smallmatrix}\right)\left(\begin{smallmatrix} 0 & -i \\ i & 0 \end{smallmatrix}\right) = \left(\begin{smallmatrix} i & 0 \\ 0 & -i \end{smallmatrix}\right) = i\sigma_3$, matching $i\varepsilon_{123}\sigma_3 = i\sigma_3$. ∎

The three vector operations — cross product, dot product, and norm — each arise from applying a standard matrix operation to the encoded vectors.

**Theorem 2.2.1** (Cross product as commutator)**.** *For vectors $\mathbf{a}, \mathbf{b} \in \mathbb{R}^3$ encoded as traceless Hermitian matrices $A, B$:*

$$\mathbf{a} \times \mathbf{b} \;\longleftrightarrow\; \frac{1}{2i}[A, B] \tag{2.13}$$

*where $[A, B] = AB - BA$ is the matrix commutator.*

*Proof.* Expand $A = \sum_j a_j\sigma_j$, $B = \sum_k b_k\sigma_k$. The commutator is:

$$[A, B] = \sum_{j,k} a_j b_k [\sigma_j, \sigma_k] = \sum_{j,k} a_j b_k \cdot 2i\varepsilon_{jkl}\sigma_l = 2i\sum_l \left(\sum_{j,k}\varepsilon_{jkl}\,a_j b_k\right)\sigma_l = 2i\,(\mathbf{a} \times \mathbf{b})\cdot\boldsymbol{\sigma}$$

Dividing by $2i$: $\frac{1}{2i}[A, B] = (\mathbf{a} \times \mathbf{b})\cdot\boldsymbol{\sigma}$, which is the traceless Hermitian matrix encoding $\mathbf{a} \times \mathbf{b}$. ∎

**Theorem 2.2.2** (Dot product from trace)**.** *For vectors encoded as traceless Hermitian matrices:*

$$\mathbf{a} \cdot \mathbf{b} = \frac{1}{2}\mathrm{tr}(AB) \tag{2.14}$$

*Proof.* $\mathrm{tr}(AB) = \sum_{j,k} a_j b_k\,\mathrm{tr}(\sigma_j\sigma_k) = \sum_{j,k} a_j b_k \cdot 2\delta_{jk} = 2\sum_j a_j b_j = 2(\mathbf{a}\cdot\mathbf{b})$. ∎

**Theorem 2.2.3** (Norm from determinant)**.** *The squared norm of a vector is:*

$$\|\mathbf{v}\|^2 = -\det(V) \tag{2.15}$$

*Proof.* $\det(V) = \det\left(\begin{smallmatrix} v_3 & v_1-iv_2 \\ v_1+iv_2 & -v_3 \end{smallmatrix}\right) = -v_3^2 - (v_1^2+v_2^2) = -(v_1^2+v_2^2+v_3^2)$. ∎

**Example 2.2.3** (Cross product via commutator)**.** *Compute $\hat{\mathbf{e}}_1 \times \hat{\mathbf{e}}_2$ using the commutator formula. Encode: $A = \sigma_1$, $B = \sigma_2$. The commutator:*

$$[A, B] = \sigma_1\sigma_2 - \sigma_2\sigma_1 = i\sigma_3 - (-i\sigma_3) = 2i\sigma_3$$

*Therefore $\frac{1}{2i}[A,B] = \sigma_3$, which encodes $\hat{\mathbf{e}}_3$. So $\hat{\mathbf{e}}_1 \times \hat{\mathbf{e}}_2 = \hat{\mathbf{e}}_3$. ✓*

**Example 2.2.4** (Rotation by conjugation)**.** *Rotate $\mathbf{v} = (1, 0, 0)$ by $90°$ about $\hat{\mathbf{e}}_3$ using $M_3(90°)$ from Example 2.2.2. The encoding is $V = \sigma_1 = \left(\begin{smallmatrix} 0 & 1 \\ 1 & 0 \end{smallmatrix}\right)$. The conjugation:*

$$V' = M_3 V M_3^\dagger = \begin{pmatrix} \frac{1+i}{\sqrt{2}} & 0 \\ 0 & \frac{1-i}{\sqrt{2}} \end{pmatrix} \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} \frac{1-i}{\sqrt{2}} & 0 \\ 0 & \frac{1+i}{\sqrt{2}} \end{pmatrix}$$

*Step 1: $M_3 V = \left(\begin{smallmatrix} 0 & \frac{1+i}{\sqrt{2}} \\ \frac{1-i}{\sqrt{2}} & 0 \end{smallmatrix}\right)$. Step 2: $(M_3 V) M_3^\dagger = \left(\begin{smallmatrix} 0 & \frac{(1+i)^2}{2} \\ \frac{(1-i)^2}{2} & 0 \end{smallmatrix}\right) = \left(\begin{smallmatrix} 0 & i \\ -i & 0 \end{smallmatrix}\right)$.*

*Since $(1+i)^2 = 2i$ and $(1-i)^2 = -2i$. Decoding: $V' = \left(\begin{smallmatrix} 0 & i \\ -i & 0 \end{smallmatrix}\right)$. The $(1,2)$ entry is $v_1' - iv_2' = i = 0 - i(-1)$, so $v_1' = 0$, $v_2' = -1$. The diagonal gives $v_3' = 0$. Result: $\mathbf{v}' = (0, -1, 0)$.*

*This is a passive rotation: the vector $(1,0,0)$ expressed in a frame rotated $+90°$ about $z$ points along $-y$. Verification: $R_3(90°)\,(1,0,0)^T = (0,-1,0)^T$. ✓*

### Rotation Properties

Conjugation by an $SU(2)$ element preserves the traceless Hermitian structure of the encoded vector and leaves its norm unchanged. These two properties confirm that conjugation is a rotation — it maps $\mathbb{R}^3$ to $\mathbb{R}^3$ and preserves distances.

**Theorem 2.2.4** (Conjugation preserves vector structure)**.** *If $V$ is traceless Hermitian and $M \in SU(2)$, then $V' = MVM^\dagger$ is traceless Hermitian.*

*Proof.* Hermitian: $(V')^\dagger = (MVM^\dagger)^\dagger = M V^\dagger M^\dagger = MVM^\dagger = V'$ since $V^\dagger = V$. Traceless: $\mathrm{tr}(V') = \mathrm{tr}(MVM^\dagger) = \mathrm{tr}(V) = 0$ by the cyclic property of the trace. ∎

**Theorem 2.2.5** (Conjugation preserves norm)**.** *For $M \in SU(2)$ and traceless Hermitian $V$:*

$$\det(MVM^\dagger) = \det(V) \tag{2.16}$$

*Therefore $\|M\mathbf{v}M^\dagger\| = \|\mathbf{v}\|$: rotation preserves the vector norm.*

*Proof.* $\det(MVM^\dagger) = \det(M)\det(V)\det(M^\dagger) = 1 \cdot \det(V) \cdot 1 = \det(V)$, since $\det(M) = 1$ for $M \in SU(2)$ and $\det(M^\dagger) = \overline{\det(M)} = 1$. ∎

### Composition and Inversion

Rotations compose by $SU(2)$ multiplication and invert by the conjugate transpose. Both properties are inherited directly from the group structure.

**Theorem 2.2.6** (Composition of rotations)**.** *If $M_{AB}$ rotates frame B to frame A, and $M_{BC}$ rotates frame C to frame B, then the composed rotation is:*

$$M_{AC} = M_{AB} \cdot M_{BC} \tag{2.17}$$

*Proof.* For any vector $V$ in frame C: $V_A = M_{AB}(M_{BC}\,V\,M_{BC}^\dagger)M_{AB}^\dagger = (M_{AB}M_{BC})\,V\,(M_{AB}M_{BC})^\dagger$, so $M_{AC} = M_{AB}M_{BC}$. Closure in $SU(2)$: $\det(M_{AC}) = \det(M_{AB})\det(M_{BC}) = 1$, and $(M_{AC})^\dagger M_{AC} = M_{BC}^\dagger M_{AB}^\dagger M_{AB} M_{BC} = I$. ∎

**Theorem 2.2.7** (Inverse rotation)**.** *For $M \in SU(2)$: $M^{-1} = M^\dagger$.*

*Proof.* Immediate from the definition $M^\dagger M = I$. ∎

**Lemma 2.2.2** (Double cover)**.** *$M$ and $-M$ produce the same rotation: $MVM^\dagger = (-M)V(-M)^\dagger$ for all traceless Hermitian $V$. Convention: choose the representative with $\mathrm{Re}(\alpha) > 0$; when $\mathrm{Re}(\alpha) = 0$ (rotation by $\pi$), choose $\mathrm{Im}(\alpha) > 0$; if both are zero, choose $\mathrm{Re}(\beta) > 0$.*

*Proof.* $(-M)V(-M)^\dagger = (-1)^2 MVM^\dagger = MVM^\dagger$. ∎

The elementary rotations (Definition 2.2.8) rotate about coordinate axes. A rotation about an arbitrary axis — needed when the rotation axis does not align with any coordinate direction — is constructed from the exponential map.

**Theorem 2.2.8** (General axis-angle rotation)**.** *A rotation by angle $\theta$ about unit axis $\hat{\mathbf{e}} = (e_1, e_2, e_3)$ is:*

$$M(\theta, \hat{\mathbf{e}}) = \cos\frac{\theta}{2}\,I + i\sin\frac{\theta}{2}\,(e_1\sigma_1 + e_2\sigma_2 + e_3\sigma_3) = \exp\!\left(\frac{i\theta}{2}\,\hat{\mathbf{e}}\cdot\boldsymbol{\sigma}\right) \tag{2.18}$$

*Proof.* Define $E = \hat{\mathbf{e}}\cdot\boldsymbol{\sigma} = e_1\sigma_1 + e_2\sigma_2 + e_3\sigma_3$. The key identity is $E^2 = I$:

$$E^2 = \sum_{j,k} e_j e_k\,\sigma_j\sigma_k = \sum_{j,k} e_j e_k\,(\delta_{jk}I + i\varepsilon_{jkl}\sigma_l) = \left(\sum_j e_j^2\right)I + i\sum_l\left(\sum_{j,k}\varepsilon_{jkl}\,e_j e_k\right)\sigma_l$$

The first sum is $\|\hat{\mathbf{e}}\|^2 = 1$. The second sum vanishes because $\varepsilon_{jkl}$ is antisymmetric while $e_j e_k$ is symmetric in $j, k$. Therefore $E^2 = I$.

With $E^2 = I$, the matrix exponential $\exp(i\frac{\theta}{2}E)$ reduces to the Euler formula for involutions: $e^{i\phi A} = \cos\phi\,I + i\sin\phi\,A$ when $A^2 = I$ (expand the Taylor series and group even/odd powers). Setting $\phi = \theta/2$ and $A = E$ gives (2.18). ∎

### Extracting the Quaternion from a Rotation Matrix

The inverse operation — recovering $M$ from a known $3 \times 3$ rotation matrix $R$ — is needed when converting from $3 \times 3$ DCM representations into the $SU(2)$ framework. The algorithm (Shepperd's method) selects the quaternion component with the largest magnitude to avoid numerical instability. The full development is in §2.9 (Theorem 2.9.1), after the adjoint bridge (§2.7) establishes the forward map $M \to R(M)$.

### Error Properties of Rotation

**Remark** (Matrix norms from scalar bounds)**.** The theorems below use the operator norm $\|A\| = \sup_{\|\mathbf{x}\|=1}\|A\mathbf{x}\|$ for $2 \times 2$ matrices. Ch 1 establishes error propagation for scalar operations. The extension to matrices is entry-by-entry: a matrix-vector product $y_j = \sum_k A_{jk} x_k$ is a sum of products, each bounded by Ch 1, Theorems 1.3.1–1.3.2. The matrix norm then bounds the total: $\|A\mathbf{x}\| \leq \|A\|\,\|\mathbf{x}\|$ by definition, and for unitary matrices $\|M\| = 1$ (all singular values are $1$).

A rotation acts on both the physical state and the error state (§2.6). Whether it amplifies errors depends on whether the rotation itself is exact. When $M$ is known exactly, conjugation is an isometry: it rotates the error vector without changing its magnitude. When $M$ carries error from the evaluation of trigonometric functions, the rotated state acquires an additional error proportional to the state magnitude.

**Theorem 2.2.9** (Isometry — exact rotation preserves error magnitude)**.** *For $M \in SU(2)$ with $\delta(M) = 0$:*

$$\delta(V') = M\,\delta(V)\,M^\dagger, \qquad \|\delta(\mathbf{v}')\| = \|\delta(\mathbf{v})\| \tag{2.19}$$

*Proof.* $V' = MVM^\dagger$ is linear in $V$ for fixed $M$, so $\delta(V') = M\,\delta(V)\,M^\dagger$ exactly (no linearization required). The error $\delta(V)$ is the difference of two traceless Hermitian matrices (the true and computed values), so $\delta(V)$ is itself traceless Hermitian. Norm preservation then follows from Theorem 2.2.5: $\|\delta(\mathbf{v}')\|^2 = -\det(M\,\delta(V)\,M^\dagger) = -\det(\delta(V)) = \|\delta(\mathbf{v})\|^2$. ∎

**Lemma 2.2.3** (Spectral norm of encoded vector)**.** *For a traceless Hermitian $V$ encoding $\mathbf{v}$: $\|V\| = \|\mathbf{v}\|$ where $\|V\|$ is the spectral norm (largest singular value).*

*Proof.* The eigenvalues of $V = \left(\begin{smallmatrix} v_3 & v_1-iv_2 \\ v_1+iv_2 & -v_3 \end{smallmatrix}\right)$ are $\pm\sqrt{v_1^2+v_2^2+v_3^2} = \pm\|\mathbf{v}\|$ (from the characteristic equation $\lambda^2 = v_1^2+v_2^2+v_3^2$). For Hermitian matrices, the spectral norm is the largest absolute eigenvalue. ∎

**Theorem 2.2.10** (Inexact rotation error)**.** *When $M$ carries entry-wise error $\delta(M)$:* [P.2.2]

$$\|\delta(V')\| \leq \|\delta(\mathbf{v})\| + 2\,\|\mathbf{v}\|\,\|\delta(M)\| \tag{2.20}$$

*The first term is the rotated input error (unchanged in magnitude by Theorem 2.2.9). The second term is proportional to the state magnitude — a larger vector suffers more from the same rotation error.*

*Proof.* The perturbed rotation is $M + \delta M$. Expanding the conjugation:

$$(M + \delta M)\,V\,(M + \delta M)^\dagger = MVM^\dagger + (\delta M)\,V\,M^\dagger + M\,V\,(\delta M)^\dagger + (\delta M)\,V\,(\delta M)^\dagger$$

The first term is the exact result $V' = MVM^\dagger$. The last term is $O(\|\delta M\|^2)$ and is dropped for the first-order bound.

Each cross term is bounded by submultiplicativity of the matrix norm. For $(\delta M)\,V\,M^\dagger$:

$$\|(\delta M)\,V\,M^\dagger\| \leq \|\delta M\| \cdot \|V\| \cdot \|M^\dagger\| = \|\delta M\| \cdot \|\mathbf{v}\| \cdot 1$$

where $\|M^\dagger\| = 1$ because $M$ is unitary (the spectral norm of a unitary matrix is $1$), and $\|V\| = \|\mathbf{v}\|$ (the spectral norm of the traceless Hermitian matrix $V$ equals the Euclidean norm of the encoded vector). The same bound holds for $M\,V\,(\delta M)^\dagger$ by the same reasoning. With two cross terms, the total rotation-induced error is $\leq 2\|\mathbf{v}\|\,\|\delta M\|$.

The input error $\delta(V)$ propagates by the isometry of Theorem 2.2.9: $\|\delta(\mathbf{v}')\|_{\mathrm{input}} = \|\delta(\mathbf{v})\|$. Adding both contributions gives (2.20).

The full bound including the second-order term is $\|\delta(V')\| \leq \|\delta(\mathbf{v})\| + 2\|\mathbf{v}\|\,\|\delta(M)\| + \|\mathbf{v}\|\,\|\delta(M)\|^2$. Equation (2.20) is the first-order approximation, valid when $\|\delta(M)\| \ll 1$ (i.e., $\mathrm{rd}(M) \geq 4$, Ch 1 §1.9). ∎

The size of $\delta(M)$ is determined by the error in the input angle, since $M_k(\alpha)$ is constructed from $\cos(\alpha/2)$ and $\sin(\alpha/2)$.

**Theorem 2.2.11** (Angle-to-rotation sensitivity)**.** *The elementary rotation $M_k(\alpha)$ satisfies:*

$$\left\|\frac{dM_k}{d\alpha}\right\| = \frac{1}{2}, \qquad \|\delta(M_k)\| \leq \frac{1}{2}\,|\delta(\alpha)| \tag{2.21}$$

*The rotation error is half the angle error, independent of $\alpha$.* [P.2.3]

*Proof.* Differentiate (2.9): $dM_k/d\alpha = -\frac{1}{2}\sin(\alpha/2)\,I + \frac{i}{2}\cos(\alpha/2)\,\sigma_k$. The eigenvalues of $dM_k/d\alpha$ are $\pm\frac{i}{2}e^{\mp i\alpha/2}$, each with modulus $\frac{1}{2}$. The spectral norm $\|dM_k/d\alpha\| = \frac{1}{2}$. The bound follows from Ch 1, Theorem 1.4.1 with Lipschitz constant $\frac{1}{2}$. ∎

When multiple elementary rotations are composed (as in the perifocal-to-TEME transform of §2.5), each contributes an angle error scaled by $\frac{1}{2}$. By Theorem 2.2.9, the error from each factor passes through subsequent rotations without amplification.

**Corollary 2.2.1** (Composed rotation error)**.** *For the 3-1-3 composition $M = M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u)$:*

$$\|\delta(M)\| \leq \frac{1}{2}\bigl(|\delta(\Omega)| + |\delta(i)| + |\delta(u)|\bigr) + 16\,\epsilon_{\mathrm{mach}} \tag{2.22}$$

*The first three terms are the angle errors, each attenuated by $\frac{1}{2}$ (Theorem 2.2.11). The rounding term $16\,\epsilon_{\mathrm{mach}}$ arises from two $SU(2)$ multiplications, each contributing at most $8\epsilon_{\mathrm{mach}}$ per entry.* [P.2.4]

*Proof.* By Theorem 2.2.11, the error in each factor is $\|\delta(M_k)\| \leq \frac{1}{2}|\delta(\text{angle})|$. The errors combine additively through the two matrix multiplications (Ch 1, Theorem 1.3.2). Each $2 \times 2$ complex matrix multiplication produces 4 output entries, each a sum of 2 complex products. Each complex product involves 4 real multiplications and 2 real additions; with round-to-nearest, the accumulated rounding per entry is bounded by $\leq 4\epsilon_{\mathrm{mach}}$ (Ch 1, Theorem 1.3.2 applied to each real operation). Taking the maximum over entries: $\|\delta(M_{\mathrm{round}})\|_{\max} \leq 8\epsilon_{\mathrm{mach}}$ per composition, and $16\epsilon_{\mathrm{mach}}$ for two compositions. ∎

---

## §2.3 Incorporating Translation: The Dual Extension

A coordinate transform may displace the origin in addition to rotating the axes. Translation is not a linear operation in the $SU(2)$ framework — adding a constant vector $T$ to the conjugation result $MVM^\dagger$ makes the action affine rather than linear. The dual number extension absorbs this addition into the algebra, so that rotation and translation compose by a single product.

### Dual Numbers

**Definition 2.3.1** (Dual numbers)**.** *A dual number is an expression $\hat{a} = a + \varepsilon\,a'$ where $a, a' \in \mathbb{C}$ and $\varepsilon$ is the dual unit satisfying $\varepsilon^2 = 0$, $\varepsilon \neq 0$. Arithmetic follows from this rule:*

$$(\hat{a})(\hat{b}) = (a + \varepsilon a')(b + \varepsilon b') = ab + \varepsilon(ab' + a'b) \tag{2.23}$$

*The dual part captures first-order perturbation: evaluating a polynomial $P$ at a dual number yields $P(a + \varepsilon\,b) = P(a) + \varepsilon\,b\,P'(a)$. This extends to any analytic function via its Taylor series: $f(a + \varepsilon\,b) = f(a) + \varepsilon\,b\,f'(a)$, since $\varepsilon^2 = 0$ kills all higher-order terms.*

### Dual $SU(2)$

**Definition 2.3.2** (Dual $SU(2)$ matrix)**.** *A dual $SU(2)$ matrix is:*

$$\hat{M} = M + \varepsilon\, D \tag{2.24}$$

*where $M \in SU(2)$ is the rotation (the non-dual part) and $D$ is a $2 \times 2$ complex matrix encoding the translation (the dual part). The pair has 8 real parameters (4 from $M$, 4 from $D$).*

**Definition 2.3.3** (Unit dual $SU(2)$ constraints)**.** *A unit dual $SU(2)$ matrix satisfies:*

1. *$M^\dagger M = I$ (the non-dual part is unitary)*
2. *$M^\dagger D + D^\dagger M = 0$ ($D$ lies in the tangent space of $SU(2)$ at $M$)*

*These two constraints reduce the 8 real parameters to 6 degrees of freedom: 3 for rotation and 3 for translation.*

**Definition 2.3.4** (Translation encoding)**.** *Given rotation $M \in SU(2)$ and translation vector $\mathbf{t} = (t_1, t_2, t_3)$, the dual part is:*

$$D = \frac{1}{2}\,T \cdot M \tag{2.25}$$

*where $T = t_1\sigma_1 + t_2\sigma_2 + t_3\sigma_3$ is the translation encoded as a traceless Hermitian matrix (Definition 2.2.4). The factor of $\frac{1}{2}$ is a convention that simplifies the composition rule (Theorem 2.3.2).*

**Theorem 2.3.1** (Translation extraction)**.** *Given a unit dual $SU(2)$ matrix $\hat{M} = M + \varepsilon D$, the translation vector is recovered by:*

$$T = 2\,D\,M^\dagger \tag{2.26}$$

*and the rotation is the non-dual part $M$.*

*Proof.* From (2.25): $D = \frac{1}{2}TM$. Multiply on the right by $M^\dagger$: $DM^\dagger = \frac{1}{2}TMM^\dagger = \frac{1}{2}T$ (since $MM^\dagger = I$). Therefore $T = 2DM^\dagger$. The result $2DM^\dagger$ is traceless Hermitian because it equals the traceless Hermitian matrix $T$ that generated $D$, recovered by the invertibility of right-multiplication by $M$. ∎

### Composition and Inverse

**Theorem 2.3.2** (Dual $SU(2)$ composition)**.** *The product of two dual $SU(2)$ matrices is:*

$$\hat{M}_1 \hat{M}_2 = (M_1 + \varepsilon D_1)(M_2 + \varepsilon D_2) = M_1 M_2 + \varepsilon(M_1 D_2 + D_1 M_2) \tag{2.27}$$

*since $\varepsilon^2 = 0$. The composed rotation is $M_1 M_2$ (Theorem 2.2.6). The composed translation, extracted via Theorem 2.3.1, is:*

$$T_{12} = 2(M_1 D_2 + D_1 M_2)(M_1 M_2)^\dagger = M_1\,T_2\,M_1^\dagger + T_1 \tag{2.28}$$

*In words: rotate the second translation by the first rotation (via conjugation), then add the first translation. This is the natural composition of rigid motions.*

*Proof.* The product follows from dual number arithmetic (2.23). For the translation: $(M_1 D_2 + D_1 M_2)(M_1 M_2)^\dagger = M_1 D_2 M_2^\dagger M_1^\dagger + D_1 M_2 M_2^\dagger M_1^\dagger = M_1 D_2 M_2^\dagger M_1^\dagger + D_1 M_1^\dagger$. Using $T_k = 2D_k M_k^\dagger$ (Theorem 2.3.1): $D_k M_k^\dagger = \frac{1}{2}T_k$. Substituting and multiplying by 2 gives $T_{12} = M_1 T_2 M_1^\dagger + T_1$. ∎

**Theorem 2.3.3** (Dual $SU(2)$ inverse)**.** *The inverse of a unit dual $SU(2)$ matrix is:*

$$\hat{M}^{-1} = M^\dagger - \varepsilon\,D^\dagger \tag{2.29}$$

*Proof.* Verify $\hat{M}\,\hat{M}^{-1} = \hat{I}$:

$$(M + \varepsilon D)(M^\dagger - \varepsilon D^\dagger) = MM^\dagger + \varepsilon(-M D^\dagger + D M^\dagger) = I + \varepsilon(-MD^\dagger + DM^\dagger)$$

The dual part is $-MD^\dagger + DM^\dagger$. Substituting $D = \frac{1}{2}TM$: $DM^\dagger = \frac{1}{2}T$ and $MD^\dagger = \frac{1}{2}MM^\dagger T^\dagger = \frac{1}{2}T$ (using $T^\dagger = T$). Therefore $-MD^\dagger + DM^\dagger = 0$. ∎

**Theorem 2.3.4** (Equivalence to $4 \times 4$ homogeneous matrix)**.** *The dual $SU(2)$ composition rule (2.28) produces the same rotation and translation as the $4 \times 4$ homogeneous matrix product. The composed rotation $R_1 R_2$ corresponds to $M_1 M_2$ (Theorem 2.2.6), and the composed translation $R_1\mathbf{t}_2 + \mathbf{t}_1$ corresponds to $M_1 T_2 M_1^\dagger + T_1$ — the same vector in the two different encodings. The connection between $SU(2)$ conjugation and $3 \times 3$ matrix-vector multiplication is established in §2.7.*

*Proof.* Standard block multiplication of two $4 \times 4$ matrices gives rotation $R_1 R_2$ and translation $R_1\mathbf{t}_2 + \mathbf{t}_1$. The equivalence $M_1 T_2 M_1^\dagger \leftrightarrow R(M_1)\,\mathbf{t}_2$ is the adjoint representation (§2.7, Theorem 2.7.1). ∎

**Example 2.3.1** (Composing rotation + translation)**.** *Compose two transforms: (1) rotate $90°$ about $\hat{\mathbf{e}}_3$, then (2) translate by $(100, 0, 0)$ km.*

*Transform 1 (rotation only): $\hat{M}_1 = M_3(\pi/2) + \varepsilon \cdot 0$, with $M_3(\pi/2) = \left(\begin{smallmatrix} \frac{1+i}{\sqrt{2}} & 0 \\ 0 & \frac{1-i}{\sqrt{2}} \end{smallmatrix}\right)$ (Example 2.2.2).*

*Transform 2 (translation only): $\hat{M}_2 = I + \varepsilon\,\frac{100}{2}\sigma_1 = I + \varepsilon\left(\begin{smallmatrix} 0 & 50 \\ 50 & 0 \end{smallmatrix}\right)$ (Theorem 2.3.7 with $\lambda = 100$, $k = 1$).*

*Composition (Theorem 2.3.2): $\hat{M}_2 \cdot \hat{M}_1 = I \cdot M_3 + \varepsilon(I \cdot 0 + 50\sigma_1 \cdot M_3) = M_3 + \varepsilon\,50\sigma_1 M_3$.*

*Extracting the translation (Theorem 2.3.1): $T = 2 \cdot 50\sigma_1 M_3 \cdot M_3^\dagger = 100\sigma_1$, encoding $\mathbf{t} = (100, 0, 0)$. ✓*

*Action on $\mathbf{r} = (7000, 0, 0)$: first rotate to $(0, -7000, 0)$ (Example 2.2.4 scaled), then translate to $(100, -7000, 0)$. The 4×4 homogeneous product gives the same result: $R_3(90°)(7000,0,0)^T + (100,0,0)^T = (100, -7000, 0)^T$. ✓*

### Linearity of the Dual Action

**Theorem 2.3.5** (The dual action is affine in the operand)**.** *For fixed $\hat{M}$, the action on a position $R_{\mathrm{pos}}$ is:*

$$R_{\mathrm{pos}}' = M\,R_{\mathrm{pos}}\,M^\dagger + T \tag{2.30}$$

*The conjugation $R_{\mathrm{pos}} \mapsto MR_{\mathrm{pos}}M^\dagger$ is linear in $R_{\mathrm{pos}}$. The addition of $T$ makes the full action affine. For the error state: the constant $T$ contributes zero error (assuming $T$ itself is exact):*

$$\delta(R_{\mathrm{pos}}') = M\,\delta(R_{\mathrm{pos}})\,M^\dagger + \delta(T) \tag{2.31}$$

*When $\delta(T) = 0$, the error transforms by the same conjugation as the physical state.*

*Proof.* Linearity of conjugation: $M(\alpha R_1 + \beta R_2)M^\dagger = \alpha MR_1M^\dagger + \beta MR_2M^\dagger$. The error of the sum $MR_{\mathrm{pos}}M^\dagger + T$ is $M\,\delta(R_{\mathrm{pos}})\,M^\dagger + \delta(T)$ by Ch 1, Theorem 1.3.1 (addition: errors add). ∎

**Remark** (Affine action in the standard pipeline)**.** All coordinate transforms in the standard orbit propagation pipeline have $T = 0$: perifocal → TEME is a pure rotation, and TEME → PEF is a pure rotation with transport. Since $T = 0$ implies $\delta(T) = 0$, the affine term vanishes and the error transforms by the same conjugation as the physical state — no distinction between the state row and the error row. The affine nature of the dual action is therefore a formal property with no operational consequence for the standard pipeline. It becomes relevant only in Generalization sections where translation is nonzero (e.g., maneuver modeling, station-to-satellite offsets).

### Position Extraction Error

**Theorem 2.3.6** (Error in position extraction)**.** *Extracting position from the dual part via $T = 2DM^\dagger$ (Theorem 2.3.1) is a bilinear operation in $D$ and $M^\dagger$. The error bound is:*

$$\|\delta(\mathbf{r})\| \leq 2\|\delta(D)\| + \|\mathbf{r}\|\,\|\delta(M)\| \tag{2.32}$$

*The first term is the dual-part error (scaled by the factor of 2 in the extraction). The second term couples the rotation error through the position magnitude: $\|D\| = \frac{1}{2}\|T\| = \frac{1}{2}\|\mathbf{r}\|$ from (2.25), and $\|M^\dagger\| = 1$ (unitarity).* [P.2.5]

*Proof.* By Ch 1, Theorem 1.3.2 applied to the product $T = 2DM^\dagger$: $\|\delta(T)\| \leq 2(\|D\|\,\|\delta(M^\dagger)\| + \|M^\dagger\|\,\|\delta(D)\|) = 2(\frac{1}{2}\|\mathbf{r}\|\,\|\delta(M)\| + \|\delta(D)\|)$, giving (2.32). ∎

**Remark.** For a LEO satellite at $\|\mathbf{r}\| \approx 7000$ km with TLE-level angle errors $\|\delta(M)\| \sim 5 \times 10^{-9}$ (Theorem 2.2.11 with $\delta(\alpha) \sim 10^{-8}$): the rotation-coupled position error is $\sim 7000 \times 5 \times 10^{-9} \approx 3.5 \times 10^{-5}$ km $\approx 0.035$ m — negligible compared to TLE position accuracy ($\sim 1$ km).

### Elementary Operations

The following theorems verify that the dual $SU(2)$ framework reproduces every elementary operation of the classical $4 \times 4$ homogeneous matrix. Each is confirmed by explicit computation.

**Theorem 2.3.7** (Elementary translation along $\hat{\mathbf{e}}_k$)**.** *A translation by $\lambda$ along coordinate axis $\hat{\mathbf{e}}_k$ (no rotation) is the dual $SU(2)$ element:*

$$\hat{M}_{\mathrm{trans}}(\lambda, k) = I + \varepsilon\,\frac{\lambda}{2}\,\sigma_k \tag{2.33}$$

*The rotation part is $M = I$. The dual part is $D = \frac{\lambda}{2}\sigma_k$, encoding the translation $\mathbf{t} = \lambda\hat{\mathbf{e}}_k$.*

*Action: $R_{\mathrm{pos}}' = R_{\mathrm{pos}} + \lambda\sigma_k$, which adds $\lambda$ to the $k$-th component of $\mathbf{r}$. Explicitly:*

*Along $\hat{\mathbf{e}}_1$:* $D = \frac{\lambda}{2}\left(\begin{smallmatrix} 0 & 1 \\ 1 & 0 \end{smallmatrix}\right)$, giving $r_1' = r_1 + \lambda$, $r_2' = r_2$, $r_3' = r_3$.

*Along $\hat{\mathbf{e}}_2$:* $D = \frac{\lambda}{2}\left(\begin{smallmatrix} 0 & -i \\ i & 0 \end{smallmatrix}\right)$, giving $r_2' = r_2 + \lambda$, with $r_1, r_3$ unchanged.

*Along $\hat{\mathbf{e}}_3$:* $D = \frac{\lambda}{2}\left(\begin{smallmatrix} 1 & 0 \\ 0 & -1 \end{smallmatrix}\right)$, giving $r_3' = r_3 + \lambda$, with $r_1, r_2$ unchanged.

*Proof.* From Definition 2.3.4: $D = \frac{1}{2}TM = \frac{1}{2}(\lambda\sigma_k)(I) = \frac{\lambda}{2}\sigma_k$. Extraction (Theorem 2.3.1): $T = 2DM^\dagger = \lambda\sigma_k$. The encoded vector is $\lambda\hat{\mathbf{e}}_k$ (Definition 2.2.5). The action: $R_{\mathrm{pos}}' = MR_{\mathrm{pos}}M^\dagger + T = R_{\mathrm{pos}} + \lambda\sigma_k$. ∎

**Theorem 2.3.8** (General translation)**.** *A translation by $\mathbf{t} = (\lambda_1, \lambda_2, \lambda_3)$ is:*

$$\hat{M}_{\mathrm{trans}}(\mathbf{t}) = I + \varepsilon\,\frac{1}{2}(\lambda_1\sigma_1 + \lambda_2\sigma_2 + \lambda_3\sigma_3) \tag{2.34}$$

*This is the superposition of the three elementary translations (Theorem 2.3.7).*

*Proof.* $T = \lambda_1\sigma_1 + \lambda_2\sigma_2 + \lambda_3\sigma_3$ encodes $\mathbf{t}$ (Definition 2.2.5). $D = \frac{1}{2}TI = \frac{1}{2}T$. The action adds $T$ to $R_{\mathrm{pos}}$. ∎

The three elementary rotations are verified by direct conjugation. The first axis is computed in full; the other two follow the same pattern with different matrix structures.

**Remark** (Sign convention)**.** Throughout §§2.2–2.5, conjugation $V' = M_k(\alpha)\,V\,M_k(\alpha)^\dagger$ implements a **passive** rotation: it expresses $\mathbf{v}$ in a frame rotated by $+\alpha$ about $\hat{\mathbf{e}}_k$. To rotate a vector **actively** by $+\alpha$, negate the angle. Example 2.2.4 demonstrated this: conjugating $(1,0,0)$ by $M_3(+\pi/2)$ gave $(0,-1,0)$ (passive), while $M_3(-\pi/2)$ gives $(0,+1,0)$ (active). The negative angles in §2.5 ($M_3(-\Omega)$, $M_1(-i)$, $M_3(-u)$) implement passive rotations that undo the orbital angles.

**Theorem 2.3.9** (Elementary rotation about $\hat{\mathbf{e}}_1$)**.** *The passive rotation $V' = M_1(\alpha)\,V\,M_1(\alpha)^\dagger$ transforms $\mathbf{v} = (v_1, v_2, v_3)$ to:*

$$v_1' = v_1, \qquad v_2' = v_2\cos\alpha + v_3\sin\alpha, \qquad v_3' = -v_2\sin\alpha + v_3\cos\alpha \tag{2.35}$$

*This is the standard passive rotation matrix $R_1(\alpha)$.*

*Proof.* Write $c = \cos\frac{\alpha}{2}$, $s = \sin\frac{\alpha}{2}$. Computing $V' = M_1 V M_1^\dagger$ entry by entry:

$(V')_{11} = c^2 v_3 + ics(v_1 + iv_2) - ics(v_1 - iv_2) - s^2 v_3 = (c^2 - s^2)v_3 - 2csv_2 = v_3\cos\alpha - v_2\sin\alpha$

$(V')_{12} = (c^2 + s^2)v_1 - i(c^2 - s^2)v_2 - 2icsv_3 = v_1 - i(v_2\cos\alpha + v_3\sin\alpha)$

By Definition 2.2.4, $(V')_{11} = v_3'$ and $(V')_{12} = v_1' - iv_2'$. Therefore $v_3' = v_3\cos\alpha - v_2\sin\alpha$, $v_1' = v_1$, $v_2' = v_2\cos\alpha + v_3\sin\alpha$, confirming (2.35). ∎

**Theorem 2.3.10** (Elementary rotation about $\hat{\mathbf{e}}_3$)**.** *The passive rotation $V' = M_3(\alpha)\,V\,M_3(\alpha)^\dagger$ transforms $\mathbf{v}$ to:*

$$v_1' = v_1\cos\alpha + v_2\sin\alpha, \qquad v_2' = -v_1\sin\alpha + v_2\cos\alpha, \qquad v_3' = v_3 \tag{2.36}$$

*This is the standard passive rotation matrix $R_3(\alpha)$.*

*Proof.* Since $M_3(\alpha) = \left(\begin{smallmatrix} e^{i\alpha/2} & 0 \\ 0 & e^{-i\alpha/2} \end{smallmatrix}\right)$ is diagonal, conjugation multiplies the off-diagonal entry $(v_1 - iv_2)$ by $e^{i\alpha}$. Therefore $v_1' - iv_2' = (v_1 - iv_2)(\cos\alpha + i\sin\alpha) = (v_1\cos\alpha + v_2\sin\alpha) - i(-v_1\sin\alpha + v_2\cos\alpha)$, and $v_3' = v_3$. ∎

The $\hat{\mathbf{e}}_2$ rotation follows the same computational pattern as $\hat{\mathbf{e}}_1$, with the roles of the components permuted.

**Theorem 2.3.11** (Elementary rotation about $\hat{\mathbf{e}}_2$)**.** *The passive rotation $V' = M_2(\alpha)\,V\,M_2(\alpha)^\dagger$ transforms $\mathbf{v}$ to:*

$$v_1' = v_1\cos\alpha - v_3\sin\alpha, \qquad v_2' = v_2, \qquad v_3' = v_1\sin\alpha + v_3\cos\alpha \tag{2.37}$$

*This is the standard passive rotation matrix $R_2(\alpha)$.*

*Proof.* Write $c = \cos\frac{\alpha}{2}$, $s = \sin\frac{\alpha}{2}$. Computing $V' = M_2 V M_2^\dagger$:

$(V')_{11} = (c^2 - s^2)v_3 + 2csv_1 = v_3\cos\alpha + v_1\sin\alpha$

$(V')_{12} = (c^2 - s^2)v_1 - i(c^2 + s^2)v_2 - 2csv_3 = (v_1\cos\alpha - v_3\sin\alpha) - iv_2$

Therefore $v_3' = v_1\sin\alpha + v_3\cos\alpha$, $v_1' = v_1\cos\alpha - v_3\sin\alpha$, $v_2' = v_2$. ∎

The remaining theorems verify that composition of dual $SU(2)$ elements reproduces the expected behavior for the three fundamental cases: rotation + translation, two translations, and two rotations about different axes.

**Theorem 2.3.12** (Composition of translation and rotation)**.** *A rotation by $\alpha$ about $\hat{\mathbf{e}}_3$ followed by a translation by $\lambda$ along $\hat{\mathbf{e}}_1$:*

$$\hat{M} = \hat{M}_{\mathrm{trans}}(\lambda, 1) \cdot \hat{M}_{\mathrm{rot}}(\alpha, 3) = M_3(\alpha) + \varepsilon\frac{\lambda}{2}\sigma_1 M_3(\alpha) \tag{2.38}$$

*Rotation part: $M_3(\alpha)$. Translation: $T = \lambda\sigma_1$, giving $\mathbf{t} = \lambda\hat{\mathbf{e}}_1$. The action: $\mathbf{r}' = R_3(\alpha)\mathbf{r} + \lambda\hat{\mathbf{e}}_1$.*

*Proof.* Dual product (Theorem 2.3.2): non-dual part $= I \cdot M_3(\alpha) = M_3(\alpha)$; dual part $= I \cdot 0 + \frac{\lambda}{2}\sigma_1 \cdot M_3(\alpha)$. Extraction: $T = 2 \cdot \frac{\lambda}{2}\sigma_1 M_3 \cdot M_3^\dagger = \lambda\sigma_1$. ∎

**Theorem 2.3.13** (Composition of two translations)**.** *Two translations compose by vector addition:*

$$\hat{M}_{\mathrm{trans}}(\mathbf{t}_1) \cdot \hat{M}_{\mathrm{trans}}(\mathbf{t}_2) = \hat{M}_{\mathrm{trans}}(\mathbf{t}_1 + \mathbf{t}_2) \tag{2.39}$$

*Proof.* Both have $M = I$. By Theorem 2.3.2: dual part $= D_1 + D_2 = \frac{1}{2}(T_1 + T_2)$. The combined translation encodes $\mathbf{t}_1 + \mathbf{t}_2$. ∎

**Theorem 2.3.14** (Composition of two rotations about different axes)**.** *A rotation by $\alpha$ about $\hat{\mathbf{e}}_1$ followed by a rotation by $\beta$ about $\hat{\mathbf{e}}_3$ (both passive):*

$$M = M_3(\beta) \cdot M_1(\alpha) = \begin{pmatrix} e^{i\beta/2}\cos\frac{\alpha}{2} & ie^{i\beta/2}\sin\frac{\alpha}{2} \\ ie^{-i\beta/2}\sin\frac{\alpha}{2} & e^{-i\beta/2}\cos\frac{\alpha}{2} \end{pmatrix} \tag{2.40}$$

*This is a single $SU(2)$ element encoding a combined rotation that is not about any coordinate axis.*

*Proof.* Direct $2 \times 2$ matrix multiplication. Verification: $\det = e^{i\beta/2}e^{-i\beta/2}(\cos^2\frac{\alpha}{2} + \sin^2\frac{\alpha}{2}) = 1$. ∎

---

## §2.4 The Complete Rigid Body State

The dual $SU(2)$ element $\hat{M}$ encodes a rigid body's configuration — where it is and how it is oriented. The velocity — both translational and angular — is encoded as a second dual element: a pure dual quaternion belonging to the Lie algebra $\mathfrak{se}(3)$. The complete state is a pair of dual quaternions $(\hat{M}, \hat{\Omega}_b)$, and the connection between them is the kinematic equation $\dot{\hat{M}} = \hat{M} \cdot \frac{1}{2}\hat{\Omega}_b$.

### State Definition

**Definition 2.4.1** (Velocity dual quaternion)**.** *The velocity of a rigid body is encoded as a pure dual quaternion:*

$$\hat{\Omega}_b = \Omega_b + \varepsilon\,V_b \tag{2.41}$$

*where $\Omega_b$ is the body-frame angular velocity (traceless Hermitian, 3 real parameters) and $V_b$ is the body-frame translational velocity (traceless Hermitian, 3 real parameters). Both parts are elements of the Lie algebra $\mathfrak{su}(2)$. Total: 6 real parameters, 6 DOF.*

*The configuration dual quaternion $\hat{M} = M + \varepsilon D$ (Definition 2.3.2) is a unit dual quaternion — an element of the Lie group $SE(3)$. The velocity dual quaternion $\hat{\Omega}_b$ is a pure dual quaternion — an element of the Lie algebra $\mathfrak{se}(3)$. The two objects have different algebraic types and different transformation laws (Theorem 2.4.2).*

**Definition 2.4.2** (Complete rigid body state)**.** *The state is the pair $(\hat{M}, \hat{\Omega}_b)$:*

| Dual quaternion | Type | Encodes | Parameters | DOF |
|-----------------|------|---------|------------|-----|
| $\hat{M} = M + \varepsilon D$ | Unit (Lie group $SE(3)$) | Orientation + position | 7 (4+3) | 6 (3+3) |
| $\hat{\Omega}_b = \Omega_b + \varepsilon V_b$ | Pure (Lie algebra $\mathfrak{se}(3)$) | Angular velocity + translational velocity | 6 (3+3) | 6 (3+3) |

*Total: 13 real parameters, 12 degrees of freedom.*

**Assumption 2.4.1** (Instantaneous angular velocity)**.** *The angular velocity $\Omega_T$ of a rotating reference frame is treated as constant during a single transformation step. For time-varying $\Omega_T(t)$, additional angular acceleration terms appear. For Earth rotation, this assumption is valid: the rotation rate varies by $< 10^{-8}$ rad/s over any practical propagation step.*

### The Kinematic Equation

**Theorem 2.4.1** (Kinematic equation)**.** *The time evolution of the configuration is:*

$$\dot{\hat{M}} = \hat{M} \cdot \tfrac{1}{2}\hat{\Omega}_b \tag{2.42}$$

*Proof.* Expand: $\hat{M} \cdot \frac{1}{2}\hat{\Omega}_b = (M + \varepsilon D) \cdot \frac{1}{2}(\Omega_b + \varepsilon V_b) = \frac{1}{2}M\Omega_b + \varepsilon\frac{1}{2}(MV_b + D\Omega_b)$.

The non-dual part: $\frac{1}{2}M\Omega_b = \dot{M}$. This is the standard $SU(2)$ kinematic equation for body-frame angular velocity.

The dual part: $\frac{1}{2}(MV_b + D\Omega_b)$. Substituting $D = \frac{1}{2}R_{\mathrm{pos}}M$ and $V_b = M^\dagger V_{\mathrm{vel}} M$ (the body-frame velocity):

$$\tfrac{1}{2}(M \cdot M^\dagger V_{\mathrm{vel}} M + \tfrac{1}{2}R_{\mathrm{pos}} M \Omega_b) = \tfrac{1}{2}V_{\mathrm{vel}} M + \tfrac{1}{4}R_{\mathrm{pos}} M\Omega_b$$

This equals $\dot{D} = \frac{d}{dt}(\frac{1}{2}R_{\mathrm{pos}} M) = \frac{1}{2}\dot{R}_{\mathrm{pos}} M + \frac{1}{2}R_{\mathrm{pos}} \dot{M} = \frac{1}{2}V_{\mathrm{vel}} M + \frac{1}{2}R_{\mathrm{pos}} \cdot \frac{1}{2}M\Omega_b$. ✓ ∎

### Spatial and Body Twists

**Definition 2.4.3** (Spatial twist)**.** *The spatial twist is the velocity expressed in the spatial (inertial) frame, obtained by the adjoint action:*

$$\hat{\Omega}_s = \hat{M}\,\hat{\Omega}_b\,\hat{M}^{-1} \tag{2.43}$$

*The non-dual part $\Omega_s = M\Omega_b M^\dagger$ is the angular velocity in the spatial frame. The dual part encodes the spatial twist velocity (see Theorem 2.4.3).*

**Theorem 2.4.2** (Frame transform of the velocity)**.** *Under a frame transform with configuration $\hat{M}_T$ and frame velocity $\hat{\Omega}_{T,s}$:*

*The configuration transforms by left multiplication (Theorem 2.3.2):*

$$\hat{M}' = \hat{M}_T \cdot \hat{M} \tag{2.44}$$

*The spatial twist transforms by the adjoint action plus the frame's own twist:*

$$\hat{\Omega}_s' = \hat{M}_T\,\hat{\Omega}_s\,\hat{M}_T^{-1} + \hat{\Omega}_{T,s} \tag{2.45}$$

*This is the complete velocity transform. For fixed $\hat{M}_T$, the adjoint action is linear in $\hat{\Omega}_s$ (Principle 1, §2.6).*

*Proof.* Equation (2.44) is Theorem 2.3.2. For (2.45): the adjoint action of a Lie group on its Lie algebra transforms algebra elements between frames. The conjugation $\hat{M}_T\,\hat{\Omega}_s\,\hat{M}_T^{-1}$ rotates and couples the twist components; the addition of $\hat{\Omega}_{T,s}$ accounts for the frame's own motion. ∎

### Spatial Twist Velocity and Material Velocity

The spatial twist velocity $V_s$ (the dual part of $\hat{\Omega}_s$) differs from the material velocity $V_{\mathrm{vel}}$ of the body's reference point. The difference is the transport term.

**Theorem 2.4.3** (Spatial twist expansion)**.** *The dual part of $\hat{\Omega}_s = \hat{M}\,\hat{\Omega}_b\,\hat{M}^{-1}$ is:*

$$V_s = V_{\mathrm{vel}} + \frac{1}{2i}[\Omega_s, R_{\mathrm{pos}}] \tag{2.46}$$

*where $\Omega_s = M\Omega_b M^\dagger$ is the spatial angular velocity and $R_{\mathrm{pos}} = 2DM^\dagger$ is the position. The second term encodes $\boldsymbol{\omega}_s \times \mathbf{r}$ — the velocity contribution from the body's rotation at its reference point.*

*Proof.* Expand $\hat{M}\,\hat{\Omega}_b\,\hat{M}^{-1}$ using $\hat{M} = M + \varepsilon D$ and $\hat{M}^{-1} = M^\dagger - \varepsilon M^\dagger D M^\dagger$ (the quaternion conjugate, not the dual conjugate):

$$\hat{M}\,\hat{\Omega}_b\,\hat{M}^{-1} = (M + \varepsilon D)(\Omega_b + \varepsilon V_b)(M^\dagger - \varepsilon M^\dagger D M^\dagger)$$

Non-dual part: $M\Omega_b M^\dagger = \Omega_s$. ✓

Dual part: $MV_b M^\dagger + D\Omega_b M^\dagger - M\Omega_b M^\dagger D M^\dagger$. Substituting $D = \frac{1}{2}R_{\mathrm{pos}}M$:

$D\Omega_b M^\dagger = \frac{1}{2}R_{\mathrm{pos}} M\Omega_b M^\dagger = \frac{1}{2}R_{\mathrm{pos}}\Omega_s$

$M\Omega_b M^\dagger D M^\dagger = \Omega_s \cdot \frac{1}{2}R_{\mathrm{pos}} M \cdot M^\dagger = \frac{1}{2}\Omega_s R_{\mathrm{pos}}$

Dual part $= MV_b M^\dagger + \frac{1}{2}(R_{\mathrm{pos}}\Omega_s - \Omega_s R_{\mathrm{pos}}) = V_{\mathrm{vel}} + \frac{1}{2}[R_{\mathrm{pos}}, \Omega_s]$

Since $[R_{\mathrm{pos}}, \Omega_s] = -[\Omega_s, R_{\mathrm{pos}}]$, this is $V_{\mathrm{vel}} - \frac{1}{2}[\Omega_s, R_{\mathrm{pos}}] = V_{\mathrm{vel}} + \frac{1}{2i} \cdot i \cdot (-[\Omega_s, R_{\mathrm{pos}}])$. Using $\frac{1}{2i}[A,B] \leftrightarrow \mathbf{a} \times \mathbf{b}$: the dual part encodes $\mathbf{v}_{\mathrm{vel}} + \boldsymbol{\omega}_s \times \mathbf{r}$. ∎

**Corollary 2.4.1** (Material velocity from spatial twist)**.** *The material velocity is recovered from the spatial twist by:*

$$V_{\mathrm{vel}} = V_s - \frac{1}{2i}[\Omega_s, R_{\mathrm{pos}}] \tag{2.47}$$

*For the point-mass case ($\Omega_b = 0$, hence $\Omega_s = 0$): $V_s = V_{\mathrm{vel}}$ — the spatial twist velocity equals the material velocity.*

### The Transport Term

**Corollary 2.4.2** (Transport term from frame transform)**.** *Under a frame transform (Theorem 2.4.2) with $T_T = 0$ (shared origin) and frame angular velocity $\Omega_T$, the material velocity transforms as:*

$$V_{\mathrm{vel}}' = M_T\,V_{\mathrm{vel}}\,M_T^\dagger - \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}'] \tag{2.48}$$

*The transport term $-\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$ encodes $-\boldsymbol{\omega}_T \times \mathbf{r}'$: transforming to a rotating frame subtracts the frame's rotational velocity at the satellite's position.*

*Verification (geostationary satellite):* $\mathbf{v}_{\mathrm{TEME}} \approx (0, 3.075, 0)$ km/s, $\mathbf{r}_{\mathrm{PEF}} = (42.68, 0, 0)$ km, $\boldsymbol{\omega}_{E} \times \mathbf{r}_{\mathrm{PEF}} = (0, 3.075, 0)$ km/s. With minus: $\mathbf{v}_{\mathrm{PEF}} = (0,0,0)$ ✓. With plus: $(0, 6.15, 0)$ — wrong.

*Proof.* The spatial twist transforms cleanly by (2.45). Converting back to material velocity via Corollary 2.4.1, the difference between the old and new $\frac{1}{2i}[\Omega_s, R_{\mathrm{pos}}]$ terms produces the transport coupling. For $T_T = 0$ and point-mass ($\Omega_b = 0$), the spatial twist velocity equals the material velocity in both frames, and the only additional term is $-\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$ from the frame's twist $\hat{\Omega}_{T,s}$. ∎

### The Transport Theorem

### Calculus Verification of the Transport Term

The transport term in Corollary 2.4.2 was derived from the Lie algebra structure — the adjoint action on twists and the spatial-to-material conversion. As verification, the same result follows from direct calculus: differentiating the $SU(2)$ conjugation $M(t)\,V\,M(t)^\dagger$ with respect to time.

**Theorem 2.4.4** (Transport from time derivative)**.** *If $M_T(t)$ is a time-varying rotation with angular velocity $\boldsymbol{\omega}_T$, the velocity of a position vector that is stationary in the rotating frame is:*

$$\frac{d}{dt}\left(M_T R_{\mathrm{pos}} M_T^\dagger\right) \;\longleftrightarrow\; \boldsymbol{\omega}_T \times \mathbf{r}' \tag{2.49}$$

*Proof.* $\dot{M}_T = \frac{i}{2}\Omega_T' M_T$ and $\dot{M}_T^\dagger = -M_T^\dagger \frac{i}{2}\Omega_T'$. By the product rule:

$$\frac{d}{dt}(M_T R_{\mathrm{pos}} M_T^\dagger) = \frac{i}{2}[\Omega_T', R_{\mathrm{pos}}'] = \frac{1}{2i}[R_{\mathrm{pos}}', \Omega_T'] \;\longleftrightarrow\; \mathbf{r}' \times \boldsymbol{\omega}_T = -\boldsymbol{\omega}_T \times \mathbf{r}'$$

The negative sign: this is the inertial-frame velocity of a point co-rotating with the frame. The material velocity transform (2.48) subtracts this to get the velocity in the rotating frame. The geostationary verification (Example 2.4.1) confirms the sign. ∎

### Properties of the Transport Term

**Theorem 2.4.5** (Transport term is linear in position)**.** *For fixed $\Omega_T$, the transport map $R_{\mathrm{pos}} \mapsto -\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ (equation 2.48) is linear:*

$$-\frac{1}{2i}[\Omega_T, \alpha R_1 + \beta R_2] = \alpha\,\bigl(-\frac{1}{2i}[\Omega_T, R_1]\bigr) + \beta\,\bigl(-\frac{1}{2i}[\Omega_T, R_2]\bigr) \tag{2.50}$$

*Proof.* The commutator is bilinear. Negation preserves linearity. ∎

**Corollary 2.4.3** (Transport term error)**.** *The transport velocity error from position error is:*

$$\delta(V_{\mathrm{transport}}) = -\frac{1}{2i}[\Omega_T, \delta(R_{\mathrm{pos}})] \tag{2.51}$$

*The error magnitude:*

$$\|\delta(\mathbf{v}_{\mathrm{transport}})\| \leq |\omega_T| \cdot \|\delta(\mathbf{r})\| \tag{2.52}$$

*Proof.* Linearity (Theorem 2.4.5) applied to $\delta(R_{\mathrm{pos}})$ gives (2.51). For the magnitude: $\|\boldsymbol{\omega} \times \delta\mathbf{r}\| \leq |\omega|\,\|\delta\mathbf{r}\|$. ∎

**Example 2.4.1** (Transport velocity from Earth rotation)**.** *A satellite at $\mathbf{r}_{\mathrm{PEF}} = (6700, 1200, 0)$ km. Earth rotation: $\omega_E = 7.292 \times 10^{-5}$ rad/s about $\hat{\mathbf{e}}_3$.*

*Cartesian: $\boldsymbol{\omega}_{E} \times \mathbf{r}_{\mathrm{PEF}} = (-\omega_E \cdot 1200,\; \omega_E \cdot 6700,\; 0) = (-87.5, 488.6, 0)$ m/s. Magnitude: $496$ m/s.*

*$SU(2)$: $\frac{1}{2i}[\omega_E\sigma_3, R_{\mathrm{pos}}]$ — verified in Theorem 2.5.2 (equation 2.59).*

*Geostationary check: $(42.68, 0, 0)$ km, $\mathbf{v}_{\mathrm{TEME}} = (0, 3.075, 0)$ km/s. Transport: $(0, 3.075, 0)$ km/s. Material velocity in PEF: $\mathbf{v}_{\mathrm{PEF}} = \mathbf{v}_{\mathrm{TEME}} - \boldsymbol{\omega}_{E} \times \mathbf{r}_{\mathrm{PEF}} = \mathbf{0}$. ✓*

*Position error $\delta r = 1$ km contributes transport velocity error $|\omega_E| \times 1 \approx 0.073$ m/s (Corollary 2.4.3).*

### Point-Mass Reduction

**Corollary 2.4.4** (Point-mass reduction)**.** *For a point-mass satellite with no body attitude tracking (the standard mode):*
- *$\Omega_b = 0$ (no body angular velocity), so $\hat{\Omega}_b = \varepsilon V_b$*
- *The spatial twist velocity equals the material velocity (Corollary 2.4.1)*
- *The state reduces to 6 active quantities: $R_{\mathrm{pos}}$ (3) and $V_{\mathrm{vel}}$ (3)*
- *The configuration transform (2.44) and material velocity transform (2.48) are the complete equations*

*Proof.* Substitute $\Omega_b = 0$ into (2.41)–(2.48). ∎

### Near-Identity Perturbations

**Theorem 2.4.6** (Near-identity composition is additive)**.** *A small perturbation to the rotation is $M_{\mathrm{pert}} = I + \frac{i}{2}\varepsilon_k\sigma_k + O(\varepsilon^2)$ where $|\varepsilon_k| \ll 1$. Composing two near-identity perturbations:*

$$(I + \frac{i}{2}A)(I + \frac{i}{2}B) = I + \frac{i}{2}(A + B) + O(\varepsilon^2) \tag{2.53}$$

*To first order, perturbations compose by addition.*

*Proof.* Expand: $I + \frac{i}{2}A + \frac{i}{2}B + \frac{i^2}{4}AB = I + \frac{i}{2}(A+B) - \frac{1}{4}AB$. The last term is $O(\varepsilon^2)$. ∎

**Remark.** Theorem 2.4.6 applies to the short-period and long-period corrections of Ch 18–19. Each correction is a near-identity perturbation, and composing them is approximately additive. The second-order cross term is $O(\varepsilon_1\varepsilon_2) \sim O(J_2^2) \sim 10^{-6}$, well below the model accuracy $\delta_a$.

---

## §2.5 The Orbit Propagation Coordinate Transform Pipeline

The orbit propagation pipeline transforms the satellite state $(\hat{M}, \hat{\Omega}_b)$ from the orbital-plane frame (perifocal) through the inertial frame (TEME) to the Earth-fixed frame (PEF). Each step acts on the pair: the configuration $\hat{M}$ transforms by left multiplication (Theorem 2.4.2, equation 2.44), and the velocity twist $\hat{\Omega}_s$ transforms by the adjoint action (equation 2.45). In the point-mass reduction ($\Omega_b = 0$), the material velocity equals the spatial twist velocity (Corollary 2.4.1), and the transforms reduce to conjugation plus the transport correction (equation 2.48).

### Frame Definitions

**Definition 2.5.1** (Perifocal frame)**.** *The orbital-plane frame with $\hat{\mathbf{x}}$ toward perigee, $\hat{\mathbf{z}}$ along the angular momentum vector, and $\hat{\mathbf{y}}$ completing the right-handed triad. Position in this frame is $(r\cos\nu, r\sin\nu, 0)$ where $r$ is the orbital radius and $\nu$ is the true anomaly.*

**Definition 2.5.2** (TEME frame)**.** *True Equator, Mean Equinox — the reference frame of the propagator (Ch 30). The $\hat{\mathbf{z}}$ axis is aligned with the Earth's instantaneous rotation axis, and $\hat{\mathbf{x}}$ points toward the mean vernal equinox.* [A.2.1]

**Definition 2.5.3** (PEF frame)**.** *Pseudo Earth-Fixed — TEME rotated by the Greenwich Mean Sidereal Time angle $\theta_{\mathrm{GMST}}$ (Ch 29) about the polar axis. Approximately Earth-fixed (neglecting polar motion).*

### Transform 1: Perifocal → TEME

The first step transforms from the orbital plane to the inertial frame. It is a pure rotation with no translation and no frame angular velocity. The frame transform dual quaternion is $\hat{M}_T = M_{\mathrm{PF \to TEME}}$ (pure rotation, $D_T = 0$), with $\hat{\Omega}_{T,s} = 0$ (inertial frame, no frame velocity).

**Theorem 2.5.1** (Perifocal → TEME rotation)**.** *The rotation from the perifocal frame to TEME is the 3-1-3 Euler composition:*

$$M_{\mathrm{PF \to TEME}} = M_3(-\Omega_{\mathrm{node}}) \cdot M_1(-i) \cdot M_3(-u) \tag{2.54}$$

*where $\Omega_{\mathrm{node}}$ is the right ascension of the ascending node, $i$ is the inclination, and $u = \omega + \nu$ is the argument of latitude. The negative signs implement passive rotation (frame rotation). Since $\hat{\Omega}_{T,s} = 0$ (inertial target frame), the adjoint action (2.45) reduces to pure conjugation, and the material velocity equals the spatial twist velocity (Corollary 2.4.1):*

$$R_{\mathrm{pos,TEME}} = M\,R_{\mathrm{pos,PF}}\,M^\dagger, \qquad V_{\mathrm{vel,TEME}} = M\,V_{\mathrm{vel,PF}}\,M^\dagger \tag{2.55}$$

*In the dual quaternion form: $\hat{M}' = M_{\mathrm{PF \to TEME}} \cdot \hat{M}$ (configuration), and $\hat{\Omega}_s' = M_{\mathrm{PF \to TEME}}\,\hat{\Omega}_s\,M_{\mathrm{PF \to TEME}}^{-1}$ (velocity twist, with no frame twist added).*

*Proof.* The perifocal frame is obtained from TEME by three successive rotations: (1) about $\hat{\mathbf{z}}$ by $\Omega_{\mathrm{node}}$, (2) about $\hat{\mathbf{x}}$ by $i$, (3) about the new $\hat{\mathbf{z}}$ by $u$. The inverse (PF → TEME) reverses these with negated angles. Composition follows from Theorem 2.2.6. With $\hat{\Omega}_{T,s} = 0$, the transport term in (2.48) vanishes. ∎

**Example 2.5.1** (ISS-like orbit: PF → TEME rotation)**.** *Orbital elements: $\Omega_{\mathrm{node}} = 30°$, $i = 51.6°$, $u = 45°$. The three $SU(2)$ factors:*

*$M_3(-u) = M_3(-45°)$: $\alpha/2 = -22.5°$, so $M_3(-45°) = \left(\begin{smallmatrix} e^{-i\pi/8} & 0 \\ 0 & e^{i\pi/8} \end{smallmatrix}\right) = \left(\begin{smallmatrix} 0.9239 - 0.3827i & 0 \\ 0 & 0.9239 + 0.3827i \end{smallmatrix}\right)$.*

*$M_1(-i) = M_1(-51.6°)$: $\alpha/2 = -25.8°$, so $M_1(-51.6°) = \left(\begin{smallmatrix} 0.9003 & -0.4352i \\ -0.4352i & 0.9003 \end{smallmatrix}\right)$.*

*$M_3(-\Omega) = M_3(-30°)$: $\alpha/2 = -15°$, so $M_3(-30°) = \left(\begin{smallmatrix} 0.9659 - 0.2588i & 0 \\ 0 & 0.9659 + 0.2588i \end{smallmatrix}\right)$.*

*Product: $M_{\mathrm{PF\to TEME}} = M_3(-30°) \cdot M_1(-51.6°) \cdot M_3(-45°)$. The two diagonal factors bracket the general factor — this is a sequence of two phase multiplications and one general $SU(2)$ multiply. The result is a single $SU(2)$ element encoding the full 3-1-3 Euler rotation, with no singularity risk (all entries well-conditioned since $i = 51.6° \neq 0$).*

**Remark** (Computational structure)**.** Since $M_3(\alpha)$ is diagonal (equation 2.10), the first and last factors in (2.54) are phase multiplications. Each $M_3$ multiplication costs 4 real multiplications instead of 16 for a general $SU(2)$ product. The two diagonal factors cannot be combined because $M_1(-i)$ lies between them.

**Proposition 2.5.1** (Singularity-free at $i = 0$)**.** *When $i = 0$ (equatorial orbit), $M_1(0) = I$ and:*

$$M_{\mathrm{PF \to TEME}} = M_3(-\Omega) \cdot M_3(-u) = M_3(-\Omega - u) \tag{2.56}$$

*The individual angles $\Omega$ and $u$ are undefined (the ascending node is arbitrary for an equatorial orbit), but their sum $\Omega + u$ is well-defined (it is the satellite's longitude). The product (2.56) is a single well-conditioned diagonal matrix. Similarly, when $e = 0$ (circular orbit), $\omega$ is undefined but $u = \omega + \nu$ is well-defined, and $M_3(-u)$ is non-singular. No conditional guards or Lyddane-type modifications are needed.*

*Proof.* $M_1(0) = I$. Then $M_3(-\Omega) \cdot I \cdot M_3(-u) = M_3(-\Omega-u)$, with entries $e^{\pm i(\Omega+u)/2}$ — well-defined for any value of $\Omega + u$. ∎

### Transform 2: TEME → PEF

The second step accounts for Earth's rotation. The frame transform has $\hat{M}_T = M_3(\theta_{\mathrm{GMST}})$ (pure rotation, $D_T = 0$) and the frame's spatial twist $\hat{\Omega}_{T,s} = \Omega_T + \varepsilon \cdot 0$ (the frame rotates at $\omega_E$ but its origin has zero velocity). The adjoint action (2.45) rotates the satellite's twist, and the addition of $\hat{\Omega}_{T,s}$ contributes the Earth rotation rate. The material velocity transform (2.48) introduces the transport term.

**Theorem 2.5.2** (TEME → PEF with transport)**.** *The rotation is $M_T = M_3(\theta_{\mathrm{GMST}})$ and the angular velocity is $\Omega_T = \omega_E\sigma_3$:*

$$M_T = \begin{pmatrix} e^{i\theta/2} & 0 \\ 0 & e^{-i\theta/2} \end{pmatrix}, \qquad \Omega_T = \begin{pmatrix} \omega_E & 0 \\ 0 & -\omega_E \end{pmatrix} \tag{2.57}$$

*No translation: $T_T = 0$. The material velocity transform (equation 2.48):*

$$V_{\mathrm{vel,PEF}} = M_T\,V_{\mathrm{vel,TEME}}\,M_T^\dagger - \frac{1}{2i}[\Omega_T, R_{\mathrm{pos,PEF}}] \tag{2.58}$$

*The transport commutator has a simple form because both $M_T$ and $\Omega_T$ are diagonal. Writing $R_{\mathrm{pos,PEF}} = \left(\begin{smallmatrix} r_3 & r_- \\ r_+ & -r_3 \end{smallmatrix}\right)$ where $r_\pm = r_1 \pm ir_2$:*

$$\frac{1}{2i}[\Omega_T, R_{\mathrm{pos,PEF}}] = \omega_E\begin{pmatrix} 0 & -ir_- \\ ir_+ & 0 \end{pmatrix} \tag{2.59}$$

*The transport velocity is $\boldsymbol{\omega}_{E} \times \mathbf{r}_{\mathrm{PEF}} = (-\omega_E r_2,\; \omega_E r_1,\; 0)$.* [M.2.1]

*Proof.* The commutator $[\sigma_3, R_{\mathrm{pos}}] = \sigma_3 R_{\mathrm{pos}} - R_{\mathrm{pos}} \sigma_3$:

$$\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}\begin{pmatrix} r_3 & r_- \\ r_+ & -r_3 \end{pmatrix} - \begin{pmatrix} r_3 & r_- \\ r_+ & -r_3 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} = \begin{pmatrix} 0 & 2r_- \\ -2r_+ & 0 \end{pmatrix}$$

Multiplying by $\frac{\omega_E}{2i}$: the $(1,2)$ entry is $\frac{\omega_E}{2i} \cdot 2r_- = -i\omega_E r_- = -i\omega_E(r_1 - ir_2) = -\omega_E r_2 - i\omega_E r_1$. By Definition 2.2.4, this encodes $v_1 - iv_2$, so $v_1 = -\omega_E r_2$ and $v_2 = \omega_E r_1$, confirming $\boldsymbol{\omega}_{E} \times \mathbf{r}_{\mathrm{PEF}}$. ∎

### Pipeline Composition

**Theorem 2.5.3** (Perifocal → PEF composition)**.** *The complete transform is:*

$$M_{\mathrm{PF \to PEF}} = M_3(\theta_{\mathrm{GMST}}) \cdot M_3(-\Omega_{\mathrm{node}}) \cdot M_1(-i) \cdot M_3(-u) \tag{2.60}$$

*The first two diagonal factors compose trivially: $M_3(\theta) \cdot M_3(-\Omega) = M_3(\theta - \Omega)$, reducing the chain to three multiplications. In the dual quaternion form: $\hat{M}_{\mathrm{PF \to PEF}} = \hat{M}_{\mathrm{TEME \to PEF}} \cdot \hat{M}_{\mathrm{PF \to TEME}}$. The material velocity includes the transport term from the TEME → PEF step only:*

$$V_{\mathrm{vel,PEF}} = M_{\mathrm{total}}\,V_{\mathrm{vel,PF}}\,M_{\mathrm{total}}^\dagger - \frac{1}{2i}[\Omega_E, R_{\mathrm{pos,PEF}}] \tag{2.61}$$

*Proof.* Composition of the two rotations gives (2.60) by Theorem 2.2.6. The velocity: Transform 1 has $\Omega_T = 0$ (no transport), so no transport term is added. The only transport contribution comes from Transform 2 (Earth rotation), acting on the already-rotated position. ∎

### Error Sources

The following table traces each error source through the $SU(2)$ construction. The symbolic bounds are precision-independent; the binary64 column evaluates them for IEEE 754 double precision.

| Source | Category | Bound | Binary64 |
|--------|----------|-------|----------|
| $\Omega_{\mathrm{node}}$ from TLE | $\sigma_m$ | $\|\delta M\| \leq \frac{1}{2}\delta\Omega$ (Thm 2.2.11) | TLE-dependent |
| $i$ from TLE | $\sigma_m$ | $\|\delta M\| \leq \frac{1}{2}\delta i$ | TLE-dependent |
| $u$ from TLE | $\sigma_m$ | $\|\delta M\| \leq \frac{1}{2}\delta u$ | TLE-dependent |
| sin/cos evaluation | $\delta_p$ | $\leq \epsilon_{\mathrm{mach}}$ per entry | $\sim 10^{-16}$ |
| $SU(2)$ multiply rounding | $\delta_p$ | $\leq 8\epsilon_{\mathrm{mach}}$ per entry | $\sim 10^{-15}$ per multiply |
| $\theta_{\mathrm{GMST}}$ polynomial | $\delta_p$ | Angle sensitivity (Thm 2.2.11) | Ch 29 |
| GMST coefficients | $\sigma_m$ | Through $M_3(\theta)$ | IAU 1982 |
| $\omega_E$ | $\sigma_m$ | Through commutator (Cor 2.4.3) | $\sim 10^{-12}$ rad/s |
| TEME frame definition | $\delta_a$ | Irreducible floor | $\sim 0.1$ arcsec [A.2.1] |

---

## §2.6 Error Propagation Through the Framework

The orbit propagation pipeline consists of many operations: rotations, translations, corrections, solvers. Each operation propagates the three errors ($\sigma_m$, $\delta_p$, $\delta_a$) from its inputs to its outputs. This section establishes two principles that classify every operation in the pipeline and determine how its errors propagate.

### The Two Principles

The operations divide cleanly into two categories, each with a distinct error propagation rule.

**Theorem 2.6.1** (Principle 1: linear operations reuse the transform)**.** *For any linear operation $\mathcal{L}$ in the framework — conjugation $MVM^\dagger$ for fixed $M$, commutator $\frac{1}{2i}[\Omega, V]$ for fixed $\Omega$, addition, or scalar multiplication — the error transforms by the same operation:*

$$\delta(\mathcal{L}(V)) = \mathcal{L}(\delta(V)) \tag{2.62}$$

*This is exact, not a first-order approximation. The same computation that transforms the physical state, applied to the error state, gives the error of the output.*

*Proof.* For a linear map: $\mathcal{L}(V + \delta V) = \mathcal{L}(V) + \mathcal{L}(\delta V)$. Therefore $\mathcal{L}(V + \delta V) - \mathcal{L}(V) = \mathcal{L}(\delta V)$, which is the error of the output. This holds exactly for all linear maps: addition and scalar multiplication (Ch 1, §1.3), conjugation for fixed $M$ (Theorem 2.2.9), the commutator for fixed $\Omega$ (Theorem 2.4.5), and the adjoint action for fixed $\hat{M}_T$ (Theorem 2.4.2). ∎

**Theorem 2.6.2** (Principle 2: nonlinear operations require rigorous sensitivity bounds)**.** *For a nonlinear operation $g$, the error satisfies:*

$$|\delta(g(V))| \leq \sup|g'| \cdot |\delta(V)| \tag{2.63}$$

*where $|g'|$ is the sensitivity (derivative magnitude) of $g$ with respect to its input, bounded over the error interval. The supremum must be a proven upper bound — Ch 1, Theorem 1.4.1 (for real-valued functions) or its complex extension (for functions of complex arguments).*

*Proof.* Ch 1, Theorem 1.4.1 (Mean Value Theorem bound) states: for a continuously differentiable function $g$, $|g(x+\epsilon) - g(x)| \leq \sup_{|\xi-x| \leq |\epsilon|}|g'(\xi)| \cdot |\epsilon|$. Setting $\epsilon = \delta(V)$ and taking the supremum over the error interval gives (2.63). ∎

**Remark.** The distinction between the two principles corresponds to a natural division in the pipeline: the *application* of a transform to the state is linear (Principle 1), while the *construction* of the transform from parameters is generally nonlinear (Principle 2). Computing $M_3(\alpha)$ from the angle $\alpha$ involves sin and cos (nonlinear). But once $M_3(\alpha)$ is known, applying it to a vector via $M_3 V M_3^\dagger$ is linear in $V$. Most pipeline operations are applications — rotating vectors, adding corrections, composing transforms. The nonlinear steps are concentrated at the interfaces where physical parameters enter the algebra.

### Pipeline Classification

**Proposition 2.6.1** (Pipeline classification)**.** *Every operation in the orbit propagation pipeline is classified for error propagation:*

**Linear operations** (Principle 1 — error row uses the same transform):

| Operation | Expression | Why linear |
|-----------|-----------|-----------|
| Rotation of any vector | $MVM^\dagger$ | Linear in $V$ for fixed $M$ (Thm 2.2.9) |
| Translation | $+ T$ | Constant addition (Ch 1, Thm 1.3.1) |
| Transport term | $-\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ | Linear in $R_{\mathrm{pos}}$ for fixed $\Omega_T$ (Thm 2.4.5) |
| Adjoint action on twist | $\hat{M}_T\,\hat{\Omega}_s\,\hat{M}_T^{-1}$ | Linear in $\hat{\Omega}_s$ for fixed $\hat{M}_T$ (Thm 2.4.2) |
| Spatial→material conversion | $V_s - \frac{1}{2i}[\Omega_s, R_{\mathrm{pos}}]$ | Linear in $(V_s, \Omega_s, R_{\mathrm{pos}})$ separately (Cor 2.4.1) |
| Configuration composition | $\hat{M}_T \cdot \hat{M}$ | Linear in $\hat{M}$ for fixed $\hat{M}_T$ |
| Secular rate additions | $M_0 + \dot{M}\Delta t$ | Linear in both terms |
| Short/long-period corrections | $+ \delta V$ | Addition (Ch 1, Thm 1.3.1) |
| Dot product | $\frac{1}{2}\mathrm{tr}(AB)$ | Bilinear; linear in each argument (Thm 2.2.2) |
| Cross product | $\frac{1}{2i}[A, B]$ | Bilinear; linear in each argument (Thm 2.2.1) |

**Remark** (Hidden assumption)**.** The entries marked "linear" are linear in one operand when the other is exact. The $SU(2)$ composition $M_1 \cdot M_2$ is linear in $M_2$ for exact $M_1$ — the normal case (a known frame rotation applied to an uncertain state). When **both** factors carry errors, the bilinear cross-term requires Ch 1, Theorem 1.3.2. The same applies to dot and cross products with two uncertain operands.

**Nonlinear operations** (Principle 2 — error row requires a rigorous sensitivity bound):

| Operation | Sensitivity bound | Source |
|-----------|------------------|--------|
| Angle $\to M_k(\alpha)$ | $\|\partial M_k/\partial\alpha\| = \frac{1}{2}$ | Thm 2.2.11 |
| Kepler equation solver | $\leq 1/(1-e)$ | Ch 9 |
| sin/cos of orbital angles | $|\cos x|\,\delta + \delta^2/2$ | Ch 1, Cors. 1.4.1–1.4.2 |
| Renormalization $M/\sqrt{\det M}$ | $\leq 1 + O(\det M - 1)$ | Thm 2.6.3 |
| Cube root (element recovery) | $\frac{2}{3}x^{-1/3}$ | Ch 32 |
| Density model $((q_0-s)/(r-s))^\tau$ | $|\tau|\rho/(r-s)$ | Ch 21 |
| Euclidean norm $\sqrt{-\det V}$ | $\leq 1$ (Lipschitz-1) | Ch 1, Cor. 1.4.8 |

**Remark** (Cross product cancellation)**.** When $\mathbf{a} \approx k\mathbf{b}$ (nearly parallel vectors), the commutator $\frac{1}{2i}[A, B] \approx 0$ and the subtraction loses digits. The condition number is $\kappa \approx \|A\|\|B\| / \|\frac{1}{2i}[A,B]\|$, detected by $\mathrm{rd}(\mathbf{a} \times \mathbf{b}) \to 0$ (Ch 1, §1.7). This arises in angular momentum $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ when the orbit is nearly rectilinear.

### Renormalization

**Theorem 2.6.3** (Renormalization sensitivity)**.** *For $M$ near-unitary with $\det M = 1 + \varepsilon$ where $|\varepsilon| \ll 1$:*

$$M_{\mathrm{renorm}} = \frac{M}{\sqrt{\det M}} \approx M(1 - \varepsilon/2) \tag{2.64}$$

*The correction is $O(\varepsilon)$ and the sensitivity of the renormalization to the input is $\leq 1 + O(\varepsilon)$. After $n$ $SU(2)$ multiplications, $|\varepsilon| \leq n \cdot 8\epsilon_{\mathrm{mach}}$. For binary64 over a full propagation year ($\sim 2 \times 10^6$ compositions): $|\varepsilon| \leq 2 \times 10^{-9}$, far below measurement error floors. Renormalization is rarely needed at double precision or wider.*

*Proof.* $\sqrt{\det M} = \sqrt{1 + \varepsilon} = 1 + \varepsilon/2 + O(\varepsilon^2)$. Division: $M/(1 + \varepsilon/2) \approx M(1 - \varepsilon/2)$. Each entry is corrected by $\leq |\varepsilon|/2$. ∎

### Linearization Validity

**Remark** (When the parallel-row principle degrades)**.** Principle 1 (Theorem 2.6.1) is exact for linear operations regardless of error magnitude. Principle 2 (Theorem 2.6.2) uses a first-order sensitivity bound, which is valid when the error is small relative to the value — specifically, when $\mathrm{rd}(v) \geq d_{\min}$ (Ch 1, Definition 1.9.1). At $d_{\min} = 4$ (relative error $\leq 10^{-4}$), the neglected second-order term in Ch 1, Theorem 1.3.2 (the cross term $\delta_x\,\delta_y$) is $\leq 10^{-8}$ of the first-order terms — four orders of magnitude smaller. [A.2.2]

When $\mathrm{rd}(v) < d_{\min}$, the linearization may produce loose bounds. In this regime, the full nonlinear propagation rules of Ch 1 (Theorems 1.3.1–1.3.3 applied to each scalar operation) must be used, and a formula-switching event should be triggered (Ch 1, §1.9).

### Three-Error Independence

**Proposition 2.6.2** (Error categories propagate independently)**.** *The three error categories $\sigma_m$, $\delta_p$, $\delta_a$ propagate through all framework operations independently. Each category uses the same functional form — the same linear map (Principle 1) or the same sensitivity bound (Principle 2) — applied to different input magnitudes.*

*Proof.* Each propagation rule has the form $\delta_{\mathrm{out}} = g(|\mathrm{values}|, \delta_{\mathrm{in}})$ where $g$ is the same function regardless of which error category $\delta$ represents. The three categories are independent because they arise from disjoint physical mechanisms (Ch 1, Proposition 1.8.1): measurement from observations, precision from computation, accuracy from the model. ∎

---

## §2.7 The Adjoint Bridge: From $SU(2)$ to $3 \times 3$ Matrices

The $SU(2)$ framework of §§2.2–2.6 works entirely within $2 \times 2$ complex matrices. The existing orbit propagation literature — and the component-wise computations of downstream chapters — uses $3 \times 3$ real rotation matrices and $7 \times 7$ state matrices. The adjoint representation provides the exact bridge between these two forms: it converts an $SU(2)$ element $M$ into the $3 \times 3$ rotation matrix $R(M)$ that performs the same rotation on $\mathbb{R}^3$ vectors.

The conversion is a group homomorphism (Theorem 2.7.3): composing in $SU(2)$ and then converting gives the same result as converting each factor and composing in $SO(3)$. This means the two forms are interchangeable — any result proved in the $SU(2)$ framework automatically holds in the matrix framework, and vice versa.

### The Adjoint Representation

**Definition 2.7.1** (Adjoint rotation matrix)**.** *For $M \in SU(2)$, the $3 \times 3$ rotation matrix $R(M)$ is defined by:*

$$R_{jk}(M) = \frac{1}{2}\mathrm{tr}(\sigma_j\,M\,\sigma_k\,M^\dagger) \tag{2.65}$$

*The $(j,k)$ entry is the $j$-th Pauli coefficient of the rotated basis vector $M\sigma_k M^\dagger$ — that is, the $j$-th component of the image of $\hat{\mathbf{e}}_k$ under the rotation.*

**Theorem 2.7.1** (Adjoint equivalence)**.** *For any vector $\mathbf{v}$ encoded as $V = \sum_k v_k \sigma_k$:*

$$MVM^\dagger = \sum_j \left(\sum_k R_{jk}(M)\,v_k\right) \sigma_j \tag{2.66}$$

*The Pauli coefficients of the conjugated matrix are the components of $R(M)\,\mathbf{v}$. Conjugation by $M$ and left-multiplication by $R(M)$ produce the same rotated vector:*

$$\mathbf{v}' = R(M)\,\mathbf{v} \tag{2.67}$$

*Proof.* Expand: $MVM^\dagger = \sum_k v_k (M\sigma_k M^\dagger)$. Each $M\sigma_k M^\dagger$ is traceless Hermitian (Theorem 2.2.4), so it has a unique Pauli expansion: $M\sigma_k M^\dagger = \sum_j R_{jk}\,\sigma_j$ where $R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j\,M\sigma_k M^\dagger)$ (by orthogonality: $\frac{1}{2}\mathrm{tr}(\sigma_j \sigma_l) = \delta_{jl}$). Substituting: $MVM^\dagger = \sum_{jk} R_{jk}\,v_k\,\sigma_j$, and the coefficient of $\sigma_j$ is $(R\mathbf{v})_j$. ∎

**Theorem 2.7.2** ($R(M) \in SO(3)$)**.** *The matrix $R(M)$ is orthogonal with determinant $+1$.*

*Proof.* Orthogonality: since conjugation by $M$ preserves the inner product $\frac{1}{2}\mathrm{tr}(AB)$ (Theorem 2.2.2) and the norm $\sqrt{-\det V}$ (Theorem 2.2.3), the map $\mathbf{v} \mapsto R(M)\mathbf{v}$ is an isometry. Every isometry of $\mathbb{R}^3$ that fixes the origin is orthogonal. Determinant: $R(I) = I_3$ (the identity rotation), and $R(M)$ varies continuously from $I_3$ as $M$ varies from $I$. Since $\det$ is continuous and $\det(I_3) = +1$, the determinant cannot jump to $-1$. ∎

**Theorem 2.7.3** (Adjoint is a group homomorphism)**.** *$R(M_1 M_2) = R(M_1)\,R(M_2)$.*

*Proof.* $(M_1 M_2)\,V\,(M_1 M_2)^\dagger = M_1(M_2 V M_2^\dagger)M_1^\dagger$. The inner conjugation applies $R(M_2)$; the outer applies $R(M_1)$ to the result. ∎

**Remark.** Theorem 2.7.3 is the reason the two-form architecture works: composing in $SU(2)$ and converting once gives the same result as converting each factor and composing as $3 \times 3$ matrices. The $SU(2)$ product is cheaper (16 real multiplications for a $2 \times 2$ complex product vs. 27 for a $3 \times 3$ real product), and the conversion (Theorem 2.7.4) is done once per composed transform, not per vector.

### Explicit Entries

**Theorem 2.7.4** (Explicit $R(M)$)**.** *For $M = \left(\begin{smallmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{smallmatrix}\right)$ with quaternion components $\alpha = w + iz$, $\beta = ix + y$:*

$$R(M) = \begin{pmatrix} 1 - 2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\ 2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\ 2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2) \end{pmatrix} \tag{2.68}$$

*Every entry is a quadratic polynomial in the quaternion components $(w,x,y,z)$. This is the standard rotation matrix used throughout the orbit propagation literature.*

*Proof.* Compute $R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j M \sigma_k M^\dagger)$ for each $(j,k)$ by direct multiplication. For example: $R_{33} = \frac{1}{2}\mathrm{tr}(\sigma_3 M \sigma_3 M^\dagger)$. Since $\sigma_3 M \sigma_3 M^\dagger$ has trace $2(|\alpha|^2 - |\beta|^2) = 2(w^2+z^2-x^2-y^2)$, we get $R_{33} = w^2+z^2-x^2-y^2 = 1 - 2(x^2+y^2)$ (using $w^2+x^2+y^2+z^2 = 1$). The remaining entries follow analogously. ∎

**Remark.** Equivalently, in terms of $(\alpha, \beta)$ directly:

$$R(M) = \begin{pmatrix} \mathrm{Re}(\alpha^2 - \beta^2) & -\mathrm{Im}(\alpha^2 + \beta^2) & 2\,\mathrm{Im}(\alpha\bar{\beta}) \\ \mathrm{Im}(\alpha^2 - \beta^2) & \mathrm{Re}(\alpha^2 + \beta^2) & -2\,\mathrm{Re}(\alpha\bar{\beta}) \\ -2\,\mathrm{Im}(\alpha\beta) & 2\,\mathrm{Re}(\alpha\beta) & |\alpha|^2 - |\beta|^2 \end{pmatrix} \tag{2.69}$$

**Example 2.7.1** (Adjoint conversion for a $45°$ rotation about $\hat{\mathbf{e}}_3$)**.** *For $M_3(45°)$: $\alpha = e^{i\pi/8} = \cos 22.5° + i\sin 22.5° \approx 0.9239 + 0.3827i$, $\beta = 0$. Applying (2.68) with $(w,x,y,z) = (0.9239, 0, 0, 0.3827)$:*

*$R_{11} = 1 - 2(y^2 + z^2) = 1 - 2(0 + 0.1465) = 0.7071$*
*$R_{12} = 2(xy - wz) = 2(0 - 0.9239 \times 0.3827) = -0.7071$*
*$R_{13} = 2(xz + wy) = 0$*
*$R_{21} = 2(xy + wz) = 0.7071$*
*$R_{22} = 1 - 2(x^2 + z^2) = 0.7071$*
*$R_{33} = 1 - 2(x^2 + y^2) = 1$*

$$R(M_3(45°)) = \begin{pmatrix} 0.7071 & -0.7071 & 0 \\ 0.7071 & 0.7071 & 0 \\ 0 & 0 & 1 \end{pmatrix} = R_3(45°) \;\checkmark$$

*This is the standard $3 \times 3$ rotation matrix for $45°$ about $\hat{\mathbf{e}}_3$. The 9 entries required 18 multiplications and 9 additions from the quaternion components — computed once, then applied to any number of vectors by the $3 \times 3$ matrix-vector product.*

### Conversion Error

**Theorem 2.7.5** (Conversion precision)**.** *Evaluating (2.68) or (2.69) in finite arithmetic introduces entry-wise error:*

$$\|\delta(R)\|_{\max} \leq 4\,\epsilon_{\mathrm{mach}} \tag{2.70}$$

*Each entry is a sum of products of $(w,x,y,z)$ with $w^2+x^2+y^2+z^2 = 1$. Each product introduces $\leq \epsilon_{\mathrm{mach}}$ rounding; the most complex entry involves two products and a subtraction, giving $\leq 4\epsilon_{\mathrm{mach}}$.* [P.2.6]

**Theorem 2.7.6** (Angle sensitivity through the bridge)**.** *For an elementary rotation $M_k(\alpha)$:*

$$\left\|\frac{dR}{d\alpha}\right\| = 1 \tag{2.71}$$

*Therefore $\|\delta(R)\| \leq |\delta(\alpha)|$: the rotation matrix error equals the angle error. The factor-of-2 difference from the $SU(2)$ sensitivity ($\|\delta M\| \leq \frac{1}{2}|\delta\alpha|$, Theorem 2.2.11) arises because $R$ is quadratic in $M$ — the two-sided conjugation doubles the effect.*

*Proof.* $\frac{dR_{jk}}{d\alpha} = \frac{1}{2}\mathrm{tr}(\sigma_j \frac{dM}{d\alpha} \sigma_k M^\dagger + \sigma_j M \sigma_k \frac{dM^\dagger}{d\alpha})$. Since $\|\frac{dM_k}{d\alpha}\| = \frac{1}{2}$ (each entry is $\cos(\alpha/2)$ or $\sin(\alpha/2)$), and the trace formula involves two terms, the total sensitivity is $2 \times \frac{1}{2} = 1$. ∎

---

## §2.8 The Linear Representation

The dual quaternion pair $(\hat{M}, \hat{\Omega}_b)$ and the traditional matrix forms are fully equivalent representations of the same state and the same transforms. This section develops the conversion: each dual quaternion maps to a $4 \times 4$ matrix, and the two $4 \times 4$ matrices assemble into the $7 \times 7$ state matrix. The conversion proceeds through the adjoint bridge (§2.7) and the homogeneous coordinate trick.

### The $4 \times 4$ Homogeneous Matrix (from $\hat{M}$)

**Definition 2.8.1** (Homogeneous matrix)**.** *A rigid-body transform with rotation $R \in SO(3)$ and translation $\mathbf{t} \in \mathbb{R}^3$ is represented by the $4 \times 4$ matrix:*

$$H = \begin{pmatrix} R & \mathbf{t} \\ \mathbf{0}^T & 1 \end{pmatrix} \tag{2.72}$$

*The action on a 4-vector $(\mathbf{x}, w_H)^T$:*

$$H \begin{pmatrix} \mathbf{x} \\ w_H \end{pmatrix} = \begin{pmatrix} R\mathbf{x} + w_H\,\mathbf{t} \\ w_H \end{pmatrix} \tag{2.73}$$

*The fourth coordinate $w_H$ controls whether translation is applied. A position (a point in space) carries $w_H = 1$; a displacement (the difference between two points) carries $w_H = 0$.*

**Remark** (Affine vs. projective interpretation)**.** In projective geometry, a zero fourth coordinate represents a "point at infinity." Here $w_H$ is not a projective weight but an affine tag: $w_H = 1$ marks a position, $w_H = 0$ marks a displacement. The distinction arises naturally — the difference of two positions $(\tilde{\mathbf{r}}, 1) - (\mathbf{r}, 1) = (\delta\mathbf{r}, 0)$ carries $w_H = 0$. Unlike projective coordinates, scalar multiples are NOT equivalent: $(\delta\mathbf{r}, 0)$ and $(2\delta\mathbf{r}, 0)$ are distinct (different magnitudes). The subscript $H$ distinguishes the homogeneous coordinate from the quaternion scalar component $w$ (Definition 2.2.6).

**Example 2.8.1** (Worked: rotation + translation on state and displacement)**.** *A $90°$ rotation about $\hat{\mathbf{e}}_3$ followed by a $100$ km translation along $\hat{\mathbf{e}}_1$. Satellite at $\mathbf{r} = (7000, 0, 0)$ km with displacement $\delta\mathbf{r} = (0.5, 0.3, 0.1)$ km.*

$$H = \begin{pmatrix} 0 & -1 & 0 & 100 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

*State ($w_H = 1$):*
$H \cdot (7000, 0, 0, 1)^T = (100, 7000, 0, 1)^T$ — rotated and translated. ✓

*Displacement ($w_H = 0$):*
$H \cdot (0.5, 0.3, 0.1, 0)^T = (-0.3, 0.5, 0.1, 0)^T$ — rotated only, translation killed by $w_H = 0$. ✓

*Verification:* Transform the true position $\tilde{\mathbf{r}} = (7000.5, 0.3, 0.1)$:
$H \cdot (7000.5, 0.3, 0.1, 1)^T = (99.7, 7000.5, 0.1, 1)^T$.
Difference: $(99.7 - 100, 7000.5 - 7000, 0.1 - 0, 1-1) = (-0.3, 0.5, 0.1, 0)$. Matches the displacement output exactly. ∎

**Remark** (Signed displacements vs. non-negative bounds)**.** Example 2.8.1 uses a signed displacement $\delta\mathbf{r} = (0.5, 0.3, 0.1)$. The rotation produces a signed output $(-0.3, 0.5, 0.1)$, which is correct for the actual perturbation. For non-negative error bounds ($\sigma_m(r_1) \geq 0$, etc.), rotation can produce negative components, which are meaningless as bounds. The correct propagation for non-negative bounds uses $|R|$ (the matrix of absolute values) instead of $R$. For scalar error norms $\delta_{\mathrm{pos}} = \|\boldsymbol{\delta}(\mathbf{r})\|$, the propagation is trivial: rotation preserves the norm. In Example 2.8.1: $\|\delta\mathbf{r}\| = \sqrt{0.35} \approx 0.592$ km before and after.

**Theorem 2.8.1** (Homogeneous composition)**.** *$H_1 H_2 = \left(\begin{smallmatrix} R_1 R_2 & R_1\mathbf{t}_2 + \mathbf{t}_1 \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$.*

*Proof.* Block multiplication of the $2 \times 2$ block structure:

$$H_1 H_2 = \begin{pmatrix} R_1 & \mathbf{t}_1 \\ \mathbf{0}^T & 1 \end{pmatrix}\begin{pmatrix} R_2 & \mathbf{t}_2 \\ \mathbf{0}^T & 1 \end{pmatrix} = \begin{pmatrix} R_1 R_2 + \mathbf{t}_1\mathbf{0}^T & R_1\mathbf{t}_2 + \mathbf{t}_1 \cdot 1 \\ \mathbf{0}^T R_2 + 1 \cdot \mathbf{0}^T & \mathbf{0}^T\mathbf{t}_2 + 1 \end{pmatrix} = \begin{pmatrix} R_1 R_2 & R_1\mathbf{t}_2 + \mathbf{t}_1 \\ \mathbf{0}^T & 1 \end{pmatrix}$$

The rotation composes as $R_1 R_2$. The translation $R_1\mathbf{t}_2 + \mathbf{t}_1$ rotates the second translation by the first rotation, then adds the first translation — the same composition rule as the dual $SU(2)$ product (Theorem 2.3.2, equation 2.28). ∎

**Theorem 2.8.2** (Homogeneous inverse)**.** *$H^{-1} = \left(\begin{smallmatrix} R^T & -R^T\mathbf{t} \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$.*

*Proof.* Verify: $R R^T = I_3$, $R(-R^T\mathbf{t}) + \mathbf{t} = \mathbf{0}$. ∎

### The $4 \times 4$ Twist Matrix (from $\hat{\Omega}_b$)

The velocity dual quaternion $\hat{\Omega}_b = \Omega_b + \varepsilon V_b$ maps to the $4 \times 4$ twist matrix — the standard matrix representation of a Lie algebra element $\mathfrak{se}(3)$:

**Definition 2.8.2a** (Twist matrix)**.** *For $\hat{\Omega}_b$ with angular velocity $\boldsymbol{\omega}_b$ (decoded from $\Omega_b$) and translational velocity $\mathbf{v}_b$ (decoded from $V_b$):*

$$\Xi(\hat{\Omega}_b) = \begin{pmatrix} [\boldsymbol{\omega}_b]_\times & \mathbf{v}_b \\ \mathbf{0}^T & 0 \end{pmatrix}$$

*where $[\boldsymbol{\omega}_b]_\times$ is the $3 \times 3$ skew-symmetric matrix (Definition 2.8.2). The twist matrix is the Cartesian equivalent of $\hat{\Omega}_b$: the angular velocity maps through $\Omega_b \to [\boldsymbol{\omega}_b]_\times$ (the Lie algebra isomorphism $\mathfrak{su}(2) \to \mathfrak{so}(3)$), and the translational velocity decodes directly.*

**Remark** (Conversion summary)**.** The two dual quaternions convert to their $4 \times 4$ linear equivalents by:

| Dual quaternion | Type | $4 \times 4$ linear form | Conversion |
|---|---|---|---|
| $\hat{M}$ (configuration) | Unit dual, $SE(3)$ | $H = \left(\begin{smallmatrix} R(M) & \mathbf{t} \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$ | $R(M)$ via Thm 2.7.4; $\mathbf{t}$ via $2DM^\dagger$ |
| $\hat{\Omega}_b$ (velocity) | Pure dual, $\mathfrak{se}(3)$ | $\Xi = \left(\begin{smallmatrix} [\boldsymbol{\omega}]_\times & \mathbf{v} \\ \mathbf{0}^T & 0 \end{smallmatrix}\right)$ | $[\boldsymbol{\omega}]_\times$ via Lie algebra map; $\mathbf{v}$ decoded directly |

### The $7 \times 7$ State Matrix

The $4 \times 4$ homogeneous matrix handles the configuration. The twist matrix handles velocity. The $7 \times 7$ state matrix combines both into a single linear operator on $(\mathbf{r}, \mathbf{v}, 1)^T$, incorporating the transport coupling between position and material velocity (Corollary 2.4.2, equation 2.48).

**Definition 2.8.2** (Skew-symmetric matrix)**.** *For $\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3)$:*

$$[\boldsymbol{\omega}]_\times = \begin{pmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{pmatrix} \tag{2.74}$$

*satisfying $[\boldsymbol{\omega}]_\times \mathbf{v} = \boldsymbol{\omega} \times \mathbf{v}$. This is the Cartesian form of the commutator $\frac{1}{2i}[\Omega, V]$ (Theorem 2.2.1).*

**Definition 2.8.3** (The $7 \times 7$ state matrix)**.** *A rotating-frame transform with rotation $R = R(M_T)$, angular velocity $\boldsymbol{\omega}_T$, translation $\mathbf{t}$, and velocity offset $\Delta\mathbf{v}$:*

$$T_7 = \begin{pmatrix} R & \mathbf{0} & \mathbf{t} \\ -[\boldsymbol{\omega}_T]_\times R & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \tag{2.75}$$

*acting on the state vector $(\mathbf{r}, \mathbf{v}, 1)^T$:*

$$\begin{pmatrix} \mathbf{r}' \\ \mathbf{v}' \\ 1 \end{pmatrix} = T_7 \begin{pmatrix} \mathbf{r} \\ \mathbf{v} \\ 1 \end{pmatrix} = \begin{pmatrix} R\mathbf{r} + \mathbf{t} \\ -[\boldsymbol{\omega}_T]_\times R\mathbf{r} + R\mathbf{v} + \Delta\mathbf{v} \\ 1 \end{pmatrix} \tag{2.76}$$

The three velocity contributions are:
- $R\mathbf{v}$: velocity rotated into the new frame
- $-[\boldsymbol{\omega}_T]_\times R\mathbf{r}$: the transport correction (subtracts the frame's rotational velocity)
- $\Delta\mathbf{v}$: velocity offset (zero in all standard pipeline transforms)

**Theorem 2.8.3** (Equivalence to dual quaternion state transform)**.** *When $\mathbf{t} = 0$ (the standard orbit propagation case), equation (2.76) is the Cartesian-component form of the material velocity transform (Corollary 2.4.2, equation 2.48). The position row is the configuration transform (equation 2.44) in Cartesian coordinates. When $\mathbf{t} \neq 0$, the transport term in (2.76) acts on $R\mathbf{r}$ rather than $\mathbf{r}' = R\mathbf{r} + \mathbf{t}$; the missing coupling $\boldsymbol{\omega}_T \times \mathbf{t}$ must be absorbed into $\Delta\mathbf{v}$. Specifically:*
- *The position row $\mathbf{r}' = R\mathbf{r} + \mathbf{t}$ is the configuration transform (2.44) in Cartesian coordinates.*
- *The velocity row $\mathbf{v}' = R\mathbf{v} - [\boldsymbol{\omega}_T]_\times R\mathbf{r} + \Delta\mathbf{v}$ is the material velocity (2.48) with $-\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$ written as $-[\boldsymbol{\omega}_T]_\times\mathbf{r}'$ (with $\mathbf{t} = 0$, giving $-[\boldsymbol{\omega}_T]_\times R\mathbf{r}$).*

*Proof.* The equivalence $MVM^\dagger \leftrightarrow R(M)\mathbf{v}$ is Theorem 2.7.1. The equivalence $\frac{1}{2i}[\Omega, V] \leftrightarrow [\boldsymbol{\omega}]_\times\mathbf{r}$ follows from Theorem 2.2.1 (cross product as commutator) written in components. ∎

**Lemma 2.8.1** (Rotation of a cross product)**.** *For $R \in SO(3)$: $R[\boldsymbol{\omega}]_\times R^T = [R\boldsymbol{\omega}]_\times$.*

*Proof.* For any $\mathbf{v}$: $R[\boldsymbol{\omega}]_\times R^T \mathbf{v} = R(\boldsymbol{\omega} \times R^T\mathbf{v}) = R\boldsymbol{\omega} \times RR^T\mathbf{v} = (R\boldsymbol{\omega}) \times \mathbf{v}$, where the second equality uses $R(\mathbf{a} \times \mathbf{b}) = (R\mathbf{a}) \times (R\mathbf{b})$ for $R \in SO(3)$ (this holds because $\det R = +1$ preserves the orientation of the cross product). Since this holds for all $\mathbf{v}$: $R[\boldsymbol{\omega}]_\times R^T = [R\boldsymbol{\omega}]_\times$. ∎

**Theorem 2.8.4** ($7 \times 7$ composition)**.** *$T_{7,1} \cdot T_{7,2}$ is a $7 \times 7$ matrix of the same form with:*

$$R_{12} = R_1 R_2, \qquad \boldsymbol{\omega}_{T,12} = \boldsymbol{\omega}_{T,1} + R_1\,\boldsymbol{\omega}_{T,2}$$

$$\mathbf{t}_{12} = R_1\mathbf{t}_2 + \mathbf{t}_1, \qquad \Delta\mathbf{v}_{12} = -[\boldsymbol{\omega}_{T,1}]_\times R_1\mathbf{t}_2 + R_1\Delta\mathbf{v}_2 + \Delta\mathbf{v}_1$$

*Proof.* Block multiplication of two matrices in the form (2.75). The lower-left block: $(-[\boldsymbol{\omega}_{T,1}]_\times R_1)R_2 + R_1(-[\boldsymbol{\omega}_{T,2}]_\times R_2) = -([\boldsymbol{\omega}_{T,1}]_\times + R_1[\boldsymbol{\omega}_{T,2}]_\times R_1^T)R_{12}$. By Lemma 2.8.1: $R_1[\boldsymbol{\omega}_{T,2}]_\times R_1^T = [R_1\boldsymbol{\omega}_{T,2}]_\times$. The lower-left block is $-[\boldsymbol{\omega}_{T,1} + R_1\boldsymbol{\omega}_{T,2}]_\times R_{12}$, confirming $\boldsymbol{\omega}_{T,12} = \boldsymbol{\omega}_{T,1} + R_1\boldsymbol{\omega}_{T,2}$ with the correct minus sign. ∎

**Remark** (Two equivalent representations)**.** The dual quaternion pair $(\hat{M}, \hat{\Omega}_b)$ and the matrix pair $(H, \Xi)$ / the $7 \times 7$ matrix $T_7$ are fully equivalent representations of the same state and transforms. The dual quaternion form is preferable for: singularity analysis (Proposition 2.5.1), compact composition (equation 2.60), formal error architecture (§2.6, Principle 1), the kinematic equation (Theorem 2.4.1), and algebraic manipulation. The matrix form is preferable for: connection to existing orbit propagation literature, component-wise computations that downstream chapters require (Ch 30), and contexts where a linear operator on a state vector is needed. Both forms produce identical results — the matrix form IS the dual quaternion form expanded in Cartesian coordinates. The code-to-theorem mapping for specific source modules is in Appendix C.

---

## §2.9 The Inverse Bridge: Shepperd's Method

The adjoint bridge (§2.7) converts $M \to R(M)$. The inverse — recovering $M$ from a known $3 \times 3$ rotation matrix $R$ — is needed when converting from $3 \times 3$ matrix representations or from the output of other propagation software.

**Theorem 2.9.1** (Shepperd's method)**.** *Given $R \in SO(3)$, the $SU(2)$ element $M = \left(\begin{smallmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{smallmatrix}\right)$ is recovered by:*

1. *Compute the four candidates $q_k^2$:*

$$q_0^2 = \tfrac{1}{4}(1 + R_{11} + R_{22} + R_{33}), \quad q_1^2 = \tfrac{1}{4}(1 + R_{11} - R_{22} - R_{33})$$

$$q_2^2 = \tfrac{1}{4}(1 - R_{11} + R_{22} - R_{33}), \quad q_3^2 = \tfrac{1}{4}(1 - R_{11} - R_{22} + R_{33}) \tag{2.77}$$

2. *Select the largest $q_k^2$ (ensuring the square root argument is $\geq \frac{1}{4}$ for numerical stability).*

3. *Compute the remaining components from the off-diagonal entries of $R$:*

$$q_0 q_1 = \tfrac{1}{4}(R_{32} - R_{23}), \quad q_0 q_2 = \tfrac{1}{4}(R_{13} - R_{31}), \quad q_0 q_3 = \tfrac{1}{4}(R_{21} - R_{12})$$

$$q_1 q_2 = \tfrac{1}{4}(R_{12} + R_{21}), \quad q_1 q_3 = \tfrac{1}{4}(R_{13} + R_{31}), \quad q_2 q_3 = \tfrac{1}{4}(R_{23} + R_{32}) \tag{2.78}$$

4. *Divide each $q_j q_k$ by the selected $q_k$ to obtain the remaining components.* [P.2.1]

*The quaternion $(q_0, q_1, q_2, q_3) = (w, x, y, z)$ maps to $M$ via $\alpha = w + iz$, $\beta = ix + y$.*

*Proof.* Invert the quadratic expressions of Theorem 2.7.4: (2.77) follows from sums and differences of diagonal entries of $R$, using $w^2+x^2+y^2+z^2 = 1$. The products (2.78) follow from sums and differences of off-diagonal entries. Selecting the largest $q_k^2$ ensures the divisor in step 4 is maximized, minimizing division error. ∎

**Remark.** Shepperd's method avoids the singularity of the classical atan2-based extraction (which fails when $q_0 \approx 0$, i.e., rotation near $180°$). The branch selection guarantees a well-conditioned square root (argument $\geq \frac{1}{4}$) for any rotation.

---

## §2.10 The Parallel-Row Error Architecture

The error propagation principles of §2.6 lead to a simple computational architecture: the physical state and the error state have the same structure, and for linear operations they are transformed by the same operation.

### The Four-Row Structure

**Definition 2.10.1** (Parallel-row state)**.** *The complete tracked state consists of four rows, each carrying the dual quaternion pair $(\hat{M}, \hat{\Omega}_b)$:*

| Row | Content | Structure | Parameters |
|-----|---------|-----------|------------|
| 1 | Physical state | $(\hat{M}, \hat{\Omega}_b)$ | $M$ (4) + $R_{\mathrm{pos}}$ (3) + $\Omega_b$ (3) + $V_b$ (3) = 13 |
| 2 | Measurement error $\sigma_m$ | $(\delta\hat{M}, \delta\hat{\Omega})_{\sigma_m}$ | 13 |
| 3 | Precision error $\delta_p$ | $(\delta\hat{M}, \delta\hat{\Omega})_{\delta_p}$ | 13 |
| 4 | Accuracy error $\delta_a$ | $(\delta\hat{M}, \delta\hat{\Omega})_{\delta_a}$ | 13 |

*Total storage: 4 $\times$ 13 = 52 real numbers. All four rows pass through the same transform pipeline. For linear operations (the majority — including the adjoint action on $\hat{\Omega}_s$, Theorem 2.4.2), the same operation processes each row. For nonlinear operations, the error rows are multiplied by the rigorous sensitivity bound (Theorem 2.6.2).*

**Remark** (Merged mode)**.** When the error category breakdown is not needed, the three error rows merge into a single total error row: $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$ per component (Ch 1, Definition 1.8.1). This reduces storage to 2 $\times$ 13 = 26 real numbers and requires one extra transform pass instead of three.

### Scalar Error Summaries

**Definition 2.10.2** (Scalar error reduction)**.** *The per-component errors in each row are reduced to four scalar summaries:*

$$\delta(\mathrm{position}) = \sqrt{\delta r_1^2 + \delta r_2^2 + \delta r_3^2} = \sqrt{-\det(\delta R_{\mathrm{pos}})} \tag{2.79}$$

$$\delta(\mathrm{velocity}) = \sqrt{\delta v_1^2 + \delta v_2^2 + \delta v_3^2} = \sqrt{-\det(\delta V_{\mathrm{vel}})} \tag{2.80}$$

$$\delta(\mathrm{orientation}) = \|\boldsymbol{\varepsilon}\| \;\text{where}\; \delta M = \tfrac{i}{2}\varepsilon_k\sigma_k \cdot M \;\text{(perturbation angle in radians)} \tag{2.81}$$

$$\delta(\mathrm{angular\,velocity}) = \sqrt{\delta\omega_1^2 + \delta\omega_2^2 + \delta\omega_3^2} = \sqrt{-\det(\delta\Omega)} \tag{2.82}$$

*The three vector quantities use the same algebraic operation: $\sqrt{-\det(\cdot)}$ applied to a traceless Hermitian error matrix (Theorem 2.2.3). The orientation error extracts the perturbation angle from the $SU(2)$ tangent vector.*

**Remark** (Connection to reliable digits)**.** The scalar summaries feed into Ch 1, Definition 1.9.1:

$$\mathrm{rd}(\mathrm{position}) = \left\lfloor -\log_{10}\!\left(\frac{\delta(\mathrm{position})}{\|\mathbf{r}\|}\right) \right\rfloor$$

This is the decision criterion for formula switching (Ch 1, §1.9) and for output quality assessment (Ch 38).

### Point-Mass Reduction

**Corollary 2.10.1** (Point-mass mode)**.** *For the standard propagator with no body attitude tracking (Corollary 2.4.4):*
- *$\hat{\Omega}_b = \varepsilon V_b$ (angular velocity zero, velocity dual reduces to pure translation)*
- *The spatial twist velocity equals the material velocity (Corollary 2.4.1)*
- *The orientation $M$ carries the frame rotation only; its error reflects TLE angle uncertainty*
- *Active error components: 4 (orientation) + 3 (position) + 3 (velocity) = 10 per error row*
- *Total storage: 13 (physical) + 3 $\times$ 10 (active errors) = 43 real numbers, or 23 in merged mode*

---

## §2.11 Summary: How Subsequent Chapters Use This Framework

Every subsequent chapter that performs a coordinate transform or state propagation includes a "State Framework" section showing the operation in both representations — the dual quaternion pair $(\hat{M}, \hat{\Omega}_b)$ and the $7 \times 7$ matrix — following this pattern:

1. Express the operation as an $SU(2)$ transform or traceless Hermitian matrix operation (§§2.2–2.4)
2. Classify it as linear or nonlinear for error propagation (§2.6, Proposition 2.6.1)
3. For linear operations: state that Principle 1 (Theorem 2.6.1) applies — the error row uses the same transform
4. For nonlinear operations: provide the rigorous sensitivity bound (Theorem 2.6.2) with proof
5. Show how the operation composes with adjacent pipeline steps (§2.4, Theorem 2.4.1)
6. When the $3 \times 3$ / $7 \times 7$ matrix form is needed: derive it via the adjoint bridge (§2.7) and the homogeneous lift (§2.8)

### Downstream Chapter Reference

| Chapter | Framework aspect used | This chapter's section |
|---------|---------------------|----------------------|
| Ch 8 (Keplerian orbit) | Perifocal frame vectors; rotation to TEME | §2.5, Thm 2.5.1 |
| Ch 9 (Kepler's equation) | Nonlinear solver sensitivity bound | §2.6, Thm 2.6.2 |
| Ch 13 (Geopotential) | Force gradient as traceless Hermitian perturbation to $\hat{\Omega}$ | §2.4, Def 2.4.1 |
| Ch 16–17 (Brouwer secular) | Secular rates as near-identity perturbations | §2.4, Thm 2.4.6 |
| Ch 18–19 (Short/long-period) | Additive corrections (linear, Principle 1) | §2.6, Thm 2.6.1 |
| Ch 20 (Osculating elements) | Complete mean-to-osculating transform | §2.5, Thm 2.5.3 |
| Ch 27 (Third-body) | Third-body force as traceless Hermitian perturbation | §2.4 |
| Ch 29 (Sidereal time) | GMST as input to TEME→PEF rotation | §2.5, Thm 2.5.2 |
| Ch 30 (Coordinate transforms) | All frame rotations and compositions | §§2.2, 2.3, 2.5, 2.7–2.8 |
| Ch 32 (Element recovery) | Cube root sensitivity (nonlinear, Principle 2) | §2.6, Thm 2.6.2 |
| Ch 34–35 (Propagation) | Full pipeline as chain of $SU(2)$ products | §2.5, Thm 2.5.3 |
| Ch 36 (Architecture) | Matched-pair compatibility with $SU(2)$ | Ch 3, §3.5 (Def. 3.5.1), §3.7 (Def. 3.7.1) |
| Ch 38 (State vector output) | Error budget from parallel-row architecture | §§2.6, 2.10 |

---

## Error Notes

**[P.2.1]** Shepperd branch selection. The square root in Theorem 2.9.1 has argument $\geq \frac{1}{4}$ (the largest of four non-negative quantities summing to $1$ is $\geq \frac{1}{4}$). The division in step 4 has denominator $\geq \frac{1}{2}$. The precision error per component is bounded by $\sqrt{2}\,\epsilon_{\mathrm{mach}}$ (Ch 1, Corollary 1.4.4 for the square root, Ch 1, Theorem 1.3.3 for the divisions). *Remedy:* use wider arithmetic type.

**[P.2.2]** Inexact rotation entry error. Each entry of $M_k(\alpha)$ is $\cos(\alpha/2)$ or $\sin(\alpha/2)$, evaluated to precision $\epsilon_{\mathrm{mach}}$ (Ch 1, Corollaries 1.4.1–1.4.2). The entry-wise error: $\|\delta(M)\|_{\max} \leq \epsilon_{\mathrm{mach}}$. *Remedy:* none needed at double precision or wider.

**[P.2.3]** Elementary rotation rounding. The half-angle parameterization provides a factor-of-2 bonus: $\|\delta M\| \leq \frac{1}{2}|\delta\alpha|$ (Theorem 2.2.11). The rounding from evaluating $\cos(\alpha/2)$ and $\sin(\alpha/2)$ is proportional to $\epsilon_{\mathrm{mach}}$, independent of $\alpha$. *Remedy:* use wider arithmetic type.

**[P.2.4]** Composed rotation rounding. Two $SU(2)$ multiplications accumulate $\leq 16\epsilon_{\mathrm{mach}}$ total rounding (Corollary 2.2.1). This is smaller than angle measurement errors by many orders of magnitude for any practical TLE precision. *Remedy:* none needed.

**[P.2.5]** Position extraction error coupling. Extracting position via $T = 2DM^\dagger$ (Theorem 2.3.6) couples rotation error through position magnitude: $\|\delta(\mathbf{r})\| \leq 2\|\delta(D)\| + \|\mathbf{r}\|\,\|\delta(M)\|$. *Remedy:* use wider arithmetic type for the $SU(2)$ entries.

**[P.2.6]** Adjoint conversion rounding. The 9 entries of $R(M)$ are quadratic in the quaternion components (Theorem 2.7.4). Each entry requires $\leq 4$ floating-point operations on values with $|w|, |x|, |y|, |z| \leq 1$. Entry-wise error: $\leq 4\epsilon_{\mathrm{mach}}$ (Theorem 2.7.5). *Remedy:* use wider arithmetic type.

**[A.2.1]** TEME frame accuracy. The TEME frame differs from precise GCRS/ITRS frames by $\sim 0.1$ arcsec ($\sim 5 \times 10^{-7}$ rad) due to the truncated precession/nutation model. This is an irreducible accuracy error. *Remedy:* apply IAU 2006/2000A corrections (beyond the standard model; see Ch 3).

**[A.2.2]** Linearization in error propagation. Principle 2 (Theorem 2.6.2) uses a first-order sensitivity bound. The neglected second-order terms are $O(\delta^2)$. At $\mathrm{rd} \geq 4$ (relative error $\leq 10^{-4}$), these terms are $\leq 10^{-8}$ of the first-order bound. When $\mathrm{rd} < 4$, the full nonlinear rules of Ch 1 must be used. *Remedy:* fall back to Ch 1 scalar propagation rules when reliable-digits criterion is not met.

**[M.2.1]** GMST angle measurement. The GMST polynomial (Ch 29) determines $\theta_{\mathrm{GMST}}$. Measurement error in its coefficients propagates through Theorem 2.2.11: $\|\delta(\mathbf{r}_{\mathrm{PEF}})\| \leq \|\delta(\mathbf{r}_{\mathrm{TEME}})\| + \|\mathbf{r}\| \cdot |\delta(\theta_{\mathrm{GMST}})|$. *Remedy:* use updated IERS values (but this may break the matched pair — see Ch 3).
