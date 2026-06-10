# Chapter 2: The State Framework

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\sigma_k$ | Pauli matrices ($k = 1, 2, 3$) | §2.2, Def. 2.2.1 |
| $I$ | The $2 \times 2$ identity matrix | §2.2 |
| $M$ | An element of $SU(2)$: a $2 \times 2$ unitary matrix with $\det M = 1$ | §2.2, Def. 2.2.3 |
| $A^\dagger$ | Conjugate transpose (Hermitian adjoint): $(A^\dagger)_{jk} = \overline{a_{kj}}$ | §2.2, Def. 2.2.1a |
| $\alpha, \beta$ | Complex parameters of $M$ satisfying $M = \left(\begin{smallmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{smallmatrix}\right)$ with $\|\alpha\|^2+\|\beta\|^2=1$ | §2.2, Def. 2.2.3 |
| $V$ | A $2 \times 2$ traceless Hermitian matrix encoding a 3D vector (Def. 2.2.2, Def. 2.2.3) | §2.2 |
| $[A, B]$ | Matrix commutator $AB - BA$ | §2.2 |
| $M_k(\alpha)$ | Elementary rotation about axis $k$ by angle $\alpha$ | §2.2, Def. 2.2.5 |
| $\hat{M}$ | A dual $SU(2)$ matrix: $\hat{M} = M + \varepsilon D$, encoding rotation + translation | §2.3, Def. 2.3.2 |
| $\varepsilon$ | The dual unit: $\varepsilon^2 = 0$, $\varepsilon \neq 0$ | §2.3, Def. 2.3.1 |
| $D$ | Dual part of $\hat{M}$, encoding translation: $\mathbf{t} = 2D M^\dagger$ (vector part) | §2.3 |
| $R_{\mathrm{pos}}$ | Position vector encoded as traceless Hermitian matrix | §2.4 |
| $V_{\mathrm{vel}}$ | Velocity vector encoded as traceless Hermitian matrix | §2.4 |
| $\Omega$ | Angular velocity encoded as traceless Hermitian matrix | §2.4 |
| $R(q)$ | The $3 \times 3$ rotation matrix derived from the SU(2) element | §2.7 |
| $\omega_E$ | Earth's rotation rate ($\approx 7.292 \times 10^{-5}$ rad/s) | §2.5 |
| $\theta_{\mathrm{GMST}}$ | Greenwich Mean Sidereal Time angle | §2.5 |
| $\Omega_{\mathrm{node}}, i, u$ | Right ascension of ascending node, inclination, argument of latitude | §2.5 |

---

## §2.1 Introduction

Orbit propagation produces a satellite's position and velocity in one reference frame. Comparing with observations, planning maneuvers, or computing ground tracks requires transforming this state into other frames. Each transformation involves rotation (the frames are oriented differently), translation (the frames have different origins), and — for rotating frames — velocity coupling (the transport theorem adds $\boldsymbol{\omega} \times \mathbf{r}$ to the observed velocity).

The traditional approach handles these three aspects separately: a $3 \times 3$ rotation matrix rotates vectors, a translation vector is added, and the transport term is computed from a cross product. Composing multiple transforms requires tracking which matrix, which offset, and which angular velocity applies at each stage — a bookkeeping exercise that grows with the number of frames in the pipeline.

This chapter develops a unified algebraic framework in which a single mathematical object encodes rotation, translation, and velocity coupling simultaneously. The framework has three properties that make it the foundation for all subsequent chapters:

1. **Singularity-free.** The representation of orientation has no coordinate singularity (no gimbal lock), eliminating the need for the Lyddane modification and similar patches at special inclinations or eccentricities.

2. **Compositional.** Two transforms compose by a single algebraic product, producing a third transform of the same type. The entire orbit propagation pipeline — from orbital-plane coordinates through inertial frame to Earth-fixed frame — is a chain of such products.

3. **Error-compatible.** Every transform in the framework is a linear operation on the state. By the composition principle of Ch. 1 (Theorem 1.3.4), the error state transforms by the same operation as the physical state, with no separate error formulas needed for linear steps. The three error categories ($\sigma_m$, $\delta_p$, $\delta_a$) propagate independently through the same transform.

The algebraic foundation is the group $SU(2)$ of $2 \times 2$ unitary matrices with unit determinant, together with its dual extension for incorporating translation. Every 3D vector — position, velocity, angular velocity — is encoded as a $2 \times 2$ traceless Hermitian matrix (a matrix equal to its own conjugate transpose with zero diagonal sum; Definition 2.2.2), and every rotation acts by matrix conjugation: $V' = M V M^\dagger$. The cross product, dot product, and norm are native operations in this matrix algebra (commutator, trace, and determinant respectively), so the transport theorem and error propagation emerge from standard matrix calculus rather than vector-specific formulas.

The traditional $3 \times 3$ rotation matrix and $7 \times 7$ state matrix are derived in §2.7 as component-wise expansions of the $SU(2)$ operations, connecting the framework to the existing orbit propagation literature.

---

## §2.2 The Algebra of Rotation

### Definitions

**Definition 2.2.1** (Pauli matrices)**.** *The three Pauli matrices are:*

$$\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} \tag{2.1}$$

*Together with the identity $I$, they form a basis for the real vector space of $2 \times 2$ Hermitian matrices.*

**Definition 2.2.1a** (Conjugate transpose)**.** *For a matrix $A \in \mathbb{C}^{m \times n}$ with entries $a_{jk}$, the conjugate transpose (Hermitian adjoint) $A^\dagger \in \mathbb{C}^{n \times m}$ is defined by:*

$$(A^\dagger)_{jk} = \overline{a_{kj}} \tag{2.1b}$$

*That is: transpose the matrix and take the complex conjugate of every entry. For a $2 \times 2$ matrix:*

$$\begin{pmatrix} a & b \\ c & d \end{pmatrix}^\dagger = \begin{pmatrix} \bar{a} & \bar{c} \\ \bar{b} & \bar{d} \end{pmatrix}$$

*Key properties: $(AB)^\dagger = B^\dagger A^\dagger$ (reversal of order), $(A^\dagger)^\dagger = A$, and $\det(A^\dagger) = \overline{\det(A)}$. For real matrices, $A^\dagger = A^T$ (ordinary transpose).*

**Definition 2.2.2** (Hermitian and traceless Hermitian matrices)**.** *A matrix $A \in \mathbb{C}^{2 \times 2}$ is Hermitian if $A = A^\dagger$ (equal to its own conjugate transpose). In components: the diagonal entries are real and the off-diagonal entries are complex conjugates of each other. A Hermitian matrix is traceless if additionally $\mathrm{tr}(A) = a_{11} + a_{22} = 0$. Explicitly, every $2 \times 2$ traceless Hermitian matrix has the form:*

$$A = \begin{pmatrix} a & b \\ \bar{b} & -a \end{pmatrix} \tag{2.1a}$$

*where $a \in \mathbb{R}$ and $b \in \mathbb{C}$. This is a 3-dimensional real vector space, with the Pauli matrices $\sigma_1, \sigma_2, \sigma_3$ as a basis.*

**Definition 2.2.3** (Vector encoding)**.** *A vector $\mathbf{v} = (v_1, v_2, v_3) \in \mathbb{R}^3$ is encoded as the traceless Hermitian matrix:*

$$V = v_1 \sigma_1 + v_2 \sigma_2 + v_3 \sigma_3 = \begin{pmatrix} v_3 & v_1 - iv_2 \\ v_1 + iv_2 & -v_3 \end{pmatrix} \tag{2.2}$$

*This is an isomorphism: every traceless Hermitian $2 \times 2$ matrix corresponds to exactly one vector in $\mathbb{R}^3$, and vice versa. The components are recovered by $v_k = \frac{1}{2}\mathrm{tr}(\sigma_k V)$.*

**Definition 2.2.4** (Rotation element)**.** *A rotation is an element $M \in SU(2)$: a $2 \times 2$ unitary matrix with unit determinant:*

$$M^\dagger M = I, \qquad \det M = 1 \tag{2.3}$$

*The general form, parameterized by $\alpha, \beta \in \mathbb{C}$ with $|\alpha|^2 + |\beta|^2 = 1$:*

$$M = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix} \tag{2.4}$$

*This is equivalent to the unit quaternion $q = (w, x, y, z)$ via the identification $\alpha = w + iz$, $\beta = ix + y$.*

**Assumption 2.2.1** (Right-handed coordinates)**.** *All coordinate systems are right-handed throughout. Rotation by angle $\theta$ about axis $\hat{\mathbf{e}}$ follows the right-hand rule.*

**Definition 2.2.5** (Rotation action)**.** *The rotation $M$ acts on a vector $V$ by conjugation:*

$$V' = M\, V\, M^\dagger \tag{2.5}$$

**Definition 2.2.6** (Elementary rotations)**.** *Rotation by angle $\alpha$ about coordinate axis $\hat{\mathbf{e}}_k$:*

$$M_k(\alpha) = \cos\frac{\alpha}{2}\,I + i\sin\frac{\alpha}{2}\,\sigma_k \tag{2.6}$$

*Explicitly:*

$$M_1(\alpha) = \begin{pmatrix} \cos\frac{\alpha}{2} & i\sin\frac{\alpha}{2} \\ i\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}, \quad M_2(\alpha) = \begin{pmatrix} \cos\frac{\alpha}{2} & \sin\frac{\alpha}{2} \\ -\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}, \quad M_3(\alpha) = \begin{pmatrix} e^{i\alpha/2} & 0 \\ 0 & e^{-i\alpha/2} \end{pmatrix} \tag{2.7}$$

*Note: $M_3(\alpha)$ is diagonal — rotation about $\hat{\mathbf{e}}_3$ is a phase multiplication.*

### Algebraic Properties

**Lemma 2.2.1** (Pauli algebra)**.** *The Pauli matrices satisfy:*

$$\sigma_j \sigma_k = \delta_{jk}\,I + i\,\varepsilon_{jkl}\,\sigma_l \tag{2.8}$$

*where $\delta_{jk}$ is the Kronecker delta and $\varepsilon_{jkl}$ is the Levi-Civita symbol. Equivalently:*
- *Anticommutator: $\{\sigma_j, \sigma_k\} = \sigma_j\sigma_k + \sigma_k\sigma_j = 2\delta_{jk}\,I$*
- *Commutator: $[\sigma_j, \sigma_k] = \sigma_j\sigma_k - \sigma_k\sigma_j = 2i\,\varepsilon_{jkl}\,\sigma_l$*
- *Trace: $\mathrm{tr}(\sigma_j\sigma_k) = 2\delta_{jk}$*

*Proof.* Direct computation from (2.1). For example, $\sigma_1\sigma_2 = \left(\begin{smallmatrix} 0 & 1 \\ 1 & 0 \end{smallmatrix}\right)\left(\begin{smallmatrix} 0 & -i \\ i & 0 \end{smallmatrix}\right) = \left(\begin{smallmatrix} i & 0 \\ 0 & -i \end{smallmatrix}\right) = i\sigma_3$, confirming $\sigma_1\sigma_2 = \delta_{12}I + i\varepsilon_{123}\sigma_3 = i\sigma_3$. The remaining products follow by cyclic permutation or by direct multiplication. ∎

**Theorem 2.2.1** (Cross product as commutator)**.** *For vectors $\mathbf{a}, \mathbf{b} \in \mathbb{R}^3$ encoded as traceless Hermitian matrices $A, B$:*

$$\mathbf{a} \times \mathbf{b} \;\longleftrightarrow\; \frac{1}{2i}[A, B] \tag{2.9}$$

*Proof.* Expand $A = \sum_j a_j\sigma_j$, $B = \sum_k b_k\sigma_k$. The commutator is:

$$[A, B] = \sum_{j,k} a_j b_k [\sigma_j, \sigma_k] = \sum_{j,k} a_j b_k \cdot 2i\varepsilon_{jkl}\sigma_l = 2i\sum_l \left(\sum_{j,k}\varepsilon_{jkl}\,a_j b_k\right)\sigma_l = 2i\,(\mathbf{a} \times \mathbf{b})\cdot\boldsymbol{\sigma}$$

Dividing by $2i$: $\frac{1}{2i}[A, B] = (\mathbf{a} \times \mathbf{b})\cdot\boldsymbol{\sigma}$, which is the traceless Hermitian matrix encoding $\mathbf{a} \times \mathbf{b}$. ∎

**Theorem 2.2.2** (Dot product from trace)**.** *For vectors encoded as traceless Hermitian matrices:*

$$\mathbf{a} \cdot \mathbf{b} = \frac{1}{2}\mathrm{tr}(AB) \tag{2.10}$$

*Proof.* $\mathrm{tr}(AB) = \sum_{j,k} a_j b_k\,\mathrm{tr}(\sigma_j\sigma_k) = \sum_{j,k} a_j b_k \cdot 2\delta_{jk} = 2\sum_j a_j b_j = 2(\mathbf{a}\cdot\mathbf{b})$. ∎

**Theorem 2.2.3** (Norm from determinant)**.** *The squared norm of a vector is:*

$$\|\mathbf{v}\|^2 = -\det(V) \tag{2.11}$$

*Proof.* $\det(V) = \det\left(\begin{smallmatrix} v_3 & v_1-iv_2 \\ v_1+iv_2 & -v_3 \end{smallmatrix}\right) = -v_3^2 - (v_1^2+v_2^2) = -(v_1^2+v_2^2+v_3^2)$. ∎

### Rotation Theorems

**Theorem 2.2.4** (Conjugation preserves vector structure)**.** *If $V$ is traceless Hermitian and $M \in SU(2)$, then $V' = MVM^\dagger$ is traceless Hermitian.*

*Proof.* Hermitian: $(V')^\dagger = (MVM^\dagger)^\dagger = M V^\dagger M^\dagger = MVM^\dagger = V'$ since $V^\dagger = V$. Traceless: $\mathrm{tr}(V') = \mathrm{tr}(MVM^\dagger) = \mathrm{tr}(V) = 0$ by the cyclic property of the trace. ∎

**Theorem 2.2.5** (Conjugation preserves norm)**.** *For $M \in SU(2)$ and traceless Hermitian $V$:*

$$\det(MVM^\dagger) = \det(V) \tag{2.12}$$

*Therefore $\|M\mathbf{v}M^\dagger\| = \|\mathbf{v}\|$: rotation preserves the vector norm.*

*Proof.* $\det(MVM^\dagger) = \det(M)\det(V)\det(M^\dagger) = 1 \cdot \det(V) \cdot 1 = \det(V)$, since $\det(M) = 1$ for $M \in SU(2)$. ∎

**Theorem 2.2.6** (Composition of rotations)**.** *If $M_{AB}$ rotates frame B to frame A, and $M_{BC}$ rotates frame C to frame B, then:*

$$M_{AC} = M_{AB} \cdot M_{BC} \tag{2.13}$$

*rotates frame C to frame A. The composition is ordinary $2 \times 2$ matrix multiplication.*

*Proof.* For any vector $V$ in frame C: $V_A = M_{AB}(M_{BC}\,V\,M_{BC}^\dagger)M_{AB}^\dagger = (M_{AB}M_{BC})\,V\,(M_{AB}M_{BC})^\dagger$, so $M_{AC} = M_{AB}M_{BC}$. Closure in $SU(2)$: $\det(M_{AC}) = \det(M_{AB})\det(M_{BC}) = 1$, and $(M_{AC})^\dagger M_{AC} = M_{BC}^\dagger M_{AB}^\dagger M_{AB} M_{BC} = I$. ∎

**Theorem 2.2.7** (Inverse rotation)**.** *For $M \in SU(2)$:*

$$M^{-1} = M^\dagger \tag{2.14}$$

*Proof.* Immediate from the definition $M^\dagger M = I$. ∎

**Lemma 2.2.2** (Double cover)**.** *$M$ and $-M$ produce the same rotation: $MVM^\dagger = (-M)V(-M)^\dagger$ for all traceless Hermitian $V$. Convention: choose the representative with $\mathrm{Re}(\alpha) \geq 0$.*

*Proof.* $(-M)V(-M)^\dagger = (-1)^2 MVM^\dagger = MVM^\dagger$. ∎

**Theorem 2.2.8** (Rotation from axis and angle)**.** *A rotation by angle $\theta$ about unit axis $\hat{\mathbf{e}} = (e_1, e_2, e_3)$ is:*

$$M(\theta, \hat{\mathbf{e}}) = \cos\frac{\theta}{2}\,I + i\sin\frac{\theta}{2}\,(e_1\sigma_1 + e_2\sigma_2 + e_3\sigma_3) \tag{2.15}$$

*Proof.* The elementary rotation (2.6) is the special case $\hat{\mathbf{e}} = \hat{\mathbf{e}}_k$. For a general axis, $M(\theta, \hat{\mathbf{e}}) = \exp(i\frac{\theta}{2}\hat{\mathbf{e}}\cdot\boldsymbol{\sigma})$. Since $(\hat{\mathbf{e}}\cdot\boldsymbol{\sigma})^2 = I$ (from Lemma 2.2.1 and $\|\hat{\mathbf{e}}\| = 1$), the matrix exponential reduces to the stated form via $e^{i\phi A} = \cos\phi\,I + i\sin\phi\,A$ when $A^2 = I$. ∎

**Theorem 2.2.9** (Extracting quaternion from rotation matrix — Shepperd's method)**.** *Given a $3 \times 3$ rotation matrix $R_{jk}$, the SU(2) element $M$ with parameters $(w, x, y, z)$ is extracted by:*

1. *Compute $T = R_{11} + R_{22} + R_{33}$ (trace) and the three diagonals $R_{11}, R_{22}, R_{33}$.*
2. *Select the largest of $\{1+T, 1+2R_{11}-T, 1+2R_{22}-T, 1+2R_{33}-T\}$.*
3. *Extract the corresponding component via a square root, then the remaining three via ratios.*

*The four cases (with $S$ denoting the square root):*
- *If $1+T$ is largest: $w = \frac{1}{2}\sqrt{1+T}$, $x = (R_{32}-R_{23})/(4w)$, $y = (R_{13}-R_{31})/(4w)$, $z = (R_{21}-R_{12})/(4w)$*
- *If $1+2R_{11}-T$ is largest: $x = \frac{1}{2}\sqrt{1+2R_{11}-T}$, then $w, y, z$ from ratios*
- *Similarly for $R_{22}$ and $R_{33}$ dominant*

*The branch selection ensures the denominator is bounded away from zero.* [P.2.1]

*Proof.* The relation between $R_{jk}$ and $(w,x,y,z)$ is given in §2.7, Theorem 2.7.1. Inverting that system: $1+T = 4w^2$, $1+2R_{11}-T = 4x^2$, etc. Choosing the largest ensures $|S| \geq \frac{1}{2}\sqrt{2}$, so the division is well-conditioned. ∎

### Error Properties of Rotation

A rotation acts on both the physical state and the error state (§2.6). The central question is: does a rotation amplify errors? The answer depends on whether the rotation itself is exact.

When the rotation matrix $M$ is known exactly — as when composing two previously computed rotations — conjugation is an isometry: it changes the direction of the error vector but not its magnitude.

**Theorem 2.2.10** (Isometry — exact rotation preserves error)**.** *For $M \in SU(2)$ with $\delta(M) = 0$:*

$$\delta(V') = M\,\delta(V)\,M^\dagger, \qquad \|\delta(\mathbf{v}')\| = \|\delta(\mathbf{v})\| \tag{2.16}$$

*Proof.* $V' = MVM^\dagger$ is linear in $V$ for fixed $M$, so $\delta(V') = M\,\delta(V)\,M^\dagger$ exactly. Norm preservation follows from Theorem 2.2.5: $\|\delta(\mathbf{v}')\|^2 = -\det(M\,\delta(V)\,M^\dagger) = -\det(\delta(V)) = \|\delta(\mathbf{v})\|^2$. ∎

In practice, $M$ is computed from orbital angles via trigonometric functions, and carries its own precision and measurement errors. These contribute an additional error term proportional to the magnitude of the vector being rotated.

**Theorem 2.2.11** (Inexact rotation error)**.** *When $M$ carries error $\delta(M)$:* [P.2.2]

$$\|\delta(V')\| \leq \|\delta(\mathbf{v})\| + 2\,\|\mathbf{v}\|\,\|\delta(M)\| \tag{2.17}$$

*Proof.* Expand $(M + \delta M)\,V\,(M + \delta M)^\dagger$. The zeroth-order term is $MVM^\dagger$. The first-order terms are $\delta M\,V\,M^\dagger + M\,V\,\delta M^\dagger$, bounded in norm by $2\|\mathbf{v}\|\,\|\delta(M)\|$ (two cross terms, each bounded by Ch 1, Theorem 1.3.2). The second-order term $\delta M\,V\,\delta M^\dagger$ is $O(\|\delta(M)\|^2)$. Combined with the isometry of Theorem 2.2.10 for the input error $\delta(V)$, the total bound is (2.17). ∎

To apply Theorem 2.2.11, one needs to know how large $\delta(M)$ is. Since $M_k(\alpha)$ is constructed from $\cos(\alpha/2)$ and $\sin(\alpha/2)$, its error is determined by the error in the input angle $\alpha$. The following theorem quantifies this relationship.

**Theorem 2.2.12** (Angle-to-rotation sensitivity)**.** *The elementary rotation $M_k(\alpha)$ satisfies:*

$$\left\|\frac{dM_k}{d\alpha}\right\| = \frac{1}{2}, \qquad \|\delta(M_k)\| \leq \frac{1}{2}\,|\delta(\alpha)| \tag{2.18}$$

*The rotation error is half the angle error, independent of $\alpha$.* [P.2.3]

*Proof.* Differentiate (2.6): $dM_k/d\alpha = -\frac{1}{2}\sin(\alpha/2)\,I + \frac{i}{2}\cos(\alpha/2)\,\sigma_k$. The eigenvalues are $\frac{1}{2}(-\sin(\alpha/2) \pm \cos(\alpha/2))$, each with modulus $\frac{1}{2}$. The bound follows from Ch 1, Theorem 1.4.1 with Lipschitz constant $\frac{1}{2}$. ∎

For the perifocal-to-inertial transform (§2.5), three elementary rotations are composed. Each contributes an angle error scaled by $\frac{1}{2}$, and by Theorem 2.2.10, each factor's contribution passes through subsequent rotations without amplification.

**Corollary 2.2.1** (Composed rotation error)**.** *For the 3-1-3 composition $M = M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u)$:*

$$\|\delta(M)\| \leq \frac{1}{2}\bigl(|\delta(\Omega)| + |\delta(i)| + |\delta(u)|\bigr) + 16\,\epsilon_{\mathrm{mach}} \tag{2.19}$$

*The rounding term $16\,\epsilon_{\mathrm{mach}}$ arises from two $SU(2)$ multiplications at $8\epsilon_{\mathrm{mach}}$ per entry each (for binary64: $16\epsilon_{\mathrm{mach}} \approx 2 \times 10^{-15}$).*

---

## §2.3 Incorporating Translation: The Dual Extension

### Dual Numbers

**Definition 2.3.1** (Dual numbers)**.** *A dual number is an expression $\hat{a} = a + \varepsilon\, a'$ where $a, a' \in \mathbb{C}$ and $\varepsilon$ is the dual unit satisfying $\varepsilon^2 = 0$, $\varepsilon \neq 0$. Arithmetic follows from this rule:*

$$(\hat{a})(\hat{b}) = (a + \varepsilon a')(b + \varepsilon b') = ab + \varepsilon(ab' + a'b) \tag{2.20}$$

*The dual part captures first-order perturbation: if $a(\epsilon) = a_0 + a_0'\epsilon$, then $\hat{a} = a_0 + \varepsilon a_0'$ encodes both the value and its derivative. This is algebraic automatic differentiation.*

### Dual SU(2)

**Definition 2.3.2** (Dual SU(2) matrix)**.** *A dual SU(2) matrix is:*

$$\hat{M} = M + \varepsilon\, D \tag{2.21}$$

*where $M \in SU(2)$ is the rotation and $D$ is a $2 \times 2$ complex matrix encoding the translation. The pair $(\hat{M})$ has 8 real parameters (4 from $M$, 4 from $D$).*

**Definition 2.3.3** (Unit dual SU(2) constraints)**.** *A unit dual SU(2) matrix satisfies:*

1. *$M^\dagger M = I$ (the non-dual part is unitary)*
2. *$M^\dagger D + D^\dagger M = 0$ (orthogonality: $D$ lies in the tangent space of $SU(2)$ at $M$)*

*These two constraints reduce the 8 real parameters to 6 degrees of freedom, matching the rigid body configuration: 3 for rotation + 3 for translation.*

**Definition 2.3.4** (Translation encoding)**.** *Given rotation $M \in SU(2)$ and translation vector $\mathbf{t} = (t_1, t_2, t_3)$, the dual part is:*

$$D = \frac{1}{2}\,T \cdot M \tag{2.22}$$

*where $T = t_1\sigma_1 + t_2\sigma_2 + t_3\sigma_3$ is the translation encoded as a traceless Hermitian matrix (Definition 2.2.3). The factor of $\frac{1}{2}$ is a convention that simplifies the composition rule.*

**Theorem 2.3.1** (Translation extraction)**.** *Given a unit dual SU(2) matrix $\hat{M} = M + \varepsilon D$, the translation vector is recovered by:*

$$T = 2\,D\,M^\dagger \tag{2.23}$$

*and the rotation is simply $M$ (the non-dual part of $\hat{M}$).*

*Proof.* From (2.22): $D = \frac{1}{2}TM$. Multiply on the right by $M^\dagger$: $DM^\dagger = \frac{1}{2}TMM^\dagger = \frac{1}{2}T$ (since $MM^\dagger = I$). Therefore $T = 2DM^\dagger$.

It remains to verify that $T = 2DM^\dagger$ is traceless Hermitian. By construction, $T$ was defined as a traceless Hermitian matrix (Definition 2.3.4), and $D = \frac{1}{2}TM$ is derived from it. The extraction formula $T = 2DM^\dagger$ recovers the original $T$ by the invertibility of right-multiplication by $M$ (which is invertible since $\det M = 1$). Therefore $2DM^\dagger$ equals the traceless Hermitian matrix $T$ that generated $D$. ∎

### Composition and Inverse

**Theorem 2.3.2** (Dual SU(2) composition)**.** *The product of two dual SU(2) matrices is:*

$$\hat{M}_1 \hat{M}_2 = (M_1 + \varepsilon D_1)(M_2 + \varepsilon D_2) = M_1 M_2 + \varepsilon(M_1 D_2 + D_1 M_2) \tag{2.24}$$

*since $\varepsilon^2 = 0$. The composed rotation is $M_1 M_2$ (Theorem 2.2.6). The composed translation, extracted via Theorem 2.3.1, is:*

$$T_{12} = 2(M_1 D_2 + D_1 M_2)(M_1 M_2)^\dagger = 2M_1 D_2 M_2^\dagger M_1^\dagger + 2D_1 M_1^\dagger$$

$$= M_1\,T_2\,M_1^\dagger + T_1 \tag{2.25}$$

*where $T_k = 2D_k M_k^\dagger$. In words: rotate $\mathbf{t}_2$ by $M_1$ (via conjugation), then add $\mathbf{t}_1$. This is the natural composition of two rigid motions — apply the second displacement, then the first.*

*Proof.* The product follows from dual number arithmetic (2.20). For the translation extraction: $(M_1 D_2 + D_1 M_2)(M_1 M_2)^\dagger = M_1 D_2 M_2^\dagger M_1^\dagger + D_1 M_2 M_2^\dagger M_1^\dagger = M_1 D_2 M_2^\dagger M_1^\dagger + D_1 M_1^\dagger$. Using $T_k = 2D_k M_k^\dagger$: $D_k M_k^\dagger = \frac{1}{2}T_k$, and $M_1(\frac{1}{2}T_2)M_1^\dagger = \frac{1}{2}M_1 T_2 M_1^\dagger$. Multiplying by 2 gives (2.25). ∎

**Theorem 2.3.3** (Dual SU(2) inverse)**.** *The inverse of a unit dual SU(2) matrix is:*

$$\hat{M}^{-1} = M^\dagger - \varepsilon\,D^\dagger \tag{2.26}$$

*Proof.* Verify $\hat{M}\,\hat{M}^{-1} = \hat{I}$:

$$(M + \varepsilon D)(M^\dagger - \varepsilon D^\dagger) = MM^\dagger + \varepsilon(- M D^\dagger + D M^\dagger) = I + \varepsilon(-MD^\dagger + DM^\dagger)$$

The dual part is $-MD^\dagger + DM^\dagger$. Substituting $D = \frac{1}{2}TM$ (Definition 2.3.4): $DM^\dagger = \frac{1}{2}TMM^\dagger = \frac{1}{2}T$ and $MD^\dagger = M(\frac{1}{2}TM)^\dagger = \frac{1}{2}MM^\dagger T^\dagger = \frac{1}{2}T$ (using $MM^\dagger = I$ and $T^\dagger = T$). Therefore $-MD^\dagger + DM^\dagger = -\frac{1}{2}T + \frac{1}{2}T = 0$. ∎

**Theorem 2.3.4** (Equivalence to $4 \times 4$ homogeneous matrix)**.** *The dual SU(2) composition rule (2.25) produces the same rotation and translation as the $4 \times 4$ homogeneous matrix product. Specifically, if $H_k = \left(\begin{smallmatrix} R_k & \mathbf{t}_k \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$ where $R_k = R(M_k)$ is the $3 \times 3$ rotation matrix from $M_k$ (derived in §2.7), then:*

$$H_1 H_2 = \begin{pmatrix} R_1 R_2 & R_1\mathbf{t}_2 + \mathbf{t}_1 \\ \mathbf{0}^T & 1 \end{pmatrix}$$

*The composed rotation $R_1 R_2$ corresponds to $M_1 M_2$ (Theorem 2.2.6), and the composed translation $R_1\mathbf{t}_2 + \mathbf{t}_1$ corresponds to $M_1 T_2 M_1^\dagger + T_1$ (equation 2.25), which is the same vector expressed in the Hermitian matrix encoding.*

*Proof.* The $4 \times 4$ product is standard block multiplication. The equivalence between $M_1 T_2 M_1^\dagger$ (SU(2) conjugation) and $R_1\mathbf{t}_2$ (matrix-vector product) is established in §2.7, Theorem 2.7.1. ∎

### Linearity of the Dual SU(2) Action

**Theorem 2.3.5** (The dual action is linear in the operand)**.** *For fixed $\hat{M}$, the map $V \mapsto MVM^\dagger + T$ (rotate and translate a vector) is an affine function of $V$. The linear part is the conjugation $V \mapsto MVM^\dagger$, and the constant part is the translation $T$.*

*In dual form: the action $\hat{V} \mapsto \hat{M}\hat{V}\hat{M}^\dagger$ on a dual-encoded point is linear in $\hat{V}$ for fixed $\hat{M}$.*

*Proof.* Conjugation $V \mapsto MVM^\dagger$ is linear in $V$ (Theorem 2.2.10 used this). Addition of the constant $T$ makes the full action affine. For the error state: $\delta(MVM^\dagger + T) = M\,\delta(V)\,M^\dagger$ — the constant $T$ drops out (Ch 1, Theorem 1.3.1: the error of a sum is the sum of the errors, and a constant contributes zero error). ∎

**Remark.** Theorem 2.3.5 is the reason the parallel-row error architecture works for rotation + translation: the translation is absorbed into the algebra, and the error row transforms by the same conjugation as the physical row. The constant translation contributes zero error (assuming the translation itself is exact). When the translation carries its own error $\delta(T)$, it adds to the position error directly: $\delta(V') = M\,\delta(V)\,M^\dagger + \delta(T)$.

### Position Extraction Error

**Theorem 2.3.6** (Error in position extraction)**.** *Extracting position from the dual part via $T = 2DM^\dagger$ (Theorem 2.3.1) is a bilinear operation in $D$ and $M^\dagger$. The error bound is:*

$$\|\delta(\mathbf{r})\| \leq \|\delta(D)\| \cdot \|M^\dagger\| + \|D\| \cdot \|\delta(M^\dagger)\| = \delta_D + \frac{1}{2}\|\mathbf{r}\|\,\delta_M \tag{2.27}$$

*where $\delta_D = 2\|\delta(D)\|$ is the dual-part error (scaled by the factor of 2 in the extraction), $\delta_M = \|\delta(M)\|$ is the rotation error, and $\|D\| = \frac{1}{2}\|T\| = \frac{1}{2}\|\mathbf{r}\|$ from the encoding (2.22). The $\|M^\dagger\| = 1$ factor (unitarity) means the dual-part error passes through unscaled.* [P.2.4]

*Proof.* $T = 2DM^\dagger$ is a product of two matrices. By Ch 1, Theorem 1.3.2 (multiplication bound, valid for complex matrix entries): $\|\delta(T)\| \leq 2(\|D\|\,\|\delta(M^\dagger)\| + \|M^\dagger\|\,\|\delta(D)\|) = 2(\frac{1}{2}\|\mathbf{r}\|\,\delta_M + 1 \cdot \|\delta(D)\|)$. Since $T$ encodes a position vector, $\|\delta(\mathbf{r})\| = \sqrt{-\det(\delta T)}$ (Theorem 2.2.3), and the bound follows from the submultiplicativity of the determinant. ∎

**Remark.** The position error has two sources: (1) the dual-part error $\delta_D$, which represents direct uncertainty in the position encoding, and (2) the rotation error $\delta_M$ coupled through the position magnitude $\|\mathbf{r}\|$. For a LEO satellite at $\|\mathbf{r}\| \approx 7000$ km with TLE-level angle errors $\delta_M \sim 5 \times 10^{-9}$ (from Theorem 2.2.12 with $\delta(\alpha) \sim 10^{-8}$): the rotation-coupled position error is $\sim 7000 \times 5 \times 10^{-9} \approx 3.5 \times 10^{-5}$ km $\approx 0.035$ m. This is negligible compared to the TLE position accuracy ($\sim 1$ km).

### Elementary Operations: Worked Examples

The following theorems demonstrate every elementary operation — the building blocks from which all orbit propagation transforms are composed. Each is verified by explicit matrix computation. Together they show that the $SU(2)$ framework reproduces the full repertoire of the classical $4 \times 4$ homogeneous matrix.

**Theorem 2.3.7** (Elementary translation along $\hat{\mathbf{e}}_k$)**.** *A translation by $\lambda$ along coordinate axis $\hat{\mathbf{e}}_k$ (no rotation) is the dual $SU(2)$ element:*

$$\hat{M}_{\mathrm{trans}}(\lambda, k) = I + \varepsilon\,\frac{\lambda}{2}\,\sigma_k \tag{2.57}$$

*The rotation part is $M = I$ (identity — no rotation). The dual part is $D = \frac{\lambda}{2}\sigma_k$, encoding the translation $\mathbf{t} = \lambda\hat{\mathbf{e}}_k$ via $D = \frac{1}{2}T \cdot I = \frac{1}{2}T$.*

*Action on a position vector $R_{\mathrm{pos}}$:*

$$R_{\mathrm{pos}}' = I \cdot R_{\mathrm{pos}} \cdot I + 2D \cdot I = R_{\mathrm{pos}} + \lambda\sigma_k$$

*which adds $\lambda$ to the $k$-th component of $\mathbf{r}$. Explicitly for each axis:*

*Along $\hat{\mathbf{e}}_1$:* $\;D = \frac{\lambda}{2}\sigma_1 = \frac{\lambda}{2}\left(\begin{smallmatrix} 0 & 1 \\ 1 & 0 \end{smallmatrix}\right)$. *The translated position has $r_1' = r_1 + \lambda$, $r_2' = r_2$, $r_3' = r_3$.*

*Along $\hat{\mathbf{e}}_2$:* $\;D = \frac{\lambda}{2}\sigma_2 = \frac{\lambda}{2}\left(\begin{smallmatrix} 0 & -i \\ i & 0 \end{smallmatrix}\right)$. *The translated position has $r_2' = r_2 + \lambda$, with $r_1, r_3$ unchanged.*

*Along $\hat{\mathbf{e}}_3$:* $\;D = \frac{\lambda}{2}\sigma_3 = \frac{\lambda}{2}\left(\begin{smallmatrix} 1 & 0 \\ 0 & -1 \end{smallmatrix}\right)$. *The translated position has $r_3' = r_3 + \lambda$, with $r_1, r_2$ unchanged.*

*Proof.* From Definition 2.3.4: $D = \frac{1}{2}TM = \frac{1}{2}(\lambda\sigma_k)(I) = \frac{\lambda}{2}\sigma_k$. Extraction (Theorem 2.3.1): $T = 2DM^\dagger = 2 \cdot \frac{\lambda}{2}\sigma_k \cdot I = \lambda\sigma_k$. The vector encoded by $\lambda\sigma_k$ is $\lambda\hat{\mathbf{e}}_k$ (by Definition 2.2.3). The action $R_{\mathrm{pos}}' = MR_{\mathrm{pos}}M^\dagger + T = R_{\mathrm{pos}} + \lambda\sigma_k$ adds $\lambda$ to the $k$-th component. ∎

**Theorem 2.3.8** (General translation)**.** *A translation by $\mathbf{t} = (\lambda_1, \lambda_2, \lambda_3)$ is:*

$$\hat{M}_{\mathrm{trans}}(\mathbf{t}) = I + \varepsilon\,\frac{1}{2}(\lambda_1\sigma_1 + \lambda_2\sigma_2 + \lambda_3\sigma_3) = I + \varepsilon\,\frac{1}{2}\begin{pmatrix} \lambda_3 & \lambda_1 - i\lambda_2 \\ \lambda_1 + i\lambda_2 & -\lambda_3 \end{pmatrix} \tag{2.58}$$

*This is the superposition of the three elementary translations (Theorem 2.3.7).*

*Proof.* $T = \lambda_1\sigma_1 + \lambda_2\sigma_2 + \lambda_3\sigma_3$ encodes $\mathbf{t}$ (Definition 2.2.3). $D = \frac{1}{2}TI = \frac{1}{2}T$. The action adds $T$ to $R_{\mathrm{pos}}$, translating each component by $\lambda_k$. ∎

**Theorem 2.3.9** (Elementary rotation about $\hat{\mathbf{e}}_1$)**.** *The passive rotation $V' = M_1(\alpha)\,V\,M_1(\alpha)^\dagger$ transforms a vector $\mathbf{v} = (v_1, v_2, v_3)$ to:*

$$v_1' = v_1, \qquad v_2' = v_2\cos\alpha + v_3\sin\alpha, \qquad v_3' = -v_2\sin\alpha + v_3\cos\alpha \tag{2.59}$$

*This is the standard passive rotation matrix $R_1(\alpha)$.*

*Proof.* Write $c = \cos\frac{\alpha}{2}$, $s = \sin\frac{\alpha}{2}$. The product $V' = M_1 V M_1^\dagger$ yields:

$$(V')_{11} = c^2 v_3 + ics(v_1 + iv_2) - ics(v_1 - iv_2) - s^2 v_3 = (c^2 - s^2)v_3 - 2csv_2 = v_3\cos\alpha - v_2\sin\alpha$$

$$(V')_{12} = (c^2 + s^2)v_1 - i(c^2 - s^2)v_2 - 2icsv_3 = v_1 - i(v_2\cos\alpha + v_3\sin\alpha)$$

By Definition 2.2.3, $(V')_{11} = v_3'$ and $(V')_{12} = v_1' - iv_2'$. Therefore $v_3' = v_3\cos\alpha - v_2\sin\alpha$, $v_1' = v_1$, and $v_2' = v_2\cos\alpha + v_3\sin\alpha$, confirming (2.59). ∎

**Theorem 2.3.10** (Elementary rotation about $\hat{\mathbf{e}}_3$)**.** *The passive rotation $V' = M_3(\alpha)\,V\,M_3(\alpha)^\dagger$ transforms a vector to:*

$$v_1' = v_1\cos\alpha + v_2\sin\alpha, \qquad v_2' = -v_1\sin\alpha + v_2\cos\alpha, \qquad v_3' = v_3 \tag{2.60}$$

*This is the standard passive rotation matrix $R_3(\alpha)$.*

*Proof.* $M_3(\alpha) = \left(\begin{smallmatrix} e^{i\alpha/2} & 0 \\ 0 & e^{-i\alpha/2} \end{smallmatrix}\right)$ is diagonal, so $M_3 V M_3^\dagger$ multiplies the off-diagonal entry $(v_1 - iv_2)$ by $e^{i\alpha/2} \cdot e^{i\alpha/2} = e^{i\alpha}$ and the diagonal entries are unchanged. Therefore $v_1' - iv_2' = (v_1 - iv_2)e^{i\alpha} = (v_1 - iv_2)(\cos\alpha + i\sin\alpha) = (v_1\cos\alpha + v_2\sin\alpha) - i(-v_1\sin\alpha + v_2\cos\alpha)$, and $v_3' = v_3$, confirming (2.60). ∎

### Sign Convention

**Remark** (Passive vs. active rotation)**.** The conjugation $V' = M_k(\alpha)\,V\,M_k(\alpha)^\dagger$ implements a **passive** rotation: it expresses the vector $\mathbf{v}$ in a frame rotated by $+\alpha$ about $\hat{\mathbf{e}}_k$. Equivalently, the vector appears to rotate by $-\alpha$ relative to the new frame. To rotate a vector **actively** by $+\alpha$ (keeping the frame fixed), negate the angle: $V' = M_k(-\alpha)\,V\,M_k(-\alpha)^\dagger$.

This is verified by a concrete test: conjugating $\hat{\mathbf{e}}_1 = (1,0,0)$ by $M_3(+\pi/2)$ gives $(0, -1, 0)$ — the $x$-axis of the old frame points along $-y$ in a frame rotated $+\pi/2$ about $z$, which is the correct passive result. Negating the angle gives $(0, +1, 0)$, which is the correct active result.

All theorems in §§2.2–2.5 use the passive convention (frame rotation), which is the standard convention in orbit propagation. The negative angles in Theorem 2.5.1 ($M_3(-\Omega)$, $M_1(-i)$, $M_3(-u)$) undo the orbital rotations, transforming from the perifocal frame to the inertial frame.

**Theorem 2.3.11** (Elementary rotation about $\hat{\mathbf{e}}_2$)**.** *The passive rotation $V' = M_2(\alpha)\,V\,M_2(\alpha)^\dagger$ transforms a vector to:*

$$v_1' = v_1\cos\alpha - v_3\sin\alpha, \qquad v_2' = v_2, \qquad v_3' = v_1\sin\alpha + v_3\cos\alpha \tag{2.61}$$

*This is the standard passive rotation matrix $R_2(\alpha)$.*

*Proof.* Write $c = \cos\frac{\alpha}{2}$, $s = \sin\frac{\alpha}{2}$. The product $V' = M_2 V M_2^\dagger$ yields:

$$(V')_{11} = (c^2 - s^2)v_3 + 2csv_1 = v_3\cos\alpha + v_1\sin\alpha$$

$$(V')_{12} = (c^2 - s^2)v_1 - i(c^2 + s^2)v_2 - 2csv_3 = (v_1\cos\alpha - v_3\sin\alpha) - iv_2$$

By Definition 2.2.3, $(V')_{11} = v_3'$ and $(V')_{12} = v_1' - iv_2'$. Therefore $v_3' = v_1\sin\alpha + v_3\cos\alpha$, $v_1' = v_1\cos\alpha - v_3\sin\alpha$, and $v_2' = v_2$, confirming (2.61). ∎

**Theorem 2.3.12** (Composition of translation and rotation)**.** *A rotation by $\alpha$ about $\hat{\mathbf{e}}_3$ followed by a translation by $\lambda$ along $\hat{\mathbf{e}}_1$ is the dual $SU(2)$ element:*

$$\hat{M} = \hat{M}_{\mathrm{trans}}(\lambda, 1) \cdot \hat{M}_{\mathrm{rot}}(\alpha, 3) = (I + \varepsilon\frac{\lambda}{2}\sigma_1) \cdot (M_3(\alpha) + \varepsilon \cdot 0) \tag{2.62}$$

$$= M_3(\alpha) + \varepsilon\frac{\lambda}{2}\sigma_1 M_3(\alpha)$$

*The rotation part is $M_3(\alpha)$ (from the rotation factor). The dual part is $D = \frac{\lambda}{2}\sigma_1 M_3(\alpha)$. Extracting the translation: $T = 2DM_3^\dagger = \lambda\sigma_1$, so $\mathbf{t} = \lambda\hat{\mathbf{e}}_1$.*

*The action on a position: first rotate the vector by $\alpha$ about $\hat{\mathbf{e}}_3$ (passive), then translate by $\lambda$ along $\hat{\mathbf{e}}_1$. This is:*

$$\mathbf{r}' = R_3(\alpha)\mathbf{r} + \lambda\hat{\mathbf{e}}_1$$

*which matches the $4 \times 4$ homogeneous result $H = \left(\begin{smallmatrix} R_3(\alpha) & \lambda\hat{\mathbf{e}}_1 \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$.*

*Proof.* Dual quaternion product (Theorem 2.3.2): non-dual part $= I \cdot M_3(\alpha) = M_3(\alpha)$; dual part $= I \cdot 0 + \frac{\lambda}{2}\sigma_1 \cdot M_3(\alpha) = \frac{\lambda}{2}\sigma_1 M_3(\alpha)$. Extraction: $T = 2 \cdot \frac{\lambda}{2}\sigma_1 M_3 \cdot M_3^\dagger = \lambda\sigma_1$. The composed action is (equation 2.29): $R_{\mathrm{pos}}' = M_3 R_{\mathrm{pos}} M_3^\dagger + \lambda\sigma_1$. ∎

**Theorem 2.3.13** (Composition of two translations)**.** *Two translations compose by vector addition:*

$$\hat{M}_{\mathrm{trans}}(\mathbf{t}_1) \cdot \hat{M}_{\mathrm{trans}}(\mathbf{t}_2) = \hat{M}_{\mathrm{trans}}(\mathbf{t}_1 + \mathbf{t}_2) \tag{2.63}$$

*Proof.* Both have $M = I$. By Theorem 2.3.2: non-dual part $= I \cdot I = I$; dual part $= I \cdot D_2 + D_1 \cdot I = D_1 + D_2 = \frac{1}{2}(T_1 + T_2)$. The combined translation is $T_1 + T_2$, encoding $\mathbf{t}_1 + \mathbf{t}_2$. ∎

**Theorem 2.3.14** (Composition of two rotations about different axes)**.** *A rotation by $\alpha$ about $\hat{\mathbf{e}}_1$ followed by a rotation by $\beta$ about $\hat{\mathbf{e}}_3$ (both passive) is:*

$$M = M_3(\beta) \cdot M_1(\alpha) = \begin{pmatrix} e^{i\beta/2} & 0 \\ 0 & e^{-i\beta/2} \end{pmatrix}\begin{pmatrix} \cos\frac{\alpha}{2} & i\sin\frac{\alpha}{2} \\ i\sin\frac{\alpha}{2} & \cos\frac{\alpha}{2} \end{pmatrix}$$

$$= \begin{pmatrix} e^{i\beta/2}\cos\frac{\alpha}{2} & ie^{i\beta/2}\sin\frac{\alpha}{2} \\ ie^{-i\beta/2}\sin\frac{\alpha}{2} & e^{-i\beta/2}\cos\frac{\alpha}{2} \end{pmatrix} \tag{2.64}$$

*This is a single $SU(2)$ element (verified: $|\alpha'|^2 + |\beta'|^2 = \cos^2\frac{\alpha}{2} + \sin^2\frac{\alpha}{2} = 1$). It encodes a combined rotation that is NOT about any coordinate axis — it is a general rotation, achieved by composing two elementary ones.*

*Proof.* Direct $2 \times 2$ matrix multiplication. The result is $SU(2)$: $\det = e^{i\beta/2}e^{-i\beta/2}(\cos^2\frac{\alpha}{2} + \sin^2\frac{\alpha}{2}) = 1$. ∎

---

## §2.4 The Complete Rigid Body State

### State Definition

**Definition 2.4.1** (Rigid body state)**.** *The complete state of a rigid body consists of four quantities, each encoded as a $2 \times 2$ matrix:*

| Component | Matrix type | Real parameters | Physical meaning |
|-----------|-----------|-----------------|-----------------|
| Orientation $M$ | $SU(2)$ (unitary, $\det = 1$) | 4 (with 1 constraint → 3 DOF) | Which way the body is pointing |
| Position $R_{\mathrm{pos}}$ | Traceless Hermitian | 3 | Where the body is |
| Velocity $V_{\mathrm{vel}}$ | Traceless Hermitian | 3 | How fast it's translating |
| Angular velocity $\Omega$ | Traceless Hermitian | 3 | How fast it's rotating |

*Total: 13 real parameters, 12 degrees of freedom. The configuration (position + orientation) is encoded as a dual SU(2) matrix $\hat{M} = M + \varepsilon D$ where $D = \frac{1}{2}R_{\mathrm{pos}} M$ (Definition 2.3.4).*

**Assumption 2.4.1** (Instantaneous angular velocity)**.** *The angular velocity $\Omega_T$ of a rotating reference frame is treated as constant during a single transformation step. For time-varying $\Omega_T(t)$, additional angular acceleration terms appear. For Earth rotation, this assumption is valid: the rotation rate varies by $< 10^{-8}$ rad/s over any practical propagation step.*

### The State Transform

**Theorem 2.4.1** (Rigid body state transform)**.** *A frame transform with rotation $M_T \in SU(2)$, translation $T_T$ (traceless Hermitian), and frame angular velocity $\Omega_T$ (traceless Hermitian) acts on the state as:*

$$M' = M_T \cdot M \tag{2.28}$$

$$R_{\mathrm{pos}}' = M_T\,R_{\mathrm{pos}}\,M_T^\dagger + T_T \tag{2.29}$$

$$V_{\mathrm{vel}}' = M_T\,V_{\mathrm{vel}}\,M_T^\dagger + \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}'] \tag{2.30}$$

$$\Omega' = M_T\,\Omega\,M_T^\dagger + \Omega_T \tag{2.31}$$

*Equation (2.28): orientation composes by SU(2) multiplication (Theorem 2.2.6).*
*Equation (2.29): position rotates and translates (Theorem 2.3.5).*
*Equation (2.30): velocity rotates and acquires the transport term $\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$, which is $\boldsymbol{\omega}_T \times \mathbf{r}'$ (Theorem 2.2.1).*
*Equation (2.31): angular velocity rotates and adds the frame rotation rate.*

*Proof.* Equations (2.28), (2.29), and (2.31) follow directly from the composition rules of §§2.2–2.3. Equation (2.30) is the transport theorem, derived in Theorem 2.4.2 below. ∎

### The Transport Theorem

**Theorem 2.4.2** (Transport theorem from time derivative)**.** *If $M_T(t)$ is a time-varying rotation (the frame rotates with angular velocity $\boldsymbol{\omega}_T$), the velocity of a position vector that is stationary in the rotating frame is:*

$$\frac{d}{dt}\left(M_T R_{\mathrm{pos}} M_T^\dagger\right) = \frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}'] \tag{2.32}$$

*where $R_{\mathrm{pos}}' = M_T R_{\mathrm{pos}} M_T^\dagger$ is the rotated position and $\Omega_T' = M_T \Omega_T M_T^\dagger$ is the angular velocity in the target frame. This is the transport velocity $\boldsymbol{\omega}_T \times \mathbf{r}'$ (Theorem 2.2.1).*

*Proof.* For a unitary matrix rotating at angular velocity $\Omega_T$ (expressed in the body frame), the time derivative satisfies $\dot{M}_T = \frac{i}{2}\Omega_T' M_T$ where $\Omega_T' = M_T \Omega_T M_T^\dagger$. Correspondingly, $\dot{M}_T^\dagger = -M_T^\dagger \frac{i}{2}\Omega_T'$. By the product rule:

$$\frac{d}{dt}(M_T R_{\mathrm{pos}} M_T^\dagger) = \dot{M}_T R_{\mathrm{pos}} M_T^\dagger + M_T R_{\mathrm{pos}} \dot{M}_T^\dagger$$

$$= \frac{i}{2}\Omega_T' M_T R_{\mathrm{pos}} M_T^\dagger - M_T R_{\mathrm{pos}} M_T^\dagger \frac{i}{2}\Omega_T' = \frac{i}{2}(\Omega_T' R_{\mathrm{pos}}' - R_{\mathrm{pos}}' \Omega_T') = \frac{i}{2}[\Omega_T', R_{\mathrm{pos}}']$$

Since $\frac{i}{2} = -\frac{1}{2i}$ (multiply numerator and denominator by $i$: $\frac{i \cdot i}{2i} = \frac{-1}{2i}$), and using the antisymmetry of the commutator $[A,B] = -[B,A]$:

$$\frac{d}{dt}(M_T R_{\mathrm{pos}} M_T^\dagger) = -\frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}'] = \frac{1}{2i}[R_{\mathrm{pos}}', \Omega_T']$$

Under the right-hand convention, $\boldsymbol{\omega} \times \mathbf{r} \leftrightarrow \frac{1}{2i}[\Omega, R_{\mathrm{pos}}]$ (Theorem 2.2.1). The derivative gives $\frac{1}{2i}[R_{\mathrm{pos}}', \Omega_T'] = -\frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}']$, which corresponds to $\mathbf{r}' \times \boldsymbol{\omega}_T = -\boldsymbol{\omega}_T \times \mathbf{r}'$. Therefore the transport velocity in the right-hand convention is:

$$V_{\mathrm{transport}} = \frac{1}{2i}[\Omega_T', R_{\mathrm{pos}}'] \;\longleftrightarrow\; \boldsymbol{\omega}_T \times \mathbf{r}'$$

which is equation (2.30)'s transport term. ∎

**Remark.** The transport theorem is the product rule of matrix differentiation — it is not an independent physical postulate. The velocity coupling $\boldsymbol{\omega} \times \mathbf{r}$ emerges algebraically from differentiating the SU(2) conjugation $M(t)VM(t)^\dagger$.

### Properties of the Transport Term

**Theorem 2.4.3** (Transport term is linear in position)**.** *For fixed $\Omega_T$, the map $R_{\mathrm{pos}} \mapsto \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ is linear:*

$$\frac{1}{2i}[\Omega_T, \alpha R_1 + \beta R_2] = \alpha\,\frac{1}{2i}[\Omega_T, R_1] + \beta\,\frac{1}{2i}[\Omega_T, R_2] \tag{2.33}$$

*Proof.* The commutator is bilinear: $[A, \alpha B_1 + \beta B_2] = \alpha[A, B_1] + \beta[A, B_2]$. ∎

**Corollary 2.4.1** (Transport term error)**.** *The transport velocity error from position error is:*

$$\delta(V_{\mathrm{transport}}) = \frac{1}{2i}[\Omega_T, \delta(R_{\mathrm{pos}})] \tag{2.34}$$

*This is the same commutator operation applied to the error — no separate formula needed (same linear map, different input). The error magnitude is:*

$$\|\delta(\mathbf{v}_{\mathrm{transport}})\| = |\omega_T| \cdot \|\delta(\mathbf{r})\| \tag{2.35}$$

*where $|\omega_T| = \sqrt{-\det(\Omega_T)}$ is the angular velocity magnitude.*

*Proof.* Linearity (Theorem 2.4.3) gives (2.34). For the magnitude: $\|\boldsymbol{\omega} \times \delta\mathbf{r}\| = |\omega|\,\|\delta\mathbf{r}\|\,|\sin\phi|$ where $\phi$ is the angle between $\boldsymbol{\omega}$ and $\delta\mathbf{r}$. The worst case $|\sin\phi| = 1$ gives (2.35). ∎

**Remark.** For Earth rotation: $|\omega_E| \approx 7.292 \times 10^{-5}$ rad/s. A position error of $\delta r = 1$ km produces a transport velocity error of $|\omega_E| \times 1 \approx 7.3 \times 10^{-5}$ km/s $= 0.073$ m/s. This is the primary mechanism by which position error couples into velocity error through the orbit propagation pipeline.

### Point-Mass Reduction

**Corollary 2.4.2** (Point-mass reduction)**.** *For a point-mass satellite with no body attitude tracking (the standard mode of Ch 36):*
- *$M = M_T$ (orientation is purely the frame rotation, not body attitude)*
- *$\Omega = 0$ (no body angular velocity)*
- *The 13-parameter state reduces to 6 active quantities: $R_{\mathrm{pos}}$ (3) and $V_{\mathrm{vel}}$ (3), with $M$ carrying the frame rotation and $\Omega = 0$*
- *The transform equations (2.28)–(2.31) reduce to:*
  - *Position: $R_{\mathrm{pos}}' = M_T R_{\mathrm{pos}} M_T^\dagger + T_T$ (rotate and translate)*
  - *Velocity: $V_{\mathrm{vel}}' = M_T V_{\mathrm{vel}} M_T^\dagger + \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$ (rotate with transport)*
  - *Angular velocity: $\Omega' = \Omega_T$ (frame rate only)*

*No special-casing is required — setting $\Omega = 0$ in the general equations produces the reduced form naturally.*

*Proof.* Substitute $\Omega = 0$ into (2.28)–(2.31). Equation (2.31) gives $\Omega' = M_T \cdot 0 \cdot M_T^\dagger + \Omega_T = \Omega_T$. The remaining equations are unchanged. ∎

### Near-Identity Perturbations

**Theorem 2.4.4** (Near-identity composition is additive)**.** *A small perturbation to the rotation is $M_{\mathrm{pert}} = I + \frac{i}{2}\varepsilon_k\sigma_k + O(\varepsilon^2)$ where $|\varepsilon_k| \ll 1$. Composing two near-identity perturbations:*

$$(I + \frac{i}{2}A)(I + \frac{i}{2}B) = I + \frac{i}{2}(A + B) + O(\varepsilon^2) \tag{2.36}$$

*where $A = \varepsilon_k^{(1)}\sigma_k$ and $B = \varepsilon_k^{(2)}\sigma_k$. To first order, perturbations compose by addition.*

*Proof.* Expand the product: $I + \frac{i}{2}A + \frac{i}{2}B + \frac{i^2}{4}AB = I + \frac{i}{2}(A+B) - \frac{1}{4}AB$. The last term is $O(\varepsilon^2)$. ∎

**Remark.** Theorem 2.4.4 is used in Ch 18–19 (short-period and long-period corrections). Each correction is a near-identity perturbation to the orbital state. Composing multiple corrections is approximately additive — the corrections sum rather than multiply. The second-order error from this approximation is $O(\varepsilon_1 \varepsilon_2)$, where $\varepsilon_k$ are the perturbation magnitudes ($\sim J_2 \sim 10^{-3}$ for the dominant term). The cross term is $\sim 10^{-6}$, well below the model accuracy $\delta_a$.

## §2.5 The Orbit Propagation Coordinate Transform Pipeline

### Frame Definitions

**Definition 2.5.1** (Perifocal frame)**.** *The orbital-plane frame with $\hat{\mathbf{x}}$ toward perigee, $\hat{\mathbf{z}}$ along the angular momentum vector, and $\hat{\mathbf{y}}$ completing the right-handed triad. Position in this frame is $(r\cos\nu, r\sin\nu, 0)$ where $r$ is the orbital radius and $\nu$ is the true anomaly.*

**Definition 2.5.2** (TEME frame)**.** *True Equator, Mean Equinox — the reference frame used by the propagator (Ch 30). The $\hat{\mathbf{z}}$ axis is aligned with the Earth's instantaneous rotation axis, and $\hat{\mathbf{x}}$ points toward the mean vernal equinox.* [A.2.1]

**Definition 2.5.3** (PEF frame)**.** *Pseudo Earth-Fixed — TEME rotated by the Greenwich Mean Sidereal Time angle $\theta_{\mathrm{GMST}}$ (Ch 29) about the polar axis. Approximately Earth-fixed (neglecting polar motion).*

### Transform 1: Perifocal → TEME

**Theorem 2.5.1** (Perifocal → TEME rotation)**.** *The rotation from the perifocal frame to TEME is the 3-1-3 Euler composition:*

$$M_{\mathrm{PF \to TEME}} = M_3(-\Omega_{\mathrm{node}}) \cdot M_1(-i) \cdot M_3(-u) \tag{2.37}$$

*where $\Omega_{\mathrm{node}}$ is the right ascension of the ascending node, $i$ is the inclination, and $u = \omega + \nu$ is the argument of latitude. This is a pure rotation — no translation, no frame angular velocity:*

$$T_T = 0, \qquad \Omega_T = 0$$

*The velocity transforms identically to position (no transport term):*

$$R_{\mathrm{pos,TEME}} = M\,R_{\mathrm{pos,PF}}\,M^\dagger, \qquad V_{\mathrm{vel,TEME}} = M\,V_{\mathrm{vel,PF}}\,M^\dagger \tag{2.38}$$

*Proof.* The perifocal frame is obtained from TEME by: (1) rotating about $\hat{\mathbf{z}}$ by $\Omega_{\mathrm{node}}$ to align with the ascending node, (2) rotating about $\hat{\mathbf{x}}$ by $i$ to tilt into the orbital plane, (3) rotating about the new $\hat{\mathbf{z}}$ by $u$ to align with the satellite's position. The inverse (PF → TEME) reverses the angles: $M_3(-u)$ undoes step 3, $M_1(-i)$ undoes step 2, $M_3(-\Omega)$ undoes step 1. Composition follows from Theorem 2.2.6. Since $\Omega_T = 0$, the transport term in (2.30) vanishes, and velocity transforms by the same conjugation as position. ∎

**Remark** (Computational note)**.** Since $M_3(\alpha) = \left(\begin{smallmatrix} e^{i\alpha/2} & 0 \\ 0 & e^{-i\alpha/2} \end{smallmatrix}\right)$ is diagonal, the first and last factors in (2.37) are phase multiplications. The product $M_3(-\Omega) \cdot M_3(-u)$ at the outer positions could be combined: but because $M_1(-i)$ is between them, all three multiplications are needed. However, the diagonal structure means each $M_3$ multiplication costs only 4 real multiplications (instead of 16 for a general $SU(2)$ product).

**Proposition 2.5.1** (Singularity-free at $i = 0$)**.** *When $i = 0$ (equatorial orbit), $M_1(0) = I$ and the composition reduces to:*

$$M_{\mathrm{PF \to TEME}} = M_3(-\Omega) \cdot I \cdot M_3(-u) = M_3(-\Omega - u) \tag{2.39}$$

*The individual angles $\Omega$ and $u$ are undefined (the ascending node is arbitrary for an equatorial orbit), but their sum $\Omega + u$ is well-defined (it is the satellite's longitude). The $SU(2)$ product (2.39) is non-singular — a single well-conditioned diagonal matrix.*

*Similarly, when $e = 0$ (circular orbit), $\omega$ is undefined but $u = \omega + \nu$ is well-defined. The rotation $M_3(-u)$ is non-singular.*

*No conditional guards or Lyddane-type modifications are needed in the $SU(2)$ representation.*

*Proof.* $M_1(0) = \cos(0)\,I + i\sin(0)\,\sigma_1 = I$. Then $M_3(-\Omega) \cdot I \cdot M_3(-u) = M_3(-\Omega) \cdot M_3(-u)$. Since both are diagonal with entries $e^{\pm i\alpha/2}$, their product is $M_3(-\Omega-u)$, which has entries $e^{\pm i(\Omega+u)/2}$. These are well-defined for any value of $\Omega + u$. ∎

### Transform 2: TEME → PEF

**Theorem 2.5.2** (TEME → PEF with transport)**.** *The transformation from TEME to PEF is a rotation by $\theta_{\mathrm{GMST}}$ about $\hat{\mathbf{e}}_3$, with the Earth rotating at angular velocity $\boldsymbol{\omega}_E = \omega_E\hat{\mathbf{e}}_3$:*

$$M_T = M_3(\theta_{\mathrm{GMST}}) = \begin{pmatrix} e^{i\theta/2} & 0 \\ 0 & e^{-i\theta/2} \end{pmatrix} \tag{2.40}$$

$$\Omega_T = \omega_E\sigma_3 = \begin{pmatrix} \omega_E & 0 \\ 0 & -\omega_E \end{pmatrix} \tag{2.41}$$

*No translation: $T_T = 0$. The velocity transform (from Theorem 2.4.1, equation 2.30):*

$$V_{\mathrm{vel,PEF}} = M_T\,V_{\mathrm{vel,TEME}}\,M_T^\dagger + \frac{1}{2i}[\Omega_T, R_{\mathrm{pos,PEF}}] \tag{2.42}$$

*The transport commutator has a particularly simple form because both $M_T$ and $\Omega_T$ are diagonal. For $R_{\mathrm{pos,PEF}} = \left(\begin{smallmatrix} r_3 & r_- \\ r_+ & -r_3 \end{smallmatrix}\right)$ where $r_\pm = r_1 \pm ir_2$:*

$$\frac{1}{2i}[\Omega_T, R_{\mathrm{pos,PEF}}] = \frac{\omega_E}{2i}\begin{pmatrix} 0 & 2r_- \\ -2r_+ & 0 \end{pmatrix} = \omega_E\begin{pmatrix} 0 & -ir_- \\ ir_+ & 0 \end{pmatrix} \tag{2.43}$$

*The transport velocity is $\boldsymbol{\omega}_E \times \mathbf{r} = (-\omega_E r_2, \omega_E r_1, 0)$.* [M.2.1]

*Proof.* The angular velocity of Earth rotation is $\Omega_T = \omega_E\sigma_3$ (Definition 2.2.3). Write $R_{\mathrm{pos}} = \left(\begin{smallmatrix} r_3 & r_- \\ r_+ & -r_3 \end{smallmatrix}\right)$ where $r_\pm = r_1 \pm ir_2$. The commutator $[\sigma_3, R_{\mathrm{pos}}]$ is:

$$\sigma_3 R_{\mathrm{pos}} - R_{\mathrm{pos}} \sigma_3 = \begin{pmatrix} r_3 & r_- \\ -r_+ & r_3 \end{pmatrix} - \begin{pmatrix} r_3 & -r_- \\ r_+ & r_3 \end{pmatrix} = \begin{pmatrix} 0 & 2r_- \\ -2r_+ & 0 \end{pmatrix}$$

Applying the factor $\frac{\omega_E}{2i}$:

$$\frac{\omega_E}{2i}\begin{pmatrix} 0 & 2r_- \\ -2r_+ & 0 \end{pmatrix} = \omega_E\begin{pmatrix} 0 & -ir_- \\ ir_+ & 0 \end{pmatrix}$$

The $(1,2)$ entry is $-i\omega_E r_- = -i\omega_E(r_1 - ir_2) = -\omega_E r_2 - i\omega_E r_1$. By Definition 2.2.3, this encodes $v_1 - iv_2$, so $v_1 = -\omega_E r_2$ and $v_2 = \omega_E r_1$. The diagonal is zero, so $v_3 = 0$. Therefore $\mathbf{v}_{\mathrm{transport}} = (-\omega_E r_2,\; \omega_E r_1,\; 0) = \boldsymbol{\omega}_E \times \mathbf{r}$. ∎

### Pipeline Composition

**Theorem 2.5.3** (Perifocal → PEF composition)**.** *The complete transform from perifocal to PEF is:*

$$M_{\mathrm{PF \to PEF}} = M_3(\theta_{\mathrm{GMST}}) \cdot M_3(-\Omega_{\mathrm{node}}) \cdot M_1(-i) \cdot M_3(-u) \tag{2.44}$$

*Since the first two factors are both $M_3$ (diagonal), they compose trivially:*

$$M_3(\theta) \cdot M_3(-\Omega) = M_3(\theta - \Omega) \tag{2.45}$$

*reducing the composition to three $SU(2)$ multiplications (one of which is a diagonal-by-general product). The velocity includes the transport term from the TEME → PEF step:*

$$V_{\mathrm{vel,PEF}} = M_{\mathrm{total}}\,V_{\mathrm{vel,PF}}\,M_{\mathrm{total}}^\dagger + \frac{1}{2i}[\Omega_E, R_{\mathrm{pos,PEF}}] \tag{2.46}$$

*The transport term acts on the already-rotated position $R_{\mathrm{pos,PEF}} = M_{\mathrm{total}}\,R_{\mathrm{pos,PF}}\,M_{\mathrm{total}}^\dagger$.*

*Proof.* Composition of the two rotations (Theorem 2.2.6) gives (2.44). The velocity transform: Apply Theorem 2.4.1 (equation 2.30) to the TEME → PEF step, with the input velocity already rotated from perifocal to TEME by Transform 1. Since Transform 1 has $\Omega_T = 0$, no transport term is added in that step. The only transport contribution comes from Transform 2 (Earth rotation). ∎

### Error Sources

**Remark** (Error source classification)**.** *The following table traces each error source through the $SU(2)$ construction. The "Symbolic bound" column is precision-independent; the "Binary64 magnitude" column evaluates it for IEEE 754 double precision ($\epsilon_{\mathrm{mach}} = 2^{-52}$). Other arithmetic precisions scale accordingly.*

| Source | Category | Symbolic bound | Propagation | Binary64 magnitude |
|--------|----------|---------------|-------------|-------------------|
| $\Omega_{\mathrm{node}}$ from TLE | $\sigma_m$ | $\|\delta M\| \leq \frac{1}{2}\delta\Omega$ (Thm 2.2.12) | Conjugation of all vectors | TLE-dependent |
| $i$ from TLE | $\sigma_m$ | $\|\delta M\| \leq \frac{1}{2}\delta i$ | Same | TLE-dependent |
| $u$ from TLE | $\sigma_m$ | $\|\delta M\| \leq \frac{1}{2}\delta u$ | Same | TLE-dependent |
| sin/cos evaluation | $\delta_p$ | $\leq \epsilon_{\mathrm{mach}}$ per entry (Ch 1, Cors. 1.4.1–1.4.2) | $M_k$ entries | $\sim 10^{-16}$ |
| $SU(2)$ multiply rounding | $\delta_p$ | $\leq 8\epsilon_{\mathrm{mach}}$ per entry [P.2.3] | Each composition | $\sim 10^{-15}$ per multiply |
| $\theta_{\mathrm{GMST}}$ polynomial | $\delta_p$ | Angle sensitivity (Thm 2.2.12) | $M_3(\theta)$ for TEME→PEF | Ch 29 |
| GMST coefficients | $\sigma_m$ | Through $M_3(\theta)$ | $\theta_{\mathrm{GMST}}$ value | IAU 1982 |
| $\omega_E$ | $\sigma_m$ | Through commutator (Cor 2.4.1) | $\Omega_T = \omega_E\sigma_3$ | $\sim 10^{-12}$ rad/s |
| TEME frame definition | $\delta_a$ | Irreducible floor on output | Entire frame | $\sim 0.1$ arcsec [A.2.1] |

## §2.6 Error Propagation Through the Framework

### The Two Principles

**Theorem 2.6.1** (Principle 1: linear operations reuse the transform)**.** *For any linear operation $\mathcal{L}: V \mapsto \mathcal{L}(V)$ in the framework (conjugation $MVM^\dagger$, commutator $\frac{1}{2i}[\Omega, V]$, addition, scalar multiplication), the error transforms by the same operation:*

$$\delta(\mathcal{L}(V)) = \mathcal{L}(\delta(V)) \tag{2.47}$$

*This is exact — not a first-order approximation. No separate error formula is needed. The same computation, applied to the error state, gives the error of the output.*

*Proof.* For a linear map $\mathcal{L}$: $\mathcal{L}(V + \delta V) = \mathcal{L}(V) + \mathcal{L}(\delta V)$. Therefore $\delta(\mathcal{L}(V)) = \mathcal{L}(V + \delta V) - \mathcal{L}(V) = \mathcal{L}(\delta V)$. This holds exactly for linear maps (Ch 1, §1.3: addition and scalar multiplication are linear; §2.2: conjugation is linear in $V$ for fixed $M$; §2.4: the commutator is linear in each argument). ∎

**Theorem 2.6.2** (Principle 2: nonlinear operations require sensitivity bounds)**.** *For a nonlinear operation $g: V \mapsto g(V)$ (e.g., constructing $M$ from an angle, computing a square root, evaluating sin/cos), the error satisfies:*

$$|\delta(g(V))| \leq \sup|g'| \cdot |\delta(V)| \tag{2.48}$$

*where $|g'|$ is the sensitivity of $g$ to its input, bounded using Ch 1, Theorem 1.4.1 (derivative bound for scalar functions, valid for both $\mathbb{R}$ and $\mathbb{C}$). The bound must be rigorous (a proven upper bound on $|g'|$), not an approximation.* [A.2.2]

*Proof.* This is Theorem 1.4.1 applied to $g$. The supremum is over the error disk $|z - V| \leq |\delta(V)|$. For holomorphic $g$ (which includes all standard functions), the line-integral proof of Theorem 1.4.1 gives the bound. ∎

**Remark.** The distinction between Principles 1 and 2 is: the APPLICATION of a transform to the state is linear (Principle 1); the CONSTRUCTION of the transform from parameters (angles, time polynomials) is generally nonlinear (Principle 2). In the orbit propagation pipeline, most operations are applications — rotating vectors, adding corrections, composing transforms. The nonlinear steps are concentrated at the interfaces: computing sin/cos from angles, solving Kepler's equation, evaluating the density model.

### Pipeline Classification

**Theorem 2.6.3** (Complete pipeline classification)**.** *Every operation in the orbit propagation pipeline is classified as linear or nonlinear for error propagation purposes:*

**Linear operations** (error row uses same transform — Principle 1):

| Operation | Expression | Why linear |
|-----------|-----------|-----------|
| Rotation of any vector | $MVM^\dagger$ | Bilinear in $V$; linear for fixed $M$ (Thm 2.2.10) |
| Translation | $+ T$ | Constant addition (Ch 1, Thm 1.3.1) |
| Transport term | $\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ | Commutator is linear in $R_{\mathrm{pos}}$ (Thm 2.4.3) |
| $SU(2)$ composition | $M_1 \cdot M_2$ | Linear in $M_2$ for fixed $M_1$ (matrix multiplication) |
| Secular rate additions | $M_0 + \dot{M}\Delta t$ | Linear in both terms |
| Short/long-period corrections | $+ \delta V$ | Addition (Ch 1, Thm 1.3.1) |
| Dot product | $\frac{1}{2}\mathrm{tr}(AB)$ | Bilinear (Thm 2.2.2) |
| Cross product | $\frac{1}{2i}[A, B]$ | Bilinear (Thm 2.2.1) |

**Nonlinear operations** (error row requires rigorous sensitivity bound — Principle 2):

| Operation | Why nonlinear | Sensitivity | Bound source |
|-----------|--------------|-------------|-------------|
| Angle → $M_k(\alpha)$ | sin/cos are nonlinear | $\|dM_k/d\alpha\| = \frac{1}{2}$ | Thm 2.2.12 |
| Kepler equation solver | Iterative, transcendental | $\leq 1/(1-e)$ | Ch 9 |
| sin/cos of orbital angles | Transcendental | $|\cos x|\delta + \delta^2/2$ | Ch 1, Cors. 1.4.1–1.4.2 (real domain) |
| Renormalization $M/\sqrt{\det M}$ | Division by norm | $\leq 1 + O(\det M - 1)$ | See below |
| Cube root (element recovery) | Nonlinear root | $\frac{2}{3}x^{-1/3}$ | Ch 32 |
| Density model $(q_0-s)/(r-s))^\tau$ | Power law | $|\tau|\rho/(r-s)$ | Ch 21 |
| Euclidean norm $\sqrt{-\det V}$ | Square root | $\leq 1$ (norm never amplifies) | Ch 1, Cor. 1.4.4 |

**Cross product cancellation risk:** When $\mathbf{a} \approx k\mathbf{b}$ (nearly parallel), $\frac{1}{2i}[A, B] \approx 0$ and the subtraction in the commutator loses digits. The condition number is $\kappa \approx \|A\|\|B\| / \|\frac{1}{2i}[A,B]\|$, detected by $\mathrm{rd}(\mathbf{a} \times \mathbf{b}) \to 0$ (Ch 1, §1.7). This occurs in angular momentum $\mathbf{h} = \mathbf{r} \times \mathbf{v}$ when the orbit is nearly rectilinear.

### Renormalization Bound

**Theorem 2.6.4** (Renormalization sensitivity)**.** *For $M$ near-unitary with $\det M = 1 + \varepsilon$ where $|\varepsilon| \ll 1$:*

$$M_{\mathrm{renorm}} = \frac{M}{\sqrt{\det M}} \approx M(1 - \varepsilon/2) \tag{2.49}$$

*The correction is $O(\varepsilon)$ and the sensitivity of the renormalization to the input is $\leq 1 + O(\varepsilon)$, which is near-unity for small drift.*

*After $n$ $SU(2)$ multiplications, $|\varepsilon| \leq n \cdot 8\epsilon_{\mathrm{mach}}$. The need for renormalization depends on the arithmetic precision and the number of compositions. For binary64 over a full propagation year ($\sim 2 \times 10^6$ multiplications at worst): $|\varepsilon| \leq 2 \times 10^{-9}$, far below typical measurement error floors. For narrower arithmetic types, renormalization may be needed sooner.*

*Proof.* $\sqrt{\det M} = \sqrt{1 + \varepsilon} = 1 + \varepsilon/2 + O(\varepsilon^2)$. Division by $1 + \varepsilon/2$ gives $M/(1 + \varepsilon/2) \approx M(1 - \varepsilon/2)$. The correction to each matrix entry is $\leq |\varepsilon|/2 \cdot |M_{ij}| \leq |\varepsilon|/2$. ∎

### Linearization Validity

**Remark** (When the parallel-row principle breaks down)**.** The parallel-row principle (Theorem 2.6.1) is exact for linear operations. For nonlinear operations, Theorem 2.6.2 uses a first-order bound. This bound is valid when the error is small relative to the value — specifically, when $\mathrm{rd}(v) \geq d_{\min}$ (Ch 1, Definition 1.9.1). The threshold $d_{\min} = 4$ means the relative error is $\leq 10^{-4}$, and the neglected second-order term (Ch 1, Theorem 1.3.2: the cross term $\delta_x \delta_y$) is $\leq 10^{-8}$ of the first-order terms.

When $\mathrm{rd}(v) < d_{\min}$: the linearization may not be valid, and full nonlinear error propagation from Ch 1 (Theorems 1.3.1–1.3.3, applied to each scalar operation individually) must be used. This should trigger a formula-switching event (Ch 1, Definition 1.9.3).

### Three-Error Independence

**Theorem 2.6.5** (Error categories propagate independently)**.** *The three error categories $\sigma_m$, $\delta_p$, $\delta_a$ propagate through all framework operations independently. Each category uses the same functional form for the bound — the same linear map (Principle 1) or the same sensitivity bound (Principle 2) — but with different input magnitudes.*

*Proof.* Each propagation rule (Theorems 2.6.1–2.6.2) has the form $\delta_{\mathrm{out}} = g(|\mathrm{values}|, \delta_{\mathrm{in}})$ where $g$ is the same function regardless of which error category $\delta$ represents. The three categories are independent because they arise from disjoint physical mechanisms (Ch 1, Proposition 1.8.1). ∎

---

## §2.7 The Traditional Matrix Form

**Theorem 2.7.1** (Rotation matrix from $SU(2)$)**.** *The $3 \times 3$ rotation matrix $R$ corresponding to $M = \left(\begin{smallmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{smallmatrix}\right)$ is:*

$$R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j\,M\,\sigma_k\,M^\dagger) \tag{2.50}$$

*With $\alpha = w + iz$, $\beta = ix + y$ (i.e., the quaternion $(w,x,y,z)$):*

$$R = \begin{pmatrix} 1 - 2(y^2+z^2) & 2(xy-wz) & 2(xz+wy) \\ 2(xy+wz) & 1-2(x^2+z^2) & 2(yz-wx) \\ 2(xz-wy) & 2(yz+wx) & 1-2(x^2+y^2) \end{pmatrix} \tag{2.51}$$

*Every entry is a quadratic polynomial in the quaternion components. This IS the standard rotation matrix.*

*Proof.* The map $V \mapsto MVM^\dagger$ is a linear map on the 3D vector space of traceless Hermitian matrices. In the Pauli basis $\{\sigma_1, \sigma_2, \sigma_3\}$, this map has the matrix representation $R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j M \sigma_k M^\dagger)$ (since the component of $MVM^\dagger$ along $\sigma_j$ is $\frac{1}{2}\mathrm{tr}(\sigma_j MVM^\dagger) = \frac{1}{2}\sum_k v_k\mathrm{tr}(\sigma_j M\sigma_k M^\dagger) = \sum_k R_{jk} v_k$). The explicit entries (2.51) follow from expanding $\mathrm{tr}(\sigma_j M\sigma_k M^\dagger)$ with $M$ in the $(\alpha, \beta)$ parameterization. ∎

**Theorem 2.7.2** (The $7 \times 7$ state matrix from $SU(2)$)**.** *Writing the state transform (Theorem 2.4.1) in Cartesian components, the complete transformation is:*

$$\begin{pmatrix} \mathbf{r}' \\ \mathbf{v}' \\ 1 \end{pmatrix} = \begin{pmatrix} R & \mathbf{0} & \mathbf{t} \\ [\boldsymbol{\omega}_T]_\times R & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \begin{pmatrix} \mathbf{r} \\ \mathbf{v} \\ 1 \end{pmatrix} \tag{2.52}$$

*where:*
- *$R = R(M_T)$ from (2.51): the $3 \times 3$ rotation matrix from the $SU(2)$ element*
- *$[\boldsymbol{\omega}_T]_\times$: the $3 \times 3$ skew-symmetric matrix with $[\boldsymbol{\omega}]_\times\mathbf{v} = \boldsymbol{\omega} \times \mathbf{v}$, corresponding to the commutator $\frac{1}{2i}[\Omega_T, \cdot]$ written in coordinates*
- *$\mathbf{t}$: the translation vector from $T_T = 2D_T M_T^\dagger$*
- *$\Delta\mathbf{v}$: the velocity offset (zero for all standard orbit propagation transforms)*

*This is exactly the $7 \times 7$ state matrix from the archived DCM chapter. The matrix form is the $SU(2)$ algebra written in Cartesian coordinates.*

*Proof.* The position row: $\mathbf{r}' = R\mathbf{r} + \mathbf{t}$, which is equation (2.29) written in components. The velocity row: $\mathbf{v}' = R\mathbf{v} + [\boldsymbol{\omega}_T]_\times R\mathbf{r} + \Delta\mathbf{v}$. The transport term: the commutator $\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$ in Cartesian components is $\boldsymbol{\omega}_T \times \mathbf{r}'$, and $\mathbf{r}' = R\mathbf{r} + \mathbf{t}$. For the standard orbit propagation transforms, $\mathbf{t} = 0$ so $\boldsymbol{\omega}_T \times \mathbf{r}' = [\boldsymbol{\omega}_T]_\times R\mathbf{r}$, which is the lower-left block $[\boldsymbol{\omega}_T]_\times R$. ∎

**Remark** (When to use which form)**.** The $SU(2)$ form (§§2.2–2.5) is preferable for: singularity analysis (Proposition 2.5.1), compact composition (equation 2.44), error architecture (§2.6 Principle 1), and the formal framework. The $3 \times 3$ / $7 \times 7$ matrix form is preferable for: connection to existing orbit propagation literature (Vallado, Hoots & Roehrich) and the component-wise computations that downstream chapters (Ch 30) require. The code-to-theorem mapping for specific source modules is in Appendix C. Both forms produce numerically identical results — the matrix form IS the $SU(2)$ form expanded in coordinates.

---

## §2.8 The Parallel-Row Error Architecture

**Definition 2.8.1** (Parallel-row state)**.** *The complete tracked state is a set of four rows, each with 13 real components and the same structure:*

| Row | Content | Components |
|-----|---------|------------|
| 1 | Physical state | $M$ (4), $R_{\mathrm{pos}}$ (3), $V_{\mathrm{vel}}$ (3), $\Omega$ (3) |
| 2 | Measurement error $\sigma_m$ | $\delta_m(M)$ (4), $\delta_m(R_{\mathrm{pos}})$ (3), $\delta_m(V)$ (3), $\delta_m(\Omega)$ (3) |
| 3 | Precision error $\delta_p$ | $\delta_p(M)$ (4), $\delta_p(R_{\mathrm{pos}})$ (3), $\delta_p(V)$ (3), $\delta_p(\Omega)$ (3) |
| 4 | Accuracy error $\delta_a$ | $\delta_a(M)$ (4), $\delta_a(R_{\mathrm{pos}})$ (3), $\delta_a(V)$ (3), $\delta_a(\Omega)$ (3) |

*Total storage: 4 × 13 = 52 real numbers.*

*All four rows pass through the same transform pipeline. For linear operations (the majority of the orbit propagation pipeline), the same code processes each row. For nonlinear operations, the error rows are multiplied by the rigorous sensitivity bound (Theorem 2.6.2).*

**Remark** (Merged mode)**.** For production propagation where the error category breakdown is not needed, the three error rows may be merged into a single total error row: $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$ per component (Ch 1, Definition 1.8.1). This reduces the storage to 2 × 13 = 26 real numbers and requires one extra transform pass instead of three.

### Scalar Error Reduction

**Definition 2.8.2** (Scalar error summaries)**.** *The per-component errors in each row are reduced to four scalar summaries using the algebraic norm (Theorem 2.2.3):*

$$\delta(\text{position}) = \sqrt{-\det(\delta R_{\mathrm{pos}})} = \sqrt{\delta r_1^2 + \delta r_2^2 + \delta r_3^2} \tag{2.53}$$

$$\delta(\text{velocity}) = \sqrt{-\det(\delta V_{\mathrm{vel}})} = \sqrt{\delta v_1^2 + \delta v_2^2 + \delta v_3^2} \tag{2.54}$$

$$\delta(\text{orientation}) = \|\boldsymbol{\varepsilon}\| \;\text{where}\; \delta M = \frac{i}{2}\varepsilon_k\sigma_k \cdot M \;\text{(full angle, in radians — see §R.3)} \tag{2.55}$$

$$\delta(\text{angular velocity}) = \sqrt{-\det(\delta\Omega)} = \sqrt{\delta\omega_1^2 + \delta\omega_2^2 + \delta\omega_3^2} \tag{2.56}$$

*All four use the same algebraic operation: $\sqrt{-\det(\cdot)}$ applied to a traceless Hermitian error matrix (for vectors), or the perturbation magnitude (for orientation).*

**Remark** (Connection to reliable digits)**.** The scalar summaries feed into Ch 1, Definition 1.9.1:

$$\mathrm{rd}(\text{position}) = \left\lfloor -\log_{10}\!\left(\frac{\delta(\text{position})}{\|\mathbf{r}\|}\right) \right\rfloor$$

This is the decision criterion for formula switching (when $\mathrm{rd} < d_{\min}$, switch to an alternative formulation) and for output quality assessment (Ch 38).

### Point-Mass Reduction

**Corollary 2.8.1.** *For point-mass mode ($M = $ frame rotation, $\Omega = 0$):*
- *The orientation error $\delta(M)$ reflects frame rotation accuracy (from TLE angles)*
- *The angular velocity error $\delta(\Omega) = 0$ (trivially)*
- *The active error components are: 4 (orientation) + 3 (position) + 3 (velocity) = 10 per error row*
- *The storage for the complete tracked state: 13 (physical) + 3 × 10 (active error components) = 43 real numbers, or 13 + 10 = 23 in merged mode*

---

## §2.9 Summary: How Subsequent Chapters Use This Framework

Every subsequent chapter with a "State Matrix Formulation" section follows this pattern:

1. Express the chapter's operation as an $SU(2)$ transform or a traceless Hermitian matrix operation (§§2.2–2.4)
2. Classify the operation as linear or nonlinear for error propagation (§2.6, Theorem 2.6.3)
3. For linear: state that Principle 1 (Theorem 2.6.1) applies — the error row uses the same transform
4. For nonlinear: provide the rigorous sensitivity bound (Theorem 2.6.2) with proof
5. Show how the operation composes with adjacent pipeline steps (§2.4, Theorem 2.4.1)
6. Identify error sources and their categories using the classification of §2.5

### Downstream Chapter Reference

| Chapter | Framework aspect used | This chapter's section |
|---------|---------------------|----------------------|
| Ch 8 (Keplerian orbit) | Perifocal frame rotation | §2.5, Thm 2.5.1 |
| Ch 13 (Geopotential) | Force gradient as traceless Hermitian perturbation | §2.4, Thm 2.4.1 |
| Ch 15 (Kaula expansion) | Inclination/eccentricity functions in the framework | §2.4 |
| Ch 16–17 (Brouwer secular) | Secular rates as near-identity perturbations | §2.4, Thm 2.4.4 |
| Ch 18 (Short-period) | Additive corrections (linear, Principle 1) | §2.6, Thm 2.6.1 |
| Ch 20 (Osculating elements) | Complete mean-to-osculating transform as composition | §2.5, Thm 2.5.3 |
| Ch 27 (Third-body) | Third-body force as traceless Hermitian perturbation | §2.4 |
| Ch 29 (Sidereal time) | GMST as input to TEME→PEF rotation | §2.5, Thm 2.5.2 |
| Ch 30 (Coordinate transforms) | All frame rotations and compositions | §§2.2, 2.3, 2.5 |
| Ch 34–35 (Propagation) | Full pipeline as chain of $SU(2)$ products | §2.5, Thm 2.5.3 |
| Ch 38 (State vector output) | Error budget from parallel-row architecture | §§2.6, 2.8 |

---

## Error Notes

**[P.2.1]** Shepperd extraction square root. The branch selection in Theorem 2.2.9 ensures the argument of the square root is $\geq \frac{1}{2}$, bounding the precision error by $\delta \leq \delta(R_{jk})/(2\sqrt{1/2}) = \sqrt{2}\,\delta(R_{jk})$ (Ch 1, Corollary 1.4.4). *Remedy:* use the largest-diagonal branch to maximize the denominator.

**[P.2.2]** Inexact rotation entry error. Each entry of $M_k(\alpha)$ is a trigonometric function of $\alpha/2$, evaluated in finite precision. The error per entry is bounded by $\frac{1}{2}|\delta(\alpha)| + \epsilon_{\mathrm{mach}}$ (Theorem 2.2.12 plus representation error). For binary64 with $|\delta(\alpha)| \sim 10^{-8}$ (TLE precision): the angle contribution dominates at $\sim 5 \times 10^{-9}$; the representation error is $\sim 10^{-16}$. *Remedy:* improve angle measurement (not a computational concern).

**[P.2.3]** Composition rounding. Each $SU(2)$ matrix multiplication involves 8 complex multiply-adds = 16 real operations. Rounding accumulation per multiply: $\leq 8\epsilon_{\mathrm{mach}}$ per entry. For 2 compositions (the 3-1-3 sequence): total rounding $\leq 16\epsilon_{\mathrm{mach}}$. For binary64: $\approx 2 \times 10^{-15}$. For single precision: $\approx 2 \times 10^{-6}$. *Remedy:* use wider arithmetic type to reduce $\epsilon_{\mathrm{mach}}$.

**[A.2.1]** TEME frame accuracy. The TEME frame differs from the precise GCRS/ITRS frames by $\sim 0.1$ arcsec ($\sim 50$ m at GEO, $\sim 0.3$ m at LEO) due to the truncated precession/nutation model. This is an accuracy error — it cannot be reduced by computational precision. *Remedy:* apply IAU 2006/2000A precession-nutation corrections (beyond the standard propagator model; breaks the matched pair — see Ch 3).

**[M.2.1]** GMST determination. The GMST polynomial coefficients (Ch 29) carry measurement error from the IAU 1982 determination (Aoki et al.). The angle $\theta_{\mathrm{GMST}}$ propagates this error into the TEME→PEF rotation via Theorem 2.2.12. *Remedy:* use updated IERS values (but this breaks the matched pair — see Ch 3).

**[A.2.2]** Linearization in error propagation. The parallel-row principle (Theorem 2.6.1) is exact for linear operations. For nonlinear operations, Theorem 2.6.2 uses a first-order sensitivity bound. The neglected second-order terms are $O(\delta^2)$, which is $\leq 10^{-8}$ of the first-order bound when $\mathrm{rd} \geq 4$. When $\mathrm{rd} < 4$, full nonlinear propagation from Ch 1 must be used. *Remedy:* reduce to Ch 1 scalar error rules at each operation when the reliable-digits criterion is not met.

**[P.2.4]** Position extraction coupling. Extracting position via $T = 2DM^\dagger$ (Theorem 2.3.6) has error $\leq \delta_D + \frac{1}{2}\|\mathbf{r}\|\delta_M$. The rotation-coupled term $\frac{1}{2}\|\mathbf{r}\|\delta_M$ scales with orbital radius and rotation error. For a representative case ($\|\mathbf{r}\| = 7000$ km, $\delta_M \sim 5 \times 10^{-9}$ from TLE angle precision): $\sim 0.035$ m. The actual magnitude depends on the input source precision. *Remedy:* improve angle measurement (the rotation error $\delta_M$ comes from input angle uncertainties).
