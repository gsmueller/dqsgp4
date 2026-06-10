# Chapter 2: The State Matrix Framework (Alternate — Two-Layer Architecture)

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\sigma_k$ | Pauli matrices ($k = 1, 2, 3$) | §2.2, Def. 2.2.1 |
| $I$ | The $2 \times 2$ identity matrix | §2.2 |
| $A^\dagger$ | Conjugate transpose (Definition 2.2.2) | §2.2, Def. 2.2.2 |
| $M$ | An element of $SU(2)$: unitary with $\det M = 1$ | §2.2, Def. 2.2.6 |
| $\alpha$, $\beta$ | Complex parameters of $M$ satisfying $\lvert\alpha\rvert^2+\lvert\beta\rvert^2=1$ | §2.2, Def. 2.2.6 |
| $V$ | A traceless Hermitian $2 \times 2$ matrix encoding a 3D vector | §2.2, Def. 2.2.5 |
| $\delta_{jk}$ | Kronecker delta: $1$ if $j=k$, else $0$ | §2.2, Def. 2.2.9 |
| $\varepsilon_{jkl}$ | Levi-Civita symbol: fully antisymmetric tensor on $\{1,2,3\}$ | §2.2, Def. 2.2.9 |
| $R(M)$ | The $3 \times 3$ rotation matrix derived from $M$ via the adjoint representation | §2.3, Def. 2.3.1 |
| $H$ | A $4 \times 4$ homogeneous matrix: rotation $R$ + translation $\mathbf{t}$ | §2.4, Def. 2.4.1 |
| $T_7$ | The $7 \times 7$ state matrix: rotation + translation + transport term | §2.5, Def. 2.5.2 |
| $[\mathbf{a}]_\times$ | The $3 \times 3$ skew-symmetric matrix: $[\mathbf{a}]_\times \mathbf{b} = \mathbf{a} \times \mathbf{b}$ | §2.5, Def. 2.5.1 |
| $\omega_E$ | Earth's rotation rate | §2.6 |
| $\varepsilon$ | Dual unit satisfying $\varepsilon^2 = 0$ (distinct from $\varepsilon_{jkl}$) | §2.8, Def. 2.8.1 |
| $\hat{M}$ | Dual $SU(2)$ element: $\hat{M} = M + \varepsilon D$ | §2.8, Def. 2.8.3 |

---

## §2.1 Introduction

The state matrix framework provides a unified representation for coordinate transforms in orbit propagation. Every transform — from the perifocal frame to inertial, from inertial to Earth-fixed — is expressed as a single matrix that acts on the state vector by left-multiplication.

The framework has two layers, each serving a distinct purpose:

**Layer 1 — Construction.** Rotations are built and composed using the group $SU(2)$ of $2 \times 2$ unitary matrices with unit determinant. This representation is singularity-free (no gimbal lock), compact (4 real parameters with 1 constraint), and composes by matrix multiplication. Every 3D vector — position, velocity, angular velocity — can be encoded as a $2 \times 2$ traceless Hermitian matrix, and the rotation acts by conjugation: $V' = MVM^\dagger$. The cross product, dot product, and norm emerge as native operations in this algebra (§2.2).

**Layer 2 — Action.** To act on state vectors, the $SU(2)$ element $M$ is converted to a $3 \times 3$ rotation matrix $R(M)$ via the adjoint representation (§2.3). This conversion is exact and done once per transform. The rotation matrix is then embedded in a $4 \times 4$ homogeneous matrix (with translation) or a $7 \times 7$ state matrix (with translation and velocity coupling). The action on state vectors is left-multiplication — truly linear, with no affine exceptions. The same matrix acts on both the physical state and the error state (§2.7).

The dual extension ($SU(2)$ with dual numbers) provides an alternative composition mechanism for rotation + translation in Layer 1, connecting to the rich theory of dual quaternions (§2.8). However, because the dual extension produces an affine action on vectors rather than a linear one, the operational pipeline uses the homogeneous matrix form of Layer 2 for the actual state propagation.

The separation of concerns is:

| | Layer 1: $SU(2)$ | Layer 2: Matrix |
|---|---|---|
| **Purpose** | Build and compose transforms | Act on state vectors |
| **Rotation** | $M_1 M_2$ (quaternion multiply) | $R(M)\,\mathbf{v}$ (matrix-vector) |
| **Translation** | Dual extension $\hat{M}$ | Homogeneous column in $H$ |
| **Velocity coupling** | Transport theorem (§2.5) | $\dot{R}$ block in $T_7$ |
| **Singularity** | None | None (inherits from $M$) |
| **Error propagation** | Conjugation is bilinear | Left-multiply is linear: same $T_7$ on both rows |
| **Conversion** | — | $M \to R(M)$: once per transform (§2.3) |

---

## §2.2 The $SU(2)$ Construction Layer

This section develops the algebraic tools for constructing and composing rotations. All results here operate in the $2 \times 2$ complex matrix world. The conversion to the operational $3 \times 3$ form is deferred to §2.3.

### Foundations

**Definition 2.2.1** (Pauli matrices)**.** *The three Pauli matrices are:*

$$\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} \tag{2.1}$$

*Together with the identity $I$, they form a basis for the real vector space of $2 \times 2$ Hermitian matrices.*

**Definition 2.2.2** (Conjugate transpose)**.** *For $A \in \mathbb{C}^{m \times n}$, the conjugate transpose $A^\dagger$ is defined by $(A^\dagger)_{jk} = \overline{a_{kj}}$. Key properties: $(AB)^\dagger = B^\dagger A^\dagger$, $(A^\dagger)^\dagger = A$, $\det(A^\dagger) = \overline{\det(A)}$.*

**Definition 2.2.3** (Hermitian matrix)**.** *A matrix $A \in \mathbb{C}^{2 \times 2}$ is Hermitian if $A = A^\dagger$. The constraint $A = A^\dagger$ forces the diagonal entries to satisfy $a_{jj} = \overline{a_{jj}}$, i.e., $\mathrm{Im}(a_{jj}) = 0$: the diagonal entries are real. The off-diagonal entries satisfy $a_{12} = \overline{a_{21}}$. A general $2 \times 2$ Hermitian matrix therefore has the form:*

$$A = \begin{pmatrix} a & b \\ \bar{b} & c \end{pmatrix} \tag{2.2}$$

*where $a, c \in \mathbb{R}$ and $b \in \mathbb{C}$. This is a 4-dimensional real vector space.*

**Definition 2.2.4** (Traceless Hermitian matrix)**.** *A Hermitian matrix $A$ is traceless if $\mathrm{tr}(A) = a + c = 0$, i.e., $c = -a$. Every traceless Hermitian $2 \times 2$ matrix has the form:*

$$A = \begin{pmatrix} a & b \\ \bar{b} & -a \end{pmatrix} \tag{2.3}$$

*where $a \in \mathbb{R}$ (forced by the Hermitian constraint) and $b \in \mathbb{C}$. This gives 3 real degrees of freedom ($a$, $\mathrm{Re}(b)$, $\mathrm{Im}(b)$), forming a 3-dimensional real vector space with basis $\{\sigma_1, \sigma_2, \sigma_3\}$.*

**Definition 2.2.5** (Vector encoding)**.** *A vector $\mathbf{v} = (v_1, v_2, v_3) \in \mathbb{R}^3$ is encoded as the traceless Hermitian matrix:*

$$V = v_1 \sigma_1 + v_2 \sigma_2 + v_3 \sigma_3 = \begin{pmatrix} v_3 & v_1 - iv_2 \\ v_1 + iv_2 & -v_3 \end{pmatrix} \tag{2.4}$$

*This is an isomorphism: every traceless Hermitian $2 \times 2$ matrix corresponds to exactly one vector in $\mathbb{R}^3$, and vice versa. Components are recovered by $v_k = \frac{1}{2}\mathrm{tr}(\sigma_k V)$.*

**Definition 2.2.6** (Rotation element)**.** *A rotation is an element $M \in SU(2)$: a $2 \times 2$ unitary matrix with unit determinant. In the $(\alpha, \beta)$ parameterization:*

$$M = \begin{pmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{pmatrix}, \qquad |\alpha|^2 + |\beta|^2 = 1 \tag{2.5}$$

**Definition 2.2.7** (Conjugation)**.** *The conjugation of a matrix $V$ by an invertible matrix $M$ is the operation $V \mapsto MVM^\dagger$. When $M \in SU(2)$ and $V$ is traceless Hermitian, the conjugation preserves both properties (Theorem 2.2.4) and the result encodes the rotated vector:*

$$V' = M\,V\,M^\dagger \tag{2.6}$$

**Definition 2.2.8** (Elementary rotations)**.** *Rotation by angle $\alpha$ about coordinate axis $\hat{\mathbf{e}}_k$:*

$$M_k(\alpha) = \cos\frac{\alpha}{2}\,I + i\sin\frac{\alpha}{2}\,\sigma_k \tag{2.7}$$

*Equivalently, via the matrix exponential:*

$$M_k(\alpha) = \exp\!\left(\frac{i\alpha}{2}\,\sigma_k\right) \tag{2.8}$$

*This is the exponential map from the Lie algebra $\mathfrak{su}(2)$ (spanned by $\{i\sigma_1/2,\, i\sigma_2/2,\, i\sigma_3/2\}$) to the Lie group $SU(2)$. The half-angle $\alpha/2$ arises because $SU(2)$ is a double cover of $SO(3)$ (Lemma 2.2.2).*

### Algebraic Properties

The following identities connect the $SU(2)$ algebra to standard vector operations. They are used in §2.3 to derive the rotation matrix $R(M)$, and again in §2.5 for the transport theorem.

**Definition 2.2.9** (Kronecker delta and Levi-Civita symbol)**.** *The Kronecker delta $\delta_{jk}$ equals $1$ when $j = k$ and $0$ otherwise. The Levi-Civita symbol $\varepsilon_{jkl}$ is the fully antisymmetric tensor on indices $\{1, 2, 3\}$:*

$$\varepsilon_{123} = \varepsilon_{231} = \varepsilon_{312} = +1, \qquad \varepsilon_{132} = \varepsilon_{213} = \varepsilon_{321} = -1 \tag{2.9}$$

*All other values (any repeated index) are zero. The Levi-Civita symbol encodes the cross product: $(\mathbf{a} \times \mathbf{b})_l = \sum_{jk} \varepsilon_{jkl}\,a_j\,b_k$.*

**Lemma 2.2.1** (Pauli algebra)**.** *The Pauli matrices satisfy:*

$$\sigma_j \sigma_k = \delta_{jk}\,I + i\,\varepsilon_{jkl}\,\sigma_l \tag{2.10}$$

*where summation over the repeated index $l$ is implied.*

*Proof.* Direct verification for each pair $(j,k)$. For example: $\sigma_1 \sigma_2 = \left(\begin{smallmatrix} 0 & 1 \\ 1 & 0 \end{smallmatrix}\right)\left(\begin{smallmatrix} 0 & -i \\ i & 0 \end{smallmatrix}\right) = \left(\begin{smallmatrix} i & 0 \\ 0 & -i \end{smallmatrix}\right) = i\sigma_3$. ∎

**Theorem 2.2.1** (Cross product as commutator)**.** *For vectors $\mathbf{a}, \mathbf{b}$ encoded as traceless Hermitian matrices $A, B$:*

$$\mathbf{a} \times \mathbf{b} \;\longleftrightarrow\; \frac{1}{2i}[A, B] \tag{2.11}$$

*where $[A, B] = AB - BA$ is the matrix commutator.*

*Proof.* Expand $[A, B] = \sum_{jk} a_j b_k [\sigma_j, \sigma_k] = \sum_{jk} a_j b_k \cdot 2i\varepsilon_{jkl}\,\sigma_l$ (from Lemma 2.2.1: $[\sigma_j, \sigma_k] = 2i\varepsilon_{jkl}\sigma_l$). Thus $\frac{1}{2i}[A,B] = \sum_l (\sum_{jk} \varepsilon_{jkl} a_j b_k)\,\sigma_l$, and the coefficient of $\sigma_l$ is $(\mathbf{a} \times \mathbf{b})_l$. ∎

**Theorem 2.2.2** (Dot product from trace)**.** *$\mathbf{a} \cdot \mathbf{b} = \frac{1}{2}\mathrm{tr}(AB)$.*

*Proof.* $\mathrm{tr}(AB) = \sum_{jk} a_j b_k \mathrm{tr}(\sigma_j \sigma_k) = \sum_{jk} a_j b_k \cdot 2\delta_{jk} = 2\sum_j a_j b_j$. ∎

**Theorem 2.2.3** (Norm from determinant)**.** *$\|\mathbf{v}\|^2 = -\det(V)$.*

*Proof.* $\det(V) = \det\left(\begin{smallmatrix} v_3 & v_1-iv_2 \\ v_1+iv_2 & -v_3 \end{smallmatrix}\right) = -v_3^2 - (v_1^2 + v_2^2) = -(v_1^2 + v_2^2 + v_3^2)$. ∎

**Theorem 2.2.4** (Conjugation preserves structure)**.** *If $V$ is traceless Hermitian and $M \in SU(2)$, then $V' = MVM^\dagger$ is traceless Hermitian, and $\det(V') = \det(V)$. That is, conjugation preserves both vector structure and norm.*

*Proof.* Hermitian: $(MVM^\dagger)^\dagger = (M^\dagger)^\dagger V^\dagger M^\dagger = MVM^\dagger$. Traceless: $\mathrm{tr}(MVM^\dagger) = \mathrm{tr}(V M^\dagger M) = \mathrm{tr}(V)= 0$. Determinant: $\det(MVM^\dagger) = \det(M)\det(V)\det(M^\dagger) = 1 \cdot \det(V) \cdot 1 = \det(V)$. ∎

### Composition and Inversion

Rotations compose by $SU(2)$ multiplication and invert by the adjoint — both inherited from the group structure.

**Theorem 2.2.5** (Composition)**.** *If $M_{AB}$ rotates frame B to frame A, and $M_{BC}$ rotates frame C to frame B, then the composed rotation is $M_{AC} = M_{AB} \cdot M_{BC}$.*

*Proof.* $V_A = M_{AB}\,V_B\,M_{AB}^\dagger = M_{AB}(M_{BC}\,V_C\,M_{BC}^\dagger)M_{AB}^\dagger = (M_{AB}M_{BC})\,V_C\,(M_{AB}M_{BC})^\dagger$. ∎

**Theorem 2.2.6** (Inverse)**.** *$M^{-1} = M^\dagger$.*

*Proof.* $M^\dagger M = I$ by definition of $SU(2)$. ∎

**Lemma 2.2.2** (Double cover)**.** *$M$ and $-M$ produce the same rotation: $MVM^\dagger = (-M)V(-M)^\dagger$ for all traceless Hermitian $V$.*

*Proof.* $(-M)V(-M)^\dagger = (-1)^2 MVM^\dagger = MVM^\dagger$. ∎

**Theorem 2.2.7** (General axis-angle rotation)**.** *A rotation by angle $\theta$ about unit axis $\hat{\mathbf{e}} = (e_1, e_2, e_3)$ is:*

$$M(\theta, \hat{\mathbf{e}}) = \cos\frac{\theta}{2}\,I + i\sin\frac{\theta}{2}(e_1\sigma_1 + e_2\sigma_2 + e_3\sigma_3) \tag{2.12}$$

*Proof.* For axis-aligned $\hat{\mathbf{e}} = \hat{\mathbf{e}}_k$, this reduces to Definition 2.2.7. For general $\hat{\mathbf{e}}$: $(\hat{\mathbf{e}} \cdot \boldsymbol{\sigma})^2 = I$ (from Lemma 2.2.1 and $\|\hat{\mathbf{e}}\| = 1$), so the matrix exponential $\exp(i\frac{\theta}{2}\hat{\mathbf{e}} \cdot \boldsymbol{\sigma}) = \cos\frac{\theta}{2}\,I + i\sin\frac{\theta}{2}\,\hat{\mathbf{e}} \cdot \boldsymbol{\sigma}$. ∎

---

## §2.3 The Adjoint Bridge: $M \to R(M)$

The conjugation $V' = MVM^\dagger$ is linear in $V$ for fixed $M$. Since $V$ encodes a 3-vector via the Pauli basis, the linear map $V \mapsto MVM^\dagger$ can be expressed as a $3 \times 3$ real matrix acting on the coefficient vector $(v_1, v_2, v_3)$. This matrix is the rotation matrix $R(M)$ — the adjoint representation of $SU(2)$.

This section derives $R(M)$ explicitly and establishes it as the bridge from the $SU(2)$ construction layer to the matrix action layer.

**Definition 2.3.1** (Adjoint rotation matrix)**.** *For $M \in SU(2)$, the $3 \times 3$ rotation matrix $R(M)$ is defined by:*

$$R_{jk}(M) = \frac{1}{2}\mathrm{tr}(\sigma_j\,M\,\sigma_k\,M^\dagger) \tag{2.13}$$

*Equivalently: the $(j,k)$ entry of $R$ is the $j$-th Pauli coefficient of $M\sigma_k M^\dagger$ — that is, the $j$-th component of the vector obtained by rotating $\hat{\mathbf{e}}_k$.*

**Theorem 2.3.1** (Adjoint representation equivalence)**.** *For any vector $\mathbf{v}$ encoded as $V$:*

$$MVM^\dagger = \sum_j \left(\sum_k R_{jk}(M)\,v_k\right) \sigma_j \tag{2.14}$$

*That is, the Pauli coefficients of the rotated matrix are $(R(M)\,\mathbf{v})_j$. The conjugation is equivalent to left-multiplication by $R(M)$:*

$$\mathbf{v}' = R(M)\,\mathbf{v} \tag{2.15}$$

*Proof.* Expand $MVM^\dagger = M(\sum_k v_k \sigma_k)M^\dagger = \sum_k v_k\,(M\sigma_k M^\dagger)$. Each $M\sigma_k M^\dagger$ is traceless Hermitian (Theorem 2.2.4), so $M\sigma_k M^\dagger = \sum_j R_{jk}\,\sigma_j$ where $R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j\,M\sigma_k M^\dagger)$ (by orthogonality of the Pauli basis: $\frac{1}{2}\mathrm{tr}(\sigma_j \sigma_l) = \delta_{jl}$). Therefore $MVM^\dagger = \sum_{jk} R_{jk}\,v_k\,\sigma_j$, and the coefficient of $\sigma_j$ is $\sum_k R_{jk}\,v_k = (R\mathbf{v})_j$. ∎

**Theorem 2.3.2** ($R(M)$ is orthogonal with determinant $+1$)**.** *$R(M) \in SO(3)$.*

*Proof.* Orthogonality: $\sum_j R_{jk} R_{jl} = \frac{1}{4}\sum_j \mathrm{tr}(\sigma_j M\sigma_k M^\dagger)\,\mathrm{tr}(\sigma_j M\sigma_l M^\dagger)$. Using the completeness relation $\sum_j (\sigma_j)_{ab}(\sigma_j)_{cd} = 2\delta_{ad}\delta_{bc} - \delta_{ab}\delta_{cd}$ (restricted to traceless part), this reduces to $\delta_{kl}$ after applying unitarity of $M$. Determinant: $\det R = +1$ because $R(I) = I_3$ (the identity rotation) and $R$ varies continuously from $I_3$ as $M$ varies continuously from $I$, so $\det R$ cannot jump from $+1$ to $-1$. ∎

**Theorem 2.3.3** (Explicit $R(M)$ in terms of $\alpha, \beta$)**.** *For $M = \left(\begin{smallmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{smallmatrix}\right)$:*

$$R(M) = \begin{pmatrix} \mathrm{Re}(\alpha^2 - \beta^2) & -\mathrm{Im}(\alpha^2 + \beta^2) & 2\,\mathrm{Im}(\alpha\bar{\beta}) \\ \mathrm{Im}(\alpha^2 - \beta^2) & \mathrm{Re}(\alpha^2 + \beta^2) & -2\,\mathrm{Re}(\alpha\bar{\beta}) \\ -2\,\mathrm{Im}(\alpha\beta) & 2\,\mathrm{Re}(\alpha\beta) & |\alpha|^2 - |\beta|^2 \end{pmatrix} \tag{2.16}$$

*Proof.* Compute $R_{jk} = \frac{1}{2}\mathrm{tr}(\sigma_j M \sigma_k M^\dagger)$ for each $(j,k)$ by direct matrix multiplication. For example, $R_{33} = \frac{1}{2}\mathrm{tr}(\sigma_3 M \sigma_3 M^\dagger)$. Since $\sigma_3 M = \left(\begin{smallmatrix} \alpha & -\bar{\beta} \\ -\beta & -\bar{\alpha} \end{smallmatrix}\right)$ and $\sigma_3 M \sigma_3 = \left(\begin{smallmatrix} \alpha & \bar{\beta} \\ -\beta & \bar{\alpha} \end{smallmatrix}\right)$, we get $\sigma_3 M \sigma_3 M^\dagger = \left(\begin{smallmatrix} |\alpha|^2+|\beta|^2 & \cdots \\ \cdots & \cdots \end{smallmatrix}\right)$; the trace gives $R_{33} = |\alpha|^2 - |\beta|^2$. The remaining entries follow analogously. ∎

**Theorem 2.3.4** (Composition respects the adjoint)**.** *$R(M_1 M_2) = R(M_1)\,R(M_2)$.*

*Proof.* $(M_1 M_2)\,V\,(M_1 M_2)^\dagger = M_1(M_2 V M_2^\dagger)M_1^\dagger$. The inner conjugation applies $R(M_2)$ to $\mathbf{v}$; the outer applies $R(M_1)$ to the result. ∎

**Remark.** Theorem 2.3.4 is the reason the two-layer architecture works: composing in SU(2) (Layer 1) and then converting once gives the same result as converting each factor and composing in SO(3) (Layer 2). The conversion $M \to R(M)$ is a group homomorphism.

### Conversion Cost and Error

The conversion from $M$ to $R(M)$ involves 9 quadratic expressions in $(\alpha, \beta)$ — 18 real multiplications and 9 additions. This is done once per transform, not per vector.

**Theorem 2.3.5** (Conversion precision)**.** *The error in $R(M)$ from evaluating (2.16) in finite arithmetic satisfies:*

$$\|\delta(R)\|_{\max} \leq 4\,\epsilon_{\mathrm{mach}} \tag{2.17}$$

*per entry, where $\epsilon_{\mathrm{mach}}$ is the machine epsilon of the arithmetic type. For binary64: $\sim 4 \times 10^{-16}$.* [P.2.1]

*Proof.* Each entry of $R$ is a sum of products of $\alpha, \beta, \bar{\alpha}, \bar{\beta}$, each with $|\alpha|, |\beta| \leq 1$. Each product introduces $\leq \epsilon_{\mathrm{mach}}$ rounding; each sum of two products introduces $\leq 2\epsilon_{\mathrm{mach}}$. The most complex entry involves two products and a subtraction: $\leq 4\epsilon_{\mathrm{mach}}$. ∎

**Theorem 2.3.6** (Angle sensitivity through the adjoint bridge)**.** *For an elementary rotation $M_k(\alpha)$, the sensitivity of $R(M_k)$ to the input angle $\alpha$ satisfies:*

$$\left\|\frac{dR}{d\alpha}\right\| = 1 \tag{2.18}$$

*Therefore $\|\delta(R)\| \leq |\delta(\alpha)|$: the rotation matrix error equals the angle error.*

*Proof.* $\frac{dR_{jk}}{d\alpha} = \frac{1}{2}\mathrm{tr}(\sigma_j \frac{dM}{d\alpha} \sigma_k M^\dagger + \sigma_j M \sigma_k \frac{dM^\dagger}{d\alpha})$. Since $\|dM/d\alpha\| = \frac{1}{2}$ (the eigenvalues of $dM_k/d\alpha$ have modulus $\frac{1}{2}$, as each entry of $M_k$ is $\cos(\alpha/2)$ or $\sin(\alpha/2)$ with derivative $\frac{1}{2}$), and the trace formula involves two such terms, the total sensitivity is $2 \times \frac{1}{2} = 1$. ∎

**Remark.** Compare with the direct SU(2) sensitivity (the current Ch 2, Theorem 2.2.12): $\|\delta(M)\| \leq \frac{1}{2}|\delta(\alpha)|$. The factor-of-2 difference arises because $R$ is quadratic in $M$ ($R$ involves $M \cdot M^\dagger$). The operational error — the error in the rotated vector — is $\|\delta(\mathbf{v}')\| \leq \|\mathbf{v}\| \cdot \|\delta(R)\| \leq \|\mathbf{v}\| \cdot |\delta(\alpha)|$, which is the same either way: the factor of $\frac{1}{2}$ in $\delta(M)$ is squared by the two-sided action, recovering the factor of 1 in $\delta(R)$.

---

## §2.4 The Homogeneous Lift: Translation as a Linear Operation

A coordinate transform may include a translation (displacement by a constant vector) in addition to a rotation. In the $SU(2)$ layer, translation is incorporated via the dual extension (§2.8), but the action on vectors is affine: $\mathbf{r}' = R(M)\mathbf{r} + \mathbf{t}$. The homogeneous lift converts this affine action into a linear one by embedding the 3-vector into a 4-dimensional space.

**Definition 2.4.1** (Homogeneous matrix)**.** *A rigid-body transform (rotation $R \in SO(3)$ and translation $\mathbf{t} \in \mathbb{R}^3$) is represented by the $4 \times 4$ matrix:*

$$H = \begin{pmatrix} R & \mathbf{t} \\ \mathbf{0}^T & 1 \end{pmatrix} \tag{2.19}$$

*The action on a 4-vector $(\mathbf{x}, w)^T$ is:*

$$H \begin{pmatrix} \mathbf{x} \\ w \end{pmatrix} = \begin{pmatrix} R\mathbf{x} + w\,\mathbf{t} \\ w \end{pmatrix} \tag{2.20}$$

*The fourth coordinate $w$ controls whether translation is applied. A position (a point in space) carries $w = 1$; a displacement (the difference between two points) carries $w = 0$.*

**Remark** (Why $w = 0$ for displacements)**.** In projective geometry, $w = 0$ represents a "point at infinity" — a direction without a location. This projective interpretation is not what we intend here. In the affine interpretation used for rigid-body transforms, $w$ is a tag that distinguishes two types of geometric objects:

- A **point** $(\mathbf{r}, 1)$ has a location in space. Translation moves it.
- A **displacement** $(\boldsymbol{\delta}, 0)$ is the difference between two points. It has direction and magnitude but no location. Translation does not affect it.

The distinction arises naturally: if two points $(\tilde{\mathbf{r}}, 1)$ and $(\mathbf{r}, 1)$ are both translated by $\mathbf{t}$, their difference is unchanged. An error — the gap between the true and computed position — is such a difference.

The following example demonstrates this concretely.

**Example 2.4.1** (Worked: rotation + translation acting on state and error)**.** *Consider a rotation of $90°$ about $\hat{\mathbf{e}}_3$ followed by a translation of $100$ km along $\hat{\mathbf{e}}_1$. The satellite is at $\mathbf{r} = (7000, 0, 0)$ km with a position error of $\delta\mathbf{r} = (0.5, 0.3, 0.1)$ km. The homogeneous matrix is:*

$$H = \begin{pmatrix} 0 & -1 & 0 & 100 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

*State row ($w = 1$):*

$$H \begin{pmatrix} 7000 \\ 0 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 + 100 \\ 7000 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 100 \\ 7000 \\ 0 \\ 1 \end{pmatrix}$$

*The position is rotated ($x \to y$) and translated ($+100$ along $x$). ✓*

*Error row ($w = 0$):*

$$H \begin{pmatrix} 0.5 \\ 0.3 \\ 0.1 \\ 0 \end{pmatrix} = \begin{pmatrix} -0.3 \\ 0.5 \\ 0.1 \\ 0 \end{pmatrix}$$

*The error is rotated ($\delta x \to \delta y$, $\delta y \to -\delta x$) but NOT translated. The translation column contributes $100 \times 0 = 0$. ✓*

*Verification:* the true position is $\tilde{\mathbf{r}} = \mathbf{r} + \delta\mathbf{r} = (7000.5, 0.3, 0.1)$. Transforming:

$$H \begin{pmatrix} 7000.5 \\ 0.3 \\ 0.1 \\ 1 \end{pmatrix} = \begin{pmatrix} -0.3 + 100 \\ 7000.5 \\ 0.1 \\ 1 \end{pmatrix} = \begin{pmatrix} 99.7 \\ 7000.5 \\ 0.1 \\ 1 \end{pmatrix}$$

*The difference:* $(99.7, 7000.5, 0.1, 1) - (100, 7000, 0, 1) = (-0.3, 0.5, 0.1, 0)$. *This matches the error row output exactly.* ∎

**Remark** (Error bounds vs. signed perturbations)**.** The example above uses a signed perturbation vector $\delta\mathbf{r} = (0.5, 0.3, 0.1)$ — the actual displacement between true and computed positions. This signed vector transforms correctly under $H$.

However, the three-error tracking framework (Ch 1) uses non-negative bounds: $\sigma_m(r_1) \geq 0$, $\delta_p(r_2) \geq 0$, etc. These represent the maximum possible error in each component, without a known sign. Applying the rotation $R$ to a vector of non-negative bounds can produce negative components (as seen above: $\delta x' = -0.3$), which are meaningless as bounds.

For non-negative component-wise bounds, the correct propagation is $\delta(r'_j) \leq \sum_k |R_{jk}|\,\delta(r_k)$ — using the absolute-value matrix $|R|$ instead of $R$. This is not the same matrix as the state transform.

For scalar error norms $\delta_{\mathrm{pos}} = \|\boldsymbol{\delta}(\mathbf{r})\|$, the propagation is trivial: $\|R\,\boldsymbol{\delta}\| = \|\boldsymbol{\delta}\|$. Rotation preserves the norm exactly, and no matrix multiplication is needed for the error. In Example 2.4.1: $\|\delta\mathbf{r}\| = \sqrt{0.25 + 0.09 + 0.01} = \sqrt{0.35} \approx 0.592$ km before and after the rotation.

The choice between these representations — signed perturbation vectors, component-wise bounds, or scalar norms — is deferred to the final error architecture decision. The homogeneous framework supports all three; the scalar norm is the simplest and avoids the sign issue entirely.

**Theorem 2.4.1** (Homogeneous composition)**.** *$H_1 H_2 = \left(\begin{smallmatrix} R_1 R_2 & R_1\mathbf{t}_2 + \mathbf{t}_1 \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$.*

*Proof.* Block matrix multiplication. ∎

**Theorem 2.4.2** (Homogeneous inverse)**.** *$H^{-1} = \left(\begin{smallmatrix} R^T & -R^T\mathbf{t} \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$.*

*Proof.* Verify $H H^{-1} = I_4$: $R R^T = I_3$, $R(-R^T\mathbf{t}) + \mathbf{t} = -\mathbf{t} + \mathbf{t} = \mathbf{0}$. ∎

**Theorem 2.4.3** (Error propagation through homogeneous transforms)**.** *For a homogeneous transform $H$ with exact entries ($\delta(H) = 0$):*

*(a) Signed perturbation vectors: $H \cdot (\boldsymbol{\delta}, 0)^T = (R\,\boldsymbol{\delta}, 0)^T$. The same matrix $H$ acts on both state and error.*

*(b) Component-wise non-negative bounds: $\delta(r'_j) \leq \sum_k |R_{jk}|\,\delta(r_k)$. The absolute-value matrix $|R|$ replaces $R$.*

*(c) Scalar norms: $\|\boldsymbol{\delta}(\mathbf{r}')\| = \|\boldsymbol{\delta}(\mathbf{r})\|$. Rotation preserves the norm; no matrix needed.*

*In all three cases, translation does not contribute to the error.*

*Proof.* (a): From (2.17) with $w = 0$. (b): From Ch 1, Theorem 1.3.2 applied entry-by-entry to the matrix-vector product, noting that each output component is a sum of products $R_{jk}\,\delta(r_k)$ and the bound requires absolute values. (c): $\|R\,\boldsymbol{\delta}\| = \|\boldsymbol{\delta}\|$ because $R$ is orthogonal. ∎

When $H$ carries error ($\delta(H) \neq 0$), the full error bound acquires additional terms:

$$\|\boldsymbol{\delta}(\mathbf{r}')\| \leq \|\boldsymbol{\delta}(\mathbf{r})\| + \|\mathbf{r}\| \cdot \|\delta(R)\| + \|\delta(\mathbf{t})\| \tag{2.21}$$

by Ch 1, Theorem 1.4.2. Since $R$ is orthogonal, $\|R\| = 1$ and the first term simplifies. The three contributions are: input error (preserved by rotation), rotation error (proportional to position magnitude), and translation error.

---

## §2.5 The $7 \times 7$ State Matrix: Incorporating Velocity

A coordinate transform between frames that rotate relative to each other affects not only the position but also the velocity. The velocity acquires a **transport term** — the contribution from the rotation of the frame itself. The $7 \times 7$ state matrix captures this coupling.

**Definition 2.5.1** (Skew-symmetric matrix)**.** *For $\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3)$:*

$$[\boldsymbol{\omega}]_\times = \begin{pmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{pmatrix} \tag{2.22}$$

*satisfying $[\boldsymbol{\omega}]_\times \mathbf{v} = \boldsymbol{\omega} \times \mathbf{v}$ for all $\mathbf{v}$.*

**Definition 2.5.2** (The $7 \times 7$ state matrix)**.** *A rotating-frame transform with rotation $R$, angular velocity $\boldsymbol{\omega}_T$ of the target frame, translation $\mathbf{t}$, and velocity offset $\Delta\mathbf{v}$ is:*

$$T_7 = \begin{pmatrix} R & \mathbf{0} & \mathbf{t} \\ [\boldsymbol{\omega}_T]_\times R & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \tag{2.23}$$

*acting on the state vector $(\mathbf{r}, \mathbf{v}, 1)^T$:*

$$\begin{pmatrix} \mathbf{r}' \\ \mathbf{v}' \\ 1 \end{pmatrix} = T_7 \begin{pmatrix} \mathbf{r} \\ \mathbf{v} \\ 1 \end{pmatrix} = \begin{pmatrix} R\mathbf{r} + \mathbf{t} \\ [\boldsymbol{\omega}_T]_\times R\mathbf{r} + R\mathbf{v} + \Delta\mathbf{v} \\ 1 \end{pmatrix} \tag{2.24}$$

The velocity transform has three contributions:
- $R\mathbf{v}$: the velocity rotated into the new frame
- $[\boldsymbol{\omega}_T]_\times R\mathbf{r}$: the transport term — the velocity arising from the rotation of the frame at the satellite's position
- $\Delta\mathbf{v}$: a velocity offset (zero for standard orbit propagation transforms)

**Theorem 2.5.1** (Transport theorem derivation)**.** *Let $\mathbf{r}(t)$ be a position vector expressed in a rotating frame with angular velocity $\boldsymbol{\omega}_T$. The velocity in the inertial frame is:*

$$\mathbf{v}_{\mathrm{inertial}} = R\,\mathbf{v}_{\mathrm{rotating}} + \boldsymbol{\omega}_T \times \mathbf{r}_{\mathrm{inertial}} \tag{2.25}$$

*where $\mathbf{r}_{\mathrm{inertial}} = R\,\mathbf{r}_{\mathrm{rotating}} + \mathbf{t}$.*

*Proof.* The inertial position is $\mathbf{r}_I(t) = R(t)\,\mathbf{r}_R(t) + \mathbf{t}(t)$. Differentiating: $\dot{\mathbf{r}}_I = \dot{R}\,\mathbf{r}_R + R\,\dot{\mathbf{r}}_R + \dot{\mathbf{t}}$. The angular velocity satisfies $\dot{R} = [\boldsymbol{\omega}_T]_\times R$. Substituting: $\mathbf{v}_I = [\boldsymbol{\omega}_T]_\times R\,\mathbf{r}_R + R\,\mathbf{v}_R + \dot{\mathbf{t}} = \boldsymbol{\omega}_T \times \mathbf{r}_I + R\,\mathbf{v}_R + (\dot{\mathbf{t}} - \boldsymbol{\omega}_T \times \mathbf{t})$. For static translation ($\dot{\mathbf{t}} = 0$) and $\mathbf{t} = 0$ (the standard orbit propagation case), this reduces to (2.24). ∎

**Remark.** In the $SU(2)$ framework, the transport term $\boldsymbol{\omega}_T \times \mathbf{r}$ corresponds to the commutator $\frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}]$ (Theorem 2.2.1). The 7×7 matrix expresses this as the matrix product $[\boldsymbol{\omega}_T]_\times R\,\mathbf{r}$ — a left-multiplication, consistent with the action layer's linear framework.

**Theorem 2.5.2** ($7 \times 7$ composition)**.** *$T_{7,1} \cdot T_{7,2}$ is a $7 \times 7$ state matrix of the same form, with:*

$$R_{12} = R_1 R_2, \qquad \boldsymbol{\omega}_{T,12} = \boldsymbol{\omega}_{T,1} + R_1\,\boldsymbol{\omega}_{T,2}$$

$$\mathbf{t}_{12} = R_1\mathbf{t}_2 + \mathbf{t}_1, \qquad \Delta\mathbf{v}_{12} = [\boldsymbol{\omega}_{T,1}]_\times R_1\mathbf{t}_2 + R_1\Delta\mathbf{v}_2 + \Delta\mathbf{v}_1$$

*Proof.* Block matrix multiplication of two $7 \times 7$ matrices in the form (2.22). The angular velocity composition follows from the lower-left block: $[\boldsymbol{\omega}_{T,1}]_\times R_1 R_2 + R_1 [\boldsymbol{\omega}_{T,2}]_\times R_2 = ([\boldsymbol{\omega}_{T,1}]_\times + R_1 [\boldsymbol{\omega}_{T,2}]_\times R_1^T) R_1 R_2 = [\boldsymbol{\omega}_{T,1} + R_1\boldsymbol{\omega}_{T,2}]_\times R_{12}$, using the identity $R[\boldsymbol{\omega}]_\times R^T = [R\boldsymbol{\omega}]_\times$. ∎

**Corollary 2.5.1** (Error linearity of the $7 \times 7$ matrix)**.** *For exact $T_7$ ($\delta(T_7) = 0$), the error propagation is:*

$$T_7 \begin{pmatrix} \boldsymbol{\delta}(\mathbf{r}) \\ \boldsymbol{\delta}(\mathbf{v}) \\ 0 \end{pmatrix} = \begin{pmatrix} R\,\boldsymbol{\delta}(\mathbf{r}) \\ [\boldsymbol{\omega}_T]_\times R\,\boldsymbol{\delta}(\mathbf{r}) + R\,\boldsymbol{\delta}(\mathbf{v}) \\ 0 \end{pmatrix} \tag{2.26}$$

*The same $T_7$ acts on both the state ($w = 1$) and the error ($w = 0$). Translation and velocity offsets are killed by $w = 0$. Position error couples into velocity error through the transport term — this is physically correct (a position uncertainty of $\delta r$ in a rotating frame produces a velocity uncertainty of $|\omega_T| \cdot \delta r$).*

---

## §2.6 The Orbit Propagation Coordinate Transforms

This section applies the framework to the specific transforms of the orbit propagation pipeline.

### Perifocal → TEME

**Definition 2.6.1** (Perifocal frame)**.** *The frame aligned with the orbit: $\hat{\mathbf{e}}_1$ toward perigee, $\hat{\mathbf{e}}_3$ along the angular momentum vector, $\hat{\mathbf{e}}_2$ completing the right-handed system.*

**Definition 2.6.2** (TEME frame)**.** *True Equator, Mean Equinox — the reference frame used by the propagator (Ch 30). $\hat{\mathbf{z}}$ along Earth's rotation axis, $\hat{\mathbf{x}}$ toward the mean vernal equinox.* [A.2.1]

The perifocal-to-TEME transform is a pure rotation (no translation, no angular velocity): the perifocal frame is not rotating relative to TEME at the instant of evaluation.

**Theorem 2.6.1** (Perifocal → TEME rotation)**.** *The rotation from perifocal to TEME is the 3-1-3 Euler sequence:*

$$M_{\mathrm{PF} \to \mathrm{TEME}} = M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u) \tag{2.27}$$

*where $\Omega$ is the right ascension of the ascending node, $i$ is the inclination, and $u = \omega + \nu$ is the argument of latitude. The negative signs implement passive rotation (frame rotation). The rotation matrix is $R_{\mathrm{PF} \to \mathrm{TEME}} = R(M_{\mathrm{PF} \to \mathrm{TEME}})$ via (2.10). The state matrix is:*

$$T_{7,\mathrm{PF} \to \mathrm{TEME}} = \begin{pmatrix} R & \mathbf{0} & \mathbf{0} \\ \mathbf{0} & R & \mathbf{0} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \tag{2.28}$$

*with $[\boldsymbol{\omega}_T]_\times = \mathbf{0}$ (no transport term) and $\mathbf{t} = \mathbf{0}$, $\Delta\mathbf{v} = \mathbf{0}$.*

### TEME → PEF

The TEME-to-PEF (Pseudo Earth-Fixed) transform accounts for Earth's rotation. This is a rotating-frame transform with a transport term.

**Theorem 2.6.2** (TEME → PEF state matrix)**.** *The transform is:*

$$T_{7,\mathrm{TEME} \to \mathrm{PEF}} = \begin{pmatrix} R_3(\theta) & \mathbf{0} & \mathbf{0} \\ [\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta) & R_3(\theta) & \mathbf{0} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{pmatrix} \tag{2.29}$$

*where $\theta = \theta_{\mathrm{GMST}}$ is the Greenwich Mean Sidereal Time (Ch 29), $R_3(\theta) = R(M_3(\theta))$, and $\omega_E$ is the Earth rotation rate. The transport term $[\omega_E\hat{\mathbf{e}}_3]_\times R_3(\theta)\,\mathbf{r}$ produces the velocity correction $\boldsymbol{\omega}_E \times \mathbf{r}_{\mathrm{PEF}}$.* [M.2.1]

### Pipeline Composition

**Theorem 2.6.3** (Full pipeline)**.** *The complete transform from perifocal to PEF is:*

$$T_{7,\mathrm{PF} \to \mathrm{PEF}} = T_{7,\mathrm{TEME} \to \mathrm{PEF}} \cdot T_{7,\mathrm{PF} \to \mathrm{TEME}} \tag{2.30}$$

*This is a single $7 \times 7$ matrix multiply. The rotation part is $R_3(\theta) \cdot R_{\mathrm{PF \to TEME}}$, which can be computed in Layer 1 as $M_3(\theta) \cdot M_3(-\Omega) \cdot M_1(-i) \cdot M_3(-u)$ and converted once via $R(M)$.*

---

## §2.7 Error Propagation: The Parallel-Row Principle

The homogeneous structure of the $7 \times 7$ state matrix enables a simple and powerful error architecture: the same matrix that propagates the physical state also propagates the error state.

### The Principle

**Theorem 2.7.1** (Parallel-row error propagation)**.** *Let $T_7$ be a state matrix with exact entries. Define:*

- *Row 1 (state): $\mathbf{s} = (\mathbf{r}, \mathbf{v}, 1)^T$*
- *Row 2 (error): $\mathbf{e} = (\boldsymbol{\delta}(\mathbf{r}), \boldsymbol{\delta}(\mathbf{v}), 0)^T$*

*Then $T_7 \mathbf{s}$ produces the transformed state and $T_7 \mathbf{e}$ produces the transformed error, using the same matrix and the same multiplication. No separate error formulas are needed for linear transforms.*

*Proof.* Immediate from (2.23) and (2.25): the $w = 1$ vs $w = 0$ coordinate controls whether translation/velocity offsets contribute. Both undergo the same rotation and transport coupling. ∎

This principle extends to the three error categories. Since $\sigma_m$, $\delta_p$, and $\delta_a$ propagate by the same rules (Ch 1, Proposition 1.8.1), each category can be carried as a separate error row:

| Row | Content | Homogeneous $w$ |
|-----|---------|----------------|
| 1 | Physical state $(\mathbf{r}, \mathbf{v})$ | 1 |
| 2 | Measurement error $(\boldsymbol{\sigma}_m(\mathbf{r}), \boldsymbol{\sigma}_m(\mathbf{v}))$ | 0 |
| 3 | Precision error $(\boldsymbol{\delta}_p(\mathbf{r}), \boldsymbol{\delta}_p(\mathbf{v}))$ | 0 |
| 4 | Accuracy error $(\boldsymbol{\delta}_a(\mathbf{r}), \boldsymbol{\delta}_a(\mathbf{v}))$ | 0 |

All four rows are transformed by the same $T_7$. The computation is performed once and applied four times.

### Inexact Transforms

When $T_7$ carries error (from computing $R(M)$ with finite-precision trigonometric functions, or from uncertain input angles), the error rows receive an additional contribution.

**Theorem 2.7.2** (Inexact transform error)**.** *When $T_7$ has entry-wise error $\delta(T_7)$, the output error for the state row is:*

$$\boldsymbol{\delta}(\mathbf{s}') \leq T_7\,\mathbf{e} + |T_7'|\,|\mathbf{s}| \tag{2.31}$$

*where $|T_7'|$ is the matrix of absolute values of the transform error, and $|\mathbf{s}|$ is the component-wise absolute value of the state. The first term is the propagated input error (same matrix); the second is the error introduced by the inexact transform.*

*Proof.* By Ch 1, Theorem 1.4.2 (multivariate bound) applied to the matrix-vector product $T_7 \mathbf{s}$: the error from perturbing $T_7$ by $\delta(T_7)$ while holding $\mathbf{s}$ fixed is $|\delta(T_7)| \cdot |\mathbf{s}|$, and the error from perturbing $\mathbf{s}$ while holding $T_7$ fixed is $|T_7| \cdot |\boldsymbol{\delta}(\mathbf{s})|$. The bound (2.30) follows. ∎

### Pipeline Classification

**Theorem 2.7.3** (Linear vs. nonlinear operations)**.** *Every operation in the orbit propagation pipeline is classified:*

| Operation | Type | Error handling |
|-----------|------|---------------|
| Rotation of state vectors | Linear | Same $T_7$ on error row (Theorem 2.7.1) |
| Translation | Linear | Same $T_7$; translation killed by $w=0$ (Theorem 2.4.3) |
| Transport term | Linear | Same $T_7$; couples $\delta(\mathbf{r})$ into $\delta(\mathbf{v})$ (Corollary 2.5.1) |
| $SU(2)$ composition | Linear | $R(M_1 M_2) = R(M_1) R(M_2)$ (Theorem 2.3.4) |
| $M \to R(M)$ conversion | Nonlinear | Bounded by Theorem 2.3.5; $\delta(R) \leq 4\epsilon_{\mathrm{mach}}$ |
| sin/cos of input angles | Nonlinear | Bounded by Ch 1, Corollaries 1.4.1–1.4.2 |
| Kepler solver | Nonlinear | Bounded by Ch 1, Theorem 1.6.1; sensitivity derived in Ch 9 |
| $SU(2)$ renormalization | Nonlinear | Bounded by Ch 1, Corollary 1.4.4 (sqrt); rarely needed |

*The linear operations (the majority of the pipeline) use the parallel-row principle directly. The nonlinear operations are concentrated at the interfaces — constructing transforms from angles and solving Kepler's equation — and each has a rigorous sensitivity bound derived in the chapter where that operation is developed.*

---

## §2.8 The Dual Extension: Composition with Translation

The homogeneous matrix (§2.4) handles translation in the action layer. The dual extension provides an alternative mechanism in the construction layer, enabling rotation + translation to be composed by a single algebraic product — without ever forming the $3 \times 3$ matrix $R(M)$.

This section develops the dual $SU(2)$ algebra for completeness and for applications where the composition structure is advantageous (e.g., formal pipeline analysis, symbolic manipulation). The operational pipeline (§2.7) uses the homogeneous form.

### Dual Numbers

**Definition 2.8.1** (Dual unit)**.** *The dual unit $\varepsilon$ satisfies $\varepsilon^2 = 0$ and $\varepsilon \neq 0$.*

**Definition 2.8.2** (Dual number)**.** *A dual number is $\hat{a} = a + \varepsilon\,a'$ where $a, a' \in \mathbb{R}$. Arithmetic:*

$$(\hat{a} + \hat{b}) = (a + b) + \varepsilon(a' + b'), \qquad \hat{a}\,\hat{b} = ab + \varepsilon(ab' + a'b) \tag{2.32}$$

**Remark.** Evaluating a polynomial $P$ at a dual number yields $P(a + \varepsilon\,b) = P(a) + \varepsilon\,b\,P'(a)$. This is a consequence of $\varepsilon^2 = 0$ killing all terms beyond first order. The result extends to any function admitting a Taylor series: $f(a + \varepsilon\,b) = f(a) + \varepsilon\,b\,f'(a)$. The dual part carries the derivative scaled by $b$ — a property that connects dual numbers to automatic differentiation and to first-order perturbation theory.

### Dual $SU(2)$

**Definition 2.8.3** (Dual $SU(2)$ element)**.** *A dual $SU(2)$ element is $\hat{M} = M + \varepsilon D$ where $M \in SU(2)$ and $D \in \mathbb{C}^{2 \times 2}$ satisfies:*

1. *$M^\dagger M = I$ (the non-dual part is unitary)*
2. *$M^\dagger D + D^\dagger M = 0$ (the dual part is "tangent" to the unitary constraint)*

*The dual part encodes a translation: $D = \frac{1}{2}TM$ where $T = t_1\sigma_1 + t_2\sigma_2 + t_3\sigma_3$ is the translation vector encoded as a traceless Hermitian matrix.*

**Theorem 2.8.1** (Dual $SU(2)$ composition)**.** *The product of two dual $SU(2)$ elements is:*

$$\hat{M}_1 \hat{M}_2 = M_1 M_2 + \varepsilon(M_1 D_2 + D_1 M_2) \tag{2.33}$$

*The composed rotation is $M_1 M_2$. The composed translation, recovered by $T_{12} = 2(M_1 D_2 + D_1 M_2)(M_1 M_2)^\dagger$, equals $M_1 T_2 M_1^\dagger + T_1$ — rotate the second translation by the first rotation, then add the first translation. This is the natural composition of rigid motions.*

*Proof.* The product follows from (2.31) with matrix entries. The translation extraction: $(M_1 D_2 + D_1 M_2)(M_1 M_2)^\dagger = M_1 D_2 M_2^\dagger M_1^\dagger + D_1 M_1^\dagger$. Using $D_k = \frac{1}{2}T_k M_k$: $D_k M_k^\dagger = \frac{1}{2}T_k$. Therefore $T_{12} = 2(M_1 \frac{1}{2}T_2 M_1^\dagger + \frac{1}{2}T_1) = M_1 T_2 M_1^\dagger + T_1$. ∎

**Theorem 2.8.2** (Equivalence to homogeneous matrix)**.** *The dual $SU(2)$ composition (2.32) produces the same rotation and translation as the $4 \times 4$ homogeneous matrix product (Theorem 2.4.1). Specifically: if $H_k = \left(\begin{smallmatrix} R(M_k) & \mathbf{t}_k \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$, then $H_1 H_2 = \left(\begin{smallmatrix} R(M_1 M_2) & R(M_1)\mathbf{t}_2 + \mathbf{t}_1 \\ \mathbf{0}^T & 1 \end{smallmatrix}\right)$, which is the matrix form of $\hat{M}_1 \hat{M}_2$.*

*Proof.* The rotation $R(M_1 M_2) = R(M_1)R(M_2)$ by Theorem 2.3.4. The translation $R(M_1)\mathbf{t}_2 + \mathbf{t}_1$ corresponds to $M_1 T_2 M_1^\dagger + T_1$ by the adjoint equivalence (Theorem 2.3.1): the conjugation $M_1 T_2 M_1^\dagger$ is the same vector as $R(M_1)\mathbf{t}_2$. ∎

**Remark.** Theorem 2.8.2 confirms that the two layers are consistent: composing in the $SU(2)$ layer (dual quaternion product) and then converting to the matrix layer gives the same result as converting first and composing in the matrix layer (homogeneous matrix product). The choice between them is one of computational convenience, not mathematical content.

### Affine Nature of the Dual Extension

The dual $SU(2)$ element $\hat{M}$ acts on a position encoded as a traceless Hermitian matrix $V$ by:

$$V' = MVM^\dagger + T \tag{2.34}$$

This action is **affine** in $V$: the rotation $MVM^\dagger$ is linear, but the translation $+T$ is a constant shift. For error propagation, this means the error of $V'$ includes $M\,\delta(V)\,M^\dagger + \delta(T)$: the state sees $+T$, the error sees $+\delta(T)$. These are different operations unless $T$ is exact.

By contrast, the homogeneous matrix acts on $(\mathbf{v}, w)^T$ by left-multiplication — truly linear. The $w$ coordinate handles the distinction automatically (§2.4).

This is the reason the operational pipeline (§2.7) uses the homogeneous form rather than the dual extension for state propagation. The dual extension remains valuable for composition (Theorem 2.8.1), formal analysis, and the connection to dual quaternion literature.

---

## §2.9 Summary and Cross-References

The two-layer architecture provides:

1. **Singularity-free construction** via $SU(2)$ — no gimbal lock in building or composing rotations
2. **Truly linear action** via the homogeneous state matrix — the same $T_7$ acts on state and error rows
3. **Rigorous error propagation** — linear operations use the parallel-row principle; nonlinear operations (sin/cos, Kepler solver) have per-operation sensitivity bounds derived in their respective chapters

The bridge between layers is the adjoint representation $M \to R(M)$ (§2.3), which is an exact group homomorphism computed once per transform.

Downstream chapters build on this framework:

- Ch 8: position and velocity from orbital elements — perifocal vectors, then $T_{7,\mathrm{PF} \to \mathrm{TEME}}$
- Ch 9–10: Kepler solver — a nonlinear step with sensitivity bound for the error rows
- Ch 18–20: perturbation corrections — near-identity $T_7$ matrices
- Ch 29: GMST computation — the angle that parameterizes $T_{7,\mathrm{TEME} \to \mathrm{PEF}}$
- Ch 30: complete coordinate transform pipeline — composition of $T_7$ matrices
- Ch 38: final error budget — the product of all $T_7$ matrices gives the end-to-end error

---

## Error Notes

**[P.2.1]** Adjoint conversion rounding. The 9 entries of $R(M)$ are quadratic in $(\alpha, \beta)$. Each entry requires $\leq 4$ floating-point operations on values with $|\alpha|, |\beta| \leq 1$, bounding the rounding error at $\leq 4\epsilon_{\mathrm{mach}}$ per entry. For binary64: $\sim 4 \times 10^{-16}$. For single precision: $\sim 5 \times 10^{-7}$. *Remedy:* use wider arithmetic type.

**[P.2.2]** Transport term precision. The cross product $\boldsymbol{\omega}_T \times \mathbf{r}$ in the transport term is computed as $[\boldsymbol{\omega}_T]_\times R\mathbf{r}$. For Earth rotation ($|\omega_E| \approx 7.3 \times 10^{-5}$ rad/s), a position error of $\delta r$ produces velocity error $|\omega_E|\,\delta r$. *Bound:* $\delta v_{\mathrm{transport}} \leq |\omega_E|\,\delta r + |\mathbf{r}|\,\delta\omega_E$.

**[M.2.1]** GMST angle uncertainty. The GMST polynomial (Ch 29) determines $\theta_{\mathrm{GMST}}$, which parameterizes the TEME→PEF rotation. Measurement error in the polynomial coefficients propagates through Theorem 2.3.6: $\delta(R) \leq |\delta(\theta)|$. *Remedy:* use updated IERS values (but this breaks the matched pair — see Ch 3).

**[A.2.1]** TEME frame accuracy. The TEME frame differs from precise GCRS/ITRS frames by $\sim 0.1$ arcsec due to the truncated precession/nutation model. This is an irreducible accuracy error. *Remedy:* apply IAU 2006/2000A corrections (beyond the standard propagator model; see Ch 3).
