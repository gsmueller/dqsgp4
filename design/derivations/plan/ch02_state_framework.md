# Draft Plan: Chapter 2 — The State Framework

**Part I: Mathematical Foundations**

**Phase:** Build (complete as-built)

---

## Objectives

1. Build the SU(2) algebra from the Pauli basis, establishing vector encoding, conjugation-as-rotation, and the algebraic identities connecting matrix operations to cross product, dot product, and norm.
2. Prove that conjugation by SU(2) is a norm-preserving, structure-preserving rotation, and derive composition and inversion rules.
3. Extend SU(2) to dual SU(2) via dual numbers, absorbing translation into the algebra so that rigid motions compose by a single product.
4. Define the complete rigid body state as two dual quaternions (M̂, Ω̂_b), derive the kinematic equation, and establish the transport term from both the Lie algebra and direct calculus routes.
5. Build the coordinate transform pipeline: perifocal → TEME → PEF, including the 3-1-3 rotation, the Earth-rotation transport term, and the composed end-to-end transform.
6. Classify every pipeline operation as linear (Principle 1) or nonlinear (Principle 2) for error propagation, with proofs.
7. Derive the adjoint bridge M → R(M) as a group homomorphism, give explicit quaternion-to-DCM formulas, and derive the inverse (Shepperd's method).
8. Assemble the 4×4 homogeneous matrix and 7×7 state matrix from the dual quaternion pair, proving equivalence to the SU(2) operations.
9. Define the parallel-row error architecture: four rows of 13 parameters, scalar error reductions, point-mass mode.
10. Summarize how each subsequent chapter uses the framework, with a downstream reference table.

---

## Notation Table

The chapter opens with a full notation table (no section number) covering four groups:

| Group | Symbols |
|-------|---------|
| SU(2) algebra | $\sigma_k$, $\boldsymbol{\sigma}$, $I$, $A^\dagger$, $M$, $\alpha$/$\beta$, $(w,x,y,z)$, $V$, $\delta_{jk}$, $\varepsilon_{jkl}$, $[A,B]$, $M_k(\alpha)$, $\hat{\mathbf{e}}_k$ |
| Dual extension | $\hat{M}$, $\varepsilon$ (dual unit), $D$, $T$ (translation matrix) |
| Rigid body state | $R_{\mathrm{pos}}$, $V_{\mathrm{vel}}$, $\Omega$, $\Omega_T$ |
| Orbital elements / frames | $\Omega_{\mathrm{node}}$, $i$, $u$, $e$, $r$, $\omega_E$, $\theta_{\mathrm{GMST}}$ |
| Matrix forms | $R(M)$, $H$, $w_H$, $T_7$, $[\mathbf{a}]_\times$, $\Delta\mathbf{v}$ |
| Error framework | $\delta(\cdot)$, $\sigma_m$/$\delta_p$/$\delta_a$, $\epsilon_{\mathrm{mach}}$, $\mathrm{rd}(v)$ |

---

## Section Structure

### §2.1 Introduction

Narrative only. Three properties of the framework (singularity-free, compositional, error-compatible). Overview of SU(2) conjugation, the adjoint bridge, and the 7×7 state matrix as a component-wise expansion. No numbered items.

---

### §2.2 The Algebra of Rotation

#### Subsection: Foundations

- **Definition 2.2.1** (Pauli matrices). The three $\sigma_k$ matrices; basis for Hermitian 2×2 matrices.
- **Definition 2.2.2** (Conjugate transpose). $(A^\dagger)_{jk} = \overline{a_{kj}}$; key properties.
- **Definition 2.2.3** (Hermitian matrix). $A = A^\dagger$; 4-dimensional real vector space.
- **Definition 2.2.4** (Traceless Hermitian matrix). $\mathrm{tr}(A) = 0$; 3-dimensional real vector space; basis $\{\sigma_1, \sigma_2, \sigma_3\}$.
- **Definition 2.2.5** (Vector encoding). $V = v_k\sigma_k$; isomorphism $\mathbb{R}^3 \leftrightarrow$ traceless Hermitian $2\times 2$; recovery formula $v_k = \frac{1}{2}\mathrm{tr}(\sigma_k V)$.
- **Example 2.2.1** (Vector encoding). $\mathbf{v} = (3, 4, 0)$; norm via determinant; component recovery. ✓
- **Definition 2.2.6** (Rotation element). $M \in SU(2)$: $M^\dagger M = I$, $\det M = 1$; general form with $\alpha, \beta \in \mathbb{C}$; quaternion correspondence $\alpha = w+iz$, $\beta = ix+y$.
- **Assumption 2.2.1** (Right-handed coordinates). All coordinate systems right-handed throughout.
- **Definition 2.2.7** (Conjugation). $V' = MVM^\dagger$; the rotation action.
- **Definition 2.2.8** (Elementary rotations). $M_k(\alpha) = \cos\frac{\alpha}{2}I + i\sin\frac{\alpha}{2}\sigma_k = \exp(i\frac{\alpha}{2}\sigma_k)$; half-angle; explicit $M_1$, $M_2$, $M_3$.
- **Example 2.2.2** (Elementary rotation matrix). $M_3(90°)$; verification of unitarity and determinant. ✓

#### Subsection: Algebraic Properties

- **Definition 2.2.9** (Kronecker delta and Levi-Civita symbol). $\delta_{jk}$, $\varepsilon_{jkl}$; cross product encoding.
- **Lemma 2.2.1** (Pauli algebra). $\sigma_j\sigma_k = \delta_{jk}I + i\varepsilon_{jkl}\sigma_l$; anticommutator, commutator, trace forms; proved by complete multiplication table.
- **Theorem 2.2.1** (Cross product as commutator). $\mathbf{a}\times\mathbf{b} \leftrightarrow \frac{1}{2i}[A,B]$; proved from Lemma 2.2.1.
- **Theorem 2.2.2** (Dot product from trace). $\mathbf{a}\cdot\mathbf{b} = \frac{1}{2}\mathrm{tr}(AB)$; proved from trace identity.
- **Theorem 2.2.3** (Norm from determinant). $\|\mathbf{v}\|^2 = -\det(V)$; proved by direct computation.
- **Example 2.2.3** (Cross product via commutator). $\hat{\mathbf{e}}_1 \times \hat{\mathbf{e}}_2 = \hat{\mathbf{e}}_3$ via $\frac{1}{2i}[\sigma_1, \sigma_2]$. ✓
- **Example 2.2.4** (Rotation by conjugation). Rotate $(1,0,0)$ by $90°$ about $\hat{\mathbf{e}}_3$; step-by-step conjugation; passive rotation interpretation. ✓

#### Subsection: Rotation Properties

- **Theorem 2.2.4** (Conjugation preserves vector structure). $MVM^\dagger$ is traceless Hermitian; proved from Hermitian and trace properties.
- **Theorem 2.2.5** (Conjugation preserves norm). $\det(MVM^\dagger) = \det(V)$; proved by multiplicativity of det.

#### Subsection: Composition and Inversion

- **Theorem 2.2.6** (Composition of rotations). $M_{AC} = M_{AB}\cdot M_{BC}$; proved by nested conjugation; closure verified.
- **Theorem 2.2.7** (Inverse rotation). $M^{-1} = M^\dagger$; immediate from definition.
- **Lemma 2.2.2** (Double cover). $M$ and $-M$ produce the same rotation; convention for canonical representative.
- **Theorem 2.2.8** (General axis-angle rotation). $M(\theta, \hat{\mathbf{e}}) = \cos\frac{\theta}{2}I + i\sin\frac{\theta}{2}(\hat{\mathbf{e}}\cdot\boldsymbol{\sigma})$; proved using $E^2 = I$ and Euler formula for involutions.

#### Subsection: Extracting the Quaternion from a Rotation Matrix

Narrative only (forward reference to §2.9 Theorem 2.9.1 and §2.7 adjoint bridge).

#### Subsection: Error Properties of Rotation

- **Remark** (Matrix norms from scalar bounds). Operator norm; extension from Ch 1 scalar propagation; $\|M\| = 1$ for unitary $M$.
- **Theorem 2.2.9** (Isometry — exact rotation preserves error magnitude). $\delta(V') = M\,\delta(V)\,M^\dagger$; $\|\delta(\mathbf{v}')\| = \|\delta(\mathbf{v})\|$; exact (no linearization). Linear in $V$ for fixed $M$.
- **Lemma 2.2.3** (Spectral norm of encoded vector). $\|V\| = \|\mathbf{v}\|$ (spectral = Euclidean); eigenvalues $\pm\|\mathbf{v}\|$.
- **Theorem 2.2.10** (Inexact rotation error). $\|\delta(V')\| \leq \|\delta(\mathbf{v})\| + 2\|\mathbf{v}\|\|\delta(M)\|$; first-order bound; two cross-term proof via submultiplicativity. **[P.2.2]**
- **Theorem 2.2.11** (Angle-to-rotation sensitivity). $\|dM_k/d\alpha\| = \frac{1}{2}$; $\|\delta(M_k)\| \leq \frac{1}{2}|\delta(\alpha)|$; proved by eigenvalue of derivative. **[P.2.3]**
- **Corollary 2.2.1** (Composed rotation error). For 3-1-3 composition: $\|\delta(M)\| \leq \frac{1}{2}(|\delta(\Omega)| + |\delta(i)| + |\delta(u)|) + 16\epsilon_{\mathrm{mach}}$. **[P.2.4]**

---

### §2.3 Incorporating Translation: The Dual Extension

#### Subsection: Dual Numbers

- **Definition 2.3.1** (Dual numbers). $\hat{a} = a + \varepsilon a'$, $\varepsilon^2 = 0$; arithmetic rule; first-order perturbation interpretation; analytic function extension.

#### Subsection: Dual SU(2)

- **Definition 2.3.2** (Dual SU(2) matrix). $\hat{M} = M + \varepsilon D$; 8 real parameters.
- **Definition 2.3.3** (Unit dual SU(2) constraints). Two constraints reducing 8 parameters to 6 DOF.
- **Definition 2.3.4** (Translation encoding). $D = \frac{1}{2}TM$ where $T = t_k\sigma_k$; factor of $\frac{1}{2}$ convention.
- **Theorem 2.3.1** (Translation extraction). $T = 2DM^\dagger$; proved by right-multiplying definition by $M^\dagger$.

#### Subsection: Composition and Inverse

- **Theorem 2.3.2** (Dual SU(2) composition). Product rule; composed translation $T_{12} = M_1 T_2 M_1^\dagger + T_1$; natural composition of rigid motions.
- **Theorem 2.3.3** (Dual SU(2) inverse). $\hat{M}^{-1} = M^\dagger - \varepsilon D^\dagger$; proved by direct verification.
- **Theorem 2.3.4** (Equivalence to 4×4 homogeneous matrix). Same rotation and translation as homogeneous block product; equivalence via adjoint representation (§2.7).
- **Example 2.3.1** (Composing rotation + translation). $90°$ about $\hat{\mathbf{e}}_3$ then translate $(100, 0, 0)$ km; dual product; extracted translation; action on $(7000, 0, 0)$; 4×4 verification. ✓

#### Subsection: Linearity of the Dual Action

- **Theorem 2.3.5** (The dual action is affine in the operand). $R_{\mathrm{pos}}' = MR_{\mathrm{pos}}M^\dagger + T$; linear in $R_{\mathrm{pos}}$; error transforms: $\delta(R_{\mathrm{pos}}') = M\,\delta(R_{\mathrm{pos}})\,M^\dagger + \delta(T)$.
- **Remark** (Affine action in the standard pipeline). $T = 0$ for all standard orbit propagation transforms; affine term vanishes; formal property only.

#### Subsection: Position Extraction Error

- **Theorem 2.3.6** (Error in position extraction). $\|\delta(\mathbf{r})\| \leq 2\|\delta(D)\| + \|\mathbf{r}\|\|\delta(M)\|$; bilinear derivation. **[P.2.5]**
- **Remark**. LEO numerical estimate: rotation-coupled position error $\approx 0.035$ m — negligible vs. TLE accuracy.

#### Subsection: Elementary Operations

- **Theorem 2.3.7** (Elementary translation along $\hat{\mathbf{e}}_k$). $\hat{M}_{\mathrm{trans}}(\lambda, k) = I + \varepsilon\frac{\lambda}{2}\sigma_k$; explicit component actions for all three axes.
- **Theorem 2.3.8** (General translation). Superposition of three elementary translations.
- **Remark** (Sign convention). Passive vs. active rotation; negative angles in §2.5 implement passive frame rotations.
- **Theorem 2.3.9** (Elementary rotation about $\hat{\mathbf{e}}_1$). Passive conjugation → standard $R_1(\alpha)$ components; proved by entry-by-entry computation.
- **Theorem 2.3.10** (Elementary rotation about $\hat{\mathbf{e}}_3$). Diagonal $M_3$ → phase multiplication of off-diagonal; $R_3(\alpha)$ components.
- **Theorem 2.3.11** (Elementary rotation about $\hat{\mathbf{e}}_2$). $R_2(\alpha)$ components; same pattern as $\hat{\mathbf{e}}_1$ with permuted components.
- **Theorem 2.3.12** (Composition of translation and rotation). Rotation by $\alpha$ about $\hat{\mathbf{e}}_3$ then translate $\lambda$ along $\hat{\mathbf{e}}_1$; explicit product.
- **Theorem 2.3.13** (Composition of two translations). Two translations → vector addition; identity rotation factors.
- **Theorem 2.3.14** (Composition of two rotations about different axes). $M_3(\beta)\cdot M_1(\alpha)$; explicit $2\times 2$ product; not about any coordinate axis.

---

### §2.4 The Complete Rigid Body State

#### Subsection: State Definition

- **Definition 2.4.1** (Velocity dual quaternion). $\hat{\Omega}_b = \Omega_b + \varepsilon V_b$; Lie algebra $\mathfrak{se}(3)$; 6 real parameters.
- **Definition 2.4.2** (Complete rigid body state). Table: $(\hat{M}, \hat{\Omega}_b)$; Lie group vs. algebra types; 13 real parameters, 12 DOF.
- **Assumption 2.4.1** (Instantaneous angular velocity). Frame $\Omega_T$ treated as constant over one step; Earth rotation justification ($<10^{-8}$ rad/s variation).

#### Subsection: The Kinematic Equation

- **Theorem 2.4.1** (Kinematic equation). $\dot{\hat{M}} = \hat{M}\cdot\frac{1}{2}\hat{\Omega}_b$; proved by expanding the dual product and matching to $\dot{M}$ and $\dot{D}$.

#### Subsection: Spatial and Body Twists

- **Definition 2.4.3** (Spatial twist). $\hat{\Omega}_s = \hat{M}\hat{\Omega}_b\hat{M}^{-1}$; spatial angular velocity as conjugation.
- **Theorem 2.4.2** (Frame transform of the velocity). Configuration: $\hat{M}' = \hat{M}_T\cdot\hat{M}$; spatial twist: $\hat{\Omega}_s' = \hat{M}_T\hat{\Omega}_s\hat{M}_T^{-1} + \hat{\Omega}_{T,s}$.

#### Subsection: Spatial Twist Velocity and Material Velocity

- **Theorem 2.4.3** (Spatial twist expansion). Dual part of $\hat{M}\hat{\Omega}_b\hat{M}^{-1}$ is $V_{\mathrm{vel}} + \frac{1}{2i}[\Omega_s, R_{\mathrm{pos}}]$ — the $\boldsymbol{\omega}_s\times\mathbf{r}$ transport term; proved by expanding adjoint action.
- **Corollary 2.4.1** (Material velocity from spatial twist). $V_{\mathrm{vel}} = V_s - \frac{1}{2i}[\Omega_s, R_{\mathrm{pos}}]$; point-mass reduction $(\Omega_b = 0)$.

#### Subsection: The Transport Term

- **Corollary 2.4.2** (Transport term from frame transform). Under $T_T = 0$, $\Omega_T$ rotation: $V_{\mathrm{vel}}' = M_T V_{\mathrm{vel}} M_T^\dagger - \frac{1}{2i}[\Omega_T, R_{\mathrm{pos}}']$; geostationary sign verification ✓.

#### Subsection: Calculus Verification of the Transport Term

- **Theorem 2.4.4** (Transport from time derivative). $\frac{d}{dt}(M_T R_{\mathrm{pos}} M_T^\dagger) \leftrightarrow \boldsymbol{\omega}_T\times\mathbf{r}'$; proved by product rule on $\dot{M}_T = \frac{i}{2}\Omega_T' M_T$.

#### Subsection: Properties of the Transport Term

- **Theorem 2.4.5** (Transport term is linear in position). Commutator bilinearity; for fixed $\Omega_T$ the transport map is linear in $R_{\mathrm{pos}}$.
- **Corollary 2.4.3** (Transport term error). $\delta(V_{\mathrm{transport}}) = -\frac{1}{2i}[\Omega_T, \delta(R_{\mathrm{pos}})]$; magnitude bound $\leq |\omega_T|\|\delta\mathbf{r}\|$.
- **Example 2.4.1** (Transport velocity from Earth rotation). $\mathbf{r}_{\mathrm{PEF}} = (6700, 1200, 0)$ km; Cartesian and SU(2) transport; geostationary check; position error → transport velocity error.

#### Subsection: Point-Mass Reduction

- **Corollary 2.4.4** (Point-mass reduction). $\Omega_b = 0$; active state reduces to 6 quantities ($R_{\mathrm{pos}}$, $V_{\mathrm{vel}}$).

#### Subsection: Near-Identity Perturbations

- **Theorem 2.4.6** (Near-identity composition is additive). To first order in $\varepsilon$: $(I + \frac{i}{2}A)(I + \frac{i}{2}B) \approx I + \frac{i}{2}(A+B)$; error $O(\varepsilon^2)$.
- **Remark**. Application to Ch 18–19 short/long-period corrections; cross-term $O(J_2^2) \sim 10^{-6}$.

---

### §2.5 The Orbit Propagation Coordinate Transform Pipeline

#### Subsection: Frame Definitions

- **Definition 2.5.1** (Perifocal frame). $\hat{\mathbf{x}}$ toward perigee, $\hat{\mathbf{z}}$ along angular momentum; position $(r\cos\nu, r\sin\nu, 0)$.
- **Definition 2.5.2** (TEME frame). True Equator, Mean Equinox; the SGP4 propagator frame. **[A.2.1]**
- **Definition 2.5.3** (PEF frame). Pseudo Earth-Fixed; TEME rotated by $\theta_{\mathrm{GMST}}$.

#### Subsection: Transform 1 — Perifocal → TEME

- **Theorem 2.5.1** (Perifocal → TEME rotation). $M_{\mathrm{PF\to TEME}} = M_3(-\Omega_{\mathrm{node}})\cdot M_1(-i)\cdot M_3(-u)$; negative angles for passive rotation; position and velocity transforms (equations 2.55); dual quaternion forms.
- **Example 2.5.1** (ISS-like orbit: PF → TEME rotation). $\Omega = 30°$, $i = 51.6°$, $u = 45°$; three SU(2) factors computed; singularity-free commentary.
- **Remark** (Computational structure). $M_3$ diagonal → phase multiplication at 4 real ops instead of 16.
- **Proposition 2.5.1** (Singularity-free at $i = 0$). $M_1(0) = I$; $M_3(-\Omega)\cdot M_3(-u) = M_3(-\Omega-u)$; well-conditioned diagonal matrix; no Lyddane modification needed.

#### Subsection: Transform 2 — TEME → PEF

- **Theorem 2.5.2** (TEME → PEF with transport). $M_T = M_3(\theta_{\mathrm{GMST}})$; $\Omega_T = \omega_E\sigma_3$; material velocity transform (equation 2.58); transport commutator (equation 2.59) with explicit encoding $(-\omega_E r_2, \omega_E r_1, 0)$. **[M.2.1]**

#### Subsection: Pipeline Composition

- **Theorem 2.5.3** (Perifocal → PEF composition). $M_{\mathrm{PF\to PEF}} = M_3(\theta)\cdot M_3(-\Omega)\cdot M_1(-i)\cdot M_3(-u)$; first two diagonal factors combine trivially; end-to-end velocity (equation 2.61).

#### Subsection: Error Sources

Table tracing all error sources (angle measurement, sin/cos rounding, SU(2) multiply, GMST polynomial, GMST coefficients, $\omega_E$, TEME frame definition) through the SU(2) construction; symbolic bounds and binary64 values.

---

### §2.6 Error Propagation Through the Framework

#### Subsection: The Two Principles

- **Theorem 2.6.1** (Principle 1: linear operations reuse the transform). $\delta(\mathcal{L}(V)) = \mathcal{L}(\delta(V))$; exact for all linear maps; application to conjugation, commutator, addition, adjoint action.
- **Theorem 2.6.2** (Principle 2: nonlinear operations require rigorous sensitivity bounds). $|\delta(g(V))| \leq \sup|g'|\cdot|\delta(V)|$; Mean Value Theorem bound from Ch 1, Theorem 1.4.1.
- **Remark**. Application vs. construction distinction in the pipeline; sin/cos nonlinear; conjugation linear.

#### Subsection: Pipeline Classification

- **Proposition 2.6.1** (Pipeline classification). Two tables: Linear operations (10 entries: rotation, translation, transport, adjoint, spatial→material, configuration composition, secular rates, short/long-period corrections, dot product, cross product); Nonlinear operations (7 entries: angle→$M_k$, Kepler equation, sin/cos, renormalization, cube root, density model, Euclidean norm).
- **Remark** (Hidden assumption). Both-factors-uncertain requires Ch 1, Theorem 1.3.2 bilinear bound.
- **Remark** (Cross product cancellation). Near-parallel vectors; condition number; $\mathrm{rd}(\mathbf{a}\times\mathbf{b})\to 0$; angular momentum near-rectilinear case.

#### Subsection: Renormalization

- **Theorem 2.6.3** (Renormalization sensitivity). $M_{\mathrm{renorm}} = M/\sqrt{\det M} \approx M(1 - \varepsilon/2)$; double-precision accumulation estimate: $|\varepsilon| \leq 2\times 10^{-9}$ after $2\times 10^6$ compositions — renormalization rarely needed.

#### Subsection: Linearization Validity

- **Remark** (When the parallel-row principle degrades). Principle 1 exact; Principle 2 valid when $\mathrm{rd}(v) \geq d_{\min}$; fall back to Ch 1 scalar rules when $\mathrm{rd}(v) < d_{\min}$. **[A.2.2]**

#### Subsection: Three-Error Independence

- **Proposition 2.6.2** (Error categories propagate independently). $\sigma_m$, $\delta_p$, $\delta_a$ use same propagation function applied to different input magnitudes; disjoint physical mechanisms (Ch 1, Proposition 1.8.1).

---

### §2.7 The Adjoint Bridge: From SU(2) to 3×3 Matrices

#### Subsection: The Adjoint Representation

- **Definition 2.7.1** (Adjoint rotation matrix). $R_{jk}(M) = \frac{1}{2}\mathrm{tr}(\sigma_j M\sigma_k M^\dagger)$; the $(j,k)$ entry as the $j$-th component of the rotated $k$-th basis vector.
- **Theorem 2.7.1** (Adjoint equivalence). Conjugation $\leftrightarrow$ matrix-vector multiply: $\mathbf{v}' = R(M)\mathbf{v}$; proved by Pauli expansion and trace orthogonality.
- **Theorem 2.7.2** ($R(M) \in SO(3)$). Orthogonality from isometry; determinant $+1$ by continuity.
- **Theorem 2.7.3** (Adjoint is a group homomorphism). $R(M_1 M_2) = R(M_1)R(M_2)$; proved by nested conjugation.
- **Remark**. Cost comparison: SU(2) product 16 real mults vs. 3×3 product 27; convert once, apply many times.

#### Subsection: Explicit Entries

- **Theorem 2.7.4** (Explicit $R(M)$). Full $3\times 3$ matrix in quaternion components $(w,x,y,z)$ (equation 2.68); proved by computing all 9 entries from Definition 2.7.1.
- **Remark**. Alternative form in $(\alpha, \beta)$ directly (equation 2.69).
- **Example 2.7.1** (Adjoint conversion for a 45° rotation about $\hat{\mathbf{e}}_3$). Quaternion $(w,x,y,z) = (0.9239, 0, 0, 0.3827)$; all non-zero entries computed; result matches $R_3(45°)$. ✓

#### Subsection: Conversion Error

- **Theorem 2.7.5** (Conversion precision). $\|\delta(R)\|_{\max} \leq 4\epsilon_{\mathrm{mach}}$; each entry at most 2 products plus subtraction. **[P.2.6]**
- **Theorem 2.7.6** (Angle sensitivity through the bridge). $\|dR/d\alpha\| = 1$; factor-of-2 vs. SU(2) sensitivity because $R$ is quadratic in $M$.

---

### §2.8 The Linear Representation

#### Subsection: The 4×4 Homogeneous Matrix (from M̂)

- **Definition 2.8.1** (Homogeneous matrix). $H = \left(\begin{smallmatrix} R & \mathbf{t} \\ \mathbf{0}^T & 1\end{smallmatrix}\right)$; action on 4-vector; $w_H = 1$ for positions, $w_H = 0$ for displacements.
- **Remark** (Affine vs. projective interpretation). $w_H$ as affine tag, not projective weight; scalar multiples not equivalent.
- **Example 2.8.1** (Worked: rotation + translation on state and displacement). $90°$ about $\hat{\mathbf{e}}_3$ + 100 km translation; satellite at $(7000, 0, 0)$; state and displacement transformed; difference reconstruction verification. ✓
- **Remark** (Signed displacements vs. non-negative bounds). Rotation of signed displacement vs. scalar norm propagation.
- **Theorem 2.8.1** (Homogeneous composition). Block multiplication gives $R_1 R_2$ and $R_1\mathbf{t}_2 + \mathbf{t}_1$; matches dual SU(2) rule.
- **Theorem 2.8.2** (Homogeneous inverse). $H^{-1} = \left(\begin{smallmatrix} R^T & -R^T\mathbf{t} \\ \mathbf{0}^T & 1\end{smallmatrix}\right)$.

#### Subsection: The 4×4 Twist Matrix (from Ω̂_b)

- **Definition 2.8.2a** (Twist matrix). $\Xi(\hat{\Omega}_b) = \left(\begin{smallmatrix} [\boldsymbol{\omega}_b]_\times & \mathbf{v}_b \\ \mathbf{0}^T & 0\end{smallmatrix}\right)$; Lie algebra element.
- **Remark** (Conversion summary). Table: two dual quaternions to their 4×4 linear equivalents.

#### Subsection: The 7×7 State Matrix

- **Definition 2.8.2** (Skew-symmetric matrix). $[\boldsymbol{\omega}]_\times$ with $[\boldsymbol{\omega}]_\times\mathbf{v} = \boldsymbol{\omega}\times\mathbf{v}$; Cartesian form of commutator.
- **Definition 2.8.3** (The 7×7 state matrix). $T_7$ block structure (equation 2.75); action on $(\mathbf{r}, \mathbf{v}, 1)^T$ (equation 2.76); three velocity contributions named.
- **Theorem 2.8.3** (Equivalence to dual quaternion state transform). $\mathbf{t} = 0$ case gives Corollary 2.4.2 in Cartesian coordinates; $\mathbf{t} \neq 0$ case requires $\Delta\mathbf{v}$ absorption.
- **Lemma 2.8.1** (Rotation of a cross product). $R[\boldsymbol{\omega}]_\times R^T = [R\boldsymbol{\omega}]_\times$; proved via $R(\mathbf{a}\times\mathbf{b}) = (R\mathbf{a})\times(R\mathbf{b})$ for $\det R = +1$.
- **Theorem 2.8.4** ($7\times 7$ composition). Block multiplication gives $R_{12}$, $\boldsymbol{\omega}_{T,12}$, $\mathbf{t}_{12}$, $\Delta\mathbf{v}_{12}$; lower-left block uses Lemma 2.8.1.
- **Remark** (Two equivalent representations). Summary of when to prefer dual quaternion vs. matrix form; they produce identical results.

---

### §2.9 The Inverse Bridge: Shepperd's Method

- **Theorem 2.9.1** (Shepperd's method). Four-step algorithm: compute $q_k^2$ from diagonal of $R$; select largest; compute off-diagonal products; divide. Branch selection guarantees denominator $\geq \frac{1}{2}$. **[P.2.1]**
- **Remark**. Avoids classical atan2 singularity near $q_0 \approx 0$ (rotation near $180°$).

---

### §2.10 The Parallel-Row Error Architecture

#### Subsection: The Four-Row Structure

- **Definition 2.10.1** (Parallel-row state). Four rows of $(\hat{M}, \hat{\Omega}_b)$: physical, $\sigma_m$, $\delta_p$, $\delta_a$; 4×13 = 52 real numbers; linear ops use same transform; nonlinear ops multiply by sensitivity bound.
- **Remark** (Merged mode). Three error rows merge to one total row; 26 real numbers; one extra pass.

#### Subsection: Scalar Error Summaries

- **Definition 2.10.2** (Scalar error reduction). Four summaries: $\delta(\mathrm{position}) = \sqrt{-\det(\delta R_{\mathrm{pos}})}$; $\delta(\mathrm{velocity}) = \sqrt{-\det(\delta V_{\mathrm{vel}})}$; $\delta(\mathrm{orientation})$ from SU(2) tangent perturbation angle; $\delta(\mathrm{angular\,velocity}) = \sqrt{-\det(\delta\Omega)}$.
- **Remark** (Connection to reliable digits). $\mathrm{rd}(\mathrm{position}) = \lfloor -\log_{10}(\delta(\mathrm{position})/\|\mathbf{r}\|)\rfloor$; feeds Ch 1, §1.9 formula switching.

#### Subsection: Point-Mass Reduction

- **Corollary 2.10.1** (Point-mass mode). $\hat{\Omega}_b = \varepsilon V_b$; 10 active components per error row; total 43 real numbers (or 23 in merged mode).

---

### §2.11 Summary: How Subsequent Chapters Use This Framework

Narrative only. Six-step protocol for each downstream chapter. Downstream reference table mapping 12 chapters (Ch 8, 9, 13, 16–19, 20, 27, 29, 30, 32, 34–35, 36, 38) to the specific section and theorem they depend on.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Defs 1.2.1–1.2.3 | §2.10 | Error category definitions |
| Ch 1, Thms 1.3.1–1.3.4 | §2.2–2.4 | Arithmetic error composition rules |
| Ch 1, Thm 1.4.1 | §2.6 | Derivative sensitivity bound |
| Ch 1, Cors 1.4.1–1.4.2 | §2.2 | Sine/cosine bounds for rotation entries |
| Ch 1, Cor 1.4.4 | §2.9 | Square root bound for Shepperd |
| Ch 1, Cor 1.4.8 | §2.6 | Euclidean norm bound |
| Ch 1, Def 1.9.1 | §2.6 | Reliable digits metric |
| Ch 1, Prop 1.8.1 | §2.10 | Error independence |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 8 | §2.5 | Perifocal frame vectors (Thm 2.5.1) |
| Ch 9 | §2.6 | Nonlinear solver sensitivity (Thm 2.6.2) |
| Ch 13 | §2.4 | Force gradient as perturbation (Def 2.4.1) |
| Ch 16–17 | §2.4 | Secular rates as near-identity perturbations (Thm 2.4.6) |
| Ch 18–19 | §2.6 | Additive corrections via linear pipeline (Thm 2.6.1) |
| Ch 20 | §2.5 | Mean-to-osculating transform (Thm 2.5.3) |
| Ch 27 | §2.4 | Third-body force perturbation |
| Ch 29 | §2.5 | GMST as TEME→PEF rotation input (Thm 2.5.2) |
| Ch 30 | §2.2–2.8 | All frame rotations and compositions |
| Ch 32 | §2.6 | Cube root sensitivity (Thm 2.6.2) |
| Ch 34–35 | §2.5 | Full pipeline as SU(2) product chain (Thm 2.5.3) |
| Ch 36 | §2.10 | Matched-pair compatibility with state framework |
| Ch 38 | §2.6, §2.10 | Error budget from parallel-row architecture |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.2.1] | P | §2.9 | Shepperd branch selection: precision per component $\leq \sqrt{2}\,\epsilon_{\mathrm{mach}}$ |
| [P.2.2] | P | §2.2 | Inexact rotation entry: $\|\delta(M)\|_{\max} \leq \epsilon_{\mathrm{mach}}$ |
| [P.2.3] | P | §2.2 | Elementary rotation rounding: $\|\delta M\| \leq \frac{1}{2}|\delta\alpha|$ |
| [P.2.4] | P | §2.2 | Composed rotation: $\leq 16\epsilon_{\mathrm{mach}}$ total |
| [P.2.5] | P | §2.3 | Position extraction couples rotation error through $\|\mathbf{r}\|$ |
| [P.2.6] | P | §2.7 | Adjoint conversion: entry-wise $\leq 4\epsilon_{\mathrm{mach}}$ |
| [A.2.1] | A | §2.5 | TEME frame: ~0.1 arcsec from truncated precession/nutation |
| [A.2.2] | A | §2.6 | Linearization: $O(\delta^2)$ terms neglected in Principle 2 |
| [M.2.1] | M | §2.5 | GMST angle measurement propagates through TEME→PEF rotation |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 26 |
| Theorems | 48 |
| Lemmas | 4 |
| Corollaries | 6 |
| Propositions | 3 |
| Examples | 9 |
| Error Notes | 9 |
| Equations | ~82 |
| Sections | 11 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §2.1 | Build | Introduction and forward references |
| §2.2 | Build | SU(2) algebra complete |
| §2.3 | Build | Dual SU(2) complete |
| §2.4 | Build | Velocity dual quaternion complete |
| §2.5 | Build | Frame definitions complete |
| §2.6 | Build | Principles 1 and 2 complete |
| §2.7 | Build | Adjoint bridge complete |
| §2.8 | Build | 7×7 matrix complete |
| §2.9 | Build | Shepperd extraction complete |
| §2.10 | Build | Parallel-row architecture complete |
| §2.11 | Build | Summary and downstream table complete |
