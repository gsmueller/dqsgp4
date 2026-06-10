# Chapter 2: The State Matrix Framework

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathbb{H}$ | The quaternion algebra | §2.2 |
| $q = q_0 + q_1 \mathbf{i} + q_2 \mathbf{j} + q_3 \mathbf{k}$ | A quaternion with scalar part $q_0$ and vector part $\mathbf{q} = (q_1, q_2, q_3)$ | §2.2, Def. 2.2.1 |
| $q^*$ | Quaternion conjugate: $q_0 - q_1\mathbf{i} - q_2\mathbf{j} - q_3\mathbf{k}$ | §2.2, Def. 2.2.2 |
| $\|q\|$ | Quaternion norm: $\sqrt{q_0^2 + q_1^2 + q_2^2 + q_3^2}$ | §2.2, Def. 2.2.2 |
| $\mathbf{v}^{\flat}$ | A vector $\mathbf{v} \in \mathbb{R}^3$ embedded as a pure quaternion $(0, v_1, v_2, v_3)$ | §2.2, Def. 2.2.3 |
| $\varepsilon$ | The dual unit satisfying $\varepsilon^2 = 0$, $\varepsilon \neq 0$ | §2.3, Def. 2.3.1 |
| $\hat{q} = q_r + \varepsilon\, q_d$ | A dual quaternion with real part $q_r$ and dual part $q_d$ | §2.3, Def. 2.3.2 |
| $\hat{q}^*$ | Dual quaternion conjugate: $q_r^* + \varepsilon\, q_d^*$ | §2.3, Def. 2.3.3 |
| $\overline{\hat{q}}$ | Dual conjugate: $q_r - \varepsilon\, q_d$ | §2.3, Def. 2.3.3 |
| $R(q)$ | The $3 \times 3$ rotation matrix derived from unit quaternion $q$ | §2.7, Thm. 2.7.1 |
| $\omega_E$ | Earth's rotation rate $\approx 7.2921159 \times 10^{-5}$ rad/s | §2.5, Thm. 2.5.2 |
| $\theta_{\mathrm{GMST}}$ | Greenwich Mean Sidereal Time angle (Ch. 29) | §2.5, Thm. 2.5.2 |
| $\Omega, i, u$ | Right ascension, inclination, argument of latitude | §2.5, Thm. 2.5.1 |

---

## §2.1 Introduction

The SGP4 propagator transforms orbital elements through a pipeline of operations — secular update, Kepler solver, perturbation corrections, coordinate rotations, frame rotations — each modifying the satellite's position and velocity. These operations involve rotations (frame changes), translations (origin offsets), and velocity coupling (the transport theorem for rotating frames).

Chapter 1 developed error propagation for scalar operations. This chapter extends that framework to the rigid body state: position, velocity, orientation, and angular velocity. The primary algebraic tool is the **dual quaternion**, which encodes rotation and translation in a single 8-component object that composes by multiplication, has no coordinate singularities, and naturally produces the transport theorem through the product rule of its time derivative.

The chapter builds from quaternion algebra (§2.2) through dual numbers and dual quaternions (§2.3) to the full rigid body state (§2.4). Concrete SGP4 transforms are expressed as dual quaternions (§2.5). Error propagation through quaternion operations is derived from Chapter 1's scalar rules (§2.6). The traditional rotation matrix form is derived as a consequence (§2.7), not assumed as a starting point. The extended state with error tracking is defined in §2.8.

For SGP4 standard mode, the satellite is a point mass — the attitude and angular velocity components are trivial (identity orientation, zero angular velocity). The framework carries these slots for enhanced models where drag depends on body orientation (Ch. 36, Generalization).

---

## §2.2 Quaternion Algebra

**Definition 2.2.1** (Quaternion)**.** *A quaternion $q \in \mathbb{H}$ is:*

$$
q = q_0 + q_1\,\mathbf{i} + q_2\,\mathbf{j} + q_3\,\mathbf{k} \tag{2.1}
$$

*where $q_0 \in \mathbb{R}$ is the scalar part, $(q_1, q_2, q_3) \in \mathbb{R}^3$ is the vector part, and $\mathbf{i}, \mathbf{j}, \mathbf{k}$ satisfy $\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{i}\mathbf{j}\mathbf{k} = -1$.*

**Definition 2.2.2** (Conjugate, norm, inverse)**.** *For $q = q_0 + q_1\mathbf{i} + q_2\mathbf{j} + q_3\mathbf{k}$:*

$$
q^* = q_0 - q_1\mathbf{i} - q_2\mathbf{j} - q_3\mathbf{k}, \qquad \|q\| = \sqrt{q\,q^*} = \sqrt{q_0^2 + q_1^2 + q_2^2 + q_3^2}, \qquad q^{-1} = \frac{q^*}{\|q\|^2}. \tag{2.2}
$$

**Definition 2.2.3** (Pure quaternion embedding)**.** *A vector $\mathbf{v} = (v_1, v_2, v_3) \in \mathbb{R}^3$ is embedded as the pure quaternion:*

$$
\mathbf{v}^\flat = 0 + v_1\,\mathbf{i} + v_2\,\mathbf{j} + v_3\,\mathbf{k}. \tag{2.3}
$$

*The scalar part is zero. Conversely, the vector part of any pure quaternion is extracted as $\mathbf{v} = \mathrm{vec}(\mathbf{v}^\flat)$.*

**Theorem 2.2.1** (Hamilton product)**.** *The product of two quaternions $p$ and $q$ is:*

$$
pq = (p_0 q_0 - \mathbf{p} \cdot \mathbf{q}) + (p_0 \mathbf{q} + q_0 \mathbf{p} + \mathbf{p} \times \mathbf{q})^\flat \tag{2.4}
$$

*where $\mathbf{p} = (p_1, p_2, p_3)$, $\mathbf{q} = (q_1, q_2, q_3)$. The product is associative but not commutative: $pq \neq qp$ in general (the cross product term changes sign).*

*Proof concept.* Expand using $\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = -1$, $\mathbf{i}\mathbf{j} = \mathbf{k}$, $\mathbf{j}\mathbf{k} = \mathbf{i}$, $\mathbf{k}\mathbf{i} = \mathbf{j}$ and collect terms. The scalar part is a dot product; the vector part is a sum of scaled vectors plus a cross product. ∎

**Theorem 2.2.2** (Quaternion rotation)**.** *Let $q$ be a unit quaternion ($\|q\| = 1$) parameterized as:*

$$
q = \cos\frac{\theta}{2} + \sin\frac{\theta}{2}\,(n_1\,\mathbf{i} + n_2\,\mathbf{j} + n_3\,\mathbf{k}) \tag{2.5}
$$

*where $\hat{\mathbf{n}} = (n_1, n_2, n_3)$ is a unit vector and $\theta$ is the rotation angle. Then the sandwich product:*

$$
\mathbf{v}'^\flat = q\,\mathbf{v}^\flat\,q^* \tag{2.6}
$$

*rotates $\mathbf{v}$ by angle $\theta$ about axis $\hat{\mathbf{n}}$. The result $\mathbf{v}'^\flat$ is a pure quaternion, and $\mathbf{v}' = \mathrm{vec}(\mathbf{v}'^\flat)$.*

*Proof concept.* Expand $q\,\mathbf{v}^\flat\,q^*$ using (2.4). Show that the scalar part vanishes (the result is pure) and the vector part equals Rodrigues' rotation formula: $\mathbf{v}' = \mathbf{v}\cos\theta + (\hat{\mathbf{n}} \times \mathbf{v})\sin\theta + \hat{\mathbf{n}}(\hat{\mathbf{n}} \cdot \mathbf{v})(1 - \cos\theta)$. ∎

**Theorem 2.2.3** (Composition of rotations)**.** *If $q_1$ rotates frame B → A and $q_2$ rotates frame C → B, then:*

$$
q_{AC} = q_1 \cdot q_2 \tag{2.7}
$$

*rotates frame C → A. Rotation composition is quaternion multiplication.*

*Proof concept.* Apply the sandwich product twice: $\mathbf{v}_A = q_1(q_2\,\mathbf{v}_C\,q_2^*)\,q_1^* = (q_1 q_2)\,\mathbf{v}_C\,(q_1 q_2)^*$. ∎

**Theorem 2.2.4** (Inverse rotation)**.** *For a unit quaternion, $q^{-1} = q^*$. The inverse rotation is the conjugate.*

*Proof concept.* $\|q\| = 1$ implies $q^{-1} = q^*/\|q\|^2 = q^*$. Verify: $q\,q^* = q^*\,q = 1$. ∎

**Lemma 2.2.1** (Double cover)**.** *$q$ and $-q$ represent the same rotation: $q\,\mathbf{v}^\flat\,q^* = (-q)\,\mathbf{v}^\flat\,(-q)^*$. The map from unit quaternions to $SO(3)$ is two-to-one.*

*Proof concept.* The signs cancel in the sandwich product. ∎

**Definition 2.2.4** (Elementary quaternion rotations)**.** *Rotation by angle $\alpha$ about coordinate axis $k$:*

$$
q_1(\alpha) = \cos\frac{\alpha}{2} + \sin\frac{\alpha}{2}\,\mathbf{i}, \qquad q_2(\alpha) = \cos\frac{\alpha}{2} + \sin\frac{\alpha}{2}\,\mathbf{j}, \qquad q_3(\alpha) = \cos\frac{\alpha}{2} + \sin\frac{\alpha}{2}\,\mathbf{k}. \tag{2.8}
$$

---

## §2.3 Dual Numbers and Dual Quaternions

**Definition 2.3.1** (Dual number)**.** *A dual number is $\hat{a} = a + \varepsilon\,b$ where $a, b \in \mathbb{R}$ and $\varepsilon$ is the dual unit satisfying $\varepsilon^2 = 0$, $\varepsilon \neq 0$. Arithmetic:*

$$
(a + \varepsilon b)(c + \varepsilon d) = ac + \varepsilon(ad + bc). \tag{2.9}
$$

*The $\varepsilon^2 = 0$ property means dual numbers naturally carry first-order perturbation information: $f(a + \varepsilon b) = f(a) + \varepsilon\,b\,f'(a)$ for any analytic function $f$.*

**Definition 2.3.2** (Dual quaternion)**.** *A dual quaternion is:*

$$
\hat{q} = q_r + \varepsilon\,q_d \tag{2.10}
$$

*where $q_r, q_d \in \mathbb{H}$ are ordinary quaternions. This is an 8-dimensional object: 4 components from $q_r$, 4 from $q_d$.*

**Definition 2.3.3** (Dual quaternion conjugates)**.** *Three conjugation operations:*

1. *Quaternion conjugate: $\hat{q}^* = q_r^* + \varepsilon\,q_d^*$ (negate vector parts)*
2. *Dual conjugate: $\overline{\hat{q}} = q_r - \varepsilon\,q_d$ (negate dual part)*
3. *Combined conjugate: $\overline{\hat{q}}^* = q_r^* - \varepsilon\,q_d^*$*

**Definition 2.3.4** (Unit dual quaternion)**.** *A unit dual quaternion satisfies:*

$$
\|q_r\| = 1, \qquad q_r \cdot q_d = q_{r,0}\,q_{d,0} + q_{r,1}\,q_{d,1} + q_{r,2}\,q_{d,2} + q_{r,3}\,q_{d,3} = 0. \tag{2.11}
$$

*Two constraints on 8 components → 6 degrees of freedom, matching a rigid body configuration (3 translational + 3 rotational).*

**Theorem 2.3.1** (Rigid body transform via dual quaternion)**.** *A unit dual quaternion $\hat{q} = q_r + \varepsilon\,q_d$ encodes rotation by $q_r$ and translation by $\mathbf{t}$, where:*

$$
q_d = \frac{1}{2}\,\mathbf{t}^\flat\,q_r. \tag{2.12}
$$

*The action on a point $\mathbf{p}$ (embedded as dual quaternion $1 + \varepsilon\,\mathbf{p}^\flat$) is:*

$$
1 + \varepsilon\,\mathbf{p}'^\flat = \hat{q}\,(1 + \varepsilon\,\mathbf{p}^\flat)\,\overline{\hat{q}}^* \tag{2.13}
$$

*yielding $\mathbf{p}' = q_r\,\mathbf{p}^\flat\,q_r^* + \mathbf{t}^\flat$ — rotation then translation.*

*Proof concept.* Expand (2.13) using the dual quaternion product. The real part gives $q_r \cdot 1 \cdot q_r^* = 1$ (identity). The dual part gives $q_r\,\mathbf{p}^\flat\,q_r^* + 2\,q_d\,q_r^*$. By (2.12), $2\,q_d\,q_r^* = \mathbf{t}^\flat\,q_r\,q_r^* = \mathbf{t}^\flat$. ∎

**Theorem 2.3.2** (Composition of dual quaternions)**.** *Let $\hat{q}_1$ transform frame B → A and $\hat{q}_2$ transform frame C → B. Then:*

$$
\hat{q}_{AC} = \hat{q}_1 \cdot \hat{q}_2 \tag{2.14}
$$

*correctly composes both rotation and translation. The combined rotation is $q_{r,1}\,q_{r,2}$ and the combined translation is $q_{r,1}\,\mathbf{t}_2^\flat\,q_{r,1}^* + \mathbf{t}_1^\flat$ (inner translation rotated into outer frame, then outer translation added).*

*Proof concept.* Multiply $\hat{q}_1 \hat{q}_2 = (q_{r,1} + \varepsilon q_{d,1})(q_{r,2} + \varepsilon q_{d,2}) = q_{r,1}q_{r,2} + \varepsilon(q_{r,1}q_{d,2} + q_{d,1}q_{r,2})$. Extract the translation from the dual part using (2.12). ∎

**Theorem 2.3.3** (Inverse of a dual quaternion)**.** *For a unit dual quaternion:*

$$
\hat{q}^{-1} = \overline{\hat{q}}^* = q_r^* - \varepsilon\,q_d^* \tag{2.15}
$$

*Proof concept.* Verify $\hat{q}\,\hat{q}^{-1} = 1$ using the unit constraints (2.11). ∎

**Theorem 2.3.4** (Equivalence to 4×4 homogeneous matrix)**.** *Every unit dual quaternion $\hat{q}$ corresponds to a unique 4×4 homogeneous matrix $H = [R(\hat{q}), \mathbf{t}(\hat{q});\, \mathbf{0}^T, 1]$, and vice versa (up to the double cover sign). Composition of dual quaternions corresponds to matrix multiplication of homogeneous matrices.*

*Proof concept.* Construct $R$ from $q_r$ via the rotation matrix formula (Theorem 2.7.1). Extract $\mathbf{t} = 2\,\mathrm{vec}(q_d\,q_r^*)$. Verify that dual quaternion multiplication gives the same rotation and translation as the 4×4 block product. ∎

---

## §2.4 The Rigid Body State

**Definition 2.4.1** (Rigid body state)**.** *The full rigid body state is:*

$$
\mathcal{S} = (\hat{q},\, \mathbf{v},\, \boldsymbol{\omega}) \tag{2.16}
$$

*where:*
- *$\hat{q} = q_r + \varepsilon\,q_d$ is a unit dual quaternion encoding position and orientation (6 DOF)*
- *$\mathbf{v} \in \mathbb{R}^3$ is the translational velocity*
- *$\boldsymbol{\omega} \in \mathbb{R}^3$ is the angular velocity of the body*

*This is 14 components $(8 + 3 + 3)$ with 2 constraints, giving 12 DOF.*

**Remark.** For SGP4 standard mode (point mass), the attitude part of $q_r$ is identity (no body rotation beyond the orbital frame rotation), and $\boldsymbol{\omega}_{\mathrm{body}} = \mathbf{0}$. The dual quaternion reduces to encoding frame rotation + satellite position only.

**Assumption 2.4.1** (Instantaneous angular velocity)**.** *The angular velocity $\boldsymbol{\omega}$ is treated as constant during a single transformation. For time-varying $\boldsymbol{\omega}$, the transform must be re-evaluated at each epoch.*

**Theorem 2.4.1** (Transport theorem from dual quaternion derivative)**.** *Let $\hat{q}(t)$ describe a time-varying rigid body transform (e.g., Earth rotation). The time derivative of the composed transform $\hat{q}_1(t) \cdot \hat{q}_2$ (where $\hat{q}_2$ is constant) is:*

$$
\frac{d}{dt}(\hat{q}_1 \hat{q}_2) = \dot{\hat{q}}_1\,\hat{q}_2 \tag{2.17}
$$

*The derivative $\dot{\hat{q}}_1$ encodes both the angular velocity and linear velocity of frame 1. For a purely rotating frame with angular velocity $\boldsymbol{\omega}$:*

$$
\dot{q}_r = \frac{1}{2}\,\boldsymbol{\omega}^\flat\,q_r \tag{2.18}
$$

*The velocity of a point $\mathbf{r}$ in the rotating frame, as seen from the inertial frame, is:*

$$
\mathbf{v}_{\mathrm{inertial}} = q_r\,\mathbf{v}_{\mathrm{body}}^\flat\,q_r^* + \boldsymbol{\omega} \times (q_r\,\mathbf{r}^\flat\,q_r^*) \tag{2.19}
$$

*The cross-product term $\boldsymbol{\omega} \times \mathbf{r}'$ is the transport theorem — it arises from the product rule applied to $\hat{q}(t)\,\mathbf{p}\,\hat{q}(t)^{-1}$.*

*Proof concept.* Differentiate the action (2.13) with respect to time. Apply the product rule to $\hat{q}(t)\,(1 + \varepsilon\,\mathbf{p}^\flat)\,\overline{\hat{q}}(t)^*$. The time derivative of $q_r(t)\,\mathbf{v}^\flat\,q_r(t)^*$ gives $\dot{q}_r\,\mathbf{v}^\flat\,q_r^* + q_r\,\mathbf{v}^\flat\,\dot{q}_r^*$. Substitute $\dot{q}_r = \frac{1}{2}\boldsymbol{\omega}^\flat q_r$ from (2.18) and simplify to obtain the cross-product term. ∎

**Theorem 2.4.2** (Composition of stateful transforms)**.** *Let $\hat{q}_1(t)$ describe a time-varying frame (with velocity/angular velocity) and $\hat{q}_2$ a static transform. The composed transform $\hat{q}_{12} = \hat{q}_1 \hat{q}_2$ has:*

- *Configuration: $\hat{q}_{12}$ (dual quaternion multiplication, Theorem 2.3.2)*
- *Angular velocity: $\boldsymbol{\omega}_{12} = \boldsymbol{\omega}_1 + q_1\,\boldsymbol{\omega}_2^\flat\,q_1^*$ (angular velocities add, with inner frame's ω rotated into outer frame)*
- *Linear velocity: $\mathbf{v}_{12} = \mathbf{v}_1 + q_1\,\mathbf{v}_2^\flat\,q_1^* + \boldsymbol{\omega}_1 \times (q_1\,\Delta\mathbf{r}_2^\flat\,q_1^*)$ (transport term appears)*

*Proof concept.* Differentiate $\hat{q}_1(t)\,\hat{q}_2(t)$ and extract the velocity components from the dual part. The transport term arises from $\dot{q}_{r,1}$ acting on the translation of $\hat{q}_2$. ∎

**Corollary 2.4.1** (Pure rotation composition)**.** *When both transforms are pure rotations (no translation, no angular velocity), composition reduces to quaternion multiplication: $q_{12} = q_1 q_2$. No transport terms arise.*

**Corollary 2.4.2** (Earth rotation transport)**.** *For Earth rotation about $\hat{\mathbf{e}}_3$ at rate $\omega_E$, the quaternion is $q_E(t) = q_3(\omega_E t)$ and $\boldsymbol{\omega} = \omega_E\,\hat{\mathbf{e}}_3$. The transport term is $\omega_E\,\hat{\mathbf{e}}_3 \times \mathbf{r}' = \omega_E(-r'_y, r'_x, 0)^T$.*

---

## §2.5 SGP4 Coordinate Transforms as Dual Quaternions

**Definition 2.5.1** (TEME frame)**.** *The True Equator, Mean Equinox (TEME) frame is the reference frame of SGP4 output. It is not a precisely defined inertial frame.* [A.2.1]

**Definition 2.5.2** (PEF frame)**.** *The Pseudo Earth-Fixed (PEF) frame is obtained from TEME by rotating about the polar axis by the Greenwich Mean Sidereal Time angle $\theta_{\mathrm{GMST}}$ (Ch. 29).*

**Definition 2.5.3** (Perifocal frame)**.** *The perifocal frame has $\hat{\mathbf{x}}$ toward perigee, $\hat{\mathbf{z}}$ along the angular momentum vector, and is defined by the orbital elements $(\Omega, i, \omega)$ with argument of latitude $u = \omega + \nu$.*

**Theorem 2.5.1** (Perifocal → TEME quaternion)**.** *The rotation from perifocal to TEME is a pure rotation quaternion:*

$$
q_{\mathrm{PF} \to \mathrm{TEME}} = q_3(-\Omega) \cdot q_1(-i) \cdot q_3(-u) \tag{2.20}
$$

*This is the composition of three elementary quaternion rotations (Definition 2.2.4). The dual quaternion has zero translation: $\hat{q} = q_{\mathrm{PF} \to \mathrm{TEME}} + \varepsilon\,0$.*

*Proof concept.* The perifocal-to-ECI rotation is the 3-1-3 Euler sequence $(\Omega, i, u)$ applied in reverse. By Theorem 2.2.3, this is the product of the three elementary quaternion rotations. Each $q_k(-\alpha) = \cos(\alpha/2) - \sin(\alpha/2)\,\mathbf{e}_k^\flat$ (note the sign). The product yields a single unit quaternion encoding the full rotation. No translation is involved (both frames share an origin at Earth's center). ∎

**Remark** (Connection to unit vectors $\hat{\mathbf{U}}$, $\hat{\mathbf{V}}$)**.** The unit vectors used in Spacetrack Report No. 3 are $\hat{\mathbf{U}} = q\,\hat{\mathbf{e}}_1^\flat\,q^*$ and $\hat{\mathbf{V}} = q\,\hat{\mathbf{e}}_2^\flat\,q^*$ — the standard basis vectors rotated by $q$. Position and velocity in TEME are $\mathbf{r} = r\,\hat{\mathbf{U}}$ and $\mathbf{v} = \dot{r}\,\hat{\mathbf{U}} + r\dot{f}\,\hat{\mathbf{V}}$.

**Remark** (Trigonometric errors)**.** The quaternion components involve $\cos(\alpha/2)$ and $\sin(\alpha/2)$ of the orbital angles. Errors bounded by Ch. 1, Corollaries 1.4.1–1.4.2. The half-angle form means the trig arguments are half those of the DCM entries — this can improve numerical behavior for angles near $\pi$. [P.2.1]

**Theorem 2.5.2** (TEME → PEF dual quaternion)**.** *The transformation from TEME to PEF is a rotating-frame dual quaternion:*

$$
\hat{q}_{\mathrm{TEME} \to \mathrm{PEF}}(t) = q_3(\theta_{\mathrm{GMST}}(t)) + \varepsilon\,0 \tag{2.21}
$$

*with angular velocity $\boldsymbol{\omega} = \omega_E\,\hat{\mathbf{e}}_3$. The transport term (Theorem 2.4.1) produces the velocity correction $\omega_E\,\hat{\mathbf{e}}_3 \times \mathbf{r}_{\mathrm{PEF}}$.* [M.2.1]

*Proof concept.* PEF is TEME rotated by $\theta_{\mathrm{GMST}}$ about $\hat{\mathbf{e}}_3$. The quaternion is $q_3(\theta)$. No translation ($\Delta\mathbf{r} = \mathbf{0}$). The frame rotates at $\omega_E$, so by Theorem 2.4.1 (Eq. 2.18), $\dot{q}_r = \frac{1}{2}\omega_E\,\hat{\mathbf{e}}_3^\flat\,q_r$, and the velocity transform includes the transport term (Eq. 2.19). ∎

**Theorem 2.5.3** (Perifocal → PEF composition)**.** *The composed transform is:*

$$
\hat{q}_{\mathrm{PF} \to \mathrm{PEF}} = \hat{q}_{\mathrm{TEME} \to \mathrm{PEF}} \cdot \hat{q}_{\mathrm{PF} \to \mathrm{TEME}} \tag{2.22}
$$

*The rotation is $q_3(\theta_{\mathrm{GMST}}) \cdot q_3(-\Omega) \cdot q_1(-i) \cdot q_3(-u)$. The transport term from Earth rotation acts on the fully rotated position.*

*Proof concept.* Apply Theorem 2.4.2. The perifocal→TEME transform is static (no angular velocity); the TEME→PEF transform has $\boldsymbol{\omega} = \omega_E\hat{\mathbf{e}}_3$. The composed angular velocity is $\omega_E\hat{\mathbf{e}}_3$ (only the outer frame rotates). The transport term is $\omega_E\hat{\mathbf{e}}_3 \times \mathbf{r}'$ where $\mathbf{r}'$ is the fully rotated position. ∎

**Proposition 2.5.1** (Near-identity perturbation)**.** *The short-period corrections of Brouwer's theory (Ch. 18) can be expressed as a near-identity dual quaternion: $\hat{q} \approx 1 + \varepsilon\,\delta\hat{q}$ where $\|\delta\hat{q}\| \sim J_2/p^2 \ll 1$. Composition of near-identity dual quaternions is approximately additive: $(1 + \varepsilon\,\delta\hat{q}_1)(1 + \varepsilon\,\delta\hat{q}_2) \approx 1 + \varepsilon(\delta\hat{q}_1 + \delta\hat{q}_2)$.*

*Proof concept.* The corrections $\Delta r$, $\Delta u$, $\Delta\Omega$, $\Delta i$ are $O(J_2/p^2) \sim 10^{-3}$. A quaternion near identity is $q \approx 1 + \frac{1}{2}\boldsymbol{\theta}^\flat$ for small rotation $\boldsymbol{\theta}$. Products of near-identity elements: $(1+a)(1+b) = 1 + a + b + ab \approx 1 + a + b$. Full development in Ch. 18. ∎

**Remark** (PEF → ECEF)**.** Polar motion corrections ($< 0.5$ arcsec) are negligible for SGP4-class accuracy. The TEME frame accuracy is limited by truncated precession/nutation [A.2.1].

**Remark** (GMST and unit conversions)**.** $\theta_{\mathrm{GMST}}$ computed from the IAU 1982 polynomial (Ch. 29). Unit conversion factors are Tier I within the matched pair (Ch. 3, Ch. 37).

---

## §2.6 Error Propagation Principles

The key result of this section is that the error state transforms by the same operations as the physical state for all linear steps in the pipeline, and by rigorous sensitivity bounds at nonlinear steps. This eliminates the need for a separate error propagation framework.

**Definition 2.6.1** (Linear operation)**.** *An operation $f$ on the state is linear if $f(\mathcal{S} + \delta\mathcal{S}) = f(\mathcal{S}) + f(\delta\mathcal{S})$ for all perturbations $\delta\mathcal{S}$. Rotation, translation, and the transport term are linear operations on the state.*

**Definition 2.6.2** (Nonlinear operation)**.** *An operation $f$ is nonlinear if the above identity does not hold. The Kepler equation solver, trigonometric evaluations, and renormalization are nonlinear operations.*

**Theorem 2.6.1** (Error propagation through linear operations)**.** *Let $f$ be a linear operation on the state. Then:*

$$
f(\mathcal{S} + \delta\mathcal{S}) - f(\mathcal{S}) = f(\delta\mathcal{S}) \tag{2.23}
$$

*exactly. The error state $\delta\mathcal{S}$ transforms by the same operation $f$ that transforms the physical state $\mathcal{S}$. No approximation is involved.*

*Proof concept.* Linearity. The rotation $q\,(\mathbf{v} + \boldsymbol{\delta})^\flat\,q^* = q\,\mathbf{v}^\flat\,q^* + q\,\boldsymbol{\delta}^\flat\,q^*$. The transport term $\boldsymbol{\omega} \times (\mathbf{r} + \boldsymbol{\delta}) = \boldsymbol{\omega} \times \mathbf{r} + \boldsymbol{\omega} \times \boldsymbol{\delta}$. Translation $(\mathbf{r} + \boldsymbol{\delta}) + \Delta\mathbf{r} = (\mathbf{r} + \Delta\mathbf{r}) + \boldsymbol{\delta}$. In each case the perturbation passes through the same operation. ∎

**Corollary 2.6.1** (Rotation preserves error magnitude)**.** *A unit quaternion rotation is an isometry: $\|q\,\boldsymbol{\delta}^\flat\,q^*\| = \|\boldsymbol{\delta}\|$. Rotation does not amplify or reduce the total error.*

*Proof concept.* Must be proved from the properties of the sandwich product and the unit-norm constraint. The sandwich product by a unit quaternion preserves the norm of pure quaternions. ∎

**Theorem 2.6.2** (Error propagation through nonlinear operations)**.** *Let $f$ be a nonlinear operation. The error transforms by:*

$$
f(\mathcal{S} + \delta\mathcal{S}) - f(\mathcal{S}) = J_f(\mathcal{S})\,\delta\mathcal{S} + O(\|\delta\mathcal{S}\|^2) \tag{2.24}
$$

*where $J_f$ is the Jacobian of $f$ evaluated at $\mathcal{S}$. A rigorous error bound requires:*

$$
\|f(\mathcal{S} + \delta\mathcal{S}) - f(\mathcal{S})\| \leq B_f(\mathcal{S})\,\|\delta\mathcal{S}\| \tag{2.25}
$$

*where $B_f(\mathcal{S})$ is a proven upper bound on the sensitivity of $f$ at $\mathcal{S}$. This bound must be derived from the specific properties of $f$ — it cannot be assumed or approximated.*

*Proof concept.* The Mean Value Theorem (Ch. 1, Theorem 1.4.1) or its multivariate generalization (Ch. 1, Theorem 1.4.2) provides the existence of $B_f$. The specific value of $B_f$ depends on the operation and must be derived in the chapter where that operation is developed. ∎

**Remark.** The pipeline decomposes into linear and nonlinear steps:

| Step | Type | Error propagation |
|------|------|-------------------|
| Frame rotation (§2.5) | Linear | Same operation |
| Translation | Linear | Same operation (adds zero to error) |
| Transport term $\boldsymbol{\omega} \times \mathbf{r}$ | Linear | Same operation |
| Secular update (Ch. 33) | Linear (time polynomial) | Same operation |
| Kepler solver (Ch. 9) | Nonlinear | Sensitivity bound $B_f$ required |
| Trigonometric evaluation | Nonlinear | Sensitivity bound $B_f$ required |
| Short-period corrections (Ch. 18) | Nonlinear | Sensitivity bound $B_f$ required |
| Renormalization | Nonlinear | Sensitivity bound $B_f$ required |

For each nonlinear step, the bound $B_f$ is derived in the chapter that develops that step. This chapter (Ch. 2) establishes only the framework — the two principles above — not the specific bounds.

**Theorem 2.6.3** (Independence of error categories)**.** *The three error categories $\sigma_m$, $\delta_p$, $\delta_a$ propagate independently through both linear and nonlinear operations. For linear operations, this follows from linearity: $f(\delta\mathcal{S})$ depends only on $\delta\mathcal{S}$, regardless of which category it represents. For nonlinear operations, the bound $B_f$ is a property of the operation and the physical state — it does not depend on which error category is being propagated.*

*Proof concept.* The three categories share the same state-space structure (same 14 components). The operations act on that structure, not on the error labels. Since the same operation (or same bound) applies regardless of category, no cross-category coupling occurs. This is the vector generalization of Ch. 1, Proposition 1.8.1. ∎

**Theorem 2.6.4** (Dot product, cross product, and norm)**.** *The dot product $\mathbf{a} \cdot \mathbf{b}$, cross product $\mathbf{a} \times \mathbf{b}$, and Euclidean norm $\|\mathbf{v}\|$ are compositions of the arithmetic operations of Ch. 1, §1.3. Rigorous error bounds for each must be derived from Ch. 1, Theorems 1.3.1–1.3.2 (addition, multiplication) applied to their component-wise expansions.*

*Key properties to be established:*

1. *Dot product: sum of products. Error bound follows from Ch. 1, Theorems 1.3.1–1.3.2.*
2. *Cross product: differences of products. Subtractive cancellation (Ch. 1, §1.7) occurs when $\mathbf{a} \approx \lambda\mathbf{b}$. The bound must quantify when the cross product loses reliability.*
3. *Euclidean norm: square root of sum of squares. The norm does not amplify the total error magnitude (to be proved via Ch. 1, Theorem 1.4.2). At $\|\mathbf{v}\| \to 0$, see Ch. 1, Definition 1.9.2.*

*Downstream need:* Ch. 8 (orbital radius, angular momentum), Ch. 20 (osculating elements), Ch. 27 (third-body angle), Ch. 38 (state vector error output).

---

## §2.7 The Rotation Matrix as a Derived Form

**Theorem 2.7.1** (Rotation matrix from quaternion)**.** *For a unit quaternion $q = (q_0, q_1, q_2, q_3)$, the equivalent $3 \times 3$ rotation matrix is:*

$$
R(q) = \begin{bmatrix} q_0^2 + q_1^2 - q_2^2 - q_3^2 & 2(q_1 q_2 - q_0 q_3) & 2(q_1 q_3 + q_0 q_2) \\ 2(q_1 q_2 + q_0 q_3) & q_0^2 - q_1^2 + q_2^2 - q_3^2 & 2(q_2 q_3 - q_0 q_1) \\ 2(q_1 q_3 - q_0 q_2) & 2(q_2 q_3 + q_0 q_1) & q_0^2 - q_1^2 - q_2^2 + q_3^2 \end{bmatrix} \tag{2.30}
$$

*This matrix satisfies $R^T R = I$, $\det R = +1$, and $R(q)\,\mathbf{v} = \mathrm{vec}(q\,\mathbf{v}^\flat\,q^*)$.*

*Proof concept.* Expand the sandwich product $q\,\mathbf{v}^\flat\,q^*$ component by component using (2.4). Collect coefficients of $v_1, v_2, v_3$ to form the matrix entries. Orthogonality and unit determinant follow from $\|q\| = 1$. ∎

**Theorem 2.7.2** (7×7 state matrix from dual quaternion)**.** *The 7×7 state matrix of the previous formulation:*

$$
T = \begin{bmatrix} R & \mathbf{0} & \Delta\mathbf{r} \\ \dot{R} & R & \Delta\mathbf{v} \\ \mathbf{0}^T & \mathbf{0}^T & 1 \end{bmatrix} \tag{2.31}
$$

*is the component-wise expansion of the dual quaternion action (Theorem 2.3.1) and its time derivative (Theorem 2.4.1). Specifically:*

- *$R = R(q_r)$ from (2.30)*
- *$\Delta\mathbf{r} = \mathrm{vec}(2\,q_d\,q_r^*)$ from (2.12)*
- *$\dot{R} = [\boldsymbol{\omega}]_\times\,R$ from (2.18), where $[\boldsymbol{\omega}]_\times$ is the skew-symmetric matrix of $\boldsymbol{\omega}$*
- *$\Delta\mathbf{v}$ from the time derivative of the translation*

*The 7×7 matrix multiply of the old formulation is algebraically identical to dual quaternion composition + velocity extraction. The matrix form is a computational expansion of the dual quaternion form, not an independent construct.*

*Proof concept.* Write out the dual quaternion action (2.13) in component form using (2.30) for the rotation part and $2\,q_d\,q_r^*$ for the translation. Differentiate with respect to time to obtain the velocity row. Assemble into the $7 \times 7$ block structure. ∎

**Remark** (When to use which form)**.** The dual quaternion form is preferred for:

- Composition of many transforms (compact, no matrix storage)
- Singularity-free attitude tracking
- Renormalization (divide by scalar norm vs. re-orthogonalize a matrix)

The matrix form is preferred for:

- Explicit coordinate computations (direct multiplication)
- Interfacing with legacy code that expects 3×3 DCMs
- Pedagogical derivations where the block structure is instructive

The error propagation through TrackedValue is representation-agnostic at the scalar level — it tracks errors through whichever arithmetic operations are actually performed.

**Remark** (Quaternion from matrix — Shepperd's method)**.** Given $R \in SO(3)$, extract $q$ by choosing the largest of $\{q_0, q_1, q_2, q_3\}$ to avoid division by near-zero. The four candidates are:

- $4q_0^2 = 1 + R_{11} + R_{22} + R_{33}$
- $4q_1^2 = 1 + R_{11} - R_{22} - R_{33}$
- $4q_2^2 = 1 - R_{11} + R_{22} - R_{33}$
- $4q_3^2 = 1 - R_{11} - R_{22} + R_{33}$

Pick the largest, take the square root, then compute the remaining three from off-diagonal entries. This is always numerically stable.

---

## §2.8 The Extended State with Error Tracking

**Definition 2.8.1** (The augmented state)**.** *The augmented state consists of four row vectors of equal dimension:*

$$
\mathbf{S} = \begin{bmatrix} \mathbf{s} \\ \boldsymbol{\sigma}_m \\ \boldsymbol{\delta}_p \\ \boldsymbol{\delta}_a \end{bmatrix} \in \mathbb{R}^{4 \times 14} \tag{2.32}
$$

*where:*
- *Row 1: $\mathbf{s} = (q_{r,0},\, q_{r,1},\, q_{r,2},\, q_{r,3},\, q_{d,0},\, q_{d,1},\, q_{d,2},\, q_{d,3},\, v_x,\, v_y,\, v_z,\, \omega_x,\, \omega_y,\, \omega_z)$ — the physical state*
- *Row 2: $\boldsymbol{\sigma}_m$ — measurement error bounds, same 14 components*
- *Row 3: $\boldsymbol{\delta}_p$ — precision error bounds, same 14 components*
- *Row 4: $\boldsymbol{\delta}_a$ — accuracy error bounds, same 14 components*

*All four rows share the same structure. The error rows carry non-negative values representing upper bounds on the deviation of each component from its true value, classified by source (Ch. 1, Definitions 1.2.1–1.2.3).*

**Theorem 2.8.1** (Independent row propagation)**.** *Let $f$ be any operation in the propagation pipeline. The augmented state transforms as:*

- *For linear $f$ (Theorem 2.6.1): $f$ is applied identically to each row of $\mathbf{S}$. The four rows do not interact.*
- *For nonlinear $f$ (Theorem 2.6.2): $f$ is applied to row 1 (the physical state). The sensitivity bound $B_f(\mathbf{s})$, evaluated at the physical state, is applied to rows 2–4. The four rows do not interact.*

*In both cases, no row depends on any other row during the transform. The rows are independent.*

*Proof concept.* For linear operations, this follows from Theorem 2.6.1: $f(\mathbf{s} + \boldsymbol{\delta}) - f(\mathbf{s}) = f(\boldsymbol{\delta})$, so the error row transforms by the same $f$. For nonlinear operations, the bound $B_f(\mathbf{s})$ depends on row 1 (the physical state at which the sensitivity is evaluated) but not on the error rows. The bound is then applied as a scalar multiplication to each error row independently. The three error rows never interact with each other (Theorem 2.6.3). ∎

**Remark.** At nonlinear steps, the bound $B_f(\mathbf{s})$ is evaluated once from the physical state and then applied to all three error rows. This is a single evaluation of the bound, not three.

**Remark.** For a sequence of $n$ operations $f_1, f_2, \ldots, f_n$, the augmented state propagates as $\mathbf{S}_{\mathrm{out}} = f_n \circ \cdots \circ f_2 \circ f_1(\mathbf{S}_{\mathrm{in}})$, where each $f_k$ acts on all four rows independently (with the sensitivity bound for nonlinear $f_k$ evaluated at the physical state after step $k-1$). The total error at the end is valid by Ch. 1, Theorem 1.3.4 (composition of error bounds).

**Definition 2.8.2** (Scalar error summaries)**.** *The 14-component error rows are summarized into physically meaningful scalar bounds by Euclidean norm over each group:*

| Group | Components | Scalar bound |
|-------|-----------|-------------|
| Translational position | columns 5–7 of $q_d$ (after extracting $\mathbf{r}$ from dual part) | $\delta(\mathbf{r})$ |
| Orientation | columns 1–3 of $q_r$ (vector part) | $\delta(\mathrm{att})$ |
| Translational velocity | columns 9–11 | $\delta(\mathbf{v})$ |
| Angular velocity | columns 12–14 | $\delta(\boldsymbol{\omega})$ |

*Each scalar bound is computed for each error category separately, giving $4 \times 3 = 12$ summary values. The total error for reliable-digits computation (Ch. 1, Definition 1.9.1) is the sum across categories: $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$ per group.*

**Remark** (Position extraction from dual quaternion)**.** The position $\mathbf{r}$ is not stored directly — it is encoded in the dual part $q_d$ via the relation $\mathbf{t} = 2\,\mathrm{vec}(q_d\,q_r^*)$ (Theorem 2.3.1). Extracting $\mathbf{r}$ from $q_d$ is a nonlinear operation (it involves multiplication by $q_r^*$). The position error group must account for this extraction. For the error summary, the sensitivity of the extraction must be bounded (it depends on $q_r$, which is a unit quaternion, so the bound is straightforward).

**Remark** (SGP4 standard mode)**.** For point-mass propagation, the attitude part of $q_r$ is identity $(1, 0, 0, 0)$ with zero error, and $\boldsymbol{\omega}_{\mathrm{body}} = \mathbf{0}$ with zero error. The corresponding columns in the error rows are zero. The augmented state has only 8 active columns (the dual quaternion encoding frame rotation + position) plus 3 for velocity. The remaining columns carry zero and propagate as zero through all operations.

**Remark** (Enhanced mode — drag-attitude coupling)**.** When drag depends on body attitude (Ch. 21, Generalization), the cross-sectional area $A$ becomes a function of the attitude columns. This couples attitude error into the drag force (a nonlinear step). The sensitivity bound $B_{\mathrm{drag}}$ at that step will have nonzero entries coupling the attitude columns to the velocity columns. The framework handles this without modification — it is simply a nonzero bound at one nonlinear step in the pipeline.

---

## §2.9 Summary and Usage Guide

This chapter established the dual quaternion as the primary algebraic framework for representing, composing, and error-analyzing all coordinate transforms and state propagation steps.

**The "State Matrix Formulation" section pattern.** Every subsequent chapter that performs a coordinate transformation or state propagation step includes a section titled "State Matrix Formulation." That section expresses the chapter's operation as a dual quaternion (or composition of dual quaternions), with the equivalent $7 \times 7$ matrix form derived via Theorem 2.7.2 when the block structure is instructive.

**The composition principle.** The SGP4 propagation pipeline is a chain of dual quaternion multiplications:

$$
\hat{q}_{\mathrm{total}} = \hat{q}_{\mathrm{coord}} \cdot (1 + \varepsilon\,\delta\hat{q}_{\mathrm{SP}}) \cdot \hat{q}_{\mathrm{Kepler}} \cdot (1 + \varepsilon\,\delta\hat{q}_{\mathrm{LP}}) \cdot \hat{q}_{\mathrm{secular}} \tag{2.33}
$$

The error budget for the entire pipeline follows from applying §2.6 to each multiplication in the chain.

**Downstream chapters that build on this framework:**

| Chapter | Transform expressed as | Key theorem |
|---------|----------------------|-------------|
| Ch. 8 | Perifocal → ECI pure rotation quaternion | Theorem 2.5.1 |
| Ch. 16–17 | Secular perturbation as state update | Definition 2.4.1 |
| Ch. 18 | Short-period corrections as near-identity dual quaternion | Proposition 2.5.1 |
| Ch. 19 | Long-period corrections as near-identity perturbation | Proposition 2.5.1 |
| Ch. 20 | Complete mean → osculating composition | Theorem 2.4.2 |
| Ch. 29 | Sidereal time in TEME → PEF rotation | Theorem 2.5.2 |
| Ch. 30 | All coordinate transforms | Theorems 2.5.1–2.5.3 |
| Ch. 34–35 | Full propagation pipeline | Theorem 2.4.2 (chain) |
| Ch. 36 | Attitude-dependent drag (enhanced mode) | Definition 2.4.1, §2.8 |
| Ch. 38 | Final state vector error budget | §2.6, Definition 2.8.2 |

---

## Error Notes

**[P.2.1]** Trigonometric evaluation in quaternion components. Computing $\cos(\alpha/2)$ and $\sin(\alpha/2)$ from orbital angles introduces representation error bounded by Ch. 1, Corollaries 1.4.1–1.4.2. The half-angle form has smaller arguments than the DCM entries, marginally improving numerical behavior. For double precision: $\delta_p \leq 1.1 \times 10^{-16}$ per trig evaluation.

**[P.2.2]** Quaternion multiplication accumulation. Each quaternion multiplication involves 16 scalar multiplications and 12 additions. Error accumulates as $O(\sqrt{n}\,\epsilon_{\mathrm{mach}})$ per component for a chain of $n$ multiplications.

**[P.2.3]** Renormalization drift. After $n$ quaternion multiplications without renormalization, the unit-norm constraint drifts by $O(n\,\epsilon_{\mathrm{mach}})$. Renormalization (Theorem 2.6.4) restores the constraint at cost of $O(\epsilon_{\mathrm{mach}})$ additional precision error. *Remedy:* renormalize after every multiplication, or after every $k$ multiplications with $k$ chosen so that drift stays below the tolerance.

**[A.2.1]** TEME frame accuracy. The TEME frame differs from the precise GCRS/ITRS frames by $\sim 0.1$ arcsec ($\sim 50$ m at GEO, $\sim 0.3$ m at LEO). This is an accuracy error. *Remedy:* apply IAU 2006/2000A precession-nutation corrections (beyond SGP4 standard).

**[A.2.2]** Linearization in error propagation. The error bounds in §2.6 use first-order (linear) approximations. Neglected second-order cross terms are $O(\delta^2)$. *Remedy:* use full nonlinear propagation from Ch. 1 when $\mathrm{rd}(v) < 4$.

**[M.2.1]** GMST determination. The GMST polynomial coefficients (Ch. 29) carry measurement error from the IAU 1982 determination. Propagates into the transport term: $\sigma_m(\theta_{\mathrm{GMST}}) \approx 0.1$ mas after one day. *Remedy:* use updated IERS values (but breaks the matched pair — see Ch. 3).
