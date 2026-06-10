# Chapter 11: Hamiltonian Mechanics and Delaunay Variables

**Part II: The Two-Body Problem**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $T$ | Kinetic energy | §11.2 |
| $V$ | Potential energy | §11.2 |
| $L$ | Lagrangian: $L = T - V$ | §11.2 |
| $H$ | Hamiltonian: $H = T + V$ | §11.2 |
| $\mathbf{q}$ | Generalized coordinates $(q_1, \ldots, q_n)$ | §11.2 |
| $\mathbf{p}$ | Conjugate momenta $(p_1, \ldots, p_n)$ | §11.2 |
| $p_i$ | Conjugate momentum: $p_i = \partial L / \partial \dot{q}_i$ | §11.2 |
| $\{f, g\}$ | Poisson bracket of $f$ and $g$ | §11.4 |
| $J_{2n}$ | Symplectic matrix of order $2n$ | §11.3 |
| $l$ | Mean anomaly (Delaunay angle variable) | §11.5 |
| $g$ | Argument of perigee (Delaunay angle variable) | §11.5 |
| $h$ | Longitude of ascending node (Delaunay angle variable) | §11.5 |
| $L$ | $L = \sqrt{\mu a}$ (Delaunay action variable) | §11.5 |
| $G$ | $G = L\sqrt{1-e^2} = \sqrt{\mu p}$ (Delaunay action variable) | §11.5 |
| $H$ | $H = G \cos i$ (Delaunay action variable) | §11.5 |
| $H_0$ | Unperturbed (Keplerian) Hamiltonian | §11.6 |
| $a$ | Semi-major axis | Ch 8 |
| $e$ | Orbital eccentricity | Ch 8 |
| $i$ | Orbital inclination | Ch 8 |
| $\eta$ | $\eta = \sqrt{1 - e^2}$ | §11.7 |
| $\theta$ | $\theta = \cos i$ | §11.7 |
| $n$ | Mean motion: $n = \sqrt{\mu / a^3}$ | Ch 8 |
| $p$ | Semi-latus rectum: $p = a(1-e^2)$ | Ch 8 |
| $\mu$ | Gravitational parameter $GM$ | Ch 8 |
| $\beta$ | Auxiliary angle used in perturbation expansions | §11.7 |

---

## §11.1 Introduction

Chapters 8–10 developed the geometry and kinematics of Keplerian motion. Those results are complete for the unperturbed two-body problem: given initial orbital elements, the position and velocity at any later time are determined by Kepler's equation. The next step — accounting for the small forces that cause the orbital elements themselves to evolve — requires a dynamical framework that is both principled and computationally tractable.

This chapter develops that framework. The Hamiltonian formulation of mechanics (§11.2) casts the equations of motion in terms of canonical coordinates and their conjugate momenta, linked by Hamilton's equations. The Legendre transform that connects the Lagrangian and Hamiltonian pictures is developed in full, both as a conceptual bridge to the Newtonian equations already established and as a prerequisite for the transformation theory of §11.3.

Canonical transformations (§11.3) are the central tool of perturbation theory: they are changes of coordinates that preserve the form of Hamilton's equations. The symplectic structure, expressed by the condition $M^T J M = J$ on the Jacobian of the transformation, is the geometric invariant that all canonical transformations share. Poisson brackets (§11.4) encode this structure algebraically and provide the machinery for computing rates of change of any dynamical quantity under a given Hamiltonian.

The Delaunay variables (§11.5) are a set of canonical coordinates tailored to Keplerian orbits. The three angle variables $(l, g, h)$ are the mean anomaly, argument of perigee, and longitude of ascending node; the three conjugate action variables $(L, G, H)$ are functions of the semi-major axis, eccentricity, and inclination. In these variables, the unperturbed Hamiltonian takes its simplest possible form (§11.6), and the secular and periodic parts of any perturbation separate cleanly under orbit-averaging (Chapter 12).

The partial derivative identities of §11.7 are the practical engine of perturbation theory. Every application of Hamilton's equations to a perturbing function requires computing $\partial R / \partial L$, $\partial R / \partial G$, $\partial R / \partial H$ — the rates at which the disturbing function changes with the Delaunay actions. These reduce, via chain rule, to derivatives of the standard orbital elements with respect to $(L, G, H)$, which are algebraic functions of $(a, e, i)$. The full catalog of identities is tabulated here for direct use in Chapters 16–17.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 8, Keplerian elements | §11.5 | Orbital elements $(a, e, i, \Omega, \omega, M)$ used to define Delaunay variables |
| Ch 8, vis-viva relation | §11.6 | Energy integral $H = -\mu/(2a)$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 12, perturbation equations | §11.5 | Canonical framework for perturbation theory |
| Ch 12, von Zeipel method | §11.6 | Generating function structure; unperturbed Hamiltonian |
| Ch 16, first-order secular rates | §11.6, §11.7 | Secular rates and partial derivative identities |
| Ch 17, second-order secular rates | §11.7 | Explicit perturbation rate computations |

---

## §11.2 Hamilton's Equations

This section derives Hamilton's equations from the Legendre transform of the Lagrangian and establishes energy conservation for the Keplerian problem.

**Definition 11.2.1** (Generalized coordinates and momenta)**.** *[stub — generalized coordinates $q_i$, generalized velocities $\dot{q}_i$, and conjugate momenta $p_i = \partial L / \partial \dot{q}_i$]*

**Definition 11.2.2** (Lagrangian)**.** *[stub — $L(\mathbf{q}, \dot{\mathbf{q}}, t) = T - V$; for point mass in gravitational field, $T = \tfrac{1}{2}m|\dot{\mathbf{r}}|^2$, $V = -\mu m/r$]*

**Theorem 11.2.1** (Euler--Lagrange equations)**.** *[stub — the stationarity condition $\delta \int L\,dt = 0$ yields $d/dt(\partial L/\partial \dot{q}_i) - \partial L/\partial q_i = 0$]* — *Proof approach: apply the calculus of variations to Hamilton's principle; integrate by parts with fixed endpoints; invoke the fundamental lemma of the calculus of variations.*

**Definition 11.2.3** (Legendre transform and Hamiltonian)**.** *[stub — $H(\mathbf{q}, \mathbf{p}, t) = \sum_i p_i \dot{q}_i - L$, where $\dot{q}_i$ is eliminated by inverting $p_i = \partial L / \partial \dot{q}_i$]*

**Theorem 11.2.2** (Hamilton's equations)**.** *[stub — from the Legendre transform, $\dot{q}_i = \partial H / \partial p_i$ and $\dot{p}_i = -\partial H / \partial q_i$; also $\partial H / \partial t = -\partial L / \partial t$, so $H$ is conserved when $L$ has no explicit time dependence]* — *Proof approach: compute the total differential $dH = \sum_i \dot{q}_i dp_i - \sum_i (\partial L/\partial q_i) dq_i + (\partial H/\partial t)dt$; identify coefficients of $dp_i$ and $dq_i$ with Hamilton's equations.*

**Theorem 11.2.3** (Energy conservation)**.** *[stub — when the Lagrangian is time-independent, $dH/dt = 0$ along any solution of Hamilton's equations; for the Keplerian problem, $H = T + V = -\mu/(2a)$]* — *Proof approach: compute $dH/dt = \sum_i (\partial H/\partial q_i)\dot{q}_i + \sum_i (\partial H/\partial p_i)\dot{p}_i + \partial H/\partial t$; substitute Hamilton's equations; the first two sums cancel.*

**Remark** (Relation to Newton's equations)**.** *[stub — for the Keplerian problem, Hamilton's equations reduce to $\ddot{\mathbf{r}} = -\mu \mathbf{r}/r^3$; the Hamiltonian formulation adds no physics here but provides the canonical framework needed for §11.3]*

**Example 11.2.1** (Keplerian Hamiltonian in Cartesian coordinates)**.** *[stub — explicit $H = |\mathbf{p}|^2/(2m) - \mu m/|\mathbf{r}|$ for point mass in inverse-square field; verify Hamilton's equations recover $\ddot{\mathbf{r}} = -\mu\mathbf{r}/r^3$]* Values: $\mu = 398600.4418$ km$^3$/s$^2$, $\mathbf{r} = (6778, 0, 0)$ km, $\mathbf{v} = (0, 7.67, 0)$ km/s; compute $H$ and $\dot{q}_i$, $\dot{p}_i$ numerically.

---

## §11.3 Canonical Transformations

This section defines canonical transformations via the symplectic condition, introduces the four types of generating functions, and proves preservation of Hamilton's equations.

**Definition 11.3.1** (Symplectic matrix)**.** *[stub — the $2n \times 2n$ matrix $J_{2n} = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}$; properties: $J^T = -J$, $J^2 = -I$, $\det J = 1$]*

**Definition 11.3.2** (Canonical transformation)**.** *[stub — a diffeomorphism $(\mathbf{q}, \mathbf{p}) \to (\mathbf{Q}, \mathbf{P})$ is canonical if its Jacobian $M$ satisfies $M^T J_{2n} M = J_{2n}$]*

**Theorem 11.3.1** (Preservation of Hamilton's equations)**.** *[stub — if $(\mathbf{Q}, \mathbf{P})$ are related to $(\mathbf{q}, \mathbf{p})$ by a canonical transformation, and $K(\mathbf{Q}, \mathbf{P}, t) = H(\mathbf{q}(\mathbf{Q}, \mathbf{P}), \mathbf{p}(\mathbf{Q}, \mathbf{P}), t) + \partial F / \partial t$ for an appropriate generating function $F$, then $\dot{Q}_i = \partial K / \partial P_i$ and $\dot{P}_i = -\partial K / \partial Q_i$]* — *Proof approach: substitute the generating function relations into the modified Hamilton's principle $\delta\int(\sum_i p_i\dot{q}_i - H)\,dt = \delta\int(\sum_i P_i\dot{Q}_i - K)\,dt + \delta\int dF/dt\,dt$; the boundary terms from $F$ vanish, yielding Hamilton's equations in the new variables.*

**Definition 11.3.3** (Generating functions, types F1–F4)**.** *[stub — four types of generating functions for canonical transformations: $F_1(\mathbf{q}, \mathbf{Q})$, $F_2(\mathbf{q}, \mathbf{P})$, $F_3(\mathbf{p}, \mathbf{Q})$, $F_4(\mathbf{p}, \mathbf{P})$; relations between old and new variables in each case]*

**Theorem 11.3.2** (Invariance of the symplectic form)**.** *[stub — canonical transformations preserve the symplectic 2-form $\omega = \sum_i dp_i \wedge dq_i$; this is the coordinate-independent statement of the symplectic condition]* — *Proof approach: compute the pullback of $\omega$ under the transformation; use the symplectic condition $M^T J M = J$ on the Jacobian to show $\sum_i dP_i \wedge dQ_i = \sum_i dp_i \wedge dq_i$.*

**Example 11.3.1** (Identity transformation via $F_2 = \mathbf{q} \cdot \mathbf{P}$)**.** *[stub — verify that $F_2 = \sum_i q_i P_i$ generates the identity; illustrate the generating function machinery]* Values: for $n = 3$ degrees of freedom, compute $p_i = \partial F_2/\partial q_i = P_i$ and $Q_i = \partial F_2/\partial P_i = q_i$, confirming the identity map.

**Remark** (Point transformations)**.** *[stub — coordinate changes $\mathbf{q} \to \mathbf{Q}(\mathbf{q})$ that depend only on position extend to canonical transformations via the induced momentum transformation; the Delaunay variables arise from a sequence of such transformations]*

---

## §11.4 Poisson Brackets

This section defines the Poisson bracket, proves its algebraic properties, and establishes the equation of motion for arbitrary dynamical quantities.

**Definition 11.4.1** (Poisson bracket)**.** *[stub — for smooth functions $f$, $g$ on phase space, $\{f, g\} = \sum_i \left( \frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i} \right)$]*

**Theorem 11.4.1** (Algebraic properties)**.** *[stub — the Poisson bracket is: (a) bilinear, (b) antisymmetric: $\{f,g\} = -\{g,f\}$, (c) Leibniz: $\{fg, h\} = f\{g,h\} + g\{f,h\}$, (d) Jacobi identity: $\{f,\{g,h\}\} + \{g,\{h,f\}\} + \{h,\{f,g\}\} = 0$]* — *Proof approach: (a)--(c) follow directly from the definition and the product rule for partial derivatives; (d) requires expanding three nested brackets and verifying cancellation of all terms (algebraic identity in second partial derivatives).*

**Theorem 11.4.2** (Fundamental Poisson brackets)**.** *[stub — $\{q_i, q_j\} = 0$, $\{p_i, p_j\} = 0$, $\{q_i, p_j\} = \delta_{ij}$; these are preserved by any canonical transformation (canonical invariant)]* — *Proof approach: direct substitution into the Poisson bracket definition; the fundamental brackets reduce to Kronecker delta by orthogonality of the coordinate partial derivatives.*

**Theorem 11.4.3** (Equations of motion via Poisson brackets)**.** *[stub — for any smooth function $f(\mathbf{q}, \mathbf{p}, t)$, $\dot{f} = \{f, H\} + \partial f/\partial t$; Hamilton's equations are $\dot{q}_i = \{q_i, H\}$, $\dot{p}_i = \{p_i, H\}$]* — *Proof approach: apply the chain rule to $df/dt = \sum_i (\partial f/\partial q_i)\dot{q}_i + \sum_i (\partial f/\partial p_i)\dot{p}_i + \partial f/\partial t$; substitute Hamilton's equations to recover the Poisson bracket.*

**Corollary 11.4.1** (Constants of motion)**.** *[stub — $f$ is a constant of motion if and only if $\{f, H\} + \partial f / \partial t = 0$; for time-independent $f$, the condition is $\{f, H\} = 0$]*

**Theorem 11.4.4** (Invariance under canonical transformations)**.** *[stub — the Poisson bracket $\{f, g\}$ has the same value whether computed in $(\mathbf{q}, \mathbf{p})$ or any canonically equivalent $(\mathbf{Q}, \mathbf{P})$ coordinates; this invariance characterizes canonical transformations]* — *Proof approach: express the Poisson bracket in matrix form as $\{f,g\} = (\nabla f)^T J (\nabla g)$; under the canonical transformation with Jacobian $M$, the bracket becomes $(\nabla' f)^T M^T J M (\nabla' g) = (\nabla' f)^T J (\nabla' g)$ by the symplectic condition.*

**Example 11.4.1** (Angular momentum Poisson brackets)**.** *[stub — verify $\{L_x, L_y\} = L_z$ cyclically for the orbital angular momentum vector $\mathbf{L} = \mathbf{r} \times \mathbf{p}$]* Values: compute each bracket by expanding $L_x = yp_z - zp_y$, $L_y = zp_x - xp_z$, $L_z = xp_y - yp_x$ in the Poisson bracket definition. Numerical spot-check at $\mathbf{r} = (6778, 0, 0)$ km, $\mathbf{p} = (0, 7.67, 1.0)$ km/s.

---

## §11.5 Delaunay Variables

This section defines the Delaunay action-angle variables, proves they are canonical, and derives their algebraic relations to the standard Keplerian elements.

**Definition 11.5.1** (Delaunay action variables)**.** *[stub — $L = \sqrt{\mu a}$, $G = \sqrt{\mu a(1-e^2)} = \sqrt{\mu p}$, $H = \sqrt{\mu a(1-e^2)} \cos i = G \cos i$]*

**Definition 11.5.2** (Delaunay angle variables)**.** *[stub — $l = M$ (mean anomaly), $g = \omega$ (argument of perigee), $h = \Omega$ (longitude of ascending node)]*

**Theorem 11.5.1** (Canonical brackets for Delaunay variables)**.** *[stub — the Delaunay variables satisfy $\{l, L\} = 1$, $\{g, G\} = 1$, $\{h, H\} = 1$, with all other pairs having zero Poisson bracket; i.e., $(l, g, h, L, G, H)$ is a canonical coordinate set]* — *Proof approach: compute the Jacobian of the transformation from Cartesian phase-space variables to Delaunay variables; verify the symplectic condition $M^T J M = J$ directly, or equivalently verify the 15 independent Poisson bracket relations by chain-rule expansion.*

**Theorem 11.5.2** (Inversion: orbital elements from Delaunay variables)**.** *[stub — explicit algebraic relations: $a = L^2/\mu$, $e = \sqrt{1 - G^2/L^2}$, $i = \arccos(H/G)$; $M = l$, $\omega = g$, $\Omega = h$]* — *Proof approach: algebraic inversion of the definitions in Def. 11.5.1; solve $L = \sqrt{\mu a}$ for $a$, eliminate $G$ to get $e$, and use $H = G\cos i$ to get $i$.*

**Remark** (Singularities)**.** *[stub — the Delaunay angle $g$ is singular when $e = 0$ (circular orbit: $G = L$) and $h$ is singular when $i = 0$ (equatorial orbit: $H = G$); the Lyddane modification (Ch 10, Ch 18) addresses these singularities in the Brouwer theory]*

**Remark** (Physical interpretation of the actions)**.** *[stub — $L$ encodes the semi-major axis (energy), $G$ encodes the magnitude of specific angular momentum, $H$ encodes the $z$-component of specific angular momentum; these are the three classical conserved quantities of Keplerian motion]*

**Example 11.5.1** (Delaunay variables for a representative LEO orbit)**.** Values: $a = 6878.14$ km, $e = 0.001$, $i = 51.6°$, $\mu = 398600.4418$ km$^3$/s$^2$. Compute $L = \sqrt{\mu a}$, $G = L\sqrt{1-e^2}$, $H = G\cos i$ to full double precision. Verify round-trip inversion recovers $a$, $e$, $i$.

---

## §11.6 The Unperturbed Hamiltonian

This section derives the Keplerian Hamiltonian $H_0 = -\mu^2/(2L^2)$ in Delaunay variables and shows that unperturbed motion reduces to linear advance of the mean anomaly.

**Theorem 11.6.1** (Keplerian Hamiltonian in Delaunay variables)**.** *[stub — $H_0 = -\mu^2 / (2L^2)$]* — *Proof approach: substitute $a = L^2/\mu$ into $H_0 = -\mu/(2a)$; simplify to $H_0 = -\mu^2/(2L^2)$; verify independence from $G$, $H$, $l$, $g$, $h$.*

**Corollary 11.6.1** (Unperturbed equations of motion)**.** *[stub — Hamilton's equations give $\dot{l} = \partial H_0 / \partial L = \mu^2/L^3 = n$, $\dot{g} = \dot{h} = 0$, $\dot{L} = \dot{G} = \dot{H} = 0$; unperturbed Keplerian motion is $l = l_0 + n(t - t_0)$ with all other Delaunay variables constant]* — *Proof approach: apply Hamilton's equations $\dot{l} = \partial H_0/\partial L$, $\dot{L} = -\partial H_0/\partial l = 0$ (since $H_0$ is independent of $l$); analogously for $(g,G)$ and $(h,H)$.*

**Remark** (Degeneracy)**.** *[stub — $H_0$ depends only on $L$, making the unperturbed system degenerate (or anisochronous): there is one non-trivial frequency $n$ and five zero frequencies; this degeneracy means that $g$ and $h$ are undetermined by the unperturbed equations and are fixed by initial conditions; perturbation theory must carefully handle this near-degeneracy]*

**Remark** (Why Delaunay variables are natural for perturbation theory)**.** *[stub — in Delaunay variables, the perturbation $R$ appears only on the right-hand side of Hamilton's equations; the unperturbed part $H_0$ is isolated and trivially integrated; the structure $H = H_0 + R$ directly motivates the averaging and von Zeipel methods of Chapter 12]*

---

## §11.7 Partial Derivative Identities

This section derives and tabulates all first partial derivatives of $(a, e, i, \eta, \theta, n, p)$ with respect to $(L, G, H)$ for use in the perturbation equations of Chs 16--17.

**Definition 11.7.1** (Standard derived quantities)**.** *[stub — $\eta = \sqrt{1-e^2} = G/L$, $\theta = \cos i = H/G$, $p = a(1-e^2) = G^2/\mu$, $n = \mu^2/L^3$; all are algebraic functions of the Delaunay actions]*

**Theorem 11.7.1** (First partial derivatives with respect to $L$)**.** *[stub — from $a = L^2/\mu$, $n = \mu^2/L^3$, $p = G^2/\mu$, and $e^2 = 1 - G^2/L^2$:*

$$\frac{\partial a}{\partial L} = \frac{2L}{\mu} = \frac{2a}{L}, \quad
\frac{\partial n}{\partial L} = -\frac{3\mu^2}{L^4} = -\frac{3n}{L}, \quad
\frac{\partial e}{\partial L} = \frac{G^2}{L^3 e} = \frac{\eta^2}{Le}, \quad
\frac{\partial p}{\partial L} = 0 \quad (\text{stub — verify sign})$$

*full table: all first-order derivatives $\partial/\partial L$, $\partial/\partial G$, $\partial/\partial H$ of $a, e, i, \eta, \theta, n, p$]* — *Proof approach: differentiate $a = L^2/\mu$ directly; for $e$, differentiate $e^2 = 1 - G^2/L^2$ implicitly to get $2e\,\partial e/\partial L = 2G^2/L^3$; for $\eta = G/L$, apply the quotient rule.*

**Theorem 11.7.2** (First partial derivatives with respect to $G$)**.** *[stub — from $G = L\eta = \sqrt{\mu p}$, $i = \arccos(H/G)$: $\partial a/\partial G = 0$, $\partial e/\partial G = -G/(L^2 e)= -\eta/(Le)$, $\partial p/\partial G = 2G/\mu$, $\partial i/\partial G = H/(G^2\sqrt{1-(H/G)^2}) = -\theta/(G\sin i)$; full table for all quantities]* — *Proof approach: $a = L^2/\mu$ is independent of $G$; differentiate $e^2 = 1 - G^2/L^2$ to get $\partial e/\partial G$; differentiate $\cos i = H/G$ implicitly with $H$ held fixed.*

**Theorem 11.7.3** (First partial derivatives with respect to $H$)**.** *[stub — from $\theta = H/G$, $i = \arccos(H/G)$: $\partial a/\partial H = \partial e/\partial H = \partial p/\partial H = 0$, $\partial i/\partial H = -1/(G\sin i)$, $\partial \theta/\partial H = 1/G$; full table]* — *Proof approach: $a$, $e$, $p$ depend only on $L$ and $G$, not on $H$; differentiate $\cos i = H/G$ with $G$ held fixed to get $\partial i/\partial H$.*

**Theorem 11.7.4** (Chain rule identities for perturbation calculations)**.** *[stub — for a disturbing function $R(a, e, i)$, express $\partial R / \partial L$, $\partial R / \partial G$, $\partial R / \partial H$ via the chain rule and the tables of Theorems 11.7.1--11.7.3; these are the quantities that enter Hamilton's equations in Ch 16--17]* — *Proof approach: apply the multivariate chain rule $\partial R/\partial L = (\partial R/\partial a)(\partial a/\partial L) + (\partial R/\partial e)(\partial e/\partial L) + (\partial R/\partial i)(\partial i/\partial L)$; substitute the derivative table entries.*

**Table 11.7.1** (Complete first-derivative table)**.** *[stub — tabulated summary of all $\partial/\partial L$, $\partial/\partial G$, $\partial/\partial H$ for $(a, e, i, \eta, \theta, n, p)$; the definitive reference for Chapters 16 and 17]*

**Example 11.7.1** (Verification for a representative orbit)**.** Values: $a = 7000$ km, $e = 0.1$, $i = 45°$, $\mu = 398600.4418$ km$^3$/s$^2$. Compute $L$, $G$, $H$; evaluate all entries in Table 11.7.1 analytically. Cross-check $\partial n/\partial L = -3n/L$ and $\partial e/\partial G = -\eta/(Le)$ against forward-difference numerical derivatives (step size $\Delta L = 10^{-6} L$).

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.11.1] | M | §11.5 | Gravitational parameter $\mu$ is a measured quantity with precision limited by the geodetic system (WGS72 or WGS84). All Delaunay action variables inherit this measurement uncertainty, since $L = \sqrt{\mu a}$, $G = \sqrt{\mu p}$. |
| [P.11.1] | P | §11.5 | Inversion $e = \sqrt{1 - G^2/L^2}$ is susceptible to subtractive cancellation when $e \ll 1$; use $e = \sqrt{(L-G)(L+G)}/L$ or Lyddane variables. |
| [P.11.2] | P | §11.5 | Inclination $i = \arccos(H/G)$ is ill-conditioned near $i = 0$ and $i = \pi$; use Lyddane modification or equatorial regularization. |
| [A.11.1] | A | §11.2 | Hamiltonian framework assumes point mass in a conservative field; drag and solar radiation pressure are non-conservative. |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 10 |
| Theorems | 16 |
| Lemmas | 0 |
| Corollaries | 2 |
| Propositions | 0 |
| Examples | 5 |
| Error Notes | 4 |
| Equations | ~30 |
| Sections | 7 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §11.1 Introduction | Draft | |
| §11.2 Hamilton's Equations | Draft | |
| §11.3 Canonical Transformations | Draft | |
| §11.4 Poisson Brackets | Draft | |
| §11.5 Delaunay Variables | Draft | |
| §11.6 The Unperturbed Hamiltonian | Draft | |
| §11.7 Partial Derivative Identities | Draft | |
