# Draft Plan: Chapter 9 — Kepler's Equation

**Part II: The Two-Body Problem**

**Phase:** Draft

**Implementation target:** `kepler.h`

---

## Objectives

1. Derive Kepler's equation $M = E - e\sin E$ from the geometry of the eccentric anomaly construction.
2. Establish the exact relations among true anomaly $\nu$, eccentric anomaly $E$, and mean anomaly $M$ (half-angle formulae and their inverses).
3. Prove existence and uniqueness of $E$ for given $M$ and $e \in [0,1)$ by the contraction mapping theorem.
4. Derive and certify Newton's method (quadratic convergence), Halley's method (cubic convergence), and Householder's method (fourth-order convergence) for solving Kepler's equation.
5. Derive the continued fraction solution as an alternative with guaranteed convergence.
6. Characterize the near-parabolic behavior as $e \to 1^-$ and identify the failure modes.
7. Specify starting value strategies with quantified convergence bounds.
8. Map every computation to the error-tracking framework of Ch 1.

---

## Notation Table

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $M$ | Mean anomaly | Ch 8, §8.5 |
| $E$ | Eccentric anomaly | §9.2 |
| $\nu$ | True anomaly | Ch 8, §8.4 |
| $e$ | Eccentricity | Ch 8, §8.4 |
| $a$ | Semi-major axis | Ch 8, §8.5 |
| $b$ | Semi-minor axis: $b = a\sqrt{1-e^2} = a\eta$ | §9.2 |
| $\eta$ | Eccentricity complement: $\eta = \sqrt{1-e^2}$ | §9.2 |
| $f(E)$ | Kepler residual: $f(E) = E - e\sin E - M$ | §9.5 |
| $E_0$ | Starting value for iteration | §9.9 |
| $E_k$ | Iterate at step $k$ | §9.5 |
| $\epsilon_k$ | Error at step $k$: $\epsilon_k = E_k - E_*$ | §9.5 |
| $\tau$ | Caller-supplied tolerance (Ch 3, Def. 3.4.1) | Ch 3 |
| $\kappa$ | Contraction constant for Newton's method near convergence | §9.5 |
| $\mathcal{V}$ | Tracked value quadruple $(v, \sigma_m, \delta_p, \delta_a)$ | Ch 1 |

---

## Section Structure

### §9.1 Introduction

Narrative: Kepler's equation is the fundamental transcendental equation of orbit propagation — given mean anomaly $M$ (time), find eccentric anomaly $E$ (geometry), from which position and velocity follow. It is solved once per propagation step and must be solved accurately and efficiently. The equation has no closed-form algebraic solution; numerical methods are essential.

Overview of the iterative hierarchy: Newton (quadratic) suffices for most orbits; Halley (cubic) and Householder (fourth-order) reduce iteration count at high eccentricity; continued fractions provide an alternative with different convergence character. The near-parabolic regime ($e$ near 1) requires special treatment.

Forward reference table:

| Section | Feeds | Role |
|---------|-------|------|
| §9.3 ($\nu$–$E$ relations) | Ch 10 | Modified Kepler uses the same half-angle forms |
| §9.5 (Newton's method) | Ch 10, Ch 25 | Same iteration pattern used for equation of center |
| §9.4 (existence/uniqueness) | Ch 10 | Contraction mapping applies to Lyddane form |
| §9.8 (continued fraction) | Ch 4 | Application of Ch 4 continued fraction theory |

---

### §9.2 Derivation from Eccentric Anomaly Geometry

This section defines the eccentric anomaly geometrically and derives Kepler's equation $M = E - e\sin E$ from the area swept by the radius vector.

- **Definition 9.2.1** (Eccentric anomaly). Geometric construction: circumscribed circle of radius $a$; eccentric anomaly $E$ as the angle at the center between the auxiliary circle point and the reference direction.
- **Definition 9.2.2** (Cartesian coordinates via $E$). $x = a\cos E$, $y = b\sin E = a\eta\sin E$; position on ellipse from auxiliary point.
- **Theorem 9.2.1** (Kepler's equation). $M = E - e\sin E$. — *Proof approach: parametrize the ellipse sector area via the auxiliary circle; the ellipse area from periapsis to point $P$ is $(b/a)$ times the corresponding circular sector area minus triangle correction; equate to $hT/(2\pi)\cdot M$ via Kepler's second law.*
- **Corollary 9.2.1** (Velocity from $E$). $\dot{x} = -(na/r)\sin E$, $\dot{y} = (na\eta/r)\cos E$ where $r = a(1 - e\cos E)$. — *Proof approach: differentiate the Cartesian parametrization with respect to time; use $dE/dt = n\cdot a/r$ from Kepler's equation.*
- **Corollary 9.2.2** (Radius from $E$). $r = a(1 - e\cos E)$. — *Proof approach: compute $r = \sqrt{x^2 + y^2}$ from the Cartesian parametrization $x = a\cos E - ae$, $y = b\sin E$; simplify.*
- **Example 9.2.1** (Geometric verification). Verify $M = E - e\sin E$ at $E = 0$, $E = \pi/2$, $E = \pi$ for $e = 0.25$. Compute $M$ at each point and verify geometric consistency (periapsis, quarter-orbit, apoapsis).
- *Error Note placeholder:* [P.9.1] the area computation is exact; no precision error in the derivation itself.

---

### §9.3 Anomaly Relations

This section establishes the exact algebraic relations among $\nu$, $E$, and $M$ in both half-angle and sine/cosine forms.

- **Theorem 9.3.1** (True-to-eccentric anomaly half-angle formulae). $\tan(E/2) = \sqrt{(1-e)/(1+e)}\tan(\nu/2)$. — *Proof approach: eliminate $r$ between $r = a(1 - e\cos E)$ and $r = p/(1 + e\cos\nu)$; apply half-angle identities to simplify the resulting relation between $\cos E$ and $\cos\nu$.*
- **Corollary 9.3.1** (Eccentric-to-true). $\tan(\nu/2) = \sqrt{(1+e)/(1-e)}\tan(E/2)$.
- **Theorem 9.3.2** (Sine and cosine forms). $\sin\nu = \eta\sin E/(1 - e\cos E)$; $\cos\nu = (\cos E - e)/(1 - e\cos E)$. — *Proof approach: express $\cos\nu$ by eliminating $r$ between orbit equation and radius-from-$E$ formula; derive $\sin\nu$ from $\sin^2\nu + \cos^2\nu = 1$ using $\eta^2 = 1-e^2$.*
- **Corollary 9.3.2** (Inverse forms). $\sin E = \eta\sin\nu/(1 + e\cos\nu)$; $\cos E = (e + \cos\nu)/(1 + e\cos\nu)$.
- **Theorem 9.3.3** (Time from eccentric anomaly). $t - t_0 = (M - M_0)/n = (E - e\sin E - M_0)/n$.
- **Remark** (atan2 evaluation). All four quadrant-resolving conversions use atan2 (Ch 7) to avoid the quadrant ambiguity of the half-angle form.
- **Example 9.3.1** (Numerical anomaly chain). Given $M = 1.2$ rad, $e = 0.15$: solve for $E$ via Newton iteration (§9.5), then compute $\nu$ via half-angle formula. Verify the round-trip $M(E(\nu)) = M$ to machine precision. Tabulate $M$, $E$, $\nu$, and $r$.
- *Error Note placeholder:* [P.9.2] precision loss in $\cos\nu = (\cos E - e)/(1 - e\cos E)$ near $E = 0$ when $e \approx 1$ (cancellation in numerator); alternative form.

---

### §9.4 Existence and Uniqueness

This section proves existence and uniqueness of the solution to Kepler's equation and establishes the contraction mapping property for Newton's method.

- **Theorem 9.4.1** (Existence). For all $M \in \mathbb{R}$ and $e \in [0,1)$, Kepler's equation $E - e\sin E = M$ has at least one solution $E_*$. — *Proof approach: apply the intermediate value theorem to $f(E) = E - e\sin E - M$ on $[M-\pi, M+\pi]$; show $f$ changes sign.*
- **Theorem 9.4.2** (Uniqueness). For $e \in [0,1)$, the solution $E_*$ is unique. — *Proof approach: compute $f'(E) = 1 - e\cos E \geq 1 - e > 0$; strict monotonicity of $f$ implies at most one zero.*
- **Corollary 9.4.1** (Uniform boundedness). $|E - M| \leq e/(1-e)$ for $e \in [0,1)$. — *Proof approach: bound $|f(M)| = |e\sin M| \leq e$; use $f'(E) \geq 1-e$ to convert residual bound to root bound via the mean value theorem.*
- **Theorem 9.4.3** (Contraction for Newton's method). Define $g(E) = E - f(E)/f'(E)$. For $e \leq e_{\max} < 1$, the Newton map $g$ is a contraction on $[E_* - \delta, E_* + \delta]$ for sufficiently small $\delta$. — *Proof approach: compute $g'(E) = f(E)f''(E)/(f'(E))^2$; evaluate $g'(E_*) = 0$; bound $|g''|$ on a neighborhood to establish the contraction radius via the Banach fixed-point theorem.*
- *Error Note placeholder:* [A.9.1] the elliptic assumption $e < 1$ is strict; behavior at $e = 1$ (parabola) addressed in §9.10.

---

### §9.5 Newton's Method

This section derives the Newton iteration formula for Kepler's equation and proves quadratic convergence with explicit bounds.

- **Definition 9.5.1** (Newton iteration for Kepler's equation). $E_{k+1} = E_k - f(E_k)/f'(E_k) = E_k - (E_k - e\sin E_k - M)/(1 - e\cos E_k)$.
- **Theorem 9.5.1** (Quadratic convergence). The Newton iteration satisfies $|\epsilon_{k+1}| \leq C|\epsilon_k|^2$ where $C = |f''(E_*)| / (2|f'(E_*)|) = e|\sin E_*| / (2(1 - e\cos E_*)^2)$. — *Proof approach: Taylor expand $f(E_k)$ about $E_*$ to second order; use $f(E_*) = 0$ and the Newton update formula to isolate the error recurrence.*
- **Corollary 9.5.1** (Digit doubling). Starting from $|\epsilon_0| < 1$, after $k$ Newton steps $|\epsilon_k| \lesssim C^{2^k - 1}|\epsilon_0|^{2^k}$. Implication: 4–5 iterations reach double-precision convergence from any reasonable starting value.
- **Proposition 9.5.1** (Refinement step). After convergence test $|E_{k+1} - E_k| < \tau$, one additional iteration reduces the residual to $O(\tau^2/|f'|)$; the excess is added to $\delta_p$.
- **Example 9.5.1** (Newton convergence table). Tabulate $E_k$, $f(E_k)$, $|\epsilon_k|$ for $M = 2.5$ rad, $e = 0.65$, starting from $E_0 = M$. Show quadratic convergence by computing the ratio $|\epsilon_{k+1}|/|\epsilon_k|^2$ at each step. Expect convergence in 5--6 iterations to double precision.
- **Remark** (Tolerance parameter). At $\tau = \tau_{\text{standard}}$ (Ch 3), the iteration reproduces the SGP4 standard truncation. At $\tau < \tau_{\text{standard}}$, additional iterations are performed.
- *Error Note placeholder:* [P.9.3] division by $f'(E_k) = 1 - e\cos E_k$ becomes small when $e \cos E_k \approx 1$, i.e., near $E \approx 0$ at high eccentricity; bound on the division amplification.

---

### §9.6 Halley's Method

This section derives Halley's cubic-convergence iteration and quantifies the cost-benefit tradeoff relative to Newton.

- **Definition 9.6.1** (Halley iteration). $E_{k+1} = E_k - \frac{f(E_k)/f'(E_k)}{1 - \tfrac{1}{2}f(E_k)f''(E_k)/(f'(E_k))^2}$; equivalently $E_{k+1} = E_k - \frac{f \cdot f'}{(f')^2 - \tfrac{1}{2}f \cdot f''}$ where $f'' = e\sin E$.
- **Theorem 9.6.1** (Cubic convergence). The Halley iteration satisfies $|\epsilon_{k+1}| \leq K|\epsilon_k|^3$. — *Proof approach: Taylor expand $f$ to third order about $E_*$; show that the Halley denominator correction cancels the $O(\epsilon^2)$ term exactly, leaving a cubic error recurrence.*
- **Corollary 9.6.1** (Comparison with Newton). Halley requires one additional evaluation of $\sin E$ per step but halves the number of steps needed. For $e > 0.5$, Halley saves at least one iteration relative to Newton.
- **Example 9.6.1** (Halley convergence table). Same parameters as Example 9.5.1 ($M = 2.5$ rad, $e = 0.65$, $E_0 = M$). Tabulate $E_k$, $|f(E_k)|$, step count to double precision. Compare iteration count with Newton (expect 1--2 fewer steps).
- *Error Note placeholder:* [P.9.4] the denominator $1 - f f''/(2(f')^2)$ can be computed with one extra multiply; precision budget.

---

### §9.7 Householder's Method

This section derives the fourth-order Householder iteration and identifies the practical eccentricity crossover where it outperforms Newton and Halley.

- **Definition 9.7.1** (Householder iteration of order 3). Fourth-order method using $f$, $f'$, $f''$, $f'''$ where $f''' = e\cos E$ (equals $-d^3f/dE^3$ up to sign). Explicit formula.
- **Theorem 9.7.1** (Fourth-order convergence). $|\epsilon_{k+1}| \leq H|\epsilon_k|^4$. — *Proof approach: the general Householder method of order $d$ eliminates Taylor terms through order $d$, yielding $(d+1)$-th order convergence; specialize to $d = 3$ with the explicit derivatives $f' = 1 - e\cos E$, $f'' = e\sin E$, $f''' = e\cos E$.*
- **Corollary 9.7.1** (Practical benefit). For $e > 0.8$, Householder reaches double-precision in 2–3 iterations from a good starting value.
- **Remark** (Diminishing returns). Higher-order Householder methods ($d \geq 4$) exist but require additional derivatives and provide diminishing returns; $d = 3$ is the practical optimum for double precision.
- *Error Note placeholder:* [P.9.5] additional arithmetic cost of Householder vs. benefit; cross-over eccentricity where it outperforms Newton.

---

### §9.8 Continued Fraction Solution

This section derives the Lagrange inversion series and its continued fraction acceleration for solving Kepler's equation, with convergence radius analysis.

- **Theorem 9.8.1** (Lagrange inversion). $E$ can be expressed as a power series in $e$ via the Lagrange inversion theorem applied to $E = M + e\sin E$: $E = M + \sum_{k=1}^{\infty} c_k(M) e^k$ where the $c_k$ are explicit trigonometric polynomials in $M$. — *Proof approach: apply the Lagrange inversion formula for implicit functions $z = w + \epsilon\phi(z)$; compute successive $c_k$ by the Faa di Bruno differentiation formula.*
- **Definition 9.8.1** (Kapteyn series). The formal power series in $e$ for $E(M, e)$ via Bessel functions $J_k$; $E = M + 2\sum_{k=1}^{\infty} J_k(ke)\sin(kM)/k$.
- **Theorem 9.8.2** (Continued fraction form). A convergent continued fraction for $E - M$ in terms of $e$ and $M$; derivation via the Euler continued fraction method (Ch 4, §4.3). — *Proof approach: apply the Euler continued fraction transformation to the Lagrange power series; verify termwise equivalence and establish the extended convergence domain.*
- **Corollary 9.8.1** (Convergence radius). The Lagrange series converges for $|e| < e_{\text{Laplace}} \approx 0.6627$; outside this radius the continued fraction extension is required.
- **Generalization.** The continued fraction representation provides convergence for $e$ values where the power series diverges; the method from Ch 4 applies directly. Connection to Ch 4, §4.3 (Euler continued fractions).
- *Error Note placeholder:* [A.9.2] Bessel function evaluation cost vs. iterative methods; the continued fraction is competitive for moderate eccentricity.

---

### §9.9 Starting Values

This section derives and compares starting value strategies with quantified convergence bounds for each eccentricity regime.

- **Proposition 9.9.1** (Simple starting value). $E_0 = M$ gives $|\epsilon_0| \leq e$ and requires at most $\lceil \log_2(\log_2(e\tau^{-1})) \rceil + 1$ Newton steps. — *Proof approach: bound $|f(M)| = |e\sin M| \leq e$; apply the quadratic convergence bound from Theorem 9.5.1 recursively to determine step count.*
- **Proposition 9.9.2** (Improved starting value for high eccentricity). $E_0 = M + e\sin M/(1 - \sin(M+e) + \sin M)$ (Markley 1995); reduces initial residual by approximately one Newton step for $e > 0.5$.
- **Proposition 9.9.3** (Near-$\pi$ starting value). For $M$ near $\pi$ and moderate $e$, $E_0 = \pi$ reduces the initial residual; characterize the convergence neighborhood.
- **Theorem 9.9.1** (Universal starting value guarantee). For any $e \in [0, 1-\delta]$ and any $M$, the iteration starting from $E_0 = M$ converges to double-precision in at most $N_{\max}(e, \delta)$ Newton steps. Explicit formula for $N_{\max}$.
- **Example 9.9.1** (Starting value comparison). Tabulate iteration counts vs. $e$ for $M = 1.0$ rad at $e \in \{0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95\}$, comparing the three starting value strategies ($E_0 = M$, Markley, near-$\pi$). Count iterations to $|f(E_k)| < 10^{-15}$.

---

### §9.10 Near-Parabolic Behavior

This section characterizes the failure of Newton's method as $e \to 1^-$ and identifies the practical failure threshold for SGP4-relevant orbits.

- **Theorem 9.10.1** (Divergence of Newton's method near $e = 1$). The Newton contraction constant $C = e|\sin E_*|/(2(1-e\cos E_*)^2)$ diverges as $e \to 1^-$ and $E_* \to 0$ simultaneously (perigee passage at high eccentricity). — *Proof approach: evaluate $C$ along the curve $E_* = 2\arcsin\sqrt{1-e}$ (the locus of perigee passage); show $C \sim 1/(4(1-e))$ as $e \to 1^-$.*
- **Proposition 9.10.1** (Practical failure threshold). For $e > e_{\text{crit}} \approx 0.99$ and $M$ near $0$ or $2\pi$, standard Newton iteration requires more than 20 steps to converge; a reformulated equation is needed.
- **Remark** (SGP4 scope). SGP4 is defined only for $e < 1$, but operational TLEs can have eccentricities up to $e \approx 0.999$ for highly elliptic orbits (Molniya-type). The failure threshold is relevant for the implementation.
- **Remark** (Alternative: universal variable). For $e$ near 1, the universal variable formulation (Battin 1987) avoids the singularity; this is noted as a generalization path but not derived here.
- **Generalization.** Universal variables and Stumpff functions provide a unified framework valid for $e \geq 1$; the tolerance parameter approach allows the library to switch formulation at a caller-specified threshold.
- *Error Note placeholder:* [A.9.3] near-parabolic threshold is context-dependent; the accuracy error of truncating to the elliptic form.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Thm 1.4.1 (convergence error) | Error tracking, convergence error as precision error, subtractive cancellation warning |
| Ch 3 | Thm 3.4.1 (tolerance parameter) | Tolerance parameter and tau_standard; matched-pair principle |
| Ch 4 | Thm 4.3 (continued fractions) | Continued fractions, Lagrange inversion, convergence analysis tools |
| Ch 5 | Geometric-tail series evaluator | Geometric-tail series evaluator for the Kapteyn series |
| Ch 7 | atan2 quadrant resolution | atan2 for quadrant resolution in anomaly conversion |
| Ch 8 | Orbital element definitions | All orbital element definitions; solves for E given elements from Ch 8 |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 10 | Modified Kepler iteration machinery | Modified Kepler equation uses the same Newton/Halley machinery on a different f |
| Ch 11 | Delaunay element evaluation | Mean-to-eccentric anomaly conversion is the first step in Delaunay element evaluation |
| Ch 20 | Osculating element pipeline | The Kepler step from mean to eccentric anomaly is the core of the osculating element pipeline |
| Ch 25 | Equation-of-center expansion | Equation-of-center expansion uses the same small-e series |
| Appendix B | kepler.h implementation | Implementation specifications |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.9.1] | P | §9.2 | Derivation of Kepler's equation via area is exact; precision error arises from finite-precision $\sin E$ and $\cos E$ in iteration |
| [P.9.2] | P | §9.3 | Cancellation in $\cos\nu = (\cos E - e)/(1 - e\cos E)$ when $e \approx 1$ and $E \approx 0$; alternative half-angle form avoids cancellation |
| [P.9.3] | P | §9.5 | Newton denominator $1 - e\cos E$ bounded below by $1 - e$; divisor amplification factor $1/(1-e)$, genuine precision amplification |
| [P.9.4] | P | §9.6 | Halley denominator precision; correction to Newton step is $O(\tau^2 e / (f')^3)$, negligible once converged |
| [P.9.5] | P | §9.7 | Householder precision budget; three extra multiplies and one divide per step; cross-over near $e \approx 0.8$ |
| [A.9.1] | A | §9.4 | Elliptic case $e < 1$ assumed; SGP4 always elliptic |
| [A.9.2] | A | §9.8 | Kapteyn/Bessel series computationally inferior to Newton; main use is error bounding and continued fraction derivation |
| [A.9.3] | A | §9.10 | Near-parabolic failure threshold $e_{\rm crit} \approx 0.99$; Molniya orbits at $e \approx 0.74$ safely below |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 6 |
| Theorems | 14 |
| Lemmas | 0 |
| Corollaries | 9 |
| Propositions | 5 |
| Examples | 5 |
| Error Notes | 8 |
| Equations | ~30 |
| Sections | 10 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §9.1 Introduction | Draft | |
| §9.2 Derivation from Eccentric Anomaly Geometry | Draft | |
| §9.3 Anomaly Relations | Draft | |
| §9.4 Existence and Uniqueness | Draft | |
| §9.5 Newton's Method | Draft | |
| §9.6 Halley's Method | Draft | |
| §9.7 Householder's Method | Draft | |
| §9.8 Continued Fraction Solution | Draft | |
| §9.9 Starting Values | Draft | |
| §9.10 Near-Parabolic Behavior | Draft | |
