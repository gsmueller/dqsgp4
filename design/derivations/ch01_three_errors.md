# Chapter 1: The Three Fundamental Errors

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathbb{F}$ | The scalar field: $\mathbb{R}$ (real) or $\mathbb{C}$ (complex). All results in §§1.3–1.4 hold for both unless noted. | §1.3 |
| $|z|$ | Scalar absolute value ($z \in \mathbb{R}$) or modulus ($z \in \mathbb{C}$): $|a + bi| = \sqrt{a^2 + b^2}$ | §1.3 |
| $\|\mathbf{v}\|$ | Euclidean norm of a vector: $\|\mathbf{v}\| = \sqrt{\sum v_i^2}$. Reserved for vectors and matrices (Ch 2). Not used for scalars in Ch 1. | Ch 2 |
| $v$ | A computed or measured scalar value (real or complex) | §1.2 |
| $\tilde{v}$ | The true (unknowable) value that $v$ approximates | §1.2 |
| $\hat{v}$ | The result of computationally evaluating a mathematical expression for $v$ | §1.2, Def. 1.2.2 |
| $v_{\mathrm{model}}$ | The exact (infinite-precision) prediction of the mathematical model | §1.2, Def. 1.2.3 |
| $\sigma_m$ | Measurement error bound (non-negative) | §1.2, Def. 1.2.1 |
| $\delta_p$ | Precision error bound (non-negative) | §1.2, Def. 1.2.2 |
| $\delta_a$ | Accuracy (model) error bound (non-negative) | §1.2, Def. 1.2.3 |
| $\delta$ | A generic error bound (any single category), used when a result applies uniformly to all three | §1.3 |
| $\delta_{\mathrm{total}}$ | Total error bound $= \sigma_m + \delta_p + \delta_a$ | §1.8, Def. 1.8.1 |
| $\delta_{\mathrm{RSS}}$ | Root-sum-square error bound $= \sqrt{\sigma_m^2 + \delta_p^2 + \delta_a^2}$ | §1.8, Def. 1.8.2 |
| $\mathcal{V}$ | A tracked value: the quadruple $(v, \sigma_m, \delta_p, \delta_a)$ | §1.2, Def. 1.2.4 |
| $\mathrm{fl}(v)$ | The IEEE 754 floating-point representation of real number $v$ | §1.2 |
| $\mathrm{ulp}(v)$ | Unit in the last place: spacing between adjacent representable floats near $v$ | §1.2 |
| $\epsilon_{\mathrm{mach}}$ | Machine epsilon: $2^{-52} \approx 1.1 \times 10^{-16}$ for binary64 | §1.2 |
| $\epsilon$ | Convergence tolerance for an iterative solver | §1.6 |
| $L$ | Lipschitz (contraction) constant satisfying $|g(x) - g(y)| \leq L|x-y|$, with $L < 1$ | §1.6, Asm. 1.6.1 |
| $\kappa$ | Condition number of an operation (ratio of relative output error to relative input error) | §1.7, Thm. 1.7.1 |
| $\mathrm{rd}(v)$ | Number of reliable decimal digits in $v$ | §1.9, Def. 1.9.1 |
| $R_N$ | Remainder (tail) of a series truncated after $N$ terms: $R_N = \sum_{k=N+1}^{\infty} a_k$ | §1.5 |
| $S_N$ | Partial sum of a series: $S_N = \sum_{k=0}^{N} a_k$ | §1.5 |
| $\nabla f$ | Gradient vector $(\partial f/\partial x_1, \ldots, \partial f/\partial x_n)$ | §1.4, Thm. 1.4.2 |
| $\boldsymbol{\epsilon}$ | Multivariate perturbation vector $(\epsilon_1, \ldots, \epsilon_n)$ | §1.4, Thm. 1.4.2 |
| $\kappa_{\mathrm{rel}}(f, v)$ | Relative condition number $= |v \cdot f'(v) / f(v)|$ | §1.9, Thm. 1.9.1 |
| $n$ (in §1.4) | Integer quotient $\lfloor x/y \rfloor$ in modular reduction (not mean motion) | §1.4, Cor. 1.4.9 |
| $p$ (in §1.6) | Order of convergence of an iterative method | §1.6, Thm. 1.6.3 |
| $C$ (in §1.6) | Convergence constant in $p$-th order methods | §1.6, Thm. 1.6.3 |
| $m$, $M$ (in §1.6) | Bounds on $|f'|$ and $|f''|$ in Newton's method | §1.6, Thm. 1.6.2 |
| **[MP]** | Matched Pair annotation (value is exact within the coupled TLE+SGP4 system) | §1.2, Ch. 3 |

---

## §1.1 Introduction

Orbit propagation transforms a set of orbital elements, determined at an epoch by fitting observations, into a predicted position and velocity at a future time. This transformation involves:

- Physical constants known only to finite accuracy (the gravitational parameter $\mu$, zonal harmonics $J_n$, atmospheric density parameters).
- A mathematical model that approximates the true equations of motion (perturbation series truncated at finite order, simplified density profiles, neglected forces).
- A computational process that evaluates the model using finite-precision arithmetic, truncated series, and iterative solvers.

Each of these three sources introduces error into the final result, and the three sources are fundamentally independent: improving the computational precision does not compensate for a truncated perturbation model, and extending the model to higher order does not reduce measurement uncertainty in the gravitational parameter.

A meaningful error budget requires that these three sources be tracked separately through every operation. This chapter develops the mathematical framework for doing so: definitions of the three error types (§1.2), propagation rules through all arithmetic and transcendental operations (§§1.3–1.4), series and iterative errors (§§1.5–1.6), cancellation diagnostics (§1.7), error combination (§1.8), and the reliable-digits decision criterion (§1.9).

The scalar results of this chapter are extended to vectors and matrices — including the 13$\times$13 state matrix framework — in Chapter 2. The interaction between error tracking and the SGP4 matched-pair constraint (why "better" constants can degrade accuracy) is developed in Chapter 3. Fast-convergence alternatives to Taylor series (continued fractions, Padé approximants) and the tolerance parameter that governs series truncation are developed in Chapter 4.

---

## §1.2 Fundamental Definitions

**Definition 1.2.1** (Measurement error)**.** *Let $\tilde{v}$ denote the true value of a physical quantity and $v_{\mathrm{meas}}$ the value obtained by measurement. The measurement error bound $\sigma_m \geq 0$ satisfies*

$$
|v_{\mathrm{meas}} - \tilde{v}| \leq \sigma_m. \tag{1.1}
$$

*Measurement error is a property of the physical input. It cannot be reduced by computation; it can only be reduced by better measurement.*

For quantities that are exact by definition (e.g., the WGS-72 equatorial radius $a_E = 6378.135$ km, which defines the coordinate system rather than measuring the Earth), $\sigma_m = 0$ exactly. For quantities determined by observation (e.g., the gravitational parameter $\mu = GM$), $\sigma_m > 0$ reflects the current state of measurement science. [M.1.1]

**Definition 1.2.2** (Precision error)**.** *Let $v$ be a mathematical expression and $\hat{v}$ the result of evaluating $v$ by a computational process (floating-point arithmetic, truncated series, iterative convergence to tolerance $\epsilon$). The precision error bound $\delta_p \geq 0$ satisfies*

$$
|\hat{v} - v| \leq \delta_p. \tag{1.2}
$$

*Precision error is a property of the computational process. It can be reduced by using wider arithmetic types, more series terms, or tighter convergence tolerances.*

Sources of precision error include:

1. **Representation error.** A real number $v$ stored in IEEE 754 binary64 (double precision) satisfies $|\mathrm{fl}(v) - v| \leq \tfrac{1}{2}\,\mathrm{ulp}(v)$, where $\mathrm{ulp}(v) = 2^{e-52}$ for $v$ in the binade $[2^e, 2^{e+1})$. [P.1.1]
2. **Rounding accumulation.** Each floating-point operation introduces a relative error of at most $\epsilon_{\mathrm{mach}} = 2^{-52} \approx 1.1 \times 10^{-16}$ for binary64. Chains of $n$ operations can accumulate $O(n \epsilon_{\mathrm{mach}})$ error in the worst case.
3. **Series truncation to finite terms.** Evaluating $\sum_{k=0}^{N} a_k$ when the exact value is $\sum_{k=0}^{\infty} a_k$ introduces $\delta_p = |\sum_{k=N+1}^{\infty} a_k|$. (See §1.5 for the distinction between precision and accuracy truncation.)
4. **Iterative convergence.** A solver producing $x_k$ with $|x_k - x_{k-1}| < \epsilon$ has $\delta_p$ depending on the convergence rate and the stopping criterion (§1.6).

**Definition 1.2.3** (Accuracy error)**.** *Let $v_{\mathrm{model}}$ denote the result of an exact evaluation of a mathematical model, and $v_{\mathrm{true}}$ the corresponding quantity in the real physical system. The accuracy error bound $\delta_a \geq 0$ satisfies*

$$
|v_{\mathrm{model}} - v_{\mathrm{true}}| \leq \delta_a. \tag{1.3}
$$

*Accuracy error is a property of the mathematical model. It can be reduced only by using a more complete model (more perturbation terms, higher-fidelity density, additional forces), never by improving computational precision.*

Examples:
- SGP4 retains $J_2$, $J_3$, $J_4$ zonal harmonics. Omitting $J_5$ through $J_{2190}$ (the EGM2008 degree) introduces $\delta_a \sim 100$ m for LEO position after one day. [A.1.1]
- The power-law density model $\rho(r) = \rho_0((q_0 - s)/(r - s))^\tau$ is a curve fit to the real atmosphere. The model-data residual is an accuracy error. [A.1.2]
- Brouwer secular rates at order $J_2^2 + J_4$. The omitted $J_2^3$ term is an accuracy error. [A.1.3]

**Remark.** The distinction between precision truncation and accuracy truncation depends on what the series represents. If a series evaluates a known mathematical function (e.g., $\sin x = \sum (-1)^k x^{2k+1}/(2k+1)!$), truncation is a precision error — the exact function is known and more terms reduce the gap. If a series represents a perturbation expansion whose higher-order terms model additional physics (e.g., the third-order Brouwer generating function captures $J_2^3$ coupling), truncation is an accuracy error — the "exact" answer requires a different model.

**Definition 1.2.4** (Tracked value)**.** *A tracked value is an ordered quadruple*

$$
\mathcal{V} = (v,\, \sigma_m,\, \delta_p,\, \delta_a) \tag{1.4}
$$

*where $v$ is the computed value and $\sigma_m, \delta_p, \delta_a \geq 0$ are the three error bounds. The true value $\tilde{v}$ satisfies*

$$
|\tilde{v} - v| \leq \sigma_m + \delta_p + \delta_a. \tag{1.5}
$$

**Definition 1.2.5** (Tier classification)**.** *Every constant and derived quantity is classified into one of four tiers based on which error categories are nonzero:*

| Tier | $\sigma_m$ | $\delta_p$ | $\delta_a$ | Description | Example |
|------|-----------|-----------|-----------|-------------|---------|
| *I* | $0$ | repr. only | $0$ | Exact by definition | $a_E = 6378.135$ km |
| *II* | $> 0$ | repr. only | $0$ | Measured constant | $\mu = 3.986008 \times 10^5$ km$^3$/s$^2$ |
| *III* | inherited | accumulated | $0$ | Computed from Tier I/II | $e^2 = 2f - f^2$ |
| *IV* | inherited | accumulated | $> 0$ | Model-dependent | Brouwer $\dot{\Omega}$ |

*The tier of a computed quantity equals the maximum tier of its inputs, except that any operation introducing model truncation promotes the result to Tier IV regardless of input tiers.*

**Remark.** Within the SGP4 matched-pair system (Chapter 3), certain Tier I constants carry $\sigma_m = 0$ by convention even though the physical quantities they approximate have nonzero measurement uncertainty. For example, $\mu = 398600.8$ km$^3$/s$^2$ **[MP]** is exact within the WGS-72 + SGP4 system but differs from the best current measurement by $\sim 0.4$ km$^3$/s$^2$. The Matched Pair Principle explains when and why it is correct to treat such values as Tier I despite their physical imprecision.

---

## §1.3 Propagation Through Arithmetic Operations

**Remark** (Field generality)**.** The results in this section hold for both real-valued ($\mathbb{F} = \mathbb{R}$) and complex-valued ($\mathbb{F} = \mathbb{C}$) operands. The proofs use only the triangle inequality $|a + b| \leq |a| + |b|$, the reverse triangle inequality $||a| - |b|| \leq |a - b|$ (used in Theorem 1.3.3), and the multiplicative property $|ab| = |a||b|$, all of which hold for the complex modulus. Wherever $|x|$ appears below, it denotes the absolute value (for real $x$) or modulus (for complex $x$). The generalization to $\mathbb{C}$ is needed for the SU(2) state framework developed in Chapter 2.

The results apply independently to each error category. The proofs are stated for a generic error bound $\delta$, representing any of $\sigma_m$, $\delta_p$, or $\delta_a$. Let $x$, $y$ be computed values in $\mathbb{F}$ with error bounds $\delta_x$, $\delta_y$ respectively, so $|\tilde{x} - x| \leq \delta_x$ and $|\tilde{y} - y| \leq \delta_y$.

**Theorem 1.3.1** (Addition and subtraction)**.** *Let $z = x \pm y$. Then*

$$
|\tilde{z} - z| \leq \delta_x + \delta_y. \tag{1.6}
$$

*Proof.* The true value satisfies $\tilde{z} = \tilde{x} \pm \tilde{y}$. Then

$$
|\tilde{z} - z| = |(\tilde{x} \pm \tilde{y}) - (x \pm y)| = |(\tilde{x} - x) \pm (\tilde{y} - y)| \leq |\tilde{x} - x| + |\tilde{y} - y| \leq \delta_x + \delta_y
$$

by the triangle inequality. ∎

**Remark.** The absolute error bound $\delta_x + \delta_y$ is the same for addition and subtraction. Subtraction does not increase absolute error. However, when $|x - y| \ll |x|$, the *relative* error $(\delta_x + \delta_y)/|x - y|$ can be much larger than the relative errors of the inputs. This phenomenon — subtractive cancellation — is analyzed in §1.7.

**Theorem 1.3.2** (Multiplication)**.** *Let $z = x \cdot y$. Then*

$$
|\tilde{z} - z| \leq |x|\,\delta_y + |y|\,\delta_x + \delta_x\,\delta_y. \tag{1.7}
$$

*Proof.* Write $\tilde{x} = x + \epsilon_x$ and $\tilde{y} = y + \epsilon_y$ where $|\epsilon_x| \leq \delta_x$, $|\epsilon_y| \leq \delta_y$. Then

$$
\tilde{z} = \tilde{x}\,\tilde{y} = (x + \epsilon_x)(y + \epsilon_y) = xy + x\epsilon_y + y\epsilon_x + \epsilon_x\epsilon_y.
$$

Therefore

$$
|\tilde{z} - z| = |x\epsilon_y + y\epsilon_x + \epsilon_x\epsilon_y| \leq |x|\,|\epsilon_y| + |y|\,|\epsilon_x| + |\epsilon_x|\,|\epsilon_y| \leq |x|\,\delta_y + |y|\,\delta_x + \delta_x\,\delta_y.
$$

∎

**Remark.** When errors are small relative to values ($\delta_x \ll |x|$, $\delta_y \ll |y|$), the cross term $\delta_x \delta_y$ is negligible and the bound is well approximated by $|x|\,\delta_y + |y|\,\delta_x$. The cross term is retained for rigor, as it matters when errors are not small (e.g., near-zero values with finite absolute error).

**Theorem 1.3.3** (Division)**.** *Let $z = x / y$ with $|y| > \delta_y$ (i.e., the denominator error does not encompass zero). Then*

$$
|\tilde{z} - z| \leq \frac{|x|\,\delta_y + |y|\,\delta_x}{|y|\,(|y| - \delta_y)}. \tag{1.8}
$$

*If $|y| \leq \delta_y$, the bound is undefined (the true denominator may be zero or change sign, and the quotient is unbounded).*

*Proof.* Write $\tilde{x} = x + \epsilon_x$, $\tilde{y} = y + \epsilon_y$ with $|\epsilon_x| \leq \delta_x$, $|\epsilon_y| \leq \delta_y$. Since $|y| > \delta_y$, the denominator $\tilde{y} = y + \epsilon_y$ satisfies $|\tilde{y}| \geq |y| - |\epsilon_y| \geq |y| - \delta_y > 0$, so $\tilde{z} = \tilde{x}/\tilde{y}$ is well defined.

$$
\tilde{z} - z = \frac{\tilde{x}}{\tilde{y}} - \frac{x}{y} = \frac{\tilde{x}\,y - x\,\tilde{y}}{y\,\tilde{y}} = \frac{(x + \epsilon_x)y - x(y + \epsilon_y)}{y\,\tilde{y}} = \frac{y\,\epsilon_x - x\,\epsilon_y}{y\,\tilde{y}}.
$$

Taking absolute values:

$$
|\tilde{z} - z| = \frac{|y\,\epsilon_x - x\,\epsilon_y|}{|y|\,|\tilde{y}|} \leq \frac{|y|\,\delta_x + |x|\,\delta_y}{|y|\,|\tilde{y}|}.
$$

The minimum possible $|\tilde{y}|$ is $|y| - \delta_y$, so

$$
|\tilde{z} - z| \leq \frac{|y|\,\delta_x + |x|\,\delta_y}{|y|\,(|y| - \delta_y)}.
$$

∎

**Remark.** The condition $|y| > \delta_y$ is not a computational guard — it is a mathematical requirement. When the error bound on a denominator encompasses zero, the quotient is genuinely unbounded. In the tracked-value framework, this condition is checked via $\mathrm{rd}(y) > 0$: if the denominator has no reliable digits, the quotient is meaningless. This condition arises concretely in:

- Kepler's equation solvers, which divide by $f'(E) = 1 - e\cos E$ (Ch. 9–10): when $e \to 1$, this denominator can approach zero
- Density ratios $(q_0 - s)/(r - s)$ where $r \to s$ at perigee (Ch. 21)
- The Lyddane modification trigger, where division by $\sin i$ becomes ill-conditioned (Ch. 18)

**Theorem 1.3.4** (Composition of error bounds)**.** *Let $\mathcal{V}_1 = (v_1, \sigma_1, \delta_{p,1}, \delta_{a,1})$ and let $f$ be an operation whose error propagation rule produces $\mathcal{V}_2 = (f(v_1), \sigma_2, \delta_{p,2}, \delta_{a,2})$ with valid bounds $\sigma_2 \geq |\tilde{f}(\tilde{v}_1) - f(v_1)|_{\sigma}$, etc. Let $g$ be a second operation producing $\mathcal{V}_3$ from $\mathcal{V}_2$ by the same propagation rules. Then the bounds in $\mathcal{V}_3$ are valid upper bounds on $|\tilde{g}(\tilde{f}(\tilde{v}_1)) - g(f(v_1))|$ in each error category.*

*Proof.* Each propagation rule (Theorems 1.3.1–1.3.3) computes output error bounds from input error bounds using inequalities of the form $\delta_{\mathrm{out}} \geq |\tilde{z} - z|$ whenever $\delta_{\mathrm{in}} \geq |\tilde{x} - x|$. This is a monotone property: if the input bounds are valid (i.e., they are genuine upper bounds on the true error), then the output bounds are also valid. The composition inherits validity by induction.

Formally: after $f$, $\delta_{p,2}$ satisfies $|\hat{f}(v_1) - f(v_1)| \leq \delta_{p,2}$ (and analogously for $\sigma_2$, $\delta_{a,2}$). The propagation rule for $g$ requires only that its input error bounds are valid upper bounds — it does not require them to be tight. Since $\delta_{p,2}$ is valid, the output $\delta_{p,3}$ from $g$'s propagation rule is also valid. The argument applies to any finite chain of operations. ∎

**Remark.** Theorem 1.3.4 guarantees that the entire SGP4 propagation pipeline — from TLE input through secular update, Kepler solver, short-period corrections, and coordinate transform (Ch. 31–38) — produces valid error bounds at the final state vector, provided each individual step uses the propagation rules of this chapter. The bounds may become increasingly conservative (loose) through long chains, but they remain rigorous upper bounds.

---

## §1.4 Propagation Through Transcendental Functions

**Theorem 1.4.1** (Derivative bound for scalar functions)**.** *Let $f: \mathbb{F} \to \mathbb{F}$ (where $\mathbb{F} = \mathbb{R}$ or $\mathbb{C}$) be continuously differentiable on a region containing the closed disk $\{z : |z - x| \leq \delta\}$. Then*

$$
|f(x + \epsilon) - f(x)| \leq \sup_{|\zeta - x| \leq \delta} |f'(\zeta)| \cdot |\epsilon| \tag{1.9}
$$

*for all $|\epsilon| \leq \delta$.*

*Proof for $\mathbb{F} = \mathbb{R}$:* By the Mean Value Theorem, there exists $\xi$ between $x$ and $x + \epsilon$ such that $f(x + \epsilon) - f(x) = f'(\xi)\,\epsilon$. Since $|\epsilon| \leq \delta$, the point $\xi$ lies in $[x - \delta, x + \delta]$, so $|f'(\xi)| \leq \sup_{|\zeta-x| \leq \delta} |f'(\zeta)|$. The result follows.

*Proof for $\mathbb{F} = \mathbb{C}$:* Let $\gamma$ be the line segment from $x$ to $x + \epsilon$. Since $f$ is holomorphic on the disk,

$$f(x + \epsilon) - f(x) = \int_\gamma f'(\zeta)\,d\zeta.$$

Taking moduli: $|f(x+\epsilon) - f(x)| \leq \int_\gamma |f'(\zeta)|\,|d\zeta| \leq \sup_{|\zeta-x|\leq\delta} |f'(\zeta)| \cdot |\epsilon|$, since the path length is $|\epsilon|$ and every point on $\gamma$ lies within the disk.

Both proofs produce the same bound (1.9). ∎

**Remark.** Theorem 1.4.1 applies to each error category independently. If $\mathcal{X} = (x, \sigma_m, \delta_p, \delta_a)$, then $\mathcal{F} = (f(x), \sigma_m^f, \delta_p^f, \delta_a^f)$ where each bound is computed from (1.9) with the corresponding $\delta$. For the real case, the "disk" $\{z : |z-x| \leq \delta\}$ reduces to the interval $[x-\delta, x+\delta]$.

**Theorem 1.4.2** (Multivariate error propagation)**.** *Let $f: \mathbb{F}^n \to \mathbb{F}$ (where $\mathbb{F} = \mathbb{R}$ or $\mathbb{C}$) be continuously differentiable on an open set containing $\{\mathbf{x} + \boldsymbol{\epsilon} : |\epsilon_i| \leq \delta_i\}$. Then*

$$
|f(\mathbf{x} + \boldsymbol{\epsilon}) - f(\mathbf{x})| \leq \sum_{i=1}^{n} \sup_{\boldsymbol{\xi}} \left|\frac{\partial f}{\partial x_i}(\boldsymbol{\xi})\right| \cdot \delta_i \tag{1.10}
$$

*where the supremum is taken over the box $\{\boldsymbol{\xi} : |\xi_i - x_i| \leq \delta_i\}$, and $|\epsilon_i| \leq \delta_i$ for each $i$.*

*Proof.* Write the change as a telescoping sum along coordinate directions:

$$
f(\mathbf{x} + \boldsymbol{\epsilon}) - f(\mathbf{x}) = \sum_{i=1}^{n} \bigl[f(\mathbf{x} + \boldsymbol{\epsilon}^{(i)}) - f(\mathbf{x} + \boldsymbol{\epsilon}^{(i-1)})\bigr]
$$

where $\boldsymbol{\epsilon}^{(i)}$ agrees with $\boldsymbol{\epsilon}$ in the first $i$ components and is zero in the remaining $n-i$ components (so $\boldsymbol{\epsilon}^{(0)} = \mathbf{0}$ and $\boldsymbol{\epsilon}^{(n)} = \boldsymbol{\epsilon}$). Each bracketed term changes only the $i$-th argument. By Theorem 1.4.1 (which holds for both $\mathbb{R}$ and $\mathbb{C}$) applied to the $i$-th component:

$$
|f(\mathbf{x} + \boldsymbol{\epsilon}^{(i)}) - f(\mathbf{x} + \boldsymbol{\epsilon}^{(i-1)})| \leq \sup_{\boldsymbol{\xi}}\left|\frac{\partial f}{\partial x_i}(\boldsymbol{\xi})\right| \cdot |\epsilon_i|
$$

where the supremum is over the box (for $\mathbb{R}$, this follows from the MVT with equality at some intermediate point; for $\mathbb{C}$, from the line-integral bound of Theorem 1.4.1). Summing:

$$
|f(\mathbf{x} + \boldsymbol{\epsilon}) - f(\mathbf{x})| \leq \sum_{i=1}^{n} \sup_{\boldsymbol{\xi}} \left|\frac{\partial f}{\partial x_i}(\boldsymbol{\xi})\right| \cdot |\epsilon_i| \leq \sum_{i=1}^{n} \sup_{\boldsymbol{\xi}} \left|\frac{\partial f}{\partial x_i}(\boldsymbol{\xi})\right| \cdot \delta_i.
$$

∎

**Remark.** Theorem 1.4.2 is the general foundation for error propagation through multi-input functions. The arithmetic operations of §1.3 are special cases (e.g., $f(x,y) = xy$ has $\partial f/\partial x = y$, $\partial f/\partial y = x$, recovering Theorem 1.3.2 without the cross term). It is used extensively in:

- Ch. 8: orbital elements $(a,e,i,\Omega,\omega,\nu) \to (x,y,z,\dot{x},\dot{y},\dot{z})$ — a 6-input, 6-output function
- Ch. 18–20: short-period and long-period corrections applied to multiple orbital elements simultaneously
- Ch. 30: coordinate transforms involving multiple rotation angles
- Ch. 38: the complete error budget from TLE inputs to final state vector

For vector-valued functions $\mathbf{f}: \mathbb{F}^n \to \mathbb{F}^m$, the bound applies component-wise; the Jacobian matrix $\partial f_j/\partial x_i$ governs the propagation. For $\mathbb{F} = \mathbb{C}$, the partial derivatives are complex (Wirtinger) derivatives when $f$ is holomorphic. Extension to matrix operations over $\mathbb{C}$, including the SU(2) state framework, is developed in Chapter 2.

For tighter bounds on specific functions, the Taylor remainder provides quadratic corrections beyond the linear term.

**Corollary 1.4.1** (Sine bound)**.** *For $|\epsilon| \leq \delta$,*

$$
|\sin(x + \epsilon) - \sin(x)| \leq |\cos x|\,\delta + \frac{\delta^2}{2}. \tag{1.11}
$$

*The bound is capped at $2$ for real arguments (the diameter of the range of $\sin$ on $\mathbb{R}$).*

*Proof.* For real arguments, by Taylor's theorem with Lagrange remainder:

$$
\sin(x + \epsilon) = \sin(x) + \cos(x)\,\epsilon - \frac{\sin(\xi)}{2}\,\epsilon^2
$$

for some $\xi$ between $x$ and $x + \epsilon$. Therefore

$$
|\sin(x + \epsilon) - \sin(x)| \leq |\cos(x)|\,|\epsilon| + \frac{|\sin(\xi)|}{2}\,\epsilon^2 \leq |\cos x|\,\delta + \frac{\delta^2}{2}
$$

since $|\sin(\xi)| \leq 1$ for real $\xi$, and $|\epsilon| \leq \delta$. The bound cannot exceed $|\sin(x + \epsilon) - \sin(x)| \leq |\sin(x+\epsilon)| + |\sin(x)| \leq 2$ for real arguments.

**Domain note ($\mathbb{C}$):** For complex arguments, $\sin$ is entire and the first-order bound $|\cos(x)| \cdot \delta$ from Theorem 1.4.1 applies via the Cauchy estimate. However, the quadratic correction $\delta^2/2$ relies on $|\sin(\xi)| \leq 1$, which fails for complex $\xi$ (e.g., $|\sin(iy)| = |\sinh(y)| \to \infty$), and the cap at 2 does not hold. For complex arguments, use Theorem 1.4.1 directly with $\sup|f'| = \sup|\cos|$ over the disk. In the SGP4 application, all trig arguments are real orbital angles, so the real bounds apply throughout. ∎

**Corollary 1.4.2** (Cosine bound)**.** *For $|\epsilon| \leq \delta$,*

$$
|\cos(x + \epsilon) - \cos(x)| \leq |\sin x|\,\delta + \frac{\delta^2}{2}. \tag{1.12}
$$

*The bound is capped at $2$ for real arguments.*

*Proof.* For real arguments, identical to Corollary 1.4.1, using the Taylor expansion of $\cos$ about $x$:

$$
\cos(x + \epsilon) = \cos(x) - \sin(x)\,\epsilon - \frac{\cos(\xi)}{2}\,\epsilon^2.
$$

Then $|\cos(x+\epsilon) - \cos(x)| \leq |\sin x|\,|\epsilon| + \frac{|\cos \xi|}{2}\,\epsilon^2 \leq |\sin x|\,\delta + \delta^2/2$, capped at $2$ for real arguments.

**Domain note ($\mathbb{C}$):** The same considerations as Corollary 1.4.1 apply: the quadratic term and cap require $|\cos(\xi)| \leq 1$, which holds only for real $\xi$. For complex arguments, use Theorem 1.4.1 directly. ∎

**Corollary 1.4.3** (Two-argument arctangent bound)**.** *Let $z = \mathrm{atan2}(y, x)$, and let $r^2 = x^2 + y^2 > 0$. For errors $\delta_y$, $\delta_x$ on $y$ and $x$ respectively,*

$$
|\mathrm{atan2}(y + \epsilon_y, x + \epsilon_x) - \mathrm{atan2}(y, x)| \leq \frac{|x|}{r^2}\,\delta_y + \frac{|y|}{r^2}\,\delta_x \tag{1.13}
$$

*provided $r^2 > (\delta_x + \delta_y)^2$. The bound is capped at $\pi$. If $r^2 \leq (\delta_x + \delta_y)^2$, the angle is indeterminate and the bound is $\pi$.*

*Proof.* The function $\mathrm{atan2}(y,x) = \arg(x + iy)$ has partial derivatives

$$
\frac{\partial}{\partial y}\mathrm{atan2}(y,x) = \frac{x}{x^2 + y^2}, \qquad \frac{\partial}{\partial x}\mathrm{atan2}(y,x) = \frac{-y}{x^2 + y^2}
$$

wherever $r^2 = x^2 + y^2 > 0$. By the multivariate Mean Value Theorem (applied component-wise),

$$
|\Delta z| \leq \left|\frac{x}{r^2}\right|\,|\epsilon_y| + \left|\frac{-y}{r^2}\right|\,|\epsilon_x| = \frac{|x|}{r^2}\,|\epsilon_y| + \frac{|y|}{r^2}\,|\epsilon_x| \leq \frac{|x|}{r^2}\,\delta_y + \frac{|y|}{r^2}\,\delta_x.
$$

The bound cannot exceed $\pi$ (the maximum change in $\mathrm{atan2}$ from any perturbation that does not cross the branch cut and back). When $r^2 \leq (\delta_x + \delta_y)^2$, the perturbed point $(x + \epsilon_x, y + \epsilon_y)$ may be in any quadrant, so the angle is entirely undetermined and the bound is $\pi$.

**Domain note:** The function $\mathrm{atan2}(y, x)$ is defined only for real arguments. It has no standard complex extension. This corollary is not used in the SU(2) framework (Ch 2), where angles are recovered from quaternion components via algebraic operations, not inverse trigonometric functions. ∎

**Corollary 1.4.4** (Square root bound)**.** *Let $x > 0$ and $\delta < x$. Then*

$$
|\sqrt{x + \epsilon} - \sqrt{x}| \leq \frac{\delta}{2\sqrt{x - \delta}} \tag{1.14}
$$

*for all $|\epsilon| \leq \delta$.*

*Proof.* The derivative of $\sqrt{\cdot}$ is $1/(2\sqrt{\cdot})$, which is positive and decreasing on $(0, \infty)$. On the interval $[x - \delta, x + \delta]$ (which lies in $(0, \infty)$ since $\delta < x$), the supremum of $|(\sqrt{\cdot})'|$ is attained at the left endpoint:

$$
\sup_{\xi \in [x-\delta, x+\delta]} \frac{1}{2\sqrt{\xi}} = \frac{1}{2\sqrt{x - \delta}}.
$$

Applying Theorem 1.4.1 gives the result. ∎

**Remark.** When $\delta \geq x$, the argument $x + \epsilon$ may be zero or negative. In this case, $\sqrt{x + \epsilon}$ may be undefined, and the entire value $\sqrt{x}$ is uncertain. The bound is then set to $\sqrt{x}$ itself (the value could range from $0$ to $\sqrt{x + \delta}$, and $\sqrt{x}$ is an upper bound on the change $|\sqrt{x+\epsilon} - \sqrt{x}|$ for $\epsilon \geq -x$).

**Domain note ($\mathbb{C}$):** For complex arguments, $\sqrt{z}$ is defined with a branch cut (conventionally along the negative real axis). The bound (1.14) holds for $z$ in a neighborhood of the positive real axis (where $\sqrt{\cdot}$ is holomorphic and $(\sqrt{\cdot})' = 1/(2\sqrt{\cdot})$). For $z$ near the branch cut, the function is discontinuous and the bound does not apply. In the SGP4 application, $\sqrt{\cdot}$ is applied to quantities like $1 - e^2$ and $\det M$, which are positive reals in the normal operating regime.

**Corollary 1.4.5** (Exponential bound)**.** *For $|\epsilon| \leq \delta$,*

$$
|e^{x+\epsilon} - e^x| \leq e^x(e^{\delta} - 1). \tag{1.15}
$$

*When $\delta \ll 1$, this simplifies to $|e^{x+\epsilon} - e^x| \lesssim e^x \cdot \delta$.*

*Proof.* The function $e^t$ is convex and increasing. On $[x - \delta, x + \delta]$, the maximum of $|e^{x+\epsilon} - e^x|$ is attained at $\epsilon = \delta$ (since $e^{x+\delta} - e^x > e^x - e^{x-\delta}$ by convexity). Then

$$
|e^{x+\epsilon} - e^x| \leq e^{x+\delta} - e^x = e^x(e^{\delta} - 1).
$$

For small $\delta$, the bound $e^{\delta} - 1 = \delta + \delta^2/2 + \cdots \approx \delta$.

**Domain note ($\mathbb{C}$):** The exponential is entire (holomorphic everywhere on $\mathbb{C}$). For complex arguments, $|e^{z+\epsilon} - e^z| \leq |e^z|(e^{|\epsilon|} - 1)$ still holds, since $|e^w| = e^{\mathrm{Re}(w)} \leq e^{|w|}$. The bound is valid for all $z \in \mathbb{C}$. ∎

**Corollary 1.4.6** (Logarithm bound)**.** *Let $x > 0$ and $\delta < x$. Then for $|\epsilon| \leq \delta$,*

$$
|\ln(x + \epsilon) - \ln x| \leq \ln\!\left(\frac{x}{x - \delta}\right) = -\ln\!\left(1 - \frac{\delta}{x}\right). \tag{1.16}
$$

*When $\delta \ll x$, this simplifies to $|\ln(x+\epsilon) - \ln x| \lesssim \delta/x$.*

*Proof.* The derivative of $\ln$ is $1/t$, which is positive and decreasing on $(0,\infty)$. On $[x - \delta, x + \delta]$, the supremum of $|1/t|$ is $1/(x - \delta)$. By Theorem 1.4.1,

$$
|\ln(x + \epsilon) - \ln x| \leq \frac{\delta}{x - \delta}.
$$

A tighter bound comes from direct evaluation: $\ln(x + \epsilon) - \ln x = \ln(1 + \epsilon/x)$. Since $\ln$ is concave, the maximum of $|\ln(1 + \epsilon/x)|$ for $|\epsilon| \leq \delta$ is attained at $\epsilon = -\delta$ (where $\ln(1 - \delta/x)$ is most negative). Thus

$$
|\ln(x + \epsilon) - \ln x| \leq |\ln(1 - \delta/x)| = -\ln(1 - \delta/x) = \ln\!\left(\frac{x}{x - \delta}\right).
$$

For $\delta/x \ll 1$: $-\ln(1 - \delta/x) = \delta/x + \delta^2/(2x^2) + \cdots \approx \delta/x$.

**Domain note ($\mathbb{C}$):** The complex logarithm has a branch cut (conventionally along the negative real axis). The bound holds for $z$ in a neighborhood of the positive real axis where $\ln$ is holomorphic. For $z$ near the branch cut, the bound does not apply. The same branch-cut consideration applies to Corollary 1.4.7 (power function), which is defined via $\exp$ and $\ln$. ∎

**Corollary 1.4.7** (Power function bound)**.** *Let $x > 0$, $\delta < x$, and $\alpha \in \mathbb{R}$. Then for $|\epsilon| \leq \delta$,*

$$
|x^{\alpha} - (x+\epsilon)^{\alpha}| \leq x^{\alpha}\left[\left(\frac{x}{x - \delta}\right)^{|\alpha|} - 1\right]. \tag{1.17}
$$

*When $\delta \ll x$, this simplifies to $|x^{\alpha} - (x+\epsilon)^{\alpha}| \lesssim |\alpha|\,x^{\alpha-1}\,\delta$.*

*Proof.* Write $x^{\alpha} = e^{\alpha \ln x}$. By composition of Corollaries 1.4.6 and 1.4.5: the logarithm maps input error $\delta$ on $x$ to output error $\leq \ln(x/(x-\delta))$ on $\ln x$; the exponential then maps this to output error $\leq e^{|\alpha| \ln(x/(x-\delta))} - 1 = (x/(x-\delta))^{|\alpha|} - 1$ as a relative factor on $x^{\alpha}$. The linear approximation follows from $|\alpha|\,\delta/x$ being the first-order term. ∎

**Remark.** Corollary 1.4.7 is used in the density model (Ch. 21–22) for $\rho(r) = \rho_0((q_0 - s)/(r-s))^\tau$ with $\tau = 4$, and in mean motion computations $n = \mu^{1/2} a^{-3/2}$ (Ch. 8).

**Domain note ($\mathbb{C}$):** The power function $z^{\alpha} = e^{\alpha \ln z}$ for complex $z$ inherits the branch cut of the logarithm (conventionally along the negative real axis). The bound holds for $z$ in a neighborhood of the positive real axis where $\ln$ is holomorphic. For $z$ near or on the branch cut, $z^{\alpha}$ is discontinuous and the bound does not apply. In the SGP4 application, $x^{\alpha}$ is applied to positive real quantities (semi-major axis, density ratios).

**Corollary 1.4.8** (Modulus bound)**.** *For any $x \in \mathbb{F}$ and $|\epsilon| \leq \delta$,*

$$
||x + \epsilon| - |x|| \leq \delta. \tag{1.18}
$$

*Proof.* By the reverse triangle inequality (valid in both $\mathbb{R}$ and $\mathbb{C}$), $||x+\epsilon| - |x|| \leq |\epsilon| \leq \delta$. Equivalently, the modulus $|\cdot|$ is Lipschitz continuous with Lipschitz constant $1$. ∎

**Corollary 1.4.9** (Modular reduction bound)**.** *Let $z = x \bmod y$ where $y > 0$ and $\delta_y < y$. Define $n = \lfloor x/y \rfloor$ (the integer quotient). Then*

$$
|\widetilde{z} - z| \leq \delta_x + |n|\,\delta_y. \tag{1.19}
$$

*Proof.* Since $z = x - ny$, and $n$ is an integer (carrying no error of its own, though it may change by $\pm 1$ when inputs are perturbed — see Remark below),

$$
|\tilde{z} - z| = |(\tilde{x} - n\tilde{y}) - (x - ny)| = |(\tilde{x} - x) - n(\tilde{y} - y)| \leq \delta_x + |n|\,\delta_y.
$$

∎

**Remark.** The bound assumes $n$ does not change under the perturbation — i.e., $x/y$ is not close to an integer. When $x/y$ is within $(\delta_x + |x/y|\delta_y)/y$ of an integer, the quotient $n$ may jump by $\pm 1$, causing a discontinuity of magnitude $\approx y$ in $z$. In angle reduction (Ch. 7), this corresponds to the value wrapping around $2\pi$, which is geometrically benign (the angle is the same) but requires that the error-tracking implementation handles the branch consistently. See Ch. 7 for the angle-specific treatment.

**Domain note:** Modular reduction requires the floor function $\lfloor x/y \rfloor$, which is defined only for real arguments. This corollary is not applicable to complex operands. In the SU(2) framework (Ch 2), the periodicity of the rotation representation eliminates the need for angle wrapping entirely (see Ch 2, §9.3).

**Corollary 1.4.10** (Arcsine bound)**.** *Let $|x| < 1$ and $\delta < 1 - |x|$. Then for $|\epsilon| \leq \delta$,*

$$
|\arcsin(x + \epsilon) - \arcsin(x)| \leq \frac{\delta}{\sqrt{1 - (|x| + \delta)^2}}. \tag{1.20}
$$

*Proof.* The derivative of $\arcsin$ is $(\arcsin)'(t) = 1/\sqrt{1-t^2}$, which is positive and increasing in $|t|$ on $(-1, 1)$. On $[x - \delta, x + \delta] \subset (-1,1)$, the supremum of $|(\arcsin)'|$ is attained where $|t|$ is largest, i.e., at $|t| = |x| + \delta$:

$$
\sup_{\xi \in [x-\delta, x+\delta]} \frac{1}{\sqrt{1 - \xi^2}} \leq \frac{1}{\sqrt{1 - (|x|+\delta)^2}}.
$$

Applying Theorem 1.4.1 gives the result. ∎

**Remark.** As $|x| \to 1$, the bound diverges: $\arcsin$ has infinite condition number at $x = \pm 1$. This is the mathematical reason that inclination recovery from $\sin i$ is ill-conditioned near $i = 0$ and $i = \pi$, motivating the use of $\mathrm{atan2}$ (Corollary 1.4.3) instead of $\arcsin$ for recovering angles from their sine and cosine. The arccosine bound is analogous, with the same divergence at $|x| \to 1$.

**Domain note:** For real $x \in (-1, 1)$, $\arcsin$ is well-defined and the bound applies directly. For complex arguments, $\arcsin(z) = -i\ln(iz + \sqrt{1-z^2})$ has branch cuts at $|z| \geq 1$ on the real axis. The bound extends to complex $z$ near the real interval $(-1, 1)$ where $\arcsin$ is holomorphic, but diverges as $z$ approaches the branch points. In the SU(2) framework, angles are not recovered via $\arcsin$ — the quaternion components encode orientation directly.

---

## §1.5 Propagation Through Series Truncation

This section addresses the classification and bounding of errors that arise when an infinite series is evaluated to finitely many terms.

**Assumption 1.5.1.** *The series $S = \sum_{k=0}^{\infty} a_k$ converges absolutely.*

**Theorem 1.5.1** (Classification of truncation error)**.** *Let $S = \sum_{k=0}^{\infty} a_k$ and $S_N = \sum_{k=0}^{N} a_k$. The remainder $R_N = S - S_N = \sum_{k=N+1}^{\infty} a_k$ is classified as follows:*

- *If $S$ evaluates a known mathematical function (e.g., $\sin x$, $\ln(1+x)$, a hypergeometric function), then $|R_N|$ is a precision error $\delta_p$, because more terms reduce the gap to the known exact value.* [P.1.2]
- *If $S$ is a perturbation expansion whose higher-order terms model additional physical effects (e.g., the Brouwer generating function at order $J_2^3$), then $|R_N|$ is an accuracy error $\delta_a$, because the "exact" sum requires a different physical model.* [A.1.4]

*Proof.* The distinction follows from Definitions 1.2.2 and 1.2.3. Precision error measures the gap between the computed value and the model's exact prediction; accuracy error measures the gap between the model and physical reality. When the series defines a mathematical function, its exact sum is the model's prediction, and truncation is a computational shortfall. When the series is a perturbation expansion, each order represents a physical approximation, and truncation drops physical content. ∎

**Theorem 1.5.2** (Leibniz alternating series bound)**.** *Let $\sum_{k=0}^{\infty} (-1)^k b_k$ be an alternating series with $b_k > 0$, $b_k \geq b_{k+1}$ for all $k$, and $b_k \to 0$. Then*

$$
|R_N| = \left|\sum_{k=N+1}^{\infty} (-1)^k b_k\right| \leq b_{N+1}. \tag{1.21}
$$

*Proof.* Group consecutive terms of the remainder: $R_N = (-1)^{N+1}(b_{N+1} - b_{N+2}) + (-1)^{N+3}(b_{N+3} - b_{N+4}) + \cdots$. Each parenthesized pair is non-negative by the monotonicity hypothesis. Therefore $|R_N| = R_N \cdot (-1)^{N+1}$ or its negative, and in either case $|R_N| \leq b_{N+1}$ because the partial sums of the remainder oscillate with decreasing amplitude, each overshoot bounded by $b_{N+1}$.

More precisely: $|R_N| = |(-1)^{N+1}b_{N+1} + (-1)^{N+2}b_{N+2} + \cdots|$. The partial sums $\sum_{k=N+1}^{M} (-1)^k b_k$ satisfy $|\sum_{k=N+1}^{M} (-1)^k b_k| \leq b_{N+1}$ for all $M > N$ (by the standard proof of the alternating series test: consecutive partial sums bracket the limit). Taking $M \to \infty$ preserves the inequality. ∎

**Theorem 1.5.3** (Ratio test remainder bound)**.** *Suppose $|a_{k+1}/a_k| \leq r < 1$ for all $k \geq N$. Then*

$$
|R_N| \leq \frac{|a_{N+1}|}{1 - r}. \tag{1.22}
$$

*Proof.* For $k > N$, the hypothesis gives $|a_k| \leq |a_{N+1}| \cdot r^{k-N-1}$. Therefore

$$
|R_N| \leq \sum_{k=N+1}^{\infty} |a_k| \leq |a_{N+1}| \sum_{j=0}^{\infty} r^j = \frac{|a_{N+1}|}{1 - r}.
$$

∎

**Remark.** Theorem 1.5.3 applies to series with geometric-like decay, which includes the $q_0$ series (Ch. 14), the eccentricity functions $G_{lpq}(e)$ (Ch. 15), and the equation of center expansion (Ch. 25). For each application, one must verify the ratio bound $r$ and the index $N$ from which it holds.

**Remark.** The choice of truncation index $N$ is governed by a caller-specified tolerance parameter, developed in Chapter 4. At the SGP4-standard tolerance, $N$ matches the original model's truncation point exactly — reproducing the Spacetrack Report No. 3 [Hoots and Roehrich 1980] result. At tighter tolerances, additional terms are computed, with the remainder bound (Theorem 1.5.2 or 1.5.3) added to $\delta_p$. This is the Matched Pair Principle (Ch. 3) applied to series: the generalized routine accepts a tolerance that, at its standard value, reproduces the original model.

---

## §1.6 Propagation Through Iterative Convergence

**Assumption 1.6.1.** *Let $g: D \to D$ be a contraction mapping on a closed subset $D \subseteq \mathbb{F}$ (where $\mathbb{F} = \mathbb{R}$ or $\mathbb{C}$) with Lipschitz constant $L < 1$, i.e., $|g(x) - g(y)| \leq L|x - y|$ for all $x, y \in D$.*

**Theorem 1.6.1** (Contraction mapping error bound)**.** *Let $x^*$ be the unique fixed point of $g$ (guaranteed by the Banach fixed-point theorem), and let $\{x_k\}$ be the iterative sequence $x_{k+1} = g(x_k)$. Then after $k$ iterations,*

$$
|x_k - x^*| \leq \frac{L}{1 - L}\,|x_k - x_{k-1}|. \tag{1.23}
$$

*The quantity $|x_k - x_{k-1}|$ (the last correction) is computable, and the factor $L/(1-L)$ is determined by the contraction rate.*

*Proof.* By the contraction property,

$$
|x_{k+1} - x_k| = |g(x_k) - g(x_{k-1})| \leq L|x_k - x_{k-1}|.
$$

Inductively, $|x_{k+j} - x_{k+j-1}| \leq L^j |x_k - x_{k-1}|$ for all $j \geq 1$. The fixed point is the limit of the sequence, so

$$
x^* - x_k = \sum_{j=1}^{\infty} (x_{k+j} - x_{k+j-1}).
$$

Taking absolute values and applying the geometric bound:

$$
|x^* - x_k| \leq \sum_{j=1}^{\infty} |x_{k+j} - x_{k+j-1}| \leq |x_k - x_{k-1}| \sum_{j=1}^{\infty} L^j = \frac{L}{1-L}\,|x_k - x_{k-1}|.
$$

∎

**Remark.** For a solver that stops when $|x_k - x_{k-1}| < \epsilon$, the precision error is bounded by $\delta_p \leq L\epsilon/(1-L)$. When $L$ is small (fast convergence), $\delta_p \approx L\epsilon \ll \epsilon$: the last correction is a conservative bound. When $L$ is close to $1$ (slow convergence), the factor $L/(1-L)$ amplifies the correction, and more iterations are needed. [P.1.3]

**Theorem 1.6.2** (Newton's method convergence)**.** *Let $f$ be twice continuously differentiable with $f(x^*) = 0$ and $f'(x^*) \neq 0$. Newton's iteration $x_{k+1} = x_k - f(x_k)/f'(x_k)$ satisfies, for $x_k$ sufficiently close to $x^*$,*

$$
|x_{k+1} - x^*| \leq \frac{M}{2m}\,|x_k - x^*|^2 \tag{1.24}
$$

*where $m = \inf|f'|$ and $M = \sup|f''|$ on a neighborhood of $x^*$.*

*Proof.* By Taylor expansion about $x_k$:

$$
0 = f(x^*) = f(x_k) + f'(x_k)(x^* - x_k) + \frac{f''(\xi)}{2}(x^* - x_k)^2
$$

for some $\xi$ between $x_k$ and $x^*$. Rearranging:

$$
x^* - x_k + \frac{f(x_k)}{f'(x_k)} = -\frac{f''(\xi)}{2f'(x_k)}(x^* - x_k)^2.
$$

The left side is $x^* - x_{k+1}$ (by definition of Newton's step). Taking absolute values:

$$
|x^* - x_{k+1}| = \frac{|f''(\xi)|}{2|f'(x_k)|}|x^* - x_k|^2 \leq \frac{M}{2m}|x^* - x_k|^2.
$$

∎

**Remark.** Quadratic convergence means the number of correct digits approximately doubles at each step. After the last Newton correction $c_k = -f(x_k)/f'(x_k)$, the error satisfies $|x_k - x^*| \leq (M/(2m))\,c_{k-1}^2 \approx c_k$ when convergence is well established. In practice, $|c_k|$ is used as the precision error bound $\delta_p$.

**Domain note:** The proof uses the Lagrange form of the Taylor remainder ("for some $\xi$ between $x_k$ and $x^*$"), which is stated for real-valued functions. For holomorphic $f: \mathbb{C} \to \mathbb{C}$, the same quadratic convergence bound holds with $M = \sup|f''|$ and $m = \inf|f'|$ taken over a disk centered at $x^*$, proved via the integral form of the remainder. The bound (1.24) has the same form in both cases.

**Theorem 1.6.3** (General $p$-th order convergence)**.** *Let $\{x_k\}$ be an iterative sequence converging to $x^*$ with order $p \geq 1$, meaning there exists $C > 0$ such that*

$$
|x_{k+1} - x^*| \leq C\,|x_k - x^*|^p \tag{1.25}
$$

*for all $k$ sufficiently large. Then once $|x_k - x^*| < 1/C^{1/(p-1)}$ (for $p > 1$), the last correction satisfies*

$$
|x_k - x^*| \leq \frac{|x_k - x_{k-1}|}{1 - C\,|x_k - x_{k-1}|^{p-1}}. \tag{1.26}
$$

*When convergence is well established ($C|x_k - x^*|^{p-1} \ll 1$), $|x_k - x^*| \approx |x_k - x_{k-1}|$ and the last correction is a valid precision bound.*

*Proof.* From the reverse triangle inequality, $|x_k - x^*| \leq |x_k - x_{k-1}| + |x_{k-1} - x^*|$. By the convergence bound, $|x_k - x^*| \leq C|x_{k-1} - x^*|^p$, so $|x_{k-1} - x^*| \leq (|x_k - x^*|/C)^{1/p}$. In the convergence regime, $|x_{k-1} - x^*| \gg |x_k - x^*|$, so $|x_k - x_{k-1}| \approx |x_{k-1} - x^*|$ and $|x_k - x^*| \leq C|x_k - x_{k-1}|^p \ll |x_k - x_{k-1}|$. The stated bound follows from the geometric-series argument of Theorem 1.6.1, generalized to $p$-th order contraction. ∎

**Remark.** The cases relevant to this project:

| Method | Order $p$ | Application |
|--------|-----------|-------------|
| Fixed-point (contraction) | 1 | General (Theorem 1.6.1) |
| Newton | 2 | Kepler's equation (Ch. 9), cube root (Ch. 32) |
| Halley | 3 | Kepler's equation (Ch. 9) |
| Householder | 4 | Kepler's equation (Ch. 9) |

For each method, the last correction $|c_k| = |x_k - x_{k-1}|$ gives a valid (conservative) precision error bound $\delta_p$. Higher-order methods converge in fewer iterations but require higher derivatives, which carry their own precision errors. The choice of method is a tradeoff between iteration count and per-iteration cost, governed by the tolerance parameter (Ch. 4).

---

## §1.7 Subtractive Cancellation and Digit Loss

**Theorem 1.7.1** (Condition number of subtraction)**.** *Let $z = x - y$ with $x \neq y$. The condition number of the subtraction, defined as the ratio of relative output error to relative input error, is*

$$
\kappa(x, y) = \frac{|x| + |y|}{|x - y|}. \tag{1.27}
$$

*Proof.* The relative error in $z = x - y$ is

$$
\frac{|\tilde{z} - z|}{|z|} = \frac{|(\tilde{x} - x) - (\tilde{y} - y)|}{|x - y|} \leq \frac{|\tilde{x} - x| + |\tilde{y} - y|}{|x - y|}.
$$

Write $|\tilde{x} - x| \leq |x| \cdot \rho_x$ and $|\tilde{y} - y| \leq |y| \cdot \rho_y$ where $\rho_x, \rho_y$ are the relative errors of $x$ and $y$. Let $\rho = \max(\rho_x, \rho_y)$. Then

$$
\frac{|\tilde{z} - z|}{|z|} \leq \frac{|x|\rho + |y|\rho}{|x - y|} = \frac{|x| + |y|}{|x - y|}\,\rho = \kappa(x,y)\,\rho.
$$

∎

**Theorem 1.7.2** (Digit loss)**.** *Under the hypotheses of Theorem 1.7.1, the number of significant digits lost in the subtraction is approximately*

$$
\text{digits lost} \approx \log_{10} \kappa(x,y) = \log_{10}\!\left(\frac{|x| + |y|}{|x - y|}\right). \tag{1.28}
$$

*Proof.* If the inputs have $d$ reliable digits (relative error $\rho \approx 10^{-d}$), the output has relative error $\kappa \cdot \rho \approx 10^{\log_{10}\kappa - d}$, corresponding to $d - \log_{10}\kappa$ reliable digits. The loss is $\log_{10}\kappa$. ∎

**Example.** Computing $\eta = \sqrt{1 - e^2}$ for $e = 0.999$: the subtraction $1 - e^2 = 1 - 0.998001 = 0.001999$ has $\kappa = (1 + 0.998001)/0.001999 \approx 999$. This loses $\log_{10}(999) \approx 3$ digits. If the inputs had 16 digits (double precision), the result $1 - e^2$ has only 13. The subsequent square root halves the relative error (Corollary 1.4.4 shows $\sqrt{\cdot}$ has condition number $1/2$), recovering half a digit: $\eta$ has approximately 13.5 reliable digits. [P.1.4]

**Remark.** When the condition number is large, the remedy is reformulation, not higher precision. For example, the Lyddane variables $(e\cos\omega, e\sin\omega)$ avoid the $\omega$-singularity at $e = 0$ not by computing $\omega$ more precisely, but by eliminating the ill-conditioned computation entirely. The reliable-digits diagnostic (§1.9) provides the trigger for such switches: when $\mathrm{rd}(z) < d_{\min}$, switch to an alternative formula.

---

## §1.8 Combining Error Bounds

The three error categories $\sigma_m$, $\delta_p$, $\delta_a$ contribute jointly to the total uncertainty in a tracked value. The question of how to combine them depends on whether the categories are treated as independent.

**Definition 1.8.1** (Worst-case combination)**.** *The total error bound for a tracked value $(v, \sigma_m, \delta_p, \delta_a)$ is*

$$
\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a. \tag{1.29}
$$

*This is a rigorous upper bound: the true value $\tilde{v}$ satisfies $|\tilde{v} - v| \leq \delta_{\mathrm{total}}$.*

*Proof.* The gap between the true value $\tilde{v}$ and the computed output $v$ decomposes into a chain of three links through two intermediate values:

1. $\tilde{v}$ — the true physical value (unknowable).
2. $v_{\mathrm{model}}$ — the result of evaluating the mathematical model at the *exact* (error-free) inputs, with infinite precision. This is the model's exact prediction. It differs from $\tilde{v}$ because the model omits physical effects (e.g., higher-order harmonics, neglected forces).
3. $v_{\mathrm{model}}^{\mathrm{meas}}$ — the result of evaluating the same model at the *measured* inputs (which carry measurement error $\sigma_m$), still with infinite precision. It differs from $v_{\mathrm{model}}$ because the inputs are uncertain.
4. $v$ — the computed output: the model evaluated at measured inputs using finite-precision arithmetic. It differs from $v_{\mathrm{model}}^{\mathrm{meas}}$ because of computational rounding, truncation, and convergence tolerance.

By the triangle inequality:

$$
|\tilde{v} - v| = |(\tilde{v} - v_{\mathrm{model}}) + (v_{\mathrm{model}} - v_{\mathrm{model}}^{\mathrm{meas}}) + (v_{\mathrm{model}}^{\mathrm{meas}} - v)|
$$

$$
\leq \underbrace{|\tilde{v} - v_{\mathrm{model}}|}_{\leq\, \delta_a} + \underbrace{|v_{\mathrm{model}} - v_{\mathrm{model}}^{\mathrm{meas}}|}_{\leq\, \sigma_m} + \underbrace{|v_{\mathrm{model}}^{\mathrm{meas}} - v|}_{\leq\, \delta_p} = \delta_a + \sigma_m + \delta_p.
$$

The three inequalities hold by Definitions 1.2.1, 1.2.2, and 1.2.3 respectively. ∎

**Definition 1.8.2** (Root-sum-square combination)**.** *When the three error sources are statistically independent, the RSS combination*

$$
\delta_{\mathrm{RSS}} = \sqrt{\sigma_m^2 + \delta_p^2 + \delta_a^2} \tag{1.30}
$$

*provides an expected-case bound.*

**Proposition 1.8.1** (Independence of the three categories)**.** *The three error categories arise from disjoint physical mechanisms:*

- *$\sigma_m$ arises from limitations of physical measurement apparatus and observation campaigns.*
- *$\delta_p$ arises from the finite wordlength and finite iteration count of the computational process.*
- *$\delta_a$ arises from the choice of mathematical model (which forces and approximations are included).*

*Because these mechanisms operate independently — changing the arithmetic precision does not affect the measurement uncertainty, and improving the measurement does not reduce the model truncation — the three contributions are independent in the probabilistic sense when treated as random variables over an ensemble of computations. The RSS combination (Definition 1.8.2) is therefore a valid expected-case estimate.*

**Remark.** Proposition 1.8.1 is a physical modeling assertion, not a mathematical theorem. It cannot be proved from axioms; it follows from the observation that the three error sources have disjoint remedies and disjoint physical origins. The worst-case bound (Definition 1.8.1) remains the rigorous upper bound used in all formal error propagation, because independence may fail for individual realizations (e.g., a measurement error that happens to align with the model error in a particular direction).

**Remark.** Within a single error category, different inputs may be correlated. For example, TLE orbital elements are fit from a common observation arc, so their measurement errors $\sigma_m$ are correlated. The propagation rules of §1.3 use worst-case (additive) bounds, which remain valid regardless of correlation — they are conservative but never incorrect. The RSS combination should be applied only across the three categories (which are independent by Proposition 1.8.1), not across correlated inputs within a single category, unless independence is separately justified.

---

## §1.9 Reliable Digits

**Definition 1.9.1** (Reliable digits)**.** *For a tracked value $(v, \sigma_m, \delta_p, \delta_a)$ with $v \neq 0$ and $\delta_{\mathrm{total}} > 0$, the number of reliable decimal digits is*

$$
\mathrm{rd}(v) = \left\lfloor -\log_{10}\!\left(\frac{\delta_{\mathrm{total}}}{|v|}\right) \right\rfloor. \tag{1.31}
$$

*When $\delta_{\mathrm{total}} = 0$ (exact value), $\mathrm{rd}(v) = +\infty$. When $\delta_{\mathrm{total}} \geq |v|$ (error dominates), $\mathrm{rd}(v) \leq 0$ (no reliable digits).*

**Theorem 1.9.1** (Reliable digits under univariate functions)**.** *Let $\mathcal{V}_1 = (v_1, \sigma_1, \delta_{p,1}, \delta_{a,1})$ with $v_1 \neq 0$, and let $f$ be a continuously differentiable function applied with the derivative bound propagation rule (Theorem 1.4.1). Define $\mathcal{V}_2 = (f(v_1), \sigma_2, \delta_{p,2}, \delta_{a,2})$. Then, provided $f(v_1) \neq 0$:*

$$
\mathrm{rd}(\mathcal{V}_2) \geq \mathrm{rd}(\mathcal{V}_1) - \left\lceil \log_{10} \kappa_{\mathrm{rel}}(f, v_1) \right\rceil \tag{1.32}
$$

*where $\kappa_{\mathrm{rel}}(f, v_1) = |v_1 \cdot f'(v_1) / f(v_1)|$ is the relative condition number of $f$ at $v_1$.*

*Proof.* By Theorem 1.4.1, the total error satisfies $\delta_{\mathrm{total},2} \leq |f'(v_1)| \cdot \delta_{\mathrm{total},1}$ (using the linear approximation, valid when $\delta_{\mathrm{total},1} \ll |v_1|$). The relative error transforms as:

$$
\frac{\delta_{\mathrm{total},2}}{|f(v_1)|} \leq \frac{|f'(v_1)|}{|f(v_1)|} \cdot \delta_{\mathrm{total},1} = \frac{|f'(v_1)| \cdot |v_1|}{|f(v_1)|} \cdot \frac{\delta_{\mathrm{total},1}}{|v_1|} = \kappa_{\mathrm{rel}} \cdot \frac{\delta_{\mathrm{total},1}}{|v_1|}.
$$

Taking $-\log_{10}$:

$$
-\log_{10}\!\left(\frac{\delta_{\mathrm{total},2}}{|f(v_1)|}\right) \geq -\log_{10}\!\left(\frac{\delta_{\mathrm{total},1}}{|v_1|}\right) - \log_{10} \kappa_{\mathrm{rel}}.
$$

Applying the floor to the left and ceiling to the subtracted term gives the result. ∎

**Remark.** The relative condition number $\kappa_{\mathrm{rel}}$ quantifies how many digits each operation costs (or gains):

| Operation | $\kappa_{\mathrm{rel}}$ | Digit change |
|-----------|------------------------|--------------|
| $x + y$ ($x, y > 0$) | $(|x|+|y|)/|x+y| \leq 1$ | gains $\leq 0$ digits (never loses) |
| $x - y$ ($x \approx y$) | $(|x|+|y|)/|x-y| \gg 1$ | loses $\log_{10}\kappa$ digits |
| $\sin(x)$ near $x = k\pi$ | $|x\cos x/\sin x| \to \infty$ | loses digits |
| $\sin(x)$ near $x = (k+\tfrac{1}{2})\pi$ | $|x\cos x/\sin x| \approx 0$ | gains digits |
| $\sqrt{x}$ | $1/2$ | gains $\approx 0.3$ digits |
| $\mathrm{atan}(x)$ | $|x/(1+x^2)\mathrm{atan}(x)| \leq 1$ | never loses |

This provides the quantitative basis for the formula-switching criterion: when the cumulative digit loss through a chain of operations brings $\mathrm{rd}$ below a minimum threshold, switch to a reformulated computation with smaller $\kappa_{\mathrm{rel}}$.

**Remark.** The reliable-digits count is the decision variable for formula switching. When $\mathrm{rd}(z) < d_{\min}$ for some threshold $d_{\min}$ (typically 2–4), the computation should switch to an alternative formulation that avoids the ill-conditioned step. This is the principled replacement for ad hoc angular thresholds: instead of "if $i < 0.01$ rad, use Lyddane variables," the criterion becomes "if $\mathrm{rd}(\sin i) < 2$, use Lyddane variables." The mathematics determines the switch, not a hardcoded constant.

**Definition 1.9.2** (Reliable digits for near-zero values)**.** *When $v = 0$ (or $|v| < \delta_{\mathrm{total}}$), the relative error is undefined or exceeds $1$, and $\mathrm{rd}(v) \leq 0$. In this regime, the absolute error $\delta_{\mathrm{total}}$ remains the meaningful uncertainty measure. The value is known to lie in the interval $[-\delta_{\mathrm{total}}, +\delta_{\mathrm{total}}]$.*

**Remark.** A near-zero value with small absolute error is not meaningless — it may be the correct physical answer (e.g., eccentricity $e = 0$ for a circular orbit, inclination $i = 0$ for an equatorial orbit). The reliable-digits metric is inappropriate for such values. Downstream chapters must use absolute error bounds, not relative, when the quantity being computed can legitimately be zero. This applies to:

- Eccentricity $e$ near zero (Ch. 10, 18): Lyddane variables avoid dividing by $e$
- Inclination $i$ near zero (Ch. 18): Lyddane variables avoid dividing by $\sin i$
- Short-period corrections that can pass through zero (Ch. 18)

**Definition 1.9.3** (Comparison reliability)**.** *A comparison $x > y$ is reliable if and only if*

$$
|x - y| > \delta_{\mathrm{total},x} + \delta_{\mathrm{total},y}. \tag{1.33}
$$

*When this condition fails, the error intervals of $x$ and $y$ overlap, and the comparison result is indeterminate — it may be true or false depending on the actual errors.*

**Remark.** Comparison reliability matters for branch decisions in the propagator:

- Perigee altitude thresholds (98, 156, 220 km) that select drag modes (Ch. 32)
- The period threshold ($\geq 225$ min) that selects near-space vs. deep-space propagation (Ch. 35)
- The eccentricity threshold ($e < 10^{-6}$) that triggers Lyddane variables (Ch. 18)

When the comparison is unreliable, both branches should be evaluated and the results compared — if they agree to within the error bounds, the branch choice is immaterial. If they diverge, the error bound on the result should reflect the worst case across both branches.

---

## Error Notes

**[M.1.1]** Measurement errors of physical constants. For the WGS-72 gravitational parameter $\mu = GM = 3.986008 \times 10^5$ km$^3$/s$^2$, the measurement uncertainty is approximately $\sigma_m \sim 8 \times 10^{-3}$ km$^3$/s$^2$ (from satellite laser ranging and lunar laser ranging). For defined constants ($a_E$, $1/f$, $\omega$), $\sigma_m = 0$ exactly. *Remedy:* improved measurement campaigns (not a computational concern).

**[P.1.1]** IEEE 754 binary64 representation error. Every real number stored in double precision satisfies $|\mathrm{fl}(v) - v| \leq \tfrac{1}{2}\,\mathrm{ulp}(v) \leq |v| \cdot 2^{-53} \approx |v| \cdot 1.1 \times 10^{-16}$. This gives approximately 15.9 decimal digits. *Remedy:* use extended precision (e.g., `long double`, 80-bit) or arbitrary precision (e.g., `boost::multiprecision::cpp_bin_float_50`).

**[P.1.2]** Series truncation as precision error. When a series evaluates a known mathematical function and is truncated at $N$ terms, the remainder $|R_N|$ is bounded by the methods of §1.5 (Theorems 1.5.2, 1.5.3). The bound is added to $\delta_p$. *Remedy:* compute more terms, or use a continued fraction or Padé approximant with faster convergence (Ch. 4).

**[P.1.3]** Iterative convergence tolerance. A solver stopping at $|x_k - x_{k-1}| < \epsilon$ has $\delta_p \leq L\epsilon/(1-L)$ where $L$ is the contraction constant (Theorem 1.6.1). For Newton's method on Kepler's equation with $e < 1$, $L \leq e < 1$, so $\delta_p \leq e\epsilon/(1-e)$. *Remedy:* tighten $\epsilon$, or use a higher-order method (Halley, Householder) for faster convergence (Ch. 9).

**[P.1.4]** Subtractive cancellation in $1 - e^2$. For near-circular orbits ($e \approx 0$), no cancellation occurs. For highly eccentric orbits ($e \to 1$), cancellation loses $\approx \log_{10}((1+e^2)/(1-e^2))$ digits. At $e = 0.999$: $\approx 3$ digits lost. At $e = 0.99999$: $\approx 5$ digits lost. *Remedy:* rewrite as $(1-e)(1+e)$ to avoid direct subtraction, or use extended precision.

**[A.1.1]** Omitted zonal harmonics. SGP4 retains $J_2$, $J_3$, $J_4$; the EGM2008 model extends to degree 2190. For a 400 km LEO orbit, omitted harmonics contribute $\delta_a \sim 50$–$200$ m in along-track position after one day. *Remedy:* include higher-degree zonals (Ch. 13 Generalization section) or numerical integration with a full gravity model.

**[A.1.2]** Power-law density model. The SGP4 density $\rho(r) = \rho_0((q_0-s)/(r-s))^4$ is a static, spherically symmetric approximation to the real atmosphere, which varies with solar activity, local time, latitude, and altitude. Typical model-data residuals are 20–50% for the density itself, propagating to comparable fractional errors in drag-induced decay rates. *Remedy:* use NRLMSISE-00 or Jacchia-Bowman 2008 density models with space weather inputs (Ch. 21 Generalization section).

**[A.1.3]** Brouwer perturbation theory truncation. The secular rate formulas at order $J_2^2 + J_4$ omit $O(J_2^3)$ terms. For a 400 km LEO orbit with $J_2 \approx 1.08 \times 10^{-3}$, the ratio $J_2^3/J_2^2 \sim 10^{-3}$, so the omitted terms affect secular rates at the $\sim 0.1$% level, contributing $\delta_a \sim 10$–$100$ m after one day. *Remedy:* compute third-order generating function (Ch. 17 Generalization section).

**[A.1.4]** Perturbation expansion as model truncation. Each order of the von Zeipel or Lie transform perturbation expansion adds physical content (coupling between perturbation sources at that order). Truncating at order $n$ omits all coupling at order $n+1$ and above. This is an accuracy error because the dropped terms represent real physical effects, not computational shortcuts. *Remedy:* extend to higher order or use numerical integration.
