# Chapter 4: Approximation Theory and Fast Convergence

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $S_N$ | Partial sum of a series: $S_N = \sum_{k=0}^{N} a_k$ (Ch 1, §1.5) | Ch 1 |
| $R_N$ | Remainder after $N$ terms: $R_N = \sum_{k=N+1}^{\infty} a_k$ (Ch 1, §1.5) | Ch 1 |
| $\tau$ | Tolerance parameter: accuracy threshold in physical units (Ch 3, Def. 3.4.1) | Ch 3 |
| $\tau_{\mathrm{standard}}$ | Standard tolerance (Ch 3, Def. 3.4.2) | Ch 3 |
| $\Gamma(z)$ | Gamma function: $\Gamma(z) = \int_0^\infty t^{z-1} e^{-t}\,dt$ for $\mathrm{Re}(z) > 0$; satisfies $\Gamma(n+1) = n!$ | §4.2 |
| $D^\alpha$ | Fractional derivative operator of order $\alpha \in \mathbb{R}$ | §4.2, Def. 4.2.1 |
| $J^\alpha$ | Riemann–Liouville fractional integral of order $\alpha > 0$ | §4.2, Def. 4.2.2 |
| ${}^C\!D^\alpha$ | Caputo fractional derivative of order $\alpha$ | §4.2, Def. 4.2.4 |
| $E_\alpha(z)$ | Mittag-Leffler function: $E_\alpha(z) = \sum_{k=0}^{\infty} z^k / \Gamma(k\alpha + 1)$ | §4.2 |
| $[p/q]_f(x)$ | Padé approximant of order $[p/q]$ to function $f$ at $x$ | §4.3, Def. 4.3.1 |
| $C_n$ | The $n$-th convergent of a continued fraction: $C_n = A_n / B_n$ | §4.4, Def. 4.4.1 |
| $a_k$, $b_k$ | Partial numerators and partial denominators of a continued fraction | §4.4, Def. 4.4.1 |
| $T_n(x)$ | Chebyshev polynomial of the first kind, degree $n$: $T_n(x) = \cos(n \arccos x)$ | §4.5, Def. 4.5.1 |
| $\binom{\alpha}{k}$ | Generalized binomial coefficient: $\alpha(\alpha-1)\cdots(\alpha-k+1)/k!$ | §4.10, Def. 4.10.1 |
| $(\alpha)_k$ | Pochhammer symbol (rising factorial): $\alpha(\alpha+1)\cdots(\alpha+k-1)$ | §4.10 |
| $(\alpha)_k^{(-)}$ | Falling factorial: $\alpha(\alpha-1)\cdots(\alpha-k+1)$ | §4.10 |
| $\epsilon_{\mathrm{mach}}$ | Machine epsilon of the arithmetic type (Ch 1) | Ch 1 |
| $\delta_p$ | Precision error bound (Ch 1, Def. 1.2.2) | Ch 1 |

---

## §4.1 Introduction

Chapter 1 established remainder bounds for truncated series (Theorems 1.5.2–1.5.3) and convergence bounds for iterative solvers (Theorems 1.6.1–1.6.3). Chapter 3 introduced the tolerance parameter $\tau$ that governs how many terms or iterations a generalized computation uses (Definition 3.4.1) and the matched-pair compatibility criterion that constrains the choice of method at the standard tolerance (Definition 3.7.1).

This chapter develops the approximation methods themselves: the mathematical toolkit for evaluating functions to a specified tolerance with rigorous error bounds. The central question is: given a function $f(x)$, what is the most efficient representation — in terms of operations per digit of accuracy — and what is the tightest computable error bound?

The chapter proceeds from the most fundamental method (Taylor series and its fractional generalization, §4.2) through increasingly powerful alternatives: Padé approximants (§4.3), continued fractions (§4.4), Chebyshev economization (§4.5), convergence acceleration (§4.6), iterative methods of arbitrary order (§4.7), argument reduction (§4.8), polynomial evaluation (§4.9), and the generalized binomial series (§4.10). The synthesis in §4.11 provides a decision framework connecting function properties to method selection and the tolerance parameter.

A cross-cutting principle governs the entire chapter: **exact rational arithmetic for coefficients**. When a series coefficient is a rational number (as in the binomial series, Brouwer polynomials, and geodetic series), it should be represented exactly as an integer ratio, with rounding deferred to the final multiplication by a power of the variable. This eliminates coefficient rounding from the precision error entirely.

This chapter resolves the forward references from Chapter 1: "continued fractions and Padé approximants with faster convergence" (§1.5, [P.1.2]) and "the tolerance parameter" governing method selection (§1.6, Remark after Theorem 1.6.3).

---

## §4.2 Taylor Series: From Integer to Fractional Order

The Taylor series is the baseline approximation method: expand $f(x)$ in powers of $(x - a)$ using derivatives at $a$. This section develops the theory as a single arc — beginning with the classical integer-order expansion, identifying its limitations, and then generalizing to fractional orders $\alpha \in \mathbb{R}$ via the Gamma function and fractional calculus operators.

### §4.2.1 Classical Taylor Series

**Theorem 4.2.1** (Taylor's theorem with Lagrange remainder)**.** *Let $f \in C^{N+1}[a, b]$. Then for any $x \in [a, b]$,*

$$
f(x) = \sum_{k=0}^{N} \frac{f^{(k)}(a)}{k!}(x - a)^k + R_N(x) \tag{4.1}
$$

*where the Lagrange remainder satisfies*

$$
R_N(x) = \frac{f^{(N+1)}(\xi)}{(N+1)!}(x - a)^{N+1} \tag{4.2}
$$

*for some $\xi$ between $a$ and $x$.*

*Proof.* Define the auxiliary function $g(t) = f(x) - \sum_{k=0}^{N} \frac{f^{(k)}(t)}{k!}(x-t)^k - C(x-t)^{N+1}$, where $C$ is chosen so that $g(a) = 0$, i.e., $C = R_N(x)/(x-a)^{N+1}$. Then $g(x) = 0$ (all terms vanish at $t = x$). By Rolle's theorem, $g'(\xi) = 0$ for some $\xi$ between $a$ and $x$. Computing $g'(t)$ and using the telescoping cancellation of the sum's derivative, this yields $-\frac{f^{(N+1)}(\xi)}{N!}(x-\xi)^N + C(N+1)(x-\xi)^N = 0$. Dividing by $(x-\xi)^N$ and solving for $C$ gives $C = f^{(N+1)}(\xi)/((N+1)!)$, establishing (4.2). ∎

**Theorem 4.2.2** (Taylor's theorem with integral remainder)**.** *Under the same hypotheses,*

$$
R_N(x) = \frac{1}{N!}\int_a^x (x - t)^N f^{(N+1)}(t)\,dt. \tag{4.3}
$$

*Proof.* Start from the fundamental theorem of calculus: $f(x) = f(a) + \int_a^x f'(t)\,dt$. Integrate by parts with $u = f'(t)$, $dv = dt$ rewritten as $u = f'(t)$, $dv = d(-(x-t))$:

$$f(x) = f(a) + f'(a)(x-a) + \int_a^x (x-t)\,f''(t)\,dt.$$

Repeating $N$ times, each step extracting one more derivative at $a$ and increasing the power of $(x-t)$ in the integrand, yields (4.3). The integral form is sharper than Lagrange when $f^{(N+1)}$ varies significantly over $[a, x]$: the Lagrange form uses $\max|f^{(N+1)}|$, while the integral form uses the actual values of $f^{(N+1)}(t)$ weighted by $(x-t)^N$.

**Remark.** The integral remainder (4.3) has the form of a convolution integral — it is $J^{N+1}[f^{(N+1)}](x)$ in the notation of Definition 4.2.2 below. This structural identity is the bridge to the fractional generalization: replacing the integer $N$ with a real $\alpha$ in the convolution kernel gives the Riemann–Liouville fractional integral. ∎

**Corollary 4.2.1** (Connection to Ch 1)**.** *The geometric remainder bound of Theorem 1.5.3 is a special case: if $|a_{k+1}/a_k| \leq r < 1$ for all $k \geq N$, then $|R_N| \leq |a_{N+1}|/(1 - r)$. The alternating series bound of Theorem 1.5.2 is likewise a special case for series with terms of strictly decreasing magnitude and alternating sign.*

**Remark.** The classical Taylor series has three fundamental limitations that motivate the rest of this chapter:

1. **Convergence disk.** The series converges only for $|x - a| < \rho$, where $\rho$ is the distance to the nearest singularity of $f$ in the complex plane. For $\arctan(x)$ centered at 0, $\rho = 1$ (poles at $\pm i$), so the series diverges for $|x| > 1$ — yet orbital mechanics requires $\arctan$ for all real arguments.

2. **Slow convergence near the boundary.** As $|x - a| \to \rho$, the geometric convergence rate degrades. For the $q_0$ series (Ch 14) with convergence ratio $e'^2 \approx 0.007$, this is not an issue ($|x|/\rho \ll 1$). But for functions with $|x|/\rho \approx 1$, many terms are needed.

3. **Integer-power basis.** The basis $\{1, (x-a), (x-a)^2, \ldots\}$ cannot represent functions with fractional power behavior at the expansion point. The function $f(x) = \sqrt{x}$ has no Taylor series at $x = 0$ because its derivatives are unbounded there.

Limitations 1–2 motivate Padé approximants (§4.3) and continued fractions (§4.4). Limitation 3 motivates the fractional generalization that follows.

### §4.2.2 Generalizing the Power Rule

The classical derivative $d^m/dx^m\, x^n = n!/(n-m)!\, x^{n-m}$ is restricted to integer orders $m$. The Gamma function $\Gamma(z)$, which satisfies $\Gamma(n+1) = n!$ for $n \in \mathbb{N}$, provides the bridge to arbitrary real order.

**Definition 4.2.1** (Fractional derivative of a monomial)**.** *For $\alpha \in \mathbb{R}$ and $n > -1$,*

$$
D^\alpha x^n = \frac{\Gamma(n+1)}{\Gamma(n - \alpha + 1)}\,x^{n - \alpha}. \tag{4.4}
$$

*Setting $\alpha \in \mathbb{N}$ recovers the integer power rule, since $\Gamma(n+1) = n!$ and $\Gamma(n - m + 1) = (n - m)!$ for integer $m \leq n$.*

**Example 4.2.1.** The half-derivative of $x$:

$$
D^{1/2} x = \frac{\Gamma(2)}{\Gamma(3/2)}\,x^{1/2} = \frac{1}{\sqrt{\pi}/2}\,x^{1/2} = \frac{2}{\sqrt{\pi}}\,\sqrt{x}. \tag{4.5}
$$

The half-derivative of a constant ($n = 0$): $D^{1/2} x^0 = \frac{\Gamma(1)}{\Gamma(1/2)}\,x^{-1/2} = \frac{1}{\sqrt{\pi x}}$.

**Example 4.2.3** (Half-derivatives of monomials)**.** Applying Definition 4.2.1 with $\alpha = 1/2$:

| $f(x)$ | $D^{1/2} f$ | Value |
|---------|-------------|-------|
| $1$ | $\Gamma(1)/\Gamma(1/2) \cdot x^{-1/2}$ | $x^{-1/2}/\sqrt{\pi}$ |
| $x$ | $\Gamma(2)/\Gamma(3/2) \cdot x^{1/2}$ | $2\sqrt{x/\pi}$ |
| $x^2$ | $\Gamma(3)/\Gamma(5/2) \cdot x^{3/2}$ | $8x^{3/2}/(3\sqrt{\pi})$ |
| $x^3$ | $\Gamma(4)/\Gamma(7/2) \cdot x^{5/2}$ | $16x^{5/2}/(5\sqrt{\pi})$ |

The half-derivative of a polynomial is not a polynomial — it introduces half-integer powers. This is the structural reason why the classical Taylor basis $\{1, x, x^2, \ldots\}$ cannot represent $D^{1/2} f$: the half-derivative maps polynomials into a different function space.

### §4.2.3 Fractional Integration and Differentiation

To define a fractional derivative for a general function $f(x)$ (not just monomials), we start with Cauchy's formula for $n$-fold repeated integration:

$$
J^n f(x) = \frac{1}{(n-1)!}\int_a^x (x - t)^{n-1} f(t)\,dt. \tag{4.6}
$$

Replacing $(n-1)!$ with $\Gamma(\alpha)$ gives the fractional integral.

**Definition 4.2.2** (Riemann–Liouville fractional integral)**.** *For $\alpha > 0$,*

$$
J^\alpha f(x) = \frac{1}{\Gamma(\alpha)} \int_a^x (x - t)^{\alpha - 1} f(t)\,dt. \tag{4.7}
$$

**Definition 4.2.3** (Riemann–Liouville fractional derivative)**.** *For $\alpha > 0$ with $n = \lceil \alpha \rceil$,*

$$
D^\alpha f(x) = \frac{d^n}{dx^n} J^{n-\alpha} f(x) = \frac{d^n}{dx^n}\left[\frac{1}{\Gamma(n-\alpha)} \int_a^x (x-t)^{n-\alpha-1} f(t)\, dt\right]. \tag{4.8}
$$

*The operator "integrates up" by the fractional amount $n - \alpha$, then "differentiates down" by the integer amount $n$.*

**Definition 4.2.4** (Caputo fractional derivative)**.** *For $\alpha > 0$ with $n = \lceil \alpha \rceil$,*

$$
{}^C\!D^\alpha f(x) = J^{n-\alpha}\frac{d^n f}{dx^n} = \frac{1}{\Gamma(n-\alpha)} \int_a^x (x-t)^{n-\alpha-1} f^{(n)}(t)\, dt. \tag{4.9}
$$

*The Caputo form reverses the order: differentiate first (integer order $n$), then integrate (fractional order $n - \alpha$).*

**Remark.** The Caputo form is preferred in physical applications because it preserves standard initial conditions: ${}^C\!D^\alpha f$ depends on $f(a), f'(a), \ldots, f^{(n-1)}(a)$, which have direct physical meaning. The Riemann–Liouville form requires fractional initial conditions $\lim_{x \to a^+} D^{\alpha - k} f(x)$, which lack physical interpretation. Both forms reduce to the standard integer derivative when $\alpha \in \mathbb{N}$.

**Remark.** The integral remainder form of Theorem 4.2.2 is structurally identical to the Riemann–Liouville integral (4.7) with $\alpha = N + 1$ and integrand $f^{(N+1)}(t)$. This is not a coincidence — the fractional calculus framework is the natural generalization of the Taylor remainder.

### §4.2.4 The Generalized Taylor Series

**Theorem 4.2.3** (Fractional Mean Value Theorem)**.** *If $f \in C[a, b]$ and $D^\alpha f$ is continuous on $(a, b]$ for $0 < \alpha \leq 1$, then there exists $\xi \in (a, x)$ such that*

$$
f(x) = f(a) + \frac{(x-a)^\alpha}{\Gamma(\alpha + 1)}\,(D^\alpha f)(\xi). \tag{4.10}
$$

*Proof.* Apply the fractional integral $J^\alpha$ to both sides of ${}^C\!D^\alpha f(x) = g(x)$ (where $g$ is continuous by hypothesis). By the composition property $J^\alpha \,{}^C\!D^\alpha f = f - f(a)$ (valid for Caputo derivatives when $f \in C[a,b]$ and $0 < \alpha \leq 1$), we obtain $f(x) - f(a) = J^\alpha g(x) = \frac{1}{\Gamma(\alpha)}\int_a^x (x-t)^{\alpha-1} g(t)\,dt$. By the integral mean value theorem, there exists $\xi \in (a, x)$ such that the integral equals $g(\xi) \cdot \frac{(x-a)^\alpha}{\alpha\,\Gamma(\alpha)} = g(\xi) \cdot \frac{(x-a)^\alpha}{\Gamma(\alpha+1)}$, using $\alpha\,\Gamma(\alpha) = \Gamma(\alpha+1)$. Since $g = {}^C\!D^\alpha f = D^\alpha f$ under these conditions, the result follows. ∎

**Theorem 4.2.4** (Generalized Taylor series — Caputo form)**.** *Under suitable smoothness hypotheses, a function $f$ admits the expansion*

$$
f(x) = \sum_{k=0}^{\infty} \frac{(x-a)^{k\alpha}}{\Gamma(k\alpha + 1)}\,({}^C\!D^{k\alpha} f)(a). \tag{4.11}
$$

*Setting $\alpha = 1$ recovers the classical Taylor series (4.1).*

**Theorem 4.2.5** (Generalized remainder bound)**.** *The remainder after $N$ terms of (4.11) satisfies*

$$
R_N(x) = \frac{(x-a)^{(N+1)\alpha}}{\Gamma((N+1)\alpha + 1)}\,(D^{(N+1)\alpha} f)(\xi) \tag{4.12}
$$

*for some $\xi \in (a, x)$. Setting $\alpha = 1$ recovers the Lagrange remainder (4.2).*

**Remark.** The remainder (4.12) has the same structural role as the integer remainder: it provides a computable upper bound on the truncation error, added to $\delta_p$. The tolerance parameter $\tau$ (Ch 3, Def. 3.4.1) governs truncation identically — add terms until $|R_N(x)| < \tau$.

### §4.2.5 Convergence: Fractional vs Classical

**Theorem 4.2.6** (Convergence rate comparison)**.** *For a smooth analytic function, a fractional series with step $\alpha$ requires $\lceil 1/\alpha \rceil$ times as many terms as the classical series ($\alpha = 1$) to reach the same polynomial degree. However, for functions with non-integer power-law behavior $f(x) \sim C(x - a)^\beta$ with $\beta \notin \mathbb{N}$ near the expansion point, the fractional series with $\alpha | \beta$ converges while the classical series does not.*

*Proof.* For the first claim: the $k$-th term of the fractional series (4.11) has degree $k\alpha$ in $(x-a)$. To reach the same polynomial degree as the $n$-th term of the classical series (degree $n$), the fractional series requires $k = \lceil n/\alpha \rceil$ terms. For $\alpha = 1/2$, this is $2n$ terms — twice as many for the same degree.

For the second claim: let $f(x) = x^\beta$ with $\beta = m\alpha$ for some positive integer $m$. By Definition 4.2.1, $D^{k\alpha} x^\beta = \frac{\Gamma(\beta+1)}{\Gamma(\beta - k\alpha + 1)}\,x^{\beta - k\alpha}$. For $k = m$: $D^{m\alpha} x^{m\alpha} = \frac{\Gamma(m\alpha+1)}{\Gamma(1)} = \Gamma(m\alpha+1)$, a constant. For $k > m$: $\beta - k\alpha = m\alpha - k\alpha < 0$, and $\Gamma(\beta - k\alpha + 1) = \Gamma(m\alpha - k\alpha + 1)$. When $m\alpha - k\alpha + 1 \leq 0$ is a non-positive integer, $\Gamma$ has a pole and $D^{k\alpha} f = 0$. The fractional series terminates at the $m$-th term. The classical Taylor series at $x = 0$ fails because $f^{(n)}(0) = \frac{\Gamma(\beta+1)}{\Gamma(\beta-n+1)}\,0^{\beta-n}$, which is infinite for $n > \lfloor\beta\rfloor$ when $\beta \notin \mathbb{N}$. ∎

**Example 4.2.2** (Singularity resolution)**.** The function $f(x) = \sqrt{x}$ has unbounded derivatives at $x = 0$: $f'(x) = \frac{1}{2}x^{-1/2} \to \infty$. Its classical Taylor series at the origin does not exist. The $\alpha = 1/2$ fractional series terminates in one term:

$$
\sqrt{x} = \frac{x^{1/2}}{\Gamma(3/2)}\,(D^{1/2}\sqrt{x})(0) \cdot \Gamma(3/2) = x^{1/2}. \tag{4.13}
$$

**Remark** (Mittag-Leffler function)**.** The function

$$
E_\alpha(z) = \sum_{k=0}^{\infty} \frac{z^k}{\Gamma(k\alpha + 1)} \tag{4.14}
$$

is the fractional generalization of the exponential: $E_1(z) = e^z$. For $0 < \alpha < 1$, $E_\alpha$ interpolates between the exponential ($\alpha = 1$) and the geometric series $1/(1-z)$ ($\alpha \to 0$). It governs fractional relaxation processes and appears as the natural basis function for fractional-order differential equations.

**Remark** (Fractional Newton's method)**.** The generalized Taylor expansion suggests a fractional Newton iteration for root-finding. For $\alpha = 1/2$, the two-term expansion $f(x_n + h) \approx f(x_n) + \frac{2}{\sqrt{\pi}}(D^{1/2}f)(x_n)\,h^{1/2} + f'(x_n)\,h$ yields a quadratic equation in $u = \sqrt{h}$:

$$
f'(x_n)\,u^2 + \frac{2}{\sqrt{\pi}}\,(D^{1/2}f)(x_n)\,u + f(x_n) = 0 \tag{4.15}
$$

solvable in closed form. This "half-Newton" step uses the half-derivative as an intermediate correction between the function value and the full derivative. For functions with simple zeros ($f'(x^*) \neq 0$), this provides no convergence advantage over standard Newton (§4.7) — the half-derivative information is redundant for smooth functions. The fractional iteration becomes non-trivially superior when the target function has a zero of fractional order: $f(x) \sim C(x - x^*)^\beta$ with $\beta \notin \mathbb{N}$, where standard Newton's convergence degrades. In the current propagation pipeline, all iterative targets have simple zeros (Ch 9, 14, 32), so the classical methods of §4.7 suffice. The fractional framework becomes relevant for enhanced-mode applications such as osculating-to-mean element inversion at high eccentricity (where $\sqrt{1 - e^2}$ introduces half-power behavior in the Jacobian) or closest-approach analysis (where the minimum-distance function has $\sqrt{t - t_{\min}}$ behavior at contact).

---

## §4.3 Padé Approximants

A Taylor polynomial of degree $p + q$ uses $p + q + 1$ coefficients in a polynomial. A Padé approximant $[p/q]$ uses the same $p + q + 1$ coefficients in a rational function $P_p(x)/Q_q(x)$, where $P_p$ has degree $p$ and $Q_q$ has degree $q$. The rational form can represent poles and branch-point behavior that no polynomial can, and often converges in larger domains.

**Definition 4.3.1** (Padé approximant)**.** *Given a formal power series $f(x) = \sum_{k=0}^{\infty} c_k x^k$, the **Padé approximant** $[p/q]_f(x)$ is the rational function $P_p(x)/Q_q(x)$ with $\deg P_p \leq p$, $\deg Q_q \leq q$, $Q_q(0) = 1$, such that*

$$
f(x) - \frac{P_p(x)}{Q_q(x)} = O(x^{p+q+1}). \tag{4.16}
$$

**Theorem 4.3.1** (Existence and uniqueness)**.** *For any formal power series $f$ and non-negative integers $p, q$, the Padé approximant $[p/q]_f$ exists and is unique.*

*Proof.* Write $Q_q(x) = 1 + q_1 x + \cdots + q_q x^q$ and $P_p(x) = p_0 + p_1 x + \cdots + p_p x^p$. The condition (4.16) requires $Q_q(x)\,f(x) - P_p(x) = O(x^{p+q+1})$, i.e., the coefficients of $x^0, x^1, \ldots, x^{p+q}$ in the product $Q_q f - P_p$ must all vanish. Matching coefficients of $x^{p+1}, \ldots, x^{p+q}$ (where $P_p$ contributes nothing, since $\deg P_p \leq p$) gives $q$ linear equations in $q_1, \ldots, q_q$:

$$c_{p+j} + q_1 c_{p+j-1} + \cdots + q_q c_{p+j-q} = 0 \quad (j = 1, \ldots, q)$$

where $c_k$ are the Taylor coefficients of $f$. This is a Toeplitz system that determines $Q_q$ uniquely (the coefficient matrix is always nonsingular for the Padé problem). The numerator coefficients are then determined by matching $x^0, \ldots, x^p$: $p_k = c_k + q_1 c_{k-1} + \cdots + q_{\min(k,q)} c_{k-\min(k,q)}$. ∎

**Theorem 4.3.2** (Padé error order)**.** *The error of the $[p/q]$ Padé approximant satisfies*

$$
|f(x) - [p/q]_f(x)| = O(|x|^{p+q+1}). \tag{4.17}
$$

*This is the same order as the Taylor polynomial of degree $p + q$, but the constant and domain of validity may differ (and are often superior for the Padé).*

**Proposition 4.3.1** (Diagonal Padé error bound)**.** *For $f$ analytic in a disk of radius $\rho$ centered at 0, the diagonal Padé $[n/n]_f$ satisfies*

$$
|f(x) - [n/n]_f(x)| \leq C \left(\frac{|x|}{\rho}\right)^{2n+1}
$$

*for $|x| < \rho$, where $C$ depends on $f$ and $\rho$ but not on $n$. The convergence is geometric with ratio $(|x|/\rho)^2$ — faster than the Taylor series, which converges with ratio $|x|/\rho$.*

*Proof sketch.* The diagonal Padé $[n/n]$ matches $2n + 1$ Taylor coefficients in a rational function of total degree $2n$. By the Montessus de Ballore theorem, for meromorphic $f$ in a disk of radius $\rho$, the $[n/n]$ Padé convergents converge uniformly on compact subsets of the disk (minus the poles). The error decays as the $(2n+1)$-th power of $|x|/\rho$, giving the stated geometric ratio $(|x|/\rho)^2$ per unit increase in $n$. For comparison, the Taylor partial sum $S_{2n}$ has error decaying as $(|x|/\rho)^{2n+1}$ — the same power, but the Padé achieves this with a rational function that also captures pole behavior outside the disk. ∎

**Example 4.3.2.** For $|x|/\rho = 0.9$ (near the convergence boundary), the Taylor convergence ratio is $0.9$ per term, while the diagonal Padé ratio is $(0.9)^2 = 0.81$ per unit increase in $n$. After 5 terms, Taylor achieves $0.9^5 \approx 0.59$ relative error; Padé achieves $0.81^5 \approx 0.35$.

**Example 4.3.1.** The $[1/1]$ Padé approximant to $e^x$ is $(1 + x/2)/(1 - x/2)$. On $[-1, 1]$, the maximum error is $|e^x - [1/1]| \leq 0.056$, compared to $|e^x - S_2(x)| \leq 0.218$ for the degree-2 Taylor polynomial — a factor-of-4 improvement with the same number of coefficients.

**Remark.** The connection to continued fractions (§4.4): the successive diagonal Padé approximants $[0/0], [1/1], [2/2], \ldots$ are the even convergents of a continued fraction representation of $f$. This is developed in §4.4, Proposition 4.4.1.

---

## §4.4 Continued Fractions

For many functions relevant to orbital mechanics, continued fractions converge in regions where the Taylor series diverges or converges slowly. The prototypical example is $\arctan(x)$: its Taylor series converges only for $|x| \leq 1$, while the continued fraction

$$
\arctan(x) = \cfrac{x}{1 + \cfrac{x^2}{3 + \cfrac{4x^2}{5 + \cfrac{9x^2}{7 + \cdots}}}} \tag{4.18}
$$

converges for all real $x$.

**Definition 4.4.1** (Generalized continued fraction)**.** *A **generalized continued fraction** is an expression of the form*

$$
b_0 + \cfrac{a_1}{b_1 + \cfrac{a_2}{b_2 + \cfrac{a_3}{b_3 + \cdots}}} \tag{4.19}
$$

*where $a_k$ are the **partial numerators** and $b_k$ are the **partial denominators**. The **$n$-th convergent** $C_n = A_n / B_n$ is the value obtained by truncating after $n$ partial fractions.*

**Theorem 4.4.1** (Wallis–Euler recurrence)**.** *The numerators $A_n$ and denominators $B_n$ of the convergents satisfy the recurrence*

$$
A_n = b_n A_{n-1} + a_n A_{n-2}, \qquad B_n = b_n B_{n-1} + a_n B_{n-2} \tag{4.20}
$$

*with initial conditions $A_{-1} = 1$, $A_0 = b_0$, $B_{-1} = 0$, $B_0 = 1$.*

*Proof.* By induction on $n$. Base case ($n = 0$): $C_0 = b_0 = A_0/B_0$ with $A_0 = b_0$, $B_0 = 1$. For $n = 1$: $C_1 = b_0 + a_1/b_1 = (b_0 b_1 + a_1)/b_1$, so $A_1 = b_1 b_0 + a_1 = b_1 A_0 + a_1 A_{-1}$ and $B_1 = b_1 = b_1 B_0 + a_1 B_{-1}$ (using $A_{-1} = 1$, $B_{-1} = 0$). Inductive step: the $(n+1)$-th convergent is obtained from the $n$-th by replacing $b_n$ with $b_n + a_{n+1}/b_{n+1}$. In the ratio $C_n = A_n/B_n$, substituting $b_n \to b_n + a_{n+1}/b_{n+1}$ and simplifying gives $A_{n+1} = b_{n+1}A_n + a_{n+1}A_{n-1}$ and $B_{n+1} = b_{n+1}B_n + a_{n+1}B_{n-1}$, completing the induction. ∎

**Theorem 4.4.2** (Śleszyński–Pringsheim convergence criterion)**.** *The continued fraction (4.19) converges if $|b_n| \geq |a_n| + 1$ for all $n \geq 1$.*

*Proof sketch.* By the Stern–Stolz theorem, the condition $|b_n| \geq |a_n| + 1$ ensures that the denominators satisfy $|B_n| \geq 1$ for all $n$ (proved by induction using the recurrence (4.20)), so $|C_n - C_{n-1}| = |a_n/(B_n B_{n-1})| \leq |a_n|$ which, combined with the convergence of $\sum |a_n / B_n B_{n-1}|$, establishes convergence. The Śleszyński–Pringsheim condition is sufficient but not necessary.

**Example 4.4.1.** For the $\arctan$ continued fraction (4.18): $a_n = n^2 x^2$, $b_n = 2n + 1$. The Śleszyński–Pringsheim condition $|b_n| \geq |a_n| + 1$ requires $2n+1 \geq n^2 x^2 + 1$, i.e., $|x| \leq \sqrt{2n/(n^2)} = \sqrt{2/n}$. At $n = 1$: $|x| \leq \sqrt{2}$. The condition fails for large $|x|$ at small $n$, but the continued fraction nevertheless converges for all real $x$ — this follows from the stronger Seidel–Stern theorem for continued fractions with positive partial numerators, which applies since $a_n = n^2 x^2 > 0$. ∎

**Theorem 4.4.3** (Truncation error bound)**.** *For a convergent continued fraction whose convergents $C_n$ alternate (i.e., $C_1 < C_3 < \cdots < f < \cdots < C_2 < C_0$ or vice versa), the truncation error satisfies*

$$
|f - C_n| \leq |C_n - C_{n-1}|. \tag{4.21}
$$

*This is computable: successive convergent differences provide a rigorous error bound, directly usable as the stopping criterion $|C_n - C_{n-1}| < \tau$ for the tolerance parameter.*

**Proposition 4.4.1** (Euler's method)**.** *A power series $f(x) = c_0 + c_1 x + c_2 x^2 + \cdots$ with all $c_k \neq 0$ can be converted into a continued fraction with*

$$
b_0 = c_0, \quad a_1 = c_1 x, \quad b_1 = 1, \quad a_{n} = -\frac{c_n}{c_{n-1}} x, \quad b_n = 1 + \frac{c_n}{c_{n-1}} x \tag{4.22}
$$

*for $n \geq 2$. The continued fraction converges at least wherever the power series converges, and often beyond.*

*Proof sketch.* The Euler equivalence transformation converts the partial sums $S_n(x) = \sum_{k=0}^{n} c_k x^k$ into a sequence of continued fraction convergents by the QD (quotient-difference) algorithm. The QD tableau computes ratios $q_k^{(n)} = c_{n+1}/c_n$ (quotient column) and differences $e_k^{(n)}$ (difference column), from which the partial numerators and denominators (4.22) are read off. The resulting continued fraction's even convergents $C_{2n}$ equal the diagonal Padé approximants $[n/n]_f$, providing the connection stated in the Remark after Example 4.3.1. ∎

**Proposition 4.4.2** (Gauss continued fraction)**.** *The ratio of contiguous hypergeometric functions ${}_2F_1(a, b+1; c+1; x) / {}_2F_1(a, b; c; x)$ admits a continued fraction with partial numerators and denominators that are linear in $n$. This provides continued fraction representations for $\arctan$, $\log(1+x)$, and many special functions arising in perturbation theory.*

*Proof sketch.* The Gauss continued fraction expresses the ratio ${}_2F_1(a, b+1; c+1; x) / {}_2F_1(a, b; c; x)$ as a continued fraction with partial numerators $a_n$ and denominators $b_n$ that are linear functions of $n$. The $\arctan$ continued fraction (4.18) is a special case: $\arctan(x) = x \cdot {}_2F_1(1/2, 1; 3/2; -x^2)$, and the Gauss formula applied to contiguous relations of this ${}_2F_1$ yields $a_n = n^2 x^2$ and $b_n = 2n+1$ as in (4.18). Similarly, $\log(1+x) = x \cdot {}_2F_1(1, 1; 2; -x)$ gives a continued fraction convergent for all $x > -1$. ∎

**Remark** (Modified Lentz algorithm)**.** The recurrence (4.20) can produce $|A_n|$ or $|B_n|$ that grow or shrink exponentially, causing overflow or underflow. The modified Lentz algorithm normalizes at each step: define $D_n = 1/B_n$ and $C_n = A_n/A_{n-1}$, then $f_n = f_{n-1} \cdot C_n \cdot D_n$. This avoids the growth problem while computing convergents incrementally. See [P.4.2] for error analysis.

---

## §4.5 Chebyshev Polynomials and Economization

When a function must be evaluated repeatedly over a bounded interval $[a, b]$, Chebyshev polynomials provide the polynomial of given degree that minimizes the maximum error over the interval.

**Definition 4.5.1** (Chebyshev polynomial of the first kind)**.** *For $x \in [-1, 1]$,*

$$
T_n(x) = \cos(n \arccos x). \tag{4.23}
$$

*Equivalently, $T_n$ satisfies the recurrence $T_0 = 1$, $T_1 = x$, $T_{n+1} = 2x\,T_n - T_{n-1}$.*

**Theorem 4.5.1** (Minimax property)**.** *Among all monic polynomials of degree $n$, $T_n(x)/2^{n-1}$ has the smallest supremum norm on $[-1, 1]$:*

$$
\max_{x \in [-1, 1]} \left|\frac{T_n(x)}{2^{n-1}}\right| = \frac{1}{2^{n-1}} \leq \max_{x \in [-1, 1]} |p(x)| \tag{4.24}
$$

*for any monic polynomial $p$ of degree $n$.*

*Proof.* The polynomial $T_n(x)/2^{n-1}$ is monic (its leading coefficient is $2^{n-1}/2^{n-1} = 1$) and achieves the values $\pm 1/2^{n-1}$ at the $n + 1$ Chebyshev nodes $x_k = \cos(k\pi/n)$ for $k = 0, 1, \ldots, n$, with alternating sign. Suppose a monic polynomial $p$ of degree $n$ satisfies $\|p\|_\infty < 1/2^{n-1}$. Then the difference $d(x) = T_n(x)/2^{n-1} - p(x)$ is a polynomial of degree $\leq n - 1$ (the leading terms cancel since both are monic). At the Chebyshev nodes, $d(x_k) = T_n(x_k)/2^{n-1} - p(x_k)$ has sign determined by $T_n(x_k)/2^{n-1} = \pm 1/2^{n-1}$, since $|p(x_k)| < 1/2^{n-1}$. Therefore $d$ alternates in sign at $n + 1$ points, giving at least $n$ zeros by the intermediate value theorem. But $\deg d \leq n - 1$, so $d$ can have at most $n - 1$ zeros — a contradiction. ∎

**Theorem 4.5.2** (Chebyshev coefficient decay)**.** *If $f$ is analytic inside an ellipse $\mathcal{E}_\rho$ with foci at $\pm 1$ and semi-axis sum $\rho > 1$, then the Chebyshev coefficients $c_n$ in $f(x) = \sum_{n=0}^{\infty} c_n T_n(x)$ satisfy*

$$
|c_n| \leq \frac{2M}{\rho^n} \tag{4.25}
$$

*where $M = \max_{z \in \mathcal{E}_\rho} |f(z)|$.*

**Corollary 4.5.1** (Chebyshev truncation error)**.** *For $x \in [-1, 1]$,*

$$
\left|f(x) - \sum_{k=0}^{N} c_k T_k(x)\right| \leq \frac{2M\,\rho^{-N}}{\rho - 1}. \tag{4.26}
$$

**Proposition 4.5.1** (Clenshaw algorithm)**.** *The sum $S(x) = \sum_{k=0}^{N} c_k T_k(x)$ can be evaluated via the backward recurrence*

$$
y_{N+2} = y_{N+1} = 0, \qquad y_k = 2x\,y_{k+1} - y_{k+2} + c_k \quad (k = N, N-1, \ldots, 1) \tag{4.27}
$$

*with $S(x) = c_0 + x\,y_1 - y_2$. This requires $O(N)$ operations and is numerically stable — the rounding error is bounded by $O(N\,\epsilon_{\mathrm{mach}}\,\max|c_k|)$, which is typically much smaller than the Horner bound for the equivalent monomial-form polynomial.*

*Proof.* Define $y_k = \sum_{j=k}^{N} c_j U_{j-k}(x)$, where $U_m$ is the Chebyshev polynomial of the second kind satisfying $U_0 = 1$, $U_1 = 2x$, $U_{m+1} = 2x\,U_m - U_{m-1}$. The recurrence $y_k = 2x\,y_{k+1} - y_{k+2} + c_k$ follows from the $U_m$ recurrence. For the final step, use the identity $T_j(x) = U_j(x) - x\,U_{j-1}(x)$ (which follows from the cosine addition formula) to obtain $S(x) = \sum_{k=0}^{N} c_k T_k(x) = c_0 + x\,y_1 - y_2$. The numerical stability follows from the fact that $|U_m(x)| \leq m + 1$ for $|x| \leq 1$, so the intermediate values $y_k$ remain bounded. ∎

**Remark** (Chebyshev economization)**.** Given a Taylor polynomial of degree $n$, express it in the Chebyshev basis: $p(x) = \sum_{k=0}^{n} d_k T_k(x)$. Since $|T_k(x)| \leq 1$ on $[-1, 1]$, dropping the highest-degree term adds at most $|d_n|$ to the maximum error. Repeatedly dropping terms produces a degree-$m$ polynomial ($m < n$) with a smaller maximum error on $[-1, 1]$ than the degree-$m$ Taylor polynomial. This is relevant to Ch 15, where the Hough cubic polynomial fits to Hansen coefficients could be replaced by Chebyshev-economized polynomials.

---

## §4.6 Convergence Acceleration

When a series converges slowly — each partial sum approaching the limit by a fixed factor per term (linear convergence) — sequence transformations can produce a faster-converging sequence from the same partial sums.

**Definition 4.6.1** (Sequence transformation)**.** *A map $T: \{s_n\} \to \{\hat{s}_n\}$ is a **convergence acceleration** for $\{s_n\} \to s$ if $\hat{s}_n \to s$ and $(\hat{s}_n - s)/(s_n - s) \to 0$ as $n \to \infty$.*

**Theorem 4.6.1** (Aitken's $\Delta^2$ method)**.** *If $s_n \to s$ with $s_n - s \approx C\lambda^n$ for $|\lambda| < 1$ (linear convergence with ratio $\lambda$), then the transformed sequence*

$$
\hat{s}_n = s_n - \frac{(\Delta s_n)^2}{\Delta^2 s_n} \tag{4.28}
$$

*where $\Delta s_n = s_{n+1} - s_n$ and $\Delta^2 s_n = s_{n+2} - 2s_{n+1} + s_n$, satisfies $(\hat{s}_n - s) = O(\lambda^{2n})$ — superlinear convergence.*

*Proof.* Assume $s_n = s + C\lambda^n + O(\lambda^{2n})$ with $|\lambda| < 1$. Then:

$$\Delta s_n = s_{n+1} - s_n = C\lambda^n(\lambda - 1) + O(\lambda^{2n})$$

$$\Delta^2 s_n = s_{n+2} - 2s_{n+1} + s_n = C\lambda^n(\lambda - 1)^2 + O(\lambda^{2n})$$

Therefore:

$$\frac{(\Delta s_n)^2}{\Delta^2 s_n} = \frac{C^2\lambda^{2n}(\lambda-1)^2 + O(\lambda^{3n})}{C\lambda^n(\lambda-1)^2 + O(\lambda^{2n})} = C\lambda^n + O(\lambda^{2n})$$

and $\hat{s}_n = s_n - (\Delta s_n)^2/\Delta^2 s_n = (s + C\lambda^n) - C\lambda^n + O(\lambda^{2n}) = s + O(\lambda^{2n})$. The convergence ratio has been squared from $|\lambda|$ to $|\lambda|^2$, confirming superlinear acceleration. ∎

**Theorem 4.6.2** (Euler transform for alternating series)**.** *For an alternating series $\sum_{k=0}^{\infty} (-1)^k a_k$ with $a_k > 0$ decreasing, the Euler-transformed series*

$$
\sum_{k=0}^{\infty} \frac{(-1)^k}{2^{k+1}} \Delta^k a_0 \tag{4.29}
$$

*converges faster when the original $a_k$ decrease slowly (e.g., $a_k \sim 1/k$). Here $\Delta^k a_0 = \sum_{j=0}^{k} \binom{k}{j} (-1)^j a_j$ is the $k$-th forward difference.*

**Remark.** The Euler transform is most effective when the original $a_k$ decrease polynomially (e.g., $a_k \sim 1/k$) rather than geometrically. For the geodetic series in this project, $a_k \sim e'^{2k}$ with $e'^2 \approx 0.007$, so convergence is already geometric with ratio $0.007$ — each term is $\sim 150\times$ smaller than the previous. The Euler transform provides negligible improvement for such rapidly convergent series. It is most valuable for slowly converging series such as $\sum (-1)^k/(2k+1)$ (Leibniz series for $\pi/4$), where the terms decrease as $1/k$.

**Remark** (Richardson extrapolation)**.** When a numerical approximation satisfies $A(h) = L + c_1 h^{p_1} + c_2 h^{p_2} + \ldots$ for a step-size parameter $h$, combining evaluations at $h$ and $h/r$ eliminates the leading error term: $\hat{A} = (r^{p_1} A(h/r) - A(h))/(r^{p_1} - 1) = L + O(h^{p_2})$. Richardson extrapolation generalizes Aitken to multi-term asymptotic expansions and is the foundation of Romberg integration.

---

## §4.7 Iterative Methods: Convergence Orders and Stability

Iterative root-finding methods — Newton, Halley, Householder — provide convergence orders from 2 (quadratic) through arbitrary $d + 1$. But convergence order alone does not guarantee convergence: the iteration must start within the **basin of attraction** of the root. This section develops both the convergence theory and the stability criteria, then analyzes the convergence basins of the specific iterative equations in the propagation pipeline.

### §4.7.1 General Theory

**Theorem 4.7.1** (Newton's method — recap from Ch 1 §1.6)**.** *The iteration $x_{n+1} = x_n - f(x_n)/f'(x_n)$ converges quadratically to a simple root $x^*$ of $f$: $|x_{n+1} - x^*| \leq C|x_n - x^*|^2$ where $C = M/(2m)$, $M = \sup|f''|$, $m = \inf|f'|$ in a neighborhood of $x^*$.*

**Theorem 4.7.2** (Newton convergence basin)**.** *The Newton map $N(x) = x - f(x)/f'(x)$ satisfies*

$$
N'(x) = \frac{f(x)\,f''(x)}{f'(x)^2}. \tag{4.30}
$$

*At a simple root $x^*$: $N'(x^*) = 0$ (confirming quadratic convergence). The **basin of attraction** is the maximal connected set containing $x^*$ on which $|N'(x)| < 1$, i.e., the set where*

$$
|f(x)\,f''(x)| < f'(x)^2. \tag{4.31}
$$

*Within this basin, Newton converges monotonically. The **fast convergence radius** is the largest $r$ such that $|N'(x)| < 1$ for all $|x - x^*| < r$. For a function with $f'(x^*) = m > 0$ and $|f''| \leq M$, the fast convergence radius satisfies $r \geq m / M$.*

*Proof.* The Newton multiplier (4.30) is $N'(x) = f(x)\,f''(x)/f'(x)^2$. Near the root $x^*$, $f(x) \approx f'(x^*)(x - x^*) = m(x - x^*)$ and $f'(x) \approx m$, so $|N'(x)| \approx |m(x - x^*) \cdot f''(x)|/m^2 \leq M|x - x^*|/m$. The condition $|N'(x)| < 1$ gives $|x - x^*| < m/M$. ∎

**Theorem 4.7.3** (Halley's method)**.** *The iteration*

$$
x_{n+1} = x_n - \frac{2f(x_n)\,f'(x_n)}{2f'(x_n)^2 - f(x_n)\,f''(x_n)} \tag{4.32}
$$

*converges cubically to a simple root: $|x_{n+1} - x^*| \leq C|x_n - x^*|^3$.*

*Proof sketch.* Halley's method can be derived as Newton's method applied to the auxiliary function $g(x) = f(x)/\sqrt{|f'(x)|}$, which has the property that $g$ has a simple zero at $x^*$ and $g'/g$ has a simpler form than $f'/f$. The Newton step $x_{n+1} = x_n - g/g'$ simplifies to (4.32) after substitution. Equivalently, (4.32) arises from the $[1/1]$ Padé approximation to $f$ about $x_n$: setting the rational approximation $(f(x_n) + f'(x_n)h)/(1 + \gamma h) = 0$ with $\gamma$ chosen to match $f''(x_n)$ yields $h = -2ff'/(2f'^2 - ff'')$. The cubic convergence follows from the three matching conditions ($f$, $f'$, $f''$). ∎

**Theorem 4.7.4** (Halley stability criterion)**.** *The Halley denominator $\Delta_H(x) = 2f'(x)^2 - f(x)\,f''(x)$ must be nonzero for the iteration to be well-defined. At a simple root, $\Delta_H(x^*) = 2f'(x^*)^2 > 0$. The iteration can become unstable when $|\Delta_H(x)|$ is small, which occurs when $|f(x)\,f''(x)| \approx 2f'(x)^2$ — i.e., when the function is large, highly curved, and nearly flat simultaneously.*

**Theorem 4.7.5** (Householder's method of order $d$)**.** *The $d$-th order Householder iteration*

$$
x_{n+1} = x_n + d\,\frac{(1/f)^{(d-1)}(x_n)}{(1/f)^{(d)}(x_n)} \tag{4.33}
$$

*converges with order $d + 1$: $|x_{n+1} - x^*| \leq C|x_n - x^*|^{d+1}$. Newton ($d = 1$) and Halley ($d = 2$) are special cases.*

**Corollary 4.7.1** (Iteration count estimate)**.** *Starting from an approximation with $d_0$ correct digits, a method of order $p$ requires*

$$
k = \left\lceil \frac{\log(d_{\mathrm{target}} / d_0)}{\log p} \right\rceil \tag{4.34}
$$

*iterations to achieve $d_{\mathrm{target}}$ correct digits.*

| Method | Order $p$ | Iterations: 3 → 16 digits | Iterations: 3 → 50 digits | Iterations: 3 → 100 digits |
|--------|-----------|--------------------------|--------------------------|---------------------------|
| Newton | 2 | 3 | 5 | 6 |
| Halley | 3 | 2 | 3 | 4 |
| Householder-3 | 4 | 2 | 3 | 3 |

**Proposition 4.7.1** (Newton for inverse functions)**.** *If $f(x)$ can be evaluated efficiently, then $f^{-1}(y)$ can be computed by applying Newton's method to $g(x) = f(x) - y$: $x_{n+1} = x_n - (f(x_n) - y)/f'(x_n)$. Each iteration requires one evaluation of $f$ and $f'$.*

### §4.7.2 Stability of Orbital Iterative Equations

The general theory of §4.7.1 guarantees convergence within the basin of attraction. For each iterative equation in the propagation pipeline, we must verify: (a) the basin contains the starter value, and (b) the stability margin — the distance from the starter to the basin boundary — is adequate.

**Proposition 4.7.2** (Kepler equation stability)**.** *For the Kepler equation $f(E) = E - e\sin E - M$ with $0 \leq e < 1$:*

*(i) Global convergence of Newton's method.* *The derivative $f'(E) = 1 - e\cos E \geq 1 - e > 0$ for all $E$ and all $e < 1$. Since $f'$ never vanishes, $f$ is strictly increasing. Newton's method converges from any starting value — the basin of attraction is all of $\mathbb{R}$.*

*(ii) Fast convergence radius.* *The fast convergence radius (Theorem 4.7.2) satisfies $r \geq f'_{\min}/|f''|_{\max} = (1-e)/e$. This shrinks with increasing eccentricity:*

| $e$ | $r_{\mathrm{fast}}$ (rad) | $r_{\mathrm{fast}}$ (deg) | Practical implication |
|-----|--------------------------|--------------------------|----------------------|
| $0.1$ | $9.0$ | $516°$ | All starters converge fast |
| $0.5$ | $1.0$ | $57°$ | Comfortable margin |
| $0.9$ | $0.11$ | $6.4°$ | Starter must be within $6°$ of root |
| $0.95$ | $0.053$ | $3.0°$ | Starter must be within $3°$ of root |
| $0.99$ | $0.010$ | $0.6°$ | Starter must be within $0.6°$ of root |
| $0.999$ | $0.001$ | $0.06°$ | Starter must be within $0.06°$ of root |

*Outside $r_{\mathrm{fast}}$, Newton still converges (part (i)) but the initial iterations may be linear rather than quadratic, requiring additional iterations.*

*(iii) Halley denominator bound.* *The Halley denominator $\Delta_H = 2(1 - e\cos E)^2 - (E - e\sin E - M)(e\sin E)$ achieves its minimum value $\sim 2(1-e)^2$ near $E = 0$. For $e = 0.99$, $\min |\Delta_H| \approx 2 \times 10^{-4}$. The Halley iteration remains well-defined but the step size can be large when $|\Delta_H|$ is small, potentially overshooting the root. Halley is therefore preferred over Newton for $e > 0.9$ (cubic convergence compensates for the narrower fast-convergence region) but requires a starter within $r_{\mathrm{fast}}$.*

*Proof.* (i) Since $f'(E) = 1 - e\cos E \geq 1 - e > 0$ for all $E$ when $e < 1$, the function $f$ is strictly increasing. The values $f(0) = -M < 0$ and $f(2\pi) = 2\pi - M > 0$ (for $M \in (0, 2\pi)$) guarantee exactly one root $E^* \in (0, 2\pi)$ by the intermediate value theorem. For Newton's method on a strictly monotone function with a unique root, the iterates form a monotone sequence bounded by $E^*$: if $E_0 > E^*$, the sequence is decreasing (since $f(E_0) > 0$ and $f' > 0$ give a negative correction); if $E_0 < E^*$, the convexity structure ensures the iterate jumps past $E^*$ on the first step and then decreases monotonically. In either case, convergence is guaranteed.

(ii) At the root, $f(E^*) = 0$, $f'(E^*) = 1 - e\cos E^*$, $f''(E^*) = e\sin E^*$. The minimum of $f'$ over all $E$ is $f'_{\min} = 1 - e$ (at $E = 0$), and $|f''| \leq e$ everywhere. By Theorem 4.7.2, $r_{\mathrm{fast}} \geq f'_{\min} / |f''|_{\max} = (1-e)/e$.

(iii) The Halley denominator at $E = 0$ with $M$ small: $\Delta_H = 2(1-e)^2 - 0 = 2(1-e)^2$. This is the global minimum because $f'$ is minimized and $|f\cdot f''|$ is minimized (both near zero) at $E \approx 0$. For $e = 0.99$: $\Delta_H \approx 2 \times 10^{-4}$. ∎

**Proposition 4.7.3** (Cube root stability)**.** *For the cube root equation $f(a) = a^3 - b$ with $b > 0$:*

*(i) Basin of attraction.* *The Newton multiplier is $N'(a) = 2/3 - 2b/(3a^3)$. The condition $|N'(a)| < 1$ gives $a > (2b/5)^{1/3} = (2/5)^{1/3}\,a^*$ where $a^* = b^{1/3}$. Since $(2/5)^{1/3} \approx 0.737$, the basin is $[0.737\,a^*,\, \infty)$.*

*(ii) Starter adequacy.* *The starter $a_0 = b$ satisfies $a_0/a^* = b^{2/3} > 1$ for $b > 1$ (which holds for all SGP4 orbital radii in Earth radii). The starter is well within the basin.*

*(iii) Convergence rate.* *At the starter $a_0 = b$: $|N'(b)| = |2/3 - 2/(3b^2)| \approx 2/3$ for $b > 1$. The first iteration reduces the error by a factor of $\sim 2/3$, after which quadratic convergence takes over. Typical convergence: 4–5 iterations for double precision.*

**Proposition 4.7.4** (Modified Kepler equation stability)**.** *For the SGP4 modified Kepler equation $f(x) = U - a_{yn}\cos x + a_{xN}\sin x - x$ (Lyddane form):*

*(i) Global convergence.* *The derivative $f'(x) = a_{yn}\sin x + a_{xN}\cos x - 1$ satisfies $|f'(x) + 1| = |a_{yn}\sin x + a_{xN}\cos x| \leq \sqrt{a_{yn}^2 + a_{xN}^2} = e$ (the eccentricity). Therefore $f'(x) \leq e - 1 < 0$ for $e < 1$ — the function is strictly decreasing. Newton converges from any starting value.*

*(ii) Eccentricity limits.* *As $e \to 0$: $f(x) \to U - x$, $f'(x) \to -1$. The equation becomes trivial (root $x^* = U$ in one step). As $e \to 1$: the fast convergence radius shrinks as $(1-e)/e$, mirroring the standard Kepler equation. The Lyddane parameterization eliminates the $\omega$-singularity at $e = 0$ but does not widen the convergence basin at $e \to 1$.*

**Proposition 4.7.5** (J₂ fixed-point iteration stability)**.** *For the J₂ refinement $e^2 = g(e^2)$ where $g(e^2) = 3J_2 + h(e^2)$ with $|h'| \leq m \approx 0.003$:*

*The iteration is a contraction mapping with Lipschitz constant $L \approx 0.003 \ll 1$ (Ch 1, Theorem 1.6.1). The basin of attraction is the entire domain $e^2 \in [0, 1)$. Convergence is geometric with ratio $L$: after $k$ iterations, $|e^2_k - e^{*2}| \leq L^k |e^2_0 - e^{*2}|$. For double precision: $k = \lceil 16 / \log_{10}(1/L) \rceil = \lceil 16 / 2.5 \rceil = 7$ iterations suffice, though in practice 3–5 iterations reach convergence because the starter $e^2_0 = 3J_2$ is already close to the root.*

### §4.7.3 Fractional-Order Zeros

**Remark.** When $f$ has a zero of fractional order $\beta$ at $x^*$ — i.e., $f(x) \sim C(x - x^*)^\beta$ with $\beta \notin \mathbb{N}$ — Newton's method degrades. For $\beta < 1$, $f'(x^*) = \infty$ (the Newton step vanishes near the root, producing linear convergence). For $\beta > 1$ and $\beta \notin \mathbb{N}$, $f'(x^*) = 0$ (the Newton step diverges near the root). The fractional Newton iteration of §4.2.5, Equation (4.15), is designed for this regime. In the current propagation pipeline, all iterative targets have simple zeros with $f'(x^*) \neq 0$ (Propositions 4.7.2–4.7.5), so the classical methods suffice. The fractional framework becomes relevant for enhanced-mode applications such as osculating-to-mean element inversion at high eccentricity, where $\sqrt{1 - e^2}$ introduces half-power behavior in the Jacobian.

---

## §4.8 Argument Reduction

Before evaluating a transcendental function by series or continued fraction, reducing the argument to a small range can dramatically reduce the number of terms needed.

**Theorem 4.8.1** (Additive argument reduction)**.** *For a function $f$ with period $P$, $f(x) = f(x')$ where $x' = x - nP$ and $n = \lfloor x/P + 1/2 \rfloor$. The reduced argument satisfies $|x'| \leq P/2$. The precision loss in the reduction is bounded by*

$$
\delta_p(x') \leq \delta_p(x) + |n| \cdot \delta_p(P) + \epsilon_{\mathrm{mach}} \cdot |nP| \tag{4.35}
$$

*where the last term arises from rounding in the subtraction $x - nP$ (Ch 1, §1.7).*

**Example 4.8.1** (Precision loss in modular reduction)**.** For $x = 10^8$ rad and $P = 2\pi$ in double precision ($\sim 16$ digits), $n = \lfloor x/(2\pi) + 1/2 \rfloor \approx 1.59 \times 10^7$. The subtraction $x' = x - n \cdot 2\pi$ involves numbers of magnitude $\sim 10^8$, so the result $x'$ (magnitude $\leq \pi \approx 3.14$) loses $\log_{10}(10^8/3.14) \approx 7.5$ digits — leaving only $\sim 8$ reliable digits in $x'$. At $x = 10^{15}$ rad (a time argument in seconds over decades), all 16 digits are lost and $x'$ is pure noise. This is why large time arguments in sidereal time computations (Ch 29) require the Payne-Hanek algorithm or multi-word arithmetic for the period $P$.

**Theorem 4.8.2** (Multiplicative argument reduction)**.** *For the exponential: $\exp(x) = (\exp(x/2^k))^{2^k}$. Choose $k$ so that $|x/2^k| < \epsilon$; then the Taylor series for $\exp(x/2^k)$ converges in $O(1)$ terms, and the result is reconstructed by $k$ squarings. For trigonometric functions, the half-angle formulas provide the analogous reduction:*

$$
\sin x = 2\sin(x/2)\cos(x/2), \qquad \cos x = 2\cos^2(x/2) - 1. \tag{4.36}
$$

**Remark** (Payne-Hanek algorithm)**.** When $|x|$ is large and the period $P$ is irrational (e.g., $P = 2\pi$), the modular reduction $x \bmod P$ can lose all significant digits. The Payne-Hanek algorithm avoids this by using a pre-computed high-precision representation of $1/P$ (stored as a long integer bit-string) to extract the fractional part of $x/P$ directly, without computing the full quotient $n$.

**Remark** (Reconstruction error)**.** Each squaring in the reconstruction $(\exp(x/2^k))^{2^k}$ doubles the relative error: after $k$ squarings, $\delta_p(\exp(x)) \leq 2^k \cdot \delta_p(\exp(x/2^k))$. Since $k = O(\log|x|)$, the total reconstruction error is $O(\log|x| \cdot \epsilon_{\mathrm{mach}})$.

---

## §4.9 Polynomial Evaluation

**Definition 4.9.1** (Horner form)**.** *The polynomial $p(x) = c_0 + c_1 x + \cdots + c_n x^n$ is evaluated as*

$$
p(x) = c_0 + x(c_1 + x(c_2 + \cdots + x \cdot c_n)). \tag{4.37}
$$

*This requires $n$ multiplications and $n$ additions — the same as naive evaluation, but with fundamentally different rounding error.*

**Theorem 4.9.1** (Horner rounding bound)**.** *The computed Horner evaluation $\hat{p}(x)$ satisfies*

$$
|\hat{p}(x) - p(x)| \leq (2n+1)\,\epsilon_{\mathrm{mach}}\,\tilde{p}(|x|) \tag{4.38}
$$

*where $\tilde{p}(|x|) = \sum_{k=0}^{n} |c_k|\,|x|^k$ is the "condition polynomial." When the polynomial is well-conditioned ($\tilde{p}(|x|) \approx |p(x)|$), this reduces to $O(\epsilon_{\mathrm{mach}}) \cdot |p(x)|$.*

*Proof.* Define the intermediate Horner values: $q_n = c_n$, $q_k = \mathrm{fl}(\mathrm{fl}(q_{k+1} \cdot x) + c_k)$ for $k = n-1, \ldots, 0$. At each step, the multiplication introduces relative error $\leq \epsilon_{\mathrm{mach}}$ and the addition introduces relative error $\leq \epsilon_{\mathrm{mach}}$, so $q_k = (q_{k+1} x + c_k)(1 + \theta_k)$ with $|\theta_k| \leq 2\epsilon_{\mathrm{mach}}$. Unrolling the recurrence: $\hat{p}(x) = q_0 = \sum_{k=0}^{n} c_k x^k \prod_{j=0}^{k} (1 + \theta_j)$ where $|\theta_j| \leq 2\epsilon_{\mathrm{mach}}$. Since $|\prod(1+\theta_j) - 1| \leq (2k+1)\epsilon_{\mathrm{mach}}$ for small $\epsilon_{\mathrm{mach}}$ (by the standard product bound), $|\hat{p}(x) - p(x)| \leq (2n+1)\epsilon_{\mathrm{mach}} \sum |c_k| |x|^k = (2n+1)\epsilon_{\mathrm{mach}}\,\tilde{p}(|x|)$. ∎

**Remark.** Contrast with naive evaluation of $c_0 + c_1 x + c_2 x^2 + \cdots$: computing $x^k$ by repeated multiplication accumulates $O(k\,\epsilon_{\mathrm{mach}})$ relative error in each power, and the subsequent multiplication by $c_k$ adds another $\epsilon_{\mathrm{mach}}$. Summing over $k = 0, \ldots, n$ gives a total bound of $O(n^2\,\epsilon_{\mathrm{mach}}) \cdot \tilde{p}(|x|)$ — a factor of $n$ worse than Horner.

**Theorem 4.9.2** (Error propagation with uncertain argument)**.** *When $x$ carries error $\delta_x$, the total error of the Horner evaluation is*

$$
|\hat{p}(x) - p(\tilde{x})| \leq |p'(x)| \cdot \delta_x + (2n+1)\,\epsilon_{\mathrm{mach}}\,\tilde{p}(|x|) \tag{4.39}
$$

*where the first term is the sensitivity to input error (Ch 1, Theorem 1.4.1) and the second is the rounding error.*

**Remark** (Compensated Horner)**.** When $\tilde{p}(|x|) \gg |p(x)|$ — which occurs for polynomials with alternating-sign coefficients evaluated near $|x| \approx 1$ — the condition number $\tilde{p}(|x|)/|p(x)|$ is large and Horner loses significant digits (Ch 1, §1.7). Compensated Horner evaluation uses error-free transformations (Dekker splitting) to compute the rounding error of each step and accumulate a correction term, effectively doubling the working precision without changing the arithmetic type.

---

## §4.10 The Generalized Binomial Series

The binomial series $(1 + x)^\alpha$ extends the finite binomial theorem to non-integer exponents $\alpha \in \mathbb{R}$. The coefficients involve falling factorials, and when $\alpha$ is rational, the coefficients are exactly rational.

**Definition 4.10.1** (Generalized binomial coefficient)**.** *For $\alpha \in \mathbb{R}$ and $k \in \mathbb{N}$,*

$$
\binom{\alpha}{k} = \frac{(\alpha)_k^{(-)}}{k!} = \frac{\alpha(\alpha-1)(\alpha-2)\cdots(\alpha - k + 1)}{k!}. \tag{4.40}
$$

**Theorem 4.10.1** (Convergence of the binomial series)**.** *For $|x| < 1$ and any $\alpha \in \mathbb{R}$,*

$$
(1 + x)^\alpha = \sum_{k=0}^{\infty} \binom{\alpha}{k} x^k \tag{4.41}
$$

*with remainder bound*

$$
|R_N(x)| \leq \left|\binom{\alpha}{N+1}\right| \frac{|x|^{N+1}}{1 - |x|} \tag{4.42}
$$

*for $|x| < 1$. The convergence rate is geometric with ratio $|x|$.*

*Proof.* The Taylor series for $f(x) = (1+x)^\alpha$ at $x = 0$ has coefficients $f^{(k)}(0)/k! = \alpha(\alpha-1)\cdots(\alpha-k+1)/k! = \binom{\alpha}{k}$, establishing (4.41). For the remainder bound: the ratio of successive terms is $|a_{k+1}/a_k| = |(\alpha - k)/(k+1)| \cdot |x|$. For $k > |\alpha|$, this ratio is bounded by $|x| \cdot (k + |\alpha|)/(k+1) < |x| \cdot (1 + |\alpha|/k)$, which approaches $|x|$ for large $k$. Applying the geometric tail bound (Ch 1, Theorem 1.5.3) with $r = |x|$ gives $|R_N| \leq |\binom{\alpha}{N+1}| \cdot |x|^{N+1}/(1 - |x|)$ for $|x| < 1$. ∎

**Theorem 4.10.2** (Exact rational coefficients)**.** *When $\alpha = p/q$ with $p, q \in \mathbb{Z}$, each coefficient $\binom{p/q}{k}$ is a rational number:*

$$
\binom{p/q}{k} = \frac{p(p-q)(p-2q)\cdots(p-(k-1)q)}{q^k \cdot k!}. \tag{4.43}
$$

*The numerator and denominator are products of integers, computable in exact arithmetic. The only rounding in the partial sum $S_N(x)$ arises from the powers $x^k$.*

**Example 4.10.1.** For $\alpha = 1/2$ (square root): $\binom{1/2}{k} = \frac{1 \cdot (-1) \cdot (-3) \cdots (3 - 2k)}{2^k \cdot k!}$. The first few: $\binom{1/2}{0} = 1$, $\binom{1/2}{1} = 1/2$, $\binom{1/2}{2} = -1/8$, $\binom{1/2}{3} = 1/16$. These are exact rationals.

**Example 4.10.2.** For $\alpha = -2$ (relevant to $(1 + e'^2\cos^2\phi)^{-2}$ in the equipotential ellipsoid, Ch 14): $\binom{-2}{k} = (-1)^k(k+1)$. The coefficients are exact integers.

**Remark** (Pochhammer connection)**.** The rising factorial $(\alpha)_k = \alpha(\alpha+1)\cdots(\alpha+k-1)$ is related to the falling factorial by $(\alpha)_k^{(-)} = (-1)^k(-\alpha)_k$. Both appear in hypergeometric series ${}_pF_q$, connecting the binomial series to the Gauss continued fraction (§4.4, Proposition 4.4.2). The binomial series is the special case ${}_1F_0(\alpha; -; -x)$.

**Remark** (Horner nesting for binomial series)**.** The partial sum $S_N(x) = \sum_{k=0}^{N} \binom{\alpha}{k} x^k$ can be evaluated in Horner form (§4.9) with exact rational coefficients. Factor out $x$: $S_N(x) = 1 + x(c_1 + x(c_2 + \cdots))$ where $c_k = \binom{\alpha}{k}$ is computed in exact rational arithmetic. This combines the rounding advantage of Horner (Theorem 4.9.1) with the coefficient exactness of Theorem 4.10.2 — the only rounding occurs in the $n$ multiplications by $x$.

---

## §4.11 Method Selection and the Tolerance Parameter

This section synthesizes the preceding methods into a decision framework governed by the tolerance parameter $\tau$ (Ch 3, Def. 3.4.1) and the matched-pair compatibility criterion (Ch 3, Def. 3.7.1).

**Table 4.11.1** (Method selection by function property):

| Function property | Preferred method | Convergence | Reference |
|-------------------|-----------------|-------------|-----------|
| Analytic, bounded domain | Chebyshev series | Geometric, rate $\rho^{-n}$ | Thm 4.5.2 |
| Analytic, unbounded domain | Continued fraction | Geometric, convergent decay | Thm 4.4.3 |
| Rational behavior / poles | Padé approximant | Rate $(|x|/\rho)^{2n}$ | Prop 4.3.1 |
| Polynomial in standard form | Horner evaluation | Exact (finite sum) | Thm 4.9.1 |
| Non-integer power-law at origin | Fractional Taylor series | Adapted to $\alpha$ | Thm 4.2.6 |
| Slowly converging alternating | Euler transform | Accelerated | Thm 4.6.2 |
| Inverse of a fast function | Newton/Halley iteration | Quadratic/cubic per step | Thm 4.7.1–4.7.2 |
| Smooth, general case | Taylor with remainder | Geometric, rate $|x|/\rho$ | Thm 4.2.1 |

**Table 4.11.2** (Precision-level strategy):

| Precision range | Primary strategy | Supporting techniques |
|-----------------|-----------------|----------------------|
| $\sim 16$ digits (double) | Horner polynomials, short Taylor/Padé | Argument reduction (§4.8) |
| $\sim 50$ digits | Argument reduction + moderate series; Halley for inverses; continued fractions | Exact rational coefficients (Thm 4.10.2) |
| $\sim 100+$ digits | AGM for $\log$/$\exp$; binary splitting for rational series | Payne-Hanek reduction; Machin-type decomposition |

**Proposition 4.11.1** (Matched-pair compatibility)**.** *Every method listed in Tables 4.11.1–4.11.2, when applied with tolerance $\tau = \tau_{\mathrm{standard}}$ (Ch 3, Def. 3.4.2), produces results satisfying the matched-pair compatibility criterion (Ch 3, Def. 3.7.1, Eq. 3.18).*

*Proof.* By Ch 3, Definition 3.7.1, a computation $g(\tau)$ is matched-pair compatible if $|g(\tau_{\mathrm{standard}}) - f_{\mathrm{SGP4}}| \leq \delta_a(\mathcal{M})$. Each method in Tables 4.11.1–4.11.2 produces a result $g(\tau)$ that converges to the exact mathematical value $f_{\mathrm{exact}}$ as $\tau \to 0$. At $\tau = \tau_{\mathrm{standard}}$, the result satisfies $|g(\tau_{\mathrm{standard}}) - f_{\mathrm{exact}}| \leq \tau_{\mathrm{standard}}$ (by definition of the tolerance parameter, Def. 3.4.1). The Spacetrack Report No. 3 computation [Hoots and Roehrich 1980] satisfies $|f_{\mathrm{SGP4}} - f_{\mathrm{exact}}| \leq \delta_p^{\mathrm{SR3}}$ (the SR3 truncation/rounding error). By the triangle inequality, $|g(\tau_{\mathrm{standard}}) - f_{\mathrm{SGP4}}| \leq \tau_{\mathrm{standard}} + \delta_p^{\mathrm{SR3}}$. Since both $\tau_{\mathrm{standard}}$ and $\delta_p^{\mathrm{SR3}}$ are precision-level quantities (far below $\delta_a(\mathcal{M}) \approx 1$ km), the bound holds. The specific numerical verification for each computation in Chapters 5–35 is required by Ch 3, Corollary 3.7.1. ∎

**Corollary 4.11.1** (Tolerance parameter unification)**.** *Regardless of which method is selected, the external interface is identical: a callable that accepts the time-varying quantities and returns a result with a computable error bound. The tolerance $\tau$ was captured at initialization (Ch 3, Def. 3.4.1) and is invisible to the caller. The method selection affects only the internal convergence rate, not the interface.*

---

## Error Notes

**[P.4.1]** Padé approximant evaluation near poles. The evaluation of a $[p/q]$ Padé approximant requires a division by $Q_q(x)$. When $|Q_q(x)|$ is small (near a pole of the approximant), the division amplifies rounding error by the condition number $|P_p(x)|/|Q_q(x)|$. *Remedy:* detect near-pole conditions and switch to Taylor or continued fraction evaluation in those regions.

**[P.4.2]** Continued fraction overflow/underflow. The Wallis–Euler recurrence (Theorem 4.4.1) can produce $|A_n|$ or $|B_n|$ that grow or shrink exponentially. *Remedy:* normalize $A_n, B_n$ by dividing both by $B_n$ at each step (modified Lentz algorithm), so that the effective denominator is always 1.

**[P.4.3]** Chebyshev coefficient computation. Computing Chebyshev coefficients from function values at Chebyshev nodes uses the discrete cosine transform (DCT) with rounding error $O(\sqrt{n}\,\epsilon_{\mathrm{mach}})$ for $n$ nodes. *Remedy:* for the small degrees typical in orbital mechanics ($n \leq 10$), this error is negligible compared to the geometric coefficient decay.

**[P.4.4]** Horner evaluation with cancellation. When a polynomial has alternating-sign coefficients and $|x| \approx 1$, the condition number $\tilde{p}(|x|)/|p(x)|$ can be large (Ch 1, §1.7). The effective precision loss is $\log_{10}(\tilde{p}/|p|)$ digits. *Remedy:* reformulate the polynomial to avoid cancellation, or use compensated Horner evaluation.

**[P.4.5]** Argument reduction precision loss. Additive reduction $x' = x - nP$ loses $\log_{10}(n)$ digits when $n$ is large (Ch 1, §1.7 — subtractive cancellation). For $x \gg P$ with limited-precision $P$, all digits can be lost. *Remedy:* use the Payne-Hanek algorithm with a high-precision representation of $1/P$.

**[P.4.6]** Binomial series near $|x| = 1$. The binomial series $(1+x)^\alpha$ converges for $|x| < 1$ but the convergence rate degrades as $|x| \to 1$. For $|x|$ near 1, many terms are needed and rounding error accumulates. *Remedy:* use the identity $(1+x)^\alpha = \exp(\alpha \ln(1+x))$ with the methods of §4.8, or apply convergence acceleration (§4.6).

**[P.4.7]** Fractional derivative evaluation. The Riemann–Liouville fractional integral (Definition 4.2.2) requires numerical quadrature when $D^\alpha f$ cannot be computed in closed form. The quadrature error is a precision error added to $\delta_p$, with magnitude depending on the quadrature rule and the smoothness of the integrand. *Remedy:* use Gauss-Jacobi quadrature adapted to the $(x-t)^{\alpha-1}$ weight function, which is exact for polynomials up to degree $2n-1$ with $n$ nodes.
