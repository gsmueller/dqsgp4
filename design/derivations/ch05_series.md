# Chapter 5: Series Evaluation and Error Control

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $S_N$ | Partial sum of a series: $S_N = \sum_{k=0}^{N} a_k$ | Ch 1, §1.5 |
| $R_N$ | Remainder after $N$ terms: $R_N = \sum_{k=N+1}^{\infty} a_k$ | Ch 1, §1.5 |
| $\tau$ | Tolerance parameter (Ch 3, Def. 3.4.1) | Ch 3 |
| $\tau_{\mathrm{standard}}$ | Standard tolerance reproducing SGP4 truncation (Ch 3, Def. 3.4.2) | Ch 3 |
| $\delta_p$ | Precision error bound (Ch 1, Def. 1.2.2) | Ch 1 |
| $\delta_a$ | Accuracy/model error bound (Ch 1, Def. 1.2.3) | Ch 1 |
| $\mathcal{V}$ | A tracked value: the quadruple $(v, \sigma_m, \delta_p, \delta_a)$ (Ch 1, Def. 1.2.4) | Ch 1 |
| $\epsilon_{\mathrm{mach}}$ | Machine epsilon of the arithmetic type | Ch 1 |
| $r$ | Convergence ratio: upper bound on $\lvert a_{k+1}/a_k \rvert$ | §5.3 |
| $q_0$ | Auxiliary quantity for the level ellipsoid: $q_0 = \tfrac{1}{2}[(1 + 3/e'^2)\arctan e' - 3/e']$ | §5.5 |
| $q_0'$ | Derivative $\mathrm{d}q_0/\mathrm{d}e'$, expanded as a series | §5.5 |
| $e'$ | Second eccentricity of the reference ellipsoid | §5.5 |
| $e^2$ | First eccentricity squared of the reference ellipsoid | §5.6 |
| $\binom{\alpha}{k}$ | Generalized binomial coefficient (Ch 4, Def. 4.10.1) | Ch 4 |
| $W_n$ | Wallis integral $\int_0^{\pi/2} \cos^n\phi\;\mathrm{d}\phi$ (formally defined in Ch 6) | Ch 6 |

---

## §5.1 Introduction

Chapter 4 developed the abstract toolkit — Taylor series with remainder bounds, continued fractions, Padé approximants, Chebyshev economization, convergence acceleration, and the generalized binomial series — together with the decision framework of §4.11 for selecting among them. This chapter turns from theory to practice: it applies those methods to the specific series that arise throughout satellite orbit propagation.

Three observations guide the chapter:

1. **Three evaluator types cover all cases.** Every series in this textbook falls into one of three categories: alternating series (§5.2), geometric-tail series (§5.3), or finite polynomials evaluated via Horner's method (§5.4). Each category has a distinct, rigorous truncation bound.

2. **Truncation error is precision error, not accuracy error.** When a series evaluates a known mathematical function (the arctan series for $q_0$, the binomial series for $(1-e^2\sin^2\phi)^{-1/2}$), the omitted terms represent a gap between the computed value and the exact mathematical result — a precision error $\delta_p$ in the classification of Theorem 1.5.1. The single exception is the ephemeris truncation of §5.7.4, where omitted trigonometric terms represent physical effects.

3. **Every evaluator returns a tracked value.** The result is a TrackedValue $\mathcal{V} = (v, \sigma_m, \delta_p, \delta_a)$ whose precision component $\delta_p$ includes both the arithmetic rounding accumulated through the partial sum and the rigorous truncation bound for the omitted tail. The three error categories remain separable: the truncation bound is added only to $\delta_p$, leaving $\sigma_m$ and $\delta_a$ to carry measurement and model errors from upstream.

The chapter proceeds from the three generic evaluators (§§5.2–5.4) through concrete applications: the $q_0$ and $q_0'$ auxiliary series for the level ellipsoid (§5.5), the family of geodetic binomial series (§5.6), a catalog of further series appearing in later chapters (§5.7), and the synthesis connecting these evaluators to the tolerance parameter and the matched-pair principle (§5.8).

---

## §5.2 Alternating Series Evaluation

An alternating series is one whose terms change sign at every step: $S = \sum_{k=0}^{\infty} (-1)^k b_k$ with $b_k > 0$. When the terms are also monotonically decreasing, the Leibniz alternating series test guarantees convergence and provides a remainder bound that is both simple and sharp. This section formalizes the evaluation procedure and its integration with the error-tracking framework.

**Definition 5.2.1** (Alternating series evaluator)**.** *An alternating series evaluator is a procedure that takes:*

- *a term function $t: \mathbb{N} \to \mathcal{V}$ producing the $k$-th term as a tracked value, where terms alternate in sign and $|t(k+1).v| \leq |t(k).v|$,*
- *a start index $k_0 \in \mathbb{N}$,*
- *a tolerance $\tau > 0$ (bare scalar, not a tracked value),*
- *a maximum term count $N_{\max}$,*

*and returns $\mathcal{V}_S = S_{N} = \sum_{k=k_0}^{N} t(k)$, where $N$ is the smallest index satisfying $|t(N).v| < \tau$, with the truncation error $|t(N).v|$ added to the precision component of $\mathcal{V}_S$.*

**Theorem 5.2.1** (Leibniz truncation bound for tracked evaluation)**.** *Let $S = \sum_{k=k_0}^{\infty} t(k)$ be an alternating series satisfying the hypotheses of Theorem 1.5.2 (terms alternate in sign, $|t(k).v|$ is monotonically decreasing, $|t(k).v| \to 0$). The evaluator of Definition 5.2.1 returns a tracked value $\mathcal{V}_S$ satisfying*

$$
|S - \mathcal{V}_S.v| \leq |t(N).v| \tag{5.1}
$$

*and the total precision error of $\mathcal{V}_S$ is*

$$
\delta_p(\mathcal{V}_S) = \delta_p^{\mathrm{arith}}(S_N) + |t(N).v| \tag{5.2}
$$

*where $\delta_p^{\mathrm{arith}}(S_N)$ is the precision error accumulated through the $N - k_0 + 1$ additions and multiplications in the partial sum.*

*Proof.* Inequality (5.1) is Theorem 1.5.2 applied to the tail $R_N = \sum_{k=N+1}^{\infty} t(k)$. For (5.2): the partial sum $S_N$ is computed by successive additions of tracked values, so $\delta_p^{\mathrm{arith}}(S_N)$ is the precision error accumulated through the arithmetic (Ch 1, §1.3, propagation through addition). The truncation remainder $|R_N| \leq |t(N).v|$ is an independent source of precision error (Theorem 1.5.1, first case). Since both are non-negative upper bounds, their sum bounds the total precision gap between $\mathcal{V}_S.v$ and the exact series value $S$. ∎

**Corollary 5.2.1** (Error category preservation)**.** *The truncation bound $|t(N).v|$ is added exclusively to the precision component $\delta_p$. The measurement error $\sigma_m$ and accuracy error $\delta_a$ of the result inherit from the input terms via the addition rules of Ch 1, §1.3. No reclassification of error occurs.*

*Proof.* Immediate from Theorem 1.5.1: the series evaluates a known mathematical function, so truncation is a precision error. ∎

**Proposition 5.2.1** (Term count estimate)**.** *Suppose the convergence ratio $r = |t(k+1).v / t(k).v|$ is approximately constant. Then the evaluator stops after approximately*

$$
N \approx k_0 + \left\lceil \frac{\ln(\tau / |t(k_0).v|)}{\ln r} \right\rceil \tag{5.3}
$$

*terms.*

*Proof.* After $N - k_0$ terms at ratio $r$, the term magnitude is approximately $|t(k_0).v| \cdot r^{N-k_0}$. Setting this below $\tau$ and solving for $N$ gives (5.3). ∎

**Example 5.2.1.** For the $2q_0$ series (Eq. 5.11) with $r \approx e'^2 \approx 0.0067$ (WGS84), $|t(1).v| = 4/15 \cdot e'^3 \approx 1.48 \times 10^{-4}$, and $\ln r \approx \ln 0.0067 \approx -5.0$, Proposition 5.2.1 gives:

| Precision level | $\tau$ | Terms needed |
|-----------------|--------|-------------|
| Double (~16 digits) | $10^{-16}$ | ~7 |
| Extended (~34 digits) | $10^{-34}$ | ~16 |
| 50-digit | $10^{-50}$ | ~23 |

**Remark** (Monotonicity verification)**.** The Leibniz bound (5.1) requires $|t(k+1).v| \leq |t(k).v|$ for all $k$ in the evaluated range. For series with exact rational coefficients times powers of a small parameter $x$ with $|x| < 1$, this holds automatically once $k$ is large enough that $|c_{k+1}/c_k| \cdot |x| < 1$, where $c_k$ are the rational coefficients. For the series in §§5.5–5.6, this condition holds from the first term. A robust implementation can verify monotonicity at runtime and flag violations. See [P.5.1].

---

## §5.3 Geometric-Tail Series Evaluation

Many series in orbit propagation are not alternating but have geometric-like decay: each term is bounded by a constant fraction $r$ of the previous term. The geometric tail bound provides a rigorous remainder estimate for such series.

**Definition 5.3.1** (Geometric-tail series evaluator)**.** *A geometric-tail series evaluator takes:*

- *a term function $t: \mathbb{N} \to \mathcal{V}$,*
- *a start index $k_0$,*
- *a tolerance $\tau > 0$,*
- *a convergence ratio bound $r \in (0, 1)$: an upper bound on $|t(k+1).v / t(k).v|$ for $k \geq k_0$,*
- *a maximum term count $N_{\max}$,*

*and returns $\mathcal{V}_S = S_N$ where $N$ is the smallest index satisfying $|t(N).v| \cdot r / (1-r) < \tau$, with the tail bound $|t(N).v| \cdot r/(1-r)$ added to $\delta_p$.*

**Theorem 5.3.1** (Geometric tail bound for tracked evaluation)**.** *Let $S = \sum_{k=k_0}^{\infty} t(k)$ with $|t(k+1).v| \leq r |t(k).v|$ for all $k \geq k_0$ and some $r \in (0,1)$. The evaluator of Definition 5.3.1 returns a tracked value $\mathcal{V}_S$ satisfying*

$$
|S - \mathcal{V}_S.v| \leq |t(N).v| \cdot \frac{r}{1-r} \tag{5.4}
$$

*and*

$$
\delta_p(\mathcal{V}_S) = \delta_p^{\mathrm{arith}}(S_N) + |t(N).v| \cdot \frac{r}{1-r}. \tag{5.5}
$$

*Proof.* Inequality (5.4) is Theorem 1.5.3 (Eq. 1.22) applied with $|a_{N+1}| \leq r|t(N).v|$. The denominator $1-r$ sums the geometric progression of the tail. The integration into $\delta_p$ follows the same reasoning as Theorem 5.2.1. ∎

**Corollary 5.3.1** (Comparison with Leibniz)**.** *For an alternating series with convergence ratio $r$, the Leibniz bound gives $|R_N| \leq |t(N).v|$, while the geometric bound gives $|R_N| \leq |t(N).v| \cdot r/(1-r)$. The geometric bound is tighter when $r < 1/2$ (since $r/(1-r) < 1$), while the Leibniz bound is tighter when $r > 1/2$.*

*Proof.* The ratio of bounds is $r/(1-r)$, which equals 1 at $r = 1/2$, is less than 1 for $r < 1/2$, and exceeds 1 for $r > 1/2$. ∎

**Remark.** For the geodetic series of §5.6 where $r \approx e^2 \approx 0.0067$ (WGS84), the geometric bound is a factor of $0.0067/(1-0.0067) \approx 0.00674$ tighter than Leibniz — nearly 150 times smaller. This is the practical case: whenever the convergence ratio is known and small, the geometric bound should be preferred.

**Theorem 5.3.2** (Ratio test convergence detection)**.** *For a power series $\sum c_k x^k$ with $\lim_{k\to\infty} |c_{k+1}/c_k| = 1/R$, the asymptotic convergence ratio is $|x|/R$. The evaluator can use an adaptive ratio estimate*

$$
r_k = |t(k).v / t(k-1).v| \tag{5.6}
$$

*at each step, with the current $r_k$ applied in Definition 5.3.1. The resulting tail bound is rigorous provided $r_k$ is indeed an upper bound on all subsequent ratios, which holds when $|c_{k+1}/c_k|$ is eventually monotone (as it is for the binomial series, the Bessel function series, and the hypergeometric series).*

*Proof.* The asymptotic convergence ratio follows from the root test: $\limsup |c_k x^k|^{1/k} = |x|/R$ when $|c_k|^{1/k} \to 1/R$. For the binomial series $(1+x)^\alpha$ (Ch 4, Theorem 4.10.1), $|c_{k+1}/c_k| = |(\alpha-k)/(k+1)|$, which is eventually monotonically decreasing, so $r_k \geq r_{k+j}$ for all $j > 0$ once $k > |\alpha|$. The same holds for any hypergeometric series ${}_pF_q$ by Gauss's formula for the ratio of successive terms. ∎

**Remark** (Adaptive versus fixed ratio)**.** When the convergence ratio $r$ is known a priori (as for the binomial series at a known argument), the fixed ratio of Definition 5.3.1 is simplest. When $r$ is uncertain or argument-dependent, the adaptive ratio (5.6) provides a tighter bound at each step but introduces a dependence on the current and previous terms. See [P.5.4] for the case where adaptive estimation fails.

---

## §5.4 Horner Evaluation of Truncated Series

Once the truncation index $N$ has been determined (by §5.2 or §5.3), the partial sum $S_N = \sum_{k=0}^{N} c_k x^k$ can be evaluated by Horner's method. This section formalizes the combined procedure: determine $N$, then evaluate $S_N$ via Horner with full error tracking.

**Theorem 5.4.1** (Horner evaluation with tracked errors)**.** *Let $p(x) = c_0 + c_1 x + c_2 x^2 + \cdots + c_N x^N$ with each $c_k$ and $x$ given as tracked values. The Horner evaluation*

$$
p(x) = c_0 + x(c_1 + x(c_2 + \cdots + x \cdot c_N)\cdots) \tag{5.7}
$$

*requires $N$ multiplications and $N$ additions. At each step, the error propagation rules of Ch 1, §1.3 (Theorems 1.3.1–1.3.2) apply: multiplication propagates errors via $|f| \cdot \delta + |g| \cdot \delta$, and addition propagates via $\delta + \delta$. The total precision error of the result satisfies*

$$
\delta_p(p(x)) \leq N \cdot \epsilon_{\mathrm{mach}} \cdot \tilde{p}(|x|) + \sum_{k=0}^{N} |x|^k \cdot \delta_p(c_k) \tag{5.8}
$$

*where $\tilde{p}(|x|) = \sum_{k=0}^{N} |c_k.v| \cdot |x.v|^k$ is the condition polynomial (Ch 4, Theorem 4.9.1).*

*Proof.* The Horner recurrence is $h_N = c_N$, $h_k = h_{k+1} \cdot x + c_k$ for $k = N-1, \ldots, 0$, with $p(x) = h_0$. Each multiplication $h_{k+1} \cdot x$ introduces relative rounding error $\epsilon_{\mathrm{mach}}$ on the product (Ch 1, §1.3), contributing at most $\epsilon_{\mathrm{mach}} |h_{k+1}.v| \cdot |x.v|$ to the absolute error. Summing over all $N$ multiplications and using $|h_k.v| \leq \tilde{p}_k(|x|)$ gives the $N \epsilon_{\mathrm{mach}} \tilde{p}(|x|)$ term. The coefficient error term arises from Ch 1, Theorem 1.4.1 (error propagation through a polynomial). ∎

**Corollary 5.4.1** (Exact rational coefficient advantage)**.** *When each coefficient $c_k$ is an exact rational number (Ch 4, Theorem 4.10.2), its tracked value has $\delta_p(c_k) = 0$ (the representation error is zero for a rational stored as an integer ratio). The second term in (5.8) vanishes, leaving*

$$
\delta_p(p(x)) \leq N \cdot \epsilon_{\mathrm{mach}} \cdot \tilde{p}(|x|). \tag{5.9}
$$

*The only source of rounding is the $N$ multiplications by $x$.*

**Proposition 5.4.1** (Combined truncation and evaluation procedure)**.** *The complete procedure for evaluating a convergent power series $S = \sum_{k=0}^{\infty} c_k x^k$ to tolerance $\tau$ is:*

1. *Determine the truncation index $N$ such that $|R_N| < \tau$, using:*
   - *Theorem 5.2.1 (Leibniz) if the series alternates and the terms decrease monotonically, or*
   - *Theorem 5.3.1 (geometric tail) if the convergence ratio $r$ is known.*
2. *Evaluate $S_N = \sum_{k=0}^{N} c_k x^k$ via Horner (Theorem 5.4.1).*
3. *Add the truncation bound to the precision component: $\delta_p = \delta_p^{\mathrm{Horner}} + |R_N|$.*

*The result is a tracked value whose precision component accounts for both computational rounding and mathematical truncation.*

**Remark** (Well-conditioning for small arguments)**.** For the series of §§5.5–5.6, the argument is $x = e'^2 \approx 0.0067$ or $x = e^2 \approx 0.0067$ (WGS84). The condition number $\tilde{p}(|x|)/|p(x)|$ is approximately 1 (the polynomial is dominated by its leading term), so Horner evaluation introduces no significant cancellation error. This remains true for any geodetic eccentricity ($e \leq 0.1$). The problematic regime — alternating coefficients with $|x| \approx 1$ — does not arise in orbit propagation. See Ch 4, [P.4.4] for the general case.

---

## §5.5 Application: The $q_0$ and $q_0'$ Auxiliary Series

The quantities $q_0$ and $q_0'$ arise in the theory of the level (equipotential) ellipsoid: $q_0$ appears in the expression for the normal potential $U_0$ on the ellipsoidal surface, and the ratio $q_0'/q_0$ governs the distribution of normal gravity between equator and pole. Their closed forms involve the combination $\arctan(e') - 3/e'$ which, for small $e'$, produces catastrophic cancellation. The series forms derived here avoid this cancellation entirely.

### §5.5.1 The $q_0$ Series

**Definition 5.5.1** (The quantity $q_0$)**.** *The quantity $q_0$ is defined by*

$$
q_0 = \frac{1}{2}\left[\left(1 + \frac{3}{e'^2}\right)\arctan e' - \frac{3}{e'}\right] \tag{5.10}
$$

*where $e' = \sqrt{e^2/(1-e^2)}$ is the second eccentricity of the reference ellipsoid. This appears in the expression for the normal potential $U_0$ on the level ellipsoid [Heiskanen and Moritz 1967, §2-7] and in the derivation of $J_{2n}$ from the ellipsoidal flattening (Ch 14).*

**Theorem 5.5.1** (Series expansion for $q_0$)**.** *For $|e'| < 1$,*

$$
2q_0 = \sum_{n=1}^{\infty} (-1)^{n+1} \frac{4n}{(2n+1)(2n+3)} e'^{2n+1}. \tag{5.11}
$$

*The series alternates in sign, the terms decrease monotonically in absolute value for $|e'| < 1$, and the convergence ratio is $r = e'^2$.*

*Proof.* Begin with the Maclaurin series for $\arctan$:

$$
\arctan e' = \sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} e'^{2k+1}, \qquad |e'| \leq 1. \tag{5.12}
$$

Multiply by $(1 + 3/e'^2)$:

$$
\left(1 + \frac{3}{e'^2}\right)\arctan e' = \sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} e'^{2k+1} + 3\sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} e'^{2k-1}. \tag{5.13}
$$

Shift the index in the second sum by setting $j = k-1$:

$$
3\sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} e'^{2k-1} = \frac{3}{e'} + 3\sum_{j=0}^{\infty} \frac{(-1)^{j+1}}{2j+3} e'^{2j+1}. \tag{5.14}
$$

Subtract $3/e'$:

$$
\left(1 + \frac{3}{e'^2}\right)\arctan e' - \frac{3}{e'} = \sum_{k=0}^{\infty} \frac{(-1)^k}{2k+1} e'^{2k+1} + 3\sum_{j=0}^{\infty} \frac{(-1)^{j+1}}{2j+3} e'^{2j+1}. \tag{5.15}
$$

Both sums start at index 0. Combining them under a single index $n$:

$$
= \sum_{n=0}^{\infty} (-1)^n \left[\frac{1}{2n+1} - \frac{3}{2n+3}\right] e'^{2n+1}. \tag{5.16}
$$

The bracket simplifies:

$$
\frac{1}{2n+1} - \frac{3}{2n+3} = \frac{(2n+3) - 3(2n+1)}{(2n+1)(2n+3)} = \frac{-4n}{(2n+1)(2n+3)}. \tag{5.17}
$$

Therefore:

$$
\left(1 + \frac{3}{e'^2}\right)\arctan e' - \frac{3}{e'} = \sum_{n=0}^{\infty} (-1)^n \cdot \frac{-4n}{(2n+1)(2n+3)} \cdot e'^{2n+1}. \tag{5.18}
$$

The $n=0$ term vanishes (numerator $-4 \cdot 0 = 0$). For $n \geq 1$, pulling the minus sign into the alternation gives $(-1)^n \cdot (-4n) = (-1)^{n+1} \cdot 4n$. Dividing by 2 (from the definition $2q_0 = \cdots$) yields (5.11).

For the convergence properties: the $n$-th term has magnitude $\frac{4n}{(2n+1)(2n+3)} e'^{2n+1}$. The ratio of successive terms is

$$
\frac{t_{n+1}}{t_n} = \frac{4(n+1)}{(2n+3)(2n+5)} \cdot \frac{(2n+1)(2n+3)}{4n} \cdot e'^2 = \frac{(n+1)(2n+1)}{n(2n+5)} \cdot e'^2. \tag{5.19}
$$

The rational factor approaches 1 as $n \to \infty$ (from above for $n \geq 1$), so the asymptotic convergence ratio is $e'^2$. For $e' < 1$, the series converges, and for $n \geq 1$ with $e'^2 < 1$, the terms decrease monotonically. ∎

**Corollary 5.5.1** (Truncation bound for $q_0$)**.** *By Theorem 5.2.1 (Leibniz bound), stopping the series (5.11) at term $N$ gives*

$$
|R_N| \leq \frac{4N}{(2N+1)(2N+3)} e'^{2N+1}. \tag{5.20}
$$

*For WGS84 ($e' \approx 0.08182$):*

| $N$ | $\lvert R_N \rvert$ bound |
|-----|--------------------------|
| 3 | $< 5 \times 10^{-9}$ |
| 5 | $< 4 \times 10^{-14}$ |
| 7 | $< 3 \times 10^{-19}$ |
| 10 | $< 2 \times 10^{-27}$ |

**Corollary 5.5.2** (Cancellation avoidance)**.** *The closed form (5.10) computes $q_0$ as the difference of two nearly equal quantities when $e'$ is small. Specifically, for WGS84:*

$$
(1 + 3/e'^2)\arctan e' \approx 36.543, \qquad 3/e' \approx 36.543, \tag{5.21}
$$

*but the two quantities agree to 5 leading digits ($36.5434\ldots$ vs $36.5433\ldots$). The difference $2q_0 \approx 1.47 \times 10^{-4}$ is five orders of magnitude smaller than either operand, so the cancellation costs approximately $\log_{10}(36.5 / 1.47 \times 10^{-4}) \approx 5.4$ significant digits. At double precision (16 digits), this leaves $\sim$10 reliable digits — adequate but wasteful. At higher precision, the absolute digit loss is the same, but it consumes a smaller fraction of the available digits.*

*The series form (5.11) avoids this cancellation entirely: the $3/e'$ and the corresponding $n=0$ term cancel algebraically in the derivation (Eq. 5.17), leaving only terms proportional to $e'^{2n+1}$ with $n \geq 1$. No catastrophic subtraction occurs.* [P.5.2]

### §5.5.2 The $q_0'$ Series

**Definition 5.5.2** (The quantity $q_0'$)**.** *The quantity $q_0'$ is the shape factor derivative arising in the expressions for normal gravity. Its closed form (Ch 14, Definition 14.4.3) is*

$$
q_0' = 3\left(1 + \frac{1}{e'^2}\right)\left(1 - \frac{1}{e'}\arctan e'\right) - 1. \tag{5.22}
$$

*This appears in the Somigliana formula for normal gravity and in the relation between $J_2$ and the flattening (Ch 14, §14.7).*

[Heiskanen and Moritz 1967, §2-7; Moritz 1980, p. 130]

**Theorem 5.5.2** (Series expansion for $q_0'$)**.** *For $|e'| < 1$,*

$$
q_0' = \sum_{n=1}^{\infty} (-1)^{n+1} \frac{6}{(2n+1)(2n+3)} e'^{2n}. \tag{5.23}
$$

*The series alternates, the terms decrease monotonically for $|e'| < 1$, and the convergence ratio approaches $e'^2$.*

*Proof.* The equivalence of the closed form (5.22) and the series (5.23) is proved in Ch 14, Theorem 14.4.3. The derivation expands $\arctan(e')/e'$ as a Maclaurin series, subtracts from 1, multiplies by $(1 + 1/e'^2)$, and collects terms at each power of $e'^{2n}$. The $n = 0$ term vanishes after the outer subtraction of 1, and the $n \geq 1$ coefficients combine to $(-1)^{n+1} \cdot 6/((2n+1)(2n+3))$ via the identity

$$
(-1)^n\left[\frac{-1}{2n+1} + \frac{1}{2n+3}\right] = (-1)^n \cdot \frac{-2}{(2n+1)(2n+3)}. \tag{5.24}
$$

Multiplying by the outer factor of 3 yields the series (5.23). The convergence ratio $|t_{n+1}/t_n| = (2n+1)/(2n+5) \cdot e'^2 \to e'^2$ follows from the ratio of successive terms. ∎

**Corollary 5.5.3** (Structural comparison)**.** *The $q_0$ and $q_0'$ series share a common structure:*

| Property | $2q_0$ (Eq. 5.11) | $q_0'$ (Eq. 5.23) |
|----------|------|------|
| General term | $(-1)^{n+1}\frac{4n}{(2n+1)(2n+3)} e'^{2n+1}$ | $(-1)^{n+1}\frac{6}{(2n+1)(2n+3)} e'^{2n}$ |
| Starts at | $n = 1$ | $n = 1$ |
| Power of $e'$ | Odd: $2n+1$ | Even: $2n$ |
| Numerator | $4n$ (grows linearly) | $6$ (constant) |
| Denominator | $(2n+1)(2n+3)$ | $(2n+1)(2n+3)$ |
| Conv. ratio | $e'^2$ | $e'^2$ |

*Both share the same denominator pattern and convergence ratio. The $q_0'$ series converges slightly faster because its numerator is constant rather than growing linearly with $n$.*

### §5.5.3 The $U_0$ Auxiliary Series

The normal potential $U_0$ on the ellipsoidal surface (Ch 14) requires the evaluation of $\arctan(e')/e'$, which is the third member of the arctan-derived family alongside $q_0$ and $q_0'$.

**Theorem 5.5.3** (Series expansion for $\arctan(e')/e'$)**.** *For $|e'| < 1$,*

$$
\frac{\arctan e'}{e'} = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n} = 1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \frac{e'^6}{7} + \cdots \tag{5.28}
$$

*The series alternates, the terms decrease monotonically, and the convergence ratio is $r = e'^2$.*

*Proof.* Divide the Maclaurin series (5.12) by $e'$. The ratio of successive terms is $|t_{n+1}/t_n| = \frac{2n+1}{2n+3} e'^2 \to e'^2$. ∎

**Remark.** The $U_0$ computation in the ellipsoid constructor uses this series in the form $U_0 = (GM/b) \cdot S + \omega^2 a^2/3$, where $S = \sum_{n=0}^{\infty} (-1)^n e'^{2n}/(2n+1)$. The coefficients $1/(2n+1)$ are simpler than those in the $q_0$ and $q_0'$ series — they are direct reciprocals of odd integers, with no $(2n+3)$ factor in the denominator. The convergence ratio is the same $e'^2$, and the Leibniz bound applies.

### §5.5.4 Horner Nesting

**Proposition 5.5.1** (Horner form for $q_0$)**.** *The series (5.11) can be written in Horner-nested form by factoring out $e'^3$:*

$$
2q_0 = e'^3 \sum_{n=1}^{\infty} (-1)^{n+1} \frac{4n}{(2n+1)(2n+3)} (e'^2)^{n-1} = e'^3 \left[\frac{4}{3 \cdot 5} - e'^2 \cdot \frac{8}{5 \cdot 7} + e'^4 \cdot \frac{12}{7 \cdot 9} - \cdots\right]. \tag{5.25}
$$

*Setting $u = -e'^2$ and $c_n = 4(n+1)/((2n+3)(2n+5))$ for $n = 0, 1, 2, \ldots$, the bracketed expression is $\sum_{n=0}^{N-1} c_n u^n$, evaluable via Horner (Theorem 5.4.1) with exact rational coefficients (Ch 4, Theorem 4.10.2):*

$$
c_0 = \frac{4}{15}, \quad c_1 = \frac{8}{35}, \quad c_2 = \frac{12}{63} = \frac{4}{21}, \quad c_3 = \frac{16}{99}. \tag{5.26}
$$

*Each $c_n$ is a ratio of small integers, representable exactly.*

*Similarly, $q_0' = e'^2[c_0' + c_1' \cdot (-e'^2) + c_2' \cdot e'^4 - \cdots]$ with $c_n' = 6/((2n+3)(2n+5))$ — again exact rationals.*

### §5.5.5 Continued Fraction Alternative

**Proposition 5.5.2** (Continued fraction for $q_0$)**.** *Since $q_0$ is expressed through $\arctan(e')$ (Definition 5.5.1), the Gauss continued fraction for $\arctan$ (Ch 4, §4.4, Proposition 4.4.2) provides an alternative evaluation path:*

$$
\arctan z = \cfrac{z}{1 + \cfrac{z^2}{3 + \cfrac{4z^2}{5 + \cfrac{9z^2}{7 + \cdots}}}}. \tag{5.27}
$$

*For the geodetic application, however, the series (5.11) is preferred over the continued fraction for three reasons:*

1. *Error bound transparency: the Leibniz bound (Theorem 5.2.1) gives a simple, sharp remainder estimate directly from the last term. The continued fraction convergence bound (Ch 4, Theorem 4.4.3) is tighter asymptotically but requires tracking the convergent denominators $B_n$.*

2. *Cancellation avoidance: the derivation of (5.11) eliminates the catastrophic cancellation in the closed form (5.10). The continued fraction evaluates $\arctan(e')$ directly, and the subtraction $(1+3/e'^2)\arctan e' - 3/e'$ would reintroduce the same cancellation.*

3. *Performance at small arguments: for $|e'| \approx 0.08$ (WGS84), the series converges in $\sim$7 terms at double precision. The continued fraction converges in $\sim$5–6 convergents, offering no significant advantage in a context where the evaluation is performed once at initialization (Ch 3, §3.4).*

*The continued fraction becomes advantageous only when $|e'|$ approaches 1 (eccentricities far beyond any terrestrial application), where the series convergence degrades while the continued fraction maintains geometric convergence. At $\tau = \tau_{\mathrm{standard}}$, both methods produce identical results to well within the tolerance.*

### §5.5.6 Implementation Correspondence

The series evaluation of $q_0$ and $q_0'$ is implemented in `EquipotentialEllipsoid` (see Appendix C). The implementation uses:

- `ratio<T>(sign * 4 * n, (2*n+1) * (2*n+3))` for exact rational coefficients matching (5.11),
- `alternating_series<T>(1, q0_term, series_tolerance)` for Leibniz-bounded evaluation matching Theorem 5.2.1,
- Repeated multiplication `ep_power = ep_power * e_prime` for the power $e'^{2n+1}$, avoiding `pow` and its associated precision loss.

The series tolerance parameter is the same bare-`T` quantity that governs all series evaluations in the ellipsoid constructor — a single control point matching the tolerance parameter architecture of Chapter 3, Definition 3.4.1.

---

## §5.6 Application: Geodetic Binomial Series

The gravitational field of the reference ellipsoid generates a family of integrals involving $(1 + e'^2\cos^2\Phi)^{\alpha}$ for half-integer exponents $\alpha$. Each such integral expands as a binomial series integrated term-by-term against a Wallis weight, producing the published geodetic coefficients. This section formalizes the three-stage derivation pattern and its error bounds.

### §5.6.1 The Geodetic Binomial Family

**Definition 5.6.1** (Geodetic binomial series)**.** *A geodetic binomial series is a power series arising from the binomial expansion*

$$
(1 + e'^2\cos^2\Phi)^{\alpha} = \sum_{k=0}^{\infty} \binom{\alpha}{k} (e'^2\cos^2\Phi)^k \tag{5.29}
$$

*evaluated at a specific half-integer exponent $\alpha$ and integrated over latitude. This is the canonical form adopted in Ch 14, Definition 14.8.1, following Moritz (1980). The binomial argument $e'^2\cos^2\Phi$ is positive, so the expansion carries no additional $(-1)^k$ factor — the alternating signs in the published coefficients arise from $\binom{\alpha}{k}$ itself when $\alpha < 0$.*

*The alternative form $(1 - e^2\sin^2\Phi)^\alpha$ is equivalent via Corollary 14.3.1 (fundamental identity) and differs by the constant factor $(1-e^2)^\alpha$.*

*The exponents and their applications are:*

| Exponent $\alpha$ | Application | Appears in |
|---|---|---|
| $-1/2$ | Radius of curvature in the prime vertical | Ch 14 |
| $-3/2$ | Meridian arc integrand | Ch 14 |
| $-2$ | Equal-area sphere integrand | Ch 14 |
| $-5/2$ | Mean gravity numerator | Ch 14 |
| $+1/2$ | Square root of a series (series algebra step) | Ch 14 |

**Theorem 5.6.1** (Convergence for geodetic parameters)**.** *For any reference ellipsoid with $e'^2 < 1$ (satisfied by all terrestrial and planetary bodies), the binomial series (5.29) converges absolutely for all latitudes $\Phi$, since $e'^2\cos^2\Phi \leq e'^2 < 1$. By Theorem 4.10.1, the remainder after $N$ terms satisfies*

$$
|R_N(\Phi)| \leq \left|\binom{\alpha}{N+1}\right| \frac{e'^{2(N+1)}}{1 - e'^2}. \tag{5.30}
$$

*Proof.* The convergence follows from Ch 4, Theorem 4.10.1 with $|x| = e'^2\cos^2\Phi \leq e'^2$. The bound (5.30) uses the worst case $\cos^2\Phi = 1$ and the geometric tail factor $1/(1-e'^2)$ from Ch 1, Theorem 1.5.3. ∎

**Table 5.6.1** (Remainder bounds at $N = 4$ for WGS84, $e'^2 \approx 0.00674$)**.**

| $\alpha$ | $\lvert\binom{\alpha}{5}\rvert$ | $|R_4|$ bound |
|---|---|---|
| $-1/2$ | $7/256$ | $< 2 \times 10^{-13}$ |
| $-3/2$ | $693/256$ | $< 2 \times 10^{-11}$ |
| $-2$ | $6$ | $< 4 \times 10^{-11}$ |
| $-5/2$ | $45045/256$ | $< 1 \times 10^{-9}$ |

At double precision, $N = 4$ or $5$ terms suffice for all exponents. At 50-digit precision, $N \approx 23$ terms are needed (Proposition 5.2.1).

### §5.6.2 The Three-Stage Coefficient Derivation

The published geodetic series coefficients — the numbers that appear in formulas for the meridian arc $Q$, the equal-area radius $R_2$, mean gravity $\bar{\gamma}$, and the gravity formula $\gamma(\Phi)$ — all arise from the same three-stage process. This subsection formalizes the pattern.

**Theorem 5.6.2** (Three-stage coefficient derivation)**.** *Let $I(\alpha, w) = \int_0^{\pi/2} (1 + e'^2\cos^2\Phi)^{\alpha} w(\Phi)\;\mathrm{d}\Phi$ where $w(\Phi)$ is a latitude weight function ($w = 1$, $w = \cos\Phi$, etc.), following the canonical form of Ch 14, Definition 14.8.1. The coefficient of $e'^{2k}$ in the series expansion of $I$ is computed in three stages:*

*Stage 1 (Binomial expansion): Expand the integrand $(1 + x)^{\alpha}$ where $x$ is a positive quantity ($e'^2\cos^2\Phi$ or $e'^2$), obtaining coefficients $\binom{\alpha}{k}$ for each power $x^k$. The alternating signs in the published coefficients arise from $\binom{\alpha}{k}$ itself when $\alpha$ is negative — no additional $(-1)^k$ factor is needed.*

*Stage 2 (Wallis integration): Integrate each term against $w(\Phi)$ over $[0, \pi/2]$, producing Wallis-type integrals (Ch 6):*
- *For $w(\Phi) = 1$: $W_{2k} = \int_0^{\pi/2}\cos^{2k}\Phi\;\mathrm{d}\Phi = \int_0^{\pi/2}\sin^{2k}\Phi\;\mathrm{d}\Phi = \frac{(2k-1)!!}{(2k)!!}\cdot\frac{\pi}{2}$ (even Wallis).*
- *For $w(\Phi) = \cos\Phi$: $W_{2k+1} = \int_0^{\pi/2}\cos^{2k+1}\Phi\;\mathrm{d}\Phi = \int_0^{\pi/2}\sin^{2k+1}\Phi\;\mathrm{d}\Phi = \frac{(2k)!!}{(2k+1)!!}$ (odd Wallis).*

*Stage 3 (Series algebra): If the physical quantity involves a composition (square root, reciprocal, quotient) of integrated series, apply the corresponding algebraic operation on the coefficient sequences.*

*At each stage, the coefficients are exact rational numbers: the binomial coefficients are rational (Ch 4, Theorem 4.10.2 for half-integer $\alpha$), the Wallis integrals are rational multiples of $\pi/2$ or are rational outright, and the series algebra operations (Cauchy products, Newton iteration for square root) preserve rationality. The only irrational quantity is $\pi$, which appears as a common factor in even-Wallis integrals and cancels in ratios.*

*Proof.* The rationality of Stage 1 coefficients follows from Ch 4, Theorem 4.10.2. The Wallis integrals are ratios of double factorials (Ch 6), which are products of integers. The Cauchy product of two sequences of rationals is rational. Newton iteration for $\sqrt{1+u}$ with rational $u$ coefficients produces rational coefficients at each step (the iteration involves only addition, multiplication, and division by 2). ∎

**Corollary 5.6.1** (Compute, don't look up)**.** *For arbitrary-precision evaluation, the geodetic coefficients must be computed from the three-stage formula rather than stored as floating-point literals. Storing precomputed double-precision coefficients limits the result to $\sim$16 digits regardless of the arithmetic type $T$.*

**Example 5.6.1.** The meridian arc integral uses the second eccentricity form $(1 + e'^2\cos^2\Phi)^{-3/2}$ with $w(\Phi) = 1$, where the binomial argument $x = e'^2\cos^2\Phi$ is positive. Stage 1 gives $\binom{-3/2}{k}$ (no additional $(-1)^k$ — the alternating signs come from the binomial coefficient itself). Stage 2 multiplies by $W_{2k}/(\pi/2)$. The combined coefficient of $e'^{2k}$ is:

$$
c_k^{(Q)} = \binom{-3/2}{k} \cdot \frac{(2k-1)!!}{(2k)!!}. \tag{5.31}
$$

*The first few values: $c_0 = 1$, $c_1 = -3/4$, $c_2 = 45/64$, $c_3 = -175/256$, $c_4 = 11025/16384$.*

**Example 5.6.2.** The equal-area radius $R_2$ uses the form $(1 + e'^2\cos^2\Phi)^{-2}$ with $w(\Phi) = \cos\Phi$, then takes a square root. Stage 1: $\binom{-2}{k} = (-1)^k(k+1)$ (the alternation is intrinsic to the binomial coefficient). Stage 2: $\int_0^{\pi/2}\cos^{2k+1}\Phi\;\mathrm{d}\Phi = (2k)!!/(2k+1)!!$. The integral is $I = \sum_{k=0}^{\infty} (-1)^k(k+1) \cdot (2k)!!/(2k+1)!! \cdot e'^{2k} = 1 - 4e'^2/3 + 8e'^4/5 - \cdots$. Stage 3: $R_2/c = \sqrt{I}$ via the binomial series with $\alpha = 1/2$, yielding $1 - 2e'^2/3 + 26e'^4/45 - \cdots$.

### §5.6.3 Integrated Series Truncation Bound

**Theorem 5.6.3** (Truncation bound for integrated binomial series)**.** *When the binomial series (5.29) is truncated at $N$ terms and then integrated against a weight $w(\Phi)$, the remainder of the integrated series satisfies*

$$
|R_N^{\mathrm{int}}| \leq \left|\binom{\alpha}{N+1}\right| e'^{2(N+1)} \cdot \begin{cases} W_{2(N+1)}/(\pi/2) & \text{if } w = 1 \\ (2N+2)!!/(2N+3)!! & \text{if } w = \cos\Phi \end{cases} \tag{5.32}
$$

*Both multipliers are at most 1, so the integrated bound is at most the raw binomial bound (5.30). For $w = \cos\Phi$, the odd Wallis factor provides substantial additional tightening.*

*Proof.* The truncation remainder of the integrated series is $\sum_{k=N+1}^{\infty} |\binom{\alpha}{k}| e'^{2k} \cdot F_k$, where $F_k$ is the integral factor. For $w = 1$: $F_k = W_{2k}/(\pi/2) = (2k-1)!!/(2k)!! \leq 1$ with equality at $k = 0$, strictly decreasing (Wallis integral properties, Ch 6). For $w = \cos\Phi$: $F_k = (2k)!!/(2k+1)!! \leq 1$, strictly decreasing. In both cases, $F_k$ tightens the bound beyond the raw binomial remainder. ∎

### §5.6.4 Series Algebra Error Bounds

When the physical quantity involves a composition of integrated series (as in Examples 5.6.1–5.6.2), the series algebra step introduces no new truncation error — it is an exact algebraic operation on finitely many coefficients. However, the arithmetic rounding in the algebra step must be tracked.

**Lemma 5.6.1** (Series square root via Newton iteration)**.** *Let $I = 1 + u_1 x + u_2 x^2 + \cdots + u_N x^N$ be a truncated power series with $|u_1 x + \cdots + u_N x^N| = |U| < 1$. The square root $\sqrt{I}$ can be computed to $N$ terms by the binomial series $\sqrt{1 + U} = \sum_{k=0}^{M} \binom{1/2}{k} U^k$ truncated at $M$ terms, or equivalently by Newton iteration $s_{j+1} = (s_j + I/s_j)/2$ starting from $s_0 = 1$. After $j$ iterations, the error in $s_j$ (viewed as a polynomial in $x$) satisfies*

$$
|s_j - \sqrt{I}| \leq \frac{|U|^{2^j}}{2^{j+1}}. \tag{5.33}
$$

*For $|U| \approx 0.007$ (WGS84 geodetic series) and $j = 4$: error $< 10^{-40}$.*

*Proof.* Newton's method for $f(s) = s^2 - I$ has quadratic convergence. The error $e_j = s_j - \sqrt{I}$ satisfies $e_{j+1} = e_j^2/(2s_j)$. Since $s_j \geq 1$ for the starting value $s_0 = 1$ and $I$ close to 1, $|e_{j+1}| \leq |e_j|^2/2$. Iterating from $|e_0| = |1 - \sqrt{I}| \leq |U|/2$ gives $|e_j| \leq |U|^{2^j}/2^{j+1}$. ∎

**Lemma 5.6.2** (Series quotient error bound)**.** *Let $N = 1 + n_1 x + \cdots + n_M x^M$ and $D = 1 + d_1 x + \cdots + d_M x^M$ with $|D.v| \geq d_{\min} > 0$. The quotient $Q = N/D$ evaluated as a tracked-value division satisfies*

$$
\delta_p(Q) \leq \frac{\delta_p(N) + |Q.v| \cdot \delta_p(D)}{|D.v|} + \epsilon_{\mathrm{mach}} |Q.v|. \tag{5.34}
$$

*The first term propagates input precision errors through the quotient rule (Ch 1, §1.3); the second is the rounding error of the division itself.*

*Proof.* Standard error propagation for $f(N, D) = N/D$: $|\partial f/\partial N| = 1/|D|$, $|\partial f/\partial D| = |N|/|D|^2 = |Q|/|D|$. Applying Ch 1, Theorem 1.4.1. ∎

---

## §5.7 Application: Further Series in Orbit Propagation

This section catalogs the principal series arising in later chapters, classifying each by evaluator type and convergence properties. Full derivations appear in the cited chapters; the purpose here is to establish that the three evaluators of §§5.2–5.4 cover all cases.

### §5.7.1 Kepler's Equation: Bessel Function Series (Ch 9)

**Proposition 5.7.1** (Bessel series for Kepler's equation)**.** *The solution of Kepler's equation $E - e\sin E = M$ can be expressed as the Bessel series*

$$
E = M + \sum_{n=1}^{\infty} \frac{2}{n} J_n(ne) \sin(nM) \tag{5.35}
$$

*where $J_n$ is the Bessel function of the first kind. The terms decay geometrically with ratio approximately $e/2$ for moderate eccentricities. The geometric-tail evaluator (§5.3) applies with $r \approx e/2$.*

*However, for computational purposes, iterative methods (Newton's method, Ch 4 Theorem 4.7.1, developed in Ch 9) are preferred over the series: they converge in 3–4 iterations at double precision for all LEO eccentricities, compared to $\sim$10 terms of the Bessel series, and they avoid evaluating Bessel functions. The series remains valuable for theoretical analysis and for deriving the equation of center (Ch 25).*

### §5.7.2 Hansen Coefficients (Ch 15)

**Proposition 5.7.2** (Hansen coefficient series)**.** *The Hansen coefficients $X_k^{n,m}(e)$ used in the Kaula expansion of the disturbing function (Ch 15) are power series in $e$ with exact rational coefficients:*

$$
X_k^{n,m}(e) = \sum_{j=0}^{\infty} c_j^{(n,m,k)} e^j \tag{5.36}
$$

*where $c_j^{(n,m,k)}$ are rational numbers computable from the Newcomb operators. The convergence ratio is $|e|$, so the geometric-tail evaluator (§5.3) applies with $r = e$. For LEO orbits ($e \leq 0.02$), $N \leq 10$ terms suffice at double precision.*

### §5.7.3 Drag Coefficient Series (Ch 22)

**Proposition 5.7.3** (Orbit-averaged drag series)**.** *The orbit-averaged drag coefficients (Ch 22) are non-alternating power series in eccentricity with polynomial coefficients in the density scale height ratio. The terms are positive and geometrically decaying with ratio $r = e$. The geometric-tail evaluator (§5.3) applies.*

### §5.7.4 Trigonometric Ephemeris Truncation (Ch 25–26)

**Proposition 5.7.4** (Ephemeris truncation as accuracy error)**.** *The solar and lunar ephemerides (Ch 25–26) are finite truncated trigonometric sums — they do not converge to an exact function as more terms are added, because each omitted term represents a distinct physical perturbation with its own amplitude and frequency. The truncation error is therefore an accuracy error $\delta_a$ (Theorem 1.5.1, second case), not a precision error.* [A.5.1]

*This is the sole exception to the pattern of §§5.2–5.6: every other series in this chapter evaluates a known mathematical function, where truncation is $\delta_p$. The distinction matters for the error budget: precision errors are reducible by computing more terms or using higher-precision arithmetic, while accuracy errors require a fundamentally different model (a higher-fidelity ephemeris).*

---

## §5.8 The Series Evaluation Interface

### §5.8.1 Interface Design

**Definition 5.8.1** (Series evaluator interface)**.** *Every series evaluator in this chapter conforms to the following interface:*

- *Input: a term function $t: \mathbb{N} \to \mathcal{V}$, a start index $k_0$, a tolerance $\tau$ (bare scalar of type $T$), and a maximum term count $N_{\max}$.*
- *Output: a tracked value $\mathcal{V}_S$ whose value component $v$ is the partial sum $S_N$, and whose precision component $\delta_p$ includes both arithmetic rounding and the truncation bound.*

*The Horner evaluator (§5.4) differs in taking a coefficient array and argument rather than a term function, but the output contract is the same: a tracked value with precision accounting for all rounding.*

*The tolerance $\tau$ is captured at initialization time (Ch 3, Def. 3.4.1) and stored in the closure that constructs the term function. The evaluator does not know or care about the matched-pair principle — it simply evaluates until $|R_N| < \tau$. The matched-pair compatibility is guaranteed by the choice of $\tau$ (Ch 3, Corollary 3.7.1).*

**Remark** (Init-time constants vs per-step evaluations)**.** The series in this chapter divide into two categories by evaluation lifetime:

1. *Init-time constants:* Series whose arguments depend only on the defining parameters of the model ($e'$, $e^2$, $f$). The $q_0$, $q_0'$, and $U_0$ series (§5.5), the geodetic binomial coefficients (§5.6), and the even zonal harmonics $J_{2n}$ all depend only on the eccentricity of the reference ellipsoid, which is fixed at initialization. These series are evaluated once in the initializer, and the resulting TrackedValues are stored as members — the closure collapses to a constant. No series evaluation occurs per propagation step.

2. *Per-step evaluations:* Series whose arguments depend on the current orbital elements. The Kepler equation (§5.7.1), Hansen coefficients at the current eccentricity (§5.7.2), and drag coefficients (§5.7.3) change at every time step. For these, the initializer captures $\tau$ in a closure and produces a callable $\lambda(e)$ (Ch 3, Example 3.4.1) that evaluates the series at the current eccentricity. The evaluator runs per step, but $\tau$ is invisible to the caller.

The evaluators in §§5.2–5.4 are general-purpose: they accept any term function and tolerance, and return a TrackedValue. They are not specific to orbit propagation. The SGP4-specific logic resides entirely in the term functions and the initialization architecture — the series library is reusable for any application that requires error-bounded series evaluation.

**Remark** (Bound selection for the library)**.** The numerical verification of §5.9.4 shows that the Leibniz bound (Theorem 5.2.1) overestimates the actual truncation error by a factor of $\sim$170 for the $q_0$ series, while the geometric tail bound (Theorem 5.3.1) overestimates by only 10–20%. For a reusable library, this suggests that `alternating_series` should accept an optional convergence ratio parameter: when the ratio is provided, the evaluator uses the tighter geometric bound (5.4) instead of the Leibniz bound (5.1). When the ratio is not provided, the Leibniz bound serves as the conservative default. This avoids reporting a $\delta_p$ that is two orders of magnitude larger than the actual error.

**Theorem 5.8.1** (Tolerance parameter integration)**.** *Let $g(\tau)$ denote the result of a series evaluator called with tolerance $\tau$. Then*

$$
|g(\tau).v - f_{\mathrm{exact}}| \leq \tau + \delta_p^{\mathrm{arith}} \tag{5.37}
$$

*where $f_{\mathrm{exact}}$ is the exact value of the infinite series. At $\tau = \tau_{\mathrm{standard}}$, the matched-pair compatibility criterion (Ch 3, Def. 3.7.1) is satisfied:*

$$
|g(\tau_{\mathrm{standard}}).v - f_{\mathrm{SGP4}}| \leq \tau_{\mathrm{standard}} + \delta_p^{\mathrm{arith}} + \delta_p^{\mathrm{SR3}} \leq \delta_a(\mathcal{M}) \tag{5.38}
$$

*where $\delta_p^{\mathrm{SR3}}$ is the precision error of the Spacetrack Report No. 3 [Hoots and Roehrich 1980] evaluation and $\delta_a(\mathcal{M})$ is the model accuracy (all on the order of km, far above the precision-level quantities).*

*Proof.* Inequality (5.37): the truncation bound is at most $\tau$ (by the stopping criterion), and the arithmetic rounding is $\delta_p^{\mathrm{arith}}$, accumulated through the partial sum. These are independent non-negative bounds. Inequality (5.38): by the triangle inequality, $|g - f_{\mathrm{SGP4}}| \leq |g - f_{\mathrm{exact}}| + |f_{\mathrm{exact}} - f_{\mathrm{SGP4}}| \leq (\tau_{\mathrm{standard}} + \delta_p^{\mathrm{arith}}) + \delta_p^{\mathrm{SR3}}$. All three terms are precision-level quantities (at most $\sim 10^{-10}$), far below $\delta_a(\mathcal{M}) \approx 1$ km, establishing the matched-pair compatibility. This is the series-specific verification required by Ch 3, Corollary 3.7.1, and the general result proved in Ch 4, Proposition 4.11.1. ∎

### §5.8.2 Method Selection by Series

**Table 5.8.1** (Method comparison for series with multiple evaluation paths)**.**

Every series is truncated by §5.2 or §5.3, then evaluated via Horner (Proposition 5.4.1) unless noted. The preferred method (marked $\checkmark$) is selected by convergence rate and error profile together.

| Series | Method | Conv. rate | Error profile | |
|--------|--------|------------|---------------|-|
| $q_0$ | Alternating series (§5.2) | Geometric, $r \to e'^2$ (Eq. 5.19) | No cancellation; Leibniz bound (Thm 5.2.1) | $\checkmark$ |
| | CF for $\arctan$ (5.27) | Geometric, ~5--6 convergents at double (§5.5.4) | Subtractive cancellation in $(1+3/e'^2)\arctan e' - 3/e'$ [P.5.2] | |
| | Closed form (5.10) | $O(1)$ | ~5.4 digits lost to cancellation (Cor. 5.5.2) | |
| $q_0'$ | Alternating series (§5.2) | Geometric, $r \to e'^2$ | No cancellation; Leibniz bound | $\checkmark$ |
| | Closed form (5.22) | $O(1)$ | Same cancellation as $q_0$ | |
| Binomial, $\alpha = -1/2$ | Series (§5.3) | Geometric, $r \to e'^2$ (Thm 5.6.1) | Exact rational coefficients (Cor. 5.4.1); tail bound (5.30) | $\checkmark$ |
| Binomial, $\alpha = -3/2$ | Series (§5.3) | Geometric, $r \to e'^2$ (Thm 5.6.1) | Same as $\alpha = -1/2$ | $\checkmark$ |
| $R_2$ sqrt | Binomial $\alpha = 1/2$ / Newton iteration | Quadratic: $\lvert e_j\rvert \leq \lvert U\rvert^{2^j}/2^{j+1}$ | Equivalent methods (Lemma 5.6.1); inherits input $\delta_p$ | $\checkmark$ |
| $\bar{\gamma}$ quotient | TrackedValue division | N/A (single operation) | $\delta_p$ via quotient rule (Eq. 5.34) | $\checkmark$ |
| Kepler equation | Newton iteration (Ch 9) | Quadratic (Thm 4.7.1) | 3--4 iterations at double; no special functions required | $\checkmark$ |
| | Bessel series (Prop. 5.7.1) | Geometric, $r \approx e/2$ | ~10 terms at double; requires Bessel function evaluation | |
| Hansen $X_k^{n,m}$ | Power series (Prop. 5.7.2) | Geometric, $r = e$ | Exact rational coefficients; arbitrary precision | $\checkmark$ |
| Drag coefficients | Geometric-tail series (Prop. 5.7.3) | Geometric, $r = e$ | Positive terms; tail bound via Thm 5.3.1 | $\checkmark$ |
| Ephemeris | Fixed trig sum (Prop. 5.7.4) | N/A (finite) | Truncation is $\delta_a$ [A.5.1], not $\delta_p$ | $\checkmark$ |

### §5.8.3 Term Counts by Tolerance Level

**Proposition 5.8.1** (Term count as a function of tolerance)**.** *For a geometric series with convergence ratio $r$ and leading term magnitude $|t_0|$, the number of terms required to reduce the truncation bound below a tolerance $\tau$ is approximately $N \approx \lceil \ln(\tau/|t_0|) / \ln r \rceil$ (Proposition 5.2.1, Eq. 5.3). The following table gives $N$ across five tolerance levels: the model accuracy floor $\delta_a(\mathcal{M}) \approx 10^{-3}$, the current implementation tolerance $10^{-12}$, and machine epsilon for double ($10^{-16}$), quad ($10^{-34}$), and 50-digit ($10^{-50}$) arithmetic.*

**Table 5.8.2** (Term counts by tolerance regime, WGS84)**.**

| Series | $r$ | Accuracy floor $10^{-3}$ | Implementation $10^{-12}$ | Double $10^{-16}$ | Quad $10^{-34}$ | 50-digit $10^{-50}$ |
|--------|-----|:--:|:--:|:--:|:--:|:--:|
| $q_0$ | $e'^2 \approx 0.0067$ | 1 | 5 | 7 | 15 | 22 |
| $q_0'$ | $e'^2 \approx 0.0067$ | 1 | 5 | 7 | 15 | 22 |
| Binomial, $\alpha = -1/2$ | $e'^2 \approx 0.0067$ | 1 | 3 | 4 | 8 | 12 |
| Binomial, $\alpha = -3/2$ | $e'^2 \approx 0.0067$ | 1 | 4 | 5 | 10 | 15 |
| Kepler Bessel | $e/2 \approx 0.01$ | 1 | 6 | 8 | 17 | 25 |
| Hansen | $e \approx 0.02$ | 1 | 7 | 10 | 20 | 30 |

*The $\tau = 10^{-3}$ column corresponds to the SGP4 model accuracy floor $\delta_a(\mathcal{M}) \approx 1$ km (Ch 3, [A.3.1]). At this tolerance, every series converges in a single term — the higher-order corrections are smaller than the model error. The $\tau = 10^{-12}$ column reflects the tolerance used in the current implementation (see Appendix C). The three rightmost columns ($10^{-16}$, $10^{-34}$, $10^{-50}$) represent the theoretical maximum useful precision for double, quad, and 50-digit arithmetic respectively.*

*Proof.* Apply Eq. (5.3) with the leading term magnitudes from §§5.5–5.7. For $q_0$: $|t_1| \approx 4/15 \cdot e'^3 \approx 7.8 \times 10^{-5}$, $\ln r \approx -5.0$; $N(\tau = 10^{-3}) = \lceil \ln(10^{-3}/7.8 \times 10^{-5})/(-5.0) \rceil = \lceil 0.5 \rceil = 1$. For the binomial $\alpha = -1/2$: $|t_1| = |\binom{-1/2}{1}| \cdot e^2 = e^2/2 \approx 0.0033$; $N(\tau = 10^{-3}) = \lceil \ln(10^{-3}/0.0033)/(-5.0) \rceil = \lceil 0.24 \rceil = 1$. The remaining entries follow analogously. ∎

**Table 5.8.3** (SR3 truncation orders and matched tolerance)**.**

The Hoots and Roehrich (1980) implementation uses fixed-order polynomial expressions, not iterative series. The following table records the truncation order in the SR3 code, the magnitude of the first omitted term, and the tolerance $\tau_{\mathrm{match}}$ at which the generic evaluator (§§5.2--5.3) produces the same truncation index — thereby reproducing the SR3 result.

| Computation | SR3 truncation | First omitted term | $\tau_{\mathrm{match}}$ |
|-------------|---------------|--------------------|-----------------------|
| $\dot{M}$ secular rate | $O(J_2^2)$, polynomial in $\cos^2 i$ to degree 4 | $O(J_2^3) \approx 10^{-9}$ rad/min | $\sim 10^{-8}$ |
| $\dot{\omega}$ secular rate | $O(J_2^2 + J_4)$, polynomial in $\cos^2 i$ to degree 4 | $O(J_2^3) \approx 10^{-9}$ rad/min | $\sim 10^{-8}$ |
| $\dot{\Omega}$ secular rate | $O(J_2^2 + J_4)$, polynomial in $\cos^2 i$ to degree 2 | $O(J_2^3) \approx 10^{-9}$ rad/min | $\sim 10^{-8}$ |
| Kepler equation | Newton iteration, $\lvert\Delta E\rvert < 10^{-12}$, max 10 iterations, damped at 0.95 | Residual at convergence $< 10^{-12}$ | $10^{-12}$ |
| $q_0$ (ellipsoid) | Not in SR3 (ellipsoid constant, not propagator) | N/A | N/A |
| Short-period $e$ correction | $O(J_2)$, 3 terms in $e$ | $O(J_2^2 \cdot e^2) \approx 10^{-9}$ | $\sim 10^{-8}$ |
| Short-period $i$ correction | $O(J_2)$, 2 terms in $\cos i$ | $O(J_2^2) \approx 10^{-6}$ | $\sim 10^{-5}$ |

At any $\tau \leq \tau_{\mathrm{match}}$, the generic evaluator computes at least as many terms as the SR3 polynomial, and the difference from the SR3 result is below the model accuracy floor $\delta_a(\mathcal{M})$. At $\tau < \tau_{\mathrm{match}}$, additional terms are computed automatically, improving precision beyond the SR3 level.

**Remark.** The SR3 secular rates (lines 1546--1553 of `SGP4.cpp` [Hoots and Roehrich 1980]) are closed-form polynomials in $J_2$, $J_4$, $\cos i$, and $1/p^2$, not iterative series. The generic evaluator reproduces these by expanding the same Brouwer perturbation expressions (Ch 16--19) as a power series in the small parameter $J_2 a_E^2/p^2$ and truncating when the remainder drops below $\tau$. At $\tau = \tau_{\mathrm{match}}$, the evaluator retains exactly the $O(J_2)$ and $O(J_2^2 + J_4)$ terms, matching the SR3 polynomial. At tighter $\tau$, it includes $O(J_2^3)$ and higher terms that the SR3 code omits.

### §5.8.4 Accuracy-Level Stopping

**Proposition 5.8.2** (Accuracy-limited evaluation)**.** *When the tolerance $\tau$ is set to the model accuracy floor $\delta_a(\mathcal{M})$ rather than to machine epsilon, the series evaluators of §§5.2–5.3 require far fewer terms (Table 5.8.2, first data column). The algorithm is unchanged — only the stopping threshold differs. The evaluator still computes the truncation bound at each step and stops when the bound drops below $\tau$. At $\tau = \delta_a(\mathcal{M})$, the truncation bound is dominated by the model error, and computing additional terms yields no operationally meaningful improvement.*

*The precision component $\delta_p$ of the returned TrackedValue records the truncation bound regardless of the tolerance level. A caller inspecting $\delta_p$ can distinguish between a result truncated at $\tau = 10^{-3}$ (where $\delta_p \approx 10^{-3}$) and one truncated at $\tau = 10^{-16}$ (where $\delta_p \approx 10^{-16}$). The error budget remains transparent.*

**Corollary 5.8.1** (Single algorithm, two regimes)**.** *The evaluators of §§5.2–5.4 serve both the matched-pair regime ($\tau = \delta_a(\mathcal{M})$, reproducing the SR3 truncation level) and the precision regime ($\tau = \epsilon_{\mathrm{mach}}$, extracting maximum precision from the arithmetic type). The transition between regimes requires no code change — only the value of $\tau$ passed at initialization (Ch 3, Def. 3.4.1). This is the tolerance parameter principle in action: the algorithm is generic, the stopping condition determines the operating point.*

**Corollary 5.8.2** (Automatic precision scaling)**.** *Changing the arithmetic type $T$ from `double` to a higher-precision type changes $\epsilon_{\mathrm{mach}}$ and therefore the precision-regime $\tau$. The evaluators automatically compute additional terms to match — no code change is required.*

---

## §5.9 Numerical Verification

This section records concrete evaluations of the series in §§5.5–5.6 at WGS84 parameters ($e^2 = 6.6944 \times 10^{-3}$, $e' = 8.2094 \times 10^{-2}$). Reference values are computed at 60-digit precision. The implementation evaluates these series via `alternating_series<T>()` with term functions constructed from `ratio<T>()` and `exact<T>()` (see Appendix C). At `T = double`, the values agree to 15–16 significant digits; at wider `T`, the tables below are reproduced to the full precision of $T$.

### §5.9.1 The $q_0$ Series Term by Term

**Table 5.9.1** (First 10 terms of the $q_0$ series, Eq. 5.11, WGS84)**.**

| $n$ | Term $t_n$ | Partial sum $S_n = 2q_0$ | Leibniz $\lvert t_n \rvert$ | $\lvert t_{n+1}/t_n \rvert$ |
|:--:|---:|---:|:--:|:--:|
| 1 | $1.4754 \times 10^{-4}$ | $1.47540 \times 10^{-4}$ | $1.48 \times 10^{-4}$ | |
| 2 | $-8.5230 \times 10^{-7}$ | $1.46688 \times 10^{-4}$ | $8.52 \times 10^{-7}$ | $0.00578$ |
| 3 | $4.7867 \times 10^{-9}$ | $1.46693 \times 10^{-4}$ | $4.79 \times 10^{-9}$ | $0.00562$ |
| 4 | $-2.7372 \times 10^{-11}$ | $1.46693 \times 10^{-4}$ | $2.74 \times 10^{-11}$ | $0.00572$ |
| 5 | $1.5964 \times 10^{-13}$ | $1.46693 \times 10^{-4}$ | $1.60 \times 10^{-13}$ | $0.00583$ |
| 6 | $-9.4679 \times 10^{-16}$ | $1.46693 \times 10^{-4}$ | $9.47 \times 10^{-16}$ | $0.00593$ |
| 7 | $5.6928 \times 10^{-18}$ | $1.46693 \times 10^{-4}$ | $5.69 \times 10^{-18}$ | $0.00601$ |

Converged value: $q_0 = 7.33462578708345207 \times 10^{-5}$ (25 terms at 60-digit precision).

### §5.9.2 Convergence Rate Verification

**Table 5.9.2** (Actual term ratios vs Eq. 5.19 prediction)**.**

The exact term ratio from Eq. (5.19) is $\frac{(n+1)(2n+1)}{n(2n+5)} \cdot e'^2$, which approaches $e'^2 = 0.006739$ from below.

| $n$ | Actual $\lvert t_{n+1}/t_n \rvert$ | Eq. (5.19) | Asymptote $e'^2$ |
|:--:|:--:|:--:|:--:|
| 1 | $0.005\,7767$ | $0.005\,7767$ | $0.006\,7395$ |
| 2 | $0.005\,6162$ | $0.005\,6162$ | $0.006\,7395$ |
| 3 | $0.005\,7184$ | $0.005\,7184$ | $0.006\,7395$ |
| 4 | $0.005\,8323$ | $0.005\,8323$ | $0.006\,7395$ |
| 5 | $0.005\,9308$ | $0.005\,9308$ | $0.006\,7395$ |
| 6 | $0.006\,0127$ | $0.006\,0127$ | $0.006\,7395$ |
| 7 | $0.006\,0807$ | $0.006\,0807$ | $0.006\,7395$ |
| 8 | $0.006\,1378$ | $0.006\,1378$ | $0.006\,7395$ |
| 9 | $0.006\,1860$ | $0.006\,1860$ | $0.006\,7395$ |

The actual ratios match the Eq. (5.19) prediction to all displayed digits, confirming the derivation. The asymptote $e'^2$ is an upper bound but not tight for small $n$ — the actual ratio at $n = 1$ is 86% of the asymptote.

### §5.9.3 Cancellation in the Closed Form

**Table 5.9.3** (Subtractive cancellation in the closed-form $q_0$, Cor. 5.5.2)**.**

| Quantity | Value |
|----------|-------|
| $(1 + 3/e'^2)\arctan e'$ | $3.654342629738111 \times 10^{1}$ |
| $3/e'$ | $3.654327960486537 \times 10^{1}$ |
| Difference $= 2q_0$ | $1.466925157416690 \times 10^{-4}$ |

The two quantities agree to 5 leading digits. At IEEE 754 double precision, the closed form yields $q_0 = 7.334625786\mathbf{7257} \times 10^{-5}$, compared to the series reference $7.334625787\mathbf{0835} \times 10^{-5}$ — a relative error of $\sim 5 \times 10^{-11}$, leaving $\sim$10 reliable digits. The series form (Table 5.9.1) avoids this cancellation entirely and agrees with the 60-digit reference to the full precision of the arithmetic type.

### §5.9.4 Error Bound Tightness

**Table 5.9.4** (Leibniz bound vs geometric tail bound vs actual remainder, $q_0$ series)**.**

| $N$ | Leibniz $\lvert t_N \rvert$ | Geometric $\lvert t_N \rvert \cdot r/(1-r)$ | Actual $\lvert R_N \rvert$ | Leibniz / actual | Geometric / actual |
|:--:|:--:|:--:|:--:|:--:|:--:|
| 1 | $1.48 \times 10^{-4}$ | $1.00 \times 10^{-6}$ | $8.48 \times 10^{-7}$ | 174 | 1.18 |
| 2 | $8.52 \times 10^{-7}$ | $5.78 \times 10^{-9}$ | $4.76 \times 10^{-9}$ | 179 | 1.22 |
| 3 | $4.79 \times 10^{-9}$ | $3.25 \times 10^{-11}$ | $2.72 \times 10^{-11}$ | 176 | 1.19 |
| 5 | $1.60 \times 10^{-13}$ | $1.08 \times 10^{-15}$ | $9.41 \times 10^{-16}$ | 170 | 1.15 |
| 7 | $5.69 \times 10^{-18}$ | $3.86 \times 10^{-20}$ | $3.44 \times 10^{-20}$ | 166 | 1.12 |
| 10 | $1.31 \times 10^{-24}$ | $8.92 \times 10^{-27}$ | $8.13 \times 10^{-27}$ | 162 | 1.10 |

The Leibniz bound overestimates the actual error by a factor of $\sim$170. The geometric tail bound with $r = e'^2$ overestimates by only 10–20%, making it the tighter bound in practice (consistent with Corollary 5.3.1, since $r \approx 0.007 \ll 1/2$). The reported $\delta_p$ from the alternating series evaluator (§5.2) uses the Leibniz bound and is therefore conservative by two orders of magnitude.

**Remark.** The conservatism of the Leibniz bound is not a defect — it is the price of a bound that requires no knowledge of $r$. When $r$ is known (as it is for all geodetic series), the geometric tail bound should be preferred. For the $q_0$ series, switching from Leibniz to geometric reduces the reported $\delta_p$ by a factor of $\sim$150 with no change in the actual accuracy.

### §5.9.5 TrackedValue Error Budget

**Table 5.9.5** (Precision error decomposition for the $2q_0$ series at double, $\tau = 10^{-16}$)**.**

The table below tracks the error budget of the raw series $2q_0$ (Eq. 5.11) before division by 2. The implementation divides by `exact<T>(2)` to obtain $q_0$, which halves the precision error: $\delta_p(q_0) \approx \delta_p(2q_0)/2$. At convergence ($n = 7$), the implementation reports $\delta_p(q_0) = 2.88 \times 10^{-18}$, consistent with $5.77 \times 10^{-18}/2$.

| $n$ | Arithmetic $\delta_p^{\mathrm{arith}}$ | Truncation bound | Total $\delta_p$ | Dominant source |
|:--:|:--:|:--:|:--:|---|
| 1 | $3.3 \times 10^{-20}$ | $1.48 \times 10^{-4}$ | $1.48 \times 10^{-4}$ | Truncation |
| 3 | $9.8 \times 10^{-20}$ | $4.79 \times 10^{-9}$ | $4.79 \times 10^{-9}$ | Truncation |
| 5 | $1.6 \times 10^{-19}$ | $1.60 \times 10^{-13}$ | $1.60 \times 10^{-13}$ | Truncation |
| 7 | $2.3 \times 10^{-19}$ | $5.69 \times 10^{-18}$ | $5.92 \times 10^{-18}$ | Truncation |
| 8 | $2.6 \times 10^{-19}$ | $3.46 \times 10^{-20}$ | $2.95 \times 10^{-19}$ | Arithmetic |

At $n = 7$ (the stopping point for $\tau = 10^{-16}$), truncation and arithmetic rounding are comparable. By $n = 8$, arithmetic rounding dominates — additional terms would reduce the truncation bound below the arithmetic noise floor. This is the natural stopping point: further terms improve the mathematical truncation but not the computational precision.

**Remark** (Verification protocol)**.** The values in Tables 5.9.1–5.9.5 are computed at 60-digit precision. The implementation reproduces these via `EquipotentialEllipsoid<T>` constructed with WGS84 parameters and varying `series_tolerance`. At `T = double` with `series_tolerance = 1e-16`, the `q0` member agrees with Table 5.9.1 to 15 significant digits. At wider `T`, agreement extends to the full precision of the type. The `TrackedValue<T>` output carries the $\delta_p$ decomposition of Table 5.9.5 in its `errors.precision` field.

---

## Error Notes

**[P.5.1]** Alternating series monotonicity failure. The Leibniz bound (Theorem 5.2.1) requires $|t(k+1).v| \leq |t(k).v|$ for all $k$ in the evaluated range. If the terms are not yet monotonically decreasing at the start index — which can happen for series with binomial coefficients when $|\alpha|$ is large and $k$ is small — the Leibniz bound does not apply.

1. *What is affected:* the truncation bound (5.1) may be invalid for the first few terms.
2. *Quantified bound:* the Leibniz bound holds from the first index $k_0$ where $|t(k_0+1).v| < |t(k_0).v|$. For the binomial series $(1+x)^\alpha$ at $|x| < 1$, this holds for all $k \geq \lceil |\alpha| \rceil$ (from Eq. 5.19 and the analogous bound for arbitrary $\alpha$).
3. *Remedy:* accumulate terms without applying the Leibniz bound until monotonicity is established, then switch to the alternating evaluator. Alternatively, use the geometric-tail evaluator (§5.3) which does not require alternating signs or monotonicity.

**[P.5.2]** Catastrophic cancellation in the closed-form $q_0$. The closed form (5.10) computes $q_0$ as the difference of two quantities that agree in their leading $\sim$5 digits (for WGS84).

1. *What is affected:* the first $\sim$5.4 significant digits of the closed-form result are lost to subtractive cancellation (Ch 1, §1.7).
2. *Quantified bound:* the relative error of the closed form is amplified by a factor $\sim 1/e'^2 \approx 150$ compared to the series form. At double precision, this costs $\log_{10}(150) \approx 2.2$ digits beyond the $\log_{10}(3/e') \approx 3.2$ digits lost in the subtraction, for a total of $\sim$5.4 digits.
3. *Remedy:* always use the series form (5.11), which avoids the cancellation algebraically. The closed form should not be used in any computational context.

**[P.5.3]** Rounding accumulation in long series. The arithmetic rounding in the partial sum grows as $O(N \cdot \epsilon_{\mathrm{mach}} \cdot |S_N|)$ by the triangle inequality over $N$ additions (Ch 1, §1.3).

1. *What is affected:* for the geodetic series with $N \leq 25$, the accumulated rounding is $\leq 25 \epsilon_{\mathrm{mach}} |S_N|$, which is negligible compared to the truncation bound.
2. *Quantified bound:* $\delta_p^{\mathrm{arith}} \leq N \epsilon_{\mathrm{mach}} |S_N|$.
3. *Remedy:* Horner evaluation (§5.4) reduces the accumulation. For very long series ($N > 1000$), compensated summation (Kahan's algorithm) should be used, but this regime does not arise in the orbit propagation series of this chapter.

**[P.5.4]** Geometric tail bound with uncertain convergence ratio. The geometric tail bound (Theorem 5.3.1) requires an upper bound $r$ on $|t(k+1)/t(k)|$ that holds for all $k \geq N$.

1. *What is affected:* if the ratio estimate $r$ is too small, the tail bound underestimates the actual remainder, making the reported $\delta_p$ non-conservative.
2. *Quantified bound:* when the adaptive estimate (5.6) is used, the ratio $r_k$ may fluctuate if the terms have non-monotone coefficient ratios. The bound is rigorous only if $r_k \geq |t(j)/t(j-1)|$ for all $j > k$.
3. *Remedy:* for the series in this chapter, the asymptotic ratio is $e^2$ or $e'^2$, and the coefficient ratio $|c_{k+1}/c_k|$ is eventually monotone. Use the fixed ratio $r = e^2$ (or $r = e'^2$) for the geodetic series, which is rigorous from the first term. For general series, use the maximum of the last several observed ratios as $r$.

**[M.5.1]** Measurement error propagation through series. The measurement errors $\sigma_m$ of the input (e.g., the eccentricity) propagate through the series evaluation via the TrackedValue arithmetic.

1. *What is affected:* the $\sigma_m$ of $e$ propagates through the powers $e^{2k}$ and the multiplications, amplified by the condition number of the series at the evaluation point.
2. *Quantified bound:* for a power series $S(e^2)$ evaluated near $e^2 = 0.0067$, the condition number is approximately 1 (dominated by the constant term), so $\sigma_m(S) \approx |S'(e^2)| \cdot \sigma_m(e^2) \approx \sigma_m(e^2)$. The amplification is negligible.
3. *Remedy:* handled automatically by TrackedValue arithmetic. No special treatment is needed.

**[A.5.1]** Ephemeris truncation as accuracy error. The solar and lunar ephemerides (Ch 25–26) are truncated trigonometric series where each omitted term represents a physical perturbation.

1. *What is affected:* the ephemeris position is shifted by the unmodeled terms, introducing an accuracy error $\delta_a$ that cannot be reduced by increasing the arithmetic precision.
2. *Quantified bound:* the SGP4 solar ephemeris retains terms with amplitudes $\geq 0.01°$ in ecliptic longitude, giving $\delta_a \approx 0.01°$ ($\sim$6 arcsec). The lunar ephemeris has $\delta_a \approx 0.1°$ ($\sim$6 arcmin). These are accuracy errors of the ephemeris model.
3. *Remedy:* use a higher-fidelity ephemeris (DE440, VSOP87) when the application requires it. The SGP4 ephemeris is a matched-pair value [MP]: replacing it with a better ephemeris without re-fitting the TLE element set may degrade overall accuracy (Ch 3, Matched Pair Principle).
