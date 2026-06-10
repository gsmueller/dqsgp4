# Fast Convergence Methods for Arbitrary Precision Function Evaluation

## Why Not Taylor Series?

A Taylor series for $\sin(x)$ centered at 0 converges at rate $O(|x|^n / n!)$. For $x = 1$ radian and 50 digits of precision, this requires ~40 terms. Each term involves a multiprecision multiplication costing $O(p \log p)$ where $p$ is the precision in bits. That's 40 expensive multiplications.

But Newton's method for computing $\sqrt{x}$ doubles the number of correct digits per iteration — 50 digits needs only ~8 iterations. That's 5x fewer expensive operations for the same result.

The question for every function we need is: **what algorithm achieves the target precision in the fewest multiprecision operations?**

## 1. Quadratically Convergent Methods (Digits Double Per Iteration)

### 1.1 Newton's Method for Inverse Functions

If we can evaluate $f(x)$, we can compute $f^{-1}(y)$ via:

$$x_{n+1} = x_n - \frac{f(x_n) - y}{f'(x_n)}$$

Each iteration doubles the number of correct digits.

| Function | Computed Via Newton On | Iterations for 50 digits |
|---|---|---|
| $\sqrt{x}$ | $f(t) = t^2 - x$ | ~8 |
| $1/x$ | $f(t) = 1/t - x$ (no division needed!) | ~8 |
| $\arcsin(x)$ | $f(t) = \sin(t) - x$ | ~8 (if sin is fast) |
| $\arctan(x)$ | $f(t) = \tan(t) - x$ | ~8 (if tan is fast) |
| $\exp(x)$ | Via $\log$ inversion | ~8 (if log is fast) |

**Key insight:** Newton's method shifts the cost to evaluating the forward function. If the forward function is cheaper than the inverse, Newton wins.

### 1.2 Arithmetic-Geometric Mean (AGM)

The AGM of $a_0$ and $b_0$ is the common limit of:

$$a_{n+1} = \frac{a_n + b_n}{2}, \qquad b_{n+1} = \sqrt{a_n b_n}$$

This converges quadratically. The AGM connects to:

$$\log(x) = \frac{\pi}{2 \cdot \text{AGM}(1, 4/s)} - m\log 2$$

where $s = x \cdot 2^m$ is scaled so $s > 2^{p/2}$ (p = precision bits).

**Through $\log$, AGM gives fast access to $\exp$, and through those to all trig functions.**

| Function | AGM Path | Convergence |
|---|---|---|
| $\log(x)$ | Direct AGM formula | Quadratic |
| $\exp(x)$ | Newton inversion of $\log$ | Quadratic |
| $\pi$ | AGM(1, $1/\sqrt{2}$) via Gauss-Legendre | Quadratic |
| $\sin(x), \cos(x)$ | Via $\exp(ix)$ or argument reduction + AGM | Quadratic |
| $\arctan(x)$ | Via $\log$ of complex argument | Quadratic |

**Cost per iteration:** One AGM iteration requires one addition, one multiplication, and one square root — all $O(p \log p)$ for $p$-bit precision. Total for $n$ digits: $O(\log n)$ iterations × $O(n \log n)$ per iteration = $O(n \log^2 n)$.

### 1.3 Halley's Method (Cubic Convergence)

$$x_{n+1} = x_n - \frac{2f(x_n)f'(x_n)}{2f'(x_n)^2 - f(x_n)f''(x_n)}$$

Triples digits per iteration. Useful when $f''$ is cheap (e.g., for $\sin$, $\cos$).

## 2. Superlinear Series Methods

### 2.1 Binary Splitting

For a hypergeometric-like series $\sum_{k=0}^{N} a_k / b_k$ where $a_k$ and $b_k$ are integers or rationals, **binary splitting** evaluates the sum in $O(M(N\log N) \log N)$ bit operations instead of $O(N \cdot M(N))$ for naive term-by-term summation.

The idea: split the sum into two halves, compute each recursively, then combine. The key gain is that intermediate results are computed with exact integer arithmetic (no precision loss), and the final division happens only once at the end.

**Applicable to:** Taylor series for $\exp$, $\sin$, $\cos$, $\arctan$, $\log(1+x)$, and our geodetic series (which are rational-coefficient power series).

### 2.2 Machin-Type Formulas for $\arctan$

Instead of computing $\arctan(e')$ directly (slow Taylor convergence for $e' \approx 0.082$), decompose into sums of $\arctan$ of smaller arguments:

$$\arctan(x) = \sum_i c_i \arctan(a_i/b_i)$$

where each $a_i/b_i$ is small (so the Taylor series converges fast) and the $c_i$ are small integers.

Example: **Machin's formula** for $\pi$:

$$\frac{\pi}{4} = 4\arctan\frac{1}{5} - \arctan\frac{1}{239}$$

For $\arctan(e')$ where $e' \approx 0.082$, the argument is already small enough that the Taylor series converges reasonably. But we could still use:

$$\arctan(e') = \arctan\frac{p}{q} + \arctan\frac{e' - p/q}{1 + e' \cdot p/q}$$

choosing $p/q$ as a rational approximation to $e'$, making the second $\arctan$ argument very small.

### 2.3 Continued Fractions

For many functions, continued fractions converge faster than Taylor series:

$$\arctan(x) = \cfrac{x}{1 + \cfrac{x^2}{3 + \cfrac{4x^2}{5 + \cfrac{9x^2}{7 + \cdots}}}}$$

General pattern: $\arctan(x) = \cfrac{x}{1+} \cfrac{(1x)^2}{3+} \cfrac{(2x)^2}{5+} \cfrac{(3x)^2}{7+} \cdots$

This converges for all real $x$ (unlike the Taylor series which only converges for $|x| \leq 1$) and generally converges faster than the Taylor series within $|x| \leq 1$.

Continued fractions can be evaluated efficiently using the **modified Lentz algorithm** or **Wallis recurrence**.

### 2.4 Padé Approximants

Rational function approximations $[L/M]$ that match the Taylor series through order $L+M$:

$$f(x) \approx \frac{P_L(x)}{Q_M(x)}$$

For the same total degree, Padé approximants converge faster than Taylor polynomials, especially near singularities. They're particularly good for functions with poles (like $\tan$).

## 3. Argument Reduction

Before evaluating any trig function, reduce the argument to a small range where the series converges fast:

### 3.1 Standard Trig Reduction

$$\sin(x) = \sin(x \bmod 2\pi)$$

Then reduce further to $[0, \pi/4]$ using symmetries:

$$\sin(x) = \cos(\pi/2 - x), \quad \sin(\pi - x) = \sin(x), \quad \text{etc.}$$

**Critical issue for arbitrary precision:** The modular reduction $x \bmod 2\pi$ loses precision when $x$ is large. The **Payne-Hanek algorithm** handles this by using a high-precision representation of $2/\pi$ to extract the correct reduced argument.

### 3.2 Repeated Halving for $\exp$

$$\exp(x) = \exp(x/2^k)^{2^k}$$

Choose $k$ so $|x/2^k| < 2^{-p/4}$, then the Taylor series for $\exp(x/2^k)$ converges in very few terms. Then square the result $k$ times.

### 3.3 Repeated Halving for $\sin$/$\cos$

$$\sin(x) = 2\sin(x/2)\cos(x/2)$$

$$\cos(x) = 2\cos^2(x/2) - 1$$

Halve the argument $k$ times until it's small, compute $\sin$ and $\cos$ of the tiny argument with a short Taylor/Padé series, then reconstruct using the double-angle formulas.

## 4. Specific Recommendations for Our Functions

### 4.1 $\arctan(e')$ in $q_0$ and $q_0'$

$e' \approx 0.082$ — small argument, Taylor converges OK but not great.

**Recommended:** Continued fraction evaluation (converges for all real $x$, faster than Taylor for this range). Or: use the series form of $q_0$ directly (which avoids computing $\arctan$ separately):

$$2q_0 = \sum_{n=1}^{\infty} \frac{4(-1)^{n+1}n}{(2n+1)(2n+3)} e'^{2n+1}$$

This alternating series converges rapidly for $e' \approx 0.082$ (each term is smaller by factor $e'^2 \approx 0.0067$) and avoids the subtractive cancellation problem of the closed form. **This may be the best approach for our case** — we already identified this in our series analysis.

### 4.2 $\sin(\varphi)$, $\cos(\varphi)$ in Somigliana and SGP4

Geodetic latitudes are bounded: $|\varphi| \leq \pi/2$.

**Recommended:** Argument reduction to $[0, \pi/4]$ + Padé approximant or short Taylor with binary splitting. For the SGP4 inner loop where sin/cos are called many times, precompute once.

Boost.Multiprecision already uses optimized algorithms internally — we should verify what it uses and whether it's sufficient, rather than reimplementing.

### 4.3 Kepler's Equation $E - e\sin E = M$

Currently solved by Newton-Raphson (quadratic convergence). Could use:

- **Halley's method** (cubic convergence) since $f'' = e\sin E$ is cheap:

$$E_{n+1} = E_n - \frac{f \cdot f'}{f'^2 - \frac{1}{2}f \cdot f''}$$

- **Householder's method** (quartic convergence) for even fewer iterations.

For 50-digit precision: Newton needs ~8 iterations, Halley needs ~6, Householder needs ~5. Each iteration costs a sin + cos evaluation.

### 4.4 Geodetic Series ($R_2$, $Q$, $\bar{\gamma}$)

These are polynomial series in $e'^2$ (or $e^2$) with rational coefficients. Since $e'^2 \approx 0.0067$, convergence is rapid — each term is ~150x smaller than the previous.

**Recommended:** Direct evaluation with our generalized binomial framework, using exact rational arithmetic for the coefficients. Binary splitting if many terms are needed (unlikely — 8-10 terms give 50+ digits for these series).

The series approach IS the right method here — the question is just efficient coefficient generation and error bounding, which we've already designed.

### 4.5 $\pi$ to Arbitrary Precision

**Recommended:** Use `boost::math::constants::pi<T>()` which uses the optimal algorithm for the precision of `T`. For `cpp_bin_float_50`, Boost uses Machin-type or AGM internally.

We should NOT reimpute $\pi$ — use the library.

## 5. What Boost.Multiprecision Already Provides

Before implementing custom algorithms, verify what Boost already does:

| Function | Boost Algorithm | Sufficient? |
|---|---|---|
| `sin`, `cos` | Argument reduction + Taylor | Likely yes for moderate precision |
| `sqrt` | Newton's method | Yes — quadratic convergence |
| `exp`, `log` | AGM-based or Taylor with argument reduction | Likely yes |
| `atan`, `atan2` | Taylor series | **May be slow** — consider continued fraction |
| `pow(x, rational)` | Via exp/log | **Slow for multiprecision** — use explicit forms |

**Key concern:** Boost's `atan` implementation may use a naive Taylor series. For our $q_0$ computation that calls $\arctan(e')$, we should benchmark Boost's `atan` against our series form of $q_0$ and against a continued fraction implementation.

## 6. Clenshaw Summation

For evaluating sums of the form $\sum c_k P_k(x)$ where $P_k$ are orthogonal polynomials (Chebyshev, Legendre), **Clenshaw's algorithm** is numerically stable and efficient:

$$y_k = 2x \cdot y_{k+1} - y_{k+2} + c_k$$

evaluated backwards from $k = N$ to $k = 0$. This is relevant for:
- Evaluating spherical harmonic series (EGM2008)
- Evaluating Legendre polynomial sums in the gravitational potential
- Chebyshev-based function approximations

## 7. Summary Table

| Function | Naive Taylor | Best Method | Speedup |
|---|---|---|---|
| $\sqrt{x}$ | Not applicable | Newton iteration | N/A |
| $\sin(x)$, $\cos(x)$ | $O(n)$ terms | Arg reduction + short series | ~4x fewer terms |
| $\arctan(x)$ | $O(n)$ terms, $|x| \leq 1$ only | Continued fraction or Machin | ~2x faster, works for all $x$ |
| $\exp(x)$ | $O(n)$ terms | Arg reduction + binary splitting | ~$O(\sqrt{n})$ effective terms |
| $\log(x)$ | $O(n)$ terms, $|x-1| < 1$ only | AGM | $O(\log n)$ iterations |
| $\pi$ | Slow | Chudnovsky (~14 digits/term) | ~100x |
| Kepler's equation | Newton (quadratic) | Halley (cubic) or Householder (quartic) | ~30% fewer iterations |
| Geodetic series | Direct summation | Binary splitting + exact rationals | Marginal for ~10 terms |
| $q_0$, $q_0'$ | Closed form (cancellation!) | Series form (no cancellation) | Avoids precision loss entirely |

## 8. Design Impact

### The `math/` Library Should Provide:

```
math/
    series/
        binary_splitting.h      — Fast evaluation of rational series
        clenshaw.h              — Clenshaw summation for orthogonal polynomial series
        continued_fraction.h    — Continued fraction evaluation (Lentz/Wallis)
        pade.h                  — Padé approximant construction and evaluation
    functions/
        arctan_cf.h             — Continued fraction arctan (faster than Taylor)
        kepler_solver.h         — Halley/Householder Kepler equation solver
        argument_reduction.h    — Payne-Hanek and standard trig reduction
    core/
        tracked_value.h         — Error-tracked arithmetic
        binomial.h              — Generalized binomial coefficients
        wallis.h                — Wallis integrals (exact rational)
        factorial.h             — Falling/rising factorial, double factorial
```

### Key Principle: Prefer Methods That Match the Precision Level

| Precision | Strategy |
|---|---|
| `double` (~16 digits) | Boost built-in functions are fine; Horner evaluation of short polynomials |
| 50 digits | Argument reduction + moderate series; Newton/Halley for inverses; continued fractions for arctan |
| 100+ digits | AGM for log/exp; binary splitting for series; Machin-type for arctan if needed |

The `TrackedValue<T>` precision bound tells us when we've converged — we don't need to know the precision level in advance. The algorithm choice can be made at compile time based on the type `T`, or at runtime based on the requested tolerance.
