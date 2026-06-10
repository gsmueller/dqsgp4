# Equation-Specific Convergence Tricks

For each equation we will implement, this document identifies the mathematical structure that enables faster convergence than naive series evaluation.

## 1. $q_0$ and $q_0'$ — The Most Precision-Critical Computation

### The Problem

Closed form: $2q_0 = \left(1 + \frac{3}{e'^2}\right)\arctan(e') - \frac{3}{e'}$

For $e' \approx 0.082$, the term $(1 + 3/e'^2) \approx 447$ multiplies $\arctan(e') \approx 0.0819$, giving $\approx 36.6$, then subtracting $3/e' \approx 36.6$. The result $2q_0 \approx 0.000147$ is obtained by subtracting two large numbers — **catastrophic cancellation**.

### Solution: Direct Series (Avoids Cancellation Entirely)

$$2q_0 = \sum_{n=1}^{\infty} \frac{4(-1)^{n+1} n}{(2n+1)(2n+3)} e'^{2n+1}$$

This is an alternating series with ratio $\sim e'^2 \approx 0.0067$. Each term is ~150× smaller than the previous. For 50 digits: ~8 terms. For 100 digits: ~16 terms.

**Convergence trick: Exact rational coefficients.**

The coefficient $\frac{4n}{(2n+1)(2n+3)}$ is an exact rational number. Compute it in exact rational arithmetic (integer numerator/denominator), then multiply by $e'^{2n+1}$. This eliminates ALL rounding in the coefficients — the only rounding is in the power of $e'$.

**Further acceleration: Horner-like nesting.**

Factor $e'^3$ out: $2q_0 = e'^3 \sum_{n=1}^{\infty} \frac{4n}{(2n+1)(2n+3)} (-e'^2)^{n-1}$

Then evaluate the inner sum via Horner's method in $(-e'^2)$:

$$S = c_N + (-e'^2)(c_{N-1} + (-e'^2)(c_{N-2} + \cdots))$$

where $c_n = \frac{4n}{(2n+1)(2n+3)}$. This minimizes the number of multiplications.

### $q_0'$ — Same Approach

$$q_0' = 3\left[\left(1 + \frac{1}{e'^2}\right)\left(1 - \frac{1}{e'}\arctan(e')\right)\right] - 1$$

**Series derivation for $q_0'$:**

Starting from $\arctan(e') = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n+1}$:

$$\frac{1}{e'}\arctan(e') = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n} = 1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \cdots$$

$$1 - \frac{1}{e'}\arctan(e') = \frac{e'^2}{3} - \frac{e'^4}{5} + \frac{e'^6}{7} - \cdots = \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n}$$

Now $(1 + 1/e'^2)$ times this sum:

$$\left(1 + \frac{1}{e'^2}\right)\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} = \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} + \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n-2}$$

Reindex the second sum ($n \to n+1$):

$$= \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} + \sum_{n=0}^{\infty} \frac{(-1)^{n+2}}{2n+3} e'^{2n} = \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} + \frac{1}{3} + \sum_{n=1}^{\infty} \frac{(-1)^{n}}{2n+3} e'^{2n}$$

So $q_0' = 3[\frac{1}{3} + \sum_{n=1}^{\infty}(-1)^{n+1}(\frac{1}{2n+1} - \frac{1}{2n+3})e'^{2n}] - 1$

$= 1 + 3\sum_{n=1}^{\infty}(-1)^{n+1}\frac{2}{(2n+1)(2n+3)}e'^{2n} - 1$

$$q_0' = \sum_{n=1}^{\infty} \frac{6(-1)^{n+1}}{(2n+1)(2n+3)} e'^{2n}$$

**Verification:** For $n=1$: $\frac{6}{3 \times 5} e'^2 = \frac{2}{5}e'^2$. From closed form with small $e'$: $q_0' \approx \frac{2}{5}e'^2 + O(e'^4)$. ✓

This is an alternating series with ratio $\sim e'^2 \approx 0.0067$. Converges rapidly. Evaluate via Horner in $(-e'^2)$ with exact rational coefficients $\frac{6}{(2n+1)(2n+3)}$.

**Error bound:** First omitted term. Automatic via `BoundedResult`.

## 2. $R_2$ — Sphere of Equal Area

### Structure

$$R_2 = c\sqrt{\int_0^{\pi/2} \frac{\cos\Phi}{(1+e'^2\cos^2\Phi)^2} d\Phi}$$

### Three-Stage Computation

**Stage 1: Binomial expansion of integrand**

$(1+e'^2\cos^2\Phi)^{-2} = \sum_{k=0}^{N} (-1)^k(k+1) e'^{2k}\cos^{2k}\Phi$

Coefficient $(-1)^k(k+1)$ is an exact integer — zero rounding.

**Stage 2: Wallis integration (exact rational)**

$\int_0^{\pi/2}\cos^{2k+1}\Phi\,d\Phi = \frac{(2k)!!}{(2k+1)!!} = \frac{4^k(k!)^2}{(2k+1)!}$

This is an exact rational number. Compute as integer ratio.

**Stage 3: Square root via Newton iteration**

Rather than expanding $\sqrt{I}$ as another binomial series (which doubles the number of terms), use **Newton's method** for the square root:

$$x_{n+1} = \frac{1}{2}\left(x_n + \frac{I}{x_n}\right)$$

Starting from $x_0 = 1$ (since $I \approx 1$), this converges quadratically. For 50 digits: ~8 iterations. Each iteration requires one division of `TrackedValue<T>`, which is cheaper than evaluating a second binomial series.

**Convergence trick:** The integral series itself converges at rate $e'^2 \approx 0.0067$ per term. The square root Newton iteration converges quadratically. Total cost: ~8 series terms + ~8 Newton iterations = ~16 multiprecision operations.

## 3. Mean Gravity $\bar{\gamma}$ — Ratio of Two Fast-Converging Series

### Structure

$$\bar{\gamma}/\gamma_e = \frac{\sum_{k=0}^{N} a_k e^{2k}}{\sum_{k=0}^{N} b_k e^{2k}}$$

where $a_k$ involves both $e^{2k}$ and $k$ terms, and $b_k = (k+1)/(2k+1)$.

### Convergence Trick: Series Division via Recurrence

Rather than computing numerator and denominator separately and dividing, use the recurrence for the quotient of two power series directly:

If $N(x) = \sum n_k x^k$ and $D(x) = \sum d_k x^k$, then $Q(x) = N(x)/D(x) = \sum q_k x^k$ where:

$$q_0 = n_0/d_0, \quad q_k = \frac{1}{d_0}\left(n_k - \sum_{j=1}^{k} d_j q_{k-j}\right)$$

This computes the quotient coefficients incrementally, stopping when $|q_k e^{2k}|$ is below tolerance. No separate division at the end.

**Error bound:** Since both series converge at rate $e^2 \approx 0.0067$, the quotient also converges at this rate. The error from truncating at term $N$ is bounded by $|q_N e^{2N}| / (1 - e^2)$.

## 4. Somigliana Normal Gravity $\gamma(\varphi)$ — Closed Form, No Series Needed

### Structure

$$\gamma = \gamma_e \frac{1 + k\sin^2\varphi}{\sqrt{1 - e^2\sin^2\varphi}}$$

This is a closed-form expression. The only transcendental operation is the square root, which converges quadratically via Newton. The $\sin^2\varphi$ is computed once.

**No series expansion needed for evaluation.** The series form (gravity formula 1980) is used only for validation, not computation.

### Convergence Trick for the Reciprocal Square Root

Instead of computing $\sqrt{1 - e^2\sin^2\varphi}$ and dividing, compute the **reciprocal square root** directly:

$$\frac{1}{\sqrt{1 - e^2\sin^2\varphi}}$$

via Newton iteration for $f(x) = 1/x^2 - (1-e^2\sin^2\varphi)$:

$$x_{n+1} = \frac{x_n}{2}(3 - (1-e^2\sin^2\varphi)x_n^2)$$

This avoids division entirely (only multiplications), which is faster for multiprecision types. Starting from $x_0 = 1$, converges quadratically.

## 5. Kepler's Equation $E - e\sin E = M$ — Higher-Order Solvers

### Current: Newton-Raphson (Quadratic)

$$E_{n+1} = E_n - \frac{E_n - e\sin E_n - M}{1 - e\cos E_n}$$

### Better: Halley's Method (Cubic Convergence)

$$f = E - e\sin E - M, \quad f' = 1 - e\cos E, \quad f'' = e\sin E$$

$$E_{n+1} = E_n - \frac{2ff'}{2f'^2 - ff''}$$

Triples correct digits per iteration. For 50 digits: ~6 iterations instead of ~8.

### Even Better: Householder's Method (Quartic Convergence)

$$E_{n+1} = E_n - \frac{f}{f'}\left(1 + \frac{ff''}{2f'^2} + \frac{f^2(f'^2 f'' - 3f'f''^2 + 3f''^3)}{6f'^4}\right)^{-1}$$

Wait, that's getting complicated. For Kepler's equation, since $f''' = -e\cos E = -(f' - 1)$, all derivatives are trivially cheap. So:

$$f''' = -(1 - e\cos E) + 1 = -f' + 1$$

The Householder iteration of order 4:

$$E_{n+1} = E_n - \frac{f}{f'}\cdot\frac{1 - \frac{ff''}{2f'^2}}{1 - \frac{ff''}{f'^2} + \frac{f^2 f'''}{6f'^3}}$$

Each iteration costs: 1 sin, 1 cos, and a few multiplications/divisions. Quartic convergence: ~5 iterations for 50 digits.

### Convergence Trick: Starter Value

A good starting value reduces iterations by 1-2. For near-circular orbits ($e < 0.3$):

$$E_0 = M + e\sin M + \frac{e^2}{2}\sin 2M$$

This is the first three terms of the Fourier-Bessel expansion of $E(M, e)$. It gives ~3-4 correct digits, saving one iteration of any method.

For higher eccentricity, use the Markley (1995) starter or the Mikkola (1987) cubic approximation.

## 6. Brouwer Secular Rates — Exact Polynomial Evaluation

### Structure

$$\dot{M} = n + \frac{3}{2}nJ_2\frac{1}{p^2}\beta(3\cos^2 i - 1) + \frac{1}{16}nJ_2^2\frac{1}{p^4}\beta(13 - 78\cos^2 i + 137\cos^4 i) + \cdots$$

These are **finite polynomials** in $\cos^2 i$, not infinite series. No convergence issue — just exact evaluation.

### Convergence Trick: Horner's Method in $\theta^2 = \cos^2 i$

$$13 - 78\theta^2 + 137\theta^4 = 13 + \theta^2(-78 + 137\theta^2)$$

Evaluates in 2 multiplications + 2 additions instead of the naive 4 multiplications + 2 additions. For $n$ terms, Horner reduces from $O(n^2)$ to $O(n)$ multiplications.

### Note on Extending the Series

If we want higher-order Brouwer terms ($J_2^3$, $J_2 J_3$, etc.), these are also finite polynomials in $\cos^2 i$ with rational coefficients. The coefficients are derivable from the Brouwer (1959) Hamiltonian perturbation theory. Each additional order adds one power of $J_2/p^2$ and extends the $\cos^2 i$ polynomial by two degrees.

## 7. Even Zonal Harmonics $J_{2n}$ — Closed Formula

$$J_{2n} = (-1)^{n+1}\frac{3e^{2n}}{(2n+1)(2n+3)}\left(1 - n + 5n\frac{J_2}{e^2}\right)$$

This is a closed-form expression for each $n$. No series convergence issue. Exact rational prefactor $\frac{3}{(2n+1)(2n+3)}$ computed in integer arithmetic.

## 8. Deep Space Integration Step — Adaptive Step Size

### Current: Fixed 720-minute step with 2nd-order integrator

$$x_{n+1} = x_n + \dot{x}\Delta t + \frac{1}{2}\ddot{x}\Delta t^2$$

### Convergence Trick: Precision-Dependent Step Size

The local truncation error of this 2nd-order integrator is $O(\Delta t^3)$. For the standard 720-minute step:

$$\epsilon_\text{step} \sim \frac{1}{6}\dddot{x}\Delta t^3$$

For higher precision, we should either:
1. **Reduce $\Delta t$** — halving the step quarters the error but doubles the number of steps
2. **Use a higher-order integrator** — 4th-order Runge-Kutta gives $O(\Delta t^5)$ error for the same step

However, since this is the standard SGP4 integrator, changing the order or step size would make our results incompatible with standard SGP4 validation. **For SGP4 compatibility, keep the 720-minute step and 2nd-order integrator.** Track the integration error via `TrackedValue` precision bounds — the caller can see when integration error dominates.

For a future "enhanced" mode beyond standard SGP4, implement adaptive step size where the step is halved whenever the `TrackedValue` precision bound exceeds the target tolerance.

## 9. $\arctan$ via Continued Fraction (for cases where we must compute it)

If for any reason we need $\arctan(x)$ directly (rather than avoiding it via series):

### Continued Fraction (Euler, 1748)

$$\arctan(x) = \cfrac{x}{1 + \cfrac{x^2}{3 + \cfrac{4x^2}{5 + \cfrac{9x^2}{7 + \cfrac{16x^2}{9 + \cdots}}}}}$$

General form: $\arctan(x) = \cfrac{x}{1+}\cfrac{(1\cdot x)^2}{3+}\cfrac{(2\cdot x)^2}{5+}\cfrac{(3\cdot x)^2}{7+}\cdots$

For $x = e' \approx 0.082$, this converges faster than the Taylor series because the "effective ratio" of successive convergents is $\sim x^2/(4n^2)$ rather than $\sim x^2$ for Taylor.

### Evaluation: Modified Lentz Algorithm

The Lentz algorithm evaluates the continued fraction in the forward direction with $O(1)$ storage per step, and naturally provides an error estimate from the ratio of successive convergents.

```
f = tiny (e.g., 1e-300)
C = f
D = 0
for n = 0, 1, 2, ...
    a_n, b_n = CF coefficients
    D = b_n + a_n * D
    C = b_n + a_n / C
    D = 1/D
    delta = C * D
    f = f * delta
    if |delta - 1| < tolerance: break
```

**Error bound:** $|f - f_\text{true}| \lesssim |f| \cdot |delta - 1|$ at convergence.

## 10. Normal Potential $U_0$ — Mixed Closed/Series

$$U_0 = \frac{GM}{E}\arctan e' + \frac{\omega^2 a^2}{3}$$

### Convergence Trick: Use Series for arctan(e') Inside U₀

$$\frac{1}{E}\arctan e' = \frac{1}{ae}\cdot\arctan\frac{e}{\sqrt{1-e^2}} = \frac{1}{ae}\sum_{n=0}^{\infty}\frac{(-1)^n}{2n+1}\frac{e^{2n+1}}{(1-e^2)^{n+1/2}}$$

Actually, it's simpler to write:

$$\frac{\arctan e'}{E} = \frac{1}{b}\sum_{n=0}^{\infty}\frac{(-1)^n}{2n+1}e'^{2n} = \frac{1}{b}\left(1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \cdots\right)$$

So:

$$U_0 = \frac{GM}{b}\left(1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \frac{e'^6}{7} + \cdots\right) + \frac{\omega^2 a^2}{3}$$

This is the series form from [M80] Sec. 3. Alternating, rapidly converging (ratio $e'^2 \approx 0.0067$), exact rational coefficients $1/(2n+1)$.

**No cancellation in this form** — the $GM/b$ and $\omega^2 a^2/3$ terms are both positive and of different magnitude, so they add without cancellation.

## Summary: Computation Cost per Equation

| Equation | Method | Operations for 50 digits | Operations for 100 digits |
|---|---|---|---|
| $q_0$ | Horner series (exact rational coefficients) | ~8 mult | ~16 mult |
| $q_0'$ | Horner series | ~8 mult | ~16 mult |
| $R_2$ | Binomial series + Newton sqrt | ~8 + ~8 = 16 | ~16 + ~9 = 25 |
| $\bar{\gamma}$ | Series quotient recurrence | ~8 | ~16 |
| $Q$ | Horner series | ~8 | ~16 |
| $\gamma(\varphi)$ | Closed form + reciprocal sqrt Newton | ~8 | ~9 |
| Kepler | Halley iteration + good starter | ~6 | ~7 |
| $U_0$ | Alternating series | ~8 | ~16 |
| $J_{2n}$ | Closed form (single evaluation) | ~3 | ~3 |
| $\arctan(x)$ | Continued fraction (Lentz) | ~12 | ~20 |
| Brouwer polynomials | Horner evaluation | ~5 | ~5 (exact polynomial) |

All series converge at rate $e'^2 \approx 0.0067$ (~2.2 digits per term), so the number of terms scales linearly with the desired digits. The Newton/Halley iterations double/triple digits per step, so they scale logarithmically.
