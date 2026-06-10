# Derivation of Series Coefficients from First Principles

## Overview

Nearly all series coefficients in the geodetic formulas arise from the same three-step pattern:

1. **Generalized binomial expansion** of $(1 \pm \epsilon\,\text{trig}^2\Phi)^{\alpha}$ for non-integer $\alpha$
2. **Wallis-type cosine/sine integrals** to evaluate each term
3. **Series algebra** (square root, reciprocal, quotient) to combine results

This document traces how each published coefficient is generated, identifying the underlying mathematical structures that our `math/` library must implement.

---

## 1. The Generalized Binomial Theorem

### 1.1 Falling Factorial

The falling factorial of $\alpha$ to $k$ terms is:

$$\alpha^{(k)} = \alpha(\alpha - 1)(\alpha - 2) \cdots (\alpha - k + 1) = \prod_{j=0}^{k-1}(\alpha - j)$$

with $\alpha^{(0)} = 1$.

For positive integer $\alpha = n$, this reduces to $n!/(n-k)!$ and terminates at $k = n$. For non-integer or negative $\alpha$, it produces an infinite sequence.

### 1.2 Rising Factorial (Pochhammer Symbol)

The rising factorial is:

$$(\alpha)_k = \alpha(\alpha + 1)(\alpha + 2) \cdots (\alpha + k - 1) = \prod_{j=0}^{k-1}(\alpha + j)$$

with $(\alpha)_0 = 1$.

**Connection:** $(\alpha)_k = (-1)^k (-\alpha)^{(k)}$

### 1.3 Generalized Binomial Coefficient

$$\binom{\alpha}{k} = \frac{\alpha^{(k)}}{k!} = \frac{\alpha(\alpha-1)(\alpha-2)\cdots(\alpha-k+1)}{k!}$$

For non-integer $\alpha$, the generalized binomial theorem gives:

$$(1 + x)^{\alpha} = \sum_{k=0}^{\infty} \binom{\alpha}{k} x^k, \qquad |x| < 1$$

### 1.4 Expansion of $(1 - x)^{-\alpha}$ (the form needed for geodesy)

Using the rising factorial:

$$(1 - x)^{-\alpha} = \sum_{k=0}^{\infty} \frac{(\alpha)_k}{k!} x^k$$

**Proof:** $\binom{-\alpha}{k}(-x)^k = (-1)^k \frac{(-\alpha)^{(k)}}{k!}(-x)^k = \frac{(\alpha)_k}{k!}x^k$

### 1.5 Specific Cases Used in Geodesy

| Exponent $\alpha$ | Application | Binomial coefficients $\binom{\alpha}{k}$ for $k = 0, 1, 2, 3, 4$ |
|---|---|---|
| $-1/2$ | $N = a(1 - e^2\sin^2\Phi)^{-1/2}$ | $1, \; -\frac{1}{2}, \; \frac{3}{8}, \; -\frac{5}{16}, \; \frac{35}{128}$ |
| $-3/2$ | Meridian quadrant $Q$ integrand | $1, \; -\frac{3}{2}, \; \frac{15}{8}, \; -\frac{35}{16}, \; \frac{315}{128}$ |
| $-2$ | $R_2$ integrand, surface area | $1, \; -2, \; 3, \; -4, \; 5$ |
| $-5/2$ | Mean gravity numerator | $1, \; -\frac{5}{2}, \; \frac{35}{8}, \; -\frac{105}{16}, \; \frac{1155}{128}$ |
| $+1/2$ | Square root of a series | $1, \; \frac{1}{2}, \; -\frac{1}{8}, \; \frac{1}{16}, \; -\frac{5}{128}$ |

Note: For $\alpha = -2$, the binomial coefficients simplify to $\binom{-2}{k} = (-1)^k(k+1)$.

---

## 2. Wallis Integrals

### 2.1 Definition

$$W_n = \int_0^{\pi/2} \cos^n\Phi \; d\Phi = \int_0^{\pi/2} \sin^n\Phi \; d\Phi$$

### 2.2 Recurrence Relation

$$W_n = \frac{n-1}{n} W_{n-2}, \qquad W_0 = \frac{\pi}{2}, \quad W_1 = 1$$

### 2.3 Closed Forms

**Even power** ($n = 2k$):

$$W_{2k} = \frac{(2k-1)!!}{(2k)!!} \cdot \frac{\pi}{2} = \frac{1 \cdot 3 \cdot 5 \cdots (2k-1)}{2 \cdot 4 \cdot 6 \cdots (2k)} \cdot \frac{\pi}{2}$$

**Odd power** ($n = 2k+1$):

$$W_{2k+1} = \frac{(2k)!!}{(2k+1)!!} = \frac{2 \cdot 4 \cdot 6 \cdots (2k)}{1 \cdot 3 \cdot 5 \cdots (2k+1)}$$

### 2.4 Key Special Integral

The integral $\int_0^{\pi/2} \sin^{2k}\Phi \cos\Phi \; d\Phi$ is NOT a Wallis integral — it's simpler:

$$\int_0^{\pi/2} \sin^{2k}\Phi \cos\Phi \; d\Phi = \frac{1}{2k+1}$$

(Substitution $u = \sin\Phi$, $du = \cos\Phi \; d\Phi$, gives $\int_0^1 u^{2k} du = 1/(2k+1)$.)

This integral appears in the mean gravity and $J_{2n}$ computations.

### 2.5 Double Factorial Notation

$$n!! = \begin{cases} n(n-2)(n-4)\cdots 4 \cdot 2 & n \text{ even} \\ n(n-2)(n-4)\cdots 3 \cdot 1 & n \text{ odd} \end{cases}$$

**Connection to factorials:** $(2k)!! = 2^k k!$ and $(2k-1)!! = \frac{(2k)!}{2^k k!}$

**Connection to Wallis integrals:** $W_{2k} = \frac{(2k)!}{4^k (k!)^2}\cdot\frac{\pi}{2}$ and $W_{2k+1} = \frac{4^k (k!)^2}{(2k+1)!}$

---

## 3. Tracing the Coefficients

### 3.1 Meridian Quadrant $Q$ — Simplest Case

**Source:** [M80] Sec. 3

$$Q = c \int_0^{\pi/2} \frac{d\Phi}{(1 + e'^2\cos^2\Phi)^{3/2}} = c\frac{\pi}{2}\left(1 - \frac{3}{4}e'^2 + \frac{45}{64}e'^4 - \frac{175}{256}e'^6 + \frac{11025}{16384}e'^8\right)$$

**Step 1:** Expand $(1 + e'^2\cos^2\Phi)^{-3/2}$ with $\alpha = -3/2$, $x = e'^2\cos^2\Phi$:

$$(1 + e'^2\cos^2\Phi)^{-3/2} = \sum_{k=0}^{\infty} \binom{-3/2}{k} e'^{2k} \cos^{2k}\Phi$$

where $\binom{-3/2}{k} = (-1)^k \frac{(2k+1)!!}{2^k \cdot k! \cdot 2^0}$... more precisely:

$$\binom{-3/2}{k} = \frac{(-3/2)(-5/2)\cdots(-(2k+1)/2)}{k!} = \frac{(-1)^k \cdot 3 \cdot 5 \cdots (2k+1)}{2^k \cdot k!}$$

For $k \geq 1$: $\binom{-3/2}{k} = (-1)^k \frac{(2k+1)!!}{2^k \cdot k! \cdot 1}$ ... let me just tabulate:

| $k$ | $\binom{-3/2}{k}$ | Simplified |
|---|---|---|
| 0 | 1 | 1 |
| 1 | $-3/2$ | $-3/2$ |
| 2 | $(-3/2)(-5/2)/2! = 15/8$ | $15/8$ |
| 3 | $(-3/2)(-5/2)(-7/2)/3! = -35/16$ | $-35/16$ |
| 4 | $(-3/2)(-5/2)(-7/2)(-9/2)/4! = 315/128$ | $315/128$ |

**Step 2:** Integrate term-by-term using the **even-power Wallis integral**:

$$\int_0^{\pi/2} \cos^{2k}\Phi \; d\Phi = W_{2k} = \frac{(2k-1)!!}{(2k)!!} \cdot \frac{\pi}{2}$$

| $k$ | $W_{2k}$ |
|---|---|
| 0 | $\pi/2$ |
| 1 | $\frac{1}{2}\cdot\frac{\pi}{2}$ |
| 2 | $\frac{3}{8}\cdot\frac{\pi}{2}$ |
| 3 | $\frac{5}{16}\cdot\frac{\pi}{2}$ |
| 4 | $\frac{35}{128}\cdot\frac{\pi}{2}$ |

**Step 3:** Multiply binomial coefficient by Wallis integral (factoring out $\pi/2$):

| $k$ | $\binom{-3/2}{k} \times \frac{W_{2k}}{\pi/2}$ | Product | Published |
|---|---|---|---|
| 0 | $1 \times 1$ | $1$ | $1$ |
| 1 | $(-3/2)(1/2)$ | $-3/4$ | $-3/4$ ✓ |
| 2 | $(15/8)(3/8)$ | $45/64$ | $45/64$ ✓ |
| 3 | $(-35/16)(5/16)$ | $-175/256$ | $-175/256$ ✓ |
| 4 | $(315/128)(35/128)$ | $11025/16384$ | $11025/16384$ ✓ |

**General coefficient of $e'^{2k}$:**

$$c_k^{(Q)} = \binom{-3/2}{k} \cdot \frac{(2k-1)!!}{(2k)!!} = (-1)^k \frac{(2k+1)!! \cdot (2k-1)!!}{4^k (k!)^2}$$

### 3.2 Sphere of Equal Area $R_2$ — Requires Square Root

**Source:** [M80] Sec. 3, [NGA] B-15

$$R_2 = c\left(\int_0^{\pi/2} \frac{\cos\Phi \; d\Phi}{(1 + e'^2\cos^2\Phi)^2}\right)^{1/2}$$

Published: $R_2 = c(1 - \frac{2}{3}e'^2 + \frac{26}{45}e'^4 - \frac{100}{189}e'^6 + \frac{7034}{14175}e'^8)$

**Step 1:** Expand $(1 + e'^2\cos^2\Phi)^{-2}$ with $\alpha = -2$:

$$\binom{-2}{k} = (-1)^k(k+1)$$

So: $(1 + e'^2\cos^2\Phi)^{-2} = \sum_{k=0}^{\infty} (-1)^k(k+1) e'^{2k}\cos^{2k}\Phi$

**Step 2:** Integrate with $\cos\Phi$ weighting — this gives **odd-power** Wallis integrals:

$$\int_0^{\pi/2} \cos^{2k+1}\Phi \; d\Phi = W_{2k+1} = \frac{(2k)!!}{(2k+1)!!}$$

| $k$ | $W_{2k+1}$ |
|---|---|
| 0 | $1$ |
| 1 | $2/3$ |
| 2 | $8/15$ |
| 3 | $48/105 = 16/35$ |
| 4 | $384/945 = 128/315$ |

**Step 2 result:** The integral $I = \sum_{k=0}^{\infty} (-1)^k(k+1)\frac{(2k)!!}{(2k+1)!!} e'^{2k}$

| $k$ | $(-1)^k(k+1) \cdot W_{2k+1}$ | Value |
|---|---|---|
| 0 | $1 \cdot 1$ | $1$ |
| 1 | $-2 \cdot 2/3$ | $-4/3$ |
| 2 | $3 \cdot 8/15$ | $8/5$ |
| 3 | $-4 \cdot 16/35$ | $-64/35$ |
| 4 | $5 \cdot 128/315$ | $128/63$ |

So $I = 1 - \frac{4}{3}e'^2 + \frac{8}{5}e'^4 - \frac{64}{35}e'^6 + \frac{128}{63}e'^8 + \ldots$

**Step 3:** Take square root $R_2/c = \sqrt{I}$ using the generalized binomial with $\alpha = 1/2$:

$$\sqrt{1 + u} = 1 + \frac{1}{2}u - \frac{1}{8}u^2 + \frac{1}{16}u^3 - \frac{5}{128}u^4 + \ldots$$

where $u = -\frac{4}{3}e'^2 + \frac{8}{5}e'^4 - \frac{64}{35}e'^6 + \ldots$

**Collecting terms by power of $e'^2$:**

Coefficient of $e'^2$: $\frac{1}{2}(-\frac{4}{3}) = -\frac{2}{3}$ ✓

Coefficient of $e'^4$: $\frac{1}{2}(\frac{8}{5}) + (-\frac{1}{8})(-\frac{4}{3})^2 = \frac{4}{5} - \frac{2}{9} = \frac{36-10}{45} = \frac{26}{45}$ ✓

Coefficient of $e'^6$: $\frac{1}{2}(-\frac{64}{35}) + (-\frac{1}{8})\cdot 2 \cdot(-\frac{4}{3})(\frac{8}{5}) + \frac{1}{16}(-\frac{4}{3})^3$

$= -\frac{32}{35} + \frac{1}{8}\cdot\frac{64}{15} - \frac{1}{16}\cdot\frac{64}{27}$

$= -\frac{32}{35} + \frac{8}{15} - \frac{4}{27}$

$= \frac{-32 \cdot 27 + 8 \cdot 63 - 4 \cdot 35}{35 \cdot 27} = \frac{-864 + 504 - 140}{945} = \frac{-500}{945} = -\frac{100}{189}$ ✓

**Key insight:** The $R_2$ coefficients arise from composing TWO generalized binomial operations: first $\alpha = -2$ (for the integrand), then $\alpha = 1/2$ (for the square root).

### 3.3 Mean Gravity $\bar{\gamma}$ — Requires Series Quotient

**Source:** [M80] Sec. 3, [NGA] B-27

$$\bar{\gamma} = \gamma_e\left(1 + \frac{1}{6}e^2 + \frac{1}{3}k + \frac{59}{360}e^4 + \frac{5}{18}e^2 k + \ldots\right)$$

This is a ratio of two integrals, each involving generalized binomial expansions:

$$\bar{\gamma}/\gamma_e = \frac{\text{numerator}}{\text{denominator}}$$

**Denominator:** $\int_0^{\pi/2}\frac{\cos\Phi}{(1-e^2\sin^2\Phi)^2} d\Phi$

Using $(1-e^2\sin^2\Phi)^{-2}$ with rising factorial $(\alpha)_k$ for $\alpha=2$: $\frac{(2)_k}{k!} = k+1$

$$\text{denom} = \sum_{k=0}^{\infty}(k+1)e^{2k}\int_0^{\pi/2}\sin^{2k}\Phi\cos\Phi\;d\Phi = \sum_{k=0}^{\infty}\frac{k+1}{2k+1}e^{2k}$$

$= 1 + \frac{2}{3}e^2 + \frac{3}{5}e^4 + \frac{4}{7}e^6 + \ldots$

**Numerator:** Using $\gamma = \gamma_e(1+k\sin^2\Phi)(1-e^2\sin^2\Phi)^{-1/2}$:

$$\text{num} = \int_0^{\pi/2}\frac{(1+k\sin^2\Phi)\cos\Phi}{(1-e^2\sin^2\Phi)^{5/2}}d\Phi$$

Using $(1-x)^{-5/2}$ with rising factorial $(5/2)_n$:

$(5/2)_0 = 1, \; (5/2)_1 = 5/2, \; (5/2)_2 = 35/4, \; (5/2)_3 = 315/8$

$$\text{num} = \sum_{n=0}^{\infty}\frac{(5/2)_n}{n!}e^{2n}\left(\frac{1}{2n+1} + \frac{k}{2n+3}\right)$$

$= 1 + (\frac{5}{6} + \frac{k}{3})e^2 + (\frac{35}{40} + \frac{5k}{10})e^4 + \ldots$

$= 1 + (\frac{5}{6} + \frac{k}{3})e^2 + (\frac{7}{8} + \frac{k}{2})e^4 + \ldots$

**Step 3: Series quotient** $\bar{\gamma}/\gamma_e = \text{num}/\text{denom}$

Let num $= 1 + n_1 e^2 + n_2 e^4 + \ldots$ and denom $= 1 + d_1 e^2 + d_2 e^4 + \ldots$

Then $\text{num}/\text{denom} = 1 + (n_1 - d_1)e^2 + (n_2 - d_2 - d_1(n_1 - d_1))e^4 + \ldots$

$n_1 - d_1 = (\frac{5}{6} + \frac{k}{3}) - \frac{2}{3} = \frac{1}{6} + \frac{k}{3}$ ✓ (matches $\frac{1}{6}e^2 + \frac{1}{3}k$)

For the $e^4$ term: $n_2 - d_2 - d_1(n_1 - d_1) = (\frac{7}{8} + \frac{k}{2}) - \frac{3}{5} - \frac{2}{3}(\frac{1}{6}+\frac{k}{3})$

$= \frac{7}{8} + \frac{k}{2} - \frac{3}{5} - \frac{1}{9} - \frac{2k}{9}$

$= (\frac{7}{8} - \frac{3}{5} - \frac{1}{9}) + (\frac{1}{2} - \frac{2}{9})k$

$= \frac{7\cdot45 - 3\cdot72 - 1\cdot40}{360} + \frac{9-4}{18}k = \frac{315 - 216 - 40}{360} + \frac{5}{18}k = \frac{59}{360} + \frac{5}{18}k$ ✓

**The coefficients arise from: rising factorial → simple integral $1/(2k+1)$ → series long division.**

### 3.4 Even Zonal Harmonics $J_{2n}$ — Direct Formula

**Source:** [M80] Sec. 3, [HM67] (2-92)

$$J_{2n} = (-1)^{n+1} \frac{3e^{2n}}{(2n+1)(2n+3)}\left(1 - n + 5n\frac{J_2}{e^2}\right)$$

These are NOT series coefficients from integration — they come from the analytical solution of Clairaut's equation for the external potential of the level ellipsoid. The pattern $(2n+1)(2n+3)$ in the denominator relates to the eigenvalues of the Legendre equation.

### 3.5 Gravity Series Coefficients $a_{2n}$

**Source:** [M80] Sec. 3 "The Gravity Formula"

$$\gamma = \gamma_e\left(1 + \sum_{n=1}^{\infty} a_{2n}\sin^{2n}\Phi\right)$$

where:

| $n$ | $a_{2n}$ | Structure |
|---|---|---|
| 1 | $\frac{1}{2}e^2 + k$ | Linear in $e^2$ and $k$ |
| 2 | $\frac{3}{8}e^4 + \frac{1}{2}e^2 k$ | |
| 3 | $\frac{5}{16}e^6 + \frac{3}{8}e^4 k$ | |
| 4 | $\frac{35}{128}e^8 + \frac{5}{16}e^6 k$ | |

**Pattern:** The $e^{2n}$ coefficient in $a_{2n}$ is $\binom{-1/2}{n} \cdot (-1)^n = \frac{(2n-1)!!}{(2n)!!} = \frac{(2n)!}{4^n(n!)^2}$ — these are the **central binomial coefficients divided by $4^n$** (equivalently, the Maclaurin coefficients of $(1-x)^{-1/2}$).

The $e^{2(n-1)}k$ coefficient is the same pattern shifted by one index.

These arise from expanding $\gamma_e(1+k\sin^2\Phi)/(1-e^2\sin^2\Phi)^{1/2}$ as a power series in $\sin^2\Phi$, then collecting powers.

---

## 4. Summary: Mathematical Structures Needed

### 4.1 Core Functions for `math/` Library

| Function | Purpose | Used In |
|---|---|---|
| `FallingFactorial<T>(alpha, k)` | $\alpha^{(k)} = \alpha(\alpha-1)\cdots(\alpha-k+1)$ | Generalized binomial |
| `RisingFactorial<T>(alpha, k)` | $(\alpha)_k = \alpha(\alpha+1)\cdots(\alpha+k-1)$ | $(1-x)^{-\alpha}$ expansion |
| `GeneralizedBinomial<T>(alpha, k)` | $\binom{\alpha}{k} = \alpha^{(k)}/k!$ | All series |
| `DoubleFactorial<T>(n)` | $n!! = n(n-2)(n-4)\cdots$ | Wallis integrals |
| `WallisIntegral<T>(n)` | $\int_0^{\pi/2}\cos^n\Phi\;d\Phi$ | All ellipsoidal integrals |

### 4.2 Series Operations for `math/` Library

| Operation | Purpose | Used In |
|---|---|---|
| `SeriesSqrt(coeffs)` | $\sqrt{1 + c_1 x + c_2 x^2 + \ldots}$ | $R_2$ |
| `SeriesReciprocal(coeffs)` | $1/(1 + c_1 x + c_2 x^2 + \ldots)$ | $\bar{\gamma}$ denominator |
| `SeriesMultiply(a, b)` | Cauchy product | Combining series |
| `BinomialExpansion(alpha, x_series, n_terms)` | $(1+u)^\alpha$ where $u$ is a series | General composition |

### 4.3 Design Principle: Compute Coefficients, Don't Hard-Code Them

For arbitrary precision, the series coefficients should be **computed** from the generalized binomial coefficients and Wallis integrals to the full precision of type `T`, not hard-coded as `double` literals. The computation sequence:

1. Compute $e^2$ (or $e'^2$) to full precision
2. For each series needed ($Q$, $R_2$, $\bar{\gamma}$, $a_{2n}$, etc.):
   a. Compute generalized binomial coefficients $\binom{\alpha}{k}$ to full precision
   b. Compute Wallis integrals $W_n$ to full precision (these are exact rationals)
   c. Multiply to get the integral's series coefficients
   d. Apply series algebra (sqrt, reciprocal, quotient) to full precision
3. Evaluate the series at the specific $e^2$ or $e'^2$ value

The Wallis integrals and binomial coefficients are exact rational numbers — they involve only integer arithmetic (factorials, double factorials). This means they can be computed to **infinite precision** without any floating-point error, then multiplied by the (irrational) powers of $e^2$.

### 4.4 Connection Between the Three Structures

```
Generalized Binomial ─── uses ──→ Falling Factorial
     │                                    │
     │ produces coefficients               │ related by (-1)^k
     │                                    │
     ▼                                    ▼
Integrand Series ──── evaluates via ──→ Rising Factorial
     │                                (for (1-x)^{-α} form)
     │ term-by-term integration
     │
     ▼
Wallis Integrals ──── evaluate as ──→ Double Factorials
     │                                (exact rationals)
     │
     ▼
Raw Series Coefficients
     │
     │ series sqrt / reciprocal / quotient
     │ (another generalized binomial application)
     ▼
Published Coefficients
```
