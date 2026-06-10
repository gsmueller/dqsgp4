# Derivation 002: Series Expansion of $q_0'$

## Goal

Derive the series form of $q_0'$ from its definition, showing every algebraic step.

## Starting Point

From Heiskanen & Moritz (1967) Eq. (2-67), with $u = b$ for the level ellipsoid:

$$q_0' = 3\left[\left(1 + \frac{1}{e'^2}\right)\left(1 - \frac{1}{e'}\arctan(e')\right)\right] - 1 \tag{002.Eq.1}$$

This also suffers from cancellation: $(1/e')\arctan(e') \approx 1$ for small $e'$.

## Step 1: Series for $\frac{1}{e'}\arctan(e')$

From Derivation 001, Step 1:

$$\arctan(e') = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n+1} \tag{002.Eq.2}$$

Dividing by $e'$:

$$\frac{1}{e'}\arctan(e') = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n} = 1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \frac{e'^6}{7} + \cdots \tag{002.Eq.3}$$

## Step 2: $1 - \frac{1}{e'}\arctan(e')$

$$1 - \frac{1}{e'}\arctan(e') = 1 - \left(1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \cdots\right) \tag{002.Eq.4}$$

$$= \frac{e'^2}{3} - \frac{e'^4}{5} + \frac{e'^6}{7} - \cdots \tag{002.Eq.5}$$

$$= \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} \tag{002.Eq.6}$$

**Verification:** For $n=1$: $\frac{(-1)^2}{3}e'^2 = \frac{e'^2}{3}$ ✓. For $n=2$: $\frac{(-1)^3}{5}e'^4 = -\frac{e'^4}{5}$ ✓.

## Step 3: Multiply by $(1 + 1/e'^2)$

$$\left(1 + \frac{1}{e'^2}\right) \times \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} \tag{002.Eq.7}$$

$$= \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} + \frac{1}{e'^2}\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} \tag{002.Eq.8}$$

$$= \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n} + \sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n-2} \tag{002.Eq.9}$$

### Reindex the second sum

Let $m = n-1$ (so $n = m+1$, and when $n=1$, $m=0$):

$$\sum_{n=1}^{\infty} \frac{(-1)^{n+1}}{2n+1} e'^{2n-2} = \sum_{m=0}^{\infty} \frac{(-1)^{m+2}}{2(m+1)+1} e'^{2m} = \sum_{m=0}^{\infty} \frac{(-1)^m}{2m+3} e'^{2m} \tag{002.Eq.10}$$

Separating the $m=0$ term:

$$= \frac{1}{3} + \sum_{m=1}^{\infty} \frac{(-1)^m}{2m+3} e'^{2m} \tag{002.Eq.11}$$

### Combine both sums

The coefficient of $e'^{2n}$ for $n \geq 1$ is:

From first sum: $\frac{(-1)^{n+1}}{2n+1}$

From second sum: $\frac{(-1)^n}{2n+3}$

Combined: $\frac{(-1)^{n+1}}{2n+1} + \frac{(-1)^n}{2n+3} = (-1)^{n+1}\left[\frac{1}{2n+1} - \frac{1}{2n+3}\right]$

$$= (-1)^{n+1} \cdot \frac{(2n+3) - (2n+1)}{(2n+1)(2n+3)} = (-1)^{n+1} \cdot \frac{2}{(2n+1)(2n+3)} \tag{002.Eq.12}$$

Plus the constant $1/3$ from the $m=0$ term.

So:

$$\left(1 + \frac{1}{e'^2}\right)\left(1 - \frac{1}{e'}\arctan(e')\right) = \frac{1}{3} + \sum_{n=1}^{\infty} (-1)^{n+1}\frac{2}{(2n+1)(2n+3)} e'^{2n} \tag{002.Eq.13}$$

## Step 4: Multiply by 3 and subtract 1

$$q_0' = 3\left[\frac{1}{3} + \sum_{n=1}^{\infty} (-1)^{n+1}\frac{2}{(2n+1)(2n+3)} e'^{2n}\right] - 1 \tag{002.Eq.14}$$

$$= 1 + \sum_{n=1}^{\infty} (-1)^{n+1}\frac{6}{(2n+1)(2n+3)} e'^{2n} - 1 \tag{002.Eq.15}$$

$$\boxed{q_0' = \sum_{n=1}^{\infty} \frac{6(-1)^{n+1}}{(2n+1)(2n+3)} e'^{2n}} \tag{002.Eq.16}$$

## Step 5: Numerical Verification

Using $e' = 0.0820944379$ (WGS84):

### Term $n=1$:

$$\frac{6}{3 \times 5} e'^2 = \frac{6}{15} e'^2 = \frac{2}{5} e'^2 \tag{002.Eq.17}$$

$e'^2 = 0.006739497$

Term $= 0.4 \times 0.006739497 = 2.69580 \times 10^{-3}$

### Term $n=2$:

$$-\frac{6}{5 \times 7} e'^4 = -\frac{6}{35} e'^4 \tag{002.Eq.18}$$

$e'^4 = 4.54208 \times 10^{-5}$

Term $= -0.17143 \times 4.54208 \times 10^{-5} = -7.78642 \times 10^{-6}$

### Term $n=3$:

$$\frac{6}{7 \times 9} e'^6 = \frac{6}{63} e'^6 = \frac{2}{21} e'^6 \tag{002.Eq.19}$$

$e'^6 = 3.06122 \times 10^{-7}$

Term $= 0.09524 \times 3.06122 \times 10^{-7} = 2.91545 \times 10^{-8}$

### Partial sum:

$$q_0' \approx 2.69580 \times 10^{-3} - 7.78642 \times 10^{-6} + 2.91545 \times 10^{-8} = 2.68804 \times 10^{-3} \tag{002.Eq.20}$$

### Published value: $q_0' = 2.688041300461 \times 10^{-3}$ ([NGA] Appendix B, Eq. B-19)

3-term sum: $2.68804 \times 10^{-3}$. Matches to 4 significant figures. ✓

## Comparison: $q_0$ vs $q_0'$ Series

| Property | $q_0$ | $q_0'$ |
|---|---|---|
| Series variable | $e'^{2n+1}$ (odd powers) | $e'^{2n}$ (even powers) |
| Rational coefficient | $\frac{4n}{(2n+1)(2n+3)}$ | $\frac{6}{(2n+1)(2n+3)}$ |
| Starts at | $n=1$ | $n=1$ |
| First term | $\frac{4}{15}e'^3$ | $\frac{2}{5}e'^2$ |
| Convergence ratio | $e'^2 \approx 0.0067$ | $e'^2 \approx 0.0067$ |
| Alternating | Yes | Yes |

Both share the denominator $(2n+1)(2n+3)$ and converge at the same rate. The denominators arise from the Wallis-type structure of the $\arctan$ expansion.

## Dependency Chain

$$\text{Defining: } 1/f \xrightarrow{f = 1/(1/f)} f \xrightarrow{e^2 = 2f-f^2} e^2 \xrightarrow{e' = e/\sqrt{1-e^2}} e' \xrightarrow{\text{this series}} q_0' \tag{002.Eq.21}$$

Every step is a computation, not a lookup. The series IS the definition of $q_0'$ in terms of $e'$.
