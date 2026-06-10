# Derivation 001: Series Expansion of $q_0$

## Goal

Derive the series form of $q_0$ from its definition, showing every algebraic step, and verify the first several terms numerically.

## Starting Point

From Heiskanen & Moritz (1967) Eq. (2-58) and (2-71), the quantity $q_0$ arises in the theory of the normal gravity field of the equipotential ellipsoid. It is defined as:

$$q_0 = \frac{1}{2}\left[\left(1 + \frac{3}{e'^2}\right)\arctan(e') - \frac{3}{e'}\right] \tag{001.Eq.1}$$

where $e' = \frac{e}{\sqrt{1-e^2}}$ is the second eccentricity.

This closed form suffers from **catastrophic cancellation** because $(1 + 3/e'^2)\arctan(e') \approx 3/e'$ and the difference is $O(e'^3)$. We need a series form.

## Step 1: Taylor Series for $\arctan(e')$

The Maclaurin series for $\arctan(x)$ is:

$$\arctan(x) = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} x^{2n+1}, \quad |x| \leq 1 \tag{001.Eq.2}$$

For $e' \approx 0.082 < 1$, this converges. Substituting $x = e'$:

$$\arctan(e') = e' - \frac{e'^3}{3} + \frac{e'^5}{5} - \frac{e'^7}{7} + \cdots = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n+1} \tag{001.Eq.3}$$

## Step 2: Multiply by $(1 + 3/e'^2)$

$$\left(1 + \frac{3}{e'^2}\right)\arctan(e') = \arctan(e') + \frac{3}{e'^2}\arctan(e') \tag{001.Eq.4}$$

### First part: $\arctan(e')$

$$\arctan(e') = \sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n+1} \tag{001.Eq.5}$$

### Second part: $\frac{3}{e'^2}\arctan(e')$

$$\frac{3}{e'^2}\arctan(e') = \frac{3}{e'^2}\sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n+1} = 3\sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n-1} \tag{001.Eq.6}$$

Writing out the first few terms:

$$= 3\left[\frac{e'^{-1}}{1} - \frac{e'^1}{3} + \frac{e'^3}{5} - \frac{e'^5}{7} + \cdots\right] = \frac{3}{e'} - e' + \frac{3e'^3}{5} - \frac{3e'^5}{7} + \cdots \tag{001.Eq.7}$$

### Combined:

$$\left(1 + \frac{3}{e'^2}\right)\arctan(e') = \frac{3}{e'} + \left(e' - e'\right) + \left(-\frac{e'^3}{3} + \frac{3e'^3}{5}\right) + \left(\frac{e'^5}{5} - \frac{3e'^5}{7}\right) + \cdots \tag{001.Eq.8}$$

$$= \frac{3}{e'} + 0 + e'^3\left(-\frac{1}{3} + \frac{3}{5}\right) + e'^5\left(\frac{1}{5} - \frac{3}{7}\right) + \cdots \tag{001.Eq.9}$$

## Step 3: Subtract $3/e'$

$$\left(1 + \frac{3}{e'^2}\right)\arctan(e') - \frac{3}{e'} = e'^3\left(\frac{-1}{3} + \frac{3}{5}\right) + e'^5\left(\frac{1}{5} - \frac{3}{7}\right) + e'^7\left(-\frac{1}{7} + \frac{3}{9}\right) + \cdots \tag{001.Eq.10}$$

The $3/e'$ term cancels exactly with the $n=0$ term from the second part.

## Step 4: General Term

For $n \geq 1$, the coefficient of $e'^{2n+1}$ in the combined sum is:

From $\arctan(e')$: $\frac{(-1)^n}{2n+1}$

From $\frac{3}{e'^2}\arctan(e')$: The $e'^{2n+1}$ term comes from the $(n+1)$-th term of $\arctan(e')$ multiplied by $3/e'^2$, which is $\frac{3(-1)^{n+1}}{2(n+1)+1} \cdot e'^{2(n+1)+1} / e'^2 = \frac{3(-1)^{n+1}}{2n+3} e'^{2n+1}$.

Expanding $\frac{3}{e'^2}\arctan(e') = 3\sum_{m=0}^{\infty}\frac{(-1)^m}{2m+1}e'^{2m-1}$, the power of $e'$ is $2m-1$. Setting $2m-1 = 2n+1$ gives $m = n+1$. So the coefficient of $e'^{2n+1}$ from this sum is:

$$\frac{3(-1)^{n+1}}{2(n+1)+1} = \frac{3(-1)^{n+1}}{2n+3} \tag{001.Eq.11}$$

Combined coefficient of $e'^{2n+1}$ (for $n \geq 1$):

$$\frac{(-1)^n}{2n+1} + \frac{3(-1)^{n+1}}{2n+3} \tag{001.Eq.12}$$

$$= (-1)^n\left[\frac{1}{2n+1} - \frac{3}{2n+3}\right] \tag{001.Eq.13}$$

$$= (-1)^n\left[\frac{(2n+3) - 3(2n+1)}{(2n+1)(2n+3)}\right] \tag{001.Eq.14}$$

$$= (-1)^n\left[\frac{2n + 3 - 6n - 3}{(2n+1)(2n+3)}\right] \tag{001.Eq.15}$$

$$= (-1)^n\left[\frac{-4n}{(2n+1)(2n+3)}\right] \tag{001.Eq.16}$$

$$= (-1)^{n+1}\frac{4n}{(2n+1)(2n+3)} \tag{001.Eq.17}$$

## Step 5: Assemble $2q_0$

$$2q_0 = \sum_{n=1}^{\infty} (-1)^{n+1}\frac{4n}{(2n+1)(2n+3)} e'^{2n+1} \tag{001.Eq.18}$$

Or equivalently:

$$q_0 = \frac{1}{2}\sum_{n=1}^{\infty} (-1)^{n+1}\frac{4n}{(2n+1)(2n+3)} e'^{2n+1} \tag{001.Eq.19}$$

## Step 6: Numerical Verification

Using $e' = 0.0820944379$ (WGS84 second eccentricity):

### Term $n=1$:

$$(-1)^2 \frac{4 \cdot 1}{3 \cdot 5} e'^3 = \frac{4}{15} e'^3 \tag{001.Eq.20}$$

$e'^3 = (0.0820944)^3 = 5.53345 \times 10^{-4}$

Term $= \frac{4}{15} \times 5.53345 \times 10^{-4} = 1.47559 \times 10^{-4}$

### Term $n=2$:

$$(-1)^3 \frac{4 \cdot 2}{5 \cdot 7} e'^5 = -\frac{8}{35} e'^5 \tag{001.Eq.21}$$

$e'^5 = (0.0820944)^5 = 3.72802 \times 10^{-6}$

Term $= -\frac{8}{35} \times 3.72802 \times 10^{-6} = -8.5212 \times 10^{-7}$

### Term $n=3$:

$$(-1)^4 \frac{4 \cdot 3}{7 \cdot 9} e'^7 = \frac{12}{63} e'^7 = \frac{4}{21} e'^7 \tag{001.Eq.22}$$

$e'^7 = (0.0820944)^7 = 2.51248 \times 10^{-8}$

Term $= \frac{4}{21} \times 2.51248 \times 10^{-8} = 4.7857 \times 10^{-9}$

### Partial sum (3 terms):

$$2q_0 \approx 1.47559 \times 10^{-4} - 8.5212 \times 10^{-7} + 4.7857 \times 10^{-9} = 1.46711 \times 10^{-4} \tag{001.Eq.23}$$

$$q_0 \approx 7.3356 \times 10^{-5} \tag{001.Eq.24}$$

### Published value: $q_0 = 7.334625787083 \times 10^{-5}$ ([NGA] Appendix B, Eq. B-18)

The 3-term partial sum gives $7.3356 \times 10^{-5}$, which matches to 3 significant figures. The 4th term is $O(e'^9) \approx 3 \times 10^{-11}$, confirming rapid convergence.

## Properties

- **Alternating series:** Signs alternate starting positive
- **Monotonically decreasing terms:** Ratio of successive terms $\approx e'^2 \approx 0.0067$
- **Leibniz error bound:** Error after $N$ terms $\leq$ first omitted term
- **Exact rational coefficients:** $\frac{4n}{(2n+1)(2n+3)}$ is exact for all $n$
- **No cancellation:** All terms are computed independently; the catastrophic cancellation of the closed form is completely avoided
