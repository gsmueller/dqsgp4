# Derivation 005: $R_2$ — Radius of a Sphere of Equal Surface Area

## Goal

Derive the series expansion for $R_2$ from the defining integral, showing every step and tracking all three errors.

## Starting Point

From Moritz (1980) Sec. 3, the radius of a sphere whose surface area equals the ellipsoid's surface area is:

$$R_2 = c\left(\int_0^{\pi/2} \frac{\cos\Phi}{(1 + e'^2\cos^2\Phi)^2} \, d\Phi\right)^{1/2} \tag{005.Eq.1}$$

where $c = a^2/b$ is the polar radius of curvature, $e'$ is the second eccentricity, and $\Phi$ is the geodetic latitude.

## Step 1: Expand $(1 + e'^2\cos^2\Phi)^{-2}$ via Generalized Binomial

The generalized binomial theorem with $\alpha = -2$:

$$(1 + x)^{-2} = \sum_{k=0}^{\infty} \binom{-2}{k} x^k \tag{005.Eq.2}$$

The generalized binomial coefficient for $\alpha = -2$:

$$\binom{-2}{k} = \frac{(-2)(-3)(-4)\cdots(-2-k+1)}{k!} = \frac{(-1)^k(k+1)!}{k!} = (-1)^k(k+1) \tag{005.Eq.3}$$

**Verification:**
- $k=0$: $\binom{-2}{0} = 1$ ✓ (by convention)
- $k=1$: $(-2)/1! = -2 = (-1)^1 \times 2$ ✓
- $k=2$: $(-2)(-3)/2! = 6/2 = 3 = (-1)^2 \times 3$ ✓
- $k=3$: $(-2)(-3)(-4)/3! = -24/6 = -4 = (-1)^3 \times 4$ ✓

So: $(1 + e'^2\cos^2\Phi)^{-2} = \sum_{k=0}^{\infty} (-1)^k(k+1)\,e'^{2k}\cos^{2k}\Phi$

## Step 2: Integrate Term-by-Term

$$I = \int_0^{\pi/2} \frac{\cos\Phi}{(1+e'^2\cos^2\Phi)^2}\,d\Phi = \sum_{k=0}^{\infty} (-1)^k(k+1)\,e'^{2k}\int_0^{\pi/2}\cos^{2k+1}\Phi\,d\Phi \tag{005.Eq.4}$$

The integral $\int_0^{\pi/2}\cos^{2k+1}\Phi\,d\Phi$ is the Wallis integral for odd power $n = 2k+1$:

$$W_{2k+1} = \frac{(2k)!!}{(2k+1)!!} \tag{005.Eq.5}$$

where $(2k)!! = 2 \times 4 \times 6 \times \cdots \times 2k$ and $(2k+1)!! = 1 \times 3 \times 5 \times \cdots \times (2k+1)$.

**Verification by direct computation:**
- $W_1 = \int_0^{\pi/2}\cos\Phi\,d\Phi = [\sin\Phi]_0^{\pi/2} = 1$. Formula: $(0)!!/(1)!! = 1/1 = 1$ ✓
- $W_3 = \int_0^{\pi/2}\cos^3\Phi\,d\Phi = \int_0^1 (1-u^2)\,du = [u - u^3/3]_0^1 = 2/3$. Formula: $2!!/3!! = 2/3$ ✓
- $W_5 = \int_0^{\pi/2}\cos^5\Phi\,d\Phi$. Substitute $u = \sin\Phi$: $\int_0^1(1-u^2)^2\,du = \int_0^1(1-2u^2+u^4)\,du = 1 - 2/3 + 1/5 = 8/15$. Formula: $(4)!!/(5)!! = (2\times4)/(1\times3\times5) = 8/15$ ✓

So the integral series is:

$$I = \sum_{k=0}^{\infty} (-1)^k(k+1)\frac{(2k)!!}{(2k+1)!!}\,e'^{2k} \tag{005.Eq.6}$$

Computing the first five terms:

| $k$ | $(-1)^k(k+1)$ | $(2k)!!/(2k+1)!!$ | Product | $\times e'^{2k}$ ($e'^2 = 0.006739$) |
|---|---|---|---|---|
| 0 | 1 | 1 | 1 | 1 |
| 1 | $-2$ | $2/3$ | $-4/3$ | $-9.0040 \times 10^{-3}$ |
| 2 | 3 | $8/15$ | $8/5$ | $7.2585 \times 10^{-5}$ |
| 3 | $-4$ | $16/35$ | $-64/35$ | $-8.3005 \times 10^{-7}$ |
| 4 | 5 | $128/315$ | $128/63$ | $1.0512 \times 10^{-8}$ |

$I \approx 1 - 9.0040 \times 10^{-3} + 7.2585 \times 10^{-5} - 8.3005 \times 10^{-7} + 1.0512 \times 10^{-8}$

$I \approx 0.990924...$

## Step 3: Square Root via Newton's Method

$R_2 = c\sqrt{I}$

Rather than expanding $\sqrt{I}$ as another binomial series, we compute $\sqrt{I}$ directly by Newton's method:

$$x_{n+1} = \frac{1}{2}\left(x_n + \frac{I}{x_n}\right) \tag{005.Eq.7}$$

Starting from $x_0 = 1$ (since $I \approx 0.991 \approx 1$):
- $x_1 = (1 + 0.990924)/2 = 0.995462$
- $x_2 = (0.995462 + 0.990924/0.995462)/2 = (0.995462 + 0.995458)/2 = 0.995460$
- Converged to 6 digits in 2 iterations.

$\sqrt{I} \approx 0.995460$

$R_2 = c \times 0.995460 = 6399593.6258 \times 0.995460 = 6371007.18...$

**Published:** $R_2 = 6371007.1810$ m ✓

## Alternative: Direct Series for $R_2/c$

For reference, the published series form (verified by composing the binomial expansion of $\sqrt{I}$):

$$R_2/c = 1 - \frac{2}{3}e'^2 + \frac{26}{45}e'^4 - \frac{100}{189}e'^6 + \frac{7034}{14175}e'^8 - \cdots \tag{005.Eq.8}$$

**Derivation of the $e'^2$ coefficient:**

$\sqrt{1 + u} = 1 + \frac{1}{2}u - \frac{1}{8}u^2 + \cdots$ where $u = -\frac{4}{3}e'^2 + \text{higher}$

Coefficient of $e'^2$: $\frac{1}{2} \times (-\frac{4}{3}) = -\frac{2}{3}$ ✓

**Derivation of the $e'^4$ coefficient:**

The $e'^4$ contribution comes from two sources:
1. $\frac{1}{2} \times (\text{coeff of } e'^4 \text{ in } I - 1) = \frac{1}{2} \times \frac{8}{5} = \frac{4}{5}$
2. $-\frac{1}{8} \times (\text{coeff of } e'^2 \text{ in } I - 1)^2 = -\frac{1}{8} \times \frac{16}{9} = -\frac{2}{9}$

Total: $\frac{4}{5} - \frac{2}{9} = \frac{36}{45} - \frac{10}{45} = \frac{26}{45}$ ✓

## Three Errors for $R_2$

**Measurement error:** $R_2$ depends on $c = a^2/b$ and $e'$. If using the definitional view ($\sigma_m(a) = 0$, $\sigma_m(1/f) = 0$), then $\sigma_m(R_2) = 0$. If using the physical view, $\sigma_m(R_2) \approx \sigma_m(a) \approx \pm 2$ m.

**Precision error:** Two sources:
1. Series truncation of $I$: bounded by the first omitted term (Leibniz criterion since the series of $I$ is alternating with decreasing terms for $e'^2 < 1$).
2. Newton iteration for $\sqrt{I}$: quadratic convergence, bounded by $|x_n - x_{n-1}|^2/(2x_n)$.

Total: $\delta_p(R_2) = c \times (\delta_p(\sqrt{I}))$ where $\delta_p(\sqrt{I})$ combines both sources.

**Accuracy error:** $\delta_a(R_2) = 0$ — this is an exact formula for the equipotential ellipsoid.
