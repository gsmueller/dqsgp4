# Derivation 006: Even Zonal Harmonics $J_{2n}$ from the Equipotential Ellipsoid

## Goal

Show how the even zonal harmonics $J_2, J_4, J_6, \ldots$ arise from the equipotential ellipsoid theory, and derive the closed formula that expresses $J_{2n}$ in terms of $e^2$, $m$, and $q_0$.

## Context: Where $J_{2n}$ Come From

The gravitational potential outside the Earth can be expanded in spherical harmonics. For an axially symmetric body (like a level ellipsoid), only zonal harmonics ($m=0$) are nonzero:

$$V = \frac{GM}{r}\left[1 - \sum_{n=2}^{\infty} J_n \left(\frac{a}{r}\right)^n P_n(\sin\varphi')\right]$$

For a **level ellipsoid** (equipotential surface with constant rotation), the odd zonals ($J_3, J_5, \ldots$) are exactly zero by symmetry. The even zonals ($J_2, J_4, J_6, \ldots$) are completely determined by the four defining parameters through the requirement that the ellipsoid surface be an equipotential.

This is the key point: **$J_{2n}$ are not independent measurements for an equipotential ellipsoid**. They are computed quantities. The ellipsoid's shape ($e^2$) and rotation ($m$) uniquely determine every $J_{2n}$.

## Step 1: The Normal Potential on the Ellipsoid Surface

From Heiskanen & Moritz (1967), the external gravitational potential of the level ellipsoid in ellipsoidal coordinates $(u, \beta)$ is:

$$V = \frac{GM}{E}\arctan\frac{E}{u} + \frac{\omega^2 a^2}{2}\frac{q}{q_0}\left(\sin^2\beta - \frac{1}{3}\right)$$

where $q$ and $q_0$ are the functions from Derivations 001-002, and $E = ae$ is the linear eccentricity.

Converting this to the spherical harmonic form (Heiskanen & Moritz, (2-90) through (2-92')):

## Step 2: Expansion of $V$ in Powers of $a/r$

The $\arctan(E/u)$ term, when expanded for $r > a$ (exterior), produces a series in $(a/r)^n$. Each power of $(a/r)^n$ multiplied by $P_n(\cos\theta)$ gives the $J_n$ contribution.

From Heiskanen & Moritz (2-92) and (2-92'), the result is:

$$J_{2n} = \frac{(-1)^{n+1} 3 e^{2n}}{(2n+1)(2n+3)}\left(1 - n + 5n\frac{J_2}{e^2}\right)$$

## Step 3: Verification That This is Self-Consistent

This formula expresses $J_{2n}$ in terms of $J_2$ and $e^2$. But $J_2$ itself is related to $e^2$ and $m$ by the fundamental equation (Derivation 001 context):

$$J_2 = \frac{e^2}{3}\left(1 - \frac{2me'}{15q_0}\right)$$

So the full chain is:

$$1/f \to e^2 \to e' \to q_0 \to J_2 \to J_{2n}$$

and $m$ enters through $GM$ and $\omega$.

### Check for $n=1$ ($J_2$ from $J_2$):

$$J_2 = \frac{(-1)^2 \cdot 3e^2}{3 \times 5}\left(1 - 1 + 5\frac{J_2}{e^2}\right) = \frac{3e^2}{15} \cdot \frac{5J_2}{e^2} = \frac{3 \cdot 5 \cdot J_2}{15} = J_2 \quad \checkmark$$

The formula is tautologically consistent for $n=1$.

### Check for $n=2$ ($J_4$):

$$J_4 = \frac{(-1)^3 \cdot 3e^4}{5 \times 7}\left(1 - 2 + \frac{10J_2}{e^2}\right) = -\frac{3e^4}{35}\left(-1 + \frac{10J_2}{e^2}\right)$$

Using $J_2/e^2 \approx 1.0826 \times 10^{-3}/6.694 \times 10^{-3} \approx 0.16172$:

Inner: $-1 + 10 \times 0.16172 = -1 + 1.6172 = 0.6172$

$J_4 = -\frac{3 \times (6.694 \times 10^{-3})^2}{35} \times 0.6172 = -\frac{3 \times 4.481 \times 10^{-5}}{35} \times 0.6172$

$= -3.841 \times 10^{-6} \times 0.6172 = -2.371 \times 10^{-6}$

**Published (Moritz 1980):** $J_4 = -0.00000237091222$. Our calculation: $-2.371 \times 10^{-6}$. ✓

### Check for $n=3$ ($J_6$):

$$J_6 = \frac{(-1)^4 \cdot 3e^6}{7 \times 9}\left(1 - 3 + \frac{15J_2}{e^2}\right) = \frac{3e^6}{63}\left(-2 + \frac{15 \times 0.16172}{1}\right) = \frac{e^6}{21}(0.4258)$$

$e^6 = (6.694 \times 10^{-3})^3 = 2.998 \times 10^{-7}$

$J_6 = \frac{2.998 \times 10^{-7}}{21} \times 0.4258 = 6.077 \times 10^{-9}$

**Published:** $J_6 = 0.00000000608347$. Our: $6.077 \times 10^{-9}$. ✓ (to 3 figures)

## Three Errors for $J_{2n}$

**Measurement error:** $J_{2n}$ depends on $J_2$, which depends on $m$, which depends on $GM$.

$$\sigma_m(J_{2n}) = \left|\frac{\partial J_{2n}}{\partial J_2}\right| \sigma_m(J_2) = \frac{(-1)^{n+1} 3e^{2n}}{(2n+1)(2n+3)} \cdot \frac{5n}{e^2} \cdot \sigma_m(J_2)$$

And $\sigma_m(J_2)$ propagates from $\sigma_m(m)$ which propagates from $\sigma_m(GM)$.

For $J_4$: $\sigma_m(J_4) \approx |3e^4/(35)| \times |10/e^2| \times \sigma_m(J_2) \approx 3.84 \times 10^{-6} \times 1494 \times \sigma_m(J_2)$

This is a significant amplification — the higher $J_{2n}$ are increasingly sensitive to measurement error in $J_2$.

**Precision error:** From the evaluation of the formula — multiplication and division of tracked values. Small.

**Accuracy error:** **Zero for the geometric ellipsoid** — this is an exact formula, not a truncation. However, the PHYSICAL Earth's $J_{2n}$ differ from the ellipsoidal $J_{2n}$ because the Earth is not a perfect level ellipsoid. The difference is what EGM2008 captures. So when comparing ellipsoidal $J_{2n}$ to measured $J_{2n}$, there's a "model accuracy" issue — but within the ellipsoid theory itself, $\delta_a = 0$.

## The Computation Module

$J_{2n}$ should NOT be stored as constants. They should be computed by a function:

```
J2n(n, e_squared, J2) → TrackedValue<T>
```

that implements the formula above. The function IS the definition. Its inputs are $e^2$ (from the geometric chain) and $J_2$ (from $e^2$, $m$, $q_0$ — the physical chain). Nothing is looked up.

## Equations to Resolve

| Equation | Status | What must be shown |
|---|---|---|
| J₂ₙ closed formula from ellipsoidal harmonic expansion | NOT DERIVED | Reproduce arctan(E/u) expansion → spherical harmonic matching |
| Ellipsoidal-to-spherical coordinate conversion | NOT DERIVED | Show u, β_geo, λ → r, φ, λ mapping |
| Self-consistency check n=1 | VERIFIED | J₂ from J₂ (tautological) |
| n=2 check (J₄) | VERIFIED | -2.371×10⁻⁶ matches Moritz |
| n=3 check (J₆) | VERIFIED | 6.077×10⁻⁹ matches published |

## Source Documents Required

| Source | Location | Availability | What it provides |
|---|---|---|---|
| Heiskanen & Moritz (1967) Ch. 2 | External | NOT in repo | Primary source for ellipsoidal potential theory |
| GRS80 Moritz (1980) | wgs84_markdown_parts/40_GRS80_Moritz.md | ✓ | Reference values for J₂, J₄, J₆ verification |
| WGS84 NGA standard | wgs84_markdown_parts/ | ✓ | Defining parameters and published zonal harmonics |
