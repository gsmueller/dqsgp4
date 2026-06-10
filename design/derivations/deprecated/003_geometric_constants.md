# Derivation 003: Geometric Constants from Defining Parameters

## Goal

Derive $f$, $e^2$, $e'$, $b$, $E$ from the two geometric defining parameters ($a$, $1/f$), showing every step and tracking all three errors through each operation.

## Inputs

| Symbol | Value | Measurement | Precision | Accuracy |
|---|---|---|---|---|
| $a$ | $6378137.0$ m | 0 (defined) | $\delta_p(a)$ = representation error of integer 6378137 in T (zero for any reasonable T) | 0 |
| $1/f$ | $298.257223563$ | 0 (defined) | $\delta_p(1/f)$ = representation error of this decimal in T | 0 |

## Step 1: $f = 1/(1/f)$ — one division

$$f = \frac{1}{1/f} \tag{003.Eq.1}$$

**Value:** $f = 1 / 298.257223563 = 3.3528106647474807...\times 10^{-3}$

**Measurement error:** $\sigma_m(f) = 0$ (both inputs are defined, $\sigma_m = 0$)

**Precision error:** Division of exact 1 by $(1/f)$:

$$\delta_p(f) = \frac{|1| \cdot \delta_p(1/f) + |1/f| \cdot \delta_p(1)}{|1/f|^2 - |1/f| \cdot \delta_p(1/f)} \tag{003.Eq.2}$$

Since $\delta_p(1) = 0$ (integer 1 is exact):

$$\delta_p(f) = \frac{\delta_p(1/f)}{|1/f|^2 - |1/f| \cdot \delta_p(1/f)} \approx \frac{\delta_p(1/f)}{(1/f)^2} \tag{003.Eq.3}$$

For `double`: $\delta_p(1/f) \approx 298.257 \times 2^{-52} \approx 6.6 \times 10^{-14}$, so $\delta_p(f) \approx 6.6\times 10^{-14} / 298.257^2 \approx 7.4 \times 10^{-19}$.

Additionally, the result of the division itself introduces a rounding error of at most one ULP of $f$: $\approx f \times \epsilon_T$.

**Accuracy error:** $\delta_a(f) = 0$ (no model involved)

**Verification:** $f = 0.003352810664747...$  Published: $3.3528106647475 \times 10^{-3}$. ✓

## Step 2: $e^2 = 2f - f^2$ — one multiply, one multiply, one subtract

Break this into sub-operations:

### Step 2a: $2f$

$$2f = 2 \times f \tag{003.Eq.4}$$

- $\sigma_m(2f) = 2 \times \sigma_m(f) = 0$
- $\delta_p(2f) = |2| \cdot \delta_p(f) + |f| \cdot \delta_p(2) + \delta_p(f)\delta_p(2) = 2\delta_p(f)$ (since $\delta_p(2) = 0$)
- $\delta_a(2f) = 0$

### Step 2b: $f^2$

$$f^2 = f \times f \tag{003.Eq.5}$$

- $\sigma_m(f^2) = 0$
- $\delta_p(f^2) = |f|\delta_p(f) + |f|\delta_p(f) + \delta_p(f)^2 = 2|f|\delta_p(f) + \delta_p(f)^2$
- $\delta_a(f^2) = 0$

### Step 2c: $e^2 = 2f - f^2$

Subtraction. Note: $2f \approx 6.706 \times 10^{-3}$ and $f^2 \approx 1.124 \times 10^{-5}$, so this is NOT a subtractive cancellation case — the terms differ by a factor of ~600.

- $\sigma_m(e^2) = \sigma_m(2f) + \sigma_m(f^2) = 0$
- $\delta_p(e^2) = \delta_p(2f) + \delta_p(f^2) = 2\delta_p(f) + 2|f|\delta_p(f) + \delta_p(f)^2$

Numerically for `double`: $\delta_p(f) \approx 10^{-18}$, so $\delta_p(e^2) \approx 2 \times 10^{-18} + 2 \times 3.35 \times 10^{-3} \times 10^{-18} + 10^{-36} \approx 2.007 \times 10^{-18}$.

Plus the rounding error of the subtraction itself.

- $\delta_a(e^2) = 0$

**Value:** $e^2 = 2(3.3528106647475 \times 10^{-3}) - (3.3528106647475 \times 10^{-3})^2$

$= 6.7056213294950 \times 10^{-3} - 1.1241336327082 \times 10^{-5}$

$= 6.69437999014 \times 10^{-3}$

**Published:** $6.694379990141 \times 10^{-3}$ ✓

**Note on alternative formula:** $e^2 = (a^2 - b^2)/a^2$ involves computing $a^2$ and $b^2$ which are large numbers (~$4.07 \times 10^{13}$), then subtracting two nearly equal large numbers. This DOES cause subtractive cancellation. The formula $e^2 = 2f - f^2$ avoids this entirely because $2f$ and $f^2$ differ by a factor of 600. **This is why Appendix B recommends Eq. B-4 over Eq. B-3.**

## Step 3: $e = \sqrt{e^2}$

- $\sigma_m(e) = 0$
- $\delta_p(e) = \delta_p(e^2) / (2\sqrt{e^2 - \delta_p(e^2)}) \approx \delta_p(e^2) / (2\sqrt{e^2})$
- $\delta_a(e) = 0$

**Value:** $e = \sqrt{0.00669437999014} = 0.0818191908426...$

**Published:** $8.1819190842622 \times 10^{-2}$ ✓

## Step 4: $e' = e / \sqrt{1 - e^2}$ — the second eccentricity

### Step 4a: $1 - e^2$

$$1 - e^2 = 1 - 0.006694380... = 0.993305620... \tag{003.Eq.6}$$

Not a cancellation problem (subtracting a small number from 1).

- All three errors: inherited from $e^2$, propagated through subtraction.

### Step 4b: $\sqrt{1 - e^2}$

$$\sqrt{0.993305620...} = 0.996647189... \tag{003.Eq.7}$$

- Errors from sqrt propagation.

### Step 4c: $e' = e / \sqrt{1-e^2}$

- $\sigma_m(e') = 0$ (all inputs defined)
- $\delta_p(e')$: accumulated from 4 operations (division of two square roots)
- $\delta_a(e') = 0$

**Value:** $e' = 0.08181919... / 0.99664719... = 0.08209443795...$

**Published:** $8.2094437949696 \times 10^{-2}$ ✓

**Alternative formula avoiding division:** $e'^2 = e^2/(1-e^2)$ or $e' = \sqrt{e^2/(1-e^2)}$. Equivalent.

## Step 5: $b = a(1 - f)$ — the semi-minor axis

$$b = a \times (1 - f) = 6378137.0 \times (1 - 0.003352810664748...) \tag{003.Eq.8}$$

$= 6378137.0 \times 0.996647189335...$

$= 6356752.31425...$

- $\sigma_m(b) = 0$ (both inputs defined)
- $\delta_p(b)$: from multiplication of $a$ (exact integer) by $(1-f)$ (has precision error from $f$)
  - $\delta_p(b) = |a| \cdot \delta_p(1-f) = a \cdot \delta_p(f)$ (since $\delta_p(1) = 0$ and subtraction propagates)
- $\delta_a(b) = 0$

**Published:** $6356752.3142$ m ✓

## Step 6: $E = ae$ — the linear eccentricity

$$E = 6378137.0 \times 0.0818191908426... = 521854.00842... \tag{003.Eq.9}$$

- $\sigma_m(E) = 0$
- $\delta_p(E) = |a| \cdot \delta_p(e)$ (since $a$ is an exact integer)
- $\delta_a(E) = 0$

**Published:** $5.2185400842339 \times 10^5$ m ✓

## Complete Dependency Graph (Geometric Constants)

```
Defining:  a (integer, exact)     1/f (decimal, representation error only)
           |                       |
           |                       v
           |                   f = 1/(1/f)
           |                   /         \
           |              2f              f²
           |               \             /
           |            e² = 2f - f²
           |           /        |       \
           |      e = √(e²)    |    1-e²
           |        /    \      |      |
           |   E = ae   e' = e/√(1-e²)
           |                    |
           v                    v
       b = a(1-f)        [feeds into q₀, q₀']
```

All measurement errors are zero. All accuracy errors are zero. Only precision errors propagate (representation + arithmetic rounding). The entire geometric constant set is derived from exactly two numbers: 6378137 and 298.257223563.
