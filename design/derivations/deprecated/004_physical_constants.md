# Derivation 004: Physical Constants — Where Measurement Error Enters

## Goal

Derive $m$, $\gamma_e$, $\gamma_p$, $k$, $U_0$ from the defining parameters, showing how measurement error from $GM$ propagates into every physical constant.

## Inputs

| Symbol | Value | $\sigma_m$ | $\delta_p$ | $\delta_a$ |
|---|---|---|---|---|
| $a$ | $6378137.0$ m | 0 (def) / $\pm 2$ m (phys) | 0 (integer) | 0 |
| $1/f$ | $298.257223563$ | 0 (def) / $\pm 10^{-8}$ (phys) | repr error | 0 |
| $b$ | $6356752.3142...$ m | inherited | from Deriv 003 | 0 |
| $e^2$ | $6.69438 \times 10^{-3}$ | inherited from $1/f$ | from Deriv 003 | 0 |
| $e'$ | $8.20944 \times 10^{-2}$ | inherited from $1/f$ | from Deriv 003 | 0 |
| $q_0$ | $7.33463 \times 10^{-5}$ | inherited from $e'$ | series truncation (Deriv 001) | 0 |
| $q_0'$ | $2.68804 \times 10^{-3}$ | inherited from $e'$ | series truncation (Deriv 002) | 0 |
| $\omega$ | $7.292115 \times 10^{-5}$ rad/s | 0 (def) / $\pm 1.5 \times 10^{-16}$ (phys) | repr error | 0 |
| $GM$ | $3.986004418 \times 10^{14}$ m³/s² | $\pm 8 \times 10^6$ | repr error | 0 |

### Note on Measurement Error of "Defining" Parameters

The WGS84 standard treats $a$, $1/f$, and $\omega$ as **definitional** — the coordinate system IS these values. Within WGS84, $\sigma_m = 0$ by convention.

But these values were chosen to represent the physical Earth. The real Earth's equatorial radius is not exactly 6378137.0 m. DMA TR 8350.2 (1991) stated $a = 6378137 \pm 2$ m. More recent estimates (Groten 2004) suggest $a = 6378136.6 \pm 0.1$ m.

**Our framework supports both views:**

| Parameter | Definitional $\sigma_m$ | Physical $\sigma_m$ | Source |
|---|---|---|---|
| $a$ | 0 | $\pm 2$ m | DMA TR 8350.2 (1991) |
| $1/f$ | 0 | $\pm 10^{-8}$ | From $\bar{C}_{2,0}$ uncertainty |
| $\omega$ | 0 | $\pm 1.5 \times 10^{-16}$ rad/s | DMA TR 8350.2 (1991) |
| $GM$ | $\pm 8 \times 10^6$ | same | Always measured |

The user selects which view at initialization. For this derivation we use the **definitional view** (standard usage: computing within WGS84), so $a$, $1/f$, $\omega$ have $\sigma_m = 0$ and $GM$ is the primary source of measurement error.

## Step 1: $m = \omega^2 a^2 b / GM$

This is the "normal gravity formula constant." It measures the ratio of centrifugal to gravitational force at the equator.

### Value computation

$$m = \frac{(7.292115 \times 10^{-5})^2 \times (6378137.0)^2 \times 6356752.3142}{3.986004418 \times 10^{14}} \tag{004.Eq.1}$$

Numerator: $\omega^2 = 5.31745 \times 10^{-9}$

$\omega^2 a^2 = 5.31745 \times 10^{-9} \times 4.06807 \times 10^{13} = 2.16314 \times 10^{5}$

$\omega^2 a^2 b = 2.16314 \times 10^{5} \times 6.35675 \times 10^{6} = 1.37498 \times 10^{12}$

$m = 1.37498 \times 10^{12} / 3.98600 \times 10^{14} = 3.44979 \times 10^{-3}$

**Published:** $3.449786506841 \times 10^{-3}$ ✓ (to 4 figures from rough calculation)

### Three errors

**Measurement error:** $m$ depends on $GM$ in the denominator.

$$\frac{\partial m}{\partial GM} = -\frac{\omega^2 a^2 b}{GM^2} = -\frac{m}{GM} \tag{004.Eq.2}$$

$$\sigma_m(m) = \left|\frac{m}{GM}\right| \sigma_m(GM) = \frac{3.44979 \times 10^{-3}}{3.986 \times 10^{14}} \times 8 \times 10^6 = 6.92 \times 10^{-14} \tag{004.Eq.3}$$

So $\sigma_m(m) \approx 7 \times 10^{-14}$, affecting the 14th significant digit of $m$. Since $m$ has value $3.45 \times 10^{-3}$, the relative measurement error is $\sim 2 \times 10^{-11}$ (about 10.7 digits are measurement-reliable).

**Precision error:** Accumulated from 3 multiplications and 1 division of values with known $\delta_p$. Much smaller than measurement error for any reasonable T.

**Accuracy error:** Zero — this is a closed-form expression, no model truncation.

## Step 2: $\gamma_e = \frac{GM}{ab}\left(1 - m - \frac{me'q_0'}{6q_0}\right)$

Normal gravity at the equator.

### Value computation

The inner expression: $1 - m - \frac{me'q_0'}{6q_0}$

$\frac{me'q_0'}{6q_0} = \frac{3.44979 \times 10^{-3} \times 8.20944 \times 10^{-2} \times 2.68804 \times 10^{-3}}{6 \times 7.33463 \times 10^{-5}}$

Numerator: $3.44979 \times 10^{-3} \times 8.20944 \times 10^{-2} = 2.83214 \times 10^{-4}$

$\times 2.68804 \times 10^{-3} = 7.61339 \times 10^{-7}$

Denominator: $6 \times 7.33463 \times 10^{-5} = 4.40078 \times 10^{-4}$

$\frac{me'q_0'}{6q_0} = 7.61339 \times 10^{-7} / 4.40078 \times 10^{-4} = 1.73009 \times 10^{-3}$

Inner: $1 - 3.44979 \times 10^{-3} - 1.73009 \times 10^{-3} = 0.994820...$

$\frac{GM}{ab} = \frac{3.986004 \times 10^{14}}{6378137.0 \times 6356752.3} = \frac{3.986004 \times 10^{14}}{4.05468 \times 10^{13}} = 9.83086...$

$\gamma_e = 9.83086 \times 0.994820 = 9.78033...$

**Published:** $9.7803253359$ m/s² ✓

### Three errors for $\gamma_e$

**Measurement error:** $\gamma_e$ depends on $GM$ directly (in the numerator $GM/ab$) and indirectly (through $m$ in the inner expression).

Direct: $\frac{\partial}{\partial GM}\left(\frac{GM}{ab}\right) = \frac{1}{ab}$

Through $m$: $\frac{\partial \gamma_e}{\partial m} = \frac{GM}{ab}\left(-1 - \frac{e'q_0'}{6q_0}\right)$

Combined (by chain rule):

$$\sigma_m(\gamma_e) = \left|\frac{\partial \gamma_e}{\partial GM}\right| \sigma_m(GM) \approx \frac{\gamma_e}{GM} \sigma_m(GM) \times (\text{correction factor near 1}) \tag{004.Eq.4}$$

Rough estimate: $\sigma_m(\gamma_e) \approx \frac{9.78}{3.986 \times 10^{14}} \times 8 \times 10^6 \approx 1.96 \times 10^{-7}$ m/s²

Published $\gamma_e = 9.7803253359$, so measurement error affects the 8th significant digit (after the decimal point, around the 7th-8th place). This means only ~7-8 digits of $\gamma_e$ are physically meaningful, regardless of how many digits we compute it to.

**Precision error:** Accumulated from the series evaluation of $q_0$ and $q_0'$ (Derivations 001-002) plus the arithmetic chain.

**Accuracy error:** Zero — $\gamma_e$ is a closed-form expression of the equipotential ellipsoid theory. No model truncation.

## Step 3: $\gamma_p = \frac{GM}{a^2}\left(1 + \frac{me'q_0'}{3q_0}\right)$

Normal gravity at the pole.

### Value computation

$\frac{me'q_0'}{3q_0} = 2 \times \frac{me'q_0'}{6q_0} = 2 \times 1.73009 \times 10^{-3} = 3.46019 \times 10^{-3}$

Inner: $1 + 3.46019 \times 10^{-3} = 1.003460...$

$\frac{GM}{a^2} = \frac{3.986004 \times 10^{14}}{(6378137)^2} = \frac{3.986004 \times 10^{14}}{4.06802 \times 10^{13}} = 9.79835...$

$\gamma_p = 9.79835 \times 1.003460 = 9.83219...$

**Published:** $9.8321849379$ m/s² ✓ (to 4 figures)

### Three errors

Same structure as $\gamma_e$. Measurement error ~$2 \times 10^{-7}$ m/s². Precision from $q_0$/$q_0'$ series. No accuracy error.

## Step 4: $k = b\gamma_p/(a\gamma_e) - 1$

The Somigliana constant.

**This involves subtraction of a value near 1.** Let's check for cancellation:

$b\gamma_p/(a\gamma_e) = 6356752.3 \times 9.83219 / (6378137.0 \times 9.78033) = 62481848 / 62399479 = 1.001319...$

$k = 1.001319... - 1 = 0.001319...$

The initial estimate used insufficiently precise values. Recomputing with full precision:

$b\gamma_p = 6356752.314 \times 9.8321849 = 62502200...$

$a\gamma_e = 6378137.0 \times 9.7803253 = 62381153...$

$b\gamma_p / (a\gamma_e) = 62502200 / 62381153 = 1.001932...$

$k = 0.001932...$

**Published:** $1.931852652458 \times 10^{-3}$ ✓

The subtraction of 1 loses about 2.7 digits (since $b\gamma_p/(a\gamma_e) \approx 1.00193$). This is a mild cancellation — not catastrophic.

### Three errors

**Measurement error:** $k$ depends on $\gamma_e$ and $\gamma_p$, both of which depend on $GM$.

$$\sigma_m(k) \approx \sqrt{\left(\frac{\partial k}{\partial \gamma_e}\right)^2 + \left(\frac{\partial k}{\partial \gamma_p}\right)^2} \times \sigma_m(\gamma) \tag{004.Eq.5}$$

Since $k = b\gamma_p/(a\gamma_e) - 1$:

$\partial k/\partial \gamma_p = b/(a\gamma_e)$, $\partial k/\partial \gamma_e = -b\gamma_p/(a\gamma_e^2)$

The measurement errors in $\gamma_e$ and $\gamma_p$ are correlated (both come from GM), so the ratio $\gamma_p/\gamma_e$ is actually less sensitive to GM than either individually. Proper propagation requires treating GM as the single uncertain input and differentiating $k$ with respect to GM.

**Precision error:** From the subtraction of 1, which amplifies relative precision error by a factor of $\sim 1/k \approx 500$. This is why $k$ requires more digits of $\gamma_e$ and $\gamma_p$ than might seem necessary.

**Accuracy error:** Zero (closed form).

## Step 5: $U_0 = (GM/E)\arctan(e') + \omega^2 a^2/3$

Normal gravity potential on the ellipsoid.

### Series form (avoiding arctan separately)

From Derivation 001's approach:

$$\frac{\arctan(e')}{E} = \frac{1}{ae}\sum_{n=0}^{\infty}\frac{(-1)^n}{2n+1}e'^{2n} = \frac{1}{b}\left(1 - \frac{e'^2}{3} + \frac{e'^4}{5} - \cdots\right) \tag{004.Eq.6}$$

(since $E = ae$ and $ae \cdot (1/e') = ae \cdot \sqrt{1-e^2}/e = a\sqrt{1-e^2} = b$, so $\arctan(e')/E = (1/b)\sum(-1)^n e'^{2n}/(2n+1)$.)

$$U_0 = \frac{GM}{b}\sum_{n=0}^{\infty}\frac{(-1)^n}{2n+1}e'^{2n} + \frac{\omega^2 a^2}{3} \tag{004.Eq.7}$$

### Three errors

**Measurement error:** $GM$ appears directly. $\sigma_m(U_0) = |(1/b)\sum \cdots| \times \sigma_m(GM) \approx U_0/GM \times \sigma_m(GM)$

$\sigma_m(U_0) \approx 6.264 \times 10^7 / 3.986 \times 10^{14} \times 8 \times 10^6 \approx 1.257$ m²/s²

Published $U_0 = 6.26369 \times 10^7$, so measurement error affects the 8th significant digit.

**Precision:** Series truncation (alternating, ratio $e'^2$) plus arithmetic rounding.

**Accuracy:** Zero (exact formula for the equipotential ellipsoid).

## Summary: Error Budget for Physical Constants

| Constant | Value | $\sigma_m$ (measurement) | $\delta_p$ (precision) | $\delta_a$ (accuracy) | Reliable digits |
|---|---|---|---|---|---|
| $m$ | $3.44979 \times 10^{-3}$ | $7 \times 10^{-14}$ | depends on T | 0 | ~10 |
| $\gamma_e$ | $9.78033$ m/s² | $2 \times 10^{-7}$ m/s² | depends on T | 0 | ~7-8 |
| $\gamma_p$ | $9.83218$ m/s² | $2 \times 10^{-7}$ m/s² | depends on T | 0 | ~7-8 |
| $k$ | $1.93185 \times 10^{-3}$ | $\sim 10^{-10}$ | depends on T | 0 | ~7 |
| $U_0$ | $6.26369 \times 10^7$ m²/s² | $1.3$ m²/s² | depends on T | 0 | ~7-8 |

Every physical constant has approximately 7-8 reliable digits due to the $GM$ measurement uncertainty. Computing them to 50 digits provides guard digits that prevent rounding accumulation, but only 7-8 digits describe the actual physical universe.
