# The Three Fundamental Errors

Every value in our system carries three independent error quantities.

## 1. Measurement Error ($\sigma_m$)

The uncertainty in the physical measurement of an input constant. This is a property of the real world, not of our computation. It cannot be reduced by using more digits or better algorithms.

**Example:** $GM = 3.986004418 \times 10^{14} \pm 8 \times 10^{6}$ m³/s². The $\pm 8 \times 10^6$ is measurement error. No computation can reduce it.

**Propagation rule:** Standard sensitivity analysis. For $z = f(x, y)$:

$$\sigma_{m,z} = \sqrt{\left(\frac{\partial f}{\partial x}\right)^2 \sigma_{m,x}^2 + \left(\frac{\partial f}{\partial y}\right)^2 \sigma_{m,y}^2}$$

(assuming uncorrelated inputs; use covariance matrix if correlated)

**For defined constants** (a, 1/f, ω): $\sigma_m = 0$ exactly. These define the coordinate system.

## 2. Precision ($\delta_p$)

The error introduced by our computational representation and operations. This includes:
- Floating-point representation error (a decimal constant in binary)
- Series truncation (stopping after N terms)
- Rounding accumulation through arithmetic chains
- Iterative convergence tolerance (Kepler solver, Newton sqrt)

**Example:** Computing $q_0$ via the series to 8 terms gives a truncation error of $\sim 10^{-20}$. Using `double` instead of `cpp_bin_float_50` gives representation error of $\sim 10^{-16}$.

**Propagation rule:** Worst-case (additive) bounds. For $z = x \cdot y$:

$$\delta_{p,z} = |x|\delta_{p,y} + |y|\delta_{p,x} + \delta_{p,x}\delta_{p,y}$$

**This can always be reduced** by using more precision, more series terms, or tighter convergence tolerances.

## 3. Accuracy ($\delta_a$)

The error from using a simplified model instead of reality. This is a property of the mathematical model, not of the measurements or the computation.

**Example:** SGP4 uses only $J_2$, $J_3$, $J_4$ zonal harmonics. The true geopotential has terms through degree 2190 (EGM2008). The omitted terms cause ~100m orbit error for LEO satellites. No amount of precision in evaluating the $J_2$-only model eliminates this.

**Example:** The Brouwer secular rate formulas are truncated at $J_2^2$ and $J_4$. Higher-order terms ($J_2^3$, $J_2 J_3$, etc.) are omitted. This is a model accuracy limitation.

**Propagation rule:** Model-specific. For analytical models, the accuracy bound is typically estimated empirically or from the magnitude of the first omitted term in the perturbation expansion.

**This can only be reduced** by using a better model (more terms, numerical integration, etc.), NOT by using more computational precision.

## How They Combine

$$\text{Total uncertainty} = \delta_p + \sigma_m + \delta_a$$

(worst-case additive; or RSS for statistical combination)

But they must be tracked separately because they respond to different remedies:
- $\delta_p$ too large → use more digits, more series terms
- $\sigma_m$ too large → better measurements (not our problem)
- $\delta_a$ too large → better model (extend perturbation series, use numerical integration)

## Which Values Have Which Errors

| Value | $\sigma_m$ | $\delta_p$ | $\delta_a$ |
|---|---|---|---|
| $a = 6378137.0$ m | 0 (defined) | representation error in T | 0 (no model involved) |
| $GM = 3.986004418 \times 10^{14}$ | $\pm 8 \times 10^6$ | representation error in T | 0 (direct measurement) |
| $f = 1/(1/f)$ | 0 (from defined) | division rounding | 0 |
| $e^2 = 2f - f^2$ | 0 | subtraction + multiply rounding | 0 |
| $q_0$ (from series) | 0 (all inputs defined) | series truncation + rounding | 0 |
| $\gamma_e$ (from formula) | $\neq 0$ (depends on GM) | formula evaluation rounding | 0 (closed-form, no model truncation) |
| $m = \omega^2 a^2 b / GM$ | $\neq 0$ (depends on GM) | arithmetic rounding | 0 |
| $J_2$ (from $e^2$, $m$, $q_0$) | $\neq 0$ (depends on GM via $m$) | formula + series rounding | 0 |
| Brouwer $\dot{M}$ | $\neq 0$ (from $J_2$, $J_4$) | polynomial evaluation rounding | $\neq 0$ (truncated at $J_2^2$, $J_4$) |
| SGP4 position | $\neq 0$ (from TLE + gravity) | accumulated rounding | $\neq 0$ (analytical model ~1km) |

## Data Structure

```
struct ThreeErrors<T> {
    T measurement;  // σ_m: from physical measurements
    T precision;    // δ_p: from computational representation
    T accuracy;     // δ_a: from model truncation/simplification

    T total() const { return measurement + precision + accuracy; }
};

struct TrackedValue<T> {
    T value;
    ThreeErrors<T> errors;
};
```

## At Every Operation

Every arithmetic operation and function call must propagate all three independently:

```
TrackedValue<T> operator*(TrackedValue<T> a, TrackedValue<T> b) {
    return {
        .value = a.value * b.value,
        .errors = {
            // measurement: sensitivity analysis (RSS or additive)
            .measurement = |a.value| * b.errors.measurement
                         + |b.value| * a.errors.measurement,
            // precision: worst-case bound
            .precision = |a.value| * b.errors.precision
                       + |b.value| * a.errors.precision
                       + a.errors.precision * b.errors.precision,
            // accuracy: worst-case propagation
            .accuracy = |a.value| * b.errors.accuracy
                      + |b.value| * a.errors.accuracy
                      + a.errors.accuracy * b.errors.accuracy
        }
    };
}
```
