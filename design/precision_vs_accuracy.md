# Precision vs. Accuracy Tracking

## Definitions

**Precision** — The number of meaningful digits in our numerical representation. Determined by:
- The floating-point type `T` (e.g., `double` ≈ 15.9 digits, `cpp_bin_float_50` = 50 digits)
- Series truncation error (how many terms we evaluated)
- Numerical cancellation (subtractive loss of significant figures)

**Accuracy** — How close our computed value is to the true physical quantity. Determined by:
- Measurement uncertainty in the input constants
- Model fidelity (e.g., SGP4 is an analytical approximation of the true orbit)
- Omitted physical effects (atmospheric drag uncertainty, unmodeled perturbations)

**Key insight:** We can compute $q_0$ to 100 digits of precision, but if $e'$ is derived from $1/f = 298.257223563$ (12 significant figures), then $q_0$ is only *accurate* to about 12 figures regardless of how many digits we carry. However, carrying extra precision prevents *accumulation* of rounding errors through long computation chains — the extra digits serve as guard digits.

## The Error Model

Every quantity carries three pieces of information:

```cpp
template<typename T>
struct TrackedValue {
    T value;                // the computed value
    T precision_bound;      // upper bound on |value - exact_computation|
                            // (series truncation + numerical noise)
    T accuracy_bound;       // upper bound on |exact_computation - physical_truth|
                            // (input uncertainty + model error)
};
```

The **total error bound** at output is:

$$|\text{output} - \text{physical truth}| \leq \text{precision\_bound} + \text{accuracy\_bound}$$

But they propagate by different rules and have different sources, so they must be tracked separately.

## Sources of Each Error Type

### Precision Sources (Computational)

| Source | Affects | Mitigation |
|---|---|---|
| Type `T` representation | All values | Use sufficient precision type |
| Series truncation | $q_0$, $q_0'$, $R_2$, $\bar{\gamma}$, $Q$, normal gravity | `BoundedResult<T>` with tolerance |
| Subtractive cancellation | $q_0$ closed form, $e^2 = 2f - f^2$ | Use series form; use $(2f-f^2)$ not $(a^2-b^2)/a^2$ |
| Iterative convergence | $e^2$ from $J_2$ (GRS80 path) | Track $|e^2_n - e^2_{n-1}|$ |
| Accumulated rounding | Long computation chains | Carry guard digits |

### Accuracy Sources (Physical)

| Source | Affects | Magnitude | How Known |
|---|---|---|---|
| Defining parameter measurement | All derived constants | Varies by parameter | Published uncertainties |
| $a = 6378137.0$ m | Exact by definition within WGS84 | 0 (definitional) | — |
| $1/f = 298.257223563$ | Exact by definition within WGS84 | 0 (definitional) | — |
| $GM = 3.986004418 \times 10^{14}$ | Defines the mass scale | $\pm 0.000000008 \times 10^{14}$ (~1 ppb) | [NGA] Ch. 3 |
| $\omega = 7.292115 \times 10^{-5}$ | Exact by convention | 0 (conventional) | — |
| TLE element precision | SGP4 orbital elements | ~$10^{-8}$ rad in angles | TLE format (8 digits) |
| SGP4 model truncation | Propagated position/velocity | ~1 km position, ~1 m/s velocity | Empirical |
| Atmospheric drag model | LEO propagation | Dominant error source for LEO | Varies by solar activity |
| Unmodeled perturbations | Deep space propagation | Tesseral harmonics, SRP, etc. | Varies by orbit |

### Critical Distinction: Definitional vs. Measured

Within the WGS84 system, $a$ and $1/f$ are **definitional** — they define the coordinate system, not measurements of the Earth. Their accuracy within WGS84 is exact (zero uncertainty). But $GM$ is a **measured** quantity with real uncertainty, and the TLE orbital elements are measurements with limited precision.

## Propagation Rules

### Precision Propagation (Standard Error Propagation)

For $z = f(x, y)$ where $x$ and $y$ have precision bounds $\delta_x$ and $\delta_y$:

$$\delta_z \leq \left|\frac{\partial f}{\partial x}\right| \delta_x + \left|\frac{\partial f}{\partial y}\right| \delta_y$$

For common operations:
- Addition: $\delta_{x+y} = \delta_x + \delta_y$
- Multiplication: $\delta_{xy} \leq |x|\delta_y + |y|\delta_x$
- Division: $\delta_{x/y} \leq (|x|\delta_y + |y|\delta_x)/y^2$
- Square root: $\delta_{\sqrt{x}} \leq \delta_x / (2\sqrt{x})$
- Series composition: track through `BoundedResult`

### Accuracy Propagation (Sensitivity Analysis)

For a derived constant $D = f(a, f, GM, \omega)$:

$$\sigma_D^2 = \left(\frac{\partial D}{\partial a}\right)^2 \sigma_a^2 + \left(\frac{\partial D}{\partial f}\right)^2 \sigma_f^2 + \left(\frac{\partial D}{\partial GM}\right)^2 \sigma_{GM}^2 + \left(\frac{\partial D}{\partial \omega}\right)^2 \sigma_\omega^2$$

Since $a$, $f$, and $\omega$ are definitional (zero uncertainty within WGS84), the accuracy of derived constants depends primarily on $\sigma_{GM}$.

For SGP4 propagation, the accuracy is dominated by:
1. TLE element precision (input)
2. SGP4 model fidelity (systematic)
3. Atmospheric drag uncertainty (for LEO)
4. Time since epoch (error grows)

## Classification of Quantities

### Tier 1: Exact by Definition (precision-limited only)

These have zero accuracy error within the WGS84 model. Their error is purely precision (rounding, truncation).

| Quantity | Source |
|---|---|
| $a$, $1/f$, $\omega$ | WGS84 defining parameters |
| $f$, $b$, $e^2$, $e'^2$, $e$, $e'$, $E$, $c$, $AR$ | Derived from $a$ and $f$ only |
| $R_1 = a(1-f/3)$ | Derived from $a$ and $f$ only |
| $R_3 = \sqrt[3]{a^2 b}$ | Derived from $a$ and $f$ only |

### Tier 2: Dependent on GM (accuracy limited by $\sigma_{GM}$)

| Quantity | Additional input |
|---|---|
| $m = \omega^2 a^2 b / GM$ | $GM$ |
| $\gamma_e$, $\gamma_p$, $k$, $f^*$ | $GM$ via $m$ |
| $U_0$ | $GM$ directly |
| $\bar{\gamma}$ | $GM$ via $\gamma_e$ and $k$ |
| $J_2$, $\bar{C}_{2,0}$ | $GM$ via $m$ |
| $M = GM/G$ | Both $GM$ and $G$ |

### Tier 3: Series-Evaluated (precision limited by truncation + accuracy limited by inputs)

| Quantity | Precision source | Accuracy source |
|---|---|---|
| $q_0$, $q_0'$ | Series truncation | $e'$ (Tier 1) |
| $R_2$ | Series truncation + sqrt | $e'$ (Tier 1) |
| $Q$ | Series truncation | $e'$ (Tier 1) |
| $\bar{\gamma}$ | Series truncation + quotient | $\gamma_e$, $k$ (Tier 2) |

### Tier 4: SGP4 Propagation Output (accuracy dominated by model)

| Quantity | Precision source | Accuracy source |
|---|---|---|
| Position (km) | Accumulated rounding | SGP4 model (~1 km), TLE precision, drag |
| Velocity (km/s) | Accumulated rounding | SGP4 model (~1 m/s), TLE precision, drag |

## Implementation Design

### TrackedValue Class

```cpp
template<typename T>
class TrackedValue {
public:
    T value;
    T precision_bound;  // |computed - exact_formula_result|
    T accuracy_bound;   // |exact_formula_result - physical_truth|

    // Total error: precision_bound + accuracy_bound
    T total_error() const { return precision_bound + accuracy_bound; }

    // The number of reliable digits
    int reliable_digits() const;

    // Arithmetic operations propagate both bounds
    friend TrackedValue operator+(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue operator*(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue operator/(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue sqrt(const TrackedValue& a);
};
```

### Factory Methods for Defining Parameters

```cpp
// Definitional constant: zero accuracy error, zero precision error (exact representation)
template<typename T>
TrackedValue<T> Definitional(const T& value) {
    return {value, T(0), T(0)};
}

// Measured constant: zero precision error, nonzero accuracy error
template<typename T>
TrackedValue<T> Measured(const T& value, const T& uncertainty) {
    return {value, T(0), uncertainty};
}

// Series result: nonzero precision error, inherits accuracy from inputs
template<typename T>
TrackedValue<T> FromSeries(const T& value, const T& truncation_error,
                           const T& inherited_accuracy) {
    return {value, truncation_error, inherited_accuracy};
}
```

### Initialization Example

```cpp
auto a     = Definitional<T>(T("6378137.0"));
auto inv_f = Definitional<T>(T("298.257223563"));
auto GM    = Measured<T>(T("3.986004418e14"), T("8e6"));  // ±8×10⁶ m³/s²
auto omega = Definitional<T>(T("7.292115e-5"));

// f is exact (derived from definitional 1/f)
auto f = Definitional<T>(T(1) / inv_f.value);

// m depends on GM, so it inherits GM's accuracy
auto m = omega * omega * a * a * b / GM;
// m.precision_bound ≈ 0 (exact arithmetic on exact inputs + one measured)
// m.accuracy_bound  ≈ |∂m/∂GM| × σ_GM = m/GM × σ_GM
```

## Output Error Estimate

At the SGP4 propagation output, the error estimate combines:

1. **Input precision** — how well we represented the TLE elements
2. **Computational precision** — accumulated rounding through the propagator
3. **TLE accuracy** — the TLE elements themselves have limited accuracy (~8 decimal digits)
4. **Model accuracy** — SGP4's analytical approximation vs. true orbit

```cpp
struct PropagationResult<T> {
    Vector3<TrackedValue<T>> position;   // km, with error tracking
    Vector3<TrackedValue<T>> velocity;   // km/s, with error tracking

    // Convenience: total position error estimate (RSS of 3 components)
    T position_error_estimate() const;
    T velocity_error_estimate() const;

    // Breakdown
    T position_precision() const;  // computational contribution
    T position_accuracy() const;   // physical/model contribution
};
```

## Relationship to the Existing `BoundedResult<T>`

`BoundedResult<T>` tracks **precision only** (series truncation). `TrackedValue<T>` extends this to also track **accuracy** (physical uncertainty). The two compose:

```cpp
// Series evaluation returns BoundedResult (precision only)
BoundedResult<T> q0_result = EvaluateQ0Series(e_prime, tolerance);

// Combine with accuracy inherited from inputs
TrackedValue<T> q0 = FromSeries(
    q0_result.value,
    q0_result.error_bound,           // precision
    /* accuracy from e' */ T(0)       // e' is Tier 1 (definitional), so 0
);
```

For Tier 2 quantities where accuracy depends on $GM$:

```cpp
TrackedValue<T> gamma_e = ComputeGammaE(a, b, GM, m, q0, q0_prime);
// gamma_e.precision_bound = f(truncation errors in q0, q0')
// gamma_e.accuracy_bound  = |∂γₑ/∂GM| × σ_GM
```
