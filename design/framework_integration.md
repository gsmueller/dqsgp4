# Framework Integration: TrackedValue Throughout the Entire Pipeline

## Principle

`TrackedValue<T>` is not an optional wrapper applied at the boundaries. It IS the numeric type used everywhere — from the initial constants through every intermediate computation to the final position/velocity output. Every function accepts `TrackedValue<T>` and returns `TrackedValue<T>`. There is no point in the pipeline where error information is discarded or unavailable.

## What This Means Concretely

### Every SGP4 Constant

Not just the WGS84 defining parameters, but every constant in the SGP4 theory — including the deep space constants that have no documentation in geodetic standards — must be categorized and tracked.

| Constant Category | Type | Example |
|---|---|---|
| WGS84 geometric definitions | `DefinedValue<T>` | $a$, $1/f$, $\omega$ |
| WGS84 measured constants | `MeasuredValue<T>` | $GM$, $G$, $M_A$, $H$ |
| SGP4 gravity model constants | `ModelValue<T>` | $J_2$, $J_3$, $J_4$, $\mu$, $r_e$ per model variant |
| SGP4 analytical theory constants | `ModelValue<T>` | ZNS, ZES, ZNL, ZEL, Q22, ROOT22, etc. |
| Solar/lunar ephemeris constants | `ModelValue<T>` | 4.5236020, 0.91375164, 5.8351514, etc. |
| Resonance polynomial coefficients | `ModelValue<T>` | g-function polynomials (3.616, -13.247, ...) |
| Mathematical constants | `MathConstant<T>` | $\pi$, $\sqrt{5}$, $2/3$ |
| Threshold/tolerance values | `Threshold<T>` | 1e-12 convergence, 0.2 rad Lyddane, 720 min step |

The new `ModelValue<T>` category captures constants from the SGP4 analytical theory that are neither measured nor defined — they are derived from a simplified physical model. Their accuracy bound represents the fidelity of that model.

### Every Intermediate Computation

The SGP4 initialization computes ~50 intermediate quantities (c1, c2, c3, c4, c5, d2, d3, d4, eta, coef, coef1, etc.). In the current implementations these are bare `double`. In our implementation, each is a `TrackedValue<T>` carrying its accumulated precision and accuracy bounds.

When the initialization computes:

```
eta = a * e * tsi
```

This is three multiplications of `TrackedValue<T>`. The precision bound accumulates through each multiply. The accuracy bound inherits from whichever inputs carry measurement or model uncertainty.

### Every Series Evaluation

The Kepler equation solver iterates:

$$E_{n+1} = E_n - \frac{E_n - e\sin E_n - M}{1 - e\cos E_n}$$

Each iteration operates on `TrackedValue<T>`. The convergence test compares the correction magnitude against the precision bound — the solver stops when the correction is smaller than the accumulated precision, meaning further iteration cannot improve the result given the input uncertainty.

### Every Trigonometric Call

When computing $\sin(i)$ where $i$ is a `TrackedValue<T>`:

- `sin(i).value` = computed sine
- `sin(i).precision_bound` = $|\cos(i)| \times i.\text{precision\_bound}$
- `sin(i).accuracy_bound` = $|\cos(i)| \times i.\text{accuracy\_bound}$

Near the inclination singularity ($i \approx 0$ or $\pi$), $|\sin(i)|$ is small but $|\cos(i)| \approx 1$, so the error in $\sin(i)$ is approximately the error in $i$ itself. The division $ph / \sin(i)$ then amplifies the error by $1/\sin(i)$ — the TrackedValue arithmetic reveals exactly when this amplification makes the result unreliable.

### The Singularity Detection Becomes Automatic

With TrackedValue propagation, we don't need hardcoded thresholds like `5.2359877e-2 rad` for the inclination singularity. Instead:

```cpp
TrackedValue<T> sin_i = sin(inclination);

// Division by sin(i) — TrackedValue propagation computes the error
TrackedValue<T> ph_over_sini = ph / sin_i;

// If the result's error exceeds its value, it's meaningless
if (ph_over_sini.total_error() > abs(ph_over_sini.value) * threshold_factor) {
    // Switch to Lyddane modification
}
```

The singularity is detected by the error bounds themselves, not by an arbitrary angular threshold. This adapts automatically to the precision level — at `double`, the switch happens around 3°; at 50-digit precision, it might happen at 0.001°.

### The Deep Space Integrator

The Euler-Maclaurin integration in dspace uses a 720-minute step:

$$x_{n+1} = x_n + \dot{x}\Delta t + \frac{1}{2}\ddot{x}\Delta t^2$$

Each step accumulates integration error. With TrackedValue, the precision bound grows with each step, and after enough steps, the precision bound exceeds the accuracy bound — telling the caller that computational error now dominates model error, and a smaller step size would help.

```cpp
// Integration step
TrackedValue<T> xli_new = xli + xldot * delt + xndot * step2;
// xli_new.precision_bound includes truncation error from the 2nd-order integrator
// This grows with each step — visible to the caller
```

### The Final Output

```cpp
struct PropagationResult {
    Vector3<TrackedValue<T>> position_km;
    Vector3<TrackedValue<T>> velocity_km_s;

    // The TrackedValue on each component tells the full story:
    // position_km.x.value            = -4400.594
    // position_km.x.precision_bound  = 1.2e-8 km  (computational)
    // position_km.x.accuracy_bound   = 0.8 km     (TLE + model)
    // position_km.x.total_error()    = 0.8 km      (accuracy-dominated)
    // position_km.x.reliable_digits()= 4
};
```

The caller sees immediately:
- The position is computed to ~8 decimal places of precision
- But only ~4 digits are physically meaningful (dominated by TLE accuracy and SGP4 model limitations)
- The extra precision serves as guard digits preventing rounding accumulation

## Accuracy Sources in the SGP4 Pipeline

### TLE Input Accuracy

TLE elements have limited precision from the format itself:

| Element | TLE Format Digits | Approximate Accuracy |
|---|---|---|
| Inclination | 8 digits (xxx.xxxx°) | ~0.0001° = 0.000002 rad |
| RAAN | 8 digits | ~0.0001° |
| Eccentricity | 7 digits (0.xxxxxxx) | ~1e-7 |
| Arg of perigee | 8 digits | ~0.0001° |
| Mean anomaly | 8 digits | ~0.0001° |
| Mean motion | 11 digits (xx.xxxxxxxx rev/day) | ~1e-8 rev/day |
| B* | 5 digits + exponent | ~10% relative |

These are `MeasuredValue<T>` inputs. Their accuracy bounds propagate through the entire computation.

### SGP4 Model Accuracy

The SGP4 analytical model itself introduces systematic error that cannot be reduced by precision:

| Effect | Approximate Error | Time Growth |
|---|---|---|
| Truncated geopotential (J2, J3, J4 only) | ~100 m | Secular |
| Simplified atmospheric drag (B*) | ~100 m for LEO | Quadratic in time |
| Analytical vs numerical integration | ~10 m near epoch | Linear |
| Simplified solar/lunar perturbations | ~10 m for deep space | Varies |
| Resonance model approximations | ~100 m for resonant orbits | Secular |

These are `ModelValue` accuracy contributions that affect the final output regardless of computational precision.

### How They Combine at Output

```
Total position error ≈ sqrt(
    TLE_precision_error²        (from MeasuredValue inputs)
  + SGP4_model_error²           (from ModelValue theory)
  + computational_precision²    (from TrackedValue accumulation)
  + epoch_staleness²            (grows with time from TLE epoch)
)
```

The `TrackedValue` framework tracks the first and third terms automatically. The second term is a model property attached to the SGP4 constants as `ModelValue` accuracy bounds. The fourth term is computed from the time since epoch and the secular error growth rates.

## New Type: ModelValue

```cpp
template<typename T>
class ModelValue {
    T value_;
    T model_uncertainty_;    // accuracy of the analytical model
    std::string provenance_; // "Spacetrack Report #3" etc.
public:
    ModelValue(const char* value, const char* uncertainty, const char* source);
    const T& value() const { return value_; }
    const T& model_uncertainty() const { return model_uncertainty_; }
    const std::string& provenance() const { return provenance_; }
};
```

Example:
```cpp
auto ZNS = ModelValue<T>("1.19459e-5", "1e-10",
    "Solar node rate, Spacetrack Report #3 (1980)");

auto J2_wgs72 = ModelValue<T>("0.001082616", "1e-9",
    "WGS72 J2, Spacetrack Report #3");
```

When a `ModelValue` participates in arithmetic with other `TrackedValue` types, its `model_uncertainty` flows into the `accuracy_bound` of the result.

## Type Hierarchy Summary

```
                    TrackedValue<T>
                   /       |       \
          DefinedValue  MeasuredValue  DerivedValue
              |              |              |
         (a, 1/f, ω)    (GM, G, H)    (all computations)
                              |
                        ModelValue
                           |
                   (J2, ZNS, Q22, ...)
                   (SGP4 theory constants)
```

All types convert to `DerivedValue` when they participate in arithmetic. The `DerivedValue` carries the accumulated precision and accuracy bounds from all inputs.
