# Defined Values vs. Derived Values

## The Fundamental Distinction

A **defined value** is not a measurement, not an approximation, not a computation. It is an exact assertion. It has no error — not small error, not negligible error — **zero error**. The digits are not "significant figures" subject to rounding; they are the complete, exact specification.

A **derived value** is the result of applying operations to other values. Even if all inputs are defined, the result may have precision error from the computation. If any input is measured, the result inherits accuracy error.

These are categorically different things and should be different types in the code.

## Categories

### 1. Defined Constants

These are exact by decree. They define the reference system. Changing even one digit would create a different reference system. Examples:

| Value | Why It's Defined |
|---|---|
| $a = 6378137.0$ m | Defines the WGS84 ellipsoid |
| $1/f = 298.257223563$ | Defines the WGS84 ellipsoid |
| $\omega = 7.292115 \times 10^{-5}$ rad/s | Defines the WGS84 rotation rate |
| $c = 299792458$ m/s | Defines the meter via the second |
| $GM_\text{GPSNAV} = 3.9860050 \times 10^{14}$ | Defines the GPS navigation interface |

**Key property:** The defined value has zero physical uncertainty but may have nonzero representation error in type `T`.

- **Integers** (e.g., $a = 6378137$): Exactly representable in binary floating point as long as they fit in the mantissa. `double` can represent integers up to $2^{53}$ exactly. $6378137 < 2^{53}$, so `double` represents it exactly. So does any wider type.

- **Terminating decimals that are NOT powers of 2** (e.g., $1/f = 298.257223563$): The decimal digits are exact by definition, but this number has an infinite repeating binary expansion. In `double` (~15.9 decimal digits), the representation error is $\sim 10^{-13}$. In `cpp_bin_float_50` (50 decimal digits), the string constructor captures all 12 source digits exactly, but if we then divide (e.g., $f = 1/(1/f)$), the result is an infinite non-repeating decimal that truncates at 50 digits.

- **Irrational results from defined inputs** (e.g., $e = \sqrt{2f - f^2}$): The sqrt of a rational is irrational. No finite representation is exact. The representation error depends entirely on the precision of `T`.

The `DefinedValue` type must track this representation error so it correctly propagates through subsequent arithmetic.

### 2. Measured Constants

These represent our best knowledge of a physical quantity, subject to measurement uncertainty. The uncertainty is a property of the measurement, not of our computation. Examples:

| Value | Uncertainty | Source |
|---|---|---|
| $GM = 3.986004418 \times 10^{14}$ m$^3$/s$^2$ | $\pm 8 \times 10^{6}$ | Satellite tracking, [NGA] |
| $G = 6.67428 \times 10^{-11}$ m$^3$/(kg s$^2$) | $\pm 0.00067 \times 10^{-11}$ | Laboratory experiments |
| $H = 3.273795 \times 10^{-3}$ | $\pm 0.000001 \times 10^{-3}$ | Precession observations |
| $M_A = 5.1480 \times 10^{18}$ kg | $\pm 0.0003 \times 10^{18}$ | Atmospheric models |

**Key property:** No amount of computational precision changes the measurement uncertainty. Computing $GM/G$ to 100 digits doesn't make $M$ accurate to 100 digits.

### 3. Derived Exact Values (from defined inputs only)

When ALL inputs are defined constants, the result is in principle exact — the only error is computational precision (rounding, truncation). Examples:

| Value | Formula | Inputs |
|---|---|---|
| $f = 1/(1/f)$ | Division | Defined |
| $b = a(1-f)$ | Multiply, subtract | Defined |
| $e^2 = 2f - f^2$ | Multiply, subtract | Defined |
| $e' = e/\sqrt{1-e^2}$ | Divide, sqrt | Defined |
| $E = ae$ | Multiply | Defined |
| $R_1 = a(1 - f/3)$ | Multiply, subtract | Defined |
| $R_3 = \sqrt[3]{a^2 b}$ | Multiply, cbrt | Defined |
| $q_0$ | Series in $e'$ | Defined (but series truncation adds precision error) |
| $R_2$ | Series in $e'$ + sqrt | Defined (but series truncation adds precision error) |

**Key property:** These have zero accuracy error. The only error is precision — how well our computation approximates the exact mathematical result. This can be driven to zero by using more precision and more series terms.

### 4. Derived Measured Values (at least one measured input)

When ANY input is measured, the result inherits measurement uncertainty that cannot be eliminated. Examples:

| Value | Formula | Measured input |
|---|---|---|
| $m = \omega^2 a^2 b / GM$ | Arithmetic | $GM$ |
| $\gamma_e$ | Closed form involving $m$, $q_0$, $q_0'$ | $GM$ via $m$ |
| $\gamma_p$ | Closed form involving $m$, $q_0$, $q_0'$ | $GM$ via $m$ |
| $k = b\gamma_p/(a\gamma_e) - 1$ | Arithmetic | $GM$ via $\gamma_e$, $\gamma_p$ |
| $U_0 = (GM/E)\arctan e' + \omega^2 a^2/3$ | Arithmetic + transcendental | $GM$ directly |
| $J_2$ (when derived from $f$ and $m$) | Formula involving $m$ | $GM$ via $m$ |
| $M = GM/G$ | Division | Both $GM$ and $G$ |

### 5. Model-Limited Values

The propagation output. Even with perfect inputs and infinite precision, the SGP4 analytical model doesn't match the true orbit. This is a different kind of error from measurement uncertainty — it's systematic model deficiency.

| Value | Model limitation |
|---|---|
| SGP4 position | ~1 km for near-epoch, grows with time |
| SGP4 velocity | ~1 m/s for near-epoch, grows with time |

## Type System Design

```cpp
// Defined value: exact by decree, but may not be exactly representable in T
template<typename T>
class DefinedValue {
    T value_;                    // best representation in type T
    T representation_error_;     // |T_value - exact_defined_value|, >= 0
    bool exactly_representable_; // true if representation_error_ == 0

public:
    // Integer-valued defined constant (always exactly representable)
    explicit DefinedValue(int exact_integer);

    // Decimal-valued defined constant
    // The constructor computes the representation error by:
    //   1. Parsing the decimal string into T
    //   2. If T has enough precision to hold all decimal digits: exact (error = 0)
    //   3. If T has fewer digits than the source string: error = T's ULP at that magnitude
    explicit DefinedValue(const char* exact_decimal, int source_digits);

    // Rational-valued defined constant (e.g., 1/298.257223563)
    // Cannot be exactly represented in ANY finite binary/decimal float
    // The representation error is bounded by T's epsilon at that magnitude
    static DefinedValue FromRatio(const T& numerator, const T& denominator);

    const T& value() const { return value_; }
    const T& representation_error() const { return representation_error_; }
    bool is_exact() const { return exactly_representable_; }
};

// Measured value: carries measurement uncertainty
template<typename T>
class MeasuredValue {
    T value_;
    T uncertainty_;  // 1-sigma measurement uncertainty
public:
    MeasuredValue(const char* best_estimate, const char* uncertainty);
    const T& value() const { return value_; }
    const T& uncertainty() const { return uncertainty_; }
};

// Derived value: carries both precision error and (possibly) accuracy error
template<typename T>
class DerivedValue {
    T value_;
    T precision_bound_;   // from computation (rounding, series truncation)
    T accuracy_bound_;    // inherited from measured inputs (0 if all inputs defined)
public:
    const T& value() const { return value_; }
    const T& precision_bound() const { return precision_bound_; }
    const T& accuracy_bound() const { return accuracy_bound_; }
    T total_error() const { return precision_bound_ + accuracy_bound_; }
    bool is_accuracy_limited() const { return accuracy_bound_ > precision_bound_; }
    bool is_precision_limited() const { return precision_bound_ > accuracy_bound_; }
};
```

## Representation Format Considerations

### Binary vs. Decimal Floating Point

| Type | Can represent `298.257223563` exactly? | Can represent `1/3` exactly? | Can represent `6378137` exactly? |
|---|---|---|---|
| `double` (binary64) | No | No | Yes (integer < $2^{53}$) |
| `cpp_bin_float_50` (binary) | No (infinite repeating binary) | No | Yes |
| `cpp_dec_float_50` (decimal) | Yes (12 decimal digits < 50) | No (infinite repeating decimal) | Yes |

**Implication:** If we use binary floating point (`cpp_bin_float`), then defined decimal values like $1/f = 298.257223563$ have a nonzero representation error at the last bit, regardless of how many bits we use. With decimal floating point (`cpp_dec_float`), the same value is exact — but then `1/3` and `sqrt(2)` have representation errors instead.

**Design choice:** Since the WGS84 defining parameters are specified as terminating decimals, `cpp_dec_float` would represent them exactly. However, `cpp_bin_float` is faster for transcendental functions (sin, cos, arctan) which dominate the SGP4 computation. The representation error for defined constants in `cpp_bin_float_50` is $\sim 10^{-50}$, which is far below any physical significance.

The `DefinedValue` class must compute the actual representation error for whichever type `T` is used, not assume it's zero.

### Exact Representability Test

A decimal value $d$ is exactly representable in binary floating point with $p$ bits of mantissa if and only if $d \times 2^k$ is an integer for some integer $k$ where $|d \times 2^k| < 2^p$. Most geodetic constants fail this test.

For `DefinedValue`, the representation error can be bounded by:

$$\text{representation\_error} \leq \frac{|d|}{2} \times 2^{1 - p_\text{effective}}$$

where $p_\text{effective}$ is the number of mantissa bits (or decimal digits × $\log_2 10$ for decimal types).

## Arithmetic Rules by Type Interaction

### Defined ⊕ Defined → DerivedValue (precision error only)

```
DefinedValue + DefinedValue → DerivedValue(precision = rounding, accuracy = 0)
DefinedValue × DefinedValue → DerivedValue(precision = rounding, accuracy = 0)
```

The result is no longer "defined" — it's derived. Even though accuracy is zero, there may be precision error from the arithmetic. For exact operations on exact representable values (e.g., integer × integer), precision error is also zero.

### Defined ⊕ Measured → DerivedValue (accuracy inherited)

```
DefinedValue × MeasuredValue → DerivedValue(
    precision = |defined| × measured.representation_rounding,
    accuracy  = |defined| × measured.uncertainty
)
```

### Derived ⊕ Derived → DerivedValue (both propagate)

Standard error propagation as documented in `tracked_value_class.md`.

### Measured ⊕ Measured → DerivedValue (uncertainties combine)

```
GM / G → DerivedValue(
    precision = rounding,
    accuracy  = sqrt((σ_GM/GM)² + (σ_G/G)²) × |GM/G|
)
```

## What This Means for the Codebase

### The EquipotentialEllipsoid constructor

```cpp
template<typename T>
class EquipotentialEllipsoid {
public:
    // Inputs are explicitly typed as Defined or Measured
    EquipotentialEllipsoid(
        DefinedValue<T> a,
        DefinedValue<T> inv_f,
        MeasuredValue<T> GM,
        DefinedValue<T> omega,
        T series_tolerance
    );

    // Tier 1: derived from defined only → precision error only, accuracy = 0
    DerivedValue<T> f() const;
    DerivedValue<T> b() const;
    DerivedValue<T> e2() const;
    DerivedValue<T> e_prime() const;
    DerivedValue<T> E_lin() const;
    DerivedValue<T> q0() const;       // series truncation → nonzero precision
    DerivedValue<T> R2() const;       // series truncation → nonzero precision

    // Tier 2: derived from defined + measured → has accuracy error
    DerivedValue<T> m() const;        // depends on GM
    DerivedValue<T> gamma_e() const;  // depends on GM
    DerivedValue<T> gamma_p() const;  // depends on GM
    DerivedValue<T> U0() const;       // depends on GM
    DerivedValue<T> k() const;        // depends on GM
    DerivedValue<T> J2() const;       // depends on GM
};
```

### Why This Matters

When someone computes $\gamma_e$ to 50 decimal places using `cpp_bin_float_50`, the result will report:

```
gamma_e:
  value            = 9.78032533590000000000000000000000000000000000000000 m/s²
  precision_bound  = 2.3e-50 (series truncation + rounding)
  accuracy_bound   = 1.9e-08 (inherited from GM uncertainty)
  reliable digits  = 8 (accuracy-limited, not precision-limited)
```

This tells the user: "We computed this to 50 digits, but only 8 are physically meaningful because $GM$ is only known to ~10 significant figures." The extra precision serves as guard digits that prevent rounding accumulation.

For a Tier 1 quantity like $e^2$:

```
e2:
  value            = 0.00669437999014131699613723340...  (50 digits)
  precision_bound  = 1.2e-51 (rounding only)
  accuracy_bound   = 0 (all inputs are defined)
  reliable digits  = 50 (precision-limited only)
```

This tells the user: "All 50 digits are meaningful because this is derived purely from defined constants."
