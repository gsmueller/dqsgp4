# TrackedValue<T> — The Core Numeric Type

## Design Principle

`TrackedValue<T>` replaces raw `T` as the numeric type throughout the computation. Every constant, every intermediate result, and every output carries precision and accuracy bounds. Arithmetic operations propagate both error types automatically.

There is no separate `constants.h` that returns bare `T` values. Instead, constants are `TrackedValue<T>` instances with their error metadata baked in.

## Class Definition

```cpp
template<typename T>
class TrackedValue {
public:
    T value;
    T precision_bound;  // |computed_value - exact_formula_result|
    T accuracy_bound;   // |exact_formula_result - physical_truth|

    // --- Constructors ---

    // Exact mathematical constant (e.g., pi, sqrt(5))
    // Precision limited only by type T; accuracy is zero (mathematical truth)
    static TrackedValue Mathematical(const T& val, const T& type_precision);

    // Definitional physical constant (e.g., a = 6378137.0 m in WGS84)
    // Exact within the model; accuracy is zero by convention
    static TrackedValue Definitional(const T& val);

    // Measured physical constant (e.g., GM = 3.986004418e14)
    // Represented exactly in T; accuracy limited by measurement
    static TrackedValue Measured(const T& val, const T& measurement_uncertainty);

    // Series result
    static TrackedValue FromSeries(const T& val, const T& truncation_error,
                                    const T& inherited_accuracy);

    // --- Error queries ---

    T total_error() const { return precision_bound + accuracy_bound; }
    int reliable_digits() const;  // floor(-log10(total_error / |value|))

    // --- Arithmetic with error propagation ---

    friend TrackedValue operator+(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue operator-(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue operator*(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue operator/(const TrackedValue& a, const TrackedValue& b);
    friend TrackedValue operator-(const TrackedValue& a);  // unary negation

    // Math functions
    friend TrackedValue sqrt(const TrackedValue& a);
    friend TrackedValue abs(const TrackedValue& a);
    friend TrackedValue sin(const TrackedValue& a);
    friend TrackedValue cos(const TrackedValue& a);
    friend TrackedValue atan(const TrackedValue& a);
    friend TrackedValue pow(const TrackedValue& base, int exponent);
};
```

## Error Propagation Rules

### Addition / Subtraction

$$z = x \pm y$$

$$z.\text{precision} = x.\text{precision} + y.\text{precision}$$

$$z.\text{accuracy} = x.\text{accuracy} + y.\text{accuracy}$$

**Note on subtraction:** When $x \approx y$ (subtractive cancellation), the precision bound stays the same in absolute terms, but relative precision degrades because $|z|$ is small. The `reliable_digits()` method captures this:

```
reliable_digits = floor(-log10((precision + accuracy) / |value|))
```

### Multiplication

$$z = x \times y$$

$$z.\text{precision} = |x.\text{value}| \cdot y.\text{precision} + |y.\text{value}| \cdot x.\text{precision} + x.\text{precision} \cdot y.\text{precision}$$

$$z.\text{accuracy} = |x.\text{value}| \cdot y.\text{accuracy} + |y.\text{value}| \cdot x.\text{accuracy} + x.\text{accuracy} \cdot y.\text{accuracy}$$

The cross-term (precision × precision, accuracy × accuracy) is usually negligible but included for rigor. For small errors, the dominant terms are $|x|\delta_y + |y|\delta_x$.

### Division

$$z = x / y$$

$$z.\text{precision} = \frac{|x.\text{value}| \cdot y.\text{precision} + |y.\text{value}| \cdot x.\text{precision}}{y.\text{value}^2 - y.\text{precision} \cdot |y.\text{value}|}$$

$$z.\text{accuracy} = \frac{|x.\text{value}| \cdot y.\text{accuracy} + |y.\text{value}| \cdot x.\text{accuracy}}{y.\text{value}^2 - y.\text{accuracy} \cdot |y.\text{value}|}$$

(Conservative bound using the minimum possible denominator.)

### Square Root

$$z = \sqrt{x}$$

$$z.\text{precision} = \frac{x.\text{precision}}{2\sqrt{x.\text{value} - x.\text{precision}}}$$

$$z.\text{accuracy} = \frac{x.\text{accuracy}}{2\sqrt{x.\text{value} - x.\text{accuracy}}}$$

### Trigonometric Functions

$$z = \sin(x)$$

$$z.\text{precision} = |\cos(x.\text{value})| \cdot x.\text{precision}$$

$$z.\text{accuracy} = |\cos(x.\text{value})| \cdot x.\text{accuracy}$$

(Since $|d\sin/dx| = |\cos x| \leq 1$, trig functions never amplify errors.)

### Integer Powers (No Error Amplification from Exponent)

$$z = x^n$$

$$z.\text{precision} = |n| \cdot |x.\text{value}|^{n-1} \cdot x.\text{precision}$$

$$z.\text{accuracy} = |n| \cdot |x.\text{value}|^{n-1} \cdot x.\text{accuracy}$$

## Constants as TrackedValue

### Mathematical Constants

```cpp
template<typename T>
TrackedValue<T> Pi() {
    T val = boost::math::constants::pi<T>();
    // boost::math::constants provides pi to the full precision of T
    // Precision bound: the last representable ULP of T
    T prec = /* machine epsilon of T */ ;
    return TrackedValue<T>::Mathematical(val, prec);
}

template<typename T>
TrackedValue<T> Sqrt5() {
    T val = sqrt(T(5));
    T prec = /* propagated from T(5) through sqrt */ ;
    return TrackedValue<T>::Mathematical(val, prec);
}
```

### WGS84 Defining Parameters

```cpp
template<typename T>
struct WGS84Constants {
    static TrackedValue<T> a() {
        return TrackedValue<T>::Definitional(T("6378137.0"));
    }
    static TrackedValue<T> inv_f() {
        return TrackedValue<T>::Definitional(T("298.257223563"));
    }
    static TrackedValue<T> GM() {
        // Value from NGA.STND.0036 Eq. 3-4
        // Uncertainty ~±8×10⁶ m³/s² based on [NGA] discussion
        return TrackedValue<T>::Measured(T("3.986004418e14"), T("8e6"));
    }
    static TrackedValue<T> omega() {
        return TrackedValue<T>::Definitional(T("7.292115e-5"));
    }
};
```

### Derived Constants (Automatically Tracked)

```cpp
// f = 1 / (1/f) — both definitional, so f is definitional
auto f = TrackedValue<T>::Definitional(T(1)) / WGS84Constants<T>::inv_f();
// f.precision_bound ≈ machine epsilon (from division)
// f.accuracy_bound = 0 (definitional / definitional)

// e² = 2f - f² — definitional only
auto e2 = T(2) * f - f * f;
// e2.accuracy_bound = 0

// m = ω²a²b / GM — inherits GM accuracy
auto m = omega * omega * a * a * b / GM;
// m.accuracy_bound ≈ m.value * GM.accuracy_bound / GM.value
//                  ≈ 3.45e-3 * 8e6 / 3.986e14 ≈ 6.9e-14
```

## Impact on Class Hierarchy

```
math/
    tracked_value.h     — TrackedValue<T> with full error propagation
    series.h            — Series evaluation returning TrackedValue<T>
    binomial.h          — Generalized binomial returning TrackedValue<T>
    wallis.h            — Wallis integrals (exact rationals → TrackedValue with 0 precision error)
    vector3.h           — Vector3<TrackedValue<T>>

geodesy/
    equipotential_ellipsoid.h  — All members are TrackedValue<T>

sgp4/
    sgp4_propagator.h  — Output is Vector3<TrackedValue<T>> with full error budget
```

## Performance Consideration

`TrackedValue<T>` carries 3× the storage of bare `T` and each arithmetic operation does 3× the work (value + precision propagation + accuracy propagation). For `double`, this is negligible. For `cpp_bin_float_100`, the overhead is still only 3×, which is acceptable since the multiprecision arithmetic itself is already 100-1000× slower than `double`.

The error tracking can be disabled for performance-critical inner loops by providing a `TrackedValue<T, false>` specialization that skips bound propagation — but the default should always track.
