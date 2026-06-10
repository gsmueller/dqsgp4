# Error-Bounded Computation Design

## Core Principle

Every series-based computation returns both the computed value and a rigorous error bound. Functions iterate until the error bound is below a caller-specified tolerance, rather than using a fixed number of terms.

## 1. Function Prototype Pattern

```cpp
template<typename T>
struct BoundedResult {
    T value;          // computed value
    T error_bound;    // rigorous upper bound on |true_value - value|
    int terms_used;   // number of terms evaluated
};

template<typename T>
BoundedResult<T> ComputeR2(const EquipotentialEllipsoid<T>& ellipsoid, const T& tolerance);
```

Every function that evaluates a truncated series follows this pattern. The caller specifies the tolerance; the function returns when the error bound is proven to be below that tolerance.

## 2. Error Bounds for Alternating Series

Many geodetic series are **alternating** (signs alternate) with **monotonically decreasing** absolute term values. For such series, the **Leibniz criterion** gives a tight error bound:

If $S = \sum_{k=0}^{\infty} (-1)^k a_k$ where $a_k > 0$ and $a_k \downarrow 0$, and we truncate at $N$ terms:

$$S_N = \sum_{k=0}^{N-1} (-1)^k a_k$$

then

$$|S - S_N| \leq a_N$$

The error is bounded by the **first omitted term**.

### 2.1 Applicability

The following geodetic series are alternating with decreasing terms (for $e'^2 < 1$):

| Series | Alternating? | Decreasing? | Error bound |
|---|---|---|---|
| $Q$ (meridian quadrant) | Yes | Yes (for $e'^2 < 1$) | First omitted term |
| $R_2$ integral before sqrt | Yes | Yes | First omitted term |
| $q_0$ series form | Yes | Yes | First omitted term |
| $U_0$ series form | Yes | Yes | First omitted term |
| $J_{2n}$ harmonics | Yes | Yes | First omitted term |

### 2.2 Non-Alternating Series

The mean gravity $\bar{\gamma}$ involves both $e^2$ and $k$ terms (both positive), so it's not strictly alternating. For such series, we need:

**Geometric bound:** If $|r_{k+1}/r_k| \leq q < 1$ for all $k \geq N$ (where $r_k$ is the $k$-th term), then:

$$|S - S_N| \leq \frac{|r_N|}{1 - q}$$

For geodetic series with $e^2 \approx 0.0067$, the ratio $q$ is very small, making convergence rapid and bounds tight.

## 3. Error Bound for Composed Operations

### 3.1 Square Root of a Series ($R_2$)

$R_2 = c\sqrt{I}$ where $I = \sum c_k e'^{2k}$

If $I$ is known to within $\delta_I$ (i.e., $|I - I_N| \leq \delta_I$), then:

$$|R_2 - R_{2,N}| = c|\sqrt{I} - \sqrt{I_N}| \leq c \frac{\delta_I}{2\sqrt{I_N}}$$

(from the mean value theorem: $|\sqrt{x} - \sqrt{y}| \leq |x-y|/(2\sqrt{\min(x,y)})$)

### 3.2 Quotient of Two Series ($\bar{\gamma}$)

$\bar{\gamma}/\gamma_e = N/D$ where both $N$ and $D$ are known with bounds $\delta_N$ and $\delta_D$.

$$\left|\frac{N}{D} - \frac{N_n}{D_n}\right| \leq \frac{|N|\delta_D + |D|\delta_N}{D \cdot D_n}$$

Simplified (since $D \approx D_n \approx 1$): $\delta_{\bar{\gamma}} \lesssim \delta_N + \delta_D$

### 3.3 Iterative $e^2$ from $J_2$ (GRS80 path)

The iteration $e^2_{n+1} = 3J_2 + \frac{4}{15}\frac{\omega^2 a^3}{GM}\frac{a^3}{2q_0(e^2_n)}$ converges contractively.

Error bound: $|e^2 - e^2_n| \leq \frac{|e^2_n - e^2_{n-1}|}{1 - L}$ where $L < 1$ is the Lipschitz constant of the iteration map. For Earth-like parameters, $L \approx 0.003$, so convergence is very rapid.

## 4. Error Bound for $q_0$ and $q_0'$ (Numerically Sensitive)

The NGA document (Appendix B, Section 3.15) warns that $q_0$ and $q_0'$ are numerically sensitive due to subtractive cancellation in the $\arctan$ computation. The series forms avoid this:

$$2q_0 = \sum_{n=1}^{\infty} \frac{4(-1)^{n+1}n}{(2n+1)(2n+3)} e'^{2n+1}$$

This is an alternating series with decreasing terms, so the Leibniz bound applies directly. The series form is **preferred** for arbitrary precision because:

1. It avoids the subtractive cancellation in $(1+3/e'^2)\arctan(e') - 3/e'$
2. It provides a natural error bound (first omitted term)
3. It uses only rational coefficients and powers of $e'$

Similarly for $q_0'$:

$$q_0' = 3\left[\left(1+\frac{1}{e'^2}\right)\left(1-\frac{1}{e'}\arctan e'\right)\right] - 1$$

This should also be evaluated via its series expansion for arbitrary precision, with error bound from the omitted terms.

## 5. Convergence Acceleration (Optional)

For very high precision requirements, standard series may converge slowly. Options:

1. **Richardson extrapolation** — accelerate convergence of the partial sums
2. **Euler-Maclaurin summation** — convert tail of series to an integral
3. **Pade approximants** — rational function approximation of the partial sums

These are enhancements for later phases; the basic Leibniz/geometric bounds are sufficient for the initial implementation.

## 6. Implementation Pattern

```cpp
template<typename T>
BoundedResult<T> EvaluateAlternatingSeries(
    std::function<T(int)> term,  // term(k) returns the k-th term
    const T& tolerance,
    int max_terms = 1000)
{
    T sum = T(0);
    T prev_term;
    for (int k = 0; k < max_terms; ++k) {
        T t = term(k);
        sum += t;
        using std::abs;
        if (k > 0 && abs(t) < tolerance) {
            return {sum, abs(t), k + 1};
        }
        prev_term = t;
    }
    return {sum, abs(prev_term), max_terms};  // did not converge
}
```

## 7. Impact on Class Design

The `EquipotentialEllipsoid<T>` constructor takes a tolerance parameter:

```cpp
template<typename T>
class EquipotentialEllipsoid {
public:
    EquipotentialEllipsoid(T a, T inv_f, T GM, T omega, T tolerance);

    // All derived constants are BoundedResult<T>
    BoundedResult<T> q0() const;
    BoundedResult<T> q0_prime() const;
    BoundedResult<T> U0() const;
    BoundedResult<T> gamma_e() const;
    BoundedResult<T> gamma_p() const;
    BoundedResult<T> R2() const;
    BoundedResult<T> meridian_quadrant() const;
    BoundedResult<T> mean_gravity() const;

    // Series-evaluated at a point
    BoundedResult<T> normal_gravity(T phi) const;
    BoundedResult<T> J2n(int n) const;

    // Exact (no series needed)
    T b() const;       // b = a*sqrt(1-e²)
    T e2() const;      // e² = 2f - f²
    T f() const;       // f = 1/(1/f)
    T E_lin() const;   // E = a*e
    // ...
};
```

Note the distinction: quantities computed from closed-form expressions (involving only `sqrt`, `+`, `*`, `/`) are returned as plain `T`. Quantities requiring series evaluation return `BoundedResult<T>` with error bounds.
