# Style Guide

## Namespaces

Match directory structure. One namespace per conceptual domain.

```
math::           — pure mathematics, no physical knowledge
geodesy::        — equipotential ellipsoid and its properties
astronomy::      — solar system orbital elements
perturbation::   — orbital perturbation theory (Brouwer, Kaula, Hansen)
tle::            — two-line element set parsing
sgp4::           — SGP4 propagator assembly
```

## Naming

| Entity | Convention | Example |
|---|---|---|
| Namespace | `lowercase` | `math::`, `perturbation::` |
| Class | `PascalCase` | `EquipotentialEllipsoid`, `TrackedValue` |
| Function | `snake_case` | `solve_kepler`, `alternating_series` |
| Member variable | `snake_case` | `gamma_e`, `e_prime` |
| Template parameter | Single uppercase | `T` |
| Constant-returning function | `snake_case` | `exact(3)`, `ratio(1, 16)` |

## Type Conventions

Every public function accepts and returns `TrackedValue<T>`. No bare `T` in any public interface.

```cpp
// CORRECT:
TrackedValue<T> normal_gravity(const TrackedValue<T>& phi) const;

// WRONG:
T normal_gravity(T phi) const;
```

The only place bare `T` appears is:
- Inside `TrackedValue<T>` implementation (operating on `.value`)
- The `series_tolerance` parameter (a control parameter, not a physical value)
- Template type deduction contexts

## Function Signature Pattern

Every computational function follows this pattern:

```cpp
/// Brief description of what this computes.
/// Derivation: reference to design/derivations/NNN_name.md
/// Equation: reference to source document equation number
template<typename T>
TrackedValue<T> function_name(
    const TrackedValue<T>& input1,   // what this represents [units]
    const TrackedValue<T>& input2,   // what this represents [units]
    const T& tolerance = T(0))       // series tolerance (if applicable)
```

## Documentation Pattern

Every function, class, and module has a comment block stating:
1. What it computes (one sentence)
2. Where the formula comes from (derivation reference or equation number)
3. What the inputs represent (with units)
4. What errors apply (which of the three error types are relevant)

```cpp
/// Compute the normal gravity formula constant m = ω²a²b/GM.
/// Derivation: design/derivations/004_physical_constants.md, Step 1
/// Source: Moritz (1980) Eq. (2-70), NGA Appendix B Eq. (B-20)
/// Errors: measurement error propagates from GM.
```

## Return Values

Every computational function returns `TrackedValue<T>` carrying all three errors. No function returns bare `T`, `bool`, or error codes.

For functions that can fail (division by near-zero, non-convergence), the failure is encoded in the error bounds, not in exceptions or return codes:

```cpp
// Division by near-zero: the precision bound becomes very large,
// and reliable_digits() returns 0 or negative.
// The caller checks the error bounds, not a special return value.
auto result = a / b;
if (result.reliable_digits() < minimum_required) {
    // handle: switch to alternative formula, Lyddane modification, etc.
}
```

For functions returning multiple values, use a named struct:

```cpp
template<typename T>
struct BrouwerSecularRates {
    TrackedValue<T> M_dot;
    TrackedValue<T> omega_dot;
    TrackedValue<T> Omega_dot;
};
```

## Parameter Conventions

**Ordering:** Inputs describing the physical state come first, control parameters come last.

```cpp
TrackedValue<T> solve_kepler(
    const TrackedValue<T>& M,          // physical: mean anomaly
    const TrackedValue<T>& e,          // physical: eccentricity
    const T& tolerance,                // control: convergence tolerance
    int max_iterations = 50);          // control: iteration limit
```

**Const reference** for all `TrackedValue<T>` parameters. Never pass by value (they carry three extra `T` fields).

**Bare `T`** only for control parameters (tolerance, step size) that are not physical values and do not carry errors.

**Integer** for counts, indices, harmonic degree/order.

## Overflow, Underflow, and Singularity Handling

No exceptions. No error codes. The three-error framework IS the error handling.

**Overflow:** If a computation produces an infinite or NaN value, the `TrackedValue` stores it. The caller detects this via `std::isinf(result.value)` or `std::isnan(result.value)`. This should never happen in correctly formulated geodetic computations — if it does, it indicates a bug in the formula, not a runtime condition.

**Underflow:** If a value rounds to zero when it shouldn't be, the precision bound captures this. The `reliable_digits()` method returns 0, signaling the result is meaningless at this precision. The caller can increase precision (use wider `T`) and retry.

**Singularity (e.g., division by sin(i) near i=0):**

```cpp
auto sin_i = sin(inclination);
auto result = ph / sin_i;

// The division propagates sin_i's error into result's error.
// Near i=0, sin_i ≈ 0 but sin_i.errors.precision ≈ |cos(i)| * i.errors.precision ≈ i.errors.precision.
// So result.errors.precision ≈ |ph| * i_error / sin_i² which blows up.
// reliable_digits() detects this automatically.

if (result.reliable_digits() < 2) {
    // Switch to Lyddane modification — the MATH told us to switch,
    // not a hardcoded angular threshold.
}
```

**Non-convergence (Kepler solver, Newton iteration):** If the solver does not converge within `max_iterations`, it returns the best estimate with the last correction magnitude added to `precision`. The caller sees a large precision error and knows the result is unreliable.

**Negative square root argument:** If `sqrt` receives a value that could be negative within its error bounds (e.g., `1 - e²` where `e² ≈ 1`), the implementation computes `sqrt(max(value, 0))` and sets the precision bound to `sqrt(|original_error|)`. This preserves the error tracking without crashing.

## Module Interface Pattern

Each module exposes a small number of classes or functions that represent the module's conceptual output. Internal helpers are in anonymous namespaces or private methods.

```cpp
// geodesy/equipotential_ellipsoid.h — ONE class, all derived constants
// astronomy/solar_system.h — ONE class, all orbital elements
// perturbation/brouwer.h — functions computing secular rates
// perturbation/kaula.h — functions computing inclination functions
```
