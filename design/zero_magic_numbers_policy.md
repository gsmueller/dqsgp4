# Zero Magic Numbers Policy

## Principle

There shall be no unexplained numeric literal in any source file. Every constant must trace to one of:

1. **A defining parameter** — identified by name, source document, and equation number
2. **A mathematical identity** — derivable from first principles with explicit steps shown
3. **A computed result** — produced by a function that implements a documented formula
4. **A design choice** — explicitly justified with rationale (e.g., integration step size, convergence tolerance)

"Magic number" means any numeric literal whose origin is not immediately traceable. **Magic numbers are bugs.**

## Implementation Rules

### Rule 1: Every Constant Has a Name and a Source

```cpp
// WRONG:
double temp = 1.5 * 0.001082616 * pinvsq * n;

// RIGHT:
const auto three_halves = T(3) / T(2);
const auto& j2 = gravity_model.j2();  // Source: [NGA] Table 3.1 / getgravconst()
auto temp1 = three_halves * j2 * pinvsq * n_recovered;
```

### Rule 2: Rational Numbers Are Exact

```cpp
// WRONG:
double coeff = 0.0625;  // what is this?

// RIGHT:
const auto one_sixteenth = T(1) / T(16);  // = 1/16, from Brouwer (1959) Eq. 38
```

### Rule 3: Polynomial Coefficients Trace to Theory

```cpp
// WRONG (current SGP4):
xmdot = n + 0.5*temp1*betao*x3thm1
      + 0.0625*temp2*betao*(13.0 - 78.0*theta2 + 137.0*theta4);

// RIGHT (our implementation):
// Brouwer (1959) secular rate of mean anomaly, J2^2 term
// P2(cos i) = (3cos²i - 1)/2 → x3thm1 = 3θ² - 1
// J2² secular polynomial in cos²i: coefficients from Brouwer Eq. 38
// 13 - 78cos²i + 137cos⁴i evaluated via Horner
const auto brouwer_j2sq_M = T(13) + theta2 * (T(-78) + T(137) * theta2);
xmdot = n_recovered
    + (T(3)/T(2)) * temp1 * betao * x3thm1
    + (T(1)/T(16)) * temp2 * betao * brouwer_j2sq_M;
// Source: Brouwer (1959) Eq. 38, via Lane & Hoots (1979) Sec. 3
```

### Rule 4: Deep Space Constants Trace to Astronomy

```cpp
// WRONG (current SGP4):
double ZNS = 1.19459E-5;

// RIGHT:
// Solar mean motion = 2π / (365.2421897 days × 1440 min/day)
// = 2π / 525949 min = 1.19459...e-5 rad/min
const auto solar_mean_motion = math::two_pi<T>()
    / (T("365.2421897") * T(1440));  // tropical year in minutes
// Source: Astronomical Almanac, tropical year length
```

### Rule 5: Curve-Fit Coefficients Are Documented Until Re-Derived

```cpp
// For the Hansen G-function polynomial fits (Phase 3 re-derivation):
// These are cubic polynomial approximations to the exact Hansen coefficient G_{lpq}(e)
// Source: Hough (1981) curve fits to Bessel-function integrals
// Distributed via Spacetrack Report #3 FORTRAN source code
// TODO(Phase 3): Replace with exact Hansen coefficient evaluation

const auto g211_coeffs = ModelValue<T>::CurveFit(
    {T("3.616"), T("-13.247"), T("16.290"), T("0.0")},
    "G_{211}(e) cubic fit, e <= 0.65",
    "Hough (1981) via Spacetrack Report #3",
    T("1e-4")  // estimated model accuracy of cubic fit
);
```

### Rule 6: Thresholds Are Precision-Dependent

```cpp
// WRONG:
if (em < 1.0e-6) em = 1.0e-6;  // why 1e-6?

// RIGHT:
// Minimum eccentricity: below this, e-dependent formulas lose significance
// For double (~16 digits): 1e-6 leaves 10 digits for e-dependent terms
// For 50-digit: could go to 1e-40, but the SGP4 model itself is meaningless there
// Use: precision_of_T / (amplification_factor_in_formula)
const auto min_eccentricity = max(
    epsilon_of<T>() * T(1e6),  // guard factor for worst-case amplification
    T("1e-6")                   // SGP4 compatibility floor
);
```

## Constant Categories and Their Treatment

| Category | Treatment | Example |
|---|---|---|
| Integer | Exact `T(n)` | `T(3)`, `T(134)`, `T(81)` |
| Simple rational | Exact `T(p)/T(q)` | `T(1)/T(16)` = 0.0625 |
| Defining parameter | `DefinedValue` from string | `T("6378137.0")` |
| Measured parameter | `MeasuredValue` with uncertainty | `GM ± σ` |
| Mathematical constant | `boost::math::constants` | `pi<T>()` |
| Derived from theory | Computed by function | `brouwer_secular_rate(J2, e2, cosi)` |
| Astronomical value | Computed from orbital elements | `solar_mean_motion = 2π/year` |
| Curve-fit coefficient | `ModelValue` with provenance | `g211_poly[0] = T("3.616")` |
| Legacy model constant | `LegacyModelValue` with provenance | `GM_sgp4_wgs84 = T("398600.5")` |
| Threshold | Precision-dependent expression | `epsilon_of<T>() * guard` |

## What This Means for the ~120 Opaque Constants

Every constant identified in `documentation_gaps.md` gets one of these treatments:

1. **Category A (Brouwer polynomials, Kaula F-functions):** Replaced by computation functions that derive the coefficients from theory. The function IS the documentation.

2. **Category B (Solar/lunar constants):** Replaced by expressions in fundamental astronomical quantities (orbital periods, obliquity, eccentricities) that are themselves `DefinedValue` or `MeasuredValue`.

3. **Category C (Hansen G-functions):** Initially `ModelValue` with full provenance string. Later replaced by exact Bessel-series evaluation.

4. **Category C (tesseral resonance):** Recomputed from EGM2008 $\bar{C}_{nm}$, $\bar{S}_{nm}$ coefficients.

**No constant remains unexplained. The code is its own documentation.**
