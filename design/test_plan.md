# SGP4 Library Comprehensive Test Plan (see plan file for current roadmap)

## Overview

This document defines the test strategy for the SGP4 arbitrary-precision propagator library. Tests are organized by module layer, from foundation (math) through assembly (propagator). Each test specifies its purpose, inputs, expected results, and the physical or mathematical principle it validates.

All tests use TrackedValue<T> to verify that error propagation works correctly alongside the computation itself.

---

## Layer 0: math/tracked_value.h

### T-001: Named Constructors
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-001a | exact_integer | `exact_integer(42)` | value=42, all errors=0 | Integers are representable exactly |
| T-001b | defined | `defined("0.1")` | value≈0.1, precision=|0.1−fl(0.1)|, measurement=accuracy=0 | Representation error from decimal→binary |
| T-001c | measured | `measured("9.81", "0.01")` | value≈9.81, measurement=0.01, precision≈repr_err | Physical measurement uncertainty |
| T-001d | defined zero precision | `defined("0.5")` | precision=0 | 0.5 is exact in binary |

### T-002: Arithmetic Error Propagation
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-002a | Addition | `measured("1.0","0.01") + measured("2.0","0.02")` | measurement = √(0.01²+0.02²) | RSS error combination |
| T-002b | Subtraction cancellation | `measured("1.000001","0.01") - measured("1.0","0.01")` | reliable_digits drops sharply | Catastrophic cancellation detection |
| T-002c | Multiplication | `measured("3.0","0.1") * measured("4.0","0.2")` | measurement ≈ 3·0.2 + 4·0.1 = 1.0 | Product rule for errors |
| T-002d | Division near zero | `exact(1) / measured("1e-15","1e-16")` | measurement ≈ 1e14 | Error amplification in division |
| T-002e | Exact arithmetic | `exact(3) * exact(7)` | value=21, all errors=0 | Exact integers produce exact results |

### T-003: Transcendental Functions
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-003a | sin near zero | `sin(measured("1e-10","1e-20"))` | accuracy near zero (sin≈x for small x) | Small-angle linear regime |
| T-003b | cos near zero | `cos(measured("1e-10","1e-20"))` | accuracy from quadratic bound (x²/2) | Quadratic regime: cos(x)≈1−x²/2 |
| T-003c | cos near pi/2 | `cos(pi()/exact(2))` | value≈0, precision reflects cancellation | cos near zero of function |
| T-003d | sqrt | `sqrt(measured("4.0","0.1"))` | value=2.0, measurement=0.025 | δ(√x) = δx/(2√x) |
| T-003e | atan2 quadrants | All four quadrants | Correct angle in [-π,π] | Quadrant-correct arctangent |

### T-004: Angle Utilities
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-004a | wrap_two_pi positive | 7.0 | 7.0 − 2π | Range [0, 2π) |
| T-004b | wrap_two_pi negative | −1.0 | −1.0 + 2π | Range [0, 2π) |
| T-004c | wrap_two_pi large | 100π + 0.5 | 0.5 | Many wraps |
| T-004d | degrees_to_radians | 180° | π rad | Conversion factor |
| T-004e | radians_to_degrees | π rad | 180° | Inverse conversion |

---

## Layer 0.5: math/series.h, factorial.h, wallis.h

### T-010: Series Evaluation
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-010a | Alternating series (ln2) | Σ(-1)^{n+1}/n, tol=1e-10 | value≈0.6931471806, precision < tol | Leibniz criterion bounds |
| T-010b | Alternating series error bound | Same as T-010a | error.precision ≤ |last term| | Leibniz truncation bound |
| T-010c | Geometric series | Σ(1/2)^n, tol=1e-12 | value≈2.0 | Ratio test convergence |
| T-010d | Horner polynomial | p(x)=3x²+2x+1 at x=2 | value=17 | Horner evaluation |

### T-011: Factorials and Binomials
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-011a | factorial(0) | 0 | 1 | Definition: 0! = 1 |
| T-011b | factorial(10) | 10 | 3628800, exact | Integer factorial |
| T-011c | double_factorial(7) | 7 | 105 = 7·5·3·1 | Odd double factorial |
| T-011d | generalized_binomial | C(−1/2, 3) | −5/16 | Non-integer upper parameter |
| T-011e | generalized_binomial | C(−5/2, 3) | −105/16 | Verified in derivation |
| T-011f | falling_factorial | (−1/2)(−3/2)(−5/2) for k=3 | −15/8 | Product formula |

### T-012: Wallis Integrals
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-012a | wallis_odd(0) | k=0 | 1 | W₁ = 1 |
| T-012b | wallis_odd(1) | k=1 | 2/3 | W₃ = 2!!/3!! |
| T-012c | wallis_odd(2) | k=2 | 8/15 | W₅ = 4!!/5!! |
| T-012d | wallis_even(1) | k=1 | π/4 | W₂ = (1!!/2!!)·(π/2) |
| T-012e | sin_power_cos_integral | k=3 | 1/7 | ∫₀^{π/2} sin⁶(φ)cos(φ)dφ |

---

## Layer 1: math/kepler.h, math/vector3.h

### T-020: Standard Kepler Solver
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-020a | Circular orbit | M=1.5, e=0.0 | E=1.5 | E=M when e=0 |
| T-020b | Low eccentricity | M=1.5, e=0.1 | E≈1.6054 (verify) | Converges in ~3 iterations |
| T-020c | Moderate eccentricity | M=π/2, e=0.5 | Verify f(E)=0 to tolerance | Newton/Halley convergence |
| T-020d | High eccentricity | M=0.5, e=0.9 | Verify f(E)=0, check iteration count | Challenging convergence |
| T-020e | Near-parabolic | M=0.1, e=0.99 | Verify f(E)=0, may need many iterations | Boundary of elliptic regime |
| T-020f | Convergence precision | Any | error.precision ≤ tolerance | Residual added to precision error |

### T-021: SGP4 Modified Kepler Solver
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-021a | Circular (axN=0, ayn=0) | U=1.5 | x=1.5 | Trivial case: x=U |
| T-021b | Standard form equivalence | axN=e, ayn=0, U=M | Same E+ω as standard solver + ω | Algebraic equivalence |
| T-021c | Halley vs Newton | Same inputs | Same result, fewer iterations (Halley) | Cubic vs quadratic convergence |

### T-022: Vector3 Operations
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-022a | Dot product orthogonal | (1,0,0)·(0,1,0) | 0 | Perpendicular = zero dot |
| T-022b | Dot product parallel | (3,0,0)·(4,0,0) | 12 | Parallel = product of magnitudes |
| T-022c | Cross product | (1,0,0)×(0,1,0) | (0,0,1) | Right-hand rule |
| T-022d | Cross product anti-commutative | a×b vs b×a | Negation | a×b = −b×a |
| T-022e | Magnitude | (3,4,0) | 5 | Pythagorean triple |
| T-022f | Error propagation through dot | Measured vectors | Errors combine per product rule | TrackedValue arithmetic chain |

---

## Layer 2: Geodesy

### T-030: WGS84 Ellipsoid (from 1/f)
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-030a | Semi-minor axis b | a=6378137, 1/f=298.257223563 | b = 6356752.3142... | b = a(1−f) |
| T-030b | First eccentricity | Same | e = 0.0818191908... | e² = 2f−f² |
| T-030c | J₂ from ellipsoid | Same | J₂ ≈ 1.08263e-3 | Level surface formula |
| T-030d | Normal gravity equator | φ=0° | γₑ = 9.7803253... m/s² | Somigliana at equator |
| T-030e | Normal gravity pole | φ=90° | γₚ = 9.8321849... m/s² | Somigliana at pole |
| T-030f | 50-digit precision | cpp_bin_float_50 | Match NGA Appendix B to 50 digits | Guard digit verification |

### T-031: WGS72 Ellipsoid (from J₂)
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-031a | J₂ roundtrip | J₂=0.001082616 | J2n(1) = 0.001082616 exactly | Defining parameter preserved |
| T-031b | Reconstructed 1/f | Same | 1/f ≈ 298.26 | Iterative solution from J₂ |
| T-031c | xke constant | Same | 0.0743669161... | √(GM/aₑ³)·60 |

---

## Layer 2: Astronomy

### T-040: Fundamental Constants
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-040a | Solar mean motion | 365.257 day period | ZNS = 1.19459e-5 rad/min | 2π/(T·1440) |
| T-040b | Lunar anomalistic rate | 27.5546 day period | ZNL = 1.58352e-4 rad/min | Anomalistic month |
| T-040c | Lunar sidereal rate | 27.3216 day period | 0.22997 rad/day | Sidereal month (distinct from anomalistic) |
| T-040d | Obliquity | cos(ε)=0.91744867 | ε ≈ 23.4441° | Axial tilt |

---

## Layer 2: Perturbation Theory

### T-050: Brouwer Secular Rates
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-050a | Equatorial orbit | i=0° | Ω̇ = −TEMP1 (max regression) | cos(0)=1 |
| T-050b | Polar orbit | i=90° | Ω̇ = 0 | cos(90°)=0 → no regression |
| T-050c | Critical inclination | i=63.43° | ω̇ = 0 | 1−5cos²i = 0 |
| T-050d | Sun-sync orbit | i=97.4°, a=7078km | Ω̇ ≈ 0.9856°/day | Required for sun-sync |
| T-050e | J₂² polynomial | Any | 13−78θ²+137θ⁴ evaluates correctly | Horner form matches expanded |
| T-050f | Accuracy error bound | Any | accuracy ≥ O(J₂³n/p⁶) | Truncation error estimate |
| T-050g | Sat 00005 rates | From TLE | Match diagnostic output to 1e-14 relative | Integration test |

### T-051: Kaula Inclination Functions
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-051a | F₂₂₀(51.6°) | ISS inclination | 1.96962... | Exact rational coefficients |
| T-051b | F₂₂₁(51.6°) | ISS inclination | 0.92239... | Exact rational coefficients |
| T-051c | F₄₄₂(51.6°) | ISS inclination | 14.8891... (uses 315/8) | No FORTRAN decimal approximation |
| T-051d | Unimplemented (l,m,p) | (10,5,3) | Returns 0 with accuracy error flag | Graceful unknown handling |

---

## Layer 2: TLE Parsing

### T-060: TLE Format
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-060a | Standard TLE | Sat 00005 from SGP4-VER.TLE | All fields parsed correctly | Column-fixed format |
| T-060b | Eccentricity leading decimal | "1859667" → 0.1859667 | Correct interpretation | Implied "0." prefix |
| T-060c | Negative BSTAR | "−11606−4" → −1.1606e-4 | Exponent notation | Implied decimal |
| T-060d | Measurement errors | Any TLE | inclination error ≈ 8.7e-7 rad | Format-derived uncertainty |
| T-060e | Alpha-5 satellite number | "A0005" | Parsed correctly | Extended numbering |

---

## Layer 3: SGP4 Propagator (Integration)

### T-070: Element Recovery
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-070a | Sat 00005 a₀″ | From TLE | 1.353899820603 ER | δ₁, δ₀ iteration |
| T-070b | Sat 00005 n₀″ | From TLE | 0.047206301559 rad/min | Recovered mean motion |
| T-070c | Perigee altitude | From TLE | 651.33 km | (a₀″(1−e₀)−1)·XKMPER |
| T-070d | Deep-space classification | Period ≥ 225 min | Throws runtime_error | Deep-space not yet implemented |

### T-071: Drag Coefficients
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-071a | C₁ through C₅ | Sat 00005 | Match diagnostic values to 1e-14 relative | Lane-Hoots drag model |
| T-071b | D₂, D₃, D₄ | Sat 00005 | Match diagnostic values | Higher-order drag polynomial |
| T-071c | Time coefficients | Sat 00005 | t2cof through t5cof match | Taylor expansion coefficients |
| T-071d | Low-perigee adjustment | Perigee < 156 km | s₄ and qoms24 adjusted | Atmospheric density switch |

### T-072: Near-Space Propagation (Reference: tcppver.out)
| ID | Test | Input | Expected | Tolerance | Principle |
|----|------|-------|----------|-----------|-----------|
| T-072a | Sat 00005, t=0 | t=0 min | x=7022.465, y=−1400.083, z=0.040 | <1e-8 km | Initial state |
| T-072b | Sat 00005, t=360 | t=360 min | x=−7154.031, y=−3783.177, z=−3536.194 | <1e-4 km | 6-hour propagation |
| T-072c | Sat 00005, t=720 | t=720 min | Reference from tcppver.out | <1e-4 km | 12-hour propagation |
| T-072d | Sat 00005, t=1080 | t=1080 min | Reference from tcppver.out | <2e-4 km | 18-hour propagation |
| T-072e | Sat 00005, t=1440 | t=1440 min | Reference from tcppver.out | <3e-4 km | 24-hour propagation |
| T-072f | All near-earth SGP4-VER.TLE | 33 cases, filter period<225 | Match tcppver.out | <1e-3 km | Full validation suite |

### T-073: Model Configuration
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-073a | sgp4_standard preset | select("sgp4_standard") | WGS72 J₂=0.001082616, J₄=−1.65597e-6 | SR3 constants |
| T-073b | sgp4_wgs84 preset | select("sgp4_wgs84") | WGS84 J₂=0.00108263, J₄=−1.61099e-6 | Vallado getgravconst |
| T-073c | config.Jn(2) | Any preset | Matches ellipsoid.J2n(1) | J₂ from ellipsoid |
| T-073d | config.Jn(3) | WGS72 | −2.53881e-6 | Odd zonal from empirical table |
| T-073e | config.Jn(4) | WGS72 | −1.65597e-6 (empirical, NOT ellipsoid) | Matched-pair J₄ |
| T-073f | config.Jn(99) odd | Any | 0 (not in table) | Missing odd returns zero |
| T-073g | config.Jn(6) even | WGS72 | Ellipsoid fallback | Even not in table → ellipsoid |
| T-073h | config.A(3,0) | WGS72 | −(−2.53881e-6) = +2.53881e-6 | A_{n,0} = −Jₙ |
| T-073i | config.A(3,1) | Any | 0 (tesseral not supported) | Graceful unimplemented |

### T-074: Kepler Solver Comparison
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-074a | Newton vs Halley | Same TLE, same tolerance | Position difference < 1e-12 km | Same math, different convergence rate |
| T-074b | Halley iteration count | High-eccentricity TLE | Fewer iterations than Newton | Cubic vs quadratic convergence |

### T-075: Error Budget Decomposition
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-075a | Position measurement error | Any propagation | Non-zero, from TLE uncertainties | Error flows from input |
| T-075b | Position precision error | Any propagation | Non-zero, grows with time | Arithmetic accumulation |
| T-075c | Position accuracy error | Any propagation | Non-zero, from Brouwer truncation | Model truncation |
| T-075d | Total ≥ max(components) | Any propagation | total_error ≥ max(meas, prec, acc) | Error combination |

### T-076: Failure Modes
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-076a | Deep-space TLE | Sat 04632 (period=1198 min) | Throws runtime_error | Explicit failure, not silent zero |
| T-076b | Sidereal time call | Call sidereal_time lambda | Throws runtime_error | Not yet implemented |
| T-076c | Drag lambda call | Call drag_coefficients lambda | Throws runtime_error | Not yet extracted |

---

## Layer 4: Multiprecision Verification

### T-080: Precision Scaling
| ID | Test | Input | Expected | Principle |
|----|------|-------|----------|-----------|
| T-080a | double vs cpp_bin_float_50 | Same TLE, same propagation | Position agrees to ~15 digits | Guard digits don't change result |
| T-080b | 50-digit WGS84 constants | cpp_bin_float_50 | Match NGA Appendix B all 50 digits | No precision loss in derivation chain |
| T-080c | Reliable digits at t=0 | cpp_bin_float_50 | ~8 reliable (limited by GM uncertainty) | Measurement error dominates |

---

## Execution Strategy

### Existing Test Executables
| Test | Executable | Build | Status |
|------|-----------|-------|--------|
| test_math | tests/test_math/main.cpp | MSBuild or manual | Covers T-001 through T-022 partially |
| test_geodesy | tests/test_geodesy/main.cpp | MSBuild | Covers T-030, T-031 |
| test_astronomy | tests/test_astronomy/main.cpp | MSBuild | Covers T-040 |
| test_perturbation | tests/test_perturbation/main.cpp | MSBuild | Covers T-050, T-051 |
| test_tle | tests/test_tle/main.cpp | MSBuild | Covers T-060 |
| test_wgs84 | tests/test_wgs84/main.cpp | MSBuild | Covers T-030f (50-digit) |
| test_sgp4 | build_test.bat | Direct cl.exe | Covers T-070 through T-076 |

### Tests Requiring New Code
| Test IDs | Gap | Priority |
|----------|-----|----------|
| T-002b,d | Subtractive cancellation, division near zero | HIGH |
| T-003a-e | Transcendental edge cases | HIGH |
| T-011a-f | Factorial/binomial unit tests | MEDIUM |
| T-020d-e | High-eccentricity Kepler | HIGH |
| T-022a-f | Vector3 unit tests | MEDIUM |
| T-050a-d | Brouwer limiting cases (equatorial, polar, critical, sun-sync) | HIGH |
| T-072f | All 33 SGP4-VER.TLE near-earth cases | HIGH |
| T-074a-b | Kepler solver comparison | MEDIUM |
| T-075a-d | Error budget output | MEDIUM |
| T-076a-c | Exception tests | HIGH |
| T-080a-c | Multiprecision comparison | LOW |

### Pass Criteria
- All `compare()` calls: relative error < 1e-6 (or absolute < 1e-15 for near-zero values)
- All propagation: position error < tolerance specified per test case
- All exceptions: must throw `std::runtime_error` with descriptive message
- Zero compiler warnings with `/W4`
- Zero BAD flags in diagnostic output
