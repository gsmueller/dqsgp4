# Library Architecture

## Principle

Every numeric value is computed from inputs. There are no magic numbers. The library is organized by the mathematical domain that generates the values, not by the consumer that uses them.

## Module Dependency

```
                    ┌─────────────────────┐
                    │   Caller provides:   │
                    │  a, 1/f (or J₂),    │
                    │  GM, ω              │
                    │  + uncertainties     │
                    │  + TLE elements      │
                    │  + precision target  │
                    └────────┬────────────┘
                             │
              ┌──────────────┼──────────────┐
              ▼              ▼              ▼
    ┌──────────────┐ ┌─────────────┐ ┌───────────────┐
    │ math         │ │ astronomy   │ │ tle_parser    │
    │              │ │             │ │               │
    │ TrackedValue │ │ Solar/lunar │ │ Parse TLE     │
    │ Series eval  │ │ orbital     │ │ text into     │
    │ Factorials   │ │ elements    │ │ numeric       │
    │ Binomial     │ │ from fund.  │ │ elements      │
    │ Wallis       │ │ periods     │ │               │
    │ Kepler solv  │ │             │ │               │
    └──────┬───────┘ └──────┬──────┘ └───────┬───────┘
           │                │                │
           ▼                │                │
    ┌──────────────┐        │                │
    │ geodesy      │        │                │
    │              │        │                │
    │ Equipotent.  │        │                │
    │ Ellipsoid    │        │                │
    │ (a,1/f,GM,ω) │        │                │
    │ → e², b, q₀  │        │                │
    │ → γₑ, γₚ, k  │        │                │
    │ → J₂ₙ, U₀   │        │                │
    │ → Somigliana │        │                │
    └──────┬───────┘        │                │
           │                │                │
           ▼                ▼                ▼
    ┌──────────────────────────────────────────────┐
    │ perturbation_theory                          │
    │                                              │
    │ Brouwer secular rates (from J₂, J₄, e², i)  │
    │ Kaula inclination functions F_lmp(i)         │
    │ Hansen eccentricity functions G_lpq(e)       │
    │ Resonance coefficients (from J_nm via EGM)   │
    │ Third-body perturbations (from astronomy)    │
    └──────────────────┬───────────────────────────┘
                       │
                       ▼
    ┌──────────────────────────────────────────────┐
    │ sgp4                                         │
    │                                              │
    │ SGP4Propagator(ellipsoid, perturbations,     │
    │                tle_elements, astronomy)       │
    │                                              │
    │ Initialize() → near-space or deep-space      │
    │ Propagate(t) → position, velocity            │
    │               with three errors              │
    └──────────────────────────────────────────────┘
```

## Module Details

### `math/`

Pure mathematics. No domain knowledge. No physical constants.

| File | Contains | Depends On |
|---|---|---|
| `tracked_value.h` | `TrackedValue<T>` with three-error propagation through every operation | Nothing |
| `series.h` | Alternating series evaluator with error bounds, Horner evaluation, series quotient | `tracked_value.h` |
| `factorial.h` | Falling factorial, rising factorial, double factorial — exact integer/rational | `tracked_value.h` |
| `binomial.h` | Generalized binomial coefficient $\binom{\alpha}{k}$ | `factorial.h` |
| `wallis.h` | Wallis integrals $W_n$ — exact rational values | `factorial.h` |
| `continued_fraction.h` | Lentz algorithm for continued fraction evaluation | `tracked_value.h` |
| `kepler.h` | Kepler equation solver (Halley/Householder) | `tracked_value.h` |
| `vector3.h` | `Vector3<TrackedValue<T>>` with arithmetic | `tracked_value.h` |

### `geodesy/`

The equipotential ellipsoid and everything derived from it. Accepts four defining parameters as input. Computes all derived constants using `math/` functions.

| File | Contains | Depends On |
|---|---|---|
| `equipotential_ellipsoid.h` | Central class. Constructed from ($a$, $1/f$ or $J_2$, $GM$, $\omega$). Computes and stores all derived constants. | `math/*` |

The `EquipotentialEllipsoid<T>` constructor does:

```
Input: a, inv_f (or J2), GM, omega — each as TrackedValue<T>

Geometric chain (Derivation 003):
  f = 1 / inv_f
  e2 = 2*f - f*f
  e = sqrt(e2)
  e_prime = e / sqrt(1 - e2)
  b = a * (1 - f)
  E_lin = a * e
  c = a*a / b

Series evaluations (Derivations 001-002):
  q0 = series in e_prime (math/series.h)
  q0_prime = series in e_prime (math/series.h)

Physical chain (Derivation 004):
  m = omega*omega * a*a * b / GM
  gamma_e = (GM / (a*b)) * (1 - m - m*e_prime*q0_prime / (6*q0))
  gamma_p = (GM / (a*a)) * (1 + m*e_prime*q0_prime / (3*q0))
  k = b*gamma_p / (a*gamma_e) - 1
  U0 = series form (Derivation 004 Step 5)

Zonal harmonics (Derivation 006):
  J2 = (e2/3) * (1 - 2*m*e_prime / (15*q0))
  J2n(n) = formula from e2, J2

Radii (Derivation 005):
  R1 = a * (1 - f/3)
  R2 = c * sqrt(integral_series)
  R3 = cbrt(a*a*b)

Every value is a TrackedValue<T> carrying three errors.
No constant is stored. Everything is computed.
```

Methods:

```
normal_gravity(phi) → TrackedValue<T>
  // Somigliana: gamma_e * (1 + k*sin²φ) / sqrt(1 - e²*sin²φ)

J2n(n) → TrackedValue<T>
  // Derivation 006 formula
```

### `astronomy/`

Computes solar/lunar orbital elements and obliquity from fundamental periods. These are the inputs that the deep space perturbation theory needs.

| File | Contains | Depends On |
|---|---|---|
| `solar_system.h` | Solar/lunar mean motions, eccentricities, obliquity — computed from fundamental orbital periods | `math/tracked_value.h` |

Example: instead of `ZNS = 1.19459E-5` as a magic number, this module computes:

```
solar_mean_motion = two_pi / (tropical_year_in_minutes)
```

where `tropical_year_in_minutes` is an input (with its uncertainty).

The obliquity of the ecliptic, the lunar node regression rate, the lunar orbital eccentricity — all computed from fundamental astronomical quantities that are themselves inputs.

### `perturbation_theory/`

Computes the perturbation coefficients that the SGP4 propagator needs. This is where the Brouwer theory, Kaula functions, and Hansen coefficients live.

| File | Contains | Depends On |
|---|---|---|
| `brouwer.h` | Secular rates ($\dot{M}$, $\dot{\omega}$, $\dot{\Omega}$) as functions of $J_2$, $J_4$, $e^2$, $\cos i$ | `math/*`, `geodesy/` |
| `kaula.h` | Inclination functions $F_{lmp}(i)$ — computed from recursion | `math/factorial.h` |
| `hansen.h` | Eccentricity functions $G_{lpq}(e)$ — Bessel series or polynomial fit | `math/series.h` |
| `third_body.h` | Solar/lunar perturbation coefficients — computed from astronomy inputs | `astronomy/`, `math/*` |
| `resonance.h` | Tesseral resonance coefficients — computed from $J_{nm}$ via EGM or from gravity model inputs | `math/*` |

The Brouwer secular rates emerge from the perturbation theory as functions of ($J_2$, $J_4$, $e^2$, $\cos^2 i$). The polynomial coefficients (13, 78, 137, etc.) are not stored — they are computed by the Brouwer formulas. The Brouwer module IS the derivation.

The accuracy error ($\delta_a$) enters here: the Brouwer theory is truncated at $J_2^2$ and $J_4$. The magnitude of the first omitted term ($J_2^3$) bounds the accuracy error. This is computed and carried through.

### `tle/`

Parses TLE text into numeric elements.

| File | Contains | Depends On |
|---|---|---|
| `tle_parser.h` | Parse TLE line 1 and line 2 into `TleData` struct | Nothing (string parsing) |
| `tle_parser.cpp` | Implementation — supports Alpha-5, 9-digit catalog numbers | |

The parsed elements become `TrackedValue<T>` with:
- $\sigma_m$ from the TLE format precision (8 decimal digits → known measurement error)
- $\delta_p$ = 0 (the parsing is exact)
- $\delta_a$ = 0 (no model at this stage)

### `sgp4/`

The propagator itself. Assembles the pieces from geodesy, perturbation theory, astronomy, and TLE parsing.

| File | Contains | Depends On |
|---|---|---|
| `sgp4_propagator.h` | `SGP4Propagator<T>` — Initialize + Propagate | Everything above |
| `state_vector.h` | `StateVector<T>` = `Vector3<TrackedValue<T>>` position + velocity | `math/vector3.h` |

The propagator does NOT contain any constants. It receives:
- An `EquipotentialEllipsoid<T>` (from geodesy)
- Perturbation coefficients (from perturbation_theory)
- Astronomical elements (from astronomy)
- TLE elements (from tle_parser)

And produces `StateVector<T>` with all three errors propagated to the output.

## File Layout

```
src/
  math/
    tracked_value.h
    series.h
    factorial.h
    binomial.h
    wallis.h
    continued_fraction.h
    kepler.h
    vector3.h
  geodesy/
    equipotential_ellipsoid.h
  astronomy/
    solar_system.h
  perturbation_theory/
    brouwer.h
    kaula.h
    hansen.h
    third_body.h
    resonance.h
  tle/
    tle_parser.h
    tle_parser.cpp
  sgp4/
    sgp4_propagator.h
    state_vector.h
  main.cpp
build.bat
```

## What Changed from the Original Plan

The original plan had `sgp4/gravity_model.h` storing constants per model variant (wgs72old, wgs72, wgs84). This is gone. Instead, the caller constructs an `EquipotentialEllipsoid<T>` with the desired inputs. The SGP4 propagator doesn't know or care which "model" it is — it just uses the ellipsoid it was given.

The original plan had constants as static functions returning `T`. This is gone. Constants are computed values stored in the ellipsoid instance, carrying three errors.

The original plan had separate `math/constants.h` for pi, two_pi, etc. Pi is now `boost::math::constants::pi<T>()` wrapped in a `TrackedValue` at point of use — a mathematical constant with $\sigma_m = 0$, $\delta_p$ = type representation error, $\delta_a = 0$.
