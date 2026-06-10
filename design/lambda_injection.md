# Lambda Injection for Model Computations

## Principle

The model initialization function accepts all computational formulas as lambda functions (or `std::function`). No formula is hardcoded in the propagator. Every computation that might change between model versions — sidereal time, ephemeris rates, atmospheric density, precession-nutation, perturbation theory — is injected at construction time.

This means:
- Switching from IAU 1982 GMST to IAU 2006 ERA: inject a different sidereal time lambda
- Switching from simplified lunar ephemeris to JPL DE: inject a different lunar position lambda
- Switching from Brouwer J₂² truncation to higher-order theory: inject a different secular rate lambda
- Using the Astronomical Almanac 2025 instead of 1980 constants: inject different ephemeris lambdas

The propagator is a GENERIC computation engine. It has NO knowledge of which specific formulas it uses. It just calls the lambdas.

## Example: Sidereal Time

```cpp
// The sidereal time function type:
// Takes epoch (as TrackedValue) and returns Greenwich sidereal time (as TrackedValue)
template<typename T>
using SiderealTimeFn = std::function<TrackedValue<T>(const TrackedValue<T>& epoch_jd)>;

// IAU 1982 GMST (standard SGP4)
template<typename T>
TrackedValue<T> gmst_aoki_1982(const TrackedValue<T>& jd) {
    // Aoki et al. (1982) polynomial
    // ...
}

// IAU 2006 ERA (modern)
template<typename T>
TrackedValue<T> era_iau_2006(const TrackedValue<T>& jd) {
    // ERA = 2*pi*(0.7790572732640 + 1.00273781191135448 * Du)
    // ...
}

// Propagator accepts whichever one the caller provides
template<typename T>
class Propagator {
    SiderealTimeFn<T> sidereal_time_fn_;
    // ...
};
```

## Example: Lunar Position

```cpp
// The lunar ephemeris function type:
// Takes time from reference epoch, returns lunar orbital elements
template<typename T>
using LunarEphemerisFn = std::function<LunarElements<T>(const TrackedValue<T>& day)>;

// Simplified Hujsak (1979) ephemeris (standard SGP4)
template<typename T>
LunarElements<T> lunar_hujsak_1979(const TrackedValue<T>& day) {
    // Linear rates from SR3 page 59
    // XNODCE = 4.5236020 - 9.2422029E-4 * day
    // ...
}

// Could be replaced with JPL DE ephemeris for enhanced mode
```

## Full Injection Interface

The propagator model initialization accepts a struct of lambdas:

```cpp
template<typename T>
struct ModelFunctions {
    // Sidereal time
    SiderealTimeFn<T> sidereal_time;

    // Lunar/solar ephemeris
    LunarEphemerisFn<T> lunar_ephemeris;
    SolarEphemerisFn<T> solar_ephemeris;

    // Brouwer secular rates
    SecularRatesFn<T> secular_rates;

    // Kaula inclination functions
    InclinationFn<T> inclination_function;

    // Hansen eccentricity functions
    EccentricityFn<T> eccentricity_function;

    // Atmospheric drag model
    DragFn<T> drag_model;

    // Kepler equation solver
    KeplerFn<T> kepler_solver;

    // Factory for standard SGP4 model
    static ModelFunctions standard_sgp4();

    // Factory for enhanced model (modern constants)
    static ModelFunctions enhanced();
};
```

## What This Enables

1. **Standard SGP4 compatibility:** `ModelFunctions::standard_sgp4()` provides all the original formulas with 1970s constants.

2. **Astronomical Almanac updates:** Create a new factory `ModelFunctions::almanac_2025()` that injects updated ephemeris rates, sidereal time formula, and precession-nutation model.

3. **Research mode:** Replace individual lambdas to test the effect of each formula change in isolation (e.g., "what happens if I use ERA instead of GMST but keep everything else the same?").

4. **Testing:** Inject mock lambdas that return known values for unit testing specific computations.

5. **No recompilation:** Different model configurations are selected at runtime by choosing which factory method to call, not by recompiling with different #defines.

## Discovery Documentation

Every lambda implementation includes a comment documenting:
- What it computes
- Which equation or derivation it implements
- What the specific numeric constants mean
- Any findings from the magic number test suite
