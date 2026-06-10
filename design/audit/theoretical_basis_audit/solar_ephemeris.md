# Theoretical Basis Audit — `src/ephemeris/solar_ephemeris.h`

## File Summary

**File**: `src/ephemeris/solar_ephemeris.h`  
**Expected function count**: 1  
**Functions**: `compute_solar_position()`  
**Lines of code**: 73

This file computes the Sun's geocentric ecliptic longitude and geocentric distance for use in SGP4 deep-space third-body perturbations. The implementation is a low-precision analytical ephemeris (accuracy ~0.01°) based on Kepler orbital mechanics with a truncated series for the equation of center.

---

## Audit Card: `compute_solar_position`

```
=== FORMULA AUDIT CARD ===
ID:                     solar_ephemeris::compute_solar_position
Location:               src/ephemeris/solar_ephemeris.h:40-72
Mathematical statement: Compute Sun's ecliptic longitude λ☉ and distance r☉ via:
                        (1) Mean anomaly: M☉(t) = M₀ + n☉ · Δt
                        (2) Equation of center: C = 2e sin(M) + (5/4)e² sin(2M)
                        (3) True anomaly: ν = M☉ + C
                        (4) Ecliptic longitude: λ☉ = ν + ω̃
                        (5) Distance [AU]: r☉ ≈ 1 − e cos(M)

THEORY
  Underlying theorem:   Kepler's equation inversion via Lagrange series.
                        The equation of center is the first two nonzero terms
                        of the Lagrange series:
                          ν − M = 2e sin(M) + (5/4)e² sin(2M) + O(e³)
                        Convergence holds for e < e_Laplace ≈ 0.6627.
                        For the Sun, e ≈ 0.01671, so convergence is rapid.
                        Distance approximation r ≈ a(1 − e cos M) is the first-
                        order Taylor in e of the exact orbit equation
                        r = a(1−e²)/(1+e cos ν).
  Primary reference:    - Meeus, J. (1998) "Astronomical Algorithms" §25.3
                        - Vallado, Crawford & Hujsa (2006) "Revisiting
                          Spacetrack Report No. 3" (SR3) §B5, page 79-80
                        - Brouwer, D. & Clemence (1961) "Methods of
                          Celestial Mechanics" Chapter 10 (Lagrange inversion)
  Domain of validity:   Mean anomaly M☉ ∈ [0, 2π]; Sun's eccentricity
                        e ∈ (0, e_Laplace); solar orbital elements valid
                        over bounded time intervals (mission duration).

METHOD
  Method declared:      Power series truncation: equation of center kept to
                        2 terms in e·sin(·). Third term (13/12)e³ sin(3M)
                        and higher omitted.
  Method implemented:   Lines 50-72:
                        - M☉ from linear propagation (line 51)
                        - sin(M), sin(2M) computed (lines 55-56)
                        - C = 2e·sin(M) + (5/4)e²·sin(2M) (line 58)
                        - ν = M☉ + C (line 61)
                        - λ☉ = ν + ω̃ (line 64)
                        - r☉ = 1 − e·cos(M) (line 69)
  Match verdict:        ✓ matched — implementation is 2-term Lagrange
                        series truncation, as cited in Meeus §25.3
                        and SR3 §B5.

ERROR BOUND
  Bound category:       precision (truncation of convergent power series)
  Bound formula:        The omitted third term is:
                          |term_3| = |(13/12) e³ sin(3M)| ≤ (13/12) e³
                        For e = 0.01671 (Sun):
                          (13/12) × (0.01671)³ ≈ 4.9 × 10⁻⁸ rad
                                              ≈ 0.01 arcsec
                        Well below stated accuracy target ~0.01° (36 arcsec).
                        Conservative bound for C: |ΔC| ≤ 5 × 10⁻⁸ rad.
  Bound implemented:    The code does not explicitly compute or add a
                        truncation bound to the result. The bound is
                        implicitly absorbed via TrackedValue error
                        propagation (sin/cos operations track precision via
                        REQ-EF-3 / AUD-EF). No explicit `result.errors.precision +=`
                        statement for the O(e³) truncation.
  Bound verdict:        ⚠ implicit — truncation bound not explicitly computed,
                        but implicitly carried by TrackedValue propagation.
                        Acceptable ONLY IF sin/cos operations are audited
                        to handle the domain |M| ∈ [0, 2π] correctly.
                        Recommendation: add a comment documenting that
                        O(e³) and higher terms are dropped, and that their
                        magnitude is below accuracy targets.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form error propagation via
                        TrackedValue arithmetic); REQ-EF-1 (all public ops
                        return TrackedValue).
  AUD-EF applies:       AUD-EF-1 (TrackedValue propagation through
                        arithmetic); AUD-EF-2 (sin/cos operations track
                        precision).
  AUD-MC applies:       n/a (not an algebra operation; used by downstream
                        third-body perturbation).
  Verification test:    - Verify against high-precision ephemeris
                        (e.g. JPL Horizons) for a range of dates.
                        - Check that |λ_code − λ_ref| < 0.01° consistently.
                        - Verify reported errors dominate actual discrepancy.

NOTES
  - The distance formula r ≈ a(1 − e cos M) is first-order in e. Exact
    formula is r = a(1−e²)/(1 + e cos ν). Error: ~0.03% at this e, acceptable
    for deep-space perturbations.
  - Code pre-computes sin(λ☉) and cos(λ☉) for efficiency in downstream
    force calculations. These inherit precision from λ☉.
  - No 2π wraparound handling; for long-duration missions, angle
    normalization may be needed (not an issue for SGP4 propagation over
    weeks to months).
  - Solar orbital elements from CelestialBody::make_sun assumed
    precomputed with full precision.
```

---

## File-Level Verdict

- **A. Error wiring**: ✓ All intermediate TrackedValue operations (sin, cos, multiply, add) propagate errors per REQ-EF-3.
- **B. Algebra axioms**: n/a (not an algebra operation; composition of scalar functions).
- **C. Theoretical basis**:
  - Method: ✓ 2-term Lagrange series truncation matches Meeus/SR3 exactly.
  - Bound: ✓ Implicit via TrackedValue plus explicit annotation block in solar_ephemeris.h:53-82 (R14, 2026-05-13). The dropped O(e³) term magnitude is documented inline (~0.01 arcsec for the Sun vs. 0.01° = 36 arcsec accuracy target — ~36000× headroom).

**File verdict: PASS** — 2-term equation of center matches cited theory exactly; truncation error documented inline and well below stated accuracy goals. R14 (2026-05-13) closed the "documentation note recommended" follow-up.

---

## References

1. **Meeus, J.** (1998). *Astronomical Algorithms*, 2nd ed. Willmann-Bell. Chapter 25: "Solar Coordinates".
2. **Vallado, D.A., Crawford, P., & Hujsa, R.** (2006). "Revisiting Spacetrack Report No. 3: Rev 2". CSSI.
3. **Brouwer, D. & Clemence, G.M.** (1961). *Methods of Celestial Mechanics*. Academic Press.
