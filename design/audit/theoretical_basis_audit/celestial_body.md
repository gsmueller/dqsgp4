# Theoretical Basis Audit — `src/ephemeris/celestial_body.h`

## Overview

File: `src/ephemeris/celestial_body.h` (lines 1–112)

This file defines two static factory methods that construct `CelestialBody<T>` structs for the Sun and Moon from a `FundamentalConstants` bundle. The functions populate geocentric orbital elements (mean anomaly, longitude rate, eccentricity, inclination, node regression) from constant-valued inputs.

**Expected function count**: 2 (`make_sun`, `make_moon`)

**Audit scope**: Theoretical basis and error propagation for each formula assignment.

---

## Card 1: `make_sun(const auto& fc)`

```
=== FORMULA AUDIT CARD ===
ID:                     celestial_body::make_sun
Location:               src/ephemeris/celestial_body.h:58-81
Mathematical statement: Construct a CelestialBody with solar geocentric orbital
                        elements from FundamentalConstants: mean anomaly rate
                        (2π / solar_anomaly_period_days), longitude rate
                        (2π / solar_anomaly_period_days), eccentricity,
                        argument of perigee, and inclination ≡ 0.

THEORY
  Underlying theorem:   Orbital mechanics: the geocentric mean longitude
                        increases linearly with time at the mean solar motion
                        (Definition 25.2 of Ch 25 solar ephemeris). The Sun's
                        heliocentric orbit, when reflected to the geocentric
                        frame, has zero ecliptic inclination by geometry
                        (Ch 25, §25.6). All elements derive from Astronomical
                        Almanac (and matched-pair principle, REQ-CP-2,
                        Ch 3 for numerical values).
  Primary reference:    Hoots & Roehrich 1980, Spacetrack Report #3, Appendix B
                        (Solar and Lunar orbital elements); Astronomical
                        Almanac, annual editions (solar period, eccentricity);
                        Vallado et al. (2006) §3.4 (geocentric solar elements).
  Domain of validity:   For third-body perturbation of LEO/GEO satellites.
                        Solar period ~365.25 days; eccentricity ~0.01671.
                        Valid over SGP4 epoch range (1957–2050+).

METHOD
  Method declared:      Closed-form assignment of constants from inputs to
                        struct fields. No iterative solver, no series
                        truncation. Each field is assigned exactly one value
                        from fc by arithmetic (division, assignment).
  Method implemented:   Lines 59–80: each field assigned once. Mean anomaly
                        rate = 2π / (period × min_per_day). Longitude rate =
                        2π / period. Eccentricity, arg perigee, inclination
                        directly from fc; inclination set to exact(0).
  Match verdict:        ✓ matched — implementation is pure closed-form
                        assignment of defined and measured constants.

ERROR BOUND
  Bound category:       precision (ULP error from arithmetic) + measurement
                        (inherited from fc constants' measurement error).
  Bound formula:        Per REQ-EF-3 (closed-form propagation), each output
                        field is a TrackedValue<T> whose error is the
                        composition of its input errors via arithmetic.
                        Specifically:
                        - mean_anomaly_rate_rad_min = 2π / (period × min_per_day)
                          → errors propagate from period.errors and
                          min_per_day.errors via division rules.
                        - orbit_eccentricity = fc.solar_eccentricity
                          → errors copied directly from fc.
                        - All other fields: direct assignment with no
                          arithmetic, so errors = input.errors.
                        See REQ-EF-3 for the closed-form formula per operation
                        type (division, multiplication, assignment).
  Bound implemented:    The function returns a CelestialBody with each field
                        populated by arithmetic operations. Each operation
                        (/, *, =) propagates input.errors to output.errors
                        automatically via TrackedValue<T> operator overloads.
                        No explicit bound is written; propagation is implicit.
                        ✓ Matches the REQ-EF-3 model.
  Bound verdict:        ✓ matched — closed-form arithmetic propagates error
                        categories automatically. No manual bound-adding needed.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form error propagation on all
                        arithmetic operations).
  AUD-EF applies:       AUD-EF-1 (all public operations return TrackedValue<T>);
                        AUD-EF-2 (all return values have errors populated).
  AUD-MC applies:       n/a (constructors, not algebra operations).
  Verification test:    tests/test_constants_provider.cpp (verify that
                        make_sun returns consistent fields; cross-check mean
                        anomaly rate against published solar period).

NOTES
  - The function is a static factory, not an instance method. It constructs
    and returns a CelestialBody struct (value semantics).
  - Line 76–78: inclination = exact<T>(0) and node_epoch = exact<T>(0),
    node_rate = exact<T>(0). These encode the geometric fact that the Sun's
    ecliptic latitude is zero by construction. The exact<T>(0) call creates
    a TrackedValue with zero error, which is correct (geometric certainty).
  - Solar mean anomaly epoch (line 66) is set to fc.solar_anomaly_epoch.
    This is a constant from FundamentalConstants, carrying its measurement
    error from the constants provider.
  - The function does not validate inputs (e.g., period > 0); validation is
    delegated to ConstantsProvider construction (REQ-CP-3).
```

---

## Card 2: `make_moon(const auto& fc)`

```
=== FORMULA AUDIT CARD ===
ID:                     celestial_body::make_moon
Location:               src/ephemeris/celestial_body.h:84-108
Mathematical statement: Construct a CelestialBody with lunar geocentric orbital
                        elements from FundamentalConstants: mean anomaly rate
                        (2π / lunar_anomalistic_period_days), sidereal
                        longitude rate (2π / lunar_sidereal_period_days),
                        eccentricity, argument of perigee, inclination ~5.145°,
                        node regression rate (−2π / lunar_node_period_days).

THEORY
  Underlying theorem:   Orbital mechanics of the Moon: geocentric mean
                        longitude and anomaly increase linearly with the
                        sidereal and anomalistic periods, respectively
                        (Definition 26.2 of Ch 26 lunar ephemeris). The
                        ascending node regresses at a negative rate determined
                        by solar perturbation (Ch 26, §26.2). All rates and
                        elements are from Astronomical Almanac and SGP4
                        matched-pair constants (Hoots & Roehrich 1980, SR3).
  Primary reference:    Hoots & Roehrich 1980, Spacetrack Report #3, Appendix B
                        (Solar and Lunar orbital elements); Astronomical
                        Almanac, annual editions (lunar periods, inclination,
                        node regression rate); Vallado et al. (2006) §3.4
                        (geocentric lunar elements).
  Domain of validity:   For third-body perturbation of LEO/GEO satellites.
                        Lunar sidereal period ~27.32 days, anomalistic ~27.55
                        days, node regression period ~18.6 years. Valid over
                        SGP4 epoch range (1957–2050+).

METHOD
  Method declared:      Closed-form assignment of constants from inputs to
                        struct fields. No iterative solver, no series
                        truncation. Each field is assigned exactly one value
                        from fc by arithmetic (division, negation, assignment).
  Method implemented:   Lines 85–107: each field assigned once. Mean anomaly
                        rate = 2π / (anomalistic_period × min_per_day);
                        sidereal longitude rate = 2π / sidereal_period.
                        Node rate = −2π / node_period (negative: regression).
                        Inclination, eccentricity, arg perigee, node epoch,
                        longitude epoch directly from fc.
  Match verdict:        ✓ matched — implementation is pure closed-form
                        assignment of defined and measured constants.

ERROR BOUND
  Bound category:       precision (ULP error from arithmetic) + measurement
                        (inherited from fc constants' measurement error).
  Bound formula:        Per REQ-EF-3 (closed-form propagation), each output
                        field is a TrackedValue<T> whose error is the
                        composition of input errors via arithmetic. Specifically:
                        - mean_anomaly_rate_rad_min = 2π / (period × min_per_day)
                          → errors propagate from period.errors and
                          min_per_day.errors via division rules.
                        - longitude_rate_rad_day = 2π / sidereal_period
                          → errors propagate via division.
                        - node_rate_rad_day = −2π / node_period
                          → errors propagate via division and unary negation.
                        - All scalar fields (eccentricity, inclination, etc.):
                          direct assignment with no arithmetic, so
                          errors = input.errors.
                        See REQ-EF-3 for the closed-form formula per operation
                        type (division, multiplication, negation, assignment).
  Bound implemented:    The function returns a CelestialBody with each field
                        populated by arithmetic operations. Each operation
                        (/, −, *, =) propagates input.errors to output.errors
                        automatically via TrackedValue<T> operator overloads.
                        No explicit bound is written; propagation is implicit.
                        ✓ Matches the REQ-EF-3 model.
  Bound verdict:        ✓ matched — closed-form arithmetic propagates error
                        categories automatically. No manual bound-adding needed.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form error propagation on all
                        arithmetic operations).
  AUD-EF applies:       AUD-EF-1 (all public operations return TrackedValue<T>);
                        AUD-EF-2 (all return values have errors populated).
  AUD-MC applies:       n/a (constructors, not algebra operations).
  Verification test:    tests/test_constants_provider.cpp (verify that
                        make_moon returns consistent fields; cross-check mean
                        anomaly rate against published lunar anomalistic and
                        sidereal periods).

NOTES
  - Line 92: mean_anomaly_epoch is set to exact<T>(0), not fc's value. This
    reflects SGP4 epoch handling (SR3, §B): the lunar mean anomaly at epoch
    is defined as zero in the SGP4 reference frame. The exact<T>(0) call is
    correct (zero measurement and precision error by design).
  - Line 95: sidereal longitude rate (not anomalistic) is the "longitude"
    rate. This is the rate at which ℓ_☾ advances (Definition 26.2.1, Ch 26).
  - Line 99: orbit_arg_perigee = fc.lunar_perigee_epoch. This is the
    argument of perigee longitude at epoch, not the time-derivative. The
    comment confirms it's a longitude, not a rate.
  - Line 105: node_rate_rad_day is negative (−2π / period). The regression
    of the lunar ascending node is a well-known astronomical fact (18.6-year
    period). The sign is correct.
  - The function does not validate inputs (e.g., periods > 0); validation is
    delegated to ConstantsProvider construction (REQ-CP-3).
  - Line 103: lunar_node_epoch is stored; this is the node longitude at the
    reference epoch.
```

---

## File-level Verdict

- **A. Error wiring**: ✓ All six return fields (`mean_anomaly_rate_rad_min`,
  `longitude_rate_rad_day`, `orbit_eccentricity`, `orbit_inclination`,
  `node_rate_rad_day`, etc.) are `TrackedValue<T>` instances populated by
  arithmetic operations that automatically propagate input errors per REQ-EF-3.
  Conform to AUD-EF-1 and AUD-EF-2.

- **B. Algebra axioms**: n/a (constructors, not algebra operations).

- **C. Theoretical basis**:
  - Card 1 (`make_sun`): ✓ All fields are linear functions of published
    astronomical constants (period, eccentricity, etc.). Closed-form
    assignment with no truncation. Theory (Almanac / SR3) matches method
    (assignment) matches bound (error propagation per REQ-EF-3).
    **PASS.**
  - Card 2 (`make_moon`): ✓ All fields are linear functions of published
    lunar orbital constants. Closed-form assignment with no truncation.
    Theory (Almanac / SR3) matches method (assignment) matches bound (error
    propagation per REQ-EF-3). Special note: mean_anomaly_epoch set to
    exact(0) per SGP4 convention; this is correct and audited.
    **PASS.**

**File verdict: PASS** — Both factory methods are closed-form assignments of
defined and measured constants. Error propagation is automatic and correct per
REQ-EF-3. No method-theory mismatch; no unsound bounds.

---

## References

1. Hoots, F. R., & Roehrich, R. L. (1980). Spacetrack Report No. 3: Models for
   Propagating the NORAD Element Set. United States Air Force Aerospace Defense
   Command (later corrected editions: 1988, 2004).
2. Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T. S. (2006). Revisiting
   Spacetrack Report #3: Rev2. In *Proceedings of the AIAA/AAS Astrodynamics
   Specialist Conference*. AIAA Paper 2006-6753.
3. The Astronomical Almanac for the Year [annual]. U.S. Naval Observatory /
   Her Majesty's Nautical Almanac Office.
4. Chapront-Touzé, M., & Chapront, J. (1988). ELP 2000-85: A semi-analytical
   lunar ephemeris adequate for historical times. *Astronomy and Astrophysics*,
   190, 342–352.
