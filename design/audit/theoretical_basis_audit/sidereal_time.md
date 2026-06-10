# Theoretical Basis Audit — `src/astronomy/sidereal_time.h`

## Overview

File: `src/astronomy/sidereal_time.h`  
Function count: 1 (`compute_gmst`)  
Formulas: 1 (Aoki et al. 1982 GMST cubic polynomial + sidereal/solar rotation + radian conversion)  
Status: PASS (closed-form Horner expansion in T; measurement bounds wired via `TrackedValue<T>`; precision bounds propagated via the closed-form composition rule)

The function computes Greenwich Mean Sidereal Time (GMST) — the hour angle of the mean vernal equinox at Greenwich — used for the TEME→ECEF frame transformation. The polynomial is the IAU 1982 / IERS-adopted Aoki form, evaluated by Horner's scheme on Julian centuries since J2000.0, with the intra-day rotation added via the sidereal/solar ratio. Final wrap to `[0, 2π)` is a closed-form `fmod` identity.

---

## Formula Audit Cards

### Card 1: Aoki (1982) GMST Polynomial — Closed-Form Horner Expansion

```
=== FORMULA AUDIT CARD ===
ID:                     sidereal_time::compute_gmst
Location:               src/astronomy/sidereal_time.h:48-93
Mathematical statement: GMST(JD_UT1) is computed as
                          T_0h   = (JD_0h − 2451545.0) / 36525
                          f_day  = JD_UT1 − JD_0h
                          s_0h   = c₀ + c₁ T_0h + c₂ T_0h² + c₃ T_0h³   [sidereal s]
                          s_full = s_0h + r · (f_day · 86400)            [sidereal s]
                          θ      = (s_full · 2π / 86400) mod 2π          [rad ∈ [0, 2π)]
                        with
                          c₀ = 24110.54841
                          c₁ = 8640184.812866
                          c₂ = 0.093104
                          c₃ = −6.2 × 10⁻⁶
                          r  = 1.00273790935            (sidereal/solar day ratio)
                          JD_0h = floor(JD_UT1 − 0.5) + 0.5
                        Horner form actually evaluated:
                          s_0h = c₀ + T_0h · (c₁ + T_0h · (c₂ + T_0h · c₃))

THEORY
  Underlying theorem:   Closed-form polynomial identity (Horner's scheme is
                        an exact rearrangement, not an approximation).
                        The polynomial itself is the IAU 1982 expression
                        derived by Aoki et al. (1982) from VLBI and lunar-
                        laser-ranging observations of Earth orientation;
                        it is a fit, not a theorem, so its coefficients
                        carry measurement uncertainty rather than truncation
                        error. The intra-day extension uses the closed-form
                        identity
                          GMST(t) = GMST(0h) + r · (t − 0h) [seconds],
                        where r = (sidereal day)/(solar day) is exact per
                        the IAU/IERS definition of mean sidereal time.
                        The wrap-to-[0, 2π) is the closed-form mod-2π
                        identity on the circle.
  Primary reference:    - Aoki, S. et al. (1982), "The new definition of
                          Universal Time", Astronomy & Astrophysics 105,
                          pp. 359-361. (Source of the four-coefficient
                          polynomial and 1.00273790935 ratio.)
                        - IERS Conventions (2010), Petit & Luzum eds.,
                          IERS Technical Note 36, Chapter 5.5
                          ("Sidereal time and UT1").
                        - Vallado, D. (2013), Fundamentals of Astrodynamics
                          and Applications, 4th ed., §3.5 "Sidereal Time",
                          Algorithm 15 (THETAG / GSTIME).
                        - Capitaine, N. et al. (2003), "Expressions for IAU
                          2000 precession quantities", A&A 412, pp 567-586,
                          for the coefficient measurement uncertainties.
                        - design/derivations/014_sidereal_time.md
  Domain of validity:   All JD_UT1 in the modern era (Aoki polynomial is
                        formally a Taylor-style series in T but is treated
                        as exact closed form within ~century-scale T;
                        IERS recommends |T| ≲ 1 century from J2000.0 to
                        keep coefficient measurement error dominant over
                        omitted higher-order terms).

METHOD
  Method declared:      Closed-form Horner evaluation of the IAU 1982
                        cubic polynomial in T_0h, plus closed-form linear
                        addition of intra-day rotation, plus closed-form
                        scalar conversion seconds → radians, plus a
                        closed-form wrap-to-[0, 2π). No series truncation,
                        no iteration, no Padé/continued-fraction or
                        other rational approximant.
  Method implemented:   Lines 54-65: closed-form arithmetic
                          T_0h    = (JD_0h - 2451545) / 36525  (line 64)
                          frac_day = JD_UT1 - JD_0h            (line 65)
                        Line 77: Horner evaluation
                          gmst_0h_sec = c0 + t_0h*(c1 + t_0h*(c2 + t_0h*c3))
                        Lines 82-85: closed-form linear addition
                          gmst_sec = gmst_0h_sec + r * (frac_day * 86400)
                        Lines 88-89: closed-form scalar conversion
                          theta = gmst_sec * (2π / 86400)
                        Line 92: closed-form wrap
                          return wrap_two_pi(theta)
                        Every arithmetic step is `TrackedValue<T>` operator
                        overload; no loops, no iteration, no series
                        truncation in code path.
  Match verdict:        ✓ matched — implementation is a literal Horner
                        evaluation of the cited Aoki polynomial. Not a
                        Taylor truncation of some other function, not a
                        Padé approximant, not a continued fraction, not
                        an iterative solve.

ERROR BOUND
  Bound category:       measurement (dominant) + precision (subordinate)
  Bound formula:        Two distinct contributions, both rigorous under
                        REQ-EF-3 (closed-form identity propagation):

                        (1) MEASUREMENT (Aoki coefficient uncertainty).
                            Each c_k is a `measured` TrackedValue with
                            sigma σ_k from VLBI re-determinations. The
                            polynomial Σ c_k T_0h^k propagates these into
                            the output by closed-form scalar product /
                            sum, giving
                              σ_meas(s_0h) ≤ Σ_k σ_k |T_0h|^k.
                            Conversion to radians multiplies by 2π/86400,
                            then wrap_two_pi is sign-preserving on
                            magnitudes. The TrackedValue composition
                            performs this exactly under operator
                            overloading.

                        (2) PRECISION (input JD_UT1 propagation + machine
                            rounding). Because the operation is closed-form
                            polynomial, by mean-value theorem the
                            sensitivity of the output to JD_UT1 is its
                            derivative
                              dθ/d(JD_UT1) = (2π/86400) · (d s_full/d JD_UT1)
                                           = (2π/86400) · (c₁/36525 + 2 c₂ T/36525
                                                          + 3 c₃ T²/36525 + r·86400)
                            At T ≈ 0 this is dominated by the r·86400 ·
                            2π/86400 = 2π · r term ≈ 6.30 rad/day —
                            i.e. Earth's rotation rate. Closed-form
                            multiplicative propagation gives the precision
                            bound by REQ-EF-3.

                        No truncation bound is needed: the chosen method
                        is closed-form, not a Taylor / Padé / series cut.
  Bound implemented:    All four coefficients constructed via
                        `TrackedValue<T>::measured(value, sigma)` (lines
                        72-75), which seeds the `errors.measurement`
                        category with the stated σ. The `defined` constants
                        2451545.0, 36525.0, 1.00273790935, 86400.0 (lines
                        54-55, 82-83) carry zero error per REQ-EF-1.
                        The composition through arithmetic / Horner /
                        multiply / wrap propagates `measurement` and
                        `precision` through `TrackedValue<T>` operator
                        overloads per REQ-EF-3.
                        `floor(...)` at line 62 is integer-stable for the
                        modern epoch (representable exactly in `T`); the
                        constructed `jd_0h` inherits `jd_ut1.errors` so
                        the input precision is not dropped (line 63).
  Bound verdict:        ✓ matched — closed-form polynomial composition;
                        measurement bounds are seeded at the source
                        constants and propagated by TrackedValue's
                        closed-form rules; no truncation bound needed
                        because no method-level truncation occurs.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (defined constants carry zero error)
                        REQ-EF-3 (closed-form identity error propagation)
                        REQ-EF-2 indirectly (the wrap_two_pi closed form
                          is sign-preserving on magnitudes)
  AUD-EF applies:       AUD-EF-2 (algebraic operations add errors correctly,
                          exercised by every + and * in the Horner chain)
  AUD-MC applies:       n/a — `compute_gmst` consumes the algebra; it is
                        not itself an algebra operation. Algebra
                        correctness is audited under AUD-MC-1..18.
  Verification test:    tests/test_astronomy/ (sidereal-time tests, if
                        present) — verify GMST at J2000.0 (Greenwich noon,
                        2000-Jan-01 12:00 UT1) is 18h 41m 50.54841s
                        (≡ 4.894961... rad), and that GMST(t+86400s_UT1)
                        − GMST(t) = 2π · 1.00273790935 rad (mod 2π).
                        Compare against Vallado (2013) Table 3-7 reference
                        values or against an independent IERS-conformant
                        implementation (e.g. SOFA `iauGmst82`).

NOTES
  - The Aoki polynomial is not a Taylor expansion of any analytic function:
    it is an empirical fit to Earth-orientation observations. The
    "coefficients" c_k therefore carry measurement σ rather than
    truncation σ, which is exactly what the code reflects by using
    `TrackedValue<T>::measured(...)` rather than `::defined(...)`.
    The Bound category split is consequently measurement (dominant in
    absolute terms) + precision (dominant in derivative terms near
    one day).
  - Higher-order terms (T⁴ and beyond) in IERS 2010 (Chapter 5.5,
    expression for GMST 2003) are absorbed into the c_k measurement σ
    here. For multi-century propagation those higher terms become
    non-negligible; for SGP4 / modern operational ranges (|T| ≲ 1
    century) they are absorbed.
  - The `jd_0h` split (line 62) avoids catastrophic cancellation by
    keeping the polynomial argument T_0h independent of the intra-day
    fraction. This is the standard Vallado §3.5 trick and is
    independent of the polynomial's theory; it improves numerical
    *precision*, not the rigor of the bound.
  - `1.00273790935` is itself a measured quantity (the ratio of mean
    sidereal day to mean solar day), but is carried as `defined` here
    per the IAU/IERS convention that fixes its value as the operational
    definition of "mean sidereal time". Any uncertainty in r is
    absorbed into the c_k coefficient uncertainties of the Aoki fit.
    If a future audit requires r to carry σ, the change is local
    (line 82) and propagates automatically via TrackedValue.
  - `wrap_two_pi` (line 92, from `math/angles.h`) is a closed-form
    `fmod`-based identity; its REQ-EF-3 bound is propagated through
    angle-arithmetic preserving magnitudes (no truncation, no
    method-level error). Audited under `theoretical_basis_audit/angles.md`.
```

---

## File-Level Verdict

**Error wiring (REQ-EF / AUD-EF)**:  
✓ The function returns `TrackedValue<T>`. All four Aoki coefficients are seeded via `::measured(...)` so their σ enters the measurement category; defined constants are seeded via `::defined(...)` with zero error per REQ-EF-1. Every arithmetic step is a TrackedValue operator overload, propagating measurement and precision per REQ-EF-3. Input `jd_ut1.errors` is preserved through the `jd_0h` reconstruction (line 63) — no dropped error.

**Algebra axioms (AUD-MC)**:  
n/a — `compute_gmst` is a closed-form scalar computation, not itself an algebra operation. The underlying `TrackedValue<T>` algebra correctness is audited under AUD-MC-1..18.

**Theoretical basis (TBA)**:  
- **Card 1 (Aoki GMST polynomial)**: ✓ Closed-form Horner evaluation of the IAU 1982 polynomial + closed-form linear sidereal-rotation extension + closed-form scalar radian conversion + closed-form mod-2π wrap. No Taylor / Padé / continued-fraction / iterative substitution anywhere in the chain. Measurement bounds are seeded at source; precision is propagated by closed-form composition. **Method matches cited theory; bound matches the rigorous closed-form-propagation formula.**

**Overall verdict: PASS**

The function is a closed-form evaluation of a published empirical polynomial. The theory (IAU 1982 Aoki polynomial, IERS-adopted sidereal/solar ratio, mod-2π wrap) is matched by the method (Horner expansion + scalar composition + `wrap_two_pi`); the error bound is the rigorous closed-form-propagation bound under REQ-EF-3, with the Aoki coefficient σ entering correctly via `::measured(...)` constructors. No C-fail.

---

## Cross-References

**Uses (backward)**:  
- `math::TrackedValue<T>` — error propagation (AUD-EF-2, REQ-EF-3, REQ-EF-1)  
- `math::two_pi<T>()` — closed-form constant (REQ-EF-1)  
- `math::wrap_two_pi(...)` — closed-form mod-2π identity (audit card: `theoretical_basis_audit/angles.md`)  
- Aoki et al. (1982), IERS Conventions (2010) §5.5, Vallado §3.5  

**Feeds (forward)**:  
- TEME → ECEF frame rotation in any propagator stage that outputs Earth-fixed coordinates (SGP4 output transform; ground-track computation).  
- Any consumer of Earth orientation requiring 0.1 ms-equivalent-angle accuracy on operational timescales (|T| ≲ 1 century from J2000.0).
