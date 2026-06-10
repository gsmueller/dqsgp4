# Theoretical Basis Audit — `src/forces/gravity_central.h`

## Overview

File `src/forces/gravity_central.h` contains two functions:
1. **`gravity_central(state, K)`** — Computes body-frame gravitational wrench from point-mass attractor.
2. **`make_gravity_central(K)`** — Factory function; returns a lambda closure wrapping `gravity_central`.

This file implements the **central (spherically symmetric) part** of Earth's gravitational field. Model truncation (zonal, tesseral, or full harmonics) is delegated to separate force lambdas; the accuracy loss from this truncation is added to `errors.accuracy` per REQ-EF-7.

---

## Card 1: Core Formula — `gravity_central::central_force`

```
=== FORMULA AUDIT CARD ===
ID:                     gravity_central::central_force
Location:               src/forces/gravity_central.h:51-70
Mathematical statement: F_world = −(GM / |r|³) · r,
                        where r ∈ ℝ³ is the body's position in the world frame,
                        then rotated to body frame: F_body = q_r* · F_world · q_r.

THEORY
  Underlying theorem:   Newton's law of universal gravitation (1687, Principia 
                        Book III, Prop. VII): The gravitational force exerted by 
                        a point mass M on a mass m at distance r is
                          F = −(GM / r²) · r̂,
                        where r̂ = r / |r| is the unit vector. Equivalently,
                          F = −(GM / |r|³) · r.
                        This is a closed-form, exact relationship under the 
                        assumptions of Newtonian mechanics (non-relativistic, 
                        point-mass approximation).
  Primary reference:    Newton, I. (1687). Philosophiae Naturalis Principia 
                        Mathematica. Book III, Prop. VII.
                        Modern treatment: Vallado, Crawford & Hujsak (2006) 
                        "Revisiting Spacetrack Report #3", §1, Two-Body Problem.
  Domain of validity:   All r ∈ ℝ³ \ {0}. The formula is singular at r = 0 
                        (undefined physically; orbital dynamics avoids this by 
                        requiring r_min ~ 6,371 km for Earth). No branch cuts 
                        or removable singularities in the valid domain.

METHOD
  Method declared:      Closed-form algebraic evaluation of the inverse-square 
                        law. No Taylor expansion, no iterative root-finding, 
                        no series truncation.
  Method implemented:   Lines 54–64:
                          r_sq = x² + y² + z²
                          r_mag = sqrt(r_sq)
                          r_cubed = r_mag · r_sq
                          F_world = −(GM / r_cubed) · [x, y, z]
                        Lines 68–69: Rotate F_world to body frame via 
                        quaternion conjugate (q_r* · F_world · q_r).
  Match verdict:        ✓ matched — implementation is the closed-form 
                        inverse-square law; no approximations.

ERROR BOUND
  Bound category:       precision (REQ-EF-3 propagation) and accuracy 
                        (REQ-EF-7 model truncation, handled separately).
  Bound formula:        Error propagation via composed TrackedValue<T> operations:
                          1. x² + y² + z²         — per-operation error via 
                                                    multiplication and addition
                          2. sqrt(r_sq)           — sqrt error via REQ-EF-3
                          3. r_mag · r_sq         — multiplication error
                          4. GM / r_cubed         — division error
                          5. Scalar · Vector      — componentwise scaling error
                          6. Quaternion rotation  — rotation-composition error 
                                                    (via quaternion algebra axioms)
                        Each step uses the mean-value theorem bound for the 
                        composed operation: e.g., for f(a) / g(b), 
                        δ(f/g) ≤ |δf/g + f·δg/g²|. The final total_error() is 
                        the sum of these per-operation contributions.
  Bound implemented:    Implicit in the TrackedValue<T> semantics: each 
                        `operator*`, `operator/`, `operator+`, `sqrt()`, and 
                        quaternion `rotate()` updates the TrackedValue's error 
                        fields (per AUD-EF-1 through AUD-EF-7). The returned 
                        Wrench<T> carries the composed error.
  Bound verdict:        ✓ matched — error propagation is via REQ-EF-3 
                        closed-form composition rules, appropriate for a 
                        closed-form algebraic formula with no truncation.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form algebraic composition), 
                        REQ-EF-7 (model truncation: central gravity is the 
                        leading-order term; J₂ and higher are added separately).
  AUD-EF applies:       AUD-EF-1 (all intermediate results are TrackedValue<T>), 
                        AUD-EF-3 (operations return TrackedValue),
                        AUD-EF-7 (wrench assembly returns Wrench<T> with 
                        consolidated error).
  AUD-MC applies:       AUD-MC-N (quaternion rotation axioms, once quaternion.h 
                        is audited — rotation distributes over vector scaling, 
                        composition is associative).
  Verification test:    tests/test_forces/test_gravity_central.cc — 
                        compareAgainstHighPrecision() at known orbital radii 
                        (LEO r ~ 6.6e6 m, GEO r ~ 4.2e7 m) and verify that 
                        reported total_error() dominates actual error.

NOTES
  - The formula is **not an approximation**; it is exact (within the 
    assumptions of Newtonian mechanics and point-mass approximation). 
    Error is not truncation error; it is machine-precision loss from 
    finite arithmetic.
  - Model truncation (J₂, J₃, …, tesseral, sectorial harmonics) is NOT 
    included here. The contribution from omitted terms is added to 
    errors.accuracy by separate force lambdas (gravity_zonal.h, etc.), 
    implementing REQ-EF-7.
  - **Rotation composition**: The quaternion rotation q_r* · F_world · q_r 
    assumes q_r is unit (normalized). The input state.pose.rotation() is 
    typically normalized; any denormalization error is carried in the 
    quaternion's TrackedValue. This is correct but assumes quaternion 
    normalization is maintained upstream (e.g., by the propagator).
```

---

## Card 2: Factory Function — `make_gravity_central`

```
=== FORMULA AUDIT CARD ===
ID:                     gravity_central::make_gravity_central
Location:               src/forces/gravity_central.h:77–82
Mathematical statement: Returns a function object (lambda) that captures K 
                        (constants provider) and delegates to 
                        gravity_central(state, K).

THEORY
  Underlying theorem:   No new numeric formula. This is a **structural wrapper** 
                        that captures state (the constants provider K) and 
                        returns a closure suitable for the propagator's 
                        force-list interface. The numeric computation is 
                        entirely delegated to gravity_central, which implements 
                        Newton's law (as in Card 1).
  Primary reference:    (same as Card 1: Newton's gravitation)
  Domain of validity:   Same as gravity_central. The factory takes K (const 
                        ref) and ensures it outlives the returned lambda. 
                        Caller responsibility.

METHOD
  Method declared:      Factory pattern: capture + delegate. No approximation, 
                        no iteration.
  Method implemented:   Lines 79–81: lambda captures K by const ref, calls 
                        gravity_central(state, K). Identity function over 
                        gravity_central's return.
  Match verdict:        ✓ matched — wrapper is transparent.

ERROR BOUND
  Bound category:       (inherited from gravity_central)
  Bound formula:        make_gravity_central does not modify the error report. 
                        It returns a callable that, when invoked, produces the 
                        same TrackedValue<T> and Wrench<T> as gravity_central 
                        directly, with the same error audit trail.
  Bound implemented:    Implicit: the returned lambda is transparent; calling 
                        it invokes gravity_central, which propagates errors 
                        per Card 1.
  Bound verdict:        ✓ matched — no error manipulation.

CROSS-AUDIT
  REQ-EF applies:       (inherited from gravity_central)
  AUD-EF applies:       AUD-EF-1, AUD-EF-3, AUD-EF-7 (via the returned lambda 
                        and its calls to gravity_central).
  AUD-MC applies:       (inherited from gravity_central)
  Verification test:    tests/test_forces/ — verify that the lambda returned 
                        by make_gravity_central(K) produces identical output 
                        to a direct call to gravity_central(state, K).

NOTES
  - This is a **factory, not a formula**. Its audit is trivial: it is a 
    transparent wrapper. The real audit is on gravity_central.
  - Caller must ensure K outlives the lambda. This is a lifetime contract, 
    not an error-budget issue. Code review (AUD-CC) should verify this 
    assumption at each callsite.
```

---

## File-level verdict

- **A. Error wiring**: ✓ All intermediate calculations in gravity_central use 
  TrackedValue<T> (lines 54–63); the returned Wrench<T> carries the composed 
  error. Conforms to AUD-EF-1, AUD-EF-3, AUD-EF-7.

- **B. Algebra axioms**: ✓ Quaternion rotation in line 69 assumes q_r is a 
  unit quaternion (or carries its denormalization error in its TrackedValue). 
  Composition with vector scalar multiplication distributes correctly 
  (once quaternion.h passes AUD-MC-12 or similar).

- **C. Theoretical basis**:
  - Card 1 (gravity_central): ✓ Method matches cited theory (Newton's 
    closed-form law). Bound is by per-operation error propagation per 
    REQ-EF-3. **PASS.**
  - Card 2 (make_gravity_central): ✓ Transparent wrapper; no independent 
    formula to audit. **PASS.**

**File verdict: PASS** — Central gravity is correctly implemented as the 
Newton inverse-square law with error propagation via composed TrackedValue 
operations. Model-truncation accuracy (J₂, …) is handled by separate audits.

---

## References

- Newton, I. (1687). *Philosophiae Naturalis Principia Mathematica.* Book III, 
  Proposition VII.
- Vallado, Crawford & Hujsak (2006). "Revisiting Spacetrack Report #3," 
  *Journal of Guidance, Control, and Dynamics* 29(2). (Two-body problem 
  background.)
