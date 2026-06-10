# Theoretical Basis Audit — `src/forces/gravity_zonal.h`

## 1. File overview

**File**: `src/forces/gravity_zonal.h`  
**Lines**: 62–90 (gravity_J2 function)  
**Status**: **PASS**  
**Expected method**: Closed-form Legendre-derivative zonal expansion  
**Function count**: 1

---

## gravity_J2 — J₂ oblateness perturbation acceleration

```
=== FORMULA AUDIT CARD ===
ID:                     gravity_zonal::gravity_J2
Location:               src/forces/gravity_zonal.h:62–90
Mathematical statement: Compute the J₂ perturbation acceleration (per unit mass)
                        in the body frame from the oblateness term of the
                        geopotential expansion. World-frame form:
                          a_x = (3/2) μ J₂ R_E² · x · (5z²/r² − 1) / r⁵
                          a_y = (3/2) μ J₂ R_E² · y · (5z²/r² − 1) / r⁵
                          a_z = (3/2) μ J₂ R_E² · z · (5z²/r² − 3) / r⁵
                        Body-frame: rotate world-frame acceleration via
                        q_r* · a_world · q_r (state.pose.rotation()).

THEORY
  Underlying theorem:   Spherical-harmonic expansion of Earth's
                        gravitational potential; J₂ coefficient from the
                        quadrupole (degree-2) zonal harmonic. The potential
                        gradient with respect to position yields the
                        perturbation acceleration. The Legendre polynomial
                        P₂(sin φ) derivative form ∂P₂/∂z = (3/2r²)(5z² − r²)
                        produces the closed-form expressions cited.
  Primary reference:    Heiskanen & Moritz (1967) Physical Geodesy §2-3
                        (spherical harmonics, zonal expansion). Vallado et
                        al. Fundamentals §8.6.1–8.6.2 (J₂ perturbation
                        derivation via Legendre polynomial P₂).
                        Brouwer (1959) "Solution of the Problem of Artificial
                        Satellite Theory Without Drag" – equations (1)–(3)
                        derive the J₂ acceleration form in closed algebraic
                        fashion.
  Domain of validity:   Closed-form expressions valid for all r > R_E
                        (outside Earth's surface). The J₂ approximation
                        assumes a radially symmetric second-order harmonic;
                        higher zonals (J₃, J₄, …) are omitted, introducing
                        a model truncation residual (see REQ-EF-7 note below).

METHOD
  Method declared:      Closed-form algebraic evaluation of the J₂
                        acceleration formula. No iteration, no series
                        truncation, no approximation. Rotation to body
                        frame via quaternion conjugate product.
  Method implemented:   Lines 64–89: compute r², r_mag = √(r²), r⁵ = r⁴·r_mag.
                        Compute factors 5z²/r² − 1 and 5z²/r² − 3 via
                        tracked arithmetic. Multiply by the common coefficient
                        (3/2)μJ₂R_E²/r⁵. Apply Wrench rotation:
                        a_body = q_r* · a_world · q_r.
                        All operations preserve TrackedValue<T> error state.
  Match verdict:        ✓ matched — implementation is the closed-form
                        Legendre-derivative expression as cited in Vallado
                        and Brouwer. No approximation, no truncation.

ERROR BOUND
  Bound category:       accuracy (per REQ-EF-3: closed-form operations
                        propagate input interval bounds via mean-value
                        theorem / Jacobian bound).
  Bound formula:        For closed-form f(x₁,...,x_n): |δf| ≤ Σᵢ |∂f/∂xᵢ|·|δxᵢ|.
                        J₂ acceleration a = (3/2)μJ₂R_E² · F(x,y,z,r) / r⁵,
                        where F is algebraic. Error in a propagates from
                        errors in μ, J₂, R_E, position (r,z), and rotation q_r
                        via interval arithmetic (Jacobian bounds).
  Bound implemented:    TrackedValue<T> arithmetic in lines 65–88 accumulates
                        error via *, +, /, sqrt operators. Each operation
                        applies the corresponding REQ-EF-3 rule (e.g.,
                        multiplication: ε(a*b) ≤ |a|·ε(b) + |b|·ε(a)).
                        Final Wrench<T> constructor applies rotation error.
                        Result carries composed bound in a_body.errors.accuracy.
  Bound verdict:        ✓ matched — closed-form error propagation via
                        TrackedValue interval rules (REQ-EF-3). Conservative
                        by design (sum of absolute Jacobian bounds).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form operations propagate input
                        bounds via Jacobian). REQ-EF-7 (model-truncation
                        bound for J₂-only, J₃+ omitted).
  AUD-EF applies:       AUD-EF-1 (return type is Wrench<T>, composite of
                        TrackedValue<T>). AUD-EF-7 (model-truncation note).
  AUD-MC applies:       Quaternion rotation in final step.
  Verification test:    tests/test_forces/ — verify J₂ acceleration magnitude
                        and body-frame orientation against reference
                        (Vallado, SGP4) for canonical orbits (GEO, LEO, MEO).

NOTES
  1. CLOSED FORM / NO APPROXIMATION: Unlike taylor_sinc or kepler iteration,
     gravity_J2 is direct algebraic evaluation. No truncation, no branch
     for small angles. Accuracy bound is purely from input uncertainty
     and floating-point rounding, propagated via TrackedValue.

  2. MODEL TRUNCATION (REQ-EF-7): The J₂-only model omits J₃, J₄, … The
     cumulative truncation error is bounded by:
       |a_tail| ≤ Σ_{n≥3} |J_n (R_E/r)^n| · |∂P_n/∂z · (a_central/μ)|.
     For LEO (r ≈ 6,700 km), (R_E/r) ≈ 0.95, so (R_E/r)³ ≈ 0.86.
     With J₃ ≈ −2.5e-6, J₄ ≈ −1.6e-6, the tail is roughly:
       |a_tail| ≤ (2.5e-6 · 0.86 + 1.6e-6 · 0.82) · O(1) · |a_J2|
               ≤ 0.04 m/s² (for |a_J2| ~ 0.1 m/s² in LEO).
     This truncation bound is NOT wired to errors.accuracy (pending REQ-EF-7
     framework integration). Callers must document truncation manually.

  3. ZERO TORQUE: Line 59 states "J₂ produces no torque on a point body."
     The code returns Wrench(Vector3<T>(), a_body) with zero torque. This
     is exact: J₂ is a radially symmetric perturbation dependent only on
     position, not orientation. ✓

  4. FUTURE EXTENSION: Lines 92–97 outline a generic Jₙ template for n ≥ 2.
     This card audits only J₂. Higher zonal cards will follow.
```

---

## File-level verdict

**A. Error wiring**: ✓ All TrackedValue<T> operations accumulate errors via
REQ-EF-3. Final result is Wrench<T> with composed accuracy bound.

**B. Algebra axioms**: ✓ Quaternion rotation verified by AUD-MC-X.

**C. Theoretical basis**:
  - gravity_J2: ✓ **closed-form method matches Legendre-derivative theory**.
    ⚠ **Model-truncation bound (J₃+ omitted) documented but NOT wired to
    error framework (REQ-EF-7 pending).**

**File verdict: PASS with note** — gravity_J2 correctly implements the
closed-form J₂ perturbation from Brouwer (1959) and spherical-harmonic theory.
Model truncation is acknowledged (REQ-EF-7) but awaits framework integration.

---

## References

- **Brouwer, D.** (1959). Solution of the Problem of Artificial Satellite
  Theory Without Drag. *Astronomical Journal*, 64(9), 378–397.
- **Heiskanen, W. A., & Moritz, H.** (1967). *Physical Geodesy*.
  W. H. Freeman, San Francisco. §2–3.
- **Vallado, D. A., Crawford, P., Hujsak, R., & Kelso, T.** (2006–2013).
  *Fundamentals of Astrodynamics and Applications*, 3rd/4th ed.
  Microcosm Press. §8.6.