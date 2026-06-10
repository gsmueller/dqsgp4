# Theoretical Basis Audit — `src/forces/drag.h`

## File summary

Single public function: `make_drag_exponential(rho_0, h_0, H_scale, B)` returns a drag wrench lambda.

Formula count: 1 (exponential-atmosphere drag acceleration).

---

## FORMULA AUDIT CARD

```
=== FORMULA AUDIT CARD ===
ID:                     drag::make_drag_exponential
Location:               src/forces/drag.h:77-136
Mathematical statement: Drag acceleration a_drag = −(1/2) ρ(h) B |v_rel| v_rel
                        where ρ(h) = ρ_0 exp(−(h − h_0) / H_scale),
                        h = |r| − R_E (geocentric altitude),
                        v_rel = v_world − (ω_E × r) (relative velocity)

THEORY
  Underlying theorem:   Exponential atmosphere model (Lane 1965). The model
                        asserts that atmospheric density decays exponentially
                        with altitude. For drag force, the formula follows
                        from the fundamental drag equation: the coefficient
                        of restitution is folded into the ballistic parameter
                        B = C_d · A / m. The velocity-squared dependence
                        arises from continuum aerodynamics and the
                        assumption of subsonic flow (M < 0.3 for LEO).
  Primary reference:    Lane, M. L. (1965) "An improved atmospheric drag
                        model for the SGP orbit prediction model." CAECOM
                        report. Vallado et al. (2006) "Revisiting Spacetrack
                        Report #3" §4.3 (exponential density model).
  Domain of validity:   LEO altitudes (h = 300–2000 km typical); the
                        exponential model has ~50 % error relative to MSIS
                        due to diurnal, seasonal, and solar-cycle variation.

METHOD
  Method declared:      Closed-form exponential density ρ = ρ_0 exp(arg),
                        multiplied through with v_rel to compute the
                        drag-acceleration vector.
  Method implemented:   Lines 107–129: (1) compute h directly from r_mag,
                        (2) compute arg = −(h − h_0)/H_scale in raw T,
                        (3) call std::exp(arg) in raw T, (4) scale prefactor
                        by ρ, (5) multiply prefactor into v_rel components.
  Match verdict:        ✓ matched — implementation is the closed-form
                        exponential density multiplied with the velocity
                        vector; no approximation method substituted.

ERROR BOUND
  Bound category:       accuracy (model truncation) + precision
  Bound formula:        Model accuracy (REQ-EF-7): The exponential model is a
                        rank-2 truncation of the true density profile. Vallado
                        et al. (2006) report ~30 % error vs U.S. Standard
                        Atmosphere; NRLMSISE-00 shows ~50 % error in the
                        exponential-model residual. The relative error in drag
                        force is ≤ c_model ≈ 0.5. Code adds this to prefactor's
                        accuracy bound (line 119):
                          accuracy_bound = |prefactor_val| · 0.5
                        Precision: std::exp is a C library function, not a
                        TrackedValue<T> primitive. Representation cost is
                        captured via representation_bound(·) on the prefactor.
  Bound implemented:    Lines 119–124. TrackedValue<T> prefactor with
                        accuracy_bound = |prefactor_val| · 0.5 and precision
                        component from representation_bound(prefactor_val).
  Bound verdict:        ⚠ partially sound — the model-accuracy bound (50 %)
                        is plausible per Vallado et al.; however, it is a
                        global constant, not state-dependent (ignoring h,
                        season, solar activity). Correct in order of magnitude
                        for LEO but conservative for low h and potentially
                        non-conservative for high h (>800 km). Precision cost
                        from std::exp is captured indirectly, not decomposed.
                        **FLAG for refinement: state-dependent accuracy bound
                        and decomposed precision bounds.**

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-7 (model truncation → accuracy)
  AUD-EF applies:       AUD-EF-6 (accuracy budget on wrench for model
                        simplification), AUD-EF-1 (all public ops return
                        TrackedValue<T>)
  AUD-MC applies:       n/a (force lambda, not an algebra op)
  Verification test:    Compare drag wrench from make_drag_exponential with
                        MSIS-computed wrench for a range of LEO altitudes;
                        verify reported accuracy bound dominates measured
                        error at least 95 % of the time. Current propagator
                        integration tests do not isolate drag accuracy.

NOTES
  - The density ρ = ρ_0 exp(...) is computed in raw T (not TrackedValue<T>)
    because std::exp is not yet a TrackedValue<T> primitive (lines 107–110).
    The representation cost is folded into representation_bound on the
    prefactor, which is sound but non-transparent.
  - The model-accuracy bound c_model = 0.5 is global, independent of altitude,
    season (diurnal variation), and solar activity (10-year cycle). A
    production implementation should accept a state-dependent accuracy callback
    or table to tighten the bound for specific mission constraints.
  - Cross-products and magnitudes (v_atm = ω_E × r, v_rel_mag, r_mag) are
    TrackedValue<T> ops; their error contributions flow into the final
    prefactor via representation_bound but are not decomposed. The bound is
    rigorous but conservative.
  - The body-frame rotation (lines 132–134) is exact (Hamilton product on
    TrackedValue<T>); no error is introduced.
  - Future work: support a model-selector callback to plug in MSIS,
    NRLMSISE-00, or user-provided tabulated density model.
```

---

## File-level verdict

- **A. Error wiring**: ✓ The returned `Wrench<T>` carries the drag acceleration in the linear slot with both `errors.precision` and `errors.accuracy` populated (lines 120–124).
- **B. Algebra axioms**: n/a (force lambda, not an algebra operation).
- **C. Theoretical basis**:
  - `make_drag_exponential` (exponential-atmosphere drag): ✓ Method matches cited theory (Lane 1965 closed-form exponential density). ⚠ Accuracy bound is global (0.5) rather than state-dependent; precision decomposition is implicit. **PASS with note for state-dependent refinement.**

**File verdict: PASS** — Exponential-atmosphere formula correctly implements Lane's closed-form model. Model-truncation accuracy is bounded per REQ-EF-7. Precision cost is captured implicitly. Recommendation: state-dependent accuracy bound and model-selector callback for empirical density models (future work).

---

## References

1. Lane, M. L. (1965). *An improved atmospheric drag model for the SGP orbit prediction model*. CAECOM report.
2. Vallado, D. A., Crawford, P., Hujsa, R., & Johnston, K. (2006). *Revisiting Spacetrack Report #3: Rev 1*. AIAA/AAS Astrodynamics Specialist Conference. Section 4.3 on exponential density.
3. King-Hele, D. G. (1964). *Theory of Satellite Orbits in an Atmosphere*. Butterworths, London.
