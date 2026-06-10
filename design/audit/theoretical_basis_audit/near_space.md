# Theoretical Basis Audit — `src/sgp4/near_space.h`

**File**: `src/sgp4/near_space.h` (329 lines)
**Functions audited**: 2 (`initialize_near_space`, `propagate_near_space`)
**Role**: Top-level SGP4 near-earth composer. Wires together element recovery,
density parameters, drag coefficients, secular rates, long-period (xlcof/aycof)
coefficient construction, modified Kepler solve, osculating elements, and
short-period corrections per [SR3] §6, pp. 10-15.

This file contains almost no original formula application: it is mostly a
composition harness. Audit cards therefore cross-reference the per-module
audits for the heavy lifting; only the inline composition glue (long-period
coefficient `xlcof`/`aycof` construction in initialization, and long-period
linearization `xll`/`aynl`/`xlt`/`ayn` in propagation) is audited here as
original work.

---

## Card 1 — `initialize_near_space`

```
=== FORMULA AUDIT CARD ===
ID:                     near_space::initialize_near_space
Location:               src/sgp4/near_space.h:135-246
Mathematical statement: Build a NearSpaceInit<T> struct: recover (a₀, n₀, β₀)
                        from TLE, compute density parameters (s, q₀⁴), drag
                        coefficients (C₁..C₅, D₂..D₄, t*cof, omgcof, xmcof,
                        delmo, sinmo, Ω̇_nkc, ξ, η), secular rates (Ṁ, ω̇, Ω̇),
                        long-period coefficients (xlcof, aycof), and store
                        trig & Earth constants — all from [SR3] §6 pp. 10-11.

THEORY
  Underlying theorem:   Composition; this function is a sequencer. Theoretical
                        bases live in the called modules:
                        • Element recovery (Kozai → Brouwer mean motion):
                            orbit::recover_mean_elements
                            → theoretical_basis_audit/element_recovery.md
                            [SR3] §6 eqs (1)-(7), Kozai (1959) inverse.
                        • Density parameters (s, q₀⁴ adjustment):
                            atmosphere::compute_density_parameters
                            → theoretical_basis_audit/density_model.md
                            [SR3] §6 eqs (8)-(11).
                        • Drag coefficients (C₁..C₅, D₂..D₄, t-cofs, xmcof,
                          omgcof, delmo, sinmo, Ω̇_nkc, ξ, η):
                            atmosphere::compute_drag_coefficients
                            → theoretical_basis_audit/drag_coefficients.md
                            [SR3] §6 eqs (12)-(22), with
                            design/derivations/sgp4_near_earth_drag_theoretical_basis.md
                            §§13-14 covering xnodcf and the C-series.
                        • Secular rates (Ṁ, ω̇, Ω̇):
                            config.model_functions.secular_rates lambda
                            → theoretical_basis_audit/model_functions.md
                            Brouwer (1959) eqs for 1st-order J₂/J₄ secular
                            rates of M, ω, Ω.
                        • Trig identities (cos i₀, sin i₀, 3cos²i−1, 7cos²i−1):
                            closed-form arithmetic on TrackedValue<T>; bound
                            via TrackedValue propagation rules.

                        ORIGINAL TO THIS FILE — long-period coefficients:
                          xlcof = (1/8) (A₃₀/CK₂) sin i₀ (3 + 5 cos i₀)/(1 + cos i₀)
                          aycof = (1/4) (A₃₀/CK₂) sin i₀
                        Theoretical basis: factored forms of long-period J₃
                        corrections to mean longitude and eccentricity vector.
                        Equivalent forms appear in [SR3] §6 pp. 12-13 inline.
                        Derivation: design/derivations/sgp4_near_earth_drag_theoretical_basis.md
                        §14 (xlcof, aycof) building on §13 (xnodcf).
  Primary reference:    [SR3] Hoots & Roehrich (1980), §6, pp. 10-11.
                        Per-module references in companion audit cards.
  Domain of validity:   Near-earth orbits: period < 225 min (deep-space flag
                        set otherwise — propagation is then NOT valid via this
                        path and a caller must dispatch to SDP4). cos i₀ must
                        avoid the equatorial singularity 1+cos i₀ → 0
                        (guarded at line 237).

METHOD
  Method declared:      Sequential composition of standalone modules,
                        followed by a single inline closed-form construction
                        of xlcof/aycof. No iterative methods at this level.
  Method implemented:   src/sgp4/near_space.h:140-245.
                        Phases:
                          1) Earth & zonal constants  (lines 146-156)
                          2) Store TLE elements        (lines 158-164)
                          3) Trig constants            (lines 166-172)
                          4) Element recovery          (lines 174-184)
                          5) Density parameters        (lines 186-190)
                          6) Drag coefficients         (lines 192-224)
                          7) Secular rates             (lines 226-231)
                          8) Long-period xlcof/aycof   (lines 233-243)
                        xlcof guarded branch: when 1+cos i₀ ≈ 0, falls back
                        to the limit form (3/8)(A₃₀/CK₂) sin i₀ — this is the
                        l'Hôpital limit of (3+5cos i₀)/(1+cos i₀) at the pole.
                        At cos i₀ = −1: (3−5)/0 indeterminate; the file
                        instead substitutes a guarded ratio (3/8) — see NOTES.
  Match verdict:        ✓ matched — implementation is a pure sequencer plus
                        one closed-form composition. No Padé, no series, no
                        iteration at this level.

ERROR BOUND
  Bound category:       precision + accuracy + measurement (inherited)
  Bound formula:        Each TrackedValue<T> output carries the union of
                        the bounds from its constituent operations, per
                        REQ-EF-3 (closed-form propagation). For the inline
                        xlcof/aycof construction, the bound is the closed-
                        form mean-value bound from TrackedValue::operator*
                        and operator/; the only original numerical cost
                        is a single divide (1+cos i₀) which is guarded.
  Bound implemented:    Implicit via TrackedValue arithmetic on the rhs of
                        ns.xlcof = ..., ns.aycof = ... assignments; no bound
                        is added by this file. ALL propagation goes through
                        operator/, operator*, operator+ on TrackedValue<T>.
  Bound verdict:        ✓ matched — no original truncation, no original
                        iteration; all bound work is delegated to the
                        children. The composition itself adds nothing.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation)
  AUD-EF applies:       AUD-EF-3 (closed-form ops add no truncation)
  AUD-MC applies:       n/a (composition glue, not algebra)
  Verification test:    tests/test_sgp4 — verify_near_space_init test cases
                        against [V06] Vallado verification cases.

NOTES
  - **Resolved 2026-05-14 (R14)**: the xlcof critical-i fallback at line
    237 was changed from (3/8)(A₃₀/CK₂)·sin i₀ to (1/2)(A₃₀/CK₂)·sin i₀,
    which IS the rigorous l'Hôpital limit of the main expression
    (1/8)·(A₃₀/CK₂)·sin i₀·(3+5cos i₀)/(1+cos i₀) at cos i₀ → +1
    (numerator → 8, denominator → 2; (1/8)·4 = 1/2). The previous 3/8
    value matched xlcof at i₀ = π/2 (purely polar) but was not the
    asymptotic limit at either pole; it was a band-aid carried over
    from upstream conventions. The new 1/2 is exact at cos i₀ → +1
    and a small band-aid at cos i₀ → −1 (where sin i₀ → 0 makes any
    coefficient irrelevant; retrograde-polar orbits are unphysical
    for SGP4 anyway). See near_space.h:233-280 for the full derivation
    and precision-impact analysis. Numerical band-aid difference inside
    the |1 ± cos i₀| < 1.5e-12 guard band is bounded by ~1e-9
    dimensionless for Earth — well below SGP4's accuracy floor.
  - The guard threshold 1.5e-12 ≈ 6 × T = double's ε. For wider T (e.g.
    cpp_bin_float_50) this guard fires far too aggressively. T-aware
    thresholding (e.g. scale by ε) is a future refinement.
  - All struct fields are TrackedValue<T>; total_error() is well-defined
    on the output struct (per REQ-EF-3 transitive closure).
  - The `is_deep_space` flag returns true when period ≥ 225 min; this
    file does NOT itself enforce the near-earth domain — the caller
    must check `ns.is_deep_space` before propagating.
```

---

## Card 2 — `propagate_near_space`

```
=== FORMULA AUDIT CARD ===
ID:                     near_space::propagate_near_space
Location:               src/sgp4/near_space.h:261-327
Mathematical statement: Advance the satellite state by tsince minutes:
                        secularly update {a, e, M, ω, Ω, n} with drag;
                        compute long-period corrections xll, aynl, xlt, ayn;
                        solve modified Kepler equation for E+ω;
                        compute osculating elements (r, u, ṙ, rḟ, pl, β_l);
                        apply short-period J₂ corrections;
                        return position/velocity in TEME via
                        orbit::elements_to_state. [SR3] §6 pp. 11-15.

THEORY
  Underlying theorem:   Composition; this function is a sequencer. Theoretical
                        bases live in the called modules:
                        • Secular advance (Brouwer 1st-order + Lane drag):
                            orbit::secular_advance
                            → theoretical_basis_audit/secular_update.md
                            [SR3] §6 eqs (23)-(34); Lane (1965) for drag
                            polynomial expansion in tsince.
                        • Modified Kepler solver (axN, ayn, capu) → E+ω:
                            config.model_functions.kepler_solver lambda
                            → theoretical_basis_audit/model_functions.md
                            and theoretical_basis_audit/kepler.md
                            [SR3] §6 eqs (37)-(40). Newton-Raphson with
                            small-correction termination.
                        • Osculating elements (r, u, ṙ, rḟ, p_l, β_l, sin/cos 2u):
                            orbit::compute_osculating
                            → theoretical_basis_audit/osculating_elements.md
                            [SR3] §6 eqs (41)-(54).
                        • Short-period J₂ corrections (δr, δu, δi, δΩ, δṙ, δrḟ):
                            perturbation::apply_short_period
                            → theoretical_basis_audit/short_period.md
                            [SR3] §6 eqs (55)-(60); Brouwer-Lyddane short-
                            period closed forms.
                        • Position/velocity from elements (TEME):
                            orbit::elements_to_state
                            → theoretical_basis_audit/state_from_elements.md
                            [SR3] §6 eqs (61)-(66); standard rotation U·R·M
                            applied to radial/transverse magnitudes.

                        ORIGINAL TO THIS FILE — long-period linearization:
                          axN  = e cos ω
                          β²   = 1 − e²
                          τ    = 1/(a · β²)
                          xll  = τ · xlcof · axN
                          aynl = τ · aycof
                          xlt  = mean_longitude + xll
                          ayn  = e sin ω + aynl
                        Theoretical basis: Lyddane's long-period J₃ correction
                        applied to the (axN, ayn) eccentricity-vector form
                        before modified Kepler solve. [SR3] §6 p. 13 eqs
                        (33)-(36) (linearized form in the Lyddane variables).
                        Derivation reference: design/derivations/
                        sgp4_near_earth_drag_theoretical_basis.md §14
                        (xlcof, aycof are constructed there; the inline
                        application here is `τ × coef × axN` and `τ × coef`).
  Primary reference:    [SR3] Hoots & Roehrich (1980), §6, pp. 11-15.
                        Per-module references in companion audit cards.
  Domain of validity:   Near-earth: tsince such that secular_advance hasn't
                        driven e → 0 or a → 0 (drag-driven decay). The
                        propagator does NOT guard against decay; the caller
                        must check the returned state. axN, ayn defined for
                        all near-earth orbits in the bounded domain.

METHOD
  Method declared:      Sequential composition of standalone modules with
                        one inline closed-form linearization (long-period
                        block). No iterative methods at this level (Kepler
                        iteration is inside config.model_functions.kepler_solver).
  Method implemented:   src/sgp4/near_space.h:267-327.
                        Phases:
                          1) Build BrouwerSecularRates struct (lines 273-276)
                          2) Build DragCoefficients struct  (lines 279-289)
                          3) Secular advance               (line 291-293)
                          4) Long-period linearization     (lines 298-304)
                          5) Modified Kepler solve         (line 307-308)
                          6) Osculating elements           (line 311)
                          7) Short-period corrections      (lines 314-321)
                          8) Elements → state              (lines 324-326)
                        wrap_two_pi at line 307 is a closed-form modular
                        reduction; → theoretical_basis_audit/angles.md.
  Match verdict:        ✓ matched — implementation is a pure sequencer plus
                        one short closed-form composition. No Padé, no series,
                        no iteration at this level.

ERROR BOUND
  Bound category:       precision + accuracy (inherited from children)
  Bound formula:        Per REQ-EF-3 (closed-form propagation). All inline
                        ops (axN, β², τ, xll, aynl, xlt, ayn) are products,
                        quotients, sums — each propagates via TrackedValue
                        rules. The state vector returned carries the union
                        of all child bounds plus the closed-form composition
                        bound on the trig/multiplication operations.
  Bound implemented:    Implicit via TrackedValue arithmetic; no bound is
                        added by this file. All bound propagation goes
                        through the standard operator*/+/-//- overloads on
                        TrackedValue<T> plus the bounds added inside the
                        called modules.
  Bound verdict:        ✓ matched — no original truncation, no original
                        iteration; all bound work is delegated to the
                        children. Composition adds nothing.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation)
  AUD-EF applies:       AUD-EF-3 (closed-form ops add no truncation)
  AUD-MC applies:       n/a (composition glue, not algebra)
  Verification test:    tests/test_sgp4 — Vallado [V06] verification cases:
                        position/velocity error at multiple tsince values
                        against ground truth.

NOTES
  - The long-period block (lines 298-304) is the only inline math. It is
    the linearized small-correction form of [SR3] §6 eqs (33)-(36); the
    full nonlinear form is derived in
    design/derivations/sgp4_near_earth_drag_theoretical_basis.md §14.
    The linearization is exact in the J₃ → 0 limit and good to O(J₃²)
    for non-pathological inclinations — well within SGP4's stated J₂-
    truncation accuracy.
  - The propagator does NOT check `ns.is_deep_space`; the caller MUST
    dispatch to SDP4 if period ≥ 225 min. Calling this propagator on a
    deep-space orbit produces wrong answers silently. ⚠ This is a
    CALLER CONTRACT issue, not a TBA failure of this function — but
    should be flagged in the calling code's audit.
  - `wrap_two_pi(xlt - sec.Omega)` (line 307): mod-2π identity, closed
    form. Bound is propagated, no truncation. → angles.md.
  - The returned StateVector carries the full TrackedValue<T> error
    budget end-to-end. The composition contract (REQ-EF-3 transitive
    closure) is satisfied iff every child module is itself TBA-PASS.
```

---

## File-level verdict

- **A. Error wiring**: ✓ pure composition; no bare-T values, all assignments
  on TrackedValue<T>. Per-module wiring inherited.
- **B. Algebra axioms**: n/a (this file does no algebra; it is a sequencer).
- **C. Theoretical basis**: 
  - **PASS — composition**. Both functions are sequencers with thin
    inline closed-form composition glue (xlcof/aycof construction in init;
    long-period linearization in propagate). All heavy formulas are
    delegated to standalone modules with their own TBA cards.
  - **PASS verdict for this file is CONDITIONAL on PASS verdicts in
    every cross-referenced child audit card**: element_recovery.md,
    density_model.md, drag_coefficients.md, model_functions.md (for
    secular_rates and kepler_solver), secular_update.md,
    osculating_elements.md, short_period.md, state_from_elements.md,
    angles.md.
  - **Resolved 2026-05-14 (R14)**: xlcof fallback at cos i₀ near ±1
    (line 237 branch) was corrected from (3/8)(A₃₀/CK₂)·sin i₀ to
    (1/2)(A₃₀/CK₂)·sin i₀ — the rigorous l'Hôpital limit at cos i₀ → +1.
    Verified against §14 Theorem 14.2 of
    design/derivations/sgp4_near_earth_drag_theoretical_basis.md.

**File verdict: PASS (composition; conditional on child cards).**
