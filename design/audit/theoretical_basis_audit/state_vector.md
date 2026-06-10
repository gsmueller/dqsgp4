# Theoretical Basis Audit — `src/sgp4/state_vector.h`

File: `src/sgp4/state_vector.h` (90 lines).
Five numeric-producing member functions on `StateVector<T>`, each combining per-component error magnitudes across the three Cartesian components (x, y, z) into a single scalar magnitude via the **RSS (root-sum-square)** operation `sqrt(a² + b² + c²)`.

Reference for the RSS-vs-bound distinction:
- `src/math/tracked_value.h:31-40` — `ThreeErrors::total()` is the triangle-inequality (rigorous) combiner; `ThreeErrors::rss()` is documented as the "Statistical estimate (RSS). Use when reporting expected error rather than worst-case bound. Only valid if sources are independent."
- `design/specifications/error_framework.md` REQ-EF-1 — three categories are independent and must not be merged.
- `design/specifications/error_framework.md` REQ-EF-2 — `total_error()` is the rigorous upper bound `|x.value − x_true| ≤ x.total_error()` via the triangle inequality.

Below, the **spatial direction** (across x/y/z components) is consistently RSS in this file. The audit asks: for each function, is RSS the correct combination for what is being reported?

---

## Card 1 — `StateVector<T>::position_error()`

```
=== FORMULA AUDIT CARD ===
ID:                     state_vector::position_error
Location:               src/sgp4/state_vector.h:37-45
Mathematical statement: ε_r = sqrt( e_x² + e_y² + e_z² )
                        where e_i = position_km.i.total_error()
                              (triangle-sum of the i-component's three categories)

THEORY
  Underlying theorem:   Reverse triangle inequality on the Euclidean norm.
                        For r̂ = r + Δr with |Δr_i| ≤ e_i per component,
                          ‖Δr‖₂ ≤ sqrt(e_x² + e_y² + e_z²),
                        and hence | ‖r̂‖ − ‖r‖ | ≤ ‖Δr‖₂ ≤ sqrt(e_x²+e_y²+e_z²).
                        This is a rigorous Euclidean-magnitude bound, distinct
                        from RSS as a "statistical estimate" because here the
                        Pythagorean inequality is *exact* for the worst-case
                        displacement (each axis is orthogonal and independent
                        by construction of the basis).
  Primary reference:    Standard inner-product-space inequality; e.g.,
                        Horn & Johnson (1985) §5.2, Cauchy–Schwarz corollary.
  Domain of validity:   Any orthonormal Cartesian basis; the TEME frame here.

METHOD
  Method declared:      RSS across three orthogonal components. Each input
                        e_i is itself the rigorous triangle-inequality
                        total_error() of the i-component (per REQ-EF-2).
  Method implemented:   `sqrt( px*px + py*py + pz*pz ).value` with
                        px,py,pz = position_km.{x,y,z}.total_error().
                        Wraps the sqrt of a TrackedValue (zero-error
                        constructor) to reuse the sqrt operator, then
                        extracts .value and packages into a TrackedValue
                        with errors.{measurement,precision,accuracy} all
                        set to T(0).
  Match verdict:        ⚠ method-vs-return-type mismatch. The spatial RSS
                        across orthogonal components is correct for a bound
                        on Euclidean-vector-magnitude error. But the result
                        is returned as `TrackedValue<T>` with zero in all
                        three error fields, which is misleading — the value
                        IS itself an error magnitude (already a bound), not
                        a tracked value with a separable error budget.
                        Cf. position_measurement_error() etc. (Cards 3-5)
                        which correctly return plain T.

ERROR BOUND
  Bound category:       n/a — the returned quantity IS itself a bound.
  Bound formula:        N/A. By REQ-EF-2 the input e_i already satisfies
                        |Δr_i| ≤ e_i; the sqrt(Σ e_i²) is then a rigorous
                        upper bound on ‖Δr‖₂.
  Bound implemented:    The TrackedValue<T> returned has errors = (0,0,0).
                        Strictly, the sqrt() call's own intrinsic precision
                        bound (from REQ-EF-3 sqrt formula) is *discarded*:
                        the code unwraps `.value` and rebuilds with zero
                        errors, dropping the few ulps of sqrt roundoff.
  Bound verdict:        ⚠ approximate — for the *Euclidean magnitude bound*
                        purpose, the dropped sqrt-roundoff bound is many
                        orders of magnitude below the reported magnitude
                        and the discard is harmless in practice. But the
                        framework contract (REQ-EF-2) says every operation
                        propagates its errors; deliberately zeroing is a
                        framework-wiring deviation. The function should
                        either (a) return plain T (as Cards 3-5 do) or
                        (b) preserve the sqrt's own propagated bound.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (three-category separation — but this
                                  function flattens into a single magnitude,
                                  so REQ-EF-1's separation is intentionally
                                  abandoned at the reporting boundary).
                        REQ-EF-2 (input e_i are rigorous upper bounds; the
                                  RSS combination preserves rigor for
                                  Euclidean-vector-magnitude error).
  AUD-EF applies:       AUD-EF-? (final-report aggregation — not in §3 of
                                  the framework table; this is a reporting
                                  surface, not a propagation step).
  AUD-MC applies:       n/a (no algebra-axiom test directly on this
                              aggregator).
  Verification test:    None located in this read-only audit. Expected
                        location: tests/test_sgp4/ would exercise
                        position_error() against a known propagated case.

NOTES
  - The doc comment "Combines measurement, precision, and accuracy into a
    single bound" is accurate provided one reads total_error() (called per
    component) as REQ-EF-2's rigorous triangle bound. The per-category
    information is intentionally lost at this boundary; that's the point
    of position_error() vs. the three category-specific functions.
  - Returning TrackedValue<T> with zero errors is type-aesthetic only;
    callers cannot decompose the result further. Recommend returning T
    for symmetry with Cards 3-5.
  - Caller-side note ("RSS is a statistical estimate, not a rigorous
    upper bound"): this caveat is *generally* true for RSS of arbitrary
    correlated sources. For the *spatial* RSS across orthogonal Cartesian
    components, RSS is a rigorous Euclidean-norm bound. The framework
    note applies to the *category* axis (where RSS would be unsound and
    triangle-sum is mandated by REQ-EF-2) — and this file correctly uses
    RSS only across spatial components and triangle-sum (via
    total_error()) across categories. **This is the right separation.**
  - The `sqrt(...).value` pattern abuses the TrackedValue<T> sqrt operator
    just to access the underlying T sqrt. A direct `using std::sqrt; sqrt(...)`
    on a plain T would be cleaner (as Cards 3-5 do).
```

---

## Card 2 — `StateVector<T>::velocity_error()`

```
=== FORMULA AUDIT CARD ===
ID:                     state_vector::velocity_error
Location:               src/sgp4/state_vector.h:48-56
Mathematical statement: ε_v = sqrt( e_vx² + e_vy² + e_vz² )
                        where e_vi = velocity_km_s.i.total_error()

THEORY
  Underlying theorem:   Same as Card 1 — reverse triangle inequality on
                        the Euclidean norm, applied to the velocity vector
                        in the TEME frame.
  Primary reference:    Same as Card 1.
  Domain of validity:   TEME orthonormal Cartesian basis.

METHOD
  Method declared:      RSS across three orthogonal velocity components;
                        each input is the rigorous total_error() per
                        component.
  Method implemented:   Same pattern as Card 1: sqrt of TrackedValue with
                        squared per-component total_errors, then unwrap
                        .value and repackage with zero error fields.
  Match verdict:        ⚠ same return-type observation as Card 1.

ERROR BOUND
  Bound category:       n/a — the returned quantity IS itself a bound on
                        Euclidean velocity-magnitude error [km/s].
  Bound formula:        |Δv| ≤ sqrt(e_vx² + e_vy² + e_vz²) by the same
                        Pythagorean argument as Card 1.
  Bound implemented:    Returns TrackedValue<T> with errors zeroed; sqrt
                        roundoff bound discarded.
  Bound verdict:        ⚠ approximate — same comment as Card 1.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-2 (same as Card 1).
  AUD-EF applies:       Same as Card 1.
  AUD-MC applies:       n/a.
  Verification test:    tests/test_sgp4/ (not located in this audit).

NOTES
  - Identical structure to position_error(); same observations apply.
  - The doc comment lacks the "combines measurement/precision/accuracy"
    line that position_error() has; minor inconsistency, not a fault.
```

---

## Card 3 — `StateVector<T>::position_measurement_error()`

```
=== FORMULA AUDIT CARD ===
ID:                     state_vector::position_measurement_error
Location:               src/sgp4/state_vector.h:60-66
Mathematical statement: ε_m = sqrt( m_x² + m_y² + m_z² )
                        where m_i = position_km.i.errors.measurement
                        (the measurement-category error of the i-component)

THEORY
  Underlying theorem:   Reverse triangle inequality applied to the
                        *measurement-only* contribution to position error.
                        Treats the measurement category as the magnitude
                        of a "measurement-error vector" in the same TEME
                        basis; its Euclidean norm is bounded by the RSS
                        of the per-axis measurement-category magnitudes.
  Primary reference:    Same Euclidean-norm inequality as Cards 1-2.
                        Conceptual basis: REQ-EF-1 keeps the measurement
                        category separable so its contribution can be
                        reported in isolation.
  Domain of validity:   TEME orthonormal basis; measurement category as
                        defined by REQ-EF-1.

METHOD
  Method declared:      RSS across three orthogonal components of the
                        measurement-category sub-error.
  Method implemented:   `using std::sqrt; sqrt(mx*mx + my*my + mz*mz)`
                        with mx = position_km.x.errors.measurement, etc.
                        Direct on plain T — no TrackedValue wrapping.
  Match verdict:        ✓ matched. The spatial RSS is the rigorous
                        Euclidean-norm bound for the measurement-only
                        contribution; the per-category isolation (no
                        mixing with precision or accuracy) is preserved
                        per REQ-EF-1.

ERROR BOUND
  Bound category:       The returned T IS the measurement-category bound
                        on Euclidean position-error magnitude.
  Bound formula:        Returns a scalar bound; no further error tracking.
                        sqrt's own roundoff is order ε_T · |result| and
                        is not propagated (this is a reporting function,
                        a leaf of the computation graph).
  Bound implemented:    Plain T return; no error budget attached.
  Bound verdict:        ✓ rigorous within type precision for the stated
                        scope (Euclidean-norm bound on measurement
                        contribution). The sqrt roundoff is far below
                        any practically relevant magnitude for a
                        reporting function.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (measurement category isolated).
  AUD-EF applies:       The framework's AUD-EF tests target propagation
                        steps; this is a reporting boundary.
  AUD-MC applies:       n/a.
  Verification test:    Not located.

NOTES
  - Return type T (not TrackedValue<T>) is consistent with the
    "reporting boundary" semantics: this is a final number, not a
    further-propagable quantity.
  - The RSS-vs-bound caveat: spatially, RSS *is* rigorous for Euclidean
    magnitude. The framework's "RSS = statistical estimate" warning
    applies to RSS *across error categories*, which this function does
    not do (it stays inside the measurement category).
```

---

## Card 4 — `StateVector<T>::position_precision_error()`

```
=== FORMULA AUDIT CARD ===
ID:                     state_vector::position_precision_error
Location:               src/sgp4/state_vector.h:70-76
Mathematical statement: ε_p = sqrt( p_x² + p_y² + p_z² )
                        where p_i = position_km.i.errors.precision

THEORY
  Underlying theorem:   Reverse triangle inequality on the precision-
                        category contribution. Same Euclidean-norm
                        argument as Card 3, restricted to the precision
                        category (floating-point roundoff and truncation
                        accumulated through the SGP4 chain).
  Primary reference:    Same as Card 3; REQ-EF-1 for category isolation.
  Domain of validity:   TEME orthonormal basis; precision category per
                        REQ-EF-1.

METHOD
  Method declared:      RSS across three orthogonal precision-category
                        sub-errors.
  Method implemented:   `using std::sqrt; sqrt(px*px + py*py + pz*pz)`
                        with px = position_km.x.errors.precision, etc.
  Match verdict:        ✓ matched, same as Card 3.

ERROR BOUND
  Bound category:       Returns the precision-category Euclidean-norm
                        bound on position-error magnitude.
  Bound formula:        Same as Card 3.
  Bound implemented:    Plain T.
  Bound verdict:        ✓ rigorous within type precision.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       Reporting boundary.
  AUD-MC applies:       n/a.
  Verification test:    Not located.

NOTES
  - Doc comment "From floating-point arithmetic accumulation through
    ~200 operations" describes the source category content but is not
    a method or bound claim — informational.
  - Same RSS-vs-rigor commentary as Card 3 applies.
```

---

## Card 5 — `StateVector<T>::position_accuracy_error()`

```
=== FORMULA AUDIT CARD ===
ID:                     state_vector::position_accuracy_error
Location:               src/sgp4/state_vector.h:80-86
Mathematical statement: ε_a = sqrt( a_x² + a_y² + a_z² )
                        where a_i = position_km.i.errors.accuracy

THEORY
  Underlying theorem:   Reverse triangle inequality on the accuracy-
                        category contribution. Same Euclidean-norm
                        argument as Cards 3-4, restricted to the accuracy
                        category (model truncation: Brouwer series order,
                        simplified drag, etc.).
  Primary reference:    Same as Card 3; REQ-EF-1.
  Domain of validity:   TEME orthonormal basis; accuracy category per
                        REQ-EF-1.

METHOD
  Method declared:      RSS across three orthogonal accuracy-category
                        sub-errors.
  Method implemented:   `using std::sqrt; sqrt(ax*ax + ay*ay + az*az)`
                        with ax = position_km.x.errors.accuracy, etc.
  Match verdict:        ✓ matched, same as Card 3.

ERROR BOUND
  Bound category:       Returns the accuracy-category Euclidean-norm
                        bound on position-error magnitude.
  Bound formula:        Same as Card 3.
  Bound implemented:    Plain T.
  Bound verdict:        ✓ rigorous within type precision.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       Reporting boundary.
  AUD-MC applies:       n/a.
  Verification test:    Not located.

NOTES
  - Doc comment "From model truncation: Brouwer series order, simplified
    drag, etc." is informational.
  - Same RSS-vs-rigor commentary as Card 3 applies.
```

---

## File-level verdict for `state_vector.h`

- **A. Error wiring (REQ-EF / AUD-EF)**:
  - Cards 1-2 return `TrackedValue<T>` with all three error fields zeroed; the wrapping is cosmetic and discards the sqrt operator's own propagated bound. ⚠ minor framework-wiring deviation — this is a reporting boundary so the dropped bound has no downstream consumer, but the pattern is inconsistent with Cards 3-5.
  - Cards 3-5 return plain `T`, which is consistent with the leaf-reporting semantics.
- **B. Algebra axioms (AUD-MC)**: n/a — these are reporting aggregators, not algebra operations.
- **C. Theoretical basis (TBA)**:
  - All five cards: method = spatial RSS across orthogonal Cartesian components. Cited theory = reverse triangle inequality on Euclidean norm, giving the rigorous bound `|‖r̂‖ − ‖r‖| ≤ ‖Δr‖₂ ≤ sqrt(Σ Δr_i²)`. The "RSS = statistical estimate" caveat (`tracked_value.h:35-36`) applies to RSS *across error categories* (where category sources are not necessarily independent and triangle-sum is mandated by REQ-EF-2); it does NOT apply to spatial RSS across orthogonal Cartesian axes. **The file correctly separates the two axes**: across categories, it uses `total_error()` (triangle sum, rigorous); across spatial components, it uses RSS (rigorous Euclidean-norm bound). ✓

**Per-function summary**:

| Card | Function                          | Method ✓/?/⚠/✗ | Bound ✓/?/⚠/✗ | Verdict     |
|------|-----------------------------------|-----------------|-----------------|-------------|
| 1    | `position_error()`                | ⚠ return-type  | ⚠ sqrt-roundoff dropped | PASS w/ note |
| 2    | `velocity_error()`                | ⚠ return-type  | ⚠ sqrt-roundoff dropped | PASS w/ note |
| 3    | `position_measurement_error()`    | ✓               | ✓               | PASS        |
| 4    | `position_precision_error()`      | ✓               | ✓               | PASS        |
| 5    | `position_accuracy_error()`       | ✓               | ✓               | PASS        |

**File verdict: PASS with notes.** The RSS-vs-rigorous-bound concern raised in the audit task is correctly handled by the code: spatial RSS across orthogonal axes is a rigorous Euclidean-norm bound (not a statistical estimate), and the file does NOT use RSS across categories — that axis goes through `total_error()` (triangle sum) per REQ-EF-2. Two minor cosmetic issues on Cards 1-2: the `TrackedValue<T>`-with-zero-errors return type is inconsistent with the plain-T return of Cards 3-5, and the sqrt operator's own propagated bound is discarded by the `.value` unwrap.

**Recommended remediations** (non-blocking):
1. Change Cards 1-2 to return plain `T` for symmetry with Cards 3-5, or preserve the sqrt-operator's own propagated bound.
2. Add a test in `tests/test_sgp4/` that constructs a `StateVector<T>` with known per-component errors and verifies all five aggregators against hand-computed RSS values.
3. Add a comment near the function bodies clarifying that spatial RSS is rigorous (Pythagorean) and distinguishing it from the across-category RSS warned about in `tracked_value.h:35-36`.
