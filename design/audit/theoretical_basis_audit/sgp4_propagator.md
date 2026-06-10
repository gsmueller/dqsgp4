# Theoretical Basis Audit — `src/sgp4/sgp4_propagator.h`

**File**: `src/sgp4/sgp4_propagator.h` (169 lines)
**Scope**: Top-level SGP4 propagator. Pure dispatch/assembly layer — contains zero hardcoded formulas. All numeric work is delegated to `near_space.h`, `deep_space.h`, and lambdas in `ModelConfiguration<T>`.
**Audit unit**: 4 cards (ctor, propagate, model_description, is_deep_space).

---

## Card 1 — `Propagator::Propagator(config, elements, tolerance)` (constructor)

```
=== FORMULA AUDIT CARD ===
ID:                     sgp4_propagator::Propagator::ctor
Location:               src/sgp4/sgp4_propagator.h:104-114, initialize() at 147-155
Mathematical statement: No closed-form formula. Dispatches:
                          ns_init_ ← initialize_near_space(config, elements, tol)
                          use_deep_space_ ← ns_init_.is_deep_space
                          if use_deep_space_:
                              ds_init_ ← initialize_deep_space(config, elements, tol)
                        Orbit-type classification is performed inside
                        initialize_near_space (criterion: 2π / n₀ ≥ 225 min,
                        equivalently a₀ ≳ 1.524 ER per SR3 §6).

THEORY
  Underlying theorem:   Hoots & Roehrich (1980) Spacetrack Report #3 §6,
                        "Selection of perturbation theory": the near-space
                        SGP4 model applies for orbital periods < 225 minutes;
                        the deep-space SDP4 model applies otherwise. This is
                        a definitional/procedural rule, not a theorem in the
                        analytic sense.
  Primary reference:    Hoots & Roehrich (1980) SR#3 §3 (initialization) and
                        §6 (model selection criterion). Vallado et al. (2006)
                        "Revisiting Spacetrack Report #3" Rev 3 §3.
  Domain of validity:   Any TLE with valid mean motion n₀ > 0 and e₀ ∈ [0, 1).
                        Threshold 225 min is a hard procedural boundary; no
                        analytic continuity is claimed across it.

METHOD
  Method declared:      Procedural dispatch — no series, no iteration, no
                        approximation at this layer. Calls free functions
                        initialize_near_space (always) and
                        initialize_deep_space (conditionally).
  Method implemented:   src/sgp4/sgp4_propagator.h:147-155. Two function
                        calls into other TUs; bool flag assignment.
  Match verdict:        ✓ matched — pure dispatch matches the procedural
                        nature of the cited rule.

ERROR BOUND
  Bound category:       n/a at this layer. All TrackedValue<T> error bookkeeping
                        is owned by initialize_near_space / initialize_deep_space
                        (cards in near_space.md, deep_space.md). The ctor adds
                        zero new error.
  Bound formula:        n/a — composition only. Outputs ns_init_ / ds_init_
                        carry their own errors per REQ-EF-3.
  Bound implemented:    n/a — no arithmetic.
  Bound verdict:        ✓ correct — dispatch layer must not invent or absorb
                        error; it composes other layers' bounds.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form composition: TrackedValue passes
                        through ctor unchanged).
  AUD-EF applies:       n/a directly; AUD-EF wiring is checked in the
                        initialize_near_space / initialize_deep_space cards.
  AUD-MC applies:       n/a (no algebra operation here).
  Verification test:    tests/test_sgp4/ — any test that constructs a
                        Propagator and inspects is_deep_space() against a
                        TLE with a known period exercises this dispatch.

NOTES
  - The "near-space init runs first to determine orbit type" pattern is
    a SR#3-faithful design: SR#3 §3 computes a₀, n₀, and other near-space
    constants regardless, then decides whether to also build the deep-space
    long-period/resonance machinery.
  - When use_deep_space_ is true, ns_init_ is still populated and reused
    by deep-space code paths (lunar-solar terms layer atop near-space
    secular rates). This is correct per SR#3 §3.
  - No threshold appears explicitly in this file; it is encoded inside
    initialize_near_space and sets ns_init_.is_deep_space.
```

---

## Card 2 — `Propagator::propagate(tsince_minutes)`

```
=== FORMULA AUDIT CARD ===
ID:                     sgp4_propagator::Propagator::propagate
Location:               src/sgp4/sgp4_propagator.h:121-127, with delegating
                        helpers propagate_near_space (158-160) and
                        propagate_deep_space (163-165).
Mathematical statement: r(tsince), v(tsince) in TEME frame, returned as
                        StateVector<T>. Dispatches:
                          if use_deep_space_: → sgp4::propagate_deep_space(...)
                          else:               → sgp4::propagate_near_space(...)
                        Both downstream calls implement SR#3 §§3-5
                        (near-space SGP4) or SR#3 §§3-5 + §§4 lunar-solar +
                        §5 resonance (deep-space SDP4).

THEORY
  Underlying theorem:   Hoots & Roehrich (1980) SR#3 §1, "The propagator
                        evaluates the analytic solution of the mean
                        equations at time t = epoch + tsince and recovers
                        osculating elements via short-period corrections."
                        The branch follows SR#3 §6.
  Primary reference:    Hoots & Roehrich (1980) SR#3 §§1, 6. Vallado et al.
                        (2006) Rev 3 §3.
  Domain of validity:   tsince ∈ ℝ (any sign; SR#3 supports backward
                        propagation). Validity for very large |tsince|
                        degrades per Lane drag model assumptions; that is
                        a downstream concern (drag.md).

METHOD
  Method declared:      Procedural dispatch on use_deep_space_; no formula
                        at this layer.
  Method implemented:   src/sgp4/sgp4_propagator.h:121-127. One bool check,
                        one of two free-function calls. Returns whatever
                        the downstream propagator returns (StateVector<T>
                        already carries TrackedValue<T> errors).
  Match verdict:        ✓ matched — dispatch matches the cited "branch on
                        orbit type" procedural rule.

ERROR BOUND
  Bound category:       n/a at this layer. The downstream propagate_*
                        functions populate StateVector<T>'s
                        precision/accuracy/measurement budget per REQ-EF.
  Bound formula:        n/a — pure pass-through. Composition under REQ-EF-3:
                        the returned StateVector<T> carries the exact bound
                        produced by the downstream propagator with no
                        addition or subtraction at this boundary.
  Bound implemented:    n/a — `return propagate_xxx_space(...)` directly.
  Bound verdict:        ✓ correct — no error injected, no error dropped.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (pass-through composition).
  AUD-EF applies:       AUD-EF wiring of the downstream propagators is the
                        meaningful test (see near_space.md, deep_space.md).
  AUD-MC applies:       n/a.
  Verification test:    tests/test_sgp4/ — any propagation regression test
                        exercises this dispatch. The Vallado test-vector
                        suite (tests/test_sgp4_vallado/) is the primary
                        end-to-end check.

NOTES
  - `const` qualifier on propagate() is correct: PropagationState (mutable
    integrator state for deep-space resonance) is owned downstream inside
    propagate_deep_space, not in this class. Lines 78-80 declare
    PropagationState but the class does not store it as a member —
    integrator state lives inside the deep-space free function. ✓
  - No tolerance is passed to propagate() at runtime; the init-time
    tolerance_ member is forwarded. Correct: tolerance is a model
    property, not a per-call parameter, matching REQ-EF philosophy.
  - The `if (use_deep_space_)` branch is a runtime data-driven dispatch;
    compile-time selection is not possible because the orbit type comes
    from TLE data parsed at construction.
```

---

## Card 3 — `Propagator::model_description() const`

```
=== FORMULA AUDIT CARD ===
ID:                     sgp4_propagator::Propagator::model_description
Location:               src/sgp4/sgp4_propagator.h:130
Mathematical statement: Returns const std::string& — the description string
                        held by the ModelConfiguration (e.g. "WGS72 + Lane
                        drag" or similar). Not a numeric computation.

THEORY
  Underlying theorem:   n/a — metadata accessor.
  Primary reference:    n/a.
  Domain of validity:   Always valid; returns by const reference to
                        config_.description, lifetime tied to the
                        Propagator instance.

METHOD
  Method declared:      Constant-time field accessor.
  Method implemented:   `return config_.description;` (line 130).
  Match verdict:        ✓ matched — trivial accessor.

ERROR BOUND
  Bound category:       n/a — non-numeric.
  Bound formula:        n/a.
  Bound implemented:    n/a.
  Bound verdict:        n/a.

CROSS-AUDIT
  REQ-EF applies:       n/a (no TrackedValue<T> involved).
  AUD-EF applies:       n/a.
  AUD-MC applies:       n/a.
  Verification test:    Any test that prints or asserts the model
                        description string; not safety-critical.

NOTES
  - Pure metadata; out of scope for TBA in any substantive sense. Card
    included only for completeness of the 4-function inventory.
```

---

## Card 4 — `Propagator::is_deep_space() const`

```
=== FORMULA AUDIT CARD ===
ID:                     sgp4_propagator::Propagator::is_deep_space
Location:               src/sgp4/sgp4_propagator.h:133
Mathematical statement: Returns the bool set at construction by the
                        orbit-type classifier inside initialize_near_space
                        (criterion: orbital period ≥ 225 min per SR#3 §6).

THEORY
  Underlying theorem:   Hoots & Roehrich (1980) SR#3 §6 — definitional
                        "deep space" threshold. Sourced from ns_init_, not
                        recomputed here.
  Primary reference:    SR#3 §6.
  Domain of validity:   Set once at construction; immutable thereafter for
                        the lifetime of the Propagator. Returns the
                        construction-time decision.

METHOD
  Method declared:      Constant-time field accessor.
  Method implemented:   `return use_deep_space_;` (line 133).
  Match verdict:        ✓ matched — trivial accessor.

ERROR BOUND
  Bound category:       n/a — bool, not TrackedValue<T>.
  Bound formula:        n/a.
  Bound implemented:    n/a.
  Bound verdict:        n/a.

CROSS-AUDIT
  REQ-EF applies:       n/a.
  AUD-EF applies:       n/a.
  AUD-MC applies:       n/a.
  Verification test:    tests/test_sgp4/ — typically asserted against
                        TLEs with known orbit periods (e.g. GEO TLE →
                        true; LEO TLE → false). The actual classification
                        logic is tested at the near_space.md layer.

NOTES
  - Returning a bool, not a TrackedValue<bool>, is appropriate: the
    threshold 225 min is a hard procedural switch, not a measurement.
    No error budget applies.
  - Whether the threshold itself is correctly encoded inside
    initialize_near_space is a question for near_space.md, not this
    file's card.
```

---

## File-level verdict for `sgp4_propagator.h`

This is a **pure dispatch / assembly layer**. The file header explicitly states "Contains ZERO hardcoded formulas. Every computation is provided via ModelFunctions lambdas from the ModelConfiguration." The audit confirms this claim:

- **A. Error wiring**: ✓ no error is invented, absorbed, or dropped at this layer. Composition under REQ-EF-3 (TrackedValue<T> pass-through) is correct.
- **B. Algebra axioms**: n/a — no algebra operations.
- **C. Theoretical basis**:
  - **Card 1 (ctor)**: ✓ matched. Procedural dispatch per SR#3 §6. **PASS.**
  - **Card 2 (propagate)**: ✓ matched. Pass-through dispatch per SR#3 §1/§6. **PASS.**
  - **Card 3 (model_description)**: n/a — metadata accessor. **PASS (trivial).**
  - **Card 4 (is_deep_space)**: n/a — bool accessor of construction-time classification. **PASS (trivial).**

**File verdict: PASS** — the file is a thin dispatcher; substantive TBA work for SGP4 lives in `near_space.md`, `deep_space.md`, `model_selector.md`, and the drag/zonal/lambda cards.

**Out-of-scope-here / handled elsewhere:**
- The 225-min threshold encoding (classification logic) — handled in `near_space.md` (or `model_selector.md`).
- Population of `NearSpaceCoefficients<T>` and `DeepSpaceCoefficients<T>` fields (c1, c4, …, sgh2, …, d2201, …) — handled in the corresponding init cards. These structs are *declared* here but *populated* by free functions in `near_space.h` / `deep_space.h`.
- `PropagationState<T>` (xli, xni, atime) is *declared* here at lines 78-80 but the Propagator class does not store it as a member. Mutable resonance integrator state lives inside `propagate_deep_space` (file: `deep_space.h`); audit the persistence-between-calls claim there, not here.

**Open questions / flags:** None for this file. The dispatch pattern is clean and SR#3-faithful.
