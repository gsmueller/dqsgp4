# Theoretical Basis Audit — `src/sgp4/model_functions.h`

**File**: `src/sgp4/model_functions.h` (287 lines)
**Audit framework**: `design/audit/theoretical_basis_audit.md` §1, §5
**Reference**: Hoots & Roehrich (1980) Spacetrack Report #3 (SR3).

## Scope

This file defines the `ModelFunctions<T>` lambda-injection bundle (5 function-type
typedefs + struct + a single factory `standard_sgp4()`). It is **predominantly a
plumbing wrapper**: the factory binds existing module implementations
(`perturbation::compute_secular_rates`, `perturbation::inclination_function`,
`astronomy::compute_gmst`) into `std::function` slots. Those implementations are
audited in their own cards (`brouwer.md`, `kaula.md`, `sidereal_time.md` — outside
this file).

Three slots, however, **define their own numerics in-line** inside the factory
and are in scope for this audit:

| Slot | In-line? | Audit card |
|---|---|---|
| `secular_rates` | no (wraps `perturbation::compute_secular_rates`) | n/a here |
| `inclination_function` | no (wraps `perturbation::inclination_function`) | n/a here |
| `kepler_solver` | **yes** (lambda body lines 234-255) | §1 below |
| `sidereal_time` | no (wraps `astronomy::compute_gmst`) | n/a here |
| `drag_coefficients` | **yes** but throws (no numerics) | §2 below |

Function count audited here: **1 formula card** (`kepler_solver`), plus one
non-numeric notice (`drag_coefficients` stub).

---

## 1 `standard_sgp4()::kepler_solver` lambda

```
=== FORMULA AUDIT CARD ===
ID:                     model_functions::standard_sgp4::kepler_solver
Location:               src/sgp4/model_functions.h:234-255
Mathematical statement: Solve the SR3 modified Kepler equation
                          x + ayn·cos(x) - axN·sin(x) = U
                        for x = E + ω, given axN = e·cos(ω),
                        ayn = e·sin(ω) + long-period correction, U = IL_T − Ω.

THEORY
  Underlying theorem:   Newton's method applied to
                          f(x) := U − ayn·cos(x) + axN·sin(x) − x
                          f'(x) = -ayn·sin(x) - axN·cos(x) + 1.
                        Newton-Kantorovich: convergence is quadratic in a
                        neighborhood of the root when f, f' are continuous
                        and f'(x*) ≠ 0. For e < 1 we have
                          |ayn·sin x − axN·cos x| ≤ √(axN²+ayn²) = e < 1,
                        so f' = 1 − (…) > 0 everywhere — root unique, Newton
                        converges quadratically.
  Primary reference:    Hoots & Roehrich (1980) SR3 page 13, eqs. for U, (E+ω),
                        and Δ(E+ω). The SR3 update is exactly the Newton step
                        Δ = (U − ayn·cos x + axN·sin x − x) / (1 − ayn·sin x
                        − axN·cos x).
                        See also: design/derivations/009_sgp4_modified_kepler.md.
  Domain of validity:   e ∈ [0, 1); any U. The Lyddane parameterisation
                        regularises the e=0 singularity of the standard
                        E - e·sin E = M form.

METHOD
  Method declared:      Doc comment line 204 says "Halley's method with cubic
                        convergence". The KeplerFn typedef doc (lines 71-96)
                        says only "modified Kepler equation … iteration", no
                        order specified. SR3 itself prescribes Newton.
  Method implemented:   Lines 240-253: fixed-point Newton iteration
                          x ← x + f/f'
                        with no f'' term. Capped at 30 iterations.
                        Halley's correction 2ff'/(2f'² − f f'') is NOT
                        present. The implementation is plain Newton.
  Match verdict:        ✗ MISMATCH at the doc-comment level. The `@par
                        Kepler Solver` block at line 203-204 advertises
                        "Halley's method with cubic convergence", but the
                        lambda body implements first-order Newton. The
                        SR3 reference (line 84, 233) and the actual algorithm
                        agree on Newton; the Halley wording in the docstring
                        is the outlier.
                        **C-fail: theory citation in code comment disagrees
                        with implemented method.** Either update the comment
                        to "Newton's method (SR3 page 13)" or replace the
                        body with a Halley step including f'' = −ayn·cos x
                        − axN·sin x.

ERROR BOUND
  Bound category:       precision
  Bound formula:        For any contraction-mapping / Newton iteration
                        converged to step |Δ_k| < τ, the Kantorovich /
                        Ostrowski bound gives
                          |x_k − x*| ≤ |Δ_k| / (1 − q)
                        where q < 1 is the local contraction factor. In the
                        quadratic-Newton regime q ≪ 1, so |x_k − x*| ≤ |Δ_k|
                        to leading order. This is the same final-correction
                        bound used for the Halley solver audited in
                        `kepler.md`, and is mandated by REQ-EF-5.
  Bound implemented:    Line 250: `x.errors.precision = x.errors.precision
                          + abs(delta.value);`
                        added only on the converged path (after
                        |Δ| < tolerance).
  Bound verdict:        ⚠ INCOMPLETE.
                        (a) On the converged path: ✓ matches |Δ_final| as
                            mandated by REQ-EF-5.
                        (b) On the **non-converged 30-iteration cap** path
                            (line 254 returns `x` without adding the residual
                            to precision): ✗ silent — no truncation bound
                            recorded, no exception thrown, no flag in
                            `errors.precision`. AUD-EF-4 violation: iterative
                            algorithms that exit without convergence must
                            either throw or add the last-step magnitude
                            (conservatively, ∞ or a large sentinel) to
                            precision. Currently a non-convergent solve
                            returns silently with the same precision budget
                            as a converged one.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (Iterative-algorithm convergence residual
                        added to precision)
  AUD-EF applies:       AUD-EF-4 (Iterative algorithms add residual; must
                        cover both converged and non-converged exits)
  AUD-MC applies:       n/a (numerical root-finder, not an algebra op)
  Verification test:    Recommended:
                          - tests/test_sgp4/ — Kepler solver convergence
                            test at e ∈ {0, 0.1, 0.5, 0.9}; verify reported
                            precision ≥ actual error.
                          - tests/test_sgp4/ — non-convergence test with
                            tolerance below T's ε and 30-iter cap; verify
                            either exception or sentinel precision.

NOTES
  - The standard Kepler form (E − e·sin E = M) is recovered when axN = e,
    ayn = 0, U = M. The modified form's only purpose is to keep the
    iteration well-defined at e → 0 and to absorb long-period corrections
    naturally — algebra-equivalent for e > 0 (good).
  - Starter x₀ = U (line 240). For e small, U ≈ M and this is a 1st-order
    Kepler starter — adequate for Newton's quadratic convergence to reach
    machine precision in ~5-8 iterations at e < 0.3. For e → 1 the starter
    is weak; the 30-iter cap is the only safeguard. Compare with
    `math/kepler.h::solve_kepler` which uses
      E₀ = M + e·sin(M) + (e²/2)·sin(2M)
    (a 2nd-order Taylor in e). The SGP4 modified solver could benefit from
    the same starter, but this is an optimisation, not a correctness issue.
  - The 30-iteration cap should be revisited: Newton at e = 0.99 in the
    modified form can take 15-20 iterations to reach 1e-15. Cap of 30 is
    likely sufficient but unverified; the silent-non-convergence path
    above makes it impossible to detect a near-cap saturation from the
    error budget.
  - The doc-comment Halley/Newton mismatch must be reconciled. The
    cheap fix is the docstring (line 204 → "Newton's method"); the more
    expensive fix is adding f'' and the Halley correction (cubic vs.
    quadratic convergence; one extra trig pair per iter). SR3 page 13
    specifies Newton, so docstring is the right place to fix.
```

---

## 2 `standard_sgp4()::drag_coefficients` lambda

```
=== FORMULA AUDIT CARD ===
ID:                     model_functions::standard_sgp4::drag_coefficients
Location:               src/sgp4/model_functions.h:268-278
Mathematical statement: (none — placeholder that throws)

THEORY                  n/a (placeholder; theory anchor would be
                        Lane (1965), Lane-Hoots (1979) — see drag.md and
                        drag_coefficients.md).
METHOD declared:        "DragCoefficientsFn not yet extracted as injectable
                        lambda."
METHOD implemented:     `throw std::runtime_error(...)`.
Match verdict:          ✓ matched — placeholder behaves exactly as the
                        comment states (no silent zero-fill; explicit
                        throw on misuse). This is the correct
                        not-yet-implemented pattern.
ERROR BOUND
  Bound category:       n/a
  Bound implemented:    n/a (throws before any TrackedValue write)
  Bound verdict:        ✓ matched (no contract to satisfy).
CROSS-AUDIT
  Verification test:    Recommended: a death-test confirming the throw,
                        so future refactors don't silently swap in a
                        no-op that would let the caller use uninitialised
                        s4/qoms24.
NOTES
  - The inline drag-coefficient computation currently lives in
    `near_space.h` (commented at lines 264-267). When extracted, the
    audit card will move to `drag_coefficients.md` (already exists) and
    this stub should be replaced by a forwarder lambda analogous to the
    other slots.
```

---

## 3 Non-numeric content (no audit cards required)

| Item | Lines | Role | Audit |
|---|---:|---|---|
| `SecularRatesFn`, `InclinationFn`, `KeplerFn`, `SiderealTimeFn`, `DragCoefficientsFn` typedefs | 45-130 | `std::function` aliases — pure types, no numerics | n/a |
| `ModelFunctions<T>` struct definition | 147-184 | Five `std::function` fields + `description` string | n/a |
| `standard_sgp4()::secular_rates` lambda | 212-221 | One-line forwarder to `perturbation::compute_secular_rates` | audited in `brouwer.md` |
| `standard_sgp4()::inclination_function` lambda | 223-229 | One-line forwarder to `perturbation::inclination_function` | audited in `kaula.md` |
| `standard_sgp4()::sidereal_time` lambda | 260-262 | One-line forwarder to `astronomy::compute_gmst` | audited in `sidereal_time.md` (out of scope here; not yet present in `theoretical_basis_audit/`) |
| `mf.description` string assignment | 280-281 | Documentation literal | n/a |

The factory itself (`standard_sgp4()`, lines 209-283) is a configuration
constructor, not a numerical computation. The expected function count of 1 in
the task brief is correct **at the factory granularity** (1 factory, returning
1 bundle); at the formula granularity inside the factory body, 1 numerical
formula card is required (`kepler_solver`, §1 above) and 1 placeholder notice
(`drag_coefficients`, §2 above).

---

## 4 File-level verdict

- **A. Error wiring**: ⚠ INCOMPLETE for `kepler_solver` — converged path adds
  |Δ_final| to `errors.precision` (✓ REQ-EF-5), but the 30-iteration
  non-convergence path (line 254 fall-through) returns silently with no
  precision-budget penalty. AUD-EF-4 violation.

- **B. Algebra axioms**: n/a — no algebraic identities live in this file.

- **C. Theoretical basis**:
  - `kepler_solver`: ✗ **C-fail on declared method**. Doc comment line 204
    says "Halley's method with cubic convergence"; the lambda body
    (lines 240-253) implements first-order Newton. SR3 page 13 prescribes
    Newton, so the implementation matches SR3 — the docstring is the wrong
    artefact. Bound is the right Kantorovich form on the converged path
    (✓) but missing on the non-converged path (⚠).
  - `drag_coefficients`: ✓ behaves as the placeholder docstring states
    (throws). No formula to validate.
  - Forwarder lambdas (`secular_rates`, `inclination_function`,
    `sidereal_time`): out of scope here; their underlying implementations
    have their own audit cards.

**File verdict: C-FAIL (low severity)** — the file is mostly plumbing, but
the `kepler_solver` lambda has (i) a docstring "Halley" claim that doesn't
match the Newton implementation, and (ii) a silent non-convergence exit path
that drops the REQ-EF-5 contract. Both are localised fixes:

1. **Reconcile docstring with code.** Change line 203-204 from "Halley's
   method with cubic convergence" to "Newton's method (SR3 page 13);
   quadratic convergence in a neighbourhood of the root for e < 1".
2. **Cover the non-convergence path.** Either `throw` after the loop, or
   add `abs(delta.value)` (the last-computed step) to `x.errors.precision`
   before returning at line 254. The latter degrades silently but keeps
   the precision budget honest; the former forces callers to handle
   non-convergence explicitly.

Once those two edits land, the file passes TBA.

---

## 5 Cross-references

- `theoretical_basis_audit.md` §1 (audit-card template), §5 (worked example),
  §6 (`math/kepler.h::solve_kepler` audit card — comparator for this file's
  in-line Kepler solver)
- `design/derivations/009_sgp4_modified_kepler.md` — derivation of the SR3
  modified form
- `kepler.md` — standard `E - e·sin E = M` solver (Halley) audit card; this
  file's SGP4 Kepler is the Newton sibling, not the same algorithm
- `brouwer.md`, `kaula.md`, `sidereal_time.md` — audit cards for the
  forwarder slots
