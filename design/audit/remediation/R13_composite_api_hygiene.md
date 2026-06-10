# R13 — Composite-API hygiene sweep

**Status:** DONE — commit 3cc618f
**Severity:** P4 (API smell / annotation hygiene)
**Estimated scope:** ~20 LoC across three files

---

## Files

**Write:**
- `src/sgp4/state_vector.h`
- `src/dynamics/derivative.h`
- `src/atmosphere/density_model.h`

**Read:**
- `design/specifications/error_framework.md` REQ-EF-1, REQ-EF-2
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §2 ERROR SOURCE A-D4 (for density_model)

**Audit cards:**
- `design/audit/theoretical_basis_audit/state_vector.md`
- `design/audit/theoretical_basis_audit/derivative.md`
- `design/audit/theoretical_basis_audit/density_model.md`

---

## Primary issues

Three independent "API hygiene" annotations, bundled because each touches one small region in a different file.

### 1. state_vector.h::position_error / velocity_error return type

`position_error()` and `velocity_error()` compute RSS over Cartesian components. They currently return `TrackedValue<T>` with all error fields zeroed — silently discarding the sqrt's own propagated bound via the `.value` unwrap.

**Fix decision:** either
- Return plain `T` (matching the per-category accessors `position_measurement_error`, `position_precision_error`, `position_accuracy_error` at the same line range), or
- Return a real `TrackedValue<T>` with the sqrt bound propagated correctly.

Recommend the former (plain `T`) for consistency with the per-category accessors and because callers typically want a scalar magnitude, not an error-budget-on-an-error-budget.

### 2. derivative.h::zero() not the additive identity

`zero()` sets `time_rate = 1` not `0`. This is NOT the group additive identity for `Derivative<T>` (which would be `(Twist::zero(), TrackedValue::exact_integer(0))`).

**Fix decision:** either
- Rename to something like `standard_time_zero_accel()` or `identity_time_rate_zero_accel()`, or
- Fix to `time_rate = 0` (and add a separate `from_acceleration()` factory if `time_rate = 1` is needed for RK4 stage combinations).

Recommend rename: the function's purpose IS "zero acceleration with default time parameterization"; renaming clarifies intent without breaking the RK4 stage-combination semantics.

### 3. density_model.h C⁰ continuity annotation

The perigee adjustment at 156 km transitions Cases A → B; `s` and `(q₀-s)⁴` are computed differently on either side. The transition is `C⁰` continuous (the function values match) but NOT `C¹` (the derivatives differ). Per `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §2 ERROR SOURCE A-D4 this introduces a ~1% drag-rate discontinuity near the boundary.

**Fix:** add an in-code comment cross-referencing A-D4 + an explicit `errors.accuracy` annotation when the boundary is straddled or near.

---

## Theory anchor

REQ-EF-1 ("Every tracked value carries three error budgets"), REQ-EF-2 ("`total_error()` is a rigorous upper bound"), `sgp4_near_earth_drag_theoretical_basis.md` §2 ERROR SOURCE A-D4.

---

## Fix scope

~20 LoC across the three files. Three independent edits.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run all existing tests — PASS (regression). The state_vector return-type change may require minor caller updates in `tests/test_sgp4`.
3. Manual inspection of `state_vector` usage in callers (probably just `test_sgp4`).
4. Add unit test that exercises the density_model 156-km boundary, observing the accuracy field annotation.

---

## References

- Audit cards: `state_vector.md`, `derivative.md`, `density_model.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` §"Open notes (PASS-with-notes)"

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
