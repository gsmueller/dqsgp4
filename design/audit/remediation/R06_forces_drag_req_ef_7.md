# R06 — forces/drag.h REQ-EF-7 atmospheric-model truncation

**Status:** DONE — commit 3cc618f
**Severity:** P2 (systematic REQ-EF-7 wiring)
**Estimated scope:** ~10 LoC

---

## Files

**Write:**
- `src/forces/drag.h` (DQ propagator path — `make_drag_exponential`)

**Note:** this is the DQ propagator's drag lambda. Do NOT touch `src/atmosphere/drag_coefficients.h` (which is the SGP4 Lane-Hoots path — separate scope under R09).

**Read:**
- `src/dynamics/wrench.h`
- `design/specifications/error_framework.md` REQ-EF-7
- `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §2 Lemma 2.2 (Lane density model accuracy)

**Audit card:**
- `design/audit/theoretical_basis_audit/drag.md`

---

## Primary issue

`make_drag_exponential` correctly implements the Lane exponential-atmosphere model with closed-form per-step density evaluation. `errors.precision` is correctly propagated.

**Missing:** the ~30% baseline accuracy floor of the Lane model is NOT annotated as `errors.accuracy` on the produced wrench. REQ-EF-7 spec violation.

---

## Theory anchor

From `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §2 Lemma 2.2 (Lane power-law approximation property):

> The form (2.1) with τ=4, ρ₀, q₀, s as above matches the U.S. Standard Atmosphere within ~30% over altitudes 100-500 km.

Per the same document §15 ERROR SOURCE A-D2:
- Baseline density-model error: ~30%
- Solar maxima: up to ~100%
- This is the dominant accuracy limit of any Lane-based drag prediction beyond a few revolutions.

REQ-EF-7: model-truncation bound added to `errors.accuracy`. Force lambdas in `src/forces/` declare their model-truncation contribution; the propagator accumulates it into the state at each step.

---

## Fix scope

~10 LoC. Add to each linear component of the produced wrench:

```cpp
// After computing drag_accel:
auto density_accuracy_floor = T(0.30);  // 30% baseline per Lane (1965)
auto truncation_x = density_accuracy_floor * abs(drag_accel.x.value);
result.linear.x.errors.accuracy = result.linear.x.errors.accuracy + truncation_x;
// same for y, z
```

Consider whether to expose the `0.30` constant as a configurable parameter (for switching between "baseline" 0.3 and "solar-max" 1.0 regimes). For now, the conservative baseline is appropriate.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_propagator.exe` — Phase 3 (drag, 5 periods) still PASS (regression).
3. Manual check: Phase 3 results' `errors.accuracy` field should now dominate `errors.precision` after ~1 day of propagation.
4. The current Phase 3 reports 709 km position offset after 5 periods; the accuracy bound on that should be ~30% × 709 km ≈ 210 km after the fix.

---

## References

- Audit card: `design/audit/theoretical_basis_audit/drag.md`
- SGP4 drag derivation: `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` §2, §15
- Primary source: Lane (1965) "The Development of an Artificial Satellite Theory Using a Power-Law Atmospheric Density Model"

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
