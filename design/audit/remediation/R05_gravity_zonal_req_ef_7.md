# R05 — gravity_zonal.h REQ-EF-7 J₃+ truncation accuracy

**Status:** DONE — commit 3cc618f
**Severity:** P2 (systematic REQ-EF-7 wiring)
**Estimated scope:** ~10 LoC

---

## Files

**Write:**
- `src/forces/gravity_zonal.h`

**Read:**
- `src/dynamics/wrench.h` (return type structure)
- `src/constants/constants_provider.h` (`K.earth.J2n(n)` lookup)
- `design/specifications/error_framework.md` REQ-EF-7

**Audit card:**
- `design/audit/theoretical_basis_audit/gravity_zonal.md`

---

## Primary issue

`gravity_J2` correctly implements the J₂ Legendre derivative for the central-body oblateness acceleration. Its `errors.precision` is correctly propagated through the closed-form arithmetic.

**Missing:** the model-truncation residual (the omitted J₃, J₄, J₅, ... contributions) is NOT declared into `errors.accuracy`. REQ-EF-7 ("Force models add the model bound to accuracy") spec violation.

---

## Theory anchor

REQ-EF-7 from `design/specifications/error_framework.md`:

> When a model is deliberately simplified (J₂-only gravity vs full spherical harmonics, exponential drag vs an empirical density model, two-body astronomy vs JPL DE ephemeris, etc.), the model-truncation bound is added to `errors.accuracy`. Force lambdas in `src/forces/` declare their model-truncation contribution; the propagator accumulates it into the state at each step.

Bound formula: the J₃ contribution magnitude at the satellite altitude. For Earth at LEO:

```
|a_J3| ≈ |J3 · μ · a_E³ / r⁵ · (5z²/r² - 3)|
      ≈ J3 / J2 · (a_E/r) · |a_J2|
      ≈ 10⁻⁶ × |a_J2|  (for LEO)
```

So the bound is `|a_J3| ≤ (J3/J2)·(a_E/r)·|a_J2|` per acceleration component.

Higher zonals (J₄, J₅, ...) contribute O(10⁻⁹) relative to J₂; absorbable into the J₃ bound with factor 1.1× safety.

---

## Fix scope

~10 LoC. Compute and add the J₃-equivalent truncation bound to each linear component of the returned wrench:

```cpp
// In gravity_J2, after computing a_J2_x, a_J2_y, a_J2_z:
auto J3_over_J2 = abs(K.earth.J2n(2).value) / abs(K.earth.J2n(1).value);
auto aE_over_r = K.earth.a.value / r_magnitude.value;
auto safety = T(1.1);

auto trunc_x = safety * J3_over_J2 * aE_over_r * abs(a_J2_x.value);
result.linear.x.errors.accuracy = result.linear.x.errors.accuracy + trunc_x;
// same for y, z
```

(Pseudo-code — adapt to actual variable names and TrackedValue API.)

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_propagator.exe` — Phases 1-3 still PASS (regression).
3. Manual check: propagate a 400 km circular orbit one revolution under J₂ only; `errors.accuracy` field on position should be ~10⁻⁶ of the position value (J₃-equivalent bound).
4. Compare against the precision contribution: accuracy should dominate after wiring is in.

---

## References

- Audit card: `design/audit/theoretical_basis_audit/gravity_zonal.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` §"Cross-cutting observations" #2
- Primary source: Heiskanen & Moritz (1967) "Physical Geodesy" §2

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
