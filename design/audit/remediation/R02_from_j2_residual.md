# R02 — from_J2 Banach fixed-point residual

**Status:** DONE — commit 3cc618f
**Severity:** P0 (blocking C-FAIL #3, but trivial scope)
**Estimated scope:** ~5 LoC, ~10 minutes

---

## Files

**Write:**
- `src/geodesy/equipotential_ellipsoid.h` (the `from_J2` static factory, around line ~256)

**Read:**
- `src/math/kepler.h` (line 79 — established pattern for iteration-residual handling)
- `src/math/series.h` (line 147 — same pattern in `series_sqrt`)
- `design/specifications/error_framework.md` REQ-EF-5

**Audit card:**
- `design/audit/theoretical_basis_audit/equipotential_ellipsoid.md`

---

## Primary issue (C-FAIL #3)

The `from_J2` static factory performs a Banach fixed-point iteration on `1/f` given `J₂`, converging to the value of `e²` that matches the geodetic relation between `J₂` and ellipsoid flattening. The iteration converges fine.

**Problem:** the final correction magnitude is **NOT** added to `e2_guess.errors.precision` before the value is forwarded to the WGS84 ctor. REQ-EF-5 ("convergence residual is added to precision") spec violation.

**Impact:** every constant derived through `from_J2` (WGS72, GRS80 paths) carries an under-counted precision bound. Magnitude is type-epsilon-level (negligible for double; matters for `cpp_bin_float_50+`). Structurally a spec violation — the rigorous-bound contract for `total_error()` is broken on these paths.

---

## Theory anchor

REQ-EF-5 from `design/specifications/error_framework.md`:

> When an iterative algorithm (Newton, Halley, fixed-point) terminates at iteration k with correction Δ_k below the caller's tolerance, the result's `errors.precision` is incremented by |Δ_k|:
>
>     result.errors.precision += abs(delta_k);

---

## Fix scope

Add the residual update at convergence in `from_J2`. The exact placement matches `kepler.h:79`:

```cpp
// After: if (abs(correction.value) < tolerance) break;
e2_guess.errors.precision = e2_guess.errors.precision + abs(correction.value);
```

Place this immediately before the function returns or breaks out of the iteration loop.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_wgs84.exe` and `test_geodesy.exe` — both PASS (regression).
3. Manual check: print `constants.earth.J2n(1).errors.precision` after `constants::ConstantsProvider<double>::wgs72(1e-15)`:
   - Before fix: 0
   - After fix: ~1e-15 (matches tolerance)

---

## References

- Audit card: `design/audit/theoretical_basis_audit/equipotential_ellipsoid.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` §"C-FAIL #3"

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
