# R12 — Modified-Kepler + secular-update boundary handling

**Status:** DONE — commit 3cc618f
**Severity:** P3 (boundary-case rigor gap)
**Estimated scope:** ~10 LoC across two files

---

## Files

**Write:**
- `src/orbit/modified_kepler.h` (both `solve_kepler_newton` and `solve_kepler_halley`)
- `src/orbit/secular_update.h` (the e<1e-6 floor branch)

**Read:**
- `src/math/kepler.h:79` (established pattern for residual addition)
- `design/specifications/error_framework.md` REQ-EF-5, REQ-EF-9

**Audit cards:**
- `design/audit/theoretical_basis_audit/modified_kepler.md`
- `design/audit/theoretical_basis_audit/secular_update.md`

---

## Primary issues

Two related "boundary case relaxes the rigorous bound" issues, both in `src/orbit/`:

### 1. modified_kepler.h cap-hit fallback

Both `solve_kepler_newton` and `solve_kepler_halley` have a 30-iteration safety cap. On cap-hit (non-convergence within 30 iters), the function returns the current best estimate **without adding `|Δ_last|` to `errors.precision`**. AUD-EF-4 violation.

This silently passes through a potentially-unconverged result as if it were converged.

### 2. secular_update.h eccentricity floor

`secular_advance` floors the eccentricity at `e = 1e-6` (to avoid division-by-zero downstream). When the floor fires, the returned `e` value carries a non-bounded systematic offset:

```
true_e ∈ [0, ~1e-6]  but  returned_e = 1e-6
```

This silently relaxes the rigorous-bound contract: `total_error()` on the returned `e` does not include this systematic offset.

---

## Why coupled

Both are "boundary case silently relaxes the rigorous bound" issues. Both live in `src/orbit/`. Both need the same kind of fix: when the guard fires, record the discrepancy in the appropriate error category.

---

## Theory anchor

REQ-EF-5: convergence residual → precision. Applies to both Kepler solvers.

REQ-EF-9: "Catastrophic loss is signaled, not silenced." Applies to the e-floor — when the floor fires, the result is in a "catastrophic" regime and must be flagged.

---

## Fix scope

~10 LoC across the two files.

### modified_kepler.h

Both solvers — add at end of iteration loop (cap-hit fallback path):

```cpp
// After the for loop completes without convergence
result.errors.precision = result.errors.precision + abs(last_correction);
return result;
```

### secular_update.h

When the floor fires:

```cpp
if (e_value < T(1e-6)) {
    e_clamped = TrackedValue<T>(T(1e-6), ...);
    // Add the offset to errors.accuracy (this is a model approximation, not numerical):
    e_clamped.errors.accuracy = e_clamped.errors.accuracy + T(1e-6);
}
```

The bound `T(1e-6)` is the full range of the offset.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_sgp4.exe`, `test_math.exe` — all PASS (regression).
3. Test the cap-hit path: construct a stressed Kepler input (e ≈ 0.99 with high mean anomaly); verify `errors.precision > 0` even on cap-hit.
4. Test the e-floor path: construct a low-eccentricity TLE; advance until e crosses 1e-6; verify `result.e.errors.accuracy >= 1e-6` after the floor fires.

---

## References

- Audit cards: `design/audit/theoretical_basis_audit/modified_kepler.md`, `secular_update.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` §"Open notes (PASS-with-notes)"

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
