# R03 — model_functions docstring + AUD-EF-4

**Status:** DONE — commit 3cc618f
**Severity:** P1
**Estimated scope:** ~10 LoC (docstring + AUD-EF-4 patch)

---

## Files

**Write:**
- `src/sgp4/model_functions.h` (lines ~204 for docstring; ~254 for AUD-EF-4 patch)

**Read:**
- `src/math/kepler.h` (line 79 — pattern for residual addition on convergence)
- `sgp4_references/Hoots_Roehrich_1980_Spacetrack_Report_No3.pdf` §6 (modified Kepler)
- `design/audit/error_framework.md` AUD-EF-4

**Audit card:**
- `design/audit/theoretical_basis_audit/model_functions.md`

---

## Primary issue

Two issues in `model_functions.h::standard_sgp4::kepler_solver` lambda:

1. **Doc-code mismatch (docstring bug):** Docstring at line ~204 advertises *"Halley's method with cubic convergence"* but the lambda body at line ~234 implements first-order Newton (no `f''` term — only `f` and `f'`).

2. **AUD-EF-4 violation:** 30-iteration non-convergence fall-through at line ~254 returns silently. AUD-EF-4 requires every iterative algorithm to record `|Δ_final|` in `errors.precision` on every exit path.

---

## Theory anchor

SR3 §6 (modified Kepler equation): the modified Kepler equation `x - ayn·cos(x) + axN·sin(x) = U` is conventionally solved by Newton-Raphson (SR3's algorithm). The current code is correct for Newton; only the docstring is misleading.

REQ-EF-5: convergence residual → precision. AUD-EF-4: every iterative algorithm must satisfy REQ-EF-5 on every exit path.

---

## Fix scope

**Decision required:** Either upgrade the lambda to Halley to match the docstring, OR update the docstring to say Newton (matches SR3 and code).

**Recommended choice:** update docstring to "Newton" — this is what SR3 §6 prescribes and what the code actually does. Upgrading to Halley introduces new test coverage burden for no observable accuracy benefit at double precision.

Also patch the cap-hit path:

```cpp
// At line ~254 (after for loop completes without convergence)
result.errors.precision = result.errors.precision + abs(last_correction);
return result;
```

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_sgp4.exe`:
   - 8/8 near-earth cases PASS at machine precision (regression — no behavior change)
   - All position errors < 1e-7 km at t=0
3. Manual: stress the cap-hit path with an artificially tight tolerance or high-eccentricity input; observe `errors.precision > 0` in the result.

---

## References

- Audit card: `design/audit/theoretical_basis_audit/model_functions.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` §"PASS-with-doc-bug"

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
