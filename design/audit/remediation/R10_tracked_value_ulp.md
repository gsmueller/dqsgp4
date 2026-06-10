# R10 — tracked_value representation_bound 0.5 ULP tightening

**Status:** ALREADY-CORRECT — file already implements 0.5 ULP at line 441 (verified by inspection; no change required)
**Severity:** P3 (rigor tightening; bound currently conservative-but-rigorous)
**Estimated scope:** ~3 LoC

---

## Files

**Write:**
- `src/math/tracked_value.h` (line ~441 in `representation_bound`)

**Read:**
- IEEE-754 §3.4 (rounding rules)
- Goldberg (1991) "What Every Computer Scientist Should Know About Floating-Point Arithmetic" §2 (rounding error)

**Audit card:**
- `design/audit/theoretical_basis_audit/tracked_value.md` (the re-dispatched audit notes this specifically)

---

## Primary issue

`representation_bound` returns **1 ULP** of the input value, but IEEE-754 round-to-nearest-even gives at most **0.5 ULP** of rounding error. The bound is conservative by a factor of 2 — rigorous (does not under-bound), but unnecessarily loose.

Quote from the audit card:
> `representation_bound` (card #40) returns 1 ULP rather than 0.5 ULP for double/float/multiprec, making it conservative by factor 2; rigorous but loose.

---

## Theory anchor

IEEE-754-2019 §3.4 (rounding rules): round-to-nearest-even guarantees that the result is within ±0.5 ULP of the true value.

Goldberg (1991) §2: "if x is a floating-point number that is a result of an arithmetic operation, then `|x - true| ≤ 0.5 ulp(x)`."

---

## Fix scope

~3 LoC. In `representation_bound`:

```cpp
// Before: returns ulp_value
// After:
return ulp_value / T(2);
```

The exact change depends on how the existing code computes ULP. Verify against the audit-flagged code path (around line 441). The 0.5 factor must be applied for all types (float, double, long double, boost multiprec).

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_math.exe` — all PASS (regression). Tests that compare bounds against actual roundoff should now show the bound is tighter but still rigorous (actual error ≤ bound).
3. Manual: `representation_bound(1.0/3.0)` for double should now return exactly `2⁻⁵³ · 1/3 ≈ 3.7e-17` (was `2⁻⁵² · 1/3 ≈ 7.4e-17`).
4. Edge case: `representation_bound(0.0)` — confirm behavior at zero is still correct (likely returns 0 or smallest subnormal).

---

## References

- Audit card: `design/audit/theoretical_basis_audit/tracked_value.md`
- Consolidated summary: `design/audit/AUD_TBA_results.md` mention of "1 ULP rather than 0.5 ULP"
- Primary sources: IEEE-754-2019, Goldberg (1991), Higham (2002) "Accuracy and Stability of Numerical Algorithms"

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
