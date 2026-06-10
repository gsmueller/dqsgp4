# R11 — taylor_half_angle_scale wide-T correction

**Status:** DONE — commit 3cc618f
**Severity:** P3 (tightening for wide T; tight for double already)
**Estimated scope:** ~5 LoC + optional T-dependent threshold

---

## Files

**Write:**
- `src/math/small_angle_series.h` (function `taylor_half_angle_scale`, around line ~80)

**Read:**
- `design/audit/theoretical_basis_audit.md` §5.2 (audit card with the analysis)
- Abramowitz-Stegun (1964) §15.1.10 (`arcsin(s)/s` series)

**Audit card:**
- `design/audit/theoretical_basis_audit/small_angle_series.md` §5.2

---

## Primary issue

The Taylor series of `arcsin(s)/s` at `s = 0`:

```
arcsin(s)/s = 1 + s²/6 + 3s⁴/40 + 5s⁶/112 + 35s⁸/1152 + ...
```

has all-**positive** coefficients (NOT alternating). The current implementation truncates after `3s⁴/40` and bounds the truncation with `5|s|⁶/112` — the magnitude of the next term.

**Problem:** Leibniz's bound (next-term magnitude) applies only to **alternating** series. For positive-coefficient series, the tail is technically larger than just the next term. The implementation under-counts by the ratio of the geometric tail to the leading term, which is `1/(1 − 2s²) − 1`:

- For `s = 10⁻⁴` (threshold): correction factor is `1 / (1 − 2×10⁻⁸) ≈ 1 + 2×10⁻⁸` — undercount is ~2×10⁻⁸ × leading-coefficient × |s|⁶ = ~10⁻³² (negligible for double)
- For wider T (`cpp_bin_float_50`): the threshold needs to scale, and the correction matters

For double the existing bound is still rigorous (the correction is below double's ε), but for arbitrary precision T the bound becomes unsound at sufficiently wide precision.

---

## Theory anchor

Abramowitz-Stegun (1964) §15.1.10. The series `arcsin(s)/s = Σ_{n≥0} ((2n)! / (4^n · (n!)² · (2n+1))) · s^{2n}` has all-positive coefficients with ratio of successive terms ≤ `2·s²`. Hence the geometric tail estimate:

```
|R_N(s)| ≤ |c_{N+1}| · s^{2(N+1)} · 1/(1 − 2s²)    for |s| < 1/√2
```

The current code uses just `|c_{N+1}| · s^{2(N+1)}` — the leading term. The fix is to multiply by `1/(1 − 2s²)`.

---

## Fix scope

~5 LoC. In `taylor_half_angle_scale`, after computing the existing `tail_bound`:

```cpp
auto s2 = qv_norm_sq;
auto correction = exact<T>(1) / (exact<T>(1) - exact<T>(2) * s2);
tail_bound = tail_bound * correction;
result.errors.precision = result.errors.precision + tail_bound;
```

Also optionally:
- Make the small-argument threshold T-dependent: `threshold = pow(T(eps_of_T), 1/6)` so the bound stays tight across types.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_math.exe`, `test_quaternion.exe`, `test_dual_quaternion.exe` — all PASS (regression).
3. Add wide-T test (per REQ-SY-7 — `T = boost::multiprecision::cpp_bin_float_50`):
   - Compute `taylor_half_angle_scale(s, s², 1)` at `s = 1e-10` for `T = cpp_bin_float_50`
   - Compare reported bound against actual error against reference high-precision value
   - Bound should be tight (≤ 2× actual error) but rigorous (≥ actual error)

---

## References

- Audit card: `design/audit/theoretical_basis_audit/small_angle_series.md` §5.2
- Audit framework: `design/audit/theoretical_basis_audit.md` §5.2 (same content)
- Primary sources: Abramowitz-Stegun (1964), Knopp (1928)

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
