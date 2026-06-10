# R08 — element_recovery F2 series tail accuracy

**Status:** DONE — commit 3cc618f
**Severity:** P2 (systematic REQ-EF-7 wiring)
**Estimated scope:** ~5 LoC

---

## Files

**Write:**
- `src/orbit/element_recovery.h`

**Read:**
- `design/derivations/deprecated/013_element_recovery.md` (if exists — Brouwer-Lyddane derivation)
- `sgp4_references/Hoots_Roehrich_1980_Spacetrack_Report_No3.pdf` §6 (element recovery)
- `design/specifications/error_framework.md` REQ-EF-7

**Audit card:**
- `design/audit/theoretical_basis_audit/element_recovery.md`

---

## Primary issue

The Brouwer-Lyddane series reversion in `recover_mean_elements` evaluates the polynomial:

```
a₀'' = a₁ · (1 - δ₁/3 - δ₁² - (134/81)·δ₁³)
```

This is the order-3 truncation of an infinite series in `δ₁`. The dropped term is O(δ₁⁴), which for realistic TLEs (δ₁ ~ 10⁻³) is ~10⁻¹². Negligible numerically but not added to `errors.accuracy` — REQ-EF-7 gap.

---

## Theory anchor

Brouwer (1959) Theorem on element recovery; SR3 §6 lines 1340-1370 (Vallado's `initl` chain).

The series is the inverse of the (a → a₀) Brouwer correction; reversion to order 3 in δ₁ retains all terms up to O(J₂³). Higher orders give corrections to a₀'' that are O(J₂⁴) ≈ O(10⁻¹³) relative.

REQ-EF-7: model-truncation bound added to `errors.accuracy`. The bound for the O(δ₁⁴) tail is:

```
|tail| ≤ C · |δ₁|⁴ · a₁    where C ≈ 1 (next-coefficient magnitude)
```

---

## Fix scope

~5 LoC. After computing `a0_pp`:

```cpp
auto del1_abs = abs(del1.value);
auto tail_bound = pow(del1_abs, 4) * a1.value;
result.a0.errors.accuracy = result.a0.errors.accuracy + tail_bound;
// Similar propagation to n0 (since n0 derives from a0 via Kepler 3rd law):
//   n0.errors.accuracy gets the corresponding bound via the Jacobian
```

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_sgp4.exe`:
   - All 8 near-earth cases still PASS at machine precision (regression)
   - Recovered `n0.errors.accuracy` is now nonzero (~10⁻¹³ relative to `n0.value`)
3. Manual: print `ns.n0.errors.accuracy / ns.n0.value` after init for Sat 00005; expect ~10⁻¹²..10⁻¹³.

---

## References

- Audit card: `design/audit/theoretical_basis_audit/element_recovery.md`
- Primary sources: Brouwer (1959); SR3 §6
- Pattern: matches the same series-truncation-→ -accuracy approach as R09 (drag coefficients)

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
