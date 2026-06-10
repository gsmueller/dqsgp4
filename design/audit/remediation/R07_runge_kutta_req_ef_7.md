# R07 — runge_kutta.h integrator-order accuracy

**Status:** DONE — commit 3cc618f
**Severity:** P2 (systematic REQ-EF-7 wiring)
**Estimated scope:** ~15 LoC

---

## Files

**Write:**
- `src/integrators/runge_kutta.h` (both `euler` and `runge_kutta_4`)

**Read:**
- `src/dynamics/state.h`, `src/dynamics/derivative.h`
- `design/specifications/error_framework.md` REQ-EF-7
- Butcher (2003) §1.2 (Euler LTE) and §2.2 (RK4 LTE) — for bound formulas

**Audit card:**
- `design/audit/theoretical_basis_audit/runge_kutta.md`

---

## Primary issue

`euler` and `runge_kutta_4` are correctly implemented in Munthe-Kaas Lie-group form. Their per-step `errors.precision` is correctly propagated through `lie_advance_pose` + field arithmetic.

**Missing:** the per-step local truncation error (LTE) is NOT added to `errors.accuracy`:
- Euler: O(δt²) — first-order method, LTE ~ ½|y''|·δt²
- RK4: O(δt⁵) — fourth-order method, LTE ~ (1/120)|y⁽⁵⁾|·δt⁵

REQ-EF-7 spec violation. The audit cards specifically flagged this for both functions.

---

## Theory anchor

Butcher (2003) "Numerical Methods for Ordinary Differential Equations":

- **§1.2** — Forward Euler local truncation: y_{n+1} = y_n + δt·f(t_n, y_n) + O(δt²). Rigorous bound: |LTE| ≤ (δt²/2)·sup|y''(ξ)| over the step.

- **§2.2 / Eq. 2.2.5** — Classical RK4 local truncation: O(δt⁵). Rigorous bound: |LTE| ≤ (δt⁵/120)·M_5, where M_5 bounds the fifth time derivative of the state over the step.

For practical bounds in orbital mechanics, the second/fifth derivatives can be estimated from the force lambda magnitudes. A conservative bound: use the previous step's acceleration magnitude as a proxy for higher derivatives.

REQ-EF-7: model-truncation bound added to `errors.accuracy` at each integrator step.

---

## Fix scope

~15 LoC across `euler` and `runge_kutta_4`. Approximate approach:

```cpp
// In euler:
auto |accel| = sqrt(acceleration.linear.x*acceleration.linear.x + ...);
auto lte_bound = T(0.5) * |accel.value| * dt.value * dt.value;
result.pose.position.x.errors.accuracy += lte_bound;
// (... apply to all state components)

// In runge_kutta_4:
auto |accel| = ...; // same magnitude estimate
auto lte_bound = (|accel.value|/T(120)) * pow(dt.value, 5);
result.pose.position.x.errors.accuracy += lte_bound;
```

The exact form of the bound estimate depends on what the agent considers sufficiently rigorous. The Butcher closed forms above are the theoretical statement; an implementable estimator approximates `M_5` by recent acceleration magnitudes.

---

## Verification

1. Rebuild: `build.bat nodocs`
2. Run `test_propagator.exe` — Phases 1-3 still PASS (regression).
3. Manual check: propagate a known closed-form orbit (Kepler) for 1000 RK4 steps; observe `errors.accuracy` field grows ≈ quartically with step count (RK4 LTE) or linearly with step count for cumulative bound.
4. Compare Euler vs RK4 accuracy growth rates — Euler should grow as O(t·δt), RK4 as O(t·δt⁴).

---

## References

- Audit card: `design/audit/theoretical_basis_audit/runge_kutta.md`
- Primary sources: Butcher (2003), Runge (1895), Kutta (1901), Munthe-Kaas (1999)

---

## Status history

- 2026-05-13 — Created from approved plan `peppy-lobster`. OPEN.
