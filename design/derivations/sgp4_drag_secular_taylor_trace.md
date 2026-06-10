# Phases 6-8 — SGP4 Secular-Drag Taylor Coefficients D₂, D₃, D₄, t2cof…t5cof

## §D.0 Scope and framework

Traces the SGP4 higher-order secular-drag coefficients D₂, D₃, D₄
(`src/atmosphere/drag_coefficients.h:178-182`) under **Standard 10**
(`simplify(D_k_derived − D_k_code) == 0`). They are the `t²,t³,t⁴` Taylor coefficients of
the drag-decay law
```
a(t) = a₀ · tempa² ,   tempa = 1 − C₁t − D₂t² − D₃t³ − D₄t⁴       (drag_coefficients.h:14)
```

**The decay ODE.** Generalising the epoch secular law (Phase 2.A §A.9, `ȧ(0) = −2a₀C₁`) to
the current `a`, the semi-major axis obeys
```
ȧ = −2 a · C₁(a) ,    C₁(a) = B* C₂(a) .                          (D.0.1)
```
**The a-dependence of C₁ (load-bearing).** To leading order in `η` (η→0), `C₂ ≈ coef1·n·a`
with `coef1 ∝ ξ⁴`, `ξ = 1/(a−s)`, the **mean motion `n ∝ a^{−3/2}`** (Kepler), and the
explicit `a`. Hence
```
C₁(a) ∝ (a−s)^{−4} · a^{−3/2} · a = (a−s)^{−4} · a^{−1/2} ,
  d(ln C₁)/da = −4ξ − 1/(2a) .                                    (D.0.2)
```
The `−1/(2a)` (from `n ∝ a^{−3/2}` and the explicit `a`) is essential — **dropping it is the
root cause of audit finding D-5** (§D.2). The η/ψ a-dependence is a higher-order
(AFGP4→SGP4) truncation, not code-used.

**Taylor match.** Writing `a(t)/a₀ = 1 + b₁t + b₂t² + b₃t³ + b₄t⁴` from (D.0.1) via the chain
rule (`a' = g(a)`, `g = −2aC₁(a)`; `a'' = g'g`, `a''' = g''g²+g'²g`, `a'''' = g'''g³+4g''g'g²+g'³g`),
and `tempa² = 1 − 2C₁t + (C₁²−2D₂)t² + (2C₁D₂−2D₃)t³ + (D₂²+2C₁D₃−2D₄)t⁴`, equate coefficients:
```
D₂ = (C₁² − b₂)/2 ,   D₃ = C₁D₂ − b₃/2 ,   D₄ = (D₂² + 2C₁D₃ − b₄)/2 .   (D.0.3)
```

## §D.2 D₂ — close D-5  *(verify_D2.m 3/3)*

With the full `C₁(a) ∝ (a−s)^{−4}a^{−1/2}` (D.0.2), the `t²` match (D.0.3) gives
```
D₂ = 4 a₀ ξ C₁²                                                  (D.2.1)  [drag_coefficients.h:178]
```
**EXACTLY** — `verify_D2.m` D2.1/D2.2, `simplify(D₂ − 4a₀ξC₁²) = 0`. The stray `C₁²` (from the
`1` in `[1 + a₀ d(ln C₁)/da]` and the `−1/2`) cancels against `tempa²`'s `C₁²`.

**D-5 resolution.** The audit's spurious residual `D₂ = −½C₁² + 4a₀ξC₁²` is reproduced
**exactly** by dropping `n(a)` (using `C₁ ∝ (a−s)^{−4}` only, `d(ln C₁)/da = −4ξ`) —
`verify_D2.m` D2.3. So D-5 was an **omission of the mean-motion a-dependence** in the audit's
re-derivation, not a derivation error or a truncation in the code. The code `D₂ = 4a₀ξC₁²` is
correct and exact (to the order SGP4 keeps).

## §D.3 D₃, D₄  *(DERIVED — verify_D3D4.m 3/3)*

The `t³,t⁴` coefficients (D.0.3) with the same `C₁(a) ∝ (a−s)^{−4}a^{−1/2}`. Computing the
chain-rule `a''', a''''` and matching `tempa²` lands, with `temp_d = D₂ξC₁/3`,
```
D₃ = (17a₀ + s)·temp_d ,    D₄ = ½·temp_d·a₀·ξ·(221a₀ + 31s)·C₁     (D.3.1) [code:179-182]
```
**EXACTLY** — `verify_D3D4.m` D3/D4, `simplify = 0`. The `(17a₀+s)` and `(221a₀+31s)` arise
purely from the higher `a`-derivatives of `(a−s)^{−4}a^{−1/2}` (the `s` carried by the
`(a−s)^{−4}` factor, the `a₀`-coefficients 17/221 by the combination with `a^{−1/2}`). This
resolves the **cascade-from-D-5** (§11): since D₂ is exact, D₃/D₄ are exact; Q-9/Q-10 (the
`temp_d` grouping) are confirmed and the `17a₀+s`/`221a₀+31s` are now *derived*, not transcribed.

## §T t2cof…t5cof — mean-longitude drag polynomial  *(close D-6; verify_tcofs.m 4/4)*

The drag-perturbed mean motion follows the decaying semi-major axis: since `a(t) = a₀ tempa²`
and `n ∝ a^{−3/2}` (Kepler), `n(t) = n₀ tempa^{−3}`. The drag-induced extra mean longitude is
`ℓ(t) = ℓ₀ + n₀·templ` with
```
templ = ∫₀^t ( n(t')/n₀ − 1 ) dt' = ∫₀^t ( tempa^{−3} − 1 ) dt' .       (T.1)
```
The four t-cofs are the `t²,t³,t⁴,t⁵` Taylor coefficients of (T.1)
(`templ = t2cof t² + t3cof t³ + t⁴(t4cof + t·t5cof)`, `drag_coefficients.h:16`). Expanding
`tempa^{−3} − 1 = 3x + 6x² + 10x³ + 15x⁴` with `x = C₁t + D₂t² + D₃t³ + D₄t⁴` and integrating
term-by-term:
```
t2cof = (3/2)C₁ ,   t3cof = D₂ + 2C₁² ,
t4cof = ¼(3D₃ + C₁(12D₂ + 10C₁²)) ,
t5cof = ⅕(3D₄ + 12C₁D₃ + 6D₂² + 15C₁²(2D₂ + C₁²)) ,                     (T.2)
```
**all four code-matched** (`verify_tcofs.m` 4/4, `simplify = 0` vs `:185-191`). **D-6 CLOSED** —
the legacy §12 stream-of-consciousness was a failure to state (T.1); the structure is a one-line
integral of the mean-motion decay.

---

**Status:** D₂, D₃, D₄, t2cof…t5cof **COMPLETE** — all from one decaying orbit
(`a = a₀ tempa²`, `n = n₀ tempa^{−3}`). `verify_D2.m` 3/3, `verify_D3D4.m` 3/3, `verify_tcofs.m`
4/4 (`simplify(_ − code) = 0` vs `:178-191`). **D-5 CLOSED**, **D-6 CLOSED**, cascade-from-D-5
(§11) resolved.
