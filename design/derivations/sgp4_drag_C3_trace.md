# Phase 5 — Symbolic Trace of the SGP4 C₃ (J₃ Eccentricity–Drag Coupling)

## §C3.0 Scope, code target, born-digital inputs

Traces the SGP4 C₃ coefficient (`src/atmosphere/drag_coefficients.h:154-159`) — the
amplitude of the **J₃-driven long-period coupling into the drag-secular argument-of-perigee
and mean-anomaly corrections** — under **Standard 10** (`simplify(C₃_derived − C₃_code) == 0`).
**Resolves audit finding D-4.**

**Code target** (`drag_coefficients.h:156`, guard `e₀ > 1e-4` else `C₃ = 0`):
```
C₃ = coef · ξ · A₃₀ · n₀ · sin i₀ / (k₂ · e₀) ,
  coef = (q₀−s)⁴ξ⁴ ,  ξ = 1/(a₀−s) ,  A₃₀ = A(3,0) = −J₃ a_E³ ,  k₂ = J₂ a_E²/2 .
```

**Convention (load-bearing, model_selector.h:178).** `A(3,0)` returns `−Jₙ` (no ½):
`A₃₀ = −J₃ a_E³` (= `−J₃` in SGP4 normalized units `a_E = 1`). Hence the **Brouwer ratio**
```
A₃₀/k₂ = −J₃ a_E³ / (J₂ a_E²/2) = −2 (J₃/J₂) a_E = −2·j3oj2·a_E      (verify C3.2)
```
(`verify_C3.m` C3.2, `simplify = 0`). The **deprecated** `021_long_period_corrections.md`
uses Brouwer's `A₃₀ = −J₃a_E³/2` (½ absorbed) — **do not** mix conventions.

**Born-digital inputs** (sealed-room dispatch ac84fb143; I verify the assembly):
- **J₃ disturbing function** `R₃ = (μJ₃a_E³/2r⁴)(5sin³i sin³u − 3 sin i sin u)`, `u = ω+f`
  (SR3 transcript Eq 4.4, `hoots_roehrich_1980_math_derivation.md:400`; `sin(lat) = sin i sin u`).
  Every term carries `sin(n·u)` ⇒ **J₃ has no secular term**; its leading effect is a
  **long-period oscillation of the eccentricity vector** `(e sin ω, e cos ω)` (Kozai 1959).
- **SGP4 long-period structure** (the aycof/xlcof of `near_space.h:280,286`, = Vallado/Rhodes
  `propagation.py:1485,1488`).
- **Vallado/Rhodes equivalent** `cc3 = −2·coef·tsi·j3oj2·n₀·sinio/ecco` (`propagation.py:1454`,
  `tsi ≡ ξ`), identical to the code form via `A₃₀/k₂ = −2 j3oj2` (`verify_C3.m` C3.4).

> **STATUS (2026-06-02 — mechanism CORRECTED; §C3.1-§C3.3 REFUTED).** The **code-match** in §C3.4 /
> `verify_C3.m` 6/6 stands (public algebra). The **physical mechanism** is now corrected: my
> §C3.1-§C3.3 "Gauss-VE factor-2 / eccentricity-vector projection" account is **REFUTED** — it
> reproduces the code form only **coincidentally**.
> - **Code-confirmed target** (`src/orbit/secular_update.h:96-105`, born-digital): `omgcof = B*·C₃·cos ω₀`
>   feeds `delomg`, which is **added to the mean anomaly and subtracted from the argument of perigee**
>   (`M += delomg+delm`, `ω −= delomg+delm`). So **C₃ is the J₃×drag amplitude of the drag-periodic
>   MEAN-ANOMALY correction** (with a sign-flipped copy to perigee) — the **same machinery** as
>   `xmcof`'s `delm` (the density part). It is **NOT** a standalone eccentricity-vector term (those are
>   the separate, **C₃-free** `aycof`/`xlcof` long-period J₃ terms — `near_space.h:280-286`).
> - The `1/e₀` is a **plain per-element normalisation** (the `e₀·δℓ → δℓ` factor; **power exactly 1**),
>   **not** an eccentricity-vector projection. The §C3.1/§C3.3 "factor-2 cancellation" is **coincidental**.
> - **CLEAN factor derivation (2026-06-03; `verify_C3.m` 8/8, partner of `verify_xmcof_theory.m`).** Now
>   that `xmcof`/`delm` is clean-derived (the Gauss `Ṁ` transverse `1/e` pole + Lane-density power-4→3
>   antiderivative — see `sgp4_drag_corrections_trace.md` §X.1), C₃'s factors follow as the J₃ partner:
>   - **`1/e₀` (power 1) = the Gauss `Ṁ` transverse SIMPLE pole** (`verify_C3.m` C3.7: residue
>     `limit_{e→0}(e·(Ṁ−n)/n) = 2 a B* sin f`, finite & nonzero ⇒ power exactly 1) — **the same `1/e`
>     that lands `xmcof`**, definitively **not** an e-vector projection.
>   - **`A₃₀/k₂ = −2 j3oj2` = pure convention** (C3.2); **`n₀` = secular rate × t** (`delomg = B*C₃cos ω₀·t`
>     is linear in `t`, so its coefficient is a rate — unlike `xmcof`'s dimensionless phase).
>   - **`(A₃₀/k₂)·sin i₀` = the J₃ long-period e-vector amplitude** `= 4·aycof` (`verify_C3.m` C3.8;
>     `aycof = ¼(A₃₀/CK2)sin i₀`, `near_space.h:286`): C₃ inherits the **J₃ long-period amplitude**, not
>     an independent e-vector projection of the drag.
>   - **`cos ω` is CLEAN — resonance-free by EXACT cancellation** (caveat-2 mitigation,
>     `verify_xmcof_C3_mitigations.m` M4-M6): `C₃·cos ω₀·t` is the **secular (linear-in-`t`) part** of the
>     J₃-modulated drag mean-anomaly drift along the apsidal drift `ω(t)=ω₀+ω̇t`. The secular part is the
>     **linear-`t` coefficient** of the accumulated rate `[sin(ω₀+ω̇t)−sin ω₀]/ω̇`, which is `cos ω₀`
>     **independent of the apsidal frequency `ω̇ ∝ (5cos²i−1)`** — so the **critical-inclination resonance
>     cancels EXACTLY** out of `C₃` (M4/M4b). The resonant `1/ω̇` piece lives only in the **oscillatory**
>     long-period block, handled **separately** by `aycof`/`xlcof` (themselves resonance-free — only the
>     `(1+cos i)` *polar* singularity, not critical inclination). So `C₃`'s `cos ω` structure is the
>     **EXACT resonance-free secular remnant**, NOT a lossy truncation.
>   - **`ξ⁵ = coef·ξ`** = the drag density at the `I^{(0,5)}` Hansen level (same family as C₂ Part B's
>     `(8+24η²+3η⁴)`), **not** a `∂ρ/∂r` gradient. **The overall coefficient `1` is now DERIVED**
>     (`verify_xmcof_C3_mitigations.m` M8): it is the **leading (η→0) Hansen density** `I^{(0,5)}(0) = 1`
>     (Phase 1: `I^{(0,m)}(0)=1` ∀m). The dropped η-corrections (the `(8+24η²+3η⁴)`-type polynomial C₂
>     Part B keeps) are the **same Lane-density truncation** as `xmcof`'s cube and C₂'s dropped AFGP4
>     terms — **one unified, bounded residual** (M7: xmcof cube is `O(η²)`, error `3η²cos²f`), not a
>     C₃-specific gap.
>
> §C3.1-§C3.3 below are **retained but REFUTED** (kept only to document the superseded attempt; the
> `verify_C3.m` C3.1/C3.3 checks are true **identities**, now **labelled `[IDENTITY, not physics]`**).

## §C3.1 The Gauss-VE factor 2 (Phase 0 Theorem 0.3.3)

The eccentricity rate `ė = (β/na)[R sin f + T(cos f + cos E)]` (Theorem 0.3.3) carries, in its
transverse term, a **leading factor of 2** — the `2hT r̂` of the eccentricity-vector rate
`de_vec/dt` (0.3.3.10), surfaced via `2h cos f + rṙ sin f = h(cos f + cos E)` (0.3.3.16). It is
visible as
```
(cos f + cos E)|_{e→0} = 2 cos f ,    cos E = (e+cos f)/(1+e cos f)      (verify C3.1)
```
(`verify_C3.m` C3.1, `simplify = 0`). **This is the 2 that resolves D-4** (it cancels the
Brouwer `A₃₀/(2k₂)` prefactor's 2).

## §C3.2 Factor-by-factor assembly

The J₃ radial/eccentricity-vector perturbation `δr_{J₃} ∝ A₃₀ sin i sin u` modulates the Lane
drag and couples into the eccentricity-vector rate. Each factor of C₃ (dispatch ac84fb143):

| Factor | Origin |
|---|---|
| `coef = (q₀−s)⁴ξ⁴` | Lane τ=4 base density `ρ = ρ₀(q₀−s)⁴/(r−s)⁴`, orbit-averaged (= C₂ Part A). |
| **extra `ξ`** | **[CORRECTED 2026-06-03 — NOT a `∂ρ/∂r` gradient]** the drag density enters at the `I^{(0,5)}` **Hansen level** (`ξ⁵=ξ⁴·ξ`, same family as C₂ Part B's `(8+24η²+3η⁴)`); the extra `ξ` is the higher density power in the J₃-coupled mean-anomaly integral, not a density-gradient linearisation. See banner. |
| `A₃₀` | amplitude of the J₃ radial perturbation `δr_{J₃} ∝ A₃₀ sin i sin u`. |
| `n₀` | the drag **rate** is a frequency × dimensionless density factor (every C-coef carries one n₀). |
| `sin i₀` | the single `sin i·sin u` harmonic of R₃ surviving the average (the `5sin³i sin³u` is higher-harmonic, dropped). |
| `1/k₂` | packaged with `A₃₀` as the Brouwer ratio `A₃₀/k₂` (J₃-relative-to-J₂ normalization). |
| `1/e₀` | the coupling lives in the eccentricity **vector** `(e sin ω, e cos ω)`; projecting the vector rate onto scalar `(e, ω)` introduces `1/e` (singular as `e→0` ⇒ the `e₀>1e-4` guard). |

The relevant average is the **single `sin u` harmonic of `δr_{J₃}` against the Gauss
eccentricity-rate kernel** (producing the `cos ω` that feeds `omgcof` and a `sin ω`), **not**
any `⟨sin²u⟩` — the resonant `sin u/cos u` projection has coefficient 1, contributing no ½.

## §C3.3 D-4 resolution — exactly how `/(k₂·e₀)` arises

The legacy proof (`sgp4_near_earth_drag_theoretical_basis.md:282-285`) wrote the prefactor
`A₃₀/(2k₂)` and then claimed "absorbing the factor of 2 from `⟨sin²u⟩ = ½`." **Wrong twice:**
(i) `⟨sin²u⟩ = ½` *multiplies* by ½, sending the denominator to `/(4k₂e₀)` — off by **4×**
(`verify_C3.m` C3.5: `4·C₃_legacy = C₃_code`); (ii) there is **no** `sin²u` in the chain.

**Correct chain.** The `A₃₀/(2k₂)` Brouwer prefactor's 2 is cancelled by the **Gauss-VE factor 2**
of the eccentricity-vector rate `de/dt` (§C3.1):
```
C₃ = coef · ξ · [A₃₀/(2k₂)] · [2 (Gauss de/dt)] · n₀ · sin i₀ / e₀
   = coef · ξ ·  A₃₀ · n₀ · sin i₀ / (k₂ · e₀)                          (verify C3.3)
```
(`verify_C3.m` C3.3, `simplify(C₃_brouwer_gauss − C₃_code) = 0`). The denominator's 2 is the
`A₃₀/(2k₂)` 2 killed by the Gauss 2 — **not** a trig identity, **not** `⟨sin²u⟩`.

## §C3.4 Code-match witness

```
simplify(C₃_derived − C₃_code) = 0 ,     C₃_code = coef·ξ·A₃₀·n₀·sin i₀/(k₂·e₀)  [drag_coefficients.h:156]
simplify(C₃_code − (−2 coef ξ j3oj2 a_E n₀ sin i₀/e₀)) = 0                         [Vallado cc3]
omgcof = B*·C₃·cos ω₀                                                              [drag_coefficients.h:194]
```
**Verified** — `verify_C3.m` 8/8 PASS. The boxed form / code is **correct**; only the legacy
*proof prose* (§7 lines 282-285) was wrong (a Phase 9 legacy-doc cleanup item — replace the
`⟨sin²u⟩` prose with the §C3.7 Gauss-`Ṁ` `1/e` account).

---

**Status:** C₃ **code-matched** (`verify_C3.m` 8/8 vs `drag_coefficients.h:156`; omgcof vs `:194`) —
**D-4 CLOSED** (the boxed form / code is correct; the legacy `⟨sin²u⟩=½` proof step was wrong).
**Mechanism CORRECTED (2026-06-02) + factors CLEAN-DERIVED (2026-06-03):** C₃ is the **J₃×drag amplitude
of the drag-periodic mean-anomaly correction** `δℓ_D` (sign-flipped to perigee; `secular_update.h:96-105`),
the **partner of `xmcof`'s density term** — **NOT** an eccentricity-vector amplitude (those are the
C₃-free `aycof`/`xlcof`). My §C3.1-§C3.3 "Gauss-VE factor-2 / e-vector `1/e` projection" account is
**REFUTED** (the `verify_C3.m` C3.1/C3.3 checks are now relabelled `[IDENTITY, not physics]`). **Clean
factors (banner, `verify_C3.m` C3.7-C3.8):** `1/e₀` = the Gauss `Ṁ` transverse **simple pole** (power 1,
the same `1/e` as `xmcof` — `verify_xmcof_theory.m`); `A₃₀/k₂ = −2 j3oj2` convention; `n₀` = secular
rate×t; `(A₃₀/k₂)sin i₀ = 4·aycof` (J₃ long-period e-vector amplitude). `ξ⁵` (drag density at `I^{(0,5)}`
Hansen level, ∥ C₂ Part B — **not** a `∂ρ/∂r` gradient). **`cos ω` is now CLEAN** — resonance-free by
**exact cancellation** (caveat-2 mitigation, `verify_xmcof_C3_mitigations.m` M4-M6): `C₃·cos ω₀·t` is the
**secular linear-`t` part** of the J₃-drag drift along the apsidal drift, whose linear coefficient is
`cos ω₀` **independent of `ω̇ ∝ (5cos²i−1)`**, so the critical-inclination resonance cancels exactly; the
resonant `1/ω̇` oscillation is the **separate** `aycof`/`xlcof` block (also resonance-free). The `ξ⁵`
**coefficient `1` is DERIVED** = the leading Hansen density `I^{(0,5)}(0)=1`; the only residual is the
**bounded Lane-density η-truncation** shared by the whole C₂-family (`verify_xmcof_C3_mitigations.m` M7-M8). **Independent TARGET confirmation (2026-06-03, cleanroom):** an independent agent verified against
the primary sources (Lane-Hoots 1979 / SR3) that C₃ is the **J₃×drag part of the one mean-anomaly
drag-periodic** (partner of `xmcof`), that `aycof`/`xlcof` are the **separate, drag-free, pure-J₃**
e-vector mechanism, and that the `1/e₀` is the **mean-anomaly-equation geometric** factor — all CONFIRMED
(a second agent scrubbed the report for OCR first; no source equations enshrined). Its "critical-
inclination resonance" caveat is **MITIGATED**: the resonance cancels exactly out of the secular `C₃` (the
code form is the clean remnant, **not** a lossy truncation). Legacy §7 prose fix remains a doc-cleanup item.
