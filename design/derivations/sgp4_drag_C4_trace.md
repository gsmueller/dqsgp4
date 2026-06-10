# Phase 3 — Symbolic Trace of the SGP4 C₄ Eccentricity-Decay Coefficient

## §C4.0 Scope and inputs

Traces the SGP4 C₄ coefficient (`src/atmosphere/drag_coefficients.h:162-171`) — the
secular eccentricity-decay rate `e(t) = e₀ − tempe`, `tempe = B*C₄ t + …` — from first
principles under **Standard 10** (`simplify(C₄_derived − C₄_code) == 0`).

**Code target** (`drag_coefficients.h:162-171`):
```
C₄ = 2 n₀ coef1 a₀ β² · [ KEPLERIAN + J2-COUPLING ] ,     β² = betao2 = 1 − e₀² ,
  KEPLERIAN   = η(2 + ½η²) + e(½ + 2η²)
  J2-COUPLING = −(2 k₂ ξ)/(a₀ ψ²) · [ −3 con41 (1 − 2eη + η²(³⁄₂ − ½eη))
                                      + ¾ sin²i₀ (2η² − eη(1+η²)) cos 2ω₀ ]
```
`con41 = 3cos²i₀ − 1`, `ψ² = |1−η²|`, `k₂ = J₂/2`, `coef1 = (q₀−s)⁴ξ⁴ψ⁻⁷`.

**Born-digital cross-check (SR3 p.11, verbatim):** the SGP4 C₄ is, with `θ = cos i₀`,
```
C₄ = 2 n₀ (q₀−s)⁴ξ⁴ a₀ β₀² ψ⁻⁷ ·( [2η(1+e₀η) + ½e₀ + ½η³]
     − 2k₂ξ/(a₀ψ²)·[ 3(1−3θ²)(1 + ³⁄₂η² − 2e₀η − ½e₀η³)
                     + ¾(1−θ²)(2η² − e₀η − e₀η³) cos 2ω₀ ] )
```
(`sgp4_references/.../SR3_pages/page_011.md:16-17`, sealed-room dispatch a1fd272ae).
This is **algebraically identical** to the code: `[2η(1+e₀η)+½e₀+½η³] = η(2+½η²)+e₀(½+2η²)`;
`3(1−3θ²) = −3 con41`; `(1−θ²) = sin²i₀`; `(2η²−e₀η−e₀η³) = 2η²−e₀η(1+η²)`.

**Inputs consumed:**
- **Phase 0-rev1 Theorem 0.3.3** — instantaneous `ė = (β/na)[R sin f + T(cos f + cos E)]`,
  `cos E = (e+cos f)/(1+e cos f)`. (`sgp4_drag_phase0_foundations/03_gauss_ve_part1_a_e.md`)
- **Phase 0-rev1 Theorem 0.4.x / 0.5.3** — drag specialization `R,T ∝ −B*ρ|v|(ṙ, rḟ)`,
  `|v| = (na/β)(1+e²+2e cos f)^{1/2}`; orbit-average `dM = β³/(1+e cos f)² df`.
- **Phase 1 Theorem 1.3.2** — Lane integrals `I^{(0,m)}` (the `p = 0` family; Phase 2.C
  proved C₄ needs no `p ≥ 1`).
- **J₂ short-period perturbations** (born-digital, dispatch a1fd272ae; 4 concordant sources
  SR3 p.14 / dnwrnr `SGP4.cc:618` / Rhodes `propagation.py:1917` / Lara 2021):
  ```
  δr_secular = − (3/2)(k₂/p²) β (3θ²−1) · r        [r-multiplier, the C₂ Part-B input]
  Δr_shortper = (k₂/2p_L) (1−θ²) cos 2u ,  u = ω+f  [coefficient of sin²i cos2u = k₂/(2p_L)]
  Δ(rḟ)       = (k₂ n/p_L)[ (1−θ²) cos 2u − (3/2)(1−3θ²) ]   [NOTE the (3θ²−1) SECULAR part]
  Δṙ          = −(k₂ n/p_L)(1−θ²) sin 2u
  ```
  `p_L = a(1−e²) = aβ²`. The geopotential origin: `V_J2 ∝ ½(3sin²(lat)−1)`,
  `sin(lat) = sin i sin u`, so `sin²i sin²u = ½sin²i(1−cos2u)` splits into the `(3θ²−1)`
  secular constant and the `sin²i cos2u` short-period harmonic.

## §C4.A Part A — Keplerian eccentricity decay  *(DERIVED — verify_C4.m C4.1-C4.4)*

**ė-bracket collapse.** Substituting the drag `R = −B*ρ|v|ṙ`, `T = −B*ρ|v|(rḟ)` and the
Kepler primitives (Phase 0) into `ė = (β/na)[R sin f + T(cos f + cos E)]`, the entire
`f`-bracket reduces:
```
e sin²f + (1 + e cos f)(cos f + cos E) = 2(e + cos f)            (C4.A.1)  [verify C4.1]
```
(`sin²f + cos²f = 1` collapses the eccentric-anomaly terms). Hence, with the Lane density
`ρ → ρ₀(q₀−s)⁴ξ⁴(1−η cos f†)^{−4}` (B* carries ρ₀) and the Jacobian, the secular rate is
```
⟨ė⟩ = −2 B* (q₀−s)⁴ξ⁴ (n a β²) · ⟨ h_C4 ⟩ ,
  h_C4 = (1−η cos f†)^{−4} (1+e²+2e cos f)^{1/2} (e + cos f)/(1 + e cos f)² .   (C4.A.2)
```
Since `tempe`'s secular part is `B*C₄ t`, `ė_sec = −B*C₄`, so `C₄ = 2 coef1 (naβ²)·ψ⁷⟨h_C4⟩`
with `n=n₀, a=a₀` → the code prefactor `2 n₀ coef1 a₀ β²`.

**`p = 0` reduction (Phase 2.C scoping).** Expanding `h_C4` in `e` at fixed `η` (exact
density, §A.3 method):
```
h0 = cos f · (1−η cos f)^{−4} ,    h1 = (1 + 3η cos f) sin²f · (1−η cos f)^{−5}  (C4.A.3)
```
— purely `(1−η cos f)` powers. Reducing via `cos f = (1−(1−η cos f))/η` and the `sin²f`
identity to `I^{(0,m)}`:
```
ψ⁷⟨h0⟩ = η(2 + ½η²) ,    ψ⁷⟨h1⟩ = ½ + 2η²                          (C4.A.4)  [verify C4.3]
```
so the **Keplerian bracket** `ψ⁷(⟨h0⟩ + e⟨h1⟩) = η(2+½η²) + e(½+2η²)` = `drag_coefficients.h:163`
**exactly** (`verify_C4.m` C4.4, `simplify = 0`). Part A **lands**.

## §C4.B Part B — J₂ coupling  *(2026-06-02: correct target identified; clean theory → order-η CORROBORATED; code order-1 is OPERATIONAL)*

**The correct target (2026-06-02).** SGP4's eccentricity coefficients are NOT a naïve osculating Gauss
`⟨ė⟩`; they are the **Lane-Cranford translation of the Brouwer-Hori angular-momentum drag rate** into
Keplerian form. The mean-element eccentricity rate is the **kinematic** combination
```
ė = (β/(e L))·(β L̇ − Ġ) ,    e = √(1 − (G/L)²) ,  G = L β .              (C4.B.1)
```
(`verify_C4C5_target.m` check 1, `simplify=0`.) So `ė` is built from the **energy** rate `L̇` (→ `a`,
`C₁/C₂`; carries **no** β) and the **angular-momentum** rate `Ġ` (→ `e`, node; carries the β/η
normalisation) — **the same energy-vs-angular-momentum distinction that RESOLVED `xnodcf`**
(`verify_xnodcf_theory.m`; `⟨Ġ/G⟩ = −C₁β²`). This is why C₄/C₅ carry the `β²` prefactor that C₁/C₂
lack, and why the `(e+cos f)` Gauss kernel of §C4.A is the **leading realisation** of `β L̇ − Ġ`
(Part A lands exactly — positive control).

**Why the isotropic Part-B term is order-η (the cancellation).** BH61's `dG″/dt` *does* have an
order-1 isotropic J₂-secular bracket — but that **same** constant bracket is shared by `L̇` and `Ġ`,
so it **CANCELS** in the combination `β L̇ − Ġ` at `e→0` (`verify_C4C5_target.m` check 2:
`(β L̇ − Ġ)|_{e=0} = 0` for a shared bracket). What survives is **order-η**. Confirmed three further
ways (`verify_C4B_theory.m`): the C2-style density gradient gives `⟨cos f(1−η cos f)^{−5}⟩ =
(I⁰⁵−I⁰⁴)/η = O(η)` (numerator `15η³+20η`, vanishes at η=0); the transverse-velocity J₂ beat is
order-η; the osculating Gauss e-vector rate `⟨dė_vec/dt⟩ = 0` at `e=0` (the J₂ perturbations are
frequency-2 in `u`, orthogonal to the frequency-1 e-vector geometry). **The legacy §8(ii)
`⟨cos2(ω+f)⟩ = cos2ω⟨cos2f⟩` sketch is wrong** (`⟨cos2f⟩=0`).

**Result (trusted, CORROBORATED).** Part B is **code-matched to operational SGP4** (`verify_C4.m`
C4.B/C4.full, `simplify=0`, 8/8 — public algebra). The clean-theory isotropic term is **order-η**
(vanishes at `e=0`); the **code is order-1** (`−3 con41·1` at `e=η=0`). An **independent** from-scratch
investigation (its own sympy + the born-digital BH61 markdown — **not** OCR) **corroborated** this:
it could find **no** trusted target/coupling that yields the code's order-1 isotropic constant. So the
disagreement is **not** an error in my derivation — it is genuine.

**Disposition.** The clean-theory C₄ Part-B isotropic term is **order-η** (corroborated). The code's
**order-1** isotropic constant is an **operational SR3/Lane-Cranford form** — code-matched for
bit-compatibility, but **NOT reproducible from trusted (theoretical or born-digital) drag theory** →
**theoretically UNRESOLVED / operational**. (The OCR-sourced "operational erratum / reduction-step"
provenance is a *plausible hypothesis* for *why* the code differs, but is **not asserted as fact** —
OCR is never fact.)

---

**Status:** Part A (Keplerian) **DERIVED** (`verify_C4.m` C4.1-C4.4). Part B **code-matched** (8/8).
**Correct target identified** — the Lane-Cranford angular-momentum rate `ė = (β/(eL))(β L̇ − Ġ)`
(`verify_C4C5_target.m`), unifying C₄/C₅ with `xnodcf`. Clean theory gives the isotropic term
**order-η** (the order-1 J₂ bracket cancels in `β L̇ − Ġ`); **independently CORROBORATED** (not OCR).
The code's **order-1** constant is an **operational form**, theoretically UNRESOLVED (not a derivation
error on my part). Code-match for bit-compatibility unaffected. `verify_C4B_theory.m` 5/5,
`verify_C4C5_target.m` 2/2.
