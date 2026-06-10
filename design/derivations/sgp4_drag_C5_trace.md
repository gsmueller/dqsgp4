# Phase 4 — SGP4 C₅ Periodic Eccentricity Coefficient

## §C5.0 Scope

Traces C₅ (`src/atmosphere/drag_coefficients.h:174-175`), the amplitude of the **periodic**
(non-secular) eccentricity correction `e(t) = e₀ − tempe`,
`tempe = B*C₄ t + B*C₅(sin M − sin M₀)`. Under **Standard 10**, `verify_C5.m` proves the code
form equals the born-digital SR3 p.11 form; the from-scratch derivation has a **documented open
gap** (the 5/4, see §C5.2).

**Code target** (`drag_coefficients.h:174-175`):
```
C₅ = 2 coef1 a₀ β₀² ( 1 + (11/4)(η² + e₀η) + e₀η³ ) ,   β₀² = 1 − e₀² , coef1 = (q₀−s)⁴ξ⁴(1−η²)⁻⁷ᐟ² .
```
**Born-digital** (SR3 `page_011.md:19`; Rhodes `propagation.py:1462`): identical
(`1 + 2.75(etasq+eeta) + eeta·etasq`, `2.75 = 11/4`). `verify_C5.m` C5.1 confirms `simplify = 0`.

## §C5.1 Mechanism — the m=1 Fourier harmonic of ė

C₄ is the **secular** rate `⟨ė⟩` (multiplies `t`); C₅ is the **m=1 (cos M) Fourier coefficient**
of the oscillatory `ė`. The periodic eccentricity is `δe = ∫(ė − ⟨ė⟩) dM/n`; with `ė` carrying a
`cos M` harmonic of amplitude `a₁`, `∫a₁\cos M\,dM = a₁(\sin M − \sin M₀)`, so `C₅ ∝ a₁`. **C₅ has
no `n₀`** (unlike C₄'s `2 n₀ coef1 a₀ β²`): it multiplies the dimensionless phase `(sin M − sin M₀)`,
not `t`; the `∫…dM` (over `dM`, not `dt`) cancels the `n₀`. `verify_C5.m` C5.3 confirms
`C₄_prefactor = n₀·C₅_prefactor`.

## §C5.2 The 11/4 bracket — clean theory gives 1+4η² (CORROBORATED); code 11/4 is OPERATIONAL

**Trusted facts:** the code value `1 + (11/4)(η²+e₀η) + e₀η³` is code-matched (`verify_C5.m` 3/3 —
public algebra). My from-scratch derivation — the m=1 (cos M) Fourier harmonic of the eccentricity-rate
kernel via the projection `cos M dM = (cos f − 2e cos2f)df` and the kernel `h_C4 = h0 + e·h1` (the
**same** kernel whose m=0 secular average lands C₄ Part-A) — gives the bracket
```
ψ⁷·a₁[h_C4] = (1 + 4η²) + e₀(−2η − 3η³) + O(e²) ,                          (C5.2.1)
```
(`verify_C5_theory.m` 3/3). The η² coefficient is **4**, not the code's **11/4** (gap `−(5/4)η²`); the
e₀-coefficient even **disagrees in sign**. So the code's `11/4` is **not** the cos-M Fourier of the
clean-theory ė kernel.

**This is corroborated, not a personal error.** An **independent** from-scratch investigation (its own
sympy + the born-digital Brouwer-Hori 1961 markdown — **not** OCR) re-derived the **same** projection
and kernel, confirmed both are correct, and reproduced **1+4η²** vs the code's **1+11/4η²**. Its
**smoking gun:** the *same* Gauss kernel gives an m=0 (secular) harmonic that **matches** C₄ Part-A but
an m=1 (periodic) harmonic that does **not** match C₅ — so **C₅ is built from a different object than
C₄-Part-A** (it is *not* the periodic partner of C₄-A under one orbit-average). The correct C₄/C₅
target is the Lane-Cranford angular-momentum rate `ė = (β/(eL))(β L̇ − Ġ)` (`verify_C4C5_target.m`,
the `xnodcf` energy-vs-angular-momentum split), which lands C₄-A but whose periodic realisation for C₅
is the operational reduction step the code uses.

The earlier "`11/4 = 3/2 + 5/4`, the 5/4 is sealed" split is **withdrawn** (it was asserted, not
derived). OCR-sourced "fictitious-longitude integrand" suggestions are **untrusted hypotheses**, not
relied on.

---

**Status:** C₅ **code-matched** (`verify_C5.m` 3/3 — public algebra). **Clean theory gives the bracket
`1 + 4η²` (e⁰) + `e₀(−2η−3η³)`** (`verify_C5_theory.m` 3/3), **independently CORROBORATED**; the code's
**`1 + 11/4η²`** is an **operational SR3/Lane-Cranford form**, NOT reproducible from trusted theory
(C₅ is not the periodic partner of C₄-A — `verify_C4C5_target.m`). Disposition: code-matched for
bit-compatibility; theoretically **UNRESOLVED / operational** (not a derivation error — corroborated).
