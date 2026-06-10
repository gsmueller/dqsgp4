# Phase 9 — SGP4 Drag Correction Terms (xmcof, delmo, sinmo, xnodcf)

Traces the remaining drag correction terms (`src/atmosphere/drag_coefficients.h:194-209`) under
**Standard 10** (`verify_corrections.m` 4/4). Closes audit **D-8** (xmcof ξ⁴) and **D-7** (xnodcf
β-power). Born-digital inputs: SR3 p.12 (`page_012.md:7,13`), Vallado/Rhodes
`propagation.py:1480-1492` (sealed-room dispatch a968ac5e).

> **STATUS (trusted-theory, 2026-06-03).** The **code-matches** below stand (`verify_corrections.m`
> 4/4 — public algebra).
> - **xmcof's `−2/3`: CLEAN-DERIVED (2026-06-03; `verify_xmcof_theory.m` 8/8).** `xmcof` is the
>   **density part** of the same drag-periodic **mean-anomaly correction** `δℓ_D` that `C₃` is the J₃
>   part of: both feed `temp_drag = delomg + delm`, applied `M += temp_drag`, `ω −= temp_drag`
>   (`secular_update.h:97-105`, born-digital) — `xmcof`/`delm` and `C₃`/`delomg` are **partners**. From
>   the trusted Phase-0 Gauss `Ṁ` VOP (Thm 0.3.6) under drag + the Lane density, the mean-anomaly drag
>   correction's antiderivative `∫sin f(1−η cos f)⁻⁴df = −(1/3η)(1−η cos f)⁻³` gives **exactly** the
>   `−2/3` (`= 2 × (−1/3)`: transverse Gauss `2` × density power-4→3 antiderivative `−1/3`), the `1/e₀`
>   (Gauss transverse pole, power 1 — **not** e-vector geometry), the `1/η`, and `coef·B*`. Two
>   **operational** code simplifications (like C₄-B/C₅): the code drops the circular-speed scale `a₀`
>   (`xmcof_theory = a₀ × xmcof_code`), and `(1−η cos f†)⁻³ → (1+η cos M)³` (exact through `O(η)`). See
>   §X.1. The **leading `−2/3` is clean-theory-EXACT** (unlike C₅'s wholly-operational `11/4`).
> - **xnodcf's `β₀²` is RESOLVED (2026-06-02): it is the ANGULAR-MOMENTUM drag normalization.**
>   My earlier "`β₀⁴` three ways" result came from wrongly forcing `L̇/L = Ġ/G = −C₁` (energy-like)
>   on *all* the Delaunay momenta. The correction — the orbit-averaged **angular-momentum** drag rate
>   is `⟨Ġ/G⟩ = −B*⟨ρ|v|⟩ = −C₁·β²` (DERIVED: `⟨ρ|v|⟩ / ((a/μ)⟨ρ|v|³⟩) = β² = 1−e²` to O(e²)),
>   whereas the **energy** rate `L̇/L = −C₁` carries no β — gives `nodecf = ½·C₁·xhdot1·(3+4β²)`,
>   which **carries the β²** and matches the code at leading order. The code's exact `(7/2)β²` applies
>   omeosq uniformly to the gradient sum (the SGP4 `β=1`-in-numerators simplification; code − theory =
>   O(e²)). See §X.2 and `verify_xnodcf_theory.m` (10/10). The angular-momentum-β² insight was an
>   independent agent's *structural hint*, **confirmed here by my own derivation** (never source fact).

## §X.1 xmcof + delmo + sinmo — Lane periodic mean-anomaly correction (D-8)

`delm = xmcof·((1+η cos M)³ − delmo)`, `delmo = (1+η cos M₀)³` (`drag_coefficients.h:202-203`),
`sinmo = sin M₀`. The `(1+η cos M)³` is the operational form of the **cubic** `(1−η cos f†)⁻³` that
the Lane-density **power-4→3 antiderivative** produces (see the clean derivation below); subtracting
`delmo` (the same cubic at epoch) enforces `Δ = 0` at `t₀` (it is the integral's lower limit).

```
xmcof = −(2/3)·coef·B*/(e₀η) ,    coef = (q₀−s)⁴ξ⁴      (drag_coefficients.h:197)
```
**D-8 resolution (ξ⁴, code-matched).** The code's `coef = qoms4·ξ⁴` (`:121`) carries the `ξ⁴` —
`verify_corrections.m` confirms `xmcof_code = (legacy box)·ξ⁴`. The legacy §15.4 *box* dropped `ξ⁴`
(a **doc-only** defect; the code and SR3 p.12 are correct). **D-8 CLOSED.**

**Mechanism + exact `−2/3` — CLEAN-DERIVED (2026-06-03; `verify_xmcof_theory.m` 8/8).** `xmcof` is the
**density part** of the drag-periodic **mean-anomaly correction** `δℓ_D` — the **partner** of `C₃` (the
J₃ part): both feed `temp_drag`, applied `M += temp_drag`, `ω −= temp_drag` (`secular_update.h:97-105`).
The clean derivation, from the **trusted Phase-0 Gauss `Ṁ` VOP** (Thm 0.3.6 / 0.4.2.5) under drag
(0.4.1) with the Lane density `coef/(1−η cos f†)⁴`:

- The drag-periodic mean-anomaly correction is `δℓ_D = ∫_{M₀}^{M} (Ṁ_drag − n)/n · dM′` (the C₅-analog
  construction; `dt = dM/n` to `O(B*)`). The `(Ṁ−n)` rate's **transverse `(β²/(n a e))` term** carries
  a **simple `1/e` pole** (`verify_xmcof_theory.m` X.1-X.2); its residue is
  `2 B* a coef · sin f /(1−η cos f)⁴` (the `2` is the transverse Gauss factor `(1+r/p)|_{e=0}=2`).
- **The KEY step — the Lane-density antiderivative (X.3, exact):**
  `∫ sin f (1−η cos f)⁻⁴ df = (1/η)∫u⁻⁴du = −(1/(3η))(1−η cos f)⁻³`, so
  `δℓ_D = −(2/3)·(B* a coef/(e₀η))·(1−η cos f)⁻³`. This **derives, exactly:** the `−2/3` [`= 2 × (−1/3)`:
  the transverse Gauss `2`, and the `−1/3` from the **density power 4→3** antiderivative], the `1/e₀`
  [Gauss transverse pole, power **1** — per-element, **not** e-vector geometry], the `1/η` [`du = η sin f df`],
  and `coef·B*` [Lane density × drag law]. Matches Phase-0 align note (c): `1/(e₀η) = 1/e [Gauss] × 1/η
  [orbit-average]`. The `−(1+η cos M₀)³` term (`delmo`) is the **lower limit** of this definite integral
  (X.8) ⇒ `Δ=0` at `t₀`.
- **Two OPERATIONAL code simplifications** (NOT clean-theory-exact, like C₄-B/C₅): (i) the code **drops
  the circular-speed scale `a₀`** — clean theory gives `−(2/3)a₀ coef B*/(e₀η)`, i.e. `xmcof_theory =
  a₀ × xmcof_code` (X.7; `a₀ ≈ 1` in ER, and `δℓ_D` is a small correction); (ii) the density factor
  `(1−η cos f†)⁻³ → (1+η cos M)³`, a **reciprocal-binomial** identity exact only **through `O(η)`** (X.6;
  `f†→M` through `O(e)`) — **bounded**: leading error `3η²cos²f ≤ 3η²` (`verify_xmcof_C3_mitigations.m`
  M7), the **same Lane-density η-truncation** the whole C₂-family uses (M8). **Unlike C₅ (whose `11/4` is
  wholly operational), xmcof's LEADING `−2/3` is clean-theory-EXACT** — only the (proven-canonical,
  bounded) `a₀` scalar and the bounded `O(η²)` cubic shape are operational.

> **Independent TARGET confirmation (2026-06-03, cleanroom 2-agent).** Per the trust model, an
> independent agent verified — against the **primary sources** (Lane-Hoots 1979, the parent theory SGP4
> simplifies; SR3 §6, which itself stubs δM/δω and points upstream) — that this derivation solves the
> **right problem**, with a second agent scrubbing the report for OCR before I read it (structure/target
> only; no source equations enshrined). **All structural claims CONFIRMED:** `delm`/`xmcof` + `delomg`/`C3`
> are two summands of **one** mean-anomaly drag-periodic mirrored equal-and-opposite onto perigee; built
> by **integrating the drag-perturbed rate** with the epoch value as the integral's lower limit; `xmcof`
> is the **density** part (carries density/ballistic factors, **no** gravitational-harmonic factor), `C3`
> the **J₃×drag** part; `aycof`/`xlcof` are a **separate, drag-free, pure-J₃** long-period e-vector
> mechanism; and the `1/e` is the **mean-anomaly-equation geometric** factor (the `e₀>1e-4` guard is its
> symptom), **not** an e-vector projection. Two caveats (not refutations): the cleanest target framing is
> the time-integrated **action** perturbation (my rate-integral is the valid equivalent reduced picture —
> a candidate action-route refinement for the operational `a₀`); and C3's clean parent-theory J₃-drag
> object is **critical-inclination-resonance-structured**, the code form being its simplified remnant
> (consistent with the §C3 "partially operational" disposition).

> **Action-route follow-up — `a₀` is INTRINSIC, NOT a rate-route artifact (2026-06-03).** Pursued the
> confirmation's action-route caveat: derived `δℓ` via the Delaunay generating function `δℓ = ∂S/∂L`
> (`δL = −∂S/∂ℓ`) from the drag action rate `L̇ = (1/n)Ė = −B*ρ|v|³/n` (which carries `μ`, **no** explicit
> `a`). Findings: (1) the action route **reproduces the functional shape** (numerical shape-correlation
> `0.9995` vs `delm`), and `∂η/∂L ≈ −∂η/∂G` at small `e` **explains the `M += temp`, `ω −= temp`
> equal-and-opposite structure**; (2) but `δℓ = ∂S/∂L` **also carries `a₀`** — `∂/∂L` brings down
> `L²/μ = a₀` via the dominant `∂η/∂L ∝ η/(Le²)` channel. So **both routes give `δℓ ∝ a₀`**. Numerically
> (`verify_xmcof_theory.m` X.9) the **full-Gauss** rate-route `δℓ/delm_code → a₀` as `η→0`
> (`1.3055` vs `a₀=1.30` at `η=0.036`). xmcof's leading `−2/3` stays clean-exact.

> **Caveat-1 MITIGATION — the `a₀` is CANONICAL, not a route-artifact (`verify_xmcof_C3_mitigations.m`
> M1-M3, 2026-06-03).** The confirmation worried the rate-integral "could create an artificial puzzle"
> vs a "cleaner action route." Settled by **canonical perturbation theory**: a generating function `W`
> (with `δℓ=∂W/∂L` — the action route) **exists only if the action 1-form `ω = δL dℓ + δG dg` is closed**,
> i.e. the perturbation is Hamiltonian. **Drag is dissipative:** `L̇,Ġ` are independent of the perigee `g`,
> yet `Ġ` varies over `ℓ`, so `∂(δL)/∂g = 0 ≠ ∂(δG)/∂ℓ` — **`ω` is not closed ⇒ no `W` exists** (M1).
> Therefore the rate route `δℓ=∫(ℓ̇−⟨ℓ̇⟩)dt` is the **unique canonical technique**; the prior follow-up's
> non-clean `~40×` "action route" was a **mis-applied (nonexistent) `W`**, not a real ambiguity. The `a₀`
> is the `|v|/n = a` circular-speed scale **the sibling C5 retains** (`cc5 = 2 coef1·a₀·β²·…`, M2) and that
> xmcof drops — a single **localized, bounded** simplification (M3, abs. `M`-error `~(a₀−1)|δℓ| ~ 10⁻⁷
> rad`), **not** a derivation error and **not** resolvable by a cleaner route (there is none). Earlier
> "action route gives `δℓ ∝ a₀`" findings stand (the rate route is the action route, done correctly).

## §X.2 xnodcf — drag×J₂ nodal coupling (D-7 RESOLVED 2026-06-02: β² = angular-momentum drag norm)

```
xnodcf = −(21/2)·n₀·k₂·cos i₀·C₁/(a₀²β₀²) = (7/2)·β₀²·xhdot1·C₁  (drag_coefficients.h:208-209)
  xhdot1 = Ω̇_{J₂} = −3 k₂ n₀ cos i₀/(a₀²β₀⁴)  (the J₂ secular nodal rate; β₀⁴ = 1/p² normalisation)
```
The algebraic code-match `(7/2)β₀²·xhdot1·C₁ = −(21/2)n₀k₂cos i₀ C₁/(a₀²β₀²)` is public algebra
(`verify_corrections.m` 4/4, `simplify=0`). The **trusted-theory derivation** (`verify_xnodcf_theory.m`
10/10) is:

**Step 1 — gradient structure.** In Delaunay momenta `Ω̇_{J₂} ∝ H·L⁻³·G⁻⁵` (`L=√(μa)`, `G=Lβ`,
`H=G cos i`; verified `= xhdot1`). Drag is **in-plane** (velocity-aligned, no out-of-plane force) so
it does not change `i` directly ⇒ `Ḣ/H = Ġ/G`, hence
```
Ω̈/Ω̇ = Ḣ/H − 3 L̇/L − 5 Ġ/G = −3 L̇/L − 4 Ġ/G .                         (X.2.1)
```

**Step 2 — the two drag rates differ (THE KEY).** Drag couples to **energy** and **angular momentum**
with *different* normalisations:
```
energy:    L̇/L = ½ȧ/a = −C₁                          (defines C₁; NO β-power)
ang. mom.: Ġ/G = −B*ρ|v|  (torque ḣ = r×F = −B*ρ|v| h) ⇒ ⟨Ġ/G⟩ = −C₁·β²   (CARRIES omeosq)
```
The `β²` is **derived**, not asserted: with `C₁ = (a/μ)B*⟨ρ|v|³⟩` (the energy rate) and `⟨Ġ/G⟩ =
−B*⟨ρ|v|⟩`, the vis-viva ratio (vis-viva: `(a/μ)|v|² = 2a/r−1`)
```
⟨Ġ/G⟩ / (−C₁) = ⟨ρ|v|⟩ / ((a/μ)⟨ρ|v|³⟩) = (1−¼e²)/(1+¾e²) = 1−e² = β²  (to O(e²)).  (X.2.2)
```
(`verify_xnodcf_theory.m` check 3c, `simplify=0`.) This is the `omeosq` the code carries. *(Honesty:
(X.2.2) is the **near-circular / constant-density leading** result — it isolates the angular-momentum
β² **origin** cleanly; the Lane-density `(1−η cos f†)⁻⁴` adds O(η²) corrections to the exact ratio.
The qualitative point — angular-momentum drag carries an omeosq the energy drag does not — is robust;
the code's exact `(1−e²)` and its uniform application are the SGP4 operational form.)*

**Step 3 — assembly.** Substituting `L̇/L=−C₁`, `Ġ/G=−C₁β²` into (X.2.1):
```
nodecf = ½Ω̈ = ½·C₁·xhdot1·(3 + 4β²) .                                   (X.2.3)
```
This **carries the β²** (refuting my earlier "β₀⁴, no omeosq" hypothesis, which had wrongly used
`Ġ/G = −C₁`). At `e→0` it equals the code's `(7/2)C₁·xhdot1`.

**Step 4 — the code is the uniform-omeosq simplification.** The code's exact `(7/2)β²·xhdot1·C₁`
applies the omeosq factor **uniformly** to the whole gradient sum `7` (the SGP4 "β=1 in numerators"
rule), rather than only to the two angular-momentum gradient terms. Hence
```
nodecf_code − nodecf_theory = −(3/2)·C₁·xhdot1·e² = O(e²)                 (X.2.4)
```
— they agree at leading order; the difference is the SGP4 simplification, **not** an unexplained
discrepancy. (`verify_xnodcf_theory.m` check 5b.)

**Provenance note.** The structural hint "the node's secular drag drift is an *angular-momentum*
effect, so it carries β² that the energy drag (C₁/C₂) does not" came from an independent verification
agent reading the Lane-Hoots method (structure only); it is **confirmed above by my own derivation**
(X.2.2) and never taken as source fact (per the trust model, OCR is never fact).

**D-7 resolution.** The earlier "`β₀⁻³` into `C₁`" guess is still **wrong** (`C₁=B*C₂` carries only
`ψ⁻⁷`, `η=a₀e₀ξ≠e₀` — Phase 0 Rem 0.6.5, `β≠ψ`). The `β₀²` is the **angular-momentum drag
normalisation** (X.2.2), and the code's exact form is its uniform-omeosq SGP4 simplification (X.2.4).

---

**Status:** xmcof, delmo, sinmo, xnodcf **code-matched** (`verify_corrections.m` 4/4,
`simplify(_ − code) = 0` vs `:197,202-205,208-209`). **D-8 CLOSED** (code/SR3 carry ξ⁴; legacy
§15.4 box dropped it — doc fix applied). **D-7 RESOLVED (2026-06-02, theory):** the `β₀²` is the
**angular-momentum drag normalisation** `⟨Ġ/G⟩ = −C₁β²` (derived, X.2.2), giving
`nodecf = ½C₁·xhdot1·(3+4β²)` which carries the β²; the code's exact `(7/2)β²` is the uniform-omeosq
SGP4 simplification (code − theory = O(e²), X.2.4). `verify_xnodcf_theory.m` 10/10. **xmcof's `−2/3`
is now CLEAN-DERIVED (2026-06-03):** the Gauss `Ṁ` transverse `1/e` pole + the Lane density power-4→3
antiderivative give `−2/3 = 2 × (−1/3)`, `1/(e₀η)`, `coef·B*` exactly (`verify_xmcof_theory.m` 8/8); the
`a₀` scalar and `(1+η cos M)³` cubic shape are operational (§X.1 banner). **D-8 closed; xmcof complete.**
