# SGP4 Drag Re-Derivation — Remaining-Phases Plan (durable)

> **Durable in-repo roadmap** for the remaining SGP4 near-earth drag coefficient
> symbolic re-derivations. Mirrors the Phase 2.A template. Survives plan-mode loss
> (the prior session's ephemeral plan file was lost; this is the canonical copy).
> Approved 2026-06-01 on branch `session/2026-04-23`. Update the **Status ledger**
> below as phases land.

## Outcome (2026-06-01b) — plan executed

**Every SGP4 drag coefficient is now code-matched** (`simplify(derived − code) = 0`): C₁, C₂, C₃,
C₄, C₅, D₂, D₃, D₄, t2cof…t5cof, omgcof, xmcof, delmo, sinmo, xnodcf. **89 verifier checks PASS**
across 10 `verify_*.m`. Regression unchanged (**23/33, 506/623**). **Findings closed:** D-1, D-3,
D-4, D-5, D-6, D-7, D-8, D-9, D-11, Q-4, cascade-from-D-5, D-3-cascade-§8.

**STATUS (2026-06-02, after independent verification — trusted-theory).** The code values are all
code-matched (verifiers PASS — public algebra; the implementation is unaffected). The **mechanisms**
I'd claimed are NOT trusted — several rested on a sub-agent reading a **scanned** source (OCR,
inherently untrusted), which is **never fact**. Trusted-theory position:
1. **C₄ Part B — correct target identified; clean theory → order-η CORROBORATED; code order-1 is
   OPERATIONAL.** The correct C₄/C₅ target is the **Lane-Cranford angular-momentum rate**
   `ė = (β/(eL))(β L̇ − Ġ)` (`verify_C4C5_target.m`; the `xnodcf` energy-vs-angular-momentum split —
   `L̇`→a/no-β, `Ġ`→e/carries-β). C₄ Part-A lands; the isotropic Part-B term is **order-η** because the
   order-1 J₂ bracket is shared by `L̇`,`Ġ` and **cancels** in `β L̇ − Ġ`. **Independently CORROBORATED**
   (sympy + born-digital BH61, not OCR). The code's **order-1** constant is an **operational SR3/
   Lane-Cranford form**, NOT clean-theory-derivable → code-matched, theoretically UNRESOLVED (not my
   error). `verify_C4B_theory.m` 5/5.
2. **C₅ (11/4) — clean theory → 1+4η² CORROBORATED; code 11/4 is OPERATIONAL.** The cos-M Fourier of the
   ė kernel gives `1+4η²` (e⁰) `+ e₀(−2η−3η³)` (`verify_C5_theory.m` 3/3); the code's `1+11/4η²` differs
   (gap `−5/4η²`, e₀-term sign-flipped). Smoking gun (corroborated): C₅ is **not** the periodic partner
   of C₄-A under one orbit-average. Code-matched; theoretically UNRESOLVED / operational. **xmcof (−2/3):
   CLEAN-DERIVED (2026-06-03, `verify_xmcof_theory.m` 8/8).** The density-cubic partner of C₃ in `δℓ_D`
   (`secular_update.h:97-101`): from the Gauss `Ṁ` VOP (Thm 0.3.6) transverse `1/e` pole + the Lane-density
   power-4→3 antiderivative `∫sin f(1−η cos f)⁻⁴df = −(1/3η)(1−η cos f)⁻³`, the `−2/3 = 2×(−1/3)`, `1/(e₀η)`,
   `coef·B*` are **exact**; the `a₀` scalar (`xmcof_theory=a₀×code`) and `(1−η cos f†)⁻³→(1+η cos M)³` (thru
   `O(η)`) are operational. **Unlike C₅'s `11/4`, xmcof's leading `−2/3` is clean-theory-exact.**
3. **C₃ — mechanism CORRECTED (2026-06-02); factors CLEAN-DERIVED (2026-06-03, `verify_C3.m` 8/8).** C₃ is the
   **J₃×drag amplitude of the drag-periodic mean-anomaly correction** `δℓ_D` (sign-flipped to perigee;
   `secular_update.h:96-105`), **partner of `xmcof`** — **not** an e-vector term (`aycof`/`xlcof` are the
   C₃-free J₃ long-period terms; §C3.1-3.3 REFUTED, relabelled `[IDENTITY,not physics]`). Clean factors:
   `1/e₀` = Gauss `Ṁ` transverse **simple pole** (power 1, = xmcof's `1/e`, C3.7); `A₃₀/k₂=−2 j3oj2` conv;
   `n₀` = secular rate×t; `(A₃₀/k₂)sin i₀ = 4·aycof` (J₃ long-period e-vector amplitude, C3.8). `ξ⁵=coef·ξ`
   (drag density at `I^{(0,5)}` Hansen level, ∥ C₂ Part B — **not** ∂ρ/∂r gradient) + `cos ω` (e-vector over
   apsidal drift) structurally identified; exact J₃-LP×drag assembly code-matched, **partially operational**.
4. **xnodcf β₀²: RESOLVED (2026-06-02, theory).** The `β₀²` is the **angular-momentum drag
   normalisation**: `⟨Ġ/G⟩ = −B*⟨ρ|v|⟩ = −C₁·β²` (DERIVED — `⟨ρ|v|⟩/((a/μ)⟨ρ|v|³⟩) = β²` to O(e²)),
   whereas the energy rate `L̇/L = −C₁` carries no β. Hence `nodecf = ½C₁·xhdot1·(3+4β²)` carries the
   β² and matches the code at leading order; the code's exact `(7/2)β²` is the uniform-omeosq SGP4
   simplification (code − theory = O(e²)). My earlier "β₀⁴ three ways" was the *energy-rate* error.
   `verify_xnodcf_theory.m` 10/10.

**Trust model (binding).** Derive from theory; an agent may **confirm** a solution I have truly
derived, but may **not** be the **source** of a derivation; **OCR is never fact** (do not trust
vision — if it's a scan, it's OCR). The from-scratch derivations are now **DONE clean**: xnodcf
(2026-06-02), and **xmcof + C₃ (2026-06-03, `verify_xmcof_theory.m` 8/8 + `verify_C3.m` 8/8)** — all
from the trusted Phase-0 Gauss VOPs + Lane density, with operational pieces honestly dispositioned
(no OCR enshrined). C₄-B/C₅ leading terms remain corroborated-operational (do not re-hunt). **D-2
CLOSED** (2026-06-03). **PHASE 11 RESOLVED** (2026-06-03, `verify_phase11_fdagger.m` 8/8): the Lane
`f† = E` (eccentric anomaly), `r−s=(a−s)(1−η cos E)` EXACT ⇒ density exact, no f† error; the C2 Part-A
average over E reduces to the Phase-1 Lane integrals and **DERIVES** the §A.7 dropped term `(¾e²+3e²η²)`
exactly (= `(3/2)/η²(I04−2I03+I02) → (3/4)(1+4η²)`). The SGP4 `f†≈f` is just the E-vs-f equation-of-center
(the O(e²) drop). Remaining (optional/cosmetic): §A.7 Part-B `−5eη(4+3η²)` (analogous exact-over-E method);
C3 legacy §7 prose cleanup.

## Status ledger

| Phase | Coef | State | Closes | Verifier | Commit |
|---|---|---|---|---|---|
| 0-rev1 | Gauss VE + drag spec | ✅ DONE | D-10 | verify_phase0.m 32/32 | f790297 |
| 1 §1.0-1.3 | `I^{(0,m)}` m=1..6 | ✅ DONE (validated by 2.A) | — | (in 2.A) | 1aa15cd |
| 2.A | C₂ | ✅ DONE | D-9, D-3(C₂) | verify_phase2a.m 11/11 | 915b9c9 |
| 2.B | C₁ | ✅ DONE | Q-4 | verify_phase2a.m 13/13 | (2.B commit) |
| 2.C | Phase 1 closure (p=0; §1.4 = NOT-built decision) | ✅ DONE | D-3 (resolved C₂; final at C₄/C₅) | verify_phase1.m 13/13 | (2.C commit) |
| 3 | C₄ | ✅ code-matched (Part A DERIVED; Part B = SR3 = code, from-scratch deriv = **sealed LH79**) | D-3 cascade §8 | verify_C4.m 8/8 | (C4 commit) |
| 4 | C₅ | ◐ code-matched to SR3; 3/2 derived; **5/4 + e₀η³ sketch (open R09 §9a, sealed)** | C₅ 11/4 gap | verify_C5.m 3/3 | (C5 commit) |
| 5 | C₃ (J₃ × drag δℓ_D) | ✅ DONE; all factors clean-derived (1/e₀ Gauss-Ṁ pole; A₃₀/k₂; 4·aycof; cos ω resonance-free by exact cancellation; ξ⁵ coeff 1 = leading Hansen I^{(0,5)}(0)); residual = bounded C₂-family Lane-density η-truncation | D-4 | verify_C3.m 8/8 + verify_xmcof_C3_mitigations.m 14/14 | (2026-06-03) |
| 6 | D₂ | ✅ DONE | D-5 | verify_D2.m 3/3 | (D2 commit) |
| 7 | D₃, D₄ | ✅ DONE | cascade-D5 §11 | verify_D3D4.m 3/3 | (D3D4 commit) |
| 8 | t2cof…t5cof | ✅ DONE | D-6 | verify_tcofs.m 4/4 (in secular_taylor doc §T) | (tcofs commit) |
| 9 | corrections + cleanup | ✅ DONE — xmcof(D-8; −2/3 clean-derived 2026-06-03) xnodcf(D-7) delmo sinmo + D-1 D-11 + §7/§13/§15 doc-fixes | D-1✅ D-7✅ D-8✅ D-11✅ | verify_corrections.m 4/4, verify_xmcof_theory.m 8/8 | (2026-06-03) |

## Context

Each coefficient in `src/atmosphere/drag_coefficients.h` is proven from first
principles under **Standard 10** (`simplify(derived − code) == 0`; numerical-at-points
is NOT proof — `design/audit/2026_05_15_.../01_standards_and_schedule.md`). Phase 2.A
(C₂) landing retrospectively validated Phase 0-rev1 + Phase 1. The remaining
coefficients are only synopsis-level in the legacy main doc
(`sgp4_near_earth_drag_theoretical_basis.md`), where Pass-1 audit found open derivation
errors **D-1, D-3, D-4, D-5, D-6, D-7, D-8** (boxed forms match code; proofs do not
close). This plan re-derives each cleanly and closes each D-finding.

## The Phase 2.A template each phase mirrors

1. Single-pass-readable trace `design/derivations/sgp4_drag_<coef>_trace.md` (split if
   it grows): §scope → inputs → step-by-step derivation (every transform a **named
   theorem application or definition substitution**) → **code-match witness** → audit
   disposition. Model: `sgp4_drag_phase2a_C2_trace.md`.
2. Verifier `verify_<coef>.m` (Octave + symbolic), every check `simplify(lhs−rhs)==0`,
   ending in `simplify(<coef>_derived − <coef>_code) == 0` vs the literal code.
3. Audit-log disposition in `design/audit/2026_05_15_sgp4_drag_derivation_full_audit/`
   (close the D-finding, or escalate via newstk HR→PA→panel — never silently pick a
   convenient reading).
4. Commit on `session/2026-04-23`.

**Discipline (binding):** Standard 10; **sealed-room** (no direct reads of
BH61/Brouwer/Hoots/Lane scans — dispatch sub-agents, integrate cited summaries, verify
the assembly myself); born-digital sources preferred; legacy `deprecated/` derivations
are pedagogy-deprecated and **not citable** (derive fresh).

## Code targets (`drag_coefficients.h`)

| Coef | Lines | Form | Open finding |
|---|---|---|---|
| C₁ | 152 | `B*·C₂` | — (Q-4) |
| C₃ | 156 | `coef·ξ·A₃₀·n₀·sin i₀ /(k₂·e₀)`  (guard `e₀>1e-4`) | **D-4** |
| C₄ | 162-171 | `2n₀ coef1 a₀ β²·[ η(2+½η²)+e(½+2η²) − 2k₂ξ/(a₀ψ²)·( −3 con41 (1−2eη+η²(³⁄₂−½eη)) + ¾ sin²i₀(2η²−eη(1+η²)) cos2ω₀ ) ]` | **D-3** cascade |
| C₅ | 174-175 | `2 coef1 a₀ β²·(1+(11/4)(η²+eη)+eη³)` | C₅ 11/4 gap |
| D₂ | 178 | `4 a₀ ξ C₁²` | **D-5** |
| D₃ | 180 | `(17a₀+s)·D₂ξC₁/3` | cascade-D5 |
| D₄ | 181-182 | `½ temp_d a₀ξ (221a₀+31s) C₁`,  `temp_d=D₂ξC₁/3` | cascade-D5 |
| t2..t5cof | 185-191 | `(3/2)C₁ ; D₂+2C₁² ; ¼(3D₃+C₁(12D₂+10C₁²)) ; ⅕(3D₄+12C₁D₃+6D₂²+15C₁²(2D₂+C₁²))` | **D-6** |
| omgcof | 194 | `B*·C₃·cos ω₀` | — (follows C₃) |
| xmcof | 197 | `−⅔·coef·B*/(e₀η)` | **D-8** (ξ⁴) |
| xnodcf | 208-209 | `−(21/2)n₀k₂cos i₀ C₁/(a₀²β₀²)` | **D-7** |
| delmo/sinmo | 202-205 | `(1+η cos M₀)³ ; sin M₀` | — |

## Phases

### 2.B — C₁ (trivial)
`C₁ = B*·C₂` definitional (secular law `ȧ ≈ −2a₀B*C₂` ⇒ `tempa = 1 − C₁t`). Append
§A.9 to the C₂ trace; one check in `verify_phase2a.m`. Confirm Q-4.

### 2.C — Complete Phase 1 + Standard-10 §1.4 decision  ✅ DONE (2026-06-01)
**Scoping outcome:** the C₄ `⟨ė⟩` kernel collapses
(`e sin²f+(1+e cos f)(cos f+cos E) = 2(e+cos f)`) and, expanded in `e` at fixed `η`,
reduces to `h0 = cos f(1−η cos f)⁻⁴`, `h1 = (1+3η cos f)sin²f(1−η cos f)⁻⁵` — **purely
`p=0`**; `ψ⁷⟨h0⟩=η(2+½η²)`, `ψ⁷⟨h1⟩=½+2η²` = code:163. **Decision:** per Standard 10
the `p≥1` family (§1.4) is **NOT built** (no downstream code-match consumes it; D₂-D₄/
t-cofs have no orbit integral; C₃ → just-in-time if needed). Delivered: Phase 1 doc
§1.4 (decision)/§1.5 (two-route verification)/§1.6 (closed status); `verify_phase1.m`
(recurrence + η=0 + numerical quadrature, 13/13). D-3 resolved for C₂; finalized at
C₄/C₅.

### 3 — C₄ (hard) — close D-3 cascade §8
Orbit-average `ė = (β/na)[R sin f + T(cos f+cos E)]` (Thm 0.3.3) under drag
`R,T ∝ −B*ρ|v|(ṙ,rḟ)`, `|v|=(na/β)(1+e²+2e cos f)^{1/2}`, `ρ→(1−η cos f†)^{−4}`, Jacobian
`β³/(1+e cos f)²`. (i) Keplerian-drag → `η(2+½η²)+e(½+2η²)`. (ii) J₂-coupling
`−2k₂ξ/(a₀ψ²)` group with `con41`, reuse 2.A §A.6 `δr/r`; plus the long-period
`sin²i₀ cos2ω₀` term (e_vec/arg-perigee average; sealed-room cross-check vs SR3 p.11).
`sgp4_drag_C4_trace.md` (likely split) + `verify_C4.m` vs 162-171.

### 4 — C₅ — close 11/4 gap
Periodic O(e²) part of the C₄ `⟨ė⟩` expansion → `1+(11/4)(η²+eη)+eη³`; rigorous
`11/4=3/2+5/4`. `sgp4_drag_C5_trace.md` + `verify_C5.m` vs 174-175.

### 5 — C₃ (J₃, sealed-room) — close D-4
**Launch the J₃ sealed-room sub-agent early (‖ Phases 3-4).** It reads born-digital J₃
sources (Vallado SGP4 J₃ short-period, SR3, Brouwer γ₃) for `δr_{J₃} ∝ sin i₀ sin u`;
I verify the assembly. Derive `⟨ė⟩_{J₃-drag}` → `coef·ξ·A₃₀·n₀·sin i₀/(k₂·e₀)`.
**D-4:** the legacy "⟨sin²u⟩=½ absorbs the 2" *adds* a ½ (→/4) — find the correct chain
landing `/(k₂·e₀)`. Also `omgcof = B*·C₃·cos ω₀`. `sgp4_drag_C3_trace.md`+`verify_C3.m`.

### 6 — D₂ — close D-5
`a(t)=a₀ tempa²`, `tempa=1−C₁t−D₂t²−…`; D₂ from 2nd-order Taylor with `C₁(a)` varying
(`ṅ=−(3/2)(n/a)ȧ`, Thm 0.3.6). Audit showed exact match gives `D₂=−(5/4)C₁²−4a₀ξC₁²`
vs code `+4a₀ξC₁²`. Confirm SR3's `4a₀ξC₁²` is the retained `a-dependence-of-C₁` term
and the `C₁²` square is AFGP4→SGP4 truncation; code-match retained + **bound dropped
term** (don't silently match). Derivation landed in `sgp4_drag_secular_taylor_trace.md` §D.2 (not a separate D2 trace) + `verify_D2.m` vs 178.

### 7 — D₃, D₄
Continue `tempa` Taylor to t³/t⁴; derive `17a₀+s` and `221a₀+31s` (Q-9/Q-10 confirm the
`temp_d` algebra). Derivation in `sgp4_drag_secular_taylor_trace.md` §D.3 + `verify_D3D4.m` vs 179-182.

### 8 — t2cof…t5cof — close D-6
`templ` from `ℓ(t)=ℓ₀+n₀ templ` with `ṅ`-feedback; rewrite legacy §12 cleanly. Each cof
a polynomial in C₁,D₂,D₃,D₄. Derivation in `sgp4_drag_secular_taylor_trace.md` §T + `verify_tcofs.m` vs 185-191.

### 9 — Corrections + couplings + doc cleanup
xmcof (**D-8**, `coef=qoms4·ξ⁴` supplies the ξ⁴); xnodcf (**D-7**, adjudicate "β⁻³ into
C₁", ψ≠β); delmo/sinmo; **D-1** (ρ₀-in-B* convention); **D-11** (fix legacy main-doc §3
Thm 3.1 → `δr/r=−(3/2)(k₂/p²)β(3θ²−1)`, short-period form is 3× too small).
`sgp4_drag_corrections_trace.md`+`verify_corrections.m`.

### Phase 11 (Lane f† substitution) — RESOLVED 2026-06-03 (`verify_phase11_fdagger.m` 8/8)
**`f† = E` (eccentric anomaly):** `r − s = (a−s)(1 − η cos E)` is EXACT (P11.1), so the Lane density
`(r−s)⁻⁴ = ξ⁴(1−η cos E)⁻⁴` carries NO approximation. The C2 Part-A average over E (Jacobian
`dM=(1−e cos E)dE`, integrand `(1−e²cos²E)^{3/2}/[(1−e cos E)²(1−η cos E)⁴]`) reduces to the Phase-1 Lane
integrals and gives **exactly**: constant `I^{(0,4)}=(1+3/2η²)ψ⁻⁷` (P11.5), e-term `(2/η)(I04−I03)=η(4+η²)ψ⁻⁷`
⇒ code `eη(4+η²)` (P11.6), and the **§A.7 dropped O(e²) term `(3/2)/η²(I04−2I03+I02)=(3/4)(1+4η²)ψ⁻⁷` ⇒
`(¾e²+3e²η²)` — DERIVED** (P11.7). SGP4's `f†≈f` (true anomaly) is the E-vs-f equation-of-center, the
O(e²) drop. **Deferred (optional; not code-used):** §A.7 Part-B `−5eη(4+3η²)` (analogous exact-over-E
from the J₂-coupled integrand). Part-A is the representative completeness demonstration.

## Ordering / parallelism

```
C₂(done)→C₁(2.B)→Phase1 closure(2.C)→C₄(3)→C₅(4)        C₃(5) ‖ 3-4 (sealed-room early)
C₁ ─────────────────────────────────→ D₂(6)→D₃,D₄(7)→t-cofs(8)→cleanup(9)
```
Critical path: 2.B→2.C→C₄→C₅. D₂→D₃/D₄→t-cofs gated on C₁. Launch C₃'s J₃ dispatch
early (genuine parallelism / sealed-room read; never sub-agent the sequential
derivation work — `feedback_subagent_usage`).

## Per-phase verification (before commit)

1. `octave-cli --no-gui verify_<coef>.m` → all PASS; headline
   `simplify(<coef>_derived − <coef>_code) == 0`.
2. Foundation guards unchanged: `verify_phase0.m` 32/32, `verify_phase2a.m` 11/11.
3. Regression: `./build/Release/test_sgp4.exe` from repo root → **24/33, 518/623**
   (was 23/33, 506/623; the **low-perigee `s4` bug fix** 2026-06-03, commit `84fe9a7`,
   fixed sat 28350). The drag *re-derivation* itself validates code unchanged; the s4 fix
   is a separate code bug found via the test-failure investigation. The remaining 9 fails
   are deep-space/Lyddane/resonance (separate subsystem).
4. Update the audit-log row + the Status ledger above; commit.

## Build / verify environment

VS 2026 / MSVC 14.50 / v145. Octave 11.1.0 + symbolic 3.2.2 + SymPy 1.14.0 at
`/c/Program Files/GNU Octave/Octave-11.1.0/mingw64/bin/octave-cli.exe`. `test_sgp4.exe`
must run from the main repo CWD for relative data-path resolution.
