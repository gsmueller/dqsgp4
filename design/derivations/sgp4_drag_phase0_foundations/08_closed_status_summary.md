## §0.7 Closed-Status Summary for Phase 0 (rev1)

Phase 0-rev1 includes the **D-10 fix**: the drag specialization in Theorem 0.4.1
and the closed-form rates in Theorem 0.4.2 now reflect the **Newtonian quadratic
drag** `F_drag = -B*·ρ·|v|·v` from main-doc Postulate 1.2, with the `|v|` factor
and the `(1 + e² + 2 e cos f)^{3/2}` exponent restored. Theorem 0.5.3's
orbit-averaged-drag setup is updated to match (the integrand carries the
`^{3/2}` factor).

**Theorems established (with full multi-pass-aligned proofs):**

| Theorem | Subject | Closes (or supersedes) |
|---|---|---|
| **0.2.2** | Conic section `r = p/(1+e cos f)` | — |
| **0.2.3** | Angular momentum `r²ḟ = h = √(μp) = n a² β` | — |
| **0.2.4** | Kepler 3rd `n² a³ = μ`, `n a² β = h` | — |
| **0.2.5** | `dM/df = β³/(1+e cos f)²` | **D-2** (main-doc Eq. 4.1 intermediate-form transcription) |
| **0.2.6** | Orbit-average operator `⟨·⟩_M` | — |
| **0.3.2** | Gauss VE `ȧ = (2/(nβ))·[e R sin f + T(1+e cos f)]` | — |
| **0.3.3** | Gauss VE `ė = (β/(na))·[R sin f + T(cos f + cos E)]` | — |
| **0.3.4** | Gauss VE `Ω̇, i̇` from `r × F` | — |
| **0.3.5** | Gauss VE `ω̇ = ... − Ω̇ cos i` | — |
| **0.3.6** | Gauss VE `Ṁ = n − 2 r R/(n a²) + (β²/(n a e))·[R cos f − T sin f (1+r/p)]` | (prerequisite for closing **D-5**) |
| **0.4.1 (rev1)** | Drag specialization `R_drag = -B*·ρ·|v|·ṙ`, `T_drag = -B*·ρ·|v|·r ḟ`, `N_drag = 0`; `|v| = (n a/β)·√(1+e²+2e cos f)` | **D-10** (linear-vs-quadratic drag form) |
| **0.4.2 (rev1)** | Closed-form drag rates of all six elements; principal result `ȧ_drag = -(2 n B*·ρ·a²/β³)·(1+e²+2e cos f)^{3/2}` | (prerequisite for Phase 2.A) |
| **0.5.1** | Energy identity `dE/dt = (μ/(2a²))·ȧ` | — |
| **0.5.2** | Delaunay momentum `L = √(μa)`, `dL/dt = (n a / 2)·ȧ` | — |
| **0.5.3 (rev1)** | Setup of orbit-averaged drag rate `⟨ȧ⟩_drag` with corrected integrand `(1+e²+2e cos f)^{3/2}·ρ(r)·β³/(1+e cos f)²` | (input to Phase 2.A) |

**Alignment-to-SGP4 remarks established:**

| Topic | Standard discharged |
|---|---|
| `n → n₀''`, `a → a₀''` Brouwer-recovered | (0.6.2) |
| VOP / osculating-elements framework | (0.2.2 alignment) |
| Two contexts of orbit-average (epoch-frozen vs instantaneous-osculating) | (0.2.6 alignment) |
| `Ω̇_drag = i̇_drag = 0` | (0.3.4 alignment + 0.4.1-rev1 / Theorem 0.4.2-rev1 Step 4) |
| Coordinate singularities at `i = 0` and `e = 0` | (0.3.4.17, 0.3.5.14) |
| `β` ≠ `ψ` — distinct Keplerian vs Lane factors | (0.6.5) — closes **D-7** in advance |
| Non-rotating-atmosphere model choice (accuracy `A-D1`) | (0.6.6) |
| Time-scale invariance of theorem forms | (0.6.1) |
| **Newtonian quadratic drag `F_drag = -B*·ρ·|v|·v`** | (§0.4 opening, 0.4.1-rev1 hypothesis (H₁); §0.4 alignment bullet (d)) — closes **D-10** |

**Lemmas established:**

| Lemma | Subject |
|---|---|
| 0.3.3.a | Eccentricity vector identity `e_vec = (v×h)/μ − r̂` |
| 0.3.5.a | Line-of-nodes rotation `(dN̂/dt) · M̂ = Ω̇ cos i` |

**Audit findings closed by Phase 0-rev1:**

- **D-2** — Theorem 0.2.5 supersedes the main-doc Eq (4.1) intermediate-form
  transcription error.
- **D-7** — Remark 0.6.5 explicitly distinguishes `β² = 1 − e²` (Keplerian) from
  `ψ² = |1 − η²|` (Lane density residue) in advance of the §13 main-doc claim
  that would have confused them.
- **D-10** — Theorem 0.4.1-rev1 + Theorem 0.4.2-rev1 + Theorem 0.5.3-rev1 restore
  the `|v|` factor for Newtonian quadratic drag; the `verify_phase0.m`
  Postulate-1.2 match check (added under Block A4 of the remediation plan)
  guards against recurrence.

**Phase 0-rev1 verifier `verify_phase0.m`:** updated under Block A4 with the
Postulate-1.2 symbolic match check (`R_drag·r̂ + T_drag·t̂ + N_drag·n̂ ≡
−B*·ρ·|v|·v`) and the corrected closed-form check for Theorem 0.4.2 (0.4.2.1)
with the `(1+e²+2e cos f)^{3/2}` factor. The verifier is run under Block A5.

**Out of scope for Phase 0:** every C-coefficient (`C₁..C₅`), every D-coefficient
(`D₂..D₄`), `t-cofs`, `xnodcf`, `xlcof`, `aycof`, `omgcof`, `xmcof`. Those are
Phases 2-10. Phase 0 produces only the primitives that Phases 1-10 consume.

**Out of scope for Phase 0 (alignment standard 9-B):** all implementation
choices (Horner-form coefficient assemblies, exact-rational arithmetic,
`ratio<T>(p, q)`, error-bound integration, etc.). These are unaffected by the
theoretical Phase 0 results and will be addressed at each downstream Phase.

**Status:** Phase 0-rev1 derivations complete. Per Standard 10, Phase 0 is
**retrospectively validated** only when Phase 2.A (Block B of the remediation
plan) successfully derives the SGP4 code's C₂ expression from these primitives
with symbolic equality at `drag_coefficients.h:146-149`.
