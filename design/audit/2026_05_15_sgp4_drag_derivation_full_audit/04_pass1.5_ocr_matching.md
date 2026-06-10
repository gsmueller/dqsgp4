## PASS 1.5 — OCR-Matching Audit (user directive 2026-05-15)

> *"You have quite a few equations that are matched to OCR documents. Please ensure
> that each of these equations has been symbolically rederived from theoretical sources."*

Per binding rules `feedback_trust_code_over_paper.md`, `feedback_derivation_first.md`,
`feedback_vet_all_claims.md`, `feedback_no_bias_toward_reference.md`, `reference_clean_sources.md`:

- OCR-scanned [B59], [BH61], [L65], [LC69], [LH79], [SR3] papers are **error-prone references**, not derivation sources.
- Born-digital sources preferred: ProofWiki / nLab / Wikipedia for canonical theorems; Lara 2021, Sneeuw 2022 for astrodynamics; the project's own BH61 cleanroom work (`sgp4_references/.../derivation/`) for the secular-rate / long-period-generator chain.
- "Matched to [SR3 p. 11]" is **not a proof**.

### Catalog of every boxed equation, with derivation status

Each row classifies the boxed equation by source:

- **DEF** — definition or convention (no derivation needed; only consistency check)
- **POST** — postulate of the physical model (e.g. Lane power-law density). Cannot be derived from anything; only its choice can be motivated.
- **PRIM** — primitive theorem from elementary mechanics / geometry / analysis. Derivable from first principles in 1-2 pages without specialist astrodynamics machinery.
- **BH61-CR** — derivable from the BH61 cleanroom work already done in `sgp4_references/.../derivation/`. The BH61 chain produces the J₂ secular rates and the long-period generator S₁* in born-digital form.
- **OCR-MATCH** — currently only "matches [primary source]"; **must be re-derived from theoretical sources**.

| Eq | § | Subject | Status | Derivation source |
|---|---|---|---|---|
| (1.1) | §1 | `B* := C_D·A/(2m)` (or with ρ_0?) | **DEF** (carries **D-1** dimensional error) | Standard ballistic-coefficient definition; cross-check Vallado §9.6 or Lane 1965. |
| (1.2) | §1 | Newton drag law `ddot r = -B*·ρ·\|v\|·v` | **POST** | Newtonian dynamics + drag postulate. |
| (2.1) | §2 | Lane power-law density `ρ(r) = ρ_0·((q-s)/(r-s))^τ` | **POST** | Lane 1965 model choice. Not derivable. |
| (3.1) | §3 | `⟨δr⟩_{J_2} = -(k_2/p)·(3θ²-1)/2` | **BH61-CR** | Brouwer 1959 long-period generator S₁ applied to r. Available in BH61 cleanroom `ch11a_long_period_generator.md` and `ch10_foundations`. **Re-derive symbolically here, citing the cleanroom result.** |
| (3.2) | §3 | `δρ_{J_2}/ρ = 2 k_2 ξ (3θ²-1)/p` for τ=4 | **PRIM** | Linear chain rule on (2.1) + (3.1). Re-derive explicitly. |
| (4.1) | §4 | `dM = β³/(1+e cos f)² df` | **PRIM** | Kepler's third law + elementary geometry. Re-derive explicitly. Replace D-2 intermediate. |
| (4.2) | §4 | Orbit-averaged `⟨ρ·(r/a)^m⟩` identity | **PRIM** + **OCR-MATCH** | Combination of (2.1), (4.1), and the Lane f† substitution. Re-derive from first principles. Fix D-3. |
| (4.3) | §4 | `I^{(p,m)}` definition | **DEF** | — |
| (4.4) | §4 | `I^{(0,4)} = (2+η²)/(1-η²)^{7/2}` | **OCR-MATCH** | Residue calculus on the unit circle. Re-derive symbolically. |
| (4.5..) | §4 | `I^{(1,4)}, I^{(2,4)}` (referenced but not boxed) | **OCR-MATCH** | Re-derive. |
| (5.1) | §5 | C₂ boxed form | **OCR-MATCH** (currently traced to `deprecated/020`; status of 020 unverified) | Direct assembly from §3 (J₂ density coupling) + §4 (Lane integrals). Re-derive from first principles. |
| (6.1) | §6 | `C₁ := B*·C₂` | **DEF** | — |
| (6.2) | §6 | `a(t) = a_0''·(1-C₁τ+O(τ²))²` | **DEF** (construction of the time-polynomial form) | Justify the squaring convention. |
| (6.3) | §6 | `⟨dE/dt⟩` | **PRIM** | Energy / momentum identity. Trivial. |
| (7.1) | §7 | `δr_{J_3}` | **BH61-CR** | Brouwer J₃ long-period generator. Available in cleanroom. |
| (7.2) | §7 | C₃ boxed form | **OCR-MATCH** | Re-derive from §7.1 + §4 + Lagrange VE for ė. Fix D-4. |
| (8.1) | §8 | `Ġ\|drag` Lagrange VE | **PRIM** | Standard variational equations (Gauss form). |
| (8.2) | §8 | C₄ boxed form | **OCR-MATCH** | Re-derive from §8.1 + §4 + §3.2. |
| (9.1) | §9 | `Ṁ\|drag` Lagrange VE | **PRIM** | Standard variational equations. |
| (9.2) | §9 | C₅ boxed form | **OCR-MATCH** | Re-derive from §9.1 + §4 + cubic Taylor in eη. (This is the R09 §9a item.) |
| (10.1) | §10 | `D₂ = 4 a_0'' ξ C₁²` | **OCR-MATCH** + **D-5 critical** | Re-derive Taylor expansion of `a(t)/a_0''` to order τ². Resolve the missing `(1/2)C₁²` term — either the SR3 formula is exact (and the doc's proof is wrong) or it is approximate (and the order should be quantified). |
| (10.2)-(10.7) | §10 | Taylor steps | **PRIM** | Elementary. |
| (11.1) | §11 | `D₃ = (1/3)(17 a_0'' + s)·D₂·ξ·C₁` | **OCR-MATCH** | Re-derive Taylor to order τ³. (This is the R09 §11a item.) |
| (11.2) | §11 | `D₄ = (1/2) D₃ ξ a_0'' (221 a_0'' + 31 s) C₁ / (17 a_0'' + s)` | **OCR-MATCH** | Re-derive Taylor to order τ⁴. (R09 §11a continued.) |
| (12.1) | §12 | `ℓ(t)` propagation form | **DEF** | — |
| (12.2) | §12 | t₂cof ... t₅cof | **OCR-MATCH** + **D-6 sufficiency** | Re-derive from binomial series `(1-x)^{-3/2}`. (R09 §12a.) |
| (13.1) | §13 | xnodcf = `-(21/2) n_0'' k_2 cos i_0 / (a_0''² β_0²) · C₁` | **OCR-MATCH** + **D-7 suspicion** | Re-derive `Ω̇(a) ↦ Ω̇(t)` Taylor; resolve "absorb β⁻³ into C₁" question. |
| (14.1) | §14 | `a_xN, a_yN, ℓ'` decomposition | **DEF** | — |
| (14.2) | §14 | aycof = `-(1/4)·(A_30/k_2)·sin i_0` | **BH61-CR** | Long-period generator S₁* projected on e_y. Available in cleanroom. |
| (14.3) | §14 | xlcof = `(1/8)·(A_30/k_2)·sin i_0·(3+5 cos i_0)/(1+cos i_0)` | **BH61-CR** | Long-period generator S₁* projected on ℓ. Available in cleanroom. |
| (15.1) | §15 | `ω(t)` propagation form | **DEF** | — |
| (15.2) | §15 | omgcof = `B*·C₃·cos ω_0` | **PRIM** | Integration of ω̇\|_{J_3-drag} from (7.4). |
| (15.4) | §15 | xmcof = `-(2/3)·(q-s)⁴·B*/(e₀η)` (per doc) — **D-8 says code has additional ξ⁴** | **OCR-MATCH** | Re-derive M̈\|_{cubic} from (9.1) + cubic Taylor. Fix D-8. |
| (15.5) | §15 | delmo, sinmo | **DEF** | — |

### Tally

- DEF: 8 boxed equations (need consistency check only)
- POST: 2 (no derivation possible)
- PRIM: 6 (elementary; 1-2 pages each)
- BH61-CR: 4 (already done in cleanroom; cite + adapt)
- **OCR-MATCH**: **14 boxed equations that need full symbolic re-derivation from theoretical sources**

These 14 are exactly the "**matched to [SR3 p. X]**" / "**matches [LH79 p. Y]**" rows. They constitute the bulk of the SGP4 drag theory machinery.

## Re-derivation queue

Order is dictated by dependency:

| Phase | What | Depends on | OCR-MATCH equations closed |
|---:|---|---|---|
| **0** | Foundations: Kepler geometry, Lagrange VE in Gauss form, Brouwer S₁ J₂/J₃ generator pull-back (cite BH61 cleanroom). | — | (4.1) [fix D-2], (3.1) [via BH61], (7.1) [via BH61] |
| **1** | Lane integrals: residue calculus for `I^{(p,m)}` at τ=4 (m=2,3,4, p=0,1,2). Fix D-3 by re-deriving (4.2). | Phase 0 | (4.2), (4.4), (4.5..) |
| **2** | C₂ assembly: §3.2 + §4 Lane integrals → (5.1). | Phases 0, 1 | (5.1) |
| **3** | C₃: J₃ density coupling + Lagrange VE for ė → (7.2). Fix D-4 (the /2 vs /4 discrepancy). | Phases 0, 1, 2 | (7.2) |
| **4** | C₄: Lagrange VE for Ġ + Lane integrals + §3.2 → (8.2). | Phases 0, 1, 2 | (8.2) |
| **5** | C₅: cubic Taylor in eη + Lagrange VE for Ṁ → (9.2). **This is R09 §9a.** | Phases 0, 1, 2 | (9.2) |
| **6** | D₂ Taylor expansion. Resolve D-5: is `D₂ = 4 a_0'' ξ C₁²` exact or approximate? If approximate, quantify the order of the dropped `(1/2)C₁²` and add it to the error catalog. | Phase 2 | (10.1) |
| **7** | D₃, D₄ Taylor expansion to orders τ³, τ⁴. **This is R09 §11a.** | Phase 6 | (11.1), (11.2) |
| **8** | t₂cof … t₅cof via binomial series `(1-x)^{-3/2}`. **This is R09 §12a.** Rewrite §12 from scratch (close D-6). | Phase 7 | (12.2) |
| **9** | xnodcf Taylor. Resolve D-7 (β^{-3} absorption claim). | Phases 0, 1, 2, 6 | (13.1) |
| **10** | xmcof from cubic Taylor + Lagrange VE for Ṁ. Fix D-8 (ξ⁴ factor). | Phase 5 | (15.4) |
| **11** | Lane f† formalization. **This is R09 §4a.** Define f† as a series in (f, e); bound truncation error. | Phase 1 | Replaces the f† caveat in (4.2). |

Phase 1 closes D-3 (the master orbit-averaging identity).
Phase 6 closes D-5 (D₂).
Phases 7-8 close D-6 (§12 sufficiency).
Phases 9-10 close D-7, D-8.

The R09 sub-items §4a, §9a, §11a, §12a are subsumed: §4a = Phase 11, §9a = Phase 5, §11a = Phase 7, §12a = Phase 8. They cannot be cleanly done without Phases 0-1-2 first because the master orbit-average identity (D-3) is broken.

