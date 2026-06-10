### §5 Fundamental drag rate C₂ (synopsis)

- **Q-3** Theorem 5.1 boxed form matches the code (`drag_coefficients.h:146-149`) under the identifications `coef1 = (q-s)^4·ξ^4·ψ^{-7}` and `half_J2 = k₂ = J₂/2`. ✓
- **N-9** Line 201 — references `deprecated/020_c2_drag_integral_derivation.md` as the canonical full C₂ derivation. The file lives under `deprecated/` — why? Either subsumed or retired. **If subsumed, the master proof of (5.1) is not in this document** — the doc is a synopsis only. Pass 2 must verify 020 still exists and is consistent; otherwise §5 has no proof.

### §6 Drag Polynomial Rate C₁

- **Q-4 — CLOSED (2026-06-01, Phase 2.B)** Definition 6.1: `C₁ := B*·C₂`. Code:152
  matches. ✓ Now **derived** (not only quote-checked) in
  `design/derivations/sgp4_drag_phase2a_C2_trace.md` §A.9: the secular-decay law
  `a(t)=a₀ tempa²` (code:14) gives `ȧ(0)=−2a₀C₁`, the master drag rate (A.1.3→A.8.1)
  factors as `⟨ȧ⟩=−2a₀B*C₂`, hence `C₁=B*C₂`. `verify_phase2a.m` checks A.9.2
  (prefactor factor-out) + A.9.3 (`simplify(C₁_derived−C₁_code)=0`), 13/13 PASS.
- **Q-5** Proposition 6.2: `a(t) = a_0''·[1 - C₁(t-t_0) + O((t-t_0)²)]² = a_0''·[1 - 2C₁(t-t_0) + O((t-t_0)²)]`.
  - Squaring identity is correct. ✓
- The energy-loss equation (6.3) is plausible but its derivation is one line of hand-wave. ✓-with-flag.

### §7 J₃ Coupling Coefficient C₃

- **Q-6** Eq (7.1) algebra: `δr_{J_3} = -(A_{3,0}/(2k_2))·(sin i/p)·sin u · p = -(A_{3,0} sin i_0/(2k_2))·sin u`. The two `p`'s cancel. ✓
- **D-4 — CLOSED (2026-06-01, Phase 5/C₃; `design/derivations/sgp4_drag_C3_trace.md`, `verify_C3.m` 6/6)**
  The boxed form (7.2) and the code (`drag_coefficients.h:156`) are **correct**; only the
  legacy *proof prose* (lines 282-285) is wrong. The error: it invoked `⟨sin²u⟩=½` to "absorb
  the factor of 2", which would give `/(4k₂e₀)` (off by **4×**) — and there is no `⟨sin²u⟩` in
  the chain at all. **Correct factor-2 cancellation:** the Brouwer `A₃₀/(2k₂)` prefactor's 2 is
  killed by the **Gauss-VE factor 2** of the eccentricity-vector rate `de/dt` (Phase 0 Thm 0.3.3:
  `(cos f + cos E)|_{e→0} = 2cos f`), giving exactly `/(k₂e₀)`. Convention (load-bearing):
  `A₃₀ = −J₃a_E³` (no ½, model_selector.h:178) ⇒ `A₃₀/k₂ = −2 j3oj2 a_E` = Vallado `cc3`.
  J₃ inputs born-digital (sealed-room dispatch). **Legacy §7 prose fix → Phase 9 cleanup.**
  Original finding retained below for the record.

- **D-4 (MAJOR) [original]** Theorem 7.1 proof sketch (lines 273-291) is inconsistent with the boxed form.

  Stated boxed form (7.2): `C₃ = (q_0-s)^4·ξ⁵·A_{3,0}·n_0''·sin i_0 / (k₂·e_0)`.
  
  Stated proof line (7.4): `⟨ė⟩_{J3-drag} ∝ B*·(q-s)^4·ξ^4 · (A_{3,0}/(2k_2)) · (ξ·sin i_0/e_0) · n_0'' · sin ω`. Collecting:
  
  `B*·(q-s)^4·ξ⁵·A_{3,0}·n_0''·sin i_0·sin ω / (2 k_2 e_0)`.
  
  Then the proof claims "absorbing the factor of 2 from the trigonometric average `⟨sin²u⟩ = 1/2`"; but ⟨sin²u⟩ = 1/2 is a factor of *1/2*, which would *introduce* a 1/2 (i.e. *increase* the denominator's 2 to a 4), not *cancel* the existing 2.
  
  So the proof, taken literally, produces `C₃ = (q-s)^4·ξ⁵·A_{3,0}·n_0''·sin i_0 / (4 k_2 e_0)`, off by a factor of 4 from the boxed form.
  
  But the code (line 156) has `C₃ = (q-s)^4·ξ⁵·A_{3,0}·n_0''·sin i_0 / (k_2·e_0)` — agreeing with the boxed form, not with the proof.
  
  **Verdict (Pass 1):** Q-FAIL — the proof and the boxed form (and the code) are inconsistent. The boxed form matches the code, so the proof needs to be re-written or re-derived.

### §8 Eccentricity Decay Rate C₄

- **Q-7** Theorem 8.1 boxed form: cross-checks against code `drag_coefficients.h:162-171` with the identifications used in §5. ✓
- **Cascade — RESOLVED (2026-06-01, Phase 3/C₄; `sgp4_drag_C4_trace.md`, `verify_C4.m` 8/8)**
  D-3 is resolved (Phase 2.A §A.2 β³ Jacobian). C₄ **Part A (Keplerian) is DERIVED** from the
  Gauss `ė` (Thm 0.3.3) via the `p=0` Lane family `I^{(0,m)}` (no `p≥1` needed): `ψ⁷⟨h0⟩=η(2+½η²)`,
  `ψ⁷⟨h1⟩=½+2η²` = code:163. **Part B (J₂ coupling) is code-matched to SR3 = operational SGP4**
  (= code:164-171). **STATUS 2026-06-02 (trusted-theory):** my from-scratch drag×J₂ derivation gives
  the isotropic term **order-η**, the code is **order-1** — these disagree and are **UNRESOLVED from
  trusted sources** (my derivation may be incomplete, or the code may differ from correct theory). An
  OCR-sourced "operational erratum" explanation is an **untrusted hypothesis**, not asserted. Match
  the code for bit-compatibility (see trace §C4.B). The §4 Lane-integral
  inheritance is sound (the master orbit-average uses the corrected β³ Jacobian).
- **Cascade [original]** §8 proof "(i) Keplerian-drag part" invokes "the Lane integrals `I^{(p,m)}` from §4". Those Lane integrals come from Lemma 4.2 / Eq (4.2) which carries **D-3**. So C₄'s derivation **inherits** the error from §4 unless Pass 2 resolves D-3 favourably.

### §9 Second-Order Mean-Anomaly Correction C₅

- **Q-8** Theorem 9.1 boxed: cross-checks against code `drag_coefficients.h:174-175`. ✓
- **N-10** Proof sketch claims `11/4 = 3/2 + 5/4` ✓ arithmetically, but the listing of "three contributions" is confusing: a separate "1/2 η(1+η²) from higher orders" is mentioned but this term does **not** appear in either 3/2 or 5/4 — it must instead correspond to the `e_0 η³` term in the boxed form (which is listed separately as the "third bullet"). The bullet structure mixes the 11/4 decomposition with separate higher-order terms.
- **Status:** R09 §9a is the planned rigour-extension. Pass 1 confirms the gap is real.

### §10 Taylor Expansion of a(t) — D₂

- **D-5 — CLOSED (2026-06-01, Phase 6/D₂; `sgp4_drag_secular_taylor_trace.md`, `verify_D2.m` 3/3)**
  The code `D₂ = 4a₀ξC₁²` (`drag_coefficients.h:178`) is **exact**, not a truncation. **Root
  cause of the audit's residual:** the re-derivation below dropped the a-dependence of the mean
  motion `n ∝ a^{−3/2}` (and the explicit `a` in C₂). The correct decay rate is `ȧ = −2a C₁(a)`
  with `C₁(a) ∝ (a−s)^{−4}a^{−1/2}`, so `d(ln C₁)/da = −4ξ − 1/(2a)`; the `−1/(2a)` makes the
  stray `C₁²` cancel against `tempa²`, giving `D₂ = 4a₀ξC₁²` exactly (`verify_D2.m` D2.1/D2.2).
  Dropping `n(a)` (`d ln C₁/da = −4ξ`) reproduces the audit's `−½C₁²` residual **exactly**
  (`verify_D2.m` D2.3) — confirming it was an omission, not a code error. Original below.

- **D-5 (MAJOR) [original]** Theorem 10.1 proof (lines 402-447) does not close cleanly.

  Doc's derivation:
  
  (10.2) `a(τ)/a_0'' = 1 - 2C₁τ - 2D₂τ² + C₁²τ² + O(τ³)`. ✓ (squaring identity)
  
  (10.3) `ȧ/a_0'' = -2C₁ - 4D₂τ + 2C₁²τ + O(τ²)`. ✓ (τ-derivative of 10.2)
  
  Re-derivation of `ȧ/a`: from `a/a_0'' = (1 + x)²` with `x = -C₁τ - D₂τ² - …`,
  
    `ln(a/a_0'') = 2 ln(1+x) = 2(x - x²/2 + x³/3 - …)`
    
    `ȧ/a = d/dτ ln(a/a_0'') = 2 dx/dτ - x·dx/dτ + …`
    
    `dx/dτ = -C₁ - 2D₂τ + O(τ²)`,  `x = -C₁τ + O(τ²)`,  `-x·dx/dτ = (C₁τ)(C₁) + O(τ²) = C₁²τ + O(τ²) … wait sign:`
    
    `-x · dx/dτ = -(-C₁τ)·(-C₁) + O(τ²) = -C₁²τ + O(τ²)`.
    
    So `ȧ/a = 2(-C₁ - 2D₂τ) + (-C₁²τ) + O(τ²) = -2C₁ - 4D₂τ - C₁²τ + O(τ²)`.
  
  Now match to `ȧ/a = -2 C₁(a(τ))` with (10.6) `C₁(a(τ)) ≈ C₁(1 - 2C₁τ(1 + 4 a_0'' ξ))`:
  
    `-2 C₁(a(τ)) = -2 C₁ + 4 C₁² τ (1 + 4 a_0'' ξ) = -2 C₁ + 4 C₁² τ + 16 a_0'' ξ C₁² τ`.
  
  Equating:
  
    `-2C₁ - 4D₂τ - C₁²τ = -2 C₁ + 4 C₁² τ + 16 a_0'' ξ C₁² τ`
    
    `-4 D₂ - C₁² = 4 C₁² + 16 a_0'' ξ C₁²`
    
    `-4 D₂ = 5 C₁² + 16 a_0'' ξ C₁²`
    
    `D₂ = -(5/4) C₁² - 4 a_0'' ξ C₁²` (i.e., **not** the doc's claimed `D₂ = 4 a_0'' ξ C₁²`).
  
  The doc's matching step (line 444) goes from `-4 D₂ + 2 C₁² = 4 C₁² + 16 a_0'' ξ C₁²` to `D₂ = 4 a_0'' ξ C₁²` by canceling 2C₁² on the LHS with one of the 4C₁² terms on the RHS, **leaving a stray (1/2)C₁²** in the doc's own working that is silently dropped ("with our sign convention…").
  
  **Even the doc's own matching is inconsistent**: −4D₂ + 2C₁² = 4C₁² + 16a₀ξC₁² gives D₂ = −(1/2)C₁² − 4a₀ξC₁², not the boxed D₂ = +4a₀ξC₁².
  
  My re-derivation produces a different (still nonzero) `−(5/4)C₁²` residual, depending on whether `ȧ/a_0''` or `ȧ/a` is being matched.
  
  **Verdict (Pass 1):** DERIVATION-ERROR (CRITICAL). §10 Theorem 10.1 proof does not produce its own boxed form. The boxed form matches `[SR3] p. 11`, so the code is right by reference; the proof is wrong by algebra.
  
  Hypothesis for Pass 2: the SR3 formula `D₂ = 4 a_0'' ξ C₁²` is itself an *approximation* that drops the additional `C₁²` term. Lane-Hoots may explicitly note this. The doc should make the approximation explicit (with order-of-magnitude bound), not silently match the formula.

### §11 Higher-Order Drag — D₃, D₄

- **Q-9** Theorem 11.1: `D₃ = (1/3)·(17 a_0'' + s)·D₂·ξ·C₁`. Cross-check against `drag_coefficients.h:179-180`: code has `temp_d = D₂·ξ·C₁/3; D₃ = (17 a_0'' + s)·temp_d` ✓.
- **Q-10** Theorem 11.2: `D₄ = (1/2)·D₃·ξ·a_0''·(221 a_0'' + 31 s)·C₁/(17 a_0'' + s)`. The second form `(2/3)·a_0''²·ξ³·(221 a_0'' + 31 s)·C₁⁴` is obtained by substituting D₃. Let me re-check algebra:
  
    `D₄ = (1/2)·[(1/3)·(17 a_0''+s)·D₂·ξ·C₁]·ξ·a_0''·(221 a_0''+31s)·C₁ / (17 a_0''+s)`
    `   = (1/6)·D₂·ξ²·a_0''·(221 a_0''+31s)·C₁²`
    `   = (1/6)·[4·a_0''·ξ·C₁²]·ξ²·a_0''·(221 a_0''+31s)·C₁²`
    `   = (2/3)·a_0''²·ξ³·(221 a_0''+31s)·C₁⁴`. ✓
  
  But code `drag_coefficients.h:181-182` reads `D₄ = (1/2)·temp_d·a_0''·ξ·(221 a_0'' + 31 s)·C₁` where `temp_d = D₂·ξ·C₁/3`. So `D₄ = (1/2)·(D₂·ξ·C₁/3)·a_0''·ξ·(221 a_0''+31 s)·C₁ = (1/6)·D₂·ξ²·a_0''·(221 a_0''+31s)·C₁²`. ✓ (matches my expansion)
  
- **Cascade** §11 proofs are "sketches" — explicitly flagged in Part VII Open Theoretical Gaps as the work R09 §11a is supposed to formalize. Pass 1 NOTED.
- **Cascade-from-D-5 — RESOLVED (2026-06-01, Phase 7; `verify_D3D4.m` 3/3)** D₂ is **exact**
  (D-5 closed), so D₃/D₄ are exact. Both are derived from one decay-ODE Taylor
  (`sgp4_drag_secular_taylor_trace.md` §D.3): `D₃ = (17a₀+s)·temp_d`,
  `D₄ = ½ temp_d a₀ ξ (221a₀+31s) C₁`, `temp_d = D₂ξC₁/3`, with `(17a₀+s)`/`(221a₀+31s)`
  emerging from the higher `a`-derivatives of `C₁(a) ∝ (a−s)^{−4}a^{−1/2}`. Q-9/Q-10 confirmed.
- **Cascade-from-D-5 [original]** §11 D₃/D₄ depend on D₂. If D₂ is approximate (per D-5 hypothesis), then D₃/D₄ inherit the approximation.

### §12 Mean-Longitude Polynomial — t2cof … t5cof

- **D-6 — CLOSED (2026-06-01, Phase 8; `sgp4_drag_secular_taylor_trace.md` §T, `verify_tcofs.m` 4/4)**
  All four t-cofs are the `t²..t⁵` Taylor coefficients of `templ = ∫₀^t(n/n₀ − 1)dt'` with
  `n = n₀ tempa^{−3}` (since `a = a₀ tempa²`, `n ∝ a^{−3/2}`). Expanding `tempa^{−3}−1 =
  3x+6x²+10x³+15x⁴` (`x = C₁t+D₂t²+D₃t³+D₄t⁴`) and integrating gives `t2cof=(3/2)C₁`,
  `t3cof=D₂+2C₁²`, `t4cof=¼(3D₃+C₁(12D₂+10C₁²))`, `t5cof=⅕(3D₄+12C₁D₃+6D₂²+15C₁²(2D₂+C₁²))` —
  **all code-matched** (`simplify=0` vs `:185-191`). The legacy §12 stream-of-consciousness was a
  failure to state the one-line integral; t3/t4/t5cof are now derived. Original below.

- **D-6 (CRITICAL) [original]** §12 proof (lines 526-548) is **incomplete and confused**.
  
  The proof text contains:
  
  - Line 526: "actually re-examining…" — the author identifies that the τ¹ term contributes to `dot M` and not to templ.
  - Line 528: "Actually, [SR3] uses a slightly different decomposition…"
  - Line 536: "this is also linear. Let me approach differently."
  - Line 544: "This is messy. Let me just compute the t_2 cof directly from the structure."
  - Lines 546-547: a heuristic re-derivation of t_2 cof using `ṅ = -(3/2)(n/a)·ȧ`.
  - Line 548: "∎ (sketch for t₃, t₄, t₅ cof: continue the same expansion to higher powers of τ, with corrections from D₂, D₃, D₄ entering at the corresponding orders)."
  
  This is **not a proof**. It is a stream-of-consciousness mid-derivation re-think with "actually" and "let me approach differently" left in the text. The t_3, t_4, t_5 cof's are not derived at all.
  
  **Verdict (Pass 1):** SUFFICIENCY-GAP + TRANSCRIPTION-ERROR. The text needs to be entirely rewritten before §12a (R09 sub-item) extends it. The exposed working "actually" / "let me approach differently" / "this is messy" should not appear in a derivation document.
  
  Severity: HIGH for documentation rigour; LOW for numerical correctness (the code matches SR3 by transcription).

### §13 Drag-Gravity RAAN Coupling — xnodcf

- **Q-11** Theorem 13.1 algebra cascades correctly. ✓
- **D-7 — REOPENED (2026-06-02, independent verification).** The β-power is **NOT closed**. The
  β≠ψ / no-C₁-absorption part IS correct: `C₁ = B*C₂`, `C₂ ∝ ψ⁻⁷=(1−η²)⁻⁷ᐟ²`, `η=a₀''e₀ξ ≠ e₀`, so
  C₁ carries no β-power (Phase 0 §0.6 Remark 0.6.5) — the original "absorb β₀⁻³ into C₁" was wrong.
  **But the 2026-06-01 "the β₀² is just Vallado's explicit `omeosq` convention" disposition is
  WITHDRAWN:** my secular node-rate Taylor reproduces the code's leading term with **β₀⁴, not β₀²**,
  by **three** independent routes (a-Taylor, L/G/H momentum, n/p). The correct derivation is the
  (L,G,H) Delaunay-momentum node update with drag **angular-momentum** decay (`G=Lβ` decays — β is
  NOT fixed); whether it yields the code's β₀² is **UNRESOLVED** (could not reproduce — a
  clean-session from-scratch task, method only). `verify_corrections.m`'s xnodcf check is a
  code↔Vallado-assembly identity (public), **not** a derivation of the β₀². Original below.

- **D-7 [original]** Line 592 — "SR3 simplifies this to `-(21/2) k_2 n_0'' cos i_0 C₁ / (a_0''² β_0²)` by absorbing β_0^{-3} into the C₁ definition". 
  
  This is suspicious. C₁ = B*·C₂, and C₂ contains `ψ^{-7}` where `ψ² = |1-η²|`. η = a·e·ξ is **not** equal to β = √(1-e²). So ψ ≠ β.
  
  Claim: β_0^{-3} ends up absorbed into C₁ via the ψ^{-7} factor. But ψ depends on η, not on β directly. The two coincide only when e=0 (then η=0, ψ=1; also β=1). For e≠0 they differ.
  
  **Verdict (Pass 1):** suspected DERIVATION-ERROR or APPLICATION-ERROR. The "absorbed into C₁" claim looks incorrect. Pass 2 against [SR3] p. 12.

### §14 Long-Period Periodics — xlcof, aycof

- **Q-12** Theorem 14.1 algebra holds. ✓
- **Q-13** Theorem 14.2 cross-checks against code with R14's bonus l'Hôpital fix: the fallback at the critical inclination `cos i_0 = -1` is `(1/2)·(A_{3,0}/k_2)·sin i_0`. ✓
- **N-11** Line 641 — "the factor 1/8 is half of aycof's 1/4, traced to a second cosine average in the ℓ-kernel" — vague. Pass 4 (BH61 cross-check) must resolve via the long-period generating function.

### §15 Drag Corrections — omgcof, xmcof, delmo, sinmo

- **Q-14** Theorem 15.1 omgcof: matches code. ✓
- **D-8 — CLOSED (2026-06-01, Phase 9; `sgp4_drag_corrections_trace.md` §X.1, `verify_corrections.m`)**
  **Doc-only defect.** The code/SR3 p.12 carry the `ξ⁴` (`coef = qoms4·ξ⁴`, `drag_coefficients.h:121,197`;
  Vallado `propagation.py:1480`); only the §15.4 *box* dropped it. §15.4 box corrected to
  `xmcof = −(2/3)(q₀−s)⁴ξ⁴B*/(e₀η)`. `verify_corrections.m` confirms `xmcof_code = (legacy box)·ξ⁴`.
  (The `−2/3` term-by-term origin is sketch-level — sealed Lane/LH79; value born-digital-confirmed.)
  Original below.

- **D-8 (MAJOR) [original]** Theorem 15.2 boxed form — `xmcof = -(2/3)·(q_0-s)^4·B*/(e_0·η)`.

  **Code line 197:** `dc.M_drag_coef = -ratio<T>(2, 3) * coef * in.bstar / eeta` where `coef = qoms4 · ξ^4` (line 121) and `eeta = e_0·η` (line 117).
  
  So code: `xmcof = -(2/3) · (q-s)^4 · ξ^4 · B* / (e_0·η)`.
  
  Doc: `xmcof = -(2/3) · (q-s)^4 · B* / (e_0·η)`.
  
  **The doc is missing the `ξ^4` factor** present in the code.
  
  **Verdict (Pass 1):** DERIVATION-ERROR or TRANSCRIPTION-ERROR.
  Severity: HIGH (the boxed form is the canonical statement of the theorem; missing ξ^4 means the doc cannot be cited as theoretical basis for code line 197 as written).

### §16-§18 Error-Source Catalog

- **N-12** §16 P-D5 says "Use Kahan summation if `tempa` becomes very small (orbit decay)." Mitigation is recommended but **not implemented** in the code. Verify there is no actual Kahan-summation code (this is a documented gap, not a derivation error). Pass 3.
- **N-13** §17 A-D2 — "Lane power-law atmosphere" baseline 30%; cite [L65]. Pass 2.

### §19-§20 Cross-Reference Index

- **N-14** §19 lists six companion documents under `deprecated/`. Verify each exists and is internally consistent. Pass 2.
- **N-15** §20 lists three BH61 cleanroom files. These are the Pass 4 cross-corpus references. Verify each exists and resolves the relevant secular-rate derivations.

---

## Pass-1 Summary (preliminary; pending Pass 2-5)

| ID | Severity | Section | Class | One-line |
|---|---|---|---|---|
| **D-1** | HIGH | §1 Def 1.1 | DERIVATION-ERROR | B* formula missing ρ_0 factor; units inconsistent. |
| **D-2** | LOW | §4 Eq (4.1) | TRANSCRIPTION | ✅ **CLOSED (2026-06-03).** Intermediate `β·r²/a²` (off by `β²`=O(e²)) corrected to `r²/(a²β)`; final form `β³/(1+e cos f)²` was already correct (`simplify=0`; matches Phase 1 Jacobian / `verify_phase1.m`). |
| **D-3** | HIGH | §4 Eq (4.2) | DERIVATION-ERROR | Doc integrand `β·(r/a)^{m+2}/(1+e cos f)²` ≠ chain result `β³·(r/a)^m/(1+e cos f)²`; differs by f-dependent factor. |
| **D-4** | HIGH | §7 Thm 7.1 | DERIVATION-ERROR | Proof's "absorbing the factor of 2" produces /4 in denominator; boxed form has /1. |
| **D-5** | CRITICAL | §10 Thm 10.1 | DERIVATION-ERROR | Proof's matching step `-4D₂ + 2C₁² = 4C₁² + 16a_0''ξC₁²` gives `D₂ = -(1/2)C₁² - 4a_0''ξC₁²`, not the boxed `D₂ = +4a_0''ξC₁²`. The doc silently drops the (1/2)C₁² term. |
| **D-6** | HIGH | §12 proof | SUFFICIENCY-GAP | Stream-of-consciousness mid-derivation re-think left in text; t₃/t₄/t₅cof not derived at all. |
| **D-7** | MED | §13 line 592 | suspected DERIVATION-ERROR | "β_0^{-3} absorbed into C₁" claim implausible (ψ ≠ β). |
| **D-8** | HIGH | §15 Thm 15.2 | DERIVATION-ERROR or TRANSCRIPTION | Doc xmcof missing ξ^4 factor present in code. |
| **N-1..N-15** | various | various | NOTATIONAL | citations, conventions, deprecated companion docs to verify. |
| **Q-1..Q-14** | various | various | quote checks | most pass; D-4 / D-8 are quote vs proof inconsistencies. |
| **Cascade** | — | §5, §8, §11 | DERIVED | inherit D-3 (via Lane integrals) and D-5 (via D₂). |

**Pre-Pass-2 status:** the existing derivation document has **at least 5 high-or-critical-severity derivation errors** and at least one critical sufficiency gap. The R09 sub-items (§4a, §9a, §11a, §12a) cannot proceed cleanly until these are resolved, because:

- §4a (Lane f†) depends on Lemma 4.2 / Eq (4.2) which carries **D-3**.
- §9a (C₅ 11/4) cascades from the C₂ formula which uses Lane integrals from §4.
- §11a (D₃/D₄) depends on §10 / D₂ which carries **D-5**.
- §12a (t-cofs) extends §12 which carries **D-6** (stream-of-consciousness incomplete proof).

---

