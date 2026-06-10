## PASS 1 — Internal consistency (in progress)

Every line of every proof in §0 … §20 read. Findings collected below in declaration order.
A finding ID convention:

- **D-N** = derivation finding (the algebra is wrong or the chain doesn't close)
- **N-N** = notational / convention finding (ambiguous but not wrong)
- **Q-N** = quote / boxed-form-vs-proof finding (proof and boxed result disagree)
- **T-N** = transcription / formatting

Severity will be revisited after Passes 2-3 because some findings turn out to be APPLICATION-ERROR (the math is right per source but is then mis-applied) rather than DERIVATION-ERROR.

---

### §0 Notation

- **N-1** Line 19 — `[LH79]` cited as "Project Space Track Report **No. 2**". Conventional citation is *Spacetrack Report No. 2* but spelling and number need cross-check against the actual NORAD/AFSPC publication record. Pass 2.
- **N-2** Line 38 — `k₄ := -(3/8) J₄ a_E⁴`. Sign of the (3/8) prefactor varies across sources (Vallado uses `+(3/8) J_4`, Brouwer uses different splitting). Verify against [B59] §1 or [SR3] p. 9. Pass 2.
- **N-3** Line 46 — `ψ² := |1 - η²|`. The absolute value is operational; for bound orbits η < 1 so 1-η² > 0 always. Flag as **N-2/notational** but harmless.
- **N-4** Symbol table omits `f†` (the Lane fictitious anomaly that appears in §4 Eq (4.2), Lemma 4.2, Note 4.3). Should be added.

### §1 Drag Acceleration

- **D-1 — CLOSED (2026-06-01, Phase 9)** Definition 1.1 corrected to
  `B* = C_D ρ₀ A/(2m)` (the `ρ₀` was missing, making it dimensionally `[L²/mass]` not `ER⁻¹`).
  The Lane density now enters as the dimensionless profile `ρ/ρ₀ = ((q₀−s)/(r−s))⁴`, `ρ₀` riding
  in `B*` — the convention the code and the code-matched Phase 2.A §A.2 use. Original below.

- **D-1 (MAJOR) [original]** Lines 65-68, Definition 1.1.

  Stated: `B* := C_D · A / (2 m)` with units `ER⁻¹`.
  
  Dimensional analysis: C_D is dimensionless, A has units [length²], m has [mass]. So A/m has [length²/mass] = [m²/kg] in SI, which is **not** [1/length]. The formula as written is dimensionally **incompatible** with the stated ER⁻¹ unit.
  
  Standard SGP4 / SR3 definition: `B* = (C_D · ρ_0 · A) / (2 m)` with ρ_0 in [mass/volume]. Now (mass/volume)·(length²/mass) = 1/length = ER⁻¹. ✓
  
  **Verdict (Pass 1):** Definition 1.1 is missing the `ρ_0` factor. This is a **DERIVATION-ERROR** unless ρ_0 has been silently absorbed elsewhere upstream.
  
  Severity: HIGH. The downstream C₁..C₅ derivations all implicitly assume B* carries the ρ_0 factor (e.g. they multiply `B* · ρ(r)` only with the spatial profile of ρ, not its scale). If ρ_0 is missing from B*, the C-coefficients are off by ρ_0.

### §2 Lane Power-Law Atmosphere

- **N-5** Line 102 — `ρ_0 = 2.461 × 10⁻⁵ XKMPER⁻¹`. The unit XKMPER⁻¹ = 1/km is not a density unit. The Lane-Hoots literature reports the SGP4 reference density in (kg / earth-radius³) or in a dimensionless-rescaled form; the doc's choice of XKMPER⁻¹ is consistent with the SGP4 source-code convention but is misleading without an explanatory note.
- **N-6** Lines 105-106, Lemma 2.2 — claims the τ=4 power-law "matches a 4-parameter exponential atmosphere exactly at two altitudes and to second order in (r-rₐ) at one of them". This claim is not proved; the four parameters of the exponential are ρ_a, r_a, H, and the cutoff; the four parameters of the Lane form are ρ_0, q_0, s, τ. A matching argument would equate four conditions: two ρ values, one ρ′, one ρ″ at chosen points. Pass 2 must verify this against [L65].
- **Q-1** Line 116, Case C — `s* = 1 + 20/XKMPER`. Cross-check against `[SR3]` p. 11-12 (the actual fit-parameter case-table). Pass 2.

### §3 J₂-Perturbed Radial Distance

- **D-11 — CLOSED (2026-06-01, Phase 9; found by Phase 2.A §A.6)** §3 Theorem 3.1's
  `⟨δr⟩_{J2} = −(k₂/p)(3θ²−1)/2` is the *short-period amplitude* — a factor ~3 too small for the
  *secular radial multiplier* `⟨δr⟩/r = −(3/2)(k₂/p²)β(3θ²−1)` that the drag density coupling
  requires (and that C₂ Part B / C₄ Part B use). Theorem 3.1 corrected to the secular form
  (born-digital: Vallado `SGP4.cc:618`, SR3 p.14); the §3.2–§5 cascade is flagged **superseded** by
  the code-matched per-coefficient trace docs (`sgp4_drag_phase2a_C2_trace.md`, `sgp4_drag_C4_trace.md`).

- **N-7** Line 129, Theorem 3.1 — cites `[B59 §3.1, Eq.(14)]` for the secular radial perturbation `⟨δr⟩ = -(k₂/p)·(3θ²-1)/2`. The Brouwer 1959 paper has the radial perturbation distributed across the short-period and long-period generators; identifying which equation in B59 produces this **secular** radial is not obvious. Pass 2 must verify [B59 Eq.(14)] is the correct citation. Memory bank notes: `feedback_born_digital_latex.md` prefers post-LaTeX digital sources; check `reference_brouwer_1959.md` for the equation map.
- **Q-2** Algebra of Eq (3.2): re-derivation matches the doc's claim `δρ/ρ = 2 k₂ ξ (3θ²-1) / p` for τ=4. ✓ (algebra closes)

### §4 Orbit Averaging — General Formulas

- **D-9 — CLOSED (2026-06-01 by Phase 2.A; see `design/derivations/sgp4_drag_phase2a_C2_trace.md` §A.4.2)**
  Phase 2.A verified symbolically that the correct `I^{(0,4)} = (2+3η²)/(2(1−η²)^{7/2})`
  (Phase 1 Thm 1.3.2.4) is **exactly** the SGP4 C₂ Part-A constant term
  (`verify_phase2a.m`, `simplify=0`). The main-doc Eq (4.4) value below is a
  **transcription error**. Original PRELIMINARY analysis retained for the record:

  Line 190, Eq (4.4) — `I^{(0,4)} = (2+η²)/(1−η²)^{7/2} + O(e⁰η⁶)`.

  **Observation.** Under the literal reading
  `I^{(0,4)}(η) := (1/(2π)) · ∫₀^{2π} df / (1 − η cos f)^4`, the recurrence
  `dI^{(0,m)}/dη = (m/η)·[I^{(0,m+1)} − I^{(0,m)}]` applied to
  `I^{(0,1)} = (1−η²)^{−1/2}` produces

      I^{(0,1)} = (1−η²)^{−1/2}
      I^{(0,2)} = (1−η²)^{−3/2}
      I^{(0,3)} = (2 + η²) / (2 (1−η²)^{5/2})
      I^{(0,4)} = (2 + 3 η²) / (2 (1−η²)^{7/2})        ← from Phase 1 Theorem 1.3.2
      I^{(0,5)} = (8 + 24 η² + 3 η⁴) / (8 (1−η²)^{9/2})

  The Phase 1 Theorem 1.3.2.4 value (2+3η²)/(2(1−η²)^{7/2}) gives 1 at η=0
  (matches integrand ≡ 1 ⇒ integral = 1 ✓), whereas the main doc's claimed
  (2+η²)/(1−η²)^{7/2} gives 2 at η=0 (✗). So under the literal reading of
  the definition, the main doc's Eq (4.4) does not match the integral.

  **However, this is NOT yet a closed finding.** Per user directive
  2026-05-15:

  > "We can't accept results that we aren't sure of. We will need to
  > understand how the theorems apply and get the correct result, or it is
  > possible that the theorems have an incorrect assumption."

  Three possibilities are not yet ruled out:

  1. **Theorem 1.3.2 is correct AND it is the integral SGP4 computes.** The
     main doc's Eq (4.4) is a transcription error. Phase 2's clean C₂
     assembly must (a) trace how Theorem 1.3.2.4 / 1.3.2.5 propagate through
     the orbit-averaging to reproduce the SGP4 code's `(8 + 24η² + 3η⁴)`
     polynomial in C₂ Part B. Preliminary observation: the polynomial matches
     the *numerator* of `I^{(0,5)}` from Theorem 1.3.2.5 — but the ψ-power
     bookkeeping and the prefactor algebra need to be checked rigorously.
  2. **Theorem 1.3.2 is correct, but it is not the integral SGP4 computes.**
     The actual SGP4 orbit-average integrand is `(r − s)^{−4} · X(f)` with
     `r − s = (a β² − s − s e cos f) / (1 + e cos f)`, which is NOT of the
     form `(1 − η cos f)^{−4}` without the Lane fictitious-anomaly
     substitution `r − s = (1/ξ)(1 − η cos f†)`. The Lane substitution is
     approximate (Note 4.3 of the main doc: O(eη) corrections). So
     `I^{(0,m)}` here may evaluate the *post-substitution* integral, which
     differs from the *true* SGP4 orbit-average by O(eη) terms.
  3. **Theorem 1.3.2 has a hidden assumption.** The recurrence's
     differentiation-under-the-integral-sign step (1.2.1.2) assumes the
     integrand has continuous η-dependence in a neighbourhood of evaluation.
     This holds for `|η| < 1` away from `η = 1`, but the SGP4 atmospheric
     model has perigee-regime case-table behaviour (Definition 0.6.3) where
     `s` itself depends on the orbit elements via `h_p`. The Theorem 1.3.2
     formulas assume `s` is a constant independent of `(a, e)`, which it is
     within each regime but NOT across regime boundaries. The doc's Eq (4.4)
     may correspond to a different parameterization in which `s` varies
     continuously and the recurrence picks up an additional term.

  **Required next step (Phase 2.A — investigate the doc-code-integral
  connection).** Before D-9 is classified, Phase 2 must:
  - Trace the SGP4 C₂ formula from the orbit-average to the closed form.
  - Identify which Lane integrals appear in the trace.
  - Verify either that Theorem 1.3.2 values produce the SGP4 code's
    `(8 + 24η² + 3η⁴) · ψ^{−7}` form (case 1 above), or that the SGP4 actual
    integral is different from Definition 1.1.1 (case 2 or 3 above).

  Until Phase 2.A is complete, the D-9 finding is **PRELIMINARY** and must
  not be cited as a closed finding. The Phase 1 Theorem 1.3.2 closed forms
  stand as correct evaluations of the literally-defined Lane integrals, but
  the connection to SGP4 is open.

- **D-10 (CRITICAL — found 2026-05-15 by Standard-10-driven Phase 2.A setup; affects Phase 0)**

  **Status: FIXED in Phase 0-rev1 (Block A of the remediation plan, 2026-05-16).**
  Theorem 0.4.1 rewritten with `F_drag = -B*·ρ·|v|·v` (A1); Theorem 0.4.2's 5
  closed-form drag rates recomputed with `(1+e²+2e cos f)^{3/2}` factor (A2);
  cascade through §0.3.1 preview + §0.4 alignment bullet (c)/(d) + §0.5
  Theorem 0.5.3 integrand + §0.7 summary (A3); `verify_phase0.m` extended with
  (0.4.1.10) algebraic identity check, (0.4.1.11) speed-scalar closed-form
  check, 3× Postulate-1.2 match checks (r̂/t̂/n̂, regression guards), and an
  updated Theorem 0.4.2 (0.4.2.1) check using the `(...)^{3/2}` form (A4);
  verifier passes 32/32 (A5). The Postulate-1.2 match check is the symbolic
  guard that would have caught D-10 originally and now prevents recurrence.

  `design/derivations/sgp4_drag_phase0_foundations.md` §0.4 Theorem 0.4.1
  (drag specialization) is **missing the `|v|` factor** required by the
  Newtonian quadratic drag law.

  **The error.** Phase 0 §0.4 (0.4.1) writes

      F_drag := −B* · ρ(r) · v             [WRONG — linear in v]

  This is dimensionally Stokes-like (force ∝ velocity), not the quadratic
  atmospheric drag that the main doc Postulate 1.2 specifies:

      F_drag = −B* · ρ · |v_rel| · v_rel    [CORRECT — quadratic, |v|·v form]

  The consequence is that the Phase 0 closed-form drag rates in Theorem 0.4.2
  (0.4.2.1)..(0.4.2.5) are all off by a `|v|` factor. Specifically:

      ȧ_drag (Phase 0)   = −(2 B*ρ a / β²) · [e² sin² f + (1+e cos f)²]
                         = −(2 B*ρ a / β²) · (1 + e² + 2e cos f)         [WRONG]

      ȧ_drag (correct) = −(2 n B* ρ a² / β³) · (1 + e² + 2e cos f)^{3/2}  [CORRECT]

  The ratio (correct / Phase-0-wrong) = |v| = (na/β)·√(1+e²+2e cos f),
  which is the speed scalar (verified via vis-viva
  v² = μ(2/r − 1/a) = (μ/(aβ²))·(1 + e² + 2e cos f)).

  **Discovery mechanism.** This finding was surfaced by Standard 10's
  code-match requirement BEFORE Phase 2.A was formally executed. While
  preparing the Phase 2.A symbolic trace from Phase 0 + Phase 1 to SGP4 code:

  1. SGP4 C₂ Part A has polynomial `(1 + (3/2)η² + 4eη + eη³)` which (per
     standard analyses of the AFGP4/SGP4 drag derivation) emerges from the
     orbit-averaged density integral with the `|v|³` factor.
  2. From the energy form `dE/dt = F·v = −B*ρ|v|³`, ȧ = (2a²/μ)·dE/dt
     produces a `(1+e²+2e cos f)^{3/2}` integrand factor.
  3. My Phase 0 Theorem 0.4.2 produces `(1+e²+2e cos f)^1` instead.
  4. The factor of `√(1+e²+2e cos f) = |v|/(na/β)` between the two is
     exactly the missing `|v|` from Theorem 0.4.1.

  **Why the verifier did not catch it.** `verify_phase0.m`'s 27/27 PASS
  checked algebraic consistency of the (wrong) Phase 0 Theorem 0.4.1 with
  itself — e.g. "R_drag / T_drag = e sin f / (1+e cos f)" (which holds
  regardless of `|v|` since it cancels in the ratio). No check tested
  against the external reference (Postulate 1.2 quadratic drag law). This
  is exactly the failure mode Standard 10 is designed to expose:
  numerical / algebraic-internal verification is insufficient.

  **Verdict (Pass 1):** DERIVATION-ERROR. CRITICAL. Phase 0 §0.4 must be
  re-derived from the correct quadratic drag form. The verifier needs an
  additional check that the drag-acceleration formula matches the physical
  postulate `F_drag = −B*ρ|v|v`. Phase 2.A cannot proceed on the existing
  Phase 0 §0.4 — the Phase 0 §0.4 fix is the immediate next task.

  **Required fix (Phase 0-rev1):**

  1. Theorem 0.4.1: change `F_drag := −B* · ρ · v` to `F_drag := −B* · ρ · |v| · v`,
     and update R_drag / T_drag to include the `|v|` factor:

         R_drag = −B* · ρ · |v| · ṙ
         T_drag = −B* · ρ · |v| · r ḟ

  2. Theorem 0.4.2 closed form: substitute the corrected R_drag / T_drag
     into the Gauss VE (Theorem 0.3.2):

         ȧ_drag = −(2 n B* ρ a² / β³) · (1 + e² + 2 e cos f)^{3/2}

     and similarly for `ė_drag`, `ω̇_drag`, `Ṁ_drag`. The `(...)^{3/2}` factor
     is the v²·|v| = |v|³ kinematic content of the quadratic drag.

  3. Add to `verify_phase0.m` an explicit symbolic check that the Theorem
     0.4.1 result agrees with the postulated F_drag (vector) under
     `|v| := √(ṙ² + (r ḟ)²) = (na/β)·√(1+e²+2e cos f)`. This is the
     code-match-style check that should have caught D-10.

  4. Re-run verifier and confirm 27+(new)/27+(new) PASS.

  5. Re-evaluate downstream impact:
     - D-9 (preliminary) becomes irrelevant for the immediate question, since
       the orbit-averaged integrand changes from `(1+e²+2e cos f)` to
       `(1+e²+2e cos f)^{3/2}`. The Lane integrals `I^{(p,m)}` of Definition 1.1.1
       still hold in their literal definition, but the actual SGP4 integral
       projects `(1+e²+2e cos f)^{3/2}` against `ρ(r)` — which is a *different*
       combination of Lane integrals than the `(1+e²+2e cos f)^1` case.
       The connection between the SGP4 code's `(8 + 24η² + 3η⁴)·ψ^{−7}` and
       the Lane-integral expansion of `(1+e²+2e cos f)^{3/2} · ρ(r)` is the
       Phase 2.A task — now well-defined.
     - All other Phase 0 alignment remarks for drag specialization
       (§0.3.1.alignment, §0.4 alignment remarks) require revision.

  **D-1 ↔ D-10 relationship.** D-1 (B* definition missing ρ_0 factor) and
  D-10 (drag specialization missing |v|) are independent. D-10 is far more
  serious for downstream Phase work because it changes the algebraic form
  of the orbit-averaged integrand. D-1 is a dimensional-bookkeeping issue
  that doesn't change algebraic structure (only units interpretation).
  Both must be fixed before Phase 2.A can proceed.

- **D-2** Line 164, Eq (4.1) — `dM = (β r²/a²) df = β³/(1+e cos f)² df`.

  The intermediate form `β r²/a²` is **wrong**. The correct Kepler relation is `dM/df = r²/(a²β)` (from `r² ḟ = na²β` ⇒ `dM/df = n/(ḟ) = r²/(a²β)`).
  
  Substituting `r/a = β²/(1+e cos f)`:  `(r²/a²) / β = (β⁴/(1+e cos f)²) / β = β³/(1+e cos f)²`. ✓ (matches final form)
  
  So the intermediate form `β r²/a²` is a TRANSCRIPTION-ERROR; the final closed form is correct, but the intermediate has β in the **numerator** instead of the denominator.
  
  Severity: LOW — the proof line below (line 167) correctly states `dM/df = r²/(a²β)`. The transcription error is just in the displayed equation (4.1).

- **D-3 (MAJOR)** Lines 169-176, Lemma 4.2 / Eq (4.2) — orbit-averaged density × r-power identity.

  Stated: `⟨ρ(r)·(r/a)^m⟩ = ρ_0·(q-s)^4·ξ^4 · (1/2π)·∫₀^{2π} [β·(r/a)^{m+2}] / [(1+e cos f)² · (1-η cos f†)^4] df`.
  
  Re-derivation. Start: `⟨ρ(r)·(r/a)^m⟩_M = (1/2π)·∫ ρ(r)·(r/a)^m dM`. Convert with dM = β³/(1+e cos f)² df (from Eq 4.1 corrected):
  
  `= (1/2π)·∫ ρ(r)·(r/a)^m · β³/(1+e cos f)² df`.
  
  With the Lane substitution `ρ(r) = ρ_0·(q-s)^4·ξ^4 / (1-η cos f†)^4`:
  
  `= ρ_0·(q-s)^4·ξ^4 · (1/2π)·∫ [β³ · (r/a)^m] / [(1+e cos f)² · (1-η cos f†)^4] df`.
  
  Compared to the doc's claim: doc has `β·(r/a)^{m+2}` in the numerator vs. the correct `β³·(r/a)^m`.
  
  Using `(r/a)² = β⁴/(1+e cos f)²`: doc-numerator = `β·(r/a)^m·(r/a)² = β⁵·(r/a)^m/(1+e cos f)²`. Correct numerator = `β³·(r/a)^m`.
  
  Ratio (doc / correct) = `β² / (1+e cos f)²` — **f-dependent**, not a constant. So the doc's Eq (4.2) is **algebraically inconsistent** with the derivation chain.
  
  **Verdict (Pass 1):** DERIVATION-ERROR.
  Severity: HIGH. Lemma 4.2 / Eq (4.2) is the *master* orbit-averaging identity that §5 (C₂), §8 (C₄), §9 (C₅) all invoke. An error here propagates to every primary drag coefficient.

  Hypothesis for Pass 2: Lane-Hoots [LH79 p. 16] may use a different starting form (e.g. averaging over `f` directly rather than over `M`, which changes the Jacobian). Verify which form LH79 uses, and reconcile.

- **N-8** Note 4.3 (line 178) — Lane's `f†` substitution is "exact only at first order in e" with corrections O(e η). This is the open theoretical gap that §4a (R09 sub-item) is supposed to formalize. Pass 1 status: NOTED.

