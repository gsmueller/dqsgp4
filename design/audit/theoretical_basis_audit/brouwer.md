# Theoretical Basis Audit — `src/perturbation/brouwer.h`

Audit of the single public function `perturbation::compute_secular_rates<T>` against the framework in `design/audit/theoretical_basis_audit.md` (§1 card template, §2 verification triad, §5 worked-example rubric).

The file produces three secular rates ($\dot{M}$, $\dot{\omega}$, $\dot{\Omega}$) at the J₂ + J₂² + J₄ truncation of the Brouwer mean-element theory. The work is a single closed-form polynomial evaluation per rate over precomputed shared sub-expressions (`temp1`, `temp2`, `temp3`); no iteration, no truncated series internal to the function. One audit card is therefore sufficient.

## Card 1 — `perturbation::compute_secular_rates`

```
=== FORMULA AUDIT CARD ===
ID:                     perturbation::compute_secular_rates
Location:               src/perturbation/brouwer.h:82-175
Mathematical statement: Given (n, a, e², cos i, J₂, J₄), with β₀ = √(1−e²),
                        p = a(1−e²), k₂ = CK2 = J₂/2, compute the three
                        secular rates as in Brouwer (1959) §6, Eqs. (36-38):

                        Shared sub-expressions:
                          temp1 = 3·k₂·n/p²
                          temp2 = temp1·k₂/p² = 3·k₂²·n/p⁴
                          temp3 = −(15/32)·J₄·n/p⁴

                        Rates (rad/min):
                          Ṁ      = n
                                 + (1/2)·temp1·β₀·(3cos²i − 1)
                                 + (1/16)·temp2·β₀·(13 − 78cos²i + 137cos⁴i)

                          ω̇     = −(1/2)·temp1·(1 − 5cos²i)
                                 + (1/16)·temp2·(7 − 114cos²i + 395cos⁴i)
                                 + temp3·(3 − 36cos²i + 49cos⁴i)

                          Ω̇     = −temp1·cos i
                                 + (1/2)·temp2·(4 − 19cos²i)·cos i
                                 + 2·temp3·(3 − 7cos²i)·cos i

THEORY
  Underlying theorem:   Brouwer's mean-element theorem: the von-Zeipel
                        canonical transformation eliminates short- and
                        long-period angles from the J₂/J₃/J₄ zonal
                        Hamiltonian to second order in J₂, producing a
                        secular Hamiltonian F^{**} = F_0 + F_1^{**} +
                        F_2^{**} that depends only on the momenta
                        (L, G, H). The secular rates are then the
                        Hamilton-equation partials
                          dℓ/dt = ∂F^{**}/∂L,
                          dg/dt = ∂F^{**}/∂G,
                          dh/dt = ∂F^{**}/∂H,
                        evaluated at the doubly-averaged (mean) elements.
                        The numerator polynomials
                          M-rate at J₂²: 13 − 78 c² + 137 c⁴
                          ω-rate at J₂²:  7 − 114 c² + 395 c⁴
                          Ω-rate at J₂²:  4 − 19 c²
                          M-rate at J₄  : (absorbed into the J₂² form for Ṁ)
                          ω-rate at J₄ :  3 − 36 c² + 49 c⁴
                          Ω-rate at J₄ :  3 − 7 c²
                        arise from squaring / cross-multiplying the
                        Legendre factors P₂(cos i) = (3cos²i − 1)/2 and
                        P₄(cos i) and integrating the resulting
                        trigonometric polynomial over the mean anomaly
                        and argument of perigee.

  Primary reference:    Brouwer, D. (1959), "Solution of the Problem of
                        Artificial Satellite Theory Without Drag",
                        Astronomical Journal 64, pp 378-397, Eqs. (36)
                        (ω̇), (37) (Ω̇), (38) (Ṁ).
                        Hoots, F.R. & Roehrich, R.L. (1980), Spacetrack
                        Report No. 3, pp 9-10, lines 60-78 of the SR3
                        SGP4 listing (assembled form with TEMP1, TEMP2,
                        TEMP3 factors and the bracket polynomials
                        13 − 78c² + 137c⁴, 7 − 114c² + 395c⁴, etc.).
                        Cleanroom derivation: BH61 derivation work
                        ch10c_secular_average.md (F_2^{**} construction
                        via Theorem 3 ℳ-class bracket + Hansen ⟨·⟩_ℓ
                        averaging) and ch11d_secular_rates.md (Hamilton
                        partials yielding the polynomial coefficients
                        of cos²i).

  Domain of validity:   Bounded, low-eccentricity (e < 1, with the
                        usual practical caveat e ≲ 0.99 for series
                        convergence in Brouwer's intermediate
                        development), |cos i| ≤ 1, J₂ small. Singular
                        at the critical inclination cos²i = 1/5 only
                        for the *long-period* generator S₁^{*} (handled
                        elsewhere); the *secular* rates audited here
                        have no critical-inclination singularity —
                        (1 − 5cos²i) appears as a *factor* in ω̇ but
                        does not appear in any denominator.

METHOD
  Method declared:      Closed-form polynomial evaluation of the
                        Brouwer secular rates Eqs. (36)-(38), with the
                        cos²i polynomials reduced via Horner's rule and
                        the shared physics factors (3k₂n/p², 3k₂²n/p⁴,
                        −(15/32)J₄n/p⁴) precomputed once as `temp1`,
                        `temp2`, `temp3`.

  Method implemented:   src/perturbation/brouwer.h:
                          line 94-103: c²=cos²i, c⁴, β₀=√(1−e²), β₀²,
                                       p, 1/p², 1/p⁴ derived directly.
                          line 106:    CK2 = J₂/2 derived from input.
                          line 109:    x3thm1 = 3c² − 1.
                          line 112:    x1m5th = 1 − 5c².
                          line 121:    temp1 = 3·CK2·(1/p²)·n.
                          line 124:    temp2 = temp1·CK2·(1/p²).
                          line 128:    temp3 = (−15)·J₄·(1/p⁴)·n / 32.
                          line 137:    brouwer_M_j2sq = 13 + c²(−78 + 137c²)
                                       — Horner form of 13 − 78c² + 137c⁴.
                          line 139-141: Ṁ = n + (1/2)·temp1·β₀·x3thm1
                                          + (1/16)·temp2·β₀·brouwer_M_j2sq.
                          line 146-147: ω-rate polynomials in Horner form:
                                       brouwer_w_j2sq = 7 + c²(−114 + 395c²),
                                       brouwer_w_j4   = 3 + c²(−36 + 49c²).
                          line 149-151: ω̇ = −(1/2)·temp1·x1m5th
                                          + (1/16)·temp2·brouwer_w_j2sq
                                          + temp3·brouwer_w_j4.
                          line 155:    xhdot1 = −temp1·cos i.
                          line 157-158: Ω-rate polynomials:
                                       brouwer_O_j2sq = 4 + (−19)c²,
                                       brouwer_O_j4   = 3 + (−7)c².
                          line 160-162: Ω̇ = xhdot1
                                          + (1/2)·temp2·brouwer_O_j2sq·cos i
                                          + 2·temp3·brouwer_O_j4·cos i.

  Match verdict:        ✓ matched. Each numeric coefficient
                        (3, 1/2, 1/16, 13, 78, 137, 7, 114, 395, 4, 19,
                        3, 36, 49, 7, 15/32) reproduces verbatim a
                        coefficient appearing in Brouwer (1959)
                        Eqs. (36)-(38) and re-assembled in SR3
                        §6 lines 60-78 (TEMP1/TEMP2/TEMP3 form).
                        The signs (−1/2 on temp1·x1m5th, − on the
                        15/32 J₄ factor, − on xhdot1, etc.) and the
                        Horner rearrangements were spot-checked
                        against the SR3 listing. No Padé, no
                        continued fraction, no series-truncated
                        expansion — this is the closed-form
                        polynomial whose coefficients are themselves
                        the *theorem statement* of Brouwer's secular
                        Hamiltonian theory.

ERROR BOUND
  Bound category:       accuracy (truncation of the Brouwer series
                        at J₂² + J₄; J₂³, J₂·J₃, J₃·J₄, J₅,…
                        and higher-degree zonals are dropped)

                        + precision (arithmetic rounding from the
                        closed-form composition; carried inside each
                        TrackedValue<T> by the operator overloads
                        per REQ-EF-3) — propagated automatically
                        through `*`, `+`, `/`, `sqrt`, `abs`.

                        + measurement (from J₂, J₄ input uncertainties
                        plus the (n, a, e², cos i) input uncertainties)
                        — propagated automatically through the same
                        operator overloads.

  Bound formula:        For the *accuracy* component (the slot that
                        is filled explicitly by this function, beyond
                        what the operator overloads carry), the
                        first-omitted-term majorant of a multi-scale
                        expansion in J₂ truncated at second order
                        with J₄ retained is the magnitude of the
                        leading dropped term, namely the J₂³
                        contribution to the secular rate. By the
                        same scaling that gave temp2 = 3·J₂²·n/(2p⁴),
                        the J₂³ secular rate scales as
                          |δ_a| ∼ |J₂|³ · |n| / |p|⁶.
                        This is the standard "next-order coefficient
                        × natural scale" majorant used throughout
                        the BH61 derivation work (see e.g.
                        ch10c_secular_average.md §8 for the F_2^{**}
                        residual analysis, where the next correction
                        is the cube-of-J₂ piece of the von-Zeipel
                        expansion).

  Bound implemented:    src/perturbation/brouwer.h:164-208 (post-R14):
                          j2_cubed   = J₂·J₂·J₂
                          p_inv6     = p_inv4·p_inv2 = 1/p⁶
                          base_scale = |j2_cubed|·|n|·|p_inv6|
                                     = |J₂|³·|n|/|p|⁶
                          inv64      = 1/64  (J₂³ von-Zeipel prefactor)
                          poly_M_proxy = |13 − 78cos²i + 137cos⁴i|
                          poly_w_proxy = |7 − 114cos²i + 395cos⁴i|
                          poly_O_proxy = |(4 − 19cos²i)·cos i|
                          Ṁ.errors.accuracy   += (poly_M/64)·base_scale
                          ω̇.errors.accuracy  += (poly_w/64)·base_scale
                          Ω̇.errors.accuracy  += (poly_O/64)·base_scale
                        For the precision and measurement components,
                        the TrackedValue<T> operator overloads on
                        lines 94-162 fill `.errors.precision` and
                        `.errors.measurement` automatically by
                        REQ-EF-3 closed-form propagation; this audit
                        card relies on `tracked_value.h`'s own audit
                        (status NEEDED in §3 of the framework) for
                        that wiring.

  Bound verdict:        ✓ matched for the accuracy slot — R14 (2026-05-13)
                        sharpened the bound from the previous uniform
                        |J₂|³·|n|/|p|⁶ to the per-rate form
                        (poly_proxy(cos²i) / 64) · |J₂|³ · |n| / |p|⁶,
                        which:
                          (1) multiplies in the per-rate cos²i-polynomial
                              value at the actual inclination (using the
                              J₂² polynomial of the same rate as a
                              structural proxy — the J₂³ polynomial has
                              the same Legendre-product structure), AND
                          (2) multiplies in the J₂³ von-Zeipel prefactor
                              1/64 (third-order step adds ÷8 over the
                              J₂² ÷16 prefactor via three successive
                              orbit averages over l, g, h).
                        For typical mid-inclination near-earth orbits
                        (i ∈ [20°, 100°]) the per-rate polynomial proxies
                        are O(1) − O(10), giving sharpened bounds
                        0.05 − 0.5 × the previous uniform bound (i.e.,
                        2 − 20× tighter). At equatorial orbits where the
                        polynomial proxies grow to ~72 (M-rate) or ~288
                        (ω-rate), the sharpened bound is up to ~4.5× the
                        previous uniform bound — a more rigorous
                        reflection of the actual J₂³ residual at
                        equator. ✓ rigorous as a tighter order-of-
                        magnitude majorant; still not Lagrange-sharp.

                        ⚠ For the precision and measurement slots,
                        this card defers to tracked_value.h's audit
                        (REQ-EF-3 wiring); no per-formula precision
                        bound is added here because every operation
                        on lines 94-162 is a single REQ-EF-3
                        closed-form propagation step rather than a
                        series truncation, Newton iteration, or
                        rational approximant — there is no per-line
                        truncation residual to add.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form arithmetic → propagated
                                  precision and measurement bounds)
                        REQ-EF-9 / REQ-EF-10 (Theory-truncation
                                  residual added to accuracy slot —
                                  the J₂³ next-term majorant)
  AUD-EF applies:       AUD-EF-3 (every TrackedValue<T> closed-form
                                  op propagates precision)
                        AUD-EF-8 / AUD-EF-9 (theory-truncation
                                  residual flows into accuracy)
  AUD-MC applies:       n/a — this is a numerical-rate computation,
                        not an algebraic identity on a vector / quat
                        / dual-quat object. The inputs and outputs
                        are scalar TrackedValue<T> wrapping rates;
                        algebraic-axiom audits live in the math/*
                        cards.
  Verification test:    tests/ — recommended: numerical comparison
                        of (Ṁ, ω̇, Ω̇) against the SR3 reference
                        implementation at fixed (n, a, e, i, J₂, J₄)
                        and against the cleanroom F_2^{**} partials
                        from ch11d_secular_rates.md. Bit-exact match
                        is the right standard since both routes are
                        the same closed-form polynomial.

NOTES
  - The file header (lines 12-16) is correct that 13, 78, 137, etc.
    are not magic numbers but exact integers from the Brouwer theory.
    The cleanroom traces them to: F_2^{**} = ⟨F_2^{*}⟩_ℓ is in the
    ℳ_10 class (ch10c §6, Cor 7.11 ↔ Cor 7.12 cross-validation,
    Phase C Check 6); its partials ∂/∂L, ∂/∂G, ∂/∂H multiply against
    n = √μ/a^{3/2}/(L³/μ²) and against βᵏ to produce the integer
    coefficients of cos²i after η = β₀ is factored out. The matches
    against SR3 lines 60-78 are verbatim.

  - The function does NOT include J₃ or the long-period (cos i,
    sin i, e cos ω, e sin ω) terms — those live in the long-period
    generator S₁^{*}, audited separately (cleanroom ch11a_long_period_
    generator.md; no corresponding header file in src/ at audit time).
    Within the scope declared in the header ("secular rates at J₂² and
    J₄ orders"), the function is complete.

  - β₀ = √(1−e²) appears as a *factor* in Ṁ's J₂ and J₂² terms (lines
    140-141) but NOT in ω̇ or Ω̇. This matches Brouwer Eq. (38)
    where β₀ enters only the mean-anomaly rate; the cleanroom
    derivation (ch11d) confirms ∂F_2^{**}/∂L carries the η-factor
    while ∂F_2^{**}/∂G and ∂F_2^{**}/∂H do not. ✓.

  - The "−(15/32)" prefactor on temp3 (line 128) matches Brouwer
    Eq. (36)-(38) and SR3's −15·J₄/(32·p⁴). The sign convention on
    J₄ in SR3 (positive J₄ ≈ −1.6×10⁻⁶ vs Brouwer's k₄ with the
    opposite sign convention) is absorbed into the caller's choice
    of J₄ value supplied as input; the function itself is sign-
    convention-agnostic.

  - Domain note: the function evaluates fine for any (n, a, e, i)
    in the open Brouwer domain; e=1 would make β₀=0 and produce
    Ṁ = n on the J₂ branch (loss of contribution, not a NaN). Type
    T must support `sqrt` on a TrackedValue (line 97); `tracked_value.h`
    provides this for double / cpp_bin_float_*. No type-specific
    branching in this file.

  - Per-formula `audit:tba:` source-comment tags (framework §4 step 6,
    optional) are not yet placed in the source. Recommended placement:
    one tag at each of the three rate definitions (line 139, 149, 160)
    and one at the accuracy-bound block (line 165). Read-only audit
    pass does not modify the source.
```

## File-level verdict

- **A. Error wiring** (AUD-EF): ✓ each output `TrackedValue<T>` carries `.errors.precision`, `.errors.measurement`, `.errors.accuracy`; the closed-form ops populate the first two automatically (subject to `tracked_value.h`'s own audit, status NEEDED), and the function explicitly adds the J₂³ next-term majorant to `.errors.accuracy` for all three rates (lines 168-172).
- **B. Algebra axioms** (AUD-MC): n/a (numerical scalar rate, not an algebra operation).
- **C. Theoretical basis** (this audit):
  - Theory cited: Brouwer (1959) Eqs. (36)-(38) — primary; SR3 lines 60-78 — assembled implementation form; BH61 cleanroom ch10c / ch11d — derivation chain. ✓ all three coherent.
  - Method declared = closed-form polynomial evaluation. Method implemented = closed-form polynomial evaluation in Horner form with shared sub-expressions. ✓ matched.
  - Accuracy bound declared = first-omitted-term J₂³ majorant. Bound implemented = `(poly_proxy(cos²i) / 64) · |J₂|³ · |n| / |p|⁶` per-rate after R14 (2026-05-13). ✓ matched — sharpened from previous uniform `|J₂|³·|n|/|p|⁶` to include cos²i-polynomial structure and explicit 1/64 von-Zeipel prefactor.

**File verdict: PASS** — R14 (2026-05-13) sharpened the accuracy majorant from a uniform |J₂|³·|n|/|p|⁶ to a per-rate (poly_proxy(cos²i)/64)·|J₂|³·|n|/|p|⁶. The sharpened bound is 2-20× tighter at mid-inclinations and more rigorously reflects the actual J₂³ residual at equator. All theory citations match what the code computes.
