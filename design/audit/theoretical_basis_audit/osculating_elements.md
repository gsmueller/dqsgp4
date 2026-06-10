# Theoretical Basis Audit — `src/orbit/osculating_elements.h`

**File**: `src/orbit/osculating_elements.h` (86 lines)
**Function count**: 1 (`compute_osculating`)
**Scope**: convert the Kepler-equation solution `(E+ω)` plus the eccentricity-vector components `(axN, ayn)` and the semi-major axis `a` into the SGP4 osculating set `{r, ṙ, rḟ, sin u, cos u, u, sin 2u, cos 2u, e², p, β}`.

---

## 1 `compute_osculating(E_plus_w, axN, ayn, a, ke)`

```
=== FORMULA AUDIT CARD ===
ID:                     osculating_elements::compute_osculating
Location:               src/orbit/osculating_elements.h:47-83
Mathematical statement: Given (E+ω, axN, ayn, a, kₑ), compute
                          e·cosE = axN·cos(E+ω) + ayn·sin(E+ω)
                          e·sinE = axN·sin(E+ω) − ayn·cos(E+ω)
                          e²     = axN² + ayn²
                          p      = a(1 − e²)
                          r      = a(1 − e·cosE)
                          ṙ      = kₑ·√a · e·sinE / r
                          rḟ     = kₑ·√p / r
                          β      = √(1 − e²)
                          τ      = 1 / (1 + β)
                          cos u  = (a/r)·(cos(E+ω) − axN + ayn·e·sinE·τ)
                          sin u  = (a/r)·(sin(E+ω) − ayn − axN·e·sinE·τ)
                          u      = atan2(sin u, cos u)
                          sin 2u = 2·sin u·cos u
                          cos 2u = 2·cos²u − 1

THEORY
  Underlying theorem:   Closed-form two-body conic identities applied to
                        the SGP4 modified eccentric-anomaly variable
                        (E+ω) where (axN, ayn) = (e·cos ω, e·sin ω + Δ_LP).
                        Each output is an algebraic identity, not an
                        approximation:
                          (a) e·cosE / e·sinE: rotation of (axN, ayn) by
                              (E+ω), i.e. evaluating the eccentricity
                              vector in the perifocal frame. The identity
                              follows from
                                e·cos(E+ω−ω) = e·cosE,
                                e·sin(E+ω−ω) = e·sinE,
                              expanded with angle-sum formulas.
                          (b) r = a(1 − e·cosE): Kepler's ellipse-radius
                              identity (Vallado Eq. 2-79 or any conic-
                              section text).
                          (c) ṙ = kₑ·√a · e·sinE / r: time-derivative of
                              r via Ė = n/(1−e·cosE) and r = a(1−e·cosE),
                              with kₑ absorbing the n·a scaling
                              (Vallado Eq. 2-87, SR3 p. 13).
                          (d) rḟ = kₑ·√p / r: angular-momentum integral
                              h = √(μp) constant, with rḟ = h/r and
                              kₑ·√p replacing √(μp) in SGP4 units
                              (Vallado Eq. 2-87, SR3 p. 13).
                          (e) cos u, sin u via the f = E + 2 arctan((β·sin E)/
                              (1−β·cos E)) half-angle identity, expanded
                              with τ = 1/(1+β) to avoid singularity at
                              e=0 and rewritten in (axN, ayn) coordinates
                              (SR3 p. 13-14, exactly the SR3 form).
                          (f) sin 2u, cos 2u: standard double-angle
                              identities sin 2u = 2 sin u cos u,
                              cos 2u = 2 cos²u − 1.
                        u itself is recovered by atan2 of the (sin u, cos u)
                        pair — the canonical 4-quadrant inverse.
  Primary reference:    Hoots & Roehrich (1980), Spacetrack Report No. 3,
                        p. 13-14 (equations for U, sin 2u, cos 2u, r,
                        ṙ, rḟ in the SGP4 algorithm — this code is a
                        direct transcription of those equations).
                        Cross-reference: Vallado, "Fundamentals of
                        Astrodynamics and Applications" (4th ed.) §8.6.1
                        for the SGP4 algorithm; §2.2 (Eqs. 2-79, 2-87)
                        for the underlying two-body identities.
  Domain of validity:   e < 1 (elliptical orbit; β = √(1−e²) real);
                        r > 0 (which holds for e < 1 since 1 − e·cosE ≥
                        1 − e > 0); a > 0; kₑ > 0. The variable
                        (E+ω) is in the range that the Kepler solver
                        produced.
                        The τ = 1/(1+β) form avoids the e → 0 limit
                        cleanly (β → 1, τ → 1/2), so the formula is
                        valid throughout 0 ≤ e < 1.

METHOD
  Method declared:      Closed-form algebraic composition. No series,
                        no iteration, no approximation. Each output is
                        the value of the cited identity evaluated on
                        the input TrackedValues using the standard
                        TrackedValue operators (+ − × ÷ √ sin cos atan2).
  Method implemented:   Line-by-line transcription of the identities
                        above. The 13 lines in the function body
                        (osculating_elements.h:59-80) map 1:1 to the
                        13 closed-form expressions in the theory block.
  Match verdict:        ✓ matched — closed-form to closed-form.

ERROR BOUND
  Bound category:       precision (and accuracy, propagated)
  Bound formula:        Every operator used here (`+`, `−`, `×`, `÷`,
                        `√`, `sin`, `cos`, `atan2`, `exact<T>`) is a
                        TrackedValue operation with its own audited
                        closed-form error rule (REQ-EF-3): the bound
                        on the output is the closed-form propagation
                        of the bounds on the inputs through the
                        operator's derivative (mean-value form). No
                        truncation bound is needed because no series
                        is truncated and no iteration is performed.
                        Formally: for any composite F(x₁, …, xₙ) built
                        from these operators on tracked inputs,
                          err(F) ≤ Σⱼ |∂F/∂xⱼ|·err(xⱼ) + roundoff(F),
                        with each Σ term assembled automatically by
                        the operator overloads. This is REQ-EF-3's
                        triangle/mean-value propagation rule.
  Bound implemented:    The function calls only TrackedValue operators
                        and does not directly write to `errors.*`.
                        Therefore the error bound of every output of
                        `compute_osculating` is whatever the underlying
                        TrackedValue operators add — which is the
                        REQ-EF-3 closed-form propagation. Soundness of
                        this card depends on soundness of the
                        TrackedValue operators (see
                        `theoretical_basis_audit/tracked_value.md`).
  Bound verdict:        ✓ matched, subject to TrackedValue's bounds
                        being sound (which is the subject of a separate
                        card). No truncation or iteration error is
                        introduced here. ⚠ The atan2 step is the one
                        non-trivial operator: its bound must be the
                        mean-value bound for atan2 (which is bounded
                        by 1/min(sin²u + cos²u terms) ≈ 1 since
                        (sin u, cos u) live on the unit circle).
                        Verified at TrackedValue level.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation per operator).
                        No REQ-EF-5 (no iteration here — Kepler iteration
                        happens upstream in kepler::solve_kepler).
                        No REQ-EF-6 (no Taylor truncation).
  AUD-EF applies:       AUD-EF-3 (closed-form error propagation
                        through composite operations).
  AUD-MC applies:       n/a — this is a coordinate / quantity extraction,
                        not an algebra operation with its own axioms.
                        Its outputs feed the short-period correction
                        chain (e_sq, sin_2u, cos_2u, r, etc.).
  Verification test:    tests/test_orbit/ — round-trip check
                        a(1 − e·cosE) ↔ r; check that (sin u, cos u)
                        from this function reproduce true anomaly via
                        the f = u − ω identity to within the propagated
                        bound; check sin²u + cos²u = 1 to within bound.

NOTES
  - The formula for (cos u, sin u) is the e=0-safe form: at e = 0,
    axN = ayn = 0 and esine = 0, so cos u = cos(E+ω), sin u = sin(E+ω),
    which is the correct circular-orbit limit.
  - The intermediate `temp_kep = 1 / (1 + β)` is not stored as an
    output but flows into (cos_u, sin_u). Its error bound is the
    REQ-EF-3 closed-form bound on 1/(1+β), which is well-conditioned
    for β ∈ [0, 1].
  - `kₑ` is declared in the docstring as √(GM/aₑ³)·60 with units
    [rad/min]; the user is responsible for supplying it in those
    units. This is a units convention, not a theoretical-basis claim.
  - The function does not check e < 1; that precondition is enforced
    upstream (Kepler solver would not converge for e ≥ 1). The
    `sqrt(1 − e_sq)` for β would NaN if e ≥ 1; this is the
    canonical "garbage in, garbage out" signal and not a soundness
    failure of this card.
  - The double-angle outputs sin 2u, cos 2u could equivalently be
    computed via sin(2·u), cos(2·u) using the TrackedValue trig
    operators on `u`. The chosen form (2·sin u·cos u,
    2·cos²u − 1) avoids a second atan2 + sin/cos round-trip and
    is the exact form used in SR3. Both are closed-form; the
    chosen form has fewer operators and therefore a tighter
    propagated bound.
```

---

## 2 File-level verdict

- **A. Error wiring**: ✓ no direct `errors.*` writes; all bounds flow through TrackedValue operators (REQ-EF-3 path).
- **B. Algebra axioms**: n/a — coordinate / quantity extraction, no algebraic structure of its own.
- **C. Theoretical basis**: ✓ every line is a closed-form two-body identity from SR3 p. 13-14; method is line-by-line transcription; bound is REQ-EF-3 propagation.

**File verdict: PASS**, conditional on TrackedValue operator bounds being sound (audited separately in `theoretical_basis_audit/tracked_value.md`). No method-theory mismatches; no truncation; no iteration; no approximation introduced inside this function.
