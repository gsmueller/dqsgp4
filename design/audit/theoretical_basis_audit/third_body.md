# Theoretical Basis Audit — `src/perturbation/third_body.h`

**File**: `src/perturbation/third_body.h` (179 lines)
**Functions audited**: 1 (`compute_third_body_rates`)
**Framework**: `design/audit/theoretical_basis_audit.md` §1, §5

---

## Card 1 — `compute_third_body_rates`

```
=== FORMULA AUDIT CARD ===
ID:                     third_body::compute_third_body_rates
Location:               src/perturbation/third_body.h:70-176
Mathematical statement: Five secular rate corrections (dΩ/dt, dω/dt, dM/dt,
                        de/dt, di/dt) of a satellite's orbital elements
                        induced by a single third body (Sun or Moon),
                        computed from the orbit-averaged P₂(cos S) term of
                        the third-body disturbing function expanded to
                        second order in (r/r₃), via Lagrange's planetary
                        equations.

THEORY
  Underlying theorem:   Lagrange's planetary equations applied to the
                        orbit-averaged disturbing function
                          R = (μ₃/r₃) · (a/r₃)² · ⟨P₂(cos S)⟩_M
                        where S is the geocentric satellite-third-body
                        angle. P₂(cos S) = (3 cos²S − 1)/2 is the second
                        Legendre polynomial; (r/r₃)² is the leading
                        nontrivial term of the (r/r₃)^k Legendre expansion
                        after the constant.
                        The ecliptic→equatorial frame transformation by
                        the obliquity ε is a closed-form rotation about
                        the vernal equinox axis (x-axis).
                        Projection of the third-body unit vector onto the
                        satellite's orbital triad (a₁, a₂, a₃) is the
                        standard Kozai/Brouwer-Clemence decomposition.
  Primary reference:    - Brouwer & Clemence (1961), "Methods of Celestial
                          Mechanics", Chapter XI (lunar/solar perturbations).
                        - Kozai (1959), "On the Effects of the Sun and Moon
                          upon the Motion of a Close Earth Satellite",
                          Smithsonian Astrophys. Obs. Spec. Rep. 22.
                        - Hoots & Roehrich (1980), SPACETRACK Report #3,
                          §B (Lunar/Solar Perturbations — Lyddane-form for
                          SDP4). The five expressions structurally match
                          SR3 Eqs. for (dM/dt)₃, (dΩ/dt)₃, (dω/dt)₃,
                          (de/dt)₃, (di/dt)₃.
                        - File-internal: design/derivations/
                          017_third_body_perturbation.md (referenced;
                          contents not verified in this card).
  Domain of validity:   - e ∈ [0, 1) with β = √(1−e²) > 0 (denominator).
                        - sin(i) ≠ 0 — equator-aligned orbits hit a
                          singularity in dω/dt and di/dt (sin(i) in
                          denominator, line 158, 167).
                        - r/r₃ ≪ 1 (low-Earth satellite vs Sun/Moon
                          distance) — the P₂ truncation drops O((r/r₃)³)
                          and higher; for LEO and r₃ = Sun/Moon, r/r₃ ≲
                          1.7e-5 (Moon) to 4.3e-8 (Sun).

METHOD
  Method declared:      (a) Closed-form ecliptic→equatorial rotation of
                            the third-body unit vector (lines 98-100):
                              x₃ =  cos β cos λ
                              y₃ =  cos β sin λ cos ε − sin β sin ε
                              z₃ =  cos β sin λ sin ε + sin β cos ε
                        (b) Closed-form projection onto orbital triad
                            (lines 114-116):
                              a₁ =  x₃ cos Ω + y₃ sin Ω
                              a₂ = −x₃ sin Ω cos i + y₃ cos Ω cos i
                                  + z₃ sin i
                              a₃·sin i = z₃
                        (c) Closed-form algebraic combinations of
                            (a₁, a₂, a₃·sin i, sin ω, cos ω, e, β, n, a)
                            feeding the five rate formulas (lines 152-173).
                        (d) No iteration, no series truncation, no Padé,
                            no continued fraction. The orbit-average over
                            mean anomaly M has ALREADY been performed
                            analytically upstream — what remains is a
                            single closed-form polynomial expression in
                            the listed primitives.
  Method implemented:   Lines 98-173 are direct closed-form arithmetic
                        on `TrackedValue<T>` operands using +, −, *, /,
                        sin, cos, sqrt. The `exact<T>(k)` and
                        `ratio<T>(p,q)` constructors supply rational
                        constants (1, 3, 1/3, 1/2).
  Match verdict:        ✓ method matches declared theory — closed-form
                        evaluation of the orbit-averaged P₂ result. No
                        Taylor branch, no iteration. The function is a
                        polynomial map of its inputs.
                        ? un-cross-checked: each of the five rate
                          expressions vs SR3 §B Eqs. The structural
                          pattern (a₁ cos ω + a₂ sin ω), (−a₁ sin ω +
                          a₂ cos ω), (a₁² + a₂² − 1/3), the β = √(1−e²)
                          placement, and the cos(i)/sin(i)² factor in
                          dω/dt are all hallmarks of Kozai-Brouwer-Clemence
                          / SR3 form, but a term-by-term symbolic match
                          to a primary source is not performed in this
                          read-only audit.
                        ⚠ derivation file referenced (017_third_body_
                          perturbation.md) was not opened — the chain
                          from R to the five rate expressions is asserted,
                          not verified here.

ERROR BOUND
  Bound category:       precision (closed-form propagation through
                        `TrackedValue<T>` ops via REQ-EF-3); accuracy
                        contribution from the P₂ truncation of the
                        Legendre series and from the orbit-average
                        truncation in M (Kozai expansion).
  Bound formula:        Per REQ-EF-3, each TrackedValue op (+, −, *, /,
                        sin, cos, sqrt) adds its closed-form propagated
                        bound to `result.errors.precision`. No truncation
                        term is added IN this function — the function
                        is purely closed-form.
                        Accuracy bound for the THEORY (not added here):
                          - P₂ truncation drops (r/r₃)³ P₃(cos S) +
                            (r/r₃)⁴ P₄(cos S) + …; magnitude ≲
                            (r/r₃) · |P₂| · (1 + ε); for r/r₃ ≲ 2e-5
                            (Moon) this is a ~2e-5 relative truncation.
                          - Orbit-average truncation: ⟨P₂⟩_M keeps only
                            mean-anomaly-secular terms; short-period
                            terms are dropped (this is the desired
                            secularization, not an error).
                        These accuracy bounds belong to the **caller**
                        (whoever supplies `perturbation_coef` and folds
                        in r₃, μ₃), not to this function.
  Bound implemented:    No explicit `errors.*` mutation in this function
                        — all bounds propagate implicitly through the
                        `TrackedValue<T>` operators.
                        ? un-verified: that every operator used (`*`,
                          `-`, `/`, `sqrt`, `sin`, `cos`) is a
                          `TrackedValue<T>` overload that does the
                          REQ-EF-3 propagation. The file `#includes`
                          only `tracked_value.h`; the use of `sin`,
                          `cos`, `sqrt` is `auto`-deduced ADL — if
                          these resolve to `std::sin` (taking `T`, not
                          `TrackedValue<T>`), the error wiring is
                          silently broken. Cross-check needed against
                          `tracked_value.h` (which must define ADL
                          overloads for `sin`, `cos`, `sqrt`).
  Bound verdict:        ✓ for the formula content: closed-form
                          arithmetic with implicit per-op REQ-EF-3
                          propagation is the correct bound family for
                          a polynomial map of TrackedValue inputs.
                        ? accuracy contribution of the P₂ truncation
                          and orbit-average truncation is not added
                          here — must be added at the caller (the file
                          that computes `perturbation_coef` and chooses
                          to use only the P₂ term).
                        ⚠ wiring of `sin`, `cos`, `sqrt` to
                          TrackedValue overloads is not confirmed in
                          this file — defer to `tracked_value.h` audit
                          card (AUD-EF-1 / REQ-EF-3 verification).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form per-op error propagation
                        for *, /, +, −, sin, cos, sqrt).
  AUD-EF applies:       AUD-EF-1 (every public op returns
                        TrackedValue<T>), AUD-EF-3 (closed-form rule
                        applied per op), AUD-EF-10 (no bare `T` leak).
                        ? whether the rate fields (dM_dt, …) being
                          struct members satisfies the "returns
                          TrackedValue<T>" requirement — `ThirdBodyRates`
                          is a struct of five TrackedValues. ✓ each
                          field is `TrackedValue<T>`.
  AUD-MC applies:       n/a — this is a perturbation formula, not an
                        algebra operation.
  Verification test:    ? no test file referenced in the framework's
                          §3 table for `src/perturbation/`. The file is
                          NOT in the §3 audit-status table (24 files
                          listed, none under `src/perturbation/`). This
                          card is being produced as an extension /
                          out-of-scope adjunct to the §3 sweep.
                        Recommended test: compare numerical output for
                        a fixed (n, a, e, i, ω, Ω, body angles) tuple
                        against an independent re-implementation of SR3
                        §B Lyddane-form Sun/Moon rates, or against
                        Vallado's `dpper` / `dscom` reference values.

NOTES
  - The function signature takes pre-computed `body_sin_lon`,
    `body_cos_lon`, `body_sin_lat`, `body_cos_lat`, `sin_obliquity`,
    `cos_obliquity` rather than the raw angles. This is consistent
    with the project convention of computing trigs once at init
    (feedback_compute_once.md) — Sun/Moon ephemeris is slow-varying,
    so caching trig values is appropriate. ✓
  - `perturbation_coef` is described as "incorporating μ₃/r₃³ and the
    orbit-averaging factor". This is a deliberate factorization: the
    callee handles the elemental rate algebra; the caller handles the
    physical scaling. The auditor cannot verify the coefficient
    formula from this file alone — must inspect the caller (e.g.
    `compute_solar_perturbation_coefs` / `compute_lunar_perturbation_coefs`
    or whatever upstream prepares it).
  - The line 115 comment `// a₃·sin(i) = z₃ (avoids division by sin(i)
    at equatorial)` is a numerical-stability improvement. The
    expressions that genuinely need sin(i) in the denominator (lines
    158, 167) still have it — the comment refers only to the local
    representation. The singularity at i=0 / i=π is a real property
    of the chosen elements, not an artifact, and the caller must
    decide how to handle equator-aligned orbits.
  - The five rate formulas are written without intermediate equation
    labels (no SR3 §B Eq.(B.k) anchors in comments). For full TBA
    rigor, each of dΩ/dt, dω/dt, de/dt, di/dt, dM/dt should carry
    a comment naming its source equation in Brouwer-Clemence Ch. XI
    or SR3 §B. **Flag for source annotation.**
  - The header docstring (lines 3-24) cites Brouwer-Clemence (1961)
    but does NOT cite Kozai (1959) or Hoots-Roehrich (1980) SR3 §B.
    Given the obvious structural inheritance from SR3 (Lyddane-form
    for SDP4 has exactly this body-direction projection algebra),
    the SR3 §B reference should be added. **Flag for citation.**
  - The "Verdict" of this card is `?` overall, NOT `✓`: the method
    matches the declared closed-form theory cleanly, but
    term-by-term equation correspondence with a primary source
    (SR3 §B or Brouwer-Clemence Ch. XI) is not performed in a
    read-only audit. A full pass requires that symbolic check —
    plausibly a Maxima / SymPy script comparing each rate against
    the SR3 §B Lyddane form.
```

---

## File-level verdict for `third_body.h`

- **A. Error wiring**: ✓ all five rate outputs are `TrackedValue<T>`; no bare `T` in `ThirdBodyRates`. Implicit per-op propagation via `TrackedValue<T>` operator overloads — assumes `sin`/`cos`/`sqrt` resolve to TrackedValue ADL overloads (must be confirmed in `tracked_value.h` audit, not in this file).
- **B. Algebra axioms**: n/a — this is a perturbation formula, not an algebra operation.
- **C. Theoretical basis**:
  - Theory citation: ⚠ Brouwer-Clemence (1961) only — missing Kozai (1959) and Hoots-Roehrich (1980) SR3 §B, which the algebra structurally inherits from.
  - Method ↔ theory match: ✓ closed-form evaluation of orbit-averaged P₂ disturbing function via Lagrange's planetary equations — no Taylor / Padé / continued-fraction mismatch.
  - Bound: ✓ closed-form propagation via TrackedValue ops (REQ-EF-3) is the correct bound family; ? P₂-truncation and orbit-average accuracy belong to the caller, not this function; ? ADL wiring of `sin`/`cos`/`sqrt` un-confirmed at this layer.
  - Term-by-term symbolic correspondence to SR3 §B equations: ? not performed.

**File verdict: PASS-with-flags** — no method-theory mismatch detected (the function does what its docstring says it does, in the way the cited theory prescribes), but three open items warrant a follow-up:

1. Add Kozai (1959) and Hoots-Roehrich (1980) SR3 §B to the header citation block.
2. Annotate each of the five rate-formula lines (lines 152-173) with its source-equation tag (Brouwer-Clemence Ch. XI Eq.(N) or SR3 §B Eq.(B.k)).
3. Perform a term-by-term symbolic cross-check of all five rate expressions vs the SR3 §B Lyddane form (or vs the file-internal derivation `017_third_body_perturbation.md`, which was not opened in this audit).

The function passes the TBA contract structurally; the open items are documentation and verification-rigor flags, not C-fails.

---

## Out-of-scope note

`src/perturbation/third_body.h` is **not** listed in the §3 audit-status table of the framework document (24 files, all under `src/math/`, `src/constants/`, `src/dynamics/`, `src/forces/`, `src/integrators/`). This card is produced as an extension to the framework sweep, on direct request. The framework's §3 table should be updated to include `src/perturbation/` files if perturbation modules are part of the in-scope audit set.
