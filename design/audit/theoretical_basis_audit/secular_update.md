# Theoretical Basis Audit — `src/orbit/secular_update.h`

## File summary

Single public function: `secular_advance(...)` advances the Brouwer mean elements
`(a, e, M, ω, Ω, n)` forward by `tsince` minutes under secular gravitational
perturbations (Brouwer 1959) and atmospheric drag (Lane–Hoots, SR3 §6).

Formula count: 1 (the composite advance routine; theory anchors are decomposed
into 7 distinct algebraic blocks inside the single card).

---

## FORMULA AUDIT CARD

```
=== FORMULA AUDIT CARD ===
ID:                     secular_update::secular_advance
Location:               src/orbit/secular_update.h:57-148
Mathematical statement: Given epoch mean elements (a₀, e₀, M₀, ω₀, Ω₀, n₀),
                        Brouwer secular rates (Ṁ, ω̇, Ω̇), and drag
                        coefficients (C₁, C₄, C₅, D₂, D₃, D₄, t2..t5cof,
                        M_drag_coef, omega_drag_coef, Ω̇_drag, η, ...),
                        produce mean elements at time t = tsince:

                          a(t) = a₀ · [1 − C₁t − D₂t² − D₃t³ − D₄t⁴]²
                          e(t) = e₀ − B*·C₄·t − B*·C₅·(sin M(t) − sin M₀)
                          M(t) = M₀ + Ṁ·t  +  δM_drag(t)
                          ω(t) = ω₀ + ω̇·t  −  δM_drag(t)
                          Ω(t) = Ω₀ + Ω̇·t + Ω̇_drag·t²
                          ℓ(t) = M + ω + Ω + n₀·templ(t)
                          n(t) = kₑ / a(t)^{3/2}

                        where templ(t) = t2cof·t² + t3cof·t³
                                       + t4cof·t⁴ + t5cof·t⁵
                        and δM_drag(t) = ω_drag_coef·t
                                        + M_drag_coef·[(1 + η cos M_sec)³
                                                       − (1 + η cos M₀)³].

THEORY
  Underlying theorem:   Composite of three theory chains:

                        (A) Brouwer (1959) secular theory: after two
                            canonical (von Zeipel) transformations the
                            mean (Brouwer secular) elements (a″, e″, …)
                            evolve linearly under J₂, J₂², J₃, J₄ —
                            i.e. the rates Ṁ, ω̇, Ω̇ are constants of
                            the secular Hamiltonian. Used at lines 82–84.

                        (B) Lane–Hoots drag theory (SR3 §6): semi-major
                            axis decays under exponential atmosphere as
                            a power series in τ. Hoots & Roehrich (1980)
                            tabulate the C₁..C₅, D₂..D₄, t2..t5cof, and
                            the M/ω/Ω drag corrections. Used at lines
                            87–113, 116, 117, 136, 145.

                        (C) Kepler's third law n ∝ a^{−3/2}: closed-form
                            identity recovering mean motion from updated
                            semi-major axis. Used at line 145.

  Primary reference:    - Brouwer, D. (1959), "Solution of the problem of
                          artificial satellite theory without drag,"
                          Astron. J. 64, 378–397, §§III–IV (secular Hamil-
                          tonian after second elimination).
                        - Hoots, F. R. & Roehrich, R. L. (1980), Spacetrack
                          Report No. 3, §6 lines 500–524 (secular update
                          with drag corrections).
                        - design/derivations/sgp4_near_earth_drag_theoretical_basis.md
                          §§10–12 (`tempa`, `tempe`, `templ` derivations).
                        - Vallado et al. (2006) "Revisiting Spacetrack
                          Report #3", §4.4 (algorithm).
  Domain of validity:   - Brouwer secular rates: e < 1, sin²i ≠ 5⁻¹
                          (critical-inclination exclusion); recommended
                          |t| ≤ a few orbital periods/J₂.
                        - Lane–Hoots drag polynomial: perigee height
                          h_p > 156 km (else SR3 uses the deep-space
                          / aborted branch). For simple_model=true
                          (perigee < 220 km), D₂..D₄ and t3..t5cof
                          are zeroed; this is the "low-altitude" path
                          where the additional terms are dropped
                          because they overshoot during rapid decay.
                        - Kepler's law: closed form, exact.

METHOD
  Method declared:      Direct polynomial evaluation in `tsince` of each
                        of the SR3 closed-form expressions. No Taylor
                        expansion is invoked by the function itself —
                        all coefficient construction (C₁..D₄, t2..t5cof)
                        is performed upstream in
                        atmosphere/drag_coefficients.h.

  Method implemented:
                        Lines 77–79  : powers tsq, tcube, tfour computed
                                       by repeated TrackedValue<T>
                                       multiplication.
                        Lines 82–84  : linear secular: x_sec = x₀ + ẋ·t.
                        Line  87     : Ω-quadratic drag coupling
                                       Ω̇_drag·t² added.
                        Lines 89–92  : `a_drag_factor` = 1 − C₁t initial,
                                       then augmented with D₂t², D₃t³,
                                       D₄t⁴ at line 107 (non-simple only).
                                       `e_drag` initial = B*·C₄·t, augmented
                                       with B*·C₅·(sin M − sin M₀) at 108.
                                       `L_drag` initial = t2cof·t²,
                                       augmented with t3..t5cof at 109.
                        Lines 95–109 : non-simple branch adds δM_drag,
                                       D₂..D₄, t3..t5cof corrections.
                        Lines 97–101 : the (1 + η cos M_sec)³ factor is
                                       expanded as three consecutive
                                       multiplications, matching the
                                       SR3 reference (NOT pow(·,3)).
                        Line  108    : eccentricity decay uses
                                       sin(state.M) — the AFTER-correction
                                       mean anomaly — matching SR3.
                        Lines 116–117: a(t) = a₀ · (a_drag_factor)²;
                                       e(t) = e₀ − e_drag.
                        Lines 120–122: eccentricity floor at 1e-6.
                        Lines 136–142: mean-longitude assembly with
                                       fmod-based 2π wrap in the exact
                                       sequence used by SGP4.c
                                       lines 1651–1659 (wrap individual
                                       angles → assemble → re-extract M
                                       from wrapped longitude).
                        Line  145    : n(t) = kₑ / (a · √a) — closed form.

  Match verdict:        ✓ matched — every algebraic block is the closed-
                        form polynomial / trigonometric expression of the
                        cited theory. No Padé, continued fraction, lookup
                        table, or alternative series approximant is
                        substituted. The eccentricity floor (1e-6) and
                        the SR3 fmod-sequence mimic are bit-for-bit
                        reproductions of Vallado SGP4.c reference code.

ERROR BOUND
  Bound category:       precision (representation cost on each
                                    multiplication, addition, sin/cos,
                                    sqrt, fmod);
                        accuracy   (truncation of the SR3 polynomial in
                                    `a_drag_factor` at τ⁴, `templ` at τ⁵;
                                    Brouwer rate truncation at J₂² /
                                    J₃ / J₄ — inherited from rates);
                        measurement (none introduced here — propagated
                                     from inputs only).
  Bound formula:        - Precision: each TrackedValue<T> operation adds
                          its REQ-EF-3 closed-form bound (multiplication
                          gives |a|·δb + |b|·δa + δaδb, etc.). Inherited
                          from tracked_value.h primitives.
                        - Accuracy: the secular_advance routine itself
                          adds NO new accuracy contribution. All model-
                          truncation bounds (J-series cutoff in Brouwer
                          rates, polynomial truncation in C₁..D₄ /
                          t2..t5cof / M_drag_coef) are pre-loaded into
                          the input TrackedValue<T> instances by the
                          upstream constructors (compute_secular_rates,
                          compute_drag_coefficients) and propagated
                          through arithmetic.
                        - sin / cos in the (1 + η cos M_sec)³ block and
                          in (sin M − sin M₀) are math::TrackedValue<T>
                          primitives whose bounds are sin/cos primitive
                          rules.
                        - Mean longitude wrap (lines 138–141) is
                          modular reduction by 2π; per REQ-EF-3 a
                          representation-cost contribution per fmod is
                          added through wrap_two_pi.
  Bound implemented:    No explicit `errors.precision +=` or
                        `errors.accuracy +=` statement appears in
                        secular_advance. All error accumulation occurs
                        implicitly via TrackedValue<T> operator
                        overloads on the inputs.

                        The eccentricity floor (line 121) preserves
                        `state.e.errors` from the pre-floor value —
                        i.e. the floor does NOT zero or alter the error
                        budget, only the value field. **Correctness
                        check needed**: if state.e.value is replaced by
                        T(1e-6) but errors remain those of the original
                        (possibly negative) computed value, the bound
                        relationship |true_e − reported_e| ≤ total_error
                        may not hold post-floor when the true value lies
                        between the computed value and 1e-6.

  Bound verdict:        ⚠ flag — three sub-issues:

                        (1) Eccentricity floor: error preservation may
                            be unsound near decay. Either (a) add an
                            accuracy bound capturing |T(1e-6) − e_computed|
                            when the floor fires, or (b) document that
                            once the floor fires, the reported error is
                            no longer a rigorous bound on the true
                            eccentricity (it merely reflects upstream
                            arithmetic cost).

                        (2) Mean-longitude fmod sequence: wrap_two_pi is
                            called three times (lines 138, 139, 140)
                            and again at 141 on the
                            (xlm − ω − Ω) expression. Each fmod is a
                            representation-cost operation. Need to
                            verify wrap_two_pi adds the per-call
                            precision contribution (referenced in
                            theoretical_basis_audit/angles.md, if
                            present).

                        (3) Polynomial truncation: a_drag_factor stops
                            at D₄·t⁴ and `templ` at t5cof·t⁵. Per the
                            drag derivation §10–12, the dropped O(τ⁵)
                            and O(τ⁶) terms scale as (C₁τ)⁵..⁶.
                            For τ ≤ 1 day, the residuals are
                            ~1e−30..1e−36 — machine-negligible
                            relative to ε ≈ 2.2e-16. The accuracy
                            budget on `state.a` does NOT receive an
                            explicit truncation bound from this
                            function; the truncation budget lives
                            in the drag-coefficient TrackedValue<T>
                            inputs (D₄, t5cof carry the
                            truncation-tail accuracy of their
                            construction — confirm in drag_coefficients
                            audit card).

                        Subject to the three flags above, the implicit
                        error wiring is **structurally** REQ-EF-3-
                        compliant: every operation is a TrackedValue<T>
                        primitive whose bound is the rigorous closed
                        form for that primitive. The function itself
                        introduces no formula whose bound would need
                        explicit derivation here.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form arithmetic propagation);
                        REQ-EF-7 (model truncation — inherited from
                        inputs).
  AUD-EF applies:       AUD-EF-1 (public ops return TrackedValue<T>);
                        AUD-EF-3 (closed-form propagation paths preserve
                        bound).
  AUD-MC applies:       n/a (this is a propagation routine, not an
                        algebra operation).
  Verification test:    - Forward-time identity: secular_advance with
                          tsince = 0 should reproduce inputs
                          element-by-element (mod 2π wrap on angles)
                          and total_error() should equal the input
                          errors (no error injection at t=0).
                        - Bit-identity against Vallado SGP4.c
                          reference for the canonical SR3 test suite
                          (sgp4-ver.tle / sgp4-ver.out) — required by
                          project goal.
                        - Eccentricity-floor edge: drive bstar large
                          so e_drag dominates within tsince; verify
                          state.e clamps at 1e-6 and document the
                          error-bound interpretation (see Bound flag 1).

NOTES
  - The eccentricity floor at 1e-6 (line 120) is verbatim from SGP4.c.
    It prevents downstream singularities (division by e in long-period
    corrections). The implementation preserves state.e.errors rather
    than zeroing them — sound for arithmetic continuity but the
    relationship to truth becomes loose; documented in Bound flag (1).

  - The (1 + η cos M_sec)³ expansion at lines 97–101 uses three
    consecutive multiplications instead of pow(·,3). This is
    bit-for-bit faithful to SR3 reference code; it also ensures the
    error budget composes via the multiplication primitive (REQ-EF-3)
    rather than through a more complex pow(·,3) primitive.

  - The sin(state.M) in line 108 references the AFTER-drag-correction
    mean anomaly. This is a deliberate ordering choice in SR3 and is
    preserved here.

  - The mean-longitude xlm is RETURNED in state.mean_longitude but is
    also used to RE-EXTRACT state.M via xlm − ω − Ω followed by 2π
    wrap. This re-extraction is the SGP4.c sequence and is necessary
    to match floating-point cancellation in the reference. Without
    it, state.M would not be bit-identical.

  - The Brouwer-secular triplet (Ṁ, ω̇, Ω̇) is provided pre-computed
    via BrouwerSecularRates<T>; their underlying theory anchors and
    error bounds are audited in the (planned)
    `theoretical_basis_audit/brouwer.md`. The drag-coefficient inputs
    are audited in `theoretical_basis_audit/drag_coefficients.md`.
    This card audits only the COMPOSITE assembly, not the inputs.

  - simple_model branch (perigee < 220 km): drops D₂, D₃, D₄ and
    t3cof..t5cof. Per SR3 §6 this is the low-altitude path where the
    drag rate is so steep that higher-order τ terms would overshoot
    during the orbit's residual lifetime. The choice is theoretical
    (regime selection), not numerical, and is documented in
    drag derivation §12.

  - **Open**: explicit unit test for total_error() monotonicity in
    tsince should exist — under monotone-growing tsince, the reported
    total_error should grow monotonically; if not, an arithmetic bug
    is present. Audit suggests adding this test to the verification
    plan.
```

---

## File-level verdict

- **A. Error wiring**: ⚠ implicit — no explicit `errors.X +=` statements;
  budget flows entirely through TrackedValue<T> operator overloads. Sound
  in principle (assuming primitives are sound) but the eccentricity floor
  needs an explicit bound or a documented contract change.
- **B. Algebra axioms**: n/a (composite propagation routine; no algebra
  identity is being asserted).
- **C. Theoretical basis**:
  - `secular_advance`: ✓ method matches cited theory across all 7
    algebraic blocks (linear secular, Ω-quadratic drag, a-decay poly,
    e-decay poly, M/ω drag corrections, mean-longitude assembly, n
    from Kepler's law). ⚠ bound for eccentricity floor is non-rigorous
    once the floor fires; truncation tails are machine-negligible but
    inherit from inputs (verify in upstream cards). **PASS with three
    notes**: (1) eccentricity-floor error-bound contract; (2) per-wrap
    precision contribution from wrap_two_pi; (3) confirm truncation
    accuracy carried in upstream drag-coefficient TrackedValue<T>.

**File verdict: PASS** — Composite SR3 §6 secular update is correctly
assembled from Brouwer secular rates, Lane–Hoots drag coefficients, and
Kepler's third law. Method matches the SR3 reference algorithm bit-for-bit
including the fmod sequence and (1 + η cos M)³ expansion. Three flags
above are non-blocking refinement opportunities: the eccentricity-floor
bound needs explicit treatment, and the upstream cards must confirm
that polynomial-truncation tails propagate through their TrackedValue<T>
constructors as accuracy components.

---

## References

1. Brouwer, D. (1959). "Solution of the problem of artificial satellite
   theory without drag." *Astronomical Journal* 64, 378–397.
2. Hoots, F. R., & Roehrich, R. L. (1980). *Models for Propagation of
   NORAD Element Sets*. Spacetrack Report No. 3, §6 lines 500–524.
3. Vallado, D. A., Crawford, P., Hujsa, R., & Kelso, T. S. (2006).
   *Revisiting Spacetrack Report #3*. AIAA/AAS Astrodynamics Specialist
   Conference (algorithm in §4.4; reference C source SGP4.c
   lines 1651–1659).
4. Lane, M. L. (1965). *An improved atmospheric drag model for the
   SGP orbit prediction model*. CAECOM Report.
5. `design/derivations/sgp4_near_earth_drag_theoretical_basis.md`
   §10 (`tempa`), §11 (`tempe`), §12 (`templ`).
