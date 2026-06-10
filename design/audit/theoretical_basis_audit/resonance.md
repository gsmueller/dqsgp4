# Theoretical Basis Audit — `src/perturbation/resonance.h`

**File**: `src/perturbation/resonance.h` (175 lines)
**Scope**: 3 functions — `detect_resonance`, `initialize_resonance`, `step_resonance`
**Domain**: Deep-space SGP4 (SDP4) tesseral resonance — 24h synchronous & 12h half-day orbits.
**Current state**: Tesseral driving coefficients and acceleration integration are **stubbed** (`TODO`). The 24 deep-space tcppver tests are currently failing — consistent with the unimplemented tesseral terms.

> Verdict legend: ✓ matched / ? unknown-or-unimplemented / ⚠ approximate-or-flagged / ✗ mismatched.

---

## Card 1 — `detect_resonance`

```
=== FORMULA AUDIT CARD ===
ID:                     resonance::detect_resonance
Location:               src/perturbation/resonance.h:63-74
Mathematical statement: ResonanceType(period_min) =
                          SYNCHRONOUS  if 1200 < T < 1800 min
                          HALF_DAY     if  600 < T <  800 min
                          NONE         otherwise

THEORY
  Underlying theorem:   Resonance classification by commensurability of orbital
                        period T with Earth's sidereal rotation period
                        T_⊕ ≈ 1436.068 min. Two windows are admitted:
                          1:1 (synchronous, GEO)  T ≈ T_⊕
                          2:1 (half-day, GPS/Molniya) T ≈ T_⊕ / 2 ≈ 718 min
  Primary reference:    Hoots & Roehrich (1980) SR3 §6 "Deep-space resonance
                        contributions". Hujsak (1979), AIAA 79-136 "A Restricted
                        Four Body Solution for Resonating Satellites Without
                        Drag".
  Domain of validity:   Period in minutes; nonnegative reals. The windows are
                        intentionally wide (~20% / ~15%) to capture
                        eccentricity-driven period drift over the propagation
                        horizon.

METHOD
  Method declared:      Threshold classifier on the scalar period.
  Method implemented:   Two `if (period_min > A && period_min < B)` branches
                        returning the enum; default NONE.
  Match verdict:        ✓ matched — boolean classifier, no numerical method
                        involved. Bounds 1200/1800 and 600/800 are project
                        constants whose origin is the SR3 "deep-space" gate.

ERROR BOUND
  Bound category:       n/a — returns an enum, not a TrackedValue.
  Bound formula:        n/a.
  Bound implemented:    No error propagation. Input `period_min` is a bare
                        `double`, not a TrackedValue.
  Bound verdict:        n/a — discrete classifier. The audit hazard is the
                        threshold itself: a satellite with period 1199 or 801
                        would be classified NONE despite being deep-space.
                        SR3 §6 uses a tighter gate (period > 225 min triggers
                        deep-space; resonance windows then carved within).
                        **Flag**: thresholds 1200/1800 and 600/800 need a
                        primary-source citation; SR3 uses specific node-rate
                        / mean-motion gates rather than period windows.

CROSS-AUDIT
  REQ-EF applies:       n/a (no TrackedValue output)
  AUD-EF applies:       n/a
  AUD-MC applies:       n/a
  Verification test:    tests/test_perturbation — none located that exercises
                        boundary periods (~1200, ~1800, ~600, ~800).

NOTES
  - Magic numbers 1200/1800/600/800 are not derived in-file. SR3 §6 uses the
    225 min gate for deep-space activation and orbits matching specific
    mean-motion conditions for resonance. The current windows are
    conservative-loose approximations.
  - The "20%" and "15%" comments do not match the implemented widths:
    600..800 around 718 is +/-~14%, 1200..1800 around 1436 is roughly
    -16%/+25%. **Flag for tightening**.
```

---

## Card 2 — `initialize_resonance`

```
=== FORMULA AUDIT CARD ===
ID:                     resonance::initialize_resonance
Location:               src/perturbation/resonance.h:93-137
Mathematical statement: Build ResonanceState rs from epoch elements:
                          xlamo (synchronous) = M₀ + Ω₀ + ω₀ − θ_g
                          xlamo (half-day)    = M₀ + 2Ω₀ − 2θ_g
                          xli   = xlamo
                          xni   = n₀
                          atime = 0
                          {d2201,d2211,d3210,d3222,d4410,d4422,
                           d5220,d5232,d5421,d5433} = ? (UNIMPLEMENTED)

THEORY
  Underlying theorem:   SR3 §6 resonance integration variable:
                          λ = M + Ω + ω − θ_g            (synchronous, m=1)
                          λ = M + 2(Ω − θ_g)             (half-day, m=2)
                        where θ_g is GMST and (m=1, m=2) are the tesseral
                        harmonic orders that resonate with the orbit's
                        commensurate frequency.
                        Tesseral driving terms d_{ℓmpq} are products of
                        Kaula inclination functions F_{ℓmp}(i) and Hansen /
                        Kaula eccentricity functions G_{ℓpq}(e) weighted by
                        the J_{ℓm} tesseral coefficient.
  Primary reference:    Hoots & Roehrich (1980) SR3 §6, equations (6.1)-(6.6)
                        for the d-coefficient list and integration variable.
                        Kaula (1966) "Theory of Satellite Geodesy" §3.3
                        (inclination functions F_{ℓmp}) and §3.4 (eccentricity
                        functions G_{ℓpq}).
                        Hujsak (1979) AIAA 79-136 — restricted four-body
                        formulation underlying SR3 resonance treatment.
  Domain of validity:   e < 1; resonance branch fires only when
                        detect_resonance(period_min) ≠ NONE.

METHOD
  Method declared:      (a) Closed-form linear combinations of input angles
                            for xlamo (per SR3 §6).
                        (b) Closed-form Kaula F_{ℓmp}(i) and G_{ℓpq}(e)
                            evaluations for the ten d-coefficients.
                        (c) Initial integrator state copy (xli ← xlamo,
                            xni ← n₀, atime ← 0).
  Method implemented:   (a) ✓ Implemented as written — linear combinations
                            of TrackedValue angles using `exact<T>(2)` for
                            the half-day case.
                        (b) ✗ NOT IMPLEMENTED — d2201..d5433 are
                            default-constructed (zero) and never assigned.
                            Code comment line 134: "TODO: Compute tesseral
                            driving coefficients from inclination functions".
                        (c) ✓ Implemented — xli, xni, atime set as declared.
  Match verdict:        ? partial — angle aggregation matches theory; the
                        Kaula F·G tesseral coefficients are stubbed and the
                        method (whether closed-form Kaula recursion, Allan's
                        tables, or direct series) is undeclared in-file.

ERROR BOUND
  Bound category:       precision (for xlamo, xli, xni, atime) when
                        TrackedValue's closed-form add/multiply propagation
                        applies (REQ-EF-3).
  Bound formula:        (a) xlamo bound: sum-of-input-precisions —
                            err(xlamo)_sync ≤ err(M₀)+err(Ω₀)+err(ω₀)+err(θ_g)
                            err(xlamo)_hd   ≤ err(M₀)+2·err(Ω₀)+2·err(θ_g)
                        (b) d-coefficient bounds: would depend on the chosen
                            F·G evaluation method (closed-form Kaula
                            polynomial vs. Hansen X-series truncation).
                            UNDEFINED while (b) is stubbed.
                        (c) xni = n0 (copy) — bound passes through unchanged.
  Bound implemented:    (a) Operator overloads on TrackedValue add the input
                            precisions per AUD-EF-1..3.
                        (b) Zero — d-coefficients default to TrackedValue with
                            zero value AND zero error. Falsely reports
                            "exact" coefficients while the truth is "unknown".
                        (c) Pass-through copy.
  Bound verdict:        (a) ✓ matched — TrackedValue arithmetic propagation.
                        (b) ✗ UNSOUND — default-constructed TrackedValue
                            implies zero error, but the true value is also
                            zero and the true error is the entire magnitude
                            of the unimplemented tesseral terms. This is a
                            silently-dropped-error C-fail per §2 of the
                            framework. Magnitude: the missing tesseral
                            acceleration is O(J_{22} n²) ≈ 10⁻⁹ rad/min²
                            integrated over deep-space horizons → kilometers
                            of position drift unbounded.
                        (c) ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation through add/mul)
                        — applied only to the angle aggregation.
                        REQ-EF unmet for the d-coefficient slot: stubs lie
                        about precision.
  AUD-EF applies:       AUD-EF-3 (arithmetic ops propagate input errors)
                        passes for the angle line. AUD-EF doesn't currently
                        catch the "stub-as-zero-error" pattern.
  AUD-MC applies:       n/a (no algebra identity exercised here).
  Verification test:    None located. The 24 deep-space tcppver tests
                        FAIL as of HEAD — consistent with the d-coefficient
                        stub.

NOTES
  - The xlamo formula for HALF_DAY is written as `M0 + 2·Ω0 − 2·gmst`,
    omitting the perigee term — consistent with SR3 §6 which combines into
    the angle (M + 2(Ω − θ_g)). Synchronous includes ω. Matches Hujsak
    (1979) formulation.
  - The "Tesseral coefficients would be computed here from the Kaula
    inclination and eccentricity functions" comment is the explicit
    placeholder. The full SR3 implementation requires:
      F_{220}, F_{2,2,1} for synchronous m=2 terms
      F_{321}, F_{322}; F_{4,4,1}, F_{4,4,2}; F_{5,2,2}, F_{5,2,3},
      F_{5,4,2}, F_{5,4,3} for the rest
      G_{210}, G_{211}, etc. for eccentricity dependence
    Citations: Kaula (1966) §3.3 (F-tables), §3.4 (G-tables); SR3 §6.
  - The unused parameters `n0`, `e0`, `i0`, `omega0` for the SYNCHRONOUS
    branch and `e0`, `i0`, `omega0` for HALF_DAY tell the same story —
    they are the F·G arguments that would be consumed once the tesseral
    coefficients are wired.
  - **Status**: stub. Per framework §3, this card is `?` rather than
    `PASS` or `C-fail` because the method for the unimplemented portion
    is undeclared. Card becomes auditable once F·G evaluation is chosen
    (closed-form Kaula vs. Hansen-series truncation, etc.).
```

---

## Card 3 — `step_resonance`

```
=== FORMULA AUDIT CARD ===
ID:                     resonance::step_resonance
Location:               src/perturbation/resonance.h:150-172
Mathematical statement: Integrate the resonance ODE
                          d²λ/dt² = Σ_{ℓm} d_{ℓm} sin(λ − φ_{ℓm}) + (drag)
                        from t = atime to t = tsince, updating
                          xli  (integrated mean longitude λ)
                          xni  (integrated mean motion n)
                          atime (last-step time)

THEORY
  Underlying theorem:   SR3 §6 resonance equation of motion: the tesseral
                        harmonic potential contributes a longitude-dependent
                        torque on the orbit that — for resonant period —
                        does not average to zero over an orbit. The
                        resulting ODE in (λ, n) is integrated separately
                        from the secular elements because its characteristic
                        time is the resonant beat period (months-years),
                        not the orbital period.
                        Vallado (2013) §11.5 outlines the SDP4 resonance
                        leapfrog/RK4 step at 720 min (half-day) and 720 min
                        (synchronous — note same value because the SR3
                        integrator uses a constant 12 h step for both per
                        Hujsak's empirical choice).
  Primary reference:    Hoots & Roehrich (1980) SR3 §6 — defines the step
                        size and the leapfrog-style integration.
                        Vallado et al. (2006) "Revisiting Spacetrack Report
                        #3" Appendix — modern SR3 source listing for the
                        DSPACE subroutine.
  Domain of validity:   Resonance branch only (NONE returns early). The
                        720-min step is empirically chosen to suppress
                        secular aliasing while remaining numerically cheap.

METHOD
  Method declared:      "Simple leapfrog/Euler scheme" at 720-min step,
                        with acceleration from tesseral driving terms.
                        (Comment, line 142-144.)
  Method implemented:   ✗ NOT IMPLEMENTED as declared. Current code:
                          rs.xli   = rs.xlamo + rs.xni * tsince;
                          rs.atime = tsince;
                        This is a **linear extrapolation** from epoch at
                        the initial mean motion — it is neither leapfrog
                        nor Euler of the SR3 ODE, and the tesseral
                        acceleration is entirely absent. `earth_rate` is
                        explicitly discarded `(void)earth_rate;`.
                        The local `T step = T(720)` is computed but never
                        used.
  Match verdict:        ✗ MISMATCHED — declared method (leapfrog/Euler over
                        tesseral acceleration) vs. implemented method
                        (linear extrapolation at constant n) are different
                        in kind. This is the canonical "theory says X,
                        code does Y" failure pattern (framework §0).

ERROR BOUND
  Bound category:       precision (would be, for any iterative integrator).
  Bound formula:        For Euler at step h over horizon H integrating
                        an ODE with Lipschitz constant L and acceleration
                        magnitude A: |error| ≤ (A·h/2)·(e^{LH}−1)/L.
                        For leapfrog (symplectic, 2nd order): error scales
                        as O(h²·H·A). At h=720 min and A ≈ d_{22}·n² with
                        d_{22} ~ 10⁻⁹ rad/min², the step truncation alone
                        is O(h²·A) ≈ 5·10⁵ · 10⁻⁹ = 5·10⁻⁴ rad per step.
                        Over a 1-year horizon (~730 steps), bound is
                        O(0.4 rad) = 23° — significant but bounded.
  Bound implemented:    No bound added. `rs.xli`, `rs.atime` are assigned
                        via TrackedValue arithmetic so multiplication/
                        addition errors propagate, but no truncation,
                        method, or unimplemented-acceleration bound is
                        added to precision.
  Bound verdict:        ✗ UNSOUND. Two failure modes:
                          (i) The linear-extrapolation method has its own
                              truncation error vs. the true SR3 ODE, which
                              is not bounded.
                          (ii) The omitted tesseral acceleration produces
                              a residual drift O(d_{22}·n²·t²/2) =
                              O(10⁻⁹ · t²) [rad] that is silently absent
                              from the precision budget.
                        This is a C-fail per framework Theorem 2.1 (M ⇒ B
                        is false because M itself doesn't follow from T).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (iterative method residual added to
                        precision) — NOT satisfied; no residual computed.
                        REQ-EF-3 (closed-form propagation) partially
                        applied to the linear extrapolation arithmetic.
  AUD-EF applies:       AUD-EF-4 (iterative algorithms add residual) —
                        would FAIL once exercised, but no test currently
                        exists.
  AUD-MC applies:       n/a (no algebra identity exercised).
  Verification test:    The 24 deep-space tcppver test cases (Vallado et
                        al. 2006 verification deck) FAIL on this branch
                        per the task description. Failure is consistent
                        with the stub: linear extrapolation cannot
                        reproduce the resonant longitude oscillation that
                        SR3 expects.

NOTES
  - The local `T step = T(720)` is dead code at present — declared but not
    used.
  - `earth_rate` is consumed only as `(void)earth_rate` — the parameter
    is reserved for the full implementation where the acceleration
    depends on (λ − θ_g(t)) and θ_g(t) = θ_g₀ + ω_⊕·t.
  - The 720-min step for SYNCHRONOUS (period ~1436 min) is suspect: SR3
    uses 720 min for HALF_DAY but ~720 min is also used for synchronous
    in some implementations. The in-file comment "720 min for
    synchronous, 720 min for half-day" is correct per SR3.
  - **Status**: stub + method mismatch. Card is `✗` for Match verdict and
    `✗` for Bound verdict. Both downgrade once tesseral acceleration and
    leapfrog/Euler step are wired per SR3 §6.
```

---

## File-level verdict

- **A. Error wiring**: ⚠ — angle aggregation in `initialize_resonance` propagates correctly via TrackedValue ops; d-coefficients are default-zero (silently zero error for unknown-magnitude quantities); `step_resonance` adds no method/truncation residual.
- **B. Algebra axioms**: n/a — no algebra identity is exercised in this file.
- **C. Theoretical basis**:
  - Card 1 `detect_resonance`: ⚠ classifier OK in form; thresholds need primary-source citation (SR3 §6 gates) — currently magic numbers.
  - Card 2 `initialize_resonance`: ? partial — angle aggregation matches theory; Kaula F·G tesseral coefficients **stubbed** (TODO line 134); no F·G method declared.
  - Card 3 `step_resonance`: ✗ — declared method (leapfrog/Euler over tesseral acceleration) does not match implementation (linear extrapolation at constant n; tesseral acceleration absent; `earth_rate` discarded). Truncation/method bound not added.

**File verdict: C-FAIL (provisional)** — `step_resonance` exhibits the exact "theory says X, code does Y" mismatch the TBA framework was designed to detect, and `initialize_resonance` silently drops error on the unimplemented d-coefficients. This is consistent with the project state (24/24 deep-space tcppver tests failing). Promotion path: wire Kaula F·G evaluation (cite Kaula 1966 §3.3-3.4), implement SR3 §6 leapfrog at 720 min, add per-step Δ as REQ-EF-5 residual, audit thresholds in `detect_resonance` against the SR3 deep-space gate.

## Suggested next steps (out of audit scope)

1. **Cite primary source for thresholds** in `detect_resonance`. SR3 §6 uses specific mean-motion gates rather than period bands; pull those constants in and cite the equation number.
2. **Choose F·G evaluation method** for `initialize_resonance`: closed-form Kaula polynomials (preferred, finite, exact in T) vs. Hansen X-series with truncation (would need REQ-EF-6 truncation bound).
3. **Implement leapfrog/Euler** in `step_resonance` per SR3 DSPACE listing; ensure step residual is added via REQ-EF-5.
4. **Add tcppver regression hook** that exercises the 24 deep-space cases once (3) lands.
