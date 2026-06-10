# Theoretical Basis Audit — `src/sgp4/deep_space.h`

**File**: `src/sgp4/deep_space.h` (329 lines)
**Functions audited**: 2
**Status**: 24/24 SDP4 tcppver tests currently FAILING → expect many ⚠ / ? entries.

This file is a **composition orchestrator** for SDP4 (Hoots-Roehrich 1980 §6).
It contains *very little* original math; almost every formula delegates to a
standalone submodule (`perturbation/third_body.h`, `perturbation/resonance.h`,
`perturbation/short_period.h`, `orbit/osculating_elements.h`,
`orbit/state_from_elements.h`). The cards below evaluate the **orchestration**:
unit conversions, third-body coefficient construction, combined-rate assembly,
and the propagation control flow.

## Card 1 — `initialize_deep_space`

```
=== FORMULA AUDIT CARD ===
ID:                     deep_space::initialize_deep_space
Location:               src/sgp4/deep_space.h:107-245
Mathematical statement: Build a DeepSpaceInit<T> struct containing all
                        precomputed constants required by SDP4 propagation:
                        recovered Brouwer elements (a₀,n₀,β₀), Brouwer
                        secular rates (Ṁ,ω̇,Ω̇), solar and lunar third-body
                        rate corrections (Ṁ,ω̇,Ω̇,ė,i̇)_⊙,☾, resonance
                        state, long-period (xlcof, aycof) coefficients,
                        obliquity sin/cos, and trig polynomials.

THEORY
  Underlying theorem:   SDP4 deep-space model (Hoots & Roehrich 1980,
                        Spacetrack Report No. 3, §6). Composed of:
                          (a) Brouwer (1959) secular rates from J₂, J₂², J₄.
                          (b) Lunar/solar third-body theory of Lane &
                              Cranford (1969) / Hujsak (1979) — orbit-
                              averaged disturbing function for an external
                              body, giving secular rates dM/dt, dω/dt,
                              dΩ/dt, de/dt, di/dt.
                          (c) Resonance detection for geosynchronous
                              (1-day) and Molniya/12-h (1/2-day) cases,
                              with tesseral harmonic m=2,3,4 (geo) and
                              m=1 (12-h) couplings — Hujsak (1979).
                          (d) Long-period periodic coefficients
                              xlcof, aycof from Brouwer J₃/J₂.
  Primary reference:    Hoots & Roehrich (1980) SR3 §6, pages 13-16
                        (DSINIT, DSCOM, DPPER, DSPER).
                        Hujsak (1979), "A Restricted Four Body Solution
                        for Resonating Satellites Without Drag".
                        Lane & Cranford (1969) AIAA paper 69-925.
  Domain of validity:   Orbital period ≥ 225 min. Eccentricity e₀ ∈ [0,1).
                        Inclination i₀ ∈ [0,π].

METHOD
  Method declared:      Composition of standalone submodel calls:
                          1. orbit::recover_mean_elements (Brouwer recovery)
                          2. config.model_functions.secular_rates (Brouwer)
                          3. ephemeris::compute_solar_position
                          4. ephemeris::compute_lunar_position
                          5. perturbation::compute_third_body_rates (×2)
                          6. perturbation::initialize_resonance
                        Plus closed-form unit conversions and trig.
  Method implemented:   Matches declared in shape. Specifics inspected:
                          • Trig constants (lines 137-142): closed-form
                            (cos, sin, ×, polynomial in cosio).
                          • Element recovery (144-151): delegated.
                          • Brouwer rates (153-159): delegated.
                          • Obliquity (165-167): closed-form cos/sin.
                          • Third-body coefficient (line 196, 204):
                              C = (3/2) · μ₃ · sec_per_min² / (r_e³ ·
                                  r₃_er³ · n₀ · a₀²)
                            built from numeric literals GM_sun, GM_moon,
                            r_sun_AU, r_moon_km via measured() constants.
                          • Third-body rate calls (207-219): delegated.
                          • Combined rates (222-226): trivial addition.
                          • Resonance detect (229-232): delegated.
                          • Long-period xlcof/aycof (234-242): closed-form
                            algebraic from A₃₀/CK₂ and sinio/cosio.
  Match verdict:        ⚠ structurally matches declared composition, but
                        FOUR concrete concerns visible at this level
                        (see NOTES). Submodule cards must vouch for the
                        individual deliveries.

ERROR BOUND
  Bound category:       precision (numerical) + accuracy (input TLE) +
                        measurement (physical constants GM_sun, GM_moon,
                        r_sun, r_moon, obliquity, etc.).
  Bound formula:        Composite. Each delegated subcall returns
                        TrackedValue<T> with its own bound. Closed-form
                        composition propagates per REQ-EF-3. The numeric
                        literals embedded here at lines 192, 194, 200, 202
                        carry **measurement** error via
                        TrackedValue::measured(value, uncertainty):
                          • GM_sun  = 1.32712440018e11 km³/s² ± 1e4
                          • GM_moon = 4902.8 km³/s² ± 0.3
                          • r_sun   = 1.495978707e8 km ± 1.0
                          • r_moon  = 384400 km ± 1
                        The constant r_sun is a mean distance, but the
                        Sun's actual distance varies ±1.7% over the year;
                        r_moon varies ±6.7% over the month. Using a mean
                        with ±1 km uncertainty under-bounds the **model**
                        error: third-body perturbation scales as 1/r³,
                        so r±1.7% → coefficient ±5%.
  Bound implemented:    Numeric literals carry the listed measurement
                        bounds. Compose via REQ-EF-3 in TrackedValue ops.
  Bound verdict:        ⚠ unsound for r_sun and r_moon: the recorded
                        ±1 km measurement uncertainty is the **fit
                        precision** of the mean distance, NOT the
                        physical variation envelope across the orbit.
                        For deep-space accuracy this is a meaningful
                        under-bound (factor ~10⁴ low for sun, ~10⁵ low
                        for moon). Should be promoted to ±2.5e6 km
                        (sun) and ±2.5e4 km (moon), or — better — the
                        third-body perturbation coefficient should be
                        rebuilt at each step from the current sun/moon
                        ephemeris distance rather than this single mean.
                        See **Concern 3** in NOTES.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation),
                        REQ-EF-1 (accuracy carried with values),
                        REQ-EF-2 (measurement carried with constants).
  AUD-EF applies:       AUD-EF-1, AUD-EF-2, AUD-EF-3.
  AUD-MC applies:       n/a (no algebra operation here; all algebra
                        lives in submodules).
  Verification test:    tests/test_sgp4/test_sdp4_tcppver.cpp (or
                        equivalent SDP4 reference comparator). **Currently
                        failing all 24 deep-space tcppver tests** —
                        flagged on the file's status line.

NOTES — concrete concerns visible at this orchestration level

  Concern 1 — Third-body coefficient: form vs SR3.
    Lines 192-205 build perturbation_coef_{sun,moon} as
        C = (3/2)·μ₃·(60²/r_e³)/(r₃_er³ · n₀ · a₀²)            (1)
    The doc comment lines 175-180 cites this as "(3/2) × (μ₃/r₃³) /
    (n × a²)" in SGP4 normalized units. SR3 §6 (variables CC, ZNS,
    ZNL etc.) tabulates the third-body coefficient as something like
        ZNS = 1.19459e-5  (rad/min, mean motion of sun)
        ZNL = 1.5835218e-4 (rad/min, mean motion of moon)
        CS  = 2.9864797e-6 / a^(3/2) etc.
    The constants in SR3 are *closed-form composites* of physical GM,
    r, n₀, a₀. The form here (1) is dimensionally a per-time
    coefficient; SR3 absorbs the n_third_body factor into the
    coefficient. **Verification needed**: does compute_third_body_rates
    consume coefficient in the same form (3/2·μ₃/(r³·n·a²)) as SR3's
    ZNS/ZNL convention, or is there a missing factor of n_sun, n_moon,
    or a different power of a₀? This is exactly the kind of
    convention mismatch that would yield wholesale tcppver failures.
    Verdict: ?

  Concern 2 — Mean-distance r₃ vs instantaneous r₃.
    Lines 192-205 use a FIXED mean distance for r_sun = 1 AU and
    r_moon = 384400 km. SR3 §6 computes third-body rates using
    mean distances (this matches Lane-Cranford / Hujsak who orbit-
    average the disturbing function). Hujsak 1979 explicitly uses
    mean distance, mean motion. So the **method choice** here (use
    of mean r₃) is consistent with the cited theory. But the
    measurement uncertainty recorded (±1 km) does NOT cover the
    intrinsic model error from the mean-distance approximation,
    which is O(eccentricity_of_third_body × coefficient).
    Verdict: method ✓, bound ⚠.

  Concern 3 — Sun ephemeris at t=0 only.
    Lines 170-172 compute sun_pos and moon_pos at t=t_zero=0 once,
    during initialization, and these are passed to
    compute_third_body_rates which presumably uses them as the mean
    angles for the third-body orbit. For deep-space propagation
    over months/years, the third-body position changes substantially.
    SR3 §6 handles this by treating the sun/moon mean elements as
    time-varying (computed each call to DSPER from the secular
    averaged-over-one-orbit-of-the-perturber theory). **If
    propagate_deep_space does not re-call compute_third_body_rates
    with updated sun/moon positions, the secular-rate corrections
    are frozen at epoch — which is correct for the **secular**
    contribution but **omits the long-period (~SR3 DPPER) lunar/
    solar terms that depend on time-varying sun/moon angles.**
    Verdict: ? — depends on what propagate_deep_space does with
    the time-varying part. See Concern 5.

  Concern 4 — Resonance period gating via double().
    Line 232 passes `double(ds.period_min.value)` — explicit
    downcast to `double` even when T is a wider type (e.g.
    cpp_bin_float_50). The resonance dispatcher receives only
    double-precision period; for wide-T users this drops precision
    irrelevantly (period gating uses ~minutes-of-period, not
    ULP-level), so it is not a correctness failure but is a
    **code consistency** wart (AUD-CC scope).
    Verdict: ✓ (gating only).

  Concern 5 — DPPER (deep-space long-period) terms absent.
    The function does NOT compute or store the SR3 DPPER long-
    period correction state (the periodic lunar/solar corrections
    to e, i, Ω, ω, M that depend on the third-body's current
    longitude). SR3 DPPER is invoked at every propagation call,
    not just init. So `initialize_deep_space` correctly does NOT
    compute it. **But there is no call to a `deep_periodics`
    submodule in `propagate_deep_space` either** (see Card 2,
    Concern A). DPPER absence is the likely culprit for the 24/24
    failures.
    Verdict: ⚠ DPPER missing pipeline-wide.

  Concern 6 — Long-period coefficient branch at cosio→±1.
    Lines 236-241 branch xlcof: the main form
        xlcof = (1/8)·(A₃₀/CK₂)·sinio·(3+5·cosio)/(1+cosio)
    becomes 0/0 as cosio→−1. The fallback at line 240
        xlcof = (3/8)·(A₃₀/CK₂)·sinio
    is the L'Hôpital / Taylor limit at cosio=−1 (where sinio→0 too,
    so the whole term should vanish — it does: sinio·finite = 0).
    Threshold 1.5e-12 is appropriate for double; for wide T this
    threshold is too coarse but the L'Hôpital limit is exact at
    cosio=±1, so the branch is rigorous.
    Verdict: ✓ for double, ⚠ threshold for wide T (cosmetic).

  Concern 7 — `axN` etc. computed during init are NOT here.
    The init function does NOT compute axN, ayn — those are at
    line 290-296 of propagate_deep_space. So this init is purely
    constants. Good separation.
    Verdict: ✓.
```

## Card 2 — `propagate_deep_space`

```
=== FORMULA AUDIT CARD ===
ID:                     deep_space::propagate_deep_space
Location:               src/sgp4/deep_space.h:260-326
Mathematical statement: Given DeepSpaceInit and time tsince [min],
                        produce a TEME StateVector (r, v) at tsince.
                        Pipeline:
                          1. Secular advance with combined rates:
                             M = M₀ + Ṁ_total · t,  ω = ω₀ + ω̇_total · t,
                             Ω = Ω₀ + Ω̇_total · t,  e = e₀ + ė_total · t,
                             i = i₀ + i̇_total · t.
                          2. e clamp: e ← max(e, 1e-6).
                          3. a = a₀  (no drag in deep-space simplified model).
                             n = xke / a^(3/2).
                          4. Mean longitude xl = M + ω + Ω.
                          5. Long-period periodics (Brouwer J₃):
                               axN = e·cos(ω)
                               β² = 1 − e²
                               temp_lp = 1/(a·β²)
                               xll = temp_lp · xlcof · axN
                               aynl = temp_lp · aycof
                               xlt = xl + xll
                               ayn = e·sin(ω) + aynl
                          6. Modified-Kepler solver:
                             capu = wrap_2π(xlt − Ω); solve for E+ω.
                          7. Osculating elements (delegated).
                          8. Short-period J₂ corrections (delegated).
                          9. State-from-elements (delegated).

THEORY
  Underlying theorem:   SDP4 propagation as Brouwer-mean-element advance
                        with J₂ short-period and J₃ long-period
                        corrections (Hoots-Roehrich 1980 SR3 §3-§7),
                        plus deep-space DPPER lunar/solar periodic
                        corrections (SR3 §6) and resonance integration
                        for tesseral commensurabilities (SR3 §6.B,
                        Hujsak 1979).
  Primary reference:    Hoots & Roehrich (1980) §3-§7, especially
                        §6 DPPER (lunar/solar periodics) and §6 DSPER
                        (resonance integration).
                        Brouwer (1959) for the J₂ secular and periodic
                        forms.
  Domain of validity:   tsince ∈ ℝ (signed). e clamped to ≥ 1e-6 to
                        avoid singular Kepler. inclination ∈ [0,π].

METHOD
  Method declared:      Time-step the satellite by:
                          • Linear secular advance over the **combined**
                            Brouwer+solar+lunar rates.
                          • Apply long-period and short-period periodics.
                          • Solve modified Kepler for E+ω.
                          • Convert to inertial state.
                        Deep-space-specific additions over near-space:
                          • Combined rates include solar/lunar
                            contributions.
                          • Resonance integration should be active when
                            ds.resonance indicates commensurability.
                          • Deep-space long-period (DPPER) corrections
                            should be applied between secular advance
                            and Kepler.
  Method implemented:   Lines 270-325 implement:
                          • Secular advance (✓, lines 270-275).
                          • e clamp (✓, lines 278-280).
                          • a = a₀, n recomputed (lines 283-284).
                          • Long-period (J₃) periodics (✓, 290-296).
                          • Kepler via injected lambda (✓, 299-300).
                          • Osculating (✓, 303).
                          • Short-period (✓, 313-320).
                          • State-from-elements (✓, 323-325).
                        **MISSING vs declared method**:
                          (A) No call to a `deep_periodics` / DPPER
                              submodule. The lunar/solar **periodic**
                              corrections (SR3 §6 DPPER) are absent
                              entirely.
                          (B) No call to a `integrate_resonance` /
                              DSPER submodule. The `ds.resonance` state
                              is initialized but never read or
                              integrated forward. For resonant
                              satellites (geo, 12-h Molniya, GPS), this
                              means the resonance-driven Δn, ΔM, ΔΩ
                              corrections that drift the satellite back
                              to its sub-orbit are entirely absent.
                          (C) The third-body rates are frozen at
                              epoch sun/moon positions (per
                              `initialize_deep_space` Card 1, Concern 3).
                              SR3 §6 has the secular rates also depend
                              on the **current** (slowly-drifting)
                              sun/moon mean longitudes — implemented
                              by re-evaluating DSCOM-like coefficients
                              each call. Here the rates are constants.
  Match verdict:        ✗ MISMATCHED — the deep-space-specific pieces
                        (DPPER long-period, DSPER resonance, time-
                        varying third-body angles) are absent. The
                        function reduces to "near-space with constant
                        solar+lunar secular rate offsets", which is
                        not SDP4. This is the structural explanation
                        for 24/24 SDP4 tcppver failures.

ERROR BOUND
  Bound category:       precision (numerical) + accuracy (model fidelity).
  Bound formula:        For a SDP4-faithful implementation, the
                        precision bound is the composition of:
                          • Each TrackedValue closed-form op (REQ-EF-3).
                          • Kepler iterative residual (REQ-EF-5).
                          • Truncation in series-evaluation sub-ops.
                        **Accuracy** bound should additionally include
                        the **omitted-term magnitudes**:
                          • DPPER long-period (O(GM_moon/GM_earth ·
                            (a/r_moon)³ · ω_moon · tsince) per
                            element, modulated by sin/cos of third-body
                            longitude) — magnitudes ~1e-5 to 1e-3 rad
                            per day, integrating to large errors over
                            propagation timescales.
                          • Resonance Δn (resonant geosynchronous
                            satellites can drift ~degrees/day if
                            uncorrected — many orders of magnitude
                            above tcppver tolerance).
                          • Time-varying sun/moon angles in the secular
                            rates (~1° per ~360 days for solar mean
                            longitude — bounded by tsince·n_sun in
                            angle, with rate-on-rate effect ~e×ZNS).
  Bound implemented:    Only the per-step closed-form propagation
                        bounds. **No explicit "model fidelity" budget
                        for the missing DPPER/DSPER/time-varying
                        contributions** is added to `accuracy`.
  Bound verdict:        ✗ UNSOUND. The reported total_error() from this
                        function omits the dominant deep-space
                        contributions. Any go/no-go decision based on
                        it would be optimistic by orders of magnitude
                        for resonant or long-time-base propagations.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3, REQ-EF-5 (Kepler iter), and a
                        **needed-but-absent** REQ on omitted-term
                        accuracy accounting.
  AUD-EF applies:       AUD-EF-3 (closed-form ops), AUD-EF-4 (Kepler
                        residual).
  AUD-MC applies:       n/a (composition only).
  Verification test:    tests/test_sgp4/test_sdp4_tcppver.cpp
                        (24/24 currently failing → matches diagnosis).

NOTES — concrete concerns visible at this orchestration level

  Concern A — DPPER long-period periodics MISSING (critical).
    SR3 §6 DPPER computes lunar/solar long-period corrections to
    e, i, ω, M, Ω that depend on sin/cos of the **current** sun
    and moon mean longitudes. These are applied after the secular
    advance and before Kepler. They are entirely absent from this
    function. For GEO/Molniya/GPS the DPPER corrections are
    ~10⁻³ to 10⁻¹ radians and dominate the deep-space error
    budget over the day-to-month tcppver propagation intervals.
    THIS IS LIKELY THE PRIMARY CAUSE of 24/24 failures.
    Verdict: ✗.

  Concern B — DSPER resonance integration MISSING (critical for
    resonant satellites).
    SR3 §6 DSPER integrates the resonance differential equations
    (Hujsak 1979): for 1-day resonance, ṅ = …terms in
    sin(2(λ−θ_g)) etc.; for 12-h, ṅ = …terms in sin(λ−θ_g) etc.
    Here `ds.resonance` is initialized but never consumed.
    Verdict: ✗ for resonant satellites. (For non-resonant deep-
    space — i.e. period ≥ 225 min but not at exact 1-day or
    12-h commensurability — DSPER is bypassed in SR3 as well,
    so absence is only failure for the resonant subset.)

  Concern C — Time-varying sun/moon angles in secular rates.
    ds.solar_rates and ds.lunar_rates are stored constants from
    initialize_deep_space. They depend on sun_pos.cos_longitude,
    sun_pos.sin_longitude (etc.) at t=0. The sun's longitude
    drifts ~1° per day; the moon's ~13°/day. Over a 30-day
    propagation the secular-rate contributions of solar/lunar
    perturbations rotate substantially. SR3 §6 DSCOM-equivalent
    re-evaluates these each call.
    Verdict: ⚠ secondary cause of tcppver failures, smaller
    effect than DPPER but still meaningful.

  Concern D — a = a₀ "no drag for deep-space in simplified model".
    Line 283 comment acknowledges the simplification. SR3 §6
    actually does include a drag-like contribution for low-
    perigee deep-space satellites (the η, ξ, etc. terms in
    DSINIT/DSPER §6.A). For high-altitude purely-geosynchronous
    satellites drag is negligible (~10⁻¹⁵ ER/min²), so the
    simplification is defensible. **However**: the documentation
    of the simplification should explicitly state the bound on
    drag-induced error, which it doesn't.
    Verdict: ⚠ correct simplification, undocumented bound.

  Concern E — Eccentricity clamp e ← max(e, 1e-6).
    Lines 278-280: protects Kepler from e=0 singularity, but does
    NOT carry a bound for the perturbation of e introduced by the
    clamp. If the secular advance drives e below 1e-6 (possible
    for near-circular deep-space orbits after long t), the clamp
    silently injects up to 1e-6 in the e value. This should be
    added to `e.errors.precision` or `e.errors.accuracy` as a
    bound on the clamp-induced perturbation.
    Verdict: ⚠ silent dropout in error budget.

  Concern F — `inc` variable computed but mostly unused.
    Line 275 computes inc = i₀ + i̇_total · tsince. It is used in
    short-period corrections (cos_i, sin_i, cos2_i — lines 306-
    311). Good. But the long-period coefficients xlcof and aycof
    were precomputed at i = i₀ (init time), and **not updated**
    for the current inclination. For long deep-space propagations
    where i drifts noticeably (i̇_total tsince can be O(degrees)),
    using xlcof(i₀) as if it were xlcof(i(t)) introduces a small
    error. SR3 standard practice is to use mean (epoch) elements
    for long-period coefficients, so this is **consistent with
    SR3 convention** — not a defect of the implementation, but a
    convention worth noting.
    Verdict: ✓ matches SR3 convention.

  Concern G — `wrap_two_pi` only at one site.
    Line 299 wraps capu but not xl, xlt, M, etc. SDP4 lambda
    accumulators can drift by many revolutions over multi-year
    propagation; wrapping only at the Kepler input is sufficient
    because cos/sin in long-period periodics handle unwrapped
    arguments correctly, and Kepler is the only iterative
    consumer. ✓.
```

## File-level verdict for `deep_space.h`

| Dimension | Verdict | Notes |
|---|---|---|
| A. Error wiring | ⚠ | Numerical literal bounds for r_sun, r_moon under-state the **physical variation** (Card 1 Bound verdict). Eccentricity clamp drops a silent precision contribution (Card 2 Concern E). |
| B. Algebra axioms | n/a | This file is composition; algebra lives in submodules. |
| C. Theoretical basis | **✗ FAIL** | `propagate_deep_space` is missing DPPER lunar/solar long-period (Concern A), DSPER resonance integration (Concern B), and time-varying third-body angles (Concern C). The cited theory (SR3 §6) is not implemented. The implementation is effectively "near-space with constant secular offsets" — which explains the 24/24 SDP4 tcppver failures. |

**Primary defect**: The deep-space-specific machinery (DPPER, DSPER, DSCOM
time-evolution) declared by the file header docstring and theory citation is
absent from the propagation step. The structure of `propagate_deep_space`
needs three additional delegated submodule calls (analogous to the existing
`perturbation::compute_third_body_rates` and `perturbation::apply_short_period`):

  1. `perturbation::apply_deep_periodics(...)` — DPPER long-period.
  2. `perturbation::integrate_resonance(ds.resonance, tsince, ...)` —
     DSPER resonance update of n, M.
  3. `perturbation::recompute_third_body_rates_at(tsince, ...)` — or fold
     time-varying sun/moon angles into the third-body submodule.

**Secondary defects**: bound under-statement for r_sun/r_moon physical
variation; silent eccentricity-clamp accuracy drop.

**Action**: file is FLAGGED for C-fail per Theorem 2.1 of the framework.
Downstream consumers cannot trust `total_error()` from
`propagate_deep_space` as a rigorous bound until DPPER/DSPER/time-varying
contributions are either implemented or explicitly accounted for in the
accuracy budget.
