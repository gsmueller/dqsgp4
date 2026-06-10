# DQSGP4 vision review — the whole codebase against the stated intent (2026-06-10)

**Mandate (user, 2026-06-10, verbatim intent):** review all of the code against the corrected
framing — (1) the frozen SGP4 is *code to test against*, **not part of the final code base**;
(2) the DQSGP4 project is an **infinite precision/accuracy upscale of a dual-quaternion based
propagator, independent of Euler angles**.

**Verdict in one paragraph.** The mathematics already matches the vision: the dual-quaternion
state and screw-exponential integration are Euler-angle-free, the entire computational substrate
is templated on T with exactly one documented double-core exception, and accuracy is dial-up
almost everywhere. Three things do NOT match: (i) the **product surface and prose** present SGP4
as a co-equal shipped propagator (umbrella exports, "two independent propagators" in README/docs/
diagram) rather than as the test oracle + TLE-epoch adapter; (ii) the **propagated accuracy
channel saturates** (the uniform per-step LTE deposit — the one place the "infinite accuracy"
claim fails structurally); (iii) the **arbitrary-precision claim for the DQ propagator itself is
unmeasured** — no gate propagates `dynamics::Propagator<cpp_bin_float_50>`, while the SGP4 side
ironically has exactly that gate (DS1). Each is fixable; none requires moving the mathematics.

---

## 1. Lens 1 — SGP4 as test-against code, not product

### 1.1 The map (measured from the include graph, 2026-06-10)

The SGP4-theory tier is **≈ 4,986 of 18,026 src lines (27.7%)**:

| Tier piece | Files | Lines |
|---|---|---|
| `src/sgp4/` (selector, near/deep space, propagator, state vector) | 6 | 2,266 |
| `src/perturbation/` (Brouwer, Kaula, resonance, short-period, SR3 third-body) | 5 | 1,152 |
| `src/orbit/` SGP4 pipeline (element_recovery, modified_kepler, osculating, secular, state_from_elements) | 5 | 591 |
| `src/atmosphere/` Lane drag (density_model, drag_coefficients) | 2 | 329 |
| `src/astronomy/solar_system.h` (SR3 lunisolar elements) | 1 | 297 |
| `src/ephemeris/` SR3-historical trio (dispositioned R4a) | 3 | 351 |

**The coupling surface is already minimal.** Exactly three non-SGP4 files include `sgp4/`:

- `src/dqsgp4.h` — the umbrella **exports** `sgp4_propagator.h` + `model_selector.h` (product framing — the finding);
- `src/dynamics/state_from_tle.h` — the TLE→epoch-state seeder (load-bearing, see 1.2);
- `src/dynamics/state_conversion.h` — the `StateVector` parity adapters (F1; test-facing in spirit).

Everything else in the DQ core is SGP4-free. `orbit/kepler_series.h` (generative Fourier–Bessel)
is DQ-side and consumed by `ephemeris/sun_meeus.h`; the rest of `orbit/` is pipeline-internal to
SGP4. The Aoki-82 GMST (`sidereal_time.h`) is genuinely shared (the TEME convention itself).

### 1.2 The physics constraint a relocation decision must respect

A TLE's elements are **SGP4 mean elements by definition** — they have no physical meaning except
through the model that fitted them. Recovering an epoch Cartesian state requires evaluating the
SGP4 pipeline at t = 0 (the short-period corrections are km-scale for LEO; treating TLE mean
elements as osculating Kepler elements injects ~10 km-class epoch error). The t = 0 evaluation
exercises nearly the whole pipeline (init, element recovery, long-period, Kepler, short-period),
so there is no smaller "epoch-only" kernel to extract.

Therefore "SGP4 not in the final code base" has exactly two coherent realizations:

- **(A) Adapter + oracle (recommended).** SGP4 leaves the PRODUCT SURFACE (umbrella, docs
  identity, "two propagators" framing) but its templated kernel remains in-tree as (1) the
  internal engine of the TLE-epoch adapter (`state_from_tle`) and (2) the verification oracle
  (OR1 golden master, 33/33 SGP4-VER). TLE ingestion stays first-class. The kernel being
  TrackedValue<T>-templated means even the seed inherits arbitrary precision (DS1 measured this).
- **(B) Full test-side exile.** The library ingests Cartesian/osculating states only; the TLE→
  state conversion becomes a test-side utility. TLE input is demoted to a fixture format. This is
  clean ideologically but costs the dominant real-world data path.

Physical relocation of the ~5 k lines (a `src_oracle/` tree) was weighed and is NOT recommended
for the same reason the clean-sheet sln was rejected (plan §7): ~80 gated projects reference the
layout; the include-graph already proves containment; moving files is configuration-control risk
with no semantic gain. Containment should be made MECHANICAL instead (a gate — §5, T2).

### 1.3 The product-framing inventory (what says "two propagators" today)

- `README.md:3` "two independent propagators", `:8` "The two propagators — pick one";
- `docs/index.html` overview card + the `sgp4_vs_dqsgp4` diagram title ("One TLE, two
  propagators"); guide §1 ("The library ships two independent propagators");
- `dqsgp4.h` umbrella PROPAGATORS section listing SGP4 first and exporting it;
- the `Propagatable` concept + gate F2 framed as "API parity" between two products (the honest
  framing under the vision: the oracle exposes the same verb so tests and adapters can drive
  both — a test-facing property).

## 2. Lens 2 — infinite PRECISION upscale (audit)

**Substrate: clean.** Every computation path is templated on T over TrackedValue<T>; comparisons
in the integrator loops use `.value` generically; no bare-double members or locals exist in the
DQ core (grep-audited) outside the two documented cases:

- `atmosphere/nrlmsise00.h` (7 casts) — the DOUBLE-core MSIS wrapper, by design and documented:
  a ~6-digit empirical fit with 15% RMS; wider T would be perceived precision (nrlmsise00.md §3).
  This is an honest accuracy-floor exception, not a leak.
- `tle/tle_parser.h` — `TleData` stores parsed fields as double. The TLE format carries ≤ 8
  significant digits, so the double ceiling (~1e-16 rel) sits ~8 orders below the format's own
  quantization, which the tracked measurement budget already records. Immaterial; noted for
  completeness. (A string→T re-parse in `TleElements<T>` would be purist-tighter but changes
  nothing observable under the existing budgets.)

**Finding P1 (the gap): the arbitrary-precision claim for the DQ propagator is unmeasured.**
No gate instantiates `dynamics::Propagator<cpp_bin_float_50>` (or the facade) and propagates.
test_precision_scaling runs bf50 on constants/derivation chains; test_srp/test_third_body use
bf50 for formula identities only. Meanwhile the SGP4 side HAS the exact gate the vision calls
for — DS1 propagates deep-space satellites at bf50 and asserts value agreement + dramatic
precision tightening. The flagship property of the flagship propagator should be measured the
same way (no perceived fidelity): a gate that integrates an arc at double and at bf50, asserts
value agreement at the double-roundoff scale, and asserts the position-precision channel
tightens by tens of orders.

## 3. Lens 3 — infinite ACCURACY upscale (audit)

**Dial-up knobs exist across the stack** (each with a tracked truncation deposit): zonal degree +
tesseral order (geopotential / GravityField), nutation `n_terms` (678-term table, amplitude-
sorted), Poisson `n_terms` (lunar 60+60), Kepler series `K` (converges for all e < 1), obliquity/
eccentricity polynomial terms, integrator order (euler/RK4/RKF7(8)) + `dt` + adaptive `tol`,
force composition (presets/extra_forces). The fixed-accuracy floors that remain are honest model
floors (Vallado-8-4 table representativeness; MSIS 15% RMS; Meeus series grades), all measured
and gate-recorded.

**Finding A1 (the structural violation): the propagated state's accuracy channel saturates.**
The integrators deposit the per-step LTE bound (RK4: dt⁵/120·max|a|; euler: dt²/2·|a|) uniformly
on all 14 state slots (W21-asserted semantics) — position-scale magnitudes stamped onto O(1)
unit-quaternion components. `position()` then propagates those bounds through r = 2·dual⊗real*,
multiplying by position scale: ONE 30 s RK4 step yields 4.3e13 m of "accuracy"; an hour saturates
to inf. Values, measurement, and precision all stay correct and informative — but the accuracy
channel of the propagated state, the very thing an "infinite accuracy upscale" dials, is
non-informative beyond one step. Under the Q-phase this was documented + EX2-pinned as measured
semantics with the redesign filed as optional (chip task_28c4798e). **Under the vision it is core
work**: a per-slot, dimensionally consistent LTE (rotation slots get the angular LTE, dual slots
the translational LTE in metres, twist slots theirs) so extracted-position accuracy ≈ the actual
integrator LTE and shrinks with dt and with method order — i.e. dialable. Own theory increment
(extend runge_kutta_lie_group.md), W21/RK1/EX2 updates, value-channel-only change.

**RESOLVED (V3, 2026-06-10).** runge_kutta_lie_group.md §6 now derives the slot map (real ← Θ/2;
dual ← P/2 + ‖q_d‖₁·Θ/2; twist parts their own; each ≤ the old uniform stamp, so no bound weakens);
`State::add_step_lte` replaces the uniform stamp; euler/rk4/leapfrog deposit split-by-part
magnitudes (rkf78's twist-only embedded estimate unchanged, documented). W21 asserts the law
directly (torque-free: real slots EXACTLY 0; dual = P/2; twist = P); EX2 measures the extracted
accuracy finite and dt-dialing where it previously saturated to inf.

## 4. Lens 4 — independence from Euler angles (audit)

**Clean.** The state is a unit dual quaternion + se(3) twist; the pose advances by ONE screw
exponential per step (Munthe-Kaas form); the attitude alternative integrates the quaternion
ambient-then-retract. There is no Euler-angle state, no per-axis angle integration, no
roll/pitch/yaw extraction anywhere in src/ (grep-audited). Angle-built rotations appear in
exactly two places, both interfaces rather than dynamics:

- `orbit/state_from_elements.h` (SGP4 tier): perifocal→TEME via R₃(−Ω)R₁(−i)R₃(−u) — inherent
  to ELEMENT inputs (orbital elements are angles); not on the DQ integration path.
- `astronomy/frames.h` + the GCRS→ITRS chain: literature-defined axis-rotation sequences for
  frame I/O (IAU conventions), applied to OUTPUT states.

Two naming caveats so the audit isn't misread: `integrators::euler` is Euler's METHOD (first-
order RK), and `dynamics/attitude.h` implements Euler's rigid-body EQUATIONS (ω̇ = I⁻¹(τ−ω×Iω));
neither involves Euler angles.

## 5. Remediation tiers (decision requested)

- **T1 — surface & prose (cheap, value-free).** Umbrella stops exporting `sgp4/*` (the facade
  keeps its documented internal opt-in); README/index/guide/diagram reframe: *the DQSGP4
  propagator* as the product, *the frozen SGP4* as its verification oracle and TLE-epoch
  semantics provider; quickstart §1 reframed accordingly ("the oracle defines what a TLE means;
  the product is the facade"); `Propagatable`/F2 described as oracle-parity (test-facing).
  Touches: dqsgp4.h doc+includes, README, gen_docs/gen_diagrams prose, EX1/EX2 includes
  (ground_track/quickstart pull sgp4 via the umbrella today — they would include the oracle
  explicitly where they demonstrate it). Gates: rebuild + DOC1 regen + EX1/EX2.
- **T2 — mechanical containment (new gate ARCH1).** A ScriptGate pinning the measured topology:
  outside `src/sgp4/` and tests, only `dynamics/state_from_tle.h` and
  `dynamics/state_conversion.h` may include `sgp4/`; the oracle-tier module pages carry a tier
  banner ("test-against / adapter tier — not the product surface"). This makes the user's intent
  un-rottable without moving a file.
- **T3 — the vision items (code).** P1: the DQ bf50 propagation gate (the DS1 analogue for the
  product). A1: the per-slot LTE redesign promoted from optional chip to core increment.

**Recommended realization: (A) + T1 + T2 + T3.** Option (B) (full exile, Cartesian-only library)
is available but trades away first-class TLE ingestion; flag if that is the actual intent.

**USER DECISION (2026-06-10): option (A) confirmed, with an explicit retention guarantee — the
arbitrary-precision SGP4 propagator is NOT deleted (nothing in the sgp4 tier is removed); the
goal is demonstrable UNPEELABILITY of the DQ solver. Realized as: the umbrella exports the solver
surface only (the oracle reachable by explicit include), ARCH1 mechanically pins the containment
(only the named adapters include sgp4/ from outside the tier), and PREC1's translation unit
includes no sgp4/tle header at all — the compile itself is the unpeel demonstration.**

## 6. Clean bills (what was checked and NOT found)

No Euler-angle state or extraction; no bare-double computation in the DQ core beyond the two
documented cases; no SGP4 include outside the three named files; no hidden SGP4 dependence in
forces/integrators/ephemeris/frames; the DQ propagator runs entirely from Cartesian state
(`Propagator<T>` has no TLE/SGP4 knowledge — the dependency enters only through the seeding
adapter and the facade's factories). Suite state at review: 81/81 full sweep, OR1 0 km,
test_sgp4 33/33 623/623, docs DOC1-fresh, 0 undocumented public symbols.
