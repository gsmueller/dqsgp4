# Explicit Runge–Kutta on SE(3): the Butcher family and the unified `rk_step` driver

**Theory note (FIRST, before code) for survey item P5** — the integrator unification. Governed by
[[feedback_theory_first_library]] (exhaustive theory precedes code; the functions *fall out*) and
[[feedback_no_perceived_fidelity]] (any value change is gated against an independent oracle to a stated
tolerance). This note derives the explicit Runge–Kutta family, its Lie-group (Munthe-Kaas) adaptation to the
second-order rigid-body system on SE(3), and shows that the three State-space integrators the library ships —
`euler`, `runge_kutta_4`, `rkf78` — are **one driver evaluated at three Butcher tableaux**. It then disposes,
honestly, of the two integrators that are *not* members of this family (`attitude_rk4_step`, `symplectic_leapfrog`).

The companion code is `src/integrators/runge_kutta.h` (primitives + `ButcherTableau` + `rk_step` + euler/rk4) and
`src/integrators/runge_kutta_fehlberg.h` (rkf78). The gate is `tests/test_butcher_tableau` plus the standing
physics gates `test_integrator_order`, `test_integrator_rkf78`, `test_propagator`, `test_integrator_symplectic`.

---

## 1. The duplication this resolves (survey P5)

The architecture-reuse survey (`design/ARCHITECTURE_REUSE_SURVEY.md` §P3) flagged that *"RK4 is written twice"* and
that *"the rkf78 loop already is this driver with a hardcoded tableau."* Reading the three State-space integrators
confirms it precisely:

- `rkf78_step` (runge_kutta_fehlberg.h) is **already the generic driver**: it loops over a 13×13 stage matrix
  `a`, node vector `c`, and weight vector `b` held in a `Rkf78Tableau<T>` struct.
- `runge_kutta_4` (runge_kutta.h) is the **same loop, hand-unrolled to four stages** with the RK4 tableau inlined.
- `euler` (runge_kutta.h) is the **one-stage** member of the same loop.

These three are the *same arithmetic re-implemented at ≥2 sites* — the binding validity criterion for a
unification (the corrected formula plan; "unify for uniformity" is churn). So the legitimate dedup is: lift the
stage loop out of `rkf78_step` into a `ButcherTableau<T,S>` + `rk_step<T,S>` driver, and re-express euler, rk4 and
rkf78 as that driver applied to their tableaux.

Two other integrators are **deliberately excluded** — they are *different algorithms*, not copies of this loop;
see §7. The boot framing called both RK4s "Munthe-Kaas schemes on different manifolds"; that is imprecise and §7
corrects it (the State RK4 is Munthe-Kaas; the attitude RK4 is ambient-then-retract).

---

## 2. Explicit Runge–Kutta on a vector space (the classical theory)

For an ODE `dy/dt = f(t, y)` on a vector space, with state `y₀` at time `t₀` and step `h`, an explicit
`S`-stage Runge–Kutta method is defined by a **Butcher tableau** — nodes `cᵢ`, a strictly-lower-triangular stage
matrix `aᵢⱼ` (j < i), and weights `bᵢ`:

```
        k_i = f( t₀ + c_i h,  y₀ + h Σ_{j<i} a_ij k_j )          i = 0 … S−1
        y₁  = y₀ + h Σ_i b_i k_i
```

`k_i` are the **stage derivatives**. Strict lower triangularity ⇒ each `k_i` depends only on earlier stages ⇒
*explicit* (no solve). Two structural identities every consistent tableau satisfies:

- **Row-sum (node consistency):** `Σ_{j<i} a_ij = c_i`. The node `c_i` is *defined* as where stage i samples
  the independent variable; the row sum must agree so the stage state is a consistent first-order extrapolation.
- **Weight normalisation:** `Σ_i b_i = 1` (the method reproduces the exact integral of `f ≡ const`).

**Order.** A method has order `p` when its local truncation error (LTE) over one step is `O(h^{p+1})`. Order `p`
is equivalent to the **quadrature/order conditions** holding up to degree `p`; the first family (sufficient to
*check* a claimed order against a transcription error) is

```
        Σ_i b_i c_i^{k} = 1/(k+1)        k = 0 … p−1.
```

(`k=0` is `Σ b_i = 1`.) Higher trees (e.g. `Σ b_i a_ij c_j = 1/6`) are needed for `p ≥ 3` but the quadrature
family already pins every coefficient enough to catch a typo — which is exactly how `test_butcher_tableau` and
`test_integrator_rkf78` validate the tableaux: **exactly, from the rationals, with no numerical reference.**

The three tableaux:

| name   | S  | `c`                                   | `b`                                                | order `p` |
|--------|----|---------------------------------------|----------------------------------------------------|-----------|
| euler  | 1  | (0)                                   | (1)                                                | 1         |
| rk4    | 4  | (0, ½, ½, 1)                          | (⅙, ⅓, ⅓, ⅙)                                       | 4         |
| rkf78  | 13 | Fehlberg 1968 (NASA TR R-287)         | 8th-order weights; embedded 7th for error control  | 8 (7)     |

rk4's stage matrix is `a₁₀=½, a₂₁=½, a₃₂=1` (all else 0); euler's is the 1×1 zero. The rkf78 coefficients are the
exact rationals already transcribed in `Rkf78Tableau` (verified there by the order conditions through `k=7`).

---

## 3. The system here is second-order on SE(3) — why naive RK is wrong

The propagator state (`dynamics/state.h`, `derivative.h`) is **not** a flat vector. It is

```
        pose   M̂ ∈ SE(3)        (unit dual quaternion: rotation ⊕ translation)
        twist  Ω̂ ∈ se(3) ≅ ℝ⁶   (body angular ω ⊕ body linear v)
        time   t ∈ ℝ
```

evolving under the **partitioned, second-order** law

```
        dM̂/dt = ½ M̂ · Ω̂_pure       (kinematics: configuration advanced by velocity — LIE-GROUP)
        dΩ̂/dt = a_body              (dynamics: velocity advanced by acceleration — VECTOR SPACE)
        dt/dt = 1.
```

The velocity `Ω̂` and the acceleration `a_body` live in the *vector space* `se(3)` (twists add component-wise —
`Twist::operator+`). But the configuration `M̂` lives on the *curved manifold* SE(3): `M̂₁ + M̂₂` is **not** a
rigid motion. Adding a scaled derivative to a unit dual quaternion (the flat-RK update `y₀ + h Σ b_i k_i`) leaves
the unit-norm manifold — the result is not a valid pose. Naive RK on `M̂` is therefore ill-defined; it only
"works" if you renormalise afterwards, and even then the *interior* stage states are off-manifold, corrupting
every `f` evaluation.

**Munthe-Kaas / Crouch–Grossman.** The fix is to advance the configuration by the **exponential map** of an
algebra element, never by ambient addition. For a rigid body the exponential is the **screw motion**
`exp_screw : se(3) → SE(3)` (`DualQuaternion::exp_screw`, half-angle convention REQ-DQ-18). The configuration
update from a base pose `M̂₀` by a twist increment `ξ ∈ se(3)` over time `h` is the right-translation

```
        advance(M̂₀, h, ξ) = ( M̂₀ · exp_screw( (h/2) ξ ) ).normalized()      ← lie_advance_pose
```

— the `½` is the half-angle of `exp_screw`; `.normalized()` is the SE(3) retraction (REQ-EF-15) that projects
sub-ULP drift back onto the unit-DQ manifold. This is the single primitive on which the whole family is built.

### 3.1 Why the velocity-combination can be exponentiated once (the key simplification)

A general Munthe-Kaas method carries the dexpinv correction because successive exponentials from *moving* base
points do not compose additively. Here the structure is simpler, and it is worth stating exactly *why* a single
exponential of the combined twist is correct, because it is what makes the driver clean:

1. Every stage and the final update advance the pose from the **same** base `M̂₀` (not from the previous stage's
   pose). So we never compose two exponentials about different base points within a step.
2. The pose-advance depends on the configuration only through the **twist** (the kinematic law `dM̂/dt = ½M̂Ω̂` is
   right-invariant: the increment is `M̂₀ · exp(·)`). The twist itself lives in the flat `se(3)`, where the RK
   stage combination `Σ a_ij W_j` is ordinary vector arithmetic.

Hence for each stage `i` the *configuration* update is `advance(M̂₀, h, Σ_{j<i} a_ij W_j)` — one exponential of an
ordinary `se(3)` combination — and the *velocity* update is the flat RK rule `W_i = Ω̂₀ + h Σ_{j<i} a_ij A_j`.
The pose stays exactly on SE(3) at **every** stage (each is `M̂₀ · exp(algebra)`), so every `f` evaluation sees a
valid rigid motion. This is the precise sense in which "the State RK4 is a Munthe-Kaas Lie-group scheme."

---

## 4. The unified driver — the functions fall out

Writing §2 with §3's two substitutions (vector RK on the velocity `W`/acceleration `A`; exponential advance on
the pose) gives the driver directly. With `A_i` the body acceleration returned by the callback at stage `i`:

```
    stage loop  (i = 0 … S−1):
        W_i        = Ω̂₀ + h · Σ_{j<i} a_ij A_j               (twist:  flat se(3) RK)
        stagePose  = advance( M̂₀, h, Σ_{j<i} a_ij W_j )       (pose:   single screw-exp from M̂₀)
        A_i        = accel_fn( State{ stagePose, W_i, t₀ + c_i h } )

    combine:
        M̂₁ = advance( M̂₀, h, Σ_i b_i W_i )                    (pose:   screw-exp of b-weighted twists)
        Ω̂₁ = Ω̂₀ + h · Σ_i b_i A_i                             (twist:  b-weighted accelerations)
        t₁ = t₀ + h.
```

This is **exactly** the loop currently inside `rkf78_step`. So the abstraction that falls out is:

```cpp
template<typename T, int S>
struct ButcherTableau {                 // c, b strictly defined; a strictly lower triangular; exact rationals
    std::array<TrackedValue<T>, S> c, b;
    std::array<std::array<TrackedValue<T>, S>, S> a;
};

template<typename T, int S>
struct RkStepResult {                   // the integrated state + the stage accelerations (for the LTE policy)
    State<T> state;
    std::array<Twist<T>, S> stage_accel;
};

template<typename T, int S>
RkStepResult<T,S> rk_step(const ButcherTableau<T,S>&, const State<T>& y0,
                          const TrackedValue<T>& dt, const AccelFn<T>&);
```

and the three entry points are one line of tableau selection each:

- `euler(y0, dt, f)`         = `rk_step(euler_tableau, …)` + Euler LTE deposit.
- `runge_kutta_4(y0, dt, f)` = `rk_step(rk4_tableau, …)`  + RK4 LTE deposit.
- `rkf78_step(y0, dt, f)`    = `rk_step(rkf78_tableau, …)` + the embedded 7(8) error.

`rk_step` performs **pure integration only** (no error book-keeping). Each entry point applies its *own* error
policy afterward (§6) — keeping the load-bearing precision-vs-accuracy choice explicit and auditable at each
call site, per the constants-initiative discipline ([[feedback_constants_generative_or_bounded]]).

### 4.1 euler and rk4 as specialisations — worked

- **euler (S=1):** `c=(0), b=(1), a=((0))`. Stage 0: `Σ_{j<0}=∅` ⇒ `W₀ = Ω̂₀ + h·0 = Ω̂₀`, `stagePose =
  advance(M̂₀, h, 0)`, `A₀ = f(stagePose, Ω̂₀, t₀)`. Combine: `M̂₁ = advance(M̂₀, h, 1·Ω̂₀)`,
  `Ω̂₁ = Ω̂₀ + h·A₀`. This is forward Euler: pose screwed by `h·Ω̂₀`, twist kicked by `h·a`. ✔
- **rk4 (S=4):** the four stages reduce to the classical `½,½,1` cascade. Note `advance(M̂₀, h, ½W)` and
  `advance(M̂₀, h/2, W)` are *identical* — both screw by `(h/4)W` — because `exp_screw` is linear in `(time ×
  twist)`. So the generic `Σ a_ij W_j` form reproduces the hand-rolled `half_dt`-scaled form. The combine
  `Σ b_i W_i = ⅙W₀+⅓W₁+⅓W₂+⅙W₃ = ⅙(W₀+2W₁+2W₂+W₃)` reproduces `twist_avg`. ✔

### 4.2 The one numerically-visible change (round-off, not bit)

The hand-rolled `runge_kutta_4` factors the combine as `⅙·(W₀ + 2W₁ + 2W₂ + W₃)` (one rational `⅙`, exact
integer `2`). The generic driver evaluates `Σ_i b_i W_i` as a fold with each `b_i` a separate rational
(`⅙, ⅓, ⅓, ⅙`). These are equal in exact arithmetic but **differ in the last ULP** in floating point (different
rounding of `⅓` vs `2/6`). Likewise the generic stage-0 pose is `advance(M̂₀, h, 0) = M̂₀.normalized()`, whereas
hand-rolled rk4 samples `f` at `M̂₀` un-renormalised; for a unit pose these differ only at round-off (and rkf78
*already* renormalises stage 0 this way). Therefore:

- **rkf78**: bit-identical (its loop *is* `rk_step` — the same code, merely factored out).
- **euler**: bit-identical up to the stage-0 retraction (`exact(1)·W` and `+h·0` are IEEE-exact identities).
- **rk4**: round-off-different (the combine reorders) — `O(ε_machine)` per step.

This is a **round-off-gated** change, not bit-exact. It is admissible because the only consumer of `runge_kutta_4`
off the frozen path is the DQ propagator, whose tests are round-off-tolerant (energy/Lz drift `1e-6` rel, pos
closure `1 km`); the tight `1e-12`/`1e-9` gates call the force functions *directly*, not through the integrator.
The frozen SGP4/SDP4 path is analytic and **never** uses these integrators, so OR1 and `test_sgp4` (33/33) are
untouched. See §8.

---

## 5. The `ButcherTableau` type and the tableau factories

`ButcherTableau<T,S>` zero-initialises `a`, `b`, `c` to `exact<T>(0)`; the factories set the nonzero entries as
exact rationals (`ratio<T>`) / exact integers (`exact<T>`) — born-digital, W18-clean (no float literals), the
same provenance discipline `Rkf78Tableau` already used. `euler_tableau<T>()` and `rk4_tableau<T>()` live with
the euler/rk4 entry points in `runge_kutta.h`; `rkf78_tableau<T>()` (the 13-stage Fehlberg coefficients, lifted
verbatim from the former `Rkf78Tableau` constructor) lives with rkf78 in `runge_kutta_fehlberg.h`. The exact
rationals are **the** correctness statement — `test_butcher_tableau` checks all three against the order
conditions of §2 with no numerical reference, so any transcription error is caught exactly.

---

## 6. Error policies — the PER-SLOT truncation deposit (V3, 2026-06-10; supersedes the uniform stamp)

### 6.1 Why the uniform stamp had to go (measured)

Until V3, euler/rk4/leapfrog deposited ONE magnitude — `C_p·h^{p+1}·maxᵢ‖Aᵢ‖₁` over the full 6-component
twist — uniformly on **all 14 state slots** (`State::add_uniform_accuracy`). The slots are dimensionally
heterogeneous: the pose's real part is a UNIT quaternion (components O(1)); its dual part is ½·t⊗q_r
(components O(|r|/2)); the twist is velocity-scale. Stamping a translational-scale magnitude on the O(1)
real slots is not merely loose — position extraction `r = 2·q_d⊗q_r*` multiplies the real-slot bound by
|q_d| ≈ |r|/2, so the extracted position's accuracy is amplified by POSITION SCALE. Measured (2026-06-10):
one 30 s RK4 step on a LEO state deposited 1.85e6 on every slot and `position()` reported **4.3e13 m**
of accuracy; an hour of steps saturated the channel to `inf` while values, measurement, and precision all
stayed correct. For an "infinite accuracy upscale" solver, the channel a user dials must stay informative;
the uniform stamp structurally prevented that.

### 6.2 The magnitudes and the slot map (the derivation)

**Magnitudes — the dimensional envelope.** The only in-step datum is the stage-acceleration bound,
split by part:

    A_lin = maxᵢ ‖Aᵢ.linear‖₁  [m/s²] ,   A_ang = maxᵢ ‖Aᵢ.angular‖₁  [rad/s²]

A true order-p LTE needs sup|y^{(p+1)}|, which |a| cannot supply; the old `C_p·h^{p+1}·A` form bought
its h^{p+1} scaling by abandoning units (m·s^{p-1} numbers stamped on m/s and unit-quaternion slots),
and — measured during V3 bring-up — the m·s³-as-m/s velocity bound COMPOUNDS through the very next
step's pose update (the tracked arithmetic faithfully propagates err_v into the pose), re-saturating
extraction to inf even with a correct slot map. The honest envelope treats h itself as the only
available in-step timescale (the derivative-cascade proxy |y^{(k+1)}| ≲ A/h^{k-1}), which collapses
every fixed-step method's per-step LTE to the SAME dimensional bounds:

    P = (h²/2)·A_lin  [m]      Θ = (h²/2)·A_ang  [rad]      (pose: position / rotation)
    V = h·A_lin       [m/s]    Ω = h·A_ang       [rad/s]    (twist: linear / angular)

These are deliberately method-order-BLIND: with |a| as the only datum, claiming an h^{p+1} bound would
be perceived accuracy. The higher-order methods' true superiority lives where it is MEASURED — the
value-based gates (RK1's empirical orders 1/4/7.996; the energy/closure tolerances), and the adaptive
path, whose embedded 7(8) difference IS an order-true per-step estimate. The envelope is conservative
by orders of magnitude (RK1 measures rk4's real error ~mm where the envelope says ~km) but it is
finite, unit-correct, h-dialing (h² per step), and a genuine bound under the stated proxy.

**The long-arc horizon (measured, honest).** The velocity envelope accumulates h·A per step, so after
k ≈ |v|/(h·A) steps (LEO, h = 30 s: k ≈ 7546/276 ≈ 27, about 14 minutes) the rigorous velocity bound
exceeds the velocity itself; the tracked arithmetic then FAITHFULLY reports the consequences — 1/r²
in the geopotential and exp(Δh/H) in drag have genuinely unbounded images over such an input interval,
so the extracted accuracy diverges (the max-sentinel → inf). This is rigor, not defect: with |a| as
the only in-step datum, no fixed-step bound can certify a long arc. The order-true long-arc accuracy
channel is `propagate_adaptive` — the embedded 7(8) difference is a measured per-step error, orders
below the envelope, and its accumulation stays informative over arcs the envelope cannot certify.
Per-STEP, the envelope is the honest, finite, h²-dialing claim (gates W21 and EX2 pin both halves).

**The slot map.** Mapping P and Θ onto the dual-quaternion slots is small-perturbation algebra on the
unit constraint:

- **Real part.** A rotation error of angle Θ perturbs a unit quaternion by δq = q⊗[cos(Θ/2)−1, sin(Θ/2)·n̂];
  every component of δq is bounded by **Θ/2** (vector part ≤ sin(Θ/2) ≤ Θ/2; scalar part ≤ Θ²/8 ≤ Θ/2).
  → real slots receive Θ/2.
- **Dual part.** q_d = ½·t⊗q_r, so δq_d = ½(δt⊗q_r + t⊗δq_r). With ‖q_r‖ = 1, |t| = 2‖q_d‖, and the
  quaternion-product component bound ‖a⊗b‖∞ ≤ ‖a‖₁·‖b‖∞ folded conservatively into magnitudes:
  ‖δq_d‖ ≤ **P/2 + ‖q_d‖₁·Θ/2**. → dual slots receive P/2 + ‖q_d‖₁·Θ/2 (the cross term vanishes when Θ = 0).
- **Twist.** The linear slots receive V, the angular slots Ω — the parts no longer cross-contaminate,
  and the velocity bound is a velocity (m/s), so its propagation into subsequent pose updates stays
  h-scaled instead of compounding a mis-united magnitude.

Consequence for extraction: δr from `r = 2·q_d⊗q_r*` is now ≈ 2(‖δq_d‖ + ‖q_d‖·‖δq_r‖) ≈ P + 2‖q_d‖₁·Θ —
the translational LTE envelope itself (plus the honest lever-arm term when rotating), with NO
position-scale amplification of a translational deposit. Torque-free orbits get real-slot accuracy
EXACTLY 0 — the true statement the uniform stamp could not make: there is no rotational truncation when
nothing rotates.

### 6.3 The policies after V3

- **euler / rk4 / leapfrog** — the SAME dimensional envelope (P, Θ, V, Ω above; A from each method's own
  stages: A₀ for euler; maxᵢ over the four rk4 stages; max(A₀, A₁) for leapfrog), deposited via
  `State::add_step_lte(P, Θ, V, Ω)`. Method order deliberately does not enter the bound (§6.2); it is
  measured by RK1 and realized in the values.
- **rkf78** — UNCHANGED: the *embedded* 7(8) estimate `(41/840)·h·‖A₀ + A₁₀ − A₁₁ − A₁₂‖` deposited to the
  6 twist components only (a velocity-space, order-true estimate; the asymmetry is pre-existing and
  intentional — adaptive stepping is the order-true accuracy path).

`Twist::l1_norm_linear()/l1_norm_angular()` are the split proxies (l1_norm = their sum remains for LTE-free
callers); `State::add_step_lte` replaces `add_uniform_accuracy` (whose only callers were these integrators).
Gate W21 asserts the slot law directly (spin-free: real slots exactly 0, dual = P/2 = (h/4)·V, twist.linear
= V); the EX2 budget section measures the extracted accuracy FINITE and h-dialing where it previously
saturated to inf.

---

## 7. What is NOT in this family (honest dispositions)

The validity criterion forbids folding in code that is not *the same arithmetic*. Two integrators look adjacent
but are different algorithms; folding them would be churn or a silent method-change.

### 7.1 `attitude_rk4_step` (dynamics/attitude.h) — a *different* RK4, kept separate

It integrates `AttitudeState = (q ∈ SO(3), ω ∈ ℝ³)` under `q̇ = ½q⊗ω`, `ω̇ = I⁻¹(τ − ω×Iω)`. Crucially it is
**not** Munthe-Kaas: it treats `(q, ω)` as one flat 7-vector, does ordinary ambient RK4 with component-wise
quaternion addition `q + Δq` at *every* stage, and retracts **once** at the end (`next.q.normalized()`). During
the interior stages `q` leaves the unit sphere — the opposite of §3's per-stage on-manifold guarantee. So:

- **Different state type** (`AttitudeState`, not `State`), **different scheme** (ambient-then-retract vs per-stage
  Lie-advance), and a **single consumer** (`test_attitude_dynamics`, the H1 coupled-Euler check).
- It shares with the State family *only the Butcher numbers* `(⅙,⅓,⅓,⅙)` — not the arithmetic. Forcing it under
  `rk_step` would require either (a) a state-concept abstraction parameterised by retraction strategy with
  exactly one extra, differently-behaving consumer (churn), or (b) rewriting it as a Munthe-Kaas scheme on SO(3)
  — a genuine *value change* needing its own fidelity justification and gate, well beyond a dedup.

**Disposition:** left as-is. The survey's "AttitudeState as a `State` view (deletes the duplicate RK4)" is a
separate P5 *state-adapter* item — make `AttitudeState` the rotational sub-object of `State` and adopt the MK
scheme — not part of this ButcherTableau dedup. Recorded so the boundary is explicit, not forgotten. The
`rk_step` driver is written concretely over `State<T>` (matching the survey's `rk_step(ButcherTableau, State, …)`
signature); the state-concept generalisation is **deferred until a second Munthe-Kaas consumer exists**
(anti-dead-code — the same discipline that shipped only the consumed subset in L2/L3).

### 7.2 `symplectic_leapfrog` (integrators/symplectic_leapfrog.h) — a splitting method, kept separate

Velocity-Verlet (kick–drift–kick) is a **partitioned/splitting** integrator, not an explicit-RK Butcher method
in the sense of §2 (its half-kick → drift-with-force-re-evaluation → half-kick structure exploits the H = T(p) +
V(q) split for symplecticity and the resulting absence of secular energy drift). It is not a copy of the §4
loop and cannot be expressed as one without losing the structure that is its entire purpose. Left separate. (Its
pose-LTE under-bounding — it formerly deposited the O(dt³) bound to the 6 twist slots only, not the 8 pose slots,
unlike euler/rk4 — was fixed in the follow-up commit: the drift advances the pose, so the bound now deposits
uniformly across all 14 components via `State::add_uniform_accuracy`, gated by `test_integrator_symplectic`
(task_731830b8).)

---

## 8. Validation & the no-perceived-fidelity analysis

This increment **adds no fidelity** — it re-expresses three existing methods through one driver. So the oracle
question is conformance, not a new physical claim. Three independent checks:

1. **Tableau correctness — exact, no numerical reference.** All three tableaux satisfy `Σ_{j} a_ij = c_i`,
   `Σ b_i = 1`, and `Σ b_i c_i^k = 1/(k+1)` to their claimed order, checked from the rationals
   (`test_butcher_tableau`; rkf78 also in `test_integrator_rkf78`). A transcription error breaks one identity.
2. **Empirical order — independent physical oracle.** The closed-form circular Kepler orbit `r(t) = R(cos ωt,
   sin ωt, 0)` is exact truth; `test_integrator_order` measures `p = log₂(E(n)/E(2n))` and asserts euler `~1`,
   rk4 `~4`, rkf78 `≥5`. These ratios are round-off-insensitive, so they survive the §4.2 combine reorder.
3. **Round-off equivalence to the prior implementations** — established standalone (vcvars64 + `cl`) before
   committing: rkf78 bit-identical, euler bit-identical up to stage-0 retraction, rk4 `O(ε)`. The standing
   physics gates (`test_propagator` closure `<1 km` / energy `<1e-6`; `test_integrator_symplectic`;
   `test_adaptive_stepping` `<1 m` vs a fine reference) then exercise the refactored entry points end-to-end and
   confirm no behavioural change beyond round-off.

**Frozen invariants preserved:** OR1 (`test_sgp4_regression`) 0 km and `test_sgp4` 33/33 — the SGP4/SDP4 path is
analytic and never routes through these integrators; the inline J2/J3/J4 + `sgp4::ZonalHarmonics` are untouched.
The full sweep (incl. the MATLAB W*.m derivation verifiers) is run at the phase boundary.

---

## 9. Summary — the abstraction

```
  AccelFn<T>          : State<T> → Twist<T>          (body acceleration callback; injected by the propagator)
  lie_advance_pose    : (Pose, h, Twist) → Pose      (M̂₀ · exp_screw((h/2)ξ), retracted — the SE(3) primitive)
  ButcherTableau<T,S> : { c, b, a }                  (exact-rational tableau; euler S=1, rk4 S=4, rkf78 S=13)
  rk_step<T,S>        : (tableau, y0, h, f) → {state, stage_accel}   (the generic Lie-group RK driver)
  euler / rk4 / rkf78 : rk_step(tableau) + per-method LTE/embedded-error policy
```

One driver, three tableaux, three error policies. `attitude_rk4_step` (ambient-then-retract on SO(3), single
consumer) and `symplectic_leapfrog` (splitting) are *not* members and stay separate — the honest boundary the
validity criterion draws.
