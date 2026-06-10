# The DQ propagator facade completes — force injection, presets, adaptive, adapters (R3)

**Design note (FIRST, before code) for replan item R3** (`design/PROFESSIONAL_LIBRARY_PLAN.md` §6) — the L5
capstone: everything L4 built becomes *reachable* through the facade, with the default model bit-unchanged.
Governed by [[feedback_theory_first_library]] (the API falls out of the stated contract) and
[[feedback_no_perceived_fidelity]] (every claim gated). Companion code: `dynamics/dq_sgp4_propagator.h`,
`dynamics/propagator.h`, `dynamics/state_from_tle.h`; gates: new FM1 (`test_force_presets`) + extensions to
AD1/F2; the whole standing suite is the bit-regression backstop.

---

## 1. The defect being fixed (the doc-lie)

`dq_sgp4_propagator.h` has promised since its creation: *"Drag and third-body forces can be layered on via the
explicit constructor."* The constructor takes **no force list** — `make_propagator` hardcodes the unified
geopotential as the only force. The gated L4 forces (`make_drag` + the R1 table, `make_third_body_force`,
`make_srp_force`) are facade-unreachable. R3a makes the comment true.

## 2. The injection contract

The underlying `dynamics::Propagator` already sums an arbitrary `std::vector<ForceFn<T>>` (wrench-additive,
REQ-PR-2). The facade exposes exactly that:

- The explicit constructor gains `std::vector<ForceFn<T>> extra_forces = {}` — **appended after** the
  gravitational core inside `make_propagator`.
- **Invariant (the regression claim): default-empty `extra_forces` is BIT-IDENTICAL to the prior facade** —
  the same force list, same order, same integrator. Existing callers recompile unchanged; G1/F2/INJ1/DS1 and
  every propagation gate re-verify it.

## 3. Named presets — an options struct, not an enum explosion

The additions compose (lunisolar × drag × SRP × …), so named enum presets would be combinatorial. The
gravitational core keeps the existing `PropagatorMode { Authentic, Boosted }` selector; the perturbation
additions become an options struct, each member opt-in:

```cpp
template<class T> struct DqForceOptions {
    bool lunisolar = false;                 // Sun + Moon third-body (TB1-gated forces)
    std::optional<TV> drag_B;               // C_d·A/m  → Vallado-8-4 table drag (ATM1)
    std::optional<TV> srp_cr_area_over_mass;// C_R·A/m  → cannonball SRP + shadow (SRP1)
};
// from_tle(td, tol, mode, dt_max, options) assembles the force list; default options = {} = today's model.
```

`drag_B` and `srp_cr_area_over_mass` are caller-supplied because they are *physical properties of the
spacecraft* — no honest default exists (the same reason `make_drag` takes B).

**The base epoch for the Sun-dependent forces** is derived inside the preset from the TLE epoch
(`TleElements::from_tle_data(td).epoch_jd`, tagged UTC — what a TLE is). Honesty on the time scale: the
ephemeris argument nominally wants TT; TT − UTC ≈ 69 s (2000s era). That offsets the Moon ≤ ~40″ and the Sun
≤ ~10″ of arc — i.e. ≤ ~2e-4 *relative* on perturbations that are themselves ~1e-7 of central gravity, a
~2e-11 contribution to the acceleration: far below every propagation tolerance (the same magnitude argument
that demoted nutation). Documented at the construction site; a real UTC→TT conversion arrives with the L1
leap-second increment when a consumer forces it.

## 4. Adaptive stepping through the facade (R3b)

`rkf78_propagate_adaptive` (AD1-gated) is standalone-only; `Propagator` exposes only fixed-cadence
`propagate_to`. The facade method falls out of composing what exists:

```cpp
Propagator<T>::propagate_adaptive(y0, t_target, dt_initial, tol, dt_min)
    = integrators::rkf78_propagate_adaptive(y0, t_target, dt_initial, [this] accel = compute_acceleration, tol, dt_min)
```

**Claim: bit-identical to calling the standalone loop with the same callback** (it *is* that call). Gate: AD1
extended with a facade-vs-standalone bit-equality check.

## 5. Adapter dedup + the `Propagatable` concept (R3c)

- `state_from_tle` re-implements `to_state`'s km→m + identity-attitude pack inline (the survey-P5 “written 3×”
  finding; the third site is `to_state_vector`, the inverse). It now **delegates** to `to_state(sv, 0)` — the
  same per-component `· exact(1000)` multiplies in the same order, so **bit-identical**; the entire seeded
  suite (G1/F2/INJ1/DS1/OR1-adjacent) is the regression gate.
- `Propagatable<P,T>`: the concept naming the verb the two propagators share —
  `p.propagate(TrackedValue<T> tsince_minutes)`. Its honest consumer is the F2 api-parity gate: a *generic*
  function templated on the concept drives both `sgp4::Propagator` and `DqSgp4Propagator` through one code
  path (the concept is consumed, not scaffolding).

## 6. Gates

| claim | gate |
|---|---|
| default options ≡ prior facade, bit-identical positions | FM1 (and the standing suite) |
| each preset toggles a real, sane-magnitude perturbation (LEO ordering: drag(B=0.5) > lunisolar > SRP(0.026)) | FM1 |
| arbitrary `extra_forces` injection works (the doc-lie fixed) | FM1 (counting lambda + ½at² displacement) |
| facade adaptive ≡ standalone adaptive, bit-exact | AD1 extension |
| `state_from_tle` ≡ `to_state` delegation | the seeded suite (bit-regression) |
| `Propagatable` consumed generically over both propagators | F2 extension |

Everything is additive; the SGP4 authentic path and the default DQ model are bit-frozen throughout.
