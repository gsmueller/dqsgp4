# Which propagator? — SGP4 vs DQSGP4

This library ships **two independent propagators**. They solve different problems, take different inputs,
and produce output in different frames and units. This guide helps you pick.

## TL;DR

- **Want authentic SGP4 output from a TLE?** Use **`sgp4::Propagator`**. It reproduces the official
  reference to machine precision (33/33 satellites, 623/623 points on the `SGP4-VER` suite).
- **Want to integrate your own force model, carry rigid-body attitude, or swap the integrator?** Use
  **`dynamics::Propagator`** (the dual-quaternion "DQSGP4" propagator).

They are **not** layered — neither calls the other. Choosing one does not involve the other.

## Side by side

| Dimension | `sgp4::Propagator<T>` (SGP4) | `dynamics::Propagator<T>` (DQSGP4) |
|---|---|---|
| **Method** | Analytical mean-element theory: secular rates + periodic corrections evaluated in closed form at each epoch | Numerical integration of the rigid-body equations of motion on SE(3), state as a dual quaternion |
| **Namespace / header** | `sgp4::` / `sgp4/sgp4_propagator.h` | `dynamics::` / `dynamics/propagator.h` |
| **Input** | A TLE → `tle::TleElements<T>` | An initial Cartesian `dynamics::State<T>` (pose + twist) |
| **Output frame / units** | TEME, **km** and **km·s⁻¹** (`StateVector<T>`) | Inertial Cartesian, **m** and **m·s⁻¹** (`State<T>`) |
| **Forces** | The fixed SGP4 model (WGS72 gravity, atmospheric drag via B\*, lunar-solar + tesseral resonance for deep space) — not user-editable | A user-supplied `std::vector<ForceFn<T>>`; ships `gravity_central`, generic zonal `gravity_zonal` (J₂…Jₙ from a shared `ZonalHarmonics` provider, `gravity_zonal.h`), and exponential `drag` |
| **Time stepping** | None — closed-form evaluation at any time-since-epoch | A user-supplied `IntegratorFn<T>`; ships `euler` and `runge_kutta_4` (RK4 on the Lie group) |
| **Attitude** | None (point in space) | Full SE(3) orientation + body angular velocity |
| **Cost** | One closed-form evaluation per time | One integration from epoch to the target time (step-count × forces) |
| **Validation** | **33/33** SGP4-VER satellites, residuals ≤1e-9 km near-earth, ≤1.5e-7 km deep-space | LEO smoke test (`test_propagator`): energy + angular-momentum conservation over one orbit; central + J₂ + drag |
| **Maturity** | Reference-grade, fully coefficient-audited | General-purpose engine; `ConstantsProvider` delegates to the shared geodesy earth-model factories and a generative `ZonalHarmonics` provider |

## Decision tree

1. **Do you have a TLE and want the standard SGP4 answer?** → `sgp4::Propagator`. Stop here.
2. **Do you need a force or perturbation SGP4 does not model, or rigid-body attitude, or to control the
   integrator?** → `dynamics::Propagator`.
3. **Do you want to start DQSGP4 from a TLE-derived state** (e.g. to compare the two, or to continue a
   TLE orbit under a custom force model)? → use the bridge `dynamics::state_from_tle(...)` to produce the
   initial `State<T>`, then drive `dynamics::Propagator` with your forces.

## Why two propagators?

SGP4 is a *specific, frozen* analytical theory whose entire value is bit-faithful reproduction of the
catalog. You cannot change its physics without ceasing to be SGP4. The DQSGP4 propagator is the opposite:
a *general* numerical engine whose value is flexibility — inject any forces, any integrator, any inertia,
and carry full attitude. Keeping them separate lets each be exactly what it should be: SGP4 stays a pure
faithful assembly of the reference; DQSGP4 stays an unconstrained experimentation platform. Both share the
same `TrackedValue<T>` three-error framework, so an error budget computed in one is comparable to the
other.

## The bridge (TLE → DQSGP4 initial state)

`dynamics::state_from_tle(td, tolerance)` produces a `dynamics::State<T>` at epoch from a parsed TLE: it
recovers the epoch position/velocity (in metres, inertial) and packs them into the dual-quaternion state
with identity attitude and zero body rate. From there you propagate with your own force list — note that
the resulting trajectory follows *your* forces, not SGP4's, so it will diverge from the SGP4 answer unless
your forces match SGP4's model.

## The easy path: `DqSgp4Propagator` (TLE → DQSGP4, one line)

For the common case — seed DQSGP4 from a TLE and propagate under a gravitational model —
`dynamics::DqSgp4Propagator<T>` wraps the bridge, a default force model, and an integrator behind the *same
verb* as SGP4, so the two propagators feel identical apart from their mathematics:

```cpp
auto p = dynamics::DqSgp4Propagator<T>::authentic(td, tol, dt_max); // one line
dynamics::State<T> s = p.propagate(tsince_minutes);                 // SGP4's verb
```

It takes minutes-since-epoch (the SGP4 convention) and returns a `State<T>`; use the F1 converter
`dynamics::to_state_vector(s)` for a km/TEME `StateVector` directly comparable to SGP4's output.

### Authentic vs boosted modes

A TLE is WGS72 by definition, so the epoch state is **always** seeded from the authentic WGS72 SGP4. What
the mode selects is the *propagation* model:

| Mode | Ellipsoid | Zonal field | Use |
|---|---|---|---|
| `authentic` | WGS72 | WGS72 J₂…J₄ | the DQ analogue of `sgp4_standard` — "SGP4 as advertised" |
| `boosted` | WGS84 | EGM2008 J₂…J₉ | a higher-fidelity continuation of the same orbit |

`DqSgp4Propagator::from_tle(td, tol, PropagatorMode::Boosted, dt)` selects the mode by name. The
authentic-seed / boosted-propagation split is explicit, not silent. Widening the numeric type `T` (e.g.
`boost::multiprecision::cpp_bin_float_50`) is an orthogonal precision boost: the constants are generative,
so their error budgets tighten with `T`, while the authentic `double` path stays bit-faithful to the
reference (guarded by the frozen-SGP4 regression oracle).

## Input formats

Both propagators accept the same `tle::TleData`, populated from either format:

- **TLE** (two/three-line) → `tle::parse(line1, line2, td)`. Classic and **Alpha-5** catalog numbers are
  decoded to `td.catalog_number`; the column-68 checksums are validated and reported.
- **OMM** (CCSDS Orbit Mean-elements Message, KVN — the modern Celestrak/Space-Track format) →
  `tle::parse_omm_kvn(text, td)`. It maps to the same `TleData`, so it drives either propagator unchanged.
