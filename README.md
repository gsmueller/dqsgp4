# DQSGP4

DQSGP4 is a C++20 satellite orbit propagator built on dual-quaternion rigid-body
dynamics. The state of the spacecraft — orientation and position together — is carried
as a unit dual quaternion on SE(3) and advanced by screw-motion integration, so the
attitude is represented by a quaternion throughout and the formulation has no
gimbal-lock singularities. Force models for the geopotential, atmospheric drag, lunar
and solar gravity, and solar radiation pressure plug into the integration as
independent components, and each has been verified against an independent reference.

Three properties distinguish the library:

- **Arbitrary precision.** Every computation is templated on the numeric type. The same
  code runs in `double` or in 50-digit `boost::multiprecision::cpp_bin_float_50`, and
  the reported precision budget tightens accordingly.
- **Tracked error budgets.** Every computed value is a `math::TrackedValue<T>` carrying
  three separate error bounds — measurement uncertainty from the input data, numerical
  precision from finite arithmetic, and model accuracy from truncated theory — and the
  bounds propagate through every arithmetic operation.
- **Adjustable fidelity.** Accuracy scales with the gravity-field degree and order, the
  ephemeris series length, the set of enabled force models, and the integrator order
  and step size.

**Documentation:** <https://gsmueller.github.io/dqsgp4/> — a task-oriented guide,
the complete API index, per-module pages with theory, usage, and test coverage, and
generated diagrams.

The library also contains a complete, independently verified implementation of the
analytical SGP4/SDP4 theory. It serves two purposes. First, it is the reference that
the numerical propagator is tested against: the test suite reproduces the published
SGP4 verification set (33 satellites, 623 evaluation points) and pins the results
bit-for-bit. Second, it recovers initial states from two-line element sets, whose
orbital elements are defined in terms of the SGP4 model and would lose kilometres of
accuracy if treated as osculating Kepler elements. The analytical implementation is
kept architecturally separate: only two adapter files connect it to the numerical
propagator, and an automated check in the test suite enforces that separation.

|  | Numerical propagator (DQSGP4) | Analytical reference (SGP4/SDP4) |
|---|---|---|
| Entry points | `dynamics::DqSgp4Propagator<T>`, `dynamics::Propagator<T>` | `sgp4::Propagator<T>` |
| Model | Dual-quaternion integration on SE(3) with pluggable force models | Mean-element analytical theory (Spacetrack Report #3) |
| Input | A Cartesian state, or a TLE/OMM through the adapters | A TLE |
| Output | SE(3) pose and body twist, in metres and m/s | TEME position and velocity, in kilometres and km/s |
| Fidelity | Selectable: gravity degree/order, lunisolar, drag, radiation pressure, integrator | Fixed by the published theory |

See [`design/propagator_choice.md`](design/propagator_choice.md) for guidance on
choosing between them. Note that the two sides use different units; convert when
comparing results.

## Propagating a TLE with the analytical model

```cpp
#include "sgp4/sgp4_propagator.h"   // not part of the umbrella header; include explicitly
using T = double;

// 1. Parse the two TLE lines.
tle::TleData td;
tle::parse(line1, line2, td);

// 2. Pick a model. The WGS72 "sgp4_standard" preset reproduces the published reference.
sgp4::ModelConfiguration<T> config = sgp4::ModelSelector<T>::select("sgp4_standard", T(1e-12));
tle::TleElements<T>          elems  = tle::TleElements<T>::from_tle_data(td);

// 3. Build the propagator once, then propagate to any time since epoch (minutes).
sgp4::Propagator<T> prop(config, elems, T(1e-12));
sgp4::StateVector<T> sv = prop.propagate(math::TrackedValue<T>::exact_integer(0));

// sv.position_km.{x,y,z} and sv.velocity_km_s.{x,y,z} are TEME components; every
// component carries .value and .errors.{measurement, precision, accuracy}.
double x_km = sv.position_km.x.value;
```

Model presets: `"sgp4_standard"` (WGS72, the reference), `"sgp4_wgs72_old"`,
`"sgp4_wgs84"`, `"modern_2020"`, and `"research_full"`. A custom configuration can be
assembled with
`sgp4::ModelSelector<T>::custom().gravity(...).astronomy(...).kepler(...).build(tol)`.

## Propagating with the numerical model

The one-line path from a TLE:

```cpp
auto p = dynamics::DqSgp4Propagator<T>::authentic(td, tol, dt_max);
dynamics::State<T> s = p.propagate(tsince_min);   // same call convention as the analytical side
```

The `authentic` mode propagates under WGS72 gravity to match the analytical model's
assumptions; `boosted` uses WGS84 with EGM2008 zonal harmonics. Perturbations beyond
gravity are enabled per propagator, and the default (gravity-only) model is unchanged
when none are requested:

```cpp
dynamics::DqForceOptions<T> opt;
opt.lunisolar = true;                     // Sun and Moon third-body attraction
opt.drag_B = tv(0.01);                    // drag with C_d·A/m, static-table atmosphere
opt.srp_cr_area_over_mass = tv(0.026);    // cannonball solar radiation pressure
auto p = dynamics::DqSgp4Propagator<T>::from_tle(td, tol, mode, dt_max, opt);
```

Additional force models can be injected as lambdas through the explicit constructor,
and adaptive RKF7(8) stepping is available as
`propagate_adaptive(y0, t_target, dt0, tol, dt_min)` on the underlying engine. The
umbrella header `#include "dqsgp4.h"` provides the whole public surface.

For full control, the engine can be assembled by hand:

```cpp
#include "dqsgp4.h"
using T = double;

// 1. Constants (reference ellipsoid). wgs84 / wgs72 / grs80 factories.
constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));

// 2. Body inertia (point_mass / uniform_sphere / diagonal).
dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(math::exact<T>(1));

// 3. The force list. Each entry maps (State, ConstantsProvider) to a Wrench.
//    The geopotential evaluates the monopole and zonal terms in a single
//    Cunningham recursion pass.
constants::GravityField<T> field(constants::ZonalHarmonics<T>::egm2008(T(1e-12)),
                                 constants::TesseralHarmonics<T>{});
std::vector<dynamics::ForceFn<T>> forces;
forces.push_back([field](const dynamics::State<T>& s,
                         const constants::ConstantsProvider<T>& KK) {
    return forces::geopotential(s, KK, field, 4, 0, math::exact<T>(0));
});

// 4. The integrator (euler / runge_kutta_4).
dynamics::IntegratorFn<T> integrator =
    [](const dynamics::State<T>& y0, const math::TrackedValue<T>& dt,
       const integrators::AccelFn<T>& accel) {
        return integrators::runge_kutta_4(y0, dt, accel);
    };

dynamics::Propagator<T> prop(K, inertia, std::move(forces), integrator);

// 5. An initial Cartesian state (metres, m/s), identity attitude, zero spin.
dynamics::State<T> s0 = dynamics::State<T>::from_kinematics(
    math::Quaternion<T>::identity(), r0 /*Vector3<T> m*/, omega0,
    v0 /*Vector3<T> m/s*/, math::exact<T>(0));

// 6. Propagate to t_target with a maximum step of dt_max.
dynamics::State<T> sf = prop.propagate_to(s0, t_target, dt_max);
// sf.position() (m) and sf.linear_velocity() (m/s) are error-tracked.
```

A complete, runnable program covering both propagators, the perturbation options,
adaptive stepping, the error budgets, and precision selection is
[`examples/quickstart.cpp`](examples/quickstart.cpp); it is compiled and executed as
part of the test suite, and its checks include the published SGP4 reference values for
satellite 00005. The generated documentation under [`docs/`](docs/index.html) includes
a task-oriented guide (`docs/guide.html`) and a complete API index (`docs/api.html`).

## Force models and their verification

| Force | Model | Verified against |
|---|---|---|
| Gravity | Monopole, zonal, and tesseral terms evaluated in a single Cunningham recursion pass | Independent per-term implementations and the closed-form J₂ acceleration; agreement at the 5e-16 level |
| Drag | Pluggable density interface; a 27-band piecewise-exponential static atmosphere (Vallado Table 8-4) | The published table values; band-to-band consistency verified at 1e-4 |
| Third body | Newtonian tidal attraction of the Sun and Moon in Battin's cancellation-free form | JPL DE430 ephemeris positions |
| Radiation pressure | Cannonball model with a cylindrical Earth shadow; the solar pressure constant is computed from the IAU nominal solar luminosity | The published 1361 W/m² total solar irradiance (the computed value implies 1361.17) |

## Earth-fixed outputs and ground tracks

The full IAU Earth-orientation chain is implemented in
[`astronomy/earth_rotation.h`](src/astronomy/earth_rotation.h):
`gcrs_to_itrs(t_tt, jd_ut1, x_p, y_p)` composes polar motion, Earth rotation,
the IAU 2000A nutation series (678 terms, truncatable with a tracked error bound),
and IAU 2006 precession. Each layer is verified against the ERFA reference
implementation; polar motion matches exactly and the full matrix agrees with
`c2t06a` to about 0.4 milliarcseconds. The example
[`examples/ground_track.cpp`](examples/ground_track.cpp) computes a ground track from
a TLE through the numerical propagator and this chain.

## Building and testing

See [`BUILDING.md`](BUILDING.md). In short: build `sgp4.sln` (Release | x64) with
MSBuild and run the test executables from `build/Release/`, or run the whole
acceptance suite with `pwsh tools/run_acceptance.ps1`. The headline regression is
`test_sgp4.exe`, which reproduces the official SGP4 verification set (33 satellites,
623 points).

## Layout

- `src/dqsgp4.h` — the umbrella header (the public surface in one include).
- `src/dynamics/`, `src/forces/`, `src/integrators/`, `src/atmosphere/` — the
  dual-quaternion propagator, its force models, and the atmosphere models.
- `src/sgp4/` — the analytical SGP4/SDP4 implementation and its model selector.
- `src/ephemeris/` — Sun and Moon series ephemerides and the GCRS position chain.
- `src/astronomy/` — epochs, reference frames, precession and nutation, sidereal
  angles, and the GCRS-to-ITRS chain.
- `src/math/` — `TrackedValue`, quaternion and dual-quaternion algebra, series
  helpers, and Kepler solvers.
- `src/geodesy/`, `src/perturbation/`, `src/tle/` — the level ellipsoid, the
  analytical perturbation terms, and TLE/OMM parsing.
- `design/` — specifications, audits, and the mathematical derivation notes behind
  each module.
- `tests/` — the per-module test programs that make up the acceptance suite.
- `examples/` — `quickstart.cpp` and `ground_track.cpp`, both built and run by the
  test suite.
- `docs/` — the generated documentation set (regenerate with
  `python tools/gen_docs.py`; a test verifies it stays current).

## References

- Hoots, F. R. & Roehrich, R. L. (1980), *Spacetrack Report #3*.
- Vallado, D. A. et al. (2006), *Revisiting Spacetrack Report #3* (AIAA 2006-6753).

## License

Copyright (c) 2026 Graham Mueller (remstadt@gmail.com).

This project is source-available under the
[PolyForm Noncommercial License 1.0.0](LICENSE.md)
(SPDX: `PolyForm-Noncommercial-1.0.0`). It is free to use, modify, and share for
noncommercial purposes — personal use, research, education, and evaluation. It is
deliberately not an OSI "open source" license; commercial use requires a separate
license. Contact remstadt@gmail.com.

Third-party data provenance and acknowledgements are recorded in
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md). The software comes with no warranty
(see the license) and is not intended for operational, safety-of-life, or
conjunction-assessment use.
