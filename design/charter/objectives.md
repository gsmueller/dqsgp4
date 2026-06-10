# Project Objectives

This file is binding. Every implementation choice must serve one or more of these
objectives. Every requirement in `design/specifications/` and every audit criterion
in `design/audit/` must trace back to an objective here. A feature that does not
advance an objective does not belong in the codebase.

Each objective has a stable identifier (`OBJ-N`) cited by downstream specs and tests.


## OBJ-1 — Dual quaternion six-DOF state propagator

Provide a propagator whose state is the pair (M̂, Ω̂) of a unit pose dual quaternion
and a pure twist dual quaternion. The propagator advances this state in time under
an injected wrench model that may include gravity, drag, third-body, radiation
pressure, control inputs, and arbitrary user-supplied terms.

Rationale. Dual quaternions encode SE(3) (position + attitude) with constant-time
composition, exact unit-norm enforcement under Lie-group retraction, and
singularity-free interpolation. They are the right representation for a six-DOF
propagator that may need to handle attitude-coupled drag, propulsive maneuvers,
gravity-gradient torque, and rigid orbital dynamics simultaneously, without an
architectural change between cases.


## OBJ-2 — Extensible accuracy and precision

The numeric precision of every quantity can be raised — to arbitrary precision via
Boost multiprecision — without recompiling client code or modifying library code.
State, force model, integrator, and constants are template-parameterized on the
underlying scalar type T; precision is selected at the call site.

Rationale. Applications range from quick LEO survey runs (`double` suffices) to
scientific propagation where 50+ decimal digits are required to compare analytical
claims against numerical baselines. The propagator must serve both without forks.


## OBJ-3 — Three-error tracking on every output

Every quantity computed inside the propagator is a `TrackedValue<T>` carrying
rigorous upper bounds on three independent error sources:

- measurement — propagated from input physical uncertainty
- precision   — accumulated numerical roundoff, truncation, residual
- accuracy    — accumulated model truncation and simplification

Rationale. A propagation result without error bounds is unverifiable. Tracking
three independent budgets lets the caller see which source dominates and where
investment will reduce total error.


## OBJ-4 — Lambda-injected force, torque, and ephemeris models

Every contribution to the body wrench, every ephemeris callback, every time-system
conversion, every Earth-orientation transformation is supplied by the caller as a
`std::function`. The propagator contains no force model, no ephemeris, and no
sidereal-time formula.

Rationale. The propagator is a generic engine. Force-model assumptions and choice
of ephemeris are research questions, not architectural ones. Swapping a J₂ zonal
for a full spherical-harmonic model must be a single-line change at the call site.


## OBJ-5 — Swappable fundamental constants

A `ConstantsProvider<T>` bundle injects the Earth ellipsoid, gravity field,
astronomy, time system, and physical constants. Factories produce WGS84, EGM2008,
JGM-3, or user-defined providers. The propagator reads constants only through this
interface; there is no global default.

Rationale. Model variability includes which constants you believe. The propagator
must not commit to a single set.


## OBJ-6 — Reusable, application-independent

The math library is general-purpose. Nothing in `src/math/`, `src/dynamics/`,
`src/constants/`, `src/forces/`, or `src/integrators/` may reference SGP4, TLE, or
any other application-specific concept. Applications consume the library; the
library does not name them.

Rationale. A generic six-DOF propagator is valuable across many missions. Coupling
to one application locks the value behind one model.


## OBJ-7 — Documentation as a first-class artifact

Every public symbol is documented in Doxygen-extractable form. Every algorithmic
file has a derivation in `design/` describing the formulas it implements, with
citations. Every requirement is testable and every test cites its requirement.

Rationale. A propagator whose error bounds, force-model assumptions, or constants
the user must reverse-engineer from source is not usable for science.


## OBJ-8 — Deterministic, reproducible, sandbox-friendly

Identical inputs produce identical outputs across compiler versions, optimization
levels, and runs. There is no randomness, no environment lookup, no timing
dependency, no platform-specific behavior in the computational path.

Rationale. Bit-for-bit reproducibility is the contract under which "rigorous
upper bound" remains meaningful when the result is recomputed.


## Traceability

| Objective | Primary realization |
|---|---|
| OBJ-1 | `dynamics/propagator.h`, `dynamics/state.h`, `math/dual_quaternion.h` |
| OBJ-2 | template T on every public type; `boost::multiprecision::cpp_bin_float_*` |
| OBJ-3 | `math/tracked_value.h`, propagated through every operator |
| OBJ-4 | `dynamics/wrench.h` + `std::function` injection at construction |
| OBJ-5 | `constants/constants_provider.h` + factories |
| OBJ-6 | absence of application identifiers in library namespaces |
| OBJ-7 | Doxyfile coverage check + `design/` content audit |
| OBJ-8 | absence of `static` mutable state, RNG, environment reads |
