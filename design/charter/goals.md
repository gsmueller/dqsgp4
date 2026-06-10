# Project Goals

Goals are measurable targets that the project commits to reaching. Each goal is
referenced by at least one requirement in `design/specifications/` and verified by
at least one test in `tests/`.

Goals differ from objectives (which describe the project's purpose) and
constraints (which describe boundaries that must hold at all times). A goal is
something you aim at and either hit or miss; a constraint is something you do not
violate.

Each goal has a stable identifier (`GOAL-N`).


## GOAL-1 — Reference LEO orbit converges to within 1 m over one period at `T = double`

Propagate a circular 400 km LEO orbit for one orbital period (≈ 5560 s) under
point-mass gravity (no oblateness) with `T = double` and RK4 at 60 s step size.
Position error — defined as the Euclidean distance between the propagated state
and the analytical closed-form Keplerian solution at the same time — must be below
1 m. The reported precision-error component on the position must be a strict upper
bound: `total_error() ≥ |actual error|` at every sample point.

Verifies: OBJ-1, OBJ-3, CON-4.


## GOAL-2 — Same orbit, `T = cpp_bin_float_50`, error < 1 mm

Repeat GOAL-1 with `T = boost::multiprecision::cpp_bin_float_50` and a tighter
integrator (RKF7(8), reltol = 1e-30). Position error must be below 1 mm.
`total_error()` must remain a strict upper bound.

Verifies: OBJ-2, OBJ-3.


## GOAL-3 — Constants-provider swap changes results consistently

Propagate the same Keplerian initial state under WGS84 and EGM2008 providers. The
first-period position difference must equal the difference predicted by the GM and
$a^2 \omega$ implied by the two providers, to within the combined precision of
both runs.

Verifies: OBJ-5, OBJ-4.


## GOAL-4 — J₂ injection produces precessing orbit

Add `gravity_J2(K)` to the wrench list. The argument of perigee precesses at the
analytical rate $\dot\omega = (3 n J_2 / 2)(R_E / a)^2 (4 - 5\sin^2 i)/(1 - e^2)^2$
to within the integrator's reported precision.

Verifies: OBJ-1, OBJ-4.


## GOAL-5 — Attitude propagation conserves angular momentum

Propagate a torque-free spinning rigid body. Angular momentum in the world frame
is conserved up to the integrator's order; precision error reported on the
momentum vector is a strict upper bound on the drift.

Verifies: OBJ-1.


## GOAL-6 — Library builds clean under MSBuild

`build.bat` returns success on a clean checkout. No new warnings beyond the
existing baseline at `/W4`.

Verifies: CON-10.


## GOAL-7 — 100 % Doxygen coverage of public API

Doxygen runs with `EXTRACT_ALL = NO`, `WARN_IF_UNDOCUMENTED = YES`,
`WARN_AS_ERROR = YES` and exits zero. Every public symbol in `src/` is documented.

Verifies: OBJ-7, CON-11.


## GOAL-8 — Every specification has at least one test

A cross-reference audit (script in `tests/audit/`) confirms every requirement ID
in `design/specifications/*.md` is cited by the name of at least one test in
`tests/`, and every test cites at least one requirement ID. No orphans in either
direction.

Verifies: CON-12.


## GOAL-9 — Smoke test in CI

Running `tests/test_smoke.exe` from a clean build returns exit code 0 on the
reference LEO orbit. Energy is conserved to 1 part in 10^7; position error within
1 m at end of one orbital period.

Verifies: OBJ-1, OBJ-3.


## GOAL-10 — Self-audit of error framework

A dedicated test program (`tests/audit/test_error_framework.cpp`) generates random
inputs at known precision, runs them through every operator on `TrackedValue<T>`,
`Quaternion<T>`, `DualQuaternion<T>`, `Vector3<T>`, and compares the reported
`total_error()` to the actual numerical error measured against a high-precision
reference. The reported bound must dominate the actual error in 100 % of trials.

Verifies: OBJ-3, CON-4, CON-5, CON-6.


## GOAL-11 — Style consistency across `src/`

A static audit (script in `tests/audit/`) confirms every file in `src/math/`,
`src/dynamics/`, `src/constants/`, `src/forces/`, `src/integrators/` satisfies the
contract in `design/audit/code_consistency.md`: naming, comment style, math
notation, lexicon, helper-method conventions, error-framework integration.

Verifies: CON-13.


## GOAL-12 — Long-haul drift bounded

Propagate a GPS-altitude orbit (a ≈ 26 600 km) for 30 days under WGS84 + J₂.
Reported `total_error()` on the final position is below 10 m. Compared against an
EGM2008-based reference run, position difference at 30 days is within the
reported bound (i.e., the bound is tight enough to capture the model-truncation
gap from J₂-only).

Verifies: OBJ-3, CON-4, CON-6.
