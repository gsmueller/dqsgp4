# Specification: System

Top-level system requirements that span multiple subsystems. Each
requirement has a stable identifier `REQ-SY-N` and either decomposes
into more specific requirements in a sub-spec or is verified directly
by an end-to-end test.

Sub-specifications:
- `error_framework.md` — REQ-EF-N (rigorous three-error tracking)
- `dual_quaternion_algebra.md` — REQ-DQ-N (algebraic laws)
- `constants_provider.md` — REQ-CP-N (constants injection)
- `propagator.md` — REQ-PR-N (propagator orchestration)
- `integrator.md` — REQ-IN-N (time-stepping)

Top-level requirements here either constrain the system as a whole or
realize an objective without a dedicated sub-spec.


## REQ-SY-1 — Inputs are explicit

Verifies: OBJ-5, CON-1, CON-2.

Every numeric input to a propagator instance arrives through the
constructor: `ConstantsProvider`, `Inertia`, force lambdas,
`IntegratorFn`, initial `State`. There are no defaults, no globals, no
environment reads in algorithmic code.

A reader can determine a propagator's complete dependency by inspecting
its constructor arguments.


## REQ-SY-2 — Outputs are explicit and bounded

Verifies: OBJ-3, CON-4, CON-5.

The propagator's outputs are the `State<T>` at requested times. Every
component of the output state carries its full three-error budget
(measurement, precision, accuracy), and `total_error()` is a rigorous
upper bound on the deviation from the value that would be obtained by
infinite-precision arithmetic under the same model.

A reader can convert any output value into a confidence statement
without consulting external information.


## REQ-SY-3 — End-to-end test exists and passes

Verifies: every OBJ.

`tests/test_propagator.cpp` exercises the full pipeline:
`ConstantsProvider → Inertia → forces → integrator → Propagator →
propagate_to`. It verifies position closure, energy conservation,
angular momentum conservation, three-error population, and
multi-lambda force composition on a 400 km LEO orbit (WGS 84 central
gravity ± J₂).

Status: build is green; every assertion passes.


## REQ-SY-4 — Library is reusable beyond SGP4

Verifies: OBJ-6, CON-9.

Code in `src/math/`, `src/dynamics/`, `src/constants/`, `src/forces/`,
`src/integrators/` does not reference SGP4, TLE, Brouwer, Hoots, or any
other application identifier. A grep over those directories for those
terms returns zero matches.

Audited by `design/audit/code_consistency.md` AUD-CC-14.


## REQ-SY-5 — Build reproducibility

Verifies: OBJ-8, CON-10.

A clean checkout on a system with Visual Studio 2026 + vcpkg-installed
Boost produces an identical set of binaries via `build.bat`, regardless
of run order or environment beyond the standard MSBuild dependencies.
No network access during build (vcpkg restore is a prerequisite); no
generated-code step beyond Doxygen documentation.


## REQ-SY-6 — Documentation tooling is integrated

Verifies: OBJ-7, CON-11.

`build.bat` invokes Doxygen as part of the standard build. Every public
symbol in `src/` is extracted into the generated HTML under
`docs/html/`. New modules added by following the existing file-header
and symbol-doc conventions (AUD-CC-3, AUD-CC-4) are picked up
automatically.

Audited by `design/audit/documentation.md` AUD-DOC-N.


## REQ-SY-7 — Numeric precision is selectable per call site

Verifies: OBJ-2.

Every public type and function in the library is template on `T`. The
caller selects `T` at the call site: `double`,
`boost::multiprecision::cpp_bin_float_50`, `cpp_bin_float_100`, or any
compatible scalar type. No source change is required to raise
precision; the same code paths run.

`test_propagator.cpp` instantiates with `T = double`. A forthcoming
`test_propagator_high_precision.cpp` will demonstrate the same code at
`T = cpp_bin_float_50` with tightened tolerances reflecting the wider
representation.


## REQ-SY-8 — Force-model swap is a single edit at the call site

Verifies: OBJ-4, REQ-PR-11.

Adding J₂ to a central-gravity propagator is one line:

    forces.push_back([](const State<T>& s,
                        const ConstantsProvider<T>& K) {
        return forces::gravity_J2(s, K);
    });

Removing it is one line. No library code changes. The propagator does
not know which forces are in its list. `tests/test_propagator.cpp`
Phase 2 validates this directly.


## REQ-SY-9 — Constants swap is a single factory selection

Verifies: OBJ-5, REQ-CP-2.

Switching between WGS 84, WGS 72, and GRS 80 is a single function call:

    constants::ConstantsProvider<T>::wgs84(tol)   // vs
    constants::ConstantsProvider<T>::wgs72(tol)   // vs
    constants::ConstantsProvider<T>::grs80(tol)

The propagator does not know which convention is in use. A forthcoming
`tests/test_constants_provider.cpp` will propagate the same initial
state under all three providers and verify the differences match the
analytical prediction from the differing `GM` / `R_E` / `J₂` / `ω`.


## REQ-SY-10 — Integrator swap is a lambda assignment

Verifies: OBJ-4, REQ-PR-12, REQ-IN-1.

Switching from RK4 to Euler (or to any future RKF7(8), DOPRI8,
symplectic scheme) is one line: assign a different `IntegratorFn<T>` to
the propagator's integrator slot. The propagator does not know which
scheme is in use; force lambdas do not know either.


## REQ-SY-11 — No application-coupled state or inheritance

Verifies: OBJ-6, CON-9.

The library types (`Pose`, `Twist`, `Wrench`, `Inertia`, `State`,
`Derivative`, `DualQuaternion`, `Quaternion`, `DualNumber`, `Vector3`,
`TrackedValue`) are self-contained value types. None inherits from an
application class; none has a virtual method bound to an application
interface; none takes an application object as a constructor argument.

This permits the library to be vendored into any application without
import cycles.


## REQ-SY-12 — Test-coverage map is explicit and complete

Verifies: CON-12, OBJ-7.

Every requirement identifier (REQ-EF-N, REQ-DQ-N, REQ-CP-N, REQ-PR-N,
REQ-IN-N, REQ-SY-N) is cited by at least one audit identifier (AUD-N)
in `design/audit/`, and every audit identifier is cited by at least
one test in `tests/`.

The full mapping is published in `design/audit/test_coverage.md`.
Adding a new requirement obligates adding an audit and a test; adding
a new test obligates citing the audit it verifies.


## REQ-SY-13 — All Doxygen output is non-warning-free

Verifies: AUD-DOC-1.

Doxygen runs with `EXTRACT_ALL = NO`, `WARN_IF_UNDOCUMENTED = YES`. A
clean build produces zero "undocumented" warnings for any public symbol
in the library. Failure indicates either a new symbol missing its
`///` block or a mis-routed Doxyfile rule.


## Cross-reference

Audited holistically by:
- `design/audit/code_consistency.md` (style)
- `design/audit/error_framework.md` (framework conformance)
- `design/audit/mathematical_correctness.md` (algebraic correctness)
- `design/audit/documentation.md` (Doxygen coverage)
- `design/audit/test_coverage.md` (test–requirement mapping)

Verified end-to-end by:
- `tests/test_dual_number.cpp` (REQ-DQ-1..3)
- `tests/test_quaternion.cpp` (REQ-DQ-4..14)
- `tests/test_dual_quaternion.cpp` (REQ-DQ-15..18)
- `tests/test_propagator.cpp` (REQ-PR + REQ-IN + REQ-EF + REQ-SY-3)
