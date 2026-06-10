# Audit: Test Coverage

Maps every audit identifier (`AUD-CC-N`, `AUD-EF-N`, `AUD-MC-N`,
`AUD-DOC-N`) to the test or tests that verify it, and conversely every
test to the audit identifiers it covers. New requirements / audits /
tests append to this table. Every audit must be cited by at least one
test; every test must cite at least one audit (REQ-SY-12).


## AUD-CC (code consistency)

| Audit | Verified by |
|---|---|
| AUD-CC-1 (namespace) | code review + `tests/audit/test_code_consistency.cpp` (planned: grep) |
| AUD-CC-2 (`#pragma once`) | `tests/audit/test_code_consistency.cpp` (planned: grep) |
| AUD-CC-3 (file-header doc) | code review + `tests/audit/test_code_consistency.cpp` (planned) |
| AUD-CC-4 (symbol-level doc) | code review + AUD-DOC-1 |
| AUD-CC-5 (math notation) | `tests/audit/test_code_consistency.cpp` (planned: grep `$...$`) |
| AUD-CC-6 (naming) | code review |
| AUD-CC-7 (lexicon) | `tests/audit/test_code_consistency.cpp` (planned: grep forbidden) |
| AUD-CC-8 (operator placement) | code review |
| AUD-CC-9 (closed form first) | code review |
| AUD-CC-10 (`exact<T>` / `ratio<T>`) | `tests/audit/test_code_consistency.cpp` (planned: grep `exact_integer` outside tracked_value.h) |
| AUD-CC-11 (`pi<T>` / `two_pi<T>`) | `tests/audit/test_code_consistency.cpp` (planned: grep boost constants outside angles/tracked_value) |
| AUD-CC-12 (section delimiters) | code review |
| AUD-CC-13 (`using std::`) | `tests/audit/test_code_consistency.cpp` (planned: grep `using namespace`) |
| AUD-CC-14 (no app coupling) | `tests/audit/test_code_consistency.cpp` (planned: grep sgp4/tle/brouwer/hoots) |
| AUD-CC-15 (no magic numbers) | code review + `tests/audit/test_code_consistency.cpp` (planned) |
| AUD-CC-16 (line length 100) | `tests/audit/test_code_consistency.cpp` (planned) |
| AUD-CC-17 (include order) | code review |
| AUD-CC-18 (const correctness) | code review |


## AUD-EF (error framework)

| Audit | Verified by | Status |
|---|---|---|
| AUD-EF-1 (tracked returns) | all 4 test programs (composite types instantiate) | ✓ |
| AUD-EF-2 (per-category derivations) | all 4 tests — every assertion bound by `total_error` | ✓ |
| AUD-EF-3 (series tail → precision) | `tests/test_propagator.cpp` (RK4 over 93 steps produces non-zero reported precision) | ✓ |
| AUD-EF-4 (iter residual → precision) | `solve_kepler` exercised by existing SGP4 tests; no iterative algorithm currently on the DQ path | ✓ (existing) |
| AUD-EF-5 (Taylor branch → precision) | `tests/test_quaternion.cpp`, `tests/test_dual_quaternion.cpp` (small-angle samples + screw round-trip) | ✓ |
| AUD-EF-6 (model truncation → accuracy) | (gap: `gravity_J2` documents the residual but does not write to `errors.accuracy`) | follow-up |
| AUD-EF-7 (composites compose) | all 4 tests (Pose, Twist, Wrench, State all instantiate) | ✓ |
| AUD-EF-8 (no silent error reduction) | `tests/audit/test_code_consistency.cpp` (planned: grep `errors.X =`) | planned |
| AUD-EF-9 (catastrophic regions) | `tests/test_propagator.cpp` (Inertia point-mass safe-div) | ✓ |
| AUD-EF-10 (property-based test exists) | this test-coverage document + `tests/test_*.cpp` | ✓ |


## AUD-MC (mathematical correctness)

| Audit | Verified by |
|---|---|
| AUD-MC-1 (ε² = 0) | `tests/test_dual_number.cpp::test_epsilon_squared` |
| AUD-MC-2 (ring axioms) | `tests/test_dual_number.cpp::test_ring_axioms` |
| AUD-MC-3 (forward AD) | `tests/test_dual_number.cpp::test_forward_ad` |
| AUD-MC-4 (Hamilton associativity) | `tests/test_quaternion.cpp::test_associativity` |
| AUD-MC-5 (identity laws) | `tests/test_quaternion.cpp::test_identity` |
| AUD-MC-6 (conjugate involutive) | `tests/test_quaternion.cpp::test_conjugate_involutive` |
| AUD-MC-7 (conjugate of product) | `tests/test_quaternion.cpp::test_conjugate_of_product` |
| AUD-MC-8 (magnitude multiplicative) | `tests/test_quaternion.cpp::test_magnitude_multiplicative` |
| AUD-MC-9 (inverse identity) | `tests/test_quaternion.cpp::test_inverse` |
| AUD-MC-10 (rotation preserves length) | `tests/test_quaternion.cpp::test_rotation_preserves_length` |
| AUD-MC-11 (rotation composition) | `tests/test_quaternion.cpp::test_rotation_composition` |
| AUD-MC-12 (exp/log round-trip) | `tests/test_quaternion.cpp::test_exp_log_roundtrip` |
| AUD-MC-13 (axis-angle round-trip) | `tests/test_quaternion.cpp::test_axis_angle_roundtrip` |
| AUD-MC-14 (normalize idempotent) | `tests/test_quaternion.cpp::test_normalize_idempotent` |
| AUD-MC-15 (DQ associativity) | `tests/test_dual_quaternion.cpp::test_dq_associativity` |
| AUD-MC-16 (conjugate dualities) | `tests/test_dual_quaternion.cpp::test_conjugate_dualities` |
| AUD-MC-17 (pose action composition) | `tests/test_dual_quaternion.cpp::test_pose_action_composition` |
| AUD-MC-18 (screw exp/log round-trip) | `tests/test_dual_quaternion.cpp::test_screw_exp_log_roundtrip` |


## REQ-PR (propagator)

Every propagator requirement is verified by `tests/test_propagator.cpp`.

| Requirement | Phase | Specific assertion |
|---|---|---|
| REQ-PR-1 (pure orchestrator) | construction | Pipeline assembles without library-side glue |
| REQ-PR-2 (force-sum / inertia) | both phases | `acceleration_from_wrench` composes correctly |
| REQ-PR-3 (delegate to integrator) | both phases | `step` calls `runge_kutta_4` |
| REQ-PR-4 (multi-step) | both phases | `propagate_to` lands at exact target time |
| REQ-PR-5 (force signature) | both phases | `(State, K) → Wrench` compiles |
| REQ-PR-6 (once-only construction) | both phases | factories called once at main entry |
| REQ-PR-7 (read-only access) | both phases | `K.earth.J2n(1)` reads through provider |
| REQ-PR-8 (three-error propagation) | both phases | reported precision > 0 |
| REQ-PR-9 (pose retraction) | both phases | `normalized()` inside `lie_advance_pose` |
| REQ-PR-10 (deterministic) | both phases | seeded-free integration produces same result |
| REQ-PR-11 (no force special-cased) | Phase 2 | `gravity_J2` appended; 123 km position effect observed |
| REQ-PR-12 (composable integrator) | both phases | lambda integrator slot |


## REQ-IN (integrator)

| Requirement | Verified by |
|---|---|
| REQ-IN-1 (signature) | `tests/test_propagator.cpp` (IntegratorFn instantiation) |
| REQ-IN-2 (Lie-group pose advance) | `tests/test_propagator.cpp` (position closure 28.8 m) |
| REQ-IN-3 (linear twist advance) | `tests/test_propagator.cpp` (energy / Lz conservation) |
| REQ-IN-4 (order claimed) | (planned: Richardson order-convergence test) |
| REQ-IN-5 (no global state) | `tests/test_propagator.cpp` (deterministic by REQ-PR-10) |
| REQ-IN-6 (adaptive feedback) | (planned: rkf78 test) |
| REQ-IN-7 (retraction every step) | `tests/test_propagator.cpp` (unit-DQ invariant maintained) |
| REQ-IN-8 (three-error propagation) | `tests/test_propagator.cpp` (precision sum > 0) |
| REQ-IN-9 (symplectic) | (planned: symplectic_leapfrog test) |
| REQ-IN-10 (composable forces) | `tests/test_propagator.cpp` Phase 1 + Phase 2 |
| REQ-IN-11 (Lie-group on unit-DQ) | `tests/test_propagator.cpp` (via DualQuaternion::exp_screw) |
| REQ-IN-12 (complexity declared) | code review of `runge_kutta.h` |


## REQ-CP (constants provider)

| Requirement | Verified by |
|---|---|
| REQ-CP-1 (single point of injection) | `tests/test_propagator.cpp` (everything reads through K) |
| REQ-CP-2 (named factory methods) | `tests/test_propagator.cpp` (uses `wgs84`) |
| REQ-CP-3 (three-error tracking) | `tests/test_propagator.cpp` (K.earth.GM has populated errors) |
| REQ-CP-4 (earth sub-bundle) | `tests/test_propagator.cpp` (reads `K.earth.GM`, `K.earth.a`, `K.earth.J2n(1)`) |
| REQ-CP-5..8 (forthcoming sub-bundles) | (planned: drag/third-body/SRP tests) |
| REQ-CP-9 (once-only construction) | `tests/test_propagator.cpp` (single `wgs84(tol)` call) |
| REQ-CP-10 (no global state) | (verified by absence of static state in `constants/`) |


## REQ-EF and REQ-DQ

Cross-referenced in their respective specifications. Each requirement
in `error_framework.md` and `dual_quaternion_algebra.md` cites an audit
ID that verifies it, and the audit table above maps that audit to a
test.


## Test-program scoreboard

| Test program | Passes (current) | Fails |
|---|---|---|
| `test_dual_number` | 113 / 113 | 0 |
| `test_quaternion` | 104 / 104 | 0 |
| `test_dual_quaternion` | 72 / 72 | 0 |
| `test_propagator` | 8 / 8 (both phases) | 0 |
| **Total** | **297 / 297** | **0** |


## Acknowledged gaps

These are documented limitations of the current coverage, not bugs.
Each has a planned remediation.

- **AUD-CC-\* automated checks**: `tests/audit/test_code_consistency.cpp`
  not yet implemented. Currently relies on code review and grep.
- **AUD-EF-6**: `gravity_J2` documents its model-truncation residual
  but does not write to `errors.accuracy`. A future pass will record
  the residual through the framework.
- **REQ-CP-5..8 sub-bundles**: gravity_field, astronomy, time, and
  fundamental_constants types are specified but not yet implemented.
  Land alongside the force lambdas that consume them (drag, third-body,
  SRP).
- **REQ-IN-4 order-convergence**: a Richardson test halving dt and
  asserting the expected ratio of position errors is planned but not
  written.
- **REQ-IN-6 adaptive feedback**: requires `rkf78` (planned).
- **REQ-IN-9 symplectic**: requires `symplectic_leapfrog` (planned).
- **REQ-SY-7 high-precision**: `test_propagator_high_precision` at
  `T = cpp_bin_float_50` is planned.
- **REQ-SY-9 constants swap**: `test_constants_provider` propagating
  the same orbit under wgs84/wgs72/grs80 is planned.
