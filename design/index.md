# Design Documentation Index

Project governance and technical-design content for the dual quaternion
propagator. The four layers below trace one direction: charter drives
specifications, specifications drive audits, audits drive tests. The IDs in
each layer are stable: a test cites an `AUD-N`, which cites a `REQ-N`, which
cites an `OBJ-N` or `CON-N`. Renumbering after first publication is forbidden.

> The coverage table at the bottom is kept honest by `tools/check_index.ps1`
> (acceptance gate **W9**): every spine-layer count must equal the actual file
> count on disk, and every specification and audit file must be linked here.


## 1. Charter — why the project exists

Stable IDs: `OBJ-N`, `CON-N`, `GOAL-N`.

- [charter/objectives.md](charter/objectives.md) — OBJ-1..8: what the project
  is for.
- [charter/constraints.md](charter/constraints.md) — CON-1..14: hard
  boundaries.
- [charter/goals.md](charter/goals.md) — GOAL-1..12: measurable targets.


## 2. Specifications — what the system must do

Stable IDs: `REQ-<area>-N` (e.g., `REQ-EF-1`, `REQ-DQ-12`). Each file opens
with a `# Specification: <area>` heading.

- [specifications/error_framework.md](specifications/error_framework.md) —
  REQ-EF-1..15: how errors propagate; the contract for OBJ-3 and CON-4/5/6.
- [specifications/dual_quaternion_algebra.md](specifications/dual_quaternion_algebra.md)
  — REQ-DQ-1..18: algebraic laws for DualNumber, Quaternion, DualQuaternion.
- [specifications/system.md](specifications/system.md) — REQ-SY-1..13:
  system-level inputs/outputs, reusability, build, precision selection.
- [specifications/constants_provider.md](specifications/constants_provider.md)
  — REQ-CP-1..10: single-injection constants bundle and named conventions.
- [specifications/propagator.md](specifications/propagator.md) — REQ-PR-1..12:
  state evolution, force-lambda composition, integrator injection.
- [specifications/integrator.md](specifications/integrator.md) — REQ-IN-1..12:
  the integrator interface and the concrete schemes.


## 3. Audit — how to verify

Stable IDs: `AUD-<area>-N` (e.g., `AUD-CC-7`, `AUD-EF-5`, `AUD-MC-12`). Each
governance audit file opens with a `# Audit: <area>` heading — the marker that
distinguishes the five audit-layer documents from the SGP4-derivation working
notes that also live under `audit/` (e.g. `AUDIT_BACKLOG.md`,
`theoretical_basis_audit/`, `remediation_plan.md`), which are not part of the
charter→spec→audit→test spine.

- [audit/code_consistency.md](audit/code_consistency.md) — AUD-CC-1..18:
  binding style and lexicon contract.
- [audit/error_framework.md](audit/error_framework.md) — AUD-EF-1..10:
  procedural checklist for error-framework conformance.
- [audit/mathematical_correctness.md](audit/mathematical_correctness.md) —
  AUD-MC-1..18: algebraic-law verification, sampling discipline.
- [audit/documentation.md](audit/documentation.md) — AUD-DOC-1..11:
  documentation and cross-citation conformance checklist.
- [audit/test_coverage.md](audit/test_coverage.md) — the audit→test coverage
  matrix: maps every `AUD-CC/EF/MC/DOC` ID to the test(s) that verify it and
  back (REQ-SY-12). Defines no new IDs of its own.


## 4. Module designs — technical-design docs

Architectural narratives for each module. This is an open-ended narrative
layer (not a closed, stable-ID governance set), so it carries no fixed file
count.

- [dual_quaternion_propagator.md](dual_quaternion_propagator.md) — top-level
  architecture: state, time-evolution, lambda-injection, integrator,
  constants provider, error propagation, module map.
- [library_architecture.md](library_architecture.md) — module dependency
  graph (existing math / geodesy / perturbation / sgp4 chain).
- [lambda_injection.md](lambda_injection.md) — `std::function` injection
  pattern (existing).
- [error_bounded_computation.md](error_bounded_computation.md) — series
  error bounds, `BoundedResult<T>`, iterate-to-tolerance (existing).
- Other existing technical-design markdown files in `design/`.


## 5. Test plan — how the audit is realized

Tests cite the audit IDs they verify. The five governance tests below are the
test programs whose source cites at least one `AUD-<area>-N` ID:

- [tests/test_dual_number/main.cpp](../tests/test_dual_number/main.cpp) —
  AUD-MC-1..3 (REQ-DQ-1..3).
- [tests/test_quaternion/main.cpp](../tests/test_quaternion/main.cpp) —
  AUD-MC-4..14 (REQ-DQ-4..14).
- [tests/test_dual_quaternion/main.cpp](../tests/test_dual_quaternion/main.cpp)
  — AUD-MC-15..18 (REQ-DQ-15..18).
- [tests/test_propagator/main.cpp](../tests/test_propagator/main.cpp) —
  AUD-EF-2/3/6 (REQ-EF-3/7): the propagator pipeline exercises tracked-error
  propagation end to end.
- [tests/test_force_models/main.cpp](../tests/test_force_models/main.cpp) —
  AUD-EF-6 (REQ-EF-6): each force's model-truncation residual reaches
  `errors.accuracy` (zero for the exact monopole, nonzero for J₂ and drag).
- [tests/test_code_consistency/main.cpp](../tests/test_code_consistency/main.cpp)
  — AUD-CC-1/2/5/10/11/13/14/16: scans the five library directories and asserts
  the grep-detectable code-consistency items.
- [tests/test_error_framework/main.cpp](../tests/test_error_framework/main.cpp)
  — AUD-EF-1/2/3/4/5/7/9/10 (REQ-EF-2): property tests asserting
  `total_error()` bounds the true error against a cpp_bin_float_50 reference.

All governance tests named in this plan are now implemented.


## Traceability conventions

- Every ID is stable. Renumbering after first publication is forbidden.
- Every requirement cites at least one objective (OBJ-) or constraint (CON-).
- Every audit cites at least one requirement (REQ-).
- Every test name cites at least one audit ID
  (e.g., `test_AUD_MC_4_hamilton_associativity`).
- New IDs append to the existing sequence in each file. Deletion is allowed
  (the ID becomes "retired"); reuse of a retired ID is not.


## Coverage status

Counts are the number of files present on disk for each spine layer (verified
by `tools/check_index.ps1`); a layer is complete when every planned document
exists. Module designs is the open narrative layer and is intentionally not
counted.

| Layer | Files | Coverage |
|---|---|---|
| Charter | 3 / 3 | objectives, constraints, goals |
| Specifications | 6 / 6 | error_framework, dual_quaternion_algebra, system, constants_provider, propagator, integrator |
| Audit | 5 / 5 | code_consistency, error_framework, mathematical_correctness, documentation, test_coverage |
| Module designs | narrative | open-ended; dual_quaternion_propagator + existing technical docs |
| Tests | 7 / 7 | dual_number, quaternion, dual_quaternion, propagator, force_models, code_consistency, error_framework |
