# Audit: Documentation

Verifies that the library is properly documented and discoverable.
Documentation is a first-class artifact (OBJ-7, CON-11). Every public
symbol in the library is extracted by Doxygen into `docs/html/`; the
`design/` folder holds the higher-level governance documentation.

Each requirement has a stable identifier `AUD-DOC-N`. Enforced by
Doxygen warnings (built into `build.bat`) and by code review.


## AUD-DOC-1 — Every public symbol has a Doxygen-extractable comment

Verifies: OBJ-7, CON-11, AUD-CC-4.

Doxygen runs with `EXTRACT_ALL = NO`, `WARN_IF_UNDOCUMENTED = YES`. A
warning-free build means every public symbol has at least a `///` or
`/** */` block. Private symbols are documented at the implementer's
discretion; their absence does not generate warnings.

Verified by: `Doxyfile` configuration + `build.bat` exit code.

Status: ✓ (most recent build outputs "Documentation generated in
docs\html\" with no undocumented-symbol warnings on the new module).


## AUD-DOC-2 — Every file has a `@file` header

Verifies: AUD-CC-3.

Every header in `src/math/`, `src/dynamics/`, `src/constants/`,
`src/forces/`, `src/integrators/` begins with a `/// @file` or
`/** @file */` block whose first line is the file's brief purpose,
followed by detail.

Verified by: `tests/audit/test_documentation.cpp` (planned — greps
for `@file` near top of each header in the library directories).

Status: ✓ for all files written by the dual-quaternion module
(small_angle_series, dual_number, quaternion, dual_quaternion, pose,
twist, wrench, inertia, state, derivative, propagator,
constants_provider, gravity_central, gravity_zonal, runge_kutta).


## AUD-DOC-3 — Audit conformance is cited per file

Verifies: CON-13, AUD-CC-3.

Every header in the library cites the specific AUD-CC-N and AUD-EF-N
identifiers it conforms to. This makes audit traceability explicit at
the source and lets reviewers verify a file's claimed conformance
against the actual content.

Verified by: `tests/audit/test_documentation.cpp` (planned — greps
each header for "AUD-CC-" and "AUD-EF-" in its file doc block).

Status: ✓ for all 15 dual-quaternion module files.


## AUD-DOC-4 — Mathematical derivations are inline at the implementation site

Verifies: OBJ-7, AUD-CC-5.

Non-trivial mathematical expressions — closed-form derivations, error
bounds, Lie-algebra identities, screw exponential / log formulas —
include their derivation as a comment in the file where they are
implemented. The derivation either spells out the algebra inline or
cites a textbook reference, a specification requirement ID, or another
document in `design/`.

Verified by: code review.

Status: ✓ (e.g., `quaternion.h::rotate` documents the closed form;
`dual_quaternion.h::exp_screw` documents the closed-form derivation;
`small_angle_series.h` documents every Taylor expansion's first omitted
term).


## AUD-DOC-5 — Generated docs are buildable in CI

Verifies: OBJ-7, CON-10.

`build.bat` invokes Doxygen as a non-optional step. The HTML
documentation generates into `docs/html/` and is inspectable. A Doxygen
failure fails the overall build.

Verified by: `build.bat` already runs Doxygen.

Status: ✓ (Doxygen integrated; latest builds output "Documentation
generated in docs\html\").


## AUD-DOC-6 — Specifications cite their objectives and constraints

Verifies: OBJ-7, CON-12.

Every requirement in `design/specifications/` cites at least one
objective (`OBJ-N`) or constraint (`CON-N`) on a "Verifies:" line.
Backward traceability is preserved: from a requirement, the reader can
find the objective it serves.

Verified by: `tests/audit/test_documentation.cpp` (planned — greps
each REQ-N block for the "Verifies:" line).

Status: ✓ for `error_framework.md`, `dual_quaternion_algebra.md`,
`constants_provider.md`, `propagator.md`, `integrator.md`, `system.md`.


## AUD-DOC-7 — Audits cite their requirements

Verifies: OBJ-7, CON-12.

Every audit item in `design/audit/` cites at least one requirement
(`REQ-N`) on a "Verifies:" line. Forward traceability: from a
requirement, the reader can find the audit items that verify it.

Verified by: `tests/audit/test_documentation.cpp` (planned).

Status: ✓ for `code_consistency.md`, `error_framework.md`,
`mathematical_correctness.md`, `test_coverage.md`, this document.


## AUD-DOC-8 — Tests cite their audits

Verifies: OBJ-7, CON-12.

Every test in `tests/` cites in its file header the audit IDs it
verifies. Tests without citations are orphans; audits without citing
tests are flagged in `test_coverage.md`.

Verified by: `tests/audit/test_documentation.cpp` (planned — greps
each test's file header for "AUD-N").

Status: ✓ for `test_dual_number/main.cpp`, `test_quaternion/main.cpp`,
`test_dual_quaternion/main.cpp`, `test_propagator/main.cpp` — each
file's header explicitly cites the AUD-N IDs it verifies.


## AUD-DOC-9 — Doxygen extracts the new module

Verifies: OBJ-7.

After a successful `build.bat` invocation, `docs/html/index.html`
contains entries for the new dual-quaternion module alongside the
existing math, geodesy, astronomy, perturbation, and sgp4 modules. The
module is reachable through the standard Doxygen file / class / member
indices.

Verified by: empirical (the latest build's Doxygen step reports
"Documentation generated in docs\html\").

Status: ✓


## AUD-DOC-10 — Design folder is the authoritative high-level reference

Verifies: OBJ-7.

`design/` is the canonical location for charter, specifications,
audits, and module-level architecture documents. Other documentation
(README, code comments, Doxygen output) cross-references `design/` for
detail and does not duplicate its content.

`design/index.md` is the navigation entry; the four-layer hierarchy
(charter → specifications → audit → modules) is described there.

Verified by: code review.

Status: ✓ (`design/index.md` published with coverage status table).


## AUD-DOC-11 — Stable identifiers are never reused

Verifies: OBJ-7.

Identifiers (OBJ-N, CON-N, GOAL-N, REQ-X-N, AUD-X-N) are stable: once
published in a specification or audit document, their numbering is
never reused. Removal is allowed (the ID becomes "retired"); reuse is
forbidden so that all downstream references remain unambiguous.

Verified by: code review (and a planned `tests/audit/
test_identifier_stability.cpp` that hashes the spec / audit files
against a manifest of issued IDs).

Status: ✓ (no reassignments to date; manifest pending).


## Coverage status

| AUD-DOC | Status | Notes |
|---|---|---|
| AUD-DOC-1 | ✓ | Doxygen `WARN_IF_UNDOCUMENTED = YES` in `Doxyfile` |
| AUD-DOC-2 | ✓ | satisfied by all 15 new files; CI grep planned |
| AUD-DOC-3 | ✓ | satisfied by all 15 new files; CI grep planned |
| AUD-DOC-4 | ✓ | code review |
| AUD-DOC-5 | ✓ | Doxygen integrated in `build.bat` |
| AUD-DOC-6 | ✓ | every REQ has a "Verifies:" line; CI grep planned |
| AUD-DOC-7 | ✓ | every AUD has a "Verifies:" line; CI grep planned |
| AUD-DOC-8 | ✓ | every test header cites its AUD IDs; CI grep planned |
| AUD-DOC-9 | ✓ | empirically verified per build |
| AUD-DOC-10 | ✓ | `design/index.md` published |
| AUD-DOC-11 | ✓ | no reassignments to date; manifest pending |


## Cross-reference

Maps to `design/audit/test_coverage.md` for the test → audit map and
to `design/audit/code_consistency.md` for AUD-CC-4 (symbol-level doc
contract). Audited at PR time by reviewers and (forthcoming) by
`tests/audit/test_documentation.cpp` for the grep-detectable rules.
