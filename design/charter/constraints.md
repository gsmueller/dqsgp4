# Project Constraints

Hard constraints. Code or designs that violate any of these must be revised or
rejected. Each constraint has a stable identifier (`CON-N`) cited by downstream
specs and tests.


## CON-1 — No magic numbers in computational paths

Every numeric constant used in a derivation, formula, or coefficient is computed
at construction time from inputs supplied by the caller (the `ConstantsProvider<T>`
chain or the call site). Hardcoded floating-point literals such as `1.0826E-3`
(J₂), `6378.137` (a [km]), or `0.7790572732640` (ERA polynomial constant) appear
only inside named-constant factory functions — never in algorithmic code.

Rationale. A propagator that bakes in a 1970s ellipsoid is locked to a 1970s
answer.


## CON-2 — No global state

All state is local to a propagator instance. No `static` mutable variables, no
singletons, no `thread_local` caches, no environment reads in algorithmic code.
Two propagators instantiated in the same process must not interfere.

Rationale. Multi-threaded ensemble propagation and deterministic replay both
require instance isolation.


## CON-3 — No hidden force model

Every force or torque contribution that appears in the integrated wrench must have
been passed in as a lambda at construction time. The propagator must not call a
named force function (e.g., `gravity_J2`, `drag_exponential`) from inside itself.

Rationale. If you can't read off the force model from the constructor arguments,
you can't reason about what you propagated.


## CON-4 — Rigorous upper bounds, not estimates

Every error reported by the library is a rigorous mathematical upper bound on the
true error, derived from the operation's analytical properties — not a statistical
estimate. `total_error()` is a guarantee.

Rationale. Scientific propagation requires "this answer is correct to within X" —
not "this answer is probably correct to within X."


## CON-5 — Three-error tracking is mandatory

Every `TrackedValue<T>` in the computational path carries all three of
(measurement, precision, accuracy). Operations that drop any error category for
performance or simplicity are not permitted.

Rationale. Knowing which error dominates is the whole point of separated budgets.
Discarding one renders the others uninterpretable.


## CON-6 — Truncations and residuals are accounted

Series truncations, Taylor remainder terms, Newton/Halley convergence residuals,
and fixed-iteration cutoffs all add their rigorous bound to the result's
`precision` error budget. Model truncations (e.g., a J₂-only gravity field used to
represent a full spherical-harmonic field) add to the `accuracy` error budget.
Approximations that compute a value without recording the approximation error are
forbidden.

Rationale. An unrecorded approximation is an unbounded error.


## CON-7 — No third-party dependencies beyond vcpkg Boost

The library depends only on the C++20 standard library and Boost (multiprecision,
math constants) via vcpkg. No Eigen, no GSL, no Ceres, no external numerical
solvers, no FFI calls.

Rationale. Portability, reproducibility, and audit-trail integrity.


## CON-8 — Header-only where reasonable

Library code is header-only by default. Compilation-unit dependencies are reserved
for code that is genuinely impractical as templates (e.g., parsers, factories of
large static tables). The existing math library is header-only and the dual
quaternion propagator must remain so.

Rationale. Template specialization on `T` requires header-only at the point of
instantiation; consistency with the existing pattern minimizes friction.


## CON-9 — No application coupling

Code in `src/math/`, `src/dynamics/`, `src/constants/`, `src/forces/`, and
`src/integrators/` must not reference application-specific identifiers (`sgp4`,
`tle`, `brouwer`, `hoots`, etc.). The library does not know its consumers.

Rationale. See OBJ-6. Reversing this coupling later is more expensive than
preventing it now.


## CON-10 — Buildable under MSBuild via `build.bat`

The library and its tests build clean under Visual Studio 2026 / MSBuild, with
`build.bat` as the canonical entry point. CMake or alternate paths are permitted
as additional options but not as substitutes.

Rationale. This is the project's developer-environment baseline.


## CON-11 — Doxygen-ingestible documentation

Every public symbol is documented in a Doxygen-recognized form (`///` for symbols,
`//` for file headers, or `/** */` blocks where structure benefits). The existing
`Doxyfile` must continue to extract a complete set of pages from the source tree
without manual intervention.

Rationale. Discoverability and maintainability scale with documentation, not code.


## CON-12 — All claims tested

Every requirement in `design/specifications/` corresponds to at least one test in
`tests/`. New requirements ship with new tests. The audit criteria in
`design/audit/` are themselves enforced by audit-test programs that exit non-zero
on violation.

Rationale. Untested specifications are aspirations.


## CON-13 — One voice across the library

Code style, comment style, naming, lexicon, and Doxygen conventions are uniform
across `src/math/`, `src/dynamics/`, `src/constants/`, `src/forces/`,
`src/integrators/`. See `design/audit/code_consistency.md` for the binding style
contract.

Rationale. A library written in many voices reads as many libraries. The audit
folder defines one voice.


## CON-14 — No data exfiltration in algorithmic code

Library code must not perform network access, file I/O outside explicit
serialization functions, telemetry emission, or environment reads. Side effects
are confined to the function output and the caller-supplied logger (if any).

Rationale. Determinism (OBJ-8), reproducibility, and operational safety.
