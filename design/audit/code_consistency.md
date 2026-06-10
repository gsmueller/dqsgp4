# Audit: Code Consistency

Binding style and lexicon contract for `src/math/`, `src/dynamics/`,
`src/constants/`, `src/forces/`, `src/integrators/`. Each item has a stable
identifier (`AUD-CC-N`) and cites the originating constraint. A file passes
the audit if every applicable item is satisfied. Enforced by manual review at
PR time and by `tests/audit/test_code_consistency.cpp` (filename / grep based)
at CI time.


## AUD-CC-1 — Namespace

Verifies: CON-9, CON-13.

| Directory | Namespace |
|---|---|
| `src/math/`        | `math` |
| `src/dynamics/`    | `dynamics` |
| `src/constants/`   | `constants` |
| `src/forces/`      | `forces` |
| `src/integrators/` | `integrators` |

No `using namespace` at file scope. No anonymous namespaces in headers.


## AUD-CC-2 — Header guards

Verifies: CON-13.

Every header begins with `#pragma once`. No `#ifndef` guards.


## AUD-CC-3 — File-header documentation

Verifies: CON-11, OBJ-7.

The first lines after `#pragma once` are a doc block describing the file's
purpose. Two permitted forms:

    /// @file <name>.h
    /// <one-line brief>
    /// <optional detail lines>

or

    /** @file <name>.h
     *  <one-line brief>
     *  <optional detail lines>
     */

The block states: what the file provides, what it depends on (the include list
itself is acceptable — narrative dependencies live in `design/`), any non-obvious
context, and the audit-conformance citations (`AUD-CC-N`, `AUD-EF-N`).

Lines do not exceed 100 characters.


## AUD-CC-4 — Symbol-level documentation

Verifies: CON-11, OBJ-7.

Every public symbol (struct, free function, member function, named constructor)
is preceded by a `///` block. First line is the brief; subsequent lines are
detail. `@brief` tag is not used — the first line is the brief by convention.

`@param`, `@return`, `@tparam` are used where they add clarity; not used when the
parameter name carries the meaning.


## AUD-CC-5 — Math notation in code comments

Verifies: CON-13.

ASCII operators (`*`, `+`, `-`, `/`, `^`) with Greek letters and Unicode math
characters where natural:

| Concept | Use |
|---|---|
| Greek letters | `alpha`, `theta` (spelled out) or `α`, `θ` (Unicode) |
| Comparisons | `≤`, `≥`, `≠`, `≈`, `±` |
| Operators | `×`, `·`, `√`, `∞`, `→`, `↦` |
| Exponents | `²`, `³`, `ⁿ` or `^2`, `^3`, `^n` |
| Sums / products | `sum_{k=0}^{N}`, `prod_{k=1}^{n}` |

LaTeX `$...$` math is **not** used in code comments. It is reserved for
`design/` markdown files where Doxygen and GitHub render it.

Equation cross-references cite IDs in `design/` (theorem numbers, proposition
IDs, requirement IDs) — never bare equation numbers that may drift.


## AUD-CC-6 — Naming

Verifies: CON-13.

| Element | Style | Example |
|---|---|---|
| Types | `PascalCase` | `TrackedValue`, `DualQuaternion`, `Pose` |
| Free functions | `snake_case` | `from_axis_angle`, `wrap_two_pi` |
| Member functions | `snake_case` | `magnitude`, `log_unit`, `conjugate` |
| Variables | `snake_case` | `half_angle`, `qv_norm`, `theta_sq` |
| Template parameters | single uppercase or `PascalCase` | `T`, `Func` |
| Named constructors | `snake_case` static methods | `identity`, `zero`, `pure` |
| Private helpers | `snake_case`, **no trailing underscore** | `taylor_sinc` (not `sinc_`) |

Constructor parameters that would shadow public single-letter fields take a
trailing underscore (`x_`, `y_`, `z_`, `w_`). Multi-letter fields use short
ctor-param abbreviations to avoid the shadow (`val` for `value`, `err` for
`errors`, `r` for `real`, `d` for `dual`).


## AUD-CC-7 — Lexicon

Verifies: CON-13.

| Concept | Use | Don't use |
|---|---|---|
| Euclidean length of a vector / quaternion / DQ | `magnitude` | `norm`, `length`, `modulus`, `abs` |
| Squared Euclidean length | `magnitude_squared` | `norm_squared`, `length2` |
| Multiplicative identity | `identity` | `one`, `unit`, `eye` |
| Additive identity | `zero` | |
| Multiplicative inverse | `inverse` | `reciprocal`, `recip` |
| Conjugate (quaternion / DQ) | `conjugate` | `conj`, `bar` |
| Renormalize to unit magnitude (const) | `normalized` | `unit`, `unitize` |
| Apply rotation to a 3-vector | `rotate(v)` | `apply`, `transform`, `act_on` |
| Compose two poses / twists / wrenches | `operator*` | `compose`, `combine` |
| Log of a unit quaternion | `log_unit` | `log`, `axis_angle` |
| Exp of a pure quaternion | `exp_pure` | `from_rotvec`, `expq` |
| Time advance | `step` | `tick`, `advance`, `update` |
| Pose dual quaternion | `Pose` | `Frame`, `Transform`, `SE3Element` |
| Twist dual quaternion | `Twist` | `Velocity`, `VelocityState` |
| Wrench dual quaternion | `Wrench` | `Force`, `ForceTorque`, `GeneralizedForce` |
| Body inertia (dual) | `Inertia` | `MomentOfInertia`, `Mass` |


## AUD-CC-8 — Operator placement

Verifies: CON-13.

Binary operators (`+`, `-`, `*`, `/`, comparisons) on a type T are declared as
hidden friends inside `struct T`:

    friend T operator+(const T& a, const T& b) { ... }

Unary operators (`-`, `!`) and instance queries (`q.conjugate()`,
`v.magnitude()`) are members. Free-function math operations on T (`sin`, `cos`,
`sqrt`, etc.) are friends so they participate in ADL.


## AUD-CC-9 — Closed-form first, series second

Verifies: REQ-EF-3, REQ-EF-4.

Algebraic operations are written in closed form. Series expansions are reserved
for cases where no closed form exists, or as small-argument stability branches
where the closed form is numerically ill-conditioned. Every series operation
conforms to REQ-EF-4. Every fixed-order Taylor branch conforms to REQ-EF-6.


## AUD-CC-10 — `exact<T>(n)` and `ratio<T>(num, den)`

Verifies: CON-13.

Integer tracked constants use `exact<T>(n)`. Rational tracked constants use
`ratio<T>(num, den)`. The long form `TrackedValue<T>::exact_integer(n)` is used
only inside `src/math/tracked_value.h` (its definition site).


## AUD-CC-11 — `pi<T>()` and `two_pi<T>()`

Verifies: CON-13.

Tracked π and 2π come from `src/math/angles.h`. Code that needs π includes
`angles.h` and calls `pi<T>()` / `two_pi<T>()`. Direct access to
`boost::math::constants::pi<T>()` is reserved for `angles.h` and
`tracked_value.h` themselves.


## AUD-CC-12 — Section delimiters

Verifies: CON-13.

Structs with more than ~10 members or methods use `// --- Section ---`
delimiters between logical groups, in this order:

    // --- Constructors ---
    // --- Named constructors ---
    // --- Component access ---
    // --- Arithmetic ---
    // --- Norm and conjugation ---       (algebras)
    // --- Lie algebra: exp and log ---   (group elements)
    // --- Action ---                     (group elements)
    // ... domain-specific groups ...
    // (then private:  and helpers)

Smaller structs may omit delimiters.


## AUD-CC-13 — `using std::` declarations

Verifies: CON-13.

`using std::abs;` etc. appear inside function bodies, never at namespace scope.
ADL on tracked types is preferred when available.


## AUD-CC-14 — No application coupling

Verifies: CON-9.

Library code does not reference SGP4, TLE, Brouwer, Hoots, or any other
application-specific identifier. A grep over the library namespaces for these
terms must return zero hits.


## AUD-CC-15 — No magic numbers

Verifies: CON-1.

Floating-point literals do not appear in algorithmic code. Permitted exceptions:

- Small integers as the argument to `exact<T>(...)`.
- Numerical thresholds for branch selection (e.g., a `1e-4` Taylor-branch
  cutoff). Each threshold is a named constant in the function or struct, with a
  derivation comment that justifies the choice — typically the truncation-bound
  formula at the threshold.
- Named-constant factory functions (`wgs84_provider`, `egm2008_provider`) are
  exempt — that is what they exist for.


## AUD-CC-16 — Line length

Verifies: CON-13.

Lines do not exceed 100 characters. Long expressions break before the operator
and indent by 4 spaces from the line they continue.


## AUD-CC-17 — `#include` order

Verifies: CON-13.

Includes are grouped with blank-line separation:

1. Module headers under `src/` (same library).
2. Boost headers.
3. C++ standard library headers.

Within each group, alphabetical.


## AUD-CC-18 — Const correctness

Verifies: CON-13.

Member functions that do not modify state are `const`. Function parameters that
are not modified are passed by `const&` for non-trivial types (including
`TrackedValue<T>`, `Vector3<T>`, and all composites). Pass scalars (`T`, `int`)
by value.


## Cross-reference

| AUD-CC-N | Grep-detectable | Manual review | Verified by |
|---|---|---|---|
| AUD-CC-1 | yes (namespace decl) | | `test_code_consistency.cpp` |
| AUD-CC-2 | yes (`#pragma once`) | | `test_code_consistency.cpp` |
| AUD-CC-3 | partial (presence) | yes (quality) | reviewer |
| AUD-CC-4 | partial | yes | reviewer |
| AUD-CC-5 | yes (`$...$` regex) | | `test_code_consistency.cpp` |
| AUD-CC-6 | partial | yes (semantic) | reviewer |
| AUD-CC-7 | yes (forbidden-name grep) | | `test_code_consistency.cpp` |
| AUD-CC-8 | partial | yes | reviewer |
| AUD-CC-9 | | yes | reviewer |
| AUD-CC-10 | yes (`exact_integer` outside `tracked_value.h`) | | `test_code_consistency.cpp` |
| AUD-CC-11 | yes (boost-constants grep) | | `test_code_consistency.cpp` |
| AUD-CC-12 | | yes | reviewer |
| AUD-CC-13 | yes (`using std::` at namespace scope) | | `test_code_consistency.cpp` |
| AUD-CC-14 | yes (forbidden-identifier grep) | | `test_code_consistency.cpp` |
| AUD-CC-15 | partial (float-literal regex) | yes | reviewer |
| AUD-CC-16 | yes | | `test_code_consistency.cpp` |
| AUD-CC-17 | partial | yes | reviewer |
| AUD-CC-18 | partial | yes | reviewer |
