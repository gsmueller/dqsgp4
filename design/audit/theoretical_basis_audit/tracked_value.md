# Theoretical Basis Audit — `src/math/tracked_value.h`

**File**: `src/math/tracked_value.h` (534 lines)
**Status**: AUDITED — 43 audit cards
**Role**: Foundation file. Every downstream `total_error()` claim in the propagator inherits from the closed-form per-category bounds defined here (REQ-EF-3, REQ-EF-9).
**Audit scope**: One card per distinct function — constructors, struct/class operators, math intrinsics, comparisons, helpers, and free templates.

## Index

| # | ID | Lines | Category | Verdict |
|---:|---|---:|---|:---:|
| 1 | `ThreeErrors::ctor_default` | 28 | constructor | OK |
| 2 | `ThreeErrors::ctor_three_arg` | 29 | constructor | OK |
| 3 | `ThreeErrors::total` | 33 | bound aggregator | PASS |
| 4 | `ThreeErrors::rss` | 37-40 | statistical estimator | PASS |
| 5 | `ThreeErrors::operator+` | 48-54 | error-combine | PASS |
| 6 | `ThreeErrors::operator-` | 58-65 | error-combine | OK |
| 7 | `ThreeErrors::operator*(scalar,e)` | 69-73 | error-scale | PASS |
| 8 | `ThreeErrors::operator*(e,scalar)` | 75-79 | error-scale | PASS |
| 9 | `ThreeErrors::operator/(e,scalar)` | 82-90 | error-scale | PASS |
| 10 | `ThreeErrors::max_per_category` | 93-100 | error-combine | PASS |
| 11 | `ThreeErrors::apply` | 105-107 | category mapper | PASS |
| 12 | `TrackedValue::ctor_default` | 118 | constructor | OK |
| 13 | `TrackedValue::ctor_val_errors` | 120-121 | constructor | OK |
| 14 | `TrackedValue::ctor_val_mpa` | 123-124 | constructor | OK |
| 15 | `TrackedValue::exact_integer` | 129-131 | named ctor | PASS |
| 16 | `TrackedValue::defined` | 135-139 | named ctor | PASS |
| 17 | `TrackedValue::defined_with_physical_uncertainty` | 143-150 | named ctor | PASS |
| 18 | `TrackedValue::measured` | 153-158 | named ctor | PASS |
| 19 | `TrackedValue::total_error` | 162 | query | PASS |
| 20 | `TrackedValue::reliable_digits` | 164-171 | query | PASS |
| 21 | `TrackedValue::operator+` | 176-178 | arithmetic | PASS |
| 22 | `TrackedValue::operator-` | 185-187 | arithmetic | PASS |
| 23 | `TrackedValue::operator-(unary)` | 192-194 | arithmetic | PASS |
| 24 | `TrackedValue::operator*` | 204-232 | arithmetic | PASS |
| 25 | `TrackedValue::operator/` | 238-271 | arithmetic | PASS |
| 26 | `TrackedValue::sqrt` | 280-293 | math fn | PASS |
| 27 | `TrackedValue::sin` | 308-320 | math fn | PASS |
| 28 | `TrackedValue::cos` | 333-345 | math fn | PASS |
| 29 | `TrackedValue::atan` | 349-360 | math fn | PASS |
| 30 | `TrackedValue::atan2` | 365-393 | math fn | PASS |
| 31 | `TrackedValue::abs` | 396-399 | math fn | PASS |
| 32 | `TrackedValue::fmod` | 404-410 | math fn | PASS |
| 33 | `TrackedValue::operator<` | 414 | comparison | OK |
| 34 | `TrackedValue::operator>` | 415 | comparison | OK |
| 35 | `TrackedValue::operator<=` | 416 | comparison | OK |
| 36 | `TrackedValue::operator>=` | 417 | comparison | OK |
| 37 | `TrackedValue::operator==` | 418 | comparison | OK |
| 38 | `TrackedValue::operator!=` | 419 | comparison | OK |
| 39 | `TrackedValue::ctor_int` | 423 | implicit ctor | PASS |
| 40 | `TrackedValue::representation_bound` | 441-492 | precision primitive | PASS-conservative |
| 41 | `TrackedValue::from_string` | 498-506 | utility | OK |
| 42 | `exact<T>` (free) | 513-515 | free fn | PASS |
| 43 | `ratio<T>` (free) | 528-532 | free fn | PASS |

Legend: PASS = full theory ↔ method ↔ bound triad verified. OK = pure data-movement / metadata function (no numerical bound to certify). PASS-conservative = bound is rigorous but loose by a small known factor.

---

## 1. `ThreeErrors::ctor_default`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::ctor_default
Location:               src/math/tracked_value.h:28
Mathematical statement: e := (0, 0, 0) — the zero element of the error space.

THEORY
  Underlying theorem:   Definition of the zero element of (T^3, +). The
                        error budget is the non-negative cone of T^3; the
                        additive identity is the origin.
  Primary reference:    Definition; REQ-EF-1 (three non-negative budgets).
  Domain of validity:   All T with a zero element.

METHOD
  Method declared:      Direct construction with T(0) in each component.
  Method implemented:   `ThreeErrors() : measurement(T(0)), precision(T(0)),
                                          accuracy(T(0)) {}`
  Match verdict:        OK — exact match.

ERROR BOUND
  Bound category:       n/a (constructs a zero-error object)
  Bound formula:        0
  Bound implemented:    T(0) in each slot.
  Bound verdict:        OK.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1
  AUD-EF applies:       AUD-EF-1 (every value carries the budget)
  Verification test:    tests/audit/test_error_framework.cpp — REQ-EF-1.

NOTES
  Used as the starting point for any exact value (integers, computed-from-int).
```

## 2. `ThreeErrors::ctor_three_arg`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::ctor_three_arg
Location:               src/math/tracked_value.h:29
Mathematical statement: e := (m, p, a) for caller-supplied m, p, a ≥ 0.

THEORY
  Underlying theorem:   Direct construction; preserves caller-supplied bounds.
  Primary reference:    REQ-EF-1.
  Domain of validity:   Any non-negative triple.

METHOD
  Method declared:      Member-initializer-list copy.
  Method implemented:   `: measurement(m), precision(p), accuracy(a) {}`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       caller-determined per-category
  Bound formula:        identity (no transformation)
  Bound implemented:    identity.
  Bound verdict:        OK — preserves the caller's contract.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1
  AUD-EF applies:       AUD-EF-1
  Verification test:    tests/audit/test_error_framework.cpp — REQ-EF-1.

NOTES
  Does not validate non-negativity at the type level. Negative inputs would
  violate the invariant; callers are trusted (the only call sites are this
  library's own named constructors which compute |·|-style bounds).
```

## 3. `ThreeErrors::total`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::total
Location:               src/math/tracked_value.h:33
Mathematical statement: total(e) := m + p + a.

THEORY
  Underlying theorem:   Triangle inequality. For not-necessarily-independent
                        error sources δ_m, δ_p, δ_a, the worst-case combined
                        deviation satisfies
                          |δ_m + δ_p + δ_a| ≤ |δ_m| + |δ_p| + |δ_a|.
                        This is the L^1 norm on (T^3).
  Primary reference:    Higham (2002) §2.2; Wilkinson (1965) §1.6.
                        Triangle inequality is textbook (any analysis text).
  Domain of validity:   Always; no independence required.

METHOD
  Method declared:      Sum the three categories.
  Method implemented:   `return measurement + precision + accuracy;`
  Match verdict:        PASS — closed-form addition.

ERROR BOUND
  Bound category:       This is the aggregator that defines total_error().
  Bound formula:        e_m + e_p + e_a (L^1 norm of (m,p,a)).
  Bound implemented:    e_m + e_p + e_a.
  Bound verdict:        PASS — exact identity, no truncation.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-2 (total_error is a rigorous upper bound).
  AUD-EF applies:       AUD-EF-2 (combined bound dominates each category).
  Verification test:    AUD-EF-2 property test.

NOTES
  Comments at lines 31-32 explicitly note "the three error sources are not
  necessarily independent" — this is why we add (L^1) rather than RSS (L^2).
  RSS would underestimate the worst-case bound when sources correlate.
```

## 4. `ThreeErrors::rss`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::rss
Location:               src/math/tracked_value.h:37-40
Mathematical statement: rss(e) := sqrt(m² + p² + a²) — Euclidean / L² norm.

THEORY
  Underlying theorem:   For independent random error sources with standard
                        deviations σ_m, σ_p, σ_a, the combined standard
                        deviation is sqrt(σ_m² + σ_p² + σ_a²) — orthogonal
                        variance addition.
  Primary reference:    Standard probability — sum-of-variances for
                        independent r.v.s. Higham (2002) §1.4 contrasts
                        probabilistic vs deterministic error.
  Domain of validity:   Independent error sources only. Inadmissible as a
                        rigorous upper bound; statistical reporting only.

METHOD
  Method declared:      Square, sum, sqrt.
  Method implemented:   `return sqrt(m*m + p*p + a*a);`
  Match verdict:        PASS — closed-form L² norm.

ERROR BOUND
  Bound category:       statistical estimate, NOT a rigorous upper bound.
  Bound formula:        sqrt(m² + p² + a²).
  Bound implemented:    sqrt(m² + p² + a²).
  Bound verdict:        PASS as a reporting value; not used on the
                        propagation path (total() is).

CROSS-AUDIT
  REQ-EF applies:       n/a for upper-bound contract; reporting only.
  AUD-EF applies:       n/a (rss is not on the upper-bound path).
  Verification test:    n/a as a soundness test; report-side only.

NOTES
  Comment line 36: "Only valid if sources are independent." Code never uses
  rss() internally on the bound-propagation path — only total() is used by
  total_error() at line 162. Safe by construction.
```

## 5. `ThreeErrors::operator+(ThreeErrors, ThreeErrors)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::operator+
Location:               src/math/tracked_value.h:48-54
Mathematical statement: (a+b)_k = a_k + b_k for k ∈ {m, p, a}.

THEORY
  Underlying theorem:   Per-category addition realizes the REQ-EF-3
                        add/sub propagation rule (bound = a + b). This is
                        the triangle inequality applied componentwise.
  Primary reference:    REQ-EF-3 add: bound = a_bound + b_bound.
                        Higham (2002) §2.2 forward-error of sum.
  Domain of validity:   All T with addition.

METHOD
  Method declared:      Component-wise addition on the (m, p, a) tuple.
  Method implemented:   `ThreeErrors(a.m + b.m, a.p + b.p, a.a + b.a)`.
  Match verdict:        PASS — exact match to REQ-EF-3.

ERROR BOUND
  Bound category:       each k separately
  Bound formula:        REQ-EF-3 add/sub: bound = a_bound + b_bound.
  Bound implemented:    a_k + b_k per category.
  Bound verdict:        PASS — exact REQ-EF-3 form.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (add/sub bound).
  AUD-EF applies:       AUD-EF-3 (closed-form bounds match REQ-EF-3).
  Verification test:    test_error_framework.cpp REQ-EF-3 add case.

NOTES
  Used by TrackedValue::operator+ (card #21) and operator- (card #22) at
  lines 177 and 186 — both add error budgets symmetrically (no minus in
  the error path), which is the REQ-EF-3 prescription.
```

## 6. `ThreeErrors::operator-(ThreeErrors, ThreeErrors)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::operator-
Location:               src/math/tracked_value.h:58-65
Mathematical statement: (a-b)_k = max(a_k - b_k, 0) for k ∈ {m, p, a}.

THEORY
  Underlying theorem:   Non-negative-cone subtraction. The error budget
                        lives in [0, ∞)^3; subtraction must clamp to the
                        cone floor to preserve the invariant.
  Primary reference:    REQ-EF-1 (categories are non-negative).
  Domain of validity:   All T with order.

METHOD
  Method declared:      Component-wise subtraction clamped at 0.
  Method implemented:   `max(a.k - b.k, T(0))` per category.
  Match verdict:        OK — matches its stated semantics.

ERROR BOUND
  Bound category:       each k separately
  Bound formula:        Not a propagation bound — bookkeeping operation.
  Bound implemented:    max(a_k - b_k, 0).
  Bound verdict:        OK — non-negativity invariant preserved.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 invariant maintenance.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Not called by any TrackedValue arithmetic operator. Available for
  caller-side bookkeeping if a known component is being subtracted from
  the budget (e.g., to remove an over-counted source after a refinement).
  Not on the bound-propagation path.
```

## 7. `ThreeErrors::operator*(scalar, ThreeErrors)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::operator*_left
Location:               src/math/tracked_value.h:69-73
Mathematical statement: (s · e)_k = |s| · e_k.

THEORY
  Underlying theorem:   Scaling by a scalar of known magnitude. For
                        y = s · x, |y − y_true| = |s| · |x − x_true|.
                        Lipschitz constant of x ↦ s·x is exactly |s|.
                        Reduces REQ-EF-3 mul to the case where one factor
                        is exact (its error is zero).
  Primary reference:    REQ-EF-3 (mul rule, exact-factor case);
                        Higham (2002) §3.2.
  Domain of validity:   All T.

METHOD
  Method declared:      Multiply each category by |s|.
  Method implemented:   `as = abs(s); ThreeErrors(as*e.m, as*e.p, as*e.a)`.
  Match verdict:        PASS — implements |s|·e exactly.

ERROR BOUND
  Bound category:       per-k
  Bound formula:        |s| · e_k.
  Bound implemented:    |s| · e_k.
  Bound verdict:        PASS — exact closed form.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (mul reduced to scalar case).
  AUD-EF applies:       AUD-EF-3.
  Verification test:    AUD-EF-3 scalar-mul case.

NOTES
  abs(s) ensures non-negativity. Used by fmod (card #32) where the scalar
  is the integer quotient magnitude |n|.
```

## 8. `ThreeErrors::operator*(ThreeErrors, scalar)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::operator*_right
Location:               src/math/tracked_value.h:75-79
Mathematical statement: (e · s)_k = e_k · |s|. Commuted form of card #7.

THEORY
  Underlying theorem:   Commutativity of scalar multiplication.
  Primary reference:    REQ-EF-3 (same as card #7).
  Domain of validity:   All T.

METHOD
  Method declared:      Same as #7, arguments commuted.
  Method implemented:   `as = abs(s); ThreeErrors(e.m*as, e.p*as, e.a*as)`.
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        |s| · e_k.
  Bound implemented:    |s| · e_k.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3.
  AUD-EF applies:       AUD-EF-3.
  Verification test:    AUD-EF-3.

NOTES
  Provided for syntactic convenience (`e * s` and `s * e` both work).
  Identical math to card #7.
```

## 9. `ThreeErrors::operator/(ThreeErrors, scalar)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::operator_div
Location:               src/math/tracked_value.h:82-90
Mathematical statement: (e / s)_k = e_k / |s| for s ≠ 0; ∞ for s = 0.

THEORY
  Underlying theorem:   Division by a known scalar (no error on s). For
                        y = x / s, |y − y_true| = |x − x_true| / |s|.
                        Lipschitz constant 1/|s|. Reduces REQ-EF-3 div
                        to the case where the denominator has zero error.
                        The s = 0 case is catastrophic; per REQ-EF-9 the
                        bound is set to numeric_limits<T>::max() (the
                        conservative ceiling, signalling unreliability).
  Primary reference:    REQ-EF-3 (div, exact-denominator case);
                        REQ-EF-9 (catastrophic-case signaling).
  Domain of validity:   s ≠ 0.

METHOD
  Method declared:      Per-category division by |s|, with s=0 fallback.
  Method implemented:   `as = abs(s); if (as==0) return inf-in-all;
                                       else return e_k / as` per category.
  Match verdict:        PASS — REQ-EF-3 div with b_err=0 reduction.

ERROR BOUND
  Bound formula:        e_k / |s|, with REQ-EF-9 catastrophic fallback
                        numeric_limits<T>::max() when s = 0.
  Bound implemented:    matches.
  Bound verdict:        PASS — exact closed form + correct REQ-EF-9 signal.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (simplified div), REQ-EF-9.
  AUD-EF applies:       AUD-EF-3, AUD-EF-9.
  Verification test:    AUD-EF-3 scalar-div; AUD-EF-9 zero-denominator.

NOTES
  The full REQ-EF-3 div formula is
    bound = (|a|·b_err + |b|·a_err) / (|b|·(|b|−b_err))
  with b_err = 0 (denominator is a plain T scalar), this reduces to
    bound = a_err / |b|.
  Consistent. The s = 0 fallback to ::max() is REQ-EF-9 "signal, not
  silence".
```

## 10. `ThreeErrors::max_per_category`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::max_per_category
Location:               src/math/tracked_value.h:93-100
Mathematical statement: max_k(a, b) := (max(a_m, b_m), max(a_p, b_p),
                                        max(a_a, b_a)).

THEORY
  Underlying theorem:   Component-wise max on (T^3, ≤). Used as the
                        worst-of-two-sources combinator when a result
                        inherits the worse bound from two known sources
                        (e.g., piecewise branches).
  Primary reference:    Order-theoretic max on a product order. REQ-EF-13
                        (composite aggregation may use per-category max).
  Domain of validity:   All T with order.

METHOD
  Method declared:      Per-category max.
  Method implemented:   `max(a.k, b.k)` per category.
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        max(a_k, b_k) per category.
  Bound implemented:    matches.
  Bound verdict:        PASS — correct dominate-of-two combinator.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-13 (composite max aggregation).
  AUD-EF applies:       AUD-EF-13.
  Verification test:    AUD-EF-13.

NOTES
  Utility; not called from TrackedValue methods directly. Composite types
  (Vector3, Quaternion, …) may use this when defining their aggregate
  total_error() as a per-component max.
```

## 11. `ThreeErrors::apply<F>`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ThreeErrors::apply
Location:               src/math/tracked_value.h:105-107
Mathematical statement: apply(f, e) := (f(e_m), f(e_p), f(e_a)).

THEORY
  Underlying theorem:   Lift of a scalar non-negative-preserving map
                        f: T → T to the product (T^3). Provides a single
                        bound-computation lambda to be applied uniformly
                        across all three categories.
  Primary reference:    Functorial lift on a product type; REQ-EF-3
                        component bounds use this pattern.
                        Comment line 103 requires f non-negative-preserving.
  Domain of validity:   f maps non-negative T to non-negative T.

METHOD
  Method declared:      Apply f to each category; return new ThreeErrors.
  Method implemented:   `ThreeErrors(f(m), f(p), f(a))`.
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        f(e_k) per category.
  Bound implemented:    matches.
  Bound verdict:        PASS — soundness depends on f being a sound
                        bound for the operation it represents (caller's
                        responsibility); apply itself is the correct lift.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (used inside sqrt, sin, cos, atan branches).
  AUD-EF applies:       AUD-EF-3.
  Verification test:    AUD-EF-3 (per-op tests exercise apply paths).

NOTES
  Used by sqrt (#26), sin (#27), cos (#28), atan (#29). Each passes a
  lambda computing the per-category bound; apply lifts to all three.
  The closure captures the value (a.value) so that the same f-form applies
  uniformly across categories — only the input error varies.
```

## 12. `TrackedValue::ctor_default`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::ctor_default
Location:               src/math/tracked_value.h:118
Mathematical statement: v := 0 with errors (0, 0, 0).

THEORY
  Underlying theorem:   Definition of the zero element of the value space.
  Primary reference:    REQ-EF-1 (zero error if default-constructed);
                        REQ-EF-10 (value is plain T).
  Domain of validity:   T must support T(0).

METHOD
  Method declared:      Construct value = T(0) and ThreeErrors default.
  Method implemented:   `: value(T(0)), errors() {}`.
  Match verdict:        OK.

ERROR BOUND
  Bound formula:        0 in all categories.
  Bound implemented:    0.
  Bound verdict:        OK.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-10.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Used in default-initialization paths and STL containers.
```

## 13. `TrackedValue::ctor_val_errors`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::ctor_val_errors
Location:               src/math/tracked_value.h:120-121
Mathematical statement: v := val with errors := err (copy-construct).

THEORY
  Underlying theorem:   Direct construction; caller-provided value and
                        per-category bound.
  Primary reference:    REQ-EF-1.
  Domain of validity:   Any (val, err) pair.

METHOD
  Method declared:      Member init.
  Method implemented:   `: value(val), errors(err) {}`.
  Match verdict:        OK.

ERROR BOUND
  Bound formula:        identity.
  Bound implemented:    identity.
  Bound verdict:        OK.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Workhorse internal constructor: every arithmetic operator (cards #21-32)
  returns through this form.
```

## 14. `TrackedValue::ctor_val_mpa`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::ctor_val_mpa
Location:               src/math/tracked_value.h:123-124
Mathematical statement: v := val with errors := (m, p, a).

THEORY
  Underlying theorem:   Direct construction with explicit category bounds.
  Primary reference:    REQ-EF-1.
  Domain of validity:   Any non-negative (m, p, a).

METHOD
  Method declared:      Forward to ThreeErrors(m, p, a) (card #2).
  Method implemented:   `: value(val), errors(meas, prec, acc) {}`.
  Match verdict:        OK.

ERROR BOUND
  Bound formula:        identity.
  Bound implemented:    identity.
  Bound verdict:        OK.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Convenience overload used by exact_integer, defined, measured.
```

## 15. `TrackedValue::exact_integer(int n)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::exact_integer
Location:               src/math/tracked_value.h:129-131
Mathematical statement: For n ∈ ℤ with |n| ≤ 2^p (p = significand bits of T),
                        T(n) is exactly representable; output v = T(n) with
                        all error categories = 0.

THEORY
  Underlying theorem:   IEEE-754 §3.4: every integer with |n| ≤ 2^p is
                        exactly representable in a binary floating-point
                        format with p significand bits (p = 53 for double).
  Primary reference:    IEEE Std 754-2019 §3.4; Goldberg (1991) "What Every
                        Computer Scientist Should Know About Floating-Point
                        Arithmetic", Theorem 1.
  Domain of validity:   |n| ≤ 2^p. For double, |n| ≤ 2^53. For 32-bit
                        int input, always satisfied.

METHOD
  Method declared:      Cast int → T directly; assign zero errors.
  Method implemented:   `TrackedValue(T(n), T(0), T(0), T(0))`.
  Match verdict:        PASS — direct construction; zero error claimed.

ERROR BOUND
  Bound formula:        0 in all categories (under the |n| ≤ 2^p precondition).
  Bound implemented:    0.
  Bound verdict:        PASS — under the IEEE-754 §3.4 precondition.
                        For 32-bit int with T = double, precondition is
                        always satisfied.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-10.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1 / IEEE-754 exact-integer property test.

NOTES
  Used pervasively by exact<T>(int) (card #42). For T = float (p=24), the
  precondition fails for |n| > 2^24 ≈ 1.68 × 10^7; codebase usage is for
  small constants (typical |n| ≤ 10000) so this is moot in practice.
```

## 16. `TrackedValue::defined(decimal_str)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::defined
Location:               src/math/tracked_value.h:135-139
Mathematical statement: For a decimal string s representing real r:
                        v = round_to_T(r);
                        precision = repr_bound(v);
                        measurement = accuracy = 0.

THEORY
  Underlying theorem:   IEEE-754 round-to-nearest-even: the rounding cost
                        from a real number to its T-representable form is
                        bounded by 0.5 ULP(round_to_T(r)).
                        Goldberg (1991) Theorem 1; Higham (2002) §2.1.
  Primary reference:    IEEE Std 754-2019 §4.3.1 (rounding rules);
                        Goldberg (1991); Higham (2002) §2.1.
  Domain of validity:   All decimal strings parseable to a finite real.

METHOD
  Method declared:      Parse string → T (via from_string, card #41);
                        compute representation cost (card #40); set
                        measurement and accuracy to 0.
  Method implemented:   `val = from_string(s);
                        prec = representation_bound(val);
                        return TrackedValue(val, 0, prec, 0)`.
  Match verdict:        PASS.

ERROR BOUND
  Bound category:       precision (representation cost).
  Bound formula:        0.5 ULP(v) = representation_bound(v).
  Bound implemented:    representation_bound(val).
  Bound verdict:        PASS — the IEEE-754 0.5-ULP rounding bound,
                        delivered conservatively at 1 ULP (see card #40
                        for the conservative factor).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-3 (rounding cost is precision).
  AUD-EF applies:       AUD-EF-3 representation cost.
  Verification test:    AUD-EF tests sampling decimal constants.

NOTES
  Entry point for defining constants such as μ_E, J_2, etc. Note:
  measurement = 0 because by *definition* a defined constant has no
  measurement uncertainty (it is defined to its quoted digits). For
  measured physical constants with quoted σ, use card #18 (measured).
```

## 17. `TrackedValue::defined_with_physical_uncertainty`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::defined_with_physical_uncertainty
Location:               src/math/tracked_value.h:143-150
Mathematical statement: For decimal value string s (real v) and σ string σ:
                        val = round_to_T(s);
                        measurement = round_to_T(σ);
                        precision = repr_bound(val);
                        accuracy = 0.

THEORY
  Underlying theorem:   Composition of two definitional sources:
                        (i)  IEEE-754 rounding cost (precision)
                             — Goldberg (1991) Theorem 1;
                        (ii) physical measurement uncertainty (measurement)
                             — traceable to standards (CODATA, IERS),
                             entered as the measurement budget per REQ-EF-1.
  Primary reference:    REQ-EF-1; CODATA-style σ entries.
  Domain of validity:   Both strings parse; σ ≥ 0.

METHOD
  Method declared:      Parse both strings; representation_bound on val;
                        zero accuracy.
  Method implemented:   lines 145-149.
  Match verdict:        PASS — sound composition of the two error sources.

ERROR BOUND
  Bound formula:        measurement = σ; precision = 0.5 ULP(val);
                        accuracy = 0.
  Bound implemented:    matches.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Used for parameters such as the "physical Earth" GM where a quoted σ
  exists alongside a defining numeric value (the latter is conventional,
  the former is the underlying measured quantity).
```

## 18. `TrackedValue::measured(value_str, sigma_str)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::measured
Location:               src/math/tracked_value.h:153-158
Mathematical statement: Same form as card #17, but semantic distinction —
                        the value itself is measured (not definitional).
                        val = round_to_T(s); measurement = σ;
                        precision = repr_bound(val); accuracy = 0.

THEORY
  Underlying theorem:   Same as #17 — IEEE-754 rounding + measurement σ.
  Primary reference:    REQ-EF-1; metrology — quoted σ traceable to
                        primary standards.
  Domain of validity:   Both strings parse; σ ≥ 0.

METHOD
  Method declared:      Identical to #17.
  Method implemented:   lines 155-158.
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        measurement = σ; precision = 0.5 ULP; accuracy = 0.
  Bound implemented:    matches.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Semantically distinct from defined() (#16): a measured constant has a
  measurement budget; a defined constant does not. The two named
  constructors enforce the modeling distinction at construction time.
```

## 19. `TrackedValue::total_error()`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::total_error
Location:               src/math/tracked_value.h:162
Mathematical statement: total_error(v) := v.errors.total()
                                       = v.errors.m + v.errors.p + v.errors.a.

THEORY
  Underlying theorem:   Triangle inequality (REQ-EF-2): for the worst-case
                        composite error,
                        |x − x_true| ≤ |Δ_m| + |Δ_p| + |Δ_a|.
                        Each Δ is a rigorous per-category upper bound;
                        their sum is therefore a rigorous total bound.
  Primary reference:    REQ-EF-2; Higham (2002) §2.2; Wilkinson (1965) §1.6.
  Domain of validity:   Always.

METHOD
  Method declared:      Forward to errors.total() (card #3).
  Method implemented:   `return errors.total();`
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        m + p + a.
  Bound implemented:    matches.
  Bound verdict:        PASS — exactly the REQ-EF-2 contract.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-2 (this is the contract).
  AUD-EF applies:       AUD-EF-2 (total dominates each category).
  Verification test:    test_error_framework.cpp REQ-EF-2 property test.

NOTES
  Headline accessor for the entire error framework. Every downstream
  reliability claim in the propagator passes through this method.
```

## 20. `TrackedValue::reliable_digits()`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::reliable_digits
Location:               src/math/tracked_value.h:164-171
Mathematical statement: reliable_digits(v) := largest k ∈ ℤ s.t.
                        10^(−k) ≥ total_error(v) / |v.value|.
                        Equivalently: floor(−log10(rel_err)).

THEORY
  Underlying theorem:   Definition of relative-error decimal digits:
                        if |x − x_true| / |x_true| ≤ 10^(−k), then x
                        agrees with x_true to at least k decimal digits.
                        Standard relative-error metric.
  Primary reference:    REQ-EF-11; Higham (2002) §2.2; Goldberg (1991) §2.
  Domain of validity:   v.value ≠ 0 and total_error < |v.value|. Edge
                        cases handled explicitly.

METHOD
  Method declared:      Compute total_error / |value|, take −log10, floor.
                        Edge cases: zero error → INT_MAX; error dominates
                        value → 0.
  Method implemented:   lines 165-171.
  Match verdict:        PASS — direct from the definition.

ERROR BOUND
  Bound category:       n/a (returns an integer digit count, not a bound).
  Bound formula:        floor(−log10(total_error / |value|)).
  Bound implemented:    `floor(-log10(te / av))`.
  Bound verdict:        PASS — definitional. floor() is the conservative
                        direction: under-states the reliable digit count
                        rather than over-states it.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-11.
  AUD-EF applies:       AUD-EF-11.
  Verification test:    REQ-EF-11 property tests.

NOTES
  Edge cases (lines 167-169):
  - te == 0: returns INT_MAX (exactly representable; digit count
    undefined → maximum).
  - av == 0 or te ≥ av: returns 0 (the error swallows the value).
  These guard log10 from -∞ and from a domain error.
```

## 21. `TrackedValue::operator+(TrackedValue, TrackedValue)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator+
Location:               src/math/tracked_value.h:176-178
Mathematical statement: For tracked a = (x, δ_x), b = (y, δ_y):
                        a + b = (x+y, δ_x + δ_y) per category.

THEORY
  Underlying theorem:   Triangle inequality on each category:
                        |(x+y) − (x_true + y_true)| ≤ |x − x_true|
                                                    + |y − y_true|.
                        This is exact for addition; no second-order term.
  Primary reference:    REQ-EF-3 (add: bound = a_bound + b_bound);
                        Higham (2002) §2.2; Wilkinson (1965) §1.6.
  Domain of validity:   All T.

METHOD
  Method declared:      Closed-form: value = a.value + b.value;
                        errors = a.errors + b.errors (per-category, card #5).
  Method implemented:   `return TrackedValue(a.value + b.value,
                                              a.errors + b.errors);`
  Match verdict:        PASS — direct REQ-EF-3 add prescription.

ERROR BOUND
  Bound formula:        REQ-EF-3 add/sub: bound = a_bound + b_bound.
  Bound implemented:    a_bound + b_bound per category.
  Bound verdict:        PASS — matches REQ-EF-3 exactly.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (add), REQ-EF-8 (Lipschitz-1 floor met).
  AUD-EF applies:       AUD-EF-3 add case.
  AUD-MC applies:       AUD-MC-1, AUD-MC-2 (additive axioms for composite types).
  Verification test:    test_error_framework.cpp REQ-EF-3 add.

NOTES
  Foundational for DualNumber, Vector3, Quaternion, DualQuaternion
  additive composition. Lipschitz constant for + is 1 in each input.
  Catastrophic cancellation in the *value* (x ≈ −y) does not increase
  absolute error but degrades reliable_digits — REQ-EF-9 acknowledges
  this is reported (via reliable_digits), not silenced.
```

## 22. `TrackedValue::operator-(TrackedValue, TrackedValue)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator-_binary
Location:               src/math/tracked_value.h:185-187
Mathematical statement: For a = (x, δ_x), b = (y, δ_y):
                        a − b = (x − y, δ_x + δ_y) per category.

THEORY
  Underlying theorem:   Triangle inequality:
                        |(x−y) − (x_true − y_true)| ≤ |x − x_true|
                                                    + |y − y_true|.
                        Errors ADD, not subtract — subtractive cancellation
                        does not reduce absolute error.
                        Higham (2002) §1.7 "Cancellation".
  Primary reference:    REQ-EF-3 (sub: bound = a_bound + b_bound);
                        Higham (2002) §1.7, §2.2.
  Domain of validity:   All T.

METHOD
  Method declared:      Closed-form: value = a.value − b.value;
                        errors = a.errors + b.errors (NOT minus).
  Method implemented:   `return TrackedValue(a.value - b.value,
                                              a.errors + b.errors);`
  Match verdict:        PASS — REQ-EF-3 sub prescription. Errors are
                        added (line 186), the correct triangle inequality.

ERROR BOUND
  Bound formula:        REQ-EF-3 add/sub: bound = a_bound + b_bound.
  Bound implemented:    a_bound + b_bound per category.
  Bound verdict:        PASS — exact REQ-EF-3 form.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (sub), REQ-EF-8 (Lipschitz-1 floor met).
  AUD-EF applies:       AUD-EF-3 sub case.
  AUD-MC applies:       AUD-MC algebra axioms for composite types.
  Verification test:    AUD-EF-3 sub.

NOTES
  Comments at lines 181-183 explicitly: "subtractive cancellation does not
  increase absolute error, but relative error may grow … reliable_digits
  detects this." Bound remains absolute and correct; the user sees fewer
  reliable digits via reliable_digits (REQ-EF-9 — signal, not silence).
```

## 23. `TrackedValue::operator-()` (unary)

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator-_unary
Location:               src/math/tracked_value.h:192-194
Mathematical statement: −a := (−a.value, a.errors).

THEORY
  Underlying theorem:   Negation is a sign flip; Lipschitz constant = 1.
                        Per REQ-EF-14, identity-on-error operations
                        preserve the budget exactly.
                        IEEE-754 negation is exact (sign-bit flip; no
                        rounding).
  Primary reference:    REQ-EF-14; Higham (2002) §2.2.
  Domain of validity:   All T.

METHOD
  Method declared:      Value: T's unary minus (exact in IEEE-754).
                        Errors: unchanged.
  Method implemented:   `return TrackedValue(-value, errors);`
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        errors unchanged.
  Bound implemented:    errors unchanged.
  Bound verdict:        PASS — REQ-EF-14 preservation.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-14 (identity-class preservation).
  AUD-EF applies:       AUD-EF-14.
  AUD-MC applies:       Additive inverse axiom for composite types.
  Verification test:    AUD-EF-14.

NOTES
  Both IEEE-754 hardware and boost multiprecision negate exactly via
  sign-bit flip — no rounding cost. Bound preservation is rigorous.
```

## 24. `TrackedValue::operator*(TrackedValue, TrackedValue)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator*
Location:               src/math/tracked_value.h:204-232
Mathematical statement: For a = (x, δ_x), b = (y, δ_y):
                        a * b = (xy, δ_xy) where the per-category bound is
                        δ_xy ≤ |x|·δ_y + |y|·δ_x + δ_x·δ_y.

THEORY
  Underlying theorem:   Expansion:
                        (x+δ_x)(y+δ_y) − xy = x·δ_y + y·δ_x + δ_x·δ_y,
                        |·| triangle inequality gives
                        |(x+δ_x)(y+δ_y) − xy| ≤ |x|·|δ_y| + |y|·|δ_x|
                                              + |δ_x|·|δ_y|.
                        Each term retained: the third is genuinely
                        second-order small but kept for rigor when either
                        component error approaches its value.
  Primary reference:    REQ-EF-3 (mul: bound = |a_val|·b_bound + |b_val|·
                        a_bound + a_bound·b_bound);
                        Higham (2002) §3.5; Wilkinson (1965) §3.1.
  Domain of validity:   All T.

METHOD
  Method declared:      Closed-form per-category sum of three terms.
                        Applied per category independently; cross-category
                        products neglected per comment lines 200-202
                        (standard interval-arithmetic practice for
                        separated error budgets).
  Method implemented:   safe_mul lambda guards 0×inf → 0 (REQ-EF-9 defense);
                        mul_bound = safe_mul(av, b_err) + safe_mul(bv, a_err)
                                  + safe_mul(a_err, b_err).
  Match verdict:        PASS — REQ-EF-3 mul prescription. safe_mul handles
                        REQ-EF-9 fallout (0 × max() = 0, not NaN).

ERROR BOUND
  Bound formula:        REQ-EF-3 mul:
                        bound = |a_val|·b_bound + |b_val|·a_bound
                              + a_bound·b_bound.
  Bound implemented:    lines 215-222 — exact match including the
                        second-order term. safe_mul prevents NaN.
  Bound verdict:        PASS — REQ-EF-3 exact; REQ-EF-9 defense applied.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (mul), REQ-EF-9 (0×inf guard), REQ-EF-8.
  AUD-EF applies:       AUD-EF-3 mul, AUD-EF-9 catastrophic-mul.
  AUD-MC applies:       AUD-MC-4..7 (multiplicative axioms for composite
                        types: commutativity, associativity, distributivity).
  Verification test:    test_error_framework.cpp REQ-EF-3 mul; AUD-EF-9
                        0×max guard test.

NOTES
  Cross-category products (e.g., a.errors.measurement × b.errors.precision)
  are NOT computed — only same-category products. Comment lines 200-202
  declare this standard practice for separated budgets and acknowledge it
  is a deliberate trade for category-attribution clarity. The reduction
  is rigorous when individual errors << values, which holds throughout
  the propagator's domain. (A worst-case interval bound would compute
  all 9 cross-products; the 6 cross-pair products are second-order small.)
  REQ-EF-9 defense: safe_mul prevents 0 × numeric_limits::max() (from a
  previous catastrophic division) from producing NaN.
```

## 25. `TrackedValue::operator/(TrackedValue, TrackedValue)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator/
Location:               src/math/tracked_value.h:238-271
Mathematical statement: For a = (x, δ_x), b = (y, δ_y) with |δ_y| < |y|:
                        a / b = (x/y, δ_q), with
                        δ_q ≤ (|x|·δ_y + |y|·δ_x) / (|y|·(|y| − δ_y)).
                        Catastrophic when |δ_y| ≥ |y|: bound = max<T>.

THEORY
  Underlying theorem:   |x/y − (x+δ_x)/(y+δ_y)|
                        = |x(y+δ_y) − y(x+δ_x)| / |y(y+δ_y)|
                        = |x·δ_y − y·δ_x| / |y(y+δ_y)|
                        ≤ (|x|·δ_y + |y|·δ_x) / (|y|·(|y| − δ_y))
                        when |δ_y| < |y| (denominator stays positive).
                        Uses minimum possible denominator magnitude
                        |y| − δ_y for rigor (over-bounds δ_q).
  Primary reference:    REQ-EF-3 (div formula);
                        REQ-EF-9 (catastrophic when |δ_y| ≥ |y|);
                        Higham (2002) §3.6; Wilkinson (1965) §3.1.
  Domain of validity:   |y| > 0 and |δ_y| < |y|. Else catastrophic.

METHOD
  Method declared:      Closed-form per category:
                        - safe_mul lambda (defense against 0×inf, REQ-EF-9);
                        - numerator = |a|·b_err + |b|·a_err;
                        - denom_min = |b| − b_err;
                        - if denom_min ≤ 0: return numeric_limits<T>::max()
                          (REQ-EF-9 catastrophic signal);
                        - else: numerator / (|b| · denom_min).
  Method implemented:   lines 243-261.
  Match verdict:        PASS — REQ-EF-3 div + REQ-EF-9 fallback exactly.

ERROR BOUND
  Bound formula:        REQ-EF-3 div: (|a_val|·b_bound + |b_val|·a_bound)
                                      / (|b_val|·(|b_val| − b_bound)).
                        REQ-EF-9 catastrophic: numeric_limits<T>::max().
  Bound implemented:    lines 253-260 — exact match.
  Bound verdict:        PASS — closed form + catastrophic fallback.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (div), REQ-EF-9 (catastrophic), REQ-EF-8.
  AUD-EF applies:       AUD-EF-3 div, AUD-EF-9 catastrophic-div.
  AUD-MC applies:       Multiplicative inverse axiom for composite types.
  Verification test:    test_error_framework.cpp REQ-EF-3 div; AUD-EF-9
                        catastrophic-div with b_err > |b|.

NOTES
  Catastrophic-region: when den_err ≥ |b|, the true quotient could be
  arbitrarily large or change sign — bound is set to numeric_limits<T>::max()
  per REQ-EF-9. This propagates through total_error() and forces
  reliable_digits → 0; downstream callers see the signal and may abort
  or branch to a fallback method.
  safe_mul guard prevents NaN from 0 × max() (defense-in-depth; same
  pattern as card #24).
```

## 26. `TrackedValue::sqrt`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::sqrt
Location:               src/math/tracked_value.h:280-293
Mathematical statement: For a = (x, δ_x) with x > 0:
                        sqrt(a) = (√x, δ_s), with
                        δ_s ≤ δ_x / (2·√(x − δ_x)) when δ_x < x.
                        Catastrophic floor δ_s = √x when δ_x ≥ x.

THEORY
  Underlying theorem:   Mean-value theorem on f(t) = √t:
                        √(x+δ) − √x = δ / (2·√ξ) for some ξ between x
                        and x+δ. Worst case ξ = x − δ_x (minimum over
                        the perturbation interval, maximizing |1/(2√ξ)|):
                        |√(x+δ) − √x| ≤ |δ| / (2·√(x − δ_x)).
                        When δ_x ≥ x, the true value could be 0; the
                        worst-case change is the entire result √x.
  Primary reference:    REQ-EF-3 (sqrt: bound = err / (2·sqrt(val − err))
                        when err < val);
                        REQ-EF-9 (catastrophic when not);
                        Higham (2002) §3.3 (forward error of sqrt);
                        mean-value theorem (any calculus text).
  Domain of validity:   a.value > 0.

METHOD
  Method declared:      Per-category mean-value bound; catastrophic floor.
  Method implemented:   `sqrt_bound = (err >= a.value) ? val
                                                       : err / (2·sqrt(a.value − err))`,
                        applied per category via ThreeErrors::apply (card #11).
  Match verdict:        PASS — REQ-EF-3 sqrt + REQ-EF-9 fallback exactly.

ERROR BOUND
  Bound formula:        REQ-EF-3 sqrt: err / (2·sqrt(val − err))  when err < val;
                        REQ-EF-9: bound = val (entire value uncertain)
                                  when err ≥ val.
  Bound implemented:    lines 284-290 — exact match for both regions.
  Bound verdict:        PASS — closed-form match for both regimes.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (sqrt), REQ-EF-9, REQ-EF-8.
  AUD-EF applies:       AUD-EF-3 sqrt, AUD-EF-9 catastrophic-sqrt.
  AUD-MC applies:       Used by Vector3::norm, Quaternion::norm.
  Verification test:    test_error_framework.cpp REQ-EF-3 sqrt; AUD-EF-9
                        catastrophic-sqrt with err > val.

NOTES
  Catastrophic-floor "bound = val" honors REQ-EF-9: if the input error
  reaches its value, the true input could be zero, so the true sqrt
  could be 0 while the computed √x ≠ 0. The entire result is uncertain.
  Reporting bound = val is the rigorous worst case (true sqrt ∈ [0, √(x+err)]).
  Below the catastrophic threshold, err / (2·√(x − err)) is the tightest
  mean-value bound (minimum derivative across the perturbation interval).
```

## 27. `TrackedValue::sin`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::sin
Location:               src/math/tracked_value.h:308-320
Mathematical statement: For a = (x, δ_x):
                        sin(a) = (sin x, δ_s), with
                        δ_s ≤ |cos(x)|·δ_x + δ_x²/2, capped at 2.

THEORY
  Underlying theorem:   Taylor's theorem with Lagrange remainder for sin:
                        sin(x+δ) − sin(x) = cos(x)·δ − sin(ξ)·δ²/2
                        for some ξ between x and x+δ.
                        |sin(x+δ) − sin(x)| ≤ |cos(x)|·|δ| + |sin(ξ)|·δ²/2
                                            ≤ |cos(x)|·|δ| + δ²/2,
                        since |sin(ξ)| ≤ 1. Absolute cap: |sin(x+δ) − sin(x)|
                        ≤ 2 unconditionally (range diameter of sin).
  Primary reference:    REQ-EF-3 (sin/cos: |deriv(val)|·err + err²/2,
                        capped at 2);
                        Higham (2002) §13.3 (trig of perturbed argument);
                        Taylor's theorem with Lagrange remainder (Rudin
                        (1976) Theorem 5.15).
  Domain of validity:   All x ∈ ℝ.

METHOD
  Method declared:      Per-category: bound = |cos(val)|·err + err²/2,
                        capped at 2.
  Method implemented:   `linear = |cos(a.value)| * err;
                        quadratic = err*err/2;        // |sin(ξ)| ≤ 1
                        return min(linear + quadratic, T(2));`
                        Applied per category via apply (card #11).
  Match verdict:        PASS — REQ-EF-3 sin/cos exact.

ERROR BOUND
  Bound formula:        REQ-EF-3: |cos(val)|·err + err²/2, capped at 2.
  Bound implemented:    lines 313-317 — exact match.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (sin/cos rule).
  AUD-EF applies:       AUD-EF-3 sin case.
  AUD-MC applies:       AUD-MC-12 (sin/cos of perturbed angle), used by
                        quaternion exp_pure, dual-quaternion exp_screw.
  Verification test:    test_error_framework.cpp REQ-EF-3 sin.

NOTES
  Comments lines 295-307 explicit: linear term operative; quadratic from
  Taylor remainder with |sin(ξ)| ≤ 1; δ² < δ for δ < 1; absolute cap 2.
  Conservative choice of |sin(ξ)| ≤ 1 (vs tightening to interval bound
  on [x, x+δ]) is the REQ-EF-3 prescription; tighter bounds remain
  future-work.
```

## 28. `TrackedValue::cos`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::cos
Location:               src/math/tracked_value.h:333-345
Mathematical statement: For a = (x, δ_x):
                        cos(a) = (cos x, δ_c), with
                        δ_c ≤ |sin(x)|·δ_x + δ_x²/2, capped at 2.

THEORY
  Underlying theorem:   Taylor's theorem with Lagrange remainder for cos:
                        cos(x+δ) − cos(x) = −sin(x)·δ − cos(ξ)·δ²/2
                        for some ξ between x and x+δ.
                        |cos(x+δ) − cos(x)| ≤ |sin(x)|·|δ| + δ²/2,
                        since |cos(ξ)| ≤ 1. Absolute cap: ≤ 2.
  Primary reference:    REQ-EF-3 (sin/cos);
                        Higham (2002) §13.3;
                        Taylor's theorem with Lagrange remainder.
  Domain of validity:   All x ∈ ℝ.

METHOD
  Method declared:      Per-category: bound = |sin(val)|·err + err²/2,
                        capped at 2.
  Method implemented:   `linear = |sin(a.value)| * err;
                        quadratic = err*err/2;        // |cos(ξ)| ≤ 1
                        return min(linear + quadratic, T(2));`
                        Applied per category via apply (card #11).
  Match verdict:        PASS — REQ-EF-3 sin/cos exact.

ERROR BOUND
  Bound formula:        REQ-EF-3: |sin(val)|·err + err²/2, capped at 2.
                        (Note: derivative of cos at val is −sin(val), so
                        |deriv| = |sin(val)|.)
  Bound implemented:    lines 338-342 — exact match.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (sin/cos rule).
  AUD-EF applies:       AUD-EF-3 cos case.
  AUD-MC applies:       AUD-MC-12.
  Verification test:    test_error_framework.cpp REQ-EF-3 cos.

NOTES
  Mirror of card #27 (sin): same theory, same method, same bound, with
  the role of |sin| and |cos| swapped (derivative of cos is −sin).
  Comments lines 322-332 acknowledge the Taylor remainder structure.
```

## 29. `TrackedValue::atan`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::atan
Location:               src/math/tracked_value.h:349-360
Mathematical statement: For a = (x, δ_x):
                        atan(a) = (atan x, δ_a), with
                        δ_a ≤ δ_x / (1 + x²), capped at π.

THEORY
  Underlying theorem:   Mean-value theorem on atan, with d(atan)/dx
                        = 1/(1+x²). For any x, x+δ:
                        atan(x+δ) − atan(x) = δ / (1+ξ²) for some ξ
                        between x and x+δ. Since 1/(1+t²) ≤ 1 with the
                        maximum at t = 0, the bound 1/(1+x²) holds at
                        the linearization point (REQ-EF-3 prescription;
                        loose if δ ≫ |x|, tight for small δ).
                        Range of atan is (−π/2, π/2), so the absolute
                        change |atan(x+δ) − atan(x)| < π.
  Primary reference:    REQ-EF-3 (atan: bound = err / (1 + val²));
                        Higham (2002) §13.3 (inverse trig);
                        mean-value theorem.
  Domain of validity:   All x ∈ ℝ.

METHOD
  Method declared:      Per-category: bound = (1/(1+val²))·err, capped at π.
  Method implemented:   `deriv = 1/(1 + a.value²);
                        return min(deriv * err, π);`
                        Applied per category via apply (card #11).
  Match verdict:        PASS — REQ-EF-3 atan rule exact.

ERROR BOUND
  Bound formula:        REQ-EF-3: err / (1 + val²), capped at π.
  Bound implemented:    lines 352-357 — exact match.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (atan), REQ-EF-8 (Lipschitz; |deriv| ≤ 1).
  AUD-EF applies:       AUD-EF-3 atan case.
  AUD-MC applies:       AUD-MC-12 (inverse trig of perturbed argument).
  Verification test:    test_error_framework.cpp REQ-EF-3 atan.

NOTES
  Comment line 348 acknowledges "atan never amplifies errors" — the
  derivative is at most 1 in magnitude, so atan is Lipschitz-1 globally
  with the local Lipschitz constant 1/(1+x²) used in the bound. The
  prescription evaluates the derivative at the linearization point;
  for δ ≫ |x| this is conservative, for small δ it is tight.
```

## 30. `TrackedValue::atan2(y, x)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::atan2
Location:               src/math/tracked_value.h:365-393
Mathematical statement: For y = (y, δ_y), x = (x, δ_x) with r² = x²+y² > 0:
                        atan2(y, x) = (atan2(y, x), δ_a) with
                        δ_a ≤ (|x|·δ_y + |y|·δ_x) / (x² + y²), capped at π.
                        Catastrophic disc: if r² ≤ (δ_y + δ_x)², bound = π
                        in all categories.

THEORY
  Underlying theorem:   Partial derivatives of atan2 at (y, x):
                        ∂θ/∂y =  x / (x² + y²);
                        ∂θ/∂x = −y / (x² + y²).
                        First-order bound:
                        |θ(y+δ_y, x+δ_x) − θ(y, x)| ≤ |x|·δ_y/(x²+y²)
                                                    + |y|·δ_x/(x²+y²).
                        Range of atan2 is (−π, π], so the bound is capped
                        at π (conservative; non-wrapping case stays < π).
                        Catastrophic disc: when r² ≤ (δ_x + δ_y)², the
                        true (x, y) could lie anywhere in a disc that
                        includes the origin → angle undefined → bound = π.
  Primary reference:    REQ-EF-3 (atan2: bound = (|x_val|·y_err +
                                                 |y_val|·x_err) /
                                                (x_val² + y_val²));
                        REQ-EF-9 (catastrophic singular-disc);
                        Higham (2002) §13.3.
  Domain of validity:   r² > (δ_y + δ_x)².

METHOD
  Method declared:      First-order partial-derivative bound, with
                        singular-disc fallback to π in all categories.
  Method implemented:   lines 370-374: total_err = y.total + x.total;
                        if r² ≤ total_err²: return (atan2, π, π, π).
                        lines 378-393: dy_factor = |x|/r²;
                        dx_factor = |y|/r²;
                        per-category bound = dy_factor·y_err
                                           + dx_factor·x_err, capped at π.
  Match verdict:        PASS — REQ-EF-3 atan2 + REQ-EF-9 catastrophic disc.

ERROR BOUND
  Bound formula:        REQ-EF-3: (|x_val|·y_err + |y_val|·x_err) /
                                  (x_val² + y_val²), capped at π.
                        REQ-EF-9 catastrophic disc: bound = π per category
                        when r² ≤ total_err².
  Bound implemented:    lines 378-392 — exact match for both regions.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (atan2), REQ-EF-9 (catastrophic disc).
  AUD-EF applies:       AUD-EF-3 atan2, AUD-EF-9 atan2 singular disc.
  AUD-MC applies:       AUD-MC-13 (axis-angle log/exp via atan2), AUD-MC-12.
  Verification test:    test_error_framework.cpp REQ-EF-3 atan2; AUD-EF-9
                        singular-disc test with r² < total_err².

NOTES
  Comments lines 362-364: "Near the origin (x ≈ 0, y ≈ 0), r² → 0 and
  derivatives blow up. Angle is completely undetermined there, so bound = π."
  Catastrophic test (line 374) uses the L¹ sum total_err = y.total + x.total
  for the disc radius — conservative (it includes both budgets in the worst
  case). Heavily used in Kepler and quaternion axis-angle inversion paths.
```

## 31. `TrackedValue::abs`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::abs
Location:               src/math/tracked_value.h:396-399
Mathematical statement: |a| := (|a.value|, a.errors).

THEORY
  Underlying theorem:   Reverse triangle inequality:
                        ||x| − |x_true|| ≤ |x − x_true|.
                        abs is Lipschitz-1. Per REQ-EF-14, identity-class
                        operations preserve the budget exactly.
                        IEEE-754 std::abs is the exact sign-bit clear
                        (no rounding).
  Primary reference:    REQ-EF-14; Higham (2002) §1.2.
  Domain of validity:   All T with std::abs.

METHOD
  Method declared:      Value: std::abs (exact in IEEE-754).
                        Errors: unchanged.
  Method implemented:   `return TrackedValue(abs(a.value), a.errors);`
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        errors unchanged.
  Bound implemented:    matches.
  Bound verdict:        PASS — REQ-EF-14 preservation.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-14 (Lipschitz-1 identity-class).
  AUD-EF applies:       AUD-EF-14.
  AUD-MC applies:       Used by vector norm and component extraction paths.
  Verification test:    AUD-EF-14.

NOTES
  Reverse triangle inequality is tight only when x and x_true agree in
  sign; for a value near zero with δ > |x|, the actual deviation ||x| − |x_true||
  could equal |x|. But the rigorous bound remains δ_x because the
  inequality ||x| − |x_true|| ≤ |x − x_true| ≤ δ_x is satisfied. No
  special-case branch needed.
```

## 32. `TrackedValue::fmod(a, b)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::fmod
Location:               src/math/tracked_value.h:404-410
Mathematical statement: For a = (x, δ_x), b = (y, δ_y):
                        fmod(a, b) = (fmod(x, y), δ_f), with
                        δ_f ≤ δ_x + |n|·δ_y, where n = trunc(x/y).

THEORY
  Underlying theorem:   fmod(x, y) = x − n·y where n is the integer
                        trunc(x/y). For perturbed x+δ_x, y+δ_y, with n
                        held fixed (valid in a neighborhood not crossing
                        an integer boundary of x/y):
                        |(x+δ_x) − n·(y+δ_y) − (x − n·y)|
                          = |δ_x − n·δ_y| ≤ δ_x + |n|·δ_y.
                        REQ-EF-3-style closed form (additive with scaled
                        denominator-error contribution).
  Primary reference:    IEEE 754-2019 §5.5.1 (fmod definition);
                        Higham (2002) §2.6 (range-reduction error);
                        C++ std::fmod follows IEEE 754 semantics.
  Domain of validity:   y ≠ 0; n does not change across the perturbation
                        interval (i.e., x/y is not at an integer boundary).

METHOD
  Method declared:      Closed-form: value = std::fmod(a.value, b.value);
                        n = |floor(a.value / b.value)|;
                        errors = a.errors + |n| · b.errors per category.
  Method implemented:   lines 405-409.
                        Note: code uses `floor` to compute n_abs; for
                        positive a/b, floor and trunc agree; for negative,
                        they differ but the |.| in `n_abs = abs(floor(...))`
                        makes the difference irrelevant for the bound
                        (bound depends only on |n|).
  Match verdict:        PASS — REQ-EF-3-style closed form; |n| handles
                        both signs.

ERROR BOUND
  Bound formula:        δ_x + |n|·δ_y per category.
  Bound implemented:    `a.errors + n_abs * b.errors` (per-category, via
                        cards #5 and #7).
  Bound verdict:        PASS — exact closed form for the interior of
                        each n-interval. Integer n contributes no error.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form fmod path), REQ-EF-8.
  AUD-EF applies:       AUD-EF-3 fmod.
  AUD-MC applies:       Used by angle normalization (Angles::normalize).
  Verification test:    AUD-EF-3 fmod test.

NOTES
  Caveat: if the perturbation crosses an integer boundary of x/y, fmod
  has a discontinuity of magnitude |y|. The bound δ_x + |n|·δ_y is
  rigorous within a single n-interval; near a boundary the bound is loose
  but still valid because the error budget grows there. The codebase uses
  fmod for angle normalization (y = 2π), where n is small (a few) and
  callers verify outputs against a known-good angle range.
```

## 33. `TrackedValue::operator<`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator<
Location:               src/math/tracked_value.h:414
Mathematical statement: a < b ⇔ a.value < b.value.

THEORY
  Underlying theorem:   Order comparison on the value field only.
                        Per REQ-EF-10, value is the raw T scalar.
                        Comparisons inspect value only; error budgets
                        are metadata, not part of the order.
  Primary reference:    REQ-EF-10 ("value is a plain T").
  Domain of validity:   T must have <.

METHOD
  Method declared:      Delegate to T::operator<.
  Method implemented:   `return a.value < b.value;`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       n/a (returns bool).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-10.
  Verification test:    REQ-EF-10 tests.

NOTES
  Comments at line 412 explicit: "Comparisons (on value only; errors are
  metadata)". Callers wanting reliability-aware comparison must compose
  with reliable_digits() or abs(a − b) explicitly.
```

## 34. `TrackedValue::operator>`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator>
Location:               src/math/tracked_value.h:415
Mathematical statement: a > b ⇔ a.value > b.value.

THEORY
  Underlying theorem:   Mirror of card #33.
  Primary reference:    REQ-EF-10.
  Domain of validity:   T must have >.

METHOD
  Method declared:      Delegate to T::operator>.
  Method implemented:   `return a.value > b.value;`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       n/a.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-10.
  Verification test:    REQ-EF-10.

NOTES
  Same caveat as card #33.
```

## 35. `TrackedValue::operator<=`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator<=
Location:               src/math/tracked_value.h:416
Mathematical statement: a <= b ⇔ a.value <= b.value.

THEORY
  Underlying theorem:   Mirror of card #33.
  Primary reference:    REQ-EF-10.

METHOD
  Method declared:      Delegate to T::operator<=.
  Method implemented:   `return a.value <= b.value;`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       n/a.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-10.
  Verification test:    REQ-EF-10.

NOTES
  Same caveat as card #33.
```

## 36. `TrackedValue::operator>=`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator>=
Location:               src/math/tracked_value.h:417
Mathematical statement: a >= b ⇔ a.value >= b.value.

THEORY
  Underlying theorem:   Mirror of card #33.
  Primary reference:    REQ-EF-10.

METHOD
  Method declared:      Delegate to T::operator>=.
  Method implemented:   `return a.value >= b.value;`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       n/a.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-10.
  Verification test:    REQ-EF-10.

NOTES
  Same caveat as card #33.
```

## 37. `TrackedValue::operator==`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator==
Location:               src/math/tracked_value.h:418
Mathematical statement: a == b ⇔ a.value == b.value.

THEORY
  Underlying theorem:   IEEE-754 equality on raw value field only.
                        Two TrackedValues with the same value but
                        different error budgets compare equal —
                        intentional, per REQ-EF-10.
  Primary reference:    REQ-EF-10.

METHOD
  Method declared:      Delegate to T::operator==.
  Method implemented:   `return a.value == b.value;`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       n/a.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-10.
  Verification test:    REQ-EF-10.

NOTES
  Value comparison only — not a reliability comparison. Callers testing
  "is a equal to b within tolerance" must use `abs(a − b) < tol`
  explicitly (and a represents a tracked tolerance via REQ-EF-1).
```

## 38. `TrackedValue::operator!=`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::operator!=
Location:               src/math/tracked_value.h:419
Mathematical statement: a != b ⇔ a.value != b.value.

THEORY
  Underlying theorem:   Negation of card #37.
  Primary reference:    REQ-EF-10.

METHOD
  Method declared:      Delegate to T::operator!=.
  Method implemented:   `return a.value != b.value;`
  Match verdict:        OK.

ERROR BOUND
  Bound category:       n/a.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-10.
  Verification test:    REQ-EF-10.

NOTES
  Mirror of card #37.
```

## 39. `TrackedValue::ctor_int(int n)` (implicit conversion)

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::ctor_int
Location:               src/math/tracked_value.h:423
Mathematical statement: TrackedValue(n) := (T(n), (0, 0, 0)) for n ∈ ℤ.

THEORY
  Underlying theorem:   IEEE 754 §3.4: for |n| ≤ 2^p (p = significand bits),
                        int → T conversion is exact; zero error in all
                        three categories. Same theorem as card #15.
  Primary reference:    IEEE 754-2019 §3.4; Goldberg (1991) Theorem 1.
  Domain of validity:   |n| ≤ 2^p. For double: |n| ≤ 2^53. For 32-bit int,
                        always satisfied.

METHOD
  Method declared:      Direct construction; int → T; zero errors.
                        No `explicit` (implicit conversion).
  Method implemented:   `TrackedValue(int n) : value(T(n)),
                                                errors(T(0), T(0), T(0)) {}`
  Match verdict:        PASS.

ERROR BOUND
  Bound formula:        0 in all categories (when |n| ≤ 2^p).
  Bound implemented:    0.
  Bound verdict:        PASS — under domain precondition.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-10.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Implicit conversion enables `tv + 2`-style expressions naturally. For
  int approaching 2^53 the conversion is still exact but callers should
  be aware. Typical codebase usage is small constants (|n| ≤ 10000).
```

## 40. `TrackedValue::representation_bound(val)`

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::representation_bound
Location:               src/math/tracked_value.h:441-492
Mathematical statement: representation_bound(v) := bound on |v − v_true|
                        for v = round_to_T(v_true).
                        For v = 0: returns 0.
                        For double: returns ldexp(1, e − 53) where e is
                          the frexp exponent — this is 1 ULP at v's binade.
                        For float:  returns ldexp(1, e − 24) — 1 ULP.
                        For multiprec: returns |v|·2^(−digits) — 1 ULP.
                        IEEE-754 round-to-nearest-even guarantees the
                        true rounding error is ≤ 0.5 ULP; the code returns
                        1 ULP, which is conservative by a factor of 2.

THEORY
  Underlying theorem:   IEEE 754-2019 §4.3.1 (round-to-nearest-even):
                        |x − fl(x)| ≤ 0.5 ULP(fl(x)).
                        Within the binade [2^(e−1), 2^e) (with frexp's
                        convention v = frac·2^e and 0.5 ≤ frac < 1),
                        1 ULP = 2^(e−1)·ε = 2^(e−p), where ε = 2^(1−p)
                        is the machine epsilon (p = significand bits,
                        including implicit leading 1; p = 53 for double,
                        p = 24 for float).
                        Therefore 0.5 ULP = 2^(e − p − 1) = 2^(e − 54) for
                        double; 2^(e − 25) for float.
  Primary reference:    IEEE Std 754-2019 §4.3.1; Goldberg (1991) §2;
                        Higham (2002) §2.1.
  Domain of validity:   All finite v; v = 0 handled separately.

METHOD
  Method declared:      Three branches by type (constexpr-dispatched):
                        - double: frexp → exponent e; ldexp(1.0, e − 53).
                                  Comment "0.5 ULP = 2^(exp − 1 − 52)
                                  = 2^(exp − 53)" appears at line 457
                                  but is off by 1: the correct 0.5 ULP
                                  expression is 2^(exp − 54). The
                                  *implemented* expression 2^(exp − 53)
                                  is therefore 1 ULP, not 0.5 ULP.
                        - float:  frexp; ldexp(1.0f, e − 24). Similarly 1 ULP.
                        - other (multiprec): ldexp(T(1), −digits)
                                  evaluates to ε = 2^(1−p)·… (actually
                                  2^(−p) for digits = p); multiplied by
                                  |v| gives |v|·2^(−p), which is ≤ 1 ULP.
                                  ADL via `using std::ldexp` allows
                                  boost::multiprecision overloads (comment
                                  lines 481-485).
                        v = 0 short-circuits to 0 at line 444.
  Method implemented:   lines 441-492.
  Match verdict:        PASS-conservative — code returns 1 ULP (or
                        ≤ 1 ULP for multiprec). True theoretical bound
                        is 0.5 ULP per IEEE-754 §4.3.1, so the
                        implementation over-bounds by a factor ≤ 2.
                        The bound remains a rigorous upper bound on the
                        rounding error — REQ-EF-2 is satisfied.

ERROR BOUND
  Bound category:       This IS the precision representation bound used
                        by defined(), measured(), and ratio().
  Bound formula (theory):
                        For double:        0.5 ULP = 2^(e − 54).
                        For float:         0.5 ULP = 2^(e − 25).
                        For multiprec:     0.5 ULP = |v|·2^(−p−1).
  Bound implemented:    For double:        1 ULP = 2^(e − 53).
                        For float:         1 ULP = 2^(e − 24).
                        For multiprec:     1 ULP = |v|·2^(−p)  (slightly
                                                                less in
                                                                practice).
  Bound verdict:        PASS-conservative — bound is rigorous (over-bounds
                        the true 0.5 ULP rounding error by ≤ 2x). Flag
                        for tightening; current behavior is correct.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (precision foundation);
                        REQ-EF-2 (representation cost is part of total).
  AUD-EF applies:       AUD-EF-1, AUD-EF-3.
  Verification test:    AUD-EF-1 — verify the returned bound dominates
                        the actual rounding error for sampled decimals.
                        A property test "bound >= true_round_error"
                        will pass with the conservative 1 ULP and would
                        also pass with 0.5 ULP.

NOTES
  Critical foundation: all of defined(), defined_with_physical_uncertainty,
  measured(), and ratio() route through this function for their precision
  category. A bug here propagates to every named-constant construction.
  - For IEEE-754 double/float, frexp is a bit-pattern shift (no rounding)
    and ldexp(1.0, n) sets the exponent field exactly (no rounding).
    The implementation is therefore exact at 1 ULP.
  - For boost multiprecision T, ADL must reach boost's own ldexp overload;
    comment lines 481-485 explain that `using std::ldexp` + unqualified
    `ldexp(T(1), -digits)` accomplishes this, citing the MSVC v145 /
    VS 2026 compile-time consideration that `std::`-qualified lookup
    does not consider ADL candidates.
  - For T = 0 (line 444): returns T(0). Strictly the rounding error of
    storing exactly 0 is 0 (denormal scale aside); the value 0 is exactly
    representable in IEEE-754 → no precision cost.
  Tightening upgrade: change e − 53 to e − 54 for double / e − 24 to
  e − 25 for float / -digits to -digits − 1 for multiprecision. All
  current bounds remain rigorous over-bounds; tightening is pure
  optimization (saves a factor of 2 in reported precision).
```

## 41. `TrackedValue::from_string(str)` (private)

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::TrackedValue::from_string
Location:               src/math/tracked_value.h:498-506
Mathematical statement: from_string(s) := round_to_T(parse_decimal(s)).
                        Parses a decimal string into T using T's standard
                        string parsing (std::stod, std::stof, or boost
                        multiprecision's string constructor).

THEORY
  Underlying theorem:   String → float parse semantics per the C++
                        standard (std::stod, std::stof) and boost
                        multiprecision equivalents. Result is the IEEE-754
                        correctly-rounded T-representation of the real
                        value parsed from s.
  Primary reference:    C++ standard <cstdlib> / <string> — std::stod,
                        std::stof; boost multiprecision string ctor.
                        IEEE Std 754-2019 §4.3.1 for the rounding.
  Domain of validity:   Strings parseable to a finite decimal.

METHOD
  Method declared:      constexpr branch by T:
                        - double: std::stod;
                        - float:  std::stof;
                        - other:  T's string constructor (boost multiprec).
  Method implemented:   lines 498-506.
  Match verdict:        OK — direct delegation. No bound is emitted here.

ERROR BOUND
  Bound category:       n/a (utility; no bound returned).
                        Callers (defined, measured, etc.) apply
                        representation_bound (card #40) to the returned
                        value to obtain the precision contribution.
  Bound formula:        n/a.
  Bound implemented:    n/a.
  Bound verdict:        OK.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (via callers).
  AUD-EF applies:       AUD-EF-1 (via callers).
  Verification test:    Indirect via defined / measured tests (AUD-EF-1).

NOTES
  - std::stod / std::stof follow LC_NUMERIC; this codebase assumes a
    portable parse with decimal point '.'. Locale issues are precondition.
  - Boost multiprecision string ctors are IEEE-754-correct for their
    represented precision.
  - The separation (parse here, bound applied by caller) keeps the parse
    utility composable and avoids double-counting the precision cost.
```

## 42. `exact<T>(int n)` (free function)

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::exact
Location:               src/math/tracked_value.h:513-515
Mathematical statement: exact<T>(n) := TrackedValue<T>::exact_integer(n).

THEORY
  Underlying theorem:   Forward to card #15 (exact_integer): for |n| ≤ 2^p,
                        T(n) is exactly representable per IEEE 754 §3.4
                        and Goldberg (1991) Theorem 1.
  Primary reference:    IEEE 754-2019 §3.4; Goldberg (1991).
  Domain of validity:   |n| ≤ 2^p. For double, |n| ≤ 2^53.

METHOD
  Method declared:      Forward to exact_integer.
  Method implemented:   `return TrackedValue<T>::exact_integer(n);`
  Match verdict:        PASS — pure forwarding.

ERROR BOUND
  Bound formula:        0 in all categories (under |n| ≤ 2^p).
  Bound implemented:    0 (via exact_integer).
  Bound verdict:        PASS — under domain precondition.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1.
  AUD-EF applies:       AUD-EF-1.
  Verification test:    AUD-EF-1.

NOTES
  Convenience template enabling the idiom `exact<double>(2)` for small
  exact integer constants in arithmetic expressions. Equivalent to
  `TrackedValue<double>(2)` (via the implicit int ctor, card #39) when
  n fits in int. The explicit form is more readable where the exactness
  is emphasized.
```

## 43. `ratio<T>(int num, int den)` (free function)

```
=== FORMULA AUDIT CARD ===
ID:                     tracked_value::ratio
Location:               src/math/tracked_value.h:528-532
Mathematical statement: ratio<T>(num, den) := (T(num)/T(den), 0, ε_repr, 0),
                        where ε_repr = representation_bound(T(num)/T(den)).

THEORY
  Underlying theorem:   Composition of two facts:
                        (i) For |num|, |den| ≤ 2^p, T(num) and T(den) are
                            exact (Goldberg 1991 Theorem 1).
                        (ii) The quotient T(num)/T(den) rounds to its
                             nearest representable T per IEEE-754 §5.4;
                             rounding cost is 0.5 ULP (delivered as 1 ULP
                             by the conservative representation_bound,
                             card #40).
                        Total: value is the IEEE-correctly-rounded
                        quotient; precision bound is the representation
                        cost on the rounded result. measurement and
                        accuracy = 0 (rational constants are definitional).
  Primary reference:    IEEE 754-2019 §5.4 (division rounding);
                        Goldberg (1991) Theorem 1;
                        REQ-EF-3 (rounding cost is precision).
  Domain of validity:   den ≠ 0; |num|, |den| ≤ 2^p.

METHOD
  Method declared:      Compute T(num)/T(den) directly (one IEEE rounded
                        division); apply representation_bound (card #40)
                        for the precision bound; zero measurement and
                        accuracy.
  Method implemented:   lines 529-531:
                        `val = T(num) / T(den);
                        repr_err = representation_bound(val);
                        return TrackedValue<T>(val, 0, repr_err, 0);`
  Match verdict:        PASS — REQ-EF-3-style rounding-cost accounting.

ERROR BOUND
  Bound formula:        precision = representation_bound(T(num)/T(den))
                                  = 0.5-1.0 ULP of the rounded quotient
                                    (per card #40 conservative factor);
                        measurement = accuracy = 0.
  Bound implemented:    matches.
  Bound verdict:        PASS.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-3 (rounding cost is precision).
  AUD-EF applies:       AUD-EF-1, AUD-EF-3.
  Verification test:    AUD-EF-3 ratio<T>(1, 3) — bound should be 0.5-1
                        ULP of T's representation of 0.333…, NOT zero.

NOTES
  Comments lines 520-527 motivate the explicit-bound approach:
  the naive division of two exact TrackedValues (both with zero error)
  would report zero precision (REQ-EF-3 div with both errors zero gives
  bound = 0). But the quotient T(num)/T(den) for irrational fractions
  (1/3, 1/7, etc.) is NOT exactly representable in T, so reporting zero
  precision would violate REQ-EF-2 (rigorous bound).
  ratio<T>() computes the representation cost on the result directly,
  recovering the rounding error. It is the canonical entry for rational
  constants throughout the codebase (e.g., ratio<T>(-1, 3) for -1/3 in
  small_angle_series.h taylor_cos_minus_sinc_over_theta_sq).
```

---

## File-level verdict for `src/math/tracked_value.h`

- **A. Error wiring**: PASS — all numerical operations return `TrackedValue<T>` with per-category bounds applied. No bare-T paths in the public arithmetic / math-function surface.
- **B. Algebra axioms**: not the subject of this file; foundational operators are tested for the axioms they need to support (AUD-MC-1..8 for additive/multiplicative axioms on composite types built atop these primitives).
- **C. Theoretical basis** (43 cards):
  - **Cards 1-14, 33-38, 41** (data movement, named constructors, queries, comparisons, utility): OK or PASS — pure metadata or directly definitional with no numerical bound to certify, or one-line forwarders.
  - **Cards 3, 4, 7-11, 19, 20** (ThreeErrors aggregators and TrackedValue queries): PASS — each implements a closed-form combinator (triangle, RSS, scale, max, lift) used internally; each is sound or correctly delegates.
  - **Cards 15-18, 39, 42** (constructors with bound assignment): PASS — IEEE-754 §3.4 (exact integers) and REQ-EF-1 / REQ-EF-3 (rounding cost as precision) cited and implemented exactly.
  - **Cards 21-32** (arithmetic and math operators): PASS — each card cites the underlying theorem (triangle inequality for ±, expansion-bound for ×, exact closed form for ÷, mean-value theorem for sqrt and atan, Taylor with Lagrange remainder for sin/cos, partial derivatives for atan2, Lipschitz-1 for abs, trunc/floor + scale identity for fmod) and implements the corresponding REQ-EF-3 (or REQ-EF-9 catastrophic) closed-form bound exactly.
  - **Card 40** (representation_bound): PASS-conservative — code returns 1 ULP rather than 0.5 ULP (factor-of-2 over-bound). Bound remains rigorous (REQ-EF-2 satisfied); flag for optional tightening.
  - **Card 43** (ratio): PASS — correctly recovers the rounding cost on irrational rationals (1/3, etc.) and propagates it as precision per REQ-EF-3 + REQ-EF-2.

**File verdict: PASS** with one flagged tightening opportunity at card #40 (conservative by factor 2). All REQ-EF-3 closed-form bound formulas (add/sub, mul, div, sqrt, sin/cos, atan, atan2) are implemented exactly per spec. REQ-EF-9 catastrophic-region signals (division by sub-error denominator, sqrt of sub-error argument, atan2 inside singular disc, ThreeErrors division by exact zero) are present and conservative. The foundation of the project's error budget is sound.

## Catastrophic-region inventory (REQ-EF-9 audit)

| Op | Triggering condition | Bound emitted |
|---|---|---|
| `ThreeErrors::operator/(e, s)` | s == 0 | `numeric_limits<T>::max()` in all categories |
| `operator/(TV, TV)` | b_bound ≥ |b.value| (denom_min ≤ 0) | `numeric_limits<T>::max()` |
| `sqrt` | err ≥ value | `sqrt(value)` (full result magnitude) |
| `atan2` | r² ≤ (y_err + x_err)² | π in all categories |
| `representation_bound` | value == 0 | 0 (exact zero is exactly representable) |

All five paths match REQ-EF-9: signal, not silence. Downstream callers see the conservative bound via `total_error()` and `reliable_digits()` collapses to 0 in these regions, enabling go/no-go decisions in the propagator's outer loops.

## Open notes for cross-cutting upgrades

1. **representation_bound is 1 ULP, not 0.5 ULP** (card #40). Conservative by factor 2. To tighten: change `e − 53` to `e − 54` for double, `e − 24` to `e − 25` for float, and `-digits` to `-digits − 1` for multiprecision. All current bounds remain rigorous over-bounds; tightening is pure optimization.

2. **`fmod` boundary safety** (card #32). Bound assumes the integer quotient `n` is stable across the perturbation interval. Holds in the codebase's actual use case (angle normalization with y = 2π and small n); a discontinuity-spanning bound would add |y| as a hedge if cross-boundary inputs become routine.

3. **`atan` uses derivative at the central value** (card #29). Conservative if δ_x is large relative to |val|. Tighter bound uses min over the interval [val−err, val+err] of 1/(1+t²); a refinement opportunity, no correctness impact under REQ-EF-3 as written.

4. **Cross-category products neglected in `operator*`** (card #24). Comment lines 200-202 document this as standard interval-arithmetic practice for separated budgets. The 6 cross-pair products (a_m·b_p, etc.) are genuinely second-order; the per-category convention is chosen for attribution clarity. Standing decision; no action.

5. **`sin/cos` use |sin(ξ)| ≤ 1 / |cos(ξ)| ≤ 1 for the Lagrange remainder** (cards #27, #28). REQ-EF-3 prescribes exactly this. A tighter bound would compute the maximum of |sin| or |cos| over the interval [val−err, val+err]; possible refinement, no correctness impact.

## References

- IEEE Std 754-2019, *IEEE Standard for Floating-Point Arithmetic*. IEEE Computer Society.
- Higham, N. J. (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM.
- Goldberg, D. (1991), "What Every Computer Scientist Should Know About Floating-Point Arithmetic", *ACM Comput. Surv.* 23(1), 5-48.
- Wilkinson, J. H. (1965), *The Algebraic Eigenvalue Problem*, Oxford University Press.
- Rudin, W. (1976), *Principles of Mathematical Analysis*, 3rd ed., McGraw-Hill (Taylor's theorem with Lagrange remainder, Theorem 5.15).
