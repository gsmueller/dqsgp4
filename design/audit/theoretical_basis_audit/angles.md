# Theoretical Basis Audit — `src/math/angles.h`

**Scope.** Six angle-manipulation entry points: `pi()`, `two_pi()`, `degrees_to_radians()`, `radians_to_degrees()`, `wrap_two_pi()`, `wrap_neg_pos_pi()`. All inputs/outputs are `TrackedValue<T>` in radians.

Companion to framework `design/audit/theoretical_basis_audit.md` (§1 card schema, §5 worked example).

---

## Card 1 — `angles::pi`

```
=== FORMULA AUDIT CARD ===
ID:                     angles::pi
Location:               src/math/angles.h:13-17
Mathematical statement: Return the constant π wrapped as a TrackedValue<T>,
                        with measurement = 0, precision = π·ε, accuracy = 0.

THEORY
  Underlying theorem:   Definition. π is the mathematical constant; the
                        Boost.Math constants library returns it to the
                        maximum representable precision of T (boost::math
                        ::constants::pi<T>() returns the rounded T-image
                        of π with relative error ≤ ε(T)/2).
  Primary reference:    Boost.Math Constants documentation — every floating-
                        point T has |fl(π) − π| ≤ |π|·ε(T)/2 by
                        IEEE round-to-nearest of an exactly-representable
                        constant table. (Higham, "Accuracy and Stability of
                        Numerical Algorithms" 2nd ed., §2.2)
  Domain of validity:   All T satisfying the boost::math::constants concept.

METHOD
  Method declared:      Closed-form table-lookup constant via
                        boost::math::constants::pi<T>().
  Method implemented:   `T val = boost::math::constants::pi<T>();`
                        wrapped as TrackedValue(val, 0, val·ε, 0).
  Match verdict:        ✓ matched — closed-form constant, no truncation, no
                        iteration. The only error source is the floating-
                        point representation of π itself.

ERROR BOUND
  Bound category:       precision
  Bound formula:        |fl(π) − π| ≤ π · ε(T) (a conservative upper bound;
                        the half-ulp bound π·ε/2 also valid but rarely
                        materially tighter at this scale).
  Bound implemented:    `val * std::numeric_limits<T>::epsilon()` i.e. π·ε.
  Bound verdict:        ✓ matched — bound is the representation bound of π
                        in type T, conservative by factor 2 over half-ulp.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form constants carry their
                        representation bound as precision).
  AUD-EF applies:       AUD-EF-1 (every public op returns TrackedValue<T>).
  AUD-MC applies:       n/a (constant, not an algebra op).
  Verification test:    None dedicated; exercised transitively via
                        tests/test_math (deg→rad conversions).

NOTES
  - The bound `val·ε` is the conservative form of the half-ulp bound; both
    are rigorous. Half-ulp `val·ε/2` would be tighter but the factor-2
    margin absorbs any rounding mode subtlety.
  - For wide T (e.g. cpp_bin_float_50) boost::math::constants returns π
    to the full type precision; the same `val·ε` form remains rigorous.
```

---

## Card 2 — `angles::two_pi`

```
=== FORMULA AUDIT CARD ===
ID:                     angles::two_pi
Location:               src/math/angles.h:20-23
Mathematical statement: Return 2π as a TrackedValue<T>.

THEORY
  Underlying theorem:   Definition + closed-form scalar multiplication.
                        2π = 2 · π; multiplication of an exact integer by a
                        TrackedValue propagates precision exactly per
                        REQ-EF-3 (the integer 2 is representable with zero
                        error in any IEEE T; the product 2·π may incur one
                        rounding ≤ ε but the precision-propagation rule
                        |2|·δπ = 2·π·ε dominates).
  Primary reference:    REQ-EF-3 (closed-form error propagation, see
                        design/specifications/error_framework.md).
                        Higham §2.2 for floating-point multiplication of an
                        exactly-representable integer by a near-exact value.
  Domain of validity:   All T.

METHOD
  Method declared:      Closed-form: exact<T>(2) * pi<T>().
  Method implemented:   `return exact<T>(2) * pi<T>();`
  Match verdict:        ✓ matched — closed-form scalar product, no
                        truncation, no iteration.

ERROR BOUND
  Bound category:       precision
  Bound formula:        δ(2π) = 2 · δπ = 2 · π · ε = 2π · ε.
                        (TrackedValue::operator* on exact_integer(2) gives
                        precision = 2·precision(π), see REQ-EF-3.)
  Bound implemented:    Inherited from `operator*` on `TrackedValue<T>`.
                        Result precision = 2 · (π·ε) = 2π·ε.
  Bound verdict:        ✓ matched — bound is the representation bound of
                        2π scaled by the precision-multiplication rule.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3.
  AUD-EF applies:       AUD-EF-1, AUD-EF-2 (multiplication error wiring).
  AUD-MC applies:       n/a.
  Verification test:    Exercised by `wrap_two_pi` callers in
                        tests/test_math.

NOTES
  - Result precision matches the half-ulp of 2π under double precision
    (2π · ε ≈ 1.4e-15), as expected for a near-exactly-representable
    multiple of π.
  - No additional rounding error from the multiplication is tracked
    explicitly; this is absorbed into the conservative `val·ε`
    representation bound already on π.
```

---

## Card 3 — `angles::degrees_to_radians`

```
=== FORMULA AUDIT CARD ===
ID:                     angles::degrees_to_radians
Location:               src/math/angles.h:26-29
Mathematical statement: rad = deg · π / 180.

THEORY
  Underlying theorem:   Closed-form unit conversion identity. The conversion
                        factor π/180 is exact in the reals; in floating-
                        point T it inherits the precision of π plus the
                        rounding of the division by the exact integer 180.
  Primary reference:    SI/ISO 80000-3 unit-conversion definition;
                        REQ-EF-3 (closed-form propagation: for
                        z = x · y / k with k an exact integer,
                        δz = (|y/k|·δx + |x/k|·δy) by triangle inequality).
  Domain of validity:   All real T; no branch cuts.

METHOD
  Method declared:      Closed-form scalar product/quotient
                        `degrees * pi<T>() / exact<T>(180)`.
  Method implemented:   `return degrees * pi<T>() / exact<T>(180);`
  Match verdict:        ✓ matched — single closed-form expression, no
                        Taylor, Padé, or iteration.

ERROR BOUND
  Bound category:       precision (representation) + propagated measurement.
  Bound formula:        δ_precision = (π/180)·δdeg_precision
                                    + (deg/180)·δπ_precision
                        δ_measurement = (π/180)·δdeg_measurement.
                        Both follow from REQ-EF-3 closed-form propagation
                        rules for `operator*` and `operator/` chained.
  Bound implemented:    Inherited from `TrackedValue<T>::operator*` and
                        `TrackedValue<T>::operator/` (see tracked_value.h
                        AUD-EF-2,3).
  Bound verdict:        ✓ matched — bound is exactly what closed-form
                        operator chaining produces under REQ-EF-3.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3.
  AUD-EF applies:       AUD-EF-1, AUD-EF-2 (operator wiring).
  AUD-MC applies:       n/a (unit conversion, not algebra-axiom).
  Verification test:    tests/test_math/main.cpp:52 — deg90 = 90°→π/2.

NOTES
  - The division by `exact<T>(180)` does not add measurement error
    (180 is integer-exact); it adds the representation bound of 1/180
    via `ratio`-style propagation already wired into operator/.
  - At T = double, output precision on deg=90 input with δdeg=1e-15 is
    π/2 · 1e-15 plus π·ε/180 ≈ 2·10⁻¹⁷ — dominated by input δdeg.
```

---

## Card 4 — `angles::radians_to_degrees`

```
=== FORMULA AUDIT CARD ===
ID:                     angles::radians_to_degrees
Location:               src/math/angles.h:32-35
Mathematical statement: deg = rad · 180 / π.

THEORY
  Underlying theorem:   Inverse of degrees_to_radians; same closed-form
                        identity. Multiplication by exact integer 180 and
                        division by π propagate precision per REQ-EF-3.
                        Division by π is well-defined since π ≠ 0 with
                        |π − fl(π)| ≪ |π|, so the propagation
                        δ(x/π) = (1/π)δx + (|x|/π²)δπ remains finite and
                        rigorous (mean-value theorem on 1/y near y=π).
  Primary reference:    Same as Card 3. Higham §2.4 for division error
                        propagation.
  Domain of validity:   All real T.

METHOD
  Method declared:      Closed-form `radians * exact<T>(180) / pi<T>()`.
  Method implemented:   `return radians * exact<T>(180) / pi<T>();`
  Match verdict:        ✓ matched — closed-form, no iteration or series.

ERROR BOUND
  Bound category:       precision + measurement (propagated).
  Bound formula:        δ_precision = (180/π)·δrad_precision
                                    + (180·rad/π²)·δπ_precision.
                        δ_measurement = (180/π)·δrad_measurement.
                        (Follows from REQ-EF-3 quotient rule applied to
                        `operator/` in tracked_value.h.)
  Bound implemented:    Inherited from `operator*` then `operator/`.
  Bound verdict:        ✓ matched — bound matches REQ-EF-3 closed-form
                        propagation for the chained product/quotient.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3.
  AUD-EF applies:       AUD-EF-1, AUD-EF-2.
  AUD-MC applies:       n/a.
  Verification test:    Round-trip with degrees_to_radians; not yet a
                        dedicated test in tests/test_math (recommend adding
                        a deg→rad→deg round-trip identity check).

NOTES
  - The amplifier 180/π ≈ 57.3 means input precision in radians is
    multiplied by ~57 in the degrees output — a known property of this
    conversion, not a bound failure.
```

---

## Card 5 — `angles::wrap_two_pi`

```
=== FORMULA AUDIT CARD ===
ID:                     angles::wrap_two_pi
Location:               src/math/angles.h:38-47
Mathematical statement: Wrap angle θ to the half-open interval [0, 2π);
                        result = θ mod 2π, with negative results shifted
                        by +2π.

THEORY
  Underlying theorem:   The mod identity: for any real θ, the value
                        θ − 2π·⌊θ/(2π)⌋ ∈ [0, 2π). Equivalently
                        fmod(θ, 2π) ∈ (−2π, 2π) per C99/IEEE-754, with a
                        single +2π correction when negative giving the
                        canonical principal value in [0, 2π).
                        fmod itself is Lipschitz-1 in its first argument
                        (the file header comment) and piecewise affine,
                        so per-category error in θ passes through as-is;
                        error in 2π scales by |⌊θ/(2π)⌋|.
  Primary reference:    Goldberg (1991) "What Every Computer Scientist
                        Should Know About Floating-Point Arithmetic",
                        Theorem 11 (fmod/remainder is exact when the
                        quotient is representable). C99 §7.12.10.1 fmod.
                        ProofWiki — Real Modulo Operation Properties.
  Domain of validity:   All finite real θ. For θ that is an exact integer
                        multiple of 2π, fmod returns ±0 (handled by the
                        sign-correction branch).

METHOD
  Method declared:      Closed-form mod-reduction via `fmod(θ, 2π)` plus
                        single conditional add of 2π if negative.
  Method implemented:   `auto result = fmod(angle, tp);`
                        `if (result.value < T(0)) result = result + tp;`
                        — fmod is the `tracked_value.h:404` overload, which
                        already encodes the bound
                        δ(fmod) ≤ δθ + |⌊θ/2π⌋|·δ(2π).
  Match verdict:        ✓ matched — closed-form mod, no Taylor, no Padé,
                        no iteration. The cited theory (mod identity +
                        Lipschitz-1 fmod) matches the implementation
                        exactly.

ERROR BOUND
  Bound category:       precision + measurement (passed-through).
  Bound formula:        Per tracked_value.h fmod overload:
                          δ(fmod(θ, 2π)) ≤ δθ + |⌊θ/(2π)⌋|·δ(2π)
                        for each error category. The +2π correction in the
                        negative branch reuses `operator+`, adding another
                        δ(2π) to the precision (only when triggered).
  Bound implemented:    Inherited from `fmod` (tracked_value.h:404-410):
                          `TrackedValue(result, a.errors + n_abs * b.errors)`
                        with n_abs = |⌊θ/(2π)⌋|. The conditional
                        `result + tp` adds δ(2π) when the sign correction
                        fires (REQ-EF-3 addition rule).
  Bound verdict:        ✓ matched — bound is the rigorous fmod bound
                        cited in tracked_value.h:401-402, exactly equal
                        to δθ + n·δ(2π) for each category. The bound
                        DEPENDS on the tracked_value::fmod bound (per the
                        task's note); that bound is the Lipschitz-1
                        composition of the mod identity, so the chain is
                        rigorous.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation under fmod).
  AUD-EF applies:       AUD-EF-2 (binary-op wiring — fmod and operator+
                        both add to errors per their category).
  AUD-MC applies:       n/a (utility, not algebra-axiom op).
  Verification test:    No dedicated wrap test currently in
                        tests/test_math; recommend adding one for
                        θ ∈ {0, π, 2π, 3π, −π, very-large multiples}.

NOTES
  - For very large |θ| (≫ 2π), |⌊θ/(2π)⌋| can amplify δ(2π) substantially.
    This is Goldberg's classical argument-reduction issue: the precision
    in 2π is amplified by the integer quotient before subtraction. The
    bound is rigorous; the loss-of-significance is correctly tracked, not
    hidden.
  - The conditional `result + tp` adds another δ(2π) to precision; this
    is conservative but rigorous (the actual rounded value of result+tp
    differs from the true principal value by at most one more ε(2π)).
  - The half-open vs closed-open boundary at θ=2π is mathematically
    [0, 2π) but fmod's behavior at exact multiples returns ±0 — the
    sign branch maps −0 to +2π via the conditional, which is the
    correct principal value at the boundary.
```

---

## Card 6 — `angles::wrap_neg_pos_pi`

```
=== FORMULA AUDIT CARD ===
ID:                     angles::wrap_neg_pos_pi
Location:               src/math/angles.h:50-55
Mathematical statement: Wrap angle θ to the half-open interval [−π, π);
                        result = wrap_two_pi(θ + π) − π.

THEORY
  Underlying theorem:   Shift-and-wrap identity:
                          ((θ + π) mod 2π) − π ∈ [−π, π)
                        for all real θ. This is the canonical reduction
                        used in numerical libraries (cf. C99
                        `remainder(x, 2π)`-style range).
                        Composition: addition is closed-form (REQ-EF-3),
                        wrap_two_pi gives [0, 2π) by Card 5, and a final
                        subtraction of π shifts to [−π, π). Each step is
                        a closed-form TrackedValue operation, so error
                        propagation chains rigorously.
  Primary reference:    ProofWiki — Real Modulo and Principal Value;
                        Goldberg (1991) Theorem 11; equivalent to
                        std::remainder(x, 2π) up to sign convention at
                        the open endpoint.
  Domain of validity:   All finite real θ.

METHOD
  Method declared:      Closed-form composition: `wrap_two_pi(θ + π) − π`.
  Method implemented:   `auto p = pi<T>();`
                        `auto result = wrap_two_pi(angle + p) - p;`
  Match verdict:        ✓ matched — exact closed-form composition of
                        REQ-EF-3 operations; no series, no iteration.

ERROR BOUND
  Bound category:       precision + measurement (propagated).
  Bound formula:        Compose the per-step REQ-EF-3 bounds:
                          δ(θ + π)   = δθ + δπ
                          δ(wrap)    ≤ δ(θ+π) + |⌊(θ+π)/2π⌋|·δ(2π)
                          δ(result)  = δ(wrap) + δπ.
                        Net: δ ≤ δθ + 2·δπ + n·δ(2π), with n = the
                        integer quotient inside fmod.
  Bound implemented:    Inherited from `operator+` (tracked_value.h),
                        `wrap_two_pi` (Card 5), and `operator-`. The
                        composition adds the appropriate per-step bound
                        to each error category — REQ-EF-3 throughout.
  Bound verdict:        ✓ matched — bound is the closed-form composition
                        of REQ-EF-3 propagation on three operations.
                        DEPENDS on the tracked_value::fmod bound (per task
                        note) since wrap_two_pi is invoked internally;
                        same rigorous Lipschitz-1 chain as Card 5.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3.
  AUD-EF applies:       AUD-EF-2 (binary-op wiring across the three steps).
  AUD-MC applies:       n/a.
  Verification test:    No dedicated test currently; recommend
                        θ ∈ {0, π, −π, 2π, −2π, 3π/2, −3π/2, very-large}.

NOTES
  - At θ = π the result is mathematically the open endpoint and falls into
    [−π, π) as +(-π) or −π depending on floating-point representation of
    π · (1 + ε). The bound δ ≥ 2·δπ ensures the reported precision
    covers this boundary uncertainty.
  - Like Card 5, large |θ| amplifies δ(2π) by the integer quotient n
    inside the internal fmod — rigorous tracking of argument-reduction
    loss-of-significance, not a bound failure.
  - The interval `[−π, π)` is the principal-value choice for unique
    representation; the alternative `(−π, π]` would require an extra
    sign test but the project standard is the half-open `[)` form.
```

---

## File-level verdict for `angles.h`

- **A. Error wiring**: ✓ all six entry points return `TrackedValue<T>` and propagate via REQ-EF-3 closed-form rules (Cards 1–6 confirm AUD-EF-1, AUD-EF-2).
- **B. Algebra axioms**: n/a — utilities, not algebra ops. (Used downstream in quaternion/dual-quaternion code where AUD-MC tests apply.)
- **C. Theoretical basis**:
  - Card 1 `pi`: ✓ all three slots matched. **PASS.**
  - Card 2 `two_pi`: ✓ all three slots matched. **PASS.**
  - Card 3 `degrees_to_radians`: ✓ all three slots matched. **PASS.**
  - Card 4 `radians_to_degrees`: ✓ all three slots matched. **PASS.**
  - Card 5 `wrap_two_pi`: ✓ all three slots matched, **conditioned on** the rigorous `tracked_value::fmod` bound at `tracked_value.h:401-410`. **PASS.**
  - Card 6 `wrap_neg_pos_pi`: ✓ all three slots matched, **conditioned on** Card 5 (transitively on tracked_value::fmod). **PASS.**

**File verdict: PASS** — every formula is closed-form (no Taylor or iterative method to mismatch), every bound is the REQ-EF-3 closed-form propagation result, and the two wrap functions cleanly compose the rigorous fmod bound supplied by `tracked_value.h`. No C-fails. Recommended follow-up: add dedicated round-trip and wrap-boundary tests in `tests/test_math/main.cpp`.
