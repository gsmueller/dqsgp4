# Theoretical Basis Audit — `src/math/wallis.h`

**File**: `src/math/wallis.h` (59 lines)
**Expected functions**: 4 (wallis_odd, wallis_even, wallis, sin_power_cos_integral)
**Audit status**: PASS
**Summary**: All four functions implement closed-form exact identities (Wallis' reduction formulas for power-of-cosine integrals and a direct rational integral). No series approximation, no iteration. Error propagation via closed-form rules (REQ-EF-3).

---

## 4.1 `wallis_odd(k)` — Odd-power Wallis integral

`
=== FORMULA AUDIT CARD ===
ID:                     wallis::wallis_odd
Location:               src/math/wallis.h:22-24
Mathematical statement: W_{2k+1} = ∫_0^{π/2} cos^{2k+1}(φ) dφ = (2k)!! / (2k+1)!!

THEORY
  Underlying theorem:   Wallis' reduction formula (Wallis 1656, *Arithmetica
                        Infinitorum*; modern treatment: Whittaker & Watson
                        (1927) §12.5 or ProofWiki "Integral of Power of
                        Cosine"). For n = 2k+1 (odd), the repeated application
                        of integration by parts reduces the power recursively,
                        yielding the closed form:
                          W_{2k+1} = (2k)!! / (2k+1)!!
                        where n!! = n × (n−2) × (n−4) × ⋯ × 1.
                        This is an exact rational number (no π).
  Primary reference:    Whittaker & Watson (1927) §12.5, or ProofWiki
                        "Integral of Power of Cosine / Odd Index Form".
  Domain of validity:   All k ≥ 0. The double-factorial form is exact.

METHOD
  Method declared:      Closed-form identity via double-factorial ratio.
  Method implemented:   `double_factorial<T>(2*k) / double_factorial<T>(2*k+1)`
  Match verdict:        ✓ matched — direct evaluation of the closed form.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Closed-form identity with no approximation. Error
                        propagates from the two `double_factorial<T>()` calls
                        via REQ-EF-3 (division of two TrackedValues):
                          errors.precision = (|1/v| · err_num + |u/v²| · err_denom)
                        where u = double_factorial(2k), v = double_factorial(2k+1),
                        err_num = u.errors.precision, err_denom = v.errors.precision.
  Bound implemented:    `operator/(TrackedValue, TrackedValue)` in
                        tracked_value.h applies REQ-EF-3 division rule.
  Bound verdict:        ✓ matched — propagated via closed-form rule.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form error propagation)
  AUD-EF applies:       AUD-EF-3 (Binary operators add per REQ-EF-3)
  AUD-MC applies:       n/a (not an algebra operation)
  Verification test:    tests/test_math — wallis_odd(k) for k=0,1,2 should
                        match known values (1, 2/3, 8/15, …).

NOTES
  - The denominator 2k+1!! is never zero (all terms are positive integers).
  - Both double_factorial calls return exact integer ratios; their composition
    is also exact rational (no accumulated truncation).
  - This is the reference implementation; no alternative method.
```

---

## 4.2 `wallis_even(k)` — Even-power Wallis integral

`
=== FORMULA AUDIT CARD ===
ID:                     wallis::wallis_even
Location:               src/math/wallis.h:29-39
Mathematical statement: W_{2k} = ∫_0^{π/2} cos^{2k}(φ) dφ = (2k−1)!! / (2k)!! × π/2

THEORY
  Underlying theorem:   Wallis' reduction formula (same as 4.1) applied to
                        even powers n = 2k. The recurrence reduces to:
                          W_{2k} = ((2k−1)!! / (2k)!!) × (π/2)
                        The π/2 factor arises from the base case W_0 = π/2.
  Primary reference:    Whittaker & Watson (1927) §12.5, or ProofWiki
                        "Integral of Power of Cosine / Even Index Form".
  Domain of validity:   All k ≥ 0. For k = 0, W_0 = π/2 exactly.

METHOD
  Method declared:      Closed-form: (2k−1)!! / (2k)!! × π/2, with base case
                        W_0 = π/2 computed from boost::math::constants::pi<T>()
                        with precision bound = π_val × ε.
  Method implemented:   For k=0: construct π/2 with precision bound.
                        For k>0: recursively call wallis_even(0) for π/2,
                        then multiply by (2k−1)!! / (2k)!!.
  Match verdict:        ✓ matched — direct closed-form evaluation.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Base case (k=0): π_val is from boost with machine
                        epsilon relative error. Precision bound assigned:
                          trunc_bound = π_val × ε
                        For k > 0: π/2 (from base case) is multiplied by the
                        double-factorial ratio. Error propagates via REQ-EF-3
                        multiplication rule.
  Bound implemented:    k=0: `T pi_prec = pi_val * std::numeric_limits<T>::epsilon();`
                        k>0: `double_factorial(2k−1) / double_factorial(2k) * pi_half`
                        propagates via `operator*(TrackedValue, TrackedValue)`.
  Bound verdict:        ✓ matched — base case uses machine epsilon as
                        conservative bound for boost's π; composition uses
                        closed-form multiplication rule.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form rules for *, /)
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    tests/test_math — wallis_even(k) for k=0,1,2 should
                        match known values (π/2, π/8, π/128, …).

NOTES
  - The recursion wallis_even<T>(0) is re-evaluated each call for k>0.
    For efficiency, the compiler may inline and cache, but the code does not
    explicitly memoize. This is acceptable: boost::math::constants::pi<T>()
    is typically instantiated once at compile time.
  - The precision bound for π (ε relative) is conservative: boost's pi<T>()
    is typically correct to full machine precision. Using ε is safe and simple.
  - No alternative method; this is the canonical form.
```

---

## 4.3 `wallis(n)` — Dispatcher

`
=== FORMULA AUDIT CARD ===
ID:                     wallis::wallis
Location:               src/math/wallis.h:43-49
Mathematical statement: W_n = { W_{2k} if n=2k, W_{2k+1} if n=2k+1 }
                        i.e. ∫_0^{π/2} cos^n(φ) dφ for any n ≥ 0.

THEORY
  Underlying theorem:   Wallis' reduction formula for arbitrary n.
                        The closed form depends on parity of n.
  Primary reference:    (same as 4.1 and 4.2)
  Domain of validity:   All n ≥ 0.

METHOD
  Method declared:      Dispatcher: if n is even, call wallis_even(n/2);
                        otherwise call wallis_odd(n/2).
  Method implemented:   `if (n % 2 == 0) return wallis_even<T>(n/2);
                        else return wallis_odd<T>(n/2);`
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision
  Bound formula:        No independent error; propagates from the called
                        function (wallis_even or wallis_odd).
  Bound implemented:    Pass-through: result.errors = called_result.errors.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       (inherited from wallis_even or wallis_odd)
  AUD-EF applies:       (inherited)
  AUD-MC applies:       n/a
  Verification test:    (inherited)

NOTES
  - Pure wrapper; no arithmetic on TrackedValue in this function.
  - Integer parity check (n % 2) is exact and adds no error.
```

---

## 4.4 `sin_power_cos_integral(k)` — Direct rational integral

`
=== FORMULA AUDIT CARD ===
ID:                     wallis::sin_power_cos_integral
Location:               src/math/wallis.h:55-57
Mathematical statement: I_k = ∫_0^{π/2} sin^{2k}(φ) cos(φ) dφ = 1/(2k+1)

THEORY
  Underlying theorem:   Direct substitution u = sin(φ), du = cos(φ) dφ gives
                          I_k = ∫_0^1 u^{2k} du = [u^{2k+1} / (2k+1)]_0^1 = 1/(2k+1)
                        This is an immediate consequence of the power-rule for
                        integration and the fundamental theorem of calculus.
                        Exact closed form; no approximation.
  Primary reference:    Any calculus textbook; direct application of u-substitution.
  Domain of validity:   All k ≥ 0. The result is an exact rational.

METHOD
  Method declared:      Closed-form rational.
  Method implemented:   `ratio<T>(1, 2*k+1)`
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Closed-form identity. Error propagates from
                        `ratio<T>(1, 2*k+1)` via REQ-EF-3 (rational construction).
  Bound implemented:    The `ratio<T>` helper constructs an exact rational
                        with relative error bound = ε (machine epsilon, or
                        conservative bound for the rounded representation).
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form rational)
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    tests/test_math — sin_power_cos_integral(k) for
                        k=0,1,2 should match 1, 1/3, 1/5, ….

NOTES
  - The denominator 2k+1 is always positive and nonzero.
  - The numerator 1 is exact in all T.
  - This integral appears in the derivation of gravity and J₂ terms;
    it is simpler than the full Wallis integrals and often preferred.
```

---

## File-level verdict for `wallis.h`

- **A. Error wiring**: ✓ All four functions return `TrackedValue<T>`. Error propagation via:
  - wallis_odd: division (REQ-EF-3)
  - wallis_even: multiplication and division (REQ-EF-3)
  - wallis: pass-through (inherited)
  - sin_power_cos_integral: rational constructor (REQ-EF-3)

- **B. Algebra axioms**: n/a (utility functions, not algebra operations).

- **C. Theoretical basis**: ✓ All four cards verified:
  - 4.1 wallis_odd: ✓ closed-form double-factorial identity matches Wallis' formula.
  - 4.2 wallis_even: ✓ closed-form (2k−1)!!/(2k)!! × π/2 matches Wallis' formula.
  - 4.3 wallis: ✓ dispatcher matches parity case split.
  - 4.4 sin_power_cos_integral: ✓ direct substitution integral matches u-sub closed form.

**File verdict: PASS** — all formulas are exact closed-form identities with no approximation error. Bounds propagate correctly via REQ-EF-3.
