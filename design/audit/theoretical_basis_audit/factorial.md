# Theoretical Basis Audit — `src/math/factorial.h`

## File summary
Five functions implementing exact and generalized factorial-like operations. All return `TrackedValue<T>` with error propagation via REQ-EF-3 (product rule). Integer inputs yield exact results (zero error); TrackedValue inputs propagate error through multiplication chains.

---

## FORMULA AUDIT CARD 1 of 5: factorial

```
=== FORMULA AUDIT CARD ===
ID:                     factorial::factorial
Location:               src/math/factorial.h:19-25
Mathematical statement: n! = ∏_{i=1}^{n} i for integer n ≥ 0

THEORY
  Underlying theorem:   Definition of factorial as the product
                        n! = 1 · 2 · 3 · ... · n.
  Primary reference:    Knuth (1997) "The Art of Computer Programming" Vol. 1, §1.2.5.
  Domain of validity:   n ∈ ℤ_≥0; for T=double, n in [0,170].

METHOD
  Method declared:      Exact integer product: iterate i from 2 to n.
  Method implemented:   result = exact<T>(1); for i in [2,n]: result *= exact<T>(i).
  Match verdict:        ✓ matched — direct definition.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Each factor is exact<T>(...), zero error. Product: zero.
  Bound implemented:    result = exact<T>(1); all multiplications exact.
  Bound verdict:        ✓ matched — zero precision error.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Product error propagation for exact inputs)
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    factorial(0)=1, factorial(5)=120, both error=0.

NOTES
  - No precision loss; pure mathematical identity.
```

---

## FORMULA AUDIT CARD 2 of 5: double_factorial

```
=== FORMULA AUDIT CARD ===
ID:                     factorial::double_factorial
Location:               src/math/factorial.h:30-37
Mathematical statement: n!! = n · (n-2) · (n-4) · ... · (2 or 1)

THEORY
  Underlying theorem:   Definition of double factorial: product of every
                        second integer down to 1 or 2.
  Primary reference:    Abramowitz & Stegun (1964) §6.1.48.
  Domain of validity:   n ∈ ℤ_≥0. Convention: 0!! = 1.

METHOD
  Method declared:      Exact product: iterate i from n down by steps of 2.
  Method implemented:   if (n ≤ 0) return exact<T>(1);
                        result = exact<T>(1);
                        for (i = n; i >= 2; i -= 2): result *= exact<T>(i).
  Match verdict:        ✓ matched — direct definition.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Each factor exact<T>(...); product zero error.
  Bound implemented:    Multiplications of exact<T> quantities.
  Bound verdict:        ✓ matched — zero precision error.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    double_factorial(5)=15, double_factorial(6)=48, error=0.

NOTES
  - Convention at n ≤ 0 returns 1; no precision loss.
```

---

## FORMULA AUDIT CARD 3 of 5: falling_factorial

```
=== FORMULA AUDIT CARD ===
ID:                     factorial::falling_factorial
Location:               src/math/factorial.h:43-50
Mathematical statement: (α)_k = α · (α−1) · (α−2) · ... · (α−k+1)

THEORY
  Underlying theorem:   Definition of falling factorial (Pochhammer negative order).
                        For integer α ≥ k: (α)_k = α! / (α−k)!.
  Primary reference:    Abramowitz & Stegun (1964) §6.1.22.
                        Graham, Knuth & Patashnik (1994) §5.4.
  Domain of validity:   k ∈ ℤ_≥0. α ∈ T or TrackedValue<T>.

METHOD
  Method declared:      Exact product: α · (α−1) · (α−2) · ... · (α−k+1).
  Method implemented:   if (k ≤ 0) return exact<T>(1);
                        result = alpha;
                        for (j = 1; j < k; ++j): result = result * (alpha - exact<T>(j)).
  Match verdict:        ✓ matched — direct definition.

ERROR BOUND
  Bound category:       precision
  Bound formula:        If α exact, each factor exact → zero error.
                        If α is TrackedValue with error δ_α, REQ-EF-3 applies:
                        error propagates through k-fold product chain.
  Bound implemented:    Each factor (alpha - exact<T>(j)) is exact subtraction;
                        product chains via *, applying REQ-EF-3.
  Bound verdict:        ✓ matched (via REQ-EF-3 multiplication operator).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Product error propagation)
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    falling_factorial(exact<T>(5), 3) = 60, error=0.
                        TrackedValue input propagates per REQ-EF-3.

NOTES
  - For integer α ≥ k, yields exact integer; for non-integer or
    TrackedValue, errors propagate multiplicatively per REQ-EF-3.
```

---

## FORMULA AUDIT CARD 4 of 5: rising_factorial

```
=== FORMULA AUDIT CARD ===
ID:                     factorial::rising_factorial
Location:               src/math/factorial.h:54-61
Mathematical statement: (α)_k^(r) = α · (α+1) · (α+2) · ... · (α+k−1)

THEORY
  Underlying theorem:   Definition of rising factorial (Pochhammer positive order).
                        Γ(α+k) / Γ(α) for α ∉ {0, −1, −2, ...}.
  Primary reference:    Abramowitz & Stegun (1964) §6.1.22.
                        Graham, Knuth & Patashnik (1994) §5.4.
  Domain of validity:   k ∈ ℤ_≥0. α ∈ T or TrackedValue<T>.

METHOD
  Method declared:      Exact product: α · (α+1) · (α+2) · ... · (α+k−1).
  Method implemented:   if (k ≤ 0) return exact<T>(1);
                        result = alpha;
                        for (j = 1; j < k; ++j): result = result * (alpha + exact<T>(j)).
  Match verdict:        ✓ matched — direct definition.

ERROR BOUND
  Bound category:       precision
  Bound formula:        If α exact, product exact (zero error).
                        If α is TrackedValue with error δ_α, REQ-EF-3:
                        error through k-fold product chain.
  Bound implemented:    Multiplication chain via TrackedValue *, REQ-EF-3.
  Bound verdict:        ✓ matched (via REQ-EF-3 multiplication).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    rising_factorial(exact<T>(2), 3) = 24, error=0.

NOTES
  - Generalization of binomial to non-integer α.
  - Used in binomial_series.h and quaternion log/exp.
```

---

## FORMULA AUDIT CARD 5 of 5: generalized_binomial

```
=== FORMULA AUDIT CARD ===
ID:                     factorial::generalized_binomial
Location:               src/math/factorial.h:65-68
Mathematical statement: C(α, k) = (α)_k / k! = falling_factorial(α, k) / k!

THEORY
  Underlying theorem:   Newton's generalized binomial coefficient.
                        Extends C(n,k) = n!/(k!(n-k)!) to non-integer α via
                        Γ(α+1) / (Γ(k+1) · Γ(α−k+1)).
  Primary reference:    Newton (1676). Modern: Abramowitz & Stegun (1964) §6.1.
                        Graham, Knuth & Patashnik (1994) §5.4.
  Domain of validity:   k ∈ ℤ_≥0. α ∈ T or TrackedValue<T>.

METHOD
  Method declared:      Composition: C(α, k) = falling_factorial(α, k) / k!.
  Method implemented:   if (k == 0) return exact<T>(1);
                        return falling_factorial(alpha, k) / factorial<T>(k).
  Match verdict:        ✓ matched — direct definition.

ERROR BOUND
  Bound category:       precision
  Bound formula:        Composition of two operations:
                        1. falling_factorial(α, k) → F with error δ_F
                        2. factorial<T>(k) → D = k! with error 0 (exact)
                        3. Division F / D: errors = |F/D| · (δ_F/|F| + 0)
                           = δ_F / |k!|
  Bound implemented:    falling_factorial(...) / factorial<T>(k) applies
                        REQ-EF-3 division after REQ-EF-3 product chain.
  Bound verdict:        ✓ matched (via composition of REQ-EF-3 operations).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Composite error: product × division)
  AUD-EF applies:       AUD-EF-3
  AUD-MC applies:       n/a
  Verification test:    C(n, k) matches classical binomial for integer n ≥ k.
                        C(α, 0) = 1 for all α, error = 0.

NOTES
  - Building block for binomial series (1+x)^α = Σ_k C(α,k) x^k.
  - Integer α ≥ k: exact classical binomial coefficient.
  - Non-integer/TrackedValue α: error through falling_factorial, then division.
```

---

## File-level verdict

- **A. Error wiring**: ✓ All five functions return `TrackedValue<T>` with REQ-EF-3 propagation (products and division).
  
- **B. Algebra axioms**: n/a (numerical helpers supporting algebra layer).

- **C. Theoretical basis**:
  - Card 1 (`factorial`): ✓ **PASS**
  - Card 2 (`double_factorial`): ✓ **PASS**
  - Card 3 (`falling_factorial`): ✓ **PASS**
  - Card 4 (`rising_factorial`): ✓ **PASS**
  - Card 5 (`generalized_binomial`): ✓ **PASS**

**File verdict: PASS** — All formulas are direct-definition products/compositions with REQ-EF-3 error propagation. No truncation, approximation, or method-theory mismatch.
