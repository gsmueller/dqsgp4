# Theoretical Basis Audit — `src/perturbation/kaula.h`

**File**: `src/perturbation/kaula.h` (176 lines)
**Functions producing `TrackedValue<T>`**: 1 (`inclination_function`)
**Expected count**: 1 ✓

The single template dispatches on `(l, m, p)` to one of 12 hand-coded closed-form polynomial branches, plus a fallback. Each branch is an independent formula and is audited as its own card. The fallback (no match) is audited as card 13.

---

## 1 `inclination_function(2, 2, 0, sin_i, cos_i)` — F_{220}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_220
Location:               src/perturbation/kaula.h:76-79
Mathematical statement: F_{220}(i) = (3/4) (1 + cos i)²

THEORY
  Underlying theorem:   Kaula's inclination function F_{lmp}(i), defined by
                        the finite double sum (Kaula 1966 Eq. 3.61):
                          F_{lmp}(i) = Σ_t (2l−2t)! / (t!(l−t)!(l−m−2t)!2^(2l−2t))
                                       · sin^(l−m−2t) i
                                       · Σ_s C(m,s) cos^s i · ...
                        At (l=2, m=2, p=0) the double sum collapses to a
                        single term yielding the exact closed-form rational
                        polynomial (3/4)(1+cos i)².
  Primary reference:    Kaula (1966) "Theory of Satellite Geodesy" §3.3
                        Table 1 (or Eq. 3.61 evaluated at l=m=2, p=0).
                        Also: Allan (1965) closed-form expansion for low l.
                        SGP4 cross-reference: Hoots & Roehrich (1980)
                        Spacetrack Report No. 3, p. 62-63 `f220` line.
  Domain of validity:   i ∈ [0, π]; sin_i, cos_i satisfy sin²+cos²=1
                        (caller responsibility — not enforced here).

METHOD
  Method declared:      Exact closed-form polynomial in cos i with rational
                        coefficient 3/4. NOT a Taylor expansion, NOT Kaula's
                        recursion, NOT a series — direct evaluation of the
                        algebraic identity (3/4)(1+cos i)².
  Method implemented:   `ratio<T>(3, 4) * one_plus_cos * one_plus_cos`
                        where `one_plus_cos = exact<T>(1) + cos_i`.
                        Matches SGP4 `f220 = 0.75 * (1 + 2 cos i + cos²i)`
                        algebraically (expanded form vs factored form).
  Match verdict:        ✓ matched — closed-form rational polynomial, exact.

ERROR BOUND
  Bound category:       precision (propagated from sin_i, cos_i inputs)
  Bound formula:        REQ-EF-3 closed-form propagation. For
                        f(x) = (3/4)(1+x)², propagated bound is
                        |3/2 (1+cos i)| · ε(cos_i) plus rounding of the
                        multiplications by `exact` / `ratio` (handled by
                        `TrackedValue::operator*`).
  Bound implemented:    Inherited from `TrackedValue::operator+` and
                        `operator*` chain (see `tracked_value.h`). No
                        explicit truncation bound is added here — none is
                        needed because the formula is exact-rational.
  Bound verdict:        ✓ matched — no method-introduced error;
                        closed-form propagation only.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form propagation)
  AUD-EF applies:       n/a (no truncation; AUD-EF-3 covers `*`/`+`)
  AUD-MC applies:       n/a (numerical helper)
  Verification test:    Magic-number cross-check vs SGP4 `f220` —
                        file header documents zero residual.

NOTES
  - The branch `if (l==2 && m==2 && p==0)` is reached recursively by the
    F_441 branch (kaula.h:119) via `inclination_function(2,2,0,...)`.
```

---

## 2 `inclination_function(2, 2, 1, sin_i, cos_i)` — F_{221}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_221
Location:               src/perturbation/kaula.h:83-85
Mathematical statement: F_{221}(i) = (3/2) sin²(i)

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 evaluated at (l=2, m=2, p=1).
  Primary reference:    Kaula (1966) §3.3 Table 1.
                        SGP4: Hoots & Roehrich (1980) p.62 `f221`.
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form (3/2) sin²(i).
  Method implemented:   `ratio<T>(3, 2) * sin2` where `sin2 = sin_i*sin_i`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3; bound = |3 sin i| · ε(sin_i).
  Bound implemented:    Inherited from `*` chain.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f221 = 1.5 * sini2` magic-number check.

NOTES
  - Exact-rational coefficient; no FORTRAN truncation artifact.
```

---

## 3 `inclination_function(3, 1, 1, ...)` — F_{311}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_311
Location:               src/perturbation/kaula.h:90-93
Mathematical statement: F_{311}(i) = (15/16) sin²(i) (1 + 3 cos i)
                                    − (3/4)(1 + cos i)

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=3, m=1, p=1).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f311`.
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form polynomial with rationals
                        15/16 and 3/4.
  Method implemented:   `ratio<T>(15,16) * sin2 * (1 + 3 cos_i)
                         - ratio<T>(3,4) * (1 + cos_i)`. Matches SGP4
                        `f311` per the code comment.
  Match verdict:        ✓ matched (R14, 2026-05-13) — the two-term form
                        is the canonical un-normalized Kaula F_{311},
                        cross-verified against Wakker (2015) Table 23.2,
                        Allan (1965) App. A, and SR3 p. 62 `f311`.
                        Inline spot-check comment in kaula.h:87-114
                        documents the algebraic equivalence + numerical
                        spot-check at i = 51.6° (ISS). Resolved from ? to ✓.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f311` magic-number check (file header
                        asserts zero residual against legacy code).

NOTES
  - **Resolved 2026-05-13 (R14)**: the stale one-term annotation was
    replaced with an explicit SPOT-CHECK comment block (kaula.h:87-114)
    citing Wakker (2015) Table 23.2, Allan (1965) App. A, and SR3 p. 62.
    The two-term SGP4 form is the canonical un-normalized Kaula F_{311}.
    Numerical spot-check at i = 51.6° documented inline.
```

---

## 4 `inclination_function(3, 2, 1, ...)` — F_{321}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_321
Location:               src/perturbation/kaula.h:97-100
Mathematical statement: F_{321}(i) = (15/8) sin(i) (1 − 2 cos i − 3 cos²i)

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=3, m=2, p=1).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f321`.
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form with rational 15/8.
  Method implemented:   `ratio<T>(15,8) * sin_i * (1 − 2 cos_i − 3 cos²i)`.
  Match verdict:        ✓ matched (algebraically identical to SGP4
                        `f321 = 1.875 * sinim * (1 − 2 cosim − 3 cosisq)`).

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f321` magic-number check.
```

---

## 5 `inclination_function(3, 2, 2, ...)` — F_{322}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_322
Location:               src/perturbation/kaula.h:104-107
Mathematical statement: F_{322}(i) = −(15/8) sin(i) (1 + 2 cos i − 3 cos²i)

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=3, m=2, p=2).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f322`.
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form with rational −15/8.
  Method implemented:   `ratio<T>(-15,8) * sin_i * (1 + 2 cos_i − 3 cos²i)`.
  Match verdict:        ✓ matched (SGP4 `f322 = -1.875 * sinim * ...`).

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f322` magic-number check.
```

---

## 6 `inclination_function(3, 3, 0, ...)` — F_{330}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_330
Location:               src/perturbation/kaula.h:111-114
Mathematical statement: F_{330}(i) = (15/8) (1 + cos i)³

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=3, m=3, p=0).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f330` (synchronous
                        resonance term).
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form (15/8)(1+cos i)³.
  Method implemented:   `ratio<T>(15,8) * opc * opc * opc` with
                        `opc = 1 + cos_i`.
  Match verdict:        ✓ matched (SGP4 `f330 = 1.875 * (1+cosim)^3`).

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f330` magic-number check.
```

---

## 7 `inclination_function(4, 4, 1, ...)` — F_{441}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_441
Location:               src/perturbation/kaula.h:118-121
Mathematical statement: F_{441}(i) = 35 sin²(i) · F_{220}(i)
                                   = 35 · sin²(i) · (3/4)(1+cos i)²
                                   = (105/4) sin²(i) (1+cos i)²

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=4, m=4, p=1). The
                        algebraic identity F_{441} = 35 sin²i · F_{220}
                        is asserted by the SGP4 code (`f441 = 35 * sini2
                        * f220`); it factors out the F_{220} sub-pattern
                        rather than expanding the Kaula sum independently.
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f441`.
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form, structured as 35·sin²i·F_{220}
                        with F_{220} obtained by recursive call.
  Method implemented:   `exact<T>(35) * sin2 * inclination_function(2,2,0,...)`
                        — recurses into card 1.
  Match verdict:        ✓ matched (SGP4 form). The recursion is purely
                        compile-time / direct call into another branch
                        of the same function; not Kaula's general
                        recursion, just code reuse.

ERROR BOUND
  Bound category:       precision (propagated through F_{220} card 1)
  Bound formula:        REQ-EF-3 chained through recursive call.
  Bound implemented:    Inherited from F_{220} return + outer `*` chain.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f441` magic-number check.

NOTES
  - The recursive call is the only branch reuse in the file. Cycle-free
    because the recursion target (l=2) cannot recurse back.
```

---

## 8 `inclination_function(4, 4, 2, ...)` — F_{442}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_442
Location:               src/perturbation/kaula.h:125-127
Mathematical statement: F_{442}(i) = (315/8) sin⁴(i)

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=4, m=4, p=2).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f442` (which uses
                        the truncated decimal 39.3750 = 315/8).
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form (315/8) sin⁴ i.
  Method implemented:   `ratio<T>(315, 8) * sin2 * sin2`.
  Match verdict:        ✓ matched; FORTRAN truncation 39.3750 is replaced
                        with the exact rational 315/8.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f442` magic-number check.

NOTES
  - 315/8 = 39.375 exactly; the SGP4 literal `39.3750` already coincides
    with the rational. This branch eliminates rounding-on-construction
    by going through `ratio<T>(315, 8)` rather than a double literal.
```

---

## 9 `inclination_function(5, 2, 2, ...)` — F_{522}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_522
Location:               src/perturbation/kaula.h:132-137
Mathematical statement: F_{522}(i) = (315/32) sin(i) · [ sin²i (1−2 cos i−5 cos²i)
                                     + (1/3)(−2 + 4 cos i + 6 cos²i) ]

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=5, m=2, p=2).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f522` (which uses
                        the truncated decimal 9.84375 = 315/32).
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form with rational 315/32 and an
                        inner 1/3 coefficient on the additive bracket.
  Method implemented:   `ratio<T>(315,32) * sin_i * ( sin2 * (1−2cos_i
                         −5 cos2) + ratio<T>(1,3) * (−2+4 cos_i+6 cos2) )`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f522` magic-number check.

NOTES
  - 9.84375 = 315/32 exactly; the `12` "trailing" suffix in some FORTRAN
    listings (per file header note) is the truncation artifact this
    module eliminates.
```

---

## 10 `inclination_function(5, 2, 3, ...)` — F_{523}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_523
Location:               src/perturbation/kaula.h:142-147
Mathematical statement: F_{523}(i) = sin(i) · [ (315/64) sin²i (−2 − 4 cos i
                                     + 10 cos²i) + (105/16)(1 + 2 cos i − 3 cos²i) ]

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=5, m=2, p=3).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f523` (decimals
                        4.92187512 ≈ 315/64 and 6.56250012 ≈ 105/16).
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form with rationals 315/64, 105/16.
  Method implemented:   `sin_i * ( ratio<T>(315,64) * sin2 * (-2 - 4 cos_i
                         + 10 cos2) + ratio<T>(105,16) * (1 + 2 cos_i
                         - 3 cos2) )`.
  Match verdict:        ✓ matched. **FORTRAN artifact eliminated**:
                        4.92187512 → 315/64 = 4.921875 (the trailing
                        "12" was a FORTRAN truncation, removed here);
                        6.56250012 → 105/16 = 6.5625 likewise.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched. The rational replacement removes a
                        constant rounding contribution of order
                        ~1e-8 / coefficient that was previously baked
                        into SGP4's f523 — a strict improvement in the
                        precision claim.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f523` magic-number check; expect the
                        residual against exact rationals to be exactly
                        the truncation of 4.92187512 vs 4.921875,
                        amplified by sin²i.
```

---

## 11 `inclination_function(5, 4, 2, ...)` — F_{542}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_542
Location:               src/perturbation/kaula.h:152-156
Mathematical statement: F_{542}(i) = (945/32) sin(i) · [2 − 8 cos i
                                     + cos²i (−12 + 8 cos i + 10 cos²i)]

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=5, m=4, p=2).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f542` (decimal
                        29.53125 = 945/32).
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form with rational 945/32.
  Method implemented:   `ratio<T>(945,32) * sin_i * (2 - 8 cos_i + cos2
                         * (-12 + 8 cos_i + 10 cos2))`. Horner-style
                        nesting in cos²i.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f542` magic-number check.

NOTES
  - 29.53125 = 945/32 exactly; SGP4 decimal happens to be exactly
    representable in double, so the residual against the rational
    form is zero by construction.
```

---

## 12 `inclination_function(5, 4, 3, ...)` — F_{543}

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::F_543
Location:               src/perturbation/kaula.h:161-165
Mathematical statement: F_{543}(i) = (945/32) sin(i) · [−2 − 8 cos i
                                     + cos²i (12 + 8 cos i − 10 cos²i)]

THEORY
  Underlying theorem:   Kaula 1966 Eq. 3.61 at (l=5, m=4, p=3).
  Primary reference:    Kaula (1966) §3.3 Table 1; SGP4 `f543` (decimal
                        29.53125 = 945/32).
  Domain of validity:   i ∈ [0, π].

METHOD
  Method declared:      Exact closed-form with rational 945/32.
  Method implemented:   `ratio<T>(945,32) * sin_i * (−2 − 8 cos_i + cos2
                         * (12 + 8 cos_i − 10 cos2))`.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       precision (propagated)
  Bound formula:        REQ-EF-3.
  Bound implemented:    Inherited.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  Verification test:    SGP4 `f543` magic-number check.
```

---

## 13 Fallback path (no `(l,m,p)` branch matches)

```
=== FORMULA AUDIT CARD ===
ID:                     kaula::inclination_function::fallback_zero
Location:               src/perturbation/kaula.h:170-172
Mathematical statement: F_{lmp}(i) returns 0 with errors.accuracy = 1
                        when (l,m,p) is outside the supported set
                        { (2,2,0), (2,2,1), (3,1,1), (3,2,1), (3,2,2),
                          (3,3,0), (4,4,1), (4,4,2), (5,2,2), (5,2,3),
                          (5,4,2), (5,4,3) }.

THEORY
  Underlying theorem:   n/a — sentinel return, not a numerical method.
                        The exact F_{lmp} for unhandled indices would
                        be computed via Kaula's general double sum
                        (Kaula 1966 Eq. 3.61) or Gooding & King (1989)
                        recursion. This branch is a "not implemented"
                        sentinel, not an approximation.
  Primary reference:    n/a.
  Domain of validity:   Outside the 12 hand-coded SGP4-required cases.

METHOD
  Method declared:      Sentinel return of 0 with accuracy=1 to signal
                        "not implemented".
  Method implemented:   `result = exact<T>(0); result.errors.accuracy = T(1);
                         return result;`.
  Match verdict:        ⚠ sentinel semantics — calling code must treat a
                        unity `accuracy` as a fatal signal, not a tight
                        bound on 0. This is not a numerical bound — it
                        is an "absent value" marker.

ERROR BOUND
  Bound category:       accuracy (set to 1)
  Bound formula:        n/a — sentinel, not a method bound.
  Bound implemented:    `errors.accuracy = T(1)` (literal 1.0).
  Bound verdict:        ⚠ semantic — the value `1` is large enough to
                        force any caller's `total_error()` check to
                        reject the computation. Acceptable as a signal,
                        but is not a quantitative bound on a quantity.

CROSS-AUDIT
  REQ-EF applies:       n/a (sentinel, not a propagated bound)
  AUD-EF applies:       caller must test `errors.accuracy >= 1` for
                        the not-implemented case.
  Verification test:    A test calling, e.g., `inclination_function(6, 0, 0, ...)`
                        should observe accuracy=1 and value=0.

NOTES
  - **Open issue (C-warning)**: this fallback departs from REQ-EF
    semantics by overloading the `accuracy` slot to mean "function
    not defined for these arguments" rather than "rigorous accuracy
    bound on the returned value". Callers must understand this
    overload. A cleaner design would be a `std::optional` return or
    a dedicated `errors.not_implemented` flag.
  - For SGP4's resonance perturbations, only the 12 cards above are
    needed; this branch should never fire under normal usage. Tests
    should still confirm the sentinel behavior.
  - Header documentation mentions extension via "Gooding & King (1989)
    recursion for larger l" — this is forward-looking, not currently
    implemented.
```

---

## File-level verdict

**A. Error wiring**: ✓ all branches return a `TrackedValue<T>` constructed via `ratio`, `exact`, and `TrackedValue` arithmetic; bounds are propagated through `tracked_value.h` operators (REQ-EF-3 closed-form path). No method introduces truncation, so no truncation bound is added.

**B. Algebra axioms**: n/a — this is a numerical helper, not an algebra operation. Its use in resonance forcing downstream is exercised by the higher-level perturbation tests.

**C. Theoretical basis**:
- All 12 numerical branches: ✓ implementation is **exact closed-form polynomial in (sin i, cos i) with exact rational coefficients**, not a Taylor expansion, not Kaula's recursion, not a series. The theory citation (Kaula 1966 Table 1) and method (closed-form polynomial) match. The bound is REQ-EF-3 closed-form propagation, which is the rigorous bound for an exact-rational closed form. **PASS** on the central audit question — no Taylor/Padé/continued-fraction confusion is possible because the operation has no truncation.
- F_{311} (card 3): ✓ resolved by R14 (2026-05-13). The two-term SGP4 form is the canonical un-normalized Kaula F_{311}, cross-verified against Wakker (2015) Table 23.2, Allan (1965) App. A, and SR3 p. 62. Numerical spot-check at i = 51.6° documented inline in kaula.h:87-114.
- Fallback (card 13): ⚠ overloads `errors.accuracy = 1` as a "not implemented" sentinel rather than a quantitative bound. Semantically clear but cross-cuts REQ-EF.

**File verdict: PASS** with one ⚠ design note (fallback sentinel semantics). The R14 ? primary-source spot-check for F_{311} is closed (2026-05-13). The file's central improvement over SGP4 — replacing FORTRAN-truncated decimals (4.92187512, 6.56250012) with exact rationals (315/64, 105/16) — strengthens, not weakens, the precision claim, and the theory→method→bound chain is sound for all 12 implemented cases.
