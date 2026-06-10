# Theoretical Basis Audit — src/sgp4/precomputed.h

**File**: src/sgp4/precomputed.h (115 lines)
**Status**: PASS
**Audit date**: 2026-05-13

## Overview

`KaulaTable<T>` is a fixed-capacity lookup table caching precomputed Kaula inclination function values F_{l,m,p}(i) keyed by integer triples (l, m, p). Because the orbital inclination is fixed at epoch, F_{lmp}(i) is constant for the life of a propagator instance and is computed once at TLE-driven construction.

Four functions are audited:

1. `set(l, m, p, val)` — container utility (insert).
2. `F(l, m, p)` — container utility (lookup with zero default).
3. `has(l, m, p)` — container utility (presence test).
4. `precompute(sin_i, cos_i)` — driver that evaluates `perturbation::inclination_function(l, m, p, sin_i, cos_i)` for the 12 (l,m,p) triples needed by SGP4 deep-space and stores the results.

Cards 1-3 are container utilities that perform no numeric computation; they store, retrieve, and test membership of values that arrive already wrapped in `TrackedValue<T>`. Their theoretical basis is structural per REQ-EF-12 (composite types compose).

Card 4 (`precompute`) is the substantive theoretical-basis card: it delegates each F_{lmp}(i) evaluation to `perturbation::inclination_function`. The Kaula (1966) theoretical basis lives in `perturbation/kaula.h` and is audited via that file. This card's responsibility is to confirm that (a) the list of (l,m,p) triples matches what the SGP4 deep-space model requires and (b) the delegation transmits errors faithfully.

---

## Function Audit Cards

### CARD 1: KaulaTable::set(l, m, p, val)

```
=== FORMULA AUDIT CARD ===
ID:                     precomputed::KaulaTable::set
Location:               src/sgp4/precomputed.h:61-65
Mathematical statement: Insert entry (l, m, p, val) at position `count`, increment count.
                        Silently drops entries beyond MAX_ENTRIES = 64.

THEORY
  Underlying theorem:   None (container operation; structural).
                        Capacity-bounded append on a flat array.
  Primary reference:    REQ-EF-12 (composite types compose). The TrackedValue<T>
                        passed in retains its own theoretical-basis card from
                        its producer (perturbation::inclination_function).
  Domain of validity:   count < MAX_ENTRIES = 64. SGP4 needs ~12 entries
                        (well below capacity).

METHOD
  Method declared:      Direct array write at index `count`, then increment.
                        Caller TrackedValue<T> copied by value.
  Method implemented:   `entries[count++] = {l, m, p, val};` guarded by
                        `if (count < MAX_ENTRIES)`.
  Match verdict:        PASS — pure container operation, no formula.

ERROR BOUND
  Bound category:       None (no computation; pass-through).
  Bound formula:        N/A — `val.errors` are copied unchanged into the
                        stored Entry.
  Bound implemented:    Aggregate-init copy preserves all error fields.
  Bound verdict:        PASS — transparent storage; no error introduced.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12 (composite types compose).
  AUD-EF applies:       AUD-EF-7 (accessor / pass-through wiring).
  AUD-MC applies:       N/A (no algebra operation).
  Verification test:    Round-trip via precompute() + F(l,m,p) (Card 4 + Card 2)
                        — the precomputed F_{220}(i) etc. must equal the
                        direct perturbation::inclination_function output.

NOTES
  - Silent drop above MAX_ENTRIES is a structural cap, not a numeric concern.
    The SGP4 model uses 12 entries; cap of 64 leaves >5x margin.
  - No bounds-check error reported (set is fail-silent). Acceptable because
    callers are constexpr lists, not runtime data.
```

---

### CARD 2: KaulaTable::F(l, m, p)

```
=== FORMULA AUDIT CARD ===
ID:                     precomputed::KaulaTable::F
Location:               src/sgp4/precomputed.h:69-76
Mathematical statement: Lookup: return the stored TrackedValue<T> matching
                        (l, m, p), or exact<T>(0) if not present.

THEORY
  Underlying theorem:   None (container lookup; structural).
                        Linear scan over a small (count <= 64) flat array.
  Primary reference:    REQ-EF-12 (composite types compose); the returned
                        TrackedValue<T> carries the error budget set by its
                        producer (perturbation::inclination_function).
  Domain of validity:   Any (l, m, p) tuple — miss returns exact zero.

METHOD
  Method declared:      Linear scan; return matching entry's value, or
                        exact<T>(0) on miss.
  Method implemented:   `for (i = 0; i < count; ++i) { if match: return
                        entries[i].value; } return math::exact<T>(0);`.
  Match verdict:        PASS — direct linear scan as declared.

ERROR BOUND
  Bound category:       None (lookup; no transformation).
  Bound formula:        Hit: returned value's errors are whatever the
                        producer stored (typically from
                        perturbation::inclination_function).
                        Miss: exact<T>(0) returns zero precision / accuracy /
                        measurement (an exact constant).
  Bound implemented:    Hit: aggregate copy preserves all error fields.
                        Miss: math::exact<T>(0) yields zero-error
                        TrackedValue.
  Bound verdict:        PASS — transparent pass-through on hit, exact zero
                        on miss.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (exact constants have zero precision error),
                        REQ-EF-12 (composite types compose).
  AUD-EF applies:       AUD-EF-3 (exact constant propagation),
                        AUD-EF-7 (accessor wiring).
  AUD-MC applies:       N/A (lookup, not algebra).
  Verification test:    F(2, 2, 0) after precompute() should equal
                        inclination_function(2, 2, 0, sin_i, cos_i)
                        identically. F(99, 99, 99) should return
                        exact zero.

NOTES
  - Miss-returns-zero is a silent default. Callers must verify presence with
    has() if a miss is semantically distinct from a true zero (in SGP4
    deep-space, missing triples are not expected — the precompute() list
    is exhaustive for the model).
  - Linear scan is O(count) <= O(64); cost is negligible vs. force
    evaluation.
```

---

### CARD 3: KaulaTable::has(l, m, p)

```
=== FORMULA AUDIT CARD ===
ID:                     precomputed::KaulaTable::has
Location:               src/sgp4/precomputed.h:79-86
Mathematical statement: Presence test: return true iff an entry matching
                        (l, m, p) is in the table.

THEORY
  Underlying theorem:   None (container membership; structural).
  Primary reference:    REQ-EF-12 (composite types compose).
  Domain of validity:   Any (l, m, p) tuple.

METHOD
  Method declared:      Linear scan; return true on match, false otherwise.
  Method implemented:   `for (i = 0; i < count; ++i) { if match: return
                        true; } return false;`.
  Match verdict:        PASS — direct linear scan as declared.

ERROR BOUND
  Bound category:       None (boolean predicate; no TrackedValue<T> result).
  Bound formula:        N/A — return type is `bool`.
  Bound implemented:    N/A.
  Bound verdict:        N/A — no error budget for bool returns.

CROSS-AUDIT
  REQ-EF applies:       N/A (bool predicate).
  AUD-EF applies:       N/A.
  AUD-MC applies:       N/A.
  Verification test:    After precompute(), has(2, 2, 0) == true;
                        has(99, 99, 99) == false.

NOTES
  - Boolean predicate; intentionally outside the TrackedValue<T> error
    framework.
  - Used to distinguish true-zero F_{lmp} from never-precomputed (since
    F() conflates them).
```

---

### CARD 4: KaulaTable::precompute(sin_i, cos_i)

```
=== FORMULA AUDIT CARD ===
ID:                     precomputed::KaulaTable::precompute
Location:               src/sgp4/precomputed.h:90-112
Mathematical statement: Evaluate F_{l,m,p}(i) at fixed inclination i for the
                        12 (l, m, p) triples
                          {(2,2,0), (2,2,1),
                           (3,1,1), (3,2,1), (3,2,2), (3,3,0),
                           (4,4,1), (4,4,2),
                           (5,2,2), (5,2,3), (5,4,2), (5,4,3)}
                        and store all 12 in a fresh KaulaTable<T>.

THEORY
  Underlying theorem:   Kaula (1966) §3.3 — F_{lmp}(i) is the inclination
                        function appearing in the spherical-harmonic
                        expansion of the Earth's gravitational potential
                        expressed in orbital elements. The closed-form
                        polynomial expressions in sin(i), cos(i) used here
                        are the standard table entries for the (l, m, p)
                        triples relevant to SGP4 deep-space resonance
                        terms (12-hour and 24-hour, l = 2..5).
  Primary reference:    Kaula, W. M. (1966), "Theory of Satellite Geodesy",
                        Ch. 3 §3.3, Table 1.
                        Hoots & Roehrich (1980), Spacetrack Report No. 3,
                        pp. 62-63.
                        Closed-form polynomials are realized in
                        src/perturbation/kaula.h — its audit card carries
                        the formula-by-formula theoretical-basis details
                        for each F_{lmp}.
  Domain of validity:   sin_i, cos_i arbitrary in [-1, 1] (consistent with
                        i in [0, pi]); the closed-form polynomials are
                        defined for all i.

METHOD
  Method declared:      For each of 12 (l, m, p) triples in a constexpr
                        list, call perturbation::inclination_function(l, m, p,
                        sin_i, cos_i) and table.set(l, m, p, val).
                        No numeric transformation of the inputs occurs in
                        this file — precompute is the driver, the formula
                        bodies live in kaula.h.
  Method implemented:   Lines 98-109: constexpr array `needed[]` of 12 LMP
                        structs; range-for loop calling
                        perturbation::inclination_function and
                        table.set(...) for each.
  Match verdict:        PASS — pure driver; formula evaluation delegated
                        verbatim to perturbation::inclination_function.

ERROR BOUND
  Bound category:       Each stored TrackedValue inherits the error
                        budget produced by perturbation::inclination_function.
                        That function composes closed-form polynomial
                        operations (mul/add/scale by ratio<T>) on the
                        TrackedValue<T> inputs sin_i, cos_i — its bounds
                        are propagated per REQ-EF-3 (closed-form
                        arithmetic propagates input errors through
                        triangle / product rules).
  Bound formula:        For each F_{lmp} = P_{lmp}(sin i, cos i):
                          err(F) <= sum over monomials of
                                    |partial P / partial sin_i| * err(sin_i)
                                    + |partial P / partial cos_i| * err(cos_i)
                                    + (rounding in mul/add at type T's
                                       precision).
                        See kaula.h audit for per-formula tightness.
  Bound implemented:    Errors flow uniformly through TrackedValue<T>
                        operators (*, +, ratio<T>, exact<T>) invoked
                        inside perturbation::inclination_function;
                        stored unchanged via set().
  Bound verdict:        PASS — bound is whatever
                        perturbation::inclination_function produces;
                        precompute itself does not add or drop error.
                        Tightness of the underlying F_{lmp} formulas is
                        audited in kaula.h.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form arithmetic propagation),
                        REQ-EF-12 (composite types compose).
  AUD-EF applies:       AUD-EF-3 (exact / closed-form constant
                        propagation), AUD-EF-7 (pass-through wiring).
  AUD-MC applies:       N/A (no algebra-axiom test at this layer;
                        kaula.h's per-F_{lmp} formulas are tested against
                        SGP4 reference values to machine precision).
  Verification test:    precompute(sin_i, cos_i).F(l, m, p) ==
                          perturbation::inclination_function(l, m, p,
                                                             sin_i, cos_i)
                        for each of the 12 triples. Numerically: legacy
                        SGP4 magic numbers (f220, f221, ..., f543) should
                        match table.F(l, m, p).value to machine precision.

NOTES
  - The 12-triple list is the complete set used by SGP4 deep-space
    perturbations (12-hour and 24-hour resonance, lunar / solar terms).
    Confirmed against Hoots & Roehrich (1980) Spacetrack Report No. 3
    pp. 62-63 and against the magic numbers in kaula.h.
  - Triple list is constexpr and resides solely in this file; if SGP4
    extensions (e.g. ?higher-order resonance) need additional (l, m, p),
    this list must be amended. has() and F() are agnostic to list
    content.
  - If perturbation::inclination_function is called with an (l, m, p)
    not in its switch (kaula.h:170), it returns exact<T>(0) with
    accuracy = 1 to flag "not implemented." All 12 triples in this
    file's list are explicitly handled in kaula.h, so no fallback path
    fires under normal use. Audit kaula.h to confirm coverage.
  - Method is delegation; no theory-method mismatch is possible at
    this layer. Any future C-fail would originate in kaula.h.
```

---

## File-Level Verdict

| Audit dimension | Result | Notes |
|---|---|---|
| **A. Error wiring** | PASS | All four functions either pass `TrackedValue<T>` through (set, F, precompute) or return a non-tracked predicate (has). No bare T leaks; no error dropped. Errors flow uniformly per REQ-EF-12. |
| **B. Algebra axioms** | N/A | KaulaTable is a container, not an algebra. The F_{lmp} formulas it caches are evaluated in `perturbation/kaula.h` and tested against SGP4 magic numbers (see kaula.h `@par Magic Number Test Results`). |
| **C. Theoretical basis** | PASS | Cards 1-3 are structural container ops (no formula). Card 4 delegates verbatim to `perturbation::inclination_function`, whose Kaula (1966) §3.3 theoretical basis is audited in kaula.h. Theory-method match is PASS at this layer; substantive correctness rests on kaula.h. |

---

## Summary

**File verdict: PASS**

`KaulaTable<T>` is a fixed-capacity lookup wrapper that caches Kaula F_{lmp}(i) values at TLE-driven construction time. The four functions split cleanly:

1. PASS — `set` — capacity-guarded array append; pure container; no theory.
2. PASS — `F` — linear-scan lookup; closed-form `exact<T>(0)` on miss; transparent pass-through on hit.
3. PASS — `has` — linear-scan boolean predicate; outside the TrackedValue<T> framework by design.
4. PASS — `precompute` — driver for the 12 (l, m, p) triples required by SGP4 deep-space; delegates each F_{lmp}(i) evaluation to `perturbation::inclination_function`. The Kaula (1966) §3.3 theoretical basis lives in `kaula.h`; this file introduces no new formula and no new error.

The substantive theoretical-basis claim — that F_{lmp}(i) is computed as the exact rational polynomial in sin(i), cos(i) — is owned by `perturbation/kaula.h`. Errors propagate from `inclination_function` through `set` and `F` unchanged. **No tightening or amendments needed at this layer.**
