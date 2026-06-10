# Theoretical Basis Audit — `src/sgp4/model_selector.h`

**Source**: `src/sgp4/model_selector.h` (631 lines)
**Framework**: `design/audit/theoretical_basis_audit.md` (§1 card schema, §5 worked example)
**Function count**: 22

## Scope summary

This file is primarily a **configuration / dispatch layer**: it builds
`ZonalHarmonics<T>`, `ModelConfiguration<T>`, and `CustomBuilder<T>`
instances and wires them to factory presets. Most public members are
container helpers (add / has / max_degree) or builder setters that
return `*this`. They produce no numeric output and have no theoretical
basis to audit (verdict: **n/a**).

The math-relevant members are:

| # | Member | Line(s) | Verdict |
|---|---|---:|:---:|
| 1 | `ZonalHarmonics::add` | 83-85 | n/a |
| 2 | `ZonalHarmonics::Jn` | 89-94 | ✓ |
| 3 | `ZonalHarmonics::J3` | 97 | ✓ |
| 4 | `ZonalHarmonics::J5` | 98 | ✓ |
| 5 | `ZonalHarmonics::J7` | 99 | ✓ |
| 6 | `ZonalHarmonics::J9` | 100 | ✓ |
| 7 | `ZonalHarmonics::max_degree` | 103-109 | n/a |
| 8 | `ZonalHarmonics::has` | 112-117 | n/a |
| 9 | `ModelConfiguration::Jn` | 138-148 | ⚠ |
| 10 | `ModelConfiguration::A` | 167-179 | ⚠ |
| 11 | `CustomBuilder::gravity` (setter) | 191 | n/a |
| 12 | `CustomBuilder::astronomy` (setter) | 192 | n/a |
| 13 | `CustomBuilder::perturbation` (setter) | 193 | n/a |
| 14 | `CustomBuilder::kepler` (setter) | 194 | n/a |
| 15 | `CustomBuilder::drag` (setter) | 195 | n/a |
| 16 | `CustomBuilder::build` | 197-214 | n/a |
| 17 | `CustomBuilder::make_zonals` | 245-304 | ⚠ |
| 18 | `CustomBuilder::make_ellipsoid` | 309-362 | ⚠ |
| 19 | `CustomBuilder::make_astronomy` | 367-395 | ? |
| 20 | `CustomBuilder::make_model_functions` | 400-474 | ✓ |
| 21 | `ModelSelector::select` | 510-622 | n/a |
| 22 | `ModelSelector::custom` | 625-627 | n/a |

Below: full audit cards for the math-relevant members (#2-6, #9-10, #17-20).
Stubs for the others (no theoretical basis, no card required).

---

## Card 1 — `ZonalHarmonics::add` (line 83-85)

n/a — container insert. No numeric output, no theory.

---

## Card 2 — `ZonalHarmonics::Jn`

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::ZonalHarmonics::Jn
Location:               src/sgp4/model_selector.h:89-94
Mathematical statement: Jₙ(degree) = lookup in empirical table; 0 if absent.

THEORY
  Underlying theorem:   Definition — table lookup over (degree, value) pairs.
                        The Jₙ values themselves come from gravity-field
                        models (see make_zonals card for source theory).
  Primary reference:    The lookup itself: no theorem; pure data dispatch.
                        For the values returned, see make_zonals (card #17).
  Domain of validity:   Any integer degree.

METHOD
  Method declared:      Linear scan of `entries` vector; return first
                        matching value; return exact 0 otherwise.
  Method implemented:   `for (e : entries) if (e.first==degree) return e.second;`
                        else `return math::exact<T>(0);`
  Match verdict:        ✓ matched — exact match of declared method.

ERROR BOUND
  Bound category:       measurement (inherited from stored TrackedValue)
  Bound formula:        Pass-through: the returned TrackedValue carries its
                        stored measurement-error (set at add() time).
                        Fallback exact<T>(0) carries 0 in all 3 categories.
  Bound implemented:    No bound added; the stored TrackedValue is returned
                        by value, preserving its error fields.
  Bound verdict:        ✓ matched — pass-through is the rigorous form for a
                        lookup.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form / pass-through)
  AUD-EF applies:       AUD-EF-1 (TrackedValue passes through)
  AUD-MC applies:       n/a
  Verification test:    tests/test_sgp4/ — round-trip add(n,v) then Jn(n)==v.

NOTES
  - Linear scan is O(n) in the table size, but tables here are ~5-10 entries.
  - The fallback `exact<T>(0)` is correct under the documented convention:
    "Odd zonals not in the table return zero" (lines 87-88, 146-147).
```

---

## Card 3-6 — `ZonalHarmonics::J3, J5, J7, J9` (lines 97-100)

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::ZonalHarmonics::J{3,5,7,9}
Location:               src/sgp4/model_selector.h:97-100
Mathematical statement: J_k = Jn(k) for k ∈ {3, 5, 7, 9}.

THEORY
  Underlying theorem:   Definition — convenience aliases for Jn(k).
  Primary reference:    see Jn card (above) and make_zonals card.
  Domain of validity:   k ∈ {3, 5, 7, 9}.

METHOD
  Method declared:      Direct forward to Jn(k).
  Method implemented:   `return Jn(k);` for each k.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       inherited from Jn.
  Bound formula:        Pass-through.
  Bound implemented:    Pass-through.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a
  Verification test:    none specific — covered by Jn test.

NOTES
  - Pure convenience; one card covers all four.
```

---

## Card 7 — `ZonalHarmonics::max_degree` (line 103-109)

n/a — integer max over keys; no numeric output with theoretical basis.

## Card 8 — `ZonalHarmonics::has` (line 112-117)

n/a — boolean membership query; no theoretical basis.

---

## Card 9 — `ModelConfiguration::Jn`

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::ModelConfiguration::Jn
Location:               src/sgp4/model_selector.h:138-148
Mathematical statement: Resolve Jₙ for n ≥ 2 by the matched-pair principle:
                          (1) empirical table   if zonals.has(n);
                          (2) ellipsoid.J2n(n/2) if n is even, n ≥ 2;
                          (3) 0                  otherwise.

THEORY
  Underlying theorem:   (1) Matched-pair principle — the same Jₙ used during
                            TLE fitting must be used during propagation
                            (otherwise the matched-pair invariant of [SR3]
                            is violated; see Hoots-Roehrich 1980 §B).
                        (2) For even zonals not in the empirical table, the
                            level-ellipsoid analytical formula provides
                            Jₙ = -3eₚ² · (1 - n/2 + 5n/2 · J₂/eₚ²) · …
                            (Heiskanen & Moritz 1967, "Physical Geodesy",
                            §2-9, eqs. 2-92 to 2-95). The
                            EquipotentialEllipsoid::J2n method implements
                            this; treated as a black box here.
                        (3) Odd zonals (n=3,5,7,…) break equatorial symmetry
                            and cannot be obtained from the level ellipsoid
                            (which is equatorially symmetric by construction);
                            they MUST come from a satellite gravity model.
                            Default = 0 when no empirical override is loaded.
  Primary reference:    - [SR3] page 59 COMMON/C1 block (XJ3, XJ4).
                        - NGA.STND.0036 (2014) Table 3.4 (EGM2008 Cₙ₀).
                        - Heiskanen & Moritz (1967) §2-9 for the level-
                          ellipsoid J₂ₙ recurrence.
  Domain of validity:   n ≥ 2 integer. n < 2 returns 0 (undocumented but
                        consistent with "not available").

METHOD
  Method declared:      Three-branch dispatch: table → even-ellipsoid → 0.
  Method implemented:   `if (zonals.has(n)) return zonals.Jn(n);
                         if (n>=2 && n%2==0) return ellipsoid.J2n(n/2);
                         return math::exact<T>(0);`
  Match verdict:        ✓ matched — exact dispatch logic.

ERROR BOUND
  Bound category:       (1) measurement (table value);
                        (2) precision + measurement (ellipsoid output);
                        (3) zero in all categories.
  Bound formula:        Pass-through from selected branch.
  Bound implemented:    Pass-through.
  Bound verdict:        ✓ matched — no bound added at dispatch.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form / pass-through dispatch)
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a
  Verification test:    tests/test_sgp4/ — Jn(2) = ellipsoid.J2 when no
                        empirical override; Jn(3) = empirical value when
                        loaded; Jn(5) = 0 in wgs72 (not in table).

NOTES
  ⚠ Dispatch verdict: PASS, but there is a subtle convention issue. The
    EquipotentialEllipsoid::J2n returns the analytical level-ellipsoid value,
    which can differ from the empirical EGM2008 value at high degree
    (the comment on lines 51-54 notes "the ellipsoid approximation loses
    accuracy at high degree"). The choice of "ellipsoid fallback" for
    even zonals not in the table means that high-degree even zonals will
    silently use the analytical value instead of zero — which may or may
    not be desired depending on whether the perturbation theory was
    fitted against the analytical value or the empirical value.
    Recommendation: add an `accuracy` error term reflecting
    |Jₙ_ellipsoid - Jₙ_empirical| at the fallback, or document the
    convention explicitly. **Flag for note.**
  - The `n < 2` case returns 0 silently (undocumented). Should probably
    assert or document.
```

---

## Card 10 — `ModelConfiguration::A`

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::ModelConfiguration::A
Location:               src/sgp4/model_selector.h:167-179
Mathematical statement: Aₙ,₀ = -Jₙ                          (zonal, m=0)
                        Aₙ,ₘ = 0                            (m ≠ 0, not yet)
                        In SGP4's normalized units (aₑ = 1 ER), this is
                        the [SR3] coefficient A_{nm} appearing in the
                        geopotential expansion.

THEORY
  Underlying theorem:   Geopotential spherical-harmonic expansion. For
                        a zonal harmonic, the SR3 perturbation coefficient
                        is A_{n,0} = -C_{n,0} · aₑⁿ, and the relation
                        between unnormalized Cₙ,₀ and Jₙ is
                            C_{n,0} = -Jₙ,
                        so A_{n,0} = Jₙ · aₑⁿ. In SGP4-normalized units
                        aₑ = 1 (Earth radii) ⇒ A_{n,0} = Jₙ. But the
                        comment on lines 156-158 says A_{n,0} = -Jₙ (with
                        opposite sign), citing SR3 directly.
  Primary reference:    - Kaula (1966) "Theory of Satellite Geodesy" §3.2.
                        - Hoots-Roehrich 1980 [SR3] §B: A₃₀ = -J₃·aₑ³,
                          A₄₀ = -J₄·aₑ⁴.
                        - The sign convention in the code (-Jₙ) matches
                          SR3 directly; the comment chain on lines
                          156-161 walks the double negation.
  Domain of validity:   n ≥ 2; only m = 0 currently implemented.

METHOD
  Method declared:      Closed-form: m≠0 → 0; m=0 → -Jₙ (with implicit
                        aₑⁿ = 1 in SGP4-normalized units).
  Method implemented:   `if (m != 0) return math::exact<T>(0);
                         return -Jn(n);`
  Match verdict:        ✓ matched to comment + SR3 sign convention.

ERROR BOUND
  Bound category:       Pass-through from Jn(n); negation does not change
                        any error category (REQ-EF-3 unary negation).
  Bound formula:        Same as Jn(n).
  Bound implemented:    `-Jn(n)` propagates errors via tracked_value's
                        unary minus.
  Bound verdict:        ✓ matched.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form unary op preserves bounds)
  AUD-EF applies:       AUD-EF-1
  AUD-MC applies:       n/a (no algebraic identity tested here)
  Verification test:    tests/test_sgp4/ — A(3,0) = -J₃; A(4,0) = -J₄;
                        A(n,m≠0) = 0.

NOTES
  ⚠ The comment on lines 158-161 says explicitly: "In SGP4's normalized
    units where aₑ = 1 Earth radius, A_{n,0} simplifies to -Jₙ. But this
    function carries aₑ explicitly so the result is correct in any unit
    system." However, the IMPLEMENTATION does NOT carry aₑ explicitly —
    it returns `-Jn(n)` directly, which is correct ONLY in normalized
    units. **Comment-code mismatch.** Either:
      (a) Update the comment to acknowledge that aₑ is implicitly 1.
      (b) Multiply by `pow(ellipsoid.a_normalized, n)` if a unit-agnostic
          A is wanted.
    Currently (a) is the safe interpretation since all callers are in
    SGP4-normalized units, but the documentation is misleading.
    **Flag for code or comment fix.**
  - The fallback at m ≠ 0 returns 0. The comment correctly notes that
    tesseral/sectoral support is future; document expectation now is 0.
```

---

## Cards 11-16 — `CustomBuilder` setters and `build()` (lines 191-214)

n/a — builder pattern setters return `*this`; `build()` aggregates factory
outputs into a `ModelConfiguration` struct. No new numeric formulas are
introduced. Theoretical basis lives in the factories (`make_zonals`,
`make_ellipsoid`, `make_astronomy`, `make_model_functions`) which are
audited below.

---

## Card 17 — `CustomBuilder::make_zonals` (lines 245-304)

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::CustomBuilder::make_zonals
Location:               src/sgp4/model_selector.h:245-304
Mathematical statement: Build a ZonalHarmonics<T> table populated with
                        empirical Jₙ values for the named gravity field.
                        Cases:
                          "wgs72_old" / "wgs72":   J₃, J₄ from [SR3] p.59
                          "wgs84_sgp4":            J₃, J₄ from Vallado (2006)
                          "wgs84_precise" / "grs80": J₃..J₉ from EGM2008
                          (default): WGS72 J₃, J₄
                        Conversion C̄ₙ₀ → Jₙ: Jₙ = -√(2n+1) · C̄ₙ₀.

THEORY
  Underlying theorem:   (1) Geophysical: odd zonals Jₙ break equatorial
                            symmetry; only an empirical gravity model
                            (EGM2008, EGM96) can supply them.
                        (2) Normalization conversion (lines 67-69, 236-237):
                            Jₙ = -√(2n+1) · C̄ₙ₀ — standard normalization
                            for fully-normalized spherical harmonics.
                        (3) The √(2n+1) factor is the Schmidt-quasi-norm
                            convention; with the additional Kronecker delta
                            δ_{m0} the full factor is √((2n+1)(n-m)!/(n+m)!)
                            which reduces to √(2n+1) for m=0.
  Primary reference:    - Hoots & Roehrich (1980) "Spacetrack Report No. 3"
                          page 59 COMMON/C1 (XJ2..XJ4, XJ3).
                        - Vallado (2006) `getgravconst` (WGS84 SGP4 values).
                        - NGA.STND.0036 (2014) Table 3.4 (EGM2008 C̄ₙ₀).
                        - IERS Conventions (2010) Table 6.2 (C̄₃₀, C̄₄₀).
                        - Heiskanen & Moritz (1967) §1-14 for the
                          normalization conversion identity.
  Domain of validity:   n = 3, 4, 5, …, 9 for "wgs84_precise"; n = 3, 4
                        for the others.

METHOD
  Method declared:      Hard-coded table of `TV::measured(value, uncertainty)`
                        entries, one switch branch per named gravity model.
  Method implemented:   Each branch calls `zh.add(n, TV::measured("…","…"))`
                        with literal strings for value and uncertainty.
                        The uncertainties for wgs84_precise are documented
                        as "formal EGM2008 errors propagated through √(2n+1)"
                        (line 288).
  Match verdict:        ✓ matched — direct table population.

ERROR BOUND
  Bound category:       measurement (the gravity-model uncertainty)
  Bound formula:        For each entry, the formal EGM2008 (or other)
                        coefficient error σ(C̄ₙ₀) is propagated through
                        the normalization conversion as
                        σ(Jₙ) = √(2n+1) · σ(C̄ₙ₀).
  Bound implemented:    The literal uncertainties in the string arguments
                        to `TV::measured` claim to be the propagated values.
                        Sample check (n=3):
                          C̄₃₀ = +0.9571612e-6 ± 4.9e-12 (per comment line 278)
                          √7 ≈ 2.64575
                          σ(J₃) = 4.9e-12 · 2.64575 ≈ 1.30e-11 ≈ 0.0000000000130
                          Code uses "0.0000000000130" ✓ matches.
                        Sample check (n=4):
                          C̄₄₀ = +0.5399659e-6 ± 4.7e-12
                          √9 = 3.0
                          σ(J₄) = 4.7e-12 · 3 ≈ 1.41e-11 ≈ 0.0000000000141
                          Code uses "0.0000000000141" ✓ matches.
                        Sample check (n=5):
                          σ(J₅) listed as 0.000000000010 = 1e-11
                          (No C̄₅₀ uncertainty given in comment; appears to
                          be a placeholder uniform 1e-11 for n ≥ 5.)
                        Sample check (wgs72 J₃):
                          Listed as "0.00000000001" = 1e-11. [SR3] does not
                          publish a formal uncertainty for J₃; this 1e-11
                          appears to be an analyst-supplied conservative
                          bound. **No primary-source justification.**
  Bound verdict:        ⚠ partial — for the IERS/EGM2008 entries n=3, 4,
                        the uncertainty arithmetic checks out against
                        the comment; for n ≥ 5 the uniform 1e-11 figure
                        is not derived in-file. For WGS72 J₃ and J₄, the
                        1e-11 figure appears to be an analyst placeholder
                        rather than a published formal error.
                        **Flag: trace n ≥ 5 EGM2008 σ(C̄ₙ₀) from primary
                        source to validate σ(Jₙ) = √(2n+1) · σ(C̄ₙ₀).**

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-2 (measured values declare measurement error)
  AUD-EF applies:       AUD-EF-2 (TV::measured constructs with declared σ)
  AUD-MC applies:       n/a
  Verification test:    tests/test_sgp4/ — make_zonals("wgs72").Jn(3).value
                        equals -2.53881e-6 ± 1e-11; etc.

NOTES
  - The comment block on lines 270-296 documents the C̄ₙ₀ values used to
    derive the listed Jₙ. For n=3, 4 the arithmetic is reproducible:
       -√7  · 0.9571612e-6 ≈ -2.532153e-6 ✓
       -√9  · 0.5399659e-6 ≈ -1.6198977e-6 ✓
       -√11 · 0.0686729e-6 ≈ -2.27733e-7  (code: -2.27735e-7, 4th sig
                                            digit differs — likely C̄₅₀
                                            given to higher precision)
       -√13 · (-0.1499537e-6) ≈ +5.40743e-7 (code: +5.40774e-7, 5th
                                              sig digit differs)
       -√15 · 0.0905123e-6 ≈ -3.50492e-7   (code: -3.50447e-7, 4th sig
                                              digit differs)
    The discrepancies at n=5,6,7 in the 4th-5th significant digit suggest
    either: (a) more decimal places of C̄ₙ₀ in the actual EGM2008 file than
    are reproduced in the comment, or (b) a small typo somewhere.
    **Flag: re-derive n=5..9 entries from raw EGM2008 file with full
    precision.** Either fix the comment to match the code, or fix the
    code to match the comment.
  - The "WGS72 J₃ = -2.53881e-6" and "wgs84_sgp4 J₃ = -2.53215306e-6"
    values do not agree (~0.3% relative difference). Both come from
    different historical fits; the matched-pair principle says use
    whichever matches the perturbation-theory fit. Documented adequately.
  - The default branch (no name match) returns WGS72 J₃/J₄. The selector
    `select("sgp4_standard")` requests "wgs72", so the default branch is
    only reachable via unrecognized gravity names; this is silent
    fallback. Consider warning or throwing.
```

---

## Card 18 — `CustomBuilder::make_ellipsoid` (lines 309-362)

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::CustomBuilder::make_ellipsoid
Location:               src/sgp4/model_selector.h:309-362
Mathematical statement: Construct an EquipotentialEllipsoid<T> from the
                        named gravity-field's four defining parameters
                        (a, J₂ or 1/f, GM, ω).

THEORY
  Underlying theorem:   The level ellipsoid is uniquely determined by any
                        four of its defining parameters (a, 1/f, GM, ω,
                        J₂); the other three are computed by inversion.
                        EquipotentialEllipsoid::from_J2 and the direct
                        constructor implement this inversion (see
                        equipotential_ellipsoid.h audit).
  Primary reference:    - Heiskanen & Moritz (1967) "Physical Geodesy"
                          Ch. 2, especially eqs. 2-73 to 2-97 for the
                          (a, J₂, GM, ω) → 1/f and J₂ₙ chain.
                        - Moritz (1980) "Geodetic Reference System 1980".
                        - NGA.STND.0036 (2014) Table 3.1 (WGS84 precise).
                        - [SR3] for WGS72 / WGS72_old values.
                        - Vallado (2006) `getgravconst` for WGS84-SGP4.
  Domain of validity:   The defining parameters provided by each branch.

METHOD
  Method declared:      Hard-coded TV::defined / TV::measured per-branch
                        construction. Default fallback to WGS72.
  Method implemented:   `auto a = TV::defined("…"); … return
                         EquipotentialEllipsoid<T>::from_J2(a, J2, GM, omega, tol);`
                        (or direct constructor for wgs84_precise with 1/f).
  Match verdict:        ✓ matched — direct dispatch to ellipsoid constructor.

ERROR BOUND
  Bound category:       Inherited from ellipsoid; per-branch:
                        - wgs72_old, wgs72, wgs84_sgp4, grs80: `TV::defined`
                          for everything → zero measurement error,
                          truncation/precision from from_J2 inversion.
                        - wgs84_precise: `TV::measured("398600.4418", "0.008")`
                          for GM → measurement error 0.008 km³/s² (= 8000
                          m³/s², matching NGA.STND.0036 quoted ±8e6 m³/s²).
                          Also `TV::measured("23.4393","0.0001")` for IAU
                          2006 obliquity. Everything else defined.
  Bound formula:        Pass-through from constructor; this function adds
                        no bound of its own.
  Bound implemented:    Pass-through.
  Bound verdict:        ✓ matched — bound is whatever the ellipsoid
                        constructor produces, plus the per-branch
                        measurement uncertainties.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (defined constants carry zero measurement),
                        REQ-EF-2 (measured values carry declared σ).
  AUD-EF applies:       AUD-EF-2, AUD-EF-3.
  AUD-MC applies:       n/a
  Verification test:    tests/test_geodesy/ + tests/test_sgp4/ — for each
                        named gravity, ellipsoid.a etc. match the expected
                        published values to the declared precision.

NOTES
  ⚠ The "wgs84_precise" branch is the only one carrying measurement
    uncertainty on GM and obliquity. For the "eleven" preset documentation
    (lines 561-617) claims "WGS84 with refined GM and full measurement
    uncertainty propagated through every derived constant" — this is
    correct only if the ellipsoid constructor propagates the GM
    uncertainty into all derived quantities (b, e², 1/f when J₂ is
    derived, etc.). Audit dependency: that propagation must be verified
    in `equipotential_ellipsoid.h` audit, not here.
  - The 0.008 km³/s² uncertainty on GM should be cross-checked against
    NGA.STND.0036 Table 3.1's quoted σ (~8e-9 km³/s² is typical; 0.008
    is much larger). **Flag: verify σ(GM) against primary source.**
    [Note: NGA gives ±5e5 m³/s² = 5e-4 km³/s² for GM_⊕ in some editions;
    if 0.008 km³/s² = 8e6 m³/s² differs, the comment "±8e6 m³/s²" at
    line 608 also disagrees with table.]
  - The "wgs84_precise" branch uses the direct constructor (with 1/f),
    not from_J2; this matches the "geometric precision" rationale on
    lines 567-569 of the eleven-preset comment.
  - The default fallback to WGS72 is silent.
```

---

## Card 19 — `CustomBuilder::make_astronomy` (lines 367-395)

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::CustomBuilder::make_astronomy
Location:               src/sgp4/model_selector.h:367-395
Mathematical statement: Construct an astronomy::FundamentalConstants<T>
                        for the named astronomy era.
                        - "sr3_1980": pristine sgp4_standard().
                        - "iau_2006" / "almanac_201[0,5]"/"almanac_2025":
                          sgp4_standard() with obliquity overridden to
                          IAU 2006 J2000.0 (23.4393° ± 0.0001°) and
                          solar_anomaly_period to sidereal year
                          (365.256363004 ± 1e-9 days).

THEORY
  Underlying theorem:   - SGP4 standard constants: [SR3] p. 59 DATA
                          statements (1970s NOA Almanac values).
                        - IAU 2006 obliquity at J2000.0: ε = 23.4393°
                          (84381.406 arcsec) from IAU 2006 precession
                          (Capitaine, Wallace, Chapront 2003 / IAU 2006).
                        - Sidereal year: 365.256363004 days from IERS
                          Conventions (2010) / VSOP87 mean.
  Primary reference:    - IAU 2006 Resolution B1 (precession model).
                        - IERS Conventions (2010) Table 1.1 (sidereal year).
                        - Hoots-Roehrich (1980) [SR3] p. 59.
  Domain of validity:   Epoch ≈ J2000.0 for the IAU 2006 numerical values.

METHOD
  Method declared:      Dispatch: name "sr3_1980" → unmodified standard;
                        IAU-2006-family names → standard with obliquity
                        and solar-period overrides.
  Method implemented:   Two branches matching declared method exactly.
                        The override block (lines 384-389) calls
                        `math::degrees_to_radians` and TV::measured.
  Match verdict:        ✓ matched.

ERROR BOUND
  Bound category:       (sr3_1980) inherited from sgp4_standard();
                        (iau_2006-family) measurement error 0.0001° on
                        obliquity (= ~1.7e-6 rad), and 1e-9 days on the
                        solar period.
  Bound formula:        For the obliquity override, the propagation
                        through degrees_to_radians multiplies by π/180.
                        σ(ε_rad) = (π/180) · 0.0001° ≈ 1.745e-6 rad.
  Bound implemented:    The TV is constructed with σ=0.0001 in degrees,
                        and degrees_to_radians applies the linear scale —
                        which itself is a closed-form unary op (REQ-EF-3).
                        So σ propagates correctly assuming degrees_to_radians
                        scales the measurement-error field. Verified in
                        `angles.h` audit.
  Bound verdict:        ? — DEPENDENT on degrees_to_radians correctness.
                        Card itself is ✓ matched; bound check defers to
                        angles.h audit.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-2 (measured TVs), REQ-EF-3 (closed-form
                        degrees_to_radians).
  AUD-EF applies:       AUD-EF-1, AUD-EF-3.
  AUD-MC applies:       n/a
  Verification test:    tests/test_sgp4/ — make_astronomy("iau_2006").
                        obliquity.value ≈ 0.40909 rad ± 1.7e-6.

NOTES
  ? The aliases "almanac_2010", "almanac_2015", "almanac_2020",
    "almanac_2025" all map to the same IAU 2006 J2000.0 obliquity. This
    is incorrect at the level of accuracy claimed by "eleven": the actual
    obliquity changes by ~47 arcsec/century, so the difference between
    2010 and 2025 is ~7 arcsec ≈ 3.4e-5 rad, larger than the declared
    σ = 1.7e-6 rad. The TODO on line 379 acknowledges:
       "populate from specific Almanac editions" — currently unhonored.
    **C-fail risk for almanac_20XX presets**: the claimed σ understates
    the time-evolution offset from J2000.0. If the propagation epoch
    is ≠ J2000.0, an additional accuracy-error term should reflect
    the precession-rate uncertainty over the elapsed interval.
    **Flag: implement obliquity drift or widen σ.**
  - The "sr3_1980" branch is a pure delegate; verdict for it inherits
    from astronomy::FundamentalConstants<T>::sgp4_standard() audit.
  - The fallback (unknown name) returns sgp4_standard() silently.
```

---

## Card 20 — `CustomBuilder::make_model_functions` (lines 400-474)

```
=== FORMULA AUDIT CARD ===
ID:                     model_selector::CustomBuilder::make_model_functions
Location:               src/sgp4/model_selector.h:400-474
Mathematical statement: Build a ModelFunctions<T>, possibly overriding the
                        kepler_solver lambda based on the kepler name:
                          "newton":      Newton-Raphson on the SGP4
                                         modified form
                                         f(x) = U - ayn·cos(x) + axN·sin(x) - x
                          "householder": Halley iteration on the same f
                                         (cubic convergence).
                          (default):     keep standard_sgp4()'s solver.
                        Perturbation override stub: only "brouwer_j2sq" is
                        recognized; it's the default so no override is
                        applied.

THEORY
  Underlying theorem:   - Newton-Raphson convergence theorem
                          (Kantorovich); quadratic on simple roots with
                          f' ≠ 0 in a neighborhood. For SGP4's modified
                          Kepler form f' = 1 - ayn·sin(x) - axN·cos(x)
                          which is positive for axN²+ayn² < 1 (a
                          consequence of e < 1).
                        - Halley's iteration: cubic-order root finder
                          using f, f', f''. Δ = 2·f·f' / (2·f'² - f·f'')
                          (Halley 1694; Conte & de Boor 1980 §3.7).
                          The label "householder" in this file is a
                          misnomer — Householder's method is a family
                          containing Halley as the d=1 case (next cubic-
                          order member); the code implements pure Halley.
                        - SGP4 modified Kepler form: derived in
                          Hoots-Roehrich (1980) [SR3] §C, where the
                          longitude-of-ascending-node frame yields
                          U = E - ayn·cos(E) + axN·sin(E)
                          for axN = e·cos(ω), ayn = e·sin(ω).
                          (Note: f below is the residual `U - … - x`,
                          NOT the orbital true anomaly.)
  Primary reference:    - Halley (1694) Phil. Trans. R. Soc. 18, 136.
                        - Conte & de Boor (1980) §3.7.
                        - Hoots-Roehrich (1980) [SR3] §C.
                        - Battin (1999) §5.3 for Kepler-equation context.
  Domain of validity:   axN, ayn from SGP4 internal state; convergence
                        guaranteed for elliptic orbits (e < 1).

METHOD
  Method declared:      - "newton": iterate x ← x + (f / f'); cap 30
                          iterations; tolerance break.
                        - "householder" (in fact Halley): iterate
                          x ← x + 2·f·f' / (2·f'² - f·f'');
                          cap 15 iterations; tolerance break.
                        - Bound at convergence: `precision += |delta|`.
  Method implemented:   The lambdas at lines 410-430 and 438-461
                        implement exactly the declared methods, including
                        the |delta| < tolerance break and the precision
                        increment.
  Match verdict:        ✓ matched — the code IS Newton or Halley as
                        declared (not a continued fraction, not a
                        Padé approximant, not a series). The "Halley"
                        formula 2ff'/(2f'² - ff'') is the textbook one;
                        the label "householder" in the dispatch name
                        is a vocabulary issue, not an implementation
                        defect. See NOTES.

ERROR BOUND
  Bound category:       precision (convergence residual)
  Bound formula:        For any iterative root-finder converging at order
                        ≥ 2, |x_∞ - x_k| ≤ |Δ_k| / (1 - q^p), with q < 1
                        and p the order. In practice |Δ_k| itself is a
                        conservative bound (REQ-EF-5; see kepler.h audit
                        card §6 of framework doc). This is the
                        Kantorovich / Ostrowski bound.
  Bound implemented:    `x.errors.precision = x.errors.precision +
                         abs(delta.value);` (line 425 and line 456).
                        Added at the iteration that crosses tolerance.
  Bound verdict:        ✓ matched — final |Δ| as the precision bound,
                        per REQ-EF-5.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (iterative residual → precision).
  AUD-EF applies:       AUD-EF-4 (iterative algorithms add residual).
  AUD-MC applies:       n/a
  Verification test:    tests/test_sgp4/ — with both "newton" and
                        "householder", verify Kepler equation residual
                        < tolerance and that precision ≥ |true error|.

NOTES
  ⚠ The dispatch name "householder" is a misnomer. The implementation
    (lines 438-461) is Halley's method:
       Δ = 2·f·f' / (2·f'² − f·f'')
    Householder's method of order d uses higher derivatives; Halley
    is the d=1 member of the Householder family. The comment on
    line 433 correctly labels it as "Halley (cubic convergence)" but
    the user-facing string says "householder". The "eleven" preset
    documentation on line 583 also misnames it: "Householder quartic
    convergence — triples then quadruples correct digits". Halley is
    cubic, not quartic; the doubling/tripling claim is incorrect.
    **Flag: rename dispatch string to "halley" OR implement true
    Householder-d=2 (quartic) using f'''. Update preset comment.**
  - The fallback iteration cap (30 for Newton, 15 for Halley) does NOT
    add an error term if the cap is reached without converging — the
    final |delta| is returned unannotated except for the standard
    precision-bound path. For pathological inputs (numerically loose
    e ≈ 1 cases), the cap could be hit; the silent return without
    flagging is a potential issue. **Note: see also kepler.h.**
  - "Perturbation theory selection" is currently no-op for the only
    branch ("brouwer_j2sq"). Future "brouwer_j2cu", "brouwer_full"
    branches will need their own audit cards.
  - The drag parameter is ignored (line 403 `const std::string& /*drag*/`).
    Future drag selection will require its own dispatch and audit.
```

---

## Cards 21-22 — `ModelSelector::select` and `::custom` (lines 510-627)

n/a — pure dispatch into `CustomBuilder`. Each preset (`"sgp4_standard"`,
`"sgp4_wgs72_old"`, `"sgp4_wgs84"`, `"modern_2020"`, `"research_full"`,
`"eleven"`) sets builder fields then calls `build()`. No new numeric
formulas are introduced beyond what `CustomBuilder::build()` and its
underlying factories produce.

One documentation discrepancy worth flagging:

- The "eleven" preset's docstring (lines 561-617) advertises:
  "Householder quartic convergence — triples then quadruples
   correct digits per iteration"
  but `make_model_functions` for `kep == "householder"` implements
  Halley (cubic). See Card 20 NOTES.
- The same docstring on line 608 says "±8e6 m³/s² measurement uncertainty"
  for GM, but `make_ellipsoid` for wgs84_precise uses
  `TV::measured("398600.4418", "0.008")` = ±0.008 km³/s² = ±8e3 m³/s².
  The docstring overstates by 10³.

---

## File-level verdict

- **A. Error wiring (EF)**: ✓ where applicable — math-relevant members
  (`Jn`, `A`, kepler lambdas) propagate TrackedValue<T> through standard
  closed-form / pass-through paths. The Kepler lambdas correctly add
  `|Δ|` to `precision` per REQ-EF-5.
- **B. Algebra axioms (MC)**: n/a — this file is configuration / dispatch.
- **C. Theoretical basis (TBA)**:
  - Container utilities (`add`, `has`, `max_degree`, setters, `build`,
    `select`, `custom`): n/a — no formulas.
  - `ZonalHarmonics::Jn`, `J3`, `J5`, `J7`, `J9`: ✓ pass-through lookup.
  - `ModelConfiguration::Jn`: ⚠ dispatch logic correct; high-degree
    ellipsoid-fallback convention should be documented or carry an
    accuracy term.
  - `ModelConfiguration::A`: ⚠ correct in SGP4-normalized units; the
    comment claims unit-agnosticism but the implementation does not
    carry aₑⁿ. Comment-code mismatch.
  - `make_zonals`: ⚠ wgs72_old / wgs84_precise n=3,4 σ values
    reproduce comment-derived numbers; n=5..9 entries have 4th-5th
    significant-digit mismatches vs. the C̄ₙ₀ values quoted in the
    comment. Trace from primary EGM2008 file required.
  - `make_ellipsoid`: ⚠ wgs84_precise GM uncertainty 0.008 km³/s²
    inconsistent with the "eleven" preset's claim of "±8e6 m³/s²"
    (which would be 8 km³/s², 10³ larger). Verify against
    NGA.STND.0036.
  - `make_astronomy`: ? all four `almanac_201[0,5,20,25]` aliases
    return the same J2000.0 obliquity. Time-evolution drift (~47
    arcsec/century) is larger than the declared σ ~ 0.0001°.
    Either widen σ or implement epoch-specific obliquity per TODO.
  - `make_model_functions`: ✓ Newton and Halley correctly implemented;
    bounds are |Δ| per REQ-EF-5. ⚠ Dispatch string "householder" is
    misnomer for Halley; "eleven" docstring incorrectly calls it
    "quartic". Rename or upgrade to true quartic.

**File verdict: PASS with notes** — no theory ↔ method mismatches in the
implementations themselves. All flagged items are either (a) documentation
inconsistencies that misrepresent the implementation, or (b) σ-derivation
traceability gaps that should be closed against primary sources.

---

## Action items (priority order)

1. **[Doc fix]** Rename Kepler dispatch string from `"householder"` to
   `"halley"` and update the "eleven" preset docstring to say "cubic
   convergence — Halley iteration". Or, alternatively, implement a true
   Householder-d=2 (quartic) iteration using f''' and keep the name.
2. **[Doc fix]** Reconcile the "eleven" preset's "±8e6 m³/s²" GM σ claim
   with the actual `TV::measured("398600.4418", "0.008")` = ±8e3 m³/s².
3. **[Data trace]** For `make_zonals("wgs84_precise")` entries n=5..9,
   re-derive Jₙ from the full-precision EGM2008 C̄ₙ₀ values and
   reconcile with the comment block.
4. **[Convention]** Decide whether high-degree even-zonal fallback in
   `ModelConfiguration::Jn` (ellipsoid analytical formula) should add
   an `accuracy` term reflecting |J_ellipsoid - J_empirical|, or be
   documented as the explicit policy.
5. **[Future]** When `almanac_20XX` presets are populated per-edition
   (TODO line 379), implement epoch-specific obliquity drift or widen
   σ to encompass the J2000.0-to-epoch offset.
6. **[Future]** Implement `brouwer_j2cu` / `brouwer_full` perturbation
   branches and write their audit cards.
7. **[Future]** Wire the `drag` parameter through `make_model_functions`
   and audit the resulting drag-model lambdas.
