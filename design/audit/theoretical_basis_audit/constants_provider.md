# Theoretical Basis Audit — `src/constants/constants_provider.h`

**Status**: drafted 2026-05-13.
**Scope**: 4 numeric-producing entry points — the bundle constructor and three named-convention factories (`wgs84`, `wgs72`, `grs80`).
**Framework**: see `design/audit/theoretical_basis_audit.md` §1 (per-formula card schema), §5 (worked example).

The file produces no "computed" numbers of its own; every value is a literal from a published convention wrapped in `TrackedValue<T>::defined(...)` or `TrackedValue<T>::measured(...)` and forwarded into a `geodesy::EquipotentialEllipsoid<T>` constructor (direct or `from_J2`). The audit therefore focuses on **convention fidelity**: does the literal match the cited document to its full published precision, is the uncertainty class (`defined` vs `measured`) the correct one per the convention's text, and does the construction path (`(a, 1/f, GM, ω, …)` vs `from_J2(a, J2, GM, ω, …)`) match the convention's *defining* parameter set?

## 1. `ConstantsProvider::ConstantsProvider(const EquipotentialEllipsoid<T>&)`

```
=== FORMULA AUDIT CARD ===
ID:                     constants_provider::ctor
Location:               src/constants/constants_provider.h:39-40
Mathematical statement: Bundle constructor — copy-store a pre-built
                        EquipotentialEllipsoid<T> into the `earth` member.
                        No numeric formula; pure aggregation.

THEORY
  Underlying theorem:   None. This is a dependency-injection wrapper; the
                        only mathematical content is that copying a
                        TrackedValue<T> preserves its value and its
                        per-category error bounds (REQ-EF-2 / REQ-EF-3).
  Primary reference:    REQ-CP-1 (single constants bundle), REQ-CP-4..8
                        (sub-bundle structure). The math content is
                        delegated to EquipotentialEllipsoid<T>'s own audit
                        card (forthcoming).
  Domain of validity:   Any well-formed EquipotentialEllipsoid<T>; the
                        ellipsoid's own constructor enforces a > 0, etc.

METHOD
  Method declared:      Member-wise copy construction.
  Method implemented:   `: earth(e) {}` — single member-initializer copy.
  Match verdict:        ✓ matched — pure copy, no transformation.

ERROR BOUND
  Bound category:       n/a (no new uncertainty introduced).
  Bound formula:        Identity — output errors equal input errors per
                        category, by definition of TrackedValue copy.
  Bound implemented:    Inherited via the copy-construction of `earth`,
                        which itself member-wise-copies each
                        TrackedValue<T>.
  Bound verdict:        ✓ matched — no error category is touched here.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-2 (per-category error preservation under
                        copy), REQ-EF-3 (closed-form propagation — trivial
                        identity here).
  AUD-EF applies:       AUD-EF-1 (TrackedValue plumbing); AUD-EF-7
                        (composition preserves bounds).
  AUD-MC applies:       n/a (not an algebra operation).
  Verification test:    tests/test_constants_provider/ (or equivalent) —
                        verify that ConstantsProvider(e).earth has the
                        same numeric value and per-category errors as `e`.

NOTES
  - The constructor is `explicit`, preventing accidental implicit
    construction from a stray ellipsoid; this is the AUD-CC-7 stylistic
    contract.
  - There is no factory-less path that mutates internal state; all
    populated bundles come from `wgs84` / `wgs72` / `grs80` factories
    below, or from an explicit `ConstantsProvider(ellipsoid)` call where
    the caller has already audited the ellipsoid.
```

## 2. `ConstantsProvider::wgs84(series_tolerance)`

```
=== FORMULA AUDIT CARD ===
ID:                     constants_provider::wgs84
Location:               src/constants/constants_provider.h:51-60
Mathematical statement: Construct the WGS 84 reference ellipsoid bundle from
                        the four defining parameters
                          a     = 6 378 137.0 m              (definitional)
                          1/f   = 298.257 223 563             (definitional)
                          GM    = 3.986 004 418 × 10^14 m³/s² (measured,
                                                              σ ≈ 8 × 10^5)
                          ω     = 7.292 115 1467 × 10^-5 rad/s (definitional)
                        and forward them to the EquipotentialEllipsoid<T>
                        primary constructor (the (a, 1/f, GM, ω, tol) path).

THEORY
  Underlying theorem:   WGS 84 convention as published. WGS 84 fixes
                        (a, 1/f, GM, ω) as its defining set; all other
                        ellipsoidal quantities (b, e², J₂, n, …) are
                        derived from these via the closed-form identities
                        of an equipotential rotating ellipsoid (Somigliana,
                        Heiskanen & Moritz). GM is the only one of the
                        four that is *measured* (it is fixed to a value
                        consistent with current observations and carries
                        an associated 1σ uncertainty); a, 1/f, ω are
                        adopted by definition for the purpose of the
                        reference system.
  Primary reference:    NIMA TR 8350.2, 3rd Edition (Amendment 1, 2000),
                        with 2014 Addendum. Table 3.1 ("WGS 84 Defining
                        Parameters") and §3.2:
                          a      = 6 378 137.0 m
                          f      = 1 / 298.257 223 563
                          GM     = 3 986 004.418 × 10^8 m³/s²
                                 = 3.986 004 418 × 10^14 m³/s²
                          ω_E    = 7 292 115.1467 × 10^-11 rad/s
                                 = 7.292 115 1467 × 10^-5 rad/s
                        and §3.4 / Table 3.3 for the GM uncertainty
                          σ(GM)  = 8 × 10^5 m³/s² (1σ, from satellite
                                                   laser ranging fits).
                        (NIMA TR 8350.2 §3.4.1 explicitly classifies GM
                        as "measured" / "adopted from observations" and
                        a, 1/f, ω as "defining" / "definitional".)
  Domain of validity:   Earth — rotating equipotential ellipsoid model.
                        Forwarded to EquipotentialEllipsoid<T> which
                        requires a > 0, 1/f > 1, GM > 0, ω real.

METHOD
  Method declared:      Direct literal binding of the four NIMA-published
                        defining values into TrackedValue<T> wrappers,
                        with `measured` used for GM (carries σ) and
                        `defined` used for a, 1/f, ω (zero measurement
                        uncertainty by convention). Then construction
                        of EquipotentialEllipsoid<T> via the
                        (a, 1/f, GM, ω, tol) constructor — this is the
                        convention's defining path, so no derivation
                        mismatch.
  Method implemented:   Lines 53-58:
                          a     = TrackedValue<T>::defined("6378137.0");
                          inv_f = TrackedValue<T>::defined("298.257223563");
                          GM    = TrackedValue<T>::measured(
                                      "3.986004418e14", "8e5");
                          omega = TrackedValue<T>::defined(
                                      "7.2921151467e-5");
                          return ConstantsProvider(
                              EquipotentialEllipsoid<T>(
                                  a, inv_f, GM, omega, series_tolerance));
                        All four literals match the NIMA values to their
                        full published precision.
  Match verdict:        ✓ matched — literals and uncertainty classes
                        match NIMA TR 8350.2 Tables 3.1 / 3.3 exactly,
                        and the construction path (a, 1/f, GM, ω) is the
                        convention's defining set.

ERROR BOUND
  Bound category:       For GM:    measurement (σ = 8 × 10^5 m³/s²,
                                   from NIMA TR 8350.2 §3.4 / Table 3.3).
                        For a, 1/f, ω: none — `defined` means the value
                                        is exact by convention (precision
                                        bound = 0; measurement bound = 0;
                                        accuracy bound = 0 at the
                                        TrackedValue construction site).
  Bound formula:        TrackedValue<T>::defined(s): all error categories
                                                     initialized to 0.
                        TrackedValue<T>::measured(v, σ): errors.measurement
                                                         initialized to σ.
                        These are constructor-defined contracts of
                        TrackedValue<T>, not derived bounds — they encode
                        the convention's own published uncertainty model.
  Bound implemented:    `defined(...)` and `measured(...)` factories of
                        TrackedValue<T> as currently implemented at
                        src/math/tracked_value.h:135 and :153. The σ
                        literal "8e5" is the NIMA-published 1σ value for
                        GM.
  Bound verdict:        ✓ matched — the σ literal equals NIMA TR 8350.2
                        Table 3.3 GM uncertainty exactly. Defining
                        parameters carry zero uncertainty as required by
                        the convention's text.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (initial bound categorization at the
                        point of construction); REQ-EF-2 (per-category
                        preservation through forwarding into the
                        ellipsoid constructor).
  AUD-EF applies:       AUD-EF-1 (TrackedValue construction); AUD-EF-7
                        (composition preservation).
  AUD-MC applies:       n/a (no algebra performed here; the ellipsoid's
                        derived quantities are audited in its own card).
  Verification test:    tests/test_constants_provider — verify
                        `wgs84(tol).earth.a().value() == 6378137.0`,
                        `wgs84(tol).earth.GM().errors.measurement ≈ 8e5`,
                        and that `defined` values report 0 measurement
                        error.

NOTES
  - **Uncertainty class of GM**: per the user's task note and NIMA TR
    8350.2 §3.4 / Table 3.3, GM is *measured* (it is fixed to be
    consistent with observations and carries a 1σ uncertainty of
    8 × 10^5 m³/s²). The implementation correctly uses
    `TrackedValue<T>::measured(...)` here, distinguishing it from the
    three defining parameters that use `defined(...)`.
  - **σ source**: NIMA TR 8350.2 Table 3.3 ("WGS 84 Parameter Values
    and Sigmas") gives σ(GM) = 8 × 10^5 m³/s² (1σ). The 2014 Addendum
    does not change this value; it adds the EGM2008 / EGM-aware GM_e
    convention (GM_e = 3.986 004 415 × 10^14 with σ ≈ 8 × 10^5) for
    Earth-orbiting satellite work, but TR 8350.2's GM remains the
    reference value at the resolution of this constants provider.
  - **1/f precision**: "298.257223563" matches NIMA TR 8350.2 to all
    9 significant figures of the published flattening. The reciprocal
    flattening was *adopted* as definitional in WGS 84 (specifically
    in the 1984 redefinition relative to WGS 72); the literal is exact.
  - **ω precision**: "7.2921151467e-5" matches the IAU 1976 nominal
    Earth rotation rate as adopted into WGS 84. Definitional.
  - **a precision**: "6378137.0" is the WGS 84 semi-major-axis
    definitional value (exact integer + .0).
  - **Construction path**: this routes through the EquipotentialEllipsoid
    primary constructor with (a, 1/f, GM, ω, tol) — the same parameter
    set NIMA TR 8350.2 publishes as defining. No `from_J2`-style
    inversion is needed (or used), so there is no inversion-derived
    rounding to budget.
```

## 3. `ConstantsProvider::wgs72(series_tolerance)`

```
=== FORMULA AUDIT CARD ===
ID:                     constants_provider::wgs72
Location:               src/constants/constants_provider.h:73-82
Mathematical statement: Construct the WGS 72 reference ellipsoid bundle
                        from the convention's four defining parameters
                          a     = 6 378 135.0 m              (definitional)
                          GM    = 3.986 008 × 10^14 m³/s²     (definitional)
                          J₂    = 1.082 616 × 10^-3            (definitional)
                          ω     = 7.292 115 1467 × 10^-5 rad/s
                        and forward them to the `from_J2` factory of
                        EquipotentialEllipsoid<T>, which inverts the
                        Heiskanen-Moritz J₂ ↔ (a, f, GM, ω) closed-form
                        relation to recover 1/f from the given J₂.

THEORY
  Underlying theorem:   WGS 72 convention as published. WGS 72 differs
                        from WGS 84 in which parameters are *defining*:
                        WGS 72 publishes (a, GM, J₂, ω) as definitional —
                        not (a, 1/f, GM, ω). All four are exact by
                        adoption (no σ); J₂ is treated as defining rather
                        than measured because WGS 72 was constructed
                        before the precision-GM era and adopted J₂ as a
                        fixed reference quantity. The flattening f is
                        then *derived* from J₂ via the equipotential
                        closed form
                          J₂ = (2/3)·(f − m/2 + …)
                        (Heiskanen & Moritz 1967 §2.10 series; the
                        EquipotentialEllipsoid::from_J2 factory implements
                        the exact closed-form inversion to caller-
                        specified tolerance — see its own audit card
                        for the inversion bound).
  Primary reference:    Hoots & Roehrich 1980, "Spacetrack Report #3:
                        Models for Propagation of NORAD Element Sets",
                        Appendix B "WGS-72 Constants", which adopts the
                        original WGS 72 publication (DMA TR 8350.2-A
                        / "World Geodetic System 1972", 1974):
                          a_E   = 6 378 135 m
                          GM    = 3.986 008 × 10^14 m³/s²
                          J₂    = 1.082 616 × 10^-3
                          ω_E   = 7.292 115 1467 × 10^-5 rad/s
                        SR3 §B explicitly lists these four as the WGS 72
                        constants used by SGP4/SDP4. ω is carried over
                        from the same IAU adoption that WGS 84 uses, so
                        the literal is identical to wgs84's.
                        Background on the equipotential J₂ ↔ f relation:
                        Heiskanen & Moritz, "Physical Geodesy" (1967)
                        Chapter 2 §2.10 ("Series Expansions for the
                        Normal Gravity Field"); Moritz (1980) "Geodetic
                        Reference System 1980" Bull. Géod. 54 — same
                        machinery, applied to GRS 80.
  Domain of validity:   Earth — rotating equipotential ellipsoid model.
                        `from_J2` requires the inversion to converge to
                        the caller's tolerance.

METHOD
  Method declared:      Direct literal binding of the four WGS 72 defining
                        values as `defined` TrackedValue<T> (no σ — all
                        four are conventional adoptions per WGS 72 / SR3
                        Appendix B). Then EquipotentialEllipsoid<T>::
                        from_J2 inverts J₂ ↔ 1/f, yielding an ellipsoid
                        whose `J2()` accessor reproduces the input J₂
                        exactly (modulo `series_tolerance`).
  Method implemented:   Lines 75-81:
                          a     = TrackedValue<T>::defined("6378135.0");
                          GM    = TrackedValue<T>::defined("3.986008e14");
                          J2    = TrackedValue<T>::defined("1.082616e-3");
                          omega = TrackedValue<T>::defined("7.2921151467e-5");
                          return ConstantsProvider(
                              EquipotentialEllipsoid<T>::from_J2(
                                  a, J2, GM, omega, series_tolerance));
                        All four literals match SR3 Appendix B exactly.
  Match verdict:        ✓ matched — literals, uncertainty class
                        (all `defined`), and construction path
                        (`from_J2`, because J₂ is the convention's
                        defining-not-derived quantity) all match the
                        convention. The header comment explicitly
                        documents WHY `from_J2` is used:
                          "Uses the `from_J2` initialization path because
                          J₂ is defining; deriving 1/f from J₂ avoids
                          the 4.8 × 10^-9 mismatch documented in
                          equipotential_ellipsoid.h 'Key Finding:
                          WGS72 J2'."
                        — this is the C-criterion compliance: method
                        (J₂-inversion) matches theory (J₂ is the
                        defining parameter).

ERROR BOUND
  Bound category:       For a, GM, J₂, ω: none at the constants-provider
                                          layer — all are `defined`, so
                                          measurement = 0 and precision = 0
                                          at construction.
                        The from_J2 inversion introduces a *precision*
                        bound (the Newton/series tolerance for inverting
                        J₂ → 1/f) which is budgeted **inside**
                        EquipotentialEllipsoid::from_J2; that bound is
                        audited in the ellipsoid's own card, not here.
                        At this layer, the constants_provider passes
                        the caller's `series_tolerance` through and is
                        not responsible for the inversion budget itself.
  Bound formula:        TrackedValue<T>::defined(s) → all errors 0.
                        Inversion bound: deferred to
                        EquipotentialEllipsoid<T>::from_J2 audit (the
                        inversion is a Newton/series iteration whose
                        bound is the final-correction magnitude per
                        REQ-EF-5, mediated by `series_tolerance`).
  Bound implemented:    `defined(...)` constructor for all four; tolerance
                        forwarded to `from_J2`. The ellipsoid factory is
                        responsible for adding the inversion residual to
                        `inv_f.errors.precision` (or its derived
                        equivalent). The constants_provider itself does
                        not add to any bound here — it is correctly
                        passive at this layer.
  Bound verdict:        ✓ matched — at the constants_provider layer, no
                        bound is computed; all bound computation is
                        delegated to the ellipsoid factory, which is the
                        correct separation of concerns. Defining params
                        carry zero uncertainty per the convention.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1 (initial categorization); REQ-EF-3
                        (closed-form propagation into the ellipsoid).
                        The inversion bound at the ellipsoid layer is
                        REQ-EF-5 (iterative tolerance), audited there.
  AUD-EF applies:       AUD-EF-1 (TrackedValue plumbing).
  AUD-MC applies:       n/a.
  Verification test:    tests/test_constants_provider — verify
                        `wgs72(tol).earth.a().value() == 6378135.0`,
                        `wgs72(tol).earth.GM().value() == 3.986008e14`,
                        `wgs72(tol).earth.J2().value() == 1.082616e-3`
                        within `tol`, and that all four `defined` values
                        report zero measurement error. Crucially, verify
                        that `from_J2` reproduces the input J₂ to within
                        `tol` (round-trip test); this is the audit
                        anchor for the "4.8 × 10^-9 mismatch avoided"
                        claim in the header comment.

NOTES
  - **No GM uncertainty**: unlike WGS 84, WGS 72 / SR3 Appendix B
    publishes GM as a *definitional* value (no σ). This reflects the
    pre-precision-laser-ranging era of WGS 72's adoption. The
    implementation correctly uses `defined` for all four parameters.
    A purist might argue an "epistemic" uncertainty band is
    appropriate, but the convention's published text treats GM as
    exact-by-adoption for WGS 72, and the audit honors the convention.
  - **Why `from_J2` and not the primary constructor**: the header
    comment on lines 71-72 is the operative audit anchor. WGS 72's
    *defining* set is (a, GM, J₂, ω) — NOT (a, 1/f, GM, ω). If we
    forward-fed an externally-computed 1/f into the primary
    constructor and then asked for `J2()`, the closed-form J₂
    derivation would yield a value ~4.8 × 10^-9 off from the
    convention's 1.082 616 × 10^-3. By inverting from J₂ instead,
    we ensure `earth.J2().value() == 1.082616e-3` exactly (to
    `series_tolerance`), which is required for SGP4/SDP4 to match
    SR3 numerical references at the documented precision.
  - **ω literal**: identical to wgs84's "7.2921151467e-5". SR3 §B
    adopts the IAU 1976 rotation rate, the same value WGS 84 later
    adopted. The header comment on line 68 notes this.
  - **ε (mean motion of the Greenwich meridian) and tu (time unit)**:
    not exposed here. SR3 §B lists additional derived quantities
    (k_e, time conversion); these are computed downstream from
    (a, GM) inside the propagator and are audited under
    `gravity_central.h` / `propagator.h`, not here.
```

## 4. `ConstantsProvider::grs80(series_tolerance)`

```
=== FORMULA AUDIT CARD ===
ID:                     constants_provider::grs80
Location:               src/constants/constants_provider.h:91-100
Mathematical statement: Construct the GRS 80 reference ellipsoid bundle
                        from its four defining parameters
                          a     = 6 378 137.0 m            (definitional)
                          GM    = 3.986 005 × 10^14 m³/s²  (definitional)
                          J₂    = 1.082 63 × 10^-3          (definitional)
                          ω     = 7.292 115 × 10^-5 rad/s   (definitional)
                        forwarded to the `from_J2` factory of
                        EquipotentialEllipsoid<T>, which inverts J₂ ↔ 1/f
                        via the equipotential closed-form relation.

THEORY
  Underlying theorem:   GRS 80 convention as published. GRS 80, like
                        WGS 72, adopts (a, GM, J₂, ω) as its defining
                        set — explicitly NOT (a, 1/f, GM, ω). Moritz's
                        Bull. Géod. 54 paper (1980) lists J₂ (not the
                        dynamical form factor of WGS 84) as the second
                        defining parameter, alongside a, GM, ω. All four
                        are adopted by international convention; the
                        flattening, eccentricity, second flattening, etc.
                        are all derived from these via the same
                        Heiskanen-Moritz equipotential machinery.
  Primary reference:    H. Moritz, "Geodetic Reference System 1980",
                        Bulletin Géodésique 54 (1980), pp. 395-405,
                        Table 1 ("Defining Parameters"):
                          a       = 6 378 137 m
                          GM      = 3.986 005 × 10^14 m³/s²
                          J₂      = 1.082 63 × 10^-3
                          ω       = 7.292 115 × 10^-5 rad/s
                        Reaffirmed in Moritz, "Geodetic Reference System
                        1980", J. Geodesy 74 (2000), pp. 128-133
                        (republished after editorial review). The
                        defining-parameter set is unchanged from the
                        1980 IAG adoption.
  Domain of validity:   Earth — rotating equipotential ellipsoid model.
                        `from_J2` inversion must converge to caller's
                        tolerance.

METHOD
  Method declared:      Direct literal binding of the four GRS 80
                        defining values as `defined` TrackedValue<T>,
                        then `EquipotentialEllipsoid<T>::from_J2`
                        inversion to recover 1/f from J₂ (same theoretical
                        machinery as wgs72, applied to GRS 80's
                        numerical adoption).
  Method implemented:   Lines 93-99:
                          a     = TrackedValue<T>::defined("6378137.0");
                          GM    = TrackedValue<T>::defined("3.986005e14");
                          J2    = TrackedValue<T>::defined("1.08263e-3");
                          omega = TrackedValue<T>::defined("7.292115e-5");
                          return ConstantsProvider(
                              EquipotentialEllipsoid<T>::from_J2(
                                  a, J2, GM, omega, series_tolerance));
                        All four literals match Moritz 1980 Table 1
                        exactly.
  Match verdict:        ✓ matched — literals, uncertainty class
                        (all `defined`), and construction path (`from_J2`,
                        because J₂ is GRS 80's defining quantity) all
                        match Moritz 1980. The shared use of `from_J2`
                        with wgs72 is principled: both conventions
                        adopt J₂ as defining, and the same closed-form
                        inversion theory governs both.

ERROR BOUND
  Bound category:       For a, GM, J₂, ω: none at the constants-provider
                                          layer — all `defined`, so all
                                          per-category errors are 0 at
                                          construction.
                        The from_J2 inversion residual is budgeted inside
                        EquipotentialEllipsoid<T>::from_J2 (REQ-EF-5,
                        iterative tolerance) — same as wgs72.
  Bound formula:        TrackedValue<T>::defined(s) → all errors 0.
                        Inversion bound: deferred to the ellipsoid
                        factory's audit.
  Bound implemented:    `defined(...)` for all four parameters; tolerance
                        passed through to `from_J2`. The constants_provider
                        layer adds nothing.
  Bound verdict:        ✓ matched — defining parameters are exact;
                        inversion budget is delegated to the ellipsoid
                        factory.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-1, REQ-EF-3; iterative-tolerance bound
                        REQ-EF-5 at the ellipsoid layer (audited there).
  AUD-EF applies:       AUD-EF-1.
  AUD-MC applies:       n/a.
  Verification test:    tests/test_constants_provider — verify literals,
                        verify `defined` zero-error contract, verify
                        `from_J2` round-trip reproduces input J₂ to
                        within `series_tolerance`.

NOTES
  - **ω literal precision**: "7.292115e-5" (seven significant figures)
    matches Moritz 1980 Table 1 exactly. This is *less* precise than
    the WGS 72 / WGS 84 ω literal "7.2921151467e-5" (eleven sig figs).
    The two conventions are using ω values from different epochs of IAU
    adoption: GRS 80 used the rounded 1980-era value, WGS 84 carries
    the more precise IAU 1976 value. Both are convention-exact ("defined")
    *within their own conventions*; the difference (~4.67 × 10^-13)
    is not a bug — it is the documented divergence of the two
    conventions' rotation-rate adoptions. Auditor confirms the literal
    matches Moritz 1980's published value to its full sig-fig count.
  - **a literal**: "6378137.0" — identical to WGS 84's. GRS 80 and
    WGS 84 share the same defining semi-major axis. Definitional.
  - **GM literal**: "3.986005e14" — six significant figures, matching
    Moritz 1980. (WGS 84 publishes a more precise GM with σ; GRS 80
    publishes its GM as a defining value at lower precision. Both are
    self-consistent within their conventions.)
  - **J₂ literal**: "1.08263e-3" — six significant figures, matching
    Moritz 1980. (WGS 72 publishes "1.082616e-3" — same precision class,
    slightly different value. The two J₂ values are different
    convention adoptions, not the same number rounded differently.)
  - **No σ for GM** in GRS 80: like WGS 72 and unlike WGS 84, GRS 80
    publishes GM as a definitional adoption without uncertainty.
    Implementation correctly uses `defined`.
  - **Construction path matches wgs72's reasoning**: GRS 80's defining
    set is (a, GM, J₂, ω) — J₂ is the second defining parameter. Using
    `from_J2` ensures `earth.J2().value()` reproduces the input J₂
    exactly (to tolerance), matching the convention's adopted value.
```

## 5. File-level verdict

**Dimensions:**

- **A. Error wiring**: ✓. Every literal is wrapped in a `TrackedValue<T>::defined(...)` or `TrackedValue<T>::measured(...)` factory before reaching the ellipsoid constructor; per-category bounds are initialized correctly at the boundary between "literal in source" and "TrackedValue in the propagator". GM in WGS 84 is the only `measured` value and correctly carries σ = 8 × 10^5 m³/s².

- **B. Algebra axioms**: n/a at this layer (no algebra is performed). The algebraic operations on these constants live in `EquipotentialEllipsoid` and the force lambdas; their audits are separate.

- **C. Theoretical basis**:
  - 4.1 ctor: ✓ identity copy. **PASS**.
  - 4.2 `wgs84`: ✓ literals, uncertainty classes, and construction path all match NIMA TR 8350.2 §3.2 / Table 3.1 / Table 3.3. **PASS**.
  - 4.3 `wgs72`: ✓ literals match SR3 §B / WGS 72 publication; correct use of `from_J2` per the convention's defining-J₂ adoption; the 4.8 × 10^-9 J₂-mismatch hazard is correctly avoided and documented in code. **PASS**.
  - 4.4 `grs80`: ✓ literals match Moritz 1980 Bull. Géod. 54 Table 1; correct use of `from_J2`. **PASS**.

**File verdict: PASS** — all four entry points correctly bind their convention's published values to the right uncertainty class (`defined` vs `measured`) and forward them through the right construction path (primary ctor vs `from_J2`) per each convention's *defining* parameter set. No C-fail surfaced. Three notes carried forward to downstream audit cards:

1. The `from_J2` inversion residual must be audited in `equipotential_ellipsoid.md` (REQ-EF-5).
2. The "WGS 72 J₂ mismatch avoided" claim in the header comment is testable by a `from_J2` round-trip unit test — recommend a verification anchor in `tests/test_constants_provider/`.
3. ω literal divergence between GRS 80 (7e-7 sig figs) and WGS 72/84 (11 sig figs) is a *convention-level* difference, not an implementation bug; documented in the GRS 80 NOTES above.
