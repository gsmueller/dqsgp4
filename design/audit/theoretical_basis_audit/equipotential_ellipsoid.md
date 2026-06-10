# Theoretical Basis Audit — `src/geodesy/equipotential_ellipsoid.h`

**File**: `src/geodesy/equipotential_ellipsoid.h` (308 lines)
**Functions audited**: 4 (1 ctor, 1 factory, 2 methods)
**Audit date**: 2026-05-13
**Status**: PASS-with-notes — all four formulas have cited theory and matched methods; bound verdicts are ✓ for series-driven slots and ⚠ for the two iterative paths (`from_J2`, ctor cube-root) where the Kantorovich-style residual is added but not all input-error propagation slots are exercised. One C-fail concern flagged for `from_J2` against REQ-EF-5.

---

## 1. `EquipotentialEllipsoid(a, 1/f, GM, ω, series_tolerance)` — WGS84 constructor

```
=== FORMULA AUDIT CARD ===
ID:                     equipotential_ellipsoid::ctor_inv_f
Location:               src/geodesy/equipotential_ellipsoid.h:89-189
Mathematical statement: Given the four defining parameters (a, 1/f, GM, ω) in
                        parameter set A (Definition 14.2.1), compute the full
                        derived constants of the equipotential ellipsoid:
                          - Geometric:  f, e², e, e', e'², b, E_lin, c
                          - Series:     q₀, q₀', u₀-series  (Ch 14, Thm 14.4.1–14.4.4)
                          - Physical:   m, γ_e, γ_p, k_som, U₀
                          - Radii:      R₁, R₃ = ∛(a²b)

THEORY
  Underlying theorem:   Composition of:
                        (a) Definition 14.3.1 of e² = 2f - f² (closed form).
                        (b) Theorem 14.4.1 (series for 2q₀) — alternating series
                            2q₀ = Σ_{n≥1} 4(-1)^{n+1}n / ((2n+1)(2n+3)) · e'^{2n+1}.
                        (c) Theorem 14.4.2 (series for q₀') — alternating series
                            q₀' = Σ_{n≥1} 6(-1)^{n+1} / ((2n+1)(2n+3)) · e'^{2n}.
                        (d) Theorem 14.4.4 (series for arctan(e')/e' driving U₀).
                        (e) Heiskanen–Moritz (1967) §2-7..2-10 closed forms for
                            γ_e (HM Eq. 2-73), γ_p (HM Eq. 2-74), Somigliana k
                            (HM Eq. 2-78), U₀ (HM Eq. 2-61).
                        (f) Newton iteration for x³ = a²b (real cube root) —
                            corresponds to fixed-point map x_{n+1} = (2x + a²b/x²)/3.
  Primary reference:    Heiskanen & Moritz (1967), Physical Geodesy §§2-7..2-10;
                        Moritz (1980), Geodetic Reference System 1980 §§2-3;
                        NGA.STND.0036 (2014) (NIMA TR 8350.2 successor) Ch 3 & App B;
                        design/derivations/ch14_equipotential_ellipsoid.md.
  Domain of validity:   1/f > 0 (so 0 < f < 1), so 0 < e² < 1, e' bounded;
                        for Earth-like inputs e'² ≈ 6.7e-3 ⇒ series geometric
                        ratio e'² ≪ 1 ⇒ rapid convergence (Cor 14.4.2).

METHOD
  Method declared:      (a) Closed form e² = 2f - f², b = a(1-f), c = a²/b,
                            e' = e/√(1-e²) — all algebraic, no truncation.
                        (b)/(c) `math::alternating_series` over the q₀-term and
                            q₀'-term lambdas (Leibniz-bound mode since
                            convergence_ratio defaults to -1).
                        (d) `math::alternating_series` for the u₀ Maclaurin
                            (Thm 14.4.4) ∑(-1)^n e'^{2n}/(2n+1).
                        (e) Closed-form algebraic combinations of GM, a, b, m,
                            e', q₀, q₀' for γ_e, γ_p, k_som, U₀.
                        (f) Newton cube-root iteration to series_tolerance,
                            capped at 60 iterations (matches series_sqrt style).
  Method implemented:   src/geodesy/equipotential_ellipsoid.h:100-188.
                        (a) Lines 101-108: pure algebraic chain, uses sqrt() from
                            TrackedValue (closed form Newton-bound elsewhere).
                        (b) Lines 112-124: q0_term lambda computes the n-th term
                            with `int sign`, `ratio<T>(...)`, and a `for` loop
                            building e'^{2n+1} by repeated multiplication. Sum
                            via `alternating_series(1, q0_term, series_tolerance)`.
                        (c) Lines 128-138: same structure for q0p_term.
                        (d) Lines 156-167: same structure for u0_term (note: this
                            series starts at n=0 — see NOTES re Theorem 14.4.4).
                        (e) Lines 141-152: m, γ_e, γ_p, k_som via closed-form
                            TrackedValue arithmetic (HM-style).
                        (f) Lines 175-188: Newton cube-root x_{n+1} = (2x + a²b/x²)/3
                            starting from x₀ = a; stops when |Δ| < tol.
  Match verdict:        ✓ matched — every leg of the construction is the method
                        cited by Ch 14 / HM 1967 / Moritz 1980:
                        - geometric chain is closed-form (no Padé, no CF);
                        - 2q₀ / q₀' / u₀ are alternating Maclaurin series, not
                          rational approximants;
                        - γ_e, γ_p, k_som are the canonical HM closed forms;
                        - R₃ is Newton-on-cubic-root, the standard iterative root.

ERROR BOUND
  Bound category:       precision (truncation + Newton residual) flows through
                        TrackedValue per REQ-EF-3 for the closed-form legs.
                        Series legs: per REQ-EF-6 (truncated alternating).
                        Newton cube root: per REQ-EF-5 (iterative residual).
  Bound formula:        (a) Closed-form steps: bound is propagated through the
                            REQ-EF-3 closed-form rules of TrackedValue ops.
                        (b)/(c)/(d) Leibniz-form bound: |R_N| ≤ |first omitted term|,
                            added inside `alternating_series` once truncation
                            magnitude < tolerance (series.h:53-56). This is the
                            rigorous Leibniz bound for an alternating series of
                            monotonically decreasing magnitudes (Cor 14.4.2).
                        (e) Closed-form REQ-EF-3 propagation.
                        (f) Newton cube root: |R| ≤ |Δ_final| at convergence, the
                            Kantorovich/Ostrowski residual for quadratically
                            convergent iterations (REQ-EF-5).
  Bound implemented:    (b)/(c)/(d) Inside `math::alternating_series`
                            (src/math/series.h:53-56): on convergence
                            `sum.errors.precision += truncation` with
                            `truncation = abs(last_term.value)` for the default
                            (Leibniz) mode used here. ✓
                        (e) TrackedValue ops propagate errors per their REQ-EF-3
                            rules. ✓ inherited.
                        (f) Lines 183-185: at convergence, the code adds
                            `correction = abs((x_new - x).value)` to
                            `x.errors.precision` via the construction
                            `math::TrackedValue<T>(correction, T(0), correction, T(0)).errors.precision`,
                            i.e. it ADDS the bound exactly to precision.
                            This matches REQ-EF-5 (iterative residual into precision).
  Bound verdict:        ✓ matched for legs (a)–(e); ⚠ Newton-cube-root residual is
                        added correctly in shape but only on the *value-magnitude*
                        of Δ; the convergence-region Kantorovich majorant
                        |Δ|/(1-q) is omitted. For the cube-root map near
                        x* ≈ a∛(1-f), q ≈ |1 - f| · O(1) ≪ 1, so the practical
                        gap is ~q · |Δ| ≪ tolerance — bound remains rigorous
                        within type precision for double. Flag for tightening
                        if used with wide T (cpp_bin_float_50).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form rules), REQ-EF-5 (Newton residual),
                        REQ-EF-6 (Taylor/Leibniz truncation).
  AUD-EF applies:       AUD-EF-4 (iterative residual to precision),
                        AUD-EF-5 (Taylor branches add truncation bound),
                        AUD-EF-7..10 (per-category closed-form propagation).
  AUD-MC applies:       n/a — this is geodetic data construction, not algebra.
  Verification test:    tests/test_geodesy/ should exercise:
                          (i) e²(WGS84) reproduces NGA.STND.0036 value to
                              type precision;
                          (ii) q₀, q₀' match Moritz (1980) Table 1 to declared
                              precision;
                          (iii) γ_e, γ_p match Moritz (1980) Table 2 to declared
                              precision;
                          (iv) R₃ from a, b matches a · (1-f)^{1/3} to declared
                              precision.

NOTES
  - The U₀ helper series starts at n=0 (line 166: `alternating_series<T>(0, ...)`),
    which matches Theorem 14.4.4 (the Maclaurin of arctan(e')/e' starts at n=0,
    not n=1). The 2q₀ and q₀' lambdas start at n=1 per Theorems 14.4.1/14.4.2.
    All three start indices match the cited theorems exactly.
  - The cube-root Newton uses `int i = 0; i < 60` matching the cap used in
    `math::series_sqrt` (series.h:135). The starting guess x₀ = a is within
    a factor (1-f)^{1/3} ≈ 0.9989 of the root for terrestrial flattening, so
    quadratic convergence kicks in immediately; typically 3-4 iterations
    suffice for double.
  - **C-pass per Theorem 2.1 of the framework**: theory cited (Ch 14 series
    representations + HM closed forms + Newton root), method implemented
    matches (no Padé, no CF), bounds match Leibniz / closed-form / Newton-
    residual as required.
```

---

## 2. `EquipotentialEllipsoid::from_J2(a, J₂, GM, ω, series_tolerance)` — GRS80/WGS72 factory

```
=== FORMULA AUDIT CARD ===
ID:                     equipotential_ellipsoid::from_J2
Location:               src/geodesy/equipotential_ellipsoid.h:214-265
Mathematical statement: Given (a, J₂, GM, ω), solve Brouwer's formula
                          J₂ = (e²/3) · (1 − 2me'/(15q₀))
                        for e² by fixed-point iteration starting from
                        e²₀ = 3J₂ and updating
                          e²_new = 3J₂ + (2/15) · m(e²) · e'(e²) · e² / q₀(e²),
                        then derive 1/f from e² and forward to the WGS84 ctor.

THEORY
  Underlying theorem:   Theorem 14.5.1 (Brouwer's formula) with
                        Corollary 14.5.1 (fixed-point form):
                          e² = 3J₂ + (2m e' e²) / (15 q₀)
                        and convergence by the Banach fixed-point theorem
                        (Ch 1, §1.6) — contraction ratio dominated by
                        2m/(15q₀) ≈ O(10⁻²) ≪ 1 for terrestrial bodies, so
                        the map is contractive on (0, 1).
                        SR3 (Hoots & Roehrich 1980) Appendix §B uses the same
                        WGS72 convention with J₂ as a direct input.
  Primary reference:    Heiskanen & Moritz (1967), §2-10 Eqs. 2-90, 2-92;
                        Moritz (1980), pp. 129-131;
                        Hoots & Roehrich (1980), SPACETRACK Report No. 3,
                        Appendix B (WGS72 constants);
                        design/derivations/ch14_equipotential_ellipsoid.md
                        §14.5 (Theorem 14.5.1, Corollary 14.5.1).
  Domain of validity:   J₂ > 0, 3J₂ < 1, and the contraction condition
                        |g'(e²)| < 1 (guaranteed for terrestrial parameters).

METHOD
  Method declared:      Fixed-point iteration of the Brouwer map g(e²) per
                        Cor 14.5.1; cap 20 iterations; stop when
                        |e²_new − e²_guess| < series_tolerance.
                        Inside each iteration:
                          - e' = √(e²/(1-e²))                  (closed form)
                          - q₀ via 2q₀ alternating series       (Thm 14.4.1)
                          - b ≈ a√(1-e²)                       (closed form)
                          - m = ω²a²b/GM                       (closed form)
  Method implemented:   src/geodesy/equipotential_ellipsoid.h:223-257:
                        Loop `for (int iter = 0; iter < 20; ++iter)` computes
                        e_prime_iter, q0_iter (via the same alternating_series
                        lambda used in the ctor), b_iter, m_iter, and assembles
                          e2_new = 3*J2_in
                                 + 2 * m_iter * e_prime_iter * e2_guess
                                   / (15 * q0_iter).
                        Note: line 250 actually computes the update as
                          (2 · m · e' · e²_guess) / (15 · q₀),
                        which matches Cor 14.5.1's fixed-point form
                          e² = 3J₂ + 2me'e²/(15q₀)  ✓
                        Stops on `correction < series_tolerance` (line 256).
                        After convergence: lines 260-264 compute f from e²,
                        then 1/f, then `return EquipotentialEllipsoid(a, inv_f,
                        GM, ω, series_tolerance)` to forward into the WGS84
                        constructor and re-derive all downstream constants
                        from 1/f.
  Match verdict:        ✓ matched — Banach fixed-point iteration of the
                        Brouwer/Moritz map, NOT a series expansion in J₂, NOT
                        a Padé or continued fraction. The inner q₀ is the
                        cited Theorem 14.4.1 series. The contraction-mapping
                        framing is the cited Cor 14.5.1 framing.

ERROR BOUND
  Bound category:       precision (Banach contraction residual + inner-series
                        truncation, per REQ-EF-5 and REQ-EF-6).
  Bound formula:        For a Banach fixed-point iteration with contraction
                        constant q < 1, the rigorous post-convergence bound is
                          |e²_k − e²_*| ≤ |Δ_k| / (1 − q)
                        where Δ_k is the final correction. For the Brouwer
                        map q ≈ 2m/(15q₀) ≈ O(10⁻²), so |Δ_k|/(1-q) ≈ |Δ_k|.
                        The bound to ADD to precision per REQ-EF-5 is the
                        final correction magnitude (or, more conservatively,
                        |Δ_k|/(1-q)).
  Bound implemented:    ✗ NO bound is added by `from_J2`. Lines 252-257 only
                        check `correction < series_tolerance` and `break`;
                        the loop body does not touch `e2_guess.errors.precision`,
                        and no precision contribution from the J₂→e² iteration
                        residual is recorded before forwarding to the ctor on
                        line 264. The downstream ctor will compute its OWN
                        precision contributions for the q₀, q₀', u₀ series and
                        for the cube-root Newton, but the iteration residual
                        between (a, J₂, GM, ω) and the derived 1/f is dropped.
  Bound verdict:        ✗ C-FAIL per Theorem 2.1 of the framework. The
                        Kantorovich-style residual |Δ_final| of the J₂→e²
                        Banach iteration is NOT added to any error category.
                        REQ-EF-5 explicitly requires this. The downstream
                        `EquipotentialEllipsoid` it returns thus carries an
                        under-counted precision bound on every derived
                        constant — including e², b, q₀, J_{2n}, γ_e, γ_p.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-5 (Convergence residual added to precision) —
                        VIOLATED by current implementation.
                        REQ-EF-6 (Taylor/Leibniz truncation) — inherited
                        correctly via the inner `alternating_series` call.
  AUD-EF applies:       AUD-EF-4 (Iterative algorithms add residual) — this
                        is the test that should catch the omission above.
  AUD-MC applies:       n/a.
  Verification test:    tests/test_geodesy/ should:
                          (i) Construct from GRS80 J₂ and confirm e² matches
                              Moritz (1980, p. 131) value 0.006 694 380 022 90
                              to type precision after ≤ 4 iterations
                              (Example 14.5.1 in Ch 14).
                          (ii) Construct from SGP4 WGS72 J₂ = 0.001082616 and
                              confirm derived 1/f differs from the WGS72-A
                              value 298.26 by the documented ~4.8e-9 gap
                              (this is the "Key Finding" in the file header).
                          (iii) Confirm `e2.errors.precision` includes a
                              contribution from the J₂→e² iteration residual.
                              CURRENTLY EXPECTED TO FAIL (see Bound verdict).

NOTES
  - The fix to close the C-fail is straightforward: after the loop, before
    forwarding to the ctor, add the final correction magnitude to
    `e2_guess.errors.precision`, propagate that to `f_derived` and
    `inv_f_derived` via the existing TrackedValue arithmetic (which already
    handles closed-form propagation per REQ-EF-3), and forward the
    resulting `inv_f_derived` with its now-correct precision to the ctor.
    Pattern: same as `series_sqrt` (series.h:147) and the cube-root Newton
    in the ctor (line 184).
  - The note in REQ-EF-5: "use Kantorovich residual or final-correction
    magnitude; both are rigorous" supports either form. The simpler
    `|Δ_final|` form is recommended here for consistency with `kepler.h`
    and `series_sqrt`.
  - The b iteration uses the CONVENTIONAL Moritz form
      b ≈ a√(1-e²),
    rather than the exact b = a(1-f), because f hasn't been determined yet.
    This is consistent with Cor 14.5.1's iteration design.
  - **C-fail per Theorem 2.1**: theory and method match, but the implemented
    bound omits the iteration residual. The function passes A (error wiring
    of inputs flows through), passes B (no algebra ops to break), but fails C.
    Downstream `total_error()` of every derived constant from a `from_J2`
    constructed ellipsoid is under-counted by ~|Δ_final| · O(downstream
    Jacobian). For Earth at series_tolerance = 1e-16, |Δ_final| ≈ 1e-16
    (~machine epsilon), so the practical magnitude of the gap is at the
    type-epsilon level — the C-fail is structural (REQ-EF-5 violation),
    not numerically catastrophic.
```

---

## 3. `normal_gravity(phi)` — Somigliana formula

```
=== FORMULA AUDIT CARD ===
ID:                     equipotential_ellipsoid::normal_gravity
Location:               src/geodesy/equipotential_ellipsoid.h:271-278
Mathematical statement: γ(φ) = γ_e · (1 + k·sin²φ) / √(1 − e²·sin²φ)
                        — normal gravity on the ellipsoid at geodetic
                        latitude φ, via the Somigliana closed form.

THEORY
  Underlying theorem:   Somigliana's formula (Heiskanen & Moritz 1967, §2-9,
                        Eq. 2-78). γ on the ellipsoid is the magnitude of the
                        gradient of the normal potential U evaluated on the
                        ellipsoid surface; reduction of the gradient via the
                        Somigliana constant k = bγ_p/(aγ_e) − 1 produces the
                        closed form above.
                        Ch 14, §14.7 (Theorem 14.7.1 / Eq. 14.7.x in our
                        derivation document).
  Primary reference:    Heiskanen & Moritz (1967), §2-9 Eq. 2-78;
                        Moritz (1980), §3 p. 138, Eq. (4.1);
                        NGA.STND.0036 (2014) §4 (gravity on ellipsoid);
                        design/derivations/ch14_equipotential_ellipsoid.md §14.7.
  Domain of validity:   φ ∈ [-π/2, π/2]; 0 < e² < 1 (true for any physical
                        ellipsoid); 1 − e²·sin²φ > 0 (always, since
                        e²·sin²φ ≤ e² < 1). Closed form, no truncation, no
                        iteration.

METHOD
  Method declared:      Direct closed-form evaluation of γ_e · N / D where
                          N = 1 + k_som · sin²φ
                          D = √(1 − e² · sin²φ),
                        with sin via TrackedValue's `sin` and √ via TrackedValue's
                        `sqrt` (each propagating closed-form bounds per REQ-EF-3).
  Method implemented:   src/geodesy/equipotential_ellipsoid.h:273-277:
                          sin_phi = sin(phi)
                          sin2_phi = sin_phi * sin_phi
                          numerator = 1 + k_som * sin2_phi
                          denominator = sqrt(1 − e2 * sin2_phi)
                          return gamma_e * numerator / denominator.
                        This is *literally* the Somigliana formula with no
                        algebraic rearrangement.
  Match verdict:        ✓ matched — pure closed-form Somigliana, no series,
                        no Padé, no rational approximant of the formula.

ERROR BOUND
  Bound category:       precision (closed-form REQ-EF-3 propagation from
                        inputs γ_e, k_som, e², and φ, plus the closed-form
                        bound of sin and sqrt via TrackedValue ops).
  Bound formula:        Composition of REQ-EF-3 closed-form propagation rules:
                          - sin(φ): δ_p propagates as |cos(φ)| · δ_p^φ (Lipschitz
                            with |cos|≤1; the TrackedValue::sin op handles this).
                          - x²: δ_p propagates as 2|x| · δ_p^x.
                          - 1 ± y: δ_p propagates additively.
                          - sqrt(z): δ_p propagates as δ_p^z / (2√z).
                          - Multiplication & division: standard product / quotient
                            rules.
                        All rules are the rigorous closed-form bounds per REQ-EF-3
                        (triangle inequality + first-order Taylor for monotonic
                        Lipschitz functions); no method-specific truncation arises.
  Bound implemented:    ✓ Bound is whatever the underlying TrackedValue operations
                        provide — the function does not add any additional
                        truncation or iteration bound (correctly, since none is
                        needed: every step is closed form). The precision bound
                        on the result is exactly the REQ-EF-3 composition of
                        precision bounds on γ_e, k_som, e², and φ.
  Bound verdict:        ✓ matched — closed-form theorem (Somigliana) is
                        implemented as a closed-form arithmetic expression, with
                        propagation entirely inherited from TrackedValue. There
                        is no "next-omitted-term" or iteration residual to add.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form propagation through standard funcs).
  AUD-EF applies:       AUD-EF-7..10 (per-category propagation rules of `sin`,
                        `sqrt`, multiplication, division).
  AUD-MC applies:       n/a (no algebra-axiom property at stake; this is a
                        scalar formula).
  Verification test:    tests/test_geodesy/ should:
                          (i) normal_gravity(0) = γ_e exactly (sin²0 = 0
                              ⇒ numerator = 1, denominator = 1);
                          (ii) normal_gravity(π/2) = γ_p exactly (uses the
                              identity 1 + k = bγ_p/(aγ_e), giving γ_p
                              after the division by √(1-e²) = b/a);
                          (iii) For WGS84, normal_gravity(φ) reproduces the
                                Moritz Table 3 / NGA.STND.0036 §4 values to
                                type precision over a grid of φ.

NOTES
  - The Somigliana formula is the canonical closed-form expression of γ on
    the ellipsoid; no series or iterative method is appropriate here.
    Theorem-method-bound triad is trivially consistent.
  - The "Lipschitz / mean-value" entry in §8's detection rule table
    (theoretical_basis_audit.md §8) classifies sin's bound rule — both the
    Lipschitz form |cos(ξ)|·δ and the mean-value form are rigorous at first
    order; the TrackedValue::sin op's choice is auditable separately as
    part of the tracked_value.h card.
  - **C-pass per Theorem 2.1**: closed-form theorem cited (Somigliana),
    closed-form method implemented, closed-form REQ-EF-3 propagation
    inherited correctly.
```

---

## 4. `J2n(n)` — Even zonal harmonic

```
=== FORMULA AUDIT CARD ===
ID:                     equipotential_ellipsoid::J2n
Location:               src/geodesy/equipotential_ellipsoid.h:282-304
Mathematical statement: For n ≥ 1,
                          J_{2n} = (-1)^{n+1} · 3 e^{2n} / ((2n+1)(2n+3))
                                   · (1 − n + 5n J₂/e²),
                        with the n=1 case computed directly from
                          J₂ = (e²/3) · (1 − 2me'/(15q₀))
                        (Theorem 14.5.1) to avoid circular reference to
                        itself in the (1 − n + 5nJ₂/e²) inner factor.

THEORY
  Underlying theorem:   Theorem 14.6.1 (Heiskanen & Moritz Eq. 2-92; Moritz
                        (1980) p. 130). J_{2n} of the level ellipsoid is the
                        spherical-harmonic coefficient determined by the
                        ellipsoidal shape and J₂ — not an independent input.
                        For n = 1, the formula reduces to the Brouwer relation
                        Theorem 14.5.1 above.
  Primary reference:    Heiskanen & Moritz (1967), §2-10 Eq. 2-92;
                        Moritz (1980), p. 130;
                        design/derivations/ch14_equipotential_ellipsoid.md
                        §14.6 (Theorem 14.6.1).
  Domain of validity:   n ≥ 1; 0 < e² < 1; finite J₂. The formula is exact
                        for the level ellipsoid (not a truncated series),
                        so domain is just "physical-ellipsoid parameters".

METHOD
  Method declared:      Direct closed-form evaluation of Eq. 14.6.1 for n ≥ 2;
                        direct closed-form evaluation of Brouwer's J₂ formula
                        for n = 1 (special-cased to avoid the circular
                        J_{2n=1} → J₂ self-reference).
                        e^{2n} computed by REPEATED MULTIPLICATION of e² (no
                        `pow`, no Taylor expansion in e).
  Method implemented:   src/geodesy/equipotential_ellipsoid.h:285-303:
                        Line 287-288: J₂ via Brouwer Theorem 14.5.1.
                        Line 290: short-circuit `if (n == 1) return J2;`.
                        Lines 295-298: e^{2n} built by an `n`-step loop
                                       multiplying e² in.
                        Line 300: `prefactor = ratio<T>(sign*3, (2n+1)(2n+3))`.
                        Line 301: `inner = 1 − n + 5n · J₂/e²` (the Eq. 14.6.1
                                  inner factor).
                        Line 303: return prefactor · e^{2n} · inner.
                        Every operation is closed-form TrackedValue arithmetic;
                        no `pow` (correctly — repeated multiplication is more
                        precise for small integer powers), no series.
  Match verdict:        ✓ matched — closed-form Theorem 14.6.1 implemented as
                        closed-form arithmetic.

ERROR BOUND
  Bound category:       precision (closed-form REQ-EF-3 propagation only).
                        For n = 1: the precision bound inherits the q₀ series
                        truncation bound (added during construction).
  Bound formula:        REQ-EF-3 composition of closed-form rules:
                          - Repeated multiplication of e² (n times): δ_p scales
                            as n · |e²|^{n-1} · δ_p^{e²} to first order.
                          - Division by integer denominators (2n+1)(2n+3): exact
                            (TrackedValue with zero precision error since it's
                            `ratio<T>`).
                          - 1 − n + 5n · J₂/e²: standard sum/difference/quotient.
                          - Outer multiplication: standard product rule.
                        No method-specific truncation or iteration arises — the
                        formula is a finite closed-form expression in the
                        stored constants (e², J₂ for n ≥ 2; q₀, m, e' for n = 1).
  Bound implemented:    ✓ Bound is inherited entirely from the TrackedValue
                        operations and from `e2`, `J2`, `q0`, `m_const`, `e_prime`
                        precision bounds set at construction. The function adds
                        no additional bound (correctly: closed-form theorem).
  Bound verdict:        ✓ matched — closed-form theorem implemented as
                        closed-form arithmetic, with REQ-EF-3 propagation
                        composed correctly. The n = 1 special case correctly
                        avoids the self-reference and inherits the q₀ series
                        truncation bound from `q0.errors.precision` (set in
                        the ctor).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (Closed-form propagation).
  AUD-EF applies:       AUD-EF-7..10 (multiplicative / additive propagation).
  AUD-MC applies:       n/a.
  Verification test:    tests/test_geodesy/ should:
                          (i) J2n(1) reproduces J₂ from the constructor's
                              stored Brouwer evaluation (and matches GRS80
                              specification value to type precision when
                              constructed via `from_J2`);
                          (ii) J2n(2), J2n(3) reproduce Moritz (1980) p. 130
                              tabulated values for GRS80 to type precision;
                          (iii) Sign alternates (-1)^{n+1} for n = 1..5.

NOTES
  - The n = 1 special case is REQUIRED, not optional: substituting n = 1 into
    the inner factor would give (1 - 1 + 5 J₂/e²) = 5 J₂/e², which combined
    with the prefactor 3 e² / (3 · 5) = e²/5 gives 5J₂/e² · e²/5 = J₂. So
    the formula would reduce to J₂ itself ALGEBRAICALLY — but in the code
    `J2` is just-computed from the Brouwer formula at line 287, so the
    "circular reference" comment is correct in spirit but algebraically a
    tautology. The short-circuit at line 290 is a CLARITY measure, not a
    correctness measure. (The algebraic identity is part of why
    Theorem 14.6.1 is consistent with Theorem 14.5.1.)
  - Repeated multiplication of e² rather than `pow(e2, n)` is the right
    choice: it preserves the closed-form REQ-EF-3 propagation through each
    individual product without invoking a transcendental `pow` bound rule.
  - **C-pass per Theorem 2.1**: closed-form theorem cited (Eq. 14.6.1),
    closed-form method implemented (no Padé, no series), bound inherited
    correctly from inputs.
```

---

## 5. File-level verdict

### Triad summary

| Function | Theory cited | Method matches theory | Bound matches method | Verdict |
|---|---|---|---|---|
| ctor `(a, 1/f, GM, ω)` | Ch 14 §§14.3-14.4, HM §§2-7..2-10, Newton root | ✓ closed-form + alt-series + Newton root | ✓ Leibniz / REQ-EF-3 / Newton residual | ✓ PASS |
| `from_J2(a, J₂, GM, ω)` | Thm 14.5.1, Cor 14.5.1, Banach fixed point | ✓ contraction-map iteration | ✗ Banach residual not added (REQ-EF-5) | ✗ C-FAIL |
| `normal_gravity(φ)` | Somigliana, HM Eq. 2-78 | ✓ pure closed form | ✓ REQ-EF-3 inherited | ✓ PASS |
| `J2n(n)` | Thm 14.6.1, HM Eq. 2-92 | ✓ pure closed form (with n=1 special case) | ✓ REQ-EF-3 inherited | ✓ PASS |

### Findings

- **3 / 4 PASS, 1 / 4 C-FAIL.**
- **C-fail**: `from_J2` does not add the Banach fixed-point residual
  `|Δ_final|` (REQ-EF-5) to the precision of the derived `e²` (or onward
  to `1/f` before forwarding to the WGS84 ctor). Per Theorem 2.1 of the
  framework, this invalidates the `total_error()` claim of every
  downstream-derived constant obtained via the `from_J2` path. The
  magnitude is type-epsilon-level at default tolerance, but the
  requirement is structural.

### REQ-EF-5 remediation note

The fix is mechanical and matches existing patterns in
`src/math/kepler.h:87` and `src/math/series.h:147-148`:

```
// After the for-loop, before deriving inv_f:
e2_guess.errors.precision = e2_guess.errors.precision + correction_final;
```

where `correction_final` is the magnitude of the final `e2_new - e2_guess`
delta. The existing `f_derived` and `inv_f_derived` computations
already use TrackedValue arithmetic, so the precision will propagate
correctly through them and into the forwarded constructor.

### No method-theory mismatches detected

None of the four functions implement a method outside its cited theory.
In particular, the §8 detection-rule table of the framework was checked:
- No `numerator / denominator` patterns suggesting Padé approximation
  (the only quotient structures are the canonical HM closed forms and
  the `inner = 1 − n + 5nJ₂/e²` factor of Eq. 14.6.1).
- No `for` loops summing terms where a closed-form identity is cited
  (the only loops are: (i) alternating-series evaluation, correctly
  cited as series; (ii) Newton iterations, correctly cited as iterative;
  (iii) the n-step e^{2n} construction, which is repeated multiplication
  by definition, not a series approximation).
- No continued fractions, no rational approximants of `arctan` or other
  transcendentals; the only transcendental that appears (atan, sin, sqrt)
  enter via TrackedValue's own ops and are audited as part of
  `tracked_value.h`.

### Cross-references

- Framework: `design/audit/theoretical_basis_audit.md` §1 (audit card schema),
  §2 (Theorem 2.1 verification triad), §5 (worked example), §8 (detection
  rules).
- Theory document: `design/derivations/ch14_equipotential_ellipsoid.md`
  §§14.2 (defining parameter sets), 14.3 (geometric chain), 14.4 (q₀, q₀',
  U₀ series), 14.5 (Brouwer formula + Cor 14.5.1 iteration), 14.6
  (Theorem 14.6.1 for J_{2n}), 14.7 (Somigliana).
- Companion source: `src/math/series.h` (`alternating_series`,
  `series_sqrt`), `src/math/tracked_value.h` (REQ-EF-3 closed-form
  propagation), `src/math/kepler.h` (canonical iterative-residual pattern
  cited by §6 of framework).
- External: Heiskanen & Moritz (1967) Physical Geodesy §§2-7..2-10;
  Moritz (1980) GRS80 §§2-3 and Tables 1-3; NGA.STND.0036 (2014)
  Ch 3 / App B (NIMA TR 8350.2 successor); Hoots & Roehrich (1980)
  SPACETRACK Report No. 3 Appendix B (WGS72 J₂ direct-input convention).

### Remediation track

1. **`from_J2` REQ-EF-5 fix** (C-FAIL → C-PASS): add Banach residual to
   `e2_guess.errors.precision` before forwarding. ~5 LoC. Update file
   status to PASS after the fix lands and a test exercises the
   precision-inclusion path.
2. **Wide-T tightening (optional)**: the cube-root Newton in the ctor
   omits the `|Δ|/(1-q)` correction factor; rigorous within type
   precision for double but should be made explicit for wide T. Flag
   only; not a current C-fail.

After (1) lands, this file moves from PASS-with-notes to **PASS**.
