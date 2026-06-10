# Theoretical Basis Audit — `src/orbit/state_from_elements.h`

**Scope.** One function in this file produces orbit state output: `elements_to_state`. One audit card.

**Companion to.** `design/audit/theoretical_basis_audit.md` (framework §1, §5).

---

## 1 `elements_to_state(r_er, u, i, Omega, rdot, rfdot, re_km)`

```
=== FORMULA AUDIT CARD ===
ID:                     state_from_elements::elements_to_state
Location:               src/orbit/state_from_elements.h:39-89
Mathematical statement: Given corrected (r, u, i, Ω, ṙ, rφ̇), produce
                        TEME position r⃗ = r·û  and velocity
                        v⃗ = ṙ·û + rφ̇·v̂ in km and km/s.
                        Orientation built from
                            û = R₃(−Ω) R₁(−i) R₃(−u) · ê_x
                            v̂ = R₃(−Ω) R₁(−i) R₃(−u) · ê_y
                        via direct sin/cos of u, i, Ω (no matrix
                        multiplication). Components:
                          xmx = −sin(Ω)·cos(i)
                          xmy =  cos(Ω)·cos(i)
                          ûx  = xmx·sin(u) + cos(Ω)·cos(u)
                          ûy  = xmy·sin(u) + sin(Ω)·cos(u)
                          ûz  =  sin(i)·sin(u)
                          v̂x  = xmx·cos(u) − cos(Ω)·sin(u)
                          v̂y  = xmy·cos(u) − sin(Ω)·sin(u)
                          v̂z  =  sin(i)·cos(u)
                        Unit conversion: position × re_km [km];
                        velocity × re_km/60 [km/s].

THEORY
  Underlying theorem:   Composition of three elementary rotations
                        R₃(−Ω) R₁(−i) R₃(−u) is the standard
                        perifocal-to-inertial transformation. The
                        column-extraction trick (û, v̂ directly as
                        rotated basis vectors) is an algebraic
                        identity, not a series.
  Primary reference:    Hoots & Roehrich (1980), Spacetrack Report
                        No. 3, pp. 14-15 (xmx, xmy, ux, uy, uz, vx,
                        vy, vz construction verbatim). Vallado §8
                        documents TEME frame definition. Battin
                        (1999) Ch. 2 covers the underlying Euler-
                        rotation algebra.
  Domain of validity:   All (u, i, Ω) ∈ ℝ³; no singularity in the
                        formula (rotations are entire). Result is
                        meaningful for r_er > 0; signs of ṙ, rφ̇
                        unrestricted.

METHOD
  Method declared:      Closed-form algebraic identity — the
                        components of the product rotation applied
                        to ê_x and ê_y, expressed entirely in
                        sin/cos of the three Euler angles. No
                        truncation, no iteration.
  Method implemented:   Six sin/cos calls, then closed-form
                        polynomial combinations matching SR3
                        formulas line-for-line:
                          line 62-63: xmx, xmy
                          line 65-67: ux, uy, uz
                          line 69-71: vx, vy, vz
                          line 74-76: x, y, z (= r·û·re_km)
                          line 80-83: xdot, ydot, zdot
                                       (= (ṙ·û + rφ̇·v̂)·re_km/60)
                        Sign conventions match SR3 page 14-15.
  Match verdict:        ✓ matched — implementation is a transcription
                        of SR3 §6 orientation-vector construction.

ERROR BOUND
  Bound category:       inherited from inputs (no new method-specific
                        bound introduced by this function)
  Bound formula:        REQ-EF-3 (closed-form propagation): every
                        arithmetic op on TrackedValue<T> (sin, cos, ×,
                        +, −) adds the appropriate per-op precision
                        bound from `tracked_value.h`. The function
                        itself introduces no truncation, no iteration,
                        no series — it only composes closed-form ops.
                        Total `errors.precision` accumulates triangle-
                        inequality-style across the dozen multiplies
                        and adds; no separate "method bound" is
                        required at this level.
  Bound implemented:    The function returns auto-typed TrackedValue
                        results from sin/cos/×/+/− chains. No
                        manual `.errors.precision += ...` call is
                        made here. Bounds flow through the operator
                        overloads of `TrackedValue<T>` as wired by
                        AUD-EF tests.
  Bound verdict:        ✓ matched — closed-form identity →
                        TrackedValue<T> arithmetic carries the bound
                        per REQ-EF-3. No method-specific bound owed.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-3 (closed-form per-op propagation).
                        No REQ-EF-5/6 obligation (no iteration,
                        no Taylor truncation).
  AUD-EF applies:       AUD-EF-1..3 (closed-form chain wiring) —
                        validated transitively through the TrackedValue
                        operator overloads used here.
  AUD-MC applies:       n/a — this is not an algebra-operator
                        implementation; it is an orientation-vector
                        composition. Verified at the integration
                        level by SGP4 end-to-end tests where the
                        propagated state must agree with SR3
                        reference vectors.
  Verification test:    tests/sgp4 end-to-end tests of state output
                        against SR3 reference; any TEME-frame
                        check at the propagator boundary exercises
                        this function.

NOTES
  - The "xmx, xmy" intermediates are precisely the SR3 names; the
    code preserves the literal SR3 variable naming. This makes the
    "transcription of SR3 §6" verification a textual diff.
  - The factor re_km/60 in vkmpersec implicitly carries the
    convention that input ṙ, rφ̇ are in [ER/min] (per the docstring).
    Units are documented but not type-checked; this is an interface
    contract, not a theoretical-basis concern.
  - The two factors `r_er * u{x,y,z} * re_km` could be re-ordered
    for minor numerical preference, but every multiplication order
    is closed-form and gives the same rigorous bound.
  - There is no small-angle branch here because no function in the
    formula has a removable singularity. sin/cos of large angles
    are evaluated by the host-T std library; their precision is
    captured by `TrackedValue<T>::sin/cos` per REQ-EF (not this
    function's concern).
  - No method-theory mismatch detected: theory says "compose three
    rotations and read off two columns"; code does exactly that.
  - One observation worth flagging for the file-level review: the
    formula in this implementation matches SR3 for the SGP4 path,
    but a dual-quaternion propagator (the user's primary engine in
    this codebase) would more naturally express the same rotation
    composition as a unit quaternion `q = q_Ω · q_i · q_u` and
    extract û, v̂ as `q · ê_{x,y} · q*`. The two are algebraically
    equivalent; the SR3-style component formula is retained here
    because this file appears in the SGP4 integration path.
```

---

## 2 File-level verdict

- **A. Error wiring**: ✓ all internal arithmetic flows through TrackedValue<T> operators (REQ-EF-3); no bare T leakage.
- **B. Algebra axioms**: n/a (orientation composition, not algebra op).
- **C. Theoretical basis**:
  - 2.1 `elements_to_state`: ✓ theory = SR3 §6 / Vallado §8 Euler composition; ✓ method = closed-form column extraction matching SR3 line-for-line; ✓ bound = inherited via REQ-EF-3.

**File verdict: PASS.** No method-theory mismatch; bound obligations satisfied at the per-op level via TrackedValue<T>.
