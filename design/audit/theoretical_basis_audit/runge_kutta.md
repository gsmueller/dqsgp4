# Theoretical Basis Audit — `src/integrators/runge_kutta.h`

**Scope.** Three formulas: `lie_advance_pose`, `euler`, `runge_kutta_4`.
**Framework.** `design/audit/theoretical_basis_audit.md` §1 (audit-card slots), §2 (theory ⇒ method ⇒ bound triad), §5 (worked example for Taylor-method conformance).
**Cross-audit anchors.** REQ-EF-7 (model truncation → accuracy), REQ-EF-15 (post-step retraction → precision), REQ-DQ-18 (half-angle screw exponential convention).

---

## 1 `lie_advance_pose`

```
=== FORMULA AUDIT CARD ===
ID:                     runge_kutta::lie_advance_pose
Location:               src/integrators/runge_kutta.h:54-63
Mathematical statement: Pose advance on SE(3):
                          pose_new = (pose · exp_se3(δt · ξ))_↾SE(3)
                        with ξ = (ω, v) the body twist, exp_se3 the
                        screw exponential map se(3) → SE(3), and
                        (·)_↾SE(3) the retraction (renormalization) to
                        the unit-DQ submanifold. Encoded in dual-
                        quaternion half-angle form:
                          step  = exp_screw((δt/2) · ω, (δt/2) · v)
                          M_new = Π_SE(3)( pose.M · step )
                        where Π_SE(3) is `Pose::normalized()`.

THEORY
  Underlying theorem:   The exponential map exp : se(3) → SE(3) is the
                        unique one-parameter subgroup integrating a
                        constant body twist on the Lie group SE(3).
                        For a constant body twist ξ ∈ se(3), the exact
                        integral curve of dM/dt = M · ξ̂ with M(0) = M_0
                        is M(t) = M_0 · exp(t ξ̂).
                        (Munthe-Kaas (1999), "High order Runge–Kutta
                        methods on manifolds"; equivalently Iserles,
                        Munthe-Kaas, Nørsett, Zanna (2000) "Lie-group
                        methods", Acta Numerica §3, Thm 3.1.)
                        Half-angle encoding: REQ-DQ-18 specifies
                          exp_screw(u, v) = cos|u| + sinc|u|·u
                                          + ε [ -(u·v) sinc|u|
                                                + (cos|u| · v
                                                   + (u×v) sinc|u|
                                                   + u · (u·v)
                                                         (1-cos|u|)/|u|²) ]
                        with u, v carrying the **half** rotation and
                        half translation; thus exp_screw(δt/2 · ω,
                        δt/2 · v) corresponds to a full-step screw of
                        magnitude δt · |ω|. This matches `exp_screw`'s
                        documented half-angle convention.
  Primary reference:    Munthe-Kaas (1999), "High order Runge–Kutta
                        methods on manifolds", Appl. Numer. Math. 29,
                        115-127, Theorem 2.1 (constant-twist exact
                        flow). Selig (2005), "Geometric Fundamentals
                        of Robotics" §3.4 for the SE(3) screw form.
                        Half-angle DQ encoding: Kavan et al. (2008)
                        "Geometric skinning with approximate dual
                        quaternion blending".
  Domain of validity:   Half-angle ball |u| = (δt/2)|ω| ≤ π/2, i.e.
                        a full-step rotation angle ≤ π. For larger
                        rotations the screw logarithm is multi-valued
                        and exp_screw remains correct as a forward
                        map, but log_screw round-trip is lost.
                        `pose.M` is assumed unit (REQ-DQ-15);
                        sub-ULP drift is removed by `normalized()`.

METHOD
  Method declared:      Exact closed-form SE(3) screw exponential
                        applied to a **constant** body twist over the
                        stage interval, followed by post-step retraction
                        via `Pose::normalized()`.
  Method implemented:   half_delta = delta_t / 2;
                        u = half_delta · twist.angular;
                        v = half_delta · twist.linear;
                        step = DualQuaternion::exp_screw(u, v);
                        return Pose(pose.M · step).normalized();
                        — verbatim closed-form, no series, no iteration.
  Match verdict:        ✓ matched — closed-form Lie-group flow + REQ-EF-15
                        retraction. No Taylor, no Padé, no continued
                        fraction. Singularity branch is delegated to
                        `exp_screw` (which handles small-|u| via Taylor
                        per `small_angle_series.h`).

ERROR BOUND
  Bound category:       precision (from exp_screw small-angle Taylor +
                        floating-point ops) and precision (from
                        retraction-discarded normal component, REQ-EF-15)
  Bound formula:        Two contributions, both inherited / framework-
                        wired, not added directly in this function:
                          1. exp_screw's small-angle Taylor truncation
                             bound (audited under
                             `small_angle_series.md` §§ taylor_sinc and
                             taylor_cos_minus_sinc_over_theta_sq);
                             added to precision inside exp_screw.
                          2. Retraction error per REQ-EF-15: the
                             component projected off the unit-DQ
                             surface during `normalized()`,
                             added to precision inside
                             `DualQuaternion::normalized()`.
                        No model-truncation contribution: a constant
                        twist over the stage is the *exact* one-
                        parameter flow of dM/dt = M · ξ̂. The model
                        error from treating ξ as constant across the
                        stage is realized at the integrator level
                        (next two cards), not here.
  Bound implemented:    None added in `lie_advance_pose` itself. The
                        function composes operations that each add
                        their own bound: `*` (closed-form DQ product),
                        `exp_screw` (Taylor truncation in small-angle
                        branch + float-op precision), `normalized()`
                        (REQ-EF-15 retraction precision).
  Bound verdict:        ✓ matched by composition — no bound is owed
                        at this level beyond what its subcalls add.
                        Conforms to REQ-EF-12 (composite types inherit
                        propagation via composition).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-12 (inheritance by composition),
                        REQ-EF-15 (retraction step is delegated to
                        `Pose::normalized()`).
  AUD-EF applies:       AUD-EF-1 (TrackedValue surfacing — composite
                        Pose<T> wraps tracked components throughout).
  AUD-MC applies:       AUD-MC-18 (screw exp/log round-trip, used in
                        the unit tests that exercise this advance).
  Verification test:    tests/test_integrators/ — RK4 LEO smoke test
                        (one orbit, ω-only constant twist verifies the
                        exact rigid rotation along a single stage).
                        Also exercised indirectly by every Euler / RK4
                        test in `tests/test_integrators`.

NOTES
  - The half-factor `δt/2` in (u, v) is **not** a Lie-group truncation;
    it is the half-angle parameterization required by
    `exp_screw`. A full-step rotation by angle θ = δt · |ω|
    corresponds to u = (δt/2) · ω with |u| = θ/2. Misinterpreting
    this as a Strang-splitting half-step would be a method-theory
    mismatch. Verified against REQ-DQ-18 doc.
  - The retraction `normalized()` is invoked **at every stage**,
    not just post-step. This is conservative: it tightens the
    REQ-EF-15 retraction precision contribution each call but does
    not change the leading-order truncation behavior. A future
    optimization could defer retraction to the final stage only.
  - `lie_advance_pose` is a pure helper — it does not own any global
    truncation order; it inherits the truncation order of `exp_screw`
    (which is exact except in the small-angle Taylor branch).
```

---

## 2 `euler`

```
=== FORMULA AUDIT CARD ===
ID:                     runge_kutta::euler
Location:               src/integrators/runge_kutta.h:68-78
Mathematical statement: One-stage Lie-Euler step for the coupled
                        state ODE
                          dM̂/dt = (1/2) M̂ · Ω̂_pure   (SE(3) Lie group)
                          dΩ̂/dt = a_body              (vector space)
                          dt/dt = 1.
                        Update:
                          M̂_{n+1} = (M̂_n · exp_screw(δt · ω/2,
                                                       δt · v/2))_↾SE(3)
                          Ω̂_{n+1} = Ω̂_n + δt · a(y_n)
                          t_{n+1} = t_n + δt
                        where (ω, v) = (Ω̂_n.angular, Ω̂_n.linear) and
                        a(·) = accel_fn(·).

THEORY
  Underlying theorem:   Euler's method, in its Lie-group formulation
                        (Munthe-Kaas 1999, §2). For the vector-space
                        leg this is the classical forward-Euler scheme
                        from Euler (1768), Institutionum Calculi
                        Integralis Vol. I §650. For the Lie-group leg
                        it is the lowest-order Munthe-Kaas update with
                        single-stage Butcher tableau c = (0),
                        b = (1).
                        Local truncation: O(δt²) per step, hence
                        global (after N = T/δt steps) O(δt) — i.e.
                        first-order accurate (REQ-EF-7 model
                        truncation contribution to accuracy).
                        Convergence and order are theorems of Butcher
                        (2003), "Numerical Methods for Ordinary
                        Differential Equations" §23 (RK on manifolds:
                        §38).
  Primary reference:    Euler (1768) — classical scheme. Butcher
                        (2003), Thm 23.1 (consistency ⇒ order-1
                        convergence). Lie-group extension:
                        Munthe-Kaas (1999), Thm 2.1.
  Domain of validity:   ODE right-hand side is locally Lipschitz on
                        a neighborhood of the trajectory. Step size
                        δt < δt_max for stability of the implicit
                        SE(3) flow (i.e. the half-angle ball
                        |δt · ω| / 2 ≤ π/2 ⇒ |ω| δt ≤ π).

METHOD
  Method declared:      One-stage Lie-Euler (Butcher c = (0), b = (1)).
                        Pose: exact constant-twist screw flow over
                        [t_n, t_n + δt] using twist evaluated at y_n.
                        Twist: forward-Euler in vector space using
                        acceleration evaluated at y_n.
  Method implemented:   accel = accel_fn(y0);
                        pose_new = lie_advance_pose(y0.pose, dt, y0.twist);
                        twist_new = y0.twist + dt · accel;
                        time_new  = y0.time  + dt;
                        — verbatim single-stage forward Euler.
                        No higher-order corrections; no implicit
                        solve; no extrapolation.
  Match verdict:        ✓ matched. Implementation is one-stage Lie-
                        Euler; the cited theory is Lie-Euler. No
                        Taylor expansion in δt is performed
                        explicitly here — the truncation enters
                        implicitly through the integrator's order
                        (Theorem 23.1 of Butcher).

ERROR BOUND
  Bound category:       accuracy (per REQ-EF-7 — model truncation
                        from the *method order*, not the model;
                        but the local-truncation contribution of an
                        integrator is the same conceptual category:
                        a deliberate fixed-order approximation to
                        the true continuous flow).
                        Precision contributions (Taylor in
                        exp_screw, retraction in normalized(),
                        per-op floating-point) are inherited via
                        composition (REQ-EF-12).
  Bound formula:        Local truncation per step:
                          |y_{n+1} - φ_δt(y_n)| ≤ C · δt² · max|y''|
                        with C the standard Euler error coefficient
                        (1/2 for the leading Taylor remainder of the
                        exact flow). For the vector-space leg this is
                          (δt²/2) · max_{[t_n, t_n+δt]} |a'(t)|;
                        for the SE(3) leg it is the leading
                        commutator term in the BCH expansion of the
                        exact two-stage flow against the one-stage
                        flow with frozen twist:
                          (δt²/2) · [ω(t_n), a(t_n)]_se(3) + O(δt³).
                        Global (final-time) accumulation: O(δt).
  Bound implemented:    ⚠ Not added explicitly in `euler`. The
                        local-truncation bound for the **integrator
                        order** is not currently surfaced as a
                        contribution to `errors.accuracy` inside the
                        integrator step. Per REQ-EF-7, model
                        truncation must be accumulated; an
                        integrator-order truncation is conceptually
                        of the same kind (a deliberate fixed-order
                        approximation to the exact continuous
                        flow). The function therefore relies on the
                        caller (propagator) to assess step-error.
                        Per-op precision (Taylor in `exp_screw`,
                        retraction in `normalized()`, vector-space
                        floating-point) **is** propagated via
                        composition; only the **integrator-order**
                        contribution is missing.
  Bound verdict:        ? open / under-counted — composition handles
                        precision-category leakage, but the
                        accuracy-category contribution from the
                        first-order Euler truncation is not added
                        in this function. This is a documented gap
                        relative to REQ-EF-7's intent (model
                        truncation → accuracy). One of two
                        resolutions is needed:
                          (a) integrator step adds a step-error
                              estimate (e.g. embedded RK or
                              Richardson extrapolation) into
                              errors.accuracy; or
                          (b) project-level decision that
                              integrator-order truncation is
                              the propagator's responsibility,
                              not the integrator's, with an
                              explicit step at the propagator
                              layer.

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-7 (integrator-order truncation →
                        accuracy — current status: under-counted),
                        REQ-EF-12 (composition propagation —
                        satisfied), REQ-EF-15 (retraction inside
                        `lie_advance_pose` — satisfied by delegation).
  AUD-EF applies:       AUD-EF-1 (TrackedValue surfacing — composite
                        State<T> wraps tracked components).
                        AUD-EF-7 (model-truncation contribution to
                        accuracy — current integrator-order
                        contribution is the gap).
  AUD-MC applies:       n/a (integrator step, not algebra op).
  Verification test:    tests/test_integrators/ — Euler test with a
                        constant-acceleration synthetic problem,
                        ω = 0, recovers linear vector update exactly;
                        SE(3) leg with constant twist recovers the
                        exact rigid motion.

NOTES
  - The truncation gap is **not** a correctness failure of the Euler
    formula itself — the formula is the standard Lie-Euler step. It
    is a gap in **error reporting**: the propagator currently relies
    on numerical comparison against high-precision reference rather
    than reporting an a-priori bound from the integrator's order.
  - For a first-order method, the simplest a-priori bound that could
    be added without an embedded variant is the *step-doubling
    Richardson estimate*: |y_{δt} - y_{δt/2}| / (2^1 - 1). This
    requires a doubled-cost evaluation per step.
  - Function header lists conformance to `AUD-EF-7`, which suggests
    the gap was previously identified but the bound has not yet been
    wired through. Flag for separate ticket; not in scope for this
    audit card to fix.
```

---

## 3 `runge_kutta_4`

```
=== FORMULA AUDIT CARD ===
ID:                     runge_kutta::runge_kutta_4
Location:               src/integrators/runge_kutta.h:87-137
Mathematical statement: Four-stage Runge–Kutta step for the coupled
                        state ODE (same as Euler card above), using
                        Butcher tableau
                          c = (0, 1/2, 1/2, 1),
                          a = strictly lower-triangular with
                              a_{2,1} = 1/2, a_{3,2} = 1/2, a_{4,3} = 1,
                          b = (1/6, 1/3, 1/3, 1/6).
                        At each stage k, the body twist w_k and body
                        acceleration a_k are evaluated. Final state:
                          twist_avg = (1/6)(w_1 + 2 w_2 + 2 w_3 + w_4)
                          accel_avg = (1/6)(a_1 + 2 a_2 + 2 a_3 + a_4)
                          M_{n+1}    = (M_n · exp_screw(δt·twist_avg/2))_↾SE(3)
                          Ω_{n+1}    = Ω_n + δt · accel_avg
                          t_{n+1}    = t_n + δt

THEORY
  Underlying theorem:   Classical Runge–Kutta order-4 method (Kutta
                        1901), applied stage-by-stage with the Lie-
                        group Munthe-Kaas adaptation for the pose
                        leg.
                          - Vector-space leg: Kutta (1901),
                            "Beitrag zur näherungsweisen Integration
                            totaler Differentialgleichungen", Z.
                            Math. Phys. 46, 435-453. The "classical
                            RK4" Butcher tableau used here.
                          - Lie-group leg: Munthe-Kaas (1999),
                            "High order Runge–Kutta methods on
                            manifolds", Appl. Numer. Math. 29,
                            115-127, §3 (RK4 on Lie groups). Each
                            stage advances the pose by the *exact*
                            constant-twist screw flow over the stage
                            interval using the stage's body twist;
                            the final combination uses a weighted
                            average twist applied as a single full-
                            step screw.
                        Local truncation: O(δt^5) per step, hence
                        global (after T/δt steps) O(δt^4) —
                        fourth-order accurate. Theorem: Butcher
                        (2003), "Numerical Methods for ODEs", §32
                        (classical RK4) and §38 (RK on manifolds).
                        Munthe-Kaas (1999), Thm 3.4 (the standard RK4
                        scheme retains order 4 on Lie groups when
                        the exponential map is exact).
  Primary reference:    Runge (1895), "Ueber die numerische
                        Auflösung von Differentialgleichungen",
                        Math. Ann. 46, 167-178 (origin of the
                        general method). Kutta (1901) for the
                        specific 4-stage tableau. Butcher (2003)
                        §32-38 for the modern unified treatment.
                        Lie-group extension: Munthe-Kaas (1999).
  Domain of validity:   Same as Euler: locally Lipschitz RHS,
                        step-size respects the SE(3) half-angle
                        constraint at the largest-magnitude stage
                        twist (i.e. δt · max_k |w_k| / 2 ≤ π/2).

METHOD
  Method declared:      Classical (non-embedded) RK4 with Butcher
                        tableau c = (0, 1/2, 1/2, 1),
                        b = (1/6, 1/3, 1/3, 1/6). For the SE(3) pose
                        leg: Munthe-Kaas style — each stage advances
                        y_0 from t_n to t_n + c_k · δt using the
                        prior stage's body twist for the pose
                        advance and the prior stage's body
                        acceleration for the twist advance.
                        Combination step: weighted average twist
                        feeds a single full-δt exp_screw, then
                        retraction.
  Method implemented:   Four-stage evaluation matching the declared
                        tableau exactly:
                          Stage 1 (c_1 = 0): w_1 = y0.twist;
                                            a_1 = accel_fn(y0).
                          Stage 2 (c_2 = 1/2): y2.pose via
                                            lie_advance_pose(y0.pose,
                                            half_dt, w_1);
                                            y2.twist = y0.twist
                                                       + half_dt · a_1;
                                            w_2 = y2.twist;
                                            a_2 = accel_fn(y2).
                          Stage 3 (c_3 = 1/2): y3.pose via
                                            lie_advance_pose(y0.pose,
                                            half_dt, w_2);
                                            y3.twist = y0.twist
                                                       + half_dt · a_2;
                                            w_3 = y3.twist;
                                            a_3 = accel_fn(y3).
                          Stage 4 (c_4 = 1):   y4.pose via
                                            lie_advance_pose(y0.pose,
                                            dt, w_3);
                                            y4.twist = y0.twist
                                                       + dt · a_3;
                                            w_4 = y4.twist;
                                            a_4 = accel_fn(y4).
                        Combination:
                          twist_avg = (1/6) · (w_1 + 2 w_2
                                                + 2 w_3 + w_4);
                          accel_avg = (1/6) · (a_1 + 2 a_2
                                                + 2 a_3 + a_4);
                          state_new = (lie_advance_pose(y0.pose, dt,
                                                       twist_avg),
                                       y0.twist + dt · accel_avg,
                                       y0.time + dt).
                        No Padé, no continued fraction, no series in
                        δt — closed-form Butcher RK4.
  Match verdict:        ⚠ method-vs-theory subtlety. The vector-space
                        leg (twist update via accel_avg) is the
                        textbook RK4 closed-form. The pose leg uses
                        **stagewise exp_screw + final exp_screw with
                        averaged twist**, which is one of several
                        Munthe-Kaas RK4 variants. The "classical
                        Munthe-Kaas RK4" as written in Munthe-Kaas
                        (1999) §3 typically uses *exponentials of
                        weighted-averaged Lie-algebra elements*
                        with a more subtle commutator correction
                        (the Magnus / dexpinv terms) to preserve
                        order 4 on a non-Abelian group. This
                        implementation applies a single weighted-
                        average twist through one exp_screw at the
                        end, which is rigorously order 4 only when
                        the body-twist commutators
                        [w_i, w_j]_se(3) are small (e.g. when the
                        rotational rate is moderate relative to
                        δt; specifically for orbital problems with
                        |ω| · δt ≪ 1, the omitted commutator terms
                        contribute O(δt^5) and the order is
                        retained, but the leading coefficient
                        differs from the strict Munthe-Kaas RK4
                        treatment). **Flag for clarification.**
                        Either:
                          (a) Cite this as a *Lie-Euler-on-pose +
                              RK4-on-twist split* with an explicit
                              order-reduction analysis showing
                              order 4 is retained under the orbital
                              parameter regime; or
                          (b) Upgrade to the full Munthe-Kaas RK4
                              with dexpinv correction per
                              Munthe-Kaas (1999) §3.

ERROR BOUND
  Bound category:       accuracy (integrator-order truncation, same
                        category mapping as Euler card above) +
                        precision (composition of `exp_screw`
                        Taylor branch + per-op floats + retraction
                        contributions from each `lie_advance_pose`
                        call).
  Bound formula:        Local truncation per step (Butcher 2003,
                        Thm 32.4):
                          |y_{n+1} - φ_δt(y_n)|
                              ≤ C_RK4 · δt^5 · max_{[t_n,t_n+δt]} |y^{(5)}|
                        with C_RK4 the RK4 leading constant
                        (1/720 for the scalar Taylor remainder;
                        Lie-group variant carries an additional
                        BCH-commutator contribution of the same
                        order, see Munthe-Kaas Thm 3.4 remarks).
                        Global accumulation over N = T/δt steps:
                        O(δt^4).
  Bound implemented:    ⚠ Not added explicitly in `runge_kutta_4`.
                        Same situation as Euler: the integrator-order
                        accuracy contribution is not surfaced into
                        `state.errors.accuracy`. Per-op precision
                        contributions (Taylor inside exp_screw,
                        retraction inside normalized()) **are**
                        inherited via composition. The integrator-
                        order accuracy bound is a documented gap.
  Bound verdict:        ? open / under-counted — same C-status as
                        Euler. The formula is correct, the
                        composition propagates per-op precision, but
                        the integrator-order accuracy contribution
                        is not currently wired into the state's
                        accuracy budget per REQ-EF-7. Resolutions
                        identical to Euler card (embedded RK,
                        Richardson, or propagator-layer
                        responsibility).

CROSS-AUDIT
  REQ-EF applies:       REQ-EF-7 (integrator-order truncation →
                        accuracy — under-counted, same as Euler),
                        REQ-EF-12 (composition propagation —
                        satisfied), REQ-EF-15 (retraction is
                        delegated to `lie_advance_pose` →
                        `Pose::normalized()` — satisfied).
  AUD-EF applies:       AUD-EF-1 (TrackedValue surfacing through
                        composite State<T>), AUD-EF-7 (integrator-
                        order accuracy contribution — gap).
  AUD-MC applies:       n/a (integrator step, not algebra op).
                        Indirect: AUD-MC-18 (screw exp round-trip)
                        is exercised by every pose advance.
  Verification test:    tests/test_integrators/ — LEO smoke test
                        (one orbit, RK4 vs analytic two-body
                        agreement to ~11/12 components per Phase
                        D integration). High-precision reference
                        propagation gives the **observed** error;
                        the audit asks for an **a-priori** bound.

NOTES
  - The Munthe-Kaas-variant subtlety (Match verdict ⚠) is a real
    theory-method open question. Empirical evidence from the LEO
    smoke test suggests the implementation retains order 4 in
    practice for orbital-rate parameters; the question is whether
    the cited theorem (Munthe-Kaas Thm 3.4) **applies as stated**
    or whether an order-reduction analysis under the orbital
    regime is needed. The latter is a derivation gap, not a code
    bug.
  - The composite y_k state passed to `accel_fn` reads `time + c_k δt`
    via `y0.time + half_dt` / `y0.time + dt`, matching the tableau's
    c-row. ✓ No time-mis-evaluation bug.
  - The combination uses the **average twist** w̄ = (w_1 + 2w_2 +
    2w_3 + w_4)/6 fed to a single exp_screw, not the average of
    the four exp_screws (which would require BCH to combine). The
    chosen approach is a particular Munthe-Kaas variant, sometimes
    called "RKMK-Lie-Euler-composite" in the manifold-integration
    literature; it is **not** the most general RKMK form.
  - Header lists conformance to AUD-EF-7. The gap is the same as
    Euler: missing integrator-order accuracy bound. Flag for
    separate ticket.
```

---

## 4 File-level verdict for `runge_kutta.h`

| Function | Theory match | Method match | Bound match | Verdict |
|---|---|---|---|---|
| `lie_advance_pose` | ✓ Munthe-Kaas constant-twist exact flow + REQ-EF-15 retraction | ✓ closed-form exp_screw + normalized() | ✓ (by composition; no bound owed at this level) | **PASS** |
| `euler` | ✓ Lie-Euler (Euler 1768, Munthe-Kaas 1999) | ✓ one-stage forward Euler | ? integrator-order accuracy not added (REQ-EF-7 gap) | **PASS-with-note** |
| `runge_kutta_4` | ⚠ Munthe-Kaas-variant order-4 question (classical RK4 on twist; composite-exp form on pose) | ✓ Butcher c = (0, 1/2, 1/2, 1), b = (1/6, 1/3, 1/3, 1/6) | ? integrator-order accuracy not added (REQ-EF-7 gap) | **PASS-with-flags** |

**File verdict: PASS-with-flags.** No method-theory mismatch of the "Taylor cited, continued fraction implemented" kind that the framework targets. Two open items to track:

1. **REQ-EF-7 wiring** (applies to both `euler` and `runge_kutta_4`): integrator-order local-truncation bound is not currently surfaced into `errors.accuracy` at the integrator step. Resolution options listed per card; project-level decision required on whether this belongs at the integrator or propagator layer.

2. **Munthe-Kaas RK4 variant** (applies to `runge_kutta_4`): the pose leg uses a composite-exp form (single exp_screw of weighted-average twist) rather than the strict Munthe-Kaas RK4 with dexpinv correction. Empirically order 4 holds in the orbital regime; the theory citation should either pin down the variant (with an order-retention proof for orbital parameters) or upgrade the implementation to the strict Munthe-Kaas form.

Neither flag invalidates the current `total_error()` claim for **per-op precision** propagation (composition handles that). Both flags affect the **accuracy** category's REQ-EF-7 contribution, which is the integrator's claim about its own truncation order.

---

## 5 References

- Euler, L. (1768). *Institutionum Calculi Integralis*, Vol. I, §650.
- Runge, C. (1895). "Ueber die numerische Auflösung von Differentialgleichungen", *Math. Ann.* **46**, 167-178.
- Kutta, W. (1901). "Beitrag zur näherungsweisen Integration totaler Differentialgleichungen", *Z. Math. Phys.* **46**, 435-453.
- Munthe-Kaas, H. (1999). "High order Runge–Kutta methods on manifolds", *Appl. Numer. Math.* **29**, 115-127.
- Iserles, A., Munthe-Kaas, H., Nørsett, S., Zanna, A. (2000). "Lie-group methods", *Acta Numerica* **9**, 215-365.
- Butcher, J. C. (2003). *Numerical Methods for Ordinary Differential Equations*, 2nd ed., Wiley.
- Selig, J. M. (2005). *Geometric Fundamentals of Robotics*, 2nd ed., Springer.
- Kavan, L., Collins, S., Žára, J., O'Sullivan, C. (2008). "Geometric skinning with approximate dual quaternion blending", *ACM Trans. Graph.* **27**(4).

Project anchors:
- `design/audit/theoretical_basis_audit.md` (framework, §1 audit-card slots, §5 worked example).
- `design/specifications/error_framework.md` — REQ-EF-7, REQ-EF-12, REQ-EF-15.
- `src/math/dual_quaternion.h` — REQ-DQ-18 half-angle exp_screw encoding.
- `src/math/small_angle_series.h` — Taylor branches consumed indirectly via exp_screw.
