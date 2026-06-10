### Theorem 0.3.4 (Gauss variational equations for `Ω̇` and `i̇`)

**Hypotheses.** As in Theorem 0.3.2.

**Conclusion.** Let `u := ω + f` be the argument of latitude. Then

```
Ω̇ = (r N sin u) / (n a² β sin i) ,                                  (0.3.4.1)
i̇ = (r N cos u) / (n a² β) .                                         (0.3.4.2)
```

In particular, both vanish when `N = 0` (purely in-plane perturbation).

**Proof.**

**Step 1 (Differentiate the angular momentum vector).** From `h_vec := r × v`,

```
dh_vec/dt = ṙ × v + r × v̇ = v × v + r × (−μ r̂/r² + F)
          = r × F                                                    (0.3.4.3)
```

(since `v × v = 0` and `r × (−μ r̂/r²) = (−μ/r²)(r × r̂) = 0` because `r ∥ r̂`).

**Step 2 (Component-decomposition of `r × F`).** With `r = r r̂` and
`F = R r̂ + T t̂ + N n̂`:

```
r × F = r r̂ × (R r̂ + T t̂ + N n̂)
      = r [R (r̂ × r̂) + T (r̂ × t̂) + N (r̂ × n̂)]
      = r [0 + T n̂ + N (−t̂)]
      = r T n̂ − r N t̂ .                                            (0.3.4.4)
```

**Step 3 (Split `dh_vec/dt` into magnitude and direction).** Write
`h_vec = h n̂` with `h = ‖h_vec‖` the scalar magnitude and `n̂` the unit orbit
normal. Differentiating:

```
dh_vec/dt = (dh/dt) n̂ + h (dn̂/dt) .                                (0.3.4.5)
```

Combining (0.3.4.3-4) and (0.3.4.5), and recalling that `dn̂/dt ⊥ n̂` (any unit
vector's derivative is perpendicular to itself):

- **n̂ component of (0.3.4.5)**: only `(dh/dt) n̂` contributes, so
  `dh/dt = r T` (matches Theorem 0.2.3 differentiated under perturbation).
- **t̂ component of (0.3.4.5)**: only `h (dn̂/dt)` contributes, equal to
  `−r N t̂` from (0.3.4.4). Therefore

```
h (dn̂/dt) = −r N t̂                                                   (0.3.4.6)
```

(the orbit-plane rotation vector lies along `t̂` ⊥ the line of nodes — see
Step 5 for the geometric decomposition).

**Step 4 (Inertial-frame parameterization of `n̂`).** Choose inertial Cartesian
axes `(x̂, ŷ, ẑ)` with `ẑ` along the celestial north pole. The orbital plane is
specified by `(i, Ω)`:

```
n̂(i, Ω) = sin i · sin Ω · x̂  −  sin i · cos Ω · ŷ  +  cos i · ẑ .  (0.3.4.7)
```

The line of nodes is `N̂ := cos Ω · x̂ + sin Ω · ŷ` (the in-equator direction
along which the orbital plane intersects the equator). The "out-of-node" in-plane
direction is

```
M̂ := cos i · (−sin Ω x̂ + cos Ω ŷ) + sin i · ẑ ,                    (0.3.4.8)
```

and the right-handed orthonormal basis on the orbital plane is `(N̂, M̂, n̂)`.

**Step 5 (Differentiate `n̂` w.r.t. time).** Using (0.3.4.7):

```
∂n̂/∂Ω = (cos Ω) sin i · x̂ + (sin Ω) sin i · ŷ = sin i · N̂ ,        (0.3.4.9)
∂n̂/∂i = cos i sin Ω · x̂ − cos i cos Ω · ŷ − sin i · ẑ = −M̂ .      (0.3.4.10)
```

By the chain rule,

```
dn̂/dt = Ω̇ · (∂n̂/∂Ω) + i̇ · (∂n̂/∂i) = Ω̇ sin i · N̂ − i̇ · M̂ .    (0.3.4.11)
```

**Step 6 (Express `t̂` in the `(N̂, M̂)` basis).** The position direction
`r̂` makes angle `u := ω + f` with the ascending node `N̂` (this is the
**argument of latitude** by definition):

```
r̂ = cos u · N̂ + sin u · M̂ ,                                       (0.3.4.12)
t̂ = −sin u · N̂ + cos u · M̂ .                                      (0.3.4.13)
```

**Step 7 (Match components in (0.3.4.6) using (0.3.4.11) and (0.3.4.13)).**

```
h · (Ω̇ sin i · N̂ − i̇ · M̂) = −r N · (−sin u N̂ + cos u M̂)
                              = r N sin u · N̂ − r N cos u · M̂ .
```

Matching coefficients of `N̂`:

```
h Ω̇ sin i = r N sin u   ⇒   Ω̇ = (r N sin u) / (h sin i) .         (0.3.4.14)
```

Matching coefficients of `M̂`:

```
−h i̇ = −r N cos u   ⇒   i̇ = (r N cos u) / h .                     (0.3.4.15)
```

**Step 8 (Substitute `h = n a² β`).** By Theorem 0.2.4 (0.2.4.2),
`h = n a² β`. Substituting into (0.3.4.14)-(0.3.4.15):

```
Ω̇ = (r N sin u) / (n a² β sin i) ,
i̇ = (r N cos u) / (n a² β) .                                        (0.3.4.16)
```

These match (0.3.4.1)-(0.3.4.2). ∎

**Remark 0.3.4.17 (Singularity at `i = 0`).** The formula (0.3.4.1) for `Ω̇`
diverges as `sin i → 0`. This is a coordinate singularity: for equatorial orbits
the line of nodes is undefined, and `Ω` itself is not a meaningful angle.
Implementations resolve this either by using equinoctial elements (which are
nonsingular at `i = 0`) or by branching code that sets `Ω̇ = 0` (and reassigns
the meaning of `ω` as longitude-of-perihelion) when `sin i` is below a tolerance.
SGP4's `near_space.h` carries the xlcof critical-inclination fallback (R14
bonus fix in the audit) for the related `cos i = −1` singularity in the
long-period periodics; the `sin i = 0` case is handled by the Brouwer secular
rates upstream.

**Alignment to SGP4.**
- (a) **Symbol bridge.** `n → n₀''`, `a → a₀''`. `u = ω + f` is the
  **argument of latitude**, used in J₂ short-period and J₃ long-period
  derivations downstream.
- (b) **For drag (`N = 0`).** Both `Ω̇_drag = 0` and `i̇_drag = 0` identically.
  Drag does not change the orbital plane orientation — only the in-plane
  elements `(a, e, ω, M)` decay.
- (c) **For J₃ short-period (`N ≠ 0`).** The J₃ zonal harmonic produces an
  out-of-plane radial perturbation `δr_{J_3} ∝ sin i sin u` (Brouwer 1959
  Eq. 19, derived in the BH61 cleanroom). This couples to drag via the density
  modification (`δρ = (∂ρ/∂r) δr_{J_3}`) and, separately, contributes to the
  orbital-plane drift via Theorem 0.3.4. The drift contribution is one of the
  reasons the long-period `xlcof / aycof` corrections in §14 are needed.
- (d) **Coordinate singularity (i → 0).** SGP4 invariants for low-inclination
  orbits are handled by `near_space.h`. Phase 0 derivations remain valid for
  `i > 0`; the Lyddane Lyddane-edge SGP4 failures noted in the audit (sats
  21897, 22674, etc.) are downstream-of-Phase-0 coordinate-handling issues.

**Alignment to implementation (deferred per Standard 9-B).** None at Phase 0.

---

### Lemma 0.3.5.a (Rotation of the line of nodes within the orbital plane)

The line of nodes `N̂(Ω) = cos Ω x̂ + sin Ω ŷ` lies in the equatorial plane and
rotates around the inertial `ẑ`-axis at rate `Ω̇`. Its projection onto the
orbital plane has angular velocity, with respect to the orbit-normal `n̂`,

```
(dN̂/dt) · M̂ = Ω̇ cos i .                                            (0.3.5.a.1)
```

In words: the line of nodes drifts through the orbital plane at rate `Ω̇ cos i`
around `n̂`.

**Proof.** Differentiating `N̂ = cos Ω x̂ + sin Ω ŷ` w.r.t. time:

```
dN̂/dt = Ω̇ (−sin Ω x̂ + cos Ω ŷ) .                                  (0.3.5.a.2)
```

Project onto `M̂` from (0.3.4.8):

```
M̂ · (dN̂/dt) = [cos i (−sin Ω x̂ + cos Ω ŷ) + sin i ẑ] · Ω̇ (−sin Ω x̂ + cos Ω ŷ)
            = Ω̇ cos i [(−sin Ω)(−sin Ω) + (cos Ω)(cos Ω)] + 0
            = Ω̇ cos i [sin²Ω + cos²Ω]
            = Ω̇ cos i .                                              (0.3.5.a.3)
```

∎

### Theorem 0.3.5 (Gauss variational equation for `ω̇`)

**Hypotheses.** As in Theorem 0.3.2.

**Conclusion.**

```
ω̇ = (β / (n a e)) · [−R cos f + T sin f (1 + r/p)] − Ω̇ cos i ,    (0.3.5.1)
```

with `Ω̇` given by Theorem 0.3.4 (0.3.4.1). For purely in-plane perturbations
(`N = 0`), `Ω̇ = 0`, so

```
ω̇\|_{N=0} = (β / (n a e)) · [−R cos f + T sin f (1 + r/p)] .         (0.3.5.2)
```

**Proof.**

**Step 1 (Geometric setup).** The eccentricity vector `e_vec = e P̂` has
magnitude `e` and direction `P̂` (toward perigee). Within the orbital plane,
`P̂` makes angle `ω` with the ascending node `N̂`. The argument of perigee `ω`
is therefore defined as

```
ω := ∠(N̂, P̂)                                                       (0.3.5.3)
```

measured in the orbital plane (right-handed around `n̂`).

**Step 2 (Decompose `ω̇` into intrinsic and nodal parts).** The rate of change
of `ω` decomposes as

```
ω̇ = (rate of P̂ around n̂, inertial) − (rate of N̂ around n̂, inertial) . (0.3.5.4)
```

The first term is the **inertial rotation rate of `e_vec` around `n̂`**, equal to
`(de_vec/dt · Q̂)/e` for `Q̂ := n̂ × P̂` (since `de_vec/dt = ė P̂ + e (rotation of P̂)`,
and the `P̂`-component contributes `ė` not rotation; the `Q̂`-component
contributes the rotation rate divided by `e`).

The second term, **by Lemma 0.3.5.a**, equals `Ω̇ cos i`. So

```
ω̇ = (de_vec/dt · Q̂) / e − Ω̇ cos i .                              (0.3.5.5)
```

**Step 3 (Compute `de_vec/dt · Q̂` from 0.3.3.10).** Recall (0.3.3.10):

```
de_vec/dt = (1/μ) [2 h T r̂ − (h R + r ṙ T) t̂ − r ṙ N n̂] .         (0.3.5.6)
```

Project onto `Q̂` using `r̂ · Q̂ = sin f`, `t̂ · Q̂ = cos f`, `n̂ · Q̂ = 0`:

```
de_vec/dt · Q̂ = (1/μ) [2 h T sin f − (h R + r ṙ T) cos f − 0]
              = (1/μ) [−h R cos f + T (2 h sin f − r ṙ cos f)] .   (0.3.5.7)
```

**Step 4 (Simplify the `T` bracket using `1 + r/p` identity).** From (0.3.3.13)
`r ṙ = e h sin f / (1 + e cos f)`. So

```
2 h sin f − r ṙ cos f = 2 h sin f − (e h sin f cos f) / (1 + e cos f)
                     = (h sin f / (1 + e cos f)) · [2 (1 + e cos f) − e cos f]
                     = (h sin f / (1 + e cos f)) · [2 + e cos f] .  (0.3.5.8)
```

Now observe that `2 + e cos f = 1 + (1 + e cos f)`, so

```
(2 + e cos f) / (1 + e cos f) = 1 + 1/(1 + e cos f) = 1 + r/p ,    (0.3.5.9)
```

using `r/p = 1/(1 + e cos f)` from Theorem 0.2.2. Therefore

```
2 h sin f − r ṙ cos f = h sin f · (1 + r/p) .                       (0.3.5.10)
```

**Step 5 (Assemble the `de_vec/dt · Q̂` formula).** Substitute (0.3.5.10) into
(0.3.5.7):

```
de_vec/dt · Q̂ = (1/μ) [−h R cos f + T h sin f (1 + r/p)]
              = (h/μ) [−R cos f + T sin f (1 + r/p)] .              (0.3.5.11)
```

By Theorem 0.2.4 and Kepler 3rd, `h/μ = (n a² β)/(n² a³) = β/(n a)`. Substituting:

```
de_vec/dt · Q̂ = (β/(n a)) · [−R cos f + T sin f (1 + r/p)] .       (0.3.5.12)
```

**Step 6 (Final form for `ω̇`).** Substitute (0.3.5.12) into (0.3.5.5):

```
ω̇ = (β/(n a e)) · [−R cos f + T sin f (1 + r/p)] − Ω̇ cos i .     (0.3.5.13)
```

This is (0.3.5.1). The drag specialization (0.3.5.2) follows by setting `N = 0`
(which gives `Ω̇ = 0` by Theorem 0.3.4 (0.3.4.1)). ∎

**Remark 0.3.5.14 (Singularity at `e → 0`).** Like (0.3.4.1) at `i = 0`, the
formula (0.3.5.1) for `ω̇` diverges as `e → 0`. This is again a coordinate
singularity: for circular orbits the perigee is undefined and `ω` is not
meaningful. SGP4's `drag_coefficients.h` branches on `e₀ > 10⁻⁴` (audited at
R03 / R12), zeroing the `C₃ / xmcof` coefficients below that threshold. The
analytical `1/e` divergence is exactly the structural reason for the threshold;
see also the **N**-component contribution analysis in Phase 3 (C₃ derivation).

**Alignment to SGP4.**
- (a) **Symbol bridge.** `n → n₀''`, `a → a₀''`, `e → e₀`. For drag (Phase 4 / 10
  derivations), `ω̇` is the drag-driven perigee-rotation rate. For J₂ secular
  drift, `ω̇_{J_2}` is handled by the Brouwer secular generator (BH61 cleanroom)
  and added separately at propagation time.
- (b) **For drag (`N = 0`).** The `−Ω̇ cos i` correction vanishes. (0.3.5.2)
  applies directly. The drag-driven perigee rate is `O(B*)` and appears in §10
  / Phase 10 (omgcof derivation).
- (c) **For J₃ short-period (`N ≠ 0`).** The Ω̇ correction is nonzero. It is
  absorbed into the long-period J₃ generator chain (xlcof / aycof / omgcof,
  Phase 7).
- (d) **Eccentricity-vector convention.** SGP4 long-period periodics
  (xlcof / aycof in §14) act on `e_vec` components `(a_xN, a_yN)`, not on
  scalar `(e, ω)`. The eccentricity-vector formulation in Lemma 0.3.3.a +
  Theorem 0.3.5 is what enables that direct connection.

**Alignment to implementation (deferred per Standard 9-B).** None at Phase 0.

---

### Theorem 0.3.6 (Gauss variational equation for `Ṁ`)

**Hypotheses.** As in Theorem 0.3.2.

**Conclusion.**

```
Ṁ = n − (2 r R) / (n a²) + (β² / (n a e)) · [R cos f − T sin f (1 + r/p)] . (0.3.6.1)
```

Equivalently, collecting the `R`-terms:

```
Ṁ = n − (β² / (n a e)) · [R (2 e r/p − cos f) + T sin f (1 + r/p)] .   (0.3.6.2)
```

The unperturbed limit `R = T = 0` recovers `Ṁ = n` (Kepler).

**Proof.**

**Step 1 (Differentiate Kepler's equation).** From `M = E − e sin E`:

```
dM/dt = (1 − e cos E) dE/dt − sin E · ė .                            (0.3.6.3)
```

Using `1 − e cos E = r/a` (this is Theorem 0.2.2 in its `(a, e, E)` form,
`r = a(1 − e cos E)`):

```
dM/dt = (r/a) · dE/dt − sin E · ė .                                  (0.3.6.4)
```

**Step 2 (Solve for `dE/dt` from `r = a(1 − e cos E)`).** Differentiate:

```
ṙ = ȧ (1 − e cos E) − a ė cos E + a e sin E · dE/dt .                (0.3.6.5)
```

Solve for `dE/dt`:

```
dE/dt = (ṙ − ȧ (r/a) + a ė cos E) / (a e sin E)                      (0.3.6.6)
```

(using `1 − e cos E = r/a` from Step 1).

**Step 3 (Substitute `dE/dt` into `Ṁ` formula).** Combining (0.3.6.4) and
(0.3.6.6):

```
Ṁ = (r/a) · (ṙ − ȧ (r/a) + a ė cos E) / (a e sin E) − ė sin E
   = (r ṙ − ȧ r²/a + r a ė cos E) / (a² e sin E) − ė sin E
   = (r ṙ) / (a² e sin E)   ⎫
     − (ȧ r²) / (a³ e sin E) ⎬  three terms                        (0.3.6.7)
     + ė · [(r cos E) / (a e sin E) − sin E] ⎭
```

The three terms are evaluated in Steps 4–6.

**Step 4 (Term 1: `(r ṙ)/(a² e sin E) = n`).** Substitute the closed forms:
- `ṙ = (n a e sin f) / β` from (0.3.2.10);
- `r = a β² / (1 + e cos f)` from Theorem 0.2.2;
- `sin E = β sin f / (1 + e cos f)`, the **eccentric-true-anomaly identity**
  (which follows from `cos E = (e + cos f)/(1 + e cos f)` and `sin²E + cos²E = 1`;
  cite [WIKI-KE]).

Compute:

```
r ṙ = (a β² / (1 + e cos f)) · (n a e sin f / β) = (n a² β e sin f) / (1 + e cos f) .
a² e sin E = a² e · (β sin f / (1 + e cos f)) .

(r ṙ) / (a² e sin E) = [(n a² β e sin f)/(1 + e cos f)] / [(a² e β sin f)/(1 + e cos f)] = n .
                                                                       (0.3.6.8)
```

This is the unperturbed Kepler rate, as expected.

**Step 5 (Term 2: `(ȧ r²)/(a³ e sin E)`).** Compute `r²/(a³ e sin E)`:

```
r²/a³ = β⁴ / (a (1 + e cos f)²) ,    sin E = β sin f / (1 + e cos f) ,
r²/(a³ sin E) = (β⁴ / (a (1 + e cos f)²)) · ((1 + e cos f) / (β sin f))
            = β³ / (a (1 + e cos f) sin f) .                          (0.3.6.9)
```

Then `(ȧ r²) / (a³ e sin E) = ȧ · β³ / (a e (1 + e cos f) sin f)`. Substitute
`ȧ` from Theorem 0.3.2 (0.3.2.15) = `(2 / (n β)) · [e R sin f + T (1 + e cos f)]`:

```
(ȧ r²)/(a³ e sin E) = (2 / (n β)) · [e R sin f + T (1 + e cos f)] · β³ / (a e (1 + e cos f) sin f)
                    = (2 β² / (n a e (1 + e cos f) sin f)) · [e R sin f + T (1 + e cos f)] .
                                                                        (0.3.6.10)
```

Expand:

```
= (2 β² e R sin f)/(n a e (1+e cos f) sin f) + (2 β² T (1+e cos f))/(n a e (1+e cos f) sin f)
= (2 β² R)/(n a (1+e cos f)) + (2 β² T)/(n a e sin f) .                (0.3.6.11)
```

Note the first piece can be simplified using `(1+e cos f) = p/r`:
`(2 β² R)/(n a · p/r) = 2 β² r R / (n a p) = 2 r R / (n a²)` (using `p = a β²`).
So

```
(ȧ r²)/(a³ e sin E) = (2 r R)/(n a²) + (2 β² T)/(n a e sin f) .       (0.3.6.12)
```

**Step 6 (Term 3: `ė · [(r cos E)/(a e sin E) − sin E]`).** Compute the bracket
factor. We have, from algebraic simplification:

```
r cos E − a e sin² E = (a β² / (1 + e cos f)²) · [(e + cos f) − e sin² f]
                    = (a β² / (1 + e cos f)²) · [e (1 − sin² f) + cos f]
                    = (a β² / (1 + e cos f)²) · [e cos² f + cos f]
                    = (a β² cos f (1 + e cos f)) / (1 + e cos f)²
                    = (a β² cos f) / (1 + e cos f) .                  (0.3.6.13)
```

Divide by `a e sin E = a e · β sin f / (1 + e cos f)`:

```
(r cos E − a e sin² E) / (a e sin E) = [(a β² cos f)/(1+e cos f)] / [(a e β sin f)/(1+e cos f)]
                                     = (β cos f) / (e sin f) .         (0.3.6.14)
```

Now substitute `ė` from Theorem 0.3.3 (0.3.3.18) = `(β / (n a)) · [R sin f + T (cos f + cos E)]`:

```
ė · (β cos f) / (e sin f) = (β / (n a)) · [R sin f + T (cos f + cos E)] · (β cos f) / (e sin f)
                          = (β² cos f) / (n a e sin f) · [R sin f + T (cos f + cos E)]
                          = (β² R cos f) / (n a e) + (β² T cos f (cos f + cos E)) / (n a e sin f) .
                                                                        (0.3.6.15)
```

**Step 7 (Combine Steps 4-6).** Substituting (0.3.6.8), (0.3.6.12), and (0.3.6.15)
into (0.3.6.7):

```
Ṁ = n − [(2 r R)/(n a²) + (2 β² T)/(n a e sin f)]
       + (β² R cos f)/(n a e) + (β² T cos f (cos f + cos E))/(n a e sin f) .
                                                                        (0.3.6.16)
```

**Step 8 (Group R- and T-coefficients).**

```
Ṁ = n + R · [−2 r/(n a²) + β² cos f /(n a e)]
       + T · [−2 β²/(n a e sin f) + β² cos f (cos f + cos E)/(n a e sin f)]
   = n + R · [−2 r/(n a²) + β² cos f /(n a e)]
       + T · (β² /(n a e sin f)) · [−2 + cos f (cos f + cos E)] .       (0.3.6.17)
```

**Step 9 (Simplify the T-bracket using `cos E` identity).** Compute
`−2 + cos f (cos f + cos E)`. From `cos E = (e + cos f)/(1 + e cos f)`:

```
cos f (cos f + cos E) = cos² f + cos f (e + cos f)/(1 + e cos f)
                     = [cos² f (1+e cos f) + e cos f + cos² f] / (1+e cos f)
                     = [2 cos² f + e cos f + e cos³ f] / (1+e cos f) .  (0.3.6.18)
```

So

```
−2 + cos f (cos f + cos E) = [−2 (1+e cos f) + 2 cos² f + e cos f + e cos³ f] / (1+e cos f)
                           = [−2 + 2 cos² f − e cos f + e cos³ f] / (1+e cos f)
                           = [−2 sin² f − e cos f (1 − cos² f)] / (1+e cos f)   ⎫ since 2 cos²f − 2 = −2 sin²f
                           = [−2 sin² f − e cos f · sin² f] / (1+e cos f)        ⎬ using 1 − cos²f = sin²f
                           = −sin² f · (2 + e cos f) / (1+e cos f) .                 ⎭   (0.3.6.19)
```

Use `(2 + e cos f)/(1 + e cos f) = 1 + r/p` from (0.3.5.9):

```
−2 + cos f (cos f + cos E) = −sin² f · (1 + r/p) .                    (0.3.6.20)
```

Therefore the T-coefficient in (0.3.6.17) becomes

```
T · (β² /(n a e sin f)) · (−sin² f · (1 + r/p)) = −T · (β² sin f / (n a e)) · (1 + r/p) .
                                                                        (0.3.6.21)
```

**Step 10 (Final form).** Substituting into (0.3.6.17):

```
Ṁ = n + R · [−2 r/(n a²) + β² cos f /(n a e)] − T · (β² sin f / (n a e)) · (1 + r/p)
   = n − (2 r R)/(n a²) + (β² / (n a e)) · [R cos f − T sin f (1 + r/p)] . (0.3.6.22)
```

This is (0.3.6.1). The equivalent form (0.3.6.2) is obtained by absorbing the
`−2 r R/(n a²)` term into the bracket using `−2 r/(n a²) = −(β²/(n a e)) · (2 e r/p)`
(verified: `(β²·2 e r)/(n a e p) = 2 β² r / (n a p) = 2 r/(n a²)` with `p = a β²`).
∎

**Remark 0.3.6.23 (Why `Ṁ` has more terms than `ȧ, ė, ω̇`).** The mean anomaly
`M` is the only orbital element whose unperturbed rate is nonzero (`= n`). The
perturbation enters at first order in `F` (the linear `R, T` terms) but also
through the perturbation-induced change in `n = √(μ/a³)` itself: as `a` decays,
`n` increases, contributing to `Ṁ` indirectly. The `−2 r R/(n a²)` term in
(0.3.6.1) captures exactly this `dn/dt = (-3/(2a)) n ȧ` contribution after some
algebra (cf. Step 5's manipulation of the `ȧ` term).

This is the structural reason `M` is often replaced in perturbation analyses by
the **mean anomaly at epoch** `l₀ := M − n(t − t₀)`, whose unperturbed rate is
zero. The Brouwer secular formulation in (BH61 cleanroom) operates on
`(l₀, g₀, h₀)` rather than `(M, ω, Ω)` to avoid this complication. Phase 6 of
the SGP4-drag re-derivation (the `D₂` Taylor expansion) explicitly handles the
`n(a(t))` back-reaction; see §10 of the main derivation document (which has
the **D-5 critical finding** that the proof there does not close cleanly — this
Phase 0 result is what makes the Phase 6 re-derivation possible).

**Alignment to SGP4.**
- (a) **Symbol bridge.** `n → n₀''`, `a → a₀''`, `e → e₀`. The Brouwer-recovered
  mean motion `n₀''` is the value that, when used in `Ṁ = n₀''`, reproduces the
  time-averaged actual mean motion including the J₂-secular correction. The
  perturbation correction in (0.3.6.1) is the **additional** drag and J-zonal
  short-period contribution beyond the Brouwer-secular drift.
- (b) **For drag (`N = 0`).** Both `R_drag` and `T_drag` contribute. The drag
  decreases `a` (so `Ṁ` increases via the `−(2 r R)/(n a²)` channel **when
  R_drag has the correct sign**) and modifies the eccentricity / argument-of-
  perigee chain (so `Ṁ` is also affected via the `T` channel). Phase 4 (C₄) and
  Phase 10 (xmcof) extract the secular and `cos M` harmonic components of
  this rate.
- (c) **Singularity at `e → 0`.** As in Theorem 0.3.5, the `1/e` factor in the
  `R cos f`-coefficient is a coordinate singularity. SGP4's `drag_coefficients.h`
  branches at `e₀ > 10⁻⁴`; the rate (0.3.6.1) is well-defined for `e > 0`. The
  drag-modified mean anomaly correction `xmcof` (R12 audit) has the explicit
  `1/(e₀ η)` form that traces to this `1/e` from (0.3.6.1) plus a `1/η` from
  the orbit-average evaluation.
- (d) **Brouwer-recovered subtraction.** When SGP4 stores `M(t₀) = M₀` and
  propagates `M(t) = M₀ + n₀''·(t − t₀) + (drag corrections)`, the
  drag corrections are the integrated form of (0.3.6.1) minus `n₀''`. This is
  the bookkeeping that the Phase 8 t-cofs derivation must close.

**Alignment to implementation (deferred per Standard 9-B).** None at Phase 0.
Implementation choices in `omgcof / xmcof / delmo / sinmo` (R12-audited code)
use specific harmonic decompositions (e.g. `(1 + η cos M)³` cubic-only); those
are Phase 10 considerations.

---

