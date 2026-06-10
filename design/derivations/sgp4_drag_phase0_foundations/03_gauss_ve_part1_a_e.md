## §0.3 Gauss-Lagrange Variational Equations

The Gauss-Lagrange variational equations express the time derivatives of the
orbital elements `(a, e, i, ω, Ω, M)` under a perturbing acceleration `F` (per
unit mass). For SGP4 drag, the perturbation has only a tangential component
(along-track), so the in-plane rates `(ȧ, ė, ω̇, Ṁ)` dominate while `(Ω̇, i̇) = 0`
identically. We derive the full set here so that the framework also covers J₃
out-of-plane couplings used in Phases 7 and 9.

The derivation strategy throughout §0.3 is to **avoid Lagrange-bracket
machinery** (which would require the Hamiltonian formulation) and work directly
from the **conservation-law derivatives** of the unperturbed integrals of motion
under the perturbation. Specifically:

- `ȧ` from `dE/dt = v · F` and `E = -μ/(2a)`.
- `i̇, Ω̇` from `dh/dt = r × F` projected onto the orbit-normal.
- `ė` from `d|e|/dt = d|v × h / μ − r̂|/dt` projected appropriately.
- `ω̇` from the projection of the eccentricity-vector rate onto the in-plane
  direction perpendicular to the eccentricity vector.
- `Ṁ` from the Lagrange relation `dM/dt = n + (correction)` with the correction
  derived from the chain rule on the elements.

This path uses only elementary vector calculus and the Keplerian primitives of
§0.2 — no Hamiltonian / canonical-conjugate / Poisson-bracket machinery and no
OCR-source dependencies.

### Definition 0.3.1 (Radial-transverse-normal decomposition of `F`)

Let `r̂ := r / ‖r‖` be the radial unit vector pointing from the central body to
the satellite. Let `t̂` be the unit vector in the orbital plane, perpendicular to
`r̂` in the direction of increasing true anomaly `f` (the *transverse* direction).
Let `n̂ := r̂ × t̂` be the orbit-normal unit vector.

The triad `(r̂, t̂, n̂)` is right-handed orthonormal. Any perturbing acceleration
`F` (per unit mass) acting on the satellite decomposes as

```
F = R · r̂ + T · t̂ + N · n̂                                          (0.3.1.1)
```

where `R, T, N` are the radial, transverse, and normal scalar components,
respectively.

**Remark 0.3.1.2 (Velocity in this frame).** Differentiating `r = r r̂` with
respect to time, and using `dr̂/dt = ḟ t̂` (the unit-vector identity for a
rotating planar frame), gives

```
v = ṙ r̂ + r ḟ t̂ .                                                    (0.3.1.2)
```

There is no `n̂` component in `v` because the orbit lies in the plane spanned by
`(r̂, t̂)`.

**Alignment to SGP4.** The `(r̂, t̂, n̂)` decomposition is frame-agnostic — it is
defined relative to the instantaneous orbital plane. In SGP4 the orbital plane
itself is parameterized by `(i, Ω)` in the TEME frame; perturbations are computed
in `(R, T, N)` and then projected back to TEME via `(i, Ω, ω)` when constructing
the position/velocity output (in `state_from_elements.h`). Phase 0 derivations
work entirely in the orbital-plane frame; the TEME projection is downstream of
Phase 0.

**Alignment to drag specialization (preview).** SGP4's drag acceleration is
the **Newtonian quadratic drag** of main-doc Postulate 1.2:

```
F_drag = −B* · ρ(r) · |v| · v ,
```

with force magnitude `∝ |v|²` (since `|F_drag| = B* ρ |v|·|v| = B* ρ |v|²`) and
direction anti-velocity (`−v̂`). The decomposition into the `(r̂, t̂, n̂)` frame
uses the velocity (0.3.1.2) `v = ṙ r̂ + r ḟ t̂`:

```
R_drag = −B* · ρ · |v| · ṙ ,
T_drag = −B* · ρ · |v| · r ḟ = −B* · ρ · |v| · h/r ,
N_drag = 0 ,
```

with the speed scalar `|v| = √(ṙ² + (r ḟ)²) = (n a / β) · √(1 + e² + 2 e cos f)`
(per Theorem 0.4.1). The radial component `R_drag ∝ ṙ` vanishes at perigee and
apogee and is `O(e)` smaller than `T_drag` on average — but it is *not* zero,
and the Lane drag derivation retains it. See §0.4 (drag specialization) for
the derivation that keeps both `R_drag` and `T_drag`, and which surfaces the
`(1 + e² + 2 e cos f)^{3/2}` factor in the `ȧ_drag` closed form (0.4.2.1) as
the structural fingerprint of quadratic drag.

### Theorem 0.3.2 (Gauss variational equation for `ȧ`)

**Hypotheses.**

- (H₁) The satellite follows the perturbed two-body equation
  `r̈ = -μ r/r³ + F` where `F` is the perturbing acceleration per unit mass.
- (H₂) The orbit is bound (E < 0) and instantaneously characterized by
  osculating Keplerian elements per the VOP framework (Theorem 0.2.2 alignment).

**Conclusion.**

```
ȧ = (2 / (n β)) · [e R sin f + T · p / r]                            (0.3.2.1)
```

Equivalently, using `(1 + e cos f) = p/r` from Theorem 0.2.2:

```
ȧ = (2 / (n β)) · [e R sin f + T (1 + e cos f)] .                    (0.3.2.2)
```

**Proof.**

**Step 1 (Specific orbital energy in terms of `a`).** For a bound Keplerian
orbit, the specific orbital energy is `E = (1/2) v² − μ/r = −μ/(2a)` (standard
two-body energy integral; cite [WIKI-OST]). Differentiating with respect to time
under perturbation:

```
dE/dt = (μ / (2 a²)) · ȧ .                                           (0.3.2.3)
```

**Step 2 (Rate of energy from work-by-perturbation).** The Keplerian central
force `-μ r/r³` is conservative, so it contributes nothing to `dE/dt`. Only the
perturbation `F` does work:

```
dE/dt = v · F .                                                       (0.3.2.4)
```

This follows from `dE/dt = v · v̇ + (μ/r³) · r · v̇₀` where `v̇₀ = -μr/r³` is the
Keplerian acceleration and `v̇ = v̇₀ + F`. The Keplerian part cancels by energy
conservation under the unperturbed flow; only `v · F` remains.

**Step 3 (Substitute the frame decomposition).** From (0.3.1.1) and (0.3.1.2):

```
v · F = (ṙ r̂ + r ḟ t̂) · (R r̂ + T t̂ + N n̂)
      = ṙ · R + r ḟ · T .                                             (0.3.2.5)
```

The `n̂` component vanishes because `v` has no `n̂` component.

**Step 4 (Express `ṙ` and `r ḟ` in elements).** From Theorem 0.2.2,
`r = p/(1 + e cos f)`. Differentiating:

```
ṙ = (∂r/∂f) · ḟ = (p · e sin f / (1 + e cos f)²) · ḟ .               (0.3.2.6)
```

From Theorem 0.2.3, `r² ḟ = h`, so `ḟ = h/r²`. Substituting into (0.3.2.6):

```
ṙ = (p · e sin f / (1 + e cos f)²) · (h / r²)
  = (p · e sin f · h) / ((1 + e cos f)² · r²) .                       (0.3.2.7)
```

By Theorem 0.2.2 `r = p/(1+e cos f)`, so `r² (1+e cos f)² = p²`. Substituting:

```
ṙ = (p · e sin f · h) / p² = (e h sin f) / p .                       (0.3.2.8)
```

For `r ḟ`, again from `r² ḟ = h`:

```
r ḟ = h / r = h (1 + e cos f) / p .                                  (0.3.2.9)
```

**Step 5 (Substitute Theorem 0.2.4 to eliminate `h` in favour of `n, a, β`).**
By (0.2.4.2), `h = n a² β`. Also `p = a β²`. So `h/p = n a² β / (a β²) = n a / β`.
Substituting:

```
ṙ = e · (n a / β) · sin f = (n a e sin f) / β ,                      (0.3.2.10)
r ḟ = (n a / β) · (1 + e cos f) .                                    (0.3.2.11)
```

**Step 6 (Assemble v · F).** Combining (0.3.2.5), (0.3.2.10), (0.3.2.11):

```
v · F = (n a e sin f / β) · R + ((n a / β) · (1 + e cos f)) · T
      = (n a / β) · [e R sin f + T (1 + e cos f)] .                  (0.3.2.12)
```

**Step 7 (Solve for `ȧ`).** Combining (0.3.2.3), (0.3.2.4), (0.3.2.12):

```
(μ / (2 a²)) · ȧ = (n a / β) · [e R sin f + T (1 + e cos f)]
```

```
ȧ = (2 a² · n a) / (μ β) · [e R sin f + T (1 + e cos f)]
  = (2 n a³ / (μ β)) · [e R sin f + T (1 + e cos f)] .               (0.3.2.13)
```

**Step 8 (Simplify the prefactor).** By Theorem 0.2.4 `n² a³ = μ`, so
`a³ = μ/n²`, and therefore `n a³ = μ/n`. Substituting in (0.3.2.13):

```
2 n a³ / (μ β) = 2 (μ/n) / (μ β) = 2 / (n β) .                       (0.3.2.14)
```

**Step 9 (Final form).** Substituting (0.3.2.14) into (0.3.2.13):

```
ȧ = (2 / (n β)) · [e R sin f + T (1 + e cos f)] .                    (0.3.2.15)
```

Using `(1 + e cos f) = p/r` (Theorem 0.2.2 rearranged), the equivalent compact
form is

```
ȧ = (2 / (n β)) · [e R sin f + T p / r] .                            (0.3.2.16)
```

(0.3.2.15) and (0.3.2.16) match the claimed forms (0.3.2.1) and (0.3.2.2). ∎

**Remark 0.3.2.17 (Sign convention and units).** `ȧ` is positive when the
perturbation does positive work on the satellite, growing the orbit; negative
when energy is removed. For drag (Step in §0.4), `T < 0` everywhere (drag opposes
motion), so `ȧ < 0` and the orbit decays — sign-consistent with physical
intuition. The units of (0.3.2.15) are `[L]/[T]`: `R, T` have units `[L]/[T]²`
(acceleration), `n` has units `[T]⁻¹`, so `R/(n β)` has units `[L]/[T]` ✓.

**Alignment to SGP4 (theorem-level).**
- (a) **Symbol bridge.** In Phase 0 generic notation `a` and `n`; in SGP4 these
  are `a₀''` and `n₀''` (Brouwer-recovered, per 0.2.4 alignment).
- (b) **VOP hypothesis (H₂).** `ȧ` here is the *instantaneous* rate; it varies
  along the orbit because `R(f)` and `T(f)` depend on `f` (through the density
  profile in drag, through the inclination in J₂, etc.). The *secular* rate
  `⟨ȧ⟩` is obtained by orbit-averaging (0.3.2.15) via Definition 0.2.6,
  context-1 — see Theorem 0.5.X (energy identity, Phase 0 §0.5) and Phase 2
  (C₂ assembly).
- (c) **No truncation has been applied.** (0.3.2.15) is exact under (H₁)-(H₂).
  Truncations enter only when (i) `R, T` themselves are approximated (e.g. the
  Lane density model truncated at `τ = 4`), or (ii) the orbit average is
  evaluated with the Lane `f†` substitution (a separate approximation, Phase 11).
  In particular: the AFGP4 → SGP4 dropping of `O(e²)` terms is *not* in the
  Gauss VE (0.3.2.15); it is in the subsequent orbit-averaging step.
- (d) **`Ω̇, i̇`-free.** The drag application of (0.3.2.15) inherits `N = 0` (§0.4),
  so out-of-plane couplings to `Ω̇, i̇` do not enter `ȧ`. They will enter the
  J₃-coupling Phase 3 / Phase 9 (where the radial out-of-phase perturbation
  `δr_{J_3} ∝ sin u` produces effective `N ≠ 0`).

**Alignment to implementation (deferred per Standard 9-B).** None at Phase 0:
(0.3.2.15) is the canonical Gauss form. Implementation choices (Horner-form
evaluation of the bracket, pre-computed `2 / (n β)` etc.) are downstream.

---

### Lemma 0.3.3.a (Eccentricity vector identity)

For a Keplerian orbit with position `r`, velocity `v`, and angular momentum
`h := r × v`, the **eccentricity vector**

```
e_vec := (v × h)/μ − r̂                                                (0.3.3.a.1)
```

has magnitude equal to the scalar eccentricity `e` and points from the focus
toward perigee.

**Proof.** Standard; cite [WIKI-OST §"Eccentricity vector"]. The two-line
verification: differentiate (0.3.3.a.1) under the unperturbed Kepler equation
`v̇ = −μ r̂/r²`. By Step 4 below applied with `F = 0`, the result is zero, so
`e_vec` is conserved. Its magnitude can be evaluated at perigee where `ṙ = 0`:
`e_vec|_perigee = (v_p × h)/μ − r̂_perigee`. With `|v_p| = √(μ(1+e)/(a(1−e)))` and
`v_p ⊥ r̂_perigee`, `|v_p × h| = |v_p| · h = √(μ(1+e)/(a(1−e))) · √(μa(1−e²))
= μ √((1+e)²/(1)) = μ(1+e)`. So `|e_vec|_perigee = (1+e) − 1 = e`. Direction:
`v_p × h` is along `r̂_perigee` (by right-hand rule with `v_p` transverse,
`h` orbit-normal), so `e_vec|_perigee = e r̂_perigee` ⇒ `e_vec` points toward
perigee. ∎

**Alignment to SGP4.** The eccentricity vector is the canonical bridge between
the Cartesian state `(r, v)` and the in-plane elements `(e, ω)`. In SGP4 the
"long-period periodic" corrections to the eccentricity vector (xlcof / aycof
in §14, derived in Phase 7) act directly on `e_vec`'s components in the orbital
plane, rather than on the scalar `(e, ω)` separately. The use of `e_vec`
throughout §0.3 makes that connection natural.

### Theorem 0.3.3 (Gauss variational equation for `ė`)

**Hypotheses.** As in Theorem 0.3.2 (perturbed two-body, VOP, osculating
elements).

**Conclusion.**

```
ė = (β / (n a)) · [R sin f + T (cos f + cos E)]                       (0.3.3.1)
```

where `E` is the eccentric anomaly and `cos E = (e + cos f)/(1 + e cos f)`.

Equivalently, eliminating `cos E`:

```
ė = (β / (n a)) · [R sin f + T · (2 cos f + e + e cos²f) / (1 + e cos f)] . (0.3.3.2)
```

**Proof.**

**Step 1 (Differentiate the eccentricity vector).** From Lemma 0.3.3.a,

```
de_vec/dt = (v̇ × h)/μ + (v × ḣ)/μ − dr̂/dt .                        (0.3.3.3)
```

**Step 2 (Apply perturbed two-body identities).** Under hypothesis (H₁),
`v̇ = −μ r̂/r² + F` and `ḣ = d(r × v)/dt = ṙ × v + r × v̇ = v × v + r × v̇
= r × (−μ r̂/r² + F) = r × F` (since `r × r̂ = 0`).

Substituting into (0.3.3.3):

```
de_vec/dt = ((−μ r̂/r² + F) × h)/μ + (v × (r × F))/μ − dr̂/dt
          = −((r̂ × h))/r² + (F × h)/μ + (v × (r × F))/μ − dr̂/dt . (0.3.3.4)
```

**Step 3 (Identify and cancel the Keplerian term).** With `h = h n̂` (orbit
normal), and the right-handed orthonormal triad `(r̂, t̂, n̂)`:
`r̂ × h = r̂ × (h n̂) = h (r̂ × n̂) = −h t̂` (right-hand rule).
So `−(r̂ × h)/r² = h t̂ / r² = ḟ t̂` (using `r² ḟ = h` from Theorem 0.2.3).

Also `dr̂/dt = ḟ t̂` (derivative of a unit vector in a rotating planar frame).
Therefore the Keplerian and `dr̂/dt` terms cancel:

```
−(r̂ × h)/r² − dr̂/dt = ḟ t̂ − ḟ t̂ = 0 .                              (0.3.3.5)
```

This confirms `e_vec` is conserved under unperturbed flow (matches Lemma 0.3.3.a).
The perturbation-driven part of (0.3.3.4) is therefore

```
de_vec/dt = (F × h)/μ + (v × (r × F))/μ .                           (0.3.3.6)
```

**Step 4 (Evaluate `F × h`).** With `h = h n̂` and `F = R r̂ + T t̂ + N n̂`:

```
F × h = (R r̂ + T t̂ + N n̂) × h n̂
      = h [R (r̂ × n̂) + T (t̂ × n̂) + N (n̂ × n̂)]
      = h [R (−t̂) + T (r̂) + 0]
      = h (T r̂ − R t̂) .                                              (0.3.3.7)
```

**Step 5 (Evaluate `v × (r × F)` using the BAC-CAB identity).**

```
v × (r × F) = r (v · F) − F (v · r) .                                (0.3.3.8)
```

From (0.3.2.5) `v · F = R ṙ + T r ḟ`. From `v = ṙ r̂ + r ḟ t̂` and `r = r r̂`,
`v · r = ṙ r`. So

```
v × (r × F) = r r̂ · (R ṙ + T r ḟ) − F · (ṙ r)
            = r (R ṙ + T r ḟ) r̂ − ṙ r (R r̂ + T t̂ + N n̂)
            = r [(R ṙ + T r ḟ) − ṙ R] r̂ − ṙ r T t̂ − ṙ r N n̂
            = r (T r ḟ) r̂ − ṙ r T t̂ − ṙ r N n̂
            = (r² ḟ) T r̂ − r ṙ T t̂ − r ṙ N n̂
            = h T r̂ − r ṙ T t̂ − r ṙ N n̂                              (0.3.3.9)
```

(using `r² ḟ = h` in the last step).

**Step 6 (Sum and collect components).** Combining (0.3.3.7) and (0.3.3.9) into
(0.3.3.6):

```
de_vec/dt = (1/μ) [h (T r̂ − R t̂) + h T r̂ − r ṙ T t̂ − r ṙ N n̂]
          = (1/μ) [2 h T · r̂ − (h R + r ṙ T) · t̂ − r ṙ N · n̂] .   (0.3.3.10)
```

**Step 7 (Project onto the perigee direction `P̂`).** Let `P̂` be the unit vector
in the orbital plane pointing from focus toward perigee (the direction of
`e_vec` itself), and `Q̂ := n̂ × P̂` be the in-plane orthogonal direction. Then
`r̂ = cos f · P̂ + sin f · Q̂` and `t̂ = −sin f · P̂ + cos f · Q̂` (standard
orbital-frame rotation by angle `f` from perigee).

The scalar eccentricity is `e = e_vec · P̂`, so `ė = de_vec/dt · P̂` provided
`P̂` itself is not rotating in the orbital plane. For purely in-plane
perturbations (`N = 0`, drag case), `P̂` is fixed in the orbital plane and this
holds. For general perturbations (`N ≠ 0`), `P̂` rotates at rate `ω̇` and an
additional term `−e ω̇ Q̂` appears in `de_vec/dt`; however, projecting onto `P̂`
extracts only the `ė` component since `Q̂ · P̂ = 0`. So in all cases

```
ė = de_vec/dt · P̂ .                                                  (0.3.3.11)
```

Project (0.3.3.10):

```
de_vec/dt · P̂ = (1/μ) [2 h T (r̂ · P̂) − (h R + r ṙ T)(t̂ · P̂) − r ṙ N (n̂ · P̂)]
              = (1/μ) [2 h T cos f − (h R + r ṙ T)(−sin f) − 0]
              = (1/μ) [2 h T cos f + (h R + r ṙ T) sin f]
              = (1/μ) [h R sin f + T (2 h cos f + r ṙ sin f)] .       (0.3.3.12)
```

**Step 8 (Simplify the `T` bracket using `cos E` identity).** Substitute
`r ṙ` from §0.3.2 Step 4. From (0.3.2.8) `ṙ = e h sin f / p`, and from
Theorem 0.2.2 `r = p/(1 + e cos f)`. Therefore

```
r ṙ = (p/(1 + e cos f)) · (e h sin f / p) = e h sin f / (1 + e cos f) . (0.3.3.13)
```

Substituting into the `T` bracket:

```
2 h cos f + r ṙ sin f = 2 h cos f + (e h sin² f)/(1 + e cos f)
                     = (h / (1 + e cos f)) · [2 cos f (1 + e cos f) + e sin² f]
                     = (h / (1 + e cos f)) · [2 cos f + 2 e cos² f + e (1 − cos² f)]
                     = (h / (1 + e cos f)) · [2 cos f + e cos² f + e]
                     = (h / (1 + e cos f)) · [2 cos f + e (1 + cos² f)] . (0.3.3.14)
```

Now invoke the **eccentric-anomaly identity** `cos E = (e + cos f)/(1 + e cos f)`
([WIKI-KE §"Relationship between eccentric and true anomaly"]). Adding `cos f`:

```
cos f + cos E = cos f + (e + cos f)/(1 + e cos f)
              = [cos f (1 + e cos f) + e + cos f]/(1 + e cos f)
              = [cos f + e cos² f + e + cos f]/(1 + e cos f)
              = [2 cos f + e (1 + cos² f)]/(1 + e cos f) .             (0.3.3.15)
```

Comparing (0.3.3.14) and (0.3.3.15), the bracket equals `h (cos f + cos E)`:

```
2 h cos f + r ṙ sin f = h (cos f + cos E) .                          (0.3.3.16)
```

**Step 9 (Final assembly).** Substitute (0.3.3.16) into (0.3.3.12):

```
ė = (1/μ) [h R sin f + T · h (cos f + cos E)]
   = (h/μ) [R sin f + T (cos f + cos E)] .                            (0.3.3.17)
```

By Theorem 0.2.4 (0.2.4.2), `h = n a² β`, and `μ = n² a³` (Kepler 3rd), so
`h/μ = (n a² β)/(n² a³) = β / (n a)`. Substituting:

```
ė = (β / (n a)) · [R sin f + T (cos f + cos E)] .                    (0.3.3.18)
```

This matches (0.3.3.1). The equivalent form (0.3.3.2) follows by substituting
`cos E = (e + cos f)/(1 + e cos f)` from (0.3.3.15) and simplifying:
`cos f + cos E = (2 cos f + e + e cos² f)/(1 + e cos f)`. ∎

**Remark 0.3.3.19 (Where `N` went).** The eccentricity vector has no `n̂`
component for in-plane perturbations, and the projection onto `P̂` zeros out the
`n̂` component for any perturbation. So `N` does not appear in `ė` at all —
out-of-plane perturbations do not change the eccentricity at first order. (They
do change `i` and `Ω`, derived in Theorem 0.3.5.)

**Alignment to SGP4.**
- (a) **Symbol bridge.** `n → n₀''`, `a → a₀''`, `e → e₀` plus long-period
  corrections from §14 / Phase 7. `E` is solved iteratively from Kepler's
  equation `M = E − e sin E` (the modified-Kepler iterator in
  `src/orbit/modified_kepler.h`, which is Halley-method per audit R03).
- (b) **VOP hypothesis.** As in Theorem 0.3.2, `ė` here is the instantaneous
  rate. The secular `⟨ė⟩` is obtained by orbit-averaging — this is the
  derivation that produces `C₄` (eccentricity decay) in Phase 4.
- (c) **No truncation.** (0.3.3.1) is exact under (H₁)-(H₂).
- (d) **`N`-independence.** `ė` depends only on `R` and `T`, not on `N`. The
  SGP4 J₂ short-period out-of-plane perturbation contributes via `R` (radial)
  not `N` (out-of-plane), so its effect on `e` is captured here. The J₃
  long-period out-of-plane perturbation has `δr_{J_3} ∝ sin u · sin i`
  primarily, which appears as `R` through the density-coupling chain in Phase 3.

**Alignment to implementation (deferred per Standard 9-B).** None at Phase 0.
Implementation choices in `C₃ / C₄` use Horner-form polynomial evaluation and
exact-rational coefficient construction (`ratio<T>(p,q)`); these are downstream.

---

