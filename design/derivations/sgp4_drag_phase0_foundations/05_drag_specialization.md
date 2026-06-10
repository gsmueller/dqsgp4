## §0.4 Drag Specialization

The Gauss VE in §0.3 are general — they hold for any perturbing acceleration
`F = R r̂ + T t̂ + N n̂`. This section specializes to SGP4's drag model:

```
F_drag := −B*(t) · ρ(r) · |v| · v ,                                   (0.4.1)
```

where `B*` is the modified ballistic coefficient (Definition 1.1 of the main
document — flagged with **D-1** for the missing `ρ_0` factor in its dimensional
form; the corrected definition is `B* = C_D · ρ_0 · A / (2 m)` with units
ER⁻¹), `ρ(r)` is the atmospheric density at geocentric distance `r` (Lane
power-law model, postulated separately in §2 of the main document), `v` is
the (non-rotating) satellite velocity vector, and `|v| := ‖v‖` is the speed
scalar.

The `|v| · v` form (main-doc Postulate 1.2's `|v_rel| · v_rel` under the Lane
non-rotating-atmosphere assumption) is the **Newtonian quadratic drag law**:
the force magnitude scales as `|v|²` (since `|F_drag| = B* · ρ · |v| · |v| =
B* · ρ · |v|²`) and the force direction is `−v̂` (anti-velocity). The
factorization `|v| · v` (rather than `v² · v̂`) is operationally cleaner for
the Gauss-VE decomposition because `v` already decomposes into `(ṙ, r ḟ)`
components per (0.3.1.2).

### Theorem 0.4.1 (Decomposition of `F_drag` in the (r̂, t̂, n̂) frame)

**Hypotheses.**
- (H₁) Perturbed two-body equation `r̈ = -μ r/r³ + F_drag` with `F_drag` from
  (0.4.1).
- (H₂) Lane non-rotating-atmosphere assumption: `v_rel = v` (main doc
  Assumption 1.3).

**Conclusion.**

```
R_drag = −B* · ρ · |v| · ṙ ,   T_drag = −B* · ρ · |v| · r ḟ ,   N_drag = 0 . (0.4.1.1)
```

The speed scalar `|v|` is

```
|v| = √(ṙ² + (r ḟ)²) ,                                                  (0.4.1.2)
```

with `ṙ` and `r ḟ` from (0.3.2.10) and (0.3.2.11). Substituting:

```
|v| = (n a / β) · √(1 + e² + 2 e cos f) .                              (0.4.1.3)
```

Equivalently in terms of orbit elements:

```
R_drag = −B* · ρ · (n² a² e sin f / β²) · √(1 + e² + 2 e cos f) ,
T_drag = −B* · ρ · (n² a² (1 + e cos f) / β²) · √(1 + e² + 2 e cos f) . (0.4.1.4)
```

**Proof.**

**Step 1 (Velocity decomposition).** From (0.3.1.2), `v = ṙ r̂ + r ḟ t̂` —
the velocity has only `r̂` and `t̂` components, no `n̂` component (since the
orbit is planar). The speed scalar is

```
|v|² = v · v = (ṙ r̂ + r ḟ t̂) · (ṙ r̂ + r ḟ t̂) = ṙ² + (r ḟ)²        (0.4.1.5)
```

(cross-term `2 ṙ · (r ḟ) (r̂ · t̂)` vanishes since `r̂ ⊥ t̂`). Taking the
positive square root (the speed is non-negative by definition) gives (0.4.1.2).

**Step 2 (Apply the drag postulate with the `|v|` factor).** From (0.4.1):

```
F_drag = (−B* · ρ · |v|) · v
       = (−B* · ρ · |v|) · (ṙ r̂ + r ḟ t̂)
       = (−B* · ρ · |v| · ṙ) r̂ + (−B* · ρ · |v| · r ḟ) t̂ + 0 · n̂ . (0.4.1.6)
```

Comparing to `F = R r̂ + T t̂ + N n̂`:

```
R_drag = −B* · ρ · |v| · ṙ ,   T_drag = −B* · ρ · |v| · r ḟ ,   N_drag = 0 . (0.4.1.7)
```

This matches the conclusion (0.4.1.1).

**Step 3 (Substitute orbit-element forms of `ṙ` and `r ḟ`).** From (0.3.2.10)
`ṙ = (n a e sin f) / β` and from (0.3.2.11) `r ḟ = (n a / β)(1 + e cos f)`.
Substituting into (0.4.1.7):

```
R_drag = −B* · ρ · |v| · (n a e sin f / β) ,
T_drag = −B* · ρ · |v| · (n a / β) (1 + e cos f) .                    (0.4.1.8)
```

**Step 4 (Closed form for `|v|`).** Combining Step 1's `|v|² = ṙ² + (r ḟ)²`
with the orbit-element forms of Step 3:

```
|v|² = (n a e sin f / β)² + (n a (1 + e cos f) / β)²
     = (n a / β)² · [e² sin² f + (1 + e cos f)²] .                     (0.4.1.9)
```

Expand the bracket:

```
e² sin² f + (1 + e cos f)² = e² sin² f + 1 + 2 e cos f + e² cos² f
                          = 1 + 2 e cos f + e² (sin² f + cos² f)
                          = 1 + e² + 2 e cos f .                       (0.4.1.10)
```

So `|v|² = (n a / β)² · (1 + e² + 2 e cos f)`, and the positive square root
gives

```
|v| = (n a / β) · √(1 + e² + 2 e cos f) ,                              (0.4.1.11)
```

matching (0.4.1.3). Substituting (0.4.1.11) into (0.4.1.8) yields (0.4.1.4). ∎

**Remark 0.4.1.6 (Relative magnitudes of R_drag and T_drag).** The ratio
`R_drag/T_drag = ṙ/(r ḟ) = (e sin f)/(1 + e cos f)`. For typical LEO orbits
`e ∼ 0.01 − 0.1`, so `|R_drag/T_drag| ≲ e ≪ 1`. The radial drag is `O(e)`
smaller than the transverse on average, but is non-zero — and crucially the
SGP4 derivation retains both. The C₂ Part-A polynomial `(1 + (3/2)η² + 4e η + e η³)`
that emerges in Phase 2 contains contributions from both R and T components;
dropping `R_drag` would be a different (and wrong) approximation.

**Alignment to SGP4 (critical).**
- (a) **The non-rotating-atmosphere assumption.** Theorem 0.4.1 uses the SGP4
  postulate `v_rel = v` (Assumption 1.3 of the main document). For Earth's actual
  atmosphere `v_rel = v − ω_⊕ × r` with `|ω_⊕ × r| ∼ 0.46 km/s` at LEO; against
  satellite speed `∼ 7.7 km/s` this is `∼ 6 %`. The drag force `F ∝ v_rel² ≈
  v² (1 ∓ 2 ω_⊕ r cos i / v)`, so the relative-velocity assumption introduces a
  systematic accuracy error of order `10⁻¹` (with inclination-dependent sign).
  This is catalogued as `A-D1` in §17 of the main document. Acceptance of this
  error is a SGP4 model choice, not a derivation error.
- (b) **The Lane density model.** `ρ(r)` is postulated as the power-law form
  (Lane 1965). Its derivation lives in §2 of the main document plus Phase 11
  (Lane f† formalization). For Phase 0 purposes, `ρ(r)` is a given function of
  `r`, and only its evaluation at the orbit position enters (0.4.1.5).
- (c) **Specialization of §0.3 Gauss VE.** Substituting `N_drag = 0` into the
  alignment remarks for Theorems 0.3.2-0.3.6:
  - `ȧ_drag` from (0.3.2.15) with full `R, T`.
  - `ė_drag` from (0.3.3.1) with full `R, T`.
  - `Ω̇_drag = 0`, `i̇_drag = 0` from Theorem 0.3.4.
  - `ω̇_drag` from (0.3.5.2) with no Ω̇ correction.
  - `Ṁ_drag` from (0.3.6.1) with full `R, T`.

### Theorem 0.4.2 (Closed-form drag rates of orbit elements)

**Hypotheses.** As in Theorem 0.4.1: corrected quadratic-drag form
`F_drag = −B*·ρ·|v|·v` with `R_drag, T_drag, N_drag` from (0.4.1.7), and the
speed-scalar closed form for `|v|` from (0.4.1.11).

**Conclusion.** Substituting (0.4.1.7) into the Gauss VE from §0.3 yields the
closed-form instantaneous drag rates:

```
ȧ_drag = −(2 n B* ρ · a² / β³) · (1 + e² + 2 e cos f)^{3/2} .          (0.4.2.1)
```

```
ė_drag = −(B* ρ |v|) · [e sin² f + (1 + e cos f) (cos f + cos E)] .    (0.4.2.2)
```

```
ω̇_drag = −(B* ρ |v| / e) · [−e sin f cos f + (1 + e cos f) sin f (1 + r/p)] . (0.4.2.3)
```

```
Ω̇_drag = 0 ,    i̇_drag = 0 .                                          (0.4.2.4)
```

```
Ṁ_drag = n − (2 r R_drag) / (n a²) + (β² / (n a e)) · [R_drag cos f − T_drag sin f (1 + r/p)] ,
                                                                       (0.4.2.5)
```

where `R_drag, T_drag` are substituted from (0.4.1.7) (each carrying the `|v|`
factor) and `|v|` from (0.4.1.11) `(n a / β) · √(1 + e² + 2 e cos f)`.

**Proof.**

**Step 1 (`ȧ_drag`, full simplification).** Substitute `R_drag, T_drag` from
(0.4.1.7) into Gauss VE Theorem 0.3.2 (0.3.2.15):

```
ȧ_drag = (2 / (n β)) · [e R_drag sin f + T_drag (1 + e cos f)]
        = (2 / (n β)) · [e · (−B* ρ |v| · n a e sin f / β) · sin f
                        + (−B* ρ |v| · n a (1 + e cos f) / β) · (1 + e cos f)]
        = (2 / (n β)) · (−B* ρ |v| · n a / β) · [e² sin² f + (1 + e cos f)²]
        = −(2 B* ρ |v| · a / β²) · [e² sin² f + (1 + e cos f)²] .      (0.4.2.6)
```

By (0.4.1.10), `e² sin² f + (1 + e cos f)² = 1 + e² + 2 e cos f`. Substituting:

```
ȧ_drag = −(2 B* ρ |v| · a / β²) · (1 + e² + 2 e cos f) .              (0.4.2.7)
```

Substitute the closed form for `|v|` from (0.4.1.11):

```
ȧ_drag = −(2 B* ρ · a / β²) · (n a / β) · √(1 + e² + 2 e cos f) · (1 + e² + 2 e cos f)
        = −(2 n B* ρ · a² / β³) · (1 + e² + 2 e cos f)^{3/2} .         (0.4.2.8)
```

This is (0.4.2.1). The `^{3/2}` exponent emerges from the combination of
`|v|`'s `^{1/2}` (the speed scalar under quadratic drag, by (0.4.1.11)) with
the bracket factor's `^1`. This is the structural fingerprint of Newtonian
quadratic drag in the orbit-averaged-rate formulation, and is exactly what was
**missing** before D-10 was fixed.

**Step 2 (`ė_drag`).** Substitute `R_drag, T_drag` from (0.4.1.7) into Gauss
VE Theorem 0.3.3 (0.3.3.18):

```
ė_drag = (β / (n a)) · [R_drag sin f + T_drag (cos f + cos E)]
        = (β / (n a)) · [(−B* ρ |v| · n a e sin f / β) · sin f
                        + (−B* ρ |v| · n a (1 + e cos f) / β) · (cos f + cos E)]
        = (β / (n a)) · (−B* ρ |v| · n a / β) · [e sin² f + (1 + e cos f)(cos f + cos E)]
        = −(B* ρ |v|) · [e sin² f + (1 + e cos f)(cos f + cos E)] .   (0.4.2.9)
```

This is (0.4.2.2). Unlike Step 1, the trigonometric bracket here does **not**
factor into `(1 + e² + 2 e cos f)`, so `|v|` remains as an explicit
multiplicative scalar. The orbit-averaging step in Phase 4 will expand
`|v| · [bracket]` and project onto the Lane density profile.

**Step 3 (`ω̇_drag`).** Substitute `R_drag, T_drag` from (0.4.1.7) into Gauss
VE Theorem 0.3.5 (0.3.5.2) (with `Ω̇ = 0` since `N_drag = 0`):

```
ω̇_drag = (β / (n a e)) · [−R_drag cos f + T_drag sin f (1 + r/p)]
         = (β / (n a e)) · [B* ρ |v| · (n a e sin f / β) · cos f
                           − B* ρ |v| · (n a / β)(1 + e cos f) · sin f · (1 + r/p)]
         = (β / (n a e)) · (B* ρ |v| · n a / β) · [e sin f cos f − (1 + e cos f) sin f (1 + r/p)]
         = (B* ρ |v| / e) · sin f · [e cos f − (1 + e cos f)(1 + r/p)]
         = −(B* ρ |v| / e) · [−e sin f cos f + (1 + e cos f) sin f (1 + r/p)] . (0.4.2.10)
```

This is (0.4.2.3). Same `|v|`-as-multiplicative-scalar observation as Step 2.

**Step 4 (`Ω̇_drag, i̇_drag`).** By Theorem 0.3.4, both `Ω̇` and `i̇` are
proportional to `N` (the out-of-plane perturbation component). From (0.4.1.7),
`N_drag = 0`, so

```
Ω̇_drag = (r · 0 · sin u) / (n a² β sin i) = 0 ,
i̇_drag = (r · 0 · cos u) / (n a² β) = 0 .                              (0.4.2.11)
```

This is (0.4.2.4). The vanishing is **independent of `|v|`** — the `|v|`
factor only enters through `R, T`, not through `N`. The drag does not rotate
the orbital plane.

**Step 5 (`Ṁ_drag`).** Substitute `R_drag, T_drag` from (0.4.1.7) into Gauss
VE Theorem 0.3.6 (0.3.6.22):

```
Ṁ_drag = n − (2 r R_drag) / (n a²) + (β² / (n a e)) · [R_drag cos f − T_drag sin f (1 + r/p)] .
                                                                       (0.4.2.12)
```

Each of `R_drag` and `T_drag` here carries the `|v|` factor per (0.4.1.7);
substituting `|v|` from (0.4.1.11) produces an expanded form involving
`√(1 + e² + 2 e cos f)` as a multiplicative factor in every drag-driven term.
For the Phase 5 / Phase 10 orbit-averaging step, the Gauss-VE structure of
(0.4.2.12) is the more useful representation (the orbit average projects the
trigonometric bracket onto specific harmonics), so we leave (0.4.2.5) in this
form rather than expanding further. This is (0.4.2.5). ∎

**Remark 0.4.2.6.** (0.4.2.1)-(0.4.2.5) are the **instantaneous** drag rates,
varying around an orbit as `(f, ρ(r(f)))` vary. They are not yet the **secular
rates** `⟨ȧ⟩`, `⟨ė⟩`, etc. — those follow by orbit-averaging via Definition 0.2.6
context-1 (the C₁..C₅ coefficient definitions of Phase 2-5).

**Alignment to SGP4 (final, Phase 0-rev1).**
- (a) The negative sign in (0.4.2.1)-(0.4.2.5) reflects energy loss: drag
  removes orbital energy, so `ȧ_drag < 0` and the orbit decays. ✓
- (b) `Ω̇_drag = 0` and `i̇_drag = 0` from Theorem 0.3.4 (and Step 4 of the
  proof above) confirm the SGP4 observed behaviour that drag does not rotate
  the orbital plane — only in-plane elements decay. This conclusion is
  **independent of `|v|`** (the `|v|` factor enters only `R, T`, never `N`).
- (c) The `(1 + e² + 2 e cos f)^{3/2}` factor in (0.4.2.1) and the `|v| · [...]`
  multiplicative structure in (0.4.2.2)-(0.4.2.3) are the leading-order
  trigonometric content of the C-coefficient integrands. The `^{3/2}` exponent
  in particular is the structural fingerprint of **Newtonian quadratic drag**
  (Postulate 1.2): it emerges from `|v|^{1/2}` (the speed scalar's structure
  via (0.4.1.11)) combined with the bracket factor's `^1`. This factor was
  **missing** prior to Phase 0-rev1 (audit finding **D-10**) and is restored
  here. Phase 2.A's symbolic trace will project this against the Lane density
  profile `ρ(r(f))` via Phase 1's Lane integrals, yielding the SGP4 C₁..C₅
  closed forms.
- (d) **D-10 closure** (Phase 0-rev1). The corrected Theorem 0.4.1 + 0.4.2
  + Theorem 0.5.3 close audit finding D-10. The Postulate-1.2 match check
  added under Block A4 to `verify_phase0.m` is the verification that would
  have caught D-10 originally and now guards against its recurrence.

---

