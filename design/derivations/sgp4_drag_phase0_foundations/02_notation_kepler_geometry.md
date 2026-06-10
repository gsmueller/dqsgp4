## §0.1 Notation

| Symbol | Definition | Units |
|---|---|---|
| `μ` | Gravitational parameter `GM_⊕` of the central body | m³·s⁻² (ER³·min⁻² in SGP4 units) |
| `a` | Semi-major axis of the orbit | m (ER) |
| `e` | Eccentricity (`0 ≤ e < 1` for bound orbits) | dimensionless |
| `i` | Inclination of the orbit plane to the equator | rad |
| `Ω` | Right ascension of the ascending node (RAAN) | rad |
| `ω` | Argument of perigee | rad |
| `f` | True anomaly (angle from perigee to satellite, measured from focus) | rad |
| `E` | Eccentric anomaly | rad |
| `M` | Mean anomaly, `M := E − e sin E` (Kepler's equation) | rad |
| `n` | Mean motion, `n := √(μ/a³)` | rad·s⁻¹ |
| `β` | `β := √(1−e²)` | dimensionless |
| `p` | Semi-latus rectum, `p := a(1−e²) = a β²` | m (ER) |
| `r` | Geocentric distance, `r := p/(1 + e cos f)` (function of f) | m (ER) |
| `h` | Specific angular momentum, `h := \|r × v\|` | m²·s⁻¹ |
| `θ` | `θ := cos i` | dimensionless |
| `u` | Argument of latitude, `u := ω + f` | rad |
| `r̂, t̂, n̂` | Radial / transverse / orbit-normal unit vectors | dimensionless |
| `R, T, N` | Components of a perturbing acceleration in the (`r̂, t̂, n̂`) frame | m·s⁻² |
| `⟨X⟩_M` | Orbit average of `X(f)` over mean anomaly: `(1/2π) ∫₀^{2π} X dM` | matches X |

**Conventions.** Vectors are bolded only in physical context; in algebraic
manipulation we treat scalar component equations.  The orbital frame `(r̂, t̂, n̂)` is
defined so that `t̂` lies in the orbital plane perpendicular to `r̂` in the direction
of motion (true-anomaly-increasing direction), and `n̂ := r̂ × t̂` is normal to the
orbit.  Times are SI seconds in §0.1–§0.4; the SGP4 minute-based convention is
re-established in §0.5.

---

## §0.2 Kepler Geometry

This section establishes the two-body Keplerian primitives needed downstream: the
conic-section radius `r(f)`, angular momentum conservation `r²ḟ = h`, the third-law
relation `n²a³ = μ`, the differential relation `dM/df`, and the orbit-average
operator.

### Definition 0.2.1 (Two-body Keplerian orbit)

Let `r ∈ ℝ³` denote the position of a particle relative to a central body of
gravitational parameter `μ > 0`, and let `v := ṙ`. The **two-body Keplerian orbit**
is the trajectory satisfying

```
r̈ = −μ r / ‖r‖³ .                                                    (0.2.1.1)
```

By a standard reduction (energy + angular-momentum integrals; see [WIKI-OE], [BATT99
§3]), any bounded solution of (0.2.1.1) lies in a fixed plane, and within that
plane the trajectory is a conic section (ellipse for `e < 1`) with one focus at the
origin.

### Theorem 0.2.2 (Conic section in polar form)

**Hypothesis.** Let `r(f)` be the radius along a bound Keplerian orbit
(Definition 0.2.1), parameterized by the true anomaly `f` measured from perigee.

**Conclusion.**

```
r(f) = p / (1 + e cos f)   where   p := a(1 − e²) = a β² .            (0.2.2.1)
```

**Proof.** Standard conic-section derivation; cite [WIKI-OE] for the elementary
two-step proof (angular-momentum conservation + Binet's equation reduces (0.2.1.1)
to the linear ODE `u″ + u = μ/h²` for `u := 1/r`, whose general solution
`u = (μ/h²)(1 + e cos(f − f₀))` with `f₀ := 0` at perigee yields (0.2.2.1) after
identifying `p = h²/μ` and using `h² = μ a (1−e²) = μ p`). ∎

**Remark 0.2.2.2 (Range of r).** From (0.2.2.1), the perigee radius is
`r_p := p/(1+e) = a(1−e)` (at `f = 0`) and the apogee radius is
`r_a := p/(1−e) = a(1+e)` (at `f = π`).

**Alignment to SGP4.** Theorem 0.2.2 is stated for an *unperturbed* Keplerian orbit
— `(a, e)` are constants of motion of (0.2.1.1). SGP4 propagates a *perturbed*
orbit where `(a, e)` are themselves time-dependent under J₂-secular drift, drag
decay, and J₃ long-period periodics. The relation `r = p/(1 + e cos f)` is applied
in SGP4 under the **variation-of-parameters (VOP) / osculating-elements**
formulation: at each instant `t`, the satellite is treated as an instantaneously
Keplerian orbit with elements `(a(t), e(t), …)` taken from the perturbation
integrators, and `r(f)` is evaluated with those instantaneous elements.

The VOP applicability requires that the elements vary on a time scale much longer
than one orbital period (the *secularity* hypothesis). For SGP4 near-earth LEO,
drag rates are `|ȧ/a| ∼ 10⁻⁷ min⁻¹` against an orbital frequency `n ∼ 10⁻¹ min⁻¹`;
the timescale separation is ~6 orders of magnitude. The hypothesis holds with
wide margin.

In SGP4 symbol bridge: `a → a₀''` (Brouwer-recovered, evaluated at epoch
+ secular drift), `e → e₀` plus the C₃-driven long-period correction in §7 / §15.
`f` is solved via the modified Kepler equation each propagation step.

### Theorem 0.2.3 (Angular momentum constant; Kepler's second law)

**Hypothesis.** As in Theorem 0.2.2.

**Conclusion.**

```
r² ḟ = h     where     h = √(μ p) = √(μ a β²) = √(μ a) · β .          (0.2.3.1)
```

**Proof.**

**Step 1 (Cross product is conserved).** From (0.2.1.1), the time derivative
`d(r × v)/dt = ṙ × v + r × v̇ = v × v + r × (−μ r/r³) = 0`. Hence `r × v` is
constant in time. Call its magnitude `h := ‖r × v‖`.

**Step 2 (Component identification).** Decompose `v = ṙ r̂ + r ḟ t̂` (radial plus
transverse). Then `r × v = r r̂ × (ṙ r̂ + r ḟ t̂) = r² ḟ (r̂ × t̂) = r² ḟ n̂`. Taking
magnitudes, `r² ḟ = h`.

**Step 3 (Identify h in terms of orbit elements).** At perigee `f = 0`, ṙ = 0 and
`r = r_p = a(1−e)`. From energy conservation
`(1/2)v² − μ/r = −μ/(2a)`, the speed at perigee is `v_p² = μ(1+e)/(a(1−e))`. So
`h = r_p · v_p = a(1−e) · √(μ(1+e)/(a(1−e))) = √(μ a (1+e)(1−e)) = √(μ a (1−e²))
= √(μ a β²) = √(μ a) · β = √(μ p)`. ∎

**Alignment to SGP4.** Theorem 0.2.3's `h = √(μ p)` is the angular-momentum
magnitude of the *unperturbed* Keplerian orbit with semi-latus rectum `p = a β²`.
Under SGP4's VOP (Theorem 0.2.2 alignment remark), `h` is the angular momentum of
the **Brouwer-recovered reference Keplerian** evaluated at the current osculating
elements: `h ≡ h(t) = √(μ · a₀''(t) · β₀(t)²)`. Drag along-track forces act on
`a` (and through it on `h`); the rate `dh/dt|_{drag}` enters §0.3 Theorem 0.3.3
(`ė` Gauss VE) and §0.4 drag specialization through the link `dh = (1/(2h))·d(μp)
= (μ/(2h))·(β²·da + 2 a β·dβ)`. The relation `r²ḟ = h` is invoked instantaneously
with the current `h(t)`.

The conserved-`h` interpretation from Step 1 is valid only for the unperturbed
problem; under drag, `h(t)` decays monotonically. The drag derivations consume
Theorem 0.2.3 in two distinct modes: (i) as an algebraic identity at each instant
of time (always valid), and (ii) as a Jacobian for orbit averaging at fixed
elements (valid under the secularity hypothesis, Theorem 0.2.5 alignment).

### Theorem 0.2.4 (Kepler's third law)

**Hypothesis.** As in Theorem 0.2.2.

**Conclusion.** The mean motion `n := 2π/T` (where T is the orbital period)
satisfies

```
n² a³ = μ ,           equivalently         n a² β = h .               (0.2.4.1)
```

**Proof.**

**Step 1 (Area-rate from Step 2 of 0.2.3).** The areal velocity is
`dA/dt = (1/2) r² ḟ = h/2` (constant). Integrating over one orbital period,
`A = h T / 2`.

**Step 2 (Total area).** For the ellipse, `A = π a b = π a · a β = π a² β`.

**Step 3 (Combine).** `π a² β = h T / 2` ⇒ `T = 2π a² β / h`. Using
`h = √(μ a) · β` from Theorem 0.2.3, `T = 2π a² β / (√(μ a) · β) = 2π a² / √(μ a)
= 2π √(a³/μ)`. Hence `n = 2π/T = √(μ/a³)` ⇒ `n² a³ = μ`.

**Step 4 (Alternative form).** Combine `n = √(μ/a³)` from Step 3 with `h = √(μ a) · β`
from Theorem 0.2.3:

```
n · a² · β = √(μ/a³) · a² · β = √((μ/a³) · a⁴) · β = √(μ a) · β = h . (0.2.4.2)
```

Equation (0.2.4.2) is the form used downstream: it expresses `h` as `n a² β`, which
appears in every Gauss-Lagrange variational equation in §0.3. ∎

**Alignment to SGP4 (critical).** "Mean motion `n`" in Theorem 0.2.4 is the
**instantaneous Keplerian** mean motion satisfying `n² a³ = μ` at the current value
of `a`. SGP4 stores and propagates the **Brouwer-recovered** mean motion `n₀''`,
which is NOT the instantaneous Keplerian mean motion — it is the value that, when
combined with the Brouwer J₂-secular rate generator, reproduces the time-averaged
true mean motion. The relation between the two is

```
n₀ = n₀''  ·  (1 + δ_n(a₀, e₀, i₀, J₂))                              (0.2.4.SGP4)
```

with `δ_n` the Brouwer recovery factor (derived in element_recovery.h and the
project's `deprecated/013` content). The TLE-stored `n_kozai` is `n₀` (Kozai-mean);
SGP4 inverts this to obtain `n₀''` and the corresponding `a₀''`.

In every downstream Phase, when a Phase-0 result uses the symbol `n`, it is
understood to be `n₀''` (and `a` is `a₀''`). The numerical recovery is performed
once at SGP4 initialization; the drag derivations consume `n₀'' / a₀''` as if they
were instantaneous-Keplerian elements. This is the "Brouwer secular elements
treated as the reference Keplerian" convention.

A consequence is that `h = n₀'' · a₀''² · β₀` from (0.2.4.2) is the angular
momentum *of the Brouwer-recovered reference Keplerian* — not the instantaneous
angular momentum of the true perturbed orbit, which contains additional short-period
J₂ wobble (the SPP corrections in `short_period.h`). The drag derivations operate
on the Brouwer reference and accept this offset as part of the SGP4 model.

### Theorem 0.2.5 (Differential relation `dM/df`)

**Hypothesis.** As in Theorem 0.2.2.

**Conclusion.**

```
dM/df = r² / (a² β) = β³ / (1 + e cos f)² .                           (0.2.5.1)
```

**Proof.**

**Step 1 (Differentiate Kepler's equation).** `M = E − e sin E` ⇒
`dM/dt = (1 − e cos E) · dE/dt`. We do not need `E` explicitly; we work via the
chain `dM/df = (dM/dt) / (df/dt)`.

**Step 2 (Apply third-law / angular-momentum identities).** By definition of mean
motion, `dM/dt = n` (M increases at the constant rate n). By Theorem 0.2.3,
`df/dt = h / r²`. Hence

```
dM/df = (dM/dt) / (df/dt) = n / (h/r²) = n r² / h .                   (0.2.5.2)
```

**Step 3 (Substitute Theorem 0.2.4).** By (0.2.4.1) `h = n a² β`. So
`n r² / h = n r² / (n a² β) = r² / (a² β)`.

**Step 4 (Substitute Theorem 0.2.2).** By (0.2.2.1) `r = a β²/(1 + e cos f)`. So
`r² / a² = β⁴ / (1 + e cos f)²`, and therefore

```
dM/df = β⁴ / [(1 + e cos f)² · β] = β³ / (1 + e cos f)² .             (0.2.5.3)
```

This establishes both forms in (0.2.5.1). ∎

**Remark 0.2.5.4 (Comparison with main-doc Eq (4.1)).** The main derivation
document at `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` line 164
states the same relation but with the intermediate form `β r²/a²` (β in the
numerator) — that intermediate is **incorrect**; the correct form is `r²/(a²β)`
(β in the denominator). The error is documented as **D-2** in the audit log. The
final closed form `β³/(1+e cos f)²` in the main doc IS correct; only the
intermediate algebra has the typo. **This proof supersedes that proof.**

**Alignment to SGP4 (critical).** The dM/df relation in Theorem 0.2.5 is derived
under fixed-element assumptions: `(a, e, β)` are the constants of motion of an
unperturbed orbit. Phase 1+ of this corpus uses (0.2.5.1) to convert `∫ X dM` into
`∫ X (β³/(1 + e cos f)²) df` for orbit-averaged drag integrals.

In SGP4, the elements `(a, e)` decay due to drag and oscillate due to J₂
short-period periodics, so they are NOT constants over a single orbit. The
substitution dM → df therefore inherits the **secularity hypothesis**:

> The integrand `X(f, a, e, …)` is averaged over `M ∈ [0, 2π)` **with the
> elements held fixed at their current osculating values**. The average is then
> assigned to the *secular rate* of whatever orbital element X is the variational
> derivative of.

This is the textbook variation-of-parameters formulation; cf. [BATT99 §10.3]
("First-Order Perturbation Theory"). It is consistent **provided** the element
drift per period is much smaller than the element value (verified for SGP4 in
§0.2.2 alignment remark).

A consequence: orbital averages of `r(f)`, `ρ(r(f))`, etc. computed via (0.2.5.1)
implicitly project the integrand onto the **secular part** of its time evolution
— short-period and long-period harmonics in `(a, e)` themselves are dropped at
this step. They are reintroduced separately as the J₂ short-period corrections
(`short_period.h`) and the J₃ long-period corrections (`xlcof / aycof`, §14 / Phase 7).
**No claim** is made by Theorem 0.2.5 that the substitution captures these
oscillatory corrections; alignment-to-SGP4 means we use (0.2.5.1) only for the
secular component and integrate the others via separate generating-function pull-backs.

### Definition 0.2.6 (Orbit average over mean anomaly)

For a measurable function `X : [0, 2π) → ℝ` along an orbit, the **orbit average
over the mean anomaly** is

```
⟨X⟩_M := (1/(2π)) ∫₀^{2π} X(M) dM .                                  (0.2.6.1)
```

Equivalently, by Theorem 0.2.5, when X is written as a function of f,

```
⟨X⟩_M = (1/(2π)) ∫₀^{2π} X(f) · (β³ / (1 + e cos f)²) df .            (0.2.6.2)
```

The two forms are interchangeable; (0.2.6.2) is the working form used throughout the
Lane / SGP4 drag derivations because the integrand factors `(1 + e cos f)^p` and
`(1 − η cos f†)^m` (the latter via the Lane substitution, Phase 11) are naturally
expressed in f.

**Remark 0.2.6.3 (Time average equals M-average).** Because `Ṁ = n` is constant,
the orbit average over physical time `t` over one period equals the M-average:
`(1/T) ∫₀^T X dt = (1/(nT)) ∫₀^{nT} X dM = (1/(2π)) ∫₀^{2π} X dM = ⟨X⟩_M`.

**Alignment to SGP4.** SGP4 invokes Definition 0.2.6 in **two distinct contexts**
that the formal definition does not distinguish:

1. **Secular-rate context** (used by Phases 1-6 in this corpus). The orbit-average
   is evaluated **once at epoch** with elements `(a₀'', e₀, i₀)` taken from the
   Brouwer-recovered reference Keplerian. The result is a *constant* secular rate
   that multiplies time `τ = t − t₀`. Subsequent element drift (drag decay, J₂
   secular drift) is propagated separately by appending more time-polynomial terms
   (the D₂-D₄ coefficients in Phases 6-7). The orbit-average operator itself is
   *not* re-evaluated; the rates are frozen at their epoch values.

2. **Instantaneous-osculating context** (used by Phases 7+ for long-period
   periodics). The orbit-average is evaluated at the *current* osculating elements
   at the propagation step, then projected onto a specific harmonic (e.g. the J₃
   long-period kernel that produces `xlcof / aycof`).

The two contexts assign different meanings to "the elements held fixed during the
average." In context 1, "fixed" means at epoch. In context 2, "fixed" means at
the current step.

**Each phase that invokes the orbit-average operator must declare which context
applies.** The corpus is wrong by silent confusion if a context-1 result is used
in a context-2 role or vice versa.

**Secularity hypothesis statement (used in both contexts).** The orbit-average
is meaningful only if the integrand `X(f, a, e, …)` is well-approximated by its
value at fixed elements over one orbital period. Numerically: the per-period
element drift `Δa, Δe` must satisfy `|Δa|·∂X/∂a ≪ |X|` and analogously for `e`.
For SGP4 LEO drag, `|Δa/a| ∼ 10⁻⁵` per period and `∂X/∂a ≲ X/a` for the
integrands of interest, so the residual is `≲ 10⁻⁵` per period — well below the
target accuracy.

---

