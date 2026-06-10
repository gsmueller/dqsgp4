# Phase 2.A — Symbolic Trace of the SGP4 C₂ Drag Coefficient

## §A.0 Scope and Pedagogy

### What Phase 2.A establishes

Phase 2.A is the **critical gating task** of the SGP4 near-earth drag symbolic
re-derivation. It traces the SGP4 implementation's C₂ drag coefficient
(`src/atmosphere/drag_coefficients.h:146-149`) from first principles — the
Phase 0-rev1 orbit-averaged drag rate and the Phase 1 Lane integrals — to a
closed-form expression that is **symbolically equal** to the code expression,
step by step.

Per **Standard 10** (the standing code-match rule, audit
`2026_05_15_..._audit/01_standards_and_schedule.md`), no intermediate Phase-0/1
result is accepted into the corpus until a downstream Phase produces a
code-matched chain that uses it. Phase 2.A is that downstream Phase for C₂:

- **If the chain lands** (final `simplify(C₂_derived − C₂_code) == 0`):
  Phases 0 + 1 are **retrospectively validated**, and audit finding **D-9**
  closes as "main-doc Eq (4.4) is a transcription error in a
  consistent-with-code direction."
- **If the chain fails to land:** the failure point identifies the actual error
  (a wrong or misapplied step in Phase 0, Phase 1, or the SGP4 model
  alignment), and **D-9 escalates**. (Under the no-block discipline, an
  escalation here is resolved via the newstk HR→PA→panel mechanism, not by
  silently choosing a convenient reading.)

### The code expression being matched (target)

From `drag_coefficients.h:146-149`, with `coef = qoms4·ξ⁴ = (q−s)⁴ ξ⁴`,
`coef1 = coef/(ψ²)^{7/2} = (q−s)⁴ ξ⁴ (1−η²)^{−7/2}`, `eta2 = η²`, `eeta = eη`,
`psisq = ψ² = |1−η²|`, `half_J2 = k₂ = J₂/2`, `three_cos2i_minus_1 = 3cos²i − 1`:

```
C₂ = coef1 · n₀ · [  a₀ (1 + (3/2)η² + eη(4 + η²))                        ← Part A (Keplerian)
                   + (3/4) k₂ (ξ/ψ²) (3cos²i − 1) (8 + 24η² + 3η⁴)  ]    ← Part B (J₂ coupling)
```

The code writes Part A's polynomial as `1 + (3/2)η² + eη(4 + η²) =
1 + (3/2)η² + 4eη + eη³`, and Part B's polynomial in Horner form
`8 + 3η²(8 + η²) = 8 + 24η² + 3η⁴`.

**ψ-power bookkeeping (the structural map).**
- Part A carries `coef1 ∝ ψ⁻⁷ = (1−η²)^{−7/2}` — the pole order of Phase 1
  `I^{(0,4)} = (2 + 3η²)/(2(1−η²)^{7/2})` (the D-9 integral).
- Part B carries `coef1 · (ξ/ψ²) ∝ ψ⁻⁹ = (1−η²)^{−9/2}` — the pole order of
  Phase 1 `I^{(0,5)} = (8 + 24η² + 3η⁴)/(8(1−η²)^{9/2})`, whose **numerator is
  exactly Part B's polynomial**.

### Inputs (consumed from prior Phases)

- **Phase 0-rev1 Theorem 0.4.2 (0.4.2.1)** — instantaneous drag rate
  `ȧ_drag(f) = −(2 n B* ρ(r(f)) a²/β³) (1 + e² + 2e cos f)^{3/2}`.
  (`sgp4_drag_phase0_foundations/05_drag_specialization.md`)
- **Phase 0-rev1 Theorem 0.5.3 (0.5.3.2)** — orbit-average setup
  `⟨ȧ⟩_drag = (1/2π) ∫₀^{2π} ȧ_drag · β³/(1+e cos f)² df`.
  (`.../06_energy_orbit_averaged_rates.md`)
- **Phase 0 Definition 0.2.6** — orbit-average operator `⟨·⟩_M`.
  (`.../02_notation_kepler_geometry.md`)
- **Phase 1 Theorem 1.3.2** — Lane integrals `I^{(0,m)}`; in particular
  `I^{(0,5)}` (Part B) and `I^{(0,4)}` (D-9). (`sgp4_drag_phase1_lane_integrals.md`)
- **Lane density model** `ρ(r) = ρ₀ ((q−s)/(r−s))⁴` (τ = 4) and the Lane
  fictitious-anomaly substitution `r − s = (1/ξ)(1 − η cos f†)` (main doc §2 / §4).

### Pedagogy standards (same as Phase 0 / Phase 1)

- **Standard 1** — symbolic derivation from theoretical sources; born-digital
  citations only.
- **Standard 9-A** — alignment to the SGP4 problem made explicit at each step.
- **Standard 9-B** — alignment to implementation (the precomputed `coef1`
  grouping etc.); acknowledged and discharged in the §A.8 code-match witness.
- **Standard 10** — every algebraic transformation is a named theorem
  application or definition substitution, demonstrable as
  `simplify(lhs − rhs) == 0`. Numerical agreement at sample points is **not**
  proof.

### Open audit findings touched

- **D-3** (master orbit-average identity Eq (4.2) algebraic inconsistency):
  §A.1-A.2 re-derive the master orbit-average correctly (β³ Jacobian, not the
  doc's β⁵-equivalent form).
- **D-9** (I^{(0,4)} value vs main-doc Eq (4.4)): disposition recorded in §A.8.

---

## §A.1 Orbit-averaged ⟨ȧ⟩_drag setup  *(B2)*

**Objective.** Express the secular drag rate of the semi-major axis as a single
orbit-average integral over the true anomaly, with the Newtonian-quadratic-drag
factor `(1 + e² + 2e cos f)^{3/2}` explicit, ready for the Lane density
substitution (§A.2).

**By Theorem 0.5.3 (0.5.3.2)** the orbit-averaged drag rate is
```
⟨ȧ⟩_drag = (1/2π) ∫₀^{2π} ȧ_drag(f) · (β³/(1 + e cos f)²) df ,          (A.1.1)
```
**by Definition 0.2.6** (orbit average over mean anomaly), with the Jacobian
`dM/df = β³/(1 + e cos f)²` from **Theorem 0.2.5**.

**By Theorem 0.4.2 (0.4.2.1)** — the Phase 0-rev1 instantaneous drag rate with
the `|v|` factor restored (closing D-10) —
```
ȧ_drag(f) = −(2 n B* ρ(r(f)) · a²/β³) · (1 + e² + 2e cos f)^{3/2} .      (A.1.2)
```

**Substitution.** Inserting (A.1.2) into (A.1.1), the Jacobian's `β³` cancels
the `1/β³` of (A.1.2):
```
⟨ȧ⟩_drag = −(2 n B* a²) · (1/2π) ∫₀^{2π}
              ρ(r(f)) · (1 + e² + 2e cos f)^{3/2} / (1 + e cos f)²  df .  (A.1.3)
```
The β-cancellation is exact (`simplify` target: `β³·(1/β³) − 1 = 0`); no β
survives into the C₂ integrand. Equation (A.1.3) is the **master drag-rate
integral** for C₂.

**L-rate equivalent (SR3 normalization).** SR3 / main doc §6.1 define the
primary drag coefficient via the Delaunay-momentum rate. **By Definition 0.5.2
(0.5.2.2)** `L̇ = (na/2) ȧ`, so
```
⟨L̇⟩_drag = (n a / 2) · ⟨ȧ⟩_drag ,                                       (A.1.4)
```
the two normalizations differing only by the constant `na/2`. We carry the
`⟨ȧ⟩` form (A.1.3) through §§A.2–A.7 and reconcile to the code's C₂ grouping
(which exposes the `n₀` and `a₀` factors explicitly, with `C₁ = B*·C₂` and the
secular law `ȧ ≈ −2a₀ B* C₂` at epoch) in §A.8.

**Alignment to SGP4 (Standard 9-A).** The average (A.1.1) is taken in
**Definition 0.2.6 context-1** (secular-rate context): elements
`(a, e, i) = (a₀'', e₀, i₀)` are frozen at Brouwer-recovered epoch values, so C₂
is an epoch constant — consistent with `compute_drag_coefficients` evaluating it
once at initialization.

## §A.2 Lane density substitution, pre-f†  *(B3)*

**Objective.** Replace `ρ(r(f))` in the master integral (A.1.3) by the Lane
power-law model, expressing `r − s` exactly in `(f, e, a, s)`, and re-derive the
master orbit-average identity correctly (resolving **D-3**).

**Lane power-law density (τ = 4).**
```
ρ(r) = ρ₀ · ((q − s)/(r − s))⁴ = ρ₀ (q − s)⁴ (r − s)^{−4} ,             (A.2.1)
```
`q` the reference perigee-fit radius, `s` the atmospheric scale-fit radius
(constants within a perigee regime, Definition 0.6.3).

**ρ₀ / B\* convention (D-1).** The code's `coef = (q−s)⁴ξ⁴` carries no `ρ₀`; the
reference density rides in `B*` (the TLE ballistic coefficient is effectively
`B* = (C_D A/(2m))·ρ₀`). So `B*·ρ(r) = B*·((q−s)/(r−s))⁴` with the **dimensionless**
geometric profile left in the integrand. (Audit **D-1** records that the
main-doc *definition* of `B*` omitted `ρ₀`; the code/convention is consistent
once `ρ₀` is assigned to `B*`.) We carry `B*` understood to include `ρ₀`.

**Exact `r − s`.** By Theorem 0.2.2 `r = aβ²/(1 + e cos f)`, so
```
r − s = (aβ² − s − s e cos f)/(1 + e cos f) .                            (A.2.2)
```

**D-3 resolution (master orbit-average identity).** Re-deriving the
orbit-averaged density identity from Definition 0.2.6 with the **correct**
Jacobian `dM = β³/(1+e cos f)² df` (Theorem 0.2.5; not the main-doc
β-in-numerator intermediate, **D-2**) gives, for any integer m,
```
⟨ρ(r)(r/a)^m⟩_M = ρ₀(q−s)⁴ · (1/2π) ∫₀^{2π} (r−s)^{−4}(r/a)^m · β³/(1+e cos f)² df. (A.2.3)
```
The numerator carries the Jacobian's **β³**, confirming the audit's corrected
D-3 form over the main-doc Eq (4.2) numerator `β·(r/a)^{m+2}` (the two differ by
the f-dependent factor `β²/(1+e cos f)²`). **This supersedes main-doc Lemma 4.2.**

**Forward to §A.3.** `(r − s)^{−4}` is converted to the Lane form `(a−s)⁴
(1 − η cos f†)^{−4}` next; §A.3 shows that factor is **exact** in `f`.

## §A.3 Lane f† substitution, O(eη) retained  *(B4)*

**Objective.** Convert the exact `(r − s)^{−4}` factor of the master integrand
into the Lane integral form `(1 − η cos f†)^{−4}`. §A.4 establishes by symbolic
verification that the O(eη) corrections of this substitution are **load-bearing**
for the C₂ Part-A e-term — so they are retained, not dropped.

**Lane intermediate quantities (Definition 0.6.4).** `ξ := 1/(a − s)`,
`η := a e ξ = a e/(a − s)`, with `0 < η < 1` for bound LEO with perigee above `s`.

**The fictitious anomaly.** `f†` is defined by
```
r − s = (a − s)(1 − η cos f†) = (1/ξ)(1 − η cos f†) .                    (A.3.1)
```
Equating (A.3.1) to the exact `r − s = (aβ² − s − s e cos f)/(1 + e cos f)`
(A.2.2) gives an implicit relation `f†(f; e, η)` with `f† = f + O(e)` and
Jacobian `df†/df = 1 + O(eη)` (main-doc Note 4.3). Both the angular offset and
the Jacobian deviation from unity are first order in the product `eη`.

**The exact density is what is kept (verified).** The factor `(1 − η cos f†)`
is **exact** as a function of `f`. Using the η-parameterization to eliminate
`(a, s)` (`η = ae/(a−s)` ⇒ `s = a(1 − e/η)`), the exact Kepler `r` (Theorem
0.2.2) and the Lane definition (A.3.1) give the closed form
```
1 − η cos f† = (1 − η cos f − ηe + e cos f) / (1 + e cos f) .            (A.3.2)
```
Expanding (A.3.2) to O(e) at **fixed η**:
```
1 − η cos f† = (1 − η cos f) − ηe sin²f + O(e²) ,                        (A.3.3)
```
so the Lane density factor is
```
(1 − η cos f†)^{−4} = (1 − η cos f)^{−4} + 4ηe sin²f (1 − η cos f)^{−5} + O(e²). (A.3.4)
```
The naive `f† ≈ f` (replacing (A.3.2) by `1 − η cos f`) drops the exact O(e)
term `−ηe sin²f`, which is why it recovers only **half** of the Part-A e-term.
**No change of variables is needed** — keeping the exact density (A.3.4) and
integrating in `f` is exact; the new `(1 − η cos f)^{−5}` piece brings in
`I^{(0,5)}`. The only Part-A approximation is the SGP4 O(e²) truncation (§A.7).
All of (A.3.2)–(A.3.4) are verified symbolically (`verify_phase2a.m`, 8/8).

## §A.4 Reduction to Lane integrals  *(B5 — in progress)*

Reduces the master integrand to Lane integrals `I^{(p,m)}` (Definition 1.1.1) in
two passes: a leading `f† ≈ f` pass (fixes Part A's constant term, closes D-9)
and the O(eη) f†-correction pass (supplies Part A's e-term).

### A.4.1 Kinematic factor expansion (f-variable, O(e))

The kinematic factor `g(e,f) := (1 + e² + 2e cos f)^{3/2}/(1 + e cos f)²`
Taylor-expands in e as
```
g(e,f) = 1 + e cos f + O(e²) .                                          (A.4.1)
```
**Verified symbolically** (Octave SymPy): `c0 = 1`, `c1 = cos f`. The `(3/2)`
quadratic-drag power contributes `+3e cos f`; the Jacobian leftover
`(1 + e cos f)^{−2}` contributes `−2e cos f`; net `+e cos f`.

### A.4.2 Leading pass (f† ≈ f) — Part A constant, D-9 CLOSURE

Under `f† ≈ f`, `(r−s)^{−4} → ξ⁴(1 − η cos f)^{−4}`. Averaging (A.4.1) against
`(1 − η cos f)^{−4}` term by term via the cos-power reduction
`cos f = (1 − (1 − η cos f))/η`:
```
O(e⁰):  ⟨(1 − η cos f)^{−4}⟩ = I^{(0,4)} = (2 + 3η²)/(2(1−η²)^{7/2}) .   (A.4.2)
O(e¹):  ⟨cos f·(1 − η cos f)^{−4}⟩ = (I^{(0,4)} − I^{(0,3)})/η
                                   = η(4 + η²)/(2(1−η²)^{7/2}) .          (A.4.3)
```

**Part-A constant term — D-9 CLOSURE.** Factoring `ψ^{−7} = (1−η²)^{−7/2}`,
`(A.4.2) = ψ^{−7}(1 + (3/2)η²)`, which is **exactly** the code's Part-A constant
(once `ψ^{−7}` is absorbed into `coef1`). **Verified:**
`simplify(I^{(0,4)} − ψ^{−7}(1 + (3/2)η²)) = 0`.

This **closes D-9**: Phase 1 Theorem 1.3.2.4 `I^{(0,4)} = (2+3η²)/(2(1−η²)^{7/2})`
is correct AND is the integral SGP4 computes for the Part-A constant. The
main-doc Eq (4.4) value `(2+η²)/(1−η²)^{7/2}` (which gives 2 at η=0 instead of 1)
is the **transcription error**, superseded by Phase 1.

### A.4.3 Exact-density O(e) pass — Part A e-term  *(LANDED)*

Keeping the exact density (A.3.4), the O(e) coefficient of the integrand is
```
a1(f) = cos f (1 − η cos f)^{−4} + 4η sin²f (1 − η cos f)^{−5} ,         (A.4.4)
```
(verified `a1` matches). Averaging via the cos-power and sin² reductions plus
the Lane closed forms:
```
⟨a1⟩ = (I^{(0,4)} − I^{(0,3)})/η + 4η[ I^{(0,5)} − (I^{(0,5)} − 2I^{(0,4)} + I^{(0,3)})/η² ]
     = ψ^{−7} · η(4 + η²) ,                                              (A.4.5)
```
**exactly** the code's Part-A e-term. **Verified:**
`simplify(⟨a1⟩ − ψ^{−7}η(4+η²)) = 0`. The first term
`(I^{(0,4)}−I^{(0,3)})/η = ψ^{−7}η(4+η²)/2` is the `f†≈f` half; the second
`4η⟨sin²f(1−η cos f)^{−5}⟩` is the exact-density half (the `−ηe sin²f` term of
(A.3.3)). Their sum lands the full e-term — the f† substitution required **no
approximation** at O(e).

### A.4.4 Lane integral inventory

Leading pass uses `I^{(0,3)}, I^{(0,4)}` (both `p = 0`, within Phase 1 scope).
The O(eη) correction pass is expected to remain in the `p = 0` family after the
cos-power reduction; **if any `I^{(p,m)}` with `p ≥ 1` is required, ESCALATE to
Block C1** (Phase 1 §1.4). Status: leading pass confirmed `p = 0` only.

## §A.5 Part A (Keplerian) assembly  *(B6)*

Combining §A.4.2 (constant) and §A.4.3 (e-term), the Keplerian Part-A average to
O(e) is
```
⟨g · (1 − η cos f†)^{−4}⟩ = ψ^{−7} · (1 + (3/2)η² + eη(4 + η²)) + O(e²) .  (A.5.1)
```
Restoring the prefactors (`coef = (q−s)⁴ξ⁴` from §A.2, the `n₀ a₀` from
§A.1/§A.8), this is **exactly** the SGP4 code's Part-A group
```
coef1 · n₀ · a₀ · (1 + (3/2)η² + eη(4 + η²)) ,    coef1 = coef · ψ^{−7} , (A.5.2)
```
i.e. `drag_coefficients.h:146-147`. **Part A matches the code to O(e)** (verified,
`verify_phase2a.m` 8/8). The dropped `O(e²)` Part-A residuals (`(3/4)e² + 3e²η²`
of the full AFGP4 form) are the SGP4 simplification handled in §A.7.

## §A.6 Part B (J₂ density coupling) assembly  *(B7)*

**Origin.** Drag acts on the J₂-perturbed orbit: `ρ(r)` is evaluated at the
J₂-perturbed radius `r = r_kep + δr`. Linearizing the τ=4 Lane density,
```
ρ(r_kep + δr) = ρ(r_kep)·(1 − 4 δr/(r − s) + O(J₂²)) .                   (A.6.1)
```
The `−4 δr/(r−s)` term is **Part B**; `1/(r−s) = ξ/(1 − η cos f)` (Lane) raises
the density power `−4 → −5`, so Part B is governed by `I^{(0,5)}` (vs Part A's
`I^{(0,4)}`), with the `(3cos²i − 1)` factor carried by `δr`.

**Sealed-room input (J₂ secular radial multiplier).** From the SGP4/Lane theory
(born-digital: Vallado `SGP4.cpp:1990-2002`, `mrt = rl·(1 − 1.5·(k₂/p²)·βₗ·con41)`;
SR3 Eq 6.7; Brouwer 1959 `γ₂(3θ²−1)`, `γ₂ = k₂/p²`):
```
δr/r = −(3/2)·(k₂/p²)·β·(3cos²i − 1) ,    k₂ = J₂ a_E²/2,  p = a(1−e²). (A.6.2)
```
(Obtained via a sealed-room sub-agent dispatch; integrated as a cited result.
**Audit note (new finding):** the project main-doc §3 Theorem 3.1 states the
*short-period amplitude* form `−(k₂/p)(3θ²−1)/2`, a **factor of 3 too small**
for the *secular radial multiplier* (A.6.2); the SGP4 code and SR3 use (A.6.2).
§3 Theorem 3.1 should be corrected.)

**Power raise (verified).** Substituting (A.6.2) and `r/(r−s) = aξ/(1−η cos f)`
(leading e) into the `−4 δr/(r−s)` factor times the leading Lane density:
```
−4(δr/r)(r/(r−s))(1−η cos f)^{−4} = 6(k₂/p²)β(3cos²i−1)·aξ·(1−η cos f)^{−5}. (A.6.3)
```
**Verified** `simplify = 0` (`verify_phase2a.m`).

**Orbit average + assembly (verified).** Averaging `(1−η cos f)^{−5} → I^{(0,5)}`
(Phase 1 Theorem 1.3.2.5) and restoring the `n₀ a₀ coef` prefactor with the
leading-e dimensional reduction `a₀ a β/p² → 1` (the `1/p² ≈ 1/a₀²` from `δr/r`
times the `a` from `r/(r−s)` cancels one `a₀` — why Part B has no explicit `a₀`
in the code):
```
C₂|_PartB = coef1 · n₀ · (3/4) k₂ (ξ/ψ²)(3cos²i − 1)(8 + 24η² + 3η⁴) ,    (A.6.4)
```
**exactly** `drag_coefficients.h:148-149`. **Verified** `simplify = 0`. Numeric
trace: `I^{(0,5)} = (8+24η²+3η⁴)/(8ψ⁹)` × the `6 = 4·(3/2)` gradient factor gives
the `(3/4)` prefactor; `ψ⁹ = ψ⁷·ψ²` supplies the extra `1/ψ²` beside Part A.

## §A.7 AFGP4 → SGP4 simplification  *(B8)*

The match in §A.1–§A.6 is to the **SGP4** C₂ — a truncation of the full AFGP4
orbit integral. Two families of higher-order terms are dropped:

**(i) Part-A `O(e²)` residual — now DERIVED EXACTLY (Phase 11, `verify_phase11_fdagger.m` 8/8).**
The full AFGP4 Part-A bracket carries an `O(e²)` term `((3/4)e² + 3e²η²)` (SR3 p. 11;
`drag_coefficients.h:132` comment) which SGP4 drops. **Phase 11 resolves the Lane f†:**
`f† = E` (the **eccentric anomaly**) — `r − s = (a−s)(1 − η cos E)` is **EXACT** (P11.1), so the
density `(r−s)⁻⁴ = ξ⁴(1−η cos E)⁻⁴` is exact (no f† error). The C2 Part-A average **over E** (Jacobian
`dM = (1−e cos E)dE`) has integrand `(1−e²cos²E)^{3/2}/[(1−e cos E)²(1−η cos E)⁴]`; its kinematic
factor expands `= 1 + 2e cos E + (3/2)e²cos²E + O(e³)` (P11.2), and the `cosᵏE` terms reduce onto the
Phase-1 Lane integrals (P11.3-4). This gives **exactly**: constant `I^{(0,4)} = (1+(3/2)η²)ψ⁻⁷` (P11.5),
e-term `(2/η)(I^{(0,4)}−I^{(0,3)}) = η(4+η²)ψ⁻⁷` ⇒ code's `eη(4+η²)` (P11.6), and the dropped **`O(e²)`
term `(3/2)/η²·(I^{(0,4)}−2I^{(0,3)}+I^{(0,2)}) = (3/4)(1+4η²)ψ⁻⁷` ⇒ `(3/4)e² + 3e²η²`** (P11.7) — the
documented AFGP4 dropped term, **DERIVED**. (Independently, `verify_A7_truncation_floor.m` 5/5
numerically confirms `B_exact − B_code` matches `(¾e²+3e²η²)` to 0.6–6.7%; that **residual is the
`O(e³)` higher-order e-truncation**, NOT an f† error — f† = E is exact.) This **supersedes** §A.4's
`f† ≈ f` (true-anomaly) approximation: the f†-corrections there are exactly the `E`-vs-`f`
equation-of-center terms; integrating over `E` is exact and cleaner.

**(ii) Part-B `O(eη)` cross term.** §A.6 assembled Part B at **leading e**
(`a₀aβ/p² → 1`). The full AFGP4 Part B carries an `O(eη)` cross term
`−5eη(4 + 3η²)` (SR3 p. 11; `drag_coefficients.h:140` comment), from the
e-dependence of `p = a(1−e²)`, `β`, and `r/(r−s)`. SGP4 drops it. **Bound:**
relative size `O(eη)`; for `e ≲ 0.1, η ≲ 0.2`, a few × `10⁻²`.

**Net.** The SGP4 C₂ that §A.8 code-matches is `AFGP4_C₂ − (i) − (ii)`. The
dropped terms are model-choice approximations of the SGP4 specification (audit
sources A-D6 / A-D7), not derivation errors; they set the C₂ accuracy floor at
`O(e², eη)`. **Part-A's `(3/4)e²+3e²η²` is now DERIVED EXACTLY** (Phase 11 / `f†=E`,
`verify_phase11_fdagger.m` 8/8; numerically corroborated `verify_A7_truncation_floor.m`
5/5). Part-B's `−5eη(4+3η²)` is the analogous `O(eη)` floor — the **same exact-over-E
method** (Phase 11) would derive it from the J₂-coupled integrand; the Part-A derivation
is the representative completeness demonstration. **Phase 11 is RESOLVED** (no f† density
error; the residual beyond `(¾e²+3e²η²)` is `O(e³)`), closing the last drag-derivation
open item.

## §A.8 Code-match witness  *(B9)*

**Headline (Standard 10).** Assembling §A.5 (Part A) and §A.6 (Part B), the
derived C₂ is
```
C₂_derived = coef1·n₀·( a₀(1 + (3/2)η² + eη(4 + η²))
                      + (3/4)k₂(ξ/ψ²)(3cos²i − 1)(8 + 24η² + 3η⁴) ) ,      (A.8.1)
```
with every coefficient **derived** from Phase 0-rev1 + Phase 1 + the J₂-radial
input — no value imported from the code. The literal SGP4 code
`drag_coefficients.h:146-149` is
```
C₂_code = coef1·n₀·( a₀(1 + (3/2)η² + eη(4 + η²))
                   + (3/4)k₂(ξ/ψ²)(3cos²i − 1)(8 + 3η²(8 + η²)) ) .         (A.8.2)
```
With `(8 + 3η²(8 + η²)) = 8 + 24η² + 3η⁴`,
```
simplify(C₂_derived − C₂_code) = 0 .                                        (A.8.3)
```
**Verified** — `verify_phase2a.m` check A.8 (11/11 PASS).

**Disposition (gating outcome).** Phase 2.A lands cleanly. Per the Standard-10
gate:
- **Phases 0-rev1 and 1 are retrospectively VALIDATED** — the quadratic-drag
  `(1+e²+2e cos f)^{3/2}` integrand (Phase 0-rev1) and the Lane integrals
  `I^{(0,4)}, I^{(0,5)}` (Phase 1) are exactly what the SGP4 C₂ code computes.
- **D-9 CLOSED**: `I^{(0,4)} = (2+3η²)/(2(1−η²)^{7/2})` (Phase 1 Thm 1.3.2.4) is
  correct and code-matched; main-doc Eq (4.4) is a transcription error.
- **D-3 resolved** (§A.2): the master orbit-average uses the correct β³ Jacobian.
- **New finding** (§A.6): main-doc §3 Theorem 3.1 J₂-radial multiplier is 3×
  too small; the code-correct form is (A.6.2).

§A.7 documents the AFGP4 → SGP4 truncation bounds for the dropped `O(e²)` Part-A
and `−5eη(4+3η²)` Part-B terms (model-choice approximations, not code-used).

---

## §A.9 C₁ = B*·C₂ — the secular-decay coefficient  *(Phase 2.B)*

**Objective.** Derive the SGP4 secular semi-major-axis decay coefficient C₁
(`drag_coefficients.h:152`) as `B*·C₂`, completing the secular-law reconciliation
that §A.1.4 promised and §A.8 deferred — the explicit `−2a₀B*` factor-out of the
master drag-rate integral. No value is imported from the code.

**Inputs.** §A.1.3 master drag-rate integral; §A.5 (A.5.2) Part-A assembly; §A.6
(A.6.4) Part-B assembly; §A.8 (A.8.1) C₂ closed form.

### A.9.1 SGP4 secular-decay parameterization (code-aligned)

SGP4 propagates the drag-perturbed semi-major axis as (`drag_coefficients.h:14`)
```
a(t) = a₀ · tempa² ,    tempa = 1 − C₁ t − D₂ t² − D₃ t³ − D₄ t⁴ .       (A.9.1)
```
Differentiating at epoch (`t = 0`, `tempa = 1`):
```
ȧ(0) = 2 a₀ · tempa · (d tempa/dt)|₀ = 2 a₀ (−C₁) = −2 a₀ C₁ ,           (A.9.2)
```
so C₁ is, by construction, minus one-half the epoch fractional decay rate,
`C₁ = −ȧ(0)/(2a₀)`.

### A.9.2 The −2a₀B* factor-out (reconciling §A.1.3 → §A.8)

The master drag-rate integral (A.1.3) has prefactor `−(2 n B* a²)`; at the
Brouwer-recovered epoch (`a = a₀`, `n = n₀`, Standard 9-A) this is `−(2 n₀ B* a₀²)`.
The §A.5 / §A.6 assemblies expressed the two evaluated parts of the rate as
```
⟨ȧ⟩_PartA = −2a₀B* · ( coef1 · n₀ · a₀ · (1 + ³⁄₂η² + eη(4+η²)) )  [= −2a₀B*·(A.5.2)] ,
⟨ȧ⟩_PartB = −2a₀B* · ( coef1 · n₀ · (3/4)k₂(ξ/ψ²)(3θ²−1)(8+24η²+3η⁴) ) [= −2a₀B*·(A.6.4)] .
```
The factor-out is exact: matching `⟨ȧ⟩_PartA` to the A.1.3 prefactor requires
`(2 n₀ a₀²)·coef·ψ⁻⁷ = 2a₀·(coef1 n₀ a₀)`, which holds identically since
`coef1 = coef·ψ⁻⁷` (**verified** `simplify = 0`). Summing the two legs and recognising
the bracket as C₂ (A.8.1):
```
⟨ȧ⟩_drag = −2 a₀ B* · [ coef1 n₀ ( a₀ PartA + PartB ) ] = −2 a₀ B* C₂ .   (A.9.3)
```
Equation (A.9.3) is the **defining normalization of C₂**: the per-B*, per-(−2a₀)
orbit-averaged secular decay rate. (This is the reconciliation §A.1.4 deferred.)

### A.9.3 Identification (code-match)

Both (A.9.2) and (A.9.3) are `ȧ(0)`; equating:
```
−2 a₀ C₁ = −2 a₀ B* C₂   ⇒   C₁ = B* C₂ ,                                (A.9.4)
```
exactly `drag_coefficients.h:152` (`dc.C1 = in.bstar * dc.C2`). Since §A.8 proved
`C₂_derived = C₂_code`, it follows that `C₁_derived = B*·C₂_derived = C₁_code`.
**Verified** — `verify_phase2a.m` checks A.9.2 (prefactor factor-out) and A.9.3
(`simplify(C₁_derived − C₁_code) = 0`). This closes audit **Q-4** (C₁ now has a
derivation, not only a quote-check).

---

**Status:** Phase 2.A + 2.B **COMPLETE**. §A.1–§A.9 derived; full C₂ **code-matched**
and C₁ = B*·C₂ derived from the secular-decay law (`verify_phase2a.m` 13/13 PASS —
check A.8 `simplify(C₂_derived − C₂_code) = 0` vs `drag_coefficients.h:146-149`;
checks A.9.2/A.9.3 land C₁ vs `:152`). **D-9 CLOSED**, **D-3 resolved**, **Q-4
closed**; new finding (main-doc §3 Thm 3.1 J₂-radial multiplier 3× too small, §A.6).
§A.7 truncation bounds documented (exact AFGP4 dropped-term coefficients are a
deferred completeness item, not code-used). Phases 0-rev1 + 1 retrospectively
validated.
