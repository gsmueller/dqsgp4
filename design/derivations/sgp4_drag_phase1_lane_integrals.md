# Phase 1 — Lane Integrals via Residue Calculus

## §1.0 Scope and Pedagogy

This is the **Phase 1** document of the symbolic re-derivation of the SGP4
near-earth drag coefficients per the multi-pass audit at
`design/audit/2026_05_15_sgp4_drag_derivation_full_audit.md`.

### What Phase 1 establishes

Phase 1 evaluates a family of orbit-averaged integrals — the **Lane integrals**

```
I^{(p,m)}(η, e) := (1/(2π)) · ∫₀^{2π} (1 + e cos f)^p / (1 − η cos f)^m  df    (1.0.1)
```

— for the parameter values needed downstream by Phases 2-10:
- `p ∈ {0, 1, 2, 3}`
- `m ∈ {1, 2, 3, 4, 5, 6}` (with `m = 4` being the principal SGP4 case for
  the Lane power-law density at τ = 4)

The derivation method is **residue calculus on the unit circle**, plus a
**recurrence relation** that propagates the result from `m` to `m+1` via
differentiation. Both approaches yield identical closed forms; we present the
recurrence approach as primary (because it is algebraically cleaner) and the
residue-calculus approach as the alternative verification path (because it is
the standard textbook construction).

### Scope (and explicitly out of scope)

**In scope:** the Lane integrals (1.0.1) for the parameter values above,
expressed as closed-form rational functions of `η, e` over polynomial powers
of `(1−η²)`.

**Out of scope (Phase 11):** the relationship between the integration
variable `f` here and Lane's *fictitious anomaly* `f†` introduced in the main
document's §4 Note 4.3. **PHASE 11 RESOLVED (2026-06-03, `verify_phase11_fdagger.m`
8/8): `f† = E`, the eccentric anomaly** — `r − s = (1/ξ)(1 − η cos f†)` is **EXACT**
with `f† = E` (since `r = a(1−e cos E)` and `η = ae/(a−s)`), so the Lane density
substitution carries **no approximation error**. The Lane integrals are therefore
`I^{(0,m)}(η) = (1/2π)∫(1−η cos E)⁻ᵐ dE` **over `E`**, and the SGP4 orbit-average
`⟨ρX⟩ = (1/2π)∫(1−η cos E)⁻⁴ X (1−e cos E) dE` reduces to them **exactly** (the C2
Part-A derivation, C2-trace §A.7). Phase 1 establishes the integrals **as integrals over
a generic angular variable** under `|η| < 1`; Phase 11 identifies that variable as `E`.
SGP4's `f† ≈ f` (true-anomaly) treatment is the source of the `O(e²)` AFGP4 dropped terms
(the `E`-vs-`f` equation-of-center), **not** a density error.

### Pedagogy

Identical to Phase 0 — three standards: (1) symbolic derivation from
theoretical sources, (9-A) alignment to SGP4-problem, (9-B) alignment to
implementation (acknowledged, deferred). Born-digital citations only —
Wikipedia, Gradshteyn-Ryzhik 8th ed. (born-digital), Battin. No OCR
primary-source cites.

### Born-digital source list for Phase 1

| Cite | Reference |
|---|---|
| **[WIKI-RES]** | Wikipedia: *Residue theorem*; *Residue (complex analysis)*. |
| **[GR8]** | Gradshteyn-Ryzhik, *Table of Integrals, Series, and Products*, 8th ed (2014). Eq. 3.616 series. |
| **[BATT99]** | Battin (1999), *Astrodynamics*. Hansen-Bessel coefficient setup. |
| **[WIKI-HC]** | Wikipedia: *Hansen coefficients* / *Bessel functions of the first kind* — `I^{(p,m)}` are Hansen coefficients in disguise. |

---

## §1.1 Notation and Definitions

| Symbol | Definition | Domain |
|---|---|---|
| `η` | Lane scaled-eccentricity parameter; orbit-element form `η = a₀'' e₀ ξ` (Phase 0 Def 0.6.4.2) | `|η| < 1` for integrability |
| `e` | True eccentricity | `0 ≤ e < 1` |
| `f` | Integration variable (true anomaly or Lane fictitious — see Phase 11 alignment) | `f ∈ [0, 2π)` |
| `I^{(p,m)}` | Lane integral (1.0.1) | depends on `(p, m, η, e)` |

### Definition 1.1.1 (Lane integral)

For non-negative integers `p, m` with `m ≥ 1`, and parameters `η, e` with
`|η| < 1` and `|e| < 1`,

```
I^{(p,m)}(η, e) := (1/(2π)) · ∫₀^{2π} (1 + e cos f)^p / (1 − η cos f)^m  df .  (1.1.1.1)
```

The integral is absolutely convergent under `|η| < 1` (since the denominator
never vanishes on the unit circle) and `|e| < ∞` for finite `p`.

**Remark 1.1.1.2 (Special structure for `e = η`).** When `e = η`, the integrand
factors as `[(1+e cos f)/(1−e cos f)]^... ` only at specific `p, m`; in general
the integral is a non-trivial function of both `e` and `η`. SGP4 has `η = a e ξ`
which is `e · (a ξ)`. With `a ξ = a/(a−s) ≈ 1.01` for typical LEO, `η ≈ 1.01 e`,
so `η/e ≈ 1.01 ≠ 1`. The two are independent parameters in (1.1.1.1).

### Definition 1.1.2 (`I^{(p,m)}` for `p = 0` — pure η-integral)

The `p = 0` case is the **pure η-integral**

```
I^{(0,m)}(η) := (1/(2π)) · ∫₀^{2π}  df / (1 − η cos f)^m .                  (1.1.2.1)
```

The `e`-dependence has dropped out; `I^{(0,m)}` is a function of `η` alone.
This is the case most directly related to the residue-calculus closed forms,
and the case that produces the `ψ^{-(2m-1)} = (1−η²)^{-(2m-1)/2}` factor
appearing in the SGP4 C-coefficient formulas.

---

## §1.2 Recurrence Relation

### Theorem 1.2.1 (Differentiation recurrence)

**Hypotheses.** As in Definition 1.1.2 (`p = 0`, `m ≥ 1`, `|η| < 1`).

**Conclusion.**

```
I^{(0,m+1)}(η) = I^{(0,m)}(η) + (η / m) · dI^{(0,m)}(η) / dη .             (1.2.1.1)
```

**Proof.**

**Step 1 (Differentiate under the integral sign).** From (1.1.2.1):

```
dI^{(0,m)}/dη = (1/(2π)) · ∫₀^{2π} d/dη [(1 − η cos f)^{−m}] df
              = (1/(2π)) · ∫₀^{2π} m cos f · (1 − η cos f)^{−(m+1)} df .   (1.2.1.2)
```

**Step 2 (Algebraic identity for `cos f`).** From `1 − η cos f = 1 − η cos f`,
rearrange:

```
cos f = (1 − (1 − η cos f)) / η .                                          (1.2.1.3)
```

**Step 3 (Substitute into the integrand of (1.2.1.2)).**

```
cos f · (1 − η cos f)^{−(m+1)} = (1/η) · [1 − (1 − η cos f)] · (1 − η cos f)^{−(m+1)}
                              = (1/η) · [(1 − η cos f)^{−(m+1)} − (1 − η cos f)^{−m}] .
                                                                         (1.2.1.4)
```

**Step 4 (Integrate term-by-term).** Applying `(1/(2π)) · ∫₀^{2π} · df` to
both sides of (1.2.1.4):

```
(1/(2π)) · ∫₀^{2π} cos f (1 − η cos f)^{−(m+1)} df = (1/η) · [I^{(0,m+1)} − I^{(0,m)}] .
                                                                         (1.2.1.5)
```

**Step 5 (Combine with (1.2.1.2)).** Substituting (1.2.1.5):

```
dI^{(0,m)}/dη = m · (1/η) · [I^{(0,m+1)} − I^{(0,m)}]
              = (m/η) · [I^{(0,m+1)} − I^{(0,m)}] .                       (1.2.1.6)
```

**Step 6 (Solve for `I^{(0,m+1)}`).** Rearranging (1.2.1.6):

```
I^{(0,m+1)} = I^{(0,m)} + (η/m) · dI^{(0,m)}/dη .                         (1.2.1.7)
```

This is (1.2.1.1). ∎

**Remark 1.2.1.8 (Why this recurrence is useful).** Given the base case
`I^{(0,1)}` (which is a classical result derivable directly by residue
calculus), the recurrence (1.2.1.7) generates `I^{(0,m)}` for all `m ≥ 2` by
repeated differentiation and addition. The resulting closed forms are rational
functions of `η` over polynomial powers of `(1 − η²)`.

**Alignment to SGP4.**
- (a) **Symbol bridge.** `η = a₀'' e₀ ξ` per Phase 0 (Definition 0.6.4.2).
  The integral is parameterized by `η` alone (for `p = 0`); downstream
  formulas use the SGP4 `η` directly.
- (b) **Use in C-coefficients.** The `ψ^{-7} = (1 − η²)^{-7/2}` factor that
  appears throughout the C₂ / C₄ / C₅ formulas (Phase 2-5) traces directly
  to `I^{(0,4)}` evaluated by repeated application of (1.2.1.7). The (1−η²)
  pole of order `m − 1/2` is the **Lane density residue**, named in §4 of
  the main doc.
- (c) **Convergence under `η → 1`.** The orbit-element form `η = a e ξ =
  ae/(a−s)`. As perigee approaches the atmospheric fitting altitude `s`,
  `ξ → ∞` and `η → 1` — so the Lane integrals diverge. This is a
  *physical* singularity (perigee crossing the atmospheric model boundary),
  not an artefact of the integration method.

**Alignment to implementation (deferred per Standard 9-B).** SGP4's
`drag_coefficients.h` (line 122) precomputes `coef1 = (q-s)^4 · ξ^4 · ψ^{-7}`
in `O(1)` arithmetic at orbit initialization. The implementation never
evaluates an integral at runtime; the closed forms derived here are inlined
into the precomputed coefficient groups.

---

## §1.3 Closed Forms for `I^{(0,m)}`

### Theorem 1.3.1 (Base case `I^{(0,1)}`)

```
I^{(0,1)}(η) = 1 / √(1 − η²) = (1 − η²)^{−1/2} .                          (1.3.1.1)
```

**Proof.** Classical result; cite [GR8 §3.616] or derive via residue
calculus.

**Residue derivation.** Substitute `z = e^{if}`, `df = dz/(iz)`,
`cos f = (z + 1/z)/2`. The integrand becomes

```
1/(1 − η · (z+1/z)/2) = 2z / (2z − η z² − η) = −2z / (η z² − 2z + η) .   (1.3.1.2)
```

The denominator `Q(z) := η z² − 2z + η` has roots

```
z_± = (1 ± √(1 − η²)) / η ,                                              (1.3.1.3)
```

with `z_− z_+ = 1` (product of roots equals constant term over leading
coefficient `= η/η = 1`). For `|η| < 1` and `η > 0`, `z_−` is inside the unit
circle (`|z_−| < 1`) and `z_+` outside (`|z_+| > 1`).

Apply the residue theorem:

```
∫₀^{2π} df / (1 − η cos f) = ∮ (−2z/(η Q(z))) · dz/(iz) = (−2/(iη)) ∮ dz/Q(z)
                          = (−2/(iη)) · 2πi · Res_{z=z_−} 1/Q(z)
                          = (−4π/η) · 1/(η (z_− − z_+))
                          = (−4π/η²) · 1/(−2√(1−η²)/η)
                          = (−4π/η²) · (−η)/(2√(1−η²))
                          = 2π / √(1 − η²) .                              (1.3.1.4)
```

Dividing by `2π` gives `I^{(0,1)}(η) = 1/√(1 − η²)`. ∎

### Theorem 1.3.2 (Closed forms `I^{(0,m)}` for `m = 1 … 6`)

**Hypotheses.** As in Definition 1.1.2.

**Conclusion.** Apply Theorem 1.2.1 recurrence starting from Theorem 1.3.1:

```
I^{(0,1)}(η) = (1 − η²)^{−1/2}                                            (1.3.2.1)

I^{(0,2)}(η) = (1 − η²)^{−3/2}                                            (1.3.2.2)

I^{(0,3)}(η) = (2 + η²) / (2 · (1 − η²)^{5/2})                            (1.3.2.3)

I^{(0,4)}(η) = (2 + 3 η²) / (2 · (1 − η²)^{7/2})                          (1.3.2.4)

I^{(0,5)}(η) = (8 + 24 η² + 3 η⁴) / (8 · (1 − η²)^{9/2})                  (1.3.2.5)

I^{(0,6)}(η) = (8 + 40 η² + 15 η⁴) / (8 · (1 − η²)^{11/2})                (1.3.2.6)
```

**Proof.** Mechanical application of (1.2.1.7) to (1.3.1.1).

**Step `m = 1 → 2`.** `I^{(0,2)} = I^{(0,1)} + η · dI^{(0,1)}/dη`.
With `dI^{(0,1)}/dη = η · (1−η²)^{−3/2}`:

```
I^{(0,2)} = (1−η²)^{−1/2} + η² (1−η²)^{−3/2}
         = (1−η²)^{−3/2} · [(1−η²) + η²]
         = (1−η²)^{−3/2} .                                                (1.3.2.7)
```

**Step `m = 2 → 3`.** `I^{(0,3)} = I^{(0,2)} + (η/2) · dI^{(0,2)}/dη`.
With `dI^{(0,2)}/dη = 3 η · (1−η²)^{−5/2}`:

```
I^{(0,3)} = (1−η²)^{−3/2} + (3 η²/2) · (1−η²)^{−5/2}
         = (1−η²)^{−5/2} · [(1−η²) + 3 η²/2]
         = (1−η²)^{−5/2} · [1 + η²/2]
         = (2 + η²) / (2 · (1−η²)^{5/2}) .                                (1.3.2.8)
```

**Step `m = 3 → 4`.** `I^{(0,4)} = I^{(0,3)} + (η/3) · dI^{(0,3)}/dη`.
Compute `dI^{(0,3)}/dη`:

```
dI^{(0,3)}/dη = (1/2) · d/dη [(2+η²) (1−η²)^{−5/2}]
             = (1/2) · [2η (1−η²)^{−5/2} + (2+η²) · 5η (1−η²)^{−7/2}]
             = η (1−η²)^{−7/2} · [(1−η²) + (5/2)(2+η²)]
             = η (1−η²)^{−7/2} · [1 − η² + 5 + 5η²/2]
             = η (1−η²)^{−7/2} · [6 + 3 η²/2]
             = (3η/2) (1−η²)^{−7/2} · (4 + η²) .                          (1.3.2.9)
```

Substituting:

```
I^{(0,4)} = (2 + η²) / (2 (1−η²)^{5/2}) + (η/3) · (3η/2)(4 + η²)/(1−η²)^{7/2}
         = (2 + η²) / (2 (1−η²)^{5/2}) + (η²/2)(4 + η²)/(1−η²)^{7/2}
         = (1/(2 (1−η²)^{7/2})) · [(2 + η²)(1 − η²) + η²(4 + η²)]
         = (1/(2 (1−η²)^{7/2})) · [2 + η² − 2η² − η⁴ + 4η² + η⁴]
         = (1/(2 (1−η²)^{7/2})) · [2 + 3η²]
         = (2 + 3 η²) / (2 (1−η²)^{7/2}) .                                (1.3.2.10)
```

This is (1.3.2.4).

**Steps `m = 4 → 5 → 6`.** Same pattern, expansion is mechanical. The
results (1.3.2.5)-(1.3.2.6) are verified by `verify_phase1.m` (companion).

For brevity the m=5, 6 algebra is omitted from this proof; the verifier
checks both. ∎

### Remark 1.3.3 (Discrepancy with main-doc Eq (4.4) — D-9 PRELIMINARY, pending Phase 2 trace)

The main derivation document `sgp4_near_earth_drag_theoretical_basis.md` Eq (4.4) claims

```
I^{(0,4)}_doc = (2 + η²) / (1 − η²)^{7/2}                                 (1.3.3.1)
```

The closed form derived in Theorem 1.3.2 (1.3.2.4) is

```
I^{(0,4)} = (2 + 3 η²) / (2 · (1 − η²)^{7/2}) .                           (1.3.3.2)
```

**Sanity check at η = 0.** Integrand `1/(1−0·cos f)^4 ≡ 1`, so any
correctly-defined `I^{(0,4)}(0) = (1/2π) · 2π · 1 = 1`. Theorem 1.3.2 gives
`2 / (2·1) = 1`. ✓ The doc value gives `2 / 1 = 2`. ✗

**However, this discrepancy does NOT yet establish that the main doc is
wrong.** Per user directive 2026-05-15:

> "We can't accept results that we aren't sure of. We will need to
> understand how the theorems apply and get the correct result, or it is
> possible that the theorems have an incorrect assumption."

**What is certain:**

1. **Theorem 1.3.2 (1.3.2.4) is the mathematically-correct closed form of
   the integral as defined in (1.0.1) / Definition 1.1.2.** The proof in
   §1.2 is rigorous, and the η=0 sanity check confirms the value.

2. **The SGP4 code in `drag_coefficients.h:146-149` produces a C₂ value that
   matches 7 independent reference SGP4 implementations** (per the
   phantom-bug audit cross-check, commit `03b620a`).

**What is NOT yet certain:**

(a) Whether the SGP4 orbit-average literally computes the integral as
defined in (1.0.1). The Lane fictitious-anomaly substitution
`r − s = (1/ξ)(1 − η cos f†)` introduces O(eη) approximation error (Note
4.3 of the main doc). The actual SGP4 integrand BEFORE the substitution
is `(r − s)^{−4} · X(f)` with `r − s = (aβ² − s − se cos f)/(1 + e cos f)`,
which is NOT of the form `(1 − η cos f)^{−4}`. The Lane substitution
transforms the integrand, but the transformation is the source of an
additional accuracy bound (Phase 11 work).

(b) Whether the main doc's Eq (4.4) is intended as the literal value of the
post-substitution integral, or as something else (e.g., a Hansen
coefficient with different normalization, or the LH79 p.16 transcription
with its own structure).

(c) Whether the SGP4 code's `(8 + 24η² + 3η⁴)` polynomial in C₂ Part B
emerges from `I^{(0,5)}` from Theorem 1.3.2 (1.3.2.5) — preliminary
observation that the polynomial matches the numerator of `I^{(0,5)}` — or
from some other Lane-integral combination, or from a parameterization in
which the doc's Eq (4.4) is consistent.

**Audit-finding status:** D-9 is **PRELIMINARY** (not closed) until Phase
2 traces the explicit connection between the Lane integrals derived here
and the SGP4 C₂ code. The Phase 1 Theorem 1.3.2 closed forms stand as
correct evaluations of the integrals as defined; the question of whether
the SGP4 code computes those integrals (or computes a slightly different
integral that the main doc characterizes by Eq (4.4)) is **open** and is
the first task of Phase 2.

**Required Phase 2 entry task (Phase 2.A):**
1. Start from the SGP4 orbit-averaging step prior to any Lane substitution.
2. Express the integrand explicitly in terms of `f, e, ξ, s, a` (no `η, f†`).
3. Apply the Lane substitution rigorously, identifying its approximation
   order in `e`.
4. Reduce the resulting integrand to a combination of Lane integrals
   `I^{(p,m)}` per Definition 1.1.1.
5. Substitute the Phase 1 Theorem 1.3.2 closed forms.
6. Compare the result to the SGP4 code's `(8 + 24η² + 3η⁴) · ψ^{−7}` form.
7. If they match: D-9 is closed as "main doc Eq (4.4) is a transcription
   error". If they don't match: D-9 is escalated, and Theorem 1.3.2 (or its
   hypotheses) must be re-examined.

**Alignment to SGP4 (Phase 1 standalone — independent of D-9 resolution).**
- (a) Theorem 1.3.2 gives the closed forms of integrals AS DEFINED in
  (1.0.1). These are mathematically valid under `|η| < 1` and
  differentiation-under-the-integral hypotheses (no atmospheric perigee
  case-table boundary in the integrand).
- (b) The connection of these integrals to the SGP4 orbit-average is
  Phase 2 work (specifically Phase 2.A); Phase 1 makes no claim about
  that connection.
- (c) The values `I^{(0,m)}` from Theorem 1.3.2 are *inputs* to Phase 2's
  C₂ derivation; if Phase 2.A determines that SGP4's integral is different
  from Definition 1.1.1, then the Lane-integral inputs to C₂ will be
  re-specified, and Phase 1's Theorem 1.3.2 still stands as the correct
  evaluation of the literally-defined integral.

**Alignment to implementation (deferred per Standard 9-B).** Same as
before — implementation choices are not exercised at the Lane-integral
level; they are downstream of Phase 2 / Phase 6 / etc.

---

## §1.4 The `p ≥ 1` family — scoped, NOT built (Standard-10 decision, Phase 2.C)

Phase 1 was originally scoped (§1.0) to cover `I^{(p,m)}` for `p ∈ {0,1,2,3}`, on the
expectation that the eccentricity-rate coefficients C₄/C₅ — whose `ė` integrands carry
`(1 + e cos f)` factors — would require the `p ≥ 1` family. **Phase 2.C scoping shows
they do not.**

The C₄ secular eccentricity rate orbit-averages the instantaneous `ė` (Phase 0
Theorem 0.3.3) under drag. After the **exact `ė`-bracket collapse**
```
e sin²f + (1 + e cos f)(cos f + cos E) = 2 (e + cos f)        [verified simplify = 0]
```
(`cos E = (e + cos f)/(1 + e cos f)`), the `f`-average kernel is
```
h_C4(e,f) = (1 − η cos f†)^{−4} · (1 + e² + 2e cos f)^{1/2} · (e + cos f) / (1 + e cos f)² .
```
Expanding in `e` at **fixed `η`** (the Phase 2.A exact-density method, §A.3) gives
```
h0 = cos f · (1 − η cos f)^{−4} ,
h1 = (1 + 3η cos f) · sin²f · (1 − η cos f)^{−5} ,
```
both containing **only** `(1 − η cos f)` powers — every `(1 + e cos f)` factor has been
expanded away. They reduce, via the cos-power identity `cos f = (1 − (1 − η cos f))/η`
(1.2.1.3) and the `sin²f` identity, to the `p = 0` family `I^{(0,m)}`, `m ≤ 6`. The
averages land the SGP4 C₄ Keplerian bracket exactly:
```
ψ⁷·⟨h0⟩ = η(2 + ½η²) ,    ψ⁷·⟨h1⟩ = ½ + 2η²        [= drag_coefficients.h:163]
```
(verified `simplify = 0`; the full C₄ assembly is Phase 3, `verify_C4.m`).

**Decision.** Per **Standard 10** — no theoretical result is accepted into the corpus
without a downstream code-match that consumes it — the `p ≥ 1` Lane family is **not
constructed**. The candidates that could plausibly need it all reduce to `p = 0`:
- **C₄, C₅** — the eccentricity-rate coefficients — reduce to `I^{(0,m)}` as shown.
- **D₂–D₄, t2cof…t5cof** — `t`-Taylor coefficients of `a(t)`/`ℓ(t)` — involve **no
  orbit integral at all** (they are polynomials in C₁, C₂, …).

Should **C₃** (Phase 5, the J₃ coupling) require a specific `I^{(p,m)}` with `p ≥ 1`,
it will be derived **just-in-time** at that phase. Phase 1 therefore **closes with the
`p = 0` family** `I^{(0,m)}`, `m = 1 … 6` (Theorem 1.3.2).

## §1.5 Independent verification of `I^{(0,m)}` — `verify_phase1.m`

The `p = 0` closed forms are established by **two independent routes** (corpus
independence requirement, `feedback_faked_verification_auto_fail`):

1. **Analytic recurrence** — Theorem 1.2.1 `I^{(0,m+1)} = I^{(0,m)} + (η/m)
   dI^{(0,m)}/dη` from the residue base case `I^{(0,1)} = (1−η²)^{−1/2}` (§1.3.1).
2. **Numerical quadrature** — high-resolution periodic-trapezoidal evaluation of
   `(1/2π)∫₀^{2π}(1−η cos f)^{−m} df` (spectrally accurate for this analytic periodic
   integrand) vs the closed form, over `η ∈ {0.1, 0.3, 0.5, 0.7, 0.85}`, `m = 1…6`
   (max relative error `≈ 9·10⁻¹⁶`, machine precision).

Both routes plus the `η = 0` sanity (`I^{(0,m)}(0) = 1`) are checked by the companion
verifier **`verify_phase1.m`** (13/13 PASS). The residue-calculus route for `m ≥ 2`
(m-th order pole at `z₋`, §1.3.1 framework) is an available third path, not ground out
in full here because the recurrence + quadrature already confirm every used value,
which is **additionally** retrospectively validated by the Phase 2.A C₂ code-match
(consuming `I^{(0,3)}, I^{(0,4)}, I^{(0,5)}`) under the Standard-10 gate.

## §1.6 Closed status

| Item | Status |
|---|---|
| `I^{(0,m)}`, m = 1…6 (Theorem 1.3.2) | ✅ recurrence + quadrature (`verify_phase1.m` 13/13); C₂ code-match (Phase 2.A) |
| `I^{(p,m)}`, p ≥ 1 (§1.4) | ⛔ not built — no downstream code-match (Standard 10); C₄/C₅ are `p = 0` |
| **D-9** (main-doc Eq 4.4) | ✅ CLOSED (Phase 2.A §A.4.2: `I^{(0,4)}` code-matched) |
| **D-3** (master orbit-average identity) | ◐ resolved for C₂ (Phase 2.A §A.2, β³ Jacobian); finalized when C₄/C₅ land (Phase 3/4) |

Phase 1 is **closed** at the `p = 0` family — the complete set the SGP4 drag
code-matches consume.

---

## §1.7 Pre-Closure Summary (this commit)

**Theorems established (with full multi-pass-aligned proofs):**

| Theorem | Subject |
|---|---|
| **1.2.1** | Lane integral recurrence `I^{(0,m+1)} = I^{(0,m)} + (η/m) · dI^{(0,m)}/dη` |
| **1.3.1** | Base case `I^{(0,1)} = (1−η²)^{−1/2}` (residue derivation provided) |
| **1.3.2** | Closed forms `I^{(0,m)}` for `m = 1 … 6` |

**Audit findings closed by this commit:**

| ID | Subject | Closure |
|---|---|---|
| **D-9** | Main-doc Eq (4.4) incorrect | Theorem 1.3.2 (1.3.2.4) supersedes |

**Closed status (updated Phase 2.C, 2026-06-01):** §1.0-§1.6 complete.
§1.4 records the Standard-10 decision that the `p ≥ 1` family is **not built**
(C₄/C₅ scoping shows they reduce to `p = 0`); §1.5 documents the two-route
independent verification; §1.6 is the closed-status summary. The companion
verifier **`verify_phase1.m` is written and PASSes 13/13** (recurrence + η=0
sanity + numerical quadrature). Phase 1 is **closed at the `p = 0` family**
`I^{(0,m)}`, m = 1…6 — the complete set the SGP4 drag code-matches consume.
D-3 is resolved for C₂ (Phase 2.A §A.2) and finalized when C₄/C₅ land.
