## §0.5 Energy Identity and Orbit-Averaged Drag Rates

This section connects the in-plane Gauss VE for `ȧ` to the energy form (and to
the Delaunay momentum `L`), and sets up the orbit-averaging framework that
Phase 1 (Lane integrals) will evaluate. The actual evaluation of the orbit
integrals is **deferred to Phase 1**; here we only define the orbit-averaged
quantities and verify the elementary algebraic identities.

### Theorem 0.5.1 (Energy-↔-`a` identity for the rate)

**Hypotheses.** As in Theorem 0.3.2 (perturbed two-body, VOP).

**Conclusion.**

```
dE/dt = (μ / (2 a²)) · ȧ ,                                            (0.5.1.1)
```

where `E = −μ/(2a)` is the specific orbital energy.

**Proof.** Direct differentiation of `E = −μ/(2a)`:

```
dE/dt = −μ · d(1/a)/dt · (1/2)
      = −μ · (−1/a²) · ȧ · (1/2)
      = (μ / (2 a²)) · ȧ .                                           (0.5.1.2)
```

∎

**Remark 0.5.1.3 (Compatibility with Theorem 0.3.2).** Combining (0.5.1.1) with
(0.3.2.4) `dE/dt = v · F` yields

```
(μ / (2 a²)) · ȧ = v · F   ⇒   ȧ = (2 a² / μ) · (v · F) .              (0.5.1.4)
```

Re-deriving via (0.3.2.5) `v · F = ṙ R + r ḟ T`, (0.3.2.10) `ṙ = (n a e sin f)/β`,
and (0.3.2.11) `r ḟ = (n a / β)(1 + e cos f)` reproduces (0.3.2.15). This is a
**consistency check** between the two derivation paths of Theorem 0.3.2 (the
direct energy derivation and the Gauss VE algebraic form).

### Definition 0.5.2 (Delaunay momentum and its rate)

The **Delaunay momentum** conjugate to the mean anomaly `M` is

```
L := √(μ a)    (units of [L]² / [T] in SI; ER² / min in SGP4 units) .  (0.5.2.1)
```

Its instantaneous rate is

```
dL/dt = (μ / (2 L)) · (ȧ / a) · a = (μ / (2 √(μa))) · ȧ
       = (1/2) · √(μ/a) · ȧ = (n a / 2) · ȧ                            (0.5.2.2)
```

(using `n = √(μ/a³)` and so `√(μ/a) = n a`).

**Alignment to SGP4.** SGP4's Lane / Lane-Hoots derivation works in `L` because
`L` is conserved under the unperturbed two-body flow (whereas `ȧ` has a
trivially-zero unperturbed rate). The C₁ definition in [SR3]'s Eq. 4 / our main
doc §6.1 is in terms of `L̇/L`:

```
⟨dL/dt⟩_drag / L = ?                                                   (0.5.2.3)
```

This is the canonical "L-rate" form that Phase 2 will close. By (0.5.2.2),
`⟨dL/dt⟩_drag = (n a / 2) · ⟨ȧ⟩_drag`, so the two forms are equivalent up to
the trivial factor `n a / 2`.

### Theorem 0.5.3 (Setup of the orbit-averaged drag rate)

**Hypotheses.** As in Theorem 0.4.2 (drag specialization, instantaneous rates
(0.4.2.1)).

**Conclusion.** The orbit-averaged secular rate of `a` driven by drag is

```
⟨ȧ⟩_drag = (1 / (2π)) · ∫_0^{2π} ȧ_drag dM .                          (0.5.3.1)
```

Converting `dM ↔ df` via Theorem 0.2.5 (`dM = β³/(1+e cos f)² df`):

```
⟨ȧ⟩_drag = (1 / (2π)) · ∫_0^{2π} ȧ_drag · (β³ / (1 + e cos f)²) df ,    (0.5.3.2)
```

with `ȧ_drag` from (0.4.2.1) — Phase 0-rev1, corrected for D-10 with the
Newtonian quadratic drag form — substituted as a function of `f`:

```
ȧ_drag(f) = −(2 n B* ρ(r(f)) · a² / β³) · (1 + e² + 2 e cos f)^{3/2} . (0.5.3.3)
```

**Proof.** Direct application of Definition 0.2.6 (orbit-average) to the
instantaneous drag rate (0.4.2.1). The Jacobian `β³/(1+e cos f)²` is from
Theorem 0.2.5. ∎

**Remark 0.5.3.4 (Why this is left as an integral).** The integrand of (0.5.3.2)
depends on `ρ(r(f))` and on the `(1 + e² + 2 e cos f)^{3/2}` factor from
(0.5.3.3). Under the Lane density model `ρ(r) = ρ₀ ((q − s)/(r − s))^τ` with
`τ = 4`, and after the Lane f† substitution `r − s = (1/ξ)(1 − η cos f†)`
(approximate at first order in `e`; full treatment is Phase 11 / Phase 2.A §A.3),
the integrand becomes a product of:

- `(1 − η cos f†)^{−4}` (Lane density power-law factor),
- `(1 + e² + 2 e cos f)^{3/2}` (quadratic-drag fingerprint, expanded as a
  polynomial in `(e, cos f)`),
- `β³/(1 + e cos f)²` (Jacobian).

The **residue-calculus evaluation** of the pure `(1 − η cos f†)^{−m}` part is
the content of **Phase 1** (Theorem 1.3.2 gives `I^{(0,m)}` closed forms for
`m = 1..6`). The further expansion of `(1 + e² + 2 e cos f)^{3/2}` as a
polynomial in `(e, η)`, and the reduction of the full integrand to a specific
combination of Lane integrals `I^{(p, m)}`, is the content of **Phase 2.A**
(the symbolic trace from Phase 0-rev1 + Phase 1 to the SGP4 code's `C₂`
expression at `drag_coefficients.h:146-149`). Phase 0-rev1 has set up the
integral; Phase 1 evaluates the `(1 − η cos f†)^{−m}` Lane integrals; Phase
2.A traces the chain to the SGP4 code expression.

**Alignment to SGP4.**
- (a) **The integrand structure (Phase 0-rev1).** The `(1 + e² + 2 e cos f)^{3/2}`
  factor in (0.5.3.3) is the **trigonometric content of Newtonian quadratic
  drag**. It combines `|v|^{1/2}` (the speed scalar's structure, since `|v|² =
  (n a/β)²·(1+e²+2e cos f)`, (0.4.1.11)) with the bracket factor `^1` produced
  by the Gauss VE for `ȧ` (Step 1 of Theorem 0.4.2). The `R_drag · ṙ` and
  `T_drag · r ḟ` contributions to (0.4.2.1) combined into the `(1 + e² + 2 e cos f)`
  bracket via (0.4.1.10); the `|v|`-substitution then promoted this to the
  `^{3/2}` form. Phase 2.A's symbolic trace projects this against `ρ(r(f))`
  via the Lane substitution + Lane integrals to produce the C₂ closed form.
- (b) **Why both R_drag and T_drag are needed.** As flagged in Remark 0.4.1.6
  and the audit-log's D-1 finding, the radial drag is `O(e)` smaller than the
  transverse but is **non-zero**. Phase 2.A's symbolic trace must produce
  contributions from both. Dropping `R_drag` would yield a different
  `(1 + e cos f)²` bracket (instead of `1 + e² + 2 e cos f`) and lose the
  `4 e η` coefficient in the eventual C₂ Part-A polynomial.
- (c) **AFGP4 vs SGP4 simplification step.** The AFGP4 result retains
  additional `O(e²)` terms (the `(3/4) e² + 3 e² η²` Part-A residuals); SGP4
  drops them. The Phase 0-rev1 result (0.5.3.2)/(0.5.3.3) is the AFGP4
  integral form before any truncation. Phase 2.A will perform the AFGP4 → SGP4
  simplification with explicit accuracy bounds for the dropped terms
  (resolving the audit's **A-D6 / A-D7** error sources).
- (d) **D-10 closure connection.** Prior to Phase 0-rev1, the integrand
  (0.5.3.3) had a `^1` exponent rather than `^{3/2}`, because the Phase 0
  §0.4 drag specialization missed the `|v|` factor. Restoring `|v|` in
  Theorems 0.4.1 / 0.4.2 (Phase 0-rev1, Blocks A1-A2) propagates here as the
  `^{3/2}` correction. Phase 2.A's symbolic trace against the SGP4 code
  expression is what retrospectively validates the correction.

**Alignment to implementation (deferred per Standard 9-B).** None at Phase 0.
The Phase 2.A closed form will involve precomputed group `coef1 = (q-s)^4 ξ^4 ψ^{-7}`
that the implementation computes once at orbit-initialization
(`drag_coefficients.h:121-122`); the algebraic equivalence is a Phase 2.A
concern under the Standard-10 code-match witness.

---

