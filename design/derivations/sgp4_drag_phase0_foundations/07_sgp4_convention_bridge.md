## §0.6 SGP4-Convention Bridge

This section pins down the symbol bridge from Phase 0's generic Keplerian
notation to SGP4's specific identifiers. Each downstream Phase will invoke this
bridge explicitly so that no symbol means two different things by silent
context.

### Definition 0.6.1 (SGP4 time scale)

SGP4 measures time in **minutes from epoch**:

```
τ := (t − t₀)    [min]                                                (0.6.1.1)
```

where `t₀` is the TLE epoch time. All rates in the SGP4 implementation
(`drag_coefficients.h`, `secular_update.h`) carry units of `[*]/min`.

**Alignment to Phase 0.** Phase 0 derivations are in SI seconds for clarity
(`n = √(μ/a³)` in `rad/s` when `μ` is in `m³/s²` and `a` in `m`). The conversion
to minutes is multiplicative: in SGP4 units `μ = 398600.4418 km³/s² × (60 s/min)²`,
or equivalently `n_min = 60 · n_sec`. The algebraic form of every Phase 0
theorem (e.g. Theorem 0.3.2 `ȧ = (2/(n β))·[...]`) is **scale-invariant**: if
`n` is in `rad/min`, then `ȧ` is in `ER/min` provided `R, T` are in `ER/min²`.
No theorem statement changes.

**Alignment to implementation (deferred per Standard 9-B).** Some numerical
constants (`XKMPER`, `KE`, etc.) are pre-computed in minute / Earth-radius
units. See `constants_provider.h` for the canonical values.

### Definition 0.6.2 (Brouwer-recovered orbit elements `a₀'', n₀''`)

The TLE-stored mean motion `n_kozai := 2π / period` is the **Kozai-mean** motion.
SGP4's initialization (per [SR3] §2 and our `element_recovery.h`) transforms
this to the **Brouwer-recovered** mean motion `n₀''` via the
J₂-secular-pre-correction:

```
a₀ := (μ / n_kozai²)^{1/3} ,                                          (0.6.2.1)
δ₀ := (3 k₂ / (a₀² β₀³)) · (3 θ² − 1) / 2 ,                            (0.6.2.2)
a₀'' := a₀ · (1 − δ₀/3 − δ₀² − …)    (Brouwer iterative recovery) ,    (0.6.2.3)
n₀'' := √(μ / (a₀'')³) .                                              (0.6.2.4)
```

(The iterative form (0.6.2.3) is in `element_recovery.h:138-150`; the
audit-tagged R08 covers its `O(δ₁⁴)` tail accuracy bound.)

**Alignment to Phase 0.** Throughout this corpus, when a Phase 0 theorem
writes `n` or `a`, it is understood that the SGP4 substitution is
`n → n₀''`, `a → a₀''`. The Phase 0 derivations make no statement about
the recovery transformation itself — that is the BH61 cleanroom result
(see Phase 0a, below).

**Alignment to implementation.** `element_recovery.h` is the SGP4 init step
that produces `(a₀'', n₀'')`. Phase 0 theorems can be cited as the source for
formulas that consume `(a₀'', n₀'')` downstream (e.g. Phase 2's `C₂` formula).

### Definition 0.6.3 (SGP4 atmospheric parameters)

SGP4's atmospheric model is the **Lane (1965) power-law density** with the
following parameters:

```
ρ(r) := ρ₀ · ((q_0 − s) / (r − s))^τ     for r > s ,                  (0.6.3.1)
τ := 4 ,                                                              (0.6.3.2)
q_0 := 1 + 120 km / XKMPER     (ER) ,                                 (0.6.3.3)
s := 1 + 78 km / XKMPER        (ER, nominal — adjusted for low perigee) , (0.6.3.4)
ρ_0 := 2.461 × 10⁻⁵ · XKMPER⁻¹ · (rescaled density unit) .           (0.6.3.5)
```

The `s` parameter has three regimes depending on perigee altitude `h_p`
(per [SR3] line 467-469):

- **Regime I** `h_p > 156 km`: nominal `s = 1 + 78/XKMPER`.
- **Regime II** `98 < h_p ≤ 156 km`: `s* = 1 + (h_p − 78)/XKMPER`,
  `(q − s*)^4 = [(q_0 − s)_{nom} + (s − s*)]^4`.
- **Regime III** `h_p ≤ 98 km`: `s* = 1 + 20/XKMPER`; analogous adjustment.

The case boundaries maintain `C⁰` continuity of `ρ(perigee)` but not `C¹`.

**Alignment to Phase 0.** The density model is a **postulate** of the SGP4
problem, not a derived result. Theorem 0.4.1 (drag specialization) consumes it
as a given. Phase 1 (Lane integrals) evaluates orbit averages of `ρ(r)·X(f)`
under this specific functional form.

**Alignment to implementation.** `density_model.h` carries the case-table for
the perigee regimes. The `C⁰`-only discontinuity at the 156 km boundary is
catalogued as `A-D4` error source.

### Definition 0.6.4 (Lane intermediate quantities)

The orbit-averaging of `ρ(r(f))·X(f)` for `τ = 4` produces the natural
combinations

```
ξ := 1 / (a₀'' − s)              (units: ER⁻¹) ,                       (0.6.4.1)
η := a₀'' · e₀ · ξ                (dimensionless) ,                    (0.6.4.2)
ψ² := |1 − η²|                   (dimensionless) ,                     (0.6.4.3)
β₀ := √(1 − e₀²) , β₀² := 1 − e₀² (dimensionless) .                  (0.6.4.4)
```

The `(q_0 − s)^4 · ξ^4` and `ψ^{−7}` factors that appear throughout the
C-coefficient formulas trace to the residue-calculus evaluation of the Lane
integrals at `τ = 4` (Phase 1).

### Remark 0.6.5 (β and ψ are NOT the same)

A frequent source of confusion in the Lane-Hoots literature (and a contributing
factor to audit finding **D-7**): `β² = 1 − e²` (the Keplerian factor from
`p = a β²`) and `ψ² = |1 − η²|` (the Lane density-residue factor from
orbit-averaging `(r − s)^{−4}`) are **distinct**.

- `β` depends only on `e`: `β = √(1 − e²)`.
- `ψ` depends on `η = a e ξ` (eccentricity × scale-height ratio): `ψ = √(1 − η²)`.

For typical LEO with `e₀ = 0.01` and `a₀'' ≈ 1 + 800/XKMPER`,
`s ≈ 1 + 78/XKMPER`, `ξ ≈ XKMPER / 722 km ≈ 8.83 ER⁻¹`, `η ≈ a₀ e₀ ξ ≈ 1.01 · 0.01 · 8.83 ≈ 0.089`.
So `β ≈ 0.99995`, `ψ ≈ 0.996`. Numerically close but **algebraically distinct**.

Claims like "absorb β^{−3} into C₁ via the ψ^{−7} factor" (which appeared in the
main doc's §13 proof and triggered audit finding **D-7**) are **wrong**: `β` and
`ψ` are not interconvertible. Phase 9 (xnodcf re-derivation) must avoid this
identification.

**Alignment to Phase 0.** All Phase 0 theorem statements use only `β`, never
`ψ`. The `ψ` factor enters in Phase 1 (Lane integrals) as a separate
residue-calculus output, not as a Keplerian-geometry quantity.

### Remark 0.6.6 (Non-rotating-atmosphere assumption)

Phase 0 §0.4 used the Lane non-rotating-atmosphere assumption (`v_rel = v`,
Assumption 1.3 of the main document). This is a SGP4 model choice with
a `~10⁻¹` accuracy cost catalogued as `A-D1`. Phase 0 theorems are valid
under any choice of `v_rel`; the drag specialization (Theorem 0.4.1) is what
fixes `v_rel = v`.

For a generalized preset that includes Earth rotation, replace (0.4.1) with
`F_drag = −B* ρ |v_rel| v_rel` where `v_rel = v − ω_⊕ × r`, and the
specialization in Theorem 0.4.1 picks up additional `ω_⊕ × r` terms in
`(R_drag, T_drag, N_drag)`. The `N_drag = 0` conclusion no longer holds.
This generalization is **out of scope for the current SGP4 corpus** but is
flagged here for future extension.

### Sub-phase 0a (BH61 cleanroom dispatch — pending)

The Phase 3 (C₃, J₃ coupling) and Phase 9 (xnodcf, J₂-RAAN coupling)
derivations require, respectively, the **Brouwer J₃ long-period generator**
result `δr_{J_3}` (main doc Eq. 7.1) and the **Brouwer J₂ secular RAAN rate**
`Ω̇_{J_2} = -3 k₂ n cos i / (p² β)` (main doc Eq. 13.1 setup).

These two results live in the BH61 cleanroom at
`sgp4_references/vallado_celestrak/documentation/SGP4/Brouwer_Hori 1961/derivation/`.
Per binding rule `feedback_consolidator_no_BH61_reads.md`, the main session
cannot directly read those files. The Phase 0a sub-task dispatches a sub-agent
to extract the two results and return:

- **Cite-key**: which BH61 cleanroom chapter/theorem produces the result.
- **Closed form**: the result expression in our notation `(a, e, i, p, …)` plus
  the chain of named theorems / lemmas the BH61 cleanroom proof invokes.
- **Verifier handle**: which verifier file in the BH61 cleanroom mechanically
  checks the result.
- **No quoted equations** from the cleanroom (per binding rule).

The Phase 0a dispatch is **deferred until Phase 3 / Phase 9 begin** to avoid
loading unnecessary context. Phase 0a is listed in the audit log's Phase 1.5
Re-derivation Queue and the todo list.

---

