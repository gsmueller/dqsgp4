# Phase 0 — Kepler Geometry and Gauss-Lagrange Variational Equations

## §0.0 Scope and Pedagogy

This is the **Phase 0** foundational document of the symbolic re-derivation of the
SGP4 near-earth drag coefficients, executed per the multi-pass audit at
`design/audit/2026_05_15_sgp4_drag_derivation_full_audit.md`.

### Why this document exists

The pre-existing `design/derivations/sgp4_near_earth_drag_theoretical_basis.md` and
its companion `design/derivations/deprecated/020_c2_drag_integral_derivation.md` are
**OCR-transcriptions** of [LH79] / [SR3] with cross-checks to code, **not** symbolic
re-derivations from first-principles theoretical sources. The 2026-05-15 audit
identified **8 derivation findings (D-1 through D-8)** in the main document, plus
the explicit user verdict that the entire `deprecated/` folder is
*pedagogy-deprecated*.

Phase 0 begins the replacement corpus. It establishes the **foundational primitives**
that every subsequent phase (Lane integrals, C-coefficients, D-coefficients, t-cofs,
xnodcf, xlcof / aycof, omgcof / xmcof, Lane f†) consumes. All Phase 0 content
is derivable from elementary mechanics + Kepler's laws — no OCR sources, no BH61
sub-agent dispatch required.

### Pedagogy standard

Per `feedback_textbook_theorem_rigor.md`, `feedback_solver_derives_from_theory.md`,
`feedback_derivation_first.md`, `feedback_no_TD_labels.md`,
`feedback_born_digital_latex.md`, `feedback_no_bias_toward_reference.md`,
`feedback_no_convention_for_incorrectness.md`:

1. **Theorem / Definition / Lemma / Proposition** are the only structural blocks. No
   "(T)" / "(D)" / "(U)" labels.
2. **Hypotheses are explicit and labeled** `(H₁), (H₂), …`.
3. **Every step in every proof is numbered and labeled** with a bold descriptive
   step label, e.g. `**Step 1 (Kepler 2nd law).**`. No steps may be skipped.
4. **Equation numbers** use the form `(N.X.E.k)` — section N, theorem X, equation k.
5. **Citations are bold-markup** referring to named results:
   `**By Theorem 0.2.3**`, `**By Definition 0.2.6**`.
   Citations to external sources may be made **only** to born-digital references:
   Wikipedia / ProofWiki / nLab for canonical theorems; Battin (1999) / Roy (2005) /
   Vallado-LaTeX (2013) / Lara (2021) / Sneeuw (2022) for astrodynamics; never to OCR'd
   Lane-Hoots / Brouwer / SR3 primary papers.
6. **Every theorem closes** with `∎`.
7. **A mechanical verifier** (Octave / SymPy) accompanies each phase and reproduces
   each closed-form result. Verifier returns 0 on PASS.
8. **Every approximation** is named and bounded: dropped terms must be identified
   with an explicit error-order estimate added to the precision / accuracy / uncertainty
   catalog.
9. **Alignment-to-SGP4 check after each theorem.** Born-digital sources solve the
   *general* problem; SGP4 solves a *specific* problem with its own conventions
   (osculating elements vs constants; Brouwer-recovered `n₀'' / a₀''` vs textbook
   `n / a`; minute-scale time; Lane density model; AFGP4 → SGP4 simplifications;
   the TEME frame). Each theorem in this corpus is followed by an explicit
   `**Alignment to SGP4.**` remark that pins down (a) which SGP4-side identifier the
   theorem-side symbol maps to, (b) what hidden assumption the SGP4 application
   carries, and (c) whether any specialization or approximation has been applied
   beyond what the general theorem states.

### The alignment principle

> "Sources may have features that appear correct, but are misaligned on scrutiny."
> — user directive, 2026-05-15.

Concrete alignment risks for the Phase-0 layer:

- **Osculating vs constant elements.** Textbook Kepler theory treats `(a, e, i, ω, Ω, M)`
  as constants of motion of the unperturbed two-body problem. SGP4 treats them as
  *osculating* elements under the **variation-of-parameters (VOP)** formulation —
  they vary slowly, and the Keplerian relations are applied instantaneously. Every
  use of `r(f)`, `h = √(μp)`, `dM/df` in this corpus assumes VOP.
- **Brouwer-recovered vs textbook elements.** SGP4 stores Brouwer-recovered
  `n₀'' / a₀''` (J₂-secular-pre-corrected); these are NOT the textbook
  instantaneous-Keplerian `n / a`. They are related by the Brouwer recovery
  transformation (see element_recovery.h / `deprecated/013` content). In the
  drag derivations, every `n` and `a` is implicitly `n₀''` and `a₀''`.
- **Time scale.** Textbook Kepler uses SI seconds; SGP4 uses minutes from epoch.
  All rates here are time-scale-invariant in form but the numerical implementation
  must be consistent.
- **Frame.** Textbook orbital relations are frame-agnostic (orbital-plane geometry);
  SGP4 outputs in the TEME frame. The drag derivations all live within the orbital
  plane and are therefore frame-independent at the theorem level. Frame matters only
  at the position/velocity-output level (handled in `state_from_elements.h`).
- **Lane density model.** The atmospheric density profile in §2 (Phase 1+) is a
  SGP4-specific model choice — `ρ(r) = ρ₀·((q₀-s)/(r-s))^τ` with `τ = 4`. No
  textbook canonicalises this; alignment-to-SGP4 means we postulate this specific
  form and propagate consequences.
- **AFGP4 → SGP4 simplifications.** Lane-Hoots 1979 retains higher-order eccentricity
  terms that SGP4 drops (the `O(e²)` Part A terms; the `-5eη(4+3η²)` Part B term).
  Each dropped term must be flagged with an explicit accuracy bound; the SGP4
  formula is then identified as an approximation of AFGP4 with named, bounded
  residuals.

When any **Alignment to SGP4** remark cannot be discharged (e.g., the textbook
result is genuinely different from what SGP4 uses), the theorem must be
re-specialized within this corpus — not silently adopted under a hope that
"close enough" matches the problem.

### The implementation-alignment principle (Standard 9)

> "Even with correct theory, there is a risk. For example, if the theory derives
> Taylor, but we need a faster convergence, then the coefficients may need
> specialized alterations. Worry about one problem at a time though."
> — user directive, 2026-05-15.

Even when a derivation is **theoretically correct** (Standard 1) AND **aligned to
SGP4's specific problem** (Standard 9 part A), the **implementation** may use a
different concrete form than the theoretical canonical one:

- **Series re-summation.** The theory may yield a Taylor series; the implementation
  may use a Padé approximant, Chebyshev rearrangement, or asymptotic re-expansion
  for faster convergence on a specific argument range. The literal Taylor
  coefficients are not what the code computes; the coefficient transformation must
  itself be derived, error-bounded, and documented.
- **Range reduction.** The theory may state a result for `θ ∈ [0, 2π)`; the
  implementation may reduce to `θ ∈ [-π, π]` or `[0, π/2]` with sign-tracking, and
  the resulting half-/quarter-angle formulae change the coefficient layout.
- **Pre-computed constants.** Per `feedback_compute_once.md` and
  `feedback_lambda_precompute.md`, the implementation may pre-compute coefficient
  groups at orbit-initialization time that the theory presents as nested formulas;
  alignment requires checking that the precomputed groups are algebraically
  equivalent to the theoretical form.
- **Convergence acceleration.** A series whose theoretical convergence is
  geometric in `e²` may need an Euler / Shanks / Wynn ε-acceleration to be
  efficient for `e ≳ 0.5`; the accelerated coefficients are functions of the
  raw theoretical coefficients but are *not* equal to them.
- **Iterative solvers.** Kepler's equation has Newton-style and Halley-style
  iterators; their convergence orders (quadratic vs cubic) come from different
  series expansions of the same underlying inversion. The "coefficients" in the
  implementation are the Newton increment vs Halley increment, not anything
  appearing in the canonical statement of Kepler's equation.

When an implementation diverges from the canonical theoretical form, a separate
remark `**Alignment to implementation.**` flags the divergence in the relevant
Phase, with: (a) the implementation form actually computed, (b) the
coefficient-transformation derivation (which is itself subject to Standards 1-2),
and (c) the convergence / accuracy consequence.

**Status for Phase 0:** none of the Phase 0 theorems (Kepler geometry, Gauss-Lagrange
VE in symbolic form) make implementation choices; they are pure theory. Standard 9
is recorded here for completeness and will be exercised in later Phases (notably
Phase 1's Lane-integral residue evaluation, Phase 6-8's Taylor-vs-binomial-series
choice for D-coefficients and t-cofs, and any modified-Kepler iterator analysis
that downstream surfaces). Per the user directive: "Worry about one problem at a
time."

### The code-matching principle (Standard 10, standing rule 2026-05-15)

> *"The solver shall not use theoretical results that can not be matched to the
> code. Numerical testing at points is not acceptable."*
> — user standing directive, 2026-05-15.

A standing rule (formally recorded in
`design/audit/2026_05_15_sgp4_drag_derivation_full_audit.md` STANDING BINDING
RULE) that all theoretical results consumed by the SGP4 corpus must be
symbolically matched to a corresponding code expression by a chain of named
theorem applications and definition-based substitutions, with step-by-step
equality demonstrable in the cleanroom proof format. Numerical agreement at
sample points (e.g. η = 0, e = 0.01) is **insufficient** as proof of equality
— it can serve only as a sanity check.

For Phase 0 specifically:
- The Gauss-Lagrange VE in Theorems 0.3.2-0.3.6 are **intermediate**
  theoretical primitives. They have no direct code match in
  `drag_coefficients.h` or other SGP4 source files. Their validity is not
  established standalone; it will be established **retrospectively** when
  Phase 2 / Phase 4 / etc. symbolically derive a code-matched C-coefficient
  expression that consumes them.
- If a downstream Phase fails to symbolically match the code, the failure
  identifies an error in some step of the chain — possibly in Phase 0, Phase 1,
  the substitution between them, or the orbit-averaging step.
- The Phase 0 mechanical verifier `verify_phase0.m` (27/27 PASS) confirms
  internal algebraic consistency, but does NOT establish that any Phase 0
  theorem is code-matched. Code-matching for Phase 0 theorems is intrinsically
  retrospective.

**Pedagogy implication:** every Phase ≥ 2 document must end with a
**Code-match witness** section that explicitly shows symbolic equality
between the derived closed form and the corresponding code expression, with
references to the source file and line numbers. Phase 0 / Phase 1 documents do
not have such a section because they produce intermediate results, not
code-matched closed forms.

### Born-digital source list for Phase 0

| Cite | Reference | URL or canonical source |
|---|---|---|
| **[WIKI-OE]** | Wikipedia: *Orbit equation* | https://en.wikipedia.org/wiki/Orbit_equation |
| **[WIKI-KL]** | Wikipedia: *Kepler's laws of planetary motion* | https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion |
| **[WIKI-KE]** | Wikipedia: *Kepler's equation* | https://en.wikipedia.org/wiki/Kepler%27s_equation |
| **[WIKI-GVE]** | Wikipedia: *Orbital perturbation analysis (variation of parameters method) / Gauss planetary equations* | https://en.wikipedia.org/wiki/Orbital_perturbation_analysis |
| **[WIKI-OST]** | Wikipedia: *Orbital state vectors* | https://en.wikipedia.org/wiki/Orbital_state_vectors |
| **[BATT99]** | Battin, R.H. (1999), *An Introduction to the Mathematics and Methods of Astrodynamics, Revised Edition*. AIAA Education Series. (Born-digital edition; chapters 3, 9, 10.) |
| **[ROY05]** | Roy, A.E. (2005), *Orbital Motion*, 4th ed. Institute of Physics. (Born-digital edition; chapter 6.) |
| **[VAL13L]** | Vallado, D.A. (2013), *Fundamentals of Astrodynamics and Applications*, 4th ed. (LaTeX-typeset edition; chapter 9.) |

Each named theorem below cites the canonical born-digital reference for the
classical result it specializes; the proof in this document is then independently
written without quoting the OCR sources.

---

