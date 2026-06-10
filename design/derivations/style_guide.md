# Style Guide: Mathematical Derivations

This document is the single authoritative reference for all chapter files under `design/derivations/`. It governs structure, formatting, mathematical typesetting, and quality standards. It is separate from `design/style_guide.md`, which governs C++ code conventions.

---

## 1. Chapter File Structure

Every chapter file follows this layout, top to bottom:

```
# Chapter N: Title

**Part X: Part Title**

## Notation

[Notation table]

---

## §N.1 Introduction

[Opening, roadmap, forward-reference table, maturity, model error]

---

## §N.2 First Body Section

[Definitions, Theorems, Proofs, Examples, Remarks]

---

...

---

## §N.k Summary and Downstream Usage

[Usage pattern, downstream reference table]

---

## Error Notes

[All [M.N.x], [P.N.x], [A.N.x] notes collected here]
```

Use `---` (horizontal rule) between top-level sections for visual separation.

---

## 2. The Introduction Section (§N.1)

Every chapter begins with §N.1 containing five elements in this order.

### 2a. Opening Paragraph

One paragraph stating the chapter's purpose and what it establishes. Present tense, third person.

### 2b. Roadmap

A brief mapping of the section structure:

```markdown
The chapter proceeds from [topic] (§N.2) through [topic] (§N.3), ...,
and concludes with [topic] (§N.k).
```

### 2c. Forward-Reference Table

A table showing which later chapters depend on this chapter's results:

```markdown
**Forward references.**

| Section | Feeds | Role |
|---------|-------|------|
| §14.3–4 ($q_0$, $q_0'$, $\arctan(e')/e'$) | Ch 5 | Series evaluation with error bounds |
| §14.5 ($J_2 \to e^2$ iteration) | Ch 1, §1.6 | Application of contraction mapping theory |
| §14.6 ($J_{2n}$) | Ch 7 | Zonal perturbation theory |
```

### 2d. Maturity Declaration

Required for every chapter. Identifies which sections are at which phase of the Draft/Develop/Build cycle:

```markdown
**Maturity.** Sections §§14.3, 14.4, 14.8 are at Build level (complete proofs).
All other sections are at Develop level (proof approaches and numerical examples,
but no full proofs).
```

If the entire chapter is at one level: `**Maturity.** All sections are at Build level.`

### 2e. Model Error Context

For chapters that introduce physical models (Ch 8, 13, 14, 21, etc.), include a paragraph on model error — what the model assumes and where reality departs:

```markdown
**Model error.** The equipotential ellipsoid is a mathematical model, not a
description of the real Earth. Every derived constant inherits an accuracy error
(model error) from the gap between the model assumptions — hydrostatic equilibrium,
exact ellipsoidal shape, time-invariant parameters — and the physical Earth, which
deviates from equilibrium, has a geoid that differs from the ellipsoid by up to
~100 m, and whose parameters (notably $J_2$) change secularly. The model error for
each parameter class is characterized in the Error Notes (§14.10) and tracked
through the tier classification.
```

This element is omitted for pure-mathematics chapters (Ch 1, 4, 5, 6, 7).

---

## 3. Notation Tables

Placed immediately after the part marker, before §N.1. Three columns: **Symbol | Meaning | Introduced**.

### Grouping

For chapters with 15+ symbols, group by semantic category using bold subheader rows:

```markdown
## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| | **Geometric constants** | |
| $a$ | Semi-major axis of the reference ellipsoid | §14.2 |
| $b$ | Semi-minor axis | §14.3 |
| | **Physical constants** | |
| $GM$ | Geocentric gravitational constant | §14.2 |
| $\omega$ | Angular velocity of rotation | §14.2 |
| | **Derived quantities** | |
| $e$ | First eccentricity | §14.3 |
```

For chapters with fewer than 15 symbols, a flat (ungrouped) table is acceptable.

### Cross-references

External symbols include their source: `Ch 1, Def. 1.2.1` or `Ch 14, §14.3`.

### Matched-pair markers

Append **[MP]** to any symbol whose value is constrained by the matched-pair principle (Ch 3):

```markdown
| $\mu$ | Gravitational parameter $= GM$ **[MP]** | §8.2 |
```

---

## 4. Formal Items

### Format

All formal items use the same structure:

```markdown
**ItemType N.k.j** (Name)**.** *Statement in italics.*
```

The double-bold terminator `**.**` separates the label from the italic statement. The name in parentheses is plain text.

### Verbatim example (Definition)

```markdown
**Definition 1.2.1** (Measurement error)**.** *Let $\tilde{v}$ denote the true
value of a physical quantity and $v_{\mathrm{meas}}$ the value obtained by
measurement. The measurement error bound $\sigma_m \geq 0$ satisfies*

$$|v_{\mathrm{meas}} - \tilde{v}| \leq \sigma_m. \tag{1.1}$$

*Measurement error is a property of the physical input. It cannot be reduced
by computation; it can only be reduced by better measurement.*
```

### Verbatim example (Theorem with proof)

```markdown
**Theorem 1.4.1** (Derivative bound for scalar functions)**.** *Let
$f: \mathbb{F} \to \mathbb{F}$ be continuously differentiable on a region
containing the closed disk $\{z : |z - x| \leq \delta\}$. Then*

$$|f(x + \epsilon) - f(x)| \leq \sup_{|\zeta - x| \leq \delta}
|f'(\zeta)| \cdot |\epsilon| \tag{1.9}$$

*for all $|\epsilon| \leq \delta$.*

*Proof.* By the Mean Value Theorem, there exists $\xi$ between $x$ and
$x + \epsilon$ such that $f(x + \epsilon) - f(x) = f'(\xi)\,\epsilon$.
Since $|\epsilon| \leq \delta$, the point $\xi$ lies in
$[x - \delta, x + \delta]$, so $|f'(\xi)| \leq \sup |f'(\zeta)|$.
The result follows. ∎
```

### Item types

| Item Type | Use for | Proof required? |
|-----------|---------|-----------------|
| **Definition** | Introducing a named concept, quantity, or notation | No |
| **Theorem** | A result requiring proof with broad applicability | Yes (at Build) |
| **Lemma** | A stepping stone to a theorem | Yes (at Build) |
| **Corollary** | A direct consequence of a theorem or lemma | Yes (usually brief) |
| **Proposition** | A narrower result, or a design principle | Optional |
| **Assumption** | A condition required by subsequent theorems | No |
| **Remark** | Commentary, motivation, convention choice, or connection | No |
| **Example** | A numerical or worked illustration | No |

### Numbering

Items are numbered by section: Definition N.k.1, Theorem N.k.2, etc. Each item type has its own counter within the section.

**Isolation principle.** Inserting a new theorem in §N.3 renumbers nothing in §N.4 or any other chapter. Each section is an independent numbering scope.

**Remarks** are unnumbered: use `**Remark.**` or `**Remark** (topic name)**.**`.

### Conditions

When a theorem has specific conditions (convergence domain, parameter range), state them explicitly after the main statement:

```markdown
*Conditions:* $|e'| < 1$ (satisfied by all terrestrial and planetary bodies).
```

---

## 5. Proofs

Three maturity levels, each with a distinct format.

### Build level — complete proof

```markdown
*Proof.* [Full proof text. May span multiple paragraphs with labeled steps.]

**Step 1.** [First step with justification.]

**Step 2.** [Second step.]

∎
```

The halmos symbol ∎ appears only after complete proofs. Never use `*Proof.*` without completing the proof and closing with ∎.

### Develop level — proof approach

```markdown
*Proof approach.* [High-level strategy: which technique, which key identity,
which reference. Enough for a reader to reconstruct the proof, but not the
proof itself. No ∎ marker.]
```

Verbatim example:

```markdown
*Proof approach.* Derive from the boundary condition $U = U_0$ on the ellipsoid
surface. Express the external gravitational potential in ellipsoidal harmonics
(§14.4.1), then convert to a spherical harmonic expansion in powers of
$(a/r)^{2n} P_{2n}(\sin\phi)$. The coefficient of the $P_2(\sin\phi)$ term, when
identified with the standard spherical harmonic definition of $J_2$, yields the
Brouwer formula. [Heiskanen and Moritz 1967, Eqs. 2-90, 2-92; Moritz 1980, p. 129]
```

### Draft level — proof sketch or absent

At Draft level, a proof may be entirely absent (the theorem is stated as a stub) or given as a one-sentence sketch:

```markdown
*Proof sketch.* Follows by integration by parts and the Wallis recurrence.
```

---

## 6. Examples

Examples follow the item they illustrate. They contain concrete numerical values, not just symbolic assertions.

```markdown
**Example 14.5.1** (Iteration for GRS80)**.** *Starting from
$e^2_0 = 3J_2 = 0.003\,247\,89$:*

| Iteration | $e^2$ | Digits correct |
|-----------|-------|----------------|
| 0 | $0.003\,247\,89$ | 1 |
| 1 | $0.006\,694\,4$ | 5 |
| 2 | $0.006\,694\,380\,02$ | 9 |
| 3 | $0.006\,694\,380\,022\,90$ | 13+ |

*Convergence to the GRS80 value $e^2 = 0.006\,694\,380\,022\,90$
[Moritz 1980, p. 131] is achieved in 3 iterations to double-precision
accuracy.*
```

**Requirements:**
- Every numerical value must cite its source: `[Moritz 1980, p. 131]` or `[NGA.STND.0036, 2014]`.
- Every example should be independently verifiable (a reader with a calculator, or the `verify_chN.py` script, can reproduce every number).
- Examples are numbered per section: Example N.k.j.

---

## 7. Equations

### Display equations

All display equations are numbered sequentially within the chapter:

```markdown
$$E = -\frac{\mu}{2a} \tag{8.1}$$
```

Tag format: `(ch.eq)` — equation 1 of chapter 8 is `(8.1)`.

### Inline equations

Use `$...$` for inline math. No numbering for inline equations.

### Referencing

Use **Eq.** consistently (not "Equation" in full):

```markdown
Substituting (14.7) into Eq. (14.5) yields...
```

When context is clear, the bare number `(8.3)` suffices.

---

## 8. Tables

### Data tables

```markdown
**Table N.k.j** (Descriptive title)**.**

| Column 1 | Column 2 | Column 3 |
|----------|----------|----------|
| data | data | data |
```

Center-align numeric columns where possible. Include units in column headers.

### Summary tables

Used in the Summary section to map results to downstream chapters. Follow the format from Ch 2 §2.11:

```markdown
| Chapter | Result used | This chapter's section |
|---------|------------|----------------------|
| Ch 8 | Theorem N.k.j | §N.k |
```

---

## 9. Cross-References

### Within the chapter

Reference by item number: `Theorem 14.3.1`, `Definition 14.4.2`, `Eq. (14.5)`.

### To other chapters

Use the format `Ch X, §X.Y` or `Ch X, Theorem X.Y.Z`. On first reference in a section, include a brief role:

```markdown
...the contraction mapping theorem (Ch 1, Theorem 1.6.1)...
```

### Formal cross-reference notation

- To a formal statement: `(Theorem 1.3.2)`, `(Definition 4.2.1)`, `(Lemma 7.5.1)`
- To a section: `(Ch 2, §2.3)` or `(§2.3)` when within the same chapter
- To an equation: `Eq. (8.3)` or `(8.3)`
- To an error note: `[A.5.2]` — the chapter number is embedded in the tag

### Forward references

When a later chapter's result is needed, state it as an assumption and cite the chapter. Do not forward-reference undefined concepts.

---

## 10. Summary Section (§N.last)

Every chapter ends with a summary section before Error Notes. It contains:

1. **Usage pattern**: How downstream chapters use this chapter's results (1–2 sentences).
2. **Downstream chapter reference table**:

```markdown
## §N.k Summary and Downstream Usage

Every subsequent chapter that performs a coordinate transform includes a
"State Framework" section showing the operation in both representations,
following this pattern: [brief description of the pattern].

**Downstream Chapter Reference**

| Chapter | Framework aspect used | This chapter's section |
|---------|---------------------|----------------------|
| Ch 8 (Keplerian orbit) | Perifocal frame vectors | §2.5, Thm 2.5.1 |
| Ch 9 (Kepler's equation) | Nonlinear solver sensitivity | §2.6, Thm 2.6.2 |
| Ch 13 (Geopotential) | Force gradient perturbation | §2.4, Def 2.4.1 |
```

For short chapters (under 5 sections), the summary may be a concluding Remark instead of a full section.

---

## 11. Error Notes

Placed at the end of the chapter, after the summary section and a `---` separator.

### Tag system

- **[M.N.k]** — Measurement error ($\sigma_m$) note $k$ in chapter $N$
- **[P.N.k]** — Precision error ($\delta_p$) note $k$ in chapter $N$
- **[A.N.k]** — Accuracy/model error ($\delta_a$) note $k$ in chapter $N$

Tags appear **inline** in the text at the point where the error source is introduced: `...this subtraction costs approximately 5.4 digits [P.14.2]...`

### End-of-chapter format

Each note collected at the end uses this structure:

```markdown
## Error Notes

**[M.1.1]** Measurement errors of physical constants. For the WGS-72
gravitational parameter $\mu = GM = 3.986008 \times 10^5$ km$^3$/s$^2$,
the measurement uncertainty is approximately $\sigma_m \sim 8 \times 10^{-3}$
km$^3$/s$^2$ (from satellite laser ranging). For defined constants
($a_E$, $1/f$, $\omega$), $\sigma_m = 0$ exactly.
*Remedy:* improved measurement campaigns (not a computational concern).
```

### Required sub-fields by maturity level

| Level | Issue | Magnitude | Remedy |
|-------|-------|-----------|--------|
| Build | Required | Required (quantified bound or formula) | Required |
| Develop | Required | Required (order-of-magnitude estimate) | May be placeholder |
| Draft | One-line description suffices | Optional | Optional |

### Placeholder marking

If an error note is not yet elaborated, mark it explicitly:

```markdown
**[P.N.k]** *Placeholder: precision error in the iterative solution for $e^2$
from $J_2$ (§14.5).*
```

Never leave a tag referenced inline without at least a placeholder in the Error Notes section.

---

## 12. Special Sections

Certain sections appear across multiple chapters with consistent formatting.

### State Framework

Chapters that perform coordinate transforms or state propagation include:

```markdown
### §N.k State Framework
```

This section expresses the chapter's operation in both representations — the dual quaternion pair $(\hat{M}, \hat{\Omega}_b)$ and the $7 \times 7$ matrix on $(r, v, 1)^T$ — following the framework of Ch 2. The section opens with a cross-reference to Ch 2 and follows the pattern:

1. Express as SU(2) transform or traceless Hermitian operation (Ch 2, §§2.2–2.4)
2. Classify as linear or nonlinear for error propagation (Ch 2, §2.6)
3. Apply Principle 1 (linear) or Principle 2 (nonlinear) for error rows
4. Show composition with adjacent pipeline steps

### Generalization

```markdown
### §N.k Generalization
```

Documents what additional terms, higher-order corrections, or alternative methods extend the standard SGP4 result. Connected to the tolerance parameter (Ch 3, Definition 3.4.1): every generalized routine accepts a tolerance that reproduces the original model at the standard value and computes additional terms at tighter tolerance.

### Full Accuracy Recovery

```markdown
### §N.k Full Accuracy Recovery
```

Documents: (a) the exact untruncated expression, (b) where and why SGP4 truncates, (c) the computational technique (continued fraction, higher-order series, numerical quadrature) for the full expression to arbitrary precision.

---

## 13. Mathematical Typesetting

### LaTeX conventions

- Inline math: `$...$`
- Display equations: `$$...$$` with `\tag{ch.eq}` for numbered equations
- Vectors: bold lowercase `$\mathbf{r}$`, `$\mathbf{v}$`
- Matrices: bold uppercase `$\mathbf{M}$`, `$\mathbf{R}$`
- Operators: upright `$\sin$`, `$\cos$`, `$\log$`, `$\exp$`
- Differentials: upright d: `$\mathrm{d}x$`
- Absolute value: `$|x|$`; norms: `$\|x\|$`
- Partial derivatives: `$\partial f / \partial x$` inline, `$\frac{\partial f}{\partial x}$` display
- Thousand separators in numerical values: `$6\,378\,137$` (thin space `\,`)

### Units

State units in square brackets after the first occurrence of a dimensional quantity: $a$ [m], $\mu$ [km$^3$/s$^2$]. Derivations use consistent units; conversions appear only at interfaces (TLE parsing, output formatting).

---

## 14. Global Notation Conventions

The following notation is consistent across all chapters:

| Symbol | Meaning |
|--------|---------|
| $a$ | Semi-major axis |
| $e$ | Eccentricity |
| $i$ | Inclination |
| $\Omega$ | Right ascension of the ascending node |
| $\omega$ | Argument of perigee |
| $M$ | Mean anomaly |
| $E$ | Eccentric anomaly |
| $\nu$ | True anomaly |
| $n$ | Mean motion |
| $p = a(1-e^2)$ | Semi-latus rectum |
| $\theta = \cos i$ | Cosine of inclination |
| $\eta = \sqrt{1-e^2}$ | Eccentricity function |
| subscript $_o$ | Value at epoch (letter "o", not zero) |
| $\sigma_m$ | Measurement error bound |
| $\delta_p$ | Precision error bound |
| $\delta_a$ | Accuracy/model error bound |
| $\mu = GM$ | Gravitational parameter |
| $a_E$ | Earth's equatorial radius |
| $J_n$ | Zonal harmonic coefficient of degree $n$ |
| $k_2 = \tfrac{1}{2} J_2 a_E^2$ | Gravitational coefficient |
| $(l, g, h, L, G, H)$ | Delaunay variables |
| $\hat{M}$ | Configuration dual quaternion ($SU(2)$ + translation) |
| $\hat{\Omega}_b$ | Velocity dual quaternion (angular velocity + translational velocity) |

---

## 15. Derivation Standard

**Every coefficient, formula, and bound must be DERIVED from first principles.** No result is accepted because it appears in reference code, matches a test case, or is stated without proof in a historical document.

The derivation chain is:
1. State assumptions (physical laws, mathematical hypotheses)
2. Apply standard mathematical techniques (algebra, calculus, perturbation theory)
3. Arrive at the result

If a standard result is used (e.g., Taylor's theorem, the contraction mapping theorem), cite it by name. It need not be re-proved, but its hypotheses must be verified for the application at hand.

**Vallado test cases** are cross-checks, not proofs. They may appear in verification sections but never in derivations.

---

## 16. Source Citations

Citations use the format `[Author Year, Eq. N]` or `[Author Year, §M]`:

- [Brouwer 1959, Eq. 29]
- [Lara 2021, §3.2]
- [Sneeuw 2022, Eq. 4.17]
- [Heiskanen and Moritz 1967, Eq. 2-70]
- [Moritz 1980, p. 131]

### Source hierarchy

**Primary sources** (clean digital, post-LaTeX typesetting):
- Lara, M. (2021) — Brouwer theory
- Na, K. S. et al. (2012) — SGP4 analysis
- Sneeuw, N. (2022) — Dynamic Satellite Geodesy lecture notes
- Moritz, H. (1980) — Geodetic Reference System 1980

**Historical sources** (cross-verification only, not primary derivation references):
- Brouwer, D. (1959) — original perturbation theory
- Lane, M. H. and Hoots, F. R. (1979) — AFGP4/SGP4 formulation
- Hoots, F. R. and Roehrich, R. L. (1980) — Spacetrack Report No. 3

**Specification documents:**
- NGA.STND.0036_1.0.0_WGS84 (2014) — WGS84 parameters
- DMA (1974), DTIC ADA110165 — WGS72 parameters

---

## 17. Tier Classification

Every constant and derived quantity receives a tier classification:

| Tier | Description | $\sigma_m$ | $\delta_p$ | $\delta_a$ | Example |
|------|-------------|-----------|-----------|-----------|---------|
| I | Exact by definition | 0 | repr. only | [A.N.k] | $a_E = 6378.135$ km (WGS-72) |
| II | Algebraic function of Tier I | 0 or inherited | 0 (algebraic) | Inherited | $e^2 = 2f - f^2$ |
| III | Computed via series/iteration | inherited | accumulated | inherited | $q_0$, $\gamma_e$, $R_2$ |
| IV | Model-dependent | inherited | accumulated | nonzero | Brouwer secular rates |

The tier is stated at the point of introduction, usually within or immediately after the defining equation.

---

## 18. Matched Pair Annotations

When a constant or truncation level is a **matched-pair value** — meaning it is "wrong" in an absolute sense but correct within the coupled TLE + SGP4 system — annotate it:

```markdown
$\mu = 398600.8$ km$^3$/s$^2$ **[MP]**
```

The **[MP]** tag indicates that this value must not be replaced with a "more accurate" constant without re-fitting the element set. The Matched Pair Principle is developed in Ch 3.

---

## 19. Writing Style

- Third person, present tense for mathematical statements: "The bound satisfies..." not "We compute..."
- Active voice for derivation steps: "Expanding the product yields..." not "The product was expanded..."
- No colloquialisms, no rhetorical questions
- Define every symbol before or upon first use
- No forward references to undefined concepts — if a later chapter's result is needed, state it as an assumption and cite the chapter
- No code snippets in chapter files (the code-to-theorem mapping is in Appendix C)
- Proofs are self-contained: a reader should not need external references to follow the argument, though citations indicate where a standard technique originates
- Section openings state the section's objective in one sentence: "This section derives [what] from [what]. The objective is to establish [result]."

---

## 20. File Naming

Chapter files: `ch01_three_errors.md` through `ch38_state_vector_output.md`

The two-digit chapter number ensures lexicographic ordering matches logical ordering. The suffix is a short snake_case descriptor.

Appendix files: `app_a_constants.md`, `app_b_notation.md`, `app_c_code_map.md`, `app_d_sources.md`

Plan files: `design/derivations/plan/ch01_three_errors.md`, etc. (one plan per chapter).

---

## 21. Plan File Format

Every chapter has a corresponding plan file in `design/derivations/plan/`. The plan file uses this structure:

```markdown
# Draft Plan: Chapter N — Title

## Objectives
[Numbered list of 5–10 specific, measurable objectives]

## Section Structure
### §N.1 Introduction
  - Forward-reference table entries: ...
  - Maturity: ...

### §N.2 Section Title
  - Definitions: N.2.1 (Name), N.2.2 (Name)
  - Theorems: N.2.1 (Name) — proof strategy: [brief]
  - Examples: N.2.1 (Name) — numerical content: [brief]
  - Error Notes: [P.N.1] — [brief description]

[... repeat for each section ...]

## Cross-References
**Uses (backward):**
| Source | Section | Role |
|--------|---------|------|

**Feeds (forward):**
| Target | Section | Role |
|--------|---------|------|

## Error Notes
| Tag | Type | § | Description |
|-----|------|---|-------------|

## Estimated Count
| Item | Count |
|------|-------|
| Definitions | X |
| Theorems | X |
| Lemmas | X |
| Corollaries | X |
| Propositions | X |
| Examples | X |
| Error Notes | X |
| Equations | ~X |
| Sections | X |

## Maturity
| Section | Level | Notes |
|---------|-------|-------|
```

---

## 22. Quality Checklist

Before declaring a chapter complete at any maturity level, verify:

### Draft

- [ ] §N.1 has opening paragraph, roadmap, and forward-reference table
- [ ] Maturity declaration present
- [ ] All sections have objective statements
- [ ] All Definitions are numbered and named (statement may be placeholder)
- [ ] All Theorems are numbered and named (proof not required)
- [ ] Error note tags are placed inline; placeholders exist in Error Notes section
- [ ] Notation table exists with all symbols

### Develop (all of Draft, plus)

- [ ] All Theorem statements are complete (full mathematical content)
- [ ] All Theorems have proof approaches
- [ ] All Examples have numerical values with source citations
- [ ] Error notes have Issue and Magnitude sub-fields
- [ ] Cross-references to other chapters are specific (§X.Y, Theorem X.Y.Z)
- [ ] Maturity declaration is accurate

### Build (all of Develop, plus)

- [ ] All proofs are complete with ∎
- [ ] Error notes have Issue, Magnitude, and Remedy sub-fields
- [ ] Summary section with downstream chapter reference table
- [ ] All numerical values independently verified (verification script exists)
- [ ] No placeholder markers remain
