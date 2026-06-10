# Draft Plan: Chapter 3 — The Matched Pair Principle

**Part I: Mathematical Foundations**

## Objectives

1. Formally define the matched pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ — the coupled system of a propagation model and the orbital elements fitted to it — and prove that substituting "better" components without re-fitting can degrade accuracy (Theorem 3.2.1).
2. Introduce **LegacyModelValue** as a fourth value category (alongside DefinedValue, MeasuredValue, ModelValue) for constants that are frozen-wrong but specification-exact within the matched pair; derive their error classification (Theorem 3.3.1).
3. Define the tolerance parameter $\tau$ and the closure/lambda architecture: every generalized computation accepts $\tau$ at initialization, captures it in a closure, and iterates internally until its remainder bound falls below $\tau$. Define $\tau_{\mathrm{standard}}$ as agreement with the Hoots and Roehrich (1980) result to within $\delta_a(\mathcal{M})$.
4. Define two operating modes — standard (preserves matched pair) and enhanced (breaks it by design) — and characterize their respective error budgets.
5. Resolve the forward references [A.2.1] and [M.2.1] from Chapter 2 (TEME frame and GMST coefficients), and the **[MP]** and $\tau$ previews from Chapter 1.
6. State the downstream verification requirement for all Chapters 4–35 computations that replace SR3 computations.

---

## Section Structure

### §3.1 Introduction

- Narrative orientation only; no formal items.
- Previews the four main results: Matched Pair (§3.2), LegacyModelValue (§3.3), Tolerance Parameter (§3.4), Two Operating Modes (§3.5).
- Resolves forward references from Ch 1 (§1.2 Remark after Def. 1.2.5; §1.5 Remark after Thm. 1.5.3) and Ch 2 ([A.2.1], [M.2.1]).
- Key statement: "The results of this chapter are not algebraic derivations. They are constraints derived from the structure of the least-squares fitting process."

---

### §3.2 The Matched Pair

- **Definitions:**
  - **Definition 3.2.1** (Matched pair) — Propagation model $\mathcal{M}$ as triple $(\{c_i\}, \{N_j\}, \{f_k\})$; element set $\mathbf{e}_{\mathcal{M}}$ as least-squares minimizer (Eq. 3.1); matched pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$. Formally defines the **[MP]** annotation.
- **Theorems:**
  - **Theorem 3.2.1** (Degradation by substitution) — Propagating $\mathbf{e}_{\mathcal{M}}$ under modified model $\mathcal{M}'$ can increase accuracy error (Eqs. 3.2–3.3). Proof by 1-D analogy: compensating bias $\beta$ (Eqs. 3.4–3.5); matched-pair residual is $O(|c_0 - c_{\mathrm{true}}|^2)$ (Eq. 3.6); substitution residual is first-order (Eq. 3.7). Existence result: only triggers for components the fitting process compensated for.
- **Remarks:**
  - Remark 1 (after Def. 3.2.1): **[MP]** tag formally grounded in Def. 3.2.1.
  - Remark 2 (after Def. 3.2.1): Minimization (3.1) is schematic; actual fitting involves differential correction with weighting.
- **Examples:**
  - **Example 3.2.1** (Gravitational parameter mismatch) — $\mu_{\mathrm{SGP4}} = 398600.5$ vs. $\mu_{\mathrm{WGS84}} = 398600.4418$; $\Delta\mu = 0.0582$ km$^3$/s$^2$; radial bias formula (Eqs. 3.10–3.11); table with LEO/MEO/GEO effects; matched-pair double-correction mechanism.
  - **Example 3.2.2** (G-function polynomial approximations) — Hough cubic Hansen coefficient approximations; combined compensation for $\epsilon_G + \epsilon_P$; eccentricity dependence of effect; reference to Ch 15 for full analysis.
  - **Example 3.2.3** (Brouwer secular rate truncation) — SGP4 truncated at $J_2^2 + J_4$; omitted terms $\sim J_2 \approx 10^{-3}$; linear secular drift formula (Eq. 3.12); degradation grows with propagation time.

---

### §3.3 LegacyModelValue

- **Definitions:**
  - **Definition 3.3.1** (LegacyModelValue) — Three-condition definition: component of matched-pair model; differs from best measurement; retained because elements were fitted to it. **[MP]** annotation and Tier I classification within the matched-pair system.
- **Theorems:**
  - **Theorem 3.3.1** (Error classification of LegacyModelValue) — Three-part: (i) $\sigma_m(c) = 0$; (ii) $\delta_p(c) \leq |c_0| \cdot \epsilon_{\mathrm{mach}}$; (iii) $\delta_a(c) = \delta_a(\mathcal{M})$. Proof by argument from structure: (i) specification vs. measurement; (ii) IEEE 754 representation; (iii) avoidance of double-counting.
- **Remarks:**
  - Remark 1 (after Def. 3.3.1): Contrast with DefinedValue — "exactness" from matching, not geodetic definition.
  - Remark 2 (after Thm. 3.3.1): $\Delta\mu = 0.0582$ km$^3$/s$^2$ is real; the theorem classifies it, not negates it.
  - Remark 3 (Contrast with value categories): Comparison table of all four value categories with error signatures and examples.
  - Remark 4 (after Table 3.3.1): $\omega_E$ is a DefinedValue, not a LegacyModelValue — no discrepancy at relevant precision.
  - Remark 5 (after Table 3.3.1): SGP4 $J_2$ in WGS-84 mode is a specification bypassing ellipsoid derivation.
- **Tables:**
  - **Table 3.3.1** (LegacyModelValues in SGP4) — Seven rows: $\mu$ (WGS-72), $\mu$ (WGS-84), $J_2$, $J_3$, $J_4$, $k_e$, GMST coefficients, $q_0$/$s$. Columns: SGP4 value [MP], modern best, $|c_0 - c_{\mathrm{true}}|$, radial effect, source.

---

### §3.4 The Tolerance Parameter Principle

- **Definitions:**
  - **Definition 3.4.1** (Tolerance parameter) — $\tau > 0$ in physical units; initializer captures $\tau$ in a closure; $\lambda$ iterates internally until remainder bound $< \tau$; caller never sees $\tau$.
  - **Definition 3.4.2** (Standard tolerance) — $\tau_{\mathrm{standard}}$: agreement with SR3 computation to within $\delta_a(\mathcal{M})$ (Eq. 3.14). Explicitly does NOT promise bitwise reproduction.
- **Propositions:**
  - **Proposition 3.4.1** (Tolerance parameter across computational patterns) — Uniform application to: series with remainder bounds (stopping at $|R_N| < \tau$, refs. Thm. 1.5.2–1.5.3 → Ch 5); continued fractions/Padé ($|C_n - C_{n-1}| < \tau$, → Ch 4); iterative solvers (contraction mapping bound $< \tau$, Thm. 1.6.1 → Ch 9). Table format. Proof by case: SR3 truncation point analysis for each pattern.
- **Remarks:**
  - Remark 1 (after Def. 3.4.1): Mathematical vs. architectural separation of concerns.
  - Remark 2 (after Def. 3.4.2): Bitwise reproduction caveat — dual quaternions and arbitrary-precision arithmetic differ from SR3 at precision level.
  - Remark 3 (Not a term count): Algorithm self-monitors remainder; reports actual remainder as $\delta_p$. Back-reference to Ch 1 §1.5.
- **Examples:**
  - **Example 3.4.1** (Closure architecture) — Concrete geometric-series $\lambda$ construction; $N = \min\{n : |a_{n+1}|/(1-r) < \tau\}$ (Eq. 3.13).

---

### §3.5 Two Operating Modes

- **Definitions:**
  - **Definition 3.5.1** (Standard mode) — LegacyModelValues + $\tau_{\mathrm{standard}}$; preserves matched pair; error budget table ($\sigma_m$, $\delta_p$, $\delta_a \approx 1$ km per [A.3.1]).
  - **Definition 3.5.2** (Enhanced mode) — MeasuredValues + $\tau < \tau_{\mathrm{standard}}$; breaks matched pair by design; error budget table; output will NOT match standard SGP4 test cases.
- **Remarks:**
  - Remark 1 (after Def. 3.5.1, Representation caveat): Not bitwise identical to SR3; dual quaternion and arbitrary-precision differences are at $\delta_p$ level, below accuracy floor.
  - Remark 2 (after Def. 3.5.2): Same callable interfaces in both modes; mode selection is purely an initialization choice; reference to Ch 36 injectable lambda architecture.
  - Remark 3 (after Def. 3.5.2): Standard mode is default; enhanced mode is opt-in for comparison studies, research, and validation.

---

### §3.6 Breaking the Matched Pair Safely

- **Propositions:**
  - **Proposition 3.6.1** (Safe modification) — $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ can be replaced by $(\mathcal{M}', \mathbf{e}_{\mathcal{M}'})$ if elements are re-fitted to $\mathcal{M}'$. Proof: re-fitting (3.1) under $\mathcal{M}'$ produces valid matched pair; if $\mathcal{M}'$ has fewer systematic errors, $\delta_a(\mathcal{M}', \mathbf{e}_{\mathcal{M}'}) \leq \delta_a(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ (Eq. 3.3).
  - **Proposition 3.6.2** (Unsafe modification) — Substituting a single component while keeping $\mathbf{e}_{\mathcal{M}}$ breaks the pair; governed by Theorem 3.2.1. Proof: direct application of Thm. 3.2.1 (Eqs. 3.6–3.7).
- **Remarks:**
  - Remark 1 (after Prop. 3.6.1): In practice, re-fitting requires 18 SDS infrastructure; individual users cannot re-fit TLE elements.
- **Examples:**
  - **Example 3.6.1** (IAU 2006 nutation — resolves [A.2.1]) — TEME frame accuracy $\sim 0.1$ arcsec; positional effect at orbital radius formula (Eq. 3.15); table with LEO/GPS/GEO effects (3.4, 13, 21 m); all below 1 km accuracy floor but move prediction in wrong direction; enhanced-mode remedy.
  - **Example 3.6.2** (IERS Earth rotation — resolves [M.2.1]) — GMST polynomial (Eq. 3.16); error propagation through frame rotation (Eq. 3.17, ref. Ch 2 Thm. 2.2.11); IAU 1982 vs. IERS 2010 difference $\sim 10^{-6}$ rad for few-day arcs; enhanced-mode remedy.

---

### §3.7 The Matched Pair as a Constraint on Approximation

- **Definitions:**
  - **Definition 3.7.1** (Matched-pair compatibility) — $g(\tau)$ is matched-pair compatible with $f_{\mathrm{SGP4}}$ iff $|g(\tau_{\mathrm{standard}}) - f_{\mathrm{SGP4}}| \leq \delta_a(\mathcal{M})$ (Eq. 3.18). No constraint at tighter $\tau$.
- **Corollaries:**
  - **Corollary 3.7.1** (Downstream verification requirement) — Every Ch 4–35 generalized computation must verify Eq. (3.18) on representative orbital elements ($e \in [0, 0.9]$, $i \in [0°, 180°]$). Table of directly constrained chapters: Ch 4 (Padé), Ch 5 (series), Ch 9 (Kepler solver), Ch 15 (Hansen coefficients), Ch 25 (equation of center).
- **Remarks:**
  - Remark 1 (after Def. 3.7.1): Constraint is one-sided — agree with SR3 at $\tau_{\mathrm{standard}}$; improve beyond it for tighter $\tau$.
  - Remark 2 (Why this matters): Continued fraction convergents are different numbers from Taylor partial sums; faster convergence ≠ matched-pair compatible without verification.
  - Remark 3 (after Cor. 3.7.1): Constraint is "a design discipline": stated, checked per computation, documented — but almost always satisfied in practice.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Def 1.2.5 | §3.2 | [MP] annotation on Tier I constants |
| Ch 1, Thm 1.5.2–1.5.3 | §3.4 | Series remainder bounds for tolerance parameter |
| Ch 1, Thm 1.6.1 | §3.4 | Iterative convergence for tolerance parameter |
| Ch 2, Def 2.5.2 | §3.6 | TEME frame definition (resolved by Example 3.6.1) |
| Ch 2, Thm 2.2.11 | §3.6 | Angle sensitivity (resolved by Example 3.6.2) |
| Ch 2, §2.3 | §3.5 | Dual quaternion representation caveat |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 4, 5, 9, 15, 25 | §3.7 | Matched-pair compatibility verification (Cor 3.7.1) |
| Ch 15 | §3.2 | Hansen coefficient degradation example |
| Ch 16–17 | §3.2 | Brouwer secular rate degradation example |
| Ch 29 | §3.6 | GMST polynomial matched-pair constraint |
| Ch 36 | §3.4–3.5 | Injectable lambda architecture, standard/enhanced modes |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.3.1] | A | §3.5, §3.7 | SGP4 accuracy floor (~1 km position, ~1 m/s velocity); combined effect of [A.1.1]–[A.1.4] |
| [M.3.1] | M | §3.3 | LegacyModelValue convention: $\sigma_m = 0$, discrepancy absorbed into $\delta_a$ |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 7 |
| Theorems | 2 |
| Lemmas | 0 |
| Corollaries | 1 |
| Propositions | 3 |
| Examples | 6 |
| Error Notes | 2 |
| Equations | ~18 |
| Sections | 7 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §3.1 | Develop | Introduction narrative |
| §3.2 | Develop | Matched pair definition and degradation theorem |
| §3.3 | Develop | LegacyModelValue classification |
| §3.4 | Develop | Tolerance parameter framework |
| §3.5 | Develop | Standard and enhanced modes |
| §3.6 | Develop | Safe/unsafe modification examples |
| §3.7 | Develop | Matched-pair compatibility constraint |
