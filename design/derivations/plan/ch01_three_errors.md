# Draft Plan: Chapter 1 — The Three Fundamental Errors

**Part I: Mathematical Foundations**

**Phase:** Build (chapter is complete)

**Implementation target:** `tracked_value.h` (error-tracking scaffold for all subsequent chapters)

---

## Objectives

1. Define the three fundamental and independent error categories — measurement error ($\sigma_m$), precision error ($\delta_p$), and accuracy error ($\delta_a$) — and establish that they cannot substitute for one another.
2. Formalize the tracked-value quadruple $\mathcal{V} = (v, \sigma_m, \delta_p, \delta_a)$ and the tier classification of constants and derived quantities.
3. Derive error propagation rules for all arithmetic operations (addition, subtraction, multiplication, division, composition) over both $\mathbb{R}$ and $\mathbb{C}$.
4. Derive error propagation rules for all transcendental functions needed in SGP4: sine, cosine, atan2, square root, exponential, logarithm, power, modulus, modular reduction, and arcsine.
5. Classify series truncation error as either precision or accuracy depending on whether the series evaluates a known function or represents a perturbation expansion.
6. Derive precision error bounds for iterative solvers: contraction mapping (linear convergence), Newton's method (quadratic), and general $p$-th order methods.
7. Quantify subtractive cancellation via the condition number of subtraction and the digit-loss formula.
8. Define worst-case (additive) and RSS error combination and prove independence of the three categories.
9. Define the reliable-digits metric $\mathrm{rd}(v)$ and the comparison reliability criterion as the decision variables for formula-switching and branch-taking throughout the propagator.

---

## Notation Table

*(As written in the chapter — reproduced for completeness.)*

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathbb{F}$ | Scalar field: $\mathbb{R}$ or $\mathbb{C}$ | §1.3 |
| $|z|$ | Absolute value or modulus | §1.3 |
| $\|\mathbf{v}\|$ | Euclidean norm (vectors/matrices; not used in Ch 1) | Ch 2 |
| $v$ | Computed or measured scalar | §1.2 |
| $\tilde{v}$ | True (unknowable) value | §1.2 |
| $\hat{v}$ | Result of computational evaluation | §1.2, Def. 1.2.2 |
| $v_{\mathrm{model}}$ | Exact model prediction (infinite precision) | §1.2, Def. 1.2.3 |
| $\sigma_m$ | Measurement error bound | §1.2, Def. 1.2.1 |
| $\delta_p$ | Precision error bound | §1.2, Def. 1.2.2 |
| $\delta_a$ | Accuracy error bound | §1.2, Def. 1.2.3 |
| $\delta$ | Generic error bound (any category) | §1.3 |
| $\delta_{\mathrm{total}}$ | Total error bound $= \sigma_m + \delta_p + \delta_a$ | §1.8, Def. 1.8.1 |
| $\delta_{\mathrm{RSS}}$ | RSS error bound $= \sqrt{\sigma_m^2 + \delta_p^2 + \delta_a^2}$ | §1.8, Def. 1.8.2 |
| $\mathcal{V}$ | Tracked value quadruple $(v, \sigma_m, \delta_p, \delta_a)$ | §1.2, Def. 1.2.4 |
| $\mathrm{fl}(v)$ | IEEE 754 floating-point representation of $v$ | §1.2 |
| $\mathrm{ulp}(v)$ | Unit in the last place | §1.2 |
| $\epsilon_{\mathrm{mach}}$ | Machine epsilon: $2^{-52}$ for binary64 | §1.2 |
| $\epsilon$ | Convergence tolerance for an iterative solver | §1.6 |
| $L$ | Lipschitz (contraction) constant, $L < 1$ | §1.6, Asm. 1.6.1 |
| $\kappa$ | Condition number of an operation | §1.7, Thm. 1.7.1 |
| $\mathrm{rd}(v)$ | Reliable decimal digits in $v$ | §1.9, Def. 1.9.1 |
| $R_N$ | Remainder of a series truncated after $N$ terms | §1.5 |
| $S_N$ | Partial sum: $S_N = \sum_{k=0}^{N} a_k$ | §1.5 |
| $\nabla f$ | Gradient vector | §1.4, Thm. 1.4.2 |
| $\boldsymbol{\epsilon}$ | Multivariate perturbation vector | §1.4, Thm. 1.4.2 |
| $\kappa_{\mathrm{rel}}(f, v)$ | Relative condition number $= |v \cdot f'(v) / f(v)|$ | §1.9, Thm. 1.9.1 |
| $n$ (in §1.4) | Integer quotient $\lfloor x/y \rfloor$ in modular reduction | §1.4, Cor. 1.4.9 |
| $p$ (in §1.6) | Order of convergence of an iterative method | §1.6, Thm. 1.6.3 |
| $C$ (in §1.6) | Convergence constant in $p$-th order methods | §1.6, Thm. 1.6.3 |
| $m$, $M$ (in §1.6) | Bounds on $|f'|$ and $|f''|$ in Newton's method | §1.6, Thm. 1.6.2 |
| **[MP]** | Matched Pair annotation | §1.2, Ch. 3 |

---

## Section Structure

### §1.1 Introduction

Narrative only. Motivates the three-error framework from the structure of orbit propagation (physical constants, model approximations, computational process). States that the three sources are fundamentally independent and previews the chapter structure. Forward-references Ch 2 (vectors/matrices), Ch 3 (matched-pair constraint), and Ch 4 (fast-convergence series alternatives).

- No definitions, theorems, or error notes.

---

### §1.2 Fundamental Definitions

- **Definition 1.2.1** (Measurement error). Bound $\sigma_m \geq 0$ on $|v_{\mathrm{meas}} - \tilde{v}|$; equation (1.1). Followed by remark distinguishing exact-by-definition constants ($\sigma_m = 0$) from measured constants ($\sigma_m > 0$). Cites [M.1.1].
- **Definition 1.2.2** (Precision error). Bound $\delta_p \geq 0$ on $|\hat{v} - v|$; equation (1.2). Four listed sources: representation error [P.1.1], rounding accumulation, series truncation, iterative convergence.
- **Definition 1.2.3** (Accuracy error). Bound $\delta_a \geq 0$ on $|v_{\mathrm{model}} - v_{\mathrm{true}}|$; equation (1.3). Three examples citing [A.1.1], [A.1.2], [A.1.3]. Followed by remark distinguishing precision truncation (known function) from accuracy truncation (perturbation expansion).
- **Definition 1.2.4** (Tracked value). Ordered quadruple $\mathcal{V} = (v, \sigma_m, \delta_p, \delta_a)$; equations (1.4) and (1.5) (total-error bound on $|\tilde{v} - v|$).
- **Definition 1.2.5** (Tier classification). Four-tier table (Tier I–IV) based on which error categories are nonzero. Tier promotion rule for model-dependent computations. Followed by remark on the SGP4 matched-pair system and the [MP] annotation.

---

### §1.3 Propagation Through Arithmetic Operations

Opening remark on field generality: all results hold for $\mathbb{F} = \mathbb{R}$ and $\mathbb{F} = \mathbb{C}$; uses triangle inequality, reverse triangle inequality, and multiplicative modulus property. Proofs are stated for a generic $\delta$ applying to any of the three categories.

- **Theorem 1.3.1** (Addition and subtraction). $|\tilde{z} - z| \leq \delta_x + \delta_y$; equation (1.6). Proof by triangle inequality. Remark: absolute error bound is the same for addition and subtraction; relative error can be large (subtractive cancellation treated in §1.7).
- **Theorem 1.3.2** (Multiplication). $|\tilde{z} - z| \leq |x|\delta_y + |y|\delta_x + \delta_x\delta_y$; equation (1.7). Proof by expansion with $\tilde{x} = x + \epsilon_x$. Remark: cross term $\delta_x\delta_y$ is retained for rigor; negligible when errors are small.
- **Theorem 1.3.3** (Division). Bound (1.8) requiring $|y| > \delta_y$. Proof via reverse triangle inequality on $|\tilde{y}|$. Remark: condition $|y| > \delta_y$ is a mathematical requirement (not a guard); the denominator condition arises in Kepler solvers, density ratios, and Lyddane modification.
- **Theorem 1.3.4** (Composition of error bounds). Chained operations produce valid output bounds by induction, provided each intermediate step uses the propagation rules. Proof by monotonicity of the propagation rules. Remark: guarantees validity of the full SGP4 pipeline error budget.

---

### §1.4 Propagation Through Transcendental Functions

- **Theorem 1.4.1** (Derivative bound for scalar functions). $|f(x+\epsilon) - f(x)| \leq \sup|f'| \cdot |\epsilon|$; equation (1.9). Two proofs: MVT for $\mathbb{R}$, contour integral for $\mathbb{C}$. Remark: applies to each error category independently.
- **Theorem 1.4.2** (Multivariate error propagation). Bound (1.10) via telescoping-sum argument. Proof using Theorem 1.4.1 coordinate-by-coordinate. Remark: arithmetic operations are special cases; used in Ch 8, Ch 18–20, Ch 30, Ch 38.
- **Corollary 1.4.1** (Sine bound). $|\sin(x+\epsilon) - \sin(x)| \leq |\cos x|\delta + \delta^2/2$, capped at 2; equation (1.11). Taylor remainder proof for $\mathbb{R}$. Domain note for $\mathbb{C}$.
- **Corollary 1.4.2** (Cosine bound). $|\cos(x+\epsilon) - \cos(x)| \leq |\sin x|\delta + \delta^2/2$, capped at 2; equation (1.12). Analogous proof. Domain note for $\mathbb{C}$.
- **Corollary 1.4.3** (Two-argument arctangent bound). Bound (1.13) in terms of $|x|/r^2$ and $|y|/r^2$; requires $r^2 > (\delta_x + \delta_y)^2$. Proof via partial derivatives. Domain note: atan2 is real-only; not used in SU(2) framework.
- **Corollary 1.4.4** (Square root bound). $|\sqrt{x+\epsilon} - \sqrt{x}| \leq \delta/(2\sqrt{x-\delta})$; equation (1.14). Proof by monotone derivative supremum. Remark on the $\delta \geq x$ edge case. Domain note on branch cut.
- **Corollary 1.4.5** (Exponential bound). $|e^{x+\epsilon} - e^x| \leq e^x(e^\delta - 1)$; equation (1.15). Proof by convexity. Domain note: valid for all $z \in \mathbb{C}$.
- **Corollary 1.4.6** (Logarithm bound). $|\ln(x+\epsilon) - \ln x| \leq \ln(x/(x-\delta))$; equation (1.16). Proof by concavity; tighter than derivative bound. Domain note on branch cut.
- **Corollary 1.4.7** (Power function bound). Bound (1.17) via composition of Corollaries 1.4.6 and 1.4.5. Remark: used for density model ($\tau = 4$) and mean motion ($-3/2$ power). Domain note on branch cut.
- **Corollary 1.4.8** (Modulus bound). $||x+\epsilon| - |x|| \leq \delta$; equation (1.18). Proof by reverse triangle inequality; valid for both $\mathbb{R}$ and $\mathbb{C}$.
- **Corollary 1.4.9** (Modular reduction bound). $|\tilde{z} - z| \leq \delta_x + |n|\delta_y$ where $n = \lfloor x/y \rfloor$; equation (1.19). Proof by writing $z = x - ny$. Remark on the integer-jump edge case (angle wrapping). Domain note: real-only.
- **Corollary 1.4.10** (Arcsine bound). Bound (1.20); requires $\delta < 1 - |x|$. Proof by monotone derivative supremum. Remark: divergence at $|x| \to 1$ motivates atan2 over arcsin for angle recovery. Domain note.

---

### §1.5 Propagation Through Series Truncation

- **Assumption 1.5.1.** Series $S = \sum_{k=0}^\infty a_k$ converges absolutely.
- **Theorem 1.5.1** (Classification of truncation error). Remainder $R_N$ is a precision error [P.1.2] when the series evaluates a known function; an accuracy error [A.1.4] when it is a perturbation expansion. Proof by appeal to Definitions 1.2.2 and 1.2.3.
- **Theorem 1.5.2** (Leibniz alternating series bound). $|R_N| \leq b_{N+1}$; equation (1.21). Proof by bracketing argument on partial sums.
- **Theorem 1.5.3** (Ratio test remainder bound). $|R_N| \leq |a_{N+1}|/(1-r)$; equation (1.22). Proof by geometric bound. Remark: applies to $q_0$ series, eccentricity functions $G_{lpq}(e)$, equation-of-center expansion. Second remark: tolerance parameter (Ch 4) governs truncation index $N$; at standard tolerance reproduces Hoots and Roehrich 1980.

---

### §1.6 Propagation Through Iterative Convergence

- **Assumption 1.6.1.** $g: D \to D$ is a contraction mapping with Lipschitz constant $L < 1$.
- **Theorem 1.6.1** (Contraction mapping error bound). $|x_k - x^*| \leq (L/(1-L))|x_k - x_{k-1}|$; equation (1.23). Proof by geometric series. Remark [P.1.3]: for stopping criterion $|x_k - x_{k-1}| < \epsilon$, precision error bounded by $L\epsilon/(1-L)$.
- **Theorem 1.6.2** (Newton's method convergence). $|x_{k+1} - x^*| \leq (M/2m)|x_k - x^*|^2$; equation (1.24). Proof by Taylor expansion with Lagrange remainder (MVT for $\mathbb{R}$; integral form for $\mathbb{C}$). Remark: digit-doubling; last correction $|c_k|$ used as $\delta_p$.
- **Theorem 1.6.3** (General $p$-th order convergence). Bound (1.25) and (1.26). Proof via reverse triangle inequality and geometric-series argument. Remark: table of methods (fixed-point, Newton, Halley, Householder) with orders and SGP4 applications.

---

### §1.7 Subtractive Cancellation and Digit Loss

- **Theorem 1.7.1** (Condition number of subtraction). $\kappa(x,y) = (|x|+|y|)/|x-y|$; equation (1.27). Proof via relative-error decomposition.
- **Theorem 1.7.2** (Digit loss). Digits lost $\approx \log_{10}\kappa$; equation (1.28). Proof by substituting $d$-digit relative error.
- **Example** (unnamed, inline). Computing $\eta = \sqrt{1-e^2}$ for $e = 0.999$: cancellation in $1 - e^2$ loses $\approx 3$ digits; subsequent $\sqrt{\cdot}$ recovers $\approx 0.5$ digits. Cites [P.1.4]. Remark: the remedy is reformulation ($1-e^2 = (1-e)(1+e)$), not higher precision; $\mathrm{rd}(z) < d_{\min}$ triggers the switch.

---

### §1.8 Combining Error Bounds

- **Definition 1.8.1** (Worst-case combination). $\delta_{\mathrm{total}} = \sigma_m + \delta_p + \delta_a$; equation (1.29). Proof by four-link chain $\tilde{v} \to v_{\mathrm{model}} \to v_{\mathrm{model}}^{\mathrm{meas}} \to v$ and triangle inequality.
- **Definition 1.8.2** (Root-sum-square combination). $\delta_{\mathrm{RSS}} = \sqrt{\sigma_m^2 + \delta_p^2 + \delta_a^2}$; equation (1.30). Expected-case bound when categories are independent.
- **Proposition 1.8.1** (Independence of the three categories). Physical argument: the three categories have disjoint mechanisms and disjoint remedies. Not a mathematical theorem. Remark: worst-case bound remains the rigorous bound; Proposition 1.8.1 justifies RSS as expected-case only. Second remark: within-category correlations (e.g., correlated TLE measurement errors) do not invalidate the additive rule.

---

### §1.9 Reliable Digits

- **Definition 1.9.1** (Reliable digits). $\mathrm{rd}(v) = \lfloor -\log_{10}(\delta_{\mathrm{total}}/|v|) \rfloor$; equation (1.31). Edge cases: exact value ($\mathrm{rd} = +\infty$), error-dominated ($\mathrm{rd} \leq 0$).
- **Theorem 1.9.1** (Reliable digits under univariate functions). $\mathrm{rd}(\mathcal{V}_2) \geq \mathrm{rd}(\mathcal{V}_1) - \lceil \log_{10}\kappa_{\mathrm{rel}}(f, v_1) \rceil$; equation (1.32). Proof via relative-error transformation under Theorem 1.4.1. Remark: table of $\kappa_{\mathrm{rel}}$ for common operations (addition, near-cancellation subtraction, $\sin$ at zeros, $\sin$ at extrema, $\sqrt{\cdot}$, $\mathrm{atan}$); basis for formula-switching criterion.
- **Definition 1.9.2** (Reliable digits for near-zero values). When $|v| < \delta_{\mathrm{total}}$, $\mathrm{rd}(v) \leq 0$; absolute error $\delta_{\mathrm{total}}$ is the meaningful measure. Remark: near-zero is not meaningless (e.g., $e = 0$ or $i = 0$ are physically correct answers); downstream chapters must use absolute bounds when the quantity can legitimately be zero (eccentricity, inclination, short-period corrections).
- **Definition 1.9.3** (Comparison reliability). $x > y$ is reliable if and only if $|x - y| > \delta_{\mathrm{total},x} + \delta_{\mathrm{total},y}$; equation (1.33). Remark: applies to perigee altitude thresholds (98, 156, 220 km), period threshold (225 min), eccentricity threshold ($e < 10^{-6}$) — all branch decisions in the SGP4 propagator.

---

### Error Notes Section

Unnumbered section listing all error notes in full, with remedy statements.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| (none) | — | Ch 1 is the foundational chapter; no backward dependencies |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 2 | §1.3–1.4 | Extends arithmetic/transcendental error rules to vectors, matrices, dual quaternions |
| Ch 3 | §1.2 | [MP] annotation on Tier I constants; tolerance parameter |
| Ch 4 | §1.5 | Tolerance parameter governing series truncation; fast-convergence alternatives |
| Ch 7 | §1.4 | Modular reduction bound (Cor 1.4.9), atan2 bound (Cor 1.4.3) |
| Ch 8 | §1.4 | Multivariate propagation (Thm 1.4.2) for element-to-state conversion |
| Ch 9–10 | §1.6 | Convergence error bounds (Thm 1.6.1–1.6.3); anomaly quadrant (Cor 1.4.3) |
| Ch 13 | §1.5 | Series truncation classification (Thm 1.5.1) |
| Ch 17–19 | §1.5 | Perturbation expansion truncation as accuracy error (Thm 1.5.1) |
| Ch 21–22 | §1.2, §1.4 | Density model accuracy error (Def 1.2.3); power-function bound (Cor 1.4.7) |
| Ch 32 | §1.7 | Subtractive cancellation diagnostics |
| Ch 35 | §1.9 | Period threshold branch decision (Def 1.9.3) |
| Ch 38 | §1.3 | Full error budget via composition (Thm 1.3.4) |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.1.1] | M | §1.2 | Measurement errors of physical constants ($\mu$: $\sigma_m \sim 8 \times 10^{-3}$ km³/s²; defined constants: $\sigma_m = 0$) |
| [P.1.1] | P | §1.2 | IEEE 754 binary64 representation error ($\leq |v| \cdot 2^{-53}$, ~15.9 digits) |
| [P.1.2] | P | §1.5 | Series truncation as precision error (remainder bounded by Thm 1.5.2/1.5.3) |
| [P.1.3] | P | §1.6 | Iterative convergence tolerance ($\delta_p \leq L\epsilon/(1-L)$) |
| [P.1.4] | P | §1.7 | Subtractive cancellation in $1 - e^2$ (up to 5 digits lost at $e = 0.99999$) |
| [A.1.1] | A | §1.5 | Omitted zonal harmonics ($\delta_a \sim 50$–$200$ m/day for 400 km LEO) |
| [A.1.2] | A | §1.5 | Power-law density model (20–50% density residuals) |
| [A.1.3] | A | §1.5 | Brouwer perturbation truncation at $J_2^2 + J_4$ ($\delta_a \sim 10$–$100$ m/day) |
| [A.1.4] | A | §1.5 | Perturbation expansion as model truncation (drops real physics at each order) |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 10 |
| Theorems | 15 |
| Lemmas | 0 |
| Corollaries | 10 |
| Propositions | 1 |
| Examples | 1 |
| Error Notes | 9 |
| Equations | ~33 |
| Sections | 9 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §1.1 | Build | Narrative introduction |
| §1.2 | Build | Fundamental definitions complete |
| §1.3 | Build | Arithmetic propagation proofs complete |
| §1.4 | Build | Transcendental propagation proofs complete |
| §1.5 | Build | Series truncation classification complete |
| §1.6 | Build | Iterative convergence proofs complete |
| §1.7 | Build | Subtractive cancellation analysis complete |
| §1.8 | Build | Error combination rules complete |
| §1.9 | Build | Reliable digits and comparison complete |
