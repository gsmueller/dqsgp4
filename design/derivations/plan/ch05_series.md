# Draft Plan: Chapter 5 — Series Evaluation and Error Control

**Part I: Mathematical Foundations**

---

## Objectives

- Formalize the three generic series evaluators (alternating, geometric-tail, Horner) that cover all series arising in satellite orbit propagation.
- Integrate each evaluator with the TrackedValue error-tracking framework (precision vs. accuracy, measurement vs. model error).
- Apply the generic evaluators to the specific series of the project: the $q_0$, $q_0'$, and $U_0$ arctan-derived auxiliary series (§5.5) and the geodetic binomial family (§5.6).
- Catalog further orbit-propagation series (Kepler/Bessel, Hansen, drag, ephemeris) by evaluator type (§5.7).
- Connect evaluators to the tolerance parameter and matched-pair principle (§5.8).
- Provide numerical verification against 60-digit reference values (§5.9).

---

## Section Structure

### §5.1 Introduction

No formal items. Prose motivating the three-evaluator taxonomy, the precision-vs-accuracy classification of truncation error, and the TrackedValue output contract. Roadmap of §§5.2–5.8.

---

### §5.2 Alternating Series Evaluation

- **Definition 5.2.1** (Alternating series evaluator) — Term function $t:\mathbb{N}\to\mathcal{V}$, start index, tolerance $\tau$, $N_{\max}$; returns $\mathcal{V}_S$ with Leibniz truncation bound in $\delta_p$.
- **Theorem 5.2.1** (Leibniz truncation bound for tracked evaluation) — Bound (5.1): $|S - \mathcal{V}_S.v| \leq |t(N).v|$; total precision (5.2): arithmetic rounding plus truncation bound.
- **Corollary 5.2.1** (Error category preservation) — Truncation bound added exclusively to $\delta_p$; $\sigma_m$ and $\delta_a$ inherit from input terms via Ch 1 addition rules.
- **Proposition 5.2.1** (Term count estimate) — $N \approx k_0 + \lceil \ln(\tau/|t(k_0).v|)/\ln r \rceil$ (Eq. 5.3) for approximately constant convergence ratio.
- **Example 5.2.1** — Term count table for $2q_0$ series at double, extended, and 50-digit precision.
- **Remark** (Monotonicity verification) — When $|c_{k+1}/c_k| \cdot |x| < 1$, monotonicity holds automatically; runtime check recommended. Cross-reference [P.5.1].

---

### §5.3 Geometric-Tail Series Evaluation

- **Definition 5.3.1** (Geometric-tail series evaluator) — Adds convergence ratio bound $r\in(0,1)$; stopping criterion $|t(N).v| \cdot r/(1-r) < \tau$; tail bound added to $\delta_p$.
- **Theorem 5.3.1** (Geometric tail bound for tracked evaluation) — Bound (5.4): $|S - \mathcal{V}_S.v| \leq |t(N).v| \cdot r/(1-r)$; precision (5.5).
- **Corollary 5.3.1** (Comparison with Leibniz) — Geometric bound tighter when $r < 1/2$; Leibniz tighter when $r > 1/2$.
- **Remark** — For WGS84 geodetic series ($r \approx e^2 \approx 0.0067$), geometric bound is ~150 times tighter than Leibniz.
- **Theorem 5.3.2** (Ratio test convergence detection) — Adaptive ratio estimate $r_k = |t(k).v/t(k-1).v|$ (Eq. 5.6); rigorous when $|c_{k+1}/c_k|$ is eventually monotone (binomial, Bessel, hypergeometric).
- **Remark** (Adaptive versus fixed ratio) — Fixed ratio preferred when $r$ known a priori; adaptive when argument-dependent. Cross-reference [P.5.4].

---

### §5.4 Horner Evaluation of Truncated Series

- **Theorem 5.4.1** (Horner evaluation with tracked errors) — $N$ multiplications and $N$ additions; total precision bound (5.8): $N \epsilon_{\mathrm{mach}} \tilde{p}(|x|)$ plus coefficient error term; uses condition polynomial $\tilde{p}$.
- **Corollary 5.4.1** (Exact rational coefficient advantage) — When $\delta_p(c_k) = 0$ for all coefficients, simplified bound (5.9): rounding from multiplications by $x$ only.
- **Proposition 5.4.1** (Combined truncation and evaluation procedure) — Three-step protocol: (1) determine $N$ via Theorem 5.2.1 or 5.3.1, (2) evaluate $S_N$ via Horner, (3) add truncation bound to $\delta_p$.
- **Remark** (Well-conditioning for small arguments) — For geodetic series with $x = e'^2$ or $e^2 \approx 0.0067$, condition number $\approx 1$; no cancellation. Cross-reference Ch 4 [P.4.4].

---

### §5.5 Application: The $q_0$ and $q_0'$ Auxiliary Series

#### §5.5.1 The $q_0$ Series

- **Definition 5.5.1** (The quantity $q_0$) — Closed form (5.10): $q_0 = \tfrac{1}{2}[(1+3/e'^2)\arctan e' - 3/e']$; role in normal potential $U_0$ and $J_{2n}$ derivation (Ch 14).
- **Theorem 5.5.1** (Series expansion for $q_0$) — Series (5.11): $2q_0 = \sum_{n=1}^\infty (-1)^{n+1} \frac{4n}{(2n+1)(2n+3)} e'^{2n+1}$; alternating, monotone, convergence ratio $r = e'^2$. Full proof from $\arctan$ Maclaurin expansion.
- **Corollary 5.5.1** (Truncation bound for $q_0$) — Leibniz bound (5.20): $|R_N| \leq \frac{4N}{(2N+1)(2N+3)} e'^{2N+1}$; table of bounds for $N = 3,5,7,10$ at WGS84.
- **Corollary 5.5.2** (Cancellation avoidance) — Closed form loses ~5.4 significant digits (WGS84); series avoids cancellation algebraically via Eq. (5.17). Cross-reference [P.5.2].

#### §5.5.2 The $q_0'$ Series

- **Definition 5.5.2** (The quantity $q_0'$) — Closed form (5.22): $q_0' = 3(1+1/e'^2)(1-\arctan(e')/e') - 1$; role in Somigliana formula and $J_2$–flattening relation (Ch 14, §14.7).
- **Theorem 5.5.2** (Series expansion for $q_0'$) — Series (5.23): $q_0' = \sum_{n=1}^\infty (-1)^{n+1} \frac{6}{(2n+1)(2n+3)} e'^{2n}$; proof deferred to Ch 14, Theorem 14.4.3; convergence ratio $\to e'^2$.
- **Corollary 5.5.3** (Structural comparison) — Side-by-side table comparing $2q_0$ and $q_0'$: odd vs. even powers, linear vs. constant numerator, shared denominator and convergence ratio.

#### §5.5.3 The $U_0$ Auxiliary Series

- **Theorem 5.5.3** (Series expansion for $\arctan(e')/e'$) — Series (5.28): $\arctan(e')/e' = \sum_{n=0}^\infty (-1)^n e'^{2n}/(2n+1)$; alternating, monotone, ratio $r = e'^2$.
- **Remark** — Role in $U_0$ computation; simpler coefficients ($1/(2n+1)$ with no $(2n+3)$ denominator) compared to $q_0$ and $q_0'$.

#### §5.5.4 Horner Nesting

- **Proposition 5.5.1** (Horner form for $q_0$) — Factored form (5.25) with $u = -e'^2$; first four exact rational coefficients (5.26): $4/15, 8/35, 4/21, 16/99$. Analogous Horner form for $q_0'$.

#### §5.5.5 Continued Fraction Alternative

- **Proposition 5.5.2** (Continued fraction for $q_0$) — Gauss CF for $\arctan$ (Eq. 5.27); three reasons series (5.11) is preferred over CF: error bound transparency, cancellation avoidance, performance at small arguments.

#### §5.5.6 Implementation Correspondence

Prose (no formal items). Maps series (5.11) to `EquipotentialEllipsoid` implementation: `ratio<T>`, `alternating_series<T>`, power accumulation without `pow`.

---

### §5.6 Application: Geodetic Binomial Series

#### §5.6.1 The Geodetic Binomial Family

- **Definition 5.6.1** (Geodetic binomial series) — Canonical form (5.29): $(1+e'^2\cos^2\Phi)^\alpha$; sign convention (alternating signs from $\binom{\alpha}{k}$ itself, no extra $(-1)^k$); table of five exponents and their applications.
- **Theorem 5.6.1** (Convergence for geodetic parameters) — Absolute convergence for $e'^2 < 1$; remainder bound (5.30) via Ch 4, Theorem 4.10.1 and Ch 1, Theorem 1.5.3.
- **Table 5.6.1** — Remainder bounds at $N=4$ for WGS84, four exponents ($\alpha = -1/2, -3/2, -2, -5/2$).

#### §5.6.2 The Three-Stage Coefficient Derivation

- **Theorem 5.6.2** (Three-stage coefficient derivation) — Stage 1: binomial expansion yielding $\binom{\alpha}{k}$; Stage 2: Wallis integration ($W_{2k}$ for $w=1$, $(2k)!!/(2k+1)!!$ for $w=\cos\Phi$); Stage 3: series algebra (Cauchy product, Newton iteration for square root). Rationality of all coefficients at each stage.
- **Corollary 5.6.1** (Compute, don't look up) — Coefficients must be computed from the three-stage formula for arbitrary-precision use; storing double-precision literals limits result to ~16 digits.
- **Example 5.6.1** — Meridian arc $Q$: $\alpha = -3/2$, $w=1$; Stage 1–2 combined coefficient (5.31); first five values.
- **Example 5.6.2** — Equal-area radius $R_2$: $\alpha = -2$, $w=\cos\Phi$, then square root (Stage 3 Newton iteration); result $R_2/c = 1 - 2e'^2/3 + 26e'^4/45 - \cdots$.

#### §5.6.3 Integrated Series Truncation Bound

- **Theorem 5.6.3** (Truncation bound for integrated binomial series) — Bound (5.32): Wallis factor $\leq 1$ tightens the raw binomial bound (5.30); separate cases for $w=1$ and $w=\cos\Phi$.

#### §5.6.4 Series Algebra Error Bounds

- **Lemma 5.6.1** (Series square root via Newton iteration) — Quadratic convergence bound (5.33): $|s_j - \sqrt{I}| \leq |U|^{2^j}/2^{j+1}$; 4 iterations give error $< 10^{-40}$ at WGS84.
- **Lemma 5.6.2** (Series quotient error bound) — Precision bound (5.34) via quotient rule; input precision errors plus single division rounding.

---

### §5.7 Application: Further Series in Orbit Propagation

#### §5.7.1 Kepler's Equation: Bessel Function Series (Ch 9)

- **Proposition 5.7.1** (Bessel series for Kepler's equation) — Series (5.35): geometric tail with $r \approx e/2$; Newton iteration (Ch 4, Theorem 4.7.1) preferred computationally; Bessel series retained for theory and equation-of-center derivation (Ch 25).

#### §5.7.2 Hansen Coefficients (Ch 15)

- **Proposition 5.7.2** (Hansen coefficient series) — Power series in $e$ (5.36) with exact rational Newcomb-operator coefficients; geometric tail, $r = e$; $N \leq 10$ terms at double for LEO.

#### §5.7.3 Drag Coefficient Series (Ch 22)

- **Proposition 5.7.3** (Orbit-averaged drag series) — Non-alternating, positive, geometric-tail with $r = e$; geometric-tail evaluator (§5.3) applies.

#### §5.7.4 Trigonometric Ephemeris Truncation (Ch 25–26)

- **Proposition 5.7.4** (Ephemeris truncation as accuracy error) — Sole exception to the chapter pattern: each omitted term is a distinct physical perturbation; truncation is $\delta_a$, not $\delta_p$. Cross-reference [A.5.1].

---

### §5.8 The Series Evaluation Interface

#### §5.8.1 Interface Design

- **Definition 5.8.1** (Series evaluator interface) — Common input/output contract: term function, start index, tolerance $\tau$, $N_{\max}$; output TrackedValue with $\delta_p$ covering rounding and truncation.
- **Remark** (Init-time constants vs per-step evaluations) — Two evaluation lifetimes: (1) init-time constants ($q_0$, $q_0'$, $U_0$, binomial coefficients, $J_{2n}$) — evaluated once, stored as members; (2) per-step evaluations (Kepler, Hansen, drag) — callable $\lambda(e)$ capturing $\tau$ in closure.
- **Remark** (Bound selection for the library) — Leibniz overestimates by ~170×; geometric tail overestimates by 10–20%; library should accept optional $r$ for tighter bound.
- **Theorem 5.8.1** (Tolerance parameter integration) — Bound (5.37): $|g(\tau).v - f_{\mathrm{exact}}| \leq \tau + \delta_p^{\mathrm{arith}}$; matched-pair compatibility at $\tau_{\mathrm{standard}}$ (5.38): all precision-level gaps are $\ll \delta_a(\mathcal{M}) \approx 1$ km.

#### §5.8.2 Method Selection by Series

- **Table 5.8.1** (Method comparison for series with multiple evaluation paths) — $q_0$, $q_0'$, geodetic binomial (four exponents), $R_2$ sqrt, $\bar{\gamma}$ quotient, Kepler, Hansen, drag, ephemeris; preferred method checkmarked.

#### §5.8.3 Term Counts by Tolerance Level

- **Proposition 5.8.1** (Term count as a function of tolerance) — Applies Eq. (5.3) to five tolerance levels: $10^{-3}$ (accuracy floor), $10^{-12}$ (implementation), $10^{-16}$ (double), $10^{-34}$ (quad), $10^{-50}$ (50-digit).
- **Table 5.8.2** (Term counts by tolerance regime, WGS84) — Six series × five tolerance columns.
- **Table 5.8.3** (SR3 truncation orders and matched tolerance) — Seven SR3 computations: truncation order, magnitude of first omitted term, $\tau_{\mathrm{match}}$.
- **Remark** — SR3 secular rates (Hoots and Roehrich 1980, lines 1546–1553 of `SGP4.cpp`) are closed-form polynomials; generic evaluator reproduces them at $\tau = \tau_{\mathrm{match}}$.

#### §5.8.4 Accuracy-Level Stopping

- **Proposition 5.8.2** (Accuracy-limited evaluation) — At $\tau = \delta_a(\mathcal{M})$, single term suffices for all geodetic series; $\delta_p$ records the truncation bound regardless of tolerance level.
- **Corollary 5.8.1** (Single algorithm, two regimes) — Same evaluator serves matched-pair regime and precision regime; only $\tau$ value changes at initialization.
- **Corollary 5.8.2** (Automatic precision scaling) — Widening $T$ changes $\epsilon_{\mathrm{mach}}$ and therefore the precision-regime $\tau$; evaluators automatically compute additional terms.

---

### §5.9 Numerical Verification

No new formal items. Reference values computed at 60-digit precision; implementation via `EquipotentialEllipsoid<T>` with WGS84 parameters.

#### §5.9.1 The $q_0$ Series Term by Term

- **Table 5.9.1** — First 10 terms: term value, partial sum, Leibniz bound, ratio $|t_{n+1}/t_n|$; converged value $q_0 = 7.33462578708345207 \times 10^{-5}$ (25 terms, 60-digit).

#### §5.9.2 Convergence Rate Verification

- **Table 5.9.2** — Actual ratios vs. Eq. (5.19) prediction vs. asymptote $e'^2$; agreement to all displayed digits.

#### §5.9.3 Cancellation in the Closed Form

- **Table 5.9.3** — Closed-form $q_0$: operands, difference; double-precision result loses ~5 digits; series agrees to full double precision.

#### §5.9.4 Error Bound Tightness

- **Table 5.9.4** — Leibniz bound vs. geometric tail bound vs. actual remainder for $N = 1,2,3,5,7,10$; Leibniz/actual $\approx 170$; geometric/actual $\approx 1.1$–$1.2$.
- **Remark** — Leibniz conservatism is not a defect; geometric bound preferred when $r$ known.

#### §5.9.5 TrackedValue Error Budget

- **Table 5.9.5** — Arithmetic $\delta_p^{\mathrm{arith}}$ vs. truncation bound vs. total $\delta_p$ at $n = 1,3,5,7,8$ for $\tau = 10^{-16}$; crossover from truncation-dominated to arithmetic-dominated at $n = 8$.
- **Remark** (Verification protocol) — `T = double` reproduces Table 5.9.1 to 15 sig figs; `TrackedValue<T>.errors.precision` carries Table 5.9.5 decomposition.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Defs 1.2.2–1.2.4 | Notation | $\delta_p$, $\delta_a$, $\mathcal{V}$ definitions |
| Ch 1, Thms 1.3.1–1.3.2 | §5.2–5.4 | Error propagation through arithmetic |
| Ch 1, Thms 1.5.1–1.5.3 | §5.2–5.3 | Truncation classification; Leibniz and geometric bounds |
| Ch 1, §1.7 | §5.5 | Catastrophic cancellation analysis for [P.5.2] |
| Ch 3, Defs 3.4.1–3.4.2 | §5.8 | Tolerance parameter $\tau$ and $\tau_{\mathrm{standard}}$ |
| Ch 3, Cor 3.7.1 | §5.8 | Matched-pair compatibility criterion |
| Ch 4, §4.4 | §5.5.5 | Gauss CF for $\arctan$; CF convergence bound |
| Ch 4, §4.7 | §5.7.1 | Newton iteration framework |
| Ch 4, §4.9 | §5.4 | Condition polynomial for Horner |
| Ch 4, §4.10 | §5.3, §5.6 | Generalized binomial series; exact rational coefficients |
| Ch 4, §4.11 | §5.8 | Decision framework; matched-pair verification |
| Ch 6 | §5.6.2 | Wallis integrals $W_n$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 9 | §5.7.1 | Kepler equation series and iteration |
| Ch 14 | §5.5–5.6 | $q_0$, $q_0'$ series; geodetic binomial series |
| Ch 15 | §5.7.2 | Hansen coefficient series |
| Ch 22 | §5.7.3 | Drag coefficient series |
| Ch 25–26 | §5.7.4 | Ephemeris truncation analysis |
| App C | §5.5.6 | Code mapping for series evaluators |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.5.1] | P | §5.2 | Alternating series monotonicity failure: Leibniz bound invalid until monotonicity established |
| [P.5.2] | P | §5.5.1 | Catastrophic cancellation in closed-form $q_0$ (~5.4 digits lost); use series form |
| [P.5.3] | P | §5.4 | Rounding accumulation $O(N \epsilon_{\mathrm{mach}} |S_N|)$; Kahan summation for $N > 1000$ |
| [P.5.4] | P | §5.3 | Geometric tail bound non-conservative if $r_k$ underestimated; use fixed $r$ |
| [M.5.1] | M | §5.8 | Measurement error propagation through series; condition number $\approx 1$ for geodetic |
| [A.5.1] | A | §5.7.4 | Ephemeris truncation is $\delta_a$, not $\delta_p$; physical perturbation omitted |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 7 |
| Theorems | 11 |
| Lemmas | 2 |
| Corollaries | 9 |
| Propositions | 10 |
| Examples | 3 |
| Error Notes | 6 |
| Equations | ~38 |
| Sections | 9 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §5.1 | Build | Introduction and forward references |
| §5.2 | Build | Alternating series evaluator |
| §5.3 | Build | Geometric-tail evaluator |
| §5.4 | Build | Horner evaluator |
| §5.5 | Build | Equipotential series ($q_0$, $q_0'$, arctan) |
| §5.6 | Build | Geodetic binomial series |
| §5.7 | Build | Downstream applications |
| §5.8 | Build | Tolerance integration and matched-pair |
| §5.9 | Build | Verification tables |
