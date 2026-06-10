# Draft Plan: Chapter 4 — Approximation Theory and Fast Convergence

## Objectives

- Develop the complete toolkit of approximation methods for evaluating functions to a caller-specified tolerance with rigorous error bounds.
- Resolve forward references from Chapter 1: the continued-fraction / Padé methods promised in §1.5 [P.1.2] and the tolerance-parameter method selection promised after Theorem 1.6.3.
- Establish the cross-cutting principle of exact rational arithmetic for series coefficients (defer rounding to the final multiplication by a power of the variable).
- Provide a decision framework (§4.11) that maps function properties to the optimal method at any precision level, unified under the tolerance parameter $\tau$ from Chapter 3.
- Verify that every method is matched-pair compatible with SR3 at $\tau_{\mathrm{standard}}$.

---

## Section Structure

### §4.1 Introduction
- No new definitions or theorems.
- Narrative: roles of Ch 1 remainder bounds and Ch 3 tolerance parameter; overview of §4.2–§4.11 arc; exact-rational-coefficient principle; resolution of Ch 1 forward references [P.1.2].

---

### §4.2 Taylor Series: From Integer to Fractional Order

#### §4.2.1 Classical Taylor Series
- **Theorem 4.2.1** — Taylor's theorem with Lagrange remainder. Proof via auxiliary function and Rolle's theorem.
- **Theorem 4.2.2** — Taylor's theorem with integral remainder. Proof by iterated integration by parts; structural remark that integral form = $J^{N+1}[f^{(N+1)}]$.
- **Corollary 4.2.1** — Connection to Ch 1: geometric and alternating-series remainder bounds (Theorems 1.5.2–1.5.3) as special cases.
- **Remark** — Three limitations of classical Taylor: (1) convergence disk, (2) slow convergence near boundary, (3) integer-power basis. Motivates §4.3–§4.4 (lims 1–2) and fractional generalization (lim 3).

#### §4.2.2 Generalizing the Power Rule
- **Definition 4.2.1** — Fractional derivative of a monomial: $D^\alpha x^n = \frac{\Gamma(n+1)}{\Gamma(n-\alpha+1)} x^{n-\alpha}$.
- **Example 4.2.1** — Half-derivative of $x$: $D^{1/2} x = \frac{2}{\sqrt{\pi}}\sqrt{x}$; half-derivative of a constant.
- **Example 4.2.3** — Table of half-derivatives of monomials $1, x, x^2, x^3$. Structural note: half-derivative maps polynomials to a different function space.

#### §4.2.3 Fractional Integration and Differentiation
- **Definition 4.2.2** — Riemann–Liouville fractional integral $J^\alpha f$.
- **Definition 4.2.3** — Riemann–Liouville fractional derivative $D^\alpha f$ (integrate up by $n-\alpha$, differentiate down by $n$).
- **Definition 4.2.4** — Caputo fractional derivative ${}^C\!D^\alpha f$ (differentiate first, then integrate). Preferred in physical applications; preserves standard initial conditions.
- **Remark** — Caputo vs. Riemann–Liouville: physical interpretation of initial conditions.
- **Remark** — Integral remainder of Theorem 4.2.2 is structurally identical to $J^{N+1}[f^{(N+1)}]$.

#### §4.2.4 The Generalized Taylor Series
- **Theorem 4.2.3** — Fractional Mean Value Theorem. Proof via $J^\alpha \,{}^C\!D^\alpha$ composition and integral MVT.
- **Theorem 4.2.4** — Generalized Taylor series (Caputo form): $f(x) = \sum_{k=0}^\infty \frac{(x-a)^{k\alpha}}{\Gamma(k\alpha+1)}({}^C\!D^{k\alpha}f)(a)$. Reduces to classical at $\alpha=1$.
- **Theorem 4.2.5** — Generalized remainder bound (fractional Lagrange form). Reduces to (4.2) at $\alpha=1$.
- **Remark** — Remainder (4.12) has the same role as integer remainder; $\tau$ governs truncation identically.

#### §4.2.5 Convergence: Fractional vs Classical
- **Theorem 4.2.6** — Convergence rate comparison: fractional series needs $\lceil 1/\alpha \rceil$ times as many terms for smooth functions, but converges for non-integer power-law $f(x) \sim C(x-a)^\beta$ where classical fails.
- **Example 4.2.2** — Singularity resolution: $\sqrt{x}$ has no classical Taylor series at 0; $\alpha=1/2$ fractional series terminates in one term.
- **Remark** — Mittag-Leffler function $E_\alpha(z)$: fractional generalization of $e^z$.
- **Remark** — Fractional Newton's method: half-Newton iteration for zeros of fractional order; classical methods suffice for all current pipeline targets (Ch 9, 14, 32); fractional regime relevant at high-eccentricity osculating-to-mean inversion.

---

### §4.3 Padé Approximants
- **Definition 4.3.1** — Padé approximant $[p/q]_f$: rational $P_p/Q_q$ with $Q_q(0)=1$ matching $p+q+1$ Taylor coefficients to order $O(x^{p+q+1})$.
- **Theorem 4.3.1** — Existence and uniqueness. Proof via Toeplitz linear system for $Q_q$ coefficients.
- **Theorem 4.3.2** — Padé error order: $|f(x) - [p/q]_f(x)| = O(|x|^{p+q+1})$.
- **Proposition 4.3.1** — Diagonal Padé error bound: geometric convergence ratio $(|x|/\rho)^2$ vs. Taylor's $(|x|/\rho)$. Proof sketch via Montessus de Ballore theorem.
- **Example 4.3.2** — Numerical comparison at $|x|/\rho = 0.9$: Taylor ratio $0.9$, Padé ratio $0.81$.
- **Example 4.3.1** — $[1/1]$ Padé for $e^x$: max error $0.056$ on $[-1,1]$ vs. $0.218$ for degree-2 Taylor.
- **Remark** — Connection to continued fractions: diagonal Padé sequence = even convergents of CF (→ §4.4, Proposition 4.4.1).

---

### §4.4 Continued Fractions
- **Definition 4.4.1** — Generalized continued fraction; $n$-th convergent $C_n = A_n/B_n$; partial numerators $a_k$, partial denominators $b_k$.
- **Theorem 4.4.1** — Wallis–Euler recurrence for $A_n$, $B_n$. Proof by induction.
- **Theorem 4.4.2** — Śleszyński–Pringsheim convergence criterion: $|b_n| \geq |a_n| + 1$ implies convergence. Proof sketch via Stern–Stolz theorem.
- **Example 4.4.1** — $\arctan$ continued fraction: Śleszyński–Pringsheim fails for large $|x|$ at small $n$; convergence for all real $x$ follows from Seidel–Stern theorem (positive $a_n$).
- **Theorem 4.4.3** — Truncation error bound for alternating convergents: $|f - C_n| \leq |C_n - C_{n-1}|$. Directly usable as stopping criterion $|C_n - C_{n-1}| < \tau$.
- **Proposition 4.4.1** — Euler's method: power-series-to-CF conversion via QD algorithm; even convergents = diagonal Padé.
- **Proposition 4.4.2** — Gauss continued fraction: ratio of contiguous ${}_2F_1$; recovers $\arctan$, $\log(1+x)$.
- **Remark** — Modified Lentz algorithm: avoids overflow/underflow in Wallis–Euler recurrence (→ [P.4.2]).

---

### §4.5 Chebyshev Polynomials and Economization
- **Definition 4.5.1** — Chebyshev polynomial of first kind: $T_n(x) = \cos(n\arccos x)$; three-term recurrence.
- **Theorem 4.5.1** — Minimax property: $T_n/2^{n-1}$ is the monic degree-$n$ polynomial with smallest sup-norm on $[-1,1]$. Proof by contradiction using equioscillation at $n+1$ nodes.
- **Theorem 4.5.2** — Chebyshev coefficient decay: $|c_n| \leq 2M/\rho^n$ for $f$ analytic inside ellipse $\mathcal{E}_\rho$.
- **Corollary 4.5.1** — Chebyshev truncation error bound: $|f - \sum_{k=0}^N c_k T_k| \leq 2M\rho^{-N}/(\rho-1)$.
- **Proposition 4.5.1** — Clenshaw algorithm: $O(N)$ backward recurrence, numerically stable, error $O(N\,\epsilon_{\mathrm{mach}}\max|c_k|)$. Proof via Chebyshev-of-second-kind representation.
- **Remark** — Chebyshev economization: drop highest-degree terms of Taylor-to-Chebyshev reexpansion; relevant to Ch 15 Hansen coefficient polynomials.

---

### §4.6 Convergence Acceleration
- **Definition 4.6.1** — Sequence transformation / convergence acceleration: $(\hat{s}_n - s)/(s_n - s) \to 0$.
- **Theorem 4.6.1** — Aitken's $\Delta^2$ method: linear convergence $s_n - s \approx C\lambda^n$ → accelerated $\hat{s}_n - s = O(\lambda^{2n})$. Full proof by substitution.
- **Theorem 4.6.2** — Euler transform for alternating series: $\sum (-1)^k a_k \to \sum (-1)^k \Delta^k a_0 / 2^{k+1}$.
- **Remark** — Euler transform applicability: effective for polynomial decrease ($a_k \sim 1/k$); negligible improvement for geodetic series ($a_k \sim e'^{2k}$, ratio $0.007$).
- **Remark** — Richardson extrapolation: eliminates leading $h^{p_1}$ error term; foundation of Romberg integration.

---

### §4.7 Iterative Methods: Convergence Orders and Stability

#### §4.7.1 General Theory
- **Theorem 4.7.1** — Newton's method (recap from Ch 1 §1.6): quadratic convergence to simple root.
- **Theorem 4.7.2** — Newton convergence basin: multiplier $N'(x) = f f'' / f'^2$; fast convergence radius $r \geq f'_{\min}/|f''|_{\max}$.
- **Theorem 4.7.3** — Halley's method: cubic convergence. Proof sketch via $[1/1]$ Padé approximation to $f$.
- **Theorem 4.7.4** — Halley stability criterion: denominator $\Delta_H = 2f'^2 - ff''$ must be nonzero; instability when $|ff''| \approx 2f'^2$.
- **Theorem 4.7.5** — Householder's method of order $d$: convergence order $d+1$ via derivatives of $1/f$. Newton ($d=1$) and Halley ($d=2$) as special cases.
- **Corollary 4.7.1** — Iteration count estimate: $k = \lceil \log(d_{\mathrm{target}}/d_0) / \log p \rceil$. Table: Newton/Halley/Householder-3 from 3 to 16/50/100 digits.
- **Proposition 4.7.1** — Newton for inverse functions: apply Newton to $g(x) = f(x) - y$.

#### §4.7.2 Stability of Orbital Iterative Equations
- **Proposition 4.7.2** — Kepler equation stability: (i) global convergence (f' > 0 always); (ii) fast convergence radius $(1-e)/e$ with table by eccentricity; (iii) Halley denominator bound $\sim 2(1-e)^2$.
- **Proposition 4.7.3** — Cube root stability: basin $[0.737\,a^*, \infty)$; starter $a_0 = b$ adequate for $b > 1$; ~4–5 Newton iterations for double precision.
- **Proposition 4.7.4** — Modified Kepler equation (Lyddane form) stability: global convergence since $f'(x) \leq e-1 < 0$; eccentricity limits mirror standard Kepler.
- **Proposition 4.7.5** — $J_2$ fixed-point iteration stability: contraction with $L \approx 0.003$; 7 iterations suffice for double precision; 3–5 in practice.

#### §4.7.3 Fractional-Order Zeros
- **Remark** — Newton degrades at fractional-order zeros ($f \sim C(x-x^*)^\beta$, $\beta \notin \mathbb{N}$): fractional Newton (§4.2.5) for those cases; all current pipeline targets have simple zeros (Props. 4.7.2–4.7.5).

---

### §4.8 Argument Reduction
- **Theorem 4.8.1** — Additive argument reduction: periodic $f$ with period $P$; precision loss bound (4.35) involving $|n|\cdot\delta_p(P)$ and $\epsilon_{\mathrm{mach}}\cdot|nP|$.
- **Example 4.8.1** — Precision loss at large $x$: $x = 10^8$ rad loses $\sim 7.5$ digits; $x = 10^{15}$ rad loses all 16 digits. Motivates Payne-Hanek (Ch 29).
- **Theorem 4.8.2** — Multiplicative argument reduction: $\exp(x) = (\exp(x/2^k))^{2^k}$; half-angle formulas for trig.
- **Remark** — Payne-Hanek algorithm: pre-computed high-precision $1/P$ extracts fractional part without full quotient.
- **Remark** — Reconstruction error: each squaring doubles relative error; total $O(\log|x|\cdot\epsilon_{\mathrm{mach}})$ after $k$ squarings.

---

### §4.9 Polynomial Evaluation
- **Definition 4.9.1** — Horner form: $p(x) = c_0 + x(c_1 + x(c_2 + \cdots + x\cdot c_n))$; $n$ multiplications and $n$ additions.
- **Theorem 4.9.1** — Horner rounding bound: $|\hat{p}(x) - p(x)| \leq (2n+1)\epsilon_{\mathrm{mach}}\,\tilde{p}(|x|)$. Full proof by unrolling recurrence with product error bound.
- **Remark** — Naive evaluation is $O(n^2\,\epsilon_{\mathrm{mach}})$ vs. Horner's $O(n\,\epsilon_{\mathrm{mach}})$.
- **Theorem 4.9.2** — Error propagation with uncertain argument: total error = $|p'(x)|\cdot\delta_x + (2n+1)\epsilon_{\mathrm{mach}}\,\tilde{p}(|x|)$ (cf. Ch 1, Theorem 1.4.1).
- **Remark** — Compensated Horner: Dekker splitting for large condition number $\tilde{p}(|x|)/|p(x)|$; effectively doubles working precision. Addresses [P.4.4].

---

### §4.10 The Generalized Binomial Series
- **Definition 4.10.1** — Generalized binomial coefficient: $\binom{\alpha}{k} = (\alpha)_k^{(-)} / k! = \alpha(\alpha-1)\cdots(\alpha-k+1)/k!$.
- **Theorem 4.10.1** — Convergence of binomial series for $|x| < 1$: $(1+x)^\alpha = \sum_k \binom{\alpha}{k} x^k$; geometric remainder bound $|R_N| \leq |\binom{\alpha}{N+1}|\,|x|^{N+1}/(1-|x|)$. Proof via Taylor coefficients and geometric tail bound (Ch 1, Theorem 1.5.3).
- **Theorem 4.10.2** — Exact rational coefficients: when $\alpha = p/q$, each $\binom{p/q}{k}$ is exactly rational per (4.43); only rounding is in the $n$ multiplications by $x$.
- **Example 4.10.1** — $\alpha = 1/2$ (square root): first four exact rational coefficients $1, 1/2, -1/8, 1/16$.
- **Example 4.10.2** — $\alpha = -2$ (relevant to $(1 + e'^2\cos^2\phi)^{-2}$ in Ch 14): $\binom{-2}{k} = (-1)^k(k+1)$ are exact integers.
- **Remark** — Pochhammer connection: rising vs. falling factorials; binomial series as ${}_1F_0(\alpha;-;-x)$; connection to Gauss CF (§4.4, Prop. 4.4.2).
- **Remark** — Horner nesting for binomial series: combines Theorem 4.9.1 rounding advantage with Theorem 4.10.2 coefficient exactness.

---

### §4.11 Method Selection and the Tolerance Parameter
- **Table 4.11.1** — Method selection by function property (8 rows: Chebyshev, CF, Padé, Horner, fractional Taylor, Euler, Newton/Halley, Taylor).
- **Table 4.11.2** — Precision-level strategy (3 rows: ~16 digits, ~50 digits, ~100+ digits).
- **Proposition 4.11.1** — Matched-pair compatibility: every method in Tables 4.11.1–4.11.2 at $\tau = \tau_{\mathrm{standard}}$ satisfies the Ch 3 compatibility criterion (Def. 3.7.1). Proof via triangle inequality over $\tau_{\mathrm{standard}} + \delta_p^{\mathrm{SR3}} \ll \delta_a(\mathcal{M})$.
- **Corollary 4.11.1** — Tolerance parameter unification: method selection is internal; external interface is identical regardless of method.

---

## Cross-References

| Reference | Direction | Target |
|-----------|-----------|--------|
**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Thms 1.5.2–1.5.3 | §4.2, §4.10 | Series remainder bounds |
| Ch 1, Thms 1.6.1–1.6.3 | §4.7 | Iterative convergence bounds |
| Ch 1, Thm 1.4.1 | §4.9 | Derivative sensitivity for Horner evaluation |
| Ch 1, §1.7 | §4.8–4.9 | Subtractive cancellation analysis |
| Ch 3, Def 3.4.1 | §4.2, §4.4, §4.11 | Tolerance parameter $\tau$ |
| Ch 3, Def 3.4.2 | §4.11 | Standard tolerance $\tau_{\mathrm{standard}}$ |
| Ch 3, Def 3.7.1, Cor 3.7.1 | §4.11 | Matched-pair compatibility |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 5 | §4.4, §4.7, §4.9, §4.10 | CF convergence, Newton iteration, Horner form, binomial series |
| Ch 9 | §4.7 | Iterative solver framework for Kepler equation |
| Ch 14 | §4.10 | Binomial series for geodetic integrals |
| Ch 15 | §4.5 | Chebyshev economization for Hansen coefficients |
| Ch 29 | §4.8 | Payne-Hanek argument reduction for sidereal time |
| Ch 32 | §4.7 | Iterative solver for element recovery |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.4.1] | P | §4.3 | Padé near poles: division by small $Q_q(x)$ amplifies rounding |
| [P.4.2] | P | §4.4 | Wallis-Euler recurrence overflow/underflow; use modified Lentz |
| [P.4.3] | P | §4.5 | DCT rounding $O(\sqrt{n}\,\epsilon_{\mathrm{mach}})$; negligible for $n \leq 10$ |
| [P.4.4] | P | §4.9 | Horner with alternating signs near $|x| \approx 1$: large condition number |
| [P.4.5] | P | §4.8 | Additive argument reduction loses $\log_{10}(n)$ digits for large $n$ |
| [P.4.6] | P | §4.10 | Binomial series near $|x| = 1$: slow convergence, accumulating rounding |
| [P.4.7] | P | §4.2.3 | Riemann-Liouville integral may require numerical quadrature |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 10 |
| Theorems | 26 |
| Lemmas | 0 |
| Corollaries | 4 |
| Propositions | 11 |
| Examples | 9 |
| Error Notes | 7 |
| Equations | ~43 |
| Sections | 11 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §4.1 | Build | Introduction and forward references |
| §4.2 | Build | Taylor series and fractional calculus |
| §4.3 | Build | Pade approximants |
| §4.4 | Build | Continued fractions |
| §4.5 | Build | Chebyshev polynomials |
| §4.6 | Build | Convergence acceleration |
| §4.7 | Build | Iterative solvers |
| §4.8 | Build | Argument reduction |
| §4.9 | Build | Horner evaluation |
| §4.10 | Build | Generalized binomial series |
| §4.11 | Build | Decision framework |
