# Draft Plan: Chapter 6 — Combinatorics and Special Functions

## Objectives

1. Define factorial, double factorial, and their extensions to non-integer arguments via the Gamma function.
2. Establish generalized binomial coefficients C(α, k) for real α and integer k ≥ 0, including sign patterns and recurrences needed by Ch 5.
3. Define Pochhammer symbols (rising and falling factorials) and unify them with binomial coefficients.
4. Derive the Wallis integral sequence W_n = ∫₀^{π/2} cos^n φ dφ, establish the two-term recurrence, and give closed forms via double factorial.
5. Connect all to the Gamma function, enabling extension to non-integer orders.
6. All definitions are general-purpose (not application-specific).

**Implements:** `factorial.h`, `wallis.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $n!$ | Factorial: $n! = n \cdot (n-1)!$, $0! = 1$ | §6.2 |
| $n!!$ | Double factorial: $n!! = n \cdot (n-2)!!$, base cases $0!! = 1!! = 1$ | §6.2 |
| $\Gamma(z)$ | Euler Gamma function: $\Gamma(n+1) = n!$ for $n \in \mathbb{N}$ | §6.6 |
| $\binom{\alpha}{k}$ | Generalized binomial coefficient for $\alpha \in \mathbb{R}$, $k \in \mathbb{N}_0$ | §6.3 |
| $(\alpha)_k$ | Pochhammer rising factorial: $\alpha(\alpha+1)\cdots(\alpha+k-1)$ | §6.4 |
| $[\alpha]_k$ | Falling factorial: $\alpha(\alpha-1)\cdots(\alpha-k+1)$ | §6.4 |
| $W_n$ | Wallis integral: $W_n = \int_0^{\pi/2} \cos^n\varphi\,d\varphi$ | §6.5 |
| $u$ | Machine epsilon (unit roundoff) | §6.2 |

---

## Section Structure

### §6.1 Introduction

This section provides a road map of the chapter, establishes notation, and lists forward references.

- Road map, notation table, forward references to Ch 5 (binomial coefficients), Ch 14 (Wallis integrals in geodetic coefficients), Ch 15 (Pochhammer in Hansen coefficients)

### §6.2 Factorial and Double Factorial

This section defines the factorial and double factorial, relates them, and establishes the Stirling approximation for asymptotic use.

- Definition 6.2.1 (Factorial): n! = n·(n−1)!, 0! = 1
- Definition 6.2.2 (Double Factorial): n!! = n·(n−2)!!, base cases 0!! = 1, 1!! = 1, (−1)!! = 1
- Lemma 6.2.1 (Factorial via Double Factorial): (2m)! = 2^m · m! · (2m−1)!! — *Proof approach: induction on m*
- Theorem 6.2.1 (Stirling Approximation): ln(n!) = n ln n − n + ½ ln(2πn) + O(1/n) — *Proof approach: Laplace method on Gamma integral*
- Example 6.2.1: Tabulate n! and n!! for n = 0..10; verify Lemma 6.2.1 at m = 5: (10)! = 2^5 · 5! · 9!!. Source: direct integer computation.
- [M.6.1] Integer overflow boundary for 64-bit n!
- [P.6.1] Stirling accuracy for small n

### §6.3 Generalized Binomial Coefficients

This section defines generalized binomial coefficients for real upper index, establishes recurrences used in implementation, and characterizes sign patterns.

- Definition 6.3.1: C(α, k) = α(α−1)···(α−k+1)/k! for α ∈ ℝ, k ∈ ℕ₀
- Lemma 6.3.1 (Pascal Recurrence): C(α, k) = C(α−1, k−1) + C(α−1, k) — *Proof approach: direct algebraic expansion of falling factorial products*
- Lemma 6.3.2 (Recurrence in k): C(α, k) = C(α, k−1)·(α−k+1)/k [implementation recurrence] — *Proof approach: ratio of consecutive terms from Definition 6.3.1*
- Lemma 6.3.3 (Upper Negation): C(−α, k) = (−1)^k · C(α+k−1, k) — *Proof approach: rewrite falling factorial (−α)(−α−1)···(−α−k+1) and factor out (−1)^k*
- Theorem 6.3.1 (Sign Pattern): For α < 0, C(α, k) has sign (−1)^k — *Proof approach: immediate from Lemma 6.3.3 with α > 0*
- Corollary 6.3.1 (Finite Support): C(n, k) = 0 for k > n when n ∈ ℕ — *Proof approach: falling factorial contains factor (n−n) = 0*
- Example 6.3.1: Compute C(1/2, k) for k = 0..4; verify against closed form (−1)^{k+1}/(2k(2k−1)·C(2k,k)/4^k). Source: direct rational arithmetic.
- Example 6.3.2: Compute C(−1/2, k) for k = 0..4 via upper negation (Lemma 6.3.3); cross-check against Definition 6.3.1. Source: direct rational arithmetic.
- [P.6.2] Cancellation in falling factorial for α near non-negative integer

### §6.4 Pochhammer Symbols

This section defines rising and falling factorials, relates them by duality, and connects both to the Gamma function.

- Definition 6.4.1 (Rising Factorial): (α)_k = α(α+1)···(α+k−1)
- Definition 6.4.2 (Falling Factorial): [α]_k = α(α−1)···(α−k+1)
- Lemma 6.4.1 (Rising–Falling Duality): [α]_k = (−1)^k (−α)_k — *Proof approach: direct sign extraction from each factor*
- Theorem 6.4.1 (Gamma Representation): (α)_k = Γ(α+k)/Γ(α) — *Proof approach: induction on k using Γ(z+1) = zΓ(z)*
- Example 6.4.1: Compute (3/2)_k for k = 0..4; verify via Gamma ratio Γ(3/2+k)/Γ(3/2). Source: exact rational/half-integer Gamma values.

### §6.5 Wallis Integrals

This section derives the Wallis integral recurrence, gives closed forms for even and odd orders via double factorial, and proves the Wallis product formula.

- Definition 6.5.1: W_n = ∫₀^{π/2} cos^n φ dφ
- Theorem 6.5.1 (Wallis Recurrence): W_n = ((n−1)/n)·W_{n−2}, W_0 = π/2, W_1 = 1 — *Proof approach: integration by parts on cos^n φ = cos^{n−2} φ · cos² φ*
- Theorem 6.5.2 (Even Closed Form): W_{2m} = (π/2)·(2m−1)!!/(2m)!! — *Proof approach: unroll recurrence from W_0 = π/2*
- Theorem 6.5.3 (Odd Closed Form): W_{2m+1} = (2m)!!/(2m+1)!! — *Proof approach: unroll recurrence from W_1 = 1*
- Theorem 6.5.4 (Wallis Product): lim [(2m)!!/(2m−1)!!]²/(2m+1) = π/2 — *Proof approach: squeeze W_{2m+1} ≤ W_{2m} ≤ W_{2m−1} and take ratio*
- Example 6.5.1: Tabulate W_n for n = 0..10; verify each entry by recurrence, closed form, and direct numerical quadrature. Source: exact rational multiples of π.
- [P.6.3] Recurrence rounding vs direct double-factorial computation

### §6.6 Gamma Function Connection

This section introduces the Gamma function and uses it to unify factorials, binomial coefficients, Pochhammer symbols, and Wallis integrals into a single analytic framework.

- Definition 6.6.1 (Gamma Function): Γ(z) = ∫₀^∞ t^{z−1} e^{−t} dt
- Theorem 6.6.1: Γ(n+1) = n! — *Proof approach: integration by parts and induction*
- Theorem 6.6.2: Γ(1/2) = √π — *Proof approach: substitution t = u² reducing to Gaussian integral*
- Corollary 6.6.1: Γ(n+1/2) = (2n−1)!!·√π/2^n — *Proof approach: repeated application of Γ(z+1) = zΓ(z) from Γ(1/2)*
- Theorem 6.6.3 (Wallis via Gamma): W_n = (√π/2)·Γ((n+1)/2)/Γ((n+2)/2) — *Proof approach: Beta function identity B(a,b) = Γ(a)Γ(b)/Γ(a+b) with symmetry of cos^n integral*
- Theorem 6.6.4 (Reflection): Γ(z)Γ(1−z) = π/sin(πz) — *Proof approach: stated without proof (standard result); verified numerically*
- Theorem 6.6.5 (Duplication): Γ(z)Γ(z+1/2) = (√π/2^{2z−1})·Γ(2z) — *Proof approach: stated without proof (standard result); verified numerically*
- Corollary 6.6.2: C(α, k) = Γ(α+1)/(k!·Γ(α−k+1)) — *Proof approach: substitute Gamma representation of falling factorial from Definition 6.3.1*
- Example 6.6.1: Verify W_4 = 3π/16 and W_5 = 8/15 via closed form (Thm 6.5.2/6.5.3), Gamma formula (Thm 6.6.3), and recurrence (Thm 6.5.1). Source: exact rational multiples of π.
- [P.6.4] Gamma evaluation near poles
- [P.6.5] C(α, k) via Gamma ratio for large negative α

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| (none) | — | Ch 6 is foundational; no backward dependencies |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 5, §5.3 | §6.3 | Generalized binomial coefficients $C(\alpha, k)$ |
| Ch 5, §5.6.2 | §6.5 | Wallis integrals $W_n$ for geodetic series |
| Ch 14, §14.8 | §6.5 | Wallis integrals $W_{2k}$ for geodetic coefficient pattern |
| Ch 15, §15.2 | §6.4 | Pochhammer symbols for Hansen coefficients |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.6.1] | M | §6.2 | 64-bit overflow boundary for n! |
| [P.6.1] | P | §6.2 | Stirling accuracy for small n |
| [P.6.2] | P | §6.3 | Falling factorial cancellation for α near non-negative integer |
| [P.6.3] | P | §6.5 | Recurrence rounding vs direct double-factorial |
| [P.6.4] | P | §6.6 | Gamma evaluation near poles |
| [P.6.5] | P | §6.6 | Gamma ratio for large negative α |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 7 |
| Theorems | 12 |
| Lemmas | 5 |
| Corollaries | 3 |
| Propositions | 0 |
| Examples | 6 |
| Error Notes | 6 |
| Equations | ~20 |
| Sections | 6 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §6.1 | Draft | Introduction and forward references |
| §6.2 | Draft | Factorial and double factorial |
| §6.3 | Draft | Generalized binomial coefficients |
| §6.4 | Draft | Pochhammer symbols |
| §6.5 | Draft | Wallis integrals |
| §6.6 | Draft | Gamma function connection |
