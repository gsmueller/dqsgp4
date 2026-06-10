# Draft Plan: Chapter 7 — Angle Arithmetic

## Objectives

1. Establish exact degree-radian conversion and characterize floating-point error from finite π approximation.
2. Define modular reduction to [0, 2π) and [−π, π); analyze precision loss for large arguments.
3. Present Cody-Waite extended-precision argument reduction.
4. Prove Lipschitz-1 property of modular reduction; derive error propagation bounds.
5. Introduce multi-revolution angle representation for preserving precision over many orbital periods.

**Implements:** `angles.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\theta$ | Angle value (radians unless otherwise stated) | §7.1 |
| $\pi$ | Mathematical constant (irrational, $\approx 3.14159265358979$) | §7.2 |
| $\mathrm{fl}(x)$ | Floating-point representation of $x$ | §7.2 |
| $u$ | Machine epsilon (unit roundoff, $2^{-53}$ for IEEE 754 double) | §7.2 |
| $\mathrm{mod}(\theta, 2\pi)$ | Canonical reduction to $[0, 2\pi)$ | §7.3 |
| $\mathrm{wrap}(\theta)$ | Canonical reduction to $[-\pi, \pi)$ | §7.3 |
| $C_1, C_2, C_3$ | Cody-Waite extended-precision split of $2\pi$ | §7.3 |
| $\theta_{\mathrm{rev}}$ | Multi-revolution angle representation | §7.5 |

---

## Section Structure

### §7.1 Introduction

This section motivates angle arithmetic as a precision-critical subsystem and provides a road map of the chapter.

- Angle precision is a dominant error source at large elapsed times
- Application-independent results; forward references to Ch 8, 9, 16–19, 29

### §7.2 Degree-Radian Conversion

This section establishes that degree-to-radian conversion introduces an irreducible representational error and quantifies its magnitude.

- Definition 7.2.1 (Conversion Factors): deg_to_rad = π/180 (exact irrational)
- Lemma 7.2.1 (Conversion Error): fl(π/180) has relative error ε_π ≈ u(x) — *Proof approach: direct computation of fl(π)/fl(180) vs exact π/180*
- Lemma 7.2.2 (Round-Trip Error): deg→rad→deg loses at most 2u(x)·|θ| — *Proof approach: chain two multiplication error bounds*
- Theorem 7.2.1 (Non-Exactness): no binary float exactly represents π/180 — *Proof approach: π is transcendental, so π/180 is irrational, hence not in dyadic rationals*
- Example 7.2.1: Compute fl(π/180) in double precision, compare to arbitrary-precision π/180; report relative error in ulps. Source: IEEE 754 double-precision arithmetic.
- [P.7.1] Relative error in every radian value from degree conversion
- [P.7.2] Round-trip conversion error

### §7.3 Modular Reduction

This section defines modular angle reduction, quantifies the catastrophic cancellation in naive fmod, and presents the Cody-Waite extended-precision technique.

- Definition 7.3.1 (Canonical Range [0, 2π)): mod(θ, 2π) = unique φ ∈ [0, 2π)
- Definition 7.3.2 (Canonical Range [−π, π)): wrap(θ)
- Definition 7.3.3 (Cody-Waite Reduction): extended-precision splitting 2π ≈ C₁+C₂+C₃
- Lemma 7.3.1 (Naive fmod Precision Loss): |φ̂ − φ| ≲ N·u(x)·2π for N revolutions — *Proof approach: error model for fl(θ − N·fl(2π)) with N = ⌊θ/(2π)⌋*
- Theorem 7.3.1 (Cody-Waite Error Bound): |φ̂ − φ| ≤ C·u(x)^k·|θ| with k pieces — *Proof approach: telescoping cancellation in (θ − N·C₁) − N·C₂ − N·C₃ with exact high-order subtraction*
- Example 7.3.1: Naive fmod on θ = 2π·10^6 + 0.1; count bits lost. Source: IEEE 754 double-precision fmod vs arbitrary-precision reference.
- Example 7.3.2: Cody-Waite 3-piece reduction on same θ = 2π·10^6 + 0.1; compare residual error to naive. Source: Cody-Waite constants from Cody and Waite (1980).
- [P.7.3] Bit loss as function of revolution count N
- [P.7.4] Cody-Waite residual for 2-piece vs 3-piece

### §7.4 Error Propagation in Angle Arithmetic

This section proves that modular reduction and trigonometric evaluation are Lipschitz-1, bounding how angle errors propagate through the computation pipeline.

- Definition 7.4.1 (Lipschitz-1 Map)
- Theorem 7.4.1: mod(·, 2π) and wrap(·) are Lipschitz-1 — *Proof approach: case analysis on whether reduction crosses the branch cut*
- Corollary 7.4.1: Error bound after reduction (with branch-cut caveat) — *Proof approach: direct application of Lipschitz-1 property*
- Theorem 7.4.2: |sin(φ̂) − sin(φ)| ≤ |φ̂ − φ| (trig error amplification ≤ 1) — *Proof approach: mean value theorem, |cos| ≤ 1*
- Lemma 7.4.1 (Digit Loss): ⌊log₂(ωt/(2π))⌋ bits lost before reduction — *Proof approach: floating-point subtraction cancellation count*
- Example 7.4.1: Compute sin(θ) for θ = 2π·10^6 + 0.1 with naive fmod vs Cody-Waite reduction; report absolute error vs arbitrary-precision sin. Source: IEEE 754 double precision, arbitrary-precision reference via mpfr.
- [P.7.5] Digit loss formula
- [A.7.1] sin/cos accuracy after Cody-Waite vs naive for representative N

### §7.5 Multi-Revolution Representation

This section introduces the multi-revolution angle representation that preserves full floating-point precision regardless of revolution count.

- Definition 7.5.1: Multi-revolution angle (N, φ) ∈ ℤ × [0, 2π), θ = N·2π + φ
- Theorem 7.5.1 (Precision Preservation): (N, φ) with Cody-Waite gives O(u²·|θ|) error vs O(u·|θ|) naive — *Proof approach: error analysis of incremental update φ_{i+1} = CW(φ_i + Δφ) vs monolithic θ = ωt*
- Lemma 7.5.1 (Addition): (N₁,φ₁) + (N₂,φ₂) with carry k ∈ {0,1} — *Proof approach: direct algebraic manipulation with floor division*
- Corollary 7.5.1: sin/cos evaluation uses only φ; N cancels by periodicity — *Proof approach: sin(2πN + φ) = sin(φ) by periodicity*
- Theorem 7.5.2 (Precision Budget): bits remaining = max(0, 53 − ⌊log₂ N⌋) naive; full 53 with multi-rev — *Proof approach: count significand bits consumed by integer part N*
- Example 7.5.1: LEO satellite (T ≈ 90 min), 1 day propagation: compute N ≈ 16, φ; verify all 53 bits retained in multi-rev vs 49 bits naive. Source: representative LEO mean motion n ≈ 0.0011 rad/s.
- Example 7.5.2: 10-year propagation of same LEO orbit: N ≈ 58,400; quantify naive bit-loss (≈16 bits) vs full precision in multi-rev. Source: same LEO parameters.
- [P.7.6] Precision budget for LEO/MEO/GEO/deep-space
- [P.7.7] Cody-Waite improvement factor in multi-rev

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, Cor 1.4.9 | §7.4 | Modular reduction bound |
| Ch 1, Cor 1.4.1–1.4.2 | §7.4 | Sine/cosine error bounds |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 8 | §7.2–7.5 | Orbital angle conversions and precision |
| Ch 9 | §7.3–7.4 | Kepler equation angle reduction |
| Ch 16–19 | §7.3–7.5 | Secular/periodic angle accumulation precision |
| Ch 29 | §7.3 | Sidereal time argument reduction |
| Ch 33 | §7.5 | Multi-revolution mean longitude propagation |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.7.1] | P | §7.2 | Radian conversion error from finite π approximation |
| [P.7.2] | P | §7.2 | Round-trip deg→rad→deg loses at most 2u(x)·|θ| |
| [P.7.3] | P | §7.3 | Naive fmod bit loss: ~log₂(N) bits for N revolutions |
| [P.7.4] | P | §7.3 | Cody-Waite residual for 2-piece vs 3-piece |
| [P.7.5] | P | §7.4 | Digit loss formula: ⌊log₂(ωt/(2π))⌋ bits |
| [P.7.6] | P | §7.5 | Precision budget by orbit type (LEO/MEO/GEO/deep) |
| [P.7.7] | P | §7.5 | Cody-Waite improvement factor in multi-rev |
| [A.7.1] | A | §7.4 | Trig accuracy after reduction vs naive for representative N |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 6 |
| Theorems | 6 |
| Lemmas | 5 |
| Corollaries | 2 |
| Propositions | 0 |
| Examples | 6 |
| Error Notes | 8 |
| Equations | ~15 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §7.1 | Draft | Introduction and forward references |
| §7.2 | Draft | Degree-radian conversion |
| §7.3 | Draft | Modular reduction and Cody-Waite |
| §7.4 | Draft | Error propagation in angle arithmetic |
| §7.5 | Draft | Multi-revolution representation |
