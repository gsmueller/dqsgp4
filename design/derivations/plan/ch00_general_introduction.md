# Draft Plan: Chapter 0 — General Introduction

**Part:** Prefatory matter (before Part I)

## Objectives

1. Orient the reader to the scope and philosophy of the textbook without pre-empting derivations.
2. Establish the narrative arc: observable phenomena → physical theories → mathematical approximations → computational error budget.
3. Motivate the two-dual-quaternion state representation as the natural unifying structure.
4. Frame the "generalization" theme: every SGP4 shortcut is paired with the path to recover the dropped accuracy.
5. Serve as a reading guide: how the nine parts relate, what each part depends on, and how a reader with specific goals can navigate selectively.

## Section Structure

### §0.1 Introduction

This section frames the central problem, provides a roadmap of the chapter, and establishes forward references.

- Forward-reference table:

| Section | Feeds | Role |
|---------|-------|------|
| §0.2 | Ch 13–15, Ch 21–22, Ch 23–27, Ch 28 | Qualitative overview of physical phenomena |
| §0.3 | Ch 11–12, Ch 16–19 | Mathematical theories overview |
| §0.4 | Ch 1 | Informal motivation for three-error framework |
| §0.5 | Ch 2 | Informal motivation for dual-quaternion state |
| §0.6 | Ch 4, Ch 36 | Generalization philosophy |

- Maturity: Draft

### §0.2 The Central Problem and Physical Phenomena

This section surveys the physical phenomena that perturb satellite orbits, ordered by magnitude, to motivate the mathematical apparatus developed in subsequent parts.

Stub: A satellite position is needed at time $t$. The observable inputs are a TLE. The required output is a state vector in some Earth-fixed frame, with an honest account of how many digits are reliable. Describe the gap between the two and why bridging it requires the full apparatus of the textbook.

Enumerate the physical effects in order of magnitude: Keplerian motion, $J_2$ oblateness, atmospheric drag, higher zonals, tesseral resonance, third-body (Moon, Sun), solar radiation pressure (out of scope). One paragraph per effect. Figures: schematic of orbit precession, drag spiral, resonance trapping. No derivations — qualitative descriptions only. Forward-reference to the relevant Part.

### §0.3 A Map of Mathematical Theories

This section traces the chain of mathematical theories — from Hamiltonian mechanics through perturbation theory to SGP4 — identifying what each link gains and sacrifices.

Stub: Each phenomenon is handled by a mathematical theory. List the chain: Hamiltonian mechanics → perturbation theory (von Zeipel, Lie) → Brouwer's solution → Kozai's simplification → SGP4. Note what each step gains and what it sacrifices. Identify where Kaula's expansion enters (tesseral terms). Note why the power-law density model is chosen over a physical atmosphere.

### §0.4 Approximations and the Error Budget

This section informally introduces the three-error framework (measurement, precision, model) and motivates the claim that reliable digit counts can be tracked through every computation.

Stub: Every approximation in the chain has a cost. Define the three error categories informally (measurement, precision, model) and explain why they must be tracked separately. State the central claim: at the output of a carefully implemented propagator, it is possible to know, for each coordinate, how many decimal digits are reliable. The formal machinery is in Ch 1; this section is the motivation.

### §0.5 Why Two Dual Quaternions

This section motivates the two-dual-quaternion state representation (M-hat, Omega-hat) as the natural algebraic structure that makes propagation compositional.

Stub: Explain the configuration dual quaternion M̂ = M + εD (orientation paired with position via the dual part D = ½R_pos·M) and the velocity dual quaternion Ω̂_b = Ω_b + εV_b. State informally why a single algebraic object encoding both geometric and dynamic quantities makes the entire propagation pipeline compositional: every frame change, every perturbation step, every output transform is a multiplication or adjoint action on the pair (M̂, Ω̂). The linear equivalent — the 7×7 homogeneous matrix acting on (r, v, 1)^T — is fully equivalent and equally developed. Cross-reference §2.3 (configuration dual quaternion), §2.4 (velocity dual quaternion), §2.8 (7×7 conversion).

### §0.6 The Generalization Philosophy

This section establishes the three-register presentation strategy — exact expression, SGP4 truncation with error tag, and generalized evaluation to caller-specified tolerance — that recurs throughout the textbook.

Stub: SGP4 was designed for speed on 1970s hardware. Every truncation made sense in that context. This textbook presents every result in three registers: (1) the exact, untruncated expression; (2) the SGP4 truncation with its error tag; (3) the technique — continued fraction, Padé approximant, higher-order series — to evaluate the full expression to caller-specified tolerance. The code architecture (Ch 36, injectable lambdas) makes it possible to substitute the full computation for the truncated one at runtime. This section is a compact reading guide to the Generalization sections that appear throughout the textbook.

### §0.7 A Reading Guide

This section provides a prerequisite map of the nine parts and recommended reading paths for three reader profiles.

Stub: Table organizing the nine Parts by prerequisite structure. Paths for three reader profiles: (a) orbit propagation user who only wants to understand SGP4 — read Parts I, IV, VIII, IX; (b) geodesist interested in the gravity field — read Parts I, III, IV; (c) mathematician interested in perturbation theory — read Parts I, II, IV, V, VI.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| (none) | — | Chapter 0 is the entry point; no backward dependencies |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 1 | §0.4 | Motivates the three-error framework |
| Ch 2 | §0.5 | Motivates the dual-quaternion state representation |
| Ch 3 | §0.5 | Motivates the matched-pair principle |
| Ch 4 | §0.6 | Motivates the approximation/generalization toolkit |
| Ch 36 | §0.6 | Motivates the injectable-lambda architecture |
| App B | §0.7 | Reading guide references notation index |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.0.1] | Accuracy | §0.2 | Qualitative magnitudes are representative for 400 km circular; vary by orbit type |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 0 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 0 |
| Examples | 0 |
| Error Notes | 1 |
| Equations | 0 |
| Sections | 7 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §0.1 | Draft | Forward-reference table stub present |
| §0.2 | Draft | Qualitative narrative, no formal items |
| §0.3 | Draft | Qualitative narrative, no formal items |
| §0.4 | Draft | Informal preview of Ch 1 framework |
| §0.5 | Draft | Informal preview of Ch 2 framework |
| §0.6 | Draft | Generalization philosophy narrative |
| §0.7 | Draft | Reading guide with table |
