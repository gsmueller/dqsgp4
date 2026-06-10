# Draft Plan: Chapter 36 — The Propagator Architecture

**Part IX: The SGP4 Propagator** | Implementation files: `sgp4_propagator.h`, `model_functions.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathcal{M}$ | A ModelFunctions instance: a named collection of injectable lambdas | §36.2 |
| $\Lambda_{\mathrm{Kepler}}$ | Injectable lambda for the Kepler equation solver | §36.2 |
| $\Lambda_{\mathrm{drag}}$ | Injectable lambda for the drag coefficient computation | §36.2 |
| $\Lambda_{\mathrm{series}}$ | Injectable lambda for series evaluation to tolerance $\tau$ | §36.2 |
| $\tau$ | Tolerance parameter (Ch 3, Def. 3.4.1) | Ch 3 |
| $\tau_{\mathrm{standard}}$ | Standard tolerance: reproduces Hoots and Roehrich (1980) truncation | Ch 3 |

---

## Objectives

1. Motivate the injectable-lambda architecture as the implementation of the matched-pair and generalization principles.
2. Define the ModelFunctions interface: which computations are injectable, with what signatures.
3. Specify the six model presets and their physical justification.
4. Explain why computation is separated from theory: the propagator is a fixed control flow; the injected functions provide the mathematical content at each stage.
5. Show that setting all $\tau = \tau_{\mathrm{standard}}$ and all lambdas to their SGP4-standard implementations reproduces the Hoots and Roehrich (1980) model.

## Section Structure

### §36.1 Introduction: Separation of Theory and Implementation

This section motivates the injectable-lambda architecture as the computational realization of the matched-pair and generalization principles from Ch 3 and Ch 4.

Stub: The SGP4 propagator is a pipeline of mathematical operations. Each operation can be implemented at different accuracy levels — truncated (fast, compatible with SR3) or full (slower, accurate to tolerance $\tau$). The injectable-lambda architecture separates the control flow (which operations to apply, in which order) from the mathematical content of each operation (how accurately). This is the computational realization of the generalization principle (Ch 4, §4.9; Ch 3, §3.5).

### §36.2 The ModelFunctions Interface

This section defines the ModelFunctions interface, specifying the type signature and closure semantics of each injectable lambda.

**Definition 36.2.1** (ModelFunctions): a struct holding typed lambda functions for each injectable computation stage.

**Example 36.2.1** (Standard SGP4 lambda set): For the Standard SGP4 preset, $\Lambda_{\mathrm{Kepler}}$ is the Newton iteration of Ch 9 with 3 iterations (tolerance $\tau_{\mathrm{standard}}$), $\Lambda_{\mathrm{drag}}$ evaluates the Lane expansion at $O(e^4)$, and $\Lambda_{\mathrm{series}}$ truncates the generalized binomial series at 4 terms. Show that the composition reproduces the Hoots and Roehrich (1980) results to within floating-point equivalence. Source: Vallado et al. (2006) SGP4 verification TLE set.

Stub: Definition 36.2.1 (ModelFunctions): a struct holding typed lambda functions for each injectable computation stage. Required slots: (1) Kepler equation solver $\Lambda_{\mathrm{Kepler}}(M, e, \omega; \tau)$; (2) drag coefficient evaluator $\Lambda_{\mathrm{drag}}(a, e, i, B^*; \tau)$; (3) series evaluator $\Lambda_{\mathrm{series}}(x, \alpha, \tau)$ for generalized binomial series; (4) gravity field order $N_{\max}$; (5) atmosphere model selector (power-law vs. numerical). Each lambda captures precomputed constants in its closure (Ch 37 precomputed quantities). [A.36.1] The interface must remain stable across model presets; adding a new preset must not require changing the propagator control flow.

### §36.3 The Six Model Presets

This section specifies the six model presets by their $\tau$ settings and lambda implementations, and derives the expected accuracy gain for each step up.

Stub: Table 36.3.1: the six presets, their names, $\tau$ settings, and physical justification. Standard SGP4 preset: all $\tau = \tau_{\mathrm{standard}}$, all lambdas reproduce Hoots and Roehrich (1980) exactly. SGP4+: tighter $\tau$ in Kepler solver, full continued-fraction series, standard gravity coefficients. Full gravity: $N_{\max}$ up to 70×70. Full atmosphere: NRLMSISE-00 density instead of power law. High precision: all $\tau = 10^{-14}$, full gravity, full atmosphere. Research: all lambdas user-injectable for experimentation. Derive the expected accuracy gain for each step up in preset level, using the error-bound framework of Ch 1 and Ch 5.

### §36.4 Lambda Closures and Precomputed Constants

This section specifies the closure semantics of injectable lambdas and proves that capturing only const precomputed data makes the lambda thread-safe.

Stub: Each lambda is constructed at initialization and captures its precomputed constants via closure. The capture semantics: the lambda holds a const reference to the PrecomputedConstants struct (Ch 37); the struct is initialized once and never modified after construction. This is the implementation of the "compute once" principle.

**Proposition 36.4.1** (Thread safety): a lambda that captures only const data from PrecomputedConstants and per-step orbital-element arguments is thread-safe when called from separate threads with distinct element arguments. — *Proof approach: const data is read-only; each thread operates on a distinct copy of the orbital-element arguments; by definition, read-only shared state and thread-local mutable state produce no data races.*

### §36.5 Why Computation is Separated from Theory

This section explains how the injectable-lambda design enforces matched-pair consistency and enables independent upgrading of individual pipeline stages.

Stub: The propagator chapter (this chapter) specifies only the control flow and the interface. Chapters 31–35 and 37–38 each derive the mathematical content of one stage. The separation means: (1) any stage can be upgraded independently; (2) the theoretical derivation of each stage stands alone; (3) the matched-pair relationship is enforced by the preset system, not by mixing derived values from different presets. Cross-reference Ch 3 (matched pair), Ch 37 (initialization boundary).

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 3 | Matched pair principle, tolerance parameter | Preset system enforces matched-pair consistency |
| Ch 4 | Convergence acceleration | Tolerance-driven series evaluation |
| Ch 5 | Series evaluation with tolerance | Lambda signature for series evaluator |
| Ch 9–10 | Kepler equation solvers | Injectable Kepler solver implementations |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 37 | Precomputed constants | Lambda closures capture precomputed data |
| App C | Code-to-theorem mapping | `model_functions.h` audit reference |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.36.1] | A | §36.2 | Interface stability: any preset change must preserve propagator control flow |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 2 |
| Examples | 2 |
| Error Notes | 1 |
| Equations | ~10 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §36.1 Introduction: Separation of Theory and Implementation | Draft | |
| §36.2 The ModelFunctions Interface | Draft | |
| §36.3 The Six Model Presets | Draft | |
| §36.4 Lambda Closures and Precomputed Constants | Draft | |
| §36.5 Why Computation is Separated from Theory | Draft | |
