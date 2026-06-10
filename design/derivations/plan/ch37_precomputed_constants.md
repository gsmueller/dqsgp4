# Draft Plan: Chapter 37 — Precomputed Constants

**Part IX: The SGP4 Propagator** | Implementation file: `precomputed.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathcal{P}_{\mathrm{sat}}$ | Satellite-specific precomputed constants (depend on TLE) | §37.2 |
| $\mathcal{P}_{\mathrm{model}}$ | Model-specific precomputed constants (independent of TLE) | §37.3 |
| $t_{\mathrm{init}}$ | Initialization time: when precomputed constants are computed | §37.1 |

---

## Objectives

1. Enumerate all constants that must be computed exactly once at initialization, organized by category.
2. Establish the initialization-vs-propagation boundary: which quantities depend on TLE only ($\mathcal{P}_{\mathrm{sat}}$) vs. on model only ($\mathcal{P}_{\mathrm{model}}$) vs. on $t$ (per-step).
3. Show that this organization implements the "compute once" principle without sacrificing generality.
4. Document why each constant cannot be deferred to propagation time without either redundant computation or incorrect results.

## Section Structure

### §37.1 The Initialization Principle

This section formally defines precomputed constants, proves the correctness criterion for precomputation, and establishes the initialization-vs-propagation boundary.

**Definition 37.1.1** (precomputed constant): a quantity whose value depends only on the TLE epoch elements and the model parameters, not on the propagation time $t$. All precomputed constants are computed once at initialization time $t_{\mathrm{init}}$ and stored in a const struct.

**Theorem 37.1.1** (correctness of precomputation): a quantity $f(e_o, i_o, a_o'', B^*)$ is correctly precomputed if and only if it does not depend on $t$ through $a_t, e_t$, or any angle updated during secular propagation. — *Proof approach: inspect the secular update formulas (Ch 33) and the long/short-period correction formulas (Ch 18–19); these depend on $t$ only through $\Delta t$ and the time-evolved $a_t, e_t$ derived from it; any formula depending only on the epoch values $(a_o'', e_o, i_o, B^*)$ is therefore independent of $t$.*

**Example 37.1.1** (correctness check): Verify that the Brouwer secular rate $\dot{\Omega}_0 = -\tfrac{3}{2}k_2 n_o'' \cos i_o / (a_o''{}^2 (1-e_o^2)^2)$ satisfies the precomputation criterion (depends only on epoch elements). Verify that $a_t = a_o''[1 - C_1\Delta t - \ldots]^2$ does not (it depends on $\Delta t$). Source: Ch 33, §33.2 and §33.4 formulas.

### §37.2 Satellite-Specific Constants

This section provides a complete annotated table of $\mathcal{P}_{\mathrm{sat}}$ with each entry's derivation chapter, formula, and simplified-drag-mode dependency.

Stub: Complete annotated table of $\mathcal{P}_{\mathrm{sat}}$: for each entry, name, formula, chapter where derived, and whether it depends on simplified-drag mode. Groups: element recovery ($n_o'', a_o'', s^*, (q_0-s^*)^4$); drag coefficients ($C_1, C_2, C_3, C_4, C_5, D_2, D_3, D_4$); Brouwer secular rates ($\dot{M}_0, \dot{\omega}_0, \dot{\Omega}_0$); Brouwer short-period coefficient groups (from Ch 18); long-period coefficient groups (from Ch 19); deep-space initialization outputs (from Ch 35, §35.2, when applicable).

### §37.3 Model-Specific Constants

This section provides a complete annotated table of $\mathcal{P}_{\mathrm{model}}$, distinguishing universal physical constants from model-dependent parameters.

Stub: Complete annotated table of $\mathcal{P}_{\mathrm{model}}$: physical constants from Appendix A ($k_e, k_2, k_4, A_{3,0}, a_E$); Kaula inclination function table $F_{lmp}(i_o)$ for $l \leq 5$ evaluated at the satellite's inclination (satellite-dependent but stored with model for efficiency); WGS72 harmonic coefficients $C_{nm}, S_{nm}$ for resonance computation.

### §37.4 The Initialization-vs-Propagation Boundary

This section gives the formal statement of the boundary and identifies the special case of lambda closures that hold const references to $\mathcal{P}_{\mathrm{sat}}$.

Stub: Formal statement of the boundary: $\mathcal{P}_{\mathrm{sat}} \cup \mathcal{P}_{\mathrm{model}}$ is everything computed before the first call to propagate$(t)$. Everything depending on $t$ (including $\theta_{\mathrm{GMST}}(t)$, which depends on the propagation epoch through $T_{\mathrm{UT1}}$) is per-step. Special case: the closure captures (Ch 36, §36.4) mean the lambda functions carry references into $\mathcal{P}_{\mathrm{sat}}$; they must not modify it.

### §37.5 Implementation Notes

This section documents the initialization order and the simplified-drag-mode branching that determines which members of $\mathcal{P}_{\mathrm{sat}}$ are populated.

Stub: The PrecomputedConstants struct is initialized in a single pass through the initialization function. The initialization function calls, in order: TLE parsing (Ch 31), element recovery (Ch 32), drag coefficient computation (Ch 22), secular rate computation (Ch 16–17), short-period coefficient computation (Ch 18), long-period coefficient computation (Ch 19), and — if deep-space — DPINIT (Ch 35). Cross-reference each initialization step to the chapter that derives the quantities. Note: simplified-drag mode changes which members of $\mathcal{P}_{\mathrm{sat}}$ are populated.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Part IX chapters (Ch 31–35) | Various | Each contributes entries to $\mathcal{P}_{\mathrm{sat}}$ |
| Ch 16–19 | Secular rates, short/long-period coefficients | Initialization-order dependencies |
| Ch 22 | Drag coefficients | Drag coefficient group in $\mathcal{P}_{\mathrm{sat}}$ |
| App A | Physical constants | Model-specific constants $\mathcal{P}_{\mathrm{model}}$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 36 | Lambda closures | Closures capture const references to $\mathcal{P}_{\mathrm{sat}}$ |
| App C | Code mapping | `precomputed.h` audit reference |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.37.1] | A | §37.2 | SGP4 must use WGS-72 gravitational constants ($\mu$, $J_2$, $J_3$, $J_4$, $a_E$, $k_e$), not WGS-84; the TLE-fitting process used WGS-72 parameters (matched-pair principle, Ch 3); substituting WGS-84 constants produces erroneous positions; documented by Vallado et al. (2006) AIAA 2006-6753 |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 1 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 0 |
| Examples | 1 |
| Error Notes | 1 |
| Equations | ~15 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §37.1 The Initialization Principle | Draft | |
| §37.2 Satellite-Specific Constants | Draft | |
| §37.3 Model-Specific Constants | Draft | |
| §37.4 The Initialization-vs-Propagation Boundary | Draft | |
| §37.5 Implementation Notes | Draft | |
