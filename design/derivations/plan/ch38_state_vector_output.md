# Draft Plan: Chapter 38 — The State Vector Output

**Part IX: The SGP4 Propagator** | Implementation file: `state_vector.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $\mathbf{r}_{b,\mathrm{TEME}}$ | Position in TEME frame (km) | Ch 30 |
| $\mathbf{v}_{b,\mathrm{TEME}}$ | Velocity in TEME frame (km/s) | Ch 30 |
| $\sigma_m^{\mathrm{out}}$ | Measurement error at output | §38.2 |
| $\delta_p^{\mathrm{out}}$ | Precision error at output | §38.2 |
| $\delta_a^{\mathrm{out}}$ | Accuracy/model error at output | §38.2 |
| $N_r$ | Reliable digits in position: $-\log_{10}(\delta_{\mathrm{total}}/|r|)$ | §38.3 |
| $N_v$ | Reliable digits in velocity: $-\log_{10}(\delta_{\mathrm{total}}/|v|)$ | §38.3 |

---

## Objectives

1. Define the state vector output: TEME position (km) and velocity (km/s) with their complete TrackedValues.
2. Perform the three-error decomposition at output: characterize $\sigma_m, \delta_p, \delta_a$ independently.
3. Define and compute reliable digit counts $N_r, N_v$.
4. Construct the complete end-to-end error budget from TLE input to state vector output.
5. Identify the dominant error source at each propagation time scale.

## Section Structure

### §38.1 Output Definition

This section defines the state vector output as a pair of TrackedValue 3-vectors in the TEME frame, specifying units and frame limitations.

**Definition 38.1.1** (state vector output): the pair $(\mathbf{r}_{b,\mathrm{TEME}}, \mathbf{v}_{b,\mathrm{TEME}})$, each a 3-vector of TrackedValues in the TEME frame as defined in Ch 30.

Stub: Definition 38.1.1 (state vector output): the pair $(\mathbf{r}_{b,\mathrm{TEME}}, \mathbf{v}_{b,\mathrm{TEME}})$, each a 3-vector of TrackedValues, in the TEME frame as defined in Ch 30. Units: km and km/s (not SI meters) to match the TLE convention. Note the frame limitation: TEME is not standard ECI; Ch 30 gives the chain to ECEF or J2000 ECI.

### §38.2 Three-Error Decomposition

This section decomposes the output error into its three independent components and derives each component's growth rate with propagation time.

**Theorem 38.2.1** (output error decomposition): At output, the three error components $(\sigma_m^{\mathrm{out}}, \delta_p^{\mathrm{out}}, \delta_a^{\mathrm{out}})$ are each bounded by explicit expressions in terms of per-stage contributions and propagation time $|\Delta t|$. — *Proof approach: apply the error propagation rules of Ch 1 (TrackedValue composition, Thm 1.3.1) stage by stage through the SGP4 pipeline (Ch 31–35), carrying each error component independently; the triangle inequality bounds the cumulative total.*

**Example 38.2.1** (LEO error at 1 day): ISS TLE propagated 24 hours: compute $\sigma_m^{\mathrm{out}} \approx 10$ m (from TLE digit precision), $\delta_p^{\mathrm{out}} \approx 0.1$ m (from arithmetic rounding), $\delta_a^{\mathrm{out}} \approx 1$–$10$ km (from Brouwer model truncation and drag uncertainty). Total position uncertainty $\approx \delta_a^{\mathrm{out}}$. Source: SR3 validation data; Vallado et al. (2006) comparison against precision ephemerides.

Stub: Theorem 38.2.1 (output error decomposition). At output, the three error components are: (1) $\sigma_m^{\mathrm{out}}$: the measurement error propagated from the TLE through the full pipeline; this component reflects the finite precision of TLE encoding and is bounded by the initial $\sigma_m$ scaled by the condition number of the pipeline. (2) $\delta_p^{\mathrm{out}}$: the arithmetic precision error accumulated through all computations; dominated by the stage with the most subtractive cancellation (the short-period correction near $e \to 0$ or $i \to 0$). (3) $\delta_a^{\mathrm{out}}$: the model error from all accuracy approximations made in the pipeline; at $t = 0$ this is zero (the SGP4 model is self-consistent at epoch); it grows with $|\Delta t|$ as the physical orbit deviates from the model. Proof: by applying the error propagation rules of Ch 1 through each stage of the pipeline, carrying each error component separately. ∎ [A.38.1] The dominant error source for $|\Delta t| < 1$ day is the model error $\delta_a$ from the first-order Brouwer theory. For $|\Delta t| > 7$ days, drag model error dominates for LEO; tesseral resonance error dominates for GEO.

### §38.3 Reliable Digits

This section defines the reliable digit count for position and velocity and evaluates it for representative LEO and GEO orbits at several propagation times.

**Definition 38.3.1** (reliable digits): $N_r = -\log_{10}(\delta_{\mathrm{total},r} / |\mathbf{r}|)$ where $\delta_{\mathrm{total},r} = \sigma_m^{\mathrm{out}} + \delta_p^{\mathrm{out}} + \delta_a^{\mathrm{out}}$.

**Example 38.3.1** (reliable digits for LEO/GEO): (a) 400 km LEO at $t = 0$: $\delta_a = 0$, $\sigma_m \approx 10$ m, $|\mathbf{r}| \approx 6778$ km → $N_r \approx 9.8$. (b) 400 km LEO at $t = 1$ day: $\delta_a \approx 1$ km → $N_r \approx 3.8$. (c) GEO at $t = 1$ day: $\delta_a \approx 0.1$ km → $N_r \approx 5.6$. Source: SGP4 accuracy estimates from Vallado et al. (2006).

Stub: Definition 38.3.1 (reliable digits): $N_r = -\log_{10}(\delta_{\mathrm{total},r} / |\mathbf{r}|)$ where $\delta_{\mathrm{total},r} = \sigma_m^{\mathrm{out}} + \delta_p^{\mathrm{out}} + \delta_a^{\mathrm{out}}$ for position. Analogously for velocity. Note that $N_r$ and $N_v$ are scalar summaries; they can be computed per-component or for the full vector norm. Example: for a 400 km LEO at $t = 0$, typical values are $N_r \approx 10$ (limited by TLE measurement precision) and $\delta_a^{\mathrm{out}} \approx 0$. At $t = 1$ day: $N_r \approx 8$ (model error ~10 m out of ~7000 km position magnitude). [P.38.1] The reliable-digit formula uses $\log_{10}$, which amplifies near-zero errors when $|\mathbf{r}|$ is large; report $N_r$ with a stated scale rather than as an absolute.

### §38.4 The Complete Error Budget

This section assembles a complete end-to-end error budget tracing each component from its source in TLE parsing through every pipeline stage to the final state vector.

Stub: Table 38.4.1: full end-to-end error budget tracing each error component from its source to the output. Rows: TLE parsing ($\sigma_m$ from digit count), element recovery (precision from cube root), secular update (precision from polynomial accumulation, model from Brouwer truncation), long-period correction (model from $J_3$ first-order), short-period correction (model and precision near singularities), Kepler step (precision from Newton convergence), coordinate transform (precision from trig, model from TEME approximation). Columns: chapter, error type, typical magnitude for a 400 km LEO at $t=0$ and at $t=1$ day. [A.38.2] The table entries are order-of-magnitude estimates; exact magnitudes are orbit-dependent and must be computed per TrackedValue at runtime.

### §38.5 Dominant Error Sources by Time Scale

This section qualitatively characterizes which error source dominates at each propagation time scale and identifies the regime boundaries for LEO, MEO, and GEO orbit types.

Stub: Qualitative discussion: At $t = 0$: all error from TLE measurement ($\sigma_m$) plus negligible precision; model error identically zero by matched-pair construction. For $0 < |\Delta t| < 1$ day: short-period model error grows; drag model error begins. For $1 < |\Delta t| < 7$ days: drag model error dominates for LEO; tesseral resonance error for GEO/GPS. For $|\Delta t| > 7$ days: non-modeled drag variations (solar flux), higher-order gravity, relativity, solar radiation pressure. [A.38.3] SR3 validation studies show ~1 km position error at 1 day for average LEO; this is dominated by $\delta_a$, not $\delta_p$.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Thm 1.3.1 (TrackedValue), three error categories | Error decomposition framework |
| Ch 30 | TEME frame definition | Output coordinate frame |
| Ch 33–35 | Pipeline stages (secular, near-space, deep-space) | Error contributions traced through pipeline |
| Ch 37 | Precomputed constants | Initialized error-free at $t = t_o$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| App B | Error symbol notation | $\sigma_m^{\mathrm{out}}, \delta_p^{\mathrm{out}}, \delta_a^{\mathrm{out}}$ definitions |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.38.1] | P | §38.3 | Reliable-digit formula amplifies near-zero errors when $|\mathbf{r}|$ large |
| [A.38.1] | A | §38.2 | Dominant error source transitions with $|\Delta t|$ |
| [A.38.2] | A | §38.4 | Error budget entries order-of-magnitude; exact values orbit-dependent |
| [A.38.3] | A | §38.5 | SR3 validation data provides $\delta_a$ benchmark |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 3 |
| Theorems | 1 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 0 |
| Examples | 3 |
| Error Notes | 4 |
| Equations | ~15 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §38.1 Output Definition | Draft | |
| §38.2 Three-Error Decomposition | Draft | |
| §38.3 Reliable Digits | Draft | |
| §38.4 The Complete Error Budget | Draft | |
| §38.5 Dominant Error Sources by Time Scale | Draft | |
