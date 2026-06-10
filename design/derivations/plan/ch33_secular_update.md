# Draft Plan: Chapter 33 — Secular Update

**Part IX: The SGP4 Propagator** | Implementation file: `secular_update.h`

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $t_o$ | TLE epoch | Ch 31 |
| $\Delta t$ | Propagation time interval: $t - t_o$ | §33.1 |
| $M_{\mathrm{DF}}, \omega_{\mathrm{DF}}, \Omega_{\mathrm{DF}}$ | Drift-corrected mean elements at $t$ | §33.2 |
| $\dot{M}_0, \dot{\omega}_0, \dot{\Omega}_0$ | Brouwer secular rates at epoch | Ch 16–17 |
| $\delta\omega$ | Drag correction to argument of perigee | §33.3 |
| $\delta M$ | Drag correction to mean anomaly | §33.3 |
| $C_1, C_2, C_3, C_4, C_5$ | Drag-related precomputed coefficients | Ch 22 |
| $D_2, D_3, D_4$ | Time-polynomial coefficients for $a$ decay | Ch 22 |
| $e_t$ | Eccentricity at time $t$ | §33.4 |
| $a_t$ | Semi-major axis at time $t$ | §33.4 |
| $M_p$ | Mean anomaly at time $t$ before final correction | §33.5 |
| $L_{\mathrm{total}}$ | Mean longitude accumulation polynomial | §33.5 |

---

## Objectives

1. Derive the secular mean-element update formulas ($M_{\mathrm{DF}}, \omega_{\mathrm{DF}}, \Omega_{\mathrm{DF}}$) as linear-in-time applications of Brouwer rates.
2. Derive the drag corrections $\delta\omega$ and $\delta M$ from the drag coefficients $C_3, C_5$ (Ch 22).
3. Derive the $\Omega$ drag-gravity coupling (the 21/2 coefficient) from the cross-term between $J_2$ and drag.
4. Derive the eccentricity and semi-major axis decay formulas.
5. Derive the mean longitude polynomial through $(\Delta t)^5$.
6. Characterize the simplified drag mode (perigee < 220 km) and what is truncated.

## Section Structure

### §33.1 Introduction

This section describes the role of the secular update in the SGP4 pipeline and provides a road map of the chapter.

Stub: The secular update advances the mean orbital elements from epoch $t_o$ to time $t$ using rates precomputed in Ch 32 and Ch 22. No numerical integration is required; all corrections are polynomial in $\Delta t$. The pipeline: this chapter computes the updated mean elements; Ch 34 (or Ch 35 for deep space) then applies long-period and short-period corrections.

### §33.2 Mean Element Advance: M_DF, ω_DF, Ω_DF

This section derives the linear-in-time advance formulas for the three mean angular elements using the Brouwer secular rates as precomputed constants.

**Example 33.2.1** (ISS 24-hour propagation): For ISS ($n_o'' \approx 1.131\times10^{-3}$ rad/s, $i \approx 51.6°$), compute $\dot{M}_0$, $\dot{\omega}_0$, $\dot{\Omega}_0$ from Ch 16–17 formulas, then evaluate $M_{\mathrm{DF}}$, $\omega_{\mathrm{DF}}$, $\Omega_{\mathrm{DF}}$ at $\Delta t = 1440$ min. Expected: $\dot{\Omega}_0 \approx -5.06°$/day, $\Delta\Omega \approx -5.06°$. Source: Vallado et al. (2006) SGP4 verification test case NORAD 25544.

Stub: $M_{\mathrm{DF}} = M_o + \dot{M}_0 \Delta t$, $\omega_{\mathrm{DF}} = \omega_o + \dot{\omega}_0 \Delta t$, $\Omega_{\mathrm{DF}} = \Omega_o + \dot{\Omega}_0 \Delta t$. Here $\dot{M}_0, \dot{\omega}_0, \dot{\Omega}_0$ are the total first- and second-order Brouwer secular rates (Ch 16–17). Definition 33.2.1 (secular rates as precomputed constants): these are computed once from the epoch elements and reused at every propagation step. Cross-reference Ch 37 (precomputed constants). [P.33.1] $\dot{M}_0 \Delta t$ accumulates precision error for large $\Delta t$; compensated summation or multi-revolution representation recommended for $|\Delta t| > 10^6$ s.

### §33.3 Drag Corrections δω and δM

This section derives the drag corrections to the argument of perigee and mean anomaly from the $J_3$-eccentricity coupling coefficient $C_3$.

Stub: $\delta\omega = B^* C_3 \cos\omega_o \cdot \Delta t$. Derive from the $J_3$-eccentricity coupling (Ch 22, §22.5, coefficient $C_3$). $\delta M = M_{\mathrm{DF}} + B^* C_5[(\sin M_p - \sin M_o)] + \text{lower order}$ — derive the functional form. The $\cos\omega_o$ factor in $\delta\omega$: its physical origin is the $J_3$ odd-zonal coupling that selects $\sin\omega$ in the eccentricity perturbation (Ch 19). [A.33.1] Simplified drag mode: when perigee < 220 km, $C_3 = 0$ and $\delta M$ is reduced to the $C_1$-only term; derive the justification from the atmospheric model parameters.

### §33.4 Eccentricity and Semi-Major Axis Decay

This section derives the polynomial expressions for eccentricity and semi-major axis decay and justifies the squared form of the $a_t$ formula.

Stub: $e_t = e_o - B^* C_4 \Delta t - B^* C_5(\sin M_p - \sin M_o)$. Derive from the drag variational equation (Ch 22, §22.4). $a_t = a_o''[1 - C_1\Delta t - D_2(\Delta t)^2 - D_3(\Delta t)^3 - D_4(\Delta t)^4]^2$. Derive the squared form: it arises because $a$ decays as $(1 - \epsilon)^2 \approx 1 - 2\epsilon$ to first order, but the exact form is kept to reduce error for large $\Delta t$. Derive $D_2, D_3, D_4$ from the Taylor expansion of the drag-averaged semi-major axis decay (Ch 22, §22.6). [A.33.2] The $(\Delta t)^5$ term is dropped; its magnitude is bounded by $D_4(\Delta t)^4 \cdot C_1\Delta t$.

### §33.5 The Ω Drag-Gravity Coupling

This section derives the 21/2 coefficient in the $\Omega$ drag-gravity coupling by differentiating the Brouwer node regression rate with respect to semi-major axis and multiplying by the mean drag rate.

Stub: The 21/2 coefficient appears in $\Omega_{\mathrm{DF}}$ as a drag-induced node precession correction. Derive it from the cross-product of $J_2$ secular node regression and the drag-induced $\dot{a}$: $\frac{\partial\dot{\Omega}}{\partial a}\cdot\dot{a}_{\mathrm{drag}} \cdot \Delta t$ produces a term $-\tfrac{21}{2}\frac{k_2^2 n \Delta t \cos i}{p^2}\cdot\{$drag factor$\}$. Full derivation: express $\dot\Omega$ as a function of $a$ from Brouwer (Ch 16), differentiate with respect to $a$, multiply by the mean drag rate. Derive the exact prefactor 21/2 from this chain. [A.33.3] This is a first-order coupling; the second-order term is dropped.

### §33.6 Mean Longitude Polynomial

This section derives the fifth-order polynomial for the mean anomaly at time $t$ by integrating the drag-decaying mean motion and collecting coefficients by power of $\Delta t$.

Stub: The complete expression for the mean anomaly argument at time $t$ through fifth order in $\Delta t$: $L = M_{\mathrm{DF}} + \delta M + $ polynomial correction. Derive each coefficient from the corresponding $D_k$ and $C_1$ terms by integrating the semi-major axis decay formula. Show that the polynomial arises from the Taylor expansion of $n(a(t))$, where $a(t)$ is the decayed semi-major axis. Table 33.6.1: the polynomial coefficients with their derivation cross-references. [A.33.4] Fifth-order truncation; sixth-order term bounded by $D_4^2 (\Delta t)^8$ — negligible for propagation intervals less than several weeks.

### §33.7 Simplified Drag Mode

This section defines the simplified drag mode flag, states which pipeline stages are affected, and derives a bound on the accuracy cost.

Stub: When perigee altitude $< 220$ km: drop $D_2, D_3, D_4$ (all higher-order drag terms) and retain only $C_1$. Rationale: at very low altitudes, the orbit decays so rapidly that higher-order terms in $\Delta t$ provide negligible improvement over the span before re-entry. Formal bound on the dropped terms for perigee between 156 and 220 km. Definition 33.7.1 (simplified mode flag): set at initialization, governs §33.4 and §33.6.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 16–17, Brouwer secular rates | §33.2 | $\dot{M}_0, \dot{\omega}_0, \dot{\Omega}_0$ for mean element advance |
| Ch 19, $J_3$ long-period coupling | §33.3 | Origin of drag correction $\delta\omega$ |
| Ch 22, drag coefficients | §33.3, §33.4 | $C_1, \ldots, C_5, D_2, D_3, D_4$ derivations |
| Ch 32, element recovery | §33.2, §33.4 | Epoch elements $a_o'', n_o'', e_o, B^*$ |
| Ch 37, precomputed constants | §33.2 | Secular rates computed once at initialization |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 34, near-space pipeline | §33.2, §33.4 | Updated mean elements at propagation time $t$ |
| Ch 35, deep-space pipeline | §33.2, §33.4 | Updated mean elements at propagation time $t$ |
| Ch 38, output error budget | §33.6 | Polynomial truncation errors for total budget |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [P.33.1] | P | §33.2 | Mean anomaly accumulation for large $|\Delta t|$; precision error grows |
| [A.33.1] | A | §33.3 | Simplified drag mode drops $C_3, C_5$ corrections |
| [A.33.2] | A | §33.4 | Fifth-order polynomial truncation in $a(t)$ |
| [A.33.3] | A | §33.5 | First-order drag-gravity coupling only; 21/2 coefficient |
| [A.33.4] | A | §33.6 | Mean longitude polynomial truncated at fifth order |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 4 |
| Theorems | 1 |
| Lemmas | 1 |
| Corollaries | 0 |
| Propositions | 1 |
| Examples | 1 |
| Error Notes | 5 |
| Equations | ~15 |
| Sections | 7 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §33.1 | Draft | Introduction |
| §33.2 | Draft | Mean element advance: $M_{\mathrm{DF}}, \omega_{\mathrm{DF}}, \Omega_{\mathrm{DF}}$ |
| §33.3 | Draft | Drag corrections $\delta\omega$ and $\delta M$ |
| §33.4 | Draft | Eccentricity and semi-major axis decay |
| §33.5 | Draft | The $\Omega$ drag-gravity coupling |
| §33.6 | Draft | Mean longitude polynomial |
| §33.7 | Draft | Simplified drag mode |
