# Draft Plan: Chapter 22 — Orbit-Averaged Drag Coefficients

**Part V: Atmospheric Drag**

Target header: `drag_coefficients.h`

---

## Objectives

1. State the drag force law and the Gauss variational equation for the semi-major axis decay under drag.
2. Derive the orbit-averaged drag integral for a general power-law exponent $\tau$.
3. Carry out the combinatorial expansion (Lane 1965) that converts the orbit-average into a sum of binomial coefficients and Wallis-type integrals.
4. Evaluate the general result at $\tau = 4$ to obtain the AFGP4 expression for the mean-motion derivative $C_2$.
5. Derive the $J_2$-density coupling term and the decomposition $C_2 = C_2^A + C_2^B$.
6. Derive the $3/8\,J_2$ coefficient from first principles.
7. Identify terms dropped in the AFGP4 $\to$ SGP4 transition and bound their magnitude.
8. Derive the higher-order coefficients $C_3$, $C_4$, $C_5$ and time-polynomial coefficients $D_2$, $D_3$, $D_4$.
9. State the generalization to arbitrary $\tau$, restored $O(e^2)$ terms, and co-rotating atmosphere.

---

## Section Structure

### §22.1 Introduction

This section states the chapter scope: derivation of all orbit-averaged drag coefficients from the Lane combinatorial method.

*Stub.* This chapter derives the analytic orbit-averaged drag coefficients that govern secular orbital decay. The derivation follows the Lane (1965) combinatorial method and the AFGP4 elaboration (Lane, Cranford 1969; Hoots, Roehrich 1980). Every coefficient is derived algebraically; no values are accepted from reference code without proof.

The chapter is the longest in Part V: it contains the principal quantitative machinery of the atmospheric drag model.

---

### §22.2 The Drag Force and Gauss Variational Equation

This section derives the instantaneous semi-major axis decay rate from the drag acceleration using the Gauss variational equation and introduces the ballistic coefficient $B^*$.

*Stub.* State the drag acceleration:

$$\vec{a}_{\rm drag} = -\frac{1}{2}C_D\frac{A}{m}\rho v^2 \hat{v} \tag{22.1}$$

where $\hat{v} = \vec{v}/v$ is the velocity unit vector, $C_D$ is the drag coefficient, $A/m$ is the area-to-mass ratio. Define the ballistic coefficient $B^* = C_D A/(2m)$ (SGP4 convention) and the matched-pair form (Ch 3).

Derive the Gauss variational equation for $\dot{a}$ under a purely tangential perturbing force $f_T$:

$$\dot{a} = \frac{2a^2 v}{\mu}\,f_T \tag{22.2}$$

Substitute $f_T = -B^* \rho v^2$ to obtain the instantaneous $\dot{a}$ due to drag. State that orbit-averaging (Ch 12) converts this to a secular mean-motion change.

**Definition 22.2.1** (Ballistic coefficient $B^*$). *Stub — formal definition in SGP4 convention.*

**Notation table:**

| Symbol | Meaning |
|--------|---------|
| $B^*$ | SGP4 ballistic coefficient $= C_D A/(2m)$ |
| $\tau$ | Power-law density exponent |
| $\rho_0$ | Reference density |
| $q_0 - s$ | Reference altitude above shifted surface |
| $\xi$ | Intermediate quantity: $(q_0 - s)/a$ |
| $\eta$ | $= a\sqrt{1-e^2}/(q_0 - s)$ (dimensionless semi-latus rectum) |
| $\psi$ | $= (q_0 - s)^2/(a^2(1-e^2))$; $\psi = \xi^2/\eta^2$ |
| $C_1$ | COEF1: leading factor in $\dot{n}$ expansion |
| $C_2$ | Principal secular mean-motion decay coefficient |
| $C_3$ | $J_3$-eccentricity coupling coefficient |
| $C_4$ | Eccentricity decay coefficient |
| $C_5$ | Mean anomaly correction coefficient |
| $D_2, D_3, D_4$ | Time-polynomial coefficients for $B^*$ time variation |

---

### §22.3 The Orbit-Averaged Drag Integral

This section sets up the orbit-averaged drag integral by substituting the power-law density and the Keplerian orbit geometry, reducing the problem to a definite integral over eccentric anomaly.

*Stub.* The secular mean-motion rate from drag is:

$$\left\langle\frac{dn}{dt}\right\rangle = \frac{1}{2\pi}\int_0^{2\pi} \frac{dn}{dE}\,dE \tag{22.3}$$

Substitute the density model $\rho = \rho_0[(q_0-s)/(r-s)]^\tau$ and $r = a(1-e\cos E)$, approximate $r - s \approx r$ (valid when $r \gg s$) to get:

$$\rho(r) \approx \rho_0\left(\frac{q_0-s}{a}\right)^\tau \frac{1}{(1-e\cos E)^\tau} \tag{22.4}$$

The velocity is $v^2 = \mu(2/r - 1/a)$. Expand in powers of $e\cos E$ using the binomial series (Ch 5, §5.6).

**Theorem 22.3.1** (Orbit-averaged drag integral). *Stub — exact form before Lane expansion.* — *Proof approach: substitute the power-law density and $r = a(1 - e\cos E)$ into the Gauss variational equation for $\dot{a}$, then orbit-average by integrating over $E \in [0, 2\pi]$ using $dt = (r/na) dE$.*

**Error Note [A.22.1]:** Approximating $r - s \approx r$ introduces fractional error $s/r \approx 6371/7000 \approx 0.9$ for LEO altitudes. This is the motivation for using $s$ as a fitting parameter absorbed into $\rho_0$ rather than a physical surface offset. [Tier III for $s \ne 0$; Tier I when $s = 0$.]

---

### §22.4 The Lane Combinatorial Expansion

This section applies the Lane (1965) combinatorial method to reduce the orbit-averaged drag integral to a sum of binomial coefficients and Wallis integrals.

*Stub.* The integrand contains $(1 - e\cos E)^{-(\tau + k)}$ for various integer $k$. Lane (1965) expands these as binomial series in $e\cos E$:

$$(1 - e\cos E)^{-n} = \sum_{j=0}^\infty \binom{n+j-1}{j} e^j \cos^j E \tag{22.5}$$

After term-by-term orbit averaging (exchanging sum and integral), each $\cos^j E$ integral becomes a Wallis integral $W_j$ (Ch 6). The result is an infinite series in powers of $e^2$ multiplied by binomial coefficients times Wallis values.

**Theorem 22.4.1** (Lane expansion). *Stub — express the orbit-averaged density integral as a double sum over binomial coefficients and Wallis integrals.* — *Proof approach: expand $(1 - e\cos E)^{-(\tau+k)}$ by the binomial series (Ch 5, generalized binomial expansion), interchange sum and integral, then evaluate each $\int_0^{2\pi}\cos^j E\,dE/(2\pi)$ using Wallis integrals (Ch 6, Wallis recurrence and closed forms).*

**Proposition 22.4.1** (Convergence). *Stub — the Lane series converges for $e < 1$ with geometric-tail bound from Ch 5, §5.3.* — *Proof approach: apply the Ch 5 geometric-tail bound to the binomial series with $x = e\cos E$; the dominant tail term at order $j$ is $O(e^j)$, so the series converges absolutely for $e < 1$ with geometric-tail remainder.*

**Error Note [P.22.1]:** Truncating the Lane series at finite order introduces precision error with a geometric-tail bound. SGP4 truncates at effective $O(e^4)$ for $C_2$; this bound is $O(e^6)$.

---

### §22.5 Evaluation at $\tau = 4$: The AFGP4 Result

This section specializes the Lane expansion to $\tau = 4$ and derives the full AFGP4 expression for $C_2$ through $O(e^4)$.

*Stub.* Set $\tau = 4$. Carry out the Lane expansion to obtain the full AFGP4 expression for $C_2$. Retain all terms through $O(e^4)$ (or higher, clearly labeled). Define the intermediate quantities $\xi$, $\eta$, $\psi$, COEF1 and state their role.

The full AFGP4 $C_2$ through $O(e^4)$ is:

$$C_2 = {\rm COEF1} \cdot (n_0\xi^{-7/2})\left[\text{(Lane terms in } e^2, e^4)\right] + (\text{J}_2\text{-coupling term}) \tag{22.6}$$

*Stub — state the complete algebraic expression with all coefficients derived, not matched.*

**Theorem 22.5.1** (AFGP4 $C_2$ coefficient). *Stub — full expression at $\tau = 4$.* — *Proof approach: substitute $\tau = 4$ into Theorem 22.4.1, evaluate all binomial coefficients and Wallis integral values in closed form (Ch 6), collect terms by powers of $e^2$, and group into the COEF1 prefactor and the bracketed polynomial.*

---

### §22.6 The $J_2$-Density Coupling and $C_2$ Decomposition

This section derives the $J_2$-density coupling correction by linearizing the density about the Keplerian radius and orbit-averaging the correction term.

*Stub.* The $J_2$ gravitational perturbation produces a radial displacement $\delta r_{J_2}$ at each point on the orbit (Ch 18, short-period correction). This shifts the density evaluation:

$$\rho(r + \delta r_{J_2}) \approx \rho(r)\left(1 + \tau\frac{\delta r_{J_2}}{r - s}\right) \tag{22.7}$$

Orbit-averaging the correction term produces an additional contribution to $\langle\dot{n}\rangle$. Define:

- **Part A ($C_2^A$):** Keplerian drag — integral with unperturbed density (§22.5).
- **Part B ($C_2^B$):** $J_2$-density coupling — integral of the linearized $J_2$ density correction.

**Theorem 22.6.1** ($C_2$ decomposition). *Stub — $C_2 = C_2^A + C_2^B$ with explicit expressions for both parts.* — *Proof approach: expand density as $\rho(r + \delta r_{J_2}) \approx \rho(r)(1 + \tau\,\delta r_{J_2}/(r-s))$, substitute into the orbit-averaged integral, and identify the unperturbed term $C_2^A$ (from §22.5) and the $J_2$ correction term $C_2^B$ separately.*

---

### §22.7 The $3/8\, J_2$ Coefficient

This section traces the $3/8\,J_2$ factor through three multiplicative steps, establishing that each factor arises from a distinct averaging or potential-definition step.

*Stub.* The $J_2$-density coupling term in $C_2^B$ contains the secular part of the $J_2$ radial correction orbit-averaged with the power-law density. Derive that the relevant factor is:

$$\frac{3}{2}\frac{J_2}{2}\frac{3\cos^2 i - 1}{2} = \frac{3}{8}J_2(3\cos^2 i - 1) \tag{22.8}$$

Trace the $3/8$ from: factor of $3/2$ (orbit average of $P_2(\sin\phi)$ over the orbit); factor of $1/2$ from the $J_2$ potential definition; another factor of $1/2$ from the radial perturbation formula. Each factor must be derived, not asserted.

**Theorem 22.7.1** (Origin of the $3/8\,J_2$ factor). *Stub — algebraic derivation from the orbit average of the $J_2$ radial displacement.* — *Proof approach: compute the orbit-average of $\delta r_{J_2}(E)$ (the $J_2$ radial displacement from Ch 18, short-period correction) by integrating over eccentric anomaly; collect all numerical prefactors step by step to show the net factor is $3/8$ times $J_2(3\cos^2 i - 1)$.*

---

### §22.8 Terms Dropped: AFGP4 to SGP4

This section identifies every term present in the AFGP4 derivation that was dropped in the SGP4 implementation and bounds the magnitude of each truncation.

*Stub.* Identify and bound every term present in the full AFGP4 derivation that was dropped in the SGP4 implementation:

1. $O(e^4)$ and higher terms in $C_2$.
2. The eccentricity dependence of $\eta$ in the COEF1 factor.
3. Sub-leading terms in the $J_2$-density coupling.
4. Atmospheric oblateness correction.

For each dropped term, give a representative magnitude for a typical LEO orbit ($a \approx 7000$ km, $e \approx 0.001$, $i \approx 51.6°$).

**Error Note [A.22.2]:** The AFGP4-to-SGP4 truncation introduces accuracy error of order $O(e^4 C_2)$. For typical LEO orbits ($e \approx 0.001$), this is $\sim 10^{-12}$ per revolution — negligible for propagation over days but potentially relevant for multi-year predictions. [Tier I for $e < 0.01$; Tier II for $e > 0.1$.]

---

### §22.9 Higher-Order Coefficients $C_3$, $C_4$, $C_5$

This section derives the $J_3$ coupling coefficient $C_3$, the eccentricity decay coefficient $C_4$, and the mean anomaly correction $C_5$ from their respective Gauss variational equations.

*Stub.*

**$C_3$:** Arises from the $J_3$ odd-zonal coupling with eccentricity. Derive from the Gauss variational equation for $\dot{e}$ under $J_3$; the orbit-averaged result involves $\sin\omega$.

**$C_4$:** Governs the secular eccentricity decay under drag. Derive from the Gauss variational equation for $\dot{e}$ with the power-law density and the Lane expansion.

**$C_5$:** Mean anomaly correction coefficient. Arises from the periodic component of the drag-induced anomaly shift.

State each coefficient in fully explicit algebraic form with derivation stubs.

**Theorem 22.9.1** ($C_3$ from $J_3$ coupling). *Stub.* — *Proof approach: apply the Gauss variational equation for $\dot{e}$ under the $J_3$ disturbing function; orbit-average and collect the $\sin\omega$ term; identify the coefficient as $C_3$.*

**Theorem 22.9.2** ($C_4$ eccentricity decay coefficient). *Stub.* — *Proof approach: apply the Gauss variational equation for $\dot{e}$ under the power-law drag force; orbit-average using the Lane expansion through $O(e^2)$; extract the secular eccentricity decay rate.*

**Theorem 22.9.3** ($C_5$ mean anomaly correction). *Stub.* — *Proof approach: compute the integral of the tangential drag perturbation over one orbit to obtain the accumulated anomaly shift; the orbit-averaged coefficient of the periodic (non-secular) part is $C_5$.*

---

### §22.10 Time-Polynomial Coefficients $D_2$, $D_3$, $D_4$

This section derives the time-polynomial coefficients by repeated differentiation of the drag equation, establishing that $D_k$ are fully determined by lower-order coefficients.

*Stub.* The ballistic coefficient $B^*(t)$ is treated as constant in the matched-pair model. However, the secular variation of $n$ and $a$ produces a slow variation in the effective $B^*$ when expressed as a time polynomial. Define $D_k$ as the $k$-th order Taylor coefficient in the expansion of $n(t) - n_0$ in powers of $t$:

$$n(t) = n_0 + C_1 t + D_2 t^2 + D_3 t^3 + D_4 t^4 + \cdots \tag{22.9}$$

Derive each $D_k$ by repeated differentiation of the drag equation and substitution of lower-order results. Note that $D_k$ are derived quantities, not free parameters.

**Theorem 22.10.1** (Time-polynomial coefficients). *Stub — derive $D_2$, $D_3$, $D_4$ explicitly.* — *Proof approach: differentiate $\dot{n} = C_1 + 2C_2 n(t)^{2/3}B^*$ with respect to time; substitute $\dot{n}$ back to obtain $\ddot{n}$ in terms of lower-order quantities; continue to third and fourth derivatives; collect coefficients by power of $t$ to read off $D_2$, $D_3$, $D_4$.*

**Error Note [P.22.2]:** The time polynomial is an asymptotic expansion; it is accurate for $t \ll t_{\rm decay}$ where $t_{\rm decay} = n_0/|\dot{n}|$. For orbital lifetimes longer than a few weeks, higher-order terms may be needed. [Tier II for $t > 10$ days.]

---

### §22.11 Generalization

This section extends the drag coefficient derivations to general $\tau$, restored $O(e^2)$ terms, and a co-rotating atmosphere.

*Stub.*

**22.11.1 Restoring $O(e^2)$ terms.** The Lane expansion through $O(e^6)$ recovers two additional terms in $C_2$. For moderate eccentricity ($e \approx 0.1$), these contribute at the $10^{-4}$ level relative to $C_2$.

**22.11.2 General $\tau$.** The Lane expansion for arbitrary integer $\tau \geq 1$ produces $C_2^{(\tau)}$ as a polynomial in $e^2$ with binomial coefficient prefactors. The derivation generalizes straightforwardly; provide the general formula.

**22.11.3 Co-rotating atmosphere.** Replacing $v^2$ with $v_{\rm rel}^2 = v^2 - 2\vec{v}\cdot(\vec{\omega}_E\times\vec{r}) + |\vec{\omega}_E\times\vec{r}|^2$ introduces two additional orbit-averaged terms. The cross term $-2\vec{v}\cdot(\vec{\omega}_E\times\vec{r})$ averaged over the orbit gives a correction proportional to $\omega_E^2 p / \mu$; the squared term gives $O(\omega_E^2 r^2/v^2) \approx 0.4\%$ correction.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | error propagation bounds | Error propagation through the coefficient expressions |
| Ch 3 | matched-pair principle | Matched-pair principle: $B^*$ absorbs density model error |
| Ch 5 | generalized binomial expansion | Binomial series used in Lane expansion |
| Ch 6 | Wallis integrals | Wallis integrals $W_n$ used in orbit averages |
| Ch 12 | orbit-averaging principle | Orbit-averaging principle; Gauss variational equations |
| Ch 13 | zonal harmonic definitions | $J_2$, $J_3$ definitions and sign conventions |
| Ch 16 | first-order secular rates | First-order secular rates (needed to express $p$ in corrections) |
| Ch 18 | short-period $J_2$ radial correction | $J_2$ radial displacement used in $C_2^B$ derivation |
| Ch 21 | power-law density model | Power-law density model $\rho(r)$ |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 33 | secular update | $C_1$, $C_2$, $D_2$–$D_4$ used in mean-motion secular propagation |
| Ch 34 | near-space pipeline | $C_3$, $C_4$, $C_5$ used in near-space eccentricity and anomaly updates |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.22.1] | A | §22.3 | $r - s \approx r$ approximation: error $s/r \approx 91\%$ — must use $s$ as fitted parameter; Tier III |
| [A.22.2] | A | §22.8 | AFGP4-to-SGP4 truncation: $O(e^4 C_2)$; Tier I for $e < 0.01$; Tier II for $e > 0.1$ |
| [A.22.3] | A | §22.11 | Non-rotating atmosphere in drag integral: ~5% in secular $\dot{n}$; inherited from [A.21.3] |
| [P.22.1] | P | §22.4 | Lane series truncation: geometric-tail bound $O(e^{2N+2})$ |
| [P.22.2] | P | §22.10 | Time-polynomial asymptotic expansion valid only for $t \ll t_{\rm decay}$; Tier II for $t > 10$ days |
| [A.22.4] | A | §22.4 | Lane (1965) original drag coefficient derivation contains errors in the $C_1$–$C_5$ formulas corrected by Lane & Hoots (1979); do not use Lane (1965) directly; derive from Lane & Hoots (1979) |

---

## Errata: BH61 Normalization Paradox and INPE-2746 Misattribution

**Date:** 2026-04-05

**Scope:** This errata applies to any chapter referencing BH61 Eq(14) or the INPE "spurious Poisson terms" claim (Ch 12, Ch 17, Ch 22).

### Erratum 1: BH61 Eq(14) structural form verified — cos(2f+2g) terms are present in D(Brouwer's form)

**Error (in project documents `theorem_dp1.md`, `cleanroom_phase_b_results.md`, `cleanroom_phase_a3_results.md`):** These documents concluded that BH61 Eq(14)'s cos(2f+2g) terms are spurious, based on a derivation showing D(∂S₁/∂l) = 6ΓB₀(a/r) with no cos(2f+2g) term.

**Correction:** The derivation applied D to a DIFFERENT function (using Γ = μ²/(a³η³) normalization) than BH61's (using n = μ²/L³ normalization). These two forms of ∂S₁/∂l differ by a velocity-dependent factor -(μ/a)^(3/2), which equals -1 only when a = μ — the normalized test case (μ=1, a=1) used in all prior numerical verification.

- D(our form) = 6ΓB₀(a/r) — **no cos(2f+2g)** (81/81 numerical match at μ=1, a=1)
- D(Brouwer's form) matches BH61 Eq(14) structural form — **WITH cos(2f+2g)** (81/81 numerical match at μ=1, a=1)

Both D computations are numerically verified for their respective input functions. The cos(2f+2g) terms in BH61 Eq(14) arise from D(ρ³) multiplying cos(2f+2g), not from D(cos(2f+2g)). Individual coefficients in BH61 Eq(14) still require re-derivation to verify the transcription.

**Verification:** `verify_normalization_ratio.py`, `verify_D_brouwer_form.py`, `verify_breakthrough_algebra.py`

### Erratum 2: INPE-2746 "spurious Poisson terms" misattributed

**Error:** INPE-2746-PRE/322 (Fitzgibbon, De Moraes, Lobão, 1983) was cited as independent corroboration that BH61 Eq(14) contains spurious terms.

**Correction:** Reading the INPE paper reveals it never mentions Eq(14), the D operator, δp₁, or cos(2f+2g). The "spurious Poisson terms" refer to artifacts of BH61's **Method of Successive Approximations during orbit integration**. This is a well-known issue where successive approximation applied to oscillatory systems generates artificial mixed secular-periodic (Poisson) terms. Fitzgibbon's alternative (Variation of Parameters with Bessel function density expansion) avoids these integration artifacts by construction.

The INPE claim and our normalization issue are about completely different things:
- **INPE:** Integration method artifacts (Method of Successive Approximations)
- **Our issue:** D-operator normalization (Γ vs n prefactor)

**Source:** OCR translation from Portuguese at `C:\Users\graha\Desktop\INPE-2746 en.md` (use for logic only, not symbol-level derivation)

### Note: Fitzgibbon's Variation of Parameters approach (INPE-2746 §4)

INPE-2746 describes an alternative integration framework for the drag-oblateness coupling:

**Method:** Lagrange Variation of Parameters (VoP) using Brouwer (1959) as the unperturbed solution, with the density expanded in modified Bessel functions $I_n(\alpha a e)$ instead of Taylor series in $\alpha a e$.

**Claimed advantages over BH61 (per INPE-2746):**
1. INPE-2746 states the VoP solution "does not contain Poisson terms" that BH61's does. What those terms are is not specified in the conference paper.
2. Bessel function density expansion converges uniformly for all $\alpha a e$ (vs Taylor divergence for $\alpha a e > 1$)

**Limitations:**
1. The full derivation is in Fitzgibbon's Master's thesis (ITA, 1982, ref [5]), not in the 10-page conference paper
2. Requires algebraic manipulator for the trigonometric series
3. Uses exponential density model — not directly comparable to SGP4's power-law model

**Relevance to this chapter:** The Fitzgibbon VoP approach is NOT needed for deriving SGP4's C₁-C₅ coefficients (which come from Lane's power-law theory). However, it provides:
- An independent validation pathway for BH61's orbit-averaged drag rates (if the full thesis can be obtained)
- A better density expansion technique (Bessel functions) if the exponential model is ever needed

**Recommendation:** Do not pursue the Fitzgibbon approach for the SGP4 derivation chain. The Lane (1965) ��� Lane & Hoots (1979) pathway is the correct one for Ch 22. The Fitzgibbon Bessel density expansion is documented in Ch 21 §21.8.3 as an alternative technique for exponential atmospheres.

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 8 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 2 |
| Examples | 2 |
| Error Notes | 6 |
| Equations | ~20 |
| Sections | 11 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §22.1 | Draft | Stub |
| §22.2 | Draft | Stub |
| §22.3 | Draft | Stub |
| §22.4 | Draft | Stub |
| §22.5 | Draft | Stub |
| §22.6 | Draft | Stub |
| §22.7 | Draft | Stub |
| §22.8 | Draft | Stub |
| §22.9 | Draft | Stub |
| §22.10 | Draft | Stub |
| §22.11 | Draft | Stub |
