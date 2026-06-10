# Chapter 3: The Matched Pair Principle

**Part I: Mathematical Foundations**

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| **[MP]** | Matched Pair annotation: value is exact within the coupled TLE+SGP4 system (Ch 1, §1.2) | §3.2, Def. 3.2.1 |
| $\tau$ | Tolerance parameter: an accuracy threshold in physical units (km, km/s, rad) | §3.4, Def. 3.4.1 |
| $\tau_{\mathrm{standard}}$ | The standard tolerance: chosen so that the result lands within the SGP4 accuracy floor | §3.4, Def. 3.4.2 |
| $\mathcal{M}$ | A propagation model: the complete specification of constants, truncations, and formulas | §3.2, Def. 3.2.1 |
| $\mathbf{e}_{\mathcal{M}}$ | An element set fitted to model $\mathcal{M}$ by least-squares observation fitting | §3.2, Def. 3.2.1 |
| $\boldsymbol{\beta}_{\mathcal{M}}$ | Compensating bias vector: the shift in fitted elements that absorbs model deficiencies | §3.2, Thm. 3.2.1 |
| $\sigma_m$, $\delta_p$, $\delta_a$ | Measurement, precision, accuracy error bounds (Ch 1, Defs. 1.2.1–1.2.3) | Ch 1 |
| $R_N$ | Remainder of a series truncated after $N$ terms (Ch 1, §1.5) | Ch 1 |

---

## §3.1 Introduction

Orbit propagation transforms a set of orbital elements into a predicted position and velocity. But the orbital elements themselves are not raw observations — they are parameters fitted by a least-squares process to a specific propagation model, with its specific physical constants, truncation orders, and computational formulas. The elements encode the model's biases: when the model uses a gravitational parameter that differs from the true value by $\Delta\mu$, the fitted semi-major axis shifts by an amount that compensates for $\Delta\mu$, so that the predicted orbit still matches the observations.

This compensation is not an accident. It is a structural consequence of least-squares fitting. The fitted parameters minimize the observation residuals under the model's assumptions. Change those assumptions — by substituting a "better" constant, extending a truncated series, or refining a polynomial approximation — and the previously fitted parameters no longer minimize the residuals under the new model. The compensating biases, which were helpful under the old model, become net errors under the new one.

This chapter develops the mathematical framework for this constraint. The central results are:

- **The Matched Pair** (§3.2): TLE elements and the SGP4 model form a coupled system. Substituting a "better" component without re-fitting can degrade accuracy (Theorem 3.2.1).
- **LegacyModelValue** (§3.3): a new value category for constants that are "wrong" by modern standards but correct within the coupled system. Their measurement error is zero by convention (Theorem 3.3.1).
- **The Tolerance Parameter** (§3.4): every generalized computation accepts an accuracy threshold $\tau$ at initialization. At $\tau_{\mathrm{standard}}$, the result agrees with the SGP4 computation of Hoots and Roehrich (1980) to within the model's accuracy floor. At tighter $\tau$, additional terms are computed automatically.
- **Two Operating Modes** (§3.5): standard mode preserves the matched pair; enhanced mode uses modern constants and tighter tolerances, with independently validated accuracy.

The results of this chapter are not algebraic derivations. They are constraints derived from the structure of the least-squares fitting process and the definition of the SGP4 model. The "proofs" are arguments from structure, not from equations. Every downstream chapter must respect these constraints.

This chapter resolves the forward references from Chapter 1 — the **[MP]** annotation in the Tier I classification (§1.2, Remark after Definition 1.2.5) and the tolerance parameter preview (§1.5, Remark after Theorem 1.5.3) — and from Chapter 2 — the TEME frame accuracy note [A.2.1] and the GMST measurement note [M.2.1].

---

## §3.2 The Matched Pair

**Definition 3.2.1** (Matched pair)**.** *A **propagation model** $\mathcal{M}$ is the complete specification of:*

1. *Physical constants $\{c_i\}$ (gravitational parameter, zonal harmonics, Earth rotation rate, etc.).*
2. *Truncation orders $\{N_j\}$ (perturbation series order, polynomial degree, iteration count, etc.).*
3. *Computational formulas $\{f_k\}$ (secular rate expressions, periodic correction formulas, coordinate transformations, etc.).*

*An **element set** $\mathbf{e}_{\mathcal{M}}$ is a set of orbital elements (mean motion, eccentricity, inclination, RAAN, argument of perigee, mean anomaly) determined by least-squares minimization of observation residuals against $\mathcal{M}$:*

$$
\mathbf{e}_{\mathcal{M}} = \arg\min_{\mathbf{e}} \sum_{k} \|\mathbf{r}_{\mathrm{obs},k} - \mathbf{r}_{\mathcal{M}}(\mathbf{e}, t_k)\|^2 \tag{3.1}
$$

*where $\mathbf{r}_{\mathrm{obs},k}$ are observation vectors and $\mathbf{r}_{\mathcal{M}}(\mathbf{e}, t)$ is the position predicted by $\mathcal{M}$ at time $t$ with elements $\mathbf{e}$. The pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ is a **matched pair**. The elements $\mathbf{e}_{\mathcal{M}}$ are meaningful only when propagated by $\mathcal{M}$.*

**Remark.** The **[MP]** annotation introduced in Ch 1 (§1.2, Remark after Definition 1.2.5) is formally defined by this definition: a constant or value carries the **[MP]** tag if and only if it is a component of $\mathcal{M}$ in a matched pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$.

**Remark.** The minimization (3.1) is schematic. The actual TLE fitting process (performed by the 18th Space Defense Squadron) involves differential correction with weighting, outlier rejection, and multi-arc constraints. The essential point is that the fitted elements $\mathbf{e}_{\mathcal{M}}$ depend on the model $\mathcal{M}$: change $\mathcal{M}$, and the same observations produce different optimal elements.

**Theorem 3.2.1** (Degradation by substitution)**.** *Let $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ be a matched pair. Let $\mathcal{M}'$ be a modified model obtained by replacing one or more components of $\mathcal{M}$ with "improved" versions — more accurate constants, higher-order truncations, or more precise formulas. If the original elements $\mathbf{e}_{\mathcal{M}}$ are propagated by $\mathcal{M}'$ without re-fitting, the total propagation error can increase:*

$$
\delta_a(\mathcal{M}', \mathbf{e}_{\mathcal{M}}) \geq \delta_a(\mathcal{M}, \mathbf{e}_{\mathcal{M}}) \tag{3.2}
$$

*in general, even though $\mathcal{M}'$ is objectively a better model in the sense that*

$$
\delta_a(\mathcal{M}', \mathbf{e}_{\mathcal{M}'}) \leq \delta_a(\mathcal{M}, \mathbf{e}_{\mathcal{M}}) \tag{3.3}
$$

*when the elements are re-fitted to $\mathcal{M}'$.*

*Proof.* Consider a one-dimensional analogy that captures the essential mechanism. Suppose the model predicts position as $r = f(\theta; c)$, where $\theta$ is an orbital parameter (element) and $c$ is a model constant. The true position is $r_{\mathrm{true}} = f(\theta_{\mathrm{true}}; c_{\mathrm{true}})$. The model uses an inexact constant $c_0 \neq c_{\mathrm{true}}$, so the model prediction is $f(\theta; c_0)$.

The least-squares fit finds the element value $\hat{\theta}$ that minimizes $|r_{\mathrm{true}} - f(\hat{\theta}; c_0)|^2$. Linearizing $f$ about the true values:

$$
f(\theta; c) \approx f(\theta_{\mathrm{true}}; c_{\mathrm{true}}) + \frac{\partial f}{\partial \theta}\,(\theta - \theta_{\mathrm{true}}) + \frac{\partial f}{\partial c}\,(c - c_{\mathrm{true}}). \tag{3.4}
$$

Setting $c = c_0$ and minimizing over $\theta$:

$$
0 = \frac{\partial}{\partial \theta} |r_{\mathrm{true}} - f(\theta; c_0)|^2 \implies \hat{\theta} - \theta_{\mathrm{true}} = -\frac{\partial f / \partial c}{\partial f / \partial \theta}\,(c_0 - c_{\mathrm{true}}). \tag{3.5}
$$

The compensating bias is $\beta = \hat{\theta} - \theta_{\mathrm{true}}$. This bias is proportional to the constant error $c_0 - c_{\mathrm{true}}$ and is precisely the amount needed to cancel the model deficiency. At the fitted element, the residual is:

$$
r_{\mathrm{true}} - f(\hat{\theta}; c_0) = O(|c_0 - c_{\mathrm{true}}|^2) \tag{3.6}
$$

where the leading-order terms cancel by construction of $\hat{\theta}$.

Now substitute the "better" constant $c_{\mathrm{true}}$ while keeping $\hat{\theta}$:

$$
r_{\mathrm{true}} - f(\hat{\theta}; c_{\mathrm{true}}) \approx \frac{\partial f}{\partial \theta}\,\beta = -\frac{\partial f}{\partial c}\,(c_0 - c_{\mathrm{true}}). \tag{3.7}
$$

The residual is first-order in $|c_0 - c_{\mathrm{true}}|$ — larger than the second-order residual (3.6) of the matched pair. The "improvement" in the constant has degraded the prediction because the compensating bias $\beta$ is still present in the elements.

The argument generalizes to multiple elements and multiple constants. The fitted element vector $\hat{\mathbf{e}}$ absorbs a compensating bias vector $\boldsymbol{\beta}_{\mathcal{M}}$ that is a function of all model deficiencies jointly. Correcting some deficiencies (by substitution) while leaving $\boldsymbol{\beta}_{\mathcal{M}}$ unchanged creates a mismatch between the remaining compensation and the partially corrected model.

This is an existence result: degradation occurs when the substituted component is one that the fitting process has compensated for. Substitutions that do not interact with the fitted elements through the least-squares normal equations — for example, fixing an implementation bug that produced incorrect output for given inputs — do not trigger this mechanism. ∎

The following three examples instantiate Theorem 3.2.1 with specific numerical values from the SGP4 codebase.

---

**Example 3.2.1** (Gravitational parameter mismatch)**.**

The SGP4 "wgs84" gravity model uses:

$$
\mu_{\mathrm{SGP4}} = 398600.5 \text{ km}^3/\text{s}^2 \quad \textbf{[MP]} \tag{3.8}
$$

This is $GM_{\mathrm{GPSNAV}}$, the GPS navigation value retained for receiver compatibility ([NGA] Table 3.4, Eq 3-12; `model_selector.h` line 335: `TV::defined("398600.5")`). The refined WGS84 value is:

$$
\mu_{\mathrm{WGS84}} = 398600.4418 \pm 0.008 \text{ km}^3/\text{s}^2 \tag{3.9}
$$

([NGA] Table 3.1, Eq 3-4; `model_selector.h` line 345: `TV::measured("398600.4418", "0.008")`). The difference is $\Delta\mu = 0.0582$ km$^3$/s$^2$.

For a circular orbit with mean motion $n$, the orbital radius is $r = (\mu/n^2)^{1/3}$. Differentiating with respect to $\mu$ at fixed $n$:

$$
\frac{\partial r}{\partial \mu} = \frac{1}{3}\frac{r}{\mu}. \tag{3.10}
$$

The radial bias from using $\mu_{\mathrm{SGP4}}$ instead of $\mu_{\mathrm{WGS84}}$ is:

$$
\Delta r = \frac{1}{3}\frac{r}{\mu}\,\Delta\mu. \tag{3.11}
$$

| Orbit | $r$ (km) | $\Delta r$ (m) | Context |
|-------|----------|-----------------|---------|
| LEO (400 km) | $6778$ | $0.33$ | ISS, most operational satellites |
| MEO (20200 km) | $26560$ | $1.29$ | GPS constellation |
| GEO (35786 km) | $42164$ | $2.05$ | Geostationary satellites |

For the GPS orbit, $\Delta r \approx 1.3$ m. This is the radial bias documented in [NGA] Sec 3.7.1, which was corrected in the Operational Control Segment (OCS) orbit estimation process in 1994 but retained in the broadcast navigation message.

Within the matched pair: the TLE fitting process used $\mu_{\mathrm{SGP4}}$, so the fitted mean motion $n$ absorbed the $\Delta\mu$ bias. The fitted $n$ is slightly different from what it would be if $\mu_{\mathrm{WGS84}}$ had been used. Substituting $\mu_{\mathrm{WGS84}}$ without re-fitting shifts the predicted radius by $\Delta r$ in the wrong direction — because the elements already compensate for the "wrong" $\mu$, producing a double correction (Theorem 3.2.1, Equation 3.7).

---

**Example 3.2.2** (G-function polynomial approximations)**.**

The SGP4 deep-space module evaluates Hansen coefficients $G_{lpq}(e)$ using cubic polynomial fits due to Hough. These are approximations — the exact Hansen coefficients are expressible as Bessel function series (Ch 15). The rest of the deep-space perturbation theory (secular and long-period terms from lunar and solar gravitational effects) was derived assuming these specific polynomial approximations.

Denote the polynomial approximation error $\epsilon_G(e) = G_{lpq}^{\mathrm{Hough}}(e) - G_{lpq}^{\mathrm{exact}}(e)$ and the perturbation theory truncation error $\epsilon_P$ (from omitted higher-order coupling terms). Within the matched pair, the fitting process adjusts the orbital elements to minimize the combined effect. The fitted elements carry a bias that compensates for $\epsilon_G + \epsilon_P$ jointly.

Replacing the Hough polynomials with exact Bessel-series Hansen coefficients eliminates $\epsilon_G$ but disrupts the compensation. The perturbation theory still carries $\epsilon_P$, and the fitted elements still carry the joint compensating bias. The result is an overcompensation of magnitude $\sim |\epsilon_G|$.

For low eccentricity ($e < 0.2$, covering most operational satellites), the Hough polynomial fits are accurate to $\sim 10^{-4}$ relative to the exact Hansen coefficients, making the matched-pair effect negligible compared to the $\sim 1$ km accuracy floor. For high eccentricity ($e > 0.6$, GTO/Molniya orbits), the polynomial approximation error grows and the error cancellation between $\epsilon_G$ and $\epsilon_P$ becomes more significant. The full analysis of $\epsilon_G(e)$ is developed in Ch 15, where the exact Hansen coefficients are derived and the polynomial approximation error is bounded as a function of eccentricity.

---

**Example 3.2.3** (Brouwer secular rate truncation)**.**

The Brouwer perturbation theory (Ch 16–17) computes secular rates of the orbital elements — the long-term drift in $\Omega$ (RAAN), $\omega$ (argument of perigee), and $M_0$ (mean anomaly) — as functions of $J_2$, $J_3$, $J_4$, eccentricity $e$, and inclination $i$. SGP4 uses the Brouwer theory truncated at order $J_2^2 + J_4$, omitting terms of order $J_2^3$, $J_2 J_3$, and higher.

The omitted terms have magnitude $\sim J_2^3 / J_2^2 = J_2 \approx 10^{-3}$ relative to the retained terms (Ch 1, [A.1.3]). For a 400 km LEO orbit, this contributes $\delta_a \sim 10$–$100$ m to along-track position after one day of propagation.

The TLE fitting process, using the truncated secular rates, produces mean elements $\mathbf{e}_{\mathcal{M}}$ whose mean motion $n$, eccentricity $e$, and inclination $i$ absorb the effect of the omitted higher-order secular terms. The fitted mean motion includes a bias $\delta n$ that compensates for the missing $J_2^3$ contribution to the mean motion rate.

Restoring the $J_2^3$ terms (using the full Brouwer series from Ch 17) while keeping $\mathbf{e}_{\mathcal{M}}$ produces a double-counting effect: the elements carry $\delta n$ (which compensated for the omitted terms), and the extended model also computes those terms directly. The net effect is a secular drift:

$$
\Delta_{\mathrm{along-track}}(t) \approx 2\,r\,\delta n \cdot t \tag{3.12}
$$

in along-track position, where $r$ is the orbital radius and $t$ is the propagation time. This drift grows linearly with $t$, eventually exceeding the original truncation error for propagation intervals beyond $\sim 1$ day. The degradation is worse for longer propagation intervals — precisely where accuracy matters most.

---

## §3.3 LegacyModelValue

The three-error framework of Chapter 1 classifies every constant into one of four tiers based on which error categories are nonzero (Definition 1.2.5). The matched pair introduces a complication: the gravitational parameter $\mu = 398600.5$ km$^3$/s$^2$ is not the best available measurement ($398600.4418 \pm 0.008$), yet within the matched pair it must be treated as exact. It is not a DefinedValue (no standard defines it to be exact), not a MeasuredValue (we are not using it as a measurement), and not a ModelValue (it is not derived from a model computation). It is a specification: a value that is what it is because the model says so.

**Definition 3.3.1** (LegacyModelValue)**.** *A **LegacyModelValue** is a constant $c$ that satisfies all of the following:*

1. *$c$ is a component of a matched-pair model $\mathcal{M}$ (Definition 3.2.1).*
2. *$c$ differs from the best current measurement of the corresponding physical quantity.*
3. *$c$ is retained because TLE element sets $\mathbf{e}_{\mathcal{M}}$ were fitted assuming this specific value.*

*A LegacyModelValue is annotated with **[MP]** in all documentation and carries Tier I classification within the matched-pair system.*

**Remark.** The key distinction from DefinedValue: a DefinedValue (e.g., $a_E = 6378.135$ km in WGS-72) is exact because the geodetic standard defines it to be so. A LegacyModelValue (e.g., $\mu = 398600.5$ km$^3$/s$^2$) is treated as exact within the coupled system even though a better measurement exists. The "exactness" is a consequence of the matching (Definition 3.2.1), not of the geodetic definition.

**Theorem 3.3.1** (Error classification of LegacyModelValue)**.** *For a LegacyModelValue $c$ with specified value $c_0$ in the matched-pair model $\mathcal{M}$:*

*(i)* $\sigma_m(c) = 0$.

*(ii)* $\delta_p(c) = |c_0 - \mathrm{fl}(c_0)| \leq |c_0| \cdot \epsilon_{\mathrm{mach}}$.

*(iii)* $\delta_a(c) = \delta_a(\mathcal{M})$.

*The discrepancy between the specified value and the true physical value — i.e., $|c_0 - c_{\mathrm{true}}|$ — is absorbed into the model accuracy $\delta_a(\mathcal{M})$, not into $\sigma_m$.*

*Proof.* Part (i): By Definition 3.3.1, $c_0$ is a component of the model specification $\mathcal{M}$, not a measurement of a physical quantity. The concept of measurement uncertainty ($\sigma_m$, Definition 1.2.1) applies to values obtained by physical measurement processes — sensor readings, observation campaigns, laboratory experiments. A specification has no measurement process and therefore no measurement uncertainty. This is analogous to the SI definition of the speed of light: $c = 299792458$ m/s has $\sigma_m = 0$ not because the measurement is perfect, but because the SI defines the metre in terms of $c$. Within the matched pair, $c_0$ is the value that the model specifies; it is not an approximation to a measurement.

Part (ii): The value $c_0$ is exact within the model specification. The only source of precision error is the representation of $c_0$ in finite-precision arithmetic. By Ch 1, Definition 1.2.2, the representation error satisfies $|c_0 - \mathrm{fl}(c_0)| \leq \frac{1}{2}\,\mathrm{ulp}(c_0) \leq |c_0| \cdot \epsilon_{\mathrm{mach}}$. This is the standard Tier I precision error.

Part (iii): The model accuracy $\delta_a(\mathcal{M})$ is conventionally estimated by comparing SGP4 predictions (using $\mathcal{M}$ with its specific constants, including $c_0$) against high-fidelity numerical integration or dense tracking data. These comparisons inherently include the effect of using $c_0$ instead of $c_{\mathrm{true}}$: the numerical integrator uses the best available constants, while SGP4 uses $c_0$. The resulting accuracy estimate $\delta_a(\mathcal{M}) \approx 1$ km already incorporates the impact of all constant mismatches.

Attributing the discrepancy $|c_0 - c_{\mathrm{true}}|$ separately to $\sigma_m(c)$ would double-count: the effect of using $c_0$ instead of $c_{\mathrm{true}}$ would appear once in $\sigma_m$ (as a "measurement error" of the constant) and again in $\delta_a(\mathcal{M})$ (as part of the model's overall inaccuracy, which was estimated with $c_0$ in place). The correct classification assigns the full discrepancy to $\delta_a(\mathcal{M})$ and sets $\sigma_m = 0$, avoiding the double-count. ∎

**Remark.** Theorem 3.3.1 does not claim that the discrepancy $|c_0 - c_{\mathrm{true}}|$ is negligible. For $\mu$, it is $0.0582$ km$^3$/s$^2$, producing a $\sim 1.3$ m radial effect at GPS altitude (Example 3.2.1). The theorem claims that this discrepancy is correctly classified as part of $\delta_a(\mathcal{M})$, not as $\sigma_m$. The classification determines how the error propagates through downstream computations: $\sigma_m = 0$ means the constant contributes no measurement uncertainty to any quantity derived from it.

**Remark** (Contrast with value categories)**.** The four value categories and their error signatures:

| Category | $\sigma_m$ | $\delta_p$ | $\delta_a$ | Example |
|----------|-----------|-----------|-----------|---------|
| DefinedValue | $0$ | repr. | $0$ | $a_E = 6378.135$ km (WGS-72) |
| MeasuredValue | $> 0$ | repr. | $0$ | $\mu = 398600.4418 \pm 0.008$ (WGS-84 refined) |
| ModelValue | inherited | accumulated | $> 0$ | $G_{211}(e) \approx 3.616 + \ldots$ (Hough polynomial) |
| LegacyModelValue | $0$ | repr. | $\delta_a(\mathcal{M})$ | $\mu = 398600.5$ **[MP]** (SGP4 WGS-84) |

The LegacyModelValue occupies a unique position: it has $\sigma_m = 0$ like a DefinedValue, but nonzero $\delta_a$ like a ModelValue. This reflects its dual nature — exact within the specification, but part of an approximate model.

**Table 3.3.1** (LegacyModelValues in SGP4)**.**

| Constant | SGP4 value **[MP]** | Modern best | $\|c_0 - c_{\mathrm{true}}\|$ | Radial effect | Source |
|----------|---------------------|-------------|-------------------------------|---------------|--------|
| $\mu$ (km$^3$/s$^2$), WGS-72 | $398600.8$ | $398600.4418 \pm 0.008$ | $0.358$ | $\sim 0.6$ m (LEO) | `model_selector.h:325` |
| $\mu$ (km$^3$/s$^2$), WGS-84 | $398600.5$ | $398600.4418 \pm 0.008$ | $0.058$ | $\sim 1.3$ m (GPS) | `model_selector.h:335` |
| $J_2$, WGS-72 | $0.001082616$ | $0.00108262998905$ | $\sim 10^{-9}$ | $< 0.01$ m | `model_selector.h:318` |
| $J_3$ | $-2.53881 \times 10^{-6}$ | $-2.53215 \times 10^{-6}$ | $6.7 \times 10^{-9}$ | $< 0.1$ m | `model_selector.h:299` |
| $J_4$ | $-1.65597 \times 10^{-6}$ | $-1.61990 \times 10^{-6}$ | $3.6 \times 10^{-8}$ | $< 0.1$ m | `model_selector.h:300` |
| $k_e$ | derived from $\mu$ | — | — | — | derived |
| GMST coefficients | IAU 1982 | IERS 2010 | varies by term | $< 1$ m (few-day arc) | Ch 29 |
| $q_0$, $s$ (drag params) | model-specific | empirical | N/A | drag-dependent | Ch 21 |

**Remark.** The Earth rotation rate $\omega_E = 7.292115 \times 10^{-5}$ rad/s is used identically in both SGP4 and modern references at the stated precision, so it does not qualify as a LegacyModelValue under Definition 3.3.1 (condition 2 is not satisfied — there is no discrepancy at the relevant precision). It is a DefinedValue within the WGS system.

**Remark.** For the WGS-84/SGP4 mode, $J_2 = 0.00108262998905$ is specified directly as a model input (`model_selector.h` line 334: `TV::defined("0.00108262998905")`), bypassing the derivation from $1/f$ that the WGS-84 precise mode uses (`equipotential_ellipsoid.h` lines 27–32). This makes the SGP4 $J_2$ a LegacyModelValue: it is a specification, not derived from the defining parameters.

---

## §3.4 The Tolerance Parameter Principle

Chapter 1 established that series truncation (§1.5, Theorems 1.5.2–1.5.3) and iterative convergence (§1.6, Theorem 1.6.1) both produce precision errors with computable bounds. The key question for the matched pair is: how many terms does the generalized computation retain? If it retains the same number as the Hoots and Roehrich (1980) implementation, it reproduces that result. If it retains more, it produces a more precise result that may violate the matched pair.

The tolerance parameter provides the answer: the algorithm does not count terms — it monitors its own remainder bound and stops when the bound drops below a caller-specified threshold $\tau$. At $\tau_{\mathrm{standard}}$, this produces a result within the SGP4 accuracy floor. At tighter $\tau$, more terms are computed automatically.

**Definition 3.4.1** (Tolerance parameter)**.** *A **tolerance parameter** $\tau > 0$ is an accuracy threshold expressed in the physical units of the quantity being computed (km for position, km/s for velocity, rad for angles). Every generalized computation in the propagation pipeline is constructed by an **initializer** that accepts $\tau$ (along with any constants that remain fixed for the propagation arc). The initializer captures $\tau$ in a closure and produces a callable $\lambda$ whose parameter signature contains only the quantities that change per propagation step — orbital elements, time, epoch, etc. The algorithm within $\lambda$ internally iterates or extends until its remainder bound (Theorem 1.5.2 or 1.5.3 for series; Theorem 1.6.1 for iterative methods) is below the captured $\tau$.*

*The caller of $\lambda$ never sees $\tau$. It was fixed at initialization.*

**Remark.** The separation of concerns is both mathematical and architectural:

- **Mathematical:** $\tau$ controls the accuracy/precision tradeoff. At large $\tau$ (loose tolerance), fewer terms suffice, and the result agrees with the SR3 truncation. At small $\tau$ (tight tolerance), more terms are computed, improving precision beyond the SR3 level.
- **Architectural:** The propagation loop calls $\lambda$ with time-varying arguments (orbital elements, epoch). It does not know or care how many terms $\lambda$ uses internally. The tolerance was decided once, at initialization, and is invisible to the caller. This is the closure/lambda pattern developed in Ch 36.

**Example 3.4.1** (Closure architecture)**.**

Consider a series $S(e) = \sum_{k=0}^{\infty} a_k(e)$ with geometric remainder bound $|R_N| \leq |a_{N+1}|/(1-r)$ (Theorem 1.5.3). The initializer receives $\tau$ and the ratio bound $r$, and produces:

$$
\lambda(e) = S_N(e) \quad \text{where } N = \min\{n : |a_{n+1}(e)|/(1-r) < \tau\}. \tag{3.13}
$$

The propagation loop calls $\lambda(e)$ at each time step with the current eccentricity. The value of $\tau$ is invisible — it was captured by the initializer and determines $N$ internally. At $\tau_{\mathrm{standard}}$, the series truncates at the same index as the SR3 implementation. At tighter $\tau$, additional terms are included automatically.

**Definition 3.4.2** (Standard tolerance)**.** *The **standard tolerance** $\tau_{\mathrm{standard}}$ is the value of $\tau$ at which the generalized computation $g$ agrees with the corresponding SR3 computation $f_{\mathrm{SGP4}}$ [Hoots and Roehrich 1980] to within the model's accuracy floor:*

$$
|g(\tau_{\mathrm{standard}}) - f_{\mathrm{SGP4}}| \leq \delta_a(\mathcal{M}). \tag{3.14}
$$

**Remark.** Equation (3.14) does NOT promise bitwise reproduction of SGP4. The change in underlying representation — dual quaternions (Ch 2) and arbitrary-precision arithmetic (Ch 1) — means the output will differ from the SR3 implementation at the precision level. What is preserved is agreement within the $\sim 1$ km accuracy floor that defines the matched pair. This is the correct guarantee: the SGP4 accuracy floor is the dominant error source, so differences below this floor are operationally meaningless.

**Remark** (Not a term count)**.** The tolerance $\tau$ is an accuracy threshold, not a term count. The user of $\lambda$ does not need to know the algorithm's convergence rate or the specific number of terms required. The algorithm — whether a truncated series, a continued fraction, or an iterative solver — is responsible for:

1. Computing its own remainder bound at each step.
2. Stopping when the bound drops below $\tau$.
3. Reporting the actual remainder as part of the precision error $\delta_p$.

This is the principle previewed in Ch 1 §1.5 (Remark after Theorem 1.5.3): "The choice of truncation index $N$ is governed by a caller-specified tolerance parameter."

**Proposition 3.4.1** (Tolerance parameter across computational patterns)**.** *The tolerance parameter $\tau$ applies uniformly to the three main computational patterns in the propagation pipeline:*

| Pattern | $\tau$ controls | Stopping criterion | Reference |
|---------|-----------------|-------------------|-----------|
| Series with remainder bounds | Truncation index $N$ | $|R_N| < \tau$ | Thm 1.5.2, 1.5.3 $\to$ Ch 5 |
| Continued fractions / Padé | Number of convergents $n$ | $|C_n - C_{n-1}| < \tau$ | Ch 4 |
| Iterative solvers | Convergence tolerance $\epsilon$ | $|x_k - x_{k-1}| \cdot L/(1-L) < \tau$ | Thm 1.6.1 $\to$ Ch 9 |

*In each case, tightening $\tau$ automatically produces more terms or iterations. At $\tau_{\mathrm{standard}}$, the computation matches the SR3 result to within $\delta_a(\mathcal{M})$.*

*Proof.* For the series case: let $N_0$ be the SR3 truncation point. The remainder $|R_{N_0}|$ is the precision error of the SR3 computation at that truncation. Setting $\tau_{\mathrm{standard}} \geq |R_{N_0}|$ ensures that the generalized computation, which adds terms until $|R_N| < \tau$, stops at $N \leq N_0$ or produces a result within $|R_{N_0}|$ of the SR3 value. Since $|R_{N_0}|$ is a precision-level quantity (far below $\delta_a(\mathcal{M})$), Equation (3.14) is satisfied.

For the iterative solver case: the SR3 Kepler equation solver uses a fixed iteration count (typically 10 Newton iterations). After $k$ iterations, the contraction mapping bound (Theorem 1.6.1) gives $|x_k - x^*| \leq L^k |x_0 - x^*| / (1-L)$. Setting $\tau_{\mathrm{standard}}$ equal to (or larger than) this bound at $k = 10$ ensures the generalized solver converges to the same level.

For the continued fraction case: the argument is analogous. The successive convergent difference $|C_n - C_{n-1}|$ decreases monotonically for well-conditioned continued fractions. Setting $\tau_{\mathrm{standard}}$ above the difference at the convergent that corresponds to the SR3 Taylor truncation ensures compatibility. ∎

---

## §3.5 Two Operating Modes

The matched pair (§3.2), the LegacyModelValue category (§3.3), and the tolerance parameter (§3.4) together define two natural operating modes.

**Definition 3.5.1** (Standard mode)**.** *In **standard mode**, the propagator is initialized with:*

- *All physical constants set to their LegacyModelValues (Definition 3.3.1, Table 3.3.1).*
- *The tolerance parameter set to $\tau_{\mathrm{standard}}$ (Definition 3.4.2).*

*The resulting propagation preserves the matched pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$. The error budget is:*

| Error category | Bound | Source |
|----------------|-------|--------|
| $\sigma_m$ | Inherited from TLE element fitting | Observation campaign |
| $\delta_p$ | $\delta_p(\text{repr.}) + \delta_p(\text{truncation at } \tau_{\mathrm{standard}})$ | Ch 1, §§1.3–1.6 |
| $\delta_a$ | $\delta_a(\mathcal{M}) \approx 1$ km (position) | [A.3.1] |

**Remark** (Representation caveat)**.** The implementation uses dual quaternions (Ch 2) for rotations and translations, and arbitrary-precision arithmetic (Ch 1) for intermediate calculations. These representation changes mean the output will NOT be bitwise identical to the SR3 implementation (which uses $3 \times 3$ rotation matrices and IEEE 754 binary64 arithmetic). However, the differences are at the precision level ($\delta_p$) and are far below the accuracy floor $\delta_a(\mathcal{M}) \approx 1$ km. Standard mode preserves the matched pair in the accuracy sense, not the bit-pattern sense.

**Definition 3.5.2** (Enhanced mode)**.** *In **enhanced mode**, the propagator is initialized with:*

- *Physical constants set to their MeasuredValues with nonzero $\sigma_m$ (e.g., $\mu = 398600.4418 \pm 0.008$ km$^3$/s$^2$).*
- *The tolerance parameter set to $\tau < \tau_{\mathrm{standard}}$ (tighter, computing additional terms).*
- *Re-derived mathematics where available (exact Hansen coefficients via Bessel series, higher-order Brouwer secular rates, etc.).*

*The resulting propagation breaks the matched pair by design. The error budget is:*

| Error category | Bound | Source |
|----------------|-------|--------|
| $\sigma_m$ | Measurement uncertainties of modern constants | [NGA], IERS |
| $\delta_p$ | $\delta_p(\text{repr.}) + \delta_p(\text{truncation at tighter } \tau)$ | Ch 1, §§1.3–1.6 |
| $\delta_a$ | Validated independently | Numerical integration |

*Output will NOT match standard SGP4 test cases — by design.*

**Remark.** The two modes share the same callable interfaces — only the captured constants and tolerance differ. The propagation loop code is identical in both modes. The mode selection happens entirely at initialization, through the choice of which constants and which $\tau$ the initializer captures. This is the injectable lambda architecture developed in Ch 36.

**Remark.** Standard mode is the default. Enhanced mode is opt-in and is intended for:

- Comparison studies: quantifying the cost of the matched-pair constraint.
- Research applications with independently fitted elements.
- Validation of the standard mode against an independent reference.

---

## §3.6 Breaking the Matched Pair Safely

The matched pair is a constraint on how existing TLE elements may be propagated. It is not a permanent limitation on the propagator's capability. This section identifies the conditions under which the constraint can be relaxed.

**Proposition 3.6.1** (Safe modification)**.** *The matched pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ can be safely replaced by a new pair $(\mathcal{M}', \mathbf{e}_{\mathcal{M}'})$ — with a different model $\mathcal{M}'$ using improved constants, extended truncations, or refined formulas — provided that the elements $\mathbf{e}_{\mathcal{M}'}$ are re-fitted to $\mathcal{M}'$ from the same (or comparable) observation data.*

*Proof.* The re-fitting process (3.1) applied with model $\mathcal{M}'$ produces new optimal elements $\mathbf{e}_{\mathcal{M}'}$ that absorb the biases of $\mathcal{M}'$. By Definition 3.2.1, the new pair $(\mathcal{M}', \mathbf{e}_{\mathcal{M}'})$ is a valid matched pair. If $\mathcal{M}'$ has fewer systematic errors than $\mathcal{M}$ (smaller biases to absorb), the fitting residuals are smaller and $\delta_a(\mathcal{M}', \mathbf{e}_{\mathcal{M}'}) \leq \delta_a(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$ (Equation 3.3). ∎

**Remark.** In practice, re-fitting requires access to the observation data and the element generation infrastructure. For the public TLE system, this means the 18th Space Defense Squadron (18 SDS) would need to adopt the enhanced model $\mathcal{M}'$ for their operational fitting process. Individual users cannot re-fit TLE elements — they can only propagate them with the matched model.

**Proposition 3.6.2** (Unsafe modification)**.** *Substituting a single component of $\mathcal{M}$ — one constant, one truncation order, or one formula — while keeping TLE elements $\mathbf{e}_{\mathcal{M}}$ generated by the original model, breaks the matched pair. The resulting propagation error is governed by Theorem 3.2.1: it can increase even though the substituted component is individually "better."*

*Proof.* This is a direct application of Theorem 3.2.1. The modified model $\mathcal{M}'$ differs from $\mathcal{M}$ in at least one component. The elements $\mathbf{e}_{\mathcal{M}}$ were not fitted to $\mathcal{M}'$, so the compensating bias $\boldsymbol{\beta}_{\mathcal{M}}$ is mismatched. By the mechanism of Theorem 3.2.1 (Equations 3.6–3.7), the residual under $(\mathcal{M}', \mathbf{e}_{\mathcal{M}})$ can exceed the residual under $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$. ∎

The following examples resolve the forward references from Chapter 2.

---

**Example 3.6.1** (IAU 2006 nutation — resolves [A.2.1])**.**

The TEME frame (Ch 2, Definition 2.5.2) is defined by the True Equator and Mean Equinox, using a simplified precession/nutation model. This frame differs from the precise GCRS/ITRS frames (as defined by the IAU 2006/2000A precession-nutation theory) by $\sim 0.1$ arcsec ($\sim 5 \times 10^{-7}$ rad). This is the irreducible accuracy error documented in [A.2.1].

Applying the IAU 2006/2000A corrections to the frame transformation would improve the frame accuracy. But TLE elements were generated assuming the simplified TEME definition. The fitted right ascension of the ascending node $\Omega_{\mathrm{node}}$ and inclination $i$ absorb the frame discrepancy — their values compensate for the simplified nutation model.

Using the corrected nutation with the old fitted elements produces a net angular error of $\sim 0.1$ arcsec (the magnitude of the correction itself), because the elements carry a compensating bias that is no longer appropriate. The positional effect at orbital radius $r$ is:

$$
\Delta_{\mathrm{position}} \approx r \cdot \delta\theta \tag{3.15}
$$

| Orbit | $r$ (km) | $\Delta_{\mathrm{position}}$ (m) |
|-------|----------|----------------------------------|
| LEO (400 km) | $6778$ | $\sim 3.4$ |
| GPS (20200 km) | $26560$ | $\sim 13$ |
| GEO (35786 km) | $42164$ | $\sim 21$ |

All values are below the $\sim 1$ km accuracy floor, so the IAU nutation correction would not by itself push the propagation out of the matched-pair tolerance. But the correction introduces error rather than removing it: the position moves farther from truth, not closer, because the fitted elements already compensate for the simplified model (Proposition 3.6.2).

In enhanced mode (§3.5), the IAU 2006/2000A corrections would be applied consistently with re-fitted elements, creating a new matched pair in the GCRS/ITRS frame (Proposition 3.6.1).

---

**Example 3.6.2** (IERS Earth rotation — resolves [M.2.1])**.**

The GMST polynomial (Ch 29) converts from the TEME frame to the PEF (Pseudo Earth-Fixed) frame via:

$$
R_{\mathrm{TEME} \to \mathrm{PEF}} = M_3(-\theta_{\mathrm{GMST}}) \tag{3.16}
$$

where $\theta_{\mathrm{GMST}}$ is computed from a polynomial in Julian centuries whose coefficients date from the IAU 1982 theory. The IERS 2010 Conventions provide updated coefficients.

The measurement error in the GMST coefficients propagates through the frame rotation (Ch 2, Theorem 2.2.11):

$$
\|\delta(\mathbf{r}_{\mathrm{PEF}})\| \leq \|\delta(\mathbf{r}_{\mathrm{TEME}})\| + \|\mathbf{r}\| \cdot |\delta(\theta_{\mathrm{GMST}})| \tag{3.17}
$$

as documented in [M.2.1]. Using the updated IERS values for $\theta_{\mathrm{GMST}}$ while keeping TLE elements fitted to the old polynomial breaks the matched pair in the TEME$\to$PEF transformation. The fitted RAAN $\Omega_{\mathrm{node}}$ absorbs the difference between the old and new sidereal time models — its value includes a compensating bias for the old GMST polynomial.

The IAU 1982 and IERS 2010 GMST polynomials differ primarily in their higher-order terms (quadratic and cubic in Julian centuries from J2000). For TLE propagation intervals of a few days, the difference in $\theta_{\mathrm{GMST}}$ is small ($\sim 10^{-6}$ rad or less), making the positional effect via Equation (3.17) well below the accuracy floor. The principle nevertheless applies: the substitution moves the prediction in the wrong direction (Proposition 3.6.2).

In enhanced mode, a consistent set of GMST coefficients and re-fitted elements would be used together, forming a new matched pair in the updated frame (Proposition 3.6.1).

---

## §3.7 The Matched Pair as a Constraint on Approximation

The preceding sections established the matched pair (§3.2), classified the constants involved (§3.3), and introduced the tolerance parameter that mediates between SR3-compatible and enhanced computation (§3.4). This section draws the consequence for all downstream chapters: the choice of approximation method is constrained by the matched pair.

**Definition 3.7.1** (Matched-pair compatibility)**.** *A generalized computation $g(\tau)$ that accepts a tolerance parameter $\tau$ is **matched-pair compatible** with the corresponding SR3 computation $f_{\mathrm{SGP4}}$ [Hoots and Roehrich 1980] if and only if:*

$$
|g(\tau_{\mathrm{standard}}) - f_{\mathrm{SGP4}}| \leq \delta_a(\mathcal{M}) \tag{3.18}
$$

*at the standard tolerance $\tau_{\mathrm{standard}}$ (Definition 3.4.2). No constraint is placed on $g(\tau)$ for $\tau < \tau_{\mathrm{standard}}$.*

**Remark.** Equation (3.18) constrains $g$ only at the standard tolerance. At tighter tolerances, the generalized computation is free to diverge from the SR3 result — indeed, it should, because tighter tolerance means computing additional terms that the SR3 computation omits. The constraint is one-sided: agree with SR3 at $\tau_{\mathrm{standard}}$; improve beyond it for $\tau < \tau_{\mathrm{standard}}$.

**Remark** (Why this matters)**.** A continued fraction that converges faster than a Taylor series is mathematically superior — it may reach a given accuracy in fewer terms. But the partial convergents of a continued fraction are different numbers from the partial sums of a Taylor series. If the continued fraction's convergent at $\tau_{\mathrm{standard}}$ differs from the Taylor partial sum at the SR3 truncation point by more than $\delta_a(\mathcal{M})$, the continued fraction violates Equation (3.18) even though it is more accurate in absolute terms.

In practice, this constraint is rarely binding. For most SGP4 computations, the SR3 truncation retains 3–5 terms of a well-behaved series. Any reasonable approximation method — Taylor, Padé, continued fraction — will agree with the SR3 truncation to far better than 1 km at the standard tolerance. The constraint is a design discipline: it must be stated, checked for each computation, and documented — but it will almost always be satisfied.

**Corollary 3.7.1** (Downstream verification requirement)**.** *Every generalized computation developed in Chapters 4–35 that replaces an SR3 computation must verify Equation (3.18) as part of its validation. The verification consists of evaluating $g(\tau_{\mathrm{standard}})$ and $f_{\mathrm{SGP4}}$ on a representative set of orbital elements (eccentricities $0$ to $0.9$, inclinations $0°$ to $180°$) and confirming that the difference is below $\delta_a(\mathcal{M})$.*

The following chapters are directly constrained:

| Chapter | Computation | SR3 form | Verification |
|---------|------------|----------|--------------|
| Ch 4 | Continued fractions / Padé | Taylor truncation | Convergent at $\tau_{\mathrm{standard}}$ agrees with partial sum |
| Ch 5 | Series with remainder bounds | Fixed $N$-term truncation | $N$ at $\tau_{\mathrm{standard}}$ matches SR3 truncation |
| Ch 9 | Kepler's equation solver | Fixed-iteration Newton | Residual at $\tau_{\mathrm{standard}}$ matches SR3 convergence |
| Ch 15 | Hansen coefficients $G_{lpq}(e)$ | Hough cubic polynomials | Bessel series at $\tau_{\mathrm{standard}}$ agrees with polynomials |
| Ch 25 | Equation of center | Truncated Fourier series | Extended series at $\tau_{\mathrm{standard}}$ agrees with SR3 |

---

## Error Notes

**[A.3.1]** SGP4 model accuracy floor. The SGP4 model, operating within the matched pair $(\mathcal{M}, \mathbf{e}_{\mathcal{M}})$, has an irreducible accuracy error of $\delta_a \approx 1$ km for position and $\sim 1$ m/s for velocity over prediction intervals of a few days. This is the combined effect of omitted zonal harmonics ([A.1.1]), simplified atmospheric density ([A.1.2]), Brouwer theory truncation ([A.1.3]), and perturbation expansion truncation ([A.1.4]). The accuracy floor defines the matched pair: any computation that agrees with the SR3 result to within this floor is matched-pair compatible (Definition 3.7.1). *Remedy:* enhanced mode (Definition 3.5.2) with re-derived mathematics, modern constants, and independently validated error budget.

**[M.3.1]** LegacyModelValue measurement convention. LegacyModelValues carry $\sigma_m = 0$ by convention (Theorem 3.3.1), even though the physical quantities they approximate have nonzero measurement uncertainty. The discrepancy $|c_0 - c_{\mathrm{true}}|$ is absorbed into $\delta_a(\mathcal{M})$, avoiding double-counting (Theorem 3.3.1, part iii). This is the correct classification within the matched-pair system: $c_0$ is a specification, not a measurement (Definition 3.3.1). *Remedy:* in enhanced mode (Definition 3.5.2), replace LegacyModelValues with MeasuredValues carrying their true $\sigma_m$.
