# Draft Plan: Chapter 14 — The Equipotential Ellipsoid

## Objectives

1. Establish the four defining parameters of the level ellipsoid and the two canonical parameter sets (A and B) used in geodetic practice (§14.2).
2. Derive all geometric constants algebraically from the defining parameters (§14.3).
3. Derive the normal gravitational potential and its shape factors $q_0$, $q_0'$ (§14.4) — including closed forms, power-series equivalents, and cancellation analysis.
4. Establish Brouwer's formula relating $J_2$ to $e^2$ and develop the iterative solution when $J_2$ is the defining parameter (§14.5).
5. Derive the even zonal harmonics $J_{2n}$ as determined quantities of the normal field (§14.6).
6. Derive normal gravity on the ellipsoid surface via the Somigliana formula and its series expansion (§14.7).
7. Identify the canonical geodetic integral form and derive the meridian arc, equal-area radius, and mean gravity via the three-stage coefficient pattern (§14.8).
8. Define and derive the three standard mean radii $R_1$, $R_2$, $R_3$ (§14.9).
9. Classify all derived constants by error tier with complete error-type breakdown (§14.10).

---

## Section Structure

### §14.1 Introduction
**Maturity:** Develop level (narrative, forward-reference table, maturity statement — no proofs).

- **No formal numbered items.**
- Forward-reference table links §14.3–4 → Ch 5, §14.5 → Ch 1 §1.6, §14.6 → Ch 7, §14.7 → Ch 6/8, §14.8 → Ch 5/6.
- Maturity declaration: §§14.3, 14.4, 14.8 are Build level; all others are Develop level.
- Model error framing: every derived constant carries accuracy error [A.14.1] from model–reality gap.

---

### §14.2 Defining Parameters
**Maturity:** Develop level (definitions and remarks complete; no proofs required by nature of the section).

- **Definition 14.2.1** (Defining parameter set A): $(a, 1/f, GM, \omega)$ — used by WGS84, WGS72, GRS1967.
- **Definition 14.2.2** (Defining parameter set B): $(a, J_2, GM, \omega)$ — used by GRS80.
- **Remark** (Non-interchangeability): set A gives $e^2$ algebraically; set B requires iteration to recover $e^2$.
- **Remark** (Historical note on WGS84): original DMA 1987 used $\bar{C}_{2,0}$ (set B); late 1990s revision adopted $1/f$ as exact (set A); $\Delta(1/f) \approx 1.46 \times 10^{-6}$ relative to GRS80, ~0.1 mm in $b$.
- **Table 14.2.1** (Defining parameters for WGS84, WGS72, GRS80): values, parameter set, source.
- **Remark** (Tier I and model error): all defining parameters are Tier I with $\sigma_m = 0$, $\delta_p = 0$; accuracy error [A.14.1] characterised by inter-system evolution.
- **Remark** (Atmospheric mass): $GM$ values include atmosphere ($GM_\text{atm} \approx 3.5 \times 10^8$ m³s⁻²); solid-Earth applications must subtract.

---

### §14.3 Derived Geometric Constants
**Maturity:** Build level (Theorem 14.3.1 and Corollary 14.3.1 carry complete proofs; remainder are definitional).

- **Definition 14.3.1** (First eccentricity): $e^2 = 2f - f^2 = (a^2-b^2)/a^2$.
- **Definition 14.3.2** (Second eccentricity): $e'^2 = e^2/(1-e^2) = (a^2-b^2)/b^2$.
- **Definition 14.3.3** (Linear eccentricity): $E = ae = be' = \sqrt{a^2-b^2}$.
- **Theorem 14.3.1** (Fundamental identity): $1 - e^2\sin^2\Phi = (1-e^2)(1 + e'^2\cos^2\Phi)$. *Complete proof* (algebraic, 4 lines, uses Definition 14.3.2).
- **Corollary 14.3.1** (Power-law form): $(1 - e^2\sin^2\Phi)^\alpha = (1-e^2)^\alpha(1 + e'^2\cos^2\Phi)^\alpha$. *Complete proof* (one line: apply exponent $\alpha$).
- **Remark** (Convention choice): second eccentricity form adopted as canonical per Moritz (1980); reasons: positive binomial argument, $(-1)^k$ absorbed into $\binom{\alpha}{k}$, matches published tables. Alternative first-eccentricity form equivalent via Corollary 14.3.1.
- **Definition 14.3.4** (Polar radius of curvature): $c = a^2/b$.
- **Remark** (Tier): all geometric constants are Tier III (computed from Tier I).

---

### §14.4 The Normal Gravitational Potential
**Maturity:** Build level (Theorems 14.4.1, 14.4.3, 14.4.4 carry complete proofs; Theorem 14.4.2 stated with proof deferred).

#### §14.4.1 The External Potential
- **Definition 14.4.1** (Normal potential): $U = V + V_c$ where $V_c = -\tfrac{1}{2}\omega^2(x^2+y^2)$.
- **Remark**: equipotential condition $U = U_0$ constrains shape to give rise to $q_0$ and $q_0'$.

#### §14.4.2 The Shape Factor $q_0$
- **Definition 14.4.2** (Shape factor $q_0$): closed form $q_0 = \tfrac{1}{2}[(1 + 3/e'^2)\arctan e' - 3/e']$. Sourced to Heiskanen & Moritz (1967) Eq. 2-70; Moritz (1980) p. 130.
- **Theorem 14.4.1** (Series representation of $q_0$): $2q_0 = \sum_{n=1}^{\infty} \frac{4(-1)^{n+1}n}{(2n+1)(2n+3)} e'^{2n+1}$, for $|e'| < 1$. *Complete proof* (5 steps: Maclaurin expansion of $\arctan e'$, multiply by $(1 + 3/e'^2)$, separate $k=0$ term, index shift, algebraic cancellation of $3/e'$, simplify bracket).
- **Corollary 14.4.1** (Cancellation avoidance): closed form loses ~5.4 decimal digits for WGS84 due to $O(1/e')$ subtraction; series avoids this because $3/e'$ cancels algebraically in proof Step 3. *Complete proof* (one sentence: leading term is $O(e'^3)$).
- **Corollary 14.4.2** (Convergence rate): ratio of successive terms is $\frac{(n+1)(2n+1)}{n(2n+5)} \cdot e'^2 \to e'^2$; approximately $0.0058$ at $n=1$ and $0.0067$ asymptotically for WGS84. *Complete proof* (direct computation of ratio).
- **Example 14.4.1**: first four terms of $2q_0$ at WGS84 — $t_1 = 1.475 \times 10^{-4}$, $t_2 = -8.523 \times 10^{-7}$, $t_3 = 4.787 \times 10^{-9}$, $t_4 = -2.737 \times 10^{-11}$; 7 terms suffice for double precision.
- **Remark**: error-bounded evaluation deferred to Ch 5, §5.5.1.

#### §14.4.3 The Derivative $q_0'$
- **Definition 14.4.3** (Shape factor derivative $q_0'$): closed form $q_0' = 3(1 + 1/e'^2)(1 - \arctan(e')/e') - 1$. Sourced to Heiskanen & Moritz (1967) §2-7; Moritz (1980) p. 130.
- **Remark** (Distinguishing from superficially similar expressions): the incorrect form $3[\arctan(e')/e' - 1/(1+e'^2)] - 1 \approx -0.987$ for WGS84; the correct form evaluates to approximately $+0.00269$.
- **Theorem 14.4.2** (Series representation of $q_0'$): $q_0' = \sum_{n=1}^{\infty} \frac{6(-1)^{n+1}}{(2n+1)(2n+3)} e'^{2n}$, for $|e'| < 1$. *Stated; proof given as Theorem 14.4.3.*
- **Theorem 14.4.3** (Equivalence of closed form and series for $q_0'$): closed form equals the series. *Complete proof* (5 steps: Maclaurin expansion of $\arctan(e')/e'$, subtract from 1, multiply by $(1+1/e'^2)$, index shift, combine — constant term vanishes, $n \geq 1$ terms yield the series).
- **Corollary 14.4.3** (Structural comparison with $q_0$): $q_0'$ and $2q_0$ share denominator $(2n+1)(2n+3)$ and convergence ratio $e'^2$; differ in numerator ($6$ vs $4n$) and power parity (even $e'^{2n}$ vs odd $e'^{2n+1}$); $q_0'$ converges slightly faster. *Complete proof* (comparison table of general terms plus ratio analysis).
- **Example 14.4.2**: first three terms of $q_0'$ at WGS84 — $t_1 = 2.696 \times 10^{-3}$, $t_2 = -7.786 \times 10^{-6}$, $t_3 = 2.915 \times 10^{-8}$; first term gives 3 significant digits.
- **Remark**: error-bounded evaluation deferred to Ch 5, §5.5.2.

#### §14.4.4 The Normal Potential $U_0$
- **Definition 14.4.4** (Normal potential on the ellipsoid): $U_0 = (GM/E)\arctan e' + \tfrac{1}{3}\omega^2 a^2$.
- **Theorem 14.4.4** (Series representation of $U_0$): $(GM/E)\arctan e' = (GM/b)\sum_{n=0}^{\infty} \frac{(-1)^n}{2n+1} e'^{2n}$. *Complete proof* (substitute $E = be'$, apply Maclaurin expansion of $\arctan(e')/e'$).
- **Remark**: this is the $\arctan(e')/e'$ series; evaluation with error bounds deferred to Ch 5, Theorem 5.5.3; coefficients $1/(2n+1)$ simpler than $q_0$ and $q_0'$.

---

### §14.5 Brouwer's Formula: $J_2$ from Eccentricity
**Maturity:** Develop level (proof approaches only; Example 14.5.1 complete with numerical table).

- **Theorem 14.5.1** (Brouwer's formula): $J_2 = (e^2/3)(1 - 2me'/(15q_0))$ where $m = \omega^2 a^2 b / GM$. Sourced to Heiskanen & Moritz (1967) Eqs. 2-90, 2-92; Moritz (1980) p. 129. *Proof approach only*: derive from boundary condition $U = U_0$; express external potential in ellipsoidal harmonics; convert to spherical harmonic expansion; identify coefficient of $P_2$.
- **Corollary 14.5.1** (Iterative solution for $e^2$): when $J_2$ is the defining parameter, solve $e^2 = 3J_2 + 2me'e^2/(15q_0)$ iteratively from $e^2 \approx 3J_2$. *Proof approach only*: rearrange to fixed-point form $e^2 = g(e^2)$; contraction ratio $O(m/q_0) \approx O(10^{-2})$; Banach fixed-point theorem (Ch 1, §1.6) guarantees convergence.
- **Example 14.5.1** (Iteration for GRS80): table of 4 iterations converging to $e^2 = 0.006\,694\,380\,022\,90$; contraction ratio ~$0.02$; 3 iterations for double-precision accuracy.
- **Remark** (Model error): $J_2$ changes secularly at $\dot{J}_2 \approx -2.6 \times 10^{-11}$ yr⁻¹; sensitivity $\partial e^2/\partial J_2 \approx 3$; $\Delta J_2 = 10^{-11}$ produces $\Delta e^2 \approx 3 \times 10^{-11}$.
- **Remark**: convergence governed by Ch 1, §1.6.

---

### §14.6 Even Zonal Harmonics $J_{2n}$
**Maturity:** Develop level (proof approach only; Example 14.6.1 complete with numerical table; error note fully elaborated).

- **Theorem 14.6.1** (Even zonal harmonic formula): $J_{2n} = (-1)^{n+1} \frac{3e^{2n}}{(2n+1)(2n+3)}\bigl(1 - n + 5n J_2/e^2\bigr)$ for $n \geq 1$. Sourced to Heiskanen & Moritz (1967) Eq. 2-92; Moritz (1980) p. 130. *Proof approach only*: expand ellipsoidal potential in ellipsoidal harmonics; convert to spherical harmonics; match $P_{2n}$ coefficients; denominator $(2n+1)(2n+3)$ traces to Legendre eigenvalues (same pattern as $q_0$, $q_0'$ series and geodetic integral denominators).
- **Remark** (Consistency at $n=1$): substituting $n=1$ into Theorem 14.6.1 gives $J_2 \equiv J_2$ identically — formula self-consistent at lowest order; non-trivial content begins at $n=2$.
- **Example 14.6.1** (Even zonal harmonics for GRS80): table of $J_4$, $J_6$, $J_8$ with comparison to Moritz (1980) p. 131 values; magnitudes decrease by factor ~$e^2 \approx 1/150$ per order.
- **Remark** (Model error [A.14.2]): $J_{2n}$ from Theorem 14.6.1 are for the normal field; observed Earth harmonics (EGM2008) differ due to non-hydrostatic density; anomalous zonal $\Delta J_{2n} = J_{2n}^\text{observed} - J_{2n}^\text{ellipsoid}$ encodes mantle convection, post-glacial rebound, CMB topography; for perturbation theory the observed values enter the disturbing function.

---

### §14.7 The Somigliana Normal Gravity Formula
**Maturity:** Develop level (all proof approaches only; numerical examples complete).

#### §14.7.1 Gravity at the Equator and Pole
- **Theorem 14.7.1** (Equatorial and polar gravity): $\gamma_e = (GM/ab)(1 - m - me'q_0'/(6q_0))$, $\gamma_p = (GM/a^2)(1 + me'q_0'/(3q_0))$. Sourced to Moritz (1980) p. 130 Eqs. (2-73), (2-74); Heiskanen & Moritz (1967) p. 69. *Proof approach only*: differentiate $U$ with respect to outward normal; $q_0'$ arises from differentiating the $\arctan$-based potential; $m$ is centrifugal-to-gravitational ratio.
- **Definition 14.7.2** (Centrifugal ratio): $m = \omega^2 a^2 b / GM$. Note: $m \neq \omega^2 a/\gamma_e$. For GRS80, $m = 0.003\,449\,786\,003\,08$.

#### §14.7.2 The Somigliana Formula
- **Definition 14.7.1** (Somigliana constant): $k = b\gamma_p/(a\gamma_e) - 1$.
- **Theorem 14.7.2** (Somigliana normal gravity formula): $\gamma = \gamma_e(1 + k\sin^2\Phi)/\sqrt{1 - e^2\sin^2\Phi}$. Sourced to Moritz (1980) p. 131. *Proof approach only*: gradient of $U$ on ellipsoid in ellipsoidal coordinates; numerator linear in $\sin^2\Phi$; denominator is metric coefficient $a\sqrt{1-e^2\sin^2\Phi}$; evaluate at poles to pin $\gamma_e$, $\gamma_p$; Somigliana (1929) original derivation.

#### §14.7.3 The Gravity Series
- **Theorem 14.7.3** (Gravity series expansion): $\gamma = \gamma_e(1 + \sum_{n=1}^{\infty} a_{2n}\sin^{2n}\Phi)$ where $a_2 = \tfrac{1}{2}e^2 + k$, $a_4 = -\tfrac{3}{8}e^4 + \tfrac{1}{2}e^2 k$, $a_6 = \tfrac{5}{16}e^6 + \tfrac{3}{8}e^4 k$, $a_8 = \tfrac{35}{128}e^8 + \tfrac{5}{16}e^6 k$. *Proof approach only*: expand $(1-e^2\sin^2\Phi)^{-1/2}$ via generalized binomial series (Ch 5, Theorem 5.3.1) with $\alpha=-1/2$; multiply by $(1+k\sin^2\Phi)$; coefficient of $\sin^{2n}\Phi$ is $a_{2n} = \binom{-1/2}{n}(-e^2)^n + k\binom{-1/2}{n-1}(-e^2)^{n-1}$.
- **Remark** (Clairaut's theorem): consistency check $f + f^* = (\omega^2 b/\gamma_e)(1 + e'q_0'/(2q_0))$; GRS80 values: $f + f^* = 0.003\,352\,811 + 0.005\,302\,440 = 0.008\,655\,251$. Sourced to Heiskanen & Moritz (1967) Eq. 2-75; Moritz (1980) p. 131.
- **Example 14.7.1** (GRS80 and WGS84 gravity parameters): two tables — GRS80 values of $m$, $\gamma_e$, $\gamma_p$, $k$, $f^*$; WGS84 values of $\gamma_e$, $\gamma_p$, $k$, $e^2$; GRS80–WGS84 difference $\Delta\gamma_e \approx 1.4 \times 10^{-5}$ m s⁻² arising entirely from $\Delta(1/f)$.
- **Remark** (Model error [A.14.1]): $\gamma$ is the model gravity; true $g = \gamma + \delta g$ where $\delta g$ is gravity anomaly ($\pm 10$–$300$ mGal); Somigliana formula exact for model; model–reality gap is the accuracy error; for orbit propagation the normal gravity is the reference field for perturbations.

---

### §14.8 Geodetic Integrals
**Maturity:** Build level for §14.8.1 (definitions and remarks complete with full proofs where applicable); §14.8.2 and §14.8.3 Build level (series stated fully); §14.8.4 Develop level (definition and remark only, no derivation).

#### §14.8.1 The Canonical Form
- **Definition 14.8.1** (Geodetic integral): $I(\alpha, w) = \int_0^{\pi/2} (1 + e'^2\cos^2\Phi)^\alpha w(\Phi)\,\mathrm{d}\Phi$ where $\alpha$ is a half-integer and $w(\Phi)$ is a latitude weight. Conditions: $e'^2 > 0$, convergence absolute since $e'^2\cos^2\Phi < 1$.
- **Remark** (Convention): second eccentricity form consistent with §14.3 convention; no $(-1)^k$ in binomial argument; alternating signs in published tables come from $\binom{\alpha}{k}$ alone for $\alpha < 0$.
- **Remark** (Relationship to alternative form): $(1-e^2\sin^2\Phi)^\alpha = (1-e^2)^\alpha(1+e'^2\cos^2\Phi)^\alpha$ via Corollary 14.3.1; first-eccentricity form has explicit $(-1)^k$, all-positive combined coefficients, and separate prefactor $(1-e^2)^\alpha$.
- **Table 14.8.1** (Exponents and applications): $\alpha = -3/2$ (meridian arc, $w=1$); $\alpha = -2$ (equal-area sphere, $w=\cos\Phi$); $\alpha = -5/2$ (mean gravity numerator, $w=\cos\Phi$); $\alpha = -2$ (mean gravity denominator, $w=\cos\Phi$); $\alpha = -1/2$ (Somigliana pointwise).

#### §14.8.2 The Meridian Arc
- **Definition 14.8.2** (Meridian quadrant): $Q = c\int_0^{\pi/2}(1 + e'^2\cos^2\Phi)^{-3/2}\,\mathrm{d}\Phi$.
- **Theorem 14.8.1** (Meridian arc series): $Q = c(\pi/2)(1 - \tfrac{3}{4}e'^2 + \tfrac{45}{64}e'^4 - \tfrac{175}{256}e'^6 + \tfrac{11025}{16384}e'^8 + \cdots)$.
- **Theorem 14.8.2** (Three-stage coefficient derivation): coefficients of $e'^{2k}$ in any geodetic integral arise from (1) binomial expansion giving $\binom{\alpha}{k}$, (2) Wallis-type integration of $\cos^{2k}\Phi \cdot w(\Phi)$ over $[0,\pi/2]$, (3) series algebra for composed quantities. All coefficients are exact rational numbers at each stage.

#### §14.8.3 Equal-Area Radius
- **Definition 14.8.3** (Equal-area radius $R_2$): $R_2 = c(\int_0^{\pi/2}\cos\Phi\,(1+e'^2\cos^2\Phi)^{-2}\,\mathrm{d}\Phi)^{1/2}$.
- **Theorem 14.8.3** (Equal-area radius series): $R_2 = c(1 - \tfrac{2}{3}e'^2 + \tfrac{26}{45}e'^4 - \tfrac{100}{189}e'^6 + \tfrac{7034}{14175}e'^8 + \cdots)$.

#### §14.8.4 Mean Gravity
- **Definition 14.8.4** (Mean gravity): average of normal gravity over the ellipsoid surface.
- **Remark**: involves quotient of two geodetic integrals; requires series quotient (stage 3 of Theorem 14.8.2). No further derivation present.

---

### §14.9 Mean Radii
**Maturity:** Develop level ($R_1$ and $R_3$ have algebraic proof approaches; $R_2$ deferred to §14.8.3; example and tier remarks complete).

- **Definition 14.9.1** (Arithmetic mean radius): $R_1 = (2a+b)/3 = a(1-f/3)$. *Proof approach*: arithmetic mean of semi-axes $(a,a,b)$; substitute $b = a(1-f)$; algebraic identity, no precision error.
- **Definition 14.9.2** (Equal-volume sphere radius): $R_3 = (a^2 b)^{1/3}$. *Proof approach*: equate ellipsoid volume $\tfrac{4}{3}\pi a^2 b$ to sphere volume $\tfrac{4}{3}\pi R_3^3$; cancel; cube root. Purely algebraic.
- **Example 14.9.1** (Mean radii for GRS80): table of $R_1 = 6\,371\,008.7714$ m, $R_2 = 6\,371\,007.1810$ m, $R_3 = 6\,371\,000.7900$ m; ordering $R_3 < R_2 < R_1 < a$; all three radii differ by less than 8 m.
- **Remark** (Tier classification): $R_1$, $R_3$ are Tier II (algebraic, $\delta_p = 0$); $R_2$ is Tier III (geodetic integral with series truncation error).
- **Remark** (Model error): mean radii describe the reference ellipsoid, not the real Earth; geoid undulations $N$ up to ~100 m; mean geoid undulation $\bar{N} \approx -30$ m in current models; $R_1$, $R_3$ have no model error beyond inherited from defining parameters.

---

### §14.10 Tier Classification Summary
**Maturity:** Develop level (classification table and remarks complete; no proofs required).

- **Table 14.10.1** (Tier classification): full table covering all derived constants — Tier I (defining parameters), Tier II (algebraic constants), Tier II with [P.14.1] ($e^2$ from $J_2$ iteration), Tier III ($q_0$, $q_0'$, $U_0$, $J_{2n}$ for $n \geq 2$, $\gamma_e$, $\gamma_p$, $k$, $m$, $f^*$, $a_{2n}$, $R_2$, $Q$, $\bar{\gamma}$), Tier II ($R_1$, $R_3$).
- **Remark**: $\sigma_m = 0$ for all rows — within a geodetic system, defining parameters are exact specifications; inter-system differences (e.g., $\Delta a = 2$ m between WGS72 and WGS84) evidence accuracy error.
- **Remark**: Tier II algebraic constants have $\delta_p = 0$; $e^2$ obtained iteratively from $J_2$ carries [P.14.1] and must be iterated to convergence.

---

### Error Notes Section (standalone, after §14.10)
**Maturity:** Develop level (placeholder notes plus full elaboration of accuracy errors).

- **[M.14.1]** (Placeholder): measurement error propagation. Within a geodetic system $\sigma_m = 0$; inter-system differences reveal accuracy evolution, each system internally self-consistent.
- **[P.14.1]** (Placeholder): precision error for $e^2$-from-$J_2$ iteration. Contraction ratio ~$0.02$; $|e^2_k - e^2_*| \leq (0.02)^k |e^2_0 - e^2_*|$; 3 iterations for double-precision.
- **[P.14.2]** (Placeholder): precision error in series evaluation of $q_0$, $q_0'$, $U_0$ (developed in Ch 5). All alternating with convergence ratio $e'^2 \approx 0.0067$; Leibniz bound applies; 7 terms suffice for double precision.
- **[A.14.1]** (Model error — equipotential ellipsoid vs. real Earth): five failure modes fully elaborated — (1) hydrostatic equilibrium (observed $f = 1/298.257$ vs. hydrostatic $f_h \approx 1/299.5$, ~0.4% excess flattening from post-glacial rebound and mantle convection); (2) time-invariant parameters ($\dot{J}_2 \approx -2.6 \times 10^{-11}$ yr⁻¹ from LAGEOS; $\dot{\omega}/\omega \approx -5 \times 10^{-22}$ s⁻¹ from tidal friction); (3) ellipsoidal shape (geoid undulations up to ~100 m); (4) normal gravity field (free-air anomalies ~$-300$ to $+300$ mGal); (5) atmospheric mass ($GM_\text{atm} \approx 3.5 \times 10^8$ m³s⁻²). Perturbation theory framing: model errors define the reference; perturbation terms measure deviations.
- **[A.14.2]** (Ellipsoidal vs. observed zonal harmonics): $J_{2n}$ from Theorem 14.6.1 are the normal-field harmonics; observed harmonics (EGM2008) differ; $\Delta J_{2n}$ encodes non-hydrostatic density structure; for perturbation theory the observed $J_{2n}$ enter the disturbing function.

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1, §1.6 (contraction mapping) | §14.5 | Banach fixed-point for $J_2 \to e^2$ iteration |
| Ch 5, §5.5.1–5.5.3 | §14.4 | Series evaluation of $q_0$, $q_0'$, $\arctan(e')/e'$ |
| Ch 5, Thm 5.3.1 | §14.7.3 | Generalized binomial series with $\alpha = -1/2$ |
| Ch 5, §5.6 | §14.8 | Geodetic binomial series framework |
| Ch 6 | §14.7–14.8 | Wallis integrals and binomial coefficients |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 7 | §14.6 | Zonal harmonic values $J_{2n}$ for perturbation theory |
| Ch 8 | §14.7 | Normal gravity for atmospheric drag corrections |
| Ch 13 | §14.6 | Reference ellipsoid $J_{2n}$ for disturbing function |
| Ch 16 | §14.2–14.3 | Defining parameters and geometric constants |
| App A | §14.2 | WGS72/WGS84 parameter tables |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.14.1] | M | §14.10 | Measurement error: $\sigma_m = 0$ within a system; inter-system differences reveal accuracy evolution |
| [P.14.1] | P | §14.5 | $e^2$-from-$J_2$ iteration: contraction ratio ~0.02; 3 iterations for double precision |
| [P.14.2] | P | §14.4 | Series evaluation of $q_0$, $q_0'$, $U_0$: convergence ratio $e'^2 \approx 0.0067$; 7 terms suffice |
| [A.14.1] | A | §14.2 | Equipotential model vs. real Earth: 5 failure modes (hydrostatic, time-varying, shape, gravity, atmosphere) |
| [A.14.2] | A | §14.6 | Ellipsoidal vs. observed $J_{2n}$: $\Delta J_{2n}$ encodes non-hydrostatic density structure |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §14.1 | Develop | Narrative only |
| §14.2 | Develop | Definitional; no proofs needed |
| §14.3 | Build | Theorem 14.3.1 and Corollary 14.3.1 fully proved |
| §14.4.1 | Develop | Definitional |
| §14.4.2 | Build | Theorem 14.4.1, Corollaries 14.4.1–2 fully proved |
| §14.4.3 | Build | Theorem 14.4.3 fully proved |
| §14.4.4 | Build | Theorem 14.4.4 fully proved |
| §14.5 | Develop | Proof approaches only |
| §14.6 | Develop | Proof approach only |
| §14.7 | Develop | Proof approaches only |
| §14.8 | Build | Series and three-stage theorem stated |
| §14.9 | Develop | Proof approaches for $R_1$, $R_3$ |
| §14.10 | Develop | Tier table complete |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 18 |
| Theorems | 13 |
| Lemmas | 0 |
| Corollaries | 5 |
| Propositions | 0 |
| Examples | 7 |
| Error Notes | 5 |
| Equations | ~10 |
| Sections | 10 |
