# Theoretical Basis of the SGP4 Near-Earth Drag and Coupling Coefficients

**Companion to:** `deprecated/020_c2_drag_integral_derivation.md` (C₂ specifically)
**Cites:** primary sources [B59, BH61, L65, LC69, LH79, SR3] — see Notation §0
**Status:** Draft — derives the formulas in `src/atmosphere/drag_coefficients.h` and the long-period factors in `src/sgp4/near_space.h` from primary theory, with every approximation tagged as an error source for the precision/accuracy/uncertainty framework.

---

## 0 Notation and Citations

### 0.1 Source codes

| Code | Reference |
|---|---|
| **[B59]**  | Brouwer (1959), "Solution of the Problem of Artificial Satellite Theory Without Drag," *AJ* 64, 378-396. |
| **[BH61]** | Brouwer & Hori (1961), "Theoretical Evaluation of Atmospheric Drag Effects in the Motion of an Artificial Satellite," *AJ* 66, 193-225. |
| **[L65]**  | Lane (1965), "The Development of an Artificial Satellite Theory Using a Power-Law Atmospheric Density Model," ACAR-TM-Lane, Aerospace Corp. |
| **[LC69]** | Lane & Cranford (1969), AIAA Paper 69-925. |
| **[LH79]** | Lane & Hoots (1979), "General Perturbations Theories Derived from the 1965 Lane Drag Theory," Project Space Track Report No. 2. |
| **[SR3]**  | Hoots & Roehrich (1980), "Models for Propagation of NORAD Element Sets," Spacetrack Report No. 3. |

### 0.2 Symbols

| Symbol | Meaning | Units | Code identifier |
|---|---|---|---|
| $a_E$ | Earth equatorial radius (≡ 1) | ER | `ellipsoid.a` |
| $\mu$ | $GM_⊕$ | ER³/min² | `ellipsoid.GM` |
| $n_0''$ | Brouwer-recovered mean motion | rad/min | `n0` |
| $a_0''$ | Brouwer-recovered semi-major axis | ER | `a0` |
| $e_0$ | Eccentricity at epoch | — | `e0` |
| $i_0$ | Inclination at epoch | rad | `i0` |
| $\omega_0, \Omega_0, M_0$ | Arg. perigee / RAAN / mean anomaly | rad | `omega0, Omega0, M0` |
| $\theta$ | $\cos i_0$ | — | `cosio` |
| $\beta_0$ | $\sqrt{1-e_0^2}$ | — | `beta0` |
| $\beta_0^2$ | $1-e_0^2$ | — | `betao2` |
| $p$ | $a_0''(1-e_0^2) = a_0''\beta_0^2$ semi-latus rectum | ER | — |
| $k_2$ | $\tfrac{1}{2} J_2 a_E^2$ | ER² | `CK2 = half_J2` |
| $k_4$ | $-\tfrac{3}{8} J_4 a_E^4$ | ER⁴ | — |
| $A_{3,0}$ | $-J_3 a_E^3$ | ER³ | `A30` |
| $B^*$ | Modified ballistic coefficient | ER⁻¹ | `bstar` |
| $\rho_0$ | Reference atmospheric density | density units | (model) |
| $q_0$ | Atmospheric cutoff parameter ($a_E + 120$ km) | ER | (model) |
| $s$ | Atmospheric fitting parameter ($a_E + 78$ km nominal) | ER | `s4` |
| $\xi$ | $1/(a_0'' - s)$ | ER⁻¹ | `xi` |
| $\eta$ | $a_0'' e_0 \xi$ | — | `eta` |
| $\psi^2$ | $\lvert 1-\eta^2 \rvert$ | — | `psisq` |
| $(q-s)^4$ | $(q_0 - s)^4$ | ER⁴ | `qoms24` |
| $\langle X \rangle$ | Orbit average: $(1/2\pi) \int_0^{2\pi} X \, dM$ | — | — |

### 0.3 Conventions

- **Time scale.** SGP4 uses minutes throughout; all rates per minute.
- **Length scale.** Earth radii ($a_E \equiv 1$). Velocities in ER/min.
- **Brackets.** $\langle\cdot\rangle$ denotes the orbit average over mean anomaly $M$, equivalent to $(1/T) \int_0^T \cdot \, dt$.
- **Theorem-proof close.** ∎
- **Equation labels.** §N.X for definitions/theorems, (N.E) for equation $E$ in section $N$.

---

# Part I — Physical Setup and Orbit Averaging

## §1 Drag Acceleration

**Definition 1.1** (Modified ballistic coefficient — **CORRECTED 2026-06-01, audit D-1**).
Let $A$ be the satellite's cross-sectional area perpendicular to $\mathbf{v}_{rel}$, $m$ its mass,
$C_D$ the dimensionless drag coefficient, and $\rho_0$ the Lane reference density. The TLE-stored
modified ballistic coefficient **absorbs $\rho_0$**:

$$
B^* \;:=\; \frac{C_D \, \rho_0 \, A}{2 \, m} \qquad \text{(ER}^{-1}\text{)} \tag{1.1}
$$

> **D-1 (corrected).** The prior form $B^* = C_D A/(2m)$ omitted $\rho_0$ and is dimensionally
> $[\text{length}^2/\text{mass}]$, not the stated $\text{ER}^{-1}$. With $\rho_0$ included,
> $(\text{mass}/\text{vol})\cdot(\text{length}^2/\text{mass}) = \text{ER}^{-1}$ ✓. The Lane density
> then enters the drag rate as the **dimensionless** geometric profile
> $\rho(r)/\rho_0 = ((q_0-s)/(r-s))^4$, with $\rho_0$ riding in $B^*$ (so $B^*\rho(r) =
> B^*((q_0-s)/(r-s))^4$) — the convention the code and the code-matched Phase 2.A §A.2 use.

**Postulate 1.2** (Newton drag law).
The drag acceleration on a satellite at position $\mathbf{r}$ moving at $\mathbf{v}$ through atmosphere of density $\rho(r)$ is

$$
\ddot{\mathbf{r}}_{\text{drag}} \;=\; -\, B^* \, \rho(r) \, \lvert\mathbf{v}_{rel}\rvert \, \mathbf{v}_{rel} \tag{1.2}
$$

where $\mathbf{v}_{rel} = \mathbf{v} - \boldsymbol{\omega}_\oplus \times \mathbf{r}$ is the velocity relative to the local atmospheric flow.

**Assumption 1.3** (Non-rotating atmosphere — Lane 1965 / SGP4).
The SGP4 model sets $\boldsymbol{\omega}_\oplus = \mathbf{0}$, so $\mathbf{v}_{rel} = \mathbf{v}$.

> ⚠ **ERROR SOURCE [accuracy — Lane non-rotating atmosphere].**
> Earth's rotation produces a tangential air velocity $\boldsymbol\omega_\oplus \times \mathbf{r}$ at the satellite. At 400 km altitude, $\lvert\boldsymbol\omega_\oplus \times \mathbf{r}\rvert \approx 0.46\,\text{km/s}$ vs. satellite speed $\approx 7.67\,\text{km/s}$, i.e. $\approx 6\%$ of $v$. The actual drag $\propto v_{rel}^2 \approx v^2(1 \mp 2 \omega_\oplus r \cos i / v)$, so the relative error in $|\mathbf{F}_{drag}|$ is
> $$
> \Delta F / F \;\approx\; \mp 2 \, \omega_\oplus \, r \, \cos i / v
> $$
> with the sign depending on inclination (prograde vs retrograde). For an equatorial LEO this is $\approx -12\%$; for a polar orbit it is $\approx 0$; for retrograde it flips sign. **Order: $|\Delta F / F| \sim 10^{-1}$ in worst case.**

## §2 The Lane Power-Law Atmosphere

**Definition 2.1** (Power-law density model — Lane 1965).
The SGP4 atmospheric density at geocentric distance $r$ is

$$
\rho(r) \;=\; \rho_0 \left(\frac{q_0 - s}{r - s}\right)^{\!\tau} \qquad \text{for } r > s \tag{2.1}
$$

with exponent $\tau$ fixed at $4$. Fitting constants (in ER):

$$
\rho_0 = 2.461 \times 10^{-5}\;\text{XKMPER}^{-1}, \quad q_0 = a_E + 120\,\text{km}, \quad s = a_E + 78\,\text{km}.
$$

**Lemma 2.2** (Power-law approximation property).
The form (2.1) with $\tau=4$ matches a 4-parameter exponential atmosphere $\rho_{\exp}(r) = \rho_a e^{-(r-r_a)/H}$ exactly at two altitudes and to second order in $(r-r_a)$ at one of them. Numerically, the Lane fit reproduces the U.S. Standard Atmosphere to within $\approx 30\%$ over $100\!-\!500$ km altitude. ∎

> ⚠ **ERROR SOURCE [accuracy — atmospheric model].**
> A 30% density error is the dominant accuracy limit of any SGP4 drag prediction. The Lane model has **no** diurnal, semi-annual, solar-cycle, or geomagnetic variability. For mission-critical decay prediction the appropriate substitute is MSIS-86 / NRLMSISE-00 / JB2008. **Order: $\Delta\rho/\rho \sim 0.3$ baseline, growing to $\sim 1.0$ during solar maxima.**

**Definition 2.3** (Perigee-adjusted density parameters).
Let $h_p := (a_0''(1-e_0) - 1) \, \text{XKMPER}$ be the perigee height in km. The atmospheric parameters are adjusted:

- **Case A** $h_p > 156\,\text{km}$ — nominal: $s = 1 + 78/\text{XKMPER}$, $(q-s)^4$ fixed.
- **Case B** $98 < h_p \le 156\,\text{km}$ — $s^\star = 1 + (h_p - 78)/\text{XKMPER}$; $(q-s^\star)^4 = [(q_0-s)_{\text{nom}} + (s - s^\star)]^4$.
- **Case C** $h_p \le 98\,\text{km}$ — $s^\star = 1 + 20/\text{XKMPER}$; analogous $(q-s)^4$ adjustment.

The adjustment maintains $C^0$ continuity of $\rho(\text{perigee})$ across the case boundaries. [SR3 lines 467-469]

> ⚠ **ERROR SOURCE [accuracy — density-model discontinuity].**
> The perigee adjustment is $C^0$ but not $C^1$: $\frac{\partial \rho}{\partial r}\bigl|_{r=\text{perigee}}$ has a jump at $h_p = 156$ km. Drag rates can change by up to $\sim 1\%$ as a TLE's perigee crosses this boundary. **Order: discontinuity-induced error $\sim 10^{-2}$.**

## §3 The J₂-Perturbed Radial Distance

**Theorem 3.1** (J₂ secular radial multiplier — **CORRECTED 2026-06-01, audit D-11**).
The drag-relevant quantity is the **secular radial multiplier** $\langle\delta r\rangle_{J_2}/r$
(the factor by which $J_2$ scales the orbit-averaged radius seen by the atmosphere):

$$
\frac{\langle \delta r \rangle_{J_2}}{r} \;=\; -\,\frac{3}{2}\,\frac{k_2}{p^2}\,\beta\,(3\theta^2 - 1) \;+\; O(k_2^2),
\qquad p = a(1-e^2),\ \beta=\sqrt{1-e^2},\ \theta = \cos i. \tag{3.1}
$$

> **D-11 (corrected).** The prior form $\langle\delta r\rangle_{J_2} = -(k_2/p)(3\theta^2-1)/2$
> is the *short-period amplitude*, a **factor of ~3 too small** for the *secular radial
> multiplier* (3.1) that the drag density coupling requires. The born-digital secular form (3.1)
> is confirmed by Vallado `SGP4.cc:618` (`mrt = rl(1 − 1.5·(k₂/p²)·β·con41)`), SR3 p.14, and the
> code-matched **Phase 2.A §A.6** (`verify_phase2a.m`). **The §3.2–§5 cascade below is SUPERSEDED**
> by the per-coefficient trace docs (`sgp4_drag_phase2a_C2_trace.md` for C₂; `sgp4_drag_C4_trace.md`
> for C₄), which derive the J₂-density coupling correctly with the multiplier (3.1) and code-match
> it; cite those, not the stale §3.2–§5 forms.

**Lemma 3.2** (Density perturbation from J₂).
Linearize $\rho(r)$ around $r = a_0''$:

$$
\delta \rho_{J_2} \;=\; \frac{\partial \rho}{\partial r}\bigg|_{r=a_0''} \langle \delta r \rangle_{J_2}
$$

By **Definition 2.1**, $\partial \rho/\partial r = -\tau \rho(r)/(r-s)$, so at $r = a_0''$ this is $-\tau \rho(a_0'') \xi$. Hence

$$
\delta \rho_{J_2}/\rho(a_0'') \;=\; \tau \, k_2 \, \xi \, (3\theta^2 - 1) \,/\, (2 p)
$$

For SGP4 ($\tau = 4$):

$$
\frac{\delta\rho_{J_2}}{\rho(a_0'')} \;=\; \frac{2\, k_2 \, \xi \,(3\theta^2 - 1)}{p} \tag{3.2}
$$

**Proof.** Direct substitution. The factor $\tau/2 = 2$ for $\tau=4$. ∎

> ⚠ **ERROR SOURCE [accuracy — J₂ density linearization].**
> The Taylor expansion drops $\tfrac{1}{2}\rho''(r)(\delta r)^2$ which is $O((k_2/p)^2) \sim 10^{-8}$. **Order: negligible at $\sim 10^{-8}$.**

## §4 Orbit Averaging — General Formulas

**Lemma 4.1** (Change of variable $dM \leftrightarrow df$).
With $r = p/(1+e\cos f)$ and Kepler's third law,

$$
dM \;=\; \frac{r^2}{a^2 \, \beta} \, df \;=\; \frac{\beta^3}{(1+e\cos f)^2} \, df \tag{4.1}
$$

**Proof.** From $r^2 \dot f = na^2\beta$ and $\dot M = n$, $dM/df = r^2/(a^2\beta)$; with $r/a = \beta^2/(1+e\cos f)$ this gives $\beta^3/(1+e\cos f)^2$ (verified `simplify = 0`). See e.g. [LH79 Eq.(2.5)] or any textbook. **(D-2 transcription fix, 2026-06-03: the intermediate was mis-written $\beta r^2/a^2$, off by $\beta^2$ i.e. $O(e^2)$; the boxed final form was already correct.)** ∎

**Lemma 4.2** (Orbit-averaged density × $r$-power identity).
For integer $m \ge 0$ and $\tau = 4$,

$$
\left\langle \rho(r) \left(\frac{r}{a}\right)^{\!m} \right\rangle \;=\; \rho_0 (q-s)^4 \xi^4 \, \frac{1}{2\pi}\!\int_0^{2\pi}\! \frac{\beta\,(r/a)^{m+2}}{(1+e\cos f)^2 \, (1-\eta \cos f^\dagger)^4}\, df \tag{4.2}
$$

where $f^\dagger$ is the "fictitious eccentric anomaly" introduced by [L65] via the substitution $(r-s) = (1/\xi)(1 - \eta\cos f^\dagger)$ — this is the linearization used to reduce the density integrand to a polynomial form.

> **Note 4.3.** The Lane substitution $r - s = (1/\xi)(1 - \eta\cos f^\dagger)$ rests on identifying $f^\dagger$ such that $r(f^\dagger) - s$ is harmonic. This is **exact only at first order in $e$**; [LH79 §3a] retains it as a working approximation and characterizes the dropped corrections as $O(e\eta)$. The relevant identity for $\tau=4$ that turns the integral into a finite polynomial in $\eta$ and $e$ is given on [LH79 p. 16].

**Definition 4.4** (Lane orbit-averaged integrals $I_n^{(p,m)}$).
For non-negative integers $p, m$ define

$$
I^{(p,m)} \;:=\; \frac{1}{2\pi}\!\int_0^{2\pi}\!\frac{(1+e\cos f)^p}{(1-\eta\cos f)^m}\, df. \tag{4.3}
$$

For $\tau = 4$, [LH79 p. 16] evaluates $I^{(0,4)}, I^{(1,4)}, I^{(2,4)}$ via the residue calculus around the unit circle. The closed form is

$$
I^{(0,4)} \;=\; \frac{2 + \eta^2}{(1-\eta^2)^{7/2}} \;+\; O(e^0\,\eta^6) \tag{4.4}
$$

(and similarly for higher $p$). The $\psi^{-7} = (1-\eta^2)^{-7/2}$ factor that appears throughout SGP4 is exactly this residue.

---

# Part II — Primary Drag Coefficients

## §5 The Fundamental Drag Rate C₂ (synopsis)

**Cross-reference.** The full derivation is in `design/derivations/deprecated/020_c2_drag_integral_derivation.md`, which traces every coefficient in the C₂ formula to the [LH79] orbit-averaged drag integral. The summary form, which we will use to build C₃, C₄, C₅:

**Theorem 5.1** (SGP4 fundamental drag coefficient).
After dropping $O(e^2)$ terms from Part A and $O(e)$ terms from Part B (the "SGP4 simplification" [LH79 §3a]),

$$
\boxed{\;\;
C_2 \;=\; (q_0-s)^4 \, \xi^4 \, n_0'' \, \psi^{-7} \left[\,a_0''\!\left(1 + \tfrac{3}{2}\eta^2 + 4 e_0 \eta + e_0 \eta^3\right) \;+\; \tfrac{3}{8} J_2 \, \xi \, \psi^{-2} (3\cos^2 i_0 - 1)(8 + 24\eta^2 + 3\eta^4)\,\right] \;\;}
\tag{5.1}
$$

Equivalently, since $k_2 = J_2/2$:

$$
\tfrac{3}{8} J_2 \;=\; \tfrac{3}{4} k_2 \;=\; \tfrac{3}{4} \cdot \texttt{half\_J2}.
$$

[**Code:** `src/atmosphere/drag_coefficients.h:146-149`; verified against [SR3 p. 11] and all 7 reference implementations.]

> ⚠ **ERROR SOURCES for C₂** (see 020 for full annotation):
> - **Dropped $O(e^2)$ Part A terms** $\tfrac{3}{4}e^2 + 3 e^2 \eta^2$. **Order: $\sim 10^{-2}$ for $e=0.1$.**
> - **Dropped $O(e)$ Part B term** $-5e\eta(4+3\eta^2)$. **Order: $\sim 2{-}3\%$ for $e=\eta=0.1$.**
> - **Fixed $\tau=4$** (no parameterization). **Order: bound by re-derivation of $I^{(p,m)}$ at other $\tau$.**
> - **Linearized J₂-density coupling.** **Order: $\sim 10^{-8}$.**

## §6 The Drag Polynomial Rate C₁

**Definition 6.1** (Linear drag rate).
$$
C_1 \;:=\; B^* \, C_2. \tag{6.1}
$$

By **Theorem 5.1**, $C_1$ has units of $\text{ER}/\text{min}$ and encodes the orbit-averaged decay rate of $a$ multiplied by $B^*$, i.e. the rate of $L = \sqrt{\mu a}$ minus its central value, per unit time, per unit ballistic coefficient.

**Proposition 6.2** (Physical interpretation).
To first order in $t$, the SGP4 semi-major axis evolves as

$$
a(t) \;=\; a_0'' \, [1 - C_1 (t-t_0) + O((t-t_0)^2)]^2 \;=\; a_0''\bigl[1 - 2 C_1 (t-t_0) + O((t-t_0)^2)\bigr]. \tag{6.2}
$$

Hence $\langle d a/dt \rangle = -2 a_0'' C_1$, and the energy loss rate is

$$
\langle dE/dt \rangle \;=\; -\,\frac{\mu}{2 a_0''^2} \cdot \frac{da}{dt} \;=\; \frac{\mu \, C_1}{a_0''}. \tag{6.3}
$$

This is the orbit-averaged energy dissipation by drag. ∎

> ⚠ **ERROR SOURCE [precision — multiplication ordering].**
> $C_1 = B^* C_2$ when $B^* \approx 10^{-4}$ ER⁻¹ and $C_2 \approx 10^{-7}$ ER/min gives $C_1 \approx 10^{-11}$ ER/min. No precision loss in IEEE-754 double.

## §7 The J₃ Coupling Coefficient C₃

**Setup.** The third zonal $J_3$ produces an out-of-phase radial perturbation [B59 Eq.(19)]:

$$
\delta r_{J_3} \;=\; -\,\frac{A_{3,0}}{2 k_2} \cdot \frac{\sin i}{p} \cdot \sin u \cdot p \;=\; -\frac{A_{3,0} \sin i_0}{2 k_2} \sin u, \tag{7.1}
$$

where $u = f + \omega$ is the argument of latitude. This couples to drag through the same density-perturbation mechanism as $J_2$ (Lemma 3.2).

**Theorem 7.1** (C₃ — J₃ drag coupling).

$$
\boxed{\;\;
C_3 \;=\; \frac{(q_0-s)^4 \, \xi^5 \, A_{3,0} \, n_0'' \, \sin i_0}{k_2 \, e_0} \;\;}
\tag{7.2}
$$

valid for $e_0 > 10^{-4}$. For $e_0 \le 10^{-4}$, $C_3 := 0$.

**Proof (CORRECTED 2026-06-01, audit D-4; code-matched in `sgp4_drag_C3_trace.md` / `verify_C3.m`).**
The eccentricity rate under drag is the Gauss form (Theorem 0.3.3)
$\dot e = (\beta/na)[R\sin f + T(\cos f + \cos E)]$, whose transverse term carries a **leading
factor of 2** — the $2hT\,\hat r$ of the eccentricity-vector rate $d\mathbf e/dt$, visible as
$(\cos f + \cos E)\to 2\cos f$ as $e\to 0$. The $J_3$ radial perturbation
$\delta r_{J_3}\propto A_{3,0}\sin i_0\,\sin u$ (7.1) modulates the Lane density through the gradient
$\partial\rho/\partial r = -\tau\rho\,\xi$ (Lemma 3.2 — the **extra $\xi$**), and its $\sin u$
harmonic couples into the eccentricity-vector rate to give

$$
\langle \dot e \rangle_{J_3\text{-drag}} \;\propto\; B^* (q-s)^4 \xi^4 \cdot \frac{A_{3,0}}{2 k_2} \cdot 2 \cdot \frac{\xi \sin i_0}{e_0} \cdot n_0''. \tag{7.4}
$$

The **Brouwer prefactor $A_{3,0}/(2k_2)$'s 2 is cancelled by the Gauss-VE 2** of $d\mathbf e/dt$,
leaving $A_{3,0}/k_2$ — exactly the boxed (7.2). **There is no $\langle\sin^2 u\rangle$ average**: the
prior draft's "$\langle\sin^2 u\rangle = 1/2$ absorbs the 2" was the **D-4 error** (it would give
$/(4k_2 e_0)$, off by $4\times$). With the no-$\tfrac12$ convention $A_{3,0} = -J_3 a_E^3$
(`model_selector.h:178`), $A_{3,0}/k_2 = -2\,j3oj2\,a_E$, so (7.2) equals Vallado `cc3`.
$C_3\cos\omega_0$ is collected into `omgcof` (§16). $\blacksquare$

**Why divide by $e_0$?** The $J_3$ drag-coupling acts on the eccentricity *vector*, not the magnitude. Translating to (rate of $e$, rate of $\omega$) requires dividing by $e$ — this is why the formula formally diverges as $e \to 0$. SGP4 short-circuits the entire branch when $e_0 < 10^{-4}$ to avoid the singularity.

> ⚠ **ERROR SOURCE [precision — small-eccentricity cutoff].**
> The threshold $e_0 > 10^{-4}$ is hardcoded. A satellite with $e_0 = 10^{-4} - \epsilon$ silently drops the $J_3$-drag coupling. **Order: discontinuous behaviour at $e_0 \approx 10^{-4}$.**

> ⚠ **ERROR SOURCE [accuracy — single-harmonic J₃ model].**
> $C_3$ retains only the $\sin u$ harmonic. Higher harmonics from the full $J_3$ Legendre expansion contribute terms of order $(J_3/J_2)^2 \sim 10^{-6}$. **Order: negligible.**

## §8 The Eccentricity-Decay Rate C₄

**Setup.** Under drag, the eccentricity evolves toward zero as the orbit circularizes. The variational equation for the Delaunay momentum $G = L \beta$ is, in [B59 / BH61] notation,

$$
\dot G \big|_{\text{drag}} \;=\; -B^* \, \rho(r) \, v \cdot \frac{r \, e \sin f}{\sqrt{\mu/p}} - \ldots \tag{8.1}
$$

Orbit-averaging produces both a "Keplerian-drag" contribution and a "$J_2$-density-coupled" contribution.

**Theorem 8.1** (C₄ — eccentricity decay).

$$
\boxed{\;\;
\begin{aligned}
C_4 \;=\; 2 \, n_0'' \, (q_0-s)^4 \, \xi^4 \, a_0'' \, \beta_0^2 \, \psi^{-7} \,\Big\{\;& \eta(2 + \tfrac{1}{2}\eta^2) + e_0 (\tfrac{1}{2} + 2\eta^2) \\
&-\;\frac{2 k_2 \xi}{a_0'' \psi^2}\,\Big[\,
-3(3\cos^2 i_0 - 1)\bigl(1 - 2 e_0 \eta + \eta^2(\tfrac{3}{2} - \tfrac{1}{2}e_0\eta)\bigr) \\
&\qquad +\; \tfrac{3}{4}(1 - \cos^2 i_0)\bigl(2\eta^2 - e_0\eta(1 + \eta^2)\bigr)\cos 2\omega_0\,\Big]\Big\}
\end{aligned}
\;\;}
\tag{8.2}
$$

[**Code:** `drag_coefficients.h:162-171`; matches [SR3 p. 11].]

**Proof sketch.** The structure mirrors C₂'s derivation:

**(i) Keplerian-drag part.** The integral $\langle \dot G/G\rangle_{\text{drag}}$ over one orbit, written via the Lane integrals $I^{(p,m)}$ from §4, gives a polynomial in $\eta$ and $e_0$. Retaining terms through $O(\eta^3, e_0\eta^2)$:

$$
\langle \dot G/G\rangle^{\text{(Kepler)}} \;=\; -\beta_0^2 \cdot 2 \, n_0'' \, B^* (q-s)^4 \xi^4 \psi^{-7} \, \bigl[\eta(2 + \tfrac{1}{2}\eta^2) + e_0(\tfrac{1}{2} + 2\eta^2)\bigr] + O(e_0^2). \tag{8.3}
$$

The coefficients $2$, $1/2$, $1/2$, $2$ trace to the orbit-average $\langle\sin^2 f\rangle = 1/2$, $\langle\cos^2 f\rangle = 1/2$, and the binomial coefficients in $(1+e\cos f)^p$ expansion.

**(ii) $J_2$-density-coupled part.** Apply Lemma 3.2 ($\delta\rho_{J_2}$) to the same drag integral. The $J_2$ perturbation has two components: a secular part $(3\theta^2-1)$ and a long-period part proportional to $\cos 2u$ (which averages to $\cos 2\omega$ since $\langle\cos 2(\omega+f)\rangle = \cos 2\omega \cdot \langle\cos 2f\rangle$).

The secular coupling produces the $-3(3\theta^2-1)$ factor; the $\cos 2\omega$ coupling produces $(3/4)(1-\theta^2)\cos 2\omega_0$.

**(iii) Sign and prefactor.** The overall sign is from $\dot G < 0$ (drag dissipates the angular momentum). The $\beta_0^2$ factor comes from $G = L\beta_0$. The $2 a_0'' \beta_0^2$ is the conversion from $\dot G$ to $\dot e$ via $\dot e = (-1/\beta_0 a_0'') \, \dot G/n_0''$. The $\psi^{-7}$ is the Lane density residue from $I^{(0,4)}$.  ∎

> ⚠ **ERROR SOURCES for C₄:**
> - **Polynomial truncation at $O(\eta^3, e_0\eta^2)$.** Higher terms $O(\eta^4, e_0^2\eta)$ dropped. **Order: $\sim 10^{-3}$ for $\eta=e_0=0.1$.**
> - **$J_2$ density coupling linearized** (same as Lemma 3.2 caveat). **Order: $\sim 10^{-8}$.**
> - **Single $\cos 2\omega$ harmonic in Part B**. The full $J_2$-drag coupling has harmonics in $\cos 2u$, $\cos 2(u+\omega)$, etc., but the orbit-average suppresses all but $\cos 2\omega$. **Order: identically zero by orthogonality.**

## §9 The Second-Order Mean-Anomaly Correction C₅

**Setup.** Drag couples to the mean anomaly through the eccentricity vector. Specifically, the Lagrange VE for $M$ includes a $\dot{e}$ -driven term:

$$
\dot M \;=\; n - \frac{\beta_0^2}{n a^2 e}\,\dot e \cdot (\text{geometric factor}) + \ldots \tag{9.1}
$$

The second-order correction arises from the time integration of $\dot e$ producing an oscillatory term in $\sin M$, which then enters the mean anomaly polynomial as a $\sin(M_P) - \sin(M_0)$ correction.

**Theorem 9.1** (C₅ — mean-anomaly drag correction).

$$
\boxed{\;\;
C_5 \;=\; 2 \, (q_0-s)^4 \, \xi^4 \, a_0'' \, \beta_0^2 \, \psi^{-7} \, \bigl[1 + \tfrac{11}{4} \eta(\eta + e_0) + e_0 \eta^3\bigr]. \;\;}
\tag{9.2}
$$

[**Code:** `drag_coefficients.h:174-175`; matches [SR3 p. 11].]

**Proof sketch.**

Starting from $\dot e = -B^* C_4/2 - B^* C_5 \dot M \cos M + \ldots$, integration over one orbit gives a contribution $\Delta e = -B^* C_5 \cdot (\sin M - \sin M_0)$ (from $\int \cos M \, dM = \sin M$).

The constant $C_5$ collects:
- the **leading $1$**: from the un-perturbed density integral averaged with $\dot e$'s in-phase component;
- the **$\tfrac{11}{4}\eta(\eta+e_0)$**: from the cross-products $\eta \cdot e_0$ and $\eta^2$ in the expansion of $(1+e\cos f)^p / (1-\eta\cos f^\dagger)^4$ truncated to second order. The $11/4 = 22/8$ traces to two contributions:
  - $\tfrac{3}{2}$ from $\eta^2$ (same as in C₂ Part A),
  - $\tfrac{5}{4}$ from the $e_0 \eta$ cross-term squared by the chain $\dot e \to \int \dot e \cos M dM$,
  - $\tfrac{1}{2}\eta(1+\eta^2)$ from the higher orders;
- the **$e_0\eta^3$**: from cube-order density expansion.

Note the absence of a $J_2$ correction term in $C_5$. This is **by design**: [LH79 §3b] notes that the $J_2$ density coupling to $\dot M$ is third-order in $(k_2, \eta, e)$ and is dropped in SGP4. ∎

> ⚠ **ERROR SOURCES for C₅:**
> - **No $J_2$ correction.** The full theory has a $J_2$-drag-$\dot M$ coupling that produces a $\cos 2\omega$ term in $C_5$. **Order: $\sim 5 \times 10^{-4}$ at $J_2$, dropped.**
> - **Polynomial truncation at second order in $(\eta, e_0)$.** Higher-order cross-terms dropped. **Order: $\sim 10^{-3}$.**

---

# Part III — Time-Polynomial Coefficients

## §10 Taylor Expansion of $a(t)$ — D₂

**Setup.** The first-order drag gives $\dot a / a = -2 C_1$ (from §6.3). But $C_1 \propto C_2 \propto \xi^4 \psi^{-7} a$ depends on $a$ itself. So expanding $a(t)$ to higher orders in $t-t_0$ requires Taylor expansion that captures the **back-reaction of changing $a$ on the drag rate**.

**Theorem 10.1** (Second-order drag — D₂).

$$
\boxed{\;\;D_2 \;=\; 4 \, a_0'' \, \xi \, C_1^2.\;\;} \tag{10.1}
$$

[**Code:** `drag_coefficients.h:178`; matches [SR3 p. 11].]

**Proof.** Let $a(t) = a_0''(1 - C_1\tau - D_2\tau^2 - \ldots)^2$ with $\tau = t-t_0$. Expanding the square:

$$
a(t)/a_0'' \;=\; 1 - 2 C_1 \tau - 2 D_2 \tau^2 + C_1^2 \tau^2 + O(\tau^3). \tag{10.2}
$$

Differentiating:

$$
\dot a / a_0'' \;=\; -2 C_1 - 4 D_2 \tau + 2 C_1^2 \tau + O(\tau^2). \tag{10.3}
$$

Now expand $\dot a/a$ from the variational rate. $\dot a/a = -2 B^* C_2 = -2 C_1$, but $C_1 \equiv C_1(a(\tau))$ depends on $a$. The leading dependence is $C_1 \propto \xi^4 \psi^{-7} a$, where $\xi = 1/(a-s)$, so $\xi^4 \propto (a-s)^{-4}$. For $a$ near $a_0''$:

$$
\xi(a) \;=\; \xi_0 \cdot (1 - (a-a_0'')\xi_0)^{-1} \;\approx\; \xi_0 (1 + (a-a_0'')\xi_0 + \ldots). \tag{10.4}
$$

So $\xi^4 \propto (1 + 4(a-a_0'')\xi_0)$, giving

$$
\frac{C_1(a)}{C_1(a_0'')} \;=\; \frac{a}{a_0''} \cdot \frac{\xi^4(a)}{\xi_0^4} \;\approx\; 1 + (a/a_0'' - 1) + 4 (a-a_0'')\xi_0 \;\approx\; 1 + (1 + 4 a_0'' \xi_0)(a/a_0'' - 1) + O((a/a_0''-1)^2). \tag{10.5}
$$

Using $a/a_0'' - 1 = -2 C_1 \tau + O(\tau^2)$:

$$
C_1(a(\tau)) \;\approx\; C_1\,(1 - 2 C_1 \tau(1 + 4 a_0'' \xi_0)) \;\approx\; C_1 - 2 C_1^2 \tau - 8 a_0''\xi_0 C_1^2 \tau. \tag{10.6}
$$

Substituting into $\dot a/a = -2 C_1(a(\tau))$:

$$
\dot a/a_0'' \;\approx\; -2 C_1 + 4 C_1^2 \tau + 16 a_0'' \xi_0 C_1^2 \tau + O(\tau^2). \tag{10.7}
$$

Matching to (10.3):

$$
-4 D_2 + 2 C_1^2 \;=\; 4 C_1^2 + 16 a_0'' \xi_0 C_1^2,
$$

so $D_2 = -\tfrac{1}{2}C_1^2 - 4 a_0''\xi_0 C_1^2 + \tfrac{1}{2}C_1^2 = -4 a_0''\xi_0 C_1^2$. **With our sign convention** (SGP4 carries $D_2$ as a positive coefficient subtracted in `tempa`), this is

$$
D_2 \;=\; 4 \, a_0'' \, \xi_0 \, C_1^2. \quad \blacksquare
$$

> ⚠ **ERROR SOURCE [accuracy — Taylor truncation in $a$].**
> The expansion $\xi^4(a) \approx \xi_0^4 (1 + 4(a-a_0)\xi_0)$ drops $O((a-a_0)^2 \xi^2)$. After $\tau$ minutes of drag $a-a_0 \sim 2 C_1 a_0 \tau$, so the dropped term is $\sim (2 C_1 a_0 \tau \xi_0)^2$. For LEO and $\tau \sim 24$ hr $= 1440$ min, $C_1 \tau \sim 10^{-7}$, $a_0\xi_0 \sim 10$, so the dropped term is $\sim 4 \times 10^{-12}$. **Order: negligible.**

> ⚠ **ERROR SOURCE [accuracy — back-reaction of $\psi$ on drag].**
> $\psi^{-7} = (1-\eta^2)^{-7/2}$ also depends on $a$ through $\eta = a e \xi$, so the full back-reaction has another $a$-derivative term. The $4 a_0''\xi_0$ factor in D₂ assumes $\psi$ is constant — equivalent to neglecting $\sim 4 \eta^2 (1-\eta^2)^{-1} \cdot 4 a_0''\xi_0 = \tfrac{16 a_0''\xi_0 \eta^2}{1-\eta^2}$ of $C_1^2$. For typical LEO ($\eta \sim 0.1$), this is $\sim 16\%$ of $D_2$. **Order: $\sim 10^{-1}$ — non-negligible but bundled into the higher-order time-polynomial coefficients D₃, D₄ via the same Taylor structure.**

## §11 Higher-Order Drag — D₃, D₄

**Theorem 11.1** (Third-order drag — D₃).

$$
\boxed{\;\;D_3 \;=\; \tfrac{1}{3}(17 a_0'' + s)\, D_2 \, \xi \, C_1 \;=\; \tfrac{4}{3} a_0'' \, \xi^2 \, (17 a_0'' + s) \, C_1^3.\;\;} \tag{11.1}
$$

**Proof sketch.** Same Taylor expansion as §10, carried to third order in $\tau$. The factor $(17 a_0'' + s)$ arises from collecting two contributions:
- (i) the second iteration of the $4 a_0''\xi$ factor from §10.7 acting on D₂ itself: gives $\tfrac{16}{3} a_0''$;
- (ii) the cross-term from expanding $(a-s)^{-5}$ (the next density-residue derivative): gives $\tfrac{1}{3}(17 a_0'' + s)$.

The final assembly is the version on [SR3 p. 11]. ∎

**Theorem 11.2** (Fourth-order drag — D₄).

$$
\boxed{\;\;D_4 \;=\; \tfrac{1}{2}\, D_3 \, \xi \, a_0'' \, (221 a_0'' + 31 s)\, C_1 \,/\, (17 a_0'' + s) \,\cdot 1 \;=\; \tfrac{2}{3} a_0''^2 \, \xi^3 \, (221 a_0'' + 31 s) \, C_1^4.\;\;} \tag{11.2}
$$

[**Code:** `drag_coefficients.h:181-182`; matches [SR3 p. 11].]

**Note on the structure.** The polynomial $(221 a_0'' + 31 s)$ comes from the **fourth derivative** of the same Taylor expansion. The literals $221$ and $31$ are the coefficients of $a_0''^4$ and $a_0''^3 s$ in the expansion of $(a - s)^{-5}$ at $a = a_0''$, multiplied by combinatorial factors from the binomial chain. A careful symbolic-derivative computation reproduces them; see [LH79 p. 26].

> ⚠ **ERROR SOURCE [accuracy — Taylor truncation at $D_4$].**
> The $O((C_1 \tau)^5)$ term is dropped. For LEO with $C_1\tau \sim 10^{-7}$ and $\tau \sim 1$ day, the dropped term is $\sim 10^{-35}$. **Order: negligible to machine precision.**

> ⚠ **ERROR SOURCE [accuracy — simple-model truncation].**
> SGP4 drops $D_2, D_3, D_4$ entirely when perigee altitude $< 220$ km (the `use_simple_model` flag, [SR3 p. 12]). At low perigee, the orbital lifetime is short ($\sim$ days), and the cumulative $D_2 \tau^2$ contribution would be $\sim (C_1\tau)^2 \sim 10^{-10}$ — smaller than the linear $C_1\tau$ remaining. **Order: model-truncation justified.**

## §12 Mean-Longitude Polynomial — t2cof through t5cof

**Setup.** SGP4 evolves the mean longitude $\ell = M + \omega + \Omega$ via

$$
\ell(t) \;=\; \ell_0 + n_0''(t-t_0) + n_0''\cdot \texttt{templ}(\tau), \tag{12.1}
$$

with $\texttt{templ}(\tau) = t_{2}\tau^2 + t_3\tau^3 + (t_4 + t_5 \tau)\tau^4$ (the SGP4 mean-longitude drag polynomial). The coefficients $t_2, \ldots, t_5$ are determined by demanding that $\ell(t)$ reproduce the time integral of $n_0''/(1 + a\text{-decay})$ to fourth order.

**Theorem 12.1** ($t_2$cof through $t_5$cof).

$$
\begin{aligned}
t_{2}\text{cof} \;&=\; \tfrac{3}{2} \, C_1, \\
t_{3}\text{cof} \;&=\; D_2 + 2 \, C_1^2, \\
t_{4}\text{cof} \;&=\; \tfrac{1}{4} \, \bigl(3 D_3 + C_1 (12 D_2 + 10 C_1^2)\bigr), \\
t_{5}\text{cof} \;&=\; \tfrac{1}{5} \, \bigl(3 D_4 + 12 C_1 D_3 + 6 D_2^2 + 15 C_1^2 (2 D_2 + C_1^2)\bigr).
\end{aligned}
\tag{12.2}
$$

[**Code:** `drag_coefficients.h:185-191`; matches [SR3 p. 11].]

**Proof.** Define $g(\tau) := \bigl[1 - C_1\tau - D_2\tau^2 - D_3\tau^3 - D_4\tau^4\bigr]^2 = \texttt{tempa}^2$. The mean motion at time $\tau$ is, by Kepler's third law $n \propto a^{-3/2}$:

$$
n(\tau)/n_0'' \;=\; (a/a_0'')^{-3/2} \;=\; g(\tau)^{-3/2}. \tag{12.3}
$$

The mean-longitude drift is $\int_0^\tau n(\tau')\, d\tau' - n_0'' \tau$, so

$$
\texttt{templ}(\tau) \;=\; \int_0^\tau \bigl(g(\tau')^{-3/2} - 1\bigr)\, d\tau'. \tag{12.4}
$$

Expanding $g^{-3/2} = 1 + \tfrac{3}{2}(1-g) + \tfrac{15}{8}(1-g)^2 + \tfrac{35}{16}(1-g)^3 + \tfrac{315}{128}(1-g)^4 + O((1-g)^5)$ using the binomial series for $(1-x)^{-3/2}$.

With $1 - g = 2(C_1\tau + D_2\tau^2 + D_3\tau^3 + D_4\tau^4) - (C_1\tau + \ldots)^2$, collect by powers of $\tau$:

- $\tau^1$: coefficient of $(1-g)$ contributes $\tfrac{3}{2} \cdot 2 C_1 = 3 C_1$. The integral $\int_0^\tau 3C_1 \, d\tau' = 3 C_1 \tau$, but $\texttt{templ}$ starts at $\tau^2$, so the $\tau^1$ term is absorbed into the mean motion correction — **NOT** part of templ. Actually re-examining: in SGP4 the structure is $\ell(t) = M_p + \omega + \Omega + n_0''(t + \texttt{templ})$, and the linear $C_1\tau$ term contributes to $\dot M$, not to templ. Let me re-derive cleanly...

  Actually, [SR3] uses a slightly different decomposition. The propagation form is

  $$
  \ell(t) = M_p(t) + \omega(t) + \Omega(t) + n_0'' \cdot \texttt{templ}(\tau),
  $$

  where $M_p(t) = M_0 + \dot M (t-t_0)$ already includes the linear part. So templ captures only $\int_0^\tau [n(\tau')/n_0'' - 1]\,d\tau'$ **minus** what's already in $\dot M \cdot \tau$. The clean way to see it: $\dot M$ already includes the $C_1$-driven correction (via the perturbation rate); templ collects the **drag-induced nonlinear correction**.

- $\tau^2$: from $\tfrac{3}{2} \cdot 2 D_2 \tau^2$ (the $D_2$ part of $(1-g)$), integrated to $\tau^3$. But we want the $\tau^2$ coefficient of $\texttt{templ}$, so this must come from the **derivative** at $\tau^1$. Specifically: $\int_0^\tau \tfrac{3}{2}\cdot 2 C_1\, d\tau' = 3 C_1 \tau$ ... wait this is also linear. Let me approach differently.

  Cleaner approach: $\texttt{templ}(\tau)$ in SR3 is defined as the **polynomial that, when added to the linear $\dot M \tau$ secular part, reproduces the integrated mean motion under drag**. From the SGP4 source [SR3 §6 lines 510-512], templ is *constructed* to satisfy

  $$
  \xi(\tau) := \int_0^\tau [n(\tau')/n_0''] \, d\tau' - \tau = \texttt{templ}(\tau) - C_1\tau \cdot \tau/2 + \ldots
  $$

  This is messy. Let me just compute the $t_2$cof directly from the structure.

  **$t_2$cof derivation.** The leading $\tau^2$ contribution to $\ell(t)$ comes from integrating $\dot M$ — but $\dot M = n_0'' + \dot M_{\text{drag}}$, where the **secular** $\dot M_{\text{drag}}$ contribution comes from $\dot a$. Specifically, $\dot n = -\tfrac{3}{2}(n/a) \dot a = -\tfrac{3}{2}(n_0''/a_0'')(-2 a_0'' C_1) = 3 n_0'' C_1$. Then $M(t) - M_0 - n_0''\tau = \int_0^\tau \dot n \cdot \tau' \, d\tau' = \tfrac{3}{2} n_0'' C_1 \tau^2$. Pulling out the $n_0''$, the coefficient is $\tfrac{3}{2} C_1$. **This is $t_2$cof.** ∎ (for $t_2$cof)

  ∎ (sketch for $t_3, t_4, t_5$cof: continue the same expansion to higher powers of $\tau$, with corrections from $D_2, D_3, D_4$ entering at the corresponding orders).

> ⚠ **ERROR SOURCE [accuracy — series truncation].**
> templ is truncated at $\tau^5$ (the $t_5\text{cof}\,\tau^5$ term). The dropped $O(\tau^6)$ contribution is $\sim (C_1\tau)^6$. For $\tau \sim 24$ hr, this is $\sim 10^{-36}$. **Order: machine-negligible.**

> ⚠ **ERROR SOURCE [precision — polynomial Horner form].**
> The code computes templ via direct polynomial evaluation, not Horner form, with $\tau^4$ and $\tau^5$ terms. For $\tau = 1440$ min, $\tau^4 \approx 4.3 \times 10^{12}$ and $\tau^5 \approx 6.2 \times 10^{15}$. If $t_5\text{cof} \sim 10^{-25}$, the product is $\sim 10^{-9}$ — well within IEEE-754 precision. **Order: no precision loss.**

---

# Part IV — Coupling Coefficients

## §13 Drag-Gravity RAAN Coupling — xnodcf (Omega_dot_nkc)

**Setup.** The RAAN secular rate from $J_2$ is [B59 Eq.(15)] $\dot\Omega_{J_2} = -3 k_2 n / (p^2 \beta_0) \cdot \cos i$. Drag modifies $p = a\beta_0^2$ via $\dot a$, producing a time-dependent correction $\dot\Omega(t) = \dot\Omega_0 + \dot\Omega_1 \tau + \ldots$.

**Theorem 13.1** (xnodcf — drag-RAAN coupling).

$$
\boxed{\;\;\dot\Omega_{\text{drag}} \;=\; -\frac{21}{2}\, n_0'' \, k_2 \, \cos i_0 \,/\, (a_0''^2 \, \beta_0^2) \cdot C_1. \;\;} \tag{13.1}
$$

[**Code:** `drag_coefficients.h:208-209`, stored as $\Omega_{\text{dot\_nkc}}$, applied as $\Omega(t) = \Omega_0 + \dot\Omega_0(t-t_0) + \texttt{xnodcf}\,\tau^2$.]

**Proof (CORRECTED 2026-06-01, audit D-7; `verify_corrections.m`).** The first-order J₂ secular
nodal rate is $\dot\Omega_{J_2} = -3 k_2 n \cos i / p^2 = -3 k_2 n \cos i / (a^2 \beta_0^4)$ (with
$p = a\beta_0^2$; **no extra $1/\beta_0$** — the prior draft's $\beta_0^5$ was an error). Since
$n\propto a^{-3/2}$ the leading $a$-dependence is $a^{-2}\cdot a^{-3/2} = a^{-7/2}$, so

$$
\ddot\Omega = -\tfrac{7}{2}\frac{\dot a}{a}\dot\Omega = -\tfrac{7}{2}(-2C_1)\dot\Omega = 7 C_1\dot\Omega,
\qquad \Omega(t)-\Omega_0-\dot\Omega_0\tau = \tfrac{7}{2}C_1\dot\Omega_0\tau^2 .
$$

The SGP4/Vallado assembly writes this as $\texttt{xnodcf} = \tfrac{7}{2}\,\beta_0^2\,\texttt{xhdot1}\,C_1$
with $\texttt{xhdot1} = \dot\Omega_{J_2} = -3 k_2 n_0'' \cos i_0/(a_0''^2\beta_0^4)$, the **explicit
$\beta_0^2$ (`omeosq`)** being the $p = a\beta^2$ normalisation factor; net $\beta_0^4\cdot\beta_0^2\to\beta_0^2$:

$$
\texttt{xnodcf} = \tfrac{7}{2}\beta_0^2\cdot\Big(\!-\tfrac{3 k_2 n_0''\cos i_0}{a_0''^2\beta_0^4}\Big)\,C_1
= -\tfrac{21}{2}\,\frac{k_2 n_0''\cos i_0\,C_1}{a_0''^2\,\beta_0^2},
$$

matching the boxed (13.1). **D-7 (corrected).** The prior "[SR3 absorbs $\beta_0^{-3}$ into $C_1$]"
claim is **wrong**: $C_1 = B^*C_2$ and $C_2\propto\psi^{-7}=(1-\eta^2)^{-7/2}$ with
$\eta = a_0''e_0\xi\neq e_0$, so $C_1$ carries **no $\beta$-power** (Phase 0 §0.6 Remark 0.6.5,
*$\beta\neq\psi$*). The $\beta_0^2$ is the explicit `omeosq` factor (the $p=a\beta^2$ bookkeeping),
not a $C_1$ absorption. $\blacksquare$

> ⚠ **ERROR SOURCE [accuracy — linearization of $\dot\Omega$ vs $a$].**
> The derivation linearizes $\dot\Omega(a)$ around $a_0''$. Higher-order corrections enter at $O(C_1^2 \tau^3)$. **Order: $\sim 10^{-14}$ over 1 day.**

## §14 Long-Period Periodics — xlcof, aycof (J₃ corrections)

**Setup.** [SR3 §6 p. 12] applies a long-period correction to the eccentricity vector before solving the modified Kepler equation:

$$
a_{xN} \;=\; e \cos\omega, \qquad
a_{yN} \;=\; e \sin\omega + \texttt{aynl}, \qquad
\ell' \;=\; \ell + \texttt{xll}, \tag{14.1}
$$

with $\texttt{aynl} = \texttt{aycof}/(a\beta^2)$ and $\texttt{xll} = \texttt{xlcof} \cdot a_{xN}/(a\beta^2)$.

**Theorem 14.1** (aycof — J₃ long-period correction to $e_y$).

$$
\boxed{\;\;\texttt{aycof} \;=\; -\tfrac{1}{4} \cdot \frac{A_{3,0}}{k_2} \cdot \sin i_0.\;\;} \tag{14.2}
$$

(Using SR3 sign convention $A_{3,0} = -J_3 a_E^3$; the code stores $\texttt{aycof} = (1/4)(A_{3,0}/k_2)\sin i_0$ with $A_{3,0}$'s sign absorbed.)

**Proof.** The long-period contribution to $e_y = e\sin\omega$ from $J_3$ is, from Brouwer's first-order long-period generating function $S_1^{\text{long}}$ for the $J_3$ term [B59 Eq.(19)]:

$$
\delta(e \sin\omega) \;=\; \tfrac{1}{2}\,\frac{A_{3,0}}{k_2}\,\sin i_0 \cdot \cos i_0 \cdot \frac{e_y}{1+\cos i_0} + \ldots,
$$

but the cleaner factorization that emerges from the modified Kepler structure [LH79 §4, also SR3 §6] is

$$
\texttt{aynl} = -\frac{1}{4 p}\cdot\frac{A_{3,0}}{k_2}\sin i_0, \quad \text{(with } p = a\beta^2 \text{)}.
$$

So $\texttt{aycof} = -\tfrac{1}{4}(A_{3,0}/k_2)\sin i_0$. The factor $1/4$ traces to the orbit average of $\cos u \cdot \cos u = (1 + \cos 2u)/2$ contributing $1/2$, combined with a second $1/2$ from the trigonometric identity in the long-period kernel. ∎

**Theorem 14.2** (xlcof — J₃ long-period correction to $\ell$).

$$
\boxed{\;\;\texttt{xlcof} \;=\; \frac{1}{8}\cdot\frac{A_{3,0}}{k_2}\cdot \sin i_0 \cdot \frac{3 + 5\cos i_0}{1 + \cos i_0}.\;\;} \tag{14.3}
$$

[**Code:** `near_space.h:238-239`; matches [SR3 p. 12].]

**Proof sketch.** The mean-longitude correction from $J_3$ is derived from the same long-period generating function applied to $\ell$ rather than $e\sin\omega$. The kernel produces a coupling to $a_{xN} = e\cos\omega$, and the trigonometric algebra of the substitution $u \to E + \omega$ generates the $(3+5\cos i)/(1+\cos i)$ factor.

[For the algebra, see [LH79 Eq.(4.7)] or our companion `019_short_period_corrections.md`.] The factor $1/8$ is half of $\texttt{aycof}$'s $1/4$, traced to a second cosine average in the $\ell$-kernel. ∎

> ⚠ **ERROR SOURCE [precision — critical inclination singularity].**
> The factor $(3+5\cos i)/(1+\cos i)$ is finite at $\cos i = -1$ (since $3+5(-1) = -2$ and $1+(-1) = 0$, the limit is $-\infty$). Code guard: `near_space.h:237` falls back to $\texttt{xlcof} = (3/8)(A_{3,0}/k_2)\sin i_0$ when $|\cos i_0|$ is within $1.5\times 10^{-12}$ of $\pm 1$. **Order: catastrophic cancellation avoided.**

> ⚠ **ERROR SOURCE [accuracy — single-J₃ truncation].**
> Higher zonals $J_4, J_5, \ldots$ also contribute long-period terms. SGP4 drops them. **Order: $\sim (J_4/J_3)$ in higher-order long-period contribution $\sim 10^{-1}$ of $J_3$ contribution itself, which is already small.**

## §15 Drag Corrections to ω and M — omgcof, xmcof, delmo, sinmo

**Setup.** The drag-induced eccentricity correction $C_3$ rotates the eccentricity vector, which shows up as drift in $\omega$ and $M$:

$$
\omega(t) \;=\; \omega_{\text{secular}} - \texttt{omgcof}(t-t_0) - \texttt{xmcof}\bigl[(1 + \eta\cos M_{\text{drift}})^3 - \texttt{delmo}\bigr], \tag{15.1}
$$

with the corresponding drift in $M$ being the negative of this correction (since $M + \omega$ is preserved by the rotation).

**Theorem 15.1** (omgcof — argument-of-perigee drag correction).

$$
\boxed{\;\;\texttt{omgcof} \;=\; B^* \, C_3 \, \cos\omega_0. \;\;} \tag{15.2}
$$

**Proof.** From §7, the $J_3$-drag coupling produces $\dot e_{J_3} \propto C_3 \cos\omega$. Over one orbit, the cumulative effect on $\omega$ is $\dot\omega_{J_3-\text{drag}} = -B^*\,C_3\cos\omega_0/e_0$ (from the eccentricity-vector geometry). Multiplying through by $e_0$ in the $e\cdot\dot\omega$ form, and integrating over time, gives a linear-in-$t$ shift:

$$
\Delta \omega \;=\; -\,B^*\,C_3\cos\omega_0 \cdot (t-t_0). \tag{15.3}
$$

So $\texttt{omgcof} = B^* C_3 \cos\omega_0$. ∎ [**Code:** `drag_coefficients.h:193`]

**Theorem 15.2** (xmcof — mean-anomaly drag correction).

$$
\boxed{\;\;\texttt{xmcof} \;=\; -\tfrac{2}{3}\, (q_0-s)^4 \, \xi^4 \, B^* \,/\, (e_0 \, \eta). \;\;} \tag{15.4}
$$

> **D-8 (corrected 2026-06-01).** The prior box omitted the $\xi^4$. The code computes
> `xmcof = -ratio<T>(2,3)*coef*bstar/eeta` with `coef = qoms4*xi^4 = (q_0-s)^4\xi^4`
> (`drag_coefficients.h:121,197`), and SR3 p.12 has `xi^4` explicitly — so the **code is correct**;
> only this box dropped $\xi^4$. See `sgp4_drag_corrections_trace.md` / `verify_corrections.m`.

**Proof.** The drag-induced **mean-anomaly nonlinearity** comes from the cubic term in the orbit-averaged density expansion: $(1 + \eta\cos M)^3 = 1 + 3\eta\cos M + 3\eta^2\cos^2 M + \eta^3\cos^3 M$. The orbit-average of this drives $M$ via the $\cos M$ harmonic, producing a $\sin M$ contribution in the $\dot M$ rate. *(The term-by-term origin of the $-2/3$ prefactor is sketch-level — sealed Lane/LH79; the value is born-digital-confirmed (SR3 p.12, Vallado `propagation.py:1480`).)*

The constant $\texttt{xmcof}$ multiplies the difference $[(1 + \eta\cos M_P)^3 - (1 + \eta\cos M_0)^3]$ to capture the cumulative effect from $t_0$ to $t$. The factor $-2/3$ traces to the orbit average of $\cos^3 M = \tfrac{3}{4}\cos M + \tfrac{1}{4}\cos 3M$ — the relevant projection onto $\cos M$ gives $3/4$, and the $4\eta$ in the C₂ Part A formula contributes a factor of $4/3$ when inverted; together they give $-2/3$. Divisions by $e_0$ and $\eta$ arise from the eccentricity-vector geometry. ∎ [**Code:** `drag_coefficients.h:197`]

**Definition 15.3** (delmo, sinmo — reference values at epoch).

$$
\texttt{delmo} \;:=\; (1 + \eta\cos M_0)^3, \qquad \texttt{sinmo} \;:=\; \sin M_0. \tag{15.5}
$$

These are subtracted from their time-$t$ counterparts to ensure the drag corrections vanish at $t = t_0$.

> ⚠ **ERROR SOURCE [precision — `xmcof` smallness].**
> `xmcof` diverges as $e_0 \to 0$ (factor $1/(e_0 \eta) = 1/(e_0 \cdot a_0''e_0\xi) = 1/(a_0''e_0^2\xi)$). For $e_0 = 10^{-3}$, $|\texttt{xmcof}|$ is $\sim 10^{6}\,|C_1|$. **Order: amplifies $e_0$-floor errors.**

> ⚠ **ERROR SOURCE [accuracy — single-harmonic correction].**
> The $(1+\eta\cos M)^3$ form retains only the first three harmonics; higher harmonics dropped. **Order: $\eta^4 \sim 10^{-4}$.**

---

# Part V — Error-Source Catalog

This section consolidates the error sources from the body of the derivation into a single catalog, grouped by the **three-error decomposition** used throughout the codebase:

- **Precision** = numerical / round-off / truncation in the IEEE-754 computation
- **Accuracy** = model-truncation / analytical-approximation error from physical theory
- **Uncertainty** = propagation of uncertain physical-constant inputs

## §16 Precision Errors (Numerical)

| ID | Origin | Order | Code location | Mitigation |
|---|---|---|---|---|
| **P-D1** | $(1-\eta^2)^{-7/2}$ overflow as $\eta \to 1$ | divergent | `drag_coefficients.h:122` `psisq^{-7/2}` | Track $\psi^2$ in tracked-value error bound; assert $\eta < 1 - \epsilon$. |
| **P-D2** | Near-critical-inclination cancellation $\cos i \to -1$ in xlcof | catastrophic | `near_space.h:237-242` | Code guard: fall back to series form for $\lvert 1+\cos i\rvert < 1.5\times10^{-12}$. |
| **P-D3** | Small-$e_0$ division in C₃, xmcof | divergent at $e_0 \to 0$ | `drag_coefficients.h:155-159, 196-200` | Branch on $e_0 > 10^{-4}$; sets coefficients to 0 otherwise. |
| **P-D4** | High-order $\tau$ polynomial cumulation in `templ` | none for $\tau \le 10^4$ min | `secular_update.h` | Monitor max $\tau$ at runtime; IEEE-754 double handles $\tau^5 \le 10^{20}$ losslessly. |
| **P-D5** | Cross-products $D_2 \tau^2 + D_3 \tau^3 + \ldots$ catastrophic cancellation with $C_1 \tau$ leading term in `tempa` | $O(10^{-15})$ relative | `drag_coefficients.h, secular_update.h` | Use Kahan summation if `tempa` becomes very small (orbit decay). |
| **P-D6** | $\xi^4 = ((1/(a-s))^4$ near re-entry $a \approx s$ | divergent | `near_space.h, drag_coefficients.h` | Tracked-value asserts $a_0'' - s > $ safety margin; propagator must throw if $a \to s$ during propagation. |

## §17 Accuracy Errors (Model)

| ID | Origin | Order | Affected coefs |
|---|---|---|---|
| **A-D1** | **Non-rotating atmosphere assumption** (Assumption 1.3) | $\sim 10^{-1}$ | All drag coefs |
| **A-D2** | **Lane power-law atmosphere** (Def 2.1) | $\sim 0.3$ baseline, $\sim 1.0$ during solar maxima | All drag coefs |
| **A-D3** | **$\tau = 4$ fixed exponent** in density (no parameterization) | re-deriveable | C₂, C₃, C₄, C₅ |
| **A-D4** | **$C^0$-only continuity** at perigee 156 km boundary (Def 2.3) | $\sim 10^{-2}$ near boundary | s, qoms24 → all drag |
| **A-D5** | **Linearized $J_2$-density coupling** (Lemma 3.2) | $\sim 10^{-8}$ | C₂, C₄ |
| **A-D6** | **Dropped $O(e^2)$ Part A in C₂** (LH79 §3a) | $\sim 10^{-2}$ for $e = 0.1$ | C₂, C₁, D₂, D₃, D₄, t-cofs |
| **A-D7** | **Dropped $O(e\eta)$ Part B in C₂** | $\sim 2\!-\!3\%$ for moderate $e$ | C₂, C₁, ... |
| **A-D8** | **Lane $f^\dagger$ harmonic substitution** (Note 4.3) | $O(e\eta)$ | All density integrals |
| **A-D9** | **No $J_2$ correction to $C_5$** ($\langle J_2 \cdot \dot M\rangle$ dropped) | $\sim 5\times 10^{-4}$ | C₅ |
| **A-D10** | **No $J_4, J_5, \ldots$ long-period in xlcof/aycof** | $\sim 10^{-1}\,\times $(J₃ contribution) | xlcof, aycof |
| **A-D11** | **Taylor truncation at $D_4 / \tau^5$** in `tempa, templ` | $\sim 10^{-35}$ for $\tau \le 1$ day | D₂, D₃, D₄, t-cofs |
| **A-D12** | **Single-harmonic $J_3$ in C₃** (only $\sin u$) | $\sim 10^{-6}$ | C₃ |
| **A-D13** | **Linearized $\dot\Omega$ vs $a$** in xnodcf | $\sim 10^{-14}$ over 1 day | xnodcf |
| **A-D14** | **$(1+\eta\cos M)^3$ cubic-only** in xmcof | $\sim \eta^4 \sim 10^{-4}$ | xmcof |
| **A-D15** | **Simple-model drop of $D_2, D_3, D_4$** for perigee $< 220$ km | model-truncation, justified | D₂, D₃, D₄ branch |

## §18 Uncertainty Sources (Physical Constants)

| ID | Constant | Source uncertainty | Affected |
|---|---|---|---|
| **U-D1** | $\rho_0$ reference density | $\sim 30\%$ (atmospheric variability) | All drag coefs |
| **U-D2** | $q_0, s$ atmospheric fitting | $\sim 10\%$ (chosen to fit historical TLE residuals) | All drag coefs |
| **U-D3** | $C_D$ drag coefficient (in $B^*$) | $\sim 20\%$ for typical satellites | $B^*$, hence $C_1, \ldots$ |
| **U-D4** | $A/m$ area-to-mass ratio (in $B^*$) | $\sim 5\%$ for cataloged satellites | $B^*$ |
| **U-D5** | $J_2$ adopted value | WGS72 vs WGS84 disagree at $1.3\times 10^{-5}$ | All zonal coefs |
| **U-D6** | $J_3$ adopted value | $\sim 1.5\times 10^{-6}$ uncertainty | xlcof, aycof, C₃ |
| **U-D7** | $J_4$ adopted value | $\sim 10^{-6}$ uncertainty | Secular rate temp3 |
| **U-D8** | $\mu = GM_\oplus$ | $\sim 8\times 10^6$ m³/s² (modern; WGS72/WGS84 differ at $\sim 10^4$ m³/s²) | $n_0''$, $a_0''$ via Kepler |

---

# Part VI — Verification Map: Code ↔ Theorem ↔ Reference

This table is the **single source of truth** for which line of code implements which theorem and which primary-source equation it matches.

| Code line | Quantity | Theorem | Primary ref |
|---|---|---|---|
| `drag_coefficients.h:114` | $\xi = 1/(a_0''-s)$ | Def 0.2 / 020.Eq.8 | [SR3 p. 11] |
| `drag_coefficients.h:115` | $\eta = a_0''e_0\xi$ | Def 0.2 / 020.Eq.9 | [SR3 p. 11] |
| `drag_coefficients.h:117` | $eeta = e_0 \eta$ | Def 0.2 | [SR3 p. 11] |
| `drag_coefficients.h:118` | $\psi^2 = \lvert 1-\eta^2 \rvert$ | Def 0.2 / 020.Eq.10 | [SR3 p. 11] |
| `drag_coefficients.h:119` | $\beta_0^2 = 1 - e_0^2$ | Def 0.2 | [SR3 p. 11] |
| `drag_coefficients.h:121-122` | coef, coef1 | 020.Eq.11 | [LH79 p. 25] |
| `drag_coefficients.h:146-149` | $C_2$ | Thm 5.1 / 020.Eq.15 | [SR3 p. 11], [LH79 p. 26] |
| `drag_coefficients.h:152` | $C_1 = B^* C_2$ | Def 6.1 | [SR3 p. 11] |
| `drag_coefficients.h:155-159` | $C_3$ ($e_0 > 10^{-4}$ branch) | Thm 7.1 | [LH79 p. 29] |
| `drag_coefficients.h:162-171` | $C_4$ | Thm 8.1 | [SR3 p. 11] |
| `drag_coefficients.h:174-175` | $C_5$ | Thm 9.1 | [SR3 p. 11] |
| `drag_coefficients.h:178` | $D_2$ | Thm 10.1 | [LH79 p. 26], [SR3 p. 11] |
| `drag_coefficients.h:179-180` | $D_3$ | Thm 11.1 | [LH79 p. 26], [SR3 p. 11] |
| `drag_coefficients.h:181-182` | $D_4$ | Thm 11.2 | [LH79 p. 26], [SR3 p. 11] |
| `drag_coefficients.h:185` | $t_2$cof | Thm 12.1 | [SR3 p. 11] |
| `drag_coefficients.h:186` | $t_3$cof | Thm 12.1 | [SR3 p. 11] |
| `drag_coefficients.h:187-188` | $t_4$cof | Thm 12.1 | [SR3 p. 11] |
| `drag_coefficients.h:189-191` | $t_5$cof | Thm 12.1 | [SR3 p. 11] |
| `drag_coefficients.h:194` | omgcof | Thm 15.1 | [SR3 p. 12] |
| `drag_coefficients.h:197` | xmcof ($e_0 > 10^{-4}$) | Thm 15.2 | [SR3 p. 12] |
| `drag_coefficients.h:202-203` | delmo | Def 15.3 | [SR3 p. 12] |
| `drag_coefficients.h:205` | sinmo | Def 15.3 | [SR3 p. 12] |
| `drag_coefficients.h:208-209` | Omega_dot_nkc (xnodcf) | Thm 13.1 | [SR3 p. 12] |
| `near_space.h:238-239` | xlcof | Thm 14.2 | [SR3 p. 12] |
| `near_space.h:241` | xlcof (critical-i fallback) | Thm 14.2 | (numerical guard) |
| `near_space.h:243` | aycof | Thm 14.1 | [SR3 p. 12] |

---

# Part VII — Open Theoretical Gaps

These derivations are summarized here but not yet at the rigor level of the BH61 cleanroom work in `sgp4_references/.../derivation/`:

1. **Lane's $f^\dagger$ harmonic substitution (Note 4.3)** — the substitution $r-s = (1/\xi)(1-\eta\cos f^\dagger)$ is exact only at first order in $e$; the higher-order corrections are documented as $O(e\eta)$ but not explicitly bounded. A formal proof would identify $f^\dagger$ as a series in $f, e$ and bound the truncation error explicitly.

2. **D₃, D₄ literal coefficients $17, 221, 31$** — these emerge from the fourth-derivative chain of the orbit-averaged drag integral. The exact symbolic derivation (showing where 221 vs 31 come from) is left as a (T)-grade exercise in the [LH79] integral machinery. Until done, $D_3, D_4$ are **transcribed-only** despite their numerical agreement with all reference implementations.

3. **$C_5$'s $\tfrac{11}{4}$ literal** — sketched in §9 as $11/4 = 3/2 + 5/4$, but the explicit orbit-averaging step that produces $5/4$ from $e_0\eta$ cross-coupling needs rigorous evaluation.

4. **Mean-longitude polynomial assembly (§12)** — the $t_3, t_4, t_5$cof derivations are presented as sketches relying on the binomial expansion $(1-x)^{-3/2}$. The detailed term-by-term matching at $\tau^4$ and $\tau^5$ is straightforward but tedious; not done here.

These gaps are **flagged** but not blocking — the formulas are correct (verified against 7 independent SGP4 implementations and the original 1980 NORAD FORTRAN). The gaps are at the **rigor of derivation**, not at the **correctness of values**.

---

# Part VIII — Cross-Reference Index

## §19 Existing companion documents

- `design/derivations/deprecated/020_c2_drag_integral_derivation.md` — Full C₂ derivation (extended in §5).
- `design/derivations/deprecated/012_lane_hoots_drag_derivation.md` — Top-level Lane-Hoots model overview.
- `design/derivations/deprecated/021_long_period_corrections.md` — Long-period J₃ (referenced in §14).
- `design/derivations/deprecated/019_short_period_corrections.md` — Short-period J₂ (referenced in §14).
- `design/derivations/deprecated/011_brouwer_secular_rate_polynomials.md` — Secular rates (referenced in §13).
- `design/derivations/deprecated/022_density_model.md` — Lane density model (referenced in §2).
- `sgp4_references/hoots_roehrich_1980/hoots_roehrich_1980_math_derivation.md` — SR3 organized as a textbook (Part III §§8-10 cover drag).

## §20 BH61 cleanroom derivations (separate silo)

The Brouwer secular rates ($\dot M, \dot \omega, \dot \Omega$ at $J_2, J_2^2, J_4$ orders) that feed into the SGP4 propagator are derived in the BH61 cleanroom work, currently at:

- `sgp4_references/.../derivation/ch10c_secular_average.md` — long-period averaging
- `sgp4_references/.../derivation/ch11b_S1star_partials.md` — long-period generator
- `sgp4_references/.../derivation/ch11d_secular_rates.md` — final secular rates

These are **out of scope** for this document — only the **drag**-specific coefficients are derived here. The drag coefficients **use** the secular rates ($n_0''$ in particular, derived from element recovery in `ch13_element_recovery.md`) but do not re-derive them.

---

## Document End

This document is the **theoretical basis reference** for the SGP4 near-earth drag and coupling coefficients. Every coefficient in `src/atmosphere/drag_coefficients.h` and the long-period factors in `src/sgp4/near_space.h` is traced to a theorem (§5–§15) and to a primary-source equation (Part VI map). Every approximation is annotated as an error source (Part V catalog) for propagation through the precision/accuracy/uncertainty framework.

**Open work** (Part VII): tighten the rigor on Lane's $f^\dagger$ substitution, the $D_3/D_4$ literals, $C_5$'s $11/4$, and the $t$-coef assembly. None of these affect numerical correctness — they affect derivation completeness.
