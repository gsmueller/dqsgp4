# Derivation 013: Mean Element Recovery (Un-Kozai)

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Derive the iteration that converts TLE mean motion $n_o$ (which includes
secular J₂ effects in the Kozai/Brouwer formulation) to the Brouwer
double-primed mean motion $n_o''$ and semi-major axis $a_o''$ used by all
subsequent SGP4 computations. Every secular rate, drag coefficient, and
short-period correction in SGP4 takes `no_unkozai` ($= n_o''$) as input.

This derivation must show:
1. Why the TLE mean motion differs from the Keplerian mean motion
2. The first-order correction $\delta_1$ and its origin in Brouwer's theory
3. The series reversion that produces the $1 - \delta/3 - \delta^2 - (134/81)\delta^3$ polynomial
4. Why the 134/81 coefficient is exactly that rational number

**Sources:**
- Brouwer (1959) Sections 3-6 and Section 9, pp. 393-395.
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Brouwer_1959_ADS.pdf`
  (cited as [B59])
- Kozai (1959), "The Motion of a Close Earth Satellite," *Astronomical Journal* 64, pp. 367-377.
  (cited as [K59]) — **NOT in local repo, may need to obtain**
- Hoots & Roehrich (1980), Spacetrack Report No. 3, Section 6, p. 10.
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Spacetrack_Report_No3_Hoots_Roehrich_1980.pdf`
  (cited as [SR3])
- Vallado SGP4.cpp, `initl()` function, lines 1242-1251 (the un-Kozai block).

**Code:** `src/orbit/element_recovery.h`

---

## TRANSLATE: Notation

### Symbol Table

| Symbol | Type | Definition | Units | Code Identifier |
|--------|------|-----------|-------|-----------------|
| $n_o$ | scalar | TLE mean motion (as read from element set) | rad/min | `n_input`, Vallado: `no_kozai` |
| $n_o''$ | scalar | Recovered Brouwer mean motion | rad/min | `n0`, Vallado: `no_unkozai` |
| $a_1$ | scalar | Keplerian SMA from $n_o$: $a_1 = (k_e/n_o)^{2/3}$ | ER | — (intermediate) |
| $a_o$ | scalar | First-corrected SMA | ER | — (intermediate) |
| $a_o''$ | scalar | Recovered Brouwer mean SMA | ER | `a0` |
| $\delta_1$ | scalar | J₂ correction evaluated at $a_1$ | — | `del1` |
| $\delta_o$ | scalar | J₂ correction evaluated at $a_o$ | — | `del0` |
| $k_2$ | scalar | $J_2/2$ (= CK2) | — | `half_J2` |
| $k_e$ | scalar | $\sqrt{GM/a_E^3} \times 60$ | rad/min | `ke` |
| $\theta$ | scalar | $\cos i_o$ | — | `cos_i0` |
| $\beta_o$ | scalar | $\sqrt{1-e_o^2}$ | — | `betao` |

**Key distinction:** In Vallado's naming, `no_kozai` is the TLE input and `no_unkozai`
is the recovered Brouwer value. The term "Kozai" here refers to the fact that the
TLE mean motion includes secular J₂ effects (as in Kozai's 1959 formulation).
The "un-Kozai" process removes this to recover the Brouwer mean motion that the
SGP4 analytical theory requires.

---

## BUILD Step 1: The Relationship Between TLE and Brouwer Mean Motions

### 1.1: What the TLE mean motion represents

The TLE mean motion $n_o$ is not a pure Keplerian mean motion. It includes
the secular effect of J₂ on the mean longitude rate. Specifically, Kozai (1959)
showed that the "mean" mean motion (averaged over one orbital period) for a
satellite in an oblate potential is:

$$n_{Kozai} = n_{Keplerian}\left(1 + \delta\right) \tag{013.Eq.1}$$

where $\delta$ is the J₂ secular correction to the mean motion.

**TODO — DERIVE:** Show from [B59.Eq.39] or [K59] how the secular mean motion
correction $\delta$ relates to the Brouwer first-order rate. The key equation is
the Brouwer secular rate for $\dot{l}''$ [B59, p. 393]:

$$\dot{l}'' = n_0\left\{1 + \frac{3}{2}\gamma_2'\eta(-1+3\theta^2) + O(\gamma_2'^2)\right\}$$

At first order, $\delta = \frac{3}{2}\gamma_2'\eta(-1+3\theta^2)$. The relationship
between $\gamma_2'$ evaluated at $a_1$ vs $a_o$ determines the iteration structure.

### 1.2: The δ correction in SGP4 notation

From [SR3] p. 10:

$$\delta_1 = \frac{3}{2}\frac{k_2(3\cos^2 i_o - 1)}{a_1^2(1-e_o^2)^{3/2}} \tag{013.Eq.2}$$

**TODO — DERIVE:** Show that this is exactly $(3/2)k_2(3\theta^2-1)/(a^2\beta^3)$,
which is the first-order Brouwer Ṁ correction from 010.Eq.52 divided by $n$.
Specifically, from 010.Eq.52:

$$\dot{M} = n + \frac{1}{2}T_1\beta(3\theta^2-1) = n\left(1 + \frac{3k_2(3\theta^2-1)}{2a^2\beta^4}\cdot\beta\right) = n(1+\delta)$$

where $\delta = \frac{3k_2(3\theta^2-1)}{2a^2\beta^3}$. This matches (013.Eq.2). ✓

But note: [SR3] p. 10 (SGP model, Section 5) uses $\delta_1 = \frac{3}{4}J_2\frac{a_E^2}{a_1^2}\frac{(3\cos^2 i_o-1)}{(1-e_o^2)^{3/2}}$,
while [SR3] p. 10 (SGP4 model, Section 6) uses $\delta_1 = \frac{3}{2}\frac{k_2(3\cos^2 i_o-1)}{a_1^2(1-e_o^2)^{3/2}}$.
Since $k_2 = J_2 a_E^2/2$ (in Brouwer's convention with dimensions) or $k_2 = J_2/2$ (dimensionless, with $a_E$ factored into $a$), these are equivalent when $a$ is in Earth radii.

**TODO — VERIFY:** Confirm the dimensional analysis. In SGP4, $a$ is in Earth radii
(so $a_E = 1$), and $k_2 = J_2/2$ (dimensionless). Then $(3/2)(J_2/2)(3\theta^2-1)/(a^2\beta^3) = (3/4)J_2(3\theta^2-1)/(a^2\beta^3)$. The [SR3] SGP form has $(3/4)J_2 a_E^2/a^2$, so with $a_E = 1$: identical. ✓

---

## BUILD Step 2: The Series Reversion — Deriving 1/3, 1, 134/81

### 2.1: The inversion problem

We know $n_o = n_{Brouwer}(1 + \delta(a))$ where $a$ depends on $n_{Brouwer}$ via
Kepler's third law: $a = (k_e/n)^{2/3}$.

Starting from the TLE $n_o$, we compute $a_1 = (k_e/n_o)^{2/3}$ and
$\delta_1 = \delta(a_1)$. But $a_1$ is not the true Brouwer SMA — it's
contaminated by the $\delta$ correction. We need to invert: find $a_o$
such that when we compute $n_o = (k_e/a_o^{3/2})(1+\delta(a_o))$, we
recover the observed $n_o$.

The correction formula from [SR3] p. 10:

$$a_o = a_1\left(1 - \frac{1}{3}\delta_1 - \delta_1^2 - \frac{134}{81}\delta_1^3\right) \tag{013.Eq.3}$$

**TODO — DERIVE:** This is a series reversion. The argument proceeds as follows:

From Kepler's third law, $n = k_e/a^{3/2}$, so $a = (k_e/n)^{2/3}$.

If $n_{true} = n_o/(1+\delta)$ and $a_{true} = (k_e/n_{true})^{2/3}$, then:

$$a_{true} = \left(\frac{k_e}{n_o/(1+\delta)}\right)^{2/3} = \left(\frac{k_e(1+\delta)}{n_o}\right)^{2/3} = a_1(1+\delta)^{2/3} \tag{013.Eq.4}$$

But $\delta$ itself depends on $a_{true}$: $\delta = c/a_{true}^2$ where $c = (3/2)k_2(3\theta^2-1)/\beta^3$.

So $\delta_1 = c/a_1^2$ and $\delta_{true} = c/a_{true}^2 = \delta_1(a_1/a_{true})^2$.

Substituting (013.Eq.4): $a_{true} = a_1(1+\delta_{true})^{2/3}$, so
$\delta_{true} = \delta_1/(1+\delta_{true})^{4/3}$.

This is a fixed-point equation. The series reversion expands $a_{true}/a_1$
in powers of $\delta_1$:

$$\frac{a_{true}}{a_1} = 1 + c_1\delta_1 + c_2\delta_1^2 + c_3\delta_1^3 + \ldots \tag{013.Eq.5}$$

**TODO — DERIVE $c_1, c_2, c_3$:**

From (013.Eq.4): $a_{true}/a_1 = (1+\delta_{true})^{2/3}$.

And $\delta_{true} = \delta_1(a_1/a_{true})^2$.

Let $x = a_{true}/a_1$ and $\delta_{true} = \delta_1/x^2$. Then:

$$x = (1 + \delta_1/x^2)^{2/3} \tag{013.Eq.6}$$

Expand $(1+u)^{2/3} = 1 + \frac{2}{3}u - \frac{1}{9}u^2 + \frac{4}{81}u^3 - \ldots$ where $u = \delta_1/x^2$:

$$x = 1 + \frac{2}{3}\frac{\delta_1}{x^2} - \frac{1}{9}\frac{\delta_1^2}{x^4} + \frac{4}{81}\frac{\delta_1^3}{x^6} - \ldots \tag{013.Eq.7}$$

Now substitute $x = 1 + c_1\delta_1 + c_2\delta_1^2 + c_3\delta_1^3 + \ldots$ and
$1/x^2 = 1 - 2c_1\delta_1 + (3c_1^2 - 2c_2)\delta_1^2 + \ldots$:

**At $O(\delta_1)$:** $c_1 = 2/3 \cdot 1 = 2/3$

Hmm — but [SR3] has $-1/3$, not $+2/3$.

**TODO — RESOLVE SIGN:** The discrepancy may be because $\delta_1$ in [SR3]
has the opposite sign from $\delta$ in our formulation, or because the
expansion should be of $a_o/a_1$ not $a_{true}/a_1$. The [SR3] formula uses
$a_o = a_1(1-\delta_1/3-\delta_1^2-134\delta_1^3/81)$, and then computes
$\delta_0$ at $a_o$ and sets $n_o'' = n_o/(1+\delta_0)$. So $a_o$ is NOT
$a_{true}$ — it's a one-step approximation that gets refined by recomputing
$\delta$ at the new $a$.

The key issue is: **what quantity is being expanded?** Two possible readings:

**(A)** Direct inversion: $n_{true} = n_o/(1+\delta)$ where $\delta = f(a_{true})$
and $a_{true} = (k_e/n_{true})^{2/3}$. This gives a self-consistent equation.

**(B)** Perturbative correction: $a_o = a_1(1 + \text{correction})$ where the
correction removes the $\delta_1$ bias. Since $a_1 = (k_e/n_o)^{2/3}$ includes
the $\delta$ effect, we need $a_o < a_1$ (the true SMA is smaller because the
mean motion was enhanced by J₂). So the correction should be negative, consistent
with $-\delta_1/3$.

**TODO — COMPLETE THE SERIES REVERSION:**

The approach: From $n_o = n_{Brouwer}(1+\delta)$ and $a_1 = (k_e/n_o)^{2/3}$,
$a_{Brouwer} = (k_e/n_{Brouwer})^{2/3}$:

$$a_1 = \left(\frac{k_e}{n_{Brouwer}(1+\delta)}\right)^{2/3} = a_{Brouwer}(1+\delta)^{-2/3} \tag{013.Eq.8}$$

So: $a_{Brouwer} = a_1(1+\delta)^{2/3}$

$(1+\delta)^{2/3} = 1 + \frac{2}{3}\delta - \frac{1}{9}\delta^2 + \frac{4}{81}\delta^3 - \ldots$

But $\delta = \delta(a_{Brouwer}) \neq \delta_1 = \delta(a_1)$.

Since $\delta \propto 1/a^2$: $\delta = \delta_1(a_1/a_{Brouwer})^2 = \delta_1/(1+\delta)^{4/3}$.

At zeroth order: $\delta \approx \delta_1$.
At first order: $\delta \approx \delta_1(1 - \frac{4}{3}\delta_1)$.
At second order: $\delta \approx \delta_1(1 - \frac{4}{3}\delta_1 + \frac{22}{9}\delta_1^2)$.

Substituting back:
$$a_{Brouwer}/a_1 = 1 + \frac{2}{3}\delta_1 - \frac{2}{3}\cdot\frac{4}{3}\delta_1^2 + \ldots$$

This gives $+2/3$ at first order, not $-1/3$. The sign depends on whether
$\delta_1$ is defined positive or negative. From (013.Eq.2):
$\delta_1 = (3/2)k_2(3\cos^2 i-1)/(a_1^2\beta^3)$.

For $i < 54.7°$ (most orbits): $3\cos^2 i - 1 > 0$, so $\delta_1 > 0$.
The Brouwer mean motion $n''$ is LARGER than Keplerian (equatorial orbits
move faster due to oblateness), so $n_o = n'' > n_{Keplerian}$ and
$a_1 = (k_e/n_o)^{2/3} < a_{Keplerian}$. The Brouwer $a''$ should be LARGER
than $a_1$ since $n'' = k_e/a''^{3/2}$ and $n'' < n_o$ (wait, $n'' = n_o/(1+\delta)$
which is $< n_o$ if $\delta > 0$, so $a'' > a_1$). ✓

So $a'' = a_1(1 + 2\delta_1/3 + \ldots)$, but [SR3] has $a_o = a_1(1 - \delta_1/3 + \ldots)$.

**TODO — RESOLVE:** This means [SR3]'s $a_o$ is NOT $a_{Brouwer}$ directly.
Read the SGP4 code more carefully — the iteration structure suggests $a_o$ is
an intermediate that gets refined by recomputing $\delta_0$ at $a_o$, and then
$n_o'' = n_o/(1+\delta_0)$ and $a_o'' = (k_e/n_o'')^{2/3}$. The series reversion
may be for a DIFFERENT quantity than $a_{Brouwer}/a_1$.

Possible resolution: The [SR3] formula may be inverting $n_o = n_{Keplerian}(1+\delta)$
to find $n_{Keplerian}$ (not $n_{Brouwer}$), where $n_{Keplerian} = k_e/a^{3/2}$.
Then $a_o$ is the Keplerian SMA corresponding to the UNPERTURBED mean motion.

---

## BUILD Step 3: The 134/81 Coefficient

### 3.1: Binomial expansion of $(1+\delta)^{-2/3}$

**TODO — DERIVE:** The coefficient 134/81 arises from the third-order term in
the series reversion. The binomial expansion:

$$(1+\delta)^{-2/3} = 1 - \frac{2}{3}\delta + \frac{2\cdot 5}{3\cdot 6}\delta^2 - \frac{2\cdot 5\cdot 8}{3\cdot 6\cdot 9}\delta^3 + \ldots$$

$$= 1 - \frac{2}{3}\delta + \frac{5}{9}\delta^2 - \frac{40}{81}\delta^3 + \ldots$$

But this gives 40/81, not 134/81. The discrepancy is because the series
reversion involves the SELF-CONSISTENT equation $\delta = f(a(\delta))$, not a
simple binomial expansion. The 134/81 includes cross-terms from substituting
$\delta = \delta_1(a_1/a_o)^2$ back into the expansion.

**TODO — COMPUTE:** Carry out the full series reversion to third order:

Let $a_o/a_1 = 1 + \alpha_1\delta_1 + \alpha_2\delta_1^2 + \alpha_3\delta_1^3$.

From $a_o = a_1(1+\delta(a_o))^{-2/3}$ (or the appropriate relation):
- Express $\delta(a_o) = \delta_1(a_1/a_o)^2 = \delta_1/(1+\alpha_1\delta_1+\ldots)^2$
- Expand $(1+\delta(a_o))^{-2/3}$ in powers of $\delta_1$
- Match coefficients order by order
- Verify $\alpha_1 = -1/3$, $\alpha_2 = -1$, $\alpha_3 = -134/81$

### 3.2: Cross-check

**TODO — VERIFY:** The integers in 134/81:
- $134 = 2 + 132 = 2 + 4\times33$ ... unclear pattern
- $81 = 3^4$
- The coefficient should emerge from combining binomial coefficients with
  the self-consistency iteration. Specifically, the generalized binomial
  coefficient $\binom{-2/3}{3}$ combined with the feedback from $\delta \propto a^{-2}$.

**TODO — VERIFY NUMERICALLY:** For a test orbit with known $\delta_1$, verify
that the cubic formula gives the same $a_o$ as direct numerical iteration.

---

## BUILD Step 4: The Two-Step Iteration Structure

### 4.1: Why two steps suffice

[SR3] p. 10 computes $a_o$ from the cubic, then recomputes $\delta_0$ at $a_o$,
then sets $n_o'' = n_o/(1+\delta_0)$.

**TODO — DERIVE:** Show that this two-step process converges to $O(\delta^4)$
accuracy. The cubic gives $a_o$ accurate to $O(\delta_1^4)$. Recomputing
$\delta_0$ at this improved $a_o$ absorbs the $O(\delta_1^4)$ error into a
new correction that is $O(\delta_1^4)$ itself, so the final $n_o''$ is accurate
to $O(\delta_1^4)$.

For LEO with $\delta_1 \sim 10^{-3}$: the residual is $O(10^{-12})$, well below
any physical measurement capability. ✓

### 4.2: The final step — $a_o''$ from $n_o''$

The code computes $a_o'' = (k_e/n_o'')^{2/3}$ via Kepler's third law (NOT
from $a_o/(1-\delta_0)$, which would differ at $O(\delta^2)$). This is noted
in `element_recovery.h` line 123.

**TODO — DERIVE:** Show the $O(\delta^2)$ difference between
$a_o'' = (k_e/n_o'')^{2/3}$ and the approximation $a_o/(1-\delta_0)$. The code
comment states this was a source of systematic position errors at $t > 0$.

---

## PAUSE: Evaluation

### Correctness
- **TODO:** Verify $-1/3$ coefficient at first order
- **TODO:** Verify $-1$ coefficient at second order
- **TODO:** Verify $-134/81$ coefficient at third order
- **TODO:** Numerical cross-check against direct iteration for 3+ test orbits

### Accuracy-limiting assumptions
1. **Truncation at $O(\delta^3)$:** The next term would be $O(\delta^4) \sim O(J_2^4/a^8)$.
   For LEO: $\delta \sim 10^{-3}$, so $\delta^4 \sim 10^{-12}$ — sub-picometer effect.
2. **$e$ and $i$ treated as unaffected by the Kozai-Brouwer transformation:**
   In reality, the canonical transformation also shifts $e$ and $i$ at $O(J_2)$,
   but this effect is absorbed into the short-period corrections applied later.

### Precision improvements
1. **The cubic polynomial evaluation:** Use Horner form:
   $1 - \delta_1(1/3 + \delta_1(1 + (134/81)\delta_1))$ — 3 multiplies, 3 adds.
   Already in code. ✓
2. **Cube root via Newton iteration:** The code uses Newton's method for
   $(k_e/n)^{2/3}$ instead of `pow()`. This preserves error tracking through
   the TrackedValue framework. ✓

### Generalizations
1. **Higher-order reversion:** A fourth-order term could be derived for
   extremely low orbits where $\delta \sim 10^{-2}$, but is not needed for SGP4.
2. **Full canonical transformation:** Instead of the perturbative inversion,
   one could apply the full Brouwer canonical transformation (mean → osculating)
   and its inverse. This is the "enhanced preset" approach.

---

## Equations to Resolve

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| 013.Eq.1 | TEMPLATE | $n_{Kozai} = n_{Brouwer}(1+\delta)$ — relationship to [B59.Eq.39] Ṁ secular rate |
| 013.Eq.2 | LINKED to 010.Eq.52 | $\delta_1 = (3/2)k_2(3\theta^2-1)/(a^2\beta^3)$ — already derived as Ṁ/n-1 |
| 013.Eq.3 | **CRITICAL — NOT DERIVED** | $a_o = a_1(1-\delta_1/3-\delta_1^2-134\delta_1^3/81)$ — series reversion |
| 013.Eq.4-8 | TEMPLATE | Setup for the series reversion derivation |
| 134/81 | **CRITICAL — NOT DERIVED** | Must emerge from the third-order self-consistent expansion |
| Two-step convergence | TEMPLATE | Show $O(\delta^4)$ residual after two-step process |
| $a_o''$ vs $a_o/(1-\delta_0)$ | TEMPLATE | Show $O(\delta^2)$ difference and why Kepler's law is preferred |

## Source Documents Required

| Source | Location | What it provides |
|--------|----------|-----------------|
| [B59] Brouwer (1959) | `Brouwer_1959_ADS.pdf` ✓ | Secular rate $\dot{l}''$ (Eq. 39), generating function structure |
| [SR3] Spacetrack Report 3 | `Spacetrack_Report_No3...pdf` ✓ | The element recovery formulas (p. 10), Section 5 (SGP/Kozai) and Section 6 (SGP4/Brouwer) |
| [K59] Kozai (1959) | **NOT in local repo** | Original Kozai mean motion formulation. May not be strictly needed — the relationship is implicit in [B59] and [SR3] |
| [LH79] Lane & Hoots (1979) | `Lane_Hoots_1979...pdf` ✓ | Section 5 (geopotential simplification) shows how SGP4 simplified the element recovery from AFGP4 |

**Assessment:** The derivation is fully achievable with existing local sources ([B59] and [SR3]).
The series reversion is self-contained algebra once the starting relationship $n_o = n_{Brouwer}(1+\delta(a))$
is established. [K59] would add historical context but is not mathematically required.
