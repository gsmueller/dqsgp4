# Derivation 019: J₂ Short-Period Corrections

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Derive the six short-period corrections applied after the Kepler equation solve:
corrections to radius ($r_k$), argument of latitude ($u_k$), inclination ($i_k$),
RAAN ($\Omega_k$), radial velocity ($\dot{r}_k$), and transverse velocity ($r\dot{f}_k$).

These corrections remove the short-period (one-orbit) oscillations caused by J₂
that were averaged out when computing the secular rates. They are the final step
before the coordinate transformation to Cartesian position/velocity.

**Code:** `src/perturbation/short_period.h` — all 6 corrections on lines 99-114.

**Sources:**
- Brouwer (1959), Eqs. 17-22 — short-period terms from generating function $S_1$
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Brouwer_1959_ADS.pdf`
  (cited as [B59])
- Lyddane (1963), "Small Eccentricities or Inclinations in the Brouwer Theory"
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Lyddane_1963_Small_Eccentricities_Brouwer.pdf`
  (cited as [L63])
- Hoots & Roehrich (1980), Spacetrack Report No. 3, pp. 6, 14 (SGP formulas and SGP4 application)
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Spacetrack_Report_No3_Hoots_Roehrich_1980.pdf`
  (cited as [SR3])

---

## TRANSLATE: Notation

### Symbol Table

| Symbol | Type | Definition | Units | Code Identifier |
|--------|------|-----------|-------|-----------------|
| $S_1$ | function | First-order generating function (from 011.Eq.10) | action×angle | — |
| $r$ | scalar | Osculating radius (from Kepler solve) | ER | `r_osc` |
| $u$ | scalar | Osculating argument of latitude ($= f + \omega$) | rad | `u_osc` |
| $r_k$ | scalar | Corrected radius | ER | `corr.r` |
| $u_k$ | scalar | Corrected argument of latitude | rad | `corr.u` |
| $i_k$ | scalar | Corrected inclination | rad | `corr.i` |
| $\Omega_k$ | scalar | Corrected RAAN | rad | `corr.Omega` |
| $p_L$ | scalar | Semi-latus rectum at current time $= a(1-e_L^2)$ | ER | `pl` |
| $\beta_L$ | scalar | $\sqrt{1-e_L^2}$ at current time | — | `beta_l` |

---

## BUILD Step 1: Origin in the Brouwer Generating Function

The generating function $S_1$ (derived in 011.Eq.10) transforms from mean
(primed) to osculating (unprimed) elements via [B59.Eq.6]:

$$L = L' + \frac{\partial S_1}{\partial l'}, \quad l = l' - \frac{\partial S_1}{\partial L'} \tag{019.Eq.1}$$

$$G = G' + \frac{\partial S_1}{\partial g'}, \quad g = g' - \frac{\partial S_1}{\partial G'} \tag{019.Eq.2}$$

$$H = H', \quad h = h' - \frac{\partial S_1}{\partial H'} \tag{019.Eq.3}$$

**TODO — DERIVE:** Compute each partial derivative of $S_1$ from 011.Eq.10,
then convert from Delaunay variables $(L, G, H, l, g, h)$ to Keplerian elements
$(a, e, i, M, \omega, \Omega)$ using the chain rule relationships from 010.Eq.20-25.

The results should be the short-period oscillations [B59.Eq.17-22]:

For the SGP/SGP4 simplified form, these reduce to the corrections on [SR3] p. 6/14.

---

## BUILD Step 2: The Six Corrections

### 2.1: Radius correction

From [SR3] p. 6 (SGP) and p. 14 (SGP4):

$$r_k = r + \frac{1}{4}J_2\frac{a_E^2}{p_L}\sin^2 i_o \cos 2u \tag{019.Eq.4}$$

**TODO — DERIVE from $S_1$:** The radius correction comes from $\Delta a$ (from
$\partial S_1/\partial l$) converted to $\Delta r$ using the vis-viva relation.
The code (`short_period.h:100-101`) also has a multiplicative correction
$r(1 - \frac{3}{2}(k_2/p^2)\beta(3\cos^2 i-1))$. Show where this comes from.

The coefficient 1/4 must emerge from the Legendre polynomial expansion. ✓/✗

### 2.2: Argument of latitude correction

$$u_k = u - \frac{1}{8}J_2\frac{a_E^2}{p_L^2}(7\cos^2 i_o - 1)\sin 2u \tag{019.Eq.5}$$

**TODO — DERIVE:** The 1/8 coefficient and the $(7\cos^2 i - 1)$ polynomial.
The "7" arises from combining $\partial S_1/\partial G'$ with $\partial S_1/\partial L'$
when converting from $\Delta g$ and $\Delta l$ to $\Delta u = \Delta(f+\omega)$.

### 2.3: RAAN correction

$$\Omega_k = \Omega_{so} + \frac{3}{2}J_2\frac{a_E^2}{p_L^2}\cos i_o \sin 2u \tag{019.Eq.6}$$

**TODO — DERIVE:** From $h - h' = -\partial S_1/\partial H'$. The 3/2 coefficient
and the $\cos i$ factor must emerge from differentiating $S_1$ with respect to $H$.

### 2.4: Inclination correction

$$i_k = i_o + \frac{3}{2}J_2\frac{a_E^2}{p_L^2}\sin i_o \cos i_o \cos 2u \tag{019.Eq.7}$$

**TODO — DERIVE:** From $\Delta i$ via the $G$ and $H$ corrections.
The $\sin i\cos i$ factor comes from $\partial(\cos i)/\partial G$ and $\partial(\cos i)/\partial H$.

### 2.5: Velocity corrections

$$\dot{r}_k = \dot{r} - nJ_2\frac{a_E^2}{2p}\sin^2 i_o \sin 2u \tag{019.Eq.8}$$

$$r\dot{f}_k = r\dot{f} + nJ_2\frac{a_E^2}{2p}\left(\sin^2 i_o \cos 2u + \frac{3}{2}(3\cos^2 i_o - 1)\right) \tag{019.Eq.9}$$

**TODO — DERIVE:** These come from differentiating the position corrections
with respect to time, or equivalently from the velocity components of the
generating function transformation.

---

## PAUSE: Evaluation

### Correctness
- **TODO:** Each coefficient (1/4, 1/8, 3/2, 7) traced to specific $S_1$ partial derivative
- **TODO:** Verify at equatorial ($i=0$): $u_k$ correction $\to -(3/4)(J_2/p^2)\sin 2u$ (from $7\times 1-1=6$, times $1/8 = 6/8 = 3/4$)
- **TODO:** Verify at polar ($i=90°$): inclination correction vanishes ($\sin i\cos i = 0$) ✓

### Accuracy-limiting assumptions
1. First-order in J₂ only — J₂² short-period terms omitted (~7 mm for LEO)
2. J₄ short-period terms omitted
3. Inclination $i_o$ treated as constant (no short-period variation of $i$ fed back)

### Precision improvements
1. The $(7\cos^2 i-1)$ expression: evaluate as $7\cos^2 i - 1 = 6\cos^2 i + (\cos^2 i - 1) = 6\cos^2 i - \sin^2 i$ to avoid cancellation near $\cos^2 i = 1/7$
2. The $\sin 2u$, $\cos 2u$ should be computed via $2\sin u\cos u$ and $\cos^2 u - \sin^2 u$ (already done in code)

---

## Equations to Resolve

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| 019.Eq.1-3 | TEMPLATE | Canonical transformation via $S_1$ partial derivatives |
| 019.Eq.4 | **NOT DERIVED** | Radius correction: 1/4 coefficient from $\partial S_1/\partial l$ → $\Delta a$ → $\Delta r$ |
| 019.Eq.5 | **NOT DERIVED** | Argument of latitude: 1/8 coefficient and $(7\cos^2 i - 1)$ polynomial |
| 019.Eq.6 | **NOT DERIVED** | RAAN: 3/2 coefficient from $-\partial S_1/\partial H$ |
| 019.Eq.7 | **NOT DERIVED** | Inclination: 3/2 coefficient from $\Delta i$ via $G$, $H$ corrections |
| 019.Eq.8-9 | **NOT DERIVED** | Velocity corrections from time-differentiation of position corrections |

## Source Documents Required

| Source | Location | What it provides |
|--------|----------|-----------------|
| [B59] Brouwer (1959) | `Brouwer_1959_ADS.pdf` ✓ | Eqs. 17-22: short-period corrections in Delaunay variables |
| [L63] Lyddane (1963) | `Lyddane_1963_Small_Eccentricities_Brouwer.pdf` ✓ | Reformulation for small eccentricity/inclination |
| [SR3] pp. 6, 14 | `Spacetrack_Report_No3...pdf` ✓ | SGP and SGP4 short-period formulas in Keplerian elements |
| Derivation 011 | `011_brouwer_secular_rate_polynomials.md` ✓ | $S_1$ generating function (011.Eq.10) |

**Assessment:** Fully achievable with existing local sources. The $S_1$ generating function
is already derived in 011. The missing step is computing its 6 partial derivatives
and converting from Delaunay to Keplerian elements.
