# Derivation 021: J₃ Long-Period Corrections (Lyddane Modification)

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Derive the long-period corrections to the eccentricity vector and mean longitude
that arise from the J₃ (odd) zonal harmonic. These are the `xlcof`, `aycof`
coefficients and their application via $a_{yNSL}$ and $\delta L$ in the SGP4
propagation step.

In the Brouwer theory, J₃ produces NO secular terms (because $P_3$ is odd in
$\sin u$, its orbit average vanishes), but it DOES produce long-period terms
that vary with $\omega$ (the argument of perigee). These long-period terms
are incorporated into the modified Kepler equation by adjusting the eccentricity
vector component $a_{yN}$ and the mean longitude $L$.

The Lyddane (1963) modification reformulates these corrections to avoid the
singularity at $e = 0$ that appears in the standard Brouwer formulation.

**Code:**
- `src/sgp4/near_space.h` lines 233-243 — `xlcof`, `aycof` initialization
- `src/sgp4/near_space.h` lines 300-303 — `xll`, `aynl` application during propagation
- `src/sgp4/deep_space.h` lines 236-242, 292-294 — same for deep-space

**Sources:**
- Brouwer (1959), Section 8 (third and fifth harmonics), pp. 390-391
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Brouwer_1959_ADS.pdf`
  (cited as [B59])
- Lyddane (1963), "Small Eccentricities or Inclinations in the Brouwer Theory"
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Lyddane_1963_Small_Eccentricities_Brouwer.pdf`
  (cited as [L63])
- Hoots & Roehrich (1980), Spacetrack Report No. 3, pp. 12-13
  PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Spacetrack_Report_No3_Hoots_Roehrich_1980.pdf`
  (cited as [SR3])

---

## TRANSLATE: Notation

### Symbol Table

| Symbol | Type | Definition | Units | Code Identifier |
|--------|------|-----------|-------|-----------------|
| $A_{3,0}$ | scalar | $-J_3 a_E^3/2$ — Brouwer's J₃ coefficient | length³ | `A_30`, `A30` |
| $k_2$ | scalar | $J_2/2$ | — | `CK2`, `half_J2` |
| $A_{3,0}/k_2$ | scalar | Ratio used in long-period formulas = $-J_3 a_E/(J_2)$ | length | Vallado: `j3oj2` |
| $a_{xN}$ | scalar | $e\cos\omega$ — x-component of eccentricity vector | — | `axN` |
| $a_{yN}$ | scalar | $e\sin\omega$ — y-component of eccentricity vector | — | — |
| $a_{yNSL}$ | scalar | $a_{yN}$ + long-period correction from J₃ | — | `ayn` in code |
| $\delta L$ | scalar | Long-period correction to mean longitude from J₃ | rad | `xll` |
| `xlcof` | scalar | Precomputed coefficient for $\delta L$ | — | `xlcof` |
| `aycof` | scalar | Precomputed coefficient for $a_{yNSL}$ correction | — | `aycof` |

---

## BUILD Step 1: J₃ Disturbing Function and Long-Period Terms

### 1.1: The J₃ potential

The J₃ term in the geopotential is [B59, p. 390]:

$$U_3 = \frac{\mu A_{3,0}}{r^4}\left(\frac{5}{2}\sin^3\beta - \frac{3}{2}\sin\beta\right) \tag{021.Eq.1}$$

Using $\sin\beta = \sin i\sin(g+f)$ and the identity $\sin^3 u = \frac{3}{4}\sin u - \frac{1}{4}\sin 3u$:

$$U_3 = \frac{\mu A_{3,0}}{r^4}\sin i\left[\left(-\frac{3}{2}+\frac{15}{8}\sin^2 i\right)\sin(g+f) - \frac{5}{8}\sin^3 i\sin(3g+3f)\right] \tag{021.Eq.2}$$

**TODO — DERIVE:** Show the full expansion from $P_3(\sin\beta)$ substituted with
$\sin\beta = \sin i\sin u$, then separate into harmonics of $u = g + f$.

### 1.2: Why J₃ produces no secular terms

The orbit average of $U_3$ vanishes because every term contains $\sin(ng + mf)$
with $m \neq 0$, and $\int_0^{2\pi}\sin(mf)\,dM = 0$ for integer $m \geq 1$.

**TODO — VERIFY:** Show explicitly that $\langle U_3 \rangle = 0$.

### 1.3: The long-period generating function $\Delta_3 S_1^*$

Brouwer [B59, p. 391] gives the J₃ contribution to the long-period
generating function. The secular part of $F$ gains no J₃ term, but
the long-period generating function $S_1^*$ (which removes terms
depending on $g'$ but not $l'$) gains the J₃ correction:

$$\frac{\partial}{\partial g'}\Delta_3 S_1^* = -\frac{1}{4}\frac{\mu A_{3,0}}{L'^2 k_2}\frac{L'^2}{G''^2}G''e''\sin I''\sin g' \tag{021.Eq.3}$$

Integrating:

$$\Delta_3 S_1^* = +\frac{1}{4}\frac{\mu A_{3,0}}{L'^2 k_2}\frac{L'^2}{G''^2}G''e''\sin I''\cos g' \tag{021.Eq.4}$$

**TODO — DERIVE:** Trace this from the homological equation for the long-period
transformation: $\frac{\partial F_{1s}^*}{\partial G'}\frac{\partial\Delta_3 S_1^*}{\partial g'} + \Delta_3 F_{2p}^* = 0$.

### 1.4: Effect on orbital elements

From the canonical transformation [B59, pp. 390-391]:

$$\Delta_3 G' = -\frac{\partial\Delta_3 S_1^*}{\partial g''} = \frac{1}{4}\frac{A_{3,0}}{k_2}\frac{L'^2}{G''^2}G''e''\sin I''\sin g'' \tag{021.Eq.5}$$

$$\Delta_3 g' = +\frac{\partial\Delta_3 S_1^*}{\partial G''} \tag{021.Eq.6}$$

$$\Delta_3 h' = +\frac{\partial\Delta_3 S_1^*}{\partial H''} \tag{021.Eq.7}$$

**TODO — DERIVE:** Compute all three partial derivatives explicitly, convert
to Keplerian elements using the chain rule.

---

## BUILD Step 2: The SGP4 Long-Period Formulas

### 2.1: The $a_{yNSL}$ correction [SR3, p. 12]

$$a_{yNSL} = e\sin\omega_{so} - \frac{1}{2}\frac{J_3}{J_2}\frac{a_E}{p}\sin i_o \tag{021.Eq.8}$$

The second term is the long-period J₃ correction. In code notation:

`aycof` = $\frac{1}{4}\frac{A_{3,0}}{k_2}\sin i_o$

since $\frac{1}{2}\frac{J_3}{J_2}\frac{a_E}{p} = \frac{1}{2}\frac{-2A_{3,0}/a_E^2}{2k_2/a_E^2}\frac{a_E}{p} = -\frac{A_{3,0}}{2k_2 p}$
and the correction becomes $+\frac{A_{3,0}}{2k_2 p}\sin i = \frac{1}{4}\frac{A_{3,0}}{k_2}\sin i \times \frac{2}{p}$...

**TODO — DERIVE:** Trace the exact path from (021.Eq.5) to the `aycof` formula.
Show the dimensional analysis carefully (Brouwer's $k_2$ has dimensions length²
while the code's `half_J2` is dimensionless with $a_E$ factored out).

### 2.2: The mean longitude correction [SR3, p. 12]

$$L = L_s - \frac{1}{4}\frac{J_3}{J_2}\frac{a_E}{p}a_{xNSL}\sin i_o\left[\frac{3+5\cos i_o}{1+\cos i_o}\right] \tag{021.Eq.9}$$

The factor $(3+5\cos i)/(1+\cos i)$ is the key coefficient. In code:

`xlcof` = $\frac{1}{8}\frac{A_{3,0}}{k_2}\sin i_o\frac{(3+5\cos i_o)}{(1+\cos i_o)}$

**TODO — DERIVE:** Show where $(3+5\cos i)/(1+\cos i)$ comes from. This factor
arises from combining the $\Delta_3 l'$ and $\Delta_3 g'$ corrections when computing
$\Delta(l + g) = \Delta M + \Delta\omega$. The $(1+\cos i)$ in the denominator comes
from the Lyddane reformulation that avoids the $e$ divisor.

**TODO — DERIVE the Lyddane singularity avoidance:** The standard Brouwer
J₃ correction to $\omega$ contains $1/e$, which is singular for circular orbits.
Lyddane (1963) showed that the correction to $e\sin\omega$ (= $a_{yN}$) is
well-behaved, and the correction to $l + g$ (mean longitude) avoids the $1/e$
singularity. Show this explicitly.

### 2.3: The near-retrograde singularity

The $(1+\cos i)$ denominator is singular at $i = 180°$ (retrograde). The code
handles this with a branch (`near_space.h:237-241`):

```
if (cosio > -1 + 1.5e-12 && cosio < 1 - 1.5e-12) {
    xlcof = (1/8) * A30/CK2 * sinio * (3 + 5*cosio) / (1 + cosio);
} else {
    xlcof = (3/8) * A30/CK2 * sinio;  // limit as cosio → ±1
}
```

**TODO — DERIVE:** Show that $\lim_{i\to\pi}(3+5\cos i)/(1+\cos i) = (3-5)/(1-1)$
is indeterminate, and evaluate via L'Hôpital or Taylor expansion. The code uses
the value $3/8 \times (A_{3,0}/k_2)\sin i$ as the limiting form. Verify this.

---

## PAUSE: Evaluation

### Correctness
- **TODO:** Trace each coefficient (1/4, 1/8, 3, 5) to specific Brouwer equation
- **TODO:** Verify $(3+5\cos i)/(1+\cos i)$ simplifies correctly at $i = 0$ (equatorial): $(3+5)/2 = 4$
- **TODO:** Verify the near-retrograde limit

### Accuracy-limiting assumptions
1. J₃ only — no J₅ long-period terms (J₅ is ~1000× smaller than J₃)
2. The J₃ correction is applied as a constant offset to $a_{yN}$ — no time variation
   within the propagation (the variation comes through $\omega$ changing secularly)

### Precision improvements
1. The $(3+5\cos i)/(1+\cos i)$ factor: for $i$ near $180°$, use the identity
   $(3+5c)/(1+c) = 5 - 2/(1+c)$ which shows the limit is $5 - 2/0^+ = -\infty$.
   The code's branch at $\cos i \approx -1$ is necessary.
2. For the enhanced preset: include J₅ long-period terms [B59, pp. 391-393]

---

## Equations to Resolve

| Equation | Status | What must be shown |
|----------|--------|--------------------|
| 021.Eq.1-2 | TEMPLATE | J₃ potential expansion in orbital elements |
| 021.Eq.3-4 | **NOT DERIVED** | Long-period generating function from [B59] p. 391 |
| 021.Eq.5-7 | **NOT DERIVED** | Canonical transformation: $\Delta G'$, $\Delta g'$, $\Delta h'$ from $S_1^*$ |
| 021.Eq.8 | **NOT DERIVED** | $a_{yNSL}$ formula — `aycof` coefficient origin |
| 021.Eq.9 | **NOT DERIVED** | Mean longitude correction — `xlcof` with $(3+5\cos i)/(1+\cos i)$ |
| Lyddane singularity | **NOT DERIVED** | Why $e\sin\omega$ correction is well-behaved while $\omega$ correction has $1/e$ |
| Near-retrograde limit | **NOT DERIVED** | L'Hôpital evaluation of $(3+5c)/(1+c)$ at $c = -1$ |

## Source Documents Required

| Source | Location | What it provides |
|--------|----------|-----------------|
| [B59] Section 8 | `Brouwer_1959_ADS.pdf` pp. 390-391 ✓ | J₃ disturbing function, $\Delta_3 S_1^*$, canonical corrections |
| [B59] Section 9 | `Brouwer_1959_ADS.pdf` pp. 394-395 ✓ | Computational formulas for J₃ long-period terms |
| [L63] Lyddane (1963) | `Lyddane_1963...pdf` ✓ | Singularity-free reformulation for small $e$ |
| [SR3] pp. 12-13 | `Spacetrack_Report_No3...pdf` ✓ | SGP4 implementation of long-period corrections |

**Assessment:** Fully achievable with existing local sources. All three required
papers ([B59], [L63], [SR3]) are available in the repository. The derivation
traces from the J₃ potential → generating function → canonical transformation
→ Keplerian corrections → Lyddane reformulation → code coefficients.
