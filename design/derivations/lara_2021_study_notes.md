# Study Notes: Lara (2021) "Brouwer's satellite solution redux"

arXiv:2009.10665v2, November 3, 2020
Source PDF: `sgp4_references/vallado_celestrak/documentation/SGP4/Lara_2021_Brouwer_Redux_arXiv_2009.10665v2.pdf`

## Purpose
Unlock the derivation chain behind Brouwer (1959) and Brouwer & Hori (1961).
Specifically: understand how S₁ is constructed, how δp_j/δq_j arise, and how
the Fourier-eccentricity expansion machinery works.

## Notation Map: Lara ↔ Brouwer & Hori (1961)

| Lara | Brouwer & Hori | Meaning |
|------|----------------|---------|
| $\mathcal{W}_1$ | $S_1$ | First-order determining/generating function |
| $\mathcal{C}_1$ | $S_1^*$ | Integration constant (long-period part) |
| $\epsilon = -C_{2,0} = J_2$ | $k_2$ (note: $k_2 = -\frac{1}{2}C_{2,0}R_\oplus^2$) | Small parameter |
| $s = \sin I$ | $\sin I$ | Sine of inclination |
| $\theta = \cos I$ | $\theta = H/G$ | Cosine of inclination |
| $B_0 = 1-\frac{3}{2}s^2$ | $(-\frac{1}{2}+\frac{3}{2}\theta^2)$ | Inclination polynomial (same thing) |
| $B_1 = \frac{3}{4}s^2$ | $(\frac{3}{2}-\frac{3}{2}\theta^2) \cdot \frac{1}{2}$ | Inclination polynomial |
| $\phi = f - \ell$ | equation of the center | True minus mean anomaly |
| $p = a\eta^2$ | $a(1-e^2)$ | Semi-latus rectum |
| $\Delta\xi = J_2\{\xi, \mathcal{W}_1\}$ | $\delta p_j, \delta q_j$ via operator $D$ | First-order periodic corrections |

## Derivation Chain

### Step 1: The Hamiltonian (Lara Eq 1, p.3)

$$\mathcal{H} = -\frac{\mu}{2a} + \frac{\mu R_\oplus^2}{r \cdot r^2} \frac{1}{2} C_{2,0}\left[1 - \frac{3}{2}s^2 + \frac{3}{2}s^2\cos(2f+2\omega)\right]$$

This is the gravitational potential of an oblate Earth in Delaunay variables.
- First term: Kepler (two-body)
- Second term: J₂ perturbation, split into:
  - Secular part: $(1-\frac{3}{2}s^2)$ — does not depend on $f$ or $g$
  - Short-period part: $\frac{3}{2}s^2\cos(2f+2g)$ — depends on true anomaly $f$

**Connection to BH61**: This is the force function $U$ in BH61 Eq (3), written in
the specific form for J₂ only.

### Step 2: Eccentricity function expansion (Lara Eq 3-4, p.5)

The key identity that converts powers of $1/r$ into Fourier series:

$$\frac{1}{r^j} = \frac{1}{r^2}\left(\frac{1+e\cos f}{p}\right)^{j-2}$$

This means every $r^{-j}$ term can be expanded as a polynomial in $e\cos f$,
which then expands into $\cos(kf)$ harmonics. Combined with $\cos(2f+2g)$,
this produces terms $\cos(jf+2ig)$.

The Hamiltonian perturbation becomes (Lara Eq 4):

$$\tilde{\mathcal{H}}_{0,1} = \mathcal{H}_{1,0} = \frac{R_\oplus^2}{r^2}\frac{1}{\eta^2}\sum_{i=0}^{1} B_i(s) \sum_{j=i}^{2i+1}(2-j^*)^i e^{|j-2i|}\cos(jf+2ig)$$

where $j^* = j \mod 2$.

**Connection to BH61**: This is what BH61 calls $\partial S_1/\partial l_1''$ after differentiation —
the perturbation expressed as a Fourier series in $f$ and $g$ with eccentricity
coefficients. The $B_i(s)$ are inclination functions.

### Step 3: The homological equation — constructing S₁ (Lara Eq 2, p.5)

$$\mathcal{W}_m = \frac{1}{n}\int(\tilde{\mathcal{H}}_{0,m} - \mathcal{H}_{0,m})d\ell + \mathcal{C}_m$$

**This is the central equation.** The generating function $\mathcal{W}_1$ (= Brouwer's $S_1$)
is obtained by:
1. Take the perturbation Hamiltonian $\tilde{\mathcal{H}}_{0,1}$
2. Subtract its average over $\ell$ (the new Hamiltonian $\mathcal{H}_{0,1}$)
3. Integrate over $\ell$
4. Add an integration constant $\mathcal{C}_1$ (chosen to eliminate long-period terms)

The average $\mathcal{H}_{0,1} = \langle\tilde{\mathcal{H}}_{0,1}\rangle_\ell$ removes all terms that depend on $\ell$
(i.e., all short-period terms), leaving:

$$\mathcal{H}_{0,1} = \mathcal{H}_{0,0}\frac{R_\oplus^2}{p^2}\eta\left(1-\frac{3}{2}s^2\right)$$

**Connection to BH61**: $\mathcal{H}_{0,1}$ is Brouwer's $F^{**}$ (the averaged Hamiltonian).
The secular frequencies $\partial F^{**}/\partial L_j''$ in BH61 Eq (6) come from this.

### Step 4: Integration via equation of the center (Lara Eq 5-6, p.6)

The integration $\int ... d\ell$ is converted to $\int ... df$ using:
$$\eta a^2 d\ell = r^2 df$$

This is the key trick. The result for $\mathcal{W}_1$ (Lara Eq 6):

$$\mathcal{W}_1 = -G\frac{R_\oplus^2}{p^2}\frac{1}{2}\left[B_0\phi + \sum_{i=0}^{1}B_i\sum_{j=\max(i,1)}^{2i+1}\frac{(2-j^*)^i}{j}e^{|j-2i|}\sin(jf+2ig)\right] + \mathcal{C}_1$$

where $\phi = f - \ell$ is the equation of the center.

**Key observation**: The $B_0\phi$ term contains secular growth (since $\phi$ grows with time
at a rate different from the mean anomaly). This is the term that produces
the secular perturbation in $\ell$. The sin terms are purely periodic.

**Connection to BH61**: This IS Brouwer's $S_1$ (Eq 9 of BH61, referenced in Eq 11
and used to compute $\delta p_j$ via operator $D$). The $\phi$ term corresponds to
Brouwer's $(f-l)$ terms in Eq (14), and the sin terms correspond to the
$\sin(2g+f)$, $\sin(2g+3f)$ etc. terms.

### Step 5: The integration constant C₁ (Lara Eq 13, p.9)

$$\mathcal{C}_1 = G\frac{R_\oplus^2}{p^2}\frac{15s^2-14}{32(5s^2-4)}s^2e^2\sin 2g$$

**This is Brouwer's $S_1^*$.** It's chosen specifically to cancel the long-period
terms that would otherwise appear at second order. The $(5s^2-4) = (1-5\cos^2 I)$
in the denominator is the critical inclination divisor.

**Connection to BH61**: $S_1^*$ appears in BH61 Eqs (9) and (11). It's the piece that
Brouwer separates into the "elimination of the perigee" transformation.
The $\sin 2g$ dependence produces the long-period terms in $e$ and $I$.

### Step 6: First-order periodic corrections (Lara Eq 15, p.10)

The periodic corrections to any element $\xi$ are:

$$\Delta\xi = J_2\{\xi, \mathcal{W}_1\}$$

where $\{,\}$ denotes the Poisson bracket. For example, $\Delta a$:

$$\Delta a = a\frac{R_\oplus^2}{p^2}\frac{1}{4\eta^2}\sum_{i=0}^{1}B_i(s)\sum_{j=-i}^{3+2i}A_{i,j}(\eta)e^{|j-2i|}\cos(jf+2ig)$$

with explicit $A_{i,j}$ coefficients listed.

**Connection to BH61**: These are the $\delta a$, $\delta e$, $\delta I$, $\delta l$ formulas
on pp.200-201 of BH61. The Poisson bracket $\{\xi, \mathcal{W}_1\}$ is equivalent to
BH61's operator $D$ acting on $\partial(S_1+S_1^*)/\partial l_j''$ (Eq 11).

## Pages 6-10 Detailed Notes

### Page 6: Explicit W₁ computation

Eq (5) rearranges the homological equation using $\phi = f - \ell$:
$$\mathcal{W}_1 = \frac{1}{n}\left[\mathcal{H}_{0,1}\phi + \int\left(\tilde{\mathcal{H}}_{0,1}\frac{r^2}{a^2\eta} - \mathcal{H}_{0,1}\right)df\right] + \mathcal{C}_1$$

The integrand in Eq (5) contains ONLY periodic functions of $f$ (the constant
parts cancel by construction). This means the integral produces only $\sin$ terms.

Eq (6) gives the explicit result. The first term $B_0\phi$ is Brouwer's "equation
of the center" contribution. The sum gives the short-period oscillations.

**Key insight**: The $B_0\phi$ term in $\mathcal{W}_1$ does NOT produce a short-period
correction itself — but when $D$ acts on it (BH61 language), the derivative
$d\phi/dl$ produces secular + constant terms. This is why short-period terms
in $S_1$ produce secular/long-period terms in $\delta p_j$ (BH61 p.198).

### Page 7: Second-order computation setup

Eq (7) gives $\tilde{\mathcal{H}}_{0,2}'$ — the second-order known terms from Poisson brackets
$\{\mathcal{H}_{1,0}, \mathcal{W}_1\}$. This has three structural blocks:
1. Terms free from $\phi$ and factored by $a^2/r^2$ — trivial integration in $f$
2. Terms free from both $\phi$ and $r$ — integration of $\cos$ functions in elliptic motion
3. Terms of the form $(p/r)^2\phi\sin(mf+\alpha)$ — integrated by parts using Eq (8)

Table 1 lists all the inclination polynomials $B_{i,j,k}(s)$ needed.

Eq (8) is the integration-by-parts formula:
$$\frac{m}{\eta^3}\int\frac{p^2}{r^2}\phi\sin mf\,d\ell = -\phi\cos mf + \frac{\sin mf}{m} - \int\cos mf\,d\ell$$

**Connection to BH61**: This machinery is what BH61 uses implicitly when computing
the variation terms $(p_j)_x$, $(p_j)_e$ in Sections I-VI. The integration-by-parts
produces the mixed secular terms ($nl \cdot \cos kl$, $nl \cdot \sin kl$) that appear in
BH61's equations for $dL''/dt$ etc.

### Page 8: Integration constant and long-period elimination

Eq (9) gives $\tilde{\mathcal{H}}_{0,2}^*$ — the Poisson bracket $\{\mathcal{H}_{1,0}+\mathcal{H}_{0,1}, \mathcal{C}_1\}$.
This contains the derivatives $\partial\mathcal{C}_1/\partial g$, $\partial\mathcal{C}_1/\partial G$, $\partial\mathcal{C}_1/\partial L$.

Table 2 lists inclination polynomials $b_{i,j,k}(s)$ and eccentricity polynomials $q_{i,j}(e)$.

### Page 9: Brouwer's second-order Hamiltonian recovered

Eq (10) gives the averaging rule:
$$\frac{1}{2\pi}\int_0^{2\pi}\cos(mf+\alpha)\,d\ell = \left(\frac{-e}{1+\eta}\right)^m(1+m\eta)\cos\alpha$$

This is the key formula for averaging products of $\cos(mf)$ over the mean anomaly.

Eq (11): After averaging, $\langle\tilde{\mathcal{H}}_{0,2}'\rangle_\ell$ gives Brouwer's second-order
Hamiltonian (secular + long-period parts). The secular part matches Brouwer exactly.

Eq (12): $\langle\tilde{\mathcal{H}}_{0,2}^*\rangle_\ell$ gives the long-period contribution from $\mathcal{C}_1$.

Eq (13): $\mathcal{C}_1$ is CHOSEN so that Eq (12) cancels the long-period part of Eq (11).
This determines $\mathcal{C}_1$ uniquely (up to a pure constant).

Eq (14): The complete second-order Hamiltonian $\mathcal{H}_{0,2}$ matches Brouwer's result
exactly, confirming the equivalence of the single-transformation approach.

**Connection to BH61**: The $\mathcal{C}_1 \sim \sin 2g$ dependence is what produces the
long-period perturbations in $e$, $I$, and $g$ through the Poisson bracket
$\{\xi, \mathcal{C}_1\}$. This is the $S_1^*$ in BH61's Eq (11): $\delta q_j = D(\partial(S_1+S_1^*)/\partial L_j'')$.

### Page 10: First-order periodic corrections

Eq (15): $\Delta a$ — first-order periodic correction to semi-major axis.
This is computed from $\Delta a = J_2\{a, \mathcal{W}_1\} = J_2 \cdot 2a/L \cdot \partial\mathcal{W}_1/\partial\ell$.

The $A_{i,j}(\eta)$ coefficients are listed explicitly:
- $A_{0,0} = 10-6\eta^2-4\eta^3$
- $A_{0,1} = A_{1,1} = A_{1,3} = 15-3\eta^2$
- etc.

**Connection to BH61**: $\Delta a$ from Lara Eq (15) should match BH61's $\delta a$ formula
on p.200. The Lara notation uses $B_i(s)$ and $A_{i,j}(\eta)$ while BH61 uses
$(-1/2+3\theta^2/2)$ and closed-form $(a^3/r^3)$ expressions. These are equivalent
but in different representations (Fourier-expanded vs. closed-form).

## Summary: What This Unlocks

1. **How S₁ is built**: Integrate (perturbation − its average) over mean anomaly.
   The integration uses $\eta a^2 d\ell = r^2 df$ to convert to true anomaly integrals.

2. **How S₁* (C₁) is determined**: Choose it to cancel long-period terms at
   second order. This fixes the $\sin 2g$ coefficient uniquely.

3. **How δp_j, δq_j arise**: Poisson brackets $\{\xi, \mathcal{W}_1\}$. The operator $D$
   in BH61 is equivalent to the specific Poisson bracket structure.

4. **Why short-period → secular**: The equation of the center $\phi = f-\ell$
   in $\mathcal{W}_1$ has a derivative $d\phi/d\ell$ that produces constant terms.
   When $D$ acts on $\phi\cos(2g)$ terms, it generates secular contributions.

5. **The eccentricity function machinery**: Powers of $1/r$ are expanded using
   $(1+e\cos f)^n$, which gives explicit $e^k\cos(jf)$ terms. Combined with
   Kepler's equation, this produces the Fourier-eccentricity double series.

## Pages 11-17: Second-order and third-order secular terms

### Page 11: Table 3 — inclination polynomials for second-order generating function
Massive table of $\beta_{i,j,k}(s)$ polynomials in $s = \sin I$, some with terms up to $s^{10}$.
These feed into the second-order generating function $\mathcal{V}_2$ (Lara Eq 16).

The third-order Hamiltonian computation follows the same pattern:
$\tilde{\mathcal{H}}_{0,3} = \{\mathcal{H}_{0,2}+\mathcal{H}_{1,1}, \mathcal{W}_1\} + \{\mathcal{H}_{0,1}+2\mathcal{H}_{1,0}, \mathcal{W}_2\}$

### Page 12: Integration constant $\mathcal{C}_2$ and third-order Hamiltonian

Eq (17): $\langle\tilde{\mathcal{H}}_{0,3}'\rangle_\ell$ — the averaged third-order terms, with inclination
polynomials from Table 4 (up to $s^{10}$, $\eta^k$ terms).

Eq (18): $\langle\tilde{\mathcal{H}}_{0,3}^*\rangle_\ell = -\mathcal{H}_{0,0}\frac{R_\oplus^2}{p^2}\frac{9}{2}(5s^2-4)\frac{1}{L}\frac{\partial\mathcal{C}_2}{\partial g}$

Eq (19): $\mathcal{C}_2$ chosen to cancel long-period terms, involves $\sin 2ig/(2i)$ structure.

Eq (20): The completely reduced third-order Hamiltonian $\mathcal{H}_{0,3}$.

### Page 13: Second-order periodic corrections (Eq 21)

$$\delta a = a\frac{R_\oplus^4}{p^4}\frac{1}{4^4\eta^4}\left[24\eta^7(5s^4+8s^2-8)+48\eta^5(15^2-14)s^2e^2\cos 2g + \sum\sum\sum A_{i,j,k}\ldots\cos(jf+2ig)\right]$$

Table 5 lists the $A_{i,j,k}$ coefficients — hundreds of integer values, some with
$(3s^2-2)$ factors. This is the second-order $\delta a$ in closed eccentricity form.

**Key insight**: At second order, direct and inverse transformations differ.
The inverse (osculating → mean) uses $\delta\xi = \{\Delta\xi, \mathcal{W}_1\} + \{\xi, -\mathcal{W}_2\}$
evaluated in original variables, not prime variables.

### Page 14: Initialization problem (Section 4)

**Critical practical issue**: Brouwer's periodic corrections are first-order,
but secular frequencies are computed through the reduced Hamiltonian to
second or third order. This mismatch causes initialization errors that
grow secularly.

Eq (22): The energy conservation equation $E_0 = -\mu/(2L'^2) + \sum J_2^m \mathcal{H}_{0,m}/m!$

Eq (23): Breakwell and Vagners' calibration — replace $L'$ with $\hat{L}$ computed
by solving the energy equation. This makes secular frequency initialization
consistent to arbitrary order without computing higher-order periodic corrections.

**Connection to BH61**: This is NOT addressed in BH61. The drag problem compounds
the initialization issue because the mean elements change secularly.

### Page 15-16: Performance results (Figure 1)

{1:2:1} = first-order inverse, second-order secular, first-order direct: ~2.5 km/month error
{1+:2:1} = calibrated initialization: ~20 m/month
{2:2:2} = full second-order: ~1 m/month
{2+:3:2} = third-order secular calibration: ~few cm/month

Single transformation evaluates periodic corrections 2x faster than
the three-transformation (parallax→perigee→Delaunay) approach.

### Page 17: Conclusions

Brouwer's splitting into three transformations is conventional, not necessary.
A single transformation achieves the same result with better compiler optimization
because common subexpressions like $(5s^2-4)$ appear ~300 times (vs 73 in the split version).

## Assessment: Does This Paper Unlock BH61?

### What it DOES unlock:
1. **How S₁ is constructed** — Eq (2) homological equation + Eq (5)-(6) explicit integration ✓
2. **How S₁* (C₁) is determined** — Eq (13), cancel long-period terms at next order ✓
3. **How periodic corrections arise** — Poisson brackets $\{\xi, \mathcal{W}_1\}$ = Eq (15) ✓
4. **The eccentricity function expansion** — Eq (3)-(4), converting $r^{-j}$ to Fourier ✓
5. **The averaging rule** — Eq (10), mean of $\cos(mf+\alpha)$ over $\ell$ ✓
6. **The integration-by-parts trick** — Eq (8), handling $\phi\sin(mf)$ integrals ✓
7. **Why the equation of center produces secular terms** — the $B_0\phi$ term in $\mathcal{W}_1$ ✓

### What it does NOT cover (gaps remaining for BH61):
1. **The drag problem** — Lara treats only conservative J₂. BH61 adds the
   non-conservative drag force requiring the $V\exp(-\alpha r)$ expansion (Eqs 16-19)
   and the drag-oblateness coupling terms (Sections I-VI).

2. **Higher harmonics J₃, J₄** — Lara focuses on J₂. The Δ₃ and Δ₄ terms in
   BH61 Eq (14') come from extending the Legendre expansion to odd ($P_3$)
   and even ($P_4$) harmonics. The machinery is the same but the inclination
   functions differ.

3. **Integration of non-conservative equations** — BH61 Eqs (23)-(32) integrate
   $dL''/dt$, $dl''/dt$ etc. as ODEs. In the conservative problem this step
   doesn't exist (the transformation IS the solution). The drag problem
   requires actual time integration of the mean element equations.

4. **The operator D** — Lara uses Poisson brackets directly. BH61's operator
   $D = -\sum \xi_j'' \partial/\partial\xi_j''$ is an equivalent but different formulation.
   The connection is: $D$ acting on a Kepler function equals the Poisson
   bracket with $\mathcal{W}_1$ for the specific case of velocity-proportional drag.

## Next Steps

1. **Use Lara's Eq (6)** to write out $\mathcal{W}_1$ (= $S_1$) explicitly, then compute
   $\partial\mathcal{W}_1/\partial g$ and $\partial\mathcal{W}_1/\partial G$ to verify BH61's $\delta p_j$, $\delta q_j$ in Eq (14).

2. **Use Lara's Eq (15)** to verify BH61's $\delta a$, $\delta e$, $\delta I$, $\delta l$ on pp.200-201.
   The $A_{i,j}(\eta)$ coefficients from Lara should map to BH61's closed-form expressions.

3. **Extend to J₃, J₄** by adding $P_3$ and $P_4$ Legendre terms to the Hamiltonian.
   This would let us verify the Δ₃ and Δ₄ terms in BH61 Eq (14'), including
   the disputed signs.

4. **Consider acquiring Lara's textbook** (De Gruyter 2021) for the full
   pedagogical treatment including worked examples and the drag problem context.
