# Spherical Harmonics and Their Relationship to the J Terms

## 1. The Gravitational Potential Expansion

The gravitational potential of the Earth at a point $(r, \varphi', \lambda)$ outside the Earth is:

$$V(r, \varphi', \lambda) = \frac{GM}{r}\left[1 + \sum_{n=2}^{\infty}\left(\frac{a}{r}\right)^n \sum_{m=0}^{n} (\bar{C}_{nm}\cos m\lambda + \bar{S}_{nm}\sin m\lambda)\bar{P}_{nm}(\sin\varphi')\right]$$

where:
- $(r, \varphi', \lambda)$ = geocentric spherical coordinates (radius, latitude, longitude)
- $\bar{C}_{nm}$, $\bar{S}_{nm}$ = fully normalized spherical harmonic coefficients
- $\bar{P}_{nm}(\sin\varphi')$ = fully normalized associated Legendre functions
- $n$ = degree, $m$ = order

## 2. Building Blocks: Legendre Polynomials

### 2.1 Legendre Polynomials $P_n(x)$

The Legendre polynomial of degree $n$ is defined by Rodrigues' formula:

$$P_n(x) = \frac{1}{2^n n!}\frac{d^n}{dx^n}(x^2 - 1)^n$$

First few:

$$P_0(x) = 1, \quad P_1(x) = x, \quad P_2(x) = \frac{3x^2 - 1}{2}, \quad P_3(x) = \frac{5x^3 - 3x}{2}$$

$$P_4(x) = \frac{35x^4 - 30x^2 + 3}{8}$$

**Recurrence (for computation):**

$$P_{n+1}(x) = \frac{(2n+1)xP_n(x) - nP_{n-1}(x)}{n+1}$$

Starting from $P_0 = 1$, $P_1 = x$.

### 2.2 Associated Legendre Functions $P_{nm}(x)$

$$P_{nm}(x) = (1-x^2)^{m/2}\frac{d^m}{dx^m}P_n(x)$$

or equivalently:

$$P_{nm}(x) = \frac{(1-x^2)^{m/2}}{2^n n!}\frac{d^{n+m}}{dx^{n+m}}(x^2-1)^n$$

**Key property:** $P_{n0}(x) = P_n(x)$ (the associated function reduces to the Legendre polynomial when $m = 0$).

**Recurrence (for computation):**

Sectoral ($m = n$): $P_{mm}(\sin\varphi) = (2m-1)!!\cos^m\varphi$

Then build up $n$ from $m$:

$$P_{n+1,m}(\sin\varphi) = \frac{(2n+1)\sin\varphi \cdot P_{nm}(\sin\varphi) - (n+m)P_{n-1,m}(\sin\varphi)}{n-m+1}$$

### 2.3 Normalization

The **fully normalized** associated Legendre functions (geodetic convention) are:

$$\bar{P}_{nm}(x) = \sqrt{\frac{(n-m)!(2n+1)k}{(n+m)!}} P_{nm}(x)$$

where $k = 1$ if $m = 0$, $k = 2$ if $m > 0$.

The normalization ensures:

$$\frac{1}{4\pi}\int_0^{2\pi}\int_{-\pi/2}^{\pi/2} [\bar{P}_{nm}(\sin\varphi)\cos m\lambda]^2 \cos\varphi\,d\varphi\,d\lambda = 1$$

The **unnormalized** (conventional) coefficients relate to normalized by:

$$\bar{C}_{nm} = \sqrt{\frac{(n+m)!}{(n-m)!(2n+1)k}} C_{nm}$$

## 3. Zonal Harmonics: The $J_n$ Terms

### 3.1 Definition

When $m = 0$, the spherical harmonic depends only on latitude (no longitude dependence). These are called **zonal harmonics**. The gravitational potential from zonals alone is:

$$V_\text{zonal} = \frac{GM}{r}\left[1 - \sum_{n=2}^{\infty}J_n\left(\frac{a}{r}\right)^n P_n(\sin\varphi')\right]$$

where the **conventional (unnormalized) zonal coefficients** are:

$$J_n = -C_{n0}$$

and the relationship to the **normalized** coefficients is:

$$\bar{C}_{n0} = \frac{C_{n0}}{\sqrt{2n+1}} = \frac{-J_n}{\sqrt{2n+1}}$$

Or equivalently, solving for $J_n$:

$$J_n = -\bar{C}_{n0} \cdot \sqrt{2n+1}$$

**Derivation:** The potential term for degree $n$, order $m=0$ must be the same whether using normalized or unnormalized quantities:

$$\bar{C}_{n0} \cdot \bar{P}_{n0}(\sin\varphi') = C_{n0} \cdot P_{n}(\sin\varphi')$$

Since $\bar{P}_{n0} = \sqrt{2n+1} \cdot P_n$ (for $m=0$, $k=1$), we get $\bar{C}_{n0} \cdot \sqrt{2n+1} = C_{n0}$, so $\bar{C}_{n0} = C_{n0}/\sqrt{2n+1}$.

**Numerical verification for $n=2$:**

$$J_2 = -\bar{C}_{20} \cdot \sqrt{5} = -(-4.84166774985 \times 10^{-4}) \times 2.2360679\ldots = 1.082629821313 \times 10^{-3} \; \checkmark$$

### 3.2 Specific Values

| Coefficient | $\bar{C}_{n0}$ (normalized) | $J_n$ (conventional) | Relationship |
|---|---|---|---|
| $n=2$ | $\bar{C}_{20} = -4.84166774985 \times 10^{-4}$ | $J_2 = 1.082629821313 \times 10^{-3}$ | $J_2 = -\bar{C}_{20} \cdot \sqrt{5}$ |
| $n=3$ | $\bar{C}_{30}$ | $J_3 = -C_{30}$ | $J_3 = -\bar{C}_{30} \cdot \sqrt{7}$ |
| $n=4$ | $\bar{C}_{40}$ | $J_4$ | $J_4 = -\bar{C}_{40} \cdot 3$ |

### 3.3 Physical Meaning

Each $J_n$ describes a specific latitude-dependent mass distribution:

- $J_2 \approx 1.08 \times 10^{-3}$: **Oblateness** — the Earth is fatter at the equator. By far the dominant term. Causes orbital precession.
- $J_3 \approx -2.5 \times 10^{-6}$: **Pear shape** — slight north-south asymmetry. Causes long-period perturbations in eccentricity and argument of perigee.
- $J_4 \approx -1.6 \times 10^{-6}$: **Octupole correction** to oblateness. Causes secular perturbations.
- Higher $J_n$: increasingly fine latitude structure, decreasing magnitude (~$J_n \sim J_2^{n/2}$).

### 3.4 The $J_n$ for an Equipotential Ellipsoid

For a level ellipsoid (like GRS80/WGS84), ALL even zonal harmonics $J_{2n}$ are determined by just $J_2$ (or equivalently $e^2$) and $m = \omega^2 a^2 b/GM$:

$$J_{2n} = (-1)^{n+1}\frac{3e^{2n}}{(2n+1)(2n+3)}\left(1 - n + 5n\frac{J_2}{e^2}\right)$$

The odd zonals ($J_3$, $J_5$, ...) are **zero** for a symmetric ellipsoid. Any nonzero odd $J_n$ indicates the real Earth's deviation from a perfect ellipsoid of revolution.

## 4. Tesseral and Sectoral Harmonics

### 4.1 Classification

| Type | Condition | Longitude Dependence | Example |
|---|---|---|---|
| **Zonal** | $m = 0$ | None (latitude bands) | $J_2$, $J_3$, $J_4$ |
| **Sectoral** | $m = n$ | $\cos(n\lambda)$ or $\sin(n\lambda)$ (longitude sectors) | $C_{22}$, $S_{22}$ |
| **Tesseral** | $0 < m < n$ | Both lat and lon (checkerboard pattern) | $C_{21}$, $S_{31}$ |

### 4.2 The $C_{22}$, $S_{22}$ Terms (Earth's Ellipticity)

The second-degree sectoral harmonics describe the non-axisymmetric part of the geopotential at degree 2:

$$V_{22} = \frac{GM}{r}\left(\frac{a}{r}\right)^2 (C_{22}\cos 2\lambda + S_{22}\sin 2\lambda)\bar{P}_{22}(\sin\varphi')$$

These create the "tesseral bulge" that causes resonance effects on 12-hour and 24-hour satellites. The SGP4 resonance constants Q22, ROOT22, etc. derive from these coefficients.

### 4.3 Relationship to SGP4 Resonance Constants

The SGP4 resonance mechanism considers the interaction between a satellite's orbital frequency and the Earth's rotation-dependent tesseral harmonics. For a satellite with mean motion $n$ near $k$ revolutions per sidereal day:

$$Q_{nm} \propto \sqrt{C_{nm}^2 + S_{nm}^2}$$

Specifically:
- $Q_{22} \propto J_{2,2} = \sqrt{C_{22}^2 + S_{22}^2}$ — the 24-hour (synchronous) resonance
- $Q_{31}$, $Q_{33}$ — 12-hour resonance contributions
- ROOT values are $\sqrt{Q_{nm} \times Q_{n'm'}}$ cross-coupling terms

## 5. Generating Spherical Harmonics: Computation

### 5.1 The Recurrence for $\bar{P}_{nm}$

For arbitrary-precision computation, the **stable** recurrence is:

**Sectoral seed ($m = n$):**

$$\bar{P}_{nn}(\sin\varphi) = \cos\varphi \cdot \sqrt{\frac{2n+1}{2n}} \cdot \bar{P}_{n-1,n-1}(\sin\varphi)$$

starting from $\bar{P}_{00} = 1$.

**Tesseral recurrence ($n > m$):**

$$\bar{P}_{nm}(\sin\varphi) = \sin\varphi\cdot a_{nm}\cdot\bar{P}_{n-1,m}(\sin\varphi) - b_{nm}\cdot\bar{P}_{n-2,m}(\sin\varphi)$$

where:

$$a_{nm} = \sqrt{\frac{(2n-1)(2n+1)}{(n-m)(n+m)}}, \qquad b_{nm} = \sqrt{\frac{(2n+1)(n+m-1)(n-m-1)}{(2n-3)(n-m)(n+m)}}$$

These recurrence coefficients are **exact algebraic expressions** in integers under the square root. For arbitrary precision, compute them to full precision from the integer formula.

### 5.2 Stability

The "column-first" recurrence (increasing $n$ for fixed $m$) is numerically stable. The "row-first" approach is unstable for large $n$. Always use column-first.

### 5.3 Derivative for Gravity Computation

For computing gravity (the gradient of $V$), we need $d\bar{P}_{nm}/d\varphi$:

$$\frac{d\bar{P}_{nm}}{d\varphi} = n\sin\varphi\cdot\bar{P}_{nm} - \sqrt{(n^2-m^2)(2n+1)/(2n-1)}\cdot\bar{P}_{n-1,m}$$

divided by $\cos\varphi$ (with special handling at the poles).

## 6. Connection to Our Design

### Where Spherical Harmonics Appear

1. **EquipotentialEllipsoid** — the $J_{2n}$ even zonal coefficients define the ellipsoid's gravity field
2. **GravitationalModel** — stores $\bar{C}_{nm}$, $\bar{S}_{nm}$ and evaluates $V$ via the spherical harmonic sum
3. **GeoidModel** — computes geoid undulation from the spherical harmonic expansion
4. **SGP4 Resonance** — the tesseral $Q_{nm}$ values come from $\bar{C}_{nm}$, $\bar{S}_{nm}$

### What We Need to Implement

| Component | Needed For | Priority |
|---|---|---|
| Normalized Legendre recurrence | Gravitational potential evaluation | P1 |
| $J_n$ from $\bar{C}_{n0}$ conversion | Ellipsoid-to-SGP4 constant bridge | P0 |
| $Q_{nm}$ from $\bar{C}_{nm}$, $\bar{S}_{nm}$ | Re-deriving SGP4 resonance constants | P2 |
| Clenshaw summation for Legendre series | Efficient spherical harmonic evaluation | P1 |
| Full gradient computation | Gravity vector from EGM | P2 |
