# Documentation Gaps: Undocumented Formulas and Constants in SGP4

## Coverage Summary

| Category | Documented In | Coverage |
|---|---|---|
| WGS84 defining parameters | [NGA] Ch. 3 | Complete |
| WGS84 derived geometric constants | [NGA] App. B, [M80] Sec. 3 | Complete |
| WGS84 derived physical constants | [NGA] App. B, [M80] Sec. 3 | Complete |
| Somigliana gravity formula | [NGA] Ch. 4, [M80] Sec. 3, [HM67] | Complete |
| ECI-ECEF transformation | [NGA] App. A, [IERS10] | Complete |
| Brouwer secular perturbation polynomials | Brouwer (1959) | Need separate reference |
| Kaula inclination functions (F_lmp) | Kaula (1966) | Need separate reference |
| Hansen/Kaula eccentricity functions (G_lpq) | Hough curve fits | **OPAQUE — ~100 magic numbers** |
| Tesseral resonance constants (Q22, ROOT22...) | Unknown 1970s gravity model | **OPAQUE** |
| Resonance phase angles (G22, FASX2...) | Epoch-frozen, unknown derivation | **OPAQUE** |
| Solar/lunar ephemeris constants | ~1970s simplified ephemeris | **PARTIALLY OPAQUE** |
| Deep space secular coefficients | Third-body perturbation theory | Partially traceable |
| GMST formula | Aoki et al. (1982) | Documented but SUPERSEDED by IERS 2010 |

## Additional References Needed

To close the documentation gaps, we need these additional source documents:

| Document | What It Covers |
|---|---|
| **Brouwer (1959)** "Solution of the Problem of Artificial Satellite Theory Without Drag" *Astron. J.* 64, 378-397 | Secular rate polynomials (13, 78, 137, etc.) and short-period correction formulas |
| **Lyddane (1963)** "Small Eccentricities or Inclinations in the Brouwer Theory of the Artificial Satellite" *Astron. J.* 68, 555-558 | The Lyddane modification for small inclination singularity |
| **Kaula (1966)** *Theory of Satellite Geodesy* | Inclination functions F_lmp and eccentricity functions G_lpq (exact forms) |
| **Hoots & Roehrich (1980)** Spacetrack Report #3 | The distributed FORTRAN source, including the deep space code |
| **Hough (1981)** Deep space theory extension | The G-function polynomial curve fits (if available) |
| **Lane & Hoots (1979)** "General Perturbations Theories Derived from the 1965 Lane Drag Theory" | Atmospheric drag model used in SGP4 |
| **Aoki et al. (1982)** "The New Definition of Universal Time" *Astron. & Astrophys.* 105 | The GMST formula used in the SGP4 sidereal time computation |

## Opaque Constants: Strategy for Our Implementation

### Option A: Re-derive from First Principles

For the G-function polynomials (Hansen coefficients) and F-functions (Kaula inclination functions):

1. The exact Hansen coefficients $G_{lpq}(e)$ are defined as:
   $$G_{lpq}(e) = \frac{1}{2\pi}\int_0^{2\pi} \left(\frac{r}{a}\right)^l \cos[(l-2p+q)v - (l-2p)E] \, dM$$
   where $v$ = true anomaly, $E$ = eccentric anomaly, $M$ = mean anomaly.
   These can be evaluated to arbitrary precision using Bessel function series.

2. The exact Kaula inclination functions $F_{lmp}(i)$ are:
   $$F_{lmp}(i) = \sum_t \binom{2l-2t}{l}\binom{l}{t}\binom{l-t}{m+l-2p-2t} \frac{(\sin i)^{l-m-2t} \cdot (\cos i)^{...}}{2^l \cdot l!} \cdots$$
   (Kaula's recursion). These are exact algebraic expressions in $\sin i$ and $\cos i$.

**Benefit:** Full traceability, arbitrary precision, proper error bounds.
**Cost:** Significant mathematical development effort.

### Option B: Adopt Published Polynomial Fits with ModelValue Tracking

Keep the existing polynomial approximations but:
1. Document them as `ModelValue` with provenance "Hough curve fits to Hansen coefficients"
2. Assign accuracy bounds based on the degree of the polynomial approximation
3. Flag them in the error budget as "model-limited, not precision-limited"

**Benefit:** Quick to implement, compatible with standard SGP4 validation.
**Cost:** Cannot improve accuracy beyond the polynomial approximation; opaque provenance remains.

### Option C: Hybrid

Use Option B initially for SGP4 compatibility, then implement Option A as a future enhancement phase (Phase 7 from the original plan: "Research Higher-Precision References").

**This is the recommended approach.** It lets us build and validate the framework now while planning a principled re-derivation later.

## Specific Magic Numbers Decoded

### F-Coefficient Rational Values

| Code Value | Exact Rational | Verified? |
|---|---|---|
| 0.75 | $3/4$ | Yes |
| 1.5 | $3/2$ | Yes |
| 1.875 | $15/8$ | Yes |
| 0.9375 | $15/16$ | Yes |
| 35.0 | $35$ | Yes |
| 39.3750 | $315/8$ | Yes |
| 9.84375 | $315/32$ | Yes |
| 4.92187512 | $315/64 = 4.921875$ exactly | **Trailing "12" is FORTRAN artifact** |
| 6.56250012 | $105/16 = 6.5625$ exactly | **Trailing "12" is FORTRAN artifact** |
| 29.53125 | $945/32$ | Yes |
| 0.33333333 | $1/3$ | Truncated decimal of irrational |

The trailing "12" digits on 4.92187512 and 6.56250012 are artifacts from the original FORTRAN code where these rational numbers were entered as decimal approximations. The exact values are $315/64$ and $105/16$. **Our implementation should use the exact rational forms.**

### Brouwer Polynomial Coefficients

| Code Value | Context | Origin |
|---|---|---|
| 13, 78, 137 | $J_2^2$ secular term in $\dot{M}$ | Brouwer (1959) Eq. 38 |
| 7, 114, 395 | $J_2^2$ secular term in $\dot{\omega}$ | Brouwer (1959) Eq. 36 |
| 3, 36, 49 | $J_4$ secular term in $\dot{\omega}$ | Brouwer (1959) Eq. 36 |
| 4, 19 | $J_2^2$ secular term in $\dot{\Omega}$ | Brouwer (1959) Eq. 37 |
| 3, 7 | $J_4$ secular term in $\dot{\Omega}$ | Brouwer (1959) Eq. 37 |
| 134/81 | 3rd-order element recovery | Brouwer (1959) — series reversion |
| 0.0625 = 1/16 | Scaling factor for $J_2^2$ terms | From the expansion normalization |
| -0.46875 = -15/32 | $J_4$ coefficient in temp3 | Brouwer (1959) |

These should all be tracked as exact rationals in our implementation, not as floating-point approximations.

### Synchronous Resonance Constants

| Code Value | Exact? | Context |
|---|---|---|
| 0.8125 = 13/16 | Exact | g200 polynomial |
| 6.60937 ≈ 423/64 = 6.609375 | **Truncated** | g300 polynomial — should be exact 423/64 |

### Deep Space Integration

| Code Value | Meaning |
|---|---|
| 720.0 minutes | Integration step = 12 hours |
| 259200.0 | $720^2 / 2 = 259200$ — precomputed for efficiency |

These are design choices, not physical constants.
