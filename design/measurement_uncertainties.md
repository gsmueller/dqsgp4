# Measurement Uncertainties for MeasuredValue Constants

Source: DMA TR 8350.2 Second Edition (1 September 1991), Table 3.1

## WGS84 Defining Parameters with 1-sigma Accuracy Estimates

| Parameter | Value | Accuracy (1σ) | Status in 2014 NGA Standard |
|---|---|---|---|
| $a$ | $6378137$ m | $\pm 2$ m | **Definitional** (no uncertainty by convention) |
| $\bar{C}_{2,0}$ | $-484.16685 \times 10^{-6}$ | $\pm 1.30 \times 10^{-6}$ | Replaced by $1/f$ as defining param |
| $\omega$ | $7292115 \times 10^{-11}$ rad/s | $\pm 0.1500 \times 10^{-11}$ rad/s | **Definitional** (no uncertainty by convention) |
| $GM$ | $3986005 \times 10^{8}$ m³/s² | $\pm 0.6 \times 10^{8}$ m³/s² | Refined to $3.986004418 \times 10^{14}$ (still measured) |

## Special Parameters with Accuracy Estimates

| Parameter | Value | Accuracy (1σ) | Notes |
|---|---|---|---|
| $GM'$ (atm excluded) | $3986001.5 \times 10^{8}$ | $\pm 0.6 \times 10^{8}$ | Inherits $GM$ uncertainty |
| $\omega^*$ (precessing) | $(7292115.8553 \times 10^{-11} + 4.3 \times 10^{-15} T_U)$ | $\pm 0.1500 \times 10^{-11}$ | Inherits $\omega$ uncertainty |

## Other Constants with Accuracy (Table 3.4)

| Constant | Value | Accuracy | Notes |
|---|---|---|---|
| $c$ (speed of light) | $299792458$ m/s | $\pm 1.2$ m/s | Now defined exactly (1983) |
| $H$ (dynamical ellipticity) | $1/305.4413$ | $\pm 0.0005$ (on $1/H$) | $H = 3.273795 \times 10^{-3}$ |
| $G$ (gravitational constant) | $6.673 \times 10^{-11}$ m³/kg/s² | Not stated | Updated to $6.67428 \times 10^{-11}$ in 2014 |
| $GM_A$ (atmosphere) | $3.5 \times 10^{8}$ m³/s² | Not stated | Updated to $3.4359 \times 10^{8}$ in 2014 |

## Evolution of Defining Parameters

| Edition | Fourth Defining Parameter | $GM$ Value |
|---|---|---|
| DMA TR 8350.2 (1987, 1991) | $\bar{C}_{2,0} = -484.16685 \times 10^{-6}$ | $3986005 \times 10^{8}$ (original) |
| NGA.STND.0036 (2014) | $1/f = 298.257223563$ | $3.986004418 \times 10^{14}$ (refined 1994) |

## Implications for TrackedValue Design

The 1991 DMA document provides the original measurement uncertainties that were later "adopted away" when $a$ and $\omega$ became definitional in 2014. For the `MeasuredValue` type:

- $a$: **Definitional** in WGS84 — use `DefinedValue` with representation error only
- $\omega$: **Definitional** in WGS84 — use `DefinedValue` with representation error only
- $GM$: **Measured** — use `MeasuredValue` with $\sigma_{GM} \approx 8 \times 10^{6}$ m³/s² (from the refined value's uncertainty, which is about $0.008 \times 10^{8}$, tighter than the original $\pm 0.6 \times 10^8$)
- $H$: **Measured** — use `MeasuredValue` with $\sigma_H \approx 5 \times 10^{-7}$ (from $\sigma_{1/H} = 0.0005$)
- $G$: **Measured** — use `MeasuredValue` with $\sigma_G \approx 0.00067 \times 10^{-11}$ (2014 value)

The key insight is that the 2014 standard made $a$, $1/f$, and $\omega$ definitional, eliminating their measurement uncertainty BY CONVENTION — not because they're known perfectly, but because changing them would change the coordinate system itself.
