# Derivation 015: Low-Precision Solar Ephemeris

## Purpose

Compute the Sun's geocentric position (ecliptic longitude and distance)
at a given Julian date. This is needed for the deep-space third-body
perturbation — we need the Sun's position to compute its gravitational
effect on the satellite.

We do NOT need JPL DE440 precision. The solar perturbation on the satellite
is a small correction (meters to tens of meters per orbit), so the Sun's
position only needs to be known to ~0.01° (~1 arcminute). A low-precision
analytical ephemeris is sufficient and much faster than table interpolation.

## The Sun's Apparent Orbit

From Earth's perspective, the Sun appears to orbit Earth on the ecliptic
plane. This is just Earth's orbit viewed from the geocentric frame.

The Sun's ecliptic longitude λ☉ advances from the vernal equinox at a
rate of ~360°/365.25 days = 0.9856°/day. The orbit is slightly elliptical
(e = 0.01671), so the angular rate varies through the year.

## Step 1: Mean Anomaly

The Sun's mean anomaly M☉ advances at a constant rate:

$$M_☉ = M_0 + n_☉ (t - t_0)$$

where:
- n☉ = 2π / T☉ (mean motion, T☉ = 365.256363004 days = tropical year)
- M₀ = mean anomaly at epoch
- t₀ = reference epoch (typically J2000.0)

In SGP4, the solar rates are already in FundamentalConstants:
- ZNS = 1.19459e-5 rad/min (solar mean anomaly rate)
- Solar anomaly epoch: ZMOS from [SR3] page 59

So: M☉(t) = ZMOS + ZNS × (t − t_epoch) where t is in minutes.

## Step 2: Equation of Center

The true anomaly ν differs from the mean anomaly by the equation of center:

$$\nu = M + 2e\sin M + \frac{5}{4}e^2\sin 2M + \frac{13}{12}e^3\sin 3M + \ldots$$

For the Sun with e = 0.01671, the series converges very rapidly:
- First term: 2e sin M ≈ 0.0334 rad ≈ 1.915° (max)
- Second term: (5/4)e² sin 2M ≈ 0.00035 rad ≈ 0.020°
- Third term: negligible (~0.0003°)

Two terms give 0.02° accuracy, more than sufficient for SGP4 deep-space.

## Step 3: Ecliptic Longitude

$$\lambda_☉ = \nu + \tilde{\omega}$$

where ω̃ is the longitude of perihelion. This advances slowly due to
planetary perturbations (~0.00005°/year), but for SGP4 purposes it can
be treated as a constant (~282.94° at J2000.0).

In the SGP4 formulation, the solar position is expressed through
its ecliptic longitude directly, computed from the mean longitude and
equation of center.

Mean longitude: L☉ = M☉ + ω̃

The solar ecliptic longitude: λ☉ = L☉ + equation_of_center

## Step 4: Geocentric Distance

From the orbit equation:

$$r_☉ = \frac{a(1-e^2)}{1 + e\cos\nu}$$

where a = 1 AU = 149597870.7 km, e = 0.01671.

For the perturbation calculation, we need r in Earth radii or AU.
The ratio r/r☉ that appears in the third-body disturbing function uses
a☉ = 1 AU as the normalizing distance.

## Step 5: Ecliptic to Equatorial Transformation

The perturbation acts in 3D, so we need the Sun's position in the
equatorial frame (or at least its declination/right ascension).

The transformation from ecliptic (λ, β=0) to equatorial (α, δ) is:

$$\sin\delta = \sin\epsilon\sin\lambda$$
$$\cos\alpha\cos\delta = \cos\lambda$$
$$\sin\alpha\cos\delta = \cos\epsilon\sin\lambda$$

where ε ≈ 23.44° is the obliquity.

For the SGP4 deep-space model, the Sun's position is expressed in terms
of sin(λ), cos(λ), combined with sin(ε), cos(ε) (which are the ZSINIS
and ZCOSIS constants already in FundamentalConstants).

## Summary

The solar ephemeris for SGP4 needs:
1. Mean anomaly: M☉(t) = ZMOS + ZNS × Δt
2. Equation of center: C = 2e sin M + (5/4)e² sin 2M
3. Ecliptic longitude: λ = M + ω̃ + C
4. Distance: r = a(1-e²)/(1+e cos ν) (or approximate as a(1-e cos M))
5. Convert to equatorial using obliquity ε

The FundamentalConstants struct already holds ZNS, ZES (eccentricity),
ZCOSIS, ZSINIS (obliquity), and the solar argument of perigee.

## Equations to Resolve

| Equation | Status | What must be shown |
|---|---|---|
| Equation of center coefficient 2e | NOT DERIVED | First-order Lagrange inversion of Kepler equation |
| Equation of center coefficient (5/4)e² | NOT DERIVED | Second-order Lagrange inversion |
| Equation of center coefficient (13/12)e³ | NOT DERIVED | Third-order Lagrange inversion (if needed for accuracy) |
| Distance formula r ≈ a(1-e·cos M) | CITED | First-order Kepler approximation |
| Convergence radius (Laplace limit e < 0.6627) | NOT DERIVED | Singularity of Lagrange inversion |

## Source Documents Required

| Source | Location | Availability | What it provides |
|---|---|---|---|
| Self-contained math (Lagrange inversion theorem) | N/A | ✓ | No external source strictly needed |
| Meeus (1998) Ch. 25 | External | NOT in repo | Cross-check reference |
| Brouwer & Clemence (1961) | External | NOT in repo | Lagrange inversion in celestial mechanics context |
