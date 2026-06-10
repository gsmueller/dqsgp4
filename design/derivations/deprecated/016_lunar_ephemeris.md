# Derivation 016: Low-Precision Lunar Ephemeris

## Purpose

Compute the Moon's geocentric position for SGP4 deep-space perturbations.
The Moon is a stronger perturber than the Sun for high-altitude satellites
(because it's much closer), so this ephemeris is critical.

As with the solar ephemeris, we need ~0.1° accuracy, not DE440 precision.

## The Moon's Orbit

The Moon's orbit is far more complex than the Sun's apparent orbit:

1. **Mean longitude** (L) advances at the sidereal rate: ~13.176°/day
   (one orbit in 27.322 days = sidereal month)

2. **Mean anomaly** (l) advances at the anomalistic rate: ~13.065°/day
   (perigee-to-perigee in 27.555 days = anomalistic month)
   These differ because the Moon's perigee advances.

3. **Argument of latitude** (F) = L − Ω, where Ω is the longitude of
   the ascending node

4. **Mean elongation** (D) = L − L☉, the angular separation from the Sun

The SGP4 FundamentalConstants already hold the relevant rates:
- ZNL = 1.5835218e-4 rad/min (anomalistic rate)
- Lunar sidereal rate: 0.22997 rad/day
- Node regression rate: −2π/6798.38 days
- Perigee advance rate: +2π/3231.50 days

## Step 1: Mean Elements

At time t (minutes from epoch):

- Mean anomaly: l = l₀ + ZNL × t
- Mean longitude: L = L₀ + n_sidereal × t
- Node longitude: Ω = Ω₀ + Ω̇ × t (regression)
- Perigee longitude: ω̃ = ω̃₀ + ω̃̇ × t (advance)

## Step 2: Principal Perturbations

The Moon's ecliptic longitude has several large perturbations beyond
the basic Keplerian orbit:

1. **Equation of center**: 6.289° sin(l) (largest term)
2. **Evection**: 1.274° sin(2D − l) (Sun's effect on eccentricity)
3. **Variation**: 0.658° sin(2D) (Sun's effect on longitude)
4. **Annual equation**: 0.214° sin(M☉) (solar distance effect)

For SGP4 (0.1° accuracy sufficient), only the equation of center and
possibly evection matter. The SGP4 formulation uses a simplified lunar
model that captures the dominant terms.

## Step 3: Latitude

The Moon's ecliptic latitude β is dominated by:

$$\beta \approx 5.128° \sin(F)$$

where F = L − Ω is the argument of latitude. The Moon's orbit is
inclined ~5.145° to the ecliptic.

## Step 4: Distance

The Moon's geocentric distance varies between ~356,000 and ~407,000 km
(mean ~384,400 km). The variation is primarily from eccentricity
(e ≈ 0.0549):

$$r = a(1 - e\cos l) \approx 384400(1 - 0.0549\cos l) \text{ km}$$

## Summary for SGP4

The SGP4 deep-space model needs:
1. Lunar mean anomaly: l(t) from ZNL rate
2. Lunar mean longitude: L(t) from sidereal rate
3. Lunar node longitude: Ω(t) from node regression rate
4. Equation of center: ~6.3° sin(l)
5. Ecliptic longitude: λ = L + equation_of_center (simplified)
6. Ecliptic latitude: β ≈ 5.13° sin(F)
7. sin(λ), cos(λ), sin(β), cos(β) for the perturbation integrals

## Equations to Resolve

| Equation | Status | What must be shown |
|---|---|---|
| Equation of center 2e_moon·sin(l) ≈ 6.289° | NOT DERIVED | Lagrange inversion with e_moon = 0.0549 |
| Evection coefficient 1.274° | NOT DERIVED | Solar disturbing function R☉ coupling with eccentricity |
| Variation coefficient 0.658° | NOT DERIVED | Direct solar tidal effect on lunar longitude |
| Annual equation coefficient 0.214° | NOT DERIVED | Solar distance modulation of perturbation |
| Ecliptic latitude β = i·sin(F) | CITED | Small-angle approximation of lunar inclination 5.145° |

## Source Documents Required

| Source | Location | Availability | What it provides |
|---|---|---|---|
| Meeus (1998) Ch. 47 | External | NOT in repo | CRITICAL for evection/variation/annual equation coefficients |
| Brown (1896) lunar theory | External | NOT in repo | Original source |
| Chapront-Touzé & Chapront (1988) | External | NOT in repo | Modern lunar theory |
| Self-contained math for equation of center | N/A | ✓ | Lagrange inversion derivation only |
