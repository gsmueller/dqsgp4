# Derivation 018: Orbital Resonance Detection and Modeling

## Purpose

When a satellite's orbital period is commensurate with Earth's rotation
period, tesseral harmonics (m ≠ 0) produce persistent perturbations that
don't average out over many orbits. This module detects and models these
resonances.

## Resonance Types

### Synchronous (period ≈ 1436 minutes)
- GEO satellites (geostationary/geosynchronous)
- The satellite returns to the same longitude every orbit
- J₂₂ (Earth's equatorial ellipticity) drives a libration in mean longitude
- Period of libration: ~800 days for typical GEO

### Half-day (period ≈ 718 minutes)
- GPS-like orbits, Molniya-type
- The satellite passes over the same ground track every 2 revolutions
- J₂₂, J₃₂, J₄₄, J₅₂, J₅₄ all contribute
- More complex dynamics than synchronous

### None (all other periods)
- No resonance — tesseral terms average out
- Standard secular/periodic treatment sufficient

## Detection Criterion

The period threshold is empirical:
- Period > 1200 min AND period < 1800 min → likely synchronous
- Period > 600 min AND period < 800 min → likely half-day
- Otherwise → no resonance

More precisely, check if the mean motion n₀ is close to a rational
multiple of Earth's rotation rate ωₑ:
- n ≈ ωₑ → synchronous (1:1)
- n ≈ 2ωₑ → half-day (2:1)

## Resonance Equations

The resonance variables are the mean longitude relative to Earth's
rotation (the "resonance angle"). The equations of motion for this
angle involve the tesseral harmonic coefficients G_{lmpq} which depend
on the specific (l,m,p) combinations.

For synchronous resonance, the dominant term is:
$$\ddot{\lambda}_r \propto J_{22} \sin(2\lambda_r - 2\lambda_{22})$$

This is a pendulum equation — the satellite librates about one of
two stable longitudes separated by 180°.

## Implementation

The resonance module:
1. Classifies the orbit (none, synchronous, half-day) from the period
2. If resonant: computes the tesseral harmonic driving coefficients
3. Provides an integrator that steps the resonance ODEs forward in time
4. The step size is ~720 minutes (one orbit for synchronous, half orbit for half-day)

## Equations to Resolve

| Equation | Status | What must be shown |
|---|---|---|
| Synchronous resonance condition n ≈ ω_E | NOT DERIVED | From tesseral harmonic potential V₂₂ |
| Half-day resonance condition n ≈ 2ω_E | NOT DERIVED | From V₃₂, V₃₃ tesseral harmonics |
| Pendulum equation φ̈ ∝ J₂₂ sin(2φ) | NOT DERIVED | Derive from tesseral disturbing function |
| Kaula G-function eccentricity polynomials G_{lpq}(e) | NOT DERIVED | ~100 polynomial coefficients from Hansen coefficient integrals |
| Tesseral driving coefficients d2201..d5433 | NOT DERIVED | Kaula F×G products for specific (l,m,p,q) combinations |
| Resonance integration (linear → quadratic) | NOT DERIVED | ODE integration scheme for libration |
| Phase angles G22, G32, G44, G52, G54 | NOT DERIVED | From EGM C_{nm}, S_{nm} via λ_{nm} = (1/m)arctan(S/C) |
| FASX2, FASX4, FASX6 | NOT DERIVED | Synchronous resonance phase combinations |

## Source Documents Required

| Source | Location | Availability | What it provides |
|---|---|---|---|
| Sneeuw (2022) "Dynamic Satellite Geodesy" | `sgp4_references/.../Sneeuw_2022_...LectureNotes.pdf` | ✓ IN REPO | Eq. 6.9: Wigner d-coefficients (closed form); Eq. 6.13: inclination function F_lmk; Eq. 6.15: eccentricity function G_lkq; Eqs. 7.8-7.10: resonance condition and repeat orbit theory; Section 7.3: orbit perturbation spectrum |
| Na et al. (2012) "Detailed Re-derivation..." | `sgp4_references/.../Detailed_Re-derivation_...Kaul.pdf` | ✓ IN REPO | Eqs. 3-13a-f: Lagrange planetary equations; NOTE: corrects one error in Kaula (1966) — must verify |
| Kaula (1966) images | `sgp4_references/kaula_1966/` | ✓ IN REPO (images) | Tables of F_lmp and G_lpq values for cross-checking — UNVERIFIED transcription |
| Hujsak (1979) AIAA 79-136 | `sgp4_references/.../Hujsak_1979_...pdf` | ✓ IN REPO | Tesseral resonance equations |
| EGM2008 coefficients | `wgs84_markdown_parts/` | ✓ IN REPO | For phase angles from C_{nm}, S_{nm} |
| Spacetrack Report No. 3 Section 10 pp. 58-68 | In repo | ✓ | FORTRAN implementation |
