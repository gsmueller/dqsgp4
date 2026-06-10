# Derivation 017: Third-Body Gravitational Perturbation

## Purpose

Derive the secular rate corrections to a satellite's orbital elements
from the gravitational attraction of a third body (Sun or Moon).

## The Disturbing Function

For a satellite at position **r** relative to Earth, and a third body at
position **r₃** relative to Earth, the disturbing function is:

$$R_3 = \mu_3 \left(\frac{1}{|\mathbf{r} - \mathbf{r}_3|} - \frac{\mathbf{r} \cdot \mathbf{r}_3}{|\mathbf{r}_3|^3}\right)$$

The first term is the direct attraction of the third body on the satellite.
The second term is the indirect effect (the third body's attraction on Earth,
which we subtract because we work in the geocentric frame).

## Step 1: Legendre Expansion

Since |r| << |r₃| (the satellite is much closer to Earth than to the Sun/Moon),
expand in Legendre polynomials:

$$\frac{1}{|\mathbf{r} - \mathbf{r}_3|} = \frac{1}{r_3} \sum_{n=0}^{\infty} \left(\frac{r}{r_3}\right)^n P_n(\cos S)$$

where S is the angle between **r** and **r₃**.

The n=0 term is constant (doesn't affect the orbit). The n=1 term cancels
with the indirect term. So the disturbing function starts at n=2:

$$R_3 = \frac{\mu_3}{r_3} \sum_{n=2}^{\infty} \left(\frac{r}{r_3}\right)^n P_n(\cos S)$$

For SGP4, we truncate at n=2 (the quadrupole term):

$$R_3 \approx \frac{\mu_3 r^2}{r_3^3} P_2(\cos S) = \frac{\mu_3 r^2}{2 r_3^3}(3\cos^2 S - 1)$$

## Step 2: Express cos(S) in Orbital Elements

The angle S between the satellite and the third body depends on:
- The satellite's argument of latitude u = ω + f
- The satellite's orbital plane orientation (i, Ω)
- The third body's ecliptic longitude λ₃ and latitude β₃

After expressing both position vectors in the equatorial frame and
taking the dot product, cos(S) contains terms involving:
- sin(u)sin(λ₃), cos(u)cos(λ₃), etc.
- Combined with the inclination i, node Ω, and obliquity ε

## Step 3: Orbit Average

Average R₃ over one satellite revolution (integrate over mean anomaly M
from 0 to 2π). The orbit-averaged disturbing function depends only on
the slow variables (a, e, i, ω, Ω) and the third body's position.

The orbit average of r² P₂(cos S) produces terms that are:
- **Secular**: independent of the satellite's argument of latitude
- **Long-period**: dependent on ω (argument of perigee)
- **Short-period**: dependent on u (averaged out)

## Step 4: Secular Rate Extraction

From the orbit-averaged disturbing function, apply Lagrange's planetary
equations to extract the secular rates:

$$\dot{\Omega}_3 = -\frac{1}{na^2\beta\sin i}\frac{\partial\langle R_3\rangle}{\partial i}$$

$$\dot{\omega}_3 = \frac{\beta}{na^2 e}\frac{\partial\langle R_3\rangle}{\partial e} - \frac{\cos i}{na^2\beta\sin i}\frac{\partial\langle R_3\rangle}{\partial i}$$

$$\dot{M}_3 = -\frac{2}{na}\frac{\partial\langle R_3\rangle}{\partial a}$$

$$\dot{e}_3 = -\frac{\beta}{na^2 e}\frac{\partial\langle R_3\rangle}{\partial\omega}$$  (from the e equation)

$$\dot{i}_3 = \frac{1}{na^2\beta\sin i}\left(\frac{\partial\langle R_3\rangle}{\partial\omega}\cos i - \frac{\partial\langle R_3\rangle}{\partial\Omega}\right)$$

## Step 5: SGP4 Implementation

The SGP4 deep-space subroutine computes these rates for both the Sun and Moon,
using the simplified ephemeris from phases 9-10. The rates are then added to
the Brouwer rates to get the total secular rates:

$$\dot{M}_{total} = \dot{M}_{Brouwer} + \dot{M}_{Sun} + \dot{M}_{Moon}$$

and similarly for the other elements.

The deep-space model also includes perturbations in eccentricity (de/dt)
and inclination (di/dt), which are zero in the near-earth Brouwer-only model.

## Key Parameters

The strength of the perturbation scales as μ₃/r₃³:
- Sun: μ☉/r☉³ ≈ 1.32e20 / (1.496e11)³ ≈ 3.96e-14 /s² (in SI)
- Moon: μ☽/r☽³ ≈ 4.90e12 / (3.844e8)³ ≈ 8.63e-14 /s² (in SI)

The Moon's perturbation is ~2.2× stronger than the Sun's because it's
much closer, even though the Sun is much more massive. This is why
the lunar ephemeris matters more than the solar for perturbation accuracy.

## Summary

The third-body perturbation module:
1. Takes the satellite's orbital elements + a third body's position
2. Computes cos(S) from the geometry
3. Evaluates the orbit-averaged P₂(cos S) term
4. Extracts five secular rate corrections via Lagrange's equations
5. Returns corrections to M, ω, Ω, e, i

The same function works for both Sun and Moon — just pass different
ephemeris inputs and μ₃/r₃³ scaling.

## Equations to Resolve

| Equation | Status | What must be shown |
|---|---|---|
| cos(S) decomposition via rotation matrices | NOT DERIVED | Express cos(S) = A sin u + B cos u with explicit A, B |
| Orbit average ⟨P₂(cos S)⟩ = (3(A²+B²)/2 - 1)/2 | NOT DERIVED | Show averaging integrals explicitly |
| Secular rate dΩ/dt from ∂⟨R₃⟩/∂i | NOT DERIVED | Lagrange equation applied to averaged disturbing function |
| Secular rate dω/dt from ∂⟨R₃⟩/∂e, ∂⟨R₃⟩/∂i | NOT DERIVED | Both partial derivatives needed |
| Secular rate dM/dt from ∂⟨R₃⟩/∂a | NOT DERIVED | Mean motion correction |
| Secular rate de/dt from ∂⟨R₃⟩/∂ω | NOT DERIVED | Eccentricity secular change |
| Secular rate di/dt from ∂⟨R₃⟩/∂Ω | NOT DERIVED | Inclination secular change |
| Moon/Sun scaling: Moon 2.2× stronger | CITED | Ratio of (GM/r³) for Moon vs Sun |

## Source Documents Required

| Source | Location | Availability | What it provides |
|---|---|---|---|
| Self-contained orbital mechanics math | N/A | ✓ | Legendre expansion, orbit averaging, Lagrange equations |
| Sneeuw (2022) "Dynamic Satellite Geodesy" | `sgp4_references/.../Sneeuw_2022_...LectureNotes.pdf` | ✓ IN REPO | Eqs. 6.13-6.17: Kaula expansion with F_lmk, G_lkq; Eqs. 7.5a-f: LPE with Kaula functions; Eqs. 7.1a-c: partial derivatives of R |
| Na et al. (2012) "Detailed Re-derivation..." | `sgp4_references/.../Detailed_Re-derivation_...Kaul.pdf` | ✓ IN REPO | Eqs. 3-13a-f: Lagrange planetary equations (clean digital); NOTE: corrects one error in Kaula (1966) — must verify |
| Lara (2021) "Brouwer's satellite solution redux" | `sgp4_references/.../Lara_2021_...arXiv_2009.10665v2.pdf` | ✓ IN REPO | Eq. 14: Second-order secular Hamiltonian (cross-verified) |
