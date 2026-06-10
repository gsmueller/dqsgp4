# Reference Frames — Theory Note (L2)

**Module L2 of the professional-library re-architecture** (`design/PROFESSIONAL_LIBRARY_PLAN.md`). THEORY
FIRST: this note develops the frame theory exhaustively; the `Frame` type and the rotation transforms (§5)
are its direct consequence. No frame code is written or changed before this note stands on its own
([[feedback_theory_first_library]]). Fidelity is verified against an independent oracle —
ERFA (IAU SOFA) / astropy plus analytic identities — never asserted ([[feedback_no_perceived_fidelity]]).
Builds on L1 (`Epoch`, `centuries_since_j2000` — [[time_scales_and_epochs]]).

## 0. Why this module exists (grounding: the existing state)

The library has exactly **one** frame transform, and it is inline and incomplete:
- `gravity_tesseral.h` (D2) rotates the satellite position ECI→ECEF by `Rz(−GMST)` about +z and back by
  `Rz(+GMST)`, hand-coded componentwise (`cg = cos(gmst)`, `sg = sin(gmst)`; `x' = cg·x + sg·y`, …). GMST
  comes from L1 (`sidereal_time.h`).
- The SGP4 output is **TEME** (True Equator, Mean Equinox of date), and the TEME→PEF(≈ECEF) rotation is taken
  as `Rz(GMST)` — the standard SGP4 convention, which omits the equation of equinoxes and polar motion.
- The mean obliquity ε(t) exists as a generator (`obliquity.h`, SC1, oracle-verified) but is not wired into a
  frame rotation. `Matrix3<T>` exists (`matrix3.h`, H1) but the one rotation above is hand-coded, not a matrix.
- There is **no precession, no nutation, no `Frame` abstraction**. The lunar/solar ephemeris (L3) emits
  mean-ecliptic-of-date λ/β with no typed path into the equatorial propagation frame.

So frame handling is a single ad-hoc rotation with the frame identities living only in comments — the same
"missing theory" smell L1 cured for epochs. Writing the rotation theory down yields a reusable, composable,
oracle-verified `Frame` layer; the inline `Rz(GMST)` becomes one named transform within it.

## 1. The frames and the rotation that defines each (conventions pinned)

A *reference frame* is an orthonormal triad; relating two is a proper rotation (∈ SO(3), det = +1). The chain
from the celestial frame to the Earth-fixed frame, each step a named rotation evaluated at an `Epoch`:

| Frame | Definition | Rotation from previous |
|---|---|---|
| **GCRF / ICRF** (≈ J2000) | quasi-inertial celestial frame; mean equator & equinox of J2000 | — (root) |
| **MOD** (mean of date) | mean equator & equinox of date | **precession** `P(t)` |
| **TOD** (true of date) | true equator & equinox of date | **nutation** `N(t)` |
| **TEME** | true equator, *mean* equinox of date (SGP4) | `Rz(−Eq.Equinoxes)` from TOD (small) |
| **Ecliptic of date** | mean ecliptic & equinox of date | `Rx(+ε(t))` from MOD (mean obliquity, SC1) |
| **PEF / ECEF (ITRF)** | Earth-fixed | sidereal `Rz(GMST/GAST)` then polar motion `W(t)` |

The full GCRF→ECEF map is `W·R·N·P`. **SGP4 works in TEME and uses TEME→ECEF ≈ `Rz(GMST)`** (no
equation-of-equinoxes, no polar motion) — the existing `gravity_tesseral.h` transform.

## 2. The rotation algebra (primitives)

Elementary rotations as `Matrix3<T>` (right-handed, active about a body axis):

```
        ⎡1    0      0  ⎤        ⎡ cosθ 0 −sinθ⎤        ⎡ cosθ  sinθ 0⎤
 Rx(θ)= ⎢0  cosθ  sinθ ⎥  Ry(θ)=⎢ 0    1   0  ⎥  Rz(θ)=⎢−sinθ  cosθ 0⎥
        ⎣0 −sinθ  cosθ ⎦        ⎣ sinθ 0  cosθ⎦        ⎣  0     0   1⎦
```

Each is orthonormal (`RᵀR = I`), `det = +1`, `R(θ)⁻¹ = R(−θ) = R(θ)ᵀ`. Frame chains are products; the
inverse of a chain is the reversed product of transposes. (The `gravity_tesseral.h` componentwise `Rz` is
exactly `Rz(−gmst)` applied to a vector — it becomes `sidereal_rotation(gmst)` below, value-preserving.)

## 3. The component transforms (theory + provenance)

- **Precession `P(t)`** — IAU 2006 (Capitaine) equatorial precession, the angles (ζ_A, θ_A, z_A) or the
  4-rotation (γ, φ, ψ, ε) formulation, **polynomials in `T = centuries_since_j2000` (L1)**. Each coefficient
  a finite-digit model fit → `model_coefficient` (digits → accuracy, storage → T-precision); the series order
  is a truncation → accuracy. Born-digital (IERS Conv. 2010 / IAU 2006).
- **Nutation `N(t)`** — IAU 1980 or 2000 series: `Δψ` (longitude), `Δε` (obliquity), each a sum of periodic
  terms in the five Delaunay arguments (l, l′, F, D, Ω). The amplitudes are born-digital; keeping the largest
  K terms is a truncation whose bound (Σ omitted amplitudes) is the **accuracy** — the same generator pattern
  as the lunar periodics. (For SGP4-grade work, nutation is often dropped; then `N = I` with a tracked ≈ arcsec
  accuracy.)
- **Obliquity `ε(t)`** — already a generator (`obliquity_iau2006`, SC1), oracle-verified; the rotation
  equatorial↔ecliptic is `Rx(±ε)`.
- **Sidereal rotation `R`** — `Rz(GMST)` (mean) or `Rz(GAST)` (true; GAST = GMST + Δψ·cos ε, the equation of
  equinoxes). GMST is L1's `sidereal_time`. The SGP4 TEME→ECEF uses the GMST form.
- **Polar motion `W(t)`** — `Rx(−yp)·Ry(−xp)`, with (xp, yp) IERS-observed (sub-arcsec); a `measured` input,
  default identity for SGP4-grade with a tracked bound.

## 4. The three-error budget for a transform

- Rotation primitives: exact functions of the (tracked) angle → carry the angle's precision; orthonormality is
  exact (identity-verified).
- Precession/nutation: **accuracy** = series truncation (omitted terms) + the coefficient digit floors;
  **precision** scales with `T`. The dominant *modelling* choice — dropping nutation, or treating TEME as TOD
  — is carried as an explicit accuracy (~arcsec), not hidden.
- TEME convention: the equation-of-equinoxes between TEME and TOD (~arcsec) is an accuracy if a TEME vector is
  consumed as TOD; zero when the SGP4 `Rz(GMST)` convention is used end-to-end (value-preserving).

## 5. The abstraction that falls out (the functions, as consequences)

```cpp
enum class Frame { GCRF, MOD, TOD, TEME, EclipticOfDate, ECEF };

// Primitives — Matrix3<T> (H1), exact orthonormal rotations:
template<typename T> math::Matrix3<T> rot_x(const math::TrackedValue<T>& a);
template<typename T> math::Matrix3<T> rot_y(const math::TrackedValue<T>& a);
template<typename T> math::Matrix3<T> rot_z(const math::TrackedValue<T>& a);

// Component transforms — Epoch-parameterised (L1), provenance-tagged (§3):
template<typename T> math::Matrix3<T> precession (const Epoch<T>&);          // P(t)
template<typename T> math::Matrix3<T> nutation  (const Epoch<T>&, int n_terms); // N(t), truncation→accuracy
template<typename T> math::Matrix3<T> obliquity_rotation(const Epoch<T>&);   // Rx(ε), uses SC1
template<typename T> math::Matrix3<T> sidereal_rotation (const math::TrackedValue<T>& gmst); // Rz(GMST)

// Composed: the rotation that takes a vector from `to` into `from` coordinates.
template<typename T> math::Matrix3<T> frame_transform(Frame from, Frame to, const Epoch<T>&);
```

All `Matrix3<T>` (composable, invertible via transpose), `Epoch`-parameterised, `TrackedValue`-tracked. The
existing `gravity_tesseral.h` rotation is `sidereal_rotation(gmst)` applied to the position — the inline `Rz`
disappears into the named, tested primitive.

## 6. Oracle & verification (independent truth, never self-consistency)

| Claim | Independent oracle / check |
|---|---|
| `Rx/Ry/Rz` orthonormal, det = +1, `R(θ)R(−θ) = I`, `R(0) = I`, axis images at π/2 | analytic identities |
| frame chain round-trip `GCRF→ECEF→GCRF = I` | identity (composition of orthonormal) |
| precession / nutation / obliquity / GMST matrices | **ERFA (IAU SOFA) / astropy** to a stated tolerance |
| `sidereal_rotation(gmst)` vs the inline `gravity_tesseral.h` `Rz(±gmst)` | **bit-exact** (value-preserving migration) |
| nutation truncation bound | Σ omitted amplitudes vs the full ERFA series |
| SGP4 path | **OR1 0 km** + 33/33 — TEME/`Rz(GMST)` convention reproduced exactly |

ERFA is the IAU-standard implementation, independent of any in-repo code, so a match is proof. The TEME→ECEF
`Rz(GMST)` convention reproduces today's behavior bit-for-bit.

## 7. Migration (value-preserving; SGP4 path stays OR1-frozen)

`gravity_tesseral.h`'s inline `Rz(±gmst)` → `sidereal_rotation(gmst)` (a `Matrix3` applied to the vector),
**bit-exact** (same cos/sin componentwise). Precession/nutation/obliquity rotations are **new** capability
(consumed by L3's ephemeris frame-consistency and any GCRF↔ECEF need), each oracle-gated; they do not touch
the SGP4 path. Implementation proceeds only after this note is accepted, bottom-up, each unit isolated and
oracle-gated (a new `test_reference_frames` gate).

## References (born-digital)

- IERS Conventions (2010), Ch. 5 — precession-nutation, the rotation sequence, polar motion.
- IAU 2006 / Capitaine et al. — precession; IAU 1980/2000 nutation series.
- Vallado, *Fundamentals of Astrodynamics*, Ch. 3 — TEME, the GCRF→ECEF chain, the SGP4 convention.
- ERFA (IAU SOFA) / astropy `coordinates` — the verification oracle.
