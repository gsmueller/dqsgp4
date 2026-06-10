# The Geopotential Acceleration — One Field, One Recursion

*Theory note preceding `constants/gravity_field.h` + `forces/geopotential.h` (corrected
formula-layer Stage 2). Per the theory-first mandate, the acceleration for every order —
monopole, zonal (m=0), tesseral (m≥1) — is derived here so the evaluator FALLS OUT of a
single gradient, and per the no-perceived-fidelity mandate the independent oracle is fixed
before code. Sources, all born-digital or in-repo: Cunningham (1970) "On the computation of
the spherical harmonic terms…", Cel. Mech. 2; Montenbruck & Gill, *Satellite Orbits* §3.2
(eqns 3.29–3.33); Heiskanen & Moritz, *Physical Geodesy* §1; Kaula, *Theory of Satellite
Geodesy* §1.3.*

## 1. The potential

The exterior gravitational potential of a body, in body-fixed spherical coordinates
(r, φ = geocentric latitude, λ = longitude), is the solid-harmonic series

    U(r,φ,λ) = (GM/r) Σ_{n=0}^∞ Σ_{m=0}^n (R/r)^n P_nm(sin φ)
                                  · [ C_nm cos(mλ) + S_nm sin(mλ) ],          (1)

with R the reference radius, P_nm the (unnormalized) associated Legendre functions, and
C_nm, S_nm the dimensionless Stokes coefficients. The decomposition by order m is exactly
the zonal / tesseral / sectoral split:

  - **monopole** n=m=0: C_00 = 1, and the n=0 term is (GM/r) — the point-mass potential;
  - **zonal** m=0, n≥2: longitude-independent (no λ dependence), C_n0 = −J_n, S_n0 = 0.
    The minus sign is the J_n convention (Kaula §1.3): U_zonal = −(GM/r) Σ J_n (R/r)^n P_n;
  - **tesseral/sectoral** m≥1: the longitude-dependent C_nm, S_nm.

There is exactly ONE coefficient table here; "zonal harmonics" and "tesseral harmonics" are
the m=0 column and the m≥1 block of it. This is what `GravityField` makes literal:

    C(0,0) = 1                         (monopole)
    C(n,0) = −J_n,  S(n,0) = 0         (zonal column; J_n from the ellipsoid/EGM generators)
    C(n,m) = C_nm,  S(n,m) = S_nm      (tesseral block, m≥1, denormalized Stokes)

No coefficient string or honesty tag changes: `GravityField` simply *presents* the existing
`ZonalHarmonics` Jₙ (as −Jₙ at m=0) and `TesseralHarmonics` C_nm/S_nm behind one accessor.
The Kaula denormalization N_nm·C̄_nm at m=0 reduces to N_n0 = √(2n+1), recovering exactly the
`jn_from_cbar` relation Jₙ = −√(2n+1)·C̄ₙ₀ already in `zonal_harmonics.h` — the two halves are
one table, confirming the fusion adds no new numerics.

## 2. Cunningham's singularity-free harmonics

Writing the gradient of (1) directly in (r,φ,λ) introduces 1/cos φ poles. Cunningham (1970)
removes them by working in body-fixed Cartesian (x,y,z) with the two real solid harmonics

    V_nm = (R/r)^{n+1} P_nm(sin φ) cos(mλ),   W_nm = (R/r)^{n+1} P_nm(sin φ) sin(mλ).   (2)

These satisfy the pure recursions (no Earth constants) already implemented in
`math/spherical_harmonics.h` (`cunningham_vw`) — Montenbruck & Gill (3.29)–(3.31):

    sectoral   V_mm = (2m−1)[ (xR/r²) V_{m-1,m-1} − (yR/r²) W_{m-1,m-1} ],  (W likewise),
    vertical   V_nm = ((2n−1)/(n−m))(zR/r²) V_{n-1,m} − ((n+m−1)/(n−m))(R²/r²) V_{n-2,m},

seeded by V_00 = R/r, W_00 = 0. Note V_n0 = (R/r)^{n+1} P_n(sin φ): **the m=0 column of V is
exactly the zonal harmonic** — the same physical term `gravity_zonal` builds from a
Legendre-in-u recurrence. Same physics, different arithmetic (this is the genuine
duplication the unification removes from the runtime path; §5).

## 3. The acceleration — one formula, three cases

The acceleration is a = ∇U. Cunningham's key identity is that ∂/∂{x,y,z} of a degree-n
solid harmonic is a linear combination of degree-(n+1) harmonics, so the gradient of (1)
re-sums into V/W of degree n+1. With the scale κ ≡ GM/R², Montenbruck & Gill (3.33):

**m ≥ 1** (tesseral; the existing `gravity_tesseral` sum):

    ẍ_nm = κ·½[ (−C_nm V_{n+1,m+1} − S_nm W_{n+1,m+1})
                 + (n−m+2)(n−m+1)( C_nm V_{n+1,m-1} + S_nm W_{n+1,m-1}) ]
    ÿ_nm = κ·½[ (−C_nm W_{n+1,m+1} + S_nm V_{n+1,m+1})
                 + (n−m+2)(n−m+1)(−C_nm W_{n+1,m-1} + S_nm V_{n+1,m-1}) ]

**m = 0** (zonal + monopole). The ½-formula cannot be used verbatim: it needs the
order m−1 = −1 harmonics. Negative orders fold back via the reflection (Cunningham 1970)

    V_{n,-1} = −((n−1)!/(n+1)!) V_{n,1},   W_{n,-1} = −((n−1)!/(n+1)!) W_{n,1},

i.e. V_{n,-1} = −V_{n,1}/[(n+1)n], W likewise. Substituting m=0 into the ½-formula, the
(n+2)(n+1) factor multiplies the −1 terms = −V_{n+1,1}/[(n+2)(n+1)], and the two halves
COMBINE to a single clean term (factor 1, not ½):

    ẍ_n0 = κ·( −C_n0 V_{n+1,1} ),     ÿ_n0 = κ·( −C_n0 W_{n+1,1} ).            (3)

**all m** (the z-component is uniform):

    z̈_nm = κ·(n−m+1)( −C_nm V_{n+1,m} − S_nm W_{n+1,m} ).

**The monopole falls out.** Put n=m=0, C_00=1 in (3) and the z-row. With V_11 = xR²/r³,
W_11 = yR²/r³ (from the seed V_00=R/r through the sectoral step), V_10 = zR²/r³:

    ẍ = κ(−V_11) = (GM/R²)(−xR²/r³) = −GM x/r³,   and likewise ÿ, z̈,

which is **exactly point-mass gravity −GM r/r³**. So a single loop over 0 ≤ m ≤ n ≤ N with
the (3)/½/z-row cases reproduces monopole + zonal + tesseral — no separate central or zonal
code path. This is the evaluator `geopotential_accel_ecef`; `gravity_J2` becomes the (2,0)
term and `max_m=0` selects zonal-only.

## 4. Truncation — the Kaula tail (accuracy channel)

Summing to degree N omits n>N. Kaula's rule of thumb bounds the unnormalized coefficient
magnitude by |C_nm|,|S_nm| ≲ 10⁻⁵/n², and (R/r)^{n} decays geometrically for r>R, so the
omitted acceleration is bounded by a convergent tail

    ‖a − a_N‖ ≤ κ · Σ_{n>N} (n+1)(R/r)^{n+1} · c_n,   c_n = Kaula bound at degree n,

deposited in the **accuracy** channel via `add_bound(..., ErrorChannel::accuracy)` — the
model-fidelity floor of a finite field, T-independent. (The current `gravity_tesseral`
lacks this tracked tail; the unified evaluator gives both halves one honest bound.) The
**precision** channel keeps coming, as always, from the TrackedValue arithmetic of the V/W
recursion and the coefficient generators (tightens with wider T).

## 5. Independent oracle — round-off, NOT bit (no perceived fidelity)

The unified evaluator changes the *arithmetic* of the zonal term: Cunningham V/W (column
m=0) vs `gravity_zonal`'s Legendre-in-u recurrence. These are algebraically equal but differ
at the last ULPs — so the claim is **round-off agreement, never bit-exact** (selling it as
bit-exact would be a seam). The gate `test_geopotential` pins three independent checks:

  1. **Summed-legacy oracle.** geopotential(field, N, M) ≈ gravity_central + gravity_zonal(N)
     + gravity_tesseral(N) at a stated round-off tolerance (≈1e-12 relative), over a spread
     of sample positions (LEO/MEO/GEO, high/low latitude). The legacy trio is independently
     gated (W8/D2/GAL1) and uses *different* arithmetic for the zonal/central halves — a
     genuine independent reference, not the same upstream.
  2. **Closed-form J₂.** With only (2,0) active, the evaluator must match the textbook
     closed J₂ acceleration a_J2 = −(3/2)J₂(GM/r²)(R/r)²·[ (1−5(z/r)²)x̂/r, …, (3−5(z/r)²)ẑ/r ]
     (Montenbruck & Gill 3.34) to round-off — an oracle external to the recursion itself.
  3. **Monopole identity.** With only (0,0), the evaluator must equal −GM r/r³ to round-off.

Only after all three pass does the DQ-side force list switch to the evaluator. The frozen
SGP4 inline J₂/J₃/J₄ (`near_space.h`/`deep_space.h`) and `sgp4::ZonalHarmonics` are NOT
routed through it — they stay bit-frozen (OR1, 33/33).

## 6. What falls out

  - `constants/gravity_field.h` — `GravityField<T>`: one C(n,m)/S(n,m)/J accessor fusing the
    existing `ZonalHarmonics` (m=0 column, as −Jₙ) and `TesseralHarmonics` (m≥1), honesty
    tags threaded VERBATIM; `egm2008()/wgs72()` factories delegate to the existing two.
  - `forces/geopotential.h` — `geopotential_accel_ecef(r_ecef, mu, Re, field, max_n, max_m)`:
    one `cunningham_vw` pass + the §3 loop; monopole+zonal+tesseral; Kaula tail in accuracy;
    plus the shared `force_wrench_from_world_accel(state, a_world)` body-frame epilogue (the
    4-copy rotation-to-body that ends every force lambda).
  - `tests/test_geopotential` — the §5 round-off gate (new ExeGate; takes 69→70 gates).

The evaluator subsumes `gravity_central`/`gravity_zonal`/`gravity_tesseral`/`gravity_J2` on
the runtime path; the legacy functions remain as the gate's independent oracle.
