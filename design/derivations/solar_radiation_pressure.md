# Solar radiation pressure — the cannonball force, generatively sourced

**Theory note (FIRST, before code) for replan item R2** (`design/PROFESSIONAL_LIBRARY_PLAN.md` §6) — the force
that completes the L4 layer's classical set (geopotential ✓, drag ✓ R1, third-body ✓, SRP here). Governed by
[[feedback_theory_first_library]], [[feedback_no_perceived_fidelity]], and the generative-constants mandate
([[feedback_constants_generative_or_bounded]]). Companion code `src/forces/srp.h`; gate `tests/test_srp` (SRP1).

---

## 1. The physics — momentum flux of sunlight

A photon of energy E carries momentum E/c. A surface intercepting radiative flux S [W/m²] therefore absorbs
momentum at rate S/c per unit area — a **pressure**

```
  P(d) = S(d)/c ,        S(d) = L☉ / (4π d²)        (isotropic source at distance d)
```

For the **cannonball model** (the standard first-order SRP treatment: an effective area A, mass m, and a
radiation-pressure coefficient C_R that absorbs the surface's reflectance — C_R = 1 perfect absorber, → 2 for
a perfect normal specular reflector; typically ~1.2–1.5 for spacecraft):

```
  ┌──────────────────────────────────────────────────────────────────┐
  │  a_srp = ν · P₁ᴬᵁ · (AU/d)² · C_R · (A/m) · d̂ ,   d⃗ = r − s     │
  └──────────────────────────────────────────────────────────────────┘
```

with `r` the satellite and `s` the Sun (same frame, metres), `d̂` pointing **away from the Sun** (radiation
pushes outward), and `ν ∈ {0,1}` the shadow factor (§3). The attitude-dependence of a real spacecraft (flat
panels, thermal re-emission) is *not* modeled — that is what "cannonball" means, stated openly; C_R·(A/m) is
the caller's single lumped parameter, exactly as drag's B = C_d·A/m.

## 2. P₁ᴬᵁ — generated, not copied (the constants treatment)

The 1-AU pressure is **generated from three exact-by-convention constants** (never the copied "4.56e-6"):

```
  P₁ᴬᵁ = L☉ / (4π · AU² · c)
```

| constant | value | provenance |
|---|---|---|
| L☉ | 3.828 × 10²⁶ W | IAU 2015 Resolution B3 **nominal** solar luminosity — exact by convention → `defined` (+ CR1B allowlist) |
| AU | 149 597 870 700 m | IAU 2012 B2, exact → `defined` (allowlisted since the third-body increment) |
| c | 299 792 458 m/s | SI defining constant, exact → `defined` (+ allowlist) |

Numerically: S(1 AU) = L☉/(4π AU²) = **1361.2 W/m²** (the modern total solar irradiance the B3 nominal was
chosen to round) and P₁ᴬᵁ = S/c = **4.5398 × 10⁻⁶ N/m²**. The value is thus *precision-only* through the
generator — but honesty requires the **nominal-vs-real** layer stated: the actual TSI varies ~±0.1 % over the
solar cycle about the nominal, so a 1e-3 relative **representativeness accuracy** is deposited on P₁ᴬᵁ (the
same two-layer honesty as R1's static atmosphere: fidelity to the *named convention* is exact; the convention's
distance from the variable Sun is declared, not buried).

## 3. Shadow — the cylindrical umbra (consumed subset)

Sunlight is occulted by the Earth. The **cylindrical** shadow model (the standard first-order treatment):
with ŝ the unit geocenter→Sun direction, the satellite is **sunlit** iff

```
  r·ŝ ≥ 0      (on the day side of the terminator plane)        OR
  |r − (r·ŝ)ŝ| ≥ R_E      (outside the shadow cylinder of radius R_E)
```

else ν = 0 (umbra). The refinements deliberately *not* modeled, documented with their scale: the **conical**
umbra/penumbra structure (the Sun's finite angular size ~0.53° makes the true umbra a cone shorter than the
cylinder and adds a penumbral transition; at LEO the penumbral arc lasts seconds of orbit), Earth oblateness
(~1/298 of R_E on the cylinder radius), and atmospheric refraction. These bias eclipse-boundary timing by
seconds — far below the cannonball model's own C_R uncertainty — and the boolean ν is the seam a future
illumination-fraction model replaces. Branch selection is by `.value` comparison (the established
density-regime precedent).

## 4. Frames, epoch, and the shared third-body plumbing

SRP needs exactly what the third-body force needs: the Sun's position at the state's absolute epoch, in the
propagator's TEME frame. That arithmetic — `state.time → Epoch` (two-part) then the dominant erfa-gated
IAU-2006 precession rotation, with nutation+Eqeq (≤ ~24″) omitted — is **already implemented and gated** in
`third_body_force` (`third_body_perturbation.md` §9). Two sites re-implementing it would violate the validity
criterion at birth, so it is factored into a shared helper `third_body_position_teme(body, base, elapsed)` that
both forces call (`third_body_force` refactors onto it bit-identically — same operations, same order). The
omitted-nutation direction error (~30″ bound) enters SRP as a ≤1.5e-4 relative acceleration error — deposited
in the accuracy channel like the third-body's.

The `ThirdBody<T>` descriptor (μ + position closure) is reused as the Sun handle; SRP consumes only
`.position` (μ is the gravitational member's). One descriptor, two forces — the body-genericity intended.

## 5. The abstraction — what falls out

```cpp
TrackedValue<T> solar_radiation_pressure_1au();              // P₁ᴬᵁ generated per §2 (+ TSI band)
Vector3<T>     srp_accel(r_sat, s_sun, cr_area_over_mass);   // §1 acceleration, §3 shadow, pure geometry
Wrench<T>      srp_force(state, K, sun, base, cr_am);        // the ForceFn: epoch map + TEME step + epilogue
make_srp_force(sun, base, cr_am)                             // factory for the propagator force list (opt-in)
```

`srp_accel` is the pure, oracle-checkable unit; `srp_force` composes it with the shared §4 plumbing and the
`force_wrench_from_world_accel` epilogue. Like the third-body and R1 models, it is **NOT in the default
DqSgp4 force list** — the gate is its consumer until the R3 presets wire it.

## 6. Validation (SRP1) — what is claimed, and against what

SRP adds three claims, each verified at its own layer; nothing else is asserted:

1. **The generated P₁ᴬᵁ** — S(1 AU) from the generator must land on the published TSI basis of the IAU
   nominal: |S − 1361 W/m²| < 1 W/m² (the cross-check that the generative algebra and the three defined values
   compose correctly), and P₁ᴬᵁ = S/c recomputed independently in doubles at machine-eps.
2. **The geometry** — analytic identities, no numerical reference: inverse-square (a(2d)/a(d) = 1/4 to
   round-off); direction exactly anti-sunward; shadow truth-table (sub-solar point lit, anti-solar LEO point
   dark, anti-solar point at |r⊥| > R_E lit — the cylinder boundary); ν = 0 ⇒ exactly zero wrench.
3. **Magnitude sanity** — C_R = 1.3, A/m = 0.02 m²/kg ⇒ |a| ≈ 1.18e-7 m/s² at 1 AU (the GEO-class textbook
   scale, ~1e-8 of central gravity at LEO).
4. **The budget** — the 1e-3 TSI band present on P (and majorizing the |nominal − measured-TSI| offset);
   precision tightens with wider T; the refactored third-body force stays bit-identical through the shared
   helper (its TB1 gate re-run unchanged).

The Sun position itself inherits TB1's DE430 gating — SRP adds no new ephemeris fidelity claim.
