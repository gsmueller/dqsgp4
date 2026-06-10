# Third-body gravitational perturbation — the Cartesian force (Sun / Moon)

**Theory note (FIRST, before code) for roadmap step L4 / survey item (a)** — wire a third-body lunisolar
perturbation that consumes the L3 ephemeris. Governed by [[feedback_theory_first_library]] (theory precedes code)
and [[feedback_no_perceived_fidelity]] (every value gated against an independent oracle). This note derives the
Newtonian third-body perturbing acceleration, its numerically-stable evaluation (Battin's `f(q)`), the frame and
unit reconciliation, and the two-layer JPL DE430 oracle. The companion code is `src/forces/third_body.h`; the gate
is `tests/test_third_body`.

---

## 1. Scope — why a *Cartesian force*, not the secular-rate orphan

The library already holds a third-body routine: `perturbation/third_body.h::compute_third_body_rates`. It is the
**orbit-averaged, P₂(cos S) secular-rate** form — it returns dΩ/dt, dω/dt, dM/dt, de/dt, di/dt via Lagrange's
planetary equations, the same averaged-element approach SR3 uses inline (DPPER/DPINIT). It is **unwired** (no
`src/` consumer). It belongs to an *averaged-element* propagator, and wiring it would brush the frozen SR3 path.

The DQ propagator does **not** integrate averaged elements — it integrates the **Cartesian Newtonian state**
(`gravity_central`, `gravity_J2`, `drag` are all Cartesian `Force`s returning a body-frame `Wrench`). So the
member that fits — and the one the architecture survey called *"the single biggest latent integration: everything
exists, nothing is wired"* — is a **Cartesian third-body acceleration `Force`**, built from the already-shipped,
DE430-gated `ephemeris/body_position_gcrs.h` + textbook Newtonian physics. This **never touches SR3**: it is a new
additive force, exactly like `geopotential`, gated against an independent oracle. The secular-rate orphan stays as
it is (it serves a different, averaged consumer — anti-dead-code, not this increment's concern).

---

## 2. The Newtonian third-body perturbing acceleration

Let (Earth-centred inertial):

```
  r   = satellite geocentric position
  s   = third-body (Sun or Moon) geocentric position
  μ₃  = GM of the third body
```

The satellite and the Earth are **both** attracted by the third body. In the Earth-centred (non-inertial) frame
the satellite's equation of motion carries the third body's *direct* pull on the satellite minus its pull on the
Earth (the indirect/indirect-oblateness term):

```
  a_3body  =  μ₃ [ (s − r) / ‖s − r‖³   −   s / ‖s‖³ ]
                 └─ direct: body pulls satellite ─┘   └─ indirect: body pulls Earth ─┘
```

The bracketed difference is the **tidal** perturbation: a uniform pull on the whole Earth–satellite system
produces no relative motion; only the *gradient* of the body's field across the system perturbs the orbit. The
magnitude scales as the tidal field `μ₃‖r‖/‖s‖³` — the Moon (μ₃ small, s small) beats the Sun (μ₃ huge, s huge)
by ~2.2× because the `1/s³` tidal falloff rewards proximity over mass.

---

## 3. Numerical stability — the `f(q)` reformulation (the heart of the note)

Equation §2 evaluated literally is a **catastrophic-cancellation trap**: for a LEO satellite the two vectors
`(s−r)/‖s−r‖³` and `s/‖s‖³` agree to ~6–8 significant figures (‖r‖/‖s‖ ~ 7000 km / 1.5e8 km ~ 5e-5 for the Sun),
so their difference loses most of the mantissa to round-off. A perturbation that is itself ~1e-7 of central
gravity cannot be computed by subtracting two nearly-equal O(1) vectors. This is precisely the "perceived
fidelity" failure mode — the numbers would *look* like a third-body acceleration while being dominated by
subtraction noise. Battin's reformulation removes the cancellation analytically.

Decompose, writing `ρ = s − r`:

```
  a_3body = μ₃[ ρ/‖ρ‖³ − s/‖s‖³ ]
          = μ₃[ s/‖ρ‖³ − r/‖ρ‖³ − s/‖s‖³ ]
          = μ₃[ −r/‖ρ‖³ + s·(1/‖ρ‖³ − 1/‖s‖³) ].
```

The `−r/‖ρ‖³` term is harmless (it is O(‖r‖), no cancellation). All the cancellation lives in
`(1/‖ρ‖³ − 1/‖s‖³)`. Define the dimensionless

```
  q = (r · (r − 2s)) / ‖s‖²          (= (‖r‖² − 2 r·s)/‖s‖²,  so ‖ρ‖² = ‖s‖²(1+q))
```

Then `1/‖ρ‖³ − 1/‖s‖³ = ‖s‖⁻³[(1+q)^(−3/2) − 1] = −f(q)/‖s‖³`, with `f(q) = 1 − (1+q)^(−3/2)`. Substituting and
collecting over `‖s‖³` (using `‖ρ‖³ = ‖s‖³(1+q)^(3/2)`):

```
  ┌─────────────────────────────────────────────────────────────────────┐
  │  a_3body = −(μ₃ / ‖s‖³) · [ r·(1+q)^(−3/2)  +  f(q)·s ]              │
  └─────────────────────────────────────────────────────────────────────┘
```

with `f(q)` evaluated in the **cancellation-free** form (rationalising `(1+q)^(3/2) − 1` via the difference of
cubes of `(1+q)^(1/2)`):

```
  f(q) = q (3 + 3q + q²) / [ (1+q)^(3/2) · ((1+q)^(3/2) + 1) ]
```

Now both pieces are well-conditioned: `f(q) ≈ 3q/2` for small `q`, and `(1+q)^(−3/2) ≈ 1`, so the two terms of
the bracket are each O(‖r‖) and *add* (no subtraction of near-equal quantities). The perturbation is recovered to
full machine precision.

**Sign / magnitude check** (satellite on the Earth–body line, `r = ε s`, `ε = ‖r‖/‖s‖ ≪ 1`): `q = ε(ε−2) ≈ −2ε`,
`f(q) ≈ 3q/2 ≈ −3ε`, so `a ≈ −(μ₃/‖s‖³)[εs − 3ε s] = +2ε μ₃ s/‖s‖³`. The acceleration points **toward** the body
with magnitude `2 μ₃‖r‖/‖s‖³` — the expected tidal stretch along the Earth–body line. ✔ (A satellite closer to
the body than Earth is feels a stronger pull toward it; this is the leading tidal term.)

---

## 4. Frame and units

- **Body position.** `ephemeris/body_position_gcrs(EclipticState, Epoch)` returns the body's GCRS Cartesian
  position from the L3 `sun_meeus` / `moon_meeus` ecliptic state through the L2 frame chain (`ecliptic_to_gcrs =
  Pᵀ·Rx(−ε_A)`), DE430-gated end-to-end in `test_frame_chain` (≤3.7″ Moon, ≤25″ Sun). `radius` units: **AU** for
  the Sun, **km** for the Moon (pinned per instance in `ecliptic_state.h`).
- **Units → metres.** The DQ propagator works in SI (position m, μ in m³/s²). The third-body position is converted
  to metres: ×AU (Sun) / ×1000 (Moon). AU is the IAU 2012 *defining* constant (`defined`, exact). μ₃ are
  born-digital (§5).
- **Frame consistency (gate vs wiring).** The acceleration in §3 requires `r` and `s` in the **same** frame. The
  standalone force + gate use a single consistent frame (GCRS). Wiring into the DQ propagator — whose state is
  **TEME** (SGP4-seeded, TEME→ECEF = Rz(GMST)) — additionally needs the TEME↔GCRS rotation (precession + nutation,
  the L2 chain); that frame step rides with the *wiring*, which is deferred (§7). For a perturbation at the 1e-7
  level the residual frame misalignment is itself an accuracy contribution, bounded there.

---

## 5. The constants (born-digital, honest)

Three new physical constants, each with provenance per the Constants Initiative
([[feedback_constants_generative_or_bounded]]):

| constant | value | provenance | error budget |
|---|---|---|---|
| AU | 149 597 870 700 m | `defined` (IAU 2012 Res. B2 — exact by definition) | exact; no accuracy error |
| GM☉ | 1.327 124 400 18 × 10²⁰ m³/s² | `measured`/`model_coefficient` (IAU 2015 / DE430 heliocentric) | adopted; last-digit truncation |
| GM☾ | 4.902 800 066 × 10¹² m³/s² | `measured`/`model_coefficient` (DE430 lunar, = GM⊕/81.30057) | adopted; last-digit truncation |

GM☉/GM☾ are *adopted* values (finite published digits) → `model_coefficient` (digits → accuracy, binary storage
→ T-precision), exactly as `sun_meeus`'s VSOP coefficients are tagged. AU is `defined` (exact since 2012). These
sit with the force/`ThirdBody` descriptor as born-digital `model_coefficient`/`defined` values (W18-clean: no
float literals); promoting them into `ConstantsProvider` with a CR1B-allowlist entry is deferred until a provider
consumer needs them (anti-dead-code).

---

## 6. The abstraction — the functions fall out

```cpp
// A third body = its GM + a position-of-epoch in GCRS metres (the ephemeris instance, unit-converted).
template<class T> struct ThirdBody {
    TrackedValue<T> mu;                                   // μ₃  [m³/s²], born-digital
    std::function<Vector3<T>(const Epoch<T>&)> position;  // GCRS metres at an epoch (sun_meeus / moon_meeus)
};

// Newtonian perturbing acceleration in the satellite's (GCRS) frame — the §3 stable form.
template<class T> Vector3<T> third_body_accel(const Vector3<T>& r_sat, const Vector3<T>& s, const TrackedValue<T>& mu);

// The propagator Force: body-frame wrench. Captures the body descriptor + base epoch; maps state.time → epoch,
// evaluates the ephemeris, and hands the world accel to the shared epilogue.
template<class T> Wrench<T> third_body_force(const State<T>&, const ConstantsProvider<T>&, const ThirdBody<T>&, const Epoch<T>& base);
```

`third_body_accel` is the §3 formula (pure, the unit of the machine-eps oracle). `third_body_force` composes it
with `body_position_gcrs` and closes with `force_wrench_from_world_accel` (the shared 4-force epilogue,
`force_common.h`). Generic over the body: Sun and Moon are two `ThirdBody` instances of one force — the same
body-genericity `sun_meeus`/`moon_meeus` already realise for the ephemeris.

---

## 7. Validation — the two-layer DE430 oracle (no perceived fidelity)

The fidelity claim factors into two independent checks, mirroring the frame-chain gate:

1. **The formula — machine-eps, independent re-derivation.** `third_body_accel(r, s, μ)` (the §3 Battin form) vs a
   naive direct evaluation of §2 `μ[(s−r)/‖s−r‖³ − s/‖s‖³]` at *high precision* (`cpp_bin_float_50`, where the
   cancellation is harmless), over LEO/GEO satellites × Sun/Moon geometries. They must agree to ≪ the perturbation
   magnitude — this proves the stable reformulation is *algebraically* the perturbation, and (run in `double`)
   exhibits the cancellation the high-precision path avoids. The naive form is the independent oracle for the
   *formula*; it is never used in production (it is the thing §3 exists to replace).
2. **The ephemeris→acceleration map — vs JPL DE430.** Feed the force the **model** body position
   (`sun_meeus`/`moon_meeus` → `body_position_gcrs`) and compare to the acceleration from the **DE430** body
   position (the independent *numerical* ephemeris, NOT the analytical theory the model realises — the
   independence lesson). The residual is the body-position truncation (Meeus ≤25″ Sun / ≤3.7″ Moon) propagated
   through the acceleration, which the model's tracked accuracy must **majorize**. `tools/gen_third_body_oracle.py`
   emits the DE430 reference; `test_third_body` embeds it.

Physical sanity also checked: |a☾| ≈ 2.2 |a☉| at a fixed satellite; |a_3body| ~ 1e-6 m/s² (≈1e-7 of central
gravity at LEO); the acceleration points along the Earth–body line for an on-line satellite (§3 sign check).

**Frozen invariants.** This is a new additive force; OR1 (`test_sgp4_regression`) and `test_sgp4` (33/33) are
analytic SGP4 and untouched. The force is **not** wired into the DQ propagator's default list in this increment —
the **gate is its consumer** (exactly as `geopotential` was gated standalone before its later DQ switch); wiring
it (with the TEME↔GCRS frame step of §4) is a separate, gated increment once a propagation consumer wants it.

---

## 8. Summary

```
  q          = (r·(r − 2s))/‖s‖²
  f(q)       = q(3+3q+q²) / [ (1+q)^{3/2}((1+q)^{3/2}+1) ]          (cancellation-free)
  a_3body    = −(μ₃/‖s‖³)[ r(1+q)^{−3/2} + f(q) s ]                  (Battin, machine-precise)
  ThirdBody  = {μ₃ born-digital, position = body_position_gcrs(sun/moon_meeus)·unit}
  Force      = third_body_accel → force_wrench_from_world_accel       (Sun, Moon = two instances)
  oracle     = (1) naive §2 form @ bf50 [formula]  ⊕  (2) DE430 body position [ephemeris→accel]
```

A new Cartesian perturbation force, the second real consumer of `body_position_gcrs`, resolving the survey's
biggest latent integration — built from textbook physics + the gated ephemeris, with SR3 untouched and every
value DE430-gated.

---

## 9. Propagator wiring — the State→Wrench force (the precession-frame step)

The deferred wiring of §7, landed as an OPT-IN force (`third_body_force` / `make_third_body_force`), not in the
default DqSgp4 list — so the DQ propagator's existing value path is unchanged and the wiring's gate is its consumer.
Three pieces:

1. **Epoch map.** The DQ propagator's `State::time` is elapsed seconds since the propagation epoch. The absolute
   epoch is `base ⊕ time` formed two-part — `Epoch::from_jd_two_part(base.jd_day, base.jd_frac + time/86400)` —
   preserving the day part (L1 §3). The ephemeris is evaluated there.

2. **The frame step — GCRS → TEME by precession alone (the honest approximation).** The DQ propagator's inertial
   frame is **TEME** (SGP4-seeded; TEME→ECEF = Rz(GMST) in the geopotential force). `body_position_gcrs` returns
   GCRS. The exact map is `GCRS→TEME = R3(−Eqeq)·N·P` (precession, then nutation N, then the equinox offset Eqeq).
   Its terms separate cleanly by size: **precession P ~ 0.3°** at 2024 (the dominant, accumulating ~1.4°/century),
   then **nutation N ~ 17″** and **the equation of equinoxes Eqeq ~ 17″** (both arcsec). We apply the dominant,
   erfa-gated `precession_iau2006` (`s_teme = P·s_gcrs`; FRAME2 gates P vs erfa.pmat06 at machine-ε) and **omit
   N + Eqeq**. `third_body_accel` is rotation-covariant, so rotating only the body position is equivalent to
   rotating the whole problem — the third-body physics stays in its TB1-gated form.

3. **The omitted-frame bound (no perceived fidelity).** Omitting N + Eqeq is a ~√(17²+17²) ≈ 24″ direction error
   on the body, hence a ~24″ relative error on the acceleration. It is **declared, not hidden**: a conservative 30″
   bound is deposited into the acceleration's `accuracy` channel (`add_bound`). The end-to-end fidelity is then:
   **Sun** — ephemeris-limited (25″ Meeus ≳ 24″ frame), so precession-only loses nothing; **Moon** — *frame-limited*
   at ~24″ (its ephemeris is ~4–10″), an honest, tracked degradation. Tightening the Moon below 24″ is the concrete
   consumer that would justify building **nutation N(t)** + Eqeq (the IAU 2000B series, erfa-gated) — deferred until
   that need is real (anti-dead-code; the consumer now exists in principle but is not yet forcing). For the DQ
   propagator dynamics the third-body force is ~1e-7 of central gravity, so a 24″ frame error perturbs the
   integrated state by ~1e-7·1e-4 — far below the round-off-tolerant propagator gates.
