# Architecture Reuse / Generalization Survey → Advanced API (2026-06-07)

A whole-codebase survey of how the ~69 `src/` functions can be **reused, integrated, and generalized** while
performing the *same* tasks, consolidated into a proposed advanced API. Conducted by four parallel readers over
the four architectural clusters (core-math; astronomy/ephemeris/orbit; forces/constants/geodesy;
propagation/dynamics/io), each reading every function in its cluster. This note synthesizes their findings into
the cross-cutting themes and a sequenced plan.

**Governing constraint — the frozen firewall.** The authentic SGP4/SDP4 path is behavior-locked (69 acceptance
gates + the 0-km OR1 regression + test_sgp4 33/33). Every proposal below is **additive**: a new generic API may
*wrap* a frozen function, but the frozen value path stays bit-identical. The survey's most useful structural
finding is that the frozen surface is *narrower than it looks*: only the **inline** `sgp4/` arithmetic
(`deep_space.h`, `near_space.h`, `model_*.h`), `state_from_elements`/`secular_update`/`osculating_elements`/
`modified_kepler`, `compute_gmst`/`sidereal_rotation`, the TLE/OMM parsers, and the `StateVector` accessors are
locked. The `ephemeris/*` struct family, `solar_system.h::DerivedOrbitalElements`, `kaula::inclination_function`,
and `compute_third_body_rates` are **unwired** (dead-stored / no `src/` consumer) — they can be reshaped freely.

## 1. The five cross-cutting patterns (the real reuse, not the cosmetic)

The numeric **substrate is excellent and must not be refactored**: `TrackedValue<T>` carries all three-error
logic, and `Vector3`/`Matrix3`/`Quaternion`/`DualNumber`/`DualQuaternion` compose it with *zero* error code of
their own. The duplication lives one layer up, in five patterns each re-implemented 3–5×:

### P1 — Truncated series with a tracked tail bound (the #1 win; spans math ↔ astronomy ↔ forces)
The same idea — *sum a truncated series of TrackedValues; bound the omitted tail rigorously; deposit the bound
in a chosen error channel* — is hand-rolled at least six times, with a **dead** `math/series.h::horner`:
- polynomial-in-t + Σ|cₖ||t|ᵏ tail → **accuracy**: `obliquity.h:57-78`, `earth_orbit.h:66-88` (line-for-line
  the same modulo coefficients), and the flat L0/M polynomials in `sun_meeus.h`, `sidereal_time.h`.
- Bessel harmonic sum + ρ(e) remainder → **accuracy**: `orbit/kepler_series.h:50-63` (new).
- Poisson sum + Σ omitted |Aₖ| → **accuracy**: `ephemeris/poisson_series.h` (new).
- Leibniz/geometric term-functor → **precision**: `math/series.h::alternating_series` (used only by geodesy).
- 18 hand-written `errors.{precision,accuracy} += bound` splice sites (the channel choice is load-bearing and
  currently un-auditable).

The two new files I just shipped (`kepler_series.h`, `poisson_series.h`) are *fresh instances of this very
pattern* — confirming it's the central abstraction. **The unifier is `TrackedPolynomial<T>` + a generic
`tracked_series_sum(term, tail_bound, ErrorChannel)`**, with the precision-vs-accuracy channel an *explicit*
parameter (the constants-initiative distinction must never be defaulted). Behavior-exact: identical Horner,
identical triangle-inequality bound; gated unchanged by `test_series_constants`/`test_earth_constants`.

### P2 — The ephemeris model: mean elements + equation of centre + frame (astronomy ↔ ephemeris ↔ orbit)
`x(t)=x₀+ẋ·Δt` appears 5× (`solar_ephemeris`, `lunar_ephemeris` ×3, `celestial_body`, `solar_system`); the
2-term equation of centre 3× (two now obsolete vs `orbit::equation_of_center`); the orbit radius 4× at three
fidelity grades. `FundamentalConstants`/`DerivedOrbitalElements` and `CelestialBody` are **two representations
of identical physics** (`celestial_body.h` is already `template<class FC>` over the former). These collapse to
one **`OrbitalModel<Body, Source>`** (body-generic, element-source-pluggable) whose `position_of_date` makes the
Sun the `i=0, node_rate=0` special case of the Moon — exactly the L3 re-architecture in progress (`sun_meeus.h`
+ `moon_meeus.h` are the first two instances; `poisson_series` and `equation_of_center` are their shared series
engines). The capstone is **`body_position(Epoch, Frame)`** = `position_of_date` → spherical→Cartesian →
`rot_x(obliquity_iau2006(t))` → precession/nutation → `sidereal_rotation(GMST)`, which finally gives
`obliquity.h` (currently zero consumers) and `frames.h::rot_x/y` their justifying consumer.

### P3 — The force-and-integrate loop (forces ↔ dynamics)
- The body-frame epilogue `a_body = pose.conjugate().rotate(a_world); return Wrench({}, a_body)` is verbatim in
  **4** force lambdas. → one `force_wrench_from_world_accel(state, a_world)`.
- The `r,r²,r³,1/r,sinφ` radial preamble (3×) and the `(R/r)ⁿ` ladder (2×: zonal Legendre = the m=0 column of
  the tesseral Cunningham `Vₙₘ/Wₙₘ` recursion) → **one `Geopotential<T>` that computes monopole + zonal +
  tesseral in a single Cunningham pass** (the existing `gravity_tesseral_accel_ecef` minus its `m<1` skip,
  plus the monopole). `gravity_J2` is a dead specialization of `gravity_zonal(…,2)`.
- The scalar-accuracy broadcast across pose/twist components is hand-stamped **46 times** across the
  integrators; the L1-norm of a `Twist` is re-rolled 4×; RK4 is written **twice** (`runge_kutta.h` over `State`,
  `attitude.h` over `AttitudeState`) for the same Butcher tableau. → `Vector3::add_accuracy` / `State::
  add_uniform_accuracy` / `Twist::l1_norm`, and **one `ButcherTableau<T,S>` + `rk_step` driver** (euler S=1,
  rk4 S=4, rkf78 S=13 — the rkf78 loop already *is* this driver with a hardcoded tableau).

### P4 — Coefficient providers (constants ↔ forces)
`ZonalHarmonics` and `TesseralHarmonics` are the m=0 and m≥1 halves of one table; the Kaula denormalization
`Nₙₘ·model_coefficient(C̄ₙₘ)` reduces to `−√(2n+1)` at m=0 — three near-identical sites. There are **two**
`ZonalHarmonics` types (sgp4 vs constants) encoding the same WGS72 J₃ independently. → one **`GravityField<T>`**
with `C(n,m)/S(n,m)` + `egm2008/wgs72` factories, threading each coefficient's existing honesty tag unchanged;
folded into `ConstantsProvider.gravity` so the force signature collapses back to the uniform `(State,K)→Wrench`.

### P5 — State representations & the propagator facade (dynamics ↔ io)
Three state types (`dynamics::State` DQ+twist, `AttitudeState` q+ω, `sgp4::StateVector` r+v) with no common
interface; `AttitudeState` is literally the rotational sub-object of `State`. The km↔m + identity-attitude pack
is written 3×. `dynamics::Propagator` *already is* the Force-list+Integrator+State engine; adaptive stepping
just isn't reachable through it. → `to_state`/`to_state_vector` as the one adapter module, `AttitudeState` as a
`State` view (deletes the duplicate RK4), a **`Propagatable` concept** unifying `sgp4::Propagator::propagate` /
`DqSgp4Propagator::propagate`, and `Propagator::propagate_adaptive` wiring the existing standalone RKF78.

## 2. The advanced API (the unified abstractions)

```cpp
// P1 — the truncated-series family (math/)
template<class T> struct TrackedPolynomial { /* coeffs + |c| mags */ TrackedValue<T> eval(x, n_terms, ErrorChannel); };
template<class T,class Term,class Tail> TrackedValue<T> tracked_series_sum(Term, Tail, ErrorChannel, tol, max);
// (orbit::equation_of_center [Bessel] and ephemeris::poisson_series [Delaunay] become its two realised instances)

// P2 — the ephemeris (astronomy/ephemeris/, the L3 re-architecture)
template<class T> struct OrbitalModel<Body,Source> { EclipticState<T> position_of_date(const Epoch<T>&) const; };
template<class T> Vector3<T> body_position(const OrbitalModel<T>&, const Epoch<T>&, Frame);   // composes frames+obliquity

// P3 — forces & integration
template<class F,class T> concept Force = requires(F f, State<T> s, ConstantsProvider<T> K){ {f(s,K)}->Wrench<T>; };
template<class T> class Geopotential { Vector3<T> accel_ecef(r,mu,Re) const; ForceFn<T> as_force(K, Epoch<T>) const; };
template<class T,size_t S> StepResult<T> rk_step(const ButcherTableau<T,S>&, const State<T>&, dt, const AccelFn<T>&);

// P4 — coefficients
template<class T> class GravityField { TrackedValue<T> C(n,m), S(n,m); static GravityField egm2008(tol), wgs72(tol); };

// P5 — state & propagation
template<class P,class T> concept Propagatable = requires(P p, TrackedValue<T> t){ p.propagate(t); };
```

Every one is additive and behavior-preserving: the frozen functions are wrapped (e.g. `Geopotential::as_force`
reuses `sidereal_rotation`; `OrbitalModel<Sun,SR3>` reproduces the constant-rate model without touching the
inline `deep_space.h` SR3 arithmetic). Where a *value* would change (swapping the 2-term EoC for the Bessel
series), it is done only as the deliberate, astropy/DE430-gated L3 rebuild — never sold as a refactor.

## 3. Orphans waiting for a consumer (anti-dead-code resolution)
- `compute_third_body_rates` (zero consumers) → a **`ThirdBody<T>{ephemeris, μ}`** adapter that feeds the
  freshly-shipped `sun_meeus`/`moon_meeus` (λ,β) + `obliquity_iau2006` into both a Gauss-rate form and a
  Cartesian-acceleration `Force`. This is the single biggest *latent* integration — everything exists, nothing
  is wired.
- `obliquity_iau2006` (zero consumers) → consumed by `body_position` (P2) and `ThirdBody` (replacing the
  open-coded ecliptic→equatorial rotation in `third_body.h:98-100` with `frames::rot_x(ε)`).
- `kaula::inclination_function` (the de-duplicated home for the F_lmp polynomials `resonance.h` inlines) →
  grow the Gooding–King recursion for arbitrary (l,m,p).
- `attitude.h` (orphan island: 2nd RK4 + 3rd state type + 2nd inverse-inertia) → fold into `integrators/` +
  `dynamics/` (P3/P5), the highest-leverage flexible-side cleanup.

## 4. Recommended sequencing (highest value / lowest risk first)
1. **P1 `TrackedPolynomial` + `add_bound`** — exact, gated, retires 6+ copies + the dead `series.h::horner`;
   lands in the L1–L3 layers currently being written (low retrofit). *Do first.*
2. **P2 finish the ephemeris** — `OrbitalModel<Body>` over the existing `sun_meeus`/`moon_meeus` instances, then
   `body_position` to wire `obliquity`+`frames` (resolves two orphans). *In progress.*
3. **P4 `GravityField` → P3 `Geopotential` + `Force` concept** — modern-side, gated against the summed legacy
   forces (round-off tolerance, not bit — different recurrence). Wires `ThirdBody` (orphan).
4. **P5 `ButcherTableau` + state adapters + `Propagatable`** — deletes `attitude.h`'s duplicate RK4, wires
   adaptive stepping into `Propagator`.

## 5. What NOT to touch (skeptic's bottom line)
- The `TrackedValue` substrate and the algebraic wrappers — already optimally factored.
- The per-function transcendental error bounds in `tracked_value.h` — each is individually analytically
  verified; the line-savings of a bound-combinator do not justify perturbing 20 gated bounds.
- `Matrix3`→`Matrix<N>` genericity — no N≠3 consumer; an LU inverse would change the error budget.
- Anything on the frozen SGP4 path's value computation — wrap, never rewrite; size any unification gate against
  an *independent* oracle (the summed legacy forces / DE430 / ERFA), per the no-perceived-fidelity mandate.

*Full per-cluster evidence (file:line) in the four survey-agent transcripts of session 2026-06-06b.*
