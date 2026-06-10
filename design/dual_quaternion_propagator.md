# Dual Quaternion Propagator

## Goal

A general-purpose six-DOF (position + attitude) state propagator built on the
dual quaternion representation of SE(3). The propagator is generic over:

- **Numeric precision** — `double`, `cpp_bin_float_50`, or any Boost
  multiprecision type, all wrapped in `TrackedValue<T>` for rigorous
  three-error tracking (measurement / precision / accuracy).
- **Force and torque models** — every contribution to the body wrench
  (point-mass gravity, $J_2$/$J_4$/$J_n$ zonals, full spherical harmonics,
  third-body, atmospheric drag, solar radiation pressure, magnetic torques,
  control inputs) is injected as a `std::function` lambda.
- **Fundamental constants** — Earth parameters, time system, astronomical
  ephemerides, and physical constants are bundled in a `ConstantsProvider<T>`
  and injected at construction.
- **Integrator** — RK4, adaptive RKF7(8), Dormand–Prince 8(7), or any user
  scheme can be injected as a function object.

The propagator itself contains **no constants and no hardcoded formulas**.
It is a state-evolution engine.

## State

The state of a single rigid body is:

| Symbol | Type | Meaning |
|---|---|---|
| $\hat{M}$ | unit dual quaternion | Pose: position + attitude |
| $\hat{\Omega}$ | pure dual quaternion | Twist: angular + linear velocity, body frame |
| $t$ | TrackedValue | Time since epoch |

The unit dual quaternion $\hat{M} = q_r + \varepsilon\, q_d$ encodes:

- $q_r$ — rotation quaternion from body to world: a world-frame vector is
  $v_{\text{world}} = q_r\, v_{\text{body}}\, q_r^*$.
- $q_d = \tfrac{1}{2}\, t\, q_r$ where $t$ is the pure quaternion of the
  body's position (in the world frame).

The twist $\hat{\Omega} = \omega + \varepsilon\, v$ encodes:

- $\omega$ — body angular velocity as a pure quaternion.
- $v$ — body linear velocity as a pure quaternion; the world-frame velocity
  is $v_{\text{world}} = q_r\, v\, q_r^*$.

The (pose, twist) pair represents the same content as a $4\times 4$ SE(3)
matrix plus a 6-vector twist, but with constant-time composition, exact
unit-norm enforcement via Lie-group retraction, and smooth interpolation
via Sclerp.

## Time evolution

For a rigid body with mass $m$, principal inertia tensor $I$ (body frame),
and applied wrench $\hat{W} = \tau + \varepsilon\, F$ (body frame):

$$\dot{\hat{M}} = \tfrac{1}{2}\, \hat{M}\, \hat{\Omega}$$

$$\dot{\hat{\Omega}} = \hat{I}^{-1}\!\left(\hat{W} - \hat{\Omega} \times \hat{I}\, \hat{\Omega}\right)$$

where $\hat{I} = I + \varepsilon\, m\, \mathbb{1}$ is the dual inertia and
$\times$ denotes the dual quaternion cross product (commutator).

For a point-mass orbital body, the attitude dynamics decouple from the
translation: the orbit is governed entirely by the linear part of $\hat{W}$.
But the same engine handles spinning satellites, propulsive maneuvers, and
attitude-coupled drag without architectural change — pose and twist evolve
together via a single ODE.

## Wrench model (lambda injection)

Every force / torque contribution is a lambda of fixed signature:

```cpp
template<typename T>
using WrenchFn = std::function<
    DualQuaternion<T>(
        const Pose<T>&             M,
        const Twist<T>&            Omega,
        const TrackedValue<T>&     t,
        const ConstantsProvider<T>& K
    )>;
```

The propagator sums contributions before each integrator stage:

```cpp
DualQuaternion<T> W_total = DualQuaternion<T>::zero();
for (const auto& fn : wrenches) W_total = W_total + fn(M, Omega, t, K);
```

A caller composes a force model as a list:

- `gravity_central(GM)` — Newtonian point-mass
- `gravity_zonal_J2(GM, R_E, J2)` — added oblateness
- `gravity_zonal_Jn({J2, J3, J4, ...})`
- `gravity_spherical_harmonics(N_max, M_max, coeffs)`
- `drag_exponential(rho0, h0, H_scale, C_d, A, m)`
- `drag_density_model(density_fn, C_d, A, m)`
- `third_body(ephemeris_fn, GM_third)`
- `radiation_pressure(C_R, A, m, sun_position_fn)`
- `attitude_control(K_p, K_d, target_pose_fn)`

No combination is hardcoded.

## Constants provider

```cpp
template<typename T>
struct ConstantsProvider {
    EquipotentialEllipsoid<T> earth;       // a, 1/f, GM, omega + derived
    GravityField<T>           gravity;     // J_n harmonics, C_nm/S_nm
    SolarSystem<T>            astronomy;   // Sun/Moon ephemerides
    TimeSystem<T>             time;        // TAI/UTC/UT1/TT/GPS
    FundamentalConstants<T>   fundamentals;// c, G, ...

    static ConstantsProvider wgs84();
    static ConstantsProvider egm2008();
    static ConstantsProvider jgm3();
    static ConstantsProvider custom(/* user inputs */);
};
```

Each field carries the three-error budget. Swapping providers swaps every
downstream propagation result by construction — no recompile, no global
state, no hidden constants.

## Integrator

```cpp
template<typename T>
using Integrator = std::function<
    State<T>(
        const State<T>&        y0,
        const TrackedValue<T>& dt,
        const Derivative<T>&   f
    )>;
```

Standard implementations provided:

- `runge_kutta_4<T>` — fixed-step RK4 with unit-DQ renormalization
  (Lie-group retraction) after each step.
- `rkf78<T>` — adaptive Runge–Kutta–Fehlberg 7(8) with embedded error
  estimate driving step-size control.
- `dormand_prince_8<T>` — Dormand–Prince 8(7) for very high accuracy.
- `symplectic_leapfrog<T>` — for energy-preservation in conservative
  problems (no drag).

The unit-norm constraint on $\hat{M}$ is enforced by retraction after each
accepted step: $\hat{M} \leftarrow \hat{M}\,/\,\|\hat{M}\|$. This adds at
most $O(\varepsilon_{\text{machine}})$ drift per step and keeps the state
exactly on the SE(3) manifold to machine precision.

## Error propagation

Every quantity ($\hat{M}$, $\hat{\Omega}$, $\hat{W}$, $t$) is built from
`TrackedValue<T>` components. The three errors propagate through every
arithmetic operation. After propagation, the caller can query any output:

```cpp
auto pos = state.M.position();          // Vector3<TrackedValue<T>>
pos.x.errors.measurement;               // propagated input uncertainty
pos.x.errors.precision;                 // numerical roundoff accumulated
pos.x.errors.accuracy;                  // model truncation (harmonics, etc.)
pos.x.reliable_digits();                // meaningful decimal digits left
```

This is the primary mechanism by which "extension of accuracy and precision"
is delivered: the same propagator code, instantiated on a wider `T` and a
more complete force model, will report tighter errors automatically. The
caller can prove (not just hope) that a given prediction is accurate to a
stated tolerance.

## Module map

```
src/
  math/
    quaternion.h            NEW   unit + general quaternions
    dual_number.h           NEW   a + εb, ε² = 0
    dual_quaternion.h       NEW   SE(3) representation
  dynamics/                 NEW
    pose.h                        wraps unit dual quaternion
    twist.h                       wraps pure dual quaternion
    inertia.h                     dual inertia tensor
    wrench.h                      force/torque dual quaternion
    state.h                       (pose, twist, time, metadata)
    derivative.h                  d/dt of state
    propagator.h                  top-level integration loop
    forces/
      gravity_central.h
      gravity_zonal.h
      gravity_spherical_harmonics.h
      drag.h
      third_body.h
      radiation_pressure.h
    integrators/
      runge_kutta.h               RK4
      rkf78.h                     adaptive 7(8)
      dormand_prince.h            DOPRI8(7)
  constants/                NEW
    constants_provider.h
    wgs84_provider.h
    egm2008_provider.h
    gravity_field.h
    fundamental_constants.h
    time_system.h
```

## Coexistence with existing code

The existing `src/sgp4/`, `src/perturbation/`, etc. are unchanged. The new
dual quaternion propagator shares only the `math/` foundation (TrackedValue,
Vector3, factorials, etc.). The existing `EquipotentialEllipsoid<T>` is
reused as the `earth` field of `ConstantsProvider<T>` — its computation of
$e^2$, $b$, $q_0$, $J_{2n}$, gravities from $(a, 1/f, GM, \omega)$ is
exactly what the new gravity-field models need as input.

## Implementation phases

| # | Phase | Deliverable |
|---|---|---|
| 1 | Algebra | `quaternion.h`, `dual_number.h`, `dual_quaternion.h` + tests |
| 2 | State | `pose.h`, `twist.h`, `state.h`, `wrench.h`, `inertia.h` + tests |
| 3 | Constants | `constants_provider.h`, `wgs84_provider.h` + tests |
| 4 | Forces | `gravity_central.h`, `gravity_zonal.h` + tests |
| 5 | Integrator | `runge_kutta.h` (RK4) + tests |
| 6 | Propagator | `propagator.h` end-to-end + smoke test |
| 7 | Extensions | RKF78, full spherical harmonics, drag, third-body, SRP |
| 8 | Docs | Doxygen `@page` for the module, examples directory |

## References

- E. Pennestrì, P. P. Valentini, *Linear dual algebra algorithms and their
  application to kinematics*, Multibody Dyn. (2009).
- B. Kenwright, *A beginners guide to dual-quaternions*, WSCG 2012.
- Y. Wang, C. Chen, *Dual-quaternion-based satellite pose tracking with
  finite-time convergence*, IEEE Trans. Aerospace (2017).
- G. S. Chirikjian, *Stochastic Models, Information Theory, and Lie Groups,
  Vol. II*, Birkhäuser (2012) — Lie-group integrators on SE(3).
