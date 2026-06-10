#pragma once

/// @file propagator.h
/// Propagator: orchestrator that advances a `State<T>` through time under
/// a list of injected force lambdas and an injected integrator.
///
/// The propagator owns the constants provider, the body inertia, the
/// force list, and the integrator at construction. It exposes:
///
///   step(state, dt)                       — single integrator step
///   propagate_to(state, t_target, dt_max) — multi-step to a target time
///   compute_acceleration(state)           — body acceleration callback
///
/// Each step:
///
///   1. `compute_acceleration` sums every force lambda evaluated at the
///      stage state and divides by the inertia to obtain a body twist
///      of accelerations (REQ-PR-2).
///   2. The integrator advances pose, twist, and time using the
///      body-acceleration callback (REQ-PR-3).
///   3. The integrator's `lie_advance_pose` applies the SE(3) retraction
///      that keeps the pose on the unit-DQ manifold (REQ-EF-15,
///      REQ-PR-9).
///
/// No hardcoded force model. No hardcoded constants. No hardcoded
/// integrator scheme. Every dependency arrives through the constructor
/// (REQ-PR-1; OBJ-1, OBJ-4, OBJ-5, CON-3).
///
/// Force lambda signature is `(State, ConstantsProvider) → Wrench`. The
/// propagator passes its own `K_` member to each lambda at evaluation
/// time, so force lambdas do not capture `K` themselves (REQ-PR-5).
///
/// @par Usage
/// @code
/// constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
/// dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(math::exact<T>(1));
///
/// std::vector<dynamics::ForceFn<T>> forces;
/// forces.push_back([](const dynamics::State<T>& s, const constants::ConstantsProvider<T>& KK) {
///     return forces::gravity_central(s, KK);   // + gravity_J2, drag, ...
/// });
///
/// dynamics::IntegratorFn<T> integrator =
///     [](const dynamics::State<T>& y0, const math::TrackedValue<T>& dt,
///        const integrators::AccelFn<T>& accel) {
///         return integrators::runge_kutta_4(y0, dt, accel);
///     };
///
/// dynamics::Propagator<T> prop(K, inertia, std::move(forces), integrator);
/// dynamics::State<T> sf = prop.propagate_to(state0, t_target, dt_max);
/// @endcode
/// A full runnable example (with the TLE->state bridge) is in
/// `examples/quickstart.cpp`; see also README.md and design/propagator_choice.md.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-6, AUD-EF-7.

#include "../constants/constants_provider.h"
#include "../integrators/runge_kutta.h"
#include "../integrators/runge_kutta_fehlberg.h"
#include "../math/tracked_value.h"
#include "inertia.h"
#include "state.h"
#include "twist.h"
#include "wrench.h"

#include <functional>
#include <utility>
#include <vector>

namespace dynamics {

/// Force lambda signature: takes a State and a ConstantsProvider, returns
/// the body-frame Wrench contributed by this force at that state.
template<typename T>
using ForceFn = std::function<
    Wrench<T>(const State<T>&,
              const constants::ConstantsProvider<T>&)>;

/// Integrator signature: takes a state, a step size, and a
/// body-acceleration callback; returns the advanced state.
template<typename T>
using IntegratorFn = std::function<
    State<T>(const State<T>&,
             const math::TrackedValue<T>&,
             const integrators::AccelFn<T>&)>;


/// Propagator: holds the model and advances states through time.
///
/// @tparam T  Underlying numeric type.
template<typename T>
class Propagator {
public:
    /// Construct a propagator with all dependencies passed in.
    ///
    /// @param K          Constants provider (moved into the propagator).
    /// @param inertia    Body inertia (moved into the propagator).
    /// @param forces     List of force lambdas (moved into the propagator).
    /// @param integrator Integrator function (moved into the propagator).
    Propagator(constants::ConstantsProvider<T> K,
               Inertia<T> inertia,
               std::vector<ForceFn<T>> forces,
               IntegratorFn<T> integrator)
        : K_(std::move(K)),
          inertia_(std::move(inertia)),
          forces_(std::move(forces)),
          integrator_(std::move(integrator)) {}

    // --- Read-only accessors ---

    /// The injected constants provider.
    const constants::ConstantsProvider<T>& constants() const { return K_; }

    /// The injected body inertia.
    const Inertia<T>& inertia() const { return inertia_; }

    /// The injected force list (read-only).
    const std::vector<ForceFn<T>>& forces() const { return forces_; }

    // --- Acceleration callback ---

    /// Sum every force lambda evaluated at `state` and divide by the
    /// inertia to obtain the body acceleration twist (REQ-PR-2).
    Twist<T> compute_acceleration(const State<T>& state) const {
        Wrench<T> total = Wrench<T>::zero();
        for (const ForceFn<T>& f : forces_) {
            total = total + f(state, K_);
        }
        return inertia_.acceleration_from_wrench(total);
    }

    // --- Stepping ---

    /// Advance the state by a single integrator step of size `dt`
    /// (REQ-PR-3).
    State<T> step(const State<T>& y0,
                  const math::TrackedValue<T>& dt) const {
        integrators::AccelFn<T> accel_fn =
            [this](const State<T>& s) -> Twist<T> {
                return this->compute_acceleration(s);
            };
        return integrator_(y0, dt, accel_fn);
    }

    /// Multi-step propagation from `y0` to time `t_target` using steps of
    /// size at most `dt_max`. The final step is shortened so the elapsed
    /// time lands on `t_target` (REQ-PR-4).
    ///
    /// Precondition: `t_target.value > y0.time.value`.
    State<T> propagate_to(const State<T>& y0,
                          const math::TrackedValue<T>& t_target,
                          const math::TrackedValue<T>& dt_max) const {
        State<T> y = y0;
        while (y.time.value < t_target.value) {
            math::TrackedValue<T> remaining = t_target - y.time;
            math::TrackedValue<T> step_size =
                (remaining.value < dt_max.value) ? remaining : dt_max;
            y = step(y, step_size);
        }
        return y;
    }

    /// Adaptive RKF7(8) propagation through the facade (REQ-IN-6; R3b,
    /// design/derivations/dq_propagator_facade.md §4): delegates to the
    /// AD1-gated `integrators::rkf78_propagate_adaptive` with this propagator's
    /// `compute_acceleration` as the callback — BIT-IDENTICAL to calling the
    /// standalone loop with the same callback (it IS that call). Returns the
    /// final state plus accepted/rejected step diagnostics.
    integrators::AdaptiveResult<T> propagate_adaptive(
        const State<T>& y0,
        const math::TrackedValue<T>& t_target,
        const math::TrackedValue<T>& dt_initial,
        const T& tol,
        const T& dt_min) const {
        integrators::AccelFn<T> accel_fn =
            [this](const State<T>& s) -> Twist<T> {
                return this->compute_acceleration(s);
            };
        return integrators::rkf78_propagate_adaptive<T>(
            y0, t_target, dt_initial, accel_fn, tol, dt_min);
    }

private:
    constants::ConstantsProvider<T> K_;
    Inertia<T> inertia_;
    std::vector<ForceFn<T>> forces_;
    IntegratorFn<T> integrator_;
};

} // namespace dynamics
