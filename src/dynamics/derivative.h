#pragma once

/// @file derivative.h
/// Derivative: the time-derivative of `State<T>` produced by a force model.
///
/// In the dual quaternion formulation, the time derivative of the state
/// (pose, twist, time) breaks into three pieces:
///
///   dM̂/dt = (1/2) M̂ · Ω̂_pure        (governed by the current twist)
///   dΩ̂/dt = body_acceleration         (governed by the applied wrench)
///   dt/dt = 1.
///
/// The pose-derivative is determined by the current twist itself, so it
/// is not stored as a separate field here — the integrator advances the
/// pose using the current state's twist via the screw exponential. What a
/// force model produces is the body acceleration (ω_dot, v_dot), packaged
/// as a Twist of accelerations.
///
/// Arithmetic. Derivatives form a vector space under addition and scalar
/// multiplication; this is what RK4 / RKF7(8) and similar integrators
/// combine across stages.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../math/tracked_value.h"
#include "twist.h"

namespace dynamics {

/// Time derivative of a State<T>.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct Derivative {
    Twist<T> acceleration;             ///< dΩ̂/dt (body angular + linear).
    math::TrackedValue<T> time_rate;   ///< dt/dt (= 1 in the standard
                                       ///< parameterization).

    // --- Constructors ---

    /// Zero acceleration with the standard time parameterization (dt/dt = 1).
    Derivative() : acceleration(), time_rate(math::exact<T>(1)) {}

    /// From an explicit acceleration twist and time rate.
    Derivative(const Twist<T>& a, const math::TrackedValue<T>& dt_dt)
        : acceleration(a), time_rate(dt_dt) {}

    // --- Named constructors ---

    /// Build from a body acceleration with the standard parameterization
    /// dt/dt = 1.
    static Derivative from_acceleration(const Twist<T>& a) {
        return Derivative(a, math::exact<T>(1));
    }

    /// Zero acceleration with the standard time parameterization dt/dt = 1.
    ///
    /// NOTE: this is NOT the additive identity of the underlying direct-sum
    /// vector space se(3) × ℝ (which would have time_rate = 0). The name
    /// reflects the function's actual purpose: "no acceleration applied,
    /// time still advances at the standard rate." Renamed from `zero()`
    /// to make the non-identity nature explicit; cf. audit card
    /// `design/audit/theoretical_basis_audit/derivative.md` Card 4.
    static Derivative standard_acceleration_zero() {
        return Derivative(Twist<T>::zero(), math::exact<T>(1));
    }

    // --- Arithmetic ---

    /// Component-wise addition (used by RK stage combiners).
    friend Derivative operator+(const Derivative& a, const Derivative& b) {
        return Derivative(
            a.acceleration + b.acceleration,
            a.time_rate + b.time_rate
        );
    }

    /// Component-wise subtraction.
    friend Derivative operator-(const Derivative& a, const Derivative& b) {
        return Derivative(
            a.acceleration - b.acceleration,
            a.time_rate - b.time_rate
        );
    }

    /// Scalar-derivative multiplication (scalar on left).
    friend Derivative operator*(const math::TrackedValue<T>& s,
                                const Derivative& d) {
        return Derivative(s * d.acceleration, s * d.time_rate);
    }

    /// Scalar-derivative multiplication (scalar on right).
    friend Derivative operator*(const Derivative& d,
                                const math::TrackedValue<T>& s) {
        return Derivative(d.acceleration * s, d.time_rate * s);
    }
};

} // namespace dynamics
