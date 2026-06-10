#pragma once

/// @file twist.h
/// Twist: body-frame angular and linear velocity.
///
/// A twist Ω̂ = ω + ε v is the pure dual quaternion representation of a
/// body-frame velocity: angular component ω is the body angular velocity,
/// linear component v is the body linear velocity. Twists drive pose
/// evolution under
///
///   dM̂/dt = (1/2) M̂ Ω̂_pure,
///
/// where Ω̂_pure is the pure dual quaternion form returned by
/// `as_dual_quaternion()`.
///
/// Storage. The two 3-vectors are stored directly rather than as a pure
/// dual quaternion (which would have two unused scalar slots). Arithmetic
/// is element-wise on (angular, linear).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../math/dual_quaternion.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <cmath>

namespace dynamics {

/// Twist: body-frame angular + linear velocity as a pair of 3-vectors.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct Twist {
    math::Vector3<T> angular;  ///< Body angular velocity ω.
    math::Vector3<T> linear;   ///< Body linear velocity v.

    // --- Constructors ---

    /// Zero twist.
    Twist() : angular(), linear() {}

    /// From explicit angular and linear body velocities.
    Twist(const math::Vector3<T>& omega, const math::Vector3<T>& v)
        : angular(omega), linear(v) {}

    // --- Named constructors ---

    /// Zero twist (no motion).
    static Twist zero() { return Twist(); }

    /// Pure angular twist (no linear motion).
    static Twist pure_angular(const math::Vector3<T>& omega) {
        return Twist(omega, math::Vector3<T>());
    }

    /// Pure linear twist (no rotation).
    static Twist pure_linear(const math::Vector3<T>& v) {
        return Twist(math::Vector3<T>(), v);
    }

    // --- Conversion ---

    /// Lift to a pure dual quaternion (scalar parts zero). Used when an
    /// algebraic operation requires the dual quaternion form, e.g., the
    /// screw exponential exp_screw(this).
    math::DualQuaternion<T> as_dual_quaternion() const {
        return math::DualQuaternion<T>::from_screw(angular, linear);
    }

    // --- Arithmetic ---

    /// Component-wise addition of two twists.
    friend Twist operator+(const Twist& a, const Twist& b) {
        return Twist(a.angular + b.angular, a.linear + b.linear);
    }

    /// Component-wise subtraction of two twists.
    friend Twist operator-(const Twist& a, const Twist& b) {
        return Twist(a.angular - b.angular, a.linear - b.linear);
    }

    /// Scalar-twist multiplication (scalar on left).
    friend Twist operator*(const math::TrackedValue<T>& s, const Twist& w) {
        return Twist(s * w.angular, s * w.linear);
    }

    /// Scalar-twist multiplication (scalar on right).
    friend Twist operator*(const Twist& w, const math::TrackedValue<T>& s) {
        return Twist(w.angular * s, w.linear * s);
    }

    // --- Magnitude ---

    /// L1 norm of the six numeric values, Σ|componentᵢ.value| — the magnitude
    /// proxy the integrators use for their truncation-error bounds. The flat
    /// 6-term order (angular x,y,z then linear x,y,z) reproduces the integrators'
    /// prior inline abs-sum BIT-for-bit (used as a |y⁽ⁿ⁾| LTE proxy).
    T l1_norm() const {
        using std::abs;
        return abs(angular.x.value) + abs(angular.y.value) + abs(angular.z.value)
             + abs(linear.x.value)  + abs(linear.y.value)  + abs(linear.z.value);
    }

    /// L1 magnitude of the ANGULAR part alone — the rotational LTE proxy
    /// (runge_kutta_lie_group.md §6.2; exactly 0 for torque-free motion).
    T l1_norm_angular() const {
        using std::abs;
        return abs(angular.x.value) + abs(angular.y.value) + abs(angular.z.value);
    }

    /// L1 magnitude of the LINEAR part alone — the translational LTE proxy.
    T l1_norm_linear() const {
        using std::abs;
        return abs(linear.x.value) + abs(linear.y.value) + abs(linear.z.value);
    }
};

} // namespace dynamics
