#pragma once

/// @file wrench.h
/// Wrench: body-frame torque and force.
///
/// A wrench Ŵ = τ + ε F is the pure dual quaternion representation of a
/// body-frame generalized force: angular component τ is the body torque,
/// linear component F is the body force. Wrenches drive twist evolution
/// under the SE(3) Newton–Euler equation
///
///   dΩ̂/dt = Î⁻¹ (Ŵ − Ω̂ × Î Ω̂),
///
/// where Î is the body inertia (real = inertia tensor, dual = m · 𝟙) and
/// × is the dual quaternion commutator.
///
/// Storage. The two 3-vectors are stored directly. Arithmetic is
/// element-wise on (torque, force). Total wrench at a body is the sum of
/// contributions from each force lambda passed to the propagator
/// (REQ-PR-2; implemented in dynamics::Propagator::compute_acceleration).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../math/dual_quaternion.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

namespace dynamics {

/// Wrench: body-frame torque + force as a pair of 3-vectors.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct Wrench {
    math::Vector3<T> torque;  ///< Body torque τ.
    math::Vector3<T> force;   ///< Body force F.

    // --- Constructors ---

    /// Zero wrench.
    Wrench() : torque(), force() {}

    /// From explicit body torque and force.
    Wrench(const math::Vector3<T>& tau, const math::Vector3<T>& F)
        : torque(tau), force(F) {}

    // --- Named constructors ---

    /// Zero wrench (no force, no torque).
    static Wrench zero() { return Wrench(); }

    /// Pure torque (no force).
    static Wrench pure_torque(const math::Vector3<T>& tau) {
        return Wrench(tau, math::Vector3<T>());
    }

    /// Pure force (no torque).
    static Wrench pure_force(const math::Vector3<T>& F) {
        return Wrench(math::Vector3<T>(), F);
    }

    // --- Conversion ---

    /// Lift to a pure dual quaternion (scalar parts zero). Used when an
    /// algebraic operation requires the dual quaternion form.
    math::DualQuaternion<T> as_dual_quaternion() const {
        return math::DualQuaternion<T>::from_screw(torque, force);
    }

    // --- Arithmetic ---

    /// Component-wise addition of two wrenches.
    friend Wrench operator+(const Wrench& a, const Wrench& b) {
        return Wrench(a.torque + b.torque, a.force + b.force);
    }

    /// Component-wise subtraction of two wrenches.
    friend Wrench operator-(const Wrench& a, const Wrench& b) {
        return Wrench(a.torque - b.torque, a.force - b.force);
    }

    /// Scalar-wrench multiplication (scalar on left).
    friend Wrench operator*(const math::TrackedValue<T>& s, const Wrench& w) {
        return Wrench(s * w.torque, s * w.force);
    }

    /// Scalar-wrench multiplication (scalar on right).
    friend Wrench operator*(const Wrench& w, const math::TrackedValue<T>& s) {
        return Wrench(w.torque * s, w.force * s);
    }
};

} // namespace dynamics
