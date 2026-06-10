#pragma once

/// @file inertia.h
/// Rigid-body inertia: mass + diagonal body-frame inertia tensor.
///
/// In the dual quaternion formulation of rigid-body dynamics, the dual
/// inertia operator H maps a body twist Ω̂ = ω + ε v to the body
/// generalized momentum:
///
///   H Ω̂ = (I_body · ω) + ε (m · v).
///
/// This implementation supports a diagonal body inertia tensor with
/// principal moments (I_xx, I_yy, I_zz). The full off-diagonal case
/// requires a Matrix3<T> infrastructure and is deferred to a future phase
/// (no existing application in the math library needs it yet).
///
/// Point-mass orbital propagation. Setting `principal_moments` to zero
/// degenerates the angular dynamics; the propagator must in that case be
/// configured to skip attitude integration. `acceleration_from_wrench` is
/// not valid on a point-mass inertia (division by zero); use the named
/// constructor `point_mass(m)` only with a propagator that will not call
/// the inverse operation.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../math/tracked_value.h"
#include "../math/vector3.h"
#include "twist.h"
#include "wrench.h"

namespace dynamics {

/// Rigid-body inertia: mass and diagonal principal moments.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct Inertia {
    math::TrackedValue<T> mass;          ///< Body mass m.
    math::Vector3<T> principal_moments;  ///< Diagonal body inertia
                                         ///< (I_xx, I_yy, I_zz).

    // --- Constructors ---

    /// Zero mass and moments (an uninitialized body; prefer the factories).
    Inertia() : mass(), principal_moments() {}

    /// From an explicit mass and diagonal principal moments.
    Inertia(const math::TrackedValue<T>& m, const math::Vector3<T>& I)
        : mass(m), principal_moments(I) {}

    // --- Named constructors ---

    /// Point mass: zero rotational inertia. Use only with a propagator
    /// configured to skip attitude dynamics.
    static Inertia point_mass(const math::TrackedValue<T>& m) {
        return Inertia(m, math::Vector3<T>());
    }

    /// Uniform solid sphere of radius r: I = (2/5) m r² on each principal
    /// axis.
    static Inertia uniform_sphere(const math::TrackedValue<T>& m,
                                  const math::TrackedValue<T>& r) {
        math::TrackedValue<T> two_fifths_m_r_sq =
            math::ratio<T>(2, 5) * m * r * r;
        return Inertia(m, math::Vector3<T>(two_fifths_m_r_sq,
                                           two_fifths_m_r_sq,
                                           two_fifths_m_r_sq));
    }

    /// Diagonal body inertia from three explicit principal moments.
    static Inertia diagonal(const math::TrackedValue<T>& m,
                            const math::TrackedValue<T>& I_xx,
                            const math::TrackedValue<T>& I_yy,
                            const math::TrackedValue<T>& I_zz) {
        return Inertia(m, math::Vector3<T>(I_xx, I_yy, I_zz));
    }

    // --- Operations ---

    /// Compute the body generalized momentum (L, p) of a given twist:
    ///   L_i = I_ii · ω_i,    p = m · v.
    /// Packed in Wrench-shape (torque slot holds L, force slot holds p).
    /// The Wrench type is reused for shape; momentum is semantically dual
    /// to twist, not the same kind of object as an applied force.
    Wrench<T> momentum_of(const Twist<T>& w) const {
        math::Vector3<T> L(
            principal_moments.x * w.angular.x,
            principal_moments.y * w.angular.y,
            principal_moments.z * w.angular.z
        );
        math::Vector3<T> p(
            mass * w.linear.x,
            mass * w.linear.y,
            mass * w.linear.z
        );
        return Wrench<T>(L, p);
    }

    /// Apply the inverse inertia: given a body wrench (τ, F), return the
    /// body acceleration (ω_dot, v_dot):
    ///   ω_dot_i = τ_i / I_ii
    ///   v_dot   = F / m.
    ///
    /// Zero-moment case (point-mass mode). When `principal_moments.i` is
    /// exactly zero, the body has no rotational inertia about axis i.
    /// With zero torque along that axis (gravity-only orbital propagation
    /// has this property), the angular acceleration is taken to be zero.
    /// Nonzero torque with zero moment is physically undefined (infinite
    /// α); the result is set to a max-bound `TrackedValue` to signal
    /// catastrophe per REQ-EF-9 / AUD-EF-9.
    ///
    /// Precondition: `mass.value > 0`.
    Twist<T> acceleration_from_wrench(const Wrench<T>& f) const {
        math::TrackedValue<T> (*safe_div)(
            const math::TrackedValue<T>&, const math::TrackedValue<T>&) =
            [](const math::TrackedValue<T>& num,
               const math::TrackedValue<T>& den) -> math::TrackedValue<T> {
            if (den.value == T(0)) {
                if (num.value == T(0)) {
                    // Documented zero-moment / zero-torque (point mass).
                    return math::exact<T>(0);
                }
                // Catastrophic regime: nonzero torque on zero moment.
                T inf = std::numeric_limits<T>::max();
                return math::TrackedValue<T>(inf, T(0), inf, T(0));
            }
            return num / den;
        };
        math::Vector3<T> alpha(
            safe_div(f.torque.x, principal_moments.x),
            safe_div(f.torque.y, principal_moments.y),
            safe_div(f.torque.z, principal_moments.z)
        );
        math::Vector3<T> a(
            f.force.x / mass,
            f.force.y / mass,
            f.force.z / mass
        );
        return Twist<T>(alpha, a);
    }

    /// Apply the inverse mass only: linear acceleration F/m. Use when the
    /// propagator skips attitude dynamics (point-mass mode).
    math::Vector3<T> linear_acceleration_from_force(
        const math::Vector3<T>& F) const {
        return math::Vector3<T>(F.x / mass, F.y / mass, F.z / mass);
    }
};

} // namespace dynamics
