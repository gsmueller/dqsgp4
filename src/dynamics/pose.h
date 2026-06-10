#pragma once

/// @file pose.h
/// Pose: a unit dual quaternion representing an element of SE(3).
///
/// Wraps math::DualQuaternion<T> with the unit-norm semantic constraint
/// (REQ-DQ-15). The constraint is not enforced at runtime — callers use
/// `normalized()` after operations that could cause drift, typically once
/// per integrator step (REQ-EF-15).
///
/// A pose encodes position and orientation jointly: rotation in the real
/// quaternion q_r, translation-coupled term in the dual quaternion q_d.
/// The translation 3-vector is recovered as t = 2 (q_d · q_r*).vector().
///
/// Composition. Pose composition is the dual quaternion product, with
/// right-to-left semantics: `(a * b)` applies b first, then a.
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

/// Pose: rigid-body position + orientation, encoded as a unit dual
/// quaternion.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct Pose {
    math::DualQuaternion<T> M;  ///< Underlying unit dual quaternion.

    // --- Constructors ---

    /// Default-construct as the SE(3) identity (origin, identity rotation).
    Pose() : M(math::DualQuaternion<T>::identity()) {}

    /// Wrap an existing dual quaternion. Caller is responsible for unit
    /// constraint; no renormalization is performed.
    explicit Pose(const math::DualQuaternion<T>& dq) : M(dq) {}

    // --- Named constructors ---

    /// SE(3) identity: zero translation, identity rotation.
    static Pose identity() {
        return Pose(math::DualQuaternion<T>::identity());
    }

    /// Build a pose from a unit rotation quaternion q_r and a translation
    /// 3-vector t. Caller is responsible for the unit-norm of q_r.
    static Pose from_rotation_translation(const math::Quaternion<T>& q_r,
                                          const math::Vector3<T>& t) {
        return Pose(math::DualQuaternion<T>::from_pose(q_r, t));
    }

    // --- Component access ---

    /// Rotation quaternion (q_r).
    math::Quaternion<T> rotation() const { return M.rotation(); }

    /// Translation 3-vector recovered from the pose.
    math::Vector3<T> translation() const { return M.translation(); }

    /// Underlying dual quaternion (read-only access).
    const math::DualQuaternion<T>& as_dual_quaternion() const { return M; }

    // --- Composition ---

    /// Compose two poses: `(a * b)` applies b first, then a.
    friend Pose operator*(const Pose& a, const Pose& b) {
        return Pose(a.M * b.M);
    }

    /// Multiplicative inverse. Precondition: this is a unit pose.
    Pose inverse() const { return Pose(M.inverse()); }

    /// SE(3) retraction back onto the manifold (REQ-EF-15). Used after
    /// integration to project drift back to the unit-pose surface.
    Pose normalized() const { return Pose(M.normalized()); }

    // --- Action ---

    /// Apply this pose to a 3-point: p' = R p + t.
    /// Precondition: this is a unit pose.
    math::Vector3<T> apply(const math::Vector3<T>& p) const {
        return M.apply(p);
    }

    /// Apply only the rotation component to a 3-direction (no
    /// translation). Precondition: this is a unit pose.
    math::Vector3<T> apply_direction(const math::Vector3<T>& v) const {
        return M.apply_direction(v);
    }

    // --- Interpolation ---

    /// Screw linear interpolation between two unit poses at parameter
    /// t ∈ [0, 1]: `sclerp(a, b, 0) == a`, `sclerp(a, b, 1) == b`. Traverses
    /// the screw geodesic at constant rate.
    static Pose sclerp(const Pose& a, const Pose& b,
                       const math::TrackedValue<T>& t) {
        return Pose(math::DualQuaternion<T>::sclerp(a.M, b.M, t));
    }
};

} // namespace dynamics
