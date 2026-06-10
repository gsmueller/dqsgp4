#pragma once

/// @file quaternion.h
/// Quaternion algebra over TrackedValue<T> components.
///
/// A quaternion q = w + x i + y j + z k is stored as four TrackedValue<T>
/// components. Arithmetic propagates the three-error budget (measurement,
/// precision, accuracy) through every operation by composing the underlying
/// TrackedValue<T> operator overloads. Closed-form algebra is used
/// throughout (REQ-EF-3). The Lie-algebra exp/log maps delegate the
/// small-argument Taylor branches to `small_angle_series.h`, where the
/// truncation bound is added to precision (REQ-EF-6).
///
/// Conventions.
///   - Hamilton product. The composition q1 q2 applied to a vector v as
///     (q1 q2) v (q1 q2)* rotates v first by q2 and then by q1. Matches
///     Boost.Geometry and Eigen.
///   - Active right-handed rotation. A unit quaternion rotates a 3-vector
///     by v' = q v q* with v promoted to a pure quaternion.
///   - Half-angle. A rotation by angle theta about unit axis n̂ is
///     q = cos(theta/2) + sin(theta/2) (n_x i + n_y j + n_z k).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-15, AUD-CC-16,
///   AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-2, AUD-EF-5, AUD-EF-7.

#include "small_angle_series.h"
#include "tracked_value.h"
#include "vector3.h"

#include <cmath>

namespace math {

/// Quaternion q = w + x i + y j + z k with error-tracked components.
///
/// No runtime unit-norm constraint is enforced; routines that require a unit
/// quaternion (rotate, log_unit, exp_pure) state so in their documentation.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct Quaternion {
    TrackedValue<T> w;  ///< Scalar component.
    TrackedValue<T> x;  ///< i component.
    TrackedValue<T> y;  ///< j component.
    TrackedValue<T> z;  ///< k component.

    // --- Constructors ---

    /// Zero quaternion (all four components zero).
    Quaternion() : w(), x(), y(), z() {}

    /// From explicit (w, x, y, z) components.
    Quaternion(const TrackedValue<T>& w_, const TrackedValue<T>& x_,
               const TrackedValue<T>& y_, const TrackedValue<T>& z_)
        : w(w_), x(x_), y(y_), z(z_) {}

    // --- Named constructors ---

    /// Multiplicative identity q = 1 + 0 i + 0 j + 0 k.
    static Quaternion identity() {
        return Quaternion(exact<T>(1), exact<T>(0), exact<T>(0), exact<T>(0));
    }

    /// Additive identity q = 0.
    static Quaternion zero() { return Quaternion(); }

    /// Pure quaternion: scalar part is zero, vector part is v.
    static Quaternion pure(const Vector3<T>& v) {
        return Quaternion(exact<T>(0), v.x, v.y, v.z);
    }

    /// Axis-angle construction. The axis must already be unit length; this
    /// constructor does not renormalize.
    ///
    ///   q = cos(theta/2) + sin(theta/2) (n_x i + n_y j + n_z k)
    static Quaternion from_axis_angle(const Vector3<T>& axis_unit,
                                      const TrackedValue<T>& angle) {
        TrackedValue<T> half = angle / exact<T>(2);
        TrackedValue<T> c = cos(half);
        TrackedValue<T> s = sin(half);
        return Quaternion(c, s * axis_unit.x, s * axis_unit.y, s * axis_unit.z);
    }

    // --- Component access ---

    /// Scalar component.
    TrackedValue<T> scalar() const { return w; }

    /// Vector component as a Vector3.
    Vector3<T> vector() const { return Vector3<T>(x, y, z); }

    // --- Arithmetic ---

    /// Component-wise addition.
    friend Quaternion operator+(const Quaternion& a, const Quaternion& b) {
        return Quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);
    }

    /// Component-wise subtraction.
    friend Quaternion operator-(const Quaternion& a, const Quaternion& b) {
        return Quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z);
    }

    /// Unary negation.
    Quaternion operator-() const {
        return Quaternion(-w, -x, -y, -z);
    }

    /// Hamilton product (non-commutative).
    ///
    ///   (w1 + v1)(w2 + v2)
    ///     = w1 w2 − v1 · v2
    ///       + w1 v2 + w2 v1
    ///       + v1 × v2
    ///
    /// Right-to-left composition: (q1 q2) applies q2 first, then q1.
    friend Quaternion operator*(const Quaternion& a, const Quaternion& b) {
        return Quaternion(
            a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
            a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
        );
    }

    /// Scalar-quaternion multiplication (scalar on left).
    friend Quaternion operator*(const TrackedValue<T>& s, const Quaternion& q) {
        return Quaternion(s*q.w, s*q.x, s*q.y, s*q.z);
    }

    /// Scalar-quaternion multiplication (scalar on right).
    friend Quaternion operator*(const Quaternion& q, const TrackedValue<T>& s) {
        return Quaternion(q.w*s, q.x*s, q.y*s, q.z*s);
    }

    // --- Magnitude and conjugation ---

    /// Quaternion conjugate q* = w − x i − y j − z k.
    Quaternion conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }

    /// Squared Euclidean magnitude: w² + x² + y² + z².
    TrackedValue<T> magnitude_squared() const {
        return w*w + x*x + y*y + z*z;
    }

    /// Euclidean magnitude: √(w² + x² + y² + z²).
    TrackedValue<T> magnitude() const {
        return sqrt(magnitude_squared());
    }

    /// Multiplicative inverse q⁻¹ = q* / |q|².
    /// For a unit quaternion this equals the conjugate; using `conjugate()`
    /// directly skips the squared-magnitude division.
    Quaternion inverse() const {
        TrackedValue<T> m2 = magnitude_squared();
        return Quaternion(w/m2, -x/m2, -y/m2, -z/m2);
    }

    /// Renormalize to unit magnitude. Used as a Lie-group retraction back
    /// onto the unit sphere after a numerical integration step (REQ-EF-15).
    Quaternion normalized() const {
        TrackedValue<T> m = magnitude();
        return Quaternion(w/m, x/m, y/m, z/m);
    }

    // --- Action ---

    /// Apply this (assumed unit) quaternion as a rotation to a 3-vector:
    ///   v' = q v q*  in pure-quaternion algebra.
    ///
    /// Implemented via the algebraic identity
    ///   v' = v + 2 w (q_v × v) + 2 q_v × (q_v × v)
    /// which uses fewer multiplications than the explicit triple product
    /// and keeps the result as a Vector3<T> for further composition.
    Vector3<T> rotate(const Vector3<T>& v) const {
        Vector3<T> qv(x, y, z);
        Vector3<T> qv_cross_v = qv.cross(v);
        TrackedValue<T> two = exact<T>(2);
        Vector3<T> two_qv_cross_v = two * qv_cross_v;
        return v + w * two_qv_cross_v + qv.cross(two_qv_cross_v);
    }

    // --- Lie algebra: exp and log ---

    /// Exponential of a pure quaternion u = (0, v):
    ///
    ///   exp(0 + v) = cos|v| + sinc|v| · v.
    ///
    /// Used for unit-quaternion integration: if dq/dt = (1/2) q ω with ω a
    /// pure body angular-velocity quaternion, the Lie-group step is
    ///   q ← q · exp_pure((dt/2) ω_body)
    /// which keeps q exactly on the unit sphere up to floating-point error.
    ///
    /// Inverse of `log_unit` on the half-angle ball (|v| ≤ π/2); see
    /// REQ-DQ-12.
    static Quaternion exp_pure(const Vector3<T>& v) {
        TrackedValue<T> theta_sq = v.x*v.x + v.y*v.y + v.z*v.z;
        TrackedValue<T> theta = sqrt(theta_sq);
        TrackedValue<T> sinc_th = taylor_sinc(theta, theta_sq);
        TrackedValue<T> c = cos(theta);
        return Quaternion(c, sinc_th*v.x, sinc_th*v.y, sinc_th*v.z);
    }

    /// Logarithm of a unit quaternion, returned as the half-angle rotation
    /// vector (the vector part of the pure-quaternion logarithm).
    ///
    /// If q = cos(theta/2) + sin(theta/2) · n̂, returns (theta/2) · n̂.
    /// Doubling the result yields the Rodrigues rotation vector
    /// theta · n̂.
    ///
    /// Shortest-path convention. q and −q represent the same SO(3)
    /// rotation; when w.value < 0, the implementation operates on −q so
    /// the reported rotation is the shorter of the two equivalent angles.
    /// Inverse of `exp_pure` on the half-angle ball; see REQ-DQ-12.
    Vector3<T> log_unit() const {
        TrackedValue<T> w_e, x_e, y_e, z_e;
        if (w.value < T(0)) {
            w_e = -w; x_e = -x; y_e = -y; z_e = -z;
        } else {
            w_e = w; x_e = x; y_e = y; z_e = z;
        }
        TrackedValue<T> qv_norm_sq = x_e*x_e + y_e*y_e + z_e*z_e;
        TrackedValue<T> qv_norm = sqrt(qv_norm_sq);
        TrackedValue<T> scale =
            taylor_half_angle_scale(qv_norm, qv_norm_sq, w_e);
        return Vector3<T>(scale*x_e, scale*y_e, scale*z_e);
    }
};

} // namespace math
