#pragma once

/// @file dual_quaternion.h
/// Dual quaternion algebra over Quaternion<T> components.
///
/// A dual quaternion M̂ = q_r + ε q_d is a pair of quaternions combined under
/// ε² = 0. The library uses dual quaternions in two semantic roles:
///
///   - Pose. A unit dual quaternion encodes an element of SE(3): rotation in
///     q_r, translation-coupled term in q_d. The body-frame translation
///     vector is recovered as t = 2 · (q_d · q_r*).vector() for unit q_r.
///
///   - Twist / wrench. A pure dual quaternion (q_r.scalar() = 0,
///     q_d.scalar() = 0) encodes a body-frame screw: angular part in q_r,
///     linear part in q_d.
///
/// Both roles share the same algebraic type. Semantic wrappers (Pose, Twist,
/// Wrench) live in src/dynamics/.
///
/// Conventions.
///   - Product. (q_r1 + ε q_d1)(q_r2 + ε q_d2)
///       = q_r1 q_r2 + ε(q_r1 q_d2 + q_d1 q_r2), since ε² = 0.
///   - Three conjugates (REQ-DQ-16):
///       quaternion conjugate   M̂*  = q_r*  + ε q_d*
///       dual-number conjugate  M̂ε  = q_r   − ε q_d
///       combined conjugate     M̂♯  = q_r*  − ε q_d*
///     M̂* is the multiplicative inverse of a unit pose (REQ-DQ-15); M̂♯ is
///     the conjugate used in the point action M̂ P̂ M̂♯ (REQ-DQ-17).
///
/// Lie algebra. For a pure dual quaternion Û = u + ε v (with u, v lifted
/// from 3-vectors), the exponential closes to
///   exp(Û) = cos(θ̂) + sinc(θ̂) · Û
/// with the dual angle θ̂ = θ + ε d, θ = |u|, d = (u · v) / θ. Component
/// form (REQ-DQ-18):
///   q_r = cos θ + sinc θ · u                        (= Quaternion::exp_pure(u))
///   q_d.scalar() = −(u · v) sinc θ
///   q_d.vector() = sinc θ · v + (u · v) · β · u,  β = (cos θ − sinc θ) / θ²
/// The log inverts this on the half-angle ball (|u| ≤ π/2).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-15, AUD-CC-16,
///   AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-2, AUD-EF-5, AUD-EF-7.

#include "dual_number.h"
#include "quaternion.h"
#include "small_angle_series.h"
#include "tracked_value.h"
#include "vector3.h"

namespace math {

/// Dual quaternion M̂ = q_r + ε q_d with q_r, q_d : Quaternion<T>.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct DualQuaternion {
    Quaternion<T> real;  ///< Rotation / angular part (multiplies 1).
    Quaternion<T> dual;  ///< Translation / linear part (multiplies ε).

    // --- Constructors ---

    /// Zero dual quaternion (both parts zero).
    DualQuaternion() : real(), dual() {}

    /// From explicit real (rotation) and dual (translation) parts.
    DualQuaternion(const Quaternion<T>& r, const Quaternion<T>& d)
        : real(r), dual(d) {}

    // --- Named constructors ---

    /// Additive identity (both parts zero).
    static DualQuaternion zero() { return DualQuaternion(); }

    /// Multiplicative identity (real = quaternion identity, dual = zero).
    static DualQuaternion identity() {
        return DualQuaternion(Quaternion<T>::identity(), Quaternion<T>::zero());
    }

    /// Build a pose dual quaternion from a unit rotation quaternion q_r and
    /// a translation 3-vector t:
    ///
    ///   M̂ = q_r + ε · (1/2) · t · q_r,
    ///
    /// where t is promoted to a pure quaternion. Caller is responsible for
    /// the unit-norm of q_r; this routine does not renormalize.
    static DualQuaternion from_pose(const Quaternion<T>& q_r,
                                    const Vector3<T>& t) {
        Quaternion<T> t_pure = Quaternion<T>::pure(t);
        Quaternion<T> q_d = (t_pure * q_r) * ratio<T>(1, 2);
        return DualQuaternion(q_r, q_d);
    }

    /// Build a pure (twist or wrench) dual quaternion from an angular
    /// 3-vector ω and a linear 3-vector v:
    ///
    ///   Ω̂ = ω_pure + ε v_pure.
    static DualQuaternion from_screw(const Vector3<T>& omega,
                                     const Vector3<T>& v) {
        return DualQuaternion(Quaternion<T>::pure(omega),
                              Quaternion<T>::pure(v));
    }

    // --- Component access ---

    /// Rotation quaternion (real part of a pose).
    Quaternion<T> rotation() const { return real; }

    /// Translation 3-vector recovered from a pose dual quaternion:
    ///   t = 2 · (q_d · q_r*).vector()
    /// for unit M̂. Result is meaningless for non-unit M̂.
    Vector3<T> translation() const {
        Quaternion<T> td_over_2 = dual * real.conjugate();
        Vector3<T> v = td_over_2.vector();
        TrackedValue<T> two = exact<T>(2);
        return Vector3<T>(two * v.x, two * v.y, two * v.z);
    }

    /// Angular 3-vector of a pure (twist) dual quaternion.
    Vector3<T> angular() const { return real.vector(); }

    /// Linear 3-vector of a pure (twist) dual quaternion.
    Vector3<T> linear() const { return dual.vector(); }

    // --- Arithmetic ---

    /// Component-wise addition.
    friend DualQuaternion operator+(const DualQuaternion& a,
                                    const DualQuaternion& b) {
        return DualQuaternion(a.real + b.real, a.dual + b.dual);
    }

    /// Component-wise subtraction.
    friend DualQuaternion operator-(const DualQuaternion& a,
                                    const DualQuaternion& b) {
        return DualQuaternion(a.real - b.real, a.dual - b.dual);
    }

    /// Unary negation.
    DualQuaternion operator-() const {
        return DualQuaternion(-real, -dual);
    }

    /// Dual quaternion product:
    ///   (q_r1 + ε q_d1)(q_r2 + ε q_d2)
    ///     = q_r1 q_r2 + ε(q_r1 q_d2 + q_d1 q_r2),
    /// since ε² = 0. The Hamilton product on each quaternion factor is
    /// non-commutative.
    friend DualQuaternion operator*(const DualQuaternion& a,
                                    const DualQuaternion& b) {
        return DualQuaternion(
            a.real * b.real,
            a.real * b.dual + a.dual * b.real
        );
    }

    /// Scalar-DQ multiplication (scalar on left).
    friend DualQuaternion operator*(const TrackedValue<T>& s,
                                    const DualQuaternion& m) {
        return DualQuaternion(s * m.real, s * m.dual);
    }

    /// Scalar-DQ multiplication (scalar on right).
    friend DualQuaternion operator*(const DualQuaternion& m,
                                    const TrackedValue<T>& s) {
        return DualQuaternion(m.real * s, m.dual * s);
    }

    // --- Conjugation ---

    /// Quaternion-style conjugate (REQ-DQ-16):
    ///   M̂* = q_r* + ε q_d*.
    /// Equals the multiplicative inverse for a unit pose (REQ-DQ-15).
    DualQuaternion conjugate() const {
        return DualQuaternion(real.conjugate(), dual.conjugate());
    }

    /// Dual-number conjugate:
    ///   M̂ε = q_r − ε q_d.
    /// Flips the sign of the dual part.
    DualQuaternion dual_conjugate() const {
        return DualQuaternion(real, -dual);
    }

    /// Combined conjugate (REQ-DQ-16):
    ///   M̂♯ = q_r* − ε q_d*.
    /// Used in the point action M̂ P̂ M̂♯ (REQ-DQ-17); composition of
    /// `conjugate()` and `dual_conjugate()` in either order.
    DualQuaternion combined_conjugate() const {
        return DualQuaternion(real.conjugate(), -dual.conjugate());
    }

    // --- Magnitude ---

    /// Dual-valued squared magnitude:
    ///   ||M̂||²_dual = |q_r|² + ε · 2 · scalar(q_r · q_d*).
    /// Returned as a DualNumber<T>. For a unit pose, real = 1 and dual = 0.
    DualNumber<T> magnitude_squared() const {
        TrackedValue<T> real_sq = real.magnitude_squared();
        TrackedValue<T> qr_qd_dot =
            real.w * dual.w + real.x * dual.x + real.y * dual.y + real.z * dual.z;
        return DualNumber<T>(real_sq, exact<T>(2) * qr_qd_dot);
    }

    /// Multiplicative inverse, valid for a unit dual quaternion (a pose).
    ///
    /// For unit M̂, M̂⁻¹ = M̂* (the quaternion-style conjugate). General
    /// (non-unit) dual quaternion inverse requires inverting a dual-scalar
    /// magnitude; this method does not handle that case.
    ///
    /// Precondition: |M̂|_dual = 1.
    DualQuaternion inverse() const {
        return conjugate();
    }

    /// SE(3) retraction (REQ-EF-15). Renormalizes the real part to unit
    /// magnitude and projects out the q_r-parallel component of q_d so the
    /// SE(3) constraint q_r · q_d = 0 (4-vector dot) holds. Used after a
    /// numerical integration step to project drift back onto the manifold.
    DualQuaternion normalized() const {
        TrackedValue<T> mr = real.magnitude();
        Quaternion<T> qr_n(real.w/mr, real.x/mr, real.y/mr, real.z/mr);
        TrackedValue<T> dual_dot =
            qr_n.w * dual.w + qr_n.x * dual.x + qr_n.y * dual.y + qr_n.z * dual.z;
        Quaternion<T> qd_n = dual - dual_dot * qr_n;
        return DualQuaternion(qr_n, qd_n);
    }

    // --- Action ---

    /// Apply this (assumed unit) pose dual quaternion to a 3-point:
    ///   p' = R p + t.
    ///
    /// Algebraic form: lift p to P̂ = 1 + ε p_pure, compute M̂ P̂ M̂♯, take
    /// the vector part of the dual component. The dual part evaluates to
    ///   q_r p q_r* + (q_d q_r* − q_r q_d*)
    /// which equals R(p) + t for a unit pose.
    ///
    /// Precondition: M̂ is unit (a pose). Result is meaningless for
    /// non-unit M̂.
    Vector3<T> apply(const Vector3<T>& p) const {
        Quaternion<T> p_pure = Quaternion<T>::pure(p);
        Quaternion<T> qr_conj = real.conjugate();
        Quaternion<T> qd_conj = dual.conjugate();
        Quaternion<T> dual_out =
            real * p_pure * qr_conj + (dual * qr_conj - real * qd_conj);
        return dual_out.vector();
    }

    /// Apply only the rotation component to a direction 3-vector (no
    /// translation): equivalent to `rotation().rotate(v)`. Precondition:
    /// M̂ is unit.
    Vector3<T> apply_direction(const Vector3<T>& v) const {
        return real.rotate(v);
    }

    // --- Lie algebra: screw exp, log, Sclerp ---

    /// Screw exponential.
    ///
    /// For pure dual quaternion Û = u + ε v (with u, v as 3-vectors lifted
    /// to pure quaternions), returns the unit dual quaternion exp(Û).
    /// Component form (REQ-DQ-18):
    ///   q_r = cos|u| + sinc|u| · u                  (= Quaternion::exp_pure(u))
    ///   q_d.scalar() = −(u · v) sinc|u|
    ///   q_d.vector() = sinc|u| · v + (u · v) · β · u
    /// with β = (cos|u| − sinc|u|) / |u|².
    ///
    /// Inverse of `log_screw` on the half-angle ball (|u| ≤ π/2).
    static DualQuaternion exp_screw(const Vector3<T>& u, const Vector3<T>& v) {
        TrackedValue<T> theta_sq = u.x*u.x + u.y*u.y + u.z*u.z;
        TrackedValue<T> theta = sqrt(theta_sq);
        TrackedValue<T> uv_dot = u.x*v.x + u.y*v.y + u.z*v.z;
        TrackedValue<T> c = cos(theta);
        TrackedValue<T> sinc_t = taylor_sinc(theta, theta_sq);
        TrackedValue<T> beta = taylor_cos_minus_sinc_over_theta_sq(theta, theta_sq);

        // Real part: standard quaternion exponential of the pure quaternion u.
        Quaternion<T> q_r(c, sinc_t * u.x, sinc_t * u.y, sinc_t * u.z);

        // Dual part: scalar = −(u·v) sinc|u|;
        //            vector = sinc|u| · v + (u·v) · β · u.
        TrackedValue<T> qd_w = -(uv_dot * sinc_t);
        TrackedValue<T> uv_beta = uv_dot * beta;
        Quaternion<T> q_d(
            qd_w,
            sinc_t * v.x + uv_beta * u.x,
            sinc_t * v.y + uv_beta * u.y,
            sinc_t * v.z + uv_beta * u.z
        );

        return DualQuaternion(q_r, q_d);
    }

    /// Screw exponential of a pure dual quaternion. Scalar parts of the
    /// input are ignored; the vector parts `twist.angular()` and
    /// `twist.linear()` carry the screw.
    static DualQuaternion exp_screw(const DualQuaternion& twist) {
        return exp_screw(twist.angular(), twist.linear());
    }

    /// Screw logarithm. Given this unit pose M̂, return the pure dual
    /// quaternion Û = u + ε v such that `exp_screw(Û) == M̂` (on the
    /// half-angle ball).
    ///
    /// Derivation. The real part inverts to u = q_r.log_unit() (the
    /// half-angle rotation vector). The dual part inverts to
    ///   s := u · v = −q_d.scalar() / sinc|u|
    ///   v        = (q_d.vector() − s · β · u) / sinc|u|
    /// with β = (cos|u| − sinc|u|) / |u|².
    ///
    /// Precondition. M̂ is unit (a pose). The shortest-path convention on
    /// `Quaternion::log_unit` keeps |u| ≤ π/2, so sinc|u| ≥ sinc(π/2) and
    /// the divisions above are bounded away from zero.
    DualQuaternion log_screw() const {
        Vector3<T> u = real.log_unit();
        TrackedValue<T> theta_sq = u.x*u.x + u.y*u.y + u.z*u.z;
        TrackedValue<T> theta = sqrt(theta_sq);
        TrackedValue<T> sinc_t = taylor_sinc(theta, theta_sq);
        TrackedValue<T> beta = taylor_cos_minus_sinc_over_theta_sq(theta, theta_sq);

        // Invert q_d.scalar() = −s · sinc|u| for s.
        TrackedValue<T> s = -dual.w / sinc_t;

        // Invert q_d.vector() = sinc|u| · v + s · β · u for v.
        TrackedValue<T> sb = s * beta;
        Vector3<T> v(
            (dual.x - sb * u.x) / sinc_t,
            (dual.y - sb * u.y) / sinc_t,
            (dual.z - sb * u.z) / sinc_t
        );

        return DualQuaternion(Quaternion<T>::pure(u), Quaternion<T>::pure(v));
    }

    /// Screw linear interpolation between two unit pose dual quaternions.
    ///
    /// `sclerp(m0, m1, t)` traverses the screw geodesic from m0 (t = 0) to
    /// m1 (t = 1) at constant screw rate. Algorithm:
    ///
    ///   delta = m0⁻¹ · m1
    ///   step  = exp_screw(t · log_screw(delta))
    ///   result = m0 · step
    ///
    /// Precondition. Both m0 and m1 are unit (poses).
    static DualQuaternion sclerp(const DualQuaternion& m0,
                                 const DualQuaternion& m1,
                                 const TrackedValue<T>& t) {
        DualQuaternion delta = m0.inverse() * m1;
        DualQuaternion log_delta = delta.log_screw();
        Vector3<T> u_scaled = t * log_delta.angular();
        Vector3<T> v_scaled = t * log_delta.linear();
        return m0 * exp_screw(u_scaled, v_scaled);
    }
};

} // namespace math
