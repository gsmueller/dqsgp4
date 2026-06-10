#pragma once

/**
 * @file vector3.h
 * @brief 3D vector with TrackedValue components and full error propagation.
 *
 * Every arithmetic operation (add, subtract, scale, dot, cross, magnitude)
 * propagates all three error components (measurement, precision, accuracy)
 * through the TrackedValue arithmetic operators.
 */

#include "tracked_value.h"

namespace math {

/**
 * @brief Three-dimensional vector of TrackedValue components.
 *
 * Used for position [km] and velocity [km/s] in the TEME frame.
 * Each component carries its own three-error budget. Vector operations
 * compose the component errors via the TrackedValue operator overloads.
 *
 * @tparam T  The underlying numeric type (double, cpp_bin_float_50, etc.)
 */
template<typename T>
struct Vector3 {
    TrackedValue<T> x, y, z;  ///< components, each with its own budget

    /// Zero vector.
    Vector3() : x(), y(), z() {}
    /// From explicit (x, y, z) components.
    Vector3(const TrackedValue<T>& x_, const TrackedValue<T>& y_, const TrackedValue<T>& z_)
        : x(x_), y(y_), z(z_) {}

    /// Component-wise addition.
    friend Vector3 operator+(const Vector3& a, const Vector3& b) {
        return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    /// Component-wise subtraction.
    friend Vector3 operator-(const Vector3& a, const Vector3& b) {
        return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    /// Scalar-vector multiplication (scalar on left).
    friend Vector3 operator*(const TrackedValue<T>& s, const Vector3& v) {
        return Vector3(s * v.x, s * v.y, s * v.z);
    }

    /// Scalar-vector multiplication (scalar on right).
    friend Vector3 operator*(const Vector3& v, const TrackedValue<T>& s) {
        return Vector3(v.x * s, v.y * s, v.z * s);
    }

    /// Dot product: returns a scalar TrackedValue with propagated errors.
    TrackedValue<T> dot(const Vector3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    /// Cross product: returns a Vector3 with propagated errors.
    Vector3 cross(const Vector3& other) const {
        return Vector3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    /// Euclidean magnitude: √(x² + y² + z²) with propagated errors.
    TrackedValue<T> magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }

    /// Scalar division (each component ÷ scalar), full error propagation via
    /// the TrackedValue division operator (which carries its own minimum-
    /// denominator degenerate guard).
    friend Vector3 operator/(const Vector3& v, const TrackedValue<T>& s) {
        return Vector3(v.x / s, v.y / s, v.z / s);
    }

    /// Unit vector v / |v|, with propagated errors. When |v| is within its own
    /// error of zero the direction is undetermined; the TrackedValue division's
    /// degenerate guard (denominator magnitude ≤ error ⇒ max bound) flags it.
    Vector3 normalize() const {
        return *this / magnitude();
    }
};

} // namespace math
