#pragma once

/**
 * @file matrix3.h
 * @brief 3×3 matrix of TrackedValue components with full error propagation.
 *
 * Provides the general (off-diagonal) tensor infrastructure issue H1 needs for
 * a rigid-body inertia tensor with products of inertia and its inverse — the
 * capability the diagonal-only dynamics/inertia.h lacked. Every operation
 * (matrix-vector product, determinant, inverse) propagates all three error
 * components through the TrackedValue arithmetic.
 */

#include "tracked_value.h"
#include "vector3.h"

namespace math {

/// 3×3 matrix, row-major `m[row][col]`, of TrackedValue components.
template<typename T>
struct Matrix3 {
    TrackedValue<T> m[3][3];  ///< entries; default-constructed TrackedValue == 0

    /// Zero matrix (every entry the exact-zero TrackedValue).
    Matrix3() {}

    /// Diagonal matrix from a vector of diagonal entries.
    static Matrix3 diagonal(const Vector3<T>& d) {
        Matrix3 M;
        M.m[0][0] = d.x;
        M.m[1][1] = d.y;
        M.m[2][2] = d.z;
        return M;
    }

    /// Symmetric matrix (e.g. an inertia tensor) from the three principal
    /// moments and the three products of inertia.
    static Matrix3 symmetric(const TrackedValue<T>& Ixx, const TrackedValue<T>& Iyy,
                             const TrackedValue<T>& Izz, const TrackedValue<T>& Ixy,
                             const TrackedValue<T>& Ixz, const TrackedValue<T>& Iyz) {
        Matrix3 M;
        M.m[0][0] = Ixx; M.m[1][1] = Iyy; M.m[2][2] = Izz;
        M.m[0][1] = Ixy; M.m[1][0] = Ixy;
        M.m[0][2] = Ixz; M.m[2][0] = Ixz;
        M.m[1][2] = Iyz; M.m[2][1] = Iyz;
        return M;
    }

    /// Matrix-vector product M·v.
    friend Vector3<T> operator*(const Matrix3& M, const Vector3<T>& v) {
        return Vector3<T>(
            M.m[0][0] * v.x + M.m[0][1] * v.y + M.m[0][2] * v.z,
            M.m[1][0] * v.x + M.m[1][1] * v.y + M.m[1][2] * v.z,
            M.m[2][0] * v.x + M.m[2][1] * v.y + M.m[2][2] * v.z);
    }

    /// Matrix-matrix product A·B (row i of A dotted with column j of B). Composes
    /// frame rotations; (A·B)·v = A·(B·v), so right-most acts first on a vector.
    friend Matrix3 operator*(const Matrix3& A, const Matrix3& B) {
        Matrix3 M;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                M.m[i][j] = A.m[i][0] * B.m[0][j]
                          + A.m[i][1] * B.m[1][j]
                          + A.m[i][2] * B.m[2][j];
            }
        }
        return M;
    }

    /// Transpose.
    Matrix3 transpose() const {
        Matrix3 M;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) M.m[i][j] = m[j][i];
        return M;
    }

    /// Determinant (cofactor expansion along row 0).
    TrackedValue<T> determinant() const {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    /// Inverse via the adjugate / determinant. The TrackedValue division
    /// carries its own minimum-denominator degenerate guard (a singular matrix
    /// flags a max-bound error rather than producing a silent NaN).
    Matrix3 inverse() const {
        TrackedValue<T> det = determinant();
        Matrix3 inv;
        inv.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
        inv.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
        inv.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
        inv.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
        inv.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
        inv.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
        inv.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
        inv.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
        inv.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
        return inv;
    }
};

} // namespace math
