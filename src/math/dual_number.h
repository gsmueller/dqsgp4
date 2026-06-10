#pragma once

/// @file dual_number.h
/// Dual numbers a + ε b with ε² = 0, error-tracked.
///
/// Dual numbers extend the reals by adjoining a nilpotent ε satisfying ε² = 0.
/// Two uses in this library:
///
///   1. Forward-mode automatic differentiation. f(a + εb) = f(a) + ε b f'(a)
///      exactly (no truncation), since ε² = 0 truncates every power series at
///      the linear term in ε. Setting b = 1 gives both f(a) and f'(a) with no
///      finite-difference error.
///
///   2. SE(3) screw-motion calculus. A unit dual quaternion encodes a
///      rigid-body pose; sin and cos of a dual angle θ̂ = θ + ε d give the
///      screw exponential in closed form:
///        sin(θ + εd) = sin θ + ε d cos θ
///        cos(θ + εd) = cos θ − ε d sin θ.
///
/// Every operation in this header is exact in closed form — ε² = 0 truncates
/// the relevant Taylor series at the linear ε term, so no series error is
/// introduced. The result's three-error budget is propagated entirely through
/// the underlying TrackedValue<T> operators (REQ-EF-3).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-13, AUD-CC-15, AUD-CC-16, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-2, AUD-EF-7.

#include "tracked_value.h"

namespace math {

/// Dual number a + ε b with ε² = 0. Both components are error-tracked.
///
/// @tparam T  Underlying numeric type (double, cpp_bin_float_50, etc.).
template<typename T>
struct DualNumber {
    TrackedValue<T> real;  ///< Real part a.
    TrackedValue<T> dual;  ///< Coefficient of ε.

    // --- Constructors ---

    /// Zero dual number (real = dual = 0).
    DualNumber() : real(), dual() {}

    /// Promote a real number to a dual with zero dual part.
    DualNumber(const TrackedValue<T>& r) : real(r), dual() {}

    /// Construct from explicit real and dual parts.
    DualNumber(const TrackedValue<T>& r, const TrackedValue<T>& d)
        : real(r), dual(d) {}

    // --- Named constructors ---

    /// Additive identity: 0 + ε · 0.
    static DualNumber zero() { return DualNumber(); }

    /// Multiplicative identity: 1 + ε · 0.
    static DualNumber identity() {
        return DualNumber(exact<T>(1), exact<T>(0));
    }

    /// The nilpotent unit ε itself: 0 + ε · 1.
    static DualNumber epsilon() {
        return DualNumber(exact<T>(0), exact<T>(1));
    }

    // --- Arithmetic ---

    /// (a + εb) + (c + εd) = (a + c) + ε(b + d).
    friend DualNumber operator+(const DualNumber& a, const DualNumber& b) {
        return DualNumber(a.real + b.real, a.dual + b.dual);
    }

    /// (a + εb) − (c + εd) = (a − c) + ε(b − d).
    friend DualNumber operator-(const DualNumber& a, const DualNumber& b) {
        return DualNumber(a.real - b.real, a.dual - b.dual);
    }

    /// Unary negation.
    DualNumber operator-() const {
        return DualNumber(-real, -dual);
    }

    /// Multiplication: (a + εb)(c + εd) = ac + ε(ad + bc), since ε² = 0.
    friend DualNumber operator*(const DualNumber& a, const DualNumber& b) {
        return DualNumber(
            a.real * b.real,
            a.real * b.dual + a.dual * b.real
        );
    }

    /// Scalar multiplication (scalar on left).
    friend DualNumber operator*(const TrackedValue<T>& s, const DualNumber& a) {
        return DualNumber(s * a.real, s * a.dual);
    }

    /// Scalar multiplication (scalar on right).
    friend DualNumber operator*(const DualNumber& a, const TrackedValue<T>& s) {
        return DualNumber(a.real * s, a.dual * s);
    }

    /// Division.
    ///
    /// Derivation. 1/(c + εd) is the Taylor expansion of 1/x at x = c
    /// evaluated at the infinitesimal shift εd. Because ε² = 0 the series
    /// terminates at the linear ε term:
    ///   1/(c + εd) = 1/c − ε d/c².
    /// Multiplying by (a + εb) and dropping ε²:
    ///   (a + εb) / (c + εd) = a/c + ε (b c − a d)/c².
    /// Error in each category propagates through the underlying TrackedValue<T>
    /// division and multiplication (REQ-EF-3).
    friend DualNumber operator/(const DualNumber& a, const DualNumber& b) {
        TrackedValue<T> c2 = b.real * b.real;
        return DualNumber(
            a.real / b.real,
            (a.dual * b.real - a.real * b.dual) / c2
        );
    }

    // --- Elementary functions (closed form via ε² = 0) ---

    /// sqrt(a + εb) = √a + ε b/(2√a).
    /// f(x) = √x has f'(x) = 1/(2√x); applying f(a + εb) = f(a) + ε b f'(a).
    friend DualNumber sqrt(const DualNumber& a) {
        TrackedValue<T> sr = sqrt(a.real);
        return DualNumber(sr, a.dual / (exact<T>(2) * sr));
    }

    /// sin(a + εb) = sin a + ε b cos a. f'(x) = cos x.
    friend DualNumber sin(const DualNumber& a) {
        return DualNumber(sin(a.real), a.dual * cos(a.real));
    }

    /// cos(a + εb) = cos a − ε b sin a. f'(x) = −sin x.
    friend DualNumber cos(const DualNumber& a) {
        return DualNumber(cos(a.real), -(a.dual * sin(a.real)));
    }
};

} // namespace math
