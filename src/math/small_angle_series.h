#pragma once

/// @file small_angle_series.h
/// Small-argument Taylor expansions for trigonometric helpers, with
/// truncation bounds added to errors.precision per REQ-EF-6 / AUD-EF-5.
///
/// These helpers are shared between Quaternion<T> (exp_pure, log_unit) and
/// DualQuaternion<T> (exp_screw, log_screw, sclerp). Centralizing them
/// keeps the Taylor-branch logic in one place and preserves "one voice"
/// across modules (CON-13, AUD-CC-13).
///
/// Every helper:
///   - accepts theta (and where useful theta_sq) to avoid recomputation;
///   - returns a TrackedValue<T> whose precision component includes the
///     rigorous Taylor-truncation bound when the small-argument branch
///     fires (AUD-EF-5);
///   - otherwise falls through to the closed-form expression with normal
///     per-category propagation (REQ-EF-3).
///
/// Threshold (T-dependent, B3.2). Each helper branches at |argument| <
/// taylor_branch_threshold<T>(C), the point where the next-term Taylor
/// truncation balances the closed form's representation/cancellation error.
/// Below it the Taylor branch is the more accurate of the two; above it the
/// closed form is. The crossover is tau = (C·eps)^(1/6) (derived per helper
/// below). A wider T (smaller eps) SHRINKS tau, so the closed form is used over
/// more of the domain exactly where its extra digits make it the better choice
/// — the fixed 1e-4 was T-independent and, for wide T, kept the Taylor branch
/// (with its fixed-order truncation) live well past the point the closed form
/// had become more accurate.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-8, AUD-CC-9,
///   AUD-CC-10, AUD-CC-13, AUD-CC-15, AUD-CC-16, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-2, AUD-EF-5.

#include "tracked_value.h"

#include <cmath>
#include <limits>

namespace math {

/// T-dependent small-argument branch threshold (B3.2). Returns tau =
/// (balance_constant · eps)^(1/6), where eps = machine epsilon of T. The
/// exponent 1/6 is common to all three helpers: each balances a next-Taylor-
/// term truncation against the closed form's error (a representation floor ~eps
/// for sinc / half-angle, or a ~eps/theta^2 cancellation for cos-minus-sinc),
/// and both forms of that balance reduce to theta^6 = C·eps. A larger eps gives
/// a larger tau (use Taylor more — double); a tiny eps gives a tiny tau (use the
/// closed form more — wide T), which is the correctness fix over the old 1e-4.
template<typename T>
T taylor_branch_threshold(const T& balance_constant) {
    using std::pow;
    T eps = std::numeric_limits<T>::epsilon();
    return pow(balance_constant * eps, T(1) / T(6));
}

/// sinc(theta) = sin(theta) / theta with a small-argument Taylor branch
/// to avoid 0/0 near theta = 0.
///
/// Branch (|theta| below the B3.2 T-dependent threshold (5040·eps)^(1/6)):
///   sinc(theta) = 1 − theta²/6 + theta⁴/120 + R(theta),
///   |R(theta)| ≤ |theta|⁶ / 5040  (magnitude of the next alternating
///                                  Taylor term).
/// The bound |theta_sq|³ / 5040 is added to result.errors.precision.
template<typename T>
TrackedValue<T> taylor_sinc(const TrackedValue<T>& theta,
                            const TrackedValue<T>& theta_sq) {
    using std::abs;
    // Next dropped term is theta^6/5040; balanced against the ~eps closed-form
    // floor this gives tau = (5040·eps)^(1/6) (B3.2).
    const T threshold = taylor_branch_threshold<T>(T(5040));
    if (abs(theta.value) < threshold) {
        TrackedValue<T> t4 = theta_sq * theta_sq;
        TrackedValue<T> result =
            exact<T>(1) - theta_sq / exact<T>(6) + t4 / exact<T>(120);
        T ts_abs = abs(theta_sq.value);
        T trunc_bound = ts_abs * ts_abs * ts_abs / T(5040);
        result.errors.precision = result.errors.precision + trunc_bound;
        return result;
    }
    return sin(theta) / theta;
}

/// Half-angle scale factor for Quaternion::log_unit:
///   scale = atan2(qv_norm, w_pos) / qv_norm
/// with a small-qv_norm Taylor branch.
///
/// For a unit quaternion with w_pos ≥ 0, qv_norm = sin(angle/2). The scale
/// reduces to arcsin(qv_norm) / qv_norm. Its Taylor expansion is
///   arcsin(s)/s = 1 + s²/6 + 3 s⁴/40 + 5 s⁶/112 + …
/// The branch truncates after the s⁴ term and adds 5 |s_sq|³ / 112 to
/// result.errors.precision (magnitude of the next term).
///
/// Precondition. w_pos ≥ 0 (caller has flipped q → −q if necessary).
template<typename T>
TrackedValue<T> taylor_half_angle_scale(const TrackedValue<T>& qv_norm,
                                        const TrackedValue<T>& qv_norm_sq,
                                        const TrackedValue<T>& w_pos) {
    using std::abs;
    // Next dropped term is 5·s^6/112; balanced against the ~eps closed-form
    // floor this gives tau = ((112/5)·eps)^(1/6) (B3.2).
    const T threshold = taylor_branch_threshold<T>(T(112) / T(5));
    if (abs(qv_norm.value) < threshold) {
        TrackedValue<T> s2 = qv_norm_sq;
        TrackedValue<T> s4 = s2 * s2;
        TrackedValue<T> result =
            exact<T>(1) + s2 / exact<T>(6) + exact<T>(3) * s4 / exact<T>(40);
        T s2_abs = abs(qv_norm_sq.value);
        // R11 / Abramowitz-Stegun §15.1.10: arcsin(s)/s has POSITIVE
        // coefficients, so the next-term magnitude alone is not the
        // rigorous Leibniz bound. The geometric correction factor
        // 1/(1 - 2 s²) accounts for the tail of the positive series
        // (successive-term ratio ≤ 2 s² for |s| < 1/√2).
        T leading = T(5) * s2_abs * s2_abs * s2_abs / T(112);
        T correction_denom = T(1) - T(2) * s2_abs;  // > 0 for s² < 0.5
        T trunc_bound = leading / correction_denom;
        result.errors.precision = result.errors.precision + trunc_bound;
        return result;
    }
    return atan2(qv_norm, w_pos) / qv_norm;
}

/// (cos(theta) − sinc(theta)) / theta², used by DualQuaternion::exp_screw
/// and DualQuaternion::log_screw to evaluate the coefficient β of the
/// u-direction in the screw dual-part vector.
///
/// Branch (|theta| below the B3.2 T-dependent threshold (840·eps)^(1/6)):
///   (cos θ − sinc θ) / θ² = −1/3 + θ²/30 − θ⁴/840 + …
/// The branch truncates after the θ²/30 term and adds |theta_sq|² / 840 to
/// result.errors.precision (magnitude of the next term).
///
/// Outside the branch, the closed-form expression
///   (cos(theta) − sin(theta)/theta) / theta_sq
/// is used directly. The subtraction has bounded condition number above the
/// threshold because cos and sinc separate noticeably past that point (the
/// threshold is exactly where the ~eps/θ² cancellation stops dominating).
template<typename T>
TrackedValue<T> taylor_cos_minus_sinc_over_theta_sq(
    const TrackedValue<T>& theta,
    const TrackedValue<T>& theta_sq) {
    using std::abs;
    // Next dropped term is theta^4/840; balanced against the ~eps/theta^2
    // closed-form cancellation this gives tau = (840·eps)^(1/6) (B3.2).
    const T threshold = taylor_branch_threshold<T>(T(840));
    if (abs(theta.value) < threshold) {
        TrackedValue<T> result = ratio<T>(-1, 3) + theta_sq / exact<T>(30);
        T ts_abs = abs(theta_sq.value);
        T trunc_bound = ts_abs * ts_abs / T(840);
        result.errors.precision = result.errors.precision + trunc_bound;
        return result;
    }
    return (cos(theta) - sin(theta) / theta) / theta_sq;
}

} // namespace math
