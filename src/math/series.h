#pragma once

// Series evaluation with rigorous error bounds.
//
// Evaluates truncated series to a caller-specified tolerance,
// returning the value with the truncation error added to the
// precision bound of the result.

#include "tracked_value.h"
#include <functional>

namespace math {

/// Evaluate an alternating series: S = sum_{k=start}^{inf} term(k)
/// where terms alternate in sign and decrease in absolute value.
///
/// Stops when the truncation bound drops below tolerance.
/// Returns TrackedValue with truncation error added to precision bound.
///
/// Two bound modes:
///   convergence_ratio < 0 (default): Leibniz bound |term(N)|.
///     Rigorous, requires no knowledge of the ratio. Conservative.
///   convergence_ratio >= 0: geometric tail bound |term(N)| * r/(1-r).
///     Much tighter when r is small (e.g., r = e'^2 ~ 0.007 gives
///     a bound ~150x tighter than Leibniz). See Ch 5, §5.9.4.
template<typename T>
TrackedValue<T> alternating_series(
    int start,
    // term(k) returns the k-th term as a TrackedValue.
    // The term function must produce alternating signs and decreasing |value|.
    std::function<TrackedValue<T>(int)> term,
    const T& tolerance,
    int max_terms = 10000,
    const T& convergence_ratio = T(-1))
{
    TrackedValue<T> sum = exact<T>(0);
    TrackedValue<T> last_term;
    bool use_geometric = (convergence_ratio >= T(0) && convergence_ratio < T(1));

    for (int k = start; k < start + max_terms; ++k) {
        last_term = term(k);
        sum = sum + last_term;

        using std::abs;
        T term_mag = abs(last_term.value);
        T truncation;
        if (use_geometric) {
            truncation = term_mag * convergence_ratio / (T(1) - convergence_ratio);
        } else {
            truncation = term_mag;
        }

        if (truncation < tolerance) {
            sum.errors.precision = sum.errors.precision + truncation;
            return sum;
        }
    }

    // Did not converge within max_terms
    using std::abs;
    T term_mag = abs(last_term.value);
    if (use_geometric) {
        sum.errors.precision = sum.errors.precision +
            term_mag * convergence_ratio / (T(1) - convergence_ratio);
    } else {
        sum.errors.precision = sum.errors.precision + term_mag;
    }
    return sum;
}

/// Evaluate a general series with geometric-type convergence.
/// For non-alternating series where we can bound the tail by
/// |tail| <= |last_term| * ratio / (1 - ratio).
///
/// The caller provides an estimate of the convergence ratio.
template<typename T>
TrackedValue<T> geometric_series(
    int start,
    std::function<TrackedValue<T>(int)> term,
    const T& tolerance,
    const T& convergence_ratio,  // upper bound on |term(k+1)/term(k)|
    int max_terms = 10000)
{
    TrackedValue<T> sum = exact<T>(0);
    TrackedValue<T> last_term;

    for (int k = start; k < start + max_terms; ++k) {
        last_term = term(k);
        sum = sum + last_term;

        using std::abs;
        T tail_bound = abs(last_term.value) * convergence_ratio / (T(1) - convergence_ratio);
        if (tail_bound < tolerance) {
            sum.errors.precision = sum.errors.precision + tail_bound;
            return sum;
        }
    }

    using std::abs;
    sum.errors.precision = sum.errors.precision + abs(last_term.value);
    return sum;
}

/// Newton iteration for sqrt(I) where I is close to 1.
///
/// Computes s such that s^2 = I, starting from s_0 = 1.
/// Iteration: s_{j+1} = (s_j + I/s_j) / 2
///
/// Quadratic convergence: |e_j| <= |U|^{2^j} / 2^{j+1}
/// where U = I.value - 1. For |U| ~ 0.007 (WGS84 geodetic series),
/// 4 iterations reach delta_p < 1e-40. See Ch 5, Lemma 5.6.1.
///
/// The result inherits the three-error budget of I through the
/// TrackedValue arithmetic in each iteration step.
template<typename T>
TrackedValue<T> series_sqrt(
    const TrackedValue<T>& I,
    const T& tolerance,
    int max_iterations = 60)
{
    TrackedValue<T> s = exact<T>(1);

    for (int j = 0; j < max_iterations; ++j) {
        TrackedValue<T> s_new = (s + I / s) / exact<T>(2);

        using std::abs;
        T correction = abs((s_new - s).value);
        s = s_new;

        if (correction < tolerance) {
            s.errors.precision = s.errors.precision + correction;
            return s;
        }
    }

    // Did not converge — add last correction as precision bound
    using std::abs;
    T correction = abs((I / s - s).value) / T(2);
    s.errors.precision = s.errors.precision + correction;
    return s;
}

} // namespace math
