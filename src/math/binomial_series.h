#pragma once

// Generalized binomial series (1+x)^alpha with error-bounded evaluation.
//
// Provides two levels of abstraction:
//
// 1. make_binomial_evaluator — returns a closure that evaluates (1+x)^alpha
//    for any x with |x| < 1. The closure captures alpha and the tolerance,
//    and computes exact rational binomial coefficients C(alpha, k) on the fly.
//    This is the general-purpose building block.
//
// 2. geodetic_binomial_coefficient — computes the k-th coefficient of an
//    integrated geodetic binomial series (Stage 1 x Stage 2 of the
//    three-stage pattern in Ch 5, Theorem 5.6.2). This is the specialized
//    function for ellipsoidal geometry.
//
// All coefficients are exact rationals (Ch 4, Theorem 4.10.2) — no
// floating-point literals. Rounding occurs only in the final multiplication
// by powers of the argument.
//
// References:
//   Ch 4, Theorem 4.10.1 — binomial series convergence and remainder
//   Ch 4, Theorem 4.10.2 — exact rational coefficients for alpha = p/q
//   Ch 5, §5.6 — geodetic binomial series and three-stage derivation
//   Ch 5, Lemma 5.6.1 — series_sqrt (in series.h)

#include "tracked_value.h"
#include "factorial.h"
#include "wallis.h"
#include "series.h"
#include <functional>
#include <vector>

namespace math {

/// Evaluate (1+x)^alpha to the specified tolerance using the binomial series
/// with exact rational coefficients and geometric tail bound.
///
/// Returns a closure: given x (as TrackedValue), produces (1+x)^alpha.
/// The closure captures alpha and tolerance at initialization time.
/// At each call, it evaluates the series to the number of terms needed
/// for the given |x|, using the geometric tail bound (Ch 1, Thm 1.5.3):
///   |R_N| <= |C(alpha, N+1)| * |x|^{N+1} / (1 - |x|)
///
/// This is the init-time closure pattern from Ch 3, §3.4:
/// alpha and tolerance are fixed; x varies per call.
template<typename T>
std::function<TrackedValue<T>(const TrackedValue<T>&)>
make_binomial_evaluator(
    const TrackedValue<T>& alpha,
    const T& tolerance)
{
    // Capture alpha and tolerance in the closure.
    // Binomial coefficients are computed per call because the number
    // of terms needed depends on |x|.
    return [alpha, tolerance](const TrackedValue<T>& x) -> TrackedValue<T> {
        using std::abs;

        T abs_x = abs(x.value);

        // For |x| >= 1, the binomial series does not converge.
        // Return a value with large precision error to signal the failure.
        if (abs_x >= T(1)) {
            TrackedValue<T> result(T(0), T(0), abs_x, T(0));
            return result;
        }

        // Evaluate sum_{k=0}^{N} C(alpha, k) * x^k via Horner-like accumulation.
        // Coefficients C(alpha, k) are exact rationals computed incrementally:
        //   C(alpha, 0) = 1
        //   C(alpha, k) = C(alpha, k-1) * (alpha - k + 1) / k
        TrackedValue<T> sum = exact<T>(1);       // k=0 term
        TrackedValue<T> x_power = exact<T>(1);   // x^k
        TrackedValue<T> coeff = exact<T>(1);      // C(alpha, k)

        for (int k = 1; k < 10000; ++k) {
            // Incremental coefficient: C(alpha, k) = C(alpha, k-1) * (alpha-k+1)/k
            coeff = coeff * (alpha - exact<T>(k - 1)) / exact<T>(k);
            x_power = x_power * x;

            TrackedValue<T> term = coeff * x_power;
            sum = sum + term;

            // Geometric tail bound: |C(alpha, N+1)| * |x|^{N+1} / (1 - |x|)
            T tail_bound = abs(term.value) * abs_x / (T(1) - abs_x);
            if (tail_bound < tolerance) {
                sum.errors.precision = sum.errors.precision + tail_bound;
                return sum;
            }
        }

        // Did not converge — add conservative bound
        sum.errors.precision = sum.errors.precision + abs(x_power.value);
        return sum;
    };
}

/// Compute the k-th coefficient of an integrated geodetic binomial series.
///
/// This implements the three-stage coefficient derivation (Ch 5, Thm 5.6.2):
///   Stage 1: Generalized binomial coefficient C(alpha, k) * (-1)^k
///   Stage 2: Wallis integration factor
///     - cosine_weight = false (w=1): even Wallis ratio (2k-1)!!/(2k)!!
///     - cosine_weight = true  (w=cos phi): simple integral 1/(2k+1)
///
/// The result is an exact rational number (no floating-point error in the
/// coefficient itself). The only precision error comes from the representation
/// of the rational as type T.
///
/// Usage:
///   TrackedValue<T> c_k = geodetic_binomial_coefficient<T>(ratio<T>(-3, 2), k, false);
///   // Returns the k-th meridian arc coefficient (alpha = -3/2, w = 1)
template<typename T>
TrackedValue<T> geodetic_binomial_coefficient(
    const TrackedValue<T>& alpha,
    int k,
    bool cosine_weight = false)
{
    // Stage 1: C(alpha, k)
    // The (-1)^k factor from the argument (-e^2 sin^2 phi)^k is NOT
    // included here — it is part of the series evaluation, not the
    // coefficient. The caller multiplies by (-e^2)^k or (e^2)^k as
    // appropriate for the specific series.
    TrackedValue<T> binom = generalized_binomial(alpha, k);

    // Stage 2: Wallis integration factor
    if (cosine_weight) {
        // integral of sin^{2k}(phi) * cos(phi) dphi = 1/(2k+1)
        return binom * sin_power_cos_integral<T>(k);
    } else {
        // integral of sin^{2k}(phi) dphi = W_{2k}
        // Normalized: W_{2k} / (pi/2) = (2k-1)!! / (2k)!!
        if (k == 0) return binom;  // W_0 / (pi/2) = 1
        return binom * double_factorial<T>(2 * k - 1) / double_factorial<T>(2 * k);
    }
}

} // namespace math
