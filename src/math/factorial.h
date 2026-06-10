#pragma once

// Exact factorial, double factorial, falling factorial, rising factorial.
//
// These return TrackedValue with:
//   measurement = 0 (mathematical, not physical)
//   precision = 0 when inputs are exact integers (the result is an exact integer)
//   accuracy = 0 (no model)
//
// When inputs are TrackedValue (non-integer alpha for generalized binomial),
// errors propagate through the multiplications.

#include "tracked_value.h"

namespace math {

/// n! as exact integer TrackedValue
template<typename T>
TrackedValue<T> factorial(int n) {
    TrackedValue<T> result = exact<T>(1);
    for (int i = 2; i <= n; ++i) {
        result = result * exact<T>(i);
    }
    return result;
}

/// n!! (double factorial) as exact integer TrackedValue
/// n!! = n * (n-2) * (n-4) * ... * (2 or 1)
template<typename T>
TrackedValue<T> double_factorial(int n) {
    if (n <= 0) return exact<T>(1);
    TrackedValue<T> result = exact<T>(1);
    for (int i = n; i >= 2; i -= 2) {
        result = result * exact<T>(i);
    }
    return result;
}

/// Falling factorial: alpha * (alpha-1) * (alpha-2) * ... * (alpha-k+1)
/// For integer alpha and k <= alpha, the result is exact: alpha! / (alpha-k)!
/// For TrackedValue alpha (non-integer), errors propagate through k multiplications.
template<typename T>
TrackedValue<T> falling_factorial(const TrackedValue<T>& alpha, int k) {
    if (k <= 0) return exact<T>(1);
    TrackedValue<T> result = alpha;
    for (int j = 1; j < k; ++j) {
        result = result * (alpha - exact<T>(j));
    }
    return result;
}

/// Rising factorial (Pochhammer): alpha * (alpha+1) * (alpha+2) * ... * (alpha+k-1)
template<typename T>
TrackedValue<T> rising_factorial(const TrackedValue<T>& alpha, int k) {
    if (k <= 0) return exact<T>(1);
    TrackedValue<T> result = alpha;
    for (int j = 1; j < k; ++j) {
        result = result * (alpha + exact<T>(j));
    }
    return result;
}

/// Generalized binomial coefficient: C(alpha, k) = falling_factorial(alpha, k) / k!
template<typename T>
TrackedValue<T> generalized_binomial(const TrackedValue<T>& alpha, int k) {
    if (k == 0) return exact<T>(1);
    return falling_factorial(alpha, k) / factorial<T>(k);
}

} // namespace math
