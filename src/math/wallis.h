#pragma once

// Wallis integrals: W_n = integral from 0 to pi/2 of cos^n(phi) dphi
//
// These are exact rational numbers computed from double factorials.
// W_{2k}   = (2k-1)!! / (2k)!! * pi/2
// W_{2k+1} = (2k)!! / (2k+1)!!
//
// For the geodetic series we need only the odd-power Wallis integrals
// (which are rational, no pi factor) and the even-power ones (which
// involve pi).

#include "angles.h"
#include "factorial.h"
#include "tracked_value.h"

namespace math {

/// Wallis integral for odd power n = 2k+1 (exact rational, no pi)
/// W_{2k+1} = (2k)!! / (2k+1)!!
template<typename T>
TrackedValue<T> wallis_odd(int k) {
    return double_factorial<T>(2 * k) / double_factorial<T>(2 * k + 1);
}

/// Wallis integral for even power n = 2k (involves pi)
/// W_{2k} = (2k-1)!! / (2k)!! * pi/2
template<typename T>
TrackedValue<T> wallis_even(int k) {
    if (k == 0) {
        // W_0 = pi/2. pi<T>() (from angles.h) is the canonical tracked π
        // (AUD-CC-11); it builds the same (value, 0, value·ε, 0) TrackedValue.
        return pi<T>() / exact<T>(2);
    }
    TrackedValue<T> pi_half = wallis_even<T>(0);
    return double_factorial<T>(2 * k - 1) / double_factorial<T>(2 * k) * pi_half;
}

/// General Wallis integral W_n
template<typename T>
TrackedValue<T> wallis(int n) {
    if (n % 2 == 0) {
        return wallis_even<T>(n / 2);
    } else {
        return wallis_odd<T>(n / 2);
    }
}

/// The simpler integral used in mean gravity and J2n derivations:
/// integral from 0 to pi/2 of sin^{2k}(phi) * cos(phi) dphi = 1/(2k+1)
/// This is exact rational, computed trivially.
template<typename T>
TrackedValue<T> sin_power_cos_integral(int k) {
    return ratio<T>(1, 2 * k + 1);
}

} // namespace math
