#pragma once

/// @file bessel.h
/// Tracked Bessel function of the first kind J_n(x), integer order n ≥ 0.
///
/// The value is computed in the numeric type T via boost::math::cyl_bessel_j, so
/// a wider T sharpens the tracked PRECISION (the calling card). The input error δ
/// is propagated through the derivative
///   J_n'(x) = (J_{n−1}(x) − J_{n+1}(x))/2     (J_0'(x) = −J_1(x)),
/// with a δ²/2 second-order remainder (|J_n''| ≤ 1) — the same linear+quadratic
/// forward bound the tracked sin/cos use — and capped at 1 (|J_n| ≤ 1 ∀ real x).
///
/// Used by the Fourier–Bessel solution of Kepler's equation
/// (orbit/kepler_series.h, design/derivations/ephemeris_series.md).

#include "tracked_value.h"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>

namespace math {

/// Bessel function of the first kind J_n(x), tracked: the value via
/// boost::math::cyl_bessel_j evaluated in T, the input-error propagation via
/// the derivative identity J_n'(x) = (J_{n-1}(x) − J_{n+1}(x))/2 per category.
template<typename T>
TrackedValue<T> bessel_j(int n, const TrackedValue<T>& x) {
    using std::abs;
    using std::min;
    using boost::math::cyl_bessel_j;

    T xv = x.value;
    T jn = cyl_bessel_j(n, xv);
    // J_n'(x) = (J_{n−1} − J_{n+1})/2, with the n = 0 special case J_0' = −J_1.
    T jp = (n == 0) ? -cyl_bessel_j(1, xv)
                    : (cyl_bessel_j(n - 1, xv) - cyl_bessel_j(n + 1, xv)) / T(2);
    T absjp = abs(jp);

    return TrackedValue<T>(jn, x.errors.apply([&](T err) -> T {
        // |J_n(x+δ) − J_n(x)| ≤ |J_n'(x)|·δ + δ²/2  (|J_n''| ≤ 1), and ≤ 1 always.
        return min(absjp * err + err * err / T(2), T(1));
    }));
}

} // namespace math
