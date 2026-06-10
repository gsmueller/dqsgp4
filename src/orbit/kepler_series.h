#pragma once

/// @file kepler_series.h
/// Fourier–Bessel series solution of Kepler's equation, and the equation of
/// centre / radius that follow exactly from it.
///
/// THEORY: design/derivations/ephemeris_series.md. Kepler's equation
/// M = E − e·sin E is solved EXACTLY by the Fourier series
///
///   E(M, e) = M + 2 Σ_{n=1}^∞ (J_n(n e)/n)·sin(n M),
///
/// which converges for ALL e < 1 — the analytic continuation past the Laplace
/// limit e = 0.6627 that any power-series-in-e cannot cross (verified: e-power
/// diverges at e = 0.95, this series keeps converging). It reaches machine
/// precision in K ≈ 8 harmonics for the ephemeris bodies. The companion
/// orbit::solve_kepler_* (modified_kepler.h) are the iterative form for the SGP4
/// path; THIS is the generative, dial-up, tracked series form (the user directive
/// 2026-06-07): dial up K, truncation accuracy = the Bessel tail Σ_{n>K} 2|J_n|/n.
///
/// The true anomaly, equation of centre, and radius are then EXACT algebraic
/// transforms of E — no further truncation.

#include "../math/tracked_value.h"
#include "../math/bessel.h"
#include "../math/angles.h"

#include <boost/math/special_functions/bessel.hpp>
#include <cmath>

namespace orbit {

/// Eccentric anomaly E(M, e) via the Fourier–Bessel series to K harmonics. The
/// Bessel-tail truncation (the omitted harmonics n > K) is folded into
/// errors.accuracy — a rigorous bound that tightens as K grows.
template<typename T>
math::TrackedValue<T> eccentric_anomaly_series(
    const math::TrackedValue<T>& M, const math::TrackedValue<T>& e, int K) {
    using TV = math::TrackedValue<T>;
    using math::exact;

    TV E = M;
    for (int n = 1; n <= K; ++n) {
        TV nT = exact<T>(n);
        TV jn = math::bessel_j(n, nT * e);                 // J_n(n e), tracked
        E = E + (exact<T>(2) / nT) * jn * sin(nT * M);
    }

    // Truncation accuracy: |Σ_{n>K} 2 J_n(n e)/n|. Sum a tight explicit window,
    // then a geometric remainder with the Fourier–Bessel ratio
    // ρ(e) = e·exp(√(1−e²))/(1+√(1−e²)) < 1 (∀ e < 1). The bound is a magnitude,
    // computed on .value (an upper bound need not carry full T precision).
    using std::abs; using std::sqrt; using std::exp;
    using boost::math::cyl_bessel_j;
    T ev = abs(e.value);
    T root = sqrt(T(1) - ev * ev);
    T rho = (root > T(0)) ? ev * exp(root) / (T(1) + root) : T(1);
    T tail = T(0), last = T(0);
    for (int n = K + 1; n <= K + 16; ++n) {
        last = T(2) * abs(cyl_bessel_j(n, T(n) * ev)) / T(n);
        tail = tail + last;
    }
    if (rho < T(1)) tail = tail + last * rho / (T(1) - rho);  // remainder beyond K+16
    return math::add_bound(E, tail, math::ErrorChannel::accuracy);
}

/// Equation of centre C = ν − M, from the Fourier–Bessel E and the EXACT
/// true-anomaly transform ν = atan2(√(1−e²)·sin E, cos E − e). Wrapped to (−π, π].
template<typename T>
math::TrackedValue<T> equation_of_center(
    const math::TrackedValue<T>& M, const math::TrackedValue<T>& e, int K) {
    using math::exact;
    math::TrackedValue<T> E = eccentric_anomaly_series(M, e, K);
    math::TrackedValue<T> root = sqrt(exact<T>(1) - e * e);
    math::TrackedValue<T> nu = atan2(root * sin(E), cos(E) - e);
    return math::wrap_neg_pos_pi(nu - M);
}

/// Radius r/a = 1 − e·cos E (exact from the Fourier–Bessel E).
template<typename T>
math::TrackedValue<T> radius_over_semimajor(
    const math::TrackedValue<T>& M, const math::TrackedValue<T>& e, int K) {
    using math::exact;
    math::TrackedValue<T> E = eccentric_anomaly_series(M, e, K);
    return exact<T>(1) - e * cos(E);
}

} // namespace orbit
