#pragma once

/// Kepler equation solver: E - e*sin(E) = M
///
/// Uses Halley's method (cubic convergence).
/// The caller specifies a tolerance. The solver computes the number of
/// iterations needed from the convergence rate and the starting error,
/// then iterates until the tolerance is met.
///
/// All operations on TrackedValue<T> — errors propagate through every step.
/// The convergence residual is added to the precision error of the result.
///
/// Theoretical basis: design/audit/theoretical_basis_audit/kepler.md (the
/// audited Halley-iteration derivation). The generative Fourier–Bessel
/// companion (dial-up series, converges for all e < 1) is
/// design/derivations/ephemeris_series.md / orbit/kepler_series.h.

#include "tracked_value.h"

namespace math {

/// Solve E - e*sin(E) = M for E.
///
/// Halley's method converges cubically: correct digits triple each iteration.
/// Starting from a 3-4 digit approximation, the iteration count to reach
/// the target tolerance is ceil(log3(target_digits / starter_digits)).
///
/// For 50 digits: ~6 iterations. For 100 digits: ~7. For double: ~3.
///
/// @param M         mean anomaly [rad]
/// @param e         eccentricity (0 <= e < 1)
/// @param tolerance convergence tolerance — iteration stops when |correction| < tolerance
/// @return          eccentric anomaly E [rad] with precision error from convergence
template<typename T>
TrackedValue<T> solve_kepler(
    const TrackedValue<T>& M,
    const TrackedValue<T>& e,
    const T& tolerance)
{
    // Starter: E_0 = M + e*sin(M) + (e²/2)*sin(2M)
    // Gives ~3-4 correct digits for e < 0.3
    // For high e, this still converges but needs more iterations
    TrackedValue<T> E = M + e * sin(M)
                        + e * e * sin(M + M) / exact<T>(2);

    // Estimate initial error: |E_0 - E_true| ≈ e³ (rough bound)
    // For e = 0.1, initial error ≈ 0.001 → ~3 correct digits
    // Halley triples digits: 3 → 9 → 27 → 81 ...
    // Needed iterations = ceil(log(target_digits / 3) / log(3))
    // But we don't precompute — we just iterate until convergence

    // Safety limit: even for 1000-digit types, 20 iterations of Halley
    // gives 3^20 ≈ 3.5 billion digits — far more than any practical type
    constexpr int safety_limit = 20;

    TrackedValue<T> last_correction;

    for (int iter = 0; iter < safety_limit; ++iter) {
        TrackedValue<T> sinE = sin(E);
        TrackedValue<T> cosE = cos(E);

        // f = E - e*sin(E) - M
        TrackedValue<T> f  = E - e * sinE - M;
        // f' = 1 - e*cos(E)
        TrackedValue<T> fp = exact<T>(1) - e * cosE;
        // f'' = e*sin(E)
        TrackedValue<T> fpp = e * sinE;

        // Halley correction: delta = 2*f*f' / (2*f'² - f*f'')
        TrackedValue<T> numerator = exact<T>(2) * f * fp;
        TrackedValue<T> denominator = exact<T>(2) * fp * fp - f * fpp;
        TrackedValue<T> delta = numerator / denominator;

        E = E - delta;
        last_correction = delta;

        using std::abs;
        if (abs(delta.value) < tolerance) {
            // Convergence achieved.
            // The true error is bounded by |delta|³ / (something) due to
            // cubic convergence, but conservatively we use |delta| itself.
            E.errors.precision = E.errors.precision + abs(delta.value);
            return E;
        }
    }

    // If we reach here, something is wrong (shouldn't happen for e < 1).
    // Return best estimate with last correction as precision error.
    using std::abs;
    E.errors.precision = E.errors.precision + abs(last_correction.value);
    return E;
}

} // namespace math
