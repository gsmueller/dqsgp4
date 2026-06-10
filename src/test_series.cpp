/// Series Evaluation Verification — Ch 5, §5.9
///
/// Verifies all series evaluators and computed values against
/// 60-digit reference values. Run at T = double; reference values
/// from arbitrary-precision computation.
///
/// Tests:
///   1. q0 series value at WGS84 (Table 5.9.1)
///   2. q0' series value
///   3. U0 series value (arctan(e')/e')
///   4. Convergence rate verification (Table 5.9.2, Eq. 5.19)
///   5. Cancellation comparison: closed form vs series (Table 5.9.3)
///   6. Error bound tightness: Leibniz vs geometric (Table 5.9.4)
///   7. TrackedValue error budget decomposition (Table 5.9.5)
///   8. series_sqrt on a value near 1
///   9. geodetic_binomial_coefficient for alpha = -3/2
///  10. make_binomial_evaluator closure for (1+x)^{-1/2}

#include <iostream>
#include <iomanip>
#include <cmath>

#include "math/tracked_value.h"
#include "math/series.h"
#include "math/binomial_series.h"
#include "geodesy/equipotential_ellipsoid.h"

using T = double;
using TV = math::TrackedValue<T>;

// ----------------------------------------------------------------
// Test infrastructure
// ----------------------------------------------------------------

int tests_passed = 0;
int tests_failed = 0;

void check(const char* name, T actual, T expected, T rel_tol = 1e-12) {
    T err = std::abs(actual - expected);
    T scale = std::abs(expected);
    bool ok;
    if (scale > T(0)) {
        ok = (err / scale) < rel_tol;
    } else {
        ok = err < rel_tol;
    }
    if (ok) {
        ++tests_passed;
    } else {
        ++tests_failed;
        std::cerr << "  FAIL: " << name
                  << "  actual=" << std::setprecision(17) << actual
                  << "  expected=" << expected
                  << "  rel_err=" << (scale > T(0) ? err / scale : err) << "\n";
    }
}

// ----------------------------------------------------------------
// WGS84 parameters
// ----------------------------------------------------------------

const T a_wgs84   = 6378137.0;
const T inv_f_wgs84 = 298.257223563;
const T GM_wgs84  = 3.986004418e14;
const T omega_wgs84 = 7.292115e-5;

const T f_wgs84  = 1.0 / inv_f_wgs84;
const T e2_wgs84 = 2.0 * f_wgs84 - f_wgs84 * f_wgs84;
const T ep2_wgs84 = e2_wgs84 / (1.0 - e2_wgs84);
const T ep_wgs84  = std::sqrt(ep2_wgs84);

// Reference values from 60-digit computation
const T q0_ref    = 7.33462578708345207e-05;
const T q0p_ref   = 2.68804130046088833e-03;
const T atan_ep_over_ep_ref = 9.97762541746654685e-01;

// ----------------------------------------------------------------
// Test 1: q0 series value
// ----------------------------------------------------------------

void test_q0_series() {
    std::cout << "--- Test 1: q0 series (Thm 5.5.1) ---\n";

    // Build ellipsoid — q0 is computed in constructor via alternating_series
    // Delegate to the WGS84 library factory — ONE honest source (item F3):
    // no re-typed constants, and GM carries its published sigma (measured)
    // rather than a wrong defined()=exact. q0/q0' depend only on the
    // flattening, so the series values these tests check are unchanged.
    geodesy::EquipotentialEllipsoid<T> ell = geodesy::EquipotentialEllipsoid<T>::wgs84(1e-16);

    std::cout << "  q0 = " << std::setprecision(17) << ell.q0.value << "\n";
    std::cout << "  q0.delta_p = " << ell.q0.errors.precision << "\n";
    check("q0 value", ell.q0.value, q0_ref, 1e-14);
}

// ----------------------------------------------------------------
// Test 2: q0' series value
// ----------------------------------------------------------------

void test_q0p_series() {
    std::cout << "--- Test 2: q0' series (Thm 5.5.2) ---\n";

    // Delegate to the WGS84 library factory — ONE honest source (item F3):
    // no re-typed constants, and GM carries its published sigma (measured)
    // rather than a wrong defined()=exact. q0/q0' depend only on the
    // flattening, so the series values these tests check are unchanged.
    geodesy::EquipotentialEllipsoid<T> ell = geodesy::EquipotentialEllipsoid<T>::wgs84(1e-16);

    std::cout << "  q0' = " << std::setprecision(17) << ell.q0p.value << "\n";
    std::cout << "  q0'.delta_p = " << ell.q0p.errors.precision << "\n";
    check("q0' value", ell.q0p.value, q0p_ref, 1e-14);
}

// ----------------------------------------------------------------
// Test 3: U0 series (arctan(e')/e')
// ----------------------------------------------------------------

void test_u0_series() {
    std::cout << "--- Test 3: U0 auxiliary series (Thm 5.5.3) ---\n";

    // Evaluate arctan(e')/e' directly via alternating_series
    TV ep = math::exact<T>(0) + math::TrackedValue<T>(ep_wgs84, 0.0, 0.0, 1e-15);

    std::function<TV(int)> u0_term = [&](int n) -> TV {
        int sign = (n % 2 == 0) ? 1 : -1;
        TV coeff = math::ratio<T>(sign, 2 * n + 1);
        TV ep_power = math::exact<T>(1);
        for (int i = 0; i < 2 * n; ++i) {
            ep_power = ep_power * ep;
        }
        return coeff * ep_power;
    };

    TV result = math::alternating_series<T>(0, u0_term, 1e-16);
    std::cout << "  arctan(e')/e' = " << std::setprecision(17) << result.value << "\n";
    std::cout << "  delta_p = " << result.errors.precision << "\n";
    check("arctan(e')/e'", result.value, atan_ep_over_ep_ref, 1e-14);
}

// ----------------------------------------------------------------
// Test 4: Convergence rate verification (Eq. 5.19)
// ----------------------------------------------------------------

void test_convergence_rate() {
    std::cout << "--- Test 4: Convergence rate (Eq. 5.19) ---\n";

    // Compute first 10 q0 terms and verify ratios match Eq. 5.19
    std::vector<T> term_values;
    for (int n = 1; n <= 10; ++n) {
        int sign = (n % 2 == 1) ? 1 : -1;
        T coeff = sign * 4.0 * n / ((2.0*n+1) * (2.0*n+3));
        T term = coeff * std::pow(ep_wgs84, 2*n+1);
        term_values.push_back(term);
    }

    std::cout << std::setprecision(12);
    for (int n = 1; n < 10; ++n) {
        T actual_ratio = std::abs(term_values[n] / term_values[n-1]);
        // Eq. 5.19: ((n+1)*(2n+1)) / (n*(2n+5)) * e'^2
        T predicted = ((n+1.0)*(2.0*n+1)) / (n*(2.0*n+5)) * ep2_wgs84;
        std::cout << "  n=" << n << ": actual=" << actual_ratio
                  << "  Eq.5.19=" << predicted
                  << "  e'^2=" << ep2_wgs84 << "\n";
        check("ratio at n", actual_ratio, predicted, 1e-12);
    }
}

// ----------------------------------------------------------------
// Test 5: Cancellation in closed form (Cor. 5.5.2)
// ----------------------------------------------------------------

void test_cancellation() {
    std::cout << "--- Test 5: Cancellation in closed form ---\n";

    T term_a = (1.0 + 3.0 / ep2_wgs84) * std::atan(ep_wgs84);
    T term_b = 3.0 / ep_wgs84;
    T q0_closed = 0.5 * (term_a - term_b);

    std::cout << std::setprecision(17);
    std::cout << "  (1+3/e'^2)*atan(e') = " << term_a << "\n";
    std::cout << "  3/e'                = " << term_b << "\n";
    std::cout << "  q0 closed form      = " << q0_closed << "\n";
    std::cout << "  q0 series (ref)     = " << q0_ref << "\n";

    T rel_err = std::abs(q0_closed - q0_ref) / q0_ref;
    T digits_reliable = -std::log10(rel_err);
    std::cout << "  Reliable digits     = " << std::setprecision(1)
              << digits_reliable << " (expect ~12-13)\n";

    // The closed form should agree to ~10-12 digits, not 15+
    check("closed form agrees roughly", q0_closed, q0_ref, 1e-9);
}

// ----------------------------------------------------------------
// Test 6: Error bound tightness (Table 5.9.4)
// ----------------------------------------------------------------

void test_bound_tightness() {
    std::cout << "--- Test 6: Bound tightness (Leibniz vs geometric) ---\n";

    // Compute q0 series incrementally, track Leibniz and geometric bounds
    // vs actual remainder (computed from high-N sum)
    T ref_2q0 = 2.0 * q0_ref;

    T partial = 0.0;
    std::cout << std::setprecision(3);
    std::cout << "  N   Leibniz      Geometric    Actual       L/A     G/A\n";

    for (int n = 1; n <= 10; ++n) {
        int sign = (n % 2 == 1) ? 1 : -1;
        T coeff = sign * 4.0 * n / ((2.0*n+1) * (2.0*n+3));
        T term = coeff * std::pow(ep_wgs84, 2*n+1);
        partial += term;

        T leibniz = std::abs(term);
        T geom = std::abs(term) * ep2_wgs84 / (1.0 - ep2_wgs84);
        T actual = std::abs(ref_2q0 - partial);

        if (actual > 0) {
            std::cout << "  " << n
                      << "   " << std::scientific << leibniz
                      << "  " << geom
                      << "  " << actual
                      << "  " << std::fixed << std::setprecision(0)
                      << leibniz / actual
                      << "      " << std::setprecision(2) << geom / actual
                      << "\n";
        }
    }

    // Verify geometric bound is always >= actual (rigorous)
    // Limit to N=5 because at double precision the partial sum
    // accumulates ~1e-19 arithmetic rounding. For N > 5, the
    // mathematical remainder is comparable to this noise floor,
    // so "actual remainder" computed from ref - partial is unreliable.
    partial = 0.0;
    for (int n = 1; n <= 5; ++n) {
        int sign = (n % 2 == 1) ? 1 : -1;
        T coeff = sign * 4.0 * n / ((2.0*n+1) * (2.0*n+3));
        T term = coeff * std::pow(ep_wgs84, 2*n+1);
        partial += term;
        T geom = std::abs(term) * ep2_wgs84 / (1.0 - ep2_wgs84);
        T actual = std::abs(ref_2q0 - partial);
        check("geom >= actual", (geom >= actual) ? 1.0 : 0.0, 1.0, 1e-15);
    }
}

// ----------------------------------------------------------------
// Test 7: TrackedValue error budget (Table 5.9.5)
// ----------------------------------------------------------------

void test_error_budget() {
    std::cout << "--- Test 7: TrackedValue error budget ---\n";

    // Delegate to the WGS84 library factory — ONE honest source (item F3):
    // no re-typed constants, and GM carries its published sigma (measured)
    // rather than a wrong defined()=exact. q0/q0' depend only on the
    // flattening, so the series values these tests check are unchanged.
    geodesy::EquipotentialEllipsoid<T> ell = geodesy::EquipotentialEllipsoid<T>::wgs84(1e-16);

    std::cout << "  q0.value       = " << std::setprecision(17) << ell.q0.value << "\n";
    std::cout << "  q0.delta_p     = " << std::scientific << ell.q0.errors.precision << "\n";
    std::cout << "  q0.sigma_m     = " << ell.q0.errors.measurement << "\n";
    std::cout << "  q0.delta_a     = " << ell.q0.errors.accuracy << "\n";

    // delta_p should be small (truncation + arithmetic rounding)
    check("q0 delta_p < 1e-15", (ell.q0.errors.precision < 1e-15) ? 1.0 : 0.0, 1.0, 1e-15);
}

// ----------------------------------------------------------------
// Test 8: series_sqrt (Lemma 5.6.1)
// ----------------------------------------------------------------

void test_series_sqrt() {
    std::cout << "--- Test 8: series_sqrt (Lemma 5.6.1) ---\n";

    // Compute sqrt(1 + U) where U = 0.007 (WGS84-like)
    TV I = math::exact<T>(1) + math::TrackedValue<T>(0.007, 0.0, 0.0, 1e-15);
    TV result = math::series_sqrt(I, 1e-16);

    T expected = std::sqrt(1.007);
    std::cout << "  series_sqrt(1.007) = " << std::setprecision(17) << result.value << "\n";
    std::cout << "  std::sqrt(1.007)   = " << expected << "\n";
    std::cout << "  delta_p            = " << std::scientific << result.errors.precision << "\n";
    check("series_sqrt", result.value, expected, 1e-15);
}

// ----------------------------------------------------------------
// Test 9: geodetic_binomial_coefficient (Thm 5.6.2)
// ----------------------------------------------------------------

void test_geodetic_binomial() {
    std::cout << "--- Test 9: geodetic_binomial_coefficient ---\n";

    // Meridian arc: alpha = -3/2, w = 1
    // Coefficients: C(-3/2, k) * (2k-1)!!/(2k)!!
    // k=0: 1, k=1: -3/4, k=2: 45/64, k=3: -175/256, k=4: 11025/16384
    TV alpha = math::ratio<T>(-3, 2);

    T expected[] = { 1.0, -3.0/4.0, 45.0/64.0, -175.0/256.0, 11025.0/16384.0 };

    for (int k = 0; k <= 4; ++k) {
        TV c_k = math::geodetic_binomial_coefficient(alpha, k, false);
        std::cout << "  c_" << k << " = " << std::setprecision(12) << c_k.value
                  << "  (expected " << expected[k] << ")\n";
        check("meridian arc coeff", c_k.value, expected[k], 1e-14);
    }
}

// ----------------------------------------------------------------
// Test 10: make_binomial_evaluator closure
// ----------------------------------------------------------------

void test_binomial_evaluator() {
    std::cout << "--- Test 10: make_binomial_evaluator ---\n";

    // (1 + x)^{-1/2} evaluated at x = -e^2*sin^2(phi) for phi = pi/4
    TV alpha = math::ratio<T>(-1, 2);
    std::function<TV(const TV&)> eval = math::make_binomial_evaluator<T>(alpha, 1e-16);

    T sin2_phi = 0.5;  // sin^2(pi/4)
    TV x = math::exact<T>(0) - math::TrackedValue<T>(e2_wgs84 * sin2_phi, 0.0, 0.0, 1e-15);

    TV result = eval(x);
    T expected = 1.0 / std::sqrt(1.0 - e2_wgs84 * sin2_phi);

    std::cout << "  (1-e^2*0.5)^{-1/2} series = " << std::setprecision(17) << result.value << "\n";
    std::cout << "  direct                     = " << expected << "\n";
    std::cout << "  delta_p                    = " << std::scientific << result.errors.precision << "\n";
    check("binomial evaluator", result.value, expected, 1e-14);
}

// ----------------------------------------------------------------
// Test 11: alternating_series with convergence_ratio (tighter bound)
// ----------------------------------------------------------------

void test_alternating_with_ratio() {
    std::cout << "--- Test 11: alternating_series with convergence_ratio ---\n";

    // Evaluate q0 with Leibniz (default) and with geometric bound (r = e'^2)
    TV ep = math::TrackedValue<T>(ep_wgs84, 0.0, 0.0, 1e-15);

    std::function<TV(int)> q0_term = [&](int n) -> TV {
        int sign = (n % 2 == 1) ? 1 : -1;
        TV coeff = math::ratio<T>(sign * 4 * n, (2*n+1) * (2*n+3));
        TV ep_power = ep;
        for (int i = 1; i < 2*n+1; ++i) {
            ep_power = ep_power * ep;
        }
        return coeff * ep_power;
    };

    // Leibniz (default)
    TV result_leibniz = math::alternating_series<T>(1, q0_term, 1e-16);
    T dp_leibniz = result_leibniz.errors.precision;

    // Geometric (tighter)
    TV result_geom = math::alternating_series<T>(1, q0_term, 1e-16, 10000, ep2_wgs84);
    T dp_geom = result_geom.errors.precision;

    std::cout << "  Leibniz delta_p  = " << std::scientific << dp_leibniz << "\n";
    std::cout << "  Geometric delta_p = " << dp_geom << "\n";
    std::cout << "  Ratio (Leibniz/Geometric) = " << std::fixed << std::setprecision(0)
              << dp_leibniz / dp_geom << "\n";

    // Values should agree to high precision (may differ by ~1 term's worth)
    check("same value", result_leibniz.value, result_geom.value, 1e-12);

    // Both bounds should be well below tolerance
    check("leibniz dp < 1e-15", (dp_leibniz < 1e-15) ? 1.0 : 0.0, 1.0, 1e-15);
    check("geom dp < 1e-15", (dp_geom < 1e-15) ? 1.0 : 0.0, 1.0, 1e-15);
}

// ----------------------------------------------------------------

int main() {
    std::cout << std::setprecision(17);
    std::cout << "WGS84: e^2 = " << e2_wgs84 << ", e'^2 = " << ep2_wgs84
              << ", e' = " << ep_wgs84 << "\n\n";

    test_q0_series();
    test_q0p_series();
    test_u0_series();
    test_convergence_rate();
    test_cancellation();
    test_bound_tightness();
    test_error_budget();
    test_series_sqrt();
    test_geodetic_binomial();
    test_binomial_evaluator();
    test_alternating_with_ratio();

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << tests_passed << "  Failed: " << tests_failed << "\n";
    std::cout << "========================================\n";

    return tests_failed > 0 ? 1 : 0;
}
