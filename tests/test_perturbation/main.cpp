/// @file test_perturbation/main.cpp
/// @brief Tests for the perturbation/ module: Brouwer rates, Kaula functions.

#include <iostream>
#include <iomanip>
#include <cmath>
#include "math/tracked_value.h"
#include "perturbation/brouwer.h"
#include "perturbation/kaula.h"

using T = double;

int main() {
    std::cout << "=== Perturbation Module Tests ===\n\n";

    // Test Kaula inclination functions at ISS inclination
    double inc_rad = 51.6443 * 3.14159265358979 / 180.0;
    math::TrackedValue<T> sin_i = math::TrackedValue<T>(std::sin(inc_rad), 0.0, 1e-15, 0.0);
    math::TrackedValue<T> cos_i = math::TrackedValue<T>(std::cos(inc_rad), 0.0, 1e-15, 0.0);

    void (*check)(const char*, double, double, double) = [](const char* name, double computed, double legacy, double tol) {
        double diff = std::abs(computed - legacy);
        bool ok = diff <= tol;
        std::cout << (ok ? "PASS" : "FAIL") << " " << name
                  << ": " << std::setprecision(12) << computed
                  << " diff=" << std::scientific << diff << "\n";
    };

    math::TrackedValue<T> f220 = perturbation::inclination_function(2, 2, 0, sin_i, cos_i);
    check("F_220", f220.value, 0.75 * std::pow(1.0 + std::cos(inc_rad), 2), 1e-12);

    math::TrackedValue<T> f221 = perturbation::inclination_function(2, 2, 1, sin_i, cos_i);
    check("F_221", f221.value, 1.5 * std::pow(std::sin(inc_rad), 2), 1e-12);

    math::TrackedValue<T> f442 = perturbation::inclination_function(4, 4, 2, sin_i, cos_i);
    check("F_442", f442.value, 39.375 * std::pow(std::sin(inc_rad), 4), 1e-10);

    // Test Brouwer secular rates
    math::TrackedValue<T> n = math::TrackedValue<T>(0.0472, 0.0, 1e-15, 0.0);  // ~133 min period
    math::TrackedValue<T> a = math::TrackedValue<T>(1.1, 0.0, 1e-15, 0.0);
    math::TrackedValue<T> e2 = math::TrackedValue<T>(0.01, 0.0, 1e-15, 0.0);
    math::TrackedValue<T> J2 = math::TrackedValue<T>(0.001082616, 0.0, 1e-15, 0.0);
    math::TrackedValue<T> J4 = math::TrackedValue<T>(-0.00000165597, 0.0, 1e-15, 0.0);

    perturbation::BrouwerSecularRates<T> rates = perturbation::compute_secular_rates(n, a, e2, cos_i, J2, J4);
    std::cout << "\nBrouwer rates (ISS-like orbit):\n";
    std::cout << "  M_dot     = " << std::setprecision(12) << rates.M_dot.value << " rad/min\n";
    std::cout << "  omega_dot = " << rates.omega_dot.value << " rad/min\n";
    std::cout << "  Omega_dot = " << rates.Omega_dot.value << " rad/min\n";
    std::cout << "  accuracy  = " << std::scientific << rates.M_dot.errors.accuracy << "\n";

    std::cout << "\n=== Perturbation Module Tests Done ===\n";
    return 0;
}
