/// @file test_math/main.cpp
/// @brief Tests for the math/ module: TrackedValue, series, factorial, wallis, kepler, angles.

#include <iostream>
#include <cmath>
#include "math/tracked_value.h"
#include "math/series.h"
#include "math/factorial.h"
#include "math/wallis.h"
#include "math/kepler.h"
#include "math/angles.h"
#include "math/vector3.h"

using T = double;

int main() {
    std::cout << "=== Math Module Tests ===\n\n";

    // TrackedValue basic arithmetic
    math::TrackedValue<T> a = math::TrackedValue<T>::defined("3.14159265358979");
    math::TrackedValue<T> b = math::exact<T>(2);
    math::TrackedValue<T> c = a * b;
    std::cout << "pi * 2 = " << c.value << " (precision=" << c.errors.precision << ")\n";

    // Series: evaluate q0 series
    math::TrackedValue<T> e_prime = math::TrackedValue<T>(0.0820944379, 0.0, 1e-15, 0.0);
    std::function<math::TrackedValue<T>(int)> q0_term = [&](int n) -> math::TrackedValue<T> {
        int sign = (n % 2 == 1) ? 1 : -1;
        math::TrackedValue<T> coeff = math::ratio<T>(sign * 4 * n, (2*n+1) * (2*n+3));
        math::TrackedValue<T> ep_power = e_prime;
        for (int i = 1; i < 2*n+1; ++i) ep_power = ep_power * e_prime;
        return coeff * ep_power;
    };
    math::TrackedValue<T> two_q0 = math::alternating_series<T>(1, q0_term, 1e-15);
    math::TrackedValue<T> q0 = two_q0 / math::exact<T>(2);
    std::cout << "q0 = " << q0.value << " (published: 7.334625787083e-05)\n";

    // Kepler solver
    math::TrackedValue<T> M = math::TrackedValue<T>(1.5, 0.0, 1e-15, 0.0);
    math::TrackedValue<T> e = math::TrackedValue<T>(0.1, 0.0, 1e-15, 0.0);
    math::TrackedValue<T> E = math::solve_kepler(M, e, 1e-14);
    std::cout << "Kepler E(" << M.value << ", " << e.value << ") = " << E.value << "\n";
    // Verify: E - e*sin(E) should = M
    double residual = E.value - 0.1 * std::sin(E.value) - 1.5;
    std::cout << "  residual: " << residual << "\n";

    // Wallis integrals
    math::TrackedValue<T> w3 = math::wallis_odd<T>(1);
    std::cout << "W_3 (Wallis odd k=1) = " << w3.value << " (expected: 2/3 = 0.6667)\n";

    // Angles
    math::TrackedValue<T> deg90 = math::degrees_to_radians(math::TrackedValue<T>(90.0, 0.0, 1e-15, 0.0));
    std::cout << "90 deg in rad = " << deg90.value << " (expected: 1.5708)\n";

    std::cout << "\n=== Math Module Tests Done ===\n";
    return 0;
}
