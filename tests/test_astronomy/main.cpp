/// @file test_astronomy/main.cpp
/// @brief Tests for the astronomy/ module: FundamentalConstants, DerivedOrbitalElements.

#include <iostream>
#include <iomanip>
#include <cmath>
#include "math/tracked_value.h"
#include "astronomy/solar_system.h"

using T = double;

int main() {
    std::cout << "=== Astronomy Module Tests ===\n\n";

    astronomy::FundamentalConstants<T> fc = astronomy::FundamentalConstants<T>::sgp4_standard();
    astronomy::DerivedOrbitalElements<T> astro = astronomy::DerivedOrbitalElements<T>::compute(fc);

    // Compare against SGP4 deep-space constants from SR3 page 59
    void (*check)(const char*, double, double, double) = [](const char* name, double computed, double legacy, double tol) {
        double diff = std::abs(computed - legacy);
        bool ok = diff <= tol;
        std::cout << (ok ? "PASS" : "FAIL") << " " << name
                  << ": " << std::setprecision(10) << computed
                  << " vs " << legacy
                  << " diff=" << std::scientific << diff << "\n";
    };

    check("ZNS", astro.solar_mean_motion.value, 1.19459e-5, 1e-9);
    check("ZES", astro.solar_eccentricity.value, 0.01675, 1e-10);
    check("ZSINIS", astro.sin_obliquity.value, 0.39785416, 1e-6);
    check("ZCOSIS", astro.cos_obliquity.value, 0.91744867, 1e-6);
    check("ZNL", astro.lunar_mean_motion.value, 1.5835218e-4, 1e-9);
    check("ZEL", astro.lunar_eccentricity.value, 0.05490, 1e-10);

    std::cout << "\n=== Astronomy Module Tests Done ===\n";
    return 0;
}
