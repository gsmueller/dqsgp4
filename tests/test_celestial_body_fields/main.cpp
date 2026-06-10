/// test_celestial_body_fields — BUG1 (DQSGP4 Completion Roadmap, issue BUG1).
///
/// Verifies that ephemeris::CelestialBody::make_sun / make_moon populate every
/// field with a sensible value — no value-initialized traps. The never-assigned
/// `mu_over_r3` member was removed (it returned 0 if ever read, could not be
/// computed from FundamentalConstants, and was redundant with the
/// `perturbation_coef` passed to perturbation::compute_third_body_rates).
///
/// ExeGate BUG1: nonzero exit on any failed check.

#include "astronomy/solar_system.h"
#include "ephemeris/celestial_body.h"

#include <cmath>
#include <iostream>

namespace {

using T = double;

int passed = 0;
int failed = 0;

void check(const char* name, bool ok) {
    if (ok) {
        ++passed;
        std::cout << "  PASS: " << name << "\n";
    } else {
        ++failed;
        std::cout << "  FAIL: " << name << "\n";
    }
}

bool close_d(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

} // namespace

int main() {
    astronomy::FundamentalConstants<T> fc =
        astronomy::FundamentalConstants<T>::sgp4_standard();
    ephemeris::CelestialBody<T> sun = ephemeris::CelestialBody<T>::make_sun(fc);
    ephemeris::CelestialBody<T> moon = ephemeris::CelestialBody<T>::make_moon(fc);

    // Sun: every field populated and sensible.
    check("sun name is Sun", sun.name == "Sun");
    check("sun mean-anomaly rate > 0", sun.mean_anomaly_rate_rad_min.value > 0.0);
    check("sun longitude rate > 0", sun.longitude_rate_rad_day.value > 0.0);
    check("sun eccentricity matches fc",
          close_d(sun.orbit_eccentricity.value, fc.solar_eccentricity.value, 1e-15));
    check("sun inclination == 0 (on the ecliptic)", sun.orbit_inclination.value == 0.0);

    // Moon: every field populated and sensible.
    check("moon name is Moon", moon.name == "Moon");
    check("moon mean-anomaly rate > 0", moon.mean_anomaly_rate_rad_min.value > 0.0);
    check("moon longitude rate > 0", moon.longitude_rate_rad_day.value > 0.0);
    check("moon eccentricity matches fc",
          close_d(moon.orbit_eccentricity.value, fc.lunar_eccentricity.value, 1e-15));
    check("moon inclination > 0 (~5.145 deg)", moon.orbit_inclination.value > 0.0);
    check("moon node rate < 0 (regresses)", moon.node_rate_rad_day.value < 0.0);

    std::cout << "\n  celestial body fields: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
