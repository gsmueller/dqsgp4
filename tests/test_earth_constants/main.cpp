/// test_earth_constants — C1 (DQSGP4 Completion Roadmap, issue C1).
///
/// The modern (IAU 2006 / Almanac) astronomy preset must encode the unambiguous
/// J2000 quantities from authoritative born-digital sources, per CR1:
///   - obliquity = 84381.406" (IAU 2006 B1 / IAU 2009 NSFA), via the generative
///     arcsec * pi/648000 form -> 0.40909280422 rad, with mas-level precision;
///   - Earth orbital (solar) eccentricity 0.016708634 (J2000 VSOP, the leading
///     term of the Meeus 25.4 secular series — see test_series_constants/SC2);
///   - lunar mean eccentricity 0.0549006;
/// and these must differ from the SR3 (1970s) preset.
///
/// ExeGate C1: nonzero exit on any failed check.

#include "sgp4/model_selector.h"

#include <cmath>
#include <iostream>

namespace {

using T = double;
const double PI = 3.14159265358979323846;

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
    const T tol = T(1e-12);
    sgp4::ModelConfiguration<T> modern = sgp4::ModelSelector<T>::select("research_full", tol);
    sgp4::ModelConfiguration<T> sr3 = sgp4::ModelSelector<T>::select("sgp4_standard", tol);

    // Obliquity: 84381.406" * pi/648000 = 0.40909280422 rad (IAU 2006 B1).
    const double eps = modern.astro_constants.obliquity.value;
    const double eps_expected = 84381.406 * PI / 648000.0;
    check("obliquity == IAU 2006 84381.406 arcsec (0.40909280 rad)",
          close_d(eps, eps_expected, 1e-12));

    const double eps_sr3 = sr3.astro_constants.obliquity.value;
    check("modern obliquity differs from the SR3 23.4441 deg value",
          std::abs(eps - eps_sr3) > 1e-5);

    // The generative arcsec form is precise to the mas adoption (~2.4e-9 rad),
    // far tighter than a 6-figure degree truncation (~8.7e-7 rad).
    check("obliquity precision <= 1e-8 rad (generative, mas-bounded)",
          static_cast<double>(modern.astro_constants.obliquity.errors.precision) <= 1e-8);

    // Series-generated (astronomy::earth_eccentricity) at t=0: the born-digital
    // leading coefficient 0.016708634 (Meeus 25.4 / VSOP87), sharper than the
    // former 7-figure 0.0167086 stamp. The full secular series is gated in SC2.
    check("solar eccentricity == 0.016708634 (J2000 VSOP leading term)",
          close_d(modern.astro_constants.solar_eccentricity.value, 0.016708634, 1e-12));

    check("lunar eccentricity == 0.0549006",
          close_d(modern.astro_constants.lunar_eccentricity.value, 0.0549006, 1e-12));

    std::cout << "\n  earth/almanac constants: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
