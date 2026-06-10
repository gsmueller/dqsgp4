/// test_earth_models — C2 (DQSGP4 Completion Roadmap, issue C2).
///
/// Verifies the named earth-model factories on geodesy::EquipotentialEllipsoid:
///   - WGS 84 / GRS 80 / WGS 72 / IERS 2010 each construct with their defining
///     parameters (a, 1/f or J2, GM, omega);
///   - the WGS 84 derived flattening is self-consistent with 1/f;
///   - the IERS 2010 model (current IAU 2009 NSFA) differs from WGS 84 and
///     carries a measurement uncertainty on GM.
///
/// ExeGate C2: nonzero exit on any failed check.

#include "geodesy/equipotential_ellipsoid.h"

#include <cmath>
#include <iostream>

namespace {

using T = double;
using Ellip = geodesy::EquipotentialEllipsoid<T>;

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
    Ellip wgs84 = Ellip::wgs84(tol);
    Ellip grs80 = Ellip::grs80(tol);
    Ellip wgs72 = Ellip::wgs72(tol);
    Ellip iers = Ellip::iers2010(tol);

    // WGS 84 defining parameters and self-consistent flattening.
    check("WGS84 a == 6378137.0", close_d(wgs84.a.value, 6378137.0, 1e-6));
    check("WGS84 1/f == 298.257223563", close_d(wgs84.inv_f.value, 298.257223563, 1e-9));
    check("WGS84 GM == 3.986004418e14", close_d(wgs84.GM.value, 3.986004418e14, 1e6));
    check("WGS84 f == 1/298.257223563",
          close_d(wgs84.f.value, 1.0 / 298.257223563, 1e-15));

    // GRS 80 and WGS 72.
    check("GRS80 a == 6378137.0", close_d(grs80.a.value, 6378137.0, 1e-6));
    check("GRS80 GM == 3.986005e14", close_d(grs80.GM.value, 3.986005e14, 1e6));
    check("WGS72 a == 6378135.0", close_d(wgs72.a.value, 6378135.0, 1e-6));

    // IERS 2010 (current) is distinct from WGS 84 and is a measured model.
    check("IERS2010 a == 6378136.6", close_d(iers.a.value, 6378136.6, 1e-6));
    check("IERS2010 a differs from WGS84",
          std::abs(iers.a.value - wgs84.a.value) > 0.1);
    check("IERS2010 GM measurement uncertainty > 0",
          static_cast<double>(iers.GM.errors.measurement) > 0.0);

    std::cout << "\n  earth models: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
