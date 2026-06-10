/// test_constants_single_source (gate F3) — the WGS-convention constants have
/// ONE source of truth (the geodesy EquipotentialEllipsoid factories, C2), and
/// every consumer must stay bound to it:
///
///   • the DQ `ConstantsProvider` bundle now *delegates* to that source, so its
///     ellipsoid must be BYTE-IDENTICAL (F3-a — the former inline copy is gone);
///
///   • the authentic SGP4 model carries the SAME physical ellipsoid expressed in
///     its load-bearing km working unit, so it must stay SCALE-CONSISTENT with
///     the SI source (a·1000 = aₘ, GM·1e9 = GMₘ, dimensionless parts equal).
///
/// This gate fails if anyone edits one representation's defining numbers without
/// the other — the duplication drift risk the F3 item names. The authentic SGP4
/// construction is deliberately untouched (it is the DQSGP4 test oracle); this
/// test only *observes* it.

#include "constants/constants_provider.h"
#include "geodesy/equipotential_ellipsoid.h"
#include "math/tracked_value.h"
#include "sgp4/model_selector.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

using T = double;
using TV = math::TrackedValue<T>;

int failures = 0;

/// Two encodings of the same constant in the same unit must be bit-identical.
void exact_eq(const std::string& what, const TV& x, const TV& y) {
    if (x.value != y.value) {
        std::cout << "  FAIL " << what << ": " << x.value << " != " << y.value << "\n";
        ++failures;
    } else {
        std::cout << "  ok   " << what << " = " << x.value << "\n";
    }
}

/// A km-unit value times its unit scale must match the SI value to rounding
/// (derived quantities may differ by a few ULP between the km and m chains).
void scale_eq(const std::string& what, const TV& km, double scale, const TV& si,
              double rel_tol) {
    double k = km.value * scale;
    double s = si.value;
    double rel = std::abs(k - s) / std::max(std::abs(s), 1.0);
    if (rel > rel_tol) {
        std::cout << "  FAIL " << what << ": km*scale=" << k << " vs si=" << s
                  << " (rel " << rel << " > " << rel_tol << ")\n";
        ++failures;
    } else {
        std::cout << "  ok   " << what << " (rel " << rel << ")\n";
    }
}

} // namespace

int main() {
    const T tol = T(1e-15);

    // Single source (SI / metres).
    geodesy::EquipotentialEllipsoid<T> si72 = geodesy::EquipotentialEllipsoid<T>::wgs72(tol);
    geodesy::EquipotentialEllipsoid<T> si84 = geodesy::EquipotentialEllipsoid<T>::wgs84(tol);
    geodesy::EquipotentialEllipsoid<T> si80 = geodesy::EquipotentialEllipsoid<T>::grs80(tol);

    // Consumer 1 — DQ ConstantsProvider: must be byte-identical (delegation).
    std::cout << "DQ ConstantsProvider == geodesy single source (byte-identical delegation):\n";
    geodesy::EquipotentialEllipsoid<T> dq72 = constants::ConstantsProvider<T>::wgs72(tol).earth;
    geodesy::EquipotentialEllipsoid<T> dq84 = constants::ConstantsProvider<T>::wgs84(tol).earth;
    geodesy::EquipotentialEllipsoid<T> dq80 = constants::ConstantsProvider<T>::grs80(tol).earth;
    exact_eq("wgs72.a", dq72.a, si72.a);
    exact_eq("wgs72.GM", dq72.GM, si72.GM);
    exact_eq("wgs72.inv_f", dq72.inv_f, si72.inv_f);
    exact_eq("wgs72.omega", dq72.omega, si72.omega);
    exact_eq("wgs72.e2", dq72.e2, si72.e2);
    exact_eq("wgs84.a", dq84.a, si84.a);
    exact_eq("wgs84.GM", dq84.GM, si84.GM);
    exact_eq("wgs84.inv_f", dq84.inv_f, si84.inv_f);
    exact_eq("wgs84.e2", dq84.e2, si84.e2);
    exact_eq("grs80.a", dq80.a, si80.a);
    exact_eq("grs80.GM", dq80.GM, si80.GM);
    exact_eq("grs80.inv_f", dq80.inv_f, si80.inv_f);

    // Consumer 2 — authentic SGP4 km WGS72: must be scale-consistent with SI.
    std::cout << "\nSGP4 km WGS72 <-> geodesy SI WGS72 (scale-consistent):\n";
    geodesy::EquipotentialEllipsoid<T> km72 =
        sgp4::ModelSelector<T>::select("sgp4_standard", tol).ellipsoid;
    std::cout << "  (SGP4 ellipsoid a = " << km72.a.value << " km, GM = "
              << km72.GM.value << " km^3/s^2)\n";
    scale_eq("a  (km->m)", km72.a, 1000.0, si72.a, 1e-12);
    scale_eq("GM (km^3->m^3)", km72.GM, 1e9, si72.GM, 1e-12);
    scale_eq("omega", km72.omega, 1.0, si72.omega, 1e-12);
    scale_eq("inv_f", km72.inv_f, 1.0, si72.inv_f, 1e-12);
    scale_eq("e2", km72.e2, 1.0, si72.e2, 1e-12);

    std::cout << "\n" << (failures == 0 ? "PASS" : "FAIL") << " — " << failures
              << " divergence(s) from the single source\n";
    return failures == 0 ? 0 : 1;
}
