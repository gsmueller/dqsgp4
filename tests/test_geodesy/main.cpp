/// @file test_geodesy/main.cpp
/// @brief Tests for the geodesy/ module: EquipotentialEllipsoid.

#include <iostream>
#include <iomanip>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "math/tracked_value.h"
#include "geodesy/equipotential_ellipsoid.h"

using T = boost::multiprecision::cpp_bin_float_50;

int main() {
    std::cout << "=== Geodesy Module Tests (50-digit) ===\n\n";

    math::TrackedValue<T> a     = math::TrackedValue<T>::defined("6378137.0");
    math::TrackedValue<T> inv_f = math::TrackedValue<T>::defined("298.257223563");
    math::TrackedValue<T> GM    = math::TrackedValue<T>::measured("3.986004418e14", "8e6");
    math::TrackedValue<T> omega = math::TrackedValue<T>::defined("7.292115e-5");

    geodesy::EquipotentialEllipsoid<T> wgs84(a, inv_f, GM, omega, T("1e-45"));

    std::cout << std::setprecision(15);
    std::cout << "f     = " << wgs84.f.value << " (pub: 3.3528106647475e-03)\n";
    std::cout << "e2    = " << wgs84.e2.value << " (pub: 6.694379990141e-03)\n";
    std::cout << "b     = " << wgs84.b.value << " (pub: 6356752.3142)\n";
    std::cout << "q0    = " << wgs84.q0.value << " (pub: 7.334625787083e-05)\n";
    std::cout << "gamma_e = " << wgs84.gamma_e.value << " (pub: 9.7803253359)\n";
    std::cout << "J2    = " << wgs84.J2n(1).value << " (pub: 1.082629821313e-03)\n";

    // Test J2 init path (WGS72)
    math::TrackedValue<T> a72 = math::TrackedValue<T>::defined("6378.135");
    math::TrackedValue<T> J2_72 = math::TrackedValue<T>::defined("0.001082616");
    math::TrackedValue<T> GM72 = math::TrackedValue<T>::defined("398600.8");
    math::TrackedValue<T> om72 = math::TrackedValue<T>::defined("7.292115e-5");
    geodesy::EquipotentialEllipsoid<T> wgs72 = geodesy::EquipotentialEllipsoid<T>::from_J2(a72, J2_72, GM72, om72, T("1e-45"));
    std::cout << "\nWGS72 (from J2):\n";
    std::cout << "J2    = " << wgs72.J2n(1).value << " (expected: 0.001082616)\n";

    std::cout << "\n=== Geodesy Module Tests Done ===\n";
    return 0;
}
