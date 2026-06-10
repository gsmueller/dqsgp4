/// @file test_tle/main.cpp
/// @brief Tests for the tle/ module: TLE parsing.

#include <iostream>
#include <iomanip>
#include "math/tracked_value.h"
#include "tle/tle_parser.h"

using T = double;

int main() {
    std::cout << "=== TLE Module Tests ===\n\n";

    // ISS TLE
    std::string line1 = "1 25544U 98067A   20045.18587073  .00000950  00000-0  25302-4 0  9990";
    std::string line2 = "2 25544  51.6443 242.7440 0004885 264.6060 219.4814 15.49165514212180";

    tle::TleData td;
    if (tle::parse("ISS (ZARYA)", line1, line2, td)) {
        std::cout << "PASS: Parsed ISS TLE\n";
        std::cout << "  sat_id: " << td.satellite_id << " (expected: 25544)\n";
        std::cout << "  incl: " << td.inclination_deg << " (expected: 51.6443)\n";
        std::cout << "  ecc: " << std::setprecision(7) << td.eccentricity << " (expected: 0.0004885)\n";
        std::cout << "  n: " << std::setprecision(11) << td.mean_motion_rev_day << " (expected: 15.49165514)\n";

        tle::TleElements<T> elems = tle::TleElements<T>::from_tle_data(td);
        std::cout << "  incl [rad]: " << elems.inclination.value << "\n";
        std::cout << "  incl measurement error: " << elems.inclination.errors.measurement << "\n";
        std::cout << "  incl reliable digits: " << elems.inclination.reliable_digits() << "\n";
    } else {
        std::cout << "FAIL: Could not parse ISS TLE\n";
    }

    // Test bad input
    tle::TleData td_bad;
    if (!tle::parse("X", "bad", td_bad)) {
        std::cout << "PASS: Rejected bad TLE\n";
    }

    std::cout << "\n=== TLE Module Tests Done ===\n";
    return 0;
}
