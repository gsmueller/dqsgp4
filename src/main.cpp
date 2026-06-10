// Test driver for the geodesy library.
//
// Constructs a WGS84 EquipotentialEllipsoid from the four defining parameters,
// then prints all derived constants with their three error bounds.
// Verifies against published values from NGA.STND.0036 Appendix B.

#include <iostream>
#include <iomanip>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>

#include "math/tracked_value.h"
#include "geodesy/equipotential_ellipsoid.h"
#include "tle/tle_parser.h"

using T = boost::multiprecision::cpp_bin_float_50;

template<typename ValT>
void print_value(const char* name, const math::TrackedValue<ValT>& tv,
                 const char* published = nullptr)
{
    std::cout << std::setprecision(15) << std::scientific;
    std::cout << name << ":\n";
    std::cout << "  value       = " << tv.value << "\n";
    std::cout << "  measurement = " << tv.errors.measurement << "\n";
    std::cout << "  precision   = " << tv.errors.precision << "\n";
    std::cout << "  accuracy    = " << tv.errors.accuracy << "\n";
    std::cout << "  total error = " << tv.total_error() << "\n";
    std::cout << "  reliable digits = " << tv.reliable_digits() << "\n";
    if (published) {
        std::cout << "  published   = " << published << "\n";
    }
    std::cout << "\n";
}

int main()
{
    std::cout << "=== WGS84 Equipotential Ellipsoid Test ===\n\n";

    // Four defining parameters as TrackedValue
    // a and 1/f are definitional: measurement error = 0
    // GM is measured: measurement error = +/- 8e6 m^3/s^2
    // omega is definitional: measurement error = 0

    math::TrackedValue<T> a     = math::TrackedValue<T>::defined("6378137.0");
    math::TrackedValue<T> inv_f = math::TrackedValue<T>::defined("298.257223563");
    math::TrackedValue<T> GM    = math::TrackedValue<T>::measured("3.986004418e14", "8e6");
    math::TrackedValue<T> omega = math::TrackedValue<T>::defined("7.292115e-5");

    // Series tolerance: evaluate until truncation error < 1e-45
    // (well within 50-digit precision of cpp_bin_float_50)
    T tolerance("1e-45");

    std::cout << "Constructing ellipsoid from:\n";
    print_value("a", a);
    print_value("1/f", inv_f);
    print_value("GM", GM);
    print_value("omega", omega);

    std::cout << "Series tolerance: " << tolerance << "\n\n";
    std::cout << "Computing derived constants...\n\n";

    geodesy::EquipotentialEllipsoid<T> wgs84(a, inv_f, GM, omega, tolerance);

    // Print all derived geometric constants (Derivation 003)
    std::cout << "=== Geometric Constants (measurement error = 0) ===\n\n";

    print_value("f (flattening)", wgs84.f,
                "3.3528106647475e-03");
    print_value("e^2 (first eccentricity squared)", wgs84.e2,
                "6.694379990141e-03");
    print_value("e (first eccentricity)", wgs84.e,
                "8.1819190842622e-02");
    print_value("e' (second eccentricity)", wgs84.e_prime,
                "8.2094437949696e-02");
    print_value("b (semi-minor axis)", wgs84.b,
                "6356752.3142 m");
    print_value("E (linear eccentricity)", wgs84.E_lin,
                "5.2185400842339e+05 m");
    print_value("c (polar radius of curvature)", wgs84.c,
                "6399593.6258 m");

    // Series-evaluated terms (Derivations 001-002)
    std::cout << "=== Series-Evaluated Terms ===\n\n";

    print_value("q0", wgs84.q0,
                "7.334625787083e-05");
    print_value("q0'", wgs84.q0p,
                "2.688041300461e-03");

    // Physical constants (Derivation 004) — measurement error propagates from GM
    std::cout << "=== Physical Constants (measurement error from GM) ===\n\n";

    print_value("m (gravity formula constant)", wgs84.m_const,
                "3.449786506841e-03");
    print_value("gamma_e (equator gravity)", wgs84.gamma_e,
                "9.7803253359 m/s^2");
    print_value("gamma_p (pole gravity)", wgs84.gamma_p,
                "9.8321849379 m/s^2");
    print_value("k (Somigliana constant)", wgs84.k_som,
                "1.931852652458e-03");
    print_value("U0 (normal potential)", wgs84.U0,
                "6.26368517146e+07 m^2/s^2");

    // Radii (Derivation 005)
    std::cout << "=== Radii ===\n\n";

    print_value("R1 (mean radius)", wgs84.R1,
                "6371008.7714 m");
    print_value("R3 (equal volume sphere)", wgs84.R3,
                "6371000.7900 m");

    // Zonal harmonics (Derivation 006)
    std::cout << "=== Zonal Harmonics ===\n\n";

    print_value("J2", wgs84.J2n(1),
                "1.082629821313e-03");
    print_value("J4", wgs84.J2n(2),
                "-2.37091222e-06");
    print_value("J6", wgs84.J2n(3),
                "6.08347e-09");

    // Normal gravity at sample latitudes
    std::cout << "=== Normal Gravity ===\n\n";

    math::TrackedValue<T> phi_0 = math::TrackedValue<T>::exact_integer(0);
    // pi/4 rad: GENERATED (boost quarter_pi to full T), not stamped as a decimal
    // -- a transcendental is exact-by-definition but its decimal is infinite, so
    // defined("0.7853981633974483") was a 16-digit truncation (wrong beyond
    // double, precision frozen at the decimal). Matches the phi_90 form below.
    math::TrackedValue<T> phi_45 = math::TrackedValue<T>(
        boost::math::constants::quarter_pi<T>(), T(0),
        boost::math::constants::quarter_pi<T>() * std::numeric_limits<T>::epsilon(), T(0));
    math::TrackedValue<T> phi_90 = math::TrackedValue<T>(
        boost::math::constants::half_pi<T>(), T(0),
        boost::math::constants::half_pi<T>() * std::numeric_limits<T>::epsilon(), T(0));

    print_value("gamma(0 deg) [should = gamma_e]", wgs84.normal_gravity(phi_0),
                "9.7803253359 m/s^2");
    print_value("gamma(45 deg)", wgs84.normal_gravity(phi_45),
                "9.8061999 m/s^2");
    print_value("gamma(90 deg) [should = gamma_p]", wgs84.normal_gravity(phi_90),
                "9.8321849379 m/s^2");

    // === TLE Parsing Test ===
    std::cout << "=== TLE Parsing Test ===\n\n";

    // ISS (ZARYA) — a well-known test case
    std::string tle_line1 = "1 25544U 98067A   20045.18587073  .00000950  00000-0  25302-4 0  9990";
    std::string tle_line2 = "2 25544  51.6443 242.7440 0004885 264.6060 219.4814 15.49165514212180";

    tle::TleData td;
    if (tle::parse("ISS (ZARYA)", tle_line1, tle_line2, td)) {
        std::cout << "Parsed TLE for: " << td.name << "\n";
        std::cout << "  Satellite ID:      " << td.satellite_id << "\n";
        std::cout << "  Classification:    " << td.classification << "\n";
        std::cout << "  Intl Designator:   " << td.intl_designator << "\n";
        std::cout << "  Epoch Year:        " << td.epoch_year << "\n";
        std::cout << "  Epoch Day:         " << std::setprecision(11) << td.epoch_day << "\n";
        std::cout << "  Mean Motion dt/2:  " << td.mean_motion_dt2 << "\n";
        std::cout << "  Mean Motion ddt/6: " << td.mean_motion_ddt6 << "\n";
        std::cout << "  B*:                " << td.bstar << "\n";
        std::cout << "  Ephemeris Type:    " << td.ephemeris_type << "\n";
        std::cout << "  Inclination:       " << td.inclination_deg << " deg\n";
        std::cout << "  RAAN:              " << td.raan_deg << " deg\n";
        std::cout << "  Eccentricity:      " << std::setprecision(8) << td.eccentricity << "\n";
        std::cout << "  Arg Perigee:       " << td.arg_perigee_deg << " deg\n";
        std::cout << "  Mean Anomaly:      " << td.mean_anomaly_deg << " deg\n";
        std::cout << "  Mean Motion:       " << std::setprecision(11) << td.mean_motion_rev_day << " rev/day\n";
        std::cout << "  Revolution #:      " << td.revolution_number << "\n";

        // Convert to TrackedValue elements
        tle::TleElements<T> elements = tle::TleElements<T>::from_tle_data(td);
        std::cout << "\n  TrackedValue conversions:\n";
        print_value("  inclination [rad]", elements.inclination);
        print_value("  eccentricity", elements.eccentricity);
        print_value("  mean motion [rad/min]", elements.mean_motion);
        print_value("  B*", elements.bstar);
    } else {
        std::cout << "ERROR: Failed to parse TLE\n";
    }

    std::cout << "\n=== Done ===\n";
    return 0;
}
