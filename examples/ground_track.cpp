/// examples/ground_track.cpp — a satellite ground track through the full chain
/// (gate EX1, so this example can never rot):
///
///   TLE  →  DqSgp4Propagator (boosted + lunisolar preset)  →  TEME state
///        →  GCRS (inverse IAU2006 precession)  →  ITRS (gcrs_to_itrs: nutation
///           + GAST + polar motion)  →  geocentric subpoint (lat, lon, alt).
///
/// Uses ONLY the umbrella header. Honesty notes, matching the library's grades:
///  - TEME→GCRS here is the inverse precession (the dominant term); omitted
///    nutation/equinox ≈ 24″ ≈ 0.8 km on the ground — fine for a ground track.
///  - The TLE epoch (UTC) is used for TT and UT1 alike: |TT−UTC| ≈ 69 s shifts
///    the SUBPOINT longitude by ~0.29° (Earth rotates 15″/s); pass real ΔAT/ΔUT1
///    for timing-grade tracks. Polar motion x_p = y_p = 0 here (≤ ~0.3″ ≈ 9 m).
///  - The subpoint is GEOCENTRIC (lat = asin(z/r)); geodetic conversion adds up
///    to ~0.19° at mid-latitudes.
///
/// Build (from the repo root, after vcvars64):
///   cl /std:c++20 /EHsc /O2 /I src /I vcpkg_installed\x64-windows\include ^
///      examples\ground_track.cpp src\tle\tle_parser.cpp /Fe:ground_track.exe

#include "dqsgp4.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using T = double;
using TV = math::TrackedValue<T>;

int main() {
    std::cout << "ground_track: TLE -> DQSGP4 -> ITRS -> subpoint\n\n";
    const T tol = T(1e-12);
    const double PI_D = 3.14159265358979323846;

    // SGP4-VER satellite 00005 (i = 34.27 deg, e = 0.186).
    std::string l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    std::string l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tle::TleData td;
    if (!tle::parse(l1, l2, td)) {
        std::cerr << "TLE parse failed\n";
        return 1;
    }

    // Boosted gravity + the lunisolar third-body preset, one line each.
    dynamics::DqForceOptions<T> opt;
    opt.lunisolar = true;
    auto prop = dynamics::DqSgp4Propagator<T>::from_tle(
        td, tol, dynamics::PropagatorMode::Boosted, math::exact<T>(30), opt);

    // The TLE epoch as the absolute time base (UTC; see the honesty note above).
    tle::TleElements<T> el = tle::TleElements<T>::from_tle_data(td);
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(tol);
    const double Re = K.earth.a.value;

    std::cout << "  t [min]   lat [deg]   lon [deg]   alt [km]\n";
    int bad = 0;
    for (int minutes = 0; minutes <= 120; minutes += 15) {
        // [guide:itrs]
        dynamics::State<T> s = prop.propagate(math::exact<T>(minutes));

        // Absolute epoch of this sample (two-part; day part preserved).
        TV jd = el.epoch_jd + math::exact<T>(minutes) / math::exact<T>(1440);
        TV t_tt = astronomy::centuries_since_j2000<T>(jd);

        // TEME -> GCRS by the inverse precession (dominant term; note above).
        math::Vector3<T> r_teme = s.position();
        math::Vector3<T> r_gcrs = astronomy::precession_iau2006<T>(t_tt).transpose() * r_teme;

        // GCRS -> ITRS: nutation + GAST + polar motion (x_p = y_p = 0 demo EOP).
        math::Matrix3<T> M = astronomy::gcrs_to_itrs<T>(
            t_tt, jd, math::exact<T>(0), math::exact<T>(0));
        math::Vector3<T> r_itrs = M * r_gcrs;
        // [guide:end]

        const double x = r_itrs.x.value, y = r_itrs.y.value, z = r_itrs.z.value;
        const double r = std::sqrt(x * x + y * y + z * z);
        const double lat = std::asin(z / r) * 180.0 / PI_D;     // geocentric
        const double lon = std::atan2(y, x) * 180.0 / PI_D;
        const double alt_km = (r - Re) / 1000.0;

        std::cout << std::setw(7) << minutes << "   " << std::fixed
                  << std::setprecision(4) << std::setw(9) << lat << "   "
                  << std::setw(9) << lon << "   " << std::setprecision(1)
                  << std::setw(8) << alt_km << "\n" << std::defaultfloat;

        // Physical gates (EX1): |geocentric lat| <= inclination; altitude within
        // the orbit's perigee/apogee band; everything finite.
        if (!(std::isfinite(lat) && std::isfinite(lon) && std::isfinite(alt_km))) ++bad;
        if (std::abs(lat) > 34.5) ++bad;                 // i = 34.2682 deg
        if (alt_km < 400.0 || alt_km > 4500.0) ++bad;    // perigee ~640, apogee ~3850
    }

    if (bad > 0) {
        std::cerr << "\nFAIL — " << bad << " sample(s) violated the physical bounds\n";
        return 1;
    }
    std::cout << "\nPASS — all samples within the orbit's physical bounds\n";
    return 0;
}
