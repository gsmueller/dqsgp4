/// test_geopotential — corrected formula-layer Stage 2 (geopotential unification).
///
/// The unified forces::geopotential_accel_ecef evaluates monopole + zonal +
/// tesseral in ONE Cunningham V/W pass (design/derivations/geopotential.md §3).
/// It changes the ARITHMETIC of the zonal half (Cunningham V/W vs the Legendre-
/// in-u recurrence of gravity_zonal), so it is NOT bit-identical to the summed
/// legacy forces — it is gated at ROUND-OFF against three INDEPENDENT oracles
/// (§5), never sold as bit-exact (no perceived fidelity):
///
///   1. Summed-legacy: geopotential ≈ gravity_central + gravity_zonal +
///      gravity_tesseral, over LEO/MEO/GEO × latitudes, at round-off. The legacy
///      trio uses different arithmetic for the central/zonal halves (independent).
///   2. Closed-form J₂ (Montenbruck-Gill 3.34): the (2,0) term must match the
///      textbook a_J2 — an oracle external to the recursion.
///   3. Monopole identity: the (0,0) term must equal −GM·r/r³ to round-off.
///
/// Plus the TrackedValue contracts: precision tightens with a wider T, and the
/// omitted-degree Kaula tail is booked (> 0) in the accuracy channel.
///
/// ExeGate: nonzero exit code on any failed check.

#include "forces/geopotential.h"
#include "forces/gravity_central.h"
#include "forces/gravity_zonal.h"
#include "forces/gravity_tesseral.h"
#include "constants/gravity_field.h"
#include "constants/zonal_harmonics.h"
#include "constants/tesseral_harmonics.h"
#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_bin_float_50;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0),
        math::TrackedValue<T>::representation_bound(val), T(0));
}

// ECEF point from radius, geocentric latitude, longitude (radians).
template<typename T>
math::Vector3<T> point(double r, double lat, double lon) {
    return math::Vector3<T>(tv<T>(r * std::cos(lat) * std::cos(lon)),
                            tv<T>(r * std::cos(lat) * std::sin(lon)),
                            tv<T>(r * std::sin(lat)));
}

template<typename T>
dynamics::State<T> state_at(const math::Vector3<T>& r_eci) {
    math::Quaternion<T> q = math::Quaternion<T>::identity();
    math::Vector3<T> omega, v(tv<T>(0.0), tv<T>(0.0), tv<T>(0.0));
    return dynamics::State<T>::from_kinematics(q, r_eci, omega, v, tv<T>(0.0));
}

double mag3(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z); }

struct Site { const char* name; double r, lat, lon; };

} // namespace

int main() {
    using T = double;
    const double tol_legacy = 1e-9;   // round-off vs summed legacy (relative)
    const double tol_j2     = 1e-9;   // vs closed-form J2 (relative)
    const double tol_mono   = 1e-12;  // vs -GM r/r3 (relative)

    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    constants::GravityField<T> field = constants::GravityField<T>::egm2008(T(1e-12));
    constants::ZonalHarmonics<T> zh = constants::ZonalHarmonics<T>::egm2008(T(1e-12));
    constants::TesseralHarmonics<T> th = constants::TesseralHarmonics<T>::egm2008();
    const math::TrackedValue<T> mu = K.earth.GM, Re = K.earth.a;
    const double mu_v = mu.value, Re_v = Re.value;
    const int N = field.max_degree();

    const Site sites[] = {
        {"LEO equator",   7.0e6, 0.0,        0.3},
        {"LEO mid-lat",   7.0e6, 0.7853982,  1.1},
        {"LEO near-pole", 7.0e6, 1.4,       -0.5},
        {"MEO",           2.0e7, 0.5,        2.0},
        {"GEO",           4.2164e7, 0.1,    -1.7},
    };

    // ---- Oracle 1: summed-legacy at round-off -----------------------------
    double worst_legacy = 0.0;
    for (const Site& st : sites) {
        math::Vector3<T> r = point<T>(st.r, st.lat, st.lon);
        dynamics::State<T> s = state_at<T>(r);

        math::Vector3<T> au = forces::geopotential_accel_ecef<T>(r, mu, Re, field, N, N);

        dynamics::Wrench<T> wc = forces::gravity_central<T>(s, K);
        dynamics::Wrench<T> wz = forces::gravity_zonal<T>(s, K, zh, N);
        math::Vector3<T> at = forces::gravity_tesseral_accel_ecef<T>(r, mu, Re, th, N);

        double lx = wc.force.x.value + wz.force.x.value + at.x.value;
        double ly = wc.force.y.value + wz.force.y.value + at.y.value;
        double lz = wc.force.z.value + wz.force.z.value + at.z.value;
        double d = mag3(au.x.value - lx, au.y.value - ly, au.z.value - lz);
        double rel = d / mag3(lx, ly, lz);
        if (rel > worst_legacy) worst_legacy = rel;
        std::cout << "    " << st.name << ": |a|=" << mag3(lx, ly, lz)
                  << "  rel.diff=" << rel << "\n";
    }
    check("unified == central+zonal+tesseral at round-off (oracle 1)",
          worst_legacy < tol_legacy);

    // ---- Oracle 2: closed-form J2 (Montenbruck-Gill 3.34) -----------------
    // The (2,0) term must match a_J2 = -(3/2) J2 mu R^2 / r^5 *
    //   [ x(1-5(z/r)^2), y(1-5(z/r)^2), z(3-5(z/r)^2) ], independent of the
    // recursion. Compare geopotential(2,0) to monopole_closed + J2_closed.
    double worst_j2 = 0.0;
    double J2 = zh.Jn(2).value;
    for (const Site& st : sites) {
        math::Vector3<T> r = point<T>(st.r, st.lat, st.lon);
        double x = r.x.value, y = r.y.value, z = r.z.value;
        double r2 = x * x + y * y + z * z, rm = std::sqrt(r2);
        double zr2 = (z * z) / r2;
        double pref = -1.5 * J2 * mu_v * Re_v * Re_v / (r2 * r2 * rm);  // -(3/2)J2 mu R^2/r^5
        double jx = pref * x * (1.0 - 5.0 * zr2);
        double jy = pref * y * (1.0 - 5.0 * zr2);
        double jz = pref * z * (3.0 - 5.0 * zr2);
        double mp = -mu_v / (r2 * rm);                                 // -GM/r^3
        double ox = mp * x + jx, oy = mp * y + jy, oz = mp * z + jz;

        math::Vector3<T> au = forces::geopotential_accel_ecef<T>(r, mu, Re, field, 2, 0);
        double d = mag3(au.x.value - ox, au.y.value - oy, au.z.value - oz);
        double rel = d / mag3(ox, oy, oz);
        if (rel > worst_j2) worst_j2 = rel;
    }
    check("unified (2,0) == monopole + closed-form J2 (oracle 2)",
          worst_j2 < tol_j2);
    std::cout << "    worst J2 rel.diff=" << worst_j2 << "\n";

    // ---- Oracle 3: monopole identity --------------------------------------
    double worst_mono = 0.0;
    for (const Site& st : sites) {
        math::Vector3<T> r = point<T>(st.r, st.lat, st.lon);
        double x = r.x.value, y = r.y.value, z = r.z.value;
        double r2 = x * x + y * y + z * z, rm = std::sqrt(r2);
        double mp = -mu_v / (r2 * rm);
        math::Vector3<T> au = forces::geopotential_accel_ecef<T>(r, mu, Re, field, 0, 0);
        double d = mag3(au.x.value - mp * x, au.y.value - mp * y, au.z.value - mp * z);
        double rel = d / mag3(mp * x, mp * y, mp * z);
        if (rel > worst_mono) worst_mono = rel;
    }
    check("unified (0,0) == -GM r/r^3 (oracle 3)", worst_mono < tol_mono);

    // ---- TrackedValue contracts -------------------------------------------
    // Kaula truncation tail booked (> 0) in accuracy.
    math::Vector3<T> a_acc = forces::geopotential_accel_ecef<T>(
        point<T>(7.0e6, 0.4, 0.9), mu, Re, field, N, N);
    check("Kaula truncation tail tracked (> 0) in accuracy",
          a_acc.x.errors.accuracy > 0.0 && a_acc.z.errors.accuracy > 0.0);

    // Precision tightens with a wider numeric type T (the calling card).
    double pd = a_acc.x.errors.precision;
    constants::ConstantsProvider<cpp_bin_float_50> Kb =
        constants::ConstantsProvider<cpp_bin_float_50>::wgs84(cpp_bin_float_50("1e-30"));
    constants::GravityField<cpp_bin_float_50> fb =
        constants::GravityField<cpp_bin_float_50>::egm2008(cpp_bin_float_50("1e-30"));
    math::Vector3<cpp_bin_float_50> ab = forces::geopotential_accel_ecef<cpp_bin_float_50>(
        point<cpp_bin_float_50>(7.0e6, 0.4, 0.9), Kb.earth.GM, Kb.earth.a, fb,
        fb.max_degree(), fb.max_degree());
    double pb = static_cast<double>(ab.x.errors.precision);
    // The precision is dominated by the zonal Jₙ generator's series tolerance
    // (1e-12 at double, 1e-30 at bin_float) — convergent math that tightens with a
    // wider T by ~the tolerance ratio. So precision tightens by many orders (here
    // ~1e-16), the calling card, though not the representation-only ~1e-34 of the
    // model_coefficient-only tesseral path.
    std::cout << "    precision double=" << pd << "  bin_float=" << pb << "\n";
    check("geopotential precision tracked (> 0)", pd > 0.0 && pb > 0.0);
    check("geopotential precision tightens with wider T", pb < pd * 1e-10);

    std::cout << "\n  geopotential: " << passed << " passed, " << failed
              << " failed\n";
    std::cout << "    worst legacy rel.diff=" << worst_legacy
              << "  worst mono rel.diff=" << worst_mono << "\n";
    return failed == 0 ? 0 : 1;
}
