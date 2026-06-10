/// test_srp — solar radiation pressure, cannonball + cylindrical shadow (SRP1, R2).
///
/// Verifies src/forces/srp.h per design/derivations/solar_radiation_pressure.md §6:
///
///   1. THE GENERATED CONSTANT — P₁ᴬᵁ = L☉/(4π·AU²·c) from three exact-by-
///      convention values: the implied S(1 AU) = L☉/(4π·AU²) must land on the
///      published TSI basis of the IAU nominal (|S − 1361| < 1 W/m²), and the
///      tracked P must equal an independent double recomputation at machine-eps.
///   2. GEOMETRY IDENTITIES (no numerical reference): inverse-square
///      (|a(AU)|/|a(2AU)| = 4 to round-off); direction exactly anti-sunward;
///      the shadow truth table (sub-solar lit / anti-solar LEO dark with an
///      EXACT-ZERO acceleration / outside the cylinder lit).
///   3. Magnitude sanity: C_R·A/m = 0.026 m²/kg ⇒ |a| ≈ 1.18e-7·(AU/d)² m/s².
///   4. Budget: the 1e-3 TSI representativeness band present and majorizing the
///      |nominal − 1361| offset; precision tightens with wider T.
///   5. Wiring: srp_force in umbra returns an exactly-zero wrench; the ForceFn
///      runs in a DQ Propagator alongside central gravity (finite, perturbs).
///
/// The Sun position inherits TB1's DE430 gating (no new ephemeris claim here).
/// Exit 0 iff every check passes.

#include "astronomy/epoch.h"
#include "constants/constants_provider.h"
#include "dynamics/inertia.h"
#include "dynamics/propagator.h"
#include "dynamics/state.h"
#include "forces/gravity_central.h"
#include "forces/srp.h"
#include "forces/third_body.h"
#include "integrators/runge_kutta.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using boost::multiprecision::cpp_bin_float_50;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

template<typename T>
math::Vector3<T> vec3(double x, double y, double z) {
    return math::Vector3<T>(tv<T>(x), tv<T>(y), tv<T>(z));
}

template<typename T>
double vmag(const math::Vector3<T>& a) {
    return std::sqrt(static_cast<double>(
        (a.x * a.x + a.y * a.y + a.z * a.z).value));
}

const double AU = 149597870700.0;
const double L_SUN = 3.828e26;
const double C_LIGHT = 299792458.0;
const double PI_D = 3.14159265358979323846;

} // namespace

int main() {
    std::cout << std::setprecision(8);
    std::cout << "test_srp: cannonball SRP + cylindrical shadow (SRP1, R2)\n\n";
    using T = double;

    // ---- 1. the generated P_1AU ----
    std::cout << "=== generated constant ===\n";
    math::TrackedValue<T> P = forces::solar_radiation_pressure_1au<T>();
    {
        double S = L_SUN / (4.0 * PI_D * AU * AU);             // implied TSI [W/m^2]
        double P_ref = S / C_LIGHT;                            // independent recompute
        std::cout << "  S(1AU) = " << S << " W/m^2,  P_1AU = " << P.value << " N/m^2\n";
        check("implied TSI lands on the published 1361 W/m^2 basis (|dS| < 1)",
              std::abs(S - 1361.0) < 1.0);
        check("tracked P == independent recompute at machine-eps",
              std::abs(P.value - P_ref) / P_ref < 1e-14);
    }

    // ---- 2. geometry identities ----
    std::cout << "\n=== geometry identities ===\n";
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    const double RE = K.earth.a.value;
    math::TrackedValue<T> cr_am = tv<T>(0.026);   // C_R=1.3, A/m=0.02
    {
        math::Vector3<T> origin = vec3<T>(0.0, 0.0, 0.0);
        math::Vector3<T> s1 = vec3<T>(AU, 0.0, 0.0);
        math::Vector3<T> s2 = vec3<T>(2.0 * AU, 0.0, 0.0);
        math::Vector3<T> a1 = forces::srp_accel<T>(origin, s1, cr_am, K.earth.a);
        math::Vector3<T> a2 = forces::srp_accel<T>(origin, s2, cr_am, K.earth.a);
        double ratio = vmag(a1) / vmag(a2);
        std::cout << "  |a(1AU)|/|a(2AU)| = " << std::setprecision(17) << ratio
                  << std::setprecision(8) << "\n";
        check("inverse-square: ratio == 4 to round-off", std::abs(ratio - 4.0) < 1e-12);
        check("direction exactly anti-sunward (a_x < 0, a_y = a_z = 0)",
              a1.x.value < 0.0 && a1.y.value == 0.0 && a1.z.value == 0.0);
    }
    {
        // Shadow truth table with the Sun along +x.
        math::Vector3<T> s = vec3<T>(AU, 0.0, 0.0);
        double h = RE + 400000.0;
        math::Vector3<T> sub_solar = vec3<T>(h, 0.0, 0.0);
        math::Vector3<T> anti_solar = vec3<T>(-h, 0.0, 0.0);
        math::Vector3<T> outside_cyl = vec3<T>(-h, RE + 1000000.0, 0.0);
        math::Vector3<T> a_lit = forces::srp_accel<T>(sub_solar, s, cr_am, K.earth.a);
        math::Vector3<T> a_dark = forces::srp_accel<T>(anti_solar, s, cr_am, K.earth.a);
        math::Vector3<T> a_out = forces::srp_accel<T>(outside_cyl, s, cr_am, K.earth.a);
        check("sub-solar point is lit (|a| > 0)", vmag(a_lit) > 0.0);
        check("anti-solar LEO point is in umbra (a EXACTLY zero)",
              a_dark.x.value == 0.0 && a_dark.y.value == 0.0 && a_dark.z.value == 0.0);
        check("night-side point outside the cylinder is lit", vmag(a_out) > 0.0);
    }

    // ---- 3. magnitude sanity at the real Sun ----
    std::cout << "\n=== magnitude at the real Sun (J2000 epoch) ===\n";
    forces::ThirdBody<T> sun = forces::sun_third_body<T>();
    astronomy::Epoch<T> base =
        astronomy::Epoch<T>::from_jd(tv<T>(2451545.0), astronomy::TimeScale::TT);
    math::Vector3<T> s_teme = forces::third_body_position_teme<T>(sun, base, tv<T>(0.0));
    {
        double smag = vmag(s_teme);
        double h = RE + 400000.0;
        math::Vector3<T> r_lit(tv<T>(s_teme.x.value / smag * h),
                               tv<T>(s_teme.y.value / smag * h),
                               tv<T>(s_teme.z.value / smag * h));
        math::Vector3<T> a = forces::srp_accel<T>(r_lit, s_teme, cr_am, K.earth.a);
        double amag = vmag(a);
        std::cout << "  d_sun = " << smag / AU << " AU,  |a_srp| = " << amag << " m/s^2\n";
        check("|a_srp| ~ 1.2e-7 m/s^2 for C_R*A/m = 0.026 (1.0e-7..1.4e-7)",
              amag > 1.0e-7 && amag < 1.4e-7);
    }

    // ---- 4. error budget ----
    std::cout << "\n=== error budget ===\n";
    {
        double S = L_SUN / (4.0 * PI_D * AU * AU);
        double nominal_offset = std::abs(S - 1361.0) / 1361.0;
        double band_rel = static_cast<double>(P.errors.accuracy) / P.value;
        std::cout << "  band " << band_rel << " vs nominal-TSI offset " << nominal_offset << "\n";
        check("1e-3 TSI band present and majorizes the nominal-vs-1361 offset",
              band_rel >= 0.99e-3 && band_rel > nominal_offset);

        math::TrackedValue<cpp_bin_float_50> Pb =
            forces::solar_radiation_pressure_1au<cpp_bin_float_50>();
        double pd = static_cast<double>(P.errors.precision);
        double pb = static_cast<double>(Pb.errors.precision);
        std::cout << "  precision: double = " << pd << "  bf50 = " << pb << "\n";
        check("precision > 0 (framework alive)", pd > 0.0 && pb > 0.0);
        check("precision tightens with wider T", pb < pd);
    }

    // ---- 5. wiring ----
    std::cout << "\n=== wiring (wrench + propagator) ===\n";
    {
        // In umbra the full force path must return an exactly-zero wrench.
        double smag = vmag(s_teme);
        double h = RE + 400000.0;
        math::Vector3<T> r_dark(tv<T>(-s_teme.x.value / smag * h),
                                tv<T>(-s_teme.y.value / smag * h),
                                tv<T>(-s_teme.z.value / smag * h));
        math::Vector3<T> zero;
        dynamics::State<T> s_sh = dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), r_dark, zero, zero, tv<T>(0.0));
        dynamics::Wrench<T> w = forces::srp_force<T>(s_sh, K, sun, base, cr_am);
        check("umbra wrench force is exactly zero",
              w.force.x.value == 0.0 && w.force.y.value == 0.0 && w.force.z.value == 0.0);

        // The ForceFn runs in a DQ Propagator alongside central gravity.
        dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(tv<T>(1.0));
        const double r_orbit = RE + 400000.0;
        const double v_orbit = std::sqrt(K.earth.GM.value / r_orbit);
        dynamics::State<T> orb0 = dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), vec3<T>(r_orbit, 0.0, 0.0), zero,
            vec3<T>(0.0, v_orbit, 0.0), tv<T>(0.0));
        dynamics::IntegratorFn<T> integ =
            [](const dynamics::State<T>& y, const math::TrackedValue<T>& dt,
               const integrators::AccelFn<T>& f) { return integrators::runge_kutta_4(y, dt, f); };
        auto central = [](const dynamics::State<T>& st,
                          const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_central(st, KK);
        };
        std::vector<dynamics::ForceFn<T>> fc, fs;
        fc.push_back(central);
        fs.push_back(central);
        fs.push_back(forces::make_srp_force<T>(sun, base, cr_am));
        dynamics::Propagator<T> pc(K, inertia, std::move(fc), integ);
        dynamics::Propagator<T> ps(K, inertia, std::move(fs), integ);
        dynamics::State<T> ec = pc.propagate_to(orb0, tv<T>(600.0), tv<T>(60.0));
        dynamics::State<T> es = ps.propagate_to(orb0, tv<T>(600.0), tv<T>(60.0));
        double diff = std::sqrt(
            std::pow(es.position().x.value - ec.position().x.value, 2) +
            std::pow(es.position().y.value - ec.position().y.value, 2) +
            std::pow(es.position().z.value - ec.position().z.value, 2));
        bool finite = std::isfinite(es.position().x.value) &&
                      std::isfinite(es.position().y.value);
        std::cout << "  central-vs-(central+SRP) position diff over 600 s: " << diff << " m\n";
        check("SRP force usable in the DQ propagator (finite + perturbs)",
              finite && diff > 0.0);
    }

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
