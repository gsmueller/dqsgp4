/// test_tesseral — D2 (DQSGP4 Completion Roadmap, register D2).
///
/// Verifies the tesseral/sectoral (longitude-dependent) geopotential:
///   - the EGM2008 C̄_22/S̄_22 denormalize to the expected unnormalized values;
///   - the sectoral acceleration is LONGITUDE-DEPENDENT — the defining tesseral
///     signature (a zonal field would be identical at every longitude);
///   - the GMST coupling is real: the SAME inertial position feels a different
///     tesseral pull as the Earth rotates beneath it (gmst = 0 vs gmst = 90°);
///   - the GMST rotation is the identity at gmst = 0;
///   - precision tightens with a wider numeric type T.
///
/// ExeGate D2: nonzero exit code on any failed check.

#include "forces/gravity_tesseral.h"
#include "constants/tesseral_harmonics.h"
#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>
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
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

template<typename T>
math::Vector3<T> eq_point(double r, double lambda) {  // equatorial point at longitude λ
    return math::Vector3<T>(tv<T>(r * std::cos(lambda)), tv<T>(r * std::sin(lambda)), tv<T>(0.0));
}

double vmag(const math::Vector3<double>& a) {
    return std::sqrt(a.x.value * a.x.value + a.y.value * a.y.value + a.z.value * a.z.value);
}

double vdiff(const math::Vector3<double>& a, const math::Vector3<double>& b) {
    double dx = a.x.value - b.x.value, dy = a.y.value - b.y.value, dz = a.z.value - b.z.value;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void test_denormalization() {
    using T = double;
    auto th = constants::TesseralHarmonics<T>::egm2008();
    math::TrackedValue<T> C, S;
    th.get(2, 2, C, S);
    // N_22 = sqrt(5/12) = 0.6454972; C_22 = N_22 * 2.43938e-6 = 1.5746e-6.
    double N22 = std::sqrt(5.0 / 12.0);
    check("C22 denormalized = N22 * Cbar22", std::abs(C.value - N22 * 2.43938e-6) < 1e-13);
    check("S22 denormalized = N22 * Sbar22", std::abs(S.value - N22 * (-1.40027e-6)) < 1e-13);
    check("tesseral table degree = 2", th.max_degree() == 2);
}

void test_longitude_dependence() {
    using T = double;
    auto th = constants::TesseralHarmonics<T>::egm2008();
    auto K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    const double r = 7.0e6;

    math::Vector3<T> a0  = forces::gravity_tesseral_accel_ecef<T>(
        eq_point<T>(r, 0.0), K.earth.GM, K.earth.a, th, 2);
    math::Vector3<T> a90 = forces::gravity_tesseral_accel_ecef<T>(
        eq_point<T>(r, boost::math::constants::pi<T>() / 2.0), K.earth.GM, K.earth.a, th, 2);

    check("tesseral acceleration is nonzero", vmag(a0) > 1e-9);
    // The sectoral C22 term has a cos(2λ)/sin(2λ) structure: at λ=0 vs λ=90° the
    // 2λ argument flips by 180°, so the acceleration must differ markedly. A
    // purely zonal field would give identical accelerations at both longitudes.
    check("sectoral acceleration is LONGITUDE-dependent (D2 core)",
          vdiff(a0, a90) > 1e-8);
    std::cout << "    |a(λ=0)|=" << vmag(a0) << "  |a(λ=90°)|=" << vmag(a90)
              << "  |Δ|=" << vdiff(a0, a90) << " m/s²\n";
}

void test_gmst_rotation() {
    using T = double;
    auto th = constants::TesseralHarmonics<T>::egm2008();
    auto K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    const double r = 7.0e6;

    math::Vector3<T> r_eci = eq_point<T>(r, 0.0);  // (r, 0, 0)
    math::Quaternion<T> q = math::Quaternion<T>::identity();
    math::Vector3<T> omega, v(tv<T>(0.0), tv<T>(7500.0), tv<T>(0.0));
    dynamics::State<T> s = dynamics::State<T>::from_kinematics(q, r_eci, omega, v, tv<T>(0.0));

    // gmst = 0: ECEF == ECI, identity pose => body force == direct ECEF accel.
    dynamics::Wrench<T> w0 = forces::gravity_tesseral<T>(s, K, th, 2, tv<T>(0.0));
    math::Vector3<T> direct = forces::gravity_tesseral_accel_ecef<T>(r_eci, K.earth.GM, K.earth.a, th, 2);
    check("gmst=0 rotation is identity",
          std::abs(w0.force.x.value - direct.x.value) < 1e-18
          && std::abs(w0.force.y.value - direct.y.value) < 1e-18
          && std::abs(w0.force.z.value - direct.z.value) < 1e-18);

    // GMST coupling: the SAME ECI position feels a different tesseral pull as
    // the Earth rotates 90° beneath it.
    dynamics::Wrench<T> w90 = forces::gravity_tesseral<T>(
        s, K, th, 2, tv<T>(boost::math::constants::pi<T>() / 2.0));
    check("Earth rotation (GMST) changes the inertial tesseral pull",
          vdiff(w0.force, w90.force) > 1e-8);
}

void test_precision_and_accuracy() {
    // C̄_22/S̄_22 are gravity-MODEL coefficients encoded with model_coefficient
    // (2026-06-05 panel ruling, uniform by-nature categorization): their binary
    // STORAGE is a computational precision that TIGHTENS with a wider T, while
    // the finite EGM2008 source digits are a model-fidelity ACCURACY floor that
    // is T-INDEPENDENT. So the tesseral acceleration's precision scales with T
    // (the calling card reaches the geopotential too), and accuracy carries the
    // honest source floor.
    auto thd = constants::TesseralHarmonics<double>::egm2008();
    auto Kd = constants::ConstantsProvider<double>::wgs84(1e-12);
    auto ad = forces::gravity_tesseral_accel_ecef<double>(
        eq_point<double>(7.0e6, 0.3), Kd.earth.GM, Kd.earth.a, thd, 2).x;
    double pd = ad.errors.precision, accd = ad.errors.accuracy;

    auto thb = constants::TesseralHarmonics<cpp_bin_float_50>::egm2008();
    auto Kb = constants::ConstantsProvider<cpp_bin_float_50>::wgs84(cpp_bin_float_50("1e-12"));
    auto ab = forces::gravity_tesseral_accel_ecef<cpp_bin_float_50>(
        eq_point<cpp_bin_float_50>(7.0e6, 0.3), Kb.earth.GM, Kb.earth.a, thb, 2).x;
    double pb = static_cast<double>(ab.errors.precision);
    double accb = static_cast<double>(ab.errors.accuracy);

    check("tesseral precision budget tracked (> 0)", pd > 0.0 && pb > 0.0);
    check("tesseral precision tightens with wider T (computational)", pb < pd * 1e-20);
    check("tesseral accuracy budget tracked (> 0, the model-digit floor)", accd > 0.0 && accb > 0.0);
    check("tesseral accuracy source-floored by EGM2008 digits (T-independent)",
          std::abs(accd - accb) / accd < 1e-6);
}

} // namespace

int main() {
    test_denormalization();
    test_longitude_dependence();
    test_gmst_rotation();
    test_precision_and_accuracy();

    std::cout << "\n  tesseral: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
