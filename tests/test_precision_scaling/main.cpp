/// test_precision_scaling — selectable precision (OBJ-2 / REQ-SY-7).
///
/// The library is template-parameterized on the scalar type T, so the SAME
/// code runs at `double` (~16 digits) or `cpp_bin_float_50` (~50 digits) with
/// no source changes, and the reported error bounds tighten accordingly. This
/// test runs identical computations at both precisions and verifies:
///   - the values agree to double precision (cpp_bin_float_50 is the truncation
///     target of double);
///   - the cpp_bin_float_50 reported total_error is dramatically smaller than
///     double's (the precision bound scales with the type's machine epsilon);
///   - a real force-model call (gravity_central) instantiates and runs at
///     cpp_bin_float_50 with the same tightening.
///
/// Exit 0 iff every check passes.

#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "forces/gravity_central.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>

using D = double;
using HP = boost::multiprecision::cpp_bin_float_50;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

/// A representative tracked computation: |(r, r, r)| = sqrt(3) r, with r given
/// its binary-representation precision. Returns the result tracked value.
template<typename T>
math::TrackedValue<T> magnitude_chain() {
    T r = T(7000000);  // ~LEO radius scale [m]
    math::TrackedValue<T> x(r, T(0), math::TrackedValue<T>::representation_bound(r), T(0));
    math::TrackedValue<T> sumsq = x * x + x * x + x * x;
    return sqrt(sumsq);
}

/// gravity_central force x-component at a fixed LEO position, tracked.
template<typename T>
math::TrackedValue<T> gravity_fx(const T& tol) {
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(tol);
    T c = (K.earth.a.value + T(400000)) / sqrt(T(3));
    auto tvv = [](const T& v) {
        return math::TrackedValue<T>(v, T(0), math::TrackedValue<T>::representation_bound(v), T(0));
    };
    math::Vector3<T> r(tvv(c), tvv(c), tvv(c));
    math::Vector3<T> v(tvv(T(0)), tvv(T(0)), tvv(T(0)));
    math::Vector3<T> w0;
    dynamics::State<T> s = dynamics::State<T>::from_kinematics(
        math::Quaternion<T>::identity(), r, w0, v, tvv(T(0)));
    return forces::gravity_central(s, K).force.x;
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(6);
    std::cout << "test_precision_scaling: double vs cpp_bin_float_50 (REQ-SY-7)\n\n";

    // ---- tracked magnitude chain ----
    std::cout << "=== magnitude chain |(r,r,r)| = sqrt(3) r ===\n";
    math::TrackedValue<D> md = magnitude_chain<D>();
    math::TrackedValue<HP> mh = magnitude_chain<HP>();
    const double err_d = md.total_error();
    const double err_h = static_cast<double>(mh.total_error());
    const double val_gap = std::abs(static_cast<double>(mh.value) - md.value);
    std::cout << "  value            : double=" << md.value
              << "  hp=" << static_cast<double>(mh.value) << "\n";
    std::cout << "  total_error      : double=" << err_d << "  hp=" << err_h << "\n";
    std::cout << "  |value gap|      : " << val_gap << "\n";
    check("values agree to double precision", val_gap < 1e-6);
    check("double error bound is nonzero", err_d > 0.0);
    check("cpp_bin_float_50 error bound is far smaller (< 1e-20 of double's)",
          err_h < 1e-20 * err_d);

    // ---- gravity_central at both precisions ----
    std::cout << "\n=== gravity_central force.x at LEO ===\n";
    math::TrackedValue<D> gd = gravity_fx<D>(1e-12);
    math::TrackedValue<HP> gh = gravity_fx<HP>(HP("1e-40"));
    const double gerr_d = gd.errors.precision;
    const double gerr_h = static_cast<double>(gh.errors.precision);
    const double gval_gap = std::abs(static_cast<double>(gh.value) - gd.value);
    std::cout << "  force.x          : double=" << gd.value
              << "  hp=" << static_cast<double>(gh.value) << "\n";
    std::cout << "  precision bound  : double=" << gerr_d << "  hp=" << gerr_h << "\n";
    check("gravity_central runs at cpp_bin_float_50 and agrees with double",
          gval_gap < 1e-9 * std::abs(gd.value));
    check("cpp_bin_float_50 force precision is far tighter (< 1e-20 of double's)",
          gerr_h < 1e-20 * gerr_d && gerr_d > 0.0);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
