/// test_drag_density — DRAG1 (DQSGP4 Completion Roadmap, register DRAG1).
///
/// Verifies the pluggable atmospheric-density model added to forces/drag.h:
///   - DensityModel<T> is a real seam: a custom callback drives make_drag;
///   - the exponential (Lane) model computes rho_0 * exp(-(h-h_0)/H_scale);
///   - the drag-acceleration accuracy is MODEL-DERIVED — it propagates from the
///     density model's declared fractional accuracy, so doubling the model's
///     rel_accuracy doubles the drag accuracy (it is no longer a flat 30%
///     stamped on the output);
///   - the NRLMSISE-00 hook stub returns a callable DensityModel (documented
///     fallback, wider 50% band);
///   - density precision tightens with a wider numeric type T;
///   - make_drag_exponential == make_drag(exponential_density_model(...)).
///
/// ExeGate DRAG1: nonzero exit code on any failed check.

#include "forces/drag.h"
#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "dynamics/wrench.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <functional>
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

// A LEO state at 400 km, identity orientation (body frame == world frame).
template<typename T>
dynamics::State<T> leo_state(const constants::ConstantsProvider<T>& K) {
    const T r_mag = K.earth.a.value + T(400000.0);
    const T comp = r_mag / sqrt(T(3));
    math::Vector3<T> r(tv<T>(static_cast<double>(comp)),
                       tv<T>(static_cast<double>(comp)),
                       tv<T>(static_cast<double>(comp)));
    math::Vector3<T> v(tv<T>(0.0), tv<T>(7600.0), tv<T>(1000.0));
    math::Vector3<T> omega_body;
    math::Quaternion<T> q = math::Quaternion<T>::identity();
    return dynamics::State<T>::from_kinematics(q, r, omega_body, v, tv<T>(0.0));
}

void test_exponential_value() {
    using T = double;
    forces::DensityModel<T> model =
        forces::exponential_density_model<T>(tv<T>(2.789e-10), tv<T>(200000.0), tv<T>(50000.0));
    // At h = 300 km: rho = 2.789e-10 * exp(-(300000-200000)/50000) = 2.789e-10*e^-2.
    math::TrackedValue<T> rho = model(tv<T>(300000.0));
    double expect = 2.789e-10 * std::exp(-2.0);
    check("exponential density value", std::abs(rho.value - expect) <= 1e-24);
    check("exponential density accuracy = 30% (default)",
          std::abs(rho.errors.accuracy - 0.30 * rho.value) <= 1e-24);
}

void test_pluggable_callback() {
    using T = double;
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    dynamics::State<T> s = leo_state<T>(K);

    // A custom density model: constant 1e-11 kg/m^3 regardless of altitude.
    forces::DensityModel<T> constant_rho =
        [](const math::TrackedValue<T>&) { return tv<T>(1e-11); };
    auto drag = forces::make_drag<T>(constant_rho, tv<T>(0.5));
    dynamics::Wrench<T> w = drag(s, K);
    double fmag = std::sqrt(w.force.x.value * w.force.x.value
                          + w.force.y.value * w.force.y.value
                          + w.force.z.value * w.force.z.value);
    check("custom density model drives make_drag (nonzero drag)", fmag > 0.0);

    // Zero-density model => zero drag (the callback genuinely controls density).
    forces::DensityModel<T> zero_rho =
        [](const math::TrackedValue<T>&) { return tv<T>(0.0); };
    dynamics::Wrench<T> w0 = forces::make_drag<T>(zero_rho, tv<T>(0.5))(s, K);
    check("zero-density model => zero drag",
          w0.force.x.value == 0.0 && w0.force.y.value == 0.0 && w0.force.z.value == 0.0);
}

void test_model_derived_accuracy() {
    using T = double;
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    dynamics::State<T> s = leo_state<T>(K);

    // Two exponential models identical except for declared accuracy: 10% vs 30%.
    auto drag10 = forces::make_drag<T>(
        forces::exponential_density_model<T>(tv<T>(2.789e-10), tv<T>(200000.0), tv<T>(50000.0), 0.10),
        tv<T>(0.5));
    auto drag30 = forces::make_drag<T>(
        forces::exponential_density_model<T>(tv<T>(2.789e-10), tv<T>(200000.0), tv<T>(50000.0), 0.30),
        tv<T>(0.5));
    dynamics::Wrench<T> w10 = drag10(s, K);
    dynamics::Wrench<T> w30 = drag30(s, K);

    // The force VALUE is identical (same density value)...
    check("force value independent of declared accuracy",
          std::abs(w10.force.x.value - w30.force.x.value) <= 1e-30);
    // ...but the accuracy is MODEL-DERIVED: 3x the declared accuracy => 3x the
    // drag accuracy. (A flat hardcoded bound could not do this.)
    check("drag accuracy > 0 (model-derived)", w10.force.x.errors.accuracy > 0.0);
    double ratio = w30.force.x.errors.accuracy / w10.force.x.errors.accuracy;
    check("drag accuracy scales with the density model's declared accuracy (3x)",
          std::abs(ratio - 3.0) < 1e-6);
}

void test_nrlmsise_stub() {
    using T = double;
    forces::DensityModel<T> msis =
        forces::nrlmsise00_density_model_stub<T>(tv<T>(2.789e-10), tv<T>(200000.0), tv<T>(50000.0));
    math::TrackedValue<T> rho = msis(tv<T>(250000.0));
    check("NRLMSISE-00 stub is callable & finite", std::isfinite(rho.value) && rho.value > 0.0);
    // Documented fallback: a wider 50% accuracy band.
    check("NRLMSISE-00 stub tags a 50% accuracy band",
          std::abs(rho.errors.accuracy - 0.50 * rho.value) <= 1e-24);
}

void test_precision_scales() {
    // Density precision tightens with a wider T (the budget propagates through
    // the tracked exp; the model is computed in T, not raw double).
    forces::DensityModel<double> md =
        forces::exponential_density_model<double>(tv<double>(2.789e-10), tv<double>(200000.0), tv<double>(50000.0));
    double pd = md(tv<double>(300000.0)).errors.precision;

    forces::DensityModel<cpp_bin_float_50> mb =
        forces::exponential_density_model<cpp_bin_float_50>(
            tv<cpp_bin_float_50>(2.789e-10), tv<cpp_bin_float_50>(200000.0), tv<cpp_bin_float_50>(50000.0));
    double pb = static_cast<double>(mb(tv<cpp_bin_float_50>(300000.0)).errors.precision);

    check("density precision > 0", pd > 0.0 && pb > 0.0);
    check("density precision tightens with wider T", pb < pd * 1e-20);
}

void test_convenience_equivalence() {
    using T = double;
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    dynamics::State<T> s = leo_state<T>(K);

    auto via_conv = forces::make_drag_exponential<T>(
        tv<T>(2.789e-10), tv<T>(200000.0), tv<T>(50000.0), tv<T>(0.5));
    auto via_generic = forces::make_drag<T>(
        forces::exponential_density_model<T>(tv<T>(2.789e-10), tv<T>(200000.0), tv<T>(50000.0)),
        tv<T>(0.5));
    dynamics::Wrench<T> a = via_conv(s, K);
    dynamics::Wrench<T> b = via_generic(s, K);
    check("make_drag_exponential == make_drag(exponential_density_model)",
          a.force.x.value == b.force.x.value && a.force.y.value == b.force.y.value
          && a.force.z.value == b.force.z.value
          && std::abs(a.force.x.errors.accuracy - b.force.x.errors.accuracy) <= 1e-30);
}

} // namespace

int main() {
    test_exponential_value();
    test_pluggable_callback();
    test_model_derived_accuracy();
    test_nrlmsise_stub();
    test_precision_scales();
    test_convenience_equivalence();

    std::cout << "\n  drag density: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
