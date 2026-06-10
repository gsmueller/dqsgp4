/// test_integrator_order — empirical convergence order of the integrators
/// (REQ-IN-4: "local truncation error is O(dt^p) where p is the claimed
/// order").
///
/// A circular Kepler orbit has the exact solution r(t) = R (cos wt, sin wt, 0)
/// with w = sqrt(GM / R^3). Integrating to a fixed time T at step dt = T/n and
/// comparing to that analytic position gives the global error E(n). Halving dt
/// (doubling n) scales the error by 2^p, so the empirical order is
/// p = log2(E(n) / E(2n)). We assert each integrator hits its claimed order:
/// euler ~ 1, runge_kutta_4 ~ 4, rkf78 high (>= 5; the measurable order before
/// the round-off floor).
///
/// Exit 0 iff every measured order matches the claim.

#include "constants/constants_provider.h"
#include "dynamics/inertia.h"
#include "dynamics/state.h"
#include "forces/gravity_central.h"
#include "integrators/runge_kutta.h"
#include "integrators/runge_kutta_fehlberg.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

using T = double;
using TV = math::TrackedValue<T>;
using math::Quaternion;
using math::Vector3;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool ok, double measured) {
    if (ok) {
        ++passed;
        std::cout << "  PASS: " << name << "  (measured order " << measured << ")\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "  (measured order " << measured << ")\n";
    }
}

TV tv(double v) { return TV(v, T(0), TV::representation_bound(v), T(0)); }

} // anonymous namespace

int main() {
    std::cout << std::setprecision(6);
    std::cout << "test_integrator_order: empirical convergence order (REQ-IN-4)\n\n";

    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    const T GM = K.earth.GM.value;
    const T R = K.earth.a.value + 400000.0;
    const T v_orbit = std::sqrt(GM / R);
    const T w = v_orbit / R;                       // angular rate, rad/s
    const T pi = boost::math::constants::pi<T>();
    const T period = T(2) * pi / w;
    const T Tend = period / 8.0;                   // 1/8 orbit (theta = pi/4)

    // analytic position at Tend
    const T theta = w * Tend;
    const T ax = R * std::cos(theta);
    const T ay = R * std::sin(theta);

    Vector3<T> r0(tv(R), tv(0.0), tv(0.0));
    Vector3<T> v0(tv(0.0), tv(v_orbit), tv(0.0));
    Vector3<T> omega0;
    dynamics::State<T> state0 = dynamics::State<T>::from_kinematics(
        Quaternion<T>::identity(), r0, omega0, v0, tv(0.0));
    dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(tv(1.0));

    integrators::AccelFn<T> accel = [&](const dynamics::State<T>& s) {
        return inertia.acceleration_from_wrench(forces::gravity_central(s, K));
    };

    using Integ = std::function<dynamics::State<T>(
        const dynamics::State<T>&, const TV&, const integrators::AccelFn<T>&)>;

    auto err_at = [&](const Integ& integ, int n) {
        const T dt = Tend / n;
        dynamics::State<T> y = state0;
        for (int i = 0; i < n; ++i) y = integ(y, tv(dt), accel);
        Vector3<T> r = y.position();
        return std::sqrt(std::pow(r.x.value - ax, 2) + std::pow(r.y.value - ay, 2)
                                                      + std::pow(r.z.value, 2));
    };

    auto order = [&](const char* name, const Integ& integ, int n) {
        const T e1 = err_at(integ, n);
        const T e2 = err_at(integ, 2 * n);
        const double p = std::log2(e1 / e2);
        std::cout << "  " << name << ": E(" << n << ")=" << e1
                  << "  E(" << 2 * n << ")=" << e2 << "  -> order " << p << "\n";
        return p;
    };

    Integ euler = [](const dynamics::State<T>& y, const TV& h,
                     const integrators::AccelFn<T>& f) { return integrators::euler(y, h, f); };
    Integ rk4 = [](const dynamics::State<T>& y, const TV& h,
                   const integrators::AccelFn<T>& f) { return integrators::runge_kutta_4(y, h, f); };
    Integ rkf = [](const dynamics::State<T>& y, const TV& h,
                   const integrators::AccelFn<T>& f) { return integrators::rkf78(y, h, f); };

    std::cout << "=== measured convergence orders ===\n";
    const double p_euler = order("euler ", euler, 400);
    const double p_rk4   = order("rk4   ", rk4, 24);
    const double p_rkf   = order("rkf78 ", rkf, 5);

    std::cout << "\n";
    check("euler is first order (0.8 <= p <= 1.5)", p_euler >= 0.8 && p_euler <= 1.5, p_euler);
    check("runge_kutta_4 is fourth order (3.5 <= p <= 4.5)", p_rk4 >= 3.5 && p_rk4 <= 4.5, p_rk4);
    check("rkf78 is high order (p >= 5)", p_rkf >= 5.0, p_rkf);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
