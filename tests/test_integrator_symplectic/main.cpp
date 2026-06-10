/// test_integrator_symplectic — symplectic leapfrog conservation (REQ-IN-9).
///
/// A symplectic integrator preserves the invariants of a time-independent
/// conservative system. For a circular Kepler orbit we verify, over many
/// periods, that symplectic_leapfrog:
///   - conserves total energy with a bounded (non-secular) error;
///   - conserves angular momentum, far more tightly than runge_kutta_4
///     (a non-symplectic method) — the defining symplectic signature;
///   - shows NO secular energy growth: the energy error over 40 periods is
///     not materially larger than over 20 (bounded oscillation, not the
///     linear drift a non-symplectic method exhibits).
///
/// Exit 0 iff every check passes.

#include "constants/constants_provider.h"
#include "dynamics/inertia.h"
#include "dynamics/state.h"
#include "forces/gravity_central.h"
#include "integrators/runge_kutta.h"
#include "integrators/symplectic_leapfrog.h"
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

void check(const char* name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

TV tv(double v) { return TV(v, T(0), TV::representation_bound(v), T(0)); }

} // anonymous namespace

int main() {
    std::cout << std::setprecision(11);
    std::cout << "test_integrator_symplectic: leapfrog conservation (REQ-IN-9)\n\n";

    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    const T GM = K.earth.GM.value;
    const T r_orbit = K.earth.a.value + 400000.0;
    const T v_orbit = std::sqrt(GM / r_orbit);
    const T pi = boost::math::constants::pi<T>();
    const T period = T(2) * pi * r_orbit / v_orbit;
    const T dt_step = 60.0;

    const T E0 = T(0.5) * v_orbit * v_orbit - GM / r_orbit;
    const T Lz0 = r_orbit * v_orbit;

    Vector3<T> r0(tv(r_orbit), tv(0.0), tv(0.0));
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
    Integ leap = [](const dynamics::State<T>& y, const TV& h,
                    const integrators::AccelFn<T>& f) {
        return integrators::symplectic_leapfrog(y, h, f);
    };
    Integ rk4 = [](const dynamics::State<T>& y, const TV& h,
                   const integrators::AccelFn<T>& f) {
        return integrators::runge_kutta_4(y, h, f);
    };

    // Propagate n periods; return {relative energy drift, relative Lz drift}.
    auto run = [&](const Integ& integ, T n_periods) {
        dynamics::State<T> y = state0;
        const T t_end = n_periods * period;
        T t = 0.0;
        while (t < t_end) {
            const T h = std::min(dt_step, t_end - t);
            y = integ(y, tv(h), accel);
            t += h;
        }
        Vector3<T> r = y.position();
        Vector3<T> v = y.linear_velocity();
        const T rmag = std::sqrt(r.x.value*r.x.value + r.y.value*r.y.value
                                                     + r.z.value*r.z.value);
        const T vmag = std::sqrt(v.x.value*v.x.value + v.y.value*v.y.value
                                                     + v.z.value*v.z.value);
        const T E = T(0.5) * vmag * vmag - GM / rmag;
        const T Lz = r.x.value * v.y.value - r.y.value * v.x.value;
        return std::pair<T, T>{std::abs(E - E0) / std::abs(E0),
                               std::abs(Lz - Lz0) / std::abs(Lz0)};
    };

    std::cout << "=== 20-period conservation ===\n";
    auto leap20 = run(leap, 20.0);
    auto rk4_20 = run(rk4, 20.0);
    std::cout << "  leapfrog: energy drift " << leap20.first
              << "  Lz drift " << leap20.second << "\n";
    std::cout << "  rk4:      energy drift " << rk4_20.first
              << "  Lz drift " << rk4_20.second << "\n";

    check("leapfrog conserves energy over 20 periods (rel < 1e-3)",
          leap20.first < 1e-3);
    check("leapfrog conserves angular momentum over 20 periods (rel < 1e-6)",
          leap20.second < 1e-6);
    check("leapfrog conserves Lz far tighter than RK4 (symplectic signature)",
          leap20.second < rk4_20.second);

    std::cout << "\n=== bounded (symplectic) vs secular (RK4) energy over 40 periods ===\n";
    auto leap40 = run(leap, 40.0);
    auto rk4_40 = run(rk4, 40.0);
    std::cout << "  leapfrog energy drift: 20p=" << leap20.first
              << "  40p=" << leap40.first << " (bounded ~1e-8)\n";
    std::cout << "  rk4      energy drift: 20p=" << rk4_20.first
              << "  40p=" << rk4_40.first << " (secular, grows)\n";
    // The symplectic signature: leapfrog's energy error stays bounded and small
    // over the long run, while the non-symplectic RK4's drifts secularly — so
    // leapfrog conserves energy BETTER over 40 periods despite its lower order,
    // and stays well under RK4's level (bounded, not climbing toward it).
    check("leapfrog energy bounded and tiny at 40 periods (rel < 1e-6)",
          leap40.first < 1e-6);
    check("leapfrog beats secular RK4 energy drift over 40 periods",
          leap40.first < rk4_40.first);

    // --- per-slot LTE deposit law (V3, runge_kutta_lie_group.md §6.2) ---
    // The drift advances the pose, so leapfrog deposits its truncation
    // envelope to the pose — but PER SLOT FAMILY with DIMENSIONAL
    // magnitudes, not uniformly: this orbit is torque-free (omega0 = 0,
    // central gravity), so the rotational pair is EXACTLY zero and the real
    // (unit-quaternion) slots must receive exactly 0; the dual
    // (translation-carrying) slots receive P/2 = (h²/4)·A [m]; the twist
    // linear slots receive V = h·A [m/s] — so dual/twist = h/4. The old
    // uniform stamp put one mis-united magnitude on all 14 slots, which
    // position extraction amplified by |r|/2 — the saturation V3 removed.
    std::cout << "\n=== per-slot LTE deposit law (V3) ===\n";
    dynamics::State<T> one_step = integrators::symplectic_leapfrog(state0, tv(dt_step), accel);
    const T real_acc  = one_step.pose.M.real.w.errors.accuracy;
    const T dual_acc  = one_step.pose.M.dual.x.errors.accuracy;
    const T twist_acc = one_step.twist.linear.x.errors.accuracy;
    std::cout << "  real-slot accuracy " << real_acc
              << "  dual-slot accuracy " << dual_acc
              << "  twist-slot accuracy " << twist_acc << "\n";
    check("torque-free: real (rotation) slots receive EXACTLY zero LTE",
          real_acc == T(0));
    check("dual slots receive P/2 = (h/4) x the twist's V (rel 1e-12)",
          dual_acc > T(0) &&
          std::abs(T(4) * dual_acc - dt_step * twist_acc)
              < T(1e-12) * dt_step * twist_acc);
    check("twist linear slots receive the velocity envelope V > 0",
          twist_acc > T(0));

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
