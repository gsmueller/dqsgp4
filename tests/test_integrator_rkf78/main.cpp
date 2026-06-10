/// test_integrator_rkf78 — adaptive RKF7(8) integrator (REQ-IN-4/6).
///
/// Three layers of verification:
///   1. Tableau correctness, checked EXACTLY from the coefficients (no
///      numerical reference): every stage row sums to its node
///      (sum_j a[i][j] = c[i]); the 8th-order weights sum to 1; and the
///      quadrature order conditions sum_i b[i] c[i]^k = 1/(k+1) hold for
///      k = 0..7 (necessary for order 8). A transcription error in any
///      coefficient breaks one of these.
///   2. Integration accuracy: propagate a circular Kepler orbit one period
///      and confirm rkf78 closes far tighter than runge_kutta_4 at the same
///      step (the higher-order method is dramatically more accurate).
///   3. Adaptive feedback (REQ-IN-6): rkf78_step returns a finite, positive,
///      small embedded local-error estimate.
///
/// Exit 0 iff every check passes.

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
    std::cout << "test_integrator_rkf78: RKF7(8) tableau + orbit + adaptive error\n\n";

    // ---- 1. Tableau order conditions (exact, no reference) ----
    std::cout << "=== tableau correctness ===\n";
    auto tab = integrators::rkf78_tableau<T>();
    const int N = integrators::kRkf78Stages;

    bool rows_ok = true;
    for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) s += tab.a[i][j].value;
        if (std::abs(s - tab.c[i].value) > 1e-10) rows_ok = false;
    }
    check("each stage row sums to its node (sum_j a[i][j] = c[i])", rows_ok);

    double bsum = 0.0;
    for (int i = 0; i < N; ++i) bsum += tab.b[i].value;
    check("8th-order weights sum to 1", std::abs(bsum - 1.0) < 1e-12);

    bool quad_ok = true;
    for (int k = 0; k <= 7; ++k) {
        double s = 0.0;
        for (int i = 0; i < N; ++i) s += tab.b[i].value * std::pow(tab.c[i].value, k);
        if (std::abs(s - 1.0 / (k + 1)) > 1e-10) quad_ok = false;
    }
    check("quadrature order conditions sum_i b_i c_i^k = 1/(k+1), k=0..7", quad_ok);

    // ---- shared circular-orbit setup (mirrors test_propagator Phase 1) ----
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    const T GM = K.earth.GM.value;
    const T r_orbit = K.earth.a.value + 400000.0;
    const T v_orbit = std::sqrt(GM / r_orbit);
    const T pi = boost::math::constants::pi<T>();
    const T period = T(2) * pi * r_orbit / v_orbit;
    const T dt_step = 60.0;

    Vector3<T> r0(tv(r_orbit), tv(0.0), tv(0.0));
    Vector3<T> v0(tv(0.0), tv(v_orbit), tv(0.0));
    Vector3<T> omega0;
    dynamics::State<T> state0 = dynamics::State<T>::from_kinematics(
        Quaternion<T>::identity(), r0, omega0, v0, tv(0.0));
    dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(tv(1.0));

    integrators::AccelFn<T> accel = [&](const dynamics::State<T>& s) {
        return inertia.acceleration_from_wrench(forces::gravity_central(s, K));
    };

    auto closure = [&](bool use_rkf78) {
        dynamics::State<T> y = state0;
        T t = 0.0;
        while (t < period) {
            const T h = std::min(dt_step, period - t);   // exact final partial step
            y = use_rkf78 ? integrators::rkf78(y, tv(h), accel)
                          : integrators::runge_kutta_4(y, tv(h), accel);
            t += h;
        }
        Vector3<T> rf = y.position();
        return std::sqrt(std::pow(rf.x.value - r0.x.value, 2) +
                         std::pow(rf.y.value - r0.y.value, 2) +
                         std::pow(rf.z.value - r0.z.value, 2));
    };

    // ---- 2. Integration accuracy ----
    std::cout << "\n=== orbit integration (one period, dt = 60 s) ===\n";
    const T rk4_err = closure(false);
    const T rkf_err = closure(true);
    std::cout << "  RK4    position closure: " << rk4_err << " m\n";
    std::cout << "  RKF78  position closure: " << rkf_err << " m\n";
    check("rkf78 integrates the orbit (closure < 1 m)", rkf_err < 1.0);
    check("rkf78 is more accurate than RK4 at the same step", rkf_err < rk4_err);

    // ---- 3. Adaptive embedded error estimate (REQ-IN-6) ----
    std::cout << "\n=== adaptive error feedback ===\n";
    integrators::StepResult<T> step = integrators::rkf78_step(state0, tv(dt_step), accel);
    const T e = step.error.value;
    std::cout << "  embedded local error estimate: " << e << "\n";
    check("error estimate is finite", std::isfinite(e));
    check("error estimate is positive", e > 0.0);
    check("error estimate is small (< 1e-3 of orbit scale)", e < 1e-3 * r_orbit);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
