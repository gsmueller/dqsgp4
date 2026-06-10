/// test_adaptive_stepping (gate AD1) — closed-loop RKF7(8) step-size control.
///
///   A. the adaptive driver lands exactly on the target time;
///   B. it reproduces a fixed fine-step reference to a small tolerance;
///   C. a tighter local-error tolerance takes more steps and yields a smaller
///      error vs the reference (the loop genuinely adapts);
///   D. the controller rejects/retries when needed and never stalls.

#include "constants/constants_provider.h"
#include "dynamics/inertia.h"
#include "dynamics/propagator.h"
#include "dynamics/state.h"
#include "forces/gravity_central.h"
#include "integrators/runge_kutta.h"
#include "integrators/runge_kutta_fehlberg.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using T = double;
using TV = math::TrackedValue<T>;

namespace {

int failed = 0;

TV tv(double v) { return TV(static_cast<T>(v), T(0), T(0), T(0)); }

void check(const char* name, bool ok, double detail) {
    if (ok) {
        std::cout << "  PASS: " << name << "  (" << std::setprecision(6) << detail << ")\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "  (" << std::setprecision(17) << detail << ")\n";
    }
}

double dist(const math::Vector3<T>& a, const math::Vector3<T>& b) {
    double dx = a.x.value - b.x.value;
    double dy = a.y.value - b.y.value;
    double dz = a.z.value - b.z.value;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

int main() {
    std::cout << "test_adaptive_stepping (AD1): RKF7(8) step-size control\n\n";

    // A central-gravity propagator supplies the acceleration callback.
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(math::exact<T>(1));
    std::vector<dynamics::ForceFn<T>> forces;
    forces.push_back([](const dynamics::State<T>& s, const constants::ConstantsProvider<T>& KK) {
        return forces::gravity_central(s, KK);
    });
    dynamics::IntegratorFn<T> integ =
        [](const dynamics::State<T>& y0, const TV& dt, const integrators::AccelFn<T>& a) {
            return integrators::rkf78(y0, dt, a);
        };
    dynamics::Propagator<T> prop(K, inertia, std::move(forces), integ);
    integrators::AccelFn<T> accel =
        [&prop](const dynamics::State<T>& s) { return prop.compute_acceleration(s); };

    // Circular LEO at 7000 km: v = sqrt(mu/r) ~ 7546 m/s.
    math::Vector3<T> r(tv(7.0e6), tv(0.0), tv(0.0));
    math::Vector3<T> v(tv(0.0), tv(7546.0), tv(0.0));
    math::Vector3<T> omega_body;
    dynamics::State<T> y0 = dynamics::State<T>::from_kinematics(
        math::Quaternion<T>::identity(), r, omega_body, v, tv(0.0));

    const TV t_target = math::exact<T>(600);   // 600 s
    const T dt_min = T(1) / T(100);            // 0.01 s floor

    // Reference: fixed fine RKF7(8) — 600 one-second steps.
    dynamics::State<T> ref = y0;
    for (int i = 0; i < 600; ++i) ref = integrators::rkf78(ref, math::exact<T>(1), accel);

    // Adaptive runs at three tolerances.
    integrators::AdaptiveResult<T> loose =
        integrators::rkf78_propagate_adaptive(y0, t_target, math::exact<T>(10), accel, T(1e-2), dt_min);
    integrators::AdaptiveResult<T> med =
        integrators::rkf78_propagate_adaptive(y0, t_target, math::exact<T>(10), accel, T(1e-6), dt_min);
    integrators::AdaptiveResult<T> tight =
        integrators::rkf78_propagate_adaptive(y0, t_target, math::exact<T>(10), accel, T(1e-10), dt_min);

    double e_loose = dist(loose.state.position(), ref.position());
    double e_med = dist(med.state.position(), ref.position());
    double e_tight = dist(tight.state.position(), ref.position());
    std::cout << "  [diag] loose(1e-2): " << loose.accepted_steps << " steps, "
              << loose.rejected_steps << " rej, err " << e_loose << " m\n";
    std::cout << "  [diag] med  (1e-6): " << med.accepted_steps << " steps, "
              << med.rejected_steps << " rej, err " << e_med << " m\n";
    std::cout << "  [diag] tight(1e-10): " << tight.accepted_steps << " steps, "
              << tight.rejected_steps << " rej, err " << e_tight << " m\n\n";

    // A. Lands exactly on the target time.
    check("lands on t_target", std::abs(med.state.time.value - 600.0) < 1e-9,
          med.state.time.value);

    // B. Reproduces the fine reference.
    check("med (1e-6) matches fine reference (< 1 m)", e_med < 1.0, e_med);

    // C. Tighter tolerance -> more steps and smaller error.
    check("tighter tol -> more steps", tight.accepted_steps > loose.accepted_steps,
          static_cast<double>(tight.accepted_steps - loose.accepted_steps));
    check("tighter tol -> smaller error", e_tight <= e_loose, e_tight - e_loose);

    // D. The loop ran (steps taken) and stayed bounded.
    check("accepted steps > 0", med.accepted_steps > 0,
          static_cast<double>(med.accepted_steps));
    check("step count bounded (< 100000)", med.accepted_steps < 100000,
          static_cast<double>(med.accepted_steps));

    // E. R3b — the facade method is the standalone loop with the same callback,
    // so the result must be BIT-IDENTICAL (dq_propagator_facade.md section 4).
    integrators::AdaptiveResult<T> fac =
        prop.propagate_adaptive(y0, t_target, math::exact<T>(10), T(1e-6), dt_min);
    check("facade propagate_adaptive position BIT-IDENTICAL to standalone",
          fac.state.position().x.value == med.state.position().x.value &&
          fac.state.position().y.value == med.state.position().y.value &&
          fac.state.position().z.value == med.state.position().z.value,
          dist(fac.state.position(), med.state.position()));
    check("facade step diagnostics identical",
          fac.accepted_steps == med.accepted_steps &&
          fac.rejected_steps == med.rejected_steps,
          static_cast<double>(fac.accepted_steps));

    std::cout << "\n" << (failed == 0 ? "PASS" : "FAIL") << " — " << failed
              << " failure(s)\n";
    return failed == 0 ? 0 : 1;
}
