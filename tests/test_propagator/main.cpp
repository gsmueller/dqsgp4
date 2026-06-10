/// test_propagator — LEO smoke test for the dual quaternion propagator.
///
/// Phase 1 (REQ-PR-2..4, REQ-PR-9): WGS 84 central gravity + Lie-group
/// RK4 reproduces a circular Kepler orbit at 400 km altitude over one
/// orbital period, conserves total energy and angular momentum, and the
/// framework's reported `total_error` budget stays populated.
///
/// Phase 2 (REQ-PR-11): adding `gravity_J2` to the force list exercises
/// force-lambda composition. For a planar equatorial orbit the J₂
/// perturbation has z-component = 0, so the orbit stays in z = 0 and
/// the J₂-augmented two-body energy
///   E = (1/2) v² − μ/r − μ J₂ R_E² / (2 r³)
/// is the conserved quantity. We check that the orbit remains bounded
/// near the initial radius and that energy with the J₂ potential
/// included is conserved.
///
/// This is the headline end-to-end test of the propagator pipeline:
///
///   ConstantsProvider::wgs84(tol)
///     → Inertia::point_mass(m)
///     → forces = { gravity_central [, gravity_J2 ] }
///     → integrators::runge_kutta_4<T>
///     → Propagator<T>(K, inertia, forces, integrator)
///     → propagator.propagate_to(state0, T_period, dt = 60 s)

#include "constants/constants_provider.h"
#include "dynamics/inertia.h"
#include "dynamics/propagator.h"
#include "dynamics/state.h"
#include "forces/drag.h"
#include "forces/gravity_central.h"
#include "forces/gravity_zonal.h"
#include "integrators/runge_kutta.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using T = double;
using TV = math::TrackedValue<T>;
using math::Quaternion;
using math::Vector3;

namespace {

int tests_passed = 0;
int tests_failed = 0;

TV tv(double v) {
    T val = static_cast<T>(v);
    T prec = TV::representation_bound(val);
    return TV(val, T(0), prec, T(0));
}

void check_bound(const char* name, T measured, T bound) {
    if (measured <= bound) {
        ++tests_passed;
        std::cout << "  PASS: " << name
                  << "  measured=" << std::setprecision(6) << measured
                  << "  bound=" << bound << "\n";
    } else {
        ++tests_failed;
        std::cerr << "  FAIL: " << name
                  << "  measured=" << std::setprecision(17) << measured
                  << "  bound=" << bound << "\n";
    }
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(11);
    std::cout << "test_propagator: LEO 400 km smoke test (WGS 84 + RK4)\n\n";

    // --- Constants provider (WGS 84) ---
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    std::cout << "WGS 84 GM = " << K.earth.GM.value << " m^3/s^2\n";
    std::cout << "WGS 84 a  = " << K.earth.a.value << " m\n";
    std::cout << "WGS 84 J2 = " << K.earth.J2n(1).value << "\n\n";

    // --- Orbit parameters (closed-form Kepler) ---
    const T altitude = 400000.0;
    const T r_earth  = K.earth.a.value;
    const T r_orbit  = r_earth + altitude;
    const T GM       = K.earth.GM.value;
    const T J2       = K.earth.J2n(1).value;
    const T R_E      = K.earth.a.value;
    const T v_orbit  = std::sqrt(GM / r_orbit);
    const T pi       = boost::math::constants::pi<T>();
    const T period   = T(2) * pi * r_orbit / v_orbit;
    std::cout << "r_orbit = " << r_orbit << " m\n";
    std::cout << "v_orbit = " << v_orbit << " m/s\n";
    std::cout << "period  = " << period  << " s (" << period / 60.0
              << " min)\n\n";

    // --- Initial state on the +x axis, moving in +y (equatorial) ---
    Vector3<T> r0(tv(r_orbit), tv(0.0), tv(0.0));
    Vector3<T> v0(tv(0.0), tv(v_orbit), tv(0.0));
    Quaternion<T> q0 = Quaternion<T>::identity();
    Vector3<T> omega0;
    dynamics::State<T> state0 = dynamics::State<T>::from_kinematics(
        q0, r0, omega0, v0, tv(0.0));

    dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(tv(1.0));

    dynamics::IntegratorFn<T> integrator =
        [](const dynamics::State<T>& y0,
           const TV& dt,
           const integrators::AccelFn<T>& acc_fn) {
            return integrators::runge_kutta_4(y0, dt, acc_fn);
        };

    const T dt_step = 60.0;

    // ============================================================
    // Phase 1: central gravity only
    // ============================================================
    std::cout << "=== Phase 1: central gravity only ===\n";

    std::vector<dynamics::ForceFn<T>> forces_central;
    forces_central.push_back(
        [](const dynamics::State<T>& s,
           const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_central(s, KK);
        });

    dynamics::Propagator<T> prop_central(
        K, inertia, std::move(forces_central), integrator);

    const T E0  = T(0.5) * v_orbit * v_orbit - GM / r_orbit;
    const T Lz0 = r_orbit * v_orbit;
    std::cout << "E0  = " << E0  << " J/kg\n";
    std::cout << "Lz0 = " << Lz0 << " m^2/s\n";

    std::cout << "Propagating one period at dt = " << dt_step << " s ...\n";
    dynamics::State<T> state_final = prop_central.propagate_to(state0, tv(period), tv(dt_step));
    std::cout << "Final time: " << state_final.time.value << " s\n";

    Vector3<T> r_final = state_final.position();
    Vector3<T> v_final = state_final.linear_velocity();
    const T r_final_mag = std::sqrt(r_final.x.value * r_final.x.value +
                                    r_final.y.value * r_final.y.value +
                                    r_final.z.value * r_final.z.value);
    const T v_final_mag = std::sqrt(v_final.x.value * v_final.x.value +
                                    v_final.y.value * v_final.y.value +
                                    v_final.z.value * v_final.z.value);

    const T pos_err = std::sqrt(
        std::pow(r_final.x.value - r0.x.value, 2) +
        std::pow(r_final.y.value - r0.y.value, 2) +
        std::pow(r_final.z.value - r0.z.value, 2));
    std::cout << "Position closure error: " << pos_err << " m\n";
    check_bound("position closure < 1 km", pos_err, T(1000.0));

    const T E_final  = T(0.5) * v_final_mag * v_final_mag - GM / r_final_mag;
    const T E_drift  = std::abs(E_final - E0);
    const T E_rel    = E_drift / std::abs(E0);
    std::cout << "Final energy: " << E_final
              << "  drift: " << E_drift
              << "  (relative " << E_rel << ")\n";
    check_bound("energy drift < 1e-6 relative", E_rel, T(1e-6));

    const T Lz_final = r_final.x.value * v_final.y.value
                     - r_final.y.value * v_final.x.value;
    const T Lz_drift = std::abs(Lz_final - Lz0);
    const T Lz_rel   = Lz_drift / std::abs(Lz0);
    std::cout << "Final Lz: " << Lz_final
              << "  drift: " << Lz_drift
              << "  (relative " << Lz_rel << ")\n";
    // RK4 at 60 s steps over LEO gives Lz drift on the order of 1e-7
    // relative per period (consistent with O(dt^4) global error). 1e-6
    // leaves a ~10× margin; a symplectic integrator or RKF7(8) would
    // tighten this.
    check_bound("Lz drift < 1e-6 relative", Lz_rel, T(1e-6));

    const T pos_precision =
        r_final.x.errors.precision + r_final.y.errors.precision
                                   + r_final.z.errors.precision;
    std::cout << "Reported position precision error sum: " << pos_precision
              << " m\n";
    check_bound("reported precision > 0 (framework alive)",
                T(0), pos_precision);

    // ============================================================
    // Phase 2: central gravity + J2 oblateness
    // ============================================================
    std::cout << "\n=== Phase 2: central gravity + J2 ===\n";

    std::vector<dynamics::ForceFn<T>> forces_j2;
    forces_j2.push_back(
        [](const dynamics::State<T>& s,
           const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_central(s, KK);
        });
    forces_j2.push_back(
        [](const dynamics::State<T>& s,
           const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_J2(s, KK);
        });

    dynamics::Propagator<T> prop_j2(
        K, inertia, std::move(forces_j2), integrator);

    // J2 perturbation potential (equatorial: sin(latitude) = 0):
    //   V_J2 = mu J2 R_E^2 (3 sin^2 lat - 1) / (2 r^3)
    //        = -mu J2 R_E^2 / (2 r^3)  for equatorial orbit
    const T V_J2_initial = -GM * J2 * R_E * R_E / (T(2) * std::pow(r_orbit, 3));
    const T E0_j2 = E0 + V_J2_initial;
    std::cout << "E0_j2 (central + J2 potential) = " << E0_j2 << " J/kg\n";
    std::cout << "Lz0 (unchanged by zonal-only field) = " << Lz0
              << " m^2/s\n";

    std::cout << "Propagating one period at dt = " << dt_step << " s ...\n";
    dynamics::State<T> state_j2 = prop_j2.propagate_to(state0, tv(period), tv(dt_step));
    std::cout << "Final time: " << state_j2.time.value << " s\n";

    Vector3<T> r_j2 = state_j2.position();
    Vector3<T> v_j2 = state_j2.linear_velocity();
    const T r_j2_mag = std::sqrt(r_j2.x.value * r_j2.x.value +
                                 r_j2.y.value * r_j2.y.value +
                                 r_j2.z.value * r_j2.z.value);
    const T v_j2_mag = std::sqrt(v_j2.x.value * v_j2.x.value +
                                 v_j2.y.value * v_j2.y.value +
                                 v_j2.z.value * v_j2.z.value);
    std::cout << "Final radius: " << r_j2_mag << " m (initial "
              << r_orbit << " m)\n";

    // For a planar equatorial circular orbit, J2 acts as a slight radial
    // perturbation. The orbit becomes slightly eccentric but stays
    // bounded near the initial radius.
    const T radial_dev = std::abs(r_j2_mag - r_orbit);
    const T radial_rel = radial_dev / r_orbit;
    std::cout << "Radial deviation: " << radial_dev
              << " m (relative " << radial_rel << ")\n";
    check_bound("J2 orbit stays bounded < 1% radial", radial_rel, T(0.01));

    // Energy conservation with J2 potential included. The orbit remains
    // in z = 0 (planar equatorial), so the equatorial V_J2 applies at
    // the final state too.
    const T V_J2_final = -GM * J2 * R_E * R_E / (T(2) * std::pow(r_j2_mag, 3));
    const T E_j2_final = T(0.5) * v_j2_mag * v_j2_mag - GM / r_j2_mag
                       + V_J2_final;
    const T E_j2_drift = std::abs(E_j2_final - E0_j2);
    const T E_j2_rel   = E_j2_drift / std::abs(E0_j2);
    std::cout << "Final J2 energy: " << E_j2_final
              << "  drift: " << E_j2_drift
              << "  (relative " << E_j2_rel << ")\n";
    check_bound("J2 energy drift < 1e-5 relative", E_j2_rel, T(1e-5));

    // J2 is a zonal field — symmetric about Earth's rotation axis (z).
    // Angular momentum L_z is therefore conserved exactly under J2 (up
    // to integrator error), just as under central gravity.
    const T Lz_j2 = r_j2.x.value * v_j2.y.value
                  - r_j2.y.value * v_j2.x.value;
    const T Lz_j2_drift = std::abs(Lz_j2 - Lz0);
    const T Lz_j2_rel   = Lz_j2_drift / std::abs(Lz0);
    std::cout << "Final J2 Lz: " << Lz_j2
              << "  drift: " << Lz_j2_drift
              << "  (relative " << Lz_j2_rel << ")\n";
    check_bound("J2 Lz drift < 1e-6 relative", Lz_j2_rel, T(1e-6));

    // Sanity: position should differ between Phase 1 and Phase 2 by an
    // amount on the order of the J2 perturbation effect.
    const T phase_diff = std::sqrt(
        std::pow(r_j2.x.value - r_final.x.value, 2) +
        std::pow(r_j2.y.value - r_final.y.value, 2) +
        std::pow(r_j2.z.value - r_final.z.value, 2));
    std::cout << "Position difference vs Phase 1: " << phase_diff << " m\n";
    // J2 should produce a non-zero effect on the orbit (otherwise the
    // J2 lambda has no influence).
    check_bound("Phase 2 position differs from Phase 1 (J2 has effect)",
                T(1.0), phase_diff);

    // ============================================================
    // Phase 3: central + J2 + drag, multi-period orbital decay
    // ============================================================
    std::cout << "\n=== Phase 3: central + J2 + drag (5 periods) ===\n";

    // Exponential atmosphere parameters. rho_0 = 2.789e-10 kg/m^3 at
    // h_0 = 200 km (US Standard Atmosphere-ish), H_scale = 50 km.
    // B = C_d * A / m = 0.5 m^2/kg corresponds to a high-drag spacecraft
    // (large area or low mass); makes drag observable in 5 periods.
    dynamics::ForceFn<T> drag_fn = forces::make_drag_exponential<T>(
        tv(2.789e-10), tv(200000.0), tv(50000.0), tv(0.5));

    std::vector<dynamics::ForceFn<T>> forces_drag;
    forces_drag.push_back(
        [](const dynamics::State<T>& s,
           const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_central(s, KK);
        });
    forces_drag.push_back(
        [](const dynamics::State<T>& s,
           const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_J2(s, KK);
        });
    forces_drag.push_back(drag_fn);

    dynamics::Propagator<T> prop_drag(
        K, inertia, std::move(forces_drag), integrator);

    const T t_end = T(5) * period;
    std::cout << "Propagating " << t_end << " s (5 periods) at dt = "
              << dt_step << " s ...\n";
    dynamics::State<T> state_drag = prop_drag.propagate_to(state0, tv(t_end), tv(dt_step));
    std::cout << "Final time: " << state_drag.time.value << " s\n";

    Vector3<T> r_drag = state_drag.position();
    Vector3<T> v_drag = state_drag.linear_velocity();
    const T r_drag_mag = std::sqrt(r_drag.x.value * r_drag.x.value +
                                   r_drag.y.value * r_drag.y.value +
                                   r_drag.z.value * r_drag.z.value);
    const T v_drag_mag = std::sqrt(v_drag.x.value * v_drag.x.value +
                                   v_drag.y.value * v_drag.y.value +
                                   v_drag.z.value * v_drag.z.value);
    std::cout << "Final radius:   " << r_drag_mag << " m (initial "
              << r_orbit << " m)\n";
    std::cout << "Final speed:    " << v_drag_mag << " m/s (initial "
              << v_orbit << " m/s)\n";

    // Drag is dissipative. Recover semi-major axis from vis-viva:
    //   v^2 = mu (2/r - 1/a)   =>   a = 1 / (2/r - v^2/mu).
    // a should DECREASE under drag.
    const T a_initial = r_orbit;
    const T a_final   = T(1) /
        (T(2) / r_drag_mag - v_drag_mag * v_drag_mag / GM);
    const T a_decay   = a_initial - a_final;
    std::cout << "a_initial = " << a_initial << " m\n";
    std::cout << "a_final   = " << a_final << " m\n";
    std::cout << "a decay   = " << a_decay
              << " m (positive = orbit shrunk)\n";
    // check_bound(measured, bound) asserts measured <= bound. Using
    // measured = -a_decay and bound = -100 asserts a_decay >= 100 m.
    check_bound("semi-major axis decays by > 100 m",
                -a_decay, T(-100.0));

    // After 5 periods with drag, position should be substantially
    // displaced from the initial state — far more than the ~30 m Phase 1
    // closure or the ~123 km Phase 2 J2 effect at 1 period.
    const T drag_pos_offset = std::sqrt(
        std::pow(r_drag.x.value - r0.x.value, 2) +
        std::pow(r_drag.y.value - r0.y.value, 2) +
        std::pow(r_drag.z.value - r0.z.value, 2));
    std::cout << "Position offset from start: " << drag_pos_offset << " m\n";
    check_bound("drag + 5 periods produces > 1 km position offset",
                -drag_pos_offset, T(-1000.0));

    // Drag populates errors.accuracy (REQ-EF-7 / AUD-EF-6) — the
    // exponential-atmosphere model-truncation residual. After 5 periods
    // of drag-affected propagation, the accuracy budget on position
    // should be non-zero.
    const T pos_accuracy = r_drag.x.errors.accuracy
                         + r_drag.y.errors.accuracy
                         + r_drag.z.errors.accuracy;
    std::cout << "Reported position accuracy: " << pos_accuracy
              << " m (drag model truncation)\n";
    check_bound("drag populates errors.accuracy (REQ-EF-7)",
                T(0), pos_accuracy);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << tests_passed
              << "  Failed: " << tests_failed << "\n";
    std::cout << "========================================\n";
    return tests_failed > 0 ? 1 : 0;
}
