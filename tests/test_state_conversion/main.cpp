/// test_state_conversion — F1 (DQSGP4 Completion Roadmap, issue F1).
///
/// Verifies the State <-> StateVector converters that give the SGP4 and DQSGP4
/// APIs a shared state vocabulary:
///   - unit scaling km <-> m;
///   - round-trip StateVector -> State -> StateVector (identity attitude);
///   - the DQ body-frame velocity is rotated into the world frame (magnitude
///     preserved, direction changed under a non-identity attitude);
///   - the new per-category velocity error accessors (F1-c).
///
/// ExeGate F1: nonzero exit on any failed check.

#include "dynamics/state_conversion.h"
#include "math/quaternion.h"
#include "math/vector3.h"

#include <cmath>
#include <iostream>

namespace {

using T = double;
using TV = math::TrackedValue<T>;
using V = math::Vector3<T>;
const double PI = 3.14159265358979323846;

int passed = 0;
int failed = 0;

void check(const char* name, bool ok) {
    if (ok) {
        ++passed;
        std::cout << "  PASS: " << name << "\n";
    } else {
        ++failed;
        std::cout << "  FAIL: " << name << "\n";
    }
}

bool close_d(double a, double b, double tol) {
    return std::abs(a - b) <= tol;
}

// A clean input TrackedValue: only the binary-representation precision.
TV mk(double v) {
    return TV(T(v), T(0), TV::representation_bound(T(v)), T(0));
}

} // namespace

int main() {
    // A StateVector in km / TEME: 7000 km radial, 7.5 km/s cross-track.
    sgp4::StateVector<T> sv;
    sv.position_km = V(mk(7000.0), mk(0.0), mk(0.0));
    sv.velocity_km_s = V(mk(0.0), mk(7.5), mk(0.0));

    // km -> m, identity attitude.
    dynamics::State<T> st = dynamics::to_state(sv);
    check("to_state position x == 7e6 m", close_d(st.position().x.value, 7.0e6, 1e-3));
    check("to_state velocity y == 7500 m/s", close_d(st.linear_velocity().y.value, 7500.0, 1e-6));
    check("to_state identity attitude (w==1)", close_d(st.orientation().w.value, 1.0, 1e-12));

    // Round-trip back to km / TEME.
    sgp4::StateVector<T> sv2 = dynamics::to_state_vector(st);
    check("round-trip position x == 7000 km", close_d(sv2.position_km.x.value, 7000.0, 1e-9));
    check("round-trip velocity y == 7.5 km/s", close_d(sv2.velocity_km_s.y.value, 7.5, 1e-12));

    // F1-c: per-category velocity accessors exist and are non-negative.
    check("velocity_measurement_error >= 0", sv2.velocity_measurement_error() >= 0.0);
    check("velocity_precision_error >= 0", sv2.velocity_precision_error() >= 0.0);
    check("velocity_accuracy_error >= 0", sv2.velocity_accuracy_error() >= 0.0);

    // Non-identity attitude: the body-frame velocity must be rotated into the
    // world frame. A 90-degree rotation about +z moves a body +y velocity onto
    // the world x-axis; magnitude is preserved.
    math::Quaternion<T> q = math::Quaternion<T>::from_axis_angle(
        V(mk(0.0), mk(0.0), mk(1.0)), mk(PI / 2.0));
    V body_vel(mk(0.0), mk(7500.0), mk(0.0));  // m/s, body frame
    dynamics::State<T> rot = dynamics::State<T>::from_kinematics(
        q, V(mk(7.0e6), mk(0.0), mk(0.0)), V(), body_vel, mk(0.0));
    sgp4::StateVector<T> rsv = dynamics::to_state_vector(rot);
    check("rotated velocity magnitude preserved (7.5 km/s)",
          close_d(rsv.velocity_km_s.magnitude().value, 7.5, 1e-9));
    check("rotated velocity now along x (|x| ~ 7.5)",
          close_d(std::abs(rsv.velocity_km_s.x.value), 7.5, 1e-9));
    check("rotated velocity y ~ 0", close_d(rsv.velocity_km_s.y.value, 0.0, 1e-9));

    std::cout << "\n  state conversion: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
