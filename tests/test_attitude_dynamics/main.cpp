/// test_attitude_dynamics — H1 (DQSGP4 Completion Roadmap, register H1 / DEFERRED).
///
/// Verifies the coupled-axis rigid-body attitude dynamics built in H1:
///   - math::Matrix3<T> inverse is correct for a general (off-diagonal) inertia
///     tensor with products of inertia: I·(I⁻¹·v) = v;
///   - the full Euler equation ω̇ = I⁻¹(τ − ω×Iω) reproduces the ANALYTIC
///     free symmetric-top precession: the transverse angular velocity rotates
///     at Ω = ω_z(C−A)/A in the body frame;
///   - torque-free motion conserves rotational energy ½ω·Iω and the angular
///     momentum magnitude |I·ω|;
///   - the integrated attitude quaternion stays unit-norm;
///   - the point-mass torque restriction is lifted: a nonzero torque on a body
///     with a real inertia tensor gives a finite ω̇ (no zero-moment catastrophe)
///     and the gyroscopic coupling cross-feeds the axes;
///   - precision tightens with a wider numeric type T.
///
/// ExeGate H1: nonzero exit code on any failed check.

#include "dynamics/attitude.h"
#include "math/matrix3.h"
#include "math/quaternion.h"
#include "math/vector3.h"
#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
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
math::Vector3<T> vec(double x, double y, double z) {
    return math::Vector3<T>(tv<T>(x), tv<T>(y), tv<T>(z));
}

void test_matrix3_inverse() {
    using T = double;
    // General symmetric inertia tensor with products of inertia (Matrix3 path).
    math::Matrix3<T> I = math::Matrix3<T>::symmetric(
        tv<T>(2.0), tv<T>(3.0), tv<T>(4.0), tv<T>(0.5), tv<T>(0.3), tv<T>(0.2));
    math::Matrix3<T> Iinv = I.inverse();
    math::Vector3<T> v = vec<T>(1.0, -2.0, 3.0);
    math::Vector3<T> roundtrip = I * (Iinv * v);  // I·(I⁻¹·v) should be v
    check("Matrix3 inverse: I·(I⁻¹·v) = v (off-diagonal tensor)",
          std::abs(roundtrip.x.value - 1.0) < 1e-12
          && std::abs(roundtrip.y.value + 2.0) < 1e-12
          && std::abs(roundtrip.z.value - 3.0) < 1e-12);
}

void test_free_symmetric_top() {
    using T = double;
    // Symmetric top: I = diag(A, A, C), A=2, C=1; torque-free; ω0 = (0.1, 0, 1).
    const double A = 2.0, C = 1.0, wz = 1.0, wt = 0.1;
    math::Matrix3<T> I = math::Matrix3<T>::diagonal(vec<T>(A, A, C));
    math::Matrix3<T> Iinv = I.inverse();
    math::Vector3<T> torque;  // zero

    dynamics::AttitudeState<T> s{math::Quaternion<T>::identity(), vec<T>(wt, 0.0, wz)};
    double E0 = dynamics::rotational_energy(I, s.omega).value;
    double L0 = dynamics::angular_momentum(I, s.omega).magnitude().value;

    const int N = 1000;
    math::TrackedValue<T> dt = tv<T>(0.001);  // total t = 1 s
    for (int i = 0; i < N; ++i) s = dynamics::attitude_rk4_step(s, I, Iinv, torque, dt);

    // Analytic free-top body precession rate Ω = ω_z(C−A)/A; after t = 1 s the
    // transverse (ω_x, ω_y) has rotated by Ω·t. Here Ω = 1·(1−2)/2 = −0.5 rad/s.
    double Omega = wz * (C - A) / A;
    double angle = std::atan2(s.omega.y.value, s.omega.x.value);
    check("free symmetric-top precession matches analytic Ω·t",
          std::abs(angle - Omega * 1.0) < 1e-4);
    check("spin-axis ω_z conserved", std::abs(s.omega.z.value - wz) < 1e-9);
    double wxy = std::sqrt(s.omega.x.value * s.omega.x.value + s.omega.y.value * s.omega.y.value);
    check("transverse |ω_xy| conserved", std::abs(wxy - wt) < 1e-6);

    double E1 = dynamics::rotational_energy(I, s.omega).value;
    double L1 = dynamics::angular_momentum(I, s.omega).magnitude().value;
    check("rotational energy conserved (torque-free)", std::abs(E1 - E0) < 1e-9);
    check("angular-momentum magnitude conserved", std::abs(L1 - L0) < 1e-9);

    double qn = std::sqrt(s.q.w.value * s.q.w.value + s.q.x.value * s.q.x.value
                        + s.q.y.value * s.q.y.value + s.q.z.value * s.q.z.value);
    check("attitude quaternion stays unit-norm", std::abs(qn - 1.0) < 1e-9);
}

void test_torque_and_coupling() {
    using T = double;
    // Full asymmetric inertia; a nonzero torque must give a FINITE ω̇ (the
    // point-mass zero-moment catastrophe does not apply to a real tensor).
    math::Matrix3<T> I = math::Matrix3<T>::diagonal(vec<T>(2.0, 3.0, 1.0));
    math::Matrix3<T> Iinv = I.inverse();
    // ω spans the x and z axes (different moments Ix=2, Iz=1), so the
    // gyroscopic term ω×(Iω) has a nonzero y-component (Ix−Iz)ω_zω_x.
    math::Vector3<T> omega = vec<T>(1.0, 0.0, 5.0);
    math::Vector3<T> torque = vec<T>(1.0, 0.0, 0.0);  // torque about x
    math::Vector3<T> wdot = dynamics::euler_angular_acceleration(I, Iinv, omega, torque);
    bool finite = std::isfinite(wdot.x.value) && std::isfinite(wdot.y.value)
               && std::isfinite(wdot.z.value);
    check("torque on a real inertia gives finite ω̇ (point-mass restriction lifted)", finite);
    // Gyroscopic coupling: the −ω×(Iω) term cross-couples the axes, so ω̇_y is
    // driven by the x/z spin even though no torque acts about y. A decoupled
    // τ/I model (the old diagonal Inertia path) would give ω̇_y = 0 exactly.
    // Here gyro_y = (Ix−Iz)ω_zω_x = (2−1)·5·1 = 5 ⇒ ω̇_y = −5/Iy = −5/3.
    check("gyroscopic coupling cross-feeds the axes (ω̇_y ≠ 0)",
          std::abs(wdot.y.value - (-5.0 / 3.0)) < 1e-9);
}

void test_precision_scales() {
    auto run = [](auto tag) -> double {
        using T = decltype(tag);
        math::Matrix3<T> I = math::Matrix3<T>::diagonal(vec<T>(2.0, 2.0, 1.0));
        math::Matrix3<T> Iinv = I.inverse();
        math::Vector3<T> wdot = dynamics::euler_angular_acceleration(
            I, Iinv, vec<T>(0.1, 0.0, 1.0), math::Vector3<T>());
        return static_cast<double>(wdot.y.errors.precision);
    };
    double pd = run(double{});
    double pb = run(cpp_bin_float_50{});
    check("angular-acceleration precision > 0", pd > 0.0 && pb > 0.0);
    check("angular-acceleration precision tightens with wider T", pb < pd * 1e-20);
}

} // namespace

int main() {
    test_matrix3_inverse();
    test_free_symmetric_top();
    test_torque_and_coupling();
    test_precision_scales();

    std::cout << "\n  attitude dynamics: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
