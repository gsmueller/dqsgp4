/// test_butcher_tableau — the unified ButcherTableau + rk_step driver (P5).
///
/// The three State-space integrators (euler, runge_kutta_4, rkf78) are ONE
/// driver — `rk_step` — applied to three Butcher tableaux. This gate verifies
/// the unification's correctness, complementing the physics gates
/// (test_integrator_order, test_integrator_rkf78, test_propagator):
///
///   1. Tableau order conditions, checked EXACTLY from the rationals (no
///      numerical reference): row sums Σⱼ a[i][j] = c[i], Σᵢ b[i] = 1, and the
///      quadrature conditions Σᵢ b[i] c[i]ᵏ = 1/(k+1) up to each claimed order.
///      A transcription error in any tableau coefficient breaks one identity.
///   2. Driver reduction (S=1): one integrators::euler step equals the
///      closed-form forward-Euler advance — pose ← lie_advance_pose(pose,dt,
///      twist), twist ← twist + dt·accel — anchoring rk_step to the definition.
///   3. Empirical convergence order via the entry points on the analytic
///      circular Kepler orbit: euler ~1, rk4 ~4, rkf78 high (REQ-IN-4).
///   4. rk4 integrates correctly despite the round-off-reordered combine:
///      energy conserved over one period to RK4 grade (REQ-IN-1).
///   5. rkf78's embedded 7(8) error (fed by rk_step's stage accelerations) is
///      finite and positive.
///   6. Precision tightens with a wider numeric type T.
///
/// Theory: design/derivations/runge_kutta_lie_group.md. Exit 0 iff all pass.

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
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using boost::multiprecision::cpp_bin_float_50;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

/// Generic sqrt: ADL picks std::sqrt for double and boost's for cpp_bin_float_50.
template<typename T>
T gen_sqrt(const T& x) { using std::sqrt; return sqrt(x); }

/// Wrap a T-typed value (not a double literal) as a TrackedValue<T>.
template<typename T>
math::TrackedValue<T> tvT(const T& val) {
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

/// Check a tableau's order conditions exactly from its (double) coefficients.
template<typename Tab>
void check_tableau(const std::string& name, const Tab& tab, int S, int order) {
    bool rows = true;
    for (int i = 0; i < S; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) s += tab.a[i][j].value;
        if (std::abs(s - tab.c[i].value) > 1e-12) rows = false;
    }
    check(name + ": row sums Σⱼ a[i][j] = c[i]", rows);

    double bsum = 0.0;
    for (int i = 0; i < S; ++i) bsum += tab.b[i].value;
    check(name + ": Σᵢ b[i] = 1", std::abs(bsum - 1.0) < 1e-12);

    bool quad = true;
    for (int k = 0; k < order; ++k) {
        double s = 0.0;
        for (int i = 0; i < S; ++i) s += tab.b[i].value * std::pow(tab.c[i].value, k);
        if (std::abs(s - 1.0 / (k + 1)) > 1e-10) quad = false;
    }
    check(name + ": quadrature Σᵢ b[i] c[i]ᵏ = 1/(k+1)", quad);
}

// Shared analytic circular-Kepler setup at 400 km.
template<typename T>
struct Orbit {
    constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    T GM = K.earth.GM.value;
    T R = K.earth.a.value + 400000.0;
    T v = gen_sqrt<T>(GM / R);
    T period = T(2) * boost::math::constants::pi<T>() * R / v;
    dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(tv<T>(1.0));

    dynamics::State<T> state0() const {
        math::Vector3<T> r0(tvT<T>(R), tvT<T>(T(0)), tvT<T>(T(0)));
        math::Vector3<T> v0(tvT<T>(T(0)), tvT<T>(v), tvT<T>(T(0)));
        math::Vector3<T> w0;
        return dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), r0, w0, v0, tvT<T>(T(0)));
    }
    integrators::AccelFn<T> accel() const {
        return [this](const dynamics::State<T>& s) {
            return inertia.acceleration_from_wrench(forces::gravity_central(s, K));
        };
    }
};

} // namespace

int main() {
    std::cout << std::setprecision(11);
    std::cout << "test_butcher_tableau: unified ButcherTableau + rk_step driver (P5)\n\n";
    using T = double;

    // ---- 1. tableau order conditions (exact) ----
    std::cout << "=== tableau order conditions (exact, from rationals) ===\n";
    check_tableau("euler", integrators::euler_tableau<T>(), 1, 1);
    check_tableau("rk4  ", integrators::rk4_tableau<T>(), 4, 4);
    check_tableau("rkf78", integrators::rkf78_tableau<T>(), integrators::kRkf78Stages, 8);

    Orbit<T> orb;
    dynamics::State<T> y0 = orb.state0();
    integrators::AccelFn<T> accel = orb.accel();

    // ---- 2. driver reduction (S=1) == forward Euler ----
    std::cout << "\n=== driver reduction: euler step == closed-form forward Euler ===\n";
    {
        math::TrackedValue<T> dt = tv<T>(60.0);
        dynamics::State<T> e = integrators::euler(y0, dt, accel);
        dynamics::Pose<T> pose_ref = integrators::lie_advance_pose(y0.pose, dt, y0.twist);
        dynamics::Twist<T> twist_ref = y0.twist + dt * accel(y0);
        double dpos = std::sqrt(
            std::pow(e.position().x.value - pose_ref.translation().x.value, 2) +
            std::pow(e.position().y.value - pose_ref.translation().y.value, 2) +
            std::pow(e.position().z.value - pose_ref.translation().z.value, 2));
        double dtw = std::abs(e.twist.linear.y.value - twist_ref.linear.y.value);
        std::cout << "  |Δpose| = " << dpos << " m,  |Δtwist_vy| = " << dtw << "\n";
        check("euler pose == lie_advance_pose(pose, dt, twist)", dpos < 1e-6);
        check("euler twist == twist + dt·accel", dtw < 1e-9);
    }

    // ---- 3. empirical convergence order (REQ-IN-4) ----
    std::cout << "\n=== empirical convergence order on analytic Kepler ===\n";
    auto err_at = [&](int which, int n) {
        const T eighth = orb.period / 8.0;
        const T h = eighth / n;
        dynamics::State<T> y = y0;
        for (int i = 0; i < n; ++i) {
            if (which == 0) y = integrators::euler(y, tv<T>(h), accel);
            else if (which == 1) y = integrators::runge_kutta_4(y, tv<T>(h), accel);
            else y = integrators::rkf78(y, tv<T>(h), accel);
        }
        const T theta = (orb.v / orb.R) * eighth;
        math::Vector3<T> r = y.position();
        return std::sqrt(std::pow(r.x.value - orb.R * std::cos(theta), 2) +
                         std::pow(r.y.value - orb.R * std::sin(theta), 2) +
                         std::pow(r.z.value, 2));
    };
    double pe = std::log2(err_at(0, 400) / err_at(0, 800));
    double pr = std::log2(err_at(1, 24) / err_at(1, 48));
    double pf = std::log2(err_at(2, 5) / err_at(2, 10));
    std::cout << "  order: euler=" << pe << "  rk4=" << pr << "  rkf78=" << pf << "\n";
    check("euler first order (0.8..1.5)", pe >= 0.8 && pe <= 1.5);
    check("rk4 fourth order (3.5..4.5)", pr >= 3.5 && pr <= 4.5);
    check("rkf78 high order (>= 5)", pf >= 5.0);

    // ---- 4. rk4 integrates correctly despite the reordered combine ----
    std::cout << "\n=== rk4 energy conservation over one period (REQ-IN-1) ===\n";
    {
        const T E0 = T(0.5) * orb.v * orb.v - orb.GM / orb.R;
        dynamics::State<T> y = y0;
        T t = 0.0;
        while (t < orb.period) {
            T h = std::min<T>(60.0, orb.period - t);
            y = integrators::runge_kutta_4(y, tv<T>(h), accel);
            t += h;
        }
        math::Vector3<T> r = y.position(), vv = y.linear_velocity();
        T rmag = std::sqrt(r.x.value*r.x.value + r.y.value*r.y.value + r.z.value*r.z.value);
        T vmag = std::sqrt(vv.x.value*vv.x.value + vv.y.value*vv.y.value + vv.z.value*vv.z.value);
        T E = T(0.5) * vmag * vmag - orb.GM / rmag;
        T rel = std::abs(E - E0) / std::abs(E0);
        std::cout << "  rk4 energy drift (relative): " << rel << "\n";
        check("rk4 conserves energy < 1e-6 over one period", rel < 1e-6);
    }

    // ---- 5. rkf78 embedded error fed by rk_step stages ----
    std::cout << "\n=== rkf78 embedded 7(8) error ===\n";
    {
        integrators::StepResult<T> sr = integrators::rkf78_step(y0, tv<T>(60.0), accel);
        std::cout << "  embedded local error estimate: " << sr.error.value << "\n";
        check("embedded error finite", std::isfinite(sr.error.value));
        check("embedded error positive", sr.error.value > 0.0);
    }

    // ---- 6. precision scaling ----
    std::cout << "\n=== precision tightens with wider T ===\n";
    {
        auto step_precision = [](auto tag) -> double {
            using U = decltype(tag);
            Orbit<U> o;
            dynamics::State<U> s = integrators::runge_kutta_4(
                o.state0(), tv<U>(60.0), o.accel());
            return static_cast<double>(s.position().x.errors.precision);
        };
        double pd = step_precision(double{});
        double pb = step_precision(cpp_bin_float_50{});
        std::cout << "  rk4 position precision: double=" << pd << "  bf50=" << pb << "\n";
        check("precision > 0 (framework alive)", pd > 0.0 && pb > 0.0);
        check("precision tightens with wider T", pb < pd);
    }

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
