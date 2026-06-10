/// test_gravity_zonal_jn (gate D1) — the generic gravity_zonal force.
///
/// Validates the closed-form Jₙ zonal sum that consumes the shared generative
/// ZonalHarmonics provider (GAL1), so the geopotential is no longer pinned at
/// J₂-only:
///   A. max_n = 2 reproduces gravity_J2 (the regression anchor);
///   B. J₃ breaks north–south symmetry — a NONZERO a_z at the equator, where
///      the J₂-only field is exactly zero;
///   C. the truncation residual SHRINKS monotonically as more zonals are
///      modelled (it is the real omitted-tail magnitude);
///   D. a full J₂…J₉ (EGM2008) sum is finite and nonzero.
///
/// Pass/fail convention mirrors tests/test_force_models/main.cpp.

#include "constants/constants_provider.h"
#include "constants/zonal_harmonics.h"
#include "dynamics/state.h"
#include "dynamics/wrench.h"
#include "forces/gravity_zonal.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

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

/// TrackedValue with zero measurement/accuracy (only the representation bound),
/// so any accuracy in a force output is attributable to the force model.
TV tv(double v) {
    T val = static_cast<T>(v);
    return TV(val, T(0), TV::representation_bound(val), T(0));
}

void check(const char* name, bool ok, double detail) {
    if (ok) {
        ++passed;
        std::cout << "  PASS: " << name << "  (" << std::setprecision(6) << detail << ")\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "  (" << std::setprecision(17) << detail << ")\n";
    }
}

dynamics::State<T> make_state(const Vector3<T>& r) {
    Vector3<T> v(tv(0.0), tv(7000.0), tv(1000.0));
    Vector3<T> omega_body;
    Quaternion<T> q = Quaternion<T>::identity();
    return dynamics::State<T>::from_kinematics(q, r, omega_body, v, tv(0.0));
}

double rel(double a, double b) {
    return std::abs(a - b) / std::max(std::abs(b), 1e-300);
}

} // namespace

int main() {
    std::cout << std::setprecision(11)
              << "test_gravity_zonal_jn (D1): generic Jn zonal gravity\n\n";

    constants::ConstantsProvider<T> K72 = constants::ConstantsProvider<T>::wgs72(T(1e-12));
    constants::ZonalHarmonics<T> zh72 = constants::ZonalHarmonics<T>::wgs72(T(1e-12));
    const T r_mag = K72.earth.a.value + T(700000.0);   // 700 km altitude
    const T c = r_mag / std::sqrt(T(3));

    // A. max_n = 2 reproduces gravity_J2 (force value).
    std::cout << "=== A. gravity_zonal(.,2) reproduces gravity_J2 ===\n";
    {
        dynamics::State<T> s = make_state(Vector3<T>(tv(c), tv(c), tv(c)));
        Vector3<T> g2 = forces::gravity_J2(s, K72).force;
        Vector3<T> gz = forces::gravity_zonal(s, K72, zh72, 2).force;
        check("a_x matches gravity_J2 (rel < 1e-10)", rel(gz.x.value, g2.x.value) < 1e-10,
              rel(gz.x.value, g2.x.value));
        check("a_z matches gravity_J2 (rel < 1e-10)", rel(gz.z.value, g2.z.value) < 1e-10,
              rel(gz.z.value, g2.z.value));
    }

    // B. J₃ north–south asymmetry: nonzero a_z at the equator (z = 0).
    std::cout << "\n=== B. J3 north-south asymmetry (equatorial a_z) ===\n";
    {
        dynamics::State<T> s = make_state(Vector3<T>(tv(r_mag), tv(0.0), tv(0.0)));
        T az2 = forces::gravity_zonal(s, K72, zh72, 2).force.z.value;   // J₂ only
        T az3 = forces::gravity_zonal(s, K72, zh72, 3).force.z.value;   // + J₃
        check("J2-only a_z == 0 at equator", std::abs(az2) < 1e-12, az2);
        check("J2+J3 a_z != 0 at equator (J3 breaks symmetry)", std::abs(az3) > 1e-9, az3);
    }

    // C. Truncation residual shrinks monotonically with max_n (EGM2008: J₂…J₉).
    std::cout << "\n=== C. truncation residual shrinks with max_n ===\n";
    {
        constants::ConstantsProvider<T> K84 = constants::ConstantsProvider<T>::wgs84(T(1e-12));
        constants::ZonalHarmonics<T> zh = constants::ZonalHarmonics<T>::egm2008(T(1e-12));
        dynamics::State<T> s = make_state(Vector3<T>(tv(c), tv(c), tv(c)));
        T acc2 = forces::gravity_zonal(s, K84, zh, 2).force.x.errors.accuracy;
        T acc4 = forces::gravity_zonal(s, K84, zh, 4).force.x.errors.accuracy;
        T acc8 = forces::gravity_zonal(s, K84, zh, 8).force.x.errors.accuracy;
        check("residual(4) < residual(2)", acc4 < acc2, acc4 / acc2);
        check("residual(8) < residual(4)", acc8 < acc4, acc8 / acc4);
        check("residual(8) > 0", acc8 > T(0), acc8);
    }

    // D. Full J₂…J₉ sum is finite and nonzero.
    std::cout << "\n=== D. full J2..J9 (EGM2008) sum finite & nonzero ===\n";
    {
        constants::ConstantsProvider<T> K84 = constants::ConstantsProvider<T>::wgs84(T(1e-12));
        constants::ZonalHarmonics<T> zh = constants::ZonalHarmonics<T>::egm2008(T(1e-12));
        dynamics::State<T> s = make_state(Vector3<T>(tv(c), tv(c), tv(c)));
        Vector3<T> g = forces::gravity_zonal(s, K84, zh, 9).force;
        double m = std::sqrt(g.x.value * g.x.value + g.y.value * g.y.value
                                                   + g.z.value * g.z.value);
        check("|a(J2..J9)| finite", std::isfinite(m), m);
        check("|a(J2..J9)| > 0", m > T(0), m);
    }

    std::cout << "\n" << (failed == 0 ? "PASS" : "FAIL") << " — " << passed
              << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
