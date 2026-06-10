/// test_constants_swap — swappable constants (OBJ-5 / REQ-SY-9).
///
/// The propagator and force models read every fundamental constant through the
/// injected ConstantsProvider; there is no global default. Swapping the
/// provider (WGS84 / WGS72 / GRS80) must therefore change the computed result,
/// and the change must track the swapped constant exactly. This test:
///   - confirms the three named conventions carry distinct defining constants
///     (Earth radius, GM, J2);
///   - confirms gravity_central reads GM from the provider — at a fixed
///     position the force ratio equals the GM ratio exactly (the force is
///     linear in the provider's GM, with nothing hardcoded);
///   - confirms gravity_J2 reads J2 from the provider (the J2 force differs);
///   - confirms the providers are independent values with no shared global
///     state (CON-2): building one does not perturb another.
///
/// Exit 0 iff every check passes.

#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "forces/gravity_central.h"
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

void check(const char* name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

TV tv(double v) { return TV(v, T(0), TV::representation_bound(v), T(0)); }

// A fixed position, independent of any provider, so a force comparison isolates
// the swapped constant.
dynamics::State<T> fixed_state() {
    Vector3<T> r(tv(7000000.0), tv(100000.0), tv(50000.0));
    Vector3<T> v(tv(0.0), tv(0.0), tv(0.0));
    Vector3<T> w0;
    return dynamics::State<T>::from_kinematics(
        Quaternion<T>::identity(), r, w0, v, tv(0.0));
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(12);
    std::cout << "test_constants_swap: swappable ConstantsProvider (REQ-SY-9)\n\n";

    constants::ConstantsProvider<T> wgs84 = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    constants::ConstantsProvider<T> wgs72 = constants::ConstantsProvider<T>::wgs72(T(1e-12));
    constants::ConstantsProvider<T> grs80 = constants::ConstantsProvider<T>::grs80(T(1e-12));

    const T a84 = wgs84.earth.a.value, a72 = wgs72.earth.a.value;
    const T gm84 = wgs84.earth.GM.value, gm72 = wgs72.earth.GM.value, gm80 = grs80.earth.GM.value;
    const T j2_84 = wgs84.earth.J2n(1).value, j2_72 = wgs72.earth.J2n(1).value;

    std::cout << "=== defining constants ===\n";
    std::cout << "  a  : wgs84=" << a84 << "  wgs72=" << a72 << "\n";
    std::cout << "  GM : wgs84=" << gm84 << "  wgs72=" << gm72 << "  grs80=" << gm80 << "\n";
    std::cout << "  J2 : wgs84=" << j2_84 << "  wgs72=" << j2_72 << "\n";

    check("WGS84 and WGS72 Earth radius differ (6378137 vs 6378135)",
          std::abs(a84 - 6378137.0) < 1e-6 && std::abs(a72 - 6378135.0) < 1e-6 && a84 != a72);
    check("WGS84 / WGS72 / GRS80 GM are all distinct",
          gm84 != gm72 && gm72 != gm80 && gm84 != gm80);
    check("WGS84 and WGS72 J2 differ", j2_84 != j2_72);

    // --- gravity_central reads GM from the provider ---
    std::cout << "\n=== gravity_central reads GM from the provider ===\n";
    dynamics::State<T> s = fixed_state();
    const T f84 = forces::gravity_central(s, wgs84).force.x.value;
    const T f72 = forces::gravity_central(s, wgs72).force.x.value;
    const T force_ratio = f72 / f84;
    const T gm_ratio = gm72 / gm84;
    std::cout << "  force.x: wgs84=" << f84 << "  wgs72=" << f72 << "\n";
    std::cout << "  force ratio=" << force_ratio << "  GM ratio=" << gm_ratio << "\n";
    check("central-gravity force changes when the provider changes", f72 != f84);
    check("force ratio equals the GM ratio exactly (GM flows through linearly)",
          std::abs(force_ratio - gm_ratio) < 1e-12);

    // --- gravity_J2 reads J2 from the provider ---
    std::cout << "\n=== gravity_J2 reads J2 from the provider ===\n";
    const T j84 = forces::gravity_J2(s, wgs84).force.x.value;
    const T j72 = forces::gravity_J2(s, wgs72).force.x.value;
    std::cout << "  J2 force.x: wgs84=" << j84 << "  wgs72=" << j72 << "\n";
    check("J2 force changes when the provider changes", j84 != j72);

    // --- no global state: providers are independent values ---
    std::cout << "\n=== no global state (CON-2) ===\n";
    constants::ConstantsProvider<T> wgs84_again = constants::ConstantsProvider<T>::wgs84(T(1e-12));
    check("rebuilding a provider reproduces identical constants (pure, no global state)",
          wgs84_again.earth.GM.value == gm84 && wgs84.earth.GM.value == gm84);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
