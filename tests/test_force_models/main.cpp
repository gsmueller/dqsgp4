/// test_force_models — AUD-EF-6 / REQ-EF-6 behavioral test of the force
/// lambdas' model-truncation accuracy budget.
///
/// Contract under test (AUD-EF-6 / REQ-EF-6). Every force lambda records its
/// model-truncation residual in the `errors.accuracy` category of the wrench
/// it returns:
///
///   gravity_central — the spherically symmetric monopole −GM·r/|r|³ is EXACT
///                     (Newton's shell theorem). Its model-truncation residual
///                     is identically ZERO and is recorded as such.
///   gravity_J2      — omits the J₃, J₄, … zonals, so it carries a NONZERO
///                     residual (a fraction of |a_J2|, per gravity_zonal.h).
///   drag            — the exponential-atmosphere (Lane 1965) density model has
///                     a NONZERO residual (~30 % baseline of |a_drag|).
///
/// This is a direct unit test of the three force functions (no propagator):
/// it builds one LEO `State<T>`, evaluates each force, and asserts the
/// `accuracy` component of the returned `Wrench`'s force vector. Inputs are
/// constructed with zero accuracy (WGS 84 `GM` is `measured()` → accuracy 0;
/// the position is built with accuracy 0), and the orientation is identity so
/// the body frame coincides with the world frame and the accuracy budget
/// carries through the inverse rotation unchanged. Hence `gravity_central`'s
/// residual is provably *exactly* zero, not merely small, while the J₂ and
/// drag lambdas add their residual explicitly to every component.
///
/// Pass/fail and exit convention mirror tests/test_propagator/main.cpp:
/// per-check counters, nonzero exit code on any failure.

#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "dynamics/wrench.h"
#include "forces/drag.h"
#include "forces/gravity_central.h"
#include "forces/gravity_zonal.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

using T = double;
using TV = math::TrackedValue<T>;
using math::Quaternion;
using math::Vector3;

namespace {

int tests_passed = 0;
int tests_failed = 0;

/// Build a TrackedValue with measurement = 0 and accuracy = 0 (only the
/// unavoidable binary-representation precision bound). A clean input: any
/// accuracy in the output is therefore attributable to the force model.
TV tv(double v) {
    T val = static_cast<T>(v);
    T prec = TV::representation_bound(val);
    return TV(val, T(0), prec, T(0));
}

/// AUD-EF-6: assert a force component's model-truncation residual is exactly
/// zero (the exact-monopole case).
void check_zero_accuracy(const char* name, T accuracy) {
    if (accuracy == T(0)) {
        ++tests_passed;
        std::cout << "  PASS: " << name
                  << "  accuracy=" << std::setprecision(6) << accuracy << "\n";
    } else {
        ++tests_failed;
        std::cerr << "  FAIL: " << name
                  << "  accuracy=" << std::setprecision(17) << accuracy
                  << "  (expected exactly 0)\n";
    }
}

/// AUD-EF-6: assert a force component records a strictly positive
/// model-truncation residual (the truncated/approximate-model case).
void check_positive_accuracy(const char* name, T accuracy) {
    if (accuracy > T(0)) {
        ++tests_passed;
        std::cout << "  PASS: " << name
                  << "  accuracy=" << std::setprecision(6) << accuracy << "\n";
    } else {
        ++tests_failed;
        std::cerr << "  FAIL: " << name
                  << "  accuracy=" << std::setprecision(17) << accuracy
                  << "  (expected > 0)\n";
    }
}

/// Sanity guard: assert the force itself is nonzero, so the accuracy checks
/// above are not trivially satisfied by an all-zero wrench.
void check_force_nonzero(const char* name, const Vector3<T>& f) {
    T mag = std::sqrt(f.x.value * f.x.value + f.y.value * f.y.value
                                            + f.z.value * f.z.value);
    if (mag > T(0)) {
        ++tests_passed;
        std::cout << "  PASS: " << name
                  << "  |F|=" << std::setprecision(6) << mag << "\n";
    } else {
        ++tests_failed;
        std::cerr << "  FAIL: " << name << "  |F|=0 (force vanished)\n";
    }
}

} // anonymous namespace

int main() {
    std::cout << std::setprecision(11);
    std::cout << "test_force_models: AUD-EF-6 / REQ-EF-6 accuracy-budget contract\n\n";

    // --- Constants provider (WGS 84) ---
    // GM is measured() → accuracy 0, so gravity_central's residual is exactly
    // zero rather than merely small.
    constants::ConstantsProvider<T> K =
        constants::ConstantsProvider<T>::wgs84(T(1e-12));

    // --- LEO state at 400 km altitude ---
    // Position with all three components equal and nonzero (r/√3 each) so the
    // J₂ and drag accelerations are nonzero in every component; the recorded
    // residual is then exercised on x, y, and z independently. Identity
    // orientation ⇒ body frame == world frame ⇒ the inverse rotation carries
    // the accuracy budget through unchanged.
    const T r_mag = K.earth.a.value + T(400000.0);
    const T comp  = r_mag / std::sqrt(T(3.0));
    Vector3<T> r(tv(comp), tv(comp), tv(comp));
    // Any velocity giving nonzero v_rel = v − ω_E×r in all components; the
    // exact value is immaterial to the residual contract.
    Vector3<T> v(tv(0.0), tv(7600.0), tv(1000.0));
    Vector3<T> omega_body;  // body angular velocity: irrelevant to these forces
    Quaternion<T> q = Quaternion<T>::identity();
    dynamics::State<T> state =
        dynamics::State<T>::from_kinematics(q, r, omega_body, v, tv(0.0));

    std::cout << "LEO state: |r| = " << r_mag << " m (400 km alt), "
              << "identity orientation\n\n";

    // ============================================================
    // gravity_central — exact monopole → ZERO model-truncation residual
    // ============================================================
    std::cout << "=== gravity_central (exact monopole) ===\n";
    dynamics::Wrench<T> w_cen = forces::gravity_central(state, K);
    check_force_nonzero("gravity_central produces a nonzero force", w_cen.force);
    check_zero_accuracy("AUD-EF-6 gravity_central force.x accuracy == 0",
                        w_cen.force.x.errors.accuracy);
    check_zero_accuracy("AUD-EF-6 gravity_central force.y accuracy == 0",
                        w_cen.force.y.errors.accuracy);
    check_zero_accuracy("AUD-EF-6 gravity_central force.z accuracy == 0",
                        w_cen.force.z.errors.accuracy);

    // ============================================================
    // gravity_J2 — omits J₃+ → NONZERO model-truncation residual
    // ============================================================
    std::cout << "\n=== gravity_J2 (omits J3, J4, ...) ===\n";
    dynamics::Wrench<T> w_j2 = forces::gravity_J2(state, K);
    check_force_nonzero("gravity_J2 produces a nonzero force", w_j2.force);
    check_positive_accuracy("AUD-EF-6 gravity_J2 force.x accuracy > 0",
                            w_j2.force.x.errors.accuracy);
    check_positive_accuracy("AUD-EF-6 gravity_J2 force.y accuracy > 0",
                            w_j2.force.y.errors.accuracy);
    check_positive_accuracy("AUD-EF-6 gravity_J2 force.z accuracy > 0",
                            w_j2.force.z.errors.accuracy);

    // ============================================================
    // drag — exponential-atmosphere model → NONZERO residual
    // ============================================================
    std::cout << "\n=== drag (exponential atmosphere, Lane 1965) ===\n";
    // rho_0 = 2.789e-10 kg/m³ at h_0 = 200 km, H_scale = 50 km,
    // B = C_d·A/m = 0.5 m²/kg — same parameters as test_propagator Phase 3.
    std::function<dynamics::Wrench<T>(const dynamics::State<T>&,
                                      const constants::ConstantsProvider<T>&)>
        drag_fn = forces::make_drag_exponential<T>(
            tv(2.789e-10), tv(200000.0), tv(50000.0), tv(0.5));
    dynamics::Wrench<T> w_drag = drag_fn(state, K);
    check_force_nonzero("drag produces a nonzero force", w_drag.force);
    check_positive_accuracy("AUD-EF-6 drag force.x accuracy > 0",
                            w_drag.force.x.errors.accuracy);
    check_positive_accuracy("AUD-EF-6 drag force.y accuracy > 0",
                            w_drag.force.y.errors.accuracy);
    check_positive_accuracy("AUD-EF-6 drag force.z accuracy > 0",
                            w_drag.force.z.errors.accuracy);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << tests_passed
              << "  Failed: " << tests_failed << "\n";
    std::cout << "========================================\n";
    return tests_failed > 0 ? 1 : 0;
}
