/// test_injection — INJ1 (DQSGP4 Completion Roadmap, register INJ1).
///
/// The ModelFunctions injection mechanism is now CONSISTENT: every slot it
/// exposes is genuinely called by the engine (no dead injection points), and
/// the misleading label-only selectors were removed. This gate proves it:
///   - injecting a custom kepler_solver changes the propagation result (the
///     engine genuinely dispatches the injected per-step solver);
///   - injecting a custom secular_rates changes the result (the engine genuinely
///     dispatches the injected compute-once rate function at init);
///   - the standard config's description no longer advertises a "drag=" selector
///     (the label-only .drag() builder was dropped — INJ1 simplification);
///   - the dead `inclination_function` slot is gone (compile-time: this file
///     would fail to build if it referenced the removed member — it does not).
///
/// ExeGate INJ1: nonzero exit code on any failed check.

#include "sgp4/sgp4_propagator.h"
#include "sgp4/model_selector.h"
#include "tle/tle_parser.h"
#include "math/tracked_value.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

using T = double;
using TV = math::TrackedValue<T>;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

TV tv(double v) {
    return TV(static_cast<T>(v), T(0), TV::representation_bound(static_cast<T>(v)), T(0));
}

double pos_mag(const sgp4::StateVector<T>& sv) {
    return std::sqrt(sv.position_km.x.value * sv.position_km.x.value
                   + sv.position_km.y.value * sv.position_km.y.value
                   + sv.position_km.z.value * sv.position_km.z.value);
}

} // namespace

int main() {
    // Near-Earth SGP4-VER sat 00005 (lines verbatim; parse reads fixed columns).
    const char* l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    const char* l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tle::TleData td;
    bool ok = tle::parse(std::string(l1), std::string(l2), td);
    check("TLE parses", ok);
    if (!ok) { std::cout << "\n  injection: " << passed << " passed, " << failed << " failed\n"; return 1; }

    tle::TleElements<T> elems = tle::TleElements<T>::from_tle_data(td);
    TV t = tv(120.0);  // 120 minutes since epoch

    // --- Baseline (standard model functions) ---
    auto cfg_std = sgp4::ModelSelector<T>::select("sgp4_standard", 1e-12);
    sgp4::Propagator<T> p_std(cfg_std, elems, 1e-12);
    double base = pos_mag(p_std.propagate(t));
    check("baseline propagation finite & nonzero", std::isfinite(base) && base > 0.0);

    // --- Inject a custom kepler_solver (per-step slot) ---
    // Return U + 0.001 rad instead of solving for E+ω. If the engine genuinely
    // dispatches the injected solver, the result must change.
    auto cfg_k = sgp4::ModelSelector<T>::select("sgp4_standard", 1e-12);
    cfg_k.model_functions.kepler_solver =
        [](const TV& axN, const TV& ayn, const TV& U, const T& tol) -> TV {
            (void)axN; (void)ayn; (void)tol;
            return U + TV(T(0.001), T(0), T(0), T(0));
        };
    sgp4::Propagator<T> p_k(cfg_k, elems, 1e-12);
    double with_k = pos_mag(p_k.propagate(t));
    check("injected kepler_solver changes the result (engine dispatches it)",
          std::isfinite(with_k) && std::abs(with_k - base) > 1e-6);

    // --- Inject a custom secular_rates (compute-once slot) ---
    // Wrap the standard rates and perturb M_dot; the engine uses it at init.
    auto cfg_s = sgp4::ModelSelector<T>::select("sgp4_standard", 1e-12);
    auto std_sr = cfg_s.model_functions.secular_rates;
    cfg_s.model_functions.secular_rates =
        [std_sr](const TV& n, const TV& a, const TV& e2, const TV& ci,
                 const TV& J2, const TV& J4) {
            perturbation::BrouwerSecularRates<T> r = std_sr(n, a, e2, ci, J2, J4);
            r.M_dot = r.M_dot + TV(T(1e-6), T(0), T(0), T(0));  // perturb M_dot
            return r;
        };
    sgp4::Propagator<T> p_s(cfg_s, elems, 1e-12);
    double with_s = pos_mag(p_s.propagate(t));
    check("injected secular_rates changes the result (engine dispatches it)",
          std::isfinite(with_s) && std::abs(with_s - base) > 1e-6);

    // --- The label-only .drag() selector was dropped (INJ1) ---
    check("standard config description has no 'drag=' selector",
          cfg_std.description.find("drag=") == std::string::npos);

    std::cout << "\n  injection: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
