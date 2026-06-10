/// test_force_presets (gate FM1, R3a) — the DQ facade's force injection + presets.
///
/// Verifies design/derivations/dq_propagator_facade.md §2/§3/§6:
///   1. REGRESSION: from_tle with default DqForceOptions is BIT-IDENTICAL to the
///      prior gravitational-only facade (authentic()) over a full propagation.
///   2. The seed (epoch state) is option-independent.
///   3. Each preset toggles a real perturbation: lunisolar / Vallado-8-4 drag /
///      SRP each diverge from the default model; lunisolar dominates SRP (the
///      tidal accelerations ~2e-6 vs the cannonball ~1.2e-7 m/s² here).
///   4. THE DOC-LIE FIXED: arbitrary extra_forces injection through the explicit
///      constructor — a counting lambda is invoked exactly steps × 4 RK4 stages
///      times, and a constant 1e-6 m/s² test force displaces the orbit by the
///      ½at² textbook value.
///   5. Composed options (all three presets at once) propagate finite.
///
/// Exit 0 iff every check passes.

#include "dynamics/dq_sgp4_propagator.h"
#include "math/tracked_value.h"
#include "math/vector3.h"
#include "tle/tle_parser.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using T = double;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

double dist(const math::Vector3<T>& a, const math::Vector3<T>& b) {
    double dx = a.x.value - b.x.value;
    double dy = a.y.value - b.y.value;
    double dz = a.z.value - b.z.value;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

int main() {
    std::cout << std::setprecision(6);
    std::cout << "test_force_presets (FM1): DQ facade force injection + presets\n\n";
    const T tol = T(1e-12);
    const math::TrackedValue<T> dt = math::exact<T>(30);
    const math::TrackedValue<T> t60 = math::exact<T>(60);   // minutes

    std::string l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    std::string l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tle::TleData td;
    if (!tle::parse(l1, l2, td)) {
        std::cerr << "TLE parse failed\n";
        return 1;
    }
    using DQ = dynamics::DqSgp4Propagator<T>;
    using Opt = dynamics::DqForceOptions<T>;
    const auto mode = dynamics::PropagatorMode::Authentic;

    // ---- 1. default options == the prior facade, bit-identical ----
    std::cout << "=== regression: default options == prior facade ===\n";
    DQ prior = DQ::authentic(td, tol, dt);
    DQ deflt = DQ::from_tle(td, tol, mode, dt, Opt{});
    dynamics::State<T> s_prior = prior.propagate(t60);
    dynamics::State<T> s_deflt = deflt.propagate(t60);
    {
        bool bit = s_prior.position().x.value == s_deflt.position().x.value &&
                   s_prior.position().y.value == s_deflt.position().y.value &&
                   s_prior.position().z.value == s_deflt.position().z.value;
        std::cout << "  |d| = " << dist(s_prior.position(), s_deflt.position()) << " m\n";
        check("default-options propagation BIT-IDENTICAL to authentic()", bit);
        check("seed state option-independent",
              prior.epoch_state().position().x.value ==
              deflt.epoch_state().position().x.value);
    }

    // ---- 2/3. each preset toggles a real perturbation ----
    std::cout << "\n=== presets perturb ===\n";
    {
        Opt o_luni; o_luni.lunisolar = true;
        Opt o_drag; o_drag.drag_B = math::TrackedValue<T>(0.5, T(0), T(0), T(0));
        Opt o_srp;  o_srp.srp_cr_area_over_mass =
            math::TrackedValue<T>(0.026, T(0), T(0), T(0));

        dynamics::State<T> s_luni = DQ::from_tle(td, tol, mode, dt, o_luni).propagate(t60);
        dynamics::State<T> s_drag = DQ::from_tle(td, tol, mode, dt, o_drag).propagate(t60);
        dynamics::State<T> s_srp  = DQ::from_tle(td, tol, mode, dt, o_srp).propagate(t60);

        double d_luni = dist(s_luni.position(), s_deflt.position());
        double d_drag = dist(s_drag.position(), s_deflt.position());
        double d_srp  = dist(s_srp.position(),  s_deflt.position());
        std::cout << "  divergence over 60 min:  lunisolar " << d_luni
                  << " m   drag(B=0.5) " << d_drag << " m   SRP(0.026) " << d_srp << " m\n";
        check("lunisolar preset perturbs (> 0.05 m / 60 min)", d_luni > 0.05);
        check("drag preset perturbs (> 0.05 m / 60 min)", d_drag > 0.05);
        check("SRP preset perturbs (> 0.01 m / 60 min)", d_srp > 0.01);
        check("lunisolar dominates SRP here (tidal ~2e-6 vs cannonball ~1.2e-7)",
              d_luni > d_srp);
    }

    // ---- 4. arbitrary extra_forces injection (the doc-lie fixed) ----
    std::cout << "\n=== extra_forces injection ===\n";
    {
        // Counting lambda: zero wrench, counts invocations. 60 s at dt_max=30 s
        // = 2 RK4 steps x 4 stages = 8 acceleration evaluations.
        auto count = std::make_shared<int>(0);
        std::vector<dynamics::ForceFn<T>> probe;
        probe.push_back([count](const dynamics::State<T>&,
                                const constants::ConstantsProvider<T>&) {
            ++(*count);
            return dynamics::Wrench<T>::zero();
        });
        dynamics::State<T> seed = dynamics::state_from_tle<T>(td, tol);
        DQ counted(seed, constants::ConstantsProvider<T>::wgs72(tol),
                   constants::ZonalHarmonics<T>::wgs72(tol), 4, dt, std::move(probe));
        counted.propagate(math::exact<T>(1));   // 1 minute = 60 s
        std::cout << "  counting force invoked " << *count << " times (expect 8)\n";
        check("injected force evaluated at every RK4 stage (2 steps x 4)", *count == 8);

        // Constant 1e-6 m/s^2 along body-x for 60 s: displacement ~ half*a*t^2.
        std::vector<dynamics::ForceFn<T>> push;
        push.push_back([](const dynamics::State<T>&,
                          const constants::ConstantsProvider<T>&) {
            math::Vector3<T> f(math::TrackedValue<T>(1e-6, T(0), T(0), T(0)),
                               math::exact<T>(0), math::exact<T>(0));
            return dynamics::Wrench<T>(math::Vector3<T>(), f);
        });
        DQ pushed(seed, constants::ConstantsProvider<T>::wgs72(tol),
                  constants::ZonalHarmonics<T>::wgs72(tol), 4, dt, std::move(push));
        DQ plain(seed, constants::ConstantsProvider<T>::wgs72(tol),
                 constants::ZonalHarmonics<T>::wgs72(tol), 4, dt);
        dynamics::State<T> sp = pushed.propagate(math::exact<T>(1));
        dynamics::State<T> s0 = plain.propagate(math::exact<T>(1));
        double d = dist(sp.position(), s0.position());
        double expect = 0.5 * 1e-6 * 60.0 * 60.0;   // half*a*t^2 = 1.8e-3 m
        std::cout << "  constant-force displacement " << d << " m (half*a*t^2 = "
                  << expect << ")\n";
        check("constant test force gives the half*a*t^2 displacement (+-15%)",
              std::abs(d - expect) < 0.15 * expect);
    }

    // ---- 5. composed options ----
    std::cout << "\n=== composed presets ===\n";
    {
        Opt all;
        all.lunisolar = true;
        all.drag_B = math::TrackedValue<T>(0.5, T(0), T(0), T(0));
        all.srp_cr_area_over_mass = math::TrackedValue<T>(0.026, T(0), T(0), T(0));
        dynamics::State<T> s_all = DQ::from_tle(td, tol, mode, dt, all).propagate(t60);
        double d = dist(s_all.position(), s_deflt.position());
        bool finite = std::isfinite(s_all.position().x.value) &&
                      std::isfinite(s_all.position().y.value) &&
                      std::isfinite(s_all.position().z.value);
        std::cout << "  all-presets divergence " << d << " m\n";
        check("composed presets propagate finite and perturb", finite && d > 0.05);
    }

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
