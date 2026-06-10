/// test_api_parity (gate F2) — the DQSGP4 propagator presents the SAME interface
/// as the analytical SGP4 propagator: a TLE factory and a propagate(tsince_min)
/// verb, so the two models differ only in their mathematics, not their shape.
///
///   SGP4:    Propagator p(cfg, elems, tol);           sv = p.propagate(t_min);
///   DQSGP4:  auto p = DqSgp4Propagator::authentic(td, tol, dt);  s = p.propagate(t_min);
///
/// Checks:
///   A. the DQ epoch state reproduces SGP4 at t=0 (the seed IS the analytical
///      state, round-tripped through the F1 state bridge);
///   B. propagate(t>0) advances and stays finite (forward integration);
///   C. authentic vs boosted are both finite and distinct (the G1 mode duality);
///   D. the verb is symmetric — one tsince drives both propagators.

#include "dynamics/dq_sgp4_propagator.h"
#include "dynamics/propagatable.h"
#include "dynamics/state_conversion.h"
#include "sgp4/model_selector.h"
#include "sgp4/sgp4_propagator.h"
#include "tle/tle_parser.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using T = double;

namespace {

int failed = 0;

/// R3c — the Propagatable concept's generic consumer: ONE function drives any
/// propagator (dq_propagator_facade.md section 5). Returns the native state.
template<typename P>
    requires dynamics::Propagatable<P, T>
auto propagate_generic(const P& p, const math::TrackedValue<T>& tsince) {
    return p.propagate(tsince);
}

void check(const char* name, bool ok, double d) {
    if (ok) {
        std::cout << "  PASS: " << name << "  (" << std::setprecision(6) << d << ")\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "  (" << std::setprecision(17) << d << ")\n";
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
    std::cout << std::setprecision(11)
              << "test_api_parity (F2): SGP4 <-> DQSGP4 interface parity\n\n";
    const T tol = T(1e-12);

    // SGP4-VER satellite 00005.
    std::string l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    std::string l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tle::TleData td;
    if (!tle::parse(l1, l2, td)) {
        std::cerr << "TLE parse failed\n";
        return 1;
    }

    // SGP4 reference path (3-line setup).
    sgp4::ModelConfiguration<T> cfg = sgp4::ModelSelector<T>::select("sgp4_standard", tol);
    tle::TleElements<T> elems = tle::TleElements<T>::from_tle_data(td);
    sgp4::Propagator<T> sgp4_prop(cfg, elems, tol);

    // DQSGP4 path (1-line setup) — same verb.
    dynamics::DqSgp4Propagator<T> dq =
        dynamics::DqSgp4Propagator<T>::authentic(td, tol, math::exact<T>(30));

    // A. DQ epoch (via the F1 bridge) reproduces SGP4 at t=0.
    std::cout << "=== A. DQ epoch == SGP4 t=0 ===\n";
    sgp4::StateVector<T> sv0 = sgp4_prop.propagate(math::exact<T>(0));
    dynamics::State<T> dq0 = dq.propagate(math::exact<T>(0));
    sgp4::StateVector<T> dq0_km = dynamics::to_state_vector(dq0);
    double pos_err = dist(dq0_km.position_km, sv0.position_km);
    double vel_err = dist(dq0_km.velocity_km_s, sv0.velocity_km_s);
    check("epoch position matches SGP4 t=0 (km)", pos_err < 1e-6, pos_err);
    check("epoch velocity matches SGP4 t=0 (km/s)", vel_err < 1e-9, vel_err);

    // B. propagate(10 min) advances and stays finite.
    std::cout << "\n=== B. propagate(10 min) advances ===\n";
    dynamics::State<T> dq10 = dq.propagate(math::exact<T>(10));
    double moved_km = dist(dq10.position(), dq0.position()) / 1000.0;
    check("position advanced > 1 km", moved_km > 1.0, moved_km);
    check("position finite", std::isfinite(dq10.position().x.value),
          dq10.position().x.value);

    // C. authentic vs boosted: both finite, distinct (G1 mode duality).
    std::cout << "\n=== C. authentic vs boosted (G1 modes) ===\n";
    dynamics::DqSgp4Propagator<T> dqb =
        dynamics::DqSgp4Propagator<T>::boosted(td, tol, math::exact<T>(30));
    dynamics::State<T> dqb10 = dqb.propagate(math::exact<T>(10));
    check("boosted finite", std::isfinite(dqb10.position().x.value),
          dqb10.position().x.value);
    double mode_diff_m = dist(dqb10.position(), dq10.position());
    check("boosted differs from authentic (m)", mode_diff_m > 0.0, mode_diff_m);

    // D. Verb symmetry: one tsince drives both propagators to finite states.
    std::cout << "\n=== D. verb symmetry ===\n";
    math::TrackedValue<T> tsince = math::exact<T>(5);
    sgp4::StateVector<T> a = sgp4_prop.propagate(tsince);
    dynamics::State<T> b = dq.propagate(tsince);
    check("both propagate(tsince) return finite states",
          std::isfinite(a.position_km.x.value) && std::isfinite(b.position().x.value),
          0.0);

    // E. R3c — the Propagatable concept, consumed generically: the SAME generic
    // function drives BOTH propagator families and reproduces the direct calls
    // bit-for-bit (the concept names a real shared verb, not scaffolding).
    std::cout << "\n=== E. Propagatable concept (generic driver) ===\n";
    static_assert(dynamics::Propagatable<sgp4::Propagator<T>, T>);
    static_assert(dynamics::Propagatable<dynamics::DqSgp4Propagator<T>, T>);
    sgp4::StateVector<T> ga = propagate_generic(sgp4_prop, tsince);
    dynamics::State<T> gb = propagate_generic(dq, tsince);
    check("generic driver == direct SGP4 call (bit)",
          ga.position_km.x.value == a.position_km.x.value &&
          ga.position_km.y.value == a.position_km.y.value &&
          ga.position_km.z.value == a.position_km.z.value, 0.0);
    check("generic driver == direct DQSGP4 call (bit)",
          gb.position().x.value == b.position().x.value &&
          gb.position().y.value == b.position().y.value &&
          gb.position().z.value == b.position().z.value, 0.0);

    std::cout << "\n" << (failed == 0 ? "PASS" : "FAIL") << " — " << failed
              << " failure(s)\n";
    return failed == 0 ? 0 : 1;
}
