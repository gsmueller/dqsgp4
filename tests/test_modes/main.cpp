/// test_modes (gate G1) — the authentic vs boosted fidelity-mode duality.
///
///   - DQ authentic uses the WGS72 ellipsoid (a = 6378135 m); boosted uses WGS84
///     (a = 6378137 m);
///   - BOTH seed from the authentic WGS72 SGP4, so the epoch state is identical
///     (a TLE is WGS72 by definition — only the propagation model differs);
///   - the two modes diverge under propagation;
///   - from_tle(mode) dispatches to the correct model by a named parameter;
///   - the SGP4 authentic preset (sgp4_standard) is WGS72.

#include "dynamics/dq_sgp4_propagator.h"
#include "sgp4/model_selector.h"
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

void check(const char* name, bool ok, double detail) {
    if (ok) {
        std::cout << "  PASS: " << name << "  (" << std::setprecision(6) << detail << ")\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "  (" << std::setprecision(17) << detail << ")\n";
    }
}

void check_close(const char* name, double got, double want, double tol) {
    if (std::abs(got - want) <= tol) {
        std::cout << "  PASS: " << name << " = " << std::setprecision(10) << got << "\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << " = " << std::setprecision(17) << got
                  << " (want " << want << ")\n";
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
    std::cout << "test_modes (G1): authentic vs boosted mode duality\n\n";
    const T tol = T(1e-12);
    const math::TrackedValue<T> dt = math::exact<T>(30);

    std::string l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    std::string l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tle::TleData td;
    if (!tle::parse(l1, l2, td)) {
        std::cerr << "TLE parse failed\n";
        return 1;
    }

    // Named-mode factories map to the right earth model.
    std::cout << "=== mode -> earth model ===\n";
    dynamics::DqSgp4Propagator<T> A = dynamics::DqSgp4Propagator<T>::authentic(td, tol, dt);
    dynamics::DqSgp4Propagator<T> B = dynamics::DqSgp4Propagator<T>::boosted(td, tol, dt);
    check_close("authentic earth = WGS72 a (m)", A.propagator().constants().earth.a.value,
                6378135.0, 1.0);
    check_close("boosted earth = WGS84 a (m)", B.propagator().constants().earth.a.value,
                6378137.0, 1.0);

    // Both seed from the authentic WGS72 SGP4 — identical epoch state.
    std::cout << "\n=== shared authentic seed ===\n";
    double seed_diff = dist(A.epoch_state().position(), B.epoch_state().position());
    check("both modes share the authentic WGS72 seed", seed_diff == 0.0, seed_diff);

    // The modes diverge once propagated under their different models.
    std::cout << "\n=== modes diverge under propagation ===\n";
    dynamics::State<T> a20 = A.propagate(math::exact<T>(20));
    dynamics::State<T> b20 = B.propagate(math::exact<T>(20));
    double mode_diff = dist(a20.position(), b20.position());
    check("authentic vs boosted diverge (m)", mode_diff > 0.0, mode_diff);

    // from_tle(mode) dispatches by named parameter.
    std::cout << "\n=== from_tle(mode) dispatch ===\n";
    dynamics::DqSgp4Propagator<T> Fa = dynamics::DqSgp4Propagator<T>::from_tle(
        td, tol, dynamics::PropagatorMode::Authentic, dt);
    dynamics::DqSgp4Propagator<T> Fb = dynamics::DqSgp4Propagator<T>::from_tle(
        td, tol, dynamics::PropagatorMode::Boosted, dt);
    check_close("from_tle(Authentic) earth = WGS72",
                Fa.propagator().constants().earth.a.value, 6378135.0, 1.0);
    check_close("from_tle(Boosted) earth = WGS84",
                Fb.propagator().constants().earth.a.value, 6378137.0, 1.0);

    // SGP4 side: the authentic preset is WGS72.
    std::cout << "\n=== SGP4 authentic preset ===\n";
    geodesy::EquipotentialEllipsoid<T> sgp4_auth =
        sgp4::ModelSelector<T>::select("sgp4_standard", tol).ellipsoid;
    check_close("SGP4 sgp4_standard = WGS72 a (km)", sgp4_auth.a.value, 6378.135, 1e-3);

    std::cout << "\n" << (failed == 0 ? "PASS" : "FAIL") << " — " << failed
              << " failure(s)\n";
    return failed == 0 ? 0 : 1;
}
