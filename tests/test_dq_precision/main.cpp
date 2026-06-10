/// test_dq_precision — PREC1 (vision review T3, design/DQSGP4_VISION_REVIEW.md §2).
///
/// THE SOLVER'S OWN ARBITRARY-PRECISION GATE — the DS1 analogue for the
/// product. DS1 proved the SGP4 oracle's deep-space propagation tightens its
/// precision budget at cpp_bin_float_50; until this gate, nothing measured the
/// same property for the dual-quaternion solver itself (the no-perceived-
/// fidelity rule: the flagship claim needs its measurement).
///
/// Propagates dynamics::Propagator<T> (geopotential J2..J4 + RK4) over a
/// 600 s LEO arc from an EXACT-integer Cartesian seed, and asserts:
///   1. the VALUES agree at the double-roundoff scale (bf50 is the truth
///      stick: the same arithmetic carried to 50 digits);
///   2. with the SERIES/ITERATION TOLERANCE SCALED WITH T (1e-12 at double,
///      1e-40 at bf50), the position and twist PRECISION channels tighten by
///      tens of orders — the "infinite precision upscale" property. The
///      tolerance is part of the dial: the wgs72 ellipsoid solves e² from J₂
///      iteratively TO THE CALLER'S TOLERANCE, and that convergence residual
///      is honestly deposited as precision on the regenerated even zonals.
///   3. the floor is REAL and honest (measured here, pinned as a check):
///      at bf50 with the double-grade 1e-12 tolerance, the propagated
///      precision does NOT tighten — the J₂ convergence floor (~2e-11)
///      dominates identically at both T. Upscaling T without tightening the
///      tolerance buys nothing, by design, and the budget SAYS so.
///
/// UNPEEL DEMONSTRATION (gate ARCH1's companion): this translation unit
/// includes NO sgp4/ and NO tle/ header — the solver compiles, seeds, and
/// propagates without the oracle tier existing in the TU at all. The seed is
/// a plain Cartesian state (exact integers), not a TLE.
///
/// ExeGate PREC1: returns a nonzero exit code on any failed check.

#include "constants/constants_provider.h"
#include "constants/gravity_field.h"
#include "constants/tesseral_harmonics.h"
#include "constants/zonal_harmonics.h"
#include "dynamics/inertia.h"
#include "dynamics/propagator.h"
#include "dynamics/state.h"
#include "forces/geopotential.h"
#include "integrators/runge_kutta.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

using BF50 = boost::multiprecision::cpp_bin_float_50;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "") {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name
                                  << (detail.empty() ? "" : "  [" + detail + "]") << "\n"; }
}

/// Propagate the standard arc at type T with the caller's series/iteration
/// tolerance; report the final position value, its precision budget, and the
/// twist precision (all as double for comparison; the bf50 precision
/// exponents are far inside double range).
template<typename T>
void run_arc(const T& tol, double pos_out[3], double pos_prec_out[3],
             double& vel_prec_out) {
    using TV = math::TrackedValue<T>;

    constants::ConstantsProvider<T> K =
        constants::ConstantsProvider<T>::wgs72(tol);
    constants::GravityField<T> field(
        constants::ZonalHarmonics<T>::wgs72(tol),
        constants::TesseralHarmonics<T>{});

    std::vector<dynamics::ForceFn<T>> forces;
    forces.push_back([field](const dynamics::State<T>& s,
                             const constants::ConstantsProvider<T>& KK) {
        return forces::geopotential(s, KK, field, 4, 0, math::exact<T>(0));
    });
    dynamics::IntegratorFn<T> integrator =
        [](const dynamics::State<T>& y0, const math::TrackedValue<T>& dt,
           const integrators::AccelFn<T>& accel) {
            return integrators::runge_kutta_4(y0, dt, accel);
        };
    dynamics::Propagator<T> prop(K, dynamics::Inertia<T>::point_mass(math::exact<T>(1)),
                                 std::move(forces), integrator);

    // Exact-integer LEO seed (~7000 km circular-ish): the inputs carry ZERO
    // measurement error, so the output budgets are pure computation.
    math::Vector3<T> r0(math::exact<T>(7000000), math::exact<T>(0), math::exact<T>(0));
    math::Vector3<T> v0(math::exact<T>(0), math::exact<T>(7546), math::exact<T>(0));
    dynamics::State<T> s0 = dynamics::State<T>::from_kinematics(
        math::Quaternion<T>::identity(), r0, math::Vector3<T>(), v0, math::exact<T>(0));

    dynamics::State<T> sf = prop.propagate_to(s0, math::exact<T>(600), math::exact<T>(60));

    math::Vector3<T> rf = sf.position();
    pos_out[0] = static_cast<double>(rf.x.value);
    pos_out[1] = static_cast<double>(rf.y.value);
    pos_out[2] = static_cast<double>(rf.z.value);
    pos_prec_out[0] = static_cast<double>(rf.x.errors.precision);
    pos_prec_out[1] = static_cast<double>(rf.y.errors.precision);
    pos_prec_out[2] = static_cast<double>(rf.z.errors.precision);
    vel_prec_out = static_cast<double>(sf.twist.linear.y.errors.precision);
}

}  // namespace

int main() {
    std::cout << "test_dq_precision (PREC1): the solver's arbitrary-precision "
                 "upscale, measured\n";

    double pos_d[3], prec_d[3], vprec_d;
    double pos_h[3], prec_h[3], vprec_h;
    double pos_f[3], prec_f[3], vprec_f;  // bf50 at the DOUBLE-grade tolerance
    run_arc<double>(1e-12, pos_d, prec_d, vprec_d);
    run_arc<BF50>(BF50("1e-40"), pos_h, prec_h, vprec_h);
    run_arc<BF50>(BF50("1e-12"), pos_f, prec_f, vprec_f);

    std::cout << "  double: r=(" << pos_d[0] << ", " << pos_d[1] << ", " << pos_d[2]
              << ") m, |prec| x/y = " << prec_d[0] << " / " << prec_d[1] << " m\n";
    std::cout << "  bf50:   r=(" << pos_h[0] << ", " << pos_h[1] << ", " << pos_h[2]
              << ") m, |prec| x/y = " << prec_h[0] << " / " << prec_h[1] << " m\n";

    // 1. Values agree at the double-roundoff scale. ~10 RK4 steps of ~1e7 m
    //    arithmetic at 1e-16 relative roundoff -> sub-mm divergence; band 1 mm.
    const double dx = pos_d[0] - pos_h[0];
    const double dy = pos_d[1] - pos_h[1];
    const double dz = pos_d[2] - pos_h[2];
    const double dv = std::sqrt(dx * dx + dy * dy + dz * dz);
    check("values agree across T (double vs bf50 < 1 mm over 600 s)",
          dv < 1e-3, "dv=" + std::to_string(dv));

    // 2. Sanity: the arc actually went somewhere (not a degenerate run).
    check("arc propagated (moved > 4000 km along-track)",
          std::abs(pos_h[1]) > 4.0e6);

    // 3. The precision channel is REAL at double (positive, finite, sane):
    //    representation-scale per-step roundoff accumulated over the arc.
    check("double precision budget is positive and finite",
          prec_d[0] > 0.0 && prec_d[1] > 0.0 && std::isfinite(prec_d[0]) &&
          std::isfinite(prec_d[1]));

    // 4. THE UPSCALE: with (T, tolerance) dialed together, the channel
    //    tightens by tens of orders (assert >= 1e20; measured ~1e25+).
    check("bf50@1e-40 position precision tightens by >= 1e20 (x)",
          prec_h[0] > 0.0 && prec_d[0] / prec_h[0] > 1e20,
          "ratio=" + std::to_string(prec_d[0] / (prec_h[0] > 0 ? prec_h[0] : 1e-300)));
    check("bf50@1e-40 position precision tightens by >= 1e20 (y)",
          prec_h[1] > 0.0 && prec_d[1] / prec_h[1] > 1e20);
    check("bf50@1e-40 twist precision tightens by >= 1e20",
          vprec_h > 0.0 && vprec_d / vprec_h > 1e20);

    // 5. THE FLOOR, pinned: bf50 at the double-grade tolerance does NOT
    //    tighten — the J2 e²-iteration residual (deposited as precision on
    //    the regenerated even zonals) dominates identically at both T.
    //    Honest behavior: T alone is not the dial; (T, tolerance) is.
    check("bf50@1e-12 stays at the tolerance floor (ratio < 10)",
          prec_f[0] > 0.0 && prec_d[0] / prec_f[0] < 10.0,
          "ratio=" + std::to_string(prec_d[0] / (prec_f[0] > 0 ? prec_f[0] : 1e-300)));

    std::cout << (failed == 0 ? "ALL PASS" : "FAILURES") << " — " << passed
              << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
