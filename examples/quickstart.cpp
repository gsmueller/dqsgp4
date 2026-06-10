/// examples/quickstart.cpp — the guided tour of the public surface (gate EX2,
/// so this example can never rot). The docs help guide (docs/guide.html)
/// extracts its snippets VERBATIM from the `// [guide:tag]` regions below, so
/// what the guide shows is exactly what compiles and runs here.
///
/// Sections (mirroring the guide and the per-module usage snippets):
///   1. Propagate a TLE both ways — analytical SGP4 and the DQ facade
///   2. Perturbation presets — lunisolar / drag / SRP via DqForceOptions
///   3. Adaptive stepping — RKF7(8) with accept/reject diagnostics
///   4. Reading the three-error budget of any tracked quantity
///   5. Choosing T — double vs cpp_bin_float_50 (precision vs accuracy)
///   6. Assembling the low-level engine by hand (geopotential + RK4)
///   7. Constants with provenance + the level ellipsoid (normal gravity)
///   8. Sun and Moon positions (Meeus series → GCRS)
///   9. Kepler solvers and the generative equation of centre
///
/// Each section ends in PHYSICAL assertions (bands, not bits — the bit-grade
/// claims belong to the per-module gates: OR1, FM1, AD1, DS1).
///
/// Build (from the repo root, after vcvars64):
///   cl /std:c++20 /EHsc /O2 /I src /I vcpkg_installed\x64-windows\include ^
///      examples\quickstart.cpp src\tle\tle_parser.cpp /Fe:quickstart.exe

#include "dqsgp4.h"

// Opt-ins (deliberately NOT in the umbrella): the SGP4/SDP4 ORACLE surface
// (the test-against tier — included here exactly because this tour
// demonstrates it) and the generative Kepler series helpers.
#include "sgp4/model_selector.h"
#include "sgp4/sgp4_propagator.h"
#include "orbit/kepler_series.h"

#include <boost/multiprecision/cpp_bin_float.hpp>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using T = double;

namespace {
int bad = 0;
void check(const char* name, bool ok) {
    std::cout << (ok ? "  ok    " : "  FAIL  ") << name << "\n";
    if (!ok) ++bad;
}
}  // namespace

int main() {
    std::cout << "quickstart: the public surface in seven sections\n\n";
    std::cout << std::setprecision(9);

    // ---- 1. A TLE through both propagators --------------------------------
    std::cout << "[1] TLE -> the analytical model, then the numerical propagator\n";
    // [guide:sgp4]
    // The standard verification satellite 00005.
    std::string l1 = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
    std::string l2 = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
    tle::TleData td;
    if (!tle::parse(l1, l2, td)) { std::cerr << "TLE parse failed\n"; return 1; }

    // The analytical SGP4 propagator (explicit include; it is not part of the
    // umbrella header). A TLE's orbital elements are defined in terms of this
    // model, so it also provides the numerical propagator's initial state.
    sgp4::ModelConfiguration<T> config =
        sgp4::ModelSelector<T>::select("sgp4_standard", T(1e-12));
    tle::TleElements<T> elems = tle::TleElements<T>::from_tle_data(td);
    sgp4::Propagator<T> sgp4_prop(config, elems, T(1e-12));

    sgp4::StateVector<T> sv0 = sgp4_prop.propagate(math::exact<T>(0));  // minutes
    // [guide:end]
    std::cout << "  SGP4 t=0 pos (km): " << sv0.position_km.x.value << "  "
              << sv0.position_km.y.value << "  " << sv0.position_km.z.value << "\n";
    check("SGP4 t=0 matches the SGP4-VER reference (1e-3 km)",
          std::abs(sv0.position_km.x.value - 7022.465293) < 1e-3 &&
          std::abs(sv0.position_km.y.value - (-1400.082968)) < 1e-3 &&
          std::abs(sv0.position_km.z.value - 0.039952) < 1e-3);

    // [guide:facade]
    // The dual-quaternion numerical propagator, with the same call convention:
    // propagate(minutes since epoch). The epoch state is recovered through the
    // WGS72 SGP4 model at t = 0; the propagation itself is SE(3) integration.
    auto dq = dynamics::DqSgp4Propagator<T>::from_tle(
        td, T(1e-12), dynamics::PropagatorMode::Authentic, math::exact<T>(30));
    dynamics::State<T> s0 = dq.propagate(math::exact<T>(0));  // TEME metres
    // [guide:end]
    std::cout << "  DQ   t=0 pos (m):  " << s0.position().x.value << "  "
              << s0.position().y.value << "  " << s0.position().z.value << "\n";
    check("facade epoch state == SGP4 t=0 (km -> m, rel 1e-9)",
          std::abs(s0.position().x.value / (sv0.position_km.x.value * 1000.0) - 1.0) < 1e-9 &&
          std::abs(s0.position().y.value / (sv0.position_km.y.value * 1000.0) - 1.0) < 1e-9);

    // ---- 2. Perturbation presets ----------------------------------------
    std::cout << "\n[2] perturbation presets (opt-in; default stays bit-identical)\n";
    // [guide:presets]
    dynamics::DqForceOptions<T> opt;
    opt.lunisolar = true;                                       // Sun + Moon (TB1)
    opt.drag_B = math::TrackedValue<T>::measured("0.02", "0.002");  // C_d*A/m [m^2/kg] (ATM1)
    opt.srp_cr_area_over_mass =
        math::TrackedValue<T>::measured("0.02", "0.002");       // C_R*A/m [m^2/kg] (SRP1)

    auto boosted = dynamics::DqSgp4Propagator<T>::from_tle(
        td, T(1e-12), dynamics::PropagatorMode::Boosted, math::exact<T>(30), opt);
    dynamics::State<T> s60 = boosted.propagate(math::exact<T>(60));  // t = +60 min
    // [guide:end]
    {
        auto plain = dynamics::DqSgp4Propagator<T>::from_tle(
            td, T(1e-12), dynamics::PropagatorMode::Boosted, math::exact<T>(30));
        dynamics::State<T> p60 = plain.propagate(math::exact<T>(60));
        const double dx = s60.position().x.value - p60.position().x.value;
        const double dy = s60.position().y.value - p60.position().y.value;
        const double dz = s60.position().z.value - p60.position().z.value;
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        std::cout << "  presets move the 60-min endpoint by " << d << " m\n";
        check("perturbations land in the physical band (1..100 m @ 60 min)",
              d > 1.0 && d < 100.0);
    }

    // ---- 3. Adaptive stepping -------------------------------------------
    std::cout << "\n[3] adaptive stepping (RKF7(8), AD1-gated loop)\n";
    // [guide:adaptive]
    const dynamics::Propagator<T>& engine = boosted.propagator();
    math::TrackedValue<T> t_end = boosted.epoch_state().time + math::exact<T>(600);
    integrators::AdaptiveResult<T> ad = engine.propagate_adaptive(
        boosted.epoch_state(), t_end,
        math::exact<T>(60),     // initial dt [s]
        T(1e-6),                // per-step local-error tolerance
        T(1) / T(100));         // dt floor [s]
    // ad.state, ad.accepted_steps, ad.rejected_steps
    // [guide:end]
    std::cout << "  adaptive: " << ad.accepted_steps << " accepted, "
              << ad.rejected_steps << " rejected\n";
    check("adaptive lands on t_target (1e-9 s)",
          std::abs(ad.state.time.value - t_end.value) < 1e-9);
    {
        dynamics::State<T> fixed = engine.propagate_to(
            boosted.epoch_state(), t_end, math::exact<T>(60));
        const double dx = ad.state.position().x.value - fixed.position().x.value;
        const double dy = ad.state.position().y.value - fixed.position().y.value;
        const double dz = ad.state.position().z.value - fixed.position().z.value;
        check("adaptive ~ fixed-step over 600 s (sanity band 10 m)",
              std::sqrt(dx * dx + dy * dy + dz * dz) < 10.0);
    }

    // ---- 4. The three-error budget --------------------------------------
    std::cout << "\n[4] the three-error budget of any tracked quantity\n";
    // [guide:errors]
    math::TrackedValue<T> x = sv0.position_km.x;  // any computed quantity
    std::cout << "  value       " << x.value << " km\n"
              << "  measurement " << x.errors.measurement << " km (physical input sigma)\n"
              << "  precision   " << x.errors.precision << " km (representation/rounding)\n"
              << "  accuracy    " << x.errors.accuracy << " km (model truncation)\n"
              << "  total       " << x.errors.total() << " km (triangle-inequality bound)\n";
    // [guide:end]
    check("SGP4 budget channels are finite and non-negative",
          std::isfinite(x.errors.measurement) && x.errors.measurement >= 0 &&
          std::isfinite(x.errors.precision) && x.errors.precision >= 0 &&
          std::isfinite(x.errors.accuracy) && x.errors.accuracy >= 0);
    {
        // The numerically propagated state. The fixed-step integrators record
        // a conservative per-step truncation bound, assigned by component:
        // position-scale (h²/2)·A on the translation part of the pose,
        // velocity-scale h·A on the twist, and the rotational bounds — exactly
        // zero when no torque acts — on the rotational components
        // (runge_kutta_lie_group.md §6.2). Only the acceleration magnitude is
        // available to bound the higher derivatives, so the bound is loose by
        // design, and over a long fixed-step arc it accumulates past the state
        // scale and legitimately diverges through 1/r² and the drag
        // exponential. For a meaningful long-arc accuracy figure, use the
        // adaptive integrator, whose embedded 7(8) difference measures the
        // actual per-step error.
        dynamics::State<T> s1step = dq.propagate(math::exact<T>(1) / math::exact<T>(2));
        math::TrackedValue<T> p1 = s1step.position().x;
        std::cout << "  numerical, 1 step (30 s):  acc " << p1.errors.accuracy
                  << " m (the conservative per-step bound)\n";
        check("one-step accuracy bound is finite and at the expected scale (< 1e5 m)",
              std::isfinite(p1.errors.accuracy) && p1.errors.accuracy > 0 &&
              p1.errors.accuracy < 1e5);

        math::TrackedValue<T> px = s60.position().x;
        std::cout << "  numerical, 60 min:  meas " << px.errors.measurement
                  << " m, prec " << px.errors.precision
                  << " m, acc " << px.errors.accuracy
                  << " m (the fixed-step bound diverges on long arcs)\n";
        check("measurement and precision budgets stay finite through propagation",
              std::isfinite(px.errors.measurement) && px.errors.measurement > 0 &&
              std::isfinite(px.errors.precision) && px.errors.precision > 0);

        math::TrackedValue<T> ax = ad.state.position().x;
        std::cout << "  adaptive, 600 s:    acc " << ax.errors.accuracy
                  << " m (measured per-step error, accumulated)\n";
        check("adaptive accuracy stays finite over the arc",
              std::isfinite(ax.errors.accuracy) && ax.errors.accuracy >= 0);
    }

    // ---- 5. Choosing T ----------------------------------------------------
    std::cout << "\n[5] choosing T: double vs cpp_bin_float_50\n";
    // [guide:precision]
    using BF50 = boost::multiprecision::cpp_bin_float_50;
    // Everything is templated on the numeric type. Here the analytical model
    // runs at 50 digits; the test suite measures the same precision scaling
    // for the numerical propagator.
    sgp4::ModelConfiguration<BF50> config50 =
        sgp4::ModelSelector<BF50>::select("sgp4_standard", BF50(1e-12));
    sgp4::Propagator<BF50> prop50(
        config50, tle::TleElements<BF50>::from_tle_data(td), BF50(1e-12));
    sgp4::StateVector<BF50> sv50 = prop50.propagate(math::exact<BF50>(0));
    // The PRECISION channel tightens with the wider T; the ACCURACY channel
    // (model truncation) does not — wider arithmetic cannot improve a model.
    // [guide:end]
    {
        const double vd = sv0.position_km.x.value;
        const double v50 = static_cast<double>(sv50.position_km.x.value);
        const double pd = sv0.position_km.x.errors.precision;
        const double p50 = static_cast<double>(sv50.position_km.x.errors.precision);
        std::cout << "  x(double) = " << vd << " km, precision " << pd << " km\n";
        std::cout << "  x(bf50)   = " << v50 << " km, precision " << p50 << " km\n";
        check("values agree across T (1e-9 km)", std::abs(vd - v50) < 1e-9);
        check("bf50 precision is strictly tighter", p50 < pd);
    }

    // ---- 6. The engine by hand -------------------------------------------
    std::cout << "\n[7] assembling the low-level engine by hand\n";
    {
        // [guide:integrators]
        // The canonical gravity path: monopole + zonal J2..Jn in ONE
        // Cunningham pass (forces::geopotential, gate GEOPOT), an RK4
        // integrator, and the force-list engine the facade wraps.
        constants::ConstantsProvider<T> K =
            constants::ConstantsProvider<T>::wgs84(T(1e-12));
        constants::GravityField<T> field(
            constants::ZonalHarmonics<T>::egm2008(T(1e-12)),
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
        dynamics::Propagator<T> hand(K, dynamics::Inertia<T>::point_mass(math::exact<T>(1)),
                                     std::move(forces), integrator);
        dynamics::State<T> hs = hand.propagate_to(
            s0, s0.time + math::exact<T>(600), math::exact<T>(60));
        // [guide:end]
        const double r = std::sqrt(hs.position().x.value * hs.position().x.value +
                                   hs.position().y.value * hs.position().y.value +
                                   hs.position().z.value * hs.position().z.value);
        std::cout << "  |r(600 s)| = " << r / 1000.0 << " km\n";
        check("hand-assembled engine stays in the orbit's radial band",
              r > 6.9e6 && r < 1.1e7);
    }

    // ---- 7. Constants + the level ellipsoid --------------------------------
    std::cout << "\n[8] constants with provenance; Somigliana normal gravity\n";
    {
        // [guide:constants]
        // A ConstantsProvider is a derived level ellipsoid: four defining
        // parameters (provenance-tagged), everything else series-derived
        // with tracked budgets (gates C1, F3, SC1).
        constants::ConstantsProvider<T> K =
            constants::ConstantsProvider<T>::wgs84(T(1e-12));
        math::TrackedValue<T> GM = K.earth.GM;       // defining parameter
        math::TrackedValue<T> J2 = K.earth.J2n(1);   // derived even zonal
        std::cout << "  GM  = " << GM.value << " m^3/s^2 (total error "
                  << GM.errors.total() << ")\n"
                  << "  J2  = " << J2.value << " (series-derived, total error "
                  << J2.errors.total() << ")\n";
        // [guide:end]
        check("WGS-84 GM matches its defining value",
              std::abs(GM.value - 3.986004418e14) < 1e6);
        check("derived J2 lands on the published WGS-84 value (1e-7)",
              std::abs(J2.value - 1.08263e-3) < 1e-7);

        // [guide:geodesy]
        // Normal gravity on the ellipsoid (Somigliana), latitude in radians.
        math::TrackedValue<T> phi45 = math::pi<T>() / math::exact<T>(4);
        math::TrackedValue<T> g45 = K.earth.normal_gravity(phi45);
        math::TrackedValue<T> g0 = K.earth.normal_gravity(math::exact<T>(0));
        math::TrackedValue<T> g90 = K.earth.normal_gravity(math::pi<T>() / math::exact<T>(2));
        std::cout << "  gamma(0)  = " << g0.value << " m/s^2\n"
                  << "  gamma(45) = " << g45.value << " m/s^2\n"
                  << "  gamma(90) = " << g90.value << " m/s^2\n";
        // [guide:end]
        check("normal gravity at 45 deg in the textbook band",
              g45.value > 9.79 && g45.value < 9.82);
        check("gravity increases monotonically equator -> pole",
              g0.value < g45.value && g45.value < g90.value);
    }

    // ---- 8. Sun and Moon (Meeus series -> GCRS) ----------------------------
    std::cout << "\n[9] ephemeris: Sun and Moon GCRS positions at J2000\n";
    {
        // [guide:ephemeris]
        // Meeus generative series (DE430-gated, EPH/FRAME2): ecliptic-of-date
        // state -> GCRS Cartesian via the IAU2006 chain. Sun in AU, Moon in km.
        astronomy::Epoch<T> j2000 = astronomy::Epoch<T>::from_jd(
            math::exact<T>(2451545), astronomy::TimeScale::TT);
        math::Vector3<T> sun = ephemeris::body_position_gcrs<T>(
            ephemeris::sun_meeus_of_date<T>(j2000), j2000);
        math::Vector3<T> moon = ephemeris::body_position_gcrs<T>(
            ephemeris::moon_meeus_of_date<T>(j2000), j2000);
        // [guide:end]
        const double dsun = std::sqrt(sun.x.value * sun.x.value +
                                      sun.y.value * sun.y.value +
                                      sun.z.value * sun.z.value);
        const double dmoon = std::sqrt(moon.x.value * moon.x.value +
                                       moon.y.value * moon.y.value +
                                       moon.z.value * moon.z.value);
        std::cout << "  |sun|  = " << dsun << " AU\n"
                  << "  |moon| = " << dmoon << " km\n";
        check("Sun distance within the annual band (0.983..1.017 AU)",
              dsun > 0.983 && dsun < 1.017);
        check("Moon distance within the orbit band (356e3..407e3 km)",
              dmoon > 356000.0 && dmoon < 407000.0);
    }

    // ---- 9. Kepler solvers + the generative equation of centre ------------
    std::cout << "\n[10] orbit: Kepler residual and equation of centre\n";
    {
        // [guide:orbit]
        // The iterative solver (SGP4's modified-Kepler form, U = E + omega)
        // and the generative Fourier-Bessel series (converges for ALL e < 1,
        // dial-up K) — two routes to the same anomaly (gates C2, EPH).
        math::TrackedValue<T> e = math::TrackedValue<T>::from_truncated_decimal("0.1");
        math::TrackedValue<T> M = math::exact<T>(1);   // mean anomaly [rad]
        math::TrackedValue<T> E = orbit::eccentric_anomaly_series<T>(M, e, 8);
        math::TrackedValue<T> C = orbit::equation_of_center<T>(M, e, 8);
        // [guide:end]
        const double resid = std::abs(E.value - e.value * std::sin(E.value) - M.value);
        std::cout << "  E(M=1, e=0.1) = " << E.value
                  << "  Kepler residual = " << resid
                  << "  tracked accuracy = " << E.errors.accuracy << "\n"
                  << "  C = nu - M    = " << C.value << " rad\n";
        check("K=8 residual at its measured truncation grade (< 1e-8)", resid < 1e-8);
        check("tracked accuracy bound MAJORIZES the measured residual",
              E.errors.accuracy >= resid);
        check("equation of centre is positive past perigee for small e",
              C.value > 0.0 && C.value < 0.25);
    }

    std::cout << "\n" << (bad == 0 ? "PASS" : "FAIL") << " — " << bad
              << " failed check(s)\n";
    return bad == 0 ? 0 : 1;
}
