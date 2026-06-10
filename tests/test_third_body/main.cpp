/// test_third_body — Cartesian third-body perturbation force (gate TB1, L4).
///
/// Validates the Newtonian third-body acceleration (src/forces/third_body.h),
/// theory design/derivations/third_body_perturbation.md, in two independent
/// layers (no perceived fidelity):
///
///   Layer 1 — the FORMULA, machine-eps vs an independent re-derivation. The
///     Battin stable form a = −(μ/‖s‖³)[r(1+q)^(−3/2) + f(q)s] vs the naive
///     a = μ[(s−r)/‖s−r‖³ − s/‖s‖³] at cpp_bin_float_50 must agree to ≪ the
///     perturbation (they are the same quantity). At `double` the naive form
///     loses ~5 figures to catastrophic cancellation while Battin holds to ~1e-16
///     — the cancellation Battin exists to remove (and the perceived-fidelity trap).
///
///   Layer 2 — the EPHEMERIS→acceleration map vs JPL DE430. Feed the force the
///     MODEL ephemeris (sun_meeus / moon_meeus) and compare to the acceleration
///     from the in-repo DE430 body positions (independent numerical integration).
///     The residual is the body-position truncation, which the force's TRACKED
///     accuracy must majorize, for Sun and Moon, over five epochs spanning decades.
///
///   plus physical sanity (|a☾| ≈ 2.2|a☉|, |a| ~ 1e-7 of central gravity at LEO,
///   acceleration toward the body for an on-line satellite) and precision scaling.
///
/// Oracle: tools/gen_third_body_oracle.py (DE430 Sun/Moon GCRS km). Exit 0 iff all pass.

#include "astronomy/epoch.h"
#include "dynamics/inertia.h"
#include "dynamics/propagator.h"
#include "forces/gravity_central.h"
#include "forces/third_body.h"
#include "integrators/runge_kutta.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

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

template<typename T>
math::Vector3<T> vec3(double x, double y, double z) {
    return math::Vector3<T>(tv<T>(x), tv<T>(y), tv<T>(z));
}

/// Naive §2 form — the independent oracle for the FORMULA (exact algebra; suffers
/// the LEO cancellation at low precision, which is exactly what Battin removes).
template<typename T>
math::Vector3<T> naive_third_body(const math::Vector3<T>& r, const math::Vector3<T>& s,
                                  const math::TrackedValue<T>& mu) {
    using TV = math::TrackedValue<T>;
    math::Vector3<T> d(s.x - r.x, s.y - r.y, s.z - r.z);
    TV d2 = d.x * d.x + d.y * d.y + d.z * d.z;
    TV d3 = d2 * sqrt(d2);
    TV s2 = s.x * s.x + s.y * s.y + s.z * s.z;
    TV s3 = s2 * sqrt(s2);
    return math::Vector3<T>(mu * (d.x / d3 - s.x / s3),
                            mu * (d.y / d3 - s.y / s3),
                            mu * (d.z / d3 - s.z / s3));
}

template<typename T>
double dmag(const math::Vector3<T>& a) {
    using std::sqrt;
    return static_cast<double>(sqrt(a.x.value * a.x.value + a.y.value * a.y.value
                                  + a.z.value * a.z.value));
}

// Distance between two acceleration vectors (as doubles), allowing mixed T for one.
template<typename T>
double dist_to_ref(const math::Vector3<T>& a, const double ref[3]) {
    double dx = static_cast<double>(a.x.value) - ref[0];
    double dy = static_cast<double>(a.y.value) - ref[1];
    double dz = static_cast<double>(a.z.value) - ref[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

template<typename T>
double accuracy_l1(const math::Vector3<T>& a) {
    using std::abs;
    return static_cast<double>(a.x.errors.accuracy + a.y.errors.accuracy + a.z.errors.accuracy);
}

// DE430 Sun/Moon geocentric positions (J2000 mean-equatorial ≈ GCRS, km) — emitted
// by tools/gen_third_body_oracle.py from the in-repo JPL DE430 table.
struct Row { double jd; double sun[3]; double moon[3]; };
const Row ORACLE[] = {
    {2451545.0, {26499033.7568, -132757417.3548, -57556718.3992}, {-291608.3848, -266716.8335, -76102.4860}},
    {2455455.0, {-149117802.8510, 18254421.0806, 7914251.3915}, {-31395.8604, -354133.4293, -162301.4360}},
    {2458849.5, {24884971.6986, -133017487.8671, -57663412.0880}, {390185.6383, -76522.5985, -70724.6561}},
    {2460462.5, {50222120.0327, 131332676.0738, 56930496.0349}, {368305.3099, -10570.1034, -14560.1512}},
    {2449796.5, {148930123.9145, -2426286.9030, -1051657.3464}, {-289619.4202, -207427.2619, -88753.8919}},
};
const int N_ORACLE = 5;

}  // namespace

int main() {
    std::cout << std::setprecision(6);
    std::cout << "test_third_body: Cartesian third-body perturbation (TB1, L4)\n\n";
    using T = double;

    // ============================================================
    // Layer 1 — the formula: Battin vs naive (machine-eps) + cancellation demo
    // ============================================================
    std::cout << "=== Layer 1: formula vs naive re-derivation ===\n";
    {
        // A Sun-scale geometry: LEO satellite, body ~1 AU away (ratio ~2e4 ⇒ the
        // naive difference cancels ~5 figures).
        const double mu_d = 1.32712440018e20;
        // bf50 truth.
        math::Vector3<cpp_bin_float_50> r_b = vec3<cpp_bin_float_50>(7.0e6, 1.0e6, 2.0e6);
        math::Vector3<cpp_bin_float_50> s_b =
            vec3<cpp_bin_float_50>(2.65e10, -1.33e11, -5.76e10);
        math::TrackedValue<cpp_bin_float_50> mu_b = tv<cpp_bin_float_50>(mu_d);
        math::Vector3<cpp_bin_float_50> battin_b = forces::third_body_accel(r_b, s_b, mu_b);
        math::Vector3<cpp_bin_float_50> naive_b = naive_third_body(r_b, s_b, mu_b);
        double ref[3] = {static_cast<double>(battin_b.x.value),
                         static_cast<double>(battin_b.y.value),
                         static_cast<double>(battin_b.z.value)};
        double amag = dmag(battin_b);
        double algebraic = dist_to_ref(naive_b, ref) / amag;
        std::cout << "  |a| = " << amag << " m/s^2;  |Battin-naive|/|a| @bf50 = "
                  << algebraic << "\n";
        check("Battin == naive at bf50 (same quantity, <1e-30)", algebraic < 1e-30);

        // double: Battin accurate, naive degraded by cancellation.
        math::Vector3<T> r_d = vec3<T>(7.0e6, 1.0e6, 2.0e6);
        math::Vector3<T> s_d = vec3<T>(2.65e10, -1.33e11, -5.76e10);
        math::TrackedValue<T> mu_dd = tv<T>(mu_d);
        double battin_err = dist_to_ref(forces::third_body_accel(r_d, s_d, mu_dd), ref) / amag;
        double naive_err = dist_to_ref(naive_third_body(r_d, s_d, mu_dd), ref) / amag;
        std::cout << "  @double rel err:  Battin = " << battin_err
                  << "   naive = " << naive_err << "\n";
        check("Battin is accurate at double (<1e-12)", battin_err < 1e-12);
        check("naive loses figures to cancellation (>=100x Battin)",
              naive_err > 100.0 * battin_err);
    }

    // ============================================================
    // Layer 2 — ephemeris→acceleration vs JPL DE430 (Sun and Moon)
    // ============================================================
    std::cout << "\n=== Layer 2: model ephemeris vs DE430 acceleration ===\n";
    forces::ThirdBody<T> sun = forces::sun_third_body<T>();
    forces::ThirdBody<T> moon = forces::moon_third_body<T>();
    // A fixed GCRS test point (LEO, 7000 km geocentric) — the perturbation is a
    // field, evaluated at this point; not an orbit.
    math::Vector3<T> r_sat = vec3<T>(7.0e6, 0.0, 0.0);

    bool sun_ok = true, moon_ok = true;
    double worst_sun = 0.0, worst_moon = 0.0;
    for (int i = 0; i < N_ORACLE; ++i) {
        astronomy::Epoch<T> ep =
            astronomy::Epoch<T>::from_jd(tv<T>(ORACLE[i].jd), astronomy::TimeScale::TT);

        // Sun
        math::Vector3<T> a_model_s = forces::third_body_perturbation(r_sat, sun, ep);
        math::Vector3<T> s_de = vec3<T>(ORACLE[i].sun[0] * 1000.0, ORACLE[i].sun[1] * 1000.0,
                                        ORACLE[i].sun[2] * 1000.0);
        math::Vector3<T> a_de_s = forces::third_body_accel(r_sat, s_de, sun.mu);
        double ref_s[3] = {a_de_s.x.value, a_de_s.y.value, a_de_s.z.value};
        double res_s = dist_to_ref(a_model_s, ref_s);
        double acc_s = accuracy_l1(a_model_s);
        if (res_s > acc_s) sun_ok = false;
        worst_sun = std::max(worst_sun, res_s / dmag(a_de_s));

        // Moon
        math::Vector3<T> a_model_m = forces::third_body_perturbation(r_sat, moon, ep);
        math::Vector3<T> m_de = vec3<T>(ORACLE[i].moon[0] * 1000.0, ORACLE[i].moon[1] * 1000.0,
                                        ORACLE[i].moon[2] * 1000.0);
        math::Vector3<T> a_de_m = forces::third_body_accel(r_sat, m_de, moon.mu);
        double ref_m[3] = {a_de_m.x.value, a_de_m.y.value, a_de_m.z.value};
        double res_m = dist_to_ref(a_model_m, ref_m);
        double acc_m = accuracy_l1(a_model_m);
        if (res_m > acc_m) moon_ok = false;
        worst_moon = std::max(worst_moon, res_m / dmag(a_de_m));

        std::cout << "  jd " << std::setprecision(8) << ORACLE[i].jd << std::setprecision(4)
                  << "  Sun: res/|a|=" << res_s / dmag(a_de_s) << " (acc " << acc_s << " vs res " << res_s
                  << ")  Moon: res/|a|=" << res_m / dmag(a_de_m) << "\n";
    }
    std::cout << "  worst Sun rel residual " << worst_sun
              << ", worst Moon rel residual " << worst_moon << "\n";
    check("Sun: tracked accuracy majorizes DE430 residual (all epochs)", sun_ok);
    check("Moon: tracked accuracy majorizes DE430 residual (all epochs)", moon_ok);
    check("Sun residual within Meeus grade (<1e-3 rel)", worst_sun < 1e-3);
    check("Moon residual within Meeus grade (<1e-3 rel)", worst_moon < 1e-3);

    // ============================================================
    // Physical sanity
    // ============================================================
    std::cout << "\n=== physical sanity ===\n";
    {
        astronomy::Epoch<T> ep =
            astronomy::Epoch<T>::from_jd(tv<T>(ORACLE[0].jd), astronomy::TimeScale::TT);
        double a_sun = dmag(forces::third_body_perturbation(r_sat, sun, ep));
        double a_moon = dmag(forces::third_body_perturbation(r_sat, moon, ep));
        double ratio = a_moon / a_sun;
        std::cout << "  |a_sun| = " << a_sun << "  |a_moon| = " << a_moon
                  << "  ratio = " << ratio << "\n";
        // The Moon's tidal pull exceeds the Sun's by a factor of order 2 despite
        // being ~2.7e7x less massive — proximity (1/r^3) beats mass. The exact
        // ratio is geometry/epoch-dependent (textbook average ~2.2; here, at J2000
        // with the Sun near perihelion and distinct satellite-body angles, 2.72).
        check("|a_moon| > |a_sun| by a factor of order 2 (proximity beats mass)",
              ratio > 1.5 && ratio < 4.0);
        check("|a_3body| ~ 1e-6 m/s^2 at LEO (1e-7 of central g)",
              a_sun > 1e-8 && a_moon < 1e-4);

        // Sign: a satellite ON the Earth–body line, sunward of Earth, is pulled
        // toward the body (tidal stretch). Use the DE430 Sun direction.
        math::Vector3<T> s_de = vec3<T>(ORACLE[0].sun[0] * 1000.0, ORACLE[0].sun[1] * 1000.0,
                                        ORACLE[0].sun[2] * 1000.0);
        double smag = dmag(s_de);
        math::Vector3<T> r_on = vec3<T>(s_de.x.value / smag * 7.0e6, s_de.y.value / smag * 7.0e6,
                                        s_de.z.value / smag * 7.0e6);
        math::Vector3<T> a_on = forces::third_body_accel(r_on, s_de, sun.mu);
        double along = (a_on.x.value * s_de.x.value + a_on.y.value * s_de.y.value
                      + a_on.z.value * s_de.z.value) / smag;
        std::cout << "  on-line satellite: a along body direction = " << along << " m/s^2\n";
        check("on-line satellite pulled toward the body (a·ŝ > 0)", along > 0.0);
    }

    // ============================================================
    // Precision scaling
    // ============================================================
    std::cout << "\n=== precision scaling ===\n";
    {
        auto prec = [](auto tag) -> double {
            using U = decltype(tag);
            math::Vector3<U> r = vec3<U>(7.0e6, 1.0e6, 2.0e6);
            math::Vector3<U> s = vec3<U>(2.65e10, -1.33e11, -5.76e10);
            math::TrackedValue<U> mu = tv<U>(1.32712440018e20);
            math::Vector3<U> a = forces::third_body_accel(r, s, mu);
            return static_cast<double>(a.x.errors.precision);
        };
        double pd = prec(double{});
        double pb = prec(cpp_bin_float_50{});
        std::cout << "  precision: double = " << pd << "  bf50 = " << pb << "\n";
        check("precision > 0 (framework alive)", pd > 0.0 && pb > 0.0);
        check("precision tightens with wider T", pb < pd);
    }

    // ============================================================
    // Propagator wiring (third_body_perturbation.md §9)
    // ============================================================
    std::cout << "\n=== propagator wiring: State->Wrench + epoch map + precession frame ===\n";
    {
        constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
        astronomy::Epoch<T> base =
            astronomy::Epoch<T>::from_jd(tv<T>(2460462.5), astronomy::TimeScale::TT);  // 2024
        math::Vector3<T> rsat = vec3<T>(7.0e6, 0.0, 0.0);
        math::Vector3<T> zero;
        dynamics::State<T> s0 = dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), rsat, zero, zero, tv<T>(0.0));
        dynamics::Wrench<T> w = forces::third_body_force(s0, K, sun, base);

        // (a) wiring composition: identity attitude -> wrench.force == manual accel.
        math::Vector3<T> s_gcrs = sun.position(base);
        math::TrackedValue<T> tcen = astronomy::centuries_since_j2000(base);
        math::Vector3<T> s_teme = astronomy::precession_iau2006<T>(tcen) * s_gcrs;
        math::Vector3<T> a_manual = forces::third_body_accel(rsat, s_teme, sun.mu);
        double ref_a[3] = {a_manual.x.value, a_manual.y.value, a_manual.z.value};
        double dconsist = dist_to_ref(w.force, ref_a) / dmag(a_manual);
        std::cout << "  wrench vs manual rel: " << dconsist << "\n";
        check("wrench force == third_body_accel(r, P*s_gcrs) (wiring composes)", dconsist < 1e-12);

        // (b) epoch map: +1 day moves the body -> the force changes.
        dynamics::State<T> s1 = dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), rsat, zero, zero, tv<T>(86400.0));
        dynamics::Wrench<T> w1 = forces::third_body_force(s1, K, sun, base);
        double ref_w[3] = {w.force.x.value, w.force.y.value, w.force.z.value};
        double dday = dist_to_ref(w1.force, ref_w) / dmag(w.force);
        std::cout << "  force change over +1 day: " << dday << " rel\n";
        check("epoch map advances the ephemeris (force changes over a day)", dday > 1e-3);

        // (c) precession frame applied: s_teme rotated from s_gcrs by ~0.3 deg at 2024.
        double sg = dmag(s_gcrs), st = dmag(s_teme);
        double cosang = (s_gcrs.x.value * s_teme.x.value + s_gcrs.y.value * s_teme.y.value
                       + s_gcrs.z.value * s_teme.z.value) / (sg * st);
        double ang = std::acos(std::min(1.0, cosang));
        std::cout << "  precession angle GCRS->TEME at 2024: "
                  << ang * 648000.0 / 3.14159265358979 << " arcsec\n";
        check("precession frame applied (angle ~0.3 deg, not identity)", ang > 1e-3 && ang < 0.02);

        // (d) omitted nutation/Eqeq deposited as a frame accuracy bound.
        double amag = dmag(w.force);
        double acc = static_cast<double>(w.force.x.errors.accuracy);
        std::cout << "  force.x accuracy " << acc << " (frame bound ~" << 1.45e-4 * amag << ")\n";
        check("omitted nutation/Eqeq tracked as a frame accuracy bound", acc > 1e-5 * amag);

        // (e) usable in the DQ propagator: central + third-body stays finite and perturbs.
        dynamics::Inertia<T> inertia = dynamics::Inertia<T>::point_mass(tv<T>(1.0));
        double vc = std::sqrt(K.earth.GM.value / 7.0e6);
        math::Vector3<T> v0 = vec3<T>(0.0, vc, 0.0);
        dynamics::State<T> orb0 = dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), rsat, zero, v0, tv<T>(0.0));
        dynamics::IntegratorFn<T> integ =
            [](const dynamics::State<T>& y, const math::TrackedValue<T>& dt,
               const integrators::AccelFn<T>& f) { return integrators::runge_kutta_4(y, dt, f); };
        auto central = [](const dynamics::State<T>& s, const constants::ConstantsProvider<T>& KK) {
            return forces::gravity_central(s, KK);
        };
        std::vector<dynamics::ForceFn<T>> fc, ft;
        fc.push_back(central);
        ft.push_back(central);
        ft.push_back(forces::make_third_body_force<T>(sun, base));
        dynamics::Propagator<T> pc(K, inertia, std::move(fc), integ);
        dynamics::Propagator<T> pt(K, inertia, std::move(ft), integ);
        dynamics::State<T> ec = pc.propagate_to(orb0, tv<T>(600.0), tv<T>(60.0));
        dynamics::State<T> et = pt.propagate_to(orb0, tv<T>(600.0), tv<T>(60.0));
        double pos_diff = std::sqrt(
            std::pow(et.position().x.value - ec.position().x.value, 2) +
            std::pow(et.position().y.value - ec.position().y.value, 2) +
            std::pow(et.position().z.value - ec.position().z.value, 2));
        bool finite = std::isfinite(et.position().x.value) && std::isfinite(et.position().y.value);
        std::cout << "  central-vs-(central+3body) position diff over 600 s: " << pos_diff << " m\n";
        check("third-body force usable in the DQ propagator (finite + perturbs)",
              finite && pos_diff > 0.0);
    }

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
