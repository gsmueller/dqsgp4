/// SGP4 Verification Test Driver
///
/// Parses the SGP4-VER.TLE test file and (when propagator is complete)
/// compares output against the tcppver.out reference.
///
/// For now: parses TLEs, converts to TrackedValue elements, and verifies
/// the astronomy and perturbation modules against the SGP4 magic numbers.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include <boost/multiprecision/cpp_bin_float.hpp>

#include "math/tracked_value.h"
#include "math/angles.h"
#include "geodesy/equipotential_ellipsoid.h"
#include "astronomy/solar_system.h"
#include "perturbation/brouwer.h"
#include "perturbation/kaula.h"
#include "sgp4/model_selector.h"
#include "sgp4/sgp4_propagator.h"
#include "diagnostic.h"
#include "tcppver_parser.h"
#include "tle/tle_parser.h"

using T = double; // Use double for validation against legacy output

/// Parse the SGP4-VER.TLE file which has a special format:
/// - Comment lines start with #
/// - TLE line 1 and line 2 in standard format
/// - Line 2 has appended start/stop/step times after column 69
struct TestCase {
    tle::TleData tle;
    double start_min;
    double stop_min;
    double step_min;
};

std::vector<TestCase> parse_ver_file(const std::string& path) {
    std::vector<TestCase> cases;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open " << path << "\n";
        return cases;
    }

    std::string line;
    std::string pending_line1;
    bool have_line1 = false;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        if (line[0] == '1' && line.size() >= 69) {
            pending_line1 = line.substr(0, 69);
            have_line1 = true;
        }
        else if (line[0] == '2' && line.size() >= 69 && have_line1) {
            TestCase tc;
            std::string line2_tle = line.substr(0, 69);

            if (tle::parse(pending_line1, line2_tle, tc.tle)) {
                // Parse the appended start/stop/step from after column 69
                if (line.size() > 69) {
                    std::istringstream extra(line.substr(69));
                    extra >> tc.start_min >> tc.stop_min >> tc.step_min;
                } else {
                    tc.start_min = 0;
                    tc.stop_min = 1440;
                    tc.step_min = 120;
                }
                cases.push_back(tc);
            }
            have_line1 = false;
        }
    }

    return cases;
}

/// Compare a computed value against a legacy SGP4 constant.
/// Reports the discrepancy and whether it's within expected bounds.
void check_magic_number(const char* name,
                        double computed, double legacy,
                        double expected_tolerance)
{
    double diff = std::abs(computed - legacy);
    double rel = (legacy != 0.0) ? diff / std::abs(legacy) : diff;
    bool ok = diff <= expected_tolerance;

    std::cout << (ok ? "  PASS" : "  FAIL") << " " << name
              << ": computed=" << std::setprecision(12) << computed
              << " legacy=" << legacy
              << " diff=" << std::scientific << diff
              << " rel=" << rel
              << (ok ? "" : " *** EXCEEDS TOLERANCE ***")
              << "\n";
}

int main() {
    std::cout << "=== SGP4 Verification Test Suite ===\n\n";

    // --- Magic Number Tests: Astronomy Module ---
    std::cout << "--- Astronomy Module: Solar/Lunar Constants ---\n";

    astronomy::FundamentalConstants<T> astro_fc = astronomy::FundamentalConstants<T>::sgp4_standard();
    astronomy::DerivedOrbitalElements<T> astro = astronomy::DerivedOrbitalElements<T>::compute(astro_fc);

    // Compare against SGP4 deep-space constants (from code_to_documentation_xref.md)
    // ZNS = 1.19459E-5 (solar mean motion in rad/min)
    check_magic_number("ZNS (solar mean motion)",
        astro.solar_mean_motion.value, 1.19459e-5, 1e-9);

    // ZES = 0.01675 (solar eccentricity)
    check_magic_number("ZES (solar eccentricity)",
        astro.solar_eccentricity.value, 0.01675, 1e-10);

    // ZSINIS = 0.39785416 (sin obliquity)
    check_magic_number("ZSINIS (sin obliquity)",
        astro.sin_obliquity.value, 0.39785416, 1e-6);

    // ZCOSIS = 0.91744867 (cos obliquity)
    check_magic_number("ZCOSIS (cos obliquity)",
        astro.cos_obliquity.value, 0.91744867, 1e-6);

    // ZNL = 1.5835218E-4 (lunar mean motion in rad/min)
    check_magic_number("ZNL (lunar mean motion)",
        astro.lunar_mean_motion.value, 1.5835218e-4, 1e-9);

    // ZEL = 0.05490 (lunar eccentricity)
    check_magic_number("ZEL (lunar eccentricity)",
        astro.lunar_eccentricity.value, 0.05490, 1e-10);

    // --- Magic Number Tests: Geodesy Module ---
    std::cout << "\n--- Geodesy Module: WGS72 Constants (J2 init path) ---\n";

    // Build WGS72 ellipsoid using J₂ as the defining parameter
    // (the SGP4 standard gravity model specifies J₂ directly, not 1/f)
    math::TrackedValue<T> a_wgs72 = math::TrackedValue<T>::defined("6378.135");   // km
    math::TrackedValue<T> J2_wgs72 = math::TrackedValue<T>::defined("0.001082616");
    math::TrackedValue<T> GM_wgs72 = math::TrackedValue<T>::defined("398600.8");  // km³/s²
    math::TrackedValue<T> omega_wgs72 = math::TrackedValue<T>::defined("7.292115e-5");

    geodesy::EquipotentialEllipsoid<T> wgs72 = geodesy::EquipotentialEllipsoid<T>::from_J2(
        a_wgs72, J2_wgs72, GM_wgs72, omega_wgs72, 1e-15);

    // J2 (WGS72): should match 0.001082616
    check_magic_number("J2 (WGS72)",
        wgs72.J2n(1).value, 0.001082616, 1e-8);

    // CK2 = J2/2
    double ck2_computed = wgs72.J2n(1).value / 2.0;
    check_magic_number("CK2 = J2/2",
        ck2_computed, 0.5 * 0.001082616, 1e-12);

    // xke = 60 / sqrt(re³/mu)
    double re = 6378.135;
    double mu = 398600.8;
    double xke_legacy = 60.0 / std::sqrt(re*re*re / mu);
    double xke_computed = 60.0 / std::sqrt(
        a_wgs72.value * a_wgs72.value * a_wgs72.value / GM_wgs72.value);
    check_magic_number("xke = 60/sqrt(re³/mu)",
        xke_computed, xke_legacy, 1e-12);

    // --- Parse SGP4-VER.TLE ---
    std::cout << "\n--- Parsing SGP4-VER.TLE ---\n";

    std::string ver_path = "sgp4_references/aholinch_sgp4/data/SGP4-VER.TLE";
    std::vector<TestCase> test_cases = parse_ver_file(ver_path);
    std::cout << "Parsed " << test_cases.size() << " test cases\n\n";

    for (size_t i = 0; i < test_cases.size() && i < 5; ++i) {
        TestCase& tc = test_cases[i];
        std::cout << "Test case " << (i+1) << ": sat " << tc.tle.satellite_id
                  << " epoch=" << tc.tle.epoch_year << "/" << std::setprecision(8)
                  << tc.tle.epoch_day
                  << " e=" << std::setprecision(7) << tc.tle.eccentricity
                  << " i=" << tc.tle.inclination_deg << "°"
                  << " n=" << std::setprecision(11) << tc.tle.mean_motion_rev_day << " rev/day"
                  << " t=[" << tc.start_min << "," << tc.stop_min << "]"
                  << " step=" << tc.step_min << " min\n";

        // Convert to TrackedValue elements
        tle::TleElements<T> elements = tle::TleElements<T>::from_tle_data(tc.tle);

        // Compute Brouwer secular rates for this orbit
        math::TrackedValue<T> cos_i = cos(elements.inclination);
        math::TrackedValue<T> e2 = elements.eccentricity * elements.eccentricity;
        math::TrackedValue<T> J2 = wgs72.J2n(1);
        math::TrackedValue<T> J4 = wgs72.J2n(2);

        // Need semi-major axis from mean motion
        // a = (xke / n)^(2/3) — but we need n in rad/min
        double n_rad_min = elements.mean_motion.value;
        double a_er = std::pow(xke_legacy / n_rad_min, 2.0/3.0); // Earth radii

        math::TrackedValue<T> a_tv = math::TrackedValue<T>(a_er, 0.0, 1e-10, 0.0);

        perturbation::BrouwerSecularRates<T> rates = perturbation::compute_secular_rates(
            elements.mean_motion, a_tv, e2, cos_i, J2, J4);

        std::cout << "  Brouwer M_dot:     " << std::setprecision(12) << rates.M_dot.value
                  << " rad/min (accuracy=" << std::scientific << rates.M_dot.errors.accuracy << ")\n";
        std::cout << "  Brouwer omega_dot: " << std::setprecision(12) << rates.omega_dot.value
                  << " rad/min\n";
        std::cout << "  Brouwer Omega_dot: " << std::setprecision(12) << rates.Omega_dot.value
                  << " rad/min\n";

        // Period check: deep space if period >= 225 min
        double period = 2.0 * 3.14159265358979 / n_rad_min;
        std::cout << "  Period: " << std::setprecision(4) << period << " min"
                  << (period >= 225 ? " (DEEP SPACE)" : " (NEAR EARTH)") << "\n\n";
    }

    // --- Kaula Inclination Function Tests ---
    std::cout << "--- Kaula Inclination Functions ---\n";

    // Test with ISS-like inclination: 51.6443 degrees
    double inc_rad = 51.6443 * 3.14159265358979 / 180.0;
    math::TrackedValue<T> sin_i = math::TrackedValue<T>(std::sin(inc_rad), 0.0, 1e-15, 0.0);
    math::TrackedValue<T> cos_i_test = math::TrackedValue<T>(std::cos(inc_rad), 0.0, 1e-15, 0.0);

    // F_220 = 3/4 * (1 + cos i)²
    math::TrackedValue<T> f220 = perturbation::inclination_function(2, 2, 0, sin_i, cos_i_test);
    double f220_legacy = 0.75 * std::pow(1.0 + std::cos(inc_rad), 2);
    check_magic_number("F_220(51.6°)",
        f220.value, f220_legacy, 1e-12);

    // F_221 = 3/2 * sin²(i)
    math::TrackedValue<T> f221 = perturbation::inclination_function(2, 2, 1, sin_i, cos_i_test);
    double f221_legacy = 1.5 * std::pow(std::sin(inc_rad), 2);
    check_magic_number("F_221(51.6°)",
        f221.value, f221_legacy, 1e-12);

    // F_442 = 315/8 * sin⁴(i) — tests the exact rational 39.375 = 315/8
    math::TrackedValue<T> f442 = perturbation::inclination_function(4, 4, 2, sin_i, cos_i_test);
    double f442_legacy = 39.3750 * std::pow(std::sin(inc_rad), 4);
    check_magic_number("F_442(51.6°) [315/8 vs 39.3750]",
        f442.value, f442_legacy, 1e-10);

    // --- Model Selector Tests ---
    std::cout << "\n--- Model Selector ---\n";

    // Test all presets compile and produce valid configurations
    const char* presets[] = {
        "sgp4_standard", "sgp4_wgs72_old", "sgp4_wgs84",
        "modern_2020", "research_full", "eleven"
    };

    for (const char* preset : presets) {
        sgp4::ModelConfiguration<T> config = sgp4::ModelSelector<T>::select(preset, 1e-15);
        std::cout << "  " << preset << ": "
                  << config.description << "\n"
                  << "    J2=" << std::setprecision(10) << config.ellipsoid.J2n(1).value
                  << "  a=" << config.ellipsoid.a.value
                  << "  GM=" << config.ellipsoid.GM.value << "\n";
    }

    // Test custom builder
    sgp4::ModelConfiguration<T> custom = sgp4::ModelSelector<T>::custom()
        .gravity("wgs84_precise")
        .astronomy("almanac_2025")
        .kepler("householder")
        .build(1e-15);
    std::cout << "  custom: " << custom.description << "\n"
              << "    J2=" << std::setprecision(10) << custom.ellipsoid.J2n(1).value
              << "  GM=" << custom.ellipsoid.GM.value << "\n";

    // --- Diagnostic: Element Recovery for Sat 00005 ---
    if (!test_cases.empty()) {
        diagnose_sat00005<T>(test_cases[0].tle);
    }

    // --- Propagation Test: First SGP4-VER.TLE case ---
    std::cout << "\n--- Propagation Test: Sat 00005 ---\n";

    if (!test_cases.empty()) {
        TestCase& tc = test_cases[0];
        tle::TleElements<T> tle_elems = tle::TleElements<T>::from_tle_data(tc.tle);

        // Use standard SGP4 configuration
        sgp4::ModelConfiguration<T> sgp4_config = sgp4::ModelSelector<T>::select("sgp4_standard", 1e-12);

        // Create propagator
        sgp4::Propagator<T> prop(sgp4_config, tle_elems, 1e-12);

        std::cout << "  Satellite: " << tc.tle.satellite_id
                  << "  Deep space: " << (prop.is_deep_space() ? "YES" : "NO") << "\n";

        // Reference output for sat 00005 at t=0 [from tcppver.out]:
        // x=7022.46529266  y=-1400.08296755  z=0.03995155
        // vx=1.893841015  vy=6.405893759  vz=4.534807250
        double ref_x = 7022.46529266, ref_y = -1400.08296755, ref_z = 0.03995155;
        double ref_vx = 1.893841015, ref_vy = 6.405893759, ref_vz = 4.534807250;

        // Propagate at t=0
        math::TrackedValue<T> t0 = math::TrackedValue<T>::exact_integer(0);
        sgp4::StateVector<T> sv = prop.propagate(t0);

        std::cout << std::setprecision(8) << std::fixed;
        std::cout << "  t=0 min:\n";
        std::cout << "    pos: " << sv.position_km.x.value << "  "
                  << sv.position_km.y.value << "  "
                  << sv.position_km.z.value << "\n";
        std::cout << "    ref: " << ref_x << "  " << ref_y << "  " << ref_z << "\n";
        std::cout << "    vel: " << sv.velocity_km_s.x.value << "  "
                  << sv.velocity_km_s.y.value << "  "
                  << sv.velocity_km_s.z.value << "\n";
        std::cout << "    ref: " << ref_vx << "  " << ref_vy << "  " << ref_vz << "\n";

        double dx = std::abs(sv.position_km.x.value - ref_x);
        double dy = std::abs(sv.position_km.y.value - ref_y);
        double dz = std::abs(sv.position_km.z.value - ref_z);
        double pos_err = std::sqrt(dx*dx + dy*dy + dz*dz);
        std::cout << "    position error: " << std::scientific << pos_err << " km\n";

        // Multi-time-point verification against tcppver.out reference
        // Sat 00005: t, x, y, z, vx, vy, vz from reference output
        struct RefPoint { double t, x, y, z, vx, vy, vz; };
        RefPoint refs[] = {
            {    0.0,  7022.46529266, -1400.08296755,     0.03995155,  1.893841015,  6.405893759,  4.534807250},
            {  360.0, -7154.03120202, -3783.17682504, -3536.19412294,  4.741887409, -4.151817765, -2.093935425},
            {  720.0, -7134.59340119,  6531.68641334,  3260.27186483, -4.113793027, -2.911922039, -2.557327851},
            { 1080.0,  5568.53901181,  4492.06992591,  3863.87641983, -4.209106476,  5.159719888,  2.744852980},
            { 1440.0,  -938.55923943, -6268.18748831, -4294.02924751,  7.536105209, -0.427127707,  0.989878080},
        };

        std::cout << "\n  Multi-point verification:\n";
        std::cout << std::setw(10) << "t[min]" << std::setw(14) << "pos_err[km]"
                  << std::setw(14) << "vel_err[km/s]" << "\n";

        for (RefPoint& ref : refs) {
            math::TrackedValue<T> t_val = math::TrackedValue<T>(T(ref.t), T(0), T(0), T(0));
            sgp4::StateVector<T> sv_ref = prop.propagate(t_val);
            double dx2 = sv_ref.position_km.x.value - ref.x;
            double dy2 = sv_ref.position_km.y.value - ref.y;
            double dz2 = sv_ref.position_km.z.value - ref.z;
            double pe = std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
            double dvx = sv_ref.velocity_km_s.x.value - ref.vx;
            double dvy = sv_ref.velocity_km_s.y.value - ref.vy;
            double dvz = sv_ref.velocity_km_s.z.value - ref.vz;
            double ve = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
            std::cout << std::fixed << std::setprecision(1) << std::setw(10) << ref.t
                      << std::scientific << std::setprecision(4) << std::setw(14) << pe
                      << std::setw(14) << ve << "\n";
        }
    }

    // =================================================================
    // COEFFICIENT DIAGNOSTIC FOR SAT 06251
    // =================================================================
    std::cout << "\n--- Coefficient Diagnostic: Sat 06251 ---\n";
    for (size_t idx = 0; idx < test_cases.size(); ++idx) {
        if (test_cases[idx].tle.satellite_id != "06251") continue;
        TestCase& tc6 = test_cases[idx];
        tle::TleElements<T> tle6 = tle::TleElements<T>::from_tle_data(tc6.tle);
        sgp4::ModelConfiguration<T> cfg6 = sgp4::ModelSelector<T>::select("sgp4_standard", 1e-12);
        sgp4::NearSpaceInit<T> ns6 = sgp4::initialize_near_space(cfg6, tle6, 1e-12);

        // Compute reference coefficients using same formulas but raw double
        double kXKMPER = double(cfg6.ellipsoid.a.value);
        double kGM = double(cfg6.ellipsoid.GM.value);
        double kXKE = 60.0/std::sqrt(kXKMPER*kXKMPER*kXKMPER/kGM);
        double kJ2 = double(cfg6.Jn(2).value);
        double kCK2 = kJ2/2.0;
        double kJ4 = double(cfg6.Jn(4).value);

        double n0 = tc6.tle.mean_motion_rev_day * 2*3.14159265358979/1440.0;
        double e0 = tc6.tle.eccentricity;
        double i0 = tc6.tle.inclination_deg * 3.14159265358979/180.0;
        double cosio = std::cos(i0);
        double theta2 = cosio*cosio;
        double x3thm1 = 3*theta2-1;
        double betao2 = 1-e0*e0;
        double betao = std::sqrt(betao2);
        double a1 = std::pow(kXKE/n0, 2.0/3.0);
        double temp = 1.5*kCK2*x3thm1/(betao*betao2);
        double del1 = temp/(a1*a1);
        double a0 = a1*(1-del1*(1.0/3.0+del1*(1+134.0/81.0*del1)));
        double del0 = temp/(a0*a0);
        double n0pp = n0/(1+del0);
        double a0pp = a0/(1-del0);

        compare("a0pp", double(ns6.a0.value), a0pp);
        compare("n0pp", double(ns6.n0.value), n0pp);

        double s4 = 1+78.0/kXKMPER;
        double q0 = 1+120.0/kXKMPER;
        double qoms24 = std::pow(q0-s4, 4);
        double xi = 1.0/(a0pp-s4);
        double eta = a0pp*e0*xi;
        double eta2 = eta*eta;
        double eeta = e0*eta;
        double psisq = std::abs(1-eta2);
        double coef = qoms24*std::pow(xi,4);
        double coef1 = coef/std::pow(psisq, 3.5);

        compare("xi", double(ns6.xi.value), xi);
        compare("eta", double(ns6.eta.value), eta);

        // 0.75 * kCK2 = 0.375 * J2 (Vallado); cross-check dnwrnr SGP4.cc:131.
        double C2 = coef1*n0pp*(a0pp*(1+1.5*eta2+eeta*(4+eta2))
            +0.75*kCK2*xi/psisq*x3thm1*(8+3*eta2*(8+eta2)));
        double C1 = tc6.tle.bstar*C2;

        compare("C2", double(ns6.C2.value), C2);
        compare("C1", double(ns6.C1.value), C1);

        // Secular rates
        double pinvsq = 1.0/(a0pp*a0pp*betao2*betao2);
        double t1 = 3*kCK2*pinvsq*n0pp;
        double t2 = t1*kCK2*pinvsq;
        double t3 = -(15.0/32.0)*kJ4*pinvsq*pinvsq*n0pp;
        double mdot = n0pp + 0.5*t1*betao*x3thm1
            + (1.0/16.0)*t2*betao*(13+theta2*(-78+137*theta2));
        compare("M_dot", double(ns6.M_dot.value), mdot);

        // Now propagate at t=120 and compare intermediates
        double t = 120.0;
        double tsq = t*t, tcu = tsq*t, tfo = tcu*t;

        double omgdot = double(ns6.omega_dot.value);
        double nodedot = double(ns6.Omega_dot.value);
        double omega0 = tc6.tle.arg_perigee_deg*3.14159265358979/180.0;
        double Omega0 = tc6.tle.raan_deg*3.14159265358979/180.0;
        double M0 = tc6.tle.mean_anomaly_deg*3.14159265358979/180.0;

        double xmdf = M0 + mdot*t;
        double omgadf = omega0 + omgdot*t;
        double xnoddf = Omega0 + nodedot*t;
        double xnodcf = double(ns6.Omega_dot_nkc.value);
        double xnode = xnoddf + xnodcf*tsq;

        double tempa_ref = 1 - C1*t;
        double tempe_ref = tc6.tle.bstar*double(ns6.C4.value)*t;
        double templ_ref = double(ns6.t2cof.value)*tsq;

        // Non-simple model corrections
        bool simple = double(ns6.perigee_km.value) < 220.0;
        std::cout << "  simple_model=" << (simple ? "YES" : "NO")
                  << "  perigee=" << double(ns6.perigee_km.value) << " km\n";

        if (!simple) {
            double omgcof_v = double(ns6.omgcof.value);
            double xmcof_v = double(ns6.xmcof.value);
            double delmo_v = double(ns6.delmo.value);
            double eta_v = double(ns6.eta.value);
            double sinmo_v = double(ns6.sinmo.value);

            double delomg = omgcof_v*t;
            double delm = xmcof_v*(std::pow(1+eta_v*std::cos(xmdf),3) - delmo_v);
            double xmp = xmdf + delomg + delm;
            double omega = omgadf - delomg - delm;

            tempa_ref -= double(ns6.D2.value)*tsq + double(ns6.D3.value)*tcu + double(ns6.D4.value)*tfo;
            tempe_ref += tc6.tle.bstar*double(ns6.C5.value)*(std::sin(xmp)-sinmo_v);
            templ_ref += double(ns6.t3cof.value)*tcu + tfo*(double(ns6.t4cof.value)+t*double(ns6.t5cof.value));
        }

        // Our propagation
        math::TrackedValue<T> t_tv = math::TrackedValue<T>(T(120.0), T(0), T(0), T(0));
        sgp4::StateVector<T> sv6 = sgp4::propagate_near_space(ns6, cfg6, t_tv, T(1e-12));

        // Compare tempa/tempe/templ — but we can't access them directly.
        // Instead compare a and e:
        double a_ref = a0pp*tempa_ref*tempa_ref;
        double e_ref = e0 - tempe_ref;
        double xl_ref = xmdf + omgadf + xnoddf + xnodcf*tsq + n0pp*templ_ref;

        // Our a: back-compute from position (can't access directly)
        // Just show the reference values and final position error
        std::cout << "  t=120: tempa=" << std::setprecision(15) << tempa_ref
                  << "  tempe=" << tempe_ref << "  templ=" << templ_ref << "\n";
        std::cout << "  t=120: a_ref=" << a_ref << "  e_ref=" << e_ref << "\n";

        // FULL raw-double propagation at t=120 following exact FORTRAN sequence.
        // Compare every intermediate against our TrackedValue path.
        {
            double ref_x6=-3935.69800083, ref_y6=409.10980837, ref_z6=5471.33577327;

            // Our TrackedValue propagation
            math::TrackedValue<T> t_tv = math::TrackedValue<T>(T(120.0), T(0), T(0), T(0));
            sgp4::StateVector<T> sv6 = sgp4::propagate_near_space(ns6, cfg6, t_tv, T(1e-12));
            double our_x = sv6.position_km.x.value;
            double our_y = sv6.position_km.y.value;
            double our_z = sv6.position_km.z.value;
            double our_err = std::sqrt(std::pow(our_x-ref_x6,2)+std::pow(our_y-ref_y6,2)+std::pow(our_z-ref_z6,2));

            // Raw-double propagation matching FORTRAN sequence exactly
            double kXKE = double(ns6.xke.value);
            double kCK2 = double(ns6.CK2.value);
            double kXKMPER = double(ns6.re_km.value);
            double a0d = double(ns6.a0.value);
            double n0d = double(ns6.n0.value);
            double e0d = double(ns6.e0.value);
            double i0d = double(ns6.i0.value);
            double w0d = double(ns6.omega0.value);
            double O0d = double(ns6.Omega0.value);
            double M0d = double(ns6.M0.value);
            double bsd = double(ns6.bstar.value);
            double mdotd = double(ns6.M_dot.value);
            double wdotd = double(ns6.omega_dot.value);
            double Odotd = double(ns6.Omega_dot.value);

            double C1d = double(ns6.C1.value);
            double C4d = double(ns6.C4.value);
            double C5d = double(ns6.C5.value);
            double D2d = double(ns6.D2.value);
            double D3d = double(ns6.D3.value);
            double D4d = double(ns6.D4.value);
            double t2d = double(ns6.t2cof.value);
            double t3d = double(ns6.t3cof.value);
            double t4d = double(ns6.t4cof.value);
            double t5d = double(ns6.t5cof.value);
            double omgcofd = double(ns6.omgcof.value);
            double xmcofd = double(ns6.xmcof.value);
            double delmod = double(ns6.delmo.value);
            double sinmod = double(ns6.sinmo.value);
            double xnodcfd = double(ns6.Omega_dot_nkc.value);
            double etad = double(ns6.eta.value);
            double xlcofd = double(ns6.xlcof.value);
            double aycofd = double(ns6.aycof.value);
            double cosiod = double(ns6.cosio.value);
            double siniod = double(ns6.sinio.value);
            double theta2d = double(ns6.theta2.value);
            double x3thm1d = double(ns6.x3thm1.value);
            double x1mth2d = double(ns6.x1mth2.value);
            double x7thm1d = double(ns6.x7thm1.value);

            double td = 120.0;
            double tsqd = td*td, tcud = tsqd*td, tfod = tcud*td;

            // Secular
            double xmdf_d = M0d + mdotd*td;
            double omgadf_d = w0d + wdotd*td;
            double xnoddf_d = O0d + Odotd*td;
            double xnode_d = xnoddf_d + xnodcfd*tsqd;

            double tempa_d = 1.0 - C1d*td;
            double tempe_d = bsd*C4d*td;
            double templ_d = t2d*tsqd;

            // Non-simple corrections
            double delomg_d = omgcofd*td;
            double delm_d = xmcofd*(std::pow(1.0+etad*std::cos(xmdf_d),3.0)-delmod);
            double temp_drag_d = delomg_d + delm_d;
            double xmp_d = xmdf_d + temp_drag_d;
            double omega_d = omgadf_d - temp_drag_d;

            tempa_d -= D2d*tsqd + D3d*tcud + D4d*tfod;
            tempe_d += bsd*C5d*(std::sin(xmp_d)-sinmod);
            templ_d += t3d*tcud + tfod*(t4d + td*t5d);

            double a_d = a0d*tempa_d*tempa_d;
            double e_d = e0d - tempe_d;
            if (e_d < 1e-6) e_d = 1e-6;

            double xl_d = xmp_d + omega_d + xnode_d + n0d*templ_d;

            // Long-period
            double beta2_d = 1.0-e_d*e_d;
            double n_d = kXKE/(a_d*std::sqrt(a_d));
            double axN_d = e_d*std::cos(omega_d);
            double temp_lp_d = 1.0/(a_d*beta2_d);
            double xll_d = temp_lp_d*xlcofd*axN_d;
            double aynl_d = temp_lp_d*aycofd;
            double xlt_d = xl_d + xll_d;
            double ayn_d = e_d*std::sin(omega_d) + aynl_d;

            // Kepler
            double capu_d = std::fmod(xlt_d - xnode_d, 2.0*3.14159265358979323846);
            if (capu_d < 0) capu_d += 2.0*3.14159265358979323846;
            double epw_d = capu_d;
            for (int iter = 0; iter < 30; ++iter) {
                double sx = std::sin(epw_d), cx = std::cos(epw_d);
                double f = capu_d - ayn_d*cx + axN_d*sx - epw_d;
                double fp = -ayn_d*sx - axN_d*cx + 1.0;
                double delta = f/fp;
                epw_d += delta;
                if (std::abs(delta) < 1e-12) break;
            }

            // Osculating
            double sinepw_d = std::sin(epw_d), cosepw_d = std::cos(epw_d);
            double ecose_d = axN_d*cosepw_d + ayn_d*sinepw_d;
            double esine_d = axN_d*sinepw_d - ayn_d*cosepw_d;
            double el2_d = axN_d*axN_d + ayn_d*ayn_d;
            double pl_d = a_d*(1.0-el2_d);
            double r_d = a_d*(1.0-ecose_d);
            double rdot_d = kXKE*std::sqrt(a_d)*esine_d/r_d;
            double rfdot_d = kXKE*std::sqrt(pl_d)/r_d;
            double betal_d = std::sqrt(1.0-el2_d);
            double tk_d = 1.0/(1.0+betal_d);
            double cosu_d = (a_d/r_d)*(cosepw_d-axN_d+ayn_d*esine_d*tk_d);
            double sinu_d = (a_d/r_d)*(sinepw_d-ayn_d-axN_d*esine_d*tk_d);
            double u_d = std::atan2(sinu_d,cosu_d);
            double sin2u_d = 2*sinu_d*cosu_d;
            double cos2u_d = 2*cosu_d*cosu_d-1;

            // Short-period
            double sp2_d = kCK2/pl_d;
            double sp3_d = sp2_d/pl_d;
            double rk_d = r_d*(1-1.5*sp3_d*betal_d*x3thm1d)+0.5*sp2_d*x1mth2d*cos2u_d;
            double uk_d = u_d-0.25*sp3_d*x7thm1d*sin2u_d;
            double xnodek_d = xnode_d+1.5*sp3_d*cosiod*sin2u_d;
            double xinck_d = i0d+1.5*sp3_d*cosiod*siniod*cos2u_d;
            double rdotk_d = rdot_d-n_d*sp2_d*x1mth2d*sin2u_d;
            double rfdotk_d = rfdot_d+n_d*sp2_d*(x1mth2d*cos2u_d+1.5*x3thm1d);

            // Position
            double sinuk_d=std::sin(uk_d),cosuk_d=std::cos(uk_d);
            double sinik_d=std::sin(xinck_d),cosik_d=std::cos(xinck_d);
            double sinnok_d=std::sin(xnodek_d),cosnok_d=std::cos(xnodek_d);
            double xmx_d=-sinnok_d*cosik_d, xmy_d=cosnok_d*cosik_d;
            double ux_d=xmx_d*sinuk_d+cosnok_d*cosuk_d;
            double uy_d=xmy_d*sinuk_d+sinnok_d*cosuk_d;
            double uz_d=sinik_d*sinuk_d;
            double x_d=rk_d*ux_d*kXKMPER;
            double y_d=rk_d*uy_d*kXKMPER;
            double z_d=rk_d*uz_d*kXKMPER;

            double raw_err = std::sqrt(std::pow(x_d-ref_x6,2)+std::pow(y_d-ref_y6,2)+std::pow(z_d-ref_z6,2));
            double our_vs_raw = std::sqrt(std::pow(our_x-x_d,2)+std::pow(our_y-y_d,2)+std::pow(our_z-z_d,2));

            std::cout << std::scientific << std::setprecision(6);
            std::cout << "  t=120 comparison:\n";
            std::cout << "    our_err vs tcppver:  " << our_err << " km\n";
            std::cout << "    raw_err vs tcppver:  " << raw_err << " km\n";
            std::cout << "    our vs raw:          " << our_vs_raw << " km\n";
            std::cout << "    raw pos: " << std::setprecision(11) << x_d << "  " << y_d << "  " << z_d << "\n";
            std::cout << "    our pos: " << our_x << "  " << our_y << "  " << our_z << "\n";
            std::cout << "    ref pos: " << ref_x6 << "  " << ref_y6 << "  " << ref_z6 << "\n";
        }

        std::cout << "\n";
        break;
    }

    // =================================================================
    // FULL VALIDATION: All 33 SGP4-VER.TLE cases against tcppver.out
    // =================================================================
    std::cout << "\n=== Full Validation: All SGP4-VER.TLE Cases ===\n";

    std::vector<ReferenceCase> ref_cases = parse_tcppver(
        "sgp4_references/aholinch_sgp4/data/tcppver.out");
    std::cout << "Loaded " << ref_cases.size() << " reference cases from tcppver.out\n\n";

    sgp4::ModelConfiguration<T> sgp4_config = sgp4::ModelSelector<T>::select("sgp4_standard", 1e-12);

    int total_pass = 0, total_fail = 0, total_skip = 0;
    int total_points_pass = 0, total_points_fail = 0;

    // Match reference cases to TLE cases by satellite ID
    for (size_t i = 0; i < test_cases.size(); ++i) {
        TestCase& tc = test_cases[i];

        // Find matching reference case.
        // TLE stores "00005", tcppver.out stores "5".
        // Try integer comparison (handles leading zeros); fall back to string
        // for Alpha-5 IDs like "A0005" that aren't purely numeric.
        const ReferenceCase* ref = nullptr;
        bool (*try_parse_int)(const std::string&, int&) = [](const std::string& s, int& out) -> bool {
            try { out = std::stoi(s); return true; }
            catch (...) { return false; }
        };
        int tle_num = 0;
        bool tle_is_int = try_parse_int(tc.tle.satellite_id, tle_num);

        for (ReferenceCase& rc : ref_cases) {
            if (tle_is_int) {
                int ref_num = 0;
                if (try_parse_int(rc.satellite_id, ref_num) && ref_num == tle_num) {
                    ref = &rc;
                    break;
                }
            }
            // String fallback
            if (rc.satellite_id == tc.tle.satellite_id) {
                ref = &rc;
                break;
            }
        }

        tle::TleElements<T> tle_elems = tle::TleElements<T>::from_tle_data(tc.tle);

        // Try to create propagator — deep-space may work or may have issues
        try {
            sgp4::Propagator<T> prop(sgp4_config, tle_elems, 1e-12);

            double max_pos_err = 0;
            double t0_err = -1;
            int points_pass = 0, points_fail = 0;

            if (ref && !ref->points.empty()) {
                for (const ReferencePoint& pt : ref->points) {
                    // Expected-error point: the reference errored here (NaN), so
                    // the correct behavior is for our propagator to error too.
                    bool expect_error = ref->t0_error && pt.t_min == 0.0;
                    math::TrackedValue<T> t_val = math::TrackedValue<T>(T(pt.t_min), T(0), T(0), T(0));
                    try {
                        sgp4::StateVector<T> sv = prop.propagate(t_val);
                        bool we_errored = std::isnan(double(sv.position_km.x.value));
                        if (expect_error) {
                            // Pass iff we also errored (matched the reference error).
                            if (pt.t_min == 0.0) t0_err = we_errored ? 0.0 : -1.0;
                            if (we_errored) points_pass++; else points_fail++;
                            continue;
                        }
                        if (we_errored) { points_fail++; continue; }  // unexpected error
                        double dx = sv.position_km.x.value - pt.x;
                        double dy = sv.position_km.y.value - pt.y;
                        double dz = sv.position_km.z.value - pt.z;
                        double pe = std::sqrt(dx*dx + dy*dy + dz*dz);
                        if (pe > max_pos_err) max_pos_err = pe;
                        if (pt.t_min == 0.0) t0_err = pe;

                        // Tolerance: 1e-3 km for near-earth, 1e-1 km for deep-space
                        double tol = prop.is_deep_space() ? 0.1 : 0.001;
                        if (pe < tol) points_pass++; else points_fail++;
                    } catch (const std::exception& ex) {
                        // A thrown error also satisfies an expected-error point.
                        if (expect_error) points_pass++; else points_fail++;
                    }
                }
            }

            bool pass = (points_fail == 0 && ref && !ref->points.empty());
            std::cout << (pass ? "PASS" : "FAIL") << "  sat " << tc.tle.satellite_id
                      << "  " << (prop.is_deep_space() ? "DEEP" : "NEAR")
                      << std::scientific << std::setprecision(4)
                      << "  t0=" << t0_err
                      << "  max=" << max_pos_err << " km"
                      << "  bstar=" << tc.tle.bstar
                      << "  (" << points_pass << "/" << (points_pass+points_fail) << ")\n";

            // For failing cases (near-earth OR deep-space), show per-point error growth
            if (!pass && ref && max_pos_err > 0.01) {
                std::cout << "    Error growth for sat " << tc.tle.satellite_id << ":\n";
                int shown = 0;
                for (const ReferencePoint& pt : ref->points) {
                    if (shown >= 8) break;
                    math::TrackedValue<T> t_v = math::TrackedValue<T>(T(pt.t_min), T(0), T(0), T(0));
                    try {
                        sgp4::StateVector<T> sv2 = prop.propagate(t_v);
                        double pe2 = std::sqrt(
                            std::pow(sv2.position_km.x.value - pt.x, 2) +
                            std::pow(sv2.position_km.y.value - pt.y, 2) +
                            std::pow(sv2.position_km.z.value - pt.z, 2));
                        std::cout << "      t=" << std::fixed << std::setprecision(1) << pt.t_min
                                  << "  err=" << std::scientific << std::setprecision(6) << pe2
                                  << "  dx=" << sv2.position_km.x.value - pt.x
                                  << "  dy=" << sv2.position_km.y.value - pt.y
                                  << "  dz=" << sv2.position_km.z.value - pt.z << "\n";
                    } catch (...) {
                        std::cout << "      t=" << pt.t_min << "  EXCEPTION\n";
                    }
                    shown++;
                }
            }

            if (pass) total_pass++; else total_fail++;
            total_points_pass += points_pass;
            total_points_fail += points_fail;

        } catch (const std::exception& ex) {
            std::cout << "SKIP  sat " << tc.tle.satellite_id << "  " << ex.what() << "\n";
            total_skip++;
        }
    }

    std::cout << "\n--- Summary ---\n";
    std::cout << "  Satellites: " << total_pass << " pass, " << total_fail << " fail, "
              << total_skip << " skip (of " << test_cases.size() << ")\n";
    std::cout << "  Data points: " << total_points_pass << " pass, " << total_points_fail << " fail\n";

    std::cout << "\n=== Done ===\n";
    return 0;
}
