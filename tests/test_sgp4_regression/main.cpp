/// test_sgp4_regression — frozen SGP4 self-consistency oracle (TEMPORARY).
///
/// The 33/33 reference test compares our SGP4 against Vallado's tcppver.out to a
/// kilometre-scale tolerance, so a refactor could shift our outputs WITHIN that
/// tolerance and pass undetected. This oracle is the stronger guard the redesign
/// needs: it captures the CURRENT (verified, 33/33) SGP4's full-precision outputs
/// over the SGP4-VER battery into a golden file on first run, then on every later
/// run re-propagates and asserts the live SGP4 reproduces the frozen baseline to
/// a tolerance FAR tighter than the reference test. SGP4 is bit-deterministic
/// (finite outputs are identical across rebuilds), so this is a tight guard; the
/// per-point NaN error-sentinels are compared by error-state (not value).
///
/// It is the "temp version of SGP4" for regression during the common-generative-
/// astro-library redesign (notably the F3 constants unification, which touches
/// the SGP4 path): any drift fails here immediately. Refresh the golden (delete
/// it; the next run recaptures) or remove this oracle once the redesign's final
/// functions are complete.
///
/// ExeGate OR1: nonzero exit on any drift from the frozen baseline.

#include "sgp4/model_selector.h"
#include "sgp4/sgp4_propagator.h"
#include "math/tracked_value.h"
#include "tle/tle_parser.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using T = double;

struct TestCase {
    tle::TleData tle;
    double start_min;
    double stop_min;
    double step_min;
};

/// Parse SGP4-VER.TLE (start/stop/step appended after column 69). Mirrors the
/// parser in tests/test_sgp4/main.cpp so the oracle uses the same battery.
std::vector<TestCase> parse_ver_file(const std::string& path) {
    std::vector<TestCase> cases;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open " << path << "\n";
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
        } else if (line[0] == '2' && line.size() >= 69 && have_line1) {
            TestCase tc;
            std::string line2 = line.substr(0, 69);
            if (tle::parse(pending_line1, line2, tc.tle)) {
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

/// Propagate the current SGP4 over the battery, collecting full-precision
/// position/velocity at each deterministic start..stop:step time.
std::vector<std::array<double, 6>> run_battery(const std::vector<TestCase>& cases) {
    std::vector<std::array<double, 6>> out;
    for (const TestCase& tc : cases) {
        sgp4::ModelConfiguration<T> config =
            sgp4::ModelSelector<T>::select("sgp4_standard", T(1e-12));
        tle::TleElements<T> elems = tle::TleElements<T>::from_tle_data(tc.tle);
        sgp4::Propagator<T> prop(config, elems, T(1e-12));
        const double step = (tc.step_min > 0.0) ? tc.step_min : 1.0;
        for (double t = tc.start_min; t <= tc.stop_min + 1e-9; t += step) {
            math::TrackedValue<T> tsince(T(t), T(0), T(0), T(0));
            sgp4::StateVector<T> sv = prop.propagate(tsince);
            out.push_back({sv.position_km.x.value, sv.position_km.y.value,
                           sv.position_km.z.value, sv.velocity_km_s.x.value,
                           sv.velocity_km_s.y.value, sv.velocity_km_s.z.value});
            if (tc.step_min <= 0.0) break;
        }
    }
    return out;
}

/// Serialize a component: finite -> 17 significant figures (exact double
/// round-trip); non-finite (an error sentinel) -> "nan" (a token `>>` can read).
void write_comp(std::ostream& os, double v) {
    if (std::isfinite(v)) os << std::setprecision(17) << v;
    else os << "nan";
}

/// Parse a component token. A finite double never contains 'n'/'i', so any token
/// that does is a non-finite sentinel and maps to NaN.
double read_comp(const std::string& tok) {
    if (tok.find_first_of("nNiI") != std::string::npos) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::stod(tok);
}

/// Components agree if both finite and within tol, or both non-finite (both are
/// error sentinels — what matters is that the error state is reproduced).
bool agree(double a, double b, double tol) {
    bool fa = std::isfinite(a);
    bool fb = std::isfinite(b);
    if (fa && fb) return std::abs(a - b) <= tol;
    return fa == fb;
}

} // namespace

int main() {
    const std::string ver_path = "sgp4_references/aholinch_sgp4/data/SGP4-VER.TLE";
    const std::string golden_path = "tests/test_sgp4_regression/golden.txt";
    const double pos_tol = 1e-6;   // km   — 1 mm; ~1000x tighter than the 33/33 test
    const double vel_tol = 1e-9;   // km/s

    std::vector<TestCase> cases = parse_ver_file(ver_path);
    if (cases.empty()) {
        std::cerr << "no SGP4-VER cases parsed\n";
        return 1;
    }
    std::vector<std::array<double, 6>> live = run_battery(cases);

    // Capture-if-absent: the first run freezes the current SGP4 as the baseline.
    std::ifstream gin(golden_path);
    if (!gin.is_open()) {
        std::ofstream gout(golden_path);
        for (const std::array<double, 6>& o : live) {
            for (int k = 0; k < 6; ++k) {
                write_comp(gout, o[k]);
                gout << (k < 5 ? ' ' : '\n');
            }
        }
        std::cout << "CAPTURED golden baseline (" << live.size()
                  << " points) — re-baseline; PASS.\n";
        return 0;
    }

    std::size_t idx = 0;
    int drift = 0;
    double max_pos = 0.0;
    double max_vel = 0.0;
    std::string line;
    while (std::getline(gin, line) && idx < live.size()) {
        std::istringstream ss(line);
        std::array<double, 6> g{};
        for (int k = 0; k < 6; ++k) {
            std::string tok;
            ss >> tok;
            g[k] = read_comp(tok);
        }
        const std::array<double, 6>& o = live[idx];
        for (int k = 0; k < 3; ++k) {
            double d = std::abs(o[k] - g[k]);
            if (std::isfinite(d) && d > max_pos) max_pos = d;
            if (!agree(o[k], g[k], pos_tol)) ++drift;
        }
        for (int k = 3; k < 6; ++k) {
            double d = std::abs(o[k] - g[k]);
            if (std::isfinite(d) && d > max_vel) max_vel = d;
            if (!agree(o[k], g[k], vel_tol)) ++drift;
        }
        ++idx;
    }
    if (idx != live.size()) {
        std::cerr << "golden/live point-count mismatch (golden " << idx
                  << " vs live " << live.size() << ")\n";
        return 1;
    }
    std::cout << "frozen SGP4 oracle: " << live.size() << " points; max pos drift "
              << std::setprecision(6) << max_pos << " km, max vel drift " << max_vel
              << " km/s\n";
    if (drift > 0) {
        std::cout << drift << " component(s) exceeded the frozen baseline "
                  << "(tol " << pos_tol << " km / " << vel_tol << " km/s); FAIL.\n";
        return 1;
    }
    std::cout << "reproduces the frozen baseline; PASS.\n";
    return 0;
}
