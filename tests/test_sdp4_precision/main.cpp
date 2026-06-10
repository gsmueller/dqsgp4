/// test_sdp4_precision — DS1 (DQSGP4 Completion Roadmap, register §DS1).
///
/// The DQSGP4 calling card is a rigorously-tracked accuracy/precision budget
/// that propagates through the WHOLE computation and tightens with a wider
/// numeric type T. Before DS1 the SDP4 deep-space evolution ran in raw double
/// and only re-injected the epoch budget at the Kepler stage, so for wide T the
/// reported precision shrank but the computation never actually achieved it.
///
/// DS1 made the deep-space evolution (DSCOM coefficient buildup, DPPER periodic
/// corrections, the secular advance, and the resonance leapfrog) run entirely
/// in TrackedValue<T>. This gate asserts the consequence: for a deep-space
/// satellite, the propagated position-precision budget at cpp_bin_float_50 is
/// dramatically tighter than at double — AND the values agree (same model, just
/// computed to more digits). Three sats cover every deep-space path:
///   11801  — non-resonant SDP4 (DSCOM + DPPER + secular only);
///   08195  — 12h half-day resonant (the half-day leapfrog branch);
///   24208  — 24h synchronous resonant (the synchronous leapfrog branch).
///
/// ExeGate DS1: returns a nonzero exit code on any failed check.

#include "sgp4/sgp4_propagator.h"
#include "sgp4/model_selector.h"
#include "tle/tle_parser.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_bin_float_50;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

struct DsCase {
    const char* name;
    const char* l1;
    const char* l2;
    int tmin;  // minutes since epoch (integer => exact input)
};

// Position magnitude [km] and its propagated precision + accuracy budgets.
template<typename T>
void propagate(const tle::TleData& td, int tmin,
               double& pos_mag, double& pos_precision, double& pos_accuracy,
               bool& deep, bool& finite) {
    using TV = math::TrackedValue<T>;
    tle::TleElements<T> elems = tle::TleElements<T>::from_tle_data(td);
    sgp4::ModelConfiguration<T> cfg = sgp4::ModelSelector<T>::select("sgp4_standard", T(1e-12));
    sgp4::Propagator<T> prop(cfg, elems, T(1e-12));
    deep = prop.is_deep_space();

    TV t = TV::exact_integer(tmin);
    sgp4::StateVector<T> sv = prop.propagate(t);

    TV mag = sv.position_km.magnitude();
    finite = std::isfinite(static_cast<double>(mag.value));
    pos_mag = static_cast<double>(mag.value);
    pos_precision = static_cast<double>(mag.errors.precision);
    pos_accuracy = static_cast<double>(mag.errors.accuracy);
}

// Returns the precision-scaling ratio (double/bf50) for the case, or 0 on a
// parse failure, so main can assert at least one case scales strongly.
double run(const DsCase& c) {
    tle::TleData td;
    bool ok = tle::parse(std::string(c.l1), std::string(c.l2), td);
    check(std::string(c.name) + ": TLE parses", ok);
    if (!ok) return 0.0;

    double mag_d = 0, prec_d = 0, acc_d = 0, mag_b = 0, prec_b = 0, acc_b = 0;
    bool deep_d = false, deep_b = false, fin_d = false, fin_b = false;
    propagate<double>(td, c.tmin, mag_d, prec_d, acc_d, deep_d, fin_d);
    propagate<cpp_bin_float_50>(td, c.tmin, mag_b, prec_b, acc_b, deep_b, fin_b);

    check(std::string(c.name) + ": is deep space", deep_d && deep_b);
    check(std::string(c.name) + ": finite state at t", fin_d && fin_b);

    // Same model: the double and cpp_bin_float_50 position magnitudes agree to
    // well within double's own rounding-amplified error over the arc.
    double rel = std::abs(mag_d - mag_b) / std::max(mag_b, 1.0);
    check(std::string(c.name) + ": double and bf50 values agree", rel < 1e-6);

    // The propagated precision budget is real (nonzero, the WHOLE deep-space
    // evolution now runs in T) ...
    check(std::string(c.name) + ": precision budget > 0", prec_d > 0.0 && prec_b > 0.0);
    // ... and never WORSENS with a wider T (more mantissa bits cannot reduce
    // computational precision) ...
    check(std::string(c.name) + ": precision does not worsen with wider T",
          prec_b <= prec_d * (1.0 + 1e-9));
    // ... and the headline honest statement: at cpp_bin_float_50 the
    // COMPUTATIONAL precision is far below the MODEL-accuracy floor, so the
    // tracked budget proves the deep-space result is limited by the SR3 model,
    // not by the arithmetic. (Pre-DS1 the precision was the re-injected epoch
    // budget — it could not make this statement.)
    check(std::string(c.name) + ": computation far below model floor at bf50",
          prec_b * 100.0 < acc_b);

    // The accuracy budget is the SR3 model-fidelity floor: nonzero and
    // T-INDEPENDENT (the lunisolar coefficients are specified to finite digits,
    // so a wider T cannot sharpen the MODEL — only the computation). It lives in
    // accuracy, not precision, via TrackedValue::model_coefficient.
    check(std::string(c.name) + ": accuracy floor > 0 (model fidelity)", acc_d > 0.0 && acc_b > 0.0);
    double acc_rel = std::abs(acc_d - acc_b) / std::max(acc_b, 1e-30);
    check(std::string(c.name) + ": accuracy floor is T-independent", acc_rel < 1e-6);

    double ratio = prec_b > 0.0 ? prec_d / prec_b : 0.0;
    std::cout << "    " << c.name << ": |r|=" << mag_b
              << " km | precision double=" << prec_d << " bf50=" << prec_b
              << " (x" << ratio << " tighter) | accuracy=" << acc_d << " km\n";
    return ratio;
}

} // namespace

int main() {
    // Deep-space SGP4-VER cases (lines verbatim from SGP4-VER.TLE; parse reads
    // fixed columns and ignores the trailing start/stop/step fields).
    DsCase cases[] = {
        {"11801 (non-resonant SDP4)",
         "1 11801U          80230.29629788  .01431103  00000-0  14311-1      13",
         "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13      0.0      1440.0        360.00",
         1440},
        {"08195 (12h half-day resonant)",
         "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813",
         "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656      0.0      2880.0        120.00",
         2880},
        {"24208 (24h synchronous resonant)",
         "1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600",
         "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119      0.0      1440.0        120.00",
         1440},
    };

    double best_ratio = 0.0;
    for (const DsCase& c : cases) best_ratio = std::max(best_ratio, run(c));

    // The deep-space evolution genuinely propagates type-T precision: on a clean
    // (low-drag, low-J3) deep-space case the precision sharpens by >100x going
    // double -> cpp_bin_float_50. (High-ecc/high-inc cases retain a residual
    // floor from the J3 / earth-model book coefficients, which the CR1-b
    // constant-honesty sweep moves to accuracy via model_coefficient too.)
    check("at least one deep-space case: precision sharpens >100x with wider T",
          best_ratio > 100.0);

    std::cout << "\n  SDP4 precision: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
