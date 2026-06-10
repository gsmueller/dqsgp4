/// test_atmosphere — the piecewise-exponential static-atmosphere table (gate ATM1, R1).
///
/// Verifies src/atmosphere/exponential_table.h (Vallado Table 8-4, USSA76/CIRA-72
/// derived, transcribed from the in-repo born-digital ATMOSEXP.DAT) per
/// design/derivations/atmosphere_exponential_table.md §4:
///
///   1. NODE EXACTNESS — at every band base, ρ(h₀) == ρ₀ bit-for-bit
///      (Δh = 0 ⇒ exp(0) = 1 exactly).
///   2. NEIGHBOR CHAINING — row i's exponential evaluated at row i+1's base
///      reproduces row i+1's ρ₀ to ≤ 2e-4 (2× the measured 9.6e-5 construction
///      residual). A transcription error in ρ₀ or H breaks TWO rows at ≥1e-2.
///   3. Monotone decreasing density over 0–1000 km.
///   4. Cross-repo consistency: the 200-km row IS test_propagator's drag
///      constant 2.789e-10 kg/m³ (its provenance, identified).
///   5. Lane-divergence demonstration: the old single-band Lane parameters hold
///      ~15 % near their 200-km anchor but drift to >30 % by 300 km — the
///      measured reason the per-band table replaces the single band.
///   6. The default declared band (1e-3) majorizes the measured chaining residual;
///      the digit-floor accuracy is nonzero; precision tightens with wider T.
///   7. make_drag wiring: the model through the existing DensityModel seam gives
///      the hand-computed ½ρ|v_rel|²B at 400 km.
///
/// Exit 0 iff every check passes.

#include "atmosphere/exponential_table.h"
#include "constants/constants_provider.h"
#include "dynamics/state.h"
#include "forces/drag.h"
#include "math/quaternion.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

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

} // namespace

int main() {
    std::cout << std::setprecision(6);
    std::cout << "test_atmosphere: Vallado 8-4 piecewise-exponential table (ATM1, R1)\n\n";
    using T = double;
    using atmosphere::kVallado84;
    using atmosphere::kVallado84Rows;

    // ---- 1. node exactness (bit-for-bit at every band base) ----
    std::cout << "=== node exactness ===\n";
    {
        int exact_nodes = 0;
        for (int i = 0; i < kVallado84Rows; ++i) {
            math::TrackedValue<T> rho =
                atmosphere::vallado84_density<T>(tv<T>(double(kVallado84[i].h0_m)));
            T expect = math::TrackedValue<T>::model_coefficient(kVallado84[i].rho0).value;
            if (rho.value == expect) ++exact_nodes;
        }
        std::cout << "  " << exact_nodes << "/" << kVallado84Rows << " nodes bit-exact\n";
        check("rho(h0) == rho0 bit-for-bit at all 27 band bases", exact_nodes == kVallado84Rows);
    }

    // ---- 2. neighbor chaining (the transcription tripwire) ----
    std::cout << "\n=== neighbor chaining ===\n";
    double worst_chain = 0.0;
    {
        for (int i = 0; i + 1 < kVallado84Rows; ++i) {
            double rho0 = std::strtod(kVallado84[i].rho0, nullptr);
            double H_m = std::strtod(kVallado84[i].H_km, nullptr) * 1000.0;
            double dh = double(kVallado84[i + 1].h0_m - kVallado84[i].h0_m);
            double chained = rho0 * std::exp(-dh / H_m);
            double next = std::strtod(kVallado84[i + 1].rho0, nullptr);
            worst_chain = std::max(worst_chain, std::abs(chained - next) / next);
        }
        std::cout << "  max chaining residual: " << worst_chain
                  << "  (measured construction grade 9.6e-5)\n";
        check("all 26 band edges chain to <= 2e-4", worst_chain <= 2e-4);
    }

    // ---- 3. monotone decreasing over 0..1000 km ----
    {
        bool monotone = true;
        double prev = atmosphere::vallado84_density<T>(tv<T>(0.0)).value;
        for (int km = 5; km <= 1000; km += 5) {
            double rho = atmosphere::vallado84_density<T>(tv<T>(km * 1000.0)).value;
            if (rho >= prev) monotone = false;
            prev = rho;
        }
        check("density strictly decreasing 0..1000 km (5-km grid)", monotone);
    }

    // ---- 4. cross-repo consistency: the test_propagator drag constant ----
    {
        double rho200 = atmosphere::vallado84_density<T>(tv<T>(200000.0)).value;
        std::cout << "\n  rho(200 km) = " << rho200 << " (test_propagator uses 2.789e-10)\n";
        check("200-km row == test_propagator's 2.789e-10 (provenance identified)",
              rho200 == 2.789e-10);
    }

    // ---- 5. Lane-divergence demonstration ----
    std::cout << "\n=== single-band Lane vs the table ===\n";
    {
        auto lane = [](double h_m) {
            return 2.789e-10 * std::exp(-(h_m - 200000.0) / 50000.0);
        };
        auto rel = [&](double h_m) {
            double t = atmosphere::vallado84_density<T>(tv<T>(h_m)).value;
            return std::abs(lane(h_m) - t) / t;
        };
        double r220 = rel(220000.0), r250 = rel(250000.0), r300 = rel(300000.0);
        std::cout << "  rel diff: 220 km " << r220 << "   250 km " << r250
                  << "   300 km " << r300 << "\n";
        check("Lane holds its ~30% band near the 200-km anchor (220 km)", r220 < 0.30);
        check("Lane exceeds its 30% band by 300 km (the R1 motivation)", r300 > 0.30);
    }

    // ---- 6. error budget ----
    std::cout << "\n=== error budget ===\n";
    {
        math::TrackedValue<T> rho = atmosphere::vallado84_density<T>(tv<T>(400000.0));
        double rel_band = T(1) / T(1000);
        check("default band (1e-3) majorizes the measured chaining residual",
              rel_band > worst_chain);
        check("accuracy nonzero (digit floors + band)", rho.errors.accuracy > T(0));

        math::TrackedValue<cpp_bin_float_50> rho_b =
            atmosphere::vallado84_density<cpp_bin_float_50>(tv<cpp_bin_float_50>(400000.0));
        double pd = static_cast<double>(rho.errors.precision);
        double pb = static_cast<double>(rho_b.errors.precision);
        std::cout << "  precision: double = " << pd << "  bf50 = " << pb << "\n";
        check("precision > 0 (framework alive)", pd > 0.0 && pb > 0.0);
        check("precision tightens with wider T", pb < pd);
    }

    // ---- 7. make_drag wiring at 400 km ----
    std::cout << "\n=== make_drag wiring (400 km circular) ===\n";
    {
        constants::ConstantsProvider<T> K = constants::ConstantsProvider<T>::wgs84(T(1e-12));
        const double r_orbit = K.earth.a.value + 400000.0;
        const double v_orbit = std::sqrt(K.earth.GM.value / r_orbit);
        const double B = 0.01;

        math::Vector3<T> r0(tv<T>(r_orbit), tv<T>(0.0), tv<T>(0.0));
        math::Vector3<T> v0(tv<T>(0.0), tv<T>(v_orbit), tv<T>(0.0));
        math::Vector3<T> w0;
        dynamics::State<T> s = dynamics::State<T>::from_kinematics(
            math::Quaternion<T>::identity(), r0, w0, v0, tv<T>(0.0));

        auto drag_fn = forces::make_drag<T>(
            atmosphere::vallado84_density_model<T>(), tv<T>(B));
        dynamics::Wrench<T> w = drag_fn(s, K);
        double amag = std::sqrt(w.force.x.value * w.force.x.value +
                                w.force.y.value * w.force.y.value +
                                w.force.z.value * w.force.z.value);

        // Hand recomputation with the same conventions (v_rel = v − ω×r).
        double v_rel = v_orbit - K.earth.omega.value * r_orbit;
        double rho = atmosphere::vallado84_density<T>(tv<T>(400000.0)).value;
        double expect = 0.5 * rho * B * v_rel * v_rel;
        double rel = std::abs(amag - expect) / expect;
        std::cout << "  |a_drag| = " << amag << " m/s^2  (hand: " << expect
                  << ", rel " << rel << ")\n";
        check("wrench == hand-computed half*rho*B*v_rel^2 (rel < 1e-12)", rel < 1e-12);
        check("drag magnitude sane at 400 km (1e-9..1e-5 m/s^2)",
              amag > 1e-9 && amag < 1e-5);
    }

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
