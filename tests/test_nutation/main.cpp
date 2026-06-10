/// test_nutation (gate NUT1, R4b) — IAU 2000A luni-solar nutation + the
/// Earth-fixed chain GCRS→ITRS, erfa-arbitrated at every layer.
///
/// Verifies design/derivations/nutation_itrs.md §6 against the embedded
/// gen_itrs_oracle.py reference (erfa nut06a/gmst06/gst06a/pom00/c2t06a at five
/// epochs, with REAL IERS x_p/y_p/ΔUT1/ΔAT parsed from the in-repo CSSI EOP
/// file). Measured calibration (recorded in the note):
///   Δψ residual vs nut06a  0.479 mas  (the real planetary signal; the tracked
///                                      5.54 mas floor majorizes it 11×)
///   Δε residual            0.375 mas  (floor 1.90 mas)
///   GMST                   2.41 μas   (the AstroLib-form arbitration bound)
///   GAST                   0.436 mas
///   polar motion           EXACTLY 0  (bit-for-bit vs erfa.pom00)
///   full chain elements    1.85e-9    (≈ 0.4 mas, the honest mas grade)
///   dial-up n=20 deficit   1.06 mas   (tail bound 61.8 mas majorizes)
///
/// Exit 0 iff every check passes.

#include "astronomy/earth_rotation.h"
#include "astronomy/epoch.h"
#include "math/matrix3.h"
#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using boost::multiprecision::cpp_bin_float_50;
using T = double;
using TV = math::TrackedValue<T>;

namespace {

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cerr << "  FAIL: " << name << "\n"; }
}

TV tv(double v) { return TV(v, T(0), TV::representation_bound(v), T(0)); }

const double PI_D = 3.14159265358979323846;
const double AS2R = PI_D / (180.0 * 3600.0);

double ang_diff(double a, double b) {
    double d = std::fmod(a - b + 3 * PI_D, 2 * PI_D) - PI_D;
    return std::abs(d);
}

// gen_itrs_oracle.py output, embedded VERBATIM (erfa 2.0.1.5; EOP from the
// in-repo EOP-All-v1.1_2025-01-10.txt).
struct Row {
    double jd_tt, jd_ut1, xp_as, yp_as;
    double dpsi, deps, gmst, gast;
    double c2t[9];
    double pom[9];
};
const Row ORACLE[] = {
  {2451544.5007428704, 2451544.5000041141, 0.0432610, 0.3779910,
   -6.75427320233025274e-05, -2.79231974583952498e-05, 1.74479315466571738e+00, 1.74473119562538259e+00,
   {-1.73120441988564217e-01, 9.84900660997882649e-01, 2.30568391158817568e-05, -9.84900660845308584e-01, -1.73120441181193629e-01, -3.33422117460838426e-05, -2.88471562278640908e-05, -2.84809145165942360e-05, 9.99999999178339483e-01},
   {9.99999999999978018e-01, 3.11463224668065305e-15, 2.09735246584794840e-07, 3.81236130417203947e-13, 9.99999999998320899e-01, -1.83255208136168028e-06, -2.09735246584448403e-07, 1.83255208136171987e-06, 9.99999999998298916e-01}},
  {2455454.5007660184, 2455454.4999993467, 0.2249280, 0.3970510,
   8.35159973215753840e-05, 8.75642446574797483e-06, 6.17582593906126665e+00, 6.17590255153390277e+00,
   {9.93982623354567218e-01, -1.09532692176054186e-01, -1.06480694898606568e-03, 1.09532619172818868e-01, 9.93983193648654728e-01, -1.26811526013030282e-04, 1.07229021961561069e-03, 9.41735926216587919e-06, 9.99999425052333857e-01},
   {9.99999999999405476e-01, -2.43895508584864951e-11, 1.09048171664584074e-06, 2.64886818927486525e-11, 9.99999999998147260e-01, -1.92495756897989051e-06, -1.09048171659687144e-06, 1.92495756900763126e-06, 9.99999999997552735e-01}},
  {2458849.5008007409, 2458849.4999979497, 0.0766140, 0.2823090,
   -7.99657379992177106e-05, -8.25126074740569501e-06, 1.74744231300448427e+00, 1.74736895702739381e+00,
   {-1.71324162932526891e-01, 9.85214654533156931e-01, 3.40130241655677631e-04, -9.85212862035159453e-01, -1.71324500061000362e-01, 1.87940398493412837e-03, 1.90988899135240164e-03, -1.31133743154055769e-05, 9.99998176074376799e-01},
   {9.99999999999931055e-01, -4.55693717530737433e-11, 3.71435153645251372e-07, 4.60777448909353327e-11, 9.99999999999063416e-01, -1.36867265500299800e-06, -3.71435153582533929e-07, 1.36867265502001849e-06, 9.99999999998994471e-01}},
  {2460462.5008007409, 2460462.4999997634, 0.0339010, 0.4507450,
   -2.14508276354261212e-05, 4.15644255846936160e-05, 4.36281576685096972e+00, 4.36279608863737689e+00,
   {-3.47620181727700206e-01, -9.37635045564759517e-01, 8.54742284030599302e-04, 9.37632450335704748e-01, -3.47621230854622598e-01, -2.20634008693803077e-03, 2.36586835278570443e-03, 3.44657602066395102e-05, 9.99997200735606273e-01},
   {9.99999999999986455e-01, -5.56321260251076253e-11, 1.64356686032943063e-07, 5.59912903234992225e-11, 9.99999999997612243e-01, -2.18527342691540885e-06, -1.64356685910979213e-07, 2.18527342692458222e-06, 9.99999999997598699e-01}},
  {2449796.5007081483, 2449796.5000021891, -0.0087670, 0.5498250,
   5.18069305739313332e-05, -3.25770998575075201e-05, 3.09022756579974622e+00, 3.09027508999674172e+00,
   {-9.98735387114118023e-01, 5.02735595600260765e-02, -4.42419684196408832e-04, -5.02735388536778571e-02, -9.98735484474997603e-01, -5.78067422761140840e-05, -4.44766388337965114e-04, -3.54916359418608131e-05, 9.99999900461596813e-01},
   {9.99999999999999112e-01, 1.09080711213417403e-11, -4.25036154228730075e-08, -1.10213698986128668e-11, 9.99999999996447175e-01, -2.66562682215734700e-06, 4.25036153936451497e-08, 2.66562682215781329e-06, 9.99999999996446287e-01}},
};
const int N_ORACLE = 5;

} // namespace

int main() {
    std::cout << std::setprecision(4);
    std::cout << "test_nutation (NUT1, R4b): IAU 2000A nutation + GCRS->ITRS chain\n\n";

    double w_psi = 0, w_eps = 0, w_gmst = 0, w_gast = 0, w_pom = 0, w_c2t = 0;
    bool psi_majorized = true, eps_majorized = true, gast_majorized = true;
    bool ortho_ok = true;

    for (int k = 0; k < N_ORACLE; ++k) {
        const Row& r = ORACLE[k];
        TV t_tt = astronomy::centuries_since_j2000<T>(tv(r.jd_tt));

        astronomy::NutationAngles<T> ang = astronomy::nutation_angles<T>(t_tt);
        double dp = std::abs(ang.dpsi.value - r.dpsi);
        double de = std::abs(ang.deps.value - r.deps);
        w_psi = std::max(w_psi, dp);
        w_eps = std::max(w_eps, de);
        if (ang.dpsi.errors.accuracy <= dp) psi_majorized = false;
        if (ang.deps.errors.accuracy <= de) eps_majorized = false;

        double gmst = astronomy::gmst_iau2006<T>(tv(r.jd_ut1), t_tt).value;
        TV gast_tv = astronomy::gast_iau2006<T>(tv(r.jd_ut1), t_tt);
        double dgast = ang_diff(gast_tv.value, r.gast);
        w_gmst = std::max(w_gmst, ang_diff(gmst, r.gmst));
        w_gast = std::max(w_gast, dgast);
        if (gast_tv.errors.accuracy <= dgast) gast_majorized = false;

        math::Matrix3<T> W = astronomy::polar_motion<T>(
            tv(r.xp_as * AS2R), tv(r.yp_as * AS2R), t_tt);
        math::Matrix3<T> M = astronomy::gcrs_to_itrs<T>(
            t_tt, tv(r.jd_ut1), tv(r.xp_as * AS2R), tv(r.yp_as * AS2R));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                w_pom = std::max(w_pom, std::abs(W.m[i][j].value - r.pom[3 * i + j]));
                w_c2t = std::max(w_c2t, std::abs(M.m[i][j].value - r.c2t[3 * i + j]));
            }
        }
        // Chain orthonormality: M·Mᵀ = I to round-off.
        math::Matrix3<T> MMt = M * M.transpose();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double want = (i == j) ? 1.0 : 0.0;
                if (std::abs(MMt.m[i][j].value - want) > 1e-12) ortho_ok = false;
            }
        }
    }

    std::cout << "=== angles vs erfa.nut06a (the COMPLETE model) ===\n";
    std::cout << "  worst dpsi " << w_psi / AS2R * 1e3 << " mas, deps "
              << w_eps / AS2R * 1e3 << " mas (floors 5.54 / 1.90 mas)\n";
    check("dpsi residual < 1 mas (measured 0.48; the planetary signal)",
          w_psi < 1e-3 * AS2R);
    check("deps residual < 1 mas (measured 0.38)", w_eps < 1e-3 * AS2R);
    check("tracked dpsi accuracy MAJORIZES the residual (all epochs)", psi_majorized);
    check("tracked deps accuracy MAJORIZES the residual (all epochs)", eps_majorized);

    std::cout << "\n=== rotation angles ===\n";
    std::cout << "  worst GMST " << w_gmst / AS2R * 1e6 << " uas, GAST "
              << w_gast / AS2R * 1e3 << " mas\n";
    check("GMST06 vs erfa.gmst06 < 10 uas (measured 2.4)", w_gmst < 1e-5 * AS2R);
    check("GAST vs erfa.gst06a < 1 mas (measured 0.44)", w_gast < 1e-3 * AS2R);
    check("tracked GAST accuracy MAJORIZES the residual", gast_majorized);

    std::cout << "\n=== matrices ===\n";
    std::cout << "  pom element diff " << w_pom << ", chain element diff " << w_c2t << "\n";
    check("polar motion vs erfa.pom00 at machine level (< 1e-12)", w_pom < 1e-12);
    check("full chain vs erfa.c2t06a < 5e-9 elements (measured 1.85e-9)", w_c2t < 5e-9);
    check("chain orthonormal (M*Mt = I to 1e-12)", ortho_ok);

    // Dial-up truncation: n = 20 vs full; the tracked tail majorizes the deficit.
    std::cout << "\n=== dial-up truncation ===\n";
    {
        TV t0 = astronomy::centuries_since_j2000<T>(tv(ORACLE[0].jd_tt));
        astronomy::NutationAngles<T> full = astronomy::nutation_angles<T>(t0);
        astronomy::NutationAngles<T> a20 = astronomy::nutation_angles<T>(t0, 20);
        double deficit = std::abs(a20.dpsi.value - full.dpsi.value);
        std::cout << "  n=20 deficit " << deficit / AS2R * 1e3 << " mas, tail bound "
                  << a20.dpsi.errors.accuracy / AS2R * 1e3 << " mas\n";
        check("n=20 tracked tail MAJORIZES the truncation deficit",
              a20.dpsi.errors.accuracy > deficit);
        check("dialing terms up reduces the truncation bound",
              full.dpsi.errors.accuracy < a20.dpsi.errors.accuracy);
    }

    // O4 (sidereal-ratio provenance, closed as RATIFICATION — no value change):
    // the adopted SGP4 ratio 1.00273790935 IS the Aoki-82 GMST rate, exactly:
    // 1 + c1/(86400 s/d * 36525 d/cy) with c1 = 8640184.812866 sidereal-s/cy —
    // the very coefficient sidereal_time.h carries. The generative identity
    // reproduces the adopted value within its written-digit floor.
    std::cout << "\n=== O4: sidereal-ratio provenance (generative ratification) ===\n";
    {
        double generated = 1.0 + 8640184.812866 / (86400.0 * 36525.0);
        double adopted = 1.00273790935;
        std::cout << std::setprecision(15) << "  generated " << generated
                  << " vs adopted " << adopted << std::setprecision(4) << "\n";
        check("adopted sidereal ratio == 1 + c1_aoki/(86400*36525) to its digit floor",
              std::abs(generated - adopted) < 5e-12);
    }

    // Precision tightens with wider T.
    std::cout << "\n=== precision scaling ===\n";
    {
        TV t0 = astronomy::centuries_since_j2000<T>(tv(ORACLE[0].jd_tt));
        astronomy::NutationAngles<T> ad = astronomy::nutation_angles<T>(t0, 60);
        math::TrackedValue<cpp_bin_float_50> tb =
            astronomy::centuries_since_j2000<cpp_bin_float_50>(
                tv(ORACLE[0].jd_tt).value == 0
                    ? math::TrackedValue<cpp_bin_float_50>()
                    : math::TrackedValue<cpp_bin_float_50>(
                          cpp_bin_float_50(ORACLE[0].jd_tt), cpp_bin_float_50(0),
                          math::TrackedValue<cpp_bin_float_50>::representation_bound(
                              cpp_bin_float_50(ORACLE[0].jd_tt)),
                          cpp_bin_float_50(0)));
        astronomy::NutationAngles<cpp_bin_float_50> ab =
            astronomy::nutation_angles<cpp_bin_float_50>(tb, 60);
        double pd = static_cast<double>(ad.dpsi.errors.precision);
        double pb = static_cast<double>(ab.dpsi.errors.precision);
        std::cout << "  dpsi precision: double = " << pd << "  bf50 = " << pb << "\n";
        check("precision > 0 (framework alive)", pd > 0.0 && pb > 0.0);
        check("precision tightens with wider T", pb < pd);
    }

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
