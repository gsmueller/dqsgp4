/// test_frame_chain — L3 ecliptic-of-date → GCRS frame chain (frame_chain.md §6).
///
/// Element-wise CONFORMANCE oracle: the IAU 2006 rotations must reproduce the ERFA
/// (IAU SOFA) reference to ~1e-12 per element at the four §A epochs:
///   - precession_iau2006(t)            == erfa.pmat06   (bias-precession matrix);
///   - ecliptic_to_gcrs(t)              == erfa.ecm06ᵀ   (composite ecliptic→GCRS);
///   - obliquity_iau2006(t)             == erfa.obl06     (mean obliquity ε_A).
/// These ARE the IAU analytical theory, so matching them proves the rotations
/// implement IAU 2006 — a conformance check, NOT a fidelity claim (the body-position
/// fidelity vs JPL DE430 lands with body_position_gcrs). Reference values are
/// born-digital from tools/gen_frame_oracle.py (erfa 2.0.1.5), embedded below.
///
/// Plus: the new Matrix3·Matrix3 product is exercised by an R·Rᵀ = I identity, and
/// the precision channel tightens with a wider numeric type T (the calling card).
///
/// ExeGate: nonzero exit code on any failed check.

#include "astronomy/epoch.h"
#include "astronomy/frames.h"
#include "astronomy/obliquity.h"
#include "ephemeris/body_position_gcrs.h"
#include "ephemeris/moon_meeus.h"
#include "math/matrix3.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

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

template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0),
        math::TrackedValue<T>::representation_bound(val), T(0));
}

// erfa reference (tools/gen_frame_oracle.py): JD, obl06 [rad], pmat06[9] (row-major
// GCRS→date), ecm06[9] (row-major GCRS→ecliptic). Composite ecliptic→GCRS = ecm06ᵀ.
struct Ref { double jd; double obl; double P[9]; double M[9]; };
const Ref kRefs[] = {
  { 2415020.0, 4.093196610614512898e-01,
    {9.997029465818912941e-01, 2.235135417607689995e-02, 9.717796148294214481e-03,
     -2.235135345466072662e-02, 9.997501713850852623e-01, -1.086933176862592276e-04,
     -9.717797807582641625e-03, -1.085448665461563600e-04, 9.999527751968004807e-01},
    {9.997029465818912941e-01, 2.235135417607689995e-02, 9.717796148294214481e-03,
     -2.437248888645331818e-02, 9.171194100108934766e-01, 3.978667736400593458e-01,
     -1.951829760774188849e-05, -3.979854328336227542e-01, 9.173917346865896327e-01} },
  { 2451545.0, 4.090926006005828897e-01,
    {9.999999999999941158e-01, -7.078368960971556128e-08, 8.056213977613186084e-08,
     7.078368694637676268e-08, 9.999999999999968914e-01, 3.305943735432137487e-08,
     -8.056214211620057482e-08, -3.305943169218394928e-08, 9.999999999999962252e-01},
    {9.999999999999941158e-01, -7.078368960971556128e-08, 8.056213977613186084e-08,
     3.289700407741964649e-08, 9.174821299149583664e-01, 3.977769994440479295e-01,
     -1.020704472548435540e-07, -3.977769994440430446e-01, 9.174821299149555909e-01} },
  { 2459000.5, 4.090462507949458870e-01,
    {9.999876158842090357e-01, -4.564545070165593535e-03, -1.983180909044943239e-03,
     4.564545186470075205e-03, 9.999895823993776300e-01, -4.467542036312721621e-06,
     1.983180641355248609e-03, -4.584832162479735018e-06, 9.999980334848279639e-01},
    {9.999876158842090357e-01, -4.564545070165593535e-03, -1.983180909044943239e-03,
     4.976752100178826257e-03, 9.174891972646534999e-01, 3.977295624431190779e-01,
     4.092546761373183975e-06, -3.977345067139183565e-01, 9.175005515814620427e-01} },
  { 2469807.5, 4.089790660606221206e-01,
    {9.999256843098003333e-01, -1.118167244077446741e-02, -4.857577483145684126e-03,
     1.118167289642968930e-02, 9.999374827751563721e-01, -2.706512742406408378e-05,
     4.857576434271398090e-03, -2.725272642511988863e-05, 9.999882015346350794e-01},
    {9.999256843098003333e-01, -1.118167244077446741e-02, -4.857577483145684126e-03,
     1.219121615115483989e-02, 9.174590996486576833e-01, 3.976432757140878627e-01,
     1.031180684893211186e-05, -3.976729443566680011e-01, 9.175272362280941785e-01} },
};

double t_jcen(double jd) { return (jd - 2451545.0) / 36525.0; }

} // namespace

int main() {
    using T = double;
    const double tol = 1e-12;   // element-wise vs erfa (round-off of identical theory)

    double worst_P = 0.0, worst_M = 0.0, worst_obl = 0.0;
    for (const Ref& r : kRefs) {
        math::TrackedValue<T> t = tv<T>(t_jcen(r.jd));

        // precession_iau2006 == erfa.pmat06 (element-wise).
        math::Matrix3<T> P = astronomy::precession_iau2006<T>(t);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double d = std::abs(P.m[i][j].value - r.P[i * 3 + j]);
                if (d > worst_P) worst_P = d;
            }

        // ecliptic_to_gcrs == erfa.ecm06ᵀ : my M[i][j] vs ecm06[j][i].
        math::Matrix3<T> M = astronomy::ecliptic_to_gcrs<T>(t);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double d = std::abs(M.m[i][j].value - r.M[j * 3 + i]);
                if (d > worst_M) worst_M = d;
            }

        // obliquity_iau2006 == erfa.obl06.
        double od = std::abs(astronomy::obliquity_iau2006<T>(t).value - r.obl);
        if (od > worst_obl) worst_obl = od;
    }
    std::cout << "    worst |Δ| vs erfa:  pmat06=" << worst_P
              << "  ecm06ᵀ=" << worst_M << "  obl06=" << worst_obl << "\n";
    check("precession_iau2006 == erfa.pmat06 (<1e-12 element-wise)", worst_P < tol);
    check("ecliptic_to_gcrs == erfa.ecm06ᵀ (<1e-12 element-wise)", worst_M < tol);
    check("obliquity_iau2006 == erfa.obl06 (<1e-12 rad)", worst_obl < tol);

    // Matrix3·Matrix3 (new) exercised by R·Rᵀ = I.
    math::Matrix3<T> R = astronomy::rot_x<T>(tv<T>(0.7))
                       * astronomy::rot_z<T>(tv<T>(-1.3));
    math::Matrix3<T> I = R * R.transpose();
    double worst_I = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double want = (i == j) ? 1.0 : 0.0;
            double d = std::abs(I.m[i][j].value - want);
            if (d > worst_I) worst_I = d;
        }
    check("Matrix3 product: R·Rᵀ = I (<1e-14)", worst_I < 1e-14);

    // Precision channel tightens with a wider T (the calling card), on the
    // precession off-diagonal at ~2020 (a nonzero element).
    math::TrackedValue<T> td = tv<T>(t_jcen(2459000.5));
    double pd = astronomy::precession_iau2006<double>(td).m[0][1].errors.precision;
    math::TrackedValue<cpp_bin_float_50> tb = tv<cpp_bin_float_50>(t_jcen(2459000.5));
    double pb = static_cast<double>(
        astronomy::precession_iau2006<cpp_bin_float_50>(tb).m[0][1].errors.precision);
    check("precession precision tracked (> 0)", pd > 0.0 && pb > 0.0);
    check("precession precision tightens with wider T", pb < pd * 1e-20);

    // ---- end-to-end FIDELITY: Meeus Moon → GCRS vs JPL DE430 (independent) ----
    // DE430 (in-repo) gives the Moon's GCRS (J2000 mean-equatorial) position
    // directly (km); we form moon_meeus → ecliptic Cartesian → ecliptic_to_gcrs →
    // GCRS and compare the DIRECTION. The chain is mean-CONSISTENT (Meeus is
    // mean-ecliptic-referenced, ecm06ᵀ is the mean-ecliptic→GCRS rotation), so the
    // angular residual is PURELY the Meeus truncation — the orthonormal rotation
    // adds ≤1e-12 and nutation does NOT enter (it is a mean↔true equinox frame
    // difference we never visit). DE430 is NUMERICAL — independent of the
    // analytical theory moon_meeus realizes (the no-perceived-fidelity oracle).
    struct MoonRef { double jd, x, y, z; };
    static const MoonRef kMoon[] = {
        { 2451545.000000, -291608.3848, -266716.8335,  -76102.4860 },
        { 2455455.000000,  -31395.8604, -354133.4293, -162301.4360 },
        { 2458849.500000,  390185.6383,  -76522.5985,  -70724.6561 },
        { 2460462.500000,  368305.3099,  -10570.1034,  -14560.1512 },
        { 2449796.500000, -289619.4202, -207427.2619,  -88753.8919 },
    };
    const double pi = 3.14159265358979323846;
    double worst_sep = 0.0;  // arcsec
    for (const MoonRef& mr : kMoon) {
        astronomy::Epoch<T> ep =
            astronomy::Epoch<T>::from_jd(tv<T>(mr.jd), astronomy::TimeScale::TT);
        ephemeris::EclipticState<T> ecl = ephemeris::moon_meeus_of_date<T>(ep);
        math::Vector3<T> g = ephemeris::body_position_gcrs<T>(ecl, ep);
        double gx = g.x.value, gy = g.y.value, gz = g.z.value;
        double gm = std::sqrt(gx * gx + gy * gy + gz * gz);
        double dm = std::sqrt(mr.x * mr.x + mr.y * mr.y + mr.z * mr.z);
        double cosang = (gx * mr.x + gy * mr.y + gz * mr.z) / (gm * dm);
        if (cosang > 1.0) cosang = 1.0;
        double sep = std::acos(cosang) * 648000.0 / pi;  // arcsec
        if (sep > worst_sep) worst_sep = sep;
    }
    std::cout << "    end-to-end Moon→GCRS vs DE430: worst sep = " << worst_sep
              << " arcsec\n";
    // Bound = the Meeus 60-term grade (~10″, born_digital_sources.md §B),
    // majorizing the measured ~3.7″ — and far below the degrees a broken rotation
    // would give, so it catches both a frame error and a Meeus regression.
    check("Meeus Moon → GCRS within Meeus accuracy of DE430 (fidelity)",
          worst_sep < 10.0);

    std::cout << "\n  frame chain: " << passed << " passed, " << failed
              << " failed\n";
    return failed == 0 ? 0 : 1;
}
