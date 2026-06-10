/// test_time_scales — L1 (time scales & epochs).
///
/// Verifies astronomy::Epoch / TimeScale and the derived views against:
///   - the EXACT convention bases (oracle-anchored vs astropy/ERFA 2026-06-06:
///     JD 2451545.0 TT = 2000-01-01 12:00 = J2000.0; JD 2415020.0 TT =
///     1899-12-31 12:00 = the J1900 base, exactly 36525 d = 100 Julian yr
///     before J2000);
///   - the in-repo offsets the views must reproduce bit-for-bit, so the
///     pending migration of deep_space.h / sidereal_time.h is value-preserving
///     (OR1-safe);
///   - precision tightening with a wider numeric type T (the calling card).
///
/// Theory: design/derivations/time_scales_and_epochs.md. ExeGate TIME1.

#include "astronomy/epoch.h"
#include "math/tracked_value.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_bin_float_50;
int passed = 0, failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

template<typename T>
math::TrackedValue<T> tv(double v) {
    return math::TrackedValue<T>(static_cast<T>(v), T(0),
                                 math::TrackedValue<T>::representation_bound(static_cast<T>(v)), T(0));
}

} // namespace

int main() {
    using T = double;
    using astronomy::Epoch;
    using astronomy::TimeScale;
    using astronomy::days_since_1900;
    using astronomy::centuries_since_j2000;

    // 1. Exact convention bases (oracle-anchored). Values + exactness.
    check("J2000 base == JD 2451545",
          astronomy::jd_epoch_j2000<T>().value == 2451545.0);
    check("J1900 base == JD 2415020",
          astronomy::jd_epoch_j1900<T>().value == 2415020.0);
    check("Julian century == 36525 d",
          astronomy::days_per_century<T>().value == 36525.0);
    // The key identity: the J1900 base is exactly one Julian century before J2000.
    check("J2000 - J1900 == 1 Julian century (exact)",
          (astronomy::jd_epoch_j2000<T>().value - astronomy::jd_epoch_j1900<T>().value) == 36525.0);
    // These bases are exact-by-convention -> zero accuracy AND zero precision.
    check("convention bases are EXACT (precision 0, accuracy 0)",
          astronomy::jd_epoch_j2000<T>().errors.precision == T(0) &&
          astronomy::jd_epoch_j2000<T>().errors.accuracy == T(0) &&
          astronomy::days_per_century<T>().errors.precision == T(0));

    // 2. The views reproduce the in-repo offsets BIT-FOR-BIT (migration-ready).
    //    Sample epoch: a 2024-era TLE JD.
    const double JD = 2460310.5;
    Epoch<T> e = Epoch<T>::from_jd(tv<T>(JD));
    check("days_since_1900 == JD - 2415020 (deep_space.h arithmetic)",
          days_since_1900(e).value == (JD - 2415020.0));
    check("centuries_since_j2000 == (JD - 2451545)/36525 (sidereal_time.h arithmetic)",
          centuries_since_j2000(e).value == ((JD - 2451545.0) / 36525.0));

    // 3. Epoch construction round-trips; the scale tag is preserved.
    check("from_jd round-trips the JD value", e.jd().value == JD);
    Epoch<T> e2 = Epoch<T>::from_jd_two_part(tv<T>(2460310.0), tv<T>(0.5), TimeScale::TT);
    check("from_jd_two_part: jd() == day + frac", e2.jd().value == JD);
    check("scale tag preserved", e2.scale == TimeScale::TT && e.scale == TimeScale::UTC);

    // 4. Precision of the J2000-century view tightens with a wider T (the
    //    calling card): the epoch's representation precision propagates through.
    auto cprec = [JD](auto dummy) {
        using U = decltype(dummy);
        Epoch<U> ep = Epoch<U>::from_jd(
            math::TrackedValue<U>(U(JD), U(0),
                                  math::TrackedValue<U>::representation_bound(U(JD)), U(0)));
        return static_cast<double>(centuries_since_j2000(ep).errors.precision);
    };
    double p_d = cprec(double(0));
    double p_b = cprec(cpp_bin_float_50(0));
    check("century-view precision > 0 (double)", p_d > 0.0);
    check("century-view precision tightens with wider T", p_b < p_d * 1e-20);

    std::cout << "\n  time scales: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
