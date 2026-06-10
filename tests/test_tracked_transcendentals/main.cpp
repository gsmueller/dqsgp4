/// test_tracked_transcendentals — B3 (DQSGP4 Completion Roadmap, register B3).
///
/// Verifies the TrackedValue<T> transcendental primitives added in B3.1
///   asin, acos, tan, sinh, cosh, tanh, expm1, log1p, log10, log2, hypot,
///   min, max, clamp
/// and the T-dependent small-argument Taylor threshold added in B3.2.
///
/// For each function we check four things:
///   (1) value      — matches the std:: result;
///   (2) UPPER BOUND — the reported error bound is a genuine rigorous upper
///                     bound: sampling the true function across the whole
///                     input interval [x-d, x+d], the largest deviation from
///                     f(x) never exceeds the bound (this is the whole point of
///                     the three-error budget — bounds, never estimates);
///   (3) degenerate  — the documented guards fire (domain edge, pole, d>=x);
///   (4) scaling     — the precision budget TIGHTENS with a wider numeric type
///                     T (the calling card of DQSGP4): cpp_bin_float_50 shrinks
///                     the tracked precision by many orders vs double.
///
/// ExeGate B3: returns a nonzero exit code on any failed check.

#include "math/tracked_value.h"
#include "math/small_angle_series.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>

namespace {

using boost::multiprecision::cpp_bin_float_50;

int passed = 0;
int failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

bool close_rel(double a, double b, double tol) {
    double d = std::abs(a - b);
    double s = std::max(std::abs(a), std::abs(b));
    return d <= tol * std::max(s, 1.0);
}

// (1)+(2): value matches std::f, and the measurement error bound is a genuine
// upper bound on the true excursion of f over the input interval [x-d, x+d].
// We sample the interval (not just endpoints) so the test is valid even when
// the derivative is non-monotone over the interval.
void check_unary(const std::string& name,
                 const std::function<double(double)>& f,
                 const std::function<math::TrackedValue<double>(const math::TrackedValue<double>&)>& tvf,
                 double x, double d) {
    using TV = math::TrackedValue<double>;
    TV a(x, d, TV::representation_bound(x), 0.0);  // error d lives in measurement
    TV r = tvf(a);

    check(name + " value", close_rel(r.value, f(x), 1e-13));

    const double f0 = f(x);
    double worst = 0.0;
    const int N = 41;
    for (int i = 0; i <= N; ++i) {
        double t = x - d + (2.0 * d) * (double(i) / N);
        worst = std::max(worst, std::abs(f(t) - f0));
    }
    const double bound = r.errors.measurement;
    // bound must cover the observed worst-case excursion (small slack for the
    // finite sample missing the exact extremum, and fp rounding).
    check(name + " bound is an upper bound", bound + 1e-15 >= worst * (1.0 - 1e-9));
}

// (4): precision tightens with a wider T. Returns the precision error of f(x)
// when the only input error is the representation bound of x in type T.
template<typename T>
double precision_at(const std::function<math::TrackedValue<T>(const math::TrackedValue<T>&)>& tvf,
                    double x_d) {
    using TV = math::TrackedValue<T>;
    T x = T(x_d);
    TV a(x, T(0), TV::representation_bound(x), T(0));
    return static_cast<double>(tvf(a).errors.precision);
}

void test_asin_acos() {
    using TV = math::TrackedValue<double>;
    check_unary("asin", [](double v){ return std::asin(v); },
                [](const TV& v){ return asin(v); }, 0.5, 1e-3);
    check_unary("acos", [](double v){ return std::acos(v); },
                [](const TV& v){ return acos(v); }, 0.5, 1e-3);
    check_unary("asin near 0", [](double v){ return std::asin(v); },
                [](const TV& v){ return asin(v); }, 0.01, 1e-4);

    // degenerate: |x|+d reaches the domain edge -> capped at pi.
    TV near_edge(0.999, 1e-2, TV::representation_bound(0.999), 0.0);
    TV r = asin(near_edge);
    const double pi = boost::math::constants::pi<double>();
    check("asin domain-edge guard = pi", close_rel(r.errors.measurement, pi, 1e-12));

    // scaling
    double pd = precision_at<double>([](const math::TrackedValue<double>& v){ return asin(v); }, 0.5);
    double pb = precision_at<cpp_bin_float_50>([](const math::TrackedValue<cpp_bin_float_50>& v){ return asin(v); }, 0.5);
    check("asin precision tightens with wider T", pb < pd * 1e-20 && pb > 0.0);
}

void test_tan() {
    using TV = math::TrackedValue<double>;
    check_unary("tan", [](double v){ return std::tan(v); },
                [](const TV& v){ return tan(v); }, 0.5, 1e-3);
    check_unary("tan negative", [](double v){ return std::tan(v); },
                [](const TV& v){ return tan(v); }, -1.0, 1e-3);

    // degenerate: an interval straddling pi/2 contains a pole -> +inf sentinel.
    const double pi = boost::math::constants::pi<double>();
    TV across_pole(pi / 2.0 - 1e-4, 1e-3, TV::representation_bound(pi / 2.0), 0.0);
    TV r = tan(across_pole);
    check("tan pole guard = numeric max",
          r.errors.measurement == std::numeric_limits<double>::max());

    double pd = precision_at<double>([](const math::TrackedValue<double>& v){ return tan(v); }, 0.5);
    double pb = precision_at<cpp_bin_float_50>([](const math::TrackedValue<cpp_bin_float_50>& v){ return tan(v); }, 0.5);
    check("tan precision tightens with wider T", pb < pd * 1e-20 && pb > 0.0);
}

void test_hyperbolics() {
    using TV = math::TrackedValue<double>;
    check_unary("sinh", [](double v){ return std::sinh(v); },
                [](const TV& v){ return sinh(v); }, 0.5, 1e-3);
    check_unary("cosh", [](double v){ return std::cosh(v); },
                [](const TV& v){ return cosh(v); }, 0.5, 1e-3);
    check_unary("cosh straddling 0", [](double v){ return std::cosh(v); },
                [](const TV& v){ return cosh(v); }, 0.0, 1e-3);
    check_unary("tanh", [](double v){ return std::tanh(v); },
                [](const TV& v){ return tanh(v); }, 0.5, 1e-3);
    check_unary("tanh straddling 0", [](double v){ return std::tanh(v); },
                [](const TV& v){ return tanh(v); }, 0.0005, 1e-3);

    double pd = precision_at<double>([](const math::TrackedValue<double>& v){ return sinh(v); }, 0.5);
    double pb = precision_at<cpp_bin_float_50>([](const math::TrackedValue<cpp_bin_float_50>& v){ return sinh(v); }, 0.5);
    check("sinh precision tightens with wider T", pb < pd * 1e-20 && pb > 0.0);
}

void test_expm1_logs() {
    using TV = math::TrackedValue<double>;
    check_unary("expm1", [](double v){ return std::expm1(v); },
                [](const TV& v){ return expm1(v); }, 0.5, 1e-3);
    check_unary("log1p", [](double v){ return std::log1p(v); },
                [](const TV& v){ return log1p(v); }, 0.5, 1e-3);
    check_unary("log10", [](double v){ return std::log10(v); },
                [](const TV& v){ return log10(v); }, 2.0, 1e-3);
    check_unary("log2", [](double v){ return std::log2(v); },
                [](const TV& v){ return log2(v); }, 2.0, 1e-3);

    // expm1 must NOT lose the precision budget for representation-scale d
    // (the catastrophic-cancellation trap exp/expm1 was built to avoid).
    TV a(0.5, 0.0, TV::representation_bound(0.5), 0.0);
    check("expm1 precision > 0 (no cancellation collapse)", expm1(a).errors.precision > 0.0);

    // log1p degenerate: d >= 1+x -> finite flag, value still finite.
    TV degen(-0.5, 1.0, TV::representation_bound(0.5), 0.0);  // 1+x = 0.5 <= d = 1
    TV rd = log1p(degen);
    check("log1p degenerate finite", std::isfinite(rd.errors.measurement) && rd.errors.measurement > 0.0);

    double pd = precision_at<double>([](const math::TrackedValue<double>& v){ return log10(v); }, 2.0);
    double pb = precision_at<cpp_bin_float_50>([](const math::TrackedValue<cpp_bin_float_50>& v){ return log10(v); }, 2.0);
    check("log10 precision tightens with wider T", pb < pd * 1e-20 && pb > 0.0);
}

void test_hypot() {
    using TV = math::TrackedValue<double>;
    TV x(3.0, 1e-3, TV::representation_bound(3.0), 0.0);
    TV y(4.0, 2e-3, TV::representation_bound(4.0), 0.0);
    TV r = hypot(x, y);
    check("hypot value", close_rel(r.value, 5.0, 1e-13));
    // bound = (|x|/r)dx + (|y|/r)dy = 0.6*1e-3 + 0.8*2e-3 = 2.2e-3, and it must
    // be a true upper bound on the gradient-driven excursion.
    double expect = 0.6 * 1e-3 + 0.8 * 2e-3;
    check("hypot measurement bound", close_rel(r.errors.measurement, expect, 1e-9));
    // r = 0 case: Lipschitz-1, bound = dx + dy.
    TV zx(0.0, 1e-3, 0.0, 0.0), zy(0.0, 2e-3, 0.0, 0.0);
    TV rz = hypot(zx, zy);
    check("hypot at origin = dx+dy", close_rel(rz.errors.measurement, 3e-3, 1e-12));
}

void test_min_max_clamp() {
    using TV = math::TrackedValue<double>;
    // well-separated: min/max pick the right value and carry that value's error.
    TV a(1.0, 0.01, 1e-16, 0.0), b(2.0, 0.5, 1e-16, 0.0);
    TV mn = min(a, b), mx = max(a, b);
    check("min value", close_rel(mn.value, 1.0, 1e-15));
    check("max value", close_rel(mx.value, 2.0, 1e-15));
    check("min separated carries a's error", close_rel(mn.errors.measurement, 0.01, 1e-15));
    check("max separated carries b's error", close_rel(mx.errors.measurement, 0.5, 1e-15));

    // overlapping error bars: bound is the per-category max (rigorous fallback).
    TV c(1.0, 0.4, 1e-16, 0.0), e(1.1, 0.4, 1e-16, 0.0);  // gap 0.1 < 0.8 total
    TV mn2 = min(c, e);
    check("min overlapping uses per-category max", close_rel(mn2.errors.measurement, 0.4, 1e-15));

    // clamp forces into range; value at the boundary, errors rigorous.
    TV over(1.5, 1e-3, 1e-16, 0.0);
    TV cl = clamp(over, TV::exact_integer(-1), TV::exact_integer(1));
    check("clamp clips above hi", close_rel(cl.value, 1.0, 1e-15));
    TV inside(0.3, 1e-3, 1e-16, 0.0);
    TV cl2 = clamp(inside, TV::exact_integer(-1), TV::exact_integer(1));
    check("clamp passes through interior value", close_rel(cl2.value, 0.3, 1e-15));
}

// B3.2 — the T-dependent Taylor threshold shrinks with a wider T, and the
// small-angle helpers stay accurate across the (old-threshold, new-threshold)
// band for both double and cpp_bin_float_50.
void test_taylor_threshold() {
    double th_d  = static_cast<double>(math::taylor_branch_threshold<double>(5040.0));
    double th_bf = static_cast<double>(math::taylor_branch_threshold<cpp_bin_float_50>(cpp_bin_float_50(5040)));
    check("threshold(double) ~ 1e-2", th_d > 1e-3 && th_d < 1e-1);
    check("threshold shrinks for wider T", th_bf < th_d * 1e-4);

    // sinc accuracy at a mid-band angle (1e-3): Taylor for double, closed form
    // for bf50 -- both must match the true value.
    using TVd = math::TrackedValue<double>;
    TVd theta(1e-3, 0.0, TVd::representation_bound(1e-3), 0.0);
    TVd theta_sq = theta * theta;
    TVd sinc = math::taylor_sinc(theta, theta_sq);
    double truth = std::sin(1e-3) / 1e-3;
    check("taylor_sinc value accurate at 1e-3", close_rel(sinc.value, truth, 1e-14));

    using TVb = math::TrackedValue<cpp_bin_float_50>;
    TVb thb(cpp_bin_float_50("1e-3"), cpp_bin_float_50(0),
            TVb::representation_bound(cpp_bin_float_50("1e-3")), cpp_bin_float_50(0));
    TVb thb_sq = thb * thb;
    TVb sincb = math::taylor_sinc(thb, thb_sq);
    cpp_bin_float_50 truthb = sin(cpp_bin_float_50("1e-3")) / cpp_bin_float_50("1e-3");
    check("taylor_sinc value accurate at 1e-3 (bf50)",
          close_rel(static_cast<double>(sincb.value), static_cast<double>(truthb), 1e-14));
    // The bf50 sinc precision budget must be far tighter than double's.
    check("taylor_sinc precision tightens with wider T",
          static_cast<double>(sincb.errors.precision) < static_cast<double>(sinc.errors.precision) * 1e-10);
}

// B3.3 — Vector3 normalize() and scalar division.
void test_vector3() {
    using TV = math::TrackedValue<double>;
    using V3 = math::Vector3<double>;
    V3 v(TV(3.0, 0.0, TV::representation_bound(3.0), 0.0),
         TV(0.0, 0.0, 0.0, 0.0),
         TV(4.0, 0.0, TV::representation_bound(4.0), 0.0));
    V3 u = v.normalize();
    check("normalize x = 0.6", close_rel(u.x.value, 0.6, 1e-13));
    check("normalize z = 0.8", close_rel(u.z.value, 0.8, 1e-13));
    double mag = std::sqrt(u.x.value * u.x.value + u.y.value * u.y.value + u.z.value * u.z.value);
    check("normalize unit magnitude", close_rel(mag, 1.0, 1e-13));

    V3 w = v / TV::exact_integer(2);
    check("scalar division", close_rel(w.x.value, 1.5, 1e-13) && close_rel(w.z.value, 2.0, 1e-13));
}

} // namespace

int main() {
    test_asin_acos();
    test_tan();
    test_hyperbolics();
    test_expm1_logs();
    test_hypot();
    test_min_max_clamp();
    test_taylor_threshold();
    test_vector3();

    std::cout << "\n  tracked transcendentals: " << passed << " passed, "
              << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
