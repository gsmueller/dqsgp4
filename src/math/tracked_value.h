#pragma once

/**
 * @file tracked_value.h
 * @brief TrackedValue<T> — the core numeric type with rigorous three-error
 *        propagation (the substrate everything above this layer computes in).
 *
 * Every value carries three independent error bounds:
 *   measurement  — physical measurement uncertainty of inputs
 *   precision    — computational representation and truncation error
 *   accuracy     — model fidelity (truncated perturbation series, etc.)
 *
 * All bounds are RIGOROUS UPPER BOUNDS, never estimates: each operation's
 * bound is derived from a worst-case derivative (or exact convexity) argument
 * recorded beside it. There is no bare T anywhere in the computational path
 * above this layer. Gate W12 (test_error_framework) asserts the propagation;
 * B1/B2 cover the exp and pow paths; the named constructors encode the
 * constants-honesty provenance scheme (CR1B).
 */

#include <cmath>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <boost/math/constants/constants.hpp>

namespace math {

/// The three independent error categories every TrackedValue carries. Each is
/// a non-negative rigorous upper bound; the categories propagate separately
/// so a wider T tightens precision without touching the model floor.
template<typename T>
struct ThreeErrors {
    T measurement;  ///< sigma_m: from physical input uncertainties
    T precision;    ///< delta_p: from representation, truncation, rounding
    T accuracy;     ///< delta_a: from model truncation/simplification

    /// All categories zero (an exact quantity).
    ThreeErrors() : measurement(T(0)), precision(T(0)), accuracy(T(0)) {}
    /// Explicit per-category bounds.
    ThreeErrors(const T& m, const T& p, const T& a) : measurement(m), precision(p), accuracy(a) {}

    /// Conservative upper bound (triangle inequality). Correct because the
    /// three error sources are not necessarily independent.
    T total() const { return measurement + precision + accuracy; }

    /// Statistical estimate (RSS). Use when reporting expected error rather
    /// than worst-case bound. Only valid if sources are independent.
    T rss() const {
        using std::sqrt;
        return sqrt(measurement * measurement + precision * precision + accuracy * accuracy);
    }

    // --- Per-category arithmetic ---
    // Applies the same operation to each error category independently.
    // This is the correct propagation for additive error combination
    // (addition, subtraction) and for scaling by a non-negative factor.

    /// Per-category addition: used for error combination in addition/subtraction.
    friend ThreeErrors operator+(const ThreeErrors& a, const ThreeErrors& b) {
        return ThreeErrors(
            a.measurement + b.measurement,
            a.precision + b.precision,
            a.accuracy + b.accuracy
        );
    }

    /// Per-category subtraction. Result is clamped to zero per-category
    /// because error bounds are non-negative by definition.
    friend ThreeErrors operator-(const ThreeErrors& a, const ThreeErrors& b) {
        using std::max;
        return ThreeErrors(
            max(a.measurement - b.measurement, T(0)),
            max(a.precision - b.precision, T(0)),
            max(a.accuracy - b.accuracy, T(0))
        );
    }

    /// Scale all error categories by a factor. Uses |s| to guarantee
    /// the result is non-negative regardless of the sign of s.
    friend ThreeErrors operator*(const T& s, const ThreeErrors& e) {
        using std::abs;
        T as = abs(s);
        return ThreeErrors(as * e.measurement, as * e.precision, as * e.accuracy);
    }

    friend ThreeErrors operator*(const ThreeErrors& e, const T& s) {
        using std::abs;
        T as = abs(s);
        return ThreeErrors(e.measurement * as, e.precision * as, e.accuracy * as);
    }

    /// Per-category division by a positive scalar.
    friend ThreeErrors operator/(const ThreeErrors& e, const T& s) {
        using std::abs;
        T as = abs(s);
        if (as == T(0)) {
            T inf = std::numeric_limits<T>::max();
            return ThreeErrors(inf, inf, inf);
        }
        return ThreeErrors(e.measurement / as, e.precision / as, e.accuracy / as);
    }

    /// Per-category maximum: used when the bound is the max of two sources.
    static ThreeErrors max_per_category(const ThreeErrors& a, const ThreeErrors& b) {
        using std::max;
        return ThreeErrors(
            max(a.measurement, b.measurement),
            max(a.precision, b.precision),
            max(a.accuracy, b.accuracy)
        );
    }

    /// Apply a scalar function f to each category independently.
    /// f must be non-negative-preserving (e.g., a bound computation).
    template<typename F>
    ThreeErrors apply(F f) const {
        return ThreeErrors(f(measurement), f(precision), f(accuracy));
    }
};

namespace detail {
/// Guards against 0 × overflowed-error → NaN: |x·y| is 0 when either factor is
/// exactly 0, regardless of the other's (possibly overflowed) error magnitude.
template<typename T>
inline T tv_safe_mul(T x, T y) {
    return (x == T(0) || y == T(0)) ? T(0) : x * y;
}
/// Multiplication error bound |x·y| ≤ |a|·δb + |b|·δa + δa·δb (one category).
template<typename T>
inline T tv_mul_bound(T av, T bv, T a_err, T b_err) {
    return tv_safe_mul(av, b_err) + tv_safe_mul(bv, a_err) + tv_safe_mul(a_err, b_err);
}
/// Division error bound (|num|·δden + |den|·δnum) / (|den|·(|den|−δden)), using
/// the minimum possible denominator magnitude for rigor.
template<typename T>
inline T tv_div_bound(T num_abs, T den_abs, T num_err, T den_err) {
    T numerator = tv_safe_mul(num_abs, den_err) + tv_safe_mul(den_abs, num_err);
    T denom_min = den_abs - den_err;
    if (denom_min <= T(0)) {
        return std::numeric_limits<T>::max();
    }
    return numerator / (den_abs * denom_min);
}
}  // namespace detail

/// A numeric value plus its ThreeErrors budget. Construct through the named
/// constructors (exact_integer / defined / measured / from_truncated_decimal /
/// model_coefficient) so the provenance of every constant is explicit; all
/// arithmetic and the friend math functions propagate every category with the
/// worst-case bound derived in the comment beside each operation.
template<typename T>
class TrackedValue {
public:
    T value;                 ///< the computed value
    ThreeErrors<T> errors;   ///< its three-category rigorous budget

    // --- Constructors ---

    /// Exact zero.
    TrackedValue() : value(T(0)), errors() {}

    /// A value with an explicit, already-formed budget.
    TrackedValue(const T& val, const ThreeErrors<T>& err)
        : value(val), errors(err) {}

    /// A value with explicit per-category bounds.
    TrackedValue(const T& val, const T& meas, const T& prec, const T& acc)
        : value(val), errors(meas, prec, acc) {}

    // --- Named constructors ---

    /// Exact integer: zero error in all three categories
    static TrackedValue exact_integer(int n) {
        return TrackedValue(T(n), T(0), T(0), T(0));
    }

    /// Defined decimal constant: precision error from binary representation
    /// measurement = 0 (definitional), accuracy = 0 (no model)
    static TrackedValue defined(const char* decimal_str) {
        T val = from_string(decimal_str);
        T prec = representation_bound(val);
        return TrackedValue(val, T(0), prec, T(0));
    }

    /// Defined with known physical measurement uncertainty
    /// (for the "physical Earth" view of defining parameters)
    static TrackedValue defined_with_physical_uncertainty(
        const char* decimal_str, const char* measurement_sigma)
    {
        T val = from_string(decimal_str);
        T meas = from_string(measurement_sigma);
        T prec = representation_bound(val);
        return TrackedValue(val, meas, prec, T(0));
    }

    /// Measured physical constant
    static TrackedValue measured(const char* value_str, const char* sigma_str) {
        T val = from_string(value_str);
        T meas = from_string(sigma_str);
        T prec = representation_bound(val);
        return TrackedValue(val, meas, prec, T(0));
    }

    /// Place value of the least-significant displayed digit of a written decimal
    /// string: the decimal "unit in the last place". E.g. "23.4393" -> 1e-4,
    /// "84381.406" -> 1e-3, "3.9860e14" -> 1e10. Used to derive the honest
    /// truncation bound of a finite-digit source value.
    static T decimal_ulp(const char* s) {
        int frac = 0;            // digits after the decimal point
        int exp = 0;             // explicit exponent magnitude
        int exp_sign = 1;
        bool after_dot = false;
        bool in_exp = false;
        for (const char* p = s; *p; ++p) {
            char c = *p;
            if (c == '.') { after_dot = true; }
            else if (c == 'e' || c == 'E') { in_exp = true; }
            else if (c == '-') { if (in_exp) exp_sign = -1; }
            else if (c >= '0' && c <= '9') {
                if (in_exp) { exp = exp * 10 + (c - '0'); }
                else if (after_dot) { ++frac; }
            }
        }
        int place = exp_sign * exp - frac;  // power of ten of the last digit
        T ulp = T(1);
        if (place >= 0) { for (int i = 0; i < place; ++i) ulp = ulp * T(10); }
        else { for (int i = 0; i < -place; ++i) ulp = ulp / T(10); }
        return ulp;
    }

    /// A constant known only to the digits written in `s` (a book / Almanac value
    /// with finite displayed precision). Records the honest representation bound
    /// 0.5 * decimal_ulp(s) — the "Taylor bound for the current representation" —
    /// in the precision category, as the max with the binary representation cost.
    /// Unlike defined(), this bound does NOT shrink with a wider T: the source
    /// supplies no further digits. Prefer a generative exact/series form (which
    /// does scale with T) whenever one exists.
    static TrackedValue from_truncated_decimal(const char* s) {
        using std::max;
        T val = from_string(s);
        T prec = max(decimal_ulp(s) / T(2), representation_bound(val));
        return TrackedValue(val, T(0), prec, T(0));
    }

    /// A model-DEFINING coefficient written to a finite number of decimal digits
    /// (e.g. an SR3 lunisolar DSCOM/DPPER/resonance constant, or an adopted
    /// physical eccentricity / mean-motion). It separates the two distinct error
    /// sources that from_truncated_decimal conflates into one precision number:
    ///   precision = 0.5 binary ULP  — how well the *given* decimal is stored in
    ///               type T; purely COMPUTATIONAL, so it tightens with a wider T;
    ///   accuracy  = 0.5 * decimal_ulp(s)  — the finite-digit MODEL-FIDELITY
    ///               limit (a more-digits model would differ); T-INDEPENDENT.
    /// The TOTAL error matches from_truncated_decimal (so it is no less honest),
    /// but the model floor lives in accuracy instead of masking the computational
    /// precision — which is what lets a wider T sharpen the tracked precision
    /// through a model built on finite-digit coefficients (the precision/accuracy
    /// separation of feedback_precision_accuracy). Unlike defined(), the
    /// finite-digit error is NOT dropped — it is recorded in accuracy, so the
    /// total is never under-claimed (the CR1 honesty requirement).
    static TrackedValue model_coefficient(const char* s) {
        T val = from_string(s);
        T prec = representation_bound(val);
        T acc = decimal_ulp(s) / T(2);
        return TrackedValue(val, T(0), prec, acc);
    }

    // --- Error queries ---

    /// Shorthand for errors.total() — the conservative summed bound.
    T total_error() const { return errors.total(); }

    /// Decimal digits of the value the total error leaves intact
    /// (floor(-log10(total/|value|)); 0 when the error dominates).
    int reliable_digits() const {
        using std::abs; using std::log10; using std::floor;
        T te = total_error();
        if (te == T(0)) return std::numeric_limits<int>::max(); // exact representation
        T av = abs(value);
        if (av == T(0) || te >= av) return 0; // error dominates
        return static_cast<int>(floor(-log10(te / av)));
    }

    // --- Arithmetic: addition ---

    /// Bound: |f(x+dx, y+dy) - f(x,y)| <= dx + dy (exact for addition).
    friend TrackedValue operator+(const TrackedValue& a, const TrackedValue& b) {
        return TrackedValue(a.value + b.value, a.errors + b.errors);
    }

    // --- Arithmetic: subtraction ---

    /// Same absolute error bound as addition. Subtractive cancellation does not
    /// increase absolute error, but relative error may grow because |a-b| may
    /// be small. reliable_digits() detects this.
    friend TrackedValue operator-(const TrackedValue& a, const TrackedValue& b) {
        return TrackedValue(a.value - b.value, a.errors + b.errors);
    }

    // --- Arithmetic: unary negation ---
    // No change in any error.

    TrackedValue operator-() const {
        return TrackedValue(-value, errors);
    }

    // --- Arithmetic: multiplication ---

    /// Bound: |(x+δx)(y+δy) - xy| ≤ |x|·δy + |y|·δx + δx·δy.
    ///
    /// Applied per error category. Cross-category products (e.g., a's measurement
    /// error × b's precision error) are neglected — they are second-order small
    /// when individual errors are small relative to values. This is standard
    /// interval arithmetic practice for separated error budgets.
    friend TrackedValue operator*(const TrackedValue& a, const TrackedValue& b) {
        using std::abs;
        T av = abs(a.value);
        T bv = abs(b.value);

        // Per-category bound |x·y| ≤ |a|·δb + |b|·δa + δa·δb, each product guarded
        // against 0 × overflowed-error → NaN (see math::detail::tv_safe_mul).
        return TrackedValue(
            a.value * b.value,
            ThreeErrors<T>(
                detail::tv_mul_bound(av, bv, a.errors.measurement, b.errors.measurement),
                detail::tv_mul_bound(av, bv, a.errors.precision, b.errors.precision),
                detail::tv_mul_bound(av, bv, a.errors.accuracy, b.errors.accuracy)
            )
        );
    }

    // --- Arithmetic: division ---

    /// Bound: |(x/y) - ((x+dx)/(y+dy))| ≤ (|x|·δy + |y|·δx) / (|y|·(|y| - δy)).
    /// Uses minimum possible denominator magnitude for rigor.
    friend TrackedValue operator/(const TrackedValue& a, const TrackedValue& b) {
        using std::abs;
        T av = abs(a.value);
        T bv = abs(b.value);

        // Per-category bound (|num|·δden + |den|·δnum) / (|den|·(|den|−δden)),
        // using the minimum denominator magnitude for rigor; the 0 × overflowed-
        // error guard (math::detail::tv_safe_mul) keeps a zero numerator at 0.
        return TrackedValue(
            a.value / b.value,
            ThreeErrors<T>(
                detail::tv_div_bound(av, bv, a.errors.measurement, b.errors.measurement),
                detail::tv_div_bound(av, bv, a.errors.precision, b.errors.precision),
                detail::tv_div_bound(av, bv, a.errors.accuracy, b.errors.accuracy)
            )
        );
    }

    // --- Math functions ---

    /// sqrt: |√(x+d) - √x| ≤ d / (2√(x-d))  when d < x
    ///
    /// When d ≥ x, the argument could be zero or negative, meaning √(x-d) is
    /// undefined. In that case the worst-case change is the entire result √x
    /// (the value could be anywhere from 0 to √(x+d)).
    friend TrackedValue sqrt(const TrackedValue& a) {
        using std::sqrt;
        T val = sqrt(a.value);

        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err >= a.value) {
                // Error exceeds value — result could be zero. Entire value is uncertain.
                return val;
            }
            return err / (T(2) * sqrt(a.value - err));
        }));
    }

    /// sin: The input x has error δ, so we compute sin(x + δ) instead of sin(x).
    ///
    /// Bound: |sin(x+δ) − sin(x)| ≤ |cos(x)|·δ + δ²/2
    ///
    /// The linear term |cos(x)|·δ is the operative bound — it is the derivative
    /// of sin times the input error. For representation-scale δ (~1e-16), this
    /// completely determines the output error.
    ///
    /// The δ²/2 term is from the Taylor remainder (|sin(ξ)| ≤ 1 for the unknown
    /// intermediate ξ). For δ < 1, δ² < δ so this term is always smaller than
    /// the linear term — it tightens the bound, never loosens it.
    ///
    /// Absolute cap: |sin(x+δ) − sin(x)| ≤ 2 always.
    friend TrackedValue sin(const TrackedValue& a) {
        using std::sin; using std::cos; using std::abs; using std::min;
        T s = sin(a.value);
        T abs_cos = abs(cos(a.value));

        return TrackedValue(s, a.errors.apply([&](T err) -> T {
            T linear = abs_cos * err;
            T quadratic = err * err / T(2);  // |sin(ξ)| ≤ 1
            return min(linear + quadratic, T(2));
        }));
    }

    /// cos: The input x has error δ, so we compute cos(x + δ) instead of cos(x).
    ///
    /// Bound: |cos(x+δ) − cos(x)| ≤ |sin(x)|·δ + δ²/2
    ///
    /// The linear term |sin(x)|·δ is the operative bound. For representation-
    /// scale δ (~1e-16), it completely determines the output error.
    ///
    /// The δ²/2 term is from the Taylor remainder (|cos(ξ)| ≤ 1). For δ < 1,
    /// δ² < δ so this is always smaller than the linear term.
    ///
    /// Absolute cap: |cos(x+δ) − cos(x)| ≤ 2 always.
    friend TrackedValue cos(const TrackedValue& a) {
        using std::sin; using std::cos; using std::abs; using std::min;
        T cv = cos(a.value);
        T abs_sin = abs(sin(a.value));

        return TrackedValue(cv, a.errors.apply([&](T err) -> T {
            T linear = abs_sin * err;
            T quadratic = err * err / T(2);  // |cos(ξ)| ≤ 1
            return min(linear + quadratic, T(2));
        }));
    }

    /// atan: |atan(x+d) − atan(x)| ≤ |d|/(1+x²), capped at π.
    /// The derivative 1/(1+x²) ≤ 1, so atan never amplifies errors.
    friend TrackedValue atan(const TrackedValue& a) {
        using std::atan; using std::min;
        T val = atan(a.value);
        T deriv = T(1) / (T(1) + a.value * a.value);
        T pi_val = boost::math::constants::pi<T>();

        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            return min(deriv * err, pi_val);
        }));
    }

    /// atan2: partial derivatives are ∂/∂y = x/(x²+y²), ∂/∂x = −y/(x²+y²).
    /// Near the origin (x≈0, y≈0), r²→0 and the derivatives blow up.
    /// The angle is completely undetermined there, so bound = π.
    friend TrackedValue atan2(const TrackedValue& y, const TrackedValue& x) {
        using std::atan2; using std::abs; using std::min;
        T val = atan2(y.value, x.value);
        T pi_val = boost::math::constants::pi<T>();

        T r2 = x.value * x.value + y.value * y.value;

        // Guard: if r² is tiny relative to the errors, angle is undetermined
        T total_err = y.errors.total() + x.errors.total();
        if (r2 <= total_err * total_err) {
            return TrackedValue(val, ThreeErrors<T>(pi_val, pi_val, pi_val));
        }

        T dy_factor = abs(x.value) / r2;
        T dx_factor = abs(y.value) / r2;

        return TrackedValue(
            val,
            ThreeErrors<T>(
                min(dy_factor * y.errors.measurement + dx_factor * x.errors.measurement, pi_val),
                min(dy_factor * y.errors.precision + dx_factor * x.errors.precision, pi_val),
                min(dy_factor * y.errors.accuracy + dx_factor * x.errors.accuracy, pi_val)
            )
        );
    }

    /// abs: error bounds unchanged (abs is Lipschitz-1)
    friend TrackedValue abs(const TrackedValue& a) {
        using std::abs;
        return TrackedValue(abs(a.value), a.errors);
    }

    /// fmod: fmod(x,y) = x − n·y where n = trunc(x/y).
    /// Error bound: δ(fmod) ≤ δx + |n|·δy for each category.
    /// (n is an integer, so it contributes no error itself, but it scales y's error.)
    friend TrackedValue fmod(const TrackedValue& a, const TrackedValue& b) {
        using std::fmod; using std::abs; using std::floor;
        T result = fmod(a.value, b.value);
        T n_abs = abs(floor(a.value / b.value));

        return TrackedValue(result, a.errors + n_abs * b.errors);
    }

    /// exp: exp is convex and increasing, so for an input error d the maximum
    /// deviation over [x-d, x+d] occurs at +d and equals exp(x)*(exp(d) - 1).
    /// This is the exact forward bound (no Taylor truncation) — what the
    /// mean-value theorem majorises — and it reduces to exp(x)*d for small d.
    /// expm1(d) computes exp(d)-1 without the catastrophic cancellation that
    /// rounds the naive form to zero for representation-scale d (which would
    /// silently drop the precision budget — the exact failure B1 guards against).
    /// Propagated per error category, so a wider T tightens the precision bound.
    friend TrackedValue exp(const TrackedValue& a) {
        using std::exp; using std::expm1;
        T val = exp(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            return val * expm1(err);
        }));
    }

    /// log (natural): log is concave and increasing, so for an input error d the
    /// maximum deviation over [x-d, x+d] occurs at -d, bounded via the
    /// mean-value theorem by log(x) - log(x-d) <= d / (x - d) (the derivative
    /// 1/(x-d) at the worst point). When d >= x the argument may be non-positive
    /// and log is unbounded below; we then return a large finite bound to flag
    /// the result as unreliable (cf. sqrt's degenerate guard above).
    friend TrackedValue log(const TrackedValue& a) {
        using std::log; using std::abs;
        T val = log(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err >= a.value) {
                return abs(val) + T(1);  // argument may be <= 0; result unbounded
            }
            return err / (a.value - err);
        }));
    }

    /// cbrt: cube root, total over R and monotonic increasing. Its derivative
    /// 1/(3*cbrt(t)^2) grows as |t| -> 0, so over [x-d, x+d] the worst base is the
    /// endpoint nearest zero. Written as d*(derivative), NOT as a difference of
    /// two cbrt values, to avoid the cancellation that would zero the precision
    /// budget for representation-scale d (cf. exp/expm1). When the interval
    /// straddles 0 the derivative is unbounded; we then return the full excursion
    /// as a finite large bound.
    friend TrackedValue cbrt(const TrackedValue& a) {
        using std::cbrt; using std::abs;
        T val = cbrt(a.value);
        T ax = abs(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err >= ax) {
                return abs(cbrt(ax + err));  // straddles 0; bound the full swing
            }
            T c = cbrt(ax - err);            // base nearest zero (by magnitude)
            return err / (T(3) * c * c);
        }));
    }

    /// pow(x, p) for a constant exponent p (base x > 0): derivative is p*x^(p-1).
    /// Over [x-d, x+d], |x^(p-1)| is maximised toward 0 when p<1 (since p-1<0) and
    /// away from 0 when p>1. Written as d*(derivative) to avoid cancellation. If
    /// the worst base is non-positive the result is ill-defined and we flag a
    /// large finite bound. (A tracked-exponent pow(x,y)=exp(y*log x) is deferred
    /// until a use needs it.)
    friend TrackedValue pow(const TrackedValue& a, const T& p) {
        using std::pow; using std::abs;
        T val = pow(a.value, p);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err <= T(0)) return T(0);
            T worst = (p < T(1)) ? (a.value - err) : (a.value + err);
            if (worst <= T(0)) {
                return abs(val) + T(1);  // base may be <= 0; pow ill-defined
            }
            return abs(p) * pow(worst, p - T(1)) * err;
        }));
    }

    /// Tracked-base, tracked-exponent pow: x^y = exp(y * log x), for base x > 0.
    /// The deferred case above, now that tracked exp/log exist (B-series). Rather
    /// than re-derive a two-variable forward bound, COMPOSE the already-verified
    /// tracked log, operator*, and exp -- the three-error budget then propagates
    /// rigorously through BOTH the base and the exponent, and log()'s degenerate
    /// guard (d >= x) handles a non-positive / straddling base. Not on the SGP4
    /// path (SGP4 uses std::pow / the constant-exponent overload above).
    friend TrackedValue pow(const TrackedValue& x, const TrackedValue& y) {
        return exp(y * log(x));
    }

    /// asin: domain [-1, 1], derivative 1/sqrt(1-x^2) which blows up toward
    /// |x| -> 1. Over [x-d, x+d] the worst (largest) derivative is at the
    /// endpoint nearest +-1, i.e. base magnitude |x|+d. Written as d*(derivative)
    /// to avoid the cancellation a difference-of-asin would suffer. If |x|+d
    /// reaches the domain edge the derivative is unbounded; the result range is
    /// [-pi/2, pi/2] (width pi), so we cap/flag the excursion at pi.
    friend TrackedValue asin(const TrackedValue& a) {
        using std::asin; using std::sqrt; using std::abs; using std::min;
        T val = asin(a.value);
        T ax = abs(a.value);
        T pi_val = boost::math::constants::pi<T>();
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err <= T(0)) return T(0);
            T worst = ax + err;             // nearest the +-1 singularity
            if (worst >= T(1)) return pi_val;
            return min(err / sqrt(T(1) - worst * worst), pi_val);
        }));
    }

    /// acos: domain [-1, 1], derivative -1/sqrt(1-x^2) (same magnitude as asin).
    /// Range [0, pi] (width pi), so the degenerate/cap bound is pi.
    friend TrackedValue acos(const TrackedValue& a) {
        using std::acos; using std::sqrt; using std::abs; using std::min;
        T val = acos(a.value);
        T ax = abs(a.value);
        T pi_val = boost::math::constants::pi<T>();
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err <= T(0)) return T(0);
            T worst = ax + err;
            if (worst >= T(1)) return pi_val;
            return min(err / sqrt(T(1) - worst * worst), pi_val);
        }));
    }

    /// tan: derivative sec^2(x) = 1/cos^2(x), unbounded at the poles x = pi/2+k*pi.
    /// Over [x-d, x+d] (with d < pi/2 so at most one pole can lie inside, since
    /// poles are pi apart) the largest sec^2 occurs where |cos| is smallest. With
    /// no pole inside, |cos| is minimised at an endpoint, so the worst sec^2 is
    /// 1/min(cos(x-d)^2, cos(x+d)^2). A sign change of cos across the interval
    /// (cm*cp <= 0), or d >= pi/2, means a pole is in range and tan is unbounded:
    /// the rigorous bound is then +infinity (numeric max sentinel).
    friend TrackedValue tan(const TrackedValue& a) {
        using std::tan; using std::cos; using std::abs; using std::min;
        T val = tan(a.value);
        T half_pi = boost::math::constants::pi<T>() / T(2);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err <= T(0)) return T(0);
            if (err >= half_pi) return std::numeric_limits<T>::max();
            T cm = cos(a.value - err);
            T cp = cos(a.value + err);
            if (cm * cp <= T(0)) return std::numeric_limits<T>::max();  // pole inside
            T min_abs = min(abs(cm), abs(cp));
            return err / (min_abs * min_abs);
        }));
    }

    /// sinh: derivative cosh(x) > 0, increasing in |x|, so the worst derivative on
    /// [x-d, x+d] is cosh(|x|+d). Written d*(derivative); cosh is computed at the
    /// worst base directly (no cancellation). Defined on all of R.
    friend TrackedValue sinh(const TrackedValue& a) {
        using std::sinh; using std::cosh; using std::abs;
        T val = sinh(a.value);
        T ax = abs(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            return cosh(ax + err) * err;
        }));
    }

    /// cosh: derivative sinh(x); |sinh| is increasing in |x|, worst at |x|+d.
    /// sinh(|x|+d) >= 0 since |x|+d >= 0. Defined on all of R.
    friend TrackedValue cosh(const TrackedValue& a) {
        using std::sinh; using std::abs;
        T val = cosh(a.value);
        T ax = abs(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            return sinh(ax + err) * err;
        }));
    }

    /// tanh: derivative sech^2(x) = 1 - tanh^2(x) <= 1, maximal at x = 0. Over
    /// [x-d, x+d] the worst (largest) derivative is at the base nearest 0: 0 when
    /// the interval straddles 0 (d >= |x|), else |x|-d. Range (-1, 1) width 2.
    friend TrackedValue tanh(const TrackedValue& a) {
        using std::tanh; using std::cosh; using std::abs; using std::min;
        T val = tanh(a.value);
        T ax = abs(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err <= T(0)) return T(0);
            T base = (err >= ax) ? T(0) : (ax - err);  // nearest 0 => largest sech^2
            T ch = cosh(base);
            T sech2 = T(1) / (ch * ch);
            return min(sech2 * err, T(2));
        }));
    }

    /// expm1: exp(x) - 1. Same derivative as exp, and the forward difference
    /// expm1(x+d) - expm1(x) = exp(x+d) - exp(x) = exp(x)*(exp(d) - 1) is the
    /// EXACT convex forward bound, identical to exp's. expm1(err) avoids the
    /// cancellation that would zero the bound for representation-scale err.
    friend TrackedValue expm1(const TrackedValue& a) {
        using std::exp; using std::expm1;
        T val = expm1(a.value);
        T ex = exp(a.value);
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            return ex * expm1(err);
        }));
    }

    /// log1p: log(1+x), domain x > -1. Concave (like log), so the worst case over
    /// [x-d, x+d] is at -d: log1p(x) - log1p(x-d) <= d/(1+x-d). When d >= 1+x the
    /// argument may be <= 0 and log1p is unbounded below; flag a finite bound.
    friend TrackedValue log1p(const TrackedValue& a) {
        using std::log1p; using std::abs;
        T val = log1p(a.value);
        T onep = T(1) + a.value;
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err >= onep) return abs(val) + T(1);
            return err / (onep - err);
        }));
    }

    /// log10: log(x)/ln(10). Same concave bound as log scaled by 1/ln(10).
    friend TrackedValue log10(const TrackedValue& a) {
        using std::log10; using std::log; using std::abs;
        T val = log10(a.value);
        T ln10 = log(T(10));
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err >= a.value) return abs(val) + T(1);
            return err / ((a.value - err) * ln10);
        }));
    }

    /// log2: log(x)/ln(2). Same concave bound as log scaled by 1/ln(2).
    friend TrackedValue log2(const TrackedValue& a) {
        using std::log2; using std::log; using std::abs;
        T val = log2(a.value);
        T ln2 = log(T(2));
        return TrackedValue(val, a.errors.apply([&](T err) -> T {
            if (err >= a.value) return abs(val) + T(1);
            return err / ((a.value - err) * ln2);
        }));
    }

    /// hypot: sqrt(x^2 + y^2). The gradient (x/r, y/r) has each component <= 1 in
    /// magnitude, so per category the bound is (|x|/r)*dx + (|y|/r)*dy, which is
    /// also <= dx + dy (Lipschitz-1 in each argument). At r = 0 the function is
    /// not differentiable but is still Lipschitz-1, so dx + dy is the rigorous
    /// bound there. Binary, so we combine the two inputs' categories directly.
    friend TrackedValue hypot(const TrackedValue& x, const TrackedValue& y) {
        using std::hypot; using std::abs;
        T r = hypot(x.value, y.value);
        if (r == T(0)) {
            return TrackedValue(r, x.errors + y.errors);  // |hypot(dx,dy)| <= dx+dy
        }
        T fx = abs(x.value) / r;  // <= 1
        T fy = abs(y.value) / r;  // <= 1
        return TrackedValue(r, ThreeErrors<T>(
            fx * x.errors.measurement + fy * y.errors.measurement,
            fx * x.errors.precision   + fy * y.errors.precision,
            fx * x.errors.accuracy    + fy * y.errors.accuracy));
    }

    /// min: selects the smaller value. min is jointly 1-Lipschitz in the sup norm,
    /// so |min(a*,b*) - min(a,b)| <= max(da, db). When the two values are
    /// definitively separated (gap exceeds the summed total errors) the result is
    /// unambiguously one input and carries exactly that input's errors; otherwise
    /// (the error bars overlap) the per-category max is the rigorous bound.
    friend TrackedValue min(const TrackedValue& a, const TrackedValue& b) {
        using std::abs;
        const TrackedValue& chosen = (a.value <= b.value) ? a : b;
        const TrackedValue& other  = (a.value <= b.value) ? b : a;
        T gap = abs(a.value - b.value);
        if (gap > chosen.errors.total() + other.errors.total()) {
            return chosen;  // error bars do not overlap: choice is certain
        }
        return TrackedValue(chosen.value, ThreeErrors<T>::max_per_category(a.errors, b.errors));
    }

    /// max: selects the larger value; mirror of min.
    friend TrackedValue max(const TrackedValue& a, const TrackedValue& b) {
        using std::abs;
        const TrackedValue& chosen = (a.value >= b.value) ? a : b;
        const TrackedValue& other  = (a.value >= b.value) ? b : a;
        T gap = abs(a.value - b.value);
        if (gap > chosen.errors.total() + other.errors.total()) {
            return chosen;
        }
        return TrackedValue(chosen.value, ThreeErrors<T>::max_per_category(a.errors, b.errors));
    }

    /// clamp(x, lo, hi) = max(lo, min(x, hi)); inherits min/max's rigorous error
    /// handling. Used to force an argument into a valid domain (e.g. [-1, 1]
    /// before asin/acos) while keeping a faithful error budget.
    friend TrackedValue clamp(const TrackedValue& x,
                              const TrackedValue& lo, const TrackedValue& hi) {
        return max(lo, min(x, hi));
    }

    // --- Comparisons (on value only; errors are metadata) ---

    /// Value-only comparison; the error budgets do not participate.
    friend bool operator<(const TrackedValue& a, const TrackedValue& b) {
        return a.value < b.value;
    }
    /// Value-only comparison; the error budgets do not participate.
    friend bool operator>(const TrackedValue& a, const TrackedValue& b) {
        return a.value > b.value;
    }
    /// Value-only comparison; the error budgets do not participate.
    friend bool operator<=(const TrackedValue& a, const TrackedValue& b) {
        return a.value <= b.value;
    }
    /// Value-only comparison; the error budgets do not participate.
    friend bool operator>=(const TrackedValue& a, const TrackedValue& b) {
        return a.value >= b.value;
    }
    /// Value-only equality; the error budgets do not participate.
    friend bool operator==(const TrackedValue& a, const TrackedValue& b) {
        return a.value == b.value;
    }
    /// Value-only inequality; the error budgets do not participate.
    friend bool operator!=(const TrackedValue& a, const TrackedValue& b) {
        return a.value != b.value;
    }

    /// Implicit conversion from an exact integer (zero error, all categories).
    TrackedValue(int n) : value(T(n)), errors(T(0), T(0), T(0)) {}

    /// Rigorous upper bound on the representation error of val in type T.
    ///
    /// IEEE 754: a real number x is stored as fl(x), the nearest representable
    /// float. Within the binade [2^e, 2^(e+1)), representable numbers are
    /// spaced at 1 ULP = 2^e × ε, where ε = epsilon() = 2^(1−p).
    ///
    /// Round-to-nearest-even (IEEE default) guarantees |x − fl(x)| ≤ 0.5 ULP.
    ///
    /// We compute 0.5 ULP without using floating-point division:
    /// - For double: extract the binade exponent from the IEEE 754 bit pattern
    ///   using integer operations, then use ldexp (exact power-of-2 multiply).
    /// - For multiprecision: frexp returns an integer exponent; ldexp applies it.
    ///   Both are exact operations on the exponent field, no rounding involved.
    ///
    /// epsilon() = 2^(1-p). 0.5 ULP at binade exponent e = 2^(e-1) × 2^(1-p) = 2^(e-p).
    /// We compute this as ldexp(1, e - p) where p = significand bits.
    static T representation_bound(const T& val) {
        using std::abs;
        T av = abs(val);
        if (av == T(0)) return T(0);

        // Extract binade exponent. frexp returns integer exp such that
        // val = frac × 2^exp with 0.5 ≤ frac < 1, so val ∈ [2^(exp-1), 2^exp).
        // The binade base is 2^(exp-1), and 1 ULP = 2^(exp-1) × ε.
        //
        // For double: frexp operates on the IEEE 754 bit pattern — the exponent
        // extraction is an integer shift, not a floating-point operation.
        // For multiprecision: boost provides frexp that returns the base-2 exponent.

        if constexpr (std::is_same_v<T, double>) {
            int exp;
            std::frexp(av, &exp);
            // 0.5 ULP = 2^(exp - 1 - 52) = 2^(exp - 53)
            // ldexp(1.0, n) is exact: it sets the exponent field, no rounding.
            return std::ldexp(1.0, exp - 53);  // magic-ok: 1.0 = exact ldexp mantissa (2^n)
        } else if constexpr (std::is_same_v<T, float>) {
            int exp;
            std::frexp(static_cast<float>(av), &exp);
            // float: 24-bit significand. 0.5 ULP = 2^(exp - 1 - 23) = 2^(exp - 24)
            return std::ldexp(1.0f, exp - 24);  // magic-ok: 1.0f = exact ldexp mantissa (2^n)
        } else {
            // Multiprecision types: epsilon() = 2^(1-p) where p = digits (in bits).
            // 0.5 ULP at binade exp = 2^(exp-1) × epsilon / 2 = 2^(exp-1) × 2^(-p)
            //                       = 2^(exp - 1 - p)
            //
            // For types where frexp/ldexp are available (boost multiprecision):
            // ldexp(T(1), exp - 1 - digits) gives exact 0.5 ULP.
            //
            // Fallback for types without frexp: |val| × ε / 2 is a valid upper
            // bound (overestimates by at most 2× within a binade).
            int digits = std::numeric_limits<T>::digits; // significand bits
            if (digits > 0) {
                // Use the portable formula: 0.5 ULP ≤ |val| × 2^(-digits)
                // This is |val| × epsilon / 2, which is between 0.5 and 1.0 ULP.
                // It's the tightest portable bound without frexp.
                //
                // `using std::ldexp` plus an unqualified call enables ADL so
                // boost::multiprecision's own `ldexp` overload is picked for
                // its types. Under stricter MSVC instantiation (VS 2026 /
                // v145), the `std::`-qualified form does not consider ADL
                // candidates and fails to compile for non-builtin T.
                using std::ldexp;
                T half_eps = ldexp(T(1), -digits);
                return av * half_eps;
            }
            return av * std::numeric_limits<T>::epsilon();
        }
    }

private:
    /// Convert a string to type T.
    /// For double/float: uses std::stod/stof.
    /// For multiprecision: uses T's string constructor.
    static T from_string(const char* str) {
        if constexpr (std::is_same_v<T, double>) {
            return std::stod(str);
        } else if constexpr (std::is_same_v<T, float>) {
            return std::stof(str);
        } else {
            return T(str);
        }
    }

};

// --- Convenience constructors ---

/// Free-function shorthand for TrackedValue<T>::exact_integer(n).
template<typename T>
TrackedValue<T> exact(int n) {
    return TrackedValue<T>::exact_integer(n);
}

/// Exact rational number num/den.
///
/// The mathematical value num/den is exact, but the binary floating-point
/// representation of irrational fractions (1/3, 1/7, etc.) introduces
/// precision error bounded by 0.5 ULP (IEEE round-to-nearest-even).
///
/// Note: the division operator on two exact TrackedValues would report
/// zero precision error (since both inputs have zero error). But the
/// quotient T(num)/T(den) itself may not be exactly representable.
/// We compute the representation bound of the result directly.
template<typename T>
TrackedValue<T> ratio(int num, int den) {
    T val = T(num) / T(den);
    T repr_err = TrackedValue<T>::representation_bound(val);
    return TrackedValue<T>(val, T(0), repr_err, T(0));
}

// --- Truncation/tail-bound deposit (the error-channel chokepoint) -------------

/// Which rigorous-error category a truncation / tail bound is deposited into.
/// LOAD-BEARING (constants-honesty CR1B): the PRECISION channel is a
/// convergent-MATH truncation that tightens with a wider numeric type T or more
/// terms (Leibniz / geometric / Newton residual); the ACCURACY channel is the
/// T-independent MODEL floor — a truncated-PHYSICS series (omitted obliquity
/// terms, omitted periodics, omitted gravity degrees). There is deliberately NO
/// default: every caller states the channel explicitly so the choice is auditable.
enum class ErrorChannel { precision, accuracy };

/// Deposit a non-negative rigorous `bound` into ONE explicit error category. The
/// single sanctioned chokepoint for widening a TrackedValue's error budget after
/// the fact — it replaces the hand-typed `v.errors.<field> = v.errors.<field> +
/// bound` splices so the channel choice is greppable and unit-auditable. Pure
/// (by value); the deposited sum and the untouched channels/value are
/// bit-identical to the manual form (same operands, same `+`, same field).
template<typename T>
TrackedValue<T> add_bound(TrackedValue<T> v, const T& bound, ErrorChannel ch) {
    if (ch == ErrorChannel::precision) v.errors.precision = v.errors.precision + bound;
    else                               v.errors.accuracy  = v.errors.accuracy  + bound;
    return v;
}

} // namespace math
