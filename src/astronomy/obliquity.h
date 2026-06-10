#pragma once

/// @file obliquity.h
/// Series-based mean obliquity of the ecliptic (IAU 2006 / IERS Conventions
/// 2010 Eq. 5.40; Capitaine, Wallace & Chapront 2003, adopted IAU 2006).
///
/// The obliquity is NOT true by definition — it is the leading terms of a time
/// polynomial. Per the user directive (2026-06-05): a constant that is not
/// exact-by-definition must carry SERIES-based precision AND accuracy, tracked
/// through the series evaluation rather than stamped as a decimal. Concretely:
///
///   ε_A(t) = Σ_{k=0}^{5} c_k t^k   [arcsec],   t = TT Julian centuries from J2000
///
///   - PRECISION  = the representation/rounding of the Horner evaluation, which
///                  TIGHTENS with a wider numeric type T (the calling card);
///   - ACCURACY   = the series-TRUNCATION bound Σ_{k≥n} |c_k| |t|^k (the omitted
///                  higher-order terms), which TIGHTENS as more terms n are kept
///                  — the honest physical error of using a finite expansion,
///                  NOT the typographic "how many digits were written" floor.
///
/// The coefficients are born-digital (IERS Conventions 2010 Eq. 5.40), each a
/// finite-digit model-fit coefficient (model_coefficient: its digits are a
/// model-fidelity accuracy, its binary storage a T-scaling precision). The
/// arcsec→radian conversion uses the generative π/648000 (π to full T).

#include "../math/angles.h"
#include "../math/tracked_polynomial.h"   // TrackedPolynomial: Horner + tracked tail
#include "../math/tracked_value.h"

#include <vector>

namespace astronomy {

/// IAU 2006 mean obliquity ε_A(t) [radians], t = TT Julian centuries from
/// J2000.0, summing the first `n_terms` polynomial terms (1..6). The returned
/// value carries the series-truncation accuracy (omitted terms n..5 plus a
/// conservative degree-6 tail) and the representation precision.
template<typename T>
math::TrackedValue<T> obliquity_iau2006(const math::TrackedValue<T>& t_jcen,
                                        int n_terms = 6) {
    using TV = math::TrackedValue<T>;

    // IAU 2006 coefficients in arcseconds (IERS Conventions 2010, Eq. 5.40).
    static const char* kC[6] = {
        "84381.406", "-46.836769", "-0.0001831",
        "0.00200340", "-0.000000576", "-0.0000000434"};
    // Magnitudes for the conservative truncation bound only (mirror kC; an upper
    // bound need not be carried to full T precision).
    static const double kMag[6] = {
        84381.406, 46.836769, 0.0001831, 0.00200340, 0.000000576, 0.0000000434};

    // ε_A(t) as a TrackedPolynomial: Horner over the first n_terms, with the
    // Σ_{k≥n_terms}|c_k||t|^k truncation deposited in the ACCURACY channel plus
    // one conservative degree-6 tail term (extra_tail_terms = 1) — bit-identical
    // to the prior hand-rolled Horner + tail loops. Each c_k a model_coefficient
    // (digits -> accuracy, storage -> T-scaling precision).
    static const math::TrackedPolynomial<T> kEps(
        std::vector<TV>{TV::model_coefficient(kC[0]), TV::model_coefficient(kC[1]),
                        TV::model_coefficient(kC[2]), TV::model_coefficient(kC[3]),
                        TV::model_coefficient(kC[4]), TV::model_coefficient(kC[5])},
        std::vector<T>{T(kMag[0]), T(kMag[1]), T(kMag[2]),
                       T(kMag[3]), T(kMag[4]), T(kMag[5])});

    TV eps_arcsec = kEps.eval(t_jcen, n_terms, math::ErrorChannel::accuracy, 1);

    // arcsec -> radians via the generative π/648000 (π scales with T).
    return eps_arcsec * (math::pi<T>() / math::exact<T>(648000));
}

} // namespace astronomy
