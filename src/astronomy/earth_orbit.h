#pragma once

/// @file earth_orbit.h
/// Series-based secular elements of Earth's heliocentric orbit (the Sun's
/// apparent orbit, as the lunisolar perturbation theory sees it).
///
/// Like the obliquity (astronomy/obliquity.h), Earth's orbital eccentricity is
/// NOT true by definition — it is the leading terms of a secular time series
/// driven by the planetary perturbations (Laplace–Lagrange secular theory,
/// realised numerically by VSOP). Per the user directive (2026-06-05): a
/// constant that is not exact-by-definition must carry SERIES-based precision
/// AND accuracy, tracked through the series evaluation rather than stamped as a
/// decimal. Concretely:
///
///   e(t) = Σ_{k=0}^{2} c_k t^k   [dimensionless],   t = TT Julian centuries
///                                                        from J2000.0
///
///   - PRECISION  = the representation/rounding of the Horner evaluation, which
///                  TIGHTENS with a wider numeric type T (the calling card);
///   - ACCURACY   = the series-TRUNCATION bound Σ_{k≥n} |c_k| |t|^k (the omitted
///                  higher-order secular terms), which TIGHTENS as more terms n
///                  are kept — the honest physical error of a finite expansion,
///                  NOT the typographic "how many digits were written" floor.
///
/// The coefficients are born-digital, each a finite-digit model-fit coefficient
/// (model_coefficient: its digits are a model-fidelity accuracy, its binary
/// storage a T-scaling precision). Source: Meeus, *Astronomical Algorithms*
/// (2nd ed.), Ch. 25, Eq. 25.4 (from VSOP87); realised verbatim in the pymeeus
/// reference implementation (Sun.true_longitude_coarse):
///   e = 0.016708634 - t*(0.000042037 + t*0.0000001267).
///
/// NB on the time unit: this is the PER-CENTURY form (t in Julian centuries),
/// matching the obliquity series' argument so a single satellite-epoch t feeds
/// both. The per-millennium VSOP form would scale c_1, c_2 by 10, 100 — a unit
/// trap this form deliberately avoids.

#include "../math/tracked_polynomial.h"   // TrackedPolynomial: Horner + tracked tail
#include "../math/tracked_value.h"

#include <vector>

namespace astronomy {

/// Earth's orbital (the apparent solar) eccentricity e(t) [dimensionless],
/// t = TT Julian centuries from J2000.0, summing the first `n_terms` polynomial
/// terms (1..3). The returned value carries the series-truncation accuracy
/// (omitted higher-order terms plus a conservative degree-3 tail) and the
/// representation precision.
template<typename T>
math::TrackedValue<T> earth_eccentricity(const math::TrackedValue<T>& t_jcen,
                                         int n_terms = 3) {
    using TV = math::TrackedValue<T>;

    // Meeus Eq. 25.4 (VSOP87) coefficients, dimensionless, t in Julian centuries.
    static const char* kC[3] = {
        "0.016708634", "-0.000042037", "-0.0000001267"};
    // Magnitudes for the conservative truncation bound only (mirror kC).
    static const double kMag[3] = {0.016708634, 0.000042037, 0.0000001267};

    // e(t) as a TrackedPolynomial: Horner over the first n_terms, with the
    // Σ_{k≥n_terms}|c_k||t|^k truncation in the ACCURACY channel plus one
    // conservative degree-3 tail term (extra_tail_terms = 1) — bit-identical to
    // the prior hand-rolled Horner + tail loops. Each c_k a model_coefficient
    // (digits -> accuracy, storage -> T-scaling precision).
    static const math::TrackedPolynomial<T> kEcc(
        std::vector<TV>{TV::model_coefficient(kC[0]), TV::model_coefficient(kC[1]),
                        TV::model_coefficient(kC[2])},
        std::vector<T>{T(kMag[0]), T(kMag[1]), T(kMag[2])});

    return kEcc.eval(t_jcen, n_terms, math::ErrorChannel::accuracy, 1);
}

} // namespace astronomy
