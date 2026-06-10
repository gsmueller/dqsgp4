#pragma once

/**
 * @file moon_meeus.h
 * @brief Lunar ephemeris instance — Meeus §47 / ELP-2000 Poisson series (L3).
 *
 * The Moon's geocentric position in the MEAN ecliptic & equinox OF DATE (λ, β, Δ),
 * as a function of an `Epoch<T>` (L1). Unlike the two-body Sun, the Moon's motion
 * is dominated by SOLAR PERTURBATIONS (evection, variation, annual equation, …),
 * so it is a POISSON SERIES in the Delaunay arguments (D, M, M', F) — NOT a closed
 * equation of centre. The series is GENERATED (dial up `n_terms`, tracked
 * truncation) via the reusable `poisson_series` engine over the born-digital
 * Meeus §47 tables (lunar_terms_meeus.h). Theory: design/derivations/
 * ephemeris_series.md §3.
 *
 * @par Fidelity (no perceived fidelity — verified vs an INDEPENDENT oracle)
 * The amplitudes are born-digital (pymeeus's verbatim Meeus §47); the assembled
 * position is gated against **JPL DE430** (the in-repo numerical-integration
 * ephemeris, via erfa.eqec06) — NOT against astropy's `builtin` Moon, which is
 * itself the Meeus/ELP truncation (erfa.moon98) and so would only check the code
 * reproduces the theory. Measured residual (full 60-term): ≤ 3.7″ in λ, ≤ 1.5″ in
 * β, ≤ 5 km in Δ over 1995–2024 — within Meeus §47's documented ~10″/~4″ grade,
 * which is carried as the model-fidelity accuracy floor. Generator/oracle:
 * tools/gen_lunar_terms.py (data) + tools/gen_lunar_oracle.py (DE430 reference).
 *
 * @par References (born-digital)
 * - Meeus, Astronomical Algorithms (1998), Ch.47 — Tables 47.A/47.B (ELP-2000).
 * - JPL DE430 (in-repo sunmooneph_430t12.txt); ERFA eqec06 — the oracle.
 */

#include "../astronomy/epoch.h"          // Epoch, centuries_since_j2000 (L1)
#include "../math/angles.h"
#include "../math/tracked_value.h"
#include "ecliptic_state.h"
#include "lunar_terms_meeus.h"           // born-digital kLunarLR / kLunarB tables
#include "poisson_series.h"              // the reusable Poisson engine

#include <boost/math/constants/constants.hpp>
#include <vector>

namespace ephemeris {

/// Moon's geocentric position in the mean ecliptic & equinox of date, Meeus §47.
/// `n_terms` dials the periodic series up (1..60); the omitted-term truncation is
/// tracked in errors.accuracy. `epoch` is interpreted on the TT scale.
template<typename T>
EclipticState<T> moon_meeus_of_date(const astronomy::Epoch<T>& epoch, int n_terms = 60) {
    using TV = math::TrackedValue<T>;
    using math::exact;
    using math::ratio;
    using math::degrees_to_radians;
    using math::wrap_two_pi;

    TV t = astronomy::centuries_since_j2000(epoch.jd());
    TV t2 = t * t, t3 = t2 * t, t4 = t3 * t;

    // Mean elements [deg] (Meeus §47): decimal coeffs are finite-digit fits
    // (model_coefficient); the 1/N tail coeffs are exact rationals (ratio).
    TV Lp = TV::model_coefficient("218.3164477") + TV::model_coefficient("481267.88123421") * t
          + TV::model_coefficient("-0.0015786") * t2 + ratio<T>(1, 538841) * t3 - ratio<T>(1, 65194000) * t4;
    TV D = TV::model_coefficient("297.8501921") + TV::model_coefficient("445267.1114034") * t
          + TV::model_coefficient("-0.0018819") * t2 + ratio<T>(1, 545868) * t3 - ratio<T>(1, 113065000) * t4;
    TV Msun = TV::model_coefficient("357.5291092") + TV::model_coefficient("35999.0502909") * t
          + TV::model_coefficient("-0.0001536") * t2 + ratio<T>(1, 24490000) * t3;
    TV Mp = TV::model_coefficient("134.9633964") + TV::model_coefficient("477198.8675055") * t
          + TV::model_coefficient("0.0087414") * t2 + ratio<T>(1, 69699) * t3 - ratio<T>(1, 14712000) * t4;
    TV F = TV::model_coefficient("93.2720950") + TV::model_coefficient("483202.0175233") * t
          + TV::model_coefficient("-0.0036539") * t2 - ratio<T>(1, 3526000) * t3 + ratio<T>(1, 863310000) * t4;

    // Eccentricity correction E (applied as E^|cM| to terms in the Sun's anomaly).
    TV E = exact<T>(1) - TV::model_coefficient("0.002516") * t - TV::model_coefficient("0.0000074") * t2;
    TV E2 = E * E;

    TV Dr = degrees_to_radians(D), Mr = degrees_to_radians(Msun),
       Mpr = degrees_to_radians(Mp), Fr = degrees_to_radians(F), Lpr = degrees_to_radians(Lp);

    auto e_factor = [&](int cM) -> TV {
        int a = cM < 0 ? -cM : cM;
        return (a == 0) ? exact<T>(1) : (a == 1 ? E : E2);  // Meeus §47: |cM| ≤ 2
    };

    // Longitude (Σl, 1e-6°) and distance (Σr, 1e-3 km) share the (D,M,M',F) phases.
    const int NLR = static_cast<int>(sizeof(kLunarLR) / sizeof(kLunarLR[0]));
    std::vector<TV> arg_lr(NLR), amp_l(NLR), amp_r(NLR);
    for (int k = 0; k < NLR; ++k) {
        const LunarTermLR& r = kLunarLR[k];
        arg_lr[k] = exact<T>(r.cD) * Dr + exact<T>(r.cM) * Mr + exact<T>(r.cMp) * Mpr + exact<T>(r.cF) * Fr;
        TV e = e_factor(r.cM);
        amp_l[k] = exact<T>(static_cast<int>(r.sig_l)) * e;   // exact-integer amplitude
        amp_r[k] = exact<T>(static_cast<int>(r.sig_r)) * e;
    }
    TV sum_l = poisson_series(arg_lr, amp_l, n_terms, false);  // Σl [1e-6 deg]
    TV sum_r = poisson_series(arg_lr, amp_r, n_terms, true);   // Σr [1e-3 km]

    // Latitude (Σb, 1e-6°).
    const int NB = static_cast<int>(sizeof(kLunarB) / sizeof(kLunarB[0]));
    std::vector<TV> arg_b(NB), amp_b(NB);
    for (int k = 0; k < NB; ++k) {
        const LunarTermB& r = kLunarB[k];
        arg_b[k] = exact<T>(r.cD) * Dr + exact<T>(r.cM) * Mr + exact<T>(r.cMp) * Mpr + exact<T>(r.cF) * Fr;
        amp_b[k] = exact<T>(static_cast<int>(r.sig_b)) * e_factor(r.cM);
    }
    TV sum_b = poisson_series(arg_b, amp_b, n_terms, false);   // Σb [1e-6 deg]

    // Additive (Venus / Jupiter / Earth-flattening) terms, Meeus §47 — always on.
    TV A1 = degrees_to_radians(TV::model_coefficient("119.75") + TV::model_coefficient("131.849") * t);
    TV A2 = degrees_to_radians(TV::model_coefficient("53.09") + TV::model_coefficient("479264.290") * t);
    TV A3 = degrees_to_radians(TV::model_coefficient("313.45") + TV::model_coefficient("481266.484") * t);
    sum_l = sum_l + exact<T>(3958) * sin(A1) + exact<T>(1962) * sin(Lpr - Fr) + exact<T>(318) * sin(A2);
    sum_b = sum_b + exact<T>(-2235) * sin(Lpr) + exact<T>(382) * sin(A3)
                  + exact<T>(175) * sin(A1 - Fr) + exact<T>(175) * sin(A1 + Fr)
                  + exact<T>(127) * sin(Lpr - Mpr) + exact<T>(-115) * sin(Lpr + Mpr);

    // Assemble: Σl,Σb in 1e-6°; Σr in 1e-3 km added to the mean distance.
    TV micro_deg = degrees_to_radians(exact<T>(1)) / exact<T>(1000000);  // 1e-6° in radians
    EclipticState<T> s;
    s.lon = wrap_two_pi(Lpr + sum_l * micro_deg);
    s.lat = sum_b * micro_deg;
    s.radius = TV::model_coefficient("385000.56") + sum_r / exact<T>(1000);  // km

    // Model-fidelity floor (born-digital Meeus §47 grade ≈ 10″ λ / 4″ β; ~20 km Δ),
    // verified to MAJORIZE the measured JPL DE430 residual (≤ 3.7″/1.5″/5 km). For
    // n_terms < 60 the (larger) series-truncation accuracy dominates this floor.
    const T arcsec = boost::math::constants::pi<T>() / T(180) / T(3600);
    s.lon = math::add_bound(s.lon, T(10) * arcsec, math::ErrorChannel::accuracy);
    s.lat = math::add_bound(s.lat, T(4) * arcsec, math::ErrorChannel::accuracy);
    s.radius = math::add_bound(s.radius, T(20), math::ErrorChannel::accuracy);
    return s;
}

} // namespace ephemeris
