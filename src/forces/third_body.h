#pragma once

/// @file third_body.h
/// Cartesian third-body (Sun / Moon) gravitational perturbation — the Newtonian
/// tidal acceleration a satellite feels from a distant body, built from the L3
/// ephemeris (ephemeris/body_position_gcrs.h) and textbook physics. This is the
/// CARTESIAN force member, distinct from the orbit-averaged secular-rate form in
/// perturbation/third_body.h (which serves an averaged-element / SR3-style
/// consumer and stays unwired). Theory: design/derivations/third_body_perturbation.md.
///
///   a_3body = −(μ₃/‖s‖³) [ r·(1+q)^(−3/2) + f(q)·s ]      (Battin, cancellation-free)
///       q   = (r·(r − 2s))/‖s‖²
///       f(q)= q(3+3q+q²) / [ (1+q)^(3/2)·((1+q)^(3/2)+1) ]
///
/// r = satellite geocentric position, s = third-body geocentric position (same
/// frame, GCRS metres). The Battin form removes the catastrophic cancellation of
/// the naive μ₃[(s−r)/‖s−r‖³ − s/‖s‖³] (the two vectors agree to ~7 figures for a
/// LEO satellite — a perceived-fidelity trap). No perceived fidelity: gated in
/// tests/test_third_body against (1) the naive form at cpp_bin_float_50 [formula]
/// and (2) JPL DE430 body positions [ephemeris→acceleration].
///
/// This increment ships the CONSUMED SUBSET — the GCRS-frame acceleration. The
/// State→Wrench propagator force (with the state.time→epoch map and the TEME↔GCRS
/// frame step, §4/§7 of the note) is DEFERRED to its wiring increment, exactly as
/// geopotential was gated standalone before its DQ-propagator switch (anti-dead-code).

#include "force_common.h"

#include "../astronomy/epoch.h"
#include "../astronomy/frames.h"
#include "../constants/constants_provider.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../ephemeris/body_position_gcrs.h"
#include "../ephemeris/moon_meeus.h"
#include "../ephemeris/sun_meeus.h"
#include "../math/angles.h"
#include "../math/matrix3.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <cmath>
#include <functional>

namespace forces {

/// Newtonian third-body perturbing acceleration in the satellite's frame.
/// `r_sat` and `s` (third-body position) are geocentric, same frame, same units
/// (metres); `mu` is the body's GM [m³/s²]. Returns the perturbing acceleration
/// [m/s²]. Evaluated in the Battin `f(q)` form so the LEO cancellation never
/// occurs (third_body_perturbation.md §3); all error channels propagate from the
/// tracked inputs (the body-position accuracy carries the ephemeris truncation).
template<typename T>
math::Vector3<T> third_body_accel(const math::Vector3<T>& r_sat,
                                  const math::Vector3<T>& s,
                                  const math::TrackedValue<T>& mu) {
    using math::exact;
    using TV = math::TrackedValue<T>;

    TV s2 = s.x * s.x + s.y * s.y + s.z * s.z;     // ‖s‖²
    TV s_mag = sqrt(s2);
    TV s3 = s2 * s_mag;                            // ‖s‖³

    // q = (r·(r − 2s))/‖s‖² = (‖r‖² − 2 r·s)/‖s‖².
    TV r2 = r_sat.x * r_sat.x + r_sat.y * r_sat.y + r_sat.z * r_sat.z;
    TV r_dot_s = r_sat.x * s.x + r_sat.y * s.y + r_sat.z * s.z;
    TV q = (r2 - exact<T>(2) * r_dot_s) / s2;

    TV one = exact<T>(1);
    TV onepq = one + q;
    TV onepq32 = onepq * sqrt(onepq);             // (1+q)^(3/2)

    // f(q) = q(3 + 3q + q²) / [ (1+q)^(3/2)·((1+q)^(3/2) + 1) ]  (cancellation-free).
    TV fq_num = q * (exact<T>(3) + exact<T>(3) * q + q * q);
    TV fq = fq_num / (onepq32 * (onepq32 + one));

    TV inv_onepq32 = one / onepq32;               // (1+q)^(−3/2)
    TV coef = -(mu / s3);                         // −μ/‖s‖³

    return math::Vector3<T>(
        coef * (r_sat.x * inv_onepq32 + s.x * fq),
        coef * (r_sat.y * inv_onepq32 + s.y * fq),
        coef * (r_sat.z * inv_onepq32 + s.z * fq));
}

/// A perturbing third body: its GM and a position-of-epoch function returning the
/// body's GCRS position in METRES. The position closure wraps an L3 ephemeris
/// instance (sun_meeus / moon_meeus) through body_position_gcrs and the per-body
/// unit conversion — the body-genericity the ephemeris already realises.
template<typename T>
struct ThirdBody {
    math::TrackedValue<T> mu;                                          ///< GM [m³/s²]
    std::function<math::Vector3<T>(const astronomy::Epoch<T>&)> position;  ///< GCRS [m]
};

/// IAU 2012 astronomical unit [m] (exact by definition — `defined`, registered in
/// the CR1B convention allowlist). Converts the Sun's AU-valued ephemeris to metres.
template<typename T>
math::TrackedValue<T> astronomical_unit() {
    return math::TrackedValue<T>::defined("149597870700");
}

/// The Sun as a perturbing third body. μ☉ = IAU/DE430 heliocentric GM (adopted,
/// published-σ `measured`); position from the Meeus §25 solar ephemeris (AU → m).
template<typename T>
ThirdBody<T> sun_third_body() {
    using TV = math::TrackedValue<T>;
    TV mu = TV::measured("1.32712440018e20", "8e9");
    TV au = astronomical_unit<T>();
    auto position = [au](const astronomy::Epoch<T>& epoch) -> math::Vector3<T> {
        math::Vector3<T> g = ephemeris::body_position_gcrs<T>(
            ephemeris::sun_meeus_of_date<T>(epoch), epoch);
        return math::Vector3<T>(g.x * au, g.y * au, g.z * au);
    };
    return ThirdBody<T>{mu, position};
}

/// The Moon as a perturbing third body. μ☾ = DE430 lunar GM (adopted, published-σ
/// `measured`); position from the Meeus §47 lunar ephemeris (km → m).
template<typename T>
ThirdBody<T> moon_third_body() {
    using TV = math::TrackedValue<T>;
    TV mu = TV::measured("4.902800066e12", "1e4");
    TV km = math::exact<T>(1000);   // km → m, exact (not a defined-convention)
    auto position = [km](const astronomy::Epoch<T>& epoch) -> math::Vector3<T> {
        math::Vector3<T> g = ephemeris::body_position_gcrs<T>(
            ephemeris::moon_meeus_of_date<T>(epoch), epoch);
        return math::Vector3<T>(g.x * km, g.y * km, g.z * km);
    };
    return ThirdBody<T>{mu, position};
}

/// Perturbing acceleration on a satellite at `r_sat` (GCRS metres) from `body` at
/// absolute `epoch` — the §6 composition body_position_gcrs → third_body_accel.
template<typename T>
math::Vector3<T> third_body_perturbation(const math::Vector3<T>& r_sat,
                                         const ThirdBody<T>& body,
                                         const astronomy::Epoch<T>& epoch) {
    return third_body_accel<T>(r_sat, body.position(epoch), body.mu);
}

/// Body position in the propagator's ~TEME frame at `base` ⊕ `elapsed_s` — the
/// shared epoch-map + frame step BOTH the third-body and SRP forces use
/// (third_body_perturbation.md §9; solar_radiation_pressure.md §4). The absolute
/// epoch is formed two-part (preserving the day part); the body's GCRS position is
/// rotated by the DOMINANT erfa-gated IAU2006 precession; nutation + the equation
/// of equinoxes (≤ ~24″) are omitted — each consumer deposits its own
/// frame-fidelity accuracy bound on its acceleration.
template<typename T>
math::Vector3<T> third_body_position_teme(const ThirdBody<T>& body,
                                          const astronomy::Epoch<T>& base,
                                          const math::TrackedValue<T>& elapsed_s) {
    using TV = math::TrackedValue<T>;
    TV dt_days = elapsed_s / math::exact<T>(86400);
    astronomy::Epoch<T> epoch = astronomy::Epoch<T>::from_jd_two_part(
        base.jd_day, base.jd_frac + dt_days, base.scale);
    math::Vector3<T> s_gcrs = body.position(epoch);
    TV t = astronomy::centuries_since_j2000(epoch);
    return astronomy::precession_iau2006<T>(t) * s_gcrs;
}

/// Body-frame wrench for the propagator force list (third_body_perturbation.md §9).
/// The absolute epoch is `base` + elapsed seconds (state.time); the body's GCRS
/// position is rotated into the propagator's TEME frame by the DOMINANT IAU2006
/// precession (GCRS→mean-equatorial-of-date, erfa-gated via frames.h). Nutation +
/// the equation of equinoxes (≤ ~24″, below the Sun's Meeus grade) are omitted and
/// deposited as a frame-fidelity accuracy bound — honest, not a seam. `K` is unused
/// (the body carries its own μ); the signature matches the propagator's ForceFn.
template<typename T>
dynamics::Wrench<T> third_body_force(const dynamics::State<T>& state,
                                     const constants::ConstantsProvider<T>&,
                                     const ThirdBody<T>& body,
                                     const astronomy::Epoch<T>& base) {
    using TV = math::TrackedValue<T>;
    math::Vector3<T> s_teme = third_body_position_teme<T>(body, base, state.time);

    math::Vector3<T> a = third_body_accel<T>(state.position(), s_teme, body.mu);

    // Omitted nutation + equation-of-equinoxes frame fidelity: a conservative 30″
    // direction bound → relative acceleration bound, deposited into accuracy.
    TV frame_eps = math::exact<T>(30) * (math::pi<T>() / math::exact<T>(648000));
    TV amag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    using std::abs;
    T fb = abs((frame_eps * amag).value);
    a.x = math::add_bound(a.x, fb, math::ErrorChannel::accuracy);
    a.y = math::add_bound(a.y, fb, math::ErrorChannel::accuracy);
    a.z = math::add_bound(a.z, fb, math::ErrorChannel::accuracy);

    return force_wrench_from_world_accel<T>(state, a);
}

/// Build a propagator ForceFn for `body` from base epoch `base` (captures both).
/// The consumer adds it to a Propagator's force list — OPT-IN; it is NOT in the
/// default DqSgp4 force list, whose value path stays unchanged (the gate is the
/// consumer, as geopotential was gated before its DQ switch).
template<typename T>
std::function<dynamics::Wrench<T>(const dynamics::State<T>&,
                                  const constants::ConstantsProvider<T>&)>
make_third_body_force(const ThirdBody<T>& body, const astronomy::Epoch<T>& base) {
    return [body, base](const dynamics::State<T>& s,
                        const constants::ConstantsProvider<T>& K) -> dynamics::Wrench<T> {
        return third_body_force<T>(s, K, body, base);
    };
}

}  // namespace forces
