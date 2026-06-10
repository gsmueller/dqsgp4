#pragma once

/// @file srp.h
/// Solar radiation pressure — the cannonball force with cylindrical Earth shadow
/// (replan R2; completes the L4 classical force set). Theory:
/// design/derivations/solar_radiation_pressure.md.
///
///   a_srp = ν · P₁ᴬᵁ · (AU/d)² · C_R·(A/m) · d̂ ,   d⃗ = r − s  (away from the Sun)
///
/// P₁ᴬᵁ is GENERATED from three exact-by-convention constants (never the copied
/// "4.56e-6"): P = L☉/(4π·AU²·c) with the IAU 2015 B3 nominal L☉, the IAU 2012
/// exact AU, and the SI defining c — so the value is precision-only, and the
/// nominal-vs-real TSI variability (~±0.1 % over the solar cycle) is deposited as
/// an honest representativeness accuracy band (note §2). ν ∈ {0,1} is the
/// cylindrical-umbra shadow factor (note §3; conical/penumbra refinements are
/// documented, not modeled — they bias eclipse timing by seconds, far below the
/// cannonball C_R uncertainty).
///
/// The Sun handle is the same ThirdBody<T> descriptor the gravitational
/// third-body force uses (only `.position` is consumed), and the epoch-map +
/// TEME frame step is the SHARED third_body_position_teme helper — one
/// implementation, two forces (the validity criterion applied at birth).
///
/// NOT in the default DqSgp4 force list — the gate (SRP1) is the consumer until
/// the R3 presets wire it. The SGP4 path is untouched (OR1).

#include "force_common.h"
#include "third_body.h"

#include "../astronomy/epoch.h"
#include "../constants/constants_provider.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/angles.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <cmath>
#include <functional>

namespace forces {

/// Solar radiation pressure at 1 AU, GENERATED: P = L☉/(4π·AU²·c) ≈ 4.5398e-6
/// N/m² (note §2). L☉ = IAU 2015 B3 nominal (exact by convention), AU = IAU 2012
/// (exact), c = SI defining constant (exact) — all CR1B-allowlisted `defined`.
/// The ~±0.1 % solar-cycle TSI variability about the nominal is deposited as a
/// 1e-3 relative representativeness accuracy.
template<typename T>
math::TrackedValue<T> solar_radiation_pressure_1au() {
    using TV = math::TrackedValue<T>;
    using std::abs;
    TV L_sun = TV::defined("3.828e26");    // IAU 2015 B3 nominal solar luminosity [W]
    TV c = TV::defined("299792458");       // SI defining speed of light [m/s]
    TV au = astronomical_unit<T>();
    TV four_pi = math::exact<T>(4) * math::pi<T>();
    TV P = L_sun / (four_pi * au * au * c);
    T band = abs((math::ratio<T>(1, 1000) * P).value);
    return math::add_bound(P, band, math::ErrorChannel::accuracy);
}

/// Cannonball SRP acceleration at `r_sat` from the Sun at `s_sun` (same frame,
/// metres), with the cylindrical-umbra shadow (note §1/§3). `cr_area_over_mass`
/// is the lumped C_R·(A/m) [m²/kg]; `earth_radius` the shadow-cylinder radius
/// [m]. Returns exact zero in umbra. Shadow branch selection is by value
/// comparison (the established regime-branch precedent).
template<typename T>
math::Vector3<T> srp_accel(const math::Vector3<T>& r_sat,
                           const math::Vector3<T>& s_sun,
                           const math::TrackedValue<T>& cr_area_over_mass,
                           const math::TrackedValue<T>& earth_radius) {
    using TV = math::TrackedValue<T>;

    TV s2 = s_sun.x * s_sun.x + s_sun.y * s_sun.y + s_sun.z * s_sun.z;
    TV r_dot_s = r_sat.x * s_sun.x + r_sat.y * s_sun.y + r_sat.z * s_sun.z;

    if (r_dot_s.value < T(0)) {
        // Night side of the terminator plane: in umbra iff inside the cylinder,
        // |r⊥|² = |r|² − (r·ŝ)² < R_E²  (evaluated in values; ν is boolean).
        T r2 = (r_sat.x * r_sat.x + r_sat.y * r_sat.y + r_sat.z * r_sat.z).value;
        T along = r_dot_s.value * r_dot_s.value / s2.value;   // (r·ŝ)²
        T rperp2 = r2 - along;
        T re2 = earth_radius.value * earth_radius.value;
        if (rperp2 < re2) {
            return math::Vector3<T>();   // ν = 0: exact-zero acceleration
        }
    }

    math::Vector3<T> d(r_sat.x - s_sun.x, r_sat.y - s_sun.y, r_sat.z - s_sun.z);
    TV d2 = d.x * d.x + d.y * d.y + d.z * d.z;
    TV d_mag = sqrt(d2);
    TV au = astronomical_unit<T>();
    TV P = solar_radiation_pressure_1au<T>() * (au * au / d2);
    TV coef = P * cr_area_over_mass / d_mag;   // scalar; ·d⃗ gives the d̂ direction
    return math::Vector3<T>(coef * d.x, coef * d.y, coef * d.z);
}

/// Body-frame SRP wrench for the propagator force list (note §4/§5): the shared
/// epoch-map + TEME precession step (third_body_position_teme), the cannonball
/// acceleration, the omitted-nutation frame bound (~30″, as the third-body
/// force), and the shared body-frame epilogue. `K` supplies the shadow radius.
template<typename T>
dynamics::Wrench<T> srp_force(const dynamics::State<T>& state,
                              const constants::ConstantsProvider<T>& K,
                              const ThirdBody<T>& sun,
                              const astronomy::Epoch<T>& base,
                              const math::TrackedValue<T>& cr_area_over_mass) {
    using TV = math::TrackedValue<T>;
    using std::abs;
    math::Vector3<T> s_teme = third_body_position_teme<T>(sun, base, state.time);
    math::Vector3<T> a =
        srp_accel<T>(state.position(), s_teme, cr_area_over_mass, K.earth.a);

    TV frame_eps = math::exact<T>(30) * (math::pi<T>() / math::exact<T>(648000));
    TV amag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    T fb = abs((frame_eps * amag).value);
    a.x = math::add_bound(a.x, fb, math::ErrorChannel::accuracy);
    a.y = math::add_bound(a.y, fb, math::ErrorChannel::accuracy);
    a.z = math::add_bound(a.z, fb, math::ErrorChannel::accuracy);

    return force_wrench_from_world_accel<T>(state, a);
}

/// Build a propagator ForceFn for SRP (captures the Sun descriptor, base epoch,
/// and the lumped C_R·A/m). OPT-IN — not in the default DqSgp4 force list.
template<typename T>
std::function<dynamics::Wrench<T>(const dynamics::State<T>&,
                                  const constants::ConstantsProvider<T>&)>
make_srp_force(const ThirdBody<T>& sun, const astronomy::Epoch<T>& base,
               const math::TrackedValue<T>& cr_area_over_mass) {
    return [sun, base, cr_area_over_mass](
               const dynamics::State<T>& s,
               const constants::ConstantsProvider<T>& K) -> dynamics::Wrench<T> {
        return srp_force<T>(s, K, sun, base, cr_area_over_mass);
    };
}

}  // namespace forces
