#pragma once

/// @file drag.h
/// Atmospheric-drag force lambda with an exponential-atmosphere density
/// model.
///
/// World-frame drag acceleration (per unit body mass):
///
///   a_drag = −(1/2) ρ B |v_rel| v_rel
///
/// where:
///   ρ       — atmospheric density at the body's altitude
///   B       — inverse ballistic coefficient: C_d · A / m
///   v_rel   — body velocity relative to the rotating atmosphere:
///             v_rel = v_world − (ω_E × r)
///   |v_rel| — magnitude of the relative velocity
///
/// Density comes from a PLUGGABLE `DensityModel<T>` callback (issue DRAG1):
/// altitude [m] → density [kg/m³] as a `TrackedValue<T>` that carries its OWN
/// three-error budget, including the model's accuracy. The exponential (Lane)
/// model
///
///   ρ(h) = ρ_0 · exp(−(h − h_0) / H_scale),    h = |r| − R_E
///
/// is one implementation; an NRLMSISE-00 / MSIS / tabulated model is plugged in
/// by writing another `DensityModel` and passing it to `make_drag` — the drag
/// force itself does not change. `make_drag_exponential` is the convenience that
/// wires the Lane model in.
///
/// Model truncation (REQ-EF-7 / AUD-EF-6). The density model DECLARES its own
/// fractional accuracy on the returned density (the Lane exponential defaults to
/// ~30 % vs the U.S. Standard Atmosphere), and that bound PROPAGATES through the
/// drag arithmetic into `errors.accuracy` on the acceleration — model-derived,
/// not a flat fraction stamped on the output, so a swapped-in higher-fidelity
/// model supplies a tighter bound automatically. Density is fully tracked
/// (`TrackedValue<T> exp`, issue B1), so precision/measurement propagate too.
///
/// Convention. As with `gravity_central` and `gravity_J2`, the lambda
/// returns the body-frame `Wrench<T>` with the perturbation acceleration
/// in the linear slot. Mass-independent gravitational convention applies
/// here too — `a_drag` is acceleration regardless of the body's mass;
/// `B` already incorporates `1/m`.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-6 (partial: accuracy budget on returned wrench),
///   AUD-EF-7.

#include "../constants/constants_provider.h"
#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <cmath>
#include <functional>

namespace forces {

/// A pluggable atmospheric density model (issue DRAG1): geocentric altitude
/// above the equatorial radius [m] → density [kg/m³] as a `TrackedValue<T>`
/// carrying its OWN three-error budget — crucially including the model's
/// fractional accuracy. This is the extension seam: the Lane exponential below
/// is one implementation; an NRLMSISE-00 / MSIS / tabulated model is plugged in
/// by writing another `DensityModel<T>` and passing it to `make_drag`.
template<typename T>
using DensityModel =
    std::function<math::TrackedValue<T>(const math::TrackedValue<T>& altitude_m)>;

/// Exponential (Lane 1965) atmosphere: ρ(h) = ρ_0 · exp(−(h − h_0) / H_scale),
/// fully tracked (issue B1) so precision/measurement propagate through `exp`.
/// `rel_accuracy` is the model's fractional accuracy vs a high-fidelity
/// reference; it defaults to the Lane 0.30 baseline (the model reproduces the
/// U.S. Standard Atmosphere within ~30 %, up to ~100 % at solar max — see
/// design/derivations/sgp4_near_earth_drag_theoretical_basis.md §2 Lemma 2.2,
/// §15 ERROR SOURCE A-D2). It is recorded as `errors.accuracy` on the density,
/// so the bound PROPAGATES through the drag arithmetic instead of being stamped
/// on the acceleration; a swapped-in model passes its own (tighter) value.
template<typename T>
DensityModel<T> exponential_density_model(
    math::TrackedValue<T> rho_0,
    math::TrackedValue<T> h_0,
    math::TrackedValue<T> H_scale,
    T rel_accuracy = T(3) / T(10)) {
    return [rho_0, h_0, H_scale, rel_accuracy](
        const math::TrackedValue<T>& h) -> math::TrackedValue<T> {
        using std::abs;
        math::TrackedValue<T> arg = -(h - h_0) / H_scale;
        math::TrackedValue<T> rho = rho_0 * exp(arg);
        rho.errors.accuracy = rho.errors.accuracy + rel_accuracy * abs(rho.value);
        return rho;
    };
}

/// HISTORICAL Lane-fallback stub (issue DRAG1), retained for its test callers.
/// A higher-fidelity space-weather model plugs the SAME `DensityModel<T>`
/// interface (write the model, pass it to `make_drag`; every drag consumer
/// picks it up unchanged). This stub evaluates the Lane exponential with a
/// widened 50 % accuracy band; it is a fallback, never a fake high-fidelity
/// result. New code should not call it.
template<typename T>
DensityModel<T> nrlmsise00_density_model_stub(
    math::TrackedValue<T> rho_0,
    math::TrackedValue<T> h_0,
    math::TrackedValue<T> H_scale) {
    return exponential_density_model<T>(rho_0, h_0, H_scale, T(1) / T(2));
}

/// Build a body-frame drag wrench lambda over ANY density model (issue DRAG1).
///
/// @param density  Pluggable `DensityModel<T>` (altitude → tracked density).
/// @param B        Inverse ballistic coefficient C_d · A / m [m² / kg].
///
/// Returns a `(State, ConstantsProvider) → Wrench` lambda for the propagator's
/// force list. World-frame drag a = −½ ρ B |v_rel| v_rel, with
/// v_rel = v_world − (ω_E × r). The body's mass enters through `B` (incorporates
/// 1/m); the wrench is the drag acceleration. The density model's accuracy
/// budget propagates into the acceleration's `errors.accuracy` automatically.
///
/// Lifetime. Captured parameters are stored by value in the closure.
template<typename T>
std::function<dynamics::Wrench<T>(const dynamics::State<T>&,
                                  const constants::ConstantsProvider<T>&)>
make_drag(DensityModel<T> density, math::TrackedValue<T> B) {
    return [density, B](
        const dynamics::State<T>& s,
        const constants::ConstantsProvider<T>& K) -> dynamics::Wrench<T> {

        // World-frame position and velocity.
        math::Vector3<T> r_world = s.position();
        math::Vector3<T> v_world =
            s.pose.rotation().rotate(s.twist.linear);

        // Earth rotation rate vector: along +z, magnitude ω_E.
        math::Vector3<T> omega_vec(
            math::exact<T>(0), math::exact<T>(0), K.earth.omega);

        // Atmospheric velocity at r: v_atm = ω_E × r.
        math::Vector3<T> v_atm = omega_vec.cross(r_world);
        math::Vector3<T> v_rel = v_world - v_atm;
        math::TrackedValue<T> v_rel_mag = v_rel.magnitude();

        // Geocentric altitude above the equatorial radius.
        math::TrackedValue<T> r_mag = r_world.magnitude();
        math::TrackedValue<T> h = r_mag - K.earth.a;

        // Tracked density from the pluggable model — it carries its own accuracy
        // bound, which propagates through the products below into the wrench.
        math::TrackedValue<T> rho = density(h);

        // Drag-acceleration prefactor: −(1/2) ρ B |v_rel|. Sign on the velocity
        // itself is folded in by multiplying by v_rel below.
        math::TrackedValue<T> prefactor =
            -math::ratio<T>(1, 2) * rho * B * v_rel_mag;

        math::Vector3<T> a_world(
            prefactor * v_rel.x,
            prefactor * v_rel.y,
            prefactor * v_rel.z);

        // Inverse rotation to body frame.
        math::Vector3<T> a_body =
            s.pose.rotation().conjugate().rotate(a_world);
        return dynamics::Wrench<T>(math::Vector3<T>(), a_body);
    };
}

/// Convenience: exponential-atmosphere drag (the prior public API). Wires the
/// Lane `exponential_density_model` into the generic `make_drag`.
///
/// @param rho_0    Reference density at altitude h_0 [kg / m³].
/// @param h_0      Reference altitude above R_E [m].
/// @param H_scale  Scale height [m].
/// @param B        Inverse ballistic coefficient C_d · A / m [m² / kg].
template<typename T>
std::function<dynamics::Wrench<T>(const dynamics::State<T>&,
                                  const constants::ConstantsProvider<T>&)>
make_drag_exponential(
    math::TrackedValue<T> rho_0,
    math::TrackedValue<T> h_0,
    math::TrackedValue<T> H_scale,
    math::TrackedValue<T> B) {
    return make_drag<T>(
        exponential_density_model<T>(rho_0, h_0, H_scale), B);
}

} // namespace forces
