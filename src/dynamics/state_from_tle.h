#pragma once

/// @file state_from_tle.h
/// Bridge: a TLE → an initial `dynamics::State<T>` at epoch, so the
/// dual-quaternion (DQSGP4) propagator can be seeded from two-line element
/// sets.
///
/// The epoch position/velocity are recovered with the analytical SGP4
/// propagator at t = 0, then packed into the dual-quaternion state.
/// **Including this header opts into the SGP4 dependency** — the core
/// `dynamics` propagator (`dynamics/propagator.h`) does NOT include SGP4;
/// only this bridge does.
///
/// Units: SGP4 yields kilometres / km·s⁻¹ in the TEME frame; the returned
/// `State<T>` is in metres / m·s⁻¹ (the DQSGP4 convention), with identity
/// attitude and zero body angular velocity.

#include "state.h"
#include "state_conversion.h"
#include "../sgp4/sgp4_propagator.h"
#include "../sgp4/model_selector.h"
#include "../tle/tle_parser.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"
#include "../math/quaternion.h"
#include <string>

namespace dynamics {

/// Build an epoch `State<T>` from parsed TLE data.
///
/// @param td            Parsed TLE (see `tle::parse`).
/// @param tolerance     Series-evaluation tolerance for the SGP4 model.
/// @param model_preset  SGP4 model preset (default `"sgp4_standard"`, WGS72).
/// @return The satellite state at epoch — metres / m·s⁻¹, inertial (TEME),
///         identity attitude, zero spin.
template<typename T>
State<T> state_from_tle(const tle::TleData& td,
                        const T& tolerance,
                        const std::string& model_preset = "sgp4_standard") {
    sgp4::ModelConfiguration<T> config =
        sgp4::ModelSelector<T>::select(model_preset, tolerance);
    tle::TleElements<T> elements = tle::TleElements<T>::from_tle_data(td);
    sgp4::Propagator<T> prop(config, elements, tolerance);
    sgp4::StateVector<T> sv =
        prop.propagate(math::TrackedValue<T>::exact_integer(0));

    // km → m + identity-attitude pack via the shared F1 adapter (R3c) — the
    // same per-component ·exact(1000) the prior inline block performed, so the
    // delegation is BIT-IDENTICAL (the seeded suite is the regression gate).
    return to_state(sv, math::exact<T>(0));
}

}  // namespace dynamics
