#pragma once

/// @file ecliptic_state.h
/// A body's position in ecliptic-of-date coordinates (theory note §5,
/// design/derivations/ephemeris.md). Shared by the solar (sun_meeus.h) and lunar
/// (moon_meeus.h) instances. For the Sun `lat ≈ 0` by definition; `radius` is in
/// AU for the Sun and km for the Moon (pinned per instance).

#include "../math/tracked_value.h"

namespace ephemeris {

/// Geocentric ecliptic-of-date spherical state (λ, β, Δ) of a body.
template<typename T>
struct EclipticState {
    math::TrackedValue<T> lon;     ///< ecliptic-of-date longitude λ [rad], wrapped [0, 2π)
    math::TrackedValue<T> lat;     ///< ecliptic latitude β [rad]
    math::TrackedValue<T> radius;  ///< geocentric distance (AU for Sun, km for Moon)
};

} // namespace ephemeris
