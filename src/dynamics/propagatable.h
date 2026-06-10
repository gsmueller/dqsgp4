#pragma once

/// @file propagatable.h
/// The `Propagatable` concept (R3c, design/derivations/dq_propagator_facade.md
/// §5): the verb the two propagator families share. Both
///
///   sgp4::Propagator<T>::propagate(tsince_minutes)      → StateVector<T>
///   dynamics::DqSgp4Propagator<T>::propagate(tsince)    → State<T>
///
/// expose `propagate(TrackedValue<T>)` — same name, same time convention
/// (minutes since epoch), different state types. The concept names exactly that
/// shared verb, so generic code (model comparisons, sweep harnesses, the F2
/// api-parity gate's generic driver) is written once over any propagator. It
/// deliberately does NOT constrain the return type — the state-type adapters
/// (state_conversion.h, F1) are the bridge when a common currency is needed.

#include "../math/tracked_value.h"

#include <utility>

namespace dynamics {

/// A propagator: anything exposing `propagate(tsince_minutes)` for the numeric
/// type T. The return type is the propagator's native state.
template<typename P, typename T>
concept Propagatable = requires(const P& p, const math::TrackedValue<T>& t) {
    p.propagate(t);
};

}  // namespace dynamics
