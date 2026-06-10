#pragma once

/// Angle manipulation functions.
/// All inputs and outputs are TrackedValue<T> in radians.
/// Errors propagate through every operation.

#include "tracked_value.h"
#include <boost/math/constants/constants.hpp>

namespace math {

/// Pi as a TrackedValue: measurement = 0, precision = type epsilon, accuracy = 0.
template<typename T>
TrackedValue<T> pi() {
    T val = boost::math::constants::pi<T>();
    return TrackedValue<T>(val, T(0), val * std::numeric_limits<T>::epsilon(), T(0));
}

/// Two pi as a TrackedValue.
template<typename T>
TrackedValue<T> two_pi() {
    return exact<T>(2) * pi<T>();
}

/// Convert degrees to radians: rad = deg * pi / 180.
template<typename T>
TrackedValue<T> degrees_to_radians(const TrackedValue<T>& degrees) {
    return degrees * pi<T>() / exact<T>(180);
}

/// Convert radians to degrees: deg = rad * 180 / pi.
template<typename T>
TrackedValue<T> radians_to_degrees(const TrackedValue<T>& radians) {
    return radians * exact<T>(180) / pi<T>();
}

/// Wrap angle to [0, 2pi). Errors pass through (fmod is Lipschitz-1).
template<typename T>
TrackedValue<T> wrap_two_pi(const TrackedValue<T>& angle) {
    TrackedValue<T> tp = two_pi<T>();
    TrackedValue<T> result = fmod(angle, tp);
    // Ensure positive
    if (result.value < T(0)) {
        result = result + tp;
    }
    return result;
}

/// Wrap angle to [-pi, pi). Errors pass through.
template<typename T>
TrackedValue<T> wrap_neg_pos_pi(const TrackedValue<T>& angle) {
    TrackedValue<T> p = pi<T>();
    TrackedValue<T> result = wrap_two_pi(angle + p) - p;
    return result;
}

} // namespace math
