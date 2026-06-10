#pragma once

/**
 * @file state_vector.h
 * @brief SGP4 propagation output: position and velocity with three-error budget.
 *
 * The StateVector holds the propagated satellite state in the TEME
 * (True Equator, Mean Equinox) reference frame. Each component of
 * position and velocity is a TrackedValue carrying measurement,
 * precision, and accuracy errors propagated from the TLE elements
 * through every operation in the SGP4 chain.
 *
 * Error decomposition methods allow callers to determine how much
 * of the total uncertainty comes from each source:
 * - Measurement: TLE element uncertainties (format precision, B* estimation)
 * - Precision: Floating-point arithmetic accumulation
 * - Accuracy: Model truncation (Brouwer J₂²+J₄, simplified drag, etc.)
 */

#include "../math/tracked_value.h"
#include "../math/vector3.h"

namespace sgp4 {

/**
 * @brief Satellite state: position [km] and velocity [km/s] in TEME frame.
 *
 * @tparam T  The underlying numeric type (double, cpp_bin_float_50, etc.)
 */
template<typename T>
struct StateVector {
    math::Vector3<T> position_km;    ///< TEME position [km]
    math::Vector3<T> velocity_km_s;  ///< TEME velocity [km/s]

    /// RSS total position error across all three components [km].
    /// Combines measurement, precision, and accuracy into a single bound.
    ///
    /// Returns plain T (a scalar magnitude) for consistency with the per-
    /// category accessors below. Spatial RSS across orthogonal Cartesian
    /// components is a rigorous Euclidean-norm bound (Pythagorean), distinct
    /// from the "RSS-as-statistical-estimate" caveat in tracked_value.h:35-36
    /// which applies only to RSS across the three error *categories*.
    T position_error() const {
        T px = position_km.x.total_error();
        T py = position_km.y.total_error();
        T pz = position_km.z.total_error();
        using std::sqrt;
        return sqrt(px * px + py * py + pz * pz);
    }

    /// RSS total velocity error across all three components [km/s].
    /// Returns plain T; see position_error() for the rigor argument.
    T velocity_error() const {
        T vx = velocity_km_s.x.total_error();
        T vy = velocity_km_s.y.total_error();
        T vz = velocity_km_s.z.total_error();
        using std::sqrt;
        return sqrt(vx * vx + vy * vy + vz * vz);
    }

    /// RSS measurement error contribution to position [km].
    /// From TLE format precision and B* estimation uncertainty.
    T position_measurement_error() const {
        T mx = position_km.x.errors.measurement;
        T my = position_km.y.errors.measurement;
        T mz = position_km.z.errors.measurement;
        using std::sqrt;
        return sqrt(mx * mx + my * my + mz * mz);
    }

    /// RSS precision error contribution to position [km].
    /// From floating-point arithmetic accumulation through ~200 operations.
    T position_precision_error() const {
        T px = position_km.x.errors.precision;
        T py = position_km.y.errors.precision;
        T pz = position_km.z.errors.precision;
        using std::sqrt;
        return sqrt(px * px + py * py + pz * pz);
    }

    /// RSS accuracy error contribution to position [km].
    /// From model truncation: Brouwer series order, simplified drag, etc.
    T position_accuracy_error() const {
        T ax = position_km.x.errors.accuracy;
        T ay = position_km.y.errors.accuracy;
        T az = position_km.z.errors.accuracy;
        using std::sqrt;
        return sqrt(ax * ax + ay * ay + az * az);
    }

    /// RSS measurement error contribution to velocity [km/s] (item F1-c,
    /// symmetric with the position accessors above).
    T velocity_measurement_error() const {
        T mx = velocity_km_s.x.errors.measurement;
        T my = velocity_km_s.y.errors.measurement;
        T mz = velocity_km_s.z.errors.measurement;
        using std::sqrt;
        return sqrt(mx * mx + my * my + mz * mz);
    }

    /// RSS precision error contribution to velocity [km/s].
    T velocity_precision_error() const {
        T px = velocity_km_s.x.errors.precision;
        T py = velocity_km_s.y.errors.precision;
        T pz = velocity_km_s.z.errors.precision;
        using std::sqrt;
        return sqrt(px * px + py * py + pz * pz);
    }

    /// RSS accuracy error contribution to velocity [km/s].
    T velocity_accuracy_error() const {
        T ax = velocity_km_s.x.errors.accuracy;
        T ay = velocity_km_s.y.errors.accuracy;
        T az = velocity_km_s.z.errors.accuracy;
        using std::sqrt;
        return sqrt(ax * ax + ay * ay + az * az);
    }
};

} // namespace sgp4
