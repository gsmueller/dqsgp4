#pragma once

/**
 * @file sgp4_propagator.h
 * @brief SGP4 Propagator — generic computation engine with injectable lambdas.
 *
 * Contains ZERO hardcoded formulas. Every computation is provided via
 * ModelFunctions lambdas from the ModelConfiguration. The propagator
 * is a pure assembly of injected computations.
 *
 * @par Architecture
 * - Propagator: combines a ModelConfiguration + TLE elements. The initialized
 *   constants live in NearSpaceInit / DeepSpaceInit (near_space.h, deep_space.h);
 *   the only mutable integrator state is the deep-space resonance.
 *
 * @par References
 * - Hoots & Roehrich (1980) Spacetrack Report #3
 * - Vallado et al. (2006) "Revisiting Spacetrack Report #3" Rev 3
 * - design/lambda_injection.md
 * - design/library_architecture.md
 */

#include "state_vector.h"
#include "model_functions.h"
#include "model_selector.h"
#include "near_space.h"
#include "deep_space.h"
#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "../tle/tle_parser.h"
#include <stdexcept>

namespace sgp4 {

/**
 * @brief The SGP4 propagator.
 *
 * Accepts a ModelConfiguration (which provides the ellipsoid, astronomical
 * constants, and injectable lambdas) plus TLE elements. Propagates to
 * a requested time and returns position/velocity with three-error budget.
 *
 * All computations are performed via the ModelFunctions lambdas.
 * The propagator NEVER calls domain module functions directly.
 *
 * @tparam T The numeric type (double, cpp_bin_float_50, etc.)
 */
template<typename T>
class Propagator {
public:
    /**
     * @brief Construct from a ModelConfiguration and TLE elements.
     *
     * @param config    Complete model configuration (ellipsoid + astro + lambdas)
     * @param elements  Parsed TLE orbital elements
     * @param tolerance Series evaluation tolerance
     */
    Propagator(
        const ModelConfiguration<T>& config,
        const tle::TleElements<T>& elements,
        const T& tolerance)
        : config_(config)
        , elements_(elements)
        , tolerance_(tolerance)
        , use_deep_space_(false)
    {
        initialize();
    }

    /**
     * @brief Propagate to time tsince minutes from epoch.
     * @param tsince_minutes  Time since epoch [minutes]
     * @return Position [km] and velocity [km/s] in TEME frame with three-error budget
     */
    StateVector<T> propagate(const math::TrackedValue<T>& tsince_minutes) const {
        if (use_deep_space_) {
            return propagate_deep_space(tsince_minutes);
        } else {
            return propagate_near_space(tsince_minutes);
        }
    }

    /// Get the model configuration description.
    const std::string& model_description() const { return config_.description; }

    /// Check if this satellite uses deep-space propagation.
    bool is_deep_space() const { return use_deep_space_; }

private:
    ModelConfiguration<T> config_;
    tle::TleElements<T> elements_;
    T tolerance_;

    // Initialized constants (computed once, one of these is used based on orbit type)
    NearSpaceInit<T> ns_init_;
    DeepSpaceInit<T> ds_init_;

    bool use_deep_space_;

    /// Initialize: element recovery, secular rates, all coefficients.
    void initialize() {
        // First do near-space init to determine orbit type
        ns_init_ = sgp4::initialize_near_space(config_, elements_, tolerance_);
        use_deep_space_ = ns_init_.is_deep_space;

        if (use_deep_space_) {
            ds_init_ = sgp4::initialize_deep_space(config_, elements_, tolerance_);
        }
    }

    /// Near-space propagation — delegates to near_space.h
    StateVector<T> propagate_near_space(const math::TrackedValue<T>& tsince) const {
        return sgp4::propagate_near_space(ns_init_, config_, tsince, tolerance_);
    }

    /// Deep-space propagation — delegates to deep_space.h.
    /// Passes ns_init_ for the near-earth drag coefficients (applied to
    /// deep-space sats too, per SR3 / reference propagation.py).
    StateVector<T> propagate_deep_space(const math::TrackedValue<T>& tsince) const {
        return sgp4::propagate_deep_space(ds_init_, ns_init_, config_, tsince, tolerance_);
    }
};

} // namespace sgp4
