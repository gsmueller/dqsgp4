#pragma once

/// @file dq_sgp4_propagator.h
/// API-parity bridge: a TLE-driven dual-quaternion (DQSGP4) propagator that
/// exposes the SAME verb as the analytical SGP4 propagator —
/// `propagate(tsince_minutes)` — so the two models differ only in their
/// mathematics, not their interface (REQ-PR / F2). Compare:
///
///   SGP4:    sgp4::Propagator<T> p(config, elems, tol);
///            sgp4::StateVector<T> sv = p.propagate(t_min);
///
///   DQSGP4:  auto p = dynamics::DqSgp4Propagator<T>::authentic(td, tol, dt);
///            dynamics::State<T> s = p.propagate(t_min);
///
/// **Including this header opts into the SGP4 dependency** (the epoch state is
/// recovered through the analytical model — see `state_from_tle.h`).
///
/// Two model layers are EXPLICIT (F2-c — no silent mismatch): the epoch state is
/// seeded from the AUTHENTIC WGS72 SGP4 (a TLE is WGS72 by definition), while the
/// DQ propagation runs under a chosen earth model + zonal field. The two named
/// factories surface this as the authentic-vs-boosted mode duality (G1):
///   - `authentic` — WGS72 ellipsoid + WGS72 zonals (J₂…J₄), the DQ analogue of
///     the `sgp4_standard` model;
///   - `boosted`   — WGS84 ellipsoid + EGM2008 zonals (J₂…J₉), higher fidelity,
///     still seeded from the authentic WGS72 SGP4.
///
/// Time convention (F2-d): the verb takes minutes-since-epoch (the SGP4
/// convention); internally the DQ clock is seconds, so minutes are scaled by 60.
/// Propagation is forward-only (tsince ≤ 0 returns the epoch state).
///
/// The default force model is gravitational (central monopole + zonal Jₙ),
/// BIT-UNCHANGED. Perturbations are layered on two ways (R3,
/// design/derivations/dq_propagator_facade.md): arbitrary `extra_forces`
/// appended via the explicit constructor, or the named `DqForceOptions` presets
/// on `from_tle` — lunisolar third-body (TB1), Vallado-8-4 table drag (ATM1),
/// and cannonball SRP (SRP1), each opt-in and assembled from the TLE epoch.

#include "inertia.h"
#include "propagator.h"
#include "state.h"
#include "state_from_tle.h"

#include "../astronomy/epoch.h"
#include "../atmosphere/exponential_table.h"
#include "../constants/constants_provider.h"
#include "../constants/gravity_field.h"
#include "../constants/tesseral_harmonics.h"
#include "../constants/zonal_harmonics.h"
#include "../forces/drag.h"
#include "../forces/geopotential.h"
#include "../forces/srp.h"
#include "../forces/third_body.h"
#include "../integrators/runge_kutta.h"
#include "../math/tracked_value.h"
#include "../tle/tle_parser.h"

#include <optional>
#include <utility>
#include <vector>

namespace dynamics {

/// Propagator fidelity mode (G1): the authentic frozen SGP4 model vs a boosted,
/// higher-fidelity one. The two named factories below realise each; `from_tle`
/// takes the mode as a first-class parameter so callers select it by name.
enum class PropagatorMode {
    Authentic,  ///< WGS72 ellipsoid + WGS72 zonals (J₂…J₄) — the sgp4_standard analogue.
    Boosted,    ///< WGS84 ellipsoid + EGM2008 zonals (J₂…J₉) — higher fidelity.
};

/// Opt-in perturbation additions for the DQ force model (R3 presets,
/// dq_propagator_facade.md §3). Default-constructed = none — the gravitational-
/// only default model, bit-unchanged. The additions compose (any subset);
/// `drag_B` (C_d·A/m) and `srp_cr_area_over_mass` (C_R·A/m) are caller-supplied
/// because they are physical properties of the spacecraft — no honest default
/// exists. The Sun/Moon forces' base epoch is the TLE epoch (UTC; the TT−UTC
/// ≈ 69 s offset moves the Moon ≤ ~40″ / Sun ≤ ~10″ — ≤ ~2e-4 relative on
/// perturbations that are ~1e-7 of central gravity, below every propagation
/// tolerance; a real UTC→TT conversion arrives with the L1 leap-second work).
template<typename T>
struct DqForceOptions {
    bool lunisolar = false;                                  ///< Sun + Moon third-body.
    std::optional<math::TrackedValue<T>> drag_B;             ///< Vallado-8-4 table drag.
    std::optional<math::TrackedValue<T>> srp_cr_area_over_mass;  ///< Cannonball SRP.
};

/// TLE-seeded dual-quaternion propagator with a SGP4-symmetric interface.
///
/// @tparam T  Underlying numeric type.
template<typename T>
class DqSgp4Propagator {
public:
    /// Explicit constructor: an epoch state, the DQ propagation earth model and
    /// zonal field, the zonal truncation degree, the max integrator step, and
    /// optional extra force lambdas APPENDED after the gravitational core
    /// (wrench-additive, REQ-PR-2). Default-empty `extra_forces` is bit-identical
    /// to the prior facade. Prefer the `authentic` / `boosted` / `from_tle`
    /// factories for the common cases.
    DqSgp4Propagator(State<T> epoch_state,
                     constants::ConstantsProvider<T> dq_constants,
                     constants::ZonalHarmonics<T> zonals,
                     int zonal_degree,
                     math::TrackedValue<T> dt_max,
                     std::vector<ForceFn<T>> extra_forces = {})
        : epoch_(std::move(epoch_state)),
          dt_max_(std::move(dt_max)),
          prop_(make_propagator(std::move(dq_constants), std::move(zonals),
                                zonal_degree, std::move(extra_forces))) {}

    /// State at `tsince_minutes` after epoch — the SGP4-symmetric verb. Minutes
    /// (SGP4 convention) are converted to the DQ second-clock internally.
    /// Forward-only: tsince ≤ 0 returns the epoch state.
    State<T> propagate(const math::TrackedValue<T>& tsince_minutes) const {
        math::TrackedValue<T> t_target =
            epoch_.time + tsince_minutes * math::exact<T>(60);
        if (t_target.value <= epoch_.time.value) {
            return epoch_;
        }
        return prop_.propagate_to(epoch_, t_target, dt_max_);
    }

    /// The seed state at epoch (metres / m·s⁻¹, inertial, identity attitude).
    const State<T>& epoch_state() const { return epoch_; }

    /// The underlying DQ propagator (constants, forces, integrator).
    const Propagator<T>& propagator() const { return prop_; }

    // --- Named-mode factories (authentic vs boosted — G1) -------------------

    /// Authentic mode: WGS72 ellipsoid + WGS72 zonals (J₂…J₄), seeded from the
    /// authentic WGS72 SGP4. The DQ analogue of `sgp4_standard`.
    static DqSgp4Propagator authentic(const tle::TleData& td, const T& tolerance,
                                      const math::TrackedValue<T>& dt_max) {
        State<T> seed = state_from_tle<T>(td, tolerance, "sgp4_standard");
        return DqSgp4Propagator(std::move(seed),
                                constants::ConstantsProvider<T>::wgs72(tolerance),
                                constants::ZonalHarmonics<T>::wgs72(tolerance),
                                4, dt_max);
    }

    /// Boosted mode: WGS84 ellipsoid + EGM2008 zonals (J₂…J₉), seeded from the
    /// authentic WGS72 SGP4 (TLEs are WGS72). Higher-fidelity propagation.
    static DqSgp4Propagator boosted(const tle::TleData& td, const T& tolerance,
                                    const math::TrackedValue<T>& dt_max) {
        State<T> seed = state_from_tle<T>(td, tolerance, "sgp4_standard");
        return DqSgp4Propagator(std::move(seed),
                                constants::ConstantsProvider<T>::wgs84(tolerance),
                                constants::ZonalHarmonics<T>::egm2008(tolerance),
                                9, dt_max);
    }

    /// Build for an explicit fidelity mode — the mode as a named parameter (G1) —
    /// plus opt-in perturbation presets (R3). Default options reproduce the
    /// gravitational-only model BIT-IDENTICALLY (gate FM1).
    static DqSgp4Propagator from_tle(const tle::TleData& td, const T& tolerance,
                                     PropagatorMode mode,
                                     const math::TrackedValue<T>& dt_max,
                                     const DqForceOptions<T>& options = {}) {
        std::vector<ForceFn<T>> extra = forces_from_options(options, td);
        State<T> seed = state_from_tle<T>(td, tolerance, "sgp4_standard");
        if (mode == PropagatorMode::Boosted) {
            return DqSgp4Propagator(std::move(seed),
                                    constants::ConstantsProvider<T>::wgs84(tolerance),
                                    constants::ZonalHarmonics<T>::egm2008(tolerance),
                                    9, dt_max, std::move(extra));
        }
        return DqSgp4Propagator(std::move(seed),
                                constants::ConstantsProvider<T>::wgs72(tolerance),
                                constants::ZonalHarmonics<T>::wgs72(tolerance),
                                4, dt_max, std::move(extra));
    }

private:
    /// Assemble the preset perturbation forces from the options (facade note §3).
    /// The Sun/Moon base epoch is the TLE epoch (UTC — the documented ≈TT
    /// treatment; see DqForceOptions).
    static std::vector<ForceFn<T>> forces_from_options(const DqForceOptions<T>& opt,
                                                       const tle::TleData& td) {
        std::vector<ForceFn<T>> fs;
        if (opt.drag_B) {
            fs.push_back(forces::make_drag<T>(
                atmosphere::vallado84_density_model<T>(), *opt.drag_B));
        }
        if (opt.lunisolar || opt.srp_cr_area_over_mass) {
            tle::TleElements<T> el = tle::TleElements<T>::from_tle_data(td);
            astronomy::Epoch<T> base = astronomy::Epoch<T>::from_jd(
                el.epoch_jd, astronomy::TimeScale::UTC);
            if (opt.lunisolar) {
                fs.push_back(forces::make_third_body_force<T>(
                    forces::sun_third_body<T>(), base));
                fs.push_back(forces::make_third_body_force<T>(
                    forces::moon_third_body<T>(), base));
            }
            if (opt.srp_cr_area_over_mass) {
                fs.push_back(forces::make_srp_force<T>(
                    forces::sun_third_body<T>(), base, *opt.srp_cr_area_over_mass));
            }
        }
        return fs;
    }

    /// Assemble the DQ propagator with a gravitational force model (the unified
    /// geopotential: monopole + zonal Jₙ up to `zonal_degree` in one Cunningham
    /// V/W pass), any extra force lambdas APPENDED, and an RK4 integrator.
    static Propagator<T> make_propagator(constants::ConstantsProvider<T> K,
                                         constants::ZonalHarmonics<T> zonals,
                                         int zonal_degree,
                                         std::vector<ForceFn<T>> extra_forces) {
        // Unified geopotential: monopole + zonal (order m = 0) in ONE Cunningham
        // V/W pass, replacing the separate gravity_central + gravity_zonal
        // lambdas. Round-off-equal to the prior sum (GEOPOT oracle 1, 4.9e-16);
        // rotation is linear, so rotate(mono) + rotate(zonal) = rotate(mono+zonal).
        // The field is zonal-only (empty tesseral), so the result is longitude-
        // independent and gmst-invariant — gmst is passed as 0 (R = identity).
        constants::GravityField<T> field(std::move(zonals),
                                         constants::TesseralHarmonics<T>{});
        std::vector<ForceFn<T>> forces;
        forces.push_back(
            [field, zonal_degree](const State<T>& s,
                                  const constants::ConstantsProvider<T>& KK) {
                return forces::geopotential(s, KK, field, zonal_degree, 0,
                                            math::exact<T>(0));
            });
        for (ForceFn<T>& f : extra_forces) {
            forces.push_back(std::move(f));
        }

        IntegratorFn<T> integrator =
            [](const State<T>& y0, const math::TrackedValue<T>& dt,
               const integrators::AccelFn<T>& accel) {
                return integrators::runge_kutta_4(y0, dt, accel);
            };

        return Propagator<T>(std::move(K),
                             Inertia<T>::point_mass(math::exact<T>(1)),
                             std::move(forces), std::move(integrator));
    }

    State<T> epoch_;
    math::TrackedValue<T> dt_max_;
    Propagator<T> prop_;
};

} // namespace dynamics
