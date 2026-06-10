#pragma once

/// @file geopotential.h
/// The unified geopotential acceleration — monopole + zonal + tesseral in ONE
/// Cunningham V/W pass (design/derivations/geopotential.md §3). Subsumes
/// gravity_central / gravity_zonal / gravity_tesseral / gravity_J2 on the runtime
/// path: max_m = 0 selects zonal-only; (max_n, max_m) = (2, 0) is J₂; and the
/// monopole C(0,0) = 1 carries the point mass, so −GM·r/r³ falls out of the SAME
/// loop (no separate central path). NOT bit-identical to the summed
/// gravity_central + gravity_zonal forces (Cunningham V/W vs the Legendre-in-u
/// recurrence for the zonal half) — gated at ROUND-OFF by test_geopotential
/// against that sum + a closed-form J₂ + the monopole identity, never sold as
/// bit-exact (no perceived fidelity). The superseded per-force headers are
/// retained as the GEOPOT verification oracle.
///
/// The frozen SGP4 inline J₂/J₃/J₄ and sgp4::ZonalHarmonics are NOT routed here.

#include "force_common.h"

#include "../astronomy/frames.h"   // L2: sidereal_rotation (TEME→ECEF, Rz(GMST))
#include "../constants/constants_provider.h"
#include "../constants/gravity_field.h"
#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/matrix3.h"
#include "../math/spherical_harmonics.h"   // shared Cunningham V/W recursion
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <cmath>
#include <vector>

namespace forces {

/// Geopotential acceleration in the EARTH-FIXED frame at r_ecef. One V/W pass +
/// the Montenbruck-Gill (3.33) acceleration over 0 ≤ m ≤ n ≤ N with the m = 0
/// (zonal/monopole) and m ≥ 1 (tesseral) cases of geopotential.md §3. max_m caps
/// the order (max_m = 0 ⇒ zonal-only). Returns acceleration in the units implied
/// by mu and Re; the omitted-degree Kaula tail is booked in the accuracy channel.
template<typename T>
math::Vector3<T> geopotential_accel_ecef(
    const math::Vector3<T>& r_ecef,
    const math::TrackedValue<T>& mu,
    const math::TrackedValue<T>& Re,
    const constants::GravityField<T>& field,
    int max_n, int max_m) {
    using math::TrackedValue;
    using math::exact;
    using TV = TrackedValue<T>;

    const int fd = field.max_degree();
    const int N = (max_n < fd) ? max_n : fd;
    const int M = (max_m < N) ? max_m : N;
    const int top = N + 1;  // V/W needed up to degree+order N+1

    math::VWTable<T> vw = math::cunningham_vw(r_ecef.x, r_ecef.y, r_ecef.z, Re, top);
    std::vector<std::vector<TV>>& V = vw.V;
    std::vector<std::vector<TV>>& W = vw.W;

    // Acceleration (Montenbruck-Gill 3.33), all orders in one loop.
    TV scale = mu / (Re * Re);
    TV ax = exact<T>(0), ay = exact<T>(0), az = exact<T>(0);
    for (int n = 0; n <= N; ++n) {
        const int mhi = (n < M) ? n : M;
        for (int m = 0; m <= mhi; ++m) {
            TV C = field.C(n, m);
            TV S = field.S(n, m);
            if (m == 0) {
                // m = 0 (zonal + monopole): the order −1 terms of the ½-formula
                // fold back by the negative-order reflection to one term (§3 (3)).
                ax = ax + (-C * V[n + 1][1]);
                ay = ay + (-C * W[n + 1][1]);
            } else {
                TV f = exact<T>((n - m + 2) * (n - m + 1));  // (n−m+2)!/(n−m)!
                ax = ax + math::ratio<T>(1, 2) * (
                    (-C * V[n + 1][m + 1] - S * W[n + 1][m + 1])
                    + f * (C * V[n + 1][m - 1] + S * W[n + 1][m - 1]));
                ay = ay + math::ratio<T>(1, 2) * (
                    (-C * W[n + 1][m + 1] + S * V[n + 1][m + 1])
                    + f * (-C * W[n + 1][m - 1] + S * V[n + 1][m - 1]));
            }
            az = az + exact<T>(n - m + 1) * (-C * V[n + 1][m] - S * W[n + 1][m]);
        }
    }

    math::Vector3<T> a(scale * ax, scale * ay, scale * az);

    // Kaula truncation tail → accuracy (geopotential.md §4): omitted degrees
    // n > N have |C|,|S| ~ kaula/n² and (R/r)^{n+1} decays geometrically, so the
    // omitted acceleration is bounded by the leading term inflated by the
    // geometric sum 1/(1 − R/r) — a conservative, T-independent fidelity floor.
    TV r_mag = sqrt(r_ecef.x * r_ecef.x + r_ecef.y * r_ecef.y
                   + r_ecef.z * r_ecef.z);
    TV ror = Re / r_mag;  // R/r < 1 for exterior orbits
    TV rp = exact<T>(1);
    for (int k = 0; k < N + 2; ++k) rp = rp * ror;  // (R/r)^{N+2}
    TV kaula = TV::model_coefficient("1e-5");        // Kaula's rule-of-thumb amplitude
    TV denom = exact<T>((N + 1) * (N + 1)) * (exact<T>(1) - ror);
    TV tail_tv = scale * exact<T>(N + 2) * rp * kaula / denom;
    using std::abs;
    T tail = abs(tail_tv.value);
    a.x = math::add_bound(a.x, tail, math::ErrorChannel::accuracy);
    a.y = math::add_bound(a.y, tail, math::ErrorChannel::accuracy);
    a.z = math::add_bound(a.z, tail, math::ErrorChannel::accuracy);
    return a;
}

/// Body-frame geopotential wrench for the propagator force list. Rotates the ECI
/// (TEME) position to ECEF by R = Rz(GMST), evaluates the geopotential there,
/// rotates the acceleration back by Rᵀ, then into the body frame via the shared
/// force_wrench_from_world_accel epilogue.
template<typename T>
dynamics::Wrench<T> geopotential(
    const dynamics::State<T>& state,
    const constants::ConstantsProvider<T>& K,
    const constants::GravityField<T>& field,
    int max_n, int max_m,
    const math::TrackedValue<T>& gmst) {
    math::Vector3<T> r_eci = state.position();
    math::Matrix3<T> R = astronomy::sidereal_rotation(gmst);
    math::Vector3<T> r_ecef = R * r_eci;
    math::Vector3<T> a_ecef = geopotential_accel_ecef<T>(
        r_ecef, K.earth.GM, K.earth.a, field, max_n, max_m);
    math::Vector3<T> a_eci = R.transpose() * a_ecef;
    return force_wrench_from_world_accel<T>(state, a_eci);
}

}  // namespace forces
