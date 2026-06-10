#pragma once

/// @file gravity_zonal.h
/// Zonal-harmonic gravity perturbations: closed-form `gravity_J2` and the
/// general `gravity_zonal` (J₂…J_max_n via Legendre recursion). Superseded on
/// the runtime path by the unified forces/geopotential.h; RETAINED as the
/// independent GEOPOT verification oracle (different arithmetic: Legendre-in-u
/// vs Cunningham V/W) with direct test callers.
///
/// The Earth's gravity field expressed as a zonal expansion:
///
///   V(r, φ) = −(μ/r) · [1 − Σ_{n≥2} J_n (R_E/r)^n P_n(sin φ)]
///
/// The central −μ/r term is computed by `gravity_central`. Each zonal
/// J_n contribution is a separate force lambda so the caller composes
/// the model they need: J₂-only (Brouwer-secular / SGP4-style),
/// J₂+J₃ (asymmetric N–S), J₂+J₃+J₄ (typical LEO), or the full stack.
///
/// Convention. All lambdas follow the same convention as
/// `gravity_central`: they return a body-frame `Wrench<T>` whose linear
/// component is the perturbation acceleration. Gravitational physics is
/// mass-independent at this level, so the lambda returns acceleration
/// regardless of the body's mass. The propagator divides only the
/// angular slot by the body's principal moments; the linear slot is
/// pass-through (correct for the mass-independent gravity convention,
/// and exact when the body's mass equals 1 in the inertia).
///
/// Model truncation (REQ-EF-7 / AUD-EF-6). The residual relative to the
/// full geopotential is the sum of all higher zonals not included, and it
/// IS recorded as an upper bound in `errors.accuracy`: `gravity_J2` bounds
/// the leading omitted J₃ by Kaula's a-priori magnitude (below), and
/// `gravity_zonal` sums the real omitted tail (table terms beyond `max_n`)
/// plus a Kaula bound for the beyond-table remainder (see `res_x/y/z`).
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-8,
///   AUD-CC-9, AUD-CC-10, AUD-CC-12, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../constants/constants_provider.h"
#include "../constants/zonal_harmonics.h"
#include "../dynamics/pose.h"
#include "../dynamics/state.h"
#include "../dynamics/wrench.h"
#include "../math/quaternion.h"
#include "../math/tracked_value.h"
#include "../math/vector3.h"

#include <cmath>

namespace forces {

namespace detail {

/// Kaula's rule-of-thumb a-priori magnitude for an UNMODELLED zonal Jₙ:
///   |Jₙ| ≈ 1e-5 · √(2n+1) / n²
/// (Kaula 1966, "Theory of Satellite Geodesy", §3). Used ONLY to bound the
/// truncation residual of omitted higher zonals when no measured coefficient is
/// available; never as a modelled value.
template<typename T>
T kaula_zonal_magnitude(int n) {
    using std::sqrt;
    const T kaula_rms = T(1) / T(100000);  // 1e-5 normalized-coefficient RMS
    return kaula_rms * sqrt(T(2 * n + 1)) / T(n * n);
}

}  // namespace detail

/// J₂ (oblateness) perturbation acceleration in the body frame.
///
/// World-frame J₂ acceleration (per unit mass) derived from the J₂ term
/// of the geopotential expansion:
///
///   a_x = (3/2) μ J₂ R_E² · x · (5 z²/r² − 1) / r⁵
///   a_y = (3/2) μ J₂ R_E² · y · (5 z²/r² − 1) / r⁵
///   a_z = (3/2) μ J₂ R_E² · z · (5 z²/r² − 3) / r⁵
///
/// where r = state.pose.translation(), R_E = K.earth.a, J₂ = K.earth.J2n(1).
/// The body-frame force is obtained by inverse rotation through
/// state.pose.rotation().
///
/// Closed analytical result: J₂ produces no torque on a point body. The
/// returned wrench has zero in the torque slot.
template<typename T>
dynamics::Wrench<T> gravity_J2(const dynamics::State<T>& state,
                                const constants::ConstantsProvider<T>& K) {
    math::Vector3<T> r = state.position();
    math::TrackedValue<T> r_sq = r.x*r.x + r.y*r.y + r.z*r.z;
    math::TrackedValue<T> r_mag = sqrt(r_sq);
    math::TrackedValue<T> r_5 = r_sq * r_sq * r_mag;

    math::TrackedValue<T> z_sq = r.z * r.z;
    math::TrackedValue<T> five_z2_over_r2 = math::exact<T>(5) * z_sq / r_sq;
    math::TrackedValue<T> xy_factor = five_z2_over_r2 - math::exact<T>(1);
    math::TrackedValue<T> z_factor  = five_z2_over_r2 - math::exact<T>(3);

    math::TrackedValue<T> J2 = K.earth.J2n(1);
    math::TrackedValue<T> Re = K.earth.a;
    math::TrackedValue<T> mu = K.earth.GM;
    math::TrackedValue<T> common =
        math::ratio<T>(3, 2) * mu * J2 * Re * Re / r_5;

    math::Vector3<T> a_world(
        common * xy_factor * r.x,
        common * xy_factor * r.y,
        common * z_factor  * r.z
    );

    // REQ-EF-7: model-truncation bound for the omitted J₃+ zonals. This J₂-only
    // force has no zonal table, so it bounds the leading omitted term J₃ by
    // Kaula's a-priori magnitude (|J₃| ≈ 1e-5·√7/9 ≈ 2.9e-6 — conservative vs the
    // measured 2.5e-6): |a_J3| ≤ |J₃/J₂|·(a_E/r)·|a_J2|·safety. For a residual
    // built from the REAL higher coefficients, use `gravity_zonal` with a
    // ZonalHarmonics provider. (Supersedes the former |J₄/J₂| proxy, whose
    // 1.5e-3·1.1 ≈ 1.6e-3 actually UNDER-bounded the real |J₃/J₂| ≈ 2.3e-3.)
    {
        using std::abs;
        T j3_apriori = detail::kaula_zonal_magnitude<T>(3);
        T J3_over_J2 = j3_apriori / abs(J2.value);
        T aE_over_r = K.earth.a.value / r_mag.value;
        T safety = T(11) / T(10);
        T trunc_factor = safety * J3_over_J2 * aE_over_r;
        a_world.x.errors.accuracy =
            a_world.x.errors.accuracy + trunc_factor * abs(a_world.x.value);
        a_world.y.errors.accuracy =
            a_world.y.errors.accuracy + trunc_factor * abs(a_world.y.value);
        a_world.z.errors.accuracy =
            a_world.z.errors.accuracy + trunc_factor * abs(a_world.z.value);
    }

    // Inverse rotation to body frame: a_body = q_r* · a_world · q_r.
    math::Vector3<T> a_body =
        state.pose.rotation().conjugate().rotate(a_world);
    return dynamics::Wrench<T>(math::Vector3<T>(), a_body);
}

/// Generic Jₙ zonal perturbation — the sum of the J₂…J_{max_n} contributions,
/// each the closed-form gradient of the zonal geopotential term
///
///   a_n = μ Jₙ Rₑⁿ r^{-(n+3)} · [ (n+1) Pₙ(u) + u Pₙ'(u) ] · (x, y, ·),
///
/// with u = sin φ = z/r and the z-component carrying the (x²+y²) structure.
/// The Legendre polynomial and its derivative come from the non-singular
/// recurrence
///   Pₙ = [(2n−1) u Pₙ₋₁ − (n−1) Pₙ₋₂] / n,    Pₙ' = n Pₙ₋₁ + u Pₙ₋₁'.
///
/// The coefficients Jₙ are read from the shared generative `ZonalHarmonics`
/// provider — the common-library home for the Jₖ — so even AND odd zonals come
/// from ONE standard-tagged source (no |J₄/J₂| proxy); μ and Rₑ come from the
/// ellipsoid. The truncation residual (zonals beyond max_n) is added to
/// `errors.accuracy`: the leading omitted term is built from the table's real
/// next coefficient when present, else from Kaula's a-priori magnitude, with a
/// tail safety factor and the |Pₘ'| ≤ m(m+1)/2 Legendre envelope.
///
/// For max_n = 2 this reproduces `gravity_J2`; J₃ introduces the north–south
/// asymmetry, J₄ the next even correction, and so on.
template<typename T>
dynamics::Wrench<T> gravity_zonal(const dynamics::State<T>& state,
                                   const constants::ConstantsProvider<T>& K,
                                   const constants::ZonalHarmonics<T>& zh,
                                   int max_n) {
    using math::TrackedValue;
    using math::exact;
    using std::abs;

    math::Vector3<T> r = state.position();
    TrackedValue<T> x = r.x;
    TrackedValue<T> y = r.y;
    TrackedValue<T> z = r.z;
    TrackedValue<T> r_sq  = x * x + y * y + z * z;
    TrackedValue<T> r_mag = sqrt(r_sq);
    TrackedValue<T> inv_r = exact<T>(1) / r_mag;
    TrackedValue<T> inv_r3 = inv_r / r_sq;       // 1/r³
    TrackedValue<T> u = z * inv_r;               // sin φ
    TrackedValue<T> x2y2 = x * x + y * y;

    TrackedValue<T> mu = K.earth.GM;
    TrackedValue<T> Re = K.earth.a;
    TrackedValue<T> Re_over_r = Re * inv_r;

    // The table degree is the modelling ceiling; max_n selects how many of those
    // zonals enter the FORCE. Terms above max_n (but in the table) are summed in
    // magnitude as the truncation residual, so the residual is the real omitted
    // tail and shrinks monotonically as max_n grows.
    const int table_max = zh.max_degree();
    const int force_max = (max_n < table_max) ? max_n : table_max;

    // Legendre recurrence seeds: P₀, P₁ and P₁' (P₀' = 0 is unused).
    TrackedValue<T> P_nm2 = exact<T>(1);   // P₀
    TrackedValue<T> P_nm1 = u;             // P₁
    TrackedValue<T> dP_nm1 = exact<T>(1);  // P₁'
    TrackedValue<T> Re_r_pow = Re_over_r;  // (Re/r)¹ → (Re/r)ⁿ within the loop

    math::Vector3<T> a_world;        // zero-initialized force accumulator
    T res_x = T(0), res_y = T(0), res_z = T(0);  // omitted-term magnitude sums

    for (int n = 2; n <= table_max; ++n) {
        Re_r_pow = Re_r_pow * Re_over_r;   // (Re/r)ⁿ
        TrackedValue<T> P_n =
            (exact<T>(2 * n - 1) * u * P_nm1 - exact<T>(n - 1) * P_nm2) / exact<T>(n);
        TrackedValue<T> dP_n = exact<T>(n) * P_nm1 + u * dP_nm1;

        TrackedValue<T> Jn = zh.Jn(n);
        TrackedValue<T> common = mu * Jn * Re_r_pow * inv_r3;   // μ Jₙ Rₑⁿ / r^{n+3}
        TrackedValue<T> xy_factor = exact<T>(n + 1) * P_n + u * dP_n;

        TrackedValue<T> ax = common * xy_factor * x;
        TrackedValue<T> ay = common * xy_factor * y;
        TrackedValue<T> az = common * (exact<T>(n + 1) * z * P_n - dP_n * x2y2 * inv_r);

        if (n <= force_max) {
            a_world.x = a_world.x + ax;
            a_world.y = a_world.y + ay;
            a_world.z = a_world.z + az;
        } else {
            res_x = res_x + abs(ax.value);   // real omitted-zonal magnitude
            res_y = res_y + abs(ay.value);
            res_z = res_z + abs(az.value);
        }

        P_nm2 = P_nm1;  P_nm1 = P_n;
        dP_nm1 = dP_n;
    }

    // Beyond the table, add a Kaula a-priori bound for the (table_max+1) term,
    // ×2 for the geometrically smaller infinite tail past it; |Pₘ| ≤ 1,
    // |Pₘ'| ≤ m(m+1)/2 give the (m+1)(m+2)/2 envelope. Re_r_pow is (Re/r)^{table_max}.
    {
        const int m = table_max + 1;
        const T jm = detail::kaula_zonal_magnitude<T>(m);
        const T Re_r_m = Re_r_pow.value * Re_over_r.value;   // (Re/r)^{m}
        const T envelope = T((m + 1) * (m + 2)) / T(2);
        const T tail = T(2) * mu.value * jm * Re_r_m
                     * inv_r3.value * r_mag.value * envelope;
        res_x = res_x + tail;
        res_y = res_y + tail;
        res_z = res_z + tail;
    }

    a_world.x.errors.accuracy = a_world.x.errors.accuracy + res_x;
    a_world.y.errors.accuracy = a_world.y.errors.accuracy + res_y;
    a_world.z.errors.accuracy = a_world.z.errors.accuracy + res_z;

    math::Vector3<T> a_body =
        state.pose.rotation().conjugate().rotate(a_world);
    return dynamics::Wrench<T>(math::Vector3<T>(), a_body);
}

} // namespace forces
