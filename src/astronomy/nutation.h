#pragma once

/**
 * @file nutation.h
 * @brief IAU 2000A luni-solar nutation — Δψ, Δε series, the N(t) matrix, and the
 *        equation of the equinoxes (replan R4b).
 *
 * The periodic response of Earth's spin axis to the lunisolar torque on the
 * equatorial bulge: a Poisson series in the five Delaunay arguments
 * (l, l′, F, D, Ω), dominated by the 18.6-year −17.2″ sin Ω nodal term.
 * Theory: design/derivations/nutation_itrs.md.
 *
 * Model honesty (no perceived fidelity):
 *  - the LUNI-SOLAR IAU 2000A series (678 terms, exact-integer amplitudes in
 *    1e-7″, auto-generated from the in-repo SOFA table) with the dial-up
 *    `n_terms` truncation TRACKED as Σ omitted |amplitude| — the
 *    truncated-series-with-tracked-tail pattern's third instance;
 *  - the PLANETARY series (687 terms) is NOT modeled: its measured total
 *    amplitude (Δψ 5.54 mas, Δε 1.90 mas) is deposited as an accuracy floor;
 *  - gate NUT1 measures the angles against erfa.nut06a (the COMPLETE model) —
 *    the residual is the real planetary signal, majorized by the floors.
 *
 * Fundamental arguments: the five IERS quartic polynomials [arcsec],
 * transcribed verbatim from the in-repo AstroLib fundarg (e06 branch).
 */

#include "frames.h"                    // rot_x/rot_z, Matrix3 (SOFA-equivalent, FRAME2)
#include "nutation_terms_iau2000a.h"   // auto-generated series + planetary floors
#include "../math/angles.h"
#include "../math/tracked_polynomial.h"
#include "../math/tracked_value.h"

#include <array>
#include <cmath>
#include <vector>

namespace astronomy {

/// Nutation in longitude and obliquity [rad] with their tracked budgets.
template<typename T>
struct NutationAngles {
    math::TrackedValue<T> dpsi;
    math::TrackedValue<T> deps;
};

namespace detail {

/// Quartic fundamental-argument polynomial [arcsec] from born-digital strings
/// (AstroLib fundarg verbatim) — the fw_poly pattern at degree 4.
template<typename T>
math::TrackedPolynomial<T> delaunay_poly(const char* const c[5], const double m[5]) {
    using TV = math::TrackedValue<T>;
    return math::TrackedPolynomial<T>(
        std::vector<TV>{TV::model_coefficient(c[0]), TV::model_coefficient(c[1]),
                        TV::model_coefficient(c[2]), TV::model_coefficient(c[3]),
                        TV::model_coefficient(c[4])},
        std::vector<T>{T(m[0]), T(m[1]), T(m[2]), T(m[3]), T(m[4])});
}

} // namespace detail

/// The five Delaunay arguments {l, l′, F, D, Ω} at t [TT Julian centuries from
/// J2000], radians wrapped to [0, 2π). IERS/Simon-1994 quartics, AstroLib-verbatim.
template<typename T>
std::array<math::TrackedValue<T>, 5> delaunay_arguments(const math::TrackedValue<T>& t) {
    using TV = math::TrackedValue<T>;
    static const char* kL[5]  = {"485868.249036", "1717915923.2178", "31.8792",
                                 "0.051635", "-0.00024470"};
    static const char* kLp[5] = {"1287104.793048", "129596581.0481", "-0.5532",
                                 "0.000136", "-0.00001149"};
    static const char* kF[5]  = {"335779.526232", "1739527262.8478", "-12.7512",
                                 "-0.001037", "0.00000417"};
    static const char* kD[5]  = {"1072260.703692", "1602961601.2090", "-6.3706",
                                 "0.006593", "-0.00003169"};
    static const char* kOm[5] = {"450160.398036", "-6962890.5431", "7.4722",
                                 "0.007702", "-0.00005939"};
    static const double mL[5]  = {485868.249036, 1717915923.2178, 31.8792,
                                  0.051635, 0.00024470};
    static const double mLp[5] = {1287104.793048, 129596581.0481, 0.5532,
                                  0.000136, 0.00001149};
    static const double mF[5]  = {335779.526232, 1739527262.8478, 12.7512,
                                  0.001037, 0.00000417};
    static const double mD[5]  = {1072260.703692, 1602961601.2090, 6.3706,
                                  0.006593, 0.00003169};
    static const double mOm[5] = {450160.398036, 6962890.5431, 7.4722,
                                  0.007702, 0.00005939};
    static const math::TrackedPolynomial<T> pL  = detail::delaunay_poly<T>(kL, mL);
    static const math::TrackedPolynomial<T> pLp = detail::delaunay_poly<T>(kLp, mLp);
    static const math::TrackedPolynomial<T> pF  = detail::delaunay_poly<T>(kF, mF);
    static const math::TrackedPolynomial<T> pD  = detail::delaunay_poly<T>(kD, mD);
    static const math::TrackedPolynomial<T> pOm = detail::delaunay_poly<T>(kOm, mOm);

    const math::ErrorChannel acc = math::ErrorChannel::accuracy;
    TV a2r = math::pi<T>() / math::exact<T>(648000);
    return {math::wrap_two_pi(pL.eval(t, 5, acc, 1) * a2r),
            math::wrap_two_pi(pLp.eval(t, 5, acc, 1) * a2r),
            math::wrap_two_pi(pF.eval(t, 5, acc, 1) * a2r),
            math::wrap_two_pi(pD.eval(t, 5, acc, 1) * a2r),
            math::wrap_two_pi(pOm.eval(t, 5, acc, 1) * a2r)};
}

/// Δψ, Δε [rad] from the luni-solar IAU 2000A series at t [TT centuries].
/// `n_terms` dials the amplitude-sorted series (1..678); the omitted luni-solar
/// tail AND the unmodeled planetary floors are deposited in errors.accuracy.
template<typename T>
NutationAngles<T> nutation_angles(const math::TrackedValue<T>& t,
                                  int n_terms = kNutationTermCount) {
    using TV = math::TrackedValue<T>;
    using std::abs;

    std::array<TV, 5> args = delaunay_arguments<T>(t);
    int n = n_terms;
    if (n < 1) n = 1;
    if (n > kNutationTermCount) n = kNutationTermCount;

    TV psi = math::exact<T>(0);
    TV eps = math::exact<T>(0);
    for (int i = 0; i < n; ++i) {
        const NutationTerm& tm = kNutationIau2000a[i];
        TV theta = math::exact<T>(tm.l) * args[0] + math::exact<T>(tm.lp) * args[1]
                 + math::exact<T>(tm.f) * args[2] + math::exact<T>(tm.d) * args[3]
                 + math::exact<T>(tm.om) * args[4];
        TV s = sin(theta);
        TV c = cos(theta);
        psi = psi + (math::exact<T>(int(tm.psi_sin)) +
                     math::exact<T>(int(tm.psi_t_sin)) * t) * s
                  + math::exact<T>(int(tm.psi_cos)) * c;
        eps = eps + (math::exact<T>(int(tm.eps_cos)) +
                     math::exact<T>(int(tm.eps_t_cos)) * t) * c
                  + math::exact<T>(int(tm.eps_sin)) * s;
    }

    // Unit: 1e-7 arcsec → rad (exact rationals × the pi generator).
    TV unit = math::pi<T>() / math::exact<T>(648000) / math::exact<T>(10000000);
    psi = psi * unit;
    eps = eps * unit;

    // Truncation tail (omitted luni-solar terms) + the unmodeled planetary
    // floors [1e-7 arcsec], converted with the same unit (plain T bound).
    T at = abs(t.value);
    T tail_psi = T(0);
    T tail_eps = T(0);
    for (int i = n; i < kNutationTermCount; ++i) {
        const NutationTerm& tm = kNutationIau2000a[i];
        tail_psi += T(abs(tm.psi_sin)) + T(abs(tm.psi_t_sin)) * at + T(abs(tm.psi_cos));
        tail_eps += T(abs(tm.eps_cos)) + T(abs(tm.eps_t_cos)) * at + T(abs(tm.eps_sin));
    }
    tail_psi += T(kPlanetaryPsiFloor);
    tail_eps += T(kPlanetaryEpsFloor);
    T u = abs(unit.value);
    psi = math::add_bound(psi, tail_psi * u, math::ErrorChannel::accuracy);
    eps = math::add_bound(eps, tail_eps * u, math::ErrorChannel::accuracy);
    return NutationAngles<T>{psi, eps};
}

/// The nutation matrix N(t): mean-of-date → true-of-date,
/// N = R1(−(ε_A+Δε)) · R3(−Δψ) · R1(ε_A) (the erfa `numat` composition; our
/// rot_x/rot_z are SOFA's R1/R3, proved bit-exact in FRAME2).
template<typename T>
math::Matrix3<T> nutation_matrix(const math::TrackedValue<T>& eps_a,
                                 const NutationAngles<T>& ang) {
    return rot_x(-(eps_a + ang.deps)) * rot_z(-ang.dpsi) * rot_x(eps_a);
}

/// Equation of the equinoxes [rad]: Δψ·cos ε_A + the two classical complementary
/// terms (+0.00264″ sin Ω + 0.000063″ sin 2Ω, in-repo AstroLib verbatim). The
/// remaining ~31 eect00 complementary terms (~30 μas total) are bounded by a
/// 50 μas accuracy floor (note §3).
template<typename T>
math::TrackedValue<T> equation_of_equinoxes(const NutationAngles<T>& ang,
                                            const math::TrackedValue<T>& eps_a,
                                            const math::TrackedValue<T>& omega) {
    using TV = math::TrackedValue<T>;
    using std::abs;
    TV a2r = math::pi<T>() / math::exact<T>(648000);
    TV ee = ang.dpsi * cos(eps_a)
          + TV::model_coefficient("0.00264") * a2r * sin(omega)
          + TV::model_coefficient("0.000063") * a2r * sin(math::exact<T>(2) * omega);
    T floor = abs((TV::model_coefficient("0.00005") * a2r).value);
    return math::add_bound(ee, floor, math::ErrorChannel::accuracy);
}

} // namespace astronomy
