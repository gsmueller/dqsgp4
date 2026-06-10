#pragma once

/**
 * @file resonance.h
 * @brief Orbital resonance detection and tesseral harmonic modeling (SR3 §6).
 *
 * Detects whether a satellite's orbit is in resonance with Earth's rotation
 * (24h synchronous or 12h half-day) and integrates the resonance equations
 * of motion via the SR3 leapfrog scheme at a fixed 720-minute step.
 *
 * @par Theory
 * For a resonant satellite the tesseral harmonics of the geopotential do not
 * average to zero over an orbital revolution. The mean longitude λ and mean
 * motion n satisfy a forced pendulum-like ODE:
 *
 *     dλ/dt = n + xfact
 *     dn/dt = Σ_k d_k · sin(arg_k(λ))
 *
 * where d_k are tesseral driving coefficients (Kaula F·G products with the
 * J_{ℓm} tesseral amplitudes) and arg_k are linear combinations of λ, ω,
 * and the constant phases G22, G32, G44, G52, G54 (rad).
 *
 * For 24h synchronous orbits SR3 uses three commensurabilities (m=1,2,3),
 * encoded as δ_1, δ_2, δ_3. For 12h half-day orbits (Molniya/GPS/GLONASS)
 * SR3 uses ten d-coefficients spanning m=2,3,4,5.
 *
 * @par References
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, §6 (DSINIT/DSPER)
 * - Hujsak (1979), AIAA 79-136 "A Restricted Four Body Solution for
 *   Resonating Satellites Without Drag"
 * - Kaula (1966), "Theory of Satellite Geodesy", §3.3 (F_{ℓmp})
 * - dnwrnr libsgp4 SGP4.cc lines 676-1340 (reference implementation)
 * - design/derivations/018_resonance.md
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include <cmath>
#include <functional>

namespace perturbation {

/**
 * @brief Resonance classification.
 *
 * Bands match SR3 §6 mean-motion gates:
 *   SYNCHRONOUS:  0.0034906585 < n0 < 0.0052359877 rad/min (~1 sidereal day)
 *   HALF_DAY:     0.00826 ≤ n0 ≤ 0.00924 rad/min AND e0 ≥ 0.5 (Molniya regime)
 */
enum class ResonanceType {
    NONE = 0,
    SYNCHRONOUS = 1,
    HALF_DAY = 2
};

/**
 * @brief Resonance detection + integration state (SR3 §6).
 *
 * Stores both the SR3 d-coefficients (precomputed) and the integrator
 * variables (xli, xni, atime — mutable across propagation calls).
 */
template<typename T>
struct ResonanceState {
    ResonanceType type = ResonanceType::NONE;

    // Stored initial mean motion for restart logic.
    math::TrackedValue<T> n0;

    // Integrator initial conditions.
    math::TrackedValue<T> xlamo;
    math::TrackedValue<T> xfact;

    // 24h SYNCHRONOUS resonance coefficients (m=1,2,3).
    math::TrackedValue<T> del1, del2, del3;

    // 12h HALF_DAY resonance coefficients.
    math::TrackedValue<T> d2201, d2211;
    math::TrackedValue<T> d3210, d3222;
    math::TrackedValue<T> d4410, d4422;
    math::TrackedValue<T> d5220, d5232, d5421, d5433;

    // Integrator state (mutable during propagation).
    mutable math::TrackedValue<T> xli;
    mutable math::TrackedValue<T> xni;
    mutable math::TrackedValue<T> atime;
};

/**
 * @brief Detect resonance type from recovered mean motion + eccentricity.
 *
 * Uses the SR3 §6 mean-motion gates directly (rather than period bands).
 */
inline ResonanceType detect_resonance_from_n(double n0_rad_min, double e0) {
    if (n0_rad_min > 0.0034906585 && n0_rad_min < 0.0052359877) {
        return ResonanceType::SYNCHRONOUS;
    }
    if (n0_rad_min >= 8.26e-3 && n0_rad_min <= 9.24e-3 && e0 >= 0.5) {
        return ResonanceType::HALF_DAY;
    }
    return ResonanceType::NONE;
}

/**
 * @brief Period-signature entry point (the original interface), kept for
 *        source compatibility; prefer detect_resonance_from_n.
 */
inline ResonanceType detect_resonance(double period_min) {
    if (period_min > 1200.0 && period_min < 1800.0) {
        return ResonanceType::SYNCHRONOUS;
    }
    if (period_min > 600.0 && period_min < 800.0) {
        return ResonanceType::HALF_DAY;
    }
    return ResonanceType::NONE;
}

/**
 * @brief Initialize resonance state per SR3 §6.
 *
 * Mirrors dnwrnr SGP4.cc DeepSpaceInitialise lines 889-1080.
 */
template<typename T>
ResonanceState<T> initialize_resonance(
    const math::TrackedValue<T>& n0,
    const math::TrackedValue<T>& e0,
    const math::TrackedValue<T>& i0,
    const math::TrackedValue<T>& omega0,
    const math::TrackedValue<T>& Omega0,
    const math::TrackedValue<T>& M0,
    const math::TrackedValue<T>& gmst,
    const math::TrackedValue<T>& xmdot,
    const math::TrackedValue<T>& omgdot,
    const math::TrackedValue<T>& xnodot,
    const math::TrackedValue<T>& ssl,
    const math::TrackedValue<T>& ssg,
    const math::TrackedValue<T>& ssh,
    const math::TrackedValue<T>& theta_dot)
{
    using math::exact;
    using math::ratio;

    ResonanceState<T> rs;
    rs.type = detect_resonance_from_n(double(n0.value), double(e0.value));
    rs.n0 = n0;

    if (rs.type == ResonanceType::NONE) {
        return rs;
    }

    // 1/a in Earth radii. SR3: a = (xke/n)^(2/3) so aqnv = (n0/xke)^(2/3).
    // Uses dnwrnr xke = 0.0743669161331734132 (rad/min, WGS72 convention).
    double n0d = double(n0.value);
    double xke_d = 0.0743669161331734132;
    double aqnv_d = std::pow(n0d / xke_d, 2.0 / 3.0);
    math::TrackedValue<T> aqnv = math::TrackedValue<T>(T(aqnv_d), n0.errors);

    math::TrackedValue<T> eosq = e0 * e0;
    math::TrackedValue<T> sinio = sin(i0);
    math::TrackedValue<T> cosio = cos(i0);
    math::TrackedValue<T> theta2 = cosio * cosio;
    math::TrackedValue<T> xpidot = omgdot + xnodot;

    double bfact_d = 0.0;

    if (rs.type == ResonanceType::SYNCHRONOUS) {
        // 24h resonance — dnwrnr SGP4.cc lines 897-924.
        math::TrackedValue<T> g200 = exact<T>(1) + eosq * (ratio<T>(-5,2) + ratio<T>(13,16) * eosq);
        math::TrackedValue<T> g310 = exact<T>(1) + exact<T>(2) * eosq;
        // 6.60937 is an SR3 synchronous-resonance FIT coefficient (the g300
        // polynomial, dnwrnr/SGP4.cc), not a convention -> model_coefficient
        // (digits -> model-fidelity accuracy, binary storage -> T-scaling
        // precision), per the constants initiative (each model coefficient is a
        // truncated fit, never defined()=0). Value-preserving (same .value), so
        // OR1 stays bit-exact; g300's accuracy now propagates into del3.
        math::TrackedValue<T> g300 = exact<T>(1) + eosq * (exact<T>(-6)
                  + math::TrackedValue<T>::model_coefficient("6.60937") * eosq);

        math::TrackedValue<T> one_plus_cosio = exact<T>(1) + cosio;
        math::TrackedValue<T> f220 = ratio<T>(3,4) * one_plus_cosio * one_plus_cosio;
        math::TrackedValue<T> f311 = ratio<T>(15,16) * sinio * sinio * (exact<T>(1) + exact<T>(3) * cosio)
                  - ratio<T>(3,4) * one_plus_cosio;
        math::TrackedValue<T> f330 = ratio<T>(15,8) * one_plus_cosio * one_plus_cosio * one_plus_cosio;

        math::TrackedValue<T> Q22 = math::TrackedValue<T>::measured("1.7891679e-6", "1e-12");
        math::TrackedValue<T> Q31 = math::TrackedValue<T>::measured("2.1460748e-6", "1e-12");
        math::TrackedValue<T> Q33 = math::TrackedValue<T>::measured("2.2123015e-7", "1e-13");

        math::TrackedValue<T> del1_base = exact<T>(3) * n0 * n0 * aqnv * aqnv;
        rs.del2 = exact<T>(2) * del1_base * f220 * g200 * Q22;
        rs.del3 = exact<T>(3) * del1_base * f330 * g300 * Q33 * aqnv;
        rs.del1 = del1_base * f311 * g310 * Q31 * aqnv;

        rs.xlamo = math::wrap_two_pi(M0 + Omega0 + omega0 - gmst);
        bfact_d = double((xmdot + xpidot - theta_dot + ssl + ssg + ssh).value);

        rs.d2201 = exact<T>(0); rs.d2211 = exact<T>(0);
        rs.d3210 = exact<T>(0); rs.d3222 = exact<T>(0);
        rs.d4410 = exact<T>(0); rs.d4422 = exact<T>(0);
        rs.d5220 = exact<T>(0); rs.d5232 = exact<T>(0);
        rs.d5421 = exact<T>(0); rs.d5433 = exact<T>(0);
    }
    else {
        // 12h Molniya resonance — dnwrnr SGP4.cc lines 932-1068.
        double e_d = double(e0.value);
        std::function<math::TrackedValue<T>(double, double, double, double)> poly_cubic = [&](double a, double b, double c, double d) {
            double v = a + e_d * (b + e_d * (c + e_d * d));
            return math::TrackedValue<T>(T(v), e0.errors);
        };

        double g201_d = -0.306 - (e_d - 0.64) * 0.440;
        math::TrackedValue<T> g201 = math::TrackedValue<T>(T(g201_d), e0.errors);

        math::TrackedValue<T> g211, g310, g322, g410, g422, g520;
        if (e_d <= 0.65) {
            g211 = poly_cubic(3.616, -13.247, 16.290, 0.0);
            g310 = poly_cubic(-19.302, 117.390, -228.419, 156.591);
            g322 = poly_cubic(-18.9068, 109.7927, -214.6334, 146.5816);
            g410 = poly_cubic(-41.122, 242.694, -471.094, 313.953);
            g422 = poly_cubic(-146.407, 841.880, -1629.014, 1083.435);
            g520 = poly_cubic(-532.114, 3017.977, -5740.032, 3708.276);
        } else {
            g211 = poly_cubic(-72.099, 331.819, -508.738, 266.724);
            g310 = poly_cubic(-346.844, 1582.851, -2415.925, 1246.113);
            g322 = poly_cubic(-342.585, 1554.908, -2366.899, 1215.972);
            g410 = poly_cubic(-1052.797, 4758.686, -7193.992, 3651.957);
            g422 = poly_cubic(-3581.69, 16178.11, -24462.77, 12422.52);
            if (e_d <= 0.715) {
                g520 = poly_cubic(1464.74, -4664.75, 3763.64, 0.0);
            } else {
                g520 = poly_cubic(-5149.66, 29936.92, -54087.36, 31324.56);
            }
        }
        math::TrackedValue<T> g533, g521, g532;
        if (e_d < 0.7) {
            g533 = poly_cubic(-919.2277, 4988.61, -9064.77, 5542.21);
            g521 = poly_cubic(-822.71072, 4568.6173, -8491.4146, 5337.524);
            g532 = poly_cubic(-853.666, 4690.25, -8624.77, 5341.4);
        } else {
            g533 = poly_cubic(-37995.78, 161616.52, -229838.2, 109377.94);
            g521 = poly_cubic(-51752.104, 218913.95, -309468.16, 146349.42);
            g532 = poly_cubic(-40023.88, 170470.89, -242699.48, 115605.82);
        }

        math::TrackedValue<T> sini2 = sinio * sinio;
        math::TrackedValue<T> f220 = ratio<T>(3,4) * (exact<T>(1) + exact<T>(2) * cosio + theta2);
        math::TrackedValue<T> f221 = ratio<T>(3,2) * sini2;
        math::TrackedValue<T> f321 = ratio<T>(15,8) * sinio
                  * (exact<T>(1) - exact<T>(2) * cosio - exact<T>(3) * theta2);
        math::TrackedValue<T> f322 = -ratio<T>(15,8) * sinio
                  * (exact<T>(1) + exact<T>(2) * cosio - exact<T>(3) * theta2);
        math::TrackedValue<T> f441 = exact<T>(35) * sini2 * f220;
        math::TrackedValue<T> f442 = ratio<T>(315,8) * sini2 * sini2;
        math::TrackedValue<T> f522 = ratio<T>(315,32) * sinio * (
              sini2 * (exact<T>(1) - exact<T>(2) * cosio - exact<T>(5) * theta2)
            + ratio<T>(1,3) * (exact<T>(-2) + exact<T>(4) * cosio + exact<T>(6) * theta2)
        );
        math::TrackedValue<T> f523 = sinio * (
              ratio<T>(315,64) * sini2 * (exact<T>(-2) - exact<T>(4) * cosio + exact<T>(10) * theta2)
            + ratio<T>(105,16) * (exact<T>(1) + exact<T>(2) * cosio - exact<T>(3) * theta2)
        );
        math::TrackedValue<T> f542 = ratio<T>(945,32) * sinio * (
            exact<T>(2) - exact<T>(8) * cosio
            + theta2 * (exact<T>(-12) + exact<T>(8) * cosio + exact<T>(10) * theta2)
        );
        math::TrackedValue<T> f543 = ratio<T>(945,32) * sinio * (
            exact<T>(-2) - exact<T>(8) * cosio
            + theta2 * (exact<T>(12) + exact<T>(8) * cosio - exact<T>(10) * theta2)
        );

        math::TrackedValue<T> xno2 = n0 * n0;
        math::TrackedValue<T> ainv2 = aqnv * aqnv;
        math::TrackedValue<T> temp1 = exact<T>(3) * xno2 * ainv2;

        math::TrackedValue<T> ROOT22 = math::TrackedValue<T>::measured("1.7891679e-6", "1e-12");
        math::TrackedValue<T> ROOT32 = math::TrackedValue<T>::measured("3.7393792e-7", "1e-13");
        math::TrackedValue<T> ROOT44 = math::TrackedValue<T>::measured("7.3636953e-9", "1e-15");
        math::TrackedValue<T> ROOT52 = math::TrackedValue<T>::measured("1.1428639e-7", "1e-13");
        math::TrackedValue<T> ROOT54 = math::TrackedValue<T>::measured("2.1765803e-9", "1e-15");

        {
            math::TrackedValue<T> temp = temp1 * ROOT22;
            rs.d2201 = temp * f220 * g201;
            rs.d2211 = temp * f221 * g211;
        }
        temp1 = temp1 * aqnv;
        {
            math::TrackedValue<T> temp = temp1 * ROOT32;
            rs.d3210 = temp * f321 * g310;
            rs.d3222 = temp * f322 * g322;
        }
        temp1 = temp1 * aqnv;
        {
            math::TrackedValue<T> temp = exact<T>(2) * temp1 * ROOT44;
            rs.d4410 = temp * f441 * g410;
            rs.d4422 = temp * f442 * g422;
        }
        temp1 = temp1 * aqnv;
        {
            math::TrackedValue<T> temp = temp1 * ROOT52;
            rs.d5220 = temp * f522 * g520;
            rs.d5232 = temp * f523 * g532;
        }
        {
            math::TrackedValue<T> temp = exact<T>(2) * temp1 * ROOT54;
            rs.d5421 = temp * f542 * g521;
            rs.d5433 = temp * f543 * g533;
        }

        rs.xlamo = math::wrap_two_pi(M0 + exact<T>(2) * Omega0 - exact<T>(2) * gmst);
        bfact_d = double((xmdot + exact<T>(2) * xnodot - exact<T>(2) * theta_dot
                          + ssl + exact<T>(2) * ssh).value);

        rs.del1 = exact<T>(0); rs.del2 = exact<T>(0); rs.del3 = exact<T>(0);
    }

    rs.xfact = math::TrackedValue<T>(T(bfact_d), n0.errors) - n0;

    rs.xli = rs.xlamo;
    rs.xni = n0;
    rs.atime = exact<T>(0);

    return rs;
}

/**
 * @brief Integrate resonance to time tsince via SR3 leapfrog (720-min step).
 *
 * Mirrors dnwrnr DeepSpaceSecular lines 1241-1339.
 */
template<typename T>
void step_resonance(
    const ResonanceState<T>& rs,
    const math::TrackedValue<T>& omega0,
    const math::TrackedValue<T>& omgdot,
    const math::TrackedValue<T>& tsince,
    const math::TrackedValue<T>& theta_dot,
    const math::TrackedValue<T>& gsto,
    const math::TrackedValue<T>& Omega_target,
    const math::TrackedValue<T>& omgasm_target,
    math::TrackedValue<T>& xll_out,
    math::TrackedValue<T>& xn_out)
{
    using TV = math::TrackedValue<T>;
    using math::exact;
    using math::ratio;
    if (rs.type == ResonanceType::NONE) return;

    // SR3 §6 resonance phase constants (rad) and the fixed 720-min leapfrog
    // step. Honest decimal-truncation bounds (CR1); each .value is bit-identical
    // to the original double literal, so the leapfrog value path is preserved
    // while the integrated mean longitude/motion now carry a budget (DS1.4).
    const TV G22   = TV::model_coefficient("5.7686396");
    const TV G32   = TV::model_coefficient("0.95240898");
    const TV G44   = TV::model_coefficient("1.8014998");
    const TV G52   = TV::model_coefficient("1.0508330");
    const TV G54   = TV::model_coefficient("4.4108898");
    const TV FASX2 = TV::model_coefficient("0.13130908");
    const TV FASX4 = TV::model_coefficient("2.8843198");
    const TV FASX6 = TV::model_coefficient("0.37448087");
    const TV STEP  = exact<T>(720);
    const TV STEP2 = STEP * STEP / exact<T>(2);

    TV xli   = rs.xli;
    TV xni   = rs.xni;
    TV atime = rs.atime;
    const TV& xlamo = rs.xlamo;
    const TV& xfact = rs.xfact;
    const TV& n0    = rs.n0;

    // Restart at epoch (dnwrnr lines 1252-1263). Integrator-control decisions
    // taken on the value (matches the reference step count exactly).
    if (std::fabs(double(tsince.value)) < 720.0
        || double(tsince.value) * double(atime.value) <= 0.0
        || std::fabs(double(tsince.value)) < std::fabs(double(atime.value)))
    {
        atime = exact<T>(0);
        xni = n0;
        xli = xlamo;
    }

    bool synchronous = (rs.type == ResonanceType::SYNCHRONOUS);
    bool running = true;
    while (running) {
        TV xndot, xnddt;
        if (synchronous) {
            xndot = rs.del1 * sin(xli - FASX2)
                  + rs.del2 * sin(exact<T>(2) * (xli - FASX4))
                  + rs.del3 * sin(exact<T>(3) * (xli - FASX6));
            xnddt = rs.del1 * cos(xli - FASX2)
                  + exact<T>(2) * rs.del2 * cos(exact<T>(2) * (xli - FASX4))
                  + exact<T>(3) * rs.del3 * cos(exact<T>(3) * (xli - FASX6));
        } else {
            TV xomi = omega0 + omgdot * atime;
            TV x2omi = xomi + xomi;
            TV x2li = xli + xli;
            xndot = rs.d2201 * sin(x2omi + xli - G22)
                  + rs.d2211 * sin(xli - G22)
                  + rs.d3210 * sin(xomi + xli - G32)
                  + rs.d3222 * sin(-xomi + xli - G32)
                  + rs.d4410 * sin(x2omi + x2li - G44)
                  + rs.d4422 * sin(x2li - G44)
                  + rs.d5220 * sin(xomi + xli - G52)
                  + rs.d5232 * sin(-xomi + xli - G52)
                  + rs.d5421 * sin(xomi + x2li - G54)
                  + rs.d5433 * sin(-xomi + x2li - G54);
            xnddt = rs.d2201 * cos(x2omi + xli - G22)
                  + rs.d2211 * cos(xli - G22)
                  + rs.d3210 * cos(xomi + xli - G32)
                  + rs.d3222 * cos(-xomi + xli - G32)
                  + rs.d5220 * cos(xomi + xli - G52)
                  + rs.d5232 * cos(-xomi + xli - G52)
                  + exact<T>(2) * (rs.d4410 * cos(x2omi + x2li - G44)
                         + rs.d4422 * cos(x2li - G44)
                         + rs.d5421 * cos(xomi + x2li - G54)
                         + rs.d5433 * cos(-xomi + x2li - G54));
        }
        TV xldot = xni + xfact;
        xnddt = xnddt * xldot;

        TV ft = tsince - atime;
        if (std::fabs(double(ft.value)) >= 720.0) {
            TV delt = (double(ft.value) >= 0.0) ? STEP : -STEP;
            // Parens match the reference `xli += xldot*delt + xndot*STEP2`
            // associativity (xli + (A + B)), not (xli + A) + B — bit-exactness.
            xli   = xli   + (xldot * delt + xndot * STEP2);
            xni   = xni   + (xndot * delt + xnddt * STEP2);
            atime = atime + delt;
        } else {
            xn_out = xni + xndot * ft + xnddt * ft * ft * ratio<T>(1, 2);
            TV xl_temp = xli + xldot * ft + xndot * ft * ft * ratio<T>(1, 2);
            TV theta = fmod(gsto + tsince * theta_dot, math::two_pi<T>());
            if (theta.value < T(0)) theta = theta + math::two_pi<T>();
            if (synchronous) {
                xll_out = xl_temp + theta - Omega_target - omgasm_target;
            } else {
                xll_out = xl_temp + exact<T>(2) * (theta - Omega_target);
            }
            running = false;
        }
    }

    rs.xli   = xli;
    rs.xni   = xni;
    rs.atime = atime;
}

} // namespace perturbation
