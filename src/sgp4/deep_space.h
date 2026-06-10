#pragma once

/**
 * @file deep_space.h
 * @brief SGP4 deep-space initialization and propagation (SR3 §6).
 *
 * Implements the SDP4 deep-space propagation model: orbital period
 * ≥ 225 minutes (GEO, HEO, Molniya, GPS, etc.). Adds:
 *   1. **DPPER** — sun + moon long-period periodic corrections to
 *      e, i, ω, Ω, M (computed each call, depends on time-varying
 *      sun/moon mean anomaly).
 *   2. **DSPER** — tesseral resonance integration of (λ, n) via
 *      a 720-minute leapfrog (24h synchronous or 12h Molniya).
 *   3. **Time-varying secular rates** — sun and moon contribute
 *      constant linear-rate corrections ssl, ssg, ssh, sse, ssi to
 *      dM/dt, dω/dt, dΩ/dt, de/dt, di/dt that are evaluated once
 *      at epoch using sun/moon mean anomalies at epoch.
 *
 * The DSCOM coefficient buildup and DPPER/DSPER algorithms are ported
 * from dnwrnr libsgp4 SGP4.cc (lines 676-1340), which is the canonical
 * reference SR3 §6 implementation.
 *
 * @par References
 * - Hoots & Roehrich (1980) Spacetrack Report No. 3, §6 (pages 13-16)
 * - Hujsak (1979) AIAA 79-136 (resonance theory)
 * - Lane & Cranford (1969) AIAA 69-925 (third-body theory)
 * - dnwrnr libsgp4 SGP4.cc (reference implementation)
 */

#include "../math/tracked_value.h"
#include "../math/angles.h"
#include "../astronomy/epoch.h"   // L1: days_since_1900 view (single-sourced J1900 base)
#include "../perturbation/resonance.h"
#include "../perturbation/short_period.h"
#include "../orbit/element_recovery.h"
#include "../orbit/osculating_elements.h"
#include "../orbit/state_from_elements.h"
#include "model_selector.h"
#include "near_space.h"   // NearSpaceInit: the near-earth drag coefficients (applied to deep-space too)
#include "state_vector.h"
#include <cmath>
#include <limits>

namespace sgp4 {

/**
 * @brief Deep-space initialization constants (SR3 §6).
 *
 * The published elements (a0, n0, etc.) and the Brouwer secular rates are
 * TrackedValue<T> so three-error budgets propagate through the secular
 * advance. The SR3 sun/moon coefficients are stored as doubles since the
 * dnwrnr algorithm itself is expressed in IEEE doubles and the measurement
 * uncertainties dominate by orders of magnitude.
 */
template<typename T>
struct DeepSpaceInit {
    // Recovered Brouwer elements.
    math::TrackedValue<T> a0, n0, e0, i0, omega0, Omega0, M0, bstar;
    math::TrackedValue<T> perigee_km, period_min;
    math::TrackedValue<T> beta0;

    // Earth constants.
    math::TrackedValue<T> re_km, xke, CK2;
    math::TrackedValue<T> J3, A30;

    // Brouwer secular rates (no third-body contributions).
    math::TrackedValue<T> M_dot, omega_dot, Omega_dot;

    // SR3 sun/moon long-period linear-rate contributions. TrackedValue<T> so
    // the secular advance carries a budget that scales with T (DS1).
    math::TrackedValue<T> sse, ssi, ssl, ssg, ssh;

    // DPPER sun/moon periodic coefficients (DSCOM block).
    math::TrackedValue<T> se2, si2, sl2, sgh2, sh2;
    math::TrackedValue<T> se3, si3, sl3, sgh3, sh3;
    math::TrackedValue<T> sl4, sgh4;
    math::TrackedValue<T> ee2, e3, xi2, xi3;
    math::TrackedValue<T> xl2, xl3, xl4;
    math::TrackedValue<T> xgh2, xgh3, xgh4;
    math::TrackedValue<T> xh2, xh3;
    // Phase variables (sun and moon mean anomalies at epoch).
    math::TrackedValue<T> zmos, zmol;

    // Resonance state (mutable across propagate calls).
    mutable perturbation::ResonanceState<T> resonance;

    // Trig constants at epoch.
    math::TrackedValue<T> cosio, sinio, theta2;
    math::TrackedValue<T> x3thm1, x1mth2, x7thm1;

    // Sidereal time at epoch and earth rotation rate.
    math::TrackedValue<T> gsto;
    math::TrackedValue<T> theta_dot;

    // A3 / CK2 ratio for inclination recomputation.
    math::TrackedValue<T> a3ovk2;
};

namespace detail {

/**
 * @brief DSCOM sun/moon coefficient buildup (dnwrnr SGP4.cc lines 687-887).
 *
 * Populates DeepSpaceInit.sse/ssi/ssl/ssg/ssh and the 22 DPPER coefficients
 * (se2..sh3 solar, ee2..xh3 lunar) plus the phase variables zmos, zmol.
 */
template<typename T>
void build_dpper_coefficients(
    DeepSpaceInit<T>& ds,
    const math::TrackedValue<T>& epoch_jd_j1900,
    const math::TrackedValue<T>& eosq,
    const math::TrackedValue<T>& sinio,
    const math::TrackedValue<T>& cosio,
    const math::TrackedValue<T>& betao,
    const math::TrackedValue<T>& betao2,
    const math::TrackedValue<T>& eccen,
    const math::TrackedValue<T>& argp,
    const math::TrackedValue<T>& raan,
    const math::TrackedValue<T>& xnoi)
{
    using TV = math::TrackedValue<T>;
    using math::exact;
    using math::ratio;

    // DSCOM magic numbers as honest decimal-truncation-bounded constants (CR1):
    // each .value is bit-identical to the original double literal, so the value
    // path is preserved while the accuracy budget is no longer over-claimed.
    const TV ZNS    = TV::model_coefficient("1.19459e-5");
    const TV C1SS   = TV::model_coefficient("2.9864797e-6");
    const TV ZES    = TV::model_coefficient("0.01675");
    const TV ZNL    = TV::model_coefficient("1.5835218e-4");
    const TV C1L    = TV::model_coefficient("4.7968065e-7");
    const TV ZEL    = TV::model_coefficient("0.05490");
    const TV ZCOSIS = TV::model_coefficient("0.91744867");
    const TV ZSINI  = TV::model_coefficient("0.39785416");
    const TV ZSINGS = TV::model_coefficient("-0.98088458");
    const TV ZCOSGS = TV::model_coefficient("0.1945905");

    const TV sinq = sin(raan);
    const TV cosq = cos(raan);
    const TV sing = sin(argp);
    const TV cosg = cos(argp);

    // Lunar ascending node setup (dnwrnr lines 726-743).
    const TV xnodce = math::wrap_two_pi(
        TV::model_coefficient("4.5236020")
        - TV::model_coefficient("9.2422029e-4") * epoch_jd_j1900);
    const TV stem = sin(xnodce);
    const TV ctem = cos(xnodce);
    const TV zcosil = TV::model_coefficient("0.91375164")
                    - TV::model_coefficient("0.03568096") * ctem;
    const TV zsinil = sqrt(exact<T>(1) - zcosil * zcosil);
    const TV zsinhl = TV::model_coefficient("0.089683511") * stem / zsinil;
    const TV zcoshl = sqrt(exact<T>(1) - zsinhl * zsinhl);
    const TV cz  = TV::model_coefficient("4.7199672")
                 + TV::model_coefficient("0.22997150") * epoch_jd_j1900;
    const TV gam = TV::model_coefficient("5.8351514")
                 + TV::model_coefficient("0.0019443680") * epoch_jd_j1900;
    ds.zmol = math::wrap_two_pi(cz - gam);
    TV zx = TV::model_coefficient("0.39785416") * stem / zsinil;
    TV zy = zcoshl * ctem + TV::model_coefficient("0.91744867") * zsinhl * stem;
    zx = atan2(zx, zy);
    zx = gam + zx - xnodce;
    const TV zcosgl = cos(zx);
    const TV zsingl = sin(zx);
    ds.zmos = math::wrap_two_pi(
        TV::model_coefficient("6.2565837")
        + TV::model_coefficient("0.017201977") * epoch_jd_j1900);

    TV zcosg = ZCOSGS;
    TV zsing = ZSINGS;
    TV zcosi = ZCOSIS;
    TV zsini = ZSINI;
    TV zcosh = cosq;
    TV zsinh = sinq;
    TV cc = C1SS;
    TV zn = ZNS;
    TV ze = ZES;

    TV sse = exact<T>(0), ssi = exact<T>(0), ssl = exact<T>(0), ssg = exact<T>(0), ssh_raw = exact<T>(0);
    TV se = exact<T>(0), si = exact<T>(0), sl = exact<T>(0), sgh = exact<T>(0), shdq = exact<T>(0);

    for (int cnt = 0; cnt < 2; cnt++) {
        const TV a1  = zcosg * zcosh + zsing * zcosi * zsinh;
        const TV a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
        const TV a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
        const TV a8  = zsing * zsini;
        const TV a9  = zsing * zsinh + zcosg * zcosi * zcosh;
        const TV a10 = zcosg * zsini;
        const TV a2 = cosio * a7 + sinio * a8;
        const TV a4 = cosio * a9 + sinio * a10;
        const TV a5 = -sinio * a7 + cosio * a8;
        const TV a6 = -sinio * a9 + cosio * a10;
        const TV x1 = a1 * cosg + a2 * sing;
        const TV x2 = a3 * cosg + a4 * sing;
        const TV x3 = -a1 * sing + a2 * cosg;
        const TV x4 = -a3 * sing + a4 * cosg;
        const TV x5 = a5 * sing;
        const TV x6 = a6 * sing;
        const TV x7 = a5 * cosg;
        const TV x8 = a6 * cosg;
        const TV z31 = exact<T>(12) * x1 * x1 - exact<T>(3) * x3 * x3;
        const TV z32 = exact<T>(24) * x1 * x2 - exact<T>(6) * x3 * x4;
        const TV z33 = exact<T>(12) * x2 * x2 - exact<T>(3) * x4 * x4;
        TV z1 = exact<T>(3) * (a1 * a1 + a2 * a2) + z31 * eosq;
        TV z2 = exact<T>(6) * (a1 * a3 + a2 * a4) + z32 * eosq;
        TV z3 = exact<T>(3) * (a3 * a3 + a4 * a4) + z33 * eosq;
        const TV z11 = -exact<T>(6) * a1 * a5
                     + eosq * (-exact<T>(24) * x1 * x7 - exact<T>(6) * x3 * x5);
        const TV z12 = -exact<T>(6) * (a1 * a6 + a3 * a5)
                     + eosq * (-exact<T>(24) * (x2 * x7 + x1 * x8) - exact<T>(6) * (x3 * x6 + x4 * x5));
        const TV z13 = -exact<T>(6) * a3 * a6
                     + eosq * (-exact<T>(24) * x2 * x8 - exact<T>(6) * x4 * x6);
        const TV z21 = exact<T>(6) * a2 * a5
                     + eosq * (exact<T>(24) * x1 * x5 - exact<T>(6) * x3 * x7);
        const TV z22 = exact<T>(6) * (a4 * a5 + a2 * a6)
                     + eosq * (exact<T>(24) * (x2 * x5 + x1 * x6) - exact<T>(6) * (x4 * x7 + x3 * x8));
        const TV z23 = exact<T>(6) * a4 * a6
                     + eosq * (exact<T>(24) * x2 * x6 - exact<T>(6) * x4 * x8);
        z1 = z1 + z1 + betao2 * z31;
        z2 = z2 + z2 + betao2 * z32;
        z3 = z3 + z3 + betao2 * z33;

        const TV s3 = cc * xnoi;
        const TV s2 = -ratio<T>(1, 2) * s3 / betao;
        const TV s4 = s3 * betao;
        const TV s1 = -exact<T>(15) * eccen * s4;
        const TV s5 = x1 * x3 + x2 * x4;
        const TV s6 = x2 * x3 + x1 * x4;
        const TV s7 = x2 * x4 - x1 * x3;

        se  = s1 * zn * s5;
        si  = s2 * zn * (z11 + z13);
        sl  = -zn * s3 * (z1 + z3 - exact<T>(14) - exact<T>(6) * eosq);
        sgh = s4 * zn * (z31 + z33 - exact<T>(6));

        // Lyddane shdq with low-inclination guard (dnwrnr lines 825-833).
        // Threshold sin(3°) ≈ 0.05234 maps the 3° inclination guard. The branch
        // selects on the value (a model decision), so it matches the reference.
        if (abs(sinio).value < T(5.2335956e-2)) {
            shdq = exact<T>(0);
        } else {
            shdq = (-zn * s2 * (z21 + z23)) / sinio;
        }

        // Per-body coefficient assignment (dnwrnr lines 835-846): solar pass
        // writes into ee2/e3/xi*/xl*/xgh*/xh* slots; we then copy to s* below.
        ds.ee2  =  exact<T>(2) * s1 * s6;
        ds.e3   =  exact<T>(2) * s1 * s7;
        ds.xi2  =  exact<T>(2) * s2 * z12;
        ds.xi3  =  exact<T>(2) * s2 * (z13 - z11);
        ds.xl2  = -exact<T>(2) * s3 * z2;
        ds.xl3  = -exact<T>(2) * s3 * (z3 - z1);
        ds.xl4  = -exact<T>(2) * s3 * (-exact<T>(21) - exact<T>(9) * eosq) * ze;
        ds.xgh2 =  exact<T>(2) * s4 * z32;
        ds.xgh3 =  exact<T>(2) * s4 * (z33 - z31);
        ds.xgh4 = -exact<T>(18) * s4 * ze;
        ds.xh2  = -exact<T>(2) * s2 * z22;
        ds.xh3  = -exact<T>(2) * s2 * (z23 - z21);

        if (cnt == 1) break;

        // End of solar pass: stash solar→sse/i/l/g/h_raw, copy →s*2/3/4/h2/3.
        sse = se;
        ssi = si;
        ssl = sl;
        ssh_raw = shdq;
        ssg = sgh - cosio * shdq;

        ds.se2  = ds.ee2;
        ds.si2  = ds.xi2;
        ds.sl2  = ds.xl2;
        ds.sgh2 = ds.xgh2;
        ds.sh2  = ds.xh2;
        ds.se3  = ds.e3;
        ds.si3  = ds.xi3;
        ds.sl3  = ds.xl3;
        ds.sgh3 = ds.xgh3;
        ds.sh3  = ds.xh3;
        ds.sl4  = ds.xl4;
        ds.sgh4 = ds.xgh4;

        // Switch to lunar trig for second pass.
        zcosg = zcosgl;
        zsing = zsingl;
        zcosi = zcosil;
        zsini = zsinil;
        zcosh = zcoshl * cosq + zsinhl * sinq;
        zsinh = sinq * zcoshl - cosq * zsinhl;
        zn = ZNL;
        cc = C1L;
        ze = ZEL;
    }

    // Final aggregation (dnwrnr lines 883-887): add lunar contribution. Parens
    // match the reference's `ssg += sgh - cosio*shdq` associativity exactly
    // (ssg + (sgh - cosio*shdq)), NOT (ssg + sgh) - ..., to stay bit-exact.
    sse = sse + se;
    ssi = ssi + si;
    ssl = ssl + sl;
    ssg = ssg + (sgh - cosio * shdq);
    ssh_raw = ssh_raw + shdq;

    ds.sse = sse;
    ds.ssi = ssi;
    ds.ssl = ssl;
    ds.ssg = ssg;
    ds.ssh = ssh_raw;
}

/**
 * @brief DSPER sun/moon long-period periodic corrections (dnwrnr 1088-1208).
 *
 * Inputs: tsince [min], ds.s2/3/4 and ds.x2/3/4 coefficients, ds.zmos/zmol.
 * In-place updates: em, xinc, omgasm, xnodes, xll.
 */
template<typename T>
void apply_deep_periodics(
    const DeepSpaceInit<T>& ds,
    const math::TrackedValue<T>& tsince,
    math::TrackedValue<T>& em, math::TrackedValue<T>& xinc,
    math::TrackedValue<T>& omgasm, math::TrackedValue<T>& xnodes,
    math::TrackedValue<T>& xll)
{
    using TV = math::TrackedValue<T>;
    using math::exact;
    using math::ratio;

    const TV ZES = TV::model_coefficient("0.01675");
    const TV ZNS = TV::model_coefficient("1.19459e-5");
    const TV ZNL = TV::model_coefficient("1.5835218e-4");
    const TV ZEL = TV::model_coefficient("0.05490");
    const T kPI = boost::math::constants::pi<T>();  // branch comparison only

    // Solar terms. (DS1: tracked, so the periodic correction carries a budget.)
    TV zm = ds.zmos + ZNS * tsince;
    TV zf = zm + exact<T>(2) * ZES * sin(zm);
    TV sinzf = sin(zf);
    TV f2 = ratio<T>(1, 2) * sinzf * sinzf - ratio<T>(1, 4);
    TV f3 = -ratio<T>(1, 2) * sinzf * cos(zf);

    const TV ses  = ds.se2  * f2 + ds.se3  * f3;
    const TV sis  = ds.si2  * f2 + ds.si3  * f3;
    const TV sls  = ds.sl2  * f2 + ds.sl3  * f3 + ds.sl4  * sinzf;
    const TV sghs = ds.sgh2 * f2 + ds.sgh3 * f3 + ds.sgh4 * sinzf;
    const TV shs  = ds.sh2  * f2 + ds.sh3  * f3;

    // Lunar terms.
    zm = ds.zmol + ZNL * tsince;
    zf = zm + exact<T>(2) * ZEL * sin(zm);
    sinzf = sin(zf);
    f2 = ratio<T>(1, 2) * sinzf * sinzf - ratio<T>(1, 4);
    f3 = -ratio<T>(1, 2) * sinzf * cos(zf);

    const TV sel  = ds.ee2  * f2 + ds.e3   * f3;
    const TV sil  = ds.xi2  * f2 + ds.xi3  * f3;
    const TV sll  = ds.xl2  * f2 + ds.xl3  * f3 + ds.xl4  * sinzf;
    const TV sghl = ds.xgh2 * f2 + ds.xgh3 * f3 + ds.xgh4 * sinzf;
    const TV shl  = ds.xh2  * f2 + ds.xh3  * f3;

    const TV pe = ses + sel;
    const TV pinc = sis + sil;
    const TV pl = sls + sll;
    const TV pgh = sghs + sghl;
    const TV ph = shs + shl;

    xinc = xinc + pinc;
    em = em + pe;

    const TV sinis = sin(xinc);
    const TV cosis = cos(xinc);

    if (xinc.value >= T(0.2)) {
        omgasm = omgasm + (pgh - cosis * ph / sinis);
        xnodes = xnodes + ph / sinis;
        xll = xll + pl;
    } else {
        // Lyddane modification at low inclination.
        const TV sinok = sin(xnodes);
        const TV cosok = cos(xnodes);
        TV alfdp = sinis * sinok;
        TV betdp = sinis * cosok;
        const TV dalf = ph * cosok + pinc * cosis * sinok;
        const TV dbet = -ph * sinok + pinc * cosis * cosok;
        alfdp = alfdp + dalf;
        betdp = betdp + dbet;
        xnodes = math::wrap_two_pi(xnodes);
        TV xls = xll + omgasm + cosis * xnodes;
        TV dls = pl + pgh - pinc * xnodes * sinis;
        xls = xls + dls;
        const TV oldxnodes = xnodes;
        xnodes = atan2(alfdp, betdp);
        if (std::fabs(double(oldxnodes.value) - double(xnodes.value)) > double(kPI)) {
            if (xnodes.value < oldxnodes.value) xnodes = xnodes + math::two_pi<T>();
            else xnodes = xnodes - math::two_pi<T>();
        }
        xll = xll + pl;
        omgasm = xls - xll - cosis * xnodes;
    }
}

} // namespace detail

/**
 * @brief Initialize deep-space propagation (SR3 §6).
 */
template<typename T>
DeepSpaceInit<T> initialize_deep_space(
    const ModelConfiguration<T>& config,
    const tle::TleElements<T>& elements,
    const T& tolerance)
{
    using math::exact;
    using math::ratio;

    DeepSpaceInit<T> ds;

    // Earth constants.
    ds.re_km = config.ellipsoid.a;
    math::TrackedValue<T> mu_km = config.ellipsoid.GM;
    ds.xke = exact<T>(60) / sqrt(ds.re_km * ds.re_km * ds.re_km / mu_km);
    math::TrackedValue<T> J2 = config.Jn(2);
    ds.CK2 = J2 / exact<T>(2);
    ds.J3 = config.Jn(3);
    ds.A30 = config.A(3, 0);
    math::TrackedValue<T> J4 = config.Jn(4);

    // TLE.
    ds.e0 = elements.eccentricity;
    ds.i0 = elements.inclination;
    ds.omega0 = elements.arg_perigee;
    ds.Omega0 = elements.raan;
    ds.M0 = elements.mean_anomaly;
    ds.bstar = elements.bstar;

    // Trig.
    ds.cosio = cos(ds.i0);
    ds.sinio = sin(ds.i0);
    ds.theta2 = ds.cosio * ds.cosio;
    ds.x3thm1 = exact<T>(3) * ds.theta2 - exact<T>(1);
    ds.x1mth2 = exact<T>(1) - ds.theta2;
    ds.x7thm1 = exact<T>(7) * ds.theta2 - exact<T>(1);

    // Element recovery (Brouwer Kozai inversion).
    orbit::RecoveredElements<T> recovered = orbit::recover_mean_elements(
        elements.mean_motion, ds.e0, ds.cosio, ds.CK2, ds.xke,
        ds.re_km, tolerance);
    ds.a0 = recovered.a0;
    ds.n0 = recovered.n0;
    ds.beta0 = recovered.beta0;
    ds.perigee_km = recovered.perigee_km;
    ds.period_min = recovered.period_min;

    // Brouwer secular rates.
    math::TrackedValue<T> eosq_tv = ds.e0 * ds.e0;
    perturbation::BrouwerSecularRates<T> rates = config.model_functions.secular_rates(
        ds.n0, ds.a0, eosq_tv, ds.cosio, J2, J4);
    ds.M_dot = rates.M_dot;
    ds.omega_dot = rates.omega_dot;
    ds.Omega_dot = rates.Omega_dot;

    // Earth rotation rate (kTHDT): ellipsoid omega in rad/s → rad/min.
    ds.theta_dot = config.ellipsoid.omega * exact<T>(60);

    // GMST at epoch.
    math::TrackedValue<T> gmst_tv = config.model_functions.sidereal_time(elements.epoch_jd);
    ds.gsto = gmst_tv;

    // DSCOM coefficient buildup. Inputs are TrackedValue<T> so the sun/moon
    // secular + periodic coefficients carry a budget that propagates from the
    // epoch elements and scales with T (DS1). Each input's .value is
    // bit-identical to the prior double form (same op order).
    // J1900 day count via the single-sourced view (L1): bit-identical to the
    // prior inline `elements.epoch_jd - exact<T>(2415020)`; the magic base now
    // lives once in astronomy/epoch.h.
    math::TrackedValue<T> epoch_j1900_tv = astronomy::days_since_1900<T>(elements.epoch_jd);
    math::TrackedValue<T> betao2_tv = exact<T>(1) - eosq_tv;
    math::TrackedValue<T> betao_tv = sqrt(betao2_tv);
    math::TrackedValue<T> xnoi_tv = exact<T>(1) / ds.n0;

    detail::build_dpper_coefficients(
        ds, epoch_j1900_tv,
        eosq_tv, ds.sinio, ds.cosio,
        betao_tv, betao2_tv,
        ds.e0, ds.omega0, ds.Omega0,
        xnoi_tv);

    // Resonance initialization (SR3 §6). The ss* linear-rate terms and
    // theta_dot are now TrackedValue<T>, so they feed the resonance with a real
    // budget directly — no epoch re-wrap.
    ds.resonance = perturbation::initialize_resonance(
        ds.n0, ds.e0, ds.i0, ds.omega0, ds.Omega0, ds.M0,
        gmst_tv,
        ds.M_dot, ds.omega_dot, ds.Omega_dot,
        ds.ssl, ds.ssg, ds.ssh,
        ds.theta_dot);

    // A3/CK2 ratio for the J₃ long-period coefficients. The actual xlcof/aycof
    // are recomputed at the PERTURBED inclination inside propagate_deep_space
    // (step 8), so only this ratio is stored here. (Phase 6 survey: the former
    // epoch-inclination ds.xlcof/ds.aycof were write-only dead fields — never
    // read by propagate — and were removed.)
    math::TrackedValue<T> A30_over_CK2 = ds.A30 / ds.CK2;
    ds.a3ovk2 = A30_over_CK2;

    return ds;
}

/**
 * @brief Propagate a deep-space satellite to time tsince from epoch (SR3 §6).
 *
 * Pipeline (mirrors dnwrnr FindPositionSDP4 lines 231-358):
 *   1. Brouwer secular advance.
 *   2. DSPACE secular (sun+moon linear-rate corrections + resonance leapfrog).
 *   3. Recompute a from xn.
 *   4. DSPER long-period periodics (DPPER).
 *   5. Inclination sign correction.
 *   6. Recompute inclination-dependent coefficients.
 *   7. Long-period periodics (J₃), Kepler, osculating, short-period, state.
 */
template<typename T>
StateVector<T> propagate_deep_space(
    const DeepSpaceInit<T>& ds,
    const NearSpaceInit<T>& ns,
    const ModelConfiguration<T>& config,
    const math::TrackedValue<T>& tsince,
    const T& tolerance)
{
    using math::exact;
    using math::ratio;
    using TV = math::TrackedValue<T>;

    // Error sentinel: the reference SGP4 returns r=(NaN,NaN,NaN) with a non-zero
    // error code when a propagation guard trips (mean motion <= 0, perturbed
    // eccentricity out of range, sub-orbital). A NaN state is the unambiguous
    // analog (a silent (0,0,0) could be mistaken for a valid position). See
    // python_sgp4_rhodes propagation.py:1769/1837/1885 and tests.py:727.
    StateVector<T> (*nan_state)() = []() {
        StateVector<T> sv{};
        const T nan = std::numeric_limits<T>::quiet_NaN();
        sv.position_km.x.value = nan;
        sv.position_km.y.value = nan;
        sv.position_km.z.value = nan;
        sv.velocity_km_s.x.value = nan;
        sv.velocity_km_s.y.value = nan;
        sv.velocity_km_s.z.value = nan;
        return sv;
    };

    // The deep-space secular/periodic evolution (steps 1-9) runs in TrackedValue<T>
    // (DS1), so the three-error budget propagates through it and tightens with a
    // wider T — instead of being re-injected from the epoch at the Kepler stage.
    // Each .value matches the prior raw-double path (same op order), so OR1 and
    // the 33/33 reference are bit-exact.

    // 1) Brouwer secular advance.
    const TV t2 = tsince * tsince;
    TV xmdf   = ds.M0     + ds.M_dot     * tsince;
    TV argpdf = ds.omega0 + ds.omega_dot * tsince;
    TV nodedf = ds.Omega0 + ds.Omega_dot * tsince;

    // 2) Near-earth secular atmospheric drag — SIMPLE model only (reference
    //    propagation.py:1721-1723). Every deep-space satellite uses the simple
    //    drag model: sgp4init forces isimp=1 for period >= 225 min
    //    (propagation.py:1496-1499), independent of perigee. So only the base
    //    tempa/tempe/templ and the nodecf node advance apply here, NOT the
    //    non-simple delomg/delm/D2-4/C5/t3-5cof terms (the `if isimp != 1` block,
    //    propagation.py:1725-1743). Do NOT gate on ns.use_simple_model — that is a
    //    perigee<220 km test, so deep-space sats with perigee>220 km would wrongly
    //    receive the non-simple C5 sin(M) periodic drag. DSPACE lunar-solar adds on
    //    top (step 3); the a/e/M scalings apply after DSPACE (step 4). Coefficients
    //    come from the always-computed near-space init (ns), whose recovered
    //    elements (a0,n0,e0,i0) equal ds's.
    TV tempa = exact<T>(1) - ns.C1 * tsince;
    TV tempe = ds.bstar * ns.C4 * tsince;
    TV templ = ns.t2cof * t2;
    TV mm    = xmdf;
    TV argpm = argpdf;
    TV nodem = nodedf + ns.Omega_dot_nkc * t2;

    // 3) DSPACE secular: add sun+moon linear contributions and run resonance.
    TV xn   = ds.n0;
    TV em   = ds.e0;
    TV xinc = ds.i0;

    TV xll    = mm    + ds.ssl * tsince;
    TV omgasm = argpm + ds.ssg * tsince;
    TV xnodes = nodem + ds.ssh * tsince;
    em   = em   + ds.sse * tsince;
    xinc = xinc + ds.ssi * tsince;

    // Resonance leapfrog (no-op if rs.type == NONE). Fully TrackedValue<T>
    // (DS1.4): the integrated mean longitude xll and mean motion xn carry the
    // propagated budget. xnodes and omgasm are read-only resonance targets.
    perturbation::step_resonance(
        ds.resonance, ds.omega0, ds.omega_dot,
        tsince, ds.theta_dot, ds.gsto,
        xnodes, omgasm,
        xll, xn);

    // 4) Apply drag scalings to a, e, M (reference propagation.py:1769-1794):
    //    am = (xke/nm)^(2/3)·tempa²,  em -= tempe,  M += n0·templ.
    if (xn.value <= T(0)) {        // reference error 2: mean motion <= 0
        return nan_state();
    }
    TV a = pow(ds.xke / xn, T(2) / T(3)) * tempa * tempa;
    em = em - tempe;
    xll = xll + ds.n0 * templ;

    // 5) DSPER long-period periodics (sun+moon, time-varying).
    detail::apply_deep_periodics(ds, tsince, em, xinc, omgasm, xnodes, xll);

    // 6) Inclination sign correction.
    if (xinc.value < T(0)) {
        xinc = -xinc;
        TV pi_tv = math::pi<T>();
        xnodes = xnodes + pi_tv;
        omgasm = omgasm - pi_tv;
    }

    TV xl = xll + omgasm + xnodes;
    TV omega = omgasm;
    TV e = em;

    // 7) Perturbed-eccentricity range check (reference error 3: ep<0 or ep>1,
    //    propagation.py:1837, checked AFTER the deep-space DPPER periodics).
    if (e.value < T(0) || e.value > T(1)) {
        return nan_state();
    }
    if (e.value < T(1.0e-6))         e = TV(T(1.0e-6),         e.errors);
    if (e.value > T(1) - T(1.0e-6))  e = TV(T(1) - T(1.0e-6),  e.errors);

    // 8) Recompute perturbed long-period coefficients at current inclination.
    const TV cosi_p = cos(xinc);
    const TV sini_p = sin(xinc);
    const TV theta2_p = cosi_p * cosi_p;
    const TV x3thm1_p = exact<T>(3) * theta2_p - exact<T>(1);
    const TV x1mth2_p = exact<T>(1) - theta2_p;
    const TV x7thm1_p = exact<T>(7) * theta2_p - exact<T>(1);
    TV xlcof_p, aycof_p;
    if (std::fabs(double(cosi_p.value) + 1.0) > 1.5e-12) {
        xlcof_p = ratio<T>(1, 8) * ds.a3ovk2 * sini_p
                * (exact<T>(3) + exact<T>(5) * cosi_p) / (exact<T>(1) + cosi_p);
    } else {
        xlcof_p = ratio<T>(1, 8) * ds.a3ovk2 * sini_p
                * (exact<T>(3) + exact<T>(5) * cosi_p) / TV::from_truncated_decimal("1.5e-12");
    }
    aycof_p = ratio<T>(1, 4) * ds.a3ovk2 * sini_p;

    // 9) Long-period periodics (J₃) at current inclination.
    const TV beta2 = exact<T>(1) - e * e;
    const TV xn_actual = ds.xke / pow(a, T(3) / T(2));
    const TV axn = e * cos(omega);
    const TV temp_lp = exact<T>(1) / (a * beta2);
    const TV xll_jp = temp_lp * xlcof_p * axn;
    const TV aynl = temp_lp * aycof_p;
    const TV xlt = xl + xll_jp;
    const TV ayn = e * sin(omega) + aynl;
    const TV elsq = axn * axn + ayn * ayn;
    if (elsq.value >= T(1)) {
        return nan_state();
    }

    // 10) Modified Kepler solver via injected lambda. No epoch re-injection: the
    //     budget on axn/ayn/capu now flows from the propagated evolution.
    TV capu_tv = math::wrap_two_pi(xlt - xnodes);
    TV E_plus_w_tv = config.model_functions.kepler_solver(
        axn, ayn, capu_tv, tolerance);

    // 11) Osculating elements.
    orbit::OsculatingState<T> osc = orbit::compute_osculating(
        E_plus_w_tv, axn, ayn, a, ds.xke);

    // 12) Short-period (J₂) corrections — fed by the propagated TrackedValues.
    const TV& cos_i_tv  = cosi_p;
    const TV& sin_i_tv  = sini_p;
    const TV& cos2_i_tv = theta2_p;
    const TV& sin2_i_tv = x1mth2_p;
    const TV& x3thm1_tv = x3thm1_p;
    const TV& x7thm1_tv = x7thm1_p;
    const TV& inc_tv   = xinc;
    const TV& Omega_tv = xnodes;
    const TV& n_tv     = xn_actual;

    perturbation::CorrectedElements<T> corr = perturbation::apply_short_period(
        osc.r, osc.u, inc_tv, Omega_tv,
        osc.rdot, osc.rfdot, n_tv,
        osc.pl, osc.beta_l,
        ds.CK2, cos2_i_tv, sin2_i_tv,
        cos_i_tv, sin_i_tv,
        x3thm1_tv, x7thm1_tv,
        osc.sin_2u, osc.cos_2u);

    // 13) State.
    return orbit::elements_to_state(
        corr.r, corr.u, corr.i, corr.Omega,
        corr.rdot, corr.rfdot, ds.re_km);
}

} // namespace sgp4
