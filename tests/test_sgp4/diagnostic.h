#pragma once

/// @file diagnostic.h
/// @brief Diagnostic functions for comparing intermediate SGP4 values
///        against the Vallado/dnwrnr reference implementation.

#include <iostream>
#include <iomanip>
#include <cmath>
#include "math/tracked_value.h"
#include "geodesy/equipotential_ellipsoid.h"
#include "sgp4/near_space.h"
#include "sgp4/model_selector.h"
#include "tle/tle_parser.h"

/// Compare a computed value against a reference and report the difference.
inline void compare(const char* name, double computed, double reference) {
    double diff = std::abs(computed - reference);
    double rel = (reference != 0.0) ? diff / std::abs(reference) : diff;
    bool ok = rel < 1e-6;
    std::cout << (ok ? "  OK  " : "  BAD ")
              << std::left << std::setw(20) << name
              << std::right << std::setprecision(12) << std::setw(18) << computed
              << std::setw(18) << reference
              << std::scientific << std::setw(14) << diff
              << std::setw(14) << rel
              << std::fixed << "\n";
}

/// Run diagnostic comparison of element recovery for satellite 00005.
///
/// ALL reference constants are derived from the ModelConfiguration —
/// the same single source of truth the propagator uses. No hardcoded
/// "dnwrnr" or "Vallado" constants anywhere.
template<typename T>
void diagnose_sat00005(const tle::TleData& td) {
    using math::TrackedValue;

    std::cout << "\n=== DIAGNOSTIC: Sat 00005 Element Recovery ===\n";
    std::cout << std::setw(24) << "Name" << std::setw(18) << "Computed"
              << std::setw(18) << "Reference" << std::setw(14) << "AbsDiff"
              << std::setw(14) << "RelDiff" << "\n";
    std::cout << std::string(82, '-') << "\n";

    // --- All constants from the model configuration ---
    sgp4::ModelConfiguration<T> config = sgp4::ModelSelector<T>::select("sgp4_standard", T(1e-12));
    double kXKMPER = double(config.ellipsoid.a.value);
    double kGM     = double(config.ellipsoid.GM.value);
    double kXKE    = 60.0 / std::sqrt(kXKMPER*kXKMPER*kXKMPER / kGM);
    double kJ2     = double(config.ellipsoid.J2n(1).value);
    double kCK2    = kJ2 / 2.0;
    double kTWOTHIRD = 2.0 / 3.0;
    double our_J3  = double(config.zonals.J3().value);
    double our_J4  = double(config.zonals.Jn(4).value);

    // Input elements
    double n0_rev_day = td.mean_motion_rev_day;
    double n0 = n0_rev_day * 2.0 * 3.14159265358979 / 1440.0;
    double e0 = td.eccentricity;
    double i0 = td.inclination_deg * 3.14159265358979 / 180.0;

    compare("n0 [rad/min]", n0, n0);

    // Reference element recovery (dnwrnr style)
    double a1_ref = std::pow(kXKE / n0, kTWOTHIRD);
    double cosio = std::cos(i0);
    double theta2 = cosio * cosio;
    double x3thm1 = 3.0 * theta2 - 1.0;
    double eosq = e0 * e0;
    double betao2 = 1.0 - eosq;
    double betao = std::sqrt(betao2);
    double temp = (1.5 * kCK2) * x3thm1 / (betao * betao2);
    double del1_ref = temp / (a1_ref * a1_ref);
    double a0_ref = a1_ref * (1.0 - del1_ref * (1.0/3.0 + del1_ref * (1.0 + del1_ref * 134.0/81.0)));
    double del0_ref = temp / (a0_ref * a0_ref);
    double n0pp_ref = n0 / (1.0 + del0_ref);
    double a0pp_ref = a0_ref / (1.0 - del0_ref);
    double perigee_ref = (a0pp_ref * (1.0 - e0) - 1.0) * kXKMPER;

    // Our element recovery (config already created above)
    tle::TleElements<T> tle_elems = tle::TleElements<T>::from_tle_data(td);
    sgp4::NearSpaceInit<T> ns = sgp4::initialize_near_space(config, tle_elems, T(1e-12));

    compare("a1", double(ns.a0.value), a0pp_ref);  // our a0 is a0''
    compare("cosio", double(ns.cosio.value), cosio);
    compare("theta2", double(ns.theta2.value), theta2);
    compare("x3thm1", double(ns.x3thm1.value), x3thm1);
    compare("betao", double(ns.beta0.value), betao);
    compare("n0'' [rad/min]", double(ns.n0.value), n0pp_ref);
    compare("a0'' [ER]", double(ns.a0.value), a0pp_ref);
    compare("perigee [km]", double(ns.perigee_km.value), perigee_ref);
    compare("e0", double(ns.e0.value), e0);

    // Secular rates comparison — first-order J₂ only (for quick sanity check).
    // Our M_dot includes J₂² and J₄ terms, so it will differ by ~1e-5 rad/min.
    // The full comparison against second-order rates is in the drag coefficients section below.
    double pinvsq = 1.0 / (a0pp_ref * a0pp_ref * betao2 * betao2);
    double temp1_ref = 1.5 * kCK2 * pinvsq * n0pp_ref;
    double xmdot_j2_only = n0pp_ref + 0.5 * temp1_ref * betao * x3thm1;

    std::cout << "  M_dot (J2-only ref): " << std::setprecision(12) << xmdot_j2_only
              << "  (ours: " << double(ns.M_dot.value) << ", diff from J2^2+J4 terms)\n";

    // C1 coefficient
    compare("C1", double(ns.C1.value), 0.0);  // TODO: compute reference C1

    // xi = 1/(a0'' - s)
    double s_ref = 1.0 + 78.0/kXKMPER;
    double xi_ref = 1.0 / (a0pp_ref - s_ref);
    compare("xi", double(ns.xi.value), xi_ref);

    // eta = a0'' * e0 * xi
    double eta_ref = a0pp_ref * e0 * xi_ref;
    compare("eta", double(ns.eta.value), eta_ref);

    // --- At t=0 propagation check ---
    std::cout << "\n--- t=0 position/radius check ---\n";

    double M0 = td.mean_anomaly_deg * 3.14159265358979 / 180.0;
    compare("M0 [rad]", double(ns.M0.value), M0);

    // Reference |r| at t=0 from tcppver.out
    double ref_r = std::sqrt(7022.46529266*7022.46529266 + 1400.08296755*1400.08296755
                             + 0.03995155*0.03995155);
    std::cout << "  Reference |r| at t=0 = " << ref_r << " km\n";

    // Velocity conversion: Earth_radii/min → km/s = XKMPER/60
    // NOTE: xke is NOT part of this conversion (fixed in near_space.h)
    double vkmpersec = kXKMPER / 60.0;
    compare("vkmpersec", vkmpersec, kXKMPER / 60.0);

    // --- Step-by-step propagation at t=0 using reference formulas ---
    std::cout << "\n--- Step-by-step propagation at t=0 ---\n";

    double omega0_val = td.arg_perigee_deg * 3.14159265358979 / 180.0;
    double Omega0_val = td.raan_deg * 3.14159265358979 / 180.0;

    // At t=0, all secular terms are zero. The key values are:
    // xmdf = M0, omgadf = omega0, xnoddf = Omega0
    // tempa = 1, tempe = 0, templ = 0
    // So: a = a0'', e = e0, xl = M0 + omega0 + Omega0

    double a_t0 = a0pp_ref;  // tempa=1
    double e_t0 = e0;        // tempe=0
    double xl_t0 = M0 + omega0_val + Omega0_val;

    // axN = e*cos(omega)
    double axN = e_t0 * std::cos(omega0_val);
    // ayn = e*sin(omega) + aynl (long-period correction)
    // For t=0, the long-period correction is small
    // A30 = -J3, so A30/CK2 = -J3/CK2 (J3 from config, defined at top)
    double A30_CK2 = (-our_J3) / kCK2;
    double aycof_val = 0.25 * A30_CK2 * std::sin(i0);
    double beta2_t0 = 1.0 - e_t0 * e_t0;
    double n_t0 = kXKE / std::pow(a_t0, 1.5);
    // xlcof: factor is 1/8, not 1/4 (fixed in near_space.h, must match here)
    double xlcof_val = (1.0/8.0) * A30_CK2 * std::sin(i0) * (3.0 + 5.0*cosio) / (1.0 + cosio);
    double aynl = (1.0 / (a_t0 * beta2_t0)) * aycof_val;
    double xll = (1.0 / (a_t0 * beta2_t0)) * xlcof_val * axN;
    double xlt = xl_t0 + xll;
    double ayn = e_t0 * std::sin(omega0_val) + aynl;

    compare("axN", axN, axN);
    compare("ayn", ayn, ayn);

    // capu = (xlt - Omega0) mod 2pi
    double capu = std::fmod(xlt - Omega0_val, 2*3.14159265358979);
    if (capu < 0) capu += 2*3.14159265358979;
    compare("capu (mod 2pi)", capu, capu);

    // SGP4 modified Kepler equation: x - ayn*cos(x) + axN*sin(x) = U
    // where x = (E + omega), solved with Newton iteration starting from x = U.
    double elsq = axN*axN + ayn*ayn;

    double epw = capu;  // starting value: (E+ω)₁ = U
    for (int i = 0; i < 30; ++i) {
        double sinepw = std::sin(epw);
        double cosepw = std::cos(epw);
        double f = capu - ayn*cosepw + axN*sinepw - epw;
        double fp = -ayn*sinepw - axN*cosepw + 1.0;
        double delta = f / fp;
        epw += delta;
        if (std::abs(delta) < 1e-12) break;
    }
    std::cout << "  (E+w) = " << std::setprecision(12) << epw << " rad\n";

    // Intermediate osculating quantities from (E+ω)
    double sinepw = std::sin(epw);
    double cosepw = std::cos(epw);
    double ecose = axN*cosepw + ayn*sinepw;
    double esine = axN*sinepw - ayn*cosepw;
    double pl = a_t0 * (1.0 - elsq);
    double r_t0 = a_t0 * (1.0 - ecose);
    double rdot_t0 = kXKE * std::sqrt(a_t0) * esine / r_t0;
    double rfdot_t0 = kXKE * std::sqrt(pl) / r_t0;

    std::cout << "  r = " << std::setprecision(12) << r_t0 << " ER ("
              << r_t0 * kXKMPER << " km, pre-short-period correction)\n";
    compare("rdot", rdot_t0, rdot_t0);
    compare("rfdot", rfdot_t0, rfdot_t0);

    // u (argument of latitude) from (E+ω)
    double betal = std::sqrt(1.0 - elsq);
    double temp_kep2 = 1.0 / (1.0 + betal);
    double cosu = (a_t0/r_t0) * (cosepw - axN + ayn*esine*temp_kep2);
    double sinu = (a_t0/r_t0) * (sinepw - ayn - axN*esine*temp_kep2);
    double u_t0 = std::atan2(sinu, cosu);
    compare("u [rad]", u_t0, u_t0);

    // Short-period corrections
    double sin2u = 2*sinu*cosu;
    double cos2u = 2*cosu*cosu - 1;
    double temp_sp2 = kCK2 / pl;
    double temp_sp3 = temp_sp2 / pl;
    double rk = r_t0*(1.0 - 1.5*temp_sp3*betal*x3thm1) + 0.5*temp_sp2*(1-theta2)*cos2u;
    double uk = u_t0 - 0.25*temp_sp3*(7*theta2-1)*sin2u;
    compare("rk [ER]", rk, rk);
    compare("rk [km]", rk * kXKMPER, ref_r);

    // Orientation vectors
    double xnodek = Omega0_val + 1.5*temp_sp3*cosio*sin2u;
    double xinck = i0 + 1.5*temp_sp3*cosio*std::sin(i0)*cos2u;
    double sinuk = std::sin(uk), cosuk = std::cos(uk);
    double sinik = std::sin(xinck), cosik = std::cos(xinck);
    double sinnok = std::sin(xnodek), cosnok = std::cos(xnodek);

    double xmx = -sinnok*cosik;
    double xmy = cosnok*cosik;
    double ux = xmx*sinuk + cosnok*cosuk;
    double uy = xmy*sinuk + sinnok*cosuk;
    double uz = sinik*sinuk;

    double x_pos = rk * ux * kXKMPER;
    double y_pos = rk * uy * kXKMPER;
    double z_pos = rk * uz * kXKMPER;

    compare("x [km]", x_pos, 7022.46529266);
    compare("y [km]", y_pos, -1400.08296755);
    compare("z [km]", z_pos, 0.03995155);

    double rdotk = rdot_t0 - n_t0*temp_sp2*(1-theta2)*sin2u;
    double rfdotk = rfdot_t0 + n_t0*temp_sp2*((1-theta2)*cos2u + 1.5*x3thm1);
    double vx = xmx*cosuk - cosnok*sinuk;
    double vy = xmy*cosuk - sinnok*sinuk;
    double vz = sinik*cosuk;
    double xdot = (rdotk*ux + rfdotk*vx) * vkmpersec;
    double ydot = (rdotk*uy + rfdotk*vy) * vkmpersec;
    double zdot = (rdotk*uz + rfdotk*vz) * vkmpersec;

    compare("vx [km/s]", xdot, 1.893841015);
    compare("vy [km/s]", ydot, 6.405893759);
    compare("vz [km/s]", zdot, 4.534807250);

    std::cout << "\n";

    // =====================================================
    // Drag coefficient comparison
    // =====================================================
    std::cout << "\n=== DIAGNOSTIC: Drag Coefficients ===\n";
    std::cout << std::setw(24) << "Name" << std::setw(18) << "Computed"
              << std::setw(18) << "Reference" << std::setw(14) << "AbsDiff"
              << std::setw(14) << "RelDiff" << "\n";
    std::cout << std::string(82, '-') << "\n";

    // Compute reference drag coefficients using double arithmetic
    // Following the Lane-Hoots derivation step by step
    double s4_ref = 1.0 + 78.0/kXKMPER;
    double q0_ref = 1.0 + 120.0/kXKMPER;
    double qoms24_ref = std::pow(q0_ref - s4_ref, 4.0);

    // Perigee adjustment (if needed for sat 00005)
    if (perigee_ref < 156.0 && perigee_ref >= 98.0) {
        double s_star = a0pp_ref*(1.0-e0) - s4_ref + 1.0;
        double qoms_r4 = std::pow(qoms24_ref, 0.25);
        double new_qoms = qoms_r4 + s4_ref - (1.0 + s_star);
        qoms24_ref = std::pow(new_qoms, 4.0);
        s4_ref = 1.0 + s_star;
    } else if (perigee_ref < 98.0) {
        double s_star = 20.0/kXKMPER;
        double qoms_r4 = std::pow(qoms24_ref, 0.25);
        double new_qoms = qoms_r4 + s4_ref - (1.0 + s_star);
        qoms24_ref = std::pow(new_qoms, 4.0);
        s4_ref = 1.0 + s_star;
    }

    compare("s4", double(ns.s4.value), s4_ref);
    compare("qoms24", double(ns.qoms24.value), qoms24_ref);

    double xi_ref2 = 1.0 / (a0pp_ref - s4_ref);
    double eta_ref2 = a0pp_ref * e0 * xi_ref2;
    double eta2_ref = eta_ref2 * eta_ref2;
    double eeta_ref = e0 * eta_ref2;
    double psisq_ref = std::abs(1.0 - eta2_ref);

    compare("xi", double(ns.xi.value), xi_ref2);
    compare("eta", double(ns.eta.value), eta_ref2);

    double coef_ref = qoms24_ref * std::pow(xi_ref2, 4.0);
    double coef1_ref = coef_ref / std::pow(psisq_ref, 3.5);

    // C2 reference: 0.75 * kCK2 = 0.75 * J2/2 = 0.375 * J2 = Vallado coefficient.
    // Cross-checked vs dnwrnr SGP4.cc:131 (0.75 * kCK2 form).
    double sinio_ref = std::sin(i0);
    double x1mth2_ref = 1.0 - theta2;
    double C2_ref = coef1_ref * n0pp_ref * (a0pp_ref * (1.0 + 1.5*eta2_ref + eeta_ref*(4.0 + eta2_ref))
        + 0.75 * kCK2 * xi_ref2 / psisq_ref * x3thm1
          * (8.0 + 3.0*eta2_ref*(8.0 + eta2_ref)));

    compare("C2", double(ns.C2.value), C2_ref);

    double C1_ref = td.bstar * C2_ref;
    compare("C1", double(ns.C1.value), C1_ref);

    // C3
    // A30 = -J3 (J3 from config, defined at top of function)
    double A30_ref = -our_J3;
    double kA3OVK2 = A30_ref / kCK2;
    double C3_ref2 = (e0 > 1e-4) ?
        coef_ref * xi_ref2 * A30_ref * n0pp_ref * sinio_ref / (kCK2 * e0) : 0.0;
    compare("C3", double(ns.C3.value), C3_ref2);

    // C4
    double C4_ref = 2.0 * n0pp_ref * coef1_ref * a0pp_ref * betao2
        * (eta_ref2*(2.0 + 0.5*eta2_ref) + e0*(0.5 + 2.0*eta2_ref)
           - 2.0 * kCK2 * xi_ref2 / (a0pp_ref * psisq_ref)
             * (-3.0*x3thm1*(1.0 - 2.0*eeta_ref + eta2_ref*(1.5 - 0.5*eeta_ref))
                + 0.75*x1mth2_ref*(2.0*eta2_ref - eeta_ref*(1.0+eta2_ref))
                    *std::cos(2.0*omega0_val)));
    compare("C4", double(ns.C4.value), C4_ref);

    // C5
    double C5_ref = 2.0 * coef1_ref * a0pp_ref * betao2
        * (1.0 + 2.75*(eta2_ref + eeta_ref) + eeta_ref*eta2_ref);
    compare("C5", double(ns.C5.value), C5_ref);

    // D2, D3, D4
    double D2_ref = 4.0 * a0pp_ref * xi_ref2 * C1_ref * C1_ref;
    double temp_d_ref = D2_ref * xi_ref2 * C1_ref / 3.0;
    double D3_ref = (17.0*a0pp_ref + s4_ref) * temp_d_ref;
    double D4_ref = 0.5 * temp_d_ref * a0pp_ref * xi_ref2 * (221.0*a0pp_ref + 31.0*s4_ref) * C1_ref;
    compare("D2", double(ns.D2.value), D2_ref);
    compare("D3", double(ns.D3.value), D3_ref);
    compare("D4", double(ns.D4.value), D4_ref);

    // Time coefficients
    double t2cof_ref = 1.5 * C1_ref;
    double t3cof_ref = D2_ref + 2.0*C1_ref*C1_ref;
    double t4cof_ref = 0.25*(3.0*D3_ref + C1_ref*(12.0*D2_ref + 10.0*C1_ref*C1_ref));
    double t5cof_ref = 0.2*(3.0*D4_ref + 12.0*C1_ref*D3_ref + 6.0*D2_ref*D2_ref
        + 15.0*C1_ref*C1_ref*(2.0*D2_ref + C1_ref*C1_ref));
    compare("t2cof", double(ns.t2cof.value), t2cof_ref);
    compare("t3cof", double(ns.t3cof.value), t3cof_ref);
    compare("t4cof", double(ns.t4cof.value), t4cof_ref);
    compare("t5cof", double(ns.t5cof.value), t5cof_ref);

    // Correction terms
    double omgcof_ref = td.bstar * C3_ref2 * std::cos(omega0_val);
    double xmcof_ref = (e0 > 1e-4) ? -(2.0/3.0) * coef_ref * td.bstar / eeta_ref : 0.0;
    double delmo_ref = std::pow(1.0 + eta_ref2*std::cos(M0), 3.0);
    double sinmo_ref = std::sin(M0);
    compare("omgcof", double(ns.omgcof.value), omgcof_ref);
    compare("xmcof", double(ns.xmcof.value), xmcof_ref);
    compare("delmo", double(ns.delmo.value), delmo_ref);
    compare("sinmo", double(ns.sinmo.value), sinmo_ref);

    // RAAN correction
    double xnodcf_ref = -10.5 * n0pp_ref * kCK2 * cosio / (a0pp_ref*a0pp_ref*betao2) * C1_ref;
    compare("Omega_dot_nkc", double(ns.Omega_dot_nkc.value), xnodcf_ref);

    // Secular rates — all constants from config (defined at top of function)
    std::cout << "\n  Zonal harmonic sources (all from ModelConfiguration):\n";
    std::cout << std::setprecision(15);
    std::cout << "    J2 (ellipsoid):        " << kJ2 << "\n";
    std::cout << "    J3 (zonals table):     " << our_J3 << "\n";
    std::cout << "    J4 (zonals table):     " << our_J4 << "\n";
    std::cout << "    J4 (ellipsoid, unused): " << double(config.ellipsoid.J2n(2).value) << "\n\n" << std::fixed;

    double temp1_ref2 = 3.0 * kCK2 * pinvsq * n0pp_ref;
    double temp2_ref = temp1_ref2 * kCK2 * pinvsq;
    double temp3_ref = -(15.0/32.0) * our_J4 * pinvsq * pinvsq * n0pp_ref;

    double mdot_ref = n0pp_ref + 0.5*temp1_ref2*betao*x3thm1
        + (1.0/16.0)*temp2_ref*betao*(13.0 + theta2*(-78.0 + 137.0*theta2));
    double omgdot_ref = -0.5*temp1_ref2*(1.0-5.0*theta2)
        + (1.0/16.0)*temp2_ref*(7.0 + theta2*(-114.0 + 395.0*theta2))
        + temp3_ref*(3.0 + theta2*(-36.0 + 49.0*theta2));
    double nodedot_ref = -temp1_ref2*cosio
        + 0.5*temp2_ref*cosio*(4.0 - 19.0*theta2)
        + 2.0*temp3_ref*cosio*(3.0 - 7.0*theta2);
    compare("M_dot", double(ns.M_dot.value), mdot_ref);
    compare("omega_dot", double(ns.omega_dot.value), omgdot_ref);
    compare("Omega_dot", double(ns.Omega_dot.value), nodedot_ref);

    // =====================================================
    // Propagation at t=360 step-by-step
    // =====================================================
    std::cout << "\n=== DIAGNOSTIC: Propagation at t=360 min ===\n";
    std::cout << std::setw(24) << "Name" << std::setw(18) << "Computed"
              << std::setw(18) << "Reference" << std::setw(14) << "AbsDiff"
              << std::setw(14) << "RelDiff" << "\n";
    std::cout << std::string(82, '-') << "\n";

    double t360 = 360.0;
    double t360sq = t360*t360;
    double t360cu = t360sq*t360;
    double t360fo = t360cu*t360;

    // Secular advances
    double xmdf_360 = M0 + mdot_ref*t360;
    double omgadf_360 = omega0_val + omgdot_ref*t360;
    double xnoddf_360 = Omega0_val + nodedot_ref*t360;
    double xnode_360 = xnoddf_360 + xnodcf_ref*t360sq;

    // Drag accumulation
    double tempa_360 = 1.0 - C1_ref*t360;
    double tempe_360 = td.bstar*C4_ref*t360;
    double templ_360 = t2cof_ref*t360sq;

    // Non-simple corrections
    bool simple = (perigee_ref < 220.0);
    double xmp_360, omega_360;
    if (!simple) {
        double delomg_360 = omgcof_ref*t360;
        double delm_360 = xmcof_ref*(std::pow(1.0 + eta_ref2*std::cos(xmdf_360), 3.0) - delmo_ref);
        double temp_drag_360 = delomg_360 + delm_360;
        xmp_360 = xmdf_360 + temp_drag_360;
        omega_360 = omgadf_360 - temp_drag_360;

        tempa_360 -= D2_ref*t360sq + D3_ref*t360cu + D4_ref*t360fo;
        tempe_360 += td.bstar*C5_ref*(std::sin(xmp_360) - sinmo_ref);
        templ_360 += t3cof_ref*t360cu + t360fo*(t4cof_ref + t360*t5cof_ref);
    } else {
        xmp_360 = xmdf_360;
        omega_360 = omgadf_360;
    }

    std::cout << std::setprecision(12);
    std::cout << "  tempa@360 = " << tempa_360 << " (semi-major axis drag factor)\n";
    std::cout << "  tempe@360 = " << tempe_360 << " (eccentricity correction)\n";
    std::cout << "  templ@360 = " << templ_360 << " (mean longitude correction)\n";

    double a_360 = a0pp_ref * tempa_360 * tempa_360;
    double e_360 = e0 - tempe_360;
    if (e_360 < 1e-6) e_360 = 1e-6;

    double xl_360 = xmp_360 + omega_360 + xnode_360 + n0pp_ref*templ_360;

    std::cout << "  a@360 = " << a_360 << " ER\n";
    std::cout << "  e@360 = " << e_360 << "\n";

    // Now propagate using our code and compare intermediates
    math::TrackedValue<T> t360_tv = math::TrackedValue<T>(T(360.0), T(0), T(0), T(0));
    sgp4::StateVector<T> sv_360 = sgp4::propagate_near_space(ns, config, t360_tv, T(1e-12));

    double ref_x360 = -7154.03120202, ref_y360 = -3783.17682504, ref_z360 = -3536.19412294;
    double dx360 = sv_360.position_km.x.value - ref_x360;
    double dy360 = sv_360.position_km.y.value - ref_y360;
    double dz360 = sv_360.position_km.z.value - ref_z360;
    double pe360 = std::sqrt(dx360*dx360 + dy360*dy360 + dz360*dz360);

    std::cout << "\n  Our position at t=360: "
              << std::setprecision(8) << std::fixed
              << sv_360.position_km.x.value << "  "
              << sv_360.position_km.y.value << "  "
              << sv_360.position_km.z.value << "\n";
    std::cout << "  Reference at t=360:   "
              << ref_x360 << "  " << ref_y360 << "  " << ref_z360 << "\n";
    std::cout << "  Position error: " << std::scientific << pe360 << " km\n";

    // Reference propagation using our reference coefficients
    double beta2_360 = 1.0 - e_360*e_360;
    double n_360 = kXKE / std::pow(a_360, 1.5);
    double axN_360 = e_360 * std::cos(omega_360);
    double temp_lp_360 = 1.0/(a_360*beta2_360);
    double xlcof_ref = 0.125 * kA3OVK2 * sinio_ref * (3.0+5.0*cosio)/(1.0+cosio);
    double aycof_ref = 0.25 * kA3OVK2 * sinio_ref;
    double xll_360 = temp_lp_360 * xlcof_ref * axN_360;
    double aynl_360 = temp_lp_360 * aycof_ref;
    double xlt_360 = xl_360 + xll_360;
    double ayn_360 = e_360*std::sin(omega_360) + aynl_360;

    // Modified Kepler solve at t=360
    double capu_360 = std::fmod(xlt_360 - xnode_360, 2.0*3.14159265358979);
    if (capu_360 < 0) capu_360 += 2.0*3.14159265358979;

    double epw_360 = capu_360;
    for (int iter = 0; iter < 30; ++iter) {
        double sx = std::sin(epw_360);
        double cx = std::cos(epw_360);
        double f_kep = capu_360 - ayn_360*cx + axN_360*sx - epw_360;
        double fp_kep = -ayn_360*sx - axN_360*cx + 1.0;
        double delta_kep = f_kep / fp_kep;
        epw_360 += delta_kep;
        if (std::abs(delta_kep) < 1e-12) break;
    }

    double sinepw_360 = std::sin(epw_360);
    double cosepw_360 = std::cos(epw_360);
    double ecose_360 = axN_360*cosepw_360 + ayn_360*sinepw_360;
    double esine_360 = axN_360*sinepw_360 - ayn_360*cosepw_360;
    double el2_360 = axN_360*axN_360 + ayn_360*ayn_360;
    double pl_360 = a_360*(1.0 - el2_360);
    double r_360 = a_360*(1.0 - ecose_360);
    double rdot_360 = kXKE*std::sqrt(a_360)*esine_360/r_360;
    double rfdot_360 = kXKE*std::sqrt(pl_360)/r_360;

    double betal_360 = std::sqrt(1.0 - el2_360);
    double temp_k360 = 1.0/(1.0 + betal_360);
    double cosu_360 = (a_360/r_360)*(cosepw_360 - axN_360 + ayn_360*esine_360*temp_k360);
    double sinu_360 = (a_360/r_360)*(sinepw_360 - ayn_360 - axN_360*esine_360*temp_k360);
    double u_360 = std::atan2(sinu_360, cosu_360);
    double sin2u_360 = 2.0*sinu_360*cosu_360;
    double cos2u_360 = 2.0*cosu_360*cosu_360 - 1.0;

    // Short-period corrections
    double sp2_360 = kCK2/pl_360;
    double sp3_360 = sp2_360/pl_360;
    double rk_360 = r_360*(1.0 - 1.5*sp3_360*betal_360*x3thm1) + 0.5*sp2_360*(1.0-theta2)*cos2u_360;
    double uk_360 = u_360 - 0.25*sp3_360*(7.0*theta2-1.0)*sin2u_360;
    double xnodek_360 = xnode_360 + 1.5*sp3_360*cosio*sin2u_360;
    double xinck_360 = i0 + 1.5*sp3_360*cosio*sinio_ref*cos2u_360;
    double rdotk_360 = rdot_360 - n_360*sp2_360*(1.0-theta2)*sin2u_360;
    double rfdotk_360 = rfdot_360 + n_360*sp2_360*((1.0-theta2)*cos2u_360 + 1.5*x3thm1);

    // Position
    double sinuk_360 = std::sin(uk_360), cosuk_360 = std::cos(uk_360);
    double sinik_360 = std::sin(xinck_360), cosik_360 = std::cos(xinck_360);
    double sinnok_360 = std::sin(xnodek_360), cosnok_360 = std::cos(xnodek_360);

    double xmx_360 = -sinnok_360*cosik_360;
    double xmy_360 = cosnok_360*cosik_360;
    double ux_360 = xmx_360*sinuk_360 + cosnok_360*cosuk_360;
    double uy_360 = xmy_360*sinuk_360 + sinnok_360*cosuk_360;
    double uz_360 = sinik_360*sinuk_360;

    double x_360 = rk_360*ux_360*kXKMPER;
    double y_360 = rk_360*uy_360*kXKMPER;
    double z_360 = rk_360*uz_360*kXKMPER;

    double ref_pe = std::sqrt((x_360-ref_x360)*(x_360-ref_x360)
                             +(y_360-ref_y360)*(y_360-ref_y360)
                             +(z_360-ref_z360)*(z_360-ref_z360));

    std::cout << "\n  Reference calc at t=360: "
              << std::setprecision(8) << std::fixed
              << x_360 << "  " << y_360 << "  " << z_360 << "\n";
    std::cout << "  Ref calc error vs tcppver: " << std::scientific << ref_pe << " km\n";

    // Velocity (using XKMPER/60, NOT xke)
    double vkmpersec_ref = kXKMPER / 60.0;
    double vx_360 = xmx_360*cosuk_360 - cosnok_360*sinuk_360;
    double vy_360 = xmy_360*cosuk_360 - sinnok_360*sinuk_360;
    double vz_360 = sinik_360*cosuk_360;
    double xdot_360 = (rdotk_360*ux_360 + rfdotk_360*vx_360)*vkmpersec_ref;
    double ydot_360 = (rdotk_360*uy_360 + rfdotk_360*vy_360)*vkmpersec_ref;
    double zdot_360 = (rdotk_360*uz_360 + rfdotk_360*vz_360)*vkmpersec_ref;

    double ref_vx360 = 4.741887409, ref_vy360 = -4.151817765, ref_vz360 = -2.093935425;
    double vel_pe = std::sqrt((xdot_360-ref_vx360)*(xdot_360-ref_vx360)
                             +(ydot_360-ref_vy360)*(ydot_360-ref_vy360)
                             +(zdot_360-ref_vz360)*(zdot_360-ref_vz360));

    std::cout << "  Ref calc vel at t=360: "
              << std::setprecision(8) << std::fixed
              << xdot_360 << "  " << ydot_360 << "  " << zdot_360 << "\n";
    std::cout << "  Ref calc vel error: " << std::scientific << vel_pe << " km/s\n";
}
