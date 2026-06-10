"""
Numerical verification of Poisson brackets {p_j, S_1} in Delaunay variables.
Phase A2 cleanroom derivation verification.

Computes:
  {p1, S1} by analytical assembly and by full finite differences
  {p2, S1} = -dS1/dg (analytical) vs FD
  {p3, S1} = -dS1/dh = 0 (analytical) vs FD
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

MU = 1.0  # gravitational parameter (normalized)

# ============================================================
# Core orbital mechanics
# ============================================================

def solve_kepler(M, e, tol=1e-15):
    """Solve Kepler's equation M = E - e*sin(E) for E."""
    E = M + e * sin(M)
    for _ in range(100):
        dE = (M - E + e * sin(E)) / (1 - e * cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def true_from_eccentric(E, e):
    """Compute true anomaly from eccentric anomaly."""
    return 2 * arctan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))

def delaunay_to_orbital(L, G, H, l, g, h, mu=MU):
    """Convert Delaunay variables to orbital elements, eccentric anomaly, true anomaly, radius."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    return a, e, eta, theta, E, f, r

# ============================================================
# S1 and p_j in Delaunay variables
# ============================================================

def S1_delaunay(L, G, H, l, g, h, mu=MU):
    """Compute S1 from Delaunay variables."""
    a, e, eta, theta, E, f, r = delaunay_to_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B0 = -0.5 + 1.5 * theta**2
    s2 = 1 - theta**2
    phi = f - l
    return -Gamma * (
        B0 * (phi + e * sin(f))
        + 0.75 * s2 * sin(2*f + 2*g)
        + 0.75 * e * s2 * sin(f + 2*g)
        + 0.25 * e * s2 * sin(3*f + 2*g)
    )

def p1_delaunay(L, G, H, l, g, h, mu=MU):
    """p1 = L(2a/r - 1)."""
    a = L**2 / mu
    e = sqrt(max(1 - (G/L)**2, 0))
    E = solve_kepler(l, e)
    r = a * (1 - e * cos(E))
    return L * (2 * a / r - 1)

def p2_delaunay(L, G, H, l, g, h, mu=MU):
    return G

def p3_delaunay(L, G, H, l, g, h, mu=MU):
    return H

# ============================================================
# Poisson bracket by full finite differences
# ============================================================

def poisson_bracket_fd(func_p, func_S, L, G, H, l, g, h, eps=1e-7):
    """
    Compute {func_p, func_S} by central finite differences in Delaunay variables.
    {p, S} = sum_j (dp/dl_j * dS/dL_j - dp/dL_j * dS/dl_j)
    """
    momenta = [L, G, H]
    angles = [l, g, h]
    result = 0.0
    for j in range(3):
        # dp/dl_j
        a_p = list(momenta) + list(angles); a_m = list(momenta) + list(angles)
        a_p[3+j] += eps; a_m[3+j] -= eps
        dp_dlj = (func_p(*a_p) - func_p(*a_m)) / (2*eps)

        # dS/dL_j
        a_p = list(momenta) + list(angles); a_m = list(momenta) + list(angles)
        a_p[j] += eps; a_m[j] -= eps
        dS_dLj = (func_S(*a_p) - func_S(*a_m)) / (2*eps)

        # dp/dL_j
        a_p = list(momenta) + list(angles); a_m = list(momenta) + list(angles)
        a_p[j] += eps; a_m[j] -= eps
        dp_dLj = (func_p(*a_p) - func_p(*a_m)) / (2*eps)

        # dS/dl_j
        a_p = list(momenta) + list(angles); a_m = list(momenta) + list(angles)
        a_p[3+j] += eps; a_m[3+j] -= eps
        dS_dlj = (func_S(*a_p) - func_S(*a_m)) / (2*eps)

        result += dp_dlj * dS_dLj - dp_dLj * dS_dlj
    return result

# ============================================================
# Analytical {p1, S1} via assembled partial derivatives
# ============================================================

def bracket_p1_S1_analytical(L, G, H, l, g, h, mu=MU):
    """
    {p1, S1} = dp1/dl * dS1/dL - dp1/dL * dS1/dl - dp1/dG * dS1/dg.

    All partial derivatives computed analytically with full chain rules.
    """
    a, e, eta, theta, E, f, r = delaunay_to_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    s2 = 1 - theta**2
    phi = f - l

    # ========== Partial derivatives of p1 ==========

    # dp1/dl: p1 = L(2a/r - 1), r changes via Kepler
    dr_dl = a * e * sin(E) / (1 - e * cos(E))  # = ae*sin(f)/eta
    dp1_dl = -2 * L * a / r**2 * dr_dl

    # dp1/dL: p1 = 2L^3/(mu*r) - L
    de_dL = eta**2 / (e * L) if e > 1e-14 else 0.0
    dE_dL = de_dL * sin(E) / (1 - e * cos(E))
    da_dL = 2 * a / L
    dr_dL = da_dL * (1 - e*cos(E)) + a * (-de_dL*cos(E) + e*sin(E)*dE_dL)
    dp1_dL = (2*a/r - 1) + L * (2*da_dL/r - 2*a*dr_dL/r**2)

    # dp1/dG: a is fixed, but e changes
    de_dG = -eta / (e * L) if e > 1e-14 else 0.0
    dE_dG = de_dG * sin(E) / (1 - e * cos(E))
    dr_dG = a * (-de_dG * cos(E) + e * sin(E) * dE_dG)
    dp1_dG = -2 * L * a / r**2 * dr_dG

    # ========== Partial derivatives of S1 ==========

    # dS1/dl = (1/n)(H1 - <H1>)
    dS1_dl = mu**2 * B0 / (a**3 * eta**3) - mu**2 / r**3 * (B0 + B1p * cos(2*f + 2*g))

    # dS1/dg
    dS1_dg = -Gamma * B1p * (cos(2*f+2*g) + e*cos(f+2*g) + (e/3)*cos(3*f+2*g))

    # dS1/dL (the hardest one)
    # S1 = -Gamma * F where F = B0*(phi+e*sin(f)) + (3/4)*s2*sin(2f+2g) + (3/4)*e*s2*sin(f+2g) + (1/4)*e*s2*sin(3f+2g)

    # dGamma/dL = -3*Gamma/L
    dGamma_dL = -3 * Gamma / L

    # df/dL via chain rule: df/dE * dE/dL + df/de * de/dL
    df_dE = a * eta / r  # standard relation
    # df/de at constant E:
    if abs(sin(E)) < 1e-14 and abs(sin(f)) < 1e-14:
        df_de = 0.0
    else:
        df_de = sin(E) / (eta * (1 - e*cos(E)))
    df_dL = df_dE * dE_dL + df_de * de_dL

    # dF/dL (B0, s2 don't depend on L since theta = H/G)
    dF_dL = (B0 * (df_dL + de_dL*sin(f) + e*cos(f)*df_dL)
             + 0.75*s2 * 2*cos(2*f+2*g)*df_dL
             + 0.75*s2 * (de_dL*sin(f+2*g) + e*cos(f+2*g)*df_dL)
             + 0.25*s2 * (de_dL*sin(3*f+2*g) + 3*e*cos(3*f+2*g)*df_dL))

    F = (B0 * (phi + e*sin(f))
         + 0.75*s2*sin(2*f+2*g)
         + 0.75*e*s2*sin(f+2*g)
         + 0.25*e*s2*sin(3*f+2*g))

    dS1_dL = -dGamma_dL * F - Gamma * dF_dL

    # ========== Assemble ==========
    return dp1_dl * dS1_dL - dp1_dL * dS1_dl - dp1_dG * dS1_dg

# ============================================================
# Analytical {p2, S1} = -dS1/dg
# ============================================================

def bracket_p2_S1_analytical(L, G, H, l, g, h, mu=MU):
    """{G, S1} = -dS1/dg."""
    a, e, eta, theta, E, f, r = delaunay_to_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B1p = 1.5 * (1 - theta**2)
    return Gamma * B1p * (cos(2*f+2*g) + e*cos(f+2*g) + (e/3)*cos(3*f+2*g))

# ============================================================
# Main verification
# ============================================================

def main():
    print("=" * 90)
    print("PHASE A2: POISSON BRACKET VERIFICATION")
    print("=" * 90)

    e_vals = [0.01, 0.1, 0.3]
    I_vals = [30, 60, 85]
    g_vals = [0, 45, 90]
    l_vals = [0.5, 1.5, 3.0]

    # ---- Verify individual partial derivatives ----
    print("\n--- Spot-check: analytical partials vs FD (e=0.1, I=60, g=45, l=1.5) ---")
    e_t, I_t, g_t, l_t = 0.1, 60, 45, 1.5
    eta_t = sqrt(1 - e_t**2)
    L_t = sqrt(MU * 1.0)
    G_t = L_t * eta_t
    H_t = G_t * cos(np.radians(I_t))
    gv = np.radians(g_t)
    eps = 1e-7

    # dp1/dl
    fd = (p1_delaunay(L_t, G_t, H_t, l_t+eps, gv, 0) - p1_delaunay(L_t, G_t, H_t, l_t-eps, gv, 0)) / (2*eps)
    a, e, eta, theta, E, f, r = delaunay_to_orbital(L_t, G_t, H_t, l_t, gv, 0)
    an = -2*L_t*a/r**2 * a*e*sin(E)/(1-e*cos(E))
    print(f"  dp1/dl: an={an:+18.10e}  fd={fd:+18.10e}  err={abs(an-fd)/max(abs(an),1e-20):.2e}")

    # dp1/dL
    fd = (p1_delaunay(L_t+eps, G_t, H_t, l_t, gv, 0) - p1_delaunay(L_t-eps, G_t, H_t, l_t, gv, 0)) / (2*eps)
    de_dL = eta**2/(e*L_t)
    dE_dL = de_dL*sin(E)/(1-e*cos(E))
    da_dL = 2*a/L_t
    dr_dL = da_dL*(1-e*cos(E)) + a*(-de_dL*cos(E) + e*sin(E)*dE_dL)
    an = (2*a/r-1) + L_t*(2*da_dL/r - 2*a*dr_dL/r**2)
    print(f"  dp1/dL: an={an:+18.10e}  fd={fd:+18.10e}  err={abs(an-fd)/max(abs(an),1e-20):.2e}")

    # dp1/dG
    fd = (p1_delaunay(L_t, G_t+eps, H_t, l_t, gv, 0) - p1_delaunay(L_t, G_t-eps, H_t, l_t, gv, 0)) / (2*eps)
    de_dG = -eta/(e*L_t)
    dE_dG = de_dG*sin(E)/(1-e*cos(E))
    dr_dG = a*(-de_dG*cos(E) + e*sin(E)*dE_dG)
    an = -2*L_t*a/r**2*dr_dG
    print(f"  dp1/dG: an={an:+18.10e}  fd={fd:+18.10e}  err={abs(an-fd)/max(abs(an),1e-20):.2e}")

    # dS1/dl
    fd = (S1_delaunay(L_t, G_t, H_t, l_t+eps, gv, 0) - S1_delaunay(L_t, G_t, H_t, l_t-eps, gv, 0)) / (2*eps)
    B0 = -0.5+1.5*theta**2; B1p = 1.5*(1-theta**2)
    an = MU**2*B0/(a**3*eta**3) - MU**2/r**3*(B0 + B1p*cos(2*f+2*gv))
    print(f"  dS1/dl: an={an:+18.10e}  fd={fd:+18.10e}  err={abs(an-fd)/max(abs(an),1e-20):.2e}")

    # dS1/dg
    fd = (S1_delaunay(L_t, G_t, H_t, l_t, gv+eps, 0) - S1_delaunay(L_t, G_t, H_t, l_t, gv-eps, 0)) / (2*eps)
    Gamma = MU**2/(a**3*eta**3)
    an = -Gamma*B1p*(cos(2*f+2*gv) + e*cos(f+2*gv) + (e/3)*cos(3*f+2*gv))
    print(f"  dS1/dg: an={an:+18.10e}  fd={fd:+18.10e}  err={abs(an-fd)/max(abs(an),1e-20):.2e}")

    # dS1/dL
    fd = (S1_delaunay(L_t+eps, G_t, H_t, l_t, gv, 0) - S1_delaunay(L_t-eps, G_t, H_t, l_t, gv, 0)) / (2*eps)
    an_dS1dL = bracket_p1_S1_analytical.__wrapped_dS1_dL(L_t, G_t, H_t, l_t, gv, 0) if hasattr(bracket_p1_S1_analytical, '__wrapped_dS1_dL') else None
    print(f"  dS1/dL: fd={fd:+18.10e}")

    # dS1/dH
    fd = (S1_delaunay(L_t, G_t, H_t+eps, l_t, gv, 0) - S1_delaunay(L_t, G_t, H_t-eps, l_t, gv, 0)) / (2*eps)
    print(f"  dS1/dH: fd={fd:+18.10e}")

    # ---- {p1, S1} verification ----
    print("\n" + "=" * 90)
    print("{p1, S1}: ANALYTICAL vs FINITE DIFFERENCES")
    print("=" * 90)
    print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'Analytical':>18} {'FD':>18} {'RelErr':>12} {'Status':>6}")
    print("-" * 90)

    max_err_p1 = 0.0
    pass_p1 = 0; fail_p1 = 0
    for e_val in e_vals:
        for I_deg in I_vals:
            for g_deg in g_vals:
                for l_val in l_vals:
                    eta_v = sqrt(1-e_val**2)
                    L = sqrt(MU)
                    G = L * eta_v
                    H = G * cos(np.radians(I_deg))
                    g = np.radians(g_deg)
                    h = 0.0

                    an = bracket_p1_S1_analytical(L, G, H, l_val, g, h)
                    fd = poisson_bracket_fd(p1_delaunay, S1_delaunay, L, G, H, l_val, g, h)

                    ref = max(abs(an), abs(fd), 1e-20)
                    err = abs(an - fd) / ref
                    max_err_p1 = max(max_err_p1, err)
                    status = "PASS" if err < 1e-5 else "FAIL"
                    if status == "PASS": pass_p1 += 1
                    else: fail_p1 += 1

                    print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | {an:+18.10e} {fd:+18.10e} {err:12.2e} {status:>6}")

    print(f"\n{{p1, S1}} summary: {pass_p1} PASS, {fail_p1} FAIL, max err = {max_err_p1:.2e}")

    # ---- {p2, S1} verification ----
    print("\n" + "=" * 90)
    print("{p2, S1} = -dS1/dg: ANALYTICAL vs FINITE DIFFERENCES")
    print("=" * 90)
    print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'Analytical':>18} {'FD':>18} {'RelErr':>12} {'Status':>6}")
    print("-" * 90)

    max_err_p2 = 0.0
    pass_p2 = 0; fail_p2 = 0
    for e_val in e_vals:
        for I_deg in I_vals:
            for g_deg in g_vals:
                for l_val in l_vals:
                    eta_v = sqrt(1-e_val**2)
                    L = sqrt(MU)
                    G = L * eta_v
                    H = G * cos(np.radians(I_deg))
                    g = np.radians(g_deg)
                    h = 0.0

                    an = bracket_p2_S1_analytical(L, G, H, l_val, g, h)
                    fd = poisson_bracket_fd(p2_delaunay, S1_delaunay, L, G, H, l_val, g, h)

                    ref = max(abs(an), abs(fd), 1e-20)
                    err = abs(an - fd) / ref
                    max_err_p2 = max(max_err_p2, err)
                    status = "PASS" if err < 1e-5 else "FAIL"
                    if status == "PASS": pass_p2 += 1
                    else: fail_p2 += 1

                    print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | {an:+18.10e} {fd:+18.10e} {err:12.2e} {status:>6}")

    print(f"\n{{p2, S1}} summary: {pass_p2} PASS, {fail_p2} FAIL, max err = {max_err_p2:.2e}")

    # ---- {p3, S1} verification ----
    print("\n" + "=" * 90)
    print("{p3, S1} = -dS1/dh = 0: VERIFICATION")
    print("=" * 90)

    max_p3 = 0.0
    for e_val in e_vals:
        for I_deg in I_vals:
            for g_deg in g_vals:
                for l_val in l_vals:
                    eta_v = sqrt(1-e_val**2)
                    L = sqrt(MU)
                    G = L * eta_v
                    H = G * cos(np.radians(I_deg))
                    g = np.radians(g_deg)
                    h = 0.0
                    fd = poisson_bracket_fd(p3_delaunay, S1_delaunay, L, G, H, l_val, g, h)
                    max_p3 = max(max_p3, abs(fd))

    status = "PASS" if max_p3 < 1e-10 else "FAIL"
    print(f"Max |{{p3, S1}}| over all {len(e_vals)*len(I_vals)*len(g_vals)*len(l_vals)} test cases: {max_p3:.2e} [{status}]")

    # ---- Decompose {p1, S1} into three terms ----
    print("\n" + "=" * 90)
    print("DECOMPOSITION: {p1, S1} = T1 + T2 + T3")
    print("T1 = dp1/dl * dS1/dL,  T2 = -dp1/dL * dS1/dl,  T3 = -dp1/dG * dS1/dg")
    print("=" * 90)
    print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>4} | {'T1':>14} {'T2':>14} {'T3':>14} | {'Sum':>14}")
    print("-" * 90)

    for e_val in [0.01, 0.1, 0.3]:
        for I_deg in [30, 60, 85]:
            for g_deg in [45]:
                for l_val in [1.5]:
                    eta_v = sqrt(1-e_val**2)
                    L = sqrt(MU)
                    G = L * eta_v
                    H = G * cos(np.radians(I_deg))
                    gv = np.radians(g_deg)
                    h = 0.0

                    a, e, eta, theta, E, f, r = delaunay_to_orbital(L, G, H, l_val, gv, h)
                    Gamma = MU**2/(a**3*eta**3)
                    B0 = -0.5+1.5*theta**2
                    B1p = 1.5*(1-theta**2)
                    s2 = 1-theta**2

                    # dp1/dl
                    dp1dl = -2*L*a/r**2 * a*e*sin(E)/(1-e*cos(E))

                    # dp1/dL
                    de_dL = eta**2/(e*L) if e > 1e-14 else 0
                    dEdL = de_dL*sin(E)/(1-e*cos(E))
                    dadL = 2*a/L
                    drdL = dadL*(1-e*cos(E)) + a*(-de_dL*cos(E)+e*sin(E)*dEdL)
                    dp1dL = (2*a/r-1) + L*(2*dadL/r - 2*a*drdL/r**2)

                    # dp1/dG
                    de_dG = -eta/(e*L) if e > 1e-14 else 0
                    dEdG = de_dG*sin(E)/(1-e*cos(E))
                    drdG = a*(-de_dG*cos(E)+e*sin(E)*dEdG)
                    dp1dG = -2*L*a/r**2*drdG

                    # dS1/dl
                    dS1dl = MU**2*B0/(a**3*eta**3) - MU**2/r**3*(B0+B1p*cos(2*f+2*gv))

                    # dS1/dg
                    dS1dg = -Gamma*B1p*(cos(2*f+2*gv)+e*cos(f+2*gv)+(e/3)*cos(3*f+2*gv))

                    # dS1/dL by FD
                    dS1dL = (S1_delaunay(L+eps, G, H, l_val, gv, h) - S1_delaunay(L-eps, G, H, l_val, gv, h))/(2*eps)

                    T1 = dp1dl * dS1dL
                    T2 = -dp1dL * dS1dl
                    T3 = -dp1dG * dS1dg

                    print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:4.1f} | {T1:+14.8f} {T2:+14.8f} {T3:+14.8f} | {T1+T2+T3:+14.8f}")

    # ---- Compare with Phase A candidate ----
    print("\n" + "=" * 90)
    print("COMPARISON WITH PHASE A CANDIDATE delta_p1 = 6*mu^2*B0/(a^2*r*eta^3) - 2*mu^2*B1'/(e*r^3)*[cos(f+2g)-cos(3f+2g)]")
    print("=" * 90)
    print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>4} | {'PB_FD':>14} {'PhaseA':>14} {'Ratio':>10}")
    print("-" * 80)

    for e_val in [0.01, 0.1, 0.3]:
        for I_deg in [30, 60, 85]:
            for g_deg in [0, 45, 90]:
                for l_val in [0.5, 1.5, 3.0]:
                    eta_v = sqrt(1-e_val**2)
                    L = sqrt(MU)
                    G = L * eta_v
                    H = G * cos(np.radians(I_deg))
                    gv = np.radians(g_deg)

                    pb = poisson_bracket_fd(p1_delaunay, S1_delaunay, L, G, H, l_val, gv, 0)

                    a, e, eta, theta, E, f, r = delaunay_to_orbital(L, G, H, l_val, gv, 0)
                    B0 = -0.5+1.5*theta**2; B1p = 1.5*(1-theta**2)
                    cand = 6*MU**2*B0/(a**2*r*eta**3) - 2*MU**2*B1p/(e*r**3)*(cos(f+2*gv)-cos(3*f+2*gv))

                    ratio = pb/cand if abs(cand) > 1e-20 else float('inf')
                    print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:4.1f} | {pb:+14.8f} {cand:+14.8f} {ratio:10.6f}")

    print("\n\nDONE.")

if __name__ == "__main__":
    main()
