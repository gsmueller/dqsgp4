"""
Step 4: Three-way comparison — {p1, S1} vs D(dS1/dl) vs BH61 Eq(14).

Determines whether the Poisson bracket {p1, S1} equals D(dS1/dl),
and which (if either) matches BH61 Eq(14).
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

MU = 1.0

def solve_kepler(M, e, tol=1e-15):
    E = M + e * sin(M)
    for _ in range(50):
        dE = (M - E + e * sin(E)) / (1 - e * cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def true_from_eccentric(E, e):
    return 2 * arctan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))

def orbital_from_delaunay(L, G, H, l, g, h, mu=MU):
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    return a, e, eta, theta, E, f, r


# === Functions ===

def p1_func(L, G, H, l, g, h, mu=MU):
    """p1 = L*(2a/r - 1)."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    E = solve_kepler(l, e)
    r = a * (1 - e * cos(E))
    return L * (2*a/r - 1)

def S1_func(L, G, H, l, g, h, mu=MU):
    """S1 generating function from Phase A."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
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

def dS1_dl_ours(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    Gamma = mu**2 / (a**3 * eta**3)
    return Gamma * B0 - mu**2 * (B0 + B1p * cos(2*f + 2*g)) / r**3


# === Poisson bracket {p1, S1} ===

def poisson_bracket_p1_S1(L, G, H, l, g, h, mu=MU, eps=1e-7):
    """{p1, S1} = sum_j [dp1/dl_j * dS1/dL_j - dp1/dL_j * dS1/dl_j]."""
    vars_base = [L, G, H, l, g, h]
    result = 0.0
    for j in range(3):
        Lj_idx = j
        lj_idx = j + 3

        # dp1/dl_j
        vp = list(vars_base); vp[lj_idx] += eps
        vm = list(vars_base); vm[lj_idx] -= eps
        dp1_dlj = (p1_func(*vp, mu) - p1_func(*vm, mu)) / (2*eps)

        # dS1/dL_j
        vp = list(vars_base); vp[Lj_idx] += eps
        vm = list(vars_base); vm[Lj_idx] -= eps
        dS1_dLj = (S1_func(*vp, mu) - S1_func(*vm, mu)) / (2*eps)

        # dp1/dL_j
        dp1_dLj = (p1_func(*vp, mu) - p1_func(*vm, mu)) / (2*eps)

        # dS1/dl_j
        vp = list(vars_base); vp[lj_idx] += eps
        vm = list(vars_base); vm[lj_idx] -= eps
        dS1_dlj = (S1_func(*vp, mu) - S1_func(*vm, mu)) / (2*eps)

        result += dp1_dlj * dS1_dLj - dp1_dLj * dS1_dlj

    return result


# === D operator on dS1/dl ===

def D_of_dS1_dl(L, G, H, l, g, h, mu=MU, eps=1e-7):
    """D(dS1/dl) via finite differences."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H
    sinE = sin(E); sinf = sin(f)
    q1 = 2*e*sinE + 2*eta*sinf/e if e > 1e-10 else 0.0
    q2 = -2*sinf/e if e > 1e-10 else 0.0
    q3 = 0.0

    dFdL = (dS1_dl_ours(L+eps,G,H,l,g,h,mu) - dS1_dl_ours(L-eps,G,H,l,g,h,mu)) / (2*eps)
    dFdG = (dS1_dl_ours(L,G+eps,H,l,g,h,mu) - dS1_dl_ours(L,G-eps,H,l,g,h,mu)) / (2*eps)
    dFdH = (dS1_dl_ours(L,G,H+eps,l,g,h,mu) - dS1_dl_ours(L,G,H-eps,l,g,h,mu)) / (2*eps)
    dFdl = (dS1_dl_ours(L,G,H,l+eps,g,h,mu) - dS1_dl_ours(L,G,H,l-eps,g,h,mu)) / (2*eps)
    dFdg = (dS1_dl_ours(L,G,H,l,g+eps,h,mu) - dS1_dl_ours(L,G,H,l,g-eps,h,mu)) / (2*eps)
    dFdh = (dS1_dl_ours(L,G,H,l,g,h+eps,mu) - dS1_dl_ours(L,G,H,l,g,h-eps,mu)) / (2*eps)

    return -p1*dFdL - p2*dFdG - p3*dFdH + q1*dFdl + q2*dFdg + q3*dFdh


# === BH61 Eq(14) closed form ===

def bh61_eq14(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    n = mu**2 / L**3
    rho = a / r
    bracket = (B0 * (-1/eta**3 + rho**3*(1-2*rho))
               + B1p * rho**3*(1-2*rho)*cos(2*g+2*f))
    return 3 * n * bracket


# === Our closed form ===

def our_closed_form(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    Gamma = mu**2 / (a**3 * eta**3)
    return 6 * Gamma * B0 * (a / r)


# ============================================================
# Three-way comparison
# ============================================================
print("=" * 120)
print("THREE-WAY COMPARISON: {p1,S1} vs D(dS1/dl) vs BH61 Eq(14) vs Our closed form")
print("=" * 120)
print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>5} | {'PB{p1,S1}':>16} {'D(dS1/dl)':>16} {'BH61':>16} {'6GB0r':>16} | {'PB/D':>10} {'PB/BH61':>10} {'D/6GB0r':>10}")
print("-" * 130)

a_val = 1.0
mu = MU
n_match_pb_D = 0
n_match_pb_bh = 0
n_match_D_cf = 0
n_total = 0

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    L = sqrt(mu * a_val)
    G = L * eta_val
    for I_deg in [30, 60, 85]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for g_deg in [0, 45, 90]:
            g_val = np.radians(g_deg)
            h_val = 0.0
            for l_val in [0.5, 1.5, 3.0]:
                pb = poisson_bracket_p1_S1(L, G, H, l_val, g_val, h_val, mu)
                D_val = D_of_dS1_dl(L, G, H, l_val, g_val, h_val, mu)
                bh = bh61_eq14(L, G, H, l_val, g_val, h_val, mu)
                cf = our_closed_form(L, G, H, l_val, g_val, h_val, mu)

                ratio_pb_D = pb / D_val if abs(D_val) > 1e-30 else float('nan')
                ratio_pb_bh = pb / bh if abs(bh) > 1e-30 else float('nan')
                ratio_D_cf = D_val / cf if abs(cf) > 1e-30 else float('nan')

                d1 = max(abs(pb), abs(D_val), 1e-30)
                d2 = max(abs(pb), abs(bh), 1e-30)
                d3 = max(abs(D_val), abs(cf), 1e-30)

                if abs(pb - D_val) / d1 < 1e-3: n_match_pb_D += 1
                if abs(pb - bh) / d2 < 1e-3: n_match_pb_bh += 1
                if abs(D_val - cf) / d3 < 1e-3: n_match_D_cf += 1
                n_total += 1

                print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.1f} | "
                      f"{pb:16.8e} {D_val:16.8e} {bh:16.8e} {cf:16.8e} | "
                      f"{ratio_pb_D:10.4f} {ratio_pb_bh:10.4f} {ratio_D_cf:10.4f}")

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Total test cases: {n_total}")
print(f"  {{p1,S1}} matches D(dS1/dl):        {n_match_pb_D}/{n_total}")
print(f"  {{p1,S1}} matches BH61 Eq(14):      {n_match_pb_bh}/{n_total}")
print(f"  D(dS1/dl) matches 6*Gamma*B0*(a/r): {n_match_D_cf}/{n_total}")
print()

if n_match_pb_D == n_total:
    print(">>> {p1,S1} = D(dS1/dl)  [identity confirmed]")
elif n_match_pb_bh == n_total:
    print(">>> {p1,S1} = BH61 Eq(14)  [but != D(dS1/dl), so identity delta_p = D(dS1/dl) is wrong]")
else:
    print(">>> {p1,S1} matches neither D(dS1/dl) nor BH61 consistently")
    print("    This suggests {p1,S1} and D(dS1/dl) are fundamentally different quantities.")

if n_match_D_cf == n_total:
    print(">>> D(dS1/dl) = 6*Gamma*B0*(a/r)  [cleanroom result confirmed]")
