"""
Numerical verification of ALL BH61 Eq.(14) delta_p_j and delta_q_j expressions.

Method:
  For delta_p_j = D(dS/dl_j): use ANALYTICAL dS/dl_j, apply D by single-layer FD
  For delta_q_j = D(dS/dL_j): use mixed-partial FD stencil (no nesting)

k_2 = 1 throughout.
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

def orbital_from_delaunay(L, G, H, l, g, h):
    a = L**2 / MU
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    return a, e, eta, theta, E, f, r

# =====================================================================
# S_1 + S_1* closed forms
# =====================================================================

def S_total(L, G, H, l, g, h, alpha_star=-1.0/16):
    """S = S_1 + S_1* in Brouwer normalization (k_2=1, prefactor mu^2/G^3)."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    phi = f - l
    P = MU**2 / G**3

    S1 = P * (B0 * (phi + e * sin(f))
              + B1p * (0.5*sin(2*f+2*g) + 0.5*e*sin(f+2*g) + e/6.0*sin(3*f+2*g)))

    if abs(1 - 5*theta**2) > 1e-12:
        C = (1 - 15*theta**2) * (1 - theta**2) / (1 - 5*theta**2)
        S1star = P * alpha_star * C * e**2 * sin(2*g)
    else:
        S1star = 0.0

    return S1 + S1star

# =====================================================================
# ANALYTICAL partial derivatives of S_1 + S_1*
# =====================================================================

def dS_dl_analytical(L, G, H, l, g, h, alpha_star=-1.0/16):
    """Analytical dS/dl from Brouwer (1959) Eq(13):
    dS_1/dl = (mu^2/L^3)[B0 sigma_1 + B1' sigma_2]
    where sigma_1 = (a/r)^3 - 1/eta^3, sigma_2 = (a/r)^3 cos(2g+2f).
    S_1* does not depend on l, so dS_1*/dl = 0."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    P = MU**2 / L**3
    rho3 = (a / r) ** 3
    sigma1 = rho3 - 1.0 / eta**3
    sigma2 = rho3 * cos(2*f + 2*g)
    return P * (B0 * sigma1 + B1p * sigma2)

def dS_dg_analytical(L, G, H, l, g, h, alpha_star=-1.0/16):
    """Analytical dS/dg = (mu^2/G^3) B1' [cos(2f+2g) + e cos(f+2g) + (e/3) cos(3f+2g)]
    + 2*(mu^2/G^3)*alpha_star*C*e^2*cos(2g)."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    B1p = 1.5 * (1 - theta**2)
    P = MU**2 / G**3

    dS1_dg = P * B1p * (cos(2*f+2*g) + e*cos(f+2*g) + e/3.0*cos(3*f+2*g))

    if abs(1 - 5*theta**2) > 1e-12:
        C = (1 - 15*theta**2) * (1 - theta**2) / (1 - 5*theta**2)
        dS1star_dg = P * 2 * alpha_star * C * e**2 * cos(2*g)
    else:
        dS1star_dg = 0.0

    return dS1_dg + dS1star_dg

# =====================================================================
# Single-layer D operator (operates on an analytical function)
# =====================================================================

def D_operator(F_func, L, G, H, l, g, h, eps=1e-7):
    """D(F) = sum_j [-p_j dF/dL_j + q_j dF/dl_j] by central finite differences."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)

    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H
    q1 = 2*e*sin(E) + 2*eta*sin(f)/e if e > 1e-10 else 0.0
    q2 = -2*sin(f)/e if e > 1e-10 else 0.0
    q3 = 0.0

    def dFdv(idx, eps_v):
        args = [L, G, H, l, g, h]
        args_p = list(args); args_p[idx] += eps_v
        args_m = list(args); args_m[idx] -= eps_v
        return (F_func(*args_p) - F_func(*args_m)) / (2*eps_v)

    return (-p1*dFdv(0, eps) - p2*dFdv(1, eps) - p3*dFdv(2, eps)
            + q1*dFdv(3, eps) + q2*dFdv(4, eps) + q3*dFdv(5, eps))

# =====================================================================
# Mixed-partial D operator for delta_q_j = D(dS/dL_j)
# Uses 4-point stencil: d^2S/(dv dL_j) without nesting
# =====================================================================

def D_of_dS_dLj(j_idx, L, G, H, l, g, h, alpha_star=-1.0/16, eps=1e-4):
    """Compute D(dS/dL_j) using two-layer FD with different step sizes.

    j_idx: 0 for dS/dL, 1 for dS/dG, 2 for dS/dH

    Inner layer: dS/dL_j computed with eps_inner
    Outer layer: D operator applied with eps_outer
    """
    eps_inner = eps
    eps_outer = eps * 0.1  # 10x smaller for outer

    def S_eval(LL, GG, HH, ll, gg, hh):
        return S_total(LL, GG, HH, ll, gg, hh, alpha_star)

    def dS_dLj(LL, GG, HH, ll, gg, hh):
        args_p = [LL, GG, HH, ll, gg, hh]
        args_m = [LL, GG, HH, ll, gg, hh]
        args_p[j_idx] += eps_inner
        args_m[j_idx] -= eps_inner
        return (S_eval(*args_p) - S_eval(*args_m)) / (2 * eps_inner)

    return D_operator(dS_dLj, L, G, H, l, g, h, eps=eps_outer)

# =====================================================================
# BH61 Eq.(14) formulas
# =====================================================================

def delta_p1_BH61(L, G, H, l, g, h):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    P = MU**2 / L**3
    rho = a/r
    return 3*P*(B0*(-1/eta**3 + rho**3*(1 - 2*rho)) + B1p*rho**3*(1 - 2*rho)*cos(2*g + 2*f))

def delta_p2_BH61(L, G, H, l, g, h):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    B1p = 1.5 * (1 - theta**2)
    P = MU**2 / G**3
    if abs(1 - 5*theta**2) < 1e-12:
        return 0.0
    Phi = (1 - 16*theta**2 + 15*theta**4) / (1 - 5*theta**2)
    line1 = P * B1p * (1.0/3*cos(2*g+2*f) + e*cos(2*g+f) + e/3.0*cos(2*g+3*f))
    line2 = P * e * Phi * (e/8.0*cos(2*g) + 0.5*cos(2*g+f))
    return line1 + line2

def delta_q1_BH61(L, G, H, l, g, h):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    P_L3 = MU**2 / L**3
    rho = a/r
    rho2eta2 = rho**2 * eta**2
    B0_coeff = -1 + 3*theta**2
    B1p_coeff = 1.5 - 1.5*theta**2
    eta = G/L

    if abs(e) < 1e-12 or abs(1 - 5*theta**2) < 1e-12:
        return 0.0

    # Term 1: mu^2/(e^2 L^3 G)
    pf1 = P_L3 / (e**2 * G)
    t1 = pf1 * (B0_coeff * (rho2eta2 + rho + 1) * sin(2*f)
                + B1p_coeff * ((-rho2eta2 - rho + 1) * sin(2*g)
                              + (rho2eta2 + rho + 1.0/3) * sin(2*g + 4*f)))

    # Term 2: mu^2/(e L^3 G) * (a/r)
    pf2 = P_L3 / (e * G) * rho
    t2 = pf2 * (B0_coeff * (rho2eta2 + rho + 4) * sin(f)
                + B1p_coeff * ((-rho2eta2 - rho + 2) * sin(2*g + f)
                              + (rho2eta2 + rho + 2) * sin(2*g + 3*f)))

    # Term 3: mu^2/(e L^3 G) * inclination factor
    pf3 = P_L3 / (e * G)
    incl = 0.25*(1 - 11*theta**2) - 10*theta**4/(1 - 5*theta**2)
    t3 = pf3 * incl * ((1 - 3*rho)*e*sin(2*g) + sin(2*g + f) - sin(2*g - f))

    return t1 + t2 + t3

def delta_q2_BH61(L, G, H, l, g, h):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    rho = a/r
    rho2eta2 = rho**2 * eta**2
    B0_coeff = -1 + 3*theta**2
    B1p_coeff = 1.5 - 1.5*theta**2
    phi = f - l
    sinE = sin(E)

    if abs(e) < 1e-12 or abs(1 - 5*theta**2) < 1e-12:
        return 0.0

    P_L2G2 = MU**2 / (L**2 * G**2)
    P_G4 = MU**2 / G**4

    # Groups 1-3: L^2 G^2 denominator
    pf1 = -P_L2G2 / e**2
    g1 = pf1 * (B0_coeff * (rho2eta2 + rho + 1) * sin(2*f)
                + B1p_coeff * ((-rho2eta2 - rho + 1)*sin(2*g) + (rho2eta2 + rho + 1.0/3)*sin(2*g+4*f)))

    pf2 = -P_L2G2 / e * rho
    g2 = pf2 * (B0_coeff * (rho2eta2 + rho + 4) * sin(f)
                + B1p_coeff * ((-rho2eta2 - rho + 2)*sin(2*g+f) + (rho2eta2 + rho + 2)*sin(2*g+3*f)))

    pf3 = -P_L2G2 / e * (1 - rho)
    g3 = pf3 * (B0_coeff * (rho2eta2 + rho + 1) * sin(f)
                + B1p_coeff * ((-rho2eta2 - rho - 1)*sin(2*g+f) + (rho2eta2 + rho + 1.0/3)*sin(2*g+3*f)))

    # Groups 4-7: G^4 denominator
    c369 = (-3 + 15*theta**2)
    br369 = 2*phi + e*(sin(f) - sinE) + (1-eta)*sin(f)/e
    l369 = -P_G4 * c369 * br369

    c370 = (9 - 15*theta**2)
    br370 = 1.0/3*sin(2*g+2*f) + 0.5*e*sin(2*g+f) + e/6.0*sin(2*g+3*f)
    l370 = -P_G4 * c370 * br370

    c371 = -0.5 + 5.5*theta**2 + 20*theta**4/(1-5*theta**2)
    l371 = -P_G4 * c371 * sin(2*g)

    c372 = e*(1.0/8*(1-33*theta**2) - 25*theta**4/(1-5*theta**2) - 50*theta**6/(1-5*theta**2)**2)
    l372 = -P_G4 * c372 * (sin(2*g+f) + sin(2*g-f))

    eta2 = eta**2
    ci = (1.0/8*(1-33*theta**2)*eta2 - 3.0/8*(1-55*theta**2/3.0)
          + (-25*eta2+35)*theta**4/(1-5*theta**2)
          - 50*(eta2-1)*theta**6/(1-5*theta**2)**2)
    l373 = -P_G4 * (-ci/e) * (sin(2*g+f) - sin(2*g-f))

    return g1 + g2 + g3 + l369 + l370 + l371 + l372 + l373

def delta_q3_BH61(L, G, H, l, g, h):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h)
    phi = f - l
    sinE = sin(E)
    if abs(e) < 1e-12 or abs(1 - 5*theta**2) < 1e-12:
        return 0.0
    P = MU**2 / G**4
    br375 = (2*phi + e*(sin(f)-sinE) + (1-eta)*sin(f)/e
             - 1.0/3*sin(2*g+2*f) - 0.5*e*sin(2*g+f) - e/6.0*sin(2*g+3*f))
    l375 = 6*P*theta*br375
    c376 = 11.0/8 + 10*theta**2/(1-5*theta**2) + 25*theta**4/(1-5*theta**2)**2
    l376 = -4*P*theta*c376*e*sin(2*g+f)
    return l375 + l376

# =====================================================================
# Tests
# =====================================================================

test_cases = []
for e_val in [0.01, 0.1, 0.3]:
    for I_deg in [30, 60, 85]:
        for g_deg in [0, 45, 90]:
            for l_val in [0.5, 1.5, 3.0]:
                theta = cos(np.radians(I_deg))
                a = 1.0
                L = sqrt(MU * a)
                eta = sqrt(1 - e_val**2)
                G = L * eta
                H = G * theta
                test_cases.append((L, G, H, l_val, np.radians(g_deg), 0.0,
                                   e_val, I_deg, g_deg, l_val))

print("=" * 100)
print("BH61 Eq.(14) Numerical Verification")
print("=" * 100)
print("S_1* alpha = -1/16 (BH61-required per BH61-D004)")
print()

# ----- delta_p1 and delta_p2: use analytical dS/dl_j + single-layer D -----
for name, bh61_func, inner_func in [
    ("delta_p1", delta_p1_BH61,
     lambda L,G,H,l,g,h: dS_dl_analytical(L,G,H,l,g,h)),
    ("delta_p2", delta_p2_BH61,
     lambda L,G,H,l,g,h: dS_dg_analytical(L,G,H,l,g,h)),
]:
    print(f"\n  {name}  [D(dS/d_j) via analytical inner + FD D-operator]")
    print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>5} | {'BH61':>14} {'Numerical':>14} {'Rel err':>12}")
    print("-" * 72)
    max_rel = 0.0
    n_pass = 0; n_fail = 0
    for tc in test_cases:
        L, G, H, l, g, h, e_val, I_deg, g_deg, l_val = tc
        bh61_val = bh61_func(L, G, H, l, g, h)
        num_val = D_operator(inner_func, L, G, H, l, g, h)
        denom = max(abs(bh61_val), abs(num_val), 1e-15)
        rel = abs(bh61_val - num_val) / denom
        max_rel = max(max_rel, rel)
        if rel >= 1e-3:
            n_fail += 1
            print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.2f} | {bh61_val:14.6e} {num_val:14.6e} {rel:12.2e} ***FAIL")
        else:
            n_pass += 1
    print(f"  {n_pass} pass, {n_fail} fail. Max rel err: {max_rel:.2e}")

# ----- delta_q1, delta_q2, delta_q3: mixed-partial stencil -----
for name, bh61_func, j_idx in [
    ("delta_q1", delta_q1_BH61, 0),
    ("delta_q2", delta_q2_BH61, 1),
    ("delta_q3", delta_q3_BH61, 2),
]:
    print(f"\n  {name}  [D(dS/dL_j) via 2-layer FD, eps_inner=1e-4, eps_outer=1e-5]")
    print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>5} | {'BH61':>14} {'Numerical':>14} {'Rel err':>12}")
    print("-" * 72)
    max_rel = 0.0
    n_pass = 0; n_fail = 0
    for tc in test_cases:
        L, G, H, l, g, h, e_val, I_deg, g_deg, l_val = tc
        try:
            bh61_val = bh61_func(L, G, H, l, g, h)
            num_val = D_of_dS_dLj(j_idx, L, G, H, l, g, h)
            denom = max(abs(bh61_val), abs(num_val), 1e-15)
            rel = abs(bh61_val - num_val) / denom
            max_rel = max(max_rel, rel)
            if rel >= 5e-3:
                n_fail += 1
                print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.2f} | {bh61_val:14.6e} {num_val:14.6e} {rel:12.2e} ***FAIL")
            else:
                n_pass += 1
        except Exception as ex:
            n_fail += 1
            print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.2f} | SKIP: {ex}")
    print(f"  {n_pass} pass, {n_fail} fail. Max rel err: {max_rel:.2e}")
