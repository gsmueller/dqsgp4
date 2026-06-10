"""
Step 3: Check the breakthrough algebra term by term.

Verifies each D-action claimed in phase_b_breakthrough.md:
  D(L^{-3}) = 3*(2*rho-1)/L^3
  D(rho^3) = -6*rho^3*(2*rho-1)
  D(eta^{-3}) = -6*eta^{-3}*(rho-1)

Then assembles these via Leibniz rule and compares to BH61 Eq(14).
Also checks the key identity: D(mu^2/L^3) should equal D(Gamma*eta^3).

Runs at a single detailed test point AND the full 81-point grid.
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

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

def orbital_from_delaunay(L, G, H, l, g, h, mu):
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    return a, e, eta, theta, E, f, r

def D_of_scalar(func, L, G, H, l, g, h, mu, eps=1e-7):
    """Compute D(func) via finite differences using D = sum[-p_j dF/dL_j + q_j dF/dl_j]."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H
    sinE = sin(E); sinf = sin(f)
    q1 = 2*e*sinE + 2*eta*sinf/e if e > 1e-10 else 0.0
    q2 = -2*sinf/e if e > 1e-10 else 0.0
    q3 = 0.0

    dFdL = (func(L+eps,G,H,l,g,h,mu) - func(L-eps,G,H,l,g,h,mu)) / (2*eps)
    dFdG = (func(L,G+eps,H,l,g,h,mu) - func(L,G-eps,H,l,g,h,mu)) / (2*eps)
    dFdH = (func(L,G,H+eps,l,g,h,mu) - func(L,G,H-eps,l,g,h,mu)) / (2*eps)
    dFdl = (func(L,G,H,l+eps,g,h,mu) - func(L,G,H,l-eps,g,h,mu)) / (2*eps)
    dFdg = (func(L,G,H,l,g+eps,h,mu) - func(L,G,H,l,g-eps,h,mu)) / (2*eps)
    dFdh = (func(L,G,H,l,g,h+eps,mu) - func(L,G,H,l,g,h-eps,mu)) / (2*eps)

    return -p1*dFdL - p2*dFdG - p3*dFdH + q1*dFdl + q2*dFdg + q3*dFdh


# === Individual scalar functions to test D-actions on ===

def func_L_inv3(L, G, H, l, g, h, mu):
    return 1.0 / L**3

def func_rho3(L, G, H, l, g, h, mu):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    return (a / r)**3

def func_eta_inv3(L, G, H, l, g, h, mu):
    eta = G / L
    return 1.0 / eta**3

def func_Gamma(L, G, H, l, g, h, mu):
    a = L**2 / mu
    eta = G / L
    return mu**2 / (a**3 * eta**3)

def func_n(L, G, H, l, g, h, mu):
    """n = mu^2/L^3 = Gamma * eta^3"""
    return mu**2 / L**3

def func_Gamma_eta3(L, G, H, l, g, h, mu):
    """Gamma * eta^3 — should equal n = mu^2/L^3"""
    a = L**2 / mu
    eta = G / L
    Gamma = mu**2 / (a**3 * eta**3)
    return Gamma * eta**3


# === Claimed closed-form D-actions ===

def claimed_D_L_inv3(L, G, H, l, g, h, mu):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    rho = a / r
    return 3 * (2*rho - 1) / L**3

def claimed_D_rho3(L, G, H, l, g, h, mu):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    rho = a / r
    return -6 * rho**3 * (2*rho - 1)

def claimed_D_eta_inv3(L, G, H, l, g, h, mu):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    rho = a / r
    return -6 / eta**3 * (rho - 1)

def claimed_D_Gamma(L, G, H, l, g, h, mu):
    """D(Gamma) = 6*Gamma*(a/r) — from cleanroom_phase_b_results."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    return 6 * Gamma * (a / r)


def check_D_action(label, func, claimed, L, G, H, l, g, h, mu):
    """Compare numerical D(func) against claimed closed form."""
    numerical = D_of_scalar(func, L, G, H, l, g, h, mu)
    analytical = claimed(L, G, H, l, g, h, mu)
    denom = max(abs(numerical), abs(analytical), 1e-30)
    err = abs(numerical - analytical) / denom
    status = "PASS" if err < 1e-4 else "FAIL"
    print(f"  {label:30s}: numerical={numerical:+16.8e}  claimed={analytical:+16.8e}  err={err:.2e}  {status}")
    return err


def run_detailed_test(mu, a_val, e_val, I_deg, g_deg, l_val):
    eta_val = sqrt(1 - e_val**2)
    L = sqrt(mu * a_val)
    G = L * eta_val
    theta_val = cos(np.radians(I_deg))
    H = G * theta_val
    g_val = np.radians(g_deg)
    h_val = 0.0

    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l_val, g_val, h_val, mu)
    rho = a / r
    n = mu**2 / L**3
    Gamma = mu**2 / (a**3 * eta**3)

    print(f"\n{'='*80}")
    print(f"DETAILED TEST: mu={mu}, a={a_val}, e={e_val}, I={I_deg}, g={g_deg}, l={l_val}")
    print(f"  Derived: rho={rho:.8f}, eta={eta:.8f}, n={n:.8e}, Gamma={Gamma:.8e}")
    print(f"  Gamma*eta^3 = {Gamma*eta**3:.8e},  n = {n:.8e},  equal? {abs(Gamma*eta**3 - n) < 1e-12*abs(n)}")
    print(f"{'='*80}")

    # Check individual D-actions
    print("\n--- Individual D-actions ---")
    check_D_action("D(L^{-3})", func_L_inv3, claimed_D_L_inv3, L, G, H, l_val, g_val, h_val, mu)
    check_D_action("D(rho^3)", func_rho3, claimed_D_rho3, L, G, H, l_val, g_val, h_val, mu)
    check_D_action("D(eta^{-3})", func_eta_inv3, claimed_D_eta_inv3, L, G, H, l_val, g_val, h_val, mu)
    check_D_action("D(Gamma)", func_Gamma, claimed_D_Gamma, L, G, H, l_val, g_val, h_val, mu)

    # THE KEY CHECK: D(n) = D(mu^2/L^3) vs D(Gamma*eta^3)
    print("\n--- Key identity: D(n) = D(Gamma*eta^3) ---")
    D_n = D_of_scalar(func_n, L, G, H, l_val, g_val, h_val, mu)
    D_Ge3 = D_of_scalar(func_Gamma_eta3, L, G, H, l_val, g_val, h_val, mu)
    print(f"  D(n = mu^2/L^3) numerical:       {D_n:+16.8e}")
    print(f"  D(Gamma*eta^3) numerical:         {D_Ge3:+16.8e}")
    print(f"  Difference:                        {abs(D_n - D_Ge3):.2e}")

    # Claimed D(L^{-3}) route: D(mu^2/L^3) = mu^2 * 3*(2*rho-1)/L^3
    D_n_via_L = mu**2 * 3 * (2*rho - 1) / L**3
    print(f"  D(n) via D(L^{{-3}}):               {D_n_via_L:+16.8e}")
    print(f"  err vs numerical:                  {abs(D_n - D_n_via_L)/max(abs(D_n),1e-30):.2e}")

    # D(Gamma*eta^3) via product rule: D(Gamma)*eta^3 + Gamma*D(eta^3)
    D_Gamma_val = 6 * Gamma * rho
    D_eta3_val = 3 * eta**2 * 2 * eta * (rho - 1)  # = 6*eta^3*(rho-1)
    D_Ge3_product = D_Gamma_val * eta**3 + Gamma * D_eta3_val
    print(f"  D(Gamma*eta^3) via product rule:   {D_Ge3_product:+16.8e}")
    print(f"    D(Gamma)*eta^3 = {D_Gamma_val * eta**3:+16.8e}")
    print(f"    Gamma*D(eta^3) = {Gamma * D_eta3_val:+16.8e}")
    print(f"    Sum:             {D_Ge3_product:+16.8e}")
    print(f"  D(Gamma*eta^3) = 6*Gamma*eta^3*(2*rho-1)? {6*Gamma*eta**3*(2*rho-1):+16.8e}")
    print(f"  D(n) = 3*n*(2*rho-1)?                     {3*n*(2*rho-1):+16.8e}")
    print()
    print(f"  3*n*(2*rho-1) = {3*n*(2*rho-1):+16.8e}")
    print(f"  6*n*(2*rho-1) = {6*n*(2*rho-1):+16.8e}")
    print(f"  D(n) numerical = {D_n:+16.8e}")
    print()

    if abs(D_n - 3*n*(2*rho-1))/max(abs(D_n),1e-30) < 1e-4:
        print("  >>> D(n) = 3*n*(2*rho-1)  [consistent with D(L^{-3}) = 3*(2*rho-1)/L^3]")
    elif abs(D_n - 6*n*(2*rho-1))/max(abs(D_n),1e-30) < 1e-4:
        print("  >>> D(n) = 6*n*(2*rho-1)  [consistent with D(Gamma)*eta^3 + Gamma*D(eta^3)]")
    else:
        print("  >>> D(n) matches NEITHER 3*n nor 6*n formula!")

    # Assemble the breakthrough Leibniz result
    print("\n--- Assembling breakthrough Leibniz expansion ---")
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    sigma1 = rho**3 - 1.0/eta**3
    sigma2 = rho**3 * cos(2*f + 2*g_val)

    # Brouwer's dS1/dl = n * (B0*sigma1 + B1p*sigma2)
    # D(n*(B0*s1+B1p*s2)) = D(n)*(B0*s1+B1p*s2) + n*D(B0*s1+B1p*s2)
    # D(B0) = D(B1p) = 0 (Dtheta=0)
    # D(sigma1) = D(rho^3) - D(eta^{-3})
    # D(sigma2) = D(rho^3)*cos(2f+2g) + rho^3*D(cos(2f+2g))
    # D(cos(2f+2g)) = 0 since D(f+g)=0

    D_sigma1 = -6*rho**3*(2*rho-1) - (-6/eta**3*(rho-1))
    D_sigma2 = -6*rho**3*(2*rho-1) * cos(2*f + 2*g_val) + 0  # D(cos)=0

    D_n_val = 3*n*(2*rho-1)  # Using the D(L^{-3}) route

    assembled = D_n_val*(B0*sigma1 + B1p*sigma2) + n*(B0*D_sigma1 + B1p*D_sigma2)
    print(f"  Assembled from claimed D-actions: {assembled:+16.8e}")

    # BH61 Eq(14)
    bh61 = 3*n*(B0*(-1/eta**3 + rho**3*(1-2*rho)) + B1p*rho**3*(1-2*rho)*cos(2*g_val+2*f))
    print(f"  BH61 Eq(14):                      {bh61:+16.8e}")
    print(f"  Assembled vs BH61 err:             {abs(assembled-bh61)/max(abs(bh61),1e-30):.2e}")

    # Numerical D(Brouwer's dS1/dl)
    def dS1_dl_brouwer_local(L, G, H, l, g, h, mu):
        a2, e2, eta2, theta2, E2, f2, r2 = orbital_from_delaunay(L, G, H, l, g, h, mu)
        B0l = -0.5 + 1.5*theta2**2
        B1pl = 1.5*(1-theta2**2)
        rho2 = a2/r2
        nl = mu**2/L**3
        return nl*(B0l*(rho2**3 - 1/eta2**3) + B1pl*rho2**3*cos(2*f2+2*g))

    D_brow_num = D_of_scalar(dS1_dl_brouwer_local, L, G, H, l_val, g_val, h_val, mu)
    print(f"  D(Brouwer form) numerical:        {D_brow_num:+16.8e}")
    print(f"  err vs assembled:                  {abs(D_brow_num-assembled)/max(abs(D_brow_num),1e-30):.2e}")
    print(f"  err vs BH61:                       {abs(D_brow_num-bh61)/max(abs(D_brow_num),1e-30):.2e}")

    # Our result: -6*Gamma*B0*(a/r)
    neg_ours = -6*Gamma*B0*rho
    print(f"  -6*Gamma*B0*(a/r):                {neg_ours:+16.8e}")
    print(f"  err vs D(Brouwer) numerical:       {abs(D_brow_num-neg_ours)/max(abs(D_brow_num),1e-30):.2e}")


# ============================================================
# Run tests
# ============================================================
run_detailed_test(mu=1.0, a_val=1.0, e_val=0.1, I_deg=60, g_deg=45, l_val=1.5)
run_detailed_test(mu=1.0, a_val=2.0, e_val=0.1, I_deg=60, g_deg=45, l_val=1.5)
run_detailed_test(mu=1.0, a_val=1.0, e_val=0.3, I_deg=85, g_deg=90, l_val=0.5)

# Quick grid test for D(n) identity
print(f"\n{'='*80}")
print("GRID TEST: D(n) = D(mu^2/L^3) — is it 3*n*(2*rho-1) or 6*n*(2*rho-1)?")
print(f"{'='*80}")
print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>5} | {'D(n) FD':>16} {'3*n*(2r-1)':>16} {'err_3':>10} {'6*n*(2r-1)':>16} {'err_6':>10}")
print("-" * 90)

mu = 1.0; a_val = 1.0
for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    L = sqrt(mu * a_val)
    G = L * eta_val
    for I_deg in [30, 60, 85]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for g_deg in [0, 45, 90]:
            g_val = np.radians(g_deg)
            for l_val in [0.5, 1.5, 3.0]:
                a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l_val, g_val, 0.0, mu)
                rho = a / r
                n = mu**2 / L**3
                Dn = D_of_scalar(func_n, L, G, H, l_val, g_val, 0.0, mu)
                c3 = 3*n*(2*rho-1)
                c6 = 6*n*(2*rho-1)
                d = max(abs(Dn), 1e-30)
                print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.1f} | "
                      f"{Dn:16.8e} {c3:16.8e} {abs(Dn-c3)/d:10.2e} "
                      f"{c6:16.8e} {abs(Dn-c6)/d:10.2e}")
