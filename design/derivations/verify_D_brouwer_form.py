"""
Step 2: Decisive Test — D of Brouwer's form by finite differences.

Computes D(Brouwer's dS1/dl) numerically and compares against:
  (a) -6*Gamma*B0*(a/r) — predicted by linearity if forms are negatives
  (b) BH61 Eq(14) — what the breakthrough claims
  (c) D(our dS1/dl) numerical — for cross-check

Tests with both mu=1,a=1 and mu=1,a=2 to detect a=mu coincidences.
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


# === The two forms of dS1/dl ===

def dS1_dl_ours(L, G, H, l, g, h, mu):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    Gamma = mu**2 / (a**3 * eta**3)
    return Gamma * B0 - mu**2 * (B0 + B1p * cos(2*f + 2*g)) / r**3

def dS1_dl_brouwer(L, G, H, l, g, h, mu):
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    rho = a / r
    n = mu**2 / L**3
    return n * (B0 * (rho**3 - 1.0/eta**3) + B1p * rho**3 * cos(2*f + 2*g))


# === Closed-form candidates ===

def candidate_a(L, G, H, l, g, h, mu):
    """Predicted by linearity if forms are exact negatives: -6*Gamma*B0*(a/r)."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    Gamma = mu**2 / (a**3 * eta**3)
    return -6 * Gamma * B0 * (a / r)

def candidate_bh61(L, G, H, l, g, h, mu):
    """BH61 Eq(14) with k2=1."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    n = mu**2 / L**3
    rho = a / r
    bracket = (B0 * (-1.0/eta**3 + rho**3 * (1 - 2*rho))
               + B1p * rho**3 * (1 - 2*rho) * cos(2*g + 2*f))
    return 3 * n * bracket


# === Finite-difference D operator ===

def D_finite_diff(func, L, G, H, l, g, h, mu, eps=1e-7):
    """D(F) = sum_j [-p_j * dF/dL_j + q_j * dF/dl_j]."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)

    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H
    sinE = sin(E)
    sinf = sin(f)
    q1 = 2*e*sinE + 2*eta*sinf/e if e > 1e-10 else 0.0
    q2 = -2*sinf/e if e > 1e-10 else 0.0
    q3 = 0.0

    dFdL = (func(L+eps, G, H, l, g, h, mu) - func(L-eps, G, H, l, g, h, mu)) / (2*eps)
    dFdG = (func(L, G+eps, H, l, g, h, mu) - func(L, G-eps, H, l, g, h, mu)) / (2*eps)
    dFdH = (func(L, G, H+eps, l, g, h, mu) - func(L, G, H-eps, l, g, h, mu)) / (2*eps)
    dFdl = (func(L, G, H, l+eps, g, h, mu) - func(L, G, H, l-eps, g, h, mu)) / (2*eps)
    dFdg = (func(L, G, H, l, g+eps, h, mu) - func(L, G, H, l, g-eps, h, mu)) / (2*eps)
    dFdh = (func(L, G, H, l, g, h+eps, mu) - func(L, G, H, l, g, h-eps, mu)) / (2*eps)

    return (-p1*dFdL - p2*dFdG - p3*dFdH + q1*dFdl + q2*dFdg + q3*dFdh)


def run_decisive_test(mu, a_val, label):
    print(f"\n{'='*100}")
    print(f"DECISIVE TEST: D(Brouwer's dS1/dl) — {label}")
    print(f"{'='*100}")
    print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>5} | {'D(Brouwer) FD':>16} {'Cand(a) -6GB0r':>16} {'err_a':>10} {'BH61 Eq14':>16} {'err_bh':>10} {'D(ours) FD':>16}")
    print("-" * 120)

    max_err_a = 0.0
    max_err_bh = 0.0
    n_tests = 0

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
                    # D of Brouwer's form
                    D_brow = D_finite_diff(dS1_dl_brouwer, L, G, H, l_val, g_val, h_val, mu)

                    # D of our form (cross-check)
                    D_ours = D_finite_diff(dS1_dl_ours, L, G, H, l_val, g_val, h_val, mu)

                    # Candidates
                    ca = candidate_a(L, G, H, l_val, g_val, h_val, mu)
                    cbh = candidate_bh61(L, G, H, l_val, g_val, h_val, mu)

                    denom = max(abs(D_brow), 1e-30)
                    err_a = abs(D_brow - ca) / denom
                    err_bh = abs(D_brow - cbh) / denom

                    max_err_a = max(max_err_a, err_a)
                    max_err_bh = max(max_err_bh, err_bh)
                    n_tests += 1

                    print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.1f} | "
                          f"{D_brow:16.8e} {ca:16.8e} {err_a:10.2e} "
                          f"{cbh:16.8e} {err_bh:10.2e} {D_ours:16.8e}")

    print(f"\n{n_tests} test cases completed.")
    print(f"Max rel error vs candidate (a) [-6*Gamma*B0*(a/r)]: {max_err_a:.2e}")
    print(f"Max rel error vs candidate (b) [BH61 Eq(14)]:      {max_err_bh:.2e}")

    if max_err_a < 1e-4:
        print(">>> CANDIDATE (a) MATCHES: D(Brouwer) = -6*Gamma*B0*(a/r)")
        print("    => Linearity holds, breakthrough algebra is WRONG")
    elif max_err_bh < 1e-4:
        print(">>> CANDIDATE (b) MATCHES: D(Brouwer) = BH61 Eq(14)")
        print("    => Breakthrough is correct, forms are NOT exact negatives")
    else:
        print(">>> NEITHER candidate matches! Something else is going on.")

    # Also check: does D(ours) = -D(Brouwer)?
    print(f"\nCross-check: D(ours) vs -D(Brouwer)")
    max_cross = 0.0
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
                    D_brow = D_finite_diff(dS1_dl_brouwer, L, G, H, l_val, g_val, h_val, mu)
                    D_ours = D_finite_diff(dS1_dl_ours, L, G, H, l_val, g_val, h_val, mu)
                    denom = max(abs(D_ours), abs(D_brow), 1e-30)
                    cross_err = abs(D_ours + D_brow) / denom
                    max_cross = max(max_cross, cross_err)

    print(f"  Max |D(ours) + D(Brouwer)| / max(|D(ours)|, |D(Brouwer)|) = {max_cross:.2e}")
    if max_cross < 1e-4:
        print("  => D(ours) = -D(Brouwer) CONFIRMED (linearity holds)")
    else:
        print("  => D(ours) != -D(Brouwer) — the two forms are NOT exact negatives")


# ============================================================
# Run with multiple parameter sets
# ============================================================
run_decisive_test(mu=1.0, a_val=1.0, label="Normalized (mu=1, a=1)")
run_decisive_test(mu=1.0, a_val=2.0, label="mu=1, a=2 (breaks a=mu coincidence)")
run_decisive_test(mu=398600.4418, a_val=7000.0, label="Physical units")
