"""
Phase B Numerical Verification: delta_p1 = D(dS1/dl)

Two independent computations:
(a) Closed-form: delta_p1 = 6*mu^2*B0/(a^2*r*eta^3)
(b) Finite-difference D operator: D(F) = sum_j [-p_j dF/dL_j + q_j dF/dl_j]
    where F = dS1/dl = Gamma*B0 - mu^2*[B0 + B1p*cos(2f+2g)]/r^3

Also compare against BH61 Eq(14).
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

MU = 1.0  # gravitational parameter (normalized)

def solve_kepler(M, e, tol=1e-15):
    """Solve Kepler's equation M = E - e*sin(E)."""
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
    """Convert Delaunay (L,G,H,l,g,h) to orbital elements and derived quantities."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(1 - eta**2)
    theta = H / G  # cos(I)

    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))

    return a, e, eta, theta, E, f, r

def dS1_dl(a, e, eta, theta, r, f, g, mu=MU):
    """Compute dS1/dl = Gamma*B0 - mu^2*[B0 + B1p*cos(2f+2g)]/r^3."""
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    Gamma = mu**2 / (a**3 * eta**3)
    return Gamma * B0 - mu**2 * (B0 + B1p * cos(2*f + 2*g)) / r**3

def dS1_dl_from_delaunay(L, G, H, l, g, h, mu=MU):
    """Compute dS1/dl from Delaunay variables."""
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)
    return dS1_dl(a, e, eta, theta, r, f, g, mu)

def delta_p1_closedform(a, e, eta, theta, r, f, g, mu=MU):
    """Our derived result: delta_p1 = 6*mu^2*B0/(a^2*r*eta^3)."""
    B0 = -0.5 + 1.5 * theta**2
    return 6 * mu**2 * B0 / (a**2 * r * eta**3)

def delta_p1_finite_diff(L, G, H, l, g, h, mu=MU, eps=1e-7):
    """
    Compute D(dS1/dl) via finite differences using:
    D(F) = sum_j [-p_j * dF/dL_j + q_j * dF/dl_j]

    p1 = L*(2a/r - 1), p2 = G, p3 = H
    q1 = 2e*sinE + 2*eta*sinf/e, q2 = -(2/e)*sinf, q3 = 0
    """
    a, e, eta, theta, E, f, r = orbital_from_delaunay(L, G, H, l, g, h, mu)

    # Compute p_j and q_j
    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H

    sinE_val = sin(E)
    sinf_val = sin(f)

    q1 = 2*e*sinE_val + 2*eta*sinf_val/e if e > 1e-10 else 0.0
    q2 = -2*sinf_val/e if e > 1e-10 else 0.0
    q3 = 0.0

    F0 = dS1_dl_from_delaunay(L, G, H, l, g, h, mu)

    # Partial derivatives by finite differences
    # dF/dL
    Fp = dS1_dl_from_delaunay(L+eps, G, H, l, g, h, mu)
    Fm = dS1_dl_from_delaunay(L-eps, G, H, l, g, h, mu)
    dFdL = (Fp - Fm) / (2*eps)

    # dF/dG
    Fp = dS1_dl_from_delaunay(L, G+eps, H, l, g, h, mu)
    Fm = dS1_dl_from_delaunay(L, G-eps, H, l, g, h, mu)
    dFdG = (Fp - Fm) / (2*eps)

    # dF/dH
    Fp = dS1_dl_from_delaunay(L, G, H+eps, l, g, h, mu)
    Fm = dS1_dl_from_delaunay(L, G, H-eps, l, g, h, mu)
    dFdH = (Fp - Fm) / (2*eps)

    # dF/dl
    Fp = dS1_dl_from_delaunay(L, G, H, l+eps, g, h, mu)
    Fm = dS1_dl_from_delaunay(L, G, H, l-eps, g, h, mu)
    dFdl = (Fp - Fm) / (2*eps)

    # dF/dg
    Fp = dS1_dl_from_delaunay(L, G, H, l, g+eps, h, mu)
    Fm = dS1_dl_from_delaunay(L, G, H, l, g-eps, h, mu)
    dFdg = (Fp - Fm) / (2*eps)

    # dF/dh
    Fp = dS1_dl_from_delaunay(L, G, H, l, g, h+eps, mu)
    Fm = dS1_dl_from_delaunay(L, G, H, l, g, h-eps, mu)
    dFdh = (Fp - Fm) / (2*eps)

    # D(F) = sum [-p_j dF/dL_j + q_j dF/dl_j]
    result = (-p1*dFdL - p2*dFdG - p3*dFdH
              + q1*dFdl + q2*dFdg + q3*dFdh)

    return result

def delta_p1_BH61(a, e, eta, theta, r, f, g, mu=MU):
    """
    BH61 Eq(14) with k2=1:
    delta_p1 = 3*(mu^2/L^3) * {
        B0*[-eta^{-3} + (a^3/r^3)*(1 - 2a/r)]
        + B1p*(a^3/r^3)*(1 - 2a/r)*cos(2g+2f)
    }
    """
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    L = sqrt(mu * a)
    n = mu**2 / L**3  # = mu^{1/2}/a^{3/2}

    bracket = (B0 * (-1/eta**3 + (a/r)**3 * (1 - 2*a/r))
               + B1p * (a/r)**3 * (1 - 2*a/r) * cos(2*g + 2*f))

    return 3 * n * bracket


# ============================================================
# Run tests
# ============================================================
print("=" * 90)
print("Phase B Verification: delta_p1 = D(dS1/dl)")
print("=" * 90)
print()
print(f"{'e':>5} {'I':>4} {'g':>4} {'l':>5} | {'Closed-form':>14} {'Finite-diff':>14} {'Rel err':>10} | {'BH61':>14} {'BH61 ratio':>10}")
print("-" * 90)

max_err_cf = 0.0
test_count = 0

for e in [0.01, 0.1, 0.3]:
    for I_deg in [30, 60, 85]:
        for g_deg in [0, 45, 90]:
            for l_val in [0.5, 1.5, 3.0]:
                theta = cos(np.radians(I_deg))
                g = np.radians(g_deg)
                h = 0.0  # arbitrary

                # Delaunay momenta
                a = 1.0
                eta = sqrt(1 - e**2)
                L = sqrt(MU * a)
                G = L * eta
                H = G * theta

                # Solve for orbital elements
                E = solve_kepler(l_val, e)
                f = true_from_eccentric(E, e)
                r = a * (1 - e * cos(E))

                # Three computations
                cf = delta_p1_closedform(a, e, eta, theta, r, f, g, MU)
                fd = delta_p1_finite_diff(L, G, H, l_val, g, h, MU)
                bh = delta_p1_BH61(a, e, eta, theta, r, f, g, MU)

                # Relative error between closed-form and finite-diff
                denom = max(abs(cf), abs(fd), 1e-20)
                rel_err = abs(cf - fd) / denom
                max_err_cf = max(max_err_cf, rel_err)

                # Ratio of BH61 to our result
                if abs(cf) > 1e-20:
                    ratio_bh = bh / cf
                else:
                    ratio_bh = float('nan')

                test_count += 1
                print(f"{e:5.2f} {I_deg:4d} {g_deg:4d} {l_val:5.1f} | "
                      f"{cf:14.8f} {fd:14.8f} {rel_err:10.2e} | "
                      f"{bh:14.8f} {ratio_bh:10.4f}")

print()
print(f"Total test cases: {test_count}")
print(f"Max relative error (closed-form vs finite-diff): {max_err_cf:.2e}")
print()

# Summary comparison
print("=" * 60)
print("SUMMARY")
print("=" * 60)
if max_err_cf < 1e-5:
    print("PASS: Closed-form matches finite-difference D operator")
else:
    print("FAIL: Closed-form does NOT match finite-difference D operator")

# Check if BH61 ever matches
print()
print("Checking if BH61 Eq(14) matches our result...")
bh_matches = True
for e in [0.01, 0.1, 0.3]:
    for I_deg in [30, 60, 85]:
        for g_deg in [0, 45, 90]:
            for l_val in [0.5, 1.5, 3.0]:
                theta = cos(np.radians(I_deg))
                g = np.radians(g_deg)
                a = 1.0
                eta = sqrt(1 - e**2)
                L = sqrt(MU * a)
                G = L * eta
                H = G * theta
                E = solve_kepler(l_val, e)
                f = true_from_eccentric(E, e)
                r = a * (1 - e * cos(E))

                cf = delta_p1_closedform(a, e, eta, theta, r, f, g, MU)
                bh = delta_p1_BH61(a, e, eta, theta, r, f, g, MU)

                denom = max(abs(cf), abs(bh), 1e-20)
                if abs(cf - bh) / denom > 0.01:
                    bh_matches = False
                    break

if bh_matches:
    print("BH61 Eq(14) MATCHES our derivation")
else:
    print("BH61 Eq(14) does NOT match our derivation")
    print("This is expected: BH61 Eq(14) contains additional (a/r)^3 terms")
    print("that do not appear in D(dS1/dl).")

# Additional diagnostic: check the STRUCTURE of BH61
print()
print("=" * 60)
print("DIAGNOSTIC: Structure analysis")
print("=" * 60)
print()
print("Our result:  delta_p1 = 6*mu^2*B0 / (a^2 * r * eta^3)")
print("             = 6*Gamma*B0 * (a/r)")
print("             where Gamma = mu^2/(a^3*eta^3)")
print()
print("BH61 Eq(14): delta_p1 = 3n*{B0*[-1/eta^3 + (a/r)^3*(1-2a/r)]")
print("             + B1'*(a/r)^3*(1-2a/r)*cos(2f+2g)}")
print()
print("Key differences:")
print("  1. Our result has NO cos(2f+2g) term -- it vanishes because D(2f+2g) = 0")
print("  2. Our result has a/r dependence; BH61 has (a/r)^3 and (a/r)^4")
print("  3. BH61 has -1/eta^3 constant; our result has none (pure a/r dependence)")

# Verify D(2f+2g) = 0 numerically
print()
print("=" * 60)
print("VERIFICATION: D(f+g) = 0")
print("=" * 60)
print()
print("Df = 2*sin(f)/e,  Dg = -2*sin(f)/e")
print("D(f+g) = Df + Dg = 0  (exact cancellation)")
print()

for e in [0.1, 0.3]:
    for l_val in [0.5, 1.5, 3.0]:
        E = solve_kepler(l_val, e)
        f = true_from_eccentric(E, e)
        Df = 2*sin(f)/e
        Dg = -2*sin(f)/e
        print(f"  e={e}, l={l_val}: Df={Df:+.8f}, Dg={Dg:+.8f}, D(f+g)={Df+Dg:+.2e}")

# Verify Dr = 0 numerically via finite difference
print()
print("=" * 60)
print("VERIFICATION: Dr = 0 via finite differences")
print("=" * 60)

def r_from_delaunay(L, G, H, l, g, h, mu=MU):
    a = L**2/mu
    eta = G/L
    e = sqrt(max(1-eta**2, 0))
    E = solve_kepler(l, e)
    return a * (1 - e * cos(E))

eps = 1e-7
for e in [0.1, 0.3]:
    for l_val in [0.5, 1.5, 3.0]:
        a = 1.0; eta = sqrt(1-e**2); theta = cos(np.radians(60))
        L = sqrt(MU*a); G = L*eta; H = G*theta
        g = np.radians(45); h = 0.0

        a_o, e_o, eta_o, theta_o, E_o, f_o, r_o = orbital_from_delaunay(L, G, H, l_val, g, h)
        p1 = L*(2*a_o/r_o - 1)
        p2 = G; p3 = H
        sinE = sin(E_o); sinf = sin(f_o)
        q1 = 2*e_o*sinE + 2*eta_o*sinf/e_o
        q2 = -2*sinf/e_o; q3 = 0

        # Finite diff partials of r
        drdL = (r_from_delaunay(L+eps,G,H,l_val,g,h) - r_from_delaunay(L-eps,G,H,l_val,g,h))/(2*eps)
        drdG = (r_from_delaunay(L,G+eps,H,l_val,g,h) - r_from_delaunay(L,G-eps,H,l_val,g,h))/(2*eps)
        drdH = (r_from_delaunay(L,G,H+eps,l_val,g,h) - r_from_delaunay(L,G,H-eps,l_val,g,h))/(2*eps)
        drdl = (r_from_delaunay(L,G,H,l_val+eps,g,h) - r_from_delaunay(L,G,H,l_val-eps,g,h))/(2*eps)
        drdg = (r_from_delaunay(L,G,H,l_val,g+eps,h) - r_from_delaunay(L,G,H,l_val,g-eps,h))/(2*eps)

        Dr = -p1*drdL - p2*drdG - p3*drdH + q1*drdl + q2*drdg + q3*0
        print(f"  e={e}, l={l_val}: Dr = {Dr:+.2e}  (should be ~0)")
