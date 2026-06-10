"""
Compute dp2 = D(dS1/dg) and verify numerically.

dS1/dg = -Gamma*B1p*[cos(2f+2g) + e*cos(f+2g) + (e/3)*cos(3f+2g)]

D acts on this. Key: D(2f+2g)=0 but D(f+2g)=-2sinf/e, D(3f+2g)=+2sinf/e.
So dp2 SHOULD have nontrivial harmonic structure.
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

MU = 1.0

def solve_kepler(M, e, tol=1e-15):
    E = M + e * sin(M)
    for _ in range(100):
        dE = (M - E + e * sin(E)) / (1 - e * cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def true_from_eccentric(E, e):
    return 2 * arctan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))

def get_orbital(L, G, H, l, g, h, mu=MU):
    a = L**2 / mu
    eta = G / L
    e = sqrt(1 - eta**2)
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    return a, e, eta, theta, f, r, E

def D_numerical(func, L, G, H, l, g, h, eps=1e-7, mu=MU):
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)
    p1 = L * (2*a/r - 1); p2 = G; p3 = H
    q1 = 2*e*sin(E) + 2*eta*sin(f)/e if e > 1e-14 else 0.0
    q2 = -2*sin(f)/e if e > 1e-14 else 0.0
    D_acts = [-p1, -p2, -p3, q1, q2, 0.0]
    base = [L, G, H, l, g, h]
    result = 0.0
    for i in range(6):
        if abs(D_acts[i]) < 1e-20: continue
        bp = list(base); bp[i] += eps
        bm = list(base); bm[i] -= eps
        result += D_acts[i] * (func(*bp) - func(*bm)) / (2*eps)
    return result

def dS1_dg(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B1p = 1.5 * (1 - theta**2)
    return -Gamma * B1p * (cos(2*f+2*g) + e*cos(f+2*g) + (e/3.0)*cos(3*f+2*g))

# Closed form derivation:
# D(dS1/dg) = D[-Gamma*B1p*(cos(2f+2g) + e*cos(f+2g) + (e/3)*cos(3f+2g))]
#
# = -D(Gamma)*B1p*[...] - Gamma*B1p*D[cos(2f+2g) + e*cos(f+2g) + (e/3)*cos(3f+2g)]
#
# D(Gamma) = 6*Gamma*rho
# D[cos(2f+2g)] = 0  (since D(2f+2g) = 0)
# D[e*cos(f+2g)] = De*cos(f+2g) + e*D[cos(f+2g)]
#   De = -2(e+cosf)
#   D[cos(f+2g)] = -sin(f+2g)*D(f+2g) = -sin(f+2g)*(-2sinf/e) = 2sinf*sin(f+2g)/e
#   So D[e*cos(f+2g)] = -2(e+cosf)*cos(f+2g) + 2sinf*sin(f+2g)
#                      = -2e*cos(f+2g) - 2cosf*cos(f+2g) + 2sinf*sin(f+2g)
#                      = -2e*cos(f+2g) - 2cos(2f+2g)   [using cos(A)cos(B)-sin(A)sin(B)=cos(A+B)]
#
# D[(e/3)*cos(3f+2g)] = (De/3)*cos(3f+2g) + (e/3)*D[cos(3f+2g)]
#   D[cos(3f+2g)] = -sin(3f+2g)*D(3f+2g) = -sin(3f+2g)*(2sinf/e) = -2sinf*sin(3f+2g)/e
#   So D[(e/3)*cos(3f+2g)] = -2(e+cosf)*cos(3f+2g)/3 - 2sinf*sin(3f+2g)/3
#                           = -(2e/3)*cos(3f+2g) - (2/3)[cosf*cos(3f+2g) + sinf*sin(3f+2g)]
#                           = -(2e/3)*cos(3f+2g) - (2/3)*cos(2f+2g)
#
# Total D[e*cos(f+2g) + (e/3)*cos(3f+2g)]:
# = -2e*cos(f+2g) - 2cos(2f+2g) - (2e/3)*cos(3f+2g) - (2/3)*cos(2f+2g)
# = -2e*cos(f+2g) - (8/3)*cos(2f+2g) - (2e/3)*cos(3f+2g)
#
# So D(dS1/dg) = -6*Gamma*rho*B1p*[cos(2f+2g) + e*cos(f+2g) + (e/3)*cos(3f+2g)]
#                -Gamma*B1p*[-2e*cos(f+2g) - (8/3)*cos(2f+2g) - (2e/3)*cos(3f+2g)]
# = -Gamma*B1p*[6rho*cos(2f+2g) + 6e*rho*cos(f+2g) + (6e*rho/3)*cos(3f+2g)
#               - 2e*cos(f+2g) - (8/3)*cos(2f+2g) - (2e/3)*cos(3f+2g)]
# = -Gamma*B1p*[(6rho - 8/3)*cos(2f+2g) + (6e*rho - 2e)*cos(f+2g) + (2e*rho - 2e/3)*cos(3f+2g)]
# = -Gamma*B1p*[(6rho - 8/3)*cos(2f+2g) + 2e*(3rho-1)*cos(f+2g) + (2e/3)*(3rho-1)*cos(3f+2g)]
#
# Factor: 3rho-1 = 3a/r - 1 = (3+3ecosf-eta^2)/eta^2 = (3+3ecosf-1+e^2)/eta^2 = (2+e^2+3ecosf)/eta^2

def dp2_closed(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B1p = 1.5 * (1 - theta**2)
    rho = a / r
    return -Gamma * B1p * (
        (6*rho - 8.0/3.0) * cos(2*f+2*g)
        + 2*e*(3*rho - 1) * cos(f+2*g)
        + (2*e/3.0)*(3*rho - 1) * cos(3*f+2*g)
    )

# Verify
a_val = 1.0
L_val = sqrt(MU * a_val)

print("=" * 90)
print("VERIFICATION: dp2 = D(dS1/dg)")
print("=" * 90)
print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'D(dS1/dg) FD':>18} {'Closed form':>18} {'RelErr':>12} {'Status'}")
print("-" * 90)

max_err = 0.0
pass_count = 0
fail_count = 0

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L_val * eta_val
    for I_deg in [30, 60, 85]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for g_deg in [0, 45, 90]:
            gv = np.radians(g_deg)
            for l_val in [0.5, 1.5, 3.0]:
                D_num = D_numerical(dS1_dg, L_val, G, H, l_val, gv, 0.0)
                cf = dp2_closed(L_val, G, H, l_val, gv, 0.0)
                denom = max(abs(D_num), abs(cf), 1e-15)
                rel_err = abs(D_num - cf) / denom
                max_err = max(max_err, rel_err)
                status = "PASS" if rel_err < 1e-4 else "FAIL"
                if status == "PASS": pass_count += 1
                else: fail_count += 1
                print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | "
                      f"{D_num:18.10e} {cf:18.10e} {rel_err:12.2e} {status}")

print("=" * 90)
print(f"dp2 verification: {pass_count} PASS, {fail_count} FAIL, max err = {max_err:.2e}")
print()
print("DERIVED FORMULA:")
print("dp2 = D(dS1/dg) = -Gamma*B1p*[(6*rho - 8/3)*cos(2f+2g)")
print("                                + 2e*(3*rho-1)*cos(f+2g)")
print("                                + (2e/3)*(3*rho-1)*cos(3f+2g)]")
print()
print("This has the SAME harmonic structure as BH61 Eq(14) for dp2:")
print("  cos(2f+2g), cos(f+2g), cos(3f+2g) terms with rho-dependent coefficients")
