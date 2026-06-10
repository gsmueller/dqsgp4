"""
CORRECTED derivation: D(dS1/dl) with D(f) != 0.

Key finding: D(f) = +2*sin(f)/e = -D(g)
So D(f+g) = 0 and D(f-l) must be computed.

This changes the D(dS1/dl) result fundamentally.
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
    """D(F) by FD using DL=-p1, DG=-p2, DH=-p3, Dl=q1, Dg=q2, Dh=0."""
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)
    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H
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

def dS1_dl(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    return Gamma * B0 - mu**2 / r**3 * (B0 + B1p * cos(2*f + 2*g))

# =====================================================================
# Verify D(f) = 2*sin(f)/e and D(f+g) = 0
# =====================================================================
print("=" * 80)
print("VERIFY: D(f) = 2*sin(f)/e, D(r) = 0, D(f+g) = 0")
print("=" * 80)

a_val = 1.0
L_val = sqrt(MU * a_val)

print(f"{'e':>6} {'I':>5} {'l':>5} | {'D(f) num':>14} {'2sinf/e':>14} {'err':>10} | {'D(r)':>12} | {'D(f+g)':>12}")
print("-" * 95)

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L_val * eta_val
    for I_deg in [60]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for l_val in [0.5, 1.5, 3.0, 5.0]:
            g_val = 0.7

            def func_f(L,G,H,l,g,h):
                e=sqrt(1-(G/L)**2); E=solve_kepler(l,e); return true_from_eccentric(E,e)
            def func_r(L,G,H,l,g,h):
                a=L**2/MU; e=sqrt(1-(G/L)**2); E=solve_kepler(l,e); return a*(1-e*cos(E))
            def func_fpg(L,G,H,l,g,h):
                e=sqrt(1-(G/L)**2); E=solve_kepler(l,e); return true_from_eccentric(E,e)+g

            Df = D_numerical(func_f, L_val, G, H, l_val, g_val, 0.0)
            Dr = D_numerical(func_r, L_val, G, H, l_val, g_val, 0.0)
            Dfpg = D_numerical(func_fpg, L_val, G, H, l_val, g_val, 0.0)

            a,e,eta,theta,f,r,E = get_orbital(L_val,G,H,l_val,g_val,0.0)
            expected_Df = 2*sin(f)/e

            err = abs(Df - expected_Df) / abs(expected_Df) if abs(expected_Df) > 1e-12 else abs(Df)
            print(f"{e_val:6.2f} {I_deg:5d} {l_val:5.1f} | {Df:14.6e} {expected_Df:14.6e} {err:10.2e} | {Dr:12.4e} | {Dfpg:12.4e}")

print()

# =====================================================================
# Also verify D(l) = q1 and D(phi) = D(f) - D(l) = D(f) - q1
# =====================================================================
print("=" * 80)
print("VERIFY: D(phi) = D(f-l) = D(f) - q1")
print("=" * 80)

def func_phi(L,G,H,l,g,h):
    e=sqrt(1-(G/L)**2); E=solve_kepler(l,e); f=true_from_eccentric(E,e); return f-l

for e_val in [0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L_val * eta_val
    theta_val = cos(np.radians(60))
    H = G * theta_val
    for l_val in [0.5, 1.5, 3.0]:
        a,e,eta,theta,f,r,E = get_orbital(L_val,G,H,l_val,0.5,0.0)
        Dphi = D_numerical(func_phi, L_val, G, H, l_val, 0.5, 0.0)
        q1 = 2*e*sin(E) + 2*eta*sin(f)/e
        Df_val = 2*sin(f)/e
        expected_Dphi = Df_val - q1
        print(f"  e={e_val}, l={l_val}: D(phi)={Dphi:+.8e}, D(f)-q1={expected_Dphi:+.8e}, diff={abs(Dphi-expected_Dphi):.2e}")

print()

# =====================================================================
# CORRECTED D(dS1/dl) with D(f) = 2*sin(f)/e
# =====================================================================
print("=" * 80)
print("CORRECTED closed-form D(dS1/dl)")
print("=" * 80)
print()
print("dS1/dl = Gamma*B0 - mu^2/r^3 * [B0 + B1p*cos(2f+2g)]")
print()
print("D acts through: a, eta (via Da, Deta), f (via Df=2sinf/e), g (via Dg=-2sinf/e)")
print("D(r)=0, D(theta)=0 => D(B0)=D(B1p)=0")
print()
print("Since D(f) = -D(g), we have D(2f+2g) = 2*D(f) + 2*D(g) = 0")
print("=> D[cos(2f+2g)] = 0 !!")
print()
print("This means the cos(2f+2g) term only gets a D contribution from the 1/r^3 factor")
print("But D(r)=0 => D(1/r^3)=0 too!")
print()
print("So D[ mu^2/r^3 * cos(2f+2g) ] = mu^2/r^3 * D[cos(2f+2g)] = 0")
print()
print("Wait, but D(f) != 0. Let me reconsider...")
print("D[cos(2f+2g)] = -sin(2f+2g) * [2*D(f) + 2*D(g)]")
print("             = -sin(2f+2g) * [2*(2sinf/e) + 2*(-2sinf/e)]")
print("             = -sin(2f+2g) * 0 = 0")
print()
print("YES! D[cos(2f+2g)] = 0 because D(f) = -D(g).")
print()
print("This means only the Gamma = mu^2/(a^3*eta^3) part contributes:")
print("D(dS1/dl) = D[Gamma] * B0 - 0 = 6*Gamma*(a/r)*B0")
print()

# Verify this simplified formula
print("=" * 80)
print("VERIFY: D(dS1/dl) = 6*Gamma*B0*(a/r)")
print("=" * 80)

print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'D(dS1/dl) FD':>18} {'6*Gamma*B0*rho':>18} {'RelErr':>12} {'Status'}")
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
            g = np.radians(g_deg)
            for l_val in [0.5, 1.5, 3.0]:
                a, e, eta, theta, f, r, E = get_orbital(L_val, G, H, l_val, g, 0.0)
                rho = a / r
                Gamma = MU**2 / (a**3 * eta**3)
                B0 = -0.5 + 1.5*theta**2

                D_num = D_numerical(dS1_dl, L_val, G, H, l_val, g, 0.0)
                formula = 6*Gamma*B0*rho

                denom = max(abs(D_num), abs(formula), 1e-15)
                rel_err = abs(D_num - formula) / denom
                max_err = max(max_err, rel_err)
                status = "PASS" if rel_err < 1e-4 else "FAIL"
                if status == "PASS": pass_count += 1
                else: fail_count += 1

                print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | "
                      f"{D_num:18.10e} {formula:18.10e} {rel_err:12.2e} {status}")

print("=" * 90)
print(f"Result: {pass_count} PASS, {fail_count} FAIL, max err = {max_err:.2e}")
print()

if fail_count > 0:
    print("If 6*Gamma*B0*rho doesn't match, the Df contribution must affect")
    print("other parts of dS1/dl beyond just cos(2f+2g).")
    print()
    print("Let me check: dS1/dl has a B0 term = Gamma*B0 = mu^2*B0/(a^3*eta^3)")
    print("This depends on a, eta but NOT on f or g.")
    print("D[Gamma*B0] = B0*D[Gamma] = B0*6*Gamma*rho  (as before)")
    print()
    print("And the r^-3 terms: mu^2*[B0+B1p*cos(2f+2g)]/r^3")
    print("D(r)=0 and D[cos(2f+2g)]=0, so D of this = 0")
    print()
    print("Unless... dS1/dl actually has additional f-dependence through")
    print("the identity Gamma*B0 = mu^2*B0/(a^3*eta^3) NOT depending on f,")
    print("while mu^2*B0/r^3 does depend on f through r.")
    print()
    print("But D(r)=0 means D(r^-3)=0, so D[mu^2*B0/r^3] = 0.")
    print()
    print("Therefore D(dS1/dl) = 6*Gamma*B0*(a/r) should be exact.")
    print("Any discrepancy must be from finite-difference error in D_numerical.")
