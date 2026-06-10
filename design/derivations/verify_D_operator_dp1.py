"""
CLEANROOM DERIVATION: dp1 = D(dS1/dl) via the velocity-homogeneity operator.

D is defined by its action on Delaunay variables (BH61 Eq 10):
  DL = -p1 = -L(2a/r - 1)
  DG = -p2 = -G
  DH = -p3 = -H
  Dl = q1 = 2e*sin(E) + 2*eta*sin(f)/e
  Dg = q2 = -(2/e)*sin(f)
  Dh = q3 = 0

For any function F(L, G, H, l, g, h):
  D(F) = DL*dF/dL + DG*dF/dG + DH*dF/dH + Dl*dF/dl + Dg*dF/dg + Dh*dF/dh

Step 1: Verify D-action identities (Da, De, Dr, Df, Dg, etc.)
Step 2: Compute D(dS1/dl) two ways:
  (a) Directly from D acting on Delaunay partials of dS1/dl
  (b) Via chain rule through orbital elements using D-action identities
Step 3: Compare both against the Phase A closed-form result
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

def dS1_dl(L, G, H, l, g, h, mu=MU):
    """Compute dS1/dl from the homological equation."""
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)
    Gamma = mu**2 / (a**3 * eta**3)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    return Gamma * B0 - mu**2 / r**3 * (B0 + B1p * cos(2*f + 2*g))


# =====================================================================
# STEP 1: Verify D-action identities
# =====================================================================
print("=" * 80)
print("STEP 1: Verify D-action identities")
print("=" * 80)

def D_numerical(func, L, G, H, l, g, h, eps=1e-7, mu=MU):
    """Compute D(F) by finite differences using the D-operator definition.

    D(F) = DL*dF/dL + DG*dF/dG + DH*dF/dH + Dl*dF/dl + Dg*dF/dg

    where DL=-p1, DG=-p2=-G, DH=-p3=-H, Dl=q1, Dg=q2, Dh=0.
    """
    a, e, eta, theta, f, r, E = get_orbital(L, G, H, l, g, h, mu)

    # D-actions on Delaunay variables (BH61 Eq 5, 10)
    p1 = L * (2*a/r - 1)
    p2 = G
    p3 = H
    q1 = 2*e*sin(E) + 2*eta*sin(f)/e if e > 1e-14 else 0.0
    q2 = -2*sin(f)/e if e > 1e-14 else 0.0
    # q3 = 0

    DL = -p1
    DG = -p2
    DH = -p3
    Dl = q1
    Dg = q2
    # Dh = 0

    # Compute partial derivatives of F by finite differences
    base = [L, G, H, l, g, h]
    f0 = func(*base)

    partials = []
    D_actions = [DL, DG, DH, Dl, Dg, 0.0]
    for i in range(6):
        bp = list(base); bp[i] += eps
        bm = list(base); bm[i] -= eps
        partials.append((func(*bp) - func(*bm)) / (2*eps))

    result = sum(d * p for d, p in zip(D_actions, partials))
    return result


# Test case
e_val = 0.1
a_val = 1.0
L_val = sqrt(MU * a_val)
eta_val = sqrt(1 - e_val**2)
G_val = L_val * eta_val
I_deg = 60.0
theta_val = cos(np.radians(I_deg))
H_val = G_val * theta_val
l_val = 1.5
g_val = np.radians(45.0)
h_val = 0.0

a_t, e_t, eta_t, theta_t, f_t, r_t, E_t = get_orbital(L_val, G_val, H_val, l_val, g_val, h_val)
rho_t = a_t / r_t

print(f"\nTest point: e={e_val}, I={I_deg}, g=45, l=1.5")
print(f"  a={a_t:.6f}, r={r_t:.6f}, rho=a/r={rho_t:.6f}")
print(f"  f={f_t:.6f} rad = {np.degrees(f_t):.2f} deg")
print(f"  E={E_t:.6f} rad")
print()

# Verify D(a), D(e), D(eta), D(r), D(f), D(g)
def func_a(L, G, H, l, g, h): return L**2 / MU
def func_e(L, G, H, l, g, h): return sqrt(1 - (G/L)**2)
def func_eta(L, G, H, l, g, h): return G / L
def func_theta(L, G, H, l, g, h): return H / G
def func_r(L, G, H, l, g, h):
    a = L**2/MU; e = sqrt(1-(G/L)**2); E = solve_kepler(l, e); return a*(1-e*cos(E))
def func_f(L, G, H, l, g, h):
    e = sqrt(1-(G/L)**2); E = solve_kepler(l, e); return true_from_eccentric(E, e)

# Expected D-actions from Phase A:
Da_expected = -2*a_t*(2*a_t/r_t - 1)
De_expected = -2*(e_t + cos(f_t))
Deta_expected = 2*eta_t*(a_t/r_t - 1)
Dtheta_expected = 0.0
Dr_expected = 0.0  # Phase A claim
Df_expected = 0.0  # Phase A claim
Dg_expected = -2*sin(f_t)/e_t

Da_num = D_numerical(func_a, L_val, G_val, H_val, l_val, g_val, h_val)
De_num = D_numerical(func_e, L_val, G_val, H_val, l_val, g_val, h_val)
Deta_num = D_numerical(func_eta, L_val, G_val, H_val, l_val, g_val, h_val)
Dtheta_num = D_numerical(func_theta, L_val, G_val, H_val, l_val, g_val, h_val)
Dr_num = D_numerical(func_r, L_val, G_val, H_val, l_val, g_val, h_val)
Df_num = D_numerical(func_f, L_val, G_val, H_val, l_val, g_val, h_val)

print("D-action identities:")
print(f"  D(a):     expected={Da_expected:+.8e}, numerical={Da_num:+.8e}, match={abs(Da_expected-Da_num)<1e-4*abs(Da_expected)}")
print(f"  D(e):     expected={De_expected:+.8e}, numerical={De_num:+.8e}, match={abs(De_expected-De_num)<1e-4*abs(De_expected)}")
print(f"  D(eta):   expected={Deta_expected:+.8e}, numerical={Deta_num:+.8e}, match={abs(Deta_expected-Deta_num)<1e-4*abs(Deta_expected) if abs(Deta_expected)>1e-12 else abs(Deta_num)<1e-6}")
print(f"  D(theta): expected={Dtheta_expected:+.8e}, numerical={Dtheta_num:+.8e}, match={abs(Dtheta_num)<1e-6}")
print(f"  D(r):     expected={Dr_expected:+.8e}, numerical={Dr_num:+.8e}, match={abs(Dr_num)<1e-6}")
print(f"  D(f):     expected={Df_expected:+.8e}, numerical={Df_num:+.8e}, match={abs(Df_num)<1e-6}")
print(f"  D(g) [q2]:expected={Dg_expected:+.8e}, (input, not computed)")

# =====================================================================
# STEP 2: Compute D(dS1/dl) numerically
# =====================================================================
print()
print("=" * 80)
print("STEP 2: D(dS1/dl) - numerical computation")
print("=" * 80)

D_dS1dl_num = D_numerical(dS1_dl, L_val, G_val, H_val, l_val, g_val, h_val)
print(f"\nD(dS1/dl) [numerical, FD on Delaunay] = {D_dS1dl_num:+.10e}")

# =====================================================================
# STEP 3: Phase A closed-form result
# =====================================================================
Gamma_t = MU**2 / (a_t**3 * eta_t**3)
B0_t = -0.5 + 1.5*theta_t**2
B1p_t = 1.5*(1-theta_t**2)

# Phase A formula: dp1 = 6*Gamma*B0*rho - 4*Gamma*eta^3*rho^3*B1p*sinf*sin(2f+2g)/e
phase_a = 6*Gamma_t*B0_t*rho_t - 4*Gamma_t*eta_t**3*rho_t**3*B1p_t*sin(f_t)*sin(2*f_t+2*g_val)/e_t

print(f"Phase A closed form              = {phase_a:+.10e}")
print(f"Ratio numerical/Phase_A          = {D_dS1dl_num/phase_a:.8f}")

# =====================================================================
# STEP 4: Full grid verification
# =====================================================================
print()
print("=" * 80)
print("STEP 4: Full grid verification of D(dS1/dl)")
print("=" * 80)
print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'D(dS1/dl) FD':>18} {'Phase A':>18} {'RelErr':>12} {'Status'}")
print("-" * 85)

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
                B1p = 1.5*(1-theta**2)

                # Numerical
                D_num = D_numerical(dS1_dl, L_val, G, H, l_val, g, 0.0)

                # Phase A closed form
                pa = 6*Gamma*B0*rho - 4*Gamma*eta**3*rho**3*B1p*sin(f)*sin(2*f+2*g)/e

                denom = max(abs(D_num), abs(pa), 1e-15)
                rel_err = abs(D_num - pa) / denom
                max_err = max(max_err, rel_err)
                status = "PASS" if rel_err < 1e-4 else "FAIL"
                if status == "PASS": pass_count += 1
                else: fail_count += 1

                print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | "
                      f"{D_num:18.10e} {pa:18.10e} {rel_err:12.2e} {status}")

print("=" * 85)
print(f"D(dS1/dl) verification: {pass_count} PASS, {fail_count} FAIL, max err = {max_err:.2e}")
print()

# =====================================================================
# STEP 5: If Dr or Df are nonzero, compute the correction
# =====================================================================
print("=" * 80)
print("STEP 5: Check if D(r)=0 and D(f)=0 hold across test grid")
print("=" * 80)

max_Dr = 0.0
max_Df = 0.0
for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L_val * eta_val
    for I_deg in [30, 60, 85]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for l_val in [0.5, 1.5, 3.0, 5.0]:
            Dr = D_numerical(func_r, L_val, G, H, l_val, 0.5, 0.0)
            Df = D_numerical(func_f, L_val, G, H, l_val, 0.5, 0.0)
            max_Dr = max(max_Dr, abs(Dr))
            max_Df = max(max_Df, abs(Df))

print(f"  Max |D(r)| across grid: {max_Dr:.6e}")
print(f"  Max |D(f)| across grid: {max_Df:.6e}")
if max_Dr < 1e-5 and max_Df < 1e-5:
    print("  => D(r)=0 and D(f)=0 CONFIRMED")
else:
    print("  => D(r) and/or D(f) are NONZERO!")
    print("     Phase A assumption is WRONG. Need to include Dr, Df terms.")
