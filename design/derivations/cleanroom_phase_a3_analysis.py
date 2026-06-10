"""
CLEANROOM PHASE A3: Complete analysis of δp₁ from first principles.

No assumptions carried forward. Every claim verified numerically.

QUESTION: What is the correct first-order oblateness correction δp₁
that appears in BH61 Eq(15) for the drag-oblateness coupling?

Two candidate computations:
  (A) Poisson bracket: -{p₁, S₁}
  (B) D operator: D(∂S₁/∂l)

We verify both numerically, compare them, and determine which (if either)
matches BH61 Eq(14).

All orbital mechanics from scratch. No reference to previous session results.
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

MU = 1.0

# ============================================================
# Core orbital mechanics (verified building blocks)
# ============================================================

def solve_kepler(M, e, tol=1e-15):
    E = M + e * sin(M)
    for _ in range(100):
        dE = (M - E + e * sin(E)) / (1 - e * cos(E))
        E += dE
        if abs(dE) < tol: break
    return E

def true_from_eccentric(E, e):
    return 2 * arctan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2))

def get_state(L, G, H, l, g, h, mu=MU):
    """Full orbital state from Delaunay variables."""
    a = L**2/mu
    eta = G/L
    e = sqrt(max(1-eta**2, 0))
    theta = H/G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a*(1 - e*cos(E))
    return a, e, eta, theta, E, f, r

# ============================================================
# S₁ and its partial derivatives (from Phase A, 72/72 verified)
# ============================================================

def S1(L, G, H, l, g, h, mu=MU):
    a, e, eta, theta, E, f, r = get_state(L, G, H, l, g, h, mu)
    Gamma = mu**2/(a**3 * eta**3)
    B0 = -0.5 + 1.5*theta**2
    s2 = 1 - theta**2
    phi = f - l
    return -Gamma*(B0*(phi + e*sin(f))
                   + 0.75*s2*sin(2*f+2*g)
                   + 0.75*e*s2*sin(f+2*g)
                   + 0.25*e*s2*sin(3*f+2*g))

def dS1_dl(L, G, H, l, g, h, mu=MU):
    """∂S₁/∂l from the homological equation (analytically exact)."""
    a, e, eta, theta, E, f, r = get_state(L, G, H, l, g, h, mu)
    Gamma = mu**2/(a**3 * eta**3)
    B0 = -0.5 + 1.5*theta**2
    B1p = 1.5*(1 - theta**2)
    return Gamma*B0 - mu**2/r**3 * (B0 + B1p*cos(2*f+2*g))

# ============================================================
# p₁ function
# ============================================================

def p1_func(L, G, H, l, g, h, mu=MU):
    a = L**2/mu
    e = sqrt(max(1-(G/L)**2, 0))
    E = solve_kepler(l, e)
    r = a*(1 - e*cos(E))
    return L*(2*a/r - 1)

# ============================================================
# BH61 drag functions p_k, q_k (Eq 5)
# ============================================================

def get_pq(L, G, H, l, g, h, mu=MU):
    """Return (p1, p2, p3, q1, q2, q3) from BH61 Eq(5)."""
    a, e, eta, theta, E, f, r = get_state(L, G, H, l, g, h, mu)
    p1 = L*(2*a/r - 1)
    p2 = G
    p3 = H
    if e > 1e-14:
        q1 = 2*e*sin(E) + 2*eta*sin(f)/e
        q2 = -2*sin(f)/e
    else:
        q1 = 0.0
        q2 = 0.0
    q3 = 0.0
    return p1, p2, p3, q1, q2, q3

# ============================================================
# COMPUTATION A: Poisson bracket {p₁, S₁} by finite differences
# ============================================================

def poisson_bracket_p1_S1(L, G, H, l, g, h, eps=1e-7):
    """
    {p₁, S₁} = Σ_j [∂p₁/∂l_j · ∂S₁/∂L_j - ∂p₁/∂L_j · ∂S₁/∂l_j]
    Computed entirely by central finite differences.
    """
    base = [L, G, H, l, g, h]
    result = 0.0
    for j in range(3):  # pairs (L,l), (G,g), (H,h)
        Lj, lj = j, j+3
        # ∂p₁/∂l_j
        bp = list(base); bp[lj] += eps
        bm = list(base); bm[lj] -= eps
        dp1_dlj = (p1_func(*bp) - p1_func(*bm))/(2*eps)
        # ∂S₁/∂L_j
        bp = list(base); bp[Lj] += eps
        bm = list(base); bm[Lj] -= eps
        dS1_dLj = (S1(*bp) - S1(*bm))/(2*eps)
        # ∂p₁/∂L_j
        dp1_dLj = (p1_func(*bp) - p1_func(*bm))/(2*eps)
        # ∂S₁/∂l_j
        bp = list(base); bp[lj] += eps
        bm = list(base); bm[lj] -= eps
        dS1_dlj = (S1(*bp) - S1(*bm))/(2*eps)

        result += dp1_dlj * dS1_dLj - dp1_dLj * dS1_dlj
    return result

# ============================================================
# COMPUTATION B: D(∂S₁/∂l) by finite differences on Delaunay
# Using D defined by DL_j = -p_j, Dl_j = q_j (BH61 Eq 10)
# ============================================================

def D_of_dS1dl_delaunay(L, G, H, l, g, h, eps=1e-7):
    """
    D(∂S₁/∂l) = Σ_j [DL_j · ∂(∂S₁/∂l)/∂L_j + Dl_j · ∂(∂S₁/∂l)/∂l_j]
    where DL_j = -p_j, Dl_j = q_j from BH61 Eq(5).
    """
    p1, p2, p3, q1, q2, q3 = get_pq(L, G, H, l, g, h)
    D_actions = [-p1, -p2, -p3, q1, q2, q3]
    base = [L, G, H, l, g, h]
    result = 0.0
    for i in range(6):
        if abs(D_actions[i]) < 1e-20: continue
        bp = list(base); bp[i] += eps
        bm = list(base); bm[i] -= eps
        deriv = (dS1_dl(*bp) - dS1_dl(*bm))/(2*eps)
        result += D_actions[i] * deriv
    return result

# ============================================================
# COMPUTATION C: D(∂S₁/∂l) via Cartesian velocity derivatives
# D = -Σ_i v_i · ∂/∂v_i  (fundamental definition)
# ============================================================

def delaunay_to_cartesian(L, G, H, l, g, h, mu=MU):
    """Delaunay → Cartesian (x,y,z,vx,vy,vz)."""
    a = L**2/mu; eta = G/L; e = sqrt(max(1-eta**2, 0))
    theta = H/G; I = np.arccos(np.clip(theta, -1, 1))
    E = solve_kepler(l, e); f = true_from_eccentric(E, e)
    r = a*(1 - e*cos(E)); p = a*(1-e**2); h_am = sqrt(mu*p)

    x_pf = r*cos(f); y_pf = r*sin(f)
    vx_pf = -mu/h_am * sin(f); vy_pf = mu/h_am * (e + cos(f))

    cO, sO = cos(h), sin(h); cI, sI = cos(I), sin(I)
    cw, sw = cos(g), sin(g)

    Px = cO*cw - sO*sw*cI; Py = -cO*sw - sO*cw*cI
    Qx = sO*cw + cO*sw*cI; Qy = -sO*sw + cO*cw*cI
    Wx = sw*sI;             Wy = cw*sI

    x  = Px*x_pf + Py*y_pf;  y  = Qx*x_pf + Qy*y_pf;  z  = Wx*x_pf + Wy*y_pf
    vx = Px*vx_pf + Py*vy_pf; vy = Qx*vx_pf + Qy*vy_pf; vz = Wx*vx_pf + Wy*vy_pf
    return np.array([x, y, z, vx, vy, vz])

def cartesian_to_delaunay(state, mu=MU):
    """Cartesian → Delaunay."""
    x, y, z, vx, vy, vz = state
    rv = np.array([x,y,z]); vv = np.array([vx,vy,vz])
    r = np.linalg.norm(rv); v = np.linalg.norm(vv)
    a = -mu/(v**2 - 2*mu/r)
    hv = np.cross(rv, vv); h = np.linalg.norm(hv)
    I = np.arccos(np.clip(hv[2]/h, -1, 1))
    nv = np.cross([0,0,1], hv); n = np.linalg.norm(nv)
    ev = np.cross(vv, hv)/mu - rv/r; e = np.linalg.norm(ev)

    if n > 1e-14:
        Om = np.arccos(np.clip(nv[0]/n, -1, 1))
        if nv[1] < 0: Om = 2*pi - Om
    else: Om = 0.0

    if n > 1e-14 and e > 1e-14:
        w = np.arccos(np.clip(np.dot(nv, ev)/(n*e), -1, 1))
        if ev[2] < 0: w = 2*pi - w
    else: w = 0.0

    if e > 1e-14:
        f = np.arccos(np.clip(np.dot(ev, rv)/(e*r), -1, 1))
        if np.dot(rv, vv) < 0: f = 2*pi - f
    else: f = 0.0

    E = 2*np.arctan2(sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2))
    M = E - e*sin(E)
    if M < 0: M += 2*pi

    L = sqrt(mu*a); G = L*sqrt(max(1-e**2, 0)); H = G*cos(I)
    return L, G, H, M, w, Om

def dS1_dl_from_cart(state, mu=MU):
    L, G, H, l, g, h = cartesian_to_delaunay(state, mu)
    return dS1_dl(L, G, H, l, g, h, mu)

def D_of_dS1dl_cartesian(L, G, H, l, g, h, eps=1e-7):
    """D(∂S₁/∂l) using Cartesian definition: D = -Σ v_i · ∂/∂v_i."""
    state = delaunay_to_cartesian(L, G, H, l, g, h)
    result = 0.0
    for i in range(3, 6):  # velocity components only
        sp = state.copy(); sp[i] += eps
        sm = state.copy(); sm[i] -= eps
        dF_dvi = (dS1_dl_from_cart(sp) - dS1_dl_from_cart(sm))/(2*eps)
        result += -state[i] * dF_dvi
    return result

# ============================================================
# BH61 Eq(14) formula for δp₁ (as transcribed)
# ============================================================

def bh61_eq14_dp1(L, G, H, l, g, h, mu=MU):
    """BH61 Eq(14) δp₁ with k₂=1 factored out.
    δp₁ = 3(μ²/L³){B₀[-η⁻³ + ρ³(1-2ρ)] + B₁ρ³(1-2ρ)cos(2g+2f)}
    """
    a, e, eta, theta, E, f, r = get_state(L, G, H, l, g, h, mu)
    B0 = -0.5 + 1.5*theta**2
    B1 = 1.5 - 1.5*theta**2
    rho = a/r
    prefactor = 3.0*mu**2/L**3
    return prefactor*(B0*(-1/eta**3 + rho**3*(1-2*rho))
                     + B1*rho**3*(1-2*rho)*cos(2*g+2*f))

# ============================================================
# PHASE 1: Verify D-action identities
# ============================================================

def run_phase1():
    print("=" * 90)
    print("PHASE 1: D-ACTION IDENTITIES")
    print("Verify DL=-p₁, DG=-G, DH=-H, Dl=q₁, Dg=q₂ produce correct Da, De, Dr, Df, Dg")
    print("=" * 90)

    def D_of(func, L, G, H, l, g, h, eps=1e-7):
        p1, p2, p3, q1, q2, q3 = get_pq(L, G, H, l, g, h)
        D_acts = [-p1, -p2, -p3, q1, q2, q3]
        base = [L, G, H, l, g, h]
        result = 0.0
        for i in range(6):
            if abs(D_acts[i]) < 1e-20: continue
            bp = list(base); bp[i] += eps
            bm = list(base); bm[i] -= eps
            result += D_acts[i] * (func(*bp) - func(*bm))/(2*eps)
        return result

    def fa(L,G,H,l,g,h): return L**2/MU
    def fe(L,G,H,l,g,h): return sqrt(max(1-(G/L)**2, 0))
    def feta(L,G,H,l,g,h): return G/L
    def ftheta(L,G,H,l,g,h): return H/G
    def fr(L,G,H,l,g,h):
        a=L**2/MU; e=sqrt(max(1-(G/L)**2,0)); E=solve_kepler(l,e); return a*(1-e*cos(E))
    def ff(L,G,H,l,g,h):
        e=sqrt(max(1-(G/L)**2,0)); E=solve_kepler(l,e); return true_from_eccentric(E,e)

    L0 = 1.0
    print(f"\n{'e':>6} {'I':>5} {'l':>5} | {'Da num':>12} {'Da exp':>12} | {'De num':>12} {'De exp':>12} | {'Dr':>10} | {'Df num':>12} {'2sf/e':>12} | {'Df+Dg':>10}")
    print("-" * 120)

    for e_val in [0.01, 0.1, 0.3]:
        eta_val = sqrt(1-e_val**2); G0 = L0*eta_val
        for I_deg in [60]:
            theta_val = cos(np.radians(I_deg)); H0 = G0*theta_val
            for l_val in [0.5, 1.5, 3.0]:
                g0 = 0.7
                a,e,eta,theta,E,f,r = get_state(L0,G0,H0,l_val,g0,0.0)

                Da = D_of(fa, L0,G0,H0,l_val,g0,0.0)
                De = D_of(fe, L0,G0,H0,l_val,g0,0.0)
                Dr = D_of(fr, L0,G0,H0,l_val,g0,0.0)
                Df = D_of(ff, L0,G0,H0,l_val,g0,0.0)

                Da_exp = -2*a*(2*a/r - 1)
                De_exp = -2*(e + cos(f))
                Df_exp = 2*sin(f)/e if e > 1e-14 else 0.0
                Dg_val = -2*sin(f)/e if e > 1e-14 else 0.0  # = q₂ by definition

                print(f"{e_val:6.2f} {I_deg:5d} {l_val:5.1f} | "
                      f"{Da:12.6f} {Da_exp:12.6f} | "
                      f"{De:12.6f} {De_exp:12.6f} | "
                      f"{Dr:10.2e} | "
                      f"{Df:12.6f} {Df_exp:12.6f} | "
                      f"{Df+Dg_val:10.2e}")

    print("\nKey findings:")
    print("  Da, De: match expected formulas")
    print("  Dr: should be ~0")
    print("  Df: should be 2*sin(f)/e (NOT zero)")
    print("  Df + Dg: should be ~0 (the critical identity)")

# ============================================================
# PHASE 2: Compute both candidates across full grid
# ============================================================

def run_phase2():
    print("\n" + "=" * 90)
    print("PHASE 2: COMPARE THREE COMPUTATIONS OF dp₁")
    print("  (A) -{p₁, S₁}  (negative Poisson bracket)")
    print("  (B) D(∂S₁/∂l)  (D operator, Delaunay FD)")
    print("  (C) D(∂S₁/∂l)  (D operator, Cartesian FD)")
    print("  (D) BH61 Eq(14)")
    print("=" * 90)

    L0 = 1.0
    print(f"\n{'e':>5} {'I':>4} {'g':>4} {'l':>4} | {'-{p1,S1}':>14} {'D_Delaunay':>14} {'D_Cartesian':>14} {'BH61':>14} | {'A/B':>8} {'B/C':>8} {'B/D':>8}")
    print("-" * 115)

    for e_val in [0.01, 0.1, 0.3]:
        eta_val = sqrt(1-e_val**2); G = L0*eta_val
        for I_deg in [30, 60, 85]:
            theta_val = cos(np.radians(I_deg)); H = G*theta_val
            for g_deg in [0, 45, 90]:
                gv = np.radians(g_deg)
                for l_val in [0.5, 1.5, 3.0]:
                    A = -poisson_bracket_p1_S1(L0, G, H, l_val, gv, 0.0)
                    B = D_of_dS1dl_delaunay(L0, G, H, l_val, gv, 0.0)
                    C = D_of_dS1dl_cartesian(L0, G, H, l_val, gv, 0.0)
                    D = bh61_eq14_dp1(L0, G, H, l_val, gv, 0.0)

                    rAB = A/B if abs(B) > 1e-15 else float('inf')
                    rBC = B/C if abs(C) > 1e-15 else float('inf')
                    rBD = B/D if abs(D) > 1e-15 else float('inf')

                    print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:4.1f} | "
                          f"{A:14.8f} {B:14.8f} {C:14.8f} {D:14.8f} | "
                          f"{rAB:8.4f} {rBC:8.4f} {rBD:8.4f}")

    print()
    print("INTERPRETATION:")
    print("  A/B ≈ 1.0: Poisson bracket and D operator agree (same quantity)")
    print("  A/B ≠ 1.0: They compute different things")
    print("  B/C ≈ 1.0: Delaunay D and Cartesian D agree (operator defined correctly)")
    print("  B/C ≠ 1.0: D operator definition error")
    print("  B/D ≈ 1.0: D operator matches BH61 Eq(14)")
    print("  B/D ≠ 1.0: BH61 Eq(14) is wrong")

# ============================================================
# PHASE 3: Closed-form candidates and verification
# ============================================================

def run_phase3():
    print("\n" + "=" * 90)
    print("PHASE 3: TEST CLOSED-FORM CANDIDATES")
    print("=" * 90)

    L0 = 1.0

    # Candidate 1: 6*Gamma*B0*(a/r) — from D-operator with Dr=0, D(2f+2g)=0
    # Candidate 2: More complex expression if Candidate 1 fails

    print(f"\n{'e':>5} {'I':>4} {'g':>4} {'l':>4} | {'D_Delaunay':>14} {'6*G*B0*rho':>14} {'RelErr':>10} {'Status'}")
    print("-" * 75)

    max_err = 0.0
    pass_ct = fail_ct = 0
    for e_val in [0.01, 0.1, 0.3]:
        eta_val = sqrt(1-e_val**2); G = L0*eta_val
        for I_deg in [30, 60, 85]:
            theta_val = cos(np.radians(I_deg)); H = G*theta_val
            for g_deg in [0, 45, 90]:
                gv = np.radians(g_deg)
                for l_val in [0.5, 1.5, 3.0]:
                    a, e, eta, theta, E, f, r = get_state(L0, G, H, l_val, gv, 0.0)
                    Gamma = MU**2/(a**3*eta**3)
                    B0 = -0.5 + 1.5*theta**2
                    rho = a/r

                    B_num = D_of_dS1dl_delaunay(L0, G, H, l_val, gv, 0.0)
                    cand = 6*Gamma*B0*rho

                    rel = abs(B_num - cand)/max(abs(B_num), 1e-15)
                    max_err = max(max_err, rel)
                    st = "PASS" if rel < 1e-4 else "FAIL"
                    if st == "PASS": pass_ct += 1
                    else: fail_ct += 1

                    print(f"{e_val:5.2f} {I_deg:4d} {g_deg:4d} {l_val:4.1f} | "
                          f"{B_num:14.8f} {cand:14.8f} {rel:10.2e} {st}")

    print(f"\n6*Gamma*B0*rho: {pass_ct} PASS, {fail_ct} FAIL, max err = {max_err:.2e}")

# ============================================================
# PHASE 4: Harmonic decomposition of each quantity
# ============================================================

def run_phase4():
    print("\n" + "=" * 90)
    print("PHASE 4: HARMONIC STRUCTURE IN g")
    print("For each quantity, decompose into Fourier modes in g")
    print("to check for cos(2g) long-period content")
    print("=" * 90)

    L0 = 1.0
    N_g = 64

    for e_val in [0.1, 0.3]:
        eta_val = sqrt(1-e_val**2); G = L0*eta_val
        for I_deg in [60]:
            theta_val = cos(np.radians(I_deg)); H = G*theta_val
            l_val = 1.5

            g_vals = np.linspace(0, 2*pi, N_g, endpoint=False)

            vals_A = []; vals_B = []; vals_D = []
            for gv in g_vals:
                vals_A.append(-poisson_bracket_p1_S1(L0, G, H, l_val, gv, 0.0))
                vals_B.append(D_of_dS1dl_delaunay(L0, G, H, l_val, gv, 0.0))
                vals_D.append(bh61_eq14_dp1(L0, G, H, l_val, gv, 0.0))

            vals_A = np.array(vals_A)
            vals_B = np.array(vals_B)
            vals_D = np.array(vals_D)

            print(f"\ne={e_val}, I={I_deg}, l={l_val}:")
            print(f"  {'':>20} {'mean':>12} {'|cos(2g)|':>12} {'|sin(2g)|':>12} {'|cos(g)|':>12}")
            for name, vals in [("-{p1,S1}", vals_A), ("D_Delaunay", vals_B), ("BH61 Eq14", vals_D)]:
                mean = np.mean(vals)
                c2 = 2*np.mean(vals*cos(2*g_vals))
                s2 = 2*np.mean(vals*sin(2*g_vals))
                c1 = 2*np.mean(vals*cos(g_vals))
                print(f"  {name:>20} {mean:12.6f} {abs(c2):12.6f} {abs(s2):12.6f} {abs(c1):12.6f}")

# ============================================================
# RUN ALL PHASES
# ============================================================

if __name__ == "__main__":
    run_phase1()
    run_phase2()
    run_phase3()
    run_phase4()
