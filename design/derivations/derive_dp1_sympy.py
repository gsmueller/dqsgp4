"""
CLEANROOM DERIVATION: Closed-form expression for dp1 = -{p1, S1}

Algebraically reduce the verified Poisson bracket to a closed form
in terms of (a/r, eta, theta, f, g).

All inputs are our independently derived and numerically verified
partial derivatives. No reference to BH61 coefficients.
"""
import sympy as sp
from sympy import symbols, sqrt, cos, sin, Rational, simplify, expand_trig
from sympy import collect, factor, trigsimp, cancel, together, apart

# Symbolic variables
rho = symbols('rho', positive=True)     # a/r
e, eta = symbols('e eta', positive=True)
theta = symbols('theta')                 # cos I
f, g = symbols('f g', real=True)
L, G, mu = symbols('L G mu', positive=True)
Gamma = symbols('Gamma', positive=True)  # mu^2/(a^3 eta^3)

# Inclination functions
B0 = Rational(-1, 2) + Rational(3, 2) * theta**2
B1p = Rational(3, 2) * (1 - theta**2)
s2 = 1 - theta**2

# Key relation: 1 + e*cos(f) = eta^2 * rho
ecosf = eta**2 * rho - 1  # e*cos(f)

print("=" * 70)
print("DERIVATION OF dp1 = -{p1, S1}")
print("=" * 70)
print()

# ================================================================
# Partial derivatives of p1
# ================================================================
print("--- Partial derivatives of p1 ---")

# dp1/dl = -2*L*rho^2*e*sin(f)/eta
# dp1/dL = 2*rho + 2*rho^2*eta^2*cos(f)/e - 1
#         = 2*rho + 2*rho^2*(eta^2*rho - 1)/e - 1   [using ecosf = eta^2*rho - 1... no]
# Actually: cos(f) appears, but e*cos(f) = eta^2*rho - 1, so cos(f) = (eta^2*rho - 1)/e

dp1_dl_over_L = -2 * rho**2 * e * sin(f) / eta

cosf_expr = (eta**2 * rho - 1) / e

dp1_dL = 2*rho + 2*rho**2*eta**2*cosf_expr/e - 1
dp1_dL = sp.expand(dp1_dL)
print(f"dp1/dL = {dp1_dL}")

dp1_dG = -2*rho**2*eta*cosf_expr/e
dp1_dG = sp.expand(dp1_dG)
print(f"dp1/dG = {dp1_dG}")

# ================================================================
# Partial derivatives of S1
# ================================================================
print("\n--- Partial derivatives of S1 ---")

# dS1/dl = Gamma * [B0 - eta^3*rho^3*(B0 + B1p*cos(2f+2g))]
dS1_dl_over_Gamma = B0 - eta**3 * rho**3 * (B0 + B1p*cos(2*f+2*g))
print(f"dS1/dl / Gamma = {dS1_dl_over_Gamma}")

# dS1/dg = -Gamma * B1p * [cos(2f+2g) + e*cos(f+2g) + (e/3)*cos(3f+2g)]
dS1_dg_over_Gamma = -B1p * (cos(2*f+2*g) + e*cos(f+2*g) + Rational(1,3)*e*cos(3*f+2*g))
print(f"dS1/dg / Gamma = {dS1_dg_over_Gamma}")

# dS1/dL: this is the hard one
# dS1/dL = (Gamma/L) * [3*B_bracket - sin(f)*(2+e*cos(f))/e * Phi1 - eta^2/e * Phi2]
#
# where:
#   B_bracket = B0*(phi+e*sin(f)) + (B1p/2)*sin(2f+2g) + (e*B1p/2)*sin(f+2g) + (e*B1p/6)*sin(3f+2g)
#   Phi1 = eta^2*rho*(B0 + B1p*cos(2f+2g))  [= (1+ecosf)*(B0 + B1p*cos(2f+2g))]
#   Phi2 = B0*sin(f) + (B1p/2)*sin(f+2g) + (B1p/6)*sin(3f+2g)

Phi1 = eta**2 * rho * (B0 + B1p*cos(2*f+2*g))
Phi2 = B0*sin(f) + Rational(1,2)*B1p*sin(f+2*g) + Rational(1,6)*B1p*sin(3*f+2*g)

# B_bracket (the curly bracket in S1 / (-Gamma))
phi = symbols('phi')  # phi = f - l, but appears only in B_bracket which we don't need directly
B_bracket = B0*(phi + e*sin(f)) + Rational(1,2)*B1p*sin(2*f+2*g) + Rational(1,2)*e*B1p*sin(f+2*g) + Rational(1,6)*e*B1p*sin(3*f+2*g)

# The full dS1/dL / (Gamma/L) is:
# 3*B_bracket - (sin(f)*(2+ecosf)/e)*Phi1 - (eta^2/e)*Phi2
# Using 2 + ecosf = 2 + eta^2*rho - 1 = 1 + eta^2*rho:
two_plus_ecosf = 1 + eta**2 * rho

dS1_dL_core = 3*B_bracket - sin(f)*two_plus_ecosf/e * Phi1 - eta**2/e * Phi2

# ================================================================
# Assemble the Poisson bracket terms
# ================================================================
print("\n--- Assembling dp1 = -{p1, S1} ---")
print("dp1 = dp1/dL * dS1/dl - dp1/dl * dS1/dL + dp1/dG * dS1/dg")
print()

# All three pieces, divided by Gamma:
# Piece 1: (dp1/dL) * (dS1/dl / Gamma) = dp1_dL * dS1_dl_over_Gamma
piece1 = sp.expand(dp1_dL * dS1_dl_over_Gamma)

# Piece 2: -(dp1/dl) * (dS1/dL / Gamma)
# dp1/dl = L * dp1_dl_over_L
# dS1/dL = (Gamma/L) * dS1_dL_core
# Product: L * dp1_dl_over_L * (Gamma/L) * dS1_dL_core = Gamma * dp1_dl_over_L * dS1_dL_core
# So piece2 / Gamma = -dp1_dl_over_L * dS1_dL_core / L... no wait.
#
# Let me be more careful.
# dp1 / Gamma = dp1_dL * dS1_dl/Gamma - (dp1_dl) * (dS1_dL) / Gamma + dp1_dG * dS1_dg/Gamma
#
# For piece 2: dp1_dl * dS1_dL / Gamma = (L * dp1_dl_over_L) * (Gamma/L * dS1_dL_core) / Gamma
#            = dp1_dl_over_L * dS1_dL_core
piece2 = -dp1_dl_over_L * dS1_dL_core

# Piece 3: dp1/dG * dS1/dg / Gamma = dp1_dG * dS1_dg_over_Gamma
piece3 = sp.expand(dp1_dG * dS1_dg_over_Gamma)

print("Piece 1 (dp1/dL * dS1/dl / Gamma):")
print(f"  = {piece1}")
print()
print("Piece 3 (dp1/dG * dS1/dg / Gamma):")
print(f"  = {piece3}")
print()

# ================================================================
# Focus on Pieces 1 and 3 first (no phi dependence)
# ================================================================
print("=" * 70)
print("SIMPLIFYING PIECES 1 AND 3 (closed form, no phi)")
print("=" * 70)

# Piece 1 expanded (substitute cosf = (eta^2*rho-1)/e where it appears in dp1_dL)
# dp1_dL = 2*rho + 2*rho^2*(eta^2*rho-1)^2/(e^2) ... no wait
# dp1_dL already expanded is in terms of rho and eta

# Let me re-derive dp1_dL more carefully
# dp1_dL = 2*rho + 2*rho^2*eta^2*cosf/e - 1
# With cosf = (eta^2*rho-1)/e:
# = 2*rho + 2*rho^2*(eta^2*rho-1) - 1... wait no
# 2*rho^2*eta^2*cosf/e = 2*rho^2*eta^2*(eta^2*rho-1)/(e*e) = 2*rho^2*(eta^2*rho-1)*eta^2/e^2
# Hmm, but we also have eta^2 = 1-e^2

# Actually, the cleaner approach: since ecosf = eta^2*rho - 1:
# dp1_dL = 2*rho + 2*rho^2*(ecosf)*eta/e ... no

# Let me just keep cosf as (eta^2*rho-1)/e and see what Piece 1 gives.

# Piece 1 = dp1_dL * [B0 - eta^3*rho^3*(B0 + B1p*cos(2f+2g))]
# = dp1_dL * B0 - dp1_dL * eta^3*rho^3*B0 - dp1_dL * eta^3*rho^3*B1p*cos(2f+2g)

# The first two give the B0 part, the third gives the B1p*cos(2f+2g) part.

# Let me compute dp1_dL:
dp1_dL_v2 = 2*rho + 2*rho**2*(eta**2*rho - 1)/e * eta**2/e - 1
# Wait: dp1_dL = 2*rho + 2*rho^2*eta^2*cosf/e - 1
# = 2*rho + 2*rho^2*eta^2*(eta^2*rho-1)/(e*e) - 1
# Hmm that introduces e^2 in denominator. Let me use eta^2 = 1-e^2.

# Actually, simpler: keep cosf symbolic and use ecosf = eta^2*rho-1 only at the end.
# The product dp1_dL * dS1_dl involves:
# [2rho + 2rho^2*eta^2*cosf/e - 1] * [B0 - eta^3*rho^3*(B0+B1p*cos(2f+2g))]

# Let me use a substitution approach.
# Define c = cos(f), s = sin(f), c2 = cos(2f+2g), s2g = sin(2f+2g), etc.

c, s = symbols('c s')  # cos f, sin f
c2g = cos(2*f+2*g)
s2g = sin(2*f+2*g)
c1g = cos(f+2*g)
s1g = sin(f+2*g)
c3g = cos(3*f+2*g)
s3g = sin(3*f+2*g)

# Rewrite dp1_dL in terms of c = cosf
dp1_dL_sym = 2*rho + 2*rho**2*eta**2*c/e - 1
dp1_dG_sym = -2*rho**2*eta*c/e
dp1_dl_sym = -2*rho**2*e*s/eta  # this is dp1/dl / L

# dS1/dl / Gamma
dS1_dl_sym = B0 - eta**3*rho**3*(B0 + B1p*c2g)

# dS1/dg / Gamma
dS1_dg_sym = -B1p*(c2g + e*c1g + Rational(1,3)*e*c3g)

# For dS1/dL / (Gamma/L):
Phi1_sym = eta**2*rho*(B0 + B1p*c2g)
Phi2_sym = B0*s + Rational(1,2)*B1p*s1g + Rational(1,6)*B1p*s3g
two_plus_ec = 1 + eta**2*rho  # 2 + ecosf = 2 + (eta^2*rho - 1) = 1 + eta^2*rho

# dS1_dL_core = 3*B_bracket - s*two_plus_ec/e * Phi1_sym - eta^2/e * Phi2_sym
# B_bracket contains phi which we don't want. But note:
# Piece 2 = -dp1_dl_sym * dS1_dL_core (with appropriate factors)
# = -(dp1_dl_sym) * [3*B_bracket - s*two_plus_ec/e*Phi1 - eta^2/e*Phi2]
# = -dp1_dl_sym*3*B_bracket + dp1_dl_sym*s*two_plus_ec/e*Phi1 + dp1_dl_sym*eta^2/e*Phi2

# The B_bracket term involves phi = f-l. Let's separate it:
# -dp1_dl_sym * 3 * B_bracket = (-2rho^2*e*s/eta) * 3 * B_bracket
# This contains phi. But ultimately, phi doesn't appear in the final δp1
# because it cancels. Let me check this by computing the phi-dependent part.

# B_bracket = B0*(phi + e*s) + (B1p/2)*s2g + (e*B1p/2)*s1g + (e*B1p/6)*s3g

# The phi-dependent part of Piece 2:
# -dp1_dl_sym * 3 * B0 * phi = (-2rho^2*e*s/eta) * 3*B0*phi = -6*B0*rho^2*e*s*phi/eta

# But wait, does Piece 1 also have phi dependence? No - dS1/dl has no phi.
# And Piece 3 has no phi either.
# So the phi-dependent part is ONLY from Piece 2, and it must vanish.

# Actually, phi appears in the FULL S1, not in its derivatives w.r.t. l or g.
# The derivative dS1/dL involves d(B_bracket)/dL, which has a dphi/dL = df/dL term.
# But the B_bracket itself appears multiplied by dGamma/dL.

# Let me reconsider. The full expression is:
# dS1/dL = -dGamma/dL * B - Gamma * dB/dL
#         = (3Gamma/L)*B - Gamma*(Phi1*F + Phi2*D)
# where F = df/dL = s*(2+ec)/(eL), D = de/dL = eta^2/(eL)

# So dS1/dL / (Gamma/L) = 3*B - L*Phi1*F - L*Phi2*D
#                        = 3*B - s*(1+eta^2*rho)/e * Phi1 - eta^2/e * Phi2

# The B_bracket part 3*B contains phi. So:
# Piece 2 / Gamma = -dp1_dl_sym * dS1_dL/(Gamma/L) * (1/L... wait, let me get the factors right.

# dp1 / Gamma = dp1_dL * dS1_dl/Gamma - dp1_dl * dS1_dL/Gamma + dp1_dG * dS1_dg/Gamma
# dp1_dl = L * dp1_dl_sym
# dS1_dL = (Gamma/L) * dS1_dL_core
# dp1_dl * dS1_dL / Gamma = L * dp1_dl_sym * (Gamma/L) * dS1_dL_core / Gamma = dp1_dl_sym * dS1_dL_core

# So:
# dp1/Gamma = dp1_dL_sym * dS1_dl_sym - dp1_dl_sym * dS1_dL_core + dp1_dG_sym * dS1_dg_sym

# And dS1_dL_core = 3*B_bracket - s*(1+eta^2*rho)/e*Phi1 - eta^2/e*Phi2

# Separate the B_bracket: B_bracket = B0*phi + B0*e*s + (B1p/2)*s2g + (eB1p/2)*s1g + (eB1p/6)*s3g
# = B0*phi + B_rest

B_rest = B0*e*s + Rational(1,2)*B1p*s2g + Rational(1,2)*e*B1p*s1g + Rational(1,6)*e*B1p*s3g

# dS1_dL_core = 3*B0*phi + 3*B_rest - s*(1+eta^2*rho)/e*Phi1 - eta^2/e*Phi2

# Piece 2 = -dp1_dl_sym * dS1_dL_core
#          = -dp1_dl_sym * 3*B0*phi - dp1_dl_sym*(3*B_rest - s*(1+eta^2*rho)/e*Phi1 - eta^2/e*Phi2)

# The phi-dependent part: -dp1_dl_sym * 3*B0*phi = (2rho^2*e*s/eta) * 3*B0*phi

# This must cancel against something. But pieces 1 and 3 don't have phi!
# Unless... we need to account for the fact that phi itself has an l-dependence
# that produces terms when we take ∂/∂L (since df/dL != 0).
# Wait, phi = f - l, and dphi/dL = df/dL (since l is held constant).
# This is already accounted for in dB/dL = Phi1*F + Phi2*D.
# And dGamma/dL * B appears as 3*Gamma/L * B, which has the full B including phi.

# So the phi term in Piece 2 is: -dp1_dl_sym * 3*B0*phi = (6*rho^2*e*s*B0*phi)/eta

# For this to cancel in the total dp1, there must be a compensating phi term
# from Piece 1 or Piece 3. But Piece 1 involves dS1/dl which has no phi,
# and Piece 3 involves dS1/dg which also has no phi.

# Therefore, the phi term (6*rho^2*e*s*B0*phi/eta) MUST BE the secular/long-period
# contribution to dp1! It doesn't cancel - it's the "mixed secular term" that
# BH61 talks about. When averaged over l, phi = f-l has a nonzero average
# proportional to e. This is the coupling term between drag and oblateness.

# For the SHORT-PERIOD part, we can ignore phi (it's zero on average or produces
# secular terms). Let's compute the short-period part.

print("Phi-dependent part of dp1:")
print(f"  = (6*rho^2*e*sin(f)*B0/eta) * phi")
print("This is the SECULAR/LONG-PERIOD part (mixed drag-oblateness coupling).")
print("BH61 mentions this explicitly on p.198.")
print()

# ================================================================
# SHORT-PERIOD PART (phi-independent)
# ================================================================
print("=" * 70)
print("COMPUTING SHORT-PERIOD PART")
print("=" * 70)
print()

# Total dp1/Gamma (short-period) = Piece1 + Piece2_nophi + Piece3
# where Piece2_nophi = -dp1_dl_sym * (3*B_rest - s*(1+eta^2*rho)/e*Phi1 - eta^2/e*Phi2)

Piece2_nophi = -dp1_dl_sym * (3*B_rest - s*two_plus_ec/e * Phi1_sym - eta**2/e * Phi2_sym)

# dp1_sp / Gamma = Piece1 + Piece2_nophi + Piece3
total_sp = piece1 + sp.expand(Piece2_nophi) + piece3

# This is a trigonometric expression in f and g.
# Let me expand trig products and collect by harmonics.

print("Expanding trigonometric products...")

# Use expand_trig to expand all trig products
total_expanded = sp.expand(sp.expand_trig(total_sp))

# This will be very messy. Let me try a different approach:
# evaluate numerically at specific points and fit the coefficients.

print("\nDirect symbolic simplification is complex.")
print("Switching to numerical coefficient extraction...\n")

# ================================================================
# NUMERICAL COEFFICIENT EXTRACTION
# ================================================================
# We know dp1 = Gamma * [sum of terms with different rho^n * cos(jf+2kg)]
# Since only 2g harmonics appear, the general form is:
#
# dp1/Gamma = Sum_{n} a_n * rho^n  (B0-proportional part, no g dependence)
#           + Sum_{n} [b_n*cos(2f+2g) + c_n*sin(2f+2g)] * rho^n  (B1p-proportional part)
#           + constant_term * B0  (the -eta^{-3} part in BH61)
#
# Plus the phi-dependent secular term.

import numpy as np
from numpy import pi as PI

MU_NUM = 1.0

def solve_kepler_num(M, ecc, tol=1e-15):
    E = M + ecc * np.sin(M)
    for _ in range(100):
        dE = (M - E + ecc * np.sin(E)) / (1 - ecc * np.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def true_from_eccentric_num(E, ecc):
    return 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E / 2), np.sqrt(1 - ecc) * np.cos(E / 2))

def compute_dp1_numerical(L_val, G_val, H_val, l_val, g_val, h_val=0.0, eps=1e-7):
    """Compute dp1 = -{p1, S1} by finite differences."""
    mu = MU_NUM

    def p1_func(Lv, Gv, Hv, lv, gv, hv):
        av = Lv**2 / mu
        ev = np.sqrt(1 - (Gv/Lv)**2)
        Ev = solve_kepler_num(lv, ev)
        rv = av * (1 - ev * np.cos(Ev))
        return Lv * (2*av/rv - 1)

    def S1_func(Lv, Gv, Hv, lv, gv, hv):
        av = Lv**2 / mu
        etav = Gv / Lv
        ev = np.sqrt(1 - etav**2)
        thetav = Hv / Gv
        Ev = solve_kepler_num(lv, ev)
        fv = true_from_eccentric_num(Ev, ev)
        Gammav = mu**2 / (av**3 * etav**3)
        B0v = -0.5 + 1.5 * thetav**2
        s2v = 1 - thetav**2
        phiv = fv - lv
        return -Gammav * (
            B0v * (phiv + ev * np.sin(fv))
            + 0.75 * s2v * np.sin(2*fv + 2*gv)
            + 0.75 * ev * s2v * np.sin(fv + 2*gv)
            + 0.25 * ev * s2v * np.sin(3*fv + 2*gv)
        )

    vb = [L_val, G_val, H_val, l_val, g_val, h_val]
    result = 0.0
    for j in range(3):
        Lj, lj = j, j+3
        v_p = list(vb); v_p[lj] += eps
        v_m = list(vb); v_m[lj] -= eps
        dp1_dlj = (p1_func(*v_p) - p1_func(*v_m)) / (2*eps)

        v_p = list(vb); v_p[Lj] += eps
        v_m = list(vb); v_m[Lj] -= eps
        dS1_dLj = (S1_func(*v_p) - S1_func(*v_m)) / (2*eps)
        dp1_dLj = (p1_func(*v_p) - p1_func(*v_m)) / (2*eps)

        v_p = list(vb); v_p[lj] += eps
        v_m = list(vb); v_m[lj] -= eps
        dS1_dlj = (S1_func(*v_p) - S1_func(*v_m)) / (2*eps)

        result += dp1_dlj * dS1_dLj - dp1_dLj * dS1_dlj
    return -result  # Note the MINUS sign: dp1 = -{p1, S1}


def compute_dp1_pieces(e_val, theta_val, f_val, g_val, l_val):
    """Compute individual pieces of dp1/Gamma symbolically evaluated numerically."""
    mu = MU_NUM
    eta_val = np.sqrt(1 - e_val**2)
    a_val = 1.0
    L_val = np.sqrt(mu * a_val)
    G_val = L_val * eta_val
    H_val = G_val * theta_val
    r_val = a_val * eta_val**2 / (1 + e_val * np.cos(f_val))
    rho_val = a_val / r_val
    Gamma_val = mu**2 / (a_val**3 * eta_val**3)

    dp1 = compute_dp1_numerical(L_val, G_val, H_val, l_val, g_val)
    return dp1 / Gamma_val, rho_val


# Now let's extract the functional form by computing at many points
# and fitting to a polynomial in rho times trig functions of f and g.

print("=" * 70)
print("EXTRACTING FUNCTIONAL FORM OF dp1/Gamma")
print("=" * 70)
print()

# Strategy: at fixed (e, theta, g), compute dp1/Gamma for many values of l (hence f, rho).
# Then decompose as:
#   dp1/Gamma = F_0(rho) + F_c(rho)*cos(2f+2g) + F_s(rho)*sin(2f+2g) + secular*phi
#
# To separate these, compute at g=0 and g=pi/4:
#   dp1(g=0)/Gamma = F_0(rho) + F_c(rho)*cos(2f) + F_s(rho)*sin(2f)
#   dp1(g=pi/4)/Gamma = F_0(rho) - F_s(rho)*cos(2f) + F_c(rho)*sin(2f)... no wait
#
# Actually, easier: compute at three g values and solve for F_0, F_c, F_s.
# At g, dp1/Gamma = F_0 + F_c*cos(2f+2g) + F_s*sin(2f+2g)
# At g+pi/3, dp1/Gamma = F_0 + F_c*cos(2f+2g+2pi/3) + F_s*sin(2f+2g+2pi/3)
# At g+2pi/3, dp1/Gamma = F_0 + F_c*cos(2f+2g+4pi/3) + F_s*sin(2f+2g+4pi/3)
# Average of these three = F_0 (since cos/sin of equally spaced phases sum to 0)

for e_val in [0.1]:
    eta_val = np.sqrt(1 - e_val**2)
    a_val = 1.0
    L_val = np.sqrt(MU_NUM * a_val)
    G_val = L_val * eta_val

    for I_deg in [60]:
        theta_val = np.cos(np.radians(I_deg))
        H_val = G_val * theta_val
        B0_val = -0.5 + 1.5 * theta_val**2
        B1p_val = 1.5 * (1 - theta_val**2)
        Gamma_val = MU_NUM**2 / (a_val**3 * eta_val**3)

        print(f"e={e_val}, I={I_deg}, eta={eta_val:.6f}, B0={B0_val:.4f}, B1p={B1p_val:.4f}")
        print(f"Gamma = {Gamma_val:.6f}")
        print()

        # Compute at many l values, three g values each
        N_l = 200
        l_vals = np.linspace(0.01, 2*PI-0.01, N_l)

        print(f"{'l':>8} {'f':>8} {'rho':>8} {'F_0/B0':>12} {'F_c/B1p':>12} {'F_s/B1p':>12} {'rho^3(1-2rho)':>14}")
        print("-" * 80)

        for i_l in range(0, N_l, 20):  # Print every 20th point
            l_val = l_vals[i_l]
            E_val = solve_kepler_num(l_val, e_val)
            f_val = true_from_eccentric_num(E_val, e_val)
            r_val = a_val * (1 - e_val * np.cos(E_val))
            rho_val = a_val / r_val
            phi_val = f_val - l_val

            # Compute at three g values to separate F_0, F_c, F_s
            g_offsets = [0, 2*PI/3, 4*PI/3]
            g_base = 0.5  # arbitrary base g
            dp1_vals = []
            for dg in g_offsets:
                g_test = g_base + dg
                dp1_v = compute_dp1_numerical(L_val, G_val, H_val, l_val, g_test)
                dp1_vals.append(dp1_v / Gamma_val)

            F_0 = np.mean(dp1_vals)
            # F_c*cos(2f+2g) + F_s*sin(2f+2g)
            # At g_base+dg: residual = F_c*cos(2f+2g_base+2dg) + F_s*sin(2f+2g_base+2dg)
            angle_base = 2*f_val + 2*g_base
            resid = [dp1_vals[i] - F_0 for i in range(3)]
            # resid[i] = F_c*cos(angle_base+2*g_offsets[i]) + F_s*sin(angle_base+2*g_offsets[i])
            # Solve 2x2 system using first two residuals:
            c0 = np.cos(angle_base + 2*g_offsets[0])
            s0 = np.sin(angle_base + 2*g_offsets[0])
            c1 = np.cos(angle_base + 2*g_offsets[1])
            s1 = np.sin(angle_base + 2*g_offsets[1])
            det = c0*s1 - c1*s0
            F_c = (resid[0]*s1 - resid[1]*s0) / det
            F_s = (c0*resid[1] - c1*resid[0]) / det

            # BH61's prediction for comparison
            rho3_1m2rho = rho_val**3 * (1 - 2*rho_val)

            # Normalize: F_0 should be proportional to B0, F_c to B1p
            F_0_over_B0 = F_0 / B0_val if abs(B0_val) > 1e-10 else 0
            F_c_over_B1p = F_c / B1p_val if abs(B1p_val) > 1e-10 else 0
            F_s_over_B1p = F_s / B1p_val if abs(B1p_val) > 1e-10 else 0

            print(f"{l_val:8.4f} {f_val:8.4f} {rho_val:8.4f} "
                  f"{F_0_over_B0:12.6f} {F_c_over_B1p:12.6f} {F_s_over_B1p:12.6f} "
                  f"{rho3_1m2rho:14.6f}")

        print()
        print("If BH61 is right:")
        print("  F_0/B0 should equal -eta^{-3} + rho^3(1-2*rho) = " +
              f"-{1/eta_val**3:.6f} + rho^3(1-2*rho)")
        print("  F_c/B1p should equal rho^3(1-2*rho)")
        print("  F_s/B1p should be zero")
        print()

        # Now let's check: what IS the rho-dependence?
        # Plot F_0/B0 + eta^{-3} vs rho^3(1-2rho)
        print("Testing factored form...")
        print(f"{'rho':>8} {'F_0/B0+eta^-3':>14} {'rho^3(1-2rho)':>14} {'ratio':>10} {'F_c/B1p':>14} {'ratio2':>10}")
        print("-" * 80)

        for i_l in range(0, N_l, 20):
            l_val = l_vals[i_l]
            E_val = solve_kepler_num(l_val, e_val)
            f_val = true_from_eccentric_num(E_val, e_val)
            r_val = a_val * (1 - e_val * np.cos(E_val))
            rho_val = a_val / r_val

            g_offsets = [0, 2*PI/3, 4*PI/3]
            g_base = 0.5
            dp1_vals = []
            for dg in g_offsets:
                dp1_v = compute_dp1_numerical(L_val, G_val, H_val, l_val, g_base+dg)
                dp1_vals.append(dp1_v / Gamma_val)

            F_0 = np.mean(dp1_vals)
            angle_base = 2*f_val + 2*g_base
            resid = [dp1_vals[i] - F_0 for i in range(3)]
            c0 = np.cos(angle_base + 2*g_offsets[0])
            s0 = np.sin(angle_base + 2*g_offsets[0])
            c1 = np.cos(angle_base + 2*g_offsets[1])
            s1 = np.sin(angle_base + 2*g_offsets[1])
            det = c0*s1 - c1*s0
            F_c = (resid[0]*s1 - resid[1]*s0) / det
            F_s = (c0*resid[1] - c1*resid[0]) / det

            F_0_over_B0 = F_0 / B0_val
            F_c_over_B1p = F_c / B1p_val
            shortperiod = F_0_over_B0 + 1/eta_val**3
            rho3m2rho4 = rho_val**3 * (1 - 2*rho_val)

            r1 = shortperiod / rho3m2rho4 if abs(rho3m2rho4) > 1e-10 else float('inf')
            r2 = F_c_over_B1p / rho3m2rho4 if abs(rho3m2rho4) > 1e-10 else float('inf')

            print(f"{rho_val:8.4f} {shortperiod:14.6f} {rho3m2rho4:14.6f} {r1:10.4f} "
                  f"{F_c_over_B1p:14.6f} {r2:10.4f}")

        print()
        print("If BH61's factored form is correct, both ratio columns should be ~1.0")
        print("If they're not, BH61's rho-dependence is wrong.")
        print()

        # Let me also try other candidate rho-dependencies
        print("Testing alternative forms: rho^3*(-1 + 2*rho) and 3*rho^4...")
        print(f"{'rho':>8} {'F0sp':>12} {'rho^3(-1+2rho)':>14} {'r1':>8} {'3rho^4':>10} {'r2':>8} {'Fc/B1p':>12} {'rho3(-1+2rho)':>14} {'r3':>8}")
        print("-" * 110)

        for i_l in range(0, N_l, 20):
            l_val = l_vals[i_l]
            E_val = solve_kepler_num(l_val, e_val)
            f_val = true_from_eccentric_num(E_val, e_val)
            r_val = a_val * (1 - e_val * np.cos(E_val))
            rho_val = a_val / r_val

            g_offsets = [0, 2*PI/3, 4*PI/3]
            g_base = 0.5
            dp1_vals = []
            for dg in g_offsets:
                dp1_v = compute_dp1_numerical(L_val, G_val, H_val, l_val, g_base+dg)
                dp1_vals.append(dp1_v / Gamma_val)

            F_0 = np.mean(dp1_vals)
            angle_base = 2*f_val + 2*g_base
            resid = [dp1_vals[i] - F_0 for i in range(3)]
            c0v = np.cos(angle_base + 2*g_offsets[0])
            s0v = np.sin(angle_base + 2*g_offsets[0])
            c1v = np.cos(angle_base + 2*g_offsets[1])
            s1v = np.sin(angle_base + 2*g_offsets[1])
            det = c0v*s1v - c1v*s0v
            F_c = (resid[0]*s1v - resid[1]*s0v) / det

            F_0_over_B0 = F_0 / B0_val
            F_c_over_B1p = F_c / B1p_val
            shortperiod = F_0_over_B0 + 1/eta_val**3  # remove constant part

            cand1 = rho_val**3 * (-1 + 2*rho_val)  # flipped sign from BH61
            cand2 = 3 * rho_val**4
            cand3 = rho_val**3 * (-1 + 2*rho_val)

            r1 = shortperiod / cand1 if abs(cand1) > 1e-10 else float('inf')
            r2 = shortperiod / cand2 if abs(cand2) > 1e-10 else float('inf')
            r3 = F_c_over_B1p / cand3 if abs(cand3) > 1e-10 else float('inf')

            print(f"{rho_val:8.4f} {shortperiod:12.6f} {cand1:14.6f} {r1:8.4f} "
                  f"{cand2:10.6f} {r2:8.4f} {F_c_over_B1p:12.6f} {cand3:14.6f} {r3:8.4f}")

        print()

        # Try fitting to polynomial in rho
        print("Fitting F_0/B0 + eta^{-3} to a*rho^3 + b*rho^4 + c*rho^5...")
        rhos = []
        F0sps = []
        Fcb1s = []
        for i_l in range(N_l):
            l_val = l_vals[i_l]
            E_val = solve_kepler_num(l_val, e_val)
            f_val = true_from_eccentric_num(E_val, e_val)
            r_val = a_val * (1 - e_val * np.cos(E_val))
            rho_val = a_val / r_val

            g_offsets = [0, 2*PI/3, 4*PI/3]
            g_base = 0.5
            dp1_vals = []
            for dg in g_offsets:
                dp1_v = compute_dp1_numerical(L_val, G_val, H_val, l_val, g_base+dg)
                dp1_vals.append(dp1_v / Gamma_val)

            F_0 = np.mean(dp1_vals)
            angle_base = 2*f_val + 2*g_base
            resid = [dp1_vals[i] - F_0 for i in range(3)]
            c0v = np.cos(angle_base + 2*g_offsets[0])
            s0v = np.sin(angle_base + 2*g_offsets[0])
            c1v = np.cos(angle_base + 2*g_offsets[1])
            s1v = np.sin(angle_base + 2*g_offsets[1])
            det = c0v*s1v - c1v*s0v
            F_c = (resid[0]*s1v - resid[1]*s0v) / det

            rhos.append(rho_val)
            F0sps.append(F_0 / B0_val + 1/eta_val**3)
            Fcb1s.append(F_c / B1p_val)

        rhos = np.array(rhos)
        F0sps = np.array(F0sps)
        Fcb1s = np.array(Fcb1s)

        # Fit F_0/B0 + eta^{-3} = a*rho^3 + b*rho^4 + c*rho^5
        A_mat = np.column_stack([rhos**3, rhos**4, rhos**5])
        coeffs_F0, res_F0, _, _ = np.linalg.lstsq(A_mat, F0sps, rcond=None)
        print(f"  F_0/B0 + eta^-3 = {coeffs_F0[0]:.6f}*rho^3 + {coeffs_F0[1]:.6f}*rho^4 + {coeffs_F0[2]:.6f}*rho^5")
        print(f"  Residual norm: {np.sqrt(res_F0[0]) if len(res_F0)>0 else 'N/A'}")

        # Fit F_c/B1p = a*rho^3 + b*rho^4 + c*rho^5
        coeffs_Fc, res_Fc, _, _ = np.linalg.lstsq(A_mat, Fcb1s, rcond=None)
        print(f"  F_c/B1p = {coeffs_Fc[0]:.6f}*rho^3 + {coeffs_Fc[1]:.6f}*rho^4 + {coeffs_Fc[2]:.6f}*rho^5")
        print(f"  Residual norm: {np.sqrt(res_Fc[0]) if len(res_Fc)>0 else 'N/A'}")

        print()
        print(f"BH61 claims: 1*rho^3 + (-2)*rho^4 + 0*rho^5 for BOTH F_0/B0 and F_c/B1p")
        print(f"Our result:  F_0/B0: {coeffs_F0[0]:.4f}*rho^3 + {coeffs_F0[1]:.4f}*rho^4 + {coeffs_F0[2]:.4f}*rho^5")
        print(f"             F_c/B1p: {coeffs_Fc[0]:.4f}*rho^3 + {coeffs_Fc[1]:.4f}*rho^4 + {coeffs_Fc[2]:.4f}*rho^5")

        # Also check for sin(2f+2g) component
        Fss = []
        for i_l in range(N_l):
            l_val = l_vals[i_l]
            E_val = solve_kepler_num(l_val, e_val)
            f_val = true_from_eccentric_num(E_val, e_val)
            r_val = a_val * (1 - e_val * np.cos(E_val))
            rho_val = a_val / r_val

            g_offsets = [0, 2*PI/3, 4*PI/3]
            g_base = 0.5
            dp1_vals = []
            for dg in g_offsets:
                dp1_v = compute_dp1_numerical(L_val, G_val, H_val, l_val, g_base+dg)
                dp1_vals.append(dp1_v / Gamma_val)

            F_0 = np.mean(dp1_vals)
            angle_base = 2*f_val + 2*g_base
            resid = [dp1_vals[i] - F_0 for i in range(3)]
            c0v = np.cos(angle_base + 2*g_offsets[0])
            s0v = np.sin(angle_base + 2*g_offsets[0])
            c1v = np.cos(angle_base + 2*g_offsets[1])
            s1v = np.sin(angle_base + 2*g_offsets[1])
            det = c0v*s1v - c1v*s0v
            F_s = (c0v*resid[1] - c1v*resid[0]) / det
            Fss.append(F_s / B1p_val)

        Fss = np.array(Fss)
        print(f"\n  Max |F_s/B1p| (sin(2f+2g) component): {np.max(np.abs(Fss)):.6e}")
        print(f"  Mean |F_s/B1p|: {np.mean(np.abs(Fss)):.6e}")
        if np.max(np.abs(Fss)) < 1e-3 * np.max(np.abs(Fcb1s)):
            print("  => sin(2f+2g) component is NEGLIGIBLE (< 0.1% of cos component)")
        else:
            print("  => sin(2f+2g) component is NON-NEGLIGIBLE!")
            # Fit it too
            coeffs_Fs, _, _, _ = np.linalg.lstsq(A_mat, Fss, rcond=None)
            print(f"  F_s/B1p = {coeffs_Fs[0]:.6f}*rho^3 + {coeffs_Fs[1]:.6f}*rho^4 + {coeffs_Fs[2]:.6f}*rho^5")
