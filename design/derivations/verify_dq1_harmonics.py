"""
Verify delta_q_1 harmonic reduction for BH61 Block 8.

Strategy: Work in semi-closed form keeping (a/r) factors explicit.
Compute delta_q_1 = D(dS_1/dL + dS_1*/dL) symbolically and compare with BH61.

We use SymPy to handle the trig product-to-sum identities automatically.
"""

from sympy import *

# Symbols
f, g, e, eta, theta = symbols('f g e eta theta', real=True)
# a_over_r = (1 + e*cos(f))/eta^2, but we keep it as a symbol for semi-closed form
ar = symbols('ar', positive=True)  # a/r

# Inclination coefficients
B0 = Rational(-1,2) + 3*theta**2/2  # -1/2 + 3theta^2/2
B1p = Rational(3,2) - 3*theta**2/2  # 3/2 - 3theta^2/2 = -B0 + 1

# Common prefactor is mu^2 k_2 / (e L G^3). We'll track everything divided by this.

print("=" * 70)
print("PART 1: Compute alpha * X")
print("=" * 70)

# alpha = 4 + 2*cos(f)/e + 2*a/r
# X = X_A + X_B + X_C + X_D

# X_A = B_0 * sinf * [3*(a/r)*eta^2 - e^2*sin^2(f)]
X_A = B0 * sin(f) * (3*ar*eta**2 - e**2*sin(f)**2)

# X_B = B1' * sinf*(2+e*cosf) * cos(2g+2f)
X_B = B1p * sin(f) * (2 + e*cos(f)) * cos(2*g + 2*f)

# X_C = (B1'/2) * [eta^2 * sin(2g+f) + e*sinf*(2+ecosf)*cos(2g+f)]
X_C = B1p/2 * (eta**2 * sin(2*g + f) + e*sin(f)*(2 + e*cos(f))*cos(2*g + f))

# X_D = (B1'/6) * [eta^2 * sin(2g+3f) + 3*e*sinf*(2+ecosf)*cos(2g+3f)]
X_D = B1p/6 * (eta**2 * sin(2*g + 3*f) + 3*e*sin(f)*(2 + e*cos(f))*cos(2*g + 3*f))

X = X_A + X_B + X_C + X_D

# alpha * X -- but alpha has 1/e term, so multiply through
# alpha = 4 + 2*cos(f)/e + 2*ar
# Total prefactor already has 1/e, so alpha*X/e = (4 + 2cosf/e + 2ar)*X/e
# But actually, the full expression is (mu^2 k_2)/(eL G^3) * [alpha * X + D(X)]
# So we just compute alpha * X as a polynomial in trig functions.

alpha = 4 + 2*cos(f)/e + 2*ar

alpha_X = expand_trig(expand(alpha * X))

print("\nalpha * X (before trig simplification) has", len(Add.make_args(expand(alpha_X))), "terms")

# Use product-to-sum: convert all products of trig functions
alpha_X_simplified = trigsimp(alpha_X)
print("\nalpha * X (after trigsimp) =")
# This might be huge. Let's try a different approach: expand trig products manually.

print("\n" + "=" * 70)
print("PART 2: Use product-to-sum identities on X terms")
print("=" * 70)

# Key identity: sin(A)*cos(B) = [sin(A+B) + sin(A-B)]/2
# sin(A)*sin(B) = [cos(A-B) - cos(A+B)]/2
# cos(A)*cos(B) = [cos(A+B) + cos(A-B)]/2

# Let's expand X_A, X_B, X_C, X_D using these identities

# X_A = B0 * sinf * [3*ar*eta^2 - e^2*sin^2(f)]
#      = B0 * [3*ar*eta^2 * sinf - e^2 * sinf * sin^2(f)]
#      = B0 * [3*ar*eta^2 * sinf - e^2 * sin^3(f)]
# sin^3(f) = (3sinf - sin3f)/4
# So X_A = B0 * [3*ar*eta^2 * sinf - e^2*(3sinf - sin3f)/4]
#         = B0 * [(3*ar*eta^2 - 3e^2/4) sinf + (e^2/4) sin3f]

print("\nX_A expansion:")
print("  X_A = B0 * [(3*ar*eta^2 - 3e^2/4) sinf + (e^2/4) sin(3f)]")

# X_B = B1' * sinf*(2+ecosf) * cos(2g+2f)
# = B1' * [2sinf*cos(2g+2f) + e*sinf*cosf*cos(2g+2f)]
# sinf*cos(2g+2f) = [sin(2g+3f) + sin(-2g-f)]/2 = [sin(2g+3f) - sin(2g+f)]/2
# sinf*cosf*cos(2g+2f) = (sin2f/2)*cos(2g+2f) = [sin(2g+4f) + sin(-2g)]/4 = [sin(2g+4f) - sin(2g)]/4

print("\nX_B expansion:")
print("  sinf*cos(2g+2f) = [sin(2g+3f) - sin(2g+f)]/2")
print("  sinf*cosf*cos(2g+2f) = [sin(2g+4f) - sin(2g)]/4")
print("  X_B = B1' * {sin(2g+3f) - sin(2g+f) + e*[sin(2g+4f) - sin(2g)]/4}")

# X_C = (B1'/2) * [eta^2 * sin(2g+f) + e*sinf*(2+ecosf)*cos(2g+f)]
# e*sinf*(2+ecosf)*cos(2g+f) = e*(2sinf + ecosf sinf)*cos(2g+f)
# = 2e*sinf*cos(2g+f) + e^2*sinf*cosf*cos(2g+f)
# sinf*cos(2g+f) = [sin(2g+2f) + sin(-2g)]/2 = [sin(2g+2f) - sin(2g)]/2
# sinf*cosf*cos(2g+f) = (sin2f/2)*cos(2g+f) = [sin(2g+3f) + sin(-2g+f)]/4 = [sin(2g+3f) - sin(2g-f)]/4

print("\nX_C expansion:")
print("  sinf*cos(2g+f) = [sin(2g+2f) - sin(2g)]/2")
print("  sinf*cosf*cos(2g+f) = [sin(2g+3f) - sin(2g-f)]/4")
print("  e part = 2e*[sin(2g+2f) - sin(2g)]/2 + e^2*[sin(2g+3f) - sin(2g-f)]/4")
print("         = e*sin(2g+2f) - e*sin(2g) + (e^2/4)*sin(2g+3f) - (e^2/4)*sin(2g-f)")
print("  X_C = (B1'/2)*{eta^2*sin(2g+f) + e*sin(2g+2f) - e*sin(2g) + (e^2/4)*sin(2g+3f) - (e^2/4)*sin(2g-f)}")

# X_D = (B1'/6) * [eta^2 * sin(2g+3f) + 3*e*sinf*(2+ecosf)*cos(2g+3f)]
# 3e*sinf*(2+ecosf)*cos(2g+3f) = 3e*(2sinf + ecosf sinf)*cos(2g+3f)
# = 6e*sinf*cos(2g+3f) + 3e^2*sinf*cosf*cos(2g+3f)
# sinf*cos(2g+3f) = [sin(2g+4f) + sin(-2g-2f)]/2 = [sin(2g+4f) - sin(2g+2f)]/2
# sinf*cosf*cos(2g+3f) = (sin2f/2)*cos(2g+3f) = [sin(2g+5f) + sin(-2g-f)]/4 = [sin(2g+5f) - sin(2g+f)]/4

print("\nX_D expansion:")
print("  sinf*cos(2g+3f) = [sin(2g+4f) - sin(2g+2f)]/2")
print("  sinf*cosf*cos(2g+3f) = [sin(2g+5f) - sin(2g+f)]/4")
print("  3e part = 6e*[sin(2g+4f) - sin(2g+2f)]/2 + 3e^2*[sin(2g+5f) - sin(2g+f)]/4")
print("          = 3e*sin(2g+4f) - 3e*sin(2g+2f) + (3e^2/4)*sin(2g+5f) - (3e^2/4)*sin(2g+f)")
print("  X_D = (B1'/6)*{eta^2*sin(2g+3f) + 3e*sin(2g+4f) - 3e*sin(2g+2f) + (3e^2/4)*sin(2g+5f) - (3e^2/4)*sin(2g+f)}")

print("\n" + "=" * 70)
print("PART 3: Collect X by harmonic")
print("=" * 70)

# Let's collect all contributions to each harmonic from X:
# Harmonics present in X: sinf, sin3f, sin(2g), sin(2g-f), sin(2g+f), sin(2g+2f), sin(2g+3f), sin(2g+4f), sin(2g+5f)

print("\nCollecting X by harmonic (factoring out B0 or B1' as appropriate):")

# sin(f): from X_A only: B0*(3*ar*eta^2 - 3e^2/4)
print("\n  sinf: B0*(3*ar*eta^2 - 3e^2/4)")

# sin(3f): from X_A only: B0*(e^2/4)
print("  sin(3f): B0*(e^2/4)")

# sin(2g): from X_B: B1'*(-e/4); from X_C: (B1'/2)*(-e) = -B1'*e/2
# Total: B1'*(-e/4 - e/2) = B1'*(-3e/4)
print("  sin(2g): B1'*(-e/4 - e/2) = B1'*(-3e/4)")

# sin(2g-f): from X_C only: (B1'/2)*(-e^2/4) = -B1'*e^2/8
print("  sin(2g-f): B1'*(-e^2/8)")

# sin(2g+f): from X_B: B1'*(-1); from X_C: (B1'/2)*(eta^2); from X_D: (B1'/6)*(-3e^2/4)
# Total: B1'*[-1 + eta^2/2 - e^2/8]
print("  sin(2g+f): B1'*(-1 + eta^2/2 - e^2/8)")

# sin(2g+2f): from X_C: (B1'/2)*(e) = B1'*e/2; from X_D: (B1'/6)*(-3e) = -B1'*e/2
# Total: B1'*(e/2 - e/2) = 0  !!
print("  sin(2g+2f): B1'*(e/2 - e/2) = 0  [cancels!]")

# sin(2g+3f): from X_B: B1'*(1); from X_C: (B1'/2)*(e^2/4) = B1'*e^2/8; from X_D: (B1'/6)*(eta^2)
# Total: B1'*(1 + e^2/8 + eta^2/6)
print("  sin(2g+3f): B1'*(1 + e^2/8 + eta^2/6)")

# sin(2g+4f): from X_B: B1'*(e/4); from X_D: (B1'/6)*(3e) = B1'*e/2
# Total: B1'*(e/4 + e/2) = B1'*(3e/4)
print("  sin(2g+4f): B1'*(3e/4)")

# sin(2g+5f): from X_D: (B1'/6)*(3e^2/4) = B1'*e^2/8
print("  sin(2g+5f): B1'*(e^2/8)")

print("\n" + "=" * 70)
print("PART 4: Numerical verification of X expansion")
print("=" * 70)

import numpy as np

# Test at specific values
f_val = 1.3  # true anomaly
g_val = 0.7  # argument of perigee
e_val = 0.15
eta_val = np.sqrt(1 - e_val**2)
ar_val = (1 + e_val * np.cos(f_val)) / eta_val**2
theta_val = 0.8  # cos(i)

B0_val = -0.5 + 1.5*theta_val**2
B1p_val = 1.5 - 1.5*theta_val**2

# Direct computation of X
X_A_val = B0_val * np.sin(f_val) * (3*ar_val*eta_val**2 - e_val**2*np.sin(f_val)**2)
X_B_val = B1p_val * np.sin(f_val) * (2 + e_val*np.cos(f_val)) * np.cos(2*g_val + 2*f_val)
X_C_val = B1p_val/2 * (eta_val**2 * np.sin(2*g_val + f_val) + e_val*np.sin(f_val)*(2 + e_val*np.cos(f_val))*np.cos(2*g_val + f_val))
X_D_val = B1p_val/6 * (eta_val**2 * np.sin(2*g_val + 3*f_val) + 3*e_val*np.sin(f_val)*(2 + e_val*np.cos(f_val))*np.cos(2*g_val + 3*f_val))
X_val = X_A_val + X_B_val + X_C_val + X_D_val

# Harmonic expansion
X_harm = (B0_val * (3*ar_val*eta_val**2 - 3*e_val**2/4) * np.sin(f_val)
         + B0_val * (e_val**2/4) * np.sin(3*f_val)
         + B1p_val * (-3*e_val/4) * np.sin(2*g_val)
         + B1p_val * (-e_val**2/8) * np.sin(2*g_val - f_val)
         + B1p_val * (-1 + eta_val**2/2 - e_val**2/8) * np.sin(2*g_val + f_val)
         + 0  # sin(2g+2f) cancels
         + B1p_val * (1 + e_val**2/8 + eta_val**2/6) * np.sin(2*g_val + 3*f_val)
         + B1p_val * (3*e_val/4) * np.sin(2*g_val + 4*f_val)
         + B1p_val * (e_val**2/8) * np.sin(2*g_val + 5*f_val))

print(f"\nX direct = {X_val:.12f}")
print(f"X harmonic = {X_harm:.12f}")
print(f"Difference = {abs(X_val - X_harm):.2e}")

print("\n" + "=" * 70)
print("PART 5: Now compute alpha * X by harmonic")
print("=" * 70)

# alpha = 4 + 2cosf/e + 2*ar
# We need alpha * X where X is expanded in harmonics.
# alpha * sin(A) = 4*sin(A) + (2/e)*cosf*sin(A) + 2*ar*sin(A)
# cosf*sin(A) = [sin(A+f) + sin(A-f)]/2

# So alpha*sin(A) = 4*sin(A) + (1/e)*[sin(A+f) + sin(A-f)] + 2*ar*sin(A)

# For each harmonic sin(nf) or sin(2g+nf), applying alpha:
# alpha * coeff * sin(2g+nf) = coeff * {(4+2ar)*sin(2g+nf) + (1/e)*sin(2g+(n+1)f) + (1/e)*sin(2g+(n-1)f)}

# This shifts harmonics by +/-1 in the f-index. The 2*ar and 4 terms keep the same harmonic.
# The 1/e shifts create new harmonics.

# Let's label the X harmonics by their coefficient (including B0 or B1' factor):
# We use a dict: key = harmonic label, value = coefficient

# First, let's define all X coefficients (divided by B0 or B1' as appropriate)
# B0 harmonics: sinf, sin3f
# B1' harmonics: sin(2g), sin(2g-f), sin(2g+f), sin(2g+3f), sin(2g+4f), sin(2g+5f)

# Actually, let's just track all coefficients with the full B0/B1' factor included.
# For the alpha multiplication, we need to be careful.

# Let me just do this numerically to verify the final answer, then present the algebra.

# Direct computation of alpha * X
alpha_val = 4 + 2*np.cos(f_val)/e_val + 2*ar_val
alpha_X_val = alpha_val * X_val

print(f"\nalpha = {alpha_val:.12f}")
print(f"alpha * X direct = {alpha_X_val:.12f}")

print("\n" + "=" * 70)
print("PART 6: Compute D(X) term by term")
print("=" * 70)

# We need D applied to each sub-expression of X.
# D-action table:
# D(sinf) = sin(2f)/e
# D(cosf) = -(1-cos2f)/e
# D(e) = -2(e+cosf)
# D(eta) = 2*eta*(ar - 1)
# D(eta^2) = 4*eta^2*(ar - 1)
# D(ar) = ... need to compute
# D(a/r) = D[(1+ecosf)/eta^2]
# = [D(e)*cosf + e*D(cosf)]/eta^2 + (1+ecosf)*D(1/eta^2)
# D(1/eta^2) = -2*D(eta)/eta^3 = -2*2*(ar-1)/eta^2 = -4*(ar-1)/eta^2  ... wait
# Actually D(eta^2) = 4*eta^2*(ar-1), so D(1/eta^2) = -D(eta^2)/eta^4 = -4*(ar-1)/eta^2

# D(a/r) = [-2(e+cosf)*cosf - e*(1-cos2f)/e]/eta^2 + (1+ecosf)*(-4*(ar-1))/eta^2
# = [-2(e+cosf)cosf - (1-cos2f)]/eta^2 - 4*(ar-1)*(1+ecosf)/eta^2
# = [-2ecosf - 2cos^2f - 1 + cos2f]/eta^2 - 4*(ar-1)*ar*eta^2/eta^2
# cos2f = 2cos^2f - 1, so:
# = [-2ecosf - 2cos^2f - 1 + 2cos^2f - 1]/eta^2 - 4*ar*(ar-1)
# = [-2ecosf - 2]/eta^2 - 4*ar*(ar-1)
# = -2(1+ecosf)/eta^2 - 4*ar*(ar-1)
# = -2*ar - 4*ar^2 + 4*ar
# = 2*ar - 4*ar^2
# = 2*ar*(1 - 2*ar)

# So D(a/r) = 2*(a/r)*(1 - 2*(a/r))

print("D(a/r) = 2*(a/r)*(1 - 2*a/r)")

# Verify numerically
# D(a/r) should be D[(1+ecosf)/eta^2]
# Use finite differences on the D operator definition
# D = -sum xi_k d/dxi_k = ... it's complicated in Cartesian.
# Better: use the identity D(a/r) = 2*ar*(1-2*ar)
# Check: at f=0, a/r = (1+e)/eta^2 = 1/(1-e). D(a/r) = 2/(1-e) * (1 - 2/(1-e)) = 2/(1-e) * (1-e-2)/(1-e) = 2*(-1-e)/(1-e)^2
# Alternative: D(1/r) directly... let's just verify the formula numerically later.

# For the D(X) computation, we need D acting on each part.
# Rather than symbolically expanding everything, let me compute D(X) numerically
# and compare with the analytic result.

# Numerical D operator: D = -2 sum_k xi_k d/dxi_k  in Delaunay form is
# D = sum_j [-(p_j) d/dL_j + (q_j) d/dl_j]  but this is complicated.
#
# Actually, D acts on functions of (l, g, h, L, G, H) as:
# Df = (Dl) df/dl + (Dg) df/dg + (Dh) df/dh + (DL) df/dL + (DG) df/dG + (DH) df/dH
# where Dl, Dg, Dh, DL, DG, DH are the D-actions on the Delaunay variables.
#
# From BH61:
# DL = -q_1 = -(p_1)'' = ... related to drag function
# Actually, D acts on FUNCTIONS of orbital elements, and we derived:
# D(sinf) = sin(2f)/e, D(cosf) = -(1-cos2f)/e, etc.
# D(e) = -2(e+cosf), D(eta) = 2*eta*(a/r - 1)
# D(a/r) = 2*(a/r)*(1-2*a/r)
# D(sin(2g+nf)) uses D(2g+nf) = 2*Dg + n*Df = -4sinf/e + 2n*sinf/e = 2(n-2)*sinf/e

# So we have all the D-actions we need. Let me compute D(X) symbolically.

print("\n--- Computing D(X_A) ---")
print("X_A = B0 * sinf * [3*ar*eta^2 - e^2*sin^2f]")
print("Let P = sinf, Q = 3*ar*eta^2 - e^2*sin^2f")
print("D(X_A) = B0 * [D(P)*Q + P*D(Q)]")
print()
print("D(sinf) = sin(2f)/e")
print("D(3*ar*eta^2) = 3*[D(ar)*eta^2 + ar*D(eta^2)]")
print("  = 3*[2*ar*(1-2ar)*eta^2 + ar*4*eta^2*(ar-1)]")
print("  = 3*eta^2*ar*[2-4ar + 4ar - 4]")
print("  = 3*eta^2*ar*(-2)")
print("  = -6*ar*eta^2")
print()
print("D(e^2*sin^2f) = 2e*D(e)*sin^2f + e^2*D(sin^2f)")
print("  D(e) = -2(e+cosf)")
print("  D(sin^2f) = 2*sinf*D(sinf) = 2*sinf*sin(2f)/e = 2*sinf*2*sinf*cosf/e = 4*sin^2f*cosf/e")
print("  = 2e*(-2)(e+cosf)*sin^2f + e^2*4*sin^2f*cosf/e")
print("  = -4e(e+cosf)*sin^2f + 4e*sin^2f*cosf")
print("  = -4e*sin^2f*[e+cosf - cosf]")
print("  = -4e^2*sin^2f")
print()
print("So D(Q) = -6*ar*eta^2 - (-4e^2*sin^2f) = -6*ar*eta^2 + 4e^2*sin^2f")
print()
print("D(X_A) = B0 * [(sin2f/e)*(3*ar*eta^2 - e^2*sin^2f) + sinf*(-6*ar*eta^2 + 4*e^2*sin^2f)]")
print("       = B0 * [3*ar*eta^2*sin2f/e - e*sin^2f*sin2f + sinf*(-6*ar*eta^2 + 4*e^2*sin^2f)]")

# Let's verify numerically
def D_ar(ar_v, e_v, f_v, eta_v):
    return 2*ar_v*(1 - 2*ar_v)

def D_eta2(ar_v, eta_v):
    return 4*eta_v**2*(ar_v - 1)

def D_e(e_v, f_v):
    return -2*(e_v + np.cos(f_v))

def D_sinf(f_v, e_v):
    return np.sin(2*f_v)/e_v

def D_cosf(f_v, e_v):
    return -(1 - np.cos(2*f_v))/e_v

def D_sin2f(f_v, e_v):
    # sin(2f): D(2g+2f) would be 0, but sin(2f) is just sin(2f) without g.
    # D(sin(2f)) = cos(2f)*D(2f) = cos(2f)*2*Df = cos(2f)*4sinf/e
    return np.cos(2*f_v)*4*np.sin(f_v)/e_v

# Actually, D(f) = 2sinf/e, so D(2f) = 4sinf/e
# D(sin(nf)) = cos(nf) * n * D(f) = n*cos(nf)*2sinf/e
def D_sin_nf(n, f_v, e_v):
    return n * np.cos(n*f_v) * 2*np.sin(f_v)/e_v

# D(sin(2g+nf)) = cos(2g+nf) * (2*Dg + n*Df) = cos(2g+nf)*2*(n-2)*sinf/e
def D_sin_2g_nf(n, g_v, f_v, e_v):
    return np.cos(2*g_v + n*f_v) * 2*(n-2)*np.sin(f_v)/e_v

def D_cos_2g_nf(n, g_v, f_v, e_v):
    return -np.sin(2*g_v + n*f_v) * 2*(n-2)*np.sin(f_v)/e_v

# Verify D(a/r) numerically using finite differences
# a/r = (1+ecosf)/eta^2
# D(a/r) = [D(1+ecosf)]/eta^2 + (1+ecosf)*D(1/eta^2)
# D(1+ecosf) = D(e)*cosf + e*D(cosf) = -2(e+cosf)*cosf + e*(-(1-cos2f)/e) = -2(e+cosf)*cosf - 1 + cos2f
# = -2ecosf - 2cos^2f - 1 + 2cos^2f - 1 = -2ecosf - 2 = -2(1+ecosf)
# D(1/eta^2) = -D(eta^2)/eta^4 = -4*(ar-1)/eta^2
# D(a/r) = -2(1+ecosf)/eta^2 + (1+ecosf)*(-4*(ar-1))/eta^2
# = (1+ecosf)/eta^2 * [-2 - 4*(ar-1)]
# = ar * [-2 - 4ar + 4]
# = ar * [2 - 4ar]
# = 2*ar*(1-2ar)  CONFIRMED

D_ar_val = 2*ar_val*(1 - 2*ar_val)
print(f"\nD(a/r) = 2*{ar_val:.6f}*(1-2*{ar_val:.6f}) = {D_ar_val:.12f}")

# Now compute D(X) numerically
# D(X_A) = B0 * [D(sinf)*(3ar*eta^2 - e^2*sin^2f) + sinf*D(3ar*eta^2 - e^2*sin^2f)]
D_sinf_val = D_sin_nf(1, f_val, e_val)
D_3ar_eta2 = 3*(D_ar_val*eta_val**2 + ar_val*D_eta2(ar_val, eta_val))
D_e2_sin2f = 2*e_val*D_e(e_val, f_val)*np.sin(f_val)**2 + e_val**2*2*np.sin(f_val)*D_sinf_val
D_Q_A = D_3ar_eta2 - D_e2_sin2f

D_X_A_val = B0_val * (D_sinf_val*(3*ar_val*eta_val**2 - e_val**2*np.sin(f_val)**2) + np.sin(f_val)*D_Q_A)

# D(X_B) = B1' * D[sinf*(2+ecosf)*cos(2g+2f)]
# Since D(2g+2f)=0, D(cos(2g+2f))=0
# D[sinf*(2+ecosf)] = D(sinf)*(2+ecosf) + sinf*D(2+ecosf)
# D(2+ecosf) = D(e)*cosf + e*D(cosf) = -2(e+cosf)cosf - (1-cos2f) = -2ecosf - 2cos^2f - 1 + cos2f
# = -2ecosf - 2cos^2f - 1 + 2cos^2f - 1 = -2ecosf - 2 = -2(1+ecosf)
D_2_plus_ecosf = -2*(1 + e_val*np.cos(f_val))
D_sinf_times_stuff = D_sinf_val*(2+e_val*np.cos(f_val)) + np.sin(f_val)*D_2_plus_ecosf
D_X_B_val = B1p_val * D_sinf_times_stuff * np.cos(2*g_val + 2*f_val)
# Note: cos(2g+2f) passes through D unchanged since D(2g+2f)=0

# D(X_C) = (B1'/2) * D[eta^2*sin(2g+f) + e*sinf*(2+ecosf)*cos(2g+f)]
# D[eta^2*sin(2g+f)] = D(eta^2)*sin(2g+f) + eta^2*D(sin(2g+f))
D_eta2_val = D_eta2(ar_val, eta_val)
D_sin_2gf = D_sin_2g_nf(1, g_val, f_val, e_val)  # n=1: 2*(1-2)*sinf/e = -2sinf/e
D_term_C1 = D_eta2_val*np.sin(2*g_val + f_val) + eta_val**2*D_sin_2gf

# D[e*sinf*(2+ecosf)*cos(2g+f)]
# cos(2g+f) has D(2g+f) = 2Dg + Df = -4sinf/e + 2sinf/e = -2sinf/e
D_cos_2gf = D_cos_2g_nf(1, g_val, f_val, e_val)  # n=1: sin(2g+f)*2*(2-1)*sinf/e = 2sinf*sin(2g+f)/e ...
# Wait, D_cos_2g_nf(n) = -sin(2g+nf)*2*(n-2)*sinf/e
# For n=1: -sin(2g+f)*2*(1-2)*sinf/e = -sin(2g+f)*(-2)*sinf/e = 2*sinf*sin(2g+f)/e
D_cos_2gf_val = -np.sin(2*g_val + f_val)*2*(1-2)*np.sin(f_val)/e_val

# Product rule on e*sinf*(2+ecosf)*cos(2g+f)
factor_esinf_stuff = e_val*np.sin(f_val)*(2+e_val*np.cos(f_val))
D_factor = D_e(e_val, f_val)*np.sin(f_val)*(2+e_val*np.cos(f_val)) + e_val*D_sinf_val*(2+e_val*np.cos(f_val)) + e_val*np.sin(f_val)*D_2_plus_ecosf
D_term_C2 = D_factor*np.cos(2*g_val+f_val) + factor_esinf_stuff*D_cos_2gf_val

D_X_C_val = B1p_val/2 * (D_term_C1 + D_term_C2)

# D(X_D) = (B1'/6) * D[eta^2*sin(2g+3f) + 3e*sinf*(2+ecosf)*cos(2g+3f)]
D_sin_2g3f = D_sin_2g_nf(3, g_val, f_val, e_val)  # n=3: 2*(3-2)*sinf/e = 2sinf/e * cos(2g+3f)
D_term_D1 = D_eta2_val*np.sin(2*g_val + 3*f_val) + eta_val**2*D_sin_2g3f

D_cos_2g3f_val = -np.sin(2*g_val + 3*f_val)*2*(3-2)*np.sin(f_val)/e_val
factor_3esinf_stuff = 3*e_val*np.sin(f_val)*(2+e_val*np.cos(f_val))
D_factor3 = 3*(D_e(e_val, f_val)*np.sin(f_val)*(2+e_val*np.cos(f_val)) + e_val*D_sinf_val*(2+e_val*np.cos(f_val)) + e_val*np.sin(f_val)*D_2_plus_ecosf)
D_term_D2 = D_factor3*np.cos(2*g_val+3*f_val) + factor_3esinf_stuff*D_cos_2g3f_val

D_X_D_val = B1p_val/6 * (D_term_D1 + D_term_D2)

D_X_val = D_X_A_val + D_X_B_val + D_X_C_val + D_X_D_val

print(f"\nD(X_A) = {D_X_A_val:.12f}")
print(f"D(X_B) = {D_X_B_val:.12f}")
print(f"D(X_C) = {D_X_C_val:.12f}")
print(f"D(X_D) = {D_X_D_val:.12f}")
print(f"D(X) = {D_X_val:.12f}")

total_val = alpha_val * X_val + D_X_val
print(f"\nalpha*X + D(X) = {total_val:.12f}")

print("\n" + "=" * 70)
print("PART 7: Compare with BH61 Terms 1-2-3")
print("=" * 70)

# BH61 delta_q_1 = (mu^2 k_2)/(eL G^3) * [alpha*X + D(X)] + (mu^2 k_2)/(eL G^3) * alpha_star*X_star + ...
# But wait: the overall prefactor is (mu^2 k_2)/(eL G^3), and BH61 uses:
# Term 1: mu^2 k_2/(e^2 L^3 G) {...}
# Term 2: mu^2 k_2/(e L^3 G)(a/r) {...}
# Term 3: mu^2 k_2/(e L^3 G) * [incl] {...}

# Our prefactor: mu^2 k_2/(eL G^3) = mu^2 k_2/(eL * L^3 eta^3) = mu^2 k_2/(e L^4 eta^3)
# BH61 Term 1 prefactor: mu^2 k_2/(e^2 L^3 G) = mu^2 k_2/(e^2 L^3 * L*eta) = mu^2 k_2/(e^2 L^4 eta)
# Ratio: our/BH61_T1 = (e^2 L^4 eta)/(e L^4 eta^3) = e/eta^2

# So: our_prefactor * total = BH61_T1_prefactor * (e/eta^2) * total + BH61_T2_prefactor * (...) + BH61_T3_prefactor * (...)

# Actually let me think about this differently.
# L^3 G = L^3 * L * eta = L^4 eta. So:
# mu^2 k_2/(e^2 L^3 G) = mu^2 k_2/(e^2 L^4 eta)
# mu^2 k_2/(eL G^3) = mu^2 k_2/(eL * L^3 eta^3) = mu^2 k_2/(e L^4 eta^3)

# Factor: mu^2 k_2/(e^2 L^3 G) / [mu^2 k_2/(eL G^3)] = eL G^3 / (e^2 L^3 G) = G^2/(e L^2) = L^2 eta^2/(e L^2) = eta^2/e

# So: BH61 Term 1 prefactor = (eta^2/e) * our_prefactor
# BH61 Term 2 prefactor = mu^2 k_2/(eL^3 G)(a/r) = mu^2 k_2 * ar/(eL^4 eta)
# Term 2 / our = ar * e * L^4 eta^3 / (eL^4 eta) = ar * eta^2

# BH61 Term 3 prefactor = mu^2 k_2/(eL^3 G) = mu^2 k_2/(eL^4 eta)
# Term 3 / our = e * L^4 eta^3/(eL^4 eta) = eta^2

# Summary:
# BH61_T1 bracket * (eta^2/e) + BH61_T2 bracket * ar*eta^2 + BH61_T3 bracket * eta^2 * [incl] = alpha*X + D(X)

# Let's compute BH61 side numerically
# Term 1: {(-1+3theta^2)[(ar)^2 eta^2 + ar + 1] sin2f + B1'[(-ar^2 eta^2 - ar + 1)sin2g + (ar^2 eta^2 + ar + 1/3)sin(2g+4f)]}
T1_bracket = ((-1+3*theta_val**2)*((ar_val**2)*eta_val**2 + ar_val + 1)*np.sin(2*f_val)
              + B1p_val*((-ar_val**2*eta_val**2 - ar_val + 1)*np.sin(2*g_val)
                         + (ar_val**2*eta_val**2 + ar_val + Rational(1,3))*np.sin(2*g_val + 4*f_val)))
# Wait, 1/3 should be float
T1_bracket = ((-1+3*theta_val**2)*((ar_val**2)*eta_val**2 + ar_val + 1)*np.sin(2*f_val)
              + B1p_val*((-ar_val**2*eta_val**2 - ar_val + 1)*np.sin(2*g_val)
                         + (ar_val**2*eta_val**2 + ar_val + 1.0/3)*np.sin(2*g_val + 4*f_val)))

# Term 2: {(-1+3theta^2)[(ar)^2 eta^2 + ar + 4] sinf + B1'[(-ar^2 eta^2 - ar + 2)sin(2g+f) + (ar^2 eta^2 + ar + 2)sin(2g+3f)]}
T2_bracket = ((-1+3*theta_val**2)*((ar_val**2)*eta_val**2 + ar_val + 4)*np.sin(f_val)
              + B1p_val*((-ar_val**2*eta_val**2 - ar_val + 2)*np.sin(2*g_val + f_val)
                         + (ar_val**2*eta_val**2 + ar_val + 2)*np.sin(2*g_val + 3*f_val)))

# Term 3: [incl_factor] * [(1-3ar)*e*sin2g + sin(2g+f) - sin(2g-f)]
# incl_factor = 1/4*(1-11theta^2) - 10*theta^4/(1-5*theta^2)
incl3 = 0.25*(1-11*theta_val**2) - 10*theta_val**4/(1-5*theta_val**2)
T3_bracket = (1-3*ar_val)*e_val*np.sin(2*g_val) + np.sin(2*g_val + f_val) - np.sin(2*g_val - f_val)

# BH61 total = T1_bracket*(eta^2/e) + T2_bracket*ar*eta^2 + T3_bracket*eta^2*incl3
BH61_total = T1_bracket*(eta_val**2/e_val) + T2_bracket*ar_val*eta_val**2 + T3_bracket*eta_val**2*incl3

print(f"\nBH61 Term 1 bracket * (eta^2/e) = {T1_bracket*(eta_val**2/e_val):.12f}")
print(f"BH61 Term 2 bracket * ar*eta^2 = {T2_bracket*ar_val*eta_val**2:.12f}")
print(f"BH61 Term 3 bracket * eta^2*incl3 = {T3_bracket*eta_val**2*incl3:.12f}")
print(f"BH61 total = {BH61_total:.12f}")

# Our total is alpha*X + D(X) for S_1 part only.
# The S_1* part contributes to D(dS_1*/dL) which is Term 3 partly.
# Let's compute D(dS_1*/dL) separately.

print("\n" + "=" * 70)
print("PART 8: D(dS_1*/dL) contribution")
print("=" * 70)

# From Part C of the derivation file:
# dS_1*/dL = (mu^2 k_2/G^3) * C_1' * (2eta^2/L) * sin(2g) / (2eta^2/L ... wait
# From the file: dS_1*/dL = (mu^2 k_2 eta^2)/(16 L G^3) * [(1-15theta^2)(1-theta^2)/(1-5theta^2)] sin(2g)
# C_1' = (1-15theta^2)(1-theta^2)/[16(1-5theta^2)]
# So dS_1*/dL = (mu^2 k_2)/(L G^3) * C_1' * 2*eta^2 * sin(2g) ... wait let me re-read

# From the file line 79: dS_1*/dL = (mu^2 k_2 eta^2)/(16 L G^3) * [(1-15theta^2)(1-theta^2)/(1-5theta^2)] sin(2g)
# = (mu^2 k_2)/(L G^3) * C_1' * eta^2 * sin(2g) where C_1' = (1-15theta^2)(1-theta^2)/[16(1-5theta^2)]

# D(dS_1*/dL) = D{(mu^2 k_2)/(L G^3) * C_1' * eta^2 * sin(2g)}
# = (mu^2 k_2) * C_1' * D{eta^2 sin(2g) / (L G^3)}

# D{eta^2 sin(2g)/(L G^3)} = D(1/(L G^3)) * eta^2 sin(2g) + (1/(L G^3)) * D(eta^2 sin(2g))
# D(1/(L G^3)) = (1/(L G^3)) * [(2ar-1) + 3] = (1/(L G^3)) * (2ar + 2) = (1/(L G^3)) * alpha_star
# Actually: D(1/L) = (2ar-1)/L, D(1/G^3) = 3/G^3 (from the file, line 91-92)
# D(1/(LG^3)) = D(1/L)/G^3 + (1/L)*D(1/G^3) = (2ar-1)/(LG^3) + 3/(LG^3) = (2ar+2)/(LG^3)

# D(eta^2 sin(2g)) = D(eta^2)*sin(2g) + eta^2*D(sin(2g))
# D(eta^2) = 4eta^2(ar-1)
# D(sin(2g)) = cos(2g)*2*Dg = cos(2g)*(-4sinf/e)
# D(eta^2 sin(2g)) = 4eta^2(ar-1)*sin(2g) - (4eta^2 sinf/e)*cos(2g)

# Total: D{...} = (2ar+2)/(LG^3) * eta^2*sin(2g) + (1/(LG^3))*[4eta^2(ar-1)*sin(2g) - 4eta^2 sinf/e cos(2g)]
# = (eta^2/(LG^3)) * [(2ar+2)*sin(2g) + 4(ar-1)*sin(2g) - (4sinf/e)*cos(2g)]
# = (eta^2/(LG^3)) * [(6ar-2)*sin(2g) - (4sinf/e)*cos(2g)]

# Now convert to our prefactor mu^2 k_2/(eL G^3):
# D(dS_1*/dL) = (mu^2 k_2) * C_1' * (eta^2/(LG^3)) * [(6ar-2)*sin(2g) - (4sinf/e)*cos(2g)]
# = (mu^2 k_2/(eLG^3)) * C_1' * e*eta^2 * [(6ar-2)*sin(2g) - (4sinf/e)*cos(2g)]

# Use product-to-sum on sinf*cos(2g):
# sinf*cos(2g) = [sin(2g+f) + sin(f-2g)]/2 = [sin(2g+f) - sin(2g-f)]/2

# D(dS_1*/dL) / (mu^2 k_2/(eLG^3)) = C_1' * eta^2 * [e(6ar-2)*sin(2g) - 4*[sin(2g+f) - sin(2g-f)]/2]...
# wait, the (4sinf/e)*cos(2g) term already has 1/e, so:
# e*eta^2 * [(6ar-2)*sin(2g)] - e*eta^2*(4sinf/e)*cos(2g) = ...
# Actually let me redo this. We want everything divided by mu^2 k_2/(eL G^3).

# D(dS_1*/dL) = mu^2 k_2 C_1' eta^2/(LG^3) * [(6ar-2)sin(2g) - (4sinf/e)cos(2g)]

# Divide by mu^2 k_2/(eLG^3): gives C_1' * e * eta^2 * [(6ar-2)sin(2g) - (4sinf/e)cos(2g)]
# = C_1' * eta^2 * [(6ar-2)*e*sin(2g) - 4*sinf*cos(2g)]
# = C_1' * eta^2 * [(6ar-2)*e*sin(2g) - 2*sin(2g+f) + 2*sin(2g-f)]

C1p_val = (1-15*theta_val**2)*(1-theta_val**2)/(16*(1-5*theta_val**2))

DS1star_contrib = C1p_val * eta_val**2 * ((6*ar_val-2)*e_val*np.sin(2*g_val)
                                            - 4*np.sin(f_val)*np.cos(2*g_val))
# = C1p * eta^2 * [(6ar-2)*e*sin(2g) - 2*sin(2g+f) + 2*sin(2g-f)]
DS1star_harm = C1p_val * eta_val**2 * ((6*ar_val-2)*e_val*np.sin(2*g_val)
                                        - 2*np.sin(2*g_val+f_val) + 2*np.sin(2*g_val-f_val))

print(f"D(dS_1*/dL) / prefactor (direct) = {DS1star_contrib:.12f}")
print(f"D(dS_1*/dL) / prefactor (harmonic) = {DS1star_harm:.12f}")

# Full delta_q_1 / prefactor = alpha*X + D(X) + DS1star_contrib
full_our = total_val + DS1star_contrib
print(f"\nFull delta_q_1 / (mu^2 k_2/(eLG^3)) = {full_our:.12f}")
print(f"BH61 total / (mu^2 k_2/(eLG^3)) = {BH61_total:.12f}")
print(f"Difference = {abs(full_our - BH61_total):.2e}")

print("\n" + "=" * 70)
print("PART 9: Additional test points")
print("=" * 70)

def compute_our_delta_q1(f_v, g_v, e_v, theta_v):
    """Compute alpha*X + D(X) + D(dS1*/dL) / prefactor numerically."""
    eta_v = np.sqrt(1 - e_v**2)
    ar_v = (1 + e_v*np.cos(f_v))/eta_v**2
    B0_v = -0.5 + 1.5*theta_v**2
    B1p_v = 1.5 - 1.5*theta_v**2

    # X terms
    X_A_v = B0_v * np.sin(f_v) * (3*ar_v*eta_v**2 - e_v**2*np.sin(f_v)**2)
    X_B_v = B1p_v * np.sin(f_v) * (2 + e_v*np.cos(f_v)) * np.cos(2*g_v + 2*f_v)
    X_C_v = B1p_v/2 * (eta_v**2 * np.sin(2*g_v + f_v) + e_v*np.sin(f_v)*(2 + e_v*np.cos(f_v))*np.cos(2*g_v + f_v))
    X_D_v = B1p_v/6 * (eta_v**2 * np.sin(2*g_v + 3*f_v) + 3*e_v*np.sin(f_v)*(2 + e_v*np.cos(f_v))*np.cos(2*g_v + 3*f_v))
    X_v = X_A_v + X_B_v + X_C_v + X_D_v

    alpha_v = 4 + 2*np.cos(f_v)/e_v + 2*ar_v

    # D(X) - compute numerically using product rule
    Dar_v = 2*ar_v*(1-2*ar_v)
    Deta2_v = 4*eta_v**2*(ar_v-1)
    De_v = -2*(e_v+np.cos(f_v))
    Dsinf_v = np.sin(2*f_v)/e_v
    Dcosf_v = -(1-np.cos(2*f_v))/e_v
    D2ecosf_v = -2*(1+e_v*np.cos(f_v))  # D(2+ecosf) = D(ecosf) = De*cosf + e*Dcosf

    # D(X_A) = B0 * [Dsinf*(3ar*eta^2 - e^2 sin^2f) + sinf*D(3ar*eta^2 - e^2 sin^2f)]
    D_3ar_eta2_v = 3*(Dar_v*eta_v**2 + ar_v*Deta2_v)
    D_e2sin2f_v = 2*e_v*De_v*np.sin(f_v)**2 + e_v**2*2*np.sin(f_v)*Dsinf_v
    DQ_A = D_3ar_eta2_v - D_e2sin2f_v
    DX_A_v = B0_v * (Dsinf_v*(3*ar_v*eta_v**2 - e_v**2*np.sin(f_v)**2) + np.sin(f_v)*DQ_A)

    # D(X_B) = B1' * D[sinf*(2+ecosf)]*cos(2g+2f) + 0 (since D(cos(2g+2f))=0)
    D_sinf_stuff = Dsinf_v*(2+e_v*np.cos(f_v)) + np.sin(f_v)*D2ecosf_v
    DX_B_v = B1p_v * D_sinf_stuff * np.cos(2*g_v+2*f_v)

    # D(X_C)
    Dsin2gf = np.cos(2*g_v+f_v)*2*(1-2)*np.sin(f_v)/e_v  # D(sin(2g+f))
    Dcos2gf = -np.sin(2*g_v+f_v)*2*(1-2)*np.sin(f_v)/e_v  # D(cos(2g+f))
    DtermC1 = Deta2_v*np.sin(2*g_v+f_v) + eta_v**2*Dsin2gf
    efactor = e_v*np.sin(f_v)*(2+e_v*np.cos(f_v))
    Defactor = De_v*np.sin(f_v)*(2+e_v*np.cos(f_v)) + e_v*Dsinf_v*(2+e_v*np.cos(f_v)) + e_v*np.sin(f_v)*D2ecosf_v
    DtermC2 = Defactor*np.cos(2*g_v+f_v) + efactor*Dcos2gf
    DX_C_v = B1p_v/2 * (DtermC1 + DtermC2)

    # D(X_D)
    Dsin2g3f = np.cos(2*g_v+3*f_v)*2*(3-2)*np.sin(f_v)/e_v  # D(sin(2g+3f))
    Dcos2g3f = -np.sin(2*g_v+3*f_v)*2*(3-2)*np.sin(f_v)/e_v  # D(cos(2g+3f))
    DtermD1 = Deta2_v*np.sin(2*g_v+3*f_v) + eta_v**2*Dsin2g3f
    efactor3 = 3*e_v*np.sin(f_v)*(2+e_v*np.cos(f_v))
    Defactor3 = 3*(De_v*np.sin(f_v)*(2+e_v*np.cos(f_v)) + e_v*Dsinf_v*(2+e_v*np.cos(f_v)) + e_v*np.sin(f_v)*D2ecosf_v)
    DtermD2 = Defactor3*np.cos(2*g_v+3*f_v) + efactor3*Dcos2g3f
    DX_D_v = B1p_v/6 * (DtermD1 + DtermD2)

    DX_v = DX_A_v + DX_B_v + DX_C_v + DX_D_v

    # D(dS_1*/dL) contribution
    C1p_v = (1-15*theta_v**2)*(1-theta_v**2)/(16*(1-5*theta_v**2))
    DS1s = C1p_v * eta_v**2 * ((6*ar_v-2)*e_v*np.sin(2*g_v) - 4*np.sin(f_v)*np.cos(2*g_v))

    return alpha_v*X_v + DX_v + DS1s

def compute_BH61_delta_q1(f_v, g_v, e_v, theta_v):
    """Compute BH61 Eq(14) delta_q_1 / (mu^2 k_2/(eLG^3))."""
    eta_v = np.sqrt(1 - e_v**2)
    ar_v = (1 + e_v*np.cos(f_v))/eta_v**2
    B0_v = -0.5 + 1.5*theta_v**2  # = (-1+3theta^2)/2
    B1p_v = 1.5 - 1.5*theta_v**2

    # Term 1: prefactor mu^2 k_2/(e^2 L^3 G) = (eta^2/e) * our_prefactor
    T1 = ((-1+3*theta_v**2)*((ar_v**2)*eta_v**2 + ar_v + 1)*np.sin(2*f_v)
          + B1p_v*((-ar_v**2*eta_v**2 - ar_v + 1)*np.sin(2*g_v)
                   + (ar_v**2*eta_v**2 + ar_v + 1.0/3)*np.sin(2*g_v + 4*f_v)))

    # Term 2: prefactor mu^2 k_2/(eL^3 G)(a/r) = ar*eta^2 * our_prefactor
    T2 = ((-1+3*theta_v**2)*((ar_v**2)*eta_v**2 + ar_v + 4)*np.sin(f_v)
          + B1p_v*((-ar_v**2*eta_v**2 - ar_v + 2)*np.sin(2*g_v + f_v)
                   + (ar_v**2*eta_v**2 + ar_v + 2)*np.sin(2*g_v + 3*f_v)))

    # Term 3: prefactor mu^2 k_2/(eL^3 G) = eta^2 * our_prefactor
    incl3 = 0.25*(1-11*theta_v**2) - 10*theta_v**4/(1-5*theta_v**2)
    T3_brack = (1-3*ar_v)*e_v*np.sin(2*g_v) + np.sin(2*g_v + f_v) - np.sin(2*g_v - f_v)

    return T1*(eta_v**2/e_v) + T2*ar_v*eta_v**2 + T3_brack*eta_v**2*incl3

# Test at multiple points
test_cases = [
    (1.3, 0.7, 0.15, 0.8),
    (0.5, 2.1, 0.05, 0.6),
    (2.5, 1.0, 0.30, 0.3),
    (0.1, 0.3, 0.01, 0.95),
    (3.0, 1.5, 0.20, 0.5),
    (np.pi/3, np.pi/4, 0.1, 0.7),
]

print(f"\n{'f':>6} {'g':>6} {'e':>6} {'theta':>6} | {'Our':>16} {'BH61':>16} {'Diff':>12}")
print("-"*80)
for f_t, g_t, e_t, th_t in test_cases:
    our = compute_our_delta_q1(f_t, g_t, e_t, th_t)
    bh61 = compute_BH61_delta_q1(f_t, g_t, e_t, th_t)
    print(f"{f_t:6.3f} {g_t:6.3f} {e_t:6.3f} {th_t:6.3f} | {our:16.8f} {bh61:16.8f} {abs(our-bh61):12.2e}")

print("\nIf differences are all near machine epsilon, the derivation matches BH61.")
