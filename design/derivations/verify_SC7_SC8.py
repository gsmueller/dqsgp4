"""
Verify SC-7: D(X_C) and SC-8: D(X_D) for the delta_q_1 stepwise derivation.

We verify by:
1. Computing D(X_C) symbolically using D-action rules
2. Expanding to standard harmonics sin(2g+nf)
3. Collecting coefficients

D-action rules (all verified):
  D(e) = -2(e + cosf)
  D(eta^2) = 4*eta^2*(a/r - 1)
  D(sinf) = sin2f / e
  D(cosf) = -(1-cos2f) / e
  D(sin(2g+nf)) = 2(n-2)*sinf*cos(2g+nf)/e
  D(cos(2g+nf)) = -2(n-2)*sinf*sin(2g+nf)/e
  D(B_1') = 0

Strategy: work purely in terms of trig functions of f, then use product-to-sum
to reduce everything to sin(2g + kf) harmonics.
"""

from sympy import *

# Symbols
f, g, e = symbols('f g e', real=True)
ar = symbols('ar', positive=True)  # a/r
eta2 = 1 - e**2  # eta^2

# Helper: product-to-sum expansion
def expand_trig_products(expr):
    """Expand all trig products to sums using TR8 and related identities."""
    return TR8(expand_trig(expr))

# ============================================================
# SC-7: D(X_C) where X_C = (B_1'/2)[eta^2 sin(2g+f) + e sinf(2+ecosf) cos(2g+f)]
# ============================================================
# We compute D[...] inside the brackets, since D(B_1'/2) = 0.
# Let F_C = eta^2 sin(2g+f) + e sinf(2+ecosf) cos(2g+f)

print("="*80)
print("SC-7: D(X_C)")
print("="*80)

# Define the arguments
arg_2g_f = 2*g + f
arg_2g_3f = 2*g + 3*f

# ---- Sub-term 1: D[eta^2 sin(2g+f)] ----
# = D(eta^2)*sin(2g+f) + eta^2*D(sin(2g+f))
# D(eta^2) = 4*eta^2*(ar - 1)
# D(sin(2g+f)) = 2*(1-2)*sinf*cos(2g+f)/e = -2*sinf*cos(2g+f)/e

ST1 = 4*eta2*(ar - 1)*sin(arg_2g_f) + eta2*(-2*sin(f)*cos(arg_2g_f)/e)

print("\nSub-term 1: D[eta^2 sin(2g+f)]")
ST1_expanded = expand_trig_products(expand(ST1))
print(f"  = {ST1_expanded}")

# ---- Sub-term 2: D[e * sinf*(2+ecosf) * cos(2g+f)] ----
# Expand sinf*(2+ecosf) = 2sinf + e*sinf*cosf = 2sinf + (e/2)*sin2f
# So sub-term 2 = D[2e*sinf*cos(2g+f) + (e^2/2)*sin2f*cos(2g+f)]

# Part 2a: D[2e*sinf*cos(2g+f)]
# Product of three: A=2e, B=sinf, C=cos(2g+f)
# D(A)*B*C + A*D(B)*C + A*B*D(C)
# D(2e) = -4(e+cosf)
# D(sinf) = sin2f/e
# D(cos(2g+f)) = -2*(1-2)*sinf*sin(2g+f)/e = 2*sinf*sin(2g+f)/e

P2a_term1 = (-4*(e + cos(f))) * sin(f) * cos(arg_2g_f)
P2a_term2 = 2*e * (sin(2*f)/e) * cos(arg_2g_f)
P2a_term3 = 2*e * sin(f) * (2*sin(f)*sin(arg_2g_f)/e)

P2a = P2a_term1 + P2a_term2 + P2a_term3
print("\nPart 2a: D[2e sinf cos(2g+f)]")
P2a_expanded = expand_trig_products(expand(P2a))
print(f"  = {P2a_expanded}")

# Part 2b: D[(e^2/2)*sin2f*cos(2g+f)]
# Product of three: A=e^2/2, B=sin2f, C=cos(2g+f)
# D(A)*B*C + A*D(B)*C + A*B*D(C)
# D(e^2/2) = (1/2)*(-4e(e+cosf)) = -2e(e+cosf)
# D(sin2f) = cos2f*D(2f) = cos2f*(4sinf/e)  [D(2f) = 2*D(f) = 2*(2sinf/e) = 4sinf/e]
# D(cos(2g+f)) = 2sinf sin(2g+f)/e

P2b_term1 = (-2*e*(e + cos(f))) * sin(2*f) * cos(arg_2g_f)
P2b_term2 = (e**2/2) * (cos(2*f) * 4*sin(f)/e) * cos(arg_2g_f)
P2b_term3 = (e**2/2) * sin(2*f) * (2*sin(f)*sin(arg_2g_f)/e)

P2b = P2b_term1 + P2b_term2 + P2b_term3
print("\nPart 2b: D[(e^2/2) sin2f cos(2g+f)]")
P2b_expanded = expand_trig_products(expand(P2b))
print(f"  = {P2b_expanded}")

# Total D[F_C] = ST1 + P2a + P2b
D_FC = ST1 + P2a + P2b
D_FC_expanded = expand_trig_products(expand(D_FC))
print("\n--- Total D[F_C] (expanded to trig sums) ---")
print(f"  = {D_FC_expanded}")

# Now collect coefficients of each harmonic sin(2g+kf)
# We need to express everything in terms of sin(2g), sin(2g+f), sin(2g+2f), etc.
# The ar-dependent terms keep ar explicit.

# Let's extract coefficients by substituting specific values of g
# and solving a linear system. Each harmonic sin(2g+kf) has a specific
# dependence on g. We can do this numerically for a grid of (f, g, e, ar) values.

import numpy as np

def eval_expr(expr, f_val, g_val, e_val, ar_val):
    return float(expr.subs([(f, f_val), (g, g_val), (e, e_val), (ar, ar_val)]))

# Test at multiple parameter values
test_cases = [
    (0.3, 0.7, 0.1, 1.05),
    (1.2, 2.1, 0.3, 1.15),
    (0.8, 0.5, 0.05, 1.02),
    (2.0, 1.0, 0.5, 1.3),
]

print("\n" + "="*80)
print("Numerical verification of SC-7 harmonic decomposition")
print("="*80)

# The harmonics we expect: sin(2g+kf) for k = -2, -1, 0, 1, 2, 3, 4
# Plus the ar-dependent term: 4*eta^2*ar*sin(2g+f)
# We'll verify the full expression matches.

# Let's collect terms symbolically instead.
# Replace sin/cos products with sum forms and collect.

# Actually let's just verify numerically that the hand-derived coefficients are correct.
# From the task description, the expected result structure is:
# D(F_C) = 4*eta^2*ar*sin(2g+f) + C1*sin(2g+f) + C2*sin(2g) + C3*sin(2g+2f) + ...

# Let me define the target harmonics and fit coefficients
from numpy.linalg import lstsq

harmonics_labels = ['sin(2g-2f)', 'sin(2g-f)', 'sin(2g)', 'sin(2g+f)', 'sin(2g+2f)',
                    'sin(2g+3f)', 'sin(2g+4f)', 'ar*sin(2g+f)']

def harmonic_basis(f_v, g_v, e_v, ar_v):
    eta2_v = 1 - e_v**2
    return [
        np.sin(2*g_v - 2*f_v),       # sin(2g-2f)
        np.sin(2*g_v - f_v),          # sin(2g-f)
        np.sin(2*g_v),                # sin(2g)
        np.sin(2*g_v + f_v),          # sin(2g+f)
        np.sin(2*g_v + 2*f_v),        # sin(2g+2f)
        np.sin(2*g_v + 3*f_v),        # sin(2g+3f)
        np.sin(2*g_v + 4*f_v),        # sin(2g+4f)
        ar_v * np.sin(2*g_v + f_v),   # ar * sin(2g+f)
    ]

# Generate many random test points
np.random.seed(42)
N = 200
f_vals = np.random.uniform(0.1, 2.5, N)
g_vals = np.random.uniform(0.1, 3.0, N)
e_vals = np.random.uniform(0.01, 0.5, N)
ar_vals = 1 + e_vals * np.cos(f_vals) / (1 - e_vals**2)  # consistent a/r

# But we need ar to be a free variable since the expression has ar as independent.
# Actually ar IS (1+ecosf)/eta^2, but in our symbolic expression it appears as a separate symbol.
# Let's use consistent values where ar = (1+ecosf)/(1-e^2).

A_matrix = np.zeros((N, len(harmonics_labels)))
b_vector = np.zeros(N)

for i in range(N):
    fv, gv, ev, arv = f_vals[i], g_vals[i], e_vals[i], ar_vals[i]
    A_matrix[i, :] = harmonic_basis(fv, gv, ev, arv)
    b_vector[i] = eval_expr(D_FC_expanded, fv, gv, ev, arv)

coeffs, residuals, rank, sv = lstsq(A_matrix, b_vector, rcond=None)

print("\nFitted harmonic coefficients (should be functions of e and eta^2):")
for label, c in zip(harmonics_labels, coeffs):
    print(f"  {label}: {c:.10f}")

if len(residuals) > 0 and residuals[0] < 1e-10:
    print(f"\n  Residual: {residuals[0]:.2e} -- EXCELLENT FIT")
else:
    print(f"\n  Residual: {residuals[0] if len(residuals)>0 else 'N/A'}")
    print("  WARNING: residual not small, may need more harmonics or ar-dependent terms")

# The above won't work well because coefficients themselves depend on e.
# Let's try a different approach: fix e and ar, vary f and g only.
print("\n" + "="*80)
print("Approach 2: Fix e, fit coefficients as functions of e")
print("="*80)

for e_test in [0.1, 0.3, 0.5]:
    eta2_test = 1 - e_test**2

    # Harmonics without ar-dependence in basis; instead separate ar*sin terms
    # Actually let's just use many (f,g) pairs with ar = (1+ecosf)/eta^2

    N2 = 500
    f2 = np.random.uniform(0.1, 6.0, N2)
    g2 = np.random.uniform(0.1, 6.0, N2)
    ar2 = (1 + e_test*np.cos(f2)) / eta2_test

    # Basis: sin(2g+kf) and ar*sin(2g+kf) for k = -2..4
    basis_labels2 = []
    for k in range(-2, 5):
        basis_labels2.append(f'sin(2g+{k}f)')
    for k in range(-2, 5):
        basis_labels2.append(f'ar*sin(2g+{k}f)')

    A2 = np.zeros((N2, len(basis_labels2)))
    b2 = np.zeros(N2)

    for i in range(N2):
        for j, k in enumerate(range(-2, 5)):
            A2[i, j] = np.sin(2*g2[i] + k*f2[i])
            A2[i, j+7] = ar2[i] * np.sin(2*g2[i] + k*f2[i])
        b2[i] = eval_expr(D_FC_expanded, f2[i], g2[i], e_test, ar2[i])

    c2, res2, _, _ = lstsq(A2, b2, rcond=None)

    print(f"\ne = {e_test}, eta^2 = {eta2_test:.4f}")
    print(f"  Residual: {res2[0] if len(res2)>0 else 'N/A':.2e}")
    for label, cv in zip(basis_labels2, c2):
        if abs(cv) > 1e-10:
            print(f"    {label}: {cv:.10f}")

# ============================================================
# SC-8: D(X_D) where X_D = (B_1'/6)[eta^2 sin(2g+3f) + 3e sinf(2+ecosf) cos(2g+3f)]
# ============================================================
print("\n" + "="*80)
print("SC-8: D(X_D)")
print("="*80)

# D[F_D] where F_D = eta^2 sin(2g+3f) + 3e sinf(2+ecosf) cos(2g+3f)

# Sub-term 1: D[eta^2 sin(2g+3f)]
# D(eta^2)*sin(2g+3f) + eta^2*D(sin(2g+3f))
# = 4*eta^2*(ar-1)*sin(2g+3f) + eta^2*(2*sinf*cos(2g+3f)/e)  [n=3: 2*(3-2)=2]

ST1_D = 4*eta2*(ar - 1)*sin(arg_2g_3f) + eta2*(2*sin(f)*cos(arg_2g_3f)/e)

print("\nSub-term 1: D[eta^2 sin(2g+3f)]")
ST1_D_exp = expand_trig_products(expand(ST1_D))
print(f"  = {ST1_D_exp}")

# Sub-term 2: D[3e sinf(2+ecosf) cos(2g+3f)]
# = D[6e sinf cos(2g+3f) + (3e^2/2) sin2f cos(2g+3f)]

# Part 2a: D[6e sinf cos(2g+3f)]
# A=6e, B=sinf, C=cos(2g+3f)
# D(6e) = -12(e+cosf)
# D(sinf) = sin2f/e
# D(cos(2g+3f)) = -2*(3-2)*sinf*sin(2g+3f)/e = -2*sinf*sin(2g+3f)/e

P2a_D_t1 = (-12*(e + cos(f))) * sin(f) * cos(arg_2g_3f)
P2a_D_t2 = 6*e * (sin(2*f)/e) * cos(arg_2g_3f)
P2a_D_t3 = 6*e * sin(f) * (-2*sin(f)*sin(arg_2g_3f)/e)

P2a_D = P2a_D_t1 + P2a_D_t2 + P2a_D_t3
print("\nPart 2a: D[6e sinf cos(2g+3f)]")
P2a_D_exp = expand_trig_products(expand(P2a_D))
print(f"  = {P2a_D_exp}")

# Part 2b: D[(3e^2/2) sin2f cos(2g+3f)]
# A=3e^2/2, B=sin2f, C=cos(2g+3f)
# D(3e^2/2) = (3/2)*(-4e(e+cosf)) = -6e(e+cosf)
# D(sin2f) = cos2f*4sinf/e
# D(cos(2g+3f)) = -2sinf sin(2g+3f)/e

P2b_D_t1 = (-6*e*(e + cos(f))) * sin(2*f) * cos(arg_2g_3f)
P2b_D_t2 = (3*e**2/2) * (cos(2*f) * 4*sin(f)/e) * cos(arg_2g_3f)
P2b_D_t3 = (3*e**2/2) * sin(2*f) * (-2*sin(f)*sin(arg_2g_3f)/e)

P2b_D = P2b_D_t1 + P2b_D_t2 + P2b_D_t3
print("\nPart 2b: D[(3e^2/2) sin2f cos(2g+3f)]")
P2b_D_exp = expand_trig_products(expand(P2b_D))
print(f"  = {P2b_D_exp}")

# Total D[F_D]
D_FD = ST1_D + P2a_D + P2b_D
D_FD_expanded = expand_trig_products(expand(D_FD))
print("\n--- Total D[F_D] (expanded) ---")
print(f"  = {D_FD_expanded}")

# Numerical fitting for D(X_D) as well
print("\n" + "="*80)
print("Numerical coefficient extraction for D(X_D)")
print("="*80)

for e_test in [0.1, 0.3, 0.5]:
    eta2_test = 1 - e_test**2

    N2 = 500
    f2 = np.random.uniform(0.1, 6.0, N2)
    g2 = np.random.uniform(0.1, 6.0, N2)
    ar2 = (1 + e_test*np.cos(f2)) / eta2_test

    basis_labels2 = []
    for k in range(-2, 7):
        basis_labels2.append(f'sin(2g+{k}f)')
    for k in range(-2, 7):
        basis_labels2.append(f'ar*sin(2g+{k}f)')

    n_basis = len(basis_labels2)
    A2 = np.zeros((N2, n_basis))
    b2 = np.zeros(N2)

    for i in range(N2):
        for j, k in enumerate(range(-2, 7)):
            A2[i, j] = np.sin(2*g2[i] + k*f2[i])
            A2[i, j+9] = ar2[i] * np.sin(2*g2[i] + k*f2[i])
        b2[i] = eval_expr(D_FD_expanded, f2[i], g2[i], e_test, ar2[i])

    c2, res2, _, _ = lstsq(A2, b2, rcond=None)

    print(f"\ne = {e_test}, eta^2 = {eta2_test:.4f}")
    print(f"  Residual: {res2[0] if len(res2)>0 else 'N/A':.2e}")
    for label, cv in zip(basis_labels2, c2):
        if abs(cv) > 1e-10:
            print(f"    {label}: {cv:.10f}")


# ============================================================
# Direct symbolic coefficient extraction
# ============================================================
print("\n" + "="*80)
print("Direct symbolic extraction via orthogonality")
print("="*80)

# For a function F(f, g) that is a sum of sin(2g+kf) harmonics,
# we can extract the coefficient of sin(2g+kf) by:
# C_k = (2/pi^2) * integral from 0 to 2pi dg integral from 0 to 2pi df F(f,g) sin(2g+kf)
# But this is expensive symbolically.

# Instead, let's just verify the hand computation numerically at specific e values.
# We'll compute D(X_C)/(B_1'/2) at many (f,g) points and verify it matches
# the claimed harmonic decomposition.

print("\nVerifying SC-7 hand-derived coefficients:")
print("D(X_C)/(B_1'/2) = 4*eta^2*(a/r)*sin(2g+f)")
print("  + (-4*eta^2 + 2)*sin(2g+f)")
print("  + (eta^2/e + 2e + e + 5e/4)*sin(2g)")
print("  ... etc.")

# Actually, let me just verify the full symbolic result by plugging in.
# The cleanest approach: compute D_FC symbolically, then for each harmonic,
# integrate out the harmonic and check.

# Let's use a cleaner method: evaluate at specific f values with g chosen
# to isolate harmonics.

# For extracting coefficient of sin(2g+kf), set g such that 2g+kf has a known value.
# E.g., set 2g = pi/2 - kf, then sin(2g+kf) = 1.
# Then evaluate: all other sin(2g+mf) = sin(pi/2 + (m-k)f) = cos((m-k)f)
# This doesn't isolate perfectly. Better to use the least-squares approach above.

# Let me verify by computing the symbolic integrals for specific harmonics.
# Coefficient of sin(2g+kf) in F(f,g):
# C_k(f) such that F = sum_k C_k(f) sin(2g+kf)
# This means C_k(f) = (1/pi) * int_0^{2pi} F(f,g) sin(2g+kf) dg

print("\nExtracting harmonic coefficients symbolically via integration over g:")

for expr, name in [(D_FC, "D(F_C)"), (D_FD, "D(F_D)")]:
    print(f"\n--- {name} ---")
    for k in range(-2, 7):
        # Coefficient of sin(2g+kf): integrate expr*sin(2g+kf) dg from 0 to 2pi, divide by pi
        integrand = expr * sin(2*g + k*f)
        coeff_integral = integrate(integrand, (g, 0, 2*pi))
        coeff_k = simplify(coeff_integral / pi)
        if coeff_k != 0:
            # Factor out common terms and simplify
            coeff_k = collect(expand(coeff_k), [ar, e, sin(f), cos(f)])
            print(f"  sin(2g+{k}f): {coeff_k}")
