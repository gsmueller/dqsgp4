"""
DERIVATION (not fitting) of Brouwer & Hori (1961) Section I: (p₁).

Starting from:
1. Kepler's equation: l = u - e sin u
2. x = r/a - 1 = -e cos u
3. V·exp(-αr)·p₁ / [μ·exp(-αa)] = (1-x)^{3/2}·(1+x)^{-3/2}·exp(-αax)
                                   = Σ c_k x^k

We DERIVE the expansion of c_k·x^k in terms of e and cos(nl).
"""

import sympy as sp
from sympy import (symbols, cos, sin, Rational, series, expand,
                   simplify, collect, trigsimp, sqrt, exp, Symbol,
                   Function, pi, Integer, Poly, degree)

print("=" * 70)
print("ALGEBRAIC DERIVATION OF (p₁) FOURIER EXPANSION")
print("=" * 70)

# ============================================================
# STEP 1: Kepler inversion - derive u(l,e) to O(e⁴)
# ============================================================
print("\n--- STEP 1: Kepler equation inversion ---")

e, l = symbols('e l')

# Kepler: l = u - e sin u
# Solve iteratively: u₀ = l, u₁ = l + e sin u₀ = l + e sin l, ...
# u = l + Σ_{n=1} (1/n!) (d/dl)^{n-1} [sin^n l] · e^n   [Lagrange inversion]

# Order by order:
# u₁ = e sin l
# u₂ = e sin(l + e sin l) - e sin l at O(e²) = e² sin l cos l = (e²/2) sin 2l
# More carefully using Lagrange:
# u = l + e sin l + (e²/2) sin l cos l + (e³/8)(3 sin l cos²l - sin l) + ...

# Let me derive this properly using Lagrange inversion theorem.
# If w = l + e·φ(w), then for any f:
# f(w) = f(l) + Σ_{n=1}^∞ (e^n/n!) · (d/dl)^{n-1} [φ(l)^n · f'(l)]

# Here w = u, φ(u) = sin u, so φ(l) = sin l at zeroth order.
# For f(u) = u: f'(l) = 1
# u = l + Σ (e^n/n!) · (d/dl)^{n-1} [sin^n l]

# n=1: e sin l
# n=2: (e²/2) d/dl[sin²l] = (e²/2)·2 sin l cos l = (e²/2) sin 2l
# But wait, we need d/dl[sin²l] = 2 sin l cos l = sin 2l
# So n=2: (e²/2) sin 2l  ✓

# n=3: (e³/6) (d/dl)²[sin³l]
# sin³l = (3 sin l - sin 3l)/4
# d/dl: (3 cos l - 3 cos 3l)/4
# d²/dl²: (-3 sin l + 9 sin 3l)/4
# So n=3: (e³/6)(-3 sin l + 9 sin 3l)/4 = (e³/24)(-3 sin l + 9 sin 3l)
#        = (e³/8)(-sin l + 3 sin 3l)  ✓

# n=4: (e⁴/24)(d/dl)³[sin⁴l]
# sin⁴l = (3 - 4cos 2l + cos 4l)/8
# d/dl: (8 sin 2l - 4 sin 4l)/8 = sin 2l - sin 4l/2
# d²/dl²: 2 cos 2l - 2 cos 4l
# d³/dl³: -4 sin 2l + 8 sin 4l
# So n=4: (e⁴/24)(-4 sin 2l + 8 sin 4l) = (e⁴/6)(-sin 2l + 2 sin 4l)
#        = (e⁴/3)(-sin 2l/2 + sin 4l)

# Summary:
# u = l + e sin l + (e²/2) sin 2l + (e³/8)(-sin l + 3 sin 3l)
#     + (e⁴/3)(-sin 2l/2 + sin 4l) + O(e⁵)

print("Lagrange inversion of Kepler's equation:")
print("  u = l + e sin l + (e²/2) sin 2l")
print("      + (e³/8)(-sin l + 3 sin 3l)")
print("      + (e⁴/3)(-sin 2l/2 + sin 4l)")
print("      + O(e⁵)")

# ============================================================
# STEP 2: Derive cos u as Fourier series in l
# ============================================================
print("\n--- STEP 2: cos u in terms of e, l ---")

# cos u = cos(l + δu) where δu = e sin l + ...
# = cos l cos(δu) - sin l sin(δu)
# cos(δu) = 1 - δu²/2 + δu⁴/24 - ...
# sin(δu) = δu - δu³/6 + ...

# δu = e sin l + (e²/2) sin 2l + (e³/8)(-sin l + 3 sin 3l)
#      + (e⁴/3)(-sin 2l/2 + sin 4l)

# Let me use sympy for this.
u_sym = symbols('u')
# Define u as series in e
du = e*sin(l) + (e**2/2)*sin(2*l) + (e**3/8)*(-sin(l) + 3*sin(3*l)) + \
     (e**4/3)*(-sin(2*l)/2 + sin(4*l))

# cos(l + du) expanded to O(e⁴)
# Use Taylor: cos(l+du) = cos l - sin l·du - cos l·du²/2 + sin l·du³/6 + cos l·du⁴/24

# du to needed orders
du1 = e*sin(l)  # O(e)
du2 = (e**2/2)*sin(2*l)  # O(e²)
du3 = (e**3/8)*(-sin(l) + 3*sin(3*l))  # O(e³)
du4 = (e**4/3)*(-sin(2*l)/2 + sin(4*l))  # O(e⁴)

# du² to O(e⁴):
# du² = (du1+du2+du3)² = du1² + 2·du1·du2 + du2² + 2·du1·du3 + ...
du_sq_e2 = expand(du1**2)  # e² sin²l
du_sq_e3 = expand(2*du1*du2)  # e³ · 2·sin l·sin 2l/2 = e³ sin l sin 2l
du_sq_e4 = expand(du2**2 + 2*du1*du3)  # O(e⁴)

# du³ to O(e³): du1³ = e³ sin³l
du_cu_e3 = expand(du1**3)  # e³ sin³l

# du⁴ to O(e⁴): du1⁴ = e⁴ sin⁴l
du_qu_e4 = expand(du1**4)  # e⁴ sin⁴l

# cos u = cos l - sin l·(du1+du2+du3+du4)
#         - cos l/2·(du²_e2 + du²_e3 + du²_e4)
#         + sin l/6·du³_e3
#         + cos l/24·du⁴_e4

# Let me collect by powers of e
cos_u_e0 = cos(l)
cos_u_e1 = -sin(l)*du1  # -e sin²l
cos_u_e2 = -sin(l)*du2 - cos(l)/2*du_sq_e2
cos_u_e3 = -sin(l)*du3 - cos(l)/2*du_sq_e3 + sin(l)/6*du_cu_e3
cos_u_e4 = -sin(l)*du4 - cos(l)/2*du_sq_e4 + cos(l)/24*du_qu_e4

# Simplify using trig identities
# e¹: -e sin l sin l = -e sin²l = -e(1-cos 2l)/2
c_e1 = sp.trigsimp(expand(cos_u_e1/e))
print(f"  cos u at e¹: {c_e1} · e")

# e²: -sin l · (e²/2)sin 2l - cos l/2 · e²sin²l
c_e2 = sp.trigsimp(expand(cos_u_e2/e**2))
print(f"  cos u at e²: {c_e2} · e²")

# e³:
c_e3 = sp.trigsimp(expand(cos_u_e3/e**3))
print(f"  cos u at e³: {c_e3} · e³")

# e⁴:
c_e4 = sp.trigsimp(expand(cos_u_e4/e**4))
print(f"  cos u at e⁴: {c_e4} · e⁴")

# Full cos u
print("\n  cos u = cos l + (e¹ terms)·e + (e² terms)·e² + ...")

# x = -e cos u
print("\n--- STEP 3: x = -e cos u expanded ---")
# x_e1 = -e cos l
# x_e2 = -e² · (cos u at e¹) = -e²·c_e1 = e²(1-cos 2l)/2
# x_e3 = -e³ · (cos u at e²) = -e³·c_e2
# x_e4 = -e⁴ · (cos u at e³)

x_e1 = sp.trigsimp(-e*cos(l))
x_e2 = sp.trigsimp(-e**2*c_e1)
x_e3 = sp.trigsimp(-e**3*c_e2)
x_e4 = sp.trigsimp(-e**4*c_e3)

print(f"  x at e¹: {x_e1}")
print(f"  x at e²: {sp.trigsimp(x_e2)}")
print(f"  x at e³: {sp.trigsimp(x_e3)}")
print(f"  x at e⁴: {sp.trigsimp(x_e4)}")

# ============================================================
# STEP 4: Compute powers of x as Fourier series
# ============================================================
print("\n--- STEP 4: x^k as Fourier series in l ---")

# x = -e cos l + (e²/2)(1-cos 2l) + higher
# x = -e cos l + e²/2 - e² cos 2l/2 + O(e³)
# Denote x = x₁ + x₂ + x₃ + x₄ where xₖ = O(eᵏ)
# x₁ = -e cos l
# x₂ = e²/2 - e² cos(2l)/2 = (e²/2)(1-cos 2l)

x1 = -e*cos(l)
x2 = sp.trigsimp(x_e2)
x3 = sp.trigsimp(x_e3)
x4 = sp.trigsimp(x_e4)

print(f"  x₁ = {x1}")
print(f"  x₂ = {x2}")
print(f"  x₃ = {x3}")
print(f"  x₄ = {x4}")

# x² to O(e⁴):
# x² = x₁² + 2x₁x₂ + 2x₁x₃ + x₂²  [keeping terms to e⁴]
x_sq_e2 = sp.trigsimp(expand(x1**2))  # e² cos²l
x_sq_e3 = sp.trigsimp(expand(2*x1*x2))
x_sq_e4 = sp.trigsimp(expand(2*x1*x3 + x2**2))

print(f"\n  x² at e²: {x_sq_e2}")
print(f"  x² at e³: {x_sq_e3}")
print(f"  x² at e⁴: {x_sq_e4}")

# x³ to O(e⁴):
x_cu_e3 = sp.trigsimp(expand(x1**3))
x_cu_e4 = sp.trigsimp(expand(3*x1**2*x2))

print(f"  x³ at e³: {x_cu_e3}")
print(f"  x³ at e⁴: {x_cu_e4}")

# x⁴ to O(e⁴):
x_qu_e4 = sp.trigsimp(expand(x1**4))
print(f"  x⁴ at e⁴: {x_qu_e4}")

# ============================================================
# STEP 5: Combine with the c_k coefficients
# ============================================================
print("\n--- STEP 5: (p₁) = Σ cₖ xᵏ ---")

aa = symbols('aa')

# The correct c_k from Eq (19) derivation:
c = [
    1,
    -(3 + aa),
    Rational(9,2) + 3*aa + aa**2/2,
    -(Rational(11,2) + Rational(9,2)*aa + Rational(3,2)*aa**2 + aa**3/6),
    Rational(51,8) + Rational(11,2)*aa + Rational(9,4)*aa**2 + aa**3/2 + aa**4/24,
]

# (p₁) = c₀ + c₁x + c₂x² + c₃x³ + c₄x⁴ + O(x⁵)
# Expanded in powers of e:

# e⁰ term: c₀ = 1
print("  e⁰ term: 1")

# e¹ term: c₁ · x₁ = c₁ · (-e cos l)
term_e1 = expand(c[1] * (-cos(l)))
print(f"  e¹ cos l: {term_e1}")

# e² terms:
# c₁·x₂ + c₂·(x²)_e2
# x₂ = (e²/2)(1-cos 2l)
# (x²)_e2 = e²cos²l = (e²/2)(1+cos 2l)

# c₁·x₂/e² = c₁·(1-cos 2l)/2
# c₂·(x²)_e2/e² = c₂·(1+cos 2l)/2
# constant part: (c₁+c₂)/2
# cos 2l part: (-c₁+c₂)/2

const_e2 = expand((c[1] + c[2])/2)
cos2l_e2 = expand((-c[1] + c[2])/2)
print(f"\n  e² constant: {const_e2}")
print(f"  e² cos 2l:   {cos2l_e2}")

# Check against paper:
paper_const_e2 = Rational(3,4) + aa + aa**2/4
paper_cos2l_e2 = Rational(15,4) + 2*aa + aa**2/2  # paper value

print(f"\n  Paper e² constant: {paper_const_e2}")
print(f"  Match: {simplify(const_e2 - paper_const_e2) == 0}")

print(f"\n  Paper e² cos 2l:   {paper_cos2l_e2}")
our_cos2l_e2 = expand(cos2l_e2)
print(f"  Our e² cos 2l:     {our_cos2l_e2}")
print(f"  Difference:        {simplify(our_cos2l_e2 - paper_cos2l_e2)}")

# Now do e³ terms
print("\n  e³ terms:")
# c₁·x₃ + c₂·(x²)_e3 + c₃·(x³)_e3
# x₃/e³ = the e³ part of x divided by e³
x3_norm = sp.trigsimp(expand(x3/e**3))
x_sq_e3_norm = sp.trigsimp(expand(x_sq_e3/e**3))
x_cu_e3_norm = sp.trigsimp(expand(x_cu_e3/e**3))

term_e3 = sp.trigsimp(expand(c[1]*x3_norm + c[2]*x_sq_e3_norm + c[3]*x_cu_e3_norm))
print(f"  Full e³ expression: {term_e3}")

# Collect cos l and cos 3l
# term_e3 should have cos l and cos 3l terms
cos1_e3 = term_e3.coeff(cos(l))
cos3_e3 = term_e3.coeff(cos(3*l))
print(f"  e³ cos l:  {cos1_e3}")
print(f"  e³ cos 3l: {cos3_e3}")

paper_cos1_e3 = Rational(3,8) + Rational(3,8)*aa + Rational(7,8)*aa**2 + aa**3/8  # from paper's e³ cos l
# Wait, paper says (3/8+3αa/8+7α²a²/8+α³a³/8) ... let me re-read.
# Paper: [(3+αa)e + (3/8+3αa/8+7α²a²/8+α³a³/8)e³] cos l
# So the e³ coefficient of cos l is (3/8+3αa/8+7α²a²/8+α³a³/8)
# Hmm, that doesn't look right. Let me re-check from the markdown.

# From line 544 of the markdown:
# (p₁) = μ exp(-αa) {1 + (3/4+αa+α²a²/4)e² + (21/64+3αa/64+9α²a²/32+α³a³/8+α⁴a⁴/64)e⁴
# + [(3+αa)e + (3/8+3αa/8+7α²a²/8+α³a³/8)e³] cos l
# ...}

# So paper's e³ cos l coeff = 3/8 + 3αa/8 + 7α²a²/8 + α³a³/8
paper_cos1_e3_v = Rational(3,8) + 3*aa/8 + 7*aa**2/8 + aa**3/8

print(f"\n  Paper e³ cos l:  {paper_cos1_e3_v}")
print(f"  Our e³ cos l:    {expand(cos1_e3)}")
print(f"  Difference:      {simplify(expand(cos1_e3) - paper_cos1_e3_v)}")

# cos 3l
paper_cos3_e3 = Rational(19,4) + 3*aa + 5*aa**2/8 + aa**3/24
print(f"\n  Paper e³ cos 3l: {paper_cos3_e3}")
print(f"  Our e³ cos 3l:   {expand(cos3_e3)}")
print(f"  Difference:      {simplify(expand(cos3_e3) - paper_cos3_e3)}")

# ============================================================
# Now do e⁴ terms
# ============================================================
print("\n  e⁴ terms:")
# c₁·x₄ + c₂·(x²)_e4 + c₃·(x³)_e4 + c₄·(x⁴)_e4
x4_norm = sp.trigsimp(expand(x4/e**4))
x_sq_e4_norm = sp.trigsimp(expand(x_sq_e4/e**4))
x_cu_e4_norm = sp.trigsimp(expand(x_cu_e4/e**4))
x_qu_e4_norm = sp.trigsimp(expand(x_qu_e4/e**4))

term_e4 = sp.trigsimp(expand(
    c[1]*x4_norm + c[2]*x_sq_e4_norm + c[3]*x_cu_e4_norm + c[4]*x_qu_e4_norm
))
print(f"  Full e⁴ expression: {term_e4}")

# Collect terms
const_e4 = term_e4.as_independent(cos(l), cos(2*l), cos(3*l), cos(4*l))[0]
# Hmm, need to be more careful. Let me collect cos harmonics.

# Use a substitution approach: evaluate at specific l values to extract harmonics
# Or use sympy's Fourier tools.

# Actually, let me just collect manually.
# term_e4 should contain: constant, cos 2l, cos 4l terms
# (cos l and cos 3l would be odd-order in the x expansion, but e⁴ is even...)
# Actually no: x has both even and odd parts, so x^k can produce all harmonics.
# But for (p₁), which is an even function of l (since V²·p₁² is even),
# only cos kl terms appear, and e⁴ contributes to constant, cos 2l, cos 4l.

# Let me evaluate term_e4 at specific l values to extract:
import sympy
l_val = symbols('l_val')

# Substitute l→l_val and expand
term_e4_explicit = term_e4.rewrite(cos)
print(f"\n  Rewritten: {sp.trigsimp(term_e4_explicit)}")

# Let me try a direct approach: substitute numerical l values
def eval_at_l(expr, l_num):
    return float(expr.subs([(l, l_num), (aa, 0)]))

# At αa=0, extract harmonics numerically from the symbolic expression
import numpy as np_
l_test = np_.linspace(0, 2*np_.pi, 1024, endpoint=False)

for aa_val_test in [0, 1]:
    f_vals = np_.array([float(term_e4.subs([(l, lv), (aa, aa_val_test)])) for lv in l_test])
    const_num = np_.mean(f_vals)
    cos2_num = 2*np_.mean(f_vals * np_.cos(2*l_test))
    cos4_num = 2*np_.mean(f_vals * np_.cos(4*l_test))
    cos1_num = 2*np_.mean(f_vals * np_.cos(l_test))

    print(f"\n  αa={aa_val_test}: e⁴ harmonics (from exact symbolic expression):")
    print(f"    constant: {const_num:.10f}")
    print(f"    cos l:    {cos1_num:.10f}")
    print(f"    cos 2l:   {cos2_num:.10f}")
    print(f"    cos 4l:   {cos4_num:.10f}")

    if aa_val_test == 0:
        paper_const = 21/64
        paper_cos2 = -1/16
        paper_cos4 = 391/64
        print(f"    Paper constant: {paper_const:.10f}, diff: {const_num-paper_const:.2e}")
        print(f"    Paper cos 2l:   {paper_cos2:.10f}, diff: {cos2_num-paper_cos2:.2e}")
        print(f"    Paper cos 4l:   {paper_cos4:.10f}, diff: {cos4_num-paper_cos4:.2e}")
    elif aa_val_test == 1:
        paper_const = 21/64 + 3/64 + 9/32 + 1/8 + 1/64
        paper_cos2 = -1/16 + 11/12 + 7/8 + 1/4 + 1/48
        paper_cos4 = 391/64 + 101/24 + 35/32 + 1/8 + 1/192
        print(f"    Paper constant: {paper_const:.10f}, diff: {const_num-paper_const:.2e}")
        print(f"    Paper cos 2l:   {paper_cos2:.10f}, diff: {cos2_num-paper_cos2:.2e}")
        print(f"    Paper cos 4l:   {paper_cos4:.10f}, diff: {cos4_num-paper_cos4:.2e}")

print("\n" + "=" * 70)
print("DERIVATION COMPLETE")
print("=" * 70)
