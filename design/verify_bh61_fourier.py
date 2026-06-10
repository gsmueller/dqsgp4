"""
Verification of Brouwer & Hori (1961) Section I: (p₁) expansion
in powers of e and multiples of l (mean anomaly).

The key identity: x = r/a - 1 = -e cos u, where u is related to l
by Kepler's equation: l = u - e sin u.

We need x^k as Fourier series in l, then combine with the
V·exp(-αr)·p₁ series from Eq (19).
"""

import sympy as sp
from sympy import (symbols, sqrt, cos, sin, Rational, series,
                   expand, simplify, pi, trigsimp, Symbol, I as Im,
                   exp as sp_exp, collect, Add, Mul, Pow, Integer)
import numpy as np
from fractions import Fraction

print("=" * 70)
print("VERIFICATION: Section I - (p₁) Fourier expansion")
print("=" * 70)

# ============================================================
# Kepler expansion: x = r/a - 1 as Fourier series in l
# ============================================================
# x = -e cos u, where u - e sin u = l
# Standard Hansen coefficients / Kepler expansion:
#
# r/a = 1 - e cos u = 1 + Σ_{n=1}^∞ aₙ cos(nl)
# where aₙ = -(2/n) Jₙ'(ne) ... (Bessel function derivatives)
#
# For small e, expand to e⁴:
# cos u = cos l + e(cos 2l - 1)/... Actually let me use the
# classical expansion.
#
# From Kepler's equation, the eccentric anomaly u as function of l and e:
# u = l + e sin l + (e²/2) sin 2l + (e³/8)(3 sin 3l - sin l) + ...
#
# Then x = -e cos u. Let's compute this directly.

e, l = symbols('e l')

# Kepler equation inversion: u = l + Σ aₖ(e) sin(kl)
# u = l + e sin l + e²/2 sin 2l + e³(3 sin 3l - sin l)/8
#     + e⁴(2 sin 2l·cos 2l... actually let me use the standard series.

# Standard Kepler inversion to O(e⁴):
# u = l + e sin l + (e²/2) sin 2l + (e³/8)(3 sin 3l - sin l)
#     + (e⁴/6)(2 sin 4l/4 - sin 2l) ... hmm, let me compute properly.

# Actually, the standard result is:
# u = l + Σ (2/n) Jₙ(ne) sin(nl)
# For small e:
# u = l + e sin l + e²/2 · sin 2l + e³/8 · (3 sin 3l - sin l)
#     + e⁴/192 · (32 sin 4l - 24 sin 2l) + O(e⁵)
# = l + e sin l + e²/2 · sin 2l + e³(3 sin 3l - sin l)/8
#   + e⁴(4 sin 4l - 3 sin 2l)/24 + O(e⁵)

# Wait, let me just compute numerically and use the series properly.
# cos u expanded in l:
# cos u = cos l - e(1 - cos 2l)/2 + ... etc.

# Better approach: expand cos u and sin u in Fourier series of l.
# Standard results (see e.g., Brouwer & Clemence, or any celestial mech text):
#
# cos u = -e/2 + Σ_{k=1}^∞ (2/k) J'_k(ke) cos(kl) ... (Hansen)
#
# For e⁴ accuracy, the standard Fourier-Bessel expansion gives:
# r/a = 1 - e cos u = 1 + e²/2 - e cos l - e²/2 cos 2l
#       - 3e³/8 cos l + 3e³/8 cos 3l + ... (to be careful)
#
# Actually, r/a = 1 - e cos u, and the standard expansion is:
# r/a = 1 + e²/2 + (higher even powers of e)
#       - (1 + ...) e cos l
#       + ... cos 2l + ...
#
# Let me use a direct numerical approach to verify.

def kepler_solve(M, ecc, tol=1e-15):
    """Solve Kepler's equation E - e sin E = M."""
    E = M
    for _ in range(50):
        dE = (M - E + ecc * np.sin(E)) / (1 - ecc * np.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def compute_fourier_coeffs(func_of_l, e_val, n_terms=5, n_points=1024):
    """Compute Fourier coefficients of a function of l numerically."""
    l_vals = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    f_vals = np.array([func_of_l(lv, e_val) for lv in l_vals])

    coeffs = {}
    # a_0 (constant term)
    coeffs[('const', 0)] = np.mean(f_vals)
    for k in range(1, n_terms + 1):
        coeffs[('cos', k)] = 2 * np.mean(f_vals * np.cos(k * l_vals))
        coeffs[('sin', k)] = 2 * np.mean(f_vals * np.sin(k * l_vals))
    return coeffs

def x_of_l(l_val, e_val):
    """x = r/a - 1 as function of l and e."""
    u = kepler_solve(l_val, e_val)
    return -e_val * np.cos(u)  # r/a - 1 = -e cos u

def xn_of_l(l_val, e_val, n):
    """x^n as function of l and e."""
    return x_of_l(l_val, e_val)**n

# ============================================================
# Verify (p₁) = μ·exp(-αa) · Σ of terms with cos(kl)
# ============================================================
# From the paper (Section I), (p₁) has the form:
# (p₁) = μ exp(-αa) {1 + (3/4+αa+α²a²/4)e² + (...)e⁴
#         + [(3+αa)e + (...)e³] cos l
#         + [(15/4+2αa+α²a²/2)e² + (...)e⁴] cos 2l
#         + (...)e³ cos 3l + (...)e⁴ cos 4l}

# The (p₁) expression comes from:
# V·exp(-αr)·p₁ evaluated at doubly-primed variables
# = μ·exp(-αa) · (1-x)^(3/2)/(1+x)^(3/2) · exp(-αax)
# where we need to substitute x = -e cos u and expand in e, l.

# Let me verify the constant term and cos l coefficient numerically.
print("\nNumerical verification of (p₁) Fourier coefficients:")
print("Using exact Kepler equation solution + numerical Fourier analysis\n")

# Test at several e values to extract polynomial dependence
e_test_vals = [0.001, 0.01, 0.05, 0.1]
aa_val = 1.0  # α·a = 1 for testing

def vp1_of_l(l_val, e_val, aa_v=1.0):
    """V·exp(-αr)·p₁ / [μ·exp(-αa)] as function of l, e, αa."""
    xv = x_of_l(l_val, e_val)
    return ((1-xv)/(1+xv))**1.5 * np.exp(-aa_v * xv)

# Compute Fourier coefficients at small e to verify e² terms
print("Testing (p₁) constant term vs e²:")
for e_val in [0.01, 0.05, 0.1, 0.2]:
    coeffs = compute_fourier_coeffs(lambda l, ev: vp1_of_l(l, ev, 1.0),
                                     e_val, n_terms=5)
    const = coeffs[('const', 0)]
    # Paper says const = 1 + (3/4 + αa + α²a²/4)e² + O(e⁴)
    # At αa=1: const = 1 + (3/4+1+1/4)e² = 1 + 2e²
    paper_const = 1 + 2*e_val**2
    # But we computed earlier that the x² constant in (1-x)^(3/2)/(1+x)^(3/2)
    # is 9/2 not 3. This would change the e² constant term.
    # Let me check what our independent computation gives.
    print(f"  e={e_val:.3f}: numerical={const:.10f}, paper(αa=1)={paper_const:.10f}, "
          f"diff={(const-paper_const):.2e}")

# Now let me compute what the CORRECT e² constant should be.
# The (p₁) constant term (averaged over l) at αa=1:
# From the series: μ·exp(-αa) · <(1-x)^(3/2)·(1+x)^(-3/2)·exp(-x)>
# where x = -e cos u and <...> denotes average over l.
#
# x = -e cos u, x² = e² cos² u, <cos² u> = (1+η)/2... hmm.
# Actually, <cos² u> needs careful computation.
# <cos u> = -e/2 (standard), <cos² u> = ½ + ...

print("\nDirect computation of <x^k> (average of x^k over l):")
for k in range(6):
    for e_val in [0.1]:
        avg = compute_fourier_coeffs(lambda l, ev: xn_of_l(l, ev, k),
                                      e_val, n_terms=0)[('const', 0)]
        print(f"  <x^{k}> at e={e_val}: {avg:.10f}")
        if k == 0: print(f"    Expected: 1")
        elif k == 1: print(f"    Expected: e²/2 = {e_val**2/2:.10f}")
        elif k == 2: print(f"    Expected: ~{(3*e_val**4/8 + e_val**2/2):.10f}")

# Actually, <x> = <-e cos u> = -e<cos u> = -e·(-e/2) = e²/2
# <x²> = e²<cos²u> = e²·(½ + ...)
# Standard results: <cos u> = -e/2, <cos²u> = ½(1 + e²/2 + ...)
# <x> = e²/2, <x²> = e²/2 + 3e⁴/8 + ...

print("\nThe (p₁) constant term involves:")
print("  Σ cₖ <x^k> where cₖ are from Eq (19) series")
print("  c₀=1, c₁=-(3+αa), c₂=(9/2+3αa+α²a²/2), ...")
print("  <x⁰>=1, <x¹>=e²/2, <x²>=e²/2+3e⁴/8+O(e⁶)")
print()
print("  Constant to O(e²):")
print("  = c₀ + c₁·<x> + c₂·<x²> + ...")
print("  = 1 + (-(3+αa))·(e²/2) + (9/2+3αa+α²a²/2)·(e²/2) + O(e⁴)")
print("  = 1 + e²/2·[-(3+αa) + (9/2+3αa+α²a²/2)] + O(e⁴)")
print("  = 1 + e²/2·[9/2+3αa+α²a²/2-3-αa] + O(e⁴)")
print("  = 1 + e²/2·[3/2+2αa+α²a²/2] + O(e⁴)")
print("  = 1 + (3/4+αa+α²a²/4)e² + O(e⁴)")
print()
print("  Paper says: 1 + (3/4+αa+α²a²/4)e²")
print("  ✓ This MATCHES despite the x² coefficient being 9/2 (not 3)!")
print("  The reason: only c₁·<x> matters at e² (since <x²> starts at e²).")
print("  Wait, c₂·<x²> also contributes at e² since <x²> has an e² term!")
print()

# More careful:
# <x> = e²/2 + O(e⁴)
# <x²> = e²·<cos²u> where <cos²u> = 1/2 + O(e²)
# So <x²> = e²/2 + O(e⁴)
# Therefore to O(e²):
# const = 1 + c₁·(e²/2) + c₂·(e²/2) + O(e⁴)
# = 1 + (e²/2)[c₁ + c₂]
# = 1 + (e²/2)[-(3+αa) + (c₂_const + 3αa + α²a²/2)]
# If c₂_const = 9/2: [-(3+αa) + 9/2+3αa+α²a²/2] = 3/2+2αa+α²a²/2
# If c₂_const = 3:   [-(3+αa) + 3+3αa+α²a²/2] = 2αa+α²a²/2

# Now paper says constant = 1 + (3/4+αa+α²a²/4)e²
# With c₂=9/2: (e²/2)[3/2+2αa+α²a²/2] = (3/4+αa+α²a²/4)e²  ✓
# With c₂=3:   (e²/2)[2αa+α²a²/2] = (αa+α²a²/4)e²  ✗ (missing 3/4)

print("  CRITICAL CHECK:")
print("  With c₂=9/2 (our computation): constant = 1+(3/4+αa+α²a²/4)e² ✓")
print("  With c₂=3 (paper's Eq 19):     constant = 1+(αa+α²a²/4)e² ✗")
print()
print("  CONCLUSION: The (p₁) constant in Section I is CONSISTENT with c₂=9/2.")
print("  This CONFIRMS that Eq (19) has an error: the x² constant should be 9/2.")
print("  The Section I results (if correct) implicitly used c₂=9/2.")

# Verify numerically
print("\n  Numerical verification at αa=1:")
for e_val in [0.01, 0.05, 0.1]:
    coeffs = compute_fourier_coeffs(lambda l, ev: vp1_of_l(l, ev, 1.0),
                                     e_val, n_terms=5)
    const = coeffs[('const', 0)]
    correct = 1 + (Rational(3,4)+1+Rational(1,4))*e_val**2
    wrong = 1 + (1+Rational(1,4))*e_val**2
    print(f"  e={e_val:.3f}: numerical={const:.10f}")
    print(f"           correct(c₂=9/2)={float(correct):.10f}")
    print(f"           wrong(c₂=3)    ={float(wrong):.10f}")

# ============================================================
# Verify cos l coefficient
# ============================================================
print("\n--- cos l coefficient of (p₁) ---")
print("Paper says: [(3+αa)e + (...)e³] cos l")

for e_val in [0.01, 0.05, 0.1]:
    coeffs = compute_fourier_coeffs(lambda l, ev: vp1_of_l(l, ev, 1.0),
                                     e_val, n_terms=5)
    cos1 = coeffs[('cos', 1)]
    paper_cos1 = (3+1)*e_val  # At αa=1: (3+αa)e = 4e
    print(f"  e={e_val:.3f}: numerical cos l = {cos1:.10f}, paper 4e = {4*e_val:.10f}, "
          f"ratio={cos1/(4*e_val):.6f}")

print("  (Ratio should approach 1 as e→0 if leading term is correct)")

# ============================================================
# Full e²-order verification
# ============================================================
print("\n--- Full verification of (p₁) to e² order ---")
print("Using symbolic Kepler expansion:\n")

# x = -e cos u, and from Kepler: u = l + e sin l + O(e²)
# So cos u = cos(l + e sin l + ...) = cos l - e sin l · sin l + ...
#          = cos l - e sin²l + O(e²) = cos l - e(1-cos 2l)/2 + O(e²)
# x = -e[cos l - e(1-cos 2l)/2] = -e cos l + e²(1-cos 2l)/2

# x² = e² cos²l + O(e³) = e²(1+cos 2l)/2 + O(e³)

# Now (1-x)^(3/2)/(1+x)^(3/2) · exp(-αax):
# Expand using the series coefficients c₀=1, c₁=-(3+αa), c₂=(9/2+3αa+α²a²/2):
# = 1 + c₁x + c₂x² + ...
# = 1 + c₁[-e cos l + e²(1-cos 2l)/2] + c₂[e²(1+cos 2l)/2] + O(e³)
# = 1 - c₁e cos l + (c₁/2)e²(1-cos 2l) + (c₂/2)e²(1+cos 2l) + O(e³)
# = 1 + [(c₁+c₂)/2]e² - c₁e cos l + [(c₂-c₁)/2]e² cos 2l + O(e³)

aa_sym = symbols('aa', positive=True)
c1 = -(3 + aa_sym)
c2_correct = Rational(9,2) + 3*aa_sym + aa_sym**2/2

const_e2 = (c1 + c2_correct)/2
cos1_e1 = -c1
cos2_e2 = (c2_correct - c1)/2

print(f"  Constant to e²: 1 + {expand(const_e2)}·e²")
print(f"  cos l to e¹:    {expand(cos1_e1)}·e")
print(f"  cos 2l to e²:   {expand(cos2_e2)}·e²")
print()

# Paper Section I says:
# const: 1 + (3/4+αa+α²a²/4)e²
# cos l: (3+αa)e
# cos 2l: (15/4+2αa+α²a²/2)e²  ... wait, that's a different coefficient.

# Let me check: (c₂-c₁)/2 = (9/2+3αa+α²a²/2 - (-(3+αa)))/2
# = (9/2+3αa+α²a²/2+3+αa)/2 = (15/2+4αa+α²a²/2)/2 = 15/4+2αa+α²a²/4

print(f"  Paper's constant: 1 + (3/4+αa+α²a²/4)e²")
print(f"  Our constant:     1 + {expand(const_e2)}e²")
check1 = simplify(const_e2 - (Rational(3,4) + aa_sym + aa_sym**2/4))
print(f"  Match: {check1 == 0} (diff = {check1})")

print(f"\n  Paper's cos l: (3+αa)e")
print(f"  Our cos l:     {expand(cos1_e1)}e")
check2 = simplify(cos1_e1 - (3 + aa_sym))
print(f"  Match: {check2 == 0}")

print(f"\n  Paper's cos 2l: (15/4+2αa+α²a²/2)e²")
print(f"  Our cos 2l:     {expand(cos2_e2)}e²")
paper_cos2 = Rational(15,4) + 2*aa_sym + aa_sym**2/2
check3 = simplify(cos2_e2 - paper_cos2)
print(f"  Match: {check3 == 0} (diff = {check3})")

# Hmm, our cos 2l gives (15/4+2αa+α²a²/4), paper says (15/4+2αa+α²a²/2)
# Difference is α²a²/4. Let me recheck.

print(f"\n  Detailed: (c₂-c₁)/2 = ({expand(c2_correct)} - ({expand(c1)}))/2")
diff_c2c1 = expand(c2_correct - c1)
print(f"  c₂-c₁ = {diff_c2c1}")
print(f"  (c₂-c₁)/2 = {expand(diff_c2c1/2)}")

# c₂ - c₁ = (9/2+3αa+α²a²/2) - (-(3+αa)) = 9/2+3αa+α²a²/2+3+αa = 15/2+4αa+α²a²/2
# (c₂-c₁)/2 = 15/4+2αa+α²a²/4

# Paper says 15/4+2αa+α²a²/2. So paper has α²a²/2 but we get α²a²/4.
# BUT WAIT - my e² expansion of cos 2l is oversimplified.
# I only used x = -e cos l + O(e²) and x² = e² cos²l + O(e³)
# But the e² correction to cos u matters!
#
# cos u = cos l - e(1-cos 2l)/2 + O(e²)
# x = -e cos u = -e cos l + e²(1-cos 2l)/2 + O(e³)
#
# This is correct. And c₂·x² gives c₂·e²cos²l/2·(1+cos 2l)/2 at leading order.
#
# Wait, but I also need the c₃ term? No, c₃·x³ is O(e³).
#
# Hmm, but I'm missing the exp(-αax) contribution properly.
# Actually my c₁, c₂ already include the exp factor since
# c_k are coefficients of the COMBINED series from Eq (18)×(ratio^3/2).
# Wait no - I defined c₀, c₁, c₂ as coefficients of (1-x)^(3/2)/(1+x)^(3/2)·exp(-αax),
# which is the full Eq(19) divided by μ·exp(-αa).
#
# So the expansion 1 + c₁x + c₂x² + ... is correct and my algebra should be right.
# But the cos 2l coefficient I computed is 15/4+2αa+α²a²/4 while paper says
# 15/4+2αa+α²a²/2.
#
# The difference α²a²/4 suggests I'm missing a contribution.
# Let me think... The x = -e cos l + e²(1-cos 2l)/2 expansion means:
# x² = e²cos²l - 2e³cos l·(1-cos 2l)/2 + ...
# At e² order: x² = e²cos²l = e²(1+cos 2l)/2
# And x at e² order contributes: c₁·e²(1-cos 2l)/2 to the cos 2l term.
# This gives -c₁/2·e² cos 2l.
#
# So total cos 2l at e²:
# From c₁·x at e²: c₁·e²(-cos 2l)/2 = -(c₁/2)e² cos 2l
# From c₂·x² at e²: c₂·e²(cos 2l)/2 = (c₂/2)e² cos 2l
# = (c₂-c₁)/2·e² cos 2l  ... which is what I had. So no missing term.

print("\n  The α²a² discrepancy in cos 2l coefficient:")
print(f"  Our result: 15/4+2αa+α²a²/4")
print(f"  Paper:      15/4+2αa+α²a²/2")
print(f"  Difference: α²a²/4")
print()
print("  BUT: if the paper used c₂=3 (the wrong value) instead of 9/2,")
c2_paper = 3 + 3*aa_sym + aa_sym**2/2
cos2_paper = expand((c2_paper - c1)/2)
print(f"  then cos 2l = {cos2_paper}")
# (3+3αa+α²a²/2 + 3+αa)/2 = (6+4αa+α²a²/2)/2 = 3+2αa+α²a²/4
print(f"  = 3+2αa+α²a²/4")
print(f"  Paper says 15/4+2αa+α²a²/2. Neither matches!")
print()
print("  This suggests the paper's Section I was computed independently")
print("  of Eq (19), or there are compensating errors.")

# Let me just verify numerically
print("\n  Numerical check of cos 2l coefficient at αa=1:")
for e_val in [0.01, 0.05]:
    coeffs = compute_fourier_coeffs(lambda l, ev: vp1_of_l(l, ev, 1.0),
                                     e_val, n_terms=5)
    cos2 = coeffs[('cos', 2)]
    our_val = (Rational(15,4)+2+Rational(1,4))*e_val**2  # 15/4+2+1/4 = 6
    paper_val = (Rational(15,4)+2+Rational(1,2))*e_val**2  # 15/4+2+1/2 = 6.25
    print(f"  e={e_val}: numerical cos 2l = {cos2:.12f}")
    print(f"            our formula = {float(our_val):.12f}")
    print(f"            paper formula = {float(paper_val):.12f}")
    print(f"            ratio to e²: {cos2/e_val**2:.6f}")
    # At αa=1, the correct value should be computable exactly

print("\n  At αa=1, the cos 2l coefficient / e² should be:")
print("  (from full numerical Fourier) ↑")
print("  Our prediction: 15/4+2+1/4 = 6.0")
print("  Paper prediction: 15/4+2+1/2 = 6.25")
print("  Check which one the numerics support.\n")

print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
The independent derivation reveals:
1. Eq (19) x² constant in V·exp(-αr)·p₁ is 9/2 (not 3)
2. The Section I (p₁) constant term (3/4+αa+α²a²/4)e² is CONSISTENT
   with c₂=9/2, NOT with c₂=3
3. This means Eq (19) in the markdown has a transcription error
4. The cos 2l coefficient needs further investigation

The key discriminant is the numerical Fourier analysis:
the true value confirms which set of coefficients is correct.
""")
