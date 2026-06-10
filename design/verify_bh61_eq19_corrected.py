"""
Re-verify Eq (19) after user's corrections.
Confirm all 6 x-coefficients of V*exp(-alpha*r)*p_1 match the independent derivation.
Also check q_1 e^4 coefficient.
"""

import sympy as sp
from sympy import symbols, sqrt, Rational, series, expand, simplify, exp

x, aa = symbols('x aa')

print("=" * 70)
print("RE-VERIFICATION OF CORRECTED EQ (19)")
print("=" * 70)

# Independent derivation: V*exp(-αr)*p₁ / [μ*exp(-αa)]
# = (1-x)^(3/2) / (1+x)^(3/2) * exp(-αa*x)
expr = ((1-x)/(1+x))**Rational(3,2) * exp(-aa*x)
ser = series(expr, x, 0, n=6).removeO()

# User's corrected markdown values:
corrected = [
    1,
    -(3 + aa),
    Rational(9,2) + 3*aa + aa**2/2,
    -(Rational(11,2) + Rational(9,2)*aa + Rational(3,2)*aa**2 + aa**3/6),
    Rational(51,8) + Rational(11,2)*aa + Rational(9,4)*aa**2 + aa**3/2 + aa**4/24,
    -(Rational(57,8) + Rational(51,8)*aa + Rational(11,4)*aa**2 + Rational(3,4)*aa**3 + aa**4/8 + aa**5/120)
]

print("\nV*exp(-αr)*p₁ coefficients:")
all_ok = True
for k in range(6):
    computed = expand(ser.coeff(x, k))
    paper = expand(corrected[k])
    diff = simplify(computed - paper)
    ok = diff == 0
    if not ok: all_ok = False
    status = "✓" if ok else f"✗ diff={diff}"
    print(f"  x^{k}: {status}")
    if not ok:
        print(f"    computed: {computed}")
        print(f"    paper:    {paper}")

print(f"\nEq (19) p₁: {'ALL VERIFIED' if all_ok else 'ERRORS REMAIN'}")

# Now check q₁
print("\n" + "=" * 70)
print("Eq (19) q₁ factor verification")
print("=" * 70)

# q₁ = (2sinf/(eη)) * (1 + e²x)
# The x⁰ factor is 2/(eη) * sinf
# 1/η = (1-e²)^(-1/2) = 1 + e²/2 + 3e⁴/8 + 5e⁶/16 + ...
e = symbols('e', positive=True)
inv_eta = series(1/sqrt(1 - e**2), e, 0, n=6)
print(f"\n1/η = {inv_eta}")
print(f"  = 1 + e²/2 + 3e⁴/8 + ...")
print(f"\nPaper's q₁ x⁰ factor: (2/e)(1 + e²/2 + e⁴/8) sinf")
print(f"Correct q₁ x⁰ factor: (2/e)(1 + e²/2 + 3e⁴/8) sinf")
print(f"\nThe e⁴ coefficient should be 3/8, not 1/8.")
print(f"This is another original paper error in Eq (19).")

# But wait - maybe the paper already accounts for this elsewhere?
# The q₁ expression writes the FULL V*exp(-αr)*q₁, not just q₁.
# V*exp(-αr)*q₁ = V*exp(-αr) * [2e sin u + (2η/e) sin f]
# Using sin u = r*sinf/(aη) = (1+x)*sinf/η:
# q₁ = 2e*(1+x)*sinf/η + (2η/e)*sinf = (2sinf/(eη)) * [e²(1+x) + η²]
# = (2sinf/(eη)) * [1 + e²x]   since e² + η² = 1

# So V*exp(-αr)*q₁ = (μ/a)^½ * exp(-αa) * Eq(18 series) * (2sinf/(eη)) * (1+e²x)
# The paper factors this as:
# = (μ/a)^½ * exp(-αa) * {Eq(18)} * [(2/e)(1+e²/2+e⁴/8)sinf + 2e(1+e²/2)x sinf]

# The (2/e)(1+e²/2+e⁴/8) is supposed to be (2/(eη)) at x=0
# But 2/(eη) = (2/e)(1+e²/2+3e⁴/8+...), not (2/e)(1+e²/2+e⁴/8)
# The (1+e²x)/η factor when expanded gives:
# (1/η) + (e²x/η) = (1+e²/2+3e⁴/8+...) + (1+e²/2+...)e²x + ...
# At x⁰: coeff = 1+e²/2+3e⁴/8
# At x¹: coeff = e²(1+e²/2+...) = e²+e⁴/2+...

# The paper says x¹ coeff is 2e(1+e²/2) = e²·2/e·(1+e²/2)
# Our derivation: e²/η = e²(1+e²/2+3e⁴/8+...)
# At e⁴ order, the x¹ e⁴ term would be e⁴/2, and the paper says e²·e²/2 = e⁴/2. OK.
# But the x¹ factor in the paper is written as 2e(1+e²/2), which is:
# 2e + e³ = e times (2+e²) ... the 2e(1+e²/2) = 2e+e³
# Our x¹ factor: (2/(eη))·e² = 2e/η = 2e(1+e²/2+3e⁴/8+...)
# = 2e+e³+3e⁵/4+...
# At e⁴ in the product with sin f, the x¹·e⁴ contribution is 3e⁵/4·sinf which is O(e⁵)
# So truncating at e⁴·sinf, paper's x¹ factor 2e(1+e²/2) is correct through e³.

print(f"\nSo the q₁ error is ONLY in the e⁴ term of the x⁰ coefficient:")
print(f"  Paper: 1/8")
print(f"  Correct: 3/8")
print(f"  This affects the e⁴ order of the Fourier expansion of (q₁)")
