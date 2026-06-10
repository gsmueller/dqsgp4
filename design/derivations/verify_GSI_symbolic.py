"""
Symbolic verification of the General Structural Identity (GSI) in
ABSTRACT form.

Theorem (GSI).  Let S(L, G, H, l, g) = G^{-alpha} * F(theta, e, l, g)
where theta = H/G,  e = sqrt(1 - G^2/L^2),  and F is smooth.  Then:

    dS/dG|_{L,H,l,g}  =  -(alpha/G) S  -  theta * dS/dH|_{L,G,l,g}
                                      -  (1/eta) * dS/dL|_{G,H,l,g}

where eta = sqrt(1 - e^2) = G/L.

This script verifies the theorem WITHOUT specializing F.  Treating F as
an abstract SymPy Function of (theta, e, l, g), we compute both sides
of GSI using SymPy's chain-rule machinery and show the difference is 0.

This is load-bearing for Chapters 9b, 10, 11, 12 of the BH61 derivation.
"""

import sympy as sp

# Positive reals
L, G, H = sp.symbols('L G H', positive=True)
l, g = sp.symbols('l g', real=True)
alpha = sp.symbols('alpha', real=True)

# Delaunay-derived geometric variables
theta = H / G
e = sp.sqrt(1 - G**2 / L**2)
eta = sp.sqrt(1 - e**2)           # equivalent to G/L; we let SymPy deduce this

# Abstract function F(theta, e, l, g) -- no specific form assumed.
F = sp.Function('F')
# S = G^{-alpha} * F(theta(G,H), e(L,G), l, g)
S = G**(-alpha) * F(theta, e, l, g)

print("=" * 72)
print("Symbolic verification of the General Structural Identity (GSI)")
print("in abstract form.")
print("=" * 72)
print()
print(f"S = G^(-alpha) * F(theta, e, l, g)")
print(f"  theta = H/G")
print(f"  e = sqrt(1 - G^2/L^2)")
print(f"  eta = sqrt(1 - e^2)")
print()

# Chain-rule partials.  SymPy handles the chain rule through the substituted
# theta(G,H) and e(L,G) automatically.
dS_dL = sp.diff(S, L)
dS_dG = sp.diff(S, G)
dS_dH = sp.diff(S, H)

# Predicted RHS of GSI
RHS = -(alpha / G) * S - theta * dS_dH - (1 / eta) * dS_dL

# Difference: should be identically 0
diff = sp.simplify(dS_dG - RHS)

# Simplify using eta = G/L (valid because eta^2 = 1 - e^2 = G^2/L^2 and G, L > 0).
# This step is necessary because SymPy may carry sqrt(1 - (1 - G^2/L^2)) unsimplified.
diff_substituted = diff.rewrite(sp.sqrt).subs(sp.sqrt(1 - (1 - G**2/L**2)), G/L)
diff_substituted = sp.simplify(diff_substituted)

print("Computing dS/dG - [-(alpha/G) S - theta dS/dH - (1/eta) dS/dL]...")
print()
print(f"Raw difference:")
print(f"  {diff}")
print()
print(f"After substituting eta = G/L and simplifying:")
print(f"  {diff_substituted}")
print()

if diff_substituted == 0:
    print("*** PASS: GSI verified symbolically for abstract F and alpha. ***")
    print()
    print("This establishes that GSI holds for ANY smooth F(theta, e, l, g)")
    print("and ANY real alpha.  Load-bearing theorem for Chs 9b, 10, 11, 12.")
else:
    print("*** FAIL: GSI did not reduce to 0. ***")
    print("   Further investigation required.")

print("=" * 72)
