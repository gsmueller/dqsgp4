"""
Pre-computation and symbolic verification of the second partials of
S_1 needed by Chapter 12 (delta_q_j = D(partial S_1 / partial L_j)).

This script de-risks Chapter 12 by establishing -- AHEAD OF TIME --
which second partials are independent and which follow from GSI.

The D operator action on partial S_1 / partial X for X in {L, G, H}:
  D(partial S_1 / partial X)
    = (DL) partial^2 S_1 / (partial L partial X)
      + (DG) partial^2 S_1 / (partial G partial X)
      + (DH) partial^2 S_1 / (partial H partial X)
      + (Dl) partial^2 S_1 / (partial l partial X)
      + (Dg) partial^2 S_1 / (partial g partial X)
      + (Dh) partial^2 S_1 / (partial h partial X)

With Dh = 0 (from BH61 Eq. 10), the last term vanishes for any F.
So we need partial^2 S_1 / (partial X_1 partial X_2) for X_1 in {L, G, H}
and X_2 in {L, G, H, l, g} -- that's 15 second partials.

GSI reduces these: partial S_1 / partial G = -(3/G) S_1 - theta (partial S_1 / partial H)
                                             - (1/eta) (partial S_1 / partial L)

Applied to the second partials:
  partial^2 S_1 / (partial L partial G) = partial/partial L [RHS of GSI]
  partial^2 S_1 / (partial G partial G) = partial/partial G [RHS of GSI]
  partial^2 S_1 / (partial G partial H) = partial/partial H [RHS of GSI]
  partial^2 S_1 / (partial G partial l) = partial/partial l [RHS of GSI]
  partial^2 S_1 / (partial G partial g) = partial/partial g [RHS of GSI]

So of the 15, the 5 involving partial G reduce to GSI-derivative form.
The other 10 (involving only L and H on one side) are independent but
many are zero (e.g., partial^2 S_1 / partial H^2 is computable from ch09c).

This script: abstractly verifies the GSI-reduction of second partials,
using SymPy's Function-of-variables to treat S_1 as a black box.
"""

import sympy as sp

# Set up
L, G, H, l, g = sp.symbols('L G H l g', positive=True)
alpha = sp.symbols('alpha', real=True)

theta = H / G
e = sp.sqrt(1 - G**2 / L**2)
eta = sp.sqrt(1 - e**2)

F_fn = sp.Function('F')
F = F_fn(theta, e, l, g)
S = G**(-alpha) * F

print("=" * 72)
print("Pre-computation of second partials via GSI (Ch 12 de-risk)")
print("=" * 72)
print()

# Compute S and its first partials
dS_dL = sp.diff(S, L)
dS_dG_raw = sp.diff(S, G)
dS_dH = sp.diff(S, H)
dS_dl = sp.diff(S, l)
dS_dg = sp.diff(S, g)

# GSI-form of dS/dG (should equal dS_dG_raw)
dS_dG_GSI = -(alpha / G) * S - theta * dS_dH - (1 / eta) * dS_dL

# Simplify
def sym_simplify(expr):
    """Helper: substitute eta = G/L to kill the sqrt(1 - (1 - G^2/L^2)) = G/L
    form, then simplify."""
    e1 = expr.rewrite(sp.sqrt).subs(sp.sqrt(1 - (1 - G**2/L**2)), G/L)
    return sp.simplify(e1)

diff_GSI = sym_simplify(dS_dG_raw - dS_dG_GSI)
print(f"Check 0: GSI holds for first partial?")
print(f"  dS/dG_raw - dS/dG_GSI (simplified) = {diff_GSI}")
assert diff_GSI == 0, "GSI first-partial check FAILED"
print(f"  PASS")
print()

# Now verify GSI-reduction of the 5 second partials involving G.
# Each should equal  partial/partial X [RHS of GSI_first_partial].

checks = [
    ("d^2 S / dL dG", L),
    ("d^2 S / dG^2",  G),
    ("d^2 S / dG dH", H),
    ("d^2 S / dG dl", l),
    ("d^2 S / dG dg", g),
]

print("Check 1-5: GSI reduction of second partials involving G")
for name, var in checks:
    # Direct second partial
    direct = sp.diff(dS_dG_raw, var)
    # GSI-derived: differentiate the GSI form with respect to var
    via_GSI = sp.diff(dS_dG_GSI, var)
    # Difference
    diff = sym_simplify(direct - via_GSI)
    print(f"  {name}:  direct - via_GSI = {diff}")
    assert diff == 0, f"Second-partial GSI reduction FAILED for {name}"

print(f"  All 5 PASS")
print()

# Summary: for Chapter 12, the 15 second partials reduce to 10 "independent"
# (involving only L, H on the X_1 side) plus 5 "GSI-derived" (involving G).
# Of the 10 independent, many are symmetric by mixed-partial equality:
#   d^2/dLdH = d^2/dHdL, etc.
# So truly independent ones are:
#   d^2/dL^2,  d^2/dL dH,  d^2/dL dl,  d^2/dL dg
#   d^2/dH^2,  d^2/dH dl,  d^2/dH dg
#   d^2/dl^2,  d^2/dl dg
#   d^2/dg^2  <- may be trivial if sin(kg) only appears at low orders
# Total ~10 independent second partials, each computable directly.

# For S_1 specifically:
# - d^2/dH^2 = 0 beyond trivial (H enters only through theta = H/G; d theta/dH = 1/G).
# - d^2/dl dg, d^2/dg^2 come from S_1's sin(jf+2g) and sin(2g) harmonics.
# - d^2/dL^2, d^2/dLdH, etc. use ch09a's Lemma E.1.1 (prefactor L-independent).

print("=" * 72)
print("GSI reduction of second partials verified for abstract F.")
print("Chapter 12 needs only ~10 independent second partials of S_1;")
print("the 5 involving partial G follow algebraically from GSI + first partials.")
print()
print("Full enumeration of required independent second partials:")
print("  d^2 S_1 / dL^2          <- from ch09a's Prop E.2 by differentiation")
print("  d^2 S_1 / dL dH         <- (Prop E.2 w.r.t. H)")
print("  d^2 S_1 / dL dl         <- (Prop E.2 w.r.t. l)")
print("  d^2 S_1 / dL dg         <- (Prop E.2 w.r.t. g)")
print("  d^2 S_1 / dH^2          <- from ch09c by differentiation w.r.t. H (likely zero)")
print("  d^2 S_1 / dH dl         <- from ch09c w.r.t. l")
print("  d^2 S_1 / dH dg         <- from ch09c w.r.t. g")
print("  d^2 S_1 / dl^2          <- from ch09e by differentiation w.r.t. l")
print("  d^2 S_1 / dl dg         <- from ch09e w.r.t. g")
print("  d^2 S_1 / dg^2          <- from ch09e w.r.t. g")
print()
print("These 10 will be explicitly computed in ch12 setup, NOT Chapter 9.")
print("=" * 72)
