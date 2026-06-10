"""
Independent clean-room derivation of M_1(e), M_2(e), M_3(e) and the
linear combination S(e) = 3 M_2 + 3 e M_1 + e M_3.  PURE SYMBOLIC --
no numerical approximations anywhere.

(Used because Octave's symbolic toolbox calls into SymPy under the hood
but is significantly slower per-operation; this script is the raw SymPy
equivalent of independent_X0_octave.m.)

This script does NOT use the framework's claimed closed forms.  It uses
only:
  - Kepler's equation:        l = E - e sin E  =>  dl = (1 - e cos E) dE
  - Orbit equation in E:      r/a = kappa = 1 - e cos E
  - True-anomaly identities:  cos f = (cos E - e)/kappa
                              sin f = (eta sin E)/kappa
  - eta := sqrt(1 - e^2)

Method: residue calculus on the unit circle |z| = 1 via z = e^{iE}.
This is standard complex analysis (Chapter 0.D in our textbook), not
the framework's machinery.  The framework's claimed closed forms for
X_0^{0,m} are NOT used.

Output: closed forms as rational functions of e and eta = sqrt(1 - e^2).
"""

import sympy as sp

e   = sp.Symbol('e',   positive=True)
eta = sp.Symbol('eta', positive=True)
E   = sp.Symbol('E',   real=True)
z   = sp.Symbol('z')

I = sp.I

# Build cos f, sin f, kappa as functions of (cos E, sin E, e, eta).
def build_integrand(j_val):
    """Return kappa * cos(j f) as a polynomial in (cos E, sin E, e, eta)."""
    cE, sE = sp.cos(E), sp.sin(E)
    kappa = 1 - e * cE
    cos_f = (cE - e) / kappa
    sin_f = (eta * sE) / kappa

    # cos(j f) = T_j(cos f), Chebyshev T_j of the first kind.
    # T_0 = 1, T_1 = x, T_{j+1} = 2 x T_j - T_{j-1}.
    T = [sp.S(1), cos_f]
    for k in range(1, j_val):
        T.append(2 * cos_f * T[-1] - T[-2])
    cos_jf = T[j_val]

    integrand = kappa * cos_jf
    integrand = sp.expand(integrand)
    return integrand

# Convert E-integral to z-integral over |z|=1 via z = e^{iE}.
# cos E = (z + 1/z)/2, sin E = (z - 1/z)/(2 i), dE = dz/(i z).
def to_z(integrand_E):
    cE_sub = (z + 1/z) / 2
    sE_sub = (z - 1/z) / (2*I)
    F = integrand_E.subs({sp.cos(E): cE_sub, sp.sin(E): sE_sub})
    F = sp.together(F)
    F = sp.cancel(F)
    return F

# The contour integral is sum of residues at poles inside |z|=1.
# For our integrands (after the substitutions), the relevant poles are:
#   z = 0 (pole of order varying with denominator)
#   z = z_- := (1 - eta)/e  (root of e z^2 - 2 z + e inside |z|=1).
# We sum both.
z_minus = (1 - eta) / e

print("=" * 72)
print("Independent symbolic derivation of M_1, M_2, M_3.")
print("Method: residue calculus on |z| = 1, z = e^{iE}.")
print("PURE SYMBOLIC.  No numerical approximations.")
print("=" * 72)

closed = {}
for j in (1, 2, 3):
    print(f"\n--- M_{j}(e) ---")
    F_E = build_integrand(j)
    print(f"  Integrand (in E)  : {sp.simplify(F_E)}")

    F_z = to_z(F_E) / (I * z)        # multiply by dE/dz factor
    F_z = sp.cancel(F_z)
    print(f"  Integrand (in z)  : prepared (rational function of z)")

    # Sum residues at z = 0 and z = z_-
    res0    = sp.residue(F_z, z, 0)
    res_zm  = sp.residue(F_z, z, z_minus)
    res_sum = sp.simplify(res0 + res_zm)

    # Contour integral = 2 pi i * sum(residues); divide by 2 pi for the average.
    Mj = sp.simplify(I * res_sum)
    Mj = sp.cancel(Mj)
    Mj = sp.simplify(Mj)
    print(f"  M_{j}(e) (raw)    : {Mj}")

    # Reduce eta^2 -> 1 - e^2 only where it would over-simplify;
    # keep eta in the answer.
    Mj_clean = sp.simplify(Mj)
    print(f"  M_{j}(e) (in eta) : {Mj_clean}")
    closed[j] = Mj_clean

print("\n" + "=" * 72)
print("Linear combination  S(e) = 3 M_2 + 3 e M_1 + e M_3")
print("=" * 72)

S = 3 * closed[2] + 3 * e * closed[1] + e * closed[3]
S = sp.together(S)
S = sp.cancel(S)
S = sp.simplify(S)
print(f"S(e) raw : {S}")

# Compare with -M_2 (the proposed identity).
diff = sp.simplify(S - (-closed[2]))
# Reduce via eta^2 = 1 - e^2 to confirm 0 in the quotient ring.
diff_reduced = sp.expand(diff)
for _ in range(8):
    prev = diff_reduced
    diff_reduced = sp.expand(diff_reduced.subs(eta**2, 1 - e**2))
    if diff_reduced == prev:
        break
diff_final = sp.simplify(diff_reduced)
print(f"\nTest: S(e) + M_2(e) reduced mod eta^2 = 1 - e^2:  {diff_final}")

if diff_final == 0:
    print("  CONCLUSION: S(e) = -M_2(e)  IDENTICALLY (after eta^2 = 1 - e^2).")
else:
    print(f"  CONCLUSION: S(e) is NOT equal to -M_2(e); residual = {diff_final}.")

# Also report S(e) in pure-e form, for cross-checking against a downstream
# claim that the cancellation produces a specific multiple of M_2.
S_in_e = sp.expand(S)
for _ in range(8):
    prev = S_in_e
    S_in_e = sp.expand(S_in_e.subs(eta**2, 1 - e**2))
    if S_in_e == prev:
        break
S_in_e = sp.simplify(S_in_e)
print(f"\nS(e) reduced (eta odd powers may remain) : {S_in_e}")
print(f"-M_2(e) (for comparison)                 : {sp.simplify(-closed[2])}")
