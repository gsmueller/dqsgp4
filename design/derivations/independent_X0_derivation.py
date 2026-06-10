"""
Independent clean-room derivation of three orbital integrals plus a linear
combination.

For 0 < e < 1 and the Keplerian relations
    l = E - e sin E             (Kepler's equation)
    r/a = 1 - e cos E =: kappa  (radius / semi-major axis)
    cos f = (cos E - e)/kappa   (true anomaly via eccentric anomaly)
    sin f = (eta sin E)/kappa,   eta := sqrt(1 - e^2)
we want
    M_j(e) := (1/(2 pi)) * integral_0^{2 pi} cos(j f) dl,    j = 1, 2, 3.

Differentiating Kepler's equation gives dl/dE = 1 - e cos E = kappa, so
    M_j(e) = (1/(2 pi)) * integral_0^{2 pi} kappa * cos(j f(E)) dE.

Method: residue calculus on the unit circle.

Let z = e^{iE}. Then
    cos E = (z + 1/z)/2,    sin E = (z - 1/z)/(2 i),    dE = dz / (i z).
The integral over [0, 2 pi] becomes a contour integral on |z| = 1, computed
by summing residues at poles inside the unit disk.

For j = 1:
    kappa * cos f = (1 - e cos E) * (cos E - e)/(1 - e cos E) = cos E - e,
so M_1(e) = -e directly (no contour integration needed).

For j = 2, 3: cos(j f) is built via the binomial expansion of
Re((cos f + i sin f)^j). Substituting cos f, sin f, kappa in terms of z gives
a rational function of z whose contour integral we evaluate by residues.
The only finite pole inside |z| = 1 is the root of (1 - e cos E) inside the
disk, namely
    1 - e (z + 1/z)/2 = 0   <=>   e z^2 - 2 z + e = 0
    => z_+ = (1 + eta)/e   (outside disk),    z_- = (1 - eta)/e   (inside disk),
where eta = sqrt(1 - e^2). We may also have a pole at z = 0; both are handled
by SymPy's residue() routine on the rational integrand in z.

Verification:
  (a) Series in e at small e for the closed form vs an independent
      term-by-term series integration in E (no closed form used).
  (b) High-precision mpmath quadrature of the original l-integral at
      e = 0.3, 0.7 to ~15 sig figs.
"""

import sympy as sp
from sympy import (
    symbols, sqrt, cos, sin, pi, simplify, series, Rational, S,
    Poly, expand, factor, lambdify, I, residue, together,
)
from mpmath import mp, mpf, quad, cos as mpcos, sin as mpsin, sqrt as mpsqrt


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
e, E, z = symbols('e E z', positive=True)
eta_sym = symbols('eta', positive=True)
eta = sqrt(1 - e**2)
kappa = 1 - e * cos(E)

cos_f = (cos(E) - e) / kappa
sin_f = (eta * sin(E)) / kappa


def cos_jf(jval):
    """Build cos(j f) as the real part of (cos f + i sin f)^j."""
    expr = S.Zero
    for k in range(jval + 1):
        coeff = sp.binomial(jval, k)
        # (i)^k contributes (-1)^(k/2) when k even, else imaginary -> drop
        if k % 4 == 0:
            sign = 1
        elif k % 4 == 2:
            sign = -1
        else:
            continue
        expr += sign * coeff * cos_f**(jval - k) * sin_f**k
    return expr


def integrand_in_E(jval):
    """kappa * cos(j f), as a function of E and e."""
    return sp.simplify(kappa * cos_jf(jval))


# ---------------------------------------------------------------------------
# Compute M_j by residues
# ---------------------------------------------------------------------------
def compute_Mj_residues(jval):
    """
    Convert the E-integral over [0, 2 pi] to a contour integral on |z|=1
    via z = e^{iE}, then evaluate by SymPy's residue() routine.
    """
    F = integrand_in_E(jval)
    # Substitute cos E -> (z + 1/z)/2, sin E -> (z - 1/z)/(2 i).
    F_z = F.subs({cos(E): (z + 1/z)/2, sin(E): (z - 1/z)/(2*I)})
    # dE = dz / (i z); integrand for the contour integral is F_z * dz/(i z).
    G = sp.together(F_z / (I * z))
    # G is a rational function in z (with parameter e).
    G = sp.cancel(G)

    # Identify the poles inside |z| = 1. The denominator factors as a
    # power of z and a power of (e z^2 - 2 z + e).  Roots of e z^2 - 2 z + e
    # are z_- = (1 - eta)/e (inside) and z_+ = (1 + eta)/e (outside).
    # Pole at z = 0 may also occur for j >= 2.
    z_minus = (1 - eta) / e

    # Sum of residues at z = 0 and z = z_minus
    res_sum = residue(G, z, 0) + residue(G, z, z_minus)

    # Contour integral = 2 pi i * sum(residues), then divide by (2 pi).
    Mj = sp.simplify(I * res_sum)
    Mj = sp.cancel(Mj)
    Mj = sp.simplify(Mj)
    return Mj


def to_eta_form(expr):
    """Express in terms of eta = sqrt(1-e^2) for compactness."""
    out = expr
    out = out.subs(sqrt(1 - e**2), eta_sym)
    out = out.subs((1 - e**2)**Rational(1, 2), eta_sym)
    out = sp.simplify(out)
    # also try eta^2 -> 1-e^2 collapsed both ways
    return out


# ---------------------------------------------------------------------------
# Compute M_1, M_2, M_3
# ---------------------------------------------------------------------------
print("=" * 72)
print("Computing M_1, M_2, M_3 via residue calculus on |z| = 1")
print("=" * 72)

# For M_1, the integrand reduces algebraically before we ever need residues:
F1 = sp.simplify(integrand_in_E(1))
print(f"Integrand kappa*cos(f) reduces to: {F1}")
M1_direct = sp.integrate(F1, (E, 0, 2*pi)) / (2*pi)
M1 = sp.simplify(M1_direct)
print(f"M_1(e) = {M1}    (direct, no residues needed)\n")

print("Computing M_2 by residues...")
M2 = compute_Mj_residues(2)
print(f"M_2(e) = {M2}\n")

print("Computing M_3 by residues...")
M3 = compute_Mj_residues(3)
print(f"M_3(e) = {M3}\n")

# Also express in terms of eta
M1_eta = to_eta_form(M1)
M2_eta = to_eta_form(M2)
M3_eta = to_eta_form(M3)

print("In terms of eta = sqrt(1 - e^2):")
print(f"  M_1 = {M1_eta}")
print(f"  M_2 = {M2_eta}")
print(f"  M_3 = {M3_eta}")
print()


# ---------------------------------------------------------------------------
# S(e) = 3 M_2 + 3 e M_1 + e M_3
# ---------------------------------------------------------------------------
print("=" * 72)
print("S(e) = 3 M_2 + 3 e M_1 + e M_3")
print("=" * 72)

S_expr = 3 * M2 + 3 * e * M1 + e * M3
S_simp = sp.simplify(sp.together(S_expr))
S_fact = sp.factor(S_simp)
S_eta  = to_eta_form(S_simp)

print(f"S(e)        = {S_simp}")
print(f"S(e) factored = {S_fact}")
print(f"S(e) in eta = {S_eta}")
print()

print("Proportionality tests:")
for label, Mj in [("M_1", M1), ("M_2", M2), ("M_3", M3)]:
    ratio = sp.simplify(S_simp / Mj)
    print(f"  S / {label} = {ratio}")
print()


# ---------------------------------------------------------------------------
# Verification (a): independent low-order series in e
# ---------------------------------------------------------------------------
print("=" * 72)
print("Verification (a): series in e at order 6, both routes")
print("=" * 72)

print("Series of closed forms:")
for label, expr in [("M_1", M1), ("M_2", M2), ("M_3", M3), ("S",  S_simp)]:
    s = sp.series(expr, e, 0, 7).removeO()
    print(f"  {label}(e) = {sp.expand(s)}")
print()

print("Independent series from term-by-term E-integration of integrand:")
order = 6
for jval in [1, 2, 3]:
    F = integrand_in_E(jval)
    F_series = sp.series(F, e, 0, order + 1).removeO()
    F_series = sp.expand(F_series)
    Mj_series = sp.integrate(F_series, (E, 0, 2*pi)) / (2*pi)
    Mj_series = sp.expand(sp.simplify(Mj_series))
    print(f"  M_{jval}(e) [series] = {Mj_series}")
print()


# ---------------------------------------------------------------------------
# Verification (b): high-precision quadrature at e = 0.3, 0.7
# ---------------------------------------------------------------------------
print("=" * 72)
print("Verification (b): mpmath quadrature vs closed form")
print("=" * 72)

mp.dps = 30

def Mj_numeric(jval, eval_):
    """Numerical (1/(2pi)) integral over l in [0,2pi] of cos(j f(l))."""
    e_v = mpf(eval_)
    eta_v = mpsqrt(1 - e_v * e_v)

    def kepler_solve(l):
        Ek = mpf(l)
        for _ in range(100):
            f_val = Ek - e_v * mpsin(Ek) - l
            fp = 1 - e_v * mpcos(Ek)
            dE = f_val / fp
            Ek = Ek - dE
            if abs(dE) < mpf('1e-28'):
                break
        return Ek

    def integrand(l):
        Ek = kepler_solve(l)
        cE = mpcos(Ek)
        sE = mpsin(Ek)
        kp = 1 - e_v * cE
        cF = (cE - e_v) / kp
        # cos(j f) by Chebyshev recursion in cF
        if jval == 1:
            return cF
        if jval == 2:
            return 2 * cF * cF - 1
        if jval == 3:
            return 4 * cF**3 - 3 * cF
        raise ValueError("jval not supported")

    return quad(integrand, [0, 2 * mp.pi]) / (2 * mp.pi)


def closed_value(expr, eval_):
    f = lambdify(e, expr, modules='mpmath')
    return f(mpf(eval_))


for eval_ in [mpf('0.3'), mpf('0.7')]:
    print(f"\n--- e = {eval_} ---")
    for label, Mj in [("M_1", M1), ("M_2", M2), ("M_3", M3)]:
        cf = closed_value(Mj, eval_)
        nq = Mj_numeric(int(label.split('_')[1]), eval_)
        print(f"  {label}: closed = {mp.nstr(cf, 18)}")
        print(f"        numeric = {mp.nstr(nq, 18)}")
        print(f"        diff    = {mp.nstr(cf - nq, 6)}")
    cf_S = closed_value(S_simp, eval_)
    nq_S = (3 * Mj_numeric(2, eval_)
            + 3 * eval_ * Mj_numeric(1, eval_)
            + eval_ * Mj_numeric(3, eval_))
    print(f"  S  : closed = {mp.nstr(cf_S, 18)}")
    print(f"        numeric = {mp.nstr(nq_S, 18)}")
    print(f"        diff    = {mp.nstr(cf_S - nq_S, 6)}")

print("\n" + "=" * 72)
print("Done.")
print("=" * 72)
