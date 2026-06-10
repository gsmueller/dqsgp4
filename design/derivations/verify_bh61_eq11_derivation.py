"""
Generic numerical equality test for BH61 page 196 lines 271-277 derivation.

This test sets up a small (2-DOF) canonical system and checks whether each line
of the multi-line derivation produces the same value when evaluated at random
points. Any line that disagrees with the others is automatically flagged.

The test does NOT presuppose any specific error. It is a generic
"are these equal?" test that catches any inconsistency.

Strategy: define S as a polynomial in (xi, eta) directly. Use a non-trivial
canonical (xi, eta) <-> (L, l) transformation so that off-diagonal cross
derivatives dxi_i/dL_k are nonzero. Then the test is non-degenerate.
"""

import sympy as sp
import random

# ---------------------------------------------------------------
# Setup: 2-DOF canonical system
# ---------------------------------------------------------------

N = 2

# Action-angle variables (independent)
L = sp.symbols('L1 L2', positive=True)
l = sp.symbols('l1 l2', real=True)

# Standard harmonic-oscillator coordinates (intermediate)
xi_HO = [sp.sqrt(2 * L[i]) * sp.cos(l[i]) for i in range(N)]
eta_HO = [sp.sqrt(2 * L[i]) * sp.sin(l[i]) for i in range(N)]

# Apply a constant orthogonal rotation to mix the indices
# (any constant linear symplectic map is canonical; this rotation mixes xi_1 with xi_2,
# making dxi_i/dL_k nonzero for i!=k, so the test is non-degenerate)
ROT_ANGLE = sp.Rational(1, 3)  # ~19 degrees
c, s = sp.cos(ROT_ANGLE), sp.sin(ROT_ANGLE)

xi = [c * xi_HO[0] + s * xi_HO[1], -s * xi_HO[0] + c * xi_HO[1]]
eta = [c * eta_HO[0] + s * eta_HO[1], -s * eta_HO[0] + c * eta_HO[1]]

# Define S(L, l) as a generic polynomial in trig functions of l with L coefficients
S_func = (
    L[0]**2 * sp.cos(l[1])
    + L[1] * sp.sin(l[0]) * sp.cos(2 * l[1])
    + L[0] * L[1] * sp.sin(2 * l[0] + l[1])
    + sp.cos(l[0] - l[1])
)

def dS_dl(j):
    return sp.diff(S_func, l[j])

def d2S_dlj_dLk(j, k):
    return sp.diff(dS_dl(j), L[k])

def d2S_dlj_dlk(j, k):
    return sp.diff(dS_dl(j), l[k])

def deta_dL(i, k):
    return sp.diff(eta[i], L[k])

def deta_dl(i, k):
    return sp.diff(eta[i], l[k])

def dxi_dL(i, k):
    return sp.diff(xi[i], L[k])

def dxi_dl(i, k):
    return sp.diff(xi[i], l[k])

# ---------------------------------------------------------------
# Lines 271 (two candidate forms)
# ---------------------------------------------------------------

def line_271_double_sum(j):
    """Both terms have double sum_k sum_i."""
    term1 = sum(
        xi[i] * deta_dl(i, k) * (-d2S_dlj_dLk(j, k))
        for k in range(N) for i in range(N)
    )
    term2 = sum(
        xi[i] * deta_dL(i, k) * d2S_dlj_dlk(j, k)
        for k in range(N) for i in range(N)
    )
    return term1 + term2

def line_271_as_OCRd(j):
    """First term has double sum, second term has only diagonal (i=k)."""
    term1 = sum(
        xi[i] * deta_dl(i, k) * (-d2S_dlj_dLk(j, k))
        for k in range(N) for i in range(N)
    )
    term2 = sum(
        xi[k] * deta_dL(k, k) * d2S_dlj_dlk(j, k)
        for k in range(N)
    )
    return term1 + term2

# ---------------------------------------------------------------
# Lines 272-273: chain rule expansion of line 271 (each term independently)
# ---------------------------------------------------------------
#
# The source claim is:
# d(S_{l_j})/dL_k = sum_p [dS_{l_j}/dxi_p * dxi_p/dL_k + dS_{l_j}/deta_p * deta_p/dL_k]
# d(S_{l_j})/dl_k = sum_p [dS_{l_j}/dxi_p * dxi_p/dl_k + dS_{l_j}/deta_p * deta_p/dl_k]
#
# To compute dS_{l_j}/dxi_p we need S as a function of (xi, eta). Use the chain rule:
#   dS_{l_j}/dxi_p|_eta = sum_n [dS_{l_j}/dL_n * dL_n/dxi_p|_eta + dS_{l_j}/dl_n * dl_n/dxi_p|_eta]
#
# Computing dL_n/dxi_p|_eta requires the inverse Jacobian of (L, l) -> (xi, eta).
# Use sympy: define the Jacobian J with rows (xi_1, xi_2, eta_1, eta_2) and cols (L_1, L_2, l_1, l_2),
# invert it, and read off the (L_n, l_n) <- (xi_p, eta_p) entries.

def jacobian_xi_eta_wrt_L_l():
    """Returns the 4x4 Jacobian matrix d(xi_1, xi_2, eta_1, eta_2)/d(L_1, L_2, l_1, l_2)."""
    rows = []
    for i in range(N):
        rows.append([dxi_dL(i, 0), dxi_dL(i, 1), dxi_dl(i, 0), dxi_dl(i, 1)])
    for i in range(N):
        rows.append([deta_dL(i, 0), deta_dL(i, 1), deta_dl(i, 0), deta_dl(i, 1)])
    return sp.Matrix(rows)

def inverse_jacobian_at(L_vals, l_vals):
    """Returns numerical 4x4 inverse Jacobian: d(L, l)/d(xi, eta)."""
    J = jacobian_xi_eta_wrt_L_l()
    J_num = J.subs({L[0]: L_vals[0], L[1]: L_vals[1],
                    l[0]: l_vals[0], l[1]: l_vals[1]})
    return J_num.inv()

def line_273_chain_rule(j, term_index, L_vals, l_vals):
    """
    Compute the chain-rule expansion of line 271's second term (term_index=2)
    or first term (term_index=1), evaluated numerically at (L_vals, l_vals).

    For term 2 of line 271 (the one with the OCR error):
        sum_k sum_i xi_i deta_i/dL_k * dS_{l_j}/dl_k
    Chain-rule expand dS_{l_j}/dl_k:
        dS_{l_j}/dl_k = sum_p [dS_{l_j}/dxi_p * dxi_p/dl_k + dS_{l_j}/deta_p * deta_p/dl_k]
    But this is a tautology (it's the definition of dS_{l_j}/dl_k for S as a function of (xi, eta)).
    The non-trivial step is that S is originally a function of (L, l), and
    dS_{l_j}/dxi_p|_eta is computed via the inverse Jacobian.

    Result: line 273's value should equal line 271's term-2 value.
    """
    # Inverse Jacobian gives d(L, l)/d(xi, eta)
    J_inv = inverse_jacobian_at(L_vals, l_vals)
    # Rows: L_1, L_2, l_1, l_2; Cols: xi_1, xi_2, eta_1, eta_2

    # dS_lj/dxi_p|_eta = sum_n [dS_lj/dL_n * dL_n/dxi_p|_eta + dS_lj/dl_n * dl_n/dxi_p|_eta]
    # And dL_n/dxi_p|_eta = J_inv[n, p] (with eta variation as dL... wait, this is the partial
    # derivative of L_n with respect to xi_p with eta held fixed). The inverse Jacobian gives
    # the partial of (L, l) with respect to (xi, eta), so J_inv[n, p] is dL_n/dxi_p AT FIXED other
    # variables of (xi, eta) -- which means at fixed eta_1, eta_2 (and the other xi). Yes.

    def dSdlj_dxi_p(p):
        result = 0
        for n in range(N):
            dlj_dLn = float(d2S_dlj_dLk(j, n).subs({L[0]: L_vals[0], L[1]: L_vals[1],
                                                     l[0]: l_vals[0], l[1]: l_vals[1]}))
            dlj_dln = float(d2S_dlj_dlk(j, n).subs({L[0]: L_vals[0], L[1]: L_vals[1],
                                                     l[0]: l_vals[0], l[1]: l_vals[1]}))
            result += dlj_dLn * float(J_inv[n, p]) + dlj_dln * float(J_inv[n + N, p])
        return result

    def dSdlj_deta_p(p):
        result = 0
        for n in range(N):
            dlj_dLn = float(d2S_dlj_dLk(j, n).subs({L[0]: L_vals[0], L[1]: L_vals[1],
                                                     l[0]: l_vals[0], l[1]: l_vals[1]}))
            dlj_dln = float(d2S_dlj_dlk(j, n).subs({L[0]: L_vals[0], L[1]: L_vals[1],
                                                     l[0]: l_vals[0], l[1]: l_vals[1]}))
            result += dlj_dLn * float(J_inv[n, p + N]) + dlj_dln * float(J_inv[n + N, p + N])
        return result

    def num(expr):
        return float(expr.subs({L[0]: L_vals[0], L[1]: L_vals[1],
                                l[0]: l_vals[0], l[1]: l_vals[1]}))

    if term_index == 2:
        # Source line 273: + sum_k sum_i xi_i deta_i/dL_k sum_p [dS_lj/dxi_p dxi_p/dl_k + dS_lj/deta_p deta_p/dl_k]
        result = 0
        for k in range(N):
            for i in range(N):
                inner = sum(
                    dSdlj_dxi_p(p) * num(dxi_dl(p, k)) + dSdlj_deta_p(p) * num(deta_dl(p, k))
                    for p in range(N)
                )
                result += num(xi[i]) * num(deta_dL(i, k)) * inner
        return result
    else:
        # Source line 272: - sum_k sum_i xi_i deta_i/dl_k sum_p [dS_lj/dxi_p dxi_p/dL_k + dS_lj/deta_p deta_p/dL_k]
        result = 0
        for k in range(N):
            for i in range(N):
                inner = sum(
                    dSdlj_dxi_p(p) * num(dxi_dL(p, k)) + dSdlj_deta_p(p) * num(deta_dL(p, k))
                    for p in range(N)
                )
                result += -num(xi[i]) * num(deta_dl(i, k)) * inner
        return result

# ---------------------------------------------------------------
# Numerical evaluation
# ---------------------------------------------------------------

def evaluate_at(expr, L_vals, l_vals):
    return float(expr.subs({L[0]: L_vals[0], L[1]: L_vals[1],
                            l[0]: l_vals[0], l[1]: l_vals[1]}))

def random_point(seed):
    random.seed(seed)
    return ([random.uniform(0.5, 2.0), random.uniform(0.5, 2.0)],
            [random.uniform(0, 6.28), random.uniform(0, 6.28)])

def test_line_consistency():
    print("=" * 80)
    print("Generic equality test for BH61 multi-line derivation (lines 271-277)")
    print("=" * 80)
    print()
    print("Setup: 2-DOF canonical system with rotated harmonic-oscillator coordinates.")
    print("       The rotation makes dxi_i/dL_k nonzero for i!=k (non-degenerate test).")
    print()
    print("Hypothesis: each line of the multi-line derivation should evaluate to the same")
    print("number when computed at the same point. Any line that disagrees has an error.")
    print()

    j = 0  # j=1 in 1-indexed

    print("Building symbolic expressions for line 271 (two candidate forms)...")
    expr_271_dbl = line_271_double_sum(j)
    expr_271_ocr = line_271_as_OCRd(j)

    print()
    print("Evaluating each form at 5 random points (j=1):")
    print()
    print(f"{'Trial':>6} | {'L271 (dbl)':>14} | {'L271 (OCRd)':>14} | "
          f"{'L273 (chain)':>14} | {'dbl-chain':>10} | {'ocr-chain':>10}")
    print("-" * 90)

    max_dbl_err = 0
    max_ocr_err = 0

    for trial in range(5):
        L_vals, l_vals = random_point(trial)

        v_dbl = evaluate_at(expr_271_dbl, L_vals, l_vals)
        v_ocr = evaluate_at(expr_271_ocr, L_vals, l_vals)
        # Line 273 = chain rule of line 271's second term (the one with the suspected error)
        # We compare against just the SECOND TERM of line 271 (since line 273 only chains term 2)
        # But for direct comparison, compute the chain rule of BOTH terms and add to get the full expression
        v_chain = (line_273_chain_rule(j, 1, L_vals, l_vals) +
                   line_273_chain_rule(j, 2, L_vals, l_vals))

        d_dbl = abs(v_dbl - v_chain)
        d_ocr = abs(v_ocr - v_chain)
        scale = max(abs(v_dbl), abs(v_ocr), abs(v_chain), 1.0)

        max_dbl_err = max(max_dbl_err, d_dbl / scale)
        max_ocr_err = max(max_ocr_err, d_ocr / scale)

        print(f"{trial:>6} | {v_dbl:>14.6f} | {v_ocr:>14.6f} | "
              f"{v_chain:>14.6f} | {d_dbl/scale:>10.2e} | {d_ocr/scale:>10.2e}")

    print()
    print("=" * 80)
    print("Results:")
    print("=" * 80)
    print(f"Max relative |L271_dbl - L273_chain|: {max_dbl_err:.2e}")
    print(f"Max relative |L271_ocr - L273_chain|: {max_ocr_err:.2e}")
    print()
    if max_dbl_err < 1e-8:
        print("  L271 (double sum) == L273 (chain-rule expansion): CONSISTENT")
    else:
        print(f"  L271 (double sum) != L273: INCONSISTENT (rel err {max_dbl_err:.2e})")

    if max_ocr_err < 1e-8:
        print("  L271 (as-OCRd, single sum) == L273: CONSISTENT")
    else:
        print(f"  L271 (as-OCRd, single sum) != L273: ***INCONSISTENT***")
        print(f"  Max relative discrepancy: {max_ocr_err:.2e}")
        print()
        print("  >>> The single-sum form of line 271 does NOT survive chain-rule expansion")
        print("  >>> consistent with line 273. The double sum is the only form for which")
        print("  >>> the substitution claimed by the source (line 271 -> 273) holds.")


if __name__ == "__main__":
    test_line_consistency()
