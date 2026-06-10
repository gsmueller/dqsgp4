"""
CRITICAL TEST: Does our verified Poisson bracket {p1, S1} reproduce BH61 Eq (14) δp₁?

BH61 Eq (14) gives:
  δp₁ = 3(μ²k₂/L³) { (-1/2 + 3θ²/2)[-η⁻³ + (a/r)³(1 - 2a/r)]
                     + (3/2 - 3θ²/2)(a/r)³(1 - 2a/r)cos(2g+2f) }

Our Poisson bracket gives {p₁, S₁} computed from first principles with all
partial derivatives through Kepler's equation.

If they match, we can derive BH61 from math foundations.
If they don't, either BH61 has an error in Eq (14) or we're computing
something different.

NOTE: BH61's δp₁ = D(∂S₁/∂l₁'') where D = -Σ ξ_j ∂/∂ξ_j.
The Poisson bracket {p₁, S₁} should equal this because at first order
in the Lie transform, the correction to any function F is {F, S₁}.
For F = p₁ = L(2a/r - 1), this gives the first-order oblateness correction.

The k₂ factor: BH61 includes k₂ explicitly. Our S₁ was derived with k₂=1
(the perturbation parameter factored out). So {p₁, S₁} = δp₁/k₂.
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

MU = 1.0

def solve_kepler(M, e, tol=1e-15):
    E = M + e * sin(M)
    for _ in range(100):
        dE = (M - E + e * sin(E)) / (1 - e * cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def true_from_eccentric(E, e):
    return 2 * arctan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))

def delaunay_to_orbital(L, G, H, l, g, h, mu=MU):
    a = L**2 / mu
    eta = G / L
    e = sqrt(1 - eta**2)
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    return a, e, eta, theta, f, r, E


def bh61_eq14_dp1(L, G, H, l, g, h, mu=MU):
    """
    BH61 Eq (14) δp₁, with k₂ = 1 (factored out).

    δp₁/k₂ = 3(μ²/L³) { (-1/2 + 3θ²/2)[-η⁻³ + (a/r)³(1 - 2a/r)]
                        + (3/2 - 3θ²/2)(a/r)³(1 - 2a/r)cos(2g+2f) }
    """
    a, e, eta, theta, f, r, E = delaunay_to_orbital(L, G, H, l, g, h, mu)

    B0 = -0.5 + 1.5 * theta**2
    B1 = 1.5 - 1.5 * theta**2  # = 3/2 sin²I

    ar = a / r  # a/r ratio
    ar3 = ar**3  # (a/r)³

    prefactor = 3.0 * mu**2 / L**3

    term1 = B0 * (-1.0 / eta**3 + ar3 * (1.0 - 2.0 * ar))
    term2 = B1 * ar3 * (1.0 - 2.0 * ar) * cos(2*g + 2*f)

    return prefactor * (term1 + term2)


def poisson_bracket_numerical(L, G, H, l, g, h, eps=1e-7, mu=MU):
    """Compute {p1, S1} by numerical finite differences."""
    def p1(L, G, H, l, g, h):
        a = L**2 / mu
        E = solve_kepler(l, sqrt(1 - (G/L)**2))
        r = a * (1 - sqrt(1-(G/L)**2) * cos(E))  # wrong
        # Actually: r = a*(1 - e*cos(E))
        e = sqrt(1 - (G/L)**2)
        E = solve_kepler(l, e)
        r = a * (1 - e * cos(E))
        return L * (2*a/r - 1)

    def S1(L, G, H, l, g, h):
        a = L**2 / mu
        eta = G / L
        e = sqrt(1 - eta**2)
        theta = H / G
        E = solve_kepler(l, e)
        f = true_from_eccentric(E, e)
        Gamma = mu**2 / (a**3 * eta**3)
        B0 = -0.5 + 1.5 * theta**2
        s2 = 1 - theta**2
        phi = f - l
        return -Gamma * (
            B0 * (phi + e * sin(f))
            + 0.75 * s2 * sin(2*f + 2*g)
            + 0.75 * e * s2 * sin(f + 2*g)
            + 0.25 * e * s2 * sin(3*f + 2*g)
        )

    vars_base = [L, G, H, l, g, h]
    result = 0.0
    for j in range(3):
        Lj_idx = j
        lj_idx = j + 3
        v_p = list(vars_base); v_p[lj_idx] += eps
        v_m = list(vars_base); v_m[lj_idx] -= eps
        dp1_dlj = (p1(*v_p) - p1(*v_m)) / (2*eps)

        v_p = list(vars_base); v_p[Lj_idx] += eps
        v_m = list(vars_base); v_m[Lj_idx] -= eps
        dS1_dLj = (S1(*v_p) - S1(*v_m)) / (2*eps)
        dp1_dLj = (p1(*v_p) - p1(*v_m)) / (2*eps)

        v_p = list(vars_base); v_p[lj_idx] += eps
        v_m = list(vars_base); v_m[lj_idx] -= eps
        dS1_dlj = (S1(*v_p) - S1(*v_m)) / (2*eps)

        result += dp1_dlj * dS1_dLj - dp1_dLj * dS1_dlj
    return result


# ==========================================================================
# TEST: Compare BH61 Eq(14) δp₁ against our verified Poisson bracket {p₁,S₁}
# ==========================================================================
print("=" * 95)
print("CRITICAL TEST: BH61 Eq(14) δp₁  vs  Verified Poisson bracket {p₁, S₁}")
print("=" * 95)
print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'BH61 Eq(14)':>18} {'PB (FD)':>18} {'Ratio':>12} {'RelErr':>12}")
print("-" * 95)

mu = MU
a = 1.0
L = sqrt(mu * a)

match_count = 0
mismatch_count = 0

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L * eta_val
    for I_deg in [30, 60, 85]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for g_deg in [0, 45, 90]:
            g_val = np.radians(g_deg)
            for l_val in [0.5, 1.5, 3.0]:
                h_val = 0.0

                bh61 = bh61_eq14_dp1(L, G, H, l_val, g_val, h_val, mu)
                pb = poisson_bracket_numerical(L, G, H, l_val, g_val, h_val, eps=1e-7, mu=mu)

                denom = max(abs(bh61), abs(pb), 1e-15)
                rel_err = abs(bh61 - pb) / denom
                ratio = bh61 / pb if abs(pb) > 1e-15 else float('inf')

                status = "MATCH" if rel_err < 1e-4 else "MISMATCH"
                if status == "MATCH":
                    match_count += 1
                else:
                    mismatch_count += 1

                print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | "
                      f"{bh61:18.10e} {pb:18.10e} {ratio:12.6f} {rel_err:12.2e} {status}")

print("=" * 95)
print(f"Results: {match_count} MATCH, {mismatch_count} MISMATCH out of 81 test cases")
print("=" * 95)

# If they don't match, also check: is BH61's expression perhaps {p1, S1} + something?
# Or is it related by a constant factor?
if mismatch_count > 0:
    print("\nAnalyzing the ratio pattern...")
    print("If BH61 Eq(14) = c * {p1,S1} for constant c, then ratios should be constant.")
    print("If ratios vary, the expressions are structurally different.")
