"""
Structural analysis: Can {p₁, S₁} be reduced to BH61 Eq (14) form?

Instead of comparing against the (OCR-error-prone) transcription, we ask:
Does our verified Poisson bracket have the STRUCTURE that BH61 claims?

BH61 Eq (14) claims δp₁ depends on (a,e,η,θ,f,g) through:
  - A constant (in l) term proportional to B₀ η⁻³
  - Terms proportional to (a/r)³ and (a/r)⁴
  - cos(2g+2f) dependence only (no sin terms, no other harmonics of g)

If our {p₁, S₁} has this structure, we can derive BH61 from foundations.

Approach: Numerically decompose {p��, S₁} into its harmonic content in f
at fixed (a, e, θ, g) to see what Fourier structure emerges.
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

def poisson_bracket_numerical(L, G, H, l, g, h, eps=1e-7, mu=MU):
    """Compute {p1, S1} by numerical finite differences on Delaunay variables."""
    def p1_func(L, G, H, l, g, h):
        a = L**2 / mu
        e = sqrt(1 - (G/L)**2)
        E = solve_kepler(l, e)
        r = a * (1 - e * cos(E))
        return L * (2*a/r - 1)

    def S1_func(L, G, H, l, g, h):
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
        dp1_dlj = (p1_func(*v_p) - p1_func(*v_m)) / (2*eps)

        v_p = list(vars_base); v_p[Lj_idx] += eps
        v_m = list(vars_base); v_m[Lj_idx] -= eps
        dS1_dLj = (S1_func(*v_p) - S1_func(*v_m)) / (2*eps)
        dp1_dLj = (p1_func(*v_p) - p1_func(*v_m)) / (2*eps)

        v_p = list(vars_base); v_p[lj_idx] += eps
        v_m = list(vars_base); v_m[lj_idx] -= eps
        dS1_dlj = (S1_func(*v_p) - S1_func(*v_m)) / (2*eps)

        result += dp1_dlj * dS1_dLj - dp1_dLj * dS1_dlj
    return result


# ==========================================================================
# ANALYSIS 1: Fourier decomposition of {p₁, S₁} in mean anomaly l
# ==========================================================================
print("=" * 80)
print("ANALYSIS 1: Fourier decomposition of {p1, S1} in l (mean anomaly)")
print("=" * 80)
print()
print("If BH61 Eq(14) is correct, {p₁,S₁} should decompose as:")
print("  a₀ + Σ [aₙ cos(nl) + bₙ sin(nl)]")
print("where the harmonics relate to (a/r)^k cos(jf+2g) through Kepler.")
print()

mu = MU
a_val = 1.0
L = sqrt(mu * a_val)

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L * eta_val
    for I_deg in [60]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        for g_deg in [45]:
            g_val = np.radians(g_deg)

            # Sample {p₁, S₁} over a full period in l
            N = 256
            l_vals = np.linspace(0.001, 2*pi - 0.001, N)  # avoid l=0 singularity for low e
            pb_vals = np.array([
                poisson_bracket_numerical(L, G, H, lv, g_val, 0.0, eps=1e-7)
                for lv in l_vals
            ])

            # Fourier decomposition
            # a_n = (2/N) Σ pb(l_k) cos(n l_k), b_n = (2/N) Σ pb(l_k) sin(n l_k)
            a0 = np.mean(pb_vals)
            print(f"e={e_val}, I={I_deg}°, g={g_deg}°:")
            print(f"  Mean (a₀) = {a0:.10e}")

            for n in range(1, 7):
                an = 2 * np.mean(pb_vals * cos(n * l_vals))
                bn = 2 * np.mean(pb_vals * sin(n * l_vals))
                amp = sqrt(an**2 + bn**2)
                if amp > 1e-12:
                    print(f"  n={n}: a_{n}={an:+.8e}, b_{n}={bn:+.8e}, amp={amp:.8e}, amp/a₀={amp/abs(a0):.4e}")
                else:
                    print(f"  n={n}: amplitude < 1e-12 (negligible)")
            print()


# ==========================================================================
# ANALYSIS 2: Dependence on g — extract cos(2g) structure
# ==========================================================================
print("=" * 80)
print("ANALYSIS 2: Fourier decomposition of {p₁, S₁} in g")
print("=" * 80)
print()
print("BH61 Eq(14) predicts only cos(2g+2f) dependence on g,")
print("which after averaging over l gives a cos(2g) long-period term.")
print()

for e_val in [0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L * eta_val
    for I_deg in [60]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val

        # At fixed l, vary g and decompose
        l_fixed = 1.5
        N_g = 128
        g_vals = np.linspace(0, 2*pi, N_g, endpoint=False)
        pb_vs_g = np.array([
            poisson_bracket_numerical(L, G, H, l_fixed, gv, 0.0, eps=1e-7)
            for gv in g_vals
        ])

        g0 = np.mean(pb_vs_g)
        print(f"e={e_val}, I={I_deg}°, l={l_fixed}:")
        print(f"  Mean over g = {g0:.10e}")
        for m in range(1, 5):
            am = 2 * np.mean(pb_vs_g * cos(m * g_vals))
            bm = 2 * np.mean(pb_vs_g * sin(m * g_vals))
            amp = sqrt(am**2 + bm**2)
            print(f"  m={m}: cos coeff={am:+.8e}, sin coeff={bm:+.8e}, amp={amp:.8e}")
        print()


# ==========================================================================
# ANALYSIS 3: Does {p₁,S₁} factor as BH61 claims?
# BH61: δp₁ = (3μ²/L³)[B₀·F₀(a/r,η) + B���·F₁(a/r)·cos(2g+2f)]
# where F₀ = -η⁻³ + (a/r)³(1-2a/r), F₁ = (a/r)³(1-2a/r)
# Note F₁ = F₀ + η⁻³ ... so F₀ and F₁ share the same (a/r)³(1-2a/r) piece
# ==========================================================================
print("=" * 80)
print("ANALYSIS 3: Testing BH61's claimed factored form")
print("=" * 80)
print()
print("BH61 claims: δp₁/(3μ²/L³) = B₀[-η⁻³ + (a/r)³(1-2a/r)]")
print("                             + B₁(a/r)³(1-2a/r)cos(2g+2f)")
print()
print("This means at g=0: δp₁/(3μ²/L³) = [B₀+B₁cos(2f)]·(a/r)³(1-2a/r) - B₀η⁻³")
print()
print("Test: compute R(f) = [{p₁,S₁}/(3μ²/L³) + B₀η⁻³] / [(a/r)³(1-2a/r)]")
print("If BH61 is right, R(f) = B₀ + B₁cos(2f+2g)")
print()

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G = L * eta_val
    for I_deg in [60]:
        theta_val = cos(np.radians(I_deg))
        H = G * theta_val
        B0 = -0.5 + 1.5*theta_val**2
        B1 = 1.5 - 1.5*theta_val**2

        for g_deg in [0, 45, 90]:
            g_val = np.radians(g_deg)
            print(f"  e={e_val}, I={I_deg}°, g={g_deg}°:  B₀={B0:.4f}, B₁={B1:.4f}")

            for l_val in [0.5, 1.5, 3.0, 5.0]:
                E = solve_kepler(l_val, e_val)
                f = true_from_eccentric(E, e_val)
                r = a_val * (1 - e_val * cos(E))
                ar = a_val / r

                pb = poisson_bracket_numerical(L, G, H, l_val, g_val, 0.0, eps=1e-7)
                prefactor = 3.0 * mu**2 / L**3

                # What BH61 predicts
                bh61_pred = prefactor * (B0 * (-1/eta_val**3 + ar**3*(1-2*ar))
                                        + B1 * ar**3*(1-2*ar)*cos(2*g_val+2*f))

                # Residual
                resid = pb - bh61_pred
                if abs(pb) > 1e-12:
                    rel = resid / pb
                else:
                    rel = resid

                # Also try the factored ratio
                denom = ar**3 * (1 - 2*ar)
                if abs(denom) > 1e-10:
                    R = (pb / prefactor + B0 / eta_val**3) / denom
                    expected_R = B0 + B1 * cos(2*g_val + 2*f)
                    print(f"    l={l_val:.1f}: R={R:+.6f}, B₀+B₁cos(2f+2g)={expected_R:+.6f}, "
                          f"diff={R-expected_R:+.4e}, PB={pb:+.8e}, BH61={bh61_pred:+.8e}, rel_resid={rel:+.4e}")
                else:
                    print(f"    l={l_val:.1f}: denom≈0 (a/r≈0.5), skipping ratio test. "
                          f"PB={pb:+.8e}, BH61={bh61_pred:+.8e}, rel_resid={rel:+.4e}")
            print()

print("\nIf R ≈ B₀ + B₁cos(2f+2g) consistently, then BH61 Eq(14) structure is confirmed")
print("and we CAN derive BH61 from our Poisson bracket foundations.")
