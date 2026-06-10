"""
Step 1: Ratio Test — Are our dS1/dl and Brouwer's dS1/dl exact negatives?

Our form:     dS1/dl = Gamma*B0 - mu^2*(B0 + B1p*cos(2f+2g))/r^3
Brouwer form: dS1/dl = (mu^2/L^3)*(B0*(rho^3 - eta^{-3}) + B1p*rho^3*cos(2f+2g))

If ratio = -1 everywhere, D (linear operator) must give D(Brouwer) = -D(ours).
Tests with both mu=1,a=1 (normalized) and mu=398600,a=7000 (physical units).
"""
import numpy as np
from numpy import pi, sqrt, sin, cos, arctan2

def solve_kepler(M, e, tol=1e-15):
    E = M + e * sin(M)
    for _ in range(50):
        dE = (M - E + e * sin(E)) / (1 - e * cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def true_from_eccentric(E, e):
    return 2 * arctan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))

def dS1_dl_ours(L, G, H, l, g, h, mu):
    """Our form: Gamma*B0 - mu^2*(B0 + B1p*cos(2f+2g))/r^3."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))

    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    Gamma = mu**2 / (a**3 * eta**3)

    return Gamma * B0 - mu**2 * (B0 + B1p * cos(2*f + 2*g)) / r**3

def dS1_dl_brouwer(L, G, H, l, g, h, mu):
    """Brouwer form: (mu^2/L^3)*(B0*(rho^3 - eta^{-3}) + B1p*rho^3*cos(2f+2g))."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))

    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    rho = a / r
    n = mu**2 / L**3

    sigma1 = rho**3 - 1.0 / eta**3
    sigma2 = rho**3 * cos(2*f + 2*g)

    return n * (B0 * sigma1 + B1p * sigma2)


def run_test(mu, a_val, label):
    print(f"\n{'='*80}")
    print(f"RATIO TEST: {label}  (mu={mu}, a={a_val})")
    print(f"{'='*80}")
    print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'Ours':>18} {'Brouwer':>18} {'Ratio':>14}")
    print("-" * 80)

    ratios = []
    for e_val in [0.01, 0.1, 0.3]:
        eta_val = sqrt(1 - e_val**2)
        L = sqrt(mu * a_val)
        G = L * eta_val
        for I_deg in [30, 60, 85]:
            theta_val = cos(np.radians(I_deg))
            H = G * theta_val
            for g_deg in [0, 45, 90]:
                g_val = np.radians(g_deg)
                h_val = 0.0
                for l_val in [0.5, 1.5, 3.0]:
                    ours = dS1_dl_ours(L, G, H, l_val, g_val, h_val, mu)
                    brow = dS1_dl_brouwer(L, G, H, l_val, g_val, h_val, mu)

                    if abs(brow) > 1e-30:
                        ratio = ours / brow
                    else:
                        ratio = float('nan')
                    ratios.append(ratio)

                    print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | "
                          f"{ours:18.10e} {brow:18.10e} {ratio:14.10f}")

    ratios = np.array([r for r in ratios if np.isfinite(r)])
    print(f"\nRatio statistics:")
    print(f"  Mean:   {np.mean(ratios):.15f}")
    print(f"  Std:    {np.std(ratios):.2e}")
    print(f"  Min:    {np.min(ratios):.15f}")
    print(f"  Max:    {np.max(ratios):.15f}")

    if np.all(np.abs(ratios - (-1.0)) < 1e-12):
        print(f"  RESULT: ratio = -1 CONFIRMED to machine precision")
    else:
        print(f"  RESULT: ratio is NOT -1 everywhere!")
        print(f"  Max deviation from -1: {np.max(np.abs(ratios - (-1.0))):.2e}")

    return ratios


# ============================================================
# Algebraic verification
# ============================================================
print("=" * 80)
print("ALGEBRAIC CHECK")
print("=" * 80)
print()
print("Our form:     Gamma*B0 - mu^2*(B0 + B1p*cos(2f+2g))/r^3")
print("Brouwer form: n*(B0*(rho^3 - 1/eta^3) + B1p*rho^3*cos(2f+2g))")
print()
print("Using mu^2/r^3 = mu^2*(a/r)^3/a^3 = Gamma*eta^3*rho^3:")
print("  Our = Gamma*B0 - Gamma*eta^3*rho^3*(B0 + B1p*cos)")
print("      = Gamma*(B0 - eta^3*rho^3*B0 - eta^3*rho^3*B1p*cos)")
print("      = Gamma*(B0*(1 - eta^3*rho^3) - B1p*eta^3*rho^3*cos)")
print()
print("Brouwer = n*(B0*rho^3 - B0/eta^3 + B1p*rho^3*cos)")
print()
print("If ratio=-1: Gamma*(B0*(1 - eta^3*rho^3) - B1p*eta^3*rho^3*cos)")
print("           = -n*(B0*rho^3 - B0/eta^3 + B1p*rho^3*cos)")
print()
print("Matching constant (no rho, no cos): Gamma*B0 = n*B0/eta^3")
print("  => Gamma = n/eta^3")
print("  => mu^2/(a^3*eta^3) = mu^2/(L^3*eta^3)")
print("  => a^3 = L^3 = (mu*a)^{3/2}")
print("  => a^3 = mu^{3/2} * a^{3/2}")
print("  => a^{3/2} = mu^{3/2}")
print("  => a = mu")
print()
print("This holds ONLY when a = mu (i.e., a=1 with mu=1).")
print("For general a,mu the ratio should NOT be -1!")

# ============================================================
# Numerical tests
# ============================================================
ratios_norm = run_test(mu=1.0, a_val=1.0, label="Normalized (mu=1, a=1)")
ratios_phys = run_test(mu=398600.4418, a_val=7000.0, label="Physical units (mu=398600, a=7000)")

# Also test with mu=1, a=2 to break the a=mu coincidence
ratios_a2 = run_test(mu=1.0, a_val=2.0, label="Normalized (mu=1, a=2)")

print("\n" + "=" * 80)
print("CONCLUSIONS")
print("=" * 80)
print()
print("If ratio = -1 only for a=mu, then the 'paradox' was an artifact of testing")
print("with normalized units a=1, mu=1. The two forms are NOT the same function")
print("in general, and the D operator legitimately gives different results.")
print()
print("The relationship is: our dS1/dl = -(Gamma/n) * Brouwer's dS1/dl")
print("where Gamma/n = mu^2/(a^3*eta^3) / (mu^2/L^3) = L^3/(a^3*eta^3)")
print("     = (mu*a)^{3/2} / (a^3*eta^3) = mu^{3/2} / (a^{3/2}*eta^3)")
print()
print("This is NOT constant -- it depends on a (velocity-dependent).")
print("So D(our form) != -(Gamma/n) * D(Brouwer form) in general.")
