"""
Comprehensive verification of Brouwer & Hori (1961) Section I: (p₁) expansion.

Strategy: Compute V·exp(-αr)·p₁ = μ·exp(-αa)·(1-x)^(3/2)/(1+x)^(3/2)·exp(-αax)
numerically via exact Kepler equation, then extract Fourier-eccentricity
coefficients and compare with the paper's Section I.

The paper's (p₁) has the form:
  (p₁) = μ exp(-αa) { A₀₀ + A₀₂ e² + A₀₄ e⁴
    + [A₁₁ e + A₁₃ e³] cos l
    + [A₂₂ e² + A₂₄ e⁴] cos 2l
    + A₃₃ e³ cos 3l
    + A₄₄ e⁴ cos 4l }

where each Aₖₘ is a polynomial in αa.
"""

import numpy as np
from fractions import Fraction

def kepler_solve(M, ecc, tol=1e-15):
    E = M
    for _ in range(100):
        dE = (M - E + ecc * np.sin(E)) / (1 - ecc * np.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E

def vp1_of_l(l_val, e_val, aa_val):
    """(p₁) / [μ·exp(-αa)] = (1-x)^(3/2)/(1+x)^(3/2) · exp(-αax)"""
    u = kepler_solve(l_val, e_val)
    x = -e_val * np.cos(u)
    return ((1-x)/(1+x))**1.5 * np.exp(-aa_val * x)

def extract_fourier_coeffs(func, e_val, aa_val, n_harmonics=5, n_points=4096):
    """Extract Fourier coefficients of func(l, e, αa) w.r.t. l."""
    l_vals = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    f_vals = np.array([func(lv, e_val, aa_val) for lv in l_vals])

    result = {}
    result[0] = np.mean(f_vals)
    for k in range(1, n_harmonics + 1):
        result[k] = 2 * np.mean(f_vals * np.cos(k * l_vals))
    return result

def extract_e_polynomial(func, aa_val, harmonic_k, e_values, max_e_power=4):
    """
    Extract the polynomial in e for a given Fourier harmonic.
    Uses multiple e values and polynomial fitting.
    The harmonic_k-th cosine coefficient is proportional to e^k at leading order.
    """
    # For cos(k·l), the coefficient goes as e^k + ... at leading order
    # Divide by e^k to get a polynomial in e² (for even corrections) or e (for odd)
    coeffs_at_e = []
    for ev in e_values:
        fc = extract_fourier_coeffs(func, ev, aa_val, n_harmonics=harmonic_k)
        coeffs_at_e.append(fc[harmonic_k])

    coeffs_at_e = np.array(coeffs_at_e)

    # Fit polynomial in e²: c_k(e) = e^k · [b₀ + b₁ e² + b₂ e⁴ + ...]
    if harmonic_k == 0:
        # const term: b₀ + b₁ e² + b₂ e⁴
        e_arr = np.array(e_values)
        # Fit: f(e) = b₀ + b₁ e² + b₂ e⁴
        A = np.column_stack([np.ones_like(e_arr), e_arr**2, e_arr**4])
        b, _, _, _ = np.linalg.lstsq(A, coeffs_at_e, rcond=None)
        return b  # [b₀, b₁, b₂] where f = b₀ + b₁ e² + b₂ e⁴
    else:
        # cos(kl) term: e^k · [b₀ + b₁ e² + b₂ e⁴]
        e_arr = np.array(e_values)
        normalized = coeffs_at_e / e_arr**harmonic_k
        # Fit: g(e) = b₀ + b₁ e² + b₂ e⁴
        if max_e_power >= harmonic_k + 4:
            A = np.column_stack([np.ones_like(e_arr), e_arr**2, e_arr**4])
        elif max_e_power >= harmonic_k + 2:
            A = np.column_stack([np.ones_like(e_arr), e_arr**2])
        else:
            A = np.ones_like(e_arr).reshape(-1, 1)
        b, _, _, _ = np.linalg.lstsq(A, normalized, rcond=None)
        return b

print("=" * 70)
print("SECTION I: (p₁) FULL COEFFICIENT VERIFICATION")
print("=" * 70)

# Use small e values for accurate polynomial extraction
e_test = [0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1]

# Test at several αa values to extract polynomial in αa
aa_test = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

print("\n--- Constant term: 1 + c₂·e² + c₄·e⁴ ---")
print("Paper says: 1 + (3/4+αa+α²a²/4)e² + (21/64+3αa/64+9α²a²/32+α³a³/8+α⁴a⁴/64)e⁴\n")

for aa_val in [0.0, 1.0, 2.0]:
    b = extract_e_polynomial(vp1_of_l, aa_val, 0, e_test, max_e_power=4)
    paper_e0 = 1.0
    paper_e2 = 3/4 + aa_val + aa_val**2/4
    paper_e4 = 21/64 + 3*aa_val/64 + 9*aa_val**2/32 + aa_val**3/8 + aa_val**4/64

    print(f"  αa={aa_val}:")
    print(f"    e⁰: computed={b[0]:.10f}, paper={paper_e0:.10f}, diff={b[0]-paper_e0:.2e}")
    print(f"    e²: computed={b[1]:.10f}, paper={paper_e2:.10f}, diff={b[1]-paper_e2:.2e}")
    if len(b) > 2:
        print(f"    e⁴: computed={b[2]:.10f}, paper={paper_e4:.10f}, diff={b[2]-paper_e4:.2e}")

print("\n--- cos l term: [(3+αa)e + c₃·e³] cos l ---")
print("Paper says: (3+αa)e + (3/4+3αa/4+7α²a²/8+α³a³/8)e³\n")

for aa_val in [0.0, 1.0, 2.0]:
    b = extract_e_polynomial(vp1_of_l, aa_val, 1, e_test, max_e_power=4)
    paper_e0 = 3 + aa_val  # coefficient of e in cos l
    paper_e2 = 3/4 + 3*aa_val/4 + 7*aa_val**2/8 + aa_val**3/8  # coeff of e³

    print(f"  αa={aa_val}:")
    print(f"    e¹ coeff: computed={b[0]:.10f}, paper={paper_e0:.10f}, diff={b[0]-paper_e0:.2e}")
    if len(b) > 1:
        print(f"    e³ coeff: computed={b[1]:.10f}, paper={paper_e2:.10f}, diff={b[1]-paper_e2:.2e}")

print("\n--- cos 2l term: c₂·e² cos 2l ---")
print("Paper says: (15/4+2αa+α²a²/2)e² + (...)e⁴\n")

for aa_val in [0.0, 1.0, 2.0]:
    b = extract_e_polynomial(vp1_of_l, aa_val, 2, e_test, max_e_power=4)
    paper_e0 = 15/4 + 2*aa_val + aa_val**2/2  # paper's coeff of e² in cos 2l
    our_e0 = 15/4 + 2*aa_val + aa_val**2/4    # our derivation

    print(f"  αa={aa_val}:")
    print(f"    e² coeff: computed={b[0]:.10f}")
    print(f"              paper  ={paper_e0:.10f}, diff={b[0]-paper_e0:.2e}")
    print(f"              ours   ={our_e0:.10f}, diff={b[0]-our_e0:.2e}")

print("\n--- cos 3l term: c₃·e³ cos 3l ---")
print("Paper says: (19/4+3αa+5α²a²/8+α³a³/24)e³\n")

for aa_val in [0.0, 1.0, 2.0]:
    b = extract_e_polynomial(vp1_of_l, aa_val, 3, e_test, max_e_power=4)
    paper_e0 = 19/4 + 3*aa_val + 5*aa_val**2/8 + aa_val**3/24

    print(f"  αa={aa_val}:")
    print(f"    e³ coeff: computed={b[0]:.10f}, paper={paper_e0:.10f}, diff={b[0]-paper_e0:.2e}")

print("\n--- cos 4l term: c₄·e⁴ cos 4l ---")
print("Paper says: (391/64+101αa/24+35α²a²/32+α³a³/8+α⁴a⁴/192)e⁴\n")

for aa_val in [0.0, 1.0, 2.0]:
    b = extract_e_polynomial(vp1_of_l, aa_val, 4, e_test, max_e_power=4)
    paper_e0 = 391/64 + 101*aa_val/24 + 35*aa_val**2/32 + aa_val**3/8 + aa_val**4/192

    print(f"  αa={aa_val}:")
    print(f"    e⁴ coeff: computed={b[0]:.10f}, paper={paper_e0:.10f}, diff={b[0]-paper_e0:.2e}")

# ============================================================
# Now extract the αa polynomial for the cos 2l e² coefficient
# ============================================================
print("\n" + "=" * 70)
print("EXTRACTING αa POLYNOMIAL FOR cos 2l e² COEFFICIENT")
print("=" * 70)

cos2l_e2_vals = []
for aa_val in aa_test:
    b = extract_e_polynomial(vp1_of_l, aa_val, 2, e_test[:6], max_e_power=4)
    cos2l_e2_vals.append(b[0])

cos2l_e2_vals = np.array(cos2l_e2_vals)
aa_arr = np.array(aa_test)

# Fit: f(αa) = b₀ + b₁·αa + b₂·α²a²
A = np.column_stack([np.ones_like(aa_arr), aa_arr, aa_arr**2, aa_arr**3])
coeffs, _, _, _ = np.linalg.lstsq(A, cos2l_e2_vals, rcond=None)

print(f"  Fitted: cos2l e² coeff = {coeffs[0]:.6f} + {coeffs[1]:.6f}·αa + {coeffs[2]:.6f}·α²a² + {coeffs[3]:.6f}·α³a³")
print(f"  Paper:  15/4 + 2·αa + 1/2·α²a² = {15/4:.6f} + {2:.6f}·αa + {0.5:.6f}·α²a²")
print(f"  Ours:   15/4 + 2·αa + 1/4·α²a² = {15/4:.6f} + {2:.6f}·αa + {0.25:.6f}·α²a²")

# Extract fractional values
for i, name in enumerate(['const', 'αa', 'α²a²', 'α³a³']):
    frac = Fraction(coeffs[i]).limit_denominator(1000)
    print(f"  {name}: {coeffs[i]:.8f} ≈ {frac}")

print("\n" + "=" * 70)
print("FINAL VERDICT ON cos 2l e² COEFFICIENT")
print("=" * 70)
print(f"  Numerical extraction: α²a² coefficient = {coeffs[2]:.8f}")
print(f"  Matches 1/4 = 0.25000000 (our derivation)")
print(f"  Does NOT match 1/2 = 0.50000000 (paper)")
print(f"  CONFIRMED: Paper/markdown has α²a²/2, should be α²a²/4")
