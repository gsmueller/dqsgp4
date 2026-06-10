"""
Numerical verification for Chapter 5: Series Evaluation and Error Control.

Computes concrete values for all series discussed in Ch 5,
verifying convergence rates, term counts, cancellation losses,
and error bound tightness.
"""

from mpmath import mp, mpf, atan, sqrt, log10, fabs, pi, power, binomial, gamma, factorial

# Set high precision for reference values
mp.dps = 60

# ============================================================
# WGS84 parameters
# ============================================================
a_wgs84 = mpf('6378137.0')          # semi-major axis [m]
inv_f    = mpf('298.257223563')
f        = 1 / inv_f
e2       = 2*f - f*f
e        = sqrt(e2)
ep2      = e2 / (1 - e2)
ep       = sqrt(ep2)

print("=" * 72)
print("WGS84 Parameters")
print("=" * 72)
print(f"  e^2   = {float(e2):.15e}")
print(f"  e'^2  = {float(ep2):.15e}")
print(f"  e'    = {float(ep):.15e}")
print()

# ============================================================
# 1. q0 series (Eq. 5.11): 2*q0 = sum_{n=1}^{inf} (-1)^{n+1} * 4n/((2n+1)(2n+3)) * e'^{2n+1}
# ============================================================
print("=" * 72)
print("1. q0 SERIES EVALUATION (Eq. 5.11)")
print("=" * 72)

# Reference value via closed form at 60-digit precision (no cancellation issue at this precision)
q0_closed = mpf('0.5') * ((1 + 3/ep2) * atan(ep) - 3/ep)

# Series term-by-term
print(f"\n{'n':>3} {'term':>25} {'partial sum (2*q0)':>25} {'Leibniz bound':>15} {'|t_{n+1}/t_n|':>15}")
print("-" * 90)

partial = mpf(0)
terms = []
for n in range(1, 26):
    sign = (-1)**(n+1)
    coeff = mpf(4*n) / ((2*n+1) * (2*n+3))
    term = sign * coeff * power(ep, 2*n+1)
    terms.append(term)
    partial += term
    leibniz = fabs(term)

    ratio_str = ""
    if n >= 2:
        ratio = fabs(term / terms[-2])
        ratio_str = f"{float(ratio):.10f}"

    if n <= 10 or n == 25:
        print(f"{n:3d} {float(term):>25.18e} {float(partial):>25.18e} {float(leibniz):>15.6e} {ratio_str:>15}")

q0_series = partial / 2
print(f"\nq0 (series, 25 terms)  = {float(q0_series):.20e}")
print(f"q0 (closed form, 60dp) = {float(q0_closed):.20e}")
print(f"Difference             = {float(fabs(q0_series - q0_closed)):.6e}")

# Cancellation demonstration
print(f"\n--- Cancellation in closed form (Cor. 5.5.2) ---")
term_a = (1 + 3/ep2) * atan(ep)
term_b = 3/ep
print(f"(1 + 3/e'^2)*arctan(e') = {float(term_a):.18e}")
print(f"3/e'                    = {float(term_b):.18e}")
print(f"Difference              = {float(term_a - term_b):.18e}")
print(f"2*q0                    = {float(2*q0_closed):.18e}")

# At double precision, simulate cancellation
mp_save = mp.dps
mp.dps = 16  # simulate double precision
q0_closed_dbl = mpf('0.5') * ((1 + 3/ep2) * atan(ep) - 3/ep)
mp.dps = mp_save
print(f"\nq0 at 16-digit precision (closed) = {float(q0_closed_dbl):.18e}")
print(f"q0 at 60-digit precision (closed) = {float(q0_closed):.18e}")
digits_lost = -float(log10(fabs(q0_closed_dbl - q0_closed) / fabs(q0_closed)))
print(f"Reliable digits in 16dp closed form: {digits_lost:.1f}")

# Convergence rate verification (Eq. 5.19)
print(f"\n--- Convergence rate verification (Eq. 5.19) ---")
print(f"{'n':>3} {'actual |t_{n+1}/t_n|':>22} {'Eq.5.19 prediction':>22} {'asymptote ep2':>15}")
for n in range(1, 10):
    actual = fabs(terms[n] / terms[n-1])  # terms[n] is term at n+1
    predicted = mpf((n+1)*(2*n+1)) / (n*(2*n+5)) * ep2
    print(f"{n:3d} {float(actual):>22.15f} {float(predicted):>22.15f} {float(ep2):>15.12f}")

# Term count verification at various tolerances
print(f"\n--- Term counts at various tolerances ---")
for tau_exp in [-3, -12, -16, -34, -50]:
    tau = mpf(10) ** tau_exp
    N = 0
    partial_check = mpf(0)
    for n in range(1, 200):
        sign = (-1)**(n+1)
        coeff = mpf(4*n) / ((2*n+1) * (2*n+3))
        term = sign * coeff * power(ep, 2*n+1)
        partial_check += term
        if fabs(term) < tau:
            N = n
            break
    print(f"  tau = 10^{tau_exp:3d}: N = {N:3d} terms, |last term| = {float(fabs(term)):.3e}")

# ============================================================
# 2. q0' series (Eq. 5.23)
# ============================================================
print(f"\n{'=' * 72}")
print("2. q0' SERIES EVALUATION (Eq. 5.23)")
print("=" * 72)

partial_qp = mpf(0)
terms_qp = []
header_qp = "partial sum (q0')"
print(f"\n{'n':>3} {'term':>25} {header_qp:>25} {'Leibniz bound':>15}")
print("-" * 75)
for n in range(1, 12):
    sign = (-1)**(n+1)
    coeff = mpf(6) / ((2*n+1) * (2*n+3))
    term = sign * coeff * power(ep, 2*n)
    terms_qp.append(term)
    partial_qp += term
    print(f"{n:3d} {float(term):>25.18e} {float(partial_qp):>25.18e} {float(fabs(term)):>15.6e}")

print(f"\nq0' (series, 11 terms) = {float(partial_qp):.20e}")

# ============================================================
# 3. U0 series (arctan(e')/e' = 1 - e'^2/3 + e'^4/5 - ...)
# ============================================================
print(f"\n{'=' * 72}")
print("3. U0 AUXILIARY SERIES (arctan(e')/e')")
print("=" * 72)

partial_u0 = mpf(0)
print(f"\n{'n':>3} {'term':>25} {'partial sum':>25} {'Leibniz bound':>15}")
print("-" * 75)
for n in range(0, 12):
    sign = (-1)**n
    coeff = mpf(sign) / (2*n + 1)
    term = coeff * power(ep, 2*n)
    partial_u0 += term
    print(f"{n:3d} {float(term):>25.18e} {float(partial_u0):>25.18e} {float(fabs(term)):>15.6e}")

ref_u0 = atan(ep) / ep
print(f"\narctan(e')/e' (series) = {float(partial_u0):.20e}")
print(f"arctan(e')/e' (direct) = {float(ref_u0):.20e}")

# ============================================================
# 4. Geodetic binomial: (1 - e^2*sin^2(phi))^{-1/2} at phi = pi/4
# ============================================================
print(f"\n{'=' * 72}")
print("4. GEODETIC BINOMIAL alpha=-1/2 at phi=pi/4")
print("=" * 72)

phi = pi/4
sin2phi = mpf('0.5')  # sin^2(pi/4) = 1/2
x = -e2 * sin2phi
alpha = mpf('-0.5')

# Reference
ref_binom = power(1 + x, alpha)

partial_binom = mpf(0)
print(f"\n{'k':>3} {'C(alpha,k)*x^k':>25} {'partial sum':>25} {'geom tail bound':>15}")
print("-" * 75)
for k in range(0, 12):
    # Generalized binomial coefficient
    c = mpf(1)
    for j in range(k):
        c *= (alpha - j) / (j + 1)
    term = c * power(x, k)
    partial_binom += term

    # Geometric tail bound: |term| * r/(1-r) where r = |x|
    r = fabs(x)
    if r < 1:
        tail = fabs(term) * r / (1 - r)
    else:
        tail = mpf('inf')

    print(f"{k:3d} {float(term):>25.18e} {float(partial_binom):>25.18e} {float(tail):>15.6e}")

print(f"\nSeries (11 terms)  = {float(partial_binom):.20e}")
print(f"Direct evaluation  = {float(ref_binom):.20e}")
print(f"Difference         = {float(fabs(partial_binom - ref_binom)):.6e}")

# ============================================================
# 5. Error bound tightness
# ============================================================
print(f"\n{'=' * 72}")
print("5. ERROR BOUND TIGHTNESS (q0 series)")
print("=" * 72)

# Compute q0 series to 50 terms as "exact" reference
mp.dps = 60
ref_2q0 = mpf(0)
for n in range(1, 51):
    sign = (-1)**(n+1)
    coeff = mpf(4*n) / ((2*n+1) * (2*n+3))
    ref_2q0 += sign * coeff * power(ep, 2*n+1)

print(f"\n{'N':>3} {'Leibniz |t(N)|':>18} {'Geom tail':>18} {'Actual |R_N|':>18} {'Leibniz/actual':>15} {'Geom/actual':>15}")
print("-" * 100)
partial_check = mpf(0)
for n in range(1, 16):
    sign = (-1)**(n+1)
    coeff = mpf(4*n) / ((2*n+1) * (2*n+3))
    term = sign * coeff * power(ep, 2*n+1)
    partial_check += term

    leibniz = fabs(term)
    r_ratio = ep2
    geom = fabs(term) * r_ratio / (1 - r_ratio)
    actual = fabs(ref_2q0 - partial_check)

    if actual > 0:
        leib_ratio = float(leibniz / actual)
        geom_ratio = float(geom / actual)
    else:
        leib_ratio = float('inf')
        geom_ratio = float('inf')

    print(f"{n:3d} {float(leibniz):>18.6e} {float(geom):>18.6e} {float(actual):>18.6e} {leib_ratio:>15.1f} {geom_ratio:>15.4f}")

# ============================================================
# 6. TrackedValue walkthrough for q0
# ============================================================
print(f"\n{'=' * 72}")
print("6. TrackedValue WALKTHROUGH (q0 at double precision)")
print("=" * 72)

eps_mach = mpf(2)**(-52)
print(f"\neps_mach = {float(eps_mach):.6e}")
print(f"\nStep-by-step accumulation of arithmetic delta_p:")
print(f"{'n':>3} {'|term|':>15} {'|partial|':>15} {'arith delta_p':>15} {'trunc bound':>15} {'total delta_p':>15}")
print("-" * 85)

arith_dp = mpf(0)
partial_tv = mpf(0)
for n in range(1, 9):
    sign = (-1)**(n+1)
    coeff = mpf(4*n) / ((2*n+1) * (2*n+3))
    term = sign * coeff * power(ep, 2*n+1)

    # Each addition: rounding error ~ eps_mach * |partial + term|
    partial_tv += term
    arith_dp += eps_mach * fabs(partial_tv)

    trunc = fabs(term)  # Leibniz bound
    total = arith_dp + trunc

    print(f"{n:3d} {float(fabs(term)):>15.6e} {float(fabs(partial_tv)):>15.6e} {float(arith_dp):>15.6e} {float(trunc):>15.6e} {float(total):>15.6e}")

print(f"\nAt convergence (n=7, tau=1e-16):")
print(f"  Arithmetic delta_p = {float(arith_dp):.6e}")
print(f"  Truncation bound   = {float(trunc):.6e}")
print(f"  Total delta_p      = {float(total):.6e}")
print(f"  Dominated by: {'arithmetic rounding' if arith_dp > trunc else 'truncation'}")

print(f"\n{'=' * 72}")
print("VERIFICATION COMPLETE")
print("=" * 72)
