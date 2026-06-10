"""
Independent verification: compute D(dS1/dl) using the CARTESIAN definition
D = -sum_i xi_i * d/d(xi_i)
where xi_i = dx_i/dt are Cartesian velocity components.

This bypasses the Delaunay D operator entirely and tests from first principles.
If the Cartesian and Delaunay D operators agree, the definition is correct.
If they disagree, there's an error in BH61's Eq(5) or Eq(10).
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

def delaunay_to_cartesian(L, G, H, l, g, h, mu=MU):
    """Convert Delaunay variables to Cartesian (x, y, z, vx, vy, vz)."""
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G  # cos I
    sinI = sqrt(max(1 - theta**2, 0))
    I = np.arccos(theta)
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)

    # Position and velocity in orbital frame
    r = a * (1 - e * cos(E))
    p = a * (1 - e**2)  # semi-latus rectum
    h_am = sqrt(mu * p)  # specific angular momentum

    # In perifocal frame (x along eccentricity vector, y perpendicular in orbit plane)
    x_pf = r * cos(f)
    y_pf = r * sin(f)
    vx_pf = -mu / h_am * sin(f)
    vy_pf = mu / h_am * (e + cos(f))

    # Rotation from perifocal to inertial: R = R3(-Omega) R1(-I) R3(-omega)
    Omega = h  # longitude of ascending node
    omega = g  # argument of periapsis

    cO, sO = cos(Omega), sin(Omega)
    cI, sI = cos(I), sin(I)
    cw, sw = cos(omega), sin(omega)

    # Rotation matrix columns
    Px = cO*cw - sO*sw*cI
    Py = -cO*sw - sO*cw*cI
    Qx = sO*cw + cO*sw*cI
    Qy = -sO*sw + cO*cw*cI
    Wx = sw*sI
    Wy = cw*sI

    x = Px*x_pf + Py*y_pf
    y = Qx*x_pf + Qy*y_pf
    z = Wx*x_pf + Wy*y_pf

    vx = Px*vx_pf + Py*vy_pf
    vy = Qx*vx_pf + Qy*vy_pf
    vz = Wx*vx_pf + Wy*vy_pf

    return np.array([x, y, z, vx, vy, vz])

def cartesian_to_delaunay(state, mu=MU):
    """Convert Cartesian state to Delaunay variables."""
    x, y, z, vx, vy, vz = state
    r_vec = np.array([x, y, z])
    v_vec = np.array([vx, vy, vz])
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    # Specific energy and semi-major axis
    energy = v**2 / 2 - mu / r
    a = -mu / (2 * energy)

    # Angular momentum
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)

    # Inclination
    I = np.arccos(h_vec[2] / h)

    # Node vector
    n_vec = np.cross(np.array([0, 0, 1]), h_vec)
    n = np.linalg.norm(n_vec)

    # Eccentricity vector
    e_vec = np.cross(v_vec, h_vec) / mu - r_vec / r
    e = np.linalg.norm(e_vec)

    # Longitude of ascending node
    if n > 1e-14:
        Omega = np.arccos(n_vec[0] / n)
        if n_vec[1] < 0:
            Omega = 2*pi - Omega
    else:
        Omega = 0.0

    # Argument of periapsis
    if n > 1e-14 and e > 1e-14:
        omega = np.arccos(np.clip(np.dot(n_vec, e_vec) / (n * e), -1, 1))
        if e_vec[2] < 0:
            omega = 2*pi - omega
    else:
        omega = 0.0

    # True anomaly
    if e > 1e-14:
        f = np.arccos(np.clip(np.dot(e_vec, r_vec) / (e * r), -1, 1))
        if np.dot(r_vec, v_vec) < 0:
            f = 2*pi - f
    else:
        f = 0.0

    # Eccentric anomaly
    E = 2 * np.arctan2(sqrt(1 - e) * sin(f / 2), sqrt(1 + e) * cos(f / 2))

    # Mean anomaly
    M = E - e * sin(E)
    if M < 0:
        M += 2*pi

    # Delaunay variables
    L = sqrt(mu * a)
    G = L * sqrt(1 - e**2)
    H = G * cos(I)
    l = M
    g = omega
    h = Omega

    return L, G, H, l, g, h

def dS1_dl_from_delaunay(L, G, H, l, g, h, mu=MU):
    a = L**2 / mu
    eta = G / L
    e = sqrt(max(1 - eta**2, 0))
    theta = H / G
    E = solve_kepler(l, e)
    f = true_from_eccentric(E, e)
    r = a * (1 - e * cos(E))
    Gamma = mu**2 / (a**3 * eta**3)
    B0 = -0.5 + 1.5 * theta**2
    B1p = 1.5 * (1 - theta**2)
    return Gamma * B0 - mu**2 / r**3 * (B0 + B1p * cos(2*f + 2*g))

def dS1_dl_from_cartesian(state, mu=MU):
    """Compute dS1/dl by converting Cartesian state to Delaunay."""
    L, G, H, l, g, h = cartesian_to_delaunay(state, mu)
    return dS1_dl_from_delaunay(L, G, H, l, g, h, mu)


def D_cartesian(func_cart, state, mu=MU, eps=1e-7):
    """Compute D(F) using the CARTESIAN definition: D = -sum_i v_i * dF/dv_i.

    Only velocity components (indices 3,4,5) are perturbed.
    """
    result = 0.0
    for i in range(3, 6):  # vx, vy, vz
        sp = state.copy(); sp[i] += eps
        sm = state.copy(); sm[i] -= eps
        dF_dvi = (func_cart(sp) - func_cart(sm)) / (2 * eps)
        result += -state[i] * dF_dvi  # D = -v_i * d/dv_i
    return result


def D_delaunay(L, G, H, l, g, h, mu=MU, eps=1e-7):
    """Compute D(dS1/dl) using Delaunay definition."""
    a = L**2/mu; e = sqrt(1-(G/L)**2); eta = G/L
    E = solve_kepler(l, e); f = true_from_eccentric(E, e); r = a*(1-e*cos(E))
    p1 = L*(2*a/r - 1); p2 = G; p3 = H
    q1 = 2*e*sin(E) + 2*eta*sin(f)/e if e > 1e-14 else 0.0
    q2 = -2*sin(f)/e if e > 1e-14 else 0.0
    D_acts = [-p1, -p2, -p3, q1, q2, 0.0]
    base = [L, G, H, l, g, h]
    result = 0.0
    for i in range(6):
        if abs(D_acts[i]) < 1e-20: continue
        bp = list(base); bp[i] += eps
        bm = list(base); bm[i] -= eps
        result += D_acts[i] * (dS1_dl_from_delaunay(*bp) - dS1_dl_from_delaunay(*bm)) / (2*eps)
    return result


# =====================================================================
# Test: Compare Cartesian D vs Delaunay D
# =====================================================================
print("=" * 90)
print("CRITICAL TEST: Cartesian D  vs  Delaunay D  on dS1/dl")
print("=" * 90)
print()

a_val = 1.0
L_val = sqrt(MU * a_val)

print(f"{'e':>6} {'I':>5} {'g':>5} {'l':>5} | {'D_Cartesian':>18} {'D_Delaunay':>18} {'6*G*B0*rho':>18} {'Cart/Del':>10}")
print("-" * 95)

for e_val in [0.01, 0.1, 0.3]:
    eta_val = sqrt(1 - e_val**2)
    G_val = L_val * eta_val
    for I_deg in [30, 60, 85]:
        theta_val = cos(np.radians(I_deg))
        H_val = G_val * theta_val
        for g_deg in [0, 45, 90]:
            g_val = np.radians(g_deg)
            for l_val in [0.5, 1.5, 3.0]:
                # Get Cartesian state
                state = delaunay_to_cartesian(L_val, G_val, H_val, l_val, g_val, 0.0)

                # Cartesian D
                D_cart = D_cartesian(dS1_dl_from_cartesian, state, eps=1e-7)

                # Delaunay D
                D_del = D_delaunay(L_val, G_val, H_val, l_val, g_val, 0.0, eps=1e-7)

                # Closed form
                a, e, eta, theta, f, r, E = (a_val, e_val, eta_val, theta_val,
                    true_from_eccentric(solve_kepler(l_val, e_val), e_val),
                    a_val*(1-e_val*cos(solve_kepler(l_val, e_val))),
                    solve_kepler(l_val, e_val))
                rho = a / r
                Gamma = MU**2 / (a**3 * eta**3)
                B0 = -0.5 + 1.5*theta**2
                formula = 6*Gamma*B0*rho

                ratio = D_cart / D_del if abs(D_del) > 1e-15 else float('inf')

                print(f"{e_val:6.2f} {I_deg:5d} {g_deg:5d} {l_val:5.1f} | "
                      f"{D_cart:18.10e} {D_del:18.10e} {formula:18.10e} {ratio:10.6f}")

print()
print("If Cart/Del ratio = 1.0 everywhere, the D operator definitions agree.")
print("If they disagree, BH61 Eq(5)/(10) has an error in p_j or q_j.")
