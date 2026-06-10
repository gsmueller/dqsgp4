import numpy as np

def solve_kepler(l, e, tol=1e-15):
    """Solves Kepler's equation l = E - e*sin(E) for E using Newton's method."""
    E = l if e < 0.8 else np.pi
    for _ in range(100):
        f_val = E - e * np.sin(E) - l
        f_prime = 1 - e * np.cos(E)
        delta = f_val / f_prime
        E -= delta
        if abs(delta) < tol:
            break
    return E

def get_f_from_E(E, e):
    """Computes the true anomaly f from the eccentric anomaly E."""
    half_f = np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                        np.sqrt(1 - e) * np.cos(E / 2))
    return 2 * half_f

def evaluate_H(f, g, e, theta):
    """Evaluates the integrand factor H."""
    A = (3 * theta**2 - 1) / 2
    B = 3 * (1 - theta**2) / 2
    cos_2fg = np.cos(2 * (f + g))
    term_bracket = A + B * cos_2fg
    return (1 + e * np.cos(f))**3 * cos_2fg * term_bracket

def compute_J(e, theta, g, N=4096):
    """Computes the l-average J using the periodic trapezoidal rule."""
    l_nodes = np.linspace(0, 2 * np.pi, N, endpoint=False)
    h_sum = 0.0
    for l in l_nodes:
        E = solve_kepler(l, e)
        f = get_f_from_E(E, e)
        h_sum += evaluate_H(f, g, e, theta)
    return h_sum / N

samples = [
    (0.10,  0.5, 0.30),
    (0.30,  0.7, 1.10),
    (0.50, -0.4, 2.40),
    (0.70,  0.6, 4.70),
    (0.85, -0.8, 0.00),
]

gemini_vals = [0.36444840, 0.23724300, 0.44437500, 0.29760000, 0.16593414]

print(f"{'Sample':<7}{'e':<7}{'theta':<8}{'g':<7}{'J (numerical)':<22}{'(1/2) B eta^3':<22}{'Gemini':<14}{'Gemini/analytical':<10}")
print("-" * 100)

for i, (e, theta, g) in enumerate(samples, 1):
    val = compute_J(e, theta, g)
    eta = np.sqrt(1 - e**2)
    B = 3 * (1 - theta**2) / 2
    analytical = 0.5 * B * eta**3
    ratio = gemini_vals[i-1] / analytical
    print(f"{i:<7}{e:<7.2f}{theta:<8.2f}{g:<7.2f}{val:<22.14f}{analytical:<22.14f}{gemini_vals[i-1]:<14.8f}{ratio:<10.4f}")
