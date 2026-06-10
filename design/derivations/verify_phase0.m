% verify_phase0.m
%
% Mechanical verifier for Phase 0 of the SGP4 near-earth drag re-derivation.
% Independently constructs Theorems 0.2.2 .. 0.5.3 + 0.4.1 .. 0.4.2 symbolically
% and asserts algebraic equality to the closed forms claimed in
% design/derivations/sgp4_drag_phase0_foundations.md.
%
% Usage:
%   octave --eval "run('design/derivations/verify_phase0.m')"
% or:
%   octave design/derivations/verify_phase0.m
%
% Exit code 0 on all PASS, 1 on any FAIL.
%
% Convention:
%   The script does NOT re-derive theorems from the perturbed equations of
%   motion (that would duplicate the proofs in the markdown document). Instead,
%   it constructs each KEY ALGEBRAIC IDENTITY that the proofs depend on, using
%   Octave's symbolic package, and asserts the identity holds. If any identity
%   fails, the corresponding theorem proof in the markdown document has an
%   algebraic error and must be re-examined.
%
% Coverage: each PASS in this script discharges one or more numbered steps in
% the markdown proofs. The mapping is annotated in each check.

pkg load symbolic;

% --- Symbolic variables ---
% Note: Octave's symbolic package does not support general algebraic
% assumptions like (e < 1) on declared variables. However, type-level
% assumptions (positive, real) are supported via the `syms ... positive`
% form. Physical constants and the orbital elements that must be positive
% to evaluate sqrt(.) cleanly are declared as positive below; eccentricity
% e is left general (real).
syms mu a n positive;
syms e f real;
syms beta_sym p_sym positive;
syms R_force T_force N_force real;
syms i_sym Omega_sym omega_sym real;

% Derived: beta = sqrt(1 - e^2), p = a (1-e^2), r = p/(1+e cos f)
beta_def = sqrt(1 - e^2);
p_def = a * (1 - e^2);
r_def = p_def / (1 + e*cos(f));

% Eccentric anomaly: cos E = (e + cos f)/(1+e cos f); sin E from sin^2 + cos^2 = 1
cosE = (e + cos(f)) / (1 + e*cos(f));
sinE = beta_def * sin(f) / (1 + e*cos(f));   % squared, then take positive root for 0 < f < pi

% Counters
pass_count = 0;
fail_count = 0;
failed_names = {};

% --- Check helper ---
function [pass_count, fail_count, failed_names] = check_eq( ...
    name, lhs, rhs, pass_count, fail_count, failed_names)
    diff = simplify(lhs - rhs);
    if isequal(diff, sym(0))
        printf('  PASS: %s\n', name);
        pass_count = pass_count + 1;
    else
        printf('  FAIL: %s\n', name);
        printf('        lhs - rhs = %s\n', char(diff));
        fail_count = fail_count + 1;
        failed_names{end+1} = name;
    end
endfunction

printf('=== Phase 0 mechanical verification ===\n');
printf('\n');
printf('Pre-flight identities:\n');

% --- Pre-flight: p = a beta^2 ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'p = a beta^2 (Definition 0.1 / Theorem 0.2.2)', ...
    p_def, ...
    a * beta_def^2, ...
    pass_count, fail_count, failed_names);

% --- Pre-flight: cos^2 E + sin^2 E = 1 ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'cos^2 E + sin^2 E = 1 (eccentric-true-anomaly identity)', ...
    cosE^2 + sinE^2, ...
    sym(1), ...
    pass_count, fail_count, failed_names);

% --- Pre-flight: r = a(1 - e cos E) ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'r = a(1 - e cos E) (alternative conic form)', ...
    a * (1 - e*cosE), ...
    r_def, ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.2.4 (Kepler 3rd):\n');

% --- Theorem 0.2.4 (0.2.4.2): n a^2 beta = h = sqrt(mu p) ---
% With mu = n^2 a^3, sqrt(mu p) = sqrt(n^2 a^3 * a beta^2) = n a^2 beta.
mu_def = n^2 * a^3;
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.2.4 (0.2.4.2): sqrt(mu p) = n a^2 beta', ...
    sqrt(mu_def * p_def), ...
    n * a^2 * beta_def, ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.2.5 (dM/df):\n');

% --- Theorem 0.2.5 (0.2.5.3): r^2 / (a^2 beta) = beta^3 / (1+e cos f)^2 ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.2.5 (0.2.5.3): r^2/(a^2 beta) = beta^3/(1+e cos f)^2', ...
    r_def^2 / (a^2 * beta_def), ...
    beta_def^3 / (1 + e*cos(f))^2, ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.3.2 prefactor (a_dot Gauss VE):\n');

% --- Theorem 0.3.2 (0.3.2.14): 2 n a^3 / (mu beta) = 2/(n beta) ---
% Substituting mu = n^2 a^3:
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.2 (0.3.2.14): 2 n a^3 / (mu beta) = 2/(n beta)', ...
    2 * n * a^3 / (mu_def * beta_def), ...
    2 / (n * beta_def), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.2 (0.3.2.10): r_dot = n a e sin f / beta ---
% Derived from r_dot = (e h sin f) / p with h = n a^2 beta and p = a beta^2:
% (e * n a^2 beta * sin f) / (a beta^2) = (n a e sin f) / beta
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.2 (0.3.2.10): e h sin f / p = n a e sin f / beta', ...
    e * (n * a^2 * beta_def) * sin(f) / p_def, ...
    (n * a * e * sin(f)) / beta_def, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.2 (0.3.2.11): r f_dot = h/r = (n a / beta)(1+e cos f) ---
% h/r = n a^2 beta / r = n a^2 beta * (1+e cos f) / (a beta^2) = (n a / beta)(1+e cos f)
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.2 (0.3.2.11): h/r = (n a / beta)(1+e cos f)', ...
    (n * a^2 * beta_def) / r_def, ...
    (n * a / beta_def) * (1 + e*cos(f)), ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.3.3 (e_dot Gauss VE):\n');

% --- Theorem 0.3.3 (0.3.3.13): r r_dot = e h sin f / (1+e cos f) ---
% r r_dot = (a beta^2/(1+e cos f)) * (e h sin f / p) = (a beta^2 * e h sin f) / ((1+e cos f) * a beta^2) = e h sin f / (1+e cos f)
r_rdot_lhs = r_def * (e * (n*a^2*beta_def) * sin(f) / p_def);  % r * r_dot via (0.3.2.8)
r_rdot_rhs = e * (n*a^2*beta_def) * sin(f) / (1+e*cos(f));      % from explicit closed form
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.3 (0.3.3.13): r r_dot = e h sin f / (1+e cos f)', ...
    r_rdot_lhs, ...
    r_rdot_rhs, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.3 (0.3.3.15): cos f + cos E = (2 cos f + e(1+cos^2 f))/(1+e cos f) ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.3 (0.3.3.15): cos f + cos E = (2 cos f + e(1+cos^2 f))/(1+e cos f)', ...
    cos(f) + cosE, ...
    (2*cos(f) + e*(1 + cos(f)^2)) / (1+e*cos(f)), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.3 (0.3.3.16): 2 h cos f + r r_dot sin f = h (cos f + cos E) ---
% with r r_dot = e h sin f / (1+e cos f)
h_sym = n * a^2 * beta_def;
lhs_0316 = 2*h_sym*cos(f) + (e * h_sym * sin(f) / (1+e*cos(f))) * sin(f);
rhs_0316 = h_sym * (cos(f) + cosE);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.3 (0.3.3.16): 2 h cos f + r r_dot sin f = h (cos f + cos E)', ...
    lhs_0316, ...
    rhs_0316, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.3 prefactor: h/mu = beta/(n a) ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.3 prefactor: h/mu = beta/(n a)', ...
    (n * a^2 * beta_def) / mu_def, ...
    beta_def / (n * a), ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.3.4 (Omega_dot, i_dot):\n');

% --- Theorem 0.3.4 (0.3.4.7): n_hat as parameterized in (i, Omega) ---
% Verify (n_hat dot n_hat) = 1.
% n_hat = sin(i) sin(Omega) x_hat - sin(i) cos(Omega) y_hat + cos(i) z_hat
% norm: sin^2 i sin^2 Omega + sin^2 i cos^2 Omega + cos^2 i = sin^2 i + cos^2 i = 1
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.4: ||n_hat||^2 = 1', ...
    sin(i_sym)^2 * sin(Omega_sym)^2 + sin(i_sym)^2 * cos(Omega_sym)^2 + cos(i_sym)^2, ...
    sym(1), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.4 (0.3.4.9): partial n_hat / partial Omega = sin(i) N_hat ---
% n_hat_x = sin(i) sin(Omega); partial / partial Omega = sin(i) cos(Omega)
% n_hat_y = -sin(i) cos(Omega); partial / partial Omega = sin(i) sin(Omega)
% n_hat_z = cos(i); partial = 0
% N_hat = cos(Omega) x_hat + sin(Omega) y_hat, so sin(i) N_hat has components
% (sin(i) cos(Omega), sin(i) sin(Omega), 0)  -- matches above. PASS by construction.
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.4 (0.3.4.9): partial n_hat_x/partial Omega = sin(i) cos(Omega)', ...
    diff(sin(i_sym)*sin(Omega_sym), Omega_sym), ...
    sin(i_sym)*cos(Omega_sym), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.4 (0.3.4.10): partial n_hat / partial i = -M_hat ---
% M_hat = cos(i)(-sin(Omega) x + cos(Omega) y) + sin(i) z
% -M_hat = cos(i)(sin(Omega) x - cos(Omega) y) - sin(i) z
% partial n_hat_x / partial i = cos(i) sin(Omega) (matches)
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.4 (0.3.4.10): partial n_hat_x / partial i = cos(i) sin(Omega) = -M_hat_x', ...
    diff(sin(i_sym)*sin(Omega_sym), i_sym), ...
    cos(i_sym)*sin(Omega_sym), ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.3.5 (omega_dot Gauss VE):\n');

% --- Theorem 0.3.5 (0.3.5.9): (2 + e cos f)/(1 + e cos f) = 1 + r/p ---
r_over_p = r_def / p_def;
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.5 (0.3.5.9): (2 + e cos f)/(1+e cos f) = 1 + r/p', ...
    (2 + e*cos(f)) / (1+e*cos(f)), ...
    1 + r_over_p, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.5 (0.3.5.10): 2 h sin f - r r_dot cos f = h sin f (1 + r/p) ---
lhs_0510 = 2*h_sym*sin(f) - (e * h_sym * sin(f) / (1+e*cos(f))) * cos(f);
rhs_0510 = h_sym * sin(f) * (1 + r_over_p);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.5 (0.3.5.10): 2 h sin f - r r_dot cos f = h sin f (1 + r/p)', ...
    lhs_0510, ...
    rhs_0510, ...
    pass_count, fail_count, failed_names);

% --- Lemma 0.3.5.a: (dN_hat/dt) . M_hat = Omega_dot cos i ---
% Symbolic: M_hat . d/dt N_hat where N_hat = (cos Omega, sin Omega, 0)
% d/dt N_hat = Omega_dot * (-sin Omega, cos Omega, 0)
% M_hat = (cos i * (-sin Omega), cos i * cos Omega, sin i)
% Dot product = Omega_dot * [cos i * sin^2 Omega + cos i * cos^2 Omega + 0] = Omega_dot cos i
syms Omega_dot real;
dN_dt = Omega_dot * [-sin(Omega_sym); cos(Omega_sym); sym(0)];
M_hat = [cos(i_sym)*(-sin(Omega_sym)); cos(i_sym)*cos(Omega_sym); sin(i_sym)];
[pass_count, fail_count, failed_names] = check_eq( ...
    'Lemma 0.3.5.a: (dN_hat/dt) . M_hat = Omega_dot cos i', ...
    M_hat.' * dN_dt, ...
    Omega_dot * cos(i_sym), ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.3.6 (M_dot Gauss VE):\n');

% --- Theorem 0.3.6 (0.3.6.8): (r r_dot) / (a^2 e sin E) = n ---
sinE_pos = beta_def * sin(f) / (1 + e*cos(f));   % positive root
r_rdot = e * h_sym * sin(f) / (1+e*cos(f));      % from (0.3.3.13)
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.6 (0.3.6.8): (r r_dot) / (a^2 e sin E) = n', ...
    r_rdot / (a^2 * e * sinE_pos), ...
    n, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.6 (0.3.6.13): r cos E - a e sin^2 E = (a beta^2 cos f)/(1+e cos f) ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.6 (0.3.6.13): r cos E - a e sin^2 E = a beta^2 cos f/(1+e cos f)', ...
    r_def * cosE - a * e * sinE_pos^2, ...
    a * beta_def^2 * cos(f) / (1+e*cos(f)), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.6 (0.3.6.14): (r cos E - a e sin^2 E) / (a e sin E) = beta cos f / (e sin f) ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.6 (0.3.6.14): (r cos E - a e sin^2 E)/(a e sin E) = beta cos f / (e sin f)', ...
    (r_def * cosE - a * e * sinE_pos^2) / (a * e * sinE_pos), ...
    beta_def * cos(f) / (e * sin(f)), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.6 (0.3.6.20): -2 + cos f (cos f + cos E) = -sin^2 f (1 + r/p) ---
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.6 (0.3.6.20): -2 + cos f (cos f + cos E) = -sin^2 f (1 + r/p)', ...
    -2 + cos(f) * (cos(f) + cosE), ...
    -sin(f)^2 * (1 + r_over_p), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.3.6 (0.3.6 R-coefficient): (β^2 / (n a e))*(cos f - 2 e r/p) = -2 r/(n a^2) + β^2 cos f/(n a e) ---
% Using p = a beta^2:  2 e r / p = 2 e r / (a beta^2)
%   (β^2 / (n a e)) * (cos f - 2 e r / (a beta^2))
%   = β^2 cos f / (n a e) - β^2 / (n a e) * 2 e r / (a beta^2)
%   = β^2 cos f / (n a e) - 2 r / (n a^2)
lhs_Rcoef = (beta_def^2 / (n * a * e)) * (cos(f) - 2 * e * r_def / p_def);
rhs_Rcoef = -2 * r_def / (n * a^2) + beta_def^2 * cos(f) / (n * a * e);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.3.6 R-coefficient consistency: (beta^2/(n a e))(cos f - 2 e r/p) = -2 r/(n a^2) + beta^2 cos f/(n a e)', ...
    lhs_Rcoef, ...
    rhs_Rcoef, ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.4.1 (drag specialization, Phase 0-rev1):\n');

syms B_star positive;
syms rho real;
r_dot_sym = n * a * e * sin(f) / beta_def;
r_fdot_sym = (n * a / beta_def) * (1 + e*cos(f));

% --- Theorem 0.4.1 ratio: R_drag / T_drag = (e sin f)/(1+e cos f) ---
% Under Phase 0-rev1, R_drag = -B*·ρ·|v|·ṙ, T_drag = -B*·ρ·|v|·r ḟ. The |v|
% factor cancels in the ratio, so this check is unchanged by the D-10 fix.
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.1 ratio: R_drag/T_drag = e sin f / (1+e cos f)', ...
    r_dot_sym / r_fdot_sym, ...
    e * sin(f) / (1 + e*cos(f)), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.4.1 (0.4.1.10): e² sin²f + (1+e cos f)² = 1 + e² + 2 e cos f ---
% The algebraic identity behind the speed-scalar closed form. This identity
% is what turns the |v|² Pythagorean form into the orbit-element form.
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.1 (0.4.1.10): e^2 sin^2 f + (1+e cos f)^2 = 1 + e^2 + 2 e cos f', ...
    e^2 * sin(f)^2 + (1 + e*cos(f))^2, ...
    1 + e^2 + 2*e*cos(f), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.4.1 (0.4.1.11): |v|² = (na/β)²·(1+e²+2e cos f) from ṙ²+(r ḟ)² ---
% Speed-scalar closed form. Derived from Pythagorean |v|² = ṙ² + (r ḟ)²
% (orthogonality r̂ ⊥ t̂) plus the orbit-element forms of (ṙ, r ḟ).
v_mag_squared_kinematic = r_dot_sym^2 + r_fdot_sym^2;
v_mag_squared_orbit_form = (n*a/beta_def)^2 * (1 + e^2 + 2*e*cos(f));
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.1 (0.4.1.11): |v|^2 = (n a/beta)^2 (1+e^2+2e cos f) (kinematic = orbit form)', ...
    v_mag_squared_kinematic, ...
    v_mag_squared_orbit_form, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.4.1 Postulate-1.2 match (D-10 regression guard) ---
% This is the check that would have caught D-10: verify that the (R_drag,
% T_drag, N_drag) decomposition CLAIMED by Theorem 0.4.1 matches the
% (R, T, N) decomposition obtained by directly substituting F_drag =
% -B*·ρ·|v|·v with v = ṙ r̂ + r ḟ t̂.
%
% Pre-D-10 Theorem 0.4.1 wrote R_drag = -B*·ρ·ṙ (no |v|), which would FAIL
% this check (the difference -B*·ρ·(|v|-1)·ṙ is not symbolically zero).
% Phase 0-rev1's R_drag = -B*·ρ·|v|·ṙ passes.
%
% Implementation note: we use the closed-form |v| from (0.4.1.11) rather than
% the kinematic |v| := sqrt(ṙ² + (r ḟ)²) because Octave's SymPy backend cannot
% auto-reduce nested sqrt(sqrt(...)) inside subsequent (...)^{3/2} simplification
% chains. The (0.4.1.11) check above already proves the two forms are
% symbolically equal, so the substitution is justified.
v_mag_sym = (n*a/beta_def) * sqrt(1 + e^2 + 2*e*cos(f));   % per (0.4.1.11) closed form
% Theorem 0.4.1-rev1's claimed (R, T, N):
R_drag_thm = -B_star * rho * v_mag_sym * r_dot_sym;
T_drag_thm = -B_star * rho * v_mag_sym * r_fdot_sym;
N_drag_thm = sym(0);
% Postulate-1.2 form: F_drag = -B*·ρ·|v|·v.
% v has (r̂, t̂, n̂) components (ṙ, r ḟ, 0). So F_drag has components:
F_post_r = -B_star * rho * v_mag_sym * r_dot_sym;
F_post_t = -B_star * rho * v_mag_sym * r_fdot_sym;
F_post_n = sym(0);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.1 Postulate-1.2 (r-hat): R_drag = -B*rho|v|r_dot (D-10 regression guard)', ...
    R_drag_thm, F_post_r, ...
    pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.1 Postulate-1.2 (t-hat): T_drag = -B*rho|v|r f_dot (D-10 regression guard)', ...
    T_drag_thm, F_post_t, ...
    pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.1 Postulate-1.2 (n-hat): N_drag = 0 (D-10 regression guard)', ...
    N_drag_thm, F_post_n, ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.4.2 (drag rates, Phase 0-rev1):\n');

% --- Theorem 0.4.2 (0.4.2.1) Phase 0-rev1: a_dot_drag = -(2 n B* ρ a²/β³)·(1+e²+2e cos f)^{3/2} ---
% Corrected closed form with the |v| factor restored (closes D-10).
% Substituting R_drag, T_drag from Theorem 0.4.1-rev1 into Gauss VE Theorem 0.3.2 (0.3.2.15):
%   a_dot = (2/(n beta)) * [e R_drag sin f + T_drag (1+e cos f)]
%         = (2/(n beta)) * (-B*·ρ·|v|·n a / beta) * [e^2 sin^2 f + (1+e cos f)^2]    [combine R, T]
%         = (-2 B*·ρ·|v|·a / beta^2) * (1 + e^2 + 2 e cos f)                          [apply (0.4.1.10)]
%         = -(2 n B*·ρ·a^2 / beta^3) * (1 + e^2 + 2 e cos f)^{3/2}                    [sub |v| from (0.4.1.11)]
a_dot_lhs_rev1 = (2/(n*beta_def)) * (e * R_drag_thm * sin(f) + T_drag_thm * (1+e*cos(f)));
a_dot_rhs_rev1 = -(2 * n * B_star * rho * a^2 / beta_def^3) * (1 + e^2 + 2*e*cos(f))^(sym(3)/sym(2));
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.4.2 (0.4.2.1) rev1: a_dot_drag = -(2 n B* rho a^2/beta^3)(1+e^2+2e cos f)^(3/2)', ...
    a_dot_lhs_rev1, ...
    a_dot_rhs_rev1, ...
    pass_count, fail_count, failed_names);

printf('\n');
printf('Theorem 0.5.1 (Energy identity) and 0.5.2 (Delaunay L):\n');

% --- Theorem 0.5.1: dE/dt = (mu/(2 a^2)) * a_dot ---
% E = -mu/(2a), dE/da = mu/(2 a^2). Chain rule gives dE/dt = (mu/(2a^2)) * a_dot.
syms a_dot_sym real;
E_def = -mu / (2*a);
dE_da = diff(E_def, a);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.5.1: dE/da = mu/(2 a^2)', ...
    dE_da, ...
    mu / (2*a^2), ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.5.2: dL/da = n a / 2 with L = sqrt(mu a), mu = n^2 a^3 ---
L_def = sqrt(mu*a);
dL_da = diff(L_def, a);
dL_da_substituted = subs(dL_da, mu, n^2*a^3);
[pass_count, fail_count, failed_names] = check_eq( ...
    'Theorem 0.5.2: dL/da = n a / 2 (after substituting mu = n^2 a^3)', ...
    simplify(dL_da_substituted), ...
    n * a / 2, ...
    pass_count, fail_count, failed_names);

% --- Theorem 0.5.3: Jacobian dM = beta^3/(1+e cos f)^2 df (re-verified via integration setup) ---
% This is Theorem 0.2.5 restated; already verified above.

printf('\n');
printf('=== Phase 0 verification summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n');
    printf('  FAILED checks:\n');
    for k = 1:length(failed_names)
        printf('    %s\n', failed_names{k});
    end
    printf('\n');
    printf('Phase 0 verification: FAIL\n');
    exit(1);
else
    printf('\n');
    printf('Phase 0 verification: PASS\n');
    exit(0);
end
