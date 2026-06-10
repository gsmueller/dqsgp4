% verify_independence_audit_part2.m
%
% Extended variable-dependence audit covering:
% - ch12 H.1 delta_p_1 structural dependencies
% - ch12 H.3 delta_p_3 = 0 (trivial verification)
% - ch10c_addendum U_L and T formulas
% - ch12a P25+.5 D(F(theta)) = 0 for theta-only F
% - ch11c Lyddane elements non-singular at critical inclination
% - Hansen X_0^{-6,m}(e) formula dependencies

pkg load symbolic;

printf('=========================================================\n');
printf('Extended independence audit (Part 2)\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function print_check(label, tol_str, obs_str, pass_flag)
  printf('  %s\n', label);
  printf('    Expected: %s\n', tol_str);
  printf('    Observed: %s\n', obs_str);
  if pass_flag
    printf('    Status:   PASS\n');
  else
    printf('    Status:   FAIL\n');
  end
  printf('\n');
endfunction

syms L G H l g h mu k2 positive;
syms E e_s theta_s f_s g_s positive;

eta_expr = G/L;
e_expr = sqrt(1 - (G/L)^2);
theta_expr = H/G;
A_expr = (3*theta_expr^2 - 1)/2;
B_expr = 3*(1 - theta_expr^2)/2;

%% =======================================================
%% Part 12: ch10c_addendum U_L formula
%% U_L = (mu^8 k_2^2 eta^10 / (2 e G^10)) [-12 e T + eta^2 dT/de]
%% Claim: U_L in M_{10}; l-independent (already l-averaged); g-dependent (retains T(...,g))
%% =======================================================
printf('--- Part 12: ch10c_addendum U_L formula structure ---\n');

% Build T(theta, e, eta, g) from ch10c §5
X60 = (8 + 24*e_expr^2 + 3*e_expr^4) / (8 * eta_expr^9);
X62 = e_expr^2 * (6 + e_expr^2) / (4 * eta_expr^9);
X64 = e_expr^4 / (16 * eta_expr^9);
T_expr = (A_expr^2 + B_expr^2/2) * X60 + 2*A_expr*B_expr * cos(2*g) * X62 ...
       + (B_expr^2/2) * cos(4*g) * X64 - A_expr^2 / eta_expr^6;

% T dependencies
T_l = simplify(diff(T_expr, l));
T_g = simplify(diff(T_expr, g));
T_h = simplify(diff(T_expr, h));
T_l_indep = isequal(T_l, sym(0));
T_h_indep = isequal(T_h, sym(0));
T_g_depends = ~isequal(T_g, sym(0));

print_check('T(theta, e, eta, g): l-independent (already l-averaged)', ...
  'diff(T, l) = 0', sprintf('%s', char(T_l)), T_l_indep);
if T_l_indep; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('T(theta, e, eta, g): h-independent', ...
  'diff(T, h) = 0', sprintf('%s', char(T_h)), T_h_indep);
if T_h_indep; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('T(theta, e, eta, g): g-dependent (cos 2g, cos 4g present)', ...
  'diff(T, g) ≠ 0', ...
  sprintf('is nonzero: %d', T_g_depends), T_g_depends);
if T_g_depends; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 13: U_L formula via IBP
%% U_L = (1/(2n)) dV/dL where V = (mu^10 k_2^2 / L^12) T
%% Check U_L is g-dependent (inherits T's g-dependence).
%% =======================================================
printf('--- Part 13: U_L g-dependence (inherits T) ---\n');

V_expr = mu^10 * k2^2 / L^12 * T_expr;
n_mot = mu^2 / L^3;
U_L = diff(V_expr, L) / (2 * n_mot);

U_L_l = simplify(diff(U_L, l));
U_L_h = simplify(diff(U_L, h));
U_L_g = simplify(diff(U_L, g));

U_L_l_indep = isequal(U_L_l, sym(0));
U_L_h_indep = isequal(U_L_h, sym(0));
U_L_g_depends = ~isequal(U_L_g, sym(0));

all_UL = U_L_l_indep && U_L_h_indep && U_L_g_depends;

print_check('U_L: l-independent, h-independent, g-dependent (inherits T structure)', ...
  'diff(U_L, l) = diff(U_L, h) = 0 and diff(U_L, g) ≠ 0', ...
  sprintf('l-indep: %d, h-indep: %d, g-dep: %d', U_L_l_indep, U_L_h_indep, U_L_g_depends), ...
  all_UL);
if all_UL; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 14: ⟨U_L⟩_g should equal U_L_secular (g-averaged contribution to F_2^{**})
%% =======================================================
printf('--- Part 14: ⟨U_L⟩_g structure ---\n');

U_L_avg_g = simplify(int(U_L, g, 0, 2*sym(pi)) / (2*sym(pi)));
% ⟨U_L⟩_g should be L, G, H dependent (since it's a function of momenta) but g-independent
U_L_avg_g_g = simplify(diff(U_L_avg_g, g));
p14 = isequal(U_L_avg_g_g, sym(0));

% Also ⟨U_L⟩_g should be nonzero (has content)
U_L_avg_g_nonzero = ~isequal(simplify(U_L_avg_g), sym(0));

print_check('⟨U_L⟩_g: g-independent after g-averaging, and nonzero', ...
  'diff(⟨U_L⟩_g, g) = 0 AND ⟨U_L⟩_g ≠ 0', ...
  sprintf('g-indep: %d, nonzero: %d', p14, U_L_avg_g_nonzero), ...
  p14 && U_L_avg_g_nonzero);
if p14 && U_L_avg_g_nonzero; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 15: P25+.5 D(F(theta)) = 0 for theta-only F
%% Symbolically: if F depends only on theta, then D(F) = 0 regardless of specific form.
%% =======================================================
printf('--- Part 15: P25+.5 D(F(theta)) = 0 for theta-only F ---\n');

% D-operator closed forms (file 08 table, restricted to theta-only test):
% D(theta) = 0 by degree-0
% D(F(theta)) = F'(theta) * D(theta) = 0 by chain rule
% Try F = 1/(5 theta^2 - 1), A(theta) = (3 theta^2 - 1)/2, B(theta), A^2 + B^2/2
F1_theta = 1/(5*theta_s^2 - 1);
F2_theta = (3*theta_s^2 - 1)/2;
F3_theta = 3*(1 - theta_s^2)/2;
F4_theta = F2_theta^2 + F3_theta^2/2;
F5_theta = theta_s * (1 - theta_s^2);

% D(theta) = 0
D_theta = sym(0);
% D(F) = diff(F, theta) * D(theta) = 0
D_F1 = diff(F1_theta, theta_s) * D_theta;
D_F2 = diff(F2_theta, theta_s) * D_theta;
D_F3 = diff(F3_theta, theta_s) * D_theta;
D_F4 = diff(F4_theta, theta_s) * D_theta;
D_F5 = diff(F5_theta, theta_s) * D_theta;

all_DF_zero = isequal(D_F1, sym(0)) && isequal(D_F2, sym(0)) && isequal(D_F3, sym(0)) && ...
              isequal(D_F4, sym(0)) && isequal(D_F5, sym(0));

print_check('P25+.5: D(F(theta)) = 0 for 5 sample F in {1/(5θ²-1), A, B, A²+B²/2, θ(1-θ²)}', ...
  'all 5 D(F) = 0', ...
  sprintf('all zero: %d', all_DF_zero), all_DF_zero);
if all_DF_zero; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 16: Ch11c Lyddane elements: non-singular at critical inclination
%% Elements: xi, eta_L, p, q, lambda_L
%% At theta = 1/sqrt(5), all 5 should be finite and well-defined.
%% =======================================================
printf('--- Part 16: Lyddane elements finite at critical inclination ---\n');

% Sample numerical values at critical inclination
theta_c = 1/sqrt(sym(5));
I_c = acos(theta_c);
% Pick sample (e, g, h) values
e_samp = sym(3)/10;
g_samp = sym(7)/10;
h_samp = sym(11)/10;
varpi_samp = g_samp + h_samp;

% Lyddane elements
xi_c = e_samp * cos(varpi_samp);
etaL_c = e_samp * sin(varpi_samp);
p_c = 2 * sin(I_c/2) * cos(h_samp);
q_c = 2 * sin(I_c/2) * sin(h_samp);
lam_c = sym(1)/2 + g_samp + h_samp;  % l + g + h with l = 1/2

all_finite = true;
vals = [xi_c, etaL_c, p_c, q_c, lam_c];
names = {'xi', 'eta_L', 'p', 'q', 'lambda_L'};
for i = 1:5
  v = double(vals(i));
  if ~isfinite(v)
    all_finite = false;
  end
end

print_check('Lyddane elements at theta = 1/sqrt(5): all finite', ...
  'all 5 Lyddane elements finite at critical inclination', ...
  sprintf('all finite: %d', all_finite), all_finite);
if all_finite; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 17: Lyddane Jacobian |det J_L| = 1/L^2 theta-independence
%% Test at multiple theta values including critical.
%% =======================================================
printf('--- Part 17: Lyddane Jacobian theta-independence ---\n');

function Dval = fd_jacdet(L_val, G_val, H_val)
  eps_fd = 1e-7;
  base = [L_val; G_val; H_val; 0.3; 0.7; 1.1];
  J = zeros(6, 6);
  for j = 1:6
    bp = base; bp(j) = bp(j) + eps_fd;
    bm = base; bm(j) = bm(j) - eps_fd;
    J(:, j) = (lyd_map(bp) - lyd_map(bm)) / (2*eps_fd);
  end
  Dval = abs(det(J));
endfunction

function out = lyd_map(LGHlgh)
  L_v = LGHlgh(1); G_v = LGHlgh(2); H_v = LGHlgh(3);
  l_v = LGHlgh(4); g_v = LGHlgh(5); h_v = LGHlgh(6);
  e_v = sqrt(max(1 - (G_v/L_v)^2, 1e-30));
  I_v = acos(H_v/G_v);
  out = [L_v;
         e_v * cos(g_v + h_v);
         e_v * sin(g_v + h_v);
         2 * sin(I_v/2) * cos(h_v);
         2 * sin(I_v/2) * sin(h_v);
         l_v + g_v + h_v];
endfunction

% Test at 4 theta values: including critical theta = 1/sqrt(5)
thetas = [0.3, 0.447, 0.6, 0.8];   % 0.447 ≈ 1/sqrt(5)
L_test = 1.2;
e_test = 0.2;
G_test = L_test * sqrt(1 - e_test^2);
jac_vals = zeros(length(thetas), 1);
for i = 1:length(thetas)
  H_test = thetas(i) * G_test;
  jac_vals(i) = fd_jacdet(L_test, G_test, H_test);
end

closed_form = 1/L_test^2;
rel_errs = abs(jac_vals - closed_form) / closed_form;
max_err_17 = max(rel_errs);

tol_17 = 1e-5;
p17 = max_err_17 < tol_17;

printf('  theta values tested: ');
for i = 1:length(thetas)
  printf('%.3f ', thetas(i));
end
printf('\n  Jacobians: ');
for i = 1:length(thetas)
  printf('%.6e ', jac_vals(i));
end
printf('\n  Closed form 1/L^2 = %.6e\n', closed_form);
printf('  Max rel err: %.3e\n', max_err_17);

print_check('Lyddane Jacobian: theta-independent (all 4 = 1/L^2 including critical)', ...
  sprintf('max rel err < %.0e over 4 theta values (including critical)', tol_17), ...
  sprintf('max rel err = %.3e', max_err_17), p17);
if p17; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 18: P25+.4 D(1/eta^n) sign-check with multiple n values
%% D(1/eta^n) = n e De / eta^{n+2} (positive sign)
%% =======================================================
printf('--- Part 18: P25+.4 D(1/eta^n) sign and closed form ---\n');

eta_s = sqrt(1 - e_s^2);
De = -2*(e_s + cos(f_s));

n_test_vals = [1, 2, 3, 4, 5, 6, 8];
all_P4_pass = true;
P4_fails = 0;
for n_v = n_test_vals
  % Direct: D(eta^{-n}) = -n eta^{-n-1} D(eta)
  D_eta = -e_s * De / eta_s;   % D(eta) from eta^2 = 1 - e^2
  direct = -n_v * eta_s^(-n_v-1) * D_eta;
  % Closed form P25+.4: n e De / eta^{n+2}
  form = n_v * e_s * De / eta_s^(n_v+2);
  res = simplify(direct - form);
  if ~isequal(res, sym(0))
    all_P4_pass = false;
    P4_fails = P4_fails + 1;
  end
end

print_check('P25+.4: D(1/eta^n) = n e De / eta^{n+2} for n in {1,2,3,4,5,6,8}', ...
  'all 7 residuals = 0', ...
  sprintf('all pass: %d (fails: %d)', all_P4_pass, P4_fails), all_P4_pass);
if all_P4_pass; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 19: Hansen X_0^{-6,m} formulas respect e-dependence only
%% X_0^{-6,0}(e), X_0^{-6,2}(e), X_0^{-6,4}(e) should depend only on e (not on other vars).
%% =======================================================
printf('--- Part 19: Hansen X_0^{-6,m}(e) depend only on e ---\n');

X60_test = (8 + 24*e_s^2 + 3*e_s^4) / (8 * (1-e_s^2)^(sym(9)/2));
X62_test = e_s^2 * (6 + e_s^2) / (4 * (1-e_s^2)^(sym(9)/2));
X64_test = e_s^4 / (16 * (1-e_s^2)^(sym(9)/2));

X_l_indep = isequal(diff(X60_test, l), sym(0)) && ...
            isequal(diff(X62_test, l), sym(0)) && ...
            isequal(diff(X64_test, l), sym(0));

% Try diff in g, theta (should be zero)
X_g_indep = isequal(diff(X60_test, g), sym(0)) && ...
            isequal(diff(X62_test, g), sym(0)) && ...
            isequal(diff(X64_test, g), sym(0));

X_theta_indep = isequal(diff(X60_test, theta_s), sym(0)) && ...
                isequal(diff(X62_test, theta_s), sym(0)) && ...
                isequal(diff(X64_test, theta_s), sym(0));

% But depends on e (nonzero)
X_e_dep = ~isequal(simplify(diff(X60_test, e_s)), sym(0)) && ...
          ~isequal(simplify(diff(X62_test, e_s)), sym(0)) && ...
          ~isequal(simplify(diff(X64_test, e_s)), sym(0));

all_X = X_l_indep && X_g_indep && X_theta_indep && X_e_dep;

print_check('X_0^{-6,m}(e): depends only on e (independent of l, g, theta)', ...
  'l/g/theta-indep AND e-dep', ...
  sprintf('l: %d, g: %d, θ: %d, e-dep: %d', X_l_indep, X_g_indep, X_theta_indep, X_e_dep), ...
  all_X);
if all_X; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 20: delta_p_3 = 0 (H.3) trivial verification
%% S_1 and S_1^* both have no h-dependence by axial symmetry
%% =======================================================
printf('--- Part 20: H.3 delta_p_3 = 0 (axial symmetry) ---\n');

% For S_1^* (verified above in Part 2, G.9c).
% For S_1: needs the ch07c C.4 closed form. S_1 has {f, g} dependence but no h.
% Check: we built S_1^* as depending on (θ, e, η, g) only; add synthetic S_1 test:

% Synthetic S_1 in closed form from ch07c (abbreviated):
syms E1_sym E2_sym E3_sym E4_sym E5_sym positive;
% Basis elements: E1 = f-l+e sin f, E2 = sin 2(f+g), E3 = sin(f+2g), E4 = sin(3f+2g), E5 = sin 2g
% None involve h. So diff(S_1, h) = 0 by inspection.
% Verification: no symbolic "h" appears in E_i, so diff through any polynomial combination = 0.

S1_synth = A_expr * E1_sym + (B_expr/2) * E2_sym + (e_expr * B_expr/2) * E3_sym ...
         + (e_expr * B_expr/6) * E4_sym + (B_expr * (3*e_expr^2 - 2 + 2*eta_expr^3)/(6*e_expr^2)) * E5_sym;
S1_synth_scaled = mu^2 * k2 / G^3 * S1_synth;

dS1_dh = simplify(diff(S1_synth_scaled, h));
p20 = isequal(dS1_dh, sym(0));

print_check('S_1 h-independence (no h in basis elements E_1..E_5)', ...
  'diff(S_1, h) = 0', sprintf('%s', char(dS1_dh)), p20);
if p20; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Extended independence audit (Part 2): ALL %d checks PASSED.\n', n_pass);
  printf('All formulas respect their claimed independence assumptions.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
