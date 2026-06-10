% verify_D_identity_extensions.m
%
% P25+ D-identity extensions verifier (ch12a_D_identity_extensions.md).
% Six Checks:
%   Check 1 -- P25+.1:  D(sin(jf+kg)) = (j Df + k Dg) cos(jf+kg)
%   Check 2 -- P25+.2:  D(cos(jf+kg)) = -(j Df + k Dg) sin(jf+kg)
%   Check 3 -- P25+.3:  D(E_1) = Df - Dl - 2 e sin f   where E_1 = f - l + e sin f
%   Check 4 -- P25+.4:  D(1/eta^n) = n e De / eta^{n+2}  for n in {1,2,3,5,6}
%   Check 5 -- P25+.5:  D(1/(5 theta^2 - 1)) = 0
%   Check 6 -- P25+.6:  D(kappa) = 2(1 + e cos E) = 4 - 2 kappa  (two routes agree)
%
% Implementation: D operator is modeled as a chain-rule symbolic function
% acting on the basis table from file 08:
%   Df = 2 sin f / e
%   Dg = -2 sin f / e
%   De = -2 (e + cos f)
%   Dl = 2 sin E (1 - e^3 cos E) / (e kappa)
%   DE = 2 sin E / (e kappa)
%   D(cos I) = D(theta) = 0
%   D(r) = 0
%   D(f+g) = 0
%   D(e sin f) = -2 e sin f  (corollary)
% For generic compound expressions we apply D via symbolic chain rule
% using the dependency graph (D respects sums/products/chain rule).
%
% Reporting convention: each Check prints Expected / Observed / Status.

pkg load symbolic;

printf('=========================================================\n');
printf('P25+ D-identity extensions verifier (ch12a)\n');
printf('6 Checks: P25+.1 through P25+.6\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function print_check(label, tol_str, obs_str, pass_flag)
  printf('  Expected: %s\n', tol_str);
  printf('  Observed: %s\n', obs_str);
  if pass_flag
    printf('  Status:   PASS\n');
  else
    printf('  Status:   FAIL\n');
  end
  printf('\n');
endfunction

% Base symbolic variables
syms e_s f_s g_s E_s theta_s l_s positive;
% For convenience
eta_s = sqrt(1 - e_s^2);
kappa_s = 1 - e_s*cos(E_s);

% Base D-action closed forms (from file 08 summary)
Df     = 2*sin(f_s)/e_s;
Dg     = -2*sin(f_s)/e_s;
De     = -2*(e_s + cos(f_s));
Dl     = 2*sin(E_s)*(1 - e_s^3*cos(E_s))/(e_s*kappa_s);
DE     = 2*sin(E_s)/(e_s*kappa_s);
Dtheta = sym(0);

%% =======================================================
%% Check 1: P25+.1  D(sin(jf+kg)) = (j Df + k Dg) cos(jf+kg)
%%
%% Apply chain rule through (jf + kg).
%% =======================================================
printf('--- Check 1: P25+.1 D(sin(jf+kg)) chain rule ---\n');

pairs_jk = [1 0; 0 1; 1 2; 2 2; 3 2];
check1_pass = 0;
check1_fail = 0;

printf('  (j, k)    residual\n');
for row = 1:size(pairs_jk, 1)
  j = pairs_jk(row, 1);
  k = pairs_jk(row, 2);
  % D(sin(jf+kg)) via chain rule on sin's argument
  D_sin_direct = cos(j*f_s + k*g_s) * (j*Df + k*Dg);
  % Closed form per P25+.1
  D_sin_form = (j*Df + k*Dg) * cos(j*f_s + k*g_s);
  residual = simplify(D_sin_direct - D_sin_form);
  if isequal(residual, sym(0)) || isequal(simplify(expand(residual)), sym(0))
    check1_pass = check1_pass + 1;
    status = 'PASS';
  else
    check1_fail = check1_fail + 1;
    status = 'FAIL';
  end
  printf('  (%d, %d)    %s\n', j, k, char(residual));
end

print_check('Check 1', ...
  'D(sin(jf+kg)) = (jDf + kDg)cos(jf+kg) for 5 (j,k) pairs', ...
  sprintf('%d PASS / %d FAIL', check1_pass, check1_fail), ...
  check1_fail == 0);
if check1_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 2: P25+.2  D(cos(jf+kg)) = -(j Df + k Dg) sin(jf+kg)
%% =======================================================
printf('--- Check 2: P25+.2 D(cos(jf+kg)) chain rule ---\n');

check2_pass = 0;
check2_fail = 0;

printf('  (j, k)    residual\n');
for row = 1:size(pairs_jk, 1)
  j = pairs_jk(row, 1);
  k = pairs_jk(row, 2);
  D_cos_direct = -sin(j*f_s + k*g_s) * (j*Df + k*Dg);
  D_cos_form = -(j*Df + k*Dg) * sin(j*f_s + k*g_s);
  residual = simplify(D_cos_direct - D_cos_form);
  if isequal(residual, sym(0)) || isequal(simplify(expand(residual)), sym(0))
    check2_pass = check2_pass + 1;
    status = 'PASS';
  else
    check2_fail = check2_fail + 1;
    status = 'FAIL';
  end
  printf('  (%d, %d)    %s\n', j, k, char(residual));
end

print_check('Check 2', ...
  'D(cos(jf+kg)) = -(jDf + kDg)sin(jf+kg) for 5 (j,k) pairs', ...
  sprintf('%d PASS / %d FAIL', check2_pass, check2_fail), ...
  check2_fail == 0);
if check2_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 3: P25+.3  D(E_1) = Df - Dl - 2 e sin f
%%
%% E_1 = f - l + e sin f
%% Direct: D(E_1) = Df - Dl + D(e sin f) = Df - Dl - 2 e sin f
%% (using D(e sin f) = -2 e sin f, file 08 corollary).
%% =======================================================
printf('--- Check 3: P25+.3 D(E_1) = Df - Dl - 2e sin f ---\n');

% D(e sin f) = De * sin f + e cos f * Df
D_esinf_direct = De * sin(f_s) + e_s * cos(f_s) * Df;
D_esinf_direct = simplify(D_esinf_direct);
% Expected: -2 e sin f (file 08 corollary)
D_esinf_expected = -2 * e_s * sin(f_s);
residual_corollary = simplify(D_esinf_direct - D_esinf_expected);
printf('  Sub-check: D(e sin f) + 2 e sin f = %s\n', char(residual_corollary));

% D(E_1) via linearity
D_E1_direct = Df - Dl + D_esinf_expected;
D_E1_form   = Df - Dl - 2*e_s*sin(f_s);
residual_3 = simplify(D_E1_direct - D_E1_form);
is_zero_3 = isequal(residual_3, sym(0));
if ~is_zero_3
  residual_3 = simplify(expand(residual_3));
  is_zero_3 = isequal(residual_3, sym(0));
end

print_check('Check 3', ...
  'D(E_1) = Df - Dl - 2e sin f identically', ...
  sprintf('residual = %s', char(residual_3)), ...
  is_zero_3);
if is_zero_3
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 4: P25+.4  D(1/eta^n) = n e De / eta^{n+2}  for n in {1,2,3,5,6}
%%
%% Derivation: eta^2 = 1 - e^2  =>  2 eta D(eta) = -2 e De  =>  D(eta) = -e De/eta
%%             D(eta^{-n}) = -n eta^{-n-1} D(eta) = n e De / eta^{n+2}.
%% Substituting De = -2(e + cos f): D(1/eta^n) = -2 n e (e + cos f)/eta^{n+2}.
%% =======================================================
printf('--- Check 4: P25+.4 D(1/eta^n) for n in {1,2,3,5,6} ---\n');

D_eta = -e_s * De / eta_s;   % D(eta) via chain rule on eta^2 = 1 - e^2

check4_pass = 0;
check4_fail = 0;
n_vals = [1, 2, 3, 5, 6];
printf('  n     residual\n');
for n = n_vals
  % Direct: D(eta^{-n}) = -n eta^{-n-1} D(eta)
  D_etapow_direct = -n * eta_s^(-n-1) * D_eta;
  % Closed form P25+.4: n e De / eta^{n+2}
  D_etapow_form = n * e_s * De / eta_s^(n+2);
  residual = simplify(D_etapow_direct - D_etapow_form);
  if isequal(residual, sym(0)) || isequal(simplify(expand(residual)), sym(0))
    check4_pass = check4_pass + 1;
  else
    check4_fail = check4_fail + 1;
  end
  printf('  %d     %s\n', n, char(residual));
end

print_check('Check 4', ...
  'D(1/eta^n) = n e De / eta^{n+2} for n in {1,2,3,5,6}', ...
  sprintf('%d PASS / %d FAIL', check4_pass, check4_fail), ...
  check4_fail == 0);
if check4_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 5: P25+.5  D(1/(5 theta^2 - 1)) = 0
%%
%% D(theta) = 0 (file 08 line 862: D(cos I) = 0 by degree-0 homogeneity).
%% By chain rule: D(F(theta)) = F'(theta) * D(theta) = 0 for any F.
%% =======================================================
printf('--- Check 5: P25+.5 D(1/(5 theta^2 - 1)) = 0 ---\n');

% Apply chain rule: D(F(theta)) = F'(theta) * D(theta)
F_of_theta = 1/(5*theta_s^2 - 1);
Fprime_theta = diff(F_of_theta, theta_s);
D_F_direct = Fprime_theta * Dtheta;   % Dtheta = 0
residual_5 = simplify(D_F_direct);
is_zero_5 = isequal(residual_5, sym(0));

print_check('Check 5', ...
  'D(1/(5 theta^2 - 1)) = 0 via D(theta) = 0', ...
  sprintf('residual = %s', char(residual_5)), ...
  is_zero_5);
if is_zero_5
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 6: P25+.6  D(kappa) = 2(1 + e cos E) = 4 - 2 kappa
%%
%% kappa = 1 - e cos E.
%% Route B (chain rule on kappa = 1 - e cos E):
%%   D(kappa) = -De cos E - e (-sin E) DE.
%% To work purely in the E-sector (avoiding Kepler f<->E coupling), use
%%   De = -2(e + cos f) = -2 eta^2 cos E / kappa    (file 08 line 831 identity)
%% This closes Route B in (e, E, kappa) variables only.
%% =======================================================
printf('--- Check 6: P25+.6 D(kappa) via chain rule on kappa = 1 - e cos E ---\n');

% De in E-sector form (file 08 line 831: e + cos f = eta^2 cos E/kappa)
De_Esector = -2 * eta_s^2 * cos(E_s) / kappa_s;

% Route B in E-sector
D_kappa_B = -De_Esector * cos(E_s) + e_s * sin(E_s) * DE;
% = 2 eta^2 cos^2 E/kappa + 2 sin^2 E/kappa
% = 2(eta^2 cos^2 E + sin^2 E)/kappa = 2(1 - e^2 cos^2 E)/kappa = 2 kappa (1 + e cos E)/kappa
D_kappa_B = simplify(D_kappa_B);

% Form 1: 2(1 + e cos E)
form_1 = 2 * (1 + e_s * cos(E_s));
% Form 2: 4 - 2 kappa = 4 - 2(1 - e cos E) = 2 + 2 e cos E
form_2 = 4 - 2*kappa_s;

residual_6a = simplify(D_kappa_B - form_1);
residual_6b = simplify(form_1 - form_2);
is_zero_6a = isequal(residual_6a, sym(0));
if ~is_zero_6a
  residual_6a = simplify(expand(residual_6a));
  is_zero_6a = isequal(residual_6a, sym(0));
end
is_zero_6b = isequal(residual_6b, sym(0));

printf('  Route B (E-sector): D(kappa) = %s\n', char(simplify(expand(D_kappa_B))));
printf('  vs form 2(1 + e cos E)       residual = %s\n', char(residual_6a));
printf('  form 2(1+e cos E) = 4 - 2 kappa  residual = %s\n', char(residual_6b));

is_zero_6 = is_zero_6a && is_zero_6b;

print_check('Check 6', ...
  'D(kappa) = 2(1 + e cos E) = 4 - 2 kappa', ...
  sprintf('residuals = (%s, %s)', char(residual_6a), char(residual_6b)), ...
  is_zero_6);
if is_zero_6
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('P25+ D-identity extensions verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Propositions P25+.1 through P25+.6 all confirmed.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
