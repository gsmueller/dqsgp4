% verify_ch06e_secular_rates.m
%
% Chapter 06e verifier for first-order secular-rate partials
% B.5.3, B.5.4, B.5.5 per CH_SECULAR_RATES_PLAN.md Part A.
%
% Each Check: symbolic SymPy diff(F_1^*, X) vs closed-form expression.
% Expected: simplify(routeA - routeB) = 0.

pkg load symbolic;

printf('=========================================================\n');
printf('Chapter 06e verifier — first-order secular rates B.5.3-5\n');
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

% Symbolic setup: F_1^* via ch06d B.5.1 closed form = μ⁴ k_2 A η³ / G⁶
% with A = (3θ²-1)/2, θ = H/G, η = G/L, e = sqrt(1 - G²/L²)
syms L_s G_s H_s mu_s k2_s positive;
eta_s = G_s / L_s;
theta_s = H_s / G_s;
A_s = (3*theta_s^2 - 1) / 2;
F1star_s = mu_s^4 * k2_s * A_s * eta_s^3 / G_s^6;

%% Check 1: B.5.3 ∂F_1*/∂L vs -3μ⁴ k_2 A η⁴/G⁷
printf('--- Check 1: B.5.3 ∂F_1*/∂L closed form ---\n');
dF1star_dL = diff(F1star_s, L_s);
B53_form = -3 * mu_s^4 * k2_s * A_s * eta_s^4 / G_s^7;
residual_1 = simplify(dF1star_dL - B53_form);
is_zero_1 = isequal(residual_1, sym(0));
if ~is_zero_1
  residual_1 = simplify(expand(residual_1));
  is_zero_1 = isequal(residual_1, sym(0));
end
print_check('Check 1', ...
  'simplify(diff(F_1^*, L) - B.5.3 form) = 0', ...
  sprintf('residual = %s', char(residual_1)), ...
  is_zero_1);
if is_zero_1
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% Check 2: B.5.4 ∂F_1*/∂G vs -3μ⁴ k_2 η³ (5θ²-1)/(2G⁷)
printf('--- Check 2: B.5.4 ∂F_1*/∂G closed form ---\n');
dF1star_dG = diff(F1star_s, G_s);
B54_form = -3 * mu_s^4 * k2_s * eta_s^3 * (5*theta_s^2 - 1) / (2 * G_s^7);
residual_2 = simplify(dF1star_dG - B54_form);
is_zero_2 = isequal(residual_2, sym(0));
if ~is_zero_2
  residual_2 = simplify(expand(residual_2));
  is_zero_2 = isequal(residual_2, sym(0));
end
print_check('Check 2', ...
  'simplify(diff(F_1^*, G) - B.5.4 form) = 0', ...
  sprintf('residual = %s', char(residual_2)), ...
  is_zero_2);
if is_zero_2
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% Check 3: B.5.5 ∂F_1*/∂H vs +3μ⁴ k_2 θ η³/G⁷
printf('--- Check 3: B.5.5 ∂F_1*/∂H closed form ---\n');
dF1star_dH = diff(F1star_s, H_s);
B55_form = 3 * mu_s^4 * k2_s * theta_s * eta_s^3 / G_s^7;
residual_3 = simplify(dF1star_dH - B55_form);
is_zero_3 = isequal(residual_3, sym(0));
if ~is_zero_3
  residual_3 = simplify(expand(residual_3));
  is_zero_3 = isequal(residual_3, sym(0));
end
print_check('Check 3', ...
  'simplify(diff(F_1^*, H) - B.5.5 form) = 0', ...
  sprintf('residual = %s', char(residual_3)), ...
  is_zero_3);
if is_zero_3
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% Summary
printf('=========================================================\n');
if n_fail == 0
  printf('Chapter 06e verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Propositions B.5.3, B.5.4, B.5.5 all confirmed.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
