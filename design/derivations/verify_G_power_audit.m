% verify_G_power_audit.m
%
% Independent numerical audit of G-exponents (ℳ_α class membership).
% Method: for S ∈ ℳ_α, S should satisfy S(λL, λG, λH) = λ^{-α} · S(L, G, H)
% (θ = H/G and e = √(1-G²/L²) are both scaling-invariant under uniform λ rescale).
%
% Extract α by computing log(S(λ)/S(1)) / log(λ) at a test sample.

pkg load symbolic;

printf('=========================================================\n');
printf('G-power (ℳ_α class) independent numerical audit\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function alpha_num = extract_alpha(func_handle, L0, G0_to_L0, theta0, e0, g_val, l_val)
  lambdas = [1.0, 1.5, 2.0, 0.7];
  vals = zeros(length(lambdas), 1);
  for i = 1:length(lambdas)
    lam = lambdas(i);
    L_v = lam * L0;
    G_v = G0_to_L0 * L_v;
    H_v = theta0 * G_v;
    l_v = l_val; g_v = g_val;
    vals(i) = func_handle(L_v, G_v, H_v, l_v, g_v);
  end
  ratios = vals(2:end) ./ vals(1);
  log_ratios = log(ratios) ./ log(lambdas(2:end)');
  alpha_num = -mean(log_ratios);
endfunction

function print_result(label, expected_alpha, actual_alpha, n_pass_in, n_fail_in)
  pass_flag = abs(expected_alpha - actual_alpha) < 0.01;
  printf('  %s\n', label);
  printf('    Expected α: %d\n', expected_alpha);
  printf('    Derived α:  %.6f\n', actual_alpha);
  if pass_flag
    printf('    Status: PASS\n');
    n_pass_in = n_pass_in + 1;
  else
    printf('    Status: FAIL\n');
    n_fail_in = n_fail_in + 1;
  end
  printf('\n');
  assignin('caller', 'n_pass', n_pass_in);
  assignin('caller', 'n_fail', n_fail_in);
endfunction

% Test parameters
L0 = 1.2; G_to_L = 0.9; theta0 = 0.4; e0 = sqrt(1 - G_to_L^2);
g_val = 0.7; l_val = 0.5;

% Kepler solver
function [E, f] = kep(l_v, e_v)
  E = l_v;
  for i = 1:100
    d = (E - e_v*sin(E) - l_v)/(1 - e_v*cos(E));
    E = E - d;
    if abs(d) < 1e-14; break; end
  end
  eta = sqrt(1 - e_v^2);
  f = 2*atan2(sqrt(1+e_v)*sin(E/2), sqrt(1-e_v)*cos(E/2));
endfunction

%% Test F_1 ∈ ℳ_6
function v = F1_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  [E, f] = kep(l, e);
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  v = (1+e*cos(f))^3 * (A + B_*cos(2*(f+g))) / G^6;
endfunction
a = extract_alpha(@F1_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('F_1 ∈ ℳ_6', 6, a, n_pass, n_fail);

%% Test F_1^* ∈ ℳ_6
function v = F1star_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2;
  v = A * eta^3 / G^6;
endfunction
a = extract_alpha(@F1star_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('F_1^* ∈ ℳ_6', 6, a, n_pass, n_fail);

%% Test dF_1^*/dL ∈ ℳ_7 (B.5.3)
function v = B53_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2;
  v = -3 * A * eta^4 / G^7;
endfunction
a = extract_alpha(@B53_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('B.5.3 dF_1^*/dL ∈ ℳ_7', 7, a, n_pass, n_fail);

%% Test dF_1^*/dG ∈ ℳ_7 (B.5.4)
function v = B54_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  v = -3 * eta^3 * (5*theta^2 - 1) / (2 * G^7);
endfunction
a = extract_alpha(@B54_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('B.5.4 dF_1^*/dG ∈ ℳ_7', 7, a, n_pass, n_fail);

%% Test dF_1^*/dH ∈ ℳ_7 (B.5.5)
function v = B55_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  v = 3 * theta * eta^3 / G^7;
endfunction
a = extract_alpha(@B55_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('B.5.5 dF_1^*/dH ∈ ℳ_7', 7, a, n_pass, n_fail);

%% Test dS_1/dl ∈ ℳ_3
function v = dS1_dl_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  [E, f] = kep(l, e);
  kap = 1 - e*cos(E);
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  v = ((A + B_*cos(2*(f+g)))/kap^3 - A/eta^3) / L^3;
endfunction
a = extract_alpha(@dS1_dl_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('E.6a dS_1/dl ∈ ℳ_3', 3, a, n_pass, n_fail);

%% Test S_1^{*,(T)} ∈ ℳ_3 (ch11a G.4)
function v = S1star_T_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  inner = A*B_*e^2*(6+e^2)/4*sin(2*g) + B_^2*e^4/128*sin(4*g);
  v = inner / ((5*theta^2 - 1) * eta^2 * G^3);
endfunction
a = extract_alpha(@S1star_T_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('G.4 S_1^{*,(T)} ∈ ℳ_3', 3, a, n_pass, n_fail);

%% Test c_2^{(T)} ∈ ℳ_{10}
function v = c2T_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  v = 3*A*B_*e^2*(6+e^2)*eta / (4 * G^10);
endfunction
a = extract_alpha(@c2T_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('c_2^{(T)} ∈ ℳ_{10}', 10, a, n_pass, n_fail);

%% Test c_0^{(T)} ∈ ℳ_{10} (F.15)
function v = c0T_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  Q = A^2 + B_^2/2;
  R_sym = A^2;
  Pe = 8 + 24*e^2 + 3*e^4;
  v = 3*eta / (2*G^10) * (Q*Pe/8 - R_sym*eta^3);
endfunction
a = extract_alpha(@c0T_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('F.15 c_0^{(T)} ∈ ℳ_{10}', 10, a, n_pass, n_fail);

%% Test dF_2^{**,(T)}/dH ∈ ℳ_{11}
function v = dF2ss_dH_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  Pe = 8 + 24*e^2 + 3*e^4;
  v = 9*theta / (2*G^11) * ((9*theta^2-5)*Pe*eta/16 - (3*theta^2-1)*eta^4);
endfunction
a = extract_alpha(@dF2ss_dH_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('G.10 dF_2^{**,(T)}/dH ∈ ℳ_{11}', 11, a, n_pass, n_fail);

%% Test dS_1^{*,(T)}/dg ∈ ℳ_3
function v = dS1s_dg_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  inner = A*B_*e^2*(6+e^2)/2*cos(2*g) + B_^2*e^4/32*cos(4*g);
  v = inner / ((5*theta^2 - 1) * eta^2 * G^3);
endfunction
a = extract_alpha(@dS1s_dg_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('G.9b dS_1^{*,(T)}/dg ∈ ℳ_3', 3, a, n_pass, n_fail);

%% Test dS_1^*/dL via G.6 form (T1.L at α=3), should be in ℳ_4
function v = dS1s_dL_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  % dPhi/de (treating eta = sqrt(1-e^2) as function of e)
  factor1 = 1/((5*theta^2-1) * eta^2);
  inner = A*B_*e^2*(6+e^2)/4*sin(2*g) + B_^2*e^4/128*sin(4*g);
  dfactor_de = 2*e/((5*theta^2-1) * eta^4);
  dinner_de = A*B_*(12*e + 4*e^3)/4*sin(2*g) + B_^2*4*e^3/128*sin(4*g);
  dPhi_de = dfactor_de * inner + factor1 * dinner_de;
  v = (eta^3/e) * dPhi_de / G^4;
endfunction
a = extract_alpha(@dS1s_dL_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('G.6 dS_1^*/dL ∈ ℳ_4', 4, a, n_pass, n_fail);

%% Test {F_1, S_1}_mod ∈ ℳ_{10} (via numerical_route_I structure)
% Use existing infrastructure via F_2^* = -(1/2) <{F_1, S_1}>_l
% Simpler: just test that g-averaged F_2^* is in ℳ_{10} via the c_0 + c_2 cos 2g + c_4 cos 4g
function v = F2_star_test_num(L, G, H, l, g)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B_ = 3*(1 - theta^2)/2;
  Q = A^2 + B_^2/2;
  R_sym = A^2;
  Pe = 8 + 24*e^2 + 3*e^4;
  c0 = 3*eta / (2*G^10) * (Q*Pe/8 - R_sym*eta^3);
  c2 = 3*A*B_*e^2*(6+e^2)*eta / (4 * G^10);
  c4 = 3*B_^2*e^4*eta / (64 * G^10);
  v = c0 + c2 * cos(2*g) + c4 * cos(4*g);
endfunction
a = extract_alpha(@F2_star_test_num, L0, G_to_L, theta0, e0, g_val, l_val);
print_result('F_2^{*,(T)} (via c_0+c_2 cos2g + c_4 cos4g) ∈ ℳ_{10}', 10, a, n_pass, n_fail);

printf('=========================================================\n');
if n_fail == 0
  printf('G-exponent audit: ALL %d checks PASSED.\n', n_pass);
  printf('All formulas scale correctly as G^{-α} under (L,G,H) → λ(L,G,H).\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
