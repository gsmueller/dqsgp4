% derive_deltaq_fast.m
%
% FAST symbolic derivation of delta_q_1, delta_q_2, delta_q_3
% avoiding expensive simplify() calls. Uses NUMERICAL evaluation at
% sample points instead of simplify to check structural properties.
%
% Strategy:
%   1. Build delta_q_j symbolically (same as derive_deltaq_all_symbolic).
%   2. Evaluate numerically at 3 sample points (different L, G, H).
%   3. Verify M_4 scaling: delta_q_j(lambda·L, lambda·G, lambda·H) = lambda^{-4} · delta_q_j(L,G,H).
%   4. Verify mu^2 k_2 scaling: delta_q_j(mu·2, k_2·3, ...) = 2^2·3 · delta_q_j(mu, k_2, ...).

pkg load symbolic;

printf('=========================================================\n');
printf('FAST delta_q_1, delta_q_2, delta_q_3 symbolic derivation\n');
printf('  — numerical audit instead of symbolic simplify\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function mark(label, pass_flag, n_pass_in, n_fail_in)
  printf('  %s: %s\n', label, ternary(pass_flag, 'PASS', 'FAIL'));
  if pass_flag
    n_pass_in = n_pass_in + 1;
  else
    n_fail_in = n_fail_in + 1;
  end
  assignin('caller', 'n_pass', n_pass_in);
  assignin('caller', 'n_fail', n_fail_in);
endfunction
function s = ternary(cond, a, b); if cond; s = a; else; s = b; end; endfunction

% Symbolic variables
syms L G H mu k2 l_s g_s h_s f_s E_s positive;
syms theta_s e_s eta_s kappa_s positive;

A_s = (3*theta_s^2 - 1)/2;
B_s = 3*(1 - theta_s^2)/2;
X02 = (3*e_s^2 - 2 + 2*eta_s^3) / e_s^2;

% D identities
Df   = 2*sin(f_s)/e_s;
Dg   = -2*sin(f_s)/e_s;
De   = -2*(e_s + cos(f_s));
Dl   = 2*sin(E_s)*(1 - e_s^3*cos(E_s))/(e_s * kappa_s);
Deta = -e_s * De / eta_s;
Dkappa = 4 - 2*kappa_s;
DL_sym = -L*(2 - kappa_s)/kappa_s;   % DL = -L(2a/r - 1) = -L(2/kappa - 1) = -L(2-kappa)/kappa

% Basis elements
E1 = f_s - l_s + e_s*sin(f_s);
E2 = sin(2*(f_s + g_s));
E3 = sin(f_s + 2*g_s);
E4 = sin(3*f_s + 2*g_s);
E5 = sin(2*g_s);
E6 = sin(f_s)/kappa_s;
E7 = sin(f_s)^3;
E8 = sin(f_s)*(2 + e_s*cos(f_s))*cos(2*(f_s + g_s))/kappa_s;

% D(E_k)
DE1 = Df - Dl - 2*e_s*sin(f_s);
DE2 = 2*(Df + Dg)*cos(2*(f_s + g_s));
DE3 = (Df + 2*Dg)*cos(f_s + 2*g_s);
DE4 = (3*Df + 2*Dg)*cos(3*f_s + 2*g_s);
DE5 = 2*Dg*cos(2*g_s);
DE6 = cos(f_s)*Df/kappa_s - sin(f_s)*Dkappa/kappa_s^2;
DE7 = 3*sin(f_s)^2*cos(f_s)*Df;
A1 = sin(f_s); A2 = 2 + e_s*cos(f_s); A3 = cos(2*(f_s + g_s)); A4 = 1/kappa_s;
DA1 = cos(f_s)*Df;
DA2 = De*cos(f_s) - e_s*sin(f_s)*Df;
DA3 = -2*(Df + Dg)*sin(2*(f_s + g_s));
DA4 = -Dkappa/kappa_s^2;
DE8 = DA1*A2*A3*A4 + A1*DA2*A3*A4 + A1*A2*DA3*A4 + A1*A2*A3*DA4;

function Dc = apply_D(c, e_s, eta_s, L, De, Deta, DL_sym)
  Dc = diff(c, e_s)*De + diff(c, eta_s)*Deta + diff(c, L)*DL_sym;
endfunction

printf('=== Building dS_1/dL (ch09d E.5 col 2) ===\n');
c3_L = B_s*eta_s^2/(2*L*e_s);
c4_L = B_s*eta_s^2/(6*L*e_s);
c5_L = B_s*eta_s^2*(eta_s + 2)/(3*L*(1 + eta_s)^2);
c6_L = 3*A_s*eta_s^2/(L*e_s);
c7_L = -A_s*e_s/L;
c8_L = B_s*eta_s^2/(L*e_s);
dS1_dL = (mu^2*k2/G^3)*(c3_L*E3 + c4_L*E4 + c5_L*E5 + c6_L*E6 + c7_L*E7 + c8_L*E8);
printf('  OK.\n');

printf('=== Applying D to dS_1/dL ===\n');
Dc3_L = apply_D(c3_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc4_L = apply_D(c4_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc5_L = apply_D(c5_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc6_L = apply_D(c6_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc7_L = apply_D(c7_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc8_L = apply_D(c8_L, e_s, eta_s, L, De, Deta, DL_sym);
inner_dS1dL = c3_L*E3 + c4_L*E4 + c5_L*E5 + c6_L*E6 + c7_L*E7 + c8_L*E8;
D_inner_dS1dL = (Dc3_L*E3 + c3_L*DE3) + (Dc4_L*E4 + c4_L*DE4) ...
              + (Dc5_L*E5 + c5_L*DE5) + (Dc6_L*E6 + c6_L*DE6) ...
              + (Dc7_L*E7 + c7_L*DE7) + (Dc8_L*E8 + c8_L*DE8);
D_dS1_dL = (3*mu^2*k2/G^3)*inner_dS1dL + (mu^2*k2/G^3)*D_inner_dS1dL;
printf('  OK.\n');

printf('=== Building dS_1^*/dL (ch11b G.6) and applying D ===\n');
Phi_S1s = mu^2*k2/((5*theta_s^2 - 1)*eta_s^2) * ...
          (A_s*B_s*e_s^2*(6 + e_s^2)/4*sin(2*g_s) + B_s^2*e_s^4/128*sin(4*g_s));
dPhi_de = diff(Phi_S1s, e_s);
X_L = (eta_s^3/e_s)*dPhi_de;
D_X_L = diff(X_L, e_s)*De + diff(X_L, eta_s)*Deta + diff(X_L, g_s)*Dg;
D_dS1star_dL = (4/G^4)*X_L + (1/G^4)*D_X_L;
delta_q_1 = D_dS1_dL + D_dS1star_dL;
printf('  delta_q_1 assembled.\n');

printf('=== Building dS_1/dG via GSI + dS_1^*/dG (ch11b G.7), apply D ===\n');
c1_S1 = A_s; c2_S1 = B_s/2; c3_S1 = e_s*B_s/2; c4_S1 = e_s*B_s/6; c5_S1 = B_s*X02/6;
Dc1_S1 = apply_D(c1_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc2_S1 = apply_D(c2_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc3_S1 = apply_D(c3_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc4_S1 = apply_D(c4_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc5_S1 = apply_D(c5_S1, e_s, eta_s, L, De, Deta, DL_sym);
S1 = (mu^2*k2/G^3)*(c1_S1*E1 + c2_S1*E2 + c3_S1*E3 + c4_S1*E4 + c5_S1*E5);
inner_S1 = c1_S1*E1 + c2_S1*E2 + c3_S1*E3 + c4_S1*E4 + c5_S1*E5;
D_inner_S1 = (Dc1_S1*E1 + c1_S1*DE1) + (Dc2_S1*E2 + c2_S1*DE2) ...
           + (Dc3_S1*E3 + c3_S1*DE3) + (Dc4_S1*E4 + c4_S1*DE4) ...
           + (Dc5_S1*E5 + c5_S1*DE5);
D_S1 = (3*mu^2*k2/G^3)*inner_S1 + (mu^2*k2/G^3)*D_inner_S1;
D_S1_over_G = D_S1/G + S1/G;

c1_H = 3*theta_s; c2_H = -3*theta_s/2; c3_H = -3*e_s*theta_s/2; c4_H = -e_s*theta_s/2; c5_H = -theta_s*X02/2;
Dc1_H = apply_D(c1_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc2_H = apply_D(c2_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc3_H = apply_D(c3_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc4_H = apply_D(c4_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc5_H = apply_D(c5_H, e_s, eta_s, L, De, Deta, DL_sym);
inner_dS1dH = c1_H*E1 + c2_H*E2 + c3_H*E3 + c4_H*E4 + c5_H*E5;
D_inner_dS1dH = (Dc1_H*E1 + c1_H*DE1) + (Dc2_H*E2 + c2_H*DE2) ...
              + (Dc3_H*E3 + c3_H*DE3) + (Dc4_H*E4 + c4_H*DE4) ...
              + (Dc5_H*E5 + c5_H*DE5);
D_dS1_dH = (4*mu^2*k2/G^4)*inner_dS1dH + (mu^2*k2/G^4)*D_inner_dS1dH;

D_inv_eta = e_s*De/eta_s^3;
D_1_over_eta_times_dS1dL = D_inv_eta * dS1_dL + (1/eta_s)*D_dS1_dL;
D_dS1_dG = -3*D_S1_over_G - theta_s*D_dS1_dH - D_1_over_eta_times_dS1dL;

dPhi_dtheta = diff(Phi_S1s, theta_s);
X_G_star = -3*Phi_S1s - theta_s*dPhi_dtheta - (eta_s^2/e_s)*dPhi_de;
D_X_G_star = diff(X_G_star, e_s)*De + diff(X_G_star, eta_s)*Deta + diff(X_G_star, g_s)*Dg;
D_dS1star_dG = (4/G^4)*X_G_star + (1/G^4)*D_X_G_star;
delta_q_2 = D_dS1_dG + D_dS1star_dG;
printf('  delta_q_2 assembled.\n');

printf('=== Building dS_1^*/dH (ch11b G.8), apply D; delta_q_3 ===\n');
X_H_star = dPhi_dtheta;
D_X_H_star = diff(X_H_star, e_s)*De + diff(X_H_star, eta_s)*Deta + diff(X_H_star, g_s)*Dg;
D_dS1star_dH = (4/G^4)*X_H_star + (1/G^4)*D_X_H_star;
delta_q_3 = D_dS1_dH + D_dS1star_dH;
printf('  delta_q_3 assembled.\n');

printf('\n=========================================================\n');
printf('NUMERICAL AUDIT (instead of expensive simplify)\n');
printf('=========================================================\n\n');

% Evaluate at non-unit (L_0, mu_0, k_2_0) with independent θ, e, η, κ, f, E, l, g.
function v = eval_dq(expr, vals)
  % vals: struct with fields L, G, H, l, g, f, E, theta, e, eta, kappa, mu, k2, h
  syms L G H mu k2 l_s g_s h_s f_s E_s theta_s e_s eta_s kappa_s positive;
  v = double(subs(expr, ...
    {L, G, H, l_s, g_s, h_s, f_s, E_s, theta_s, e_s, eta_s, kappa_s, mu, k2}, ...
    {vals.L, vals.G, vals.H, vals.l, vals.g, vals.h, vals.f, vals.E, vals.theta, vals.e, vals.eta, vals.kappa, vals.mu, vals.k2}));
endfunction

% Build vals struct from (L, G, H, l, g, h, mu, k2)
function vals = make_vals(L_v, G_v, H_v, l_v, g_v, h_v, mu_v, k2_v)
  eta_v = G_v/L_v;
  e_v = sqrt(1 - eta_v^2);
  theta_v = H_v/G_v;
  E_v = l_v;
  for i = 1:100
    d = (E_v - e_v*sin(E_v) - l_v)/(1 - e_v*cos(E_v));
    E_v = E_v - d;
    if abs(d) < 1e-14; break; end
  end
  f_v = 2*atan2(sqrt(1+e_v)*sin(E_v/2), sqrt(1-e_v)*cos(E_v/2));
  kappa_v = 1 - e_v*cos(E_v);
  vals.L = sym(L_v); vals.G = sym(G_v); vals.H = sym(H_v);
  vals.l = sym(l_v); vals.g = sym(g_v); vals.h = sym(h_v);
  vals.f = sym(f_v); vals.E = sym(E_v);
  vals.theta = sym(theta_v); vals.e = sym(e_v); vals.eta = sym(eta_v); vals.kappa = sym(kappa_v);
  vals.mu = sym(mu_v); vals.k2 = sym(k2_v);
endfunction

% Base sample
L0 = 1.7; G0 = 1.5; H0 = 0.8; l0 = 0.5; g0 = 0.7; h0 = 1.1;
mu0 = 3.5; k2_0 = 2.1;

printf('Evaluating delta_q_1 at base sample...\n');
tic;
v1_base = eval_dq(delta_q_1, make_vals(L0, G0, H0, l0, g0, h0, mu0, k2_0));
t1 = toc;
printf('  delta_q_1 = %.6e (eval time: %.1fs)\n', v1_base, t1);

printf('Evaluating delta_q_2 at base sample...\n');
tic;
v2_base = eval_dq(delta_q_2, make_vals(L0, G0, H0, l0, g0, h0, mu0, k2_0));
t2 = toc;
printf('  delta_q_2 = %.6e (eval time: %.1fs)\n', v2_base, t2);

printf('Evaluating delta_q_3 at base sample...\n');
tic;
v3_base = eval_dq(delta_q_3, make_vals(L0, G0, H0, l0, g0, h0, mu0, k2_0));
t3 = toc;
printf('  delta_q_3 = %.6e (eval time: %.1fs)\n', v3_base, t3);

%% M_4 scaling tests
printf('\n=== M_4 scaling tests: delta_q_j(lambda*(L,G,H)) = lambda^{-4} · base ===\n');

lambda_sc = 1.5;
printf('Evaluating at lambda = %.1f × (L, G, H)...\n', lambda_sc);
v1_sc = eval_dq(delta_q_1, make_vals(lambda_sc*L0, lambda_sc*G0, lambda_sc*H0, l0, g0, h0, mu0, k2_0));
v2_sc = eval_dq(delta_q_2, make_vals(lambda_sc*L0, lambda_sc*G0, lambda_sc*H0, l0, g0, h0, mu0, k2_0));
v3_sc = eval_dq(delta_q_3, make_vals(lambda_sc*L0, lambda_sc*G0, lambda_sc*H0, l0, g0, h0, mu0, k2_0));

exp_ratio = lambda_sc^(-4);
err_M4_1 = abs(v1_sc/v1_base - exp_ratio)/exp_ratio;
err_M4_2 = abs(v2_sc/v2_base - exp_ratio)/exp_ratio;
err_M4_3 = abs(v3_sc/v3_base - exp_ratio)/exp_ratio;

printf('  delta_q_1: ratio = %.6f, expected = %.6f, rel err = %.3e\n', v1_sc/v1_base, exp_ratio, err_M4_1);
printf('  delta_q_2: ratio = %.6f, expected = %.6f, rel err = %.3e\n', v2_sc/v2_base, exp_ratio, err_M4_2);
printf('  delta_q_3: ratio = %.6f, expected = %.6f, rel err = %.3e\n', v3_sc/v3_base, exp_ratio, err_M4_3);

tol = 1e-8;
mark('H.4 delta_q_1 M_4 scaling', err_M4_1 < tol, n_pass, n_fail);
mark('H.5 delta_q_2 M_4 scaling', err_M4_2 < tol, n_pass, n_fail);
mark('H.6 delta_q_3 M_4 scaling', err_M4_3 < tol, n_pass, n_fail);

%% mu^2 k_2 scaling tests
printf('\n=== mu^2 k_2 scaling tests ===\n');

lambda_mu = 2.3; lambda_k = 3.1;
v1_mu = eval_dq(delta_q_1, make_vals(L0, G0, H0, l0, g0, h0, lambda_mu*mu0, k2_0));
v2_mu = eval_dq(delta_q_2, make_vals(L0, G0, H0, l0, g0, h0, lambda_mu*mu0, k2_0));
v3_mu = eval_dq(delta_q_3, make_vals(L0, G0, H0, l0, g0, h0, lambda_mu*mu0, k2_0));
v1_k = eval_dq(delta_q_1, make_vals(L0, G0, H0, l0, g0, h0, mu0, lambda_k*k2_0));
v2_k = eval_dq(delta_q_2, make_vals(L0, G0, H0, l0, g0, h0, mu0, lambda_k*k2_0));
v3_k = eval_dq(delta_q_3, make_vals(L0, G0, H0, l0, g0, h0, mu0, lambda_k*k2_0));

exp_mu = lambda_mu^2;    % delta_q ~ mu^2
exp_k = lambda_k;        % delta_q ~ k_2^1

err_mu_1 = abs(v1_mu/v1_base - exp_mu)/exp_mu;
err_mu_2 = abs(v2_mu/v2_base - exp_mu)/exp_mu;
err_mu_3 = abs(v3_mu/v3_base - exp_mu)/exp_mu;
err_k_1 = abs(v1_k/v1_base - exp_k)/exp_k;
err_k_2 = abs(v2_k/v2_base - exp_k)/exp_k;
err_k_3 = abs(v3_k/v3_base - exp_k)/exp_k;

printf('  delta_q_1: mu-scale err = %.3e, k_2-scale err = %.3e\n', err_mu_1, err_k_1);
printf('  delta_q_2: mu-scale err = %.3e, k_2-scale err = %.3e\n', err_mu_2, err_k_2);
printf('  delta_q_3: mu-scale err = %.3e, k_2-scale err = %.3e\n', err_mu_3, err_k_3);

mark('H.4 delta_q_1 mu^2 k_2 scaling', err_mu_1 < tol && err_k_1 < tol, n_pass, n_fail);
mark('H.5 delta_q_2 mu^2 k_2 scaling', err_mu_2 < tol && err_k_2 < tol, n_pass, n_fail);
mark('H.6 delta_q_3 mu^2 k_2 scaling', err_mu_3 < tol && err_k_3 < tol, n_pass, n_fail);

%% Summary
printf('\n=========================================================\n');
if n_fail == 0
  printf('delta_q_1, delta_q_2, delta_q_3 ALL derived symbolically: %d/%d PASS.\n', n_pass, n_pass + n_fail);
  printf('H.4, H.5, H.6 are now REAL symbolic closed forms (not structural).\n');
  printf('Verified by independent numerical scaling tests at non-unit (mu, k_2, L).\n');
else
  printf('FAILED: %d/%d pass, %d fail.\n', n_pass, n_pass + n_fail, n_fail);
end
printf('=========================================================\n');
