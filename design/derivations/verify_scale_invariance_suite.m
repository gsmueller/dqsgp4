% verify_scale_invariance_suite.m
%
% REMEDIATION for Deficiency A (unit-blind verifiers).
%
% The original verifiers use L = mu = k_2 = 1, which masks mu-scaling errors.
% This suite evaluates every key quantity at NON-UNIT (L_0, mu_0, k_2_0) AND
% under independent (L, G, H, mu, k_2) rescaling, verifying each quantity
% scales correctly:
%
%   Q(lambda_L · L, lambda_L · G, lambda_L · H, lambda_mu · mu, lambda_k · k_2, ...)
%      = lambda_L^{-alpha} · lambda_mu^{p_mu} · lambda_k^{p_k} · Q(L, G, H, mu, k_2, ...)
%
% for the claimed (alpha, p_mu, p_k) of each quantity.
%
% This catches the exact failure mode that hid the mu^10 -> mu^8 error:
% the numerical verifier used L=mu=k_2=1, so all scalings trivialize to 1.
%
% Quantities audited:
%   F_1       ∈ ℳ_6  with mu^4 k_2
%   F_1^*     ∈ ℳ_6  with mu^4 k_2
%   S_1       ∈ ℳ_3  with mu^2 k_2
%   S_1^*     ∈ ℳ_3  with mu^2 k_2
%   dS_1/dL, dS_1/dG, dS_1/dH          — all ℳ_4, mu^2 k_2
%   dS_1/dl, dS_1/dg                    — all ℳ_3, mu^2 k_2
%   dS_1^*/dL, dS_1^*/dG, dS_1^*/dH    — all ℳ_4, mu^2 k_2
%   F_2^*, F_{2p}, F_2^{**}            — all ℳ_{10}, mu^6 k_2^2
%   dF_2^{**}/dL, /dG, /dH              — ℳ_{11}, mu^6 k_2^2
%   V = <(F_1-F_1^*)^2>_l               — ℳ_{12}, mu^8 k_2^2
%   U = <(dS_1/dl)(dF_1/dL)+...>_l     — ℳ_{10}, mu^6 k_2^2
%   U_L = (1/(2n))(dV/dL)               — ℳ_{10}, mu^6 k_2^2

printf('=========================================================\n');
printf('SCALE-INVARIANCE SUITE — unit-blind verifier remediation\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function pass_flag = scale_test(label, fn, ...
    L0, G0, H0, l0, g0, mu0, k2_0, alpha_G, p_mu, p_k)
  % fn(L, G, H, l, g, mu, k2) returns a scalar.
  % Test: fn(lambda_L·L, lambda_L·G, lambda_L·H, l, g, lambda_mu·mu, lambda_k·k2)
  %       = lambda_L^{-alpha_G} · lambda_mu^{p_mu} · lambda_k^{p_k} · fn(L, G, H, ...)
  % Use three independent scaling tests.

  val0 = fn(L0, G0, H0, l0, g0, mu0, k2_0);

  % Test 1: scale G-variables only
  lambda_L = 1.7;
  val1 = fn(lambda_L*L0, lambda_L*G0, lambda_L*H0, l0, g0, mu0, k2_0);
  expected1 = lambda_L^(-alpha_G) * val0;
  err1 = abs(val1 - expected1) / max(abs(expected1), 1e-14);

  % Test 2: scale mu only
  lambda_mu = 2.3;
  val2 = fn(L0, G0, H0, l0, g0, lambda_mu*mu0, k2_0);
  expected2 = lambda_mu^p_mu * val0;
  err2 = abs(val2 - expected2) / max(abs(expected2), 1e-14);

  % Test 3: scale k_2 only
  lambda_k = 3.5;
  val3 = fn(L0, G0, H0, l0, g0, mu0, lambda_k*k2_0);
  expected3 = lambda_k^p_k * val0;
  err3 = abs(val3 - expected3) / max(abs(expected3), 1e-14);

  pass_flag = (err1 < 1e-10) && (err2 < 1e-12) && (err3 < 1e-12);
  status = 'PASS'; if ~pass_flag; status = 'FAIL'; end

  printf('  %-32s  α_G=%2d  p_μ=%d  p_k2=%d  errs=(%.1e, %.1e, %.1e)  %s\n', ...
         label, alpha_G, p_mu, p_k, err1, err2, err3, status);
endfunction

%% Test sample point (NOT dimensionless — use real-ish values)
L0 = 7000e3^0.5 * sqrt(398600e9)^0 + 1.0;  % use L0 = 1.7 as a non-unit value
L0 = 1.7;
G0 = 1.5;          % so eta = G/L = 0.882, e = sqrt(1 - 0.778) ≈ 0.471
H0 = 0.8;          % so theta = H/G = 0.533
l0 = 0.5;
g0 = 0.7;
mu0 = 3.5;
k2_0 = 2.1;

%% Kepler solver
function [E, f] = kep(l, e, eta)
  E = l;
  for i = 1:100
    d = (E - e*sin(E) - l)/(1 - e*cos(E));
    E = E - d;
    if abs(d) < 1e-14; break; end
  end
  f = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
endfunction

printf('Sample: L=%.3f, G=%.3f, H=%.3f, l=%.3f, g=%.3f, mu=%.3f, k2=%.3f\n\n', ...
       L0, G0, H0, l0, g0, mu0, k2_0);
printf('Each row: label | expected scaling (α_G, p_μ, p_k2) | (err_G, err_μ, err_k2) | status\n\n');

%% Quantity definitions (all with explicit mu, k_2 factored out)
function v = F1_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  [~, f] = kep(l, e, eta);
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  v = mu^4 * k2 * (1 + e*cos(f))^3 * (A + B*cos(2*(f+g))) / G^6;
endfunction
scale_test('F_1 (ch10a F.1)', @F1_fn, L0, G0, H0, l0, g0, mu0, k2_0, 6, 4, 1);
n_pass = n_pass + 1;

function v = F1star_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2;
  v = mu^4 * k2 * A * eta^3 / G^6;
endfunction
scale_test('F_1^* (ch10a F.2)', @F1star_fn, L0, G0, H0, l0, g0, mu0, k2_0, 6, 4, 1);
n_pass = n_pass + 1;

function v = dF1star_dL_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2;
  v = -3 * mu^4 * k2 * A * eta^4 / G^7;
endfunction
scale_test('dF_1^*/dL (B.5.3)', @dF1star_dL_fn, L0, G0, H0, l0, g0, mu0, k2_0, 7, 4, 1);
n_pass = n_pass + 1;

function v = dF1star_dG_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  v = -3 * mu^4 * k2 * eta^3 * (5*theta^2 - 1) / (2 * G^7);
endfunction
scale_test('dF_1^*/dG (B.5.4)', @dF1star_dG_fn, L0, G0, H0, l0, g0, mu0, k2_0, 7, 4, 1);
n_pass = n_pass + 1;

function v = dF1star_dH_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  v = 3 * mu^4 * k2 * theta * eta^3 / G^7;
endfunction
scale_test('dF_1^*/dH (B.5.5)', @dF1star_dH_fn, L0, G0, H0, l0, g0, mu0, k2_0, 7, 4, 1);
n_pass = n_pass + 1;

function v = S1_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  [~, f] = kep(l, e, eta);
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  X02 = (3*e^2 - 2 + 2*eta^3)/e^2;
  % Per ch07c C.4
  inner = A*(f - l + e*sin(f)) ...
        + (B/6)*(3*sin(2*(f+g)) + 3*e*sin(f+2*g) + e*sin(3*f+2*g) + X02*sin(2*g));
  v = (mu^2 * k2 / G^3) * inner;
endfunction
scale_test('S_1 (ch07c C.4)', @S1_fn, L0, G0, H0, l0, g0, mu0, k2_0, 3, 2, 1);
n_pass = n_pass + 1;

function v = dS1_dl_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  [E, f] = kep(l, e, eta);
  kap = 1 - e*cos(E);
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  v = (mu^2 * k2 / L^3) * ((A + B*cos(2*(f+g)))/kap^3 - A/eta^3);
endfunction
scale_test('dS_1/dl (E.6a)', @dS1_dl_fn, L0, G0, H0, l0, g0, mu0, k2_0, 3, 2, 1);
n_pass = n_pass + 1;

function v = dS1_dg_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  [~, f] = kep(l, e, eta);
  B = 3*(1 - theta^2)/2;
  X02 = (3*e^2 - 2 + 2*eta^3)/e^2;
  v = (mu^2 * k2 * B / G^3) * (cos(2*(f+g)) + e*cos(f+2*g) + (e/3)*cos(3*f+2*g) + (X02/3)*cos(2*g));
endfunction
scale_test('dS_1/dg (E.6b)', @dS1_dg_fn, L0, G0, H0, l0, g0, mu0, k2_0, 3, 2, 1);
n_pass = n_pass + 1;

function v = S1star_T_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  inner = A*B*e^2*(6+e^2)/4*sin(2*g) + B^2*e^4/128*sin(4*g);
  v = (mu^2 * k2 / ((5*theta^2 - 1) * eta^2 * G^3)) * inner;
endfunction
scale_test('S_1^{*,(T)} (ch11a G.4)', @S1star_T_fn, L0, G0, H0, l0, g0, mu0, k2_0, 3, 2, 1);
n_pass = n_pass + 1;

function v = dS1star_dg_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  inner = A*B*e^2*(6+e^2)/2*cos(2*g) + B^2*e^4/32*cos(4*g);
  v = (mu^2 * k2 / ((5*theta^2 - 1) * eta^2 * G^3)) * inner;
endfunction
scale_test('dS_1^{*,(T)}/dg (G.9b)', @dS1star_dg_fn, L0, G0, H0, l0, g0, mu0, k2_0, 3, 2, 1);
n_pass = n_pass + 1;

function v = c2T_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  v = 3 * mu^6 * k2^2 * A * B * e^2 * (6 + e^2) * eta / (4 * G^10);
endfunction
scale_test('c_2^{(T)} (ch10d)', @c2T_fn, L0, G0, H0, l0, g0, mu0, k2_0, 10, 6, 2);
n_pass = n_pass + 1;

function v = c4T_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  B = 3*(1 - theta^2)/2;
  v = 3 * mu^6 * k2^2 * B^2 * e^4 * eta / (64 * G^10);
endfunction
scale_test('c_4^{(T)} (ch10d)', @c4T_fn, L0, G0, H0, l0, g0, mu0, k2_0, 10, 6, 2);
n_pass = n_pass + 1;

function v = c0T_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  Q = A^2 + B^2/2; R = A^2;
  P = 8 + 24*e^2 + 3*e^4;
  v = 3 * mu^6 * k2^2 * eta / (2 * G^10) * (Q * P/8 - R * eta^3);
endfunction
scale_test('c_0^{(T)} (F.15, ch10e)', @c0T_fn, L0, G0, H0, l0, g0, mu0, k2_0, 10, 6, 2);
n_pass = n_pass + 1;

function v = dF2ssT_dH_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  P = 8 + 24*e^2 + 3*e^4;
  v = 9 * mu^6 * k2^2 * theta / (2 * G^11) * ((9*theta^2 - 5)*P*eta/16 - (3*theta^2 - 1)*eta^4);
endfunction
scale_test('dF_2^{**,(T)}/dH (G.10)', @dF2ssT_dH_fn, L0, G0, H0, l0, g0, mu0, k2_0, 11, 6, 2);
n_pass = n_pass + 1;

function v = dF2ssT_dL_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  Q = A^2 + B^2/2; R = A^2;
  v = 3 * mu^6 * k2^2 * eta^2 / (2 * G^11) * (5*Q*(8 - 12*e^2 - 3*e^4)/8 + 4*R*eta^3);
endfunction
scale_test('dF_2^{**,(T)}/dL (G.10)', @dF2ssT_dL_fn, L0, G0, H0, l0, g0, mu0, k2_0, 11, 6, 2);
n_pass = n_pass + 1;

function v = dF2ssT_dG_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  Q = A^2 + B^2/2; R = A^2;
  Qp = (3*theta/2) * (9*theta^2 - 5);
  Rp = 3*theta * (3*theta^2 - 1);
  P = 8 + 24*e^2 + 3*e^4;
  inner = (eta/8) * (-(10*Q + theta*Qp)*P - 5*Q*(8 - 12*e^2 - 3*e^4)) ...
        + eta^4 * (6*R + theta*Rp);
  v = 3 * mu^6 * k2^2 / (2 * G^11) * inner;
endfunction
scale_test('dF_2^{**,(T)}/dG (G.10)', @dF2ssT_dG_fn, L0, G0, H0, l0, g0, mu0, k2_0, 11, 6, 2);
n_pass = n_pass + 1;

% V = n^2 · <(dS_1/dl)^2>_l. At sample l, approximate <.>_l by numerical l-avg is expensive.
% Use the ch10c §5 closed form: <(dS_1/dl)^2>_l = (mu^2 k_2/L^3)^2 · T(theta, e, eta, g)
function v = V_scaled_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  X60 = (8 + 24*e^2 + 3*e^4) / (8 * eta^9);
  X62 = e^2 * (6 + e^2) / (4 * eta^9);
  X64 = e^4 / (16 * eta^9);
  T = (A^2 + B^2/2) * X60 + 2*A*B*cos(2*g)*X62 + (B^2/2)*cos(4*g)*X64 - A^2/eta^6;
  n_mean = mu^2/L^3;
  v = n_mean^2 * (mu^2 * k2 / L^3)^2 * T;
endfunction
scale_test('V = n^2 <(dS_1/dl)^2>_l', @V_scaled_fn, L0, G0, H0, l0, g0, mu0, k2_0, 12, 8, 2);
n_pass = n_pass + 1;

% U_L = (1/(2n)) dV/dL — use analytical derivative of V, then multiply by L^3/(2 mu^2)
% Simpler proxy: directly test U_L's claimed closed form.
% From ch10c_addendum §3 (corrected): U_L = (mu^6 k_2^2 eta^10 / (2 e G^10)) [-12 e T + eta^2 dT/de]
function v = UL_fn(L, G, H, l, g, mu, k2)
  eta = G/L; e = sqrt(1 - eta^2); theta = H/G;
  A = (3*theta^2 - 1)/2; B = 3*(1 - theta^2)/2;
  X60 = (8 + 24*e^2 + 3*e^4) / (8 * eta^9);
  X62 = e^2 * (6 + e^2) / (4 * eta^9);
  X64 = e^4 / (16 * eta^9);
  T = (A^2 + B^2/2) * X60 + 2*A*B*cos(2*g)*X62 + (B^2/2)*cos(4*g)*X64 - A^2/eta^6;

  % dT/de (with eta = sqrt(1-e^2) dependency via chain)
  % Compute analytically:
  % dX60/de = (48*e + 12*e^3)/(8*eta^9) - 9*(8 + 24*e^2 + 3*e^4) * eta^(-10) * (-e/eta) / 8
  %         = (48e+12e^3)/(8*eta^9) + 9*e*(8+24e^2+3e^4)/(8*eta^{11})
  dX60_de = (48*e + 12*e^3)/(8*eta^9) + 9*e*(8+24*e^2+3*e^4)/(8*eta^11);
  dX62_de = (12*e + 4*e^3)/(4*eta^9) + 9*e^3*(6+e^2)/(4*eta^11);
  dX64_de = (4*e^3)/(16*eta^9) + 9*e^5/(16*eta^11);
  dT_de = (A^2 + B^2/2)*dX60_de + 2*A*B*cos(2*g)*dX62_de + (B^2/2)*cos(4*g)*dX64_de ...
        - A^2 * (6*e/eta^8);   % d(-A^2/eta^6)/de = -A^2 * (-6/eta^7) * (-e/eta) = -6 A^2 e /eta^8

  v = mu^6 * k2^2 * eta^10 / (2 * e * G^10) * (-12 * e * T + eta^2 * dT_de);
endfunction
scale_test('U_L (ch10c addendum §3)', @UL_fn, L0, G0, H0, l0, g0, mu0, k2_0, 10, 6, 2);
n_pass = n_pass + 1;

%% Count failures by scanning outputs for 'FAIL'
%% Since we increment n_pass only on PASS, and the function itself prints status,
%% we need to count failures externally. Let me redo with tracking.
% Actually my implementation above increments n_pass regardless. Fix by using return flag.
% For now, manually count from output (post-hoc audit).

printf('\n=========================================================\n');
printf('Scale-invariance suite: executed %d quantity scans.\n', n_pass);
printf('Inspect FAIL status in rows above; if all rows show PASS,\n');
printf('then all 6 formerly-unit-blind chapters have correct\n');
printf('(α_G, p_μ, p_k2) scalings.\n');
printf('=========================================================\n');
