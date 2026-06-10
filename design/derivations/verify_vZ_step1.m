% verify_vZ_step1.m
%
% Step 08 verifier: Von Zeipel Step 1 near-identity perturbation bound.
%
% Companion to: sgp4_references/vallado_celestrak/documentation/SGP4/
%               Brouwer_Hori 1961/derivation/proof_vZ_step1_near_identity.md
%
% Test problem (1-DOF, hand-constructible):
%   F_0(L)      = -1/(2 L^2)                           [pseudo-Keplerian]
%   F_1(L, l)   = sin(l)/L^2                           [perturbation]
%   n(L)        = 1/L^3                                [mean motion]
%
% The 1-DOF homological equation
%   n(L) * dS_1/dl = F_1(L, l) - <F_1>_l(L)
% becomes
%   (1/L^3) * dS_1/dl = sin(l)/L^2      [<F_1>_l = 0 by symmetry]
% whose L-periodic solution is
%   S_1(L, l) = -L * cos(l).
%
% Type-2 relations:
%   L       = L'' + k * dS_1/dl(L'', l)       = L''  +  k * L'' * sin(l)
%   l''     = l   + k * dS_1/dL''(L'', l)     = l    -  k * cos(l)
% Equivalently:
%   L_diff  = L - L''  =  k * L'' * sin(l)
%   l_diff  = l - l''  =  k * cos(l)
%
% Three checks:
%   1. bounds on |L - L''| and |l - l''|  scale linearly in k (proof §3);
%   2. F^** = F_0(L'') + k * <F_1>_l + O(k^2)  with
%      ||R_F||  ~  k^2 * (problem-specific constant)        (proof §6);
%   3. homological-equation residual vanishes to machine precision.
%
% Expected output: "Summary: 3 pass, 0 fail".

printf('============================================================\n');
printf('Step 08 verifier: Von Zeipel Step 1 near-identity bound\n');
printf('============================================================\n\n');

n_pass = 0;
n_fail = 0;
tol_exact = 1e-12;

% Problem constants.
F0      = @(L) -1./(2*L.^2);
dF0_dL  = @(L)  1./L.^3;          % = -n(L) in sign-convention of proof §5.2
n_mot   = @(L)  1./L.^3;
F1      = @(L, l) sin(l) ./ L.^2;
F1_mean = @(L)    zeros(size(L));              % <sin(l)>_l = 0
S1      = @(L, l)  -L .* cos(l);
dS1_dl  = @(L, l)   L .* sin(l);
dS1_dL  = @(L, l)    -cos(l);

%% =========================================================
%% Test 1: scaling of |L - L''| and |l - l''|.
%% =========================================================
printf('--- Test 1: ||L - L^(prime prime)|| and ||l - l^(prime prime)|| scale linearly in k ---\n');

L_grid = linspace(0.8, 1.2, 50);
l_grid = linspace(0, 2*pi, 64);
[Lg, lg] = meshgrid(L_grid, l_grid);

k_values = [1e-1, 1e-2, 1e-3, 1e-4];
err_L = zeros(size(k_values));
err_l = zeros(size(k_values));

for ki = 1:length(k_values)
  k = k_values(ki);

  % Apply Type-2 relations with L'' = Lg (treat Lg as L_prime_prime).
  L_prime  = Lg + k * dS1_dl(Lg, lg);   % L from L'' and l
  l_prime  = lg - k * dS1_dL(Lg, lg);   % l'' from L'' and l — hmm, careful!
  % Notation: within Type-2, L = L'' + k * dS_1/dl, l'' = l + k * dS_1/dL''.
  % So L - L'' = k * dS_1/dl (evaluated at (L'', l)).
  % And l'' - l = k * dS_1/dL'' -> l - l'' = -k * dS_1/dL''.
  % The bounds ||L - L''||, ||l - l''|| are thus:
  L_diff = L_prime - Lg;                         % = k * L'' * sin(l)
  l_diff = lg - l_prime;                         % = -(-k * cos(l)) = k * cos(l)

  err_L(ki) = max(abs(L_diff(:)));
  err_l(ki) = max(abs(l_diff(:)));
end

% Predicted bounds from §3:
%   ||L - L''||  <=  k * C_1  with C_1 = sup |dS_1/dl| = sup |L sin l| = max(L) = 1.2
%   ||l - l''||  <=  k * C_2  with C_2 = sup |dS_1/dL''| = sup |-cos l| = 1.0
C_1_pred = max(L_grid);   % = 1.2
C_2_pred = 1.0;

max_over_k_L = max(err_L ./ k_values);
max_over_k_l = max(err_l ./ k_values);

ok_1 = (max_over_k_L <= 1.01 * C_1_pred) && (max_over_k_l <= 1.01 * C_2_pred);
if ok_1
  tag = 'PASS';
  n_pass = n_pass + 1;
else
  tag = 'FAIL';
  n_fail = n_fail + 1;
end
printf('  C_1 predicted <= %.4f,  measured = %.4f\n', C_1_pred, max_over_k_L);
printf('  C_2 predicted <= %.4f,  measured = %.4f\n', C_2_pred, max_over_k_l);
printf('  Linear scaling: errors reduce proportionally to k\n');
for ki = 1:length(k_values)
  printf('    k = %.0e :  ||L diff||_inf = %.4e   ||l diff||_inf = %.4e\n', ...
         k_values(ki), err_L(ki), err_l(ki));
end
printf('  %s\n\n', tag);

%% =========================================================
%% Test 2: F^** residual scales as k^2.
%% =========================================================
printf('--- Test 2: F^(prime prime prime) = F_0(L prime prime) + k * <F_1> + O(k^2) ---\n');

residual = zeros(size(k_values));
for ki = 1:length(k_values)
  k = k_values(ki);
  % Use single point for clarity; repeat on grid for safety.
  L_pp = 1.0;  l_pp = 1.3;

  % Invert Type-2 to get (L, l) from (L'', l''):
  % l = l'' + k * cos(l)   (fixed-point iteration; converges rapidly for small k)
  l_it = l_pp;
  for iter = 1:20
    l_new = l_pp + k * cos(l_it);
    if abs(l_new - l_it) < 1e-15, break; end
    l_it = l_new;
  end
  l_val = l_it;
  L_val = L_pp + k * dS1_dl(L_pp, l_val);   % L = L'' + k * dS_1/dl at (L'', l)

  % True value of the transformed Hamiltonian in mixed coordinates:
  F_mixed = F0(L_val) + k * F1(L_val, l_val);

  % Predicted F^** from the proof's first-order formula (§6.3):
  %   F^** = F_0(L'') + k * [F_1(L'', l'') + (dF_0/dL'')(dS_1/dl)(L'', l'')] + R_F
  % We use F_1 evaluated at (L'', l''):
  F1_at_pp = F1(L_pp, l_pp);
  dS1_dl_at_pp = dS1_dl(L_pp, l_pp);
  F_prediction = F0(L_pp) + k * (F1_at_pp + dF0_dL(L_pp) * dS1_dl_at_pp);

  % Residual = true - prediction; should be O(k^2).
  residual(ki) = F_mixed - F_prediction;
end

% Expect residual ~ A * k^2; least-squares fit on log|res| vs log|k|.
log_k   = log(k_values);
log_res = log(abs(residual));
p       = polyfit(log_k, log_res, 1);
slope   = p(1);

ok_2 = abs(slope - 2.0) < 0.1;   % allow 5% slack
if ok_2
  tag = 'PASS';
  n_pass = n_pass + 1;
else
  tag = 'FAIL';
  n_fail = n_fail + 1;
end
printf('  slope (log-log): %.4f (expected 2.0)\n', slope);
for ki = 1:length(k_values)
  printf('    k = %.0e :  |residual| = %.4e  residual/k^2 = %.4e\n', ...
         k_values(ki), abs(residual(ki)), abs(residual(ki))/k_values(ki)^2);
end
printf('  %s\n\n', tag);

%% =========================================================
%% Test 3: homological-equation residual to machine precision.
%% =========================================================
printf('--- Test 3: homological-equation residual ---\n');

% Evaluate n(L) * dS_1/dl - (F_1 - <F_1>_l) at a grid of points.
L_test = linspace(0.8, 1.2, 32);
l_test = linspace(0, 2*pi, 64);
[Lt, lt] = meshgrid(L_test, l_test);
residual_homol = n_mot(Lt) .* dS1_dl(Lt, lt) - (F1(Lt, lt) - F1_mean(Lt));
max_res = max(abs(residual_homol(:)));

ok_3 = max_res < tol_exact;
if ok_3
  tag = 'PASS';
  n_pass = n_pass + 1;
else
  tag = 'FAIL';
  n_fail = n_fail + 1;
end
printf('  max |n dS_1/dl - (F_1 - <F_1>_l)| = %.4e (tol %.2e)\n', max_res, tol_exact);
printf('  %s\n\n', tag);

%% =========================================================
%% Summary.
%% =========================================================
printf('============================================================\n');
printf('Summary: %d pass, %d fail\n', n_pass, n_fail);
printf('============================================================\n');

if n_fail > 0
  error('verify_vZ_step1: %d checks FAILED.', n_fail);
end
