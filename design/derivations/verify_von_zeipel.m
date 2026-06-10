%% verify_von_zeipel.m
%% Numerical verification that the von Zeipel normalization works:
%% constructing S_1 from the homological equation produces F** that is
%% independent of the angle l'' to O(k_2^2).
%%
%% Uses a concrete 1-DOF Hamiltonian:
%%   F(L, l) = F_0(L) + k_2 * F_1(L, l)
%% where F_0(L) = mu^2/(2*L^2) and F_1(L, l) = (mu/L^2)^2 * cos(l)
%% (a simplified oblateness-like perturbation).

1;  % Octave: mark this file as a SCRIPT so the leading local function parses.
    % Without a statement before the first `function`, Octave treats the whole
    % file as a function file and errors executing the body ('k2' undefined).

%% Helper function (must appear before use in Octave scripts)
function [dLdt, dldt] = osc_eom(L, l, mu, k2)
  % Equations of motion for F(L, l) = mu^2/(2*L^2) + k2*(mu/L^2)^2*cos(l)
  % dL/dt = dF/dl = -k2*(mu/L^2)^2*sin(l)
  % dl/dt = -dF/dL = mu^2/L^3 + 2*k2*(mu^4/L^5)*cos(l)
  dLdt = -k2 * (mu/L^2)^2 * sin(l);
  dldt = mu^2/L^3 + 2*k2 * (mu^4/L^5) * cos(l);
end

printf('============================================================\n');
printf('VON ZEIPEL NORMALIZATION VERIFICATION\n');
printf('============================================================\n\n');

mu = 1.0;  % gravitational parameter (normalized)
n_tests_passed = 0;
n_tests_total = 0;

%% ================================================================
%% TEST 1: Construct S_1 from the homological equation and verify
%%          that F** is independent of l'' at first order.
%% ================================================================

printf('--- TEST 1: F** independence from l'''' ---\n\n');

L_pp = 1.5;  % mean action (held fixed)
k2 = 0.001;  % small oblateness parameter

% F_0(L) = mu^2/(2*L^2)
F0 = @(L) mu^2 ./ (2*L.^2);
dF0_dL = @(L) -mu^2 ./ L.^3;

% F_1(L, l) = (mu/L^2)^2 * cos(l)  [a simple periodic perturbation]
F_1 = @(L, l) (mu./L.^2).^2 .* cos(l);

% Mean of F_1 over l: <F_1> = (mu/L^2)^2 * <cos(l)> = 0
% So F_1_tilde = F_1 (the entire perturbation is oscillatory in this example)
F_1_mean = 0;
F_1_tilde = @(L, l) F_1(L, l) - F_1_mean;

% Homological equation: dS_1/dl = (L^3/mu^2) * F_1_tilde(L, l)
% = (L^3/mu^2) * (mu/L^2)^2 * cos(l)
% = (mu^2/L) * cos(l)  [after simplification]
%
% Integrate: S_1 = (mu^2/L) * sin(l)
S_1 = @(L, l) (mu^2./L) .* sin(l);
dS1_dl = @(L, l) (mu^2./L) .* cos(l);
dS1_dL = @(L, l) -(mu^2./L.^2) .* sin(l);

% Verify the homological equation is satisfied:
% F_1_tilde + (dF0/dL) * (dS1/dl) should = 0
n_tests_total = n_tests_total + 1;
l_test = linspace(0, 2*pi, 100);
residual = F_1_tilde(L_pp, l_test) + dF0_dL(L_pp) .* dS1_dl(L_pp, l_test);
max_residual = max(abs(residual));
printf('  Homological equation residual: %.2e (should be ~0)\n', max_residual);
if max_residual < 1e-14
  printf('  PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  FAIL\n\n');
end

% Now compute F** at many values of l'' and check it's constant
n_tests_total = n_tests_total + 1;
l_pp_vals = linspace(0, 2*pi, 200);
F_star_star = zeros(size(l_pp_vals));

for i = 1:length(l_pp_vals)
  lpp = l_pp_vals(i);

  % Inverse transformation: l = l'' - k2*dS1/dL'' + O(k2^2)
  l_osc = lpp - k2 * dS1_dL(L_pp, lpp);

  % L = L'' + k2*dS1/dl
  L_osc = L_pp + k2 * dS1_dl(L_pp, l_osc);

  % Evaluate the full Hamiltonian at osculating variables
  F_star_star(i) = F0(L_osc) + k2 * F_1(L_osc, l_osc);
end

F_variation = max(F_star_star) - min(F_star_star);
printf('  F** variation over l'''' in [0, 2pi]: %.6e\n', F_variation);
printf('  k_2^2 = %.6e\n', k2^2);
printf('  Ratio (variation / k2^2): %.4f\n', F_variation / k2^2);

if F_variation < 10 * k2^2
  printf('  PASS: F** variation is O(k_2^2) as expected\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  FAIL: F** variation is too large\n\n');
end

%% ================================================================
%% TEST 2: Verify the truncation error scales as k_2^2
%% ================================================================

printf('--- TEST 2: Truncation error scaling ---\n\n');

n_tests_total = n_tests_total + 1;
k2_values = [1e-2, 1e-3, 1e-4, 1e-5];
variations = zeros(size(k2_values));

for ik = 1:length(k2_values)
  k2_test = k2_values(ik);
  F_vals = zeros(size(l_pp_vals));

  for i = 1:length(l_pp_vals)
    lpp = l_pp_vals(i);
    l_osc = lpp - k2_test * dS1_dL(L_pp, lpp);
    L_osc = L_pp + k2_test * dS1_dl(L_pp, l_osc);
    F_vals(i) = F0(L_osc) + k2_test * F_1(L_osc, l_osc);
  end

  variations(ik) = max(F_vals) - min(F_vals);
end

printf('  k_2        | variation(F**)  | variation/k_2^2\n');
printf('  -----------|-----------------|----------------\n');
ratios = zeros(size(k2_values));
for ik = 1:length(k2_values)
  ratios(ik) = variations(ik) / k2_values(ik)^2;
  printf('  %.1e   | %.6e    | %.6f\n', k2_values(ik), variations(ik), ratios(ik));
end

% Check that the ratio is approximately constant (O(k2^2) scaling)
ratio_spread = max(ratios) / min(ratios);
printf('\n  Ratio spread (max/min): %.4f (should be ~1 for O(k2^2) scaling)\n', ratio_spread);

if ratio_spread < 2.0
  printf('  PASS: truncation error scales as k_2^2\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  FAIL: scaling is not quadratic\n\n');
end

%% ================================================================
%% TEST 3: Verify L'' is conserved under time evolution with F**
%% ================================================================

printf('--- TEST 3: L'''' conservation under F** dynamics ---\n\n');

n_tests_total = n_tests_total + 1;

% F** = F_0(L'') + k_2 * <F_1>(L'') + O(k_2^2)
% In our example, <F_1> = 0, so F** = F_0(L'') = mu^2/(2*L''^2)
% The equations of motion are:
%   dL''/dt = dF**/dl'' = 0   (exact, no l'' dependence)
%   dl''/dt = -dF**/dL'' = mu^2/L''^3 = n

% Integrate the FULL (not averaged) equations numerically
% using the osculating Hamiltonian F(L, l) = F0(L) + k2*F1(L, l)
% but in the mean-element variables via the generating function.

% Hamilton's equations in osculating variables:
%   dL/dt = dF/dl = k2 * dF1/dl = -k2 * (mu/L^2)^2 * sin(l)
%   dl/dt = -dF/dL = mu^2/L^3 - k2 * (-2*mu^2/L^5) * cos(l)
%         = mu^2/L^3 + 2*k2*(mu^2/L^5)*cos(l)

k2_dyn = 0.01;  % larger k2 to see dynamics
dt = 0.01;
n_steps = 10000;
T_total = n_steps * dt;

% Initial conditions
L_init = 1.5;
l_init = 0.5;

% Integrate osculating equations (RK4)
L_traj = zeros(1, n_steps+1);
l_traj = zeros(1, n_steps+1);
L_traj(1) = L_init;
l_traj(1) = l_init;

for step = 1:n_steps
  L_cur = L_traj(step);
  l_cur = l_traj(step);

  % RK4 for (L, l)
  [kL1, kl1] = osc_eom(L_cur, l_cur, mu, k2_dyn);
  [kL2, kl2] = osc_eom(L_cur + 0.5*dt*kL1, l_cur + 0.5*dt*kl1, mu, k2_dyn);
  [kL3, kl3] = osc_eom(L_cur + 0.5*dt*kL2, l_cur + 0.5*dt*kl2, mu, k2_dyn);
  [kL4, kl4] = osc_eom(L_cur + dt*kL3, l_cur + dt*kl3, mu, k2_dyn);

  L_traj(step+1) = L_cur + (dt/6)*(kL1 + 2*kL2 + 2*kL3 + kL4);
  l_traj(step+1) = l_cur + (dt/6)*(kl1 + 2*kl2 + 2*kl3 + kl4);
end

% Transform to mean elements using the inverse generating function
% L'' = L - k2*dS1/dl,  l'' = l + k2*dS1/dL''
% (iterative: L'' appears in dS1, so iterate)
L_pp_traj = zeros(size(L_traj));
l_pp_traj = zeros(size(l_traj));

for i = 1:length(L_traj)
  L_cur = L_traj(i);
  l_cur = l_traj(i);

  % Solve L = L'' + k2*dS1(L'', l)/dl for L''
  % Iterate: L''_new = L - k2*dS1(L''_old, l)/dl
  Lpp_est = L_cur;
  for iter = 1:20
    Lpp_est = L_cur - k2_dyn * dS1_dl(Lpp_est, l_cur);
  end
  L_pp_traj(i) = Lpp_est;

  % l'' = l + k2*dS1(L'', l)/dL''
  % dS1/dL = -(mu^2/L^2)*sin(l)
  l_pp_traj(i) = l_cur + k2_dyn * dS1_dL(Lpp_est, l_cur);
end

L_pp_variation = max(L_pp_traj) - min(L_pp_traj);
L_osc_variation = max(L_traj) - min(L_traj);

printf('  Integration: %d steps, dt = %.4f, T = %.1f\n', n_steps, dt, T_total);
printf('  k_2 = %.4f\n', k2_dyn);
printf('  Osculating L variation:   %.6e (oscillates with l)\n', L_osc_variation);
printf('  Mean L'''' variation:       %.6e (should be O(k_2^2))\n', L_pp_variation);
printf('  Ratio L_pp_var/k2^2:     %.4f\n', L_pp_variation / k2_dyn^2);

if L_pp_variation < 100 * k2_dyn^2
  printf('  PASS: L'''' is conserved to O(k_2^2)\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  FAIL: L'''' variation is too large\n\n');
end

%% ================================================================
%% TEST 4: Verify S_1 is periodic (zero-mean antiderivative condition)
%% ================================================================

printf('--- TEST 4: S_1 periodicity ---\n\n');

n_tests_total = n_tests_total + 1;

% S_1(L'', l) should be periodic in l with period 2*pi
S1_at_0 = S_1(L_pp, 0);
S1_at_2pi = S_1(L_pp, 2*pi);
periodicity_error = abs(S1_at_2pi - S1_at_0);

printf('  S_1(L'''', 0)    = %.10f\n', S1_at_0);
printf('  S_1(L'''', 2*pi) = %.10f\n', S1_at_2pi);
printf('  |S_1(2*pi) - S_1(0)| = %.2e\n', periodicity_error);

if periodicity_error < 1e-14
  printf('  PASS: S_1 is periodic\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  FAIL: S_1 is not periodic\n\n');
end

%% ================================================================
%% TEST 5: Verify that <F_1_tilde> = 0 (zero-mean condition)
%% ================================================================

printf('--- TEST 5: Zero-mean condition ---\n\n');

n_tests_total = n_tests_total + 1;

% Numerically integrate F_1_tilde over one period
N_quad = 10000;
l_quad = linspace(0, 2*pi, N_quad+1);
l_quad = l_quad(1:end-1);  % exclude endpoint for periodicity
F1_vals = F_1_tilde(L_pp, l_quad);
mean_F1 = mean(F1_vals);

printf('  <F_1_tilde> = %.6e (should be ~0)\n', mean_F1);

if abs(mean_F1) < 1e-14
  printf('  PASS: F_1_tilde has zero mean\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  FAIL: F_1_tilde has nonzero mean\n\n');
end

%% ================================================================
%% SUMMARY
%% ================================================================

printf('============================================================\n');
printf('SUMMARY: %d / %d tests PASSED\n', n_tests_passed, n_tests_total);
printf('============================================================\n');

if n_tests_passed == n_tests_total
  printf('ALL TESTS PASS\n');
else
  printf('SOME TESTS FAILED\n');
end

%% (Helper function osc_eom is defined at top of file.)
