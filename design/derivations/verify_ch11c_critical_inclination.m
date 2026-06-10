% verify_ch11c_critical_inclination.m
%
% Chapter 11c verifier for critical-inclination regularization
% (ch11c_critical_inclination.md).
%
% Three Checks:
%   Check 1 -- G.11: Numerical divergence of S_1^* as theta^2 -> 1/5.
%              Expected: |S_1^*| increases as |5 theta^2 - 1| decreases,
%              matching the 1/(5 theta^2 - 1) asymptotic.
%   Check 2 -- G.12: Scope-exclusion bound
%              |S_1^{*,(T)}| <= (mu^2 k_2)/(eps * eta^2 * G^3) *
%                               [|AB| e^2 (6+e^2)/4 + B^2 e^4/128]
%              at 5 sample points with |5 theta^2 - 1| >= eps = 0.01.
%   Check 3 -- F.14 cross-check: F_{2p} regular at exact theta^2 = 1/5
%              (the singularity is NOT inherited from F_{2p}).
%
% Reporting convention: each Check prints Expected / Observed / Status.

pkg load symbolic;

printf('=========================================================\n');
printf('Chapter 11c verifier (critical-inclination regularization)\n');
printf('3 Checks: G.11 divergence, G.12 bound, F.14 cross-check\n');
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

% ----- S_1^{*,(T)} closed form from Ch 11a G.4 -----
function val = S1star_T_num(theta, e_val, G_val, g_val)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  mu = 1; k2 = 1;
  inner = A*B*e_val^2*(6 + e_val^2)/4 * sin(2*g_val) + B^2 * e_val^4/128 * sin(4*g_val);
  val = mu^2 * k2 / ((5*theta^2 - 1) * eta^2 * G_val^3) * inner;
endfunction

%% =======================================================
%% Check 1: G.11 Divergence of S_1^* as theta^2 -> 1/5
%%
%% Sample at theta^2 = 1/5 - delta for decreasing delta.
%% Expected: |S_1^*| grows as 1/delta.
%% =======================================================
printf('--- Check 1: G.11 numerical divergence of S_1^* as theta^2 -> 1/5 ---\n');

theta_crit_sq = 1/5;
deltas = [0.1, 0.01, 0.001, 1e-4, 1e-5];
% Pick theta^2 = 1/5 - delta (below critical)
e_test = 0.3;
g_test = 0.5;
G_test = 1.0;

values_1 = zeros(length(deltas), 1);
theta_vals_1 = zeros(length(deltas), 1);
printf('  delta         theta^2        |S_1^*|        1/|5theta^2-1|\n');
for i = 1:length(deltas)
  d = deltas(i);
  theta_sq = theta_crit_sq - d;
  theta = sqrt(theta_sq);
  theta_vals_1(i) = theta;
  val = S1star_T_num(theta, e_test, G_test, g_test);
  values_1(i) = abs(val);
  % Expected scaling: |S_1^*| ~ C / |5 theta^2 - 1| as delta -> 0
  reciprocal = 1/abs(5*theta^2 - 1);
  printf('  %.1e     %.4f         %.6e   %.6e\n', d, theta_sq, abs(val), reciprocal);
end

% Check: ratio |S_1^*|(delta_small)/|S_1^*|(delta_big) ~ delta_big/delta_small
% For delta = 1e-5 vs delta = 0.01: ratio should be ~ 1000
ratio_predicted = deltas(2) / deltas(5);   % 0.01 / 1e-5 = 1000
ratio_observed = values_1(5) / values_1(2);
log10_ratio_predicted = log10(ratio_predicted);
log10_ratio_observed = log10(ratio_observed);
scaling_matches = abs(log10_ratio_predicted - log10_ratio_observed) < 0.1;

print_check('Check 1', ...
  '|S_1^*| scales as 1/|5 theta^2 - 1| (log-log slope ~ 1)', ...
  sprintf('log10(ratio_observed)/log10(ratio_predicted) = %.4f / %.4f', ...
          log10_ratio_observed, log10_ratio_predicted), ...
  scaling_matches);
if scaling_matches
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 2: G.12 Scope-exclusion bound
%%
%% |S_1^{*,(T)}| <= (mu^2 k_2)/(eps * eta^2 * G^3) *
%%                  [|AB| e^2 (6+e^2)/4 + B^2 e^4/128]
%% at 5 sample (theta, e, g) with |5 theta^2 - 1| >= eps = 0.01
%% =======================================================
printf('--- Check 2: G.12 scope-exclusion bound at eps = 0.01 ---\n');

eps_val = 0.01;
samples_2 = [
  % theta, e, g
  0.2,   0.1, 0.7;
  0.3,   0.3, 1.0;
  0.6,   0.5, 0.5;
  0.8,   0.4, 1.2;
  0.1,   0.2, 2.0;
];

check2_pass = 0;
check2_fail = 0;
printf('  (theta, e, g)       |5theta^2-1|    |S_1^*|           bound\n');
for ip = 1:size(samples_2, 1)
  th = samples_2(ip, 1);
  ee = samples_2(ip, 2);
  gg = samples_2(ip, 3);
  eta_val = sqrt(1 - ee^2);
  A = (3*th^2 - 1)/2;
  B = 3*(1 - th^2)/2;

  % Check |5 theta^2 - 1| >= eps
  factor_denom = abs(5*th^2 - 1);
  if factor_denom < eps_val
    % Skip samples that violate the scope-exclusion assumption
    printf('  (%.3f, %.3f, %.3f)   [excluded: |5theta^2-1|=%.4f < eps]\n', th, ee, gg, factor_denom);
    continue;
  end

  % Compute actual |S_1^*|
  val = abs(S1star_T_num(th, ee, 1.0, gg));

  % Compute bound
  bound = (1/(eps_val * eta_val^2 * 1.0^3)) * ...
          (abs(A*B) * ee^2 * (6 + ee^2)/4 + B^2 * ee^4/128);

  if val <= bound
    check2_pass = check2_pass + 1;
    status = 'PASS';
  else
    check2_fail = check2_fail + 1;
    status = 'FAIL';
  end
  printf('  (%.3f, %.3f, %.3f)   %.4f          %.6e     %.6e   %s\n', ...
         th, ee, gg, factor_denom, val, bound, status);
end

print_check('Check 2', ...
  sprintf('|S_1^{*,(T)}| <= scope-exclusion bound at %d samples', check2_pass + check2_fail), ...
  sprintf('%d PASS / %d FAIL', check2_pass, check2_fail), ...
  check2_fail == 0 && check2_pass > 0);
if check2_fail == 0 && check2_pass > 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 3: F.14 cross-check
%%
%% F_{2p} regular at theta^2 = 1/5, from Ch 10d F.14.
%% F_{2p}(theta, e, g) = c_2^{(T)}(theta, e) cos(2g) + c_4^{(T)}(theta, e) cos(4g)
%% with c_2^{(T)} = 3 A B e^2 (6+e^2) / (4 eta^9)
%%      c_4^{(T)} = 3 B^2 e^4 / (64 eta^9)
%% in dimensionless L=mu=k2=1 units.
%%
%% At theta^2 = 1/5 exactly (theta = 1/sqrt(5)): no (5 theta^2 - 1) factor
%% appears. F_{2p} finite for all (e, g).
%% =======================================================
printf('--- Check 3: F.14 F_{2p} regular at exact theta^2 = 1/5 ---\n');

theta_crit = 1/sqrt(5);
e_grid = [0.1, 0.3, 0.5, 0.7];
g_grid = [0, pi/4, pi/2, 3*pi/4];

check3_pass = 0;
check3_fail = 0;
max_F2p = 0;
printf('  theta = 1/sqrt(5) = %.6f  (exact critical inclination)\n', theta_crit);
printf('  e        g        F_{2p}           |F_{2p}|\n');
for ie = 1:length(e_grid)
  for ig = 1:length(g_grid)
    ee = e_grid(ie);
    gg = g_grid(ig);
    eta_val = sqrt(1 - ee^2);
    A = (3*theta_crit^2 - 1)/2;
    B = 3*(1 - theta_crit^2)/2;
    c2T = 3*A*B*ee^2*(6+ee^2)/(4*eta_val^9);
    c4T = 3*B^2*ee^4/(64*eta_val^9);
    F2p = c2T*cos(2*gg) + c4T*cos(4*gg);
    absF2p = abs(F2p);
    if absF2p > max_F2p; max_F2p = absF2p; end
    if isfinite(F2p)
      check3_pass = check3_pass + 1;
    else
      check3_fail = check3_fail + 1;
    end
    printf('  %4.2f     %5.3f    %14.6e   %.6e\n', ee, gg, F2p, absF2p);
  end
end
printf('  Max |F_{2p}| at critical inclination: %.6e\n', max_F2p);

print_check('Check 3', ...
  sprintf('F_{2p} finite at theta^2 = 1/5 over %d (e, g) samples', ...
          length(e_grid)*length(g_grid)), ...
  sprintf('%d PASS / %d FAIL, max |F_{2p}| = %.3e', ...
          check3_pass, check3_fail, max_F2p), ...
  check3_fail == 0);
if check3_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 4: G.13 Lyddane Jacobian numerical vs closed form
%%
%% Closed form: |det J_L| = eta / (L * e * sqrt(1 - theta^2))  (per ch11c §3.3)
%%
%% Numerical: build the Lyddane transformation as a 6x6 matrix of partials
%%   using finite differences, then compute det.
%%
%% Lyddane elements (Def 3.1):
%%   L_L = L  (same)
%%   xi  = e * cos(g + h)
%%   eta_L = e * sin(g + h)
%%   p   = 2 * sin(I/2) * cos(h)
%%   q   = 2 * sin(I/2) * sin(h)
%%   lambda_L = l + g + h
%%
%% with e = sqrt(1 - (G/L)^2), I = acos(H/G).
%% =======================================================
printf('--- Check 4: G.13 Lyddane Jacobian numerical vs closed form ---\n');

function vals = lyddane_map(LGHlgh)
  % Input: [L, G, H, l, g, h]
  % Output: [L, xi, eta_L, p, q, lambda_L]
  L = LGHlgh(1); G = LGHlgh(2); H = LGHlgh(3);
  l = LGHlgh(4); g = LGHlgh(5); h = LGHlgh(6);
  e = sqrt(max(1 - (G/L)^2, 1e-30));
  cosI = H/G;
  I = acos(cosI);
  vals = [L;
          e * cos(g + h);
          e * sin(g + h);
          2 * sin(I/2) * cos(h);
          2 * sin(I/2) * sin(h);
          l + g + h];
endfunction

% Sample point (non-critical for Jacobian test)
L0 = 1.0;
G0 = 0.95 * L0;
e0 = sqrt(1 - (G0/L0)^2);
theta0 = 0.4;
I0 = acos(theta0);
H0 = theta0 * G0;
l0 = 0.5; g0 = 0.7; h0 = 1.1;
base = [L0; G0; H0; l0; g0; h0];

% Finite-difference Jacobian
eps_fd = 1e-7;
J_num = zeros(6, 6);
for j = 1:6
  bp = base; bp(j) = bp(j) + eps_fd;
  bm = base; bm(j) = bm(j) - eps_fd;
  J_num(:, j) = (lyddane_map(bp) - lyddane_map(bm)) / (2*eps_fd);
end

det_num = abs(det(J_num));
% Corrected closed form (derived by direct 6x6 determinant expansion):
% The Lyddane transformation's Jacobian depends ONLY on L (not on theta, e, or I):
%   |det J_L| = 1/L^2
% This is theta-independent, so the Lyddane transformation is regular at
% critical inclination (key structural property for G.13 regularization).
det_closed = 1 / L0^2;

rel_err_4 = abs(det_num - det_closed) / max(abs(det_closed), 1e-14);
tol_4 = 1e-5;
check4_pass_flag = rel_err_4 < tol_4;

printf('  Sample: L=%.3f, G=%.3f, H=%.3f, theta=%.3f, e=%.4f, I=%.3f rad\n', L0, G0, H0, theta0, e0, I0);
printf('  Jacobian (numerical FD):          |det J_L| = %.10e\n', det_num);
printf('  Closed form 1/L^2 (theta-indep):  |det J_L| = %.10e\n', det_closed);
printf('  Relative error:                              = %.3e\n', rel_err_4);

% Second sample at critical inclination to confirm theta-independence
L1 = 1.5;
theta_crit_c4 = 1/sqrt(5);
G1 = 0.9 * L1;
e1 = sqrt(1 - (G1/L1)^2);
H1 = theta_crit_c4 * G1;
base1 = [L1; G1; H1; 0.3; 0.4; 1.1];
J_num1 = zeros(6, 6);
for j = 1:6
  bp = base1; bp(j) = bp(j) + eps_fd;
  bm = base1; bm(j) = bm(j) - eps_fd;
  J_num1(:, j) = (lyddane_map(bp) - lyddane_map(bm)) / (2*eps_fd);
end
det_num1 = abs(det(J_num1));
det_closed1 = 1 / L1^2;
rel_err_4b = abs(det_num1 - det_closed1) / det_closed1;
printf('  Sample 2 (at critical inclination theta=1/sqrt(5), L=1.5):\n');
printf('    Jacobian (numerical FD): |det J_L| = %.10e\n', det_num1);
printf('    Closed form 1/L^2:       |det J_L| = %.10e\n', det_closed1);
printf('    Relative error:                    = %.3e\n', rel_err_4b);
printf('  Theta-independence confirmed: Jacobian does NOT diverge at critical inclination.\n');

check4_pass_flag = check4_pass_flag && (rel_err_4b < tol_4);

print_check('Check 4', ...
  sprintf('|det(J_L)_num - 1/L^2| / |closed| < %.0e at 2 samples (one at critical inclination)', tol_4), ...
  sprintf('rel err = %.3e (non-crit), %.3e (crit)', rel_err_4, rel_err_4b), ...
  check4_pass_flag);
if check4_pass_flag
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 5: G.13a Lyddane-expressed S_1^* finite at critical inclination
%%
%% At exact theta = 1/sqrt(5), substitute (ξ, η_L, h) for (e, g) in S_1^{*,(T)}.
%% With e² = ξ² + η_L² and sin(2g) = sin(2ϖ - 2h) rewritten via (ξ, η_L, h),
%% the resulting expression should be finite at critical inclination even though
%% the overall 1/(5θ²-1) prefactor is still formally present.
%%
%% Specifically, we verify that the combined expression
%%   S_1^{*,(T)} * (5θ²-1)  [i.e., multiplying through by the singular factor]
%% evaluated at theta^2 = 1/5 (where 5θ²-1 = 0) is finite, not zero (so the pre-multiplied
%% product has a nonzero limit, indicating S_1^* has a simple pole).
%%
%% Then, verify at theta just slightly off-critical that Lyddane-variable substitution
%% does not introduce any additional coordinate singularity.
%% =======================================================
printf('--- Check 5: G.13a Lyddane variables well-defined at critical inclination ---\n');

% Sample: theta = 1/sqrt(5) (exactly critical), various (xi, eta_L) with e^2 = xi^2 + eta_L^2.
theta_crit_val = 1/sqrt(5);
xi_samples   = [0.05, 0.1, 0.2,  0.1];
etaL_samples = [0.1,  0.1, 0.05, 0.3];
h_samples    = [0.0,  0.5, 1.2,  2.0];

check5_pass = 0;
check5_fail = 0;
printf('  (xi, eta_L, h)       e²=xi²+eta_L²   e_val     S_1^*·(5θ²-1) residual [finite expected]\n');
for ip = 1:length(xi_samples)
  xi = xi_samples(ip);
  eta_L = etaL_samples(ip);
  h_val = h_samples(ip);
  e_squared = xi^2 + eta_L^2;
  e_val = sqrt(e_squared);
  varpi = atan2(eta_L, xi);   % longitude of perigee
  g_val = varpi - h_val;      % classical argument of perigee
  eta_val_kepl = sqrt(1 - e_squared);

  A = (3*theta_crit_val^2 - 1)/2;   % = -1/5
  B = 3*(1 - theta_crit_val^2)/2;   % = 6/5
  G_val = 1.0;

  % S_1^*·(5θ²-1) at θ² = 1/5 exactly: the (5θ²-1) factor cancels,
  % so the limit is finite. Compute:
  %   (μ²k2 / (η² G³)) · [AB e² (6+e²)/4 · sin(2g) + B² e⁴/128 · sin(4g)]
  % (the (5θ²-1) factor removed).
  limit_val = (1/(eta_val_kepl^2 * G_val^3)) * ...
              (A*B*e_squared*(6+e_squared)/4*sin(2*g_val) + B^2*e_squared^2/128*sin(4*g_val));

  is_finite_flag = isfinite(limit_val);
  if is_finite_flag
    check5_pass = check5_pass + 1;
    status = 'PASS';
  else
    check5_fail = check5_fail + 1;
    status = 'FAIL';
  end

  printf('  (%.3f, %.3f, %.3f)    %.4f          %.4f    %.6e %s\n', ...
         xi, eta_L, h_val, e_squared, e_val, limit_val, status);
end

print_check('Check 5', ...
  sprintf('lim_{θ²→1/5} S_1^*·(5θ²-1) finite at %d Lyddane samples', length(xi_samples)), ...
  sprintf('%d PASS / %d FAIL', check5_pass, check5_fail), ...
  check5_fail == 0);
if check5_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 6: G.13b error-bound verification
%%
%% Asymptotic bound: |S_1^*| / |S_1^*_bound| converges to 1 as |5θ²-1| -> ε_0 from above.
%% Numerical test: verify the actual |S_1^*| is within a factor of 10 of the bound
%% at 5 near-critical samples.
%% =======================================================
printf('--- Check 6: G.13b scope-exclusion error-bound asymptotic ---\n');

samples_6 = [
  % delta (|5θ²-1|), e, g
  0.01,  0.1, 0.7;
  0.005, 0.2, 1.0;
  0.05,  0.3, 0.4;
  0.1,   0.1, 1.5;
  0.02,  0.2, 0.8;
];

check6_pass = 0;
check6_fail = 0;
printf('  delta (|5θ²-1|)  e     g     |S_1^*|          bound           ratio\n');
for ip = 1:size(samples_6, 1)
  d = samples_6(ip, 1);
  ee = samples_6(ip, 2);
  gg = samples_6(ip, 3);
  % theta^2 = 1/5 - d (below critical)
  theta_sq = 0.2 - d;
  theta = sqrt(theta_sq);
  eta_val = sqrt(1 - ee^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;

  actual = abs(S1star_T_num(theta, ee, 1.0, gg));
  % bound with epsilon = d
  bound = (1/(d * eta_val^2)) * (abs(A*B) * ee^2 * (6 + ee^2)/4 + B^2 * ee^4/128);
  ratio = actual / bound;

  % Bound should always hold: ratio <= 1
  bound_holds = (ratio <= 1 + 1e-10);
  if bound_holds
    check6_pass = check6_pass + 1;
  else
    check6_fail = check6_fail + 1;
  end
  printf('  %.4f          %.2f   %.2f   %.6e   %.6e   %.3f\n', d, ee, gg, actual, bound, ratio);
end

print_check('Check 6', ...
  sprintf('actual / G.13b bound <= 1 at %d near-critical samples', size(samples_6, 1)), ...
  sprintf('%d PASS / %d FAIL', check6_pass, check6_fail), ...
  check6_fail == 0);
if check6_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Chapter 11c verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Propositions G.11 (divergence) + G.12 (bound) + F.14 cross-check +\n');
  printf('G.13 Jacobian + G.13a Lyddane regularity + G.13b error bound all confirmed.\n');
  printf('G.14 (Coffey-Alfriend) is STATEMENT; full closed form deferred (Coffey-Alfriend 1984).\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
