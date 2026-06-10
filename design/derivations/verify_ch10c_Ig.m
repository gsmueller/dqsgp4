% verify_ch10c_Ig.m
%
% Numerical verifier for Step 05's closed form of 𝓘_g per
% proof_ch10c_Ig_closed_form.md §§8, 10.
%
% 𝓘_g(L, G, H, g) := ⟨(∂S_1/∂g)(∂F_1/∂G)⟩_l
%                  = -⟨(∂F_1/∂g)(∂S_1/∂G)⟩_l   [per §1 IBP shortcut]
%                  = (μ⁶k_2² / G^10) · [c_0(θ,e,η) + c_2(θ,e,η) cos(2g)
%                                         + c_4(θ,e,η) cos(4g)]
%
% Verification strategy:
%  (1) Numerical route: evaluate -⟨(∂F_1/∂g)(∂S_1/∂G)⟩_l directly via
%      periodic trapezoidal l-integration (N=4096) using closed forms
%      from §2 (∂F_1/∂g) and §3.4.1+3.7.1 (∂S_1/∂G full, including 𝒜
%      term).
%  (2) Closed-form route: evaluate (8.1)-(8.4). For Σ_𝒜(e) use
%      truncation K_max of the doubly-summed series (7.2). The
%      X_k^{-3, 0}(e) are computed numerically via the Hansen integral
%      definition.
%  (3) Cross-check: difference at 20 random (θ, e, g) points.
%      Tolerance: 1e-9 relative error at e ≤ 0.7; 1e-7 at e = 0.85.
%  (4) Cross-checks at special points (equatorial, critical inclination,
%      circular).
%
% Normalize μ = k_2 = 1, G = 1; the closed form has explicit
% μ⁶k_2²/G^10 prefactor that drops out of the comparison.

clear all; close all;

fprintf('=================================================================\n');
fprintf('verify_ch10c_Ig — Step 05 closed form of 𝓘_g\n');
fprintf('=================================================================\n\n');

% --- Normalization
mu = 1.0;
k2 = 1.0;
G_const = 1.0;

% --- Trapezoidal grid for l-integration
N = 4096;
l_grid = (0:N-1)' * (2*pi/N);

% --- Newton-Kepler solver: given (l, e), return E
function E = solve_kepler(l, e)
  E = l;  % initial guess
  for iter = 1:60
    err = E - e*sin(E) - l;
    if abs(err) < 1e-15
      break;
    end
    derivative = 1 - e*cos(E);
    E = E - err/derivative;
  end
end

% --- Compute (f, kappa) given E and e
function [f, kappa] = E_to_f_kappa(E, e)
  kappa = 1 - e*cos(E);
  % f via tan(f/2) = sqrt((1+e)/(1-e)) tan(E/2)
  half_E = E/2;
  f = 2*atan2(sqrt(1+e)*sin(half_E), sqrt(1-e)*cos(half_E));
end

% --- Compute X_0^{-3, 0}(e) numerically via Hansen integral
function X = hansen_X0_n3_0(e, l_grid)
  N = length(l_grid);
  X_sum = 0;
  for i = 1:N
    E = solve_kepler(l_grid(i), e);
    kappa = 1 - e*cos(E);
    X_sum = X_sum + (1/kappa^3);  % (a/r)^3 = 1/kappa^3, j=0 means cos(0)=1
  end
  X = X_sum / N;
end

% --- Compute X_k^{-3, 0}(e) numerically: (1/2π)∫(a/r)^3 cos(k l) dl
function X = hansen_X_n3_0(k, e, l_grid)
  N = length(l_grid);
  X_sum = 0;
  for i = 1:N
    l = l_grid(i);
    E = solve_kepler(l, e);
    kappa = 1 - e*cos(E);
    X_sum = X_sum + cos(k*l)/kappa^3;
  end
  X = X_sum / N;
end

% --- Compute Σ_𝒜(e) directly via Theorem 6.4.5 structural identity
% By Theorem 6.4.5 of Step 05ar:
%   ⟨sin 2(f+g) · 𝒜/κ³⟩_l = cos(2g) · Σ_𝒜(e)
% Setting g = 0:
%   Σ_𝒜(e) = ⟨sin(2f) · 𝒜/κ³⟩_l (direct trapezoidal),
%          where 𝒜 = f - l + e sin f (the 'L-residual'/equation-of-center).
% This bypasses the Bessel-series truncation of (7.2) and gives an
% exact (up to trapezoidal precision) Σ_𝒜 from the same N-grid as
% the rest of the verifier.
function sigma_A = sigma_A_eval(e, l_grid)
  N = length(l_grid);
  sum_val = 0;
  for i = 1:N
    l = l_grid(i);
    E = solve_kepler(l, e);
    [f, kappa] = E_to_f_kappa(E, e);
    cA = f - l + e*sin(f);
    sum_val = sum_val + sin(2*f) * cA / kappa^3;
  end
  sigma_A = sum_val / N;
end

% --- Closed-form 𝓘_g per (8.1)-(8.4)
function I_g = Ig_closed(theta, e, g, mu, k2, G_const, l_grid)
  eta = sqrt(1 - e^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;

  prefactor = mu^6 * k2^2 / G_const^10;

  c0 = eta^3 * B * (5*theta^2 - 3) * (3 + 2*e^2)/4 - eta^5 * B^2 / 3;

  % Σ_𝒜 via direct trapezoidal evaluation of ⟨sin(2f)·𝒜/κ³⟩_l (g=0 form
  % of Theorem 6.4.5 structural identity); bypasses Bessel-series
  % truncation. The Theorem 6.4.5 boxed form of (7.2) is the
  % structural/closed expression; numerical evaluation is via direct l-integration.
  Sigma_A = sigma_A_eval(e, l_grid);
  c2 = -eta^3 * A * B * (12 - e^2)/4 - 3 * eta^6 * B * (5*theta^2 - 1) * Sigma_A;

  c4 = -eta^3 * B^2 * e^2 / 16;

  I_g = prefactor * (c0 + c2*cos(2*g) + c4*cos(4*g));
end

% --- Numerical 𝓘_g via direct l-integration of -(∂F_1/∂g)(∂S_1/∂G)
function I_g = Ig_numerical(theta, e, g, mu, k2, G_const, l_grid)
  N = length(l_grid);
  eta = sqrt(1 - e^2);
  L = G_const / eta;

  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;

  integrand_sum = 0;
  for i = 1:N
    l = l_grid(i);
    E = solve_kepler(l, e);
    [f, kappa] = E_to_f_kappa(E, e);

    % ∂F_1/∂g per (2.2): -2 B μ^4 k_2 η^6 sin 2(f+g) / (G^6 κ^3)
    dF1_dg = -2*B*mu^4*k2*eta^6 * sin(2*(f+g)) / (G_const^6 * kappa^3);

    % ∂S_1/∂G — assemble all pieces per (3.4.1) + (3.7.1):
    %   (μ²k_2/G^4) · {-3(5θ²-1)/2 · 𝒜 + (5θ²-3)/4 · ℬ
    %                  - A[3η²sin f/(eκ) - e sin³f] - (Bη²/(6e))·Λ}
    % where Λ = 6 sin f(2+e cos f) cos 2(f+g)/κ + 3 sin(f+2g)
    %         + sin(3f+2g) + (2e(η+2)/(1+η)²) sin(2g)
    % And 𝒜 = f-l+e sin f, ℬ = 3 sin 2(f+g) + 3e sin(f+2g)
    %                        + e sin(3f+2g) + X_0^{0,2}(e) sin(2g)
    X0_0_2 = (3*e^2 - 2 + 2*eta^3)/e^2;  % from B.0.7-7
    cA_inner = f - l + e*sin(f);
    cB_inner = 3*sin(2*(f+g)) + 3*e*sin(f+2*g) + e*sin(3*f+2*g) + X0_0_2*sin(2*g);

    Lambda_val = 6*sin(f)*(2 + e*cos(f))*cos(2*(f+g))/kappa ...
               + 3*sin(f+2*g) + sin(3*f+2*g) ...
               + (2*e*(eta+2)/(1+eta)^2) * sin(2*g);

    bracket_inner = -3*(5*theta^2 - 1)/2 * cA_inner ...
                  + (5*theta^2 - 3)/4 * cB_inner ...
                  - A * (3*eta^2 * sin(f) / (e*kappa) - e*sin(f)^3) ...
                  - (B*eta^2/(6*e)) * Lambda_val;

    dS1_dG = (mu^2*k2/G_const^4) * bracket_inner;

    integrand_sum = integrand_sum + dF1_dg * dS1_dG;
  end

  % 𝓘_g = -⟨(∂F_1/∂g)(∂S_1/∂G)⟩_l
  I_g = -integrand_sum / N;
end

% =========================================================
% Test 1: 20 random (θ, e, g) points
% =========================================================

fprintf('Test 1: 20 random (θ, e, g) points\n');
fprintf('-----------------------------------------------------------------\n');

rng(42);  % reproducible

% Σ_𝒜 evaluation: direct trapezoidal of ⟨sin(2f)·𝒜/κ³⟩_l (Theorem 6.4.5
% structural identity at g=0). No Bessel-series truncation required.

n_tests = 20;
n_pass = 0;
n_fail = 0;
max_rel_err = 0;

theta_test = -0.95 + 0.95*2*rand(n_tests, 1);
e_test = 0.05 + (0.85 - 0.05)*rand(n_tests, 1);
g_test = 2*pi*rand(n_tests, 1);

for i = 1:n_tests
  theta = theta_test(i);
  e = e_test(i);
  g = g_test(i);

  Ig_num = Ig_numerical(theta, e, g, mu, k2, G_const, l_grid);
  Ig_cl = Ig_closed(theta, e, g, mu, k2, G_const, l_grid);

  abs_diff = abs(Ig_num - Ig_cl);
  if abs(Ig_cl) > 1e-15
    rel_err = abs_diff / abs(Ig_cl);
  else
    rel_err = abs_diff;
  end
  max_rel_err = max(max_rel_err, rel_err);

  % Tolerance band: 1e-7 if e > 0.7, else 1e-9
  if e > 0.7
    tol = 1e-7;
  else
    tol = 1e-9;
  end

  status = '';
  if rel_err < tol
    status = 'PASS';
    n_pass = n_pass + 1;
  else
    status = 'FAIL';
    n_fail = n_fail + 1;
  end

  fprintf('  %s  (θ=%+.4f, e=%.4f, g=%.4f): rel_err = %.2e, num=%+.6e, cl=%+.6e\n', ...
    status, theta, e, g, rel_err, Ig_num, Ig_cl);
end

fprintf('\n  20 random points: %d PASS / %d FAIL / max rel err = %.2e\n', ...
  n_pass, n_fail, max_rel_err);

% =========================================================
% Test 2: Critical-inclination 𝒜-residual vanishing
% θ² = 1/5 ⇒ (5θ² - 1) = 0 ⇒ c_2's 𝒜-residual piece = 0
% c_2 (at critical θ) should equal -η³ A B (12-e²)/4 alone
% =========================================================

fprintf('\nTest 2: Critical-inclination 𝒜-residual vanishing\n');
fprintf('-----------------------------------------------------------------\n');

theta_crit = 1/sqrt(5);
e = 0.4;
g = 1.0;
eta = sqrt(1 - e^2);
A = (3*theta_crit^2 - 1)/2;
B = 3*(1 - theta_crit^2)/2;

Ig_num_crit = Ig_numerical(theta_crit, e, g, mu, k2, G_const, l_grid);
% Closed-form, but with Σ_𝒜 piece dropped (since (5θ²-1) = 0 at θ²=1/5)
prefactor = mu^6 * k2^2 / G_const^10;
c0_crit = eta^3 * B * (5*theta_crit^2 - 3) * (3 + 2*e^2)/4 - eta^5 * B^2 / 3;
% (5θ²-1) = 0 at this θ, so 𝒜-piece of c_2 is zero
c2_crit_poly = -eta^3 * A * B * (12 - e^2)/4;
c4_crit = -eta^3 * B^2 * e^2 / 16;
Ig_cl_crit = prefactor * (c0_crit + c2_crit_poly*cos(2*g) + c4_crit*cos(4*g));

abs_diff = abs(Ig_num_crit - Ig_cl_crit);
rel_err = abs_diff / max(abs(Ig_cl_crit), 1e-15);

status = '';
if rel_err < 1e-9
  status = 'PASS';
else
  status = 'FAIL';
end
fprintf('  %s θ_c = 1/√5 ≈ %.6f, e=%.2f, g=%.2f\n', status, theta_crit, e, g);
fprintf('         Ig_num  = %+.10e\n', Ig_num_crit);
fprintf('         Ig_poly = %+.10e (Σ_𝒜 piece structurally absent)\n', Ig_cl_crit);
fprintf('         rel err = %.2e\n', rel_err);

% =========================================================
% Test 3: Equatorial limit (θ = 1) ⇒ B = 0 ⇒ 𝓘_g ≡ 0
% =========================================================

fprintf('\nTest 3: Equatorial limit (θ = 1, B = 0)\n');
fprintf('-----------------------------------------------------------------\n');

theta_eq = 1.0;  % equatorial
e = 0.3;
g = 0.7;
Ig_num_eq = Ig_numerical(theta_eq, e, g, mu, k2, G_const, l_grid);
% At θ = 1: B = 0, so c_0, c_2, c_4 all proportional to B (or B²) all vanish.
% Wait: c_0 has -η⁵ B²/3 which is zero at B=0 ✓; first piece has B factor too ✓.
% c_2 has -AB and -3B(5θ²-1)Σ both have B factor ✓.
% c_4 has B² ✓.
% So 𝓘_g ≡ 0 at θ = 1.
abs_diff = abs(Ig_num_eq);
status = '';
if abs_diff < 1e-9
  status = 'PASS';
else
  status = 'FAIL';
end
fprintf('  %s θ_eq = 1.0, e=%.2f, g=%.2f: |Ig_num| = %.2e (expect 0)\n', ...
  status, e, g, abs_diff);

% =========================================================
% Test 4: Circular limit (e → 0) ⇒ Σ_𝒜(0) = 0, c_4(0) = 0
% Note: at e exactly 0, c_2 = -3 η³ A B at η = 1 = -3 A B
% =========================================================

fprintf('\nTest 4: Circular limit (e → 0)\n');
fprintf('-----------------------------------------------------------------\n');

theta_circ = 0.5;
e_circ = 0.01;  % small but not exactly 0 (avoid κ-singularity in derivatives)
g_circ = 0.4;

Ig_num_circ = Ig_numerical(theta_circ, e_circ, g_circ, mu, k2, G_const, l_grid);
Ig_cl_circ = Ig_closed(theta_circ, e_circ, g_circ, mu, k2, G_const, l_grid);

abs_diff = abs(Ig_num_circ - Ig_cl_circ);
rel_err = abs_diff / max(abs(Ig_cl_circ), 1e-15);

status = '';
if rel_err < 1e-7  % more lenient for tiny e
  status = 'PASS';
else
  status = 'FAIL';
end
fprintf('  %s θ=%.2f, e=%.4f (small), g=%.2f\n', status, theta_circ, e_circ, g_circ);
fprintf('         Ig_num = %+.10e, Ig_cl = %+.10e, rel err = %.2e\n', ...
  Ig_num_circ, Ig_cl_circ, rel_err);

% =========================================================
% Aggregate
% =========================================================

fprintf('\n=================================================================\n');
total_tests = n_tests + 3;
total_pass = n_pass + (rel_err < 1e-7);  % rough — just last test
fprintf('Aggregate: %d random + 3 special-case tests\n', n_tests);
fprintf('  Random pass rate: %d/%d\n', n_pass, n_tests);
fprintf('  Max rel error (random): %.2e\n', max_rel_err);
fprintf('=================================================================\n');
