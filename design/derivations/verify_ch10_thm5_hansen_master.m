%% verify_ch10_thm5_hansen_master.m
%% Numerical verification of Chapter 10 Theorem 5 (Hansen master formula).
%%
%% Theorem 5 (revised form — see ch10_foundations_thm5.md):
%%   For S in M_alpha with harmonic decomposition
%%     S = G^{-alpha} * sum c_{jkm}(theta, e, eta) * cos(j*f + k*g) * kappa^{-m}
%%                   + G^{-alpha} * sum s_{jkm}(theta, e, eta) * sin(j*f + k*g) * kappa^{-m},
%%   the l-average satisfies, harmonic-by-harmonic:
%%     <cos(j*f + k*g) * kappa^{-m}>_l = cos(k*g) * X_0^{-m, j}(e),
%%     <sin(j*f + k*g) * kappa^{-m}>_l = sin(k*g) * X_0^{-m, j}(e),
%%   where X_0^{n, m}(e) = (1/2pi) \int_0^{2pi} (r/a)^n * exp(i*m*f) dl
%%   (ch06a Definition B.0.1), and X_0^{n, m} is real (ch06a Theorem B.0.2(i)).
%%
%% Special case: pure k = 0 sine harmonics vanish because sin(0*g) = 0.
%%
%% Strategy: numerically evaluate each l-average by converting
%%     dl = kappa * dE,       with kappa = 1 - e*cos(E),
%% and using the geometric identities
%%     r/a = kappa,
%%     cos(f) = (cos(E) - e) / kappa,
%%     sin(f) = eta * sin(E) / kappa,   with eta = sqrt(1 - e^2).
%%
%% Each test compares the numerical l-average against the Hansen-coefficient
%% closed form from ch06 Corollary B.0.7.

printf('=================================================================\n');
printf('CH10 THEOREM 5 (HANSEN MASTER FORMULA) — VERIFICATION\n');
printf('=================================================================\n\n');

%% Test parameters
e_val = 0.3;                  % eccentricity
g_val = 1.0;                  % argument of perigee (radians)
eta_val = sqrt(1 - e_val^2);

printf('Test parameters: e = %.4f, g = %.4f rad, eta = %.6f\n\n', ...
       e_val, g_val, eta_val);

%% Hansen-coefficient closed forms (from ch06 Corollary B.0.7):
%%   X_0^{0, 0}    = 1                          (trivial)
%%   X_0^{-3, 0}   = 1 / eta^3                  (B.0.7-1)
%%   X_0^{-3, 2}   = 0                          (B.0.7-2)
%%   X_0^{0, 1}    = -e                         (B.0.7-6)
%%   X_0^{0, 2}    = (3*e^2 - 2 + 2*eta^3)/e^2   (B.0.7-7)
%%   X_0^{0, 3}    = (8 - 12*e^2 + 3*e^4 - 8*eta^3)/e^3   (ch07c Lemma C.4.9)

X_0_00 = 1;
X_0_m30 = 1 / eta_val^3;
X_0_m32 = 0;
X_0_01 = -e_val;
X_0_02 = (3*e_val^2 - 2 + 2*eta_val^3) / e_val^2;
X_0_03 = (8 - 12*e_val^2 + 3*e_val^4 - 8*eta_val^3) / e_val^3;

printf('Hansen coefficients at e = %.4f:\n', e_val);
printf('  X_0^{0,0}   = %.10f\n', X_0_00);
printf('  X_0^{-3,0}  = %.10f   (= 1/eta^3)\n', X_0_m30);
printf('  X_0^{-3,2}  = %.10e   (= 0)\n', X_0_m32);
printf('  X_0^{0,1}   = %.10f   (= -e)\n', X_0_01);
printf('  X_0^{0,2}   = %.10f\n', X_0_02);
printf('  X_0^{0,3}   = %.10f\n', X_0_03);
printf('\n');

%% Numerical averaging over l, parameterized by E:
%%   <F>_l = (1/(2pi)) * int_0^{2pi} F dl
%%         = (1/(2pi)) * int_0^{2pi} F * kappa dE
%% where F = F(E) is the integrand expressed in E via the geometric identities.

function avg = l_avg_via_E(e, integrand_of_E)
  %% integrand_of_E is a function handle F(E) that must accept vector E.
  %% Returns (1/(2pi)) * int_0^{2pi} F(E) * kappa(E) dE.
  %%
  %% Uses the periodic trapezoidal rule on a uniform grid E_k = 2*pi*k/N for
  %% k = 0, ..., N-1. For smooth 2pi-periodic integrands this converges
  %% exponentially; N = 4096 easily achieves machine precision.
  N = 4096;
  E_grid = (0:(N-1)) * (2*pi / N);
  kappa_vals = 1 - e * cos(E_grid);
  F_vals = integrand_of_E(E_grid);
  integrand_vals = F_vals .* kappa_vals;
  %% Periodic trapezoid: mean of function values times the period, divided by 2pi.
  avg = sum(integrand_vals) / N;
end

%% Auxiliaries: express cos(j*f) and sin(j*f) in terms of E.
%%   cos(f) = (cos(E) - e) / kappa
%%   sin(f) = eta * sin(E) / kappa
%% and iterate for j = 1, 2, 3 via the multiple-angle formulas.

function val = f_of_E(E, e)
  eta = sqrt(1 - e^2);
  kappa = 1 - e * cos(E);
  cf = (cos(E) - e) ./ kappa;
  sf = eta * sin(E) ./ kappa;
  val = atan2(sf, cf);
end

function val = cos_jf(E, e, j)
  val = cos(j * f_of_E(E, e));
end

function val = sin_jf(E, e, j)
  val = sin(j * f_of_E(E, e));
end

function val = kappa_of_E(E, e)
  val = 1 - e * cos(E);
end

%% Threshold for PASS
tol = 1e-9;
num_pass = 0;
num_fail = 0;

printf('-----------------------------------------------------------------\n');
printf('TEST 1: <1>_l = 1                                [X_0^{0,0} = 1]\n');
printf('-----------------------------------------------------------------\n');
num1 = l_avg_via_E(e_val, @(E) ones(size(E)));
exp1 = 1;                 % cos(0*g) * X_0^{0,0} = 1
err1 = abs(num1 - exp1);
printf('  numerical: %.15f\n', num1);
printf('  expected:  %.15f\n', exp1);
printf('  error:     %.3e\n', err1);
if err1 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 2: <cos(2g)>_l = cos(2g)                   [l-indep input]\n');
printf('-----------------------------------------------------------------\n');
num2 = l_avg_via_E(e_val, @(E) cos(2*g_val) * ones(size(E)));
exp2 = cos(2 * g_val) * X_0_00;
err2 = abs(num2 - exp2);
printf('  numerical: %.15f\n', num2);
printf('  expected:  %.15f\n', exp2);
printf('  error:     %.3e\n', err2);
if err2 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 3: <cos(2f + 2g)>_l = cos(2g) * X_0^{0,2}\n');
printf('-----------------------------------------------------------------\n');
integrand3 = @(E) cos(2 * f_of_E(E, e_val) + 2 * g_val);
num3 = l_avg_via_E(e_val, integrand3);
exp3 = cos(2 * g_val) * X_0_02;
err3 = abs(num3 - exp3);
printf('  numerical: %.15f\n', num3);
printf('  expected:  %.15f   (= cos(2g) * (3e^2-2+2eta^3)/e^2)\n', exp3);
printf('  error:     %.3e\n', err3);
if err3 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 4: <kappa^{-3}>_l = X_0^{-3,0} = 1/eta^3\n');
printf('-----------------------------------------------------------------\n');
integrand4 = @(E) 1 ./ kappa_of_E(E, e_val).^3;
num4 = l_avg_via_E(e_val, integrand4);
exp4 = X_0_m30;
err4 = abs(num4 - exp4);
printf('  numerical: %.15f\n', num4);
printf('  expected:  %.15f   (= 1/eta^3)\n', exp4);
printf('  error:     %.3e\n', err4);
if err4 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 5: <kappa^{-3} * cos(2f + 2g)>_l = cos(2g) * X_0^{-3,2} = 0\n');
printf('-----------------------------------------------------------------\n');
integrand5 = @(E) (1 ./ kappa_of_E(E, e_val).^3) .* cos(2 * f_of_E(E, e_val) + 2 * g_val);
num5 = l_avg_via_E(e_val, integrand5);
exp5 = cos(2 * g_val) * X_0_m32;
err5 = abs(num5 - exp5);
printf('  numerical: %.15e\n', num5);
printf('  expected:  %.15e   (= 0)\n', exp5);
printf('  error:     %.3e\n', err5);
if err5 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 6: <kappa^{-3} * sin(2f + 2g)>_l = sin(2g) * X_0^{-3,2} = 0\n');
printf('-----------------------------------------------------------------\n');
integrand6 = @(E) (1 ./ kappa_of_E(E, e_val).^3) .* sin(2 * f_of_E(E, e_val) + 2 * g_val);
num6 = l_avg_via_E(e_val, integrand6);
exp6 = sin(2 * g_val) * X_0_m32;
err6 = abs(num6 - exp6);
printf('  numerical: %.15e\n', num6);
printf('  expected:  %.15e   (= 0)\n', exp6);
printf('  error:     %.3e\n', err6);
if err6 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 7: <sin(3f)>_l = 0                        [pure k=0 sine]\n');
printf('-----------------------------------------------------------------\n');
integrand7 = @(E) sin(3 * f_of_E(E, e_val));
num7 = l_avg_via_E(e_val, integrand7);
exp7 = 0;                 % sin(0*g) * X_0^{0,3} = 0
err7 = abs(num7 - exp7);
printf('  numerical: %.15e\n', num7);
printf('  expected:  %.15e   (= 0 by reality of X_0^{0,3})\n', exp7);
printf('  error:     %.3e\n', err7);
if err7 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 8a: <sin(f + 2g)>_l = sin(2g) * X_0^{0,1} = sin(2g)*(-e)\n');
printf('         [ch07c Lemma C.4.8.1 at j=1]\n');
printf('-----------------------------------------------------------------\n');
integrand8a = @(E) sin(f_of_E(E, e_val) + 2 * g_val);
num8a = l_avg_via_E(e_val, integrand8a);
exp8a = sin(2 * g_val) * X_0_01;
err8a = abs(num8a - exp8a);
printf('  numerical: %.15f\n', num8a);
printf('  expected:  %.15f\n', exp8a);
printf('  error:     %.3e\n', err8a);
if err8a < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 8b: <sin(2f + 2g)>_l = sin(2g) * X_0^{0,2}\n');
printf('         [ch07c Lemma C.4.8.1 at j=2]\n');
printf('-----------------------------------------------------------------\n');
integrand8b = @(E) sin(2 * f_of_E(E, e_val) + 2 * g_val);
num8b = l_avg_via_E(e_val, integrand8b);
exp8b = sin(2 * g_val) * X_0_02;
err8b = abs(num8b - exp8b);
printf('  numerical: %.15f\n', num8b);
printf('  expected:  %.15f\n', exp8b);
printf('  error:     %.3e\n', err8b);
if err8b < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 9: <cos(2f + 2g)>_l sine-check: sine part sin(2g)*0 absent\n');
printf('        since <sin(2f)>_l = 0 and amplitude sum is cos(2g)*X_0^{0,2}\n');
printf('        (Verifies the angle-addition decomposition in Step 6.)\n');
printf('-----------------------------------------------------------------\n');
%% Evaluate both pieces:
%%   <cos(2f)>_l should equal X_0^{0,2}
%%   <sin(2f)>_l should equal 0
num9a = l_avg_via_E(e_val, @(E) cos(2 * f_of_E(E, e_val)));
num9b = l_avg_via_E(e_val, @(E) sin(2 * f_of_E(E, e_val)));
printf('  <cos(2f)>_l = %.15f  (expected X_0^{0,2} = %.15f)\n', ...
       num9a, X_0_02);
printf('  <sin(2f)>_l = %.15e  (expected 0)\n', num9b);
err9a = abs(num9a - X_0_02);
err9b = abs(num9b - 0);
printf('  errors: %.3e, %.3e\n', err9a, err9b);
if err9a < tol && err9b < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

printf('-----------------------------------------------------------------\n');
printf('TEST 10: <kappa^{-2} * cos(f)>_l = X_0^{-2,1} = 0   [B.0.7-4]\n');
printf('-----------------------------------------------------------------\n');
integrand10 = @(E) (1 ./ kappa_of_E(E, e_val).^2) .* cos(f_of_E(E, e_val));
num10 = l_avg_via_E(e_val, integrand10);
exp10 = 0;                %% cos(0*g) * X_0^{-2,1} = 0
err10 = abs(num10 - exp10);
printf('  numerical: %.15e\n', num10);
printf('  expected:  %.15e   (X_0^{-2,1} = 0 by B.0.7-4)\n', exp10);
printf('  error:     %.3e\n', err10);
if err10 < tol
  printf('  PASS\n\n');
  num_pass = num_pass + 1;
else
  printf('  FAIL\n\n');
  num_fail = num_fail + 1;
end

%% ================================================================
%% SUMMARY
%% ================================================================

printf('=================================================================\n');
printf('SUMMARY\n');
printf('=================================================================\n');
printf('Total tests:  %d\n', num_pass + num_fail);
printf('Passed:       %d\n', num_pass);
printf('Failed:       %d\n', num_fail);
if num_fail == 0
  printf('\n*** ALL TESTS PASSED — Theorem 5 (T5-c), (T5-s) verified ***\n');
else
  printf('\n*** %d TESTS FAILED — review Theorem 5 derivation ***\n', num_fail);
end
printf('=================================================================\n');
