% verify_ch07c_E1_fourier.m
% Numerical verification of Proposition C.4.8.3:
%
%   f - l + e*sin f = 2 * eta^3 * sum_{k>=1} (X_k^{-3, 0}(e) / k) * sin(k*l)
%
% (ch07c_mean_subtraction.md Proposition C.4.8.3, added per Ring 2b of the
%  Phase C corrective cycle.)
%
% Method:
%   LHS per sample (e, l): solve Kepler's equation l = E - e*sin(E) via
%     Newton iteration to tolerance 1e-14; compute f from
%     tan(f/2) = sqrt((1+e)/(1-e)) * tan(E/2); evaluate f - l + e*sin(f).
%
%   RHS per sample: evaluate 2 * eta^3 * sum_{k=1}^{K_trunc(e)} (X_k / k) * sin(k*l)
%     where X_k^{-3, 0}(e) is computed numerically from ch06a Definition B.0.1:
%       X_k^{-3, 0}(e) = (1/2pi) * integral_0^{2pi} kap^{-2} * cos(k * l') dE'
%     via periodic trapezoidal rule over E' at fixed e, with N_quad large enough
%     to avoid Nyquist aliasing (N_quad > 2 * K_trunc).
%
%     K_trunc(e) is chosen ADAPTIVELY: extended until the tail contribution
%     |X_k / k| * 2 * eta^3 is below 1e-13 for ~20 consecutive k. This is
%     necessary because the decay rate of X_k^{-3, 0}(e) is controlled by the
%     analyticity strip width of (a/r)^3 as a function of l, which shrinks as
%     e -> 1. For e = 0.8, ~K_trunc = 300 is needed to reach 1e-10 tolerance;
%     for e = 0.2, K_trunc = 30 suffices.
%
% Pass criterion: |LHS - RHS| < 1e-10 at all sample (e, l) points.
%
% Prerequisite theorems used (all (T)-admissible in-project results):
%   ch06a Definition B.0.1 (Hansen coefficient at general k, m, n).
%   ch06a Theorem B.0.2(i) (reality of Hansen coefficients).
%   ch06 Corollary B.0.7-1 (X_0^{-3, 0} = 1/eta^3).
%   ch06b Lemma B.1.2 (df/dl = eta/kap^2) and Kepler orbit equation.
%   ch07c Lemma C.4.8.2 (<f - l + e*sin f>_l = 0).
%
% (Per the Plan Audit Rule: this verifier is a cross-check, NOT a proof.
%  The proof is the six (T)/(D)-labeled steps in ch07c Proposition C.4.8.3.)

clear; close all;

fprintf('=== Proposition C.4.8.3 Verification ===\n');
fprintf('  f - l + e*sin f = 2 * eta^3 * sum_{k>=1} (X_k^{-3,0}(e)/k) * sin(k*l)\n');
fprintf('  Quadrature for X_k^{-3,0}: periodic trapezoidal in E (Kepler measure);\n');
fprintf('    N_quad chosen > 2 * K_trunc to avoid Nyquist aliasing.\n');
fprintf('  K_trunc(e) adaptive: extended until tail contribution < 1e-13.\n');
fprintf('  Pass tolerance: |LHS - RHS| < 1e-10.\n\n');

tol = 1e-10;
tail_threshold = 1e-13;
tail_consecutive = 20;  % consecutive small-tail steps required to truncate
K_max = 2000;  % hard upper bound
n_fail = 0;

% Sample grid
e_vals = [0.05, 0.2, 0.4, 0.6, 0.8];
l_samples = pi * [1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 15/8];

fprintf('  e      K_used   l              LHS                 RHS                 |Diff|     Status\n');

for ie = 1:length(e_vals)
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);

  % Build Hansen coefficients X_k^{-3, 0}(e) adaptively in k.
  %   X_k^{-3, 0}(e) = (1/2pi) * int_0^{2pi} kap^{-2} * cos(k*l') dE'
  % (l' = E' - e*sin E'; dl' = kap dE'; (a/r)^3 = 1/kap^3; the kap from dl'
  %  converts kap^{-3} (a/r)^3 integrand to kap^{-2} integrand in E' measure.)
  %
  % Quadrature: set N_quad = 4 * K_max to stay comfortably above Nyquist.
  N_quad = 4 * K_max;
  Eg = linspace(0, 2*pi, N_quad+1);
  Eg = Eg(1:end-1);
  dE = 2*pi / N_quad;
  kap_g = 1 - ecc*cos(Eg);
  l_g = Eg - ecc*sin(Eg);
  integrand_base = kap_g.^(-2);

  X = zeros(1, K_max);
  K_used = K_max;
  tail_count = 0;
  for k = 1:K_max
    X(k) = sum(integrand_base .* cos(k*l_g)) * dE / (2*pi);
    tail_contribution = abs(X(k) / k) * 2 * eta^3;
    if tail_contribution < tail_threshold
      tail_count = tail_count + 1;
      if tail_count >= tail_consecutive && k >= 30
        K_used = k;
        break;
      end
    else
      tail_count = 0;
    end
  end

  % Sanity check: X_0^{-3, 0} should equal 1/eta^3 (ch06 B.0.7-1).
  X0 = sum(integrand_base) * dE / (2*pi);
  X0_exact = 1/eta^3;
  if abs(X0 - X0_exact)/abs(X0_exact) > 1e-8
    fprintf('  WARNING: X_0^{-3,0} quadrature check failed at e=%.3f (got %.6e, exact %.6e)\n', ...
            ecc, X0, X0_exact);
    n_fail = n_fail + 1;
  end

  for il = 1:length(l_samples)
    l_s = l_samples(il);

    % LHS: solve Kepler via Newton for E, then convert to f.
    Es = l_s;  % initial guess
    for iter = 1:50
      num = Es - ecc*sin(Es) - l_s;
      den = 1 - ecc*cos(Es);
      delta = num / den;
      Es = Es - delta;
      if abs(delta) < 1e-14
        break;
      end
    end
    fs = 2 * atan2(sqrt(1+ecc)*sin(Es/2), sqrt(1-ecc)*cos(Es/2));

    LHS = fs - l_s + ecc*sin(fs);

    % RHS: 2 * eta^3 * sum_{k=1}^{K_used} (X_k / k) * sin(k * l_s)
    RHS = 0;
    for k = 1:K_used
      RHS = RHS + (X(k) / k) * sin(k * l_s);
    end
    RHS = 2 * eta^3 * RHS;

    adiff = abs(LHS - RHS);
    if adiff > tol
      n_fail = n_fail + 1;
      status = 'FAIL';
    else
      status = 'PASS';
    end
    fprintf('  %.3f  %6d   %10.6f  %18.12e  %18.12e  %8.2e  %s\n', ...
            ecc, K_used, l_s, LHS, RHS, adiff, status);
  end
end
fprintf('\n');

% Supplementary sanity check: zero-mean property (ch07c Lemma C.4.8.2).
fprintf('  Sanity check: <f - l + e*sin f>_l = 0 (ch07c Lemma C.4.8.2)\n');
for ie = 1:length(e_vals)
  ecc = e_vals(ie);
  % Evaluate f - l + e*sin f at a dense l-grid and average.
  Np = 4096;
  Egp = linspace(0, 2*pi, Np+1);
  Egp = Egp(1:end-1);
  lp = Egp - ecc*sin(Egp);
  fp = 2 * atan2(sqrt(1+ecc)*sin(Egp/2), sqrt(1-ecc)*cos(Egp/2));
  integrand_E = (fp - lp + ecc*sin(fp));
  % Convert E-grid average to l-grid average: <h>_l = (1/2pi) int h dl = (1/2pi) int h kap dE.
  kap_p = 1 - ecc*cos(Egp);
  mean_l = sum(integrand_E .* kap_p) * (2*pi/Np) / (2*pi);
  if abs(mean_l) > 1e-10
    n_fail = n_fail + 1;
    status = 'FAIL';
  else
    status = 'PASS';
  end
  fprintf('    e=%.3f   <f - l + e sin f>_l = %8.2e  %s\n', ecc, mean_l, status);
end
fprintf('\n');

% =============================================================
%  Summary
% =============================================================
if n_fail == 0
  fprintf('=== ALL TESTS PASSED ===\n');
  fprintf('  Proposition C.4.8.3 numerically verified at %d sample (e, l) points\n', ...
          length(e_vals) * length(l_samples));
  fprintf('  plus %d zero-mean sanity checks.\n', length(e_vals));
else
  fprintf('=== %d TEST(S) FAILED ===\n', n_fail);
end
