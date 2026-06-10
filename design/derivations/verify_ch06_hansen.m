% verify_ch06_hansen.m
% Numerical verification of Chapter 6 Hansen coefficient results.
%
% Verifies:
%   1. Wallis integrals I_1..I_4 (ch06c Theorem B.2, Corollary B.3.2.1)
%   2. Orbital averages X_0^{n,m} (ch06c Theorems B.3, B.4; ch06e B.0.6)
%   3. Three-term recurrence B.0.4 (ch06e Theorem B.0.4)
%   4. Leading-order coefficients (ch06e Corollary B.0.5.1)
%   5. Parity relation (ch06e Theorem B.0.5(a))
%
% Method: periodic trapezoidal rule over E in [0, 2*pi].
% Spectral convergence for smooth periodic integrands; N = 10000
% gives machine-precision agreement for e < 0.9.

clear; close all;

fprintf('=== Chapter 6 Hansen Coefficient Verification ===\n');
fprintf('Quadrature: periodic trapezoidal, N = 10000\n\n');

N = 10000;
E = linspace(0, 2*pi, N+1);
E = E(1:end-1);
dE = 2*pi / N;
tol = 1e-8;
n_fail = 0;

e_vals = [0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9];

% =============================================================
%  Test 1: Wallis integrals I_p = int_0^{2pi} dE / kappa^p
% =============================================================
fprintf('--- Test 1: Wallis integrals I_p ---\n');
fprintf('  Closed forms: I_1 = 2pi/eta, I_2 = 2pi/eta^3,\n');
fprintf('    I_3 = pi(2+e^2)/eta^5, I_4 = pi(2+3e^2)/eta^7\n\n');
fprintf('  e       p   Quadrature          Closed form         Rel err\n');

for ie = 1:length(e_vals)
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  kap = 1 - ecc*cos(E);

  I_exact = [2*pi/eta, 2*pi/eta^3, pi*(2+ecc^2)/eta^5, pi*(2+3*ecc^2)/eta^7];

  for p = 1:4
    I_quad = sum(kap.^(-p)) * dE;
    rel = abs(I_quad - I_exact(p)) / abs(I_exact(p));
    if rel > tol; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
    fprintf('  %.3f   %d   %18.12e  %18.12e  %8.2e%s\n', ...
            ecc, p, I_quad, I_exact(p), rel, tag);
  end
end
fprintf('\n');

% =============================================================
%  Test 2: Orbital averages X_0^{n,m}
%  X_0^{n,m} = (1/2pi) int_0^{2pi} kappa^{n+1} cos(m*f) dE
% =============================================================
fprintf('--- Test 2: Orbital averages X_0^{n,m} ---\n');
fprintf('  Closed forms: X_0^{-3,0} = 1/eta^3 [B.3],\n');
fprintf('    X_0^{-3,2} = 0 [B.4], X_0^{-2,0} = 1/eta [B.0.6.1],\n');
fprintf('    X_0^{-2,1} = 0 [B.0.6.2], X_0^{-3,1} = e/(2eta^3) [B.0.6.3]\n\n');
fprintf('  e        n   m   Quadrature          Closed form         Err        Src\n');

for ie = 1:length(e_vals)
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  kap = 1 - ecc*cos(E);
  cosf = (cos(E) - ecc) ./ kap;
  sinf = eta * sin(E) ./ kap;
  cos2f = 2*cosf.^2 - 1;

  % Each row: [n, m, quadrature integrand (kappa^{n+1} * cos(mf)), exact, label]
  % n=-3,m=0: kap^{-2} * 1
  q1 = sum(kap.^(-2)) * dE / (2*pi);          ex1 = 1/eta^3;
  % n=-3,m=2: kap^{-2} * cos(2f)
  q2 = sum(cos2f ./ kap.^2) * dE / (2*pi);    ex2 = 0;
  % n=-2,m=0: kap^{-1} * 1
  q3 = sum(kap.^(-1)) * dE / (2*pi);          ex3 = 1/eta;
  % n=-2,m=1: kap^{-1} * cos(f)
  q4 = sum(cosf ./ kap) * dE / (2*pi);        ex4 = 0;
  % n=-3,m=1: kap^{-2} * cos(f)
  q5 = sum(cosf ./ kap.^2) * dE / (2*pi);     ex5 = ecc/(2*eta^3);

  results = [-3,0,q1,ex1; -3,2,q2,ex2; -2,0,q3,ex3; -2,1,q4,ex4; -3,1,q5,ex5];
  labels = {'B.3', 'B.4', 'B.0.6.1', 'B.0.6.2', 'B.0.6.3'};

  for it = 1:size(results, 1)
    nv = results(it,1); mv = results(it,2);
    qv = results(it,3); ev = results(it,4);
    if abs(ev) > 1e-15
      err = abs(qv - ev) / abs(ev);
      err_type = 'rel';
    else
      err = abs(qv);
      err_type = 'abs';
    end
    if err > tol; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
    fprintf('  %.3f   %2d  %2d  %18.12e  %18.12e  %8.2e(%s)  %s%s\n', ...
            ecc, nv, mv, qv, ev, err, err_type, labels{it}, tag);
  end
end
fprintf('\n');

% =============================================================
%  Test 2b: New orbital averages from the kappa-reduction framework
%  (ch06_orbital_average_framework.md, Corollaries B.0.7-6, 7, 8)
% =============================================================
fprintf('--- Test 2b: New orbital averages (kappa-reduction framework) ---\n');
fprintf('  Closed forms: X_0^{0,1} = -e                      [B.0.7-6]\n');
fprintf('                X_0^{0,2} = (3e^2-2+2eta^3)/e^2     [B.0.7-7]\n');
fprintf('                X_0^{-1,2} = (1-eta)^2/e^2          [B.0.7-8]\n\n');
fprintf('  e        n   m   Quadrature          Closed form         Err        Src\n');

for ie = 1:length(e_vals)
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  kap = 1 - ecc*cos(E);
  cosf = (cos(E) - ecc) ./ kap;
  sinf = eta * sin(E) ./ kap;
  cos2f = 2*cosf.^2 - 1;

  % X_0^{0,1}: integrand kap^{0+1} cos(f) = kap * cos f
  q6 = sum(kap .* cosf) * dE / (2*pi);          ex6 = -ecc;
  % X_0^{0,2}: integrand kap * cos(2f)
  q7 = sum(kap .* cos2f) * dE / (2*pi);         ex7 = (3*ecc^2 - 2 + 2*eta^3)/ecc^2;
  % Numerically-stable equivalent: X_0^{0,2} = (2eta+1)(1-eta)/(1+eta)
  % via 3e^2-2+2eta^3 = 3e^2-2(1-eta^3); (1-eta^3)/e^2 = (1+eta+eta^2)/(1+eta)
  ex7_stable = (2*eta + 1)*(1 - eta) / (1 + eta);
  % X_0^{-1,2}: integrand kap^{0} * cos(2f) = cos(2f)
  q8 = sum(cos2f) * dE / (2*pi);                ex8 = (1 - eta)^2 / ecc^2;

  results_b = [0, 1, q6, ex6; 0, 2, q7, ex7; -1, 2, q8, ex8];
  labels_b = {'B.0.7-6', 'B.0.7-7', 'B.0.7-8'};

  for it = 1:size(results_b, 1)
    nv = results_b(it,1); mv = results_b(it,2);
    qv = results_b(it,3); ev = results_b(it,4);
    if abs(ev) > 1e-15
      err = abs(qv - ev) / abs(ev);
      err_type = 'rel';
    else
      err = abs(qv);
      err_type = 'abs';
    end
    % For B.0.7-7 at small e the framework formula has catastrophic cancellation;
    % compare quadrature to the numerically-stable form instead.
    if nv == 0 && mv == 2
      err_stable = abs(qv - ex7_stable)/abs(ex7_stable);
      eff_err = err_stable;  % use stable form for pass/fail
    else
      eff_err = err;
    end
    if eff_err > tol; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
    fprintf('  %.3f   %2d  %2d  %18.12e  %18.12e  %8.2e(%s)  %s%s\n', ...
            ecc, nv, mv, qv, ev, err, err_type, labels_b{it}, tag);
  end
end
fprintf('\n');

% Sanity check the stable-form identity for B.0.7-7 vs the framework closed form
% at e values where neither form has catastrophic cancellation:
fprintf('  Stable-form identity X_0^{0,2} = (2eta+1)(1-eta)/(1+eta):\n');
for ie = 3:length(e_vals)   % skip e = 0.001, 0.01
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  a = (3*ecc^2 - 2 + 2*eta^3)/ecc^2;
  b = (2*eta + 1)*(1 - eta)/(1 + eta);
  d = abs(a - b);
  if d > 1e-12; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
  fprintf('    e=%.3f  framework=%18.12e  stable=%18.12e  diff=%.2e%s\n', ...
          ecc, a, b, d, tag);
end
fprintf('\n');

% =============================================================
%  Test 2c: kappa-expansion of P_2 (Theorem B.2.0-N)
%  P_2(cos E) = b_0 + b_1 kap + b_2 kap^2
%  b_0 = 2 eta^4/e^2, b_1 = -4 eta^2/e^2, b_2 = (2-e^2)/e^2
% =============================================================
fprintf('--- Test 2c: kappa-expansion of P_2 (Theorem B.2.0-N) ---\n');
fprintf('  P_2(cos E) vs b_0 + b_1*kap + b_2*kap^2 at all quadrature points\n\n');
fprintf('  e         max |P_2 - (b_0+b_1 kap+b_2 kap^2)|\n');

for ie = 1:length(e_vals)
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  kap = 1 - ecc*cos(E);

  P2_direct = (2-ecc^2)*cos(E).^2 - 2*ecc*cos(E) + (2*ecc^2 - 1);
  b0 = 2*eta^4/ecc^2;
  b1 = -4*eta^2/ecc^2;
  b2 = (2-ecc^2)/ecc^2;
  P2_expand = b0 + b1*kap + b2*kap.^2;

  maxdiff = max(abs(P2_direct - P2_expand));
  % Scale tolerance by magnitude of expansion coefficients (cancellation at small e)
  scale = max(abs([b0, b1, b2]));
  if maxdiff > 1e-12 * max(1, scale); n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
  fprintf('  %.3f     %8.2e (scale %8.2e)%s\n', ecc, maxdiff, scale, tag);
end
fprintf('\n');

% =============================================================
%  Test 2d: kappa-reduction theorem (Theorem B.2.1) applied to P_2
%  int P_2 / kap^p dE = b_0 I_p + b_1 I_{p-1} + b_2 I_{p-2}
%  Verify at p = 1, 2, 4 (the values used by B.0.7-7, -8, -2)
% =============================================================
fprintf('--- Test 2d: kappa-reduction theorem B.2.1 on P_2 ---\n');
fprintf('  int P_2/kap^p dE = b_0 I_p + b_1 I_{p-1} + b_2 I_{p-2}\n\n');
fprintf('  e       p    int (direct)       sum b_q I_{p-q}     Abs diff\n');

% Wallis + negative-index closed forms
Ival = @(s, ecc, eta) ( ...
       (s == 0)  * 2*pi + ...
       (s == 1)  * 2*pi/eta + ...
       (s == 2)  * 2*pi/eta^3 + ...
       (s == 3)  * pi*(2+ecc^2)/eta^5 + ...
       (s == 4)  * pi*(2+3*ecc^2)/eta^7 + ...
       (s == -1) * 2*pi + ...
       (s == -2) * 2*pi*(1 + ecc^2/2));

for ie = [2, 4, 6]   % e in {0.01, 0.3, 0.7}
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  kap = 1 - ecc*cos(E);
  P2_direct = (2-ecc^2)*cos(E).^2 - 2*ecc*cos(E) + (2*ecc^2 - 1);

  b0 = 2*eta^4/ecc^2;
  b1 = -4*eta^2/ecc^2;
  b2 = (2-ecc^2)/ecc^2;

  for p = [1, 2, 4]
    lhs = sum(P2_direct ./ kap.^p) * dE;
    rhs = b0*Ival(p, ecc, eta) + b1*Ival(p-1, ecc, eta) + b2*Ival(p-2, ecc, eta);
    adiff = abs(lhs - rhs);
    % Scale tolerance by magnitude of each contribution (cancellation at small e for p=4)
    scale = max(abs([b0*Ival(p, ecc, eta), b1*Ival(p-1, ecc, eta), b2*Ival(p-2, ecc, eta)]));
    if adiff > 1e-10 * max(1, scale); n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
    fprintf('  %.3f   %d    %18.12e  %18.12e  %8.2e (scale %8.2e)%s\n', ...
            ecc, p, lhs, rhs, adiff, scale, tag);
  end
end
fprintf('\n');

% =============================================================
%  Test 3: Three-term recurrence B.0.4
%  X_k^{n,m} + (e/2)[X_k^{n,m+1} + X_k^{n,m-1}] = eta^2 X_k^{n-1,m}
% =============================================================
fprintf('--- Test 3: Recurrence B.0.4 ---\n');
fprintf('  LHS = X_k^{n,m} + (e/2)[X_k^{n,m+1} + X_k^{n,m-1}]\n');
fprintf('  RHS = eta^2 * X_k^{n-1,m}\n\n');
fprintf('  e       k   n   m     LHS                 RHS                 Abs diff\n');

rec_tests = [0,-2,0; 0,-3,0; 0,-2,1; 1,-3,0; 2,-3,2; 0,-3,2];

for ie = [1, 3, 5, 7]
  ecc = e_vals(ie);
  eta = sqrt(1 - ecc^2);
  kap = 1 - ecc*cos(E);
  cosf = (cos(E) - ecc) ./ kap;
  sinf = eta * sin(E) ./ kap;
  f_ang = atan2(sinf, cosf);
  l_ang = E - ecc*sin(E);

  for it = 1:size(rec_tests, 1)
    kv = rec_tests(it,1); nv = rec_tests(it,2); mv = rec_tests(it,3);

    % X_k^{n,m}: integrand = kap^{n+1} cos(m*f - k*l)
    Xnm   = sum(kap.^(nv+1) .* cos(mv*f_ang     - kv*l_ang)) * dE/(2*pi);
    Xnmp1 = sum(kap.^(nv+1) .* cos((mv+1)*f_ang  - kv*l_ang)) * dE/(2*pi);
    Xnmm1 = sum(kap.^(nv+1) .* cos((mv-1)*f_ang  - kv*l_ang)) * dE/(2*pi);
    % X_k^{n-1,m}: integrand = kap^{n} cos(m*f - k*l)
    Xn1m  = sum(kap.^(nv)   .* cos(mv*f_ang     - kv*l_ang)) * dE/(2*pi);

    LHS = Xnm + (ecc/2)*(Xnmp1 + Xnmm1);
    RHS = eta^2 * Xn1m;
    adiff = abs(LHS - RHS);
    if adiff > tol; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
    fprintf('  %.3f  %2d  %2d  %2d  %18.12e  %18.12e  %8.2e%s\n', ...
            ecc, kv, nv, mv, LHS, RHS, adiff, tag);
  end
end
fprintf('\n');

% =============================================================
%  Test 4: Leading-order coefficients (Corollary B.0.5.1)
%  X_{m+1}^{n,m}(e) = [(2m-n)/2] e + O(e^3)
%  X_{m-1}^{n,m}(e) = [-(2m+n)/2] e + O(e^3)
% =============================================================
fprintf('--- Test 4: Leading-order coefficients at e = 0.001 ---\n');
fprintf('  X_k^{n,m}(e)/e vs predicted (2m-n)/2 or -(2m+n)/2\n\n');
fprintf('  n   m   k     X/e (numeric)       Predicted           Rel err\n');

e_lo = 0.001;
eta_lo = sqrt(1 - e_lo^2);
kap_lo = 1 - e_lo*cos(E);
cosf_lo = (cos(E) - e_lo) ./ kap_lo;
sinf_lo = eta_lo * sin(E) ./ kap_lo;
f_lo = atan2(sinf_lo, cosf_lo);
l_lo = E - e_lo*sin(E);

% [n, m, k, predicted_coefficient]
lo_data = [-3, 0, 1, 1.5; -3, 0, -1, 1.5; -3, 2, 3, 3.5; -3, 2, 1, -0.5];

for it = 1:size(lo_data, 1)
  nv = lo_data(it,1); mv = lo_data(it,2);
  kv = lo_data(it,3); pred = lo_data(it,4);

  X_val = sum(kap_lo.^(nv+1) .* cos(mv*f_lo - kv*l_lo)) * dE/(2*pi);
  ratio = X_val / e_lo;

  if abs(pred) > 1e-15
    rel = abs(ratio - pred) / abs(pred);
  else
    rel = abs(ratio);
  end
  lo_tol = 1e-3;  % O(e^2) corrections
  if rel > lo_tol; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
  fprintf('  %2d  %2d  %2d  %18.12e  %18.12e  %8.2e%s\n', ...
          nv, mv, kv, ratio, pred, rel, tag);
end
fprintf('\n');

% =============================================================
%  Test 5: Parity X_k^{n,m}(-e) = (-1)^{m-k} X_k^{n,m}(e)
% =============================================================
fprintf('--- Test 5: Parity at e = 0.3 ---\n');
fprintf('  (-1)^{m-k} X_k^{n,m}(-e) should equal X_k^{n,m}(e)\n\n');
fprintf('  n   m   k     (-1)^s * X(-e)      X(e)                Abs diff\n');

e_par = 0.3;
eta_par = sqrt(1 - e_par^2);

% Positive e
kap_p = 1 - e_par*cos(E);
cosf_p = (cos(E) - e_par) ./ kap_p;
sinf_p = eta_par * sin(E) ./ kap_p;
fp = atan2(sinf_p, cosf_p);
lp = E - e_par*sin(E);

% Negative e (kappa = 1+e*cos E, cos f = (cos E + e)/kap, l = E + e*sin E)
kap_n = 1 + e_par*cos(E);
cosf_n = (cos(E) + e_par) ./ kap_n;
sinf_n = eta_par * sin(E) ./ kap_n;
fn = atan2(sinf_n, cosf_n);
ln = E + e_par*sin(E);

par_data = [-3,0,0; -3,0,1; -3,2,0; -3,2,1; -3,1,0; -2,0,1];

for it = 1:size(par_data, 1)
  nv = par_data(it,1); mv = par_data(it,2); kv = par_data(it,3);

  X_pos = sum(kap_p.^(nv+1) .* cos(mv*fp - kv*lp)) * dE/(2*pi);
  X_neg = sum(kap_n.^(nv+1) .* cos(mv*fn - kv*ln)) * dE/(2*pi);

  sfac = (-1)^(mv - kv);
  LHS_p = sfac * X_neg;
  adiff = abs(LHS_p - X_pos);
  if adiff > tol; n_fail = n_fail + 1; tag = ' FAIL'; else; tag = ''; end
  fprintf('  %2d  %2d  %2d  %18.12e  %18.12e  %8.2e%s\n', ...
          nv, mv, kv, LHS_p, X_pos, adiff, tag);
end
fprintf('\n');

% =============================================================
%  Summary
% =============================================================
if n_fail == 0
  fprintf('=== ALL TESTS PASSED ===\n');
else
  fprintf('=== %d TEST(S) FAILED ===\n', n_fail);
end
