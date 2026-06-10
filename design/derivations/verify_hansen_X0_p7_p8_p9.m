% verify_hansen_X0_p7_p8_p9.m
%
% Step 04 verifier: Hansen coefficients X_0^{-p, m} for p in {7, 8, 9},
%                   m in {0, 2, 4}  (9 pairs total).
%
% Companion to: sgp4_references/vallado_celestrak/documentation/SGP4/
%               Brouwer_Hori 1961/derivation/proof_hansen_X0_p7_p8_p9.md
%
% Five parts:
%   Part 1 - Wallis extension I_10..I_13 via TWO routes (parametric
%     differentiation of J_p(a, b), and the three-term recurrence
%     (p-1) eta^2 I_p = (2p-3) I_{p-1} - (p-2) I_{p-2} of proof_03 §2.6),
%     compared to closed forms of proof §2.5. Zero symbolic residual.
%
%   Part 2 - Symbolic assembly of each X_0^{-p, m} via the Wallis-route
%     bracket structure of proof §§3-5, compared to the boxed closed
%     forms of §§3.1-5.3. Zero residual (9 checks).
%
%   Part 3 - Auxiliary X_0^{-p, 1} for p in {6, 7, 8, 9} (proof §6):
%     derivation X_0^{-p, 1} = (1/(2 pi e))[eta^2 I_p - I_{p-1}] compared
%     to closed form. Zero residual (4 checks).
%
%   Part 4 - ch06e Theorem B.0.4 three-term recurrence at m = 0:
%     X_0^{-p, 0} + e X_0^{-p, 1} = eta^2 X_0^{-(p+1), 0}
%     for p in {6, 7, 8} (proof §7). Zero residual (3 checks).
%
%   Part 5 - Numerical cross-check at e in {0.01, 0.1, 0.3, 0.5, 0.7, 0.9}
%     for all 9 (p, m) pairs via periodic trapezoidal rule at N = 4096
%     (spectrally convergent on smooth 2pi-periodic integrands; see proof
%     §7 and OCTAVE_VERIFICATION §27-28). Relative tolerance 1e-12.
%     54 checks (9 pairs * 6 e-values).
%
% Expected total: 8 (Wallis x 2 routes at p = 10..13) + 9 (symbolic X_0)
%                 + 4 (auxiliary X_0^{-p, 1}) + 3 (B.0.4) + 54 (numerical)
%               = 78 PASSES.

pkg load symbolic;

printf('============================================================\n');
printf('Step 04 verifier: Hansen X_0^{-p, m}, p in {7,8,9}, m in {0,2,4}\n');
printf('============================================================\n\n');

n_pass = 0;
n_fail = 0;

%% =========================================================
%% Part 1: Wallis extension I_10..I_13.
%% =========================================================
printf('--- Part 1a: I_10..I_13 via parametric differentiation ---\n');

syms a b real positive;

J = cell(14, 1);
J{1} = 2 * sym(pi) / sqrt(a^2 - b^2);
J{2} = 2 * sym(pi) * a / (a^2 - b^2)^(sym(3)/2);
for p = 2:12
  J{p+1} = simplify(-diff(J{p}, a) / p);
end

syms e real positive;
eta = sqrt(1 - e^2);

I_pd = cell(14, 1);
for p = 1:13
  I_pd{p} = simplify(subs(J{p}, [a, b], [sym(1), e]));
end

printf('--- Part 1b: I_10..I_13 via three-term recurrence ---\n');

I_rec = cell(14, 1);
I_rec{1} = 2 * sym(pi) / eta;
I_rec{2} = 2 * sym(pi) / eta^3;
for p = 3:13
  I_rec{p} = simplify( ((2*p - 3) * I_rec{p-1} - (p - 2) * I_rec{p-2}) ...
                       / ((p - 1) * eta^2) );
end

printf('--- Part 1c: closed-form comparison for I_10..I_13 ---\n');

I_closed = cell(14, 1);
I_closed{5}  = sym(pi) * (8 + 24*e^2 + 3*e^4) / (4 * eta^9);
I_closed{6}  = sym(pi) * (8 + 40*e^2 + 15*e^4) / (4 * eta^11);
I_closed{7}  = sym(pi) * (16 + 120*e^2 + 90*e^4 + 5*e^6) / (8 * eta^13);
I_closed{8}  = sym(pi) * (16 + 168*e^2 + 210*e^4 + 35*e^6) / (8 * eta^15);
I_closed{9}  = sym(pi) * (128 + 1792*e^2 + 3360*e^4 + 1120*e^6 + 35*e^8) ...
               / (64 * eta^17);
I_closed{10} = sym(pi) * (128 + 2304*e^2 + 6048*e^4 + 3360*e^6 + 315*e^8) ...
               / (64 * eta^19);
I_closed{11} = sym(pi) * (256 + 5760*e^2 + 20160*e^4 + 16800*e^6 ...
                          + 3150*e^8 + 63*e^10) / (128 * eta^21);
I_closed{12} = sym(pi) * (256 + 7040*e^2 + 31680*e^4 + 36960*e^6 ...
                          + 11550*e^8 + 693*e^10) / (128 * eta^23);
I_closed{13} = sym(pi) * (1024 + 33792*e^2 + 190080*e^4 + 295680*e^6 ...
                          + 138600*e^8 + 16632*e^10 + 231*e^12) ...
               / (512 * eta^25);

for p = 10:13
  d_pd  = simplify(I_pd{p}  - I_closed{p});
  d_rec = simplify(I_rec{p} - I_closed{p});
  ok_pd  = logical(d_pd  == 0);
  ok_rec = logical(d_rec == 0);
  if ok_pd
    s_pd = 'PASS';
    n_pass = n_pass + 1;
  else
    s_pd = 'FAIL';
    n_fail = n_fail + 1;
  end
  if ok_rec
    s_rec = 'PASS';
    n_pass = n_pass + 1;
  else
    s_rec = 'FAIL';
    n_fail = n_fail + 1;
  end
  printf('  I_%-2d : param-diff %s   recurrence %s\n', p, s_pd, s_rec);
end
printf('\n');

%% =========================================================
%% Part 2: Symbolic X_0^{-p, m} assembly via Wallis route.
%% =========================================================
printf('--- Part 2: Symbolic X_0^{-p, m} assembly ---\n');

X_closed = cell(3, 3);
% rows: p = 7, 8, 9. cols: m = 0, 2, 4.
X_closed{1, 1} = (8 + 40*e^2 + 15*e^4) / (8 * eta^11);
X_closed{1, 2} = 5*e^2 * (2 + e^2) / (4 * eta^11);
X_closed{1, 3} = 5*e^4 / (16 * eta^11);
X_closed{2, 1} = (16 + 120*e^2 + 90*e^4 + 5*e^6) / (16 * eta^13);
X_closed{2, 2} = 15*e^2 * (16 + 16*e^2 + e^4) / (64 * eta^13);
X_closed{2, 3} = 3*e^4 * (10 + e^2) / (32 * eta^13);
X_closed{3, 1} = (16 + 168*e^2 + 210*e^4 + 35*e^6) / (16 * eta^15);
X_closed{3, 2} = 7*e^2 * (48 + 80*e^2 + 15*e^4) / (64 * eta^15);
X_closed{3, 3} = 7*e^4 * (10 + 3*e^2) / (32 * eta^15);

for pi_idx = 1:3
  p = 6 + pi_idx;
  for mi_idx = 1:3
    m = (mi_idx - 1) * 2;
    % Build integral via Wallis route (proof §§3-5).
    if m == 0
      intg = I_closed{p-1};
    elseif m == 2
      intg = (2*eta^4/e^2) * I_closed{p+1} ...
           - (4*eta^2/e^2) * I_closed{p} ...
           + ((2 - e^2)/e^2) * I_closed{p-1};
    elseif m == 4
      intg = (8/e^4) * (eta^8 * I_closed{p+3} - 4*eta^6 * I_closed{p+2} ...
                        + 6*eta^4 * I_closed{p+1} - 4*eta^2 * I_closed{p} ...
                        + I_closed{p-1}) ...
           - (8/e^2) * (eta^4 * I_closed{p+1} - 2*eta^2 * I_closed{p} ...
                        + I_closed{p-1}) ...
           + I_closed{p-1};
    end
    X_Wallis = simplify(intg / (2 * sym(pi)));
    d = simplify(X_Wallis - X_closed{pi_idx, mi_idx});
    if logical(d == 0)
      tag = 'PASS';
      n_pass = n_pass + 1;
    else
      tag = 'FAIL';
      n_fail = n_fail + 1;
    end
    printf('  X_0^{-%d, %d}  %s\n', p, m, tag);
  end
end
printf('\n');

%% =========================================================
%% Part 3: Auxiliary X_0^{-p, 1} (for B.0.4 cross-check).
%% =========================================================
printf('--- Part 3: Auxiliary X_0^{-p, 1} for p in {6, 7, 8, 9} ---\n');

X1_closed = cell(4, 1);
X1_closed{1} = e * (4 + 3*e^2) / (2 * eta^9);                          % p = 6
X1_closed{2} = 5*e * (8 + 12*e^2 + e^4) / (16 * eta^11);               % p = 7
X1_closed{3} = 3*e * (8 + 20*e^2 + 5*e^4) / (8 * eta^13);              % p = 8
X1_closed{4} = 7*e * (64 + 240*e^2 + 120*e^4 + 5*e^6) / (128 * eta^15);% p = 9

for pi_idx = 1:4
  p = 5 + pi_idx;  % p = 6, 7, 8, 9
  % Derivation: X_0^{-p, 1} = (1/(2 pi e))[eta^2 I_p - I_{p-1}]  (proof §6.1).
  intg1 = (eta^2/e) * I_closed{p} - (1/e) * I_closed{p-1};
  X1_derived = simplify(intg1 / (2 * sym(pi)));
  d = simplify(X1_derived - X1_closed{pi_idx});
  if logical(d == 0)
    tag = 'PASS';
    n_pass = n_pass + 1;
  else
    tag = 'FAIL';
    n_fail = n_fail + 1;
  end
  printf('  X_0^{-%d, 1}  %s\n', p, tag);
end
printf('\n');

%% =========================================================
%% Part 4: B.0.4 cross-check at m = 0 for p in {6, 7, 8}.
%% =========================================================
printf('--- Part 4: B.0.4 (m=0) at p in {6, 7, 8} ---\n');

X0_closed = cell(4, 1);
X0_closed{1} = (8 + 24*e^2 + 3*e^4) / (8 * eta^9);                     % p = 6
X0_closed{2} = X_closed{1, 1};                                         % p = 7
X0_closed{3} = X_closed{2, 1};                                         % p = 8
X0_closed{4} = X_closed{3, 1};                                         % p = 9

for pi_idx = 1:3
  p = 5 + pi_idx;  % p = 6, 7, 8
  LHS = X0_closed{pi_idx} + e * X1_closed{pi_idx};
  RHS = eta^2 * X0_closed{pi_idx + 1};
  d = simplify(LHS - RHS);
  if logical(d == 0)
    tag = 'PASS';
    n_pass = n_pass + 1;
  else
    tag = 'FAIL';
    n_fail = n_fail + 1;
  end
  printf('  B.0.4 at (n=-%d, m=0):  X_0^{-%d,0} + e X_0^{-%d,1} == eta^2 X_0^{-%d,0}  %s\n', ...
         p, p, p, p+1, tag);
end
printf('\n');

%% =========================================================
%% Part 5: Numerical cross-check.
%% =========================================================
printf('--- Part 5: Numerical cross-check (rel tol 1e-12, N=4096 trap) ---\n');

e_vals = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9];
m_vals = [0, 2, 4];
p_vals = [7, 8, 9];
tol = 1e-12;
N = 4096;
dE = 2 * pi / N;
E_nodes = (0:(N-1)) * dE;

for pi_idx = 1:3
  p = p_vals(pi_idx);
  for mi_idx = 1:3
    m = m_vals(mi_idx);
    for ei_idx = 1:6
      e_num = e_vals(ei_idx);
      eta_num = sqrt(1 - e_num^2);

      kap_n = 1 - e_num * cos(E_nodes);
      cosf  = (cos(E_nodes) - e_num) ./ kap_n;
      if m == 0
        cosmf = ones(size(E_nodes));
      elseif m == 2
        cosmf = 2 * cosf.^2 - 1;
      elseif m == 4
        cosmf = 8 * cosf.^4 - 8 * cosf.^2 + 1;
      end
      integrand_vals = kap_n.^(1 - p) .* cosmf;  % kap^{n+1} where n = -p.
      I_num = dE * sum(integrand_vals);
      X_num = I_num / (2 * pi);

      % Closed form (numerical).
      if p == 7 && m == 0
        X_exp = (8 + 40*e_num^2 + 15*e_num^4) / (8 * eta_num^11);
      elseif p == 7 && m == 2
        X_exp = 5*e_num^2 * (2 + e_num^2) / (4 * eta_num^11);
      elseif p == 7 && m == 4
        X_exp = 5 * e_num^4 / (16 * eta_num^11);
      elseif p == 8 && m == 0
        X_exp = (16 + 120*e_num^2 + 90*e_num^4 + 5*e_num^6) / (16 * eta_num^13);
      elseif p == 8 && m == 2
        X_exp = 15 * e_num^2 * (16 + 16*e_num^2 + e_num^4) / (64 * eta_num^13);
      elseif p == 8 && m == 4
        X_exp = 3 * e_num^4 * (10 + e_num^2) / (32 * eta_num^13);
      elseif p == 9 && m == 0
        X_exp = (16 + 168*e_num^2 + 210*e_num^4 + 35*e_num^6) / (16 * eta_num^15);
      elseif p == 9 && m == 2
        X_exp = 7 * e_num^2 * (48 + 80*e_num^2 + 15*e_num^4) / (64 * eta_num^15);
      elseif p == 9 && m == 4
        X_exp = 7 * e_num^4 * (10 + 3*e_num^2) / (32 * eta_num^15);
      end

      rel_err = abs(X_num - X_exp) / max(1, abs(X_exp));

      if rel_err < tol
        tag = 'PASS';
        n_pass = n_pass + 1;
      else
        tag = 'FAIL';
        n_fail = n_fail + 1;
      end
      printf('  p=%d m=%d e=%.2f: X_num=%.10e X_exp=%.10e rel_err=%.2e %s\n', ...
             p, m, e_num, X_num, X_exp, rel_err, tag);
    end
  end
end
printf('\n');

%% =========================================================
%% Summary.
%% =========================================================
printf('============================================================\n');
printf('Summary: %d pass, %d fail\n', n_pass, n_fail);
printf('============================================================\n');

if n_fail > 0
  error('verify_hansen_X0_p7_p8_p9: %d checks FAILED.', n_fail);
end
