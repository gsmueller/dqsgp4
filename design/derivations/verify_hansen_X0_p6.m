% verify_hansen_X0_p6.m
%
% Step 03 verifier: Hansen coefficients X_0^{-6, m} for m in {0, 2, 4}.
%
% Companion to: sgp4_references/vallado_celestrak/documentation/SGP4/
%               Brouwer_Hori 1961/derivation/proof_hansen_X0_p6.md
%
% Three parts:
%   Part 1 - Symbolic Wallis-chain verification:
%     (a) Build I_1 .. I_9 via parametric differentiation of J_p(a, b)
%         following ch06c Proposition B.3.2: J_{p+1} = -(1/p) dJ_p/da,
%         then specialize to (a, b) = (1, e) (this is the primary route of
%         the proof file §2.1 - §2.5).
%     (b) Independently build I_3 .. I_9 via the three-term recurrence
%         (p-1) eta^2 I_p = (2p-3) I_{p-1} - (p-2) I_{p-2}
%         (proof_hansen_X0_p6.md §2.6-R).
%     (c) Compare each of (a), (b) output against the closed-form table
%         of proof §2.6; zero symbolic residual required at each p.
%
%   Part 2 - Symbolic X_0^{-6, m} assembly:
%     For each m in {0, 2, 4}: build the Hansen-coefficient integral
%     as a linear combination of I_p (per proof §§3-5) and compare
%     against the boxed closed form; zero symbolic residual required.
%
%   Part 3 - Numerical quadrature cross-check (periodic trapezoidal):
%     For each m in {0, 2, 4} and each e in {0.01, 0.1, 0.3, 0.5, 0.8}:
%     evaluate the Hansen-integral definition
%       X_0^{-6, m}(e) = (1/2pi) * int_0^{2pi} kap^{-5} cos(m f(E, e)) dE
%     via the periodic trapezoidal rule at N = 4096 equispaced nodes.
%     For smooth 2pi-periodic integrands this rule is spectrally
%     convergent (Trefethen & Weideman 2014, "The exponentially convergent
%     trapezoidal rule"), reaching machine precision by N ~ 64-256 on
%     these integrands. Compare to closed form; require relative error
%     < 1e-12. `quadgk` is avoided here because adaptive subdivision
%     misdiagnoses tolerance on oscillatory periodic integrands (Gauss-
%     Kronrod's error estimator is not reliable for them).
%
% cos(m f) is evaluated by Chebyshev/double-angle recurrence from cos f =
% (cos E - e)/kap (wrap-immune trig-periodic use; see OCTAVE_VERIFICATION
% skill §28: cos(m theta) is wrap-immune, so no atan2 reconstruction of f
% is needed).
%
% Expected output: "Summary: 36 pass, 0 fail".
%   Part 1: 18 (9 I_p values x 2 routes).
%   Part 2:  3 (one per m).
%   Part 3: 15 (3 m-values x 5 e-values).

pkg load symbolic;

printf('========================================================\n');
printf('Step 03 verifier: Hansen X_0^{-6, m} (m = 0, 2, 4)\n');
printf('========================================================\n\n');

n_pass = 0;
n_fail = 0;

%% =========================================================
%% Part 1: Wallis-integral extension I_5 .. I_9.
%% =========================================================

%% ---- Part 1a: parametric-differentiation route. -------------
printf('--- Part 1a: I_1..I_9 via parametric differentiation ---\n');

syms a b real positive;

% J_1(a, b) and J_2(a, b) from ch06c Corollary B.3.2.1.
J = cell(9, 1);
J{1} = 2 * sym(pi) / sqrt(a^2 - b^2);
J{2} = 2 * sym(pi) * a / (a^2 - b^2)^(sym(3)/2);

% Iterate J_{p+1} = -(1/p) dJ_p/da.
for p = 2:8
  J{p+1} = simplify(-diff(J{p}, a) / p);
end

% Specialize to (a, b) = (1, e).
syms e real positive;
eta = sqrt(1 - e^2);

I_pd = cell(9, 1);
for p = 1:9
  I_pd{p} = simplify(subs(J{p}, [a, b], [sym(1), e]));
end

%% ---- Part 1b: three-term-recurrence route. -------------
printf('--- Part 1b: I_3..I_9 via three-term recurrence ---\n');

I_rec = cell(9, 1);
I_rec{1} = 2 * sym(pi) / eta;
I_rec{2} = 2 * sym(pi) / eta^3;
for p = 3:9
  I_rec{p} = simplify( ((2*p - 3) * I_rec{p-1} - (p - 2) * I_rec{p-2}) ...
                       / ((p - 1) * eta^2) );
end

%% ---- Part 1c: closed-form comparison. -------------
printf('--- Part 1c: closed-form comparison for I_1..I_9 ---\n');

I_closed = cell(9, 1);
I_closed{1} = 2 * sym(pi) / eta;
I_closed{2} = 2 * sym(pi) / eta^3;
I_closed{3} = sym(pi) * (2 + e^2) / eta^5;
I_closed{4} = sym(pi) * (2 + 3*e^2) / eta^7;
I_closed{5} = sym(pi) * (8 + 24*e^2 + 3*e^4) / (4 * eta^9);
I_closed{6} = sym(pi) * (8 + 40*e^2 + 15*e^4) / (4 * eta^11);
I_closed{7} = sym(pi) * (16 + 120*e^2 + 90*e^4 + 5*e^6) / (8 * eta^13);
I_closed{8} = sym(pi) * (16 + 168*e^2 + 210*e^4 + 35*e^6) / (8 * eta^15);
I_closed{9} = sym(pi) * (128 + 1792*e^2 + 3360*e^4 + 1120*e^6 + 35*e^8) ...
              / (64 * eta^17);

for p = 1:9
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
  printf('  I_%d : param-diff %s   recurrence %s\n', p, s_pd, s_rec);
end
printf('\n');

%% =========================================================
%% Part 2: X_0^{-6, m} assembly via the Wallis route.
%% =========================================================
printf('--- Part 2: Symbolic X_0^{-6, m} assembly ---\n');

I5 = I_closed{5};
I6 = I_closed{6};
I7 = I_closed{7};
I8 = I_closed{8};
I9 = I_closed{9};

% ----- m = 0 -----
% Integrand: P_0 kap^{-5} = kap^{-5}; integral = I_5.
X60_Wallis = I5 / (2 * sym(pi));
X60_closed = (8 + 24*e^2 + 3*e^4) / (8 * eta^9);
d60 = simplify(X60_Wallis - X60_closed);
if logical(d60 == 0)
  printf('  m = 0: X_0^{-6,0} = (8+24 e^2+3 e^4)/(8 eta^9)   PASS\n');
  n_pass = n_pass + 1;
else
  printf('  m = 0: FAIL, diff = %s\n', char(d60));
  n_fail = n_fail + 1;
end

% ----- m = 2 -----
% From proof §4 Step 2: int P_2/kap^7 dE =
%   (2 eta^4/e^2) I_7 - (4 eta^2/e^2) I_6 + ((2 - e^2)/e^2) I_5.
intP2_k7 = (2*eta^4/e^2) * I7 ...
         - (4*eta^2/e^2) * I6 ...
         + ((2 - e^2)/e^2) * I5;
X62_Wallis = intP2_k7 / (2 * sym(pi));
X62_closed = e^2 * (6 + e^2) / (4 * eta^9);
d62 = simplify(X62_Wallis - X62_closed);
if logical(d62 == 0)
  printf('  m = 2: X_0^{-6,2} = e^2 (6+e^2)/(4 eta^9)         PASS\n');
  n_pass = n_pass + 1;
else
  printf('  m = 2: FAIL, diff = %s\n', char(d62));
  n_fail = n_fail + 1;
end

% ----- m = 4 -----
% From proof §5 Step 6:
%   int P_4 kap^{-9} dE = (8/e^4) Sig_a - (8/e^2) Sig_b + I_5,
% where Sig_a = eta^8 I_9 - 4 eta^6 I_8 + 6 eta^4 I_7 - 4 eta^2 I_6 + I_5,
%       Sig_b = eta^4 I_7 - 2 eta^2 I_6 + I_5.
Sig_a = eta^8*I9 - 4*eta^6*I8 + 6*eta^4*I7 - 4*eta^2*I6 + I5;
Sig_b = eta^4*I7 - 2*eta^2*I6 + I5;
intP4_k9 = (8/e^4) * Sig_a - (8/e^2) * Sig_b + I5;
X64_Wallis = intP4_k9 / (2 * sym(pi));
X64_closed = e^4 / (16 * eta^9);
d64 = simplify(X64_Wallis - X64_closed);
if logical(d64 == 0)
  printf('  m = 4: X_0^{-6,4} = e^4/(16 eta^9)                PASS\n');
  n_pass = n_pass + 1;
else
  printf('  m = 4: FAIL, diff = %s\n', char(d64));
  n_fail = n_fail + 1;
end
printf('\n');

%% =========================================================
%% Part 3: Numerical adaptive quadrature cross-check.
%% =========================================================
printf('--- Part 3: Numerical cross-check (rel tol 1e-12, N=4096 trap) ---\n');

e_vals = [0.01, 0.1, 0.3, 0.5, 0.8];
m_vals = [0, 2, 4];
tol = 1e-12;
N = 4096;   % spectrally convergent at ~N = 64 for these integrands;
            % 4096 leaves ample headroom for the round-off floor.
dE = 2 * pi / N;
E_nodes = (0:(N-1)) * dE;  % skip duplicate endpoint for periodic rule.

for mi = 1:length(m_vals)
  m = m_vals(mi);
  for ei = 1:length(e_vals)
    e_num = e_vals(ei);
    eta_num = sqrt(1 - e_num^2);

    % Integrand kap^{-5} * cos(m f). cos(m f) via double-angle recurrence
    % from cos f = (cos E - e)/kap (wrap-immune trig-periodic use; no
    % atan2 reconstruction of f is needed). See OCTAVE_VERIFICATION §28.
    kap_n   = 1 - e_num * cos(E_nodes);
    cosf    = (cos(E_nodes) - e_num) ./ kap_n;
    if m == 0
      cosmf = ones(size(E_nodes));
    elseif m == 2
      cosmf = 2 * cosf.^2 - 1;
    elseif m == 4
      cosmf = 8 * cosf.^4 - 8 * cosf.^2 + 1;
    end
    integrand_vals = kap_n.^(-5) .* cosmf;

    I_num = dE * sum(integrand_vals);  % periodic trapezoidal rule.
    X_num = I_num / (2 * pi);

    % Closed form.
    if m == 0
      X_exp = (8 + 24*e_num^2 + 3*e_num^4) / (8 * eta_num^9);
    elseif m == 2
      X_exp = e_num^2 * (6 + e_num^2) / (4 * eta_num^9);
    elseif m == 4
      X_exp = e_num^4 / (16 * eta_num^9);
    end

    rel_err = abs(X_num - X_exp) / max(1, abs(X_exp));

    if rel_err < tol
      tag = 'PASS';
      n_pass = n_pass + 1;
    else
      tag = 'FAIL';
      n_fail = n_fail + 1;
    end
    printf('  m=%d, e=%.2f: X_num=%.16e  X_exp=%.16e  rel_err=%.2e  %s\n', ...
           m, e_num, X_num, X_exp, rel_err, tag);
  end
end
printf('\n');

%% =========================================================
%% Summary.
%% =========================================================
printf('========================================================\n');
printf('Summary: %d pass, %d fail\n', n_pass, n_fail);
printf('========================================================\n');

if n_fail > 0
  error('verify_hansen_X0_p6: %d checks FAILED.', n_fail);
end
