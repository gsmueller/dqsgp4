% verify_ch10c_hansen_extensions.m
%
% Derivation and verification of three new Hansen-coefficient extensions
% needed for Phase C's computation of <(partial S_1/partial l)^2>_l:
%
%   B.0.7-9:   X_0^{-6, 0}(e) = (8 + 24 e^2 + 3 e^4) / (8 eta^9)
%   B.0.7-10:  X_0^{-6, 2}(e) = e^2 (6 + e^2) / (4 eta^9)
%   B.0.7-11:  X_0^{-6, 4}(e) = e^4 / (16 eta^9)
%
% and the symbolic assembly
%
%   <(partial S_1/partial l)^2>_l
%     = (mu^2 k_2 / L^3)^2 * [
%         (A^2 + B^2/2) X_0^{-6, 0}(e)
%         + 2 A B cos(2g) X_0^{-6, 2}(e)
%         + (B^2/2) cos(4g) X_0^{-6, 4}(e)
%         - A^2 / eta^6
%       ]
%
% where A = A(theta) = (3 theta^2 - 1)/2 and B = B(theta) = 3(1 - theta^2)/2.
%
% Derivation strategy: compute each X_0^{-6, m} by combining the Wallis
% integrals I_p := int_0^{2pi} dE/kap^p (ch06c) via the Chebyshev
% expansion cos(m f) = T_m((cos E - e)/kap). The key simplification
%   (cos E - e) = (eta^2 - kap)/e
% lets the integrand factor cleanly.
%
% Source theorems (each labelled (T) for theorem-application / (D) for
% definition-based substitution) are traced in the script comments.

pkg load symbolic;

printf('=======================================================\n');
printf('Hansen extensions (B.0.7-9/-10/-11) + ch10c C.4.e\n');
printf('=======================================================\n\n');

n_pass = 0;
n_fail = 0;

% --- Wallis integrals I_p = int_0^{2pi} dE / kap^p (ch06c + ch06e recurrence).
% --- Seeds: I_1 = 2pi/eta, I_2 = 2pi/eta^3 (ch06c).
% --- (T) ch06e Theorem B.0.4 three-term recurrence:
%         (p-1) eta^2 I_p = (2p-3) I_{p-1} - (p-2) I_{p-2}      [derived in ch06e]
% --- Build I_3 .. I_9 symbolically from the seeds and the recurrence.
syms e real positive;
eta = sqrt(1 - e^2);

I = cell(10, 1);
I{1} = 2*sym(pi)/eta;
I{2} = 2*sym(pi)/eta^3;
for p = 3:9
  I{p} = simplify(((2*p - 3)*I{p-1} - (p - 2)*I{p-2}) / ((p - 1)*eta^2));
end

% Check the known closed forms (ch06c Theorems B.3, B.4 and seed).
expected = { 2*sym(pi)/eta, ...
             2*sym(pi)/eta^3, ...
             sym(pi)*(2 + e^2)/eta^5, ...
             sym(pi)*(2 + 3*e^2)/eta^7, ...
             sym(pi)*(8 + 24*e^2 + 3*e^4)/(4*eta^9), ...
             sym(pi)*(8 + 40*e^2 + 15*e^4)/(4*eta^11), ...
             sym(pi)*(16 + 120*e^2 + 90*e^4 + 5*e^6)/(8*eta^13), ...
             sym(pi)*(16 + 168*e^2 + 210*e^4 + 35*e^6)/(8*eta^15), ...
             sym(pi)*(128 + 1792*e^2 + 3360*e^4 + 1120*e^6 + 35*e^8)/(64*eta^17) };
printf('--- Wallis integrals I_1 .. I_9 match recurrence ---\n');
for p = 1:9
  d = simplify(I{p} - expected{p});
  if logical(d == 0)
    printf('  I_%d PASS\n', p);
    n_pass = n_pass + 1;
  else
    printf('  I_%d FAIL: diff = %s\n', p, char(d));
    n_fail = n_fail + 1;
  end
end
printf('\n');

%% =======================================================
%% X_0^{-6, 0}(e) via direct decomposition.
%%
%% Definition (ch06a B.0.1 at k = 0, m = 0, n = -6; reality T B.0.2(i)):
%%   X_0^{-6, 0}(e) = (1/2pi) int_0^{2pi} (r/a)^{-6} dl
%%                  = (1/2pi) int_0^{2pi} kap^{-5} dE      [dl = kap dE]
%%                  = I_5 / (2pi).
%% =======================================================
printf('--- X_0^{-6, 0}(e) derivation + cross-check ---\n');
X_60 = I{5} / (2*sym(pi));
X_60 = simplify(X_60);

X_60_expected = (8 + 24*e^2 + 3*e^4) / (8*eta^9);
d = simplify(X_60 - X_60_expected);
if logical(d == 0)
  printf('  PASS: X_0^{-6, 0}(e) = (8 + 24 e^2 + 3 e^4)/(8 eta^9).\n');
  n_pass = n_pass + 1;
else
  printf('  FAIL: diff = %s\n', char(d));
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% X_0^{-6, 2}(e) via the identity
%%   cos(2f) = 2 cos^2 f - 1,   cos f = (cos E - e)/kap.
%% so  cos(2f) = 2 (cos E - e)^2 / kap^2  - 1.
%%
%% Writing  cos E - e = -(kap - eta^2)/e  (from  cos E = (1-kap)/e  and
%%   1 - e^2 = eta^2),
%% therefore  (cos E - e)^2 = (kap - eta^2)^2/e^2 = (kap^2 - 2 eta^2 kap + eta^4)/e^2.
%%
%% 2pi X_0^{-6, 2}(e) = int kap^{-5} cos(2f) dE
%%                    = 2 int (cos E - e)^2 / kap^7 dE   -   I_5
%%                    = (2/e^2) int (kap^2 - 2 eta^2 kap + eta^4) / kap^7 dE  -  I_5
%%                    = (2/e^2) [I_5 - 2 eta^2 I_6 + eta^4 I_7]              -  I_5.
%% =======================================================
printf('--- X_0^{-6, 2}(e) derivation + cross-check ---\n');
numerator_62 = (eta^4 * I{7} - 2*eta^2 * I{6} + I{5});
X_62 = (2/e^2) * numerator_62 / (2*sym(pi)) - I{5}/(2*sym(pi));
X_62 = simplify(X_62);

X_62_expected = e^2 * (6 + e^2) / (4*eta^9);
d = simplify(X_62 - X_62_expected);
if logical(d == 0)
  printf('  PASS: X_0^{-6, 2}(e) = e^2 (6 + e^2)/(4 eta^9).\n');
  n_pass = n_pass + 1;
else
  printf('  FAIL: diff = %s\n', char(d));
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% X_0^{-6, 4}(e) via cos(4f) = 8 cos^4 f - 8 cos^2 f + 1
%%               = 8 (cos E - e)^4 / kap^4  -  8 (cos E - e)^2 / kap^2  +  1.
%%
%% Using (cos E - e) = (eta^2 - kap)/e:
%%   (cos E - e)^4 = (kap - eta^2)^4 / e^4
%%                 = (eta^8 - 4 eta^6 kap + 6 eta^4 kap^2 - 4 eta^2 kap^3 + kap^4) / e^4.
%%
%% 2pi X_0^{-6, 4}(e) = 8 int (cos E - e)^4/kap^9 dE  -  8 int (cos E - e)^2/kap^7 dE  +  I_5
%%                    = (8/e^4) [eta^8 I_9 - 4 eta^6 I_8 + 6 eta^4 I_7 - 4 eta^2 I_6 + I_5]
%%                      - (8/e^2) [eta^4 I_7 - 2 eta^2 I_6 + I_5]
%%                      + I_5.
%% =======================================================
printf('--- X_0^{-6, 4}(e) derivation + cross-check ---\n');
num_64_A = (eta^8 * I{9} - 4*eta^6 * I{8} + 6*eta^4 * I{7} - 4*eta^2 * I{6} + I{5});
num_64_B = (eta^4 * I{7} - 2*eta^2 * I{6} + I{5});
X_64 = ( (8/e^4)*num_64_A - (8/e^2)*num_64_B + I{5} ) / (2*sym(pi));
X_64 = simplify(X_64);

X_64_expected = e^4 / (16*eta^9);
d = simplify(X_64 - X_64_expected);
if logical(d == 0)
  printf('  PASS: X_0^{-6, 4}(e) = e^4/(16 eta^9).\n');
  n_pass = n_pass + 1;
else
  printf('  FAIL: diff = %s\n', char(d));
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Numerical cross-check: compare the closed forms above against
%% direct periodic-trapezoidal integration of the integrand at several
%% e values.  (T) ch06a Definition B.0.1 gives the integral to evaluate.
%% =======================================================
printf('--- Numerical cross-check against trapezoidal quadrature ---\n');
e_vals = [0.05, 0.2, 0.4, 0.6, 0.8];
N_E = 8192;
Eg = linspace(0, 2*pi, N_E+1);
Eg = Eg(1:end-1);
dE_num = 2*pi / N_E;

printf('  e       X_0^{-6, 0}  (quad vs closed)    X_0^{-6, 2}  (quad vs closed)    X_0^{-6, 4}  (quad vs closed)\n');
for k = 1:length(e_vals)
  ee = e_vals(k);
  et = sqrt(1 - ee^2);
  kap_n = 1 - ee*cos(Eg);
  cosE_n = cos(Eg);
  cosf_n = (cosE_n - ee) ./ kap_n;
  cos2f_n = 2*cosf_n.^2 - 1;
  cos4f_n = 8*cosf_n.^4 - 8*cosf_n.^2 + 1;

  X60_quad = sum(kap_n.^(-5)) * dE_num / (2*pi);          % dl = kap dE, (r/a)^{-n} = kap^{-n}
  X62_quad = sum(kap_n.^(-5) .* cos2f_n) * dE_num / (2*pi);
  X64_quad = sum(kap_n.^(-5) .* cos4f_n) * dE_num / (2*pi);

  X60_closed = (8 + 24*ee^2 + 3*ee^4) / (8*et^9);
  X62_closed = ee^2 * (6 + ee^2) / (4*et^9);
  X64_closed = ee^4 / (16*et^9);

  % Use absolute + relative tolerance: at small e, X_0^{-6, 4} is O(e^4) ~ 4e-7
  % at e=0.05, so a tight relative tolerance amplifies machine-precision errors.
  % abs(diff) <= atol + rtol * abs(closed) passes a combined criterion.
  atol = 1e-14;  rtol = 1e-9;
  abs_60 = abs(X60_quad - X60_closed);
  abs_62 = abs(X62_quad - X62_closed);
  abs_64 = abs(X64_quad - X64_closed);
  e60 = abs_60 / max(abs(X60_closed), 1e-14);
  e62 = abs_62 / max(abs(X62_closed), 1e-14);
  e64 = abs_64 / max(abs(X64_closed), 1e-14);

  ok_60 = abs_60 <= atol + rtol*abs(X60_closed);
  ok_62 = abs_62 <= atol + rtol*abs(X62_closed);
  ok_64 = abs_64 <= atol + rtol*abs(X64_closed);
  for ok_check = [ok_60, ok_62, ok_64]
    if ok_check
      n_pass = n_pass + 1;
    else
      n_fail = n_fail + 1;
    end
  end
  printf('  %.2f   %11.4e vs %11.4e (rel %.2e)  %11.4e vs %11.4e (rel %.2e)  %11.4e vs %11.4e (rel %.2e)\n', ...
         ee, X60_quad, X60_closed, e60, X62_quad, X62_closed, e62, X64_quad, X64_closed, e64);
end
printf('\n');

%% =======================================================
%% <(partial S_1/partial l)^2>_l symbolic assembly.
%%
%% From ch09e Proposition E.6(a):
%%   partial S_1/partial l = (mu^2 k_2 / L^3) * [ (A + B cos 2(f+g))/kap^3 - A/eta^3 ].
%%
%% Squaring and expanding (A + B cos 2(f+g))^2 = A^2 + 2AB cos 2(f+g) + B^2 cos^2 2(f+g),
%% and using cos^2 x = (1 + cos 2x)/2:
%%
%%   (partial S_1/partial l)^2 = (mu^2 k_2/L^3)^2 * [
%%       (A^2 + B^2/2)/kap^6
%%     + 2AB cos 2(f+g)/kap^6
%%     + (B^2/2) cos 4(f+g)/kap^6
%%     - 2A^2/(eta^3 kap^3)
%%     - 2AB cos 2(f+g)/(eta^3 kap^3)
%%     + A^2/eta^6
%%   ].
%%
%% Applying (T) Theorem 5 (ch10_foundations_thm5.md) per harmonic:
%%   <cos(jf + kg) * kap^{-m}>_l = cos(k g) X_0^{-m, j}(e),  (T5-c).
%%
%% Substituting X_0^{-3, 0} = 1/eta^3 (B.0.7-1), X_0^{-3, 2} = 0 (B.0.7-2),
%% and the three new corollaries X_0^{-6, m} above:
%% =======================================================
printf('--- <(partial S_1/partial l)^2>_l symbolic assembly ---\n');
syms theta g_sym mu k2 L real;
A = (3*theta^2 - 1)/2;
B = 3*(1 - theta^2)/2;

inner = (A^2 + B^2/2)*X_60_expected ...
      + 2*A*B*cos(2*g_sym)*X_62_expected ...
      + (B^2/2)*cos(4*g_sym)*X_64_expected ...
      - A^2/eta^6;
dS1_dl_sq_avg = (mu^2 * k2 / L^3)^2 * inner;

printf('  Symbolic form:\n');
printf('    <(partial S_1/partial l)^2>_l = (mu^2 k_2/L^3)^2 * { (A^2+B^2/2) X_0^{-6,0}\n');
printf('                                                      + 2AB cos(2g) X_0^{-6,2}\n');
printf('                                                      + (B^2/2) cos(4g) X_0^{-6,4}\n');
printf('                                                      - A^2/eta^6 }\n');

%% =======================================================
%% Sanity check: at e -> 0, partial S_1/partial l -> (mu^2 k_2/L^3) * B cos 2(l+g)
%% (since A-term = 0 and the only f-dependence reduces to cos 2(l+g) at circular orbit).
%% Squared and averaged: (mu^2 k_2/L^3)^2 B^2/2.
%% Our formula at e = 0, eta = 1:
%%   inner = (A^2 + B^2/2)*1 + 0 + 0 - A^2 = B^2/2. <-- PASSES directly.
%% =======================================================
printf('\n  Circular-limit sanity check (e -> 0):\n');
inner_e0 = subs(inner, e, 0);
inner_e0 = simplify(inner_e0);
expected_e0 = B^2/2;
d_e0 = simplify(inner_e0 - expected_e0);
if logical(d_e0 == 0)
  printf('    PASS: inner(e=0) = B^2/2 = (3(1-theta^2)/2)^2/2.\n');
  n_pass = n_pass + 1;
else
  printf('    FAIL: diff = %s\n', char(d_e0));
  n_fail = n_fail + 1;
end

%% =======================================================
%% Numerical cross-check: compute <(partial S_1/partial l)^2>_l directly
%% by trapezoidal quadrature of (partial S_1/partial l)^2 and compare to
%% the symbolic closed form at a few (theta, e, g) points.
%% =======================================================
printf('\n  Numerical cross-check at sample (theta, e, g) points:\n');
printf('    theta    e     g        symbolic               trapezoidal           rel diff   Status\n');
sample = [0.2, 0.3, 0.5;
          0.4, 0.4, pi/3;
          0.447, 0.3, pi/4;
          0.6, 0.6, 2*pi/3;
          0.8, 0.5, pi];
N_l = 4096;
dl_n = 2*pi/N_l;
for k = 1:size(sample, 1)
  th = sample(k, 1);
  ee = sample(k, 2);
  gg = sample(k, 3);
  et = sqrt(1 - ee^2);
  A_num = (3*th^2 - 1)/2;
  B_num = 3*(1 - th^2)/2;

  % Symbolic value, evaluated numerically.
  X60 = (8 + 24*ee^2 + 3*ee^4) / (8*et^9);
  X62 = ee^2 * (6 + ee^2) / (4*et^9);
  X64 = ee^4 / (16*et^9);
  % Use L = mu = k_2 = 1 dimensionless so the prefactor (mu^2 k_2/L^3)^2 = 1.
  sym_val = (A_num^2 + B_num^2/2)*X60 + 2*A_num*B_num*cos(2*gg)*X62 + (B_num^2/2)*cos(4*gg)*X64 - A_num^2/et^6;

  % Trapezoidal quadrature of (partial S_1/partial l)^2 in l.
  quad_sum = 0;
  for iL = 0:N_l-1
    l_iL = iL * dl_n;
    % Kepler solve (Newton):
    E_k = l_iL;
    for iter = 1:50
      delta = (E_k - ee*sin(E_k) - l_iL) / (1 - ee*cos(E_k));
      E_k = E_k - delta;
      if abs(delta) < 1e-14; break; end
    end
    kap_k = 1 - ee*cos(E_k);
    cosf_k = (cos(E_k) - ee)/kap_k;
    sinf_k = et*sin(E_k)/kap_k;
    f_k = 2 * atan2(sqrt(1+ee)*sin(E_k/2), sqrt(1-ee)*cos(E_k/2));

    cos2fg = cos(2*(f_k + gg));
    dS1_dl = ((A_num + B_num*cos2fg)/kap_k^3 - A_num/et^3);   % (mu^2 k_2 / L^3) = 1
    quad_sum = quad_sum + dS1_dl^2;
  end
  quad_val = quad_sum / N_l;

  rel_diff = abs(sym_val - quad_val) / max(abs(sym_val), 1e-14);
  if rel_diff < 1e-10
    status = 'PASS';
    n_pass = n_pass + 1;
  else
    status = 'FAIL';
    n_fail = n_fail + 1;
  end
  printf('    %5.3f  %4.2f  %6.3f  %20.12e  %20.12e  %8.2e  %s\n', ...
         th, ee, gg, sym_val, quad_val, rel_diff, status);
end
printf('\n');

%% =======================================================
%% Summary
%% =======================================================
printf('=======================================================\n');
if n_fail == 0
  printf('ALL %d CHECKS PASSED.\n', n_pass);
else
  printf('FAILED: %d pass, %d fail.\n', n_pass, n_fail);
end
printf('=======================================================\n');
