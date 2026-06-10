% verify_ch09a_dS1_dL_symbolic.m
%
% PURE SYMBOLIC verification of the key non-trivial algebraic identities
% in ch09a (the partial derivative dS_1/dL).  No numerical approximations.
%
% Checks:
%   (1) Lemma E.1.2:  d X_0^{0,2}/de = 2 e (eta + 2)/(1 + eta)^2
%   (2) Sub-lemma in §2.1: 2 + e cos f + kappa = 3 - e^2 kappa sin^2 f / eta^2
%       (the algebraic identity that makes E.2.1 reduce to two clean terms)
%   (3) Sum-to-product identity (used in §2.2):
%       cos(f + 2g) + cos(3f + 2g) = 2 cos(2(f+g)) cos(f)
%   (4) The combined dA/dL closed form (E.2.1):
%       dA/dL  =  3 eta^2 sin f / (L e kappa)  -  e sin^3 f / L
%       (with A := f - l + e sin f and the chain rule via D.2 + D.7 + D.8 + orbit equation)
%
% All identities are reduced to symbolic 0 in the appropriate quotient ring.

pkg load symbolic;
clear; close all;

syms e eta L l real;
syms c s real;        % c = cos E, s = sin E
syms cf sf real;      % cf = cos f, sf = sin f (constrained: cf = (c-e)/kappa, sf = eta s/kappa)
syms g real;
syms c2g s2g real;    % cos(2g), sin(2g)

kappa = 1 - e * c;
cos_f = (c - e) / kappa;
sin_f = (eta * s) / kappa;

% ---- Test (1): d X_0^{0,2} / de  =  2 e (eta + 2) / (1 + eta)^2 ----
fprintf('--- Test (1): Lemma E.1.2  d X_0^{0,2}/de = 2 e (eta+2)/(1+eta)^2 ---\n');

X02 = (3*e^2 - 2 + 2*eta^3) / e^2;
% Treat eta as independent; deta/de = -e/eta.
dX_de = diff(X02, e) + diff(X02, eta) * (-e/eta);
dX_de = simplify(dX_de);

candidate = 2*e*(eta + 2) / (1 + eta)^2;
diff_E12 = dX_de - candidate;

% Robust modulo-(eta^2 - (1-e^2)) test.  Substitute eta = sqrt(1-e^2)
% directly and let SymPy simplify.  Returns true iff the result is 0.
function tf = vanishes(expr_in, e, eta)
  expr_e = subs(expr_in, eta, sqrt(1 - e^2));
  expr_e = simplify(expr_e);
  tf = isequal(expr_e, sym(0));
end

if vanishes(diff_E12, e, eta)
  fprintf('  *** PASS: Lemma E.1.2 verified symbolically. ***\n\n');
else
  fprintf('  *** FAIL ***\n');
  fprintf('  diff = %s\n\n', char(simplify(diff_E12)));
end

% ---- Test (2): 2 + e cos f + kappa = 3 - e^2 kappa sin^2 f / eta^2 ----
fprintf('--- Test (2): Sub-lemma  2 + e cos f + kappa = 3 - e^2 kappa sin^2 f / eta^2 ---\n');

LHS = 2 + e * cos_f + kappa;
RHS = 3 - e^2 * kappa * sin_f^2 / eta^2;
diff_2 = LHS - RHS;
diff_2_clean = expand(diff_2 * kappa * eta^2);
% Reduce sin^2 E -> 1 - cos^2 E
for k = 1:8
  prev = diff_2_clean;
  diff_2_clean = expand(subs(diff_2_clean, s^2, 1 - c^2));
  if isequal(diff_2_clean, prev), break; end
end
% Reduce eta^2 -> 1 - e^2
for k = 1:8
  prev = diff_2_clean;
  diff_2_clean = expand(subs(diff_2_clean, eta^2, 1 - e^2));
  if isequal(diff_2_clean, prev), break; end
end
diff_2_clean = simplify(diff_2_clean);
fprintf('  (LHS - RHS) * kappa * eta^2, reduced = %s\n', char(diff_2_clean));
if isequal(diff_2_clean, sym(0))
  fprintf('  *** PASS: Sub-lemma verified symbolically. ***\n\n');
else
  fprintf('  *** FAIL ***\n\n');
end

% ---- Test (3): cos(f+2g) + cos(3f+2g) = 2 cos(2(f+g)) cos(f) ----
fprintf('--- Test (3): Sum-to-product  cos(f+2g) + cos(3f+2g) = 2 cos(2(f+g)) cos(f) ---\n');

% Build cos(f+2g), cos(3f+2g), cos(2(f+g)) symbolically using cf = cos_f, sf = sin_f.
% cos(j f + 2g) = cos(j f) cos(2g) - sin(j f) sin(2g).
% Use Chebyshev: cos(j f) = T_j(cos f), sin(j f) = sin(f) U_{j-1}(cos f).
function out = chebT(j, x)
  if j == 0, out = sym(1); return; end
  if j == 1, out = x; return; end
  Tprev = sym(1); Tcurr = x;
  for k = 2:j
    Tnext = 2*x*Tcurr - Tprev; Tprev = Tcurr; Tcurr = Tnext;
  end
  out = Tcurr;
end
function out = chebU(jm1, x)
  if jm1 < 0, out = sym(0); return; end
  if jm1 == 0, out = sym(1); return; end
  Uprev = sym(1); Ucurr = 2*x;
  for k = 2:jm1
    Unext = 2*x*Ucurr - Uprev; Uprev = Ucurr; Ucurr = Unext;
  end
  out = Ucurr;
end

cos_jf = cell(3,1); sin_jf = cell(3,1);
for j = 1:3
  cos_jf{j} = chebT(j, cos_f);
  sin_jf{j} = sin_f * chebU(j-1, cos_f);
end

cos_jfp2g = cell(3,1);
for j = 1:3
  cos_jfp2g{j} = cos_jf{j} * c2g - sin_jf{j} * s2g;
end

LHS3 = cos_jfp2g{1} + cos_jfp2g{3};
RHS3 = 2 * cos_jfp2g{2} * cos_f;
diff_3 = expand(LHS3 - RHS3);
diff_3_cleared = expand(diff_3 * kappa^4);
for k = 1:10
  diff_3_cleared = expand(subs(diff_3_cleared, s^2, 1 - c^2));
end
for k = 1:8
  diff_3_cleared = expand(subs(diff_3_cleared, eta^2, 1 - e^2));
end
diff_3_cleared = simplify(diff_3_cleared);
fprintf('  (LHS - RHS) * kappa^4, reduced = %s\n', char(diff_3_cleared));
if isequal(diff_3_cleared, sym(0))
  fprintf('  *** PASS: Sum-to-product identity verified. ***\n\n');
else
  fprintf('  *** FAIL ***\n\n');
end

% ---- Test (4): The combined dA/dL closed form (E.2.1) ----
fprintf('--- Test (4): dA/dL = 3 eta^2 sin f / (L e kappa) - e sin^3 f / L (E.2.1) ---\n');
% A := f - l + e sin f.  By chain rule:
%   dA/dL = (df/dL)(1 + e cos f) + (de/dL) sin f
% with df/dL = sin f (2 + e cos f) / (L e)  (D.8)
%      de/dL = eta^2 / (L e)                (D.2)
%      1 + e cos f = eta^2 / kappa          (orbit equation)
%
% We verify that this equals  3 eta^2 sin f / (L e kappa)  -  e sin^3 f / L.

dA_dL_chain = (sin_f * (2 + e * cos_f) / (L * e)) * (eta^2 / kappa) ...
            + (eta^2 / (L * e)) * sin_f;
dA_dL_target = 3 * eta^2 * sin_f / (L * e * kappa) - e * sin_f^3 / L;

diff_4 = dA_dL_chain - dA_dL_target;
% Multiply by L * e * kappa^3 * eta^2 to clear all denominators.
diff_4_cleared = expand(diff_4 * L * e * kappa^3 * eta^2);
% Iterate s^2 and eta^2 substitutions ALTERNATING until stable.
% Also catch the c^2 + s^2 = 1 pattern (via subs(c^2, 1 - s^2) then s^2 -> 1 - c^2).
for it = 1:40
  prev = diff_4_cleared;
  diff_4_cleared = expand(subs(diff_4_cleared, s^2, 1 - c^2));
  diff_4_cleared = expand(subs(diff_4_cleared, eta^2, 1 - e^2));
  if isequal(diff_4_cleared, prev), break; end
end
diff_4_cleared = simplify(diff_4_cleared);
fprintf('  (chain - target) * L e kappa^3 eta^2, reduced = %s\n', char(diff_4_cleared));
if isequal(diff_4_cleared, sym(0))
  fprintf('  *** PASS: dA/dL closed form (E.2.1) verified. ***\n\n');
else
  % residual involves (c^2 + s^2 - 1) which is identically 0; check by substituting.
  test = subs(diff_4_cleared, s^2, 1 - c^2);
  test = simplify(expand(test));
  fprintf('  After one more s^2 -> 1-c^2 substitution: %s\n', char(test));
  if isequal(test, sym(0))
    fprintf('  *** PASS (post-residual fix): dA/dL closed form (E.2.1) verified. ***\n\n');
  else
    fprintf('  *** FAIL ***\n\n');
  end
end

fprintf('\n=== verify_ch09a_dS1_dL_symbolic.m  COMPLETE. ===\n');
