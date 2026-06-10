% independent_X0_octave.m
% Clean-room independent symbolic derivation of three orbital integrals
%   M_j(e) := (1/(2 pi)) * integral_0^{2 pi} cos(j f) dl,  j = 1, 2, 3
% and the linear combination S(e) := 3 M_2 + 3 e M_1 + e M_3.
%
% Inputs (no project-specific framework values):
%   - Hansen integral definition; reality of M_j (M_j is real because
%     the integrand cos(jf) is real).
%   - Kepler measure relation: dl = (1 - e cos E) dE  =:  kappa dE.
%   - True-anomaly identities cos f = (cos E - e)/kappa, sin f = eta sin E/kappa.
%   - eta := sqrt(1 - e^2)  with the relation eta^2 = 1 - e^2.
%   - Standard Wallis integrals I_p := int_0^{2 pi} dE / kappa^p:
%        I_1 = 2 pi / eta
%        I_2 = 2 pi / eta^3
%        I_3 = pi (2 + e^2)   / eta^5
%        I_4 = pi (2 + 3 e^2) / eta^7
%        I_0 = 2 pi
%        I_{-1} = 2 pi
%        I_{-2} = 2 pi (1 + e^2 / 2)
%     These are textbook (chapter 0.D) results from residue calculus on
%     |z| = 1, NOT framework-specific.
%
% Method (PURE SYMBOLIC, no numerical approximation anywhere):
%
%   (1) For each j, build P_j(cos E) := Re(z^j), where
%         z = (cos E - e) + i eta sin E,
%       directly by binomial expansion.  Re(z^j) is a polynomial in cos E
%       (after the substitution sin^2 E = 1 - cos^2 E) of degree exactly j.
%   (2) Note the algebraic identity  kappa^j cos(j f) = Re(z^j) = P_j(cos E)
%       (a direct consequence of  kappa(cos f + i sin f) = z).
%       Therefore the integrand reduces to
%         kappa * cos(j f) = P_j(cos E) / kappa^{j-1}.
%   (3) Substitute  cos E = (1 - kappa)/e  (orbit equation in E) so that
%       P_j becomes a polynomial in kappa with coefficients in (e, eta).
%   (4) Integrate term by term using the standard Wallis I_p above.
%
% Output: closed forms for M_1, M_2, M_3 and a symbolic check of the
% cancellation identity  S(e) =? -M_2(e).

pkg load symbolic;
clear; close all;

syms e eta E real;
kappa = 1 - e * cos(E);

% Standard Wallis closed forms (textbook, NOT framework)
function out = I_p(p, e, eta)
  pi_s = sym(pi);
  switch p
    case -2,  out = 2*pi_s + pi_s * e^2;
    case -1,  out = 2*pi_s;
    case 0,   out = 2*pi_s;
    case 1,   out = 2*pi_s / eta;
    case 2,   out = 2*pi_s / eta^3;
    case 3,   out = pi_s * (2 + e^2) / eta^5;
    case 4,   out = pi_s * (2 + 3*e^2) / eta^7;
    otherwise
      error('I_p index %d not supplied', p);
  end
end

% Build P_j(cos E) := Re(z^j) for z = (cos E - e) + i eta sin E
% via the binomial theorem.  The result is a polynomial in cos E and (eta, e)
% after sin^2 E -> 1 - cos^2 E.
function out = build_Pj(j, e, eta, E)
  a = cos(E) - e;
  b = eta * sin(E);
  Re_zj = sym(0);
  for k = 0:j
    if mod(k, 4) == 0
      sgn =  1;
    elseif mod(k, 4) == 2
      sgn = -1;
    else
      continue;   % odd k contributes only to imaginary part
    end
    Re_zj = Re_zj + sgn * nchoosek(j, k) * a^(j - k) * b^k;
  end
  Re_zj = expand(Re_zj);
  % Eliminate sin^{2k} E iteratively.
  for it = 1:6
    prev = Re_zj;
    Re_zj = subs(Re_zj, sin(E)^2, 1 - cos(E)^2);
    Re_zj = expand(Re_zj);
    if isequal(Re_zj, prev)
      break;
    end
  end
  out = Re_zj;
end

function out = reduce_eta(expr_in, e, eta)
  out = expand(expr_in);
  for it = 1:8
    prev = out;
    out = subs(out, eta^2, 1 - e^2);
    out = expand(out);
    if isequal(out, prev)
      break;
    end
  end
  out = simplify(out);
end

% Robust equivalence test: returns true iff expr equals 0 in the quotient
% ring R(e)[eta] / (eta^2 - (1 - e^2)).
%
% Method: clear all denominators by multiplying by enough powers of e and eta;
% write the resulting numerator as A(e) + eta*B(e) (using eta^(2k) -> (1-e^2)^k);
% A and B must each vanish identically as polynomials in e.
function tf = vanishes_modulo_eta_sq(expr_in, e, eta)
  % Multiply by enough powers of e and eta to clear all denominators.
  cleared = expand(expr_in * e^10 * eta^8);
  % Reduce eta^(even) iteratively.
  for it = 1:12
    prev = cleared;
    cleared = subs(cleared, eta^2, 1 - e^2);
    cleared = expand(cleared);
    if isequal(cleared, prev)
      break;
    end
  end
  % After this, eta only appears to power 0 or 1.
  a_part = simplify(subs(cleared, eta, 0));
  b_part = simplify(subs(diff(cleared, eta), eta, 0));
  tf = isequal(a_part, sym(0)) && isequal(b_part, sym(0));
end

fprintf('=== Independent symbolic derivation: M_1, M_2, M_3 ===\n');
fprintf('Method: P_j := Re(z^j) by binomial; kappa-reduction; standard Wallis I_p\n');
fprintf('(Framework-claimed X_0^{0,m} values are NOT used.)\n\n');

results = cell(3, 1);
b_storage = cell(3, 1);

for j = 1:3
  fprintf('--- M_%d(e) ---\n', j);

  P_j = build_Pj(j, e, eta, E);
  fprintf('  P_%d(cos E) := Re(z^%d), reduced to polynomial in cos E:\n', j, j);
  fprintf('     %s\n', char(P_j));

  % Integrand = kappa * cos(j f) = P_j(cos E) / kappa^{j - 1}
  k_denom = j - 1;
  fprintf('  Integrand kappa * cos(%d f) = P_%d / kappa^%d\n', j, j, k_denom);

  % Substitute cos E -> (1 - kappa)/e to get polynomial in kappa
  syms kap;
  P_in_kap = expand(subs(P_j, cos(E), (1 - kap)/e));
  fprintf('  After cos E = (1 - kappa)/e:\n     %s\n', char(P_in_kap));

  % Extract coefficients b_q of kap^q via repeated differentiation at kap=0
  d_max = 0;
  current = P_in_kap;
  for trial = 0:10
    val0 = simplify(subs(current, kap, 0));
    if ~isequal(val0, sym(0))
      d_max = trial;
    end
    current = diff(current, kap);
  end

  fprintf('  b_q coefficients of kap^q (q = 0..%d):\n', d_max);
  b = cell(d_max + 1, 1);
  current = P_in_kap;
  for q = 0:d_max
    b{q+1} = simplify(subs(current, kap, 0) / factorial(q));
    fprintf('     b_%d = %s\n', q, char(b{q+1}));
    current = diff(current, kap);
  end
  b_storage{j} = b;

  % M_j = (1/(2 pi)) sum_q b_q I_{k_denom - q}
  Mj_int = sym(0);
  for q = 0:d_max
    Mj_int = Mj_int + b{q+1} * I_p(k_denom - q, e, eta);
  end
  Mj = simplify(Mj_int / (2 * sym(pi)));
  fprintf('  M_%d(e) = %s\n\n', j, char(Mj));

  results{j} = Mj;
end

M1 = results{1};
M2 = results{2};
M3 = results{3};

fprintf('===  Closed forms (in e and eta = sqrt(1 - e^2))  ===\n');
fprintf('  M_1(e) = %s\n', char(M1));
fprintf('  M_2(e) = %s\n', char(M2));
fprintf('  M_3(e) = %s\n\n', char(M3));

% --- Cancellation identity test ---
S_expr = simplify(3*M2 + 3*e*M1 + e*M3);
fprintf('===  S(e) := 3 M_2 + 3 e M_1 + e M_3  ===\n');
fprintf('  S(e) = %s\n\n', char(S_expr));

diff_S_negM2 = simplify(S_expr + M2);
fprintf('===  Cancellation identity test:  is  S(e) = -M_2(e) ?  ===\n');
fprintf('  S(e) + M_2(e) (raw)  = %s\n', char(diff_S_negM2));
if vanishes_modulo_eta_sq(diff_S_negM2, e, eta)
  fprintf('  *** PASS: S(e) = -M_2(e) identically (modulo eta^2 = 1 - e^2). ***\n');
  fprintf('  Independent verification of the cancellation identity:  CONFIRMED.\n');
else
  num_e = expand(numden(together(diff_S_negM2 + sym(0))));
  for it = 1:10
    num_e = expand(subs(num_e, eta^2, 1 - e^2));
  end
  fprintf('  FAIL: residual does not vanish.  Reduced numerator = %s\n', char(num_e));
end

% --- Independent check of two candidate closed forms ---
% These are NOT used in the derivation itself; they are tested for
% equivalence to the script's independently-derived M_2 and M_3.

candidate_M2 = (3*e^2 - 2 + 2*eta^3) / e^2;
diff_M2 = M2 - candidate_M2;
fprintf('\n===  Equivalence test:  M_2 == (3 e^2 - 2 + 2 eta^3)/e^2 ?  ===\n');
fprintf('  Difference (raw)  = %s\n', char(simplify(diff_M2)));
if vanishes_modulo_eta_sq(diff_M2, e, eta)
  fprintf('  PASS:  Independently-derived M_2 is equivalent to (3e^2 - 2 + 2 eta^3)/e^2.\n');
else
  fprintf('  FAIL:  Independent M_2 differs from this candidate.\n');
end

candidate_M3 = (8 - 12*e^2 + 3*e^4 - 8*eta^3) / e^3;
diff_M3 = M3 - candidate_M3;
fprintf('\n===  Equivalence test:  M_3 == (8 - 12e^2 + 3e^4 - 8 eta^3)/e^3 ?  ===\n');
fprintf('  Difference (raw)  = %s\n', char(simplify(diff_M3)));
if vanishes_modulo_eta_sq(diff_M3, e, eta)
  fprintf('  PASS:  Independently-derived M_3 is equivalent to (8 - 12e^2 + 3e^4 - 8 eta^3)/e^3.\n');
else
  fprintf('  FAIL:  Independent M_3 differs from this candidate.\n');
end
