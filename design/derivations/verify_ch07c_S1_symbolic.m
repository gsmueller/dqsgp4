% verify_ch07c_S1_symbolic.m
%
% PURE SYMBOLIC verification of Theorem C.4 (closed form of S_1) and
% Theorem C.5 (zero mean + homological equation satisfaction).
%
% No numerical approximations.  Every identity is reduced to symbolic 0
% in the polynomial quotient ring
%   Q[e, eta, c, s, c2g, s2g] / { eta^2 - (1 - e^2),  s^2 - (1 - c^2),
%                                  c2g^2 + s2g^2 - 1 }
% where c = cos E, s = sin E, c2g = cos(2g), s2g = sin(2g).
%
% Theorem C.4 (claim, modulo prefactor mu*k_2/(n*a^3*eta^3)):
%   X(E, g; e, eta) = A * (f - l + e sin f)
%                   + (B/6) * [ 3 sin(2f+2g) + 3e sin(f+2g) + e sin(3f+2g)
%                              + X02 * sin(2g) ]
% where X02 = (3 e^2 - 2 + 2 eta^3)/e^2.
%
% C.5(b) homological equation:
%   d X / d l  ==  A * ((a/r)^3 eta^3 - 1)  +  B * (a/r)^3 eta^3 cos(2(f+g))
%
% C.5(a) zero mean:  <X>_l = 0.
%
% Strategy:
%   - Use cos(jf) = T_j(cos f), sin(jf) = sin(f) U_{j-1}(cos f) with
%     cos f = (cos E - e)/kappa, sin f = (eta sin E)/kappa.
%     Substituting these gives ALL trig-of-f in terms of (cos E, sin E, e, eta).
%   - All d/dl operations reduce via d/dl = (1/kappa) d/dE.
%   - Express the difference as a single rational function in
%     (cos E, sin E, e, eta, c2g, s2g), clear denominators, reduce.

pkg load symbolic;
clear; close all;

% Use lowercase placeholders so SymPy diff() won't be asked to differentiate
% w.r.t. sin(E) -- we'll substitute back at the end if needed.
syms e eta c s c2g s2g real;     % c = cos(E), s = sin(E), c2g = cos(2g), s2g = sin(2g)
syms kap real;                    % stand-in for kappa (we will substitute)

% Direct definitions in (c, s, e, eta) -- no SymPy trig functions involved.
kappa = 1 - e * c;
cos_f = (c - e) / kappa;
sin_f = (eta * s) / kappa;

% Chebyshev T_j(x) and U_{jm1}(x)
function out = chebT(j, x)
  Tprev = sym(1);  Tcurr = x;
  if j == 0,  out = Tprev;  return;  end
  if j == 1,  out = Tcurr;  return;  end
  for k = 2:j
    Tnext = 2*x*Tcurr - Tprev;  Tprev = Tcurr;  Tcurr = Tnext;
  end
  out = Tcurr;
end

function out = chebU(jm1, x)
  if jm1 < 0,  out = sym(0);  return;  end
  Uprev = sym(1);  Ucurr = 2 * x;
  if jm1 == 0,  out = Uprev;  return;  end
  for k = 2:jm1
    Unext = 2*x*Ucurr - Uprev;  Uprev = Ucurr;  Ucurr = Unext;
  end
  out = Ucurr;
end

% Reduce s^2 -> 1 - c^2 iteratively
function out = reduce_sin(expr_in, s, c)
  out = expand(expr_in);
  for k = 1:10
    prev = out;
    out = expand(subs(out, s^2, 1 - c^2));
    if isequal(out, prev),  break;  end
  end
end

% Reduce eta^2 -> 1 - e^2 iteratively
function out = reduce_eta(expr_in, eta, e)
  out = expand(expr_in);
  for k = 1:10
    prev = out;
    out = expand(subs(out, eta^2, 1 - e^2));
    if isequal(out, prev),  break;  end
  end
end

% Reduce c2g^2 -> 1 - s2g^2 iteratively
function out = reduce_g(expr_in, c2g, s2g)
  out = expand(expr_in);
  for k = 1:6
    prev = out;
    out = expand(subs(out, c2g^2, 1 - s2g^2));
    if isequal(out, prev),  break;  end
  end
end

% Robust full-reduction: clear denominators, reduce s^2, eta^2, c2g^2, then
% extract coefficients of (eta^0 vs eta^1) * (s^0 vs s^1) * (c2g^0 vs c2g^1).
% (8 components total.)  All must vanish identically as polynomials in (e, c, s2g).
function tf = vanishes_full(expr_in, e, eta, c, s, c2g, s2g)
  cleared = expand(expr_in * e^15 * eta^10 * (1 - e * c)^15);
  for it = 1:25
    prev = cleared;
    cleared = expand(subs(cleared, s^2, 1 - c^2));
    cleared = expand(subs(cleared, eta^2, 1 - e^2));
    cleared = expand(subs(cleared, c2g^2, 1 - s2g^2));
    if isequal(cleared, prev),  break;  end
  end
  % Decompose into components.  Iterate over (eta_pow, s_pow, c2g_pow) in {0, 1}.
  tf = true;
  for ea = 0:1
    for sa = 0:1
      for ga = 0:1
        % Project onto component eta^ea * s^sa * c2g^ga
        if ea == 0
          xe = subs(cleared, eta, 0);
        else
          xe = subs(diff(cleared, eta), eta, 0);
        end
        if sa == 0
          xs = subs(xe, s, 0);
        else
          xs = subs(diff(xe, s), s, 0);
        end
        if ga == 0
          xg = subs(xs, c2g, 0);
        else
          xg = subs(diff(xs, c2g), c2g, 0);
        end
        coef = simplify(xg);
        if ~isequal(coef, sym(0))
          tf = false;
          fprintf('    component (eta^%d s^%d c2g^%d) = %s\n', ea, sa, ga, char(coef));
        end
      end
    end
  end
end

fprintf('=== Theorems C.4 + C.5 verification (PURE SYMBOLIC) ===\n');
fprintf('Quotient ring Q[e, eta, c, s, c2g, s2g] / {eta^2-(1-e^2), s^2-(1-c^2), c2g^2-(1-s2g^2)}\n\n');

syms A B;

% Build cos(jf), sin(jf), cos(jf+2g), sin(jf+2g) for j = 1, 2, 3.
cos_jf = cell(3, 1);  sin_jf = cell(3, 1);
cos_jfp2g = cell(3, 1);  sin_jfp2g = cell(3, 1);
for j = 1:3
  cos_jf{j} = chebT(j, cos_f);
  sin_jf{j} = sin_f * chebU(j - 1, cos_f);
  cos_jfp2g{j} = cos_jf{j} * c2g - sin_jf{j} * s2g;
  sin_jfp2g{j} = sin_jf{j} * c2g + cos_jf{j} * s2g;
end

% =================================================================
% C.5(b) Homological equation -- B-piece (the nontrivial check).
% =================================================================
% d/dl (sin(j f + 2 g)) = j (df/dl) cos(j f + 2 g) = j (eta/kappa^2) cos(j f + 2 g)
% Sum:
%   d/dl [3 sin(2f+2g) + 3e sin(f+2g) + e sin(3f+2g)]
% should equal:
%   6 (eta^3/kappa^3) cos(2(f+g)).
%
% (We verify the coefficient by 6 because the (B/6) prefactor cancels.)

df_dl = eta / kappa^2;
dXB_brackets = 3 * 2 * cos_jfp2g{2} * df_dl ...
             + 3 * e * 1 * cos_jfp2g{1} * df_dl ...
             + e * 3 * cos_jfp2g{3} * df_dl;
% This is the d/dl of [3 sin(2f+2g) + 3e sin(f+2g) + e sin(3f+2g)].

target_brackets = 6 * (eta^3 / kappa^3) * cos_jfp2g{2};

diff_brackets = dXB_brackets - target_brackets;

fprintf('--- C.5(b) Homological equation (B-piece) ---\n');
fprintf('  Verifying:  d/dl [3 sin(2f+2g) + 3e sin(f+2g) + e sin(3f+2g)]\n');
fprintf('           ==  6 (eta^3 / kappa^3) cos(2(f+g))\n');
fprintf('  Reducing difference symbolically...\n');
if vanishes_full(diff_brackets, e, eta, c, s, c2g, s2g)
  fprintf('  *** PASS: B-piece of homological equation verified. ***\n\n');
else
  fprintf('  *** FAIL: B-piece does not vanish. ***\n\n');
end

% =================================================================
% C.5(b) Homological equation -- A-piece.
% =================================================================
% d/dl (f - l + e sin f) = (df/dl)(1 + e cos f) - 1
%                       = (eta/kappa^2)(1 + e cos f) - 1
% Target: (a/r)^3 eta^3 - 1 = eta^3/kappa^3 - 1.
% Need:  (eta/kappa^2)(1 + e cos f) = eta^3/kappa^3
%   <=>  (1 + e cos f) = eta^2/kappa
% This is the orbit equation in true anomaly.

ddl_A_inner_minus_target = (df_dl * (1 + e * cos_f) - 1) - (eta^3 / kappa^3 - 1);

fprintf('--- C.5(b) Homological equation (A-piece) ---\n');
fprintf('  Verifying:  d/dl (f - l + e sin f)  ==  (a/r)^3 eta^3 - 1\n');
fprintf('  i.e., (df/dl)(1 + e cos f) == eta^3 / kappa^3\n');
if vanishes_full(ddl_A_inner_minus_target, e, eta, c, s, c2g, s2g)
  fprintf('  *** PASS: A-piece of homological equation verified. ***\n\n');
else
  fprintf('  *** FAIL: A-piece does not vanish. ***\n\n');
end

% =================================================================
% C.5(a) Zero mean of S_1.
% =================================================================
% By construction S_1 is the l-mean-subtracted form.  We verify by
% reassembling Theorem C.4 from M_1, M_2, M_3 and checking <S_1>_l = 0
% directly.  Use the kappa-reduction approach to compute M_1, M_2, M_3
% as in independent_X0_octave.m.

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

function Mj = compute_Mj(j, e, eta, c, s)
  syms kap;
  a_z = c - e;
  b_z = eta * s;
  Re_zj = sym(0);
  for k = 0:j
    if mod(k, 4) == 0
      sgn = 1;
    elseif mod(k, 4) == 2
      sgn = -1;
    else
      continue;
    end
    Re_zj = Re_zj + sgn * nchoosek(j, k) * a_z^(j - k) * b_z^k;
  end
  Re_zj = expand(Re_zj);
  for it = 1:8
    Re_zj = expand(subs(Re_zj, s^2, 1 - c^2));
  end
  P_in_kap = expand(subs(Re_zj, c, (1 - kap)/e));
  Mj_int = sym(0);
  d_max = 0;
  current = P_in_kap;
  for trial = 0:10
    if ~isequal(simplify(subs(current, kap, 0)), sym(0))
      d_max = trial;
    end
    current = diff(current, kap);
  end
  current = P_in_kap;
  for q = 0:d_max
    bq = simplify(subs(current, kap, 0) / factorial(q));
    Mj_int = Mj_int + bq * I_p((j - 1) - q, e, eta);
    current = diff(current, kap);
  end
  Mj = simplify(Mj_int / (2 * sym(pi)));
end

M1 = compute_Mj(1, e, eta, c, s);
M2 = compute_Mj(2, e, eta, c, s);
M3 = compute_Mj(3, e, eta, c, s);

fprintf('--- C.5(a) Zero mean ---\n');
fprintf('  Recomputed orbital averages (independent kappa-reduction):\n');
fprintf('    M_1 = %s\n', char(M1));
fprintf('    M_2 = %s\n', char(M2));
fprintf('    M_3 = %s\n', char(M3));

% In Theorem C.4, the sin(2g) coefficient (after mean subtraction) is
%   B M_2 / 6  *  (1/eta^3)
% (in the suppressed-prefactor form, the 1/eta^3 is also suppressed).
% Using the form X(E, g) of S_1 modulo prefactor mu*k_2/(n*a^3*eta^3):
%
%   X = A (f - l + e sin f)
%     + (B/2) sin(2f+2g) + (eB/2) sin(f+2g) + (eB/6) sin(3f+2g)
%     + (B M_2 / 6) sin(2g)
%
% l-mean (using <f - l + e sin f>_l = 0 from ch07b C.4.3):
%   <X>_l = (B/2) s2g M_2  +  (eB/2) s2g M_1  +  (eB/6) s2g M_3  +  (B M_2 / 6) s2g
%        = s2g * (B/6) * [3 M_2 + 3 e M_1 + e M_3 + M_2]
% By Lemma C.4.10 (independently verified in independent_X0_octave.m):
%        3 M_2 + 3 e M_1 + e M_3 = -M_2
% So:    <X>_l = s2g * (B/6) * [-M_2 + M_2] = 0.

mean_X_coef_s2g = (B/2) * M2 + (e*B/2) * M1 + (e*B/6) * M3 + (B * M2 / 6);
fprintf('  Coefficient of sin(2g) in <X>_l (must reduce to 0):\n');
fprintf('    raw = %s\n', char(mean_X_coef_s2g));
if vanishes_full(mean_X_coef_s2g, e, eta, c, s, c2g, s2g)
  fprintf('  *** PASS: <S_1>_l = 0 verified symbolically. ***\n\n');
else
  fprintf('  *** FAIL: <S_1>_l does not vanish. ***\n\n');
end

% =================================================================
% Final cross-checks: M_2 and M_3 against the framework / plan values.
% =================================================================
candidate_M2 = (3*e^2 - 2 + 2*eta^3) / e^2;
fprintf('--- Cross-check: M_2 vs Cor B.0.7-7 ((3e^2 - 2 + 2 eta^3)/e^2) ---\n');
if vanishes_full(M2 - candidate_M2, e, eta, c, s, c2g, s2g)
  fprintf('  *** PASS. ***\n\n');
else
  fprintf('  *** FAIL. ***\n\n');
end

candidate_M3 = (8 - 12*e^2 + 3*e^4 - 8*eta^3) / e^3;
fprintf('--- Cross-check: M_3 vs recurrence-derived ((8 - 12e^2 + 3e^4 - 8 eta^3)/e^3) ---\n');
if vanishes_full(M3 - candidate_M3, e, eta, c, s, c2g, s2g)
  fprintf('  *** PASS. ***\n\n');
else
  fprintf('  *** FAIL. ***\n\n');
end

fprintf('\n=== verify_ch07c_S1_symbolic.m  COMPLETE. ===\n');
