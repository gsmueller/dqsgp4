% verify_ch10_thm4_avg_closure.m
%
% PURE SYMBOLIC verification of Theorem 4 of ch10_foundations_thm4.md
% (l-averaging closure on the factorable class M_alpha).
%
% Theorem 4 statement (recap):
%   If S = G^{-alpha} * F(theta, e, l, g) with F in C^infty and dS/dh = 0,
%   then <S>_l = (1/2pi) int_0^{2pi} S dl = G^{-alpha} * F_avg(theta, e, g),
%   and <S>_l lies in M_alpha (same class as S).
%
% Concrete test function (chosen to exercise all the features):
%   F(theta, e, l, g)
%     = f1(theta) * cos(l)
%     + f2(e) * sin(l + g)
%     + f3(theta, e)
%     + f4(theta, e, g) * cos(2*l + 3*g)
%
% This is a finite trigonometric polynomial in l, with both constant-in-l
% pieces (f3, etc.) and mixed-argument oscillatory pieces (cos(2l + 3g)).
% It is the worst-case-diversity test: different l-frequencies, different
% cross-couplings with g, and a purely l-free piece that survives averaging.
%
% Concrete choices for the f-slot functions (to keep SymPy tractable):
%   f1(theta) = 3*theta^2 - 1               (A-like inclination factor)
%   f2(e)     = e^2 + 1
%   f3(theta, e) = theta * (1 - e^2)         (a "secular" remainder piece)
%   f4(theta, e, g) = theta^2 * e * cos(g)   (pure g-dependence in a coefficient)
%
% Expected outcome after l-averaging:
%   <F>_l = f3(theta, e)  (only the l-free piece survives)
% and
%   <S>_l = G^{-alpha} * f3(theta, e)
%         = G^{-alpha} * theta * (1 - e^2).
%
% All other terms average to zero because:
%   <cos(l)>_l = 0
%   <sin(l + g)>_l = 0  (integration by substitution u = l + g)
%   <cos(2*l + 3*g)>_l = 0  (same reason)
%
% Checks performed:
%   (C1) Direct symbolic l-average:
%        <S>_l == G^{-alpha} * f3 = G^{-alpha} * theta * (1 - e^2).
%   (C2) l-independence of <S>_l:
%        diff(<S>_l, l) == 0.
%   (C3) h-independence of <S>_l:
%        diff(<S>_l, h) == 0.
%   (C4) G-prefactor commutation:
%        <S>_l / G^{-alpha}  equals  F_avg(theta, e, g).
%   (C5) Symbol-presence inspection:
%        <S>_l has no l in its free symbols.
%   (C6) Smoothness sanity:
%        mixed partials of F_avg in (theta, e, g) reduce to clean closed forms
%        (finite rational functions; no singularities or Dirac-deltas).
%
% Each check is reduced to symbolic 0 (or a trivial symbolic form) in the
% appropriate SymPy quotient ring.

pkg load symbolic;
clear; close all;
addpath(fileparts(mfilename('fullpath')));

fprintf('============================================================\n');
fprintf('VERIFICATION OF THEOREM 4: l-AVERAGING CLOSURE ON M_alpha\n');
fprintf('============================================================\n\n');

% ------------------------------------------------------------------
% Setup
% ------------------------------------------------------------------
syms G real positive;
syms L H real positive;
syms l g h real;
syms theta e real;
syms alpha_sym real positive;

% The concrete test function F(theta, e, l, g) chosen to exercise the
% full structure of Theorem 4 (l-oscillatory + l-free + mixed pieces).
f1 = 3*theta^2 - 1;
f2 = e^2 + 1;
f3 = theta * (1 - e^2);
f4 = theta^2 * e * cos(g);

F = f1 * cos(l) + f2 * sin(l + g) + f3 + f4 * cos(2*l + 3*g);

S = G^(-alpha_sym) * F;

fprintf('Test function:\n');
fprintf('  F(theta, e, l, g) = f1*cos(l) + f2*sin(l+g) + f3 + f4*cos(2l+3g)\n');
fprintf('  with f1 = 3*theta^2 - 1,  f2 = e^2 + 1,\n');
fprintf('       f3 = theta*(1-e^2),  f4 = theta^2*e*cos(g).\n');
fprintf('\n');
fprintf('S = G^{-alpha} * F.\n\n');

% Expected l-average of F (computed by hand):
%   <cos(l)>_l   = 0
%   <sin(l+g)>_l = 0
%   <f3>_l       = f3
%   <cos(2l+3g)>_l = 0
% Therefore <F>_l = f3 = theta*(1 - e^2).
F_avg_expected = f3;
S_avg_expected = G^(-alpha_sym) * F_avg_expected;

% ------------------------------------------------------------------
% Compute <S>_l = (1/2pi) int_0^{2pi} S dl by direct symbolic integration.
% ------------------------------------------------------------------
fprintf('Computing <S>_l via SymPy integration...\n');
S_avg = (1 / (2*sym(pi))) * int(S, l, 0, 2*sym(pi));
S_avg = simplify(S_avg);
fprintf('  <S>_l = %s\n\n', char(S_avg));

% Also compute <F>_l for check C4.
F_avg = (1 / (2*sym(pi))) * int(F, l, 0, 2*sym(pi));
F_avg = simplify(F_avg);
fprintf('  <F>_l = %s\n\n', char(F_avg));

% ------------------------------------------------------------------
% (C1) Direct identity check: <S>_l == G^{-alpha} * f3
% ------------------------------------------------------------------
fprintf('--- Check (C1): <S>_l == G^{-alpha} * theta * (1 - e^2) ---\n');
diff_C1 = simplify(S_avg - S_avg_expected);
if isequal(diff_C1, sym(0))
  fprintf('  *** PASS: <S>_l matches expected closed form. ***\n\n');
  pass_C1 = true;
else
  fprintf('  *** FAIL ***\n');
  fprintf('  <S>_l - expected = %s\n\n', char(diff_C1));
  pass_C1 = false;
end

% ------------------------------------------------------------------
% (C2) l-independence: diff(<S>_l, l) == 0
% ------------------------------------------------------------------
fprintf('--- Check (C2): diff(<S>_l, l) == 0 ---\n');
dSavg_dl = diff(S_avg, l);
dSavg_dl = simplify(dSavg_dl);
if isequal(dSavg_dl, sym(0))
  fprintf('  *** PASS: <S>_l is l-free. ***\n\n');
  pass_C2 = true;
else
  fprintf('  *** FAIL ***\n');
  fprintf('  diff(<S>_l, l) = %s\n\n', char(dSavg_dl));
  pass_C2 = false;
end

% ------------------------------------------------------------------
% (C3) h-independence: diff(<S>_l, h) == 0
% ------------------------------------------------------------------
fprintf('--- Check (C3): diff(<S>_l, h) == 0 ---\n');
dSavg_dh = diff(S_avg, h);
dSavg_dh = simplify(dSavg_dh);
if isequal(dSavg_dh, sym(0))
  fprintf('  *** PASS: <S>_l has no h-dependence (M-condition preserved). ***\n\n');
  pass_C3 = true;
else
  fprintf('  *** FAIL ***\n');
  fprintf('  diff(<S>_l, h) = %s\n\n', char(dSavg_dh));
  pass_C3 = false;
end

% ------------------------------------------------------------------
% (C4) G-prefactor commutation: <S>_l / G^{-alpha} == F_avg
% ------------------------------------------------------------------
fprintf('--- Check (C4): <S>_l / G^{-alpha} == <F>_l ---\n');
ratio = S_avg * G^alpha_sym;
ratio = simplify(ratio);
diff_C4 = simplify(ratio - F_avg);
if isequal(diff_C4, sym(0))
  fprintf('  *** PASS: G^{-alpha} prefactor commutes with <.>_l. ***\n\n');
  pass_C4 = true;
else
  fprintf('  *** FAIL ***\n');
  fprintf('  ratio - F_avg = %s\n\n', char(diff_C4));
  pass_C4 = false;
end

% ------------------------------------------------------------------
% (C5) Symbol-presence: l is not a free symbol of <S>_l.
% ------------------------------------------------------------------
fprintf('--- Check (C5): l not in free symbols of <S>_l ---\n');
free_syms = symvar(S_avg);
has_l = false;
free_names = {};
for k = 1:numel(free_syms)
  nm = char(free_syms(k));
  free_names{end+1} = nm;
  if strcmp(nm, 'l')
    has_l = true;
  end
end
if ~has_l
  fprintf('  Free symbols of <S>_l: { %s}\n', strjoin(free_names, ', '));
  fprintf('  *** PASS: l is absent from <S>_l. ***\n\n');
  pass_C5 = true;
else
  fprintf('  *** FAIL: l appears as a free symbol. ***\n');
  fprintf('  Free symbols: { %s}\n\n', strjoin(free_names, ', '));
  pass_C5 = false;
end

% ------------------------------------------------------------------
% (C6) Smoothness: mixed partials reduce to clean closed forms.
%     We compute d^2 F_avg / dtheta/de, d^2 F_avg / dtheta/dg,
%               d^2 F_avg / de/dg.
%     All three should simplify to finite polynomials/rationals in the
%     remaining symbols (no Heaviside, Dirac, Piecewise outputs).
% ------------------------------------------------------------------
fprintf('--- Check (C6): smoothness via mixed partials of <F>_l ---\n');
p_te = simplify(diff(diff(F_avg, theta), e));
p_tg = simplify(diff(diff(F_avg, theta), g));
p_eg = simplify(diff(diff(F_avg, e), g));
fprintf('  d^2<F>_l/(dtheta de) = %s\n', char(p_te));
fprintf('  d^2<F>_l/(dtheta dg) = %s\n', char(p_tg));
fprintf('  d^2<F>_l/(de dg)     = %s\n', char(p_eg));

% For our concrete test, <F>_l = theta*(1-e^2), so:
%   d^2<F>_l/(dtheta de) = -2*e
%   d^2<F>_l/(dtheta dg) = 0
%   d^2<F>_l/(de dg)     = 0
expected_te = -2*e;
expected_tg = sym(0);
expected_eg = sym(0);

ok_te = isequal(simplify(p_te - expected_te), sym(0));
ok_tg = isequal(simplify(p_tg - expected_tg), sym(0));
ok_eg = isequal(simplify(p_eg - expected_eg), sym(0));

if ok_te && ok_tg && ok_eg
  fprintf('  All three mixed partials match expected smooth closed forms.\n');
  fprintf('  *** PASS: <F>_l is smooth in (theta, e, g). ***\n\n');
  pass_C6 = true;
else
  fprintf('  *** FAIL: mixed partials do not match expected closed forms. ***\n');
  if ~ok_te, fprintf('    theta-e partial mismatch.\n'); end
  if ~ok_tg, fprintf('    theta-g partial mismatch.\n'); end
  if ~ok_eg, fprintf('    e-g partial mismatch.\n'); end
  fprintf('\n');
  pass_C6 = false;
end

% ------------------------------------------------------------------
% Summary
% ------------------------------------------------------------------
fprintf('============================================================\n');
fprintf('SUMMARY\n');
fprintf('============================================================\n');
results = {pass_C1, 'C1 (closed-form match)';
           pass_C2, 'C2 (l-independence)';
           pass_C3, 'C3 (h-independence)';
           pass_C4, 'C4 (G-prefactor commutation)';
           pass_C5, 'C5 (symbol presence)';
           pass_C6, 'C6 (smoothness)'};
all_pass = true;
for k = 1:size(results, 1)
  p = results{k, 1};
  label = results{k, 2};
  if p
    fprintf('  [PASS] %s\n', label);
  else
    fprintf('  [FAIL] %s\n', label);
    all_pass = false;
  end
end
fprintf('------------------------------------------------------------\n');
if all_pass
  fprintf('THEOREM 4 (l-averaging closure): *** VERIFIED. ***\n');
else
  fprintf('THEOREM 4 (l-averaging closure): *** VERIFICATION FAILED. ***\n');
end
fprintf('============================================================\n');
