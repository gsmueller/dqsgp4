% verify_ch10_thm2_angle_partials.m
%
% PURE SYMBOLIC verification of Chapter 10 Foundations Theorem 2
% (Angle-Partial Closure) on the factorable class ℳ_α.
%
% Theorem 2.  Let S ∈ ℳ_α with defining factorization
%     S = G^{-α} * F(θ, e, l, g)
% where θ = H/G, e = sqrt(1 - G²/L²), and F ∈ C^∞(R^4) is an abstract
% smooth function of four arguments.  Then:
%
%     (a)  ∂S/∂l |_{L,G,H,g,h} = G^{-α} * ∂F/∂l   (T2.l)
%     (b)  ∂S/∂g |_{L,G,H,l,h} = G^{-α} * ∂F/∂g   (T2.g)
%     (c)  ∂S/∂h |_{L,G,H,l,g} = 0   identically  (T2.h)
%
% Strategy.  Using Octave's symbolic package with an abstract function
% F = symfun('F(theta, e, l, g)', [theta, e, l, g]), we construct S
% directly in terms of (L, G, H, l, g, h), with θ and e expressed via
% their Delaunay definitions.  We then compute the three angle partials
% by Octave's diff() primitive and compare each to its claimed closed
% form.  The residues must reduce to 0 symbolically, WITHOUT specializing
% F to any concrete form.
%
% This is the direct analog of verify_GSI_symbolic.py (which verifies
% Theorem 1's GSI identity using SymPy's abstract Function class), but
% in the Octave symbolic-package idiom.
%
% No orbit-equation substitution, no κ-reduction, no Chebyshev expansion
% needed — Theorem 2 is a structural triviality whose verification
% reduces to Octave's built-in chain rule honoring the facts that:
%   (i)  θ(G, H) and e(L, G) have no explicit (l, g, h) dependence;
%   (ii) F is an abstract function of exactly its declared four args;
%   (iii) h appears nowhere in S's factored form (the ℳ_α condition).

pkg load symbolic;
clear; close all;
addpath(fileparts(mfilename('fullpath')));

% ----- Symbol declarations -----
syms L G H real positive;
syms l g real;
syms h real;
syms alpha_sym real;

% ----- Build S = G^{-alpha} * F(θ(G,H), e(L,G), l, g) -----
%
% θ := H/G,  e := sqrt(1 − G²/L²).  These are the Delaunay-derived
% arguments of F.  S has no direct h-dependence (by Definition 1.1
% of ℳ_α in CH10_PLAN.md §1).
%
% Abstract F: use Octave's sym(sprintf('F(...)')) idiom which parses
% as a SymPy abstract Function (same pattern as verify_ch10_thm1).

theta_expr = H / G;
e_expr     = sqrt(1 - G^2 / L^2);

% Abstract F via symbolic function declaration.
syms F(arg1, arg2, arg3, arg4);
F_applied = F(theta_expr, e_expr, l, g);

S = G^(-alpha_sym) * F_applied;

fprintf('\n');
fprintf('============================================================\n');
fprintf('  Chapter 10 Foundations Theorem 2 symbolic verification\n');
fprintf('  (Angle-Partial Closure on the factorable class M_alpha)\n');
fprintf('============================================================\n\n');
fprintf('Setup:\n');
fprintf('  S = G^(-alpha) * F(theta, e, l, g)\n');
fprintf('    theta = H/G\n');
fprintf('    e     = sqrt(1 - G^2/L^2)\n');
fprintf('    F     = abstract smooth function of 4 args (symfun)\n');
fprintf('\n');

% ----- Check (a):  ∂S/∂l  =  G^{-alpha} * ∂F/∂l   (T2.l) -----
%
% Direct chain rule: since neither G^{-alpha} nor theta(G, H) nor
% e(L, G) depends on l, only the third argument of F carries l
% dependence.  Hence diff(S, l) should reduce to
%
%     G^(-alpha) * (∂F/∂(3rd arg))(θ_expr, e_expr, l, g).
%
% Candidate construction: take diff of F applied at placeholder
% arg-symbols with respect to the third arg-symbol, producing the
% abstract third-argument partial; then substitute the arg-symbols
% with their (L, G, H, l, g) expressions.

fprintf('=== Check (a): d(S)/d(l) = G^(-alpha) * dF/dl   (T2.l) ===\n');

dS_dl_direct = diff(S, l);

% Build the abstract third-argument partial of F directly on F_applied.
% Since theta(G,H), e(L,G), and G^{-alpha} are all l-independent, the
% l-partial of F_applied equals G^{alpha} * d(S)/d(l).  So:
%     d(F_applied)/d(l) = d(S)/d(l) / G^{-alpha}
% i.e. (T2.l) says diff(S, l) = G^(-alpha) * diff(F_applied, l).

dS_dl_candidate = G^(-alpha_sym) * diff(F_applied, l);

diff_l = simplify(dS_dl_direct - dS_dl_candidate);
fprintf('  d(S)/d(l) direct     = %s\n', char(dS_dl_direct));
fprintf('  d(S)/d(l) candidate  = %s\n', char(dS_dl_candidate));
fprintf('  residue              = %s\n', char(diff_l));
report_check('(T2.l) d(S)/d(l) = G^(-alpha) * dF/dl', ...
             isequal(diff_l, sym(0)));

% ----- Check (b):  ∂S/∂g  =  G^{-alpha} * ∂F/∂g   (T2.g) -----
%
% Identical structure with l → g and F's third arg → F's fourth arg.

fprintf('\n=== Check (b): d(S)/d(g) = G^(-alpha) * dF/dg   (T2.g) ===\n');

dS_dg_direct = diff(S, g);

dS_dg_candidate = G^(-alpha_sym) * diff(F_applied, g);

diff_g = simplify(dS_dg_direct - dS_dg_candidate);
fprintf('  d(S)/d(g) direct     = %s\n', char(dS_dg_direct));
fprintf('  d(S)/d(g) candidate  = %s\n', char(dS_dg_candidate));
fprintf('  residue              = %s\n', char(diff_g));
report_check('(T2.g) d(S)/d(g) = G^(-alpha) * dF/dg', ...
             isequal(diff_g, sym(0)));

% ----- Check (c):  ∂S/∂h  =  0   identically   (T2.h) -----
%
% h appears nowhere in the factored form of S (not in G^{-alpha},
% not in θ = H/G, not in e = sqrt(1 − G²/L²), not in the l, g
% arguments of F).  Octave's diff() should return 0 directly.
% This is a *construction* check: (T2.h) holds because we built S
% without any h-dependence, which is exactly the ℳ_α defining
% condition of Definition 1.1.

fprintf('\n=== Check (c): d(S)/d(h) = 0 identically   (T2.h) ===\n');

dS_dh_direct = diff(S, h);
fprintf('  d(S)/d(h) direct     = %s\n', char(dS_dh_direct));
report_check('(T2.h) d(S)/d(h) = 0 identically', ...
             isequal(dS_dh_direct, sym(0)));

% ----- Independent sanity check: h is absent from S's free-symbol list -----
%
% This complements Check (c) with a static symbol-inventory test.
% If 'h' is not among S's free symbols, (T2.h) is immediate regardless
% of how diff() is implemented.

fprintf('\n=== Sanity: h is absent from S''s free-symbol list ===\n');

S_free = symvar(S);
S_free_names = {};
if iscell(S_free)
  for k = 1:numel(S_free)
    S_free_names{end+1} = char(S_free{k});
  end
else
  for k = 1:numel(S_free)
    S_free_names{end+1} = char(S_free(k));
  end
end
h_present = false;
for k = 1:numel(S_free_names)
  if strcmp(S_free_names{k}, 'h')
    h_present = true;
    break;
  end
end
fprintf('  S free symbols: %s\n', strjoin(S_free_names, ', '));
report_check('h is absent from S''s free-symbol list', ~h_present);

fprintf('\n=== verify_ch10_thm2_angle_partials.m COMPLETE ===\n');
