% verify_ch10_thm1_momentum_partials.m
%
% PURE SYMBOLIC verification of Chapter 10 Theorem 1
% (Momentum-partial closure -- Extended GSI).
%
% Theorem (T1).  For any S = G^(-alpha) * F(theta, e, l, g) with
%
%     theta = H / G,
%     e     = sqrt(1 - G^2/L^2),
%     eta   = sqrt(1 - e^2) = G/L,      (Delaunay identity)
%
% the three Delaunay momentum partials admit closed forms
%
%   (T1.L)  dS/dL|_{G,H,l,g} = G^(-(alpha+1)) * (eta^3 / e) * dF/de,
%
%   (T1.G)  dS/dG|_{L,H,l,g} = G^(-(alpha+1)) * [ -alpha F
%                                                 - theta * dF/dtheta
%                                                 - (eta^2 / e) * dF/de ],
%
%   (T1.H)  dS/dH|_{L,G,l,g} = G^(-(alpha+1)) * dF/dtheta.
%
% This verifier instantiates F as an ABSTRACT symbolic function of
% (theta, e, l, g) using Octave's sym('F(...)') idiom (the Octave
% analogue of SymPy's sp.Function('F')).  This is the same pattern
% as verify_GSI_symbolic.py (which does this in pure Python).
%
% For each claim (T1.L), (T1.G), (T1.H), the verifier:
%   (i)  computes the LHS via Octave `diff()` applied to the fully
%        symbolic S = G^(-alpha) * F(theta(G,H), e(L,G), l, g);
%   (ii) constructs the RHS using the extracted abstract partials
%        F_theta_at := G * dS/dH * G^(alpha),    (from T1.H inversion)
%        F_e_at     := dS/dL * G^(alpha) * (e/eta^3);
%        but critically, F_theta_at and F_e_at are ALSO obtainable
%        directly from diff(F, H)*G and diff(F, L)/(de/dL), i.e.,
%        from independent chain-rule extractions.  We use the latter
%        independent extractions so the test is NOT circular.
%  (iii) simplifies (LHS - RHS) and checks it reduces to 0 after
%        substituting eta = G/L (via eta^2 = 1 - e^2 = G^2/L^2, G,L > 0).
%
% PASS means: the three differences are identically 0 in SymPy's
% symbolic ring.  Since F is abstract, PASS establishes Theorem 1 for
% ANY smooth F and any real alpha.
%
% Reference: ch10_foundations_thm1.md (Section 3 Proof, Section 5
%            Symbolic verification note).
% Analogue: design/derivations/verify_GSI_symbolic.py (Python, for
%           ch09b Theorem E.3 / original GSI).

pkg load symbolic;
clear; close all;
addpath(fileparts(mfilename('fullpath')));

printf('============================================================\n');
printf('Ch 10 Theorem 1 -- Momentum-partial closure (Extended GSI)\n');
printf('Symbolic verification of (T1.L), (T1.G), (T1.H)\n');
printf('for abstract F(theta, e, l, g) and real alpha.\n');
printf('============================================================\n\n');

% Declare real symbols.  L, G, H positive so that eta = G/L > 0 and
% e = sqrt(1 - G^2/L^2) is well-defined in the elliptical regime.
syms L G H real positive;
syms l g real;
syms alpha real;

% Delaunay-derived geometric variables.
theta_expr = H / G;
e_expr     = sqrt(1 - G^2 / L^2);
eta_expr   = sqrt(1 - e_expr^2);     % SymPy resolves eta = G/L via simplify

% -------------------------------------------------------------------
% Abstract F(theta, e, l, g).
%
% Octave's symbolic package wraps SymPy.  The call
%     sym(str_with_parens)
% where str_with_parens names an un-defined function parses as
% SymPy's abstract Function application, which supports chain-rule
% differentiation via Octave `diff()`.  This is the recipe used
% by verify_GSI_symbolic.py (in pure Python: sp.Function('F')(a,b,c,d)).
% -------------------------------------------------------------------

% Abstract F via symbolic function declaration.
syms F(arg1, arg2, arg3, arg4);
F_applied = F(theta_expr, e_expr, l, g);

% Build S = G^(-alpha) * F(theta, e, l, g).
S = G^(-alpha) * F_applied;

printf('S := G^(-alpha) * F(theta, e, l, g)\n');
printf('    theta = H/G\n');
printf('    e     = sqrt(1 - G^2/L^2)\n');
printf('    eta   = sqrt(1 - e^2)  (= G/L, Delaunay identity)\n\n');

% Compute the three Delaunay momentum partials via Octave diff().
% SymPy chain-rule expands theta(G,H), e(L,G), and the abstract F
% automatically.
dS_dL = diff(S, L);
dS_dG = diff(S, G);
dS_dH = diff(S, H);

% -------------------------------------------------------------------
% Extract abstract partials F_theta_at, F_e_at from independent
% chain-rule calls on F_applied directly (NOT via the Theorem 1
% predictions, so the test remains non-circular).
%
% Reasoning:
%   diff(F_applied, H) = d/dH [F(theta(H,G), e(L,G), l, g)]
%                      = dF/dtheta * d(theta)/dH
%                      = dF/dtheta * (1/G)
%   =>  dF/dtheta = G * diff(F_applied, H).
%
%   diff(F_applied, L) = dF/de * de/dL
%   =>  dF/de = diff(F_applied, L) / (de/dL).
% -------------------------------------------------------------------
F_theta_at = G * diff(F_applied, H);
de_dL      = diff(e_expr, L);
F_e_at     = diff(F_applied, L) / de_dL;

% -------------------------------------------------------------------
% CHECK 1 : (T1.H)  dS/dH = G^(-(alpha+1)) * dF/dtheta.
%
% LHS = dS/dH.
% RHS = G^(-(alpha+1)) * F_theta_at.
% Expected: LHS - RHS = 0 exactly (no sqrt reduction needed because
% dS/dH and F_theta_at both contain exactly the same Derivative
% object times prefactors that cancel exactly).
% -------------------------------------------------------------------
printf('--- Check 1: (T1.H)  dS/dH = G^(-(alpha+1)) * dF/dtheta ---\n\n');

RHS_H  = G^(-(alpha+1)) * F_theta_at;
diff_H = simplify(dS_dH - RHS_H);

printf('  dS/dH - RHS_H = %s\n', char(diff_H));
pass_H = isequal(diff_H, sym(0));
if pass_H
  printf('  [PASS] (T1.H) verified: dS/dH = G^(-(alpha+1)) * dF/dtheta\n\n');
else
  printf('  [FAIL] (T1.H): dS/dH - RHS_H did not reduce to 0.\n\n');
end

% -------------------------------------------------------------------
% CHECK 2 : (T1.L)  dS/dL = G^(-(alpha+1)) * (eta^3/e) * dF/de.
%
% LHS = dS/dL = G^(-alpha) * dF/de * de/dL = G^(-alpha) * (dF/de) * eta^2/(L*e).
% RHS = G^(-(alpha+1)) * (eta^3/e) * dF/de.
% Difference (factoring dF/de out):
%    G^(-alpha) * eta^2/(L*e) - G^(-(alpha+1)) * eta^3/e
%  = G^(-alpha)/e * [ eta^2/L - eta^3/G ]
% Under the Delaunay identity 1/L = eta/G (equivalently eta = G/L),
% eta^2/L = eta^3/G, so the difference is 0.
% To force SymPy to apply this reduction, we substitute sqrt(1-e^2) -> G/L
% or equivalently eta_expr = sqrt(1 - (1 - G^2/L^2)) -> G/L.
% -------------------------------------------------------------------
printf('--- Check 2: (T1.L)  dS/dL = G^(-(alpha+1)) * (eta^3/e) * dF/de ---\n\n');

RHS_L      = G^(-(alpha+1)) * (eta_expr^3 / e_expr) * F_e_at;
diff_L_raw = dS_dL - RHS_L;

% Force the eta = G/L reduction.  SymPy may carry sqrt(1 - (1 - G^2/L^2))
% unsimplified, or it may pre-normalize to sqrt(G^2/L^2) -- try both
% patterns.  (With G, L declared positive, sqrt(G^2/L^2) -> G/L
% automatically, but we make it explicit for robustness.)
diff_L_sub = diff_L_raw;
try
  diff_L_sub = subs(diff_L_sub, sqrt(1 - (1 - G^2/L^2)), G/L);
catch
end
try
  diff_L_sub = subs(diff_L_sub, sqrt(G^2/L^2), G/L);
catch
end
diff_L_sub = simplify(diff_L_sub);

printf('  (dS/dL - RHS_L) after eta = G/L substitution = %s\n', char(diff_L_sub));
pass_L = isequal(diff_L_sub, sym(0));
if pass_L
  printf('  [PASS] (T1.L) verified: dS/dL = G^(-(alpha+1)) * (eta^3/e) * dF/de\n\n');
else
  printf('  [FAIL] (T1.L): difference did not reduce to 0.\n\n');
end

% -------------------------------------------------------------------
% CHECK 3 : (T1.G)  dS/dG = G^(-(alpha+1)) * [-alpha F - theta dF/dtheta - (eta^2/e) dF/de].
%
% LHS = dS/dG expanded via chain rule through the explicit G^(-alpha)
%       prefactor, theta(G,H), and e(L,G).
% RHS = three-term closed form.
% The difference must reduce to 0 after eta = G/L substitution.
% -------------------------------------------------------------------
printf('--- Check 3: (T1.G)  dS/dG = G^(-(alpha+1)) * [-alpha F - theta dF/dtheta - (eta^2/e) dF/de] ---\n\n');

term_a = -alpha * F_applied * G^(-(alpha+1));
term_b = -theta_expr * F_theta_at * G^(-(alpha+1));
term_c = -(eta_expr^2 / e_expr) * F_e_at * G^(-(alpha+1));
RHS_G  = term_a + term_b + term_c;

diff_G_raw = dS_dG - RHS_G;
diff_G_sub = diff_G_raw;
try
  diff_G_sub = subs(diff_G_sub, sqrt(1 - (1 - G^2/L^2)), G/L);
catch
end
try
  diff_G_sub = subs(diff_G_sub, sqrt(G^2/L^2), G/L);
catch
end
diff_G_sub = simplify(diff_G_sub);

printf('  (dS/dG - RHS_G) after eta = G/L substitution = %s\n', char(diff_G_sub));
pass_G = isequal(diff_G_sub, sym(0));
if pass_G
  printf('  [PASS] (T1.G) verified: dS/dG closed form\n\n');
else
  printf('  [FAIL] (T1.G): difference did not reduce to 0.\n\n');
end

% -------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------
printf('============================================================\n');
printf('Ch 10 Theorem 1 verification summary\n');
printf('============================================================\n');
if pass_L, s_L = 'PASS'; else s_L = 'FAIL'; end
if pass_G, s_G = 'PASS'; else s_G = 'FAIL'; end
if pass_H, s_H = 'PASS'; else s_H = 'FAIL'; end
printf('  (T1.L)  dS/dL closed form:  %s\n', s_L);
printf('  (T1.G)  dS/dG closed form:  %s\n', s_G);
printf('  (T1.H)  dS/dH closed form:  %s\n', s_H);
printf('\n');
if pass_L && pass_G && pass_H
  printf('*** ALL THREE PASS: Theorem 1 verified for abstract F and any alpha. ***\n');
  printf('\n');
  printf('This establishes that for ANY smooth F(theta, e, l, g) and any real alpha,\n');
  printf('the three momentum partials of S = G^(-alpha) F(theta, e, l, g) lie in\n');
  printf('M_{alpha+1} with the explicit closed forms (T1.L), (T1.G), (T1.H).\n');
else
  printf('*** ONE OR MORE CHECKS FAILED.  Theorem 1 not confirmed. ***\n');
end
printf('============================================================\n');
