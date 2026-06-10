% verify_ch09d_basis_completeness.m
%
% PURE SYMBOLIC verification of Proposition E.5 (per-harmonic
% decomposition of S_1 and its three momentum partials in an 8-element
% basis).
%
% Basis ordering (Definition E.5.1) -- S_1-originating elements first,
% chain-rule-induced elements second:
%
%   E_1 = f - l + e sin(f)                           (alpha; S_1, dS_1/dH)
%   E_2 = sin(2(f+g))                                (S_1, dS_1/dH)
%   E_3 = sin(f + 2g)                                (all four partials)
%   E_4 = sin(3f + 2g)                               (all four partials)
%   E_5 = sin(2g)                                    (all four partials)
%   E_6 = sin(f)/kappa                               (dS_1/dL only + GSI transport)
%   E_7 = sin^3(f)                                   (dS_1/dL only + GSI transport)
%   E_8 = sin(f)(2 + e cos f) cos(2(f+g))/kappa      (dS_1/dL only + GSI transport)
%
% Checks (each reduces to symbolic 0 in an appropriate quotient ring):
%
%   1.  S_1  --  {E_1, E_2, E_3, E_4, E_5}      (5 directions)
%         Theorem C.4 closed form equals the basis-table reconstruction.
%
%   2.  dS_1/dL  --  {E_3, E_4, E_5, E_6, E_7, E_8}   (6 directions)
%         DIRECT CHAIN RULE on Theorem C.4's closed form (computing
%         diff(inner, e) + diff(inner, eta)*(-e/eta) + diff(inner, f)*df/dL)
%         equals the basis-table reconstruction.  This catches errors in
%         the coefficient table AND in Prop E.2 independently, because
%         the direct chain rule does not consult Prop E.2.
%
%   3.  dS_1/dH  --  {E_1, E_2, E_3, E_4, E_5}  (5 directions)
%         DIRECT CHAIN RULE (through theta only) on Theorem C.4's closed
%         form equals the basis-table reconstruction.  H enters only
%         through theta = H/G, so the chain rule is a single path.
%
%   4.  dS_1/dG  --  {E_1, ..., E_8}            (8 directions)
%         DIRECT CHAIN RULE through all paths (prefactor, theta, e, f) on
%         Theorem C.4's closed form equals the GSI-inherited basis
%         expansion dS_1/dG = -(3/G) S_1 - theta dS_1/dH - (1/eta) dS_1/dL.
%         This is a concrete cross-check of the GSI applied to S_1
%         (redundant with verify_GSI_symbolic.py but worthwhile as an
%         instantiation, and independently tests the Part-(d) coefficient
%         identity of Proposition E.5).
%
% Reductions used (in this order):
%   (R1) kappa -> eta^2 / (1 + e cos(f_sym))   [orbit equation; eliminates
%        the kappa symbol so basis form and chain-rule form live in the
%        same ring]
%   (R2) sin^2(f_sym) -> 1 - cos^2(f_sym)
%   (R3) eta^2 -> 1 - e^2
%
% Note on the eta-chain.  The basis form contains X_0^{0,2}(e) =
% (3 e^2 - 2 + 2 eta^3)/e^2, which is explicitly a function of both e
% and eta.  The direct chain-rule derivatives must account for the
% implicit eta-dependence through eta^2 = 1 - e^2 by including
% diff(inner, eta) * (-e/eta).  Omission produces spurious residues.
% (See the audit note in verify_ch09bc_partials_symbolic.m for the
% cancellation that makes omission safe for the abstract GSI test only;
% here Checks 2 and 4 require the eta-chain explicitly.)

pkg load symbolic;
clear; close all;
addpath(fileparts(mfilename('fullpath')));  % vanishes_*, report_check

syms L G H real positive;
syms l g real;
syms e real positive;
syms eta real positive;     % eta = sqrt(1 - e^2); reduced iteratively
syms f_sym real;
syms kappa real positive;   % kappa = 1 - e cos E; eliminated via orbit equation
syms mu k_2 real positive;

theta = H / G;
A_expr = (3 * theta^2 - 1) / 2;
B_expr = 3 * (1 - theta^2) / 2;
X02 = (3*e^2 - 2 + 2*eta^3) / e^2;

% ---- 8 basis elements (Version A ordering: S_1-originating first) ----
E1 = f_sym - l + e*sin(f_sym);                                     % alpha
E2 = sin(2*(f_sym + g));                                           % sin(2(f+g))
E3 = sin(f_sym + 2*g);                                             % sin(f+2g)
E4 = sin(3*f_sym + 2*g);                                           % sin(3f+2g)
E5 = sin(2*g);                                                     % sin(2g)
E6 = sin(f_sym) / kappa;                                           % sin f / kappa
E7 = sin(f_sym)^3;                                                 % sin^3 f
E8 = sin(f_sym) * (2 + e*cos(f_sym)) * cos(2*(f_sym+g)) / kappa;   % sin f (2+e cos f) cos(2(f+g))/kappa

pref3 = mu^2 * k_2 / G^3;      % S_1 and dS_1/dL prefactor
pref4 = mu^2 * k_2 / G^4;      % dS_1/dH prefactor (one more 1/G from theta chain)

% ---- Inner bracket of Theorem C.4 ----
alpha_inner = f_sym - l + e*sin(f_sym);
beta_inner = 3*sin(2*(f_sym+g)) + 3*e*sin(f_sym+2*g) ...
           + e*sin(3*f_sym+2*g) + X02*sin(2*g);
inner = A_expr * alpha_inner + (B_expr/6) * beta_inner;

% Orbit-equation substitution for kappa (used to eliminate kappa_sym from
% basis forms so that Direct and Basis forms share the same ring).
kappa_sub = eta^2 / (1 + e*cos(f_sym));

%% =========================================================================
%% Check 1:  S_1 = basis expansion (5 directions)
%% =========================================================================

fprintf('\n=== Check 1: S_1 = basis expansion (5 directions) ===\n');

% Coefficients of S_1 in the basis (multiply by pref3):
cS1_1 = A_expr;                % E_1
cS1_2 = B_expr / 2;            % E_2   [(B/6)*3]
cS1_3 = e * B_expr / 2;        % E_3   [(B/6)*3e]
cS1_4 = e * B_expr / 6;        % E_4   [(B/6)*e]
cS1_5 = B_expr * X02 / 6;      % E_5   [(B/6)*X02]
% E_6, E_7, E_8 coefficients are 0 in S_1.

S1_basis = pref3 * (cS1_1*E1 + cS1_2*E2 + cS1_3*E3 + cS1_4*E4 + cS1_5*E5);

% Direct form from Theorem C.4 (ch07c)
S1_direct = pref3 * inner;

diff_1 = simplify(expand(S1_direct - S1_basis));
fprintf('  S_1_direct - S_1_basis = %s\n', char(diff_1));
report_check('Check 1: S_1 in 5-direction basis', isequal(diff_1, sym(0)));

%% =========================================================================
%% Check 2:  dS_1/dL = basis expansion (6 directions)
%%
%% Direct form: chain rule via (de/dL, df/dL) with eta-chain tied.
%% Basis form: table of ch09d Proposition E.5 column 2.
%% =========================================================================

fprintf('\n=== Check 2: dS_1/dL = basis expansion (6 directions) ===\n');

% Coefficients of dS_1/dL in the basis (multiply by pref3):
cL_3 = B_expr * eta^2 / (2 * L * e);                              % E_3
cL_4 = B_expr * eta^2 / (6 * L * e);                              % E_4
cL_5 = B_expr * eta^2 * (eta + 2) / (3 * L * (1 + eta)^2);        % E_5
cL_6 = 3 * A_expr * eta^2 / (L * e);                              % E_6
cL_7 = -A_expr * e / L;                                            % E_7
cL_8 = B_expr * eta^2 / (L * e);                                  % E_8
% E_1, E_2 coefficients are 0 in dS_1/dL.

dS1dL_basis = pref3 * (cL_3*E3 + cL_4*E4 + cL_5*E5 ...
                     + cL_6*E6 + cL_7*E7 + cL_8*E8);

% Direct form: apply the chain rule through (e, f) to Theorem C.4's
% closed form.  Prefactor is L-free by Lemma E.1.1.
%
% From ch08:  de/dL|_G = eta^2/(L e)     (Theorem D.2)
%             df/dL|_{l,G} = sin f (2 + e cos f)/(L e)    (Theorem D.8)
% The e-chain must also tie eta via deta/de = -e/eta, because X02
% contains eta^3 explicitly; otherwise a spurious residue appears.
de_dL = eta^2 / (L * e);
df_dL = sin(f_sym) * (2 + e*cos(f_sym)) / (L * e);

dinner_de = diff(inner, e) + diff(inner, eta) * (-e/eta);
dinner_df = diff(inner, f_sym);

dS1dL_direct = pref3 * (dinner_de * de_dL + dinner_df * df_dL);

% Eliminate kappa from the basis form via the orbit equation, so both
% forms live in the same ring Q[f_sym, e, eta, ...].
dS1dL_basis_sub = subs(dS1dL_basis, kappa, kappa_sub);

diff_L = dS1dL_direct - dS1dL_basis_sub;

% Iterate: reduce sin^2 f_sym -> 1 - cos^2 f_sym, eta^2 -> 1 - e^2.
diff_L_red = expand(diff_L);
for k_iter = 1:30
  prev = diff_L_red;
  diff_L_red = expand(subs(diff_L_red, sin(f_sym)^2, 1 - cos(f_sym)^2));
  diff_L_red = expand(subs(diff_L_red, eta^2, 1 - e^2));
  if isequal(diff_L_red, prev), break; end
end
diff_L_red = simplify(diff_L_red);

fprintf('  dS_1/dL_direct - dS_1/dL_basis (reduced) = %s\n', char(diff_L_red));
report_check('Check 2: dS_1/dL in 6-direction basis (direct chain rule)', ...
             isequal(diff_L_red, sym(0)));

%% =========================================================================
%% Check 3:  dS_1/dH = basis expansion (5 directions)
%%
%% Direct form: chain rule through theta only.  H enters S_1 only
%% through theta = H/G, so dS_1/dH = (1/G) dS_1/dtheta.
%% Basis form: table of ch09d Proposition E.5 column 3.
%% =========================================================================

fprintf('\n=== Check 3: dS_1/dH = basis expansion (5 directions) ===\n');

% Coefficients of dS_1/dH in the basis (in native units
% pref4 = mu^2 k_2 / G^4, so the 3 theta factor is absorbed into the
% prefactor column):
cH_1 = 3 * theta;              % E_1
cH_2 = -3 * theta / 2;         % E_2   [(-theta/2)*3]
cH_3 = -3 * e * theta / 2;     % E_3   [(-theta/2)*3e]
cH_4 = -e * theta / 2;         % E_4   [(-theta/2)*e]
cH_5 = -theta * X02 / 2;       % E_5   [(-theta/2)*X02]
% E_6, E_7, E_8 coefficients are 0 in dS_1/dH.

dS1dH_basis = pref4 * (cH_1*E1 + cH_2*E2 + cH_3*E3 + cH_4*E4 + cH_5*E5);

% Direct form: chain rule through theta alone.  A, B depend on theta;
% alpha, beta, X02 depend on (l, g, e, f); none of {alpha, beta, X02}
% depend on H.
dtheta_dH = 1 / G;
dA_dtheta = 3 * theta;          % d((3 theta^2 - 1)/2)/dtheta
dB_dtheta = -3 * theta;         % d(3(1 - theta^2)/2)/dtheta
dinner_dtheta = dA_dtheta * alpha_inner + (dB_dtheta/6) * beta_inner;
dS1dH_direct = pref3 * dinner_dtheta * dtheta_dH;
% Note: pref3 * (1/G) = pref4.

diff_H = dS1dH_direct - dS1dH_basis;
diff_H_red = expand(diff_H);
for k_iter = 1:15
  prev = diff_H_red;
  diff_H_red = expand(subs(diff_H_red, eta^2, 1 - e^2));
  if isequal(diff_H_red, prev), break; end
end
diff_H_red = simplify(diff_H_red);

fprintf('  dS_1/dH_direct - dS_1/dH_basis (reduced) = %s\n', char(diff_H_red));
report_check('Check 3: dS_1/dH in 5-direction basis (direct chain rule)', ...
             isequal(diff_H_red, sym(0)));

%% =========================================================================
%% Check 4:  dS_1/dG = GSI-inherited basis expansion
%%
%% Direct form: chain rule through all G-dependencies (prefactor, theta,
%%   e, f) of Theorem C.4's closed form.
%% GSI form:  dS_1/dG = -(3/G) S_1 - theta dS_1/dH - (1/eta) dS_1/dL,
%%   built from the three basis reconstructions.
%%
%% This is a CONCRETE cross-check of the abstract GSI (proved in
%% verify_GSI_symbolic.py) specialized to S = S_1 at alpha = 3.  It also
%% tests the Part-(d) coefficient identity of Proposition E.5:
%%   c_k[dS_1/dG] = -(3/G) c_k[S_1] - (theta/G) c_k[dS_1/dH]
%%                - (1/eta) c_k[dS_1/dL]   (in native S_1 column units,
%%   where the theta/G factor converts the dS_1/dH column's native
%%   mu^2 k_2/G^4 prefactor to the dS_1/dG column's mu^2 k_2/G^3 via the
%%   extra 1/G),
%% since both the LHS (via direct chain rule) and the RHS (via the GSI
%% linear combination) must agree coefficient-by-coefficient.
%% =========================================================================

fprintf('\n=== Check 4: dS_1/dG = GSI-inherited basis expansion ===\n');

% GSI-inherited basis form (linear combination of the three reconstructions):
dS1dG_basis_GSI = -(3/G) * S1_basis - theta * dS1dH_basis - (1/eta) * dS1dL_basis;

% Direct form via chain rule on Theorem C.4.  All four chain paths:
%   (i)   Prefactor:  d/dG [mu^2 k_2 / G^3] * inner
%   (ii)  Theta:      pref3 * dA/dtheta * alpha * dtheta/dG
%                   + pref3 * (1/6) dB/dtheta * beta * dtheta/dG
%   (iii) e:          pref3 * d(inner)/d(e) * de/dG       (with eta-chain)
%   (iv)  f:          pref3 * d(inner)/d(f) * df/dG
%
% From ch08:  de/dG|_L = -eta/(L e)             (Theorem D.2)
%             df/dG|_{l,L} = -sin f (2+e cos f)/(L e eta)   (Theorem D.9)
%             dtheta/dG|_H = -theta/G
% And d(prefactor)/dG = -3 mu^2 k_2 / G^4.

de_dG = -eta / (L * e);
df_dG = -sin(f_sym) * (2 + e*cos(f_sym)) / (L * e * eta);
dprefactor_dG = diff(pref3, G);        % = -3 mu^2 k_2 / G^4
dtheta_dG = -theta / G;

dinner_dtheta_th = dA_dtheta * alpha_inner + (dB_dtheta/6) * beta_inner;

dS1dG_direct = dprefactor_dG * inner ...
             + pref3 * ( dinner_dtheta_th * dtheta_dG ...
                       + dinner_de * de_dG ...
                       + dinner_df * df_dG );

% Eliminate kappa from the GSI-inherited basis form (dS1dL_basis is the
% only term carrying kappa, and it appears through E_6, E_7, E_8).
dS1dG_basis_GSI_sub = subs(dS1dG_basis_GSI, kappa, kappa_sub);

diff_G = dS1dG_direct - dS1dG_basis_GSI_sub;
diff_G_red = expand(diff_G);
for k_iter = 1:30
  prev = diff_G_red;
  diff_G_red = expand(subs(diff_G_red, sin(f_sym)^2, 1 - cos(f_sym)^2));
  diff_G_red = expand(subs(diff_G_red, eta^2, 1 - e^2));
  if isequal(diff_G_red, prev), break; end
end
diff_G_red = simplify(diff_G_red);

fprintf('  dS_1/dG_direct - dS_1/dG_basis_GSI (reduced) = %s\n', char(diff_G_red));
report_check('Check 4: dS_1/dG in 8-direction basis (direct chain rule vs GSI)', ...
             isequal(diff_G_red, sym(0)));

fprintf('\n=== verify_ch09d_basis_completeness.m COMPLETE ===\n');
