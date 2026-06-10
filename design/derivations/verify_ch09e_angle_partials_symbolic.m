% verify_ch09e_angle_partials_symbolic.m
%
% PURE SYMBOLIC verification of Proposition E.6 (the three angle partials
% of S_1:  dS_1/dl, dS_1/dg, dS_1/dh).
%
% Claims verified:
%
%   (E.6a)  dS_1/dl = (mu^2 k_2 / L^3) [ (A + B cos(2(f+g)))/kappa^3
%                                         - A/eta^3 ]
%           (L^3-prefactor primary form).  Equivalently, using
%           L^3 = G^3/eta^3 (since G = L eta):
%              dS_1/dl = (mu^2 k_2 / G^3) [ eta^3 (A + B cos(2(f+g)))/kappa^3
%                                            - A ]
%              = (mu^2 k_2 / G^3) { A (eta^3/kappa^3 - 1)
%                                   + B (eta^3/kappa^3) cos(2(f+g)) }.
%
%   (E.6b)  dS_1/dg = (mu^2 k_2 B / G^3) [ cos(2(f+g)) + e cos(f+2g)
%                                           + (e/3) cos(3f+2g)
%                                           + (X02/3) cos(2g) ]
%
%   (E.6c)  dS_1/dh = 0   (h does not appear in S_1; axial symmetry of J_2).
%
% Strategy.  The verifier treats cos(f), sin(f), cos(2g), sin(2g) as
% pure SymPy trig functions.  For (E.6a), the key identity is the
% orbit equation 1 + e cos f = eta^2/kappa, i.e., cos f = (eta^2/kappa - 1)/e.
% We Chebyshev-expand cos(j f + 2g), sin(j f + 2g) into polynomials in
% (cf, sf, c2g, s2g), substitute cos f via the orbit equation, reduce
% sf^2 = 1 - cf^2, and reduce eta^2 = 1 - e^2.  The residue vanishes
% after multiplication by kappa^3 (to clear rational denominators) and
% the three reductions.
%
% Check structure:
%   (a) Direct chain-rule derivative dS_1/dl (via df/dl = eta/kappa^2)
%       equals the L^3-prefactor candidate (E.6a).  Deepest check.
%   (b) Direct SymPy diff(S_1, g) equals the candidate (E.6b).
%   (c) Direct SymPy diff(S_1, h) equals 0.
%   (d) Homological-equation cross-check: n dS_1/dl = F_1 - F_1^*
%       (Theorem C.5(b)), using F_1 from ch05d and F_1^* from ch06d.
%   (e) Symbol-presence inspection: h is absent from the free-symbol
%       list of S_1's closed form.  Independent of (c).
%
% Reductions:  (R1) orbit equation cf = (eta^2 - kappa)/(e kappa),
%              (R2) sf^2 -> 1 - cf^2,  (R3) eta^2 -> 1 - e^2.

pkg load symbolic;
clear; close all;
addpath(fileparts(mfilename('fullpath')));

syms L G H real positive;
syms l g real;
syms h real;
syms e real positive;
syms eta real positive;
syms f_sym real;
syms kappa real positive;
syms mu k_2 real positive;

theta = H / G;
A_expr = (3 * theta^2 - 1) / 2;
B_expr = 3 * (1 - theta^2) / 2;
X02 = (3*e^2 - 2 + 2*eta^3) / e^2;

alpha_inner = f_sym - l + e*sin(f_sym);
beta_inner = 3*sin(2*(f_sym+g)) + 3*e*sin(f_sym+2*g) ...
           + e*sin(3*f_sym+2*g) + X02*sin(2*g);
inner = A_expr*alpha_inner + (B_expr/6)*beta_inner;
pref3 = mu^2 * k_2 / G^3;
S1_direct = pref3 * inner;

% Chebyshev helpers (used to convert cos(j f + 2g), sin(j f + 2g) into
% polynomials in cf, sf, c2g, s2g).  T_j(cos f) = cos(j f);
% sin f * U_{j-1}(cos f) = sin(j f).
function out = chebT(j, x)
  if j == 0, out = sym(1); return; end
  if j == 1, out = x; return; end
  Tp = sym(1); Tc = x;
  for k = 2:j
    Tn = 2*x*Tc - Tp; Tp = Tc; Tc = Tn;
  end
  out = Tc;
end
function out = chebU(jm1, x)
  if jm1 < 0, out = sym(0); return; end
  if jm1 == 0, out = sym(1); return; end
  Up = sym(1); Uc = 2*x;
  for k = 2:jm1
    Un = 2*x*Uc - Up; Up = Uc; Uc = Un;
  end
  out = Uc;
end

syms cf sf c2g s2g real;

% Expansions of cos(j f + 2g), sin(j f + 2g) for j = 1, 2, 3.
cos_jfp2g = cell(3,1);
sin_jfp2g = cell(3,1);
for j = 1:3
  cjf = chebT(j, cf);
  sjf = sf * chebU(j-1, cf);
  cos_jfp2g{j} = cjf * c2g - sjf * s2g;
  sin_jfp2g{j} = sjf * c2g + cjf * s2g;
end
cos_2fp2g = cos_jfp2g{2};   % cos(2(f+g))
sin_2fp2g = sin_jfp2g{2};   % sin(2(f+g))

%% =========================================================================
%% Check (a):  dS_1/dl closed form (Proposition E.6(a))
%%
%% Direct chain rule on Theorem C.4's closed form:
%%    dS_1/dl = (explicit l-derivative) + (d/df) * (df/dl)
%% with df/dl|_e = eta/kappa^2 = (1+e cos f)^2/eta^3 (via orbit equation
%% 1 + e cos f = eta^2/kappa).  Only alpha_inner has an explicit l
%% (the -l term).
%%
%% Candidate: primary L^3-prefactor form
%%    dS_1/dl = (mu^2 k_2 / L^3) [ (A + B cos(2(f+g)))/kappa^3 - A/eta^3 ].
%% Compare to the direct chain rule after Chebyshev expansion, orbit-
%% equation substitution, and (sf^2, eta^2) reductions.
%% =========================================================================

fprintf('\n=== Check (a): dS_1/dl closed form (Proposition E.6(a), L^3 form) ===\n');

% Direct chain rule.
df_dl = eta / kappa^2;
d_inner_dl_explicit = -A_expr;              % d(-l)/dl = -1; coefficient A
d_inner_df = diff(inner, f_sym);
dS1_dl_direct = pref3 * (d_inner_dl_explicit + d_inner_df * df_dl);

% Candidate: primary L^3-prefactor form (E.6a).
pref_L3 = mu^2 * k_2 / L^3;
dS1_dl_candidate = pref_L3 * ( (A_expr + B_expr * cos(2*(f_sym+g))) / kappa^3 ...
                             - A_expr / eta^3 );

% Substitute L = G/eta (from eta = G/L) to put both forms on the same G
% prefactor.  This is an algebraic identity, not an independence-breaker,
% and it is necessary to compare (since the chain-rule form naturally
% lives on 1/G^3).
dS1_dl_candidate_sub = subs(dS1_dl_candidate, L, G / eta);

diff_l = dS1_dl_direct - dS1_dl_candidate_sub;

% Convert cos(j f + 2g), sin(j f + 2g), cos(2g), sin(2g), cos(f),
% sin(f) into the atomic-symbol ring via Chebyshev and direct
% substitution.
diff_l_cheb = diff_l;
diff_l_cheb = subs(diff_l_cheb, cos(f_sym + 2*g),     cos_jfp2g{1});
diff_l_cheb = subs(diff_l_cheb, cos(3*f_sym + 2*g),   cos_jfp2g{3});
diff_l_cheb = subs(diff_l_cheb, cos(2*f_sym + 2*g),   cos_2fp2g);
diff_l_cheb = subs(diff_l_cheb, sin(f_sym + 2*g),     sin_jfp2g{1});
diff_l_cheb = subs(diff_l_cheb, sin(3*f_sym + 2*g),   sin_jfp2g{3});
diff_l_cheb = subs(diff_l_cheb, sin(2*f_sym + 2*g),   sin_2fp2g);
diff_l_cheb = subs(diff_l_cheb, cos(2*g),             c2g);
diff_l_cheb = subs(diff_l_cheb, sin(2*g),             s2g);
diff_l_cheb = subs(diff_l_cheb, cos(f_sym),           cf);
diff_l_cheb = subs(diff_l_cheb, sin(f_sym),           sf);

% Multiply by kappa^3 to clear the dominant denominator.
diff_l_cleared = expand(diff_l_cheb * kappa^3);

% Apply the orbit equation cf = (eta^2 - kappa)/(e kappa), which is
% exactly cos f = (eta^2/kappa - 1)/e.  This substitutes cos f in terms
% of (eta, kappa, e) so that kappa no longer competes with cf in the
% denominator.
cos_f_sub = (eta^2 - kappa) / (e * kappa);
diff_l_cleared = subs(diff_l_cleared, cf, cos_f_sub);
diff_l_cleared = expand(diff_l_cleared);

% Iterate: reduce sf^2 -> 1 - cf^2 (then cf -> orbit-equation form),
% eta^2 -> 1 - e^2.
for k_iter = 1:30
  prev = diff_l_cleared;
  diff_l_cleared = expand(subs(diff_l_cleared, sf^2, 1 - cf^2));
  diff_l_cleared = expand(subs(diff_l_cleared, cf, cos_f_sub));
  diff_l_cleared = expand(subs(diff_l_cleared, eta^2, 1 - e^2));
  if isequal(diff_l_cleared, prev), break; end
end
diff_l_cleared = simplify(diff_l_cleared);

fprintf('  (dS_1/dl_direct - dS_1/dl_candidate)*kappa^3, fully reduced = %s\n', ...
        char(diff_l_cleared));
report_check('(E.6a) dS_1/dl closed form (L^3-prefactor)', ...
             isequal(diff_l_cleared, sym(0)));

%% =========================================================================
%% Check (b):  dS_1/dg closed form (Proposition E.6(b))
%%
%% Direct SymPy diff(S_1, g) equals the four-cosine candidate.
%% No orbit-equation substitution needed (f has no g-dependence; g
%% enters only through the explicit sin(... + 2g) terms).
%% =========================================================================

fprintf('\n=== Check (b): dS_1/dg closed form (Proposition E.6(b)) ===\n');

d_inner_dg = diff(inner, g);
dS1_dg_direct = pref3 * d_inner_dg;

% Candidate:
pref3B = pref3 * B_expr;
dS1_dg_candidate = pref3B * ( cos(2*(f_sym+g)) + e*cos(f_sym+2*g) ...
                             + (e/3)*cos(3*f_sym+2*g) + (X02/3)*cos(2*g) );

diff_g = expand(dS1_dg_direct - dS1_dg_candidate);
diff_g = simplify(diff_g);
fprintf('  dS_1/dg_direct - dS_1/dg_candidate = %s\n', char(diff_g));
report_check('(E.6b) dS_1/dg closed form', isequal(diff_g, sym(0)));

%% =========================================================================
%% Check (c):  dS_1/dh = 0 (Proposition E.6(c))
%%
%% h appears nowhere in S_1 (axial symmetry of J_2); diff(S_1, h) = 0.
%% =========================================================================

fprintf('\n=== Check (c): dS_1/dh = 0 (Proposition E.6(c)) ===\n');

dS1_dh_direct = diff(S1_direct, h);
fprintf('  dS_1/dh_direct = %s\n', char(dS1_dh_direct));
report_check('(E.6c) dS_1/dh = 0', isequal(dS1_dh_direct, sym(0)));

%% =========================================================================
%% Check (d):  Homological-equation cross-check
%%
%%   n dS_1/dl = F_1 - F_1^*   (Theorem C.5(b), ch07c).
%%
%% F_1 = (mu k_2 / a^3) (a/r)^3 (A + B cos(2(f+g)))      (ch05d Theorem A.17)
%%     = (mu k_2 / a^3) (1/kappa^3) (A + B cos(2(f+g)))
%% F_1^* = mu k_2 A / (a^3 eta^3)                        (ch06d Theorem B.5.1)
%% n = mu^2 / L^3,  a = L^2/mu,  so a^3 = L^6/mu^3 and n a^3 = L^3/mu.
%%
%% IMPORTANT: To make this a non-tautological test of Theorem C.5(b), we use
%% the CHAIN-RULE-DERIVED form (dS1_dl_direct, from Theorem C.4 via chain
%% rule), NOT the candidate (E.6a).  The candidate (E.6a) was constructed
%% (in derivation I, ch09e §2) as (F_1 - F_1^*)/n, so using it on the LHS
%% would reduce Check (d) to the trivial identity n * ((F_1 - F_1^*)/n) =
%% F_1 - F_1^*.
%%
%% By using dS1_dl_direct, Check (d) becomes an independent consistency test
%% between ch07c Theorem C.4 (closed form of S_1) and Theorem C.5(b)
%% (homological equation satisfaction):  does chain-rule differentiation of
%% Theorem C.4's closed form actually produce (F_1 - F_1^*)/n?
%%
%% Note: this is NOT circular with Check (a).  Check (a) tests
%% dS1_dl_direct = dS1_dl_candidate.  Check (d) tests
%% n * dS1_dl_direct = F_1 - F_1^* using independently-sourced F_1 and F_1^*
%% from ch05d/ch06d.  The combination of Checks (a) and (d) pins down the
%% closed form (E.6a) from two independent routes.
%% =========================================================================

fprintf('\n=== Check (d): n dS_1/dl = F_1 - F_1^* (Theorem C.5(b) cross-check) ===\n');

a_cubed = L^6 / mu^3;
cos2fpg_expr = cos(2*(f_sym + g));
F1 = (mu * k_2 / a_cubed) * (1/kappa^3) * (A_expr + B_expr * cos2fpg_expr);
F1_star = mu * k_2 * A_expr / (a_cubed * eta^3);
n_expr = mu^2 / L^3;

% Use the chain-rule-derived form (not the candidate), so Check (d) is an
% independent test of Theorem C.5(b) rather than a tautology.
% dS1_dl_direct lives on the G^3 prefactor; F_1 - F_1^* lives on the L^3
% prefactor (via a^3 = L^6/mu^3).  Unify by substituting L = G/eta (from
% eta = G/L), matching the approach of Check (a).
LHS_C5b = n_expr * dS1_dl_direct;         % n * (chain-rule dS_1/dl)
RHS_C5b = F1 - F1_star;
LHS_C5b = subs(LHS_C5b, L, G / eta);
RHS_C5b = subs(RHS_C5b, L, G / eta);

% Clear the common denominator (kappa^3) and reduce via orbit equation
% (cos f = (eta^2 - kappa)/(e kappa)), sf^2 -> 1 - cf^2, eta^2 -> 1 - e^2.
% Using the same atomic-symbol substitutions as Check (a); dS1_dl_direct
% contains cos(j f + 2g), sin(j f + 2g) for j = 1, 2, 3 from the B-bracket
% in dinner/df (via the Chebyshev harmonics).
diff_d_raw = LHS_C5b - RHS_C5b;
diff_d_cheb = diff_d_raw;
diff_d_cheb = subs(diff_d_cheb, cos(f_sym + 2*g),     cos_jfp2g{1});
diff_d_cheb = subs(diff_d_cheb, cos(3*f_sym + 2*g),   cos_jfp2g{3});
diff_d_cheb = subs(diff_d_cheb, cos(2*f_sym + 2*g),   cos_2fp2g);
diff_d_cheb = subs(diff_d_cheb, sin(f_sym + 2*g),     sin_jfp2g{1});
diff_d_cheb = subs(diff_d_cheb, sin(3*f_sym + 2*g),   sin_jfp2g{3});
diff_d_cheb = subs(diff_d_cheb, sin(2*f_sym + 2*g),   sin_2fp2g);
diff_d_cheb = subs(diff_d_cheb, cos(2*g),             c2g);
diff_d_cheb = subs(diff_d_cheb, sin(2*g),             s2g);
diff_d_cheb = subs(diff_d_cheb, cos(f_sym),           cf);
diff_d_cheb = subs(diff_d_cheb, sin(f_sym),           sf);
diff_d_cleared = expand(diff_d_cheb * kappa^3);
diff_d_cleared = subs(diff_d_cleared, cf, cos_f_sub);
diff_d_cleared = expand(diff_d_cleared);
for k_iter = 1:30
  prev = diff_d_cleared;
  diff_d_cleared = expand(subs(diff_d_cleared, sf^2, 1 - cf^2));
  diff_d_cleared = expand(subs(diff_d_cleared, cf, cos_f_sub));
  diff_d_cleared = expand(subs(diff_d_cleared, eta^2, 1 - e^2));
  if isequal(diff_d_cleared, prev), break; end
end
diff_d_cleared = simplify(diff_d_cleared);

fprintf('  (n dS_1/dl_direct - (F_1 - F_1^*))*kappa^3, fully reduced = %s\n', ...
        char(diff_d_cleared));
report_check('(d) Theorem C.5(b) consistency: n dS_1/dl_direct = F_1 - F_1^*', ...
             isequal(diff_d_cleared, sym(0)));

%% =========================================================================
%% Check (e):  h is absent from S_1 (independent of Check (c))
%%
%% Inspect the free-symbol list of S_1's closed form and confirm h is
%% not among them.  This is a different kind of verification: (c) is a
%% computational derivative test; (e) is a static symbol-presence
%% inspection.
%% =========================================================================

fprintf('\n=== Check (e): h does not appear in S_1''s closed form ===\n');

S1_free_symbols = symvar(S1_direct);
S1_free_names = {};
if iscell(S1_free_symbols)
  for k = 1:numel(S1_free_symbols)
    S1_free_names{end+1} = char(S1_free_symbols{k});
  end
else
  for k = 1:numel(S1_free_symbols)
    S1_free_names{end+1} = char(S1_free_symbols(k));
  end
end

h_present = false;
for k = 1:numel(S1_free_names)
  if strcmp(S1_free_names{k}, 'h')
    h_present = true;
    break;
  end
end
fprintf('  S_1 free symbols: %s\n', strjoin(S1_free_names, ', '));
report_check('(e) h is absent from S_1 free-symbol list', ~h_present);

fprintf('\n=== verify_ch09e_angle_partials_symbolic.m COMPLETE ===\n');
