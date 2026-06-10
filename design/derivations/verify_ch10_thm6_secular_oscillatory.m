%% verify_ch10_thm6_secular_oscillatory.m
%%
%% Chapter 10 -- Theorem 6 (Secular-oscillatory bracket lemma) verification.
%%
%% Claim: If F_s in M_alpha is secular (dF_s/dl = dF_s/dg = 0) and
%%        S_o in M_beta is l-mean-zero (<S_o>_l = 0), then
%%           <{F_s, S_o}>_l = 0.
%%
%% Verification strategy (concrete F_s and S_o satisfying the hypotheses):
%%   F_s = G^(-alpha) * theta^2 * eta^3  (secular: depends only on theta, e)
%%   S_o = G^(-beta) * [c1(theta, e) cos(l) + c2(theta, e) cos(2l + 3g)
%%                      + c3(theta, e) sin(l + 2g)]
%%   where c_k are chosen so that <S_o>_l = 0 (the l-averaged sum is zero).
%%
%% Because each term of S_o is either a pure l-oscillatory cosine or a
%% cos/sin of (jl + kg) with j != 0, the l-average is zero for each term
%% individually.  Hence <S_o>_l = 0.
%%
%% Compute {F_s, S_o}_modern = sum_j [(dF_s/dl_j)(dS_o/dL_j) - (dF_s/dL_j)(dS_o/dl_j)].
%% Take l-average symbolically.
%% Verify result = 0.

pkg load symbolic;

printf('============================================================\n');
printf('Ch 10 Theorem 6 -- Secular-oscillatory bracket lemma\n');
printf('Verify <{F_s, S_o}>_l = 0 for secular F_s, l-mean-zero S_o.\n');
printf('============================================================\n\n');

syms L G H real positive;
syms l g real;
syms alpha_sym beta_sym real;

theta_expr = H / G;
e_expr     = sqrt(1 - G^2 / L^2);
eta_expr   = sqrt(1 - e_expr^2);  % = G/L with G, L positive

% ----- F_s: secular (no l, g, h dependence) -----
% Concrete choice: F_s = G^(-alpha) * theta^2 * eta^3.
F_s = G^(-alpha_sym) * theta_expr^2 * eta_expr^3;

printf('F_s := G^(-alpha) * theta^2 * eta^3  (secular)\n');
printf('  dF_s/dl = %s\n', char(simplify(diff(F_s, l))));  % should be 0
printf('  dF_s/dg = %s\n', char(simplify(diff(F_s, g))));  % should be 0
printf('  dF_s/dh is omitted (H-only dependence; no h in S).\n\n');

% ----- S_o: explicit l-mean-zero form -----
% Pick: S_o = G^(-beta) * [sin(l) + sin(l + 2g) * theta^2 + cos(2l + 3g) * e]
% Each term has l-mean zero individually (integrals of sin/cos
% over [0, 2pi] at nonzero frequency vanish).
% So <S_o>_l = 0 by linearity.
S_o = G^(-beta_sym) * ( sin(l) ...
                       + sin(l + 2*g) * theta_expr^2 ...
                       + cos(2*l + 3*g) * e_expr );

% Verify <S_o>_l = 0.
S_o_avg = (1/(2*sym(pi))) * int(S_o, l, 0, 2*sym(pi));
S_o_avg_simp = simplify(S_o_avg);
printf('S_o := G^(-beta) * [sin(l) + sin(l+2g)*theta^2 + cos(2l+3g)*e]\n');
printf('<S_o>_l = %s\n', char(S_o_avg_simp));
if isequal(S_o_avg_simp, sym(0))
  printf('  [PASS] <S_o>_l = 0 (S_o is l-mean-zero by construction)\n\n');
else
  printf('  [FAIL] S_o has nonzero l-mean; test setup invalid.\n\n');
end

% ----- Compute {F_s, S_o}_{modern} -----
% Bracket: sum_j [(dF/dl_j)(dG/dL_j) - (dF/dL_j)(dG/dl_j)]
%   j=1 (l, L): (dF/dl)(dS/dL) - (dF/dL)(dS/dl)
%   j=2 (g, G): (dF/dg)(dS/dG) - (dF/dG)(dS/dg)
%   j=3 (h, H): dF/dh = 0 (F is H-dep only); dS/dh = 0 (S in M_beta).

dF_s_dl = diff(F_s, l);
dF_s_dg = diff(F_s, g);
dF_s_dL = diff(F_s, L);
dF_s_dG = diff(F_s, G);
dF_s_dH = diff(F_s, H);

dS_o_dl = diff(S_o, l);
dS_o_dg = diff(S_o, g);
dS_o_dL = diff(S_o, L);
dS_o_dG = diff(S_o, G);

PB = (dF_s_dl * dS_o_dL - dF_s_dL * dS_o_dl) ...
   + (dF_s_dg * dS_o_dG - dF_s_dG * dS_o_dg);

% (h, H) pair is zero: dF/dh = 0 by F_s being secular, dS/dh = 0 by M_beta.

printf('{F_s, S_o}_{modern} computed (l-dependent expression).\n\n');

% ----- Take l-average -----
% Integrate over l in [0, 2*pi] and divide by 2*pi.
PB_avg = (1/(2*sym(pi))) * int(PB, l, 0, 2*sym(pi));
PB_avg_simp = simplify(PB_avg);

printf('<{F_s, S_o}>_l = %s\n\n', char(PB_avg_simp));

if isequal(PB_avg_simp, sym(0))
  printf('[PASS] <{F_s, S_o}>_l = 0 (Theorem 6 verified).\n');
else
  printf('[FAIL] <{F_s, S_o}>_l is nonzero -- Theorem 6 does not hold here.\n');
  printf('       Check secular-vs-oscillatory hypothesis carefully.\n');
end

%% ============================================================
%% Additional test: abstract F_s (function of L, G, H only) and abstract S_o
%% ============================================================

printf('\n--- Test 2: abstract secular F_s (function of L, G, H only) ---\n');

syms F_abs_s(m1, m2, m3);
F_s_abs = F_abs_s(L, G, H);
% F_s_abs depends on all three momenta but NOT on l, g, h. So it's secular.

% Use the same S_o as before.
dFs_dL_abs = diff(F_s_abs, L);
dFs_dG_abs = diff(F_s_abs, G);
dFs_dH_abs = diff(F_s_abs, H);
% Angle partials are zero:
% dFs/dl = 0, dFs/dg = 0, dFs/dh = 0.

PB_abs = 0 * dS_o_dL ...        % dF_s/dl * dS_o/dL
       - dFs_dL_abs * dS_o_dl ... % - dF_s/dL * dS_o/dl
       + 0 * dS_o_dG ...        % dF_s/dg * dS_o/dG
       - dFs_dG_abs * dS_o_dg;    % - dF_s/dG * dS_o/dg
% h, H pair: 0 since both angle partials are zero by construction.

PB_abs_avg = (1/(2*sym(pi))) * int(PB_abs, l, 0, 2*sym(pi));
PB_abs_avg_simp = simplify(PB_abs_avg);

printf('<{F_s(L,G,H), S_o}>_l = %s\n', char(PB_abs_avg_simp));

if isequal(PB_abs_avg_simp, sym(0))
  printf('[PASS] Abstract-F_s case: Theorem 6 holds.\n');
else
  printf('[FAIL] Abstract-F_s case fails.\n');
end

printf('\n============================================================\n');
printf('Theorem 6 verification summary\n');
printf('============================================================\n');
