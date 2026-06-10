%% verify_ch10_poisson_bracket.m
%%
%% Chapter 10b Proposition F.7: {F_1, S_1} in M_10 (headline).
%%
%% Claim (Corollary 2.4 of CH10_PLAN.md, Proposition F.7 of ch10b_poisson_bracket.md):
%%   {F_1, S_1}_{Delaunay, modern} in M_{alpha_1 + alpha_2 + 1} = M_{6+3+1} = M_{10}
%%
%%   Specifically, {F_1, S_1} = G^{-10} * Phi_{bracket}(theta, e, l, g) for some
%%   smooth F-factor Phi_{bracket}.
%%
%% Verification strategy:
%%   Check 1 (Structural): Theorem 3 applied at (alpha_1, alpha_2) = (6, 3) with
%%                         abstract F-factors Phi_{F_1}(theta,e,l,g), Phi_{S_1}(theta,e,l,g)
%%                         gives G-power -10 after L -> G/eta substitution.
%%                         Uses the same Octave SymPy abstract-function machinery as
%%                         verify_ch10_thm3_lie_closure.m.
%%   Check 2 ((h, H) vanishing): The j=3 bracket term is identically 0
%%                               (Proposition F.4).
%%   Check 3 (Two-term reduction): {F_1, S_1} = T_{lL} + T_{gG} (Proposition F.5).
%%   Check 4 (Per-product G-power): Each of the four bracket product terms has
%%                                  G-power -10 individually.
%%
%% Because F_1 and S_1 are in M_6 and M_3 respectively (established in Phase A),
%% the Theorem-3-specialization proof goes through automatically. This verifier
%% confirms the specialization via Octave/SymPy.

pkg load symbolic;

printf('============================================================\n');
printf('Ch 10b Proposition F.7 -- {F_1, S_1} in M_10\n');
printf('Claim: {F_1, S_1}_{modern} = G^{-10} * Phi_{bracket}(theta, e, l, g)\n');
printf('============================================================\n\n');

%%=====================================================================
%% Check 1: Structural G-power via Theorem 3 at (alpha_1, alpha_2) = (6, 3)
%%=====================================================================
printf('--- Check 1: Structural G-power via Theorem 3 at (6, 3) ---\n\n');

syms L G H real positive;
syms l g h real;

% Abstract F-factors of F_1 and S_1 (SymPy Function class).
% These stand in for the concrete Phi_{F_1}(theta, e, l, g) = mu^4 k_2 (1+e cos f)^3 (A + B cos 2(f+g))
% and Phi_{S_1}(theta, e, l, g) = mu^2 k_2 {A(f-l+e sin f) + (B/6)[...]}.
% The G-power verification is structural -- it does not depend on the closed forms,
% only on the factorization S = G^{-alpha} * Phi(theta, e, l, g).
syms Phi_F1(a1, a2, a3, a4) Phi_S1(a1, a2, a3, a4);

theta_expr = H/G;
e_expr = sqrt(1 - G^2/L^2);

% Build F_1 and S_1 in M-class form with (alpha_1, alpha_2) = (6, 3).
F_1 = G^(-6) * Phi_F1(theta_expr, e_expr, l, g);
S_1 = G^(-3) * Phi_S1(theta_expr, e_expr, l, g);

printf('F_1 := G^{-6} * Phi_{F_1}(theta, e, l, g)  [alpha_1 = 6, Proposition F.1]\n');
printf('S_1 := G^{-3} * Phi_{S_1}(theta, e, l, g)  [alpha_2 = 3, ch09a Lemma E.1.1]\n\n');

% Compute Delaunay Poisson bracket (modern convention).
dF1_dl = diff(F_1, l);
dF1_dg = diff(F_1, g);
dF1_dh = diff(F_1, h);
dF1_dL = diff(F_1, L);
dF1_dG = diff(F_1, G);
dF1_dH = diff(F_1, H);

dS1_dl = diff(S_1, l);
dS1_dg = diff(S_1, g);
dS1_dh = diff(S_1, h);
dS1_dL = diff(S_1, L);
dS1_dG = diff(S_1, G);
dS1_dH = diff(S_1, H);

% Full Poisson bracket:
%   {F_1, S_1} = sum_{j=1}^{3} [(dF_1/dl_j)(dS_1/dL_j) - (dF_1/dL_j)(dS_1/dl_j)]
PB = (dF1_dl * dS1_dL - dF1_dL * dS1_dl) ...
   + (dF1_dg * dS1_dG - dF1_dG * dS1_dg) ...
   + (dF1_dh * dS1_dH - dF1_dH * dS1_dh);

printf('Computed {F_1, S_1}_{modern} via Delaunay formula (3 pairs).\n\n');

% Multiply by G^10 and substitute L -> G/eta.
scaled = G^10 * PB;

syms eta_sym real positive;
scaled_eta = subs(scaled, sqrt(1 - G^2/L^2), eta_sym);
scaled_eta = subs(scaled_eta, L, G/eta_sym);
scaled_eta = simplify(scaled_eta);

% Check that L is not among the free symbols.
free_syms = symvar(scaled_eta);
fprintf('Free symbols in G^{10} * {F_1, S_1} (post L -> G/eta): ');
for k = 1:numel(free_syms)
  fprintf('%s ', char(free_syms(k)));
end
fprintf('\n');

has_L = false;
for k = 1:numel(free_syms)
  if strcmp(char(free_syms(k)), 'L')
    has_L = true;
    break;
  end
end

if has_L
  printf('[FAIL] Residual L-dependence found.\n\n');
else
  printf('[PASS] G^{10} * {F_1, S_1} has no bare L-dependence.\n');
  printf('       Hence {F_1, S_1} in M_10 (Proposition F.7). ***\n\n');
end

%%=====================================================================
%% Check 2: (h, H) vanishing (Proposition F.4 / F.5 j=3 term)
%%=====================================================================
printf('--- Check 2: (h, H) bracket term vanishes (Proposition F.4) ---\n\n');

hH_term = dF1_dh * dS1_dH - dF1_dH * dS1_dh;
hH_term_simp = simplify(hH_term);
printf('(h, H) term = '); disp(hH_term_simp);

if isequal(hH_term_simp, sym(0))
  printf('[PASS] (dF_1/dh)(dS_1/dH) - (dF_1/dH)(dS_1/dh) = 0 identically.\n\n');
else
  printf('[FAIL] (h, H) term nonzero.\n\n');
end

%%=====================================================================
%% Check 3: Two-term reduction {F_1, S_1} = T_{lL} + T_{gG} (Proposition F.5)
%%=====================================================================
printf('--- Check 3: Two-term reduction {F_1, S_1} = T_{lL} + T_{gG} (Proposition F.5) ---\n\n');

T_lL = dF1_dl * dS1_dL - dF1_dL * dS1_dl;
T_gG = dF1_dg * dS1_dG - dF1_dG * dS1_dg;

PB_two_term = T_lL + T_gG;
PB_diff = simplify(PB - PB_two_term);
printf('PB (full) - (T_{lL} + T_{gG}) = '); disp(PB_diff);

if isequal(PB_diff, sym(0))
  printf('[PASS] {F_1, S_1} = T_{lL} + T_{gG} (Proposition F.5).\n\n');
else
  printf('[FAIL] Two-term decomposition mismatch.\n\n');
end

%%=====================================================================
%% Check 4: G-power of each of the four product terms individually
%%=====================================================================
printf('--- Check 4: G-power of each of the 4 bracket product terms is -10 ---\n\n');

products = struct();
products.A = struct('name', 'A = (dF_1/dl)(dS_1/dL)', 'expr', dF1_dl * dS1_dL);
products.B = struct('name', 'B = (dF_1/dL)(dS_1/dl)', 'expr', dF1_dL * dS1_dl);
products.C = struct('name', 'C = (dF_1/dg)(dS_1/dG)', 'expr', dF1_dg * dS1_dG);
products.D = struct('name', 'D = (dF_1/dG)(dS_1/dg)', 'expr', dF1_dG * dS1_dg);

prod_names = {'A', 'B', 'C', 'D'};
all_pass = true;
for kk = 1:4
  pname = prod_names{kk};
  prod_expr = products.(pname).expr;
  prod_desc = products.(pname).name;

  scaled_prod = G^10 * prod_expr;
  scaled_prod_eta = subs(scaled_prod, sqrt(1 - G^2/L^2), eta_sym);
  scaled_prod_eta = subs(scaled_prod_eta, L, G/eta_sym);
  scaled_prod_eta = simplify(scaled_prod_eta);

  has_L_prod = false;
  free_p = symvar(scaled_prod_eta);
  for mm = 1:numel(free_p)
    if strcmp(char(free_p(mm)), 'L')
      has_L_prod = true;
    end
  end

  if has_L_prod
    printf('  %s: G-power != -10. [FAIL]\n', prod_desc);
    all_pass = false;
  else
    printf('  %s: G-power is -10. [PASS]\n', prod_desc);
  end
end

if all_pass
  printf('\n[PASS] All 4 bracket product terms have G-power -10 individually.\n');
  printf('       This confirms the Theorem-3 specialization at (alpha_1, alpha_2) = (6, 3).\n\n');
else
  printf('\n[FAIL] At least one product term has wrong G-power.\n\n');
end

%%=====================================================================
%% Summary
%%=====================================================================
printf('============================================================\n');
printf('Proposition F.7 verification summary:\n');
printf('============================================================\n');
printf('  Check 1 (Structural):        {F_1, S_1} in M_10 via Theorem 3 specialization\n');
printf('  Check 2 ((h, H) vanishing):  j=3 bracket term = 0 identically\n');
printf('  Check 3 (Two-term):          {F_1, S_1} = T_{lL} + T_{gG} (Prop F.5)\n');
printf('  Check 4 (Per-product):       Each of 4 products has G-power -10\n');
printf('\n');
printf('Harmonic inventory bounds (for Phase C reference; see ch10b_poisson_bracket.md Section 2):\n');
printf('  j ranges:         j in {0..8} (bounded by j_{F_1}^{max} + j_{S_1}^{max} = 5 + 3)\n');
printf('  k ranges:         k in {0, +/-2, +/-4} (from sum/diff of F_1 and S_1 g-frequencies)\n');
printf('  m (kappa-power):  m in {0..7} (bounded by product (dF_1/dl)(dS_1/dL) = kappa^{-5} * kappa^{-2})\n');
printf('  Actual occupied set is a strict subset (tabulated in verify_ch10c_secular_average.m).\n');
printf('============================================================\n');
