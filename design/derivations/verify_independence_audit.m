% verify_independence_audit.m
%
% Variable-dependence audit for all formulas in ch06e, ch10e, ch11b, ch11d,
% ch10c_addendum, and related files. Verifies each formula's explicit
% variable dependencies match the claimed M_alpha class membership and the
% independence assertions (l-independent, g-independent, g-dependent, etc.).
%
% Method: symbolically differentiate each formula with respect to each
% primitive Delaunay variable (L, G, H, l, g, h) and check whether the
% result is identically zero where the formula is claimed to be independent,
% and non-zero where it is claimed to be dependent.
%
% Special care: eta = sqrt(1-e^2) = G/L. When we write formulas with eta
% explicitly, differentiation chains through both G and L.

pkg load symbolic;

printf('=========================================================\n');
printf('Variable-dependence audit for new-chapter formulas\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function print_check(label, tol_str, obs_str, pass_flag)
  printf('  %s\n', label);
  printf('    Expected: %s\n', tol_str);
  printf('    Observed: %s\n', obs_str);
  if pass_flag
    printf('    Status:   PASS\n');
  else
    printf('    Status:   FAIL\n');
  end
  printf('\n');
endfunction

% ---------------------------------------------------------------
% Primitive Delaunay symbols
% ---------------------------------------------------------------
syms L G H l g h mu k2 positive;
% Derived dimensionless variables
eta_expr = G/L;
e_expr = sqrt(1 - (G/L)^2);
theta_expr = H/G;
A_expr = (3*theta_expr^2 - 1)/2;
B_expr = 3*(1 - theta_expr^2)/2;

%% =======================================================
%% Part 1: ch06e (B.5.3-5) - First-order secular-rate partials of F_1^*
%% F_1^* = mu^4 k_2 A(theta) eta^3 / G^6 in M_6
%% Claim: l-independent, g-independent, h-independent
%% =======================================================
printf('--- Part 1: ch06e F_1^* partials (l, g, h independence) ---\n');

F1star = mu^4 * k2 * A_expr * eta_expr^3 / G^6;

% B.5.3: ∂F_1^*/∂L = -3 mu^4 k_2 A eta^4 / G^7
B53 = -3 * mu^4 * k2 * A_expr * eta_expr^4 / G^7;
% B.5.4: ∂F_1^*/∂G = -3 mu^4 k_2 eta^3 (5 theta^2 - 1) / (2 G^7)
B54 = -3 * mu^4 * k2 * eta_expr^3 * (5*theta_expr^2 - 1) / (2 * G^7);
% B.5.5: ∂F_1^*/∂H = 3 mu^4 k_2 theta eta^3 / G^7
B55 = 3 * mu^4 * k2 * theta_expr * eta_expr^3 / G^7;

% Check B.5.3 via direct diff
R1a = simplify(diff(F1star, L) - B53);
R1b = simplify(diff(F1star, G) - B54);
R1c = simplify(diff(F1star, H) - B55);
p1a = isequal(R1a, sym(0));
p1b = isequal(R1b, sym(0));
p1c = isequal(R1c, sym(0));

% l, g, h independence checks
R_B53_l = simplify(diff(B53, l));
R_B53_g = simplify(diff(B53, g));
R_B53_h = simplify(diff(B53, h));
R_B54_l = simplify(diff(B54, l));
R_B54_g = simplify(diff(B54, g));
R_B54_h = simplify(diff(B54, h));
R_B55_l = simplify(diff(B55, l));
R_B55_g = simplify(diff(B55, g));
R_B55_h = simplify(diff(B55, h));

all_B_zeros = isequal(R_B53_l, sym(0)) && isequal(R_B53_g, sym(0)) && isequal(R_B53_h, sym(0)) && ...
              isequal(R_B54_l, sym(0)) && isequal(R_B54_g, sym(0)) && isequal(R_B54_h, sym(0)) && ...
              isequal(R_B55_l, sym(0)) && isequal(R_B55_g, sym(0)) && isequal(R_B55_h, sym(0));

print_check('B.5.3 formula correctness: diff(F_1^*, L) - B53 = 0', ...
  'residual = 0', sprintf('%s', char(R1a)), p1a);
if p1a; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('B.5.4 formula correctness: diff(F_1^*, G) - B54 = 0', ...
  'residual = 0', sprintf('%s', char(R1b)), p1b);
if p1b; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('B.5.5 formula correctness: diff(F_1^*, H) - B55 = 0', ...
  'residual = 0', sprintf('%s', char(R1c)), p1c);
if p1c; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('B.5.3-5 independence: no (l, g, h) dependence', ...
  'diff(B, x) = 0 for all B in {B53, B54, B55}, x in {l, g, h}', ...
  sprintf('all 9 residuals = 0: %d', all_B_zeros), all_B_zeros);
if all_B_zeros; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 2: ch11b (G.6-G.9) - S_1^* partials
%% S_1^* = G^{-3} * Phi_{S_1^*}(theta, e, eta, g)
%% With Phi = (mu^2 k_2 / [(5 theta^2 - 1) eta^2]) [AB e^2(6+e^2)/4 sin(2g) + B^2 e^4/128 sin(4g)]
%% Claim: l-independent, h-independent, g-dependent
%% =======================================================
printf('--- Part 2: ch11b S_1^* partials (l, h independence; g dependence) ---\n');

% Build Phi_{S_1^*} using derived e, eta, theta expressions (so diff through L, G, H works)
Phi_S1star = mu^2 * k2 / ((5*theta_expr^2 - 1) * eta_expr^2) * ...
             (A_expr * B_expr * e_expr^2 * (6 + e_expr^2)/4 * sin(2*g) + ...
              B_expr^2 * e_expr^4 / 128 * sin(4*g));
S1star = Phi_S1star / G^3;

% l-independence: ∂S_1^*/∂l should be 0 identically (G.9a)
R_S1s_l = simplify(diff(S1star, l));
p_l_indep = isequal(R_S1s_l, sym(0));

% h-independence: ∂S_1^*/∂h should be 0 identically (G.9c)
R_S1s_h = simplify(diff(S1star, h));
p_h_indep = isequal(R_S1s_h, sym(0));

% g-dependence: ∂S_1^*/∂g should be NONZERO
R_S1s_g = simplify(diff(S1star, g));
g_depends = ~isequal(R_S1s_g, sym(0));

print_check('G.9a: ∂S_1^*/∂l = 0 (l-independence)', ...
  'residual = 0', sprintf('%s', char(R_S1s_l)), p_l_indep);
if p_l_indep; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('G.9c: ∂S_1^*/∂h = 0 (h-independence, axial symmetry)', ...
  'residual = 0', sprintf('%s', char(R_S1s_h)), p_h_indep);
if p_h_indep; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

print_check('G.9b: ∂S_1^*/∂g ≠ 0 (g-dependence)', ...
  'nonzero (g-harmonic content present)', ...
  sprintf('is nonzero: %d', g_depends), g_depends);
if g_depends; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

% Verify G.9b closed form directly
dS1s_dg_form = mu^2 * k2 / ((5*theta_expr^2 - 1) * eta_expr^2 * G^3) * ...
               (A_expr * B_expr * e_expr^2 * (6 + e_expr^2)/2 * cos(2*g) + ...
                B_expr^2 * e_expr^4 / 32 * cos(4*g));
R_G9b = simplify(R_S1s_g - dS1s_dg_form);
p_G9b = isequal(R_G9b, sym(0));
print_check('G.9b closed form match: diff(S_1^*, g) - boxed form = 0', ...
  'residual = 0', sprintf('%s', char(R_G9b)), p_G9b);
if p_G9b; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 3: ch10e (F.15) - F_2^{**} closed form
%% c_0^{(T)} = (3 mu^6 k_2^2 eta / (2 G^10)) [(A^2 + B^2/2)(8 + 24e^2 + 3e^4)/8 - A^2 eta^3]
%% Claim: l-independent, g-independent, h-independent (fully secular)
%% =======================================================
printf('--- Part 3: ch10e F_2^{**} c_0^{(T)} (l, g, h independence) ---\n');

Q_expr = A_expr^2 + B_expr^2/2;
R_expr = A_expr^2;
P_expr = 8 + 24*e_expr^2 + 3*e_expr^4;
c0T = (3 * mu^6 * k2^2 * eta_expr / (2 * G^10)) * ...
      (Q_expr * P_expr / 8 - R_expr * eta_expr^3);

% l, g, h independence
R_c0T_l = simplify(diff(c0T, l));
R_c0T_g = simplify(diff(c0T, g));
R_c0T_h = simplify(diff(c0T, h));
c0T_secular = isequal(R_c0T_l, sym(0)) && isequal(R_c0T_g, sym(0)) && isequal(R_c0T_h, sym(0));

print_check('F.15 c_0^{(T)}: l, g, h independence (fully secular)', ...
  'all 3 residuals = 0', ...
  sprintf('all zero: %d', c0T_secular), c0T_secular);
if c0T_secular; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

% Dependence on L, G, H: should be nonzero (it's a function of these)
c0T_L_depends = ~isequal(simplify(diff(c0T, L)), sym(0));
c0T_G_depends = ~isequal(simplify(diff(c0T, G)), sym(0));
c0T_H_depends = ~isequal(simplify(diff(c0T, H)), sym(0));
LGH_deps = c0T_L_depends && c0T_G_depends && c0T_H_depends;

print_check('F.15 c_0^{(T)}: L, G, H dependence (required for secular rates)', ...
  'all 3 dependencies nonzero', ...
  sprintf('all nonzero: %d', LGH_deps), LGH_deps);
if LGH_deps; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 4: ch11d (G.10) - Second-order secular-rate partials
%% Φ = G^10 * F_2^{**,(T)} = (3 mu^6 k_2^2 / 2) [Q(theta) P(e) eta/8 - R(theta) eta^4]
%% ∂F_2^{**,(T)}/∂L, ∂G, ∂H should be l, g, h-independent.
%% =======================================================
printf('--- Part 4: ch11d ∂F_2^{**,(T)}/∂(L,G,H) (secular, l,g,h-independent) ---\n');

F2ssT = c0T;  % same as c_0^{(T)}

% Closed-form partials from ch11d §2
dF2_dH_form = (9 * mu^6 * k2^2 * theta_expr / (2 * G^11)) * ...
              ((9*theta_expr^2 - 5) * P_expr * eta_expr / 16 - (3*theta_expr^2 - 1) * eta_expr^4);
dF2_dL_form = (3 * mu^6 * k2^2 * eta_expr^2 / (2 * G^11)) * ...
              (5 * Q_expr * (8 - 12*e_expr^2 - 3*e_expr^4) / 8 + 4 * R_expr * eta_expr^3);
Qp = (3*theta_expr/2) * (9*theta_expr^2 - 5);
Rp = 3*theta_expr * (3*theta_expr^2 - 1);
dF2_dG_inner = (eta_expr/8) * (-(10*Q_expr + theta_expr*Qp)*P_expr - 5*Q_expr*(8 - 12*e_expr^2 - 3*e_expr^4)) ...
             + eta_expr^4 * (6*R_expr + theta_expr*Rp);
dF2_dG_form = (3 * mu^6 * k2^2 / (2 * G^11)) * dF2_dG_inner;

% Direct diff
dF2_dL_direct = diff(F2ssT, L);
dF2_dG_direct = diff(F2ssT, G);
dF2_dH_direct = diff(F2ssT, H);

R_G10L = simplify(dF2_dL_direct - dF2_dL_form);
R_G10G = simplify(dF2_dG_direct - dF2_dG_form);
R_G10H = simplify(dF2_dH_direct - dF2_dH_form);

G10_correct = isequal(R_G10L, sym(0)) && isequal(R_G10G, sym(0)) && isequal(R_G10H, sym(0));
print_check('G.10 T-component closed-form correctness (all 3 partials)', ...
  'residuals all 0', sprintf('all zero: %d', G10_correct), G10_correct);
if G10_correct; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

% Secular property: l, g, h independence of the partials
dF2_dL_l = simplify(diff(dF2_dL_direct, l));
dF2_dL_g = simplify(diff(dF2_dL_direct, g));
dF2_dL_h = simplify(diff(dF2_dL_direct, h));
G10_secular = isequal(dF2_dL_l, sym(0)) && isequal(dF2_dL_g, sym(0)) && isequal(dF2_dL_h, sym(0));

print_check('G.10 T-component secular: ∂F_2^{**,(T)}/∂L has no l,g,h dependence', ...
  'all 3 residuals = 0', ...
  sprintf('all zero: %d', G10_secular), G10_secular);
if G10_secular; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 5: GSI (General Structural Identity) at alpha=3 for S_1^*
%% -(3/G) S_1^* - theta ∂S_1^*/∂H - (1/eta) ∂S_1^*/∂L = ∂S_1^*/∂G
%% =======================================================
printf('--- Part 5: GSI at alpha=3 for S_1^* ---\n');

dS1s_dL = diff(S1star, L);
dS1s_dG = diff(S1star, G);
dS1s_dH = diff(S1star, H);

GSI_LHS = dS1s_dG;
GSI_RHS = -(3/G) * S1star - theta_expr * dS1s_dH - (1/eta_expr) * dS1s_dL;
R_GSI = simplify(GSI_LHS - GSI_RHS);
p_GSI = isequal(R_GSI, sym(0));
if ~p_GSI
  R_GSI = simplify(expand(R_GSI));
  p_GSI = isequal(R_GSI, sym(0));
end

print_check('GSI at alpha=3: ∂S_1^*/∂G + 3 S_1^*/G + theta ∂S_1^*/∂H + (1/eta) ∂S_1^*/∂L = 0', ...
  'residual = 0', sprintf('%s', char(R_GSI)), p_GSI);
if p_GSI; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 6: GSI at alpha=6 for F_1^*
%% Should satisfy the GSI too.
%% =======================================================
printf('--- Part 6: GSI at alpha=6 for F_1^* ---\n');

dF1s_dL = diff(F1star, L);
dF1s_dG = diff(F1star, G);
dF1s_dH = diff(F1star, H);

GSI6_LHS = dF1s_dG;
GSI6_RHS = -(6/G) * F1star - theta_expr * dF1s_dH - (1/eta_expr) * dF1s_dL;
R_GSI6 = simplify(GSI6_LHS - GSI6_RHS);
p_GSI6 = isequal(R_GSI6, sym(0));

print_check('GSI at alpha=6: F_1^* GSI consistency', ...
  'residual = 0', sprintf('%s', char(R_GSI6)), p_GSI6);
if p_GSI6; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 7: GSI at alpha=10 for F_2^{**,(T)}
%% Should satisfy the GSI too.
%% =======================================================
printf('--- Part 7: GSI at alpha=10 for F_2^{**,(T)} ---\n');

GSI10_LHS = dF2_dG_direct;
GSI10_RHS = -(10/G) * F2ssT - theta_expr * dF2_dH_direct - (1/eta_expr) * dF2_dL_direct;
R_GSI10 = simplify(GSI10_LHS - GSI10_RHS);
p_GSI10 = isequal(R_GSI10, sym(0));
if ~p_GSI10
  R_GSI10 = simplify(expand(R_GSI10));
  p_GSI10 = isequal(R_GSI10, sym(0));
end

print_check('GSI at alpha=10: F_2^{**,(T)} GSI consistency', ...
  'residual = 0', sprintf('%s', char(R_GSI10)), p_GSI10);
if p_GSI10; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 8: Zero-mean check ⟨S_1^*⟩_g = 0
%% S_1^* has only sin(2g), sin(4g) harmonics, so g-integral over [0, 2pi] = 0.
%% =======================================================
printf('--- Part 8: Zero-mean ⟨S_1^*⟩_g = 0 ---\n');

S1s_avg_g = int(S1star, g, 0, 2*sym(pi)) / (2*sym(pi));
S1s_avg_g = simplify(S1s_avg_g);
p_zeromean = isequal(S1s_avg_g, sym(0));

print_check('⟨S_1^*⟩_g = 0 (long-period generator has zero g-mean)', ...
  'residual = 0', sprintf('%s', char(S1s_avg_g)), p_zeromean);
if p_zeromean; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 9: ⟨∂S_1^*/∂(L,G,H)⟩_g = 0 (g-average of momentum partials)
%% Since S_1^* has zero g-mean, differentiating in L, G, H commutes with ⟨·⟩_g.
%% So ⟨∂S_1^*/∂L⟩_g = ∂⟨S_1^*⟩_g/∂L = 0.
%% =======================================================
printf('--- Part 9: Zero g-mean of momentum partials of S_1^* ---\n');

dS1s_dL_avg_g = int(dS1s_dL, g, 0, 2*sym(pi)) / (2*sym(pi));
dS1s_dG_avg_g = int(dS1s_dG, g, 0, 2*sym(pi)) / (2*sym(pi));
dS1s_dH_avg_g = int(dS1s_dH, g, 0, 2*sym(pi)) / (2*sym(pi));

p_dL = isequal(simplify(dS1s_dL_avg_g), sym(0));
p_dG = isequal(simplify(dS1s_dG_avg_g), sym(0));
p_dH = isequal(simplify(dS1s_dH_avg_g), sym(0));
all_zmean_partials = p_dL && p_dG && p_dH;

print_check('⟨∂S_1^*/∂L⟩_g = ⟨∂S_1^*/∂G⟩_g = ⟨∂S_1^*/∂H⟩_g = 0', ...
  'all 3 residuals = 0', ...
  sprintf('all zero: %d (dL: %d, dG: %d, dH: %d)', all_zmean_partials, p_dL, p_dG, p_dH), ...
  all_zmean_partials);
if all_zmean_partials; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 10: F_2^{**,(T)} from g-averaging of F_{2p}^{(T)} + c_0^{(T)}
%% F_2^{*,(T)} = c_0^{(T)} + c_2^{(T)} cos(2g) + c_4^{(T)} cos(4g)
%% ⟨F_2^{*,(T)}⟩_g = c_0^{(T)} = F_2^{**,(T)}
%% =======================================================
printf('--- Part 10: F_2^{**} = ⟨F_2^*⟩_g consistency ---\n');

% c_2^{(T)} and c_4^{(T)} from ch10d §3.1
c2T = 3 * mu^6 * k2^2 * A_expr * B_expr * e_expr^2 * (6 + e_expr^2) * eta_expr / (4 * G^10);
c4T = 3 * mu^6 * k2^2 * B_expr^2 * e_expr^4 * eta_expr / (64 * G^10);

F2star_T = c0T + c2T * cos(2*g) + c4T * cos(4*g);
F2star_T_avg_g = simplify(int(F2star_T, g, 0, 2*sym(pi)) / (2*sym(pi)));
R_P10 = simplify(F2star_T_avg_g - c0T);
p_P10 = isequal(R_P10, sym(0));

print_check('F_2^{**,(T)} = ⟨c_0^{(T)} + c_2^{(T)} cos 2g + c_4^{(T)} cos 4g⟩_g = c_0^{(T)}', ...
  'residual = 0', sprintf('%s', char(R_P10)), p_P10);
if p_P10; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Part 11: Homological equation consistency
%% (∂F_1^*/∂G) (∂S_1^*/∂g) = -F_{2p} (ch11a G.1)
%% where F_{2p} = c_2^{(T)} cos(2g) + c_4^{(T)} cos(4g) (ch10d F.13, T-component)
%% =======================================================
printf('--- Part 11: Homological equation G.1 (T-component) ---\n');

F2p_T = c2T * cos(2*g) + c4T * cos(4*g);
dF1s_dG_local = diff(F1star, G);
dS1s_dg = diff(S1star, g);
LHS_G1 = dF1s_dG_local * dS1s_dg;
RHS_G1 = -F2p_T;
R_G1 = simplify(LHS_G1 - RHS_G1);
p_G1 = isequal(R_G1, sym(0));
if ~p_G1
  R_G1 = simplify(expand(R_G1));
  p_G1 = isequal(R_G1, sym(0));
end

print_check('G.1 homological: (∂F_1^*/∂G)(∂S_1^*/∂g) + F_{2p} = 0', ...
  'residual = 0', sprintf('%s', char(R_G1)), p_G1);
if p_G1; n_pass = n_pass + 1; else; n_fail = n_fail + 1; end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Variable-dependence audit: ALL %d checks PASSED.\n', n_pass);
  printf('All formulas respect their claimed independence assumptions.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
