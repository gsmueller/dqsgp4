% derive_deltaq_all_symbolic.m
%
% Full symbolic derivation of delta_q_1, delta_q_2, delta_q_3
% = D(∂(S_1 + S_1^*)/∂L), D(∂(S_1 + S_1^*)/∂G), D(∂(S_1 + S_1^*)/∂H).
%
% Replaces the "RESOLVED (structural)" placeholders in ch12 H.4, H.5, H.6
% with genuine symbolic closed forms derived via P25+ D-identity library
% + file 08 base identities.
%
% Method (same for all three):
%   1. Build ∂S_1/∂X and ∂S_1^*/∂X from the project's closed forms
%      (ch09d E.5 table + ch11b G.6/G.7/G.8).
%   2. Apply D symbolically with (θ, e, η, κ, f, E, l, g) as INDEPENDENT
%      symbols (avoid SymPy's issue with sqrt differentiation).
%   3. Run three audits on each: μ-exponent = 2, k_2-exponent = 1, ℳ_4.

pkg load symbolic;

printf('=========================================================\n');
printf('Full symbolic derivation of delta_q_1, delta_q_2, delta_q_3\n');
printf('  H.4, H.5, H.6 closed forms — replaces structural placeholders\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function mark_check(label, pass_flag, n_pass_in, n_fail_in)
  printf('  %s: %s\n', label, ternary(pass_flag, 'PASS', 'FAIL'));
  if pass_flag
    n_pass_in = n_pass_in + 1;
  else
    n_fail_in = n_fail_in + 1;
  end
  assignin('caller', 'n_pass', n_pass_in);
  assignin('caller', 'n_fail', n_fail_in);
endfunction

function s = ternary(cond, a, b)
  if cond; s = a; else; s = b; end
endfunction

% Independent symbols
syms L G H mu k2 l_s g_s h_s f_s E_s positive;
syms theta_s e_s eta_s kappa_s positive;

A_s = (3*theta_s^2 - 1)/2;
B_s = 3*(1 - theta_s^2)/2;

% Hansen X_0^{0,2}(e, eta) per ch06 B.0.7-7
X02 = (3*e_s^2 - 2 + 2*eta_s^3) / e_s^2;

% D-operator base identities (file 08 summary + file 08 line 823 for Dl)
Df   = 2*sin(f_s)/e_s;
Dg   = -2*sin(f_s)/e_s;
De   = -2*(e_s + cos(f_s));
Dtheta = sym(0);
Dl   = 2*sin(E_s)*(1 - e_s^3*cos(E_s))/(e_s * kappa_s);
DE   = 2*sin(E_s)/(e_s*kappa_s);
Deta = -e_s * De / eta_s;        % from eta^2 = 1 - e^2 => 2 eta Deta = -2e De
Dkappa = 4 - 2*kappa_s;          % from P25+.6
% D(1/L), D(1/G), D(1/e), etc.
% DL = -L(2a/r - 1). Using 1/kappa = a/r: 2a/r - 1 = 2/kappa - 1 = (2-kappa)/kappa.
% So DL = -L(2-kappa)/kappa.
DL_over_L = -(2 - kappa_s)/kappa_s;  % DL/L = -(2a/r - 1) = -(2/kappa - 1)
% D(1/G^n) = n/G^n (from DG^n = -n G^n Euler)
% D(1/L^n): chain through DL. D(L^{-n}) = -n L^{-n-1} DL = -n L^{-n-1} · L · DL_over_L = -n L^{-n} DL_over_L = n L^{-n} (2-kappa)/kappa

% Basis elements (ch09d Def E.5.1)
E1 = f_s - l_s + e_s*sin(f_s);
E2 = sin(2*(f_s + g_s));
E3 = sin(f_s + 2*g_s);
E4 = sin(3*f_s + 2*g_s);
E5 = sin(2*g_s);
E6 = sin(f_s)/kappa_s;
E7 = sin(f_s)^3;
E8 = sin(f_s)*(2 + e_s*cos(f_s))*cos(2*(f_s + g_s))/kappa_s;

% D(E_k):
DE1 = Df - Dl - 2*e_s*sin(f_s);             % P25+.3
DE2 = 2*(Df + Dg)*cos(2*(f_s + g_s));       % = 0 since Df+Dg=0
DE3 = (Df + 2*Dg)*cos(f_s + 2*g_s);         % P25+.1 j=1, k=2
DE4 = (3*Df + 2*Dg)*cos(3*f_s + 2*g_s);
DE5 = 2*Dg*cos(2*g_s);                       % P25+.1 j=0, k=2
% DE6 = D(sin f/kappa) = (cos f · Df)/kappa - sin f · Dkappa/kappa^2
DE6 = cos(f_s)*Df/kappa_s - sin(f_s)*Dkappa/kappa_s^2;
% DE7 = D(sin^3 f) = 3 sin^2 f · cos f · Df
DE7 = 3*sin(f_s)^2*cos(f_s)*Df;
% DE8: product of sin f, (2 + e cos f), cos(2(f+g)), 1/kappa — apply product rule
%  Let A1 = sin f, A2 = 2 + e cos f, A3 = cos 2(f+g), A4 = 1/kappa
%  D(A1 A2 A3 A4) = DA1·A2 A3 A4 + A1·DA2·A3 A4 + A1 A2·DA3·A4 + A1 A2 A3·DA4
DA1 = cos(f_s)*Df;
DA2 = De*cos(f_s) + e_s*(-sin(f_s))*Df;
DA3 = -2*(Df + Dg)*sin(2*(f_s + g_s));      % = 0 since Df+Dg=0
DA4 = -Dkappa/kappa_s^2;
A1 = sin(f_s); A2 = 2 + e_s*cos(f_s); A3 = cos(2*(f_s + g_s)); A4 = 1/kappa_s;
DE8 = DA1*A2*A3*A4 + A1*DA2*A3*A4 + A1*A2*DA3*A4 + A1*A2*A3*DA4;

printf('=========================================================\n');
printf('PART I: delta_q_1 = D(dS_1/dL + dS_1^*/dL) (H.4)\n');
printf('=========================================================\n\n');

% dS_1/dL from ch09d E.5 Part (b), native prefactor mu^2 k_2 / G^3
% Coefficients (column 2): c_1 = c_2 = 0; c_3..c_8 = ...
c3_L = B_s*eta_s^2/(2*L*e_s);
c4_L = B_s*eta_s^2/(6*L*e_s);
c5_L = B_s*eta_s^2*(eta_s + 2)/(3*L*(1 + eta_s)^2);
c6_L = 3*A_s*eta_s^2/(L*e_s);
c7_L = -A_s*e_s/L;
c8_L = B_s*eta_s^2/(L*e_s);

dS1_dL = (mu^2*k2/G^3)*(c3_L*E3 + c4_L*E4 + c5_L*E5 + c6_L*E6 + c7_L*E7 + c8_L*E8);

% D(dS_1/dL): use D applied to the whole expression.
% Treatment: L appears in the coefficients via the 1/L prefactor, and also in
% e, eta as e=sqrt(1-G^2/L^2), eta=G/L. But we keep (e, eta) INDEPENDENT here,
% so diff(..., L) will only pick up the 1/L in coefficients.
% For the D operator, we must combine: D acts through all variables.
%
% Simpler: observe dS_1/dL = (mu^2 k_2/G^3) · [sum], where [sum] has 1/L in coeffs.
% D((mu^2 k_2/G^3) [sum]) = (3 mu^2 k_2/G^3) [sum] + (mu^2 k_2/G^3) D([sum])
% D([sum]) = sum of D(c_k E_k) = sum of (D(c_k) E_k + c_k D(E_k)).
%
% D(c_k) = chain through (e, eta, L):
%   c_k = f_k(B, eta, e) · (1/L) typically, so Dc_k = Df_k · (1/L) + f_k · D(1/L)
%   where Df_k = (∂f_k/∂e)·De + (∂f_k/∂eta)·Deta (B has Dtheta-dependence which is 0).
%
% To compute D(c_k) fully symbolically with L-dependence: treat c_k as a function
% of (e_s, eta_s, L). Use symbolic partial derivatives.

% D acting on coefficient c(e, eta, L):
% D(c) = ∂c/∂e · De + ∂c/∂eta · Deta + ∂c/∂L · DL
% where DL = L · DL_over_L = -L(2-kappa)/kappa, so ∂c/∂L · DL.
DL_sym = L * DL_over_L;

function Dc = apply_D_to_coef(c, e_s, eta_s, L, De, Deta, DL_sym)
  Dc = diff(c, e_s)*De + diff(c, eta_s)*Deta + diff(c, L)*DL_sym;
endfunction

Dc3_L = apply_D_to_coef(c3_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc4_L = apply_D_to_coef(c4_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc5_L = apply_D_to_coef(c5_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc6_L = apply_D_to_coef(c6_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc7_L = apply_D_to_coef(c7_L, e_s, eta_s, L, De, Deta, DL_sym);
Dc8_L = apply_D_to_coef(c8_L, e_s, eta_s, L, De, Deta, DL_sym);

% D(dS_1/dL) = D(mu^2 k_2/G^3) · inner + (mu^2 k_2/G^3) · D(inner)
%   D(mu^2 k_2/G^3) = 3 mu^2 k_2/G^3 (Euler)
inner_dS1dL = c3_L*E3 + c4_L*E4 + c5_L*E5 + c6_L*E6 + c7_L*E7 + c8_L*E8;
D_inner_dS1dL = (Dc3_L*E3 + c3_L*DE3) + (Dc4_L*E4 + c4_L*DE4) ...
              + (Dc5_L*E5 + c5_L*DE5) + (Dc6_L*E6 + c6_L*DE6) ...
              + (Dc7_L*E7 + c7_L*DE7) + (Dc8_L*E8 + c8_L*DE8);

D_dS1_dL = (3*mu^2*k2/G^3)*inner_dS1dL + (mu^2*k2/G^3)*D_inner_dS1dL;

% dS_1^*/dL from ch11b G.6: G^{-4} * (eta^3/e) * ∂Phi/∂e
Phi_S1star = mu^2*k2/((5*theta_s^2 - 1)*eta_s^2) * ...
             (A_s*B_s*e_s^2*(6 + e_s^2)/4*sin(2*g_s) + B_s^2*e_s^4/128*sin(4*g_s));
dPhi_de = diff(Phi_S1star, e_s);
dS1star_dL = (eta_s^3/e_s)*dPhi_de / G^4;

% D(dS_1^*/dL). Use same structure. Note G^{-4} prefactor:
% D(G^{-4}·X) = (4/G^4)·X + (1/G^4)·D(X), where X = (eta^3/e)·∂Phi/∂e.
X_L = (eta_s^3/e_s)*dPhi_de;
% D(X_L) = ∂X_L/∂e · De + ∂X_L/∂eta · Deta + ∂X_L/∂theta · Dtheta + ∂X_L/∂g · Dg
%                                                                                   (Dtheta=0)
D_X_L = diff(X_L, e_s)*De + diff(X_L, eta_s)*Deta + diff(X_L, g_s)*Dg;
D_dS1star_dL = (4/G^4)*X_L + (1/G^4)*D_X_L;

delta_q_1 = D_dS1_dL + D_dS1star_dL;
printf('  delta_q_1 symbolic form assembled.\n');

% Dimensional audits
logderiv_mu = simplify(mu * diff(delta_q_1, mu) / delta_q_1);
mu_exp = double(subs(logderiv_mu, mu, sym(1)));
logderiv_k2 = simplify(k2 * diff(delta_q_1, k2) / delta_q_1);
k2_exp = double(subs(logderiv_k2, k2, sym(1)));
printf('  mu-exponent: %d (expected 2); k_2-exponent: %d (expected 1)\n', round(mu_exp), round(k2_exp));
mark_check('Task H.4 dim-audit', round(mu_exp) == 2 && round(k2_exp) == 1, n_pass, n_fail);

% M_4 scaling: verify diff(G^4 · delta_q_1, G) = 0 at fixed (theta, e, eta, ...)
G4_dq1 = simplify(G^4 * delta_q_1);
d_G4_dq1_dG = simplify(diff(G4_dq1, G));
mark_check('Task H.4 M_4 (diff G^4·delta_q_1, G) = 0', isequal(d_G4_dq1_dG, sym(0)), n_pass, n_fail);

% L-dependence test: G^4·delta_q_1 should depend on L only via DL_over_L (which involves
% kappa). With kappa treated as independent, L should not appear bare.
% Actually L appears in DL_sym = L·DL_over_L, and in coefficients via 1/L.
% So G^4·delta_q_1 has 1/L and L from DL combining. The "bare" L-dependence
% AFTER full simplification should vanish ONLY at the level of the coefficient
% structure. For generic kappa, we expect residual L-dependence through DL.
%
% Here the test passes if we check L-dependence AFTER substituting the DL identity
% fully, which gives back a function of (theta, e, eta, kappa, f, E, l, g).
% Since the δq_j is supposed to be a function solely of these dimensionless variables
% (and μ, k_2), no L should appear at all in G^4 · δq_1.

d_G4_dq1_dL = simplify(diff(G4_dq1, L));
is_L_indep = isequal(d_G4_dq1_dL, sym(0));
if ~is_L_indep
  d_G4_dq1_dL = simplify(expand(d_G4_dq1_dL));
  is_L_indep = isequal(d_G4_dq1_dL, sym(0));
end
mark_check('Task H.4 L-indep (after D substitution)', is_L_indep, n_pass, n_fail);

printf('\n');

printf('=========================================================\n');
printf('PART II: delta_q_2 = D(dS_1/dG + dS_1^*/dG) (H.5)\n');
printf('=========================================================\n\n');

% ∂S_1/∂G via GSI at alpha=3 (ch09b Theorem E.3):
%   ∂S_1/∂G = -(3/G) S_1 - theta · ∂S_1/∂H - (1/eta) · ∂S_1/∂L
%
% Since Dtheta = 0, we apply D directly:
% D(∂S_1/∂G) = -3·D(S_1/G) - theta·D(∂S_1/∂H) - D((1/eta)·∂S_1/∂L)

% S_1 closed form (ch07c C.4 / ch09d E.5 col 1):
c1_S1 = A_s;
c2_S1 = B_s/2;
c3_S1 = e_s*B_s/2;
c4_S1 = e_s*B_s/6;
c5_S1 = B_s*X02/6;
S1 = (mu^2*k2/G^3)*(c1_S1*E1 + c2_S1*E2 + c3_S1*E3 + c4_S1*E4 + c5_S1*E5);

% D(S_1):
Dc1_S1 = apply_D_to_coef(c1_S1, e_s, eta_s, L, De, Deta, DL_sym);  % = 0 since A only depends on theta
Dc2_S1 = apply_D_to_coef(c2_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc3_S1 = apply_D_to_coef(c3_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc4_S1 = apply_D_to_coef(c4_S1, e_s, eta_s, L, De, Deta, DL_sym);
Dc5_S1 = apply_D_to_coef(c5_S1, e_s, eta_s, L, De, Deta, DL_sym);

inner_S1 = c1_S1*E1 + c2_S1*E2 + c3_S1*E3 + c4_S1*E4 + c5_S1*E5;
D_inner_S1 = (Dc1_S1*E1 + c1_S1*DE1) + (Dc2_S1*E2 + c2_S1*DE2) ...
           + (Dc3_S1*E3 + c3_S1*DE3) + (Dc4_S1*E4 + c4_S1*DE4) ...
           + (Dc5_S1*E5 + c5_S1*DE5);

% D(S_1) = D(mu^2 k_2/G^3) · inner + (mu^2 k_2/G^3) · D(inner) = (3/G^3) · inner + (1/G^3) · D(inner)
% where D(1/G^3) = 3/G^3 (Euler)
D_S1 = (3*mu^2*k2/G^3)*inner_S1 + (mu^2*k2/G^3)*D_inner_S1;

% D(S_1/G) = D(S_1)/G + S_1·D(1/G) = D(S_1)/G + S_1/G
D_S1_over_G = D_S1/G + S1/G;

% ∂S_1/∂H (column 3) and its D
c1_H = 3*theta_s;
c2_H = -3*theta_s/2;
c3_H = -3*e_s*theta_s/2;
c4_H = -e_s*theta_s/2;
c5_H = -theta_s*X02/2;
Dc1_H = apply_D_to_coef(c1_H, e_s, eta_s, L, De, Deta, DL_sym);  % = 0
Dc2_H = apply_D_to_coef(c2_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc3_H = apply_D_to_coef(c3_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc4_H = apply_D_to_coef(c4_H, e_s, eta_s, L, De, Deta, DL_sym);
Dc5_H = apply_D_to_coef(c5_H, e_s, eta_s, L, De, Deta, DL_sym);
inner_dS1dH = c1_H*E1 + c2_H*E2 + c3_H*E3 + c4_H*E4 + c5_H*E5;
D_inner_dS1dH = (Dc1_H*E1 + c1_H*DE1) + (Dc2_H*E2 + c2_H*DE2) ...
              + (Dc3_H*E3 + c3_H*DE3) + (Dc4_H*E4 + c4_H*DE4) ...
              + (Dc5_H*E5 + c5_H*DE5);
D_dS1_dH = (4*mu^2*k2/G^4)*inner_dS1dH + (mu^2*k2/G^4)*D_inner_dS1dH;

% D((1/eta) · ∂S_1/∂L) = D(1/eta) · ∂S_1/∂L + (1/eta) · D(∂S_1/∂L)
dS1_dL_val = dS1_dL;  % already computed in Part I
D_inv_eta = e_s*De/eta_s^3;   % D(1/eta) = e De/eta^3 (from P25+.4 at n=1 but with eta not 1/eta)
% Actually D(eta^{-1}) = -eta^{-2} Deta = -(1/eta^2) · (-e De/eta) = e De/eta^3
D_1_over_eta_times_dS1dL = D_inv_eta * dS1_dL_val + (1/eta_s)*D_dS1_dL;

% D(∂S_1/∂G) via GSI:
D_dS1_dG = -3*D_S1_over_G - theta_s*D_dS1_dH - D_1_over_eta_times_dS1dL;

% ∂S_1^*/∂G via ch11b G.7: G^{-4} [-3 Phi - theta · ∂Phi/∂theta - (eta^2/e) · ∂Phi/∂e]
dPhi_dtheta = diff(Phi_S1star, theta_s);
dS1star_dG_inner = -3*Phi_S1star - theta_s*dPhi_dtheta - (eta_s^2/e_s)*dPhi_de;
dS1star_dG = dS1star_dG_inner / G^4;

% D(∂S_1^*/∂G):
% Apply D to G^{-4} · X where X = -3 Phi - theta ∂Phi/∂theta - (eta^2/e) ∂Phi/∂e
X_G_star = dS1star_dG_inner;
D_X_G_star = diff(X_G_star, e_s)*De + diff(X_G_star, eta_s)*Deta + diff(X_G_star, g_s)*Dg;
% (theta-dependence has Dtheta=0 so no theta term)
D_dS1star_dG = (4/G^4)*X_G_star + (1/G^4)*D_X_G_star;

delta_q_2 = D_dS1_dG + D_dS1star_dG;
printf('  delta_q_2 symbolic form assembled.\n');

logderiv_mu = simplify(mu * diff(delta_q_2, mu) / delta_q_2);
mu_exp = double(subs(logderiv_mu, mu, sym(1)));
logderiv_k2 = simplify(k2 * diff(delta_q_2, k2) / delta_q_2);
k2_exp = double(subs(logderiv_k2, k2, sym(1)));
printf('  mu-exponent: %d (expected 2); k_2-exponent: %d (expected 1)\n', round(mu_exp), round(k2_exp));
mark_check('Task H.5 dim-audit', round(mu_exp) == 2 && round(k2_exp) == 1, n_pass, n_fail);

G4_dq2 = simplify(G^4 * delta_q_2);
d_G4_dq2_dG = simplify(diff(G4_dq2, G));
if ~isequal(d_G4_dq2_dG, sym(0)); d_G4_dq2_dG = simplify(expand(d_G4_dq2_dG)); end
mark_check('Task H.5 M_4 (diff G^4·delta_q_2, G) = 0', isequal(d_G4_dq2_dG, sym(0)), n_pass, n_fail);

d_G4_dq2_dL = simplify(diff(G4_dq2, L));
if ~isequal(d_G4_dq2_dL, sym(0)); d_G4_dq2_dL = simplify(expand(d_G4_dq2_dL)); end
mark_check('Task H.5 L-indep', isequal(d_G4_dq2_dL, sym(0)), n_pass, n_fail);

printf('\n');

printf('=========================================================\n');
printf('PART III: delta_q_3 = D(dS_1/dH + dS_1^*/dH) (H.6)\n');
printf('  [cross-check of earlier derive_deltaq3_symbolic.m]\n');
printf('=========================================================\n\n');

% ∂S_1^*/∂H via ch11b G.8: G^{-4} · ∂Phi/∂theta
dS1star_dH = dPhi_dtheta / G^4;
X_H_star = dPhi_dtheta;
D_X_H_star = diff(X_H_star, e_s)*De + diff(X_H_star, eta_s)*Deta + diff(X_H_star, g_s)*Dg;
D_dS1star_dH = (4/G^4)*X_H_star + (1/G^4)*D_X_H_star;

delta_q_3 = D_dS1_dH + D_dS1star_dH;
printf('  delta_q_3 symbolic form assembled.\n');

logderiv_mu = simplify(mu * diff(delta_q_3, mu) / delta_q_3);
mu_exp = double(subs(logderiv_mu, mu, sym(1)));
logderiv_k2 = simplify(k2 * diff(delta_q_3, k2) / delta_q_3);
k2_exp = double(subs(logderiv_k2, k2, sym(1)));
printf('  mu-exponent: %d (expected 2); k_2-exponent: %d (expected 1)\n', round(mu_exp), round(k2_exp));
mark_check('Task H.6 dim-audit', round(mu_exp) == 2 && round(k2_exp) == 1, n_pass, n_fail);

G4_dq3 = simplify(G^4 * delta_q_3);
d_G4_dq3_dG = simplify(diff(G4_dq3, G));
if ~isequal(d_G4_dq3_dG, sym(0)); d_G4_dq3_dG = simplify(expand(d_G4_dq3_dG)); end
mark_check('Task H.6 M_4 (diff G^4·delta_q_3, G) = 0', isequal(d_G4_dq3_dG, sym(0)), n_pass, n_fail);

d_G4_dq3_dL = simplify(diff(G4_dq3, L));
if ~isequal(d_G4_dq3_dL, sym(0)); d_G4_dq3_dL = simplify(expand(d_G4_dq3_dL)); end
mark_check('Task H.6 L-indep', isequal(d_G4_dq3_dL, sym(0)), n_pass, n_fail);

%% Summary
printf('\n=========================================================\n');
if n_fail == 0
  printf('All three delta_q_j derived symbolically: %d/%d checks PASS.\n', n_pass, n_pass + n_fail);
  printf('H.4, H.5, H.6 ALL have genuine symbolic closed forms now.\n');
  printf('Dimensional audits: mu^2 k_2 prefactor + M_4 class verified symbolically.\n');
else
  printf('PARTIAL: %d/%d pass; %d fails.\n', n_pass, n_pass + n_fail, n_fail);
end
printf('=========================================================\n');
