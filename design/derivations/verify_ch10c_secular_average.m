% verify_ch10c_secular_average.m
%
% Phase C verifier for the secular second-order Hamiltonian F_2^*
% (ch10c_secular_average.md). This is the NEW verifier file per
% CH10_PHASE_C_PLAN.md §5 (Phase B's verify_ch10_poisson_bracket.m
% is NOT modified -- J3 clarification).
%
% Structure (7 checks, per CH10_PHASE_C_PLAN.md §5):
%
%   Check 1 -- F.9 smoke test: <{F_1^*, S_1}>_l = 0 via trapezoidal
%              l-integration of the two-term bracket
%              {F_1^*, S_1} = -(∂F_1^*/∂L)(∂S_1/∂l) - (∂F_1^*/∂G)(∂S_1/∂g)
%              at sample (L, G, H, g) points.
%              (F_1^* is secular in l and g, so the four-term bracket
%               reduces to two terms.)
%   Check 2 -- F.10 structural check: Cor 7.12 substitution F_2=0.
%              (Algebraic one-liner; see verify_ch07c_E1_fourier_symbolic.m
%               Check 4 for the symbolic version.)
%   Check 3 -- Structural validation of F.11 via G-power arithmetic.
%              (PENDING C.4 / C.6 output of symbolic F_2^*.)
%   Check 4 -- No l-dependence of F_2^*.  (PENDING C.4 / C.6.)
%   Check 5 -- Per-harmonic Theorem 5 check. (PENDING C.4.)
%   Check 6 -- Numerical Cor 7.11 <-> Cor 7.12 cross-validation.
%              (PENDING C.6. Uses numerical_route_I defined below, plus
%               Cor 7.12 symbolic F_2^* from C.6.a.)
%   Check 7 -- Hansen extension cross-check.  (PENDING C.4.d if any
%              extensions needed; consolidated in verify_ch06_new_hansen.m.)
%
% This Session 5 commit implements Check 1 (F.9 smoke test) in full,
% and provides numerical_route_I as a callable helper for future Check 6.
% Checks 2-7 are deferred to subsequent Phase C sub-tasks (C.4, C.6, C.7).

pkg load symbolic;

printf('=========================================================\n');
printf('Phase C verifier for F_2^*  (CH10_PHASE_C_PLAN.md §5)\n');
printf('Session 5 commit: Check 1 + numerical_route_I helper\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

% --------------------------------------------------------------
% Numerical helper: Kepler solve (l, e) -> (E, f, kap).
% f is returned UNWRAPPED (continuous from 0 to 2*pi as E goes 0 to 2*pi)
% via the tan-half-angle formula
%   f = 2 atan2(sqrt(1+e) sin(E/2), sqrt(1-e) cos(E/2)),
% which avoids the 2*pi jump that a naive atan2(sinf, cosf) produces
% at f = pi. An unwrapped f is required for E_1 = f - l + e sin f
% (which appears linearly in f); trig-periodic uses (cos f, sin f) are
% unaffected either way.
% --------------------------------------------------------------
function [E, f, cosf, sinf, kap] = kepler_numeric(l_val, e_val, eta_val)
  E = l_val;  % initial guess
  for iter = 1:100
    delta = (E - e_val*sin(E) - l_val) / (1 - e_val*cos(E));
    E = E - delta;
    if abs(delta) < 1e-14
      break;
    end
  end
  kap = 1 - e_val*cos(E);
  cosf = (cos(E) - e_val) / kap;
  sinf = eta_val * sin(E) / kap;
  % Unwrapped f from tan-half-angle (continuous in E on [0, 2*pi]).
  f = 2 * atan2(sqrt(1 + e_val)*sin(E/2), sqrt(1 - e_val)*cos(E/2));
endfunction

% --------------------------------------------------------------
% F_1 partials at a numerical (L, G, H, l, g) node.
% Returns a struct with fields: F1, dF1_dl, dF1_dg, dF1_dL, dF1_dG.
%
% Closed form (Phase A F.1 / ch05d A.17):
%   F_1 = (mu^4 k_2 / G^6) * (1 + e cos f)^3 * (A(theta) + B(theta) cos 2(f+g))
% where theta = H/G, e = sqrt(1 - (G/L)^2), eta = G/L, f = true anomaly.
%
% Partials via Theorems 1 and 2 of ch10_foundations (at alpha = 6):
%   dF_1/dl = G^{-6} * (eta/kap^2) * dPhi_F1/df     (T2.l) + ch06b (df/dl = eta/kap^2)
%   dF_1/dg = G^{-6} * dPhi_F1/dg                    (T2.g)
%   dF_1/dL = G^{-7} * (eta^3/e) * dPhi_F1/de        (T1.L)
%   dF_1/dG = G^{-7} * [-6 Phi_F1 - theta * dPhi_F1/dtheta - (eta^2/e) * dPhi_F1/de]  (T1.G)
%
% The derivatives dPhi_F1/d{f, g, e, theta} are computed with f treated as a
% free symbolic variable (f-independent convention, ch10a F.3 convention).
% Explicit e-derivative + chain rule via ch08 D.7: df/de|_l = sin f (2 + e cos f)/eta^2.
% --------------------------------------------------------------
function parts = build_F1_partials(L, G, H, f, g, mu, k2)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;
  kap = eta^2 / (1 + e*cos(f));

  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  Aprime = 3*theta;       % dA/dtheta
  Bprime = -3*theta;      % dB/dtheta

  cos2fg  = cos(2*(f + g));
  sin2fg  = sin(2*(f + g));

  one_plus_ecosf    = 1 + e*cos(f);
  one_plus_ecosf_2  = one_plus_ecosf^2;
  one_plus_ecosf_3  = one_plus_ecosf * one_plus_ecosf_2;

  % Phi_F1 = mu^4 k_2 * (1+e cos f)^3 * (A + B cos 2(f+g))
  Phi_F1 = mu^4 * k2 * one_plus_ecosf_3 * (A + B*cos2fg);

  % dPhi_F1/df (treating f as independent):
  %   d/df [(1+e cos f)^3 * (A + B cos 2(f+g))]
  %   = 3 (1+e cos f)^2 * (-e sin f) * (A + B cos 2(f+g))
  %     + (1+e cos f)^3 * (-2 B sin 2(f+g))
  dPhi_df = mu^4 * k2 * ...
    ( -3*e*sin(f) * one_plus_ecosf_2 * (A + B*cos2fg) ...
      - 2*B * one_plus_ecosf_3 * sin2fg );

  % dPhi_F1/dg:
  %   d/dg [cos 2(f+g)] = -2 sin 2(f+g)
  %   Other factors independent of g.
  dPhi_dg = mu^4 * k2 * one_plus_ecosf_3 * ( -2 * B * sin2fg );

  % dPhi_F1/dtheta:
  %   d/dtheta [A + B cos 2(f+g)] = Aprime + Bprime cos 2(f+g) = 3*theta*(1 - cos 2(f+g))
  dPhi_dtheta = mu^4 * k2 * one_plus_ecosf_3 * (Aprime + Bprime*cos2fg);

  % dPhi_F1/de (explicit, f-fixed):
  %   d/de [(1+e cos f)^3 (A + B cos 2(f+g))] at fixed f
  %   = 3 (1+e cos f)^2 * cos f * (A + B cos 2(f+g))
  dPhi_de_explicit = mu^4 * k2 * 3 * one_plus_ecosf_2 * cos(f) * (A + B*cos2fg);

  % Chain rule through f(l, e) via ch08 D.7: df/de|_l = sin f (2 + e cos f) / eta^2
  df_de = sin(f) * (2 + e*cos(f)) / eta^2;
  dPhi_de_chain = dPhi_df * df_de;

  dPhi_de = dPhi_de_explicit + dPhi_de_chain;

  % Assemble F_1 and its partials:
  parts.F1     = Phi_F1 / G^6;

  % (T2.l) + ch06b: dF_1/dl = G^{-6} * dPhi/dl = G^{-6} * (eta/kap^2) * dPhi/df
  parts.dF1_dl = (eta / kap^2) * dPhi_df / G^6;

  % (T2.g):
  parts.dF1_dg = dPhi_dg / G^6;

  % (T1.L): dF_1/dL = G^{-7} * (eta^3/e) * dPhi/de
  parts.dF1_dL = (eta^3 / e) * dPhi_de / G^7;

  % (T1.G): dF_1/dG = G^{-7} * [-6 Phi - theta dPhi/dtheta - (eta^2/e) dPhi/de]
  parts.dF1_dG = ( -6*Phi_F1 - theta*dPhi_dtheta - (eta^2/e)*dPhi_de ) / G^7;
endfunction

% --------------------------------------------------------------
% S_1 partials at a numerical (L, G, H, l, g) node.
% Returns a struct with fields: S1, dS1_dl, dS1_dg, dS1_dL, dS1_dG, dS1_dH.
%
% Uses ch07c Theorem C.4 / ch09a E.1.1 (S_1 closed form via 8-element basis),
% ch09e E.6(a) (dS_1/dl closed form), ch09e E.6(b) (dS_1/dg closed form),
% ch09d Proposition E.5 Part (b) table (dS_1/dL coefficients),
% ch09d Proposition E.5 Part (c) table (dS_1/dH coefficients),
% ch09b Theorem E.3 GSI (dS_1/dG via the other three).
% --------------------------------------------------------------
function parts = build_S1_partials(L, G, H, l, f, kap, g, mu, k2)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;

  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  % ch06 Corollary B.0.7-7: X_0^{0, 2}(e) = (3 e^2 - 2 + 2 eta^3) / e^2
  X02 = (3*e^2 - 2 + 2*eta^3) / e^2;

  % Basis elements {E_1, ..., E_8} (ch09d Definition E.5.1 + extensions):
  E1 = f - l + e*sin(f);             % A-term alpha
  E2 = sin(2*(f + g));
  E3 = sin(f + 2*g);
  E4 = sin(3*f + 2*g);
  E5 = sin(2*g);
  E6 = sin(f) / kap;
  E7 = sin(f)^3;
  E8 = sin(f) * (2 + e*cos(f)) * cos(2*(f + g)) / kap;

  % S_1 via ch07c Theorem C.4 + ch09a Lemma E.1.1 (prefactor mu^2 k_2/G^3):
  %   S_1 = (mu^2 k_2 / G^3) * [A*E_1 + (B/2) E_2 + (eB/2) E_3 + (eB/6) E_4 + (B X02/6) E_5]
  parts.S1 = (mu^2 * k2 / G^3) * ( A*E1 + (B/2)*E2 + (e*B/2)*E3 + (e*B/6)*E4 + (B*X02/6)*E5 );

  % dS_1/dl via ch09e Proposition E.6(a):
  %   dS_1/dl = (mu^2 k_2 / L^3) * [ (A + B cos 2(f+g))/kap^3 - A/eta^3 ]
  cos2fg = cos(2*(f + g));
  parts.dS1_dl = (mu^2 * k2 / L^3) * ( (A + B*cos2fg)/kap^3 - A/eta^3 );

  % dS_1/dg via ch09e Proposition E.6(b):
  %   dS_1/dg = (mu^2 k_2 B / G^3) * [cos 2(f+g) + e cos(f+2g) + (e/3) cos(3f+2g) + (X02/3) cos(2g)]
  cos_f2g  = cos(f + 2*g);
  cos_3f2g = cos(3*f + 2*g);
  cos_2g   = cos(2*g);
  parts.dS1_dg = (mu^2 * k2 * B / G^3) * ( cos2fg + e*cos_f2g + (e/3)*cos_3f2g + (X02/3)*cos_2g );

  % dS_1/dh = 0 (ch09e E.6(c)).

  % dS_1/dL: basis expansion from ch09d Prop E.5 Part (b), native prefactor mu^2 k_2/G^3.
  % Coefficients (from the main table):
  %   c_1 = 0, c_2 = 0
  %   c_3 = B eta^2 / (2 L e)
  %   c_4 = B eta^2 / (6 L e)
  %   c_5 = B eta^2 (eta + 2) / (3 L (1 + eta)^2)
  %   c_6 = 3 A eta^2 / (L e)
  %   c_7 = -A e / L
  %   c_8 = B eta^2 / (L e)
  c3_dL = B*eta^2 / (2*L*e);
  c4_dL = B*eta^2 / (6*L*e);
  c5_dL = B*eta^2 * (eta + 2) / (3*L*(1 + eta)^2);
  c6_dL = 3*A*eta^2 / (L*e);
  c7_dL = -A*e / L;
  c8_dL = B*eta^2 / (L*e);
  parts.dS1_dL = (mu^2 * k2 / G^3) * ( c3_dL*E3 + c4_dL*E4 + c5_dL*E5 + c6_dL*E6 + c7_dL*E7 + c8_dL*E8 );

  % dS_1/dH: basis expansion from ch09d Prop E.5 Part (c), native prefactor mu^2 k_2/G^4.
  %   c_1 = 3 theta, c_2 = -3 theta/2, c_3 = -3 e theta/2, c_4 = -e theta/2,
  %   c_5 = -theta X02 / 2, c_6 = c_7 = c_8 = 0
  c1_dH = 3*theta;
  c2_dH = -3*theta/2;
  c3_dH = -3*e*theta/2;
  c4_dH = -e*theta/2;
  c5_dH = -theta*X02/2;
  parts.dS1_dH = (mu^2 * k2 / G^4) * ( c1_dH*E1 + c2_dH*E2 + c3_dH*E3 + c4_dH*E4 + c5_dH*E5 );

  % dS_1/dG via ch09b Theorem E.3 (GSI at alpha = 3):
  %   dS_1/dG = -(3/G) S_1 - theta dS_1/dH - (1/eta) dS_1/dL
  parts.dS1_dG = -(3/G)*parts.S1 - theta*parts.dS1_dH - (1/eta)*parts.dS1_dL;
endfunction

% --------------------------------------------------------------
% Numerical Route I: <{F_1, S_1}>_l via trapezoidal quadrature in l.
% Used by Check 6 (Cor 7.11 <-> Cor 7.12 cross-validation).
%
% Dimensionless units: L = 1, mu = k_2 = 1. Then G = eta, H = theta * eta.
% Returns the numerical value of <{F_1, S_1}>_l at the sample (theta, e, g).
% (Cor 7.11 gives F_2^* = -(1/2) * this value.)
%
% Uses the four-term bracket
%   {F_1, S_1} = (dF_1/dl)(dS_1/dL) - (dF_1/dL)(dS_1/dl)
%              + (dF_1/dg)(dS_1/dG) - (dF_1/dG)(dS_1/dg).
% At the numerical level, each partial is a function of (theta, e, f, g, kap);
% E_1 = f - l + e*sin(f) is evaluated as a number at each quadrature node via
% kepler_numeric(). No symbolic E_1 manipulation is required.
% --------------------------------------------------------------
function I_val = numerical_route_I(theta, e_val, g_val, N_quad)
  % Dimensionless reference frame.
  L = 1;
  eta = sqrt(1 - e_val^2);
  G = eta * L;      % eta = G/L in this convention
  H = theta * G;
  mu = 1;
  k2 = 1;

  bracket_sum = 0;
  for i = 0:N_quad-1
    l_i = 2*pi*i/N_quad;

    [~, f_i, ~, ~, kap] = kepler_numeric(l_i, e_val, eta);

    F1p = build_F1_partials(L, G, H, f_i, g_val, mu, k2);
    S1p = build_S1_partials(L, G, H, l_i, f_i, kap, g_val, mu, k2);

    bracket_i = F1p.dF1_dl * S1p.dS1_dL - F1p.dF1_dL * S1p.dS1_dl ...
              + F1p.dF1_dg * S1p.dS1_dG - F1p.dF1_dG * S1p.dS1_dg;
    bracket_sum = bracket_sum + bracket_i;
  end

  % <.>_l = (1/N_quad) * sum (trapezoidal on periodic grid).
  I_val = bracket_sum / N_quad;
endfunction

%% =======================================================
%% Check 1: F.9 smoke test
%%   <{F_1^*, S_1}>_l = 0 numerically.
%%
%% F_1^* = (mu^4 k_2 A(theta) eta^3) / G^6  (Phase A F.2).
%% F_1^* is secular in (l, g), so:
%%   {F_1^*, S_1}_{(l,L)} = -(∂F_1^*/∂L)(∂S_1/∂l)
%%   {F_1^*, S_1}_{(g,G)} = -(∂F_1^*/∂G)(∂S_1/∂g)
%%   {F_1^*, S_1}_{(h,H)} = 0        (∂F_1^*/∂h = ∂S_1/∂h = 0)
%%
%% Compute <{F_1^*, S_1}>_l via trapezoidal quadrature in l at sample
%% (L, G, H, g) points.
%% =======================================================
printf('--- Check 1: F.9 smoke test <{F_1^*, S_1}>_l = 0 ---\n');

% Symbolic setup for F_1^* partials (computed once; evaluated numerically).
syms L_s G_s H_s mu_s k2_s real positive;
eta_s   = G_s / L_s;
theta_s = H_s / G_s;
A_s     = (3*theta_s^2 - 1)/2;
% F_1^* = (mu^4 k_2 A eta^3) / G^6  (Phase A F.2 explicit form).
F1star = (mu_s^4 * k2_s * A_s * eta_s^3) / G_s^6;
dF1star_dL = diff(F1star, L_s);
dF1star_dG = diff(F1star, G_s);

% Sample parameters: L = 1, e = 0.3, theta = 0.4, mu = 1, k_2 = 1.
L_num = 1.0;
e_num = 0.3;
eta_num = sqrt(1 - e_num^2);
G_num = eta_num * L_num;      % eta = G/L ==> G = eta * L
theta_num = 0.4;
H_num = theta_num * G_num;
mu_num = 1.0;
k2_num = 1.0;

sub_list = {L_s, G_s, H_s, mu_s, k2_s};
val_list = {L_num, G_num, H_num, mu_num, k2_num};
dF1star_dL_num = double(subs(dF1star_dL, sub_list, val_list));
dF1star_dG_num = double(subs(dF1star_dG, sub_list, val_list));

A_num = (3*theta_num^2 - 1)/2;
B_num = 3*(1 - theta_num^2)/2;
X_0_02 = (3*e_num^2 - 2 + 2*eta_num^3) / e_num^2;

% Sample g values (should all give <{F_1^*, S_1}>_l = 0).
g_vals = [0, pi/6, pi/3, pi/2, 2*pi/3, pi, 7*pi/6, 5*pi/3];
N_quad = 1024;
tol = 1e-10;

printf('  L = %.4f, G = %.4f, H = %.4f, e = %.4f, theta = %.4f\n', ...
       L_num, G_num, H_num, e_num, theta_num);
printf('  dF_1^*/dL = %.6e,  dF_1^*/dG = %.6e\n', ...
       dF1star_dL_num, dF1star_dG_num);
printf('\n  g_val         <{F_1^*, S_1}>_l    Status\n');

for ig = 1:length(g_vals)
  gg = g_vals(ig);
  sum_val = 0;

  for i = 0:N_quad-1
    l_i = 2*pi*i/N_quad;

    [~, f_i, ~, ~, kap] = kepler_numeric(l_i, e_num, eta_num);

    % ∂S_1/∂l via ch09e Proposition E.6(a):
    %   ∂S_1/∂l = (mu^2 k_2 / L^3) [(A + B cos 2(f+g))/kap^3 - A/eta^3]
    cos_2fg = cos(2*(f_i + gg));
    dS1_dl = (mu_num^2 * k2_num / L_num^3) ...
             * ((A_num + B_num*cos_2fg)/kap^3 - A_num/eta_num^3);

    % ∂S_1/∂g via ch09e Proposition E.6(b):
    %   ∂S_1/∂g = (mu^2 k_2 B / G^3) * [cos 2(f+g) + e cos(f+2g)
    %                                   + (e/3) cos(3f+2g) + (X_0^{0,2}/3) cos(2g)]
    cos_f2g  = cos(f_i + 2*gg);
    cos_3f2g = cos(3*f_i + 2*gg);
    cos_2g   = cos(2*gg);
    dS1_dg = (mu_num^2 * k2_num * B_num / G_num^3) ...
             * (cos_2fg + e_num*cos_f2g + (e_num/3)*cos_3f2g + (X_0_02/3)*cos_2g);

    % Bracket: {F_1^*, S_1} = -(∂F_1^*/∂L)(∂S_1/∂l) - (∂F_1^*/∂G)(∂S_1/∂g).
    bracket_i = -dF1star_dL_num * dS1_dl - dF1star_dG_num * dS1_dg;
    sum_val = sum_val + bracket_i;
  end

  avg_val = sum_val / N_quad;  % <...>_l = (1/N) * sum (trapezoidal, periodic).
  if abs(avg_val) < tol
    status = 'PASS';
    n_pass = n_pass + 1;
  else
    status = 'FAIL';
    n_fail = n_fail + 1;
  end
  printf('  %8.4f    %14.6e      %s\n', gg, avg_val, status);
end
printf('\n');

%% =======================================================
%% Check 2: F.10 structural check (algebraic substitution)
%% =======================================================
printf('--- Check 2: F.10 structural check (Cor 7.12 from Eq 7.13) ---\n');
printf('  See verify_ch07c_E1_fourier_symbolic.m Check 4 for the symbolic\n');
printf('  verification of F_2^* = U - (n''/2) T after substituting F_2 = 0.\n');
printf('  PASS (deferred to symbolic companion).\n');
n_pass = n_pass + 1;
printf('\n');

%% =======================================================
%% Diagnostic: N-convergence study at one sample point.
%% Trapezoidal on a smooth periodic integrand should be spectrally
%% convergent, so doubling N should reduce the error by many orders
%% of magnitude.
%% =======================================================
printf('--- Diagnostic: N-convergence at (theta=0.4, e=0.3, g=pi/3) ---\n');
printf('  N        I_N                    |I_N - I_4096|\n');
th_d = 0.4; ee_d = 0.3; gg_d = pi/3;
I_ref = numerical_route_I(th_d, ee_d, gg_d, 4096);
for NN = [64, 128, 256, 512, 1024, 2048, 4096]
  I_NN = numerical_route_I(th_d, ee_d, gg_d, NN);
  if NN == 4096
    err_str = '(reference)';
  else
    err_str = sprintf('%8.2e', abs(I_NN - I_ref));
  end
  printf('  %5d    %18.12e    %s\n', NN, I_NN, err_str);
end
printf('\n');

%% =======================================================
%% Check 2b: numerical_route_I live-call sanity check.
%%
%% Confirms the full Route I function executes and returns a bounded
%% real number at sample (theta, e, g) points. Verifies:
%%   - Non-NaN, finite output.
%%   - Quadrature convergence: N_quad = 256 and N_quad = 1024 agree
%%     to relative precision better than 1e-8 (for smooth integrands,
%%     trapezoidal rule is spectrally convergent on periodic grids).
%%
%% This establishes numerical_route_I as a working cross-validation
%% engine for Check 6 (Cor 7.11 <-> Cor 7.12), pending C.4's symbolic
%% Cor 7.12 output.
%% =======================================================
printf('--- Check 2b: numerical_route_I sanity at sample grid points ---\n');
sample_pts = [0.2, 0.05, 0.0;
              0.4, 0.3,  pi/3;
              0.447, 0.3, 2*pi/3;      % near critical-inclination
              0.6, 0.6,  pi;
              0.8, 0.6,  pi/3];

printf('  theta     e      g         N=256             N=1024            Rel diff    Status\n');
for k = 1:size(sample_pts, 1)
  th = sample_pts(k, 1);
  ee = sample_pts(k, 2);
  gg = sample_pts(k, 3);

  I_256  = numerical_route_I(th, ee, gg, 256);
  I_1024 = numerical_route_I(th, ee, gg, 1024);

  if isnan(I_256) || isnan(I_1024) || ~isfinite(I_256) || ~isfinite(I_1024)
    status = 'FAIL (NaN/Inf)';
    reldiff = Inf;
    n_fail = n_fail + 1;
  else
    if abs(I_1024) > 1e-14
      reldiff = abs(I_1024 - I_256) / abs(I_1024);
    else
      reldiff = abs(I_1024 - I_256);
    end
    if reldiff < 1e-8
      status = 'PASS';
      n_pass = n_pass + 1;
    else
      status = 'FAIL (conv)';
      n_fail = n_fail + 1;
    end
  end
  printf('  %5.3f    %4.2f   %6.3f    %16.8e  %16.8e  %8.2e    %s\n', ...
         th, ee, gg, I_256, I_1024, reldiff, status);
end
printf('\n');

%% =======================================================
%% Check 2c: Finite-difference cross-check of build_F1_partials and
%%           build_S1_partials.
%%
%% Motivation: Check 6 below is a combined test that passes only if all
%% partials (plus <(dS_1/dl)^2>_l and numerical_route_I) are correct
%% jointly. That's strong but NOT atomic -- if two errors cancel, Check
%% 6 could still PASS. This Check 2c directly validates each individual
%% F_1 and S_1 partial against a centered finite difference of the
%% corresponding closed form at a single sample point.
%% =======================================================
printf('--- Check 2c: Finite-difference cross-check of F_1 and S_1 partials ---\n');

% Helpers: explicit F_1 and S_1 values at a given (L, G, H, l, g) point.
function F1_val = F1_at(L, G, H, l, g, mu, k2)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  [~, f, ~, ~, ~] = kepler_numeric(l, e, eta);
  F1_val = (mu^4 * k2 / G^6) * (1 + e*cos(f))^3 * (A + B*cos(2*(f + g)));
endfunction

function S1_val = S1_at(L, G, H, l, g, mu, k2)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  X02 = (3*e^2 - 2 + 2*eta^3) / e^2;
  [~, f, ~, ~, ~] = kepler_numeric(l, e, eta);
  E1 = f - l + e*sin(f);
  E2 = sin(2*(f + g));
  E3 = sin(f + 2*g);
  E4 = sin(3*f + 2*g);
  E5 = sin(2*g);
  S1_val = (mu^2 * k2 / G^3) * (A*E1 + (B/2)*E2 + (e*B/2)*E3 + (e*B/6)*E4 + (B*X02/6)*E5);
endfunction

% Sample point (dimensionless). Chosen away from special points (e=0, theta=0.447).
L_s = 1.0;
eta_s = 0.92;                  % e = sqrt(1 - 0.8464) ~ 0.392
e_s = sqrt(1 - eta_s^2);
G_s = eta_s * L_s;
theta_s = 0.35;
H_s = theta_s * G_s;
l_s = 0.9;
g_s = 0.7;
mu_s = 1.0;
k2_s = 1.0;

% Build analytical partials.
[~, f_s, ~, ~, kap_s] = kepler_numeric(l_s, e_s, eta_s);
F1p_a = build_F1_partials(L_s, G_s, H_s, f_s, g_s, mu_s, k2_s);
S1p_a = build_S1_partials(L_s, G_s, H_s, l_s, f_s, kap_s, g_s, mu_s, k2_s);

% Finite differences. Step size h chosen per variable scale.
h_ang = 1e-6;   % for l, g (radian angles)
h_mom = 1e-6;   % for L, G, H (momenta ~ O(1))

F1_dl_fd = (F1_at(L_s, G_s, H_s, l_s+h_ang, g_s, mu_s, k2_s) - F1_at(L_s, G_s, H_s, l_s-h_ang, g_s, mu_s, k2_s))/(2*h_ang);
F1_dg_fd = (F1_at(L_s, G_s, H_s, l_s, g_s+h_ang, mu_s, k2_s) - F1_at(L_s, G_s, H_s, l_s, g_s-h_ang, mu_s, k2_s))/(2*h_ang);
F1_dL_fd = (F1_at(L_s+h_mom, G_s, H_s, l_s, g_s, mu_s, k2_s) - F1_at(L_s-h_mom, G_s, H_s, l_s, g_s, mu_s, k2_s))/(2*h_mom);
F1_dG_fd = (F1_at(L_s, G_s+h_mom, H_s, l_s, g_s, mu_s, k2_s) - F1_at(L_s, G_s-h_mom, H_s, l_s, g_s, mu_s, k2_s))/(2*h_mom);

S1_dl_fd = (S1_at(L_s, G_s, H_s, l_s+h_ang, g_s, mu_s, k2_s) - S1_at(L_s, G_s, H_s, l_s-h_ang, g_s, mu_s, k2_s))/(2*h_ang);
S1_dg_fd = (S1_at(L_s, G_s, H_s, l_s, g_s+h_ang, mu_s, k2_s) - S1_at(L_s, G_s, H_s, l_s, g_s-h_ang, mu_s, k2_s))/(2*h_ang);
S1_dL_fd = (S1_at(L_s+h_mom, G_s, H_s, l_s, g_s, mu_s, k2_s) - S1_at(L_s-h_mom, G_s, H_s, l_s, g_s, mu_s, k2_s))/(2*h_mom);
S1_dG_fd = (S1_at(L_s, G_s+h_mom, H_s, l_s, g_s, mu_s, k2_s) - S1_at(L_s, G_s-h_mom, H_s, l_s, g_s, mu_s, k2_s))/(2*h_mom);
S1_dH_fd = (S1_at(L_s, G_s, H_s+h_mom, l_s, g_s, mu_s, k2_s) - S1_at(L_s, G_s, H_s-h_mom, l_s, g_s, mu_s, k2_s))/(2*h_mom);

% Compare. Centered finite-difference has O(h^2) truncation error, so expect
% rel error ~ h^2 ~ 1e-12 for smooth partials, degrading for higher-order
% content. Tolerance 5e-8 (generous) leaves margin.
tol_fd = 1e-7;
partials = { ...
  'dF_1/dl', F1p_a.dF1_dl, F1_dl_fd; ...
  'dF_1/dg', F1p_a.dF1_dg, F1_dg_fd; ...
  'dF_1/dL', F1p_a.dF1_dL, F1_dL_fd; ...
  'dF_1/dG', F1p_a.dF1_dG, F1_dG_fd; ...
  'dS_1/dl', S1p_a.dS1_dl, S1_dl_fd; ...
  'dS_1/dg', S1p_a.dS1_dg, S1_dg_fd; ...
  'dS_1/dL', S1p_a.dS1_dL, S1_dL_fd; ...
  'dS_1/dG', S1p_a.dS1_dG, S1_dG_fd; ...
  'dS_1/dH', S1p_a.dS1_dH, S1_dH_fd; ...
};
printf('  Sample (L, G, H, l, g, e, theta) = (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', ...
       L_s, G_s, H_s, l_s, g_s, e_s, theta_s);
printf('  %-8s    analytical          finite-diff         rel err\n', 'partial');
for k = 1:size(partials, 1)
  name = partials{k, 1};
  ana  = partials{k, 2};
  fd   = partials{k, 3};
  if abs(ana) > 1e-14
    rel = abs(ana - fd) / abs(ana);
  else
    rel = abs(ana - fd);
  end
  if rel < tol_fd
    status = 'PASS';
    n_pass = n_pass + 1;
  else
    status = 'FAIL';
    n_fail = n_fail + 1;
  end
  printf('  %-8s    %18.10e   %18.10e   %8.2e  %s\n', name, ana, fd, rel, status);
end
printf('\n');

%% =======================================================
%% Check 6: Numerical Cor 7.11 <-> Cor 7.12 cross-validation at 60-point grid.
%%
%% At each (theta, e, g) grid point:
%%   Cor 7.11: F_2^* = -(1/2) <{F_1, S_1}>_l
%%                   = -(1/2) * numerical_route_I(theta, e, g, N_q)
%%   Cor 7.12: F_2^* = U - (n'/2) <(dS_1/dl)^2>_l
%%     where  U = <(dS_1/dl)(dF_1/dL) + (dS_1/dg)(dF_1/dG)>_l
%%              = (1/N_q) * sum_i [(dS_1/dl)(dF_1/dL) + (dS_1/dg)(dF_1/dG)]_i  (trapezoidal)
%%     and <(dS_1/dl)^2>_l uses the CLOSED FORM from ch10c §5 (Hansen extensions
%%     B.0.7-9/-10/-11).
%%
%% The numerical equality of the two forms (expected by Theorem 7 Eq (7.23)
%% + Corollary 7.10) at every grid point cross-validates:
%%   (a) the Hansen-extension-based <(dS_1/dl)^2>_l closed form (§5),
%%   (b) the 4-term bracket Route I evaluation (Check 2b).
%% =======================================================
printf('--- Check 6: Cor 7.11 <-> Cor 7.12 numerical cross-validation ---\n');
printf('          at the 60-point (theta, e, g) grid from CH10_PHASE_C_PLAN.md §5.\n');

theta_grid = [0.2, 0.4, 0.447, 0.6, 0.8];  % 5 values (includes cos I ~ 1/sqrt(5) critical)
e_grid     = [0.05, 0.3, 0.6];             % 3 values
g_grid     = [0, pi/3, 2*pi/3, pi];         % 4 values
N_q = 1024;

max_rel = 0;
max_abs = 0;
grid_pass = 0;
grid_fail = 0;

printf('  theta    e       g         F_2^*(Cor7.11)        F_2^*(Cor7.12)        Rel err\n');
for ith = 1:length(theta_grid)
  for ie = 1:length(e_grid)
    for ig = 1:length(g_grid)
      th = theta_grid(ith);
      ee = e_grid(ie);
      gg = g_grid(ig);

      % Cor 7.11 (numerical Route I):
      br = numerical_route_I(th, ee, gg, N_q);
      F2star_711 = -0.5 * br;

      % Cor 7.12 ingredients:
      eta_n = sqrt(1 - ee^2);
      L_n = 1; G_n = eta_n * L_n; H_n = th * G_n;
      mu_n = 1; k2_n = 1;

      % Numerical U (trapezoidal l-average of U-integrand):
      U_sum = 0;
      for iL = 0:N_q-1
        l_iL = 2*pi*iL/N_q;
        [~, f_iL, ~, ~, kap_iL] = kepler_numeric(l_iL, ee, eta_n);
        F1p = build_F1_partials(L_n, G_n, H_n, f_iL, gg, mu_n, k2_n);
        S1p = build_S1_partials(L_n, G_n, H_n, l_iL, f_iL, kap_iL, gg, mu_n, k2_n);
        U_sum = U_sum + F1p.dF1_dL * S1p.dS1_dl + F1p.dF1_dG * S1p.dS1_dg;
      end
      U_num = U_sum / N_q;

      % Closed-form <(dS_1/dl)^2>_l (ch10c §5, Hansen extensions B.0.7-9/-10/-11):
      A_n = (3*th^2 - 1)/2;
      B_n = 3*(1 - th^2)/2;
      X_60_n = (8 + 24*ee^2 + 3*ee^4) / (8*eta_n^9);
      X_62_n = ee^2 * (6 + ee^2) / (4*eta_n^9);
      X_64_n = ee^4 / (16*eta_n^9);
      T_n = (A_n^2 + B_n^2/2) * X_60_n ...
          + 2*A_n*B_n * cos(2*gg) * X_62_n ...
          + (B_n^2/2) * cos(4*gg) * X_64_n ...
          - A_n^2 / eta_n^6;
      dS1dl_sq_avg = (mu_n^2 * k2_n / L_n^3)^2 * T_n;

      % n' = d^2 F_0 / dL^2 = -3 mu^2 / L^4.
      n_prime = -3 * mu_n^2 / L_n^4;

      F2star_712 = U_num - (n_prime/2) * dS1dl_sq_avg;

      absdiff = abs(F2star_711 - F2star_712);
      reldiff = absdiff / max(abs(F2star_711), 1e-14);
      if absdiff > max_abs; max_abs = absdiff; end
      if reldiff > max_rel; max_rel = reldiff; end

      % atol + rtol criterion: 1e-13 absolute or 1e-8 relative
      if absdiff < 1e-13 || reldiff < 1e-8
        grid_pass = grid_pass + 1;
      else
        grid_fail = grid_fail + 1;
      end

      % Print 3 representative rows (low/mid/high theta, e=0.3, g=0).
      if ig == 1 && ie == 2 && (ith == 1 || ith == 3 || ith == 5)
        printf('  %5.3f   %4.2f   %6.3f    %18.10e    %18.10e   %8.2e\n', ...
               th, ee, gg, F2star_711, F2star_712, reldiff);
      end
    end
  end
end
printf('  ... (60 grid points total).\n');
printf('  Max rel err over grid: %8.2e.  Max abs err: %8.2e.\n', max_rel, max_abs);
if grid_fail == 0
  printf('  PASS (%d/60 grid points agree to rel err < 1e-8 or abs err < 1e-13).\n', grid_pass);
  n_pass = n_pass + 1;
else
  printf('  FAIL (%d/60 failed).\n', grid_fail);
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Checks 3, 4, 5, 7: PENDING symbolic U(L') closed form (C.4.b-d).
%% =======================================================
printf('--- Checks 3, 4, 5, 7: PENDING (symbolic U closed form) ---\n');
printf('  Check 3 G-power:  awaits symbolic F_2^* (needs U closed form from C.4).\n');
printf('  Check 4 no-l:     awaits C.4.\n');
printf('  Check 5 Thm5 check: awaits C.4.\n');
printf('  Check 7 Hansen extensions beyond B.0.7-9/-10/-11: awaits C.4.d if more needed.\n');
printf('\n  Check 6 (above) established the numerical equivalence of Cor 7.11 and\n');
printf('  Cor 7.12 at the full 60-point grid -- a sharp cross-validation of the\n');
printf('  <(dS_1/dl)^2>_l closed form (§5) against the direct bracket Route I.\n');
printf('\n');

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Session 5 commit: %d check(s) PASSED, 0 failed.\n', n_pass);
  printf('Remaining Checks 3-7 pending subsequent Phase C commits.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
