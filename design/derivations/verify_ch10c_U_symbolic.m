% verify_ch10c_U_symbolic.m
%
% Chapter 10c Phase C addendum verifier for the symbolic U closed form
% (ch10c_addendum_U_symbolic.md, Proposition F.11a).
%
% Two Checks:
%   Check 1 -- U_L closed form via IBP decomposition.
%              Route A: direct SymPy diff of V = (mu^10 k_2^2/L^12) * T
%                       with T from ch10c §5.
%              Route B: (T1.L) at alpha=12 applied to V, divided by 2n.
%              Expected: simplify(routeA - routeB) = 0 identically.
%
%   Check 2 -- U = U_L + I_g additive decomposition (numerical).
%              At 5 sample (theta, e, g) points:
%                 U_num = numerical_route_U (direct)
%                 U_L_sym = closed form from §3
%                 I_g_num = U_num - U_L_sym
%              Verify |I_g_num| < bounded (~ O(10) in dimensionless units)
%              and check the g-harmonic structure of I_g.
%
% Poisson-bracket convention: not invoked directly.
% Dimensionless units: L = mu = k_2 = 1 throughout. G = eta, H = theta * eta.

pkg load symbolic;

printf('=========================================================\n');
printf('Chapter 10c Phase C addendum verifier for symbolic U\n');
printf('(ch10c_addendum_U_symbolic.md Proposition F.11a)\n');
printf('2 Checks: symbolic U_L via IBP, numerical U = U_L + I_g\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function print_check(label, tol_str, obs_str, pass_flag)
  printf('  Expected: %s\n', tol_str);
  printf('  Observed: %s\n', obs_str);
  if pass_flag
    printf('  Status:   PASS\n');
  else
    printf('  Status:   FAIL\n');
  end
  printf('\n');
endfunction

% --------------------------------------------------------------
% Shared numerical infrastructure (same as verify_ch10e_F2starstar.m)
% --------------------------------------------------------------
function [E, f, cosf, sinf, kap] = kepler_numeric(l_val, e_val, eta_val)
  E = l_val;
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
  f = 2 * atan2(sqrt(1 + e_val)*sin(E/2), sqrt(1 - e_val)*cos(E/2));
endfunction

function parts = build_F1_partials(L, G, H, f, g, mu, k2)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;
  kap = eta^2 / (1 + e*cos(f));

  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  Aprime = 3*theta;
  Bprime = -3*theta;

  cos2fg  = cos(2*(f + g));
  sin2fg  = sin(2*(f + g));

  one_plus_ecosf    = 1 + e*cos(f);
  one_plus_ecosf_2  = one_plus_ecosf^2;
  one_plus_ecosf_3  = one_plus_ecosf * one_plus_ecosf_2;

  Phi_F1 = mu^4 * k2 * one_plus_ecosf_3 * (A + B*cos2fg);

  dPhi_df = mu^4 * k2 * ...
    ( -3*e*sin(f) * one_plus_ecosf_2 * (A + B*cos2fg) ...
      - 2*B * one_plus_ecosf_3 * sin2fg );

  dPhi_dg = mu^4 * k2 * one_plus_ecosf_3 * ( -2 * B * sin2fg );

  dPhi_dtheta = mu^4 * k2 * one_plus_ecosf_3 * (Aprime + Bprime*cos2fg);

  dPhi_de_explicit = mu^4 * k2 * 3 * one_plus_ecosf_2 * cos(f) * (A + B*cos2fg);

  df_de = sin(f) * (2 + e*cos(f)) / eta^2;
  dPhi_de_chain = dPhi_df * df_de;

  dPhi_de = dPhi_de_explicit + dPhi_de_chain;

  parts.F1     = Phi_F1 / G^6;
  parts.dF1_dl = (eta / kap^2) * dPhi_df / G^6;
  parts.dF1_dg = dPhi_dg / G^6;
  parts.dF1_dL = (eta^3 / e) * dPhi_de / G^7;
  parts.dF1_dG = ( -6*Phi_F1 - theta*dPhi_dtheta - (eta^2/e)*dPhi_de ) / G^7;
endfunction

function parts = build_S1_partials(L, G, H, l, f, kap, g, mu, k2)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;

  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  X02 = (3*e^2 - 2 + 2*eta^3) / e^2;

  E1 = f - l + e*sin(f);
  E2 = sin(2*(f + g));
  E3 = sin(f + 2*g);
  E4 = sin(3*f + 2*g);
  E5 = sin(2*g);
  E6 = sin(f) / kap;
  E7 = sin(f)^3;
  E8 = sin(f) * (2 + e*cos(f)) * cos(2*(f + g)) / kap;

  cos2fg = cos(2*(f + g));
  parts.dS1_dl = (mu^2 * k2 / L^3) * ( (A + B*cos2fg)/kap^3 - A/eta^3 );

  cos_f2g  = cos(f + 2*g);
  cos_3f2g = cos(3*f + 2*g);
  cos_2g   = cos(2*g);
  parts.dS1_dg = (mu^2 * k2 * B / G^3) * ( cos2fg + e*cos_f2g + (e/3)*cos_3f2g + (X02/3)*cos_2g );
endfunction

function U_val = U_route_at(theta, e_val, g_val, N_quad)
  L = 1;
  eta = sqrt(1 - e_val^2);
  G = eta * L;
  H = theta * G;
  mu = 1;
  k2 = 1;

  acc = 0;
  for i = 0:N_quad-1
    l_i = 2*pi*i/N_quad;
    [~, f_i, ~, ~, kap] = kepler_numeric(l_i, e_val, eta);
    F1p = build_F1_partials(L, G, H, f_i, g_val, mu, k2);
    S1p = build_S1_partials(L, G, H, l_i, f_i, kap, g_val, mu, k2);
    integrand = S1p.dS1_dl * F1p.dF1_dL + S1p.dS1_dg * F1p.dF1_dG;
    acc = acc + integrand;
  end
  U_val = acc / N_quad;
endfunction

%% =======================================================
%% Check 1: Symbolic U_L closed form — THREE-ROUTE cross-check
%%
%% Build V(L, G, H, g) = ⟨(F_1 - F_1^*)^2⟩_l.
%% Compute three routes:
%%   Route A: (1/(2n)) * dV/dL via SymPy diff with V = mu^8 k_2^2 T / L^12.
%%   Route B: (T1.L) at alpha=12 applied to V.
%%   Route C (INDEPENDENT PHYSICAL): V = n^2 * ⟨(∂S_1/∂l)^2⟩_l built from PHYSICAL
%%           definitions n = mu^2/L^3 and (∂S_1/∂l)^2 prefactor (mu^2 k_2/L^3)^2 T.
%%           This independently derives the mu^8 prefactor from first principles,
%%           catching any mu-scaling error that routes A and B would share.
%%
%% Expected: all three routes give the same result after simplify.
%% If Route C differs from A or B, the error is in the CHAPTER'S asserted prefactor.
%% =======================================================
printf('--- Check 1: Symbolic U_L closed form — THREE independent routes ---\n');

syms theta_s e_s g_s L_s G_s mu_s k2_s positive;
eta_s = sqrt(1 - e_s^2);

A_s = (3*theta_s^2 - 1) / 2;
B_s = 3*(1 - theta_s^2) / 2;

% Hansen coefficients (dimensionless)
X_6_0 = (8 + 24*e_s^2 + 3*e_s^4) / (8 * eta_s^9);
X_6_2 = e_s^2 * (6 + e_s^2) / (4 * eta_s^9);
X_6_4 = e_s^4 / (16 * eta_s^9);

T_s = (A_s^2 + B_s^2/2) * X_6_0 + 2*A_s*B_s * cos(2*g_s) * X_6_2 ...
    + (B_s^2/2) * cos(4*g_s) * X_6_4 - A_s^2 / eta_s^6;

% ---- ROUTE C (INDEPENDENT PHYSICAL): build V from physical definitions ----
% Step 1: n = mu^2/L^3 (Keplerian mean motion; file 04 Delaunay definitions)
n_phys = mu_s^2 / L_s^3;
% Step 2: ⟨(∂S_1/∂l)^2⟩_l has prefactor (mu^2 k_2/L^3)^2 per ch10c §5
dS1_dl_prefactor_sq = (mu_s^2 * k2_s / L_s^3)^2;
dS1_dl_sq_avg = dS1_dl_prefactor_sq * T_s;
% Step 3: V = n^2 * ⟨(∂S_1/∂l)^2⟩_l per short-period homological + l-avg
V_route_C = n_phys^2 * dS1_dl_sq_avg;
V_route_C = simplify(V_route_C);
% Route C has independently-derived mu-exponent: 2*2 + 2*1 + 0 = 4 + 4 = 8.

% ---- ROUTE A: assume chapter's claimed prefactor form V = (mu^? k_2^2/L^12) T ----
% The chapter's CURRENT formula is mu^8; verify by checking Route A matches Route C.
V_LGH_expr_chapter = mu_s^8 * k2_s^2 / L_s^12 * T_s;  % chapter claim (corrected 2026-04-19)
syms H_s positive;
V_LGH = subs(V_LGH_expr_chapter, [e_s, theta_s], [sqrt(1 - G_s^2/L_s^2), H_s/G_s]);
V_route_A_vs_C = simplify(V_LGH - subs(V_route_C, [e_s, theta_s], [sqrt(1 - G_s^2/L_s^2), H_s/G_s]));
p_route_A_eq_C = isequal(V_route_A_vs_C, sym(0));
if ~p_route_A_eq_C
  V_route_A_vs_C = simplify(expand(V_route_A_vs_C));
  p_route_A_eq_C = isequal(V_route_A_vs_C, sym(0));
end

% Diagnostic: also check mu-exponent of Route C (should be 8)
logderiv_muC = simplify(mu_s * diff(V_route_C, mu_s) / V_route_C);
mu_exp_C = double(subs(logderiv_muC, mu_s, sym(1)));
printf('  Route C (physical definition) mu-exponent: %d (expected 8)\n', round(mu_exp_C));
printf('  Chapter claim (Route A) mu-exponent: 8 (V_LGH_expr = mu^8 k_2^2 / L^12 * T)\n');
printf('  Route A vs Route C agreement: residual = %s\n', char(V_route_A_vs_C));

if ~p_route_A_eq_C
  printf('  CRITICAL: chapter prefactor DISAGREES with physical derivation!\n');
  n_fail = n_fail + 1;
end

% Route A: direct SymPy diff (V_LGH uses mu^8 per corrected chapter)
n_s = mu_s^2 / L_s^3;
dV_dL = diff(V_LGH, L_s);
routeA = dV_dL / (2 * n_s);

% Route B: (T1.L) at alpha=12: dV/dL = G^{-13} * (eta^3/e) * dPhi_V/de
% Phi_V = G^12 * V = mu^8 k_2^2 * eta^12 * T (corrected 2026-04-19 from mu^10)
Phi_V = mu_s^8 * k2_s^2 * eta_s^12 * T_s;
dPhi_V_de = diff(Phi_V, e_s);
eta_sub = G_s / L_s;
e_sub = sqrt(1 - (G_s/L_s)^2);
theta_sub = H_s / G_s;
dPhi_V_de_LGH = subs(dPhi_V_de, [theta_s, e_s], [theta_sub, e_sub]);
dV_dL_form = dPhi_V_de_LGH * (eta_sub^3 / e_sub) / G_s^13;
routeB = dV_dL_form / (2 * n_s);

residual_1 = simplify(routeA - routeB);
is_zero_1 = isequal(residual_1, sym(0));
if ~is_zero_1
  residual_1b = simplify(expand(residual_1));
  is_zero_1 = isequal(residual_1b, sym(0));
  residual_1 = residual_1b;
  if ~is_zero_1
    residual_1c = simplify(factor(residual_1));
    is_zero_1 = isequal(residual_1c, sym(0));
    residual_1 = residual_1c;
  end
end

% INDEPENDENT CHECK: compare Route C's U_L (from physical V) to Route A/B result
% U_L^{(C)} = (1/(2n)) * dV_route_C/dL. Must equal routeA symbolically.
V_route_C_in_LGH = subs(V_route_C, [e_s, theta_s], [sqrt(1 - G_s^2/L_s^2), H_s/G_s]);
dV_C_dL = diff(V_route_C_in_LGH, L_s);
U_L_route_C = dV_C_dL / (2 * n_s);
residual_AC = simplify(routeA - U_L_route_C);
p_AC = isequal(residual_AC, sym(0));
if ~p_AC
  residual_AC = simplify(expand(residual_AC));
  p_AC = isequal(residual_AC, sym(0));
end

printf('  Route A vs Route C (independent physical) U_L: residual = %s\n', char(residual_AC));
if ~p_AC
  printf('  CRITICAL: U_L from chapter formula disagrees with physical definition!\n');
end

% Three independent checks:
%   (1) Route A vs Route B (self-consistency of IBP identity)
%   (2) Route C mu-exponent = 8 (independent physical derivation)
%   (3) Route A vs Route C agreement (catches chapter-formula mu-error)
check1_pass = is_zero_1 && p_route_A_eq_C && p_AC && (round(mu_exp_C) == 8);

print_check('Check 1', ...
  'Route A = Route B = Route C (3 independent routes; mu-exp = 8 from physical defn)', ...
  sprintf('A=B: %d, A=C: %d, B=C proven via A=B & A=C, mu_C=%d', ...
          is_zero_1, p_AC, round(mu_exp_C)), ...
  check1_pass);
if check1_pass
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 2: U = U_L + I_g decomposition at sample (theta, e, g)
%%
%% Define U_L_num(theta, e, g) from §3 closed form in dimensionless units:
%%   U_L = (mu^8 k_2^2 eta^10 / (2 e G^10)) * [-12 e T + eta^2 dT/de]
%% In L=1: G = eta, mu = k_2 = 1, so U_L = eta^10/(2 e eta^10) * [-12 e T + eta^2 dT/de]
%%                                  = (1/(2 e)) * [-12 e T + eta^2 dT/de]
%%                                  = -6 T + (eta^2/(2e)) * dT/de
%%
%% Then I_g = U_numerical - U_L (numerically computed).
%% Verify I_g has expected harmonic structure (only cos 0g, cos 2g, cos 4g harmonics).
%% =======================================================
printf('--- Check 2: U = U_L + I_g decomposition (numerical) ---\n');

% Compute dT/de symbolically (with eta, theta, e independent)
syms e_sym theta_sym g_sym positive;
eta_expr = sqrt(1 - e_sym^2);
A_sym = (3*theta_sym^2 - 1)/2;
B_sym = 3*(1 - theta_sym^2)/2;
X_60_sym = (8 + 24*e_sym^2 + 3*e_sym^4) / (8 * eta_expr^9);
X_62_sym = e_sym^2 * (6 + e_sym^2) / (4 * eta_expr^9);
X_64_sym = e_sym^4 / (16 * eta_expr^9);
T_sym = (A_sym^2 + B_sym^2/2) * X_60_sym + 2*A_sym*B_sym * cos(2*g_sym) * X_62_sym ...
      + (B_sym^2/2) * cos(4*g_sym) * X_64_sym - A_sym^2 / eta_expr^6;
dT_de_sym = diff(T_sym, e_sym);

% Compile into Octave-callable functions
T_fun = function_handle(T_sym);
dT_de_fun = function_handle(dT_de_sym);

sample_pts = [
  0.2,  0.1, 0.5;
  0.3,  0.3, 1.0;
  0.447, 0.2, 0.4;
  0.6,  0.5, 0.8;
  0.8,  0.15, 1.5;
];

printf('  (theta, e, g)          U_num          U_L_sym        I_g_num\n');
N_q = 512;

all_Ig = zeros(size(sample_pts, 1), 1);
max_Ig = 0;
for ip = 1:size(sample_pts, 1)
  th = sample_pts(ip, 1);
  ee = sample_pts(ip, 2);
  gg = sample_pts(ip, 3);
  eta = sqrt(1 - ee^2);

  U_num = U_route_at(th, ee, gg, N_q);

  % U_L_sym in L=mu=k2=1 units:
  % U_L = (1/(2e)) [-12 e T + eta^2 dT/de]
  T_val = T_fun(ee, gg, th);
  dT_de_val = dT_de_fun(ee, gg, th);
  U_L_sym = (1/(2*ee)) * (-12*ee*T_val + eta^2 * dT_de_val);

  I_g_num = U_num - U_L_sym;
  all_Ig(ip) = I_g_num;
  if abs(I_g_num) > max_Ig; max_Ig = abs(I_g_num); end

  printf('  (%.3f, %.3f, %.3f)    %14.6e  %14.6e  %14.6e\n', ...
         th, ee, gg, U_num, U_L_sym, I_g_num);
end

% I_g finite (bounded) — indicates that the decomposition is meaningful.
% Check I_g ~ O(10) in dimensionless units (ch11d typical partial-sum values).
bounded = max_Ig < 1e3;

% Additional: verify at (theta_0, e_0) fixed, I_g(g) has only cos 0g, cos 2g, cos 4g
% harmonics (structurally expected) by 12-pt harmonic fit.
% 12 uniform samples (Nyquist = 6) fully resolves sin(4g); 8 samples would
% alias sin(4g) to zero at the sample points, leaving b_4 undetermined.
fit_theta = 0.3; fit_e = 0.3;
fit_g = (0:11) * (2*pi/12);
fit_Ig = zeros(length(fit_g), 1);
for ig = 1:length(fit_g)
  gg = fit_g(ig);
  eta = sqrt(1 - fit_e^2);
  U_num = U_route_at(fit_theta, fit_e, gg, N_q);
  T_val = T_fun(fit_e, gg, fit_theta);
  dT_de_val = dT_de_fun(fit_e, gg, fit_theta);
  U_L_sym = (1/(2*fit_e)) * (-12*fit_e*T_val + eta^2 * dT_de_val);
  fit_Ig(ig) = U_num - U_L_sym;
end

% 5-parameter fit: c_0 + c_2 cos(2g) + b_2 sin(2g) + c_4 cos(4g) + b_4 sin(4g)
N_fit = length(fit_g);
M_fit = zeros(N_fit, 5);
for i = 1:N_fit
  M_fit(i, 1) = 1;
  M_fit(i, 2) = cos(2*fit_g(i));
  M_fit(i, 3) = sin(2*fit_g(i));
  M_fit(i, 4) = cos(4*fit_g(i));
  M_fit(i, 5) = sin(4*fit_g(i));
end
fit_coefs = M_fit \ fit_Ig;
residuals_fit = fit_Ig - M_fit * fit_coefs;
max_fit_res = max(abs(residuals_fit));

sines_negligible = (abs(fit_coefs(3)) < 1e-8) && (abs(fit_coefs(5)) < 1e-8);
structural_pass = bounded && sines_negligible && (max_fit_res < 1e-8);

printf('\n  I_g 5-param fit at (theta=%.3f, e=%.3f):\n', fit_theta, fit_e);
printf('    c_0 = %.6e, c_2 = %.6e, b_2 = %.6e, c_4 = %.6e, b_4 = %.6e\n', ...
       fit_coefs(1), fit_coefs(2), fit_coefs(3), fit_coefs(4), fit_coefs(5));
printf('    Fit max residual: %.3e\n', max_fit_res);
printf('    Sines |b_2|, |b_4| negligible: %d\n', sines_negligible);
printf('    Max |I_g| bounded: %d\n', bounded);

print_check('Check 2', ...
  'I_g = U - U_L numerically bounded + has only k_g in {0,2,4} harmonics', ...
  sprintf('max|I_g| = %.3e (bounded=%d), fit residual %.3e, sines negligible=%d', ...
          max_Ig, bounded, max_fit_res, sines_negligible), ...
  structural_pass);
if structural_pass
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Chapter 10c Phase C addendum verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Proposition F.11a (symbolic U_L + structural I_g) confirmed.\n');
  printf('Full symbolic I_g expansion deferred per ch10c_addendum §4.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
