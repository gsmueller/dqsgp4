% verify_ch10e_F2starstar.m
%
% Chapter 10 Phase E verifier for the F_2^{**} closed form
% (ch10e_F2starstar_closed_form.md, Proposition F.15).
%
% Two Checks:
%   Check 1 -- Symbolic c_0^{(T)} closed form vs direct g-averaging of
%              -(n'/2) * (mu^2 k_2/L^3)^2 * T(theta, e, eta, g).
%              Expected: simplify(routeA - routeB) = 0 identically.
%   Check 2 -- Numerical <F_2^*>_g (via F2star_at bracket route Cor 7.11)
%              vs c_0^{(T)} (closed form) + c_0^{(U)} (numerical U-route
%              average via Cor 7.12 U-component direct evaluator).
%              Expected: rel err < 1e-10 at 5 sample (theta, e) pts.
%
% Poisson-bracket convention: F_2^* was built via {·,·}_{mod} (ch10c §2
% F.10, Cor 7.11 form). This chapter inherits that convention; no
% brackets are invoked here directly.
%
% Dimensionless units: L = mu = k_2 = 1 throughout. Then G = eta,
% H = theta*eta, L^{-10} = eta^{-10}/(G^{-10}) = 1 (since G = eta and L=1).

pkg load symbolic;

printf('=========================================================\n');
printf('Chapter 10 Phase E verifier for F_2^{**}  (ch10e §1 F.15)\n');
printf('2 Checks: symbolic c_0^{(T)} and numerical c_0 = c_0^{(T)} + c_0^{(U)}\n');
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
% Reused infrastructure (same as verify_ch10c_secular_average.m and
% verify_ch10d_harmonics.m). Kept inline so this verifier runs
% standalone.
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

  parts.S1 = (mu^2 * k2 / G^3) * ( A*E1 + (B/2)*E2 + (e*B/2)*E3 + (e*B/6)*E4 + (B*X02/6)*E5 );

  cos2fg = cos(2*(f + g));
  parts.dS1_dl = (mu^2 * k2 / L^3) * ( (A + B*cos2fg)/kap^3 - A/eta^3 );

  cos_f2g  = cos(f + 2*g);
  cos_3f2g = cos(3*f + 2*g);
  cos_2g   = cos(2*g);
  parts.dS1_dg = (mu^2 * k2 * B / G^3) * ( cos2fg + e*cos_f2g + (e/3)*cos_3f2g + (X02/3)*cos_2g );

  c3_dL = B*eta^2 / (2*L*e);
  c4_dL = B*eta^2 / (6*L*e);
  c5_dL = B*eta^2 * (eta + 2) / (3*L*(1 + eta)^2);
  c6_dL = 3*A*eta^2 / (L*e);
  c7_dL = -A*e / L;
  c8_dL = B*eta^2 / (L*e);
  parts.dS1_dL = (mu^2 * k2 / G^3) * ( c3_dL*E3 + c4_dL*E4 + c5_dL*E5 + c6_dL*E6 + c7_dL*E7 + c8_dL*E8 );

  c1_dH = 3*theta;
  c2_dH = -3*theta/2;
  c3_dH = -3*e*theta/2;
  c4_dH = -e*theta/2;
  c5_dH = -theta*X02/2;
  parts.dS1_dH = (mu^2 * k2 / G^4) * ( c1_dH*E1 + c2_dH*E2 + c3_dH*E3 + c4_dH*E4 + c5_dH*E5 );

  parts.dS1_dG = -(3/G)*parts.S1 - theta*parts.dS1_dH - (1/eta)*parts.dS1_dL;
endfunction

% Numerical route I: <{F_1, S_1}>_l  (Cor 7.11 bracket form)
function I_val = numerical_route_I(theta, e_val, g_val, N_quad)
  L = 1;
  eta = sqrt(1 - e_val^2);
  G = eta * L;
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
  I_val = bracket_sum / N_quad;
endfunction

% F_2^* via Cor 7.11:  F_2^* = -(1/2)*<{F_1, S_1}>_l
function F2s = F2star_at(theta, e_val, g_val, N_quad)
  F2s = -0.5 * numerical_route_I(theta, e_val, g_val, N_quad);
endfunction

% U-component evaluator (Cor 7.12 U-route):
%   U(theta, e, g) = <(dS_1/dl)*(dF_1/dL) + (dS_1/dg)*(dF_1/dG)>_l
% Uniform trapezoidal in l with N_quad samples.
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

% Closed-form c_0^{(T)} (dimensionless L=mu=k2=1, so G=eta, L^{10}=1 but
% the closed form uses G^{-10}*eta = eta^{-10}*eta = eta^{-9}).
%   c_0^{(T)} = (3 mu^6 k2^2 eta)/(2 G^{10}) *
%               [ (A^2 + B^2/2)*(8 + 24 e^2 + 3 e^4)/8  -  A^2 eta^3 ]
function val = c0T_of(theta, e_val)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  % In L=mu=k2=1 units: G = eta, so G^{-10} = eta^{-10}. Prefactor = (3*eta)/(2*eta^{10}) = 3/(2*eta^9).
  prefactor = 3 / (2 * eta^9);
  bracket  = (A^2 + B^2/2)*(8 + 24*e_val^2 + 3*e_val^4)/8 - A^2 * eta^3;
  val = prefactor * bracket;
endfunction

%% =======================================================
%% Check 1: Symbolic c_0^{(T)} closed form vs direct g-averaging of
%%          -(n'/2) * (mu^2 k_2 / L^3)^2 * T(theta, e, eta, g)
%%
%% Route A (direct): Build T(theta, e, g) symbolically using A, B, and
%%   Hansen coefficients X_0^{-6,0}, X_0^{-6,2}, X_0^{-6,4}; g-average
%%   by SymPy integration; multiply by (3 mu^2)/(2 L^4) * (mu^2 k_2/L^3)^2.
%%
%% Route B (closed form): c_0^{(T)} from §1 F.15 closed form.
%%
%% Expected: simplify(routeA - routeB) = 0 identically.
%% =======================================================
printf('--- Check 1: Symbolic c_0^{(T)} closed form cross-check ---\n');

syms theta_s e_s g_s L_s G_s mu_s k2_s positive;
eta_s = sqrt(1 - e_s^2);

A_s = (3*theta_s^2 - 1) / 2;
B_s = 3*(1 - theta_s^2) / 2;

% Hansen coefficients (ch06 Corollaries B.0.7-9/-10/-11)
X_6_0 = (8 + 24*e_s^2 + 3*e_s^4) / (8 * eta_s^9);   % X_0^{-6,0}
X_6_2 = e_s^2 * (6 + e_s^2)     / (4 * eta_s^9);   % X_0^{-6,2}
X_6_4 = e_s^4                   / (16 * eta_s^9);  % X_0^{-6,4}

% T(theta, e, g) from ch10c §5 (with -A^2/eta^6 subtraction)
T_s = (A_s^2 + B_s^2/2) * X_6_0 ...
    + 2*A_s*B_s * cos(2*g_s) * X_6_2 ...
    + (B_s^2/2) * cos(4*g_s) * X_6_4 ...
    - A_s^2 / eta_s^6;

% g-average via SymPy integration
T_avg = int(T_s, g_s, 0, 2*sym(pi)) / (2*sym(pi));
T_avg = simplify(T_avg);

% Route A: -(n'/2) * (mu^2 k_2 / L^3)^2 * <T>_g with n' = -3*mu^2/L^4
% so -(n'/2) = 3*mu^2/(2*L^4)
prefactor_symbolic = (3*mu_s^2 / (2*L_s^4)) * (mu_s^2 * k2_s / L_s^3)^2;
route_A = prefactor_symbolic * T_avg;
% Substitute L = G/eta (Delaunay identity) so G^{-10} prefactor emerges.
route_A_in_G = subs(route_A, L_s, G_s / eta_s);
route_A_in_G = simplify(route_A_in_G);

% Route B: Closed-form c_0^{(T)} from §1 F.15
route_B = (3 * mu_s^6 * k2_s^2 * eta_s / (2 * G_s^10)) * ...
          ((A_s^2 + B_s^2/2) * (8 + 24*e_s^2 + 3*e_s^4) / 8 - A_s^2 * eta_s^3);
route_B = simplify(route_B);

residual_1 = simplify(route_A_in_G - route_B);
is_zero_1 = isequal(residual_1, sym(0));
if ~is_zero_1
  residual_1b = simplify(expand(residual_1));
  is_zero_1 = isequal(residual_1b, sym(0));
  if ~is_zero_1
    residual_1c = simplify(factor(residual_1b));
    is_zero_1 = isequal(residual_1c, sym(0));
    residual_1 = residual_1c;
  else
    residual_1 = residual_1b;
  end
end

print_check('Check 1', ...
  'simplify(routeA - routeB) = 0 identically', ...
  sprintf('residual = %s', char(residual_1)), ...
  is_zero_1);
if is_zero_1
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 2: Numerical <F_2^*>_g vs c_0^{(T)} + c_0^{(U)}
%%
%% At 5 sample (theta, e) points, compute:
%%   LHS = <F_2^*>_g via 8-point trapezoidal in g with F2star_at (Cor 7.11)
%%   c_T = c_0^{(T)} closed form
%%   c_U = <U>_g via 8-point trapezoidal in g with U_route_at (Cor 7.12)
%%   RHS = c_T + c_U
%%
%% Expected: max |LHS - RHS| / |LHS| < 1e-10 over 5 samples.
%% (trapezoidal truncation level; F_2^* has only k_g in {0, 2, 4} so
%%  8-point is exact in g; the numerical_route_I l-quadrature with
%%  N_quad = 512 gives rel err ~ 1e-12, well under tolerance.)
%% =======================================================
printf('--- Check 2: Numerical <F_2^*>_g = c_0^{(T)} + c_0^{(U)} ---\n');

N_q_l = 512;        % l-quadrature sample count
N_q_g = 8;          % g-quadrature: 8-point exact for k_g <= 4
g_samples = (0:N_q_g-1) * (2*pi/N_q_g);

% Sample (theta, e) grid
sample_pts = [ ...
  0.2,  0.1;
  0.3,  0.3;
  0.447, 0.2;   % near critical inclination for robustness check
  0.6,  0.5;
  0.8,  0.15;
];

printf('  theta    e      <F_2^*>_g       c_0^{(T)}       c_0^{(U)}       RHS             |LHS-RHS|/|LHS|\n');

check2_pass = 0;
check2_fail = 0;
max_rel_err_2 = 0;
tol_2 = 1e-10;

for ip = 1:size(sample_pts, 1)
  th = sample_pts(ip, 1);
  ee = sample_pts(ip, 2);

  % LHS: numerical <F_2^*>_g
  LHS_sum = 0;
  for ig = 1:N_q_g
    LHS_sum = LHS_sum + F2star_at(th, ee, g_samples(ig), N_q_l);
  end
  LHS = LHS_sum / N_q_g;

  % c_T: closed form
  c_T = c0T_of(th, ee);

  % c_U: numerical <U>_g
  cU_sum = 0;
  for ig = 1:N_q_g
    cU_sum = cU_sum + U_route_at(th, ee, g_samples(ig), N_q_l);
  end
  c_U = cU_sum / N_q_g;

  RHS = c_T + c_U;

  rel_err = abs(LHS - RHS) / max(abs(LHS), 1e-14);
  if rel_err > max_rel_err_2
    max_rel_err_2 = rel_err;
  end
  if rel_err < tol_2
    check2_pass = check2_pass + 1;
  else
    check2_fail = check2_fail + 1;
  end

  printf('  %5.3f   %4.2f   %13.5e   %13.5e   %13.5e   %13.5e   %.3e\n', ...
         th, ee, LHS, c_T, c_U, RHS, rel_err);
end

print_check('Check 2', ...
  sprintf('max |LHS - RHS|/|LHS| < %.0e over 5 samples', tol_2), ...
  sprintf('max rel err = %.3e (%d PASS / %d FAIL)', max_rel_err_2, check2_pass, check2_fail), ...
  check2_fail == 0);
if check2_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Chapter 10 Phase E verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Proposition F.15 (F_2^{**} closed form) confirmed.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
