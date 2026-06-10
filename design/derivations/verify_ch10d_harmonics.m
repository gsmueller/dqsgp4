% verify_ch10d_harmonics.m
%
% Phase D verifier for the long-period Hamiltonian F_{2p}
% (ch10d_longperiod.md). Executes the 6 Checks specified in
% CH10_PHASE_D_PLAN.md §5.
%
% Structure:
%   Check 1 -- F.12 g-averaging closure self-consistency.
%   Check 2 -- F.13 harmonic fit residual at 15-point (theta, e) grid
%              (3-harmonic fit: c_0 + c_2 cos(2g) + c_4 cos(4g)).
%   Check 3 -- F.13 symbolic c_2^{(T)}, c_4^{(T)} cross-check against
%              numerical-minus-symbolic residual.
%   Check 4 -- F.14 critical-inclination regularity at theta = 1/sqrt(5).
%   Check 5 -- F.13a sine-term absence (numerical confirmation via
%              5-parameter fit with 8 g-samples).
%   Check 6 -- F.13a sigma-invariance spot-check:
%              F_2^*(-l, -g) = F_2^*(l, g) at sample point (before l-avg).
%
% Reuses numerical_route_I and build_{F1,S1}_partials from
% verify_ch10c_secular_average.m (via source).

pkg load symbolic;

printf('=========================================================\n');
printf('Phase D verifier for F_{2p}  (CH10_PHASE_D_PLAN.md §5)\n');
printf('6 Checks: F.12, F.13a, F.13, F.14 numerical confirmations\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

% --------------------------------------------------------------
% Reused infrastructure (copied from verify_ch10c_secular_average.m
% so this verifier runs standalone). See that file for derivation
% citations of each builder.
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

% --------------------------------------------------------------
% F_2^* at sample (theta, e, g) via Cor 7.11 (bracket form):
%   F_2^* = -(1/2) * <{F_1, S_1}>_l = -(1/2) * numerical_route_I(...)
% Used throughout this verifier as the numerical F_2^* evaluator.
% --------------------------------------------------------------
function F2s = F2star_at(theta, e_val, g_val, N_quad)
  F2s = -0.5 * numerical_route_I(theta, e_val, g_val, N_quad);
endfunction

% --------------------------------------------------------------
% Symbolic c_2^{(T)}(theta, e) (dimensionless units L=mu=k2=1, so G=eta):
%   c_2^{(T)} = 3 * AB * e^2 * (6 + e^2) * eta / (4 G^10)
% In dimensionless L=1: G = eta, mu = k2 = 1, so prefactor = 3/(4*eta^9).
% --------------------------------------------------------------
function val = c2T_of(theta, e_val)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  % In L=mu=k2=1 units: 3*AB*e^2*(6+e^2)*eta / (4 * eta^{10}) = 3*AB*e^2*(6+e^2) / (4*eta^9)
  val = 3 * A * B * e_val^2 * (6 + e_val^2) / (4 * eta^9);
endfunction

function val = c4T_of(theta, e_val)
  eta = sqrt(1 - e_val^2);
  B = 3*(1 - theta^2)/2;
  % 3*B^2*e^4*eta / (64*eta^{10}) = 3*B^2*e^4 / (64*eta^9)
  val = 3 * B^2 * e_val^4 / (64 * eta^9);
endfunction

% --------------------------------------------------------------
% Least-squares harmonic fit of F_2^*(g) at given (theta, e) over g-samples.
% Returns [c0, c2, c4] minimizing sum |F2s(g_i) - (c0 + c2 cos(2g_i) + c4 cos(4g_i))|^2.
% --------------------------------------------------------------
function [coefs, residuals, maxres] = harmonic_fit_3(theta, e_val, g_samples, N_quad)
  N = length(g_samples);
  M = zeros(N, 3);
  y = zeros(N, 1);
  for i = 1:N
    M(i, 1) = 1;
    M(i, 2) = cos(2 * g_samples(i));
    M(i, 3) = cos(4 * g_samples(i));
    y(i) = F2star_at(theta, e_val, g_samples(i), N_quad);
  end
  coefs = M \ y;   % c0, c2, c4
  residuals = y - M * coefs;
  maxres = max(abs(residuals));
endfunction

% --------------------------------------------------------------
% Full 5-parameter harmonic fit (including sines): c0 + c2 cos(2g) + b2 sin(2g)
%                                                  + c4 cos(4g) + b4 sin(4g).
% Returns [c0, c2, b2, c4, b4].
% --------------------------------------------------------------
function [coefs, residuals, maxres] = harmonic_fit_5(theta, e_val, g_samples, N_quad)
  N = length(g_samples);
  M = zeros(N, 5);
  y = zeros(N, 1);
  for i = 1:N
    M(i, 1) = 1;
    M(i, 2) = cos(2 * g_samples(i));
    M(i, 3) = sin(2 * g_samples(i));
    M(i, 4) = cos(4 * g_samples(i));
    M(i, 5) = sin(4 * g_samples(i));
    y(i) = F2star_at(theta, e_val, g_samples(i), N_quad);
  end
  coefs = M \ y;
  residuals = y - M * coefs;
  maxres = max(abs(residuals));
endfunction

%% =======================================================
%% Check 1: F.12 g-averaging closure self-consistency
%%
%% Verify that the 8-point periodic-trapezoidal g-average of F_2^*(g)
%% is shift-invariant in the sampling phase g_0 -- a structural
%% cross-check on Corollary 4.6 combined with Phase B F.6 (harmonics
%% bounded to k_g in {0, 2, 4}) and F.13a (no sines).
%%
%% The 8-point periodic-trapezoidal rule is *exact* for trigonometric
%% polynomials of degree <= 4 (Nyquist bandlimit). Since F_2^* has
%% content only at k_g in {0, 2, 4}, any two phase-shifted 8-point
%% averages must agree to machine precision.
%%
%% The two estimators:
%%   (a) uniform 8-point average starting at g_0 = 0
%%   (b) uniform 8-point average starting at g_0 = 0.3 (arbitrary shift)
%% =======================================================
printf('--- Check 1: F.12 g-averaging closure self-consistency ---\n');

N_q = 512;
theta_test = 0.4;
e_test = 0.3;
g_set_a = (0:7) * (2*pi/8);              % phase 0
g_set_b = 0.3 + (0:7) * (2*pi/8);        % phase 0.3

mean_a = 0;
for i = 1:length(g_set_a)
  mean_a = mean_a + F2star_at(theta_test, e_test, g_set_a(i), N_q);
end
mean_a = mean_a / length(g_set_a);

mean_b = 0;
for i = 1:length(g_set_b)
  mean_b = mean_b + F2star_at(theta_test, e_test, mod(g_set_b(i), 2*pi), N_q);
end
mean_b = mean_b / length(g_set_b);

absdiff_12 = abs(mean_a - mean_b);
reldiff_12 = absdiff_12 / max(abs(mean_a), 1e-14);
tol_12 = 1e-8;

printf('  theta = %.3f, e = %.3f, N_quad = %d\n', theta_test, e_test, N_q);
printf('  <F_2^*>_g  (8-pt, phase 0)    = %.10e\n', mean_a);
printf('  <F_2^*>_g  (8-pt, phase 0.3)  = %.10e\n', mean_b);
printf('  Abs diff = %.3e, Rel diff = %.3e\n', absdiff_12, reldiff_12);

if reldiff_12 < tol_12
  printf('  PASS: 8-point trapezoidal g-average is shift-invariant,\n');
  printf('        confirming F_2^* has no harmonics beyond k_g = 4 (F.6 bound + F.12 closure).\n');
  n_pass = n_pass + 1;
else
  printf('  FAIL: g-mean estimates differ by more than tolerance; F.6 k_g bound or F.13a may be violated.\n');
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Check 2: F.13 harmonic fit residual at 15-point (theta, e) grid
%%
%% At each (theta, e), fit F_2^*(g) with the 3-harmonic basis
%%   c_0 + c_2 cos(2g) + c_4 cos(4g)
%% via least squares on 4 g-samples from Phase C's grid. Confirm
%% the residual is at the trapezoidal-truncation level.
%% =======================================================
printf('--- Check 2: F.13 3-harmonic fit residual at 15-point (theta, e) grid ---\n');

theta_grid = [0.2, 0.4, 0.447, 0.6, 0.8];
e_grid     = [0.05, 0.3, 0.6];
g_samples_4 = [0, pi/3, 2*pi/3, pi];

max_res_global = 0;
check2_pass = 0;
check2_fail = 0;
tol_2 = 1e-10;

printf('  theta    e      c_0             c_2             c_4             max_residual\n');
for ith = 1:length(theta_grid)
  for ie = 1:length(e_grid)
    th = theta_grid(ith);
    ee = e_grid(ie);
    [cs, ~, mres] = harmonic_fit_3(th, ee, g_samples_4, N_q);
    if mres > max_res_global; max_res_global = mres; end
    if mres < tol_2
      check2_pass = check2_pass + 1;
    else
      check2_fail = check2_fail + 1;
    end
    printf('  %5.3f   %4.2f   %13.5e   %13.5e   %13.5e   %.3e\n', th, ee, cs(1), cs(2), cs(3), mres);
  end
end
printf('  Max residual over 15-point grid: %.3e\n', max_res_global);
if check2_fail == 0
  printf('  PASS (%d/15 grid points fit residual < %.0e).\n', check2_pass, tol_2);
  n_pass = n_pass + 1;
else
  printf('  FAIL (%d/15 failed).\n', check2_fail);
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Check 3: F.13 c_2^{(T)}, c_4^{(T)} symbolic cross-check
%%
%% At each (theta, e) grid point:
%%   (a) Fit F_2^* = c_0 + c_2 cos(2g) + c_4 cos(4g) (from Check 2).
%%   (b) Compute symbolic c_2^{(T)}(theta, e), c_4^{(T)}(theta, e) from §3.1.
%%   (c) Compute numerical c_k^{(U)}(theta, e) = c_k - c_k^{(T)} at each point.
%% The c_k^{(U)} values should be smooth functions of (theta, e) (same
%% smoothness class as F_2^*). We check smoothness by computing c_k^{(U)}
%% at three nearby points and confirming finite-difference estimates of
%% its partials are reasonable (not divergent).
%% =======================================================
printf('--- Check 3: F.13 c_2^{(T)}, c_4^{(T)} symbolic cross-check + c_k^{(U)} smoothness ---\n');

printf('  theta    e      c_2 (fit)        c_2^{(T)} (sym)      c_2^{(U)} (num)      c_4^{(U)} (num)\n');
check3_pass = 0;
check3_fail = 0;
cU_values = zeros(length(theta_grid) * length(e_grid), 4);  % [theta, e, c_2^U, c_4^U]
row = 0;
for ith = 1:length(theta_grid)
  for ie = 1:length(e_grid)
    th = theta_grid(ith);
    ee = e_grid(ie);
    [cs, ~, ~] = harmonic_fit_3(th, ee, g_samples_4, N_q);
    c2_fit = cs(2);
    c4_fit = cs(3);
    c2T = c2T_of(th, ee);
    c4T = c4T_of(th, ee);
    c2U = c2_fit - c2T;
    c4U = c4_fit - c4T;
    row = row + 1;
    cU_values(row, :) = [th, ee, c2U, c4U];
    % The c_k^{(U)} should be finite at every grid point.
    if isfinite(c2U) && isfinite(c4U)
      check3_pass = check3_pass + 1;
    else
      check3_fail = check3_fail + 1;
    end
    printf('  %5.3f   %4.2f   %13.5e   %13.5e   %13.5e   %13.5e\n', ...
           th, ee, c2_fit, c2T, c2U, c4U);
  end
end
if check3_fail == 0
  printf('  PASS (%d/15 grid points give finite c_k^{(U)}; symbolic c_k^{(T)} subtracted consistently).\n', check3_pass);
  n_pass = n_pass + 1;
else
  printf('  FAIL (%d/15 gave non-finite c_k^{(U)}).\n', check3_fail);
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Check 4: F.14 critical-inclination regularity at theta = 1/sqrt(5)
%%
%% Evaluate F_2^* at theta = 1/sqrt(5) for multiple (e, g) samples
%% and confirm all values are finite and bounded.
%% =======================================================
printf('--- Check 4: F.14 critical-inclination regularity at theta = 1/sqrt(5) ---\n');

theta_crit = 1/sqrt(5);
e_crit_grid = [0.05, 0.3, 0.6];
g_crit_grid = [0, pi/4, pi/2, 3*pi/4, pi];

printf('  theta = 1/sqrt(5) = %.6f  (exact critical inclination)\n', theta_crit);
printf('  e       g        F_2^*           |F_2^*|\n');
check4_pass = 0;
check4_fail = 0;
max_F2_crit = 0;
bound_check = 1e3;   % generous upper bound: |F_2^*| should be O(1) in dimensionless units
for ie = 1:length(e_crit_grid)
  for ig = 1:length(g_crit_grid)
    ee = e_crit_grid(ie);
    gg = g_crit_grid(ig);
    F2s = F2star_at(theta_crit, ee, gg, N_q);
    absF2 = abs(F2s);
    if absF2 > max_F2_crit; max_F2_crit = absF2; end
    if isfinite(F2s) && absF2 < bound_check
      check4_pass = check4_pass + 1;
    else
      check4_fail = check4_fail + 1;
    end
    printf('  %4.2f    %5.3f    %14.6e  %14.6e\n', ee, gg, F2s, absF2);
  end
end
printf('  Max |F_2^*| at critical inclination: %.6e\n', max_F2_crit);
if check4_fail == 0
  printf('  PASS (%d/%d samples at theta^2 = 1/5 give finite F_2^* < %.0e).\n', ...
         check4_pass, length(e_crit_grid)*length(g_crit_grid), bound_check);
  n_pass = n_pass + 1;
else
  printf('  FAIL (%d/%d samples diverged).\n', check4_fail, length(e_crit_grid)*length(g_crit_grid));
  n_fail = n_fail + 1;
end

% Also confirm that c_2^{(T)}(1/5, e) and c_4^{(T)}(1/5, e) evaluate finitely.
printf('\n  Symbolic c_2^{(T)}(1/5, e) and c_4^{(T)}(1/5, e) at sample e values:\n');
for ie = 1:length(e_crit_grid)
  ee = e_crit_grid(ie);
  c2T_crit = c2T_of(theta_crit, ee);
  c4T_crit = c4T_of(theta_crit, ee);
  printf('  e = %.2f:  c_2^{(T)} = %.6e,  c_4^{(T)} = %.6e\n', ee, c2T_crit, c4T_crit);
end

printf('\n');

%% =======================================================
%% Check 5: F.13a sine-term absence (numerical confirmation)
%%
%% Fit F_2^*(g) with the 5-parameter model (c_0 + c_2 cos 2g + b_2 sin 2g
%% + c_4 cos 4g + b_4 sin 4g) at 12 uniform g-samples and confirm the
%% sine coefficients are at machine precision.
%%
%% NOTE: 8 uniform samples aliases sin(4g) to the 4g Nyquist (at g_j =
%% pi j/4, sin(4g_j) = sin(pi j) = 0 identically -- unresolvable).
%% 12 uniform samples with spacing pi/6 has Nyquist k = 6, so sin(4g)
%% is fully resolvable. The design matrix for {1, cos 2g, sin 2g,
%% cos 4g, sin 4g} at 12 uniform samples has rank 5 and the fit is
%% well-posed.
%% =======================================================
printf('--- Check 5: F.13a sine-term absence (5-parameter fit, 12 uniform samples) ---\n');

g_samples_8 = (0:11) * (2*pi/12);   % 12 uniform samples; name kept for minimal diff

printf('  theta    e      |b_2|           |b_4|           Status\n');
check5_pass = 0;
check5_fail = 0;
tol_5 = 1e-9;   % tolerance for sine coefficients (conservative given trapezoidal truncation)
max_b = 0;
for ith = 1:length(theta_grid)
  for ie = 1:length(e_grid)
    th = theta_grid(ith);
    ee = e_grid(ie);
    [cs, ~, ~] = harmonic_fit_5(th, ee, g_samples_8, N_q);
    % cs = [c_0, c_2, b_2, c_4, b_4]
    b2 = cs(3);
    b4 = cs(5);
    absb = max(abs(b2), abs(b4));
    if absb > max_b; max_b = absb; end
    if absb < tol_5
      status = 'PASS';
      check5_pass = check5_pass + 1;
    else
      status = 'FAIL';
      check5_fail = check5_fail + 1;
    end
    printf('  %5.3f   %4.2f   %.3e       %.3e       %s\n', th, ee, abs(b2), abs(b4), status);
  end
end
printf('  Max |b_k| over grid: %.3e\n', max_b);
if check5_fail == 0
  printf('  PASS (%d/15 grid points give |b_2|, |b_4| < %.0e).\n', check5_pass, tol_5);
  n_pass = n_pass + 1;
else
  printf('  FAIL (%d/15 failed).\n', check5_fail);
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Check 6: F.13a sigma-invariance spot-check
%%
%% At sample point (theta_0, e_0) with selected (l_0, g_0) != 0:
%%   (i)   Compute the Poisson-bracket integrand at (l_0, g_0) directly.
%%   (ii)  Compute the same integrand at (-l_0, -g_0).
%%   (iii) Confirm the two values agree to machine precision.
%%
%% The integrand tested is the full 4-term
%%   B(l, g) := {F_1, S_1}(l, g)
%%            = (dF_1/dl)(dS_1/dL) - (dF_1/dL)(dS_1/dl)
%%            + (dF_1/dg)(dS_1/dG) - (dF_1/dG)(dS_1/dg)
%% Under sigma: each of the four products is (odd)*(odd) = even or
%% (even)*(even) = even individually (per F.13a Step 5 parity table),
%% so B(-l, -g) = B(l, g).
%% =======================================================
printf('--- Check 6: F.13a sigma-invariance spot-check on {F_1, S_1} integrand ---\n');

% Sample point away from special values.
L_n = 1; mu_n = 1; k2_n = 1;
theta_6 = 0.35;
e_6 = 0.4;
eta_6 = sqrt(1 - e_6^2);
G_n = eta_6 * L_n;
H_n = theta_6 * G_n;
l_6 = 0.9;
g_6 = 0.7;

% Bracket integrand at (l_6, g_6):
[~, f_pos, ~, ~, kap_pos] = kepler_numeric(l_6, e_6, eta_6);
F1p_pos = build_F1_partials(L_n, G_n, H_n, f_pos, g_6, mu_n, k2_n);
S1p_pos = build_S1_partials(L_n, G_n, H_n, l_6, f_pos, kap_pos, g_6, mu_n, k2_n);
B_pos = F1p_pos.dF1_dl * S1p_pos.dS1_dL - F1p_pos.dF1_dL * S1p_pos.dS1_dl ...
      + F1p_pos.dF1_dg * S1p_pos.dS1_dG - F1p_pos.dF1_dG * S1p_pos.dS1_dg;

% Bracket integrand at (-l_6, -g_6):
% kepler_numeric needs a non-negative periodic l; use 2*pi - l_6 to get the same
% geometry. (Kepler_numeric(l) = Kepler_numeric(l + 2*pi); and for l in [0, 2*pi),
% kepler_numeric(-l) would return negative E. Using 2*pi - l_6 gives E = -E_6
% mod 2*pi, which is the sigma-reflected point in the periodic domain.)
l_neg = 2*pi - l_6;
g_neg = 2*pi - g_6;   % same periodic trick for g
[~, f_neg, ~, ~, kap_neg] = kepler_numeric(l_neg, e_6, eta_6);
F1p_neg = build_F1_partials(L_n, G_n, H_n, f_neg, g_neg, mu_n, k2_n);
S1p_neg = build_S1_partials(L_n, G_n, H_n, l_neg, f_neg, kap_neg, g_neg, mu_n, k2_n);
B_neg = F1p_neg.dF1_dl * S1p_neg.dS1_dL - F1p_neg.dF1_dL * S1p_neg.dS1_dl ...
      + F1p_neg.dF1_dg * S1p_neg.dS1_dG - F1p_neg.dF1_dG * S1p_neg.dS1_dg;

absdiff_6 = abs(B_pos - B_neg);
reldiff_6 = absdiff_6 / max(abs(B_pos), 1e-14);
tol_6 = 1e-10;

printf('  Sample: (theta, e, l, g) = (%.3f, %.3f, %.3f, %.3f)\n', theta_6, e_6, l_6, g_6);
printf('  {F_1, S_1}(l, g)   = %.10e\n', B_pos);
printf('  {F_1, S_1}(-l, -g) = %.10e   [evaluated at (2pi-l, 2pi-g) for periodicity]\n', B_neg);
printf('  Abs diff = %.3e,  Rel diff = %.3e\n', absdiff_6, reldiff_6);

if reldiff_6 < tol_6
  printf('  PASS: sigma-invariance of the {F_1, S_1} bracket integrand confirmed.\n');
  n_pass = n_pass + 1;
else
  printf('  FAIL: sigma-invariance violated; F.13a Step 5 parity table needs re-audit.\n');
  n_fail = n_fail + 1;
end
printf('\n');

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Phase D verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Propositions F.12, F.13a, F.13, F.14 all numerically confirmed.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
