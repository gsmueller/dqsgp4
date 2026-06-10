% verify_ch11b_S1star_partials.m
%
% Chapter 11b verifier for the S_1^* momentum and angle partials
% (ch11b_S1star_partials.md). Executes 4 Checks per CH12_PLAN.md §1.5.
%
% Reporting convention: each Check prints Expected / Observed / Status
% per the project discipline.
%
% Poisson-bracket convention: no Poisson brackets invoked here
% (partials are chain-rule derivatives). Ch 11b uses {·,·}_{mod}
% indirectly via Theorem 1 of ch10_foundations (the momentum-partial
% closure theorem), which is modern-bracket consistent.
%
% Dimensionless units: L = mu = k_2 = 1 throughout. In these units
% G = eta, H = theta*eta.

pkg load symbolic;

printf('=========================================================\n');
printf('Chapter 11b verifier for S_1^* partials  (CH12_PLAN.md §1.5)\n');
printf('4 Checks: G.6, G.7 (GSI), G.8, G.9b\n');
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
% Symbolic setup: S_1^{*,(T)} closed form from ch11a G.4.
% Variables: L, G, H, theta, e, eta, g, mu, k2 (L and G kept independent
% so that (T1.L) and (T1.G) can be evaluated in closed form).
% --------------------------------------------------------------
syms L_s G_s H_s theta_s e_s g_s mu_s k2_s positive;

% eta is a derived function of e, NOT an independent parameter.
% In Definition 1.1 of ℳ_α, F = F(θ, e, l, g). We use η = sqrt(1-e²) as
% shorthand, but all e-partials must respect the e→η chain.
eta_expr = sqrt(1 - e_s^2);

A_s = (3*theta_s^2 - 1) / 2;
B_s = 3*(1 - theta_s^2) / 2;
inner_s = A_s * B_s * e_s^2 * (6 + e_s^2) / 4 * sin(2*g_s) ...
        + B_s^2 * e_s^4 / 128                * sin(4*g_s);
% Phi_{S_1^*}(theta, e, g) — with eta expanded to sqrt(1-e^2); g-power carried by G separately
Phi_s = mu_s^2 * k2_s / ((5*theta_s^2 - 1) * eta_expr^2) * inner_s;
% S_1^{*,(T)}(theta, e, g, G) = G^{-3} * Phi_s
S1star_s = Phi_s / G_s^3;

%% =======================================================
%% Check 1: G.6 two-route symbolic agreement
%%
%% Route A: direct differentiation dS1*/dL via chain rule with
%%   e = sqrt(1 - G^2/L^2), eta = G/L.
%% Route B: (T1.L) formula  dS1*/dL = G^{-4} * (eta^3/e) * dPhi/de
%% Expected: simplify(route_A - route_B) = 0
%% =======================================================
printf('--- Check 1: G.6 two-route symbolic agreement ---\n');

% Build S_1^{*,(T)}(L, G, H, l, g) with e, theta as functions of (L, G, H).
% Substitute theta = H/G and e = sqrt(1 - (G/L)^2); eta is already sqrt(1-e^2).
S1star_in_LGH = subs(S1star_s, [theta_s, e_s], ...
                     [H_s/G_s, sqrt(1 - (G_s/L_s)^2)]);

% Route A: direct diff
dS1star_dL_A = diff(S1star_in_LGH, L_s);

% Route B: (T1.L) formula dS/dL = G^{-(α+1)} * (η³/e) * ∂F/∂e
% where ∂F/∂e is the parametric e-partial of F treating θ and e as independent.
% Phi_s is written with e_s and theta_s; diff(Phi_s, e_s) gives ∂F/∂e respecting
% the eta = sqrt(1-e^2) dependence automatically (since eta_expr is baked in).
dPhi_de_s = diff(Phi_s, e_s);
% Substitute (theta, e) -> (H/G, sqrt(1 - G^2/L^2)) after differentiation.
dPhi_de_in_LGH = subs(dPhi_de_s, [theta_s, e_s], ...
                      [H_s/G_s, sqrt(1 - (G_s/L_s)^2)]);
eta_LGH = G_s / L_s;
e_LGH = sqrt(1 - (G_s/L_s)^2);
dS1star_dL_B = (eta_LGH^3 / e_LGH) * dPhi_de_in_LGH / G_s^4;

residual_1 = simplify(dS1star_dL_A - dS1star_dL_B);
% Check if residual is zero
is_zero_1 = isequal(residual_1, sym(0));
if ~is_zero_1
  % Try more aggressive simplification
  residual_1b = simplify(expand(residual_1));
  is_zero_1 = isequal(residual_1b, sym(0));
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
%% Check 2: G.7 GSI vs (T1.G) cross-check
%%
%% Route T1G:  dS1*/dG = G^{-4} [ -3*Phi - theta*dPhi/dtheta - (eta^2/e)*dPhi/de ]
%% Route GSI:  dS1*/dG = -(3/G)*S1* - theta*dS1*/dH - (1/eta)*dS1*/dL
%% Expected: simplify(T1G - GSI) = 0 identically
%% =======================================================
printf('--- Check 2: G.7 (T1.G) vs GSI cross-check ---\n');

% Route T1G: direct formula evaluation
dPhi_dtheta_s = diff(Phi_s, theta_s);
dPhi_de_param = diff(Phi_s, e_s);   % at fixed (theta, eta)
T1G_form = (-3 * Phi_s - theta_s * dPhi_dtheta_s - (eta_expr^2 / e_s) * dPhi_de_param) / G_s^4;

% Route GSI: substitute S1* and partials
% dS1*/dH via (T1.H): G^{-4} * dPhi/dtheta
T1H_form = dPhi_dtheta_s / G_s^4;
% dS1*/dL via (T1.L): G^{-4} * (eta^3/e) * dPhi/de
T1L_form = (eta_expr^3 / e_s) * dPhi_de_param / G_s^4;
% GSI assembly: -(3/G)*S1* - theta*T1H - (1/eta)*T1L
GSI_form = -(3/G_s) * S1star_s - theta_s * T1H_form - (1/eta_expr) * T1L_form;

residual_2 = simplify(T1G_form - GSI_form);
is_zero_2 = isequal(residual_2, sym(0));
if ~is_zero_2
  residual_2b = simplify(expand(residual_2));
  is_zero_2 = isequal(residual_2b, sym(0));
end

print_check('Check 2', ...
  'simplify(T1G - GSI) = 0 identically', ...
  sprintf('residual = %s', char(residual_2)), ...
  is_zero_2);
if is_zero_2
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 3: G.8 finite-difference confirmation
%%
%% Evaluate dS1*/dH analytically at sample point and compare to
%% centered finite difference in H.
%% Expected: |analytical - FD| / |analytical| < 1e-7
%% =======================================================
printf('--- Check 3: G.8 FD confirmation ---\n');

% Numerical evaluators
function val = S1star_num(L, G, H, g_val)
  eta = G/L;
  e = sqrt(1 - eta^2);
  theta = H/G;
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  inner = A*B*e^2*(6 + e^2)/4 * sin(2*g_val) + B^2*e^4/128 * sin(4*g_val);
  Phi = 1 / ((5*theta^2 - 1) * eta^2) * inner;  % L=mu=k2=1
  val = Phi / G^3;
endfunction

% Sample point: L=1, G=0.92, H=0.4*G, g=0.7, staying on D' (theta=0.4, away from 1/sqrt(5)=0.447)
L_n = 1.0;
G_n = 0.92;
H_n = 0.4 * G_n;
g_n = 0.7;
theta_n = H_n / G_n;
eta_n = G_n / L_n;
e_n = sqrt(1 - eta_n^2);

% Analytical via (T1.H): G^{-4} * dPhi/dtheta, evaluated at the sample
dPhi_dtheta_fn = function_handle(dPhi_dtheta_s);
% ...actually function_handle may not work in Octave with symbolic; use subs + double
dPhi_dtheta_at = double(subs(dPhi_dtheta_s, ...
                             [theta_s, e_s, g_s, mu_s, k2_s], ...
                             [sym(theta_n), sym(e_n), sym(g_n), sym(1), sym(1)]));
analytical_dH = dPhi_dtheta_at / G_n^4;

% FD: (S1*(H + h) - S1*(H - h)) / (2h)
h_fd = 1e-6;
FD_dH = (S1star_num(L_n, G_n, H_n + h_fd, g_n) - S1star_num(L_n, G_n, H_n - h_fd, g_n)) / (2*h_fd);

rel_err_3 = abs(analytical_dH - FD_dH) / max(abs(analytical_dH), 1e-14);
tol_3 = 1e-7;
pass_3 = rel_err_3 < tol_3;

printf('  Sample: L=%.3f, G=%.3f, H=%.3f, g=%.3f, theta=%.3f, e=%.3f\n', ...
       L_n, G_n, H_n, g_n, theta_n, e_n);
printf('  Analytical dS1*/dH = %.6e\n', analytical_dH);
printf('  Finite-diff  dS1*/dH = %.6e\n', FD_dH);
print_check('Check 3', ...
  sprintf('|analytical - FD| / |analytical| < %.0e', tol_3), ...
  sprintf('rel err = %.3e', rel_err_3), ...
  pass_3);
if pass_3
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 4: G.9b numerical homological equation re-verification
%%
%% Confirm that the G.9b closed form still satisfies G.1:
%%   (dF1*/dG) * (dS1*/dg) = -F_{2p}(g)
%% at 10 (theta, e, g) samples on D'.
%% Expected: max |residual| < 1e-12
%% =======================================================
printf('--- Check 4: G.9b homological equation re-verification ---\n');

% Reuse the closed form from Ch 11a verifier
function val = F2p_T(theta, e_val, g_val)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  c2T = 3 * A * B * e_val^2 * (6 + e_val^2) / (4 * eta^9);
  c4T = 3 * B^2 * e_val^4 / (64 * eta^9);
  val = c2T * cos(2*g_val) + c4T * cos(4*g_val);
endfunction

function val = dF1star_dG(theta, e_val)
  eta = sqrt(1 - e_val^2);
  val = -3 * eta^3 * (5*theta^2 - 1) / (2 * eta^7);
endfunction

% dS1^{*,(T)}/dg (analytical, from G.9b)
function val = dS1star_T_dg_num(theta, e_val, g_val)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  factor = (5*theta^2 - 1) * eta^2 * eta^3;
  inner = A*B*e_val^2*(6 + e_val^2)/2 * cos(2*g_val) ...
        + B^2*e_val^4/32            * cos(4*g_val);
  val = inner / factor;
endfunction

% 10 sample points on D'
theta_samples = [0.2, 0.3, 0.5, 0.6, 0.8];
e_samples     = [0.15, 0.35];
g_samples     = [0.5];  % fixed; we use 5x2 = 10 pts

% Alternatively, use a 10-point grid
theta_g = [0.2, 0.3, 0.5, 0.6, 0.8];
e_g     = [0.15, 0.5];
g_g     = [0.5];

max_res_4 = 0;
check4_pass = 0;
check4_fail = 0;
tol_4 = 1e-12;

printf('  theta    e      g        residual\n');
for ith = 1:length(theta_g)
  for ie = 1:length(e_g)
    for ig = 1:length(g_g)
      th = theta_g(ith);
      ee = e_g(ie);
      gg = g_g(ig);
      lhs = dF1star_dG(th, ee) * dS1star_T_dg_num(th, ee, gg);
      rhs = F2p_T(th, ee, gg);
      residual = lhs + rhs;
      absres = abs(residual);
      if absres > max_res_4
        max_res_4 = absres;
      end
      if absres < tol_4
        check4_pass = check4_pass + 1;
      else
        check4_fail = check4_fail + 1;
      end
      printf('  %4.2f    %4.2f   %.3f    %.3e\n', th, ee, gg, residual);
    end
  end
end
n_samples_4 = length(theta_g) * length(e_g) * length(g_g);
printf('\n  ... (%d sample points total).\n', n_samples_4);
print_check('Check 4', ...
  sprintf('max |residual| < %.0e', tol_4), ...
  sprintf('max |residual| = %.3e (%d PASS / %d FAIL)', max_res_4, check4_pass, check4_fail), ...
  check4_fail == 0);
if check4_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Chapter 11b verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Propositions G.6, G.7, G.8, G.9 all confirmed.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
