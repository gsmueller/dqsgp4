% verify_ch11d_secular_rates.m
%
% Chapter 11d verifier for second-order secular-rate partials
% (Proposition G.10 T-component) per CH_SECULAR_RATES_PLAN.md Part C.
%
% Five Checks:
%   Check 1 -- G.10 (T1.L): dF2**(T)/dL closed form vs direct SymPy diff.
%   Check 2 -- G.10 (T1.G): dF2**(T)/dG closed form vs direct SymPy diff.
%   Check 3 -- G.10 (T1.H): dF2**(T)/dH closed form vs direct SymPy diff.
%   Check 4 -- GSI cross-check: (T1.G) = -10/G F_2**(T) - theta dF/dH - (1/eta) dF/dL.
%   Check 5 -- Critical-inclination regularity: all three partials finite at theta^2 = 1/5.
%
% Reporting convention: each Check prints Expected / Observed / Status.

pkg load symbolic;

printf('=========================================================\n');
printf('Chapter 11d verifier for second-order secular-rate partials\n');
printf('(G.10 T-component, CH_SECULAR_RATES_PLAN.md Part C)\n');
printf('5 Checks: (T1.L), (T1.G), (T1.H), GSI, critical inclination\n');
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
% Symbolic setup: F_2^{**,(T)}(theta, e, g=avg, eta) from ch10e F.15
%   c_0^{(T)} = (3 mu^6 k_2^2 eta)/(2 G^10) * [(A^2 + B^2/2)(8 + 24e^2 + 3e^4)/8 - A^2 eta^3]
% Independent symbols: L, G, H, mu, k_2.
% Derived: eta = G/L, theta = H/G, e = sqrt(1 - G^2/L^2) = sqrt(1 - eta^2).
% --------------------------------------------------------------
syms L_s G_s H_s mu_s k2_s positive;
syms e_s theta_s positive;
% We work with two formulations:
%   Formulation 1: F_2^{**,(T)} expressed as function of (L, G, H) via Delaunay identities.
%     Used for direct SymPy diff (routeA).
%   Formulation 2: F_2^{**,(T)} expressed as F-factor Phi(theta, e, eta)
%     with eta = sqrt(1 - e^2) kept symbolic; (T1.*) applied formally.
%     Used for closed-form (routeB).
%
% For symbolic-simplification robustness, use Formulation 2 as the
% ground truth, and check that routeB from §2 matches routeA.

% ----- Formulation 1: F_2^{**,(T)}(L, G, H) -----
eta_LGH = G_s / L_s;
e_LGH = sqrt(1 - eta_LGH^2);
theta_LGH = H_s / G_s;

A_LGH = (3*theta_LGH^2 - 1) / 2;
B_LGH = 3*(1 - theta_LGH^2) / 2;
Q_LGH = A_LGH^2 + B_LGH^2/2;
R_LGH = A_LGH^2;
P_LGH = 8 + 24*e_LGH^2 + 3*e_LGH^4;

F2ssT_LGH = (3*mu_s^6*k2_s^2*eta_LGH / (2*G_s^10)) * ...
            (Q_LGH * P_LGH / 8 - R_LGH * eta_LGH^3);

%% =======================================================
%% Check 1: (T1.L): dF_2^{**,(T)}/dL
%%
%% Route A: direct SymPy diff of F_2^{**,(T)}(L, G, H) with respect to L.
%% Route B: closed form from ch11d §2.2:
%%   dF/dL = (3 mu^6 k_2^2 eta^2)/(2 G^{11}) * [5 Q(8 - 12e^2 - 3e^4)/8 + 4 R eta^3]
%% =======================================================
printf('--- Check 1: G.10 (T1.L) dF_2^{**,(T)}/dL ---\n');

routeA_dL = diff(F2ssT_LGH, L_s);

% Route B: express in (L, G, H) via substitution
eta_sub = G_s / L_s;
e_sub = sqrt(1 - (G_s/L_s)^2);
theta_sub = H_s / G_s;
A_sub = (3*theta_sub^2 - 1)/2;
B_sub = 3*(1 - theta_sub^2)/2;
Q_sub = A_sub^2 + B_sub^2/2;
R_sub = A_sub^2;

routeB_dL = (3*mu_s^6*k2_s^2 * eta_sub^2 / (2 * G_s^11)) * ...
            (5*Q_sub*(8 - 12*e_sub^2 - 3*e_sub^4)/8 + 4*R_sub * eta_sub^3);

residual_1 = simplify(routeA_dL - routeB_dL);
is_zero_1 = isequal(residual_1, sym(0));
if ~is_zero_1
  residual_1 = simplify(expand(residual_1));
  is_zero_1 = isequal(residual_1, sym(0));
  if ~is_zero_1
    residual_1 = simplify(factor(residual_1));
    is_zero_1 = isequal(residual_1, sym(0));
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
%% Check 2: (T1.G): dF_2^{**,(T)}/dG
%%
%% Route A: direct SymPy diff.
%% Route B: closed form from ch11d §2.3:
%%   dF/dG = (3 mu^6 k_2^2)/(2 G^{11}) * [
%%            (eta/8)(-(10Q + theta Q')P - 5Q(8 - 12e^2 - 3e^4))
%%            + eta^4 (6R + theta R')
%%          ]
%% =======================================================
printf('--- Check 2: G.10 (T1.G) dF_2^{**,(T)}/dG ---\n');

routeA_dG = diff(F2ssT_LGH, G_s);

% Route B: compute Q'(theta) and R'(theta) analytically, substitute
Qprime_sub = (3*theta_sub/2) * (9*theta_sub^2 - 5);
Rprime_sub = 3*theta_sub * (3*theta_sub^2 - 1);
P_sub = 8 + 24*e_sub^2 + 3*e_sub^4;

routeB_dG_inner = (eta_sub/8) * (-(10*Q_sub + theta_sub*Qprime_sub) * P_sub ...
                                 - 5*Q_sub*(8 - 12*e_sub^2 - 3*e_sub^4)) ...
                 + eta_sub^4 * (6*R_sub + theta_sub*Rprime_sub);
routeB_dG = (3*mu_s^6*k2_s^2 / (2 * G_s^11)) * routeB_dG_inner;

residual_2 = simplify(routeA_dG - routeB_dG);
is_zero_2 = isequal(residual_2, sym(0));
if ~is_zero_2
  residual_2 = simplify(expand(residual_2));
  is_zero_2 = isequal(residual_2, sym(0));
  if ~is_zero_2
    residual_2 = simplify(factor(residual_2));
    is_zero_2 = isequal(residual_2, sym(0));
  end
end

print_check('Check 2', ...
  'simplify(routeA - routeB) = 0 identically', ...
  sprintf('residual = %s', char(residual_2)), ...
  is_zero_2);
if is_zero_2
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 3: (T1.H): dF_2^{**,(T)}/dH
%%
%% Route A: direct SymPy diff.
%% Route B: closed form from ch11d §2.1:
%%   dF/dH = (9 mu^6 k_2^2 theta)/(2 G^{11}) * [(9 theta^2 - 5)(8 + 24 e^2 + 3 e^4) eta/16 - (3 theta^2 - 1) eta^4]
%% =======================================================
printf('--- Check 3: G.10 (T1.H) dF_2^{**,(T)}/dH ---\n');

routeA_dH = diff(F2ssT_LGH, H_s);

routeB_dH = (9 * mu_s^6 * k2_s^2 * theta_sub / (2 * G_s^11)) * ...
            ((9*theta_sub^2 - 5) * P_sub * eta_sub / 16 - (3*theta_sub^2 - 1) * eta_sub^4);

residual_3 = simplify(routeA_dH - routeB_dH);
is_zero_3 = isequal(residual_3, sym(0));
if ~is_zero_3
  residual_3 = simplify(expand(residual_3));
  is_zero_3 = isequal(residual_3, sym(0));
  if ~is_zero_3
    residual_3 = simplify(factor(residual_3));
    is_zero_3 = isequal(residual_3, sym(0));
  end
end

print_check('Check 3', ...
  'simplify(routeA - routeB) = 0 identically', ...
  sprintf('residual = %s', char(residual_3)), ...
  is_zero_3);
if is_zero_3
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 4: GSI cross-check
%%
%%   GSI at alpha=10:
%%     dF/dG = -(10/G) * F - theta * dF/dH - (1/eta) * dF/dL
%%
%% Both sides evaluated on F_2^{**,(T)}; expected: simplify(LHS - RHS) = 0.
%% =======================================================
printf('--- Check 4: GSI cross-check at alpha=10 ---\n');

GSI_RHS = -(10/G_s) * F2ssT_LGH - theta_sub * routeA_dH - (1/eta_sub) * routeA_dL;

residual_4 = simplify(routeA_dG - GSI_RHS);
is_zero_4 = isequal(residual_4, sym(0));
if ~is_zero_4
  residual_4 = simplify(expand(residual_4));
  is_zero_4 = isequal(residual_4, sym(0));
  if ~is_zero_4
    residual_4 = simplify(factor(residual_4));
    is_zero_4 = isequal(residual_4, sym(0));
  end
end

print_check('Check 4', ...
  'simplify((T1.G) form - GSI form) = 0 identically', ...
  sprintf('residual = %s', char(residual_4)), ...
  is_zero_4);
if is_zero_4
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 5: Critical-inclination regularity
%%
%% Evaluate all three partials at theta = 1/sqrt(5), several (e, L, G sub) samples,
%% confirm all values are finite. No (5 theta^2 - 1) denominator should appear.
%% =======================================================
printf('--- Check 5: Critical-inclination regularity (theta = 1/sqrt(5)) ---\n');

% Numerical evaluators
function val = F2ssT_dH_num(theta, e_val, G_val, mu, k2)
  eta = sqrt(1 - e_val^2);
  P = 8 + 24*e_val^2 + 3*e_val^4;
  val = (9 * mu^6 * k2^2 * theta / (2 * G_val^11)) * ...
        ((9*theta^2 - 5) * P * eta / 16 - (3*theta^2 - 1) * eta^4);
endfunction

function val = F2ssT_dL_num(theta, e_val, G_val, mu, k2)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  Q = A^2 + B^2/2;
  R = A^2;
  val = (3 * mu^6 * k2^2 * eta^2 / (2 * G_val^11)) * ...
        (5*Q*(8 - 12*e_val^2 - 3*e_val^4)/8 + 4*R * eta^3);
endfunction

function val = F2ssT_dG_num(theta, e_val, G_val, mu, k2)
  eta = sqrt(1 - e_val^2);
  A = (3*theta^2 - 1)/2;
  B = 3*(1 - theta^2)/2;
  Q = A^2 + B^2/2;
  R = A^2;
  Qp = (3*theta/2) * (9*theta^2 - 5);
  Rp = 3*theta * (3*theta^2 - 1);
  P = 8 + 24*e_val^2 + 3*e_val^4;
  inner = (eta/8) * (-(10*Q + theta*Qp)*P - 5*Q*(8 - 12*e_val^2 - 3*e_val^4)) ...
        + eta^4 * (6*R + theta*Rp);
  val = (3 * mu^6 * k2^2 / (2 * G_val^11)) * inner;
endfunction

theta_crit = 1/sqrt(5);
e_grid_crit = [0.05, 0.2, 0.5, 0.7];
G_test = 1.0;
mu_test = 1.0;
k2_test = 1.0;

printf('  theta = 1/sqrt(5) = %.6f  (exact critical inclination)\n', theta_crit);
printf('  e        dF/dL             dF/dG             dF/dH             all finite?\n');
check5_pass = 0;
check5_fail = 0;
for ie = 1:length(e_grid_crit)
  ee = e_grid_crit(ie);
  dL = F2ssT_dL_num(theta_crit, ee, G_test, mu_test, k2_test);
  dG = F2ssT_dG_num(theta_crit, ee, G_test, mu_test, k2_test);
  dH = F2ssT_dH_num(theta_crit, ee, G_test, mu_test, k2_test);
  all_finite = isfinite(dL) && isfinite(dG) && isfinite(dH);
  if all_finite
    check5_pass = check5_pass + 1;
    status = 'yes';
  else
    check5_fail = check5_fail + 1;
    status = 'NO';
  end
  printf('  %4.2f     %14.6e  %14.6e  %14.6e  %s\n', ee, dL, dG, dH, status);
end

print_check('Check 5', ...
  sprintf('all 3 partials finite at theta^2 = 1/5 over %d e-samples', length(e_grid_crit)), ...
  sprintf('%d PASS / %d FAIL', check5_pass, check5_fail), ...
  check5_fail == 0);
if check5_fail == 0
  n_pass = n_pass + 1;
else
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
if n_fail == 0
  printf('Chapter 11d verifier: ALL %d checks PASSED.\n', n_pass);
  printf('Proposition G.10 (T-component) confirmed.\n');
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
