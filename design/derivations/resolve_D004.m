%% resolve_D004.m
%% Resolve the D004 discrepancy by computing D(dS1*/dg) in harmonic form
%% and comparing coefficient by coefficient with BH61 Eq. (14).
%%
%% Our exact D(dS1*/dg) was computed in derive_delta_p2_S1star.m.
%% Here we convert it to the BH61 harmonic basis: cos(2g), cos(2g+f), cos(2g-f).
%%
%% The conversion requires:
%%   a/r = (1+e*cos(f))/(1-e^2)
%%   sin(f)*sin(2g) = (1/2)[cos(2g-f) - cos(2g+f)]

pkg load symbolic;

printf('============================================================\n');
printf('RESOLVING D004: D(dS1*/dg) in harmonic form\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms f_s g_s e_s a_s r_s real;

theta_s = H_s/G_s;

% dS1*/dg = mu^2*k2 * C_theta * [1/(G*L^2) - 1/G^3] * cos(2g)
% where C_theta = (1-16*theta^2+15*theta^4) / (8*(1-5*theta^2))
%
% D(dS1*/dg) = mu^2*k2*C_theta * [DA*cos(2g) + A*DB]
% where A = 1/(G*L^2) - 1/G^3
%       DA = (4a/r-1)/(G*L^2) - 3/G^3
%       DB = D(cos(2g)) = 4*sin(f)*sin(2g)/e

C_theta = (1 - 16*theta_s^2 + 15*theta_s^4) / (8*(1 - 5*theta_s^2));

A_expr = 1/(G_s*L_s^2) - 1/G_s^3;
DA_expr = (4*a_s/r_s - 1)/(G_s*L_s^2) - 3/G_s^3;
DB_expr = 4*sin(f_s)*sin(2*g_s)/e_s;

D_dS1star_dg = mu^2*k2*C_theta * (DA_expr*cos(2*g_s) + A_expr*DB_expr);

%% ================================================================
%% Step 1: Substitute a/r = (1+e*cos(f))/(1-e^2)
%% and 1-e^2 = G^2/L^2, so a/r = L^2*(1+e*cos(f))/G^2
%% ================================================================

printf('--- Step 1: Substitute a/r ---\n\n');

% a = L^2/mu, r = a*(1-e^2)/(1+e*cos(f)), so a/r = (1+e*cos(f))/(1-e^2)
% In Delaunay: 1-e^2 = G^2/L^2, so a/r = L^2*(1+e*cos(f))/G^2

D_sub = subs(D_dS1star_dg, a_s/r_s, L_s^2*(1+e_s*cos(f_s))/G_s^2);
% Also need to handle 4*a/r explicitly:
D_sub = subs(D_sub, a_s, L_s^2/mu);
D_sub = subs(D_sub, r_s, L_s^2*(1-e_s^2)/(mu*(1+e_s*cos(f_s))));
D_sub = simplify(D_sub);

printf('After substituting a/r:\n'); disp(D_sub);

%% ================================================================
%% Step 2: Convert sin(f)*sin(2g) to harmonics
%% sin(f)*sin(2g) = (1/2)[cos(2g-f) - cos(2g+f)]
%% ================================================================

printf('\n--- Step 2: Convert to harmonics ---\n\n');

% Replace sin(f)*sin(2g) → (1/2)[cos(2g-f) - cos(2g+f)]
% SymPy might not do this automatically. Let me expand trigonometrically.
D_expanded = expand(D_sub);

% Use trigexpand to get sum-of-products form
% Actually let's use a substitution approach.
% Factor out mu^2*k2/G^3 and the C_theta factor, then examine what's left.

prefactor = mu^2*k2/G_s^3 * (1-16*theta_s^2+15*theta_s^4)/(1-5*theta_s^2);
D_reduced = simplify(D_dS1star_dg / prefactor);

printf('D(dS1*/dg) / [mu^2*k2/G^3 * (1-16th^2+15th^4)/(1-5th^2)] =\n');
disp(D_reduced);

% Now substitute a/r in the reduced expression:
D_reduced_sub = subs(D_reduced, a_s/r_s, L_s^2*(1+e_s*cos(f_s))/G_s^2);
D_reduced_sub = subs(D_reduced_sub, a_s, L_s^2/mu);
D_reduced_sub = subs(D_reduced_sub, r_s, L_s^2*(1-e_s^2)/(mu*(1+e_s*cos(f_s))));
D_reduced_sub = simplify(D_reduced_sub);

printf('After a/r substitution:\n'); disp(D_reduced_sub);

% Substitute 1-e^2 = G^2/L^2 (so L^2 = G^2/(1-e^2)):
% Actually, keep e_s as is. Let's see the harmonic structure.

% Express in terms of e only (eliminate L, G via eta = G/L):
syms eta_s real positive;
D_reduced_eta = subs(D_reduced_sub, G_s, L_s*eta_s);
D_reduced_eta = simplify(D_reduced_eta);
printf('In terms of eta = G/L:\n'); disp(D_reduced_eta);

% Now also use eta^2 = 1-e^2:
D_reduced_e = subs(D_reduced_eta, eta_s, sqrt(1-e_s^2));
D_reduced_e = simplify(D_reduced_e);
printf('In terms of e:\n'); disp(D_reduced_e);

%% ================================================================
%% Step 3: Compare with BH61's claimed form
%%
%% BH61 claims the S1* part of delta_p_2 is:
%%   (mu^2 k2/G^3) * e * (1-16*theta^2+15*theta^4)*(1-5*theta^2)^{-1}
%%     * [e/8 * cos(2g) + 1/2 * cos(2g+f)]
%%
%% Factor out the common prefactor to compare:
%%   BH61_reduced = e * [e/8 * cos(2g) + 1/2 * cos(2g+f)]
%%                = e^2/8 * cos(2g) + e/2 * cos(2g+f)
%% ================================================================

printf('\n--- Step 3: Compare with BH61 ---\n\n');

BH61_reduced = e_s^2/8 * cos(2*g_s) + e_s/2 * cos(2*g_s + f_s);
printf('BH61 reduced form: e^2/8*cos(2g) + e/2*cos(2g+f)\n');
printf('  = '); disp(BH61_reduced);

% Our reduced form should equal BH61_reduced if correct.
% D_reduced_e should be expressible as BH61_reduced.
% Let me check the difference:

% Actually I realize my factoring above might not be clean.
% Let me just compute the full thing and compare directly.

% Full BH61 S1* contribution:
BH61_S1star = mu^2*k2/G_s^3 * e_s * (1-16*theta_s^2+15*theta_s^4)*(1-5*theta_s^2)^(-1) ...
              * (e_s/8*cos(2*g_s) + sym(1)/2*cos(2*g_s+f_s));

% Our full D(dS1*/dg) - substitute a/r and simplify:
ours_full = subs(D_dS1star_dg, a_s, L_s^2/mu);
ours_full = subs(ours_full, r_s, L_s^2*(1-e_s^2)/(mu*(1+e_s*cos(f_s))));
ours_full = simplify(ours_full);

% Use G = L*sqrt(1-e^2):
ours_full_e = subs(ours_full, G_s, L_s*sqrt(1-e_s^2));
ours_full_e = simplify(ours_full_e);

BH61_e = subs(BH61_S1star, G_s, L_s*sqrt(1-e_s^2));
BH61_e = simplify(BH61_e);

diff_e = simplify(ours_full_e - BH61_e);
printf('Our result - BH61 (all in terms of e, theta, L, f, g):\n');
disp(diff_e);

% Factor out common factors to see the residual:
if ~isequal(diff_e, sym(0))
  diff_factored = simplify(diff_e / (mu^2*k2/L_s^3));
  printf('Residual / (mu^2*k2/L^3):\n');
  disp(diff_factored);

  % Try trigexpand:
  diff_trig = expand(diff_e);
  printf('Residual expanded:\n');
  disp(diff_trig);
end

printf('\n============================================================\n');
printf('If nonzero: BH61 Eq. (14) has an error in the S1* contribution to dp2.\n');
printf('Our derivation is ground truth (SymPy, from homological equation).\n');
printf('============================================================\n');
