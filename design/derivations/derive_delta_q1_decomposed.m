%% derive_delta_q1_decomposed.m
%% Derive delta_q_1 = D(dS1/dL) + D(dS1*/dL) by decomposing into
%% individual harmonic terms and applying D to each separately.
%%
%% Key principle: keep a/r as a symbol throughout (do NOT substitute
%% the orbit equation). BH61 Eq. (14) is written in terms of a/r,
%% so the comparison is direct.
%%
%% dS1/dL from Brouwer (1959) Eq. (20) has three harmonic terms:
%%   Term_A:  ~ sin(f)
%%   Term_B1: ~ sin(2g+f)
%%   Term_B2: ~ sin(2g+3f)
%%
%% Each is processed independently: apply D, simplify, then compare
%% with the corresponding group in BH61 Eq. (14).

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING delta_q_1 BY HARMONIC DECOMPOSITION\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms f_s g_s e_s a_s r_s real;

theta_s = H_s/G_s;
eta_s = G_s/L_s;

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

% D-operator function: applies D to expr treating all variables as independent
DL_val = -L_s*(2*a_s/r_s - 1);
DG_val = -G_s;
DH_val = -H_s;
De_val = -2*(e_s + cos(f_s));
Df_val = 2*sin(f_s)/e_s;
Dg_val = -2*sin(f_s)/e_s;
Da_val = -2*a_s*(2*a_s/r_s - 1);
Dr_val = sym(0);

D_op = @(expr) diff(expr, L_s)*DL_val + diff(expr, G_s)*DG_val + ...
               diff(expr, H_s)*DH_val + diff(expr, e_s)*De_val + ...
               diff(expr, f_s)*Df_val + diff(expr, g_s)*Dg_val + ...
               diff(expr, a_s)*Da_val + diff(expr, r_s)*Dr_val;

%% ================================================================
%% dS1/dL decomposition
%%
%% dS1/dL = -(mu^2 k2/(L^3 e G)) * {
%%   A*(a^2*eta^2/r^2 + a/r + 1)*sin(f)                          [Term_A]
%%   + (B/2)*(-a^2*eta^2/r^2 - a/r + 1)*sin(2g+f)               [Term_B1]
%%   + (B/2)*(a^2*eta^2/r^2 + a/r + 1/3)*sin(2g+3f)             [Term_B2]
%% }
%%
%% where eta^2 = G^2/L^2, a^2*eta^2/r^2 = a^2*G^2/(L^2*r^2).
%% ================================================================

printf('=== PROCESSING Term_A: sin(f) group ===\n\n');

% Use z = a^2*G^2/(L^2*r^2) = a^2*eta^2/r^2 as shorthand
z = a_s^2*G_s^2/(L_s^2*r_s^2);
w = a_s/r_s;

Term_A = -(mu^2*k2/(L_s^3*e_s*G_s)) * A_coeff * (z + w + 1) * sin(f_s);

printf('Term_A = -(mu^2 k2/(L^3 e G)) * A * (z+w+1) * sin(f)\n');
printf('  where z = a^2*eta^2/r^2, w = a/r\n\n');

D_Term_A = D_op(Term_A);
D_Term_A = expand(D_Term_A);
D_Term_A = simplify(D_Term_A);
printf('D(Term_A) = '); disp(D_Term_A);

% Factor out mu^2*k2/(L^3*e*G) to see the structure:
D_A_reduced = simplify(D_Term_A / (-(mu^2*k2/(L_s^3*e_s*G_s))));
printf('D(Term_A) / [-(mu^2 k2/(L^3 e G))] = '); disp(D_A_reduced);

%% ================================================================
printf('\n=== PROCESSING Term_B1: sin(2g+f) group ===\n\n');

Term_B1 = -(mu^2*k2/(L_s^3*e_s*G_s)) * (B_coeff/2) * (-z - w + 1) * sin(2*g_s + f_s);

printf('Term_B1 = -(mu^2 k2/(L^3 e G)) * (B/2) * (-z-w+1) * sin(2g+f)\n\n');

D_Term_B1 = D_op(Term_B1);
D_Term_B1 = expand(D_Term_B1);
D_Term_B1 = simplify(D_Term_B1);
printf('D(Term_B1) = '); disp(D_Term_B1);

D_B1_reduced = simplify(D_Term_B1 / (-(mu^2*k2/(L_s^3*e_s*G_s))));
printf('D(Term_B1) / [-(mu^2 k2/(L^3 e G))] = '); disp(D_B1_reduced);

%% ================================================================
printf('\n=== PROCESSING Term_B2: sin(2g+3f) group ===\n\n');

Term_B2 = -(mu^2*k2/(L_s^3*e_s*G_s)) * (B_coeff/2) * (z + w + sym(1)/3) * sin(2*g_s + 3*f_s);

printf('Term_B2 = -(mu^2 k2/(L^3 e G)) * (B/2) * (z+w+1/3) * sin(2g+3f)\n\n');

D_Term_B2 = D_op(Term_B2);
D_Term_B2 = expand(D_Term_B2);
D_Term_B2 = simplify(D_Term_B2);
printf('D(Term_B2) = '); disp(D_Term_B2);

D_B2_reduced = simplify(D_Term_B2 / (-(mu^2*k2/(L_s^3*e_s*G_s))));
printf('D(Term_B2) / [-(mu^2 k2/(L^3 e G))] = '); disp(D_B2_reduced);

%% ================================================================
printf('\n=== PROCESSING S1* contribution: D(dS1*/dL) ===\n\n');

% dS1*/dL from derive_S1star.m:
% = mu^2*k2*(-G^4+16*G^2*H^2-15*H^4)*sin(2g) / (8*G^3*L^3*(G^2-5*H^2))

dS1star_dL = mu^2*k2*(-G_s^4+16*G_s^2*H_s^2-15*H_s^4)*sin(2*g_s) / ...
             (8*G_s^3*L_s^3*(G_s^2-5*H_s^2));

D_S1star_dL = D_op(dS1star_dL);
D_S1star_dL = expand(D_S1star_dL);
D_S1star_dL = simplify(D_S1star_dL);
printf('D(dS1*/dL) = '); disp(D_S1star_dL);

D_S1star_reduced = simplify(D_S1star_dL / (mu^2*k2/(L_s^3*e_s*G_s)));
printf('D(dS1*/dL) / [mu^2 k2/(L^3 e G)] = '); disp(D_S1star_reduced);

%% ================================================================
printf('\n=== ASSEMBLING delta_q_1 ===\n\n');

delta_q1 = D_Term_A + D_Term_B1 + D_Term_B2 + D_S1star_dL;

% Don't try to simplify the whole thing. Instead, compare piece by piece.
% BH61 Eq. (14) for delta_q_1 (reading the structure, not the coefficients):
%
% Line 1 (prefactor mu^2 k2/(e^2 L^3 G)):
%   A-group * sin(2f) + B-group * [sin(2g) + sin(2g+4f)]
%
% Line 2 (prefactor mu^2 k2/(e L^3 G) * a/r):
%   A-group * sin(f) + B-group * [sin(2g+f) + sin(2g+3f)]
%
% Line 3 (prefactor mu^2 k2/(e L^3 G)):
%   mixed theta terms * [e*sin(2g) + sin(2g+f) - sin(2g-f)]
%
% Plus the continuation on page 197 with more terms.

printf('delta_q_1 assembled. Individual pieces printed above.\n');
printf('Compare each harmonic group with BH61 Eq. (14) separately.\n');

printf('\n============================================================\n');
