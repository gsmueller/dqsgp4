%% compare_delta_q1_BH61.m
%% Compare our D(dS1/dL + dS1*/dL) with BH61 Eq. (14) delta_q_1.
%%
%% Strategy: compute both expressions, substitute all orbital relations
%% (a = L^2/mu, r = a(1-e^2)/(1+ecosf), G = L*sqrt(1-e^2), H = G*theta)
%% to reduce to (L, e, f, g, theta), then check the difference is zero.
%%
%% This is the same approach that worked for delta_p_2 (resolve_D004.m).

pkg load symbolic;

printf('============================================================\n');
printf('COMPARING delta_q_1 WITH BH61 Eq. (14)\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms f_s g_s e_s a_s r_s theta_s real;

% D-operator function
DL_val = -L_s*(2*a_s/r_s - 1);
DG_val = -G_s;
DH_val = -H_s;
De_val = -2*(e_s + cos(f_s));
Df_val = 2*sin(f_s)/e_s;
Dg_val = -2*sin(f_s)/e_s;
Da_val = -2*a_s*(2*a_s/r_s - 1);

D_op = @(expr) diff(expr, L_s)*DL_val + diff(expr, G_s)*DG_val + ...
               diff(expr, H_s)*DH_val + diff(expr, e_s)*De_val + ...
               diff(expr, f_s)*Df_val + diff(expr, g_s)*Dg_val + ...
               diff(expr, a_s)*Da_val + diff(expr, r_s)*sym(0);

%% ================================================================
%% OUR RESULT: D(dS1/dL + dS1*/dL)
%% ================================================================

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

% dS1/dL from Brouwer Eq. (20):
z = a_s^2*G_s^2/(L_s^2*r_s^2);
w = a_s/r_s;

dS1_dL = -(mu^2*k2/(L_s^3*e_s*G_s)) * ( ...
  A_coeff * (z + w + 1)*sin(f_s) + ...
  B_coeff/2 * ((-z - w + 1)*sin(2*g_s+f_s) + (z + w + sym(1)/3)*sin(2*g_s+3*f_s)) ...
);

% dS1*/dL from derive_S1star.m:
dS1star_dL = mu^2*k2*(-G_s^4+16*G_s^2*H_s^2-15*H_s^4)*sin(2*g_s) / ...
             (8*G_s^3*L_s^3*(G_s^2-5*H_s^2));

printf('Computing D(dS1/dL)...\n');
D_dS1_dL = D_op(dS1_dL);
D_dS1_dL = expand(D_dS1_dL);
printf('Done.\n');

printf('Computing D(dS1*/dL)...\n');
D_dS1star_dL = D_op(dS1star_dL);
D_dS1star_dL = expand(D_dS1star_dL);
printf('Done.\n');

our_delta_q1 = D_dS1_dL + D_dS1star_dL;

%% ================================================================
%% BH61 Eq. (14) delta_q_1
%%
%% BH61 uses eta = G/L = sqrt(1-e^2) and theta = H/G.
%% The expression spans lines 358-363 + 369-372 of the repaired document.
%%
%% I read the STRUCTURE but let SymPy build it.
%% ================================================================

printf('\nBuilding BH61 Eq. (14) delta_q_1...\n');

% Line 358 (prefactor mu^2 k2/(e^2 L^3 G)):
% {(-1+3*theta^2)*(a^2*eta^2/r^2 + a/r + 1)*sin(2f)
%  +(3/2-3/2*theta^2)*[(-a^2*eta^2/r^2 - a/r + 1)*sin(2g)
%                       +(a^2*eta^2/r^2 + a/r + 1/3)*sin(2g+4f)]}

BH61_line1 = mu^2*k2/(e_s^2*L_s^3*G_s) * ( ...
  (-1+3*theta_s^2)*(z+w+1)*sin(2*f_s) + ...
  (sym(3)/2-sym(3)/2*theta_s^2)*((-z-w+1)*sin(2*g_s) + (z+w+sym(1)/3)*sin(2*g_s+4*f_s)) ...
);

% Line 359 (prefactor mu^2 k2/(e L^3 G) * a/r):
% {(-1+3*theta^2)*(a^2*eta^2/r^2 + a/r + 4)*sin(f)
%  +(3/2-3/2*theta^2)*[(-a^2*eta^2/r^2 - a/r + 2)*sin(2g+f)
%                       +(a^2*eta^2/r^2 + a/r + 2)*sin(2g+3f)]}

BH61_line2 = mu^2*k2/(e_s*L_s^3*G_s) * (a_s/r_s) * ( ...
  (-1+3*theta_s^2)*(z+w+4)*sin(f_s) + ...
  (sym(3)/2-sym(3)/2*theta_s^2)*((-z-w+2)*sin(2*g_s+f_s) + (z+w+2)*sin(2*g_s+3*f_s)) ...
);

% Line 360 (prefactor mu^2 k2/(e L^3 G)):
% [1/4*(1-11*theta^2) - 10*theta^4*(1-5*theta^2)^{-1}]
%   * [(1-3*a/r)*e*sin(2g) + sin(2g+f) - sin(2g-f)]

C_mixed = sym(1)/4*(1-11*theta_s^2) - 10*theta_s^4*(1-5*theta_s^2)^(-1);

BH61_line3 = mu^2*k2/(e_s*L_s^3*G_s) * C_mixed * ( ...
  (1-3*w)*e_s*sin(2*g_s) + sin(2*g_s+f_s) - sin(2*g_s-f_s) ...
);

BH61_delta_q1 = BH61_line1 + BH61_line2 + BH61_line3;

%% ================================================================
%% COMPARE: substitute orbital relations and check difference
%% ================================================================

printf('Substituting orbital relations...\n');

% Substitute: a = L^2/mu, r = L^2*(1-e^2)/(mu*(1+e*cos(f))),
% G = L*sqrt(1-e^2), H = G*theta = L*sqrt(1-e^2)*theta

subs_list_from = {a_s, r_s, G_s, H_s};
subs_list_to = {L_s^2/mu, L_s^2*(1-e_s^2)/(mu*(1+e_s*cos(f_s))), ...
                L_s*sqrt(1-e_s^2), L_s*sqrt(1-e_s^2)*theta_s};

our_sub = our_delta_q1;
bh61_sub = BH61_delta_q1;

for i = 1:4
  our_sub = subs(our_sub, subs_list_from{i}, subs_list_to{i});
  bh61_sub = subs(bh61_sub, subs_list_from{i}, subs_list_to{i});
end

printf('Simplifying difference...\n');
diff_q1 = simplify(expand(our_sub - bh61_sub));

printf('\ndelta_q_1 (ours) - BH61 Eq.(14) = '); disp(diff_q1);
printf('(Should be 0 if they match)\n');

printf('\n============================================================\n');
