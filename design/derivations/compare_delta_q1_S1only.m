%% compare_delta_q1_S1only.m
%% Check if the S1-only contribution to delta_q_1 matches BH61 lines 358-359.

pkg load symbolic;

printf('============================================================\n');
printf('COMPARING delta_q_1 (S1-ONLY) WITH BH61 LINES 358-359\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms f_s g_s e_s a_s r_s theta_s real;

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

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

z = a_s^2*G_s^2/(L_s^2*r_s^2);
w = a_s/r_s;

% dS1/dL from Brouwer Eq. (20):
dS1_dL = -(mu^2*k2/(L_s^3*e_s*G_s)) * ( ...
  A_coeff * (z + w + 1)*sin(f_s) + ...
  B_coeff/2 * ((-z - w + 1)*sin(2*g_s+f_s) + (z + w + sym(1)/3)*sin(2*g_s+3*f_s)) ...
);

printf('Computing D(dS1/dL)...\n');
D_dS1_dL = D_op(dS1_dL);
D_dS1_dL = expand(D_dS1_dL);
printf('Done.\n\n');

% BH61 Eq. (14) lines 358-359 ONLY (S1 contribution to delta_q_1):
% Line 358: mu^2 k2/(e^2 L^3 G) * { 2A*(z+w+1)*sin(2f) + B*[(-z-w+1)*sin(2g)+(z+w+1/3)*sin(2g+4f)] }
% Line 359: mu^2 k2/(e L^3 G)*(a/r)*{ 2A*(z+w+4)*sin(f) + B*[(-z-w+2)*sin(2g+f)+(z+w+2)*sin(2g+3f)] }
%
% Wait: BH61 writes (-1+3*theta^2), not A. And A = -1/2+3/2*theta^2 = (1/2)*(-1+3*theta^2).
% So (-1+3*theta^2) = 2*A.
% And (3/2-3/2*theta^2) = B.

BH61_S1_line1 = mu^2*k2/(e_s^2*L_s^3*G_s) * ( ...
  2*A_coeff*(z+w+1)*sin(2*f_s) + ...
  B_coeff*((-z-w+1)*sin(2*g_s) + (z+w+sym(1)/3)*sin(2*g_s+4*f_s)) ...
);

BH61_S1_line2 = mu^2*k2/(e_s*L_s^3*G_s)*(a_s/r_s) * ( ...
  2*A_coeff*(z+w+4)*sin(f_s) + ...
  B_coeff*((-z-w+2)*sin(2*g_s+f_s) + (z+w+2)*sin(2*g_s+3*f_s)) ...
);

BH61_S1_only = BH61_S1_line1 + BH61_S1_line2;

% Substitute orbital relations:
subs_vars = {a_s, r_s, G_s, H_s};
subs_vals = {L_s^2/mu, L_s^2*(1-e_s^2)/(mu*(1+e_s*cos(f_s))), ...
             L_s*sqrt(1-e_s^2), L_s*sqrt(1-e_s^2)*theta_s};

our_sub = D_dS1_dL;
bh61_sub = BH61_S1_only;

for i = 1:4
  our_sub = subs(our_sub, subs_vars{i}, subs_vals{i});
  bh61_sub = subs(bh61_sub, subs_vars{i}, subs_vals{i});
end

printf('Simplifying S1-only difference...\n');
diff_S1 = simplify(expand(our_sub - bh61_sub));
printf('\nD(dS1/dL) - BH61 lines 358-359 = '); disp(diff_S1);
printf('(Should be 0 if S1 part matches)\n');

printf('\n============================================================\n');
