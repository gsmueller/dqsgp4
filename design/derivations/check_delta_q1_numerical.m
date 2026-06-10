%% check_delta_q1_numerical.m
%% Numerically evaluate (our delta_q_1 S1 part) - (BH61 lines 358-359)
%% at self-consistent orbital test points.

pkg load symbolic;

syms mu k2 real positive;
syms L G H real positive;
syms a_s r_s real positive;
syms l_sym g e f real;

theta = H/G;
A_expr = -sym(1)/2 + sym(3)/2 * theta^2;
B_expr =  sym(3)/2 - sym(3)/2 * theta^2;

% S1 pieces:
P0 = (mu^2*k2/G^3) * A_expr*(f - l_sym + e*sin(f));
Q1 = (mu^2*k2/G^3) * B_expr*(sym(1)/2)*sin(2*g+2*f);
Q2 = (mu^2*k2/G^3) * B_expr*(e/2)*sin(2*g+f);
Q3 = (mu^2*k2/G^3) * B_expr*(e/6)*sin(2*g+3*f);

% Kepler chain:
de_dL   = G^2/(L^3*e);
df_dL   = G^2 * sin(f)*(2+e*cos(f)) / (L^3*e*(1-e^2));
dX_dL = @(X) diff(X, e)*de_dL + diff(X, f)*df_dL;

dP0_dL = dX_dL(P0);
dQ1_dL = dX_dL(Q1);
dQ2_dL = dX_dL(Q2);
dQ3_dL = dX_dL(Q3);

% D operator:
DL_val = -L*(2*a_s/r_s - 1);
DG_val = -G;
DH_val = -H;
De_val = -2*(e + cos(f));
Df_val =  2*sin(f)/e;
Dg_val = -2*sin(f)/e;
Da_val = -2*a_s*(2*a_s/r_s - 1);
D_op = @(X) diff(X, L)*DL_val + diff(X, G)*DG_val + diff(X, H)*DH_val + ...
            diff(X, e)*De_val + diff(X, f)*Df_val + diff(X, g)*Dg_val + ...
            diff(X, a_s)*Da_val;

Dq1_P0 = D_op(dP0_dL);
Dq1_Q1 = D_op(dQ1_dL);
Dq1_Q2 = D_op(dQ2_dL);
Dq1_Q3 = D_op(dQ3_dL);

delta_q1_S1 = Dq1_P0 + Dq1_Q1 + Dq1_Q2 + Dq1_Q3;

% BH61 lines 358-359:
z = a_s^2*G^2/(L^2*r_s^2);
w = a_s/r_s;

BH61_line1 = mu^2*k2/(e^2*L^3*G) * ( ...
  2*A_expr*(z+w+1)*sin(2*f) + ...
  B_expr*((-z-w+1)*sin(2*g) + (z+w+sym(1)/3)*sin(2*g+4*f)) );

BH61_line2 = mu^2*k2/(e*L^3*G)*(a_s/r_s) * ( ...
  2*A_expr*(z+w+4)*sin(f) + ...
  B_expr*((-z-w+2)*sin(2*g+f) + (z+w+2)*sin(2*g+3*f)) );

BH61_S1 = BH61_line1 + BH61_line2;

diff_expr = delta_q1_S1 - BH61_S1;

% Self-consistent test points:
function v = mktest(mu_v, k2_v, L_v, e_v, theta_v, f_v, g_v)
  v.mu = mu_v; v.k2 = k2_v; v.L = L_v; v.e = e_v;
  v.G = L_v * sqrt(1 - e_v^2);
  v.H = v.G * theta_v;
  v.a = L_v^2/mu_v;
  v.r = v.a * (1 - e_v^2)/(1 + e_v*cos(f_v));
  v.f = f_v; v.g = g_v;
endfunction

tests = { mktest(1,    0.001, 1.2, 0.3, 0.8, 0.7, 1.5), ...
          mktest(2,    0.005, 1.5, 0.1, 0.3, 2.1, 0.7), ...
          mktest(3.98e5, -2.63e-3, 7000, 0.05, 0.6, 1.1, 2.3) };

printf('Comparing delta_q_1 (S1 part) vs BH61 lines 358-359:\n');
args = {mu, k2, L, G, H, e, f, g, a_s, r_s};
for i = 1:length(tests)
  v = tests{i};
  vals = {v.mu, v.k2, v.L, v.G, v.H, v.e, v.f, v.g, v.a, v.r};
  ours = double(subs(delta_q1_S1, args, vals));
  bh61 = double(subs(BH61_S1,     args, vals));
  d    = ours - bh61;
  printf('  Test %d: ours=%14.8g  bh61=%14.8g  diff=%14.8g  rel=%g\n', ...
    i, ours, bh61, d, abs(d)/max(abs(ours),abs(bh61)+eps));
end
