%% derive_delta_p2_S1star.m
%% Compute D(dS1*/dg) — the S1* contribution to delta_p_2.
%%
%% From derive_S1star.m, dS1*/dg' was derived from the homological equation:
%%   dS1*/dg = G*gamma2*(L^2/G^2 - L^4/G^4)*(1-16*theta^2+15*theta^4)
%%             *(1-5*theta^2)^{-1} * cos(2g) / 8
%%
%% where gamma2 = mu^2*k2/L^4, theta = H/G.
%%
%% To apply D, use the chain rule:
%%   D(f) = -p_1*df/dL - G*df/dG - H*df/dH + q_1*df/dl + q_2*df/dg
%%
%% Since dS1*/dg depends on L, G, H, and g (not l), and q_3 = 0:
%%   D(dS1*/dg) = -p_1*d(dS1*/dg)/dL - G*d(dS1*/dg)/dG - H*d(dS1*/dg)/dH
%%                + q_2*d(dS1*/dg)/dg
%%
%% But wait — this requires q_2, which we haven't independently verified.
%% Also, d(dS1*/dg)/dg involves sin(2g), and q_2 involves sin(E)/...
%%
%% SIMPLER APPROACH: dS1*/dg depends on L, G, H, g only (no l, no f, no r, no E).
%% Therefore:
%%   - The L, G, H dependence: D acts through DL = -(p_1), DG = -G, DH = -H.
%%   - The g dependence: D(cos(2g)) = -sin(2g)*2*Dg = -sin(2g)*(-4sinf/e) = 4sinf*sin(2g)/e
%%
%% Factor structure of dS1*/dg:
%%   = [mu^2*k2 * (1-16*theta^2+15*theta^4) / (8*(1-5*theta^2))] * [1/(G*L^2) - 1/G^3] * cos(2g)
%%
%% where theta = H/G is degree-0, so (1-16*theta^2+15*theta^4)/(1-5*theta^2) is degree-0.
%% The velocity-dependent factors are: 1/(G*L^2) and 1/G^3 and cos(2g).
%%
%% D acts on [1/(G*L^2) - 1/G^3] * cos(2g) = product of two factors:
%%   Factor A: 1/(G*L^2) - 1/G^3
%%   Factor B: cos(2g)
%%
%% D(A*B) = D(A)*B + A*D(B)
%%
%% D(B) = D(cos(2g)) = 4*sinf*sin(2g)/e  [from Dg = -2sinf/e]
%%
%% D(A) = D(1/(G*L^2)) - D(1/G^3)
%%   D(1/(G*L^2)) = (4a/r-1)/(G*L^2)  [derived in derive_delta_p2.m]
%%   D(1/G^3) = 3/G^3                  [Euler, degree -3]
%%   D(A) = (4a/r-1)/(G*L^2) - 3/G^3

pkg load symbolic;

printf('============================================================\n');
printf('COMPUTING D(dS1*/dg) — S1* CONTRIBUTION TO delta_p_2\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms f_s g_s e_s a_s r_s real;

theta_s = H_s/G_s;
gamma2 = mu^2*k2/L_s^4;

% dS1*/dg from derive_S1star.m:
% = G*gamma2*(L^2/G^2 - L^4/G^4)*(1-16*theta^2+15*theta^4)/(8*(1-5*theta^2)) * cos(2g)
%
% Rewrite: G*gamma2*(L^2/G^2 - L^4/G^4) = G*(mu^2*k2/L^4)*(L^2/G^2 - L^4/G^4)
%   = mu^2*k2*(1/(G*L^2) - 1/G^3)
%
% So dS1*/dg = mu^2*k2*(1/(G*L^2) - 1/G^3)*C_theta*cos(2g)
% where C_theta = (1-16*theta^2+15*theta^4)/(8*(1-5*theta^2))

C_theta = (1 - 16*theta_s^2 + 15*theta_s^4) / (8*(1 - 5*theta_s^2));

% Factor A = 1/(G*L^2) - 1/G^3
% Factor B = cos(2g)
% dS1*/dg = mu^2*k2*C_theta * A * B

printf('dS1*/dg = mu^2*k2*C_theta * [1/(G*L^2) - 1/G^3] * cos(2g)\n');
printf('where C_theta = (1-16*theta^2+15*theta^4)/(8*(1-5*theta^2))\n\n');

% D(dS1*/dg) = mu^2*k2*C_theta * [D(A)*B + A*D(B)]
% (mu^2*k2 and C_theta are degree-0 constants w.r.t. D)

% D(A):
% D(1/(G*L^2)):
%   = D(G^{-1})*L^{-2} + G^{-1}*D(L^{-2})
%   = (1/G)*L^{-2} + G^{-1}*2*L^{-2}*(2a/r-1)
%   = 1/(G*L^2) * (1 + 2*(2a/r-1))
%   = 1/(G*L^2) * (4a/r - 1)

printf('D(1/(G*L^2)):\n');
printf('  D(G^{-1}) = (-1)*G^{-2}*DG = (-1)*G^{-2}*(-G) = 1/G\n');
printf('  D(L^{-2}) = (-2)*L^{-3}*DL = (-2)*L^{-3}*(-L*(2a/r-1)) = 2*L^{-2}*(2a/r-1)\n');
printf('  D(1/(G*L^2)) = L^{-2}/G + G^{-1}*2*L^{-2}*(2a/r-1)\n');
printf('               = 1/(G*L^2) * [1 + 2*(2a/r-1)] = (4a/r-1)/(G*L^2)\n\n');

% D(1/G^3) = 3/G^3  [Euler]
printf('D(1/G^3) = 3/G^3\n\n');

DA = (4*a_s/r_s - 1)/(G_s*L_s^2) - 3/G_s^3;
printf('D(A) = (4a/r-1)/(G*L^2) - 3/G^3\n\n');

% D(B) = D(cos(2g)) = -sin(2g)*2*Dg = -sin(2g)*2*(-2sinf/e) = 4sinf*sin(2g)/e
DB = 4*sin(f_s)*sin(2*g_s)/e_s;
printf('D(cos(2g)) = 4*sin(f)*sin(2g)/e\n\n');

% A = 1/(G*L^2) - 1/G^3
A_expr = 1/(G_s*L_s^2) - 1/G_s^3;

% Full D(dS1*/dg) = mu^2*k2*C_theta * [DA*cos(2g) + A*DB]
D_dS1star_dg = mu^2*k2*C_theta * (DA*cos(2*g_s) + A_expr*DB);
D_dS1star_dg = expand(D_dS1star_dg);

printf('D(dS1*/dg) (expanded):\n'); disp(D_dS1star_dg);

% Now compare with BH61 Eq. (14), delta_p_2 second line:
% (mu^2 k2/G^3)*e*(1-16*theta^2+15*theta^4)*(1-5*theta^2)^{-1}
%   * [e/8*cos(2g) + 1/2*cos(2g+f)]

BH61_dp2_S1star = mu^2*k2/G_s^3 * e_s * (1-16*theta_s^2+15*theta_s^4)*(1-5*theta_s^2)^(-1) ...
                  * (e_s/8*cos(2*g_s) + sym(1)/2*cos(2*g_s+f_s));

diff_S1star = simplify(expand(D_dS1star_dg - BH61_dp2_S1star));
printf('\nD(dS1*/dg) - BH61 second line = '); disp(diff_S1star);
printf('(Should be 0 if BH61 is correct)\n\n');

% If not zero, try to identify what the correct expression is.
% Factor out mu^2*k2/G^3 * C_coeff_expanded:
D_dS1star_dg_factored = simplify(D_dS1star_dg / (mu^2*k2));
printf('D(dS1*/dg) / (mu^2*k2) = '); disp(D_dS1star_dg_factored);

printf('\n============================================================\n');
