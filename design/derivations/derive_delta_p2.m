%% derive_delta_p2.m
%% Derive delta_p_2 = D(dS1/dg + dS1*/dg) symbolically.
%%
%% From Eq. (11a): delta_p_2 = D(dS1/dl_2) = D(dS1/dg + dS1*/dg).
%%
%% dS1/dg from Brouwer (1959) Eq. (16):
%%   dS1/dg = (mu^2 k2/G^3) * B * [cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]
%%
%% dS1*/dg from Brouwer (1959) Eq. (33)/line 466:
%%   dS1*/dg = G*gamma2*(L^2/G^2 - L^4/G^4) *
%%             [1/8*(1-16*theta^2+15*theta^4)] * (1-5*theta^2)^{-1} * cos(2g)
%%   where gamma2 = mu^2*k2/L^4.

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING delta_p_2 = D(dS1/dg + dS1*/dg)\n');
printf('============================================================\n\n');

syms mu k2 a_s e_s f_s g_s r_s theta_s L_s G_s real;

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

%% ================================================================
%% Part 1: D(dS1/dg)
%%
%% dS1/dg = (mu^2 k2 / G^3) * B * [cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]
%%
%% Factor structure: (mu^2 k2 * B) is constant+degree-0. 1/G^3 is degree -3.
%% The bracket [...] contains trig functions of (2g+2f), (2g+f), (2g+3f)
%% and factors of e (degree 0).
%%
%% D(2g+2f) = 2*Dg + 2*Df = -4sinf/e + 4sinf/e = 0
%% D(2g+f) = 2*Dg + Df = -4sinf/e + 2sinf/e = -2sinf/e
%% D(2g+3f) = 2*Dg + 3*Df = -4sinf/e + 6sinf/e = 2sinf/e
%%
%% D(cos(2g+2f)) = 0
%% D(cos(2g+f)) = -sin(2g+f)*D(2g+f) = -sin(2g+f)*(-2sinf/e) = 2sin(2g+f)*sinf/e
%% D(cos(2g+3f)) = -sin(2g+3f)*D(2g+3f) = -sin(2g+3f)*(2sinf/e) = -2sin(2g+3f)*sinf/e
%%
%% D(e*cos(2g+f)) = De*cos(2g+f) + e*D(cos(2g+f))
%%   = -2(e+cosf)*cos(2g+f) + e*2sin(2g+f)*sinf/e
%%   = -2(e+cosf)*cos(2g+f) + 2sinf*sin(2g+f)
%%
%% D((e/3)*cos(2g+3f)) = (De/3)*cos(2g+3f) + (e/3)*D(cos(2g+3f))
%%   = (-2(e+cosf)/3)*cos(2g+3f) + (e/3)*(-2sin(2g+3f)*sinf/e)
%%   = (-2(e+cosf)/3)*cos(2g+3f) - (2sinf/3)*sin(2g+3f)
%% ================================================================

printf('--- Part 1: D(dS1/dg) ---\n\n');

% The bracket F_bracket = cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)
F_bracket = cos(2*g_s+2*f_s) + e_s*cos(2*g_s+f_s) + (e_s/3)*cos(2*g_s+3*f_s);

% dS1/dg = (mu^2 k2 / G^3) * B * F_bracket
% D(dS1/dg) = (mu^2 k2 * B) * [D(1/G^3) * F_bracket + (1/G^3) * D(F_bracket)]
%           = (mu^2 k2 * B) * [3/G^3 * F_bracket + (1/G^3) * D(F_bracket)]
%           = (mu^2 k2 * B / G^3) * [3*F_bracket + D(F_bracket)]

% Compute D(F_bracket) term by term:
% D(cos(2g+2f)) = 0
D_cos_2g2f = sym(0);

% D(e*cos(2g+f)):
D_e_cos_2gf = -2*(e_s+cos(f_s))*cos(2*g_s+f_s) + 2*sin(f_s)*sin(2*g_s+f_s);

% D((e/3)*cos(2g+3f)):
D_e3_cos_2g3f = -2*(e_s+cos(f_s))/3*cos(2*g_s+3*f_s) - 2*sin(f_s)/3*sin(2*g_s+3*f_s);

D_F_bracket = D_cos_2g2f + D_e_cos_2gf + D_e3_cos_2g3f;
D_F_bracket = expand(D_F_bracket);

printf('D(F_bracket) = '); disp(D_F_bracket);

% Full D(dS1/dg):
D_dS1_dg = mu^2*k2*B_coeff/G_s^3 * (3*F_bracket + D_F_bracket);
D_dS1_dg = expand(D_dS1_dg);

printf('D(dS1/dg) (expanded) = \n'); disp(D_dS1_dg);

% Try to collect into cos(2g+2f), cos(2g+f), cos(2g+3f) terms:
% Use trig identities to simplify.
% sin(f)*sin(2g+f) = (1/2)[cos(2g) - cos(2g+2f)]  ... product-to-sum
% sin(f)*sin(2g+3f) = (1/2)[cos(2g+2f) - cos(2g+4f)]

% Let SymPy try to simplify with trigsimp:
D_dS1_dg_simplified = simplify(D_dS1_dg);
printf('D(dS1/dg) (simplified) = \n'); disp(D_dS1_dg_simplified);

%% ================================================================
%% Part 2: D(dS1*/dg)
%%
%% dS1*/dg from Brouwer (1959):
%%   dS1*/dg = G*gamma2*(L^2/G^2 - L^4/G^4) *
%%             (1/8)*(1-16*theta^2+15*theta^4)*(1-5*theta^2)^{-1} * cos(2g)
%%   where gamma2 = mu^2*k2/L^4.
%%
%% Substitute gamma2:
%%   dS1*/dg = (mu^2 k2 / L^4) * G * (L^2/G^2 - L^4/G^4) *
%%             (1/8)*(1-16*theta^2+15*theta^4)*(1-5*theta^2)^{-1} * cos(2g)
%%
%%   = (mu^2 k2) * (1/G - L^2/G^3) *   ... wait, let me expand:
%%     G*(L^2/G^2 - L^4/G^4)/L^4 = G*L^2/(G^2*L^4) - G*L^4/(G^4*L^4)
%%                                = 1/(G*L^2) - 1/G^3
%%
%% So dS1*/dg = mu^2 k2 * [1/(G*L^2) - 1/G^3] * C_coeff * cos(2g)
%% where C_coeff = (1/8)*(1-16*theta^2+15*theta^4)*(1-5*theta^2)^{-1}
%%
%% Factor analysis:
%%   mu^2*k2: constant
%%   C_coeff: degree-0 (depends on theta = H/G, degree 0)
%%   cos(2g): D(cos(2g)) = -sin(2g)*2*Dg = -sin(2g)*(-4sinf/e) = 4sinf*sin(2g)/e
%%   1/(G*L^2): need D(G^{-1}*L^{-2})
%%   1/G^3: D(G^{-3}) = 3/G^3
%%
%% D(G^{-1}*L^{-2}):
%%   = D(G^{-1})*L^{-2} + G^{-1}*D(L^{-2})
%%   = (1/G)*L^{-2} + G^{-1}*2*L^{-2}*(2a/r-1)
%%       [D(G^{-1}) = 1/G, D(L^{-2}) = 2*L^{-2}*(2a/r-1)]
%%   Wait: D(L^{-2}) = (-2)*L^{-3}*DL = (-2)*L^{-3}*(-L*(2a/r-1)) = 2*L^{-2}*(2a/r-1).
%%   And D(G^{-1}) = (-1)*G^{-2}*DG = (-1)*G^{-2}*(-G) = 1/G. Yes.
%%   So D(1/(G*L^2)) = L^{-2}/G + 2*(2a/r-1)/(G*L^2)
%%                    = 1/(G*L^2) * (1 + 2*(2a/r-1))
%%                    = 1/(G*L^2) * (4a/r - 1)
%% ================================================================

printf('\n--- Part 2: D(dS1*/dg) ---\n\n');

C_coeff = sym(1)/8 * (1 - 16*theta_s^2 + 15*theta_s^4) * (1 - 5*theta_s^2)^(-1);

% dS1*/dg = mu^2 k2 * [1/(G*L^2) - 1/G^3] * C_coeff * cos(2g)
% Split into two terms:
% U1 = mu^2 k2 * C_coeff * cos(2g) / (G*L^2)
% U2 = -mu^2 k2 * C_coeff * cos(2g) / G^3

% D(U1) = mu^2 k2 * C_coeff * [D(cos(2g))/(G*L^2) + cos(2g)*D(1/(G*L^2))]
%   D(cos(2g)) = -sin(2g)*2*Dg = -sin(2g)*2*(-2sinf/e) = 4sinf*sin(2g)/e
%   D(1/(G*L^2)) = (4a/r-1)/(G*L^2)  [computed above]

D_cos2g = 4*sin(f_s)*sin(2*g_s)/e_s;

D_U1 = mu^2*k2*C_coeff*(D_cos2g/(G_s*L_s^2) + cos(2*g_s)*(4*a_s/r_s-1)/(G_s*L_s^2));

% D(U2) = -mu^2 k2 * C_coeff * [D(cos(2g))/G^3 + cos(2g)*D(1/G^3)]
%   D(1/G^3) = 3/G^3
D_U2 = -mu^2*k2*C_coeff*(D_cos2g/G_s^3 + cos(2*g_s)*3/G_s^3);

D_dS1star_dg = D_U1 + D_U2;
D_dS1star_dg = expand(D_dS1star_dg);

printf('D(dS1*/dg) = \n'); disp(D_dS1star_dg);

% Factor mu^2*k2/G^3 and simplify:
D_dS1star_dg_factored = simplify(D_dS1star_dg / (mu^2*k2/G_s^3));
printf('D(dS1*/dg) / (mu^2*k2/G^3) = \n'); disp(D_dS1star_dg_factored);

%% ================================================================
%% Part 3: Assemble delta_p_2 and compare with BH61 Eq. (14)
%% ================================================================

printf('\n--- Assembling delta_p_2 ---\n\n');

delta_p2 = D_dS1_dg + D_dS1star_dg;
delta_p2 = expand(delta_p2);

printf('delta_p_2 (full, expanded) = \n'); disp(delta_p2);

% BH61 Eq. (14), delta_p_2:
% delta_p_2 = (mu^2 k2/G^3)*(3/2-3/2*theta^2)*[1/3*cos(2g+2f)+e*cos(2g+f)+(e/3)*cos(2g+3f)]
%   + (mu^2 k2/G^3)*e*(1-16*theta^2+15*theta^4)*(1-5*theta^2)^{-1}
%     * [1/8*e*cos(2g) + 1/2*cos(2g+f)]

BH61_dp2_line1 = mu^2*k2/G_s^3 * B_coeff * ...
  (sym(1)/3*cos(2*g_s+2*f_s) + e_s*cos(2*g_s+f_s) + e_s/3*cos(2*g_s+3*f_s));

BH61_dp2_line2 = mu^2*k2/G_s^3 * e_s * (1-16*theta_s^2+15*theta_s^4)*(1-5*theta_s^2)^(-1) * ...
  (e_s/8*cos(2*g_s) + sym(1)/2*cos(2*g_s+f_s));

BH61_dp2 = BH61_dp2_line1 + BH61_dp2_line2;

diff_p2 = simplify(expand(delta_p2 - BH61_dp2));
printf('delta_p_2 (ours) - BH61 Eq.(14) = ');
disp(diff_p2);
printf('(Should be 0 if they match)\n');

printf('\n============================================================\n');
