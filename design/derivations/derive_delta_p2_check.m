%% derive_delta_p2_check.m
%% Check: does the S1 contribution to delta_p_2 match BH61 by itself?

pkg load symbolic;

syms mu k2 a_s e_s f_s g_s r_s theta_s L_s G_s real;

B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

% D(dS1/dg) from the main script:
% dS1/dg = (mu^2 k2/G^3) * B * [cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]
% D acts on 1/G^3 (gives factor 3) and on the bracket contents.
% The bracket's D-value was computed in the main script.

% Alternatively, the S1-only part of BH61 Eq.(14) for dp2 is:
% (mu^2 k2/G^3)*B*[1/3*cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]
%
% Note: BH61 has 1/3*cos(2g+2f), not cos(2g+2f). The factor went from 1/2 to 1/3.
% This suggests D acting on the 1/2 factor of cos(2g+2f) in S1 changed it.
% Actually no — dS1/dg differentiates the sin terms in S1 to cos terms.
% dS1/dg = (mu^2 k2/G^3)*B*[cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]
% The prefactors are 1, e, e/3 (from differentiating sin(2g+2f), etc.)

% BH61's S1-only dp2 result has 1/3 in front of cos(2g+2f).
% So the D operator changed the coefficient from 1 to 1/3.
% From D(dS1/dg): the cos(2g+2f) term gets D(1/G^3)*cos(2g+2f) = 3/G^3*cos(2g+2f)
% from the 1/G^3 factor, plus D(cos(2g+2f)) = 0. So the coefficient of cos(2g+2f)
% goes from 1 to 3*1 = 3? That's 3, not 1/3.
%
% Wait, D(dS1/dg) = (mu^2 k2 B / G^3) * (3*F_bracket + D(F_bracket)).
% The coefficient of cos(2g+2f) in 3*F_bracket is 3*1 = 3.
% D(F_bracket) has no cos(2g+2f) term (since D(cos(2g+2f)) = 0).
% So D(dS1/dg) has 3*cos(2g+2f) coefficient (times mu^2 k2 B/G^3).
%
% But BH61 has 1/3*cos(2g+2f). Factor of 9 discrepancy? That can't be right.
%
% Let me re-examine. BH61 Eq.(14) for delta_p_2:
% Line 1: (mu^2 k2/G^3)*(3/2-3/2*theta^2)*[1/3*cos(2g+2f)+e*cos(2g+f)+(e/3)*cos(2g+3f)]
%
% My result from D acting on dS1/dg should give delta_p_2(S1 part).
% Let me compute just the S1 part and compare with line 1 of BH61.

% D(dS1/dg):
% dS1/dg = (mu^2 k2/G^3) * B * bracket
% bracket = cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)
%
% D acts via product rule on [1/G^3] * [bracket]:
% D(1/G^3 * bracket) = D(1/G^3)*bracket + (1/G^3)*D(bracket)
%                     = (3/G^3)*bracket + (1/G^3)*D(bracket)

% D(bracket):
% D(cos(2g+2f)) = 0  [since D(f+g) = 0]
% D(e*cos(2g+f)) = De*cos(2g+f) + e*(-sin(2g+f))*D(2g+f)
%   D(2g+f) = 2*Dg + Df = -4sinf/e + 2sinf/e = -2sinf/e
%   = -2(e+cosf)*cos(2g+f) + e*sin(2g+f)*2sinf/e
%   = -2(e+cosf)*cos(2g+f) + 2sinf*sin(2g+f)
%   Use product-to-sum: sinf*sin(2g+f) = (1/2)[cos(2g) - cos(2g+2f)]
%   and (e+cosf)*cos(2g+f) = e*cos(2g+f) + cosf*cos(2g+f)
%   cosf*cos(2g+f) = (1/2)[cos(2g) + cos(2g+2f)]
%   So: -2e*cos(2g+f) - 2*(1/2)*cos(2g) - 2*(1/2)*cos(2g+2f) + 2*(1/2)*cos(2g) - 2*(1/2)*cos(2g+2f)
%     = -2e*cos(2g+f) - cos(2g) - cos(2g+2f) + cos(2g) - cos(2g+2f)
%     = -2e*cos(2g+f) - 2cos(2g+2f)

% D((e/3)*cos(2g+3f)) = (De/3)*cos(2g+3f) + (e/3)*(-sin(2g+3f))*D(2g+3f)
%   D(2g+3f) = 2Dg + 3Df = -4sinf/e + 6sinf/e = 2sinf/e
%   = -(2/3)(e+cosf)*cos(2g+3f) - (2sinf/3)*e*sin(2g+3f)/e
%   Wait: (e/3)*(-sin(2g+3f))*(2sinf/e) = -(2sinf/3)*sin(2g+3f)
%   And De/3 = -(2/3)*(e+cosf)
%   So: -(2/3)(e+cosf)*cos(2g+3f) - (2sinf/3)*sin(2g+3f)
%   Use: (e+cosf)*cos(2g+3f) = e*cos(2g+3f) + cosf*cos(2g+3f)
%   cosf*cos(2g+3f) = (1/2)[cos(2g+2f) + cos(2g+4f)]
%   sinf*sin(2g+3f) = (1/2)[cos(2g+2f) - cos(2g+4f)]
%   So: -(2e/3)*cos(2g+3f) - (1/3)[cos(2g+2f)+cos(2g+4f)] - (1/3)[cos(2g+2f)-cos(2g+4f)]
%     = -(2e/3)*cos(2g+3f) - (2/3)*cos(2g+2f)

% Total D(bracket):
% 0 + [-2e*cos(2g+f) - 2cos(2g+2f)] + [-(2e/3)*cos(2g+3f) - (2/3)*cos(2g+2f)]
% = -2e*cos(2g+f) - (2+2/3)*cos(2g+2f) - (2e/3)*cos(2g+3f)
% = -2e*cos(2g+f) - (8/3)*cos(2g+2f) - (2e/3)*cos(2g+3f)

% So D(dS1/dg) = (mu^2 k2 B / G^3) * [3*bracket + D(bracket)]
% = (mu^2 k2 B / G^3) * [3*cos(2g+2f) + 3e*cos(2g+f) + e*cos(2g+3f)
%    -2e*cos(2g+f) - (8/3)*cos(2g+2f) - (2e/3)*cos(2g+3f)]
% = (mu^2 k2 B / G^3) * [(3-8/3)*cos(2g+2f) + (3-2)*e*cos(2g+f) + (1-2/3)*e*cos(2g+3f)]
% = (mu^2 k2 B / G^3) * [(1/3)*cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]

printf('D(dS1/dg) S1-only contribution:\n');
printf('= (mu^2 k2 B / G^3) * [(1/3)*cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)]\n\n');

printf('BH61 Eq.(14), delta_p_2 first line:\n');
printf('= (mu^2 k2/G^3)*(3/2-3/2*theta^2)*[1/3*cos(2g+2f)+e*cos(2g+f)+(e/3)*cos(2g+3f)]\n\n');

printf('MATCH. The S1 contribution to delta_p_2 is correct.\n\n');

% Verify with SymPy:
bracket = cos(2*g_s+2*f_s) + e_s*cos(2*g_s+f_s) + (e_s/3)*cos(2*g_s+3*f_s);
D_bracket = -2*e_s*cos(2*g_s+f_s) - sym(8)/3*cos(2*g_s+2*f_s) - (2*e_s/3)*cos(2*g_s+3*f_s);
combined = 3*bracket + D_bracket;
combined = expand(combined);

% Expected: (1/3)*cos(2g+2f) + e*cos(2g+f) + (e/3)*cos(2g+3f)
expected = sym(1)/3*cos(2*g_s+2*f_s) + e_s*cos(2*g_s+f_s) + (e_s/3)*cos(2*g_s+3*f_s);
diff = simplify(combined - expected);
printf('3*bracket + D(bracket) - expected = '); disp(diff);
printf('(Should be 0)\n');

printf('\nThe discrepancy in delta_p_2 is entirely in the S1* contribution.\n');
printf('This is the territory of open discrepancy D004.\n');
