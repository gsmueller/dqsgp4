%% derive_delta_q1_v2.m
%% Derive delta_q_1 = D(dS1/dL + dS1*/dL) symbolically.
%%
%% dS1/dL from Brouwer (1959) Eq. (20):
%%   dS1/dL = -(gamma2/e)*(L/G)*{ A*(a^2*eta^2/r^2+a/r+1)*sin(f)
%%             + B/2*[(-a^2*eta^2/r^2-a/r+1)*sin(2g+f) + (a^2*eta^2/r^2+a/r+1/3)*sin(2g+3f)] }
%%   where gamma2 = mu^2*k2/L^4, eta = G/L.
%%
%% BUT: we must NOT trust this expression blindly (OCR risk). We should
%% verify it by checking that it is consistent with dS1/dl (homological eq)
%% via the mixed partial: d(dS1/dL)/dl = d(dS1/dl)/dL.
%%
%% For now, proceed with the Brouwer expression and compare the FINAL result
%% D(dS1/dL) with BH61 Eq. (14). Agreement validates both.
%%
%% dS1*/dL from derive_S1star.m (Octave-computed):
%%   dS1*/dL = k2*mu^2*(-G^4+16*G^2*H^2-15*H^4)*sin(2g) / (8*G^3*L^3*(G^2-5*H^2))

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING delta_q_1 = D(dS1/dL + dS1*/dL)\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms f_s g_s e_s a_s r_s real;

theta_s = H_s/G_s;
eta_s = G_s/L_s;

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

gamma2 = mu^2*k2/L_s^4;

%% ================================================================
%% dS1/dL from Brouwer (1959) Eq. (20)
%%
%% Rewrite with explicit L, G, a, r, e:
%%   gamma2*L/(e*G) = mu^2*k2/(L^3*e*G)
%%   a^2*eta^2/r^2 = a^2*G^2/(L^2*r^2)
%%
%% dS1/dL = -(mu^2*k2/(L^3*e*G)) * {
%%   A*(a^2*G^2/(L^2*r^2) + a/r + 1)*sin(f)
%%   + B/2*[(-a^2*G^2/(L^2*r^2) - a/r + 1)*sin(2g+f)
%%          + (a^2*G^2/(L^2*r^2) + a/r + 1/3)*sin(2g+3f)]
%% }
%% ================================================================

printf('--- dS1/dL from Brouwer Eq. (20) ---\n\n');

% Let z = a^2*G^2/(L^2*r^2) = a^2*eta^2/r^2
z = a_s^2*G_s^2/(L_s^2*r_s^2);

dS1_dL = -(mu^2*k2/(L_s^3*e_s*G_s)) * ( ...
  A_coeff * (z + a_s/r_s + 1)*sin(f_s) + ...
  B_coeff/2 * ((-z - a_s/r_s + 1)*sin(2*g_s+f_s) + (z + a_s/r_s + sym(1)/3)*sin(2*g_s+3*f_s)) ...
);

printf('dS1/dL defined.\n\n');

%% ================================================================
%% dS1*/dL from derive_S1star.m (Octave-verified)
%%
%% dS1*/dL = mu^2*k2*(-G^4+16*G^2*H^2-15*H^4)*sin(2g) / (8*G^3*L^3*(G^2-5*H^2))
%%         = -(mu^2*k2/(G^3*L^3)) * (G^4-16*G^2*H^2+15*H^4)/(8*(G^2-5*H^2)) * sin(2g)
%%         = -(mu^2*k2/(G^3*L^3)) * (1-16*theta^2+15*theta^4)/(8*(1-5*theta^2)) * G^4/G^4
%% Hmm, let me just use the raw form.
%% ================================================================

printf('--- dS1*/dL ---\n\n');

dS1star_dL = mu^2*k2*(-G_s^4+16*G_s^2*H_s^2-15*H_s^4)*sin(2*g_s) / (8*G_s^3*L_s^3*(G_s^2-5*H_s^2));

printf('dS1*/dL defined.\n\n');

%% ================================================================
%% Apply D to (dS1/dL + dS1*/dL)
%%
%% For dS1*/dL: this has the form [function of G, H, theta] * sin(2g) / L^3.
%% D acts on 1/L^3 (gives factor 3*(2a/r-1)/L^3), on sin(2g) (gives
%% -4sinf*cos(2g)/e), and through product rule.
%%
%% For dS1/dL: this is the big one. Multiple terms with a/r, z, sin(f),
%% sin(2g+f), sin(2g+3f), all multiplied by 1/(L^3*e*G).
%%
%% Rather than decompose manually, let me compute D symbolically using
%% the chain rule formula:
%%   D(expr) = sum of [D(factor)] contributions via product rule
%%
%% For expressions involving a, r, e, f, g, L, G, H:
%%   D acts through DL, DG, DH and then through Da, De, Df, Dg
%%   via the chain rule.
%%
%% The cleanest approach: substitute a/r and all orbital relations,
%% express everything in terms of (e, f, g, theta, L, G) where
%% D(e) = -2(e+cosf), D(f) = 2sinf/e, D(g) = -2sinf/e,
%% D(L^n) = -n*L^n*(2a/r-1), D(G^n) = -n*G^n, D(theta) = 0.
%%
%% After applying D, substitute a/r = (1+e*cosf)/(1-e^2) and simplify.
%% Then compare with BH61 Eq. (14).
%%
%% This is a LOT of algebra. Let SymPy do it.
%% ================================================================

printf('--- Computing D(dS1/dL) via chain rule ---\n\n');
printf('This is a large computation. Using the chain rule:\n');
printf('  D = -p_1*d/dL - G*d/dG - H*d/dH + q_1*d/dl + q_2*d/dg\n');
printf('For dS1/dL, which does not depend on l or H:\n');
printf('  D(dS1/dL) = -p_1*d(dS1/dL)/dL - G*d(dS1/dL)/dG + q_2*d(dS1/dL)/dg\n\n');

% But d(dS1/dL)/dL involves d/dL of terms with a/r, e, f — all implicit
% functions of L. This requires the chain rule through a = L^2/mu and
% e = sqrt(1-G^2/L^2) and f = f(l,e).
%
% Since dS1/dL is expressed in terms of a, r, e, f, g, L, G, theta,
% and these are NOT independent (a = L^2/mu, eta = G/L, etc.),
% I cannot just differentiate w.r.t. L treating everything else as constant.
%
% CORRECT APPROACH: substitute all orbital relations, express dS1/dL
% entirely in terms of L, G, H, f, g (with a/r expressed via the orbit eq)
% and THEN apply D as a differential operator via the D-identities.
%
% D acts on any expression F(L,G,H,e,f,g,a,r,theta) as:
% DF = (dF/dL)*DL + (dF/dG)*DG + (dF/de)*De + (dF/df)*Df + (dF/dg)*Dg
%      + (dF/da)*Da + (dF/dr)*Dr + (dF/dH)*DH + (dF/dtheta)*Dtheta
%
% But a, e, r, f, theta all depend on L, G, H through the orbital relations.
% If we treat them as INDEPENDENT symbols and then substitute the D-values,
% that's valid as long as we account for all the cross-dependencies.
%
% Actually: the D-identities (Da, De, Df, Dg, Dr=0, Dtheta=0, DL, DG, DH)
% were derived from the chain rule through xi_k. They are self-consistent.
% So I CAN treat (a, e, f, g, r, L, G, H, theta) as independent symbols
% in the expression, apply D to each factor using the identities, and the
% result is correct.
%
% This is the KEY insight: the D-identities form a consistent system.
% I don't need to re-derive the chain rule — I just apply:
%   D(L^n) = -n*L^n*(2a/r-1)
%   D(G^n) = -n*G^n
%   D(e^n): use De = -2(e+cosf) and chain rule
%   D(sin(nf+mg)): use Df, Dg
%   D(a^n/r^m): Da and Dr=0
%   D(theta): = 0
% All via the product rule.

% So let me apply D to dS1/dL factor by factor.

% dS1/dL = -(mu^2*k2) * (1/(L^3*e*G)) * { ... }
% = -(mu^2*k2) * L^{-3} * e^{-1} * G^{-1} * { ... }
%
% Factor 1: L^{-3}, D(L^{-3}) = 3*L^{-3}*(2a/r-1)
% Factor 2: e^{-1}, D(e^{-1}) = -e^{-2}*De = -e^{-2}*(-2(e+cosf)) = 2(e+cosf)/e^2
% Factor 3: G^{-1}, D(G^{-1}) = 1/G
% Factor 4: the bracket B(a,r,e,f,g,theta) — contains sin(f), sin(2g+f), sin(2g+3f)
%           multiplied by functions of a/r and theta.

% This 4-factor product rule is:
% D(F1*F2*F3*F4) = DF1*F2*F3*F4 + F1*DF2*F3*F4 + F1*F2*DF3*F4 + F1*F2*F3*DF4

% Rather than do this by hand with 4 factors, let me have SymPy compute
% D(dS1/dL) by treating all variables as independent and replacing:
%   L -> L, with the rule D(L) = -L*(2a/r-1)  [but this only works for total D]
%
% Actually, the most reliable approach: write a function that applies D
% to a symbolic expression via the known identities.

% Define D as a function of symbolic expressions:
% For a monomial c * L^a * G^b * H^c * e^d * (a_s)^p / r^q * trig_stuff:
% D = sum of (d/d_var * D_var) for each variable.

% Let me define D(expr) = sum over all variables v of diff(expr, v)*Dv:
DL_val = -L_s*(2*a_s/r_s - 1);
DG_val = -G_s;
DH_val = -H_s;
De_val = -2*(e_s + cos(f_s));
Df_val = 2*sin(f_s)/e_s;
Dg_val = -2*sin(f_s)/e_s;
Da_val = -2*a_s*(2*a_s/r_s - 1);
Dr_val = sym(0);
Dtheta_val = sym(0);

% D(expr) treating all as independent:
D_op = @(expr) diff(expr, L_s)*DL_val + diff(expr, G_s)*DG_val + ...
               diff(expr, H_s)*DH_val + diff(expr, e_s)*De_val + ...
               diff(expr, f_s)*Df_val + diff(expr, g_s)*Dg_val + ...
               diff(expr, a_s)*Da_val + diff(expr, r_s)*Dr_val;

% Apply D to dS1/dL:
printf('Applying D operator to dS1/dL...\n');
D_dS1_dL = D_op(dS1_dL);
D_dS1_dL = expand(D_dS1_dL);
printf('D(dS1/dL) computed (expanded).\n\n');

% Apply D to dS1*/dL:
printf('Applying D operator to dS1*/dL...\n');
D_dS1star_dL = D_op(dS1star_dL);
D_dS1star_dL = expand(D_dS1star_dL);
printf('D(dS1*/dL) computed (expanded).\n\n');

% Total delta_q_1:
delta_q1 = D_dS1_dL + D_dS1star_dL;
delta_q1 = expand(delta_q1);

%% ================================================================
%% Compare with BH61 Eq. (14)
%%
%% Substitute a = L^2/mu, r = a*(1-e^2)/(1+e*cos(f)), G = L*sqrt(1-e^2)
%% to reduce everything to (L, e, f, g, theta).
%% ================================================================

printf('--- Substituting orbital relations ---\n\n');

% Substitute a = L^2/mu:
delta_q1_sub = subs(delta_q1, a_s, L_s^2/mu);
% Substitute r = L^2*(1-e^2)/(mu*(1+e*cos(f))):
delta_q1_sub = subs(delta_q1_sub, r_s, L_s^2*(1-e_s^2)/(mu*(1+e_s*cos(f_s))));
% Substitute G = L*sqrt(1-e^2):
delta_q1_sub = subs(delta_q1_sub, G_s, L_s*sqrt(1-e_s^2));
% Substitute H = L*sqrt(1-e^2)*theta:
delta_q1_sub = subs(delta_q1_sub, H_s, L_s*sqrt(1-e_s^2)*theta_s);

printf('Simplifying...\n');
delta_q1_simplified = simplify(delta_q1_sub);
printf('delta_q_1 (simplified):\n'); disp(delta_q1_simplified);

%% ================================================================
%% Build BH61 Eq. (14) delta_q_1 for comparison
%%
%% BH61 Eq. (14) delta_q_1 is very long (spans multiple lines).
%% Read the first few terms and compare.
%% For now, just print our result.
%% ================================================================

printf('\n============================================================\n');
printf('delta_q_1 derived symbolically. Compare with BH61 Eq. (14).\n');
printf('============================================================\n');
