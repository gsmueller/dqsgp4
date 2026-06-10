%% derive_S1star.m
%% Derive S1* from scratch via the long-period homological equation.
%%
%% The long-period homological equation (Brouwer 1959, Eq. 33):
%%   (dF1*/dG'') * (dS1*/dg') + F_{2p}* = 0
%%
%% Therefore:
%%   dS1*/dg' = -F_{2p}* / (dF1*/dG'')
%%
%% Then integrate:
%%   S1* = integral of (dS1*/dg') dg'
%%
%% Inputs (all from Brouwer 1959, derived in earlier sections):
%%   F1* = mu^4 k2 A / (L'^3 G'^3)     [Eq. 13, the orbit-averaged perturbation]
%%     where A = -1/2 + 3/2*H^2/G^2
%%   F_{2p}* = the cos(2g') part of F2* [Eq. 29, second bracket]
%%
%% F2* was derived by Brouwer through the second-order perturbation equation (Eq. 11).
%% We take F_{2p}* as given from Eq. (29) — it is the result of a long computation
%% involving products of S1 derivatives and F1 derivatives. Deriving F2* from scratch
%% would require re-deriving all of Brouwer's Section 4 (pp. 20-40 of the paper).
%%
%% Instead, we verify our final S1* against the CONSEQUENCES: the delta_p and delta_q
%% expressions must satisfy dimensional checks, limiting cases, and cross-consistency.

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING S1* FROM THE LONG-PERIOD HOMOLOGICAL EQUATION\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms Lp Gpp H real positive;  % L', G'', H (Brouwer's notation)
syms gp real;                  % g'

% gamma2 = mu^2 * k2 / L'^4
gamma2 = mu^2 * k2 / Lp^4;

% theta = H/G'' (= cos I)
theta = H/Gpp;

%% ================================================================
%% Step 1: Compute dF1*/dG''
%%
%% F1* = mu^4 k2 / (L'^3 G''^3) * A
%% where A = -1/2 + 3/2*H^2/G''^2 = -1/2 + 3/2*theta^2
%%
%% F1* = mu^4 k2 / (L'^3 G''^3) * (-1/2 + 3/2*H^2/G''^2)
%%     = mu^4 k2 / L'^3 * (-1/(2*G''^3) + 3*H^2/(2*G''^5))
%%
%% dF1*/dG'' = mu^4 k2 / L'^3 * (3/(2*G''^4) - 15*H^2/(2*G''^6))
%%           = (3/2) * mu^4 k2 / (L'^3 * G''^4) * (1 - 5*H^2/G''^2)
%%           = (3/2) * mu^4 k2 / (L'^3 * G''^4) * (1 - 5*theta^2)
%% ================================================================

printf('--- Step 1: dF1*/dG'''' ---\n\n');

A_coeff = -sym(1)/2 + sym(3)/2 * H^2/Gpp^2;
F1star = mu^4 * k2 * A_coeff / (Lp^3 * Gpp^3);

dF1star_dGpp = diff(F1star, Gpp);
dF1star_dGpp = simplify(dF1star_dGpp);
printf('dF1*/dG'''' = '); disp(dF1star_dGpp);

% Verify it matches (3/2)*mu^4*k2/(L'^3*G''^4)*(1-5*theta^2):
expected_dF1 = sym(3)/2 * mu^4 * k2 / (Lp^3 * Gpp^4) * (1 - 5*H^2/Gpp^2);
check1 = simplify(dF1star_dGpp - expected_dF1);
printf('Check: dF1*/dG'''' - (3/2)*mu^4*k2*(1-5*theta^2)/(L''^3*G''''^4) = '); disp(check1);
printf('\n');

%% ================================================================
%% Step 2: F_{2p}* — the g'-dependent part of F2*
%%
%% From Brouwer (1959) Eq. (29), F2* has two parts: F_{2s}* and F_{2p}*.
%% F_{2p}* is the cos(2g') term:
%%
%% F_{2p}* = (mu^6 k2^2 / L'^{10}) * [-3/16*(L'^5/G''^5 - L'^7/G''^7)
%%            * (1 - 16*H^2/G''^2 + 15*H^4/G''^4)] * cos(2g')
%%
%% Simplify the prefactor:
%% mu^6 k2^2 / L'^{10} * L'^5/G''^5 = mu^6 k2^2 / (L'^5 * G''^5)
%% mu^6 k2^2 / L'^{10} * L'^7/G''^7 = mu^6 k2^2 / (L'^3 * G''^7)
%%
%% So F_{2p}* = -3/16 * mu^6 k2^2 * (1/(L'^5*G''^5) - 1/(L'^3*G''^7))
%%              * (1 - 16*theta^2 + 15*theta^4) * cos(2g')
%% ================================================================

printf('--- Step 2: F_{2p}* ---\n\n');

F2p_star = mu^6*k2^2/Lp^10 * (-sym(3)/16) * (Lp^5/Gpp^5 - Lp^7/Gpp^7) ...
           * (1 - 16*H^2/Gpp^2 + 15*H^4/Gpp^4) * cos(2*gp);

F2p_star = simplify(F2p_star);
printf('F_{2p}* = '); disp(F2p_star);

%% ================================================================
%% Step 3: Compute dS1*/dg' = -F_{2p}* / (dF1*/dG'')
%% ================================================================

printf('\n--- Step 3: dS1*/dg'' = -F_{2p}* / (dF1*/dG'''') ---\n\n');

dS1star_dg = -F2p_star / dF1star_dGpp;
dS1star_dg = simplify(dS1star_dg);
printf('dS1*/dg'' = '); disp(dS1star_dg);

% Try to express in terms of gamma2 = mu^2*k2/L'^4 and theta = H/G'':
% Factor out G'' * gamma2:
dS1star_dg_over_Ggamma2 = simplify(dS1star_dg / (Gpp * gamma2));
printf('dS1*/dg'' / (G''''*gamma2) = '); disp(dS1star_dg_over_Ggamma2);

% Factor out cos(2g'):
dS1star_dg_over_cos2g = simplify(dS1star_dg / cos(2*gp));
printf('dS1*/dg'' / cos(2g'') = '); disp(dS1star_dg_over_cos2g);

%% ================================================================
%% Step 4: Compare with Brouwer (1959) line 466
%%
%% Brouwer claims:
%% dS1*/dg' = G''*gamma2 * [1/8*(L'^2/G''^2 - L'^4/G''^4)
%%             *(1-16*H^2/G''^2+15*H^4/G''^4)] * (1-5*H^2/G''^2)^{-1} * cos(2g')
%% ================================================================

printf('\n--- Step 4: Compare with Brouwer''s claimed result ---\n\n');

brouwer_claimed = Gpp * gamma2 * sym(1)/8 * (Lp^2/Gpp^2 - Lp^4/Gpp^4) ...
                  * (1 - 16*H^2/Gpp^2 + 15*H^4/Gpp^4) ...
                  * (1 - 5*H^2/Gpp^2)^(-1) * cos(2*gp);

diff_brouwer = simplify(dS1star_dg - brouwer_claimed);
printf('Our dS1*/dg'' - Brouwer''s claim = '); disp(diff_brouwer);
printf('(Should be 0 if Brouwer is correct)\n\n');

%% ================================================================
%% Step 5: Integrate to get S1*
%%
%% S1* = integral of (dS1*/dg') dg'
%% Since dS1*/dg' = [...]*cos(2g'), the integral is [...]*sin(2g')/2.
%% ================================================================

printf('--- Step 5: Integrate to get S1* ---\n\n');

S1star = int(dS1star_dg, gp);
S1star = simplify(S1star);
printf('S1* = '); disp(S1star);

% Compare with Brouwer (1959) line 478:
% S1* = G''*gamma2*(L'^2/G''^2 - L'^4/G''^4)
%       * [1/16*(1-11*H^2/G''^2) - 5/2*H^4/G''^4*(1-5*H^2/G''^2)^{-1}] * sin(2g'')

brouwer_S1star = Gpp * gamma2 * (Lp^2/Gpp^2 - Lp^4/Gpp^4) ...
                 * (sym(1)/16*(1-11*H^2/Gpp^2) - sym(5)/2*H^4/Gpp^4*(1-5*H^2/Gpp^2)^(-1)) ...
                 * sin(2*gp);

diff_S1star = simplify(S1star - brouwer_S1star);
printf('Our S1* - Brouwer''s S1* = '); disp(diff_S1star);
printf('(Should be 0 if consistent)\n\n');

%% ================================================================
%% Step 6: Compute partial derivatives of S1* needed for Eq. (14)
%%
%% dS1*/dL', dS1*/dG'', dS1*/dH
%% ================================================================

printf('--- Step 6: Partial derivatives of S1* ---\n\n');

dS1star_dLp = diff(S1star, Lp);
dS1star_dLp = simplify(dS1star_dLp);
printf('dS1*/dL'' = '); disp(dS1star_dLp);

dS1star_dGpp_full = diff(S1star, Gpp);
dS1star_dGpp_full = simplify(dS1star_dGpp_full);
printf('dS1*/dG'''' = '); disp(dS1star_dGpp_full);

dS1star_dH = diff(S1star, H);
dS1star_dH = simplify(dS1star_dH);
printf('dS1*/dH = '); disp(dS1star_dH);

printf('\n============================================================\n');
printf('DERIVATION COMPLETE\n');
printf('============================================================\n');
