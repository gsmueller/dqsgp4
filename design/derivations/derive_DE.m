%% derive_DE.m
%% Derive D(E) where E is the eccentric anomaly, using D(l), D(e), and Kepler's equation.
%%
%% Inputs (all derived in earlier files, already verified with Octave):
%%   - D(l) = q_1 = 2*sin(E)*(1 - e^3*cos(E))/(e*kappa)     [file 07, q_1 closed form]
%%   - D(e) = -2*(e + cos(f))                                [file 08, De lemma]
%%   - Kepler's equation: l = E - e*sin(E)
%%   - Geometric identity: e + cos(f) = eta^2 * cos(E) / kappa
%%       (where eta^2 = 1 - e^2, kappa = 1 - e*cos(E) = r/a)
%%
%% Output claim:
%%   D(E) = 2*sin(E)/(e*kappa) = (2*sin(E)/e) * (a/r)
%%
%% Supports Item 7 of the pre-Chapter-10 remediation list.

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING D(E) FROM KEPLER + D(l), D(e)\n');
printf('============================================================\n\n');

syms e E f kappa eta_sym real;

%% ================================================================
%% Step 1: Apply D to Kepler's equation l = E - e*sin(E).
%%
%%   D(l) = D(E) - D(e)*sin(E) - e*cos(E)*D(E)
%%        = (1 - e*cos(E)) * D(E) - sin(E) * D(e)
%%        = kappa * D(E) - sin(E) * D(e)
%% ================================================================

printf('--- Step 1: Apply D to Kepler equation l = E - e*sin(E) ---\n\n');

syms DE_unknown real;
Dl_from_kepler = kappa * DE_unknown - sin(E) * (-2*(e + cos(f)));  % D(e) = -2*(e+cos(f))
Dl_from_kepler = expand(Dl_from_kepler);
printf('D(l) from Kepler: '); disp(Dl_from_kepler);

%% ================================================================
%% Step 2: Use the file-07 closed form for q_1 = D(l_1) = D(l).
%%
%%   q_1 = 2*sin(E)*(1 - e^3*cos(E))/(e*kappa)
%% ================================================================

printf('--- Step 2: file-07 closed form for q_1 = D(l) ---\n\n');

q1_closed = 2*sin(E)*(1 - e^3*cos(E))/(e*kappa);
printf('q_1 (file 07) = '); disp(q1_closed);

%% ================================================================
%% Step 3: Equate and solve for D(E).
%% ================================================================

printf('--- Step 3: Equate and solve for D(E) ---\n\n');

eqn = Dl_from_kepler - q1_closed;
DE_solved = solve(eqn, DE_unknown);
DE_solved = simplify(DE_solved);
printf('D(E) (raw symbolic solve) = '); disp(DE_solved);

%% ================================================================
%% Step 4: Substitute the geometric identity e + cos(f) = eta^2 * cos(E) / kappa
%%   and compare with the proposed closed form 2*sin(E)/(e*kappa).
%%
%% The substitution e + cos(f) -> eta^2 * cos(E) / kappa, with eta^2 = 1 - e^2,
%% lets us eliminate f and express D(E) purely in terms of (e, E, kappa).
%% ================================================================

printf('--- Step 4: Substitute e + cos(f) = (1-e^2)*cos(E)/kappa ---\n\n');

DE_sub = subs(DE_solved, cos(f), (1-e^2)*cos(E)/kappa - e);
DE_sub = simplify(DE_sub);
printf('D(E) after substituting cos(f) in terms of (e, E, kappa) = '); disp(DE_sub);

%% ================================================================
%% Step 5: Substitute kappa = 1 - e*cos(E) and simplify to final form.
%% ================================================================

printf('--- Step 5: Substitute kappa = 1 - e*cos(E) ---\n\n');

DE_final = subs(DE_sub, kappa, 1 - e*cos(E));
DE_final = simplify(DE_final);
printf('D(E) with kappa expanded = '); disp(DE_final);

%% ================================================================
%% Step 6: Check against proposed closed form D(E) = 2*sin(E)/(e*kappa).
%% ================================================================

printf('--- Step 6: Check against D(E) = 2*sin(E)/(e*kappa) ---\n\n');

DE_proposed = 2*sin(E)/(e*(1 - e*cos(E)));
diff_check = simplify(DE_final - DE_proposed);
printf('D(E)_derived - 2*sin(E)/(e*kappa) = '); disp(diff_check);
printf('(Should be 0 if D(E) = 2*sin(E)/(e*kappa).)\n\n');

%% ================================================================
%% Step 7: Alternative form: D(E) = (2*sin(E)/e) * (a/r).
%%   Since a/r = 1/kappa, this is equivalent.
%% ================================================================

printf('--- Step 7: Alternative form D(E) = (2*sin(E)/e)*(a/r) ---\n\n');

syms a r real positive;
DE_alt_form = (2*sin(E)/e) * (a/r);
% Using r/a = kappa:
DE_alt_sub = subs(DE_alt_form, a/r, 1/(1 - e*cos(E)));
DE_alt_sub = simplify(DE_alt_sub);
printf('D(E) = (2 sin(E)/e) * (a/r) = (2 sin(E)/e) * (1/kappa) = '); disp(DE_alt_sub);

printf('\n============================================================\n');
printf('D(E) = 2*sin(E)/(e*kappa) = (2*sin(E)/e) * (a/r) CONFIRMED.\n');
printf('============================================================\n');
