%% derive_anomaly_partials.m
%% Derive df/de|_l, df/dL|_l,G, df/dG|_l,L symbolically.
%%
%% Strategy: work with (cos E, sin E) as SYMBOLS cE, sE. Use the Kepler
%% equation l = E - e sin E (fixed l) to get dE/de implicitly, then the
%% orbit-geometry identities
%%    cos f = (cos E - e)/(1-e cos E)   [Eq. G1]
%%    sin f = eta sin E/(1-e cos E)     [Eq. G2]
%% to express df/de|_l. All partials stored in (e, f, eta) form suitable
%% for substituting into S1.
%%
%% Primitives (each proven in-script):
%%   [P1] de/dL|_G = G^2/(L^3 e)          (from e = sqrt(1-G^2/L^2))
%%   [P2] de/dG|_L = -G/(L^2 e)
%%   [P3] dE/de|_l = sin E/(1-e cos E)    (implicit diff of Kepler)
%%   [P4] df/de|_l expressed via Eqs G1 or G2 with Chain rule from [P3]

pkg load symbolic;

printf('============================================================\n');
printf('ANOMALY PARTIALS: df/de|_l, df/dL|_l, df/dG|_l (symbolic)\n');
printf('============================================================\n\n');

syms L G e eta real positive;
syms f real;
syms sE cE real;        % sin(E), cos(E) as symbols
% Constraint (used for substitution, not as a sym equation): sE^2 + cE^2 = 1

% ----------------------------------------------------------------
% [P1], [P2]: de/dL and de/dG
% ----------------------------------------------------------------
printf('[P1] e = sqrt(1-G^2/L^2)\n');
e_of_LG = sqrt(1 - G^2/L^2);
de_dL_raw = diff(e_of_LG, L);
de_dG_raw = diff(e_of_LG, G);
de_dL = simplify(subs(de_dL_raw, sqrt(1-G^2/L^2), e));
de_dG = simplify(subs(de_dG_raw, sqrt(1-G^2/L^2), e));
printf('  de/dL|_G = '); disp(de_dL);
printf('  de/dG|_L = '); disp(de_dG);

printf('\n[P2] Cross-checks:\n');
printf('  L^3 e (de/dL) - G^2 = '); disp(simplify(L^3*e*de_dL - G^2));
printf('  L^2 e (de/dG) + G   = '); disp(simplify(L^2*e*de_dG + G));

% ----------------------------------------------------------------
% [P3]: dE/de|_l from implicit differentiation of Kepler's equation.
%
%   l = E - e sin E, fix l.
%   0 = dE/de - sin E - e cos E dE/de
%   => dE/de|_l = sin E / (1 - e cos E).
%
% Store as a value assigned to cE, sE symbols:
% ----------------------------------------------------------------
printf('\n[P3] dE/de|_l = sin E / (1-e cos E)\n');
dE_de_l = sE / (1 - e*cE);

% ----------------------------------------------------------------
% [P4]: df/de|_l via the identity cos f = (cos E - e)/(1-e cos E).
%
% Differentiate both sides at fixed l:
%   -sin(f) df/de = d/de[(cE - e)/(1-e cE)]|_l
% RHS: apply quotient rule with cE = cE(e) satisfying dcE/de = -sE dE/de.
% ----------------------------------------------------------------
printf('\n[P4] df/de|_l via d/de[cos f] = -sin f * df/de|_l\n');

% Expression for cos f in terms of symbols:
cos_f_expr = (cE - e)/(1 - e*cE);

% Total derivative of cos_f_expr at fixed l:
% Partial wrt e (holding cE fixed) + partial wrt cE * dcE/de|_l.
d_cosf_de_partial = diff(cos_f_expr, e);
d_cosf_de_via_cE  = diff(cos_f_expr, cE);
dcE_de_l = -sE * dE_de_l;  % chain through E
d_cosf_de = d_cosf_de_partial + d_cosf_de_via_cE * dcE_de_l;
d_cosf_de = simplify(d_cosf_de);
printf('  d(cos f)/de|_l = '); disp(d_cosf_de);

% sin(f) from [G2]: sin(f) = eta*sE/(1-e cE). So:
sin_f_expr = eta*sE/(1 - e*cE);

df_de_l_sym = simplify(-d_cosf_de / sin_f_expr);
printf('  df/de|_l (in cE, sE, e, eta) = '); disp(df_de_l_sym);

% Substitute sE^2 -> 1-cE^2 to simplify:
df_de_l_sym2 = simplify(subs(df_de_l_sym, sE^2, 1-cE^2));
printf('  df/de|_l after sE^2->1-cE^2 = '); disp(df_de_l_sym2);

% ----------------------------------------------------------------
% Convert to (e,f) form:
%   From G1: cos E = (cos f + e)/(1+e cos f)
%   From G2: sin E = eta sin f/(1+e cos f)      [derivable from G2 inverse]
%   1 - e cos E = (1-e^2)/(1+e cos f) = eta^2/(1+e cos f)
% ----------------------------------------------------------------
printf('\n[P4] Converting to (e,f):\n');
cE_in_f = (cos(f) + e)/(1 + e*cos(f));
sE_in_f = eta*sin(f)/(1 + e*cos(f));

df_de_l_ef = subs(df_de_l_sym2, cE, cE_in_f);
df_de_l_ef = subs(df_de_l_ef, sE, sE_in_f);
df_de_l_ef = simplify(df_de_l_ef);
printf('  df/de|_l in (e,f,eta) = '); disp(df_de_l_ef);

% Expand with eta^2 = 1-e^2 for comparison to classical form:
df_de_l_expanded = simplify(subs(df_de_l_ef, eta^2, 1-e^2));
% In case it's in 1/eta form, try subs eta -> sqrt(1-e^2):
df_de_l_via_e = simplify(subs(df_de_l_ef, eta, sqrt(1-e^2)));
printf('  df/de|_l with eta->sqrt(1-e^2) = '); disp(df_de_l_via_e);

% Compare with classical: df/de|_l = sin(f)(2 + e cos f)/(1-e^2)
classical = sin(f)*(2 + e*cos(f))/(1-e^2);
diff1 = simplify(df_de_l_via_e - classical);
printf('  [check] df/de|_l - sin(f)(2+e cos f)/(1-e^2) = '); disp(diff1);

% ----------------------------------------------------------------
% [P5]: df/dL|_l,G and df/dG|_l,L by chain rule.
%   df/dL|_l = (df/de|_l) * (de/dL|_G)
%   df/dG|_l = (df/de|_l) * (de/dG|_L)
% Using results [P1], [P4].
% ----------------------------------------------------------------
printf('\n[P5] Chain: df/dL|_l,G = (df/de|_l)(de/dL), df/dG|_l = (df/de|_l)(de/dG)\n');

df_dL = simplify(df_de_l_ef * de_dL);
df_dG = simplify(df_de_l_ef * de_dG);
printf('  df/dL|_l,G = '); disp(df_dL);
printf('  df/dG|_l,L = '); disp(df_dG);

% Save results to a .mat file for downstream scripts? Octave symbolic
% doesn't play nicely with save. Just print the final values clearly.
printf('\n============================================================\n');
printf('FINAL RESULTS (for use in derive_dS1_dLGH.m)\n');
printf('============================================================\n');
printf('de/dL|_G    = G^2/(L^3*e)\n');
printf('de/dG|_L    = -G/(L^2*e)\n');
printf('df/de|_l    = sin(f)*(2+e cos(f))/(1-e^2)    [classical form, verified]\n');
printf('df/dL|_l,G  = G^2/(L^3*e) * sin(f)*(2+e cos(f))/(1-e^2)\n');
printf('            = sin(f)*(2+e cos(f))*G^2 / (L^3 e eta^2)\n');
printf('              with eta^2 = 1-e^2, G = L*eta, so G^2/(L^3 eta^2) = 1/L\n');
printf('            = sin(f)*(2+e cos(f))/(L e)     [clean form]\n');
printf('df/dG|_l,L  = -G/(L^2*e) * sin(f)*(2+e cos(f))/(1-e^2)\n');
printf('            = -sin(f)*(2+e cos(f))/(L eta e)  using G=L eta\n');

% Verify the clean forms:
printf('\n[check] df/dL with G=L*eta, eta^2=1-e^2:\n');
df_dL_clean = sin(f)*(2+e*cos(f))/(L*e);
df_dL_from_chain = subs(df_dL, G, L*eta);
df_dL_from_chain = subs(df_dL_from_chain, eta^2, 1-e^2);
df_dL_from_chain = subs(df_dL_from_chain, eta, sqrt(1-e^2));
df_dL_from_chain = simplify(df_dL_from_chain);
printf('  df/dL (chain result, G=L*eta): '); disp(df_dL_from_chain);
printf('  df/dL (clean form)           : '); disp(df_dL_clean);
printf('  diff = '); disp(simplify(df_dL_from_chain - df_dL_clean));

printf('\n[check] df/dG with G=L*eta:\n');
df_dG_clean = -sin(f)*(2+e*cos(f))/(L*eta*e);
df_dG_from_chain = subs(df_dG, G, L*eta);
df_dG_from_chain = subs(df_dG_from_chain, eta^2, 1-e^2);
df_dG_from_chain = subs(df_dG_from_chain, eta, sqrt(1-e^2));
df_dG_from_chain_alt = subs(df_dG_clean, eta, sqrt(1-e^2));
df_dG_from_chain = simplify(df_dG_from_chain);
df_dG_from_chain_alt = simplify(df_dG_from_chain_alt);
printf('  df/dG (chain result):     '); disp(df_dG_from_chain);
printf('  df/dG (clean * sqrt):     '); disp(df_dG_from_chain_alt);
printf('  diff = '); disp(simplify(df_dG_from_chain - df_dG_from_chain_alt));

printf('\n============================================================\n');
