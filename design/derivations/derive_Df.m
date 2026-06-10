%% derive_Df.m
%% Derive Df = D(true anomaly) from first principles using Octave symbolic.
%%
%% Starting point: the orbit equation r = a(1-e^2)/(1+e*cos(f))
%% and the fact that Dr = 0 (Eq. 13b: r = |eta| has no xi dependence).
%%
%% The strategy: Dr = 0 gives a constraint relating Da, De, and Df.
%% We derive Da and De from the chain rule through the Delaunay variables,
%% then solve for Df.
%%
%% Inputs (all derived in earlier files, verified with Octave):
%%   - DL_j = -(p_j) from Eq. (12), where:
%%     p_1 = L*(2a/r - 1)   [file 07]
%%     p_2 = G              [file 07]
%%   - a = L^2/mu            [file 04]
%%   - e = sqrt(1-G^2/L^2)   [file 04]
%%   - r = a(1-e^2)/(1+e*cos(f))  [Kepler orbit equation, file 03]

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING Df FROM Dr = 0\n');
printf('============================================================\n\n');

syms a_s e_s f_s r_s real;

%% ================================================================
%% Step 1: Derive Da from the chain rule
%%
%% a = L^2/mu. Therefore:
%%   Da = d(L^2/mu)/dL * DL + d(L^2/mu)/dG * DG
%%      = (2L/mu)*DL + 0
%%      = (2L/mu)*(-(p_1))
%%      = (2L/mu)*(-L*(2a/r-1))
%%      = -2*L^2/mu * (2a/r - 1)
%%      = -2*a*(2a/r - 1)
%% ================================================================

printf('--- Step 1: Da ---\n\n');
printf('a = L^2/mu depends only on L.\n');
printf('Da = (2L/mu)*DL = (2L/mu)*(-(p_1)) = (2L/mu)*(-L*(2a/r-1))\n');
printf('   = -2*(L^2/mu)*(2a/r-1) = -2*a*(2a/r-1)\n\n');

Da = -2*a_s*(2*a_s/r_s - 1);
printf('Da = '); disp(Da);

%% ================================================================
%% Step 2: Derive De from the chain rule
%%
%% e = sqrt(1 - G^2/L^2). Therefore:
%%   de/dL = G^2/(L^3*e) = (1-e^2)/(L*e)     [using G^2 = L^2(1-e^2)]
%%   de/dG = -G/(L^2*e) = -sqrt(1-e^2)/(L*e)  [using G = L*sqrt(1-e^2)]
%%
%%   De = (de/dL)*DL + (de/dG)*DG
%%      = [(1-e^2)/(L*e)] * [-(p_1)] + [-sqrt(1-e^2)/(L*e)] * [-(p_2)]
%%      = [(1-e^2)/(L*e)] * [-L*(2a/r-1)] + [-sqrt(1-e^2)/(L*e)] * [-G]
%%      = -(1-e^2)*(2a/r-1)/e + G*sqrt(1-e^2)/(L*e)
%%      = -(1-e^2)*(2a/r-1)/e + (1-e^2)/e     [using G/L = sqrt(1-e^2)]
%%      = (1-e^2)/e * [-(2a/r-1) + 1]
%%      = (1-e^2)/e * (2 - 2a/r)
%%      = -2*(1-e^2)*(a/r - 1)/e
%% ================================================================

printf('--- Step 2: De ---\n\n');
printf('e = sqrt(1-G^2/L^2).\n');
printf('de/dL = (1-e^2)/(L*e),  de/dG = -sqrt(1-e^2)/(L*e)\n');
printf('De = de/dL * DL + de/dG * DG\n');
printf('   = [(1-e^2)/(Le)] * [-L(2a/r-1)] + [-sqrt(1-e^2)/(Le)] * [-G]\n');
printf('   = -(1-e^2)(2a/r-1)/e + (1-e^2)/e\n');
printf('   = (1-e^2)/e * [1 - (2a/r-1)]\n');
printf('   = -2*(1-e^2)*(a/r - 1)/e\n\n');

De = -2*(1-e_s^2)*(a_s/r_s - 1)/e_s;
printf('De = '); disp(De);

% Simplify using orbit equation: a/r = (1+e*cos(f))/(1-e^2)
% a/r - 1 = (1+e*cos(f) - 1+e^2)/(1-e^2) = (e^2+e*cos(f))/(1-e^2) = e*(e+cos(f))/(1-e^2)
% De = -2*(1-e^2)*e*(e+cos(f)) / [e*(1-e^2)] = -2*(e+cos(f))

De_simplified = subs(De, a_s/r_s, (1+e_s*cos(f_s))/(1-e_s^2));
De_simplified = simplify(De_simplified);
printf('De (using a/r = (1+e*cos(f))/(1-e^2)):\n');
printf('De = '); disp(De_simplified);

% Let SymPy verify: should be -2*(e+cos(f))
De_check = simplify(De_simplified + 2*(e_s + cos(f_s)));
printf('De + 2*(e+cos(f)) = '); disp(De_check);
printf('(Should be 0 if De = -2*(e+cos(f)))\n\n');

%% ================================================================
%% Step 3: Derive Df from Dr = 0
%%
%% r = a*(1-e^2)/(1+e*cos(f)). Since Dr = 0:
%%   0 = D[a*(1-e^2)/(1+e*cos(f))]
%%
%% Let P = a*(1-e^2), Q = 1+e*cos(f). Then r = P/Q.
%%   0 = D(P/Q) = [D(P)*Q - P*D(Q)] / Q^2
%%   => D(P)*Q = P*D(Q)
%%
%% D(P) = Da*(1-e^2) + a*(-2*e)*De
%% D(Q) = De*cos(f) + e*(-sin(f))*Df = De*cos(f) - e*sin(f)*Df
%%
%% Solve D(P)*(1+e*cos(f)) = a*(1-e^2)*[De*cos(f) - e*sin(f)*Df] for Df.
%% ================================================================

printf('--- Step 3: Df from Dr = 0 ---\n\n');

syms Df_var real;

DP = Da*(1-e_s^2) + a_s*(-2*e_s)*De;
DQ = De*cos(f_s) - e_s*sin(f_s)*Df_var;

% Equation: DP*Q = P*DQ where Q = 1+e*cos(f), P = a*(1-e^2)
Q_expr = 1 + e_s*cos(f_s);
P_expr = a_s*(1-e_s^2);

equation = DP*Q_expr - P_expr*DQ;
equation = expand(equation);

Df_solved = solve(equation, Df_var);
Df_solved = simplify(Df_solved);

printf('Df (raw symbolic solution) = '); disp(Df_solved);

% Substitute a/r = (1+e*cos(f))/(1-e^2) to eliminate r:
Df_sub = subs(Df_solved, r_s, a_s*(1-e_s^2)/(1+e_s*cos(f_s)));
Df_sub = simplify(Df_sub);
printf('Df (with orbit equation substituted) = '); disp(Df_sub);

% Check if this equals 2*sin(f)/e:
Df_check = simplify(Df_sub - 2*sin(f_s)/e_s);
printf('Df - 2*sin(f)/e = '); disp(Df_check);
printf('(Should be 0 if Df = 2*sin(f)/e)\n\n');

%% ================================================================
%% Step 4: Derive Dg from the chain rule
%%
%% g = l_2. From Eq. (12): Dl_j = +(q_j).
%% So Dg = Dl_2 = +(q_2).
%%
%% From file 07: q_2 involves the Delaunay-to-Cartesian partials.
%% But we can also derive Dg from the conservation of (f+g):
%%
%% (f+g) is the argument of latitude, which depends on the position
%% direction relative to the ascending node. Specifically:
%%   sin(f+g) = eta_3 / (r*sin(I))
%%   cos(f+g) = [eta_1*cos(h) + eta_2*sin(h)] / r
%%
%% Both of these depend only on position (eta) and geometry (I, h).
%% Since I is degree 0 (Eq. 13a) and h is degree 0 (l_3 = h depends
%% on the angular momentum direction, which is degree 0 as a ratio),
%% (f+g) depends only on degree-0 quantities.
%%
%% Therefore D(f+g) = 0, which gives Dg = -Df = -2*sin(f)/e.
%% ================================================================

printf('--- Step 4: Dg ---\n\n');
printf('(f+g) = argument of latitude = atan2(eta_3/(r*sin(I)), ...)\n');
printf('This depends on position and angular momentum direction only.\n');
printf('All components are degree 0 in xi.\n');
printf('Therefore D(f+g) = 0.\n');
printf('Dg = -Df = -2*sin(f)/e\n\n');

% Verify symbolically: f+g is a function of eta only (through the rotation
% matrix and r). Since D annihilates functions of eta alone, D(f+g) = 0.

Dg_result = -2*sin(f_s)/e_s;
printf('Dg = '); disp(Dg_result);

%% ================================================================
%% Step 5: Derive D(e*sin(f)) as a consistency check
%%
%% D(e*sin(f)) = De*sin(f) + e*cos(f)*Df
%%             = -2*(e+cos(f))*sin(f) + e*cos(f)*2*sin(f)/e
%%             = -2*e*sin(f) - 2*sin(f)*cos(f) + 2*sin(f)*cos(f)
%%             = -2*e*sin(f)
%% ================================================================

printf('--- Step 5: D(e*sin(f)) consistency check ---\n\n');

D_esinf = De_simplified*sin(f_s) + e_s*cos(f_s)*2*sin(f_s)/e_s;
D_esinf = simplify(D_esinf);
printf('D(e*sin(f)) = De*sin(f) + e*cos(f)*Df = '); disp(D_esinf);

D_esinf_check = simplify(D_esinf + 2*e_s*sin(f_s));
printf('D(e*sin(f)) + 2*e*sin(f) = '); disp(D_esinf_check);
printf('(Should be 0 if D(e*sin(f)) = -2*e*sin(f))\n\n');

%% ================================================================
%% SUMMARY
%% ================================================================

printf('============================================================\n');
printf('DERIVED D-OPERATOR IDENTITIES\n');
printf('============================================================\n\n');
printf('Da = -2*a*(2a/r - 1)         [from Da = (2L/mu)*DL, p_1 = L(2a/r-1)]\n');
printf('De = -2*(e + cos(f))          [from chain rule through L, G]\n');
printf('Df = +2*sin(f)/e              [from Dr = 0 and the orbit equation]\n');
printf('Dg = -2*sin(f)/e              [from D(f+g) = 0]\n');
printf('D(e*sin(f)) = -2*e*sin(f)     [product rule consistency check]\n');
printf('D(f+g) = 0                    [argument of latitude is degree 0]\n');
