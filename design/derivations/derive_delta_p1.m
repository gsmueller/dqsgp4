%% derive_delta_p1.m
%% Derive delta_p_1 = D(dS1/dl) symbolically using D-operator identities.
%%
%% Starting point: dS1/dl from the homological equation (Brouwer 1959, Eq. 13):
%%   dS1/dl = (mu^2 k2 / L^3) * (A*sigma1 + B*sigma2)
%%   where sigma1 = a^3/r^3 - L^3/G^3
%%         sigma2 = (a^3/r^3)*cos(2g+2f)
%%         A = -1/2 + 3/2*theta^2
%%         B = 3/2 - 3/2*theta^2
%%         theta = H/G = cos(I)
%%
%% D-identities used (all derived in derive_Df.m):
%%   Da = -2*a*(2a/r-1)
%%   D(f+g) = 0  =>  D(cos(2g+2f)) = 0
%%   D(theta) = 0  =>  D(A) = D(B) = 0
%%   D(r) = 0
%%   D(G^n) = -n*G^n  (Euler, degree 1)
%%   D(L^n): NOT -n*L^n. Must use chain rule: D(L^n) = n*L^{n-1}*DL = -n*L^n*(2a/r-1)
%%           since DL = -(p_1) = -L*(2a/r-1).

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING delta_p_1 = D(dS1/dl)\n');
printf('============================================================\n\n');

syms mu k2 a_s e_s f_s g_s r_s theta_s L_s G_s real;

% Shorthand for 2a/r - 1 (appears everywhere through p_1)
syms alpha_ar real;  % alpha_ar := 2*a/r - 1

%% ================================================================
%% Step 1: Write dS1/dl in terms of a, r, L, G, theta, f, g
%%
%% dS1/dl = (mu^2 k2 / L^3) * [A*(a^3/r^3 - L^3/G^3) + B*(a^3/r^3)*cos(2g+2f)]
%%
%% Expand into three terms:
%%   T1 = (mu^2 k2 / L^3) * A * a^3/r^3
%%   T2 = (mu^2 k2 / L^3) * A * (-L^3/G^3)  =  -mu^2 k2 A / G^3
%%   T3 = (mu^2 k2 / L^3) * B * (a^3/r^3) * cos(2g+2f)
%% ================================================================

printf('--- dS1/dl expanded into three terms ---\n\n');

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

T1 = mu^2 * k2 * A_coeff * a_s^3 / (L_s^3 * r_s^3);
T2 = -mu^2 * k2 * A_coeff / G_s^3;
T3 = mu^2 * k2 * B_coeff * a_s^3 * cos(2*g_s + 2*f_s) / (L_s^3 * r_s^3);

printf('T1 = (mu^2 k2 A / L^3) * a^3/r^3\n');
printf('T2 = -mu^2 k2 A / G^3\n');
printf('T3 = (mu^2 k2 B / L^3) * (a^3/r^3) * cos(2g+2f)\n\n');

%% ================================================================
%% Step 2: Apply D to each term using the product rule.
%%
%% For each term, identify the factors and their D-values:
%%   D(mu^2 k2) = 0  (constants)
%%   D(A) = D(-1/2 + 3/2*theta^2) = 0  (theta = H/G is degree-0)
%%   D(B) = 0  (same reason)
%%   D(cos(2g+2f)) = 0  (D(f+g) = 0, so D(2g+2f) = 0)
%%   D(1/r^3) = 0  (r = |eta| is degree-0)
%%   D(1/G^3) = 3/G^3  (Euler: D(G^{-3}) = -(-3)*G^{-3} = 3*G^{-3})
%%   D(a^3/L^3) = ?  (this is the non-trivial piece)
%%
%% Key: a = L^2/mu, so a^3/L^3 = L^6/(mu^3 * L^3) = L^3/mu^3.
%% Wait — that's just L^3/mu^3.
%% And D(L^3) = 3*L^2*DL = 3*L^2*(-L*(2a/r-1)) = -3*L^3*(2a/r-1).
%% So D(a^3/L^3) = D(L^3/mu^3) = -3*L^3*(2a/r-1)/mu^3 = -3*(a^3/L^3)*(2a/r-1).
%%
%% Alternatively, D(a^3) = 3*a^2*Da = 3*a^2*(-2a*(2a/r-1)) = -6*a^3*(2a/r-1).
%% And D(1/L^3) = D(L^{-3}) = -(-3)*L^{-3}*(2a/r-1)... NO.
%% D(L^n) = n*L^{n-1}*DL = n*L^{n-1}*(-L*(2a/r-1)) = -n*L^n*(2a/r-1).
%% So D(L^{-3}) = -(-3)*L^{-3}*(2a/r-1)... NO!
%% D(L^n) = n*L^{n-1}*DL. For n = -3:
%%   D(L^{-3}) = (-3)*L^{-4}*DL = (-3)*L^{-4}*(-L*(2a/r-1)) = 3*L^{-3}*(2a/r-1).
%%
%% So D(a^3/L^3) = D(a^3)*L^{-3} + a^3*D(L^{-3})
%%               = -6*a^3*(2a/r-1)*L^{-3} + a^3*3*L^{-3}*(2a/r-1)
%%               = (a^3/L^3)*(2a/r-1)*(-6+3)
%%               = -3*(a^3/L^3)*(2a/r-1)
%%
%% Consistent. Good. (Or just use a^3/L^3 = L^3/mu^3 directly.)
%% ================================================================

printf('--- Applying D to each term ---\n\n');

% D(T1):
% T1 = [mu^2 k2 A / r^3] * [a^3/L^3]
% D(T1) = [mu^2 k2 A / r^3] * D(a^3/L^3)     [constants and degree-0 pass through]
%       = [mu^2 k2 A / r^3] * (-3)*(a^3/L^3)*(2a/r-1)
%       = -3*(2a/r-1) * T1

printf('D(T1) = -3*(2a/r-1) * T1\n');
printf('  (because T1 = [degree-0 stuff] * a^3/L^3, and D(a^3/L^3) = D(L^3/mu^3) = -3*L^3*(2a/r-1)/mu^3)\n\n');

% D(T2):
% T2 = [-mu^2 k2 A] * G^{-3}
% D(T2) = [-mu^2 k2 A] * D(G^{-3})     [A is degree-0]
%       = [-mu^2 k2 A] * 3*G^{-3}       [Euler: D(G^{-3}) = 3*G^{-3}]
%       = 3 * T2

printf('D(T2) = 3 * T2\n');
printf('  (because D(G^{-3}) = 3*G^{-3}, Euler degree-1 homogeneity of G)\n\n');

% D(T3):
% T3 = [mu^2 k2 B * cos(2g+2f) / r^3] * [a^3/L^3]
% D(cos(2g+2f)) = 0 (since D(f+g) = 0)
% So T3 has the same structure as T1 with an extra degree-0 factor cos(2g+2f).
% D(T3) = -3*(2a/r-1) * T3

printf('D(T3) = -3*(2a/r-1) * T3\n');
printf('  (same as T1: cos(2g+2f) is degree-0 since D(f+g) = 0)\n\n');

%% ================================================================
%% Step 3: Assemble delta_p_1 = D(T1 + T2 + T3)
%%
%% delta_p_1 = -3*(2a/r-1)*T1 + 3*T2 - 3*(2a/r-1)*T3
%%           = -3*(2a/r-1)*(T1+T3) + 3*T2
%% ================================================================

printf('--- Assembling delta_p_1 ---\n\n');

% Expand T1, T2, T3 back:
% T1 = mu^2 k2 A a^3 / (L^3 r^3)
% T2 = -mu^2 k2 A / G^3
% T3 = mu^2 k2 B a^3 cos(2g+2f) / (L^3 r^3)

% delta_p_1 = -3*(2a/r-1)*[mu^2 k2 a^3/(L^3 r^3)]*[A + B*cos(2g+2f)]
%           + 3*[-mu^2 k2 A / G^3]

% Factor out 3*mu^2*k2/L^3:
% First group: -3*(2a/r-1)*a^3/(L^3*r^3) = -3*(mu^2 k2/L^3)*(2a/r-1)*a^3/r^3 / (mu^2 k2/L^3)
%   ... actually let me factor from the full expression.

% delta_p_1 = 3*mu^2*k2 * { -(2a/r-1)*A*a^3/(L^3*r^3) - A/G^3 -(2a/r-1)*B*a^3*cos(2g+2f)/(L^3*r^3) }

% Factor mu^2*k2/L^3:
% First term: -(2a/r-1)*A*a^3/(L^3*r^3) = (mu^2 k2/L^3) * [-(2a/r-1)*A*a^3/r^3] / (mu^2 k2)
% Hmm, let me just do this with SymPy.

delta_p1 = -3*(2*a_s/r_s - 1)*T1 + 3*T2 - 3*(2*a_s/r_s - 1)*T3;
delta_p1 = expand(delta_p1);

% Factor out 3*mu^2*k2/L^3:
prefactor = 3*mu^2*k2/L_s^3;
delta_p1_reduced = simplify(delta_p1 / prefactor);

printf('delta_p_1 / (3*mu^2*k2/L^3) = \n');
disp(delta_p1_reduced);

% Now substitute a = L^2/mu to express a^3/r^3 = L^6/(mu^3*r^3):
% a^3/(L^3*r^3) = L^6/(mu^3*L^3*r^3) = L^3/(mu^3*r^3)
% But we want to keep a/r as a unit. So a^3/r^3 = (a/r)^3.
% And a/r = L^2/(mu*r).

% Let me substitute differently. Use eta = G/L (= sqrt(1-e^2)):
% L^3/G^3 = 1/eta^3 = eta^{-3}.
% And a^3/r^3 stays as is.

% Expected BH61 result (Eq. 14, first line):
% delta_p_1 = 3*(mu^2 k2/L^3) * { A*[-eta^{-3} + (a^3/r^3)*(1-2a/r)]
%              + B*(a^3/r^3)*(1-2a/r)*cos(2g+2f) }

% My result:
% delta_p_1 = 3*mu^2*k2 * { -(2a/r-1)*A*a^3/(L^3*r^3) - A/G^3 - (2a/r-1)*B*a^3*cos(2g+2f)/(L^3*r^3) }
% = 3*(mu^2*k2/L^3) * { -(2a/r-1)*A*(a^3/r^3) - A*(L^3/G^3) - (2a/r-1)*B*(a^3/r^3)*cos(2g+2f) }
%
% The A terms: -(2a/r-1)*A*(a^3/r^3) - A*(L^3/G^3)
%            = A*[-(2a/r-1)*(a^3/r^3) - L^3/G^3]
%            = A*[(1-2a/r)*(a^3/r^3) - eta^{-3}]
%            = A*[-eta^{-3} + (a^3/r^3)*(1-2a/r)]
%
% The B term: -(2a/r-1)*B*(a^3/r^3)*cos(2g+2f) = B*(1-2a/r)*(a^3/r^3)*cos(2g+2f)
%
% So delta_p_1 = 3*(mu^2*k2/L^3) * { A*[-eta^{-3} + (a^3/r^3)*(1-2a/r)]
%                 + B*(a^3/r^3)*(1-2a/r)*cos(2g+2f) }

printf('\nRearranging into BH61 form:\n');
printf('delta_p_1 = 3*(mu^2*k2/L^3) * {\n');
printf('  A*[-eta^{-3} + (a^3/r^3)*(1-2a/r)]\n');
printf('  + B*(a^3/r^3)*(1-2a/r)*cos(2g+2f)\n');
printf('}\n\n');

% Verify symbolically: build the BH61 expression and check equality.
syms eta_s real positive;

BH61_dp1 = 3*mu^2*k2/L_s^3 * ( ...
  A_coeff * (-1/eta_s^3 + a_s^3/r_s^3*(1 - 2*a_s/r_s)) + ...
  B_coeff * a_s^3/r_s^3*(1 - 2*a_s/r_s)*cos(2*g_s + 2*f_s) ...
);

% Substitute eta = G/L in our result:
delta_p1_with_eta = subs(delta_p1, G_s, L_s*eta_s);

diff_check = simplify(expand(delta_p1_with_eta - BH61_dp1));
printf('Our result - BH61 Eq.(14) first line = ');
disp(diff_check);
printf('(Should be 0 if they match)\n\n');

%% ================================================================
%% Comparison with BH61 Eq. (14) using theta notation
%% ================================================================

printf('--- Translating to BH61 notation ---\n\n');
printf('BH61 writes A = -1/2+3/2*theta^2 and uses eta = G/L = sqrt(1-e^2).\n');
printf('Our result matches BH61 Eq. (14), line 1:\n');
printf('  delta_p_1 = 3*(mu^2*k2/L^3)*{(-1/2+3/2*theta^2)*[-eta^{-3}+(a/r)^3*(1-2a/r)]\n');
printf('              +(3/2-3/2*theta^2)*(a/r)^3*(1-2a/r)*cos(2g+2f)}\n');

printf('\n============================================================\n');
printf('delta_p_1 DERIVED. No S1* contribution (S1* does not depend on l_1).\n');
printf('============================================================\n');
