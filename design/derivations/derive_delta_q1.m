%% derive_delta_q1.m
%% Derive delta_q_1 = D(dS1/dL + dS1*/dL) symbolically.
%%
%% This is the hardest of the six delta quantities because dS1/dL involves
%% implicit differentiation through the Kepler equation (f and r depend on L).
%%
%% Strategy: compute dS1/dL symbolically using SymPy, noting that:
%%   S1 = (mu^2 k2/G^3) * {A*(f-l+e*sin(f)) + B*[sin terms]}
%%   a = L^2/mu, e = sqrt(1-G^2/L^2), and f depends on L through e and l.
%%
%% Then apply D using the D-identities.
%%
%% For dS1*/dL: already computed in derive_S1star.m.

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING delta_q_1 = D(dS1/dL + dS1*/dL)\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L_s G_s H_s real positive;
syms l_s g_s f_s E_s e_s a_s r_s real;

theta_s = H_s/G_s;

%% ================================================================
%% Part 1: dS1/dL
%%
%% S1 = (mu^2 k2/G^3)*{A*(f-l+e*sin(f)) + B*[1/2*sin(2g+2f)+e/2*sin(2g+f)+e/6*sin(2g+3f)]}
%%
%% A = -1/2+3/2*theta^2 depends on G and H only, not L.
%% B = 3/2-3/2*theta^2 depends on G and H only, not L.
%% 1/G^3 does not depend on L.
%%
%% So dS1/dL = (mu^2 k2/G^3)*{A*d(f-l+e*sin(f))/dL + B*d[...]/dL}
%%
%% The L-dependence enters through e = sqrt(1-G^2/L^2) and through f(l,e).
%% de/dL = (1-e^2)/(L*e)  [computed in derive_Df.m]
%%
%% d(f)/dL at fixed l, G, H:
%%   df/dL = df/de * de/dL  (through the Kepler equation)
%%   This was computed implicitly in derive_Df.m.
%%
%% Rather than recompute all the chain rule derivatives manually,
%% let me use SymPy to differentiate S1 expressed as a function of
%% (L, G, H, f, e, l, g) with the constraint e = e(L,G) and f = f(l,e).
%%
%% Actually, the cleanest approach: express S1 in terms of L, G, H, l, g
%% using E as the intermediate variable (through Kepler's equation),
%% and let SymPy handle the implicit differentiation.
%%
%% But SymPy can't handle implicit differentiation through Kepler's equation
%% directly. Instead, use the result:
%%   d(anything)/dL|_{l,G,H,g} = (partial/partial e)*de/dL + (partial/partial f)*df/dL
%%   where df/dL was computed from the chain through E.
%%
%% From derive_Df.m:
%%   df/dL at fixed l = [df/dE * dE/dL] + [df/de|_E * de/dL]
%%
%% This is getting complex. Let me use a DIFFERENT approach.
%%
%% APPROACH: Use the identity from Eq. (9) of Brouwer (1959).
%% The generating function relations say:
%%   L = L' + dS1/dl  (the homological equation, Eq. 13)
%%   l' = l + dS1/dL' (Eq. 9 relation for coordinates)
%%
%% Wait — Brouwer (1959) Eq. (17) gives L in terms of L' directly:
%%   L = L'{1 + gamma2*[A*(a^3/r^3 - L'^3/G'^3) + B*(a^3/r^3)*cos(2g+2f)]}
%%
%% And Eq. (36) gives l' = l'' - dS1*/dL'.
%%
%% But what I need is dS1/dL', not these relations.
%%
%% Let me just compute dS1/dL by explicit differentiation of S1 w.r.t. L,
%% treating e and f as dependent on L through:
%%   e = sqrt(1-G^2/L^2)
%%   f = f(l, e) determined by Kepler's equation l = E - e*sin(E) and
%%       the true anomaly relation.
%%
%% Key partial derivatives:
%%   de/dL = (1-e^2)/(L*e)
%%   d(e*sin(f))/dL = de/dL * sin(f) + e*cos(f)*df/dL
%%   d(sin(2g+nf))/dL = n*cos(2g+nf)*df/dL
%%   d(f)/dL = df/de * de/dL  (at fixed l)
%%
%% df/de at fixed l: this is what we need. From Kepler + true anomaly:
%%   From the orbit equation r = a(1-e*cos(E)):
%%   From r*sin(f) = a*sqrt(1-e^2)*sin(E), r*cos(f) = a*(cos(E)-e):
%%   f = atan2(sqrt(1-e^2)*sin(E), cos(E)-e)
%%
%% df/de at fixed l involves dE/de at fixed l, which from Kepler:
%%   0 = dE/de - sin(E) - e*cos(E)*dE/de = (1-e*cos(E))*dE/de - sin(E)
%%   dE/de|_l = sin(E)/(1-e*cos(E)) = sin(E)*a/r
%%
%% Then df/de|_l = df/dE * dE/de|_l + df/de|_E
%%
%% This is extremely messy. Let me try a cleaner symbolic approach.
%% Express EVERYTHING in terms of E and e (eliminating f):
%%   sin(f) = sqrt(1-e^2)*sin(E)/(1-e*cos(E))
%%   cos(f) = (cos(E)-e)/(1-e*cos(E))
%%   r/a = 1-e*cos(E)
%%   l = E - e*sin(E)
%%
%% Then S1 as a function of (E, e, G, H, g):
%% ================================================================

printf('Expressing S1 in terms of E, e, G, H, g...\n\n');

% Define sin(f) and cos(f) in terms of E and e:
D_kep = 1 - e_s*cos(E_s);  % = r/a
sinf = sqrt(1-e_s^2)*sin(E_s)/D_kep;
cosf = (cos(E_s)-e_s)/D_kep;

A_coeff = -sym(1)/2 + sym(3)/2 * theta_s^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta_s^2;

% f - l: from Kepler l = E - e*sin(E), so f - l = f - E + e*sin(E).
% f = atan2(sinf_expr, cosf_expr), but we need f - l symbolically.
% Use: f - l = (f - E) + e*sin(E)
% And: tan((f-E)/2) = ... this is getting complicated.

% Actually, e*sin(f) is cleaner:
% e*sin(f) = e*sqrt(1-e^2)*sin(E)/D_kep
esinf = e_s*sqrt(1-e_s^2)*sin(E_s)/D_kep;

% S1 = (mu^2 k2/G^3) * {A*(f-l+e*sin(f)) + B*[...]}
% The term (f-l+e*sin(f)):
% f - l + e*sin(f) = f - (E-e*sin(E)) + e*sin(f) = (f-E) + e*sin(E) + e*sin(f)
%
% This still involves f explicitly. Let me instead use the KNOWN result
% that dS1/dl was already computed (= the homological equation RHS), and
% use the RELATION between dS1/dL and dS1/dl.
%
% From the generating function structure, there may be a simpler route.

% SIMPLEST APPROACH: Use SymPy to differentiate the CLOSED FORM of S1
% with respect to L, where everything is expressed through E and e.
% But f appears in S1, and f = f(E, e) is not elementary to substitute.

% ALTERNATIVE: Brouwer (1959) provides the result in Eq. (17):
%   L = L' + dS1/dl  => dS1/dl is known (homological equation)
% And Eq. (36) provides l' - l'' = -dS1*/dL'.
% But there's no explicit equation for dS1/dL in Brouwer (1959).
%
% However, BH61 Eq. (9) says L - L'' = +d(S1+S1*)/dl'' (not dL).
% And l - l'' = -d(S1+S1*)/dL''.
% So dS1/dL'' is related to (l - l'') via the first Eq. (9) relation.
%
% From Brouwer (1959) Eqs. (17)-(22), the coordinate transformations
% l' - l, g' - g, h' - h are given explicitly. These are precisely
% -dS1/dL', -dS1/dG', -dS1/dH'.
%
% So dS1/dL = -(l' - l). Brouwer (1959) Eq. (17) through the l' relation...
% Actually, Eq. (17) gives L, not l'. Let me find the l' equation.

printf('Looking for dS1/dL in Brouwer (1959)...\n');
printf('From the type-2 relation: l'' = l + dS1/dL''\n');
printf('=> dS1/dL = l'' - l = negative of the l-correction.\n\n');

% From Brouwer (1959), after Eq. (17), the l' equation should be given.
% I need to read it from the file.

printf('PAUSING: need to read Brouwer (1959) for the explicit l''-l expression.\n');
printf('This gives dS1/dL directly via the type-2 relation.\n');
