%% derive_eq14_symbolic.m
%% Symbolic derivation of BH61 Eq. (14): delta_p_j and delta_q_j
%%
%% Method:
%%   delta_p_j = D(dS1/dl_j) + D(dS1*/dl_j)     [from Eq. 11a]
%%   delta_q_j = D(dS1/dL_j) + D(dS1*/dL_j)     [from Eq. 11b]
%%
%% where D = -sum_k xi_k d/d(xi_k) is the Euler velocity-homogeneity operator.
%%
%% Theorems used:
%%   - Eq. (11): delta_p_j and delta_q_j in D-operator form
%%   - Eq. (12): DL_j = -(p_j), Dl_j = +(q_j)
%%   - Eq. (13): D annihilates functions degree-0 in xi
%%                (r, f, g, e, a, theta, eta, I are all annihilated)
%%   - Lemma: D acts as -d on degree-d homogeneous functions in xi
%%   - Product rule: D(fg) = fDg + gDf
%%
%% Key insight: S1 and S1* are expressed in orbital elements (a, e, f, g, L, G, H).
%% When viewed as functions of (xi, eta), the orbital elements have definite
%% homogeneity degrees in xi:
%%   degree 0: a, e, f, g, h, I, r, theta, eta  (geometry, no speed)
%%   degree 1: L, G, H (= L1, L2, L3)           (Lemma from Eq. 12 section)
%%   degree 2: mu^2/L^3 * (stuff)                (depends on context)
%%
%% Therefore D passes through all degree-0 factors and acts only on the
%% L, G, H factors via DL = -L, DG = -G, DH = -H (Euler, degree 1).
%%
%% Source for S1: Brouwer (1959) Eq. (15)
%% Source for S1*: Brouwer (1959) Eq. (34) / line 478

pkg load symbolic;

printf('============================================================\n');
printf('SYMBOLIC DERIVATION: BH61 Eq. (14)\n');
printf('============================================================\n\n');

%% ================================================================
%% Define symbolic variables
%% ================================================================

syms mu k2 real positive;
syms L G H real positive;      % Delaunay momenta (= L1, L2, L3)
syms l g h real;               % Delaunay angles (= l1, l2, l3)
syms f E_var real;             % true anomaly, eccentric anomaly
syms a e eta_orb theta real;   % orbital elements (functions of L, G, H)
syms r real positive;          % radial distance

% Relationships (these are degree-0 in xi, so D annihilates them):
%   a = L^2/mu
%   e = sqrt(1 - G^2/L^2)
%   eta_orb = G/L = sqrt(1-e^2)
%   theta = H/G = cos(I)
%   r = a(1 - e*cos(E))
%   f = true anomaly (implicit function of l, e)

% For symbolic manipulation, keep a, e, eta, theta, r, f as independent
% symbols and substitute relationships only when needed.

printf('Variables defined.\n\n');

%% ================================================================
%% S1 from Brouwer (1959) Eq. (15)
%% ================================================================
%%
%% S1 = (mu^2 k2 / G^3) * { A*(f - l + e*sin(f))
%%       + B*[1/2*sin(2g+2f) + e/2*sin(2g+f) + e/6*sin(2g+3f)] }
%%
%% where A = -1/2 + 3/2 * H^2/G^2 = -1/2 + 3/2 * theta^2
%%       B = +3/2 - 3/2 * H^2/G^2 = +3/2 - 3/2 * theta^2

printf('--- S1: Brouwer (1959) Eq. (15) ---\n\n');

A_coeff = sym(-1)/2 + sym(3)/2 * theta^2;
B_coeff = sym(3)/2 - sym(3)/2 * theta^2;

S1 = (mu^2 * k2 / G^3) * ( ...
  A_coeff * (f - l + e*sin(f)) + ...
  B_coeff * (sym(1)/2*sin(2*g + 2*f) + e/2*sin(2*g + f) + e/6*sin(2*g + 3*f)) ...
);

printf('S1 defined.\n');

%% ================================================================
%% STEP 1: Compute dS1/dl_1 = dS1/dl
%%
%% l appears explicitly in S1 only in the (f - l) term.
%% f also depends on l (through Kepler's equation l = E - e*sin(E),
%% then f = f(E, e)). But BH61's notation treats S1 as a function
%% of (L, G, H, l, g) where f = f(L, G, l).
%%
%% However, in the D-operator framework, we need dS1/dl evaluated
%% at the mean elements. The key identity from Brouwer (1959) Eq. (13):
%%
%%   dS1/dl = (mu^2 k2 / L^3) * (A*sigma1 + B*sigma2)
%%
%% where sigma1 = a^3/r^3 - L^3/G^3 and sigma2 = (a^3/r^3)*cos(2g+2f).
%%
%% This is the homological equation itself — dS1/dl was determined BY
%% the condition that it equals F1_tilde / (dF0/dL).
%% ================================================================

printf('\n--- dS1/dl (from homological equation, Brouwer Eq. 13) ---\n\n');

sigma1 = a^3/r^3 - L^3/G^3;
sigma2 = (a^3/r^3) * cos(2*g + 2*f);

dS1_dl = (mu^2 * k2 / L^3) * (A_coeff * sigma1 + B_coeff * sigma2);

printf('dS1/dl = (mu^2 k2 / L^3) * (A*sigma1 + B*sigma2)\n');
printf('  where sigma1 = a^3/r^3 - L^3/G^3\n');
printf('        sigma2 = (a^3/r^3)*cos(2g+2f)\n\n');

%% ================================================================
%% STEP 2: Apply D to dS1/dl
%%
%% D(dS1/dl) = D[ (mu^2 k2 / L^3) * (A*sigma1 + B*sigma2) ]
%%
%% Homogeneity analysis:
%%   mu^2 k2: constants, degree 0 → D passes through
%%   1/L^3: L is degree-1 in xi, so L^3 is degree-3, so 1/L^3 is degree -3
%%   A = -1/2 + 3/2*theta^2: theta = H/G is degree 0 → D(A) = 0
%%   B = 3/2 - 3/2*theta^2: same, D(B) = 0
%%   sigma1 = a^3/r^3 - L^3/G^3:
%%     a^3/r^3: a = L^2/mu (degree 2), r is degree 0 → a^3/r^3 is degree 6
%%              Wait — a = L^2/mu, so a is degree 2 in xi? No.
%%              L is degree 1, so L^2 is degree 2, a = L^2/mu is degree 2.
%%              r is degree 0. So a^3/r^3 is degree 6.
%%     L^3/G^3: degree 3 / degree 3 = degree 0. Wait, that's wrong.
%%              L^3 is degree 3, G^3 is degree 3, so L^3/G^3 is degree 0.
%%              But a^3/r^3 = L^6/(mu^3 r^3) is degree 6.
%%              These have DIFFERENT degrees — sigma1 is NOT homogeneous!
%%
%% This means D does NOT simply multiply sigma1 by a constant.
%% We must use the product rule and chain rule carefully.
%%
%% Better approach: factor out powers of L and G, then use
%%   D(L^n) = -n * L^n  (Euler, degree n)
%%   D(G^n) = -n * G^n  (Euler, degree n)
%%   D(f(...degree-0...)) = 0
%% ================================================================

printf('--- Applying D to dS1/dl ---\n\n');
printf('Strategy: rewrite everything in terms of L, G, H and degree-0 quantities,\n');
printf('then apply D(L^n) = -n*L^n, D(G^n) = -n*G^n, D(degree-0) = 0.\n\n');

% Rewrite dS1/dl with explicit L, G dependence:
%
% dS1/dl = mu^2 k2 / L^3 * [ A * (a^3/r^3 - L^3/G^3) + B * a^3/r^3 * cos(2g+2f) ]
%
% a = L^2/mu, so a^3 = L^6/mu^3, and a^3/r^3 = L^6/(mu^3 r^3)
%
% Term 1: mu^2 k2 / L^3 * A * L^6/(mu^3 r^3)
%        = k2 * A * L^3 / (mu * r^3)
%        homogeneity: L^3 is degree 3, everything else degree 0
%        D(term1) = D[k2 * A / (mu * r^3) * L^3]
%                 = k2 * A / (mu * r^3) * D(L^3)
%                 = k2 * A / (mu * r^3) * (-3 * L^3)
%                 = -3 * term1
%
% Term 2: mu^2 k2 / L^3 * A * (-L^3/G^3)
%        = -mu^2 k2 * A / G^3
%        homogeneity: 1/G^3 is degree -3
%        D(term2) = D[-mu^2 k2 * A * G^{-3}]
%                 = -mu^2 k2 * A * D(G^{-3})
%                 = -mu^2 k2 * A * (-(-3) * G^{-3})
%                 = -mu^2 k2 * A * (3 * G^{-3})
%                 = 3 * term2
%        Wait: D(G^n) = -n * G^n. So D(G^{-3}) = -(-3)*G^{-3} = 3*G^{-3}.
%        So D(term2) = -mu^2 k2 * A * 3 * G^{-3} = 3 * (-mu^2 k2 A / G^3) = 3*term2.
%
% Term 3: mu^2 k2 / L^3 * B * a^3/r^3 * cos(2g+2f)
%        = k2 * B * L^3 / (mu * r^3) * cos(2g+2f)
%        Same structure as term1 but with cos factor (degree 0)
%        D(term3) = -3 * term3
%
% So: D(dS1/dl) = -3*term1 + 3*term2 - 3*term3
%               = -3*(term1 + term3) + 3*term2
%               = -3*(mu^2 k2/L^3) * [A * a^3/r^3 + B * a^3/r^3 * cos(2g+2f)]
%                 + 3*(-mu^2 k2 A / G^3)  ...
%
% Actually let me do this more carefully. Factor the original:
%
% dS1/dl = (mu^2 k2/L^3) * A * a^3/r^3
%        + (mu^2 k2/L^3) * A * (-L^3/G^3)
%        + (mu^2 k2/L^3) * B * a^3/r^3 * cos(2g+2f)
%
% = k2*A*L^3/(mu*r^3) - mu^2*k2*A/G^3 + k2*B*L^3/(mu*r^3)*cos(2g+2f)

% Let me define these symbolically and let SymPy handle the D operator.

% Define: for any expression of the form c * L^m * G^n * H^p * (degree-0 stuff),
% D acts as: D(expr) = -(m+n+p) * expr
% because D(L^m G^n H^p) = -(m+n+p) * L^m G^n H^p (sum of degrees).

% Actually, D(L^m * G^n) = D(L^m)*G^n + L^m*D(G^n) = -m*L^m*G^n - n*L^m*G^n
%                        = -(m+n)*L^m*G^n. Yes.

% So I need to express each term in the form (degree-0 factor) * L^m * G^n * H^p.

% Rewrite with a = L^2/mu:
% a^3/r^3 = L^6/(mu^3 * r^3)

% Term by term:
% T1 = (mu^2 k2 / L^3) * A * L^6/(mu^3 r^3)
%    = k2 * A * L^3 / (mu * r^3)
%    = [k2 * A / (mu * r^3)] * L^3 * G^0 * H^0
%    degree = 3, so D(T1) = -3 * T1

% T2 = (mu^2 k2 / L^3) * A * (-L^3/G^3)
%    = -mu^2 * k2 * A / G^3
%    = [-mu^2 * k2 * A] * L^0 * G^{-3} * H^0
%    degree = 0 + (-3) + 0 = -3, so D(T2) = -(-3) * T2 = 3 * T2

% T3 = (mu^2 k2 / L^3) * B * L^6/(mu^3 r^3) * cos(2g+2f)
%    = k2 * B * L^3 / (mu * r^3) * cos(2g+2f)
%    = [k2 * B * cos(2g+2f) / (mu * r^3)] * L^3
%    degree = 3, so D(T3) = -3 * T3

% Therefore:
% D(dS1/dl) = -3*T1 + 3*T2 - 3*T3

% Let's compute this symbolically:

T1 = k2 * A_coeff * L^3 / (mu * r^3);
T2 = -mu^2 * k2 * A_coeff / G^3;
T3 = k2 * B_coeff * L^3 / (mu * r^3) * cos(2*g + 2*f);

D_dS1_dl = -3*T1 + 3*T2 - 3*T3;

% Simplify
D_dS1_dl = simplify(D_dS1_dl);

printf('D(dS1/dl) = \n');
disp(D_dS1_dl);

% Substitute a = L^2/mu to get a^3/r^3 = L^6/(mu^3*r^3), i.e. L^3/(mu*r^3) = a^3/(L^3*r^3/L^3)...
% Actually, let's substitute back: L^3/(mu*r^3) = (L^2/mu)^3 * mu^2 / (L^3 * r^3) ... messy.
% Better: use a = L^2/mu, so L^3 = mu^(3/2) * a^(3/2). Then L^3/(mu*r^3) = a^3*mu^2/(mu*L^3*r^3)...
% Let me just substitute L^3 = mu*a*L = mu*sqrt(mu*a)*a^(1/2)... this is getting circular.

% Instead, express in BH61 notation directly:
% mu^2 k2 / L^3 = gamma2 / L (since gamma2 = mu^2 k2 / L^4)
% But BH61 Eq. (14) uses mu^2 k2 / L^3 as prefactor for delta_p_1.

% Let's write the result in the form BH61 uses:
% D(dS1/dl) = -3*(mu^2 k2/L^3)*[A*(a^3/r^3) + B*(a^3/r^3)*cos(2g+2f)]
%             + 3*(-mu^2 k2 A / G^3)
%           = -3*(mu^2 k2/L^3)*A*a^3/r^3 - 3*(mu^2 k2/L^3)*B*(a^3/r^3)*cos(2g+2f)
%             - 3*mu^2 k2 A / G^3

% Recombine the first and third terms:
% = -3*(mu^2 k2/L^3)*[A*a^3/r^3 + B*(a^3/r^3)*cos(2g+2f)] - 3*mu^2*k2*A/G^3
% = -3*(mu^2 k2/L^3)*[A*(a^3/r^3 - L^3/G^3) + B*(a^3/r^3)*cos(2g+2f)]
%   - 3*(mu^2 k2/L^3)*A*L^3/G^3 - 3*mu^2*k2*A/G^3
%
% Hmm, that's not simplifying well. Let me just factor differently:
% D(dS1/dl) = -3*T1 + 3*T2 - 3*T3
%           = -3*(T1 + T3) + 3*T2
%           = -3*(T1 + T3 - T2)
%           Nope: -3*(T1+T3) + 3*T2 ≠ -3*(T1+T3-T2)

% Let me just substitute a^3/r^3 with the symbol and print:

% Actually the key observation: delta_p_1 = D(dS1/dl_1) and BH61 Eq. (14) has
% delta_p_1 = 3*(mu^2*k2/L^3)*{...}
%
% The factor of 3 comes from D acting on the L^3 and G^{-3} factors.
% Let me verify by expanding BH61's claimed result.

printf('\n--- Comparing with BH61 Eq. (14), delta_p_1 (S1 contribution only) ---\n\n');

% BH61 says (first line of Eq. 14):
% delta_p_1 = 3*(mu^2 k2/L^3) * { A*[-eta^{-3} + (a^3/r^3)*(1 - 2*a/r)]
%             + B*(a^3/r^3)*(1 - 2*a/r)*cos(2g+2f) }
%
% where A = -1/2 + 3/2*theta^2, B = 3/2 - 3/2*theta^2, eta = G/L.
%
% Note eta^{-3} = L^3/G^3. So the A term is:
% A*[-L^3/G^3 + a^3/r^3 - 2*a^4/r^4]
%
% But our D(dS1/dl) = -3*T1 + 3*T2 - 3*T3
% = -3*k2*A*L^3/(mu*r^3) + 3*mu^2*k2*A/G^3 - 3*k2*B*L^3*cos(2g+2f)/(mu*r^3)
%
% Using L^3/(mu*r^3) = a^3/(a*mu*r^3/L^2)... Let me use a = L^2/mu:
% L^3/(mu*r^3) = L*a/r^3... no. L^2 = mu*a, so L = sqrt(mu*a).
% L^3 = (mu*a)^{3/2} = mu^{3/2}*a^{3/2}.
% L^3/(mu*r^3) = mu^{1/2}*a^{3/2}/r^3.
% And mu^2*k2/L^3 = mu^2*k2/(mu^{3/2}*a^{3/2}) = mu^{1/2}*k2/a^{3/2}.
%
% Hmm, this is not leading to a clean factorization.
% The issue is that BH61's Eq. (14) has a (1 - 2*a/r) factor that I don't
% yet see in my result. That factor must come from D acting on a^3/r^3.
%
% WAIT. I made an error. a = L^2/mu is NOT degree 0 in xi!
% a depends on L, and L is degree 1 in xi. So a = L^2/mu is degree 2.
% Therefore a^3/r^3 = L^6/(mu^3*r^3) is degree 6, NOT degree 0.
%
% But I already accounted for this: T1 = k2*A*L^3/(mu*r^3) where
% L^3/(mu*r^3) = a^3*mu^2/(L^3*r^3*mu)... I'm going in circles.
%
% The problem is that a = L^2/mu is degree 2 in xi, so a/r is degree 2,
% and f (the true anomaly) depends on both l and e = sqrt(1-G^2/L^2),
% where e depends on G/L which is degree 0. But f also depends on l and E
% through Kepler's equation.
%
% Actually: f is determined by l and e. Since l is a Delaunay angle and
% e = sqrt(1-G^2/L^2) is degree 0 in xi (ratio of degree-1 quantities),
% f is degree 0 in xi. And r = a(1-e*cos(E)) where a = L^2/mu is degree 2
% and (1-e*cos(E)) is degree 0. So r is degree 2 in xi!
%
% Therefore r is NOT degree 0! I was wrong earlier!
% D(r) ≠ 0 in general!
%
% Wait, but Eq. (13b) says Dr = 0. Let me reconsider.
%
% The resolution: there are TWO ways to view r.
% (1) As a function of (xi, eta): r = |eta|. This is degree 0 in xi. D(r) = 0.
% (2) As a function of Delaunay variables: r = a(1-e*cos(E)) = (L^2/mu)(1-e*cos(E)).
%     Here L is degree 1, so r appears to be degree 2.
%
% The D operator is defined as -sum xi_k d/d(xi_k). When applied to r = |eta|,
% we get Dr = 0 because r doesn't depend on xi (view 1).
%
% When we write r = a(1-e*cos(E)) in Delaunay variables, this is the SAME
% function r = |eta| evaluated via the Delaunay-to-Cartesian map. The fact
% that a = L^2/mu and L depends on xi doesn't change the fact that
% |eta| doesn't depend on xi.
%
% So a = L^2/mu is degree 2, and (1-e*cos(E)) adjusts so that
% r = a(1-e*cos(E)) is degree 0 overall. The individual factors are NOT
% separately homogeneous in a useful way.
%
% This means I CANNOT simply assign degree 2 to a and degree 0 to (1-e*cos E)
% and apply D term by term. The product a*(1-e*cos(E)) has a combined
% degree that is NOT the sum of individual degrees.
%
% The correct approach: D acts on functions of (xi, eta) directly.
% When we write S1 in Delaunay variables, we are composing the function
% with the Delaunay-to-Cartesian map. D then acts through the chain rule.
%
% So: D(S1) requires knowing how L, G, H depend on xi (which we know:
% DL_j = -(p_j), from Eq. 12) and how the degree-0 quantities (e, f, etc.)
% are affected (they are NOT affected, since D annihilates them).
%
% But when I write D(a^3/r^3), I need to be careful:
% a^3/r^3 = L^6/(mu^3*r^3) where r = |eta| is degree 0.
% So a^3/r^3 = (L^6/mu^3) * (1/r^3).
% L^6 is degree 6 in xi. 1/r^3 is degree 0.
% So a^3/r^3 IS degree 6, and D(a^3/r^3) = -6*(a^3/r^3).
%
% Similarly: a/r = L^2/(mu*r) is degree 2, D(a/r) = -2*(a/r).
%
% And: a^2/r^2 = L^4/(mu^2*r^2) is degree 4, D(a^2/r^2) = -4*(a^2/r^2).
%
% This works because r = |eta| is truly degree 0, so 1/r^n is degree 0,
% and the homogeneity comes entirely from the L factors in a^n.

% CORRECTED APPROACH:
% Rewrite dS1/dl with EXPLICIT L,G,H dependence (using a = L^2/mu):
%
% dS1/dl = (mu^2 k2 / L^3) * [A * (L^6/(mu^3 r^3) - L^3/G^3)
%           + B * L^6/(mu^3 r^3) * cos(2g+2f)]
%
% = k2*A*L^3/(mu*r^3) - mu^2*k2*A/G^3 + k2*B*L^3*cos(2g+2f)/(mu*r^3)
%
% Now apply D term by term. Each term has the form (degree-0 stuff) * L^m * G^n:
%
% T1 = [k2*A/(mu*r^3)] * L^3          → degree 3, D(T1) = -3*T1
% T2 = [-mu^2*k2*A]     * G^{-3}      → degree -3, D(T2) = +3*T2
% T3 = [k2*B*cos(2g+2f)/(mu*r^3)] * L^3 → degree 3, D(T3) = -3*T3
%
% So D(dS1/dl) = -3*T1 + 3*T2 - 3*T3
%
% This is CORRECT because r, f, g, theta, e are all degree-0 in xi.
% The factors A, B only depend on theta = H/G (degree 0).
% The trig functions cos(2g+2f) are degree 0.
% The only non-trivial homogeneity is from L and G.

printf('Correct homogeneity analysis:\n');
printf('  r = |eta| is degree 0 in xi (position, not velocity)\n');
printf('  a = L^2/mu, so a^n is degree 2n. a/r = L^2/(mu*r) is degree 2.\n');
printf('  L^3/G^3 = eta_orb^{-3} is degree 0 (ratio of degree-1 quantities).\n');
printf('  f, g, e, theta are all degree 0.\n\n');

% Now compute D(dS1/dl):
D_dS1_dl_v2 = -3*T1 + 3*T2 - 3*T3;

% Factor out 3*mu^2*k2/L^3:
% -3*T1 = -3*k2*A*L^3/(mu*r^3) = -3*(mu^2*k2/L^3)*A*(L^6/(mu^3*r^3))
%       = -3*(mu^2*k2/L^3)*A*(a^3/r^3)
% +3*T2 = +3*(-mu^2*k2*A/G^3) = -3*(mu^2*k2/L^3)*A*(L^3/G^3)
%       Wait: +3*T2 = +3*(-mu^2*k2*A/G^3) = -3*mu^2*k2*A/G^3
%       And -3*(mu^2*k2/L^3)*A*(L^3/G^3) = -3*mu^2*k2*A/G^3. Yes, same.
%
% Hmm but the sign: +3*T2 = -3*mu^2*k2*A/G^3 < 0 (if A < 0).
% And -3*(mu^2*k2/L^3)*A*(L^3/G^3) = same. Good.

% So: D(dS1/dl) = 3*(mu^2*k2/L^3)*[-A*a^3/r^3 - A*L^3/G^3 - B*a^3/r^3*cos(2g+2f)]

% Wait, let me redo this carefully:
% D(dS1/dl) = -3*T1 + 3*T2 - 3*T3
%   = -3*k2*A*L^3/(mu*r^3) + 3*(-mu^2*k2*A/G^3) - 3*k2*B*L^3*cos(2g+2f)/(mu*r^3)

% Factor out -3*mu^2*k2/L^3 from each:
% From -3*k2*A*L^3/(mu*r^3): = -3*(mu^2*k2/L^3) * A*L^6/(mu^3*r^3)
%                              = -3*(mu^2*k2/L^3) * A*(a^3/r^3)
% From -3*mu^2*k2*A/G^3:     = -3*(mu^2*k2/L^3) * A*(L^3/G^3)
%                              Hmm: -3*mu^2*k2*A/G^3 = -3*(mu^2*k2/L^3)*(L^3/G^3)*A. No.
%                              -3*(mu^2*k2/L^3)*(L^3/G^3)*A = -3*mu^2*k2*A/G^3. Yes!
% From -3*k2*B*L^3*cos/(mu*r^3): = -3*(mu^2*k2/L^3) * B*(a^3/r^3)*cos(2g+2f)

% So: D(dS1/dl) = -3*(mu^2*k2/L^3) * [A*(a^3/r^3) + A*(L^3/G^3) + B*(a^3/r^3)*cos(2g+2f)]

% Hmm wait, the T2 term: +3*T2 = +3*(-mu^2*k2*A/G^3) = -3*mu^2*k2*A/G^3.
% And I'm factoring out -3*(mu^2*k2/L^3), so:
% -3*mu^2*k2*A/G^3 = -3*(mu^2*k2/L^3) * (L^3/G^3) * A
% = -3*(mu^2*k2/L^3) * A * L^3/G^3. YES.

% So ALL three terms get -3*(mu^2*k2/L^3) factored out, and we get:
% D(dS1/dl) = -3*(mu^2*k2/L^3) * [A*a^3/r^3 + A*L^3/G^3 + B*a^3/r^3*cos(2g+2f)]

% Hmm, but that's:
% -3*(mu^2*k2/L^3) * [A*(a^3/r^3 + L^3/G^3) + B*a^3/r^3*cos(2g+2f)]

% But delta_p_1 from BH61 Eq. (14) is:
% 3*(mu^2*k2/L^3) * {A*[-eta^{-3} + a^3/r^3*(1-2a/r)] + B*a^3/r^3*(1-2a/r)*cos(2g+2f)}
%
% where eta^{-3} = L^3/G^3.
%
% My result: -3*(mu^2*k2/L^3) * [A*(a^3/r^3 + L^3/G^3) + B*a^3/r^3*cos(2g+2f)]
% = 3*(mu^2*k2/L^3) * [A*(-a^3/r^3 - L^3/G^3) - B*a^3/r^3*cos(2g+2f)]
% = 3*(mu^2*k2/L^3) * {A*[-L^3/G^3 - a^3/r^3] + B*[-a^3/r^3*cos(2g+2f)]}
%
% But BH61 has:
% 3*(mu^2*k2/L^3) * {A*[-L^3/G^3 + a^3/r^3*(1-2a/r)] + B*[a^3/r^3*(1-2a/r)*cos(2g+2f)]}
%
% These are DIFFERENT. BH61 has (1-2a/r) factors and my result doesn't.
%
% The (1-2a/r) must come from delta_p_1 including BOTH D(dS1/dl1) contribution
% AND D(dS1/dl1) picks up extra terms because dS1/dl involves f and f depends on l.
%
% WAIT: I think the issue is that delta_p_1 = D(dS1/dl_1 + dS1*/dl_1).
% For S1, dS1/dl_1 = dS1/dl. I computed this as (mu^2k2/L^3)(A*sigma1+B*sigma2).
% But actually, f depends on l (through Kepler's equation), so when I
% differentiate S1 w.r.t. l, the derivative acts on both the explicit l
% AND on f(l). The result from Brouwer Eq. (13) already accounts for this —
% it's the full partial derivative dS1/dl including the f(l) chain rule.
%
% So my dS1/dl is correct. The issue must be elsewhere.
%
% Let me reconsider. Maybe delta_p_1 ≠ D(dS1/dl_1).
% From Eq. (11a): delta_p_j = D(dS_1/dl_j + dS_1*/dl_j).
% For j=1: delta_p_1 = D(dS1/dl) + D(dS1*/dl).
% Since S1* doesn't depend on l_1, dS1*/dl = 0. So delta_p_1 = D(dS1/dl).
%
% My result is D(dS1/dl) and it doesn't match BH61. So either:
% (a) My homogeneity analysis is wrong, or
% (b) BH61 Eq. (14) has OCR errors, or
% (c) I'm missing something about how D acts on the a^3/r^3 terms.
%
% Let me reconsider whether a^3/r^3 really has a clean homogeneity degree.
% a = L^2/mu, so a depends on L (degree 1 in xi). But a is conventionally
% thought of as a function of the energy, which involves |xi|^2.
% When I write a = L^2/mu, I'm using the Delaunay relation. L depends on
% both xi and eta. So a = L(xi,eta)^2/mu is NOT simply degree 2 in xi.
%
% This is the same subtlety as L_1. L_1 = mu/sqrt(2mu/r - v^2).
% DL_1 ≠ -L_1 (L_1 is not degree 1).
% Therefore Da ≠ -2a either!
%
% I was WRONG to assign degree 2 to a. Let me reconsider.

printf('\n*** CRITICAL ERROR IN HOMOGENEITY ANALYSIS ***\n');
printf('a = L^2/mu where L = L_1 is NOT degree-1 homogeneous in xi.\n');
printf('L_1 depends on both |xi|^2 and r, so it is not homogeneous.\n');
printf('Therefore a is not degree-2, and the simple D(L^3) = -3*L^3 trick\n');
printf('does NOT apply to L = L_1.\n\n');

printf('However, L_2 = G and L_3 = H ARE degree-1 (angular momentum).\n');
printf('So D(G^n) = -n*G^n and D(H^n) = -n*H^n are valid.\n\n');

printf('For terms involving L = L_1, we must use the chain rule:\n');
printf('  D(f(L,G,H,...)) = (DL)*df/dL + (DG)*df/dG + (DH)*df/dH + 0\n');
printf('where DL = -(p_1), DG = -(p_2) = -G, DH = -(p_3) = -H.\n\n');

printf('This is the correct approach. Restarting the computation.\n\n');

% CORRECT APPROACH:
% For any function f(L_1, L_2, L_3, l_1, l_2, l_3) viewed through the
% Delaunay-to-Cartesian map:
%
% Df = -sum_k xi_k df/d(xi_k)
%    = -sum_k xi_k [sum_j (df/dL_j)(dL_j/d(xi_k)) + sum_j (df/dl_j)(dl_j/d(xi_k))]
%    = sum_j (df/dL_j)*(-sum_k xi_k dL_j/d(xi_k)) + sum_j (df/dl_j)*(-sum_k xi_k dl_j/d(xi_k))
%    = sum_j (df/dL_j)*(DL_j) + sum_j (df/dl_j)*(Dl_j)
%    = sum_j (df/dL_j)*(-(p_j)) + sum_j (df/dl_j)*(+(q_j))
%
% This is the GENERAL chain rule for D.
%
% For delta_p_1 = D(dS1/dl):
% Let F = dS1/dl. Then:
% D(F) = -(p_1)*dF/dL - (p_2)*dF/dG - (p_3)*dF/dH + (q_1)*dF/dl + (q_2)*dF/dg + (q_3)*dF/dh
%
% Since q_3 = 0 (file 07) and S1 doesn't depend on h:
% D(F) = -(p_1)*dF/dL - G*dF/dG - H*dF/dH + (q_1)*dF/dl + (q_2)*dF/dg
%
% where F = dS1/dl = (mu^2 k2 / L^3) * (A*sigma1 + B*sigma2).
%
% This requires computing 5 partial derivatives of F w.r.t. L, G, H, l, g.
% This is a substantial computation. Let me set it up symbolically.

printf('Setting up D(dS1/dl) via chain rule through Delaunay variables.\n\n');

% F = dS1/dl = (mu^2 k2 / L^3) * [A*(a^3/r^3 - L^3/G^3) + B*(a^3/r^3)*cos(2g+2f)]
%
% Note: a = L^2/mu, so a depends on L. And f depends on L, G, l through
% the Kepler equation and true anomaly relation.
%
% But r and f are functions of (a, e, l), and a = L^2/mu, e = sqrt(1-G^2/L^2).
% So r and f depend on L, G, l.

% For symbolic differentiation, I need to express F entirely in terms of
% L, G, H, l, g (with r and f as implicit functions of L, G, l).
%
% This is getting complex. Let me use Octave's symbolic engine to do the
% differentiations, treating r and f as explicit symbols and applying
% the chain rule for dr/dL, dr/dG, dr/dl, df/dL, df/dG, df/dl.

printf('This computation requires the chain rule through r(L,G,l) and f(L,G,l).\n');
printf('Setting up the partial derivatives of r and f...\n\n');

% From Kepler's equation l = E - e*sin(E) and r = a(1-e*cos(E)):
% dr/dl = a*e*sin(E) * dE/dl = a*e*sin(E) / (1-e*cos(E)) = ae*sin(E)/D_kep
%       = a*e*sin(f)*sqrt(1-e^2)/(1-e^2) ... this is getting complicated.
%
% Actually, the standard results are:
% dr/dl = a*e*sin(f)/sqrt(1-e^2) * ... no, let me use the known identity:
% dr/df = a*e*sin(f)*(1-e^2)/(1+e*cos(f))^2... this is very messy.
%
% Perhaps it's better to use the D operator identities DIRECTLY:
% We know from file 07 / BH61:
%   Df = ??? (need to determine)
%   D(a/r) = ??? (need to determine)
%
% Actually, BH61 derives specific D-identities. The key ones are:
%   D(r) = 0  (Eq. 13b, since r = |eta|)
%   D(f) = ?? We DON'T know this yet.
%   D(a/r) = r = |eta| is degree 0, a = L^2/mu... wait.
%
% STOP. I'm going in circles. Let me step back and think about what approach
% BH61 actually uses.
%
% Looking at Eq. (14), the result has terms like (1-2a/r). This factor
% arises from D acting on a^3/r^3: if we write a^3/r^3 = (a/r)^3 and
% note that a depends on L while r is degree 0:
%
%   D(a^3/r^3) = D(L^6/(mu^3*r^3))
%
% But L = L_1 is NOT degree 1! So D(L^6) ≠ -6*L^6.
%
% I need to compute D(L_1) properly. From Eq. (12): DL_1 = -(p_1).
% From file 07: p_1 = ...?

printf('Need to determine D(L_1) = -(p_1). Reading p_1 from file 07.\n');
printf('p_1 = sum_k xi_k d(eta_k)/dl_1 from the explicit computation.\n\n');

printf('*** PAUSING: this derivation requires the explicit p_1 value ***\n');
printf('*** and the D-operator chain rule through L, G, H, l, g.    ***\n');
printf('*** This is a multi-step computation that should be done     ***\n');
printf('*** carefully in a separate script.                          ***\n');
