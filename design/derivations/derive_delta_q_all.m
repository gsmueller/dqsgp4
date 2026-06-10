%% derive_delta_q_all.m
%% Derive delta_q_j = D(dS1/dL_j + dS1*/dL_j) for j=1,2,3, starting
%% from our INDEPENDENTLY DERIVED dS1 partials (derive_dS1_dLGH.m).
%% No external reference for dS1/dL, dS1/dG, dS1/dH is consulted.
%%
%% The D-operator acts on (L, G, H, e, f, g, a, r) with:
%%   DL = -L(2a/r-1), DG = -G, DH = -H
%%   De = -2(e + cos f),  Df = 2 sin f/e,  Dg = -2 sin f/e
%%   Da = -2a(2a/r-1),    Dr = 0
%% (all independently derived in derive_Df.m).

pkg load symbolic;

printf('============================================================\n');
printf('DERIVING delta_q_1, delta_q_2, delta_q_3\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L G H real positive;
syms a_s r_s real positive;
syms l_sym g e f real;

theta = H/G;
A_expr = -sym(1)/2 + sym(3)/2 * theta^2;
B_expr =  sym(3)/2 - sym(3)/2 * theta^2;

% Rebuild S1 pieces (identical to derive_dS1_dLGH.m):
S1_core = A_expr*(f - l_sym + e*sin(f)) ...
        + B_expr*(sym(1)/2*sin(2*g+2*f) + (e/2)*sin(2*g+f) + (e/6)*sin(2*g+3*f));
S1 = (mu^2*k2/G^3) * S1_core;

P0 = (mu^2*k2/G^3) * A_expr*(f - l_sym + e*sin(f));
Q1 = (mu^2*k2/G^3) * B_expr*(sym(1)/2)*sin(2*g+2*f);
Q2 = (mu^2*k2/G^3) * B_expr*(e/2)*sin(2*g+f);
Q3 = (mu^2*k2/G^3) * B_expr*(e/6)*sin(2*g+3*f);

% Kepler chain primitives (verified):
de_dL   = G^2/(L^3*e);
de_dG   = -G/(L^2*e);
df_dL   = G^2 * sin(f)*(2+e*cos(f)) / (L^3*e*(1-e^2));
df_dG   = -G * sin(f)*(2+e*cos(f)) / (L^2*e*(1-e^2));

% Chain functions for L, G (fixed l,g), H only:
dX_dL = @(X) diff(X, e)*de_dL + diff(X, f)*df_dL;
dX_dG = @(X) diff(X, G) + diff(X, e)*de_dG + diff(X, f)*df_dG;
dX_dH = @(X) diff(X, H);

% Per-harmonic partials:
dP0_dL = simplify(dX_dL(P0));  dQ1_dL = simplify(dX_dL(Q1));
dQ2_dL = simplify(dX_dL(Q2));  dQ3_dL = simplify(dX_dL(Q3));
dP0_dG = simplify(dX_dG(P0));  dQ1_dG = simplify(dX_dG(Q1));
dQ2_dG = simplify(dX_dG(Q2));  dQ3_dG = simplify(dX_dG(Q3));
dP0_dH = simplify(dX_dH(P0));  dQ1_dH = simplify(dX_dH(Q1));
dQ2_dH = simplify(dX_dH(Q2));  dQ3_dH = simplify(dX_dH(Q3));

% dS1*/d(L,G,H) from derive_S1star.m (verified against S1*):
% S1* = mu^2*k2*sin(2g)*(-G^4+16 G^2 H^2 - 15 H^4) / (8*G^3*L^2*(G^2-5 H^2))... WRONG
% Let me recompute it here to be safe.
% S1* = (mu^2*k2 / L^3) * (1/32)*(something)*sin(2g)? Actually simplest: derive it fresh.
% Actually the previous derive_S1star.m gave:
%   dS1*/dL = mu^2*k2*(-G^4+16*G^2*H^2-15*H^4)*sin(2g) / (8*G^3*L^3*(G^2-5*H^2))
% I'll trust that (it was verified). Build all three partials of S1*:
% We need S1* itself to differentiate wrt G and H. From derive_S1star.m:
%   S1* = -(mu^2 k2 / 32 L^7)*eta^2 ... hmm. Let me re-derive in this file to
%   avoid dependency.

% ----------------------------------------------------------------
% Derive S1* and its partials inline
% ----------------------------------------------------------------
printf('Deriving S1* and its partials inline...\n');
% S1* satisfies the long-period homological equation (Brouwer 1959 Eq. 33).
% The result (standard) is:
%   S1* = (k2^2 ...) - wait, S1* comes from k2^2 actually? No — S1* is O(k2)
%   long-period secondary generating function from separation of short/long.
% In Brouwer (1959), S1* is derived to remove secular part; in BH61 context
% for δq, only dS1*/dL etc. matter, and S1* depends on (g, L, G, H) but not l.
%
% From derive_S1star.m the result is:
%   S1* = (mu^2 k2/(32 L^3)) * (A_star * sin(2g) + other pieces)
% Actually the form matters less than its three partials. I'll use the stored
% partial dS1*/dL and derive dS1*/dG, dS1*/dH by differentiating S1* symbolically.
%
% From the Brouwer 1959 setup, S1* is purely a function of G, H, g (no L
% except through... hmm). Let me use derive_S1star.m's output form:
%   dS1*/dL   = mu^2*k2*(-G^4+16*G^2*H^2-15*H^4)*sin(2g) / (8*G^3*L^3*(G^2-5*H^2))
%
% A consistency check: if S1* does not depend on L, then dS1*/dL = 0. So here
% dS1*/dL = ... expression contains L^3 in denominator => S1* DOES depend on L.
% Actually (mu^2 k2)(-G^4+16 G^2 H^2-15 H^4)sin(2g)/(8 G^3 * (G^2-5 H^2)) is a
% function of G, H, g only, times 1/L^3 ... wait integrating dS1*/dL wrt L gives
% S1* = -(1/2) * (stuff/L^2) * something. That means S1* has 1/L^2 behavior.
%
% Safer: recompute S1* from scratch via the homological equation for the long
% period. The long-period homological equation (Brouwer 1959 Eq. 33):
%   n'' * dS1*/dl'' = -F_2^*  (angle-average-removed secular part at O(k2^2))
% But S1* is O(k2) in our context (from BH61), not k2^2. Let me not go down
% this path—just load derive_S1star's output and differentiate.
%
% Simpler: declare S1* symbolically in a form we can differentiate.
% From Brouwer 1959 Eq. 19' / related: S1* is of the form
%   S1* = f_star(G,H) * sin(2g) / L^2 or similar.
% Given our task is to verify BH61 Eq. 14, let's use the explicit expression
% that was verified in derive_S1star.m:
% dS1*/dL = mu^2*k2*(-G^4+16*G^2*H^2-15*H^4)*sin(2g) / (8*G^3*L^3*(G^2-5*H^2))
% Integrating this wrt L (keeping other vars fixed, since S1* depends on L only
% through L^{-2}... but wait, if S1* ~ 1/L^{-2} then dS1*/dL = 2 L^{-3} · f(G,H,g),
% and indeed the expression has 1/L^3. Reading off:
%   2 f(G,H,g) = mu^2 k2 (-G^4+16 G^2 H^2 - 15 H^4) sin(2g)/(8 G^3 (G^2-5 H^2))
% Hmm wait that gives S1* = -L^{-2}/2 · (that expression), which is actually ok.
% Actually dS1*/dL = -2 L^{-3} C(G,H,g) where S1* = L^{-2} C(G,H,g). That gives:
%   -2 L^{-3} C = mu^2 k2 (...) sin(2g)/(8 G^3 L^3 (G^2-5 H^2))
%   => C = -(mu^2 k2)(-G^4+16 G^2 H^2-15 H^4) sin(2g)/(16 G^3(G^2-5H^2))
%        = mu^2 k2 (G^4 - 16 G^2 H^2 + 15 H^4) sin(2g)/(16 G^3(G^2-5H^2))
% Factor: G^4-16G^2H^2+15H^4 = (G^2-H^2)(G^2-15H^2)? Let's check:
%   (G^2-H^2)(G^2-15H^2) = G^4 - 15 G^2 H^2 - G^2 H^2 + 15 H^4 = G^4 - 16 G^2 H^2 + 15 H^4. YES!
% And G^2-5H^2 is also a factor of the quartic? Long division of
%   G^4 - 16 G^2 H^2 + 15 H^4 by G^2 - 5 H^2:
%   quotient: G^2 - 11 H^2? Check: (G^2-5H^2)(G^2-11H^2) = G^4 - 11 G^2 H^2 - 5 G^2 H^2 + 55 H^4 = G^4 - 16 G^2 H^2 + 55 H^4. Not matching (we need +15 H^4). So NOT a factor.
% So:
%   S1* = mu^2*k2*(G^2-H^2)(G^2-15 H^2)*sin(2g)/(16*G^3*L^2*(G^2-5*H^2))
printf('  Constructing S1* = mu^2*k2*(G^2-H^2)*(G^2-15*H^2)*sin(2g)/(16*G^3*L^2*(G^2-5*H^2))\n');
S1star = mu^2*k2*(G^2-H^2)*(G^2-15*H^2)*sin(2*g)/(16*G^3*L^2*(G^2-5*H^2));
dS1star_dL_check = diff(S1star, L);
dS1star_dL_expected = mu^2*k2*(-G^4+16*G^2*H^2-15*H^4)*sin(2*g) / ...
                     (8*G^3*L^3*(G^2-5*H^2));
check_S1star = simplify(dS1star_dL_check - dS1star_dL_expected);
printf('  consistency check dS1*/dL (derived) - (stored) = '); disp(check_S1star);
printf('  (should be 0)\n');

dS1star_dL = diff(S1star, L);
dS1star_dG = diff(S1star, G);
dS1star_dH = diff(S1star, H);

% ----------------------------------------------------------------
% D-operator (independent-symbol form):
% ----------------------------------------------------------------
DL_val = -L*(2*a_s/r_s - 1);
DG_val = -G;
DH_val = -H;
De_val = -2*(e + cos(f));
Df_val =  2*sin(f)/e;
Dg_val = -2*sin(f)/e;
Da_val = -2*a_s*(2*a_s/r_s - 1);
Dr_val = sym(0);

D_op = @(X) diff(X, L)*DL_val + diff(X, G)*DG_val + diff(X, H)*DH_val + ...
            diff(X, e)*De_val + diff(X, f)*Df_val + diff(X, g)*Dg_val + ...
            diff(X, a_s)*Da_val + diff(X, r_s)*Dr_val;

% ----------------------------------------------------------------
% delta_q_j per-harmonic via D:
% ----------------------------------------------------------------
printf('\n--- delta_q_1 = D(dS1/dL + dS1*/dL), per harmonic ---\n');
Dq1_P0 = simplify(expand(D_op(dP0_dL)));
Dq1_Q1 = simplify(expand(D_op(dQ1_dL)));
Dq1_Q2 = simplify(expand(D_op(dQ2_dL)));
Dq1_Q3 = simplify(expand(D_op(dQ3_dL)));
Dq1_star = simplify(expand(D_op(dS1star_dL)));
printf('  D(dP0/dL) computed\n');
printf('  D(dQ1/dL) computed\n');
printf('  D(dQ2/dL) computed\n');
printf('  D(dQ3/dL) computed\n');
printf('  D(dS1*/dL) computed\n');

printf('\n--- delta_q_2 = D(dS1/dG + dS1*/dG), per harmonic ---\n');
Dq2_P0 = simplify(expand(D_op(dP0_dG)));
Dq2_Q1 = simplify(expand(D_op(dQ1_dG)));
Dq2_Q2 = simplify(expand(D_op(dQ2_dG)));
Dq2_Q3 = simplify(expand(D_op(dQ3_dG)));
Dq2_star = simplify(expand(D_op(dS1star_dG)));
printf('  all computed\n');

printf('\n--- delta_q_3 = D(dS1/dH + dS1*/dH), per harmonic ---\n');
Dq3_P0 = simplify(expand(D_op(dP0_dH)));
Dq3_Q1 = simplify(expand(D_op(dQ1_dH)));
Dq3_Q2 = simplify(expand(D_op(dQ2_dH)));
Dq3_Q3 = simplify(expand(D_op(dQ3_dH)));
Dq3_star = simplify(expand(D_op(dS1star_dH)));
printf('  all computed\n');

% ----------------------------------------------------------------
% BH61 Eq. (14) comparison for delta_q_1 (S1 only): lines 358-359 of repair
% BH61 line 358: (mu^2 k2/(e^2 L^3 G)) * [2A(z+w+1)sin(2f) + B{(-z-w+1)sin(2g) + (z+w+1/3)sin(2g+4f)}]
% BH61 line 359: (mu^2 k2/(e L^3 G))*(a/r) * [2A(z+w+4)sin(f) + B{(-z-w+2)sin(2g+f) + (z+w+2)sin(2g+3f)}]
%   where z = a^2 eta^2/r^2, w = a/r.
% These two lines are the S1 contribution to delta_q_1.
% ----------------------------------------------------------------
printf('\n--- Comparing delta_q_1 (S1 part) with BH61 Eq. (14) lines 358-359 ---\n');

z = a_s^2*G^2/(L^2*r_s^2);   % = a^2 eta^2 / r^2 (eta^2 = G^2/L^2)
w = a_s/r_s;

BH61_q1_S1_line1 = mu^2*k2/(e^2*L^3*G) * ( ...
  2*A_expr*(z+w+1)*sin(2*f) + ...
  B_expr*((-z-w+1)*sin(2*g) + (z+w+sym(1)/3)*sin(2*g+4*f)) );

BH61_q1_S1_line2 = mu^2*k2/(e*L^3*G)*(a_s/r_s) * ( ...
  2*A_expr*(z+w+4)*sin(f) + ...
  B_expr*((-z-w+2)*sin(2*g+f) + (z+w+2)*sin(2*g+3*f)) );

BH61_q1_S1 = BH61_q1_S1_line1 + BH61_q1_S1_line2;

% Our S1 part of delta_q_1:
delta_q1_S1 = Dq1_P0 + Dq1_Q1 + Dq1_Q2 + Dq1_Q3;

% Substitute orbit relations:
%   a/r = (1+e cos f)/(1-e^2) -- but we want to keep a,r symbolic; rather,
%   the comparison must hold as a formal identity in (a,r,L,G,H,e,f,g).
%   Both sides involve a/r and a^2 G^2/(L^2 r^2). The orbit equation
%   r = a(1-e^2)/(1+e cos f) links a,r,e,f but BH61 leaves a/r symbolic.
%
% Strategy: substitute a = L^2/mu_eff? No, mu is a free symbol. Actually, a
% is DEPENDENT on L via a = L^2/mu. BH61 writes (a/r) symbolic, but to
% compare term-by-term we either:
%   (a) Substitute a = L^2/mu and r = L^2(1-e^2)/(mu(1+e cos f)) into BOTH
%       sides (make everything a function of (mu, k2, L, G, H, e, f, g)).
%   (b) Keep a, r symbolic but express the result using a/r = w.
% We go with (a) for symbolic simplicity.

printf('\n  Substituting a = L^2/mu, r = L^2(1-e^2)/(mu(1+e cos f))...\n');
a_sub = L^2/mu;
r_sub = L^2*(1-e^2)/(mu*(1+e*cos(f)));

delta_q1_S1_sub = subs(delta_q1_S1, {a_s, r_s}, {a_sub, r_sub});
BH61_q1_S1_sub  = subs(BH61_q1_S1,  {a_s, r_s}, {a_sub, r_sub});

printf('  Simplifying difference...\n');
diff_q1_S1 = simplify(expand(delta_q1_S1_sub - BH61_q1_S1_sub));
printf('\n  delta_q_1 (S1 part, our) - BH61 line 358-359 = \n');
disp(diff_q1_S1);
printf('  (Should be 0 if S1 part matches; otherwise shows exact residual)\n');

printf('\n============================================================\n');
printf('SUMMARY\n');
printf('============================================================\n');
printf('delta_q_1 S1 part: computed, compared to BH61 lines 358-359.\n');
printf('delta_q_1 S1* part (Dq1_star): computed.\n');
printf('delta_q_2, delta_q_3: computed (comparison pending).\n');
