%% check_delta_q_all_numerical.m
%% Numerical verification of delta_q_1, delta_q_2, delta_q_3 vs BH61 Eq. (14)
%% lines 358-376 of the repaired document.
%%
%% Strategy: compute D(dS1/dL_j) + D(dS1*/dL_j) symbolically from our
%% INDEPENDENTLY DERIVED partials, then evaluate numerically at several
%% self-consistent orbital test points. Compare against numerical
%% evaluation of BH61 expressions (as written in the paper).
%%
%% NOTE: test points satisfy G = L sqrt(1-e^2), 0 <= H <= G.

pkg load symbolic;

syms mu k2 real positive;
syms L G H real positive;
syms a_s r_s real positive;
syms l_sym g e f real;

theta = H/G;
A_expr = -sym(1)/2 + sym(3)/2 * theta^2;
B_expr =  sym(3)/2 - sym(3)/2 * theta^2;
eta_sym = G/L;

%% =================================================================
%% Build S1 and its partials (as in derive_dS1_dLGH.m)
%% =================================================================
P0 = (mu^2*k2/G^3) * A_expr*(f - l_sym + e*sin(f));
Q1 = (mu^2*k2/G^3) * B_expr*(sym(1)/2)*sin(2*g+2*f);
Q2 = (mu^2*k2/G^3) * B_expr*(e/2)*sin(2*g+f);
Q3 = (mu^2*k2/G^3) * B_expr*(e/6)*sin(2*g+3*f);
S1 = P0 + Q1 + Q2 + Q3;

de_dL = G^2/(L^3*e);       de_dG = -G/(L^2*e);
df_dL = G^2*sin(f)*(2+e*cos(f))/(L^3*e*(1-e^2));
df_dG = -G*sin(f)*(2+e*cos(f))/(L^2*e*(1-e^2));

dS1_dL = diff(S1,e)*de_dL + diff(S1,f)*df_dL;
dS1_dG = diff(S1,G)      + diff(S1,e)*de_dG + diff(S1,f)*df_dG;
dS1_dH = diff(S1,H);

%% =================================================================
%% S1* = mu^2 k2 (G^2-H^2)(G^2-15 H^2) sin(2g) / (16 G^3 L^2 (G^2-5H^2))
%% =================================================================
S1star = mu^2*k2*(G^2-H^2)*(G^2-15*H^2)*sin(2*g)/(16*G^3*L^2*(G^2-5*H^2));
dS1star_dL = diff(S1star, L);
dS1star_dG = diff(S1star, G);
dS1star_dH = diff(S1star, H);

%% =================================================================
%% D operator
%% =================================================================
DL_val = -L*(2*a_s/r_s - 1);
DG_val = -G; DH_val = -H;
De_val = -2*(e + cos(f));
Df_val =  2*sin(f)/e;
Dg_val = -2*sin(f)/e;
Da_val = -2*a_s*(2*a_s/r_s - 1);
D_op = @(X) diff(X,L)*DL_val + diff(X,G)*DG_val + diff(X,H)*DH_val + ...
            diff(X,e)*De_val + diff(X,f)*Df_val + diff(X,g)*Dg_val + ...
            diff(X,a_s)*Da_val;

delta_q1 = D_op(dS1_dL) + D_op(dS1star_dL);
delta_q2 = D_op(dS1_dG) + D_op(dS1star_dG);
delta_q3 = D_op(dS1_dH) + D_op(dS1star_dH);

delta_q1_S1only = D_op(dS1_dL);
delta_q1_S1star = D_op(dS1star_dL);

%% =================================================================
%% BH61 expressions from the repaired document
%% =================================================================
z = a_s^2*eta_sym^2/r_s^2;   % a^2 eta^2 / r^2
w = a_s/r_s;

% --- delta_q_1 ---
BH61_q1_line358 = mu^2*k2/(e^2*L^3*G) * ( ...
  (-1+3*theta^2)*(z+w+1)*sin(2*f) + ...
  (sym(3)/2-sym(3)/2*theta^2)*((-z-w+1)*sin(2*g) + (z+w+sym(1)/3)*sin(2*g+4*f)) );

BH61_q1_line359 = mu^2*k2/(e*L^3*G)*(a_s/r_s) * ( ...
  (-1+3*theta^2)*(z+w+4)*sin(f) + ...
  (sym(3)/2-sym(3)/2*theta^2)*((-z-w+2)*sin(2*g+f) + (z+w+2)*sin(2*g+3*f)) );

BH61_q1_line360 = mu^2*k2/(e*L^3*G) * ( ...
  (sym(1)/4*(1-11*theta^2) - 10*theta^4/(1-5*theta^2)) ) * ...
  ( (1 - 3*a_s/r_s)*e*sin(2*g) + sin(2*g+f) - sin(2*g-f) );

BH61_q1_full = BH61_q1_line358 + BH61_q1_line359 + BH61_q1_line360;
BH61_q1_S1   = BH61_q1_line358 + BH61_q1_line359;
BH61_q1_S1s  = BH61_q1_line360;

% --- delta_q_2: lines 361-374 (note: only the Eq. 14 body lines 361-363
% belong to the S1+S1* contribution; 369-374 are additional S1* pieces).
% BH61's delta_q_2 Eq. (14) starts at line 361 and continues through 374.
% Let me transcribe all six terms:
% Line 361: -mu^2 k2/(e^2 L^2 G^2) * { (-1+3θ²)(z+w+1) sin 2f + B*[...] }
%   (same angular structure as q1 line 358, but prefactor -mu^2 k2/(e^2 L^2 G^2))
BH61_q2_line361 = -mu^2*k2/(e^2*L^2*G^2) * ( ...
  (-1+3*theta^2)*(z+w+1)*sin(2*f) + ...
  (sym(3)/2-sym(3)/2*theta^2)*((-z-w+1)*sin(2*g) + (z+w+sym(1)/3)*sin(2*g+4*f)) );

BH61_q2_line362 = -mu^2*k2/(e*L^2*G^2)*(a_s/r_s) * ( ...
  (-1+3*theta^2)*(z+w+4)*sin(f) + ...
  (sym(3)/2-sym(3)/2*theta^2)*((-z-w+2)*sin(2*g+f) + (z+w+2)*sin(2*g+3*f)) );

% Line 363: -mu^2 k2/(e L^2 G^2) (1 - a/r) * {(-1+3θ²)(z+w+1) sin f + B*[...]}
BH61_q2_line363 = -mu^2*k2/(e*L^2*G^2)*(1 - a_s/r_s) * ( ...
  (-1+3*theta^2)*(z+w+1)*sin(f) + ...
  (sym(3)/2-sym(3)/2*theta^2)*((-z-w-1)*sin(2*g+f) + (z+w+sym(1)/3)*sin(2*g+3*f)) );

% Page 197 continuation (lines 369-374):
% Need "sin u" — u = f + g (argument of latitude, not f) per BH61 usage.
% Actually in BH61 notation u is often f+g. Let me check: "2(f-l)+e(sin f - sin u)".
% Let me tentatively try u = E (eccentric anomaly). Then sin E = (r/a) sin(f)/eta = ...
% Actually this part (P197 continuation) is for delta_q_2 S1* contributions.
% We skip the P197 piece in this first-pass check — the S1+S1* of our
% derivation does NOT include the k2^2 or (a/r)-expanded pieces that
% appear in the P197 continuation, which come from third-order or higher
% corrections. Since we want to verify the k2^1 drag-coupling delta_q_j,
% we compare against lines 361-363 only.
BH61_q2_eq14_primary = BH61_q2_line361 + BH61_q2_line362 + BH61_q2_line363;

% --- delta_q_3 (lines 375-376), where u is argument of latitude = f+g? ---
% "sin f - sin u" — with u = f+g would give "sin f - sin(f+g)", nonsense.
% Probably u = E (eccentric anomaly). Use sin E = sin(f)*eta/(1 + e cos f).
sin_E = sin(f)*eta_sym/(1 + e*cos(f));
sin_u = sin_E;  % interpret u as E

BH61_q3 = 6*mu^2*k2/G^4*theta*( ...
  2*(f - l_sym) + e*(sin(f) - sin_u) + (1-eta_sym)/e*sin(f) ...
  - sym(1)/3*sin(2*g+2*f) - (e/2)*sin(2*g+f) - (e/6)*sin(2*g+3*f) ) ...
  - 4*mu^2*k2/G^4*theta* ( sym(11)/8 + 10*theta^2/(1-5*theta^2) + 25*theta^4/(1-5*theta^2)^2 ) ...
  * e * sin(2*g+f);

%% =================================================================
%% Numerical tests
%% =================================================================
function v = mktest(mu_v, k2_v, L_v, e_v, theta_v, f_v, g_v, l_v)
  v.mu = mu_v; v.k2 = k2_v; v.L = L_v; v.e = e_v;
  v.G = L_v * sqrt(1 - e_v^2);
  v.H = v.G * theta_v;
  v.a = L_v^2/mu_v;
  v.r = v.a * (1 - e_v^2)/(1 + e_v*cos(f_v));
  v.f = f_v; v.g = g_v; v.l = l_v;
endfunction

tests = { mktest(1,    0.001, 1.2, 0.3, 0.8, 0.7, 1.5, 0.5), ...
          mktest(2,    0.005, 1.5, 0.1, 0.3, 2.1, 0.7, 1.2), ...
          mktest(3.98e5, -2.63e-3, 7000, 0.05, 0.6, 1.1, 2.3, 0.3) };

args = {mu, k2, L, G, H, e, f, g, l_sym, a_s, r_s};

function compare(label, ours_sym, bh61_sym, tests, args)
  printf('\n%s\n', label);
  for i = 1:length(tests)
    v = tests{i};
    vals = {v.mu, v.k2, v.L, v.G, v.H, v.e, v.f, v.g, v.l, v.a, v.r};
    o = double(subs(ours_sym, args, vals));
    b = double(subs(bh61_sym, args, vals));
    d = o - b;
    printf('  T%d: ours=%14.7g  bh61=%14.7g  diff=%13.3e  rel=%g\n', ...
      i, o, b, d, abs(d)/max(abs(o),abs(b)+eps));
  end
endfunction

compare('=== delta_q_1 (S1 only)  vs BH61 lines 358-359 ===', ...
        delta_q1_S1only, BH61_q1_S1, tests, args);
compare('=== delta_q_1 (S1* only) vs BH61 line 360       ===', ...
        delta_q1_S1star, BH61_q1_S1s, tests, args);
compare('=== delta_q_1 FULL       vs BH61 lines 358-360  ===', ...
        delta_q1, BH61_q1_full, tests, args);

compare('=== delta_q_2 FULL       vs BH61 lines 361-363  ===', ...
        delta_q2, BH61_q2_eq14_primary, tests, args);

compare('=== delta_q_3 FULL       vs BH61 lines 375-376  ===', ...
        delta_q3, BH61_q3, tests, args);
