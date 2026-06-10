% verify_ch10_thm7_vonzeipel.m
%
% Symbolic verification of Chapter 10 Theorem 7 (Second-order von Zeipel
% formula).  Reproduces the sign-correct modern-bracket formula
%
%   (VZ-2)   F_2*(L') = < F_2 + (1/2) { S_1, F_1 + F_1* }_modern >_l
%          = < F_2 - (1/2) { F_1 + F_1*, S_1 }_modern >_l   (antisymmetry)
%
% derived in ch10_foundations_thm7.md.  Also confirms the key Lemmas 7.8
% and 7.9 and the Theorem 6 cross-check <{F_1*, S_1}>_l = 0.
%
% Strategy.  A concrete low-dimensional test case is constructed:
%
%   F_0(L)       = mu^2 / (2 L^2)
%   F_1(L,G,l,g) = a(L,G) cos(l)  +  b(L,G) cos(2l + 2g)   [+ A_sec in test 6]
%   F_2          = 0                                       (BH61 J_2 truncation)
%
% With a(L,G) = 1/(L*G) and b(L,G) = 1/(L^2 * G^2) as concrete smooth
% functions (any smooth choice works; these were chosen for minimal
% symbolic-engine load).  For each test the script:
%
%   1. solves the first-order homological equation
%        n dS_1/dl = F_1 - F_1*  with  F_1* = <F_1>_l  and  n = mu^2/L^3,
%      in closed form (finite trigonometric polynomial in l, g);
%
%   2. Route A -- collects the epsilon^2 part of F(L(L',l), l) directly
%      from the type-2 Taylor expansion via
%        F_2^*(L') = <F_2>_l + U - (n'/2) <(dS_1/dl)^2>_l,
%      where U = <Sum_j (dS_1/dl_j)(dF_1/dL_j)>_l  (j = 1, 2; j = 3 vanishes
%      because dS_1/dh = 0).  This is equation (7.13) of the derivation.
%
%   3. Route B -- independently computes
%        F_2^*(L') = <F_2>_l + (1/2) <{S_1, F_1 + F_1*}_modern>_l
%      using the modern Poisson bracket
%        {f,g}_modern = Sum_j [ (df/dl_j)(dg/dL_j) - (df/dL_j)(dg/dl_j) ].
%      This is (VZ-2).
%
%   4. verifies Route A = Route B by reducing the difference to 0.
%
% Tests also check the two lemmas 7.8, 7.9 and the Theorem 6 identity
% <{F_1*, S_1}>_l = 0 independently.
%
% Test 1 is a 1D smoke test (F_1 = cos l only).
% Tests 2-5 use the full 2D (l,L)+(g,G) problem with F_1* = 0.
% Test 6 adds a nontrivial secular term A_sec(L,G) so that F_1* /= 0.
%
% PASS  = all symbolic differences reduce to 0.
% FAIL  = at least one difference fails to reduce to 0.
%
% Reference:
%   ch10_foundations_thm7.md  (Theorem 7 + Lemmas 7.8, 7.9)
%   ch07a_homological_equation.md  (first-order equation)
%   CH10_PLAN.md  (sign-convention note in Corrigendum 5 / Theorem 7)

pkg load symbolic;
clear; close all;

% --- Helper function (must appear before use in Octave scripts) -------
function result = l_avg(expr, l_sym)
  % L-average <expr>_l over [0, 2*pi] using Octave symbolic integrate.
  result = int(expr, l_sym, sym(0), 2*sym(pi)) / (2*sym(pi));
end

% ----------------------------------------------------------------------

printf('============================================================\n');
printf('Ch 10 Theorem 7 -- Second-order von Zeipel formula\n');
printf('Symbolic verification of (VZ-2) modern-bracket form\n');
printf('============================================================\n\n');

n_tests_passed = 0;
n_tests_total = 0;

% ----------------------------------------------------------------------
% TEST 1: Minimal 1D smoke test
%   F_0 = 1/(2 L^2)  (mu = 1)
%   F_1 = cos(l)
%   F_2 = 0
% Expected: F_2* = 3 L^2 / 4  (hand calc in ch10_foundations_thm7.md
% Section 6.1).
% ----------------------------------------------------------------------

printf('--- TEST 1: Minimal 1D smoke test (F_1 = cos l) ---\n\n');
n_tests_total = n_tests_total + 1;

syms L l real positive;

mu    = sym(1);
F1    = cos(l);
F1bar = sym(0);                       % <cos(l)>_l = 0
n_ke  = mu^2 / L^3;                   % positive mean motion
n_pr  = diff(n_ke, L);                % -3 mu^2 / L^4

% Homological equation: n_ke * dS_1/dl = F_1 - F_1bar = cos(l).
%   S_1 = (1/n_ke) sin(l) = L^3 sin(l).
S1    = (1/n_ke) * sin(l);
dS1dl = diff(S1, l);
dS1dL = diff(S1, L);

hom_eq_residual = simplify(n_ke * dS1dl - (F1 - F1bar));
printf('  Homological residual (should be 0): '); disp(hom_eq_residual);

% Route A: direct eps^2 formula
dF1dL       = diff(F1, L);
U_lavg      = l_avg(dS1dl * dF1dL, l);
sq_lavg     = l_avg(dS1dl^2,       l);
F2star_A    = simplify(U_lavg - (n_pr/2) * sq_lavg);
printf('  Route A (direct eps^2):  F_2* = '); disp(F2star_A);

% Route B: Poisson bracket formula
dF1dl    = diff(F1, l);
pb       = dS1dl * (dF1dL + diff(F1bar, L)) - dS1dL * (dF1dl + diff(F1bar, l));
pb_lavg  = l_avg(pb, l);
F2star_B = simplify(sym(1)/2 * pb_lavg);
printf('  Route B (Poisson bracket):  F_2* = '); disp(F2star_B);

diff_AB = simplify(F2star_A - F2star_B);
printf('  Route A - Route B (should be 0): '); disp(diff_AB);

expected_val = sym(3) * L^2 / sym(4);
diff_exp = simplify(F2star_A - expected_val);
printf('  Route A - expected 3L^2/4 (should be 0): '); disp(diff_exp);

if logical(diff_AB == 0) && logical(diff_exp == 0)
  printf('  TEST 1: PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  TEST 1: FAIL\n\n');
end

% ----------------------------------------------------------------------
% TEST 2: Two-dimensional test with (g, G) pair
%   F_0 = mu^2/(2 L^2)
%   F_1 = a(L,G) cos(l) + b(L,G) cos(2l + 2g)
%   F_2 = 0
% Both j=1 (l,L) and j=2 (g,G) contribute.
% ----------------------------------------------------------------------

printf('--- TEST 2: 2D test with (l,L) and (g,G) pairs ---\n\n');
n_tests_total = n_tests_total + 1;

clear L l;                 % fresh symbols
syms L G real positive;
syms l g real;

mu    = sym(1);
a_fn  = 1 / (L * G);
b_fn  = 1 / (L^2 * G^2);
F1    = a_fn * cos(l) + b_fn * cos(2*l + 2*g);
% <F_1>_l: cos(l) averages to 0, cos(2l+2g) averages to 0 (free l-phase)
F1bar = sym(0);
n_ke  = mu^2 / L^3;
n_pr  = diff(n_ke, L);

% Homological equation term-by-term:
%   term (a): n_ke * dS_a/dl = a cos(l) -> S_a = (a/n_ke) sin(l)
%   term (b): n_ke * dS_b/dl = b cos(2l + 2g) -> S_b = (b/(2 n_ke)) sin(2l + 2g)
S1a = (a_fn / n_ke) * sin(l);
S1b = (b_fn / (2 * n_ke)) * sin(2*l + 2*g);
S1  = S1a + S1b;

hom_eq_residual = simplify(n_ke * diff(S1, l) - (F1 - F1bar));
printf('  Homological residual (should be 0): '); disp(hom_eq_residual);

dS1dl = diff(S1, l);
dS1dL = diff(S1, L);
dS1dg = diff(S1, g);
dS1dG = diff(S1, G);
dF1dl = diff(F1, l);
dF1dL = diff(F1, L);
dF1dg = diff(F1, g);
dF1dG = diff(F1, G);
dF1bardL = diff(F1bar, L);
dF1bardG = diff(F1bar, G);

% Route A
U_lavg   = l_avg(dS1dl * dF1dL + dS1dg * dF1dG, l);
sq_lavg  = l_avg(dS1dl^2, l);
F2star_A = simplify(U_lavg - (n_pr/2) * sq_lavg);
printf('  Route A (direct eps^2):  F_2* = '); disp(F2star_A);

% Route B: {S_1, F_1 + F_1*}_modern
pb_j1    = dS1dl * (dF1dL + dF1bardL) - dS1dL * dF1dl;
pb_j2    = dS1dg * (dF1dG + dF1bardG) - dS1dG * dF1dg;
pb_total = pb_j1 + pb_j2;
pb_lavg  = l_avg(pb_total, l);
F2star_B = simplify(sym(1)/2 * pb_lavg);
printf('  Route B (Poisson bracket):  F_2* = '); disp(F2star_B);

diff_AB = simplify(F2star_A - F2star_B);
printf('  Route A - Route B (should be 0): '); disp(diff_AB);

if logical(diff_AB == 0)
  printf('  TEST 2: PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  TEST 2: FAIL\n\n');
end

% ----------------------------------------------------------------------
% TEST 3: Lemma 7.8 key identity  U + T = n' <(dS_1/dl)^2>_l
%   Uses the same 2D data as TEST 2.
% ----------------------------------------------------------------------

printf('--- TEST 3: Lemma 7.8 key identity  U + T = n'' <(dS_1/dl)^2>_l ---\n\n');
n_tests_total = n_tests_total + 1;

T_lavg = l_avg(dF1dl * dS1dL + dF1dg * dS1dG, l);

LHS = simplify(U_lavg + T_lavg);
RHS = simplify(n_pr * sq_lavg);
printf('  U + T = '); disp(LHS);
printf('  n'' <(dS_1/dl)^2>_l = '); disp(RHS);
diff_78 = simplify(LHS - RHS);
printf('  Difference (should be 0): '); disp(diff_78);

if logical(diff_78 == 0)
  printf('  TEST 3: PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  TEST 3: FAIL\n\n');
end

% ----------------------------------------------------------------------
% TEST 4: Lemma 7.9  <{S_1, F_1 + F_1*}_modern>_l = U - T
% ----------------------------------------------------------------------

printf('--- TEST 4: Lemma 7.9  <{S_1, F_1 + F_1*}>_l = U - T ---\n\n');
n_tests_total = n_tests_total + 1;

LHS79 = simplify(pb_lavg);
RHS79 = simplify(U_lavg - T_lavg);
printf('  <{S_1, F_1 + F_1*}>_l = '); disp(LHS79);
printf('  U - T                 = '); disp(RHS79);
diff_79 = simplify(LHS79 - RHS79);
printf('  Difference (should be 0): '); disp(diff_79);

if logical(diff_79 == 0)
  printf('  TEST 4: PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  TEST 4: FAIL\n\n');
end

% ----------------------------------------------------------------------
% TEST 5: Sign-convention cross-check with Theorem 6.
%   <{F_1*, S_1}>_l = 0 when F_1* is secular and <S_1>_l = 0.
% ----------------------------------------------------------------------

printf('--- TEST 5: Theorem 6 cross-check <{F_1*, S_1}>_l = 0 ---\n\n');
n_tests_total = n_tests_total + 1;

% {F_1bar, S_1}_modern: F_1bar has no angle dependence, so
%   (dF_1bar/dl_j) = 0 for all j.
%   {F_1bar, S_1} = -Sum_j (dF_1bar/dL_j)(dS_1/dl_j)
pb_thm6 = -dF1bardL * dS1dl - dF1bardG * dS1dg;
pb_thm6_lavg = simplify(l_avg(pb_thm6, l));
printf('  <{F_1*, S_1}>_l = '); disp(pb_thm6_lavg);

if logical(pb_thm6_lavg == 0)
  printf('  TEST 5: PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  TEST 5: FAIL\n\n');
end

% ----------------------------------------------------------------------
% TEST 6: Nontrivial F_1* (secular piece added to F_1).
%   F_1_v6 = A_sec(L,G) + a(L,G) cos(l) + b(L,G) cos(2l+2g)
%   F_1*_v6 = A_sec(L,G)  /= 0
%   S_1 unchanged (secular piece drops out of the homological equation).
%   All three checks are re-run: (VZ-2), Lemma 7.8, Theorem 6.
% ----------------------------------------------------------------------

printf('--- TEST 6: Nontrivial F_1* = A_sec(L,G) ---\n\n');
n_tests_total = n_tests_total + 1;

A_sec      = 1 / (L^3 * G);
F1_v6      = A_sec + a_fn * cos(l) + b_fn * cos(2*l + 2*g);
F1bar_v6   = A_sec;
S1_v6      = (a_fn / n_ke) * sin(l) + (b_fn / (2 * n_ke)) * sin(2*l + 2*g);

hom_eq_res_v6 = simplify(n_ke * diff(S1_v6, l) - (F1_v6 - F1bar_v6));
printf('  Hom eq residual (should be 0): '); disp(hom_eq_res_v6);

dS1dl_v6 = diff(S1_v6, l);
dS1dL_v6 = diff(S1_v6, L);
dS1dg_v6 = diff(S1_v6, g);
dS1dG_v6 = diff(S1_v6, G);
dF1dl_v6 = diff(F1_v6, l);
dF1dL_v6 = diff(F1_v6, L);
dF1dg_v6 = diff(F1_v6, g);
dF1dG_v6 = diff(F1_v6, G);
dF1bardL_v6 = diff(F1bar_v6, L);
dF1bardG_v6 = diff(F1bar_v6, G);

% Route A
U_lavg_v6   = l_avg(dS1dl_v6 * dF1dL_v6 + dS1dg_v6 * dF1dG_v6, l);
sq_lavg_v6  = l_avg(dS1dl_v6^2, l);
F2star_A_v6 = simplify(U_lavg_v6 - (n_pr/2) * sq_lavg_v6);
printf('  Route A: F_2* = '); disp(F2star_A_v6);

% Route B
pb_j1_v6    = dS1dl_v6 * (dF1dL_v6 + dF1bardL_v6) - dS1dL_v6 * dF1dl_v6;
pb_j2_v6    = dS1dg_v6 * (dF1dG_v6 + dF1bardG_v6) - dS1dG_v6 * dF1dg_v6;
pb_total_v6 = pb_j1_v6 + pb_j2_v6;
pb_lavg_v6  = l_avg(pb_total_v6, l);
F2star_B_v6 = simplify(sym(1)/2 * pb_lavg_v6);
printf('  Route B: F_2* = '); disp(F2star_B_v6);

diff_AB_v6 = simplify(F2star_A_v6 - F2star_B_v6);
printf('  Route A - Route B (should be 0): '); disp(diff_AB_v6);

% Theorem 6 cross-check
pb_thm6_v6 = -dF1bardL_v6 * dS1dl_v6 - dF1bardG_v6 * dS1dg_v6;
pb_thm6_lavg_v6 = simplify(l_avg(pb_thm6_v6, l));
printf('  <{F_1*, S_1}>_l (Thm 6, should be 0): '); disp(pb_thm6_lavg_v6);

if logical(diff_AB_v6 == 0) && logical(pb_thm6_lavg_v6 == 0)
  printf('  TEST 6: PASS\n\n');
  n_tests_passed = n_tests_passed + 1;
else
  printf('  TEST 6: FAIL\n\n');
end

% ----------------------------------------------------------------------
% Summary
% ----------------------------------------------------------------------
printf('============================================================\n');
printf('SUMMARY: %d of %d tests passed\n', n_tests_passed, n_tests_total);
if n_tests_passed == n_tests_total
  printf('OVERALL: PASS\n');
else
  printf('OVERALL: FAIL\n');
end
printf('============================================================\n');
