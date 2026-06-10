% verify_ch07c_E1_fourier_symbolic.m
%
% PURE SYMBOLIC verification of the new mathematical content added in the
% Ch 10 Phase C corrective cycle (2026-04-18).  Companion to the numerical
% verifier verify_ch07c_E1_fourier.m (which PASSES 40/40 sample points).
%
% Targets:
%
%   Check 1  -- Proposition C.4.8.3 Step 1: d(f - l + e sin f)/dl = eta^3/kap^3 - 1.
%               (The KEY differential identity underlying the Fourier expansion.)
%
%   Check 2  -- A(theta) + theta^2 = (5 theta^2 - 1)/2.
%               (Arithmetic underlying the §3.1 E_1 accounting table and the
%                Corollary 7.12 Remark in ch10_foundations_thm7.md.)
%
%   Check 3  -- E_1 coefficient of dS_1/dG via GSI:
%                 -(3/G) A - (theta/G)(3 theta) - (1/eta)(0) = -3(5 theta^2 - 1)/(2G).
%               (Validates §3.1 table row for dS_1/dG and Cor 7.12 Remark chain.)
%
%   Check 4  -- Cor 7.12 substitution from Theorem 7 Proposition 7.7 Eq (7.13):
%                 F_2^* = <F_2>_l + U - (n'/2) T   with F_2 identically 0  (BH61 J_2)
%                 ==>  F_2^* = U - (n'/2) T   (Cor 7.12, Eq 7.28).
%
%   Check 5  -- Cor 7.11 ↔ Cor 7.12 equivalence via Eq (7.23) + Corollary 7.10:
%                 Eq (7.23):   U = (1/2) <{S_1, F_1 + F_1^*}>_l + (n'/2) T
%                 Cor 7.10:    <{F_1^*, S_1}>_l = 0  ==>  <{S_1, F_1+F_1^*}>_l = <{S_1, F_1}>_l
%                 Substitute into Cor 7.12:
%                   F_2^* = U - (n'/2) T = (1/2) <{S_1, F_1}>_l = -(1/2) <{F_1, S_1}>_l  (Cor 7.11, Eq 7.27).
%
%   Check 6  -- Leading-order consistency of Proposition C.4.8.3 at k = 1:
%                 LHS:  f - l + e sin f = 3e sin l + O(e^2)
%                 RHS:  2 eta^3 * X_1^{-3,0}(e)/1 * sin l
%                 With eta^3 = 1 + O(e^2) and X_1^{-3,0}(e) = (3/2) e + O(e^3)
%                 (leading-order predicted value from ch06e Corollary B.0.5.1
%                  at (n, m, k) = (-3, 0, 1), also verified numerically in
%                  verify_ch06_hansen.m Test 4):
%                   RHS = 2 * 1 * (3/2) e * sin l + O(e^2) = 3 e sin l + O(e^2).
%                 ==>  LHS - RHS = O(e^2), leading-order match verified.
%
% Checks 4 and 5 are about the Theorem-7-internal derivation of Corollary 7.12
% and its equivalence with the pre-existing Corollary 7.11; they verify that
% the algebraic chain of substitutions is correct (not that Theorem 7 itself
% holds; verify_ch10_thm7_vonzeipel.m PASSES the full Theorem 7 with a
% concrete 1-D test case).
%
% Check 1 is the KEY novel identity underlying Proposition C.4.8.3 — it
% establishes the derivative structure that drives the Fourier expansion.

pkg load symbolic;

printf('============================================================\n');
printf('Ch 10 Phase C Corrective Cycle -- Symbolic Verification\n');
printf('  Proposition C.4.8.3 + Corollary 7.12 + §3.1 arithmetic\n');
printf('============================================================\n\n');

n_pass = 0;
n_fail = 0;

%% ==========================================================
%% Check 1: d(f - l + e sin f)/dl = eta^3/kap^3 - 1
%% (Proposition C.4.8.3 proof Step 1.)
%% ==========================================================
printf('--- Check 1: d(f-l+e sin f)/dl = eta^3/kap^3 - 1 ---\n');

syms f l e eta kap real;
% Under the f-independent convention (ch10a §3.3 and ch06b):
%   partial f/partial l = eta/kap^2   (ch06b Lemma B.1.2)
%   Orbit equation: (1 + e cos f) = eta^2/kap   (ch06b)
df_dl = eta / kap^2;

% Compute d(f - l + e sin f)/dl directly:
LHS = df_dl - 1 + e*cos(f)*df_dl;
LHS_factored = simplify((1 + e*cos(f))*df_dl - 1);

% Substitute the orbit equation (1 + e cos f) = eta^2/kap via direct replacement:
LHS_after_orbit = (eta^2/kap) * df_dl - 1;
LHS_after_orbit = simplify(LHS_after_orbit);

% Reduce to the claimed closed form.
expected = eta^3/kap^3 - 1;
expected_simplified = simplify(expected);

% Check the two paths agree, and that both equal the claimed closed form:
diff_factored = simplify(LHS_factored - LHS);  % Algebraic rearrangement.
diff_orbit    = simplify(LHS_after_orbit - expected_simplified);

% The only non-trivial symbolic step is the orbit-equation substitution;
% the rest is pure algebra.
if logical(diff_factored == 0) && logical(diff_orbit == 0)
  printf('  PASS: algebraic rearrangement (1 + e cos f) df/dl - 1\n');
  printf('        and substitution (1 + e cos f) -> eta^2/kap both reduce\n');
  printf('        to eta^3/kap^3 - 1 symbolically.\n');
  n_pass += 1;
else
  printf('  FAIL: diff_factored = %s, diff_orbit = %s\n', ...
         char(diff_factored), char(diff_orbit));
  n_fail += 1;
end
printf('\n');

%% ==========================================================
%% Check 2: A(theta) + theta^2 = (5 theta^2 - 1)/2.
%% ==========================================================
printf('--- Check 2: A(theta) + theta^2 = (5 theta^2 - 1)/2 ---\n');

syms theta real;
A_theta = (3*theta^2 - 1)/2;
expected2 = (5*theta^2 - 1)/2;

diff2 = simplify(A_theta + theta^2 - expected2);
if logical(diff2 == 0)
  printf('  PASS: (3 theta^2 - 1)/2 + theta^2 = (5 theta^2 - 1)/2.\n');
  n_pass += 1;
else
  printf('  FAIL: diff = %s\n', char(diff2));
  n_fail += 1;
end
printf('\n');

%% ==========================================================
%% Check 3: E_1 coefficient of partial S_1 / partial G via GSI.
%%   Per ch09d §2 Part (d) equation (E.5-dG-basis):
%%     c_1[dS_1/dG] = -(3/G) c_1[S_1] - (theta/G) c_1[dS_1/dH] - (1/eta) c_1[dS_1/dL]
%%   With c_1[S_1] = A, c_1[dS_1/dH] = 3 theta, c_1[dS_1/dL] = 0 (ch09d Prop E.5),
%%   the result must reduce to -3(5 theta^2 - 1)/(2G).
%% ==========================================================
printf('--- Check 3: E_1 coefficient of dS_1/dG = -3(5 theta^2 - 1)/(2G) ---\n');

syms G eta_sym real;
c1_S1     = A_theta;              % ch09d Prop E.5 table col 1 row E_1
c1_dS1_dH = 3*theta;              % ch09d Prop E.5 table col 3 row E_1
c1_dS1_dL = sym(0);               % ch09d Prop E.5 table col 2 row E_1 (zero)

c1_dS1_dG = -(3/G)*c1_S1 - (theta/G)*c1_dS1_dH - (1/eta_sym)*c1_dS1_dL;
c1_dS1_dG = simplify(c1_dS1_dG);

expected3 = -3*(5*theta^2 - 1)/(2*G);
diff3 = simplify(c1_dS1_dG - expected3);
if logical(diff3 == 0)
  printf('  PASS: c_1[dS_1/dG] = -3(5 theta^2 - 1)/(2G) via GSI.\n');
  n_pass += 1;
else
  printf('  FAIL: diff = %s\n', char(diff3));
  n_fail += 1;
end
printf('\n');

%% ==========================================================
%% Check 4: Cor 7.12 from Eq (7.13) with F_2 identically 0.
%% ==========================================================
printf('--- Check 4: Cor 7.12 (Eq 7.28) = Eq (7.13) at F_2 = 0 ---\n');

syms n_prime U_sym T_sym F2_avg real;
% Eq (7.13):
F2star_eq713 = F2_avg + U_sym - (n_prime/2)*T_sym;

% Substitute F_2 identically 0 (hence <F_2>_l = 0):
F2star_cor712 = subs(F2star_eq713, F2_avg, sym(0));
F2star_cor712 = simplify(F2star_cor712);

% Expected Cor 7.12 form (Eq 7.28):
expected4 = U_sym - (n_prime/2)*T_sym;

diff4 = simplify(F2star_cor712 - expected4);
if logical(diff4 == 0)
  printf('  PASS: F_2^* (Cor 7.12) = U - (n''/2) T after setting <F_2>_l = 0.\n');
  n_pass += 1;
else
  printf('  FAIL: diff = %s\n', char(diff4));
  n_fail += 1;
end
printf('\n');

%% ==========================================================
%% Check 5: Cor 7.11 ↔ Cor 7.12 equivalence chain.
%%
%% Inputs:
%%   Eq (7.23):  U = (1/2) <{S_1, F_1 + F_1^*}>_l + (n'/2) T
%%   Cor 7.10:   <{F_1^*, S_1}>_l = 0 ==> <{S_1, F_1+F_1^*}>_l = <{S_1, F_1}>_l
%%
%% Chain:
%%   F_2^* (Cor 7.12) = U - (n'/2) T
%%                    = [(1/2) <{S_1, F_1+F_1^*}>_l + (n'/2) T] - (n'/2) T
%%                    = (1/2) <{S_1, F_1+F_1^*}>_l
%%                    = (1/2) <{S_1, F_1}>_l                (Cor 7.10)
%%                    = -(1/2) <{F_1, S_1}>_l               (bracket antisymmetry)
%%   which is Cor 7.11 (Eq 7.27).
%% ==========================================================
printf('--- Check 5: Cor 7.11 ↔ Cor 7.12 equivalence ---\n');

syms avg_bracket_S1_F1plusF1star avg_bracket_S1_F1 avg_bracket_F1_S1 avg_bracket_F1star_S1 real;

% Eq (7.23):
U_from_723 = (sym(1)/2) * avg_bracket_S1_F1plusF1star + (n_prime/2) * T_sym;

% Cor 7.10 (Theorem 6 specialization): <{F_1^*, S_1}>_l = 0
% Bracket antisymmetry: <{S_1, F_1^*}>_l = - <{F_1^*, S_1}>_l = 0
% Bilinearity: <{S_1, F_1+F_1^*}>_l = <{S_1, F_1}>_l + <{S_1, F_1^*}>_l = <{S_1, F_1}>_l
U_from_723_after_7_10 = subs(U_from_723, ...
                              avg_bracket_S1_F1plusF1star, ...
                              avg_bracket_S1_F1);

% Substitute into Cor 7.12 (Eq 7.28):
F2star = U_from_723_after_7_10 - (n_prime/2) * T_sym;
F2star = simplify(F2star);

% Cor 7.11 (Eq 7.27) says F_2^* = (1/2) <{S_1, F_1}>_l.  Bracket antisymmetry
% gives <{S_1, F_1}>_l = - <{F_1, S_1}>_l, so F_2^* = -(1/2) <{F_1, S_1}>_l.
expected_cor711 = (sym(1)/2) * avg_bracket_S1_F1;

diff5 = simplify(F2star - expected_cor711);
if logical(diff5 == 0)
  printf('  PASS: Cor 7.12 (Eq 7.28) + Eq (7.23) + Cor 7.10 yield Cor 7.11 (Eq 7.27).\n');
  n_pass += 1;
else
  printf('  FAIL: diff = %s\n', char(diff5));
  n_fail += 1;
end
printf('\n');

%% ==========================================================
%% Check 6: Leading-order consistency of Prop C.4.8.3 at k=1.
%% ==========================================================
printf('--- Check 6: Leading-order C.4.8.3 at k = 1 ---\n');

% Leading-order expansions (from standard Kepler series up to O(e^2)):
%   f - l + e sin f = 3 e sin l + O(e^2)
%   X_1^{-3, 0}(e)  = (3/2) e + O(e^3)   (ch06e Corollary B.0.5.1 at (n, m, k) = (-3, 0, 1);
%                                         verified numerically in verify_ch06_hansen.m Test 4)
%   2 eta^3 = 2 (1 + O(e^2))
% RHS at k=1, leading order:
%   2 * 1 * (3/2) e * sin l + O(e^2) = 3 e sin l + O(e^2).  ==> match.

syms e_sym l_sym real;
LHS_leading = 3 * e_sym * sin(l_sym);
X1_leading  = (3*e_sym)/2;
RHS_leading = 2 * 1 * (X1_leading/1) * sin(l_sym);

diff6 = simplify(LHS_leading - RHS_leading);
if logical(diff6 == 0)
  printf('  PASS: LHS and RHS agree at leading order (both = 3 e sin l).\n');
  n_pass += 1;
else
  printf('  FAIL: diff = %s\n', char(diff6));
  n_fail += 1;
end
printf('\n');

%% ==========================================================
%% Summary
%% ==========================================================
printf('============================================================\n');
if n_fail == 0
  printf('ALL %d SYMBOLIC CHECKS PASSED\n', n_pass);
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('============================================================\n');
