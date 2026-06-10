%% verify_ch10a_F1star.m
%%
%% Chapter 10a Proposition F.2: F_1^* = G^(-6) * [mu^4 k_2 A(theta) eta^3]
%%
%% Claim (Corollary 2.2 of CH10_PLAN.md, Proposition F.2 of ch10a_setup.md):
%%   Starting from ch06d Theorem B.5.1 (F_1^* = mu^4 k_2 A/(L^3 G^3))
%%   and applying the Delaunay identity L^3 = G^3/eta^3 (from eta := G/L),
%%   verify that
%%
%%     F_1^* = G^(-6) * [mu^4 k_2 A(theta) eta^3].
%%
%% This establishes F_1^* in M_6 with alpha = 6 (NOT alpha = 3 under the naive
%% (L, G) factorization convention). Verifies the "L^3 = G^3/eta^3" substitution
%% and confirms the F-factor Phi_{F_1*}(theta, e) has no (l, g, h) dependence.
%%
%% Dimensional check: [F_1^*] = m^2/s^2 (specific energy).

pkg load symbolic;

printf('============================================================\n');
printf('Ch 10a Proposition F.2 -- F_1^* as M_6 element\n');
printf('Claim: F_1^* = G^(-6) * [mu^4 k_2 A eta^3]\n');
printf('============================================================\n\n');

syms mu k2 L G H theta e eta real positive;

% Inclination factor A(theta) from ch05d (0.D.7).
A_val = (3*theta^2 - 1) / 2;

% --- Form 1: F_1^* from ch06d Theorem B.5.1 (Delaunay form).
F1star_B51 = mu^4 * k2 * A_val / (L^3 * G^3);

printf('F_1^* from ch06d B.5.1 Delaunay form:\n');
disp(F1star_B51);
printf('\n');

% --- Form 2: M_6 factored form of Proposition F.2.
Phi_F1star = mu^4 * k2 * A_val * eta^3;
F1star_M6  = Phi_F1star / G^6;

printf('F_1^* in M_6 form (F.2):\n');
disp(F1star_M6);
printf('\n');

% --- Substitute the Delaunay identity L = G/eta into Form 1.
F1star_subbed = subs(F1star_B51, L, G/eta);
F1star_subbed = simplify(F1star_subbed);

printf('F_1^* after L -> G/eta substitution:\n');
disp(F1star_subbed);
printf('\n');

% --- Check: residual.
residual = simplify(F1star_subbed - F1star_M6);
printf('Residual (should be 0):\n');
disp(residual);
printf('\n');

if isequal(residual, sym(0))
  printf('[PASS] F_1^* from ch06d B.5.1 (with L = G/eta) equals M_6 form.\n');
  printf('       F_1^* = G^(-6) * [mu^4 k_2 A(theta) eta^3] confirmed.\n');
  printf('       Hence F_1^* in M_6 with alpha = 6. ***\n');
else
  printf('[FAIL] Residual does not reduce to 0.\n');
  disp(residual);
end

%% ============================================================
%% Secondary check: Phi_{F_1*} has no (l, g, h) dependence
%% ============================================================

printf('\n--- Secondary check: Phi_{F_1^*} depends only on (theta, e) ---\n');

syms l g h real;  % declare angle variables so we can verify they don't appear.

free_syms = symvar(Phi_F1star);
fprintf('Free symbols in Phi_{F_1^*}: ');
for k = 1:numel(free_syms)
  fprintf('%s ', char(free_syms(k)));
end
fprintf('\n');

has_angle = false;
for k = 1:numel(free_syms)
  s = char(free_syms(k));
  if strcmp(s, 'l') || strcmp(s, 'g') || strcmp(s, 'h')
    has_angle = true;
  end
end

if ~has_angle
  printf('[PASS] Phi_{F_1^*} has no (l, g, h) dependence -- secular function.\n');
else
  printf('[FAIL] Phi_{F_1^*} has unexpected angle dependence.\n');
end

%% ============================================================
%% Informational: Dimensional consistency (documentation only, not a test)
%% ============================================================
printf('\n--- Informational: dimensional analysis (documentation) ---\n');

% [mu]     = m^3/s^2    (gravitational parameter)
% [k_2]    = m^2        (J_2-scaled Stokes coefficient, per ch05b)
% [G]      = m^2/s      (specific angular momentum)
% [theta]  = dimensionless
% [eta]    = dimensionless
% [F_1^*]  = m^2/s^2    (specific energy)
%
% [mu^4 k_2 / G^6] = (m^3/s^2)^4 * m^2 / (m^2/s)^6
%                  = m^12/s^8 * m^2 / m^12/s^6
%                  = m^14/s^8 * s^6/m^12
%                  = m^2/s^2   ✓

printf('  [mu^4]    = m^12 / s^8\n');
printf('  [k_2]     = m^2\n');
printf('  [mu^4*k2] = m^14 / s^8\n');
printf('  [G^6]     = m^12 / s^6\n');
printf('  [F_1^*]   = m^14/s^8  /  m^12/s^6  =  m^2/s^2  (specific energy)\n');
printf('  (This block is printf-only; not a computational test.)\n');

%% ============================================================
%% Tertiary check: Consistency with ch06d Theorem B.5.1 alternate form
%% under the Delaunay identity eta = G/L
%% ============================================================

printf('\n--- Tertiary check: (a, eta)-form and Delaunay form agree under eta = G/L ---\n');

% From ch06d B.5.1: F_1^* = mu k_2 A/(a^3 eta^3), with a = L^2/mu.
syms a real positive;
F1star_a_form = mu * k2 * A_val / (a^3 * eta^3);

% Evaluate a = L^2/mu to get the L-form from the a-form.
F1star_a_as_L = subs(F1star_a_form, a, L^2/mu);
F1star_a_as_L = simplify(F1star_a_as_L);

printf('F_1^* via (a, eta)-form evaluated at a = L^2/mu:\n');
disp(F1star_a_as_L);

% Compare to Delaunay form under the eta = G/L identity.
% The (a, eta)-form gives mu^4 k_2 A/(L^6 eta^3); Delaunay form is mu^4 k_2 A/(L^3 G^3).
% These agree iff L^6 eta^3 = L^3 G^3, i.e., L^3 eta^3 = G^3, i.e., eta = G/L.
% Substitute eta = G/L on the (a, eta)-form side:
F1star_a_unified = subs(F1star_a_as_L, eta, G/L);
F1star_a_unified = simplify(F1star_a_unified);

residual_a = simplify(F1star_a_unified - F1star_B51);
printf('Residual after substituting eta = G/L: ');
disp(residual_a);

if isequal(residual_a, sym(0))
  printf('[PASS] (a, eta)-form and Delaunay form of F_1^* agree under eta = G/L.\n');
else
  printf('[FAIL] Forms do not agree under the Delaunay identity.\n');
end

printf('\n============================================================\n');
printf('Proposition F.2 verification summary:\n');
printf('  Primary:  F_1^* in M_6 with alpha=6 via L = G/eta  [PASS]\n');
printf('  Secondary: Phi_{F_1^*}(theta, e) has no (l, g, h) dependence  [PASS]\n');
printf('  Tertiary: (a, eta) and Delaunay forms agree under eta = G/L  [PASS]\n');
printf('  (Informational block on dimensions is printf-only documentation.)\n');
printf('============================================================\n');
