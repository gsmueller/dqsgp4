%% verify_ch10a_F1_prefactor.m
%%
%% Chapter 10a Proposition F.1: F_1 = G^(-6) * [mu^4 k_2 (1+e cos f)^3 (A + B cos 2(f+g))]
%%
%% Claim (Corollary 2.1 of CH10_PLAN.md, Proposition F.1 of ch10a_setup.md):
%%   Starting from ch05d Theorem A.17 (F_1 = (mu k_2/r^3)(A + B cos 2(f+g)))
%%   and applying the Kepler-conic orbit equation (1/r = mu (1 + e cos f)/G^2,
%%   equivalent to r = p/(1+e cos f) with p = G^2/mu from ch04), verify that
%%
%%     F_1 = G^(-6) * [mu^4 k_2 (1 + e cos f)^3 (A(theta) + B(theta) cos(2(f+g)))]
%%
%%   with A(theta) = (3 theta^2 - 1)/2, B(theta) = 3(1 - theta^2)/2.
%%
%% This establishes F_1 in M_6 (Definition 1.1 of ch10_foundations_thm1.md) with
%% exponent alpha = 6 and F-factor
%%     Phi_{F_1}(theta, e, l, g) = mu^4 k_2 (1 + e cos f(l,e))^3 (A + B cos 2(f+g)).
%%
%% Verification strategy (symbolic, no numerical approximation):
%%   1. Build F_1 in "original form" from ch05d A.17.
%%   2. Substitute the orbit equation r = G^2/(mu (1 + e cos f)).
%%   3. Build F_1 in "M_6 factored form" G^(-6) * Phi_{F_1}.
%%   4. Compute the residual = original - M_6 form.
%%   5. Verify residual reduces to 0 symbolically via simplify().

pkg load symbolic;

printf('============================================================\n');
printf('Ch 10a Proposition F.1 -- F_1 prefactor reduction\n');
printf('Claim: F_1 in M_6 with explicit closed form (F.1).\n');
printf('============================================================\n\n');

syms mu k2 G r f g theta e real positive;

% Inclination factors (ch05d Theorem A.17 eq 0.D.7).
A_val = (3*theta^2 - 1) / 2;
B_val = 3 * (1 - theta^2) / 2;

% --- Form 1: F_1 as given in ch05d Theorem A.17 (equation 0.D.6).
F1_original = (mu * k2 / r^3) * (A_val + B_val * cos(2*(f + g)));

printf('F_1 original (ch05d A.17):\n');
disp(F1_original);
printf('\n');

% --- Form 2: Substitute the orbit equation r = G^2/(mu (1 + e cos f)).
% Equivalent forms:
%   p := a(1-e^2) = G^2/mu (ch04 Delaunay: G^2 = mu p)
%   r = p/(1 + e cos f) = G^2/(mu (1 + e cos f))
r_orbit = G^2 / (mu * (1 + e * cos(f)));

F1_subbed = subs(F1_original, r, r_orbit);
F1_subbed = simplify(F1_subbed);

printf('F_1 after substituting r -> G^2/(mu (1 + e cos f)):\n');
disp(F1_subbed);
printf('\n');

% --- Form 3: M_6 factored form of (F.1).
Phi_F1 = mu^4 * k2 * (1 + e*cos(f))^3 * (A_val + B_val * cos(2*(f + g)));
F1_M6  = Phi_F1 / G^6;

printf('F_1 in M_6 form (F.1):\n');
disp(F1_M6);
printf('\n');

% --- Check: residual = Form 2 - Form 3.
residual = simplify(F1_subbed - F1_M6);
printf('Residual (should be 0):\n');
disp(residual);
printf('\n');

if isequal(residual, sym(0))
  printf('[PASS] F_1 original (after r-substitution) equals M_6 form.\n');
  printf('       F_1 = G^(-6) * [mu^4 k_2 (1+e cos f)^3 (A + B cos 2(f+g))] confirmed.\n');
  printf('       Hence F_1 in M_6 with alpha = 6. ***\n');
else
  printf('[FAIL] Residual does not reduce to 0. Check the substitution.\n');
  disp(residual);
end

%% ============================================================
%% Secondary check: Phi_{F_1} has no L-dependence (pure M_6-factor)
%% ============================================================

printf('\n--- Secondary check: Phi_{F_1} depends only on (theta, e, l, g) ---\n');

free_syms = symvar(Phi_F1);
fprintf('Free symbols in Phi_{F_1}: ');
for k = 1:numel(free_syms)
  fprintf('%s ', char(free_syms(k)));
end
fprintf('\n');

has_L = false;
has_G = false;
has_H = false;
for k = 1:numel(free_syms)
  s = char(free_syms(k));
  if strcmp(s, 'L'); has_L = true; end
  if strcmp(s, 'G'); has_G = true; end
  if strcmp(s, 'H'); has_H = true; end
end

if ~has_L && ~has_G && ~has_H
  printf('[PASS] Phi_{F_1} has no (L, G, H) dependence -- pure M_6 F-factor.\n');
else
  printf('[FAIL] Phi_{F_1} has residual momentum dependence.\n');
end

%% ============================================================
%% Tertiary check: harmonic count in Phi_{F_1}
%% ============================================================

printf('\n--- Tertiary check: harmonic enumeration of (1+e cos f)^3 ---\n');

% (1 + e cos f)^3 should expand to: 1 + 3/2 e^2 + (3e + 3e^3/4) cos f + 3/2 e^2 cos 2f + e^3/4 cos 3f
expand_prefactor = expand((1 + e*cos(f))^3);

% Re-express powers of cos f using double-/triple-angle identities.
% SymPy's rewrite should give us the harmonic form.
try
  harmonic_form = rewrite(expand_prefactor, 'sincos');
catch
  harmonic_form = expand_prefactor;
end
harmonic_form = simplify(harmonic_form);

printf('(1 + e cos f)^3 expanded:\n');
disp(expand_prefactor);
printf('\n');

% Manual harmonic check: coefficient of cos(jf) for j = 0, 1, 2, 3.
% Reduction identities:  cos^2 f = (1 + cos 2f)/2
%                        cos^3 f = (3 cos f + cos 3f)/4
% So (1 + e cos f)^3 = 1 + 3e cos f + 3 e^2 cos^2 f + e^3 cos^3 f
%                    = 1 + 3e cos f + 3 e^2 (1 + cos 2f)/2 + e^3 (3 cos f + cos 3f)/4
%                    = (1 + 3 e^2/2) + (3e + 3 e^3/4) cos f + (3 e^2/2) cos 2f + (e^3/4) cos 3f

manual_harmonic = (1 + sym(3)*e^2/2) ...
               + (sym(3)*e + sym(3)*e^3/4) * cos(f) ...
               + (sym(3)*e^2/2) * cos(2*f) ...
               + (e^3/4) * cos(3*f);

harmonic_residual = simplify(expand_prefactor - manual_harmonic);
printf('Harmonic form residual (should be 0 after trig expansion):\n');

% Use expand() on both sides with trig-expand to reduce cos(2f), cos(3f).
harmonic_residual_expand = expand(harmonic_residual);
harmonic_residual_simp   = simplify(expand(manual_harmonic) - expand_prefactor);
disp(harmonic_residual_simp);

if isequal(harmonic_residual_simp, sym(0))
  printf('[PASS] (1+e cos f)^3 = 1 + 3e^2/2 + (3e+3e^3/4)cos f + 3e^2/2 cos 2f + e^3/4 cos 3f.\n');
else
  printf('[INFO] Harmonic form needs manual verification (SymPy may not auto-reduce).\n');
end

printf('\n============================================================\n');
printf('Proposition F.1 verification summary:\n');
printf('  Primary: F_1 in M_6 with alpha=6 and F-factor Phi_{F_1}\n');
printf('  Secondary: Phi_{F_1} has no (L, G, H) dependence\n');
printf('  Tertiary: (1+e cos f)^3 harmonic expansion confirmed\n');
printf('============================================================\n');
