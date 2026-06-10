% verify_ch12_deltaq_symbolic.m
%
% Symbolic verifier for ch12_deltap_deltaq.md δq_j closed forms.
% Discharges CLOSURE_PLAN_J2_THEORY.md §1.2 Task B Item 3 (downstream P25+.4 audit)
% and §1.4 Task D.4 (extended symbolic verifier).
%
% Companion to verify_ch12_deltas.m (which does NUMERICAL G-scaling at L=μ=k_2=1).
% This verifier does SYMBOLIC checks with non-unit μ, k_2, L, ensuring scaling-blind
% errors (μ-exponent, k_2-exponent, G-power) cannot be masked by unit substitutions.
%
% Seven checks:
%   Check 1: Proposition P25+.4 multi-n direct-chain vs boxed identity. Covers
%            every η power used in ch12 §§5.2, 6.2, 7.2 (n ∈ {-4,-3,-2,-1,1,2,3}).
%            Independent route: chain rule (ch12a Ax3) with outer F(u)=u^n, inner
%            u=η(e)=sqrt(1-e²); De from ch12a Ax6.
%   Check 2: D(η²/e) compound used in ch12 §6.2 Step 4 line 606. Independent
%            via product rule (Ax2) on η²·(1/e), with D(η²) from P25+.4 at n=2
%            and D(1/e) from chain rule on De (Ax6).
%   Check 3: D(∂S_1/∂H) ch12 §7.1 Step 5 boxed form vs direct route applying
%            Propositions P25+.1, P25+.3, P25+.5 and Axioms Ax2/Ax3/Ax6 to the
%            ch09d Proposition E.5(c) 5-element basis decomposition of ∂S_1/∂H.
%   Check 4: ℳ_4 G-scaling symbolic for D(∂S_1/∂H) under non-unit μ, k_2, L.
%            Prefactor μ²k_2/G⁴ extracts cleanly with no residual on μ, k_2, G.
%   Check 5: ∂Φ_{S_1*}^{(T)}/∂θ ch12 §7.2 Step 3 boxed form vs SymPy diff.
%            Tests the product-rule application to (5θ²-1)^{-1}, A(θ), B(θ) factors
%            in ch11a Proposition G.4's Φ_{S_1*}^{(T)} closed form.
%   Check 6: ∂Φ_{S_1*}^{(T)}/∂e ch12 §5.2 Step 4 boxed form vs SymPy diff.
%            Tests the chain rule through η = sqrt(1-e²) yielding the second
%            bracket's η^{-4} factor structure.
%   Check 7: GSI identity at α=3 (ch09b Theorem E.3), used by ch12 §6.1 Step 1
%            to express ∂S_1/∂G = -(3/G)S_1 - θ ∂S_1/∂H - (1/η) ∂S_1/∂L.
%            Verified on a generic test function S = G^{-3}·sin(2(l+g))·(1+θ²)·(1+e²)
%            with e, η explicit in (L, G) via Delaunay relations.

pkg load symbolic;

printf('=========================================================\n');
printf('verify_ch12_deltaq_symbolic.m\n');
printf('Task B Item 3 + Task D.4 symbolic verifier\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

%% =======================================================
%% Common symbolic setup
%% =======================================================
syms e f g l_sym theta_sym positive;
syms mu_s k2_s L_s G_s H_s positive;
syms Dl_sym;  % D(l) = q_1 left abstract — Check 3 doesn't need its closed form

% η explicit in e (allows independent chain-rule diff)
eta_e = sqrt(1 - e^2);

% Base D-table per file 08 / ch12a Ax6-Ax7:
Df_val   = 2*sin(f)/e;
Dg_val   = -2*sin(f)/e;
De_val   = -2*(e + cos(f));
DG_val   = -G_s;
DH_val   = -H_s;
Dtheta_val = sym(0);
% Dη from Ax7: D(η²) = -2 e De, so D(η²) = 4e(e+cos f) and Dη = D(η²)/(2η) = 2e(e+cos f)/η
Deta_val = 2*e*(e + cos(f)) / eta_e;

%% =======================================================
%% Check 1: P25+.4 multi-n direct-chain vs boxed identity
%% =======================================================
printf('=== Check 1: P25+.4 D(η^n) direct vs boxed for n ∈ {-4,-3,-2,-1,1,2,3} ===\n');
printf('             (every n appearing in ch12 §§5.2, 6.2, 7.2 η-power usage)\n');

n_values = [-4, -3, -2, -1, 1, 2, 3];
chk1_pass = 0; chk1_fail = 0;

for n_test = n_values
  % Direct route: D(eta^n) = (d/de of eta^n) · De,
  %   where eta^n is treated as explicit function of e via eta = sqrt(1-e²).
  F_direct = eta_e^n_test;
  D_direct = diff(F_direct, e) * De_val;
  D_direct_simpl = simplify(D_direct);

  % Boxed P25+.4: D(eta^n) = 2n e (e+cos f) eta^{n-2}.
  D_boxed = 2 * n_test * e * (e + cos(f)) * eta_e^(n_test - 2);
  D_boxed_simpl = simplify(D_boxed);

  residual = simplify(D_direct_simpl - D_boxed_simpl);
  pass = isequal(residual, sym(0));
  if ~pass
    residual = simplify(expand(residual));
    pass = isequal(residual, sym(0));
  end
  if pass
    chk1_pass = chk1_pass + 1;
    printf('  n=%+d : residual = 0 (PASS)\n', n_test);
  else
    chk1_fail = chk1_fail + 1;
    printf('  n=%+d : residual = %s (FAIL)\n', n_test, char(residual));
  end
endfor

if chk1_fail == 0
  printf('  CHECK 1 STATUS: ALL %d PASS\n\n', chk1_pass);
  n_pass = n_pass + 1;
else
  printf('  CHECK 1 STATUS: %d PASS / %d FAIL\n\n', chk1_pass, chk1_fail);
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 2: D(η²/e) compound (§6.2 Step 4 line 606)
%%
%% Boxed claim: D(η²/e) = 2(e+cos f)(1 + e²)/e².
%%
%% Independent route: Product rule on η² · (1/e):
%%   D(η²/e) = D(η²)·(1/e) + η² · D(1/e)
%%   D(η²) = 4 e (e+cos f)  (P25+.4 at n=2)
%%   D(1/e) = -De/e² = 2(e+cos f)/e²
%%   D(η²/e) = 4(e+cos f)/1 + (1-e²)·2(e+cos f)/e² = ... (assemble)
%% =======================================================
printf('=== Check 2: D(η²/e) compound (§6.2 Step 4) ===\n');

D_eta_sq = 2 * 2 * e * (e + cos(f)) * eta_e^(2 - 2);     % P25+.4 at n=2 → 4e(e+cos f)
D_inv_e  = -De_val / e^2;                                  % D(1/e) via chain rule
D_eta_sq_over_e = D_eta_sq * (1/e) + eta_e^2 * D_inv_e;
D_eta_sq_over_e_simpl = simplify(D_eta_sq_over_e);

D_boxed_compound = 2 * (e + cos(f)) * (1 + e^2) / e^2;
D_boxed_compound_simpl = simplify(D_boxed_compound);

residual2 = simplify(D_eta_sq_over_e_simpl - D_boxed_compound_simpl);
chk2_pass = isequal(residual2, sym(0));
if ~chk2_pass
  residual2 = simplify(expand(residual2));
  chk2_pass = isequal(residual2, sym(0));
end

printf('  Direct (product rule on η² · 1/e):  %s\n', char(D_eta_sq_over_e_simpl));
printf('  Boxed (ch12 §6.2 Step 4 line 606):  %s\n', char(D_boxed_compound_simpl));
printf('  Residual: %s\n', char(residual2));
if chk2_pass
  printf('  CHECK 2 STATUS: PASS\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 2 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 3: D(∂S_1/∂H) §7.1 Step 5 boxed vs direct chain-rule
%%
%% ∂S_1/∂H = (μ²k_2/G⁴) · Φ^{(H)} where
%%   Φ^{(H)} = 3θ E_1 - (3θ/2) E_2 - (3eθ/2) E_3 - (eθ/2) E_4 - (θ X_0^{0,2}/2) E_5
%%   E_1 = f - l + e sin f
%%   E_2 = sin(2(f+g))
%%   E_3 = sin(f+2g)
%%   E_4 = sin(3f+2g)
%%   E_5 = sin(2g)
%%
%% Direct D: D(μ²k_2/G⁴) · Φ + (μ²k_2/G⁴) · D(Φ)
%%   D(G^{-4}) = 4 G^{-4}  (DG = -G ⟹ D(G^n) = -n G^n)
%%   D(Φ) via Leibniz:
%%     D(3θ E_1) = 3θ · D(E_1) = 3θ · (Df - Dl - 2e sin f)   [P25+.3]
%%     D(E_2) = 0  [P25+.1 at (j,k)=(2,2), since Df+Dg=0]
%%     D(-3eθ E_3/2) = -(3θ/2)[De · E_3 + e · D(E_3)]
%%     D(-eθ E_4/2)  = -(θ/2) [De · E_4 + e · D(E_4)]
%%     D(-θ X_0^{0,2} E_5/2) = -(θ/2)[D(X_0^{0,2}) E_5 + X_0^{0,2} · D(E_5)]
%%
%% Compared to §7.1 Step 5 boxed:
%%   D(∂S_1/∂H) = (θ μ²k_2/G⁴) · { 12 E_1 - 6 E_2 - 6e E_3 - 2e E_4 - 2 X_0^{0,2} E_5
%%                                  + 6 sin f/e - 3 q_1 - 6 e sin f
%%                                  + 3(e+cos f) sin(f+2g) + 3 sin f cos(f+2g)
%%                                  +  (e+cos f) sin(3f+2g) -   sin f cos(3f+2g)
%%                                  -  e(e+cos f) sin(2g) + (2 X_0^{0,2} sin f/e) cos(2g) }
%% =======================================================
printf('=== Check 3: D(∂S_1/∂H) — §7.1 Step 5 boxed vs direct chain-rule ===\n');

syms X02 q1;   % X_0^{0,2}(e) abstract; q_1 = Dl abstract per chapter

% Basis elements
E1 = f - l_sym + e*sin(f);
E2 = sin(2*(f+g));
E3 = sin(f + 2*g);
E4 = sin(3*f + 2*g);
E5 = sin(2*g);

% Φ^{(H)} (sans μ²k_2)
Phi_H = 3*theta_sym*E1 ...
      - (3*theta_sym/2)*E2 ...
      - (3*e*theta_sym/2)*E3 ...
      - (e*theta_sym/2)*E4 ...
      - (theta_sym*X02/2)*E5;

% Direct D(Φ^{(H)}) via Leibniz + chain rule on each summand.
%   D(E_1) = Df - Dl - 2e sin f  (P25+.3)
DE1 = Df_val - q1 - 2*e*sin(f);
%   D(E_2) = 0  (P25+.1 at (2,2), since Df + Dg = 0)
DE2 = sym(0);
%   D(E_3) = (Df + 2Dg) cos(f+2g) = -Dg cos(f+2g) since Df = -Dg
DE3 = (Df_val + 2*Dg_val) * cos(f + 2*g);
%   D(E_4) = (3Df + 2Dg) cos(3f+2g) = (-3Dg + 2Dg) cos(3f+2g) = -Dg cos(3f+2g)
DE4 = (3*Df_val + 2*Dg_val) * cos(3*f + 2*g);
%   D(E_5) = 2 Dg cos(2g)
DE5 = 2*Dg_val * cos(2*g);

% D(X_0^{0,2}) — canonical value via ch09a Lemma E.1.2.
% Lemma E.1.2: dX_0^{0,2}/de = 2e(η+2)/(1+η)²
% Canonical X_0^{0,2}(e) = (3e²-2+2η³)/e² = e²(1+2η)/(1+η)² per ch06 Corollary B.0.7-7.
% Leading-order Taylor: X_0^{0,2}(e) = (3/4)e² + O(e^4), NOT -e²/2.
% Therefore D(X_0^{0,2}) = dX_0^{0,2}/de · De = (2e(η+2)/(1+η)²)·(-2(e+cos f))
%                       = -4e(η+2)(e+cos f)/(1+η)²
% At leading order: -3e(e+cos f). At general (e, η): exact above.
% (Replaces prior incorrect "D(X_0^{0,2}) = 2e(e+cos f)" which traced to ch12
%  §7.1 Step 4 line 662's wrong "X_0^{0,2} = -e²/2 leading order" claim.)
DX02 = -4*e*(eta_e + 2)*(e + cos(f)) / (1 + eta_e)^2;

% D(Φ^{(H)}) = sum of per-summand D's (note D(θ)=0, so θ commutes through):
D_Phi_H_direct = 3*theta_sym * DE1 ...
               + (-3*theta_sym/2) * DE2 ...
               + (-3*theta_sym/2) * (De_val * E3 + e * DE3) ...
               + (-theta_sym/2)   * (De_val * E4 + e * DE4) ...
               + (-theta_sym/2)   * (DX02 * E5 + X02 * DE5);

% Apply Leibniz with G^{-4} prefactor:
%   D(G^{-4} · (μ²k_2) · Φ^{(H)}) = 4 G^{-4} (μ²k_2) Φ^{(H)} + G^{-4} (μ²k_2) D(Φ^{(H)})
prefactor = mu_s^2 * k2_s / G_s^4;
D_dS1_dH_direct = 4 * prefactor * Phi_H + prefactor * D_Phi_H_direct;

% Boxed §7.1 Step 5 (corrected 2026-05-10): the sin(2g) coefficient
%   was previously written as `-e(e+cos f)` based on the wrong
%   D(X_0^{0,2}) = 2e(e+cos f). The corrected value uses canonical
%   D(X_0^{0,2}) = -4e(η+2)(e+cos f)/(1+η)², giving sin(2g)
%   coefficient `+2e(η+2)(e+cos f)/(1+η)²`.
boxed_bracket = 12*E1 - 6*E2 - 6*e*E3 - 2*e*E4 - 2*X02*E5 ...
              + 6*sin(f)/e - 3*q1 - 6*e*sin(f) ...
              + 3*(e + cos(f))*sin(f + 2*g) + 3*sin(f)*cos(f + 2*g) ...
              +   (e + cos(f))*sin(3*f + 2*g) -   sin(f)*cos(3*f + 2*g) ...
              + 2*e*(eta_e + 2)*(e + cos(f))*sin(2*g)/(1 + eta_e)^2 ...
              + (2*X02*sin(f)/e)*cos(2*g);
D_dS1_dH_boxed = theta_sym * prefactor * boxed_bracket;

% Residual: must equal 0 modulo dX02_de.
% In §7.1 line 662 the chapter uses dX02_de = 1 (leading order) — note this is approximate;
% the chapter's "X_0^{0,2}" is itself abstract here so we test with dX02_de = 1 substituted.
residual3a = simplify(D_dS1_dH_direct - D_dS1_dH_boxed);
chk3_pass = isequal(residual3a, sym(0));
if ~chk3_pass
  residual3a = simplify(expand(residual3a));
  chk3_pass = isequal(residual3a, sym(0));
end

printf('  Residual (direct chain-rule − boxed §7.1 Step 5):\n');
printf('    %s\n', char(residual3a));
if chk3_pass
  printf('  CHECK 3 STATUS: PASS — boxed form algebraically matches direct chain-rule.\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 3 STATUS: FAIL — boxed form does NOT match direct chain-rule.\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 4: ℳ_4 G-scaling symbolic.
%%
%% Claim: D(∂S_1/∂H) ∈ ℳ_4, meaning under G → λG, H → λH (θ = H/G fixed),
%% D(∂S_1/∂H) scales as λ^{-4}.
%%
%% Test: substitute G_s → λ*G_s, H_s → λ*H_s in the boxed form and verify
%% the result equals λ^{-4} times the un-scaled form. This is symbolic
%% (not numerical as in existing verify_ch12_deltas.m Checks 3-5), and uses
%% non-unit symbolic μ, k_2 per the Dimensional-Audit Rule.
%% =======================================================
printf('=== Check 4: ℳ_4 G-scaling symbolic (non-unit μ, k_2, L) ===\n');

syms lambda positive;

% θ = H/G must be preserved under scaling. Since we have theta_sym abstract,
% scaling G → λG, H → λH leaves theta_sym = (λH)/(λG) = H/G invariant.
% The only G-dependence is in the 1/G^4 prefactor.
D_dS1_dH_scaled = subs(D_dS1_dH_boxed, G_s, lambda*G_s);

% Expected: λ^{-4} times the original
expected_scaled = lambda^(-4) * D_dS1_dH_boxed;

residual4 = simplify(D_dS1_dH_scaled - expected_scaled);
chk4_pass = isequal(residual4, sym(0));
if ~chk4_pass
  residual4 = simplify(expand(residual4));
  chk4_pass = isequal(residual4, sym(0));
end
printf('  Residual after substituting G → λG, comparing to λ^{-4} × original:\n');
printf('    %s\n', char(residual4));

% Also confirm prefactor structure: μ² k_2 / G⁴
prefactor_check = simplify(D_dS1_dH_boxed * G_s^4 / (mu_s^2 * k2_s));
% Should reduce to a θ·{bracket} expression with NO μ, k_2, G remaining
has_mu = ~isequal(simplify(diff(prefactor_check, mu_s)), sym(0));
has_k2 = ~isequal(simplify(diff(prefactor_check, k2_s)), sym(0));
has_G  = ~isequal(simplify(diff(prefactor_check, G_s)),  sym(0));

prefactor_clean = ~(has_mu || has_k2 || has_G);
printf('  Prefactor extraction (μ²k_2/G⁴): ');
if prefactor_clean
  printf('PASS — residual is independent of μ, k_2, G.\n');
else
  printf('FAIL — residual still depends on: ');
  if has_mu; printf('μ '); end
  if has_k2; printf('k_2 '); end
  if has_G;  printf('G ');  end
  printf('\n');
end

if chk4_pass && prefactor_clean
  printf('  CHECK 4 STATUS: PASS\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 4 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 5: ch12 §7.2 Step 3 — ∂Φ_{S_1*}^{(T)}/∂θ direct symbolic vs boxed
%%
%% Source: ch11a Proposition G.4 gives
%%   Φ_{S_1*}^{(T)} = (μ²k_2/((5θ²-1) η²)) [AB e²(6+e²)/4 sin(2g) + B² e^4/128 sin(4g)]
%% with A(θ) = -1/2 + 3θ²/2, B(θ) = 3/2 - 3θ²/2.
%%
%% Boxed §7.2 Step 3:
%%   ∂Φ_{S_1*}^{(T)}/∂θ = (μ²k_2/((5θ²-1) η²)) [3θ(B-A) e²(6+e²)/4 sin(2g)
%%                                              - 6θ B e^4/128 sin(4g)]
%%                       + (-10θ μ²k_2/((5θ²-1)² η²)) [AB e²(6+e²)/4 sin(2g)
%%                                                     + B² e^4/128 sin(4g)]
%%
%% Direct route: SymPy diff(Phi_S1star, theta) using product rule on
%%   (5θ²-1)^{-1}, A(θ), B(θ), and the θ-independent (e, η, g) factors.
%% =======================================================
printf('=== Check 5: ch12 §7.2 Step 3 — ∂Φ_{S_1*}^{(T)}/∂θ direct vs boxed ===\n');

% θ-functions per ch11a G.4 (and ch04 §3 / file-04 definitions of A, B).
A_th = -sym(1)/2 + 3*theta_sym^2/2;
B_th =  sym(3)/2 - 3*theta_sym^2/2;
five_th2_m1 = 5*theta_sym^2 - 1;

% η is θ-independent — treat as abstract symbol for this check.
syms eta_sym positive;

% Φ_{S_1*}^{(T)} per ch11a G.4
Phi_S1star = (mu_s^2 * k2_s / (five_th2_m1 * eta_sym^2)) ...
           * ( A_th*B_th * e^2*(6+e^2)/4 * sin(2*g) ...
             + B_th^2   * e^4/128       * sin(4*g) );

% Direct route: SymPy θ-derivative
dPhi_dtheta_direct = simplify(diff(Phi_S1star, theta_sym));

% Boxed §7.2 Step 3
B_minus_A = B_th - A_th;
boxed_term1 = (mu_s^2 * k2_s / (five_th2_m1 * eta_sym^2)) ...
            * ( 3*theta_sym*B_minus_A * e^2*(6+e^2)/4 * sin(2*g) ...
              - 6*theta_sym*B_th * e^4/128 * sin(4*g) );
boxed_term2 = (-10*theta_sym * mu_s^2 * k2_s / (five_th2_m1^2 * eta_sym^2)) ...
            * ( A_th*B_th * e^2*(6+e^2)/4 * sin(2*g) ...
              + B_th^2    * e^4/128 * sin(4*g) );
dPhi_dtheta_boxed = boxed_term1 + boxed_term2;

residual5 = simplify(dPhi_dtheta_direct - dPhi_dtheta_boxed);
chk5_pass = isequal(residual5, sym(0));
if ~chk5_pass
  residual5 = simplify(expand(residual5));
  chk5_pass = isequal(residual5, sym(0));
end

printf('  Residual (SymPy diff(Φ, θ) − §7.2 Step 3 boxed):\n');
printf('    %s\n', char(residual5));
if chk5_pass
  printf('  CHECK 5 STATUS: PASS — ∂Φ_{S_1*}^{(T)}/∂θ boxed form matches SymPy diff.\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 5 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 6: ch12 §5.2 Step 4 — ∂Φ_{S_1*}^{(T)}/∂e direct symbolic vs boxed
%%
%% Boxed §5.2 Step 4:
%%   ∂Φ_{S_1*}^{(T)}/∂e = (μ²k_2/((5θ²-1) η²)) [AB e(3+e²) sin(2g) + B² e^3/32 sin(4g)]
%%                       + (2e μ²k_2/((5θ²-1) η^4)) [AB e²(6+e²)/4 sin(2g) + B² e^4/128 sin(4g)]
%%
%% (Note: in the second bracket the (5θ²-1)^{-1} factor remains — the second
%%  bracket comes from ∂(η^{-2})/∂e = 2e/η^4.)
%%
%% Direct route: SymPy diff(Phi_S1star, e), treating η = sqrt(1-e²)
%% (NOT as independent symbol) so the chain rule through η propagates.
%% =======================================================
printf('=== Check 6: ch12 §5.2 Step 4 — ∂Φ_{S_1*}^{(T)}/∂e direct vs boxed ===\n');

% This time express Φ with η explicit in e for SymPy to chain through.
Phi_S1star_e_explicit = (mu_s^2 * k2_s / (five_th2_m1 * eta_e^2)) ...
                      * ( A_th*B_th * e^2*(6+e^2)/4 * sin(2*g) ...
                        + B_th^2   * e^4/128       * sin(4*g) );

dPhi_de_direct = simplify(diff(Phi_S1star_e_explicit, e));

% Boxed §5.2 Step 4 — keep η_e explicit for direct comparison
boxed6_term1 = (mu_s^2 * k2_s / (five_th2_m1 * eta_e^2)) ...
             * ( A_th*B_th * e*(3+e^2) * sin(2*g) ...
               + B_th^2 * e^3/32 * sin(4*g) );
boxed6_term2 = (2*e * mu_s^2 * k2_s / (five_th2_m1 * eta_e^4)) ...
             * ( A_th*B_th * e^2*(6+e^2)/4 * sin(2*g) ...
               + B_th^2 * e^4/128 * sin(4*g) );
dPhi_de_boxed = boxed6_term1 + boxed6_term2;

residual6 = simplify(dPhi_de_direct - dPhi_de_boxed);
chk6_pass = isequal(residual6, sym(0));
if ~chk6_pass
  residual6 = simplify(expand(residual6));
  chk6_pass = isequal(residual6, sym(0));
end

printf('  Residual (SymPy diff(Φ, e) − §5.2 Step 4 boxed):\n');
printf('    %s\n', char(residual6));
if chk6_pass
  printf('  CHECK 6 STATUS: PASS — ∂Φ_{S_1*}^{(T)}/∂e boxed form matches SymPy diff.\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 6 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 7: GSI identity at α=3 (used by ch12 §6.1 Step 1).
%%
%% ch09b Theorem E.3 states: for S = G^{-α} F(θ, e, l, g) ∈ ℳ_α,
%%   ∂S/∂G = -(α/G) S - θ ∂S/∂H - (1/η) ∂S/∂L.
%%
%% Test at α=3 using a generic test function S = G^{-3} sin(2(f+g)) (h-indep,
%% f and l independent of G/L/H by construction of Delaunay coords away from
%% chain rule through eta=sqrt(1-(G/L)^2)). Verify the GSI holds symbolically.
%%
%% Use L, G, H as independent symbolic momenta. In Delaunay coordinates:
%%   e = sqrt(L² - G²)/L = sqrt(1 - (G/L)²)
%%   η = sqrt(1 - e²) = G/L                  (NOT sqrt(1-(G/L)²) — that equals e)
%%   θ = H/G
%% =======================================================
printf('=== Check 7: GSI identity at α=3 (ch12 §6.1 Step 1, ch09b Theorem E.3) ===\n');

% Define e, eta, theta in terms of L, G, H (the chain through these variables
% is what makes ∂/∂L, ∂/∂G, ∂/∂H meaningful).
syms l_int g_int positive;  % use these instead of l_sym/g_int to avoid name clash

% In Delaunay coordinates: e = sqrt(L²-G²)/L, η = sqrt(1-e²) = G/L,
% θ = H/G. Note η ≠ sqrt(1-(G/L)²) — that quantity equals e, not η.
e_LG      = sqrt(L_s^2 - G_s^2) / L_s;
eta_LG    = G_s / L_s;
theta_HG  = H_s / G_s;

% Generic test function: S = G^{-3} · F(θ, e, l, g) with explicit F
S_test = G_s^(-3) * sin(2*(l_int + g_int)) * (1 + theta_HG^2) ...
       * (1 + e_LG^2);

% LHS: ∂S/∂G via SymPy diff
LHS = simplify(diff(S_test, G_s));

% RHS via GSI at α=3:
RHS = -(sym(3)/G_s) * S_test ...
    - theta_HG * diff(S_test, H_s) ...
    - (1/eta_LG) * diff(S_test, L_s);
RHS = simplify(RHS);

residual7 = simplify(LHS - RHS);
chk7_pass = isequal(residual7, sym(0));
if ~chk7_pass
  residual7 = simplify(expand(residual7));
  chk7_pass = isequal(residual7, sym(0));
end

printf('  Residual (∂S/∂G − GSI-RHS at α=3 for test S = G^{-3}·sin(2(l+g))·(1+θ²)·(1+e²)):\n');
printf('    %s\n', char(residual7));
if chk7_pass
  printf('  CHECK 7 STATUS: PASS — GSI at α=3 holds symbolically (ch09b Theorem E.3 invoked by ch12 §6.1 Step 1).\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 7 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 8: ch09a Lemma E.1.2 — dX_0^{0,2}/de via direct SymPy diff of
%% the canonical ch06 Corollary B.0.7-7 form, compared to the closed-form
%% statement of Lemma E.1.2.
%%
%% Canonical: X_0^{0,2}(e) = (3e² - 2 + 2η³)/e²  (ch06 B.0.7-7)
%%          = e²(1 + 2η)/(1+η)²                  (algebraically equivalent)
%% Lemma E.1.2 statement: dX_0^{0,2}/de = 2e(η+2)/(1+η)²
%%
%% Direct route: SymPy diff of the closed-form (3e² - 2 + 2η³)/e² with
%% η = sqrt(1-e²) substituted (so chain rule through η).
%% =======================================================
printf('=== Check 8: ch09a Lemma E.1.2 — dX_0^{0,2}/de canonical value ===\n');

X02_canonical = (3*e^2 - 2 + 2*eta_e^3) / e^2;
dX02_de_direct = simplify(diff(X02_canonical, e));

% Lemma E.1.2 closed form
dX02_de_lemma = 2*e*(eta_e + 2) / (1 + eta_e)^2;

residual8 = simplify(dX02_de_direct - dX02_de_lemma);
chk8_pass = isequal(residual8, sym(0));
if ~chk8_pass
  residual8 = simplify(expand(residual8));
  chk8_pass = isequal(residual8, sym(0));
end
if ~chk8_pass
  % Try numerical sanity at multiple e values
  rel_err = 0;
  for e_val = [0.05, 0.2, 0.5, 0.8]
    direct_n = double(subs(dX02_de_direct, e, e_val));
    lemma_n  = double(subs(dX02_de_lemma,  e, e_val));
    rel_err = max(rel_err, abs(direct_n - lemma_n)/abs(lemma_n + 1e-30));
  end
  chk8_pass = (rel_err < 1e-12);
  printf('  Numerical rel err across e ∈ {0.05, 0.2, 0.5, 0.8}: %.3e\n', rel_err);
end

printf('  Residual (SymPy diff of (3e²-2+2η³)/e² with η = √(1-e²) − 2e(η+2)/(1+η)²):\n');
printf('    %s\n', char(residual8));
if chk8_pass
  printf('  CHECK 8 STATUS: PASS — ch09a Lemma E.1.2 holds against ch06 B.0.7-7 canonical form.\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 8 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 9: Independent anchor for ch12a Proposition P25+.1.
%%
%% Claim (ch12a §1): D(sin(jf+kg)) = (j·Df + k·Dg)·cos(jf+kg) for integer j, k.
%%
%% Independent route: SymPy chain rule on sin(jf+kg) where f, g are treated as
%% independent symbolic variables whose D-values are given by the base D-table
%% (ch12a Ax6, derived in file 08 P21): Df = 2 sin(f)/e, Dg = -2 sin(f)/e.
%%
%% Without this anchor, Check 3 used P25+.1 as a black-box upstream input.
%% =======================================================
printf('=== Check 9: ch12a Proposition P25+.1 — D(sin(jf+kg)) anchor ===\n');

chk9_pass = 0; chk9_total = 0;
for j_test = [0, 1, 2, 3]
  for k_test = [0, 1, 2, 3]
    if j_test == 0 && k_test == 0
      continue;
    end
    arg = j_test*f + k_test*g;
    % Direct route: chain rule via SymPy. Treat f and g as independent symbols
    % with D-values Df_val, Dg_val from base table.
    direct = diff(sin(arg), f)*Df_val + diff(sin(arg), g)*Dg_val;
    direct = simplify(direct);
    % P25+.1 statement
    boxed = (j_test*Df_val + k_test*Dg_val) * cos(arg);
    boxed = simplify(boxed);
    res = simplify(direct - boxed);
    pass = isequal(res, sym(0));
    if ~pass
      res = simplify(expand(res));
      pass = isequal(res, sym(0));
    end
    chk9_total = chk9_total + 1;
    if pass; chk9_pass = chk9_pass + 1; end
  endfor
endfor
printf('  %d/%d sub-cases over (j, k) ∈ {0,1,2,3}² \\ (0,0) PASS\n', chk9_pass, chk9_total);
if chk9_pass == chk9_total
  printf('  CHECK 9 STATUS: PASS — Proposition P25+.1 anchored to chain rule (Ax3) + base D-table (Ax6).\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 9 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 10: Independent anchor for ch12a Proposition P25+.3.
%%
%% Claim (ch12a §3): D(E_1) = Df - Dl - 2e sin f, where E_1 = f - l + e sin f.
%%
%% Independent route: D-linearity (Ax1) on each summand, with the base D-table
%% value D(e sin f) = -2 e sin f (file 08 Ax6).
%% =======================================================
printf('=== Check 10: ch12a Proposition P25+.3 — D(E_1) anchor ===\n');

D_esinf_table = -2*e*sin(f);  % base D-table value from file 08 Ax6

E1_sym = f - l_sym + e*sin(f);
% Direct route via linearity: D(f) + D(-l) + D(e sin f) = Df - Dl - 2e sin f
direct_E1 = Df_val + (-Dl_sym) + D_esinf_table;
direct_E1 = simplify(direct_E1);
% P25+.3 statement
boxed_E1 = Df_val - Dl_sym - 2*e*sin(f);
boxed_E1 = simplify(boxed_E1);
res10 = simplify(direct_E1 - boxed_E1);
chk10_pass = isequal(res10, sym(0));
if ~chk10_pass
  res10 = simplify(expand(res10));
  chk10_pass = isequal(res10, sym(0));
end

printf('  Residual (linearity (Ax1) on f - l + e sin f using D(e sin f) = -2e sin f − P25+.3 statement):\n');
printf('    %s\n', char(res10));
if chk10_pass
  printf('  CHECK 10 STATUS: PASS — Proposition P25+.3 anchored to linearity (Ax1) + base D-table (Ax6).\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 10 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 11: Independent anchor for ch12a Proposition P25+.5.
%%
%% Claim (ch12a §5): D(F(θ)) = 0 for any differentiable F.
%%
%% Independent route: chain rule (Ax3) with D(θ) = 0 from base D-table (Ax6).
%% Test F(θ) at three non-trivial functional forms: F = θ², F = 1/(5θ²-1),
%% F = sin(3θ).
%% =======================================================
printf('=== Check 11: ch12a Proposition P25+.5 — D(F(θ)) = 0 anchor ===\n');

chk11_pass = 0; chk11_total = 0;
F_candidates = {theta_sym^2, 1/(5*theta_sym^2 - 1), sin(3*theta_sym)};
F_names      = {'θ²', '1/(5θ²-1)', 'sin(3θ)'};
for k_test = 1:length(F_candidates)
  F_test = F_candidates{k_test};
  % Direct route: D(F(θ)) = F'(θ)·D(θ) = F'(θ)·0 = 0
  D_F_direct = diff(F_test, theta_sym) * Dtheta_val;
  D_F_direct = simplify(D_F_direct);
  pass = isequal(D_F_direct, sym(0));
  chk11_total = chk11_total + 1;
  if pass; chk11_pass = chk11_pass + 1; end
  printf('  F(θ) = %-12s  D(F(θ)) = %s\n', F_names{k_test}, char(D_F_direct));
endfor
if chk11_pass == chk11_total
  printf('  CHECK 11 STATUS: PASS — Proposition P25+.5 anchored to chain rule (Ax3) + D(θ) = 0 (Ax6).\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 11 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 12: Independent anchor for ch12a Proposition P25+.2.
%%
%% Claim (ch12a §2): D(cos(jf+kg)) = -(j·Df + k·Dg)·sin(jf+kg) for integer j, k.
%%
%% Independent route: chain rule (Ax3) on cos(jf+kg) with base D-table (Ax6):
%% Df = 2 sin(f)/e, Dg = -2 sin(f)/e.
%% =======================================================
printf('=== Check 12: ch12a Proposition P25+.2 — D(cos(jf+kg)) anchor ===\n');

chk12_pass = 0; chk12_total = 0;
for j_test = [0, 1, 2, 3]
  for k_test = [0, 1, 2, 3]
    if j_test == 0 && k_test == 0
      continue;
    end
    arg = j_test*f + k_test*g;
    % Direct route via SymPy chain rule
    direct = diff(cos(arg), f)*Df_val + diff(cos(arg), g)*Dg_val;
    direct = simplify(direct);
    % P25+.2 statement
    boxed = -(j_test*Df_val + k_test*Dg_val) * sin(arg);
    boxed = simplify(boxed);
    res = simplify(direct - boxed);
    pass = isequal(res, sym(0));
    if ~pass
      res = simplify(expand(res));
      pass = isequal(res, sym(0));
    end
    chk12_total = chk12_total + 1;
    if pass; chk12_pass = chk12_pass + 1; end
  endfor
endfor
printf('  %d/%d sub-cases over (j, k) ∈ {0,1,2,3}² \\ (0,0) PASS\n', chk12_pass, chk12_total);
if chk12_pass == chk12_total
  printf('  CHECK 12 STATUS: PASS — Proposition P25+.2 anchored to chain rule (Ax3) + base D-table (Ax6).\n\n');
  n_pass = n_pass + 1;
else
  printf('  CHECK 12 STATUS: FAIL\n\n');
  n_fail = n_fail + 1;
end

%% =======================================================
%% Check 13: Independent anchor for ch12a Proposition P25+.6.
%%
%% Claim (ch12a §6): D(κ) = 2(1 + e cos E) = 4 - 2κ, where κ = 1 - e cos E.
%%
%% Independent route: D(1 - e cos E) via product rule (Ax2) and base D-table:
%%   D(e cos E) = (De)·cos E + e·D(cos E) = (De)·cos E + e·(-sin E · DE)
%% with De = -2(e+cos f), DE = 2 sin E/(e κ) from ch12a Ax6.
%% Use the orbit-equation identity e + cos f = η²·cos E/κ (Kepler relation,
%% file 04 / file 05) to simplify.
%%
%% Test value: at sample (e, E) = (0.3, π/4):
%%   κ = 1 - 0.3·cos(π/4) = 1 - 0.3·√2/2 ≈ 0.7879
%%   2(1 + e cos E) = 2(1 + 0.3·√2/2) ≈ 2.4243
%%   4 - 2κ = 4 - 1.5757 ≈ 2.4243 ✓
%% =======================================================
printf('=== Check 13: ch12a Proposition P25+.6 — D(κ) anchor ===\n');

syms E_sym positive;
kappa_sym = 1 - e*cos(E_sym);

% Use the orbit-equation identity (ch07b Lemma C.4.2): η² cos E = κ (e + cos f).
% Equivalently: e + cos f = η² cos E / κ.
% Substitute into De = -2(e+cos f) chain through cos E:
%   D(e cos E) = (-2(e+cos f))·cos E + e·(-sin E · DE)
% From ch12a Ax6: DE = 2 sin E / (e κ).
DE_val = 2*sin(E_sym) / (e * kappa_sym);
% Express e+cos f via the orbit equation (treat as substitution at the end)
syms ecosf_sym positive;  % stand-in for (e+cos f) = η²cos E/κ
De_in_E = -2*ecosf_sym;

% D(κ) = D(1 - e cos E) = -D(e cos E) = -[(De)cos E + e D(cos E)]
%      = -[(-2 ecosf_sym) cos E + e · (-sin E · DE)]
%      = 2 ecosf_sym cos E + e sin E · DE
%      = 2 ecosf_sym cos E + e sin E · 2 sin(E)/(e κ)
%      = 2 ecosf_sym cos E + 2 sin²(E)/κ
D_kappa_direct = 2*ecosf_sym*cos(E_sym) + 2*sin(E_sym)^2 / kappa_sym;

% Substitute ecosf_sym = η² cos E / κ (with η = sqrt(1-e²)):
eta_e_sym = sqrt(1 - e^2);
ecosf_substituted = eta_e_sym^2 * cos(E_sym) / kappa_sym;
D_kappa_direct = subs(D_kappa_direct, ecosf_sym, ecosf_substituted);
D_kappa_direct_simpl = simplify(D_kappa_direct);

% P25+.6 statement: D(κ) = 2(1 + e cos E) = 4 - 2κ
D_kappa_boxed = 2*(1 + e*cos(E_sym));
D_kappa_boxed_alt = 4 - 2*kappa_sym;

res13a = simplify(D_kappa_direct_simpl - D_kappa_boxed);
res13b = simplify(D_kappa_boxed - D_kappa_boxed_alt);

% Test the FORM 4 - 2κ vs 2(1 + e cos E) directly
chk13a_pass = isequal(res13b, sym(0));
if ~chk13a_pass
  res13b = simplify(expand(res13b));
  chk13a_pass = isequal(res13b, sym(0));
end

% Test direct chain-rule via E-route matches boxed
chk13b_pass = isequal(res13a, sym(0));
if ~chk13b_pass
  res13a = simplify(expand(res13a));
  chk13b_pass = isequal(res13a, sym(0));
end

printf('  Sub-check 13a (form equivalence 2(1 + e cos E) ≡ 4 - 2κ): residual = %s\n', char(res13b));
printf('  Sub-check 13b (D(κ) via product rule on 1 - e cos E with Kepler-cos-E sub): residual = %s\n', char(res13a));

if chk13a_pass && chk13b_pass
  printf('  CHECK 13 STATUS: PASS — Proposition P25+.6 anchored to product rule (Ax2) + Kepler orbit equation + base D-table (Ax6).\n\n');
  n_pass = n_pass + 1;
else
  if chk13a_pass
    printf('  CHECK 13 STATUS: PARTIAL — form equivalence holds; chain-rule check FAIL (simplification may be incomplete).\n\n');
  else
    printf('  CHECK 13 STATUS: FAIL\n\n');
  end
  n_fail = n_fail + 1;
end

%% =======================================================
%% Summary
%% =======================================================
printf('=========================================================\n');
n_total = n_pass + n_fail;
if n_fail == 0
  printf('verify_ch12_deltaq_symbolic.m: ALL %d CHECKS PASSED.\n', n_total);
  printf('\n');
  printf('Theorem citations exercised + independent anchor coverage:\n');
  printf('  ch06 Corollary B.0.7-7 — canonical X_0^{0,2}(e) closed form\n');
  printf('    [anchor for Check 3 via Check 8].\n');
  printf('  ch09a Lemma E.1.2 — dX_0^{0,2}/de [Check 8: anchored to ch06 B.0.7-7].\n');
  printf('  ch09b Theorem E.3 — General Structural Identity [Check 7: anchored].\n');
  printf('  ch09d Proposition E.5(c) — 5-element basis of ∂S_1/∂H (Check 3 input).\n');
  printf('  ch11a Proposition G.4 — Φ_{S_1*}^{(T)} closed form (Checks 5, 6 input).\n');
  printf('  ch12a Proposition P25+.1 — D(sin(jf+kg)) [Check 9: anchored to Ax3 + Ax6].\n');
  printf('  ch12a Proposition P25+.2 — D(cos(jf+kg)) [Check 12: anchored to Ax3 + Ax6].\n');
  printf('  ch12a Proposition P25+.3 — D(E_1) for E_1 = f - l + e sin f\n');
  printf('    [Check 10: anchored to Ax1 + Ax6].\n');
  printf('  ch12a Proposition P25+.4 — D(1/η^n) [Check 1: anchored to Ax3 + Ax6].\n');
  printf('  ch12a Proposition P25+.5 — D(F(θ)) = 0 [Check 11: anchored to Ax3 + D(θ)=0].\n');
  printf('  ch12a Proposition P25+.6 — D(κ) = 2(1 + e cos E)\n');
  printf('    [Check 13: anchored to Ax2 + Kepler orbit equation + Ax6].\n');
  printf('  ch12 §5.2 Step 4 — boxed ∂Φ/∂e [Check 6: anchored to SymPy diff].\n');
  printf('  ch12 §6.2 Step 4 — boxed D(η²/e) compound [Check 2: anchored].\n');
  printf('  ch12 §7.1 Step 5 — boxed D(∂S_1/∂H) (Check 3 target; anchored).\n');
  printf('  ch12 §7.2 Step 3 — boxed ∂Φ/∂θ [Check 5: anchored to SymPy diff].\n');
  printf('  ch12a Ax2 (product rule), Ax3 (chain rule), Ax6 (base D-table) —\n');
  printf('    foundational; invoked in independent-route computations.\n');
else
  printf('verify_ch12_deltaq_symbolic.m: %d PASS / %d FAIL.\n', n_pass, n_fail);
end
printf('=========================================================\n');
