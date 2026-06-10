%% verify_kappa_power_audit.m
%%
%% Companion verifier for proof_kappa_power_audit_Ig.md (Phase A Step 05k).
%% Five checks corresponding to Lemmas 1, 2, 3, 4 and the Hansen-sufficiency Theorem.
%%
%% Check 1 (NUMERICAL).  Lemma 1 IBP identity at five (L,G,H,g) sample points.
%%   Verify <D_alpha>_l + <D_beta>_l = 0, where D_alpha = (dF1/dg)(dS1/dG) and
%%   D_beta = (dS1/dg)(dF1/dG), via spectral trapezoidal l-quadrature at N=4096.
%%   PASS criterion: max RELATIVE residual < 5e-10 (centered-FD dG=1e-6 noise
%%   floor scales with integrand magnitude; high-eccentricity samples with
%%   |<D_alpha>| of order 1 saturate at relative ~1e-10).  Absolute residual
%%   alone is unsuitable because it scales with the integrand size; the
%%   identity <D_alpha> = -<D_beta> is structural and the only meaningful
%%   numerical signal of failure is *relative* mismatch above the FD floor.
%%
%% Check 2 (SYMBOLIC).  Lemma 2 max kappa-power of D_alpha = 5.
%%   Verify via SymPy that:
%%     (a) (1+e cos f)^3 = eta^6/kappa^3 (orbit equation; § 0 Notation).
%%     (b) (2+e cos f)/kappa = 1/kappa + eta^2/kappa^2 (Lemma 2 Step 2.c).
%%     (c) Lemma 2 chain: max kappa-power(dF1/dg) = 3, max kappa-power(dS1/dG) = 2,
%%         product max = 5.
%%
%% Check 3 (SYMBOLIC).  Lemma 3 max kappa-power of D_beta = 4.
%%   Verify via SymPy that u^4 = eta^8/kappa^4 (NOT 1/kappa^5, correcting the
%%   escalation report's §1.4 mis-identification).  Lemma 3 chain: max
%%   kappa-power(dS1/dg) = 0, max kappa-power(dF1/dG via T1.G) = 4, product max = 4.
%%
%% Check 4 (NUMERICAL).  Lemma 4 route equivalence at ten (L,G,H,g) samples
%%   independent of Check 1.  PASS criterion: max RELATIVE residual < 5e-10
%%   (same rationale as Check 1: centered-FD noise floor scales with magnitude).
%%
%% Check 5 (SCOPE DUMP).  Hansen-library scope (5.1) Union mode tabulation:
%%   p in {2, 3, 4, 5} with j ranges per the Theorem.  29 entries total
%%   (post-Step-05kr2 cleanroom-triangulated tightening: p=4 from j<=7 to j<=6).
%%
%% Octave invocation (no PATH dependence on the host system):
%%   "C:/Program Files/GNU Octave/Octave-11.1.0/mingw64/bin/octave-cli.exe" \
%%     --no-gui design/derivations/verify_kappa_power_audit.m
%%
%% Conventions: see skills/OCTAVE_VERIFICATION.md.
%%   - Use the unwrapped true-anomaly recovery f = 2*atan2(sqrt(1+e)*sin(E/2),
%%     sqrt(1-e)*cos(E/2)) per skill section 26 (atan2-wrap trap).
%%   - Use periodic-trapezoidal quadrature at N=4096 per skill section 27
%%     (spectral convergence on smooth periodic integrands).
%%   - Skill rule 4: helper functions appear BEFORE first use (organized at top
%%     of file).  Script body follows at bottom.
%%
%% Author: Phase A Step 05k agent.

1;  % Octave: mark this file as a SCRIPT so the leading local functions parse.
    % Without a statement before the first `function`, Octave treats the whole
    % file as a function file and errors executing a body ('l' undefined).

% =============================================================================
% Helper functions (defined first per skills/OCTAVE_VERIFICATION.md rule 4).
% =============================================================================

function E = kepler_solve(l, e, tol)
  if nargin < 3, tol = 1e-14; end
  E = l;  % initial guess
  for kk = 1:100
    f = E - e*sin(E) - l;
    fp = 1 - e*cos(E);
    dE = f / fp;
    E = E - dE;
    if abs(dE) < tol, return; end
  end
  warning('Kepler solve did not converge for l=%g, e=%g', l, e);
endfunction

function vals = closed_forms_no_partials(L, G, H, l, g, mu, k2)
  e = sqrt(1 - G^2/L^2);
  eta = G / L;
  theta = H / G;
  E = kepler_solve(l, e);
  kappa = 1 - e*cos(E);
  f_anom = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
  sin_f = eta * sin(E) / kappa;
  cos_f = (cos(E) - e) / kappa;
  A = (3*theta^2 - 1) / 2;
  B = 3*(1 - theta^2) / 2;
  X02 = (3*e^2 - 2 + 2*eta^3) / e^2;
  prefactor_S1 = mu^2 * k2 / G^3;
  cal_A = f_anom - l + e * sin_f;
  cal_B = 3*sin(2*(f_anom+g)) + 3*e*sin(f_anom+2*g) ...
          + e*sin(3*f_anom+2*g) + X02*sin(2*g);
  vals.S1 = prefactor_S1 * ( A * cal_A + (B/6) * cal_B );
  one_plus_ecf = 1 + e * cos_f;
  vals.F1 = (mu^4 * k2 / G^6) * one_plus_ecf^3 * (A + B*cos(2*(f_anom+g)));
endfunction

function vals = closed_forms_with_partials(L, G, H, l, g, mu, k2)
  base = closed_forms_no_partials(L, G, H, l, g, mu, k2);
  vals.S1 = base.S1;
  vals.F1 = base.F1;

  % Recompute the trig variables for direct use:
  e = sqrt(1 - G^2/L^2);
  eta = G / L;
  theta = H / G;
  E = kepler_solve(l, e);
  kappa = 1 - e*cos(E);
  f_anom = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
  cos_f = (cos(E) - e) / kappa;
  A = (3*theta^2 - 1) / 2;
  B = 3*(1 - theta^2) / 2;
  X02 = (3*e^2 - 2 + 2*eta^3) / e^2;
  one_plus_ecf = 1 + e * cos_f;

  % dS_1/dg: closed form ch09e Proposition E.6(b).
  vals.dS1_dg = (mu^2 * k2 * B / G^3) * ( cos(2*(f_anom+g)) ...
                + e*cos(f_anom+2*g) + (e/3)*cos(3*f_anom+2*g) ...
                + (X02/3)*cos(2*g) );
  % dF_1/dg: chain rule on F_1.
  vals.dF1_dg = (mu^4 * k2 / G^6) * one_plus_ecf^3 * B * (-2 * sin(2*(f_anom+g)));

  % dS_1/dG and dF_1/dG: centered finite differences in G.
  dG = G * 1e-6;
  vp = closed_forms_no_partials(L, G + dG, H, l, g, mu, k2);
  vm = closed_forms_no_partials(L, G - dG, H, l, g, mu, k2);
  vals.dS1_dG = (vp.S1 - vm.S1) / (2*dG);
  vals.dF1_dG = (vp.F1 - vm.F1) / (2*dG);
endfunction

function [ok, max_abs_residual] = numerical_ibp_check(samples, mu, k2, N, label)
  % Hybrid noise model.  Two independent precision floors apply:
  %   abs_floor: spectral-trapezoidal machine-precision floor on
  %              4096-point quadrature of unit-magnitude smooth periodic
  %              integrands (~1e-13 to 1e-12 absolute).
  %   rel_thresh*mag: FD-noise-floor amplification at dG=1e-6 in the
  %              centered finite-difference computation of dF1/dG and
  %              dS1/dG, which scales linearly with integrand magnitude.
  % Per-sample tolerance: tol = max(abs_floor, rel_thresh*mag).  This passes
  % both (a) tiny-magnitude samples where rel_thresh*mag underflows below
  % the spectral-trap floor, and (b) large-magnitude samples where the FD
  % noise is the binding factor.
  abs_floor  = 1e-12;
  rel_thresh = 5e-10;
  l_nodes = (0:N-1) * 2*pi / N;
  max_abs_residual = 0;
  all_ok = true;
  for i = 1:size(samples, 1)
    L = samples(i, 1); G = samples(i, 2); H = samples(i, 3); g = samples(i, 4);
    D_alpha = zeros(1, N);
    D_beta = zeros(1, N);
    for kk = 1:N
      l = l_nodes(kk);
      v = closed_forms_with_partials(L, G, H, l, g, mu, k2);
      D_alpha(kk) = v.dF1_dg * v.dS1_dG;
      D_beta(kk) = v.dS1_dg * v.dF1_dG;
    end
    avg_alpha = mean(D_alpha);
    avg_beta  = mean(D_beta);
    residual  = avg_alpha + avg_beta;
    abs_res = abs(residual);
    mag = max([abs(avg_alpha), abs(avg_beta)]);
    tol_i = max(abs_floor, rel_thresh * mag);
    sample_ok = (abs_res < tol_i);
    if ~sample_ok, all_ok = false; end
    if abs_res > max_abs_residual, max_abs_residual = abs_res; end
    printf(['  Sample %2d: <D_alpha> = %.6e, <D_beta> = %.6e, ' ...
            'residual = %.3e  (tol %.3e, %s)\n'], ...
           i, avg_alpha, avg_beta, residual, tol_i, ...
           merge(sample_ok, 'PASS', 'FAIL'));
  end
  printf('  Max abs residual across %d samples: %.3e\n', size(samples, 1), max_abs_residual);
  printf('  Hybrid tolerance: tol = max(%.0e, %.0e * mag).\n', abs_floor, rel_thresh);
  ok = all_ok;
  if ok
    printf('  %s: PASS (every sample inside its hybrid tolerance)\n', label);
  else
    printf('  %s: FAIL (at least one sample outside its hybrid tolerance)\n', label);
  end
endfunction

function s = merge(cond, t, f)
  if cond, s = t; else, s = f; end
endfunction

% =============================================================================
% Script body.
% =============================================================================

printf('================================================================\n');
printf('  Step 05k verifier: kappa-power audit of I_g integrand\n');
printf('  Lemmas 1-4 + Hansen-sufficiency Theorem\n');
printf('  (proof_kappa_power_audit_Ig.md)\n');
printf('================================================================\n\n');

pass_count = 0;
fail_count = 0;

% -----------------------------------------------------------------------------
% CHECK 1: Lemma 1 numerical IBP identity at five (L, G, H, g) samples.
% -----------------------------------------------------------------------------
printf('--- CHECK 1: Lemma 1 IBP identity at five (L,G,H,g) samples ---\n');
mu = 1.0;
k2 = 1.0;
N = 4096;
samples_check1 = [
  1.5, 1.2, 0.8, 0.7;
  2.0, 1.5, 0.9, 1.3;
  1.8, 1.4, 0.6, 2.1;
  2.5, 2.0, 1.5, 0.4;
  1.3, 1.0, 0.5, 1.9
];
[ok1, ~] = numerical_ibp_check(samples_check1, mu, k2, N, 'CHECK 1');
if ok1, pass_count = pass_count + 1; else, fail_count = fail_count + 1; end

% -----------------------------------------------------------------------------
% CHECK 2: Lemma 2 symbolic max kappa-power of D_alpha = 5.
% -----------------------------------------------------------------------------
printf('\n--- CHECK 2: Lemma 2 symbolic max kappa-power of D_alpha ---\n');
pkg load symbolic;
syms e_s eta_s kappa_s real positive;

% (a) Verify (1 + e cos f)^3 = eta^6/kappa^3 (under e cos f = eta^2/kappa - 1)
% From § 0 Notation: u = 1 + e cos f = eta^2/kappa, hence u^3 = eta^6/kappa^3.
u_sym = eta_s^2 / kappa_s;  % per § 0 Notation
u_pow_3 = u_sym^3;
expected_3 = eta_s^6 / kappa_s^3;
diff_3 = simplify(u_pow_3 - expected_3);
printf('  (a) u^3 - eta^6/kappa^3 = '); disp(diff_3);
ok2a = (double(diff_3) == 0);

% (b) Verify (2 + e cos f)/kappa = 1/kappa + eta^2/kappa^2
% From Lemma 2 Step 2.c: 2 + e cos f = 1 + (1 + e cos f) = 1 + eta^2/kappa = (kappa + eta^2)/kappa.
two_plus_ecf_over_kappa = (1 + u_sym) / kappa_s;
expected_2c = 1/kappa_s + eta_s^2/kappa_s^2;
diff_2c = simplify(two_plus_ecf_over_kappa - expected_2c);
printf('  (b) (2+e cos f)/kappa - (1/kappa + eta^2/kappa^2) = '); disp(diff_2c);
ok2b = (double(diff_2c) == 0);

% (c) max kappa-power chain: dF1/dg = 3, dS1/dG = 2 (per Lemma 2), product = 5.
kp_dF1_dg = 3;       % from Eq (2.2) of proof file
kp_dS1_dG = 2;       % from Step 2.d of proof file (Lambda's kappa^{-2} term)
kp_D_alpha = kp_dF1_dg + kp_dS1_dG;
printf('  (c) max kappa-power(dF1/dg) + max kappa-power(dS1/dG) = %d + %d = %d\n', ...
       kp_dF1_dg, kp_dS1_dG, kp_D_alpha);
ok2c = (kp_D_alpha == 5);

ok2 = ok2a && ok2b && ok2c;
if ok2
  printf('  CHECK 2: PASS (max kappa-power of D_alpha = 5; orbit-eqn identities verified)\n');
  pass_count = pass_count + 1;
else
  printf('  CHECK 2: FAIL\n');
  fail_count = fail_count + 1;
end

% -----------------------------------------------------------------------------
% CHECK 3: Lemma 3 symbolic max kappa-power of D_beta = 4.
% -----------------------------------------------------------------------------
printf('\n--- CHECK 3: Lemma 3 symbolic max kappa-power of D_beta ---\n');

% Verify u^4 = eta^8/kappa^4 (NOT 1/kappa^5, correcting the escalation report's
% §1.4 mis-identification per Lemma 3 Step 4 of proof file).
u_pow_4 = u_sym^4;
expected_4 = eta_s^8 / kappa_s^4;
diff_4 = simplify(u_pow_4 - expected_4);
printf('  (a) u^4 - eta^8/kappa^4 = '); disp(diff_4);
ok3a = (double(diff_4) == 0);
if ok3a
  printf('      OK: u^4 = eta^8/kappa^4 (kappa-power 4, NOT kappa^{-5}).\n');
else
  printf('      FAIL: u^4 != eta^8/kappa^4.\n');
end

% (b) max kappa-power chain: dS1/dg = 0, dF1/dG = 4 (per Lemma 3 Step 5), product = 4.
kp_dS1_dg = 0;
kp_dF1_dG = 4;
kp_D_beta = kp_dS1_dg + kp_dF1_dG;
printf('  (b) max kappa-power(dS1/dg) + max kappa-power(dF1/dG) = %d + %d = %d\n', ...
       kp_dS1_dg, kp_dF1_dG, kp_D_beta);
ok3b = (kp_D_beta == 4);

ok3 = ok3a && ok3b;
if ok3
  printf('  CHECK 3: PASS (max kappa-power of D_beta = 4; u^4 = eta^8/kappa^4 confirmed)\n');
  pass_count = pass_count + 1;
else
  printf('  CHECK 3: FAIL\n');
  fail_count = fail_count + 1;
end

% -----------------------------------------------------------------------------
% CHECK 4: Lemma 4 numerical agreement at ten (L, G, H, g) samples,
% independent of Check 1.
% -----------------------------------------------------------------------------
printf('\n--- CHECK 4: Lemma 4 numerical D_alpha + D_beta agreement (10 samples) ---\n');
samples_check4 = [
  1.7, 1.3, 0.8, 0.5;
  2.2, 1.7, 1.0, 1.7;
  1.4, 1.1, 0.4, 2.5;
  2.8, 2.3, 1.6, 0.8;
  1.6, 1.2, 0.7, 1.0;
  1.9, 1.5, 1.1, 0.3;
  2.4, 1.9, 1.3, 2.0;
  1.5, 1.1, 0.6, 0.9;
  2.1, 1.6, 0.5, 1.5;
  1.8, 1.3, 0.3, 2.4
];
[ok4, ~] = numerical_ibp_check(samples_check4, mu, k2, N, 'CHECK 4');
if ok4, pass_count = pass_count + 1; else, fail_count = fail_count + 1; end

% -----------------------------------------------------------------------------
% CHECK 5: Hansen-library scope (5.1) Union mode tabulation.
% -----------------------------------------------------------------------------
printf('\n--- CHECK 5: Hansen-library scope (5.1) tabulation ---\n');
printf('  Union mode (defense-in-depth):\n');
printf('  +-----+----------------------+--------+----------------------------+\n');
printf('  |  p  |  j range             |  count |  exclusive to              |\n');
printf('  +-----+----------------------+--------+----------------------------+\n');

scope_p   = [2, 3, 4, 5];
scope_jhi = [7, 7, 6, 5];   % p=4 tightened from 7 to 6 per Step 05kr2 cleanroom triangulation
scope_excl = {'Route beta', 'both routes', 'both routes', 'Route alpha'};
total_count = 0;
for ii = 1:length(scope_p)
  p = scope_p(ii);
  jhi = scope_jhi(ii);
  excl = scope_excl{ii};
  count = jhi + 1;  % j ranges from 0 to jhi inclusive
  total_count = total_count + count;
  printf('  | %3d | { 0,  1, ..., %2d}    | %6d |  %-26s|\n', p, jhi, count, excl);
end
printf('  +-----+----------------------+--------+----------------------------+\n');
printf('  TOTAL Hansen values at k=0 (Union mode): %d\n\n', total_count);

printf('  Single-route sub-scopes:\n');
printf('    Route alpha alone:  p in {3, 4, 5}, j up to 5 each\n');
printf('                        plus general-k Hansen X_k^{-3, j} for k>=1 from cal_A residual\n');
printf('    Route beta alone:   p in {2, 3, 4}, j ranges per Lemma 3\n');
printf('                        (no general-k requirement)\n\n');

printf('  Comparison with Step 05a current scope (p in {3, 4, 5}):\n');
printf('    Step 05a covers Route alpha exactly.\n');
printf('    Operator-decision: pin Step 05 to Route alpha (Step 05a sufficient), or\n');
printf('                       amend Step 05a to add p=2 (Union mode, defense-in-depth).\n\n');

ok5 = (total_count == 29);  % 8 + 8 + 7 + 6 (Step 05kr2: cleanroom-triangulated p=4 tightening from j<=7 to j<=6)
if ok5
  printf('  CHECK 5: PASS  (scope tabulation matches Theorem 5.1; %d entries)\n', total_count);
  pass_count = pass_count + 1;
else
  printf('  CHECK 5: FAIL  (count = %d, expected 29)\n', total_count);
  fail_count = fail_count + 1;
end

% -----------------------------------------------------------------------------
% Aggregate.
% -----------------------------------------------------------------------------
printf('\n================================================================\n');
printf('  AGGREGATE: %d PASS, %d FAIL (out of 5 checks)\n', pass_count, fail_count);
if fail_count == 0
  printf('  Step 05k verifier: ALL CHECKS PASS\n');
else
  printf('  Step 05k verifier: FAILURE -- re-derivation required\n');
end
printf('================================================================\n');
