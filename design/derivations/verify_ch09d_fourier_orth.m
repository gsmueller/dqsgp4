% verify_ch09d_fourier_orth.m
%
% Step 06 verifier: ch09d Fourier orthogonality via Hansen averages.
%
% Companion to: sgp4_references/vallado_celestrak/documentation/SGP4/
%               Brouwer_Hori 1961/derivation/proof_ch09d_fourier_orthogonality.md
%
% Closes the DERIVATION_GAP in ch09d Remark E.5.1.3 (the "Fourier
% orthogonality" claim for the pure-Fourier basis elements
%   E_2 = sin(2f+2g),  E_3 = sin(f+2g),  E_4 = sin(3f+2g),  E_5 = sin(2g)).
%
% Four parts:
%
%   Part 1 - Closed-form Hansen coefficients (Lemma O.6 of proof file):
%     For J in {1, 2, 3}, verify that the claimed closed form
%     X_0^{0,J}(e) of Lemma O.6 matches the defining integral
%     X_0^{0,m}(e) = (1/(2π)) int_0^{2π} κ cos(m f) dE (ch06b B.1.2)
%     evaluated by periodic trapezoidal quadrature at e in
%     {0.1, 0.3, 0.5, 0.7}. J = 0 is trivial (1 identically).
%
%   Part 2 - Sin-sector Gram matrix M^{(s)} (Theorem O.5 of proof file):
%     Build the 4x4 Gram matrix from the closed forms of Part 1 at
%     e in {0.1, 0.3, 0.5, 0.7} and verify positive definiteness
%     (det > 0, smallest eigenvalue > 0).
%
%   Part 3 - Quadrature cross-check of the Gram matrix entries:
%     Evaluate ⟨E_a E_b⟩_{l, g} by periodic trapezoidal quadrature on
%     the (E, g)-torus with the Kepler weight κ = 1 - e cos E (the
%     pull-back of dl = κ dE, ch06b Theorem B.1; the g-measure is
%     already uniform). Trapezoidal is spectrally convergent for
%     smooth 2π-periodic integrands (Trefethen & Weideman 2014);
%     per OCTAVE_VERIFICATION §27 and Lesson L4 from the envelope's
%     lessons log, quadgk is unreliable here. Compare each numerical
%     entry against the closed-form prediction of Part 2; require
%     max relative error < 1e-10.
%
%   Part 4 - Sin-cos decoupling (Proposition O.5 of proof file):
%     Verify that ⟨sin(jf + kg) cos(j'f + k'g)⟩_{l, g} = 0 (machine
%     precision, < 1e-12) for all 4x4 sin-cos pairs at
%     (k, k') = (2, 2) and a sampling of (j, j') pairs drawn from
%     {0, 1, 2, 3}. This confirms the parity-based block structure
%     of the full 8x8 Gram matrix (sin-sector + cos-sector).
%
% Pure numerical verification (no SymPy dependency). Part 1 uses
% N = 4096 nodes (1D); Parts 3-4 use N = 1024 per dimension on the
% (E, g)-torus (total grid 1024^2 ≈ 1e6 evaluations per pair, across
% 4 eccentricities — well under a minute of wall time). Trapezoidal
% convergence on these smooth 2π-periodic integrands is spectral.
%
% Expected output: "Summary: 64 pass, 0 fail".
%   Part 1: 12 (X_0^{0,m} for m = 1, 2, 3 vs quadrature at 4 eccentricities).
%   Part 2: 12 (symmetric / det>0 / lambda_min>0 at 4 eccentricities).
%   Part 3: 40 (10 unordered pairs {a,b} in {2,3,4,5} x 4 eccentricities).
%   Part 4:  4 (max |sin-cos| < 1e-12 at 4 eccentricities, enveloping 16 pairs).

printf('===============================================================\n');
printf('Step 06 verifier: ch09d Fourier orthogonality (Hansen averages)\n');
printf('===============================================================\n\n');

n_pass = 0;
n_fail = 0;

function [np, nf] = assert_true(label, cond, note, already_pass, already_fail)
  if cond
    printf('  PASS: %s  [%s]\n', label, note);
    np = already_pass + 1;
    nf = already_fail;
  else
    printf('  FAIL: %s  [%s]\n', label, note);
    np = already_pass;
    nf = already_fail + 1;
  end
endfunction

% Closed-form X_0^{0, |J|}(e) from Lemma O.6 of the proof file, inline.
%   J = 0: X_0^{0,0} = 1
%   J = 1: X_0^{0,1} = -e
%   J = 2: X_0^{0,2} = (3 e^2 - 2 + 2 eta^3) / e^2
%   J = 3: X_0^{0,3} = (8 - 12 e^2 + 3 e^4 - 8 eta^3) / e^3
function v = X0(J, e)
  eta = sqrt(1 - e^2);
  switch abs(J)
    case 0, v = 1;
    case 1, v = -e;
    case 2, v = (3*e^2 - 2 + 2*eta^3) / e^2;
    case 3, v = (8 - 12*e^2 + 3*e^4 - 8*eta^3) / e^3;
    otherwise
      error('X0: |J| > 3 out of scope of Lemma O.6');
  end
endfunction

%% =========================================================
%% Part 1: Closed-form X_0^{0, J}(e) vs periodic trapezoidal quadrature.
%% =========================================================
%
% Per Lemma O.6 of proof_ch09d_fourier_orthogonality.md:
%   X_0^{0, 1}(e) = -e                                   (O.6-01)
%   X_0^{0, 2}(e) = (3 e^2 - 2 + 2 eta^3) / e^2          (O.6-2)
%   X_0^{0, 3}(e) = (8 - 12 e^2 + 3 e^4 - 8 eta^3) / e^3 (O.6-3)
%
% Verification: evaluate the defining integral
%   X_0^{0, m}(e) = (1/(2π)) * int_0^{2π} κ cos(m f) dE     (B.1.2)
% at the test eccentricities by spectrally-convergent periodic
% trapezoidal quadrature (per Lesson L4 of _AGENT_ENVELOPE.md),
% compare to the closed form; require |abs_err| < 1e-12.

printf('--- Part 1: Closed-form X_0^{0, J} vs quadrature (Lemma O.6) ---\n');

e_test = [0.1, 0.3, 0.5, 0.7];
N_1d = 4096;  % for 1D quadrature in Part 1
E1d = (0:N_1d-1) * (2*pi / N_1d);

for ee = e_test
  eta_num = sqrt(1 - ee^2);
  kap_1d = 1 - ee*cos(E1d);
  cos_f1d = (cos(E1d) - ee) ./ kap_1d;
  sin_f1d = eta_num * sin(E1d) ./ kap_1d;

  % cos(m f) via Chebyshev recurrence (wrap-immune, OCTAVE_VERIFICATION §28).
  c0 = ones(size(cos_f1d));
  c1 = cos_f1d;
  c2 = 2*cos_f1d .* c1 - c0;
  c3 = 2*cos_f1d .* c2 - c1;

  for m = 1:3
    switch m
      case 1, cos_mf = c1;
      case 2, cos_mf = c2;
      case 3, cos_mf = c3;
    end
    integrand = kap_1d .* cos_mf;
    num_val = sum(integrand) / N_1d;  % periodic trapezoidal
    closed_val = X0(m, ee);
    abs_err = abs(num_val - closed_val);
    [n_pass, n_fail] = ...
      assert_true(sprintf('X_0^{0,%d} closed vs quadrature at e=%.1f', m, ee), ...
                  abs_err < 1e-12, ...
                  sprintf('closed=%+.10e, quad=%+.10e, abs_err=%.2e', ...
                          closed_val, num_val, abs_err), n_pass, n_fail);
  end
end

printf('\n');

%% =========================================================
%% Part 2: Sin-sector Gram matrix positive definiteness.
%% =========================================================

printf('--- Part 2: Gram matrix M^{(s)} positive definiteness ---\n');

% (X0 is the closed-form helper function defined at the top of the file.)

% Difference table for E_2..E_5 with j-values (2, 1, 3, 0).
j_vec = [2, 1, 3, 0];

for ee = e_test
  % Build 4x4 Gram matrix M (O.5-b form, times 1/2).
  M = zeros(4, 4);
  for a = 1:4
    for b = 1:4
      diff_j = abs(j_vec(a) - j_vec(b));
      M(a, b) = 0.5 * X0(diff_j, ee);
    end
  end
  % M must be symmetric positive definite.
  is_symmetric = norm(M - M.', 'fro') < 1e-14;
  det_M = det(M);
  eigs_M = eig(M);
  lambda_min = min(eigs_M);

  [n_pass, n_fail] = ...
    assert_true(sprintf('M symmetric at e=%.1f', ee), is_symmetric, ...
                sprintf('|M-M^T|_F = %.2e', norm(M-M.','fro')), ...
                n_pass, n_fail);
  [n_pass, n_fail] = ...
    assert_true(sprintf('det(M) > 0 at e=%.1f', ee), det_M > 0, ...
                sprintf('det = %.6e', det_M), n_pass, n_fail);
  [n_pass, n_fail] = ...
    assert_true(sprintf('lambda_min(M) > 0 at e=%.1f', ee), lambda_min > 0, ...
                sprintf('lambda_min = %.6e', lambda_min), n_pass, n_fail);
end

printf('\n');

%% =========================================================
%% Part 3: Quadrature cross-check of Gram matrix entries.
%% =========================================================

printf('--- Part 3: Quadrature vs closed-form Gram entries ---\n');

% Trapezoidal quadrature on (E, g) with Kepler weight κ dE.
% Number of nodes: N = 1024 per dimension (1024^2 ≈ 1e6 evaluations; these
% integrands are analytic in (E, g) with analyticity strip > 0, so
% trapezoidal converges geometrically fast in N — machine precision
% reached well below N = 1024 even at e = 0.7).
N = 1024;
E_grid = (0:N-1) * (2*pi / N);
g_grid = (0:N-1) * (2*pi / N);
[EE, GG] = meshgrid(E_grid, g_grid);

for ee = e_test
  kap = 1 - ee*cos(EE);
  eta_num = sqrt(1 - ee^2);

  % Build f(E) via cos(f) = (cos E - e)/κ, sin(f) = η sin E / κ.
  cos_f = (cos(EE) - ee) ./ kap;
  sin_f = eta_num * sin(EE) ./ kap;

  % Build cos(j f), sin(j f) for j = 0..3 via Chebyshev-recurrence
  % (wrap-immune; no atan2 reconstruction — OCTAVE_VERIFICATION §28).
  cos_jf = cell(4, 1);
  sin_jf = cell(4, 1);
  cos_jf{1} = ones(size(cos_f));        % j = 0
  sin_jf{1} = zeros(size(cos_f));
  cos_jf{2} = cos_f;                     % j = 1
  sin_jf{2} = sin_f;
  % Chebyshev recurrence:  cos((j+1) x) = 2 cos(x) cos(j x) - cos((j-1) x),
  %                         sin((j+1) x) = 2 cos(x) sin(j x) - sin((j-1) x).
  for j = 2:3
    cos_jf{j+1} = 2*cos_f .* cos_jf{j} - cos_jf{j-1};
    sin_jf{j+1} = 2*cos_f .* sin_jf{j} - sin_jf{j-1};
  end

  % Build E_a = sin(j_a f + k_a g) for a = 2..5.
  % angle addition: sin(jf + kg) = sin(jf) cos(kg) + cos(jf) sin(kg).
  k_vec = [2, 2, 2, 2];
  E_arr = cell(4, 1);
  for a = 1:4
    j_a = j_vec(a);
    k_a = k_vec(a);
    sin_jf_a = sin_jf{j_a + 1};
    cos_jf_a = cos_jf{j_a + 1};
    sin_kg = sin(k_a * GG);
    cos_kg = cos(k_a * GG);
    E_arr{a} = sin_jf_a .* cos_kg + cos_jf_a .* sin_kg;
  end

  % Quadrature weight: (κ / (2π)) dE * (1 / (2π)) dg; with periodic
  % trapezoidal rule, ⟨φ⟩_{l,g} ≈ (1/N^2) sum(φ · κ).
  for a = 1:4
    for b = a:4
      integrand = E_arr{a} .* E_arr{b} .* kap;
      num_val = sum(integrand(:)) / N^2;
      diff_j = abs(j_vec(a) - j_vec(b));
      closed_val = 0.5 * X0(diff_j, ee);
      rel_err = abs(num_val - closed_val) / max(abs(closed_val), 1e-16);
      abs_err = abs(num_val - closed_val);
      % Accept if abs_err < 1e-10 (tight for trapezoidal + spectrally
      % convergent integrand) OR if both values < 1e-14 (both machine zero).
      is_zero = abs(closed_val) < 1e-14 && abs(num_val) < 1e-14;
      pass_cond = is_zero || abs_err < 1e-10;
      [n_pass, n_fail] = ...
        assert_true(sprintf('⟨E_%d E_%d⟩ quadrature vs closed form, e=%.1f', ...
                           a+1, b+1, ee), pass_cond, ...
                    sprintf('num=%+.6e, closed=%+.6e, abs_err=%.2e', ...
                            num_val, closed_val, abs_err), n_pass, n_fail);
    end
  end
end

printf('\n');

%% =========================================================
%% Part 4: Sin-cos decoupling.
%% =========================================================

printf('--- Part 4: Sin-cos decoupling (Proposition O.5) ---\n');

% Test pairs: sin(j f + 2 g) vs cos(j' f + 2 g) for j, j' in {0,1,2,3}.
for ee = e_test
  kap = 1 - ee*cos(EE);
  eta_num = sqrt(1 - ee^2);
  cos_f = (cos(EE) - ee) ./ kap;
  sin_f = eta_num * sin(EE) ./ kap;

  cos_jf = cell(4, 1);
  sin_jf = cell(4, 1);
  cos_jf{1} = ones(size(cos_f));
  sin_jf{1} = zeros(size(cos_f));
  cos_jf{2} = cos_f;
  sin_jf{2} = sin_f;
  for j = 2:3
    cos_jf{j+1} = 2*cos_f .* cos_jf{j} - cos_jf{j-1};
    sin_jf{j+1} = 2*cos_f .* sin_jf{j} - sin_jf{j-1};
  end

  % sin(j f + 2 g) and cos(j' f + 2 g)
  cos_2g = cos(2*GG);
  sin_2g = sin(2*GG);

  max_abs = 0;
  for j = 0:3
    for jp = 0:3
      sinA = sin_jf{j+1}  .* cos_2g + cos_jf{j+1}  .* sin_2g;
      cosB = cos_jf{jp+1} .* cos_2g - sin_jf{jp+1} .* sin_2g;
      integrand = sinA .* cosB .* kap;
      val = sum(integrand(:)) / N^2;
      max_abs = max(max_abs, abs(val));
    end
  end
  [n_pass, n_fail] = ...
    assert_true(sprintf('max|⟨sin(jf+2g) cos(j''f+2g)⟩| < 1e-12 at e=%.1f', ee), ...
                max_abs < 1e-12, ...
                sprintf('max |val| = %.3e', max_abs), n_pass, n_fail);
end

printf('\n');

%% =========================================================
%% Summary
%% =========================================================

printf('===============================================================\n');
printf('Summary: %d pass, %d fail\n', n_pass, n_fail);
printf('===============================================================\n');

if n_fail > 0
  error('verify_ch09d_fourier_orth: %d FAILURES', n_fail);
end

%% (X0 closed-form helper is defined at the top of the file.)
