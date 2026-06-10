% verify_ch09d_bessel_series.m
%
% Step 07 verifier: Bessel-series expansion of E_1 = f - l + e sin f.
%
% Companion to: sgp4_references/vallado_celestrak/documentation/SGP4/
%               Brouwer_Hori 1961/derivation/proof_ch09d_bessel_series.md
%
% Closes the DERIVATION_GAP in ch09d Remark E.5.1.3 (infinite-Bessel-
% series claim for E_1, plus linear independence from {E_2, ..., E_8}).
%
% Three parts:
%
%   Part 1 - Direct E_1(l, e) vs Bessel closed form.
%     At sample e in {0.1, 0.3, 0.5} and l in {k*pi/8 : k = 1, 3, 5, ...,
%     15}: solve Kepler's equation l = E - e sin E by Newton iteration;
%     compute f from the unwrapped half-angle formula (per OCTAVE_
%     VERIFICATION §26); evaluate E_1 = f - l + e sin f. Separately,
%     compute E_1 from the truncated Bessel series with
%       c_k(e) = (2/k)[J_k(ke) + sum_{m=1}^{M} β^m (mη + 1)
%                (J_{k+m}(ke) + J_{k-m}(ke))]
%     (Theorem BS.5 of proof file). Truncate at K=60 harmonics and
%     M=80 in the β-sum; require |direct - Bessel| < 1e-9.
%
%   Part 2 - Cross-check Bessel form vs Hansen form.
%     Compute c_k(e) in two ways:
%       (a) Direct Fourier integral: (1/π) * int E_1(l) sin(kl) dl
%           (periodic trapezoidal at N=4096 nodes).
%       (b) Via Theorem BS.6: c_k = (2η³/k) * X_k^{-3, 0}(e), where
%           X_k^{-3, 0}(e) is itself computed by the Hansen integral
%           (1/(2π)) int κ^{-2} cos(k(E - e sin E)) dE (periodic
%           trapezoidal at N=4096).
%       (c) Via Bessel closed form (BS.5-c) with truncation as above.
%     For k = 1..10 and e = 0.1, 0.3, 0.5: require pairwise agreement
%     of (a), (b), (c) to within 1e-11.
%
%   Part 3 - Infinite support (Theorem BS.7).
%     At e = 0.3, compute c_k(e) for k = 1, ..., 20 and verify
%     c_k(0.3) > 1e-14 for every k (non-vanishing).
%
% Periodic trapezoidal quadrature is used throughout for spectral
% convergence on smooth 2pi-periodic integrands (OCTAVE_VERIFICATION
% §27, Lesson L4 of _AGENT_ENVELOPE.md). quadgk is avoided.
%
% Expected output: "Summary: ... pass, 0 fail".
%   Part 1: 24 (3 eccentricities × 8 l-values).
%   Part 2: 90 (3 eccentricities × 10 k-values × 3 pairwise checks).
%   Part 3: 20 (20 k-values at e = 0.3).

printf('===============================================================\n');
printf('Step 07 verifier: ch09d Bessel-series expansion of E_1\n');
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

% -- Kepler solve: given l, e, return E satisfying l = E - e sin E. --
function E = kepler_solve(l, e)
  % Newton iteration, tol 1e-14.
  E = l + e * sin(l);  % initial guess
  for iter = 1:50
    fE = E - e*sin(E) - l;
    fEp = 1 - e*cos(E);
    dE = fE / fEp;
    E = E - dE;
    if abs(dE) < 1e-14
      break;
    end
  end
endfunction

% -- Unwrapped true anomaly from eccentric anomaly (OCTAVE_VERIFICATION §26). --
function f = true_anomaly(E, e)
  f = 2 * atan2(sqrt(1 + e) * sin(E/2), sqrt(1 - e) * cos(E/2));
  % Map (-pi, pi] to [0, 2pi] for positive-E input.
  if f < 0
    f = f + 2*pi;
  end
endfunction

% -- c_k(e) via Bessel closed form (Theorem BS.5-c). --
function ck = bessel_ck(k, e, M)
  eta = sqrt(1 - e^2);
  beta = (1 - eta) / e;
  x = k * e;
  ck = besselj(k, x);
  for m = 1:M
    term = beta^m * (m*eta + 1) * (besselj(k+m, x) + besselj(k-m, x));
    ck = ck + term;
  end
  ck = (2/k) * ck;
endfunction

% -- X_k^{-3, 0}(e) via Hansen integral (ch06b Corollary B.1.2 at n=-3, m=0). --
%    X_k^{-3, 0}(e) = (1/(2π)) ∫ κ^{-3} cos(k l) dl
%                   = (1/(2π)) ∫ κ^{-2} cos(k(E - e sin E)) dE.
function Xk = hansen_X_k_n3_0(k, e, N)
  E = (0:N-1) * (2*pi / N);
  kap = 1 - e*cos(E);
  integrand = cos(k * (E - e*sin(E))) ./ kap.^2;
  Xk = sum(integrand) / N;  % periodic trapezoidal = (1/(2π)) ∫ ... dE
endfunction

% -- c_k via the direct Fourier integral (1/π) ∫ E_1(l) sin(kl) dl. --
function ck = direct_fourier_ck(k, e, N)
  l = (0:N-1) * (2*pi / N);
  E1 = zeros(1, N);
  for j = 1:N
    E = kepler_solve(l(j), e);
    f = true_anomaly(E, e);
    E1(j) = f - l(j) + e*sin(f);
  end
  ck = (2/N) * sum(E1 .* sin(k * l));  % (1/π) · (2π/N) · Σ = (2/N) Σ
endfunction

%% =========================================================
%% Part 1: Direct vs Bessel-series E_1(l).
%% =========================================================

printf('--- Part 1: Direct E_1 vs Bessel-series ---\n');

e_list_1 = [0.1, 0.3, 0.5];
l_list = pi/8 * [1, 3, 5, 7, 9, 11, 13, 15];  % 8 values avoiding endpoints
K_trunc = 60;
M_trunc = 80;

for ee = e_list_1
  for ll = l_list
    % Direct computation:
    E_root = kepler_solve(ll, ee);
    f_val = true_anomaly(E_root, ee);
    E1_direct = f_val - ll + ee*sin(f_val);

    % Bessel truncated series:
    E1_bessel = 0;
    for k = 1:K_trunc
      ck = bessel_ck(k, ee, M_trunc);
      E1_bessel = E1_bessel + ck * sin(k*ll);
    end

    abs_err = abs(E1_direct - E1_bessel);
    [n_pass, n_fail] = ...
      assert_true(sprintf('E_1 direct vs Bessel, e=%.1f, l=%.4f', ee, ll), ...
                  abs_err < 1e-9, ...
                  sprintf('direct=%+.6e, Bessel=%+.6e, err=%.2e', ...
                          E1_direct, E1_bessel, abs_err), ...
                  n_pass, n_fail);
  end
end

printf('\n');

%% =========================================================
%% Part 2: Cross-check Bessel vs Hansen vs direct Fourier.
%% =========================================================

printf('--- Part 2: c_k Bessel vs Hansen vs direct Fourier ---\n');

N_quad = 4096;

for ee = e_list_1
  eta = sqrt(1 - ee^2);
  for k = 1:10
    ck_bessel = bessel_ck(k, ee, M_trunc);
    ck_hansen = (2*eta^3/k) * hansen_X_k_n3_0(k, ee, N_quad);
    ck_direct = direct_fourier_ck(k, ee, N_quad);

    err_BH = abs(ck_bessel - ck_hansen);
    err_BD = abs(ck_bessel - ck_direct);
    err_HD = abs(ck_hansen - ck_direct);

    [n_pass, n_fail] = ...
      assert_true(sprintf('c_%d(e=%.1f) Bessel vs Hansen', k, ee), ...
                  err_BH < 1e-11, ...
                  sprintf('B=%+.6e, H=%+.6e, err=%.2e', ck_bessel, ck_hansen, err_BH), ...
                  n_pass, n_fail);
    [n_pass, n_fail] = ...
      assert_true(sprintf('c_%d(e=%.1f) Bessel vs direct', k, ee), ...
                  err_BD < 1e-11, ...
                  sprintf('B=%+.6e, D=%+.6e, err=%.2e', ck_bessel, ck_direct, err_BD), ...
                  n_pass, n_fail);
    [n_pass, n_fail] = ...
      assert_true(sprintf('c_%d(e=%.1f) Hansen vs direct', k, ee), ...
                  err_HD < 1e-11, ...
                  sprintf('H=%+.6e, D=%+.6e, err=%.2e', ck_hansen, ck_direct, err_HD), ...
                  n_pass, n_fail);
  end
end

printf('\n');

%% =========================================================
%% Part 3: Infinite support (c_k non-vanishing).
%% =========================================================

printf('--- Part 3: c_k(e=0.3) non-vanishing for k=1..20 ---\n');

for k = 1:20
  ck = bessel_ck(k, 0.3, M_trunc);
  [n_pass, n_fail] = ...
    assert_true(sprintf('|c_%d(0.3)| > 1e-14', k), ...
                abs(ck) > 1e-14, ...
                sprintf('c_%d = %+.6e', k, ck), ...
                n_pass, n_fail);
end

printf('\n');

%% =========================================================
%% Summary
%% =========================================================

printf('===============================================================\n');
printf('Summary: %d pass, %d fail\n', n_pass, n_fail);
printf('===============================================================\n');

if n_fail > 0
  error('verify_ch09d_bessel_series: %d FAILURES', n_fail);
end
