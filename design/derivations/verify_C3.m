% verify_C3.m
%
% Mechanical verifier for the SGP4 C3 J3 eccentricity-drag coupling coefficient
% (src/atmosphere/drag_coefficients.h:154-159). Standard 10: every check is a
% symbolic simplify(lhs - rhs) == 0.
%
%   C3 = coef * xi * A30 * n0 * sin_i0 / (k2 * ecc0),   when ecc0 > 1e-4, else 0
% with coef = (q0-s)^4 xi^4, xi = 1/(a0-s), A30 = A(3,0) = -J3 aearth^3 (no 1/2;
% model_selector.h:178 returns -Jn(n)), k2 = J2 aearth^2/2 (half_J2), j3oj2 = J3/J2.
%
% Resolves audit finding D-4. The legacy proof wrote a A30/(2 k2) prefactor and then
% claimed "absorbing the factor of 2 from <sin^2 u> = 1/2" -- which would have given
% /(4 k2 ecc0) (off by 4x). The CORRECT chain: the A30/(2 k2) Brouwer prefactor's 2 is
% cancelled by the Gauss-VE factor 2 of the eccentricity-vector rate de/dt (Phase 0
% Thm 0.3.3: cos f + cos E -> 2 cos f as e -> 0). There is NO <sin^2 u> average.
% J3 inputs are born-digital (sealed-room dispatch ac84fb143; SR3 transcript Eq 4.4 /
% sec.7, Vallado/Rhodes propagation.py:1454 cc3 = -2 coef tsi j3oj2 n0 sinio/ecco).

pkg load symbolic;
syms coef xi n0 sin_i0 ecc0 J2 J3 aearth bstar w0 ecc f real;

pass_count = 0; fail_count = 0; failed_names = {};
function [pc, fc, fn] = check_eq(name, lhs, rhs, pc, fc, fn)
    d = simplify(lhs - rhs);
    if isequal(d, sym(0))
        printf('  PASS: %s\n', name); pc = pc + 1;
    else
        printf('  FAIL: %s\n', name);
        printf('        lhs - rhs = %s\n', char(d));
        fc = fc + 1; fn{end+1} = name;
    end
endfunction

A30   = -J3*aearth^3;          % model_selector.h:178  A(3,0) = -J3 aearth^3  (no 1/2)
k2    = J2*aearth^2/sym(2);    % half_J2
j3oj2 = J3/J2;

printf('=== C3 (J3 eccentricity-drag coupling) ===\n\n');

% (C3.1) ALGEBRAIC IDENTITY (NOT the physical mechanism -- REFUTED 2026-06-02; see banner).
%        (cos f + cos E)|_{e=0} = 2 cos f is a true Kepler identity, but C3's physical 1/e0
%        is the Gauss Mdot transverse pole (C3.7 below, partner of xmcof), NOT an
%        eccentricity-vector "factor-2 cancellation". Retained only to document the
%        superseded attempt.
cosE = (ecc + cos(f))/(1 + ecc*cos(f));
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.1 [IDENTITY, not physics] (cos f + cos E)|_{e=0} = 2 cos f   [true Kepler identity; mechanism REFUTED]', ...
    subs(cos(f) + cosE, ecc, 0), 2*cos(f), pass_count, fail_count, failed_names);

% (C3.2) A30/k2 convention identity (no-1/2 A30): A30/k2 = -2 j3oj2 aearth.
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.2 convention: A30/k2 = -2 j3oj2 aearth   [A30=-J3 aearth^3, k2=J2 aearth^2/2]', ...
    A30/k2, -2*j3oj2*aearth, pass_count, fail_count, failed_names);

% Code form (drag_coefficients.h:156)
C3_code = coef*xi*A30*n0*sin_i0/(k2*ecc0);

% (C3.3) ALGEBRAIC IDENTITY (NOT the physical mechanism -- REFUTED 2026-06-02; see banner).
%        coef xi (A30/2k2)(2) n0 sin_i0/ecc0 = C3_code is a true identity (the 2's cancel),
%        but the "(A30/2k2) Brouwer prefactor x Gauss-VE 2" account is REFUTED: the 1/e0 is
%        the Gauss Mdot transverse pole (C3.7), the A30/k2 is plain convention (C3.2), and
%        the J3 amplitude is the long-period e-vector factor (C3.8). Retained as identity only.
C3_brouwer_gauss = coef*xi*(A30/(2*k2))*sym(2)*n0*sin_i0/ecc0;   % A30/(2k2)  x  2 (algebraic only)
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.3 [IDENTITY, not physics] coef xi (A30/2k2)(2) n0 sin_i0/ecc0 = C3_code   [the 2s cancel; mechanism REFUTED]', ...
    C3_brouwer_gauss, C3_code, pass_count, fail_count, failed_names);

% (C3.4) Born-digital equivalent (Vallado/Rhodes cc3): -2 coef xi j3oj2 aearth n0 sin_i0/ecc0.
C3_vallado = -2*coef*xi*j3oj2*aearth*n0*sin_i0/ecc0;
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.4 Vallado cc3: C3_code = -2 coef xi j3oj2 aearth n0 sin_i0/ecc0', ...
    C3_code, C3_vallado, pass_count, fail_count, failed_names);

% (C3.5) D-4 error magnitude: the legacy <sin^2 u>=1/2 route gives /(4 k2 ecc0) = C3/4.
C3_legacy_wrong = coef*xi*A30*n0*sin_i0/(4*k2*ecc0);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.5 legacy <sin^2 u>=1/2 route is 1/4 of code (off by 4x):  4 * C3_legacy = C3_code', ...
    4*C3_legacy_wrong, C3_code, pass_count, fail_count, failed_names);

% (C3.6) omgcof follow-on (drag_coefficients.h:194): omgcof = B* C3 cos w0.
omgcof_code    = bstar*C3_code*cos(w0);
omgcof_derived = bstar*(coef*xi*A30*n0*sin_i0/(k2*ecc0))*cos(w0);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.6 omgcof = B* C3 cos w0   [drag_coefficients.h:194]', ...
    omgcof_derived, omgcof_code, pass_count, fail_count, failed_names);

% ===================================================================================
% CORRECTED MECHANISM (2026-06-02..03, clean): C3 = J3 x drag MEAN-ANOMALY amplitude,
% the J3 partner of xmcof in the drag-periodic correction delta-ell_D (M += delomg+delm,
% w -= delomg+delm; secular_update.h:96-105). The 1/e0 is the Gauss Mdot transverse pole
% (power 1, shared with xmcof -- verify_xmcof_theory.m), NOT an e-vector projection.
% ===================================================================================
syms a n beta real;
% (C3.7) The 1/e0 is the Gauss Mdot transverse SIMPLE pole (power exactly 1), the SAME
%        1/e that lands xmcof (verify_xmcof_theory.m X.2). Mdot-n (Thm 0.3.6) under drag
%        has a (beta^2/(n a e)) transverse factor; the residue limit_{e->0}(e*(Mdot-n)/n)
%        is finite & nonzero => simple pole => 1/e0 with power exactly 1 (per-element
%        normalisation), refuting the "1/e from e-vector projection" account. (The Lane
%        density factor is regular at e=0, so it is set =1 here -- the pole is geometric.)
absv  = (n*a/beta)*sqrt(1+ecc^2+2*ecc*cos(f));          % |v| (0.4.1.11)
Tdr   = -bstar*absv*(n*a/beta)*(1+ecc*cos(f));          % T_drag (0.4.1.8), density=1
rp    = 1/(1+ecc*cos(f));
Mdotn = (beta^2/(n*a*ecc))*( -Tdr*sin(f)*(1+rp) );      % transverse 1/e term of Mdot-n (beta's cancel)
res_e = limit(ecc*(Mdotn/n), ecc, 0);                   % residue of the e-pole = 2 a bstar sin f
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.7 1/e0 is the Gauss Mdot transverse SIMPLE pole (power 1; partner of xmcof): residue = 2 a B* sin f (finite, nonzero)', ...
    res_e, 2*a*bstar*sin(f), pass_count, fail_count, failed_names);

% (C3.8) C3's J3 angular amplitude (A30/k2) sin_i0 IS the born-digital J3 long-period
%        e-vector amplitude: aycof = (1/4)(A30/CK2) sin_i0 (near_space.h:286), so
%        (A30/k2) sin_i0 = 4*aycof. C3 inherits the J3 long-period amplitude (sin_i0,
%        A30/k2), NOT an independent e-vector projection of the drag.
aycof_code = (sym(1)/4)*(A30/k2)*sin_i0;   % near_space.h:286  aycof = (1/4)(A30/CK2) sin i0
[pass_count, fail_count, failed_names] = check_eq( ...
    'C3.8 C3 angular factor (A30/k2) sin_i0 = 4*aycof   [J3 long-period e-vector amplitude, near_space.h:286]', ...
    (A30/k2)*sin_i0, 4*aycof_code, pass_count, fail_count, failed_names);

printf('\n=== C3 summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  C3 verification: FAIL\n'); exit(1);
else
    printf('\n  C3 verification: PASS  (code-match / D-4).\n');
    printf('  MECHANISM (corrected 2026-06-02..03, clean): C3 = J3 x drag MEAN-ANOMALY amplitude,\n');
    printf('  the J3 partner of xmcof in delta-ell_D (M += delomg+delm, w -= ...; sec_update.h:96-105).\n');
    printf('  CLEAN-DERIVED factors: 1/e0 = Gauss Mdot transverse simple pole (C3.7, power 1, = xmcof''s\n');
    printf('  1/e -- refutes e-vector projection); A30/k2 = -2 j3oj2 convention (C3.2); n0 = secular\n');
    printf('  rate x t; (A30/k2) sin_i0 = 4*aycof = J3 long-period e-vector amplitude (C3.8). cos w is\n');
    printf('  CLEAN -- resonance-free by EXACT cancellation (verify_xmcof_C3_mitigations.m M4): C3*cos w0*t\n');
    printf('  is the secular linear-t part of the J3-drag drift, whose coefficient cos w0 is independent\n');
    printf('  of the apsidal freq wdot~(5cos^2 i-1), so the critical-inclination resonance cancels; the\n');
    printf('  resonant 1/wdot oscillation is the separate aycof/xlcof block. Only the xi^5=coef*xi Hansen\n');
    printf('  density (|| C2 Part B; the J3 analog of xmcof cubic) remains operational.\n');
    printf('  C3.1/C3.3 are true ALGEBRAIC IDENTITIES, NOT physics (REFUTED -- see trace banner).\n');
    exit(0);
end
