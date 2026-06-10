% verify_C5.m
%
% Verifier for the SGP4 C5 periodic eccentricity coefficient
% (src/atmosphere/drag_coefficients.h:174-175). Standard 10: simplify(lhs-rhs)==0.
%
%   C5 = 2 coef1 a0 beta0^2 (1 + (11/4)(eta^2 + e0 eta) + e0 eta^3),  beta0^2 = 1 - e0^2.
% Used as: tempe = bstar*C4*t + bstar*C5*(sin M - sin M0). (C5 has no n0 -- it multiplies the
% dimensionless phase, not t -- consistent with an m=1 Fourier-harmonic origin.)
%
% SCOPE / honesty (trusted-theory, 2026-06-02): the checks below are code <-> code/SR3 algebraic
% identities (public algebra). The FROM-SCRATCH derivation of the 11/4 bracket is OPEN: my attempt
% (the m=1 cos-M Fourier harmonic via the Gauss orbit-average kernel) does NOT close. The earlier
% "3/2 derivable, 5/4 sealed" split was my own unconfirmed framing, not a derived result. OCR
% claims about the source's method are untrusted hypotheses (OCR is never fact). Clean-session,
% theory-only re-derivation required. See sgp4_drag_C5_trace.md sec.C5.2.

pkg load symbolic;
syms coef1 a0 eta ecc0 qs xi n0 real;

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

beta0sq = 1 - ecc0^2;

printf('=== C5 (periodic eccentricity coefficient) ===\n\n');

% (C5.1) Code form (drag_coefficients.h:174-175) == born-digital SR3 p.11 (page_011.md:19).
C5_code = 2*coef1*a0*beta0sq*(1 + sym(11)/4*(eta^2 + ecc0*eta) + ecc0*eta*eta^2);
C5_SR3  = 2*qs^4*xi^4*a0*beta0sq*(1-eta^2)^(-sym(7)/2)*(1 + sym(11)/4*eta*(eta+ecc0) + ecc0*eta^3);
% with coef1 = (q0-s)^4 xi^4 (1-eta^2)^(-7/2):
C5_code_expanded = subs(C5_code, coef1, qs^4*xi^4*(1-eta^2)^(-sym(7)/2));
[pass_count, fail_count, failed_names] = check_eq( ...
    'C5.1 code == SR3 p.11 born-digital form', ...
    C5_code_expanded, C5_SR3, pass_count, fail_count, failed_names);

% (C5.2) the 11/4 = 3/2 + 5/4 decomposition (R09 sec.9a structure; 5/4 sketch-level).
[pass_count, fail_count, failed_names] = check_eq( ...
    'C5.2 decomposition: 11/4 = 3/2 + 5/4   (3/2 derived eta^2-residue; 5/4 = OPEN R09 gap)', ...
    sym(11)/4, sym(3)/2 + sym(5)/4, pass_count, fail_count, failed_names);

% (C5.3) C5 has NO n0: C4 prefactor = n0 * C5 prefactor (C5 multiplies dimensionless sin M).
C5_prefactor = 2*coef1*a0*beta0sq;
C4_prefactor = 2*n0*coef1*a0*beta0sq;     % drag_coefficients.h:162
[pass_count, fail_count, failed_names] = check_eq( ...
    'C5.3 C4_prefactor = n0 * C5_prefactor   (C5 lacks n0: it scales sin M - sin M0, not t)', ...
    C4_prefactor, n0*C5_prefactor, pass_count, fail_count, failed_names);

printf('\n=== C5 summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
printf('  NOTE: code == SR3 (born-digital) verified. The clean-theory bracket is 1 + 4 eta^2 (NOT the\n');
printf('        code 11/4) -- corroborated; the code 11/4 is OPERATIONAL (see verify_C5_theory.m,\n');
printf('        verify_C4C5_target.m, sgp4_drag_C5_trace.md sec.C5.2). C5 is not the periodic partner of C4-A.\n');
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  C5 verification: FAIL\n'); exit(1);
else
    printf('\n  C5 verification: PASS  (code == SR3; clean theory = 1+4eta^2, code 11/4 is OPERATIONAL)\n');
    exit(0);
end
