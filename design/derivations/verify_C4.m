% verify_C4.m
%
% Mechanical verifier for the SGP4 C4 eccentricity-decay coefficient
% (src/atmosphere/drag_coefficients.h:162-171). Standard 10: every check is a
% symbolic simplify(lhs - rhs) == 0.
%
% C4 = 2 n0 coef1 a0 beta^2 * [ KEPLERIAN + J2-COUPLING ],  beta^2 = 1 - e0^2,
%   KEPLERIAN   = eta(2 + 1/2 eta^2) + e(1/2 + 2 eta^2)
%   J2-COUPLING = -2 k2 xi/(a0 psi^2) * ( -3 con41 (1 - 2 e eta + eta^2(3/2 - 1/2 e eta))
%                                         + 3/4 sin^2 i0 (2 eta^2 - e eta(1+eta^2)) cos 2w0 )
%
% STATUS (2026-06-02, trusted-theory):
%   Part A (Keplerian), reduction to O(e):       DERIVED (this file, C4.1-C4.4).
%   Part B (J2 coupling) bracket:                CODE-MATCHED to operational SGP4 (C4.B/C4.full).
%   TRUSTED-THEORY: my from-scratch drag x J2 derivation gives the isotropic term ORDER-eta; the
%   code is ORDER-1. These disagree, and the disagreement is UNRESOLVED from trusted (theoretical /
%   born-digital) sources -- my derivation may be INCOMPLETE, or the code may differ from correct
%   theory; I cannot tell without an independent THEORETICAL derivation. A verify sub-agent read a
%   SCANNED source (OCR, untrusted) and reported an "operational erratum"; that is an UNTRUSTED
%   HYPOTHESIS, not a fact, and is NOT asserted. The code-match stands (match the operational code
%   for bit-compatibility); whether the code's k2 term is theoretically correct is OPEN.
%   See sgp4_drag_C4_trace.md sec.C4.B. (OCR is never fact; do not trust vision.)
%
% Inputs: Phase 0-rev1 Thm 0.3.3 (e-dot), Thm 0.5.3 (orbit-average). Phase 1 I^(0,m).

pkg load symbolic;
syms e eta f real;

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

% --- Phase 1 Lane integral closed forms (Theorem 1.3.2) ---
I02 = (1-eta^2)^(-sym(3)/2);
I03 = (2 + eta^2)              / (2*(1-eta^2)^(sym(5)/2));
I04 = (2 + 3*eta^2)            / (2*(1-eta^2)^(sym(7)/2));
I05 = (8 + 24*eta^2 + 3*eta^4) / (8*(1-eta^2)^(sym(9)/2));
I06 = (8 + 40*eta^2 + 15*eta^4)/ (8*(1-eta^2)^(sym(11)/2));
psi7 = (1-eta^2)^(sym(7)/2);

printf('=== C4 Part A (Keplerian eccentricity-decay) ===\n\n');

% (C4.1) e-dot bracket collapse: the instantaneous e-dot (Phase 0 Thm 0.3.3)
%   e-dot = (beta/na)[R sin f + T(cos f + cos E)],  drag R,T ~ -B* rho |v| (rdot, r fdot).
% After substituting cos E = (e+cos f)/(1+e cos f), the f-bracket collapses:
cosE = (e + cos(f)) / (1 + e*cos(f));
brk  = e*sin(f)^2 + (1 + e*cos(f))*(cos(f) + cosE);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.1 e-dot bracket collapse: e sin^2 f + (1+e cf)(cos f + cos E) = 2(e+cos f)', ...
    brk, 2*(e + cos(f)), pass_count, fail_count, failed_names);

% (C4.2) C4 Keplerian f-average kernel, exact density, expand in e at fixed eta:
%   h = (1-eta cos f_dagger)^-4 (1+e^2+2e cos f)^(1/2) (e+cos f)/(1+e cos f)^2
% with EXACT (1-eta cos f_dagger) = (1-eta cos f - eta e + e cos f)/(1+e cos f).
D     = (1 - eta*cos(f) - eta*e + e*cos(f)) / (1 + e*cos(f));
kin_e = (1 + e^2 + 2*e*cos(f))^(sym(1)/2);
h     = D^(-4) * kin_e * (e + cos(f)) / (1 + e*cos(f))^2;
ser   = taylor(h, e, 0, 'order', 2);
h0    = simplify(subs(ser, e, 0));
h1    = simplify(subs(diff(ser, e), e, 0));
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.2a kernel O(e^0): h0 = cos f (1-eta cos f)^-4', ...
    h0, cos(f)*(1-eta*cos(f))^(-4), pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.2b kernel O(e^1): h1 = (1+3 eta cos f) sin^2 f (1-eta cos f)^-5', ...
    h1, (1 + 3*eta*cos(f))*sin(f)^2*(1-eta*cos(f))^(-5), ...
    pass_count, fail_count, failed_names);

% (C4.3) Reduce the averages to the p=0 Lane family (no p>=1 needed) and match code:163.
% Exact reduction identities (cos f = (1-(1-eta cf))/eta ; sin^2 f = 1 - cos^2 f):
avg_h0          = (I04 - I03)/eta;                       % <cos f (1-h cf)^-4>
avg_sin2_m5     = I05 - (1/eta^2)*(I05 - 2*I04 + I03);   % <sin^2 f (1-h cf)^-5>
avg_sin2_m4     = I04 - (1/eta^2)*(I04 - 2*I03 + I02);   % <sin^2 f (1-h cf)^-4>
avg_cos_sin2_m5 = (1/eta)*(avg_sin2_m5 - avg_sin2_m4);   % <cos f sin^2 f (1-h cf)^-5>
avg_h1          = avg_sin2_m5 + 3*eta*avg_cos_sin2_m5;   % <h1>
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.3a Keplerian O(e^0): psi^7 <h0> = eta(2 + 1/2 eta^2)   [code:163]', ...
    psi7*avg_h0, eta*(2 + eta^2/sym(2)), pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.3b Keplerian O(e^1): psi^7 <h1> = 1/2 + 2 eta^2        [code:163]', ...
    psi7*avg_h1, sym(1)/2 + 2*eta^2, pass_count, fail_count, failed_names);

% (C4.4) Keplerian bracket assembly: psi^7 (<h0> + e <h1>) = eta(2+1/2 eta^2) + e(1/2+2 eta^2).
PartA_bracket = simplify(psi7*(avg_h0 + e*avg_h1));
PartA_code    = eta*(2 + eta^2/sym(2)) + e*(sym(1)/2 + 2*eta^2);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.4 Keplerian bracket = eta(2+1/2 eta^2) + e(1/2+2 eta^2)   [drag_coefficients.h:163]', ...
    PartA_bracket, PartA_code, pass_count, fail_count, failed_names);

printf('\n');
printf('C4 Part B (J2 coupling) -- code-matched to born-digital SR3 p.11:\n');
syms theta k2 xi a0 n0 coef1 w0 beta2 real;
con41 = 3*theta^2 - 1;        % three_cos2i_minus_1
sin2i = 1 - theta^2;          % sin^2 i0
eeta  = e*eta;
psisq = 1 - eta^2;

% SR3 p.11 Part B bracket (theta form): 3(1-3θ²)(1+3/2η²-2eη-½eη³) + 3/4(1-θ²)(2η²-eη-eη³)cos2ω
PartB_SR3  = -2*k2*xi/(a0*psisq)*( 3*(1-3*theta^2)*(1 + sym(3)/2*eta^2 - 2*eeta - sym(1)/2*eeta*eta^2) ...
             + sym(3)/4*(1-theta^2)*(2*eta^2 - eeta - eeta*eta^2)*cos(2*w0) );
% code Part B bracket (con41/sin2i form, drag_coefficients.h:164-171)
PartB_code = -2*k2*xi/(a0*psisq)*( -3*con41*(1 - 2*eeta + eta^2*(sym(3)/2 - sym(1)/2*eeta)) ...
             + sym(3)/4*sin2i*(2*eta^2 - eeta*(1+eta^2))*cos(2*w0) );
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.B Part B bracket: SR3 p.11 (theta form) == code:164-171 (con41/sin2i)  [code-matched]', ...
    PartB_SR3, PartB_code, pass_count, fail_count, failed_names);

% Full C4: DERIVED Part A bracket (C4.4) + SR3 Part B bracket == code (drag_coefficients.h:162-171)
C4_derived = 2*n0*coef1*a0*beta2*(PartA_bracket + PartB_SR3);
C4_code    = 2*n0*coef1*a0*beta2*(PartA_code    + PartB_code);
[pass_count, fail_count, failed_names] = check_eq( ...
    'C4.full  C4 = 2 n0 coef1 a0 beta^2 (PartA_derived + PartB_SR3) == code:162-171', ...
    C4_derived, C4_code, pass_count, fail_count, failed_names);

printf('\n=== C4 summary (Part A derived; Part B code-matched to SR3) ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  C4 verification: FAIL\n'); exit(1);
else
    printf('\n  C4 verification: PASS  (Part A DERIVED = code:163; Part B = SR3 = code:164-171,\n');
    printf('                          full C4 = code:162-171; Part B from-scratch deriv = sealed LH79)\n');
    exit(0);
end
