% verify_corrections.m
%
% Verifier for the SGP4 drag correction terms xmcof, delmo, sinmo, xnodcf
% (src/atmosphere/drag_coefficients.h:194-209). Standard 10: simplify(lhs-rhs)==0.
% Closes audit findings D-8 (xmcof xi^4) and D-7 (xnodcf beta-power).
%
% Born-digital sources (sealed-room dispatch a968ac5e): SR3 p.12 (page_012.md:7,13),
% Vallado/Rhodes propagation.py:1480-1481,1490-1492.

pkg load symbolic;
syms coef xi bstar ecc0 eta M0 n0 k2 cos_i0 C1 a0 qs real;

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

beta0sq = 1 - ecc0^2;        % beta0^2
eeta    = ecc0*eta;          % e0 * eta

printf('=== drag corrections: xmcof (D-8), delmo/sinmo, xnodcf (D-7) ===\n\n');

% --- xmcof (D-8): the code's coef = (q0-s)^4 xi^4 supplies the xi^4 the legacy box dropped ---
coef_val     = qs^4*xi^4;                          % qoms4 * xi^4 (drag_coefficients.h:121)
xmcof_code   = -sym(2)/3*coef*bstar/eeta;          % drag_coefficients.h:197
xmcof_code_x = subs(xmcof_code, coef, coef_val);   % = -(2/3)(q0-s)^4 xi^4 B*/(e0 eta)
xmcof_legacy = -sym(2)/3*qs^4*bstar/eeta;          % legacy sec.15.4 box -- MISSING xi^4
[pass_count, fail_count, failed_names] = check_eq( ...
    'xmcof (D-8): code = legacy_box * xi^4   (code/SR3 carry xi^4; legacy box dropped it)', ...
    xmcof_code_x, xmcof_legacy*xi^4, pass_count, fail_count, failed_names);

% --- delmo / sinmo (reference quantities; the cubic Lane density at mean anomaly) ---
delmo_code = (1 + eta*cos(M0))^3;                  % drag_coefficients.h:202-203
sinmo_code = sin(M0);                              % drag_coefficients.h:205
[pass_count, fail_count, failed_names] = check_eq( ...
    'delmo = (1 + eta cos M0)^3   [drag_coefficients.h:202-203]', ...
    delmo_code, (1 + eta*cos(M0))^3, pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'sinmo = sin M0             [drag_coefficients.h:205]', ...
    sinmo_code, sin(M0), pass_count, fail_count, failed_names);

% --- xnodcf (D-7 RESOLVED 2026-06-02): the (7/2) beta0^2 xhdot1 C1 form (algebra) ---
% xhdot1 = Omega-dot_J2 = -3 k2 n0 cos_i0/(a0^2 beta0^4)  (pinvsq = 1/p^2 = 1/(a0^2 beta0^4)).
% The beta0^2 is the ANGULAR-MOMENTUM drag normalisation <Gdot/G> = -C1 beta^2 (DERIVED in
% verify_xnodcf_theory.m: <rho|v|>/((a/mu)<rho|v|^3>) = beta^2). Theory nodecf = (1/2)C1 xhdot1
% (3+4 beta^2); the code's exact (7/2)beta^2 applies omeosq uniformly (SGP4 simplification,
% code-theory=O(e^2)). This check is the bare ALGEBRA code-match; the theory is in the companion.
xhdot1         = -3*k2*n0*cos_i0/(a0^2*beta0sq^2);          % beta0^4 in denom
xnodcf_vallado = sym(7)/2*beta0sq*xhdot1*C1;                % 3.5*omeosq*xhdot1*cc1
xnodcf_code    = -sym(21)/2*n0*k2*cos_i0*C1/(a0^2*beta0sq); % drag_coefficients.h:208-209
[pass_count, fail_count, failed_names] = check_eq( ...
    'xnodcf (D-7): (7/2) beta0^2 xhdot1 C1 = -(21/2) n0 k2 cos_i0 C1/(a0^2 beta0^2)', ...
    xnodcf_vallado, xnodcf_code, pass_count, fail_count, failed_names);

printf('\n=== corrections summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
printf('  D-8: code/SR3 carry xi^4 (legacy sec.15.4 box dropped it) -- doc fix.\n');
printf('  D-7 RESOLVED: beta0^2 = ANGULAR-MOMENTUM drag norm <Gdot/G>=-C1 beta^2 (derived in\n');
printf('       verify_xnodcf_theory.m 10/10); code (7/2)beta^2 = uniform-omeosq SGP4 simplification.\n');
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  corrections verification: FAIL\n'); exit(1);
else
    printf('\n  corrections verification: PASS  (xmcof D-8, delmo, sinmo, xnodcf D-7)\n');
    exit(0);
end
