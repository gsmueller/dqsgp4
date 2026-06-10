% verify_D2.m
%
% Mechanical verifier for the SGP4 D2 higher-order secular-drag coefficient
% (src/atmosphere/drag_coefficients.h:178). Standard 10: simplify(lhs - rhs) == 0.
%
%   a(t) = a0 * tempa^2,  tempa = 1 - C1 t - D2 t^2 - ...   (drag_coefficients.h:14)
%   D2 = 4 a0 xi C1^2                                       (drag_coefficients.h:178)
%
% Resolves audit finding D-5. The audit's re-derivation produced a spurious residual
% D2 = -(1/2) C1^2 + 4 a0 xi C1^2 (or -(5/4)C1^2 ...). ROOT CAUSE (this file): it dropped
% the a-dependence of the mean motion n ~ a^-3/2 (and the explicit a in C2). The true
% drag rate is a-dot = -2 a C1(a) with C1(a) = B* C2(a) ~ xi^4 a^-1/2 (leading, eta->0:
%   C2 ~ coef1 n a, coef1 ~ xi^4, n ~ a^-3/2, explicit a  =>  C1 ~ (a-s)^-4 a^-1/2).
% Then d(ln C1)/da = -4 xi - 1/(2a); the -1/(2a) makes the stray C1^2 cancel against
% tempa^2, giving D2 = 4 a0 xi C1^2 EXACTLY. Dropping n(a) (using C1 ~ (a-s)^-4 only)
% reproduces the audit's -(1/2)C1^2 residual -- check D2.3 below demonstrates this.
% The eta/psi a-dependence is a higher-order (AFGP4->SGP4) truncation, not code-used.

pkg load symbolic;
syms A a0 s cc t real;

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

xi = 1/(a0 - s);

% --- Taylor of the decay ODE  a-dot = -2 a C1(a)  to t^2, matched to a0 tempa^2 ---
% Helper: given C1(A), return D2 from the t^2 Taylor match.
%   a-dot0   = -2 a0 C1(a0)
%   a-ddot0  = -2 a-dot0 ( C1(a0) + a0 C1'(a0) )
%   a(t)/a0  = 1 - 2 C1 t + (a-ddot0/(2 a0)) t^2
%   tempa^2  = 1 - 2 C1 t + (C1^2 - 2 D2) t^2   =>   D2 = (C1^2 - a-ddot0/(2 a0))/2
function D2 = d2_from_C1(C1_expr, A, a0, s)
    C1      = subs(C1_expr, A, a0);
    C1p     = subs(diff(C1_expr, A), A, a0);
    adot0   = -2*a0*C1;
    addot0  = -2*adot0*(C1 + a0*C1p);
    D2      = simplify((C1^2 - addot0/(2*a0))/2);
end

printf('=== D2 (Taylor of a(t); resolves D-5) ===\n\n');

% (D2.1) FULL C1(a) ~ (a-s)^-4 a^-1/2  (xi^4 from coef1, a^-1/2 from n~a^-3/2 * explicit a)
C1_full   = cc*(A - s)^(-4) * A^(-sym(1)/2);
C1_at_a0  = subs(C1_full, A, a0);
D2_full   = d2_from_C1(C1_full, A, a0, s);
[pass_count, fail_count, failed_names] = check_eq( ...
    'D2.1 with n(a)~a^-3/2:  D2 = 4 a0 xi C1^2   [drag_coefficients.h:178, EXACT]', ...
    D2_full, 4*a0*xi*C1_at_a0^2, pass_count, fail_count, failed_names);

% (D2.2) the C1^2 cancellation is exact: D2_full has NO residual C1^2 term beyond 4 a0 xi C1^2
[pass_count, fail_count, failed_names] = check_eq( ...
    'D2.2 no spurious residual:  D2 - 4 a0 xi C1^2 = 0', ...
    D2_full - 4*a0*xi*C1_at_a0^2, sym(0), pass_count, fail_count, failed_names);

% (D2.3) AUDIT'S ERROR reproduced: dropping n(a) (C1 ~ (a-s)^-4 only) gives the -(1/2)C1^2 residual
C1_no_n   = cc*(A - s)^(-4);
C1n_at_a0 = subs(C1_no_n, A, a0);
D2_no_n   = d2_from_C1(C1_no_n, A, a0, s);
[pass_count, fail_count, failed_names] = check_eq( ...
    'D2.3 dropping n(a): D2 = -(1/2)C1^2 + 4 a0 xi C1^2  [the audit D-5 residual, an omission]', ...
    D2_no_n, -C1n_at_a0^2/2 + 4*a0*xi*C1n_at_a0^2, pass_count, fail_count, failed_names);

printf('\n=== D2 summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  D2 verification: FAIL\n'); exit(1);
else
    printf('\n  D2 verification: PASS  (D2 = 4 a0 xi C1^2 exact; D-5 residual = dropped n(a))\n');
    exit(0);
end
