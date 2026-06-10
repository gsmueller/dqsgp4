% verify_D3D4.m
%
% Mechanical verifier for the SGP4 D3, D4 higher-order secular-drag coefficients
% (src/atmosphere/drag_coefficients.h:179-182). Standard 10: simplify(lhs-rhs)==0.
%
% Same framework as D2 (verify_D2.m / sgp4_drag_secular_taylor_trace.md): the t^3, t^4
% Taylor coefficients of a(t) = a0 tempa^2 under the decay ODE  a-dot = -2 a C1(a),
% C1(a) ~ (a-s)^-4 a^-1/2 (xi^4 from coef1, a^-1/2 from n~a^-3/2 * explicit a).
%
% Targets (Q-9/Q-10 confirm the temp_d grouping):
%   temp_d = D2 xi C1 / 3
%   D3 = (17 a0 + s) * temp_d                                   [drag_coefficients.h:179-180]
%   D4 = (1/2) * temp_d * a0 * xi * (221 a0 + 31 s) * C1        [drag_coefficients.h:181-182]
% The (17 a0 + s) / (221 a0 + 31 s) arise from the higher a-derivatives of (a-s)^-4 a^-1/2.

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

% --- decay ODE  a' = g(a),  g = -2 a C1(a),  C1(a) ~ (a-s)^-4 a^-1/2 ---
C1_expr = cc*(A - s)^(-4) * A^(-sym(1)/2);
g       = -2*A*C1_expr;
g1 = diff(g, A); g2 = diff(g, A, 2); g3 = diff(g, A, 3);
% evaluate at a0
G  = subs(g , A, a0);
G1 = subs(g1, A, a0);
G2 = subs(g2, A, a0);
G3 = subs(g3, A, a0);
C1 = subs(C1_expr, A, a0);

% a-derivatives at t=0 via the chain rule (a'=g):
%   a'   = g
%   a''  = g' g
%   a''' = g'' g^2 + g'^2 g
%   a'''' = g''' g^3 + 4 g'' g' g^2 + g'^3 g
ad1 = G;
ad2 = G1*G;
ad3 = G2*G^2 + G1^2*G;
ad4 = G3*G^3 + 4*G2*G1*G^2 + G1^3*G;

% a(t)/a0 Taylor coefficients  (a = a0 + ad1 t + ad2/2 t^2 + ad3/6 t^3 + ad4/24 t^4)
b2 = ad2/(2*a0);
b3 = ad3/(6*a0);
b4 = ad4/(24*a0);

% tempa^2 match (D.0.3):  D2=(C1^2-b2)/2 ; D3=C1 D2 - b3/2 ; D4=(D2^2+2 C1 D3 - b4)/2
D2 = simplify((C1^2 - b2)/2);
D3 = simplify(C1*D2 - b3/2);
D4 = simplify((D2^2 + 2*C1*D3 - b4)/2);

printf('=== D3, D4 (Taylor of a(t), t^3/t^4) ===\n\n');

% sanity: D2 still = 4 a0 xi C1^2
[pass_count, fail_count, failed_names] = check_eq( ...
    'D.sanity  D2 = 4 a0 xi C1^2', D2, 4*a0*xi*C1^2, pass_count, fail_count, failed_names);

temp_d = D2*xi*C1/3;     % drag_coefficients.h:179

% (D3) D3 = (17 a0 + s) temp_d
D3_code = (17*a0 + s)*temp_d;
[pass_count, fail_count, failed_names] = check_eq( ...
    'D3  = (17 a0 + s) D2 xi C1 / 3        [drag_coefficients.h:179-180]', ...
    D3, D3_code, pass_count, fail_count, failed_names);

% (D4) D4 = (1/2) temp_d a0 xi (221 a0 + 31 s) C1
D4_code = (sym(1)/2)*temp_d*a0*xi*(221*a0 + 31*s)*C1;
[pass_count, fail_count, failed_names] = check_eq( ...
    'D4  = (1/2) temp_d a0 xi (221 a0 + 31 s) C1   [drag_coefficients.h:181-182]', ...
    D4, D4_code, pass_count, fail_count, failed_names);

printf('\n=== D3/D4 summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  D3/D4 verification: FAIL\n'); exit(1);
else
    printf('\n  D3/D4 verification: PASS  (17a0+s and 221a0+31s from C1(a) a-derivatives)\n');
    exit(0);
end
