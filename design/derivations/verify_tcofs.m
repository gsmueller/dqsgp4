% verify_tcofs.m
%
% Mechanical verifier for the SGP4 mean-longitude time-polynomial coefficients
% t2cof..t5cof (src/atmosphere/drag_coefficients.h:185-191). Standard 10:
% simplify(lhs - rhs) == 0. Resolves audit finding D-6 (the legacy sec.12 proof was
% an incomplete stream-of-consciousness; t3/t4/t5cof were never derived).
%
% Derivation. The drag-perturbed mean motion follows the decaying semi-major axis:
%   a(t) = a0 tempa^2  =>  n(t) = n0 (a/a0)^(-3/2) = n0 tempa^(-3) .
% The drag-induced extra mean longitude is l(t) = l0 + n0*templ with
%   templ = integral_0^t ( n(t')/n0 - 1 ) dt' = integral_0^t ( tempa(t')^(-3) - 1 ) dt'.
% The four t-cofs are the t^2, t^3, t^4, t^5 Taylor coefficients of templ
%   ( templ = t2cof t^2 + t3cof t^3 + t^4(t4cof + t*t5cof) , drag_coefficients.h:16 ).
% C1, D2, D3, D4 are carried as independent symbols (the code computes the t-cofs from
% the already-evaluated C1,D2,D3,D4 -- e.g. t3cof = D2 + 2 C1^2).

pkg load symbolic;
syms C1 D2 D3 D4 t tau real;

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

% tempa in the integration variable tau; n/n0 - 1 = tempa^-3 - 1
tempa_tau = 1 - C1*tau - D2*tau^2 - D3*tau^3 - D4*tau^4;
integrand = tempa_tau^(-3) - 1;
integrand_ser = taylor(integrand, tau, 0, 'order', 5);   % terms up to tau^4
templ = int(integrand_ser, tau, 0, t);                   % templ(t), up to t^5

% Extract the t^k coefficients of templ (k = 2..5)
t2cof_d = subs(diff(templ, t, 2), t, 0)/factorial(2);
t3cof_d = subs(diff(templ, t, 3), t, 0)/factorial(3);
t4cof_d = subs(diff(templ, t, 4), t, 0)/factorial(4);
t5cof_d = subs(diff(templ, t, 5), t, 0)/factorial(5);

% Code forms (drag_coefficients.h:185-191)
t2cof_c = sym(3)/2*C1;
t3cof_c = D2 + 2*C1^2;
t4cof_c = sym(1)/4*(3*D3 + C1*(12*D2 + 10*C1^2));
t5cof_c = sym(1)/5*(3*D4 + 12*C1*D3 + 6*D2^2 + 15*C1^2*(2*D2 + C1^2));

printf('=== t2cof..t5cof (mean-longitude drag polynomial; resolves D-6) ===\n\n');
[pass_count, fail_count, failed_names] = check_eq( ...
    't2cof = (3/2) C1                                  [drag_coefficients.h:185]', ...
    t2cof_d, t2cof_c, pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    't3cof = D2 + 2 C1^2                               [drag_coefficients.h:186]', ...
    t3cof_d, t3cof_c, pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    't4cof = (1/4)(3 D3 + C1(12 D2 + 10 C1^2))         [drag_coefficients.h:187-188]', ...
    t4cof_d, t4cof_c, pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    't5cof = (1/5)(3 D4 + 12 C1 D3 + 6 D2^2 + 15 C1^2(2 D2 + C1^2))  [drag_coefficients.h:189-191]', ...
    t5cof_d, t5cof_c, pass_count, fail_count, failed_names);

printf('\n=== t-cofs summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  t-cofs verification: FAIL\n'); exit(1);
else
    printf('\n  t-cofs verification: PASS  (Taylor coefs of integral(tempa^-3 - 1)dt)\n');
    exit(0);
end
