% verify_phase1.m
%
% Independent verifier for the Phase 1 Lane integrals I^(0,m) (Theorem 1.3.2),
% m = 1..6 — the p=0 family the SGP4 near-earth drag C-coefficients consume.
%
% Two INDEPENDENT routes (corpus independence requirement):
%   (R) analytic recurrence  I^(0,m+1) = I^(0,m) + (eta/m) dI^(0,m)/deta   (Thm 1.2.1),
%       from the residue base case I^(0,1) = (1-eta^2)^-1/2 (Thm 1.3.1).
%   (N) numerical quadrature  (1/2pi) int_0^2pi (1-eta cos f)^-m df, via periodic-
%       trapezoidal sampling (spectrally accurate for this analytic periodic
%       integrand) vs the closed form.
% plus the eta=0 sanity I^(0,m)(0) = 1.
%
% Phase 2.C decision (see sgp4_drag_phase1_lane_integrals.md §1.4): the p>=1 family is
% NOT built — C4/C5 reduce to p=0 (scoping in verify_C4.m / Phase 3), and per Standard
% 10 no theoretical result is accepted without a downstream code-match. Phase 1 closes
% at the p=0 family.

pkg load symbolic;
syms eta f real;

pass = 0; fail = 0; names = {};
function [p, fl, n] = ck(nm, lhs, rhs, p, fl, n)
    d = simplify(lhs - rhs);
    if isequal(d, sym(0))
        printf('  PASS: %s\n', nm); p = p + 1;
    else
        printf('  FAIL: %s   (lhs-rhs = %s)\n', nm, char(d)); fl = fl + 1; n{end+1} = nm;
    end
end

% Closed forms (Theorem 1.3.2.1 .. 1.3.2.6)
I0 = { (1-eta^2)^(-sym(1)/2), ...
       (1-eta^2)^(-sym(3)/2), ...
       (2 + eta^2)             / (2*(1-eta^2)^(sym(5)/2)), ...
       (2 + 3*eta^2)           / (2*(1-eta^2)^(sym(7)/2)), ...
       (8 + 24*eta^2 + 3*eta^4)/ (8*(1-eta^2)^(sym(9)/2)), ...
       (8 + 40*eta^2 + 15*eta^4)/(8*(1-eta^2)^(sym(11)/2)) };

printf('=== Phase 1: I^(0,m) verification ===\n\n');

printf('Route R -- recurrence I^(0,m+1) = I^(0,m) + (eta/m) dI^(0,m)/deta:\n');
[pass,fail,names] = ck('base case I^(0,1) = (1-eta^2)^-1/2', I0{1}, (1-eta^2)^(-sym(1)/2), pass,fail,names);
for m = 1:5
    rhs = I0{m} + (eta/m)*diff(I0{m}, eta);
    [pass,fail,names] = ck(sprintf('recurrence m=%d -> %d', m, m+1), I0{m+1}, rhs, pass,fail,names);
end

printf('\neta=0 sanity (integrand == 1 => integral == 1):\n');
for m = 1:6
    [pass,fail,names] = ck(sprintf('I^(0,%d)(0) = 1', m), subs(I0{m}, eta, 0), sym(1), pass,fail,names);
end

printf('\nRoute N -- numerical periodic-trapezoidal quadrature vs closed form:\n');
Npts = 20000; fs = (0:Npts-1) * (2*pi/Npts);
Ih = cell(1, 6);
for m = 1:6; Ih{m} = function_handle(I0{m}); end   % numeric handles (avoid float->sym)
reltol = 1e-9; allnum = true; worst = 0;
for etav = [0.1 0.3 0.5 0.7 0.85]
    for m = 1:6
        num = mean( (1 - etav*cos(fs)).^(-m) );
        clf = Ih{m}(etav);
        rel = abs(num - clf) / abs(clf);
        worst = max(worst, rel);
        if rel > reltol
            printf('  FAIL: I^(0,%d) eta=%.2f  num=%.12g closed=%.12g relerr=%.2e\n', ...
                   m, etav, num, clf, rel);
            allnum = false; fail = fail + 1; names{end+1} = sprintf('num m=%d eta=%.2f', m, etav);
        end
    end
end
if allnum
    printf('  PASS: all I^(0,m), m=1..6, eta in {.1 .3 .5 .7 .85} match (max relerr %.2e < %.0e)\n', worst, reltol);
    pass = pass + 1;
end

printf('\n=== Phase 1 summary ===\n');
tot = pass + fail;
printf('  %d / %d checks PASS\n', pass, tot);
if fail > 0
    printf('\n  FAILED:\n');
    for k = 1:numel(names); printf('    %s\n', names{k}); end
    printf('\n  Phase 1 verification: FAIL\n'); exit(1);
else
    printf('\n  Phase 1 verification: PASS  (recurrence + eta=0 sanity + numerical quadrature)\n');
    exit(0);
end
