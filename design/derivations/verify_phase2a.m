% verify_phase2a.m
%
% Mechanical verifier for Phase 2.A — symbolic trace of the SGP4 C2 drag
% coefficient (drag_coefficients.h:146-149). Standard 10: every check is a
% symbolic simplify(lhs - rhs) == 0 (numerical-at-points is NOT proof).
%
% STATUS (2026-06-01):
%   Part A (Keplerian drag), reduction to O(e):                 COMPLETE.
%   Part B (J2 density coupling), leading e:                    COMPLETE
%     (J2-radial secular multiplier from sealed-room sub-agent dispatch;
%      born-digital: Vallado SGP4.cpp:1990-2002, SR3 Eq 6.7.)
%   Full C2 == drag_coefficients.h:146-149 (check A.8):         COMPLETE.
%   C1 = B*.C2 == drag_coefficients.h:152 (checks A.9.2/A.9.3): COMPLETE.
%   13/13 checks PASS.
%
% Inputs consumed:
%   Phase 0-rev1: Thm 0.4.2 (a_dot integrand), Thm 0.5.3 (orbit-average setup).
%   Phase 1:      Thm 1.3.2 Lane integrals I^(0,m).
% Target: src/atmosphere/drag_coefficients.h:146-149.

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
I03  = (2 + eta^2)            / (2*(1-eta^2)^(sym(5)/2));
I04  = (2 + 3*eta^2)          / (2*(1-eta^2)^(sym(7)/2));
I05  = (8 + 24*eta^2 + 3*eta^4) / (8*(1-eta^2)^(sym(9)/2));
psi7 = (1-eta^2)^(sym(7)/2);    % psi^7 = (1-eta^2)^(7/2)

printf('=== Phase 2.A verification ===\n\n');
printf('Part A (Keplerian drag) reduction to O(e):\n');

% (A.4.1) kinematic factor g = (1+e^2+2e cos f)^(3/2)/(1+e cos f)^2 = 1 + e cos f + O(e^2)
kin = (1 + e^2 + 2*e*cos(f))^(sym(3)/2) / (1 + e*cos(f))^2;
kin_ser = taylor(kin, e, 0, 'order', 2);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.4.1 kinematic O(e^0) = 1', simplify(subs(kin_ser, e, 0)), sym(1), ...
    pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.4.1 kinematic O(e^1) = cos f', simplify(subs(diff(kin_ser, e), e, 0)), cos(f), ...
    pass_count, fail_count, failed_names);

% (A.3.1) EXACT (1 - eta cos f_dagger) as a function of f, from exact Kepler
% r = a beta^2/(1+e cos f) and Lane def r-s = (a-s)(1-eta cos f_dagger), with
% eta = a e/(a-s) used to eliminate (a,s):  (held at fixed eta, expanded in e).
D = (1 - eta*cos(f) - eta*e + e*cos(f)) / (1 + e*cos(f));
integ = kin * D^(-4);
ser = taylor(integ, e, 0, 'order', 2);
a0 = simplify(subs(ser, e, 0));
a1 = simplify(subs(diff(ser, e), e, 0));
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.4 exact-density O(e^0): a0 = (1-eta cos f)^-4', ...
    a0, (1-eta*cos(f))^(-4), pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.4.3 exact-density O(e^1): a1 = cos f (1-h cf)^-4 + 4h sin^2 f (1-h cf)^-5', ...
    a1, cos(f)*(1-eta*cos(f))^(-4) + 4*eta*sin(f)^2*(1-eta*cos(f))^(-5), ...
    pass_count, fail_count, failed_names);

% cos-power and sin^2 reduction identities (exact algebra behind the Lane averaging)
[pass_count, fail_count, failed_names] = check_eq( ...
    'cos-power reduction: cos f (1-h cf)^-4 = (1/h)[(1-h cf)^-4 - (1-h cf)^-3]', ...
    cos(f)*(1-eta*cos(f))^(-4), ...
    (1/eta)*((1-eta*cos(f))^(-4) - (1-eta*cos(f))^(-3)), ...
    pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'sin^2 reduction: sin^2 f (1-h cf)^-5 = (1-h cf)^-5 - (1/h^2)[(1-h cf)^-5 -2(1-h cf)^-4 +(1-h cf)^-3]', ...
    sin(f)^2*(1-eta*cos(f))^(-5), ...
    (1-eta*cos(f))^(-5) - (1/eta^2)*((1-eta*cos(f))^(-5) - 2*(1-eta*cos(f))^(-4) + (1-eta*cos(f))^(-3)), ...
    pass_count, fail_count, failed_names);

% (A.4.2)-(A.4.3) Part-A averages vs SGP4 code Part-A polynomial coefficients.
% <a0> = <(1-h cf)^-4> = I04.
% <a1> = <cos f (1-h cf)^-4> + 4h <sin^2 f (1-h cf)^-5>
%       = (I04-I03)/h + 4h[ I05 - (1/h^2)(I05 - 2 I04 + I03) ].
avg_a0 = I04;
avg_a1 = (I04 - I03)/eta + 4*eta*I05 - (4/eta)*(I05 - 2*I04 + I03);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.4.2 Part-A const: <a0> = psi^-7 (1 + 3/2 eta^2)   [closes D-9]', ...
    avg_a0, (1 + sym(3)/2*eta^2)/psi7, pass_count, fail_count, failed_names);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.4.3 Part-A e-term: <a1> = psi^-7 eta(4 + eta^2)', ...
    avg_a1, eta*(4 + eta^2)/psi7, pass_count, fail_count, failed_names);

printf('\n');
printf('Part B (J2 density coupling) at leading e:\n');

syms k2 p a a0 xi beta con41 n0 coef real;
I05 = (8 + 24*eta^2 + 3*eta^4) / (8*(1-eta^2)^(sym(9)/2));
% Sealed-room input (sub-agent; born-digital: Vallado SGP4.cpp:1990-2002, SR3 Eq 6.7):
%   J2 secular radial multiplier  dr/r = -(3/2)(k2/p^2) beta (3cos^2 i - 1).
dr_r = -(sym(3)/2)*(k2/p^2)*beta*con41;
% Density gradient: rho(r_kep + dr) = rho(r_kep)(1 - 4 dr/(r-s)),
% 1/(r-s) = xi/(1-eta cos f); r/(r-s) = a xi/(1-eta cos f) at leading e.
r_rs = a*xi/(1-eta*cos(f));
PBint = -4*dr_r*r_rs*(1-eta*cos(f))^(-4);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.6 density-gradient power raise: -4(dr/r)(r/(r-s))(1-h cf)^-4 = 6(k2/p^2)b con41 a xi (1-h cf)^-5', ...
    PBint, 6*(k2/p^2)*beta*con41*a*xi*(1-eta*cos(f))^(-5), pass_count, fail_count, failed_names);
% Orbit average <(1-h cf)^-5> = I05 (Phase 1 Thm 1.3.2.5); assemble Part-B C2 term.
avgPB = 6*(k2/p^2)*beta*con41*a*xi*I05;
coef1 = coef*(1-eta^2)^(-sym(7)/2);
C2_PB_derived = subs(n0*a0*coef*avgPB, {a0, p, beta}, {a, a, sym(1)});  % leading-e: a0=a=p, beta=1
C2_PB_code = coef1*n0*(sym(3)/4)*k2*(xi/(1-eta^2))*con41*(8 + 24*eta^2 + 3*eta^4);
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.6 Part-B assembly: C2_PartB = (3/4) k2 (xi/psi^2) con41 (8+24h^2+3h^4)', ...
    C2_PB_derived, C2_PB_code, pass_count, fail_count, failed_names);

printf('\n');
printf('Full C2 code-match witness (A.8):\n');
% Derived Part-A polynomial from the Lane-integral averages: psi^7 * <a0 + e a1>.
PartA_poly_derived = simplify((1-eta^2)^(sym(7)/2)*(avg_a0 + e*avg_a1));
% Derived Part-B bracket (verified above).
PartB_bracket_derived = (sym(3)/4)*k2*(xi/(1-eta^2))*con41*(8 + 24*eta^2 + 3*eta^4);
C2_derived = coef1*n0*(a0*PartA_poly_derived + PartB_bracket_derived);
% Literal code (drag_coefficients.h:146-149), Horner Part-B poly 8+3h^2(8+h^2):
C2_code = coef1*n0*(a0*(1 + sym(3)/2*eta^2 + e*eta*(4+eta^2)) ...
          + (sym(3)/4)*k2*xi/(1-eta^2)*con41*(8 + 3*eta^2*(8+eta^2)));
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.8 FULL C2: simplify(C2_derived - C2_code) = 0', ...
    C2_derived, C2_code, pass_count, fail_count, failed_names);

printf('\n');
printf('C1 = B*.C2 secular-decay coefficient (A.9):\n');
syms bstar real;
% A.9.1: SGP4 parameterization a(t)=a0 tempa^2, tempa=1-C1 t-...  =>  adot(0) = -2 a0 C1.
% A.9.2: the master drag-rate prefactor -(2 n0 a0^2) B* coef psi^-7 (A.1.3 at epoch)
%        factors as -2 a0 B* (coef1 n0 a0) for Part A, exact since coef1 = coef psi^-7.
prefac_master   = 2*n0*a0^2*coef*(1-eta^2)^(-sym(7)/2);   % |A.1.3 Part-A prefactor| / B*
prefac_factored = 2*a0*(coef1*n0*a0);                      % -2 a0 B* (A.5.2 group) / B*
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.9.2 secular-rate prefactor factor-out: (2 n0 a0^2) coef psi^-7 = 2 a0 (coef1 n0 a0)', ...
    prefac_master, prefac_factored, pass_count, fail_count, failed_names);
% A.9.3: <adot> = -2 a0 B* C2 and adot(0) = -2 a0 C1  =>  C1 = B* C2 (code line 152).
% With A.8 (C2_derived = C2_code), C1_derived = B* C2_derived = C1_code.
C1_derived = bstar*C2_derived;
C1_code    = bstar*C2_code;     % drag_coefficients.h:152  dc.C1 = in.bstar * dc.C2
[pass_count, fail_count, failed_names] = check_eq( ...
    'A.9.3 C1 = B*.C2: simplify(C1_derived - C1_code) = 0   [code line 152; closes Q-4]', ...
    C1_derived, C1_code, pass_count, fail_count, failed_names);

printf('\n');
printf('=== Phase 2.A/2.B summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
if fail_count > 0
    printf('\n  FAILED checks:\n');
    for k = 1:length(failed_names)
        printf('    %s\n', failed_names{k});
    end
    printf('\n  Phase 2.A C2 verification: FAIL\n');
    exit(1);
else
    printf('\n  Phase 2.A/2.B verification: PASS  (C2 Part A + Part B + full C2; C1 = B*.C2)\n');
    exit(0);
end
