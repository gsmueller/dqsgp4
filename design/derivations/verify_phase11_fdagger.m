% verify_phase11_fdagger.m
%
% PHASE 11 RESOLUTION: the Lane "fictitious anomaly" f_dagger IS the eccentric anomaly E.
% The Lane substitution  r - s = (a-s)(1 - eta cos f_dagger)  is EXACT with f_dagger = E
% (since r = a(1-e cos E) and eta = a e/(a-s)).  Hence the Lane density (r-s)^-4 =
% xi^4 (1-eta cos E)^-4 is EXACT -- there is NO "f_dagger substitution error".
%
% Consequence: the C2 Part-A orbit average, done EXACTLY over E (Jacobian dM=(1-e cos E)dE),
% reduces to the Phase-1 Lane integrals I^(0,m) and reproduces the SGP4 code Part-A to O(e),
% AND yields the EXACT O(e^2) AFGP4-dropped term (3/4 e^2 + 3 e^2 eta^2) of §A.7 -- now DERIVED
% symbolically (previously only numerically validated, verify_A7_truncation_floor.m).
%
% This supersedes Phase 2.A's "f_dagger ~ f" (true-anomaly) approximation: the f_dagger-corrections
% there are exactly the E-vs-f equation-of-center terms; integrating over E is exact and cleaner.

pkg load symbolic;
syms a e s E eta real;
syms ecc psi real;          % ecc = e (alias), psi = sqrt(1-eta^2)

pass=0; fail=0; names={};
function [p,fl,nm]=ck(name,lhs,rhs,p,fl,nm)
  d=simplify(lhs-rhs);
  if isequal(d,sym(0)); printf('  PASS: %s\n',name); p=p+1;
  else printf('  FAIL: %s\n        lhs-rhs=%s\n',name,char(d)); fl=fl+1; nm{end+1}=name; end
endfunction

printf('=== Phase 11: the Lane f_dagger = eccentric anomaly E (exact density) ===\n\n');

% (P11.1) f_dagger = E exactly:  r - s = (a-s)(1 - eta cos E),  r = a(1-e cos E), eta = a e/(a-s).
eta_def = a*e/(a-s);
r_E = a*(1-e*cos(E));
[pass,fail,names]=ck('P11.1 r-s = (a-s)(1-eta cos E) EXACTLY  => Lane f_dagger = E (eccentric anomaly)', ...
   (r_E - s), (a-s)*(1-eta_def*cos(E)), pass,fail,names);

% (P11.2) EXACT C2 Part-A integrand over E (kinematic factor): using
%   1+e cos f = beta^2/(1-e cos E),  1+e^2+2e cos f = beta^2(1+e cos E)/(1-e cos E),
%   dM = (1-e cos E)dE,  (r-s)^-4 = xi^4 (1-eta cos E)^-4,
%   the integrand reduces to (1-e^2 cos^2 E)^{3/2}/[(1-e cos E)^2 (1-eta cos E)^4].
%   Its kinematic part expands as 1 + 2e cosE + (3/2)e^2 cos^2 E + O(e^3).
kin = (1-ecc^2*cos(E)^2)^(sym(3)/2) / (1-ecc*cos(E))^2;
kin_series = taylor(kin, ecc, 0, 'order', 3);     % through e^2
[pass,fail,names]=ck('P11.2 kinematic factor (1-e^2cos^2E)^{3/2}/(1-e cosE)^2 = 1 + 2e cosE + (3/2)e^2 cos^2E + O(e^3)', ...
   kin_series, 1 + 2*ecc*cos(E) + (sym(3)/2)*ecc^2*cos(E)^2, pass,fail,names);

% (P11.3,4) cos^k E reductions onto Lane integrals via cos E = (1-(1-eta cos E))/eta:
%   cosE   /(1-eta cosE)^4 = (1/eta) [ (.)^-4 - (.)^-3 ]
%   cos^2E /(1-eta cosE)^4 = (1/eta^2)[ (.)^-4 - 2(.)^-3 + (.)^-2 ]
u = 1-eta*cos(E);
[pass,fail,names]=ck('P11.3 cosE/(1-eta cosE)^4 = (1/eta)[(.)^-4 - (.)^-3]', ...
   cos(E)/u^4, (1/eta)*(u^(-4) - u^(-3)), pass,fail,names);
[pass,fail,names]=ck('P11.4 cos^2E/(1-eta cosE)^4 = (1/eta^2)[(.)^-4 - 2(.)^-3 + (.)^-2]', ...
   cos(E)^2/u^4, (1/eta^2)*(u^(-4) - 2*u^(-3) + u^(-2)), pass,fail,names);

% Phase-1 Lane integrals I^(0,m)(eta) = (1/2pi) int (1-eta cos E)^-m dE  (Thm 1.3.2):
I02 = (1-eta^2)^(-sym(3)/2);
I03 = (2+eta^2)/(2*(1-eta^2)^(sym(5)/2));
I04 = (2+3*eta^2)/(2*(1-eta^2)^(sym(7)/2));

% (P11.5) CONSTANT term: <1/(1-eta cosE)^4> = I^(0,4) = (1+(3/2)eta^2) psi^-7  = code Part-A constant.
[pass,fail,names]=ck('P11.5 constant: I^(0,4) = (1 + (3/2)eta^2)(1-eta^2)^{-7/2}   [code Part-A constant]', ...
   I04, (1 + (sym(3)/2)*eta^2)*(1-eta^2)^(-sym(7)/2), pass,fail,names);

% (P11.6) e-term: 2*<cosE/(1-eta cosE)^4> = (2/eta)(I04-I03) = eta(4+eta^2) psi^-7
%         => code Part-A e-term e*eta(4+eta^2).  EXACT, over E -- no f_dagger approximation.
[pass,fail,names]=ck('P11.6 e-term: (2/eta)(I^(0,4)-I^(0,3)) = eta(4+eta^2)(1-eta^2)^{-7/2}   [code e*eta(4+eta^2)]', ...
   (2/eta)*(I04-I03), eta*(4+eta^2)*(1-eta^2)^(-sym(7)/2), pass,fail,names);

% (P11.7) e^2-term (THE §A.7 DROPPED TERM, now DERIVED):
%   (3/2)*<cos^2E/(1-eta cosE)^4> = (3/2)/eta^2 (I04 - 2 I03 + I02) = (3/4)(1+4eta^2) psi^-7,
%   i.e. the AFGP4 Part-A O(e^2) bracket term is e^2 * (3/4)(1+4eta^2) = (3/4)e^2 + 3 e^2 eta^2.
e2coef = (sym(3)/2)/eta^2 * (I04 - 2*I03 + I02);
[pass,fail,names]=ck('P11.7a e^2-coef: (3/2)/eta^2 (I04-2I03+I02) = (3/4)(1+4eta^2)(1-eta^2)^{-7/2}', ...
   e2coef, (sym(3)/4)*(1+4*eta^2)*(1-eta^2)^(-sym(7)/2), pass,fail,names);
% strip psi^-7 (factored into coef1) -> the dropped BRACKET term:
[pass,fail,names]=ck('P11.7b §A.7 dropped term DERIVED: (3/4)(1+4eta^2) = 3/4 + 3 eta^2  => (3/4 e^2 + 3 e^2 eta^2)', ...
   (sym(3)/4)*(1+4*eta^2), sym(3)/4 + 3*eta^2, pass,fail,names);

printf('\n=== Phase 11 summary: %d / %d checks PASS ===\n', pass, pass+fail);
printf('  RESOLVED: f_dagger = E (eccentric anomaly); the Lane density (r-s)^-4 = xi^4(1-eta cosE)^-4 is\n');
printf('  EXACT (P11.1).  The C2 Part-A orbit average over E reduces to the Phase-1 Lane integrals and\n');
printf('  gives the code Part-A constant 1+(3/2)eta^2 (P11.5) and e-term e*eta(4+eta^2) (P11.6) EXACTLY,\n');
printf('  with the O(e^2) AFGP4-dropped term (3/4 e^2 + 3 e^2 eta^2) now DERIVED symbolically (P11.7) --\n');
printf('  upgrading §A.7 from numerically-validated to DERIVED.  There is NO f_dagger density error; the\n');
printf('  only C2 Part-A approximation is the SGP4 O(e^2) bracket TRUNCATION (a clean, exact accounting).\n');
if fail>0
  printf('\n  FAILED: '); for k=1:numel(names); printf('%s; ',names{k}); end; printf('\n');
  printf('  Phase 11: FAIL\n'); exit(1);
else
  printf('\n  Phase 11: PASS\n'); exit(0);
end
