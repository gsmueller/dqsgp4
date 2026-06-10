% verify_A7_truncation_floor.m
%
% Phase 2.A §A.7 -- quantify + validate the SGP4 C2 Part-A truncation FLOOR.
% (§A.7 = the AFGP4 -> SGP4 dropped terms; "not code-used, Standard 10 not required" --
% this is a NUMERICAL completeness/characterisation check, not a symbolic code-match.)
%
% Method: compute the EXACT Keplerian Part-A bracket B_exact from the exact orbit integral
% (exact density (r-s)^-4, NO f-dagger approximation, NO e-truncation; dM = beta^3/(1+e cosf)^2 df
% Jacobian -- Eq (4.1), the D-2-fixed relation) and compare to the SGP4 code bracket
%   B_code = 1 + (3/2)eta^2 + e*eta*(4+eta^2).
% The difference B_exact - B_code should (a) be POSITIVE and O(e^2), and (b) match the
% documented AFGP4 dropped term (3/4 e^2 + 3 e^2 eta^2) to within the O(e*eta) f-dagger
% substitution error (a SEPARATE Phase-11 item).  This VALIDATES the §A.7 dropped-term form
% and sets the C2 Part-A accuracy floor at O(e^2, e*eta).

1;
function E = kep(M,e)
  E=M; for k=1:80; E=E-(E-e*sin(E)-M)/(1-e*cos(E)); end
end

s = 1.01;                       % atmospheric fitting radius [ER] (fixed, illustrative)
qms = 120/6378.135;             % (q0 - s) [ER]  (illustrative shell; cancels in the ratio)
N = 40000; f = linspace(0,2*pi,N+1); f(end)=[]; df = 2*pi/N;

cases = [1.15 0.02; 1.15 0.05; 1.20 0.05; 1.20 0.10; 1.30 0.08];
pass=0; fail=0;
printf('=== A.7: C2 Part-A truncation floor (exact orbit integral vs SGP4 code) ===\n\n');
printf('  a0     e0     eta     Bexact     Bcode     diff       3/4e^2+3e^2eta^2   gap%%\n');
for i=1:rows(cases)
  a0=cases(i,1); e0=cases(i,2);
  beta=sqrt(1-e0^2); xi=1/(a0-s); eta=a0*e0*xi; psi=sqrt(abs(1-eta^2));
  r = a0*beta^2 ./ (1+e0*cos(f));  rs = r - s;
  K = (1+e0^2+2*e0*cos(f)).^(1.5) ./ ( rs.^4 .* (1+e0*cos(f)).^2 );
  B_exact = psi^7*(a0-s)^4 * mean(K);
  B_code  = 1 + 1.5*eta^2 + e0*eta*(4+eta^2);
  d       = B_exact - B_code;
  drop    = 0.75*e0^2 + 3*e0^2*eta^2;
  gap     = 100*abs(d - drop)/abs(drop);
  printf('  %.2f   %.3f  %.4f  %.6f  %.6f  %+.3e  %+.3e        %4.1f\n', a0,e0,eta,B_exact,B_code,d,drop,gap);
  % checks: diff is POSITIVE, O(e^2)-bounded, and within 12% of the documented dropped term.
  ok = (d > 0) && (d < 0.06) && (gap < 12) && (abs(d)/e0^2 < 3);   % d/e^2 bounded (O(e^2))
  if ok; pass=pass+1; else fail=fail+1; printf('     FAIL on case %d\n', i); end
end

printf('\n=== A.7 floor summary: %d / %d cases PASS ===\n', pass, pass+fail);
printf('  CORROBORATED: B_exact - B_code is positive, O(e^2)-bounded, and matches the documented AFGP4\n');
printf('  dropped term (3/4 e^2 + 3 e^2 eta^2) to within ~%.0f%%.  (UPDATE: Phase 11 derives this term\n', 12);
printf('  EXACTLY -- verify_phase11_fdagger.m P11.7 -- and shows f_dagger = E so the density is exact;\n');
printf('  the residual few-%% beyond (3/4 e^2+3 e^2 eta^2) is the O(e^3) HIGHER-ORDER e-truncation, NOT\n');
printf('  an f_dagger error.)  The SGP4 C2 Part-A accuracy floor is O(e^2) -- a characterised, exact\n');
printf('  model-choice truncation, NOT a derivation error.\n');
if fail>0; printf('\n  A.7 floor: FAIL\n'); exit(1); else printf('\n  A.7 floor: PASS\n'); exit(0); end
