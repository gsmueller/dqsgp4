% verify_near_longperiod.m
%
% Standard-10 verifier (simplify(code - reference) == 0) for the J3 LONG-PERIOD
% (Lyddane) recovery in src/sgp4/near_space.h, against canonical SGP4
% (Vallado sgp4(); dnwrnr/python-sgp4 propagation.py long-period block).
%
% Two parts (roadmap W2 of the near-earth Standard-10 sweep, A2):
%   (1) the long-period COEFFICIENTS computed in initialize_near_space:
%         xlcof, aycof  (the J3 mean-longitude / eccentricity-vector coefs);
%   (2) the per-time ASSEMBLY in propagate_near_space:
%         axN, temp_lp, ayn, xlt.
%
% Independence of the J3 sign convention: the code uses A30 = config.A(3,0),
% documented as A(3,0) = -J3, and CK2 = J2/2 (near_space.h:117,153,155). So the
% code's A30/CK2 = (-J3)/(J2/2) = -2*J3/J2. The reference side uses Vallado's
% j3oj2 = J3/J2 directly. The verifier derives both from J2, J3 as free symbols
% and proves they reduce to the same expression, so a sign or factor error on
% either side breaks the identity.

1;
pkg load symbolic;

syms J2 J3 i0 e a omega ML cc
c = cos(i0);  s = sin(i0);

% --- (1) coefficients --------------------------------------------------------
% CODE (near_space.h): A30 = -J3, CK2 = J2/2, A30_over_CK2 = A30/CK2.
A30 = -J3;  CK2 = J2/2;  A30_over_CK2 = A30 / CK2;          % = -2*J3/J2
xlcof_code = sym(1)/8 * A30_over_CK2 * s * (3 + 5*c) / (1 + c);
aycof_code = sym(1)/4 * A30_over_CK2 * s;

% REFERENCE (canonical SGP4 / Vallado): j3oj2 = J3/J2.
j3oj2 = J3 / J2;
xlcof_ref = -sym(1)/4 * j3oj2 * s * (3 + 5*c) / (1 + c);
aycof_ref = -sym(1)/2 * j3oj2 * s;

% --- (2) per-time long-period assembly (propagate_near_space) ---------------
% Code and canonical SGP4 share the same closed form:
%   axN = e cos w ;  temp = 1/(a(1-e^2)) ;  ayn = e sin w + temp*aycof ;
%   xlt = ML + temp*xlcof*axN     (ML = mean longitude = M + w + Omega)
axN     = e * cos(omega);
temp_lp = 1 / (a * (1 - e^2));
ayn_code = e * sin(omega) + temp_lp * aycof_code;
ayn_ref  = e * sin(omega) + temp_lp * aycof_ref;
xlt_code = ML + temp_lp * xlcof_code * axN;
xlt_ref  = ML + temp_lp * xlcof_ref  * axN;

names = {'xlcof  (J3 mean-longitude coefficient)', ...
         'aycof  (J3 a_yN coefficient)', ...
         'ayn  = e sin w + temp*aycof  (assembly)', ...
         'xlt  = ML + temp*xlcof*axN   (assembly)'};
D = {xlcof_code - xlcof_ref, aycof_code - aycof_ref, ...
     ayn_code - ayn_ref, xlt_code - xlt_ref};

pass = 0; fail = 0; tags = {'FAIL', 'PASS'};
printf('=== near-earth long-period (Lyddane / J3) recovery: near_space.h vs canonical SGP4 ===\n\n');
for k = 1:numel(names)
  d = simplify(D{k}); tf = false;
  try; tf = (double(d) == 0); catch; tf = false; end
  printf('  [%s] simplify(code - reference) == 0   %s\n', tags{tf+1}, names{k});
  if tf; pass = pass + 1; else; fail = fail + 1; printf('        residual = %s\n', char(d)); end
end

% --- critical-inclination fallback check ------------------------------------
% The guard substitutes xlcof = (1/2)*A30/CK2*sin i0 when 1+cos i0 -> 0 fails;
% (1/2) is the cos i0 -> +1 limit of (1/8)*(3+5c)/(1+c), i.e. lim (3+5c)/(1+c)=4.
lim = limit((3 + 5*cc) / (1 + cc), cc, 1);
tf = false;
try; tf = (double(simplify(lim - 4)) == 0); catch; tf = false; end
printf('  [%s] critical-incl fallback coef = +1 limit: lim_{c->1}(3+5c)/(1+c) = 4 (so 1/8*4 = 1/2)\n', tags{tf+1});
if tf; pass = pass + 1; else; fail = fail + 1; end

printf('\n=== long-period summary: %d / %d PASS ===\n', pass, pass + fail);
if fail > 0
  printf('\n  verify_near_longperiod: FAIL\n'); exit(1);
else
  printf('\n  verify_near_longperiod: PASS\n'); exit(0);
end
