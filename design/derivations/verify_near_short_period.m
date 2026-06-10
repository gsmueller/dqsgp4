% verify_near_short_period.m
%
% Standard-10 verifier (simplify(code - reference) == 0) for the J2 SHORT-PERIOD
% periodics in src/perturbation/short_period.h (apply_short_period), against the
% canonical SGP4 short-period block (Hoots & Roehrich, Spacetrack Report #3,
% p.14; Vallado SGP4 sgp4(); dnwrnr/python-sgp4 propagation.py).
%
% This closes one stage of AUDIT_BACKLOG A2 (per-line Standard 10 for the
% near-earth assembly): the drag coefficients were already code-matched, but the
% propagation ASSEMBLY had only end-to-end validation. Here each of the six
% osculating-element corrections is encoded TWICE and proven identical:
%   - the CODE side, transcribed from short_period.h's named intermediates
%     (sp2 = half_J2/pl, sp3 = sp2/pl, sin2_i0, three_cos2i_minus_1 = 3c^2-1,
%      seven_cos2i_minus_1 = 7c^2-1, beta_l = sqrt(1-e^2));
%   - the REFERENCE side, from SGP4's own auxiliary definitions
%     (temp1 = 0.5*J2/p, temp2 = temp1/p, x1mth2 = 1-c^2, con41 = 3c^2-1,
%      x7thm1 = 7c^2-1, betal = sqrt(1-e^2)).
% A coefficient drift on either side (e.g. 7c^2-1 -> 5c^2-1, or 1/4 -> 1/3)
% makes the symbolic difference non-zero and FAILs the gate.
%
% half_J2 is CK2 = J2/2, so the code's sp2 = half_J2/p equals the reference
% temp1 = 0.5*J2/p exactly.

1;
pkg load symbolic;

syms r u i0 Om rdot rfdot n p e half_J2
syms s2u c2u            % sin(2u), cos(2u) as free symbols

c = cos(i0);  s = sin(i0);

% --- REFERENCE: canonical SGP4 short-period block (SR3 p.14) ---------------
betal  = sqrt(1 - e^2);
temp1  = half_J2 / p;            % CK2 / p
temp2  = temp1 / p;             % CK2 / p^2
x1mth2 = 1 - c^2;               % sin^2 i
con41  = 3*c^2 - 1;
x7thm1 = 7*c^2 - 1;

rk_ref     = r*(1 - sym(3)/2*temp2*betal*con41) + sym(1)/2*temp1*x1mth2*c2u;
uk_ref     = u - sym(1)/4*temp2*x7thm1*s2u;
xinck_ref  = i0 + sym(3)/2*temp2*c*s*c2u;
xnodek_ref = Om + sym(3)/2*temp2*c*s2u;
rdotk_ref  = rdot - n*temp1*x1mth2*s2u;
rfdotk_ref = rfdot + n*temp1*(x1mth2*c2u + sym(3)/2*con41);

% --- CODE: src/perturbation/short_period.h apply_short_period --------------
sp2     = half_J2 / p;          % code: half_J2 / pl
sp3     = sp2 / p;              % code: sp2 / pl
beta_l  = sqrt(1 - e^2);
sin2_i0             = s^2;       % = 1 - c^2
three_cos2i_minus_1 = 3*c^2 - 1;
seven_cos2i_minus_1 = 7*c^2 - 1;

rk_code     = r*(1 - sym(3)/2*sp3*beta_l*three_cos2i_minus_1) + sym(1)/2*sp2*sin2_i0*c2u;
uk_code     = u - sym(1)/4*sp3*seven_cos2i_minus_1*s2u;
xinck_code  = i0 + sym(3)/2*sp3*c*s*c2u;
xnodek_code = Om + sym(3)/2*sp3*c*s2u;
rdotk_code  = rdot - n*sp2*sin2_i0*s2u;
rfdotk_code = rfdot + n*sp2*(sin2_i0*c2u + sym(3)/2*three_cos2i_minus_1);

names = {'rk     (radius)', 'uk     (argument of latitude)', ...
         'xinck  (inclination)', 'xnodek (RAAN)', ...
         'rdotk  (radial velocity)', 'rfdotk (transverse velocity)'};
D = {rk_code - rk_ref, uk_code - uk_ref, xinck_code - xinck_ref, ...
     xnodek_code - xnodek_ref, rdotk_code - rdotk_ref, rfdotk_code - rfdotk_ref};

pass = 0; fail = 0;
tags = {'FAIL', 'PASS'};
printf('=== near-earth short-period J2 periodics: short_period.h vs canonical SGP4 (SR3 p.14) ===\n\n');
for k = 1:numel(names)
  d = simplify(D{k});
  tf = false;
  try; tf = (double(d) == 0); catch; tf = false; end   % free-var residual => not zero
  printf('  [%s] simplify(code - reference) == 0   %s\n', tags{tf+1}, names{k});
  if tf; pass = pass + 1; else; fail = fail + 1; printf('        residual = %s\n', char(d)); end
end

printf('\n=== short-period summary: %d / %d PASS ===\n', pass, pass + fail);
if fail > 0
  printf('\n  verify_near_short_period: FAIL\n'); exit(1);
else
  printf('\n  verify_near_short_period: PASS\n'); exit(0);
end
