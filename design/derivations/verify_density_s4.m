% verify_density_s4.m
%
% Guards the low-perigee atmospheric-parameter (s4 / qoms4) fix in
% src/atmosphere/density_model.h (2026-06-03).  For 98 <= perigee < 156 km the
% adjusted s had a spurious +1 (a full Earth radius, ~6378 km), catastrophically
% corrupting the drag coefficients (sat 28350, perigee 135.75 km, failed by 1.1e4 km).
%
% Reference SGP4 (Hoots-Roehrich / Vallado / Rhodes propagation.py:1429-1437):
%   if perige < 156:
%       sfour = perige - 78        (km);  if perige < 98: sfour = 20
%       qzms24 = ((120 - sfour)/re)^4
%       sfour  = sfour/re + 1      (ER)
% with perige = (a0*(1-e0) - 1)*re  (perigee ALTITUDE, km).
%
% The (fixed) code computes, in ER:  s_default = 1 + 78/re;
%   98<=perige<156:  s_star = a0*(1-e0) - s_default;       s = 1 + s_star
%                    new_qoms = (qoms4_def)^(1/4) + s_default - s;  qoms4 = new_qoms^4
%   perige<98:       s_star = 20/re;                       s = 1 + s_star
%                    new_qoms = (qoms4_def)^(1/4) + s_default - s;  qoms4 = new_qoms^4
% This file checks code == reference across perigees (and that the OLD +1 form did NOT).

1;
re = 6378.135;
s_def = 1 + 78/re;
qoms4_def = ((120-78)/re)^4;          % default (q0-s)^4
root_def  = (qoms4_def)^(1/4);        % = 42/re

function [s,q] = ref(perige, re)
  sfour = perige - 78; if perige < 98; sfour = 20; end
  q = ((120 - sfour)/re)^4;
  s = sfour/re + 1;
end
function [s,q] = code_fixed(perige, a1me, re, s_def, root_def)
  if perige < 156 && perige >= 98
    s_star = a1me - s_def;                       % FIX: no +1
  else  % perige < 98
    s_star = 20/re;
  end
  s = 1 + s_star;
  new_q = root_def + s_def - s;
  q = new_q^4;
end
function [s,q] = code_buggy(perige, a1me, re, s_def, root_def)
  if perige < 156 && perige >= 98
    s_star = a1me - s_def + 1;                   % OLD spurious +1
  else
    s_star = 20/re;
  end
  s = 1 + s_star;
  new_q = root_def + s_def - (1 + s_star)*0 - s; %#ok  (s already = 1+s_star)
  new_q = root_def + s_def - s;
  q = new_q^4;
end

pass=0; fail=0;
printf('=== density s4 / qoms4 low-perigee fix vs reference SGP4 ===\n\n');
printf('  perige  s_ref      s_codeFIXED  ds        q_ref/q_codeFIXED-1   s_codeBUGGY\n');
for perige = [135.75 100.0 120.0 150.0 90.0 60.0]
  a1me = 1 + perige/re;                          % a0(1-e0) [ER]
  [sr,qr]   = ref(perige, re);
  [sf,qf]   = code_fixed(perige, a1me, re, s_def, root_def);
  [sb,~]    = code_buggy(perige, a1me, re, s_def, root_def);
  ds = abs(sf - sr); dq = abs(qf/qr - 1);
  ok = (ds < 1e-12) && (dq < 1e-12);
  printf('  %6.2f  %.7f  %.7f  %.1e   %.1e            %.6f\n', perige, sr, sf, ds, dq, sb);
  if ok; pass=pass+1; else fail=fail+1; end
  % also assert the OLD buggy form was wrong by ~1 ER for the 98-156 branch:
  if perige>=98 && perige<156
    if abs(sb - (sr+1)) < 1e-9; pass=pass+1; else fail=fail+1; printf('    (buggy-form check FAILED at %.2f)\n',perige); end
  end
end

printf('\n=== density s4 summary: %d / %d checks PASS ===\n', pass, pass+fail);
printf('  FIXED code matches reference SGP4 (s and qoms4) for both low-perigee branches; the OLD\n');
printf('  98<=perige<156 form was high by ~1 Earth radius in s.  (Regression: sat 28350 13/13.)\n');
if fail>0; printf('\n  density s4: FAIL\n'); exit(1); else printf('\n  density s4: PASS\n'); exit(0); end
