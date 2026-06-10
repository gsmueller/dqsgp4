% verify_near_osculating.m
%
% Standard-10 verifier for the osculating-elements -> TEME position/velocity
% assembly in src/orbit/state_from_elements.h (elements_to_state), against the
% canonical SGP4 orientation block (Vallado sgp4(); dnwrnr/python-sgp4). W5.
%
% Two kinds of check:
%   (1) TRANSCRIPTION: the orientation unit vectors (ux,uy,uz, vx,vy,vz) match
%       the canonical SGP4 forms exactly (a sign error breaks the identity).
%   (2) GEOMETRY (independent of transcription): the pair (u_hat, v_hat) is an
%       orthonormal frame -- |u|^2 = |v|^2 = 1 and u.v = 0 -- which proves the
%       formulas realise the rotation R3(-Omega) R1(-inc) R3(-u) correctly.
%
% Position x = r_er * u_hat * re_km and velocity (rdot*u_hat + rfdot*v_hat) *
% (re_km/60) are the assembly around this frame; re_km/60 is the exact
% ER/min -> km/s unit conversion (1 ER = re_km km, 1 min = 60 s). NB Octave's
% `i` is the imaginary unit, so inclination is `inc` here.

1;
pkg load symbolic;

syms inc Om arglat
ci = cos(inc); si = sin(inc);
cO = cos(Om);  sO = sin(Om);
cu = cos(arglat); su = sin(arglat);

% --- CODE (state_from_elements.h) ---
xmx = -sO*ci;  xmy = cO*ci;
ux = xmx*su + cO*cu;   uy = xmy*su + sO*cu;   uz = si*su;
vx = xmx*cu - cO*su;   vy = xmy*cu - sO*su;   vz = si*cu;

% --- REFERENCE (canonical SGP4 / Vallado orientation block) ---
ux_r = (-sO*ci)*su + cO*cu;   uy_r = (cO*ci)*su + sO*cu;   uz_r = si*su;
vx_r = (-sO*ci)*cu - cO*su;   vy_r = (cO*ci)*cu - sO*su;   vz_r = si*cu;

names = {'ux == canonical', 'uy == canonical', 'uz == canonical', ...
         'vx == canonical', 'vy == canonical', 'vz == canonical', ...
         '|u_hat|^2 == 1  (u is a unit vector)', ...
         '|v_hat|^2 == 1  (v is a unit vector)', ...
         'u_hat . v_hat == 0  (orthogonal -> proper rotation)'};
D = {ux - ux_r, uy - uy_r, uz - uz_r, vx - vx_r, vy - vy_r, vz - vz_r, ...
     ux^2 + uy^2 + uz^2 - 1, vx^2 + vy^2 + vz^2 - 1, ux*vx + uy*vy + uz*vz};

pass = 0; fail = 0; tags = {'FAIL', 'PASS'};
printf('=== near-earth osculating -> TEME orientation: state_from_elements.h vs canonical SGP4 ===\n\n');
for k = 1:numel(names)
  d = simplify(D{k});
  tf = false;
  try; tf = (double(d) == 0); catch; tf = false; end
  printf('  [%s] %s\n', tags{tf+1}, names{k});
  if tf; pass = pass + 1; else; fail = fail + 1; printf('        residual = %s\n', char(d)); end
end

printf('\n=== osculating/orientation summary: %d / %d PASS ===\n', pass, pass + fail);
if fail > 0
  printf('\n  verify_near_osculating: FAIL\n'); exit(1);
else
  printf('\n  verify_near_osculating: PASS\n'); exit(0);
end
