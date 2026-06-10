% verify_near_secular.m
%
% Standard-10 verifier (simplify(code - reference) == 0) for the secular-advance
% + drag assembly in src/orbit/secular_update.h (secular_advance), against the
% canonical SGP4 secular block (Vallado sgp4() lines ~1620-1660; dnwrnr/
% python-sgp4 propagation.py). Roadmap W4.
%
% The non-simple (full-drag) branch is verified; the simple-model branch is the
% same expressions with D2=D3=D4=C5=t3cof=t4cof=t5cof=omgcof=xmcof=0, so it is a
% strict subset. Both code and reference are transcribed from their own variable
% names so a coefficient/structure error on either side breaks the identity.

1;
pkg load symbolic;

syms t C1 C4 C5 D2 D3 D4 t2cof t3cof t4cof t5cof
syms omgcof xmcof delmo sinmo eta bstar nodecf
syms M0 w0 O0 Mdot wdot Odot a0 e0 ke

t2 = t*t;  t3 = t2*t;  t4 = t3*t;
Mdf    = M0 + Mdot*t;
argpdf = w0 + wdot*t;
nodedf = O0 + Odot*t;

% --- REFERENCE (canonical SGP4 / Vallado, non-simple) -----------------------
delomg_r = omgcof*t;
delm_r   = xmcof*((1 + eta*cos(Mdf))^3 - delmo);
mp_r     = Mdf + delomg_r + delm_r;
argpm_r  = argpdf - (delomg_r + delm_r);
nodem_r  = nodedf + nodecf*t2;
tempa_r  = 1 - C1*t - D2*t2 - D3*t3 - D4*t4;
tempe_r  = bstar*C4*t + bstar*C5*(sin(mp_r) - sinmo);
templ_r  = t2cof*t2 + t3cof*t3 + t4*(t4cof + t*t5cof);
am_r     = a0*tempa_r^2;
em_r     = e0 - tempe_r;
nm_r     = ke/(am_r*sqrt(am_r));

% --- CODE (secular_update.h, non-simple branch) -----------------------------
M_secular     = M0 + Mdot*t;
omega_secular = w0 + wdot*t;
Omega_secular = O0 + Odot*t;
delomg_c = omgcof*t;
delm_c   = xmcof*((1 + eta*cos(M_secular)) * (1 + eta*cos(M_secular)) ...
                  * (1 + eta*cos(M_secular)) - delmo);    % code uses (1+e)^3 expanded
temp_drag = delomg_c + delm_c;
M_c     = M_secular + temp_drag;            % state.M  (before n0*templ)
omega_c = omega_secular - temp_drag;        % state.omega
Omega_c = Omega_secular + nodecf*t2;        % state.Omega
a_drag_factor = 1 - C1*t - D2*t2 - D3*t3 - D4*t4;
e_drag = bstar*C4*t + bstar*C5*(sin(M_c) - sinmo);
L_drag = t2cof*t2 + t3cof*t3 + t4*(t4cof + t*t5cof);
a_c = a0*a_drag_factor*a_drag_factor;       % state.a
e_c = e0 - e_drag;                          % state.e  (before the 1e-6 floor)
n_c = ke/(a_c*sqrt(a_c));                   % state.n

names = {'tempa  (a-drag factor)', 'tempe  (e-drag)', 'templ  (L-drag)', ...
         'delomg', 'delm   (eta cubic)', 'mp     = Mdf + delomg + delm', ...
         'argpm  = argpdf - (delomg+delm)', 'nodem  = nodedf + nodecf*t^2', ...
         'am     = a0*tempa^2', 'em     = e0 - tempe', 'nm     = ke/(a*sqrt(a))'};
D = {a_drag_factor - tempa_r, e_drag - tempe_r, L_drag - templ_r, ...
     delomg_c - delomg_r, delm_c - delm_r, M_c - mp_r, omega_c - argpm_r, ...
     Omega_c - nodem_r, a_c - am_r, e_c - em_r, n_c - nm_r};

pass = 0; fail = 0; tags = {'FAIL', 'PASS'};
printf('=== near-earth secular advance + drag assembly: secular_update.h vs canonical SGP4 ===\n\n');
for k = 1:numel(names)
  d = simplify(D{k});
  tf = false;
  try; tf = (double(d) == 0); catch; tf = false; end
  printf('  [%s] simplify(code - reference) == 0   %s\n', tags{tf+1}, names{k});
  if tf; pass = pass + 1; else; fail = fail + 1; printf('        residual = %s\n', char(d)); end
end

printf('\n=== secular summary: %d / %d PASS ===\n', pass, pass + fail);
if fail > 0
  printf('\n  verify_near_secular: FAIL\n'); exit(1);
else
  printf('\n  verify_near_secular: PASS\n'); exit(0);
end
