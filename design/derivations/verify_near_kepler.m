% verify_near_kepler.m
%
% Standard-10 verifier (simplify(code - reference) == 0) for the SGP4 modified
% Kepler solver in src/orbit/modified_kepler.h (solve_kepler_newton /
% solve_kepler_halley), against the canonical SGP4 iteration
% (Vallado sgp4(); dnwrnr/python-sgp4 propagation.py Kepler loop). Roadmap W3.
%
% SGP4 solves the rotated Kepler equation for x = E + omega:
%   x + ayn*cos(x) - axN*sin(x) = U          (axN = e cos w, ayn = e sin w + lp)
% i.e. the residual f(x) = U - ayn*cos(x) + axN*sin(x) - x = 0.
%
% Vallado's Newton step is  x += (U - ayn*cos x + axN*sin x - x)
%                                / (1 - axN*cos x - ayn*sin x).
% The code computes f = U - ayn*cx + axN*sx - x and
% fp = -ayn*sx - axN*cx + 1 (= 1 - axN*cx - ayn*sx), then x += f/fp -- exactly
% Vallado's numerator/denominator. Note fp = -f'(x): the code's denominator is
% MINUS the true derivative, so "x += f/fp" equals the Newton step "x - f/f'".
% Halley uses fpp = ayn*cx - axN*sx = f''(x) and the cubic increment.
%
% Checks: f, fp (vs Vallado), fp == -f'(x), the Newton increment, fpp == f''(x),
% and the Halley increment == the canonical -2 f f' / (2 f'^2 - f f'').

1;
pkg load symbolic;

syms axN ayn U x

% --- CODE (modified_kepler.h) ---
f   = U - ayn*cos(x) + axN*sin(x) - x;
fp  = -ayn*sin(x) - axN*cos(x) + 1;
fpp = ayn*cos(x) - axN*sin(x);
dN_code = f / fp;
dH_code = 2*f*fp / (2*fp^2 - f*fpp);

% --- REFERENCE (canonical SGP4 / Vallado) ---
num_ref  = U - ayn*cos(x) + axN*sin(x) - x;     % Vallado numerator
den_ref  = 1 - axN*cos(x) - ayn*sin(x);         % Vallado denominator
step_ref = num_ref / den_ref;

fprime  = diff(f, x);        % true f'(x) = ayn sin x + axN cos x - 1
fdouble = diff(f, x, 2);     % true f''(x) = ayn cos x - axN sin x
dH_ref  = -2*f*fprime / (2*fprime^2 - f*fdouble);   % canonical Halley increment

names = {'Newton residual f == Vallado numerator', ...
         'Newton denominator fp == Vallado denominator', ...
         "fp == -f'(x)  (x += f/fp is a valid Newton step)", ...
         'Newton increment f/fp == Vallado step', ...
         "Halley fpp == f''(x)  (second derivative)", ...
         "Halley increment == canonical -2 f f'/(2 f'^2 - f f'')"};
D = {f - num_ref, fp - den_ref, fp + fprime, dN_code - step_ref, fpp - fdouble, dH_code - dH_ref};

pass = 0; fail = 0; tags = {'FAIL', 'PASS'};
printf('=== near-earth modified Kepler solver: modified_kepler.h vs canonical SGP4 (SR3 p.13) ===\n\n');
for k = 1:numel(names)
  d = simplify(D{k});
  tf = false;
  try; tf = (double(d) == 0); catch; tf = false; end
  printf('  [%s] simplify(code - reference) == 0   %s\n', tags{tf+1}, names{k});
  if tf; pass = pass + 1; else; fail = fail + 1; printf('        residual = %s\n', char(d)); end
end

printf('\n=== Kepler summary: %d / %d PASS ===\n', pass, pass + fail);
if fail > 0
  printf('\n  verify_near_kepler: FAIL\n'); exit(1);
else
  printf('\n  verify_near_kepler: PASS\n'); exit(0);
end
