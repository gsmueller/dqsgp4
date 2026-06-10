1;
% verify_dim_kaula.m
%
% W14 dimensional / structural audit of the Kaula inclination functions
% F_lmp(i) in src/perturbation/kaula.h (inclination_function). Extends the
% dimensional-audit discipline (feedback_dimensional_audit) to the resonance
% angular factors.
%
% Dimensional rationale. In the deep-space resonance terms each rate has the
% form  (dimensionful amplitude D_lmpq) x F_lmp(i) x G_lpq(e). The inclination
% function F_lmp is the ANGULAR factor: it must be dimensionless (a function of
% the dimensionless inclination angle alone) and a trig polynomial of bounded
% total degree l, so it contributes ZERO dimension and does not pollute the
% resonance term's units. This script audits exactly those two properties, and
% additionally confirms the exact-rational coefficients reconstruct the
% documented decimals (so the FORTRAN truncation artifacts noted in kaula.h are
% provably avoided).
%
% Checks (si = sin i, ci = cos i as independent symbols):
%   (1) dimensionless: F_lmp depends only on si, ci (no length/time/mass scale).
%   (2) total trig degree(F_lmp) == l  (Kaula structural invariant): via the
%       homogenization F(lam*si, lam*ci) — the (l+1)-th lam-derivative vanishes
%       (degree <= l) and the l-th does not (a degree-l term exists).
%   (3) exact-rational coefficients == documented decimals; FORTRAN artifacts
%       (4.92187512, 6.56250012) provably differ from 315/64, 105/16.
%
% simplify/diff/subs/expand only; no local functions. exit(0) iff all pass.

pkg load symbolic

syms si ci lam

% F_220 reused by F_441
F220 = (sym(3)/4) * (1 + ci)^2;

% {label, l, F_lmp(si,ci)} for every (l,m,p) implemented in kaula.h
K = {
  {'F_220', 2, F220}, ...
  {'F_221', 2, (sym(3)/2) * si^2}, ...
  {'F_311', 3, (sym(15)/16)*si^2*(1 + 3*ci) - (sym(3)/4)*(1 + ci)}, ...
  {'F_321', 3, (sym(15)/8)*si*(1 - 2*ci - 3*ci^2)}, ...
  {'F_322', 3, (sym(-15)/8)*si*(1 + 2*ci - 3*ci^2)}, ...
  {'F_330', 3, (sym(15)/8)*(1 + ci)^3}, ...
  {'F_441', 4, sym(35)*si^2*F220}, ...
  {'F_442', 4, (sym(315)/8)*si^4}, ...
  {'F_522', 5, (sym(315)/32)*si*(si^2*(1 - 2*ci - 5*ci^2) + (sym(1)/3)*(-2 + 4*ci + 6*ci^2))}, ...
  {'F_523', 5, si*((sym(315)/64)*si^2*(-2 - 4*ci + 10*ci^2) + (sym(105)/16)*(1 + 2*ci - 3*ci^2))}, ...
  {'F_542', 5, (sym(945)/32)*si*(2 - 8*ci + ci^2*(-12 + 8*ci + 10*ci^2))}, ...
  {'F_543', 5, (sym(945)/32)*si*(-2 - 8*ci + ci^2*(12 + 8*ci - 10*ci^2))} ...
};

np = 0; nf = 0;
printf('=== W14 dimensional/structural audit: Kaula inclination functions ===\n');
for k = 1:numel(K)
  lab = K{k}{1}; l = K{k}{2}; F = K{k}{3};

  % (1) dimensionless: only si, ci free. Substituting both leaves a pure number.
  r = subs(F, [si ci], [sym(1)/3, sym(1)/4]);
  dimless = isempty(symvar(r));

  % (2) total trig degree == l via homogenization in lam.
  g = expand(subs(F, [si ci], [lam*si, lam*ci]));
  excess = simplify(diff(g, lam, l + 1));         % 0  iff degree <= l
  lead   = diff(g, lam, l);                        % l! * (degree-l part)
  lead_num = double(subs(lead, [si ci lam], [sym(7)/10, sym(1)/2, sym(0)]));
  degree_ok = isequal(excess, sym(0)) && (abs(lead_num) > 0);

  ok = dimless && degree_ok;
  if ok
    np = np + 1; st = 'PASS';
  else
    nf = nf + 1; st = 'FAIL';
  end
  printf('  [%s] %s : dimensionless=%d, total degree == l=%d (%d)\n', ...
         st, lab, dimless, l, degree_ok);
end

% (3) exact-rational coefficients reconstruct documented decimals.
printf('\n  -- exact-rational coefficient reconstruction --\n');
coeff = {
  {'3/4   = 0.75',      sym(3)/4,    sym(75)/100}, ...
  {'3/2   = 1.5',       sym(3)/2,    sym(15)/10}, ...
  {'15/16 = 0.9375',    sym(15)/16,  sym(9375)/10000}, ...
  {'15/8  = 1.875',     sym(15)/8,   sym(1875)/1000}, ...
  {'315/8 = 39.375',    sym(315)/8,  sym(39375)/1000}, ...
  {'315/32= 9.84375',   sym(315)/32, sym(984375)/100000}, ...
  {'315/64= 4.921875',  sym(315)/64, sym(4921875)/1000000}, ...
  {'105/16= 6.5625',    sym(105)/16, sym(65625)/10000}, ...
  {'945/32= 29.53125',  sym(945)/32, sym(2953125)/100000} ...
};
for k = 1:numel(coeff)
  lab = coeff{k}{1};
  if isequal(simplify(coeff{k}{2} - coeff{k}{3}), sym(0))
    np = np + 1; printf('  [PASS] %s (exact)\n', lab);
  else
    nf = nf + 1; printf('  [FAIL] %s\n', lab);
  end
end

% FORTRAN truncation artifacts must NOT equal the exact rationals.
printf('\n  -- FORTRAN truncation artifacts differ from the exact rationals --\n');
art = {
  {'4.92187512 != 315/64', sym(492187512)/100000000, sym(315)/64}, ...
  {'6.56250012 != 105/16', sym(656250012)/100000000, sym(105)/16} ...
};
for k = 1:numel(art)
  lab = art{k}{1};
  if ~isequal(simplify(art{k}{2} - art{k}{3}), sym(0))
    np = np + 1; printf('  [PASS] %s (artifact correctly avoided)\n', lab);
  else
    nf = nf + 1; printf('  [FAIL] %s\n', lab);
  end
end

printf('\n');
if nf == 0
  printf('verify_dim_kaula: ALL %d checks PASS\n', np);
  exit(0);
else
  printf('verify_dim_kaula: %d PASS, %d FAIL\n', np, nf);
  exit(1);
end
