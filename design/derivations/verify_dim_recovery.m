1;
% verify_dim_recovery.m
%
% W15 dimensional audit of the mean-element recovery boxed formulas in
% src/orbit/element_recovery.h (recover_mean_elements). Extends the
% dimensional-audit discipline (feedback_dimensional_audit) to the
% Kozai->Brouwer element recovery.
%
% Unit system. SGP4 normalizes length to the Earth radius aE (so a, a1, a0 and
% the J2 corrections del1, del0 are dimensionless ratios) and keeps angular
% rates in rad/min. So the only running dimension in the recovery is TIME: a
% rate is Ts^-1, a period is Ts^+1, a length ratio / correction is Ts^0. The
% perigee output additionally carries a km length (Lkm) through re_km.
%
% Audited claim (per boxed formula):
%   a1, del1, a0_corrected, del0, a0_recovered, tail_bound : Ts^0  (dimensionless)
%   n0 (recovered mean motion)                              : Ts^-1 (a rate)
%   period_min                                              : Ts^+1 (a time)
%   perigee_km                                              : Ts^0 * Lkm^1
%   accuracy adds: a0 += tail (Ts^0); n0 += (3/2)(n0/a0) tail (Ts^-1)
% and the physical Kepler consistency  a1 * aE = (GM / n^2)^(1/3), confirming
% a1 is genuinely the (normalized) semi-major axis with k_e = sqrt(GM/aE^3).
%
% Method. Ts enters only through n = n0s/Ts and ke = ke0/Ts, which cancel or
% factor out, so every recovered quantity is a clean MONOMIAL in Ts (and at
% most Lkm^1). The Ts- and Lkm-exponents are therefore read off exactly from
% two numeric evaluations: p = log2(q(v=2)/q(v=1)). A non-integer result would
% flag a non-monomial (mixed-dimension) formula and fails. (Two-point numeric
% avoids symbolic simplify of the deeply nested fractional-power expressions.)
% exit(0) iff every formula audits clean.

pkg load symbolic

syms Ts Lkm ke0 n0s k2s e cosi positive

n  = n0s / Ts;          % rate, Ts^-1
ke = ke0 / Ts;          % rate, Ts^-1

beta02 = 1 - e^2;
beta0  = sqrt(beta02);
x3thm1 = 3*cosi^2 - 1;

a1   = (ke / n)^(sym(2)/3);                                   % Keplerian a (aE units)
del1 = (sym(3)/2) * k2s * x3thm1 / (a1^2 * beta0 * beta02);   % first J2 correction
a0c  = a1 * (1 - del1/3 - del1^2 - (sym(134)/81)*del1^3);     % corrected a (step 3)
del0 = (sym(3)/2) * k2s * x3thm1 / (a0c^2 * beta0 * beta02);  % refined correction
n0   = n / (1 + del0);                                        % recovered mean motion
a0r  = (ke / n0)^(sym(2)/3);                                  % recovered a (step 5)
tail = del1^4 * a1;                                           % O(del1^4) tail bound
a0_acc = tail;                                               % -> result.a0.accuracy
n0_acc = (sym(3)/2) * (n0 / a0r) * tail;                     % -> result.n0.accuracy
perigee = (a0r * (1 - e) - 1) * Lkm;                         % perigee, km
period  = sym(2) * sym(pi) / n0;                             % orbital period

% {label, expr, expected Ts-power, expected Lkm-power}
R = {
  {'a1   = (ke/n)^(2/3)        ~ Ts^0', a1,   0, 0}, ...
  {'del1 (first J2 correction) ~ Ts^0', del1, 0, 0}, ...
  {'a0   corrected (step 3)    ~ Ts^0', a0c,  0, 0}, ...
  {'del0 (refined correction)  ~ Ts^0', del0, 0, 0}, ...
  {'n0   recovered mean motion ~ Ts^-1 (rate)', n0, -1, 0}, ...
  {'a0   recovered (step 5)    ~ Ts^0', a0r,  0, 0}, ...
  {'tail O(del1^4) bound       ~ Ts^0', tail, 0, 0}, ...
  {'a0.accuracy += tail        ~ Ts^0', a0_acc, 0, 0}, ...
  {'n0.accuracy += (3/2)(n0/a0)tail ~ Ts^-1', n0_acc, -1, 0}, ...
  {'perigee_km                 ~ Ts^0 Lkm^1', perigee, 0, 1}, ...
  {'period_min = 2pi/n0        ~ Ts^+1 (time)', period, 1, 0} ...
};

% fixed dimensionless probe values: [ke0 n0s k2s e cosi]
V = [Ts Lkm ke0 n0s k2s e cosi];
v_ref = [1 1 1 1 1 sym(1)/10 sym(1)/2];   % Ts=1, Lkm=1
v_Ts2 = [2 1 1 1 1 sym(1)/10 sym(1)/2];   % Ts=2, Lkm=1
v_Lk2 = [1 2 1 1 1 sym(1)/10 sym(1)/2];   % Ts=1, Lkm=2

np = 0; nf = 0;
printf('=== W15 dimensional audit: mean-element recovery ===\n');
for k = 1:numel(R)
  lab = R{k}{1}; ex = R{k}{2}; eT = R{k}{3}; eL = R{k}{4};
  q0 = double(subs(ex, V, v_ref));
  qT = double(subs(ex, V, v_Ts2));
  qL = double(subs(ex, V, v_Lk2));
  Traw = log(qT / q0) / log(2);
  Lraw = log(qL / q0) / log(2);
  Tp = round(Traw); Lp = round(Lraw);
  monomial = (abs(Traw - Tp) < 1e-6) && (abs(Lraw - Lp) < 1e-6);
  ok = monomial && (Tp == eT) && (Lp == eL);
  if ok
    np = np + 1; st = 'PASS';
  else
    nf = nf + 1; st = 'FAIL';
  end
  printf('  [%s] %s   (got Ts^%d Lkm^%d)\n', st, lab, Tp, Lp);
end

% Kepler consistency: a1 * aE = (GM / n^2)^(1/3) with ke = sqrt(GM/aE^3).
% Numeric spot-check (fractional-power simplify is unreliable/slow); the
% relative residual must vanish to working precision.
printf('\n  -- Kepler-3rd-law consistency of a1 --\n');
syms GM aE nph positive
ke_phys = sqrt(GM / aE^3);
lhs = (ke_phys / nph)^(sym(2)/3) * aE;
rhs = (GM / nph^2)^(sym(1)/3);
subsvars = [GM aE nph];
subsvals = [sym(3986004418)/10, sym(6378137), sym(11)/10000];
L = double(subs(lhs, subsvars, subsvals));
Rk = double(subs(rhs, subsvars, subsvals));
kep_ok = (abs(L - Rk) <= 1e-6 * abs(Rk));
if kep_ok
  np = np + 1; printf('  [PASS] a1*aE = (GM/n^2)^(1/3)  (a1 is the semi-major axis; rel.resid %.1e)\n', abs(L-Rk)/abs(Rk));
else
  nf = nf + 1; printf('  [FAIL] a1*aE = (GM/n^2)^(1/3)  (lhs=%.6e rhs=%.6e)\n', L, Rk);
end

printf('\n');
if nf == 0
  printf('verify_dim_recovery: ALL %d checks PASS\n', np);
  exit(0);
else
  printf('verify_dim_recovery: %d PASS, %d FAIL\n', np, nf);
  exit(1);
end
