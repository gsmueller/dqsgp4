1;
% verify_dim_brouwer.m
%
% W13 dimensional audit of the Brouwer secular-rate boxed formulas in
% src/perturbation/brouwer.h (compute_secular_rates). Extends the
% drag-only dimensional discipline (feedback_dimensional_audit) to the
% Brouwer J2/J2^2/J4 secular rates.
%
% Claim audited: every additive term of M_dot, omega_dot, Omega_dot, and
% the J2^3 accuracy estimate has the dimension of a RATE, i.e. n^1 * L^0:
% exactly one factor of the mean motion n (dimension time^-1) and zero NET
% length dimension. Lengths enter only through the dimensionless ratios
% (aE/p), (aE/p)^2, ... once the Earth-radius unit aE — set to 1 in SGP4
% canonical units — is restored:
%     k2  = J2 * aE^2 / 2      (the code's CK2 = J2/2 with aE = 1)
%     J4 secular term carries aE^4
%     J2^3 accuracy term carries aE^6 (one extra (k2/p^2) beyond the J2^2 term)
%
% Method (matches verify_dimensional_audit_symbolic.m): all lengths share one
% symbol L (aE = a = p-scale = L); the L-exponent and n-exponent of each term
% are extracted by the log-derivative  v * d(log expr)/dv  and compared to the
% expected (L^0 or L^k, n^1 or n^0). simplify/diff/subs only; no local
% functions (Octave script-local functions are not reliably callable here).
%
% Octave symbolic. exit(0) iff every term audits clean.

pkg load symbolic

syms L e cosi J2 J4 n positive

% --- lengths share one dimension symbol L (Earth-radius unit restored) ---
aE = L;
a  = L;
p  = a * (1 - e^2);          % semi-latus rectum, dimension L
beta0  = sqrt(1 - e^2);      % dimensionless
x3thm1 = 3*cosi^2 - 1;       % dimensionless
x1m5th = 1 - 5*cosi^2;       % dimensionless

% --- base terms with aE restored (numerically aE = 1 in the code) ---
CK2   = J2 * aE^2 / 2;                 % k2  ~ L^2
temp1 = 3 * CK2 / p^2 * n;             % J2  base, a rate
temp2 = temp1 * CK2 / p^2;             % J2^2 base, a rate
temp3 = -15 * (J4 * aE^4) / (32 * p^4) * n;   % J4 base, a rate
xhdot1 = -temp1 * cosi;

% inclination polynomials (all dimensionless)
polyMj2 = 13 - 78*cosi^2 + 137*cosi^4;
polywj2 = 7  - 114*cosi^2 + 395*cosi^4;
polywj4 = 3  - 36*cosi^2  + 49*cosi^4;
polyOj2 = 4  - 19*cosi^2;
polyOj4 = 3  - 7*cosi^2;

% J2^3 accuracy scale: (k2/p^2)^3 ~ J2^3 aE^6 / p^6
acc_scale = J2^3 * aE^6 * n / p^6;
accM = (sym(1)/64) * polyMj2 * acc_scale;
accW = (sym(1)/64) * polywj2 * acc_scale;
accO = (sym(1)/64) * polyOj2 * cosi * acc_scale;

% --- test table: {label, expr, expected L-power, expected n-power} ---
% Intermediate quantities document the restored-aE dimensions; the rate
% terms are the boxed formulas (all must be L^0 n^1).
T = {
  {'p   (semi-latus rectum)            ~ L^1 n^0', p,    1, 0}, ...
  {'CK2 = J2 aE^2 / 2  (k2)            ~ L^2 n^0', CK2,  2, 0}, ...
  {'temp1 (J2 base)                    is a rate', temp1, 0, 1}, ...
  {'temp2 (J2^2 base)                  is a rate', temp2, 0, 1}, ...
  {'temp3 (J4 base)                    is a rate', temp3, 0, 1}, ...
  {'M_dot     term  n', n, 0, 1}, ...
  {'M_dot     term  (1/2) temp1 beta0 (3c^2-1)', (sym(1)/2)*temp1*beta0*x3thm1, 0, 1}, ...
  {'M_dot     term  (1/16) temp2 beta0 P_M', (sym(1)/16)*temp2*beta0*polyMj2, 0, 1}, ...
  {'omega_dot term  -(1/2) temp1 (1-5c^2)', -(sym(1)/2)*temp1*x1m5th, 0, 1}, ...
  {'omega_dot term  (1/16) temp2 P_w', (sym(1)/16)*temp2*polywj2, 0, 1}, ...
  {'omega_dot term  temp3 P_wJ4', temp3*polywj4, 0, 1}, ...
  {'Omega_dot term  xhdot1 = -temp1 cosi', xhdot1, 0, 1}, ...
  {'Omega_dot term  (1/2) temp2 (4-19c^2) cosi', (sym(1)/2)*temp2*polyOj2*cosi, 0, 1}, ...
  {'Omega_dot term  2 temp3 (3-7c^2) cosi', 2*temp3*polyOj4*cosi, 0, 1}, ...
  {'accuracy  M     (J2^3 carries aE^6)', accM, 0, 1}, ...
  {'accuracy  omega (J2^3 carries aE^6)', accW, 0, 1}, ...
  {'accuracy  Omega (J2^3 carries aE^6)', accO, 0, 1} ...
};

np = 0; nf = 0;
allvars = [e cosi J2 J4 n L];
ones6   = [1 1 1 1 1 1];

printf('=== W13 dimensional audit: Brouwer secular rates ===\n');
for k = 1:numel(T)
  lab = T{k}{1}; ex = T{k}{2}; eL = T{k}{3}; eN = T{k}{4};
  exs = simplify(ex);
  Lp = -99; Np = -99;
  if isequal(exs, sym(0))
    Lp = 0; Np = 0;
  else
    try
      ld = simplify(L * diff(exs, L) / exs);
      Lp = round(double(subs(ld, allvars, ones6)));
    catch
      Lp = -99;
    end
    try
      nd = simplify(n * diff(exs, n) / exs);
      Np = round(double(subs(nd, allvars, ones6)));
    catch
      Np = -99;
    end
  end
  ok = (Lp == eL) && (Np == eN);
  if ok
    np = np + 1; st = 'PASS';
  else
    nf = nf + 1; st = 'FAIL';
  end
  printf('  [%s] %s   (got L^%d n^%d)\n', st, lab, Lp, Np);
end

printf('\n');
if nf == 0
  printf('verify_dim_brouwer: ALL %d dimensional checks PASS\n', np);
  exit(0);
else
  printf('verify_dim_brouwer: %d PASS, %d FAIL\n', np, nf);
  exit(1);
end
