1;
% verify_dq_algebra.m
%
% W19 symbolic verification of the REQ-DQ algebraic laws
% (design/specifications/dual_quaternion_algebra.md) for the implementations in
% src/math/dual_number.h and src/math/quaternion.h. Where the numeric tests
% (test_dual_number, test_quaternion, test_dual_quaternion) SAMPLE these laws,
% this proves the closed-form / polynomial ones EXACTLY: expand(LHS - RHS) is
% identically zero in the symbolic components.
%
% Conventions matched to src/math/quaternion.h:
%   Hamilton product  (a.w b.w - a.x b.x - a.y b.y - a.z b.z, ...) [lines 120-124]
%   conjugate         q* = (w, -x, -y, -z)
%   magnitude^2       w^2 + x^2 + y^2 + z^2
%   rotation          R(q) v = vector part of q (0,v) q*  (active, v' = q v q*)
% and src/math/dual_number.h:  (a0 + a1 e)(b0 + b1 e) = a0 b0 + (a0 b1 + a1 b0) e.
%
% REQ-DQ-12/13/14/16/18 involve exp/log/normalize (transcendental, shortest-path
% branch) and are exercised by the numeric sample tests; this script proves the
% polynomial laws REQ-DQ-1..11 and REQ-DQ-15 (dual quaternion composition).
%
% Method: anonymous fns (callable from the script body, unlike Octave
% script-local functions); each law's LHS-RHS is expanded and tested == 0.
% exit(0) iff every law verifies clean.

pkg load symbolic

syms w1 x1 y1 z1 w2 x2 y2 z2 w3 x3 y3 z3 real
syms vx vy vz xv real
syms a0 a1 b0 b1 c0 c1 real

% --- quaternion algebra (4-vectors [w;x;y;z]); conventions per quaternion.h ---
qmul = @(a,b) [ a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4);
                a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3);
                a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2);
                a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1) ];
qconj = @(a) [a(1); -a(2); -a(3); -a(4)];
qmag2 = @(a) a(1)^2 + a(2)^2 + a(3)^2 + a(4)^2;
qrot  = @(q,V) qmul(qmul(q, V), qconj(q));         % q (0,v) q*

q1 = [w1; x1; y1; z1];
q2 = [w2; x2; y2; z2];
q3 = [w3; x3; y3; z3];
ident = [sym(1); sym(0); sym(0); sym(0)];
V = [sym(0); vx; vy; vz];

% --- dual numbers [value; derivative]; convention per dual_number.h ---
dmul = @(a,b) [ a(1)*b(1); a(1)*b(2) + a(2)*b(1) ];
da = [a0; a1]; db = [b0; b1]; dc = [c0; c1];
eps_dn = [sym(0); sym(1)];
dx = [xv; sym(1)];                                  % autodiff seed (x, 1)

% --- precomputed rotation intermediates ---
rot1     = qrot(q1, V);                              % R(q1) v (unnormalized)
rotc_lhs = qrot(qmul(q1, q2), V);                    % R(q1 q2) v
rotc_rhs = qrot(q1, qrot(q2, V));                    % R(q1) (R(q2) v)

% --- dual quaternions (real + eps dual); 24 independent symbols ---
ar = sym('ar', [4 1]); ad = sym('ad', [4 1]);
br = sym('br', [4 1]); bd = sym('bd', [4 1]);
cr = sym('cr', [4 1]); cd = sym('cd', [4 1]);
% product: (ar + e ad)(br + e bd) = ar*br + e(ar*bd + ad*br)
dq_re = @(Pr,Pd,Qr,Qd) qmul(Pr, Qr);
dq_du = @(Pr,Pd,Qr,Qd) qmul(Pr, Qd) + qmul(Pd, Qr);
% (A B) C, real and dual parts
ab_re = dq_re(ar,ad,br,bd); ab_du = dq_du(ar,ad,br,bd);
abc_re_L = dq_re(ab_re,ab_du,cr,cd); abc_du_L = dq_du(ab_re,ab_du,cr,cd);
% A (B C)
bc_re = dq_re(br,bd,cr,cd); bc_du = dq_du(br,bd,cr,cd);
abc_re_R = dq_re(ar,ad,bc_re,bc_du); abc_du_R = dq_du(ar,ad,bc_re,bc_du);

% --- law table: {label, difference-expression-that-must-vanish} ---
laws = {
  {'REQ-DQ-1  e^2 = 0', dmul(eps_dn, eps_dn)}, ...
  {'REQ-DQ-2  dual mult associativity', dmul(dmul(da,db),dc) - dmul(da,dmul(db,dc))}, ...
  {'REQ-DQ-2  dual mult commutativity', dmul(da,db) - dmul(db,da)}, ...
  {'REQ-DQ-2  dual distributivity', dmul(da, db+dc) - (dmul(da,db) + dmul(da,dc))}, ...
  {'REQ-DQ-3  autodiff d/dx x^2 = 2x', dmul(dx,dx) - [xv^2; 2*xv]}, ...
  {'REQ-DQ-3  autodiff d/dx x^3 = 3x^2', dmul(dmul(dx,dx),dx) - [xv^3; 3*xv^2]}, ...
  {'REQ-DQ-4  Hamilton associativity', qmul(qmul(q1,q2),q3) - qmul(q1,qmul(q2,q3))}, ...
  {'REQ-DQ-5  identity laws', qmul(q1,ident) - q1}, ...
  {'REQ-DQ-6  conjugate involutive', qconj(qconj(q1)) - q1}, ...
  {'REQ-DQ-7  conjugate of product', qconj(qmul(q1,q2)) - qmul(qconj(q2),qconj(q1))}, ...
  {'REQ-DQ-8  magnitude multiplicative', qmag2(qmul(q1,q2)) - qmag2(q1)*qmag2(q2)}, ...
  {'REQ-DQ-9  q q* = |q|^2 (inverse = q*/|q|^2)', qmul(q1,qconj(q1)) - [qmag2(q1);sym(0);sym(0);sym(0)]}, ...
  {'REQ-DQ-10 rotation scales length by |q|^2', qmag2(rot1) - qmag2(q1)^2*(vx^2+vy^2+vz^2)}, ...
  {'REQ-DQ-10 rotation result is pure', rot1(1)}, ...
  {'REQ-DQ-11 rotation composition', rotc_lhs - rotc_rhs}, ...
  {'REQ-DQ-15 dual quaternion assoc (real part)', abc_re_L - abc_re_R}, ...
  {'REQ-DQ-15 dual quaternion assoc (dual part)', abc_du_L - abc_du_R} ...
};

np = 0; nf = 0;
printf('=== W19 symbolic REQ-DQ algebra verification ===\n');
for k = 1:numel(laws)
  lab = laws{k}{1}; d = laws{k}{2};
  d = expand(d);
  tf = false;
  try
    tf = all(double(d(:)) == 0);
  catch
    tf = false;
  end
  if tf
    np = np + 1; st = 'PASS';
  else
    nf = nf + 1; st = 'FAIL';
  end
  printf('  [%s] %s\n', st, lab);
end

printf('\n');
if nf == 0
  printf('verify_dq_algebra: ALL %d REQ-DQ laws verified symbolically\n', np);
  exit(0);
else
  printf('verify_dq_algebra: %d PASS, %d FAIL\n', np, nf);
  exit(1);
end
