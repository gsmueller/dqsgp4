% verify_C5_theory.m
%
% THEORY-ONLY derivation of the SGP4 C5 periodic-eccentricity bracket, clean-session
% 2026-06-02. Trust model: derive from THEORY; the code is the arbiter.
% Companion to verify_C5.m (which holds the bare algebra code-match to SR3 p.11).
%
% PURPOSE: compute, RIGOROUSLY, what the "m=1 (cos M) Fourier harmonic of the Gauss
% eccentricity-rate kernel" gives for the C5 bracket, and CHARACTERISE the gap to the code.
%
% Setup (from sgp4_drag_C4_trace.md): the Keplerian drag ecc-rate kernel, expanded in the
% explicit eccentricity e at FIXED Lane eta, is
%     h_C4(f) = h0 + e*h1 + O(e^2),
%     h0 = cos f (1-eta cos f)^{-4},   h1 = (1+3 eta cos f) sin^2 f (1-eta cos f)^{-5}.
% C5's bracket = psi^7 * a1[h_C4], where a1[h] = (1/pi) Int_0^{2pi} h(f) cos M dM is the
% cos-M (mean-anomaly) Fourier coefficient. To O(e): cos M dM = (cos f - 2 e cos 2f) df, so
%     a1[h_C4] = A0 + e A1,   A0 = (1/pi)Int h0 cos f df,
%                             A1 = (1/pi)Int (h1 cos f - 2 h0 cos 2f) df.
%
% RESULT (this file proves):  psi^7 A0 = 1 + 4 eta^2   (my Gauss-kernel e^0 term),
%   but the code's eta-only part is 1 + (11/4) eta^2.  => a (5/4) eta^2 GAP at e^0.
%   So the from-scratch "cos-M Fourier of the Gauss kernel" does NOT close (documents the
%   OPEN status precisely, replacing the vague "does not close").
%
% Octave gotcha: syms e0 fails -> use ecc0 (here `e` is the explicit-eccentricity symbol -> `ee`).

pkg load symbolic;
syms eta ee real;          % eta = Lane density param; ee = explicit eccentricity (e)

pass = 0; fail = 0; names = {};
function [p, f, n] = ck(nm, l, r, p, f, n)
    d = simplify(l - r);
    if isequal(d, sym(0))
        printf('  PASS: %s\n', nm); p = p + 1;
    else
        printf('  FAIL: %s\n        l - r = %s\n', nm, char(d));
        f = f + 1; n{end+1} = nm;
    end
end

% --- Lane integrals  I0(m) = < (1-eta cos f)^{-m} >  (Phase 1; verify_phase1.m) ---
b2 = 1 - eta^2;                          % psi^2 = 1 - eta^2  (=beta_Lane^2)
I0 = @(m) lane_I0(m, eta, b2);
function v = lane_I0(m, eta, b2)
    switch m
        case 0, v = sym(1);
        case 1, v = 1/sqrt(b2);
        case 2, v = 1/b2^(sym(3)/2);
        case 3, v = (2+eta^2)/(2*b2^(sym(5)/2));
        case 4, v = (2+3*eta^2)/(2*b2^(sym(7)/2));
        case 5, v = (8+24*eta^2+3*eta^4)/(8*b2^(sym(9)/2));
        otherwise, error('I0(%d) not tabulated', m);
    end
end

% --- < cos^n f (1-eta cos f)^{-m} >  via  cos f = (1-(1-eta cos f))/eta ---
%     cos^n f = eta^{-n} sum_{k=0}^n C(n,k)(-1)^k (1-eta cos f)^k
%   => < cos^n f (1-eta cos f)^{-m} > = eta^{-n} sum_k C(n,k)(-1)^k I0(m-k).
avg_cn = @(n, m) avg_cosn_over(n, m, eta, @(mm) lane_I0(mm, eta, b2));
function v = avg_cosn_over(n, m, eta, I0f)
    v = sym(0);
    for k = 0:n
        v = v + nchoosek(sym(n), sym(k)) * (-1)^k * I0f(m - k);
    end
    v = simplify(v / eta^n);
end

printf('=== C5 THEORY: cos-M Fourier of the Gauss ecc-rate kernel (clean-session) ===\n\n');

% --- A0 = (1/pi) Int h0 cos f df = 2 < cos^2 f (1-eta cos f)^{-4} > ---
A0 = 2 * avg_cn(2, 4);
[pass,fail,names] = ck('(A0) (1/pi)Int h0 cos f df = (1+4 eta^2)/psi^7', ...
    A0, (1 + 4*eta^2)/b2^(sym(7)/2), pass, fail, names);

% e^0 term of the bracket:  psi^7 A0
br0_mine = simplify(b2^(sym(7)/2) * A0);
[pass,fail,names] = ck('(B0-mine) psi^7 A0 = 1 + 4 eta^2   (MY Gauss-kernel e^0 bracket)', ...
    br0_mine, 1 + 4*eta^2, pass, fail, names);

% code's eta-only part of the bracket is 1 + (11/4) eta^2:
br0_code = 1 + sym(11)/4*eta^2;
gap0 = simplify(br0_code - br0_mine);
[pass,fail,names] = ck('(GAP0) code(1+11/4 eta^2) - mine(1+4 eta^2) = -(5/4) eta^2', ...
    gap0, -sym(5)/4*eta^2, pass, fail, names);

% --- A1 = (1/pi) Int (h1 cos f - 2 h0 cos 2f) df  (the O(e) term) ---
% h1 cos f = (1+3 eta cos f) sin^2 f cos f (1-eta cos f)^{-5}
%          = (1+3 eta cos f)(cos f - cos^3 f)(1-eta cos f)^{-5}
% => (1/pi)Int = 2[ <cf(1-ec)^-5> - <cf^3(1-ec)^-5> + 3eta<cf^2(1-ec)^-5> - 3eta<cf^4(1-ec)^-5> ]
termA = 2*( avg_cn(1,5) - avg_cn(3,5) + 3*eta*avg_cn(2,5) - 3*eta*avg_cn(4,5) );
% 2 h0 cos 2f = 2 cos f (2cos^2 f - 1)(1-eta cos f)^{-4} = (4cos^3 f - 2 cos f)(1-eta cos f)^{-4}
% => (1/pi)Int 2 h0 cos2f df = 2*< (4cos^3 f - 2 cos f)(1-eta cos f)^{-4} >
termB = 2*( 4*avg_cn(3,4) - 2*avg_cn(1,4) );
A1 = simplify(termA - termB);

br1_mine = simplify(b2^(sym(7)/2) * A1);     % psi^7 A1 = the e^1 coefficient of MY bracket
% code's e^1 (explicit-e) coefficient: bracket = 1+(11/4)(eta^2+e*eta)+e*eta^3
%   => d(bracket)/d(e) = (11/4) eta + eta^3
br1_code = sym(11)/4*eta + eta^3;
gap1 = simplify(br1_code - br1_mine);
printf('  [info] psi^7 A1 (my e^1 bracket coeff)         = %s\n', char(br1_mine));
printf('  [info] code e^1 bracket coeff (11/4 eta+eta^3) = %s\n', char(br1_code));
printf('  [info] e^1 gap (code - mine)                   = %s\n', char(gap1));

printf('\n=== C5 theory summary ===\n');
printf('  %d / %d checks PASS\n', pass, pass+fail);
printf('  FINDING (trusted theory): the cos-M Fourier of the Gauss ecc-rate kernel gives\n');
printf('    e^0 bracket = 1 + 4 eta^2, but the code has 1 + (11/4) eta^2  => GAP = -(5/4) eta^2.\n');
printf('    CORROBORATED (independent sympy + born-digital BH61): the projection and kernel are correct;\n');
printf('    C5 is NOT the periodic partner of C4-A (verify_C4C5_target.m). The code 11/4 is an\n');
printf('    OPERATIONAL form, NOT clean-theory-derivable (code-matched; theoretically UNRESOLVED).\n');
if fail > 0
    printf('\n  FAILED:\n');
    for k = 1:numel(names); printf('    %s\n', names{k}); end
    printf('\n  C5 theory verification: FAIL\n'); exit(1);
else
    printf('\n  C5 theory verification: PASS  (gap precisely characterised)\n'); exit(0);
end
