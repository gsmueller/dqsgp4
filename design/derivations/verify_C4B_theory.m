% verify_C4B_theory.m
%
% THEORY-ONLY analysis of the SGP4 C4 "Part B" (drag x J2) isotropic term, clean-session
% 2026-06-02. Trust model: derive from THEORY; the code is the arbiter.
% Companion to verify_C4.m (which holds the bare algebra code-match to SR3 p.11).
%
% PURPOSE: show RIGOROUSLY that the SAME density-gradient mechanism that produced C2 Part B
% (sgp4_drag_phase2a_C2_trace.md sec.A.6) gives, for the ecc-rate kernel h_C4, an isotropic
% (3cos^2 i - 1) term that is ORDER-eta (vanishes at e=eta=0), whereas the code's isotropic
% term is ORDER-1 (leads with the constant 1). Characterises the OPEN discrepancy precisely.
%
% C2 mechanism (worked): rho(r_kep+dr) = rho(1 - 4 dr/(r-s)), with dr/r = -(3/2)(k2/p^2)beta con41,
%   raising the Lane density power -4 -> -5; C2's Part-A integrand is (1-eta cos f)^{-4} (ORDER-1)
%   so C2 Part B ~ <(1-eta cos f)^{-5}> = I0(5) = ORDER-1  => code's (8+24eta^2+3eta^4).  GOOD.
%
% C4 mechanism (same): C4's Part-A integrand is h_C4 = h0 + e h1 with h0 = cos f (1-eta cos f)^{-4}
%   (ORDER-eta after averaging, because of the cos f). The density gradient raises -4->-5:
%   C4 Part-B isotropic leading ~ < cos f (1-eta cos f)^{-5} > = (I0(5)-I0(4))/eta = ORDER-eta.
%   => leading CONSTANT is ABSENT; my isotropic term is order-eta, the code is order-1.
%
% Physical confirmation (e=0): with only the SECULAR J2 perturbations, at osculating e=0 the
% radial rate r-dot = 0 (delta-r-dot_secular = 0) so R = 0, and rho,|v|,r-f-dot are all CONSTANT
% in f; hence e-dot = (2/na) T cos f = const*cos f and <e-dot> = const*<cos f> = 0. The isotropic
% drag x J2 ecc-rate VANISHES at e=0 => order-eta, NOT order-1.

pkg load symbolic;
syms eta ee real;          % eta = Lane density param; ee = explicit eccentricity

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

b2 = 1 - eta^2;                          % psi^2 = 1 - eta^2
function v = lane_I0(m, eta, b2)
    switch m
        case 2, v = 1/b2^(sym(3)/2);
        case 3, v = (2+eta^2)/(2*b2^(sym(5)/2));
        case 4, v = (2+3*eta^2)/(2*b2^(sym(7)/2));
        case 5, v = (8+24*eta^2+3*eta^4)/(8*b2^(sym(9)/2));
        otherwise, error('I0(%d) not tabulated', m);
    end
end
I0 = @(m) lane_I0(m, eta, b2);
% < cos^n f (1-eta cos f)^{-m} > via cos f = (1-(1-eta cos f))/eta
function v = avg_cosn_over(n, m, eta, I0f)
    v = sym(0);
    for k = 0:n
        v = v + nchoosek(sym(n), sym(k)) * (-1)^k * I0f(m - k);
    end
    v = simplify(v / eta^n);
end
avg_cn = @(n, m) avg_cosn_over(n, m, eta, @(mm) lane_I0(mm, eta, b2));

printf('=== C4 Part B THEORY: density-gradient mechanism gives ORDER-eta (clean-session) ===\n\n');

% --- (1) C2 control: density mechanism on C2''s order-1 integrand gives the code polynomial ---
% C2 Part B isotropic numerator is I0(5)''s numerator (8+24eta^2+3eta^4), an ORDER-1 polynomial.
C2B_iso = simplify( 8 * b2^(sym(9)/2) * I0(5) );      % strip the 1/(8 psi^9) -> numerator
[pass,fail,names] = ck('(1) C2 Part B isotropic = (8+24eta^2+3eta^4) [ORDER-1, leads with 8]', ...
    C2B_iso, 8 + 24*eta^2 + 3*eta^4, pass, fail, names);
[pass,fail,names] = ck('(1b) C2 Part B isotropic at eta=0 is 8 (NONZERO => order-1)', ...
    subs(C2B_iso, eta, 0), sym(8), pass, fail, names);

% --- (2) C4: same mechanism on h0 = cos f (1-eta cos f)^{-4} raises power -4 -> -5 ---
% leading (e^0, fixed-eta) isotropic kernel-average:  < cos f (1-eta cos f)^{-5} >
C4B_iso_lead = avg_cn(1, 5);                          % = (I0(5)-I0(4))/eta
[pass,fail,names] = ck('(2) C4 Part B isotropic leading = <cos f (1-eta cos f)^-5> = (I0(5)-I0(4))/eta', ...
    C4B_iso_lead, simplify((I0(5) - I0(4))/eta), pass, fail, names);

% it is ORDER-eta: strip psi^9 and show the numerator VANISHES at eta=0 (no constant term)
C4B_iso_num = simplify( 8 * b2^(sym(9)/2) * C4B_iso_lead );    % numerator over 1/(8 psi^9)
[pass,fail,names] = ck('(2b) C4 Part B isotropic numerator at eta=0 is 0 (=> ORDER-eta, NO constant)', ...
    subs(C4B_iso_num, eta, 0), sym(0), pass, fail, names);
printf('  [info] C4 Part B isotropic numerator (8 psi^9 x leading) = %s\n', char(C4B_iso_num));
printf('         -> factor of eta present => order-eta; leading constant (the code''s 1) is ABSENT.\n');

% --- (3) the code''s isotropic polynomial leads with 1 (ORDER-1) ---
% code isotropic bracket: (1 - 2 e eta + eta^2(3/2 - (1/2) e eta)); at e=eta=0 it is 1.
P_code = 1 - 2*ee*eta + eta^2*(sym(3)/2 - sym(1)/2*ee*eta);
[pass,fail,names] = ck('(3) code C4 Part B isotropic poly at e=eta=0 is 1 (NONZERO => order-1)', ...
    subs(subs(P_code, ee, 0), eta, 0), sym(1), pass, fail, names);

printf('\n=== C4 Part B theory summary ===\n');
printf('  %d / %d checks PASS\n', pass, pass+fail);
printf('  FINDING (trusted theory): the C2 density-gradient mechanism, applied to C4''s ecc-rate\n');
printf('    kernel h_C4 (whose Part-A integrand carries the cos f of the (e+cos f) Gauss kernel),\n');
printf('    gives an isotropic (3cos^2 i-1) term of ORDER-eta (vanishes at e=eta=0; confirmed by\n');
printf('    the e=0 physical argument: r-dot=0 => <e-dot>=0). The code''s isotropic term is\n');
printf('    ORDER-1 (leads with 1). order-eta is CORROBORATED (independent sympy + born-digital BH61):\n');
printf('    the order-1 J2 bracket cancels in beta*Ldot-Gdot (verify_C4C5_target.m). The code order-1 is\n');
printf('    an OPERATIONAL form, NOT clean-theory-derivable (code-matched; theoretically UNRESOLVED).\n');
if fail > 0
    printf('\n  FAILED:\n');
    for k = 1:numel(names); printf('    %s\n', names{k}); end
    printf('\n  C4 Part B theory verification: FAIL\n'); exit(1);
else
    printf('\n  C4 Part B theory verification: PASS  (order-eta vs order-1 pinned down)\n'); exit(0);
end
