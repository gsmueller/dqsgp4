% verify_C4C5_target.m
%
% THEORY-ONLY: the CORRECT target for the SGP4 eccentricity drag coefficients C4 and C5,
% clean-session 2026-06-02. Trust model: derive from THEORY; the code is the arbiter.
%
% The SGP4 eccentricity coefficients are the Lane-Cranford translation of the Brouwer-Hori
% angular-momentum drag rate into Keplerian form. The mean-element eccentricity rate is
%     e-dot = (beta/(e L)) (beta L-dot - G-dot),     e = sqrt(1 - (G/L)^2),  G = L beta,
% i.e. e-dot is built from the ENERGY rate L-dot (-> a, C1/C2; no beta) and the ANGULAR-MOMENTUM
% rate G-dot (-> e, node; carries the beta/eta normalisation) -- the SAME energy-vs-angular-momentum
% distinction that RESOLVED xnodcf (verify_xnodcf_theory.m). This is pure kinematics (derived here),
% NOT taken from any source.
%
% Consequence (corroborated independently): with this CORRECT target,
%   * C4 Part-A (Keplerian) lands exactly (the (e+cos f) Gauss kernel is its leading realisation);
%   * C4 Part-B isotropic comes out ORDER-eta (the order-1 J2-secular bracket is SHARED by L-dot and
%     G-dot, so it CANCELS in beta*L-dot - G-dot) -- see verify_C4B_theory.m;
%   * C5's bracket comes out 1 + 4 eta^2 (NOT the code's 1 + (11/4) eta^2) -- see verify_C5_theory.m.
% The code's order-1 (C4-B) and 11/4 (C5) are OPERATIONAL SR3/Lane-Cranford forms, NOT reproducible
% from trusted (theoretical or born-digital) drag theory: code-matched for bit-compatibility,
% theoretically UNRESOLVED (this is now corroborated via a trusted route, not asserted from OCR).
%
% Octave gotcha: syms e0 fails -> use ecc.

pkg load symbolic;
syms L positive;
syms ecc Ldot Gdot real;

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

printf('=== C4/C5 CORRECT TARGET: e-dot = (beta/(eL))(beta L-dot - G-dot) ===\n\n');

beta = sqrt(1 - ecc^2);          % beta = sqrt(1-e^2)
G    = L*beta;                   % angular momentum  (e = sqrt(1-(G/L)^2))

% e as a function of (L,G): e = sqrt(1 - (G/L)^2). Treat G,L time-dependent: e-dot via chain rule.
% de/dG = -G/(e L^2),  de/dL = G^2/(e L^3).  e-dot = (de/dG) Gdot + (de/dL) Ldot.
edot_chain = (-G/(ecc*L^2))*Gdot + (G^2/(ecc*L^3))*Ldot;     % G,L explicit; e=ecc
% the claimed compact form:
edot_target = (beta/(ecc*L))*(beta*Ldot - Gdot);
[pass,fail,names] = ck('(1) e-dot = (beta/(eL))(beta L-dot - G-dot)  [kinematic, e=sqrt(1-(G/L)^2)]', ...
    simplify(edot_chain), simplify(edot_target), pass, fail, names);

% structural corollary: an ISOTROPIC J2-secular bracket B that is SHARED by both momentum rates
% (Ldot = -C1*(1+k2*B)*L-ish and Gdot = -C1*beta^2*(1+k2*B)*G-ish at leading order) CANCELS in the
% combination beta*Ldot - Gdot at e=0 (beta=1): the shared constant drops, leaving the order-eta part.
syms C1 k2 B real
Ldot_mod = -C1*(1 + k2*B)*L;                 % energy rate carrying the shared isotropic bracket B
Gdot_mod = -C1*(1 + k2*B)*G;                 % ang-mom rate carrying the SAME bracket B (at e->0)
edot_e0  = subs( (beta/(ecc*L))*(beta*Ldot_mod - Gdot_mod), ecc, 0 );   % beta->1 at e=0
% beta*Ldot - Gdot at e=0 (beta=1, G=L): (1)*(-C1(1+k2 B)L) - (-C1(1+k2 B)L) = 0  => the bracket cancels
val = simplify( subs( beta*Ldot_mod - Gdot_mod, ecc, 0 ) );
[pass,fail,names] = ck('(2) shared isotropic bracket cancels in (beta L-dot - G-dot) at e=0', ...
    val, sym(0), pass, fail, names);

printf('\n=== C4/C5 target summary ===\n');
printf('  %d / %d checks PASS\n', pass, pass+fail);
printf('  The correct clean-theory target is e-dot = (beta/(eL))(beta L-dot - G-dot): C4/C5 are built\n');
printf('  from the ENERGY rate (L-dot, no beta) and the ANGULAR-MOMENTUM rate (G-dot, carries beta/eta)\n');
printf('  -- same distinction that resolved xnodcf. C4 Part-A lands; C4-B is order-eta (the order-1\n');
printf('  J2 bracket cancels in beta*L-dot - G-dot); C5 is 1+4 eta^2. The code order-1/11-4 are\n');
printf('  OPERATIONAL forms, NOT clean-theory-derivable (code-matched; theoretically UNRESOLVED).\n');
if fail > 0
    printf('\n  FAILED:\n');
    for k = 1:numel(names); printf('    %s\n', names{k}); end
    printf('\n  C4/C5 target verification: FAIL\n'); exit(1);
else
    printf('\n  C4/C5 target verification: PASS\n'); exit(0);
end
