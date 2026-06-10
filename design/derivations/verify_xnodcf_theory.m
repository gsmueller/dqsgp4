% verify_xnodcf_theory.m
%
% THEORY-ONLY derivation of the SGP4 secular node t^2 coefficient (nodecf / xnodcf),
% clean-session 2026-06-02. Trust model: derive from THEORY; the code is the arbiter.
% Companion to verify_corrections.m (which holds the bare algebra code-match).
%
% RESOLUTION (2026-06-02): the code's omeosq = beta0^2 factor is the ANGULAR-MOMENTUM drag
% normalization. My earlier hypothesis ("beta^2 unexplained; theory gives beta^-4 three ways")
% was REFUTED: I had wrongly forced Ldot/L = Gdot/G = -C1 (energy-like) onto ALL the Delaunay
% momenta. The correction (derived below from first principles, NOT taken from any source):
%   * energy / semimajor:   Ldot/L = (1/2) adot/a = -C1            (no beta-power; defines C1)
%   * angular momentum:     Gdot/G = -B* <rho|v|> = -C1 * beta^2   (carries omeosq, to O(e^2))
% so nodecf = (1/2)Omega_dot[-3 Ldot/L - 4 Gdot/G] = (1/2) C1 xhdot1 (3 + 4 beta^2), which
% CARRIES the beta^2 and matches the code at leading order. The code's exact (7/2) beta^2
% applies omeosq UNIFORMLY to the gradient sum (the SGP4 "beta=1 in numerators" simplification);
% it differs from the term-by-term theory only at O(e^2): code - theory = -(3/2)C1 xhdot1 e^2.
%
% (The structural hint "angular momentum carries beta^2" came from an independent agent; per the
%  trust model it is confirmed HERE by my own derivation, never taken as source fact.)
%
% Octave gotcha: syms e0 fails -> use ecc.

pkg load symbolic;
syms a0 mu positive;
syms n0 ecc k2 cos_i0 C1 xhdot1 E real;

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

beta0sq = 1 - ecc^2;                  % beta0^2 = 1 - e0^2  (Vallado omeosq)
n0_kep  = sqrt(mu) / a0^(sym(3)/2);   % Kepler n0

printf('=== xnodcf THEORY derivation (RESOLVED, clean-session 2026-06-02) ===\n\n');

% ============================================================================
% (1) J2 nodal rate Omega_dot ∝ H L^-3 G^-5 = code xhdot1   (justifies exponents)
% ============================================================================
xhdot1_code = -3*k2*n0*cos_i0/(a0^2*beta0sq^2);
xhdot1_kep  = subs(xhdot1_code, n0, n0_kep);
L = sqrt(mu*a0);  G = L*sqrt(beta0sq);  H = G*cos_i0;
OmDot_del = -3*k2*mu^4 * H / (L^3 * G^5);
[pass,fail,names] = ck('(1) Delaunay -3 k2 mu^4 H L^-3 G^-5 = code xhdot1 (a0^-7/2 beta0^-4)', ...
    simplify(OmDot_del), xhdot1_kep, pass, fail, names);

% ============================================================================
% (2) gradient structure: Omega_ddot = Omega_dot(-3 Ldot/L - 4 Gdot/G)  (i const => Hdot/H=Gdot/G)
% ============================================================================
syms Ldot Gdot Hdot real
OmDdot_expand = -3*k2*mu^4*( Hdot*L^-3*G^-5 + H*(-3)*L^-4*Ldot*G^-5 + H*L^-3*(-5)*G^-6*Gdot );
[pass,fail,names] = ck('(2) general: Omega_ddot/Omega_dot = Hdot/H - 3 Ldot/L - 5 Gdot/G', ...
    simplify(OmDdot_expand/OmDot_del), simplify(Hdot/H - 3*Ldot/L - 5*Gdot/G), pass, fail, names);
bracket_inplane = subs(Hdot/H - 3*Ldot/L - 5*Gdot/G, Hdot, H*Gdot/G);   % i const => Hdot/H = Gdot/G
[pass,fail,names] = ck('(2b) in-plane (Hdot/H=Gdot/G): bracket = -3 Ldot/L - 4 Gdot/G', ...
    simplify(bracket_inplane), -3*Ldot/L - 4*Gdot/G, pass, fail, names);

% ============================================================================
% (3) THE KEY: the angular-momentum drag rate carries beta^2 = omeosq (to O(e^2)).
%     Instantaneous: Gdot/G = -B* rho|v|  (torque hdot = r x F = -B* rho|v| h);
%                    Ldot/L = -(a/mu) B* rho|v|^3,  C1 := (a/mu) B* <rho|v|^3>.
%     vis-viva: |v|^2 = (mu/a)(1+u)/(1-u),  u = e cos E ;  dM = (1-u) dE.
%     With rho = rho0 (const), the rho0 sqrt(mu/a) common factor cancels in the ratio:
%       <rho|v|>            ~ < sqrt((1+u)/(1-u)) (1-u) > = < sqrt(1-u^2) >
%       (a/mu)<rho|v|^3>    ~ < (1+u)^{3/2} (1-u)^{-1/2} >
%     => <Gdot/G> / (-C1) = <sqrt(1-u^2)> / <(1+u)^{3/2}(1-u)^{-1/2}> = 1 - e^2 = beta^2.
% ============================================================================
u = ecc*cos(E);
avgE = @(expr) simplify( int( taylor(expr, ecc, 0, 'order', 3), E, 0, 2*sym(pi) ) / (2*sym(pi)) );
avg_G = avgE( sqrt(1 - u^2) );                                  % <rho|v|>  (stripped)
avg_L = avgE( (1+u)^(sym(3)/2) * (1-u)^(-sym(1)/2) );           % (a/mu)<rho|v|^3>  (stripped)
[pass,fail,names] = ck('(3a) <rho|v|> ~ 1 - e^2/4', avg_G, 1 - ecc^2/4, pass, fail, names);
[pass,fail,names] = ck('(3b) (a/mu)<rho|v|^3> ~ 1 + 3e^2/4', avg_L, 1 + 3*ecc^2/4, pass, fail, names);
ratio_GtoL = taylor(avg_G/avg_L, ecc, 0, 'order', 3);
[pass,fail,names] = ck('(3c) KEY: <Gdot/G>/(-C1) = <rho|v|>/((a/mu)<rho|v|^3>) = beta^2 = 1-e^2', ...
    ratio_GtoL, 1 - ecc^2, pass, fail, names);

% ============================================================================
% (4) corrected theory nodecf = (1/2) Omega_dot (-3 Ldot/L - 4 Gdot/G),
%     with Ldot/L = -C1 (energy) and Gdot/G = -C1 beta^2 (angular momentum):
%     nodecf_theory = (1/2) C1 xhdot1 (3 + 4 beta^2)   -- CARRIES beta^2.
% ============================================================================
nodecf_theory = sym(1)/2 * C1 * xhdot1 * (3 + 4*beta0sq);
nodecf_theory_chk = sym(1)/2 * xhdot1 * ( -3*(-C1) - 4*(-C1*beta0sq) );
[pass,fail,names] = ck('(4) nodecf_theory = (1/2) C1 xhdot1 (3 + 4 beta^2)', ...
    nodecf_theory, nodecf_theory_chk, pass, fail, names);

% ============================================================================
% (5) the code form (uniform omeosq) and the precise relation to the theory.
% ============================================================================
nodecf_code = sym(7)/2 * C1 * beta0sq * xhdot1;     % (7/2) beta0^2 xhdot1 C1 = -(21/2)n0 k2 c C1/(a^2 b^2)
[pass,fail,names] = ck('(5a) code = theory at e0 -> 0 (both = (7/2) C1 xhdot1)', ...
    subs(nodecf_code, ecc, 0), subs(nodecf_theory, ecc, 0), pass, fail, names);
[pass,fail,names] = ck('(5b) code - theory = -(3/2) C1 xhdot1 e^2  (uniform-omeosq SGP4 simplification)', ...
    simplify(nodecf_code - nodecf_theory), -sym(3)/2*C1*xhdot1*ecc^2, pass, fail, names);

% also: verify code's (7/2)beta^2 xhdot1 C1 equals the literal drag_coefficients.h:208-209 form
nodecf_code_literal = -sym(21)/2*n0_kep*k2*cos_i0*C1/(a0^2*beta0sq);
[pass,fail,names] = ck('(5c) (7/2) beta^2 xhdot1 C1 = -(21/2) n0 k2 cos_i0 C1/(a0^2 beta0^2)', ...
    simplify(subs(nodecf_code, xhdot1, xhdot1_kep)), simplify(nodecf_code_literal), pass, fail, names);

printf('\n=== xnodcf theory summary ===\n');
printf('  %d / %d checks PASS\n', pass, pass+fail);
printf('  RESOLUTION: the code''s beta0^2 (omeosq) is the ANGULAR-MOMENTUM drag normalization.\n');
printf('    Derived: <Gdot/G> = -B*<rho|v|> = -C1 beta^2 (to O(e^2)); Ldot/L = -C1 (no beta).\n');
printf('    => nodecf_theory = (1/2) C1 xhdot1 (3 + 4 beta^2)  CARRIES the beta^2 and matches the\n');
printf('       code at leading order. The code''s exact (7/2) beta^2 applies omeosq uniformly to the\n');
printf('       gradient sum (SGP4 beta=1-in-numerators simplification); code - theory = O(e^2).\n');
printf('    My earlier "beta^2 unexplained / beta^-4 three ways" hypothesis is REFUTED.\n');
if fail > 0
    printf('\n  FAILED:\n');
    for k = 1:numel(names); printf('    %s\n', names{k}); end
    printf('\n  xnodcf theory verification: FAIL\n'); exit(1);
else
    printf('\n  xnodcf theory verification: PASS\n'); exit(0);
end
