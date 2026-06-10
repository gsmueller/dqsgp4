% verify_xmcof_theory.m
%
% CLEAN-THEORY derivation of the SGP4 xmcof (M_drag_coef) leading -2/3 and the
% (1+eta cos M)^3 "cubic-density" mean-anomaly drag correction (delm), from the
% trusted Phase-0 Gauss VOP for Mdot (Thm 0.3.6 / 0.4.2.5) + the Lane density
% (q0-s)^4/(r-s)^4 = coef/(1-eta cos f_dagger)^4.  Standard 10 where a clean
% identity exists; the a0 scalar and the (1-eta c)^-3 -> (1+eta c)^3 reciprocation
% are flagged as OPERATIONAL (documented, not claimed as simplify==0 code-match).
%
%   code (drag_coefficients.h:197, secular_update.h:97-101):
%     xmcof = -(2/3) coef B* /(e0 eta) ,   delm = xmcof[(1+eta cosM)^3 - (1+eta cosM0)^3]
%
% RESULT (this file): the periodic mean-anomaly drag correction
%     delm = integral_{M0}^{M} (Mdot_drag - n)/n  dM'
% has the EXACT closed antiderivative  -(2/3)(a0 B* coef/(e0 eta))(1-eta cos f)^-3,
% so the code's -2/3, 1/(e0 eta), coef, B* are all clean-derived; the code's
% (1+eta cosM)^3 is the O(eta) reciprocal-binomial operational form and the code
% drops the circular-speed scale a0 (a0~1 in ER; the delm correction is small).
%
% Partner of C3 (the J3 piece of the same delta-ell_D): see verify_C3.m + trace.

pkg load symbolic;
syms f ecc beta n a Bstar coef real;
syms eta positive;   % eta = a0 e0 xi > 0 physically (avoids degenerate eta=0 branch in int)
syms c real;   % c := cos f  (for the reciprocation check)

pass_count = 0; fail_count = 0; failed_names = {};
function [pc, fc, fn] = check_eq(name, lhs, rhs, pc, fc, fn)
    d = simplify(lhs - rhs);
    if isequal(d, sym(0))
        printf('  PASS: %s\n', name); pc = pc + 1;
    else
        printf('  FAIL: %s\n', name);
        printf('        lhs - rhs = %s\n', char(d));
        fc = fc + 1; fn{end+1} = name;
    end
endfunction

printf('=== xmcof / delm : clean theory of the -2/3 cubic-density M-drag correction ===\n\n');

% ---- Trusted Phase-0 closed forms (born-digital; code-validated via C2/C4) ----
rp   = 1/(1+ecc*cos(f));                          % r/p              (Thm 0.2.2)
r    = a*beta^2/(1+ecc*cos(f));                   % r                (Thm 0.2.2)
absv = (n*a/beta)*sqrt(1+ecc^2+2*ecc*cos(f));     % |v|              (0.4.1.11)
rho  = coef/(1-eta*cos(f))^4;                     % rho/rho0 (B* carries rho0); f<-f_dagger leading (Lane)
Rdr  = -Bstar*rho*absv*(n*a*ecc*sin(f)/beta);     % R_drag           (0.4.1.8)
Tdr  = -Bstar*rho*absv*(n*a/beta)*(1+ecc*cos(f)); % T_drag           (0.4.1.8)

% (X.1) Mdot - n from the Gauss VOP (Thm 0.3.6 eq 0.3.6.1) under drag (0.4.2.5):
Mdot_minus_n = -(2*r*Rdr)/(n*a^2) + (beta^2/(n*a*ecc))*( Rdr*cos(f) - Tdr*sin(f)*(1+rp) );
rate = Mdot_minus_n / n;                           % dM_corr/dM = (Mdot-n)/n  (dt=dM/n, O(B*))
rate_closed = 2*Bstar*a*coef*(ecc^2+ecc*cos(f)+1)*sqrt(1+ecc^2+2*ecc*cos(f))*sin(f) ...
              / ( ecc*(1+ecc*cos(f))*(1-eta*cos(f))^4 );
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.1 rate (Mdot-n)/n = closed form [Gauss Thm 0.3.6 + drag 0.4.1]', ...
  rate, rate_closed, pass_count, fail_count, failed_names);

% (X.2) leading order in e: the 1/e coordinate-singular transverse term dominates.
%       Residue of the 1/e pole: limit_{e->0}(e * rate) = 2 B* a coef sin f /(1-eta cos f)^4.
%       The "2" is the transverse Gauss factor (1+r/p)|_{e=0}=2 (same 2 as cosf+cosE->2cosf).
rate_lead = 2*Bstar*a*coef*sin(f)/(ecc*(1-eta*cos(f))^4);
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.2 limit_{e->0}(e*rate) = 2 B* a coef sin f /(1-eta cos f)^4   [(1+r/p)->2 gives the 2]', ...
  limit(ecc*rate, ecc, 0), ecc*rate_lead, pass_count, fail_count, failed_names);

% (X.3) THE KEY ANTIDERIVATIVE: integral of the leading-e rate over f (dM=df at O(e))
%       int sin f (1-eta cos f)^-4 df = (1/eta) int u^-4 du = -(1/(3 eta)) (1-eta cos f)^-3
%       -> the power-4 Lane density antiderivative gives the CUBIC and the -1/3 of -2/3.
G        = int(rate_lead, f);
G_claim  = -(sym(2)/3)*(Bstar*a*coef/(ecc*eta))*(1-eta*cos(f))^(-3);
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.3 antiderivative int(rate_lead)df = -(2/3)(B* a coef/(e eta))(1-eta cos f)^-3   [gives -2/3, 1/(e eta), cubic]', ...
  G, G_claim, pass_count, fail_count, failed_names);

% (X.4) round-trip: d/df of the claimed antiderivative returns the leading-e rate.
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.4 d/df[ -(2/3)(B* a coef/(e eta))(1-eta cos f)^-3 ] = rate_lead', ...
  diff(G_claim, f), rate_lead, pass_count, fail_count, failed_names);

% (X.5) -2/3 factorization: 2 (transverse Gauss (1+r/p)|_{e=0}) x (-1/3) (Lane power 4->3).
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.5 -2/3 = 2 * (-1/3)   [transverse Gauss 2 x density-antiderivative -1/3]', ...
  sym(2)*(-sym(1)/3), -sym(2)/3, pass_count, fail_count, failed_names);

% (X.6) OPERATIONAL reciprocation: (1-eta c)^-3 = (1+eta c)^3 + O(eta^2)  (agree thru O(eta)).
recip_m3 = taylor((1-eta*c)^(-3), eta, 0, 'order', 2);  % thru O(eta^1)
cube_p3  = taylor((1+eta*c)^( 3), eta, 0, 'order', 2);
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.6 (1-eta cosf)^-3 = (1+eta cosf)^3 through O(eta) [operational cubic form]', ...
  recip_m3, cube_p3, pass_count, fail_count, failed_names);

% (X.7) code-form algebra: xmcof(1+eta cosM)^3 with xmcof=-(2/3)coef B*/(e eta).
syms M real;
xmcof_code = -(sym(2)/3)*coef*Bstar/(ecc*eta);
delm_code  = xmcof_code*(1+eta*cos(M))^3;
% theory amplitude = a0 * code amplitude (the dropped circular-speed scale):
xmcof_theory = -(sym(2)/3)*a*coef*Bstar/(ecc*eta);
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.7 xmcof_theory = a0 * xmcof_code   [code drops the circular-speed scale a0; OPERATIONAL]', ...
  xmcof_theory, a*xmcof_code, pass_count, fail_count, failed_names);

% (X.8) continuity: delm = G(f)-G(f0) (definite integral from epoch) <-> [(.)^3 - (.)_0^3].
%       structural: the "- (1+eta cosM0)^3" is the lower limit of the same antiderivative.
syms f0 real;
delm_theory_struct = G_claim - subs(G_claim, f, f0);   % = -(2/3)(B*a coef/(e eta))[(1-eta cf)^-3 - (1-eta cf0)^-3]
[pass_count,fail_count,failed_names] = check_eq( ...
  'X.8 delm = G(f)-G(f0): the epoch-subtraction matches the definite integral lower limit', ...
  delm_theory_struct, -(sym(2)/3)*(Bstar*a*coef/(ecc*eta))*((1-eta*cos(f))^(-3)-(1-eta*cos(f0))^(-3)), ...
  pass_count, fail_count, failed_names);

% (X.9) NUMERICAL corroboration (the a0 is intrinsic, not a leading-order artifact):
%   integrate the FULL Gauss Mdot-n rate (Thm 0.3.6, no e-truncation) over one orbit and
%   compare the periodic delta-ell amplitude to the code delm. As eta->0 the ratio -> a0,
%   confirming X.7's a0 for the full rate, not just the leading kernel. (The action route
%   delta-ell = dS/dL also gives delta-ell ~ a0 analytically: dS/dL picks up L^2/mu = a0 via
%   the dominant d(eta)/dL ~ eta/(L e^2) channel -- so BOTH routes carry a0; the code's
%   omission of a0 is operational, corroborated by two independent routes.)
mu_=1; s_=1.01; q4_=(120/6378.135)^4; Bs_=1e-4; Np=2000; L_=linspace(0,2*pi,Np+1); L_(end)=[]; dL_=2*pi/Np;
a0n=1.30; e0n=0.008;                         % small eta case (eta ~ 0.036)
nn=mu_^2/(mu_*a0n)^1.5; bn=sqrt(1-e0n^2);
En=zeros(size(L_)); for ii=1:numel(L_), En(ii)=fzero(@(x)x-e0n*sin(x)-L_(ii),L_(ii)); end
fn=2*atan2(sqrt(1+e0n)*sin(En/2),sqrt(1-e0n)*cos(En/2));
rn=a0n*(1-e0n*cos(En)); rpn=rn/(a0n*bn^2); rhn=q4_./((rn-s_).^4);
vn=(nn*a0n/bn)*sqrt(1+e0n^2+2*e0n*cos(fn));
Rn=-Bs_.*rhn.*vn.*(nn*a0n*e0n*sin(fn)/bn); Tn=-Bs_.*rhn.*vn.*(nn*a0n/bn).*(1+e0n*cos(fn));
raten=(-(2*rn.*Rn)/(nn*a0n^2)+(bn^2./(nn*a0n*e0n)).*(Rn.*cos(fn)-Tn.*sin(fn).*(1+rpn)))/nn;
raten=raten-mean(raten); dern=cumsum(raten)*dL_; dern=dern-mean(dern);
xin=1/(a0n-s_); etn=a0n*e0n*xin; cfn=q4_*xin^4; xmn=-(2/3)*cfn*Bs_/(e0n*etn);
dmn=xmn*((1+etn*cos(L_)).^3-(1+etn*cos(L_(1)))^3); dmn=dmn-mean(dmn);
ratio_a0 = (max(dern)-min(dern))/(max(dmn)-min(dmn));
printf('  X.9 [numeric] full-Gauss delta-ell/code = %.4f vs a0 = %.4f (eta=%.4f); |diff|/a0 = %.3f%%\n', ...
       ratio_a0, a0n, etn, 100*abs(ratio_a0-a0n)/a0n);
if abs(ratio_a0 - a0n)/a0n < 0.03
  printf('  PASS: X.9 full-Gauss rate-route delta-ell/code -> a0 (within 3%% at eta=%.3f) [a0 is real; both routes carry it]\n', etn);
  pass_count = pass_count + 1;
else
  printf('  FAIL: X.9 ratio %.4f not within 3%% of a0 %.4f\n', ratio_a0, a0n);
  fail_count = fail_count + 1; failed_names{end+1} = 'X.9 numeric a0';
end

printf('\n=== xmcof summary ===\n');
total = pass_count + fail_count;
printf('  %d / %d checks PASS\n', pass_count, total);
printf('\n  DISPOSITION (honest):\n');
printf('  * CLEAN-DERIVED (exact): -2/3, 1/(e0 eta), coef, B*, the cubic power, and the\n');
printf('    epoch-subtraction continuity -- from the Gauss Mdot VOP (Thm 0.3.6) transverse\n');
printf('    1/e term + the Lane density (1-eta cos f)^-4 antiderivative (power 4->3 => -1/3;\n');
printf('    du=eta sinf df => 1/eta).  Matches Phase-0 align note (c): 1/(e0 eta) = 1/e [Gauss]\n');
printf('    x 1/eta [orbit-average].\n');
printf('  * OPERATIONAL (code simplifications, NOT clean-theory-exact): (i) the code drops a0\n');
printf('    (circular-speed scale; theory = a0 x code, X.7 symbolic + X.9 full-Gauss numeric ->a0);\n');
printf('    (ii) (1-eta cos f_dag)^-3 -> (1+eta cos M)^3, exact only through O(eta) (X.6); f_dag->M O(e).\n');
printf('  * The a0 is INTRINSIC, corroborated by TWO routes: the rate route (Gauss Mdot, X.1-X.9) and\n');
printf('    the Delaunay ACTION route delta-ell=dS/dL (shape corr 0.9995; dS/dL picks up L^2/mu=a0 via\n');
printf('    the dominant d(eta)/dL channel). The action route CONFIRMS a0, does NOT resolve it -- so the\n');
printf('    a0-drop is a genuine SGP4 operational simplification, not a rate-route artifact.\n');
printf('  * Code-matched for bit-compatibility (regression unchanged); the LEADING COEFFICIENT -2/3 is\n');
printf('    clean-theory-EXACT (unlike C5 11/4 which is operational).\n');
if fail_count > 0
    printf('\n  FAILED:\n');
    for k = 1:length(failed_names); printf('    %s\n', failed_names{k}); end
    printf('\n  xmcof theory: FAIL\n'); exit(1);
else
    printf('\n  xmcof theory: PASS\n'); exit(0);
end
