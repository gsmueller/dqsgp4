% vanishes_eta_only(expr, e, eta) — returns true iff expr == 0 in Q[e, eta]/(eta^2 - (1-e^2)),
% handling odd-eta residues via the even/odd-eta decomposition.
%
% Use when the expression has both even and odd powers of eta.  (The
% simpler vanishes_simple also works but this is faster when odd-eta
% is structural.)
function tf = vanishes_eta_only(expr, e, eta)
  cleared = expand(expr * e^10 * (1 + eta)^4 * eta^4);
  for it = 1:25
    prev = cleared;
    cleared = expand(subs(cleared, eta^2, 1 - e^2));
    if isequal(cleared, prev), break; end
  end
  cleared = simplify(cleared);
  a_part = simplify(subs(cleared, eta, 0));
  b_part = simplify(subs(diff(cleared, eta), eta, 0));
  tf = isequal(a_part, sym(0)) && isequal(b_part, sym(0));
end
