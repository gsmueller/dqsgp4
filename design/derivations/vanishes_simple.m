% vanishes_simple(expr, e, eta) — returns true iff expr == 0 in Q[e, eta]/(eta^2 - (1-e^2)).
%
% Method: substitute eta = sqrt(1 - e^2), then SymPy simplify.  Use for
% expressions in (e, eta) only (no trig).  Part of the verifier toolkit
% for the BH61 textbook derivation.
function tf = vanishes_simple(expr, e, eta)
  expr_e = subs(expr, eta, sqrt(1 - e^2));
  expr_e = simplify(expr_e);
  tf = isequal(expr_e, sym(0));
end
