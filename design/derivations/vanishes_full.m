% vanishes_full(expr, e, eta, c, s, c2g, s2g) — returns true iff expr == 0
% in the quotient ring
%   Q[e, eta, c, s, c2g, s2g] / { eta^2 - (1-e^2), s^2 - (1-c^2), c2g^2 - (1-s2g^2) }
% where c = cos E, s = sin E, c2g = cos(2g), s2g = sin(2g).
%
% Method: clear denominators, reduce iteratively, decompose into the
% 8 components of the quotient ring (eta^{0,1} x s^{0,1} x c2g^{0,1});
% all must vanish.
function tf = vanishes_full(expr, e, eta, c, s, c2g, s2g)
  cleared = expand(expr * e^15 * eta^10 * (1 + eta)^4 * (1 - e * c)^15);
  for it = 1:30
    prev = cleared;
    cleared = expand(subs(cleared, s^2, 1 - c^2));
    cleared = expand(subs(cleared, eta^2, 1 - e^2));
    cleared = expand(subs(cleared, c2g^2, 1 - s2g^2));
    if isequal(cleared, prev), break; end
  end
  tf = true;
  for ea = 0:1
    for sa = 0:1
      for ga = 0:1
        if ea == 0,  xe = subs(cleared, eta, 0);
        else,        xe = subs(diff(cleared, eta), eta, 0);  end
        if sa == 0,  xs = subs(xe, s, 0);
        else,        xs = subs(diff(xe, s), s, 0);  end
        if ga == 0,  xg = subs(xs, c2g, 0);
        else,        xg = subs(diff(xs, c2g), c2g, 0);  end
        coef = simplify(xg);
        if ~isequal(coef, sym(0))
          tf = false;
          fprintf('    component (eta^%d s^%d c2g^%d) = %s\n', ea, sa, ga, char(coef));
        end
      end
    end
  end
end
