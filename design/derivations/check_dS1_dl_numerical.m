%% check_dS1_dl_numerical.m
%% Cross-check: dS1/dl (our chain) equals the form used in delta_p1.
%% Test points MUST satisfy G = L sqrt(1-e^2).

pkg load symbolic;

syms mu k2 real positive;
syms L G H real positive;
syms l_sym g e f real;

theta = H/G;
A_expr = -sym(1)/2 + sym(3)/2 * theta^2;
B_expr =  sym(3)/2 - sym(3)/2 * theta^2;

S1_core = A_expr*(f - l_sym + e*sin(f)) ...
        + B_expr*(sym(1)/2*sin(2*g+2*f) + (e/2)*sin(2*g+f) + (e/6)*sin(2*g+3*f));
S1 = (mu^2*k2/G^3) * S1_core;

df_dl = (1+e*cos(f))^2/(1-e^2)^(sym(3)/2);

dS1_dl = diff(S1, l_sym) + diff(S1, f) * df_dl;

a_over_r = (1 + e*cos(f))/(1-e^2);
expected_dS1_dl = (mu^2*k2/L^3) * ( A_expr*a_over_r^3 - A_expr*(L^3/G^3) ...
                                   + B_expr*a_over_r^3*cos(2*g+2*f) );

diff_expr = dS1_dl - expected_dS1_dl;

% Self-consistent random test points: G = L sqrt(1-e^2), |H| <= G.
function v = mktest(mu_v, k2_v, L_v, e_v, theta_v, f_v, g_v, l_v)
  v.mu = mu_v; v.k2 = k2_v; v.L = L_v; v.e = e_v;
  v.G = L_v * sqrt(1 - e_v^2);
  v.H = v.G * theta_v;
  v.f = f_v; v.g = g_v; v.l = l_v;
endfunction

tests = { mktest(1,    0.001, 1.2, 0.3, 0.8, 0.7, 1.5, 0.5), ...
          mktest(2,    0.005, 1.5, 0.1, 0.3, 2.1, 0.7, 1.2), ...
          mktest(3.98e5, -2.63e-3, 7000, 0.05, 0.6, 1.1, 2.3, 0.3) };

printf('Cross-check: (our dS1/dl) vs (expected form from delta_p1)\n');
printf('Test points are self-consistent (G = L*sqrt(1-e^2)).\n\n');
for i = 1:length(tests)
  v = tests{i};
  args = {mu, k2, L, G, H, e, f, g, l_sym};
  vals = {v.mu, v.k2, v.L, v.G, v.H, v.e, v.f, v.g, v.l};
  dv  = double(subs(diff_expr, args, vals));
  rv  = double(subs(dS1_dl,    args, vals));
  ev  = double(subs(expected_dS1_dl, args, vals));
  printf('  Test %d: ref=%12.6g  expected=%12.6g  diff=%12.6g  rel=%g\n', ...
    i, rv, ev, dv, abs(dv)/max(abs(rv),abs(ev)+eps));
end
