% verify_ch05_F1_closed_form.m
%
% Symbolic + numerical verification of Theorem A.17 (Chapter 5 Part B, file ch05d):
%   F_1 = (mu k_2 / r^3) (A + B cos(2(f+g)))
% with A = (3 theta^2 - 1) / 2, B = 3 (1 - theta^2) / 2, theta = cos(I).
%
% Method: compute F_1 two independent ways and show their difference simplifies to zero.
%
% Chain 1 ("ours"):
%   Start from ch05b (0.B.7): F_1 = -(mu k_2 / r^3) (3 sin^2(beta) - 1).
%   Substitute ch05c (0.C.7): sin^2(beta) = sin^2(I) sin^2(f+g).
%   This gives F_1 directly in (r, I, f, g) without further manipulation.
%
% Chain 2 ("target"):
%   F_1 = (mu k_2 / r^3) (A + B cos(2(f+g))), with A, B as in (0.D.7).
%
% If Theorem A.17 is correct, simplify(F1_chain1 - F1_chain2) = 0 identically.
%
% NOTE on numerical substitution: Octave's symbolic package wraps numeric
% doubles via a `double_to_sym_heuristic` that performs rational approximations
% and can silently lose precision for non-exact values (e.g., 7e6 may get
% rationalized in a lossy way). To avoid this, we convert all numerical inputs
% to sym via exact-decimal strings, then force `double()` on the result.

pkg load symbolic;

syms r mu k2 I f g real;

disp('--- Chapter 5 Part B verification: F_1 closed form (Theorem A.17) ---');
disp('');

% Chain 1: from ch05b (0.B.7) + ch05c (0.C.7).
F1_chain1 = -(mu * k2 / r^3) * (3 * sin(I)^2 * sin(f+g)^2 - 1);
disp('F_1 via ch05b(0.B.7) + ch05c(0.C.7) [chain 1]:');
disp(F1_chain1);
disp('');

% Chain 2: Theorem A.17 (0.D.6) with A and B defined in (0.D.7).
theta = cos(I);
A = (3*theta^2 - 1) / 2;
B = 3 * (1 - theta^2) / 2;
F1_chain2 = (mu * k2 / r^3) * (A + B * cos(2*(f+g)));
disp('F_1 via Theorem A.17 (0.D.6) with A, B from (0.D.7) [chain 2]:');
disp(F1_chain2);
disp('');

% Symbolic difference.
diff_F1 = simplify(F1_chain1 - F1_chain2);
disp('simplify(F1_chain1 - F1_chain2) =');
disp(diff_F1);
disp('');

% Check invariant A + B = 1.
AB_sum = simplify(A + B);
disp('A + B simplifies to:');
disp(AB_sum);
disp('');

% Symbolic equality check.
is_zero = isequal(diff_F1, sym(0));
if is_zero
    disp('*** PASS: F1_chain1 - F1_chain2 simplifies to identically zero. ***');
else
    disp('*** FAIL: symbolic difference did NOT simplify to zero. ***');
    disp('Inspect the difference above.');
end
disp('');

% Helper: substitute with exact sym strings to avoid double_to_sym heuristic.
function v = eval_exact(expr, vars, vals_str)
    % vars: cell array of sym variables
    % vals_str: cell array of decimal strings
    vals_sym = cell(size(vals_str));
    for i = 1:numel(vals_str)
        vals_sym{i} = sym(vals_str{i});
    end
    v = double(subs(expr, vars, vals_sym));
end

% Numerical cross-check at 12 random points.
disp('Numerical check at 12 random (r, I, f, g, mu, k2) points (exact-sym subs):');
disp('--------------------------------------------------------------------------');
rng(42);
max_abs_diff = 0;
max_rel_diff = 0;
for t = 1:12
    r_n  = 1e7  * (0.5 + rand());
    I_n  = pi   * rand();
    f_n  = 2*pi * rand();
    g_n  = 2*pi * rand();
    mu_n = 3.986004418e14;
    k2_n = 0.5 * 1.0826267e-3 * (6.378137e6)^2;

    vars = {r, mu, k2, I, f, g};
    vals_str = {sprintf('%.17g', r_n),  sprintf('%.17g', mu_n), sprintf('%.17g', k2_n), ...
                sprintf('%.17g', I_n),  sprintf('%.17g', f_n),  sprintf('%.17g', g_n)};

    v1 = eval_exact(F1_chain1, vars, vals_str);
    v2 = eval_exact(F1_chain2, vars, vals_str);

    abs_diff = abs(v1 - v2);
    scale = max(abs(v1), abs(v2));
    rel_diff = abs_diff / max(scale, 1e-300);
    if abs_diff > max_abs_diff
        max_abs_diff = abs_diff;
    end
    if rel_diff > max_rel_diff
        max_rel_diff = rel_diff;
    end
    printf('Pt %2d: r=%.3e I=%.3f f=%.3f g=%.3f  v1=%+.6e  v2=%+.6e  |v1-v2|=%.3e  rel=%.3e\n', ...
           t, r_n, I_n, f_n, g_n, v1, v2, abs_diff, rel_diff);
end
disp('--------------------------------------------------------------------------');
printf('Max absolute difference across 12 points: %.3e\n', max_abs_diff);
printf('Max relative difference across 12 points: %.3e\n', max_rel_diff);
if max_rel_diff < 1e-14
    disp('*** PASS: numerical agreement to better than 1e-14 relative. ***');
else
    disp('*** WARN: numerical difference exceeds 1e-14 relative. Investigate. ***');
end
disp('');

% Sanity checks at special inclinations, using exact sym substitutions.
disp('Sanity at special inclinations (using F1_chain1 with exact-sym subs):');
disp('---------------------------------------------------------------------');
r_s_str  = '7000000';
mu_s_str = '3.986004418e14';
Re_s_str = '6378137';
J2_s_str = '1.0826267e-3';

% Compute k_2 = J_2 * R_e^2 / 2 symbolically-exactly and then evaluate.
k2_exact = double(sym(J2_s_str) * sym(Re_s_str)^2 / sym(2));
r_num  = double(sym(r_s_str));
mu_num = double(sym(mu_s_str));
k2_num = k2_exact;
reference = mu_num * k2_num / r_num^3;
printf('Reference scale mu*k2/r^3 at r=7e6 m: %.6e (units m^2/s^2)\n', reference);

% Rather than compute each corner via subs (which re-invokes sym), we use the
% closed form directly in doubles, then separately verify that chain 1 returns
% the same value via exact-sym subs.
%
% (I, f, g) test points:
test_cases = {
    {'0',         '0.3', '0.5'},  % equatorial, arbitrary f,g -> mu*k2/r^3 independent of f,g
    {'0',         '1.2', '2.0'},
    {'0',         '4.7', '0.5'},
    {'pi/2',      '0',   '0'},    % polar, equator crossing -> mu*k2/r^3
    {'pi/2',      '0',   'pi/2'}, % polar, pole passage -> -2*mu*k2/r^3
    {'pi/4',      'pi/6','pi/3'}, % generic
};

for tc = 1:numel(test_cases)
    I_str = test_cases{tc}{1};
    f_str = test_cases{tc}{2};
    g_str = test_cases{tc}{3};
    vars = {r, mu, k2, I, f, g};
    vals_str = {r_s_str, mu_s_str, sprintf('%.17g', k2_num), I_str, f_str, g_str};
    v_c1 = eval_exact(F1_chain1, vars, vals_str);
    v_c2 = eval_exact(F1_chain2, vars, vals_str);

    % Compute expected value in closed form (symbolic, then evaluate).
    I_sym = sym(I_str); fg_sym = sym(f_str) + sym(g_str);
    sin2_beta = sin(I_sym)^2 * sin(fg_sym)^2;
    expected = double( -( sym(mu_s_str) * sym(sprintf('%.17g', k2_num)) / sym(r_s_str)^3 ) * (3 * sin2_beta - 1) );

    printf('I=%-6s f=%-6s g=%-6s: chain1=%+.4e  chain2=%+.4e  expected=%+.4e\n', ...
           I_str, f_str, g_str, v_c1, v_c2, expected);
end
disp('');

disp('=== End of ch05 verification ===');
