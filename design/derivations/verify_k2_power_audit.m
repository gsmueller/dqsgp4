% verify_k2_power_audit.m
%
% Independent symbolic audit of k_2-exponents.
% Similar method to mu-audit: extract k_2-exponent via log derivative.

pkg load symbolic;

printf('=========================================================\n');
printf('k_2-exponent independent symbolic audit\n');
printf('=========================================================\n\n');

n_pass = 0;
n_fail = 0;

function p = k2_power(expr, k2_sym)
  expr_simplified = simplify(expr);
  if isequal(expr_simplified, sym(0))
    p = 0;
    return;
  end
  log_derivative = simplify(k2_sym * diff(expr_simplified, k2_sym) / expr_simplified);
  p_num = double(subs(log_derivative, k2_sym, sym(1)));
  p = round(p_num);
endfunction

function print_result(label, expected_k2, actual_k2, n_pass_in, n_fail_in)
  pass_flag = (expected_k2 == actual_k2);
  printf('  %s\n', label);
  printf('    Expected k_2-exponent: %d\n', expected_k2);
  printf('    Derived k_2-exponent:  %d\n', actual_k2);
  if pass_flag
    printf('    Status: PASS\n');
    n_pass_in = n_pass_in + 1;
  else
    printf('    Status: FAIL\n');
    n_fail_in = n_fail_in + 1;
  end
  printf('\n');
  assignin('caller', 'n_pass', n_pass_in);
  assignin('caller', 'n_fail', n_fail_in);
endfunction

syms L G H mu k2 l_s g_s f_s theta_s e_s positive;
eta_s = sqrt(1 - e_s^2);
A_s = (3*theta_s^2 - 1)/2;
B_s = 3*(1 - theta_s^2)/2;

% Definitions
F1 = mu^4 * k2 * (1 + e_s*cos(f_s))^3 * (A_s + B_s*cos(2*(f_s + g_s))) / G^6;
F1star = mu^4 * k2 * A_s * eta_s^3 / G^6;
n_mot = mu^2/L^3;
dS1_dl = (F1 - F1star) / n_mot;
V = (F1 - F1star)^2;

% Ch 10c §5:
P_e = 8 + 24*e_s^2 + 3*e_s^4;
X60 = P_e / (8 * eta_s^9);
X62 = e_s^2 * (6 + e_s^2) / (4 * eta_s^9);
X64 = e_s^4 / (16 * eta_s^9);
T_expr = (A_s^2 + B_s^2/2) * X60 + 2*A_s*B_s * cos(2*g_s) * X62 ...
       + (B_s^2/2) * cos(4*g_s) * X64 - A_s^2 / eta_s^6;

% ch10d c_2^{(T)}, c_4^{(T)}
c2T = 3 * mu^6 * k2^2 * A_s * B_s * e_s^2 * (6 + e_s^2) * eta_s / (4 * G^10);
c4T = 3 * mu^6 * k2^2 * B_s^2 * e_s^4 * eta_s / (64 * G^10);

% ch10e c_0^{(T)}
Q = A_s^2 + B_s^2/2;
R_sym = A_s^2;
c0T = 3 * mu^6 * k2^2 * eta_s / (2 * G^10) * (Q * P_e/8 - R_sym * eta_s^3);

% ch11a S_1^{*,(T)}
S1star = mu^2 * k2 / ((5*theta_s^2 - 1) * eta_s^2 * G^3) * ...
         (A_s * B_s * e_s^2 * (6 + e_s^2)/4 * sin(2*g_s) + ...
          B_s^2 * e_s^4 / 128 * sin(4*g_s));

% ch06e B.5.3-5
B53 = -3 * mu^4 * k2 * A_s * eta_s^4 / G^7;
B54 = -3 * mu^4 * k2 * eta_s^3 * (5*theta_s^2 - 1) / (2 * G^7);
B55 = 3 * mu^4 * k2 * theta_s * eta_s^3 / G^7;

% ch11d dF_2^{**,(T)}/dH
dF2_dH = (9 * mu^6 * k2^2 * theta_s / (2 * G^11)) * ...
         ((9*theta_s^2 - 5) * P_e * eta_s / 16 - (3*theta_s^2 - 1) * eta_s^4);

printf('=== k_2-exponent checks ===\n\n');
print_result('F_1',           1, k2_power(F1, k2), n_pass, n_fail);
print_result('F_1^*',         1, k2_power(F1star, k2), n_pass, n_fail);
print_result('dS_1/dl',       1, k2_power(dS1_dl, k2), n_pass, n_fail);
print_result('V = (F_1-F_1^*)^2',  2, k2_power(V, k2), n_pass, n_fail);
print_result('B.5.3 dF_1^*/dL',  1, k2_power(B53, k2), n_pass, n_fail);
print_result('B.5.4 dF_1^*/dG',  1, k2_power(B54, k2), n_pass, n_fail);
print_result('B.5.5 dF_1^*/dH',  1, k2_power(B55, k2), n_pass, n_fail);
print_result('c_2^{(T)}',     2, k2_power(c2T, k2), n_pass, n_fail);
print_result('c_4^{(T)}',     2, k2_power(c4T, k2), n_pass, n_fail);
print_result('c_0^{(T)}',     2, k2_power(c0T, k2), n_pass, n_fail);
print_result('S_1^{*,(T)}',   1, k2_power(S1star, k2), n_pass, n_fail);
print_result('dF_2^{**,(T)}/dH', 2, k2_power(dF2_dH, k2), n_pass, n_fail);

printf('=========================================================\n');
if n_fail == 0
  printf('k_2-exponent audit: ALL %d checks PASSED.\n', n_pass);
else
  printf('FAILED: %d pass, %d fail\n', n_pass, n_fail);
end
printf('=========================================================\n');
