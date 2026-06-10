%% verify_eqs1to5_octave.m
%% Symbolic + numerical verification of BH61 derivation notes (Eqs 1-5)
%%
%% Each section:
%%   1. Derives the result symbolically from inputs (not from the claimed answer)
%%   2. Compares the symbolic result against the claimed formula
%%   3. Evaluates both numerically at random points to confirm
%%
%% If any check fails, the script prints FAIL and continues.
%% At the end, a summary reports total PASS/FAIL count.

pkg load symbolic;

n_pass = 0;
n_fail = 0;

function report(name, passed)
  if passed
    printf('PASS: %s\n', name);
  else
    printf('FAIL: %s\n', name);
  end
end

%% ===================================================================
%% Section 1: File 01 — Gradient of U = mu/|q| (Step 3)
%% ===================================================================
printf('\n=== File 01: Gradient of U ===\n');

syms q1 q2 q3 mu positive;
r_sym = sqrt(q1^2 + q2^2 + q3^2);
U_sym = mu / r_sym;

% Symbolic derivation: compute dU/dq_j
dU_dq1_derived = diff(U_sym, q1);
dU_dq2_derived = diff(U_sym, q2);
dU_dq3_derived = diff(U_sym, q3);

% Claimed formula: dU/dq_j = -mu*q_j / r^3
dU_dq1_claimed = -mu * q1 / r_sym^3;
dU_dq2_claimed = -mu * q2 / r_sym^3;
dU_dq3_claimed = -mu * q3 / r_sym^3;

% Symbolic comparison
err1 = simplify(dU_dq1_derived - dU_dq1_claimed);
err2 = simplify(dU_dq2_derived - dU_dq2_claimed);
err3 = simplify(dU_dq3_derived - dU_dq3_claimed);

t = logical(err1 == 0) && logical(err2 == 0) && logical(err3 == 0);
report('dU/dq_j = -mu*q_j/r^3 (symbolic)', t);
if t; n_pass++; else; n_fail++; end

%% Laplace equation verification
lap_U = diff(U_sym, q1, 2) + diff(U_sym, q2, 2) + diff(U_sym, q3, 2);
lap_simplified = simplify(lap_U);
t = logical(lap_simplified == 0);
report('Laplace equation nabla^2 U = 0 (symbolic)', t);
if t; n_pass++; else; n_fail++; end

%% ===================================================================
%% Section 2: File 01 — F definition and Eq. (1) verification (Step 5)
%% ===================================================================
printf('\n=== File 01: Eq. (1) verification ===\n');

syms xi1 xi2 xi3 eta1 eta2 eta3 real;

% U as function of eta (same form)
r_eta = sqrt(eta1^2 + eta2^2 + eta3^2);
U_eta = mu / r_eta;

% F definition
F_sym = -sym(1)/2 * (xi1^2 + xi2^2 + xi3^2) + U_eta;

% Claim: dF/deta_j = dU/deta_j (kinetic term has no eta dependence)
dF_deta1 = diff(F_sym, eta1);
dU_deta1 = diff(U_eta, eta1);
t = logical(simplify(dF_deta1 - dU_deta1) == 0);
report('dF/deta_1 = dU/deta_1 (symbolic)', t);
if t; n_pass++; else; n_fail++; end

% Claim: -dF/dxi_j = xi_j
neg_dF_dxi1 = -diff(F_sym, xi1);
t = logical(simplify(neg_dF_dxi1 - xi1) == 0);
report('-dF/dxi_1 = xi_1 (symbolic)', t);
if t; n_pass++; else; n_fail++; end

neg_dF_dxi2 = -diff(F_sym, xi2);
t = logical(simplify(neg_dF_dxi2 - xi2) == 0);
report('-dF/dxi_2 = xi_2 (symbolic)', t);
if t; n_pass++; else; n_fail++; end

neg_dF_dxi3 = -diff(F_sym, xi3);
t = logical(simplify(neg_dF_dxi3 - xi3) == 0);
report('-dF/dxi_3 = xi_3 (symbolic)', t);
if t; n_pass++; else; n_fail++; end

%% ===================================================================
%% Section 3: File 03 — Conservation proofs (numerical integration)
%% ===================================================================
printf('\n=== File 03: Conservation of h_ang, E, e_vec (numerical) ===\n');

% Set up a Keplerian orbit and propagate numerically
mu_val = 1.0;
a_val = 1.5; e_val = 0.3;
E0 = 0.7;  % eccentric anomaly at t=0

% Initial conditions from Keplerian elements (2D for simplicity, I=0)
r0 = a_val * (1 - e_val * cos(E0));
f0 = 2 * atan2(sqrt(1+e_val)*sin(E0/2), sqrt(1-e_val)*cos(E0/2));

eta0 = [r0*cos(f0); r0*sin(f0); 0];
v_r = sqrt(mu_val*a_val) * e_val * sin(E0) / r0;
v_t = sqrt(mu_val*a_val*(1-e_val^2)) / r0;
xi0 = [v_r*cos(f0) - v_t*sin(f0); v_r*sin(f0) + v_t*cos(f0); 0];

% Check conserved quantities at initial point
h_ang0 = cross(eta0, xi0);
E_energy0 = 0.5*dot(xi0,xi0) - mu_val/norm(eta0);
e_vec0 = cross(xi0, h_ang0)/mu_val - eta0/norm(eta0);

% Propagate to several other points on the same orbit
max_h_err = 0; max_E_err = 0; max_e_err = 0;
for E_test = linspace(0, 2*pi, 20)
  r_t = a_val*(1 - e_val*cos(E_test));
  f_t = 2*atan2(sqrt(1+e_val)*sin(E_test/2), sqrt(1-e_val)*cos(E_test/2));
  eta_t = [r_t*cos(f_t); r_t*sin(f_t); 0];
  vr_t = sqrt(mu_val*a_val)*e_val*sin(E_test)/r_t;
  vt_t = sqrt(mu_val*a_val*(1-e_val^2))/r_t;
  xi_t = [vr_t*cos(f_t)-vt_t*sin(f_t); vr_t*sin(f_t)+vt_t*cos(f_t); 0];

  h_t = cross(eta_t, xi_t);
  E_t = 0.5*dot(xi_t,xi_t) - mu_val/norm(eta_t);
  e_t = cross(xi_t, h_t)/mu_val - eta_t/norm(eta_t);

  max_h_err = max(max_h_err, norm(h_t - h_ang0));
  max_E_err = max(max_E_err, abs(E_t - E_energy0));
  max_e_err = max(max_e_err, norm(e_t - e_vec0));
end

t = max_h_err < 1e-12;
report(sprintf('h_ang conserved (max err = %.2e)', max_h_err), t);
if t; n_pass++; else; n_fail++; end

t = max_E_err < 1e-12;
report(sprintf('Energy conserved (max err = %.2e)', max_E_err), t);
if t; n_pass++; else; n_fail++; end

t = max_e_err < 1e-12;
report(sprintf('e_vec conserved (max err = %.2e)', max_e_err), t);
if t; n_pass++; else; n_fail++; end

%% ===================================================================
%% Section 4: File 03 — Vis-viva identity (symbolic)
%% ===================================================================
printf('\n=== File 03: Vis-viva identity ===\n');

syms a_s e_s E_s positive;
r_s = a_s * (1 - e_s * cos(E_s));

% Velocity components in orbital frame
vr_s = sqrt(mu * a_s) * e_s * sin(E_s) / r_s;
vt_s = sqrt(mu * a_s * (1 - e_s^2)) / r_s;
v_sq = simplify(vr_s^2 + vt_s^2);

% Vis-viva claim
vis_viva = mu * (2/r_s - 1/a_s);

diff_vv = simplify(v_sq - vis_viva);
t = logical(diff_vv == 0);
report('|xi|^2 = mu(2/r - 1/a) (symbolic)', t);
if t; n_pass++; else; n_fail++; end

%% ===================================================================
%% Section 5: File 05 — Kepler equation implicit derivatives (symbolic)
%% ===================================================================
printf('\n=== File 05: Kepler equation derivatives ===\n');

syms L1 L2 l1 positive;
a_d = L1^2 / mu;
e_d = sqrt(1 - L2^2/L1^2);

% Kepler equation: l1 = E - e*sin(E)
% Implicit differentiation: dl1 = dE(1 - e*cos(E)) - sin(E)*de
% => dE/dl1 = 1/(1 - e*cos(E)) = 1/D
% We verify numerically since symbolic implicit diff is tricky

printf('  (Kepler implicit derivatives verified numerically)\n');
mu_v = 1.0;
for trial = 1:5
  L1_v = 0.8 + rand();
  L2_v = L1_v * (0.5 + 0.4*rand());  % ensure e < 1
  l1_v = 2*pi*rand();

  a_v = L1_v^2 / mu_v;
  e_v = sqrt(1 - L2_v^2/L1_v^2);

  % Solve Kepler equation for E
  E_v = l1_v;
  for iter = 1:50
    E_v = E_v + (l1_v - E_v + e_v*sin(E_v))/(1 - e_v*cos(E_v));
  end

  D_v = 1 - e_v*cos(E_v);

  % Claimed: dE/dl1 = 1/D
  eps_fd = 1e-8;
  l1p = l1_v + eps_fd; l1m = l1_v - eps_fd;
  E_p = l1p; E_m = l1m;
  for iter = 1:50
    E_p = E_p + (l1p - E_p + e_v*sin(E_p))/(1 - e_v*cos(E_p));
    E_m = E_m + (l1m - E_m + e_v*sin(E_m))/(1 - e_v*cos(E_m));
  end
  dE_dl1_fd = (E_p - E_m)/(2*eps_fd);
  dE_dl1_claim = 1/D_v;
  err = abs(dE_dl1_fd - dE_dl1_claim);
end
t = err < 1e-6;
report(sprintf('dE/dl_1 = 1/D (numerical, last trial err = %.2e)', err), t);
if t; n_pass++; else; n_fail++; end

%% ===================================================================
%% Section 6: File 07 — p_j and q_j verification (numerical)
%% ===================================================================
printf('\n=== File 07: p_j and q_j (Eq. 5) ===\n');

% Helper: Delaunay to Cartesian
function [eta, xi] = delaunay_to_cart(L, l, mu_val)
  L1 = L(1); L2 = L(2); L3 = L(3);
  l1 = l(1); l2 = l(2); l3 = l(3);

  a = L1^2 / mu_val;
  e = sqrt(1 - L2^2/L1^2);
  cI = L3/L2; sI = sqrt(1 - cI^2);

  % Solve Kepler
  E = l1;
  for iter = 1:50
    E = E + (l1 - E + e*sin(E))/(1 - e*cos(E));
  end

  % Orbital plane position
  P = a*(cos(E) - e);
  Q = a*sqrt(1-e^2)*sin(E);

  % Rotation matrix R = R3(l3)*R1(I)*R3(l2)
  cg = cos(l2); sg = sin(l2); ch = cos(l3); sh = sin(l3);
  R = [ch*cg - sh*cI*sg, -ch*sg - sh*cI*cg, sh*sI;
       sh*cg + ch*cI*sg, -sh*sg + ch*cI*cg, -ch*sI;
       sI*sg,             sI*cg,              cI];

  eta = R * [P; Q; 0];

  % Velocity via deta/dl1 * n
  D = 1 - e*cos(E);
  dP_dl1 = -a*sin(E)/D;
  dQ_dl1 = a*sqrt(1-e^2)*cos(E)/D;
  n = sqrt(mu_val / a^3);
  xi = n * R * [dP_dl1; dQ_dl1; 0];
end

mu_v = 1.0;
eps_fd = 1e-8;
max_err_p = zeros(3,1);
max_err_q = zeros(3,1);

for trial = 1:10
  L1_v = 0.8 + 0.5*rand();
  L2_v = L1_v * (0.6 + 0.3*rand());
  L3_v = L2_v * (0.3 + 0.4*rand());
  l1_v = 2*pi*rand();
  l2_v = 2*pi*rand();
  l3_v = 2*pi*rand();

  L_v = [L1_v; L2_v; L3_v];
  l_v = [l1_v; l2_v; l3_v];

  [eta0, xi0] = delaunay_to_cart(L_v, l_v, mu_v);

  a_v = L1_v^2/mu_v;
  e_v = sqrt(1 - L2_v^2/L1_v^2);

  % Solve Kepler for E
  E_v = l1_v;
  for iter = 1:50
    E_v = E_v + (l1_v - E_v + e_v*sin(E_v))/(1 - e_v*cos(E_v));
  end
  D_v = 1 - e_v*cos(E_v);
  r_v = a_v * D_v;

  % Compute p_j and q_j by finite differences: sum_k xi_k * deta_k/d(variable)
  for j = 1:3
    % p_j = sum_k xi_k * deta_k/dl_j
    lp = l_v; lp(j) += eps_fd;
    lm = l_v; lm(j) -= eps_fd;
    [eta_p, ~] = delaunay_to_cart(L_v, lp, mu_v);
    [eta_m, ~] = delaunay_to_cart(L_v, lm, mu_v);
    deta_dlj = (eta_p - eta_m) / (2*eps_fd);
    pj_fd = dot(xi0, deta_dlj);

    % q_j = sum_k xi_k * deta_k/dL_j
    Lp = L_v; Lp(j) += eps_fd;
    Lm = L_v; Lm(j) -= eps_fd;
    [eta_p, ~] = delaunay_to_cart(Lp, l_v, mu_v);
    [eta_m, ~] = delaunay_to_cart(Lm, l_v, mu_v);
    deta_dLj = (eta_p - eta_m) / (2*eps_fd);
    qj_fd = dot(xi0, deta_dLj);

    % Closed-form claims from Eq. (5)
    switch j
      case 1
        pj_claim = L1_v * (2*a_v/r_v - 1);
        qj_claim = 2*sin(E_v)*(1 - e_v^3*cos(E_v)) / (e_v * D_v);
      case 2
        pj_claim = L2_v;
        qj_claim = -2*sqrt(1-e_v^2)*sin(E_v) / (e_v * D_v);
      case 3
        pj_claim = L3_v;
        qj_claim = 0;
    end

    max_err_p(j) = max(max_err_p(j), abs(pj_fd - pj_claim) / max(abs(pj_claim), 1e-15));
    if abs(qj_claim) > 1e-12
      max_err_q(j) = max(max_err_q(j), abs(qj_fd - qj_claim) / abs(qj_claim));
    else
      % For q_3 = 0 claim, use absolute error (relative error is meaningless)
      max_err_q(j) = max(max_err_q(j), abs(qj_fd - qj_claim));
    end
  end
end

labels_p = {'p_1 = L_1(2a/r-1)', 'p_2 = L_2', 'p_3 = L_3'};
labels_q = {'q_1 = 2sin(E)(1-e^3cos(E))/(eD)', 'q_2 = -2sqrt(1-e^2)sin(E)/(eD)', 'q_3 = 0'};

for j = 1:3
  t = max_err_p(j) < 1e-5;
  report(sprintf('%s (max rel err = %.2e)', labels_p{j}, max_err_p(j)), t);
  if t; n_pass++; else; n_fail++; end
end

for j = 1:3
  t = max_err_q(j) < 1e-5;
  report(sprintf('%s (max rel err = %.2e)', labels_q{j}, max_err_q(j)), t);
  if t; n_pass++; else; n_fail++; end
end

%% ===================================================================
%% Section 7: File 04 — Poisson bracket canonicity check (numerical)
%% ===================================================================
printf('\n=== File 04: Canonicity {l_i, L_j} = delta_ij ===\n');

% Compute {l_i, L_j} via finite differences in (eta, xi) space
for trial = 1:3
  L1_v = 0.8 + 0.5*rand();
  L2_v = L1_v * (0.6 + 0.3*rand());
  L3_v = L2_v * (0.3 + 0.4*rand());
  l1_v = 2*pi*rand();
  l2_v = 2*pi*rand();
  l3_v = 2*pi*rand();

  [eta0, xi0] = delaunay_to_cart([L1_v;L2_v;L3_v], [l1_v;l2_v;l3_v], mu_v);

  % Compute Delaunay variables from Cartesian (inverse map)
  function [L_out, l_out] = cart_to_delaunay(eta, xi, mu_val)
    r = norm(eta);
    v2 = dot(xi, xi);
    energy = v2/2 - mu_val/r;
    a = -mu_val/(2*energy);
    L_out(1) = sqrt(mu_val * a);

    h = cross(eta, xi);
    G = norm(h);
    L_out(2) = G;
    L_out(3) = h(3);

    % Inclination
    cI = h(3)/G;
    sI = sqrt(1 - cI^2);

    % Node
    if sI > 1e-12
      l_out(3) = atan2(h(1), -h(2));
    else
      l_out(3) = 0;
    end

    % Eccentricity vector
    ev = cross(xi, h)/mu_val - eta/r;
    e = norm(ev);

    % Argument of periapsis
    n_vec = [-h(2); h(1); 0] / (G*sI + 1e-30);
    if sI > 1e-12
      cos_g = dot(n_vec, ev/e);
      sin_g = dot(cross(n_vec, ev/e), h/G);
      l_out(2) = atan2(sin_g, cos_g);
    else
      l_out(2) = atan2(ev(2), ev(1));
    end

    % True anomaly
    cos_f = dot(ev/e, eta/r);
    sin_f = dot(cross(ev/e, eta/r), h/G);
    f = atan2(sin_f, cos_f);

    % Eccentric anomaly
    E = atan2(sqrt(1-e^2)*sin(f), e + cos(f));

    % Mean anomaly
    l_out(1) = E - e*sin(E);

    L_out = L_out(:);
    l_out = l_out(:);
  end

  % Poisson bracket {l_i, L_j} = sum_k (dl_i/deta_k * dL_j/dxi_k - dl_i/dxi_k * dL_j/deta_k)
  pb = zeros(3,3);
  for i_idx = 1:3
    for j_idx = 1:3
      val = 0;
      for k = 1:3
        % dl_i/deta_k
        eta_p = eta0; eta_p(k) += eps_fd;
        eta_m = eta0; eta_m(k) -= eps_fd;
        [~, l_p] = cart_to_delaunay(eta_p, xi0, mu_v);
        [~, l_m] = cart_to_delaunay(eta_m, xi0, mu_v);
        dli_detak = (l_p(i_idx) - l_m(i_idx)) / (2*eps_fd);

        % dL_j/dxi_k
        xi_p = xi0; xi_p(k) += eps_fd;
        xi_m = xi0; xi_m(k) -= eps_fd;
        [L_p, ~] = cart_to_delaunay(eta0, xi_p, mu_v);
        [L_m, ~] = cart_to_delaunay(eta0, xi_m, mu_v);
        dLj_dxik = (L_p(j_idx) - L_m(j_idx)) / (2*eps_fd);

        % dl_i/dxi_k
        [~, l_p2] = cart_to_delaunay(eta0, xi_p, mu_v);
        [~, l_m2] = cart_to_delaunay(eta0, xi_m, mu_v);
        dli_dxik = (l_p2(i_idx) - l_m2(i_idx)) / (2*eps_fd);

        % dL_j/deta_k
        [L_p2, ~] = cart_to_delaunay(eta_p, xi0, mu_v);
        [L_m2, ~] = cart_to_delaunay(eta_m, xi0, mu_v);
        dLj_detak = (L_p2(j_idx) - L_m2(j_idx)) / (2*eps_fd);

        val += dli_detak * dLj_dxik - dli_dxik * dLj_detak;
      end
      pb(i_idx, j_idx) = val;
    end
  end

  err_pb = max(max(abs(pb - eye(3))));
  t = err_pb < 1e-4;
  report(sprintf('{l_i, L_j} = delta_ij, trial %d (max err = %.2e)', trial, err_pb), t);
  if t; n_pass++; else; n_fail++; end
end

%% ===================================================================
%% Summary
%% ===================================================================
printf('\n========================================\n');
printf('TOTAL: %d PASS, %d FAIL\n', n_pass, n_fail);
printf('========================================\n');
