%% verify_eq4_derivation.m
%% Derive and verify BH61 Eq. (4) from scratch using numerical computation.
%%
%% The chain:
%%   1. L_j and l_j are functions of (xi, eta)
%%   2. dL_j/dt = sum_k (dL_j/dxi_k)(dxi_k/dt) + sum_k (dL_j/deta_k)(deta_k/dt)
%%   3. Substitute Eq. (1): dxi_k/dt = dF/deta_k + X_k, deta_k/dt = -dF/dxi_k
%%   4. The conservative part gives dF/dl_j (by canonicity)
%%   5. The non-conservative part is P_j = sum_k (dL_j/dxi_k) X_k
%%   6. Use symplectic identity: dL_j/dxi_k = deta_k/dl_j
%%   7. Therefore P_j = sum_k X_k deta_k/dl_j
%%
%% We verify steps 5, 6, 7 numerically at random orbital elements.

pkg load symbolic;
warning('off', 'all');

OCTAVE = "C:/Program Files/GNU Octave/Octave-11.1.0/mingw64/bin/octave-cli.exe";
n_pass = 0;
n_fail = 0;

function report(name, passed, err)
  if passed
    printf('PASS: %s (err = %.2e)\n', name, err);
  else
    printf('FAIL: %s (err = %.2e)\n', name, err);
  end
end

%% Helper: Delaunay to Cartesian (forward map)
function [eta, xi] = del2cart(L, l, mu)
  L1 = L(1); L2 = L(2); L3 = L(3);
  l1 = l(1); l2 = l(2); l3 = l(3);
  a = L1^2/mu;
  e = sqrt(1 - L2^2/L1^2);
  cI = L3/L2; sI = sqrt(1 - cI^2);

  % Solve Kepler
  E = l1;
  for iter = 1:50
    E = E + (l1 - E + e*sin(E))/(1 - e*cos(E));
  end
  D = 1 - e*cos(E);

  % Orbital plane
  P = a*(cos(E) - e);
  Q = a*sqrt(1-e^2)*sin(E);

  % Rotation R = R3(l3)*R1(I)*R3(l2)
  cg = cos(l2); sg = sin(l2); ch = cos(l3); sh = sin(l3);
  R = [ch*cg - sh*cI*sg, -ch*sg - sh*cI*cg, sh*sI;
       sh*cg + ch*cI*sg, -sh*sg + ch*cI*cg, -ch*sI;
       sI*sg,             sI*cg,              cI];

  eta = R * [P; Q; 0];

  % Velocity
  dP = -a*sin(E)/D;
  dQ = a*sqrt(1-e^2)*cos(E)/D;
  n = sqrt(mu/a^3);
  xi = n * R * [dP; dQ; 0];
end

%% Helper: Cartesian to Delaunay (inverse map)
function [L, l] = cart2del(eta, xi, mu)
  r = norm(eta);
  v2 = dot(xi, xi);
  energy = v2/2 - mu/r;
  a = -mu/(2*energy);
  L(1) = sqrt(mu*a);

  h = cross(eta, xi);
  G = norm(h);
  L(2) = G;
  L(3) = h(3);

  cI = h(3)/G;
  sI = sqrt(1 - cI^2);

  % Node
  if sI > 1e-12
    l(3) = atan2(h(1), -h(2));
  else
    l(3) = 0;
  end

  % Eccentricity vector
  ev = cross(xi, h)/mu - eta/r;
  e = norm(ev);

  % Argument of periapsis
  if sI > 1e-12
    n_vec = [-h(2); h(1); 0]/(G*sI);
    cos_g = dot(n_vec, ev/e);
    sin_g = dot(cross(n_vec, ev/e), h/G);
    l(2) = atan2(sin_g, cos_g);
  else
    l(2) = atan2(ev(2), ev(1));
  end

  % True anomaly
  cos_f = dot(ev/e, eta/r);
  sin_f = dot(cross(ev/e, eta/r), h/G);
  f = atan2(sin_f, cos_f);

  % Eccentric anomaly
  E_anom = atan2(sqrt(1-e^2)*sin(f), e + cos(f));
  l(1) = E_anom - e*sin(E_anom);

  L = L(:); l = l(:);
end

mu_v = 1.0;
eps_fd = 1e-7;

printf('\n============================================================\n');
printf('VERIFICATION OF BH61 Eq. (4) DERIVATION\n');
printf('============================================================\n');

%% ================================================================
%% TEST 1: Verify the chain rule for dL_j/dt
%%
%% dL_j/dt = sum_k (dL_j/dxi_k)(dxi_k/dt) + sum_k (dL_j/deta_k)(deta_k/dt)
%%
%% With Eq. (1): dxi_k/dt = dF/deta_k + X_k, deta_k/dt = -dF/dxi_k = xi_k
%%
%% Conservative part: sum_k (dL_j/dxi_k)(dF/deta_k) + sum_k (dL_j/deta_k)(xi_k)
%%   = {L_j, F} (the Poisson bracket, which equals dF/dl_j by canonicity)
%%
%% Non-conservative part: sum_k (dL_j/dxi_k) X_k
%%
%% So: P_j = sum_k (dL_j/dxi_k) X_k
%%
%% Now the question: does dL_j/dxi_k = deta_k/dl_j (symplectic identity SI-1)?
%% And does this give P_j = sum_k X_k deta_k/dl_j (BH61 Eq. 4a)?
%% ================================================================

printf('\n--- TEST 1: Symplectic identity (SI-1) ---\n');
printf('    dL_j/dxi_k = deta_k/dl_j\n\n');

max_err_si1 = 0;
for trial = 1:5
  L1 = 0.8 + 0.5*rand(); L2 = L1*(0.6+0.3*rand()); L3 = L2*(0.3+0.4*rand());
  l1 = 2*pi*rand(); l2 = 2*pi*rand(); l3 = 2*pi*rand();
  Lv = [L1;L2;L3]; lv = [l1;l2;l3];

  [eta0, xi0] = del2cart(Lv, lv, mu_v);

  for j = 1:3
    for k = 1:3
      % Compute dL_j/dxi_k by finite differences on the inverse map
      xi_p = xi0; xi_p(k) += eps_fd;
      xi_m = xi0; xi_m(k) -= eps_fd;
      [Lp, ~] = cart2del(eta0, xi_p, mu_v);
      [Lm, ~] = cart2del(eta0, xi_m, mu_v);
      dLj_dxik = (Lp(j) - Lm(j)) / (2*eps_fd);

      % Compute deta_k/dl_j by finite differences on the forward map
      lp = lv; lp(j) += eps_fd;
      lm = lv; lm(j) -= eps_fd;
      [etap, ~] = del2cart(Lv, lp, mu_v);
      [etam, ~] = del2cart(Lv, lm, mu_v);
      detak_dlj = (etap(k) - etam(k)) / (2*eps_fd);

      err = abs(dLj_dxik - detak_dlj);
      max_err_si1 = max(max_err_si1, err);
    end
  end
end

t = max_err_si1 < 1e-4;
report('SI-1: dL_j/dxi_k = deta_k/dl_j (9 components, 5 trials)', t, max_err_si1);
if t; n_pass++; else; n_fail++; end

%% ================================================================
%% TEST 2: Symplectic identity (SI-2)
%%    dl_j/dxi_k = -deta_k/dL_j
%% ================================================================

printf('\n--- TEST 2: Symplectic identity (SI-2) ---\n');
printf('    dl_j/dxi_k = -deta_k/dL_j\n\n');

max_err_si2 = 0;
for trial = 1:5
  L1 = 0.8 + 0.5*rand(); L2 = L1*(0.6+0.3*rand()); L3 = L2*(0.3+0.4*rand());
  l1 = 2*pi*rand(); l2 = 2*pi*rand(); l3 = 2*pi*rand();
  Lv = [L1;L2;L3]; lv = [l1;l2;l3];

  [eta0, xi0] = del2cart(Lv, lv, mu_v);

  for j = 1:3
    for k = 1:3
      % dl_j/dxi_k
      xi_p = xi0; xi_p(k) += eps_fd;
      xi_m = xi0; xi_m(k) -= eps_fd;
      [~, lp] = cart2del(eta0, xi_p, mu_v);
      [~, lm] = cart2del(eta0, xi_m, mu_v);
      dlj_dxik = (lp(j) - lm(j)) / (2*eps_fd);

      % -deta_k/dL_j
      Lp = Lv; Lp(j) += eps_fd;
      Lm = Lv; Lm(j) -= eps_fd;
      [etap, ~] = del2cart(Lp, lv, mu_v);
      [etam, ~] = del2cart(Lm, lv, mu_v);
      neg_detak_dLj = -(etap(k) - etam(k)) / (2*eps_fd);

      err = abs(dlj_dxik - neg_detak_dLj);
      max_err_si2 = max(max_err_si2, err);
    end
  end
end

t = max_err_si2 < 1e-4;
report('SI-2: dl_j/dxi_k = -deta_k/dL_j (9 components, 5 trials)', t, max_err_si2);
if t; n_pass++; else; n_fail++; end

%% ================================================================
%% TEST 3: P_j from the chain rule vs P_j from Eq. (4a)
%%
%% Chain rule: P_j = sum_k (dL_j/dxi_k) X_k
%% Eq. (4a):  P_j = sum_k X_k (deta_k/dl_j)
%%
%% These should be equal (by SI-1).
%% Also verify both equal -A*V*exp(-alpha*r)*p_j where p_j = sum_k xi_k deta_k/dl_j
%% ================================================================

printf('\n--- TEST 3: P_j chain rule vs Eq. (4a) vs factored form ---\n\n');

A_drag = 0.1;  % drag coefficient
alpha_v = 0.5; % inverse scale height

max_err_P = 0;
for trial = 1:5
  L1 = 0.8 + 0.5*rand(); L2 = L1*(0.6+0.3*rand()); L3 = L2*(0.3+0.4*rand());
  l1 = 2*pi*rand(); l2 = 2*pi*rand(); l3 = 2*pi*rand();
  Lv = [L1;L2;L3]; lv = [l1;l2;l3];

  [eta0, xi0] = del2cart(Lv, lv, mu_v);
  r0 = norm(eta0);
  V0 = norm(xi0);

  % Drag acceleration: X_k = -A*V*exp(-alpha*r)*xi_k
  Xdrag = -A_drag * V0 * exp(-alpha_v * r0) * xi0;

  % Method A: chain rule P_j = sum_k (dL_j/dxi_k) X_k
  Pj_chain = zeros(3,1);
  for j = 1:3
    for k = 1:3
      xi_p = xi0; xi_p(k) += eps_fd;
      xi_m = xi0; xi_m(k) -= eps_fd;
      [Lp, ~] = cart2del(eta0, xi_p, mu_v);
      [Lm, ~] = cart2del(eta0, xi_m, mu_v);
      dLj_dxik = (Lp(j) - Lm(j)) / (2*eps_fd);
      Pj_chain(j) += dLj_dxik * Xdrag(k);
    end
  end

  % Method B: Eq. (4a) P_j = sum_k X_k (deta_k/dl_j)
  Pj_eq4a = zeros(3,1);
  for j = 1:3
    lp = lv; lp(j) += eps_fd;
    lm = lv; lm(j) -= eps_fd;
    [etap, ~] = del2cart(Lv, lp, mu_v);
    [etam, ~] = del2cart(Lv, lm, mu_v);
    deta_dlj = (etap - etam) / (2*eps_fd);
    Pj_eq4a(j) = dot(Xdrag, deta_dlj);
  end

  % Method C: factored form P_j = -A*V*exp(-alpha*r) * p_j
  pj = zeros(3,1);
  for j = 1:3
    lp = lv; lp(j) += eps_fd;
    lm = lv; lm(j) -= eps_fd;
    [etap, ~] = del2cart(Lv, lp, mu_v);
    [etam, ~] = del2cart(Lv, lm, mu_v);
    deta_dlj = (etap - etam) / (2*eps_fd);
    pj(j) = dot(xi0, deta_dlj);
  end
  Pj_factored = -A_drag * V0 * exp(-alpha_v * r0) * pj;

  err_AB = max(abs(Pj_chain - Pj_eq4a));
  err_AC = max(abs(Pj_chain - Pj_factored));
  max_err_P = max(max_err_P, max(err_AB, err_AC));
end

t = max_err_P < 1e-4;
report('P_j: chain rule = Eq.(4a) = factored (3 methods, 5 trials)', t, max_err_P);
if t; n_pass++; else; n_fail++; end

%% ================================================================
%% TEST 4: Q_j from the chain rule vs Q_j from Eq. (4b)
%%
%% The non-conservative part of dl_j/dt is: sum_k (dl_j/dxi_k) X_k
%% By SI-2: dl_j/dxi_k = -deta_k/dL_j
%% So: sum_k (dl_j/dxi_k) X_k = -sum_k X_k deta_k/dL_j = -Q_j
%%
%% BH61 defines Q_j = sum_k X_k deta_k/dL_j (POSITIVE, no minus)
%% and puts -Q_j in the equation of motion.
%%
%% The factored form: Q_j = -A*V*exp(-alpha*r)*q_j
%% where q_j = sum_k xi_k deta_k/dL_j
%% ================================================================

printf('\n--- TEST 4: Q_j chain rule vs Eq. (4b) vs factored form ---\n');
printf('    Also verifying the SIGN: non-cons part of dl_j/dt = -Q_j\n\n');

max_err_Q = 0;
max_err_sign = 0;
for trial = 1:5
  L1 = 0.8 + 0.5*rand(); L2 = L1*(0.6+0.3*rand()); L3 = L2*(0.3+0.4*rand());
  l1 = 2*pi*rand(); l2 = 2*pi*rand(); l3 = 2*pi*rand();
  Lv = [L1;L2;L3]; lv = [l1;l2;l3];

  [eta0, xi0] = del2cart(Lv, lv, mu_v);
  r0 = norm(eta0);
  V0 = norm(xi0);

  Xdrag = -A_drag * V0 * exp(-alpha_v * r0) * xi0;

  % Method A: chain rule — the non-conservative part of dl_j/dt
  % = sum_k (dl_j/dxi_k) X_k
  noncons_dlj = zeros(3,1);
  for j = 1:3
    for k = 1:3
      xi_p = xi0; xi_p(k) += eps_fd;
      xi_m = xi0; xi_m(k) -= eps_fd;
      [~, lp] = cart2del(eta0, xi_p, mu_v);
      [~, lm] = cart2del(eta0, xi_m, mu_v);
      dlj_dxik = (lp(j) - lm(j)) / (2*eps_fd);
      noncons_dlj(j) += dlj_dxik * Xdrag(k);
    end
  end

  % Method B: Q_j = sum_k X_k deta_k/dL_j (BH61 definition, POSITIVE)
  Qj_eq4b = zeros(3,1);
  for j = 1:3
    Lp = Lv; Lp(j) += eps_fd;
    Lm = Lv; Lm(j) -= eps_fd;
    [etap, ~] = del2cart(Lp, lv, mu_v);
    [etam, ~] = del2cart(Lm, lv, mu_v);
    deta_dLj = (etap - etam) / (2*eps_fd);
    Qj_eq4b(j) = dot(Xdrag, deta_dLj);
  end

  % Verify: non-conservative part of dl_j/dt = -Q_j
  err_sign = max(abs(noncons_dlj - (-Qj_eq4b)));
  max_err_sign = max(max_err_sign, err_sign);

  % Method C: factored Q_j = -A*V*exp(-alpha*r) * q_j
  qj = zeros(3,1);
  for j = 1:3
    Lp = Lv; Lp(j) += eps_fd;
    Lm = Lv; Lm(j) -= eps_fd;
    [etap, ~] = del2cart(Lp, lv, mu_v);
    [etam, ~] = del2cart(Lm, lv, mu_v);
    deta_dLj = (etap - etam) / (2*eps_fd);
    qj(j) = dot(xi0, deta_dLj);
  end
  Qj_factored = -A_drag * V0 * exp(-alpha_v * r0) * qj;

  err_BQ = max(abs(Qj_eq4b - Qj_factored));
  max_err_Q = max(max_err_Q, max(err_sign, err_BQ));
end

t = max_err_sign < 1e-4;
report('Non-cons dl_j/dt = -Q_j (sign verified, 5 trials)', t, max_err_sign);
if t; n_pass++; else; n_fail++; end

t = max_err_Q < 1e-4;
report('Q_j: Eq.(4b) = factored form (5 trials)', t, max_err_Q);
if t; n_pass++; else; n_fail++; end

%% ================================================================
%% TEST 5: Verify all six p_j, q_j closed forms against the definition
%% ================================================================

printf('\n--- TEST 5: Eq. (5) closed forms ---\n\n');

max_err_p = zeros(3,1);
max_err_q = zeros(3,1);

for trial = 1:10
  L1 = 0.8 + 0.5*rand(); L2 = L1*(0.6+0.3*rand()); L3 = L2*(0.3+0.4*rand());
  l1 = 2*pi*rand(); l2 = 2*pi*rand(); l3 = 2*pi*rand();
  Lv = [L1;L2;L3]; lv = [l1;l2;l3];

  [eta0, xi0] = del2cart(Lv, lv, mu_v);
  a_v = L1^2/mu_v;
  e_v = sqrt(1 - L2^2/L1^2);
  E_v = l1;
  for iter = 1:50; E_v = E_v + (l1 - E_v + e_v*sin(E_v))/(1 - e_v*cos(E_v)); end
  D_v = 1 - e_v*cos(E_v);
  r_v = a_v*D_v;

  for j = 1:3
    % p_j by finite differences
    lp = lv; lp(j) += eps_fd; lm = lv; lm(j) -= eps_fd;
    [etap,~] = del2cart(Lv, lp, mu_v); [etam,~] = del2cart(Lv, lm, mu_v);
    pj_fd = dot(xi0, (etap-etam)/(2*eps_fd));

    % q_j by finite differences
    Lp = Lv; Lp(j) += eps_fd; Lm = Lv; Lm(j) -= eps_fd;
    [etap,~] = del2cart(Lp, lv, mu_v); [etam,~] = del2cart(Lm, lv, mu_v);
    qj_fd = dot(xi0, (etap-etam)/(2*eps_fd));

    % Closed forms from Eq. (5)
    switch j
      case 1
        pj_claim = L1*(2*a_v/r_v - 1);
        qj_claim = 2*sin(E_v)*(1 - e_v^3*cos(E_v))/(e_v*D_v);
      case 2
        pj_claim = L2;
        qj_claim = -2*sqrt(1-e_v^2)*sin(E_v)/(e_v*D_v);
      case 3
        pj_claim = L3;
        qj_claim = 0;
    end

    max_err_p(j) = max(max_err_p(j), abs(pj_fd - pj_claim)/max(abs(pj_claim),1e-15));
    if abs(qj_claim) > 1e-12
      max_err_q(j) = max(max_err_q(j), abs(qj_fd - qj_claim)/abs(qj_claim));
    else
      max_err_q(j) = max(max_err_q(j), abs(qj_fd - qj_claim));
    end
  end
end

labels_p = {'p_1 = L_1(2a/r-1)', 'p_2 = L_2', 'p_3 = L_3'};
labels_q = {'q_1 = 2sinE(1-e^3cosE)/(eD)', 'q_2 = -2sqrt(1-e^2)sinE/(eD)', 'q_3 = 0'};

for j = 1:3
  t = max_err_p(j) < 1e-5;
  report(labels_p{j}, t, max_err_p(j));
  if t; n_pass++; else; n_fail++; end
end
for j = 1:3
  t = max_err_q(j) < 1e-5;
  report(labels_q{j}, t, max_err_q(j));
  if t; n_pass++; else; n_fail++; end
end

%% ================================================================
%% TEST 6: Verify the COMPLETE equation of motion
%%
%% dL_j/dt = dF/dl_j + P_j  should hold
%% dl_j/dt = -dF/dL_j - Q_j  should hold (note the MINUS on Q_j)
%%
%% Compute the full time derivative numerically and check both sides.
%% ================================================================

printf('\n--- TEST 6: Complete equation of motion (Eq. 3) ---\n');
printf('    dL_j/dt = dF/dl_j + P_j\n');
printf('    dl_j/dt = -dF/dL_j - Q_j  (note: MINUS Q_j)\n\n');

% Use Keplerian F = mu^2/(2*L1^2) for simplicity (no oblateness)
max_err_eom = 0;

for trial = 1:5
  L1 = 0.8 + 0.5*rand(); L2 = L1*(0.6+0.3*rand()); L3 = L2*(0.3+0.4*rand());
  l1 = 2*pi*rand(); l2 = 2*pi*rand(); l3 = 2*pi*rand();
  Lv = [L1;L2;L3]; lv = [l1;l2;l3];

  [eta0, xi0] = del2cart(Lv, lv, mu_v);
  r0 = norm(eta0);
  V0 = norm(xi0);

  Xdrag = -A_drag * V0 * exp(-alpha_v * r0) * xi0;

  % Compute deta/dt and dxi/dt from Eq. (1)
  deta_dt = xi0;  % = -dF/dxi (kinematic relation)
  dxi_dt = -mu_v * eta0 / r0^3 + Xdrag;  % = dF/deta + X (gravity + drag)

  % Compute dL_j/dt and dl_j/dt by chain rule through (xi, eta)
  dLdt_num = zeros(3,1);
  dldt_num = zeros(3,1);
  for j = 1:3
    for k = 1:3
      % dL_j/dxi_k and dL_j/deta_k
      xi_p = xi0; xi_p(k) += eps_fd; xi_m = xi0; xi_m(k) -= eps_fd;
      [Lp,~] = cart2del(eta0, xi_p, mu_v); [Lm,~] = cart2del(eta0, xi_m, mu_v);
      dLj_dxik = (Lp(j)-Lm(j))/(2*eps_fd);

      eta_p = eta0; eta_p(k) += eps_fd; eta_m = eta0; eta_m(k) -= eps_fd;
      [Lp,~] = cart2del(eta_p, xi0, mu_v); [Lm,~] = cart2del(eta_m, xi0, mu_v);
      dLj_detak = (Lp(j)-Lm(j))/(2*eps_fd);

      dLdt_num(j) += dLj_dxik * dxi_dt(k) + dLj_detak * deta_dt(k);

      % dl_j/dxi_k and dl_j/deta_k
      [~,lp2] = cart2del(eta0, xi_p, mu_v); [~,lm2] = cart2del(eta0, xi_m, mu_v);
      dlj_dxik = (lp2(j)-lm2(j))/(2*eps_fd);

      [~,lp2] = cart2del(eta_p, xi0, mu_v); [~,lm2] = cart2del(eta_m, xi0, mu_v);
      dlj_detak = (lp2(j)-lm2(j))/(2*eps_fd);

      dldt_num(j) += dlj_dxik * dxi_dt(k) + dlj_detak * deta_dt(k);
    end
  end

  % Compute dF/dl_j and dF/dL_j for Keplerian F = mu^2/(2*L1^2)
  dF_dlj = [0; 0; 0];  % F doesn't depend on l_j for Kepler
  dF_dLj = [-mu_v^2/L1^3; 0; 0];  % only L1 dependence

  % P_j and Q_j
  Pj = zeros(3,1); Qj = zeros(3,1);
  for j = 1:3
    lp = lv; lp(j) += eps_fd; lm = lv; lm(j) -= eps_fd;
    [etap,~] = del2cart(Lv,lp,mu_v); [etam,~] = del2cart(Lv,lm,mu_v);
    Pj(j) = dot(Xdrag, (etap-etam)/(2*eps_fd));

    Lp = Lv; Lp(j) += eps_fd; Lm = Lv; Lm(j) -= eps_fd;
    [etap,~] = del2cart(Lp,lv,mu_v); [etam,~] = del2cart(Lm,lv,mu_v);
    Qj(j) = dot(Xdrag, (etap-etam)/(2*eps_fd));
  end

  % Check: dL_j/dt = dF/dl_j + P_j
  err_L = max(abs(dLdt_num - (dF_dlj + Pj)));

  % Check: dl_j/dt = -dF/dL_j - Q_j  (MINUS Q_j)
  err_l = max(abs(dldt_num - (-dF_dLj - Qj)));

  max_err_eom = max(max_err_eom, max(err_L, err_l));
end

t = max_err_eom < 1e-3;
report('Full EoM: dL/dt = dF/dl + P, dl/dt = -dF/dL - Q (5 trials)', t, max_err_eom);
if t; n_pass++; else; n_fail++; end

%% ================================================================
%% SUMMARY
%% ================================================================
printf('\n============================================================\n');
printf('TOTAL: %d PASS, %d FAIL\n', n_pass, n_fail);
printf('============================================================\n');
