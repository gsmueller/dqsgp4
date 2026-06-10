%% verify_K6_osculating.m
%%
%% Numerical verification of Theorem K6 (proof_K6_osculating_elements.md)
%% for the osculating-element map Psi: (r, xi) -> (a, e, I, h, g, f) on D_reg.
%%
%% Checks performed on 500 random elliptic states:
%%  (1) Instantaneous-orbit identity: |xi|^2 = mu*(2/r - 1/a).
%%  (2) Angular-momentum-eccentricity relation: h_ang^2 = mu*a*(1-e^2).
%%  (3) e^2 identity (B.4.3): |e_vec|^2 = 1 + 2*E*h_ang^2/mu^2.
%%  (4) Round-trip: elements -> state (Layer-3 reconstruction) -> elements', matches.
%%
%% This verifier is a numerical sanity check for the algebraic identities of the
%% proof; it does NOT substitute for the proof, which is pure algebra.

function verify_K6_osculating()

  % Earth mu in SI (m^3/s^2); the identity is scale-invariant, but using a
  % physical value keeps dimensional-check sanity obvious.
  mu = 3.986004418e14;

  rng_seed = 20260419;    % reproducible random states
  randn('state', rng_seed);
  rand('state', rng_seed);

  N = 500;                % number of random states
  tol_identity = 1e-10;   % relative tolerance for identities
  tol_rt_pos   = 1e-8;    % round-trip position tolerance (m, then normalized)
  tol_rt_vel   = 1e-8;    % round-trip velocity tolerance (m/s, then normalized)

  n_fail_1 = 0;
  n_fail_2 = 0;
  n_fail_3 = 0;
  n_fail_4 = 0;

  max_err_1 = 0;
  max_err_2 = 0;
  max_err_3 = 0;
  max_err_4 = 0;

  printf('verify_K6_osculating: testing %d random elliptic states\n', N);
  printf('mu = %.6e m^3/s^2\n\n', mu);

  k = 0;
  while k < N
    % sample osculating elements on D_reg
    a = 7e6 + 1e6 * rand();           % 7000 km .. 8000 km
    e = 0.01 + 0.7 * rand();          % [0.01, 0.71)
    I_inc = 0.1 + (pi - 0.2) * rand();% [0.1 rad, pi-0.1 rad]
    h_ang_angle = 2*pi * rand();
    g_ang       = 2*pi * rand();
    f_ang       = 2*pi * rand();

    % Layer 3 reconstruction: (a,e,I,h,g,f) -> (r, xi)
    [r_vec, xi_vec] = elements_to_state(mu, a, e, I_inc, h_ang_angle, g_ang, f_ang);

    % Compute the state's energy; skip if accidentally non-elliptic (shouldn't happen)
    r_mag = norm(r_vec);
    xi_mag2 = dot(xi_vec, xi_vec);
    E_spec = 0.5*xi_mag2 - mu/r_mag;
    if E_spec >= 0
      continue;
    end
    k = k + 1;

    % ---- Forward: state -> elements via Psi ----
    [a2, e2, I2, h2, g2, f2, h_ang_vec, e_vec, h_ang_mag] = state_to_elements(mu, r_vec, xi_vec);

    % Check 1: |xi|^2 = mu*(2/r - 1/a)
    lhs1 = xi_mag2;
    rhs1 = mu*(2/r_mag - 1/a2);
    err1 = abs(lhs1 - rhs1) / abs(rhs1);
    max_err_1 = max(max_err_1, err1);
    if err1 > tol_identity, n_fail_1 = n_fail_1 + 1; end

    % Check 2: h_ang^2 = mu*a*(1-e^2)
    lhs2 = h_ang_mag^2;
    rhs2 = mu * a2 * (1 - e2^2);
    err2 = abs(lhs2 - rhs2) / abs(rhs2);
    max_err_2 = max(max_err_2, err2);
    if err2 > tol_identity, n_fail_2 = n_fail_2 + 1; end

    % Check 3: |e_vec|^2 = 1 + 2*E*h_ang^2/mu^2
    lhs3 = dot(e_vec, e_vec);
    rhs3 = 1 + 2 * E_spec * h_ang_mag^2 / mu^2;
    err3 = abs(lhs3 - rhs3) / max(1e-12, abs(rhs3));
    max_err_3 = max(max_err_3, err3);
    if err3 > tol_identity, n_fail_3 = n_fail_3 + 1; end

    % Check 4: Round-trip consistency.  Compare (a, e, I, h, g, f) on entry
    % with (a2, e2, I2, h2, g2, f2) out.  Angles modulo 2*pi.
    angle_wrap = @(x) mod(x + pi, 2*pi) - pi;
    err_a = abs(a2 - a) / a;
    err_e = abs(e2 - e);
    err_I = abs(angle_wrap(I2 - I_inc));
    err_h = abs(angle_wrap(h2 - h_ang_angle));
    err_g = abs(angle_wrap(g2 - g_ang));
    err_f = abs(angle_wrap(f2 - f_ang));
    err4 = max([err_a err_e err_I err_h err_g err_f]);
    max_err_4 = max(max_err_4, err4);
    if err4 > tol_identity, n_fail_4 = n_fail_4 + 1; end
  end

  printf('--- Results ---\n');
  printf('Check 1 (vis-viva identity):        max rel err = %.3e   fails: %d/%d\n', max_err_1, n_fail_1, N);
  printf('Check 2 (h^2 = mu*a*(1-e^2)):       max rel err = %.3e   fails: %d/%d\n', max_err_2, n_fail_2, N);
  printf('Check 3 (|e_vec|^2 = 1 + 2E h^2/mu^2): max rel err = %.3e   fails: %d/%d\n', max_err_3, n_fail_3, N);
  printf('Check 4 (round-trip (r,xi) -> elements): max err = %.3e   fails: %d/%d\n', max_err_4, n_fail_4, N);

  if n_fail_1 + n_fail_2 + n_fail_3 + n_fail_4 == 0
    printf('\nALL CHECKS PASS (tolerance %.1e)\n', tol_identity);
  else
    printf('\n*** FAILURES DETECTED ***\n');
  end

end

function [r_vec, xi_vec] = elements_to_state(mu, a, e, I, h, g, f)
  % Layer-1: Kepler's equation solution via the true-eccentric half-angle bijection
  E = 2 * atan2(sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2));
  kappa = 1 - e*cos(E);
  r = a * kappa;
  n_motion = sqrt(mu / a^3);        % mean motion
  E_dot = n_motion / kappa;         % dE/dt (for unperturbed Kepler; used here only to
                                    % reconstruct a Kepler-consistent velocity)
  % Layer-2: orbit-plane coordinates and velocities
  P_op  = a*(cos(E) - e);
  Q_op  = a*sqrt(1-e^2)*sin(E);
  dPdE  = -a*sin(E);
  dQdE  =  a*sqrt(1-e^2)*cos(E);
  Pdot  = dPdE * E_dot;
  Qdot  = dQdE * E_dot;

  % Layer-3: Euler rotation R = R3(h) R1(I) R3(g)
  c_g = cos(g); s_g = sin(g);
  c_I = cos(I); s_I = sin(I);
  c_h = cos(h); s_h = sin(h);

  R3g = [c_g -s_g 0; s_g c_g 0; 0 0 1];
  R1I = [1 0 0; 0 c_I -s_I; 0 s_I c_I];
  R3h = [c_h -s_h 0; s_h c_h 0; 0 0 1];
  R = R3h * R1I * R3g;

  r_vec = R * [P_op; Q_op; 0];
  xi_vec = R * [Pdot; Qdot; 0];
end

function [a, e, I, h_ang_angle, g_ang, f_ang, h_ang_vec, e_vec, h_ang_mag] = state_to_elements(mu, r_vec, xi_vec)
  % Forward osculating-element map per Theorem K6.
  r_mag = norm(r_vec);
  xi_mag = norm(xi_vec);

  % Conserved vectors (file 03 definitions)
  h_ang_vec = cross(r_vec, xi_vec);
  h_ang_mag = norm(h_ang_vec);
  E_spec = 0.5*xi_mag^2 - mu/r_mag;
  e_vec = cross(xi_vec, h_ang_vec)/mu - r_vec/r_mag;
  e = norm(e_vec);

  % Step 1: a
  a = -mu/(2*E_spec);

  % Step 2: h (magnitude already computed)

  % Step 3: e already computed

  % Step 4: I, h_ang_angle (= RAAN)
  I = acos(h_ang_vec(3) / h_ang_mag);
  h_ang_angle = atan2(h_ang_vec(1), -h_ang_vec(2));
  if h_ang_angle < 0, h_ang_angle = h_ang_angle + 2*pi; end

  % Orbital-plane frame
  e3 = [0; 0; 1];
  n_node = cross(e3, h_ang_vec);
  n_node_mag = norm(n_node);
  p_hat = n_node / n_node_mag;
  h_hat = h_ang_vec / h_ang_mag;
  q_hat = cross(h_hat, p_hat);

  % Step 5: g (argument of periapsis) via atan2 in orbital-plane frame
  ep_p = dot(e_vec, p_hat) / e;
  ep_q = dot(e_vec, q_hat) / e;
  g_ang = atan2(ep_q, ep_p);
  if g_ang < 0, g_ang = g_ang + 2*pi; end

  % Step 6: f (true anomaly) via unwrapped E→f half-angle formula.
  % Reference: skills/OCTAVE_VERIFICATION.md rules #24-26.  The naive
  % atan2(sin_f, cos_f) from orbital-plane projections wraps to (-pi, pi],
  % inserting a 2*pi jump at f = pi that silently breaks any linear-in-f
  % downstream use (f - l + e*sin(f), f - l, weighted sums, etc.).
  % Classification of uses in this verifier (skill #28): the round-trip
  % comparison (Check 4) uses mod(2*pi), so is wrap-immune; but to follow
  % best practice and match the elements_to_state direction (which uses
  % the unwrapped formula), we reconstruct f via E.
  %
  %   cos(f) = (r . e_vec) / (r_mag * e)
  %   sin(f) = ((e_vec x r) . h_hat) / (r_mag * e)
  %     (sign: in orbital-plane basis r_hat = cos(u)p + sin(u)q,
  %      e_hat = cos(g)p + sin(g)q with p x q = h_hat, one computes
  %      e_hat x r_hat = sin(u - g) h_hat = sin(f) h_hat.  The
  %      opposite ordering r x e_vec would give -sin(f) and silently
  %      flips f into the wrong quadrant on half the test domain.)
  %   cos(E) = (e + cos(f)) / (1 + e*cos(f))
  %   sin(E) = sqrt(1-e^2) * sin(f) / (1 + e*cos(f))
  %   E in [0, 2*pi) via atan2 + shift
  %   f = 2*atan2( sqrt(1+e)*sin(E/2),  sqrt(1-e)*cos(E/2) )  in [0, 2*pi].
  cos_f = dot(r_vec, e_vec) / (r_mag * e);
  sin_f = dot(cross(e_vec, r_vec), h_hat) / (r_mag * e);
  factor = 1 + e*cos_f;
  cos_E = (e + cos_f) / factor;
  sin_E = sqrt(1 - e^2) * sin_f / factor;
  E_raw = atan2(sin_E, cos_E);
  if E_raw < 0, E_raw = E_raw + 2*pi; end
  E_anom = E_raw;  % in [0, 2*pi)
  f_ang = 2 * atan2( sqrt(1+e)*sin(E_anom/2), sqrt(1-e)*cos(E_anom/2) );
  if f_ang < 0, f_ang = f_ang + 2*pi; end
end
