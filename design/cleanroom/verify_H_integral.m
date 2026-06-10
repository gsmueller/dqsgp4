%% verify_H_integral.m
%% Numerical verification of <H>_l where
%%   H = (1+e cos f)^3 * cos(2(f+g)) * [A(theta) + B(theta) cos(2(f+g))]
%% at the five sample points used in the Gemini cross-check.
%%
%% Direct integration in l (no dl/df transform): Newton-iterate Kepler at
%% each quadrature node, recover f from l, evaluate H, periodic-trapezoidal
%% sum at N=4096 nodes (spectrally convergent on smooth periodic integrands).

A_func = @(theta) (3*theta.^2 - 1)/2;
B_func = @(theta) 3*(1 - theta.^2)/2;

%% Newton iteration for Kepler: l = E - e*sin(E)
function E = kepler_solve(l, e)
    E = l;  % initial guess
    for k = 1:50
        f_val = E - e*sin(E) - l;
        f_pr  = 1 - e*cos(E);
        dE = -f_val ./ f_pr;
        E = E + dE;
        if max(abs(dE)) < 1e-15
            break;
        end
    end
end

%% f from E using unwrapped tan-half-angle (skill OCTAVE_VERIFICATION §26)
function f = f_from_E(E, e)
    f = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    %% unwrap to make f monotonic in E (and hence in l):
    %% f and E share branches at multiples of pi, so we add the appropriate
    %% multiples of 2*pi to make f track E continuously.
    f = f + 2*pi*round((E - f)/(2*pi));
end

function val = H_integrand(l, e, theta, g)
    E = kepler_solve(l, e);
    f = f_from_E(E, e);
    A = (3*theta^2 - 1)/2;
    B = 3*(1 - theta^2)/2;
    u = 1 + e*cos(f);
    val = u.^3 .* cos(2*(f+g)) .* (A + B*cos(2*(f+g)));
end

%% Periodic trapezoidal at N nodes
N = 4096;
l_nodes = (0:N-1) * (2*pi/N);

samples = [
    0.10,  0.5,  0.30;
    0.30,  0.7,  1.10;
    0.50, -0.4,  2.40;
    0.70,  0.6,  4.70;
    0.85, -0.8,  0.00
];

printf('\n  Sample |    e    |  theta  |    g    |  <H>_l (numerical)  |  (1/2) B eta^3 (analytical)  |  Gemini reported\n');
printf('  -------+---------+---------+---------+---------------------+------------------------------+--------------------\n');

gemini_vals = [0.36444840, 0.23724300, 0.44437500, 0.29760000, 0.16593414];

for i = 1:size(samples,1)
    e = samples(i,1);
    th = samples(i,2);
    g = samples(i,3);
    eta = sqrt(1 - e^2);
    B = 3*(1 - th^2)/2;
    analytical = 0.5 * B * eta^3;

    H_vals = arrayfun(@(l) H_integrand(l, e, th, g), l_nodes);
    numerical = mean(H_vals);

    printf('   %d     |  %.2f   |  %5.2f  |  %5.2f  |   %18.12f  |    %18.12f      |    %12.8f\n', ...
           i, e, th, g, numerical, analytical, gemini_vals(i));
end

printf('\n');
