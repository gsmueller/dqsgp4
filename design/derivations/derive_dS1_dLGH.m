%% derive_dS1_dLGH.m
%% Differentiate S1 (Brouwer 1959 Eq. 15) w.r.t. L, G, H symbolically.
%%
%% S1 = (mu^2 k2 / G^3) * { A*(f - l + e sin f)
%%                          + B*[1/2 sin(2g+2f) + (e/2) sin(2g+f)
%%                               + (e/6) sin(2g+3f)] }
%%   A = -1/2 + 3/2*theta^2, B = 3/2 - 3/2*theta^2, theta = H/G
%%
%% S1 depends on (L, G, H, l, g) through:
%%   - 1/G^3 (explicit in G)
%%   - A, B via theta = H/G (explicit in G, H)
%%   - e via e = sqrt(1 - G^2/L^2)  (implicit in L, G)
%%   - f via Kepler, f = f(l, e)     (implicit in L, G through e)
%%   - explicit l, g appear
%%
%% Chain rule (fixed l, g, and the other actions):
%%   dS1/dL|_{G,H,l,g} = (dS1/de)|_explicit * de/dL|_G
%%                     + (dS1/df)|_explicit * df/dL|_l,G
%%
%%   dS1/dG|_{L,H,l,g} = (dS1/dG)|_explicit         (1/G^3 and theta)
%%                     + (dS1/de)|_explicit * de/dG|_L
%%                     + (dS1/df)|_explicit * df/dG|_l,L
%%
%%   dS1/dH|_{L,G,l,g} = (dS1/dH)|_explicit         (theta only)
%%
%% Primitives from derive_anomaly_partials.m (VERIFIED):
%%   de/dL|_G    = G^2/(L^3*e)
%%   de/dG|_L    = -G/(L^2*e)
%%   df/de|_l    = sin(f)*(2+e cos f)/(1-e^2)
%%   df/dL|_l,G  = G^2 * sin(f)*(2+e cos f) / (L^3 * e * (1-e^2))
%%   df/dG|_l,L  = -G * sin(f)*(2+e cos f) / (L^2 * e * (1-e^2))

pkg load symbolic;

printf('============================================================\n');
printf('DIFFERENTIATING S1 (Brouwer Eq. 15) w.r.t. L, G, H, and l\n');
printf('============================================================\n\n');

syms mu k2 real positive;
syms L G H real positive;
syms l_sym g e f real;

theta = H/G;
A_expr = -sym(1)/2 + sym(3)/2 * theta^2;
B_expr =  sym(3)/2 - sym(3)/2 * theta^2;

% Build S1 with e, f as independent variables (chain rule applied manually):
S1_core = A_expr*(f - l_sym + e*sin(f)) ...
        + B_expr*(sym(1)/2*sin(2*g+2*f) + (e/2)*sin(2*g+f) + (e/6)*sin(2*g+3*f));
S1 = (mu^2*k2/G^3) * S1_core;

% Primitives:
de_dL   = G^2/(L^3*e);
de_dG   = -G/(L^2*e);
df_de_l = sin(f)*(2 + e*cos(f))/(1-e^2);
df_dL   = G^2 * sin(f)*(2+e*cos(f)) / (L^3*e*(1-e^2));
df_dG   = -G * sin(f)*(2+e*cos(f)) / (L^2*e*(1-e^2));
% df/dl|_e = eta*(a/r)^2; at fixed e.  With eta = sqrt(1-e^2),
% a/r = (1+e cos f)/(1-e^2). Keep a symbol-free form:
df_dl   = sqrt(1-e^2)*(1+e*cos(f))^2/(1-e^2)^2;
% Equivalent: df_dl = (1+e cos f)^2/(1-e^2)^{3/2}

% ----------------------------------------------------------------
% CROSS-CHECK: dS1/dl should equal the form used in delta_p1:
%   dS1/dl = (mu^2 k2/L^3)*[A*(a/r)^3 - A*(L^3/G^3) + B*(a/r)^3*cos(2g+2f)]
%   which in pure (e,f,G,L) form (a/r = (1+e cos f)/(1-e^2)) is:
% ----------------------------------------------------------------
printf('--- Cross-check: compute dS1/dl and compare with delta_p1 input ---\n');
dS1_dl = diff(S1, l_sym, 1, 'e', sym(0), 'f', sym(0));  % ignore, just wanting syntax
% Actually, chain rule through explicit l, and through f:
% dS1/dl = d_lexplicit S1  +  (dS1/df) * (df/dl|_e)
dS1_dl_explicit = diff(S1, l_sym);   % derivative treating e,f as indep of l
dS1_df_explicit = diff(S1, f);
dS1_dl = dS1_dl_explicit + dS1_df_explicit * df_dl;

% Expected form:
a_over_r = (1 + e*cos(f))/(1-e^2);
expected_dS1_dl = (mu^2*k2/L^3) * ( A_expr*a_over_r^3 - A_expr*(L^3/G^3) ...
                                   + B_expr*a_over_r^3*cos(2*g+2*f) );

diff_dS1_dl = simplify(expand(dS1_dl - expected_dS1_dl));
printf('  dS1/dl - [expected form used in delta_p1] = ');
disp(diff_dS1_dl);
printf('  (Should be 0 to confirm our Kepler chain is correct)\n\n');

% ----------------------------------------------------------------
% dS1/dL: only through e and f
% ----------------------------------------------------------------
printf('--- dS1/dL ---\n');
dS1_de_expl = diff(S1, e);
dS1_df_expl = diff(S1, f);

dS1_dL = dS1_de_expl * de_dL + dS1_df_expl * df_dL;

% Decompose by harmonic. The harmonic structure of S1 is:
%   constant (no angle) : A*(f - l + ...) has f explicit (non-harmonic); e sin f harmonic.
% Actually, expanding, the harmonics in S1_core (as functions of f,g):
%   sin(f)      : A*e*sin(f)
%   f, -l       : A*(f - l)   <- secular, no sin/cos of multiple of f,g
%   sin(2g+2f)  : B/2
%   sin(2g+f)   : B*e/2
%   sin(2g+3f)  : B*e/6
%
% We will NOT try to simplify the combined dS1/dL as a single blob.
% Instead, split by the harmonic they originate from:
%   Piece P0 = A*(f - l + e*sin(f)) *mu^2 k2/G^3
%   Piece Q1 = B*(1/2)*sin(2g+2f) *mu^2 k2/G^3
%   Piece Q2 = B*(e/2)*sin(2g+f) *mu^2 k2/G^3
%   Piece Q3 = B*(e/6)*sin(2g+3f) *mu^2 k2/G^3
%
P0 = (mu^2*k2/G^3) * A_expr*(f - l_sym + e*sin(f));
Q1 = (mu^2*k2/G^3) * B_expr*(sym(1)/2)*sin(2*g+2*f);
Q2 = (mu^2*k2/G^3) * B_expr*(e/2)*sin(2*g+f);
Q3 = (mu^2*k2/G^3) * B_expr*(e/6)*sin(2*g+3*f);

% Helper to compute dX/dL for any piece X:
dX_dL = @(X) diff(X, e)*de_dL + diff(X, f)*df_dL;
dX_dG = @(X) diff(X, G) + diff(X, e)*de_dG + diff(X, f)*df_dG;
dX_dH = @(X) diff(X, H);

printf('Piece P0 (secular + sin f): applying chain rule...\n');
dP0_dL = simplify(dX_dL(P0));
printf('  dP0/dL = '); disp(dP0_dL);

printf('Piece Q1 (sin(2g+2f)):\n');
dQ1_dL = simplify(dX_dL(Q1));
printf('  dQ1/dL = '); disp(dQ1_dL);

printf('Piece Q2 (sin(2g+f)):\n');
dQ2_dL = simplify(dX_dL(Q2));
printf('  dQ2/dL = '); disp(dQ2_dL);

printf('Piece Q3 (sin(2g+3f)):\n');
dQ3_dL = simplify(dX_dL(Q3));
printf('  dQ3/dL = '); disp(dQ3_dL);

dS1_dL_total = dP0_dL + dQ1_dL + dQ2_dL + dQ3_dL;

% ----------------------------------------------------------------
% dS1/dG: through 1/G^3 prefactor (explicit), theta=H/G in A,B,
%         through e, through f
% ----------------------------------------------------------------
printf('\n--- dS1/dG ---\n');
printf('Piece P0:\n');
dP0_dG = simplify(dX_dG(P0));
printf('  dP0/dG = '); disp(dP0_dG);

printf('Piece Q1:\n');
dQ1_dG = simplify(dX_dG(Q1));
printf('  dQ1/dG = '); disp(dQ1_dG);

printf('Piece Q2:\n');
dQ2_dG = simplify(dX_dG(Q2));
printf('  dQ2/dG = '); disp(dQ2_dG);

printf('Piece Q3:\n');
dQ3_dG = simplify(dX_dG(Q3));
printf('  dQ3/dG = '); disp(dQ3_dG);

dS1_dG_total = dP0_dG + dQ1_dG + dQ2_dG + dQ3_dG;

% ----------------------------------------------------------------
% dS1/dH: only through theta = H/G in A, B
% ----------------------------------------------------------------
printf('\n--- dS1/dH ---\n');
printf('Piece P0:\n');
dP0_dH = simplify(dX_dH(P0));
printf('  dP0/dH = '); disp(dP0_dH);

printf('Piece Q1:\n');
dQ1_dH = simplify(dX_dH(Q1));
printf('  dQ1/dH = '); disp(dQ1_dH);

printf('Piece Q2:\n');
dQ2_dH = simplify(dX_dH(Q2));
printf('  dQ2/dH = '); disp(dQ2_dH);

printf('Piece Q3:\n');
dQ3_dH = simplify(dX_dH(Q3));
printf('  dQ3/dH = '); disp(dQ3_dH);

dS1_dH_total = dP0_dH + dQ1_dH + dQ2_dH + dQ3_dH;

% ----------------------------------------------------------------
% Save expressions to a reusable .m file via sympref
% ----------------------------------------------------------------
printf('\n============================================================\n');
printf('SUMMARY:\n');
printf('  dS1/dL = sum of dP0/dL + dQ1/dL + dQ2/dL + dQ3/dL (per-harmonic)\n');
printf('  dS1/dG = sum of dP0/dG + dQ1/dG + dQ2/dG + dQ3/dG (per-harmonic)\n');
printf('  dS1/dH = sum of dP0/dH + dQ1/dH + dQ2/dH + dQ3/dH (per-harmonic)\n');
printf('============================================================\n');

% Export the per-harmonic expressions for downstream scripts.
% Octave symbolic does not save well to .mat; we write to a .m file.
fid = fopen('design/derivations/dS1_partials_cache.m', 'w');
fprintf(fid, '%% Auto-generated by derive_dS1_dLGH.m\n');
fprintf(fid, '%% Per-harmonic dS1 partial derivatives.\n');
fprintf(fid, '%% To use, run this after declaring the same symbols.\n\n');
fprintf(fid, 'dP0_dL = %s;\n', char(dP0_dL));
fprintf(fid, 'dQ1_dL = %s;\n', char(dQ1_dL));
fprintf(fid, 'dQ2_dL = %s;\n', char(dQ2_dL));
fprintf(fid, 'dQ3_dL = %s;\n', char(dQ3_dL));
fprintf(fid, 'dP0_dG = %s;\n', char(dP0_dG));
fprintf(fid, 'dQ1_dG = %s;\n', char(dQ1_dG));
fprintf(fid, 'dQ2_dG = %s;\n', char(dQ2_dG));
fprintf(fid, 'dQ3_dG = %s;\n', char(dQ3_dG));
fprintf(fid, 'dP0_dH = %s;\n', char(dP0_dH));
fprintf(fid, 'dQ1_dH = %s;\n', char(dQ1_dH));
fprintf(fid, 'dQ2_dH = %s;\n', char(dQ2_dH));
fprintf(fid, 'dQ3_dH = %s;\n', char(dQ3_dH));
fclose(fid);
printf('\nCached per-harmonic partials to design/derivations/dS1_partials_cache.m\n');
