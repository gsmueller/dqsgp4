% verify_utils.m
%
% *** NOTE: the utilities are now in separate function files  ***
% *** for Octave compatibility.  Use `addpath` and call them ***
% *** directly:                                                ***
%
%   addpath('design/derivations');
%   tf = vanishes_simple(expr, e, eta);
%   tf = vanishes_eta_only(expr, e, eta);
%   tf = vanishes_full(expr, e, eta, c, s, c2g, s2g);
%   report_check('name of check', tf);
%
% Background: Octave's MATLAB-compatibility requires each top-level
% function to live in its own file (named after the function).  A file
% with multiple function definitions is treated as a script with local
% functions -- those local functions are not callable from outside
% the script.  Hence the split into:
%   vanishes_simple.m, vanishes_eta_only.m, vanishes_full.m, report_check.m
%
% This stub file exists for documentation / historical reference.

fprintf('verify_utils.m is a documentation stub only.\n');
fprintf('Use the individual .m files: vanishes_simple, vanishes_eta_only, vanishes_full, report_check.\n');
