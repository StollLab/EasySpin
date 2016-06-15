function [err,data] = test(opt,olddata)

%==================================================================
% Test whether the different input format for D work
%==================================================================

% one electron, E == 0
%------------------------------------------------------------------
Sys.S = 3/2;

D = 100*rand;

Sys.D = D;
H1 = zfield(Sys);
Sys.D = [D 0];
H2 = zfield(Sys);
Sys.D = [-1,-1,2]/3*D;
H3 = zfield(Sys);

% test
threshold = 1e-7;
err(1) = ~areequal(H1,H2,threshold);
err(2) = ~areequal(H1,H3,threshold);

data = [];
