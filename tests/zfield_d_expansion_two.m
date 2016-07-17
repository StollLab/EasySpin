function [err,data] = test(opt,olddata)

%==================================================================
% Test whether the different input format for D work
%==================================================================

% one electron
%------------------------------------------------------------------
Sys.S = [3/2 1];
Sys.ee = 1;

D = 100*rand(1,2);


Sys.D = D;
H1 = zfield(Sys);
Sys.D = [D(1) 0; D(2), 0];
H2 = zfield(Sys);
Sys.D = D(:)*[-1,-1,2]/3;
H3 = zfield(Sys);

% test
threshold = 1e-7;
err(1) = ~areequal(H1,H2,threshold);
err(2) = ~areequal(H1,H3,threshold);

data = [];
