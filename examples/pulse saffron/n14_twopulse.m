% Two-pulse ESEEM of a 14N nucleus (saffron)
%==========================================================================
% example from Flanagan/Singel, J. Chem. Phys. 87(10), 5606-5616, 1987
%   https://doi.org/10.1063/1.453532
% (see Fig. 4(a)).

clear, clf, clc

% Two-pulse experiment parameters
Exp.Field = 324.9;           % mT
Exp.Sequence = '2pESEEM';
Exp.dt = 0.100;              % µs
Exp.nPoints = 501; 

% Spin system with S=1/2 and one 14N
Sys.Nucs = '14N';
nuI = larmorfrq(Sys.Nucs,Exp.Field); % MHz
Sys.A = 2*nuI;                       % MHz
Sys.Q = [4*0.1*nuI, 0.6];            % MHz

Opt.GridSize = 91;

saffron(Sys,Exp,Opt);
