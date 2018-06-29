% Two-pulse ESEEM of a 14N nucleus (saffron)
%==========================================================================
% example from Flanagan/Singel, J. Chem. Phys. 87(10), 5606-5616
% Analysis of 14N ESEEM patterns of randomly oriented solids
% (see Fig. 4(a)).

clear, clf, clc

Exp.Sequence = '2pESEEM';
Exp.Field = 324.9;
Exp.dt = 0.100;
Exp.nPoints = 1001;

Sys.Nucs = '14N';
nuI = larmorfrq(Sys.Nucs,Exp.Field);
Sys.A = 2*nuI;
Sys.Q = [4*0.1*nuI, 0.6];

Opt.nKnots = 31;

saffron(Sys,Exp,Opt);
