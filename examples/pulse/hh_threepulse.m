% Three-pulse ESEEM of two 1H nuclei (saffron)
%==========================================================================

clear, clf, clc

Exp.Sequence = '3pESEEM';
Exp.Field = 352.3;
Exp.dt = 0.010;
Exp.tau = 0.1;
Exp.T = 0.06;

Sys.Nucs = '1H,1H';
Sys.A = [2, -1; 5, 3];

saffron(Sys,Exp);
