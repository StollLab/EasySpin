% Biphenyl radical anion
%==========================================================================
% Hyperfine values taken from Gerson book, p. 114

clear, clf

Sys.g = 2;
Sys.Nucs = '1H,1H,1H';
Sys.n = [2 4 4];
Sys.lwpp = [0,0.02];  % mT
Sys.A = [-15.10 -7.59 1.10];  % MHz

Exp.mwFreq = 9.5;  % GHz

garlic(Sys,Exp);
