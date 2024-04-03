% Parafluoronitrobenzene (multinuclear system), fast-motional regime
%=============================================================
% Data from Carrington, Hudson, Luckhurst,
% Proc. Roy. Soc. A 284 (1965), 582-593

clear, clf

% Parameters
%-------------------------------------------------------------
g = [2.0032 2.0012 2.0097];
A_N = 40.40 + [24,-12,-12];
A_F = 22.51 + [34.9,-19.8,-15];
A_oH = [1 1 1]*9.69;
A_mH = [1 1 1]*3.16;
lw = 0.01;
tcorr = 8e-11;

Sys.g = g;
Sys.Nucs = '14N,19F,1H,1H,1H,1H';
Sys.A = [A_N;A_F;A_oH;A_oH;A_mH;A_mH];
Sys.lw = [0 lw];
Exp.mwFreq = 9.5;
Exp.nPoints = 1e4;

% Fast-motion regime spectrum
%-------------------------------------------------------------
Sys.tcorr = tcorr;
garlic(Sys,Exp);
