% Excited quintet state EPR of nitroxide-anthracene-nitroxide 
%===============================================================================

% Y.Teki et al.
% J.Am.Chem.Soc. 122, 5, 984-985
% https://doi.org/10.1021/ja993151e
% Figure 2

clear, clc

Sys.S = 2;
Sys.g = 2.0043;
Sys.D = 0.013*30e3;  % cm^-1 -> MHz
Sys.initState = {[0 0.35 0.3 0.35 0],'uncoupled'};
Sys.lwpp = 0.7;  % mT

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;

pepper(Sys,Exp);
