% Cu2+ with four imidazole nitrogens
%=============================================================
% Simulation of the EPR spectrum of a S=1/2 Cu2+ coordinated by
% four imidazoles, so that the couplings to four nitrogens are
% resolved.

% Parameters are taken from S. K. Buchanan, G. C. Dismukes,
% Biochemistry 1987, 26(16), 5049-5055 (see Fig.1 and its legend)
% https://doi.org/10.1021/bi00390a025

clear, clc, clf

CuN4.Nucs = 'Cu,14N,14N,14N,14N';
CuN4.g = [2.05 2.19];
A_Cu = [0.0019 0.0203]*30e3; % in cm^-1, converted to MHz
A_N = [0.00145 0.0017]*30e3; % in cm^-1, converted to MHz
CuN4.A = [A_Cu; A_N; A_N; A_N; A_N];  % MHz
CuN4.AFrame = [0 0 0; 0 -1 0; -1 -1 0; -2 -1 0; -3 -1 0]*pi/2;
CuN4.lwpp = [0 0.65];    % mT

Exp.mwFreq = 9.05;      % GHz
Exp.Range = [250 350];  % mT

Opt.Method = 'perturb';
Opt.GridSize = [61 0];

pepper(CuN4,Exp,Opt);
