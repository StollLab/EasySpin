% Cu2+ with four imidazole nitrogens
%=============================================================
% Simulation of the EPR spectrum of a S=1/2 Cu2+ coordinated by
% four imidazoles, so that the couplings to four nitrogens are
% resolved.

% parameters take from S. K. Buchanan, G. C. Dismukes,
% Biochemistry 1987, 26(16), 5049-5055 (see Fig.1 and its legend)

clear
Sys.Nucs = 'Cu,14N,14N,14N,14N';
Sys.g = [2.05 2.19];
A_Cu = [0.0019 0.0203]; % in cm^-1
A_N = [0.00145 0.0017]; % in cm^-1
Sys.A = [A_Cu;A_N;A_N;A_N;A_N]*30e3; % cm^-1 -> MHz
Sys.AFrame = [0 0 0; 0 -1 0; -1 -1 0; -2 -1 0; -3 -1 0]*pi/2;
Sys.lwpp = [0 0.65];    % mT

Exp.mwFreq = 9.05;      % GHz
Exp.Range = [250 350];  % mT

Opt.Method = 'perturb';
Opt.nKnots = [61 0];

pepper(Sys,Exp,Opt);
