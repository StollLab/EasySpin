% Pentacene radical anion and cation
%==========================================================================
% Values taken from Weil/Bolton/Wertz book, Table 9.3, page 249

clear, clc

% Spin center
%--------------------------------------------------------------------------
% Spin system
Sys.g = 2;
Sys.Nucs = '1H,1H,1H,1H';
Sys.n = [2 4 4 4];
Sys.lwpp = [0, 0.01];  % Lorentzian lines, mT

% Hyperfine couplings
A_anion =  unitconvert([0.4263 0.3032 0.0915 0.0870],'mT->MHz');
A_cation = unitconvert([0.5083 0.3554 0.0975 0.0757],'mT->MHz');

% Spectrometer settings
%--------------------------------------------------------------------------
Exp.mwFreq = 9.5;  % GHz
Exp.CenterSweep = [339.4, 4];  % mT
Exp.nPoints = 1e4;

% Simulations
%--------------------------------------------------------------------------
Sys.A = A_anion;
[B,spca] = garlic(Sys,Exp);
Sys.A = A_cation;
[B,spcc] = garlic(Sys,Exp);

% Graphical rendering
%--------------------------------------------------------------------------
plot(B,spca,B,spcc);
legend('radical anion','radical cation');
xlabel('magnetic field [mT]');
axis tight
title('pentacene radical spectra');
