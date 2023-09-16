% Gd complex: individual transitions
%=============================================================
% Simulation of the EPR spectrum of a Gd complex at W-band
% and visualization of the contributions of individual
% transitions.

clear, clc, clf

Sys.S = 7/2;
Sys.g = 1.992;
Sys.D = 1200;  % D, in MHz
Sys.DStrain = 200;  % FWHM of Gaussian distribution of D, in MHz
Sys.lwpp = 2;  % mT

Exp.mwFreq = 95;  % GHz
Exp.Range = [2900 3900];  % mT
Exp.Temperature = 80;  % K
Exp.Harmonic = 0;

% Simulation of individual transitions contributing to absorption spectrum
Opt.separate = 'transitions';
[B,spec0,info] = pepper(Sys,Exp,Opt);

% Plotting
nTransitions = size(info.Transitions,1);
%colororder(parula(nTransitions+1))
plot(B,spec0)
xlabel('magnetic field (mT)')
legend([num2str(info.Transitions(:,1)), repmat(' \leftrightarrow ',nTransitions,1),num2str(info.Transitions(:,2))])
legend boxoff
