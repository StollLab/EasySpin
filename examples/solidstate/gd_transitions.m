% Gd complex: individual transitions
%=============================================================
% Simulation of the EPR spectrum of a Gd complex at W-band
% and visualization of the contributions of individual
% transitions.

clear, clc, clf

% Parameters of Gd complex
Gd.S = 7/2;
Gd.g = 1.992;
Gd.D = [1200 0];  % D and E, in MHz
Gd.DStrain = 200;  % FWHM of Gaussian distribution of D, in MHz
Gd.lwpp = 2;  % mT

% Experimental parameters
Exp.mwFreq = 95;  % GHz
Exp.Range = [2900 3900];  % mT
Exp.Temperature = 10;  % K
Exp.Harmonic = 0;

% Simulation of individual transitions contributing to absorption spectrum
Opt.separate = 'transitions';
[B,spec0,info] = pepper(Gd,Exp,Opt);

% Plotting
tr = info.Transitions;
nTransitions = size(tr,1);
plot(B,spec0)
xlabel('magnetic field (mT)')
str = [num2str(tr(:,1)), repmat(' \leftrightarrow ',nTransitions,1), num2str(tr(:,2))];
legend(str)
legend boxoff
