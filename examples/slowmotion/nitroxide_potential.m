% Nitroxide slow-motional cw EPR spectrum with an orientational potential
%===============================================================================
clear, clf, clc

% Parameters
%-------------------------------------------------------------------------------
% Spin parameters
Nitroxide.g = [2.008,2.006,2.003];
Nitroxide.Nucs = '14N';
Nitroxide.A = [20,20,90]; % MHz

% Dynamic parameters
Nitroxide.lw = 0.3;
Nitroxide.tcorr = 5e-9;

% Experimental parameters
Experiment.mwFreq = 9.5; % GHz


% Simulation with and without orientational potential
%-------------------------------------------------------------------------------
Nitroxide.Potential = []; % no potential
[B,spc0] = chili(Nitroxide,Experiment);

Nitroxide.Potential = [2 0 0 1.2]; % lambda_(2,0,0) = +1.2 (in units of kB*T)
Options.LLMK = [4 0 2 2]; % basis set information
[B,spc1] = chili(Nitroxide,Experiment, Options);


% Plotting
%-------------------------------------------------------------------------------
plot(B,spc0,B,spc1);
legend('no potential','potential');
xlabel('magnetic field (mT)');
axis tight
