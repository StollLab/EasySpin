% Slow-motional cw EPR spectrum of a nitroxide radical with an orientational potential
%==========================================================================
clear, clf

% Parameters
%--------------------------------------------------------------------------
% Spin parameters
Nitroxide.g = [2.008,2.006,2.003];
Nitroxide.Nucs = '14N';
Nitroxide.A = [20,20,85];

% Dynamic parameters
Nitroxide.lw = 0.3;
Nitroxide.tcorr = 5e-9;

% Orientational potential
Nitroxide.Potential = [2 0 0 1.2]; % lambda_(2,0,0) = +1.2 

% Experimental parameters
Experiment.mwFreq = 9.5;


% Simulation and graphical rendering
%--------------------------------------------------------------------------
chili(Nitroxide,Experiment);
