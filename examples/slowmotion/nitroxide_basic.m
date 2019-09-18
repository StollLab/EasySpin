% Simple slow-motional cw EPR spectrum simulation of a nitroxide radical
%===============================================================================
clear, clf

% Parameters
%-------------------------------------------------------------------------------
% Static parameters
Nitroxide.g = [2.008,2.006,2.003];
Nitroxide.Nucs = '14N';
Nitroxide.A = [20,20,90]; % MHz

% Dynamic parameters
Nitroxide.lw = 0.3; % mT
Nitroxide.tcorr = 5e-9; % seconds

% Experimental parameters
Experiment.mwFreq = 9.5; % GHz


% Simulation and graphical rendering
%-------------------------------------------------------------------------------
chili(Nitroxide,Experiment);
