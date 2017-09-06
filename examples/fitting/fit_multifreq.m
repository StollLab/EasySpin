% Using a custom function for globally fitting multiple EPR spectra

clear, clc

% Normal fields (these are typically used by existing EasySpin functions)
Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([0.5, 3.6]);
Sys.logDiff = 8;

% Custom fields (these are used by the custom function chili_multifreq)
Sys.lwX = [0, 0.1];  % linewidth values for X-band
Sys.lwQ = [0, 0.1];  % Q-band

Exp.CenterSweepX = [335.6, 12];
Exp.CenterSweepQ = [1214.0, 20];
Exp.nPointsX = 1024;
Exp.nPointsQ = 1024;
Exp.mwFreqX = 9.4;
Exp.mwFreqQ = 34;

% Normal vary (these vary the normal fields)
Vary.A = mt2mhz([0.1, 0.15]);
Vary.logDiff = 0.5;

% Custom vary (these vary the custom fields)
Vary.lwX = [0, 0.05];
Vary.lwQ = [0, 0.05];

% Options for spectral simulation function
SimOpt.Verbosity = [];

% Options for esfit algorithms
FitOpt.PopulationSize = 500;
FitOpt.maxGenerations = 100;
FitOpt.nParticles = 100;

% Simulate mock X-band data
Sys.lw = Sys.lwX;
Exp.CenterSweep = Exp.CenterSweepX;
Exp.mwFreq = Exp.mwFreqX;

yX = chili(Sys, Exp, SimOpt);
yX = yX/max(yX);

Sys.dataX = yX;  % store data in custom field for rescaling within custom 
                 % simulation function

% Simulate mock Q-band data
Sys.lw = Sys.lwQ;
Exp.CenterSweep = Exp.CenterSweepQ;
Exp.mwFreq = Exp.mwFreqQ;

yQ = chili(Sys, Exp, SimOpt);
yQ = yQ/max(yQ);

Sys.dataQ = yQ; 

% Combine data sets
y = [yX, yQ];

% Add some noise
SNR = 20;
y = addnoise(y,SNR,'n');

% Perform Fitting
esfit('chili_multifreq', y, Sys, Vary, Exp, SimOpt, FitOpt)