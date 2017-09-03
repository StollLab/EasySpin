clear, clc

% Normal fields
Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.logDiff = 8;

% Custom fields
Sys.AX = mt2mhz([5, 36]/10);  % A-tensor values for X-band
Sys.AQ = mt2mhz([5, 36]/10);  % Q-band
Sys.lwX = [0, 0.1];  % linewidth values for X-band
Sys.lwQ = [0, 0.1];  % Q-band

Exp.CenterSweepX = [335.6, 12];
Exp.CenterSweepQ = [1214.0, 20];
Exp.nPointsX = 1024;
Exp.nPointsQ = 1024;
Exp.mwFreqX = 9.4;
Exp.mwFreqQ = 34;

% Normal vary
Vary.logDiff = 0.5;

% Custom vary
Vary.AX = mt2mhz([0.05, 0.15]);
Vary.AQ = mt2mhz([0.05, 0.15]);
Vary.lwX = [0, 0.05];
Vary.lwQ = [0, 0.05];

% Options for spectral simulation function
SimOpt.Verbosity = [];

% Options for esfit algorithms
FitOpt.PopulationSize = 500;
FitOpt.maxGenerations = 100;
FitOpt.nParticles = 100;

% Simulate X-band data
Sys.A = Sys.AX;
Sys.lw = Sys.lwX;
Exp.CenterSweep = Exp.CenterSweepX;
Exp.mwFreq = Exp.mwFreqX;

yX = chili(Sys, Exp, SimOpt);
yX = yX/max(yX);

Sys.dataX = yX;  % store data in custom field for rescaling within custom 
                 % fitting function

% Simulate Q-band data
Sys.A = Sys.AQ;
Sys.lw = Sys.lwQ;
Exp.CenterSweep = Exp.CenterSweepQ;
Exp.mwFreq = Exp.mwFreqQ;

yQ = chili(Sys, Exp, SimOpt);
yQ = yQ/max(yQ);

Sys.dataQ = yQ; 

% Combine data sets
y = [yX, yQ];

% Add some noise
y = y + 0.05*randn(1,numel(y));

% Perform Fitting
esfit('chili_multifreq', y, Sys, Vary, Exp, SimOpt, FitOpt)