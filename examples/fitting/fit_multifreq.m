% Using a custom function for globally fitting multiple EPR spectra

clear, clc

% Normal fields (these are typically used by existing EasySpin functions)
Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = unitconvert([0.5, 3.6],'mT->MHz');
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
% Vary.A = unitconvert([0.1, 0.15],'mT->MHz');
Vary.logDiff = 0.5;

% Custom vary (these vary the custom fields)
Vary.lwX = [0, 0.05];
Vary.lwQ = [0, 0.05];

% Simulate mock X-band data
Sys.lw = Sys.lwX;
Exp.CenterSweep = Exp.CenterSweepX;
Exp.mwFreq = Exp.mwFreqX;

[B{1},expdata{1}] = chili(Sys,Exp);

% Simulate mock Q-band data
Sys.lw = Sys.lwQ;
Exp.CenterSweep = Exp.CenterSweepQ;
Exp.mwFreq = Exp.mwFreqQ;

[B{2},expdata{2}] = chili(Sys,Exp);

% Add some noise
SNR = 20;
for i = 1:numel(expdata)
  expdata{i} = addnoise(expdata{i},SNR,'n');
end

% Perform Fitting
FitOpt.Method = 'simplex fcn';
result = esfit(expdata, @chili_multifreq, {Sys, Exp}, {Vary}, FitOpt);

% Plot results
for k = 1:numel(expdata)
  subplot(2,1,k)
  x = 1:numel(expdata{k});
  plot(B{k},expdata{k},B{k},result.fit{k});
end

% Run fitting in GUI
FitOpt.x = B;
esfit(expdata, @chili_multifreq, {Sys, Exp}, {Vary}, FitOpt);


% Custom function for simulating slow-motion X- and Q-band EPR spectra
% ---------------------------------------------------------------------------------
function y = chili_multifreq(Sys,Exp)

% X-band spectrum
Sys.lw = Sys.lwX;

Exp.CenterSweep = Exp.CenterSweepX;
Exp.nPoints = Exp.nPointsX;
Exp.mwFreq = Exp.mwFreqX;

y{1} = chili(Sys,Exp);

% Q-band spectrum
Sys.lw = Sys.lwQ;

Exp.CenterSweep = Exp.CenterSweepQ;
Exp.nPoints = Exp.nPointsQ;
Exp.mwFreq = Exp.mwFreqQ;

y{2} = chili(Sys,Exp);

end
