% Correlation time dependence of a nitroxide spectrum
%===============================================================================
clear, clf, clc

% Parameters
%-------------------------------------------------------------------------------
Nitroxide.g = [2.008, 2.006, 2.003];
Nitroxide.Nucs = '14N';
Nitroxide.A = [20, 20, 90]; % MHz
Nitroxide.lw = 0.1; % mT

Experiment.mwFreq = 9.5;   % GHz

Options.Verbosity = 1;
Options.LLKM = [24 14 6 6]; % moderately large basis size

% Simulate spectra over a range of rotational correlation times
%-------------------------------------------------------------------------------
tcorr = [0.03 0.1 0.3 1 3 10 30 100 300]*1e-9; % seconds
nSpectra = numel(tcorr);
for k = 1:nSpectra
  Nitroxide.tcorr = tcorr(k);
  [B{k},spc{k}] = chili(Nitroxide,Experiment,Options);
end

% Plotting
%-------------------------------------------------------------------------------
hold on

for k = 1:nSpectra
  spcnorm = spc{k}/max(spc{k});
  shift = k-1;
  plot(B{k},spcnorm + shift);
end

axis tight
xlabel('magnetic field [mT]');
