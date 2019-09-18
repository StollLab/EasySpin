% Frequency dependence of nitroxide slow-motional spectra
%===============================================================================
clear, clf

% Parameters
%-------------------------------------------------------------------------------
% Typical parameters of a nitroxide radical
Sys.g = [2.008, 2.006, 2.003];
Sys.Nucs = '14N';
Sys.A = [20, 20, 90]; % MHz
Sys.tcorr = 1e-9; % seconds

% Simulations
%-------------------------------------------------------------------------------
% Loop over several spectrometer frequencies and simulate field-sweep spectra
mw = [3, 9.5, 35, 95];     % GHz
for imw = 1:4
  Experiment.mwFreq = mw(imw);
  [B{imw},spc{imw}] = chili(Sys, Experiment);
end

% Plotting
%-------------------------------------------------------------------------------
for imw = 1:4
  subplot(4,1,imw);
  plot(B{imw},spc{imw});
  axis tight;
  title(sprintf('%g GHz',mw(imw)));
end
