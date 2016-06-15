% Frequency dependence of nitroxide slow-motional spectra
%==========================================================================
clear, clf

% Parameters
%--------------------------------------------------------------------------
g = [2.008, 2.006, 2.003];
Sys = struct('g', g, 'Nucs', '14N', 'A', [20, 20, 85]);
Sys.tcorr = 1e-9;

% Simulations
%--------------------------------------------------------------------------
mw = [3, 9.5, 35, 95];     % GHz
for k = 1:4
  Experiment.mwFreq = mw(k);
  [B{k},spc{k}] = chili(Sys, Experiment);
end

% Graphical rendering
%--------------------------------------------------------------------------
for k = 1:4
  subplot(4,1,k);
  plot(B{k},spc{k});
  axis tight;
end
