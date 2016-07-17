% Effect of tumbling on an axial g tensor
%===================================================================
clear, clf

% Spin system and experimental parameters
Sys.g = [2.008 2.008 2.003];
Sys.lw = 0.05;
Exp.mwFreq = 9.5;
Exp.Range = [337 340];

Opt.LLKM = [30 10 10 10]; % need a larger basis for slow tumbling

% Simulate spectra for various correlation times
tcorr = [1 3 10 30 100 300]*1e-9;
for k=1:numel(tcorr)
  Sys.tcorr = tcorr(k);
  [x,y(k,:)] =  chili(Sys,Exp,Opt);
end

% Plot simulated spectra
stackplot(x,y);
