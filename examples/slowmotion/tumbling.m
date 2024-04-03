% Effect of tumbling on EPR spectrum of axial g tensor
%===============================================================================
clear, clf, clc

% Spin system and experimental parameters
%-------------------------------------------------------------------------------
Sys.g = [2.008 2.008 2.003];
Sys.lw = 0.05;  % Gaussian broadening, mT
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [337 340];  % mT

Opt.LLMK = [14 3 6 6];  % use a larger basis for slow tumbling

% Simulate spectra for various correlation times
%-------------------------------------------------------------------------------
tcorr = [1 3 10 30 100 300]*1e-9; % seconds
for k = 1:numel(tcorr)
  Sys.tcorr = tcorr(k);
  [B,spc(k,:)] =  chili(Sys,Exp,Opt);
end

% Plot simulated spectra
%-------------------------------------------------------------------------------
stackplot(B,spc,'maxabs',2,compose('t_{corr} = %1.0f ns',tcorr*1e9));
xlabel('magnetic field (mT)');
