% Nitroxide radical, various correlation times
%--------------------------------------------------

clear, clf, clc

% Parameters
%---------------------------------------------------------------------
g = [2.0088,2.0061,2.0027];
A = mt2mhz([5.8,5.8,30.8]/10); % G -> MHz

Sys = struct('g',g,'Nucs','14N','A',A);
Sys.lw = [0 0.1]; % Lorentzian, mT

% Simulation
%---------------------------------------------------------------------
logtcorr = [-10.5 -10 -9.5 -9]; % £= log10(tcorr/seconds)

Exp.mwFreq = 9.5;
Exp.Range = [336 341];

for k = 1:numel(logtcorr)
  Sys.logtcorr = logtcorr(k);
  [BX,spcX(k,:)] = garlic(Sys,Exp);
end

% Plotting
%---------------------------------------------------------------------
stackplot(BX,spcX);
