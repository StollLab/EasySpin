% ENDOR excitation widths and orientation selectivity
%=======================================================================
clear, clf

% Axial S=1/2 system with 1 proton
Sys.S = 1/2;
Sys.g = [2, 2.2];
Sys.Nucs = '1H';
Sys.A = [-1 2]*2 + 3;
Sys.lwEndor = 0.1;

% X band experimental conditions
Exp.mwFreq = 9.5;
Exp.Field = 315;
Exp.Range = [8 19];

% Set high orientational resolution for the simulation
Opt.GridSize = 91;

% Now we loop over a set of different excitation widths
% to demonstrate the effect on the ENDOR spectra.
Widths = [1e5 1000 300 100];

for k = 1:4
  Exp.ExciteWidth = Widths(k);
  [freq,spectra(k,:)] = salt(Sys,Exp,Opt);
end

plot(freq,spectra); axis tight
xlabel('frequency (MHz)');
legend('100GHz','1GHz','300MHz','100MHz');
