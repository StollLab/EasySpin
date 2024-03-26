% Temperature effects for cw EPR
%==========================================================================
clear, clf

% S=1 spin system and an experiment at Q band
Sys.S = 1;
Sys.g = 2;
Sys.D = [-600 0];  % MHz
Sys.lw = 3;  % mT

Exp.mwFreq = 35;  % GHz
Exp.Range = [1210 1290];  % mT

% Now we vary the temperature
Temp = [5 2 1]; % temperatures in kelvin
for k = 1:length(Temp)
  Exp.Temperature = Temp(k);
  [B,spc_] = pepper(Sys,Exp);
  spc(k,:) = spc_/sum(cumsum(spc_));  % normalize to double integral
end

% Plotting
plot(B,spc);
title('Temperature effects (spectra normalized to double integral)');
xlabel('magnetic field (mT)');
legend('5 K','2 K','1 K');
legend boxoff
