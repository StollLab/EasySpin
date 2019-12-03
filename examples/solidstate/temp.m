% Temperature effects for cw EPR
%==========================================================================
clear, clf

% S=1 spin system and an experiment at Q band
Sys = struct('S',1,'g',2,'D',[-600 0],'lw',3);
Exp = struct('Range',[1210 1290],'mwFreq',35);

% Now we vary the temperature
Temp = [5 2 1]; % temperatures in kelvin
for k = 1:length(Temp)
  Exp.Temperature = Temp(k);
  [B,spc_] = pepper(Sys,Exp);
  spc(k,:) = spc_/sum(cumsum(spc_)); % normalize to double integral
end

% Plotting
plot(B,spc);
title('Temperature effects (spectra normalized to double integral)');
xlabel('magnetic field (mT)');
legend('5 K','2 K','1 K');
legend boxoff
