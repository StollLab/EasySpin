% Temperature effects for cw EPR
%==========================================================================
clear, clf

% S=1 spin system and an experiment at Q band
Sys = struct('S',1,'g',2,'D',[-600 0],'lw',3);
Exp = struct('Range',[1210 1290],'mwFreq',35);

% Now we vary the temperature, specified in the
% Experiment structure.

Temp = [5 2 1]; % temperatures in Kelvin
for k = 1:length(Temp)
  Exp.Temperature = Temp(k);
  [x,yy] = pepper(Sys,Exp);
  y(k,:) = yy/sum(cumsum(yy)); % normalize to double integral
end

% Display
plot(x,y);
title('Temperature effects (spectra normalized to double integral)');
xlabel('magnetic field [mT]');
legend('5K','2K','1K',0);
legend boxoff
