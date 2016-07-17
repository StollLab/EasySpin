% Triplet powder spectrum, with non-equilibrium populations (simple)
%==========================================================================
clear, clf, clc

% Spin system and experimental parameters
Sys = struct('S',1,'g',2,'lw',0.3,'D',[100,0]);
Exp = struct('mwFreq',9.5,'Range',[325 355],'Harmonic',0);

% Boltzmann equilibrium at room temperature
Exp.Temperture = 298;
[x,y{1}] = pepper(Sys,Exp);

% User-specified population vector
Exp.Temperature = [0.85 1 0.95];
[x,y{2}] = pepper(Sys,Exp);

% Plot simulation results
titles{1} = 'Thermal equilibrium populations';
titles{2} = 'Non-equilibrium populations';
for k = 1:2
  subplot(2,1,k);
  plot(x,y{k});
  title(titles{k});
  xlabel('magnetic field [mT]');
end
