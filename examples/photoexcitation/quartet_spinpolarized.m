% Photoexcited quartet state simulations
% ================================================================================
% 
% Simulation of a photoexcited quartet state populated by ISC between the
% 'trip-doublet' and 'trip-quartet' states
%
% see Fig. 1 in Kandrashkin, Y.E., Van Der Est, A. 
% Electron spin polarization of the excited quartet state of strongly 
% coupled triplet-doublet spin systems.
% J. Chem. Phys. 120, 4790-4799 (2004).
% https://doi.org/10.1063/1.1645773
%

clear; clc; clf;

% Spin system
Sys.S = 3/2;
Sys.g = 2.001;
Sys.lwpp = 2; % mT

Sys.Pop = [0 0 0.5 0.5];

% Experimental parameters
Exp.mwFreq = 9.75; % GHz
Exp.Range = [315 380]; % mT
Exp.Harmonic = 0;

Opt.Output = 'separate';

% E = 0
Sys.D = mt2mhz(30/3)*[1 0]; % MHz
[B,sim{1}] = pepper(Sys,Exp,Opt);
label{1} = 'E = 0';

% E = D/3
Sys.D = mt2mhz(30/3)*[1 1/3]; % MHz
[B,sim{2}] = pepper(Sys,Exp,Opt);
label{2} = 'E = 1/3 D';

% Plot individual transitions and sum
for i = 1:2
  subplot(2,1,i)
  title(label{i})
  hold on; box on;
  plot(B,sim{i})
  plot(B,sum(sim{i}),'k','LineWidth',1)
  legend('-3/2 -> -1/2','-1/2 -> +1/2','+1/2 -> +3/2','sum','Location','NorthWest')
  xlabel('B (mT)')
  xlim(B([1 end]))
end
