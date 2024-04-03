% Spin-correlated radical pairs with different exchange and 
% dipolar interactions
%===========================================================
% Illustration of how the spectrum of a singlet-born spin-
% correlated radical pair varies with changes in the
% strength of the exchange and dipolar coupling
% 
% see
%  - Kraffert and Behrends, Mol. Phys. 115, 2373–2386 (2017)
%    https://doi.org/10.1080/00268976.2016.1278479
%    (Fig. 4)
%

clear; clc; clf;

Sys.S = [1/2 1/2];
Sys.g = [2.0023; 2.0002];

J = [0 -1 -20 0 0 -10]; % MHz
d = [0 0 0 -1 -5 -10]; % MHz

Sys.lwpp = 0.04; % mT

% Set up the initial density matrix (in the uncoupled basis of EasySpin)
S = cgmatrix(Sys.S(1),Sys.S(2),0).';  % singlet state
Sys.initState = S*S';

Exp.Range = [346.5 349.5]; % mT
Exp.mwFreq = 9.75; % GHz
Exp.nPoints = 4096;
Exp.Harmonic = 0;

% Plot individual transitions separately
Opt.separate = 'transitions';  

for i = 1:numel(J)
  
  Sys.ee = J(i) + [1 1 -2]*d(i); % MHz
      
  [B,spc,info] = pepper(Sys,Exp,Opt);
  tr = info.Transitions;
  
  subplot(2,3,i)
  hold on; box on;
  title(sprintf('J = %1.0f MHz, D = %1.0f MHz',J(i),d(i)))
  plot(B,spc,'--')
  plot(B,sum(spc),'k','LineWidth',1)
  xlim(Exp.Range)
  xlabel('{\itB} (mT)')
  
  if i==1
    nElStates = prod(2*Sys.S+1);
    leg = [ num2str(tr(:,1)), repmat(' < - > ',nElStates,1),num2str(tr(:,2))];
    legend(leg)
  end
  
end

