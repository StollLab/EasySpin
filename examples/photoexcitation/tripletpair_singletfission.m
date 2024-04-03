% Triplet state and pair of strongly-exchange coupled 
% triplet states formed from singlet fission
%===========================================================
% Simulation of the EPR spectra of triplet states and strongly-exchange coupled
% pairs of triplet states (quintet state) formed from singlet fission
% 
% see
%  - Weiss, L. R., Bayliss, S. L., Kraffert, F., Thorley, K. J., Anthony, J. E.,
%    Bittl, R., Friend, R. H., Rao, A., Greenham, N. C., Behrends, J.
%    Nat. Phys. 13, 176–181 (2017)
%    https://doi.org/10.1038/nphys3908
%

clear; clc; clf;

% Single recombination triplet
% -----------------------------------------------------------
SysT.S = 1;
D = 1400; % MHz
SysT.D = D; % MHz
SysT.lwpp = 2; % mT

% Selective population of triplet mS=0 (T0) level
SysT.initState = {[0 1 0],'eigen'};

% Pair of strongly coupled triplet states
% -----------------------------------------------------------
SysTT.S = [1 1];
SysTT.D = [D 0; D 0]; % MHz
SysTT.J = 1e5;  % MHz
SysTT.lwpp = SysT.lwpp; % mT

% Selective population of quintet mS=0 level
%              S = [0  1 1 1  2  2 2 2 2]
%             mS = [0 -1 0 1 -2 -1 0 1 2]
SysTT.initState = {[0  0 0 0  0  0 1 0 0],'eigen'};

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [270 405]; % mT
Exp.Harmonic = 0;

Opt.separate = 'transitions';

% Calculate spectra
[B,spcT,infoT] = pepper(SysT,Exp,Opt);
[B,spcTT,infoTT] = pepper(SysTT,Exp,Opt);
trT = infoT.Transitions;
trTT= infoTT.Transitions;

% Plotting
subplot(2,2,1)
title('Triplet state')
hold on; box on;
plot(B,spcT)
plot(B,sum(spcT),'k')
xlim(Exp.Range)
xlabel('magnetic field (mT)')
leg = [num2str(trT(:,1)), repmat(' < - > ',size(trT,1),1),num2str(trT(:,2))];
legend(leg,'Location','NorthWest')

subplot(2,2,2)
levelsplot(SysT,'z',[0 400],Exp.mwFreq)

subplot(2,2,3)
title('Pair of exchange-coupled triplet states')
hold on; box on;
plot(B,spcTT)
plot(B,sum(spcTT),'k')
xlim(Exp.Range)
xlabel('magnetic field (mT)')
leg = [num2str(trTT(:,1)), repmat(' < - > ',size(trTT,1),1),num2str(trTT(:,2))];
legend(leg,'Location','NorthWest')

subplot(2,2,4)
levelsplot(SysTT,'z',[0 400],Exp.mwFreq)
