% Spin Correlated Radical Pair
%===========================================================
% Hamiltonian parameters based on:
% https://doi.org/10.1021/acs.jpca.8b07556

clear; clc;

Sys.S = [1/2 1/2];
Sys.g = [2.01566 2.00783 2.00306;  % tetrathiafulvalene
         2.00558 2.00598 2.00338]; % pyromellitimide
       
Sys.ee = [-4.55 -4.44 10.58]; % MHz
Sys.lwpp = 0.38; % mT


Exp.Range = [335 340]; % mT
Exp.mwFreq = 9.5; % GHz
Exp.nPoints = 4096;
Exp.Harmonic = 0;

% Set up the initial density matrix (in the uncoupled basis of EasySpin)
S = 1/sqrt(2)*[0 1 -1 0];
Sys.initState = S'*S;

Opt = [];
Opt.Output = 'separate';

[f,spc,tr] = pepper(Sys,Exp,Opt);

plot(f,spc,f,sum(spc),'k')

nElStates = prod(2*Sys.S+1);
leg = [ num2str(tr(:,1)), repmat(' < - > ',nElStates,1),num2str(tr(:,2))];
legend(leg)
xlabel('Magnetic Field (mT)')