% Powder spectrum of spin-polarized triplet (simple)
%==========================================================================

clear, clf, clc

% Spin system and experimental parameters
Triplet.S = 1;
Triplet.g = 2.0023;
Triplet.D = [200 -30]; % MHz
Triplet.lw = 0.2; % mT

% Population vector for Tx,Ty,Tz states
% (The triplet is assumed to be generated via inter-system crossing.)
Triplet.initState = {[0.48 0.52 0],'xyz'};

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [328 350];  % mT
Exp.Harmonic = 0;  % no field modulation

[B,spc_polarized] = pepper(Triplet,Exp);

plot(B,spc_polarized);
xlabel('magnetic field (mT)');
