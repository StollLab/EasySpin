% Powder spectrum of spin-polarized triplet (simple)
%==========================================================================

clear, clf, clc

% Spin system and experimental parameters
Triplet.S = 1;
Triplet.g = 2;
Triplet.D = 100; % MHz
Triplet.lw = 0.2; % MHz

% Population vector for Tx,Ty,Tz states
% (The triplet is assumed to be generated via inter-system crossing)
Triplet.Pop = [0.85 1 0.95];

Exp.mwFreq = 9.5; % GHz
Exp.Range = [330 350]; % mT
Exp.Harmonic = 0;

[x,y_polarized] = pepper(Triplet,Exp);

plot(x,y_polarized);
xlabel('magnetic field [mT]');
