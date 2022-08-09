% Powder spectrum of spin-polarized triplet with axial ZFS
%==========================================================================

clear, clf, clc

% Spin system and experimental parameters
Sys.S = 1;
Sys.g = 2.0023;
Sys.D = 2500; % MHz
Sys.lw = 1; % mT

% Population vector for zero-field states in order of increasing energy
% (The triplet is assumed to be generated via inter-system crossing.)
zfpop = [0 0 1];
ZFStates = [0 sqrt(2)/2 sqrt(2)/2; 1 0 0; 0 -sqrt(2)/2 sqrt(2)/2];

% Convert population vector to density matrix
Sys.initState = ZFStates*diag(zfpop)*ZFStates';

Exp.mwFreq = 9.5; % GHz
Exp.Range = [100 450]; % mT
Exp.Harmonic = 0;

[B,spec] = pepper(Sys,Exp);

plot(B,spec);
xlabel('magnetic field (mT)');
