% Time-resolved EPR spectrum of recombination triplet in photosystem II
%-------------------------------------------------------------------------------

% Recombination triplet in photosystem II (3P680)
%
% see
% Lendzian et al, Biochim. Biophys. Acta Bioenergetics, 1605, 35-46 (2003)
% https://doi.org/10.1016/S0005-2728(03)00062-8
% and references 13, 14 and 21 cited therein

clc, clf, clear

% Spin Hamiltonian parameters and broadening
Triplet.S = 1;
Triplet.D = [287e-4 43e-4]*30e3;  % cm^-1 -> MHz
Triplet.lwpp = 1;  % mT

% Experimental parameters
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [300 380];  % mT
Exp.Harmonic = 0;  % time-resolved EPR: no field modulation

% Shortcut for recombination triplet: exclusive population of T0 level
Triplet.initState = 'T0';

% General input: population vector in eigenbasis, in order of increasing energy
% Triplet.initState = {[0 1 0],'eigen'};

[B,spc] = pepper(Triplet,Exp);

% Plotting energy level diagrams for x, y and z orientations and powder spectrum
subplot(4,1,1)
levelsplot(Triplet,'x',Exp.Range,Exp.mwFreq);
title('B0 || x')

subplot(4,1,2)
levelsplot(Triplet,'y',Exp.Range,Exp.mwFreq);
title('B0 || y')

subplot(4,1,3)
levelsplot(Triplet,'z',Exp.Range,Exp.mwFreq);
title('B0 || z')

subplot(4,1,4)
plot(B,spc,'k');
xlabel('magnetic field (mT)');
