% Recombination triplet, time-resolved powder EPR spectrum
%-------------------------------------------------------------------------------
clc, clf, clear

% Spin Hamiltonian parameters and broadening
Triplet.S = 1;
Triplet.D = 500*[1 0.1];   % MHz
Triplet.lwpp = 1; % mT

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [310 370]; % mT
Exp.Harmonic = 0; % time-resolved EPR: no field modulation

Triplet.PopMode = 'highfield';
Triplet.Pop = [0 1 0];

[B,spc] = pepper(Triplet,Exp);

% Plotting energy level diagrams for z and x orientation sand powder spectrum
subplot(4,1,1)
levelsplot(Triplet,'x',Exp.Range,Exp.mwFreq);

subplot(4,1,2)
levelsplot(Triplet,'y',Exp.Range,Exp.mwFreq);

subplot(4,1,3)
levelsplot(Triplet,'z',Exp.Range,Exp.mwFreq);

subplot(4,1,4)
plot(B,spc,'k');
xlabel('magnetic field (mT)');
