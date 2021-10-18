% Recombination triplet, time-resolved powder EPR spectrum
%-------------------------------------------------------------------------------
clc, clf, clear

% Spin Hamiltonian parameters and broadening
Triplet.S = 1;
Triplet.D = 500*[1 0.1];   % MHz
Triplet.lwpp = 1; % mT

% Vector of populations, lowest to highest energy level
Triplet.PopMode = 'highfield';
Triplet.Pop = [0 1 0]; % vector of populations for the three levels

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [310 370]; % mT
Exp.Harmonic = 0;                % time-resolved EPR: no field modulation

Opt.Output = 'separate';         % separate transitions!

[B,spcsep,tr] = pepper(Triplet,Exp,Opt);
% Spectra in rows of spcsep correspond to transitions in rows of tr.

% EPR spectrum is sum of transition spectra
spc = sum(spcsep); % polarized spectrum

% Plotting energy level diagrams for z and x orientation sand powder spectrum
subplot(4,1,1)
levelsplot(Triplet,'x',Exp.Range,Exp.mwFreq);

subplot(4,1,2)
levelsplot(Triplet,'y',Exp.Range,Exp.mwFreq);

subplot(4,1,3)
levelsplot(Triplet,'z',Exp.Range,Exp.mwFreq);

subplot(4,1,4)
plot(B,spcsep,B,spc,'k');
legend('1-2','2-3','total');
legend boxoff
xlabel('magnetic field (mT)');
