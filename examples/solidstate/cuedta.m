% Mixture of two copper-EDTA complexes
%----------------------------------------
% Cu-EDTA at pH 7 gives an EPR spectrum with two
% components. It can be simulated with a single
% call to pepper.

clear

% Component 1
CuEDTA1.g = [2.04 2.321];
CuEDTA1.A = [15 136]*2.8;     % MHz
CuEDTA1.Nucs = 'Cu';
CuEDTA1.lwpp = 2;             % mT

% Component 2
CuEDTA2.g = [2.032 2.288];
CuEDTA2.A = [15 144]*2.8;     % MHz
CuEDTA2.Nucs = 'Cu';
CuEDTA2.lwpp = 2;             % mT

% Relative abundances
CuEDTA1.weight = 1;
CuEDTA2.weight = 0.6;

% Experimental parameters
Xband.mwFreq = 9.5;          % GHz
Xband.Range = [250 360];     % mT

% One call to pepper with both Cu1 and Cu2 gives
% directly the two-component spectrum.
pepper({CuEDTA1,CuEDTA2},Xband);
