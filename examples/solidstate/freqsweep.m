% basic frequency-swept EPR spectrum
%==========================================================================
clear, clf

% Spin system: nitroxide radical
%------------------------------------------------------------
Nx.g = [2.008 2.006 2.002];
Nx.Nucs = '14N';
Nx.A = [20 20 90]; % MHz

% For frequency sweeps, the convolutional width should be given in MHz.
Nx.lwpp = 5; % MHz

% Frequency-sweep experiment
%------------------------------------------------------------
Exp1.Field = 1100;         % static magnetic field, mT
Exp1.mwRange = [30.7 31];  % frequency range, GHz

% Simulate and plot frequency-sweep spectrum
pepper(Nx,Exp1);
