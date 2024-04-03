%
% Simulation for a coupled photoexcited triplet and doublet
%

clear; clc; clf;

% Spin system definition
% ---------------------------------------------------------------
Sys.S = [1 1/2];

Sys.g = [2.0023 2.0023 2.0023; ... % g-values triplet
         2.0090 2.0054 2.0011];    % g-values doublet

Sys.D = [2800 -500; 0 0]; % MHz, zero-field splitting triplet

% Coupling parameters
Sys.J = -4e5; % MHz
Sys.dip = 160; % MHz

Sys.lw = 2.5; % mT

% Triplet zero-field populations
popxyz = [0.3321 0.3521 0.3158];

% Construct initState for system of coupled triplet and doublet
% ---------------------------------------------------------------
% Zero-field eigenstates Tx, Ty and Tz in uncoupled basis
Tx = (1/sqrt(2))*[1;0;-1];
Ty = (1/sqrt(2))*[1;0;1];
Tz = [0;1;0];
ZFStates = [Tx Ty Tz];

% Initial density matrix for triplet state in uncoupled basis
initStateTriplet = ZFStates*diag(popxyz)*ZFStates';

% Initial density matrix for full system (with equal populations of doublet levels)
initState = kron(initStateTriplet,eye(2));

Sys.initState = {initState,'uncoupled'};

% Experimental parameters
% ---------------------------------------------------------------
Exp.mwFreq = 9.75; % GHz
Exp.Range = [270 430]; % mT
Exp.Harmonic = 0;

% Calculate and plot spectrum
% ---------------------------------------------------------------
[B,spec] = pepper(Sys,Exp);

plot(B,spec);
xlabel('magnetic field (mT)')