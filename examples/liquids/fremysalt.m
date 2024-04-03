% Fremy's salt (PADS), motional narrowing
%--------------------------------------------------
% Values taken from Goldman et al, J.Chem.Phys 56(2), 716 (1972)

clear, clf, clc

% PADS in 85% glycerol/15% H2O
Fremy.g = [2.00785 2.00590 2.00265];
Fremy.Nucs = '14N';
Fremy.A = unitconvert([5.5 5.0 28.7]/10,'mT->MHz'); % G -> MHz

% Rotational motion at 30 degree Celsius
Fremy.tcorr = 8e-10;    % seconds
Fremy.lw = [0, 0.03];   % mT

% X band experimental conditions
Exp.mwFreq = 9.5;             % GHz
Exp.CenterSweep = [338.5 5];  % mT

% Simulation
garlic(Fremy,Exp);
