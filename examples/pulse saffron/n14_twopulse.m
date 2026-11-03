% Two-pulse ESEEM of a 14N nucleus in cancellation regime (saffron)
%==========================================================================
% Reproduces Fig. 4(a) of
%   H. L. Flanagan, D. J. Singel, J. Chem. Phys. 87, 5606-5616 (1987)
%   https://doi.org/10.1063/1.453532
%
% In this example, the isotropic hyperfine coupling is twice the 14N Larmor
% frequency (A = 2*nuI). This is the exact cancellation regime: in one electron
% spin manifold, the hyperfine field cancels the external field at the nucleus,
% and in the other manifold the two add up to an effective nuclear Zeeman
% frequency of nuI + A/2 = 2*nuI (2 MHz here).
%
% Peaks in the spectrum:
% - 0.12, 0.24, 0.36 MHz: sharp NQR lines of the cancellation manifold,
%   2*K*eta and K*(3-+eta), with K = e^2qQ/(4h) = 0.1 MHz, eta = 0.6
% - 4.01 MHz: double-quantum line of the other manifold, about 2*(nuI + A/2)
% - 3.65-4.37 MHz: negative sum and difference combinations of the
%   double-quantum line with the NQR lines
% - Around 2 MHz: broad overlapping single-quantum lines of the other manifold,
%   and their negative sum and difference combinations with the NQR lines

clear, clf, clc

% Two-pulse experiment parameters
Exp.Field = 325.07;  % mT
Exp.Sequence = '2pESEEM';
Exp.dt = 0.100;  % µs
Exp.nPoints = 501; 

% Spin system with S=1/2 and one 14N
Sys.Nucs = '14N';
nuI = larmorfrq(Sys.Nucs,Exp.Field);  % MHz
Sys.A = 2*nuI;  % MHz
Sys.Q = [4*0.1*nuI, 0.6];  % MHz

Opt.GridSize = 91;

saffron(Sys,Exp,Opt);
