% Photoexcited quartet state
%===========================================================
% Photoexcitation of a PDI-TEMPO molecule leads to a spin-polarised quartet state 
% 
% see
%  - Mayländer et al, J. Am. Chem. Soc., 143, 7050–7058 (2021)
%    https://doi.org/10.1021/jacs.1c01620
%

clear; clc; clf;

% Spin system parameters
Sys.S = [1 1/2];
Sys.g = [2.0043 2.0043 2.0027;  % PDI
         2.0090 2.0054 2.0011]; % TEMPO

% ZFS interaction of PDI triplet
Sys.D = [990 -300; 0 0];
Sys.DStrain = [100 100; 0 0];

% Coupling between triplet and doublet states
Sys.J = -1.2e3;  % MHz
Sys.dip = -66;%[-65 -65 130];  % MHz

Sys.Nucs = '14N';
Sys.A = [0 0 19 96];  % MHz

Sys.lw = 0.65;  % mT

% High-field state populations (different populations of nuclear sublevels)
Pop = [0 0 0 0.50 0.63 0.77 0.13 0.00 0.80 0.85 0.54 0.48 0.04 0.96 0.04 0 0 0];
Sys.initState = {Pop,'eigen'};

% Experimental parameters
Exp.Range = [1165 1260];  % mT
Exp.mwFreq = 34;  % GHz
Exp.Harmonic = 0;

% Spectral simulation
[B,spc] = pepper(Sys,Exp);
 
% Plot result
title('PDI-TEMPO quartet state')
hold on; box on;
plot(B,spc)
axis tight
xlabel('Magnetic Field (mT)')
grid on
