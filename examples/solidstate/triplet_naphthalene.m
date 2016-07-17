% Photo-excited triplet state of naphthalene
%==========================================================================
% see
%  - Hutchison et al, J.Chem.Phys. 34, 908-922 (1961)
%    http://dx.doi.org/10.1063/1.1731693
%  - Wasserman et al, J.Chem.Phys. 41, 1763-1772 (1964)
%    http://dx.doi.org/10.1063/1.1726156
%  - Weil/Bolton, chapter 6.3.4

% Here we simulate the X-band EPR powder spectrum of the photo-excited
% triplet state of naphthalene. The populations of the three sublevels
% are in thermal equilibrium.

clear, clf, clc

convert = 100*clight/1e6; % cm^-1 -> MHz conversion factor

Naphthalene.S = 1;
Naphthalene.g = 2.003;
Naphthalene.D = [0.1003 0.0137]*convert;
Naphthalene.lwpp = 6;

Xband.mwFreq = 9.5;
Xband.Range = [0 500];

pepper(Naphthalene,Xband);
