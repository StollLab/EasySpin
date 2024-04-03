% Simple isotope mixture
%==========================================================================
% This examples shows how various isotopes mixtures can be specified via
% Sys.Nucs.

clear, clc

Sys.A = [12 50];  % isotropic hyperfine coupling values, MHz
Sys.n = [1 1];    % number of equivalent nuclei
Sys.lwpp = 0.05; % mT

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [330 348];  % mT
Exp.nPoints = 1e4;

% Natural-abundance C, pure 11B
Sys.Nucs = '12C,11B';
[B,y1] = garlic(Sys,Exp);

% Natural abundance C and natural abundance B
Sys.Nucs = 'C,B';
[B,y0] = garlic(Sys,Exp);

plot(B,y1,B,y0);
legend('natural abundance C, pure ^{11}B','natural abundance C and B');
legend boxoff
