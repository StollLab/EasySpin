% Two coupled electrons
%==========================================================================
clear, clf

% Two-electron system, both electrons with rhombic g tensor
Sys.S = [1/2 1/2];
Sys.g = [2 2.05 2.1; 2.2 2.25 2.3];
Sys.lwpp = 1;

% Electron-electron couplings
Sys.J = 50;     % isotropic coupling, in MHz
Sys.dip = 100;  % axial dipolar coupling, in MHz

% X band conditions
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 350];  % mT

pepper(Sys,Exp);

title('Two dipolar coupled orthorhombic S=1/2');
