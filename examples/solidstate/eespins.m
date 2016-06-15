% Two coupled electrons
%==========================================================================
clear, clf

% two-electron system, both electrons with rhombic g
Sys.S = [1/2 1/2];
Sys.g = [2 2.05 2.1; 2.2 2.25 2.3];
Sys.lwpp = 1;
Sys.J = 50; % isotropic exchange coupling, in MHz
Sys.eeD = [1 1 -2]*100; % dipolar electron-electron coupling, in MHz

% X band conditions
Exp.mwFreq = 9.5;
Exp.Range = [280 350];

pepper(Sys,Exp);

title('Two dipolar coupled orthorhombic S=1/2');
