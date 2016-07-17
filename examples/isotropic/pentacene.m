% Pentacene radical anion and cation
%==========================================================================
% Values taken from Weil/Bolton/Wertz book, Table 9.3, page 249

clear, clf

% Parameters
%--------------------------------------------------------------------------
% Spin system
Sys = struct('g',2,'Nucs','1H,1H,1H,1H','n',[2 4 4 4]);
Sys.lwpp = [0,0.01]; % Lorentzian lines

% Hyperfine couplings
A_anion =  mt2mhz([0.4263 0.3032 0.0915 0.0870]);
A_cation = mt2mhz([0.5083 0.3554 0.0975 0.0757]);

Exp.mwFreq = 9.5;
Exp.CenterSweep = [339.4, 4];
Exp.nPoints = 1e4;

% Simulations
%--------------------------------------------------------------------------
Sys.A = A_anion;
[xa,ya] = garlic(Sys,Exp);
Sys.A = A_cation;
[xc,yc] = garlic(Sys,Exp);

% Graphical rendering
%--------------------------------------------------------------------------
plot(xa,ya,'r',xc,yc,'b');
legend('radical anion','radical cation');
xlabel('magnetic field [mT]');
axis tight
title('pentacene radical spectra');
