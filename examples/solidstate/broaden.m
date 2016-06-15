% line broadening models for cw EPR spectra
%==========================================================================
clear, clf

% simple S=1/2 system with orthorhombic g matrix, at X band
Sys.g = [1.9,2.01,2.3];
Exp.Range = [285 365];    % mT
Exp.mwFreq = 9.5;         % GHz
Exp.Harmonic = 0;

% (1) broadening due to residual hyperfine couplings
Sys.lw = 0;
Sys.HStrain = [110 40 50];    % MHz, largest along gx
[x,y1] = pepper(Sys,Exp);

% (2) convolution broadening in the magnetic field domain
Sys.lw = 5;                   % mT
Sys.HStrain = [0 0 0];
[x,y2] = pepper(Sys,Exp);

% (3) g strain: Gaussian distribution of g principal values
Sys.lw = 0;
Sys.HStrain = [0 0 0];
Sys.gStrain = [0.05, 0.01, 0.01];
[x,y3] = pepper(Sys,Exp);

% plotting, normalizing all spectra to their integral
subplot(2,1,1); plot(x,y1,x,y2,x,y3); set(gca,'YTick',[]);
title('Different broadening models');
legend('HStrain','lw','gStrain');

x = x(1:end-1);
subplot(2,1,2);
plot(x,diff(y1),x,diff(y2),x,diff(y3)); set(gca,'YTick',[]);
legend('HStrain','lw','gStrain'); xlabel('magnetic field [mT]');
