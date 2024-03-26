% line broadening models for cw EPR spectra
%==========================================================================
clear, clc

% Simple S=1/2 system with orthorhombic g matrix
Sys.g = [1.9,2.01,2.3];

% X-band spectrum
Exp.Range = [280 370];    % mT
Exp.mwFreq = 9.5;         % GHz
Exp.Harmonic = 0;

% (1) broadening due to residual hyperfine couplings
Sys.lw = 0;
Sys.HStrain = [110 40 50];    % MHz, largest along gx
[B,spc1] = pepper(Sys,Exp);

% (2) convolution broadening in the magnetic field domain
Sys.lw = 5;                   % mT
Sys.HStrain = [0 0 0];        % MHz
[B,spc2] = pepper(Sys,Exp);

% (3) g strain: Gaussian distribution of g principal values
Sys.lw = 0;
Sys.HStrain = [0 0 0];       % MHz
Sys.gStrain = [0.05, 0.01, 0.01];
[B,spc3] = pepper(Sys,Exp);

% Plotting, normalizing all spectra to their integral
subplot(2,1,1);
plot(B,spc1,B,spc2,B,spc3);
title('Different broadening models');
legend('HStrain','lw','gStrain');

subplot(2,1,2);
plot(B,deriv(spc1),B,deriv(spc2),B,deriv(spc3));
legend('HStrain','lw','gStrain');
xlabel('magnetic field (mT)');
