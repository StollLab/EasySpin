% Low-spin heme: g strain vs other line broadening models
%==========================================================================
% Low-spin ferric heme, as in substrate-free cytochrome P450cam, has a
% rhombic g tensor. In frozen solution, its lines are broadened mostly by a
% distribution of g values (g strain) due to structural heterogeneity
% around the iron. This example compares this model with two simpler ones.

clear, clc

% Low-spin Fe(III) heme (S = 1/2), rhombic g tensor
Sys.g = [1.91 2.26 2.45];

% X-band spectrum
Exp.Range = [250 380];    % mT
Exp.mwFreq = 9.5;         % GHz
Exp.Harmonic = 0;

% (1) HStrain: orientation-dependent, frequency-independent Gaussian width
Sys.lw = 0;
Sys.HStrain = [400 150 250];  % MHz
[B,spc1] = pepper(Sys,Exp);

% (2) convolution broadening in the magnetic field domain
Sys.lw = 5;                   % mT
Sys.HStrain = [0 0 0];        % MHz
[B,spc2] = pepper(Sys,Exp);

% (3) g strain: Gaussian distribution of g principal values
Sys.lw = 0;
Sys.HStrain = [0 0 0];       % MHz
Sys.StrainPars = {'g(1)','g(2)','g(3)'};
Sys.StrainFWHM = [0.04 0.015 0.025];
[B,spc3] = pepper(Sys,Exp);

% Plotting, normalizing all spectra to their integral
subplot(2,1,1);
plot(B,spc1,B,spc2,B,spc3);
title('Different broadening models');
legend('HStrain','lw','g strain');

subplot(2,1,2);
plot(B,deriv(spc1),B,deriv(spc2),B,deriv(spc3));
legend('HStrain','lw','g strain');
xlabel('magnetic field (mT)');
