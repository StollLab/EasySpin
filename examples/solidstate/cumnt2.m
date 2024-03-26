% L-band single-crystal EPR spectrum of a Cu(II) complex
%===============================================================================
% Values from Kirmse et al., Inorg.Chem. 23, 3333-3338 (1984)
% https://doi.org/10.1021/ic00189a012

clear, clc, clf

% Spin system
Sys.Nucs = '(63,65)Cu';
Sys.Abund = [0.975 0.025];  % 97.5% 63Cu enrichment
Sys.lwpp = 0.3;  % mT
Sys.g = [2.020 2.023 2.089];
A_cm = [-43.1 -43.7 -171.4]*1e-4; % cm^-1
Sys.A = unitconvert(A_cm,'cm^-1->MHz'); % conversion cm^-1 -> MHz

% Experimental parameters
Exp.mwFreq = 1.4;     % GHz
Exp.Range = [10 90];  % mT
Exp.Harmonic = 0;

% Crystal
Exp.CrystalSymmetry = 'P-1';  % space group no. 2 (triclinic)
Exp.MolFrame = [0 pi/3 0];   % assumed (rad)

% Options
Opt.Verbosity = 0;

% Simulation of two spectra with different crystal orientations
Exp.SampleFrame = [0 0 0];
[B,spec1] = pepper(Sys,Exp,Opt);
Exp.SampleFrame = [0 pi/2 0];
[B,spec2] = pepper(Sys,Exp,Opt);

% Plotting
plot(B,spec1,B,spec2);
axis tight
legend('0^o','90^o');
xlabel('magnetic field (mT)');
title(sprintf('%g GHz single-crystal spectra of Cu(mnt)_2',Exp.mwFreq));
